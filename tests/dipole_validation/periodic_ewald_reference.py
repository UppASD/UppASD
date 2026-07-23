#!/usr/bin/env python3
"""Independent fp64 3D-periodic Ewald reference for point magnetic dipoles.

This deliberately small O(N^2 (R + K)) evaluator is a validation oracle for
the block-one macrocell limit.  It uses no UppASD routines and no GPU code.
Coordinates and cell vectors are expressed in one common length unit; the
returned field has the corresponding dipole-tensor prefactor omitted.

The convention is 3D periodic Ewald with conducting (tin-foil) exterior:
only the physical k=0 vector is omitted.  It implements
``B = Hessian(1/r) m`` and ``E = -1/2 sum_i m_i dot B_i``.

The displacement-kernel builder at the end of this file is intentionally a
second formulation of the same physics.  It constructs each target-minus-
source displacement directly, rather than obtaining it by evaluating a
many-particle system.  Keeping both formulations here makes the kernel and
the reciprocal-alias regression independent of UppASD C++ code.
"""
from __future__ import annotations

from dataclasses import dataclass
from itertools import product
import math
from typing import Mapping, Sequence


Vec = tuple[float, float, float]
Mat = tuple[Vec, Vec, Vec]  # columns are the three cell vectors
Tensor = tuple[float, float, float, float, float, float]  # xx,yy,zz,xy,xz,yz
Grid = tuple[int, int, int]
KernelKey = tuple[tuple[int, int, int], int, int]
Kernel = dict[KernelKey, Tensor]


def add(a: Vec, b: Vec) -> Vec:
    return (a[0] + b[0], a[1] + b[1], a[2] + b[2])


def sub(a: Vec, b: Vec) -> Vec:
    return (a[0] - b[0], a[1] - b[1], a[2] - b[2])


def scale(a: Vec, value: float) -> Vec:
    return (value * a[0], value * a[1], value * a[2])


def dot(a: Vec, b: Vec) -> float:
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]


def cross(a: Vec, b: Vec) -> Vec:
    return (a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2],
            a[0] * b[1] - a[1] * b[0])


def cell_vector(cell: Mat, n: tuple[float, float, float]) -> Vec:
    return add(add(scale(cell[0], n[0]), scale(cell[1], n[1])),
               scale(cell[2], n[2]))


def reciprocal_cell(cell: Mat) -> tuple[Mat, float]:
    """Return ``2*pi H^-T`` as columns and the positive cell volume."""
    volume = dot(cell[0], cross(cell[1], cell[2]))
    if not math.isfinite(volume) or volume <= 0.0:
        raise ValueError("cell must be right-handed with positive volume")
    factor = 2.0 * math.pi / volume
    return (scale(cross(cell[1], cell[2]), factor),
            scale(cross(cell[2], cell[0]), factor),
            scale(cross(cell[0], cell[1]), factor)), volume


def reciprocal_identity_error(cell: Mat, reciprocal: Mat | None = None) -> float:
    """Return ``max(abs(H.T @ Brec - 2*pi*I))`` in fp64."""
    if reciprocal is None:
        reciprocal, _ = reciprocal_cell(cell)
    return max(abs(dot(cell[row], reciprocal[column]) -
                   (2.0 * math.pi if row == column else 0.0))
               for row in range(3) for column in range(3))


def assert_reciprocal_identity(cell: Mat, tolerance: float = 2e-14) -> None:
    error = reciprocal_identity_error(cell)
    if error > tolerance:
        raise AssertionError(f"H^T*Brec identity residual {error:.3e} > {tolerance:.3e}")


def real_tensor(r: Vec, alpha: float) -> Tensor:
    """Symmetric Hessian of erfc(alpha*r)/r."""
    r2 = dot(r, r)
    if r2 == 0.0:
        raise ValueError("real Ewald tensor is undefined at r=0")
    radius = math.sqrt(r2)
    erfc = math.erfc(alpha * radius)
    gaussian = math.exp(-(alpha * alpha) * r2)
    # Hessian f = A r r^T + B I for f=erfc(alpha*r)/r.
    a = (3.0 * erfc / (radius**5) +
         6.0 * alpha * gaussian / (math.sqrt(math.pi) * radius**4) +
         4.0 * alpha**3 * gaussian / (math.sqrt(math.pi) * radius**2))
    b = (-erfc / (radius**3) -
         2.0 * alpha * gaussian / (math.sqrt(math.pi) * radius**2))
    x, y, z = r
    return (a * x * x + b, a * y * y + b, a * z * z + b,
            a * x * y, a * x * z, a * y * z)


def tensor_add(left: Tensor, right: Tensor) -> Tensor:
    return tuple(a + b for a, b in zip(left, right))  # type: ignore[return-value]


def tensor_scale(tensor: Tensor, value: float) -> Tensor:
    return tuple(value * item for item in tensor)  # type: ignore[return-value]


def tensor_times(tensor: Tensor, moment: Vec) -> Vec:
    xx, yy, zz, xy, xz, yz = tensor
    x, y, z = moment
    return (xx * x + xy * y + xz * z,
            xy * x + yy * y + yz * z,
            xz * x + yz * y + zz * z)


def reciprocal_tensor(wave: Vec, volume: float, alpha: float) -> Tensor:
    wave2 = dot(wave, wave)
    if wave2 == 0.0:
        raise ValueError("the physical reciprocal k=0 term is omitted")
    factor = -4.0 * math.pi * math.exp(-wave2 / (4.0 * alpha * alpha)) / (volume * wave2)
    x, y, z = wave
    return (factor * x * x, factor * y * y, factor * z * z,
            factor * x * y, factor * x * z, factor * y * z)


@dataclass(frozen=True)
class ConvergenceReport:
    alpha: float
    real_images: tuple[int, int, int]
    reciprocal_images: tuple[int, int, int]
    real_shell_changes: tuple[float, ...]
    reciprocal_shell_changes: tuple[float, ...]


@dataclass(frozen=True)
class EwaldResult:
    fields: tuple[Vec, ...]
    energy: float
    convergence: ConvergenceReport | None = None


def _validate_system(positions: Sequence[Vec], moments: Sequence[Vec], alpha: float,
                     real_images: tuple[int, int, int],
                     reciprocal_images: tuple[int, int, int]) -> None:
    if len(positions) != len(moments) or not positions:
        raise ValueError("positions and moments must be non-empty and have equal length")
    if (not math.isfinite(alpha) or alpha <= 0.0 or
            min(*real_images, *reciprocal_images) < 0):
        raise ValueError("alpha must be positive and image ranges non-negative")


def _zero_fields(count: int) -> list[Vec]:
    return [(0.0, 0.0, 0.0) for _ in range(count)]


def _real_fields(positions: Sequence[Vec], moments: Sequence[Vec], cell: Mat,
                 alpha: float, real_images: tuple[int, int, int]) -> tuple[Vec, ...]:
    fields = _zero_fields(len(positions))
    image_range = product(*(range(-extent, extent + 1) for extent in real_images))
    images = tuple(image_range)
    for target, target_pos in enumerate(positions):
        field = (0.0, 0.0, 0.0)
        for source, source_pos in enumerate(positions):
            for image in images:
                if target == source and image == (0, 0, 0):
                    continue
                r = sub(target_pos, add(source_pos, cell_vector(cell, image)))
                field = add(field, tensor_times(real_tensor(r, alpha), moments[source]))
        fields[target] = field
    return tuple(fields)


def _reciprocal_fields(positions: Sequence[Vec], moments: Sequence[Vec], cell: Mat,
                       alpha: float,
                       reciprocal_images: tuple[int, int, int]) -> tuple[Vec, ...]:
    reciprocal, volume = reciprocal_cell(cell)
    fields = _zero_fields(len(positions))
    for wave_index in product(*(range(-extent, extent + 1) for extent in reciprocal_images)):
        if wave_index == (0, 0, 0):
            continue  # tin-foil: omit only the physical reciprocal k=0
        wave = cell_vector(reciprocal, wave_index)
        prefactor = -4.0 * math.pi * math.exp(-dot(wave, wave) / (4.0 * alpha * alpha)) / (volume * dot(wave, wave))
        for target, target_pos in enumerate(positions):
            field = fields[target]
            for source, source_pos in enumerate(positions):
                phase = math.cos(dot(wave, sub(target_pos, source_pos)))
                field = add(field, scale(wave, prefactor * phase * dot(wave, moments[source])))
            fields[target] = field
    return tuple(fields)


def _sum_fields(*terms: Sequence[Vec]) -> tuple[Vec, ...]:
    if not terms:
        return ()
    return tuple(tuple(sum(term[index][component] for term in terms)
                       for component in range(3))  # type: ignore[misc]
                 for index in range(len(terms[0])))


def _max_field_difference(first: Sequence[Vec], second: Sequence[Vec]) -> float:
    return max((abs(left - right) for lhs, rhs in zip(first, second)
                for left, right in zip(lhs, rhs)), default=0.0)


def _max_field_component(fields: Sequence[Vec]) -> float:
    return max((abs(value) for field in fields for value in field), default=0.0)


def _converged(change: float, current: float, tolerance: float) -> bool:
    return change <= tolerance * max(1.0, current)


def _evaluate_with_terms(positions: Sequence[Vec], moments: Sequence[Vec], cell: Mat,
                         alpha: float, real_images: tuple[int, int, int],
                         reciprocal_images: tuple[int, int, int]) -> EwaldResult:
    real = _real_fields(positions, moments, cell, alpha, real_images)
    reciprocal = _reciprocal_fields(positions, moments, cell, alpha, reciprocal_images)
    self_prefactor = 4.0 * alpha**3 / (3.0 * math.sqrt(math.pi))
    self_term = tuple(scale(moment, self_prefactor) for moment in moments)
    fields = _sum_fields(real, reciprocal, self_term)
    energy = -0.5 * sum(dot(moment, field) for moment, field in zip(moments, fields))
    return EwaldResult(fields=fields, energy=energy)


def evaluate(positions: Sequence[Vec], moments: Sequence[Vec], cell: Mat, *, alpha: float,
             real_images: tuple[int, int, int],
             reciprocal_images: tuple[int, int, int]) -> EwaldResult:
    """Evaluate a truncated 3D tin-foil Ewald sum in fp64.

    The explicit ranges remain part of this low-level API so convergence
    studies can show exactly which real and reciprocal shells were included.
    Use :func:`evaluate_converged` when the caller wants automatic shell
    selection and diagnostics.
    """
    _validate_system(positions, moments, alpha, real_images, reciprocal_images)
    assert_reciprocal_identity(cell)
    return _evaluate_with_terms(positions, moments, cell, alpha, real_images, reciprocal_images)


def default_alpha(cell: Mat) -> float:
    """Choose a deterministic split parameter from the full periodic cell."""
    shortest = min(math.sqrt(dot(vector, vector)) for vector in cell)
    if shortest <= 0.0:
        raise ValueError("cell vectors must be non-zero")
    return math.sqrt(math.pi) / shortest


def _auto_field_shells(positions: Sequence[Vec], moments: Sequence[Vec], cell: Mat,
                       alpha: float, tolerance: float, max_shell: int,
                       real: bool) -> tuple[tuple[int, int, int], tuple[float, ...]]:
    """Converge one Ewald side while measuring that side only.

    The other side is deliberately absent.  This prevents real/reciprocal
    truncation errors from cancelling and passing a total-field comparison.
    """
    previous = _zero_fields(len(positions))
    changes: list[float] = []
    stable = 0
    for shell in range(max_shell + 1):
        extent = (shell, shell, shell)
        current = (_real_fields(positions, moments, cell, alpha, extent)
                   if real else _reciprocal_fields(positions, moments, cell, alpha, extent))
        change = _max_field_difference(current, previous)
        changes.append(change)
        if shell > 0 and _converged(change, _max_field_component(current), tolerance):
            stable += 1
            if stable == 2:
                return extent, tuple(changes)
        else:
            stable = 0
        previous = current
    side = "real" if real else "reciprocal"
    raise RuntimeError(f"{side} Ewald shells did not converge by {max_shell}")


def evaluate_converged(positions: Sequence[Vec], moments: Sequence[Vec], cell: Mat, *,
                       alpha: float | None = None, tolerance: float = 1e-12,
                       max_shell: int = 12) -> EwaldResult:
    """Evaluate after independently converging two consecutive shells.

    The reported shell changes are component-wise field maxima for the real
    and reciprocal sums separately.  ``alpha`` is optional only as a
    deterministic test convenience; production code must not expose it as a
    required user tuning parameter.
    """
    if not math.isfinite(tolerance) or tolerance <= 0.0:
        raise ValueError("tolerance must be finite and positive")
    if alpha is None:
        alpha = default_alpha(cell)
    _validate_system(positions, moments, alpha, (0, 0, 0), (0, 0, 0))
    assert_reciprocal_identity(cell)
    real_extent, real_changes = _auto_field_shells(
        positions, moments, cell, alpha, tolerance, max_shell, real=True)
    reciprocal_extent, reciprocal_changes = _auto_field_shells(
        positions, moments, cell, alpha, tolerance, max_shell, real=False)
    report = ConvergenceReport(alpha, real_extent, reciprocal_extent,
                               real_changes, reciprocal_changes)
    result = _evaluate_with_terms(positions, moments, cell, alpha,
                                  real_extent, reciprocal_extent)
    return EwaldResult(result.fields, result.energy, report)


def finite_difference_field(positions: Sequence[Vec], moments: Sequence[Vec], cell: Mat,
                            *, atom: int, component: int, step: float = 1e-6,
                            **evaluate_options: object) -> float:
    """Return ``-dE/dM[atom,component]`` by a central finite difference."""
    if not 0 <= atom < len(moments) or not 0 <= component < 3 or step <= 0.0:
        raise ValueError("finite-difference index or step is invalid")
    plus = [list(moment) for moment in moments]
    minus = [list(moment) for moment in moments]
    plus[atom][component] += step
    minus[atom][component] -= step
    plus_result = evaluate(positions, [tuple(moment) for moment in plus], cell, **evaluate_options)
    minus_result = evaluate(positions, [tuple(moment) for moment in minus], cell, **evaluate_options)
    return -(plus_result.energy - minus_result.energy) / (2.0 * step)


def grid_positions(grid: Grid, cell: Mat, basis_offsets: Sequence[Vec]) -> tuple[Vec, ...]:
    """Return cell-major, basis-fast positions in a full periodic cell ``cell``."""
    if len(basis_offsets) == 0 or min(grid) <= 0:
        raise ValueError("grid dimensions must be positive and basis non-empty")
    positions = []
    for index in product(*(range(size) for size in grid)):
        translation = cell_vector(cell, tuple(index[axis] / grid[axis] for axis in range(3)))
        positions.extend(add(translation, offset) for offset in basis_offsets)
    return tuple(positions)


def _kernel_component(grid: Grid, cell: Mat, basis_offsets: Sequence[Vec], alpha: float,
                     real_images: tuple[int, int, int],
                     reciprocal_images: tuple[int, int, int], *,
                     include_real: bool, include_reciprocal: bool,
                     include_self: bool) -> Kernel:
    reciprocal, volume = reciprocal_cell(cell)
    result: Kernel = {}
    for displacement in product(*(range(size) for size in grid)):
        translation = cell_vector(cell, tuple(displacement[axis] / grid[axis] for axis in range(3)))
        for target_basis, target_offset in enumerate(basis_offsets):
            target = add(translation, target_offset)
            for source_basis, source_offset in enumerate(basis_offsets):
                source = source_offset
                tensor: Tensor = (0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
                if include_real:
                    for image in product(*(range(-extent, extent + 1)
                                            for extent in real_images)):
                        r = sub(target, add(source, cell_vector(cell, image)))
                        if r == (0.0, 0.0, 0.0):
                            continue
                        tensor = tensor_add(tensor, real_tensor(r, alpha))
                if include_reciprocal:
                    displacement_vector = sub(target, source)
                    for wave_index in product(*(range(-extent, extent + 1)
                                                for extent in reciprocal_images)):
                        if wave_index == (0, 0, 0):
                            continue  # omit only physical k=0
                        wave = cell_vector(reciprocal, wave_index)
                        phase = math.cos(dot(wave, displacement_vector))
                        tensor = tensor_add(tensor, tensor_scale(
                            reciprocal_tensor(wave, volume, alpha), phase))
                if include_self and displacement == (0, 0, 0) and target_basis == source_basis:
                    tensor = tensor_add(tensor, tensor_scale(
                        (1.0, 1.0, 1.0, 0.0, 0.0, 0.0),
                        4.0 * alpha**3 / (3.0 * math.sqrt(math.pi))))
                result[(displacement, target_basis, source_basis)] = tensor
    return result


def _max_kernel_difference(first: Mapping[KernelKey, Tensor],
                           second: Mapping[KernelKey, Tensor]) -> float:
    return max(abs(left - right) for key in first for left, right in zip(first[key], second[key]))


def _max_kernel_component(kernel: Mapping[KernelKey, Tensor]) -> float:
    return max((abs(value) for tensor in kernel.values() for value in tensor), default=0.0)


def _auto_kernel_extent(grid: Grid, cell: Mat, basis_offsets: Sequence[Vec], alpha: float,
                        tolerance: float, max_shell: int, *, real: bool) -> tuple[tuple[int, int, int], tuple[float, ...]]:
    previous = _kernel_component(grid, cell, basis_offsets, alpha, (0, 0, 0), (0, 0, 0),
                                 include_real=real, include_reciprocal=not real,
                                 include_self=False)
    changes: list[float] = []
    stable = 0
    for shell in range(max_shell + 1):
        extent = (shell, shell, shell)
        current = _kernel_component(grid, cell, basis_offsets, alpha, extent, extent,
                                    include_real=real, include_reciprocal=not real,
                                    include_self=False)
        change = _max_kernel_difference(current, previous)
        changes.append(change)
        if shell > 0 and _converged(change, _max_kernel_component(current), tolerance):
            stable += 1
            if stable == 2:
                return extent, tuple(changes)
        else:
            stable = 0
        previous = current
    side = "real" if real else "reciprocal"
    raise RuntimeError(f"{side} kernel shells did not converge by {max_shell}")


def generate_displacement_kernel(grid: Grid, cell: Mat, basis_offsets: Sequence[Vec], *,
                                 alpha: float | None = None,
                                 real_images: tuple[int, int, int] | None = None,
                                 reciprocal_images: tuple[int, int, int] | None = None,
                                 tolerance: float = 1e-12,
                                 max_shell: int = 12) -> Kernel:
    """Build complete ``K(d,a,b)`` in target-minus-source convention.

    ``cell`` is the full periodic cell ``H``; grid-cell translations are
    ``H_i/grid_i``.  Basis offsets are Cartesian positions in that cell.  The
    returned dictionary is keyed by ``(d, target_basis, source_basis)`` and
    stores a packed symmetric 3x3 tensor.  If ranges are omitted, real and
    reciprocal shells are converged independently, just as in the oracle.
    """
    if alpha is None:
        alpha = default_alpha(cell)
    if real_images is None:
        real_images, _ = _auto_kernel_extent(grid, cell, basis_offsets, alpha, tolerance,
                                             max_shell, real=True)
    if reciprocal_images is None:
        reciprocal_images, _ = _auto_kernel_extent(grid, cell, basis_offsets, alpha, tolerance,
                                                   max_shell, real=False)
    return _kernel_component(grid, cell, basis_offsets, alpha, real_images, reciprocal_images,
                             include_real=True, include_reciprocal=True, include_self=True)


def apply_displacement_kernel(kernel: Mapping[KernelKey, Tensor], grid: Grid,
                              basis_count: int, moments: Sequence[Vec]) -> tuple[Vec, ...]:
    """Apply a displacement kernel to cell-major, basis-fast moments."""
    expected = math.prod(grid) * basis_count
    if len(moments) != expected:
        raise ValueError(f"expected {expected} moments, got {len(moments)}")
    fields = _zero_fields(expected)
    for target_cell in product(*(range(size) for size in grid)):
        target_cell_number = (target_cell[0] * grid[1] * grid[2] +
                              target_cell[1] * grid[2] + target_cell[2])
        for target_basis in range(basis_count):
            target_number = target_cell_number * basis_count + target_basis
            field = (0.0, 0.0, 0.0)
            for source_cell in product(*(range(size) for size in grid)):
                source_cell_number = (source_cell[0] * grid[1] * grid[2] +
                                      source_cell[1] * grid[2] + source_cell[2])
                displacement = tuple((target_cell[axis] - source_cell[axis]) % grid[axis]
                                     for axis in range(3))
                for source_basis in range(basis_count):
                    source_number = source_cell_number * basis_count + source_basis
                    field = add(field, tensor_times(
                        kernel[(displacement, target_basis, source_basis)], moments[source_number]))
            fields[target_number] = field
    return tuple(fields)


def kernel_energy(fields: Sequence[Vec], moments: Sequence[Vec]) -> float:
    return -0.5 * sum(dot(moment, field) for moment, field in zip(moments, fields))


def reciprocal_alias_tensor(grid: Grid, cell: Mat, q_index: tuple[int, int, int], *,
                            alpha: float, reciprocal_images: tuple[int, int, int]) -> Tensor:
    """Sum all reciprocal vectors aliasing one physical FFT bin."""
    reciprocal, volume = reciprocal_cell(cell)
    representative = tuple(index if index <= size // 2 else index - size
                           for index, size in zip(q_index, grid))
    result: Tensor = (0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
    for wave_index in product(*(range(-extent, extent + 1)
                                for extent in reciprocal_images)):
        if wave_index == (0, 0, 0):
            continue  # physical k=0 only
        if any((wave_index[axis] - representative[axis]) % grid[axis] != 0
               for axis in range(3)):
            continue
        result = tensor_add(result, reciprocal_tensor(
            cell_vector(reciprocal, wave_index), volume, alpha))
    return result


def single_representative_reciprocal_tensor(grid: Grid, cell: Mat,
                                            q_index: tuple[int, int, int], *,
                                            alpha: float) -> Tensor:
    """The incomplete one-vector construction used by the old FFT staging."""
    reciprocal, volume = reciprocal_cell(cell)
    representative = tuple(index if index <= size // 2 else index - size
                           for index, size in zip(q_index, grid))
    if representative == (0, 0, 0):
        return (0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
    return reciprocal_tensor(cell_vector(reciprocal, representative), volume, alpha)
