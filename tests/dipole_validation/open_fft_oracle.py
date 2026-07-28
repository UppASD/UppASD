#!/usr/bin/env python3
"""Independent finite/open point-dipole oracle for WP10.

This module deliberately has no UppASD, FFT, CUDA/HIP, or periodic-Ewald
dependency.  Geometry is evaluated from primitive vectors and basis offsets;
the tensor is dimensionless until :func:`apply_units` is called.
"""
from __future__ import annotations

from dataclasses import dataclass
import json
import math
from pathlib import Path
from typing import Iterable, Sequence


MU_B_J_PER_T = 9.2740100783e-24
MU0_OVER_4PI_T_M_PER_A = 1.0e-7
MRY_J = 2.179872325e-21

Vector = tuple[float, float, float]
Tensor = tuple[tuple[float, float, float], ...]


def cell_index(grid: tuple[int, int, int], x: int, y: int, z: int) -> int:
    """Return the x-fastest active-cell index frozen by WP10."""
    return x + grid[0] * (y + grid[1] * z)


def cell_coordinates(grid: tuple[int, int, int]) -> Iterable[tuple[int, int, int]]:
    for z in range(grid[2]):
        for y in range(grid[1]):
            for x in range(grid[0]):
                yield x, y, z


def _add(left: Vector, right: Vector) -> Vector:
    return tuple(left[i] + right[i] for i in range(3))  # type: ignore[return-value]


def _sub(left: Vector, right: Vector) -> Vector:
    return tuple(left[i] - right[i] for i in range(3))  # type: ignore[return-value]


def _scale(value: Vector, factor: float) -> Vector:
    return tuple(factor * value[i] for i in range(3))  # type: ignore[return-value]


def point_tensor(r: Vector) -> Tensor:
    """Analytic dimensionless tensor ``3 rr^T/|r|^5-I/|r|^3``."""
    r2 = sum(component * component for component in r)
    if r2 == 0.0:
        raise ValueError("point_tensor received the zero displacement")
    inv_r3 = r2 ** -1.5
    inv_r5 = inv_r3 / r2
    return tuple(
        tuple(3.0 * r[alpha] * r[beta] * inv_r5 - (inv_r3 if alpha == beta else 0.0)
               for beta in range(3))
        for alpha in range(3)
    )


def tensor_vector(tensor: Tensor, vector: Vector) -> Vector:
    return tuple(sum(tensor[row][column] * vector[column] for column in range(3))
                 for row in range(3))  # type: ignore[return-value]


def physical_prefactor(alat_m: float) -> float:
    """Tesla per dimensionless field for moments expressed in ``mu_B``."""
    if not math.isfinite(alat_m) or alat_m <= 0.0:
        raise ValueError(f"alat_m must be positive and finite, got {alat_m!r}")
    return MU0_OVER_4PI_T_M_PER_A * MU_B_J_PER_T / alat_m ** 3


@dataclass(frozen=True)
class FiniteCase:
    name: str
    grid: tuple[int, int, int]
    primitive_vectors: tuple[float, ...]  # column-major [C1 C2 C3]
    basis_offsets: tuple[Vector, ...]
    moments_mu_b: tuple[tuple[tuple[Vector, ...], ...], ...]  # ensemble/cell/basis/vector
    alat_m: float = 1.0e-9

    @property
    def active_cells(self) -> int:
        return self.grid[0] * self.grid[1] * self.grid[2]

    @property
    def basis(self) -> int:
        return len(self.basis_offsets)

    @property
    def atom_count(self) -> int:
        return self.active_cells * self.basis

    def primitive(self, axis: int) -> Vector:
        return tuple(self.primitive_vectors[3 * axis + component] for component in range(3))  # type: ignore[return-value]

    def position(self, cell: tuple[int, int, int], basis: int) -> Vector:
        result = self.basis_offsets[basis]
        for axis, coordinate in enumerate(cell):
            result = _add(result, _scale(self.primitive(axis), float(coordinate)))
        return result

    def validate(self) -> None:
        if len(self.primitive_vectors) != 9:
            raise ValueError(f"{self.name}: primitive_vectors must have 9 values")
        if any(d <= 0 for d in self.grid):
            raise ValueError(f"{self.name}: grid dimensions must be positive")
        if not self.basis_offsets:
            raise ValueError(f"{self.name}: at least one basis offset is required")
        if len(self.moments_mu_b) == 0:
            raise ValueError(f"{self.name}: at least one ensemble is required")
        expected = (self.active_cells, self.basis)
        for ensemble in self.moments_mu_b:
            if (len(ensemble), len(ensemble[0]) if ensemble else 0) != expected:
                raise ValueError(f"{self.name}: moments do not have cell/basis shape {expected}")
        det = (
            self.primitive_vectors[0] * (self.primitive_vectors[4] * self.primitive_vectors[8] - self.primitive_vectors[7] * self.primitive_vectors[5])
            - self.primitive_vectors[3] * (self.primitive_vectors[1] * self.primitive_vectors[8] - self.primitive_vectors[7] * self.primitive_vectors[2])
            + self.primitive_vectors[6] * (self.primitive_vectors[1] * self.primitive_vectors[5] - self.primitive_vectors[4] * self.primitive_vectors[2])
        )
        if not math.isfinite(det) or det <= 0.0:
            raise ValueError(f"{self.name}: primitive-vector determinant must be positive, got {det}")
        if any(not math.isfinite(value) for value in self.primitive_vectors):
            raise ValueError(f"{self.name}: primitive vectors contain a non-finite value")
        if any(not math.isfinite(value) for offset in self.basis_offsets for value in offset):
            raise ValueError(f"{self.name}: basis offsets contain a non-finite value")


def dimensionless_fields(case: FiniteCase) -> tuple[tuple[tuple[Vector, ...], ...], ...]:
    """Evaluate each ordered target/source pair once, excluding exact self."""
    case.validate()
    fields = [
        [[(0.0, 0.0, 0.0) for _ in range(case.basis)] for _ in range(case.active_cells)]
        for _ in case.moments_mu_b
    ]
    positions = [
        [case.position(cell, basis) for basis in range(case.basis)]
        for cell in cell_coordinates(case.grid)
    ]
    for ensemble, moments in enumerate(case.moments_mu_b):
        for target_cell in range(case.active_cells):
            for target_basis in range(case.basis):
                target_position = positions[target_cell][target_basis]
                for source_cell in range(case.active_cells):
                    for source_basis in range(case.basis):
                        if target_cell == source_cell and target_basis == source_basis:
                            continue
                        displacement = _sub(target_position, positions[source_cell][source_basis])
                        r2 = sum(component * component for component in displacement)
                        if r2 == 0.0:
                            raise ValueError(
                                f"{case.name}: nonself source/target overlap at "
                                f"target ({target_cell},{target_basis}), source ({source_cell},{source_basis})")
                        contribution = tensor_vector(point_tensor(displacement), moments[source_cell][source_basis])
                        fields[ensemble][target_cell][target_basis] = _add(
                            fields[ensemble][target_cell][target_basis], contribution)
    return tuple(tuple(tuple(vector for vector in cell) for cell in ensemble) for ensemble in fields)


def apply_units(case: FiniteCase, fields_dimensionless: Sequence[Sequence[Sequence[Vector]]]) -> tuple[tuple[tuple[Vector, ...], ...], ...]:
    scale = physical_prefactor(case.alat_m)
    return tuple(tuple(tuple(_scale(vector, scale) for vector in cell) for cell in ensemble)
                 for ensemble in fields_dimensionless)


def total_energy_j(case: FiniteCase, fields_dimensionless: Sequence[Sequence[Sequence[Vector]]]) -> tuple[float, ...]:
    """Return ``-1/2 mu_B * prefactor * sum(m dot B0)`` per ensemble."""
    scale = physical_prefactor(case.alat_m)
    return tuple(
        -0.5 * MU_B_J_PER_T * scale * sum(
            sum(moments[cell][basis][component] * fields_dimensionless[ensemble][cell][basis][component]
                for component in range(3))
            for cell in range(case.active_cells) for basis in range(case.basis)
        )
        for ensemble, moments in enumerate(case.moments_mu_b)
    )


def evaluate(case: FiniteCase) -> dict:
    dim = dimensionless_fields(case)
    fields_t = apply_units(case, dim)
    energies = total_energy_j(case, dim)
    return {
        "fields_dimensionless": dim,
        "fields_t": fields_t,
        "total_energy_j": energies,
        "energy_per_atom_mry": tuple(value / case.atom_count / MRY_J for value in energies),
    }


def kernel_batch(row: int, column: int, target_basis: int, source_basis: int, basis: int) -> int:
    return row + 3 * (column + 3 * (target_basis + basis * source_basis))


def field_batch(component: int, basis_channel: int, basis: int, ensemble: int) -> int:
    return component + 3 * (basis_channel + basis * ensemble)


def fft_cell(grid: tuple[int, int, int], coordinate: tuple[int, int, int]) -> int:
    return coordinate[0] + grid[0] * (coordinate[1] + grid[1] * coordinate[2])


def embed_kernel(case: FiniteCase, fft_grid: tuple[int, int, int] | None = None) -> tuple[float, ...]:
    """Embed the dimensionless finite kernel in zero-filled batch-major storage."""
    case.validate()
    padded = fft_grid or tuple(2 * size - 1 for size in case.grid)
    if any(padded[axis] < 2 * case.grid[axis] - 1 for axis in range(3)):
        raise ValueError(f"{case.name}: fft_grid is smaller than 2*active_grid-1")
    fft_cells = padded[0] * padded[1] * padded[2]
    result = [0.0] * (fft_cells * 9 * case.basis * case.basis)
    for dx in range(-(case.grid[0] - 1), case.grid[0]):
        for dy in range(-(case.grid[1] - 1), case.grid[1]):
            for dz in range(-(case.grid[2] - 1), case.grid[2]):
                displacement = (dx, dy, dz)
                q = tuple(value if value >= 0 else padded[axis] + value for axis, value in enumerate(displacement))
                qcell = fft_cell(padded, q)
                for target_basis in range(case.basis):
                    for source_basis in range(case.basis):
                        if displacement == (0, 0, 0) and target_basis == source_basis:
                            tensor: Tensor = ((0.0, 0.0, 0.0),) * 3
                        else:
                            r = _add(
                                _add(_scale(case.primitive(0), float(dx)),
                                     _add(_scale(case.primitive(1), float(dy)), _scale(case.primitive(2), float(dz)))),
                                _sub(case.basis_offsets[target_basis], case.basis_offsets[source_basis]),
                            )
                            r2 = sum(value * value for value in r)
                            if r2 == 0.0:
                                raise ValueError(f"{case.name}: nonself embedded displacement overlaps")
                            tensor = point_tensor(r)
                        for row in range(3):
                            for column in range(3):
                                batch = kernel_batch(row, column, target_basis, source_basis, case.basis)
                                result[qcell + fft_cells * batch] = tensor[row][column]
    return tuple(result)


def embed_active_moments(case: FiniteCase, ensemble: int, fft_grid: tuple[int, int, int]) -> tuple[float, ...]:
    """Pack exactly Nactive source cells into an otherwise zero FFT buffer."""
    if not 0 <= ensemble < len(case.moments_mu_b):
        raise IndexError("ensemble out of range")
    fft_cells = fft_grid[0] * fft_grid[1] * fft_grid[2]
    result = [0.0] * (fft_cells * 3 * case.basis)
    for x, y, z in cell_coordinates(case.grid):
        active = cell_index(case.grid, x, y, z)
        padded = fft_cell(fft_grid, (x, y, z))
        for basis in range(case.basis):
            for component in range(3):
                result[padded + fft_cells * field_batch(component, basis, case.basis, 0)] = \
                    case.moments_mu_b[ensemble][active][basis][component]
    return tuple(result)


def padded_direct_convolution(case: FiniteCase, ensemble: int,
                              fft_grid: tuple[int, int, int] | None = None) -> tuple[tuple[Vector, ...], ...]:
    """Apply the embedded kernel through padded cyclic indexing.

    This is a host test model for the FFT contract, not a second physical
    evaluator.  Its result must agree with ``dimensionless_fields`` while it
    deliberately retains the active/padded storage split.
    """
    padded = fft_grid or tuple(2 * size - 1 for size in case.grid)
    kernel = embed_kernel(case, padded)
    fft_cells = padded[0] * padded[1] * padded[2]
    fields = [[(0.0, 0.0, 0.0) for _ in range(case.basis)] for _ in range(case.active_cells)]
    moments = case.moments_mu_b[ensemble]
    for tx, ty, tz in cell_coordinates(case.grid):
        target_cell = cell_index(case.grid, tx, ty, tz)
        for sx, sy, sz in cell_coordinates(case.grid):
            source_cell = cell_index(case.grid, sx, sy, sz)
            displacement = (tx - sx, ty - sy, tz - sz)
            q = tuple(value % padded[axis] for axis, value in enumerate(displacement))
            qcell = fft_cell(padded, q)
            for target_basis in range(case.basis):
                for source_basis in range(case.basis):
                    for row in range(3):
                        value = sum(
                            kernel[qcell + fft_cells * kernel_batch(row, column, target_basis, source_basis, case.basis)]
                            * moments[source_cell][source_basis][column]
                            for column in range(3)
                        )
                        updated = list(fields[target_cell][target_basis])
                        updated[row] += value
                        fields[target_cell][target_basis] = tuple(updated)  # type: ignore[assignment]
    return tuple(tuple(vector for vector in cell) for cell in fields)


def _nested(values: object) -> object:
    if isinstance(values, tuple):
        return [_nested(value) for value in values]
    return values


def _moments(*ensembles: Sequence[Sequence[Vector]]) -> tuple[tuple[tuple[Vector, ...], ...], ...]:
    return tuple(tuple(tuple(tuple(vector) for vector in cell) for cell in ensemble) for ensemble in ensembles)


def _uniform(grid: tuple[int, int, int], vector: Vector, basis: int = 1) -> tuple[tuple[tuple[Vector, ...], ...], ...]:
    return _moments([[vector for _ in range(basis)] for _ in range(grid[0] * grid[1] * grid[2])])


def deterministic_cases() -> tuple[FiniteCase, ...]:
    """Cases are intentionally geometrically and algebraically non-redundant."""
    cubic = (1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0)
    cases: list[FiniteCase] = []
    cases.append(FiniteCase("single_1x1x1", (1, 1, 1), cubic, ((0.0, 0.0, 0.0),), _uniform((1, 1, 1), (0.7, -0.2, 0.4))))
    cases.append(FiniteCase("pair_2x1x1_axial", (2, 1, 1), cubic, ((0.0, 0.0, 0.0),), _uniform((2, 1, 1), (1.0, 0.0, 0.0))))
    cases.append(FiniteCase("pair_2x1x1_transverse", (2, 1, 1), cubic, ((0.0, 0.0, 0.0),), _uniform((2, 1, 1), (0.0, 0.0, 1.0))))
    moments_2x3 = tuple(tuple(((0.17 + 0.11 * cell, -0.23 + 0.07 * cell, 0.31 - 0.05 * cell),),) for cell in range(6))
    cases.append(FiniteCase("nonuniform_2x3x1", (2, 3, 1), cubic, ((0.0, 0.0, 0.0),), (moments_2x3,)))
    moments_g3 = tuple(tuple(((0.21 + 0.03 * cell, -0.19 + 0.02 * cell, 0.37 - 0.04 * cell),),) for cell in range(6))
    cases.append(FiniteCase("g3_gt_g2_2x1x3", (2, 1, 3), cubic, ((0.0, 0.0, 0.0),), (moments_g3,)))
    skew = (1.0, 0.20, 0.0, 0.10, 1.10, 0.05, 0.15, 0.08, 1.30)
    skew_moments = tuple(tuple(((0.13 + 0.08 * cell, -0.29 + 0.04 * cell, 0.23 - 0.03 * cell),),) for cell in range(4))
    cases.append(FiniteCase("skew_2x2x1", (2, 2, 1), skew, ((0.0, 0.0, 0.0),), (skew_moments,)))
    basis_offsets = ((0.0, 0.0, 0.0), (0.27, 0.13, 0.19))
    basis_moments = _moments(
        [[(0.2, -0.3, 0.4), (-0.5, 0.7, 0.1)], [(0.9, 0.2, -0.6), (0.3, -0.8, 0.5)]])
    cases.append(FiniteCase("basis_2x1x1_nonuniform", (2, 1, 1), cubic, basis_offsets, basis_moments))
    m4 = _moments(
        [[(0.11, -0.22, 0.33)] for _ in range(4)],
        [[(-0.41, 0.52, -0.63)] for _ in range(4)],
        [[(0.71, 0.12, -0.23)] for _ in range(4)],
        [[(-0.31, -0.42, 0.83)] for _ in range(4)])
    cases.append(FiniteCase("four_distinct_ensembles", (2, 2, 1), cubic, ((0.0, 0.0, 0.0),), m4))
    cases.append(FiniteCase("thin_film_in_plane", (4, 3, 1), cubic, ((0.0, 0.0, 0.0),), _uniform((4, 3, 1), (1.0, 0.0, 0.0))))
    cases.append(FiniteCase("thin_film_normal", (4, 3, 1), cubic, ((0.0, 0.0, 0.0),), _uniform((4, 3, 1), (0.0, 0.0, 1.0))))
    return tuple(cases)


def golden_document() -> dict:
    result = {"oracle": "independent finite open-boundary point-dipole sum", "constants": {
        "mu_b_j_per_t": MU_B_J_PER_T, "mu0_over_4pi": MU0_OVER_4PI_T_M_PER_A, "mry_j": MRY_J}}
    result["cases"] = []
    for case in deterministic_cases():
        evaluated = evaluate(case)
        result["cases"].append({"name": case.name, "grid": list(case.grid),
                                 "primitive_vectors": list(case.primitive_vectors),
                                 "basis_offsets": [list(offset) for offset in case.basis_offsets],
                                 "alat_m": case.alat_m, "moments_mu_b": _nested(case.moments_mu_b),
                                 **{key: _nested(value) for key, value in evaluated.items()}})
    return result


def main() -> int:
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write-goldens", type=Path)
    args = parser.parse_args()
    if not args.write_goldens:
        parser.error("--write-goldens is required")
    args.write_goldens.write_text(json.dumps(golden_document(), indent=2) + "\n")
    print(f"wrote {args.write_goldens}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
