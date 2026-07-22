#!/usr/bin/env python3
"""Independent fp64 3D-periodic Ewald reference for point magnetic dipoles.

This deliberately small O(N^2 (R + K)) evaluator is a validation oracle for
the block-one macrocell limit.  It uses no UppASD routines and no GPU code.
Coordinates and cell vectors are expressed in one common length unit; the
returned field has the corresponding dipole-tensor prefactor omitted.  Multiply
the result by mu0*muB/(4*pi*alat**3) when moments use UppASD lattice units.

The convention is 3D periodic Ewald with conducting (tin-foil) exterior:
k=0 is omitted.  It implements B = Hessian(1/r) m and
E = -1/2 sum_i m_i dot B_i.
"""
from __future__ import annotations

from dataclasses import dataclass
from itertools import product
import math
from typing import Sequence


Vec = tuple[float, float, float]
Mat = tuple[Vec, Vec, Vec]  # columns are the three cell vectors


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


def cell_vector(cell: Mat, n: tuple[int, int, int]) -> Vec:
    return add(add(scale(cell[0], n[0]), scale(cell[1], n[1])), scale(cell[2], n[2]))


def reciprocal_cell(cell: Mat) -> tuple[Mat, float]:
    """Return 2*pi H^-T as column vectors, and positive cell volume."""
    volume = dot(cell[0], cross(cell[1], cell[2]))
    if volume <= 0.0:
        raise ValueError("cell must be right-handed with positive volume")
    factor = 2.0 * math.pi / volume
    return (scale(cross(cell[1], cell[2]), factor),
            scale(cross(cell[2], cell[0]), factor),
            scale(cross(cell[0], cell[1]), factor)), volume


def real_tensor(r: Vec, alpha: float) -> tuple[float, float, float, float, float, float]:
    """Symmetric Hessian of erfc(alpha*r)/r, packed xx,yy,zz,xy,xz,yz."""
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
    b = -erfc / (radius**3) - 2.0 * alpha * gaussian / (math.sqrt(math.pi) * radius**2)
    x, y, z = r
    return (a * x * x + b, a * y * y + b, a * z * z + b,
            a * x * y, a * x * z, a * y * z)


def tensor_times(tensor: tuple[float, float, float, float, float, float], moment: Vec) -> Vec:
    xx, yy, zz, xy, xz, yz = tensor
    x, y, z = moment
    return (xx * x + xy * y + xz * z,
            xy * x + yy * y + yz * z,
            xz * x + yz * y + zz * z)


@dataclass(frozen=True)
class EwaldResult:
    fields: tuple[Vec, ...]
    energy: float


def evaluate(positions: Sequence[Vec], moments: Sequence[Vec], cell: Mat, *, alpha: float,
             real_images: tuple[int, int, int], reciprocal_images: tuple[int, int, int]) -> EwaldResult:
    """Evaluate a truncated but absolutely convergent 3D tin-foil Ewald sum.

    Increase both image ranges until the result is stable.  The rectangular
    image ranges are intentionally explicit: they make convergence studies
    transparent for skew cells rather than hiding an implementation choice.
    """
    if len(positions) != len(moments) or not positions:
        raise ValueError("positions and moments must be non-empty and have equal length")
    if alpha <= 0.0 or min(*real_images, *reciprocal_images) < 0:
        raise ValueError("alpha must be positive and image ranges non-negative")
    reciprocal, volume = reciprocal_cell(cell)
    fields = [(0.0, 0.0, 0.0) for _ in positions]

    for target, target_pos in enumerate(positions):
        field = (0.0, 0.0, 0.0)
        for source, source_pos in enumerate(positions):
            for image in product(*(range(-extent, extent + 1) for extent in real_images)):
                if target == source and image == (0, 0, 0):
                    continue
                r = sub(target_pos, add(source_pos, cell_vector(cell, image)))
                field = add(field, tensor_times(real_tensor(r, alpha), moments[source]))
        fields[target] = field

    for wave_index in product(*(range(-extent, extent + 1) for extent in reciprocal_images)):
        if wave_index == (0, 0, 0):
            continue  # conducting exterior convention
        wave = cell_vector(reciprocal, wave_index)
        wave2 = dot(wave, wave)
        prefactor = -4.0 * math.pi * math.exp(-wave2 / (4.0 * alpha * alpha)) / (volume * wave2)
        for target, target_pos in enumerate(positions):
            field = fields[target]
            for source, source_pos in enumerate(positions):
                phase = math.cos(dot(wave, sub(target_pos, source_pos)))
                field = add(field, scale(wave, prefactor * phase * dot(wave, moments[source])))
            fields[target] = field

    # Subtracting the smooth self tensor Hessian(erf(alpha*r)/r)|r=0 gives
    # +4 alpha^3/(3 sqrt(pi)) m in the B=Hessian(1/r)m convention.
    self_prefactor = 4.0 * alpha**3 / (3.0 * math.sqrt(math.pi))
    fields = tuple(add(field, scale(moment, self_prefactor)) for field, moment in zip(fields, moments))
    energy = -0.5 * sum(dot(moment, field) for moment, field in zip(moments, fields))
    return EwaldResult(fields=fields, energy=energy)
