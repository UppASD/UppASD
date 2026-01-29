"""
neighbors.py

NumPy-only neighbor and shell construction utilities.

Responsibilities:
- Minimum-image displacement vectors
- Brute-force neighbor list construction
- Deterministic distance-based shell assignment

NO physics logic.
NO UppASD coupling beyond indexing conventions.
"""

import numpy as np
from typing import Tuple, List


# ======================================================================
# Geometry helpers
# ======================================================================


def minimum_image_vector(
    dr: np.ndarray,
    cell: np.ndarray,
    inv_cell: np.ndarray,
    pbc: Tuple[bool, bool, bool],
) -> np.ndarray:
    """
    Apply minimum-image convention to a displacement vector.

    Parameters
    ----------
    dr : ndarray (3,)
        Cartesian displacement r_j - r_i
    cell : ndarray (3,3)
        Lattice vectors (rows)
    inv_cell : ndarray (3,3)
        Inverse lattice matrix
    pbc : tuple(bool, bool, bool)
        Periodic boundary flags
    """
    # fractional displacement
    df = dr @ inv_cell

    for k in range(3):
        if pbc[k]:
            df[k] -= np.round(df[k])

    return df @ cell


# ======================================================================
# Neighbor construction
# ======================================================================


def build_neighbors(
    positions: np.ndarray,
    cell: np.ndarray,
    cutoff: float,
    pbc: Tuple[bool, bool, bool] = (True, True, True),
):
    """
    Brute-force neighbor list using minimum-image convention.

    Parameters
    ----------
    positions : (N,3) ndarray
    cell : (3,3) ndarray
    cutoff : float
    pbc : tuple(bool, bool, bool)

    Returns
    -------
    pairs : list of (i+1, j+1)
        1-based atom index pairs
    vectors : list of ndarray
        Displacement vectors R_ij
    distances : list of float
        |R_ij|
    """
    positions = np.asarray(positions, dtype=float)
    cell = np.asarray(cell, dtype=float)
    inv_cell = np.linalg.inv(cell)

    natom = len(positions)

    pairs: List[Tuple[int, int]] = []
    vectors: List[np.ndarray] = []
    distances: List[float] = []

    for i in range(natom):
        for j in range(natom):
            if i == j:
                continue

            dr = positions[j] - positions[i]
            R = minimum_image_vector(dr, cell, inv_cell, pbc)
            r = np.linalg.norm(R)

            if r <= cutoff and r > 0.0:
                pairs.append((i + 1, j + 1))  # UppASD is 1-based
                vectors.append(R)
                distances.append(r)

    return pairs, vectors, distances


# ======================================================================
# Shell assignment
# ======================================================================


def assign_shells(
    distances: List[float],
    tol: float = 1e-5,
):
    """
    Assign deterministic shell indices based on distances.

    Distances within `tol` are considered equivalent.

    Parameters
    ----------
    distances : list of float
    tol : float

    Returns
    -------
    shell_indices : list of int
        Shell index for each distance (1-based)
    shell_radii : dict[int, float]
        Representative distance for each shell
    """
    distances = np.asarray(distances, dtype=float)

    # Sort unique shell radii deterministically
    unique = []
    for r in np.sort(distances):
        if not any(abs(r - r0) < tol for r0 in unique):
            unique.append(r)

    shell_radii = {i + 1: r for i, r in enumerate(unique)}

    shell_indices = []
    for r in distances:
        for i, r0 in shell_radii.items():
            if abs(r - r0) < tol:
                shell_indices.append(i)
                break

    return shell_indices, shell_radii
