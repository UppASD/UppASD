"""
neighbors.py

NumPy-only neighbor and shell construction utilities.

Responsibilities:
- Minimum-image displacement vectors
- Expansion-based neighbor list construction
- Deterministic distance-based shell assignment
"""

from typing import Dict, List, Tuple

import numpy as np

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
    df = dr @ inv_cell

    for k in range(3):
        if pbc[k]:
            df[k] -= np.round(df[k])

    return df @ cell


# ======================================================================
# Neighbor construction
# ======================================================================


def build_neighbors(
    basis_positions: np.ndarray,
    unit_cell: np.ndarray,
    cutoff: float,
    pbc: Tuple[bool, bool, bool] = (True, True, True),
) -> Tuple[List[Tuple[int, int]], List[np.ndarray], List[float]]:
    """
    Build neighbor list by expanding the unit cell to cover the cutoff radius.

    This ensures that for any atom i in the basis, all neighbors j (including
    periodic images) within 'cutoff' are found.

    Parameters
    ----------
    basis_positions : (N,3) ndarray
        Positions of atoms in the primary unit cell.
    unit_cell : (3,3) ndarray
        Unit cell lattice vectors (rows).
    cutoff : float
        Maximum search distance.
    pbc : tuple(bool, bool, bool)
        Periodic boundary flags.

    Returns
    -------
    pairs : list of (i+1, j+1)
        Directed bonds originating from the unit cell basis.
    vectors : list of ndarray
        Displacement vectors R_ij.
    distances : list of float
        |R_ij|.
    """
    basis_positions = np.asarray(basis_positions, dtype=float)
    unit_cell = np.asarray(unit_cell, dtype=float)
    inv_cell = np.linalg.inv(unit_cell)
    natom = len(basis_positions)

    # Determine how many unit cells we need to search in each direction
    # to guarantee we cover the cutoff radius.
    cell_norms = np.linalg.norm(unit_cell, axis=1)
    # n_expand determines the range [-n, n] for periodic images
    n_expand = np.ceil(cutoff / cell_norms).astype(int)

    pairs, vectors, distances = [], [], []

    for i in range(natom):
        p_i = basis_positions[i]

        # Search against all atoms in the basis...
        for j in range(natom):
            p_j_base = basis_positions[j]

            # ...and all relevant periodic images of those atoms
            # defined by the expansion range.
            for nx in range(-n_expand[0], n_expand[0] + 1) if pbc[0] else [0]:
                for ny in range(-n_expand[1], n_expand[1] + 1) if pbc[1] else [0]:
                    for nz in range(-n_expand[2], n_expand[2] + 1) if pbc[2] else [0]:

                        # Calculate the vector to the specific periodic image
                        offset = (
                            nx * unit_cell[0] + ny * unit_cell[1] + nz * unit_cell[2]
                        )
                        R = (p_j_base + offset) - p_i
                        r = np.linalg.norm(R)

                        # Distance guard: Skip self-interaction (r=0) and respect cutoff
                        if 1e-8 < r <= cutoff:
                            pairs.append((i + 1, j + 1))
                            vectors.append(R)
                            distances.append(r)

    return pairs, vectors, distances


# ======================================================================
# Shell assignment and Reporting
# ======================================================================


def assign_shells(
    distances: List[float],
    tol: float = 1e-5,
) -> Tuple[List[int], Dict[int, float]]:
    """
    Assign deterministic shell indices based on distances.
    """
    distances = np.asarray(distances, dtype=float)
    if distances.size == 0:
        return [], {}

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


def get_shell_report(
    pairs: List[Tuple[int, int]],
    vectors: List[np.ndarray],
    distances: List[float],
    tol: float = 1e-5,
) -> str:
    """
    Generate a formatted report of neighbor shells for the provided basis.
    """
    if not distances:
        return "No neighbors found."

    shell_indices, shell_radii = assign_shells(distances, tol=tol)

    lines = ["=== Shell Lookup Table ==="]
    lines.append(f"{'Index':<6} | {'Radius':<12}")
    lines.append("-" * 22)
    for idx, radius in sorted(shell_radii.items()):
        lines.append(f"{idx:<6} | {radius:<12.6f}")

    lines.append("\n=== Basis Neighbor Report ===")
    lines.append(
        f"{'Source':<8} | {'Target':<8} | {'Shell':<6} | {'Distance':<10} | {'Vector (r)'}"
    )
    lines.append("-" * 75)

    # Identify unique source atoms from the pairs
    sources = sorted(list(set(p[0] for p in pairs)))

    for src in sources:
        src_bonds = []
        for idx, (i, j) in enumerate(pairs):
            if i == src:
                src_bonds.append(
                    {
                        "target": j,
                        "shell": shell_indices[idx],
                        "dist": distances[idx],
                        "vec": vectors[idx],
                    }
                )

        src_bonds.sort(key=lambda x: (x["shell"], x["dist"]))
        for b in src_bonds:
            v = b["vec"]
            lines.append(
                f"Atom {src:<3}  | Atom {b['target']:<3}  | {b['shell']:<6} | "
                f"{b['dist']:<10.4f} | [{v[0]:>6.3f}, {v[1]:>6.3f}, {v[2]:>6.3f}]"
            )
        lines.append("-" * 75)

    return "\n".join(lines)
