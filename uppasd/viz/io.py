"""
viz/io.py

Data extraction and reshaping helpers for UppASD visualizations.

Design principles
-----------------
- NO plotting
- NO UppASD file format logic
- Operates only on ASDResults and NumPy arrays
- Explicit is better than implicit
"""

import numpy as np
from typing import Tuple, Optional


# ======================================================================
# Scalar time series
# ======================================================================

def timeseries(
    results,
    key: str,
    column: Optional[str] = None,
    step_key: str = "iter",
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Extract a scalar time series from ASDResults.

    Parameters
    ----------
    results : ASDResults
        Parsed UppASD results
    key : str
        Observable name (e.g. 'totenergy', 'averages', 'cumulants')
    column : str, optional
        Column to extract (required if multiple data columns exist)
    step_key : str
        Name of iteration column (default: 'iter')

    Returns
    -------
    step : ndarray (Nt,)
    values : ndarray (Nt,)
    """
    table = results.get(key)
    if table is None:
        raise KeyError(f"No observable '{key}' in results")

    if not isinstance(table, dict):
        raise TypeError(
            f"Observable '{key}' is not a column table "
            f"(got {type(table)})"
        )

    if step_key not in table:
        raise KeyError(
            f"Observable '{key}' has no step column '{step_key}'. "
            f"Available columns: {tuple(table.keys())}"
        )

    step = table[step_key]

    data_columns = [k for k in table.keys() if k != step_key]

    if column is None:
        if len(data_columns) != 1:
            raise ValueError(
                f"Observable '{key}' has multiple columns {data_columns}. "
                "Please specify column explicitly."
            )
        column = data_columns[0]

    if column not in table:
        raise KeyError(
            f"Column '{column}' not found in '{key}'. "
            f"Available: {tuple(table.keys())}"
        )

    return step, table[column]


# ======================================================================
# Local spin snapshots
# ======================================================================
def spin_snapshots(results, key="moment"):
    """
    Extract local spin snapshots.

    Parameters
    ----------
    results : ASDResults
    key : str
        'moment' or 'restart'

    Returns
    -------
    steps : ndarray (Nt,)
    spins : ndarray (Nt, Natoms, 3)
    """
    table = results[key]

    iters = table["iter"]
    atoms = table["iatom"]
    mx = table["mx"]
    my = table["my"]
    mz = table["mz"]

    # --- infer dimensions ---
    natom = atoms.max()
    steps = np.unique(iters)
    nstep = len(steps)

    spins = np.zeros((nstep, natom, 3))

    # --- fill snapshots ---
    for t, step in enumerate(steps):
        mask = iters == step
        idx = atoms[mask] - 1  # Fortran → Python index

        spins[t, idx, 0] = mx[mask]
        spins[t, idx, 1] = my[mask]
        spins[t, idx, 2] = mz[mask]

    return steps, spins


# ======================================================================
# Effective magnetic field snapshots
# ======================================================================
def effective_field_snapshots(
    results,
    key: str = "befftot",
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Extract time-resolved effective magnetic field snapshots.

    Returns
    -------
    steps : ndarray (Nt,)
    beff  : ndarray (Nt, Natom, 3)
    """
    table = results[key] if key in results.available else None
    if table is None:
        raise KeyError(f"No observable '{key}' in results")

    # --- iteration index ---
    iters = table["iter"]

    # --- atom index (UppASD naming variants) ---
    if "site" in table:
        atoms = table["site"]
    elif "iatom" in table:
        atoms = table["iatom"]
    elif "atom" in table:
        atoms = table["atom"]
    else:
        raise KeyError(
            "No atom index column found in befftot. "
            "Expected one of: 'site', 'iatom', 'atom'. "
            f"Available columns: {tuple(table.keys())}"
        )

    # --- vector components ---
    bx = table["bx"]
    by = table["by"]
    bz = table["bz"]

    steps = np.unique(iters)
    natom = int(atoms.max())

    beff = np.zeros((len(steps), natom, 3), dtype=float)

    for ti, step in enumerate(steps):
        mask = (iters == step)
        idx = atoms[mask].astype(int) - 1  # 1-based → 0-based

        beff[ti, idx, 0] = bx[mask]
        beff[ti, idx, 1] = by[mask]
        beff[ti, idx, 2] = bz[mask]

    return steps, beff


def site_snapshots(
    results,
    key: str,
    value_column: str | None = None,
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Extract site-resolved scalar snapshots.

    Typical use cases:
    - skyrmion density (dens_skynum)
    - local charges
    - site-resolved observables without vector components

    Expected table columns (after output.py normalization):
        iter   : timestep
        site   : atom index (1-based)
        value  : scalar observable (name varies)

    Parameters
    ----------
    results : ASDResults
    key : str
        Observable key (e.g. "dens_skynum")
    value_column : str, optional
        Name of scalar column. If None, auto-detects the only
        non-index column.

    Returns
    -------
    steps : ndarray, shape (Nt,)
    values : ndarray, shape (Nt, Natom)
    """
    table = results[key]
    if table is None:
        raise KeyError(f"No observable '{key}' in results")

    # Mandatory columns
    iters = table["iter"]
    sites = table["site"]  # 1-based atom index

    # Detect scalar column
    if value_column is None:
        candidates = [
            k for k in table.keys()
            if k not in ("iter", "site", "ens", "raw")
        ]
        if len(candidates) != 1:
            raise ValueError(
                f"Cannot auto-detect scalar column for '{key}'. "
                f"Candidates: {candidates}"
            )
        value_column = candidates[0]

    values_flat = table[value_column]

    # Identify time steps and atom count
    steps = np.unique(iters)
    natom = sites.max()

    Nt = len(steps)
    values = np.zeros((Nt, natom), dtype=float)

    step_index = {step: i for i, step in enumerate(steps)}

    for it, site, val in zip(iters, sites, values_flat):
        i = step_index[it]
        j = site - 1  # 1-based → 0-based
        values[i, j] = val

    return steps, values

# ======================================================================
# Supercell reshaping helpers
# ======================================================================

def reshape_supercell(
    data: np.ndarray,
    ncell: Tuple[int, int, int],
    reduce_z: bool = False,
):
    """
    Reshape flat supercell data into (N1, N2, N3) or (N1, N2).

    Assumes:
    - one atom per unit cell
    - UppASD ordering consistent with ncell

    Parameters
    ----------
    data : ndarray (Natom,)
    ncell : (N1, N2, N3)
    reduce_z : bool
        If True, average over z-direction

    Returns
    -------
    ndarray
        (N1, N2, N3) or (N1, N2)
    """
    n1, n2, n3 = ncell
    arr = data.reshape((n1, n2, n3))

    if reduce_z:
        return arr.mean(axis=2)

    return arr


# ======================================================================
# Space–time maps (1D systems)
# ======================================================================

def spin_time_map(
    results,
    component: int = 2,
    key: str = "moment",
) -> np.ndarray:
    """
    Build a space–time map m(x,t) for 1D systems.

    Parameters
    ----------
    results : ASDResults
    component : int
        0=x, 1=y, 2=z
    key : str
        Observable ('moment')

    Returns
    -------
    xt : ndarray (Nt, Natom)
        Each row corresponds to a time snapshot
    """
    steps, spins = spin_snapshots(results, key=key)
    return spins[..., component]

import numpy as np
from typing import Tuple


def skyrmion_com(
    results,
    key: str = "cmass_skynum",
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Extract skyrmion center-of-mass trajectory.

    Parameters
    ----------
    results : ASDResults
    key : str
        Observable name (default: 'cmass_skynum')

    Returns
    -------
    steps : ndarray, shape (Nt,)
        Iteration indices
    rsk : ndarray, shape (Nt, 3)
        Skyrmion COM positions (x, y, z)
    """
    table = results[key]

    steps = table["iter"]
    rx = table["rsk_x"]
    ry = table["rsk_y"]
    rz = table["rsk_z"]

    rsk = np.stack([rx, ry, rz], axis=-1)
    return steps, rsk
