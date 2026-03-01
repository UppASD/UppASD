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
from typing import Tuple, Optional, Dict


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


def snapshot_to_grid(
    vals: np.ndarray,
    coord: Optional[Dict] = None,
    ncell: Optional[Tuple[int, int, int]] = None,
    reduce_z: bool = True,
):
    """
    Convert a per-atom 1D array of values into a 2D grid suitable for heatmap plotting.

    Parameters
    ----------
    vals : ndarray (Natom,)
        Per-atom scalar values (e.g. mz) indexed by atom number 1..Natom.
    coord : dict, optional
        Coordinate table as returned by `read_coord` (keys: 'iatom','x','y').
    ncell : tuple (Nx,Ny,Nz), optional
        Supercell shape to reshape into. If provided, `vals` will be reshaped
        into (Nx,Ny,Nz) and reduced over z if `reduce_z` is True.
    reduce_z : bool
        If True and Nz>1, average over z-direction.

    Returns
    -------
    grid2d : ndarray (Ny, Nx)
        2D grid ready for `imshow(origin='lower')`.
    """
    vals = np.asarray(vals)

    if ncell is not None:
        Nx, Ny, Nz = ncell
        if vals.size != Nx * Ny * Nz:
            raise ValueError("vals length does not match ncell product")
        arr = vals.reshape((Nx, Ny, Nz))
        if reduce_z and Nz > 1:
            return arr.mean(axis=2)
        if Nz > 1:
            return arr
        return arr.reshape((Nx, Ny))

    if coord is None:
        raise ValueError("coord must be provided when ncell is not given")

    xs = coord["x"]
    ys = coord["y"]
    iatom = coord["iatom"].astype(int)

    ux = np.unique(xs)
    uy = np.unique(ys)
    nx = ux.size
    ny = uy.size

    grid = np.full((ny, nx), np.nan, dtype=float)

    # Build index maps using exact matches
    x_to_ix = {float(v): i for i, v in enumerate(ux)}
    y_to_iy = {float(v): i for i, v in enumerate(uy)}

    for ia, x, y in zip(iatom, xs, ys):
        ix = x_to_ix.get(float(x))
        iy = y_to_iy.get(float(y))
        if ix is None or iy is None:
            continue
        grid[iy, ix] = vals[ia - 1]

    return grid


def results_snapshot_grid(
    results,
    key: str = "restart",
    component: int = 2,
    snapshot_index: int = -1,
    idx: int | None = None,
    ncell: Optional[Tuple[int, int, int]] = None,
    reduce_z: bool = True,
):
    """
    Convenience wrapper: extract a single snapshot from `results` and return
    (iter_id, grid2d) ready for plotting.

    Parameters
    ----------
    results : ASDResults
    key : str
        'restart' or 'moment' (or any table with 'iter','iatom','mx','my','mz').
    component : int
        0=x,1=y,2=z component to extract.
    snapshot_index : int
        Index into time snapshots for 'moment' (supports negative indexing). For
        'restart' this is ignored.
    ncell : optional
        Supercell shape to reshape into.
    reduce_z : bool
        Whether to average over z when ncell indicates Nz>1.

    Returns
    -------
    (iter_id, grid2d)
    or None if data missing.
    """
    if key == "moment":
        # use spin_snapshots to get structured array
        steps, spins = spin_snapshots(results, key="moment")
        if len(steps) == 0:
            return None
        # accept `idx` as an alias for `snapshot_index` (notebook-friendly)
        if idx is not None:
            snapshot_index = idx
        sel_idx = snapshot_index
        iter_id = int(steps[sel_idx])
        vals = spins[sel_idx, :, component]

        coord = results.coord
        grid = snapshot_to_grid(vals, coord=coord, ncell=ncell, reduce_z=reduce_z)
        return iter_id, grid

    # fallback: treat as restart-like table
    table = results.get(key)
    if table is None:
        return None

    if "iter" in table:
        iters = np.unique(table["iter"])
        iter_id = int(iters[-1])
        mask = table["iter"] == iter_id
    else:
        iter_id = None
        mask = slice(None)

    iatom = table["iatom"][mask].astype(int)
    if component == 0:
        vals = table["mx"][mask]
    elif component == 1:
        vals = table["my"][mask]
    else:
        vals = table["mz"][mask]

    # Build full per-atom array
    natom = int(table["iatom"].max())
    full = np.full(natom, np.nan, dtype=float)
    for ia, v in zip(iatom, vals):
        full[int(ia) - 1] = v

    coord = results.coord
    grid = snapshot_to_grid(full, coord=coord, ncell=ncell, reduce_z=reduce_z)
    return iter_id, grid


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
