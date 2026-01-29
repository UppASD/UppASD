"""
output.py

Deterministic NumPy-based parsers for UppASD output files.

Design principles:
- No globbing
- No guessing
- Deterministic filenames: <stem>.<simid>.out
- simid defaults to Fortran default "_UppASD_"
- Parsing only (no reshaping, no physics logic)
"""

from pathlib import Path
import numpy as np
from typing import Optional, Dict


# ======================================================================
# simid handling
# ======================================================================

DEFAULT_SIMID = "_UppASD_"


def output_path(workdir, stem: str, simid: Optional[str] = None) -> Path:
    """
    Construct deterministic UppASD output filename.

    Parameters
    ----------
    workdir : str or Path
        UppASD run directory
    stem : str
        Output stem, e.g. "averages", "totenergy", "moment"
    simid : str, optional
        Simulation ID (<= 8 chars). Defaults to "_UppASD_"

    Returns
    -------
    Path
    """
    if simid is None:
        simid = DEFAULT_SIMID

    return Path(workdir) / f"{stem}.{simid}.out"


# ======================================================================
# Averaged quantities
# ======================================================================


def read_averages(workdir, simid: Optional[str] = None) -> Dict:
    """
    Read averages.<simid>.out

    Columns:
        Iter, <Mx>, <My>, <Mz>, <M>, M_std
    """
    path = output_path(workdir, "averages", simid)
    data = np.loadtxt(path, comments="#")

    return {
        "iter": data[:, 0].astype(int),
        "mx": data[:, 1],
        "my": data[:, 2],
        "mz": data[:, 3],
        "m": data[:, 4],
        "m_std": data[:, 5],
        "raw": data,
    }


# ======================================================================
# Energy
# ======================================================================


def read_totenergy(workdir, simid: Optional[str] = None) -> Dict:
    """
    Read totenergy.<simid>.out

    Columns:
        Iter, Tot, Exc, Ani, DM, PD, BiqDM, BQ, Dip, Zeeman,
        LSF, Chir, Ring, SA
    """
    path = output_path(workdir, "totenergy", simid)
    data = np.loadtxt(path, comments="#")

    return {
        "iter": data[:, 0].astype(int),
        "tot": data[:, 1],
        "exc": data[:, 2],
        "ani": data[:, 3],
        "dm": data[:, 4],
        "pd": data[:, 5],
        "biqdm": data[:, 6],
        "bq": data[:, 7],
        "dip": data[:, 8],
        "zeeman": data[:, 9],
        "lsf": data[:, 10],
        "chir": data[:, 11],
        "ring": data[:, 12],
        "sa": data[:, 13],
        "raw": data,
    }


# ======================================================================
# Cumulants
# ======================================================================


def read_cumulants(workdir, simid: Optional[str] = None) -> Dict:
    """
    Read cumulants.<simid>.out
    """
    path = output_path(workdir, "cumulants", simid)
    data = np.loadtxt(path, comments="#")

    return {
        "iter": data[:, 0].astype(int),
        "m": data[:, 1],
        "m2": data[:, 2],
        "m4": data[:, 3],
        "binder": data[:, 4],
        "chi": data[:, 5],
        "cv": data[:, 6],
        "energy": data[:, 7],
        "energy_exc": data[:, 8],
        "energy_lsf": data[:, 9],
        "raw": data,
    }


# ======================================================================
# Skyrmion number
# ======================================================================


def read_sknumber(workdir, simid: Optional[str] = None) -> Dict:
    """
    Read sknumber.<simid>.out
    """
    path = output_path(workdir, "sknumber", simid)
    data = np.loadtxt(path, comments="#")

    return {
        "iter": data[:, 0].astype(int),
        "sk": data[:, 1],
        "sk_avg": data[:, 2],
        "sk_std": data[:, 3],
        "raw": data,
    }


# ======================================================================
# Atom-resolved outputs
# ======================================================================


def read_restart(workdir, simid: Optional[str] = None) -> Dict:
    """
    Read restart.<simid>.out
    (final configuration)
    """
    path = output_path(workdir, "restart", simid)
    data = np.loadtxt(path, comments="#")

    return {
        "iter": data[:, 0].astype(int),
        "replica": data[:, 1].astype(int),
        "iatom": data[:, 2].astype(int),
        "magnitude": data[:, 3],
        "mx": data[:, 4],
        "my": data[:, 5],
        "mz": data[:, 6],
        "raw": data,
    }


def read_moment(workdir, simid: Optional[str] = None) -> Dict:
    """
    Read moment.<simid>.out
    (snapshots during simulation)
    """
    path = output_path(workdir, "moment", simid)
    data = np.loadtxt(path, comments="#")

    return {
        "iter": data[:, 0].astype(int),
        "replica": data[:, 1].astype(int),
        "iatom": data[:, 2].astype(int),
        "magnitude": data[:, 3],
        "mx": data[:, 4],
        "my": data[:, 5],
        "mz": data[:, 6],
        "raw": data,
    }


def read_befftot(workdir, simid: Optional[str] = None) -> Dict:
    """
    Read befftot.<simid>.out
    """
    path = output_path(workdir, "befftot", simid)
    data = np.loadtxt(path, comments="#")

    return {
        "iter": data[:, 0].astype(int),
        "iatom": data[:, 1].astype(int),
        "replica": data[:, 2].astype(int),
        "bx": data[:, 3],
        "by": data[:, 4],
        "bz": data[:, 5],
        "b": data[:, 6],
        "raw": data,
    }


def read_coord(workdir, simid: Optional[str] = None) -> Dict:
    """
    Read coord.<simid>.out
    (expanded coordinates for full system)
    """
    path = output_path(workdir, "coord", simid)
    data = np.loadtxt(path)

    return {
        "iatom": data[:, 0].astype(int),
        "x": data[:, 1],
        "y": data[:, 2],
        "z": data[:, 3],
        "species": data[:, 4].astype(int),
        "replica": data[:, 5].astype(int),
        "raw": data,
    }

def read_site_scalar(workdir, simid, prefix):
    """
    Read site-resolved scalar observable.

    Expected format:
        iter  site  value

    Example:
        dens_skynum.<simid>.out
    """
    import numpy as np
    from pathlib import Path

    filename = f"{prefix}.{simid}.out"
    path = Path(workdir) / filename

    data = np.loadtxt(path, comments="#")

    if data.ndim == 1:
        data = data[None, :]

    if data.shape[1] < 3:
        raise ValueError(f"{filename} must have at least 3 columns")

    return {
        "iter": data[:, 0].astype(int),
        "site": data[:, 1].astype(int),
        "value": data[:, 2],
        "raw": data,
    }
    
def read_cmass_skynum(workdir, simid):
    path = Path(workdir) / f"cmass_skynum.{simid}.out"
    data = np.loadtxt(path, comments="#")

    return {
        "iter": data[:, 0].astype(int),
        "rx": data[:, 1],
        "ry": data[:, 2],
        "rz": data[:, 3],
        "raw": data,
    }
    
def read_trajectories(workdir, simid):
    """
    Read all trajectory.<simid>.<ispin>.<iens>.out files.

    Returns
    -------
    dict[(ispin, iens)] -> table dict
    """
    traj = {}

    pattern = f"trajectory.{simid}.*.*.out"
    for path in Path(workdir).glob(pattern):
        parts = path.name.split(".")
        ispin = int(parts[-3])
        iens  = int(parts[-2])

        data = np.loadtxt(path, comments="#")

        traj[(ispin, iens)] = {
            "iter": data[:, 0].astype(int),
            "mx": data[:, 2],
            "my": data[:, 3],
            "mz": data[:, 4],
            "m":  data[:, 5],
            "raw": data,
        }

    return traj



def read_all_outputs(workdir):
    """
    Read all supported UppASD output files from a directory.

    Returns
    -------
    tables : dict
        Mapping name -> parsed table
    simid : str
        Simulation identifier used in filenames
    """
    workdir = Path(workdir)
    tables = {}

    # ------------------------------------------------------------
    # Determine simid
    # ------------------------------------------------------------
    simid = None

    # Strategy:
    # 1. Look for known output patterns
    # 2. Extract simid from filename
    for fname in workdir.iterdir():
        if not fname.is_file():
            continue
        parts = fname.name.split(".")
        if len(parts) >= 3:
            # e.g. averages.SCsurf_T.out
            simid = parts[-2]
            break

    if simid is None:
        simid = "_UppASD_"

    # ------------------------------------------------------------
    # Averaged / global outputs
    # ------------------------------------------------------------
    try:
        tables["averages"] = read_averages(workdir, simid)
    except FileNotFoundError:
        pass

    try:
        tables["cumulants"] = read_cumulants(workdir, simid)
    except FileNotFoundError:
        pass

    try:
        tables["totenergy"] = read_totenergy(workdir, simid)
    except FileNotFoundError:
        pass

    # ------------------------------------------------------------
    # Local / snapshot outputs
    # ------------------------------------------------------------
    try:
        tables["moment"] = read_moment(workdir, simid)
    except FileNotFoundError:
        pass

    try:
        tables["befftot"] = read_befftot(workdir, simid)
    except FileNotFoundError:
        pass

    try:
        tables["dens_skynum"] = read_site_scalar(workdir, simid, "dens_skynum")
    except FileNotFoundError:
        pass

    try:
        tables["restart"] = read_restart(workdir, simid)
    except FileNotFoundError:
        pass

    try:
        tables["coord"] = read_coord(workdir, simid)
    except FileNotFoundError:
        pass

    # cmass_skynum
    try:
        tables["cmass_skynum"] = read_cmass_skynum(workdir, simid)
    except FileNotFoundError:
        pass
    
    # trajectories
    traj = read_trajectories(workdir, simid)
    if traj:
        tables["trajectory"] = traj
    
    return tables, simid
