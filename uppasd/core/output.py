"""
output.py

Minimal NumPy-based parsers for UppASD output files.

Assumptions:
- Numeric, column-based files
- Comment headers starting with '#'
- Uniform row structure
"""

from pathlib import Path
import numpy as np
from typing import Dict, Optional, Union


PathLike = Union[str, Path]


# ======================================================================
# Internal helper
# ======================================================================


def _load_table(path: Path) -> np.ndarray:
    """
    Load a numeric table, handling empty files gracefully.
    """
    try:
        data = np.loadtxt(path, comments="#")
    except ValueError:
        # Empty or malformed file → treat as no data
        return np.empty((0, 0))

    if data.ndim == 1:
        data = data[None, :]

    return data


# ======================================================================
# Magnetization
# ======================================================================


def read_magnetization(
    workdir: PathLike,
    filename: str = "magnetization.out",
):
    """
    Read UppASD magnetization output.

    Typical columns:
        step  Mx  My  Mz  |M|

    Returns:
        dict with keys:
            step, mx, my, mz, m
    """
    path = Path(workdir) / filename
    data = _load_table(path)

    if data.size == 0:
        return None

    if data.shape[1] < 5:
        raise ValueError(f"{filename} must have at least 5 columns")

    return {
        "step": data[:, 0],
        "mx": data[:, 1],
        "my": data[:, 2],
        "mz": data[:, 3],
        "m": data[:, 4],
    }


# ======================================================================
# Energy
# ======================================================================


def read_energy(
    workdir: PathLike,
    filename: str = "energy.out",
):
    """
    Read UppASD energy output.

    Typical columns:
        step  E_total  [other contributions...]

    Returns:
        dict with keys:
            step, energy, raw
    """
    path = Path(workdir) / filename
    data = _load_table(path)

    if data.size == 0:
        return None

    if data.shape[1] < 2:
        raise ValueError(f"{filename} must have at least 2 columns")

    return {
        "step": data[:, 0],
        "energy": data[:, 1],
        "raw": data,
    }


# ======================================================================
# Generic helper
# ======================================================================


def read_table(
    workdir: PathLike,
    filename: str,
):
    """
    Generic numeric table reader for UppASD outputs.

    Returns:
        ndarray (N, M) or empty array if file is empty.
    """
    path = Path(workdir) / filename
    return _load_table(path)


# ======================================================================
# Combined convenience
# ======================================================================


def read_outputs(
    workdir: PathLike,
    magnetization: bool = True,
    energy: bool = True,
):
    """
    Read available UppASD outputs from a workspace.

    Returns:
        dict with optional keys:
            'magnetization'
            'energy'
    """
    outputs: Dict[str, Optional[dict]] = {}

    if magnetization:
        try:
            outputs["magnetization"] = read_magnetization(workdir)
        except OSError:
            outputs["magnetization"] = None

    if energy:
        try:
            outputs["energy"] = read_energy(workdir)
        except OSError:
            outputs["energy"] = None

    return outputs
