"""
sq.py

Static structure factor (S(q)) plotting helpers.

Provides utilities to convert flat (qx,qy,|S|) tables into a 2D grid
and plot it using the existing heatmap routines.
"""

import numpy as np
import matplotlib.pyplot as plt
from typing import Tuple, Optional

from .heatmap import plot_heatmap


def sq_to_grid(qx: np.ndarray, qy: np.ndarray, vals: np.ndarray, round_decimals: int = 8) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Convert flat qx,qy,values arrays into a regularly gridded 2D array.

    Parameters
    ----------
    qx, qy : ndarray
        Flat arrays of qx and qy coordinates.
    vals : ndarray
        Values at each (qx,qy) point (e.g. |S(q)|).
    round_decimals : int
        Number of decimals to round coordinates to for unique detection.

    Returns
    -------
    qx_u : ndarray (Nx,)
    qy_u : ndarray (Ny,)
    grid : ndarray (Ny, Nx)
        Grid arranged as (y, x) ready for `imshow` with origin='lower'.
    """
    qx_r = np.round(qx, round_decimals)
    qy_r = np.round(qy, round_decimals)

    qx_u = np.unique(qx_r)
    qy_u = np.unique(qy_r)

    nx = qx_u.size
    ny = qy_u.size

    grid = np.full((ny, nx), np.nan, dtype=float)

    # build index maps
    qx_index = {float(v): i for i, v in enumerate(qx_u)}
    qy_index = {float(v): i for i, v in enumerate(qy_u)}

    for x, y, v in zip(qx_r, qy_r, vals):
        ix = qx_index.get(float(x))
        iy = qy_index.get(float(y))
        if ix is None or iy is None:
            continue
        grid[iy, ix] = v

    return qx_u, qy_u, grid


def plot_sq(
    table,
    ax=None,
    cmap: str = "viridis",
    vmin: Optional[float] = None,
    vmax: Optional[float] = None,
    colorbar: bool = True,
    title: Optional[str] = None,
    show: bool = True,
    round_decimals: int = 8,
    fftshift: bool = False,
    **imshow_kwargs,
):
    """
    Plot static structure factor S(q) from a parsed table (as returned by
    `uppasd.core.output.read_sq`) or from raw arrays.

    Parameters
    ----------
    table : dict or tuple
        If dict, expected keys are 'qx','qy' and one of 'abs' or 'AbsS'.
        Alternatively a tuple (qx, qy, vals) may be provided.
    ax : matplotlib Axes, optional
    cmap, vmin, vmax, colorbar, title : plotting kwargs
    round_decimals : int
        Rounding applied to q coordinates when building the grid.
    **imshow_kwargs : passed to `imshow`
    """
    # Resolve inputs
    if isinstance(table, dict):
        qx = table.get("qx")
        if qx is None:
            qx = table.get("q_x")
        qx = np.asarray(qx)

        qy = table.get("qy")
        if qy is None:
            qy = table.get("q_y")
        qy = np.asarray(qy)
        vals = None
        if "abs" in table:
            vals = np.asarray(table["abs"])
        elif "AbsS" in table:
            vals = np.asarray(table["AbsS"])
        elif "Abs" in table:
            vals = np.asarray(table["Abs"])
        else:
            # fallback: try common names
            for k in ("res","resxx","resyy","reszz"):
                if k in table:
                    vals = np.asarray(table[k])
                    break
        if vals is None:
            raise KeyError("Could not find value column in sq table (expected 'abs' or 'AbsS')")
    else:
        qx, qy, vals = table

    qx_u, qy_u, grid = sq_to_grid(qx, qy, vals, round_decimals=round_decimals)

    # Optionally center the zero-frequency component (fftshift-like)
    if fftshift:
        # shift both coordinate axes and the grid for display
        qx_u = np.fft.fftshift(qx_u)
        qy_u = np.fft.fftshift(qy_u)
        grid = np.fft.fftshift(grid)

    # Use heatmap helper to plot
    ax = plot_heatmap(
        grid,
        ax=ax,
        title=title,
        xlabel="q_x",
        ylabel="q_y",
        cmap=cmap,
        vmin=vmin,
        vmax=vmax,
        colorbar=colorbar,
        colorbar_label="|S(q)|",
        aspect="auto",
        **imshow_kwargs,
    )

    # set ticks to q values (sparse labeling)
    nx = qx_u.size
    ny = qy_u.size
    if nx <= 10:
        ax.set_xticks(np.arange(nx))
        ax.set_xticklabels([f"{v:.3f}" for v in qx_u])
    if ny <= 10:
        ax.set_yticks(np.arange(ny))
        ax.set_yticklabels([f"{v:.3f}" for v in qy_u])

    if show:
        plt.show()

    return ax
