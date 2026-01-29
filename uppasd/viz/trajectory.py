"""
trajectory.py

Generic trajectory plotting utilities for UppASD.

Design philosophy
-----------------
- Geometry-based, not physics-based
- Reusable for skyrmion COM, domain walls, spin trajectories, etc.
- Minimal assumptions about origin of data
- Matplotlib-only dependency

Supported trajectories
----------------------
- 2D trajectories: (x(t), y(t))
- 3D trajectories: (x(t), y(t), z(t))
"""

from __future__ import annotations

import numpy as np
import matplotlib.pyplot as plt
from typing import Optional, Tuple


# ======================================================================
# 2D trajectories
# ======================================================================

def plot_2d_trajectory(
    x: np.ndarray,
    y: np.ndarray,
    *,
    ax: Optional[plt.Axes] = None,
    color: str = "C0",
    lw: float = 2.0,
    marker: Optional[str] = None,
    markersize: float = 4.0,
    start_marker: bool = True,
    end_marker: bool = True,
    equal_aspect: bool = True,
    xlabel: str = "x",
    ylabel: str = "y",
    title: Optional[str] = None,
):
    """
    Plot a 2D trajectory y(x).

    Parameters
    ----------
    x, y : ndarray (Nt,)
        Trajectory coordinates
    ax : matplotlib Axes, optional
        Target axes
    color : str
        Line color
    lw : float
        Line width
    marker : str, optional
        Marker style along trajectory
    markersize : float
        Marker size
    start_marker, end_marker : bool
        Mark start/end points explicitly
    equal_aspect : bool
        Enforce equal x/y scaling
    xlabel, ylabel : str
        Axis labels
    title : str, optional
        Plot title
    """
    if ax is None:
        fig, ax = plt.subplots()

    ax.plot(x, y, color=color, lw=lw, marker=marker, markersize=markersize)

    if start_marker:
        ax.scatter(x[0], y[0], color="green", s=60, zorder=3, label="start")

    if end_marker:
        ax.scatter(x[-1], y[-1], color="red", s=60, zorder=3, label="end")

    if equal_aspect:
        ax.set_aspect("equal", adjustable="box")

    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)

    if title:
        ax.set_title(title)

    if start_marker or end_marker:
        ax.legend(frameon=False)

    return ax


# ======================================================================
# 3D trajectories
# ======================================================================

def plot_3d_trajectory(
    x: np.ndarray,
    y: np.ndarray,
    z: np.ndarray,
    *,
    ax: Optional[plt.Axes] = None,
    color: str = "C0",
    lw: float = 2.0,
    show_sphere: bool = False,
    sphere_alpha: float = 0.15,
    xlabel: str = "x",
    ylabel: str = "y",
    zlabel: str = "z",
    title: Optional[str] = None,
):
    """
    Plot a 3D trajectory.

    Typical use:
    - Single-spin trajectory on unit sphere
    - 3D collective coordinate paths

    Parameters
    ----------
    x, y, z : ndarray (Nt,)
        Trajectory coordinates
    ax : mpl_toolkits.mplot3d Axes, optional
        Target axes
    color : str
        Line color
    lw : float
        Line width
    show_sphere : bool
        Draw unit sphere (useful for spin trajectories)
    sphere_alpha : float
        Transparency of sphere
    xlabel, ylabel, zlabel : str
        Axis labels
    title : str, optional
        Plot title
    """
    from mpl_toolkits.mplot3d import Axes3D  # noqa: F401

    if ax is None:
        fig = plt.figure()
        ax = fig.add_subplot(111, projection="3d")

    ax.plot(x, y, z, color=color, lw=lw)

    # Mark start/end
    ax.scatter(x[0], y[0], z[0], color="green", s=50, label="start")
    ax.scatter(x[-1], y[-1], z[-1], color="red", s=50, label="end")

    if show_sphere:
        _plot_unit_sphere(ax, alpha=sphere_alpha)

    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_zlabel(zlabel)

    if title:
        ax.set_title(title)

    ax.legend(frameon=False)

    return ax


# ======================================================================
# Helpers
# ======================================================================

def _plot_unit_sphere(ax, *, alpha: float = 0.2):
    """Draw a unit sphere for reference."""
    u = np.linspace(0, 2 * np.pi, 60)
    v = np.linspace(0, np.pi, 30)

    xs = np.outer(np.cos(u), np.sin(v))
    ys = np.outer(np.sin(u), np.sin(v))
    zs = np.outer(np.ones_like(u), np.cos(v))

    ax.plot_surface(
        xs,
        ys,
        zs,
        color="gray",
        alpha=alpha,
        linewidth=0,
        zorder=0,
    )

    ax.set_box_aspect((1, 1, 1))
