"""
style.py

Matplotlib style presets for UppASD visualizations.

These are opt-in helpers. Nothing is applied unless explicitly called.
"""

import matplotlib as mpl
import matplotlib.pyplot as plt


# ---------------------------------------------------------------------
# Core style setter
# ---------------------------------------------------------------------
def set_style(params: dict):
    """
    Apply a matplotlib rcParams dictionary.

    Parameters
    ----------
    params : dict
        Matplotlib rcParams-style dictionary
    """
    mpl.rcParams.update(params)


# ---------------------------------------------------------------------
# Presets
# ---------------------------------------------------------------------
def use_light():
    """
    Light, clean style for interactive work.
    """
    set_style({
        "figure.figsize": (6.4, 4.8),
        "axes.grid": False,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.labelsize": 12,
        "axes.titlesize": 13,
        "xtick.labelsize": 11,
        "ytick.labelsize": 11,
        "lines.linewidth": 2.0,
        "image.interpolation": "nearest",
        "image.cmap": "viridis",
    })


def use_dark():
    """
    Dark background style suitable for spin textures and heatmaps.
    """
    set_style({
        "figure.facecolor": "black",
        "axes.facecolor": "black",
        "axes.edgecolor": "white",
        "axes.labelcolor": "white",
        "axes.titlecolor": "white",
        "xtick.color": "white",
        "ytick.color": "white",
        "text.color": "white",
        "lines.linewidth": 2.0,
        "image.cmap": "inferno",
    })


def use_publication():
    """
    Publication-quality style (journal-safe).

    Conservative font sizes, no grids, monochrome-friendly.
    """
    set_style({
        "figure.figsize": (3.4, 2.6),
        "font.size": 9,
        "axes.labelsize": 9,
        "axes.titlesize": 9,
        "xtick.labelsize": 8,
        "ytick.labelsize": 8,
        "lines.linewidth": 1.5,
        "axes.linewidth": 0.8,
        "axes.grid": False,
        "image.cmap": "gray",
        "savefig.dpi": 300,
        "savefig.bbox": "tight",
    })
