"""
Visualization utilities for UppASD.

This submodule provides lightweight, matplotlib-based helpers for
visualizing UppASD results, including:

- Time series extraction
- Site-resolved snapshots
- Heatmaps for 1D / 2D / reduced 3D systems
"""

# ------------------------------------------------------------
# Data extraction / reshaping (NO plotting)
# ------------------------------------------------------------

from .io import (
    timeseries,
    spin_snapshots,
    effective_field_snapshots,
    site_snapshots,
    reshape_supercell,
    snapshot_to_grid,
    results_snapshot_grid,
)

# ------------------------------------------------------------
# Plotting helpers (matplotlib only)
# ------------------------------------------------------------

from .heatmap import (
    plot_heatmap,
)
from .sq import (
    plot_sq,
)

# ------------------------------------------------------------
# Public API
# ------------------------------------------------------------

__all__ = [
    # io
    "timeseries",
    "spin_snapshots",
    "effective_field_snapshots",
    "site_snapshots",
    "reshape_supercell",
    "snapshot_to_grid",
    "results_snapshot_grid",

    # plotting
    "plot_heatmap",
    "plot_sq",
]
