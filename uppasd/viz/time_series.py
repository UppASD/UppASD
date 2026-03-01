"""Time-series plotting utilities for UppASD.

This module is a relocation of the legacy `uppasd.plotting` helpers into
the `uppasd.viz` package for consistent discovery alongside other
visualization utilities (`heatmap`, `trajectory`).

Public API:
- `to_physical_time(time_array, timestep=None, unit='fs')`
- `plot_series(series_map, ...)`
"""
from typing import Dict, Any, Optional, Tuple
import numpy as np
import matplotlib.pyplot as plt


_TIME_UNITS = {
    's': 1.0,
    'ms': 1e-3,
    'us': 1e-6,
    'ns': 1e-9,
    'ps': 1e-12,
    'fs': 1e-15,
}


def to_physical_time(time_array, timestep: Optional[float] = None, unit: str = 'fs') -> np.ndarray:
    """Convert a numeric time/step index array to physical units.

    Parameters
    ----------
    time_array : array-like
        Numeric array of time indices (e.g. iteration counts or native time values).
    timestep : float or None
        Timestep in seconds. If provided, the returned array is `time_array * timestep` converted to `unit`.
        If ``None``, `time_array` is returned as a float numpy array unchanged.
    unit : str
        Output unit, one of 's', 'ms', 'us', 'ns', 'ps', 'fs'. Default 'fs'.

    Returns
    -------
    numpy.ndarray
        Time values in requested unit.
    """
    t = np.asarray(time_array, dtype=float)
    if timestep is None:
        return t
    denom = _TIME_UNITS.get(unit, _TIME_UNITS['fs'])
    # time in requested units = (time_array * timestep_seconds) / denom
    return t * (timestep / denom)


def plot_series(series_map: Dict[str, Dict[str, Any]], *,
                to_physical: bool = True,
                timestep: Optional[float] = None,
                time_unit: str = 'fs',
                connect: bool = False,
                show_points: bool = True,
                ax: Optional[plt.Axes] = None,
                figsize: Tuple[int, int] = (10, 5),
                title: Optional[str] = None) -> Tuple[plt.Figure, plt.Axes]:
    """Plot multiple series where each series carries its own time array.

    Parameters
    ----------
    series_map : dict
        Mapping name -> dict with keys: 'time' (array-like), 'values' (array-like).
        Optional keys: 'label', 'color', 'timestep' (per-series override).
    to_physical : bool
        If True, convert time indices to physical units using per-series 'timestep' or the
        provided `timestep` argument.
    timestep : float or None
        Default timestep (seconds) used when a series doesn't provide its own 'timestep'.
    time_unit : str
        Unit for x-axis when converting physical time: 'fs' by default.
    connect : bool
        If True, draw connecting lines between sample points. Otherwise draw only markers
        for each sample and optionally a faint connecting line.
    show_points : bool
        Whether to show sample markers.
    ax : matplotlib.axes.Axes or None
        Axis to plot on; if None a new figure/axis is created.
    figsize : tuple
        Figure size when creating a new figure.
    title : str or None
        Optional plot title.

    Returns
    -------
    (fig, ax)
        The matplotlib figure and axis used.
    """
    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)
    else:
        fig = ax.figure

    for key, s in series_map.items():
        t = np.asarray(s.get('time', []), dtype=float)
        v = np.asarray(s.get('values', []), dtype=float)
        label = s.get('label', key)
        color = s.get('color', None)
        series_timestep = s.get('timestep', None)

        if to_physical:
            used_ts = series_timestep if series_timestep is not None else timestep
            if used_ts is not None:
                t = to_physical_time(t, used_ts, unit=time_unit)

        if t.size == 0 or v.size == 0:
            # Nothing to plot for this series
            continue

        # If values are scalar (single point), broadcast if x has >1 points
        if v.size == 1 and t.size > 1:
            v = np.full(t.shape, v.item(), dtype=float)

        # Align lengths conservatively (trim longer axis)
        if v.size != t.size:
            min_len = min(v.size, t.size)
            if min_len == 0:
                continue
            v = v[:min_len]
            t = t[:min_len]

        if connect:
            ax.plot(t, v, '-', color=color, linewidth=1.5, label=label)
        else:
            # faint line for readability but emphasize samples
            ax.plot(t, v, '-', color=color, linewidth=0.8, alpha=0.4)
        if show_points:
            ax.plot(t, v, 'o', color=color, markersize=4, label=(None if connect else label))

    ax.set_xlabel(f"Time ({time_unit})" if to_physical else "Iteration")
    ax.set_ylabel('Value')
    if title:
        ax.set_title(title)
    ax.grid(True, alpha=0.3)
    ax.legend(loc='best')
    plt.tight_layout()
    return fig, ax
