
"""
viz/heatmap.py

Lightweight Matplotlib-based heatmap utilities for UppASD.

Design principles
-----------------
- Regular grids only (use viz.io for reshaping/reduction)
- No UppASD knowledge
- Matplotlib-only dependency
- kwargs forwarded for fine control
"""

import numpy as np
import matplotlib.pyplot as plt
from typing import Optional, Tuple, Sequence


def plot_heatmap(
    data,
    ax=None,
    title=None,
    xlabel=None,
    ylabel=None,
    cmap="viridis",
    vmin=None,
    vmax=None,
    colorbar=True,
    colorbar_label=None,
    aspect="equal",
    **imshow_kwargs,
):
    """
    Plot a 2D heatmap using matplotlib.

    Parameters
    ----------
    data : 2D ndarray
        Heatmap data.
    ax : matplotlib Axes, optional
        Existing axes to plot into.
    title, xlabel, ylabel : str, optional
        Axis labels.
    cmap : str
        Colormap.
    vmin, vmax : float, optional
        Color limits.
    colorbar : bool
        Whether to draw a colorbar.
    colorbar_label : str, optional
        Label for the colorbar.
    aspect : str
        Aspect ratio passed to imshow.
    **imshow_kwargs
        Additional arguments forwarded to imshow (NOT colorbar).
    """
    if ax is None:
        fig, ax = plt.subplots()

    im = ax.imshow(
        data,
        origin="lower",
        cmap=cmap,
        vmin=vmin,
        vmax=vmax,
        aspect=aspect,
        **imshow_kwargs,
    )

    if title:
        ax.set_title(title)
    if xlabel:
        ax.set_xlabel(xlabel)
    if ylabel:
        ax.set_ylabel(ylabel)

    if colorbar:
        cbar = plt.colorbar(im, ax=ax)
        if colorbar_label:
            cbar.set_label(colorbar_label)

    return ax

def plot_multi_heatmap(
    data: Optional[np.ndarray] = None,
    *,
    data_list: Optional[Sequence[np.ndarray]] = None,
    ncols: int = 2,
    titles: Optional[Sequence[str]] = None,
    cmap: str = "RdBu_r",
    **kwargs,
):
    """
    Plot multiple heatmaps arranged in a grid.

    Parameters
    ----------
    data : ndarray, optional
        Array of shape (Nplots, N1, N2). Kept for backward compatibility.
    data_list : sequence of ndarray, optional
        Preferred interface. List/tuple of 2D arrays (N1, N2).
    ncols : int
        Number of subplot columns.
    titles : sequence of str, optional
        Titles for each subplot.
    cmap : str
        Colormap.
    **kwargs
        Passed through to plot_heatmap (except 'cmap').

    Returns
    -------
    fig : matplotlib.figure.Figure
    axes : ndarray of matplotlib.axes.Axes
    """

    # ------------------------------------------------------------
    # Resolve input data
    # ------------------------------------------------------------
    if data_list is not None:
        data_arr = np.asarray(data_list)
    elif data is not None:
        data_arr = np.asarray(data)
    else:
        raise TypeError(
            "plot_multi_heatmap requires either 'data' or 'data_list'"
        )

    # ------------------------------------------------------------
    # Normalize shape
    # ------------------------------------------------------------
    if data_arr.ndim == 2:
        data_arr = data_arr[None, ...]  # single plot
    elif data_arr.ndim != 3:
        raise ValueError(
            "Expected data with shape (Nplots, N1, N2) "
            f"or list of 2D arrays, got shape {data_arr.shape}"
        )

    nplots = data_arr.shape[0]
    nrows = int(np.ceil(nplots / ncols))

    # ------------------------------------------------------------
    # Create figure
    # ------------------------------------------------------------
    fig, axes = plt.subplots(
        nrows,
        ncols,
        squeeze=False,
        figsize=(4 * ncols, 4 * nrows),
    )
    axes = axes.ravel()

    # ------------------------------------------------------------
    # Plot each heatmap
    # ------------------------------------------------------------
    for i in range(nplots):
        title = titles[i] if titles and i < len(titles) else None
        plot_heatmap(
            data_arr[i],
            ax=axes[i],
            title=title,
            cmap=cmap,
            **kwargs,
        )

    # ------------------------------------------------------------
    # Hide unused axes
    # ------------------------------------------------------------
    for ax in axes[nplots:]:
        ax.axis("off")

    return fig, axes
