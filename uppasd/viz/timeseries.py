"""
timeseries.py

Lightweight helpers for plotting time-series data.

Design principles:
- Zero assumptions about units or normalization
- No hidden behavior
- Full matplotlib control via **kwargs
"""

import matplotlib.pyplot as plt


def plot_timeseries(
    x,
    y,
    *,
    ax=None,
    xlabel=None,
    ylabel=None,
    title=None,
    show=True,
    **plot_kw,
):
    """
    Plot a 1D time series.

    Parameters
    ----------
    x, y : array-like
        X- and Y-data (e.g. step and observable)
    ax : matplotlib.axes.Axes, optional
        Existing axes to plot into
    xlabel, ylabel, title : str, optional
        Axis labels and title
    show : bool, default True
        Call plt.show()
    **plot_kw :
        Passed directly to ax.plot()

    Returns
    -------
    ax : matplotlib.axes.Axes
    """
    if ax is None:
        fig, ax = plt.subplots()

    ax.plot(x, y, **plot_kw)

    if xlabel:
        ax.set_xlabel(xlabel)
    if ylabel:
        ax.set_ylabel(ylabel)
    if title:
        ax.set_title(title)

    if show:
        plt.show()

    return ax
