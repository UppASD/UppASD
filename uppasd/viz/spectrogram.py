"""
spectrogram.py

Space-time maps for 1D systems.
"""

import matplotlib.pyplot as plt


def plot_xt_map(
    data,
    *,
    ax=None,
    cmap="inferno",
    origin="lower",
    colorbar=True,
    xlabel="x",
    ylabel="time",
    title=None,
    show=True,
    **imshow_kw,
):
    """
    Plot a space-time (x–t) map.

    Parameters
    ----------
    data : (Nt, Nx) array
    """
    if ax is None:
        fig, ax = plt.subplots()

    im = ax.imshow(
        data,
        cmap=cmap,
        origin=origin,
        aspect="auto",
        **imshow_kw,
    )

    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)

    if title:
        ax.set_title(title)

    if colorbar:
        plt.colorbar(im, ax=ax)

    if show:
        plt.show()

    return ax
