"""
sqw.py

S(q, ω) heatmap visualization.
"""

import matplotlib.pyplot as plt
import numpy as np


def plot_sqw(
    sqw,
    q,
    w,
    *,
    ax=None,
    cmap="magma",
    logscale=False,
    colorbar=True,
    title=None,
    show=True,
    **imshow_kw,
):
    """
    Plot S(q, ω).

    Parameters
    ----------
    sqw : (Nq, Nw) array
    q : array-like
    w : array-like
    """
    data = np.log(sqw + 1e-12) if logscale else sqw

    if ax is None:
        fig, ax = plt.subplots()

    im = ax.imshow(
        data.T,
        extent=[q.min(), q.max(), w.min(), w.max()],
        origin="lower",
        aspect="auto",
        cmap=cmap,
        **imshow_kw,
    )

    ax.set_xlabel("q")
    ax.set_ylabel("ω")

    if title:
        ax.set_title(title)

    if colorbar:
        plt.colorbar(im, ax=ax)

    if show:
        plt.show()

    return ax
