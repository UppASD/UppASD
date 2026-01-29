"""
bands.py

Band-structure plotting helpers.
"""

import matplotlib.pyplot as plt


def plot_bands(
    k,
    energies,
    *,
    ax=None,
    xlabel="k",
    ylabel="Energy",
    title=None,
    show=True,
    **plot_kw,
):
    """
    Plot band structure.

    Parameters
    ----------
    k : (Nk,) array
    energies : (Nk, Nb) array
    """
    if ax is None:
        fig, ax = plt.subplots()

    for ib in range(energies.shape[1]):
        ax.plot(k, energies[:, ib], **plot_kw)

    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)

    if title:
        ax.set_title(title)

    if show:
        plt.show()

    return ax
