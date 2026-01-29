"""
reductions.py

Explicit reduction utilities for UppASD spatial data.

Expected base data shape:
    data[t, ix, iy, iz, iatom, ...]
"""

import numpy as np


def reduce_ucell_mean(data):
    """
    Reduce over atoms in the unit cell AND z-direction.

    Resulting shape:
        data[t, ix, iy]

    Suitable for:
    - skyrmion lattices
    - continuum-like fields
    """
    return data.mean(axis=(3, 4))


def reduce_per_atom(data):
    """
    Reduce over z-direction only.

    Resulting shape:
        data[t, ix, iy, iatom]

    Suitable for:
    - sublattice-resolved visualization
    """
    return data.mean(axis=3)


def reduce_z_slice(data, iz):
    """
    Extract a single z-slice.

    Parameters
    ----------
    iz : int
        Z-index to extract

    Resulting shape:
        data[t, ix, iy, iatom]
    """
    return data[:, :, :, iz, :]


def reduce_time_average(data):
    """
    Average over time.

    Resulting shape depends on input.
    """
    return data.mean(axis=0)
