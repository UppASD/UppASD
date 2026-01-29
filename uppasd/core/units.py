"""
units.py

Simple time / unit conversion helpers for UppASD.

No global state, no physics assumptions.
"""

import numpy as np


# ======================================================================
# Time helpers
# ======================================================================


def steps_to_time(steps, timestep):
    """
    Convert simulation steps to time.

    Parameters
    ----------
    steps : array-like or scalar
        Simulation steps
    timestep : float
        Timestep in seconds

    Returns
    -------
    time : ndarray or float
        Time in seconds
    """
    return np.asarray(steps) * timestep


def steps_to_ps(steps, timestep):
    """
    Convert simulation steps to time in picoseconds.
    """
    return steps_to_time(steps, timestep) * 1e12


def steps_to_ns(steps, timestep):
    """
    Convert simulation steps to time in nanoseconds.
    """
    return steps_to_time(steps, timestep) * 1e9


# ======================================================================
# Result helpers
# ======================================================================


def magnetization_vs_time(results, timestep, component="mz", unit="ps"):
    """
    Convenience helper: magnetization vs time.

    Parameters
    ----------
    results : ASDResults
    timestep : float
        Timestep in seconds
    component : str
        'mx', 'my', 'mz', or 'm'
    unit : str
        's', 'ps', or 'ns'

    Returns
    -------
    t, m : ndarray, ndarray
    """
    steps = results.step
    if steps is None:
        return None, None

    if unit == "s":
        t = steps_to_time(steps, timestep)
    elif unit == "ps":
        t = steps_to_ps(steps, timestep)
    elif unit == "ns":
        t = steps_to_ns(steps, timestep)
    else:
        raise ValueError(f"Unknown time unit '{unit}'")

    m = getattr(results, component)
    return t, m
