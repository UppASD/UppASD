"""
UppASD Python interface.
Public, stable API.
"""

# --- Core data structures ---
from .core.system import SpinSystem
from .core.exchange import ExchangeShellTable, DMIShellTable
from .core.results import ASDResults

# --- Core helpers ---
from .core.neighbors import build_neighbors, assign_shells
from .core.units import (
    steps_to_time,
    steps_to_ps,
    steps_to_ns,
)

# --- Input handling ---
from .input.inputdata import ASDInput, InputBlock

# --- Execution ---
from .run.simulator import (
    ASDWorkspace,
    UppASDSimulator,
    run_relaxation,
    run_measurement,
)

__all__ = [
    "SpinSystem",
    "ExchangeShellTable",
    "DMIShellTable",
    "ASDResults",
    "ASDInput",
    "InputBlock",
    "ASDWorkspace",
    "UppASDSimulator",
]
