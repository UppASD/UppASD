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
from .input.inputdata import (
    ASDInput,
    InputBlock,
    set_temperature_token,
    validate_temperature_token,
)

# --- Execution ---
from .run.simulator import (
    ASDWorkspace,
    UppASDSimulator,
    run_relaxation,
    run_measurement,
)
from .run.sweep import run_temperature_sweep

__all__ = [
    "SpinSystem",
    "ExchangeShellTable",
    "DMIShellTable",
    "ASDResults",
    "ASDInput",
    "InputBlock",
    "set_temperature_token",
    "validate_temperature_token",
    "ASDWorkspace",
    "UppASDSimulator",
    "run_temperature_sweep",
]
