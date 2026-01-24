"""
UppASD: Uppsala Atomistic Spin Dynamics

A Python interface to the UppASD Fortran library for atomistic spin dynamics
and Monte Carlo simulations.

UppASD enables computational studies of magnetization dynamics, phase transitions,
and magnetic properties of materials through atomistic spin dynamics simulations.

Quick Start
-----------

**Interactive use (Jupyter notebook):**

>>> from uppasd import Simulator
>>> with Simulator() as sim:
...     sim.run_simulation()
...     print(f"Energy: {sim.energy}")
...     print(f"Magnetization: {sim.moments}")

**Script with error handling:**

>>> from uppasd import Simulator, StateError
>>> sim = Simulator(verbose=True)
>>> try:
...     sim.init_simulation()
...     sim.relax(mode='M', temperature=100, steps=1000)
...     print(f"Final energy: {sim.energy}")
... except StateError as e:
...     print(f"Simulation error: {e}")
... finally:
...     sim.cleanup()

**Low-level access (advanced users):**

>>> import uppasd.pyasd as asd
>>> natom, mensemble = asd.setup_all()
>>> moments = asd.relax(natom, mensemble, mode='M', steps=100)
>>> asd.cleanup()

Features
--------

- **Context manager support**: Automatic resource cleanup in Jupyter notebooks
- **State validation**: Prevents incorrect operation sequences
- **Logging**: Built-in progress tracking compatible with notebooks and CLI
- **Temperature sweeps**: Easy parameter studies
- **Callbacks**: Real-time monitoring during long simulations
- **Error handling**: Clear error messages and exception types

Documentation
--------------

For detailed documentation and examples, see:
- PYTHON_INTERFACE_INVENTORY.md - Overview of Python interface
- SIMULATOR_REFACTORING_PLAN.md - Architecture details and examples
- examples/notebooks/ - Interactive Jupyter notebooks

Version
-------

Version information and changelog: see CHANGELOG.md

"""

__version__ = "5.0.0"
__author__ = "UppASD Group"
__email__ = "uppasd@physics.uu.se"
__url__ = "https://github.com/UppASD/UppASD"

# Core classes and exceptions (available without C extension)
from uppasd.core import (
    UppASDBase,
    UppASDError,
    StateError,
    InitializationError,
    SimulationFailedError,
    SimulationState,
    logger,
)

# Lazy imports: only import when needed (avoid circular dependencies with C extension)
def __getattr__(name):
    """Lazy import for heavy modules that require C extension."""
    if name == "Simulator":
        from uppasd.simulator import Simulator
        return Simulator
    elif name == "pyasd":
        import uppasd.pyasd as pyasd
        return pyasd
    elif name == "InputData":
        from uppasd.inputdata import InputData
        return InputData
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")

# Public API - what should be imported with `from uppasd import *`
__all__ = [
    "Simulator",                 # Main user-facing class
    "UppASDBase",               # Base class (rarely used directly)
    "UppASDError",              # Base exception
    "StateError",               # Common exception
    "InitializationError",      # Initialization exception
    "SimulationFailedError",    # Runtime exception
    "SimulationState",          # Enum for state inspection
    "logger",                   # Module logger for customization
    "pyasd",                    # Low-level API (advanced)
    "InputData",                # Parameter management
]


