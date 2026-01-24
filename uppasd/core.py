"""
uppasd/core.py

Base classes for state management, error handling, and common utilities.

This module provides the foundation for all UppASD simulations, handling:
- Simulation lifecycle state machine
- Error handling and custom exceptions
- Logging configuration (compatible with Jupyter and CLI)
- Common properties and methods
"""

import logging
from enum import Enum
from typing import Optional

# Setup module logger (compatible with Jupyter and CLI)
logger = logging.getLogger("uppasd")
logger.setLevel(logging.INFO)

# Add handler only if none exist (avoid duplicate logs in notebooks)
if not logger.handlers:
    handler = logging.StreamHandler()
    formatter = logging.Formatter(
        "%(asctime)s [%(name)s] %(levelname)s: %(message)s",
        datefmt="%H:%M:%S"
    )
    handler.setFormatter(formatter)
    logger.addHandler(handler)


class SimulationState(Enum):
    """
    Simulation lifecycle states.
    
    States
    ------
    UNINITIALIZED
        Initial state after creating Simulator object.
        Allowed transitions: INITIALIZED
    INITIALIZED
        After init_simulation() called successfully.
        Allowed transitions: RUNNING, CLEANED_UP
    RUNNING
        During active simulation (run_simulation, relax, etc).
        Allowed transitions: INITIALIZED, COMPLETED, CLEANED_UP
    COMPLETED
        After simulation phase completed (not during execution).
        Allowed transitions: RUNNING, CLEANED_UP
    CLEANED_UP
        After cleanup() called. Final state.
        Allowed transitions: None (terminal state)
    """
    UNINITIALIZED = "uninitialized"
    INITIALIZED = "initialized"
    RUNNING = "running"
    COMPLETED = "completed"
    CLEANED_UP = "cleaned_up"


class UppASDError(Exception):
    """
    Base exception for UppASD-related errors.
    
    All UppASD-specific exceptions inherit from this class.
    """
    pass


class StateError(UppASDError):
    """
    Raised when an operation requires a different simulation state.
    
    Examples
    --------
    >>> sim = Simulator()
    >>> sim.moments  # StateError: must call init_simulation() first
    """
    pass


class InitializationError(UppASDError):
    """Raised when simulation initialization fails."""
    pass


class SimulationFailedError(UppASDError):
    """Raised when simulation execution fails."""
    pass


class UppASDBase:
    """
    Base class for UppASD simulator.
    
    Provides:
    - State machine for simulation lifecycle
    - Error wrapping around Fortran calls
    - Logging configuration (DEBUG, INFO, WARNING, ERROR, CRITICAL)
    - Common properties (natom, mensemble, temperature, etc.)
    - State validation via _check_state()
    
    This class should not be instantiated directly; use Simulator instead.
    
    Attributes
    ----------
    natom : int or None
        Number of atoms in system. Set during init_simulation().
    mensemble : int or None
        Number of independent spin ensembles/replicas. Set during init_simulation().
    state : str
        Current simulation state (read-only property).
    is_initialized : bool
        True if simulation has been initialized (read-only property).
    
    Parameters
    ----------
    verbose : bool, optional
        Enable verbose logging (DEBUG level). Default: False
    quiet : bool, optional
        Suppress all logging output. Default: False
    
    Examples
    --------
    Subclasses should implement specific functionality:
    
    >>> class MySimulator(UppASDBase):
    ...     def init_simulation(self):
    ...         self._check_state(SimulationState.UNINITIALIZED)
    ...         # ... initialization code ...
    ...         self._state = SimulationState.INITIALIZED
    """
    
    def __init__(self, verbose: bool = False, quiet: bool = False):
        """
        Initialize UppASD base simulator.
        
        Parameters
        ----------
        verbose : bool, optional
            Enable verbose logging (DEBUG level). Default: False
        quiet : bool, optional
            Suppress all logging. Default: False
            
        Raises
        ------
        ImportError
            If required dependencies (_uppasd C extension) not available.
        """
        self._state = SimulationState.UNINITIALIZED
        self.natom: Optional[int] = None
        self.mensemble: Optional[int] = None
        
        # Configure logging
        if quiet:
            logger.setLevel(logging.CRITICAL + 1)  # Disable (higher than CRITICAL)
        elif verbose:
            logger.setLevel(logging.DEBUG)
        else:
            logger.setLevel(logging.INFO)
        
        logger.debug(f"Initialized UppASDBase (verbose={verbose}, quiet={quiet})")
    
    def _check_state(self, *required_states: SimulationState) -> None:
        """
        Validate that current state is one of the required states.
        
        This method is used internally to enforce the state machine.
        It should be called at the start of any method that requires specific state.
        
        Parameters
        ----------
        required_states : SimulationState
            One or more acceptable states.
        
        Raises
        ------
        StateError
            If current state not in required_states, with helpful message.
        
        Examples
        --------
        >>> def my_method(self):
        ...     self._check_state(SimulationState.INITIALIZED)
        ...     # If we get here, we know _state == INITIALIZED
        """
        if self._state not in required_states:
            state_names = [s.value for s in required_states]
            msg = (
                f"Invalid state: '{self._state.value}'. "
                f"Expected one of: {state_names}"
            )
            logger.error(msg)
            raise StateError(msg)
    
    @property
    def state(self) -> str:
        """
        Get current simulation state as a string.
        
        Returns
        -------
        state : str
            One of: "uninitialized", "initialized", "running", 
            "completed", "cleaned_up"
        
        Examples
        --------
        >>> sim = Simulator()
        >>> print(sim.state)
        'uninitialized'
        >>> sim.init_simulation()
        >>> print(sim.state)
        'initialized'
        """
        return self._state.value
    
    @property
    def is_initialized(self) -> bool:
        """
        Check if simulation has been initialized.
        
        Returns
        -------
        bool
            True if state is one of: INITIALIZED, RUNNING, COMPLETED
            False if state is UNINITIALIZED or CLEANED_UP
        
        Examples
        --------
        >>> sim = Simulator()
        >>> sim.is_initialized
        False
        >>> sim.init_simulation()
        >>> sim.is_initialized
        True
        """
        return self._state in (
            SimulationState.INITIALIZED,
            SimulationState.RUNNING,
            SimulationState.COMPLETED
        )
    
    def __repr__(self) -> str:
        """String representation of simulator state."""
        status = f"state={self.state}"
        if self.is_initialized:
            status += f", natom={self.natom}, mensemble={self.mensemble}"
        return f"<{self.__class__.__name__}({status})>"
