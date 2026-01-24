"""
High-level Simulator class for UppASD atomistic spin dynamics simulations.

The Simulator class provides a clean, Pythonic interface for interactive use
in Jupyter notebooks, scripts, and CLI applications.

Features:
- Context manager support (automatic cleanup)
- State validation (prevents calling methods in wrong order)
- Properties for common simulation parameters
- Comprehensive logging and error handling
- Optional trajectory recording for analysis
- Callback support for progress tracking

Examples
--------
**Interactive/Jupyter usage:**

>>> from uppasd import Simulator
>>> with Simulator() as sim:
...     sim.init_simulation()
...     print(f"System: {sim.natom} atoms")
...     sim.run_simulation()
...     energy = sim.energy
...     moments = sim.moments

**Script usage with error handling:**

>>> from uppasd import Simulator, StateError
>>> sim = Simulator(verbose=True)
>>> try:
...     sim.init_simulation()
...     sim.relax(mode='M', temperature=100, steps=1000)
...     results = sim.get_trajectory()
... except StateError as e:
...     print(f"Simulation error: {e}")
... finally:
...     sim.cleanup()

**Parameter sweep (Jupyter):**

>>> from uppasd import Simulator
>>> temperatures = np.linspace(10, 500, 20)
>>> magnetizations = []
>>> 
>>> for T in temperatures:
...     with Simulator() as sim:
...         sim.init_simulation()
...         sim.temperature = T
...         sim.relax(mode='M', temperature=T, steps=500)
...         m = np.mean(np.linalg.norm(sim.moments, axis=0))
...         magnetizations.append(m)
"""

import numpy as np
import logging
from typing import Optional, Callable, Dict, List, Any, Tuple

from uppasd.core import UppASDBase, SimulationState, StateError, logger
import uppasd.pyasd as pyasd

# Configure logger for this module
logger = logging.getLogger("uppasd.simulator")
class Simulator(UppASDBase):
    """
    High-level interface to UppASD atomistic spin dynamics simulations.
    
    Provides clean, Pythonic API for initializing, running, and analyzing
    spin simulations in Jupyter notebooks, scripts, and CLI applications.
    
    Features:
    - Context manager support (automatic cleanup with ``with`` statement)
    - State validation (prevents calling methods in wrong order)
    - Properties for common simulation parameters
    - Logging and optional progress tracking
    - Trajectory recording for analysis
    
    Parameters
    ----------
    verbose : bool, optional
        Enable verbose (DEBUG) logging. Default: False
    quiet : bool, optional
        Suppress logging. Default: False
    record_trajectory : bool, optional
        Automatically record moments after each step.
        Default: False (use get_trajectory() for manual recording)
    
    Raises
    ------
    ImportError
        If _uppasd C extension not available
    
    Examples
    --------
    **Context manager (recommended):**
    
    >>> from uppasd import Simulator
    >>> with Simulator() as sim:
    ...     sim.relax(mode='M', temperature=100)
    ...     print(f"Energy: {sim.energy}")
    
    **Manual setup/cleanup:**
    
    >>> sim = Simulator()
    >>> try:
    ...     sim.init_simulation()
    ...     moments = sim.relax(steps=500)
    ... finally:
    ...     sim.cleanup()
    
    **Parameter sweep with callbacks:**
    
    >>> def progress(step, moments, energy):
    ...     if step % 100 == 0:
    ...         print(f"Step {step}: E={energy:.3f}")
    >>> 
    >>> with Simulator() as sim:
    ...     sim.relax(mode='M', steps=1000, callback=progress)
    """
    
    def __init__(
        self,
        verbose: bool = False,
        quiet: bool = False,
        record_trajectory: bool = False,
    ):
        """Initialize Simulator."""
        super().__init__(verbose=verbose, quiet=quiet)
        
        self._record_trajectory = record_trajectory
        self._trajectory: List[np.ndarray] = []
        self._trajectory_metadata: Dict[str, Any] = {}
        
        # Will be set by init_simulation()
        self.natom: int = 0
        self.mensemble: int = 0
        
        logger.debug(
            f"Simulator initialized: record_trajectory={record_trajectory}"
        )
    
    # =========================================================================
    # Context Manager Protocol
    # =========================================================================
    
    def __enter__(self):
        """
        Context manager entry: initialize simulation.
        
        Returns
        -------
        self : Simulator
        
        Examples
        --------
        >>> with Simulator() as sim:
        ...     sim.relax()
        ...     # Automatically calls cleanup() on exit
        """
        self.init_simulation()
        return self
    
    def __exit__(self, exc_type, exc_val, exc_tb):
        """
        Context manager exit: cleanup resources.
        
        Parameters
        ----------
        exc_type : type or None
            Exception type if error occurred
        exc_val : Exception or None
            Exception value if error occurred
        exc_tb : traceback or None
            Traceback if error occurred
        """
        self.cleanup()
        if exc_type is not None:
            logger.error(f"Simulation exited with {exc_type.__name__}: {exc_val}")
    
    def __repr__(self) -> str:
        """Return string representation."""
        state = self._state.name if hasattr(self, '_state') else "UNKNOWN"
        if self.natom > 0:
            return (
                f"Simulator({self.natom} atoms, state={state})"
            )
        return f"Simulator(uninitialized, state={state})"
    
    # =========================================================================
    # Initialization and Cleanup
    # =========================================================================
    
    def init_simulation(self) -> None:
        """
        Initialize UppASD simulation.
        
        Sets up all data structures and allocates arrays.
        Must be called once before any measurements or relaxations.
        
        Raises
        ------
        StateError
            If already initialized
        RuntimeError
            If Fortran setup fails
        
        Examples
        --------
        >>> sim = Simulator()
        >>> sim.init_simulation()
        >>> print(f"System: {sim.natom} atoms, {sim.mensemble} ensembles")
        """
        self._check_state(SimulationState.UNINITIALIZED)
        
        try:
            logger.info("Initializing UppASD simulation")
            pyasd.sanity_check()
            self.natom, self.mensemble = pyasd.setup_all()
            self._state = SimulationState.INITIALIZED
            logger.info(
                f"✓ Initialized: {self.natom} atoms, {self.mensemble} ensembles"
            )
        except Exception as e:
            logger.error(f"Initialization failed: {e}")
            raise
    
    def cleanup(self) -> None:
        """
        Clean up resources and deallocate arrays.
        
        Safe to call multiple times. Silently does nothing if already cleaned.
        
        Examples
        --------
        >>> sim = Simulator()
        >>> try:
        ...     sim.init_simulation()
        ...     # Do work
        ... finally:
        ...     sim.cleanup()
        """
        if self._state == SimulationState.CLEANED_UP:
            logger.debug("Already cleaned up")
            return
        
        try:
            logger.debug("Cleaning up simulation")
            pyasd.cleanup()
            self._state = SimulationState.CLEANED_UP
            logger.debug("✓ Cleanup complete")
        except Exception as e:
            logger.warning(f"Cleanup issue: {e}")
    
    # =========================================================================
    # Simulation Execution
    # =========================================================================
    
    def run_simulation(self) -> np.ndarray:
        """
        Execute full simulation pipeline.
        
        Runs: sanity check → setup → initial phase → measure → cleanup.
        
        Returns
        -------
        moments : ndarray
            Final spin moments. Shape: (3, natom, mensemble)
        
        Raises
        ------
        StateError
            If not initialized
        RuntimeError
            If simulation fails
        
        Examples
        --------
        >>> with Simulator() as sim:
        ...     moments = sim.run_simulation()
        ...     energy = sim.energy
        """
        self._check_state(SimulationState.INITIALIZED)
        
        try:
            logger.info("Running full simulation pipeline")
            self._state = SimulationState.RUNNING
            
            pyasd.initial_phase()
            pyasd.measure()
            
            moments = self.moments
            self._state = SimulationState.COMPLETED
            logger.info("✓ Simulation completed")
            return moments
        except Exception as e:
            logger.error(f"Simulation failed: {e}")
            raise RuntimeError(f"Simulation failed: {e}") from e
    
    def relax(
        self,
        mode: str = "M",
        temperature: float = 0.0,
        steps: int = 100,
        timestep: float = 1.0e-16,
        damping: float = 0.5,
        callback: Optional[Callable[[int, np.ndarray, float], None]] = None,
    ) -> np.ndarray:
        """
        Perform spin relaxation using specified algorithm.
        
        Relaxes the spin system toward equilibrium (or energy minimum if T=0)
        using the chosen method. Optionally records a trajectory.
        
        Parameters
        ----------
        mode : {"M", "S", "H"}, optional
            Relaxation method:
            - "M": Metropolis Monte Carlo (default, good for finding minima)
            - "S": Landau-Lifshitz-Gilbert / Spin dynamics
            - "H": Heat Bath
            Default: "M"
        temperature : float, optional
            Temperature in Kelvin. Default: 0.0 (ground state)
        steps : int, optional
            Number of relaxation steps. Default: 100
        timestep : float, optional
            Integration timestep (LLG only). Default: 1.0e-16 s
        damping : float, optional
            LLG damping parameter (0 to 1). Default: 0.5
        callback : callable, optional
            Function called after each step with signature:
            ``callback(step: int, moments: ndarray, energy: float)``
            Useful for real-time monitoring in notebooks.
            Default: None
        
        Returns
        -------
        moments : ndarray
            Final magnetic moments after relaxation.
            Shape: (3, natom, mensemble)
        
        Raises
        ------
        StateError
            If not initialized
        ValueError
            If parameters invalid
        RuntimeError
            If relaxation fails
        
        Examples
        --------
        **Simple relaxation:**
        
        >>> with Simulator() as sim:
        ...     moments = sim.relax(mode='M', temperature=100, steps=1000)
        
        **With callback for monitoring:**
        
        >>> def monitor(step, moments, energy):
        ...     if step % 100 == 0:
        ...         m = np.linalg.norm(np.mean(moments, axis=1))
        ...         print(f"Step {step}: |M| = {m:.3f}, E = {energy:.3f}")
        >>> 
        >>> with Simulator() as sim:
        ...     sim.relax(mode='M', steps=1000, callback=monitor)
        
        **Ground state search (T=0):**
        
        >>> with Simulator() as sim:
        ...     sim.relax(mode='M', temperature=0.0, steps=500)
        ...     print(f"Min energy: {sim.energy}")
        """
        # State check: can relax from INITIALIZED or COMPLETED
        self._check_state(SimulationState.INITIALIZED, SimulationState.COMPLETED)
        
        # Parameter validation
        if mode not in {'M', 'S', 'H'}:
            raise ValueError(
                f"Invalid mode '{mode}'. Must be 'M' (Metropolis), "
                f"'S' (spin dynamics), or 'H' (heat bath)"
            )
        
        if temperature < 0:
            raise ValueError(f"Temperature must be >= 0, got {temperature}")
        
        if steps < 1:
            raise ValueError(f"Steps must be >= 1, got {steps}")
        
        try:
            mode_names = {'M': 'Metropolis', 'S': 'LLG', 'H': 'Heat Bath'}
            logger.info(
                f"Starting {mode_names[mode]} relaxation: "
                f"T={temperature}K, {steps} steps"
            )
            self._state = SimulationState.RUNNING
            
            # Perform relaxation with optional callback
            for step in range(steps):
                moments = pyasd.relax(
                    self.natom,
                    self.mensemble,
                    mode=mode,
                    nstep=1,
                    temperature=temperature,
                    timestep=timestep,
                    damping=damping,
                )
                
                # Call user callback if provided
                if callback is not None:
                    try:
                        energy = pyasd.get_energy()
                        callback(step, moments, energy)
                    except Exception as e:
                        logger.warning(f"Callback error at step {step}: {e}")
                
                # Record trajectory if requested
                if self._record_trajectory:
                    self._trajectory.append(moments.copy())
            
            self._state = SimulationState.COMPLETED
            logger.info(f"✓ Relaxation completed ({steps} steps)")
            return moments
            
        except Exception as e:
            logger.error(f"Relaxation failed: {e}")
            raise RuntimeError(f"Relaxation failed: {e}") from e
    
    # =========================================================================
    # Properties: Magnetic Data
    # =========================================================================
    
    @property
    def moments(self) -> np.ndarray:
        """
        Get current magnetic moments.
        
        Returns
        -------
        moments : ndarray
            Magnetic moment vectors. Shape: (3, natom, mensemble)
            Each moment is a unit vector pointing in spin direction.
        
        Examples
        --------
        >>> with Simulator() as sim:
        ...     moments = sim.moments
        ...     print(moments.shape)
        ...     magnetization = np.linalg.norm(np.mean(moments, axis=1))
        """
        self._check_state(
            SimulationState.INITIALIZED,
            SimulationState.RUNNING,
            SimulationState.COMPLETED,
        )
        return pyasd.get_moments(self.natom, self.mensemble)
    
    @moments.setter
    def moments(self, value: np.ndarray) -> None:
        """
        Set magnetic moments.
        
        Parameters
        ----------
        value : ndarray
            Magnetic moments. Shape: (3, natom, mensemble)
        
        Examples
        --------
        >>> with Simulator() as sim:
        ...     new_moments = np.random.randn(3, sim.natom, sim.mensemble)
        ...     sim.moments = new_moments / np.linalg.norm(new_moments, axis=0, keepdims=True)
        """
        self._check_state(SimulationState.INITIALIZED, SimulationState.COMPLETED)
        pyasd.set_moments(value, self.natom, self.mensemble)
    
    @property
    def field(self) -> np.ndarray:
        """
        Get effective magnetic field.
        
        Returns
        -------
        field : ndarray
            Effective field (sum of all interactions).
            Shape: (3, natom, mensemble)
        """
        self._check_state(SimulationState.INITIALIZED, SimulationState.COMPLETED)
        return pyasd.get_field(self.natom, self.mensemble)
    
    @field.setter
    def field(self, value: np.ndarray) -> None:
        """Set effective magnetic field."""
        self._check_state(SimulationState.INITIALIZED, SimulationState.COMPLETED)
        pyasd.set_field(value, self.natom, self.mensemble)
    
    # =========================================================================
    # Properties: Observables
    # =========================================================================
    
    @property
    def energy(self) -> float:
        """
        Get total system energy.
        
        Returns
        -------
        energy : float
            Total energy in meV (or units set in UppASD)
        
        Examples
        --------
        >>> with Simulator() as sim:
        ...     sim.relax(steps=100)
        ...     print(f"Energy: {sim.energy:.4f} meV")
        """
        self._check_state(
            SimulationState.INITIALIZED,
            SimulationState.RUNNING,
            SimulationState.COMPLETED,
        )
        return pyasd.get_energy()
    
    @property
    def temperature(self) -> float:
        """
        Get current simulation temperature.
        
        Returns
        -------
        temperature : float
            Temperature in Kelvin
        
        Examples
        --------
        >>> with Simulator() as sim:
        ...     print(f"T = {sim.temperature} K")
        """
        return pyasd.get_temperature()
    
    @temperature.setter
    def temperature(self, value: float) -> None:
        """
        Set simulation temperature.
        
        Parameters
        ----------
        value : float
            Temperature in Kelvin (must be >= 0)
        
        Raises
        ------
        ValueError
            If temperature < 0
        
        Examples
        --------
        >>> with Simulator() as sim:
        ...     sim.temperature = 300  # Room temperature
        """
        if value < 0:
            raise ValueError(f"Temperature must be >= 0, got {value}")
        pyasd.set_temperature(value)
    
    @property
    def timestep(self) -> float:
        """
        Get integration timestep.
        
        Returns
        -------
        timestep : float
            Timestep in seconds
        """
        return pyasd.get_timestep()
    
    @property
    def nstep(self) -> int:
        """
        Get current step number.
        
        Returns
        -------
        nstep : int
            Number of completed steps
        """
        return pyasd.get_nstep()
    
    # =========================================================================
    # Trajectory and Analysis
    # =========================================================================
    
    def get_trajectory(self) -> np.ndarray:
        """
        Get recorded trajectory of spin moments.
        
        Only contains data if record_trajectory=True in __init__.
        
        Returns
        -------
        trajectory : ndarray
            Shape: (nsteps, 3, natom, mensemble) if recorded
            Empty array if no trajectory recorded
        
        Raises
        ------
        UserWarning
            If no trajectory data available
        
        Examples
        --------
        >>> sim = Simulator(record_trajectory=True)
        >>> with sim:
        ...     sim.relax(steps=100)
        ...     traj = sim.get_trajectory()
        ...     print(traj.shape)  # (100, 3, natom, mensemble)
        """
        if not self._trajectory:
            logger.warning("No trajectory recorded. Set record_trajectory=True")
            return np.array([])
        
        return np.array(self._trajectory)
    
    def reset_trajectory(self) -> None:
        """
        Clear recorded trajectory data.
        
        Useful between multiple relaxations if recording is enabled.
        
        Examples
        --------
        >>> sim = Simulator(record_trajectory=True)
        >>> with sim:
        ...     sim.relax(steps=100)
        ...     traj1 = sim.get_trajectory()
        ...     sim.reset_trajectory()
        ...     sim.relax(steps=100)
        ...     traj2 = sim.get_trajectory()
        """
        self._trajectory = []
        self._trajectory_metadata = {}
        logger.debug("Trajectory cleared")
    
    def compute_magnetization(self) -> float:
        """
        Compute average magnetization magnitude.
        
        Returns
        -------
        magnetization : float
            Mean magnitude of magnetic moments (0 to 1)
        
        Examples
        --------
        >>> with Simulator() as sim:
        ...     sim.relax(mode='M', steps=1000)
        ...     m = sim.compute_magnetization()
        ...     print(f"Magnetization: {m:.3f}")
        """
        moments = self.moments
        # Average over atoms and ensembles, return magnitude
        return float(np.mean(np.linalg.norm(moments, axis=0)))
    
    def compute_susceptibility(self) -> float:
        """
        Compute magnetic susceptibility (simple estimate).
        
        Computed as: χ ≈ (∂M/∂T) by finite difference
        
        Returns
        -------
        susceptibility : float
            Simple magnetic susceptibility estimate
        
        Examples
        --------
        >>> with Simulator() as sim:
        ...     chi = sim.compute_susceptibility()
        """
        # Store current state
        T0 = self.temperature
        m0 = self.compute_magnetization()
        
        # Perturb temperature slightly
        dT = 0.1  # K
        self.temperature = T0 + dT
        m1 = self.compute_magnetization()
        
        # Restore original temperature
        self.temperature = T0
        
        # Compute finite difference
        chi = (m1 - m0) / dT if dT != 0 else 0.0
        
        return float(chi)
