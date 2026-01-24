"""
Low-level Python wrapper for UppASD Fortran C extension (_uppasd).

This module provides functions for running UppASD simulations, performing relaxation
methods, and accessing measured quantities. All functions wrap corresponding Fortran
routines exposed via F2PY.

Functions are organized into three categories:

1. **Simulation Control**
   - `run_uppasd()`: Execute full simulation
   - `sanity_check()`: Validate configuration
   - `setup_all()`: Initialize all data structures
   - `initial_phase()`: Run initialization phase
   - `measure()`: Perform measurement
   - `cleanup()`: Deallocate resources

2. **Relaxation Methods**
   - `relax_montecarlo()`: MC relaxation
   - `relax_metropolis()`: Metropolis relaxation
   - `relax_heatbath()`: Heat bath relaxation
   - `relax_llg()`: LLG relaxation
   - `relax()`: General relaxation dispatcher

3. **Data Access**
   - Magnetic moments: `get_moments()`, `set_moments()`
   - Fields: `get_field()`, `set_field()`
   - Energy: `get_energy()`
   - Temperature: `get_temperature()`, `set_temperature()`
   - Timestep: `get_timestep()`
   - Step counters: `get_nstep()`, `get_mc_nstep()`

All functions include comprehensive error handling, logging, and NaN detection.

Examples
--------
>>> import uppasd.pyasd as asd
>>> asd.sanity_check()
>>> natom, mensemble = asd.setup_all()
>>> asd.initial_phase()
>>> asd.measure()
>>> moments = asd.get_moments(natom, mensemble)
>>> asd.cleanup()

Legacy function names (snake_case aliases) are provided for backward compatibility:
- `runuppasd()` → `run_uppasd()`
- `setupall()` → `setup_all()`
- `relaxmontecarlo()` → `relax_montecarlo()`
- etc.
"""

from __future__ import absolute_import, division, print_function

import logging
from typing import Tuple, Optional

import numpy as np

import _uppasd


# Configure module logger
logger = logging.getLogger(__name__)


def _check_nan(value: float, name: str = "value") -> None:
    """
    Check if a scalar value is NaN and log warning if so.
    
    Parameters
    ----------
    value : float
        Value to check
    name : str
        Name of the value for logging
    """
    if isinstance(value, (float, np.floating)) and np.isnan(value):
        logger.warning(f"NaN detected in {name}")


def _check_array_nan(array: np.ndarray, name: str = "array") -> None:
    """
    Check if an array contains NaN values and log warning if so.
    
    Parameters
    ----------
    array : ndarray
        Array to check
    name : str
        Name of the array for logging
    """
    if np.any(np.isnan(array)):
        n_nan = np.sum(np.isnan(array))
        logger.warning(f"NaN detected in {name}: {n_nan}/{array.size} elements")


# ============================================================================
# Simulation Control Routines
# ============================================================================

def run_uppasd() -> None:
    """
    Execute the UppASD simulation.
    
    Runs the main simulation loop with parameters configured via setup_all()
    and initial_phase().
    
    Raises
    ------
    RuntimeError
        If Fortran routine fails
    
    Examples
    --------
    >>> sanity_check()
    >>> setup_all()
    >>> initial_phase()
    >>> run_uppasd()
    >>> measure()
    """
    try:
        logger.info("Starting UppASD simulation")
        _uppasd.runuppasd()
        logger.debug("✓ Simulation completed")
    except Exception as e:
        logger.error(f"Failed to run UppASD: {e}")
        raise RuntimeError(f"UppASD simulation failed: {e}") from e


# Backward compatibility alias
runuppasd = run_uppasd


def sanity_check() -> None:
    """
    Perform a sanity check on the UppASD configuration.
    
    Validates that all required parameters are set and consistent.
    Should be called before setup_all().
    
    Raises
    ------
    RuntimeError
        If configuration is invalid
    
    Examples
    --------
    >>> sanity_check()  # Returns silently if OK, raises if not
    """
    try:
        logger.debug("Running sanity check")
        _uppasd.sanitycheck()
        logger.debug("✓ Sanity check passed")
    except Exception as e:
        logger.error(f"Sanity check failed: {e}")
        raise RuntimeError(f"Configuration invalid: {e}") from e


# Backward compatibility alias
sanitycheck = sanity_check


def num_procs() -> int:
    """
    Get the number of active OpenMP processors.
    
    Returns
    -------
    nprocs : int
        Number of processors used by UppASD
    
    Examples
    --------
    >>> n = num_procs()
    >>> print(f"Using {n} processors")
    """
    try:
        logger.debug("Fetching processor count")
        numprocs_fn = getattr(_uppasd, "numprocs", None)
        if numprocs_fn is None:
            logger.warning("numprocs not available in _uppasd; assuming 1")
            return 1

        nprocs = numprocs_fn()
        if not isinstance(nprocs, (int, np.integer)):
            logger.warning(f"numprocs returned unexpected type: {type(nprocs)}")
        return int(nprocs)
    except Exception as e:
        logger.error(f"Failed to get processor count: {e}")
        raise RuntimeError(f"Cannot retrieve processor count: {e}") from e


# Backward compatibility alias
numprocs = num_procs


def print_logo() -> None:
    """
    Print UppASD logo and version information.
    
    Calls Fortran routine to display startup banner.
    
    Examples
    --------
    >>> print_logo()
    # Prints UppASD ASCII art and version info to stdout
    """
    try:
        logger.debug("Printing UppASD logo")
        _uppasd.printlogo()
    except Exception as e:
        logger.warning(f"Could not print logo: {e}")


# Backward compatibility alias
printlogo = print_logo


def setup_all() -> Tuple[int, int]:
    """
    Initialize all UppASD data structures.
    
    Allocates and initializes arrays for atoms, moments, fields, etc.
    Must be called once before any measurements or relaxations.
    Should be called after reading input files.
    
    Returns
    -------
    natom : int
        Number of atoms in the system
    mensemble : int
        Number of independent ensembles (usually 1)
    
    Raises
    ------
    RuntimeError
        If initialization fails
    
    Examples
    --------
    >>> sanity_check()
    >>> natom, mensemble = setup_all()
    >>> print(f"System: {natom} atoms, {mensemble} ensembles")
    """
    try:
        logger.info("Setting up all UppASD data structures")
        natom, mensemble = _uppasd.setupall()
        natom = int(natom)
        mensemble = int(mensemble)
        logger.info(f"✓ Setup complete: {natom} atoms, {mensemble} ensembles")
        return natom, mensemble
    except Exception as e:
        logger.error(f"Failed to setup simulation: {e}")
        raise RuntimeError(f"Initialization failed: {e}") from e


# Backward compatibility alias
setupall = setup_all


def initial_phase() -> None:
    """
    Perform the initial phase of the UppASD simulation.
    
    Runs initialization including thermalization and initial measurements.
    Must be called after setup_all() and before run_uppasd().
    
    Raises
    ------
    RuntimeError
        If initialization fails
    
    Examples
    --------
    >>> setup_all()
    >>> initial_phase()
    >>> run_uppasd()
    """
    try:
        logger.info("Running initial phase")
        _uppasd.initialphase()
        logger.debug("✓ Initial phase completed")
    except Exception as e:
        logger.error(f"Initial phase failed: {e}")
        raise RuntimeError(f"Initialization phase failed: {e}") from e


# Backward compatibility alias
initialphase = initial_phase


def measure() -> None:
    """
    Perform measurement step of the UppASD simulation.
    
    Records observables (magnetization, energy, etc.) at current timestep.
    Typically called after each integration step.
    
    Raises
    ------
    RuntimeError
        If measurement fails
    
    Examples
    --------
    >>> initial_phase()
    >>> for step in range(n_steps):
    ...     # Perform integration step
    ...     measure()  # Record observables
    """
    try:
        logger.debug("Performing measurement")
        _uppasd.measure()
    except Exception as e:
        logger.error(f"Measurement failed: {e}")
        raise RuntimeError(f"Cannot perform measurement: {e}") from e


def cleanup() -> None:
    """
    Clean up UppASD resources.
    
    Deallocates all arrays and closes output files.
    Should be called at end of simulation.
    Safe to call multiple times.
    
    Examples
    --------
    >>> try:
    ...     setup_all()
    ...     initial_phase()
    ...     run_uppasd()
    ... finally:
    ...     cleanup()  # Always called, even on error
    """
    try:
        logger.debug("Cleaning up UppASD")
        _uppasd.cleanup()
        logger.debug("✓ Cleanup completed")
    except Exception as e:
        logger.warning(f"Cleanup encountered issue: {e}")
        # Don't raise — cleanup should not fail the program


# ============================================================================
# Relaxation Methods
# ============================================================================

def relax_montecarlo(natom: int, mensemble: int) -> np.ndarray:
    """
    Perform Monte Carlo relaxation.
    
    Relaxes spin system using Monte Carlo dynamics toward ground state
    or thermal equilibrium.
    
    Parameters
    ----------
    natom : int
        Number of atoms
    mensemble : int
        Number of ensembles
    
    Returns
    -------
    moments : ndarray
        Final spin moments. Shape: (3, natom, mensemble)
    
    Raises
    ------
    RuntimeError
        If relaxation fails or returns invalid data
    
    Examples
    --------
    >>> natom, mensemble = setup_all()
    >>> initial_phase()
    >>> moments = relax_montecarlo(natom, mensemble)
    """
    try:
        logger.info(f"Starting MC relaxation ({natom} atoms, {mensemble} ens.)")
        moments = _uppasd.relaxmontecarlo(natom, mensemble)
        moments = np.array(moments, dtype=np.float64, copy=True)
        _check_array_nan(moments, "MC relaxation moments")
        logger.debug(f"✓ MC relaxation complete, moments shape: {moments.shape}")
        return moments
    except Exception as e:
        logger.error(f"MC relaxation failed: {e}")
        raise RuntimeError(f"Monte Carlo relaxation failed: {e}") from e


# Backward compatibility alias
relaxmontecarlo = relax_montecarlo


def relax_metropolis(natom: int, mensemble: int) -> np.ndarray:
    """
    Perform Metropolis relaxation.
    
    Relaxes spin system using Metropolis algorithm with single-spin
    flip dynamics.
    
    Parameters
    ----------
    natom : int
        Number of atoms
    mensemble : int
        Number of ensembles
    
    Returns
    -------
    moments : ndarray
        Final spin moments. Shape: (3, natom, mensemble)
    
    Raises
    ------
    RuntimeError
        If relaxation fails or returns invalid data
    
    Examples
    --------
    >>> natom, mensemble = setup_all()
    >>> initial_phase()
    >>> moments = relax_metropolis(natom, mensemble)
    """
    try:
        logger.info(f"Starting Metropolis relaxation ({natom} atoms)")
        moments = _uppasd.relaxmetropolis(natom, mensemble)
        moments = np.array(moments, dtype=np.float64, copy=True)
        _check_array_nan(moments, "Metropolis relaxation moments")
        logger.debug(f"✓ Metropolis relaxation complete")
        return moments
    except Exception as e:
        logger.error(f"Metropolis relaxation failed: {e}")
        raise RuntimeError(f"Metropolis relaxation failed: {e}") from e


# Backward compatibility alias
relaxmetropolis = relax_metropolis


def relax_heatbath(natom: int, mensemble: int) -> np.ndarray:
    """
    Perform heat bath relaxation.
    
    Relaxes spin system using heat bath dynamics.
    
    Parameters
    ----------
    natom : int
        Number of atoms
    mensemble : int
        Number of ensembles
    
    Returns
    -------
    moments : ndarray
        Final spin moments. Shape: (3, natom, mensemble)
    
    Raises
    ------
    RuntimeError
        If relaxation fails or returns invalid data
    
    Examples
    --------
    >>> natom, mensemble = setup_all()
    >>> initial_phase()
    >>> moments = relax_heatbath(natom, mensemble)
    """
    try:
        logger.info(f"Starting heat bath relaxation ({natom} atoms)")
        moments = _uppasd.relaxheatbath(natom, mensemble)
        moments = np.array(moments, dtype=np.float64, copy=True)
        _check_array_nan(moments, "Heat bath moments")
        logger.debug(f"✓ Heat bath relaxation complete")
        return moments
    except Exception as e:
        logger.error(f"Heat bath relaxation failed: {e}")
        raise RuntimeError(f"Heat bath relaxation failed: {e}") from e


# Backward compatibility alias
relaxheatbath = relax_heatbath


def relax_llg(natom: int, mensemble: int) -> np.ndarray:
    """
    Perform LLG (Landau-Lifshitz-Gilbert) relaxation.
    
    Relaxes spin system using micromagnetic Landau-Lifshitz-Gilbert
    equation with damping.
    
    Parameters
    ----------
    natom : int
        Number of atoms
    mensemble : int
        Number of ensembles
    
    Returns
    -------
    moments : ndarray
        Final spin moments. Shape: (3, natom, mensemble)
    
    Raises
    ------
    RuntimeError
        If relaxation fails or returns invalid data
    
    Examples
    --------
    >>> natom, mensemble = setup_all()
    >>> initial_phase()
    >>> moments = relax_llg(natom, mensemble)
    """
    try:
        logger.info(f"Starting LLG relaxation ({natom} atoms)")
        moments = _uppasd.relaxllg(natom, mensemble)
        moments = np.array(moments, dtype=np.float64, copy=True)
        _check_array_nan(moments, "LLG relaxation moments")
        logger.debug(f"✓ LLG relaxation complete")
        return moments
    except Exception as e:
        logger.error(f"LLG relaxation failed: {e}")
        raise RuntimeError(f"LLG relaxation failed: {e}") from e


# Backward compatibility alias
relaxllg = relax_llg


def relax(
    natom: int,
    mensemble: int,
    mode: str = "S",
    nstep: int = 10,
    temperature: float = 0.0,
    timestep: float = 1.0e-16,
    damping: float = 0.5,
) -> np.ndarray:
    """
    Perform relaxation with specified method and parameters.
    
    General-purpose relaxation dispatcher that routes to the appropriate
    method based on ``mode`` parameter.
    
    Parameters
    ----------
    natom : int
        Number of atoms
    mensemble : int
        Number of ensembles
    mode : {'S', 'M', 'H'}, optional
        Relaxation method:
        - 'S': Spin-dynamics (LLG, default)
        - 'M': Metropolis (single-spin flip MC)
        - 'H': Heat bath
        Default: 'S'
    nstep : int, optional
        Number of relaxation steps. Default: 10
    temperature : float, optional
        Temperature in Kelvin. Default: 0.0 (ground state)
    timestep : float, optional
        Integration timestep in seconds. Default: 1.0e-16
    damping : float, optional
        Damping factor (0 to 1). Default: 0.5
    
    Returns
    -------
    moments : ndarray
        Final spin moments. Shape: (3, natom, mensemble)
    
    Raises
    ------
    ValueError
        If mode is invalid
    RuntimeError
        If relaxation fails or returns invalid data
    
    Examples
    --------
    **LLG relaxation (spin dynamics):**
    
    >>> natom, mensemble = setup_all()
    >>> initial_phase()
    >>> moments = relax(
    ...     natom, mensemble,
    ...     mode='S',
    ...     nstep=1000,
    ...     temperature=100,
    ...     timestep=1e-15,
    ...     damping=0.1
    ... )
    
    **Metropolis relaxation (MC):**
    
    >>> moments = relax(natom, mensemble, mode='M', temperature=300)
    
    **Heat bath:**
    
    >>> moments = relax(natom, mensemble, mode='H', temperature=100)
    """
    # Validate mode
    valid_modes = {'S', 'M', 'H'}
    if mode not in valid_modes:
        raise ValueError(
            f"Invalid relaxation mode '{mode}'. "
            f"Must be one of {valid_modes}"
        )
    
    method_names = {'S': 'LLG', 'M': 'Metropolis', 'H': 'Heat Bath'}
    method = method_names[mode]
    
    try:
        logger.info(
            f"Starting {method} relaxation: "
            f"{nstep} steps, T={temperature}K, "
            f"dt={timestep}s, α={damping}"
        )

        # Keep ip_mode aligned with explicit relaxation mode when possible
        if hasattr(_uppasd, "put_ipmode"):
            try:
                set_ipmode(mode)
            except Exception as sync_err:
                logger.debug(f"ip_mode sync skipped: {sync_err}")
        
        moments = _uppasd.relax(
            natom, mensemble, mode, nstep, temperature, timestep, damping
        )
        # Wrap Fortran array as numpy array with proper dtype
        # F2PY allocates this array, we just reference it
        moments = np.asarray(moments, dtype=np.float64)
        _check_array_nan(moments, f"{method} relaxation moments")
        logger.info(f"✓ {method} relaxation complete")
        return moments
        
    except Exception as e:
        logger.error(f"{method} relaxation failed: {e}")
        raise RuntimeError(f"Relaxation failed ({method}): {e}") from e

# ============================================================================
# Data Access: Magnetic Moments
# ============================================================================

def get_coords(natom: int) -> np.ndarray:
    """
    Get atomic coordinates.
    
    Parameters
    ----------
    natom : int
        Number of atoms
    
    Returns
    -------
    coords : ndarray
        Atomic coordinates. Shape: (3, natom)
    
    Examples
    --------
    >>> natom, _ = setup_all()
    >>> coords = get_coords(natom)
    >>> print(coords.shape)
    (3, N)
    """
    try:
        logger.debug("Fetching atomic coordinates")
        coords = _uppasd.get_coord(natom)
        coords = np.array(coords, dtype=np.float64, copy=True)
        _check_array_nan(coords, "coordinates")
        return coords
    except Exception as e:
        logger.error(f"Failed to get coordinates: {e}")
        raise RuntimeError(f"Cannot retrieve coordinates: {e}") from e


def get_moments(natom: int, mensemble: int) -> np.ndarray:
    """
    Get effective magnetic moments.
    
    Retrieves the magnetic moment vectors for all atoms and ensembles.
    
    Parameters
    ----------
    natom : int
        Number of atoms
    mensemble : int
        Number of ensembles
    
    Returns
    -------
    moments : ndarray
        Magnetic moments. Shape: (3, natom, mensemble)
        Each moment is a unit vector.
    
    Raises
    ------
    RuntimeError
        If retrieval fails or returns invalid data
    
    Examples
    --------
    >>> natom, mensemble = setup_all()
    >>> moments = get_moments(natom, mensemble)
    >>> print(f"Moments shape: {moments.shape}")
    >>> magnetization = np.mean(np.linalg.norm(moments, axis=0))
    """
    try:
        logger.debug("Fetching magnetic moments")
        moments = _uppasd.get_emom(natom, mensemble)
        # CRITICAL: Must use explicit copy, not asarray()
        # Fortran arrays passed to Python can cause malloc corruption if not copied
        # This prevents the "Incorrect checksum for freed object" error
        moments = np.array(moments, dtype=np.float64, copy=True)
        _check_array_nan(moments, "moments")
        return moments
    except Exception as e:
        logger.error(f"Failed to get moments: {e}")
        raise RuntimeError(f"Cannot retrieve moments: {e}") from e


# Backward compatibility alias
get_emom = get_moments


def set_moments(moments: np.ndarray, natom: int, mensemble: int) -> None:
    """
    Set effective magnetic moments.
    
    Sets the magnetic moment vectors for all atoms and ensembles.
    
    Parameters
    ----------
    moments : ndarray
        Magnetic moments. Shape: (3, natom, mensemble)
    natom : int
        Number of atoms
    mensemble : int
        Number of ensembles
    
    Raises
    ------
    RuntimeError
        If setting fails
    
    Examples
    --------
    >>> natom, mensemble = setup_all()
    >>> moments = np.random.randn(3, natom, mensemble)
    >>> moments /= np.linalg.norm(moments, axis=0, keepdims=True)  # Normalize
    >>> set_moments(moments, natom, mensemble)
    """
    try:
        moments = np.asarray(moments, dtype=np.float64)
        logger.debug(f"Setting moments, shape: {moments.shape}")
        _uppasd.put_emom(moments, natom, mensemble)
    except Exception as e:
        logger.error(f"Failed to set moments: {e}")
        raise RuntimeError(f"Cannot set moments: {e}") from e


# Backward compatibility alias
put_emom = set_moments


# ============================================================================
# Data Access: Fields
# ============================================================================

def get_field(natom: int, mensemble: int) -> np.ndarray:
    """
    Get effective magnetic field.
    
    Retrieves the effective field (sum of all interactions) for all atoms.
    
    Parameters
    ----------
    natom : int
        Number of atoms
    mensemble : int
        Number of ensembles
    
    Returns
    -------
    field : ndarray
        Effective field. Shape: (3, natom, mensemble)
    
    Raises
    ------
    RuntimeError
        If retrieval fails or returns invalid data
    
    Examples
    --------
    >>> natom, mensemble = setup_all()
    >>> field = get_field(natom, mensemble)
    >>> field_magnitude = np.linalg.norm(field, axis=0)
    """
    try:
        logger.debug("Fetching effective field")
        field = _uppasd.get_beff(natom, mensemble)
        field = np.array(field, dtype=np.float64, copy=True)
        _check_array_nan(field, "field")
        return field
    except Exception as e:
        logger.error(f"Failed to get field: {e}")
        raise RuntimeError(f"Cannot retrieve field: {e}") from e


# Backward compatibility alias
get_beff = get_field


def set_field(field: np.ndarray, natom: int, mensemble: int) -> None:
    """
    Set effective magnetic field.
    
    Sets the external magnetic field applied to the system.
    
    Parameters
    ----------
    field : ndarray
        Applied field. Shape: (3, natom, mensemble) or (3,) for uniform field
    natom : int
        Number of atoms
    mensemble : int
        Number of ensembles
    
    Raises
    ------
    RuntimeError
        If setting fails
    
    Examples
    --------
    >>> natom, mensemble = setup_all()
    >>> field = np.zeros((3, natom, mensemble))
    >>> field[2, :, :] = 1.0  # 1 Tesla in z-direction
    >>> set_field(field, natom, mensemble)
    """
    try:
        field = np.asarray(field, dtype=np.float64)
        logger.debug(f"Setting field, shape: {field.shape}")
        _uppasd.put_beff(field, natom, mensemble)
    except Exception as e:
        logger.error(f"Failed to set field: {e}")
        raise RuntimeError(f"Cannot set field: {e}") from e


# Backward compatibility alias
put_beff = set_field


def get_hfield() -> np.ndarray:
    """
    Get external magnetic field.
    
    Retrieves the applied external magnetic field.
    
    Returns
    -------
    hfield : ndarray
        External field. Shape: (3,) representing (Bx, By, Bz) in Tesla
    
    Raises
    ------
    RuntimeError
        If retrieval fails
    
    Examples
    --------
    >>> setup_all()
    >>> hfield = get_hfield()
    >>> print(f"Applied field: {hfield} T")
    """
    try:
        logger.debug("Fetching external magnetic field")
        hfield = _uppasd.get_hfield()
        hfield = np.array(hfield, dtype=np.float64, copy=True)
        _check_array_nan(hfield, "hfield")
        return hfield
    except Exception as e:
        logger.error(f"Failed to get hfield: {e}")
        raise RuntimeError(f"Cannot retrieve hfield: {e}") from e


def set_hfield(hfield: np.ndarray) -> None:
    """
    Set external magnetic field.
    
    Updates the applied external magnetic field.
    
    Parameters
    ----------
    hfield : ndarray
        External field. Shape: (3,) representing (Bx, By, Bz) in Tesla
    
    Raises
    ------
    ValueError
        If hfield has wrong shape or contains NaN
    RuntimeError
        If setting fails
    
    Examples
    --------
    >>> setup_all()
    >>> set_hfield(np.array([0.0, 0.0, 1.0]))  # 1 Tesla in z-direction
    """
    try:
        hfield = np.asarray(hfield, dtype=np.float64)
        if hfield.shape != (3,):
            raise ValueError(f"hfield must have shape (3,), got {hfield.shape}")
        _check_array_nan(hfield, "hfield")
        logger.debug(f"Setting external field to {hfield} T")
        _uppasd.put_hfield(hfield)
    except Exception as e:
        logger.error(f"Failed to set hfield: {e}")
        raise RuntimeError(f"Cannot set hfield: {e}") from e


# Backward compatibility alias
put_hfield = set_hfield


def get_iphfield() -> np.ndarray:
    """
    Get initial-phase external magnetic field.
    
    Retrieves the applied external magnetic field for the initial phase.
    
    Returns
    -------
    iphfield : ndarray
        Initial-phase external field. Shape: (3,) representing (Bx, By, Bz) in Tesla
    
    Raises
    ------
    RuntimeError
        If retrieval fails
    
    Examples
    --------
    >>> setup_all()
    >>> iphfield = get_iphfield()
    >>> print(f"Initial phase field: {iphfield} T")
    """
    try:
        logger.debug("Fetching initial-phase external magnetic field")
        iphfield = _uppasd.get_iphfield()
        iphfield = np.array(iphfield, dtype=np.float64, copy=True)
        _check_array_nan(iphfield, "iphfield")
        return iphfield
    except Exception as e:
        logger.error(f"Failed to get iphfield: {e}")
        raise RuntimeError(f"Cannot retrieve iphfield: {e}") from e


def set_iphfield(iphfield: np.ndarray) -> None:
    """
    Set initial-phase external magnetic field.
    
    Updates the applied external magnetic field for the initial phase.
    
    Parameters
    ----------
    iphfield : ndarray
        Initial-phase external field. Shape: (3,) representing (Bx, By, Bz) in Tesla
    
    Raises
    ------
    ValueError
        If iphfield has wrong shape or contains NaN
    RuntimeError
        If setting fails
    
    Examples
    --------
    >>> setup_all()
    >>> set_iphfield(np.array([0.0, 0.0, 0.5]))  # 0.5 Tesla in z-direction
    """
    try:
        iphfield = np.asarray(iphfield, dtype=np.float64)
        if iphfield.shape != (3,):
            raise ValueError(f"iphfield must have shape (3,), got {iphfield.shape}")
        _check_array_nan(iphfield, "iphfield")
        logger.debug(f"Setting initial-phase field to {iphfield} T")
        _uppasd.put_iphfield(iphfield)
    except Exception as e:
        logger.error(f"Failed to set iphfield: {e}")
        raise RuntimeError(f"Cannot set iphfield: {e}") from e


# Backward compatibility alias
put_iphfield = set_iphfield


# ============================================================================
# Data Access: Energy and Observables
# ============================================================================

def get_energy() -> float:
    """
    Get total system energy.
    
    Returns the sum of all interaction energies (exchange, anisotropy,
    external field, etc.).
    
    Returns
    -------
    energy : float
        Total energy in units of meV (as set in UppASD)
    
    Raises
    ------
    RuntimeError
        If retrieval fails
    
    Examples
    --------
    >>> setup_all()
    >>> initial_phase()
    >>> energy = get_energy()
    >>> print(f"System energy: {energy:.4f} meV")
    """
    try:
        logger.debug("Fetching total energy")
        energy = _uppasd.get_energy()
        energy = float(energy)
        _check_nan(energy, "energy")
        return energy
    except Exception as e:
        logger.error(f"Failed to get energy: {e}")
        raise RuntimeError(f"Cannot retrieve energy: {e}") from e


def get_nstep() -> int:
    """
    Get current simulation step number.
    
    Returns the number of completed integration steps since initialization.
    
    Returns
    -------
    nstep : int
        Current step number
    
    Examples
    --------
    >>> setup_all()
    >>> for _ in range(1000):
    ...     step = get_nstep()
    ...     print(f"Step {step}")
    """
    try:
        logger.debug("Fetching current step number")
        nstep = _uppasd.get_nstep()
        return int(nstep)
    except Exception as e:
        logger.error(f"Failed to get step number: {e}")
        raise RuntimeError(f"Cannot retrieve step number: {e}") from e


def set_nstep(nstep: int) -> None:
    """
    Set simulation step number.
    
    Updates the internal step counter to the specified value.
    
    Parameters
    ----------
    nstep : int
        The new step number
    
    Raises
    ------
    ValueError
        If nstep is negative
    RuntimeError
        If setting fails
    
    Examples
    --------
    >>> setup_all()
    >>> set_nstep(1000)
    """
    if not isinstance(nstep, int):
        nstep = int(nstep)
    
    if nstep < 0:
        raise ValueError(f"nstep must be non-negative, got {nstep}")
    
    try:
        if hasattr(_uppasd, "put_nstep"):
            logger.debug(f"Setting step number to {nstep}")
            _uppasd.put_nstep(nstep)
        else:
            logger.warning("put_nstep not available in _uppasd; step number not updated")
    except Exception as e:
        logger.error(f"Failed to set step number: {e}")
        raise RuntimeError(f"Cannot set step number: {e}") from e


# Backward compatibility alias
put_nstep = set_nstep


def get_mcnstep() -> int:
    """
    Get current Monte Carlo step number.
    
    Returns the number of completed Monte Carlo steps since initialization.
    
    Returns
    -------
    mcnstep : int
        Current MC step number
    
    Raises
    ------
    RuntimeError
        If retrieval fails
    
    Examples
    --------
    >>> setup_all()
    >>> for _ in range(10000):
    ...     step = get_mcnstep()
    ...     print(f"MC Step {step}")
    """
    try:
        logger.debug("Fetching current MC step number")
        mcnstep = _uppasd.get_mcnstep()
        return int(mcnstep)
    except Exception as e:
        logger.error(f"Failed to get MC step number: {e}")
        raise RuntimeError(f"Cannot retrieve MC step number: {e}") from e


def set_mcnstep(mcnstep: int) -> None:
    """
    Set Monte Carlo step number.
    
    Updates the internal Monte Carlo step counter.
    
    Parameters
    ----------
    mcnstep : int
        The new MC step number
    
    Raises
    ------
    ValueError
        If mcnstep is negative
    RuntimeError
        If setting fails
    
    Examples
    --------
    >>> setup_all()
    >>> set_mcnstep(5000)
    """
    if not isinstance(mcnstep, int):
        mcnstep = int(mcnstep)
    
    if mcnstep < 0:
        raise ValueError(f"mcnstep must be non-negative, got {mcnstep}")
    
    try:
        if hasattr(_uppasd, "put_mcnstep"):
            logger.debug(f"Setting MC step number to {mcnstep}")
            _uppasd.put_mcnstep(mcnstep)
        else:
            logger.warning("put_mcnstep not available in _uppasd; MC step number not updated")
    except Exception as e:
        logger.error(f"Failed to set MC step number: {e}")
        raise RuntimeError(f"Cannot set MC step number: {e}") from e


# Backward compatibility alias
put_mcnstep = set_mcnstep


# ============================================================================
# Data Access: Input Mode
# ============================================================================

def get_ipmode() -> str:
    """
    Get current initial-phase mode (``ip_mode``).

    Returns
    -------
    mode : str
        Two-character mode specifier (e.g., ``S``, ``M``, ``H``, ``SX``, ``Q``).

    Raises
    ------
    RuntimeError
        If the underlying extension lacks ip_mode support or retrieval fails.
    """
    if not hasattr(_uppasd, "get_ipmode"):
        raise RuntimeError(
            "_uppasd.build missing get_ipmode(); rebuild UppASD to enable ip_mode access"
        )

    try:
        raw_mode = _uppasd.get_ipmode()
        arr = np.asarray(raw_mode, dtype="U1")
        mode = "".join(arr.tolist()).strip().upper()
        return mode
    except Exception as e:
        logger.error(f"Failed to get ip_mode: {e}")
        raise RuntimeError(f"Cannot retrieve ip_mode: {e}") from e


def set_ipmode(mode: str) -> None:
    """
    Set initial-phase mode (``ip_mode``).

    Parameters
    ----------
    mode : str
        One- or two-character mode specifier (e.g., ``S``, ``M``, ``H``, ``SX``, ``Q``, ``Y``).

    Raises
    ------
    ValueError
        If ``mode`` is empty or longer than two characters.
    RuntimeError
        If the underlying extension lacks ip_mode support or update fails.
    """
    if not hasattr(_uppasd, "put_ipmode"):
        raise RuntimeError(
            "_uppasd.build missing put_ipmode(); rebuild UppASD to enable ip_mode updates"
        )

    if not isinstance(mode, str) or len(mode.strip()) == 0:
        raise ValueError("ip_mode must be a non-empty string")

    normalized = mode.strip().upper()
    if len(normalized) > 2:
        raise ValueError(f"ip_mode must be 1-2 characters, got '{mode}'")

    padded = normalized.ljust(2)
    payload = np.frombuffer(padded.encode("ascii"), dtype="S1")

    try:
        _uppasd.put_ipmode(payload)
    except Exception as e:
        logger.error(f"Failed to set ip_mode to '{normalized}': {e}")
        raise RuntimeError(f"Cannot set ip_mode: {e}") from e


# Backward compatibility alias
put_ipmode = set_ipmode


# ============================================================================
# Data Access: Temperature
# ============================================================================

def get_temperature() -> float:
    """
    Get current simulation temperature.
    
    Returns
    -------
    temperature : float
        Temperature in Kelvin
    
    Examples
    --------
    >>> setup_all()
    >>> T = get_temperature()
    >>> print(f"Temperature: {T} K")
    """
    try:
        logger.debug("Fetching temperature")
        temperature = _uppasd.get_temperature()
        temperature = float(temperature)
        _check_nan(temperature, "temperature")
        return temperature
    except Exception as e:
        logger.error(f"Failed to get temperature: {e}")
        raise RuntimeError(f"Cannot retrieve temperature: {e}") from e


def set_temperature(temperature: float) -> None:
    """
    Set simulation temperature.
    
    Parameters
    ----------
    temperature : float
        Temperature in Kelvin (must be >= 0)
    
    Raises
    ------
    ValueError
        If temperature is negative
    RuntimeError
        If setting fails
    
    Examples
    --------
    >>> setup_all()
    >>> set_temperature(300)  # Set to room temperature
    """
    if temperature < 0:
        raise ValueError(f"Temperature must be non-negative, got {temperature}")
    
    try:
        logger.info(f"Setting temperature to {temperature} K")
        _uppasd.put_temperature(temperature)
    except Exception as e:
        logger.error(f"Failed to set temperature: {e}")
        raise RuntimeError(f"Cannot set temperature: {e}") from e


# Backward compatibility alias
put_temperature = set_temperature


def set_iptemperature(temperature: float) -> None:
    """
    Set initial-phase temperature.
    
    Updates the temperature for the initial phase relaxation.
    
    Parameters
    ----------
    temperature : float
        Temperature in Kelvin (must be >= 0)
    
    Raises
    ------
    ValueError
        If temperature is negative
    RuntimeError
        If setting fails or if put_iptemperature not available
    
    Examples
    --------
    >>> setup_all()
    >>> set_iptemperature(100)  # Set initial phase to 100K
    """
    if temperature < 0:
        raise ValueError(f"Temperature must be non-negative, got {temperature}")
    
    if not hasattr(_uppasd, "put_iptemperature"):
        raise RuntimeError(
            "_uppasd.build missing put_iptemperature(); rebuild UppASD to enable iptemperature updates"
        )
    
    try:
        logger.debug(f"Setting initial-phase temperature to {temperature} K")
        _uppasd.put_iptemperature(temperature)
    except Exception as e:
        logger.error(f"Failed to set initial-phase temperature: {e}")
        raise RuntimeError(f"Cannot set initial-phase temperature: {e}") from e


# Backward compatibility alias
put_iptemperature = set_iptemperature


def get_timestep() -> float:
    """
    Get current integration timestep.
    
    Returns
    -------
    timestep : float
        Integration timestep in seconds
    
    Examples
    --------
    >>> setup_all()
    >>> dt = get_timestep()
    >>> print(f"Integration timestep: {dt:.2e} s")
    """
    try:
        logger.debug("Fetching timestep")
        timestep = _uppasd.get_delta_t()
        timestep = float(timestep)
        _check_nan(timestep, "timestep")
        return timestep
    except Exception as e:
        logger.error(f"Failed to get timestep: {e}")
        raise RuntimeError(f"Cannot retrieve timestep: {e}") from e


def set_timestep(timestep: float) -> None:
    """
    Set integration timestep.
    
    Parameters
    ----------
    timestep : float
        Integration timestep in seconds
    
    Raises
    ------
    ValueError
        If timestep is non-positive or NaN
    RuntimeError
        If setting fails
    
    Examples
    --------
    >>> setup_all()
    >>> set_timestep(1.0e-15)
    """
    timestep = float(timestep)
    _check_nan(timestep, "timestep")
    
    if timestep <= 0.0:
        raise ValueError(f"timestep must be positive, got {timestep}")
    
    try:
        if hasattr(_uppasd, "put_delta_t"):
            logger.debug(f"Setting timestep to {timestep:.2e}")
            _uppasd.put_delta_t(timestep)
        else:
            logger.warning("put_delta_t not available in _uppasd; timestep not updated")
    except Exception as e:
        logger.error(f"Failed to set timestep: {e}")
        raise RuntimeError(f"Cannot set timestep: {e}") from e


# Backward compatibility aliases
get_delta_t = get_timestep
set_delta_t = set_timestep
put_delta_t = set_timestep
