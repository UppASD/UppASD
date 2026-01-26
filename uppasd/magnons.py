"""
Linear Spin Wave Theory (LSWT) magnon calculations for UppASD.

This module provides Python interfaces to compute magnon dispersions, eigenvalues,
and related properties using linear spin wave theory as implemented in the Fortran
diamag module.

Features:
- Setup and compute magnon Hamiltonians from exchange tensors
- Extract magnon eigenvalues and eigenvectors
- Support for different magnetic structures via q-vectors
- Future integration with spglib/seekpath for automatic q-point generation

Key Functions:
- `compute_magnons()`: Main user-facing function
- `setup_q_mesh()`: Create reciprocal space q-point mesh
- `get_magnon_dispersion()`: Extract dispersion relation

Architecture:
- Low-level Fortran interface in pyasd module (via f90wrap)
- Magnon-specific Python layer with validation and convenience functions
- Notebook helpers for interactive use (in notebook.py)

Examples
--------
>>> from uppasd import Simulator
>>> from uppasd.magnons import compute_magnons
>>> with Simulator() as sim:
...     sim.init_simulation()
...     # Generate q-points (example: Gamma->X->M->Gamma path)
...     q_mesh = setup_q_mesh_path([...])
...     magnons = compute_magnons(sim, q_mesh)
...     dispersions = magnons['eigenvalues']

Notes
-----
- Requires Hamiltonian and magnetic moments to be initialized (call
  setup_all() and equilibrate the system first)
- Q-vectors should be in reciprocal lattice coordinates
- Automatic q-point generation (via spglib/seekpath) is a future enhancement
"""

import logging
from typing import Dict, Optional, Tuple, List
import numpy as np

logger = logging.getLogger(__name__)


def setup_q_mesh_path(
    path_points: List[Tuple[float, float, float]],
    points_per_segment: int = 50,
) -> np.ndarray:
    """
    Create a q-point mesh along a high-symmetry path in reciprocal space.
    
    Parameters
    ----------
    path_points : list of tuple
        List of (qx, qy, qz) points defining the path segments.
        Example: [(0, 0, 0), (0.5, 0, 0), (0.5, 0.5, 0)]
        for Gamma -> X -> M path.
    points_per_segment : int
        Number of points to generate per segment. Default: 50.
    
    Returns
    -------
    q_mesh : ndarray of shape (nq, 3)
        Q-points in reciprocal lattice coordinates.
    
    Examples
    --------
    >>> # Gamma -> X -> M -> Gamma path in cubic system
    >>> q_path = setup_q_mesh_path(
    ...     [(0, 0, 0), (0.5, 0, 0), (0.5, 0.5, 0), (0, 0, 0)],
    ...     points_per_segment=40
    ... )
    >>> print(f"Generated {len(q_path)} q-points")
    """
    if len(path_points) < 2:
        raise ValueError("Need at least 2 path points")
    
    path_points = np.array(path_points, dtype=float)
    q_mesh = []
    
    for i in range(len(path_points) - 1):
        start = path_points[i]
        end = path_points[i + 1]
        
        # Linear interpolation, include start but exclude end (except last segment)
        segment = np.linspace(start, end, points_per_segment, endpoint=(i == len(path_points) - 2))
        q_mesh.append(segment)
    
    q_mesh = np.vstack(q_mesh)
    logger.info(f"Generated q-mesh with {len(q_mesh)} points along {len(path_points)} path points")
    return q_mesh


def setup_q_mesh_grid(
    nq1: int = 10,
    nq2: int = 10,
    nq3: int = 1,
    reciprocal_lattice: Optional[np.ndarray] = None,
) -> np.ndarray:
    """
    Create a regular q-point mesh in reciprocal space.
    
    Parameters
    ----------
    nq1, nq2, nq3 : int
        Number of q-points along each reciprocal lattice vector.
        Default: (10, 10, 1) for 2D systems.
    reciprocal_lattice : ndarray of shape (3, 3), optional
        Reciprocal lattice vectors (rows). If provided, q-points are
        returned in Cartesian coordinates. Otherwise, returned in
        fractional coordinates.
    
    Returns
    -------
    q_mesh : ndarray of shape (nq, 3)
        Q-points in fractional coordinates (or Cartesian if
        reciprocal_lattice provided).
    
    Examples
    --------
    >>> q_mesh = setup_q_mesh_grid(nq1=20, nq2=20, nq3=1)
    >>> print(f"Generated {len(q_mesh)} q-points")
    """
    # Create fractional coordinates in [0, 1)
    q1 = np.linspace(0, 1, nq1, endpoint=False)
    q2 = np.linspace(0, 1, nq2, endpoint=False)
    q3 = np.linspace(0, 1, nq3, endpoint=False)
    
    q_mesh = np.array(np.meshgrid(q1, q2, q3, indexing='ij')).reshape(3, -1).T
    
    # Transform to Cartesian if reciprocal lattice provided
    if reciprocal_lattice is not None:
        reciprocal_lattice = np.array(reciprocal_lattice, dtype=float)
        if reciprocal_lattice.shape != (3, 3):
            raise ValueError("reciprocal_lattice must be shape (3, 3)")
        q_mesh = q_mesh @ reciprocal_lattice.T
    
    logger.info(f"Generated {len(q_mesh)}-point grid with nq1={nq1}, nq2={nq2}, nq3={nq3}")
    return q_mesh


def compute_magnons(
    simulator,
    q_mesh: np.ndarray,
    flag: int = 0,
) -> Dict:
    """
    Compute magnon dispersion and related properties using LSWT.
    
    This is the main user-facing function for magnon calculations. It requires:
    1. Simulator initialized and Hamiltonian mounted
    2. Magnetic moments equilibrated (from MD/MC relaxation or from file)
    3. Q-point mesh in reciprocal lattice coordinates
    
    Parameters
    ----------
    simulator : Simulator
        Initialized Simulator instance with Hamiltonian mounted.
    q_mesh : ndarray of shape (nq, 3)
        Q-points in reciprocal lattice coordinates.
    flag : int, optional
        Calculation type:
        - 0: Non-collinear AMS (default, returns nc_eval_q, nc_evec_q)
        - 1: Chern number calculation (returns nc_eval_qchern, nc_evec_qchern)
    
    Returns
    -------
    dict
        Contains:
        - 'eigenvalues': ndarray (nq, 2*NA) - magnon energies at each q
        - 'eigenvectors': ndarray (nq, 2*NA, 2*NA) - magnon modes (complex)
        - 'q_mesh': ndarray (nq, 3) - input q-points
        - 'nq': int - number of q-points
        - 'na': int - number of atoms in unit cell
        - 'nq_ext': int - extended number of q-points after Fortran processing
    
    Raises
    ------
    StateError
        If Simulator not initialized or Hamiltonian not mounted
    ValueError
        If q_mesh invalid or incompatible with system
    RuntimeError
        If Fortran calculation fails
    
    Examples
    --------
    >>> from uppasd import Simulator
    >>> from uppasd.magnons import compute_magnons, setup_q_mesh_path
    >>> 
    >>> with Simulator() as sim:
    ...     sim.init_simulation()
    ...     # Run some dynamics to equilibrate
    ...     sim.relax(mode='M', temperature=100, steps=500)
    ...     
    ...     # Define q-point path: Gamma -> X -> M
    ...     q_path = setup_q_mesh_path(
    ...         [(0, 0, 0), (0.5, 0, 0), (0.5, 0.5, 0)],
    ...         points_per_segment=30
    ...     )
    ...     
    ...     # Compute magnons
    ...     magnons = compute_magnons(sim, q_path)
    ...     print(f"Computed {magnons['nq']} q-points, 2*{magnons['na']} modes")
    ...     
    ...     # Access results
    ...     evals = magnons['eigenvalues']  # (nq_ext, 2*NA) from Fortran
    ...     evecs = magnons['eigenvectors']  # (nq_ext, 2*NA, 2*NA)
    
    Notes
    -----
    - Magnon energies are in units used by Fortran code (typically mRy)
    - Eigenvectors are complex; amplitude and phase carry physical info
    - For collinear systems, negative eigenvalues indicate instabilities
    - Fortran extends q-mesh internally (nq_ext = 6*nq typically) to handle
      -q points and phasons; returned eigenvalues have nq_ext points
    - Hamiltonian dimension hdim = 2*na (magnetic modes)
    """
    from uppasd.core import StateError
    import uppasd.pyasd as pyasd
    
    # Validate simulator state
    if not hasattr(simulator, 'state') or simulator.state not in ('initialized', 'running', 'completed'):
        raise StateError(
            "Simulator must be initialized and Hamiltonian mounted. "
            "Call sim.init_simulation() first."
        )
    
    # Validate q_mesh
    q_mesh = np.asarray(q_mesh, dtype=np.float64)
    if q_mesh.ndim != 2 or q_mesh.shape[1] != 3:
        raise ValueError(f"q_mesh must be shape (nq, 3), got {q_mesh.shape}")
    
    nq = len(q_mesh)
    natom = simulator.natom
    mensemble = simulator.mensemble
    # Use unit-cell atom count from InputData if available to avoid huge NA
    try:
        na = pyasd.get_na()
    except Exception:
        na = natom // mensemble  # fallback
    simid = getattr(simulator.inputdata, 'simid', 'sim')
    
    logger.info(f"Computing magnons for {nq} q-points, {na} atoms/cell, flag={flag}")
    
    try:
        # Get current magnetic moments from simulator
        emomm = np.asarray(simulator.moments, dtype=np.float64)
        # Ensure orientation (3, natom, mensemble)
        if emomm.shape[0] != 3 and 3 in emomm.shape:
            # If shape is (natom, mensemble, 3), transpose
            axes = emomm.shape
            if axes[-1] == 3:
                emomm = np.transpose(emomm, (2, 0, 1))
        if emomm.shape != (3, natom, mensemble):
            raise ValueError(f"moments shape {emomm.shape} != (3, {natom}, {mensemble})")

        # Moment magnitudes; safeguard against near-zero norms
        mmom = np.linalg.norm(emomm, axis=0)
        if np.any(mmom < 1e-3):
            logger.warning(
                "Moment magnitudes contain very small values (min=%.3e); "
                "clipping to 1e-3 for stability",
                float(mmom.min()),
            )
            mmom = np.clip(mmom, 1e-3, None)
        
        # Transpose q_mesh to Fortran convention: (3, nq)
        q_vect = q_mesh.T.copy()
        
        # Call Fortran wrapper via pyasd
        evals, evecs, nq_ext = pyasd.setup_tensor_hamiltonian(
            na, natom, mensemble, simid,
            emomm, mmom, q_vect, nq, flag
        )
        
        logger.info(f"Magnon eigenvalues shape: {evals.shape}, eigenvectors shape: {evecs.shape}")
        
        # Package results
        result = {
            'eigenvalues': evals,
            'eigenvectors': evecs,
            'q_mesh': q_mesh,
            'nq': nq,
            'na': na,
            'nq_ext': nq_ext,
            'flag': flag,
            'simid': simid,
        }
        
        return result
        
    except Exception as e:
        logger.error(f"Magnon calculation failed: {e}")
        raise RuntimeError(f"Failed to compute magnons: {e}") from e


def get_magnon_dispersion(magnons: Dict) -> Tuple[np.ndarray, np.ndarray]:
    """
    Extract magnon dispersion relation from magnon calculation results.
    
    Computes distance along q-path for plotting (useful for band structure plots).
    Note: Returns eigenvalues at nq_ext points (Fortran extended mesh), not just
    the input nq points. For visualization, may want to decimate or interpolate.
    
    Parameters
    ----------
    magnons : dict
        Results dictionary from compute_magnons()
    
    Returns
    -------
    q_distances : ndarray (nq_ext,)
        Distance along q-path (for x-axis in dispersion plots)
    energies : ndarray (nq_ext, 2*NA)
        Magnon energy bands
    
    Examples
    --------
    >>> magnons = compute_magnons(sim, q_mesh)
    >>> q_dist, energies = get_magnon_dispersion(magnons)
    >>> # Plot dispersion
    >>> for i in range(energies.shape[1]):
    ...     plt.plot(q_dist, energies[:, i], 'b-', lw=0.5)
    
    Notes
    -----
    - nq_ext is typically 6*nq_input (Fortran extends for -q and phasons)
    - q_distances are cumulative path distances, not meaningful in absolute sense
      when q-mesh is extended; use mainly for relative positioning
    - For original input q-points only, use energies[:nq_original, :] with
      appropriate indexing (depends on how Fortran orders the extended mesh)
    """
    evals = magnons['eigenvalues']
    nq_ext = evals.shape[0]
    
    # Create a placeholder q-distance array
    # In reality, would need mapping from extended to original q-points
    q_distances = np.linspace(0, 1, nq_ext)
    
    logger.info(f"Dispersion: {evals.shape[1]} bands, {nq_ext} q-points (extended)")
    
    return q_distances, evals


def filter_imaginary_modes(magnons: Dict, threshold: float = 1e-6) -> Dict:
    """
    Identify and filter out imaginary (unstable) magnon modes.
    
    Parameters
    ----------
    magnons : dict
        Results from compute_magnons()
    threshold : float
        Energy threshold below which modes are considered imaginary.
        Default: 1e-6 (in Fortran units, typically mRy)
    
    Returns
    -------
    dict
        Results with additional keys:
        - 'n_imaginary': int - number of imaginary modes detected
        - 'imaginary_mask': ndarray (nq, 2*NA) - True where E < threshold
    
    Examples
    --------
    >>> magnons = compute_magnons(sim, q_mesh)
    >>> filtered = filter_imaginary_modes(magnons, threshold=1e-5)
    >>> if filtered['n_imaginary'] > 0:
    ...     print(f"Warning: {filtered['n_imaginary']} imaginary modes")
    """
    evals = magnons['eigenvalues']
    mask = evals < threshold
    n_imag = np.sum(mask)
    
    result = magnons.copy()
    result['imaginary_mask'] = mask
    result['n_imaginary'] = n_imag
    
    if n_imag > 0:
        logger.warning(f"Found {n_imag} imaginary modes (E < {threshold})")
    
    return result
