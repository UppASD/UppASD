"""Magnon analysis moved into `uppasd.analysis`.

This file is a direct relocation of the previous `uppasd.magnons` module
into `uppasd.analysis.magnons` for clearer separation of physics code and
visualization helpers.
"""
import logging
from typing import Dict, Optional, Tuple, List
import numpy as np

logger = logging.getLogger(__name__)


def setup_q_mesh_path(
    path_points: List[Tuple[float, float, float]],
    points_per_segment: int = 50,
) -> np.ndarray:
    if len(path_points) < 2:
        raise ValueError("Need at least 2 path points")
    path_points = np.array(path_points, dtype=float)
    q_mesh = []
    for i in range(len(path_points) - 1):
        start = path_points[i]
        end = path_points[i + 1]
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
    q1 = np.linspace(0, 1, nq1, endpoint=False)
    q2 = np.linspace(0, 1, nq2, endpoint=False)
    q3 = np.linspace(0, 1, nq3, endpoint=False)
    q_mesh = np.array(np.meshgrid(q1, q2, q3, indexing='ij')).reshape(3, -1).T
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
    from uppasd.core import StateError
    import uppasd.pyasd as pyasd

    if not hasattr(simulator, 'state') or simulator.state not in ('initialized', 'running', 'completed'):
        raise StateError(
            "Simulator must be initialized and Hamiltonian mounted. Call sim.init_simulation() first."
        )

    q_mesh = np.asarray(q_mesh, dtype=np.float64)
    if q_mesh.ndim != 2 or q_mesh.shape[1] != 3:
        raise ValueError(f"q_mesh must be shape (nq, 3), got {q_mesh.shape}")

    nq = len(q_mesh)
    natom = simulator.natom
    mensemble = simulator.mensemble
    try:
        na = pyasd.get_na()
    except Exception:
        na = natom // mensemble
    simid = getattr(simulator.inputdata, 'simid', 'sim')

    logger.info(f"Computing magnons for {nq} q-points, {na} atoms/cell, flag={flag}")

    try:
        emomm = np.asarray(simulator.moments, dtype=np.float64)
        if emomm.shape[0] != 3 and 3 in emomm.shape:
            axes = emomm.shape
            if axes[-1] == 3:
                emomm = np.transpose(emomm, (2, 0, 1))
        if emomm.shape != (3, natom, mensemble):
            raise ValueError(f"moments shape {emomm.shape} != (3, {natom}, {mensemble})")

        mmom = np.linalg.norm(emomm, axis=0)
        if np.any(mmom < 1e-3):
            logger.warning("Moment magnitudes contain very small values; clipping to 1e-3 for stability")
            mmom = np.clip(mmom, 1e-3, None)

        q_vect = q_mesh.T.copy()

        evals, evecs, nq_ext = pyasd.setup_tensor_hamiltonian(
            na, natom, mensemble, simid,
            emomm, mmom, q_vect, nq, flag
        )

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
    evals = magnons['eigenvalues']
    nq_ext = evals.shape[0]
    q_distances = np.linspace(0, 1, nq_ext)
    logger.info(f"Dispersion: {evals.shape[1]} bands, {nq_ext} q-points (extended)")
    return q_distances, evals


def filter_imaginary_modes(magnons: Dict, threshold: float = 1e-6) -> Dict:
    evals = magnons['eigenvalues']
    mask = evals < threshold
    n_imag = np.sum(mask)
    result = magnons.copy()
    result['imaginary_mask'] = mask
    result['n_imaginary'] = n_imag
    if n_imag > 0:
        logger.warning(f"Found {n_imag} imaginary modes (E < {threshold})")
    return result
