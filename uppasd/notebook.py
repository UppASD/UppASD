"""
Notebook-oriented helper utilities for UppASD.

This module provides lightweight helpers for notebook workflows while keeping
dependencies minimal (numpy, matplotlib, logging only). Functions are organized
into categories:

**Structure Generation:**
- generate_supercell_coordinates: Build supercells from lattice and basis

**Magnetic Configuration Generators:**
- create_ferromagnetic_moments: All spins aligned
- create_antiferromagnetic_moments: Checkerboard AFM pattern
- create_vortex_moments: In-plane vortex texture
- create_skyrmion_moments: Néel skyrmion with radial profile
- perturb_moments: Add random perturbations with renormalization

**File I/O - Writing UppASD Inputs:**
- write_inpsd_file: Generate inpsd.dat from configuration dict
- write_posfile: Write posfile.dat with atomic positions
- write_momfile: Write momfile.dat with species-dependent moments
- write_jfile: Write jfile.dat with exchange interactions

**File I/O - Reading UppASD Files:**
- read_inpsd: Parse inpsd.dat metadata
- read_posfile: Parse atomic positions from posfile.dat
- read_momfile: Parse moments from momfile.dat
- read_jfile: Parse exchange interactions from jfile.dat
- load_restart_file: Load restart.*.out configurations
- read_generic_output_file: Parse text outputs with headers

**Relaxation Protocol Templates:**
- create_relaxation_protocol: Generate common relaxation configurations
- relaxation_protocol_mc_sd: MC annealing + SD (most common)
- relaxation_protocol_thermal_anneal: Gradual temperature reduction
- relaxation_protocol_quench: Rapid quench from high temperature

**Simulation Execution:**
- run_simulation_via_api: Run UppASD via Python API (no subprocess)
- run_simulation_api: Logging wrapper for simulations
- load_outputs: Load averages/totenergy outputs
- load_simulation_results: Load results by simid
- find_output_file: Locate output files by pattern

**Unit Conversions:**
- mev_to_mry, mry_to_mev: Energy unit conversions
- kelvin_to_mev, mev_to_kelvin: Temperature-energy conversions
- normalize_lattice_vectors: Normalize lattice for UppASD format
- reciprocal_lattice: Compute reciprocal lattice vectors
- compute_cell_volume: Calculate unit cell volume

**Simple Plotting (Matplotlib):**
- plot_magnetization_simple: M(t) with optional components
- plot_energy_simple: E(t) plot
- plot_moment_distribution: Histogram of moment magnitudes

**Magnon Calculations:**
- setup_magnon_q_mesh_path: Q-point paths for dispersion
- setup_magnon_q_mesh_grid: Regular Q-point meshes
- compute_magnon_dispersion: Run LSWT calculation

**Utilities:**
- iteration_to_time: Convert iteration indices to physical time
- print_system_info: Pretty-print parsed input parameters

Note: Advanced visualization (PyVista) and neighbor finding (ASE) are kept
in notebooks as examples to avoid heavy dependencies in the core package.
"""
from __future__ import annotations

import glob
import logging
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np

from uppasd.simulator import Simulator

# Type alias for q-point tuples
Qpoint = Tuple[float, float, float]

logger = logging.getLogger(__name__)


def read_inpsd(filepath: str = "inpsd.dat") -> Dict:
    """Parse UppASD input file (inpsd.dat) into a lightweight dict.

    Parameters
    ----------
    filepath : str
        Path to the inpsd.dat file

    Returns
    -------
    dict
        Parsed configuration with basic fields.
    """
    config = {
        "simid": None,
        "cell": np.eye(3, dtype=float),
        "ncell": [1, 1, 1],
        "posfile": None,
        "momfile": None,
        "exchange": None,
        "temperature": None,
        "mode": None,
        "nstep": None,
        "timestep": None,
        "delta_t": None,
        "sc_step": None,
        "sc_nstep": None,
        "mplambda1": None,
        "SDEalgh": None,
        "ipmode": None,
        "initmag": None,
        "mseed": None,
    }

    path = Path(filepath)
    if not path.exists():
        logger.warning("Input file not found: %s", filepath)
        return config

    with path.open("r", encoding="utf-8") as f:
        lines = f.readlines()

    cell_row = 0
    for line in lines:
        stripped = line.strip()
        if not stripped or stripped[0] in {"#", "%", "*", "="}:
            continue
        parts = stripped.split()
        key = parts[0].lower()

        if key == "simid" and len(parts) > 1:
            config["simid"] = parts[1]
        elif key == "cell" and len(parts) >= 4 and cell_row < 3:
            config["cell"][cell_row] = [float(x) for x in parts[1:4]]
            cell_row += 1
        elif key == "ncell" and len(parts) > 1:
            if len(parts) >= 4:
                config["ncell"] = [int(parts[1]), int(parts[2]), int(parts[3])]
            else:
                n = int(parts[1])
                config["ncell"] = [n, n, n]
        elif key == "posfile" and len(parts) > 1:
            config["posfile"] = parts[1]
        elif key == "momfile" and len(parts) > 1:
            config["momfile"] = parts[1]
        elif key == "exchange" and len(parts) > 1:
            config["exchange"] = parts[1]
        elif key == "temp" and len(parts) > 1:
            config["temperature"] = float(parts[1])
        elif key == "mode" and len(parts) > 1:
            config["mode"] = parts[1][0].upper()
        elif key == "nstep" and len(parts) > 1:
            config["nstep"] = int(parts[1])
        elif key == "timestep" and len(parts) > 1:
            config["timestep"] = float(parts[1])
            config["delta_t"] = float(parts[1])
        elif key == "mplambda1" and len(parts) > 1:
            config["mplambda1"] = float(parts[1])
        elif key == "sc_nstep" and len(parts) > 1:
            config["sc_nstep"] = int(parts[1])
        elif key == "sc_step" and len(parts) > 1:
            config["sc_step"] = int(parts[1])
        elif key == "sdealgh" and len(parts) > 1:
            config["SDEalgh"] = int(parts[1])
        elif key == "ipmode" and len(parts) > 1:
            config["ipmode"] = parts[1][0].upper()
        elif key == "initmag" and len(parts) > 1:
            config["initmag"] = int(parts[1])
        elif key == "mseed" and len(parts) > 1:
            config["mseed"] = int(parts[1])

    return config


def read_generic_output_file(filepath: str) -> Tuple[np.ndarray, List[str]]:
    """Read a generic UppASD text output (e.g., averages, totenergy)."""
    if not Path(filepath).exists():
        logger.warning("File not found: %s", filepath)
        return np.array([]), []

    with open(filepath, "r", encoding="utf-8") as f:
        lines = f.readlines()

    labels: List[str] = []
    data_lines: List[str] = []
    for line in lines:
        line = line.strip()
        if not line:
            continue
        if line.startswith("#"):
            header = line.lstrip("#").strip()
            labels = header.split()
        elif not line.startswith("%"):
            data_lines.append(line)

    data = np.array([list(map(float, ln.split())) for ln in data_lines]) if data_lines else np.array([])
    logger.info("Read %s: %s steps, %s columns", filepath, data.shape[0] if data.ndim else 0, data.shape[1] if data.ndim > 1 else 1)
    return data, labels


def run_simulation_via_api(
    mode: Optional[str] = None,
    ip_mode: Optional[str] = None,
    steps: Optional[int] = None,
    temperature: Optional[float] = None,
    timestep: Optional[float] = None,
    damping: Optional[float] = None,
    record_trajectory: bool = False,
) -> Dict:
    """Run UppASD using the Python API (no external process).

    Parameters
    ----------
    mode : str, optional
        Relaxation mode ('M', 'S', or 'H'). If None, uses input ipmode or 'M'.
    steps : int, optional
        Number of steps. If None, tries to use input nstep or defaults to 100.
    temperature : float, optional
        Temperature in Kelvin. If None, uses current backend temperature.
    timestep : float, optional
        Integration timestep (LLG). If None, uses current backend timestep.
    damping : float, optional
        LLG damping parameter. Default 0.5 if not provided.
    record_trajectory : bool, optional
        If True, record trajectory during relaxation.

    Returns
    -------
    dict
        Contains moments, energy, natom, mensemble, and optional trajectory.
    """
    results: Dict[str, object] = {}

    with Simulator(record_trajectory=record_trajectory) as sim:
        # Decide defaults from inputdata when available
        inferred_mode = mode or (getattr(sim.inputdata, "mode", None) or "M")
        inferred_ip_mode = ip_mode if ip_mode is not None else getattr(sim.inputdata, "ipmode", None)
        inferred_steps = steps or (sim.inputdata.nstep if sim.inputdata.nstep else 100)
        inferred_temp = temperature if temperature is not None else sim.temperature
        inferred_timestep = timestep if timestep is not None else sim.timestep
        inferred_damping = damping if damping is not None else 0.5

        moments = sim.relax(
            mode=inferred_mode,
            ip_mode=inferred_ip_mode,
            temperature=inferred_temp,
            steps=inferred_steps,
            timestep=inferred_timestep,
            damping=inferred_damping,
        )

        results["moments"] = moments
        results["energy"] = sim.energy
        results["natom"] = sim.natom
        results["mensemble"] = sim.mensemble
        if record_trajectory:
            results["trajectory"] = sim.get_trajectory()

    return results


def _find_output_file(pattern: str, simid: str = "*") -> Optional[str]:
    """Find first file matching pattern and simid."""
    full_pattern = pattern.replace("*", simid) if simid != "*" else pattern
    matches = glob.glob(full_pattern)
    if matches:
        return matches[0]
    return None


def find_output_file(pattern: str, simid: str = "*") -> Optional[str]:
    """Public wrapper to find first file matching pattern and simid."""
    return _find_output_file(pattern, simid)


def load_outputs(simid: Optional[str] = None) -> Dict:
    """Load averages and totenergy outputs for a simulation id.

    Parameters
    ----------
    simid : str, optional
        Simulation identifier. If None, matches any.

    Returns
    -------
    dict
        Parsed output data and time axes.
    """
    simid = simid or "*"
    results = {
        "averages_file": None,
        "averages_data": None,
        "averages_labels": None,
        "totenergy_file": None,
        "totenergy_data": None,
        "totenergy_labels": None,
        "avg_time": None,
        "ene_time": None,
        "magnetization": None,
        "energy": None,
    }

    ave_file = _find_output_file("averages.*.out", simid)
    if ave_file:
        data, labels = read_generic_output_file(ave_file)
        results["averages_file"] = ave_file
        results["averages_data"] = data
        results["averages_labels"] = labels
        if data.size > 0 and data.shape[1] >= 5:
            results["avg_time"] = data[:, 0]
            results["magnetization"] = data[:, 4]

    ene_file = _find_output_file("totenergy.*.out", simid)
    if ene_file:
        data, labels = read_generic_output_file(ene_file)
        results["totenergy_file"] = ene_file
        results["totenergy_data"] = data
        results["totenergy_labels"] = labels
        if data.size > 0 and data.shape[1] >= 2:
            results["ene_time"] = data[:, 0]
            results["energy"] = data[:, 1]

    return results


def iteration_to_time(iterations: np.ndarray, timestep: float) -> np.ndarray:
    """Convert iteration indices to physical time given a timestep."""
    return np.asarray(iterations, dtype=float) * float(timestep)


def _fmt(value, fmt: str, default="N/A"):
    """Helper to safely format possibly None values."""
    if value is None:
        return default
    try:
        return fmt.format(value)
    except Exception:
        return str(value)


def print_system_info(config: Dict) -> None:
    """Pretty-print system configuration from read_inpsd output."""
    print("=" * 70)
    print(" UppASD SIMULATION CONFIGURATION ".center(70, "="))
    print("=" * 70)

    print("\nSIMULATION PARAMETERS:")
    print(f"  - Simulation ID:        {config.get('simid')}")
    mode = config.get('mode')
    mode_name = {
        'S': 'Spin Dynamics',
        'M': 'Monte Carlo',
        'H': 'Heat Bath',
    }.get(mode, 'Unknown')
    print(f"  - Simulation Mode:      {mode} ({mode_name})")
    print(f"  - Temperature:          {_fmt(config.get('temperature'), '{:.2f}')} K")
    print(f"  - Measurement Steps:    {config.get('nstep')}")
    print(f"  - Time Step:            {_fmt(config.get('timestep'), '{:.2e}')} s")
    print(f"  - Damping Parameter:    {_fmt(config.get('mplambda1'), '{:.4f}')}")

    print("\nSTRUCTURAL PARAMETERS:")
    ncell = config.get('ncell', [None, None, None])
    print(f"  - Supercell Size:       {ncell[0]} x {ncell[1]} x {ncell[2]}")
    try:
        n_atoms = int(np.prod(ncell))
    except Exception:
        n_atoms = 'N/A'
    print(f"  - Total Atoms:          {n_atoms}")
    print("  - Lattice Vectors:")
    cell = config.get('cell', np.eye(3))
    for i, vec in enumerate(cell):
        print(f"    a{i+1} = [{_fmt(vec[0], '{:10.6f}')}, {_fmt(vec[1], '{:10.6f}')}, {_fmt(vec[2], '{:10.6f}')}] A")

    print("\nINPUT/OUTPUT FILES:")
    print(f"  - Position File:        {config.get('posfile')}")
    print(f"  - Moment File:          {config.get('momfile')}")
    print(f"  - Exchange File:        {config.get('exchange')}")

    print("\nALGORITHM SETTINGS:")
    sde_alg = {1: 'Midpoint', 2: 'Heun', 3: 'Heun3', 4: 'Heun_proper', 5: 'Depondt'}
    print(f"  - SDE Solver:           {sde_alg.get(config.get('SDEalgh'), 'Unknown')}")
    ip_mode_val = config.get('ipmode')
    ip_mode_name = {
        'M': 'Monte Carlo pre-run',
        'S': 'Spin dynamics pre-run',
        'H': 'Heat bath pre-run',
        'N': 'No initial phase',
    }.get(ip_mode_val, 'Unknown')
    print(f"  - Initial Phase Mode:   {ip_mode_val} ({ip_mode_name})")
    init_mag = {1: 'Random', 2: 'Cone', 3: 'Special', 4: 'From File'}
    print(f"  - Initial Magnetization:{init_mag.get(config.get('initmag'), 'Unknown')}")
    print(f"  - Random Seed:          {config.get('mseed')}")

    print("\nMEASUREMENT SETTINGS:")
    print(f"  - Measurement Interval: {config.get('sc_step') if 'sc_step' in config else config.get('SC_step')}")
    print(f"  - Total Measurements:   {config.get('sc_nstep') if 'sc_nstep' in config else config.get('SC_nstep')}")

    print("=" * 70 + "\n")


def run_simulation_api(
    mode: Optional[str] = None,
    ip_mode: Optional[str] = None,
    steps: Optional[int] = None,
    temperature: Optional[float] = None,
    timestep: Optional[float] = None,
    damping: Optional[float] = None,
    record_trajectory: bool = False,
) -> Dict:
    """Wrapper around run_simulation_via_api with logging."""
    results = run_simulation_via_api(
        mode=mode,
        ip_mode=ip_mode,
        steps=steps,
        temperature=temperature,
        timestep=timestep,
        damping=damping,
        record_trajectory=record_trajectory,
    )
    logger.info(
        "Simulation complete (API). natom=%s, mensemble=%s, energy=%s",
        results.get("natom"),
        results.get("mensemble"),
        results.get("energy"),
    )
    return results


def load_simulation_results(simid: Optional[str] = None, config: Optional[Dict] = None) -> Dict:
    """Load simulation outputs using load_outputs helper."""
    target_simid = simid or (config.get("simid") if config else None) or "*"
    results = load_outputs(target_simid)
    logger.info("Loaded results for simid=%s", target_simid)
    return results


def setup_magnon_q_mesh_path(
    path_points: List[Tuple[float, float, float]],
    points_per_segment: int = 50,
) -> np.ndarray:
    """
    Convenience wrapper for magnons.setup_q_mesh_path.
    
    Create a q-point mesh along a high-symmetry path in reciprocal space
    for magnon calculations.
    
    Parameters
    ----------
    path_points : list of tuple
        List of (qx, qy, qz) points defining path segments.
        Example: [(0, 0, 0), (0.5, 0, 0), (0.5, 0.5, 0)] for Gamma->X->M.
    points_per_segment : int
        Number of points per segment. Default: 50.
    
    Returns
    -------
    q_mesh : ndarray (nq, 3)
        Q-points in reciprocal lattice coordinates.
    """
    from uppasd.magnons import setup_q_mesh_path
    return setup_q_mesh_path(path_points, points_per_segment)


def setup_magnon_q_mesh_grid(
    nq1: int = 10,
    nq2: int = 10,
    nq3: int = 1,
) -> np.ndarray:
    """
    Convenience wrapper for magnons.setup_q_mesh_grid.
    
    Create a regular q-point mesh in reciprocal space.
    
    Parameters
    ----------
    nq1, nq2, nq3 : int
        Number of q-points per direction. Default: (10, 10, 1).
    
    Returns
    -------
    q_mesh : ndarray (nq, 3)
        Q-points in fractional coordinates.
    """
    from uppasd.magnons import setup_q_mesh_grid
    return setup_q_mesh_grid(nq1, nq2, nq3)


def compute_magnon_dispersion(
    simulator,
    q_mesh: np.ndarray,
) -> Dict:
    """
    Convenience wrapper for magnons.compute_magnons in notebooks.
    
    Compute magnon dispersion relation using LSWT.
    
    Parameters
    ----------
    simulator : Simulator
        Initialized Simulator with Hamiltonian mounted.
    q_mesh : ndarray (nq, 3)
        Q-points in reciprocal lattice coordinates.
    
    Returns
    -------
    dict
        Magnon eigenvalues, eigenvectors, and metadata.
    """
    from uppasd.magnons import compute_magnons, get_magnon_dispersion
    magnons = compute_magnons(simulator, q_mesh)
    q_dist, energies = get_magnon_dispersion(magnons)
    magnons['q_distances'] = q_dist
    magnons['energies'] = energies
    logger.info("Magnon calculation complete: %d q-points, %d bands", magnons['nq'], energies.shape[1])
    return magnons


# =============================================================================
# STRUCTURE GENERATION
# =============================================================================

def generate_supercell_coordinates(
    lattice_vectors: np.ndarray,
    basis_fractional: np.ndarray,
    n1: int,
    n2: int,
    n3: int,
    scale: float = 1.0,
) -> Tuple[np.ndarray, Optional[np.ndarray]]:
    """
    Generate Cartesian coordinates for a supercell from lattice and basis.

    Parameters
    ----------
    lattice_vectors : ndarray (3, 3)
        Normalized lattice vectors (unit cell)
    basis_fractional : ndarray (nbasis, 3) or (nbasis, 4)
        Fractional coordinates [x, y, z] or [x, y, z, type]
    n1, n2, n3 : int
        Supercell dimensions
    scale : float
        Lattice constant in Angstrom (alat)

    Returns
    -------
    positions : ndarray (natom, 3)
        Cartesian coordinates in Angstrom
    atom_types : ndarray (natom,) or None
        Atom type indices if basis_fractional has 4 columns, else None
    """
    nbasis = len(basis_fractional)
    natom = nbasis * n1 * n2 * n3
    
    has_types = basis_fractional.shape[1] == 4
    positions = np.zeros((natom, 3))
    atom_types = np.zeros(natom, dtype=int) if has_types else None
    
    cell_scaled = lattice_vectors * scale
    idx = 0
    
    for i in range(n1):
        for j in range(n2):
            for k in range(n3):
                for b in range(nbasis):
                    frac = basis_fractional[b, :3] + np.array([i, j, k])
                    positions[idx] = frac @ cell_scaled
                    if has_types:
                        atom_types[idx] = int(basis_fractional[b, 3])
                    idx += 1
    
    logger.info("Generated supercell: %d atoms (%dx%dx%d)", natom, n1, n2, n3)
    return positions, atom_types


# =============================================================================
# MAGNETIC CONFIGURATION GENERATORS
# =============================================================================

def create_ferromagnetic_moments(
    natom: int,
    moment_magnitude: float = 1.0,
    direction: Tuple[float, float, float] = (0.0, 0.0, 1.0),
) -> np.ndarray:
    """
    Create ferromagnetic configuration with all spins aligned.

    Parameters
    ----------
    natom : int
        Number of atoms
    moment_magnitude : float
        Magnitude of magnetic moment (μB)
    direction : tuple (3,)
        Unit vector direction (will be normalized)

    Returns
    -------
    moments : ndarray (natom, 3)
        Magnetic moment vectors
    """
    direction_arr = np.array(direction, dtype=float)
    direction_arr /= np.linalg.norm(direction_arr)
    moments = np.tile(direction_arr * moment_magnitude, (natom, 1))
    logger.info("Created FM configuration: %d atoms, M=%.2f μB", natom, moment_magnitude)
    return moments


def create_antiferromagnetic_moments(
    positions: np.ndarray,
    moment_magnitude: float = 1.0,
    axis: str = 'z',
) -> np.ndarray:
    """
    Create antiferromagnetic checkerboard pattern based on position parity.

    Parameters
    ----------
    positions : ndarray (natom, 3)
        Atomic positions in Angstrom
    moment_magnitude : float
        Magnitude of magnetic moment (μB)
    axis : str
        'x', 'y', or 'z' - direction of moment alignment

    Returns
    -------
    moments : ndarray (natom, 3)
        Magnetic moment vectors in checkerboard pattern
    """
    natom = len(positions)
    moments = np.zeros((natom, 3))
    
    axis_idx = {'x': 0, 'y': 1, 'z': 2}[axis.lower()]
    
    for i in range(natom):
        # Checkerboard based on sum of integer coordinates
        coord_sum = int(positions[i, 0]) + int(positions[i, 1]) + int(positions[i, 2])
        sign = 1 if coord_sum % 2 == 0 else -1
        moments[i, axis_idx] = sign * moment_magnitude
    
    logger.info("Created AFM configuration: %d atoms, M=%.2f μB", natom, moment_magnitude)
    return moments


def create_vortex_moments(
    positions: np.ndarray,
    moment_magnitude: float = 1.0,
    center: Optional[np.ndarray] = None,
    axis: str = 'z',
) -> np.ndarray:
    """
    Create in-plane vortex magnetic texture.

    Parameters
    ----------
    positions : ndarray (natom, 3)
        Atomic positions in Angstrom
    moment_magnitude : float
        Magnitude of magnetic moment (μB)
    center : ndarray (3,), optional
        Vortex center. If None, uses geometric center
    axis : str
        'x', 'y', or 'z' - normal to vortex plane

    Returns
    -------
    moments : ndarray (natom, 3)
        Magnetic moment vectors forming vortex
    """
    natom = len(positions)
    moments = np.zeros((natom, 3))
    
    if center is None:
        center = positions.mean(axis=0)
    
    # Define plane based on axis
    axis_idx = {'x': 0, 'y': 1, 'z': 2}[axis.lower()]
    plane_axes = [i for i in range(3) if i != axis_idx]
    
    for i in range(natom):
        r = positions[i] - center
        r_plane = np.array([r[plane_axes[0]], r[plane_axes[1]]])
        r_norm = np.linalg.norm(r_plane)
        
        if r_norm > 1e-6:
            # Tangential direction in plane
            theta = np.arctan2(r_plane[1], r_plane[0])
            moments[i, plane_axes[0]] = -np.sin(theta) * moment_magnitude
            moments[i, plane_axes[1]] = np.cos(theta) * moment_magnitude
    
    logger.info("Created vortex configuration: %d atoms, M=%.2f μB", natom, moment_magnitude)
    return moments


def create_skyrmion_moments(
    positions: np.ndarray,
    moment_magnitude: float = 1.0,
    center: Optional[np.ndarray] = None,
    radius: float = 5.0,
    profile: str = 'tanh',
    axis: str = 'z',
) -> np.ndarray:
    """
    Create Néel skyrmion with radial profile.

    Parameters
    ----------
    positions : ndarray (natom, 3)
        Atomic positions in Angstrom
    moment_magnitude : float
        Magnitude of magnetic moment (μB)
    center : ndarray (3,), optional
        Skyrmion center. If None, uses geometric center
    radius : float
        Skyrmion radius in Angstrom
    profile : str
        'tanh': θ(r) = π(1 - tanh((r-R)/2))
        'linear': Linear interpolation from 0 to π
    axis : str
        'x', 'y', or 'z' - out-of-plane direction

    Returns
    -------
    moments : ndarray (natom, 3)
        Magnetic moment vectors forming skyrmion
    """
    natom = len(positions)
    moments = np.zeros((natom, 3))
    
    if center is None:
        center = positions.mean(axis=0)
    
    axis_idx = {'x': 0, 'y': 1, 'z': 2}[axis.lower()]
    plane_axes = [i for i in range(3) if i != axis_idx]
    
    for i in range(natom):
        r_vec = positions[i] - center
        r_plane = np.array([r_vec[plane_axes[0]], r_vec[plane_axes[1]]])
        r = np.linalg.norm(r_plane)
        
        # Radial profile for theta angle
        if profile == 'tanh':
            theta = np.pi * (1.0 - np.tanh((r - radius) / 2.0))
        elif profile == 'linear':
            theta = np.pi * min(r / radius, 1.0)
        else:
            theta = 0.0
        
        # Azimuthal angle
        phi = np.arctan2(r_plane[1], r_plane[0]) if r > 1e-6 else 0.0
        
        # Néel skyrmion: radial helicity
        m_r = np.sin(theta) * np.cos(phi)
        m_phi = np.sin(theta) * np.sin(phi)
        m_z = np.cos(theta)
        
        moments[i, plane_axes[0]] = m_r * moment_magnitude
        moments[i, plane_axes[1]] = m_phi * moment_magnitude
        moments[i, axis_idx] = m_z * moment_magnitude
    
    logger.info("Created skyrmion configuration: %d atoms, R=%.1f Å", natom, radius)
    return moments


def perturb_moments(
    moments: np.ndarray,
    perturbation_scale: float = 0.01,
    seed: Optional[int] = None,
    renormalize: bool = True,
    moment_magnitudes: Optional[np.ndarray] = None,
) -> np.ndarray:
    """
    Add random perturbation to moment configuration.

    Useful for Monte Carlo simulations needing slightly disordered initial states.

    Parameters
    ----------
    moments : ndarray (natom, 3)
        Input magnetic moments
    perturbation_scale : float
        Scale of random perturbation
    seed : int, optional
        Random seed for reproducibility
    renormalize : bool
        If True, renormalize to preserve moment magnitudes
    moment_magnitudes : ndarray (natom,), optional
        Target magnitudes per atom for renormalization

    Returns
    -------
    perturbed_moments : ndarray (natom, 3)
        Perturbed magnetic moments
    """
    if seed is not None:
        np.random.seed(seed)
    
    perturbed = moments.copy()
    perturbed += np.random.randn(*moments.shape) * perturbation_scale
    
    if renormalize:
        for i in range(len(perturbed)):
            target_mag = moment_magnitudes[i] if moment_magnitudes is not None else np.linalg.norm(moments[i])
            current_norm = np.linalg.norm(perturbed[i])
            if current_norm > 1e-10:
                perturbed[i] = perturbed[i] / current_norm * target_mag
    
    logger.info("Perturbed moments: scale=%.3f, renormalize=%s", perturbation_scale, renormalize)
    return perturbed


# =============================================================================
# FILE I/O - UppASD INPUT FILES
# =============================================================================

def write_posfile(
    filepath: str,
    basis_fractional: np.ndarray,
    atom_types: Optional[np.ndarray] = None,
) -> None:
    """
    Write posfile.dat in UppASD format.

    Parameters
    ----------
    filepath : str
        Output file path
    basis_fractional : ndarray (nbasis, 3) or (nbasis, 4)
        Fractional coordinates [x, y, z] or [x, y, z, type]
    atom_types : ndarray (nbasis,), optional
        Atom types if not in basis_fractional
    """
    has_types = basis_fractional.shape[1] == 4
    
    with open(filepath, 'w') as f:
        for i, pos in enumerate(basis_fractional[:, :3], 1):
            if has_types:
                atype = int(basis_fractional[i - 1, 3])
            elif atom_types is not None:
                atype = int(atom_types[i - 1])
            else:
                atype = 1
            f.write(f"{i}  {atype}  {pos[0]:.6f}  {pos[1]:.6f}  {pos[2]:.6f}\n")
    
    logger.info("Wrote posfile: %s (%d atoms)", filepath, len(basis_fractional))


def write_momfile(
    filepath: str,
    moment_magnitudes: dict | np.ndarray,
    initial_directions: Optional[np.ndarray] = None,
    atom_types: Optional[np.ndarray] = None,
    anisotropy: Optional[dict] = None,
) -> None:
    """
    Write momfile.dat with species-dependent moments.

    Parameters
    ----------
    filepath : str
        Output file path
    moment_magnitudes : dict or ndarray
        Dict: {type: magnitude} or array of magnitudes per atom
    initial_directions : ndarray (natom, 3), optional
        Initial moment directions (default: [0, 0, 1])
    atom_types : ndarray (natom,), optional
        Atom types (required if moment_magnitudes is dict)
    anisotropy : dict, optional
        {type: (K, ex, ey, ez)} for species-dependent anisotropy
    """
    if isinstance(moment_magnitudes, dict):
        if atom_types is None:
            raise ValueError("atom_types required when moment_magnitudes is dict")
        natom = len(atom_types)
        mags = np.array([moment_magnitudes[t] for t in atom_types])
    else:
        mags = np.asarray(moment_magnitudes)
        natom = len(mags)
    
    if initial_directions is None:
        initial_directions = np.tile([0.0, 0.0, 1.0], (natom, 1))
    
    with open(filepath, 'w') as f:
        for i in range(natom):
            atype = atom_types[i] if atom_types is not None else 1
            mag = mags[i]
            direction = initial_directions[i]
            
            if anisotropy and atype in anisotropy:
                K, ex, ey, ez = anisotropy[atype]
                f.write(f"{i + 1}  1  {mag:.6f}  {K:.6f}  {ex:.6f}  {ey:.6f}  {ez:.6f}  "
                        f"{direction[0]:.6f} {direction[1]:.6f} {direction[2]:.6f}\n")
            else:
                f.write(f"{i + 1}  1  {mag:.6f}  {direction[0]:.6f} {direction[1]:.6f} {direction[2]:.6f}\n")
    
    logger.info("Wrote momfile: %s (%d atoms)", filepath, natom)


def write_jfile(
    filepath: str,
    exchange_list: List[Tuple],
) -> None:
    """
    Write jfile.dat from exchange interaction list.

    Parameters
    ----------
    filepath : str
        Output file path
    exchange_list : list of tuples
        Each tuple: (type_i, type_j, rx, ry, rz, J_ij)
    """
    with open(filepath, 'w') as f:
        for entry in exchange_list:
            if len(entry) == 6:
                type_i, type_j, rx, ry, rz, J_ij = entry
                f.write(f"{type_i} {type_j}  {rx:.6f}  {ry:.6f}  {rz:.6f}  {J_ij:.6f}\n")
            else:
                logger.warning("Invalid exchange entry: %s", entry)
    
    logger.info("Wrote jfile: %s (%d interactions)", filepath, len(exchange_list))


def write_inpsd_file(
    filepath: str,
    config: Dict,
) -> None:
    """
    Write inpsd.dat from configuration dictionary.

    Parameters
    ----------
    filepath : str
        Output file path
    config : dict
        Configuration with keys: simid, ncell, cell, alat, posfile, momfile,
        exchange, mode, temp, nstep, timestep, damping, ip_mode, ip_mcanneal, etc.

    Example
    -------
    >>> config = {
    ...     'simid': 'test',
    ...     'ncell': [4, 4, 3],
    ...     'cell': np.eye(3),
    ...     'alat': 2.5,
    ...     'posfile': './posfile.dat',
    ...     'momfile': './momfile.dat',
    ...     'exchange': './jfile.dat',
    ...     'mode': 'S',
    ...     'temp': 100,
    ...     'nstep': 5000,
    ...     'timestep': 1e-16,
    ...     'damping': 0.05,
    ... }
    >>> write_inpsd_file('inpsd.dat', config)
    """
    lines = []
    
    # Simulation ID
    lines.append(f"simid  {config.get('simid', 'simulation')}")
    
    # Supercell
    ncell = config.get('ncell', [1, 1, 1])
    lines.append(f"ncell  {ncell[0]}  {ncell[1]}  {ncell[2]}")
    
    # Boundary conditions
    bc = config.get('BC', 'P  P  P')
    lines.append(f"BC  {bc}")
    
    # Cell
    lines.append("cell")
    cell = config.get('cell', np.eye(3))
    for row in cell:
        lines.append(f"{row[0]:.6f}  {row[1]:.6f}  {row[2]:.6f}")
    
    # Symmetry and lattice constant
    lines.append(f"Sym  {config.get('Sym', 0)}")
    lines.append(f"alat  {config.get('alat', 1.0):.6f}")
    lines.append("")
    
    # Input files
    lines.append(f"posfile  {config.get('posfile', './posfile.dat')}")
    lines.append(f"momfile  {config.get('momfile', './momfile.dat')}")
    if 'exchange' in config:
        lines.append(f"exchange  {config['exchange']}")
    lines.append("")
    
    # Structure output
    if config.get('do_prnstruct', 1):
        lines.append("do_prnstruct  1")
        lines.append("")
    
    # Ensemble
    lines.append(f"Mensemble  {config.get('Mensemble', 1)}")
    lines.append(f"Initmag  {config.get('Initmag', 3)}")
    lines.append("")
    
    # Initial phase
    if 'ip_mode' in config:
        lines.append(f"ip_mode  {config['ip_mode']}")
        if 'ip_mcanneal' in config and config['ip_mcanneal']:
            lines.append("ip_mcanneal  1")
            anneal_params = config.get('ip_mcanneal_params', [1000, 300])
            lines.append(f"{anneal_params[0]} {anneal_params[1]}")
        lines.append("")
    
    # Main simulation
    lines.append(f"mode  {config.get('mode', 'S')}")
    lines.append(f"temp  {config.get('temp', 100)}")
    if 'damping' in config:
        lines.append(f"damping  {config['damping']:.6f}")
    if 'timestep' in config:
        lines.append(f"timestep  {config['timestep']:.2e}")
    lines.append(f"nstep  {config.get('nstep', 1000)}")
    lines.append("")
    
    # Measurements
    if config.get('do_avrg', 'Y') == 'Y':
        lines.append("do_avrg  Y")
        lines.append(f"avrg_step  {config.get('avrg_step', 100)}")
    
    with open(filepath, 'w') as f:
        f.write('\n'.join(lines))
    
    logger.info("Wrote inpsd: %s (mode=%s, nstep=%d)", filepath, config.get('mode'), config.get('nstep'))


# =============================================================================
# FILE I/O - READING UppASD FILES
# =============================================================================

def read_posfile(filepath: str) -> Tuple[np.ndarray, np.ndarray]:
    """
    Parse posfile.dat.

    Returns
    -------
    positions : ndarray (natom, 3)
        Fractional coordinates
    atom_types : ndarray (natom,)
        Atom type indices
    """
    positions = []
    atom_types = []
    
    with open(filepath, 'r') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            parts = line.split()
            if len(parts) >= 5:
                atom_types.append(int(parts[1]))
                positions.append([float(parts[2]), float(parts[3]), float(parts[4])])
    
    positions = np.array(positions)
    atom_types = np.array(atom_types)
    logger.info("Read posfile: %s (%d atoms)", filepath, len(positions))
    return positions, atom_types


def read_momfile(filepath: str) -> Dict:
    """
    Parse momfile.dat.

    Returns
    -------
    moments : dict
        Keys: 'magnitudes', 'directions', 'types', 'anisotropy' (if present)
    """
    magnitudes = []
    directions = []
    anisotropies = []
    
    with open(filepath, 'r') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            parts = line.split()
            if len(parts) >= 6:
                mag = float(parts[2])
                dir_vec = [float(parts[3]), float(parts[4]), float(parts[5])]
                magnitudes.append(mag)
                directions.append(dir_vec)
                
                if len(parts) >= 10:
                    K = float(parts[3])
                    aniso = [float(parts[4]), float(parts[5]), float(parts[6])]
                    anisotropies.append((K, aniso))
    
    result = {
        'magnitudes': np.array(magnitudes),
        'directions': np.array(directions),
    }
    if anisotropies:
        result['anisotropy'] = anisotropies
    
    logger.info("Read momfile: %s (%d moments)", filepath, len(magnitudes))
    return result


def read_jfile(filepath: str) -> List[Tuple]:
    """
    Parse jfile.dat.

    Returns
    -------
    exchange_list : list of tuples
        Each tuple: (type_i, type_j, rx, ry, rz, J_ij)
    """
    exchange_list = []
    
    with open(filepath, 'r') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            parts = line.split()
            if len(parts) >= 6:
                type_i = int(parts[0])
                type_j = int(parts[1])
                rx, ry, rz = float(parts[2]), float(parts[3]), float(parts[4])
                J_ij = float(parts[5])
                exchange_list.append((type_i, type_j, rx, ry, rz, J_ij))
    
    logger.info("Read jfile: %s (%d interactions)", filepath, len(exchange_list))
    return exchange_list


def load_restart_file(
    filepath: str,
    natom: int,
    mensemble: int = 1,
) -> np.ndarray:
    """
    Load restart.*.out magnetic configuration.

    Parameters
    ----------
    filepath : str
        Path to restart file
    natom : int
        Number of atoms
    mensemble : int
        Number of ensembles

    Returns
    -------
    moments : ndarray (natom, 3, mensemble)
        Magnetic moment configuration
    """
    moments = np.zeros((natom, 3, mensemble))
    
    with open(filepath, 'r') as f:
        lines = f.readlines()
    
    # Skip header lines
    data_lines = [l for l in lines if not l.strip().startswith('#') and l.strip()]
    
    for i, line in enumerate(data_lines[:natom * mensemble]):
        parts = line.split()
        if len(parts) >= 5:
            ens = int(parts[0]) - 1
            atom = int(parts[1]) - 1
            mx, my, mz = float(parts[2]), float(parts[3]), float(parts[4])
            moments[atom, :, ens] = [mx, my, mz]
    
    logger.info("Loaded restart: %s (%d atoms, %d ens)", filepath, natom, mensemble)
    return moments


# =============================================================================
# RELAXATION PROTOCOL TEMPLATES
# =============================================================================

def relaxation_protocol_mc_sd(
    mc_steps: int = 2000,
    mc_temp: float = 100,
    sd_steps: int = 2000,
    sd_temp: float = 50,
    damping: float = 0.1,
    timestep: float = 1e-16,
) -> Dict:
    """
    Monte Carlo annealing followed by Spin Dynamics relaxation.

    Most common relaxation protocol for finding ground states.

    Parameters
    ----------
    mc_steps : int
        Monte Carlo annealing steps
    mc_temp : float
        MC temperature (K)
    sd_steps : int
        Spin dynamics steps
    sd_temp : float
        SD temperature (K)
    damping : float
        LLG damping parameter
    timestep : float
        Integration timestep (s)

    Returns
    -------
    config : dict
        Configuration ready for write_inpsd_file
    """
    return {
        'ip_mode': 'M',
        'ip_mcanneal': True,
        'ip_mcanneal_params': [mc_steps, mc_temp],
        'mode': 'S',
        'temp': sd_temp,
        'damping': damping,
        'timestep': timestep,
        'nstep': sd_steps,
        'do_avrg': 'Y',
        'avrg_step': max(100, sd_steps // 20),
    }


def relaxation_protocol_thermal_anneal(
    steps_per_temp: int = 1000,
    temp_start: float = 300,
    temp_end: float = 10,
    n_temps: int = 10,
    mode: str = 'M',
) -> List[Dict]:
    """
    Gradual temperature reduction protocol.

    Returns a list of configurations for sequential temperature steps.

    Parameters
    ----------
    steps_per_temp : int
        Steps at each temperature
    temp_start : float
        Starting temperature (K)
    temp_end : float
        Final temperature (K)
    n_temps : int
        Number of temperature steps
    mode : str
        'M' (Monte Carlo) or 'S' (Spin Dynamics)

    Returns
    -------
    configs : list of dict
        List of configurations for each temperature
    """
    temps = np.linspace(temp_start, temp_end, n_temps)
    configs = []
    
    for i, temp in enumerate(temps):
        config = {
            'mode': mode,
            'temp': temp,
            'nstep': steps_per_temp,
            'do_avrg': 'Y',
            'avrg_step': max(10, steps_per_temp // 10),
        }
        if mode == 'S':
            config['damping'] = 0.5
            config['timestep'] = 1e-16
        configs.append(config)
    
    logger.info("Created thermal anneal protocol: %d temps, %.1f->%.1f K", n_temps, temp_start, temp_end)
    return configs


def relaxation_protocol_quench(
    mc_steps: int = 1000,
    temp_high: float = 500,
    sd_steps: int = 5000,
    temp_low: float = 10,
    damping: float = 0.5,
    timestep: float = 1e-16,
) -> Dict:
    """
    Rapid quench from high temperature.

    Parameters
    ----------
    mc_steps : int
        Monte Carlo equilibration steps at high T
    temp_high : float
        High temperature (K)
    sd_steps : int
        Spin dynamics quench steps
    temp_low : float
        Low temperature (K)
    damping : float
        LLG damping (higher = faster quench)
    timestep : float
        Integration timestep (s)

    Returns
    -------
    config : dict
        Configuration for quench protocol
    """
    return {
        'ip_mode': 'M',
        'ip_mcanneal': True,
        'ip_mcanneal_params': [mc_steps, temp_high],
        'mode': 'S',
        'temp': temp_low,
        'damping': damping,
        'timestep': timestep,
        'nstep': sd_steps,
        'do_avrg': 'Y',
        'avrg_step': max(100, sd_steps // 20),
    }


def create_relaxation_protocol(protocol_type: str = 'mc_sd', **kwargs) -> Dict | List[Dict]:
    """
    Generate inpsd configuration dict for common relaxation protocols.

    Parameters
    ----------
    protocol_type : str
        'mc_sd': Monte Carlo annealing + Spin Dynamics
        'mc_only': Pure Monte Carlo
        'sd_only': Pure Spin Dynamics with damping
        'thermal_anneal': Temperature ramp (returns list)
        'quench': Fast quench from high T

    **kwargs : additional parameters passed to specific protocol functions

    Returns
    -------
    config : dict or list of dict
        Configuration(s) ready for write_inpsd_file

    Examples
    --------
    >>> config = create_relaxation_protocol('mc_sd', mc_steps=2000, sd_temp=50)
    >>> configs = create_relaxation_protocol('thermal_anneal', n_temps=10)
    """
    if protocol_type == 'mc_sd':
        return relaxation_protocol_mc_sd(**kwargs)
    elif protocol_type == 'thermal_anneal':
        return relaxation_protocol_thermal_anneal(**kwargs)
    elif protocol_type == 'quench':
        return relaxation_protocol_quench(**kwargs)
    elif protocol_type == 'mc_only':
        steps = kwargs.get('steps', 10000)
        temp = kwargs.get('temp', 100)
        return {
            'mode': 'M',
            'temp': temp,
            'nstep': steps,
            'do_avrg': 'Y',
            'avrg_step': max(100, steps // 20),
        }
    elif protocol_type == 'sd_only':
        steps = kwargs.get('steps', 5000)
        temp = kwargs.get('temp', 50)
        damping = kwargs.get('damping', 0.1)
        timestep = kwargs.get('timestep', 1e-16)
        return {
            'mode': 'S',
            'temp': temp,
            'damping': damping,
            'timestep': timestep,
            'nstep': steps,
            'do_avrg': 'Y',
            'avrg_step': max(100, steps // 20),
        }
    else:
        raise ValueError(f"Unknown protocol type: {protocol_type}")


# =============================================================================
# UNIT CONVERSIONS & UTILITIES
# =============================================================================

def mev_to_mry(energy_mev: float) -> float:
    """Convert meV to mRy (1 meV ≈ 0.0735 mRy)."""
    return energy_mev * 0.07349862


def mry_to_mev(energy_mry: float) -> float:
    """Convert mRy to meV (1 mRy ≈ 13.6 meV)."""
    return energy_mry / 0.07349862


def kelvin_to_mev(temp_kelvin: float) -> float:
    """Convert Kelvin to meV (k_B = 0.08617 meV/K)."""
    return temp_kelvin * 0.08617333262


def mev_to_kelvin(energy_mev: float) -> float:
    """Convert meV to Kelvin."""
    return energy_mev / 0.08617333262


def normalize_lattice_vectors(lattice_vectors: np.ndarray) -> np.ndarray:
    """
    Ensure lattice vectors form a proper unit cell.

    UppASD expects normalized vectors with alat separate.
    This function normalizes by the determinant^(1/3).

    Parameters
    ----------
    lattice_vectors : ndarray (3, 3)
        Input lattice vectors

    Returns
    -------
    normalized : ndarray (3, 3)
        Normalized lattice vectors
    """
    det = np.linalg.det(lattice_vectors)
    if abs(det) < 1e-10:
        logger.warning("Lattice vectors are nearly singular!")
        return lattice_vectors
    
    scale = abs(det) ** (1.0 / 3.0)
    normalized = lattice_vectors / scale
    logger.info("Normalized lattice: det=%.6f, scale=%.6f", det, scale)
    return normalized


def reciprocal_lattice(lattice_vectors: np.ndarray) -> np.ndarray:
    """
    Compute reciprocal lattice vectors.

    Parameters
    ----------
    lattice_vectors : ndarray (3, 3)
        Real space lattice vectors (rows)

    Returns
    -------
    recip_vectors : ndarray (3, 3)
        Reciprocal lattice vectors: b_i = 2π * (a_j × a_k) / V
    """
    a1, a2, a3 = lattice_vectors[0], lattice_vectors[1], lattice_vectors[2]
    volume = np.dot(a1, np.cross(a2, a3))
    
    if abs(volume) < 1e-10:
        raise ValueError("Lattice vectors are singular!")
    
    b1 = 2 * np.pi * np.cross(a2, a3) / volume
    b2 = 2 * np.pi * np.cross(a3, a1) / volume
    b3 = 2 * np.pi * np.cross(a1, a2) / volume
    
    return np.array([b1, b2, b3])


def compute_cell_volume(lattice_vectors: np.ndarray, alat: float = 1.0) -> float:
    """
    Compute unit cell volume.

    Parameters
    ----------
    lattice_vectors : ndarray (3, 3)
        Lattice vectors
    alat : float
        Lattice constant (default 1.0)

    Returns
    -------
    volume : float
        Cell volume in Angstrom^3
    """
    scaled = lattice_vectors * alat
    a1, a2, a3 = scaled[0], scaled[1], scaled[2]
    volume = abs(np.dot(a1, np.cross(a2, a3)))
    return volume


# =============================================================================
# SIMPLE PLOTTING (MATPLOTLIB ONLY)
# =============================================================================

def plot_magnetization_simple(
    time: np.ndarray,
    magnetization: np.ndarray,
    components: bool = True,
    figsize: Tuple[int, int] = (10, 5),
):
    """
    Simple magnetization vs time plot using matplotlib.

    Parameters
    ----------
    time : ndarray
        Time axis
    magnetization : ndarray (nsteps,) or (nsteps, 3)
        Magnetization magnitude or [Mx, My, Mz]
    components : bool
        If True and magnetization has 3 columns, plot components separately
    figsize : tuple
        Figure size

    Returns
    -------
    fig, ax : matplotlib figure and axes
    """
    try:
        import matplotlib.pyplot as plt
    except ImportError:
        logger.error("Matplotlib not available")
        return None, None
    
    fig, ax = plt.subplots(figsize=figsize)
    
    if magnetization.ndim == 1 or not components:
        mag = np.linalg.norm(magnetization, axis=1) if magnetization.ndim > 1 else magnetization
        ax.plot(time, mag, 'b-', linewidth=1.5, label='|M|')
    else:
        ax.plot(time, magnetization[:, 0], 'r-', label='Mx', alpha=0.8)
        ax.plot(time, magnetization[:, 1], 'g-', label='My', alpha=0.8)
        ax.plot(time, magnetization[:, 2], 'b-', label='Mz', alpha=0.8)
    
    ax.set_xlabel('Time')
    ax.set_ylabel('Magnetization')
    ax.set_title('Magnetization vs Time')
    ax.legend()
    ax.grid(True, alpha=0.3)
    plt.tight_layout()
    
    return fig, ax


def plot_energy_simple(
    time: np.ndarray,
    energy: np.ndarray,
    figsize: Tuple[int, int] = (10, 5),
):
    """
    Simple energy vs time plot.

    Parameters
    ----------
    time : ndarray
        Time axis
    energy : ndarray
        Energy values
    figsize : tuple
        Figure size

    Returns
    -------
    fig, ax : matplotlib figure and axes
    """
    try:
        import matplotlib.pyplot as plt
    except ImportError:
        logger.error("Matplotlib not available")
        return None, None
    
    fig, ax = plt.subplots(figsize=figsize)
    ax.plot(time, energy, 'r-', linewidth=1.5)
    ax.set_xlabel('Time')
    ax.set_ylabel('Energy')
    ax.set_title('Total Energy vs Time')
    ax.grid(True, alpha=0.3)
    plt.tight_layout()
    
    return fig, ax


def plot_moment_distribution(
    moments: np.ndarray,
    bins: int = 50,
    figsize: Tuple[int, int] = (8, 6),
):
    """
    Histogram of moment magnitudes.

    Parameters
    ----------
    moments : ndarray (natom, 3)
        Magnetic moments
    bins : int
        Number of histogram bins
    figsize : tuple
        Figure size

    Returns
    -------
    fig, ax : matplotlib figure and axes
    """
    try:
        import matplotlib.pyplot as plt
    except ImportError:
        logger.error("Matplotlib not available")
        return None, None
    
    magnitudes = np.linalg.norm(moments, axis=1)
    
    fig, ax = plt.subplots(figsize=figsize)
    ax.hist(magnitudes, bins=bins, color='steelblue', alpha=0.7, edgecolor='black')
    ax.axvline(magnitudes.mean(), color='red', linestyle='--', 
               label=f'Mean: {magnitudes.mean():.3f}')
    ax.set_xlabel('Moment Magnitude (μB)')
    ax.set_ylabel('Count')
    ax.set_title('Distribution of Magnetic Moment Magnitudes')
    ax.legend()
    ax.grid(True, alpha=0.3)
    plt.tight_layout()
    
    return fig, ax
