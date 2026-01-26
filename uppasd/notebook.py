"""
Notebook-oriented helper utilities for UppASD.

This module keeps plotting out of the core package and provides
lightweight helpers for notebooks:
- read_inpsd: parse minimal metadata from inpsd.dat
- run_simulation_via_api: run UppASD via the Python API (no subprocess)
- load_outputs: convenience loader for averages/totenergy outputs
- iteration_to_time: convert iteration indices to physical time
- run_simulation_api: logging wrapper around run_simulation_via_api
- load_simulation_results: load outputs using a simid or config
- read_generic_output_file: parse text outputs and headers
- find_output_file: locate output files by simid/pattern
- print_system_info: pretty-print parsed input parameters
- setup_magnon_q_mesh: create q-point meshes for magnon calculations (wraps magnons module)
- compute_magnon_dispersion: run LSWT magnon calculation in notebook

Dependencies are kept minimal (numpy, glob, logging) to avoid heavy
plotting requirements inside the package.
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
