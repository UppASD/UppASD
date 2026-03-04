"""
interactive.py

Interactive / user-facing helpers for constructing UppASD systems
and interactions.

This module:
- uses core data structures
- does not implement physics algorithms
- is suitable for CLI, GUI, and notebooks
"""

import numpy as np
from typing import Callable, Dict

from uppasd.core.system import SpinSystem
from uppasd.core.exchange import ExchangeShellTable, DMIShellTable
from uppasd.core.anisotropy import AnisotropyTable
from uppasd.core.neighbors import (
    build_neighbors,
    assign_shells,
)
from uppasd.input.inputdata import ASDInput
from uppasd.run.simulator import (
    ASDWorkspace,
    UppASDSimulator,
)


# ======================================================================
# System helpers
# ======================================================================


def make_system(
    cell,
    positions,
    species,
    moments,
):
    """
    Explicit system constructor.
    """
    return SpinSystem(
        cell=cell,
        positions=positions,
        species=species,
        moments=moments,
    )


# ======================================================================
# Exchange construction (interactive-friendly)
# ======================================================================


def build_exchange_from_geometry(
    system: SpinSystem,
    cutoff: float,
    Jij_model: Callable[[float], float],
    pbc=(True, True, True),
    tol: float = 1e-5,
):
    """
    Build ExchangeShellTable from geometry using NumPy-only neighbors.

    Intended for interactive exploration.
    """
    pairs, vectors, distances = build_neighbors(
        positions=system.positions,
        cell=system.cell,
        cutoff=cutoff,
        pbc=pbc,
    )

    shell_indices, shell_radii = assign_shells(distances, tol=tol)

    J = ExchangeShellTable()
    for (i, j), R, r, shell in zip(
        pairs, vectors, distances, shell_indices
    ):
        Jij = Jij_model(r)
        J.add_bond(i, j, shell, R, Jij)

    return J, shell_radii


def build_dmi_from_geometry(
    system: SpinSystem,
    cutoff: float,
    D_model: Callable[[np.ndarray], np.ndarray],
    pbc=(True, True, True),
    tol: float = 1e-5,
):
    """
    Build DMIShellTable from geometry.

    D_model(R_vec) -> D_vec
    """
    pairs, vectors, distances = build_neighbors(
        positions=system.positions,
        cell=system.cell,
        cutoff=cutoff,
        pbc=pbc,
    )

    shell_indices, shell_radii = assign_shells(distances, tol=tol)

    D = DMIShellTable()
    for (i, j), R, shell in zip(pairs, vectors, shell_indices):
        Dvec = np.asarray(D_model(R), dtype=float)
        D.add_bond(i, j, shell, R, Dvec)

    return D, shell_radii


# ======================================================================
# Interactive manipulation helpers
# ======================================================================


def scale_exchange_shell(
    J: ExchangeShellTable,
    shell: int,
    factor: float,
    iatom: int = None,
    jatom: int = None,
):
    """
    Scale exchange constants in a shell.
    """
    J.scale_shell(shell=shell, factor=factor, iatom=iatom, jatom=jatom)


def remove_exchange_shell(
    J: ExchangeShellTable,
    iatom: int,
    jatom: int,
    shell: int,
):
    """
    Remove an exchange shell.
    """
    J.remove_shell(iatom, jatom, shell)


# ======================================================================
# Execution helpers (interactive-friendly)
# ======================================================================


def run_interactive_relaxation(
    workdir: str,
    system: SpinSystem,
    exchange: ExchangeShellTable = None,
    inp: ASDInput = None,
    dmi: DMIShellTable = None,
    ani: AnisotropyTable = None,
    runtime: Dict = None,
    clean: bool = True,
):
    """
    Interactive-friendly relaxation wrapper.
    """
    ws = ASDWorkspace(workdir, clean=clean)
    ws.prepare(
        system=system,
        inp=inp,
        exchange=exchange,
        dmi=dmi,
        ani=ani,
    )

    sim = UppASDSimulator(ws)
    sim.initialize()

    if runtime:
        sim.set_runtime_controls(**runtime)

    sim.relax()
    sim.finalize()

    return ws


def run_interactive_measurement(
    workdir: str,
    system: SpinSystem,
    exchange: ExchangeShellTable = None,
    inp: ASDInput = None,
    dmi: DMIShellTable = None,
    ani: AnisotropyTable = None,
    runtime: Dict = None,
    clean: bool = False,
):
    """
    Interactive-friendly measurement wrapper.
    """
    ws = ASDWorkspace(workdir, clean=clean)
    ws.prepare(
        system=system,
        inp=inp,
        exchange=exchange,
        dmi=dmi,
        ani=ani,
    )

    sim = UppASDSimulator(ws)
    sim.initialize()

    if runtime:
        sim.set_runtime_controls(**runtime)

    sim.measure()
    sim.finalize()

    return ws
