"""
notebook.py

Notebook helper utilities for UppASD.

This module is:
- OPTIONAL
- NOTEBOOK-ONLY
- NON-AUTHORITATIVE

It provides ergonomic helpers for:
- defining common systems
- building exchange / DMI interactions
- running standard workflows
- quick visualization hooks

Core logic lives elsewhere.
"""

import numpy as np

# Core API (authoritative)
from core.system import SpinSystem
from core.exchange import ExchangeShellTable, DMIShellTable
from inputdata import ASDInput
from simulator import run_relaxation, run_measurement

# Optional dependencies (lazy import!)
try:
    from ase.neighborlist import NeighborList
    _HAVE_ASE = True
except ImportError:
    _HAVE_ASE = False


# ======================================================================
# System helpers
# ======================================================================


def simple_cubic_system(
    a: float,
    species,
    moments,
):
    """
    Simple cubic system helper.

    species : list[str]
    moments : list[[mx,my,mz]]
    """
    positions = [[0.0, 0.0, 0.0]]
    cell = [[a, 0, 0], [0, a, 0], [0, 0, a]]

    return SpinSystem(
        cell=cell,
        positions=positions,
        species=species,
        moments=moments,
    )


def from_arrays(
    cell,
    positions,
    species,
    moments,
):
    """
    Minimal helper for quick notebook prototyping.
    """
    return SpinSystem(
        cell=cell,
        positions=positions,
        species=species,
        moments=moments,
    )


# ======================================================================
# Exchange / DMI helpers
# ======================================================================


def heisenberg_shells(
    shells,
):
    """
    Construct ExchangeShellTable from an explicit shell list.

    shells = [
        dict(
            iatom=1,
            jatom=1,
            shell=1,
            vectors=[[1,0,0],[0,1,0]],
            Jij=5.0,
        ),
        ...
    ]
    """
    J = ExchangeShellTable()
    for sh in shells:
        J.add_shell(
            iatom=sh["iatom"],
            jatom=sh["jatom"],
            shell=sh["shell"],
            vectors=sh["vectors"],
            Jij=sh["Jij"],
        )
    return J


def dmi_shells(
    shells,
):
    """
    Construct DMIShellTable from explicit shell definitions.
    """
    D = DMIShellTable()
    for sh in shells:
        D.add_shell(
            iatom=sh["iatom"],
            jatom=sh["jatom"],
            shell=sh["shell"],
            vectors=sh["vectors"],
            D_vectors=sh["D"],
        )
    return D


# ======================================================================
# ASE-based neighbor helpers (optional)
# ======================================================================


def exchange_from_ase(
    atoms,
    cutoff: float,
    shell_map,
    Jij_model,
):
    """
    Build ExchangeShellTable using ASE neighbor lists.

    shell_map: callable(distance) -> shell_index
    Jij_model: callable(distance) -> Jij
    """
    if not _HAVE_ASE:
        raise ImportError("ASE is required for exchange_from_ase")

    positions = atoms.get_positions()
    cell = atoms.get_cell()
    natom = len(atoms)

    cutoffs = [cutoff] * natom
    nl = NeighborList(cutoffs, self_interaction=False, bothways=True)
    nl.update(atoms)

    J = ExchangeShellTable()

    for i in range(natom):
        indices, offsets = nl.get_neighbors(i)
        for j, offset in zip(indices, offsets):
            R = offset @ cell
            r = np.linalg.norm(R)
            shell = shell_map(r)
            Jij = Jij_model(r)
            J.add_bond(i + 1, j + 1, shell, R, Jij)

    return J


# ======================================================================
# Input helpers
# ======================================================================


def basic_llg_input(
    timestep=1e-15,
    damping=0.01,
    nstep=10000,
    temperature=None,
):
    """
    Very common LLG input setup.
    """
    inp = ASDInput()
    sim = inp.block("simulation")
    sim.set(
        timestep=timestep,
        damping=damping,
        nstep=nstep,
    )
    if temperature is not None:
        sim.set(temperature=temperature)
    return inp


# ======================================================================
# Workflow helpers
# ======================================================================


def relax_and_measure(
    workdir,
    system,
    exchange,
    inp,
    runtime=None,
    clean=True,
):
    """
    One-liner workflow for notebooks.
    """
    run_relaxation(
        workdir=workdir,
        system=system,
        inp=inp,
        exchange=exchange,
        runtime=runtime,
        clean=clean,
    )

    run_measurement(
        workdir=workdir,
        system=system,
        inp=inp,
        exchange=exchange,
        runtime=runtime,
        clean=False,
    )


# ======================================================================
# Visualization hooks (placeholders)
# ======================================================================


def plot_spins(*args, **kwargs):
    """
    Placeholder hook for PyVista-based visualization.

    Implemented in uppasd.notebook.viz
    """
    raise NotImplementedError(
        "Visualization helpers live in uppasd.notebook.viz"
    )
