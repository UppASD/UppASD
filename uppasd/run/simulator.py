"""
simulator.py

High-level orchestration of UppASD simulations.

Responsibilities:
- Workspace management
- Writing structural input files
- Calling pyasd for execution
- NO physics logic
- NO defaults
"""

from pathlib import Path
import shutil
import os

from contextlib import contextmanager
from uppasd.input.inputdata import ASDInput
from uppasd.core.exchange import ExchangeShellTable, DMIShellTable
from uppasd.core.system import SpinSystem
from uppasd import pyasd


# ======================================================================
# Workspace
# ======================================================================


class ASDWorkspace:
    """
    A directory-based workspace for UppASD runs.
    """

    def __init__(self, path: str, clean: bool = False):
        self.path = Path(path).resolve()

        if clean and self.path.exists():
            shutil.rmtree(self.path)

        self.path.mkdir(parents=True, exist_ok=True)

    # ------------------------------------------------------------------

    def write_system(self, system: SpinSystem):
        """
        Write posfile and momfile.
        """
        system.write_posfile(self.path / "posfile")
        system.write_momfile(self.path / "momfile")

    def write_interactions(
        self,
        exchange: ExchangeShellTable = None,
        dmi: DMIShellTable = None,
    ):
        """
        Write jfile and dmfile if provided.
        """
        if exchange is not None:
            exchange.write_jfile(self.path / "jfile")

        if dmi is not None:
            dmi.write_dmfile(self.path / "dmfile")

    def write_input(self, inp: ASDInput):
        """
        Write inpsd.dat.
        """
        inp.write(self.path / "inpsd.dat")

    # ------------------------------------------------------------------

    def prepare(
        self,
        system: SpinSystem,
        inp: ASDInput,
        exchange: ExchangeShellTable = None,
        dmi: DMIShellTable = None,
    ):
        """
        Prepare a full UppASD workspace.
        """
        self.write_system(system)
        self.write_interactions(exchange, dmi)
        self.write_input(inp)


# ======================================================================
# Execution
# ======================================================================


class UppASDSimulator:
    """
    Thin execution wrapper around pyasd.
    """

    @contextmanager
    def _in_workspace(self):
        cwd = os.getcwd()
        os.chdir(self.workspace.path)
        try:
            yield
        finally:
            os.chdir(cwd)

    def __init__(self, workspace: ASDWorkspace):
        self.workspace = workspace

    # ------------------------------------------------------------------

    def initialize(self):
        """
        Initialize UppASD from files.
        self.natom: number of atoms
        self.mensemble: ensemble size
        """
        cwd = os.getcwd()
        os.chdir(self.workspace.path)
        try:
            pyasd.sanity_check()
            self.natom, self.mensemble = pyasd.setup_all()
        finally:
            os.chdir(cwd)

    # ------------------------------------------------------------------

    def set_runtime_controls(self, **kwargs):
        """
        Set runtime-controllable parameters via pyasd.

        Examples:
            temperature
            timestep
            nstep
            field
            damping
        """
        for key, value in kwargs.items():
            setter = getattr(pyasd, f"set_{key}", None)
            if setter is None:
                raise AttributeError(f"pyasd has no setter for '{key}'")
            setter(value)

    # ------------------------------------------------------------------
    def relax(self):
        with self._in_workspace():
            pyasd.relax_llg()

    def relax_mc(self):
        with self._in_workspace():
            pyasd.relax_mc()

    def measure(self):
        with self._in_workspace():
            pyasd.measure()

    def finalize(self):
        with self._in_workspace():
            pyasd.cleanup()


# ======================================================================
# Convenience workflows
# ======================================================================


def run_relaxation(
    workdir: str,
    system: SpinSystem,
    inp: ASDInput,
    exchange: ExchangeShellTable = None,
    dmi: DMIShellTable = None,
    runtime: dict = None,
    clean: bool = True,
):
    """
    Convenience wrapper for a standard relaxation run.
    """
    ws = ASDWorkspace(workdir, clean=clean)
    ws.prepare(system, inp, exchange, dmi)

    sim = UppASDSimulator(ws)
    sim.initialize()

    if runtime:
        sim.set_runtime_controls(**runtime)

    sim.relax()
    sim.finalize()

    return ws


def run_measurement(
    workdir: str,
    system: SpinSystem,
    inp: ASDInput,
    exchange: ExchangeShellTable = None,
    dmi: DMIShellTable = None,
    runtime: dict = None,
    clean: bool = True,
):
    """
    Convenience wrapper for a measurement run.
    """
    ws = ASDWorkspace(workdir, clean=clean)
    ws.prepare(system, inp, exchange, dmi)

    sim = UppASDSimulator(ws)
    sim.initialize()

    if runtime:
        sim.set_runtime_controls(**runtime)

    sim.measure()
    sim.finalize()

    return ws
