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

    def write_input(
        self,
        system: SpinSystem,
        inp: ASDInput,
        exchange: ExchangeShellTable = None,
        dmi: DMIShellTable = None,
    ):
        """
        Write inpsd.dat by merging system, interactions, initialization, and user input.
        
        Order:
        1. System block (geometry, files)
        2. Interactions block (exchange/DMI file references)
        3. Initialization block (with initmag default)
        4. User blocks (simulation, measurement, other keywords)
        """
        full = ASDInput()

        # System block: geometry and file references
        full.block("system").set(**system.input_block())

        # Interactions block: file references only
        if exchange is not None:
            full.block("interactions").set(exchange="./jfile")
        if dmi is not None:
            full.block("interactions").set(dm="./dmfile")

        # Initialization block: default initmag 3
        full.block("initialization").set(initmag=3)

        # User-provided blocks (simulation, measurement, interactions keywords, etc.)
        for name, block in inp.blocks.items():
            if name == "system":
                # Merge with system block if user provided system keys
                full.block("system").set(**block.as_dict())
            elif name == "interactions":
                # Merge with interactions block if user provided interaction keys
                full.block("interactions").set(**block.as_dict())
            elif name == "initialization":
                # Merge with initialization block, user can override initmag
                full.block("initialization").set(**block.as_dict())
            else:
                # Other blocks (simulation, measurement, etc.) pass through
                full.add_block(name, block)

        full.write(self.path / "inpsd.dat")

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
        self.write_input(system, inp, exchange, dmi)


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

    def relax(self, mode='S', nstep=10, temperature=None, timestep=None, damping=None):
        """
        Run a relaxation with specified method and parameters.
        
        Parameters
        ----------
        mode : {'S', 'M', 'H'}, optional
            Relaxation method:
            - 'S': Spin-dynamics (LLG, default)
            - 'M': Metropolis (single-spin flip MC)
            - 'H': Heat bath
            Default: 'S'
        nstep : int, optional
            Number of relaxation steps. Default: 10
        temperature : float, optional
            Temperature in Kelvin. If None, uses current Fortran state.
            Default: None
        timestep : float, optional
            Integration timestep in seconds. If None, uses current Fortran state.
            Default: None
        damping : float, optional
            Damping factor (0 to 1). If None, uses current Fortran state.
            Default: None
        """
        with self._in_workspace():
            if self.natom is None or self.mensemble is None:
                raise RuntimeError("Simulator not initialized")

            # Build kwargs with only provided parameters
            relax_kwargs = {
                'mode': mode,
                'nstep': nstep,
            }
            if temperature is not None:
                relax_kwargs['temperature'] = temperature
            if timestep is not None:
                relax_kwargs['timestep'] = timestep
            if damping is not None:
                relax_kwargs['damping'] = damping

            pyasd.relax(self.natom, self.mensemble, **relax_kwargs)

    def relax_mc(self, **kwargs):
        with self._in_workspace():
            if self.natom is None or self.mensemble is None:
                raise RuntimeError("Simulator not initialized")

            pyasd.relax(self.natom, self.mensemble, **kwargs)

    def measure(self):
        with self._in_workspace():
            pyasd.measure()

    def run_init_phase(self):
        with self._in_workspace():
            pyasd.initial_phase()

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
