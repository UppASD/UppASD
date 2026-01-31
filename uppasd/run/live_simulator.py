"""
live_simulator.py

Stateful, context-preserving live simulator for UppASD.

This class is designed for:
- interactive visualization (VTK / notebooks)
- live parameter steering
- incremental relaxation with full UppASD context

All stepping is done via pyasd.relax(), NOT algorithm-specific calls.

Author: UppASD
"""

from __future__ import annotations

import io
import os
from contextlib import contextmanager, redirect_stderr, redirect_stdout
from typing import Optional

import numpy as np

from uppasd import pyasd
from uppasd.run.simulator import ASDWorkspace, UppASDSimulator


@contextmanager
def _suppress_output():
    """Suppress OS-level and Python-level stdout/stderr for external libs.

    Redirects file descriptors 1 and 2 to /dev/null and also redirects
    Python's sys.stdout/sys.stderr so that both low-level and Python
    printed text are silenced inside the context.
    """
    devnull = os.open(os.devnull, os.O_RDWR)
    old_stdout_fd = os.dup(1)
    old_stderr_fd = os.dup(2)
    try:
        os.dup2(devnull, 1)
        os.dup2(devnull, 2)

        sio_out = io.StringIO()
        sio_err = io.StringIO()
        with redirect_stdout(sio_out), redirect_stderr(sio_err):
            yield
    finally:
        os.dup2(old_stdout_fd, 1)
        os.dup2(old_stderr_fd, 2)
        os.close(devnull)
        os.close(old_stdout_fd)
        os.close(old_stderr_fd)


class LiveSimulator:
    """
    Live, stateful UppASD simulator.

    Intended lifecycle:
        ws = ASDWorkspace(...)
        sim = LiveSimulator(ws)
        sim.initialize()
        sim.step("S", nstep=10)
        emom = sim.get_emom()
    """

    # ------------------------------------------------------------------
    # Construction
    # ------------------------------------------------------------------

    def __init__(self, workspace: ASDWorkspace):
        self.workspace = workspace

        self._sim: Optional[UppASDSimulator] = None

        self.natom: Optional[int] = None
        self.mensemble: Optional[int] = None

        self._initialized: bool = False
        self._emom0: Optional[np.ndarray] = None

    # ------------------------------------------------------------------
    # Initialization
    # ------------------------------------------------------------------

    def initialize(self):
        """
        Initialize UppASD and cache initial magnetic moments.
        """
        if self._initialized:
            return

        # Suppress stdout/stderr from underlying Fortran/Python initialization
        with _suppress_output():
            self._sim = UppASDSimulator(self.workspace)
            self._sim.initialize()

            self.natom = self._sim.natom
            self.mensemble = self._sim.mensemble

            # Cache initial moments for reset
            self._emom0 = pyasd.get_emom(self.natom, self.mensemble).copy()

            self._initialized = True

    # ------------------------------------------------------------------
    # Internal guard
    # ------------------------------------------------------------------

    def _check_initialized(self):
        if not self._initialized:
            raise RuntimeError("LiveSimulator not initialized")

    # ------------------------------------------------------------------
    # Unified stepping interface (PRIMARY API)
    # ------------------------------------------------------------------

    def step(self, mode: str, temperature=None, nstep: int = 1):
        """
        Advance the simulation by nstep steps using pyasd.relax().

        Parameters
        ----------
        mode : str
            "S" : LLG
            "M" : Metropolis Monte Carlo
            "H" : Heat-bath Monte Carlo
        nstep : int
            Number of steps to perform
        """
        self._check_initialized()

        mode = mode.upper()
        if mode not in ("S", "M", "H"):
            raise ValueError(f"Unknown relaxation mode '{mode}'")

        if nstep <= 0:
            return

        # Default temperature to current state
        if temperature is None:
            temperature = pyasd.get_iptemperature()

        # Safeguard against zero temperature in Heat Bath MC
        if mode == "H" and temperature < 1e-6:
            temperature = 1e-6

        # IMPORTANT:
        # pyasd.relax handles:
        # - ipmode
        # - initial phase
        # - correct MC/LLG kernels
        # - state continuity
        print(f"Stepping {nstep} steps via mode '{mode}'")
        pyasd.relax(
            natom=self.natom,
            mensemble=self.mensemble,
            mode=mode,
            temperature=temperature,
            nstep=nstep,
        )

    # ------------------------------------------------------------------
    # State accessors
    # ------------------------------------------------------------------

    def get_emom(self) -> np.ndarray:
        """
        Return magnetic moments.

        Shape: (3, natom, mensemble)
        """
        self._check_initialized()
        return pyasd.get_emom(self.natom, self.mensemble)

    def get_coords(self) -> np.ndarray:
        """
        Return atomic coordinates.

        Shape: (3, natom)
        """
        self._check_initialized()
        return pyasd.get_coords(self.natom)

    def reset_moments(self):
        """
        Restore cached initial magnetic moments.
        """
        self._check_initialized()
        if self._emom0 is not None:
            pyasd.put_emom(self._emom0, self.natom, self.mensemble)

    # ------------------------------------------------------------------
    # Parameters (thin wrappers)
    # ------------------------------------------------------------------

    def get_temperature(self) -> float:
        return pyasd.get_iptemperature()

    def set_temperature(self, T: float):
        pyasd.put_iptemperature(float(T))

    def get_field(self) -> np.ndarray:
        return pyasd.get_iphfield().copy()

    def set_field(self, B: np.ndarray):
        pyasd.set_iphfield(np.asarray(B, dtype=float))

    # ------------------------------------------------------------------
    # Energy
    # ------------------------------------------------------------------

    def calculate_energy(self) -> float:
        """
        Force energy evaluation and return total energy per atom.
        """
        self._check_initialized()
        pyasd.get_energy()
        return pyasd.get_energy()
