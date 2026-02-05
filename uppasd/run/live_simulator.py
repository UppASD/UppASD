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

from typing import Optional

import numpy as np

from uppasd import pyasd
from uppasd.run.simulator import ASDWorkspace, UppASDSimulator


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
    
    # The LiveSimulator encapsulates a short-lived `UppASDSimulator` instance
    # and exposes a thin, convenient API useful for interactive tools and
    # notebooks. Methods are intentionally lightweight wrappers around the
    # lower-level `pyasd` calls so that the Fortran-backed simulation code
    # remains authoritative for physics and defaults.

    # ------------------------------------------------------------------
    # Construction
    # ------------------------------------------------------------------

    def __init__(self, workspace: ASDWorkspace):
        """Create a LiveSimulator bound to an `ASDWorkspace`.

        Args:
            workspace (ASDWorkspace): Workspace describing the run directory
                and file layout used by the underlying simulator.

        Notes:
            Construction does not start or initialize the Fortran runtime;
            call `initialize()` before attempting to step or query state.
        """
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
        """Initialize the underlying simulator and cache initial state.

        This routine performs the minimal startup required to run
        incremental relaxations. It constructs an `UppASDSimulator`, calls
        its initialization sequence and queries atomic / ensemble sizes.

        The initial magnetic moments are cached so the user can restore the
        starting configuration later via `reset_moments()`.

        Returns:
            None

        Raises:
            RuntimeError: If initialization fails at the simulator layer.
        """
        if self._initialized:
            return

        # Suppress stdout/stderr from underlying Fortran/Python initialization
        self._sim = UppASDSimulator(self.workspace)
        self._sim.initialize()

        self.natom = self._sim.natom
        self.mensemble = self._sim.mensemble

        # Cache initial moments for reset
        self._emom0 = pyasd.get_emom(self.natom, self.mensemble).copy()
        pyasd.put_emom(pyasd.get_emom(self.natom, self.mensemble), self.natom, self.mensemble)
        self._hfield0 = pyasd.get_iphfield().copy()
        self._iptemperature0 = pyasd.get_iptemperature()
        

        self._initialized = True

    # ------------------------------------------------------------------
    # Internal guard
    # ------------------------------------------------------------------

    def _check_initialized(self):
        """Internal guard which raises if the simulator is not ready.

        Many public methods require the LiveSimulator to be initialized.
        This helper centralizes the check and the error message.
        """
        if not self._initialized:
            raise RuntimeError("LiveSimulator not initialized")

    # ------------------------------------------------------------------
    # Unified stepping interface (PRIMARY API)
    # ------------------------------------------------------------------

    def step(self, mode: str, temperature=1.0e-6, nstep: int = 1):
        """Advance the simulation by a small number of steps.

        This is the primary stepping API used by interactive tools. It calls
        `pyasd.relax(...)` with the currently cached sizes and the supplied
        mode/temperature. The function enforces simple validation and
        safeguards (for example a small non-zero temperature for Heat-bath
        MC) to avoid undefined behaviour in the Fortran layer.

        Args:
            mode (str): One-letter mode specifier: ``"S"`` (LLG), ``"M"``
                (Metropolis MC) or ``"H"`` (Heat-bath MC). Case-insensitive.
            temperature (float): Temperature to use for this step. If
                ``None`` the simulator's current `iptemperature` will be
                queried and used.
            nstep (int): Number of relaxation steps to perform.

        Returns:
            None

        Raises:
            ValueError: If `mode` is not one of ``"S"``, ``"M"``, ``"H"``.
        """
        self._check_initialized()
        # print('Step initialized with mode:', mode, 'and temperature:', temperature)

        mode = mode.upper()
        if mode not in ("S", "M", "H"):
            raise ValueError(f"Unknown relaxation mode '{mode}'")

        if nstep <= 0:
            return

        # Default temperature to current state
        # print('------------->', temperature, mode)
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
        # print(f"Stepping {nstep} steps via mode '{mode}'")
        # before calling pyasd.relax
        # pyasd.put_emom(self.get_emom(), self.natom, self.mensemble)

        moments = pyasd.relax(
            natom=self.natom,
            mensemble=self.mensemble,
            mode=mode,
            temperature=temperature,
            nstep=nstep,
        )[:, :, 0].astype(np.float64)
        
        return moments

    # ------------------------------------------------------------------
    # State accessors
    # ------------------------------------------------------------------

    def get_emom(self) -> np.ndarray:
        """Return the current magnetic moment array from the simulator.

        The returned array has shape ``(3, natom, mensemble)`` where the first
        axis indexes Cartesian components (x,y,z). The method delegates to
        ``pyasd.get_emom`` and therefore returns fresh data from the Fortran
        runtime on every call.

        Returns:
            numpy.ndarray: Moments in shape ``(3, natom, mensemble)``.
        """
        self._check_initialized()
        return pyasd.get_emom(self.natom, self.mensemble)

    def get_coords(self) -> np.ndarray:
        """Return atomic coordinates as a 2D NumPy array.

        Returns:
            numpy.ndarray: Coordinates with shape ``(3, natom)``.
        """
        self._check_initialized()
        return pyasd.get_coords(self.natom)

    def reset_moments(self):
        """Restore the cached initial magnetic moments saved at
        initialization.

        This writes the originally cached moment array back into the running
        simulator via ``pyasd.put_emom``. If the simulator was not
        initialized or no cache exists the call is a no-op.

        Returns:
            None
        """
        self._check_initialized()
        if self._emom0 is not None:
            pyasd.put_emom(self._emom0, self.natom, self.mensemble)

    # ------------------------------------------------------------------
    # Parameters (thin wrappers)
    # ------------------------------------------------------------------

    def get_temperature(self) -> float:
        """Return the simulator's current initial-phase temperature.

        This calls ``pyasd.get_iptemperature()`` and returns the floating
        point temperature used by the initial-phase drivers.

        Returns:
            float: Current temperature in simulation units (Kelvin-equivalent
                in the user's conceptual model).
        """
        return pyasd.get_iptemperature()

    def set_temperature(self, T: float):
        """Set the simulator's initial-phase temperature.

        This is a thin wrapper around ``pyasd.put_iptemperature``. The
        argument is cast to ``float`` before being forwarded to the Fortran
        layer.

        Args:
            T (float): Temperature value to set.

        Returns:
            None
        """
        #pyasd.put_iptemperature(float(T))
        pyasd.set_temperature(float(T))

    def get_field(self) -> np.ndarray:
        """Return the current initial-phase magnetic field vector.

        The returned value is a NumPy copy of the internal field vector and
        has shape ``(3,)``.

        Returns:
            numpy.ndarray: Field vector [Bx, By, Bz].
        """
        return pyasd.get_iphfield().copy()

    def set_field(self, B: np.ndarray):
        """Set the initial-phase magnetic field vector used by the simulator.

        The provided array is converted to a float NumPy array and forwarded
        to ``pyasd.set_iphfield``. The function accepts any sequence-like
        object of length three.

        Args:
            B (array-like): Magnetic field vector-like object convertible to
                ``np.ndarray`` with shape ``(3,)``.

        Returns:
            None
        """
        pyasd.set_iphfield(np.asarray(B, dtype=float))
        pyasd.set_hfield(np.asarray(B, dtype=float))

    # ------------------------------------------------------------------
    # Energy
    # ------------------------------------------------------------------

    def calculate_energy(self) -> float:
        """Evaluate and return the current total energy per atom.

        This forces an energy update inside the Fortran runtime (via
        ``pyasd.get_energy()``) and returns the scalar energy value. The
        returned energy is in the same units used elsewhere in the project
        (milli-Rydberg per atom in this codebase).

        Returns:
            float: Total energy per atom.
        """
        self._check_initialized()
        pyasd.get_energy()
        return pyasd.get_energy()

    @property
    def initialized(self) -> bool:
        return self._initialized
