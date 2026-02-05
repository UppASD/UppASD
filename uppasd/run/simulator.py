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
import io
from contextlib import contextmanager, redirect_stdout, redirect_stderr

from uppasd.input.inputdata import ASDInput
from uppasd.core.exchange import ExchangeShellTable, DMIShellTable
from uppasd.core.system import SpinSystem
from uppasd import pyasd

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

# ======================================================================
# Workspace
# ======================================================================


class ASDWorkspace:
    """Manage a filesystem workspace for an UppASD run.

    `ASDWorkspace` provides a small, explicit API for creating the run
    directory, writing the files expected by the Fortran backend (posfile,
    momfile, jfile, dmfile, inpsd.dat) and preparing a workspace ready for
    execution. It intentionally contains no physics logic — only filesystem
    orchestration.
    """

    def __init__(self, path: str, clean: bool = False):
        """Create or (optionally) recreate the workspace directory.

        Args:
            path (str): Path to the workspace directory. The path is resolved
                to an absolute `Path`.
            clean (bool): If True and the directory exists it will be
                removed recursively before creating a fresh workspace.
        """

        self.path = Path(path).resolve()

        if clean and self.path.exists():
            shutil.rmtree(self.path)

        self.path.mkdir(parents=True, exist_ok=True)

    # ------------------------------------------------------------------

    def write_system(self, system: SpinSystem):
        """Write system geometry and magnetic moments to workspace files.

        This writes the canonical `posfile` and `momfile` into the workspace
        by delegating to the provided `SpinSystem` helper which knows the
        expected on-disk layout.

        Args:
            system (SpinSystem): System object that implements
                `write_posfile(path)` and `write_momfile(path)`.
        """
        system.write_posfile(self.path / "posfile")
        system.write_momfile(self.path / "momfile")

    def write_interactions(
        self,
        exchange: ExchangeShellTable = None,
        dmi: DMIShellTable = None,
    ):
        """Write exchange (jfile) and DMI (dmfile) interaction tables.

        Only the files for tables that are provided are written. The method
        delegates to `ExchangeShellTable.write_jfile` and
        `DMIShellTable.write_dmfile` respectively.

        Args:
            exchange (ExchangeShellTable, optional): Exchange table to write.
            dmi (DMIShellTable, optional): DMI table to write.
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
        """Assemble and write the `inpsd.dat` input file.

        The function composes a canonical input file by starting from an
        empty `ASDInput` and merging blocks in a well-defined order so that
        structural and interaction file references are guaranteed to exist.

        Order of composition:
        1. system block (geometry and file references)
        2. interactions block (jfile/dmfile references)
        3. initialization block (default `initmag=3`)
        4. user-provided blocks from `inp` (merged or appended)

        Args:
            system (SpinSystem): System providing `input_block()` keys.
            inp (ASDInput): User-supplied input containing desired
                simulation/measurement blocks.
            exchange (ExchangeShellTable, optional): If present a reference
                to `./jfile` is inserted into the interactions block.
            dmi (DMIShellTable, optional): If present a reference to
                `./dmfile` is inserted into the interactions block.
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
        """Prepare the workspace by writing system, interactions and input.

        This is a convenience method that executes the sequence of file
        writes necessary before running the simulator: `write_system`,
        `write_interactions` and `write_input`.

        Args:
            system (SpinSystem): System definition.
            inp (ASDInput): User input blocks.
            exchange (ExchangeShellTable, optional): Exchange table.
            dmi (DMIShellTable, optional): DMI table.
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

    # The `UppASDSimulator` wraps the lower-level `pyasd` Fortran interface
    # and provides a small, explicit API for initialization, relaxation and
    # measurement. It handles temporary working-directory changes so callers
    # operate on an `ASDWorkspace` without managing CWD themselves.

    @contextmanager
    def _in_workspace(self):
        """Context manager that temporarily changes into the workspace.

        Use this to ensure all `pyasd` calls operate inside the run
        directory where input/output files are located. The original working
        directory is restored on exit even if an exception occurs.
        """

        cwd = os.getcwd()
        os.chdir(self.workspace.path)
        try:
            yield
        finally:
            os.chdir(cwd)

    def __init__(self, workspace: ASDWorkspace):
        """Create a simulator bound to the provided `ASDWorkspace`.

        The constructor does not perform Fortran-level initialization; call
        `initialize()` to set up the native simulation state.

        Args:
            workspace (ASDWorkspace): Workspace to use for file IO.
        """

        self.workspace = workspace

    # ------------------------------------------------------------------

    def initialize(self):
        """Initialize the Fortran simulator from the workspace files.

        This calls `pyasd.sanity_check()` and `pyasd.setup_all()` inside the
        workspace directory and stores `natom` and `mensemble` on the
        simulator object for later calls.

        Returns:
            None

        Raises:
            RuntimeError: If setup fails (propagates exceptions from `pyasd`).
        """

        cwd = os.getcwd()
        os.chdir(self.workspace.path)
        try:
            with _suppress_output():
                pyasd.sanity_check()
                self.natom, self.mensemble = pyasd.setup_all()
        finally:
            os.chdir(cwd)

    # ------------------------------------------------------------------

    def set_runtime_controls(self, **kwargs):
        """Set runtime-controllable parameters exposed by `pyasd`.

        This utility finds `pyasd.set_<key>` callables for each provided
        keyword and forwards the value. It raises an AttributeError if the
        expected setter is not present which helps detect typos early.

        Args:
            **kwargs: Keyword arguments where the key name corresponds to the
                runtime parameter (e.g. `temperature`, `timestep`, `damping`).

        Raises:
            AttributeError: If `pyasd` does not expose a `set_<key>` method.
        """

        for key, value in kwargs.items():
            setter = getattr(pyasd, f"set_{key}", None)
            if setter is None:
                raise AttributeError(f"pyasd has no setter for '{key}'")
            setter(value)

    # ------------------------------------------------------------------

    def relax(self, mode='S', nstep=10, temperature=None, timestep=None, damping=None):
        """Run a relaxation (LLG or Monte Carlo) and return updated moments.

        This is the main execution method used by higher-level tools. It
        validates that the simulator was initialized, assembles a compact
        kwargs dictionary containing only the provided runtime parameters
        and invokes `pyasd.relax` from within the workspace.

        Args:
            mode (str): One of 'S' (LLG), 'M' (Metropolis MC) or 'H' (Heat
                bath). Case-insensitive.
            nstep (int): Number of relaxation steps to perform.
            temperature (float, optional): Temperature to use for the run.
            timestep (float, optional): Integration timestep.
            damping (float, optional): Damping parameter.

        Returns:
            numpy.ndarray: The magnetic moments after relaxation with shape
                ``(3, natom, mensemble)``.

        Raises:
            RuntimeError: If the simulator has not been initialized.
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

            with _suppress_output():
                #     pyasd.relax(self.natom, self.mensemble, **relax_kwargs)
                moments = pyasd.relax(self.natom, self.mensemble, **relax_kwargs)

        return moments

    def relax_mc(self, **kwargs):
        """Run a Monte Carlo-style relaxation via `pyasd.relax`.

        This convenience wrapper forwards arbitrary keyword arguments to
        `pyasd.relax` and returns the updated moments. It performs the same
        initialization check and workspace switch as `relax`.

        Args:
            **kwargs: Forwarded to `pyasd.relax` (e.g. `mode='M'`,
                `nstep=100`, `temperature=10.0`).

        Returns:
            numpy.ndarray: Magnetic moments after relaxation.

        Raises:
            RuntimeError: If the simulator has not been initialized.
        """

        with self._in_workspace():
            if self.natom is None or self.mensemble is None:
                raise RuntimeError("Simulator not initialized")

            with _suppress_output():
                moments = pyasd.relax(self.natom, self.mensemble, **kwargs)

        return moments

    def measure(self):
        """Trigger measurement routines in the Fortran backend.

        This calls `pyasd.measure()` inside the workspace. Measurements write
        output files that can later be parsed by `ASDResults`.
        """

        with self._in_workspace():
            with _suppress_output():
                pyasd.measure()

    def run_init_phase(self):
        """Run the initialization phase (compute fields, prepare state).

        The initial phase typically computes magnetic fields and other
        quantities required before time evolution. It maps to the Fortran
        `initial_phase` driver.
        """

        with self._in_workspace():
            with _suppress_output():
                pyasd.initial_phase()

    def finalize(self):
        """Perform simulator cleanup and release Fortran-side resources.

        Calls `pyasd.cleanup()` in the workspace. This should be invoked after
        any full-run or when the Python process no longer needs the Fortran
        runtime state.
        """

        with self._in_workspace():
            with _suppress_output():
                pyasd.cleanup()

    def run_all(self):
        """
        Execute a complete UppASD simulation from initialization through cleanup.

        This is a convenience wrapper that performs:
        1. Initialization from input files (setup_all)
        2. Initial phase computation
        3. Main simulation loop (run_uppasd)
        4. Measurement output
        5. Cleanup

        This method encapsulates the standard full simulation workflow and
        is equivalent to calling initialize() → run_init_phase() → relax() → measure() → finalize()
        in sequence, but uses pyasd.run_uppasd() for the main integration loop.

        Returns
        -------
        None

        Raises
        ------
        RuntimeError
            If simulation fails or workspace is invalid

        Examples
        --------
        >>> workspace = ASDWorkspace('./run_dir', clean=True)
        >>> workspace.prepare(system=system, inp=inp, exchange=exchange)
        >>> simulator = UppASDSimulator(workspace)
        >>> simulator.run_all()  # Runs complete simulation
        >>> results = ASDResults('./run_dir')  # Load output
        """
        with self._in_workspace():
            try:
                with _suppress_output():
                    # Step 1: Initialization
                    pyasd.sanity_check()
                    self.natom, self.mensemble = pyasd.setup_all()

                    # Step 2: Initial phase (compute fields, etc.)
                    pyasd.initial_phase()

                    # Step 4: Measurement output
                    pyasd.measure()

            finally:
                # Step 5: Cleanup (always runs, even if error occurs)
                with _suppress_output():
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
    """Run a standard relaxation workflow in a transient workspace.

    This convenience helper creates an `ASDWorkspace` at ``workdir``, writes
    the required system and interaction files, initializes the Fortran
    simulator, applies any runtime controls provided via ``runtime`` and
    executes a relaxation (via :py:meth:`UppASDSimulator.relax`). The
    workspace is returned for downstream inspection or result parsing.

    Args:
        workdir (str): Filesystem path for the temporary run directory.
        system (SpinSystem): System describing positions and moments.
        inp (ASDInput): User input blocks describing simulation/measurement.
        exchange (ExchangeShellTable, optional): Exchange table to write.
        dmi (DMIShellTable, optional): DMI table to write.
        runtime (dict, optional): Runtime parameters forwarded to
            :py:meth:`UppASDSimulator.set_runtime_controls` (e.g.
            ``{'temperature': 10.0, 'timestep': 1e-15}``).
        clean (bool): If True the workspace is recreated; existing directory
            will be removed.

    Returns:
        ASDWorkspace: The prepared workspace containing input and output
            files produced by the run.

    Raises:
        Any exception raised by the underlying filesystem operations or the
        Fortran backend is propagated to the caller.

    Example:
        >>> ws = run_relaxation('./run_relax', system, inp, runtime={'temperature':10.0})
        >>> results = ASDResults(ws.path)
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
    """Run a measurement-only workflow and return the workspace.

    This helper mirrors :func:`run_relaxation` but only performs the
    measurement phase. It prepares the workspace, initializes the Fortran
    simulator, applies runtime controls if provided, and calls
    :py:meth:`UppASDSimulator.measure` before finalizing the Fortran state.

    Args:
        workdir (str): Path for the run directory.
        system (SpinSystem): System describing positions and moments.
        inp (ASDInput): Input blocks specifying measurements to perform.
        exchange (ExchangeShellTable, optional): Exchange table to write.
        dmi (DMIShellTable, optional): DMI table to write.
        runtime (dict, optional): Runtime parameters forwarded to
            :py:meth:`UppASDSimulator.set_runtime_controls`.
        clean (bool): If True recreate the workspace directory.

    Returns:
        ASDWorkspace: Workspace populated with measurement output files.

    Raises:
        Propagates exceptions from file IO and the Fortran backend.

    Example:
        >>> ws = run_measurement('./run_meas', system, inp)
        >>> # parse outputs from ws.path
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
