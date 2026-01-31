"""
NotebookLiveViewer
==================

Jupyter-native interactive viewer for UppASD LiveSimulator.

Design goals:
- NO VTK interactor loops
- NO physics logic
- Works in standard Jupyter / Colab / Binder
- Uses LiveSimulator stepping contract
- Safe repeated stepping (LLG / MC / future modes)

Author: Anders Bergman + ChatGPT
"""

from __future__ import annotations

import numpy as np
import matplotlib.pyplot as plt
import os
import sys
import io
from contextlib import contextmanager, redirect_stdout, redirect_stderr

from dataclasses import dataclass
from typing import Optional

import ipywidgets as widgets
from IPython.display import display, clear_output

# ----------------------------------------------------------------------
# Helper dataclass for visualization settings
# ----------------------------------------------------------------------

@dataclass
class ViewerConfig:
    scale: float = 1.0
    stride: int = 1
    cmap: str = "coolwarm"
    figsize: tuple = (6, 6)


@contextmanager
def _suppress_output():
    """Context manager that suppresses stdout and stderr at the OS level.

    This redirects file descriptors 1 and 2 to os.devnull so external
    libraries (Fortran/C) writing to stdout/stderr do not print to the
    notebook output cell.
    """
    devnull = os.open(os.devnull, os.O_RDWR)
    # Save original file descriptors
    old_stdout_fd = os.dup(1)
    old_stderr_fd = os.dup(2)
    try:
        # Redirect low-level fds
        os.dup2(devnull, 1)
        os.dup2(devnull, 2)

        # Also suppress Python-level sys.stdout/sys.stderr
        sio_out = io.StringIO()
        sio_err = io.StringIO()
        with redirect_stdout(sio_out), redirect_stderr(sio_err):
            yield
    finally:
        # Restore original fds
        os.dup2(old_stdout_fd, 1)
        os.dup2(old_stderr_fd, 2)
        os.close(devnull)
        os.close(old_stdout_fd)
        os.close(old_stderr_fd)


# ----------------------------------------------------------------------
# NotebookLiveViewer
# ----------------------------------------------------------------------

class NotebookLiveViewer:
    """
    Jupyter-based live viewer for UppASD LiveSimulator.

    Responsibilities:
    - UI controls (widgets)
    - Visualization
    - Delegating stepping to LiveSimulator

    NO physics logic here.
    """

    # ------------------------------------------------------------------
    def __init__(self, sim, config: Optional[ViewerConfig] = None):
        """
        Parameters
        ----------
        sim : LiveSimulator
            Initialized LiveSimulator instance.
        config : ViewerConfig, optional
            Visualization configuration.
        """
        self.sim = sim
        self.config = config or ViewerConfig()

        # Cache temperature and field for step calls
        self.temperature = self.sim.get_temperature()
        self.b_z = self.sim.get_field()[2]

        # Create figure upfront in non-interactive mode to prevent auto-display
        plt.ioff()
        self._fig, self._ax = plt.subplots(figsize=self.config.figsize)
        plt.ion()

        self._build_widgets()

    # ------------------------------------------------------------------
    # Public API
    # ------------------------------------------------------------------

    def show(self):
        """Display the full widget UI."""
        display(self.ui)
        self.update_plot()

    def step(self, mode: str, nstep: int):
        """Perform simulation step and update plot."""
        # Suppress Fortran / underlying library stdout/stderr while stepping
        with _suppress_output():
            self.sim.step(mode, temperature=self.temperature, nstep=nstep)
        self.update_plot()

    def reset(self):
        """Reset moments to initial state."""
        self.sim.reset_moments()
        self.update_plot()

    # ------------------------------------------------------------------
    # Widget construction
    # ------------------------------------------------------------------

    def _build_widgets(self):
        # --- Buttons ---
        self.btn_llg = widgets.Button(description="LLG step")
        self.btn_metropolis = widgets.Button(description="Metropolis MC")
        self.btn_heatbath = widgets.Button(description="Heatbath MC")
        self.btn_reset = widgets.Button(description="Reset moments")

        # --- Sliders ---
        self.slider_nstep = widgets.IntSlider(
            value=10, min=1, max=500, step=1, description="Steps"
        )

        self.slider_T = widgets.FloatSlider(
            value=self.sim.get_temperature(),
            min=0.0, max=1000.0, step=1.0,
            description="T (K)"
        )

        self.slider_Bz = widgets.FloatSlider(
            value=self.sim.get_field()[2],
            min=-500.0, max=500.0, step=10,
            description="Bz (T)"
        )

        # --- Status ---
        self.status = widgets.HTML()

        # --- Output area ---
        self.output = widgets.Output()

        # --- Wire callbacks ---
        self.btn_llg.on_click(lambda _: self._on_step("S"))
        self.btn_metropolis.on_click(lambda _: self._on_step("M"))
        self.btn_heatbath.on_click(lambda _: self._on_step("H"))
        self.btn_reset.on_click(lambda _: self.reset())

        self.slider_T.observe(self._on_temperature_change, names="value")
        self.slider_Bz.observe(self._on_field_change, names="value")

        # --- Layout ---
        controls = widgets.VBox([
            widgets.HBox([self.btn_llg, self.btn_metropolis, self.btn_heatbath]),
            self.slider_nstep,
            widgets.HBox([self.slider_T, self.slider_Bz]),
            self.btn_reset,
            self.status
        ])

        self.ui = widgets.HBox([controls, self.output])

    # ------------------------------------------------------------------
    # Callbacks
    # ------------------------------------------------------------------

    def _on_step(self, mode):
        nstep = self.slider_nstep.value
        self.step(mode, nstep)

    def _on_temperature_change(self, change):
        self.temperature = change["new"]
        self.sim.set_temperature(self.temperature)
        self._update_status()

    def _on_field_change(self, change):
        self.b_z = change["new"]
        B = self.sim.get_field().copy()
        B[2] = self.b_z
        self.sim.set_field(B)
        self._update_status()

    # ------------------------------------------------------------------
    # Plotting
    # ------------------------------------------------------------------

    def update_plot(self):
        with self.output:
            clear_output(wait=True)

            emom = self.sim.get_emom()[:, :, 0]   # (3, N)
            coords = self.sim.get_coords()        # (3, N)

            self._draw_quiver(coords, emom)
            self._update_status()

    def _draw_quiver(self, coords, emom):
        stride = self.config.stride
        x, y = coords[0, ::stride], coords[1, ::stride]
        u, v = emom[0, ::stride], emom[1, ::stride]

        self._ax.clear()
        self._ax.set_aspect("equal")

        self._ax.quiver(
            x, y, u, v,
            emom[2, ::stride],
            scale=self.config.scale,
            cmap=self.config.cmap
        )

        self._ax.set_title("Spin configuration")
        display(self._fig)

    # ------------------------------------------------------------------
    # Status display
    # ------------------------------------------------------------------

    def _update_status(self):
        E = self.sim.calculate_energy()
        T = self.sim.get_temperature()
        B = self.sim.get_field()

        self.status.value = f"""
        <b>Status</b><br>
        Energy: {E:.6f} mRy/atom<br>
        Temperature: {T:.2f} K<br>
        Field: ({B[0]:.2f}, {B[1]:.2f}, {B[2]:.2f}) T
        """
