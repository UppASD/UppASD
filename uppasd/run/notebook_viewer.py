"""
NotebookLiveViewer
==================

Notebook-oriented interactive viewer for UppASD LiveSimulator.

DESIGN CONTRACT (NON-NEGOTIABLE)
--------------------------------
- Viewer NEVER touches pyasd directly
- Viewer NEVER owns or mutates physics state
- Viewer NEVER stores magnetic moments
- All stepping goes through LiveSimulator.step(...)
- All visualization pulls fresh data from LiveSimulator
- Reset is simulator-owned, not viewer-owned

This file is intentionally conservative.
If something seems "missing", it is probably forbidden by design.

Author: Anders Bergman + ChatGPT
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Optional

import numpy as np
import matplotlib.pyplot as plt
import ipywidgets as widgets
from IPython.display import display, clear_output

from .live_simulator import LiveSimulator


# ======================================================================
# Step control (pure UI intent, no physics)
# ======================================================================

@dataclass
class StepControl:
    """
    Control parameters for simulation stepping.

    This object represents *intent*, not state.
    """
    mode: str = "S"   # 'S' = LLG, 'M' = Metropolis, 'H' = Heatbath
    nstep: int = 10


# ======================================================================
# Viewer configuration (pure visualization)
# ======================================================================

@dataclass
class ViewerConfig:
    """
    Visualization-only configuration.

    MUST NOT contain:
    - temperature
    - field
    - moments
    - any physics-related parameter
    """
    scale: float = 20.0
    stride: int = 1
    cmap: str = "coolwarm"
    figsize: tuple = (6, 6)


# ======================================================================
# NotebookLiveViewer
# ======================================================================

class NotebookLiveViewer:
    """
    Jupyter-based live viewer for UppASD LiveSimulator.

    Responsibilities:
    -----------------
    - UI controls (widgets)
    - Visualization
    - Delegating stepping to LiveSimulator

    Explicitly NOT responsible for:
    --------------------------------
    - Physics state
    - Moment storage
    - Temperature persistence inside UppASD
    """

    # ------------------------------------------------------------------
    def __init__(
        self,
        sim: LiveSimulator,
        config: Optional[ViewerConfig] = None,
    ):
        """
        Parameters
        ----------
        sim : LiveSimulator
            Initialized LiveSimulator instance.
        config : ViewerConfig, optional
            Visualization configuration.
        """
        if not getattr(sim, "initialized", False):
            raise RuntimeError(
                "LiveSimulator must be initialized before creating viewer"
            )

        self.sim = sim
        self.config = config or ViewerConfig()
        self.control = StepControl()

        # Local UI-side parameters (NOT physics state)
        self.temperature = self.sim.get_temperature()
        self.field = self.sim.get_field().copy()

        # Persistent matplotlib figure (important for notebooks)
        plt.ioff()
        self._fig, self._ax = plt.subplots(figsize=self.config.figsize)
        plt.ion()

        self._build_widgets()

    # ------------------------------------------------------------------
    # Public API
    # ------------------------------------------------------------------

    def show(self):
        """Display widget UI and initial visualization."""
        display(self.ui)
        self.render()

    def step(self):
        """
        Perform a simulation step via LiveSimulator.

        IMPORTANT:
        - No moment storage
        - No return value is cached
        """
        self.sim.step(
            mode=self.control.mode,
            nstep=self.control.nstep,
            temperature=self.temperature,
        )
        self.render()

    def reset(self):
        """
        Reset magnetic moments via LiveSimulator.

        Viewer does NOT attempt to restore moments itself.
        """
        self.sim.reset()
        self.render()

    # ------------------------------------------------------------------
    # UI construction
    # ------------------------------------------------------------------

    def _build_widgets(self):
        # --- Mode selection ---
        self.mode_select = widgets.Dropdown(
            options=[
                ("LLG", "S"),
                ("Metropolis MC", "M"),
                ("Heatbath MC", "H"),
            ],
            value="S",
            description="Mode:",
        )

        # --- Step count ---
        self.step_input = widgets.IntText(
            value=10,
            description="Steps:",
            layout=widgets.Layout(width="140px"),
        )

        # --- Buttons ---
        self.step_button = widgets.Button(
            description="Step",
            button_style="primary",
        )

        self.reset_button = widgets.Button(
            description="Reset",
            button_style="warning",
        )

        # --- Temperature ---
        self.slider_T = widgets.FloatSlider(
            value=self.temperature,
            min=0.0,
            max=1000.0,
            step=1.0,
            description="T (K)",
        )

        # --- Field ---
        self.slider_Bz = widgets.FloatSlider(
            value=self.field[2],
            min=-500.0,
            max=500.0,
            step=10.0,
            description="Bz (T)",
        )

        # --- Output ---
        self.output = widgets.Output()
        self.status = widgets.HTML()

        # --- Callbacks ---
        self.step_button.on_click(self._on_step)
        self.reset_button.on_click(self._on_reset)
        self.slider_T.observe(self._on_temperature_change, names="value")
        self.slider_Bz.observe(self._on_field_change, names="value")

        # --- Layout ---
        controls = widgets.VBox([
            widgets.HBox([
                self.mode_select,
                self.step_input,
                self.step_button,
                self.reset_button,
            ]),
            self.slider_T,
            self.slider_Bz,
            self.status,
        ])

        self.ui = widgets.HBox([controls, self.output])

    # ------------------------------------------------------------------
    # Callbacks
    # ------------------------------------------------------------------

    def _on_step(self, _):
        self.control.mode = self.mode_select.value
        self.control.nstep = self.step_input.value
        self.step()

    def _on_reset(self, _):
        self.reset()

    def _on_temperature_change(self, change):
        # Viewer-local value only; passed explicitly during step()
        self.temperature = change["new"]
        self._update_status()

    def _on_field_change(self, change):
        self.field[2] = change["new"]
        self.sim.set_field(self.field)
        self._update_status()

    # ------------------------------------------------------------------
    # Rendering
    # ------------------------------------------------------------------

    def render(self):
        """Render current simulation state."""
        with self.output:
            clear_output(wait=True)
            self._plot_spins()
            self._update_status()

    # ------------------------------------------------------------------
    # Visualization
    # ------------------------------------------------------------------

    def _plot_spins(self):
        """
        Geometry-aware quiver plot.

        - Arrow direction: (m_x, m_y)
        - Arrow color: m_z
        - Coordinates from real-space geometry
        """
        emom = self.sim.get_emom()[:, :, 0].copy()
        coords = self.sim.get_coords().copy()

        stride = self.config.stride

        x = coords[0, ::stride]
        y = coords[1, ::stride]
        u = emom[0, ::stride]
        v = emom[1, ::stride]
        mz = emom[2, ::stride]

        self._ax.clear()
        self._ax.set_aspect("equal")

        self._ax.quiver(
            x, y, u, v, mz,
            scale=self.config.scale,
            cmap=self.config.cmap,
        )

        self._ax.set_title("Spin configuration")
        self._ax.set_xlabel("x")
        self._ax.set_ylabel("y")

        display(self._fig)

    # ------------------------------------------------------------------
    # Status output
    # ------------------------------------------------------------------

    def _update_status(self):
        energy = self.sim.calculate_energy()
        T = self.temperature
        B = self.sim.get_field()

        self.status.value = f"""
        <b>Status</b><br>
        Energy: {energy:.6f} mRy/atom<br>
        Temperature: {T:.2f} K<br>
        Field: ({B[0]:.2f}, {B[1]:.2f}, {B[2]:.2f}) T
        """
