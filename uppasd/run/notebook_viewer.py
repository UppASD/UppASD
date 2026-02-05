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
    scale: float = 30.0
    stride: int = 1
    cmap: str = "coolwarm"
    figsize: tuple = (7, 5)


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
        # Cache atom coordinates and fix axis limits so the plot area remains
        # stable across updates (prevents shrinking when adding colorbar).
        coords = self.sim.get_coords().copy()
        x = coords[0, :]
        y = coords[1, :]
        xpad = (x.max() - x.min()) * 0.05 if x.max() > x.min() else 1.0
        ypad = (y.max() - y.min()) * 0.05 if y.max() > y.min() else 1.0
        self._xlim = (x.min() - xpad, x.max() + xpad)
        self._ylim = (y.min() - ypad, y.max() + ypad)
        self._ax.set_xlim(*self._xlim)
        self._ax.set_ylim(*self._ylim)

        # Colorbar axis, colorbar and quiver handles (created when plotting)
        self._cax = None
        self._cbar = None
        self._quiver = None

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
        try:
            self.sim.step(
                mode=self.control.mode,
                nstep=self.control.nstep,
                temperature=self.temperature,
            )
        except Exception as e:
            # Log and re-render status
            try:
                self.log_area.value += f"Step Error: {e}\n"
            except Exception:
                pass
        finally:
            self.render()

    def reset(self):
        """
        Reset magnetic moments via LiveSimulator.

        Viewer does NOT attempt to restore moments itself.
        """
        try:
            # LiveSimulator exposes `reset_moments()`; call that explicitly.
            self.sim.reset_moments()
            self.log_area.value += "Moments reset to initial state.\n"
        except Exception as e:
            try:
                self.log_area.value += f"Reset Error: {e}\n"
            except Exception:
                pass
        finally:
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
            layout=widgets.Layout(width='95%'),
            style={'description_width': '110px'},
        )

        # --- Step count ---
        self.step_input = widgets.IntSlider(
            value=10,
            min=1,
            max=500,
            step=1,
            description="Steps:",
            layout=widgets.Layout(width='95%'),
            style={'description_width': '110px'},
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
            layout=widgets.Layout(width='95%'),
            style={'description_width': '110px'},
        )

        # --- Field ---
        self.slider_Bz = widgets.FloatSlider(
            value=self.field[2],
            min=-500.0,
            max=500.0,
            step=10.0,
            description="Bz (T)",
            layout=widgets.Layout(width='95%'),
            style={'description_width': '110px'},
        )

        # --- Output ---
        self.output = widgets.Output(layout=widgets.Layout(width='65%'))
        self.status = widgets.HTML(layout=widgets.Layout(width='95%'))
        self.log_area = widgets.Textarea(
            value='',
            placeholder='Simulation logs will appear here...',
            description='',
            layout=widgets.Layout(height='120px', width='95%')
        )

        # --- Callbacks ---
        self.step_button.on_click(self._on_step)
        self.reset_button.on_click(self._on_reset)
        self.slider_T.observe(self._on_temperature_change, names="value")
        self.slider_Bz.observe(self._on_field_change, names="value")

        # --- Layout ---
        # Left controls column (stacked) and right plot column
        # Arrange controls: Mode, Temp, Field, Steps, Status, [buttons], Log
        controls = widgets.VBox([
            self.mode_select,
            self.slider_T,
            self.slider_Bz,
            self.step_input,
            self.status,
            widgets.HBox([self.step_button, self.reset_button], layout=widgets.Layout(justify_content='flex-start', width='95%')),
            self.log_area,
        ], layout=widgets.Layout(width='35%', align_items='stretch'))

        self.ui = widgets.HBox([controls, self.output], layout=widgets.Layout(width='100%'))

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

        # Keep axis limits stable and update aspect.
        self._ax.set_aspect("equal")
        self._ax.set_xlim(*getattr(self, '_xlim', (None, None)))
        self._ax.set_ylim(*getattr(self, '_ylim', (None, None)))

        # If we haven't created the quiver yet, create it and a colorbar.
        # Compute a dynamic scale so arrows are large but fit the axes. The
        # scale is computed so that the median arrow length is roughly
        # `arrow_frac` of the minimal axis extent.
        try:
            min_extent = min(self._xlim[1] - self._xlim[0], self._ylim[1] - self._ylim[0])
        except Exception:
            min_extent = None
        meanlen = np.mean(np.hypot(u, v)) if u.size > 0 else 0.0
        arrow_frac = 0.04
        if min_extent and meanlen > 0.0:
            scale_val = float(meanlen / (arrow_frac * min_extent))
            if not np.isfinite(scale_val) or scale_val <= 0.0:
                scale_val = float(self.config.scale)
        else:
            scale_val = float(self.config.scale)

        if self._quiver is None:
            self._quiver = self._ax.quiver(
                x, y, u, v, mz,
                scale=scale_val,
                scale_units='xy',
                cmap=self.config.cmap,
                pivot='mid',
            )
            # Create a dedicated colorbar axis matching the height of the main ax
            try:
                axpos = self._ax.get_position()
                pad = 0.01
                cbar_width = 0.03
                if self._cax is not None:
                    try:
                        self._cax.remove()
                    except Exception:
                        pass
                self._cax = self._fig.add_axes((axpos.x1 + pad, axpos.y0, cbar_width, axpos.height))
                self._cbar = self._fig.colorbar(self._quiver, cax=self._cax, label='$m_z$')
            except Exception:
                self._cax = None
                self._cbar = None
        else:
            # Update vector components and color array in-place to avoid
            # changing axes/layout and to keep the colorbar intact.
            try:
                # Update quiver data in-place; preserve colorbar axis position.
                self._quiver.set_UVC(u, v, mz)
                self._quiver.set_array(mz)
                # Update colorbar mappable and reposition cax to match ax height
                if self._cbar is not None and hasattr(self._cbar, 'mappable'):
                    try:
                        self._cbar.mappable.set_array(mz)
                    except Exception:
                        pass
                if self._cax is not None:
                    try:
                        axpos = self._ax.get_position()
                        pad = 0.01
                        cbar_width = 0.03
                        self._cax.set_position((axpos.x1 + pad, axpos.y0, cbar_width, axpos.height))
                    except Exception:
                        pass
            except Exception:
                # Fallback: recreate quiver if in-place update fails
                try:
                    if self._quiver is not None:
                        self._quiver.remove()
                except Exception:
                    pass
                self._quiver = self._ax.quiver(x, y, u, v, mz, scale=self.config.scale, cmap=self.config.cmap, pivot='mid')
                try:
                    if self._cbar is not None:
                        try:
                            self._cbar.remove()
                        except Exception:
                            pass
                    if self._cax is not None:
                        try:
                            self._cax.remove()
                        except Exception:
                            pass
                    axpos = self._ax.get_position()
                    pad = 0.01
                    cbar_width = 0.03
                    self._cax = self._fig.add_axes((axpos.x1 + pad, axpos.y0, cbar_width, axpos.height))
                    self._cbar = self._fig.colorbar(self._quiver, cax=self._cax, label='$m_z$')
                except Exception:
                    self._cbar = None

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
