"""Jupyter helper for interactive LiveSimulator sessions.

This module provides `JupyterLiveSimulator`, a thin adapter that wraps
`LiveSimulator` with Jupyter-friendly widgets and plotting helpers.

The intent is to offer a quick interactive frontend inside notebooks
for stepping, resetting and visualizing spin configurations without
exposing or duplicating simulator internals.

Usage example:

    from uppasd.run.jupyter_simulator import JupyterLiveSimulator
    jl = JupyterLiveSimulator(workspace)
    jl.show()

"""

import numpy as np
import matplotlib.pyplot as plt
from IPython.display import display, clear_output
import ipywidgets as widgets
from uppasd.run.live_simulator import LiveSimulator


class JupyterLiveSimulator:
    """Notebook-oriented wrapper around ``LiveSimulator``.

    The class exposes a small set of widget controls (mode, temperature,
    steps, Bz) and plotting utilities to visualize the current spin
    configuration via a matplotlib quiver plot. It delegates all
    physics operations (initialization, stepping, resetting) to
    ``LiveSimulator`` and intentionally does not store moment data.

    Args:
        workspace: Path-like or Workspace object accepted by
            ``LiveSimulator`` to locate simulation files and output.
    """

    def __init__(self, workspace):
        self.sim = LiveSimulator(workspace)
        try:
            self.sim.initialize()
        except (RuntimeError, ValueError, TypeError) as e:
            print(f"Initialization Error: {e}")
            return

        self.coords = self.sim.get_coords()
        self.x, self.y = self.coords[0, :], self.coords[1, :]

        self.output = widgets.Output()
        self.setup_widgets()

    def setup_widgets(self):
        """Create and configure ipywidgets used by the frontend.

        This method constructs dropdowns, sliders, buttons and output
        areas, wires callback handlers and performs an initial status
        update. It does not perform any heavy computation; that is
        done when the user triggers actions via the UI.
        """
        # 1. Inputs
        self.mode_sel = widgets.Dropdown(
            options=[
                ("LLG (S)", "S"),
                ("Metropolis MC (M)", "M"),
                ("Heat-bath MC (H)", "H"),
            ],
            value="S",
            description="Mode:",
        )
        self.temp_slider = widgets.FloatSlider(
            value=self.sim.get_temperature(),
            min=0.0001,
            max=1500.0,
            step=1.0,
            description="Temp (K):",
            readout_format=" 5.1f",
        )
        self.steps_input = widgets.IntText(value=100, description="N Steps:")
        self.bz_input = widgets.FloatText(
            value=self.sim.get_field()[2], description="Bz (T):"
        )

        # 2. Status & Log Displays
        self.status = widgets.HTML(
            value="<b>Status</b>: Ready", layout={"height": "100px", "width": "90%"}
        )
        self.log_area = widgets.Textarea(
            value="",
            placeholder="Simulation logs will appear here...",
            description="",
            layout={"height": "100px", "width": "90%"},
            # description='Log:', layout={'height': '100px', 'width': '90%'}
        )

        # 3. Buttons
        self.relax_btn = widgets.Button(
            description="Relax", button_style="success", icon="play"
        )
        self.reset_btn = widgets.Button(
            description="Reset", button_style="warning", icon="refresh"
        )

        self.relax_btn.on_click(self.on_relax_clicked)
        self.reset_btn.on_click(self.on_reset_clicked)

        # 4. Live Observers
        self.temp_slider.observe(self._on_control_change, names="value")
        self.mode_sel.observe(self._on_control_change, names="value")
        self.steps_input.observe(self._on_control_change, names="value")

        # 5. Layout
        controls = widgets.VBox(
            [
                self.mode_sel,
                self.temp_slider,
                self.steps_input,
                self.bz_input,
                self.status,
                widgets.HBox([self.relax_btn, self.reset_btn]),
                self.log_area,
            ]
        )
        self.ui = widgets.HBox([controls, self.output])
        self._update_status_display()
        self.on_reset_clicked(True)

    def _update_status_display(self, msg=None):
        """Refresh the HTML status widget with current parameters.

        The displayed information includes the selected integration mode,
        temperature, number of steps and Bz. If `msg` is provided it is
        shown as the "Last Action" field.

        Args:
            msg (str, optional): Short message describing the last action.
        """
        mode_map = {"S": "LLG", "M": "Metropolis", "H": "Heat-bath"}
        mode = mode_map.get(self.mode_sel.value, self.mode_sel.value)
        # Read Bz from the field input to reflect what will be used in runs
        try:
            bz_val = float(self.bz_input.value)
        except (ValueError, TypeError):
            # Fallback to reading from simulator field
            try:
                bz_val = float(self.sim.get_field()[2])
            except (ValueError, TypeError, IndexError):
                bz_val = 0.0

        self.status.value = (
            f"<div style='border: 1px solid #ccc; padding: 5px; "
            f"background-color: #f9f9f9;'>"
            f"<b>Current Config</b>: {mode} <br> "
            f"T= {self.temp_slider.value:.2f} K | "
            f"{self.steps_input.value} steps | "
            f"Bz={bz_val:.3f} T<br>"
            f"<b>Last Action</b>: {msg if msg else 'None'}</div>"
        )

    def _on_control_change(self, change):
        """Callback invoked when a control widget value changes.

        The purpose of this small hook is to refresh the status display
        immediately so the user receives instant feedback while
        adjusting parameters.

        Args:
            change: The ipywidgets change event dictionary (ignored).
        """
        self._update_status_display("Adjusting parameters...")

    def update_plot(self):
        """Render the current spin configuration to the notebook output.

        The method queries `LiveSimulator.get_emom()` for the current
        magnetic moments and draws a quiver plot using matplotlib. Any
        plotting errors are captured and written to the log area.
        """
        with self.output:
            clear_output(wait=True)
            try:
                emom = self.sim.get_emom()[:, :, 0]
                fig, ax = plt.subplots(figsize=(7, 5))
                q = ax.quiver(
                    self.x,
                    self.y,
                    emom[0, :],
                    emom[1, :],
                    emom[2, :],
                    cmap="RdYlBu",
                    pivot="mid",
                    scale=30,
                )
                ax.set_aspect("equal")
                plt.colorbar(q, label="$m_z$")
                display(fig)
                plt.close(fig)
            except (ValueError, IndexError, KeyError, RuntimeError) as e:
                self.log_area.value += f"Plotting Error: {e}\n"

    def on_relax_clicked(self, b):
        """Handler for the "Relax" button click.

        Reads the current widget values (mode, temperature, steps, Bz),
        updates the simulator field and runs ``LiveSimulator.step``. The
        UI is updated with the resulting plot and a brief log entry.

        Args:
            b: Button click event payload (ignored).
        """
        with self.output:
            try:
                # Local capture of values to avoid 'stale' widget reads
                _t = np.float64(self.temp_slider.value)
                _n = int(self.steps_input.value)
                _m = str(self.mode_sel.value)
                _bz = np.float64(self.bz_input.value)

                self.log_area.value += f"Starting relaxation:\n {_m} at {_t}K...\n"

                # Update Field
                new_field = np.array([0.0, 0.0, _bz], dtype=np.float64)
                self.sim.set_field(new_field)

                # Test to set the temperature again:
                # self.sim.set_temperature(_t)

                # Step (Avoid calling sim.set_temperature to prevent kernel crash)
                self.sim.step(mode=_m, temperature=_t, nstep=_n)

                self.update_plot()
                self._update_status_display(f"Finished {_n} steps of {_m}")
                self.log_area.value += "Relaxation complete.\n"

            except (ValueError, TypeError, RuntimeError, KeyError, IndexError) as e:
                self.log_area.value += f"Execution Error: {e}\n"

    def on_reset_clicked(self, b=None):
        """Handler for the "Reset" button click.

        Calls ``LiveSimulator.reset_moments()``, refreshes the plot and
        updates the status display. Any runtime errors from the
        simulator are appended to the widget log area.

        Args:
            b: Optional button event payload (ignored).
        """
        with self.output:
            try:
                self.sim.reset_moments()
                self.update_plot()
                self._update_status_display("Moments reset to initial state")
            except (ValueError, IndexError, KeyError, RuntimeError) as e:
                self.log_area.value += f"Reset Error: {e}\n"

    def show(self):
        """Display the UI and initial plot inside the notebook.

        This is a convenience wrapper that displays the widget layout
        and triggers an initial plot update.
        """
        display(self.ui)
        self.update_plot()
