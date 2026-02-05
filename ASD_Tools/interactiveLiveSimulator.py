#!/usr/bin/env python3
"""
interactiveLive.py

Interactive VTK-based UppASD visualizer using LiveSimulator.

Key bindings:
  S / s : LLG step (S = snapshot)
  M / m : Metropolis MC step (M = snapshot)
  H / h : Heat-bath MC step (H = snapshot)
  0     : Reset moments
  B/b   : +/- 1 T Bz
  N/n   : +/- 10 T Bz
  T/t   : +/- 1 K
  Y/y   : +/- 10 K
  a     : toggle atoms
  v     : toggle vectors
  c/C   : cycle colormap
  g/p   : Gouraud / PBR shading
  i/I   : static structure factor
  P     : screenshot

Author: Anders Bergman
Refactor: LiveSimulator-based
"""

import glob
import numpy as np
import vtk
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from scipy.ndimage import gaussian_filter

from vtk.util import numpy_support
from vtkmodules.vtkCommonColor import vtkColorSeries

from uppasd.run.simulator import ASDWorkspace
from uppasd.run.live_simulator import LiveSimulator


# ----------------------------------------------------------------------
# Utilities
# ----------------------------------------------------------------------

def create_lookup_table():
    """Create a two-tone lookup table used for mapping Z-colors.

    The table is built with 256 entries. The lower half interpolates from
    orange->green and the upper half from green->blue so that scalar values
    in the range [-1, 1] map to a visually distinct diverging palette.

    Returns:
        vtk.vtkLookupTable: A built lookup table with range set to [-1, 1].
    """
    lut = vtk.vtkLookupTable()
    lut.SetNumberOfTableValues(256)
    for i in range(128):
        lut.SetTableValue(i, (127 - i) / 127, i / 127, 0.0, 1.0)
    for i in range(128, 256):
        lut.SetTableValue(i, 0.0, (256 - i) / 128, (i - 128) / 128, 1.0)
    lut.SetTableRange(-1.0, 1.0)
    lut.Build()
    return lut


def cycleColorScheme(lut, colorSeries, backwards=False):
    """Cycle through the available color schemes in a `vtkColorSeries`.

    This helper advances (or reverses) the active `vtkColorSeries` index and
    updates the supplied `lut`. When the index wraps around we rebuild a
    default lookup table to ensure a deterministic mapping.

    Args:
        lut (vtk.vtkLookupTable): Existing lookup table to update.
        colorSeries (vtk.vtkColorSeries): Color series object managing schemes.
        backwards (bool): If True, step the series backwards; otherwise step
            forwards.

    Returns:
        tuple: ``(lut, colorSeries)`` where `lut` is the updated
        `vtkLookupTable` and `colorSeries` is the same object with an updated
        active scheme.

    Raises:
        TypeError: If `colorSeries` does not provide the expected methods.
    """
    cterm = -1 if backwards else 1
    colorIndex = colorSeries.GetColorScheme() + cterm
    if colorIndex <= 0:
        colorSeries.SetColorScheme(colorSeries.GetNumberOfColorSchemes() - 1)
        lut = create_lookup_table()
    elif colorIndex >= colorSeries.GetNumberOfColorSchemes():
        colorSeries.SetColorScheme(0)
        lut = create_lookup_table()
    else:
        colorSeries.SetColorScheme(colorIndex)
        colorSeries.BuildLookupTable(lut, vtkColorSeries.ORDINAL)
    return lut, colorSeries


def print_controls():
    """Display a concise help summary of keyboard controls for the viewer.

    The output is printed to stdout and intended to be displayed when the
    interactive viewer starts so users immediately see available key
    bindings and their effects. This function is side-effect free with
    respect to the visual state — it only emits text.

    Returns:
        None
    """
    controls = [
        ("S / s", "LLG step (snapshot)"),
        ("M / m", "Metropolis MC step (snapshot)"),
        ("H / h", "Heat-bath MC step (snapshot)"),
        ("0", "Reset moments to initial configuration"),
        ("B / b", "+/- 1 T on Bz"),
        ("N / n", "+/- 10 T on Bz"),
        ("T / t", "+/- 1 K temperature"),
        ("Y / y", "+/- 10 K temperature"),
        ("a", "Toggle atoms visibility"),
        ("v", "Toggle vectors visibility"),
        ("c / C", "Cycle colormap"),
        ("g / p", "Gouraud / PBR shading"),
        ("i / I", "Static structure factor (S(q))"),
        ("P", "Take screenshot"),
    ]

    print("\nInteractive Controls — keyboard bindings:")
    for key, desc in controls:
        print(f"  {key:10} - {desc}")
    print("\nTip: use lowercase for small adjustments, uppercase for snapshots/large steps.\n")


def plot_static_correlation(mz):
    """Compute and display a simple static structure factor S(q).

    The routine expects a 1D array representing the z-component of the
    magnetization on a square lattice (flattened row-major). It reshapes the
    vector to an N x N grid, computes the absolute value of the 2D FFT, and
    shows a smoothed logarithmic image of the result in a non-blocking
    matplotlib window.

    Args:
        mz (numpy.ndarray): 1D array of z-magnetization values of length N*N.

    Returns:
        None
    """
    n = int(np.sqrt(mz.size))
    cmat = mz.reshape((n, n))
    fmat = np.abs(np.fft.fft2(cmat.T))
    fmat[0, 0] = 0.0
    fmat = gaussian_filter(fmat, 1.0)
    plt.figure("S(q)")
    plt.clf()
    plt.imshow(np.log(fmat + 1.0), cmap=cm.Reds)
    plt.axis("off")
    plt.show(block=False)


# ----------------------------------------------------------------------
# Main
# ----------------------------------------------------------------------

def main():
    """Entry point for the VTK-based live visualizer.

    The `main()` function instantiates an `ASDWorkspace` and a
    `LiveSimulator`, initializes the simulator, constructs VTK data
    structures (points, vectors, glyphs), configures renderer/actors, and
    binds a keyboard-driven callback. After setup it enters the VTK
    interactor loop.

    Notes:
        - This function is intended to be run as a script; it mutates module
          level state (notably `temperature`) to maintain a simple user
          workflow during interactive sessions.

    Returns:
        None
    """

    global temperature
    # ------------------------------------------------------------
    # Live simulator
    # ------------------------------------------------------------

    workspace = ASDWorkspace(".", clean=False)
    sim = LiveSimulator(workspace)
    sim.initialize()

    natom = sim.natom
    mensemble = sim.mensemble
    temperature = sim.get_temperature()
    print("System initialized with", natom, "atoms and", mensemble, "ensembles at temperature", temperature, "K")
    # Print keyboard controls
    print_controls()

    # ------------------------------------------------------------
    # Geometry
    # ------------------------------------------------------------

    # coords = numpy_support.numpy_to_vtk(
    #     sim.get_emom()[0, :, 0].reshape(-1, 1) * 0.0,
    #     deep=True,
    #     array_type=vtk.VTK_DOUBLE,
    # )

    pos = numpy_support.numpy_to_vtk(
        np.loadtxt(glob.glob("coord.*.out")[0], usecols=(1, 2, 3)),
        deep=True,
        array_type=vtk.VTK_DOUBLE,
    )

    poly = vtk.vtkPolyData()
    pts = vtk.vtkPoints()
    pts.SetData(pos)
    poly.SetPoints(pts)

    # ------------------------------------------------------------
    # Moments → VTK arrays (DOUBLE to avoid segfaults)
    # ------------------------------------------------------------

    def update_vectors():
        """Fetch moment vectors from the simulator and write them into
        the VTK `poly` point-data arrays.

        This nested helper queries the live simulator for the current
        ensemble-average moments, converts them to double-precision NumPy
        arrays and stores both the 3-component vector and the scalar Z
        component (for lookup-table coloring) on `poly`.

        The function mutates the outer-scope `poly` object and does not
        return a value.
        """

        emom = sim.get_emom()[:, :, 0].T.astype(np.float64)

        vec = numpy_support.numpy_to_vtk(
            emom,
            deep=True,
            array_type=vtk.VTK_DOUBLE,
        )
        vec.SetNumberOfComponents(3)

        col = numpy_support.numpy_to_vtk(
            emom[:, 2],
            deep=True,
            array_type=vtk.VTK_DOUBLE,
        )

        poly.GetPointData().SetVectors(vec)
        poly.GetPointData().SetScalars(col)

    update_vectors()

    # ------------------------------------------------------------
    # Rendering
    # ------------------------------------------------------------

    ren = vtk.vtkOpenGLRenderer()
    ren.SetBackground(1, 1, 1)

    win = vtk.vtkRenderWindow()
    win.AddRenderer(ren)
    win.SetSize(1000, 700)

    iren = vtk.vtkRenderWindowInteractor()
    iren.SetRenderWindow(win)
    iren.SetInteractorStyle(vtk.vtkInteractorStyleTrackballCamera())

    # Atoms
    sphere = vtk.vtkSphereSource()
    sphere.SetRadius(0.6)

    lut = create_lookup_table()
    colorSeries = vtk.vtkColorSeries()
    # Use BREWER_QUALITATIVE_PASTEL1 as the default color scheme
    try:
        colorSeries.SetColorScheme(vtkColorSeries.BREWER_DIVERGING_SPECTRAL_11)
    except Exception:
        # Fallback: keep whatever default the VTK installation provides
        pass
    colorSeries.BuildLookupTable(lut, vtkColorSeries.ORDINAL)

    glyph_atoms = vtk.vtkGlyph3DMapper()
    glyph_atoms.SetInputData(poly)
    glyph_atoms.SetSourceConnection(sphere.GetOutputPort())
    glyph_atoms.SetScaleModeToNoDataScaling()
    glyph_atoms.SetLookupTable(lut)

    atom_actor = vtk.vtkActor()
    atom_actor.SetMapper(glyph_atoms)
    # Set atom actor material / shading defaults similar to interactiveASD
    try:
        p = atom_actor.GetProperty()
        p.SetSpecular(0.3)
        p.SetSpecularPower(60)
        p.SetAmbient(0.2)
        p.SetDiffuse(0.8)
        p.SetRoughness(0.6)
        p.SetMetallic(0.7)
        try:
            p.SetInterpolationToGouraud()
        except Exception:
            pass
        p.ShadingOn()
    except Exception:
        pass

    # Vectors
    arrow = vtk.vtkArrowSource()
    glyph_vec = vtk.vtkGlyph3DMapper()
    glyph_vec.SetInputData(poly)
    glyph_vec.SetSourceConnection(arrow.GetOutputPort())
    glyph_vec.OrientOn()
    glyph_vec.SetScaleFactor(2.0)
    glyph_vec.SetLookupTable(lut)

    vec_actor = vtk.vtkActor()
    vec_actor.SetMapper(glyph_vec)
    # Vector actor material defaults
    try:
        pv = vec_actor.GetProperty()
        pv.SetSpecular(0.0)
        pv.SetSpecularPower(10)
        pv.SetAmbient(0.3)
        pv.SetDiffuse(0.5)
        pv.SetOpacity(1.0)
        pv.SetRoughness(0.7)
        try:
            pv.SetInterpolationToGouraud()
        except Exception:
            pass
        pv.ShadingOn()
    except Exception:
        pass

    ren.AddActor(atom_actor)
    ren.AddActor(vec_actor)

    # ------------------------------------------------------------
    # Text
    # ------------------------------------------------------------

    def make_text(y):
        """Create and register a small `vtkTextActor` at the given Y position.

        Args:
            y (int): Vertical pixel position for the text actor's display
                origin. Higher values are further up the render window.

        Returns:
            vtk.vtkTextActor: The created and added text actor instance.
        """

        t = vtk.vtkTextActor()
        t.GetTextProperty().SetFontSize(28)
        t.GetTextProperty().SetColor(0, 0, 0)
        t.SetDisplayPosition(20, y)
        ren.AddActor(t)
        return t

    txt_T = make_text(620)
    txt_B = make_text(580)
    txt_E = make_text(40)

    # ------------------------------------------------------------
    # Key handler
    # ------------------------------------------------------------

    def on_keypress(obj, evt):
        """Handle keyboard events for the interactor and mutate visual state.

        This callback interprets a predefined set of key bindings to drive
        simulation steps, toggle visibility, adjust physical parameters like
        temperature and field, cycle color schemes, and trigger diagnostics
        (e.g. static structure factor). It updates displayed text actors and
        requests a re-render of the window when appropriate.

        Args:
            obj: The VTK interactor object emitting the event (provides
                `GetKeySym()` to query the pressed key).
            evt: Event name (ignored by this handler).

        Side effects:
            - Calls methods on the `sim` LiveSimulator instance to advance
              the simulation or change parameters.
            - Changes actor properties and lookup tables.
            - Updates overlay text actors and re-renders `win`.

        Returns:
            None
        """

        global temperature
        nonlocal lut, colorSeries
        key = obj.GetKeySym()

        if key in ("S", "s"):
            sim.step(mode="S", nstep=10, temperature=temperature)
        elif key in ("M", "m"):
            sim.step(mode="M", nstep=10, temperature=temperature)
        elif key in ("H", "h"):
            sim.step(mode="H", nstep=10, temperature=temperature + 1.0e-6)
        elif key == "0":
            sim.reset_moments()
        elif key == "B":
            B = sim.get_field()
            B[2] += 1.0
            sim.set_field(B)
        elif key == "b":
            B = sim.get_field()
            B[2] -= 1.0
            sim.set_field(B)
        elif key == "T":
            temperature += 1.0
        elif key == "t":
            temperature = max(temperature - 1.0, 1e-6)
        elif key == "a":
            atom_actor.SetVisibility(not atom_actor.GetVisibility())
        elif key == "v":
            vec_actor.SetVisibility(not vec_actor.GetVisibility())
        elif key in ("c", "C"):
            # Cycle colormap
            lut, colorSeries = cycleColorScheme(lut, colorSeries, backwards=(key == "c"))
            glyph_atoms.SetLookupTable(lut)
            glyph_vec.SetLookupTable(lut)
            glyph_atoms.Modified()
            glyph_vec.Modified()
        elif key in ("g", "G"):
            # Gouraud shading
            try:
                atom_actor.GetProperty().SetInterpolationToGouraud()
                vec_actor.GetProperty().SetInterpolationToGouraud()
            except Exception:
                pass
        elif key == "p":
            # PBR shading
            try:
                atom_actor.GetProperty().SetInterpolationToPBR()
                vec_actor.GetProperty().SetInterpolationToPBR()
            except Exception:
                pass
        elif key in ("i", "I"):
            mz = sim.get_emom()[2, :, 0]
            plot_static_correlation(mz)

        update_vectors()

        txt_T.SetInput(f"T = {temperature:.2f} K")
        txt_B.SetInput(f"Bz = {sim.get_field()[2]:.2f} T")
        txt_E.SetInput(f"E = {sim.calculate_energy():.6f} mRy/atom")

        win.Render()

    iren.AddObserver("KeyPressEvent", on_keypress)

    # ------------------------------------------------------------
    # Start
    # ------------------------------------------------------------

    win.Render()
    iren.Initialize()
    iren.Start()


if __name__ == "__main__":
    main()
