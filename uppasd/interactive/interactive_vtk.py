#!/usr/bin/env python3
"""
Interactive VTK visualization for UppASD (new Python API).

Design:
- No mutable Fortran state access
- File-based results (ASDResults)
- Explicit runtime controls
- Preserves original VTK interaction model
"""

import numpy as np
import vtk
from vtk.util import numpy_support
from pathlib import Path

# UppASD imports (NEW API)
from uppasd.core.system import SpinSystem
from uppasd.core.results import ASDResults
from uppasd.input.inputdata import ASDInput
from uppasd.run.simulator import ASDWorkspace, UppASDSimulator

# ----------------------------------------------------------------------
# USER SETUP (adapt as needed)
# ----------------------------------------------------------------------

WORKDIR = Path("vtk_run")
CUTOFF_RENDER = 1e6  # no filtering

# Example placeholders — replace with real system construction
def build_system():
    """
    User-provided system construction.
    Replace with POSCAR/ASE/etc if desired.
    """
    # Minimal dummy example
    cell = np.eye(3) * 10.0
    positions = np.array([[0, 0, 0]])
    species = [1]
    moments = np.array([[0, 0, 1.0]])

    return SpinSystem(
        cell=cell,
        positions=positions,
        species=species,
        moments=moments,
    )


def build_input():
    inp = ASDInput()
    inp.block("simulation").set(
        nstep=10,
        temperature=0.0,
    )
    return inp


# ----------------------------------------------------------------------
# VTK HELPERS
# ----------------------------------------------------------------------

def numpy_to_vtk_vectors(vectors):
    return numpy_support.numpy_to_vtk(vectors.astype(np.float32), deep=False)


def numpy_to_vtk_scalars(scalars):
    return numpy_support.numpy_to_vtk(scalars.astype(np.float32), deep=False)


# ----------------------------------------------------------------------
# MAIN APPLICATION
# ----------------------------------------------------------------------

def main():
    # ------------------------------------------------------------
    # Setup UppASD
    # ------------------------------------------------------------
    system = build_system()
    inp = build_input()

    ws = ASDWorkspace(WORKDIR, clean=True)
    ws.prepare(system=system, inp=inp)

    sim = UppASDSimulator(ws)
    sim.initialize()

    runtime = {
        "temperature": 0.0,
    }

    # ------------------------------------------------------------
    # Initial results
    # ------------------------------------------------------------
    results = ASDResults(ws.path)

    # ------------------------------------------------------------
    # VTK setup
    # ------------------------------------------------------------
    renWin = vtk.vtkRenderWindow()
    renWin.SetSize(1000, 700)
    ren = vtk.vtkRenderer()
    ren.SetBackground(1.0, 1.0, 1.0)
    renWin.AddRenderer(ren)

    # Geometry container
    polydata = vtk.vtkPolyData()

    # Initial geometry
    pts = vtk.vtkPoints()
    pts.SetData(numpy_support.numpy_to_vtk(system.positions, deep=False))
    polydata.SetPoints(pts)

    # Dummy initial vectors
    vecs = np.zeros_like(system.positions)
    scal = np.zeros(len(system.positions))

    polydata.GetPointData().SetVectors(numpy_to_vtk_vectors(vecs))
    polydata.GetPointData().SetScalars(numpy_to_vtk_scalars(scal))

    # ------------------------------------------------------------
    # Glyphs (spins)
    # ------------------------------------------------------------
    arrow = vtk.vtkArrowSource()
    arrow.SetTipRadius(0.25)
    arrow.SetShaftRadius(0.15)

    glyph = vtk.vtkGlyph3DMapper()
    glyph.SetInputData(polydata)
    glyph.SetSourceConnection(arrow.GetOutputPort())
    glyph.OrientOn()
    glyph.SetScaleModeToNoDataScaling()
    glyph.SetScaleFactor(2.0)
    glyph.ScalarVisibilityOn()

    actor = vtk.vtkActor()
    actor.SetMapper(glyph)
    ren.AddActor(actor)

    # ------------------------------------------------------------
    # Text overlay
    # ------------------------------------------------------------
    txt = vtk.vtkTextActor()
    txt.SetInput("Press M: relax | Q: quit")
    txt.GetTextProperty().SetFontSize(28)
    txt.GetTextProperty().SetColor(0, 0, 0)
    txt.SetDisplayPosition(20, 20)
    ren.AddActor(txt)

    # ------------------------------------------------------------
    # Interaction logic
    # ------------------------------------------------------------
    def update_from_results():
        nonlocal results
        results.reload()

        mag = results.magnetization
        if mag is None:
            return

        vecs = np.column_stack([mag["mx"], mag["my"], mag["mz"]])
        scal = mag["mz"]

        polydata.GetPointData().SetVectors(numpy_to_vtk_vectors(vecs))
        polydata.GetPointData().SetScalars(numpy_to_vtk_scalars(scal))
        polydata.Modified()

    def on_keypress(obj, event):
        key = obj.GetKeySym()

        if key.lower() == "q":
            sim.finalize()
            renWin.Finalize()
            obj.TerminateApp()
            return

        if key.lower() == "m":
            sim.set_runtime_controls(**runtime)
            sim.relax()
            update_from_results()
            renWin.Render()

        if key == "T":
            runtime["temperature"] += 1.0
            print("T =", runtime["temperature"])

        if key == "t":
            runtime["temperature"] = max(runtime["temperature"] - 1.0, 0.0)
            print("T =", runtime["temperature"])

    # ------------------------------------------------------------
    # Start interaction
    # ------------------------------------------------------------
    iren = vtk.vtkRenderWindowInteractor()
    iren.SetRenderWindow(renWin)
    iren.AddObserver("KeyPressEvent", on_keypress)

    renWin.Render()
    iren.Initialize()
    iren.Start()


if __name__ == "__main__":
    main()
