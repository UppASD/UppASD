#!/usr/bin/env python3
"""
interactiveASD_live.py

VTK-based interactive UppASD visualization using LiveSimulator.

Keys
----
S / s : LLG relaxation
M / m : Metropolis MC
H / h : Heat-bath MC
0     : Reset spins
P     : Screenshot
B/b   : ±1 T field (z)
N/n   : ±10 T field (z)
T/t   : ±1 K
Y/y   : ±10 K
c/C   : Cycle colormap
g/p   : Gouraud / PBR shading
a     : Toggle atoms
v     : Toggle vectors
i/I   : Static correlation S(q)
"""

import glob
import numpy as np
import vtk
import matplotlib.pyplot as plt
import matplotlib.cm as cm

from scipy.ndimage import gaussian_filter
from vtk.util import numpy_support
from vtkmodules.vtkCommonColor import vtkColorSeries
from vtk.util.numpy_support import numpy_to_vtk

from uppasd.run.live_simulator import LiveSimulator
from uppasd.run.simulator import ASDWorkspace




# ============================================================
# Utilities
# ============================================================

number_of_screenshots = 1


def screenshot(renwin):
    global number_of_screenshots
    w2i = vtk.vtkWindowToImageFilter()
    w2i.SetInput(renwin)
    w2i.ReadFrontBufferOff()
    w2i.Update()

    writer = vtk.vtkPNGWriter()
    writer.SetFileName(f"snap{number_of_screenshots:05d}.png")
    writer.SetInputConnection(w2i.GetOutputPort())
    writer.Write()

    number_of_screenshots += 1
    print("Screenshot saved")


def create_lookup_table():
    lut = vtk.vtkLookupTable()
    for i in range(128):
        lut.SetTableValue(i, (127 - i) / 127, i / 127, 0, 1)
    for i in range(128, 256):
        lut.SetTableValue(i, 0, (256 - i) / 128, (i - 128) / 128, 1)
    lut.SetTableRange(-1.0, 1.0)
    lut.Build()
    return lut


def plot_static_correlation(mat, do_log=False):
    if do_log:
        mat = np.log(mat + 1.0)
    else:
        mat = np.abs(mat)

    plt.figure()
    plt.imshow(mat, cmap=cm.Reds)
    plt.title("Static correlation S(q)")
    plt.axis("off")
    plt.show()


# ============================================================
# Main
# ============================================================

def main():
    # --------------------------------------------------------
    # Live simulator
    # --------------------------------------------------------
    workspace = ASDWorkspace(".", clean=False)
    sim = LiveSimulator(workspace)
    sim.initialize()

    natom = sim.natom
    mensemble = sim.mensemble

    # --------------------------------------------------------
    # VTK scene
    # --------------------------------------------------------
    ren = vtk.vtkRenderer()
    renwin = vtk.vtkRenderWindow()
    renwin.AddRenderer(ren)
    renwin.SetSize(1100, 700)
    ren.SetBackground(1, 1, 1)

    # --------------------------------------------------------
    # Geometry
    # --------------------------------------------------------
    coords = sim.workspace.path.glob("coord.*.out")
    coords = sorted(coords)
    if not coords:
        raise RuntimeError("No coord.*.out found")

    coords = np.loadtxt(coords[0])[:, 1:4]
    points = vtk.vtkPoints()
    points.SetData(numpy_support.numpy_to_vtk(coords))

    poly = vtk.vtkPolyData()
    poly.SetPoints(points)

    # --------------------------------------------------------
    # Initial spins
    # --------------------------------------------------------
    emom0 = sim.get_emom()
    spins = emom0[:, :, 0].T
    vec = numpy_support.numpy_to_vtk(spins, deep=False)
    mag = numpy_support.numpy_to_vtk(spins[:, 2], deep=False)

    poly.GetPointData().SetVectors(vec)
    poly.GetPointData().SetScalars(mag)

    # --------------------------------------------------------
    # Colormap
    # --------------------------------------------------------
    lut = create_lookup_table()
    color_series = vtkColorSeries()

    # --------------------------------------------------------
    # Glyphs
    # --------------------------------------------------------
    arrow = vtk.vtkArrowSource()
    glyph = vtk.vtkGlyph3DMapper()
    glyph.SetInputData(poly)
    glyph.SetSourceConnection(arrow.GetOutputPort())
    glyph.OrientOn()
    glyph.SetScaleFactor(2.0)
    glyph.SetScaleModeToNoDataScaling()
    glyph.SetLookupTable(lut)
    glyph.SetColorModeToMapScalars()

    vec_actor = vtk.vtkActor()
    vec_actor.SetMapper(glyph)

    # --------------------------------------------------------
    # Text overlays
    # --------------------------------------------------------
    def make_text(y):
        t = vtk.vtkTextActor()
        t.GetTextProperty().SetFontSize(36)
        t.GetTextProperty().SetColor(0, 0, 0)
        t.SetDisplayPosition(20, y)
        return t

    txt_T = make_text(620)
    txt_B = make_text(570)
    txt_E = make_text(40)

    ren.AddActor(vec_actor)
    ren.AddActor(txt_T)
    ren.AddActor(txt_B)
    ren.AddActor(txt_E)

    # --------------------------------------------------------
    # Update helpers
    # --------------------------------------------------------

    def update_vectors(moments=None):
        # moments may be provided (from sim.step) or None -> query sim
        if moments is None:
            emom = sim.get_emom()  # (3, N, M)
            vec = emom[:, :, 0].T  # (N, 3)
            print("update_vectors: queried sim.get_emom(), emom.shape=", emom.shape, flush=True)
        else:
            m = np.asarray(moments)
            print("update_vectors: received moments with shape", m.shape, flush=True)
            if m.ndim == 3:
                vec = m[:, :, 0].T
            elif m.ndim == 2:
                vec = m.T
            else:
                raise ValueError("Unexpected moments shape")

        vec = np.asarray(vec, dtype=np.float64, order="C")
        vtk_vec = numpy_to_vtk(vec, deep=True)
        poly.GetPointData().SetVectors(vtk_vec)
        poly.Modified()

    def update_text():
        T = sim.get_temperature()
        Bz = sim.get_field()[2]
        E = sim.calculate_energy()

        txt_T.SetInput(f"T = {T:.3f} K")
        txt_B.SetInput(f"Bz = {Bz:.2f} T")
        txt_E.SetInput(f"E = {E:.5f} mRy/atom" if E else "E = ---")

    # --------------------------------------------------------
    # Key handler
    # --------------------------------------------------------
    def on_keypress(obj, evt):
        key = obj.GetKeySym()
        print(f"Key pressed: {key}", flush=True)

        if key == "P":
            screenshot(renwin)

        elif key == "0":
            sim.reset_moments()
            update_vectors()

        elif key in ("S", "s"):
            print("Action: stepping LLG (S) 50 steps...", flush=True)
            before = sim.get_emom().copy()
            print("before[0:3,0:5,0] min/max:", before.min(), before.max(), flush=True)
            moments = sim.step("S", temperature=None, nstep=50)
            print("sim.step returned", None if moments is None else getattr(moments, 'shape', type(moments)), flush=True)
            after = sim.get_emom().copy()
            print("after[0:3,0:5,0] min/max:", after.min(), after.max(), flush=True)
            # show a small delta summary
            try:
                delta = np.linalg.norm(after - before)
                print("||after-before||_F =", delta, flush=True)
            except Exception:
                pass
            update_vectors(moments)

        elif key in ("M", "m"):
            print("Action: stepping Metropolis (M) 50 steps...", flush=True)
            before = sim.get_emom().copy()
            print("before min/max:", before.min(), before.max(), flush=True)
            moments = sim.step("M", temperature=None, nstep=50)
            print("sim.step returned", None if moments is None else getattr(moments, 'shape', type(moments)), flush=True)
            after = sim.get_emom().copy()
            print("after min/max:", after.min(), after.max(), flush=True)
            try:
                delta = np.linalg.norm(after - before)
                print("||after-before||_F =", delta, flush=True)
            except Exception:
                pass
            update_vectors(moments)

        elif key in ("H", "h"):
            print("Action: stepping Heat-bath (H) 50 steps...", flush=True)
            before = sim.get_emom().copy()
            print("before min/max:", before.min(), before.max(), flush=True)
            moments = sim.step("H", temperature=None, nstep=50)
            print("sim.step returned", None if moments is None else getattr(moments, 'shape', type(moments)), flush=True)
            after = sim.get_emom().copy()
            print("after min/max:", after.min(), after.max(), flush=True)
            try:
                delta = np.linalg.norm(after - before)
                print("||after-before||_F =", delta, flush=True)
            except Exception:
                pass
            update_vectors(moments)

        elif key == "B":
            sim.add_field([0, 0, 1])
        elif key == "b":
            sim.add_field([0, 0, -1])
        elif key == "N":
            sim.add_field([0, 0, 10])
        elif key == "n":
            sim.add_field([0, 0, -10])

        elif key == "T":
            sim.set_temperature(sim.get_temperature() + 1)
        elif key == "t":
            sim.set_temperature(max(sim.get_temperature() - 1, 1e-6))
        elif key == "Y":
            sim.set_temperature(sim.get_temperature() + 10)
        elif key == "y":
            sim.set_temperature(max(sim.get_temperature() - 10, 1e-6))

        elif key in ("c", "C"):
            idx = color_series.GetColorScheme()
            idx += -1 if key == "c" else 1
            color_series.SetColorScheme(idx % color_series.GetNumberOfColorSchemes())
            color_series.BuildLookupTable(lut)

        elif key in ("i", "I"):
            emom = sim.get_emom()
            mz = emom[2, :, 0]
            n = int(np.sqrt(len(mz)))
            mat = mz.reshape((n, n))
            fq = np.abs(np.fft.fftshift(np.fft.fft2(mat)))
            plot_static_correlation(fq, do_log=(key == "I"))

        update_text()
        renwin.Render()

    # --------------------------------------------------------
    # Start interaction
    # --------------------------------------------------------
    iren = vtk.vtkRenderWindowInteractor()
    iren.SetRenderWindow(renwin)
    iren.AddObserver("KeyPressEvent", on_keypress)

    update_text()
    renwin.Render()
    iren.Start()


if __name__ == "__main__":
    main()
