#!/usr/bin/env python3
"""
Script for creating a snapshot from SD-data.

This script reads atom positions and vector data from files and visualizes them
using VTK (Visualization Toolkit).
The visualization includes colored spheres representing atoms and arrows representing vectors.

The script also includes functions for creating a lookup table for coloring the vectors,
cycling through different color schemes, and taking screenshots of the visualization.

Author: Anders Bergman
"""

# standard libs
import argparse
import logging
import os

import numpy as np
# import uppasd.pyasd as asd
import uppasd.simulator as sim
try:
    import vtk
    from vtk.util import numpy_support
    from vtkmodules.vtkCommonColor import vtkColorSeries
except ImportError:
    print("VTK is not installed. Please install it using 'pip install vtk'.")
    exit()
try:
    from scipy.ndimage import gaussian_filter
except ImportError:
    print("Scipy is not installed. Please install it using 'pip install scipy'.")
    exit()
try:
    import matplotlib.cm as cm
    import matplotlib.pyplot as plt
except ImportError:
    print("Matplotlib is not installed. Please install it using 'pip install matplotlib'.")
    exit()


def close_window(iren):
    """
    Closes the window and terminates the application.

    Parameters:
    - iren: vtkRenderWindowInteractor object

    Returns:
    None
    """
    render_window = iren.GetRenderWindow()
    render_window.Finalize()
    iren.TerminateApp()


# Read Location of Atoms
def readAtoms(file, Nmax):
    """
    Read atoms from a file and return a vtkPoints object and the number of atoms.

    Parameters:
    - file: The file object to read atoms from.
    - Nmax: The maximum number of atoms to read.

    Returns:
    - points: A vtkPoints object containing the coordinates of the atoms.
    - nrAtoms: The number of atoms read from the file.
    """
    points = vtk.vtkPoints()
    nrAtoms = 0
    # Read all data defensively
    for line in file:
        if not line:
            break
        line = line.strip()
        if not line:
            continue
        parts = line.split()
        # need at least 4 columns: index, x, y, z
        if len(parts) < 4:
            # malformed line, stop or skip
            continue
        try:
            # a = int(parts[0])  # original atom index, unused here
            x, y, z = float(parts[1]), float(parts[2]), float(parts[3])
        except ValueError:
            # skip lines that do not parse
            continue
        if nrAtoms < Nmax:
            points.InsertPoint(nrAtoms, x, y, z)
        nrAtoms += 1
    return points, nrAtoms


# Read vectors
# We must know the time step and the number of atoms per time
def readVectorsData(file, time, nrAtoms, Nmax):
    """
    Read vector data from a file and return vtkFloatArray objects representing the vectors and colors.

    Parameters:
    - file: The file object to read the data from.
    - time: The time value associated with the data.
    - nrAtoms: The total number of atoms in the data.
    - Nmax: The maximum number of atoms to read.

    Returns:
    - vectors: A vtkFloatArray object representing the vectors.
    - colors: A vtkFloatArray object representing the colors.
    """
    # Create arrays which represent the vectors and scalar colors
    vectors = vtk.vtkFloatArray()
    colors = vtk.vtkFloatArray()

    # Define number of elements
    vectors.SetNumberOfComponents(3)
    colors.SetNumberOfComponents(1)

    # Skip header lines if present (defensive)
    for _ in range(7):
        header_line = file.readline()
        if not header_line:
            break

    i = 0
    # Read all data for a certain time
    while i < nrAtoms:
        line = file.readline()
        if not line:
            break
        line = line.strip()
        if not line:
            continue
        parts = line.split()
        # expect at least 7 columns based on previous code
        if len(parts) < 7:
            # malformed line, skip
            continue
        try:
            # validate integer fields exist
            int(parts[0])
            int(parts[1])
            x = float(parts[4])
            y = float(parts[5])
            z = float(parts[6])
        except ValueError:
            # skip lines that don't parse
            continue
        if i < Nmax:
            m = (z + 1.0) / 2.0
            vectors.InsertTuple3(i, x, y, z)
            colors.InsertValue(i, m)
        i += 1
    return vectors, colors


def plot_static_correlation(fsmat, do_log=False):
    """
    Plot the static correlation S(q).

    Parameters:
    - fsmat: numpy array
        The correlation matrix to be plotted.
    - do_log: bool, optional
        Whether to apply logarithmic scaling to the plot. Default is False.

    Returns:
    None
    """
    if do_log:
        fsmat = np.log(fsmat + 1.0)
    else:
        fsmat = np.abs(fsmat)

    plt.figure(1)
    plt.clf()
    plt.xticks([], [])
    plt.yticks([], [])
    plt.imshow(fsmat, cmap=cm.Reds)
    plt.title(r"Static correlation S(q)")
    plt.xlabel(r"$q_x$")
    plt.ylabel(r"$q_y$")
    plt.ion()
    plt.show()


def createLookupTable():
    """
    Creates and returns a vtkLookupTable object with a custom color mapping.

    The lookup table is created with 256 colors, where the first half of the table
    maps from red to green and the second half maps from green to blue.

    Returns:
        vtkLookupTable: The created lookup table object.
    """
    lut = vtk.vtkLookupTable()
    for i in range(0, 128, 1):
        lut.SetTableValue(i, (127.0 - i) / 127.0, i / 127.0, 0, 1)
    for i in range(128, 256, 1):
        lut.SetTableValue(i, 0, (256.0 - i) / 128.0, (i - 128.0) / 128, 1)
    lut.SetTableRange(-1.0, 1.0)
    return lut


def cycleColorScheme(lut, colorSeries, backwards=False):
    """
    Cycle through different color schemes in a color series.

    Args:
        lut (vtkLookupTable): The lookup table to be updated.
        colorSeries (vtkColorSeries): The color series to cycle through.
        backwards (bool, optional): If True, cycle backwards. Defaults to False.

    Returns:
        tuple: A tuple containing the updated lookup table and color series.
    """
    cterm = -1 if backwards else 1
    colorIndex = colorSeries.GetColorScheme() + cterm
    if colorIndex <= 0:
        colorSeries.SetColorScheme(colorSeries.GetNumberOfColorSchemes() - 1)
        lut = createLookupTable()
    elif colorIndex >= colorSeries.GetNumberOfColorSchemes():
        colorSeries.SetColorScheme(0)
        lut = createLookupTable()
    else:
        colorSeries.SetColorScheme(colorIndex)
        colorSeries.BuildLookupTable(lut, vtkColorSeries.ORDINAL)
    return lut, colorSeries


number_of_screenshots = 1


# A function that takes a renderwindow and saves its contents to a .png file
def Screenshot(renWin):
    """
    Takes a screenshot of the specified render window and saves it as a POV file and a PNG file.

    Parameters:
    - renWin (vtkRenderWindow): The render window to take a screenshot of.
    - number_of_screenshots (int): The current count of screenshots.

    Returns:
    - number_of_screenshots (int): The updated count of screenshots.
    """
    global number_of_screenshots
    win2im = vtk.vtkWindowToImageFilter()
    win2im.ReadFrontBufferOff()
    win2im.SetInput(renWin)
    #
    povexp = vtk.vtkPOVExporter()
    povexp.SetRenderWindow(renWin)
    # povexp.SetInput(renWin)
    renWin.Render()
    povexp.SetFileName("snap%.5d.pov" % number_of_screenshots)
    povexp.Write()
    #
    toPNG = vtk.vtkPNGWriter()
    toPNG.SetFileName("snap%.5d.png" % number_of_screenshots)
    # toPNG.SetInput(win2im.GetOutput())
    toPNG.SetInputConnection(win2im.GetOutputPort())
    toPNG.Write()

    number_of_screenshots += 1

    return number_of_screenshots


# ---------------------------
# Screenshot code ends here
# ---------------------------


# ---------------------------
# Timer code starts here
# ---------------------------
class vtkTimerCallback:
    """
    A class representing a timer callback for VTK.

    This class is used to define a callback function that will be executed
    when a timer event occurs in VTK.

    Attributes:
        timer_count (int): The count of timer events that have occurred.
    """

    def __init__(self):
        self.timer_count = 0


# For the interactive control. Set up a check for aborting rendering.
def CheckAbort(obj, event):
    # obj will be the object generating the event.  In this case it
    # is renWin.
    if obj.GetEventPending() != 0:
        obj.SetAbortRender(1)


def main():

    # parse command line options
    parser = argparse.ArgumentParser(description="Interactive UppASD VTK visualizer")
    parser.add_argument("--coords", help="Path to coord file (overrides simulator coords)")
    parser.add_argument("--moments", help="Path to moments file (not yet used)")
    parser.add_argument("--screenshots-dir", help="Directory to write screenshots", default=None)
    parser.add_argument("--verbose", help="Enable debug logging", action="store_true")
    args = parser.parse_args()

    # configure logging
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s: %(message)s")
    if args.verbose:
        logging.getLogger().setLevel(logging.DEBUG)


    # Instantiate the ASD simulation object
    asd = sim.Simulator()
    asd.init_simulation()

    renWin = vtk.vtkRenderWindow()
    ren = vtk.vtkOpenGLRenderer()
    renWin.AddRenderer(ren)

    # Set color of backgroung
    # ren.SetBackground(0.0, 0.0, 0.0)
    ren.SetBackground(1.0, 1.0, 1.0)
    renWin.SetSize(1048, 640)
    # renWin.SetSize(960,960)

    # Datatest=vtk.vtkUnstructuredGrid()
    Datatest = vtk.vtkPolyData()
    # DataScal=vtk.vtkUnstructuredGrid()

    DataFour = vtk.vtkPolyData()

    # ----------------------------
    # Screenshot code begins here
    # ----------------------------

    # ASD STUFF
    asd.get_moments()
    asd.get_coords()

    initmom = asd.moments[:, :, 0]
    currmom = asd.moments[:, :, 0].T.copy()
    vecz = numpy_support.numpy_to_vtk(currmom, deep=True)

    initcol = asd.moments[2, :, 0]
    currcol = asd.moments[2, :, 0].T.copy()
    colz = numpy_support.numpy_to_vtk(currcol, deep=True)

    # vecz,colz=(readVectorsData(directionsFile,0,nrAtoms,Nmax))

    # Size of system (unused in this path)

    # Viewing distance
    # decrease number to zoom out and vice versa
    # dollyA = 0.014
    # dollyB = 0.000
    # Open files
    # momfiles = glob.glob("restart.????????.out")
    # directionsFile = open(momfiles[0])
    # directionsFile = open("momentsparsed.out")
    logging.debug("Skip reading coords from file unless --coords is provided")
    # posfiles = glob.glob("coord.????????.out")
    # # print posfiles
    # atomsFile = open(posfiles[0])
    # # atomsFile = open("atomsparsed.out")

    # Read atom positions
    # If the user provided a coords file, use it; otherwise use asd.coords
    if args.coords:
        if not os.path.exists(args.coords):
            logging.error("Coords file not found: %s", args.coords)
            return
        with open(args.coords, "r") as atomsFile:
            atomDataPts, nrAtoms = readAtoms(atomsFile, 10_000_000)
        logging.info("Read %d atom positions from %s", nrAtoms, args.coords)
        Datatest.SetPoints(atomDataPts)
        DataFour.SetPoints(atomDataPts)
    else:
        atomPoints = vtk.vtkPoints()
        atomData = numpy_support.numpy_to_vtk(asd.coords.T.copy(), deep=True)
        atomPoints.SetData(atomData)
        nrAtoms = asd.natom
        logging.info("Number of atoms (from simulator): %d", nrAtoms)
        Datatest.SetPoints(atomPoints)
        DataFour.SetPoints(atomPoints)
    # atomData, nrAtoms = readAtoms(atomsFile, Nmax)
    # # print nrAtoms
    # Datatest.SetPoints(atomData)
    # # DataScal.SetPoints(atomData)
    # DataFour.SetPoints(atomData)

    # Read data for vectors
    # for i in range(0,55,1):
    # vecz,colz=(readVectorsData(directionsFile,0,nrAtoms,Nmax))
    Datatest.GetPointData().SetVectors(vecz)
    Datatest.GetPointData().SetScalars(colz)
    # DataScal.GetPointData().SetScalars(colz)

    DataFour.GetPointData().SetVectors(vecz)
    DataFour.GetPointData().SetScalars(colz)

    # create a small status text actor to indicate busy state (defined early so helpers can use it)
    status_txt = vtk.vtkTextActor()
    status_txt.SetInput("")
    status_txt_prop = status_txt.GetTextProperty()
    status_txt_prop.SetFontFamilyToArial()
    status_txt_prop.SetFontSize(28)
    status_txt_prop.SetColor(0, 0, 0)
    status_txt_prop.BoldOn()
    status_txt.SetDisplayPosition(20, 20)

    # Helper to run relax, fetch moments and update VTK arrays
    def run_relax_and_update(mode=None, temperature=None, take_snapshot=False):
        if mode is not None:
            # show busy status (guard in case status_txt isn't created yet)
            try:
                status_txt.SetInput("Running...")
                renWin.Render()
            except NameError:
                pass
            # prefer explicit float for temperature
            temp = temperature if temperature is not None else float(asd.inputdata.iptemp)
            asd.relax(mode=mode, temperature=temp)
            asd.get_moments()

        # update VTK arrays with deep copies to avoid dangling memory refs
        currmom_local = asd.moments[:, :, 0].T.copy()
        currcol_local = asd.moments[2, :, 0].T.copy()
        vec_local = numpy_support.numpy_to_vtk(currmom_local, deep=True)
        col_local = numpy_support.numpy_to_vtk(currcol_local, deep=True)
        Datatest.GetPointData().SetVectors(vec_local)
        Datatest.GetPointData().SetScalars(col_local)
        renWin.Render()
        if take_snapshot:
            Screenshot(renWin)
        # clear status (guard in case status_txt isn't created yet)
        try:
            status_txt.SetInput("")
            renWin.Render()
        except NameError:
            pass

    # Create colortable for the coloring of the vectors
    lut = vtk.vtkLookupTable()
    colorSeries = vtkColorSeries()
    # colorSeries.SetColorScheme(vtkColorSeries.BREWER_QUALITATIVE_PASTEL1)
    colorSeries.SetColorScheme(0)
    lut = createLookupTable()
    lut.Build()

    print(
        "vtkColorSeries ID: ",
        colorSeries.GetColorScheme(),
        colorSeries.GetNumberOfColorSchemes(),
    )

    # ATOMS  ##########################################
    # Set up atoms
    ball = vtk.vtkSphereSource()
    ball.SetRadius(1.00)
    ball.SetThetaResolution(24)
    ball.SetPhiResolution(24)

    balls = vtk.vtkGlyph3DMapper()
    balls.SetInputData(Datatest)
    balls.SetSourceConnection(ball.GetOutputPort())
    balls.SetScaleFactor(0.57)
    balls.SetScaleModeToNoDataScaling()
    balls.SetLookupTable(lut)
    balls.Update()

    atom = vtk.vtkLODActor()
    atom.SetMapper(balls)
    atom.GetProperty().SetOpacity(1.5)
    xmin, xmax = atom.GetXRange()
    ymin, ymax = atom.GetYRange()
    zmin, zmax = atom.GetZRange()
    xmid = (xmin + xmax) / 2
    ymid = (ymin + ymax) / 2
    zmid = (zmin + zmax) / 2
    # Good position for "z-direction"
    atom.SetPosition(-xmid, -ymid, -zmid)
    # Good position for "xy-direction"
    # atom.SetPosition(-xmid*1.4,-ymid,-zmid)

    # Set the color
    atom.GetProperty().SetSpecular(0.3)
    atom.GetProperty().SetSpecularPower(60)
    atom.GetProperty().SetAmbient(0.2)
    atom.GetProperty().SetDiffuse(0.8)
    atom.GetProperty().SetInterpolationToPBR()
    atom.GetProperty().SetRoughness(0.6)
    atom.GetProperty().SetMetallic(0.7)
    atom.GetProperty().ShadingOn()

    # ATOMS F #########################################
    # Set up atoms
    fball = vtk.vtkSphereSource()
    fball.SetRadius(1.00)
    fball.SetThetaResolution(24)
    fball.SetPhiResolution(24)

    fballs = vtk.vtkGlyph3DMapper()
    fballs.SetInputData(DataFour)
    fballs.SetSourceConnection(fball.GetOutputPort())
    fballs.SetScaleFactor(0.57)
    fballs.SetScaleModeToNoDataScaling()
    fballs.SetLookupTable(lut)
    fballs.Update()

    four = vtk.vtkLODActor()
    four.SetMapper(fballs)
    four.GetProperty().SetOpacity(1.5)
    xmin, xmax = four.GetXRange()
    ymin, ymax = four.GetYRange()
    zmin, zmax = four.GetZRange()
    xmid = (xmin + xmax) / 2
    ymid = (ymin + ymax) / 2
    zmid = (zmin + zmax) / 2
    # Good position for "z-direction"
    four.SetPosition(-xmid, -ymid, -zmid)
    # Good position for "xy-direction"
    # four.SetPosition(-xmid*1.4,-ymid,-zmid)

    # Set the color
    four.GetProperty().SetSpecular(0.3)
    four.GetProperty().SetSpecularPower(60)
    four.GetProperty().SetAmbient(0.2)
    four.GetProperty().SetDiffuse(0.8)

    # ARROWS ##########################################
    # Create vectors
    arrow = vtk.vtkArrowSource()
    arrow.SetTipRadius(0.26)
    arrow.SetShaftRadius(0.15)
    arrow.SetTipResolution(18)
    arrow.SetShaftResolution(18)

    # Calculate normals for shading
    arrownormals = vtk.vtkPolyDataNormals()
    arrownormals.SetInputConnection(arrow.GetOutputPort())

    # Calculate TCoords for texturing
    arrownormalstmap = vtk.vtkTextureMapToCylinder()
    arrownormalstmap.SetInputConnection(arrownormals.GetOutputPort())
    arrownormalstmap.PreventSeamOn()

    glyph = vtk.vtkGlyph3DMapper()
    # Set input connection to the normals
    # glyph.SetSourceConnection(arrow.GetOutputPort())
    glyph.SetSourceConnection(arrownormalstmap.GetOutputPort())
    glyph.OrientOn()
    # glyph.Update()
    # glyph.SetInputConnection(Datatest.GetProducerPort())
    glyph.SetInputData(Datatest)
    # glyph.SetInput(Datatest) # Position and direction
    glyph.SetScaleFactor(2.00)
    # Color the vectors according to magnitude
    # glyph.SetVectorModeToUseVector()
    # glyph.SetColorModeToColorByVector()
    glyph.SetScaleModeToNoDataScaling()
    # glyph.Update()

    # glyphMapper = vtk.vtkPolyDataMapper()
    # glyphMapper.SetInput(glyph.GetOutput())
    # glyphMapper.SetLookupTable(lut)
    glyph.SetLookupTable(lut)
    # glyph.ColorByArrayComponent(0,0)
    # glyph.UseLookupTableScalarRangeOn()
    glyph.SetColorModeToMapScalars()
    # glyph.SetScalarModeToUseFieldData()
    # glyph.ScalarVisibilityOn()
    glyph.Update()

    vector = vtk.vtkLODActor()
    vector.SetMapper(glyph)
    # Good position for "z-direction"
    vector.SetPosition(-xmid, -ymid, -zmid)
    # Good position for "xy-direction"
    # vector.SetPosition(-xmid*1.4,-ymid,-zmid)

    # vector.GetProperty().SetInterpolationToPBR()
    vector.GetProperty().SetSpecular(0.0)
    vector.GetProperty().SetSpecularPower(10)
    vector.GetProperty().SetAmbient(0.3)
    vector.GetProperty().SetDiffuse(0.5)
    vector.GetProperty().SetOpacity(1.0)
    vector.GetProperty().SetRoughness(0.7)
    vector.GetProperty().SetMetallic(0.4)
    vector.GetProperty().ShadingOn()
    vector.GetProperty().SetInterpolationToGouraud()
    print("Interpolation: ", vector.GetProperty().GetInterpolationAsString())

    # Bounding box
    outlineData = vtk.vtkOutlineFilter()
    # outlineData.SetInputConnection(Datatest.GetProducerPort())
    outlineData.SetInputData(Datatest)
    outlineMapper = vtk.vtkPolyDataMapper()
    outlineMapper.SetInputConnection(outlineData.GetOutputPort())
    # outlineMapper.SetInput(outlineData.GetOutput())
    outline = vtk.vtkActor()
    outline.SetMapper(outlineMapper)
    outline.GetProperty().SetColor(0, 0, 0)
    outline.GetProperty().SetLineWidth(5.0)
    outline.SetPosition(-xmid, -ymid, -zmid)

    # Text
    # create a text actor for Temperature
    temptxt = vtk.vtkTextActor()
    # temp='{:4.3f}'.format(asd.inputdata.iptemp[0])
    temp = "{:4.3f}".format(0)
    temptxt.SetInput("T = " + temp + " K")
    temptxtprop = temptxt.GetTextProperty()
    temptxtprop.SetFontFamilyToArial()
    temptxtprop.SetFontSize(56)
    temptxtprop.SetColor(0, 0, 0)
    temptxtprop.BoldOn()
    temptxt.SetDisplayPosition(20, 550)

    # create a text actor for Field
    fieldtxt = vtk.vtkTextActor()
    bz = "{:4.1f}".format(0)
    # bz='{:4.1f}'.format(asd.inputdata.iphfield[2])
    fieldtxt.SetInput("Bz = " + bz + " T")
    fieldtxtprop = fieldtxt.GetTextProperty()
    fieldtxtprop.SetFontFamilyToArial()
    fieldtxtprop.SetFontSize(56)
    fieldtxtprop.SetColor(0, 0, 0)
    fieldtxtprop.BoldOn()
    fieldtxt.SetDisplayPosition(20, 500)

    # create a text actor for Energy
    enetxt = vtk.vtkTextActor()
    ene = "{:6.4f}".format(0)
    # ene='{:6.4f}'.format(asd.simulationdata.total_energy)
    enetxt.SetInput("E= " + ene + " mRy/atom")
    enetxtprop = enetxt.GetTextProperty()
    enetxtprop.SetFontFamilyToArial()
    enetxtprop.SetFontSize(36)
    enetxtprop.SetColor(0, 0, 0)
    enetxtprop.BoldOn()
    enetxt.SetDisplayPosition(20, 50)

    # (status_txt already created earlier)

    # LIGHTS ON
    light = vtk.vtkLight()
    light.SetColor(1.0, 1.0, 1.0)
    ren.AddLight(light)

    # Reposition the camera
    # z-direction
    ren.GetActiveCamera().Azimuth(0)
    ren.GetActiveCamera().Elevation(0)
    # xy-direction (almost)
    ren.GetActiveCamera().ParallelProjectionOn()
    # ren.GetActiveCamera().ParallelProjectionOff()
    # d = ren.GetActiveCamera().GetDistance()
    # ren.GetActiveCamera().SetParallelScale(0.55*max(zmax,xmax))
    ren.GetActiveCamera().SetFocalPoint(0, 0, 0)
    ren.GetActiveCamera().SetViewUp(-0.866025403784439, 0.5, 0)
    # ren.GetActiveCamera().SetViewUp(0,1,0)
    ren.GetActiveCamera().SetParallelScale(0.55 * ymax)
    # ren.GetActiveCamera().SetPosition(0,0,-100*d)
    # print ren.GetActiveCamera().GetViewAngle()
    l_vp = max(xmax - xmin, zmax - zmin) / 2
    h_vp = l_vp / 0.26795 * 1.1
    # print l,h

    ren.GetActiveCamera().SetPosition(0, 0, h_vp)
    # ren.GetActiveCamera().SetPosition(0,0,50*d)

    # Enable antialiasing
    ren.UseFXAAOn()
    ren.GetFXAAOptions().SetUseHighQualityEndpoints(True)
    renWin.SetMultiSamples(4)
    # Enable SSAO
    ren.UseSSAOOn()
    ren.SetSSAOKernelSize(512)
    ren.SetSSAORadius(3.0)
    ren.SetSSAOBias(0.1)
    ren.SSAOBlurOff()

    renWin.AddObserver("AbortCheckEvent", CheckAbort)

    # Add the actors to the renderer, set the background and size
    # Atoms
    atom.SetVisibility(0)
    ren.AddActor(atom)
    # ren.AddActor(txt)
    # Vectors
    ren.AddActor(vector)
    # Text
    ren.AddActor(temptxt)
    ren.AddActor(fieldtxt)
    ren.AddActor(enetxt)
    # Status text (busy indicator)
    ren.AddActor(status_txt)

    # Update initial text fields from simulator so they show correct values
    try:
        # ensure simulator has calculated energy
        asd.calculate_energy()
    except Exception:
        # if calculate_energy fails, continue silently
        pass

    try:
        temp_val = float(asd.inputdata.iptemp)
    except Exception:
        temp_val = 0.0
    try:
        bz_val = float(asd.inputdata.iphfield[2])
    except Exception:
        # fallback if iphfield not available or not indexable
        try:
            bz_val = float(asd.inputdata.iphfield)
        except Exception:
            bz_val = 0.0

    temp = "{:4.3f}".format(temp_val)
    temptxt.SetInput("T = " + temp + " K")
    bz = "{:4.1f}".format(bz_val)
    fieldtxt.SetInput("Bz = " + bz + " T")
    ene = "{:6.4f}".format(getattr(asd, "energy", 0.0))
    enetxt.SetInput("E= " + ene + " mRy/atom")

    # Outline
    # ren.AddActor(outline)
    # Scalar bar
    # ren.AddActor(scalarBar)

    # a minimal keyboard interface
    def Keypress(obj, event):
        """
        Handle keypress events.

        Parameters:
        - obj: The object that triggered the event.
        - event: The keypress event.

        Returns:
        None
        """
        global number_of_screenshots
        key = obj.GetKeySym()
        if key == "P":
            Screenshot(renWin)
            print("Screenshot taken")
        if key == "0":
            print("Resetting UppASD")
            asd.put_moments(initmom[:, :])
            vecz = numpy_support.numpy_to_vtk(initmom[:, :].T.copy(), deep=True)
            colz = numpy_support.numpy_to_vtk(initcol.copy(), deep=True)
            Datatest.GetPointData().SetVectors(vecz)
            Datatest.GetPointData().SetScalars(colz)
        if key == "M" or key == "m":
            print("Running UppASD (M)")
            run_relax_and_update(mode="M", temperature=float(asd.inputdata.iptemp), take_snapshot=(key == "M"))
        if key == "H" or key == "h":
            print("Running UppASD (H)")
            run_relax_and_update(mode="H", temperature=float(asd.inputdata.iptemp) + 1.0e-6, take_snapshot=(key == "H"))
        if key == "S" or key == "s":
            print("Running UppASD (S)")
            run_relax_and_update(mode="S", temperature=float(asd.inputdata.iptemp), take_snapshot=(key == "S"))
        if key == "c" or key == "C":
            cycleColorScheme(lut, colorSeries, backwards=(key == "c"))
            lut.Build()
        # if key == "C":
        #   asd.inputdata.ipmode='Q'
        #   asd.pyasd.initialphase()
        #   currmom=asd.momentdata.emom[:,:,0].T
        #   currcol=asd.momentdata.emom[2,:,0].T
        #   vecz=numpy_support.numpy_to_vtk(currmom)
        #   colz=numpy_support.numpy_to_vtk(currcol)
        #   Datatest.GetPointData().SetVectors(vecz)
        #   Datatest.GetPointData().SetScalars(colz)
        # if key == "X":
        #   asd.inputdata.ipmode='Y'
        #   asd.pyasd.initialphase()
        #   currmom=asd.momentdata.emom[:,:,0].T
        #   currcol=asd.momentdata.emom[2,:,0].T
        #   vecz=numpy_support.numpy_to_vtk(currmom)
        #   colz=numpy_support.numpy_to_vtk(currcol)
        #   Datatest.GetPointData().SetVectors(vecz)
        #   Datatest.GetPointData().SetScalars(colz)
        if key == "B":
            asd.inputdata.iphfield[2] = asd.inputdata.iphfield[2] + 1.0
            asd.inputdata.update_iphfield()
            bz = "{:4.1f}".format(asd.inputdata.iphfield[2])
            fieldtxt.SetInput("Bz = " + bz + " T")
        if key == "b":
            asd.inputdata.iphfield[2] = asd.inputdata.iphfield[2] - 1.0
            asd.inputdata.update_iphfield()
            bz = "{:4.1f}".format(asd.inputdata.iphfield[2])
            fieldtxt.SetInput("Bz = " + bz + " T")
        if key == "N":
            asd.inputdata.iphfield[2] = asd.inputdata.iphfield[2] + 10.0
            asd.inputdata.update_iphfield()
            bz = "{:4.1f}".format(asd.inputdata.iphfield[2])
            fieldtxt.SetInput("Bz = " + bz + " T")
        if key == "n":
            asd.inputdata.iphfield[2] = asd.inputdata.iphfield[2] - 10.0
            asd.inputdata.update_iphfield()
            bz = "{:4.1f}".format(asd.inputdata.iphfield[2])
            fieldtxt.SetInput("Bz = " + bz + " T")
        if key == "T":
            asd.inputdata.iptemp = asd.inputdata.iptemp + 1.0
            asd.inputdata.update_iptemp()
            temp = "{:4.3f}".format(asd.inputdata.iptemp)
            temptxt.SetInput("T = " + temp + " K")
        if key == "t":
            asd.inputdata.iptemp = max(asd.inputdata.iptemp - 1.0, 1.0e-5)
            asd.inputdata.update_iptemp()
            temp = "{:4.3f}".format(asd.inputdata.iptemp)
            temptxt.SetInput("T = " + temp + " K")
        if key == "Y":
            asd.inputdata.iptemp = asd.inputdata.iptemp + 10.0
            asd.inputdata.update_iptemp()
            temp = "{:4.3f}".format(asd.inputdata.iptemp)
            temptxt.SetInput("T = " + temp + " K")
        if key == "y":
            asd.inputdata.iptemp = max(asd.inputdata.iptemp - 10.0, 1.0e-5)
            asd.inputdata.update_iptemp()
            temp = "{:4.3f}".format(asd.inputdata.iptemp)
            temptxt.SetInput("T = " + temp + " K")
        if key == "g":
            vector.GetProperty().SetInterpolationToGouraud()
        if key == "p":
            vector.GetProperty().SetInterpolationToPBR()
        # if key == "m":
        #   vector.SetVisibility(abs(1-vector.GetVisibility()))
        if key == "i" or key == "I":
            # print('-....-')
            # cdata = asd.momentdata.emom[2, :, 0].T
            asd.get_moments()
            cdata = asd.moments[2, :, 0].T
            xdim = int(np.sqrt(cdata.shape)[0])
            cmat = cdata.reshape((xdim, xdim))
            fmat = np.abs(np.fft.fft2(cmat.T))
            fmat[0, 0] = 0.0
            fsgau = gaussian_filter(fmat, sigma=1.0)
            fsmat = np.fft.fftshift(fsgau)
            fsmat = np.log(fsmat + 1.0)

            # Call the function
            plot_static_correlation(fsmat, do_log=(key == "i"))

        if key == "a":
            atom.SetVisibility(abs(1 - atom.GetVisibility()))

        if key == "v":
            vector.SetVisibility(abs(1 - atom.GetVisibility()))

        asd.calculate_energy()
        ene = "{:6.4f}".format(asd.energy)
        enetxt.SetInput("E= " + ene + " mRy/atom")
        renWin.Render()

    # Render scene
    iren = vtk.vtkRenderWindowInteractor()
    istyle = vtk.vtkInteractorStyleTrackballCamera()
    iren.SetInteractorStyle(istyle)
    iren.SetRenderWindow(renWin)
    iren.AddObserver("KeyPressEvent", Keypress)

    # ren.ResetCamera()
    iren.Initialize()
    #
    #
    # cb = vtkTimerCallback()
    # ##cb.AddActor = vector
    # iren.AddObserver('TimerEvent', cb.execute)
    # timerId = iren.CreateRepeatingTimer(100);
    # iren.SetStillUpdateRate(0.050)
    # iren.SetDesiredUpdateRate(0.050)

    renWin.Render()
    iren.Start()
    Screenshot(renWin)


if __name__ == "__main__":
    main()
