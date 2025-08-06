#!/usr/bin/env python3
# Script for creating a snapshot from SD-data
# Written by Anders Bergman, after template from Anders Hast
#
# This particular script is suited for an x-y plane of atoms/moments
# Coloring of magnetic moments is determined by their z-components
# to change this, modify the value m as read in ReadTimeData.py
#
# Zooming can be modified by the dollyA parameter
# Camera positioning can be changed using GetActiveCamera.Elevation, Roll, and Azimuth
#

import glob
import string
from math import acos, atan2

import vtk

# from ReadTimeData import * # Functions that let us read the files


# Read Location of Atoms
def readAtoms(file, Nmax):

    points = vtk.vtkPoints()
    nrAtoms = 0
    # Read ahead
    line = file.readline()
    data = line.split()

    # Read all data
    while line:
        if nrAtoms <= Nmax:
            a = int(data[0])
            x, y, z = float(data[1]), float(data[2]), float(data[3])
            # print "a ", a, " x ", x, " y ", y, " z ", z
            # points.InsertPoint(a, x, y, z)
            points.InsertPoint(nrAtoms, x, y, z)
        nrAtoms = nrAtoms + 1
        line = file.readline()
        data = line.split()
    return points, nrAtoms


# Read vectors
# We must know the time step and the number of atoms per time


def readVectorsData(file, time, nrAtoms, Nmax):
    # Create a Double array which represents the vectors
    vectors = vtk.vtkFloatArray()
    colors = vtk.vtkFloatArray()

    # Define number of elemnts
    vectors.SetNumberOfComponents(3)
    colors.SetNumberOfComponents(1)
    for i in range(7):
        line = file.readline()

    i = 0
    # Read all data for a certain time
    while i < nrAtoms:
        line = file.readline()
        if i < Nmax:
            data = line.split()
            t, a = int(data[0]), int(data[1])
            x, y, z = float(data[4]), float(data[5]), float(data[6])
            # x, y, z = float(data[3]), float(data[4]), float(data[5])
            # m = (z+1.0)/2.0
            # m = (x+y+2.0)/2.0
            # m = (z+y+x)/3.0**0.5
            m = atan2(x, y) / acos(-1) + 1
            # m = (atan2(x,z)/acos(-1)+1)/2
            # print m
            vectors.InsertTuple3(i, x, y, z)
            colors.InsertValue(i, m)
            i = i + 1
    return vectors, colors


renWin = vtk.vtkRenderWindow()
# win2im=vtk.vtkWindowToImageFilter()
##win2im.SetRenderWindow(renWin)
# win2im.SetInput(renWin)
# Create the RenderWindow,Renderer and Interactor
# ren = vtk.vtkRenderer()
# renWin.AddRenderer(ren)
# ren = vtk.vtkRenderer()
ren = vtk.vtkOpenGLRenderer()
renWin.AddRenderer(ren)


# Set color of backgroung
# ren.SetBackground(0.0, 0.0, 0.0)
ren.SetBackground(1.0, 1.0, 1.0)
renWin.SetSize(1048, 640)
# renWin.SetSize(960,960)


# Datatest=vtk.vtkUnstructuredGrid()
Datatest = vtk.vtkPolyData()
DataScal = vtk.vtkUnstructuredGrid()


# ----------------------------
# Screenshot code begins here
# ----------------------------
number_of_screenshots = 1
# A function that takes a renderwindow and saves its contents to a .png file


def Screenshot(renWin):
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

    number_of_screenshots = number_of_screenshots + 1
    return


# ---------------------------
# Screenshot code ends here
# ---------------------------

# ---------------------------
# Timer code starts here
# ---------------------------


class vtkTimerCallback:
    def __init__(self):
        self.timer_count = 0

    def execute(self, obj, event):
        vecz, colz = readVectorsData(directionsFile, self.timer_count, nrAtoms, Nmax)
        # print "Set data.."
        Datatest.GetPointData().SetVectors(vecz)
        Datatest.GetPointData().SetScalars(colz)
        DataScal.GetPointData().SetScalars(colz)
        iren = obj
        # iren.GetRenderWindow().Render()
        Screenshot(iren.GetRenderWindow())
        self.timer_count += 1


# Size of system
Nmax = 1000000

# Viewing distance
# decrease number to zoom out and vice versa
dollyA = 0.014
dollyB = 0.000
# Open files
momfiles = glob.glob("restart.*.out")
directionsFile = open(momfiles[0])
# directionsFile = open("momentsparsed.out")
posfiles = glob.glob("coord.*.out")
# print posfiles
atomsFile = open(posfiles[0])
# atomsFile = open("atomsparsed.out")

# Read atom positions
atomData, nrAtoms = readAtoms(atomsFile, Nmax)
# print nrAtoms
Datatest.SetPoints(atomData)
DataScal.SetPoints(atomData)

# Read data for vectors
# for i in range(0,55,1):
vecz, colz = readVectorsData(directionsFile, 0, nrAtoms, Nmax)
Datatest.GetPointData().SetVectors(vecz)
Datatest.GetPointData().SetScalars(colz)
DataScal.GetPointData().SetScalars(colz)

# Create colortable for the coloring of the vectors
lut = vtk.vtkLookupTable()
# for i in range(0,255,1):
#    lut.SetTableValue(i,i/255.0,0.6*i/255.0,0,1)
##    #lut.SetTableValue(i,255.0*(1.0-i/255.0),00.0*(1.0-i/255.0),0,1)
##    #lut.SetTableValue(i,255.0*(1.0-i/255.0),200.0*(1.0-i/255.0),0,1)
for i in range(0, 128, 1):
    lut.SetTableValue(i, (127.0 - i) / 127.0, i / 127.0, 0, 1)
for i in range(128, 256, 1):
    lut.SetTableValue(i, 0, (256.0 - i) / 128.0, (i - 128.0) / 128, 1)
# lut.SetTableRange(0.0,1.0);
lut.SetTableRange(-1.0, 1.0)
lut.Build()


# Set up atoms
ball = vtk.vtkSphereSource()
ball.SetRadius(1.00)
ball.SetThetaResolution(16)
ball.SetPhiResolution(16)

balls = vtk.vtkGlyph3DMapper()
# balls.SetInputConnection(Datatest.GetProducerPort())
balls.SetInputData(Datatest)
balls.SetSourceConnection(ball.GetOutputPort())
balls.SetScaleFactor(0.75)
# balls.SetVectorModeToUseVector()
# balls.SetColorModeToColorByVector()
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
# atom.GetProperty().SetColor(0.2,0.2,0.2)
# atom.GetProperty().SetSpecularColor(1.0, 0.2, 0.2)
atom.GetProperty().SetSpecular(0.3)
atom.GetProperty().SetSpecularPower(60)
atom.GetProperty().SetAmbient(0.2)
atom.GetProperty().SetDiffuse(0.8)

# Create vectors
arrow = vtk.vtkArrowSource()
arrow.SetTipRadius(0.25)
arrow.SetShaftRadius(0.15)
arrow.SetTipResolution(24)
arrow.SetShaftResolution(24)

glyph = vtk.vtkGlyph3DMapper()
glyph.SetSourceConnection(arrow.GetOutputPort())
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

vector.GetProperty().SetSpecular(0.3)
vector.GetProperty().SetSpecularPower(60)
vector.GetProperty().SetAmbient(0.2)
vector.GetProperty().SetDiffuse(0.8)
vector.GetProperty().SetOpacity(1.0)

# Scalar bar
scalarBar = vtk.vtkScalarBarActor()
scalarBar.SetLookupTable(lut)
scalarBar.SetOrientationToHorizontal()
scalarBar.SetNumberOfLabels(0)
scalarBar.SetPosition(0.1, 0.05)
scalarBar.SetWidth(0.85)
scalarBar.SetHeight(0.3)
scalarBar.GetLabelTextProperty().SetFontSize(8)

# Depth sorted field
dsDatatest = vtk.vtkDepthSortPolyData()
# dsDatatest.SetInputConnection(Datatest.GetProducerPort())
dsDatatest.SetInputData(Datatest)
dsDatatest.SetCamera(ren.GetActiveCamera())
dsDatatest.SortScalarsOn()
dsDatatest.Update()

# cubes
cube = vtk.vtkCubeSource()
cubes = vtk.vtkGlyph3DMapper()
# cubes.SetInputConnection(Datatest.GetProducerPort())
cubes.SetInputConnection(dsDatatest.GetOutputPort())
cubes.SetSourceConnection(cube.GetOutputPort())
cubes.SetScaleModeToNoDataScaling()
cubes.SetScaleFactor(0.995)
cubes.SetLookupTable(lut)
cubes.SetColorModeToMapScalars()
cubes.ScalarVisibilityOn()
cubes.OrientOff()
cubes.Update()
cubes.SetLookupTable(lut)

cubeActor = vtk.vtkActor()
cubeActor.SetMapper(cubes)
cubeActor.GetProperty().SetOpacity(0.05)
cubeActor.SetPosition(-xmid, -ymid, -zmid)

# hedgehog
hhog = vtk.vtkHedgeHog()
# hhog.SetInputConnection(Datatest.GetProducerPort())
hhog.SetInputData(Datatest)
hhog.SetScaleFactor(5.0)
hhogMapper = vtk.vtkPolyDataMapper()
hhogMapper.SetInputConnection(hhog.GetOutputPort())
hhogMapper.SetLookupTable(lut)
hhogMapper.ScalarVisibilityOn()
hhogActor = vtk.vtkActor()
hhogActor.SetMapper(hhogMapper)
hhogActor.SetPosition(-xmid, -ymid, -zmid)

# cut plane
plane = vtk.vtkPlane()
plane.SetOrigin(DataScal.GetCenter())
plane.SetNormal(0.0, 0.0, 1.0)
planeCut = vtk.vtkCutter()
# planeCut.SetInputConnection(DataScal.GetProducerPort())
planeCut.SetInputData(DataScal)
planeCut.SetInputArrayToProcess(
    0,
    0,
    0,
    vtk.vtkDataObject().FIELD_ASSOCIATION_POINTS,
    vtk.vtkDataSetAttributes().SCALARS,
)
planeCut.SetCutFunction(plane)
planeCut.GenerateCutScalarsOff()
planeCut.SetSortByToSortByCell()
clut = vtk.vtkLookupTable()
clut.SetHueRange(0, 0.67)
clut.Build()
planeMapper = vtk.vtkPolyDataMapper()
planeMapper.SetInputConnection(planeCut.GetOutputPort())
planeMapper.ScalarVisibilityOn()
planeMapper.SetScalarRange(DataScal.GetScalarRange())
# print DataScal.GetScalarRange()
# print planeCut
planeMapper.SetLookupTable(clut)
planeActor = vtk.vtkActor()
planeActor.SetMapper(planeMapper)
# print planeMapper

# clip plane
planeClip = vtk.vtkClipDataSet()
planeClip.SetInputData(DataScal)
planeClip.SetInputArrayToProcess(
    0,
    0,
    0,
    vtk.vtkDataObject().FIELD_ASSOCIATION_POINTS,
    vtk.vtkDataSetAttributes().SCALARS,
)
# planeClip.SetInputConnection(DataScal.GetProducerPort())
planeClip.SetClipFunction(plane)
# print planeClip
planeClip.InsideOutOn()
# planeClip.GenerateCutScalarsOff()
# planeClip.SetSortByToSortByCell()
clipMapper = vtk.vtkDataSetMapper()
clipMapper.SetInputConnection(planeCut.GetOutputPort())
clipMapper.ScalarVisibilityOn()
clipMapper.SetScalarRange(DataScal.GetScalarRange())
clipMapper.SetLookupTable(clut)
# print clipMapper
clipActor = vtk.vtkActor()
clipActor.SetMapper(clipMapper)
# print planeMapper

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
# create a text actor
txt = vtk.vtkTextActor()
txt.SetInput("T = TEMP K")
txtprop = txt.GetTextProperty()
txtprop.SetFontFamilyToArial()
txtprop.SetFontSize(36)
txtprop.SetColor(0, 0, 0)
txt.SetDisplayPosition(30, 550)

# Reposition the camera
# z-direction
ren.GetActiveCamera().Azimuth(0)
ren.GetActiveCamera().Elevation(0)
# xy-direction (almost)
# ren.GetActiveCamera().Azimuth(105)
# ren.GetActiveCamera().Roll(90)
# ren.GetActiveCamera().Dolly(dollyA)
ren.GetActiveCamera().ParallelProjectionOn()
# ren.GetActiveCamera().ParallelProjectionOff()
d = ren.GetActiveCamera().GetDistance()
# ren.GetActiveCamera().SetParallelScale(0.55*max(zmax,xmax))
ren.GetActiveCamera().SetFocalPoint(0, 0, 0)
# ren.GetActiveCamera().SetViewUp(-0.866025403784439,0.5,0)
ren.GetActiveCamera().SetViewUp(0, 1, 0)
ren.GetActiveCamera().SetParallelScale(0.60 * ymax)
# ren.GetActiveCamera().SetPosition(0,0,-100*d)
# print ren.GetActiveCamera().GetViewAngle()
l = max(xmax - xmin, zmax - zmin) / 2
h = l / 0.26795 * 1.1
# print l,h

ren.GetActiveCamera().SetPosition(0, 0, h)
# ren.GetActiveCamera().SetPosition(0,0,50*d)

# For the interactive control. Set up a check for aborting rendering.


def CheckAbort(obj, event):
    # obj will be the object generating the event.  In this case it
    # is renWin.
    if obj.GetEventPending() != 0:
        obj.SetAbortRender(1)


renWin.AddObserver("AbortCheckEvent", CheckAbort)

# Add the actors to the renderer, set the background and size
# Atoms
# ren.AddActor(atom)
# ren.AddActor(txt)
# Vectors
ren.AddActor(vector)
# Cubes
# ren.AddActor(cubeActor)
# HedgeHog
# ren.AddActor(hhogActor)
# CutPlane
# ren.AddActor(planeActor)
# ClipPlane
# ren.AddActor(clipActor)
# Outline
# ren.AddActor(outline)
# Scalar bar
# ren.AddActor(scalarBar)

# Render scene
iren = vtk.vtkRenderWindowInteractor()
istyle = vtk.vtkInteractorStyleTrackballCamera()
iren.SetInteractorStyle(istyle)
iren.SetRenderWindow(renWin)
# ren.ResetCamera()
iren.Initialize()
#
#
cb = vtkTimerCallback()
###cb.AddActor = vector
# iren.AddObserver('TimerEvent', cb.execute)
# timerId = iren.CreateRepeatingTimer(100);
# iren.SetStillUpdateRate(0.050)
# iren.SetDesiredUpdateRate(0.050)
renWin.Render()
iren.Start()
Screenshot(renWin)
