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

import vtk
from ReadTimeData import * # Functions that let us read the files

renWin = vtk.vtkRenderWindow()

Datatest=vtk.vtkUnstructuredGrid()
# Size of system
Nmax = NMAX

# Viewing distance
# decrease number to zoom out and vice versa
dollyA=0.018
dollyB=0.000
# Open files
directionsFile = open("momentsparsed.out")
atomsFile = open("atomsparsed.out")

# Read atom positions
atomData,nrAtoms=(readAtoms(atomsFile,Nmax))
Datatest.SetPoints(atomData)

#----------------------------
# Screenshot code begins here
#----------------------------
number_of_screenshots=1
#A function that takes a renderwindow and saves its contents to a .png file
def Screenshot(renWin):
    global number_of_screenshots
    win2im=vtk.vtkWindowToImageFilter()
    win2im.SetInput(renWin)
    toPNG=vtk.vtkPNGWriter()
    toPNG.SetFileName('snap%.3d.png' %number_of_screenshots)
    toPNG.SetInput(win2im.GetOutput())
    toPNG.Write()

    number_of_screenshots=number_of_screenshots+1
    return;
#---------------------------
# Screenshot code ends here
#---------------------------

# Set up atoms
ball = vtk.vtkSphereSource()
ball.SetRadius(0.15)
ball.SetThetaResolution(8)
ball.SetPhiResolution(8)

balls = vtk.vtkGlyph3D()
balls.SetInput(Datatest)
balls.SetSource(ball.GetOutput())
balls.SetScaleFactor(1.5)

ballMapper = vtk.vtkPolyDataMapper()
ballMapper.SetInput(balls.GetOutput())

atom = vtk.vtkActor()
atom.SetMapper(ballMapper)
xmin,xmax = atom.GetXRange()
ymin,ymax = atom.GetYRange()
zmin,zmax = atom.GetZRange()
xmid = (xmin+xmax)/2
ymid = (ymin+ymax)/2
zmid = (zmin+zmax)/2
#Good position for "z-direction"
atom.SetPosition(-xmid,-ymid,-zmid)
#Good position for "xy-direction"
#atom.SetPosition(-xmid*1.4,-ymid,-zmid)

# Set the color
atom.GetProperty().SetColor(0.2,0.2,0.2)
atom.GetProperty().SetSpecularColor(1.0, 0.2, 0.2)
atom.GetProperty().SetSpecular(0.3)
atom.GetProperty().SetSpecularPower(60)
atom.GetProperty().SetAmbient(0.2)
atom.GetProperty().SetDiffuse(0.8)

# Create vectors
arrow = vtk.vtkArrowSource()
arrow.SetTipRadius(0.20)
arrow.SetShaftRadius(0.15)
arrow.SetTipResolution(40)
arrow.SetShaftResolution(40)

# Read data for vectors
vecz=(readVectorsData(directionsFile,1,nrAtoms,Nmax))
Datatest.GetPointData().SetVectors(vecz)

glyph = vtk.vtkGlyph3D()
glyph.SetSource(arrow.GetOutput())
glyph.SetInput(Datatest) # Position and direction
glyph.SetScaleFactor(0.99)
# Color the vectors according to magnitude
glyph.SetVectorModeToUseVector()
glyph.SetColorModeToColorByVector()
glyph.Update()

# Create colortable for the coloring of the vectors
lut = vtk.vtkLookupTable()
for i in range(0,255,1):
    lut.SetTableValue(i,(i/255.0),0.25-(0.5-i/255.0)*(0.5-i/255.0),(1-i/255.0),1)
lut.SetTableRange(0.0, 1.00000);
lut.Build()

glyphMapper = vtk.vtkPolyDataMapper()
glyphMapper.SetInput(glyph.GetOutput())
glyphMapper.SetLookupTable(lut)

vector = vtk.vtkActor()
vector.SetMapper(glyphMapper)
#Good position for "z-direction"
vector.SetPosition(-xmid,-ymid,-zmid)
#Good position for "xy-direction"
#vector.SetPosition(-xmid*1.4,-ymid,-zmid)

vector.GetProperty().SetSpecular(0.3)
vector.GetProperty().SetSpecularPower(60)
vector.GetProperty().SetAmbient(0.2)
vector.GetProperty().SetDiffuse(0.8)

# Create the RenderWindow,Renderer and Interactor
ren = vtk.vtkRenderer()
renWin.AddRenderer(ren)
iren = vtk.vtkRenderWindowInteractor()
iren.SetRenderWindow(renWin)

# Add the actors to the renderer, set the background and size
# Atoms
#ren.AddActor(atom)
# Vectors
ren.AddActor(vector)

# Set color of backgroung 
ren.SetBackground(0.1, 0.1, 0.1)
#renWin.SetSize(640, 480)
renWin.SetSize(1280, 960)

# Reposition the camera
# z-direction
ren.GetActiveCamera().Azimuth(0)
ren.GetActiveCamera().Elevation(0)
# xy-direction (almost)
#ren.GetActiveCamera().Azimuth(105)
#ren.GetActiveCamera().Roll(90)
ren.GetActiveCamera().Dolly(dollyA)

# For the interactive control. Set up a check for aborting rendering.
def CheckAbort(obj, event):
    # obj will be the object generating the event.  In this case it
    # is renWin.    
    if obj.GetEventPending() != 0:
        obj.SetAbortRender(1)
 
renWin.AddObserver("AbortCheckEvent", CheckAbort)

# Render scene
# uncomment line below and to make more than one frame,
#for frame in range(1,999,1):
#    renWin.Render()
#    Screenshot(renWin)
#    vecz=(readVectorsData(directionsFile,frame,nrAtoms,Nmax))
#    Datatest.GetPointData().SetVectors(vecz)
# then the file momentsparsed.out needs to contain several 
# snapshots as well
# Single snapshot command
iren.Initialize()
renWin.Render()
iren.Start()
Screenshot(renWin)
# 



    
    
