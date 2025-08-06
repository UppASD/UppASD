###! /usr/bin/env vtk6python
# (c) Anders.Bergman@physics.uu.se
#
#

# Import the modules for the code
import glob
import string
import sys
import tkinter

import vtk

# import vtkTkRenderWindowInteractor
from vtk.tk.vtkTkRenderWindowInteractor import vtkTkRenderWindowInteractor


# Read Location of Atoms
def readAtoms(file):

    points = vtk.vtkPoints()
    coords = []
    atypes = vtk.vtkFloatArray()
    atypes.SetNumberOfComponents(1)
    nrAtoms = 0
    # Read ahead
    line = file.readline()
    data = line.split()

    # Read all data
    while line:
        a = int(data[0])
        x, y, z = float(data[1]), float(data[2]), float(data[3])
        at = float(data[5]) / 1
        cvec = []
        cvec.append(x)
        cvec.append(y)
        cvec.append(z)
        atypes.InsertValue(nrAtoms, at)
        points.InsertPoint(nrAtoms, x, y, z)
        nrAtoms = nrAtoms + 1
        line = file.readline()
        data = line.split()
        coords.append(cvec)
    return points, atypes, nrAtoms, coords


# Read vectors
# We must know the time step and the number of atoms per time
def readBLSData(file, nfreq, nrAtoms):
    # Create a Double array which represents the vectors
    freq_arr = []
    frequencies = vtk.vtkFloatArray()
    tempa = []
    freqlist = []

    minval = 100
    maxval = -100
    # Define number of elemnts
    frequencies.SetNumberOfComponents(nfreq)
    # frequencies.SetNumberOfComponents(nfreq)
    # line = file.readline()
    freq = 0
    # Read all data for a certain time
    while freq < nfreq:
        iat = 0
        tempa = []
        while iat < nrAtoms:
            line = file.readline()
            data = line.split()
            fr = float(data[2])
            x, y, z, tot = (
                float(data[3]),
                float(data[4]),
                float(data[5]),
                float(data[6]),
            )
            # m = (x*x+y*y)**0.5
            m = tot
            # m = (z*z+y*y)**0.5
            # m = x
            minval = min(m, minval)
            maxval = max(m, maxval)
            tempa.append(m)
            frequencies.InsertComponent(iat, freq, m)
            iat = iat + 1
        # frequencies.append(tempa)
        freq_arr.append(tempa)
        # frequencies.InsertNextTuple(tempa)
        # frequencies.InsertTuple(freq,tempa)
        freqlist.append(fr)
        freq = freq + 1
    return frequencies, freqlist, freq_arr, minval, maxval


# A function that takes a renderwindow and saves its contents to a .png file
def snapshot():
    win2im = vtk.vtkWindowToImageFilter()
    win2im.SetInput(renWindow)
    toPNG = vtk.vtkPNGWriter()
    toPNG.SetFileName("neighbourmap.png")
    # toPNG.SetInput(win2im.GetOutput())
    toPNG.SetInputConnection(win2im.GetOutputPort())
    toPNG.Write()
    return


# Start of main program

iAtom = 0
iFreq = 0.00
# Open files
# Coordinates
posfiles = glob.glob("coord.*.out")
atomsFile = open(posfiles[0])
# BLS file
blsfiles = glob.glob("bls.*.out")
blsfile = open(blsfiles[0])

# Read coordinates
atoms, atomTypes, nAtoms, coords = readAtoms(atomsFile)
nfreq = 9
# nfreq=10
nrAtoms = nAtoms
print(nfreq, nAtoms, nfreq * nAtoms)
minval = 10000.0
maxval = -10000.0
frequencies, freqlist, freq_arr, gminval, gmaxval = readBLSData(blsfile, nfreq, nrAtoms)
# print len(freq_arr)
print(gminval, gmaxval)
gminval = 0.00
gmaxval = 1.00
print(gminval, gmaxval)


# Update neighbour list according to selected atom
def SetUpdate_c(opc):
    global minval, maxval
    j = int(opc) - 1
    iFreq = float(freqlist[int(opc) - 1])
    AtomCoord["text"] = "[ {0:8.1f} GHz]".format(iFreq / 1000000000.0)
    # curr_freq.CopyComponent(0,frequencies,int(opc)+31)

    minval = 100
    maxval = -100
    # for i in range(nrAtoms):
    #    curr_freq.SetComponent(i,0,0.000)
    # print 'null',curr_freq.GetValueRange()
    # AtomGrid.GetPointData().SetScalars(curr_freq)
    for i in range(nrAtoms):
        # curr_freq.SetComponent(i,0,frequencies.GetComponent(i,j))
        curr_freq.InsertTuple(i, [frequencies.GetComponent(i, j)])
        minval = min(frequencies.GetComponent(i, j), minval)
        maxval = max(frequencies.GetComponent(i, j), maxval)
    print(curr_freq.GetValueRange(), minval, maxval)
    # print curr_freq.GetComponent(999,0)
    AtomGrid.GetPointData().SetScalars(curr_freq)
    ### print balls.GetScalarRange()
    # print balls.GetUseLookupTableScalarRange()
    print(minval, maxval)
    ### #balls.Update()
    # if balls.GetUseLookupTableScalarRange()==0:
    #    lut.SetTableRange(gminval,gmaxval)
    # else:
    #    lut.SetTableRange(minval,maxval)

    # print lut.GetTableRange()
    lut.SetTableRange(minval, maxval)
    # lut.SetTableRange(minval*0.95,maxval*1.05)
    ### #lut.SetTableRange(curr_freq.GetRange())
    lut.Build()
    ### #print opc,lut.GetRange()
    renWindow.Render()


# Toggle the visibility of all atoms
def toggle_Atoms():

    if balls.GetScaleMode() == 1:
        balls.SetScaleModeToNoDataScaling()
    else:
        balls.SetScaleModeToScaleByMagnitude()
    renWindow.Render()


def toggle_outline():
    # print outline.GetProperty().GetOpacity()
    if outline.GetProperty().GetOpacity() == 1.0:
        outline.GetProperty().SetOpacity(0.0)
    else:
        outline.GetProperty().SetOpacity(1.0)
    renWindow.Render()


def toggle_lores():
    if ball.GetThetaResolution() == 12:
        ball.SetThetaResolution(3)
        ball.SetPhiResolution(3)
    else:
        ball.SetThetaResolution(12)
        ball.SetPhiResolution(12)
    balls.Update()
    lut.Build()
    renWindow.Render()


def increase_opacity():
    if atom.GetProperty().GetOpacity() >= 0.67:
        atom.GetProperty().SetOpacity(1.0)
    else:
        atom.GetProperty().SetOpacity(1.5 * atom.GetProperty().GetOpacity())
    renWindow.Render()


def decrease_opacity():
    if atom.GetProperty().GetOpacity() <= 0.015:
        atom.GetProperty().SetOpacity(0.01)
    else:
        atom.GetProperty().SetOpacity(2.0 / 3.0 * atom.GetProperty().GetOpacity())
    renWindow.Render()


# Toggle the visibility of neighbour map
def toggle_scale():
    if balls.GetUseLookupTableScalarRange() == 0:
        balls.UseLookupTableScalarRangeOn()
        lut.SetTableRange(minval, maxval)
    else:
        balls.UseLookupTableScalarRangeOff()
        lut.SetTableRange(gminval, gmaxval)
    lut.Build()
    renWindow.Render()


# Reset position and distance of camera
def reset_cam():
    render.GetActiveCamera().SetPosition(0, 0, 200)
    render.GetActiveCamera().SetParallelScale(0.75 * xmax)
    render.GetActiveCamera().SetFocalPoint(0, 0, 0)
    render.GetActiveCamera().SetViewUp(0, 1, 0)
    renWindow.Render()


# Create colortable for the coloring of the vectors
lut = vtk.vtkLookupTable()
##for i in range(0,255,1):
##    lut.SetTableValue(i,(i/255.0),0.25-(0.5-i/255.0)*(0.5-i/255.0),(1-i/255.0),1)
for i in range(0, 127, 1):
    lut.SetTableValue(i, 0, i / 128.0, (128.0 - i) / 128.0, 1)
for i in range(128, 256, 1):
    lut.SetTableValue(i, (i - 128.0) / 128, (256.0 - i) / 128.0, 0, 1)
### for i in range(0,127,1):
###     lut.SetTableValue(i,(128.0-i)/128.0,i/128.0,0,1)
### for i in range(128,256,1):
###     lut.SetTableValue(i,0,(256.0-i)/128.0,(i-128.0)/128,1)
# lut.SetVectorModeToComponent()
# lut.SetTableRange(frequencies.GetValueRange())
# lut.SetTableRange(minval,maxval)
# lut.SetTableRange(-100,100)
# print frequencies.GetValueRange()
lut.Build()

# Set up atom grid
AtomGrid = vtk.vtkUnstructuredGrid()
AtomGrid.SetPoints(atoms)
# AtomGrid.GetPointData().SetScalars(frequencies)
# Set up atoms
ball = vtk.vtkSphereSource()
ball.SetRadius(1.00)
ball.SetThetaResolution(12)
ball.SetPhiResolution(12)

balls = vtk.vtkGlyph3DMapper()
# balls.SetInputConnection(AtomGrid.GetProducerPort())
balls.SetInputData(AtomGrid)
balls.SetSourceConnection(ball.GetOutputPort())
balls.SetScaleFactor(0.6)
# balls.SetVectorModeToUseVector()
# balls.SetColorModeToColorByVector()
balls.SetScaleModeToNoDataScaling()
balls.SetLookupTable(lut)
balls.ScalarVisibilityOn()
balls.UseLookupTableScalarRangeOn()
balls.SetColorModeToMapScalars()
balls.Update()

atom = vtk.vtkLODActor()
atom.SetMapper(balls)

#
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
atom.GetProperty().SetSpecularColor(1.0, 0.2, 0.2)
atom.GetProperty().SetSpecular(0.3)
atom.GetProperty().SetSpecularPower(60)
atom.GetProperty().SetAmbient(0.2)
atom.GetProperty().SetDiffuse(0.8)

curr_freq = vtk.vtkFloatArray()
curr_freq.SetNumberOfComponents(1)
curr_freq.SetNumberOfTuples(nrAtoms)

# Bounding box
outlineData = vtk.vtkOutlineFilter()
# outlineData.SetInputConnection(AtomGrid.GetProducerPort())
outlineData.SetInputData(AtomGrid)
outlineMapper = vtk.vtkPolyDataMapper()
outlineMapper.SetInputConnection(outlineData.GetOutputPort())
# outlineMapper.SetInput(outlineData.GetOutput())
outline = vtk.vtkActor()
outline.SetMapper(outlineMapper)
outline.GetProperty().SetColor(0, 0, 0)
outline.GetProperty().SetLineWidth(5.0)
outline.SetPosition(-xmid, -ymid, -zmid)
outline.GetProperty().SetOpacity(0.0)


# Setup for root window
root = tkinter.Tk()
root.title("UppASD BLS simulator")

# Create tkinter main frame
frame = tkinter.Frame(root)
frame.pack(fill=tkinter.BOTH, expand=1, side=tkinter.TOP)

# Create slicer for selecting atoms
slice = tkinter.Scale(
    root,
    from_=1,
    to=nfreq,
    orient="horizontal",
    command=SetUpdate_c,
    label="Frequency index",
)
slice.pack(fill="x", expand="false")

# Create and pack buttons
ctrl_buttons = tkinter.Frame(root)
ctrl_buttons.pack(side="bottom", anchor="n", fill="both", expand="false")
atogg_button = tkinter.Button(
    ctrl_buttons, text="Toggle atom scale", command=toggle_Atoms, background="gray"
)
atogg_button.pack(side="left", expand="false", anchor="e", fill="both")
ntogg_button = tkinter.Button(
    ctrl_buttons, text="Toggle color scale", command=toggle_scale, background="gray"
)
ntogg_button.pack(side="left", expand="false", anchor="e", fill="both")
# rtogg_button = tkinter.Button(ctrl_buttons, text="Toggle hi/low res.", command=toggle_lores,background="gray")
# rtogg_button.pack(side="left", expand="false", anchor="e", fill="both")
# otogg_button = tkinter.Button(ctrl_buttons, text="Toggle bounding box", command=toggle_outline,background="gray")
# otogg_button.pack(side="left", expand="false", anchor="e", fill="both")
incop_button = tkinter.Button(
    ctrl_buttons, text="Increase opacity", command=increase_opacity, background="gray"
)
incop_button.pack(side="left", expand="false", anchor="n", fill="both")
decop_button = tkinter.Button(
    ctrl_buttons, text="Decrease opacity", command=decrease_opacity, background="gray"
)
decop_button.pack(side="left", expand="false", anchor="n", fill="both")
cres_button = tkinter.Button(
    ctrl_buttons, text="Reset camera", command=reset_cam, background="gray"
)
cres_button.pack(side="left", expand="false", anchor="e", fill="both")
quit_button = tkinter.Button(ctrl_buttons, text="Quit", command=quit, background="red")
quit_button.pack(side="right", expand="false", anchor="e", fill="both")
snap_button = tkinter.Button(
    ctrl_buttons, text="Snapshot", command=snapshot, background="gray"
)
snap_button.pack(side="right", expand="false", anchor="e", fill="both")


#
# -- Create the text frame
textFrame = tkinter.Frame(root)
tlabel = tkinter.StringVar()
# Write the information about current atom coordinate and neighbour count
AtomLabel = tkinter.Label(
    textFrame, justify=tkinter.LEFT, anchor=tkinter.W, font=("bold")
)
AtomLabel["text"] = "Current frequency:"
AtomCoord = tkinter.Label(textFrame, justify=tkinter.LEFT, font=("bold"))
AtomCoord["text"] = "[ {0:12.8f} ]".format(iFreq)
AtomLabel.pack(side="left")
AtomCoord.pack(side="left")
textFrame.pack()


# Setup for renderer
render = vtk.vtkRenderer()
render.SetBackground(0.95, 0.95, 0.95)
render.ResetCameraClippingRange()
render.AddActor(atom)
render.AddActor(outline)
render.GetActiveCamera().ParallelProjectionOn()
d = render.GetActiveCamera().GetDistance()
render.GetActiveCamera().SetParallelScale(0.75 * xmax)
render.GetActiveCamera().SetFocalPoint(0, 0, 0)
render.GetActiveCamera().SetViewUp(0, 1, 0)
render.GetActiveCamera().SetPosition(0, 0, 200)


# Setup for rendering window
renWindow = vtk.vtkRenderWindow()
renWindow.AddRenderer(render)

# Setup for rendering window interactor
renWinInteract = vtkTkRenderWindowInteractor(root, rw=renWindow, width=600, height=600)
renWinInteract.Initialize()
renWinInteract.pack(side="top", fill="both", expand=1)
renWinInteract.Start()

# Begin execution by updating the renderer and
# starting the tkinter loop
renWindow.Render()

root.mainloop()
