#!/usr/bin/env python3
# Script for creating a snapshot from SD-data
# Written by Anders Bergman, after template from Anders Hast
#
# This particular script is suited for an x-y plane of atoms/moments
# Coloring of magnetic moments is determined by their z-components
# to change this, modify the value m as read in ReadTimeData.py
#
# Camera positioning can be changed using GetActiveCamera.Elevation, Roll, and Azimuth

import vtk
import numpy as np

# from scipy.ndimage import gaussian_filter
from vtk.util import numpy_support
from PyQt6.QtWidgets import QFileDialog


class InteractiveASD:
    """
    Class to handle rendering of simulation data from uppasd using vtk inside
    ASD_GUI.

    Inputs:
                    ren     :   vtkOpenGLRenderer()
                    renWin  :   QVTKRenderWindowInteractor().GetRenderWindow()
                    iren    :   QVTKRenderWindowInteractor().GetRenderWindow().GetInteractor()

    Author: Anders Bergman, after template from Anders Hast. Modified by Erik Karpelin.

    """

    def __init__(self, ren, renWin, iren, ASDsim=None):
        self.ren = ren
        self.renWin = renWin
        self.iren = iren
        self.number_of_screenshots = 0
        self.asd = ASDsim
        self.film = False

    def Launch(self):
        """
        Setup function to add all needed Actors and run uppasd setup. The function also
        iclude a simple keyboard interface.

        Todo:   Clean up function, remove unnecessary comments and functionality, such as
                        the keyboard inputs.
        """

        print("InteractiveASD launched! ----------------------")
        if self.asd is None:
            return
            # try:
            #     self.asd = sim.Simulator()
            #     print("ASDsimulation initialized.")
            # except ImportError:
            #     print("Launch: UppASD module not installed.")
            #     return
        else:
            self.asd.init_simulation()
        print("InteractiveASD launched!", self.asd.natom, self.asd.mensemble)

        self.Datatest = vtk.vtkPolyData()
        self.DataFour = vtk.vtkPolyData()

        self.initmom = self.asd.moments[:, :, 0].T
        self.currmom = self.asd.moments[:, :, 0].T
        self.vecz = numpy_support.numpy_to_vtk(self.currmom, deep=False)
        self.initcol = self.asd.moments[2, :, 0].T
        self.currcol = self.asd.moments[2, :, 0].T
        self.colz = numpy_support.numpy_to_vtk(self.currcol, deep=False)

        self.lut = vtk.vtkLookupTable()
        for i in range(0, 128, 1):
            self.lut.SetTableValue(i, (127.0 - i) / 127.0, i / 127.0, 0, 1)
        for i in range(128, 256, 1):
            self.lut.SetTableValue(i, 0, (256.0 - i) / 128.0, (i - 128.0) / 128, 1)
        self.lut.SetTableRange(-1.0, 1.0)
        self.lut.Build()

        # Size of system
        Nmax = 1000000

        # Open files
        # momfiles = glob.glob("restart.????????.out")
        # directionsFile = open(momfiles[0])
        # posfiles = glob.glob("coord.????????.out")
        # atomsFile = open(posfiles[0], encoding='utf-8')

        # Read atom positions
        atomPoints = vtk.vtkPoints()
        atomData = numpy_support.numpy_to_vtk(self.asd.coords.T, deep=False)
        atomPoints.SetData(atomData)
        nrAtoms = self.asd.natom
        # atomData, nrAtoms = self.readAtoms(atomsFile)
        print("Coordinates read: ", nrAtoms )
        print("Coordinates: ", self.asd.coords)
        print("Number of atoms: ", nrAtoms)
        self.Datatest.SetPoints(atomPoints)
        self.DataFour.SetPoints(atomPoints)

        # Read data for vectors
        self.Datatest.GetPointData().SetVectors(self.vecz)
        self.Datatest.GetPointData().SetScalars(self.colz)
        self.DataFour.GetPointData().SetVectors(self.vecz)
        self.DataFour.GetPointData().SetScalars(self.colz)

        # Create colortable for the coloring of the vectors

        # ATOMS  ##########################################
        # Set up atoms
        ball = vtk.vtkSphereSource()
        ball.SetRadius(1.00)
        ball.SetThetaResolution(16)
        ball.SetPhiResolution(16)

        balls = vtk.vtkGlyph3DMapper()
        balls.SetInputData(self.Datatest)
        balls.SetSourceConnection(ball.GetOutputPort())
        balls.SetScaleFactor(0.57)
        balls.SetScaleModeToNoDataScaling()
        balls.SetLookupTable(self.lut)
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
        atom.SetPosition(-xmid, -ymid, -zmid)

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
        fball.SetThetaResolution(16)
        fball.SetPhiResolution(16)

        fballs = vtk.vtkGlyph3DMapper()
        fballs.SetInputData(self.DataFour)
        fballs.SetSourceConnection(fball.GetOutputPort())
        fballs.SetScaleFactor(0.57)
        fballs.SetScaleModeToNoDataScaling()
        fballs.SetLookupTable(self.lut)
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
        four.SetPosition(-xmid, -ymid, -zmid)

        # Set the color
        four.GetProperty().SetSpecular(0.3)
        four.GetProperty().SetSpecularPower(60)
        four.GetProperty().SetAmbient(0.2)
        four.GetProperty().SetDiffuse(0.8)

        # ARROWS ##########################################
        # Create vectors
        arrow = vtk.vtkArrowSource()
        arrow.SetTipRadius(0.16)
        arrow.SetShaftRadius(0.09)
        arrow.SetTipResolution(24)
        arrow.SetShaftResolution(24)

        glyph = vtk.vtkGlyph3DMapper()
        glyph.SetSourceConnection(arrow.GetOutputPort())
        glyph.SetInputData(self.Datatest)
        # glyph.SetInput(Datatest) # Position and direction
        glyph.SetScaleFactor(1.00)
        # Color the vectors according to magnitude
        # glyph.SetVectorModeToUseVector()
        # glyph.SetColorModeToColorByVector()
        glyph.SetScaleModeToNoDataScaling()
        # glyph.Update()

        glyph.SetLookupTable(self.lut)
        glyph.SetColorModeToMapScalars()
        glyph.Update()

        vector = vtk.vtkLODActor()
        vector.SetMapper(glyph)
        vector.SetPosition(-xmid, -ymid, -zmid)
        vector.GetProperty().SetAmbient(0.3)
        vector.GetProperty().SetDiffuse(0.5)
        vector.GetProperty().SetOpacity(1.0)
        vector.GetProperty().SetRoughness(0.6)
        vector.GetProperty().SetMetallic(0.4)
        vector.GetProperty().ShadingOn()
        print("Interpolation: ", vector.GetProperty().GetInterpolationAsString())

        # Bounding box
        outlineData = vtk.vtkOutlineFilter()
        # outlineData.SetInputConnection(Datatest.GetProducerPort())
        outlineData.SetInputData(self.Datatest)
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
        self.temptxt = vtk.vtkTextActor()
        # temp = f"{0.0:4.3f}"
        temp='{:4.3f}'.format(self.asd.inputdata.get_temp())
        self.temptxt.SetInput("T = " + temp + " K")
        temptxtprop = self.temptxt.GetTextProperty()
        temptxtprop.SetFontFamilyToArial()
        temptxtprop.SetFontSize(30)
        temptxtprop.SetColor(0, 0, 0)
        temptxtprop.BoldOn()
        self.temptxt.SetDisplayPosition(20, 1350)

        # create a text actor for Field
        self.fieldtxt = vtk.vtkTextActor()
        Bfield = self.asd.inputdata.get_hfield()
        self.fieldtxt.SetInput(f"B = ({Bfield[0]:4.1f}, {Bfield[1]:4.1f}, {Bfield[2]:4.1f} ) T")
        # self.fieldtxt.SetInput(f"B = ({0:4.1f}, {1:4.1f}, {2:4.1f} ) T")
        fieldtxtprop = self.fieldtxt.GetTextProperty()
        fieldtxtprop.SetFontFamilyToArial()
        fieldtxtprop.SetFontSize(30)
        fieldtxtprop.SetColor(0, 0, 0)
        fieldtxtprop.BoldOn()
        self.fieldtxt.SetDisplayPosition(20, 1300)

        # create a text actor for Energy
        self.enetxt = vtk.vtkTextActor()
        self.asd.calculate_energy()
        ene = f"{self.asd.energy:6.4f}"
        self.enetxt.SetInput(f"E = {ene} mRy/atom")
        enetxtprop = self.enetxt.GetTextProperty()
        enetxtprop.SetFontFamilyToArial()
        enetxtprop.SetFontSize(20)
        enetxtprop.SetColor(0, 0, 0)
        enetxtprop.BoldOn()
        self.enetxt.SetDisplayPosition(20, 50)

        # LIGHTS ON
        light = vtk.vtkLight()
        light.SetColor(1.0, 1.0, 1.0)
        self.ren.AddLight(light)

        # Reposition the camera
        # z-direction
        self.ren.GetActiveCamera().Azimuth(0)
        self.ren.GetActiveCamera().Elevation(0)
        self.ren.GetActiveCamera().ParallelProjectionOn()
        # d = self.ren.GetActiveCamera().GetDistance()
        self.ren.GetActiveCamera().SetFocalPoint(0, 0, 0)
        self.ren.GetActiveCamera().SetViewUp(-0.866025403784439, 0.5, 0)
        self.ren.GetActiveCamera().SetParallelScale(0.55 * ymax)
        l_dist = max(xmax - xmin, zmax - zmin) / 2
        h = l_dist / 0.26795 * 1.1

        self.ren.GetActiveCamera().SetPosition(0, 0, h)

        # Add the actors to the renderer, set the background and size
        # Atoms
        atom.SetVisibility(0)
        self.ren.AddActor(atom)
        # Vectors
        self.ren.AddActor(vector)
        # Text
        self.ren.AddActor(self.temptxt)
        self.ren.AddActor(self.fieldtxt)
        self.ren.AddActor(self.enetxt)

        # self.iren.AddObserver("KeyPressEvent", Keypress)

        self.renWin.Render()
        self.iren.Start()

        # For the interactive control. Set up a check for aborting rendering.
        def CheckAbort(obj, event):
            # obj will be the object generating the event.  In this case it
            # is renWin.
            if obj.GetEventPending() != 0:
                obj.SetAbortRender(1)

    def S_Step(self):
        """Do a simulation using S-mode."""
        if not hasattr(self, "asd"):
            return
        self.asd.relax(mode="S", temperature=self.asd.inputdata.temp)
        currmom = self.asd.moments[:, :, 0].T
        currcol = self.asd.moments[2, :, 0].T
        vecz = numpy_support.numpy_to_vtk(currmom)
        colz = numpy_support.numpy_to_vtk(currcol)
        self.Datatest.GetPointData().SetVectors(vecz)
        self.Datatest.GetPointData().SetScalars(colz)
        # Update enegrgy
        self.asd.calculate_energy()
        ene = f"{self.asd.energy:6.4f}"
        self.enetxt.SetInput(f"E = {ene} mRy/atom")
        self.renWin.Render()
        if self.film:
            self.Screenshot()

    def M_step(self):
        """Do a simulation using Metropolis MC"""
        if not hasattr(self, "asd"):
            return
        self.asd.relax(mode="M", temperature=self.asd.inputdata.temp)
        currmom = self.asd.moments[:, :, 0].T
        currcol = self.asd.moments[2, :, 0].T
        vecz = numpy_support.numpy_to_vtk(currmom)
        colz = numpy_support.numpy_to_vtk(currcol)
        self.Datatest.GetPointData().SetVectors(vecz)
        self.Datatest.GetPointData().SetScalars(colz)
        # Update enegrgy
        self.asd.calculate_energy()
        ene = f"{self.asd.energy:6.4f}"
        self.enetxt.SetInput(f"E = {ene} mRy/atom")
        self.renWin.Render()
        if self.film:
            self.Screenshot()

    def H_step(self):
        """Do a simulation using Heat-bath MC"""
        if not hasattr(self, "asd"):
            return
        self.asd.relax(mode="H", temperature=self.asd.inputdata.temp+1.0e-6)
        currmom = self.asd.moments[:, :, 0].T
        currcol = self.asd.moments[2, :, 0].T
        vecz = numpy_support.numpy_to_vtk(currmom)
        colz = numpy_support.numpy_to_vtk(currcol)
        self.Datatest.GetPointData().SetVectors(vecz)
        self.Datatest.GetPointData().SetScalars(colz)
        # Update enegrgy
        self.asd.calculate_energy()
        ene = f"{self.asd.energy:6.4f}"
        self.enetxt.SetInput(f"E = {ene} mRy/atom")
        self.renWin.Render()
        if self.film:
            self.Screenshot()

    def Reset(self):
        """Reset data to initial."""
        if not hasattr(self, "asd"):
            return
        self.asd.put_moments(self.initmom.T)
        vecz = numpy_support.numpy_to_vtk(self.initmom)
        colz = numpy_support.numpy_to_vtk(self.initcol)
        self.Datatest.GetPointData().SetVectors(vecz)
        self.Datatest.GetPointData().SetScalars(colz)
        self.renWin.Render()

    def read_moments(self):

        dlg = QFileDialog()
        dlg.setFileMode(QFileDialog.FileMode.ExistingFile)
        dlg.setDirectory(".")
        if dlg.exec():
            mag_file = dlg.selectedFiles()[0]
            print('Reading moments from:', mag_file)
            try:
                moments = np.loadtxt(mag_file)[:, 4:]
                print('Moments read:', moments.shape)
                self.asd.put_moments(moments.T)
                self.currmom = moments
                self.vecz = numpy_support.numpy_to_vtk(self.currmom, deep=False)
                self.currcol = moments[:,2]
                self.colz = numpy_support.numpy_to_vtk(self.currcol, deep=False)
                self.Datatest.GetPointData().SetVectors(self.vecz)
                self.Datatest.GetPointData().SetScalars(self.colz)
                self.renWin.Render()
                print('Moments read and updated')
            except FileNotFoundError:
                print('Error reading file')
        else:
            print('No file selected')
        return

    def UpdateTemperature(self):
        """Update temperature actor."""
        if not hasattr(self, "asd"):
            return
        temp = f"{self.asd.inputdata.get_temp():4.3f}"
        self.temptxt.SetInput("T = " + temp + " K")
        self.renWin.Render()

    def UpdateBfield(self):
        """Update B-field actor."""
        if not hasattr(self, "asd"):
            return
        Bfield = self.asd.inputdata.get_hfield()
        self.fieldtxt.SetInput(
            f"B = ({Bfield[0]:4.1f}, {Bfield[1]:4.1f}, {Bfield[2]:4.1f} ) T"
        )
        self.renWin.Render()

    def close_window(self):
        if not hasattr(self, "asd"):
            return
        render_window = self.iren.GetRenderWindow()
        render_window.Finalize()
        self.iren.TerminateApp()

    # Read Location of Atoms
    def readAtoms(self, file):
        print('readAtoms called')
        points = vtk.vtkPoints()
        nrAtoms = 0
        # Read ahead
        line = file.readline()
        data = line.split()

        # Read all data
        while line:
            _ = int(data[0])
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
    def readVectorsData(file, time, nrAtoms):
        print('readVectorsData called')
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
            data = line.split()
            t, a = int(data[0]), int(data[1])
            if a == 1:
                x, y, z = float(data[4]), float(data[5]), float(data[6])
                # x, y, z = float(data[3]), float(data[4]), float(data[5])
                m = (z + 1.0) / 2.0
                # m = (x+y+2.0)/2.0
                # m = (z+y+x)/3.0**0.5
                # m = atan2(x,y)/acos(-1)+1
                # m = (atan2(x,z)/acos(-1)+1)/2
                # print m
                vectors.InsertTuple3(i, x, y, z)
                colors.InsertValue(i, m)
                i = i + 1
        print("Vectors read: ", i)
        return vectors, colors

        # A function that takes a renderwindow and saves its contents to a .png file

    def Screenshot(self):
        self.number_of_screenshots
        win2im = vtk.vtkWindowToImageFilter()
        win2im.ReadFrontBufferOff()
        win2im.SetInput(self.renWin)
        #
        toPNG = vtk.vtkPNGWriter()
        toPNG.SetFileName(f"snap{self.number_of_screenshots:05d}.png")
        # toPNG.SetInput(win2im.GetOutput())
        toPNG.SetInputConnection(win2im.GetOutputPort())
        toPNG.Write()

        self.number_of_screenshots += 1
        return

    def UpdateTextPlacement(self):
        pass
