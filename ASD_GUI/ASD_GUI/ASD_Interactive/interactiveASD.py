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

#---------------------------
# Timer code starts here
#---------------------------
# class vtkTimerCallback():
#     def __init__(self):
#         self.timer_count = 0

###     def execute(self,obj,event):
###         vecz,colz=(readVectorsData(directionsFile,self.timer_count,nrAtoms,Nmax))
###         #print "Set data.."
###         Datatest.GetPointData().SetVectors(vecz)
###         Datatest.GetPointData().SetScalars(colz)
###         DataScal.GetPointData().SetScalars(colz)
###         iren = obj
###         #iren.GetRenderWindow().Render()
###         Screenshot(iren.GetRenderWindow())
###         self.timer_count += 1

import vtk
from math import atan2, acos
from copy import copy, deepcopy
import glob
import string
import uppasd as asd
from vtk.util import numpy_support
import numpy as np
from scipy.ndimage import gaussian_filter
import matplotlib.pyplot as plt
import matplotlib.cm as cm

class InteractiveASD:
	"""
	Class to handle rendering of simulation data from uppasd using vtk inside
	ASD_GUI. 

	Inputs:
			ren		:	vtkOpenGLRenderer()
			renWin	:	QVTKRenderWindowInteractor().GetRenderWindow()
			iren	:	QVTKRenderWindowInteractor().GetRenderWindow().GetInteractor()

	Author: nders Bergman, after template from Anders Hast. Modified by Erik Karpelin. 
	
	"""

	def __init__(self, ren, renWin, iren):

		self.ren = ren
		self.renWin = renWin
		self.iren = iren
		self.number_of_screenshots = 0

	def Launch(self):
		"""
		Setup function to add all needed Actors and run uppasd setup. The function also
		iclude a simple keyboard interface.

		Todo:	Clean up function, remove unnecessary comments and functionality, such as 
				the keyboard inputs. 
		"""

		asd.pyasd.setupall()

		self.Datatest=vtk.vtkPolyData()
		self.DataFour=vtk.vtkPolyData()

		self.initmom=deepcopy(asd.momentdata.get_array_emom()[:,:,0].T)
		self.currmom=asd.momentdata.get_array_emom()[:,:,0].T
		self.vecz=numpy_support.numpy_to_vtk(self.currmom,deep=False)
		self.initcol=deepcopy(asd.momentdata.get_array_emom()[2,:,0].T)
		self.currcol=asd.momentdata.get_array_emom()[2,:,0].T
		self.colz=numpy_support.numpy_to_vtk(self.currcol,deep=False)

		self.lut = vtk.vtkLookupTable()
		for i in range(0,128,1):
				self.lut.SetTableValue(i,(127.0-i)/127.0,i/127.0,0,1)
		for i in range(128,256,1):
				self.lut.SetTableValue(i,0,(256.0-i)/128.0,(i-128.0)/128,1)
		self.lut.SetTableRange(-1.0,1.0)
		self.lut.Build()

		# Size of system
		Nmax = 1000000

		# Open files
		momfiles=glob.glob("restart.????????.out")
		directionsFile = open(momfiles[0])
		#directionsFile = open("momentsparsed.out")
		posfiles=glob.glob("coord.????????.out")
		#print posfiles
		atomsFile = open(posfiles[0])
		#atomsFile = open("atomsparsed.out")

		# Read atom positions
		atomData,nrAtoms=(self.readAtoms(atomsFile,Nmax))
		#print nrAtoms
		self.Datatest.SetPoints(atomData)
		#DataScal.SetPoints(atomData)
		self.DataFour.SetPoints(atomData)

		# Read data for vectors
		#for i in range(0,55,1):
		#vecz,colz=(readVectorsData(directionsFile,0,nrAtoms,Nmax))
		self.Datatest.GetPointData().SetVectors(self.vecz)
		self.Datatest.GetPointData().SetScalars(self.colz)
		#DataScal.GetPointData().SetScalars(colz)

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
		xmin,xmax = four.GetXRange()
		ymin,ymax = four.GetYRange()
		zmin,zmax = four.GetZRange()
		xmid = (xmin+xmax)/2
		ymid = (ymin+ymax)/2
		zmid = (zmin+zmax)/2
		#Good position for "z-direction"
		four.SetPosition(-xmid,-ymid,-zmid)
		#Good position for "xy-direction"
		#four.SetPosition(-xmid*1.4,-ymid,-zmid)

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
		#glyph.SetInputConnection(Datatest.GetProducerPort())
		glyph.SetInputData(self.Datatest)
		#glyph.SetInput(Datatest) # Position and direction
		glyph.SetScaleFactor(1.00)
		# Color the vectors according to magnitude
		#glyph.SetVectorModeToUseVector()
		#glyph.SetColorModeToColorByVector()
		glyph.SetScaleModeToNoDataScaling()
		#glyph.Update()

		#glyphMapper = vtk.vtkPolyDataMapper()
		#glyphMapper.SetInput(glyph.GetOutput())
		#glyphMapper.SetLookupTable(lut)
		glyph.SetLookupTable(self.lut)
		#glyph.ColorByArrayComponent(0,0)
		#glyph.UseLookupTableScalarRangeOn()
		glyph.SetColorModeToMapScalars()
		#glyph.SetScalarModeToUseFieldData()
		#glyph.ScalarVisibilityOn()
		glyph.Update()


		vector = vtk.vtkLODActor()
		vector.SetMapper(glyph)
		#Good position for "z-direction"
		vector.SetPosition(-xmid,-ymid,-zmid)
		#Good position for "xy-direction"
		#vector.SetPosition(-xmid*1.4,-ymid,-zmid)

		#vector.GetProperty().SetInterpolationToPBR()
		#vector.GetProperty().SetSpecular(0.0)
		#vector.GetProperty().SetSpecularPower(160)
		vector.GetProperty().SetAmbient(0.3)
		vector.GetProperty().SetDiffuse(0.5)
		vector.GetProperty().SetOpacity(1.0)
		vector.GetProperty().SetRoughness(0.6)
		vector.GetProperty().SetMetallic(0.4)
		vector.GetProperty().ShadingOn()
		print("Interpolation: ",vector.GetProperty().GetInterpolationAsString())



		# Bounding box
		outlineData = vtk.vtkOutlineFilter()
		#outlineData.SetInputConnection(Datatest.GetProducerPort())
		outlineData.SetInputData(self.Datatest)
		outlineMapper = vtk.vtkPolyDataMapper()
		outlineMapper.SetInputConnection(outlineData.GetOutputPort())
		#outlineMapper.SetInput(outlineData.GetOutput())
		outline = vtk.vtkActor()
		outline.SetMapper(outlineMapper)
		outline.GetProperty().SetColor(0, 0, 0)
		outline.GetProperty().SetLineWidth(5.0)
		outline.SetPosition(-xmid,-ymid,-zmid)

		# Text
		# create a text actor for Temperature
		self.temptxt = vtk.vtkTextActor()
		temp='{:4.3f}'.format(asd.inputdata.get_temp())
		self.temptxt.SetInput("T = " + temp + " K")
		temptxtprop=self.temptxt.GetTextProperty()
		temptxtprop.SetFontFamilyToArial()
		temptxtprop.SetFontSize(30)
		temptxtprop.SetColor(0,0,0)
		temptxtprop.BoldOn()
		self.temptxt.SetDisplayPosition(20,1350)

		# create a text actor for Field
		self.fieldtxt = vtk.vtkTextActor()
		Bfield = asd.inputdata.get_array_hfield()
		self.fieldtxt.SetInput(f"B = ({Bfield[0]:4.1f}, {Bfield[1]:4.1f}, {Bfield[2]:4.1f} ) T")
		fieldtxtprop=self.fieldtxt.GetTextProperty()
		fieldtxtprop.SetFontFamilyToArial()
		fieldtxtprop.SetFontSize(30)
		fieldtxtprop.SetColor(0,0,0)
		fieldtxtprop.BoldOn()
		self.fieldtxt.SetDisplayPosition(20,1300)


		# create a text actor for Energy
		enetxt = vtk.vtkTextActor()
		ene='{:6.4f}'.format(asd.simulationdata.get_total_energy())
		enetxt.SetInput("E= "+ene+" mRy/atom")
		enetxtprop=enetxt.GetTextProperty()
		enetxtprop.SetFontFamilyToArial()
		enetxtprop.SetFontSize(20)
		enetxtprop.SetColor(0,0,0)
		enetxtprop.BoldOn()
		enetxt.SetDisplayPosition(20, 50)

		# LIGHTS ON
		light = vtk.vtkLight()
		light.SetColor(1.0, 1.0, 1.0)
		self.ren.AddLight(light)

		# Reposition the camera
		# z-direction
		self.ren.GetActiveCamera().Azimuth(0)
		self.ren.GetActiveCamera().Elevation(0)
		# xy-direction (almost)
		self.ren.GetActiveCamera().ParallelProjectionOn()
		#self.ren.GetActiveCamera().ParallelProjectionOff()
		d = self.ren.GetActiveCamera().GetDistance()
		#self.ren.GetActiveCamera().SetParallelScale(0.55*max(zmax,xmax))
		self.ren.GetActiveCamera().SetFocalPoint(0,0,0)
		self.ren.GetActiveCamera().SetViewUp(-0.866025403784439,0.5,0)
		#self.ren.GetActiveCamera().SetViewUp(0,1,0)
		self.ren.GetActiveCamera().SetParallelScale(0.55*ymax)
		#ren.GetActiveCamera().SetPosition(0,0,-100*d)
		#print ren.GetActiveCamera().GetViewAngle()
		l = max(xmax-xmin,zmax-zmin)/2
		h = l/0.26795*1.1
		#print l,h

		self.ren.GetActiveCamera().SetPosition(0,0,h)
		#ren.GetActiveCamera().SetPosition(0,0,50*d)

		# For the interactive control. Set up a check for aborting rendering.
		def CheckAbort(obj, event):
				# obj will be the object generating the event.  In this case it
				# is renWin.    
				if obj.GetEventPending() != 0:
					obj.SetAbortRender(1)

		# a minimal keyboard interface
		def Keypress(obj, event):
			key = obj.GetKeySym()
			if key == "w":
				self.Screenshot(self.renWin)
				print("Screenshot taken")
			if key == "0":
				print("Resetting UppASD")
				asd.momentdata.get_array_emom()[:,:,0] =deepcopy(self.initmom.T)
				asd.momentdata.get_array_emom()[:,:,0]=deepcopy(self.initmom.T)
				vecz=numpy_support.numpy_to_vtk(self.initmom)
				colz=numpy_support.numpy_to_vtk(self.initcol)
				self.Datatest.GetPointData().SetVectors(vecz)
				self.Datatest.GetPointData().SetScalars(colz)
			if key == "M":
				print("Running UppASD")
				asd.inputdata.ipmode='M'
				for iter in range(10):
					asd.pyasd.initialphase()
					currmom=asd.momentdata.get_array_emom()[:,:,0].T
					currcol=asd.momentdata.get_array_emom()[2,:,0].T
					vecz=numpy_support.numpy_to_vtk(currmom)
					colz=numpy_support.numpy_to_vtk(currcol)
					self.Datatest.GetPointData().SetVectors(vecz)
					self.Datatest.GetPointData().SetScalars(colz)
					self.renWin.Render()
			if key == "H":
				print("Running UppASD")
				asd.inputdata.ipmode='H'
				for iter in range(10):
					asd.pyasd.initialphase()
					currmom=asd.momentdata.get_array_emom()[:,:,0].T
					currcol=asd.momentdata.get_array_emom()[2,:,0].T
					vecz=numpy_support.numpy_to_vtk(currmom)
					colz=numpy_support.numpy_to_vtk(currcol)
					self.Datatest.GetPointData().SetVectors(vecz)
					self.Datatest.GetPointData().SetScalars(colz)
					self.renWin.Render()
			if key == "S":
				print("Running UppASD")
				asd.inputdata.mode='S'
				for iter in range(10):
					asd.pyasd.measure()
					currmom=asd.momentdata.emoget_array_emom()[:,:,0].T
					currcol=asd.momentdata.get_array_emom()[2,:,0].T
					vecz=numpy_support.numpy_to_vtk(currmom)
					colz=numpy_support.numpy_to_vtk(currcol)
					self.Datatest.GetPointData().SetVectors(vecz)
					self.Datatest.GetPointData().SetScalars(colz)
					self.renWin.Render()
			if key == "m":
				print("Running UppASD")
				asd.inputdata.ipmode='M'
				asd.pyasd.initialphase()
				currmom=asd.momentdata.emom[:,:,0].T
				currcol=asd.momentdata.emom[2,:,0].T
				vecz=numpy_support.numpy_to_vtk(currmom)
				colz=numpy_support.numpy_to_vtk(currcol)
				self.Datatest.GetPointData().SetVectors(vecz)
				self.Datatest.GetPointData().SetScalars(colz)
			if key == "h":
				print("Running UppASD")
				asd.inputdata.ipmode='H'
				asd.pyasd.initialphase()
				currmom=asd.momentdata.emom[:,:,0].T
				currcol=asd.momentdata.emom[2,:,0].T
				vecz=numpy_support.numpy_to_vtk(currmom)
				colz=numpy_support.numpy_to_vtk(currcol)
				self.Datatest.GetPointData().SetVectors(vecz)
				self.Datatest.GetPointData().SetScalars(colz)
			if key == "s":
				print("Running UppASD")
				asd.inputdata.mode='S'
				asd.pyasd.measure()
				currmom=asd.momentdata.get_array_emom()[:,:,0].T
				currcol=asd.momentdata.get_array_emom()[2,:,0].T
				vecz=numpy_support.numpy_to_vtk(currmom)
				colz=numpy_support.numpy_to_vtk(currcol)
				self.Datatest.GetPointData().SetVectors(vecz)
				self.Datatest.GetPointData().SetScalars(colz)
				self.renWin.Render()
				
			if key == "C":
				asd.inputdata.ipmode='Q'
				asd.pyasd.initialphase()
				currmom=asd.momentdata.emom[:,:,0].T
				currcol=asd.momentdata.emom[2,:,0].T
				vecz=numpy_support.numpy_to_vtk(currmom)
				colz=numpy_support.numpy_to_vtk(currcol)
				self.Datatest.GetPointData().SetVectors(vecz)
				self.Datatest.GetPointData().SetScalars(colz)
			if key == "X":
				asd.inputdata.ipmode='Y'
				asd.pyasd.initialphase()
				currmom=asd.momentdata.emom[:,:,0].T
				currcol=asd.momentdata.emom[2,:,0].T
				vecz=numpy_support.numpy_to_vtk(currmom)
				colz=numpy_support.numpy_to_vtk(currcol)
				self.Datatest.GetPointData().SetVectors(vecz)
				self.Datatest.GetPointData().SetScalars(colz)
			if key == "B":
				asd.inputdata.iphfield[2] = asd.inputdata.iphfield[2] + 1.0
				bz='{:4.1f}'.format(asd.inputdata.iphfield[2])
				self.fieldtxt.SetInput("Bz = "+bz+" T")
				self.renWin.Render()
			if key == "b":
				asd.inputdata.iphfield[2] = asd.inputdata.iphfield[2] - 1.0
				bz='{:4.1f}'.format(asd.inputdata.iphfield[2])
				self.fieldtxt.SetInput("Bz = "+bz+" T")
			if key == "N":
				asd.inputdata.iphfield[2] = asd.inputdata.iphfield[2] + 10.0
				bz='{:4.1f}'.format(asd.inputdata.iphfield[2])
				self.fieldtxt.SetInput("Bz = "+bz+" T")
			if key == "n":
				asd.inputdata.iphfield[2] = asd.inputdata.iphfield[2] - 10.0
				bz='{:4.1f}'.format(asd.inputdata.iphfield[2])
				self.fieldtxt.SetInput("Bz = "+bz+" T")
			if key == "T":
				asd.inputdata.iptemp = asd.inputdata.iptemp + 10
				temp='{:4.3f}'.format(asd.inputdata.iptemp[0])
				self.temptxt.SetInput("T = " + temp + " K")
				self.renWin.Render()
			if key == "t":
				asd.inputdata.iptemp =  max(asd.inputdata.iptemp-1.0,1.0e-5)
				temp='{:4.3f}'.format(asd.inputdata.iptemp[0])
				self.temptxt.SetInput("T = " + temp + " K")
			if key == "Y":
				asd.inputdata.iptemp = asd.inputdata.iptemp + 10.0
				temp='{:4.3f}'.format(asd.inputdata.iptemp[0])
				self.temptxt.SetInput("T = " + temp + " K")
			if key == "y":
				asd.inputdata.iptemp =  max(asd.inputdata.iptemp-10.0,1.0e-5)
				temp='{:4.3f}'.format(asd.inputdata.iptemp[0])
				self.temptxt.SetInput("T = " + temp + " K")
			if key == 'g':
					vector.GetProperty().SetInterpolationToGouraud()
			if key == 'p':
					vector.GetProperty().SetInterpolationPBR()
			if key == "m":
				vector.SetVisibility(abs(1-vector.GetVisibility()))
			if key == "i":
				#print('-....-')
				cdata=asd.momentdata.emom[2,:,0].T
				xdim=int(np.sqrt(cdata.shape)[0])
				cmat=cdata.reshape((xdim,xdim))
				fmat=np.abs(np.fft.fft2(cmat.T))
				fmat[0,0]=0.0
				fsgau=gaussian_filter(fmat, sigma=1.0)
				fsmat=np.fft.fftshift(fsgau)

				plt.figure(1)
				plt.clf()
				plt.xticks([], [])
				plt.yticks([], [])
				plt.imshow(fsmat,cmap=cm.Reds)
				plt.title("Static correlation S(q)")
				plt.xlabel("qx")
				plt.ylabel("qy")
				plt.ion()
				plt.show()

			if key == "a":
				atom.SetVisibility(abs(1-atom.GetVisibility()))
				asd.pyasd.totalenergy()
				ene='{:6.4f}'.format(asd.simulationdata.total_energy)
				enetxt.SetInput("E= "+ene+" mRy/atom")
				self.renWin.Render()

		#renWin.AddObserver("AbortCheckEvent", CheckAbort)

		# Add the actors to the renderer, set the background and size
		# Atoms
		atom.SetVisibility(0)
		self.ren.AddActor(atom)
		#ren.AddActor(txt)
		# Vectors
		self.ren.AddActor(vector)
		# Text
		self.ren.AddActor(self.temptxt)
		self.ren.AddActor(self.fieldtxt)
		self.ren.AddActor(enetxt)

		# Outline
		#ren.AddActor(outline)
		# Scalar bar
		#ren.AddActor(scalarBar)

		# Render scene
		# iren = vtk.vtkRenderWindowInteractor()
		#istyle = vtk.vtkInteractorStyleTrackballCamera()
		#iren.SetInteractorStyle(istyle)
		#iren.SetRenderWindow(renWin)
		self.iren.AddObserver("KeyPressEvent", Keypress)

		#ren.ResetCamera()
		#iren.Initialize()
		#
		#
		#cb = vtkTimerCallback()
		###cb.AddActor = vector
		#iren.AddObserver('TimerEvent', cb.execute)
		#timerId = iren.CreateRepeatingTimer(100);
		#iren.SetStillUpdateRate(0.050)
		#iren.SetDesiredUpdateRate(0.050)


		self.renWin.Render()
		self.iren.Start()
		# Screenshot(renWin)

	def S_Step(self):
		""" Do a simulation using S-mode. """
		asd.inputdata.set_mode('S')
		asd.pyasd.measure()
		currmom=asd.momentdata.get_array_emom()[:,:,0].T
		currcol=asd.momentdata.get_array_emom()[2,:,0].T
		vecz=numpy_support.numpy_to_vtk(currmom)
		colz=numpy_support.numpy_to_vtk(currcol)
		self.Datatest.GetPointData().SetVectors(vecz)
		self.Datatest.GetPointData().SetScalars(colz)
		self.renWin.Render()

	def M_step(self):
		""" Do a simulation using M-mode. """
		asd.inputdata.set_mode('M')
		asd.pyasd.measure()
		currmom=asd.momentdata.get_array_emom()[:,:,0].T
		currcol=asd.momentdata.get_array_emom()[2,:,0].T
		vecz=numpy_support.numpy_to_vtk(currmom)
		colz=numpy_support.numpy_to_vtk(currcol)
		self.Datatest.GetPointData().SetVectors(vecz)
		self.Datatest.GetPointData().SetScalars(colz)
		self.renWin.Render()

	def Reset(self):
		""" Reset data to initial. """
		asd.momentdata.get_array_emom()[:,:,0] =deepcopy(self.initmom.T)
		asd.momentdata.get_array_emom()[:,:,0]=deepcopy(self.initmom.T)
		vecz=numpy_support.numpy_to_vtk(self.initmom)
		colz=numpy_support.numpy_to_vtk(self.initcol)
		self.Datatest.GetPointData().SetVectors(vecz)
		self.Datatest.GetPointData().SetScalars(colz)

	def UpdateTemperature(self):
		""" Update temperature actor. """
		temp='{:4.3f}'.format(asd.inputdata.get_temp())
		self.temptxt.SetInput("T = " + temp + " K")
		self.renWin.Render()

	def UpdateBfield(self):
		""" Update B-field actor. """
		Bfield = asd.inputdata.get_array_hfield()
		self.fieldtxt.SetInput(f"B = ({Bfield[0]:4.1f}, {Bfield[1]:4.1f}, {Bfield[2]:4.1f} ) T")
		self.renWin.Render()

	def close_window(self):
		render_window = self.iren.GetRenderWindow()
		render_window.Finalize()
		self.iren.TerminateApp()

	# Read Location of Atoms
	def readAtoms(self, file, Nmax):
		points = vtk.vtkPoints()
		nrAtoms=0
		# Read ahead
		line = file.readline()
		data = line.split()

		# Read all data
		while line:
			if nrAtoms <= Nmax:
				a = int(data[0])
				x, y, z = float(data[1]), float(data[2]), float(data[3])
				#print "a ", a, " x ", x, " y ", y, " z ", z  
				#points.InsertPoint(a, x, y, z)
				points.InsertPoint(nrAtoms, x, y, z)
			nrAtoms=nrAtoms+1
			line = file.readline()
			data = line.split()
		return points, nrAtoms

	# Read vectors
	# We must know the time step and the number of atoms per time
	def readVectorsData(file, time, nrAtoms,Nmax):
		# Create a Double array which represents the vectors
		vectors = vtk.vtkFloatArray()
		colors = vtk.vtkFloatArray()

		# Define number of elemnts
		vectors.SetNumberOfComponents(3)
		colors.SetNumberOfComponents(1)
		for i in range(7):
			line = file.readline()

		i=0
		# Read all data for a certain time
		while i < nrAtoms:
			line = file.readline()
			if i < Nmax:
				data = line.split()
				t, a = int(data[0]), int(data[1])
				x, y, z = float(data[4]), float(data[5]), float(data[6])
				#x, y, z = float(data[3]), float(data[4]), float(data[5])
				m = (z+1.0)/2.0
				#m = (x+y+2.0)/2.0
				#m = (z+y+x)/3.0**0.5
				#m = atan2(x,y)/acos(-1)+1
				#m = (atan2(x,z)/acos(-1)+1)/2
				#print m
				vectors.InsertTuple3(i ,x, y, z)
				colors.InsertValue(i ,m)
				i=i+1
		return vectors, colors
		
		#A function that takes a renderwindow and saves its contents to a .png file
	def Screenshot(self):
		self.number_of_screenshots
		win2im=vtk.vtkWindowToImageFilter()
		win2im.ReadFrontBufferOff()
		win2im.SetInput(self.renWin)
		#
		povexp=vtk.vtkPOVExporter()
		povexp.SetRenderWindow(self.renWin)
		#povexp.SetInput(renWin)
		self.renWin.Render()
		povexp.SetFileName('snap%.5d.pov' %self.number_of_screenshots)
		povexp.Write()
		#
		toPNG=vtk.vtkPNGWriter()
		toPNG.SetFileName('snap%.5d.png' %self.number_of_screenshots)
		#toPNG.SetInput(win2im.GetOutput())
		toPNG.SetInputConnection(win2im.GetOutputPort())
		toPNG.Write()

		self.number_of_screenshots += 1
		return

	def UpdateTextPlacement():
		pass