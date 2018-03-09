#!/usr/bin/env vtkpython
################################################################################
# CLASS: ASDMomActors
# @author Jonathan Chico (08/09/2017)
# @description
# Wrapper class to add the VTK  actors for the visualization of UppASD data in
# moment visualization mode.
# It contains the needed data to add the actors, modify them, as well as some helper
# functions to change them.
################################################################################
import vtk
import time
import numpy as np
import ASDVTKReading
from PyQt4 import QtCore, QtGui

class ASDMomActors():

	############################################################################
	# Main wrapper to add the needed actors for visualization
	############################################################################
	""" This is the routine where the actual VTK actors are created for the
	Moment and Restart mode. """
	def AddASD_actors(self,ren,renWin,mode,viz_type,iren):
		ASDMomActors.timer_count=0
		ASDMomActors.camera_pos=np.zeros(3,dtype=np.float32)
		ASDMomActors.camera_focal=np.zeros(3,dtype=np.float32)
		ASDMomActors.camera_yaw=0.0
		ASDMomActors.camera_roll=0.0
		ASDMomActors.camera_pitch=0.0
		ASDMomActors.camera_azimuth=0.0
		ASDMomActors.camera_elevation=0.0

		ASD_data=ASDVTKReading.ASDReading()
	    # Add the data structures with regards to reading the data
		ASD_data.ReadingWrapper(mode=mode,viz_type=viz_type)

		ASDMomActors.kmc_disp=ASD_data.kmc_flag
		ASDMomActors.cluster_disp=ASD_data.cluster_flag
	    ########################################################################
	    # Data structures for the generation of the smooth grid
	    ########################################################################
		ASDMomActors.glob_flag_2D=ASD_data.flag_2D
		ASDMomActors.glob_color_x=ASD_data.selected_colors_x
		ASDMomActors.glob_color_y=ASD_data.selected_colors_y
		ASDMomActors.glob_color_z=ASD_data.selected_colors_z
	    ########################################################################
	    # Look up tables for colors
	    ########################################################################
	    # This is a diverging RWB color mapping based on the work of Kenneth
	    # Moreland and with the vtk examples provided by Andrew Maclean
		if ASDMomActors.glob_flag_2D:
			self.lut = vtk.vtkLookupTable()
			num_colors = 256
			self.lut.SetNumberOfTableValues(num_colors)

			self.transfer_func = vtk.vtkColorTransferFunction()
			self.transfer_func.SetColorSpaceToDiverging()
			self.transfer_func.AddRGBPoint(0, 0.230, 0.299, 0.754)
			self.transfer_func.AddRGBPoint(1, 0.706, 0.016, 0.150)

			for ii,ss in enumerate([float(xx)/float(num_colors) for xx in range(num_colors)]):
				cc = self.transfer_func.GetColor(ss)
				self.lut.SetTableValue(ii, cc[0], cc[1], cc[2], 1.0)
			self.lut.Build()
		else:
			self.lut = vtk.vtkLookupTable()
			num_colors = 256
			self.lut.SetNumberOfTableValues(num_colors)

			self.transfer_func = vtk.vtkColorTransferFunction()
			self.transfer_func.SetColorSpaceToDiverging()
			self.transfer_func.AddRGBPoint(-1, 0.230, 0.299, 0.754)
			self.transfer_func.AddRGBPoint( 1, 0.706, 0.016, 0.150)

	        for ii,ss in enumerate([float(xx)/float(num_colors) for xx in range(num_colors)]):
	            cc = self.transfer_func.GetColor(ss)
	            self.lut.SetTableValue(ii, cc[0], cc[1], cc[2], 1.0)
	        self.lut.Build()

	    ########################################################################
	    # Data structures for the generation of the smooth grid
	    ########################################################################

		# Passing the data from the full system to the PolyData
		ASDMomActors.src=vtk.vtkPolyData()
		ASDMomActors.src.SetPoints(ASD_data.selected_points)
		ASDMomActors.src.GetPointData().SetScalars(ASD_data.selected_colors_z)
		ASDMomActors.src.GetPointData().SetVectors(ASD_data.selected_vectors)
		scalar_range = ASDMomActors.src.GetScalarRange()

	    ########################################################################
	    # Finding useful geometrical information of the sample
	    ########################################################################
	    # Finding the middle of the sample
	    # Also making sure that if the sample is 2D one has no problem with boudings
	    # this is mostly useful if splatters are used
		(ASDMomActors.xmin,ASDMomActors.xmax,ASDMomActors.ymin,ASDMomActors.ymax,ASDMomActors.zmin,ASDMomActors.zmax)= ASDMomActors.src.GetBounds()
		if ASDMomActors.xmin==ASDMomActors.xmax:
			ASDMomActors.xmin=0.0
			ASDMomActors.xmax=1.0
		if ASDMomActors.ymin==ASDMomActors.ymax:
			ASDMomActors.ymin=0.0
			ASDMomActors.ymax=1.0
		if ASDMomActors.zmin==ASDMomActors.zmax:
			ASDMomActors.zmin=0.0
			ASDMomActors.zmax=1.0
		ASDMomActors.xmid = (ASDMomActors.xmin+ASDMomActors.xmax)*0.5
		ASDMomActors.ymid = (ASDMomActors.ymin+ASDMomActors.ymax)*0.5
		ASDMomActors.zmid = (ASDMomActors.zmin+ASDMomActors.zmax)*0.5
		ASDMomActors.height=max(ASDMomActors.xmax,ASDMomActors.ymax,ASDMomActors.zmax)*1.75
		self.dist_x=np.absolute(ASDMomActors.xmax-ASDMomActors.xmin)
		self.dist_y=np.absolute(ASDMomActors.ymax-ASDMomActors.ymin)
		self.dist_z=np.absolute(ASDMomActors.zmax-ASDMomActors.zmin)
		ASDMomActors.camera_pos[0]=ASDMomActors.xmid
		ASDMomActors.camera_pos[1]=ASDMomActors.ymid
		ASDMomActors.camera_pos[2]=ASDMomActors.height
		ASDMomActors.camera_focal[0]=ASDMomActors.xmid
		ASDMomActors.camera_focal[1]=ASDMomActors.ymid
		ASDMomActors.camera_focal[2]=ASDMomActors.zmid
	    # The delaunay tesellation seems to be the best way to transform the point cloud
	    # to a surface for volume rendering, the problem is that it is too slow for large
	    # data sets, meaning that the best option is first to prune out the data to ensure
	    # that one has a manageable number of data points over which to do the construction

	    # surface reconstruction and splatter techniques also can be used to generate something
	    # akin to the kind of surfaces we want. The issue is that they transform the data to a
	    # regular mesh by default. And thus it is a problem for most kind of systems
		if ASDMomActors.glob_flag_2D:
			# Passing the data to generate a triangulation of the data
			MagDensMethod = vtk.vtkDelaunay2D()
			MagDensMethod.SetInputData(ASDMomActors.src)
			MagDensMethod.BoundingTriangulationOff()
			# Time the execution of the delaunay tessellation
			SM_timer = vtk.vtkExecutionTimer()
			SM_timer.SetFilter(MagDensMethod)
			MagDensMethod.Update()
			SM = SM_timer.GetElapsedWallClockTime()
			print ("2D Delaunay:", SM)

			# Creating the mapper for the smooth surfaces
			ASDMomActors.MagDensMap = vtk.vtkDataSetMapper()
			ASDMomActors.MagDensMap.SetScalarRange(scalar_range)
			ASDMomActors.MagDensMap.SetInputConnection(MagDensMethod.GetOutputPort())
			ASDMomActors.MagDensMap.SetLookupTable(self.lut)
			ASDMomActors.MagDensMap.SetColorModeToMapScalars()
			ASDMomActors.MagDensMap.Update()

			# Creating the actor for the smooth surfaces
			ASDMomActors.MagDensActor = vtk.vtkLODActor()
			ASDMomActors.MagDensActor.SetMapper(ASDMomActors.MagDensMap)
			ASDMomActors.MagDensActor.GetProperty().SetOpacity(0.75)
			ASDMomActors.MagDensActor.GetProperty().EdgeVisibilityOff()

		else:
	        ####################################################################
	        # Setting the parameters for the visualization of 3D structures with
	        # splatters
	        ####################################################################
			MagDensMethod = vtk.vtkShepardMethod()
			MagDensMethod.SetInputData(ASDMomActors.src)
			MagDensMethod.SetModelBounds(ASDMomActors.xmin,ASDMomActors.xmax,ASDMomActors.ymin,ASDMomActors.ymax,ASDMomActors.zmin,ASDMomActors.zmax)
	        # This should get rid of the problems when trying to map very thin structures in 2D
			if self.dist_x==min(self.dist_x,self.dist_y,self.dist_z):
				MagDensMethod.SetSampleDimensions(3,int(ASDMomActors.ymax),int(ASDMomActors.zmax))
			elif self.dist_y==min(self.dist_x,self.dist_y,self.dist_z):
				MagDensMethod.SetSampleDimensions(int(ASDMomActors.xmax),3,int(ASDMomActors.zmax))
			elif self.dist_z==min(self.dist_x,self.dist_y,self.dist_z):
				MagDensMethod.SetSampleDimensions(int(ASDMomActors.xmax),int(ASDMomActors.ymax),3)
	        # This parameter determines how far in the sample (normalized to 1) the
	        # method will look to interpolate, greatly affects performance
			MagDensMethod.SetMaximumDistance(0.1)

			# Time the execution of the checkerboard splatter
			SP_timer = vtk.vtkExecutionTimer()
			SP_timer.SetFilter(MagDensMethod)
			MagDensMethod.Update()
			SP = SP_timer.GetElapsedWallClockTime()
			print ("3D Shepard Method:", SP)

			# Mapper for the image obtained from the 3D reconstruction method
			ASDMomActors.MagDensMap = vtk.vtkSmartVolumeMapper()
			ASDMomActors.MagDensMap.SetBlendModeToComposite()
			ASDMomActors.MagDensMap.SetInputConnection(MagDensMethod.GetOutputPort())

			# Function for the opacity gradient
			volumeGradientOpacity = vtk.vtkPiecewiseFunction()
			volumeGradientOpacity.AddPoint(-1,0.25)
			volumeGradientOpacity.AddPoint(0.5,0.75)
			volumeGradientOpacity.AddPoint(1.0,1.0)

			# Properties of the volume to be rendered
			self.volumeProperty = vtk.vtkVolumeProperty()
			self.volumeProperty.SetInterpolationType(1)
			self.volumeProperty.SetColor(self.transfer_func)
			self.volumeProperty.SetAmbient(0.6)
			self.volumeProperty.SetDiffuse(0.6)
			self.volumeProperty.SetSpecular(0.1)
			self.volumeProperty.SetGradientOpacity(volumeGradientOpacity)
			#volumeProperty.ShadeOn()

			# Volume actor, this works in a different way than LOD actors
			ASDMomActors.MagDensActor = vtk.vtkVolume()
			ASDMomActors.MagDensActor.SetMapper(ASDMomActors.MagDensMap)
			ASDMomActors.MagDensActor.SetProperty(self.volumeProperty)

	        ####################################################################
	        # Alternative rendering methods
	        ####################################################################
	        # The checkerboard splatter method, much faster than the Shepard method
	        # however, it seems to change the values of the scalars embeded in the
	        # point data which results in incorrect diplay of the magnetization
	        #cbdSplatter = vtk.vtkCheckerboardSplatter()
	        #cbdSplatter.ScalarWarpingOff()
	        #cbdSplatter.SetFootprint(2)
	        #cbdSplatter.SetExponentFactor(-5)
	        #cbdSplatter.SetParallelSplatCrossover(2)
	        #cbdSplatter.SetOutputScalarTypeToDouble()
	        #cbdSplatter.CappingOn()
	        #cbdSplatter.ScalarWarpingOn()
	        #cbdSplatter.SetRadius(1)

	        # 3D delaunay method, the far superior as it conserves the shape of the
	        # sample, however it is extremely slow, with a rendering of a 3D image
	        # taking several minutes
	        #smooth_loop = vtk.vtkDelaunay3D()
	        #smooth_loop.SetInputData(self.src)
	        #smooth_loop.SetTolerance(0.01)
	        #smooth_loop.SetAlpha(2)
	        #smooth_loop.AlphaTrisOff()
	        #smooth_loop.Update()

	    ########################################################################
	    # Data structures for the spins
	    ########################################################################
	    # Passing the data from the full system to the PolyData
		ASDMomActors.src_spins=vtk.vtkPolyData()
		ASDMomActors.src_spins.SetPoints(ASD_data.selected_points)
		ASDMomActors.src_spins.GetPointData().SetScalars(ASD_data.selected_colors_z)
		ASDMomActors.src_spins.GetPointData().SetVectors(ASD_data.selected_vectors)
		scalar_range_spins = ASDMomActors.src_spins.GetScalarRange()
		########################################################################
		# Data structures for the contours
		########################################################################
		# Define the contour filters
		contours = vtk.vtkContourFilter()
		contours.SetInputConnection(MagDensMethod.GetOutputPort())
		# This generates the contours, it will do 5 between the -1 and 0.5 range
		cont_num=5
		range_cont=(-1,0.5)
		contours.GenerateValues(cont_num,range_cont)
		# Map the contours to graphical primitives
		contMapper = vtk.vtkPolyDataMapper()
		contMapper.SetInputConnection(contours.GetOutputPort())
		contMapper.SetScalarVisibility(False) # colored contours
		contMapper.SetScalarRange(scalar_range)
		# Create an actor for the contours
		ASDMomActors.contActor = vtk.vtkLODActor()
		ASDMomActors.contActor.SetMapper(contMapper)
		ASDMomActors.contActor.GetProperty().SetColor(0, 0, 0)
		ASDMomActors.contActor.GetProperty().SetLineWidth(1.0)

	    ########################################################################
	    # Data structures for the impurity cluster
	    ########################################################################
		if ASD_data.cluster_flag:

			# Passing the data from the cluster to the PolyData
			src_clus=vtk.vtkPolyData()
			src_clus.SetPoints(ASD_data.coord_c)
			src_clus.GetPointData().SetScalars(ASD_data.colors_clus)

			# Passing the data from the selected impurities
			src_imp=vtk.vtkPolyData()
			src_imp.SetPoints(ASD_data.points_clus_imp)
			src_imp.GetPointData().SetScalars(ASD_data.colors_imp)
			src_imp.Modified()

			atomSource = vtk.vtkDelaunay2D()
			atomSource.SetInputData(src_clus)
			atomSource.BoundingTriangulationOff()
			atomSource.Update()

 			smoothFilter =vtk.vtkSmoothPolyDataFilter()
			smoothFilter.SetInputConnection(atomSource.GetOutputPort())
			smoothFilter.SetNumberOfIterations(5)
			smoothFilter.SetRelaxationFactor(0.1)
			smoothFilter.FeatureEdgeSmoothingOff()
			smoothFilter.BoundarySmoothingOn()
			smoothFilter.Update()

			# Creating the mapper for the smooth surfaces
			atomMapper = vtk.vtkDataSetMapper()
			atomMapper.SetScalarRange(scalar_range)
			atomMapper.SetInputConnection(smoothFilter.GetOutputPort())
			atomMapper.SetColorMode(2)
			atomMapper.Update()

			# Creating the actor for the smooth surfaces
			ASDMomActors.atom = vtk.vtkLODActor()
			ASDMomActors.atom.SetMapper(atomMapper)
			ASDMomActors.atom.GetProperty().EdgeVisibilityOff()
			ASDMomActors.atom.GetProperty().SetSpecularPower(30)
			ASDMomActors.atom.GetProperty().SetAmbient(0.2)
			ASDMomActors.atom.GetProperty().SetDiffuse(0.8)
			ASDMomActors.atom.GetProperty().SetOpacity(0.50)

			# Set up imp sources
			atomSource_imp = vtk.vtkSphereSource()
			atomSource_imp.SetRadius(2.5)
			atomSource_imp.SetThetaResolution(20)
			atomSource_imp.SetPhiResolution(20)

			# Mapping the spheres to the actual points on the selected impurities
			atomMapper_imp = vtk.vtkGlyph3DMapper()
			atomMapper_imp.SetInputData(src_imp)
			atomMapper_imp.SetSourceConnection(atomSource_imp.GetOutputPort())
			atomMapper_imp.SetScaleFactor(0.2)
			atomMapper_imp.SetScaleModeToNoDataScaling()
			atomMapper_imp.Update()

			# Creating the selected impurity actors
			ASDMomActors.atom_imp = vtk.vtkLODActor()
			ASDMomActors.atom_imp.SetMapper(atomMapper_imp)
			ASDMomActors.atom_imp.GetProperty().SetSpecular(0.3)
			ASDMomActors.atom_imp.GetProperty().SetSpecularPower(30)
			ASDMomActors.atom_imp.GetProperty().SetAmbient(0.2)
			ASDMomActors.atom_imp.GetProperty().SetDiffuse(0.8)

	    ########################################################################
	    # Setting information of the directions
	    ########################################################################
		# Create vectors
		arrow = vtk.vtkArrowSource()
		arrow.SetTipRadius(0.20)
		arrow.SetShaftRadius(0.10)
		arrow.SetTipResolution(10)
		arrow.SetShaftResolution(10)

		# Create the mapper for the spins
		arrowMapper = vtk.vtkGlyph3DMapper()
		arrowMapper.SetSourceConnection(arrow.GetOutputPort())
		arrowMapper.SetInputData(ASDMomActors.src)
		arrowMapper.SetScaleFactor(0.50)
		arrowMapper.SetScalarVisibility(False)
		arrowMapper.SetScaleModeToNoDataScaling()
		arrowMapper.Update()

		# Define the vector actor for the spins
		ASDMomActors.vector = vtk.vtkLODActor()
		ASDMomActors.vector.SetMapper(arrowMapper)
		ASDMomActors.vector.GetProperty().SetSpecular(0.3)
		ASDMomActors.vector.GetProperty().SetSpecularPower(60)
		ASDMomActors.vector.GetProperty().SetAmbient(0.2)
		ASDMomActors.vector.GetProperty().SetDiffuse(0.8)
		ASDMomActors.vector.GetProperty().SetColor(0, 0, 0)
		ASDMomActors.vector.VisibilityOff()

		########################################################################
		# Setting information of the spins
		########################################################################
		# Create vectors
		ASDMomActors.spinarrow = vtk.vtkArrowSource()
		ASDMomActors.spinarrow.SetTipRadius(0.20)
		ASDMomActors.spinarrow.SetShaftRadius(0.10)
		ASDMomActors.spinarrow.SetTipResolution(10)
		ASDMomActors.spinarrow.SetShaftResolution(10)

		# Create the mapper for the spins
		ASDMomActors.SpinMapper = vtk.vtkGlyph3DMapper()
		ASDMomActors.SpinMapper.SetSourceConnection(ASDMomActors.spinarrow.GetOutputPort())
		ASDMomActors.SpinMapper.SetInputData(ASDMomActors.src_spins)
		ASDMomActors.SpinMapper.SetScalarRange(scalar_range_spins)
		ASDMomActors.SpinMapper.SetScaleFactor(0.50)
		ASDMomActors.SpinMapper.SetScaleModeToNoDataScaling()
		ASDMomActors.SpinMapper.SetLookupTable(self.lut)
		ASDMomActors.SpinMapper.SetColorModeToMapScalars()
		ASDMomActors.SpinMapper.Update()

		# Define the vector actor for the spins
		ASDMomActors.Spins = vtk.vtkLODActor()
		ASDMomActors.Spins.SetMapper(ASDMomActors.SpinMapper)
		ASDMomActors.Spins.GetProperty().SetSpecular(0.3)
		ASDMomActors.Spins.GetProperty().SetSpecularPower(60)
		ASDMomActors.Spins.GetProperty().SetAmbient(0.2)
		ASDMomActors.Spins.GetProperty().SetDiffuse(0.8)
		ASDMomActors.Spins.VisibilityOff()

		########################################################################
		# Creation of the data structures for the data clipping
		########################################################################
		# Right now this only can clip polydata, which is fine for 2D structures
		# however, for the 3d delaunay tesellation, the output is an unstructured
		# grid, which means that annother type of clipper is required
		ASDMomActors.plane = vtk.vtkPlane()
		ASDMomActors.plane.SetOrigin(ASDMomActors.xmin, ASDMomActors.ymid, 0)
		ASDMomActors.plane.SetNormal(1, 0, 0)
		if ASDMomActors.glob_flag_2D:
			self.clipper = vtk.vtkClipPolyData()
			self.clipper.SetInputConnection(MagDensMethod.GetOutputPort())
			self.clipper.SetClipFunction(ASDMomActors.plane)
			self.clipper.InsideOutOn()

			ASDMomActors.clipperMapper = vtk.vtkPolyDataMapper()
			ASDMomActors.clipperMapper.SetInputConnection(self.clipper.GetOutputPort())
			ASDMomActors.clipperMapper.SetLookupTable(self.lut)

			ASDMomActors.clipperActor = vtk.vtkLODActor()
			ASDMomActors.clipperActor.SetMapper(ASDMomActors.clipperMapper)
			ASDMomActors.clipperActor.VisibilityOff()
		else:
			self.clipper = vtk.vtkClipVolume()
			self.clipper.SetInputConnection(MagDensMethod.GetOutputPort())
			self.clipper.SetClipFunction(ASDMomActors.plane)
			self.clipper.InsideOutOn()

			ASDMomActors.clipperMapper = vtk.vtkDataSetMapper()
			ASDMomActors.clipperMapper.SetInputConnection(self.clipper.GetOutputPort())
			ASDMomActors.clipperMapper.SetLookupTable(self.transfer_func)

			ASDMomActors.clipperActor = vtk.vtkActor()
			ASDMomActors.clipperActor.SetMapper(ASDMomActors.clipperMapper)
			ASDMomActors.clipperActor.GetProperty().SetOpacity(1.00)

		if (ASD_data.kmc_flag):
			########################################################################
			# Setting data structures for the KMC particle visualization
			########################################################################
			self.KMC_src=vtk.vtkPolyData()
			self.KMC_src.SetPoints(ASD_data.coord_KMC)

			# Atom sphere
			KMC_part = vtk.vtkSphereSource()
			KMC_part.SetRadius(1.75)
			KMC_part.SetThetaResolution(40)
			KMC_part.SetPhiResolution(40)
			# Atom glyph
			KMC_part_mapper = vtk.vtkGlyph3DMapper()
			KMC_part_mapper.SetInputData(self.KMC_src)
			KMC_part_mapper.SetSourceConnection(KMC_part.GetOutputPort())
			KMC_part_mapper.SetScaleFactor(0.5)
			KMC_part_mapper.ClampingOn()
			KMC_part_mapper.SetScaleModeToNoDataScaling()
			KMC_part_mapper.SetColorModeToMapScalars()
			KMC_part_mapper.Update()
			# Atoms actors
			ASDMomActors.KMC_part_actor = vtk.vtkLODActor()
			ASDMomActors.KMC_part_actor.SetMapper(KMC_part_mapper)
			ASDMomActors.KMC_part_actor.GetProperty().SetOpacity(0.9)
			ASDMomActors.KMC_part_actor.GetProperty().SetColor(0.0, 0.0, 1.0)
			ASDMomActors.KMC_part_actor.GetProperty().EdgeVisibilityOn()
			ASDMomActors.KMC_part_actor.GetProperty().SetEdgeColor(0,0,0)
		########################################################################
		# Setting the information for the axes widget
		########################################################################
		# Create the axes actor
		axes = vtk.vtkAxesActor()
		axes.SetShaftTypeToCylinder()
		axes.SetCylinderRadius(0.05)
		axes.SetNormalizedShaftLength(0.9,0.9,0.9)
		axes.SetNormalizedTipLength(0.40,0.40,0.40)
		# The properties of the text can be controlled independently
		axes.GetXAxisCaptionActor2D().GetCaptionTextProperty().SetColor(0.0,0.0,0.0)
		axes.GetYAxisCaptionActor2D().GetCaptionTextProperty().SetColor(0.0,0.0,0.0)
		axes.GetZAxisCaptionActor2D().GetCaptionTextProperty().SetColor(0.0,0.0,0.0)
		axes.GetXAxisCaptionActor2D().GetCaptionTextProperty().ShadowOff()
		axes.GetYAxisCaptionActor2D().GetCaptionTextProperty().ShadowOff()
		axes.GetZAxisCaptionActor2D().GetCaptionTextProperty().ShadowOff()

		# The axes actor is then used as an orientation marker widget, the advantage
		# of setting it up as a widget is that it is interactive and one can move it
		# and that it moves as the zoom changes
		# Must make sure that the widget is part of the main class so that it can
		# be actually rendered and no segfaults occurr
		ASDMomActors.OrientMarker= vtk.vtkOrientationMarkerWidget()
		ASDMomActors.OrientMarker.SetOutlineColor(0.9300, 0.5700, 0.1300)
		ASDMomActors.OrientMarker.SetOrientationMarker(axes)
		ASDMomActors.OrientMarker.SetViewport(0.0, 0.0, 0.3, 0.3)

		########################################################################
		# Setting the information for the scalar bar widget
		########################################################################
		# Create the scalar bar actor
		ASDMomActors.scalar_bar = vtk.vtkScalarBarActor()
		if ASDMomActors.glob_flag_2D:
			ASDMomActors.scalar_bar.SetLookupTable(ASDMomActors.MagDensMap.GetLookupTable())
		else:
			ASDMomActors.scalar_bar.SetLookupTable(self.transfer_func)

		ASDMomActors.scalar_bar.GetLabelTextProperty().SetColor(0.0,0.0,0.0)
		ASDMomActors.scalar_bar.SetNumberOfLabels(5)
		ASDMomActors.scalar_bar.GetLabelTextProperty().ShadowOff()
		ASDMomActors.scalar_bar.GetLabelTextProperty().BoldOff()
		ASDMomActors.scalar_bar.GetLabelTextProperty().ItalicOff()
		ASDMomActors.scalar_bar.SetLabelFormat("%-#6.1f")
		ASDMomActors.scalar_bar.SetBarRatio(0.5)

		# Create the scalar_bar_widget
		ASDMomActors.scalar_bar_widget = vtk.vtkScalarBarWidget()
		ASDMomActors.scalar_bar_widget.SetScalarBarActor(ASDMomActors.scalar_bar)

		# Representation to actually control where the scalar bar is
		scalarBarRep = ASDMomActors.scalar_bar_widget.GetRepresentation()
		scalarBarRep.SetOrientation(0)  # 0 = Horizontal, 1 = Vertical
		scalarBarRep.GetPositionCoordinate().SetValue(0.30,0.05)
		scalarBarRep.GetPosition2Coordinate().SetValue(0.50,0.05)

		########################################################################
		# Setting information of the renderer
		########################################################################
		# Define the renderer
		# Add the actors to the scene
		if ASDMomActors.glob_flag_2D:
			ren.AddActor(ASDMomActors.MagDensActor)
		else:
			ren.AddViewProp(ASDMomActors.MagDensActor)

		ren.AddActor(ASDMomActors.Spins)
		ren.AddActor(ASDMomActors.vector)
		ren.AddActor(ASDMomActors.contActor)
		ren.AddActor(self.clipperActor)

		# If there is information about the cluster add the needed actors
		if ASD_data.cluster_flag:
			ren.AddActor(ASDMomActors.atom)
			ren.AddActor(ASDMomActors.atom_imp)

	    #If the KMC particles are present add them to the renderer
		if ASD_data.kmc_flag:
			ren.AddActor(ASDMomActors.KMC_part_actor)
	    # Defining the camera directions

		ren.GetActiveCamera().Azimuth(ASDMomActors.camera_azimuth)
		ren.GetActiveCamera().Elevation(ASDMomActors.camera_elevation)
		ren.GetActiveCamera().Yaw(ASDMomActors.camera_yaw)
		ren.GetActiveCamera().Roll(ASDMomActors.camera_roll)
		ren.GetActiveCamera().Pitch(ASDMomActors.camera_pitch)
		ren.GetActiveCamera().SetFocalPoint(self.camera_focal)
		ren.GetActiveCamera().SetPosition(self.camera_pos)
		ren.GetActiveCamera().SetViewUp(0,1,0)
		# Must make sure the widgets is called before the renderer is called
		# Scalar bar
		ASDMomActors.scalar_bar_widget.SetInteractor(iren)
		ASDMomActors.scalar_bar_widget.On()
		# Orient marker
		ASDMomActors.OrientMarker.SetInteractor(iren)
		ASDMomActors.OrientMarker.SetEnabled(1)
		########################################################################
		# Start the renderer
		########################################################################
		iren.Start()
		renWin.Render()

		return;

    ############################################################################
    # A function that takes a renderwindow and saves its contents to a .png file
    # @author Anders Bergman
    ############################################################################
	def Screenshot(self,renWin,number_of_screenshots,png_mode,pov_mode):

		win2im=vtk.vtkWindowToImageFilter()
		win2im.SetInput(renWin)
		win2im.Update()
		win2im.SetInputBufferTypeToRGBA()
		win2im.ReadFrontBufferOff()

		if pov_mode:
			povexp=vtk.vtkPOVExporter()
			povexp.SetInput(renWin)
			renWin.Render()
			povexp.SetFileName('snap%.5d.pov' %number_of_screenshots)
			povexp.Write()

		if png_mode:
			toPNG=vtk.vtkPNGWriter()
			toPNG.SetFileName('snap%.5d.png' %number_of_screenshots)
			toPNG.SetInputConnection(win2im.GetOutputPort())
			toPNG.Write()

		return;

	############################################################################
	# Reseting the camera to the initial positions
	############################################################################
	def reset_camera(self,ren,renWin):

	    # Defining the camera directions
		ren.GetActiveCamera().SetFocalPoint(ASDMomActors.xmid,ASDMomActors.ymid,ASDMomActors.zmid)
		ren.GetActiveCamera().SetPosition(ASDMomActors.xmid,ASDMomActors.ymid,ASDMomActors.height)
		ren.GetActiveCamera().Azimuth(0)
		ren.GetActiveCamera().Elevation(0)
		ren.GetActiveCamera().Yaw(0)
		ren.GetActiveCamera().Roll(0)
		ren.GetActiveCamera().Pitch(0)
		ren.GetActiveCamera().SetViewUp(0,1,0)
		renWin.Render()

		return

	############################################################################
	# Set the camera view up to be defined to the (1,0,0)
	############################################################################
	def set_Camera_x(self,ren,renWin):
		ren.GetActiveCamera().SetViewUp(1,0,0)
		renWin.Render()

	############################################################################
	# Set the camera view up to be defined to the (0,1,0)
	############################################################################
	def set_Camera_y(self,ren,renWin):
		ren.GetActiveCamera().SetViewUp(0,1,0)
		renWin.Render()

	############################################################################
	# Set the camera view up to be defined to the (0,0,1)
	############################################################################
	def set_Camera_z(self,ren,renWin):
		ren.GetActiveCamera().SetViewUp(0,0,1)
		renWin.Render()

	def ChangeSpinGlyph(self,renWin,keyword):
		if keyword=='Cubes':
			try:
				del ASDMomActors.spinarrow
			except:
				pass
			try:
				del self.spinsphere
			except:
				pass
			try:
				del self.spincone
			except:
				pass

			self.spincube=vtk.vtkCubeSource()
			self.spincube.SetXLength(1.0)
			self.spincube.SetYLength(1.0)
			self.spincube.SetZLength(1.0)
			ASDMomActors.SpinMapper.SetSourceConnection(self.spincube.GetOutputPort())
			ASDMomActors.SpinMapper.ClampingOn()
			ASDMomActors.SpinMapper.OrientOff()
			renWin.Render()
		if keyword=='Spheres':
			try:
				del ASDMomActors.spinarrow
			except:
				pass
			try:
				del self.spincube
			except:
				pass
			try:
				del self.spincone
			except:
				pass
			self.spinsphere = vtk.vtkSphereSource()
			self.spinsphere.SetRadius(1.00)
			self.spinsphere.SetThetaResolution(20)
			self.spinsphere.SetPhiResolution(20)
			ASDMomActors.SpinMapper.SetSourceConnection(self.spinsphere.GetOutputPort())
			ASDMomActors.SpinMapper.ClampingOn()
			ASDMomActors.SpinMapper.OrientOff()
			renWin.Render()
		if keyword=='Arrows':
			try:
				del self.spinsphere
			except:
				pass
			try:
				del self.spincube
			except:
				pass
			try:
				del self.spincone
			except:
				pass
			ASDMomActors.spinarrow = vtk.vtkArrowSource()
			ASDMomActors.spinarrow.SetTipRadius(0.20)
			ASDMomActors.spinarrow.SetShaftRadius(0.10)
			ASDMomActors.spinarrow.SetTipResolution(10)
			ASDMomActors.spinarrow.SetShaftResolution(10)
			ASDMomActors.SpinMapper.SetSourceConnection(ASDMomActors.spinarrow.GetOutputPort())
			ASDMomActors.SpinMapper.OrientOn()
			renWin.Render()
		if keyword=='Cones':
			self.spincones = vtk.vtkConeSource()
			self.spincones.SetRadius(0.50)
			self.spincones.SetHeight(1.00)
			self.spincones.SetResolution(10)
			ASDMomActors.SpinMapper.SetSourceConnection(self.spincones.GetOutputPort())
			ASDMomActors.SpinMapper.OrientOn()
			renWin.Render()
