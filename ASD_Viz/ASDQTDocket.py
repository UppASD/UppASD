#!/usr/bin/env vtkpython
################################################################################
# CLASS: ASDQTDockWindow
# @author Jonathan Chico (08/09/2017)
# @description
# Wrapper class that contains the needed information to generate the docket window present
# in the main window. It contains all the connection calls to modify the visuzalization,
# it has the needed actions to modify the visualization types, options and the
# camera.
################################################################################

import numpy as np
import ASDVTKReading
import ASDMomVTKActors
import ASDVTKVizOptions
import ASDNeighVTKActors
from PyQt4 import QtCore, QtGui

################################################################################
# Create the Docket widget to have the options for the visualization
# This will give options for the moment visualization and the nieghbours
################################################################################
class ASDQTDockWindow():
	""" The docket widget. This large area is setup on the right hand size of
	the main window, it is defined as a scrollable area so that new functionality
	can be added without problems to the size in the y-direction.
	Currently it includes the controls for the following options
	- Moment mode
    	- Scalable spins
    	- Magnetization density with a delaunay2D grid
    	- Spin directions (useful for the Magnetization density)
    	- Plane clipping with options between 3 different planes
    	- Choice of projection of the magnetization
    	- Contours from the magnetization density
    	- Axes widget for orientation
    	- Toggling of the actors
    	- Changing the color map used for the visualization
	- Neighbour mode
    	- Toggling of the atoms and neighbour cloud for neighbour mapping"""

	MomActors=ASDMomVTKActors.ASDMomActors()
	NeighActors=ASDNeighVTKActors.ASDNeighActors()

	def createASDDock(self,ASDMainWindow):
		self.ASDVizOptions=ASDVTKVizOptions.ASDVizOptions()

		camera_yaw       = 0
		camera_roll      = 0
		camera_pitch     = 0
		camera_azimuth   = 0
		camera_elevation = 0
		camera_pos       = np.zeros(3,dtype=np.float32)
		camera_focal     = np.zeros(3,dtype=np.float32)

		# Defining the font size, to avoid problems in other environments
		font=QtGui.QFont()
		font.setPointSize(12)
		# Creating the side dock widget
		VizDock = QtGui.QDockWidget('Viz. Opt.',ASDMainWindow)
		VizDock.setFont(font)
		VizDock.setMinimumSize(300, 500)
		VizDock.setAllowedAreas(QtCore.Qt.RightDockWidgetArea)

		Camera_Dock = QtGui.QDockWidget('Camera Opt.',ASDMainWindow)
		Camera_Dock.setFont(font)
		Camera_Dock.setMinimumSize(300, 500)
		Camera_Dock.setMaximumWidth(300)
		Camera_Dock.setAllowedAreas(QtCore.Qt.RightDockWidgetArea)

	    ########################################################################
	    # Vizualization data Dock
	    ########################################################################
		MainWidgetDock = QtGui.QWidget()
		DockLayout = QtGui.QVBoxLayout()
		MainWidgetDock.setLayout(DockLayout)
		########################################################################
		# Options for the Moments options in the Dock Widget
		########################################################################
		# Add a Group mode to encompass the options, this will be turned on if
		# the restart or moments mode are turned on
		self.MomentBox= QtGui.QGroupBox("Moments Mode")
		self.MomentBox.setCheckable(True)
		self.MomentBox.setChecked(False)
		self.MomentBox.setFont(font)

		# Set the layout for the moments visualization
		MomBoxLayout = QtGui.QVBoxLayout()
		self.MomentBox.setLayout(MomBoxLayout)

		# Add checkbox to add the spins
		self.SpinsGroup = QtGui.QGroupBox("Spins")
		self.SpinsGroup.setCheckable(True)
		self.SpinsGroup.setChecked(False)
		self.SpinsGroup.toggled.connect(self.ASDVizOptions.toggle_spins)
		self.SpinsGroup.setFont(font)

		SpinsBox=QtGui.QVBoxLayout()
		SpinButtonGroup=QtGui.QButtonGroup()
		self.spins_x=QtGui.QRadioButton("X Proj.")
		self.spins_x.toggled.connect(self.ASDVizOptions.set_spins_color_x)
		self.spins_x.setFont(font)
		self.spins_y=QtGui.QRadioButton("Y Proj.")
		self.spins_y.toggled.connect(self.ASDVizOptions.set_spins_color_y)
		self.spins_y.setFont(font)
		self.spins_z=QtGui.QRadioButton("Z Proj.")
		self.spins_z.setChecked(True)
		self.spins_z.toggled.connect(self.ASDVizOptions.set_spins_color_z)
		self.spins_z.setFont(font)

		# Definition or the slider for the size of the moments
		self.SpinSizeSL = QtGui.QSlider(QtCore.Qt.Horizontal)
		self.SpinSizeSL.Spinlabel = QtGui.QLabel()
		self.SpinSizeSL.Spinlabel.setText('Spin size')
		self.SpinSizeSL.Spinlabel.setFont(font)
		self.SpinSizeSL.setMinimum(10)
		self.SpinSizeSL.setMaximum(40)
		self.SpinSizeSL.setValue(10)
		self.SpinSizeSL.setTickPosition(QtGui.QSlider.TicksBelow)
		self.SpinSizeSL.setTickInterval(5)
		self.SpinSizeSL.valueChanged.connect(self.ASDVizOptions.ChangeSpinsSize)

		# Adding the widgets to define the color projection used
		SpinsBox.addWidget(self.spins_x)
		SpinsBox.addWidget(self.spins_y)
		SpinsBox.addWidget(self.spins_z)
		SpinsBox.addWidget(self.SpinSizeSL.Spinlabel)
		SpinsBox.addWidget(self.SpinSizeSL)
		self.SpinsGroup.setLayout(SpinsBox)

		MomBoxLayout.addWidget(self.SpinsGroup)

		# Add group box for the options regarding the magnetization density
		self.DensityGroup=QtGui.QGroupBox("Mag. Density")
		self.DensityGroup.setCheckable(True)
		self.DensityGroup.setChecked(True)
		self.DensityGroup.toggled.connect(self.ASDVizOptions.toggle_density)
		self.DensityGroup.setFont(font)

		# Options for the magnetization density projections
		DensBox=QtGui.QVBoxLayout()
		buttonGroup=QtGui.QButtonGroup()
		self.dens_x=QtGui.QRadioButton("X Proj.")
		self.dens_x.setFont(font)
		self.dens_x.toggled.connect(self.ASDVizOptions.set_color_x)
		self.dens_y=QtGui.QRadioButton("Y Proj.")
		self.dens_y.setFont(font)
		self.dens_y.toggled.connect(self.ASDVizOptions.set_color_y)
		self.dens_z=QtGui.QRadioButton("Z Proj.")
		self.dens_z.setFont(font)
		self.dens_z.setChecked(True)
		self.dens_z.toggled.connect(self.ASDVizOptions.set_color_z)

		# Add the widgets to choose the color projection for the density
		DensBox.addWidget(self.dens_x)
		DensBox.addWidget(self.dens_y)
		DensBox.addWidget(self.dens_z)
		self.DensityGroup.setLayout(DensBox)

		MomBoxLayout.addWidget(self.DensityGroup)

		# Add group box for the options regarding the magnetization density
		self.ClipperGroup=QtGui.QGroupBox("Clipping Planes")
		self.ClipperGroup.setFont(font)
		self.ClipperGroup.setCheckable(True)
		self.ClipperGroup.setChecked(False)
		self.ClipperGroup.toggled.connect(self.ASDVizOptions.toggle_clipper)
		ClipperBox=QtGui.QVBoxLayout()

		# Add Button sets for different planes
		PlaneButtonGroup=QtGui.QButtonGroup()
		self.plane_x=QtGui.QRadioButton("(1,0,0) Normal")
		self.plane_x.setFont(font)
		self.plane_x.toggled.connect(self.set_plane_x)
		self.plane_y=QtGui.QRadioButton("(0,1,0) Normal")
		self.plane_y.setFont(font)
		self.plane_y.toggled.connect(self.set_plane_y)
		self.plane_z=QtGui.QRadioButton("(0,0,1) Normal")
		self.plane_z.setFont(font)
		self.plane_z.toggled.connect(self.set_plane_z)

		# Slider for the plane clipper
		self.ClipperSL = QtGui.QSlider(QtCore.Qt.Horizontal)
		self.ClipperSL.Cliplabel = QtGui.QLabel()
		self.ClipperSL.setMinimum(0)
		self.ClipperSL.setValue(0)
		self.ClipperSL.setTickPosition(QtGui.QSlider.TicksBelow)
		self.ClipperSL.setTickInterval(1)
		self.ClipperSL.Cliplabel.setText('Clip. Plane Pos.=(0.0,0.0,0.0)')
		self.ClipperSL.Cliplabel.setFont(font)
		self.ClipperSL.valueChanged.connect(self.ClippingUpdate)

		# Adding the buttons and slider to the box
		ClipperBox.addWidget(self.plane_x)
		ClipperBox.addWidget(self.plane_y)
		ClipperBox.addWidget(self.plane_z)
		ClipperBox.addWidget(self.ClipperSL)
		ClipperBox.addWidget(self.ClipperSL.Cliplabel)

		# Set the box layout for the clipper
		self.ClipperGroup.setLayout(ClipperBox)
		MomBoxLayout.addWidget(self.ClipperGroup)

		# Creating a box for the rest of the options of the moment visualization
		MomOptGroup = QtGui.QGroupBox("Mom. Opt.")
		MomOptGroup.setFont(font)

		# Creating the layout for the rest of the options
		MomOptBox=QtGui.QVBoxLayout()

		# Add checkbox for the arrows with the magnetization direction (mostly to be used with the density)
		self.dir_check = QtGui.QCheckBox("Directions")
		self.dir_check.setFont(font)
		self.dir_check.toggle()
		self.dir_check.setChecked(False)
		self.dir_check.stateChanged.connect(self.ASDVizOptions.toggle_directions)

		# Adding the Checkbox for the Contours
		self.cont_check = QtGui.QCheckBox("Contours")
		self.cont_check.setFont(font)
		self.cont_check.toggle()
		self.cont_check.stateChanged.connect(self.ASDVizOptions.toggle_contours)

		# Adding the checkbox for the axes widget
		self.axes_check = QtGui.QCheckBox("Axes")
		self.axes_check.setFont(font)
		self.axes_check.toggle()
		self.axes_check.stateChanged.connect(self.ASDVizOptions.toggle_Axes)

		# Adding the checkbox for the axes widget
		self.ScalarBar_check = QtGui.QCheckBox("Scalar Bar")
		self.ScalarBar_check.setFont(font)
		self.ScalarBar_check.toggle()
		self.ScalarBar_check.stateChanged.connect(self.ASDVizOptions.toggle_ScalarBar)
		# Adding a checkbox for the cluster if present
		self.cluster_check = QtGui.QCheckBox("Cluster")
		self.cluster_check.setFont(font)
		self.cluster_check.toggle()
		self.cluster_check.stateChanged.connect(self.ASDVizOptions.toggle_cluster)
		self.cluster_check.setVisible(False)
		# Adding a checkbox for the KMC particles if present
		self.KMC_check = QtGui.QCheckBox("KMC particles")
		self.KMC_check.setFont(font)
		self.KMC_check.toggle()
		self.KMC_check.stateChanged.connect(self.ASDVizOptions.toggle_KMC)
		self.KMC_check.setVisible(False)

		# Adding radio buttons to select the different color scales
		ColorbuttonGroup=QtGui.QButtonGroup()
		self.color_map_Coolwarm=QtGui.QRadioButton("Coolwarm Color Map")
		self.color_map_Coolwarm.setChecked(True)
		self.color_map_Coolwarm.setFont(font)
		self.color_map_Coolwarm.toggled.connect(self.ASDVizOptions.set_Coolwarm_lut)
		self.color_map_BlackBody=QtGui.QRadioButton("Black Body Color Map")
		self.color_map_BlackBody.setFont(font)
		self.color_map_BlackBody.toggled.connect(self.ASDVizOptions.set_BlackBody_lut)
		self.color_map_RdGy=QtGui.QRadioButton("RdGy Color Map")
		self.color_map_RdGy.setFont(font)
		self.color_map_RdGy.toggled.connect(self.ASDVizOptions.set_RdGy_lut)
		self.color_map_Spectral=QtGui.QRadioButton("Spectral Color Map")
		self.color_map_Spectral.setFont(font)
		self.color_map_Spectral.toggled.connect(self.ASDVizOptions.set_Spectral_lut)

		# Adding the widgets to the created box
		MomOptBox.addWidget(self.dir_check)
		MomOptBox.addWidget(self.cont_check)
		MomOptBox.addWidget(self.axes_check)
		MomOptBox.addWidget(self.cluster_check)
		MomOptBox.addWidget(self.KMC_check)
		MomOptBox.addWidget(self.ScalarBar_check)
		MomOptBox.addWidget(self.color_map_Coolwarm)
		MomOptBox.addWidget(self.color_map_BlackBody)
		MomOptBox.addWidget(self.color_map_RdGy)
		MomOptBox.addWidget(self.color_map_Spectral)

		# Setting the layout for the options for the moments options
		MomOptGroup.setLayout(MomOptBox)
		MomBoxLayout.addWidget(MomOptGroup)

		########################################################################
		# Options for the Neighbour options in the Dock Widget
		########################################################################
		# Now the visualization for the neighbour atoms needs to be setup
		self.NeighBox = QtGui.QGroupBox("Neighbour Mode")
		self.NeighBox.setCheckable(True)
		self.NeighBox.setChecked(False)
		self.NeighBox.setFont(font)

		# Set the layout for the neighbour visualization
		NeighBoxLayout = QtGui.QGridLayout()
		self.NeighBox.setLayout(NeighBoxLayout)
		self.atom_viz = QtGui.QCheckBox("Atoms")
		self.atom_viz.setFont(font)
		self.atom_viz.toggle()
		self.atom_viz.stateChanged.connect(self.ASDVizOptions.toggle_NAtoms)
		NeighBoxLayout.addWidget(self.atom_viz,0,0,1,1)
		self.neigh_viz = QtGui.QCheckBox("Neighbour cloud")
		self.neigh_viz.setFont(font)
		self.neigh_viz.toggle()
		self.neigh_viz.stateChanged.connect(self.ASDVizOptions.toggle_Neigh)
		NeighBoxLayout.addWidget(self.neigh_viz,1,0,1,1)

		# Add the slider to select the neighbour that is being visualized
		self.NeighSL = QtGui.QSlider(QtCore.Qt.Horizontal)
		self.NeighSL.setMinimum(0)
		self.NeighSL.setValue(0)
		self.NeighSL.setTickPosition(QtGui.QSlider.TicksBelow)
		self.NeighSL.setTickInterval(1)
		self.NeighSL.valueChanged.connect(self.UpdateSlider)

		self.NeighSL.Nlabel = QtGui.QLabel()
		self.NeighSL.Nlabel.setText('Number of neighbours=0')
		self.NeighSL.Nlabel.setFont(font)

		self.NeighSL.SetAtom = QtGui.QLineEdit()
		self.NeighSL.SetAtom.setText(QtCore.QString.number(self.NeighSL.value()+1))
		self.NeighSL.SetAtom.editingFinished.connect(self.UpdateSlider)
		self.NeighSL.ALabel = QtGui.QLabel()
		self.NeighSL.ALabel.setText('Atom number=')

		NeighBoxLayout.addWidget(self.NeighSL.ALabel,2,0,1,1)
		NeighBoxLayout.addWidget(self.NeighSL.SetAtom,2,1,1,1)
		NeighBoxLayout.addWidget(self.NeighSL,3,0,1,2)
		NeighBoxLayout.addWidget(self.NeighSL.Nlabel,4,0,1,2)

		########################################################################
		# Add the main widgets to the dockets
		########################################################################
		# Add spacer item to ensure a correct layout
		spacerItem = QtGui.QSpacerItem(20, 40, QtGui.QSizePolicy.Minimum, QtGui.QSizePolicy.Expanding)
		DockLayout.addWidget(self.MomentBox)
		DockLayout.addWidget(self.NeighBox)
		DockLayout.addItem(spacerItem)

		########################################################################
		# Finish of Viz Options Dock
		########################################################################
		########################################################################
		# Camera options widget
		########################################################################
		CameraWidgetDock=QtGui.QWidget()
		CameraDockLayout=QtGui.QVBoxLayout()
		CameraWidgetDock.setLayout(CameraDockLayout)

		Screen_box_setup = QtGui.QGroupBox()
		Screen_box_setup.setAlignment(QtCore.Qt.AlignCenter)
		Screen_box_setup.setCheckable(False)
		Screen_box_setup.setChecked(True)

		Screen_box_layout = QtGui.QVBoxLayout()

		########################################################################
		# Camera angles setups
		########################################################################
		Camera_setup_box = QtGui.QGroupBox('Camera Angles')
		Camera_setup_box.setAlignment(QtCore.Qt.AlignCenter)
		Camera_setup_box.setCheckable(False)
		Camera_setup_box.setChecked(True)

		Camera_box_layout = QtGui.QGridLayout()

		set_elevation_label=QtGui.QLabel()
		set_elevation_label.setText("Elevation")
		self.set_elevation=QtGui.QLineEdit()

		set_azimuth_label=QtGui.QLabel()
		set_azimuth_label.setText("Azimuth")
		self.set_azimuth=QtGui.QLineEdit()

		set_roll_label=QtGui.QLabel()
		set_roll_label.setText("Roll")
		self.set_roll=QtGui.QLineEdit()

		set_pitch_label=QtGui.QLabel()
		set_pitch_label.setText("Pitch")
		self.set_pitch=QtGui.QLineEdit()

		set_yaw_label=QtGui.QLabel()
		set_yaw_label.setText("Yaw")
		self.set_yaw=QtGui.QLineEdit()

		Camera_box_layout.addWidget(set_elevation_label,0,0,1,1)
		Camera_box_layout.addWidget(self.set_elevation,1,0,1,1)
		Camera_box_layout.addWidget(set_azimuth_label,0,1,1,1)
		Camera_box_layout.addWidget(self.set_azimuth,1,1,1,1)
		Camera_box_layout.addWidget(set_roll_label,2,0,1,1)
		Camera_box_layout.addWidget(self.set_roll,3,0,1,1)
		Camera_box_layout.addWidget(set_pitch_label,2,1,1,1)
		Camera_box_layout.addWidget(self.set_pitch,3,1,1,1)
		Camera_box_layout.addWidget(set_yaw_label,2,2,1,1)
		Camera_box_layout.addWidget(self.set_yaw,3,2,1,1)

		Camera_setup_box.setLayout(Camera_box_layout)

		########################################################################
		# Setup the camera positions
		########################################################################
		Camera_pos_box = QtGui.QGroupBox('Camera Position')
		Camera_pos_box.setAlignment(QtCore.Qt.AlignCenter)
		Camera_pos_box.setCheckable(False)
		Camera_pos_box.setChecked(True)

		Camera_pos_layout = QtGui.QGridLayout()
		set_label_focal=QtGui.QLabel()
		set_label_focal.setText("Set Focal Point")
		set_label_focal_x=QtGui.QLabel()
		set_label_focal_x.setText("X")
		self.focal_x=QtGui.QLineEdit()
		set_label_focal_y=QtGui.QLabel()
		set_label_focal_y.setText("Y")
		self.focal_y=QtGui.QLineEdit()
		set_label_focal_z=QtGui.QLabel()
		set_label_focal_z.setText("Z")
		self.focal_z=QtGui.QLineEdit()

		set_label_postion=QtGui.QLabel()
		set_label_postion.setText("Set Camera Position")
		set_label_pos_x=QtGui.QLabel()
		set_label_pos_x.setText("X")
		self.pos_x=QtGui.QLineEdit()
		set_label_pos_y=QtGui.QLabel()
		set_label_pos_y.setText("Y")
		self.pos_y=QtGui.QLineEdit()
		set_label_pos_z=QtGui.QLabel()
		set_label_pos_z.setText("Z")
		self.pos_z=QtGui.QLineEdit()

		Camera_pos_layout.addWidget(set_label_focal,0,0,1,3)
		Camera_pos_layout.addWidget(set_label_focal_x,1,0,1,1)
		Camera_pos_layout.addWidget(self.focal_x,2,0,1,1)
		Camera_pos_layout.addWidget(set_label_focal_y,1,1,1,1)
		Camera_pos_layout.addWidget(self.focal_y,2,1,1,1)
		Camera_pos_layout.addWidget(set_label_focal_z,1,2,1,1)
		Camera_pos_layout.addWidget(self.focal_z,2,2,1,1)

		Camera_pos_layout.addWidget(set_label_postion,3,0,1,3)
		Camera_pos_layout.addWidget(set_label_pos_x,4,0,1,1)
		Camera_pos_layout.addWidget(self.pos_x,5,0,1,1)
		Camera_pos_layout.addWidget(set_label_pos_y,4,1,1,1)
		Camera_pos_layout.addWidget(self.pos_y,5,1,1,1)
		Camera_pos_layout.addWidget(set_label_pos_z,4,2,1,1)
		Camera_pos_layout.addWidget(self.pos_z,5,2,1,1)

		Camera_pos_box.setLayout(Camera_pos_layout)

		self.Projection_box = QtGui.QGroupBox('Parallel Projection')
		self.Projection_box.setAlignment(QtCore.Qt.AlignCenter)
		self.Projection_box.setCheckable(True)
		self.Projection_box.setChecked(False)
		self.Projection_box.toggled.connect(ASDMainWindow.CameraHandling)

		Projection_layout=QtGui.QGridLayout()
		self.proj_sl = QtGui.QSlider(QtCore.Qt.Horizontal)
		self.proj_sl.Proj_label = QtGui.QLabel()
		self.proj_sl.Proj_label.setText('Parallel Projection')
		self.proj_sl.Proj_label.setFont(font)
		self.proj_sl.setMinimum(0)
		self.proj_sl.setMaximum(100)
		self.proj_sl.setValue(10)
		self.proj_sl.setTickPosition(QtGui.QSlider.TicksBelow)
		self.proj_sl.setTickInterval(5)
		self.proj_sl.valueChanged.connect(ASDMainWindow.CameraHandling)

		proj_sl_label_text=QtGui.QLabel()
		proj_sl_label_text.setText("Parallel Scale")
		self.proj_sl_label=QtGui.QLineEdit()
		self.proj_sl_label.editingFinished.connect(ASDMainWindow.CameraHandling)

		Projection_layout.addWidget(self.proj_sl,0,0,1,2)
		Projection_layout.addWidget(proj_sl_label_text,1,0,1,1)
		Projection_layout.addWidget(self.proj_sl_label,1,1,1,1)
		self.Projection_box.setLayout(Projection_layout)

		########################################################################
		# Camera projection options
		########################################################################

		Camera_proj_box=QtGui.QGroupBox("Set View")
		Camera_proj_box.setFont(font)
		Camera_proj_box.setCheckable(False)
		Camera_proj_box.setChecked(True)
		Camera_proj_layout=QtGui.QHBoxLayout()

		# Add Button sets for different planes
		self.Camera_proj_x=QtGui.QPushButton("(1,0,0)")
		self.Camera_proj_x.clicked.connect(ASDMainWindow.CameraHandling)
		self.Camera_proj_y=QtGui.QPushButton("(0,1,0)")
		self.Camera_proj_y.clicked.connect(ASDMainWindow.CameraHandling)
		self.Camera_proj_z=QtGui.QPushButton("(0,0,1)")
		self.Camera_proj_z.clicked.connect(ASDMainWindow.CameraHandling)

		Camera_proj_layout.addWidget(self.Camera_proj_x)
		Camera_proj_layout.addWidget(self.Camera_proj_y)
		Camera_proj_layout.addWidget(self.Camera_proj_z)

		Camera_proj_box.setLayout(Camera_proj_layout)

		########################################################################
		# Creating the needed buttons to accept and reset
		########################################################################
		self.accept_button=QtGui.QPushButton("Set Camera")
		self.accept_button.clicked.connect(ASDMainWindow.CameraHandling)

		self.reset_button=QtGui.QPushButton("Reset Camera")
		self.reset_button.clicked.connect(ASDMainWindow.CameraHandling)
		########################################################################
		# Finally adding the widgets
		########################################################################
		Screen_box_layout.addWidget(Camera_setup_box)
		Screen_box_layout.addWidget(Camera_pos_box)
		Screen_box_setup.setLayout(Screen_box_layout)

		CameraDockLayout.addWidget(Screen_box_setup)
		CameraDockLayout.addWidget(Camera_proj_box)
		CameraDockLayout.addWidget(self.Projection_box)
		CameraDockLayout.addWidget(self.accept_button)
		CameraDockLayout.addWidget(self.reset_button)
		########################################################################
		# Finally adding the widgets to the docket
		########################################################################

		VizDock.setWidget(MainWidgetDock)
		Camera_Dock.setWidget(CameraWidgetDock)
		ASDMainWindow.addDockWidget(QtCore.Qt.RightDockWidgetArea, VizDock)
		ASDMainWindow.addDockWidget(QtCore.Qt.RightDockWidgetArea, Camera_Dock)
		ASDMainWindow.tabifyDockWidget(VizDock,Camera_Dock)
		VizDock.raise_()
		# Setup for the scroll area
		scrollArea = QtGui.QScrollArea(VizDock)
		scrollArea.setViewportMargins(0, 30, 0, 0)
		scrollArea.setMinimumSize(300,600)
		scrollArea.setWidgetResizable(True)
		scrollArea.setVerticalScrollBarPolicy(QtCore.Qt.ScrollBarAlwaysOn)
		scrollArea.setHorizontalScrollBarPolicy(QtCore.Qt.ScrollBarAlwaysOff)
		scrollArea.setSizePolicy(QtGui.QSizePolicy.Minimum, QtGui.QSizePolicy.Expanding)
		scrollArea.setWidget(MainWidgetDock)

		Camera_scrollArea = QtGui.QScrollArea(Camera_Dock)
		Camera_scrollArea.setViewportMargins(0, 30, 0, 0)
		Camera_scrollArea.setMinimumSize(300,600)
		Camera_scrollArea.setWidgetResizable(True)
		Camera_scrollArea.setVerticalScrollBarPolicy(QtCore.Qt.ScrollBarAlwaysOn)
		Camera_scrollArea.setHorizontalScrollBarPolicy(QtCore.Qt.ScrollBarAlwaysOff)
		Camera_scrollArea.setSizePolicy(QtGui.QSizePolicy.Minimum, QtGui.QSizePolicy.Expanding)
		Camera_scrollArea.setWidget(CameraWidgetDock)

	############################################################################
	# Set the clipping plane such that the normal plane is (1,0,0)
	############################################################################
	def set_plane_x (self,check):
	    if check:
	        ASDQTDockWindow.MomActors.plane.SetOrigin(ASDQTDockWindow.MomActors.xmin,ASDQTDockWindow.MomActors.ymid,ASDQTDockWindow.MomActors.zmid)
	        ASDQTDockWindow.MomActors.plane.SetNormal(1,0,0)
	        self.ClipperSL.setMinimum(ASDQTDockWindow.MomActors.xmin)
	        self.ClipperSL.setMaximum(ASDQTDockWindow.MomActors.xmax)
	        self.ClipperSL.setValue(ASDQTDockWindow.MomActors.xmin)

	############################################################################
	# Set the clipping plane such that the normal plane is (0,1,0)
	############################################################################
	def set_plane_y (self,check):
	    if check:
	        ASDQTDockWindow.MomActors.plane.SetOrigin(ASDQTDockWindow.MomActors.xmid,ASDQTDockWindow.MomActors.ymin,ASDQTDockWindow.MomActors.zmid)
	        ASDQTDockWindow.MomActors.plane.SetNormal(0,1,0)
	        self.ClipperSL.setMinimum(ASDQTDockWindow.MomActors.ymin)
	        self.ClipperSL.setMaximum(ASDQTDockWindow.MomActors.ymax)
	        self.ClipperSL.setValue(ASDQTDockWindow.MomActors.ymin)

	############################################################################
	# Set the clipping plane such that the normal plane is (0,0,1)
	############################################################################
	def set_plane_z (self,check):
	    if check:
	        ASDQTDockWindow.MomActors.plane.SetOrigin(ASDQTDockWindow.MomActors.xmid,ASDQTDockWindow.MomActors.ymid,ASDQTDockWindow.MomActors.zmin)
	        ASDQTDockWindow.MomActors.plane.SetNormal(0,0,1)
	        self.ClipperSL.setMinimum(ASDQTDockWindow.MomActors.zmin)
	        self.ClipperSL.setMaximum(ASDQTDockWindow.MomActors.zmax)
	        self.ClipperSL.setValue(ASDQTDockWindow.MomActors.zmin)

	############################################################################
	# Set the position of the clipping plane via the slider
	############################################################################
	def ClippingUpdate(self,value):
		if self.plane_x.isChecked():
			ASDQTDockWindow.MomActors.plane.SetOrigin(value, ASDQTDockWindow.MomActors.ymid, ASDQTDockWindow.MomActors.zmid)
			self.ClipperSL.Cliplabel.setText('Clip. Plane Pos.={:.1f},{:.1f},{:.1f}'.format(value,0,0))
		elif self.plane_y.isChecked():
			ASDQTDockWindow.MomActors.plane.SetOrigin(ASDQTDockWindow.MomActors.xmid, value, ASDQTDockWindow.MomActors.zmid)
			self.ClipperSL.Cliplabel.setText('Clip. Plane Pos.={:.1f},{:.1f},{:.1f}'.format(0,value,0))
		elif self.plane_z.isChecked():
			ASDQTDockWindow.MomActors.plane.SetOrigin(ASDQTDockWindow.MomActors.xmid,ASDQTDockWindow.MomActors.ymid,value)
			self.ClipperSL.Cliplabel.setText('Clip. Plane Pos.={:.1f},{:.1f},{:.1f}'.format(0,0,value))

	############################################################################
	# Update the slider parameters when either the text or slider position ware changed
	############################################################################
	def UpdateSlider(self):
		slid_value = self.NeighSL.value()
		line_value = self.NeighSL.SetAtom.text().toInt()[0]
		if slid_value!=line_value:
			if self.NeighSL.SetAtom.isModified():
				# hence one needs to add a connector to that
				ASD_data=ASDVTKReading.ASDReading()
				ASD_data.ReadingWrapper(1,'N')
				ASD_data.neighs, ASD_data.nTypes = ASD_data.setNeighbours(ASD_data.neighbours,line_value-1,ASD_data.coord)
				ASDQTDockWindow.NeighActors.NeighGrid.SetPoints(ASD_data.neighs)
				ASDQTDockWindow.NeighActors.NeighGrid.GetPointData().SetScalars(ASD_data.nTypes)
				ASDQTDockWindow.NeighActors.Neighs.Update()
				self.NeighSL.Nlabel.setText('Number of neighbours={:}'.format(ASD_data.neighbours[line_value-1][0]))
				self.NeighSL.setValue(line_value)
				self.NeighSL.SetAtom.setModified(False)
			else:
				ASD_data=ASDVTKReading.ASDReading()
				ASD_data.ReadingWrapper(1,'N')
				(ASD_data.neighs, ASD_data.nTypes) = ASD_data.setNeighbours(ASD_data.neighbours,slid_value,ASD_data.coord)
				ASDQTDockWindow.NeighActors.NeighGrid.SetPoints(ASD_data.neighs)
				ASDQTDockWindow.NeighActors.NeighGrid.GetPointData().SetScalars(ASD_data.nTypes)
				ASDQTDockWindow.NeighActors.Neighs.Update()
				self.NeighSL.Nlabel.setText('Number of neighbours={:}'.format(ASD_data.neighbours[slid_value][0]))
				self.NeighSL.SetAtom.setText(QtCore.QString.number(slid_value+1))


	############################################################################
	# Function definition to get the data from the updated Camera_Dock
	############################################################################
	def getCameraData(self):

		ASDQTDockWindow.MomActors.camera_yaw=float(self.set_yaw.text())
		ASDQTDockWindow.MomActors.camera_roll=float(self.set_roll.text())
		ASDQTDockWindow.MomActors.camera_pitch=float(self.set_pitch.text())
		ASDQTDockWindow.MomActors.camera_azimuth=float(self.set_azimuth.text())
		ASDQTDockWindow.MomActors.camera_elevation=float(self.set_elevation.text())
		ASDQTDockWindow.MomActors.camera_pos[0]=float(self.pos_x.text())
		ASDQTDockWindow.MomActors.camera_pos[1]=float(self.pos_y.text())
		ASDQTDockWindow.MomActors.camera_pos[2]=float(self.pos_z.text())
		ASDQTDockWindow.MomActors.camera_focal[0]=float(self.focal_x.text())
		ASDQTDockWindow.MomActors.camera_focal[1]=float(self.focal_y.text())
		ASDQTDockWindow.MomActors.camera_focal[2]=float(self.focal_z.text())

	############################################################################
	# Function to force the camera to update with new user defined values
	############################################################################
	def update_Camera(self,ren,renWin):
		ren.GetActiveCamera().SetFocalPoint(ASDQTDockWindow.MomActors.camera_focal)
		ren.GetActiveCamera().SetPosition(ASDQTDockWindow.MomActors.camera_pos)
		ren.GetActiveCamera().Elevation(ASDQTDockWindow.MomActors.camera_elevation)
		ren.GetActiveCamera().Azimuth(ASDQTDockWindow.MomActors.camera_azimuth)
		ren.GetActiveCamera().Pitch(ASDQTDockWindow.MomActors.camera_pitch)
		ren.GetActiveCamera().Roll(ASDQTDockWindow.MomActors.camera_roll)
		ren.GetActiveCamera().Yaw(ASDQTDockWindow.MomActors.camera_yaw)
		renWin.Render()
		return

	############################################################################
	# Function to update the values diaplayed in the Camera_Dock
	############################################################################
	def update_dock_info(self):
		self.focal_x.setText(QtCore.QString("%2").arg(float(ASDQTDockWindow.MomActors.camera_focal[0]),1,"f",2))
		self.focal_y.setText(QtCore.QString("%2").arg(float(ASDQTDockWindow.MomActors.camera_focal[1]),1,"f",2))
		self.focal_z.setText(QtCore.QString("%2").arg(float(ASDQTDockWindow.MomActors.camera_focal[2]),1,"f",2))
		self.pos_x.setText(QtCore.QString("%2").arg(float(ASDQTDockWindow.MomActors.camera_pos[0]),1,"f",2))
		self.pos_y.setText(QtCore.QString("%2").arg(float(ASDQTDockWindow.MomActors.camera_pos[1]),1,"f",2))
		self.pos_z.setText(QtCore.QString("%2").arg(float(ASDQTDockWindow.MomActors.camera_pos[2]),1,"f",2))
		self.set_yaw.setText(QtCore.QString("%2").arg(float(ASDQTDockWindow.MomActors.camera_yaw),1,"f",2))
		self.set_roll.setText(QtCore.QString("%2").arg(float(ASDQTDockWindow.MomActors.camera_roll),1,"f",2))
		self.set_pitch.setText(QtCore.QString("%2").arg(float(ASDQTDockWindow.MomActors.camera_pitch),1,"f",2))
		self.set_azimuth.setText(QtCore.QString("%2").arg(float(ASDQTDockWindow.MomActors.camera_azimuth),1,"f",2))
		self.set_elevation.setText(QtCore.QString("%2").arg(float(ASDQTDockWindow.MomActors.camera_elevation),1,"f",2))

		return

	############################################################################
	# Toggle the parallel projection for the camera
	############################################################################
	def toggle_projections(self,ren):
	        self.Projection_box.setChecked(True)
	        self.proj_sl_label.setText(QtCore.QString("%2").arg(float(self.proj_sl.value()),1,"f",2))

		return
	############################################################################
	# Toggle the parallel projection for the camera
	############################################################################
	def ChangeParallelProj(self,ren,renWin,line,slider):
	    if line:
			self.proj_sl.setValue(float(self.proj_sl_label.text()))
			ren.GetActiveCamera().SetParallelScale(float(self.proj_sl_label.text()))
			renWin.Render()
	    if slider:
			ren.GetActiveCamera().SetParallelScale(self.proj_sl.value())
			self.proj_sl_label.setText(QtCore.QString("%2").arg(float(self.proj_sl.value()),1,"f",2))
			renWin.Render()
