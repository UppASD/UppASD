#!/usr/bin/env vtkpython
################################################################################
# CLASS: ASDVizMenu
# @author Jonathan Chico (08/09/2017)
# @description
# Wrapper class that contains the needed information to generate the menus present
# in the main window. It contains all the connection calls to modify the visuzalization
# mirroring much of the options that appear in the Docket window.
################################################################################
import vtk
import sys
import time
import glob
import os.path
import ASDVTKReading
import ASDMomVTKActors
import ASDVTKVizOptions
from PyQt5.QtWidgets import QApplication
from PyQt5 import QtCore, QtWidgets

class ASDVizMenu():

    ASDOptions=ASDVTKVizOptions.ASDVizOptions()
    ASD_data=ASDVTKReading.ASDReading()

    def ASDMenuCreator(self,ASDMainWindow):

        menubar= ASDMainWindow.menuBar()
        menubar.setNativeMenuBar(False)

        self.exit = QtWidgets.QAction('Exit', ASDMainWindow)
        self.exit.setShortcut('Ctrl+Q')
        self.exit.setStatusTip('Exit application')
        self.exit.setCheckable(True)
        self.exit.setChecked(False)
        self.exit.toggled.connect(self.close_app)

        # The QtActionGroup should ensure  that the modes are mutually exclusive
        # However, when using it it seems that the actions are read twice ..
        self.mode_group=QtWidgets.QActionGroup(ASDMainWindow)
        self.mode_group.setExclusive(True)

        self.restart = QtWidgets.QAction('Restart',ASDMainWindow)
        self.restart.setShortcut('Ctrl+R')
        self.restart.setStatusTip('Uses the restart file')
        self.restart.setCheckable(True)
        self.restart.setChecked(False)
        self.restart.toggled.connect(ASDMainWindow.AddActors)

        self.moments = QtWidgets.QAction('Moments', ASDMainWindow)
        self.moments.setShortcut('Ctrl+M')
        self.moments.setStatusTip('Uses the moment file (Movie)')
        self.moments.setCheckable(True)
        self.moments.setChecked(False)
        self.moments.toggled.connect(ASDMainWindow.AddActors)

        # Add the option to visualize the neighbours
        self.neighbours = QtWidgets.QAction('Neighbours',ASDMainWindow)
        self.neighbours.setShortcut('Crtl+N')
        self.neighbours.setStatusTip('Uses the struct file to visualize the neighbour map')
        self.neighbours.setCheckable(True)
        self.neighbours.setChecked(False)
        self.neighbours.toggled.connect(ASDMainWindow.AddActors)

        self.mode_group.addAction(self.restart)
        self.mode_group.addAction(self.moments)
        self.mode_group.addAction(self.neighbours)

        # Creates a QtWidgets action that will be used to take an snapshot
        self.snap = QtWidgets.QAction('Snapshot', ASDMainWindow)
        self.snap.setShortcut('Ctrl+P')
        self.snap.setStatusTip('Snapshot of current spin configuration')
        # Adding the Play button to visualize a movie made from the moments
        self.PlayButton = QtWidgets.QToolButton()
        self.PlayButton.setCheckable(True)
        self.PlayButton.setChecked(False)
        self.PlayButton.setArrowType(QtCore.Qt.RightArrow)
        self.PlayButton.setToolButtonStyle(QtCore.Qt.ToolButtonIconOnly)
        self.PlayButton.toggled.connect(ASDMainWindow.PlayMovie)
        # Progress Bar showing the rendering of the moments
        self.ProgressBar = QtWidgets.QProgressBar(ASDMainWindow)
        self.ProgressBar.setValue(0)
        self.ProgressBar.setAlignment(QtCore.Qt.AlignCenter)
        self.ProgressBar.setTextVisible(False)

        self.ProgressLabel=QtWidgets.QLabel()
        self.ProgressLabel.setText('   {:}%'.format(int(self.ProgressBar.value())))
        # Creating the toolbar that displays several buttons in the windows
        # to visualize the data
        ASDVizToolbar = ASDMainWindow.addToolBar('ASD toolbar')
        ASDVizToolbar.addAction(self.exit)
        ASDVizToolbar.addActions(self.mode_group.actions())

        ASDVizToolbar.addAction(self.snap)
        ASDVizToolbar.addWidget(self.PlayButton)
        ASDVizToolbar.addWidget(self.ProgressBar)
        ASDVizToolbar.addWidget(self.ProgressLabel)

        ########################################################################
        # The actual menu bar
        ########################################################################
        # Creates the menu bar object that is added to the main window
        menubar= ASDMainWindow.menuBar()
        menubar.setNativeMenuBar(False)

        ########################################################################
        # File Menu
        ########################################################################
        # Adds an entry to the menu bar the file option
        fileMenu = menubar.addMenu('&File')
        # Setting up the coordinate file name from a File Dialog
        self.Select_coord = QtWidgets.QAction("Select coordinates file",fileMenu)
        self.Select_coord.triggered.connect(ASDVizMenu.ASD_data.getCoordFile)
        # Setting up the moment file name from a File Dialog
        self.Select_mom = QtWidgets.QAction("Select moments file",fileMenu)
        self.Select_mom.triggered.connect(ASDVizMenu.ASD_data.getMomentFile)
        # Setting up the restart file name from a File Dialog
        self.Select_restart = QtWidgets.QAction("Select restart file",fileMenu)
        self.Select_restart.triggered.connect(ASDVizMenu.ASD_data.getRestartFile)
        # Setting up the struct file name from a File Dialog
        self.Select_struct = QtWidgets.QAction("Select struct file",fileMenu)
        self.Select_struct.triggered.connect(ASDVizMenu.ASD_data.getStructFile)
        # Setting up the kmc_info file name from a File Dialog
        self.Select_KMC = QtWidgets.QAction("Select KMC file",fileMenu)
        self.Select_KMC.triggered.connect(ASDVizMenu.ASD_data.getKMCFile)
        fileMenu.addAction(self.Select_coord)
        fileMenu.addAction(self.Select_mom)
        fileMenu.addAction(self.Select_restart)
        fileMenu.addSeparator()

        fileMenu.addAction(self.Select_struct)
        fileMenu.addSeparator()
        fileMenu.addAction(self.Select_KMC)

        fileMenu.addSeparator()
        fileMenu.addAction(self.exit)
        ########################################################################
        # Snapshot Menu
        ########################################################################
        # Adds the Snapshot entry to the the menubar
        snapMenu=menubar.addMenu('&Snapshots')

        ASDVizMenu.png_snap = QtWidgets.QAction("Saving '.png'",snapMenu,checkable=True)
        ASDVizMenu.png_snap.setChecked(True)

        ASDVizMenu.pov_snap = QtWidgets.QAction("Saving '.pov'",snapMenu,checkable=True)
        ASDVizMenu.pov_snap.setChecked(True)

        snapMenu.addAction(ASDVizMenu.png_snap)
        snapMenu.addAction(ASDVizMenu.pov_snap)

        ########################################################################
        # Viz. Options Menu
        ########################################################################
        # Adds the Viz Options menu entry
        VizMenu=menubar.addMenu('&Viz. Options')
        # Adding a submenu for the Magnetization Density option
        Dens_menu=QtWidgets.QMenu("Magnetization Density",ASDMainWindow)
        # Adding an action for the visibility of the density
        self.Density_act = QtWidgets.QAction("Display",Dens_menu,checkable=True)
        self.Density_act.setChecked(True)
        self.Density_act.toggled.connect(ASDVizMenu.ASDOptions.toggle_density)
        # Creating an action group for the colors, of that way only one option can be selected
        Dens_proj_group=QtWidgets.QActionGroup(ASDMainWindow)
        Dens_proj_group.setExclusive(True)
        # X Projection
        self.Density_proj_x = QtWidgets.QAction("X Proj.",Dens_proj_group,checkable=True)
        self.Density_proj_x.setChecked(False)
        self.Density_proj_x.toggled.connect(ASDVizMenu.ASDOptions.set_color_x)
        # Y Projection
        self.Density_proj_y = QtWidgets.QAction("Y Proj.",Dens_proj_group,checkable=True)
        self.Density_proj_y.setChecked(False)
        self.Density_proj_y.toggled.connect(ASDVizMenu.ASDOptions.set_color_y)
        # Z Projection
        self.Density_proj_z = QtWidgets.QAction("Z Proj.",Dens_proj_group,checkable=True)
        self.Density_proj_z.setChecked(True)
        self.Density_proj_z.toggled.connect(ASDVizMenu.ASDOptions.set_color_z)
        # Adding the actors to the action group
        Dens_proj_group.addAction(self.Density_proj_x)
        Dens_proj_group.addAction(self.Density_proj_y)
        Dens_proj_group.addAction(self.Density_proj_z)
        # Adding a submenu for the Magnetization Density option
        Spins_menu=QtWidgets.QMenu("Atomic Spins",ASDMainWindow)
        # Adding an action for the visibility of the density
        self.Spins_act = QtWidgets.QAction("Display",Spins_menu,checkable=True)
        self.Spins_act.setChecked(False)
        self.Spins_act.toggled.connect(ASDVizMenu.ASDOptions.toggle_spins)
        # Creating an action group for the colors, of that way only one option can be selected
        Spins_proj_group=QtWidgets.QActionGroup(ASDMainWindow)
        Spins_proj_group.setExclusive(True)
        # X Projection
        self.Spin_proj_x = QtWidgets.QAction("X Proj.",Spins_proj_group,checkable=True)
        self.Spin_proj_x.setChecked(False)
        self.Spin_proj_x.toggled.connect(ASDVizMenu.ASDOptions.set_spins_color_x)
        # Y Projection
        self.Spin_proj_y = QtWidgets.QAction("Y Proj.",Spins_proj_group,checkable=True)
        self.Spin_proj_y.setChecked(False)
        self.Spin_proj_y.toggled.connect(ASDVizMenu.ASDOptions.set_spins_color_y)
        # Z Projection
        self.Spin_proj_z = QtWidgets.QAction("Z Proj.",Spins_proj_group,checkable=True)
        self.Spin_proj_z.setChecked(True)
        self.Spin_proj_z.toggled.connect(ASDVizMenu.ASDOptions.set_spins_color_z)
        # Adding the actors to the action group
        Spins_proj_group.addAction(self.Spin_proj_x)
        Spins_proj_group.addAction(self.Spin_proj_y)
        Spins_proj_group.addAction(self.Spin_proj_z)
        # Creating checkable menu entries for other options such as Axes, colorbar, etc.
        self.axes_menu_check = QtWidgets.QAction("Display Axes",VizMenu,checkable=True)
        self.axes_menu_check.setChecked(True)
        self.axes_menu_check.toggled.connect(ASDVizMenu.ASDOptions.toggle_Axes)
        self.ScalarBar_menu_check = QtWidgets.QAction("Display Color bar",VizMenu,checkable=True)
        self.ScalarBar_menu_check.setChecked(True)
        self.ScalarBar_menu_check.toggled.connect(ASDVizMenu.ASDOptions.toggle_ScalarBar)
        # Adding the directions of the spins to the Menu
        self.dir_menu_check = QtWidgets.QAction("Display Moments directions",VizMenu,checkable=True)
        self.dir_menu_check.setChecked(False)
        self.dir_menu_check.toggled.connect(ASDVizMenu.ASDOptions.toggle_directions)
        # Adding the contour toggling to the Menu
        self.cont_menu_check=QtWidgets.QAction("Display Contours",VizMenu,checkable=True)
        self.cont_menu_check.setChecked(False)
        self.cont_menu_check.toggled.connect(ASDVizMenu.ASDOptions.toggle_contours)
        # Adding an action group to handle the used color maps
        # Adding a submenu for the color map options
        Colors_menu=QtWidgets.QMenu("Color maps",ASDMainWindow)
        Colors_group=QtWidgets.QActionGroup(ASDMainWindow)
        Colors_group.setExclusive(True)
        # Coolwarm color map
        self.color_map_menu_Coolwarm = QtWidgets.QAction("Coolwarm Color Map",Colors_group,checkable=True)
        self.color_map_menu_Coolwarm.setChecked(True)
        self.color_map_menu_Coolwarm.toggled.connect(ASDVizMenu.ASDOptions.set_Coolwarm_lut)
        # Black Body Color map
        self.color_map_menu_BlackBody = QtWidgets.QAction("Black Body Color Map",Colors_group,checkable=True)
        self.color_map_menu_BlackBody.setChecked(False)
        self.color_map_menu_BlackBody.toggled.connect(ASDVizMenu.ASDOptions.set_BlackBody_lut)
        # RdGy Color map
        self.color_map_menu_RdGy = QtWidgets.QAction("RdGy Color Map",Colors_group,checkable=True)
        self.color_map_menu_RdGy.setChecked(False)
        self.color_map_menu_RdGy.toggled.connect(ASDVizMenu.ASDOptions.set_RdGy_lut)
        # Spectral Color map
        self.color_map_menu_Spectral = QtWidgets.QAction("Spectral Color Map",Colors_group,checkable=True)
        self.color_map_menu_Spectral.setChecked(False)
        self.color_map_menu_Spectral.toggled.connect(ASDVizMenu.ASDOptions.set_Spectral_lut)
        # Adding the actors to the action group
        Colors_group.addAction(self.color_map_menu_Coolwarm)
        Colors_group.addAction(self.color_map_menu_BlackBody)
        Colors_group.addAction(self.color_map_menu_RdGy)
        Colors_group.addAction(self.color_map_menu_Spectral)
        Colors_menu.addAction(self.color_map_menu_Coolwarm)
        Colors_menu.addAction(self.color_map_menu_BlackBody)
        Colors_menu.addAction(self.color_map_menu_RdGy)
        Colors_menu.addAction(self.color_map_menu_Spectral)
        ########################################################################
        # Different types of glyphs to visualize the moments
        ########################################################################
        SpinsGlyphs=QtWidgets.QMenu("Color maps",ASDMainWindow)
        SpinsGlyphs=QtWidgets.QActionGroup(ASDMainWindow)
        SpinsGlyphs.setExclusive(True)
        # Set the spins glyphs to be arrows
        self.spins_to_arrows = QtWidgets.QAction("Arrows",SpinsGlyphs,checkable=True)
        self.spins_to_arrows.setChecked(True)
        self.spins_to_arrows.toggled.connect(ASDMainWindow.ChangeGlyphs)
        # Set the spins glyphs to be cubes
        self.spins_to_cubes = QtWidgets.QAction("Cubes",SpinsGlyphs,checkable=True)
        self.spins_to_cubes.setChecked(False)
        self.spins_to_cubes.toggled.connect(ASDMainWindow.ChangeGlyphs)
        # Set the spins glyphs to be spheres
        self.spins_to_spheres = QtWidgets.QAction("Spheres",SpinsGlyphs,checkable=True)
        self.spins_to_spheres.setChecked(False)
        self.spins_to_spheres.toggled.connect(ASDMainWindow.ChangeGlyphs)
        # Set the spins glyphs to be cones
        self.spins_to_cones = QtWidgets.QAction("Cones",SpinsGlyphs,checkable=True)
        self.spins_to_cones.setChecked(False)
        self.spins_to_cones.toggled.connect(ASDMainWindow.ChangeGlyphs)

        SpinsGlyphs.addAction(self.spins_to_arrows)
        SpinsGlyphs.addAction(self.spins_to_cubes)
        SpinsGlyphs.addAction(self.spins_to_spheres)
        SpinsGlyphs.addAction(self.spins_to_cones)
        # Adding the relevant entries for the visualization options menu
        VizMenu.addMenu(Dens_menu)
        VizMenu.addMenu(Spins_menu)
        VizMenu.addMenu(Colors_menu)
        VizMenu.addSeparator()
        VizMenu.addAction(self.axes_menu_check)
        VizMenu.addAction(self.ScalarBar_menu_check)
        VizMenu.addAction(self.dir_menu_check)
        VizMenu.addAction(self.cont_menu_check)
        # Adding the actions to the Density submenu
        Dens_menu.addAction(self.Density_act)
        Dens_menu.addSeparator()
        Dens_menu.addAction(self.Density_proj_x)
        Dens_menu.addAction(self.Density_proj_y)
        Dens_menu.addAction(self.Density_proj_z)
        # Adding the actions to the Spins submenu
        Spins_menu.addAction(self.Spins_act)
        Spins_menu.addSeparator()
        Spins_menu.addAction(self.Spin_proj_x)
        Spins_menu.addAction(self.Spin_proj_y)
        Spins_menu.addAction(self.Spin_proj_z)
        Spins_menu.addSeparator()
        Spins_menu.addAction(self.spins_to_arrows)
        Spins_menu.addAction(self.spins_to_cubes)
        Spins_menu.addAction(self.spins_to_spheres)
        Spins_menu.addAction(self.spins_to_cones)

        ########################################################################
        # Adding a camera option menu
        ########################################################################
        CameraMenu=menubar.addMenu('&Camera Options')
        self.reset_camera_menu=QtWidgets.QAction("Reset Camera",CameraMenu)
        self.reset_camera_menu.triggered.connect(ASDMainWindow.CameraHandling)

        self.projection_camera_menu=QtWidgets.QAction("Parallel Projection",CameraMenu,checkable=True)
        self.projection_camera_menu.setChecked(False)
        self.projection_camera_menu.triggered.connect(ASDMainWindow.CameraHandling)

        CameraMenu.addAction(self.reset_camera_menu)
        CameraMenu.addAction(self.projection_camera_menu)

    def close_app(self):
        QApplication.exit()
