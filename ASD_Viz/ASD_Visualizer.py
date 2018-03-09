#!/usr/bin/env vtkpython
################################################################################
# PROGRAM: ASD_visualizer
# @author Jonathan Chico (08/09/2017)
# @description
# This is the main routine to run the ASD Visualizer a vtk/python visualiztion GUI
# and tool that allows the user to visualize several quantities obtained from the
# UppASD spin dynamics software package.
################################################################################
import vtk
import sys
import time
import glob
import os.path
import ASDQTMenu
import ASDQTDocket
import ASDVTKReading
import ASDMomVTKActors
import ASDNeighVTKActors
from PyQt4 import QtCore, QtGui
from PyQt4.QtGui import QApplication
from vtk.qt4.QVTKRenderWindowInteractor import QVTKRenderWindowInteractor

################################################################################
# The class for the creation of the main window where the rendering will take place
################################################################################
class Ui_MainWindow(object):
    """ Class that defines the main window where all the widgets and rendering take place.
        Changes in this class should be minimal as the bulk of the options are setup
        in the class ASD_Viewer """
    def setupUi(self, MainWindow):
        MainWindow.setObjectName("MainWindow")
        MainWindow.resize(1480, 860)
        self.centralWidget = QtGui.QWidget(MainWindow)
        self.gridlayout = QtGui.QGridLayout(self.centralWidget)
        self.vtkWidget = QVTKRenderWindowInteractor(self.centralWidget)
        self.gridlayout.addWidget(self.vtkWidget, 0, 0, 1, 1)
        MainWindow.setCentralWidget(self.centralWidget)

        return

################################################################################
# Class where the majority of the actions are actually defined
################################################################################
class ASD_Viewer(QtGui.QMainWindow):
    """ Class where the visualization options, actors and widgets are actually
    defined. The majority  of the work performed in the ASD_visualizer is performed here.
    This is the class that should be modified when one wishes to add a new visualization mode,
    or when one wishes to change how a visualization is actually performed."""

    def __init__(self, parent = None):
        """ This is where the overall changes to the UI are made. In here one should add
        any button to the toolbar, any entry to the menu bar or large widget such as the
        dock widget.
        Several variables key for the visualization that are required for other
        functions are also defined here."""
        QtGui.QMainWindow.__init__(self,parent)
        self.timer_count=0
        self.number_of_screenshots=0
        # Set the name of the window
        self.setWindowTitle('ASD Visualizer')
        # Create the UI
        self.ui=Ui_MainWindow()
        self.ui.setupUi(self)

        self.widget=QVTKRenderWindowInteractor(self)
        # Create the renderer
        self.ren = vtk.vtkRenderer()
        # Add the renderer to the main window
        self.ui.vtkWidget.GetRenderWindow().AddRenderer(self.ren)
        self.renWin=self.ui.vtkWidget.GetRenderWindow()
        # Set the window interactor
        self.iren = self.ui.vtkWidget.GetRenderWindow().GetInteractor()
        self.ren.SetBackground(1.0,1.0,1.0)
        # This calls the class in charge of creating all the needed menus
        self.ASDMenu=ASDQTMenu.ASDVizMenu()
        self.ASDMenu.ASDMenuCreator(ASDMainWindow=self)
        # Take a snapshot of the current window
        self.ASDMenu.snap.triggered.connect(self.Snapshot)
        # Making sure tha that the renderer is updated when the density actors change
        self.ASDMenu.Density_act.toggled.connect(self.UpdateRenderer)
        self.ASDMenu.Density_proj_x.toggled.connect(self.UpdateRenderer)
        self.ASDMenu.Density_proj_y.toggled.connect(self.UpdateRenderer)
        self.ASDMenu.Density_proj_z.toggled.connect(self.UpdateRenderer)
        self.ASDMenu.cont_menu_check.toggled.connect(self.UpdateRenderer)
        # Making sure that the renderer is updated when the spins change
        self.ASDMenu.Spins_act.toggled.connect(self.UpdateRenderer)
        self.ASDMenu.Spin_proj_x.toggled.connect(self.UpdateRenderer)
        self.ASDMenu.Spin_proj_y.toggled.connect(self.UpdateRenderer)
        self.ASDMenu.Spin_proj_z.toggled.connect(self.UpdateRenderer)
        # Making sure that the renderer is updated when the widgets are changed
        self.ASDMenu.axes_menu_check.toggled.connect(self.UpdateRenderer)
        self.ASDMenu.ScalarBar_menu_check.toggled.connect(self.UpdateRenderer)
        # Making sure the renderer is updated when the color maps are changed
        self.ASDMenu.color_map_menu_RdGy.toggled.connect(self.UpdateRenderer)
        self.ASDMenu.color_map_menu_Coolwarm.toggled.connect(self.UpdateRenderer)
        self.ASDMenu.color_map_menu_Spectral.toggled.connect(self.UpdateRenderer)
        self.ASDMenu.color_map_menu_BlackBody.toggled.connect(self.UpdateRenderer)
        # Making sure that the renderer is updated when misc options are changed
        self.ASDMenu.dir_menu_check.toggled.connect(self.UpdateRenderer)
        ########################################################################
        # Properties and calls for the Docked Window
        ########################################################################
        # Create the actual docket
        self.ASDQT=ASDQTDocket.ASDQTDockWindow()
        self.ASDQT.createASDDock(ASDMainWindow=self)
        # Making sure that the renderer is updated when the widgets are changed
        self.ASDQT.axes_check.stateChanged.connect(self.UpdateRenderer)
        self.ASDQT.ScalarBar_check.stateChanged.connect(self.UpdateRenderer)
        # Making sure that the renderer is updated when the color maps are changed
        self.ASDQT.color_map_RdGy.toggled.connect(self.UpdateRenderer)
        self.ASDQT.color_map_Coolwarm.toggled.connect(self.UpdateRenderer)
        self.ASDQT.color_map_Spectral.toggled.connect(self.UpdateRenderer)
        self.ASDQT.color_map_BlackBody.toggled.connect(self.UpdateRenderer)
        # Making sure that the renderer is updated when the spins are changed
        self.ASDQT.SpinsGroup.toggled.connect(self.UpdateRenderer)
        self.ASDQT.SpinSizeSL.valueChanged.connect(self.UpdateRenderer)
        self.ASDQT.spins_x.toggled.connect(self.UpdateRenderer)
        self.ASDQT.spins_y.toggled.connect(self.UpdateRenderer)
        self.ASDQT.spins_z.toggled.connect(self.UpdateRenderer)
        # Making sure that the renderer is updated when the density is changed
        self.ASDQT.cont_check.stateChanged.connect(self.UpdateRenderer)
        self.ASDQT.DensityGroup.toggled.connect(self.UpdateRenderer)
        self.ASDQT.dens_x.toggled.connect(self.UpdateRenderer)
        self.ASDQT.dens_y.toggled.connect(self.UpdateRenderer)
        self.ASDQT.dens_z.toggled.connect(self.UpdateRenderer)
        # Making sure that the renderer is updated when the clipper is changed
        self.ASDQT.ClipperGroup.toggled.connect(self.UpdateRenderer)
        self.ASDQT.ClipperSL.valueChanged.connect(self.UpdateRenderer)
        # Making sure the renderer is updated when misc actors are changed
        self.ASDQT.dir_check.stateChanged.connect(self.UpdateRenderer)
        self.ASDQT.KMC_check.stateChanged.connect(self.UpdateRenderer)
        self.ASDQT.cluster_check.stateChanged.connect(self.UpdateRenderer)
        # Making sure that the renderer is updated when the Neighbours are changed
        self.ASDQT.NeighSL.valueChanged.connect(self.UpdateRenderer)
        self.ASDQT.atom_viz.stateChanged.connect(self.UpdateRenderer)
        self.ASDQT.neigh_viz.stateChanged.connect(self.UpdateRenderer)
        self.ASDQT.NeighSL.SetAtom.editingFinished.connect(self.UpdateRenderer)
        ########################################################################
        # End of docket definitions
        ########################################################################
        return;

    ############################################################################
    # Wrapper function that takes care of adding the necessary actors and the
    # options for the different types of visualizations
    ############################################################################
    def AddActors(self):
        self.ren.RemoveAllViewProps()
        self.ren.ResetCamera()
        ########################################################################
        # This takes care of setting up the options for the visualization of the
        # magnetic moments obtained from the restart file
        ########################################################################
        if self.ASDMenu.restart.isChecked():
            self.mode=1
            self.viz_type='M'
            self.ASDQT.MomentBox.setCheckable(True)
            self.ASDQT.NeighBox.setChecked(False)
            self.ASDMenu.moments.setChecked(False)
            self.ASDMenu.neighbours.setChecked(False)
            self.MomActors=ASDMomVTKActors.ASDMomActors()
            self.MomActors.AddASD_actors(ren=self.ren,renWin=self.renWin,\
            mode=self.mode,viz_type=self.viz_type,iren=self.iren)
            if self.MomActors.kmc_disp:
                self.ASDQT.KMC_check.setVisible(True)
            if self.MomActors.cluster_disp:
                self.ASDQT.cluster_check.setVisible(True)
            self.ASDQT.update_dock_info()
            print 'Visualization of magnetic moments mode chosen'
            print 'Viewing the restart file'

        ########################################################################
        # This takes care of setting up the options for the visualization of the
        # magnetic moments obtained form the moment file
        ########################################################################
        if self.ASDMenu.moments.isChecked():
            self.mode=2
            self.viz_type='M'
            self.ASDQT.MomentBox.setCheckable(True)
            self.ASDQT.NeighBox.setChecked(False)
            self.ASDMenu.restart.setChecked(False)
            self.ASDMenu.neighbours.setChecked(False)
            self.MomActors=ASDMomVTKActors.ASDMomActors()
            self.MomActors.AddASD_actors(ren=self.ren,renWin=self.renWin,\
            mode=self.mode,viz_type=self.viz_type,iren=self.iren)
            if self.MomActors.kmc_disp:
                self.ASDQT.KMC_check.setVisible(True)
            if self.MomActors.cluster_disp:
                self.ASDQT.cluster_check.setVisible(True)
            self.ASDQT.update_dock_info()
            print 'Visualization of magnetic moments mode chosen'
            print 'Viewing the moment file'

        ########################################################################
        # This takes care of setting up the options for the Neighbour visualization
        ########################################################################
        if self.ASDMenu.neighbours.isChecked():
            self.mode=1
            self.viz_type='N'
            self.ASDMenu.restart.setChecked(False)
            self.ASDMenu.moments.setChecked(False)
            self.ASDQT.MomentBox.setChecked(False)
            self.ASDQT.NeighBox.setCheckable(True)
            self.NeighActors=ASDNeighVTKActors.ASDNeighActors()
            self.NeighActors.AddASDNeigh_actors(ren=self.ren,renWin=self.renWin,\
            mode=self.mode,viz_type=self.viz_type,iren=self.iren)
            print 'Visualization of the neighbour map mode chosen'
            self.ASDQT.NeighSL.setMaximum(self.NeighActors.SLMax-1)
            self.ASDQT.NeighSL.Nlabel.setText('Number of neighbours={:}'.format(self.NeighActors.NumNeigh-1))
        return

    ############################################################################
    # Function that calls the taking of a Snapshot of the current rendering window
    ############################################################################
    def Snapshot(self):
        self.MomActors=ASDMomVTKActors.ASDMomActors()
        self.MomActors.Screenshot(renWin=self.renWin,number_of_screenshots=self.number_of_screenshots,\
        png_mode=self.ASDMenu.png_snap.isChecked(),pov_mode=self.ASDMenu.pov_snap.isChecked())
        self.number_of_screenshots=self.number_of_screenshots+1

        return;

    ############################################################################
    # Set the conenctor to the play button
    ############################################################################
    def PlayMovie(self):
        self.iren.AddObserver('TimerEvent',self.execute)
        timerId = self.iren.CreateRepeatingTimer(100)
        self.iren.SetStillUpdateRate(15.0)
        self.iren.SetDesiredUpdateRate(15.0)
        self.iren.Start()
        self.renWin.Render()

        return

    ############################################################################
    # Create the information needed to update the moments
    # @author Anders Berman
    # Modifications for ASD_visualizer by Jonathan Chico
    ############################################################################
    def execute(self,obj,event):

        self.MomActors=ASDMomVTKActors.ASDMomActors()
        ASD_data=ASDVTKReading.ASDReading()

        if self.timer_count<=ASD_data.number_time_steps-1:

            # Update the magnetic moments
            (ASD_data.rest_mom,ASD_data.colors_x,ASD_data.colors_y,ASD_data.colors_z,\
            ASD_data.number_time_steps)=\
            ASD_data.readVectorsData(ASD_data.MomentsFile,self.timer_count,\
            ASD_data.nrAtoms,ASD_data.mode,ASD_data.number_time_steps)

            # If the KMC files are present updated them
            if ASD_data.kmc_flag:
                (ASD_data.coord_KMC,ASD_data.nrKMCpar)=\
                ASD_data.readKMCData(ASD_data.KMCFile,self.timer_count,ASD_data.nrKMCpar)
                self.KMC_src.SetPoints(ASD_data.coord_KMC)

            # If the data is 2D do the appropiate prunning
            if self.MomActors.glob_flag_2D:
                (ASD_data.selected_points,ASD_data.selected_vectors,ASD_data.selected_colors_x,\
                ASD_data.selected_colors_y,ASD_data.selected_colors_z,ASD_data.reduced_nrAtoms)=\
                ASD_data.create2Dsystem(ASD_data.rest_mom,ASD_data.colors_x,ASD_data.colors_y,ASD_data.colors_z,\
                ASD_data.coord,ASD_data.nrAtoms)
            else:
                ASD_data.selected_points=ASD_data.coord
                ASD_data.selected_vectors=ASD_data.rest_mom
                ASD_data.selected_colors_x=ASD_data.colors_x
                ASD_data.selected_colors_y=ASD_data.colors_y
                ASD_data.selected_colors_z=ASD_data.colors_z
                ASD_data.reduced_nrAtoms=ASD_data.nrAtoms

                # Check for the type of data needed
            if self.ASDQT.dens_x.isChecked():
                self.MomActors.src.GetPointData().SetScalars(ASD_data.selected_colors_x)
                self.MomActors.src_spins.GetPointData().SetScalars(ASD_data.selected_colors_x)
            if self.ASDQT.dens_y.isChecked():
                self.MomActors.src.GetPointData().SetScalars(ASD_data.selected_colors_y)
                self.MomActors.src_spins.GetPointData().SetScalars(ASD_data.selected_colors_y)
            if self.ASDQT.dens_z.isChecked():
                self.MomActors.src.GetPointData().SetScalars(ASD_data.selected_colors_z)
                self.MomActors.src_spins.GetPointData().SetScalars(ASD_data.selected_colors_z)

            self.MomActors.src.GetPointData().SetVectors(ASD_data.selected_vectors)
            self.MomActors.src_spins.GetPointData().SetVectors(ASD_data.selected_vectors)
            self.iren.Start()
            self.ASDMenu.ProgressBar.setValue(self.timer_count*100/(ASD_data.number_time_steps-1))
            self.ASDMenu.ProgressLabel.setText('   {:}%'.format(int(self.ASDMenu.ProgressBar.value())))
            self.renWin.Render()
            self.MomActors.Screenshot(renWin=self.renWin,number_of_screenshots=self.number_of_screenshots,\
            png_mode=self.ASDMenu.png_snap.isChecked(),pov_mode=self.ASDMenu.pov_snap.isChecked())
            self.timer_count=self.timer_count+1
            self.number_of_screenshots=self.number_of_screenshots+1
        # This makes sure the timer stops at the end of the time steps
        elif self.timer_count>=ASD_data.number_time_steps-1:
            self.ASDMenu.PlayButton.setChecked(False)
            self.iren.DestroyTimer()

        return

    ############################################################################
    # Camera updater. This wrapper function will take care of the majority of
    # the camera options
    ############################################################################
    def CameraHandling(self):
        if self.sender() == self.ASDQT.reset_button or self.sender() == self.ASDMenu.reset_camera_menu:
            self.MomActors=ASDMomVTKActors.ASDMomActors()
            self.MomActors.reset_camera(ren=self.ren,renWin=self.renWin)
        if self.sender() == self.ASDQT.Camera_proj_x:
            self.MomActors=ASDMomVTKActors.ASDMomActors()
            self.MomActors.set_Camera_x(ren=self.ren,renWin=self.renWin)
        if self.sender() == self.ASDQT.Camera_proj_y:
            self.MomActors=ASDMomVTKActors.ASDMomActors()
            self.MomActors.set_Camera_y(ren=self.ren,renWin=self.renWin)
        if self.sender() == self.ASDQT.Camera_proj_z:
            self.MomActors=ASDMomVTKActors.ASDMomActors()
            self.MomActors.set_Camera_z(ren=self.ren,renWin=self.renWin)
        if self.sender() == self.ASDQT.accept_button:
            self.ASDQT.getCameraData()
            self.ASDQT.update_Camera(ren=self.ren,renWin=self.renWin)
        if self.sender() == self.ASDQT.Projection_box:
            if self.ASDQT.Projection_box.isChecked():
                self.ren.GetActiveCamera().ParallelProjectionOn()
                self.ASDQT.toggle_projections(ren=self.ren)
                self.ASDMenu.projection_camera_menu.setChecked(True)
                self.renWin.Render()
            elif not self.ASDQT.Projection_box.isChecked():
                self.ren.GetActiveCamera().ParallelProjectionOff()
                self.ASDMenu.projection_camera_menu.setChecked(False)
                self.renWin.Render()
        if self.sender() == self.ASDMenu.projection_camera_menu:
            if self.ASDMenu.projection_camera_menu.isChecked():
                self.ren.GetActiveCamera().ParallelProjectionOn()
                self.ASDQT.toggle_projections(ren=self.ren)
                self.ASDQT.Projection_box.setChecked(True)
                self.renWin.Render()
            elif not self.ASDMenu.projection_camera_menu.isChecked():
                self.ren.GetActiveCamera().ParallelProjectionOff()
                self.ASDQT.Projection_box.setChecked(False)
                self.renWin.Render()
        if self.sender() == self.ASDQT.proj_sl_label:
            if self.ASDQT.proj_sl_label.editingFinished:
                self.ASDQT.ChangeParallelProj(ren=self.ren,renWin=self.renWin,line=True,slider=False)
        if self.sender() == self.ASDQT.proj_sl:
            self.ASDQT.ChangeParallelProj(ren=self.ren,renWin=self.renWin,line=False,slider=True)
        return

    ############################################################################
    # Function to change the type of glyphs that display the individual magnetic
    # moments
    ############################################################################
    def ChangeGlyphs(self):
        if self.sender() == self.ASDMenu.spins_to_cubes:
            if self.ASDMenu.spins_to_cubes.isChecked():
                self.MomActors=ASDMomVTKActors.ASDMomActors()
                self.MomActors.ChangeSpinGlyph(renWin=self.renWin,keyword='Cubes')
        if self.sender() == self.ASDMenu.spins_to_spheres:
            if self.ASDMenu.spins_to_spheres.isChecked():
                self.MomActors=ASDMomVTKActors.ASDMomActors()
                self.MomActors.ChangeSpinGlyph(renWin=self.renWin,keyword='Spheres')
        if self.sender() == self.ASDMenu.spins_to_arrows:
            if self.ASDMenu.spins_to_arrows.isChecked():
                self.MomActors=ASDMomVTKActors.ASDMomActors()
                self.MomActors.ChangeSpinGlyph(renWin=self.renWin,keyword='Arrows')
        if self.sender() == self.ASDMenu.spins_to_cones:
            if self.ASDMenu.spins_to_cones.isChecked():
                self.MomActors=ASDMomVTKActors.ASDMomActors()
                self.MomActors.ChangeSpinGlyph(renWin=self.renWin,keyword='Cones')


    ############################################################################
    # Function that takes care of updating the renderer and the Dock and Menu
    # everytime that an action is changed
    ############################################################################
    def UpdateRenderer(self):
        if self.sender() == self.ASDQT.axes_check:
            if self.ASDQT.axes_check.isChecked():
                self.ASDMenu.axes_menu_check.setChecked(True)
            else:
                self.ASDMenu.axes_menu_check.setChecked(False)
        if self.sender() == self.ASDMenu.axes_menu_check:
            if self.ASDMenu.axes_menu_check.isChecked():
                self.ASDQT.axes_check.setChecked(True)
            else:
                self.ASDQT.axes_check.setChecked(False)
        if self.sender() == self.ASDQT.ScalarBar_check:
            if self.ASDQT.ScalarBar_check.isChecked():
                self.ASDMenu.ScalarBar_menu_check.setChecked(True)
            else:
                self.ASDMenu.ScalarBar_menu_check.setChecked(False)
        if self.sender() == self.ASDMenu.ScalarBar_menu_check:
            if self.ASDMenu.ScalarBar_menu_check.isChecked():
                self.ASDQT.ScalarBar_check.setChecked(True)
            else:
                self.ASDQT.ScalarBar_check.setChecked(False)
        if self.sender() == self.ASDQT.color_map_RdGy:
            if self.ASDQT.color_map_RdGy.isChecked():
                self.ASDMenu.color_map_menu_RdGy.setChecked(True)
            else:
                self.ASDMenu.color_map_menu_RdGy.setChecked(False)
        if self.sender() == self.ASDMenu.color_map_menu_RdGy:
            if self.ASDMenu.color_map_menu_RdGy.isChecked():
                self.ASDQT.color_map_RdGy.setChecked(True)
            else:
                self.ASDQT.color_map_RdGy.setChecked(False)
        if self.sender() == self.ASDQT.color_map_Coolwarm:
            if self.ASDQT.color_map_Coolwarm.isChecked():
                self.ASDMenu.color_map_menu_Coolwarm.setChecked(True)
            else:
                self.ASDMenu.color_map_menu_Coolwarm.setChecked(False)
        if self.sender() == self.ASDMenu.color_map_menu_Coolwarm:
            if self.ASDMenu.color_map_menu_Coolwarm.isChecked():
                self.ASDQT.color_map_Coolwarm.setChecked(True)
            else:
                self.ASDQT.color_map_Coolwarm.setChecked(False)
        if self.sender() == self.ASDQT.color_map_BlackBody:
            if self.ASDQT.color_map_BlackBody.isChecked():
                self.ASDMenu.color_map_menu_BlackBody.setChecked(True)
            else:
                self.ASDMenu.color_map_menu_BlackBody.setChecked(False)
        if self.sender() == self.ASDMenu.color_map_menu_BlackBody:
            if self.ASDMenu.color_map_menu_BlackBody.isChecked():
                self.ASDQT.color_map_BlackBody.setChecked(True)
            else:
                self.ASDQT.color_map_BlackBody.setChecked(False)
        if self.sender() == self.ASDQT.color_map_Spectral:
            if self.ASDQT.color_map_Spectral.isChecked():
                self.ASDMenu.color_map_menu_Spectral.setChecked(True)
            else:
                self.ASDMenu.color_map_menu_Spectral.setChecked(False)
        if self.sender() == self.ASDMenu.color_map_menu_Spectral:
            if self.ASDMenu.color_map_menu_Spectral.isChecked():
                self.ASDQT.color_map_Spectral.setChecked(True)
            else:
                self.ASDQT.color_map_Spectral.setChecked(False)
        if self.sender() == self.ASDQT.SpinsGroup:
            if self.ASDQT.SpinsGroup.isChecked():
                self.ASDMenu.Spins_act.setChecked(True)
            else:
                self.ASDMenu.Spins_act.setChecked(False)
        if self.sender() == self.ASDMenu.Spins_act:
            if self.ASDMenu.Spins_act.isChecked():
                self.ASDQT.SpinsGroup.setChecked(True)
            else:
                self.ASDQT.SpinsGroup.setChecked(False)
        if self.sender() == self.ASDQT.spins_x:
            if self.ASDQT.spins_x.isChecked():
                self.ASDMenu.Spin_proj_x.setChecked(True)
            else:
                self.ASDMenu.Spin_proj_x.setChecked(False)
        if self.sender() == self.ASDMenu.Spin_proj_x:
            if self.ASDMenu.Spin_proj_x.isChecked():
                self.ASDQT.spins_x.setChecked(True)
            else:
                self.ASDQT.spins_x.setChecked(False)
        if self.sender() == self.ASDQT.spins_y:
            if self.ASDQT.spins_y.isChecked():
                self.ASDMenu.Spin_proj_y.setChecked(True)
            else:
                self.ASDMenu.Spin_proj_y.setChecked(False)
        if self.sender() == self.ASDMenu.Spin_proj_y:
            if self.ASDMenu.Spin_proj_y.isChecked():
                self.ASDQT.spins_y.setChecked(True)
            else:
                self.ASDQT.spins_y.setChecked(False)
        if self.sender() == self.ASDQT.spins_z:
            if self.ASDQT.spins_z.isChecked():
                self.ASDMenu.Spin_proj_z.setChecked(True)
            else:
                self.ASDMenu.Spin_proj_z.setChecked(False)
        if self.sender() == self.ASDMenu.Spin_proj_z:
            if self.ASDMenu.Spin_proj_z.isChecked():
                self.ASDQT.spins_z.setChecked(True)
            else:
                self.ASDQT.spins_z.setChecked(False)
        if self.sender() == self.ASDQT.DensityGroup:
            if self.ASDQT.DensityGroup.isChecked():
                self.ASDMenu.Density_act.setChecked(True)
            else:
                self.ASDMenu.Density_act.setChecked(False)
        if self.sender() == self.ASDMenu.Density_act:
            if self.ASDMenu.Density_act.isChecked():
                self.ASDQT.DensityGroup.setChecked(True)
            else:
                self.ASDQT.DensityGroup.setChecked(False)
        if self.sender() == self.ASDQT.dens_x:
            if self.ASDQT.dens_x.isChecked():
                self.ASDMenu.Density_proj_x.setChecked(True)
            else:
                self.ASDMenu.Density_proj_x.setChecked(False)
        if self.sender() == self.ASDMenu.Density_proj_x:
            if self.ASDMenu.Density_proj_x.isChecked():
                self.ASDQT.dens_x.setChecked(True)
            else:
                self.ASDQT.dens_x.setChecked(False)
        if self.sender() == self.ASDQT.dens_y:
            if self.ASDQT.dens_y.isChecked():
                self.ASDMenu.Density_proj_y.setChecked(True)
            else:
                self.ASDMenu.Density_proj_y.setChecked(False)
        if self.sender() == self.ASDMenu.Density_proj_y:
            if self.ASDMenu.Density_proj_y.isChecked():
                self.ASDQT.dens_y.setChecked(True)
            else:
                self.ASDQT.dens_y.setChecked(False)
        if self.sender() == self.ASDQT.dens_z:
            if self.ASDQT.dens_z.isChecked():
                self.ASDMenu.Density_proj_z.setChecked(True)
            else:
                self.ASDMenu.Density_proj_z.setChecked(False)
        if self.sender() == self.ASDMenu.Density_proj_z:
            if self.ASDMenu.Density_proj_z.isChecked():
                self.ASDQT.dens_z.setChecked(True)
            else:
                self.ASDQT.dens_z.setChecked(False)
        if self.sender() == self.ASDQT.dir_check:
            if self.ASDQT.dir_check.isChecked():
                self.ASDMenu.dir_menu_check.setChecked(True)
            else:
                self.ASDMenu.dir_menu_check.setChecked(False)
        if self.sender() == self.ASDMenu.dir_menu_check:
            if self.ASDMenu.dir_menu_check.isChecked():
                self.ASDQT.dir_check.setChecked(True)
            else:
                self.ASDQT.dir_check.setChecked(False)
        if self.sender() == self.ASDQT.cont_check:
            if self.ASDQT.cont_check.isChecked():
                self.ASDMenu.cont_menu_check.setChecked(True)
            else:
                self.ASDMenu.cont_menu_check.setChecked(False)
        if self.sender() == self.ASDMenu.cont_menu_check:
            if self.ASDMenu.cont_menu_check.isChecked():
                self.ASDQT.cont_check.setChecked(True)
            else:
                self.ASDQT.cont_check.setChecked(False)
        self.renWin.Render()

        return

################################################################################
# Main program starts
################################################################################
if __name__ == "__main__":

    app = QApplication(sys.argv)
    window = ASD_Viewer()
    window.show()
    window.iren.Initialize() # Need this line to actually show the render inside Qt
    sys.exit(app.exec_())
