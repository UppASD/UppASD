"""@package ASDUIDriver
Main driver for the different visualization types in UppASD. It contains most of the
wrapper functions and general functions to setup the UI.

This class acts as a wrapper to call the necessary functions to display the necessary
UI information. It creates the main widgets for the GUI and attach them to the different
pages depending on the backend being used.

Author
----------
Jonathan Chico
"""
# pylint: disable=invalid-name, no-name-in-module, no-member

import glob
import os.path as path
from enum import Enum

from matplotlib.backends.backend_qtagg import FigureCanvas
from matplotlib.figure import Figure
from PyQt6.QtCore import QSignalBlocker
from PyQt6.QtWidgets import QCheckBox, QFileDialog, QMainWindow
from vtk import vtkInteractorStyleTrackballCamera, vtkOpenGLRenderer

# from matplotlib.backends.backend_qt5agg import FigureCanvas
from vtk.qt.QVTKRenderWindowInteractor import QVTKRenderWindowInteractor

import ASD_GUI.ASD_Interactive.interactiveASD as IntASD
import ASD_GUI.Input_Creator.ASDInputGen as ASDInputGen
import ASD_GUI.UI.ASDInteractiveTab as ASDInteractiveTab
from ASD_GUI.PLOT import ASDPlots2D, ASDPlotsReading
from ASD_GUI.UI import ASDInputWindows, ASDUIInitHelper, ASDUIActorHelper, ASDUIPlottingHelper
from ASD_GUI.UI.ASDMenuToolbar import (
    UpdateUI,
)
from ASD_GUI.VTK_Viz import (
    ASDVTKGenActors,
    ASDVTKMomActors,
    ASDVTKReading,
    ASDVTKVizOptions,
    ASDVTKColor,
    ASDVTKTexture,
    ASDVTKCamera,
)

try:
    import uppasd.simulator as ASDsimulator
except ImportError:
    ASDsimulator = None
    print("Warning: uppasd module not found, interactive functions disabled")


class Backend(Enum):
    """
    Enum representing different backends for the ASD GUI.

    Attributes:
        UppASD_VTK (int): Backend using VTK.
        UppASD_MAT (int): Backend using MATLAB.
        UppASD_INP (int): Backend using INP.
        UppASD_INT (int): Backend using INT.
    """

    UppASD_VTK = 1
    UppASD_MAT = 2
    UppASD_INP = 3
    UppASD_INT = 4


##########################################################################
# @brief Class that defines the main window where all the widgets and rendering take place.
# @details It controls the actions which take place in the GUI. It defined the main window
# that allows for the following features:
# - Visual input file creation.
# - Matplotlib plotting of several key \c UppASD outputs.
# - VTK rendering of 3D \c UppASD data.
# @author Jonathan Chico
##########################################################################
class UppASDVizMainWindow(QMainWindow):
    ##########################################################################
    # @brief Class constructor for the main window
    # @details Class constructor for the main waindow. It initializes the inclusion
    # of several auxiliary classes that are used to setup the GUI functionality.
    # @author Jonathan Chico
    ##########################################################################
    def __init__(self):
        super(UppASDVizMainWindow, self).__init__()
        # -----------------------------------------------------------------------
        # Define the array containing the file names necessary for visualization
        # -----------------------------------------------------------------------
        self.file_names = [None] * 6
        self.plotfile_names = [None] * 7
        self.current_time = 0
        self.number_of_screenshots = 0
        self.VTKWidgetPresent = False
        self.IntWidgetPresent = False
        self.IntLaunched = False
        self.can_plot_ams = False
        self.can_plot_sqw = False
        self.hdrifile = []
        self.hdrifile_gotten = False
        self.bwBackground = False
        self.bwSinglecolor = False
        # -----------------------------------------------------------------------
        # Plotting global variables
        # -----------------------------------------------------------------------
        self.plot2D_cmap_indx = 0
        self.SQW_proj_indx = 0
        self.MagDirIndx = [3]  # range(4)
        self.EneIndx = [0]  # range(11)
        # -----------------------------------------------------------------------
        # Call the classes defining several needed functions for the reading of
        # data and the visualization
        # -----------------------------------------------------------------------
        self.ASDdata = ASDVTKReading.ASDReading()
        self.ASDsim = None
        try:
            if ASDsimulator is not None:
                self.ASDsim = ASDsimulator.Simulator()
            else:
                raise ImportError("ASDsimulator is None")
        except (ImportError, AttributeError):
            print(
                "ASDsimulator module not found or is None. Interactive functions disabled"
            )
        # self.ASDsim = None
        # try:
        #     self.ASDsim = ASDsimulator.Simulator()
        # except NameError:
        #     print("Warning: uppasd module not found, interactive functions disabled")

        # Early initialization of MomActors, will be overwritten
        self.MomActors = ASDVTKMomActors.ASDMomActors()
        self.ASDVizOpt = ASDVTKVizOptions.ASDVizOptions()
        self.ASDTexture = ASDVTKTexture.ASDTexture()
        self.ASDGenActors = ASDVTKGenActors.ASDGenActors()
        self.ASDPlotData = ASDPlotsReading.ReadPlotData()
        self.ASDPlots2D = ASDPlots2D.Abstract2DPlot()
        self.ASDCorrelationPlots = ASDPlots2D.Correlation_Plots()
        self.ASDInputGen = ASDInputGen.ASDInputGen()
        self.RestartWindow = ASDInputWindows.RestartWindow()
        self.PosfileWindow = ASDInputWindows.PosfileWindow()
        self.MomfileWindow = ASDInputWindows.MomfileWindow()
        self.InitPhaseWindow = ASDInputWindows.InitPhaseWindow()
        self.JfileWindow = ASDInputWindows.JfileWindow()
        self.DMfileWindow = ASDInputWindows.DMfileWindow()
        self.KfileWindow = ASDInputWindows.KfileWindow()
        self.InteractiveDockWidget = ASDInteractiveTab.InteractiveDock(self)
        self.ASDColor = ASDVTKColor.ASDVTKColor()
        # -----------------------------------------------------------------------
        # Set better font size
        # -----------------------------------------------------------------------
        self.setStyleSheet("QWidget{font-size:10pt}")
        # -----------------------------------------------------------------------
        # Adding the vtk Interactors
        # -----------------------------------------------------------------------
        self.vtkWidget = QVTKRenderWindowInteractor()
        # -----------------------------------------------------------------------
        # Create the renderer
        # -----------------------------------------------------------------------
        self.ren = vtkOpenGLRenderer()
        # -----------------------------------------------------------------------
        # Add the renderer to the main window
        # -----------------------------------------------------------------------
        self.vtkWidget.GetRenderWindow().AddRenderer(self.ren)
        self.renWin = self.vtkWidget.GetRenderWindow()
        # -----------------------------------------------------------------------
        # Set the window interactor for the VTK visualization
        # -----------------------------------------------------------------------
        self.iren = self.vtkWidget.GetRenderWindow().GetInteractor()
        self.iren.SetInteractorStyle(vtkInteractorStyleTrackballCamera())
        self.ren.SetBackground(1.0, 1.0, 1.0)
        # -----------------------------------------------------------------------
        # Set up the interactive simulation window and renderer
        # -----------------------------------------------------------------------
        self.IntVtkWidget = QVTKRenderWindowInteractor()
        self.Intren = vtkOpenGLRenderer()
        self.IntVtkWidget.GetRenderWindow().AddRenderer(self.Intren)
        self.IntrenWin = self.IntVtkWidget.GetRenderWindow()
        self.Intiren = self.IntVtkWidget.GetRenderWindow().GetInteractor()
        self.Intiren.SetInteractorStyle(vtkInteractorStyleTrackballCamera())
        self.Intren.SetBackground(1.0, 1.0, 1.0)
        # -----------------------------------------------------------------------
        # Initialize and setup the camera manager
        # -----------------------------------------------------------------------
        self.ASDCamera = ASDVTKCamera.CameraManager(self.ren.GetActiveCamera())
        # -----------------------------------------------------------------------
        # Initialize and setup the necessary structures for the UI
        # -----------------------------------------------------------------------
        ASDUIInitHelper.SetupUI(self)
        ASDUIInitHelper.InitUI(self)
        # -----------------------------------------------------------------------
        # Set the Plotting Canvas for the matplotlib plots
        # -----------------------------------------------------------------------
        ASDUIPlottingHelper.InitPlotUI(self)
        self.Plotting_Figure = Figure(figsize=(9, 7))
        self.Plotting_Figure.set_tight_layout(True)
        self.Plotting_canvas = FigureCanvas(self.Plotting_Figure)
        self.Plotting_Figure.add_subplot(111)
        self.Plotting_Figure.patch.set_alpha(0.0)

        self.Plotting_ax = self.Plotting_Figure.axes[0]
        self.Plot2DLayout.addWidget(self.Plotting_canvas)
        self.Plotting_Figure3D = Figure(figsize=(9, 7))
        self.Plotting_Figure3D.set_tight_layout(True)
        self.Plotting_canvas3D = FigureCanvas(self.Plotting_Figure3D)
        self.Plotting_Figure3D.add_subplot(111, projection="3d")
        self.Plotting_Figure3D.patch.set_alpha(0.0)
        self.Plotting_ax3D = self.Plotting_Figure3D.axes[0]
        self.Plot3DLayout.addWidget(self.Plotting_canvas3D)

        self.backend = Backend.UppASD_VTK

        # Interactive object
        self.InteractiveVtk = IntASD.InteractiveASD(
            self.Intren, self.IntrenWin, self.Intiren, self.ASDsim
        )

        return

    ##########################################################################
    # @brief Wrapper for the writing of the input file
    # @details This function first create the dictionary of key words, this is then
    # populated with the parameters from the GUI. The dictionary is then cleaned
    # and written to file.
    # @author Jonathan Chico
    ##########################################################################

    def WriteInputFile(self):
        """
        Generates and writes the input file using ASDInputGen methods.
        """
        self.ASDInputGen.ASDSetDefaults()
        self.ASDInputGen.ASDInputGatherer(self)
        self.ASDInputGen.clean_var()
        self.ASDInputGen.write_inpsd()
        return

    ##########################################################################
    # @brief Wrapper to create the restartfile
    # @author Jonathan Chico
    ##########################################################################

    def create_restart(self):
        """
        Creates a restart file using the ASDInputGen instance.
        """
        self.RestartWindow.write_restartfile(self.ASDInputGen)
        return

    ##########################################################################
    # Choose which kind of backend one will use to display VTK based visualizations,
    # matplotlib based visualizations
    ##########################################################################
    def chooseBackend(self):
        """
        Selects and configures the backend based on the current mode selected in the ModeSelector.
        """
        if self.ModeSelector.currentIndex() == 0:
            print("VTK")
            self.backend = Backend.UppASD_VTK
            if not self.VTKWidgetPresent:
                self.VTKWidget_Layout.addWidget(self.vtkWidget)
                self.VTKWidgetPresent = True
        elif self.ModeSelector.currentIndex() == 1:
            print("Plot")
            self.backend = Backend.UppASD_MAT
            self.vtkWidget.setVisible(False)
        elif self.ModeSelector.currentIndex() == 2:
            print("Inp")
            self.backend = Backend.UppASD_INP
            self.vtkWidget.setVisible(False)
        elif self.ModeSelector.currentIndex() == 3:
            if self.ASDsim is None:
                print("Interactive features disabled")
                self.ModeSelector.setCurrentIndex(self.ModeSelector.oldIndex)
                return
            self.backend = Backend.UppASD_INT
            if not self.IntWidgetPresent:
                self.InteractiveWidget_Layout.addWidget(self.IntVtkWidget)
                # self.InteractiveWidget_Layout.addWidget(self.vtkWidget)
                self.IntWidgetPresent = True

            # Rest of the code
            if self.CheckForInteractorFiles() and not self.IntLaunched:
                ASDInteractiveTab.InitializeInteractor(self)
                self.IntLaunched = True

        self.ModeSelector.oldIndex = self.ModeSelector.currentIndex()
        self.ResetUI()
        return

    ##########################################################################
    # Reset the UI to change between the VTK based visualization and the matplotlib
    # based visualization
    ##########################################################################
    def ResetUI(self):
        """
        Resets the UI elements based on the selected backend.
        """
        if self.backend == Backend.UppASD_VTK:
            self.OptionDock.setVisible(True)
            self.OptionDock.setEnabled(True)
            self.MatPlotOptions.setVisible(False)
            self.MatPlotOptions.setEnabled(False)
            self.InpDockWidget.setVisible(False)
            self.InpDockWidget.setEnabled(False)
            self.vtkWidget.setVisible(True)
            self.IntVtkWidget.setVisible(False)
            self.InteractiveDockWidget.setVisible(False)
            self.InteractiveDockWidget.setEnabled(False)
        elif self.backend == Backend.UppASD_MAT:
            self.OptionDock.setVisible(False)
            self.OptionDock.setEnabled(False)
            self.MatPlotOptions.setVisible(True)
            self.MatPlotOptions.setEnabled(True)
            self.InpDockWidget.setVisible(False)
            self.InpDockWidget.setEnabled(False)
            self.vtkWidget.setVisible(False)
            self.IntVtkWidget.setVisible(False)
            self.InteractiveDockWidget.setVisible(False)
            self.InteractiveDockWidget.setEnabled(False)
        elif self.backend == Backend.UppASD_INP:
            self.OptionDock.setVisible(False)
            self.OptionDock.setEnabled(False)
            self.MatPlotOptions.setVisible(False)
            self.InpDockWidget.setVisible(True)
            self.InpDockWidget.setEnabled(True)
            self.vtkWidget.setVisible(False)
            self.IntVtkWidget.setVisible(False)
            self.InteractiveDockWidget.setVisible(False)
            self.InteractiveDockWidget.setEnabled(False)
        elif self.backend == Backend.UppASD_INT:
            self.OptionDock.setVisible(False)
            self.OptionDock.setEnabled(False)
            self.MatPlotOptions.setVisible(False)
            self.InpDockWidget.setVisible(False)
            self.InpDockWidget.setEnabled(False)
            self.IntVtkWidget.setVisible(True)
            self.vtkWidget.setVisible(False)
            self.InteractiveDockWidget.setVisible(True)
            self.InteractiveDockWidget.setEnabled(True)
            self.IntrenWin.Render()
        return

    ##########################################################################
    # Initialization of some of the UI properties for 2D plots
    ##########################################################################
    #def InitPlotUI(self):
    #    """
    #    Initializes the plot UI by setting file names and disabling certain UI elements.
    #    """
    #    self.plotfile_names[0] = self.ASDPlotData.yamlfile
    #    self.plotfile_names[1] = self.ASDPlotData.amsfile
    #    self.plotfile_names[2] = self.ASDPlotData.sqwfile
    #    self.plotfile_names[3] = self.ASDPlotData.averages
    #    self.plotfile_names[4] = self.ASDPlotData.trajectory
    #    self.plotfile_names[5] = self.ASDPlotData.totenergy
    #    self.plotfile_names[6] = self.ASDPlotData.qfile
    #    self.SqwProjBox.setEnabled(False)
    #    self.SqwColorMapSelect.setEnabled(False)
    #    self.AveOpts.setEnabled(False)
    #    self.EneOpts.setEnabled(False)
    #    self.AMSDisplayOpts.setVisible(False)
    #    return

    ##########################################################################
    # Finding the file name for the VTK plots
    ##########################################################################
    def getFile(self):
        """
        Prompts the user to select a file using a file dialog.
        """
        self.ASDdata.getFileName(window=self)
        return

    ##########################################################################
    # Finding the file name for the matplotlib plots
    ##########################################################################
    def getPlotFile(self):
        """
        Opens a file dialog to get the plot file name.
        """
        self.ASDPlotData.getFileName(window=self)
        return

    ##########################################################################
    # Finding the file name for the input file generation
    ##########################################################################
    def getInpFile(self):
        """
        Opens a file dialog to select an input file using ASDInputGen.
        """
        self.ASDInputGen.getFileName(window=self)
        return

    ##########################################################################
    ##########################################################################

    def getInitPhase(self):
        """
        Handles the initialization phase when the InitPhaseDoneButton is pressed.
        """
        if self.sender() == self.InitPhaseWindow.InitPhaseDoneButton:
            self.init_phase_data = self.InitPhaseWindow.init_phase_data
        return

    ##########################################################################
    # @brief Open auxiliary windows for the inputfile creation GUI
    ##########################################################################

    def OpenWindow(self):
        """Wrapper function to display auxiliary windows in the Main Window. This handles
        the creation of the:
            - Posfile creation window
            - Momfile creation window
            - Restartfile creation window
            - Initial phase determination window
        Args
        ----------
            self: Parent Main window where the GUI is defined

        Author
        ----------
        Jonathan Chico
        """
        if self.sender() == self.InpInitMag4CreateButton:
            self.check_for_restart()
        if self.sender() == self.InpPosButtonCreate:
            print("Henlo friend")
            self.PosfileWindow.posfile_gotten = False
            # if self.InpCheckRandAlloy.isChecked():
            #     self.PosfileWindow.InPosBox.setEnabled(False)
            #     self.PosfileWindow.InPosBox.setVisible(False)
            #     self.PosfileWindow.InPosBoxRand.setEnabled(True)
            #     self.PosfileWindow.InPosBoxRand.setVisible(True)
            self.PosfileWindow.CheckForFile(self)
            self.PosfileWindow.InPosBox.setEnabled(True)
            self.PosfileWindow.InPosBox.setVisible(True)
            self.PosfileWindow.InPosBoxRand.setEnabled(False)
            self.PosfileWindow.InPosBoxRand.setVisible(False)
            self.PosfileWindow.show()
        if self.sender() == self.InpMomButtonCreate:
            self.MomfileWindow.momfile_gotten = False
            self.MomfileWindow.CheckForFile(self)
            self.MomfileWindow.show()
        if self.sender() == self.InpJfileButtonCreate:
            self.JfileWindow.jfile_gotten = False
            self.JfileWindow.CheckForFile(self)
            self.JfileWindow.show()
        if self.sender() == self.InpDMButtonCreate:
            self.DMfileWindow.DMfile_gotten = False
            self.DMfileWindow.CheckForFile(self)
            self.DMfileWindow.show()
        if self.sender() == self.InpKfileButtonCreate:
            self.KfileWindow.Kfile_gotten = False
            self.KfileWindow.CheckForFile(self)
            self.KfileWindow.show()
        if self.sender() == self.InpSetPhases:
            if self.InpInitLLG.isChecked():
                self.InitPhaseWindow.IpNphaseBox.setEnabled(True)
                self.InitPhaseWindow.IpNphaseBox.setVisible(True)
                self.InitPhaseWindow.MCannealBox.setEnabled(False)
                self.InitPhaseWindow.MCannealBox.setVisible(False)
                self.InitPhaseWindow.IpVPOBox.setVisible(False)
                self.InitPhaseWindow.IpVPOBox.setEnabled(False)
                self.InitPhaseWindow.InitPhaseAddButton.setEnabled(True)
                self.InitPhaseWindow.InitPhaseDelButton.setEnabled(True)
            if self.InpInitMcMet.isChecked() or self.InpInitMcHeat.isChecked():
                self.InitPhaseWindow.IpNphaseBox.setEnabled(False)
                self.InitPhaseWindow.IpNphaseBox.setVisible(False)
                self.InitPhaseWindow.MCannealBox.setEnabled(True)
                self.InitPhaseWindow.MCannealBox.setVisible(True)
                self.InitPhaseWindow.IpVPOBox.setVisible(False)
                self.InitPhaseWindow.IpVPOBox.setEnabled(False)
                self.InitPhaseWindow.InitPhaseAddButton.setEnabled(True)
                self.InitPhaseWindow.InitPhaseDelButton.setEnabled(True)
            if self.InpInitVPO.isChecked():
                self.InitPhaseWindow.IpNphaseBox.setEnabled(False)
                self.InitPhaseWindow.IpNphaseBox.setVisible(False)
                self.InitPhaseWindow.MCannealBox.setEnabled(False)
                self.InitPhaseWindow.MCannealBox.setVisible(False)
                self.InitPhaseWindow.IpVPOBox.setVisible(True)
                self.InitPhaseWindow.IpVPOBox.setEnabled(True)
                self.InitPhaseWindow.InitPhaseAddButton.setEnabled(False)
                self.InitPhaseWindow.InitPhaseDelButton.setEnabled(False)

            self.InitPhaseWindow.show()
        return

    ##########################################################################
    ##########################################################################

    def update_names(self):
        """
        Updates the file name using ASDInputGen.
        """
        self.ASDInputGen.update_file_name(window=self)
        return

    ##########################################################################
    # @brief Function to determine if the restartfile can be created.
    # @details This function will test if the lattice vectors have been defined,
    # as well as the posfile and momfile, that is everything which is necessary to
    # generate a restartfile
    # @author Jonathan Chico
    ##########################################################################
    def check_for_restart(self):
        """
        Checks if all required inputs are provided and attempts to restart the process.
        """
        everything_okay = True
        self.ASDInputGen.ASDInputGatherer(self)
        if not len(self.InpLineEditC1_x.text()) > 0:
            everything_okay = False
        if not len(self.InpLineEditC1_x.text()) > 0:
            everything_okay = False
        if not len(self.InpLineEditC1_z.text()) > 0:
            everything_okay = False
        if not len(self.InpLineEditC2_x.text()) > 0:
            everything_okay = False
        if not len(self.InpLineEditC2_y.text()) > 0:
            everything_okay = False
        if not len(self.InpLineEditC2_z.text()) > 0:
            everything_okay = False
        if not len(self.InpLineEditC3_x.text()) > 0:
            everything_okay = False
        if not len(self.InpLineEditC3_y.text()) > 0:
            everything_okay = False
        if not len(self.InpLineEditC3_z.text()) > 0:
            everything_okay = False
        # -----------------------------------------------------------------------
        if everything_okay:
            if (
                self.MomfileWindow.momfile_gotten or self.ASDInputGen.momfile_gotten
            ) and (
                self.PosfileWindow.posfile_gotten or self.ASDInputGen.posfile_gotten
            ):
                self.RestartWindow.restart_pre_generation(self.ASDInputGen)
                self.RestartWindow.show()
            else:
                self.Create_Restart_Error_Window = ASDInputWindows.ErrorWindow()
                self.Create_Restart_Error_Window.FunMsg.setText(
                    "I'm sorry, Dave. I'm afraid I can't do that."
                )
                self.Create_Restart_Error_Window.ErrorMsg.setText(
                    "Error: A posfile and a momfile must be first be defined"
                )
                self.Create_Restart_Error_Window.show()
                print("Error: A posfile and a momfile must be first be defined")
        else:
            self.Create_Restart_Error_Window = ASDInputWindows.ErrorWindow()
            self.Create_Restart_Error_Window.FunMsg.setText(
                "I'm sorry, Dave. I'm afraid I can't do that."
            )
            self.Create_Restart_Error_Window.ErrorMsg.setText(
                "Error: The unit cell vectors need to be defined"
            )
            self.Create_Restart_Error_Window.show()
            print("Error: The unit cell vectors need to be defined")
        return

    ##########################################################################
    # Function to select the appropriate data to plot
    ##########################################################################
    def PlottingSelector(self):
        """
        Handles the selection and plotting of different data types based on the sender action.
        """
        ASDUIPlottingHelper.PlottingSelector(self)
        return

    ##########################################################################
    # Select the projection of the S(q,w)
    ##########################################################################
    def SQW_Proj_Select(self):
        """
        Handles the selection of SQW projection index based on the sender.
        """
        if self.sender() == self.Sqw_x:
            self.SQW_proj_indx = 0
        if self.sender() == self.Sqw_y:
            self.SQW_proj_indx = 1
        if self.sender() == self.Sqw_z:
            self.SQW_proj_indx = 2
        if self.sender() == self.Sqw_2:
            self.SQW_proj_indx = 3
        self.PlottingWrapper()
        return

    ##########################################################################
    # Select the colormap over which the S(q,w) will be plotted
    ##########################################################################
    def Sqw_ColorMapSelector(self):
        """
        Selects the colormap for 2D plotting based on the sender of the signal.
        """
        if self.sender() == self.SqwCoolwarm:
            self.plot2D_cmap_indx = 0
        if self.sender() == self.SqwSpectral:
            self.plot2D_cmap_indx = 1
        if self.sender() == self.SqwBlackbody:
            self.plot2D_cmap_indx = 2
        self.PlottingWrapper()
        return

    ##########################################################################
    # Plotting the directions of the magnetization
    ##########################################################################
    def PlotMagDirSelector(self):
        """
        Updates the MagDirIndx list based on the checked state of plot options and calls the plotting function.
        """
        self.MagDirIndx = []
        if self.Plot_M_x.isChecked():
            self.MagDirIndx.append(0)
        if self.Plot_M_y.isChecked():
            self.MagDirIndx.append(1)
        if self.Plot_M_z.isChecked():
            self.MagDirIndx.append(2)
        if self.Plot_M_tot.isChecked():
            self.MagDirIndx.append(3)
        self.PlottingWrapper()

    ##########################################################################
    # Changing the marker size of the lines
    ##########################################################################
    def PlotLineChanger(self, value):
        """
        Adjusts the linewidth of the 2D plot and updates the plot.
        """
        self.ASDPlots2D.linewidth = value / 2.0
        self.PlottingWrapper()

    ##########################################################################
    # Changing the marker size of the lines
    ##########################################################################
    def PlotMarkerChanger(self, value):
        """
        Adjusts the marker size for 2D plots and updates the plot.
        """
        self.ASDPlots2D.markersize = value / 2.0
        self.PlottingWrapper()

    ##########################################################################
    # Changing the marker size of the lines
    ##########################################################################
    def PlotXGridToggle(self):
        """
        Toggles the visibility of the X-axis grid in the 2D plots.
        """
        self.ASDPlots2D.xgrid = not self.ASDPlots2D.xgrid
        self.PlottingWrapper()

    ##########################################################################
    # Changing the marker size of the lines
    ##########################################################################
    def PlotYGridToggle(self):
        """
        Toggles the visibility of the Y-axis grid in the 2D plots.
        """
        self.ASDPlots2D.ygrid = not self.ASDPlots2D.ygrid
        self.PlottingWrapper()

    ##########################################################################
    #  Toggling SQW grid lines on/off
    ##########################################################################
    def PlotSQWGridToggle(self):
        """
        Toggles the grid state for ASD correlation plots and updates the plot.
        """
        self.ASDCorrelationPlots.grid = not self.ASDCorrelationPlots.grid
        self.PlottingWrapper()

    ##########################################################################
    #  Toggling SQW grid lines on/off
    ##########################################################################
    def PlotAMSGridToggle(self):
        """
        Toggles the AMS grid visibility and updates the plot.
        """
        self.ASDPlots2D.amsgrid = not self.ASDPlots2D.amsgrid
        self.PlottingWrapper()

    ##########################################################################
    # Changing the width of S(q,w) plots
    ##########################################################################
    def SqwWidthChanger(self, value):
        """
        Adjusts the width parameter for the ASD correlation plots and updates the UI.
        """
        self.ASDCorrelationPlots.sigma_w = self.ASDCorrelationPlots.w_min * value
        self.ABCorrWidthTX.setText(f"{self.ASDCorrelationPlots.w_min*value:.3f}")
        self.PlottingWrapper()

    ##########################################################################
    # Plotting the components of the energy
    ##########################################################################
    def PlotEneCompSelector(self):
        """
        Updates the energy index list based on selected checkboxes and triggers plotting.
        """
        self.EneIndx = []
        if self.EneTotCheck.isChecked():
            self.EneIndx.append(0)
        if self.EneExcCheck.isChecked():
            self.EneIndx.append(1)
        if self.EneAniCheck.isChecked():
            self.EneIndx.append(2)
        if self.EneDMCheck.isChecked():
            self.EneIndx.append(3)
        if self.EnePdCheck.isChecked():
            self.EneIndx.append(4)
        if self.EneBqDMCheck.isChecked():
            self.EneIndx.append(5)
        if self.EneBqCheck.isChecked():
            self.EneIndx.append(6)
        if self.EneDipCheck.isChecked():
            self.EneIndx.append(7)
        if self.EneExtCheck.isChecked():
            self.EneIndx.append(8)
        if self.EneLSFCheck.isChecked():
            self.EneIndx.append(9)
        if self.EneChirCheck.isChecked():
            self.EneIndx.append(10)
        self.PlottingWrapper()
        return

    ##########################################################################
    ##########################################################################

    def ToggleInitPhase(self):
        """
        Toggles the initialization phase of the UI.
        """
        UpdateUI(self)
        return

    def ToggleHessians(self):
        """
        Toggles the Hessians in the UI.
        """
        UpdateUI(self)
        return

    ##########################################################################
    # @brief Function to selective plot the ams branches
    # @details Function to selectively plot the ams branches. It functions by
    # finding which of the checkboxes identifying each branch is selected
    # after this it creates a new data set that contains only the necessary data
    # @author Jonathan Chico
    ##########################################################################
    def AMS_PrunePlot(self):
        """
        Prunes and plots AMS data based on the state of checkboxes.
        """
        self.ams_data_y = []
        self.ams_data_x = []
        self.ams_label = []
        for ii in range(len(self.AMSCheckboxes)):
            name = f"ams_branch_{ii}"
            if self.AMSCheckboxes[name].isChecked():
                self.ams_data_x.append(self.ASDPlotData.ams_data_x[ii])
                self.ams_data_y.append(self.ASDPlotData.ams_data_y[ii])
                self.ams_label.append(self.ASDPlotData.ams_label[ii])
        self.PlottingWrapper()
        return

    ##########################################################################
    # @brief Wrapper function that takes care of plotting the selected plot
    # @details Wrapper function that takes care of plotting the selected plot, it allows
    # the user to choose between the following different types of plots
    #   - Spin-Spin correlation functions
    #       - S(q,w)
    #       - AMS
    #   - Magnetization averages
    #   - Single spin trajectories
    # @author Jonathan Chico
    ##########################################################################

    def PlottingWrapper(self):
        """Wrapper function that takes care of plotting the selected plot, it allows
        the user to choose between the following different types of plots:
            * Spin-Spin correlation functions:
                - S(q,w)
                - AMS
            * Magnetization averages
            * Single spin trajectories

        Author
        ----------
        Jonathan Chico
        """

        ASDUIPlottingHelper.PlottingWrapper(self)
        return

    ##########################################################################
    # @brief Function to save the current figure to file
    # @author Jonathan Chico
    ##########################################################################
    def SaveFig(self):
        """
        Save the current figure to a file with specified DPI.
        """
        ASDUIPlottingHelper.SaveFig(self)

    ##########################################################################
    # @brief Wrapper function that takes care of adding the necessary actors and the
    # options for the different types of visualizations
    # @details Wrapper function that takes care of adding the necessary actors and the
    # options for the different types of visualizations. It controls the visualization of
    #   - Restartfiles
    #   - Momentsfiles
    #   - Energy
    #   - Exchange neighbours
    #   - DM neighbours
    # @author Jonathan Chico
    ##########################################################################
    def AddActors(self):
        """
        Adds actors to the ASD UI.

        Utilizes ASDUIActorHelper to add actors to the current UI context.
        """
        ASDUIActorHelper.AddActors(self)

    ##########################################################################
    # @brief Enable rgb-values for single color
    # @author Anders Bergman
    ##########################################################################
    def toggle_singlecolor(self, check):
        """
        Enable or disable RGB color sliders based on the check value.
        """
        if check:
            self.RGBRedColorSlider.setEnabled(True)
            self.RGBGreenColorSlider.setEnabled(True)
            self.RGBBlueColorSlider.setEnabled(True)
        else:
            self.RGBRedColorSlider.setEnabled(False)
            self.RGBGreenColorSlider.setEnabled(False)
            self.RGBBlueColorSlider.setEnabled(False)

        return

    ##########################################################################
    # @brief Toggle grayscale background on/off
    # @author Anders Bergman
    ##########################################################################
    def toggle_bwSinglecolor(self, check):
        """
        Toggles the color sliders between black & white and single color mode.
        """
        self.bwSinglecolor = check
        rgb = [
            self.RGBRedColorSlider.value(),
            self.RGBGreenColorSlider.value(),
            self.RGBBlueColorSlider.value(),
        ]
        bw = int(sum(rgb) / 3)

        if check:
            self.RGBRedColorSlider.setValue(bw)
            self.RGBGreenColorSlider.setValue(bw)
            self.RGBRedColorSlider.setValue(bw)

        return

    ##########################################################################
    # @brief Toggle depth of field focus
    # @author Anders Bergman
    ##########################################################################
    def toggle_focus(self, check):
        """
        Toggles the focus state in the visualization options.
        """
        self.ASDVizOpt.toggle_Focus(check=check, ren=self.ren, renWin=self.renWin)

    ##########################################################################
    # @brief Toggle focal disk
    # @author Anders Bergman
    ##########################################################################
    def FocalDisk_control(self, value):
        """
        Controls the focal disk setting in the ASD visualization.
        """
        self.ASDVizOpt.setFocalDisk(value=value, ren=self.ren, renWin=self.renWin)

    ##########################################################################
    # @brief Toggle depth of field focus
    # @author Anders Bergman
    ##########################################################################
    def toggle_autofocus(self, check):
        """
        Toggles the autofocus feature in the visualization options.
        """
        self.ASDVizOpt.toggle_autoFocus(check=check, renWin=self.renWin)

    ##########################################################################
    # @brief Toggle grayscale background on/off
    # @author Anders Bergman
    ##########################################################################
    def toggle_bwBackground(self, check):
        """
        Toggles the background color between black and white based on the check value.
        """
        self.bwBackground = check
        rgb = [
            self.RGBRedBackgroundSlider.value(),
            self.RGBGreenBackgroundSlider.value(),
            self.RGBBlueBackgroundSlider.value(),
        ]
        bw = int(sum(rgb) / 3)

        if check:
            self.RGBRedBackgroundSlider.setValue(bw)
            self.RGBGreenBackgroundSlider.setValue(bw)
            self.RGBRedBackgroundSlider.setValue(bw)

        return

    ##########################################################################
    # @brief Update rgb-values for single color coloring
    # @author Anders Bergman
    ##########################################################################
    def set_singlecolor(self, value):
        """
        Set the single color for the RGB sliders and update the visualization.
        """
        self.ASDColor.set_singlecolor(window=self, value=value)

        return

    ##########################################################################
    # @brief Update rgb-values for the background
    # @author Anders Bergman
    ##########################################################################
    def set_background(self, value):
        """
        Sets the background color based on the provided value.
        """
        self.ASDColor.set_background(value=value, window=self)

        return

    ##########################################################################
    # @brief Set the lookup table for the actors
    # @details Set the lookup table for the actors, it also allows for the change
    # of the scale type for the plotting, either linear or logarithmic scale.
    # @author Jonathan Chico
    ##########################################################################
    def set_lut_db(self, mapnum):
        """
        Sets the lookup table (LUT) for the visualization based on the provided colormap number.
        """
        self.ASDColor.set_lut_db(mapnum=mapnum, window=self)
        self.renWin.Render()
        return

    ##########################################################################
    # @brief Set the lookup table for the actors
    # @details Set the lookup table for the actors, it also allows for the change
    # of the scale type for the plotting, either linear or logarithmic scale.
    # @author Jonathan Chico
    ##########################################################################
    def set_lut(self):
        """
        Sets the lookup table (LUT) scale based on the sender's state.
        """
        self.ASDColor.set_lut(window=self)

        return

    ##########################################################################
    # @brief Set the projection of the vectors
    # @details Set the projection of the vectors and the magnetization continuum
    # allowing one to set independent projections of the continuum visualization
    # and the spins.
    # @author Jonathan Chico
    ##########################################################################
    def set_projection(self):
        """Set the projection of the vectors and the magnetization continuum
        allowing one to set independent projections of the continuum visualization
        and the spins.

        Author
        ----------
        Jonathan Chico
        """
        if self.sender() == self.DensX and self.DensX.isChecked():
            self.MomActors.set_projection(atype="density", axis=0)
        if self.sender() == self.DensY and self.DensY.isChecked():
            self.MomActors.set_projection(atype="density", axis=1)
        if self.sender() == self.DensZ and self.DensZ.isChecked():
            self.MomActors.set_projection(atype="density", axis=2)
        if self.sender() == self.SpinX and self.SpinX.isChecked():
            self.MomActors.set_projection(atype="spins", axis=0)
        if self.sender() == self.SpinY and self.SpinY.isChecked():
            self.MomActors.set_projection(atype="spins", axis=1)
        if self.sender() == self.SpinZ and self.SpinZ.isChecked():
            self.MomActors.set_projection(atype="spins", axis=2)
        self.renWin.Render()
        return

    ##########################################################################
    # Display the different energy contributions
    ##########################################################################
    def set_energy_proj(self):
        """
        Sets the energy projection based on the sender button's state.
        """
        if self.sender() == self.TotEneButton and self.TotEneButton.isChecked():
            self.EneActors.src.GetPointData().SetScalars(self.ASDdata.energies[0])
            self.EneActors.EneMapper.SetScalarRange(self.EneActors.src.GetScalarRange())
            # self.EneActors.EneDensMap.SetScalarRange(self.EneActors.src.GetScalarRange())
            self.ASDGenActors.clipperMapper.SetScalarRange(
                self.EneActors.src.GetScalarRange()
            )
        if self.sender() == self.ExcEneButton and self.ExcEneButton.isChecked():
            self.EneActors.src.GetPointData().SetScalars(self.ASDdata.energies[1])
            self.EneActors.EneMapper.SetScalarRange(self.EneActors.src.GetScalarRange())
            # self.EneActors.EneDensMap.SetScalarRange(self.EneActors.src.GetScalarRange())
            self.ASDGenActors.clipperMapper.SetScalarRange(
                self.EneActors.src.GetScalarRange()
            )
        if self.sender() == self.DMEneButton and self.DMEneButton.isChecked():
            self.EneActors.src.GetPointData().SetScalars(self.ASDdata.energies[2])
            self.EneActors.EneMapper.SetScalarRange(self.EneActors.src.GetScalarRange())
            # self.EneActors.EneDensMap.SetScalarRange(self.EneActors.src.GetScalarRange())
            self.ASDGenActors.clipperMapper.SetScalarRange(
                self.EneActors.src.GetScalarRange()
            )
        if self.sender() == self.AniEneButton and self.AniEneButton.isChecked():
            self.EneActors.src.GetPointData().SetScalars(self.ASDdata.energies[3])
            self.EneActors.EneMapper.SetScalarRange(self.EneActors.src.GetScalarRange())
            # self.EneActors.EneDensMap.SetScalarRange(self.EneActors.src.GetScalarRange())
            self.ASDGenActors.clipperMapper.SetScalarRange(
                self.EneActors.src.GetScalarRange()
            )
        if self.sender() == self.BqEneButton and self.BqEneButton.isChecked():
            self.EneActors.src.GetPointData().SetScalars(self.ASDdata.energies[4])
            self.EneActors.EneMapper.SetScalarRange(self.EneActors.src.GetScalarRange())
            # self.EneActors.EneDensMap.SetScalarRange(self.EneActors.src.GetScalarRange())
            self.ASDGenActors.clipperMapper.SetScalarRange(
                self.EneActors.src.GetScalarRange()
            )
        if self.sender() == self.BqDMEneButton and self.BqDMEneButton.isChecked():
            self.EneActors.src.GetPointData().SetScalars(self.ASDdata.energies[5])
            self.EneActors.EneMapper.SetScalarRange(self.EneActors.src.GetScalarRange())
            # self.EneActors.EneDensMap.SetScalarRange(self.EneActors.src.GetScalarRange())
            self.ASDGenActors.clipperMapper.SetScalarRange(
                self.EneActors.src.GetScalarRange()
            )
        if self.sender() == self.PdEneButton and self.PdEneButton.isChecked():
            self.EneActors.src.GetPointData().SetScalars(self.ASDdata.energies[6])
            self.EneActors.EneMapper.SetScalarRange(self.EneActors.src.GetScalarRange())
            # self.EneActors.EneDensMap.SetScalarRange(self.EneActors.src.GetScalarRange())
            self.ASDGenActors.clipperMapper.SetScalarRange(
                self.EneActors.src.GetScalarRange()
            )
        if self.sender() == self.BextEneButton and self.BextEneButton.isChecked():
            self.EneActors.src.GetPointData().SetScalars(self.ASDdata.energies[7])
            self.EneActors.EneMapper.SetScalarRange(self.EneActors.src.GetScalarRange())
            # self.EneActors.EneDensMap.SetScalarRange(self.EneActors.src.GetScalarRange())
            self.ASDGenActors.clipperMapper.SetScalarRange(
                self.EneActors.src.GetScalarRange()
            )
        if self.sender() == self.DipEneButton and self.DipEneButton.isChecked():
            self.EneActors.src.GetPointData().SetScalars(self.ASDdata.energies[8])
            self.EneActors.EneMapper.SetScalarRange(self.EneActors.src.GetScalarRange())
            # self.EneActors.EneDensMap.SetScalarRange(self.EneActors.src.GetScalarRange())
            self.ASDGenActors.clipperMapper.SetScalarRange(
                self.EneActors.src.GetScalarRange()
            )
        if self.sender() == self.ChirEneButton and self.ChirEneButton.isChecked():
            self.EneActors.src.GetPointData().SetScalars(self.ASDdata.energies[9])
            self.EneActors.EneMapper.SetScalarRange(self.EneActors.src.GetScalarRange())
            # self.EneActors.EneDensMap.SetScalarRange(self.EneActors.src.GetScalarRange())
            self.ASDGenActors.clipperMapper.SetScalarRange(
                self.EneActors.src.GetScalarRange()
            )
        self.EneActors.EneDensMap.Modified()
        self.EneActors.EneMapper.Update()
        self.renWin.Render()
        return

    ##########################################################################
    # Function to change the type of glyphs that display the individual magnetic
    # moments
    ##########################################################################
    def ChangeGlyphs(self):
        """Function to change the type of glyphs that display the individual magnetic
        moments

        Author
        ----------
        Jonathan Chico
        """
        if self.sender() == self.SpinBarButton:
            if self.SpinBarButton.isChecked():
                self.MomActors.ChangeSpinGlyph(keyword="Bars")
                self.SpinCenterCheck.setEnabled(False)
        if self.sender() == self.SpinCubeButton:
            if self.SpinCubeButton.isChecked():
                self.MomActors.ChangeSpinGlyph(keyword="Cubes")
                self.SpinCenterCheck.setEnabled(False)
        if self.sender() == self.SpinSphereButton:
            if self.SpinSphereButton.isChecked():
                self.MomActors.ChangeSpinGlyph(keyword="Spheres")
                self.SpinCenterCheck.setEnabled(False)
        if self.sender() == self.SpinArrowButton:
            if self.SpinArrowButton.isChecked():
                self.MomActors.ChangeSpinGlyph(keyword="Arrows")
                self.SpinCenterCheck.setChecked(False)
                self.SpinCenterCheck.setEnabled(True)
        if self.sender() == self.SpinConeButton:
            if self.SpinConeButton.isChecked():
                self.MomActors.ChangeSpinGlyph(keyword="Cones")
                self.SpinCenterCheck.setEnabled(False)
        if self.sender() == self.SpinCenterCheck:
            if self.SpinCenterCheck.isChecked():
                self.MomActors.ChangeSpinGlyph(keyword="CenterOn")
            else:
                self.MomActors.ChangeSpinGlyph(keyword="CenterOff")
        self.renWin.Render()
        return

    ##########################################################################
    # Function to change the shading of the glyphs that display the individual
    # magnetic moments
    ##########################################################################
    def ChangeShading(self):
        """Function to change the type of glyphs that display the individual magnetic
        moments

        Author
        ----------
        Anders Bergman, Jonathan Chico
        """
        if self.sender() == self.FlatShadeButton:
            if self.FlatShadeButton.isChecked():
                self.MomActors.ChangeSpinShade(keyword="Flat")

        if self.sender() == self.GouraudShadeButton:
            if self.GouraudShadeButton.isChecked():
                self.MomActors.ChangeSpinShade(keyword="Gouraud")

        if self.sender() == self.PBRShadeButton:
            if self.PBRShadeButton.isChecked():
                self.MomActors.ChangeSpinShade(keyword="PBR")

        if self.sender() == self.PhongShadeButton:
            if self.PhongShadeButton.isChecked():
                self.MomActors.ChangeSpinShade(keyword="Phong")

        self.renWin.Render()
        return

    ##########################################################################
    # Update the neighbours
    ##########################################################################

    def NeighbourControl(self):
        """
        Updates the neighbour actors with the current window, ASD data,
        general actors, render window, and mode.
        """
        self.NeighActors.UpdateNeighbour(
            window=self,
            ASDdata=self.ASDdata,
            ASDGenActors=self.ASDGenActors,
            renWin=self.renWin,
            mode=self.mode,
        )
        return

    ##########################################################################
    # Wrapper function to handle the camera functions
    ##########################################################################

    def camera_handler(self):
        """
        Handles various camera operations based on the sender of the signal.

        This method performs different camera-related actions such as resetting the camera,
        setting the camera view direction, updating the camera, and controlling the parallel scale.
        The specific action is determined by the sender of the signal.

        Actions:
        - Reset the camera to the original position if the sender is CamResetButton.
        - Set the camera view direction to X, Y, or Z axis if the sender is SetXView,
             SetYView, or SetZView respectively.
        - Update the camera if the sender is SetCamButton.
        - Change the parallel projection scale based on input from ParallelScaleLineEdit
            or ParallelScaleSlider.
        - Toggle parallel projections if the sender is ParallelProjectBox.
        """
        # -----------------------------------------------------------------------
        # Reset the camera to the original position
        # -----------------------------------------------------------------------
        if self.sender() == self.CamResetButton:
            if self.viz_type == "M":
                self.ASDCamera.reset_camera(
                    ren=self.ren, renWin=self.renWin, current_Actors=self.MomActors
                )
            elif self.viz_type == "N":
                self.ASDCamera.reset_camera(
                    ren=self.ren, renWin=self.renWin, current_Actors=self.NeighActors
                )
            elif self.viz_type == "E":
                self.ASDCamera.reset_camera(
                    ren=self.ren, renWin=self.renWin, current_Actors=self.EneActors
                )
        # -----------------------------------------------------------------------
        # Controlling what is up in the camera
        # -----------------------------------------------------------------------
        if self.sender() == self.SetXView:
            self.ASDCamera.set_Camera_viewUp(
                ren=self.ren, renWin=self.renWin, rdir=(1, 0, 0)
            )
        if self.sender() == self.SetYView:
            self.ASDCamera.set_Camera_viewUp(
                ren=self.ren, renWin=self.renWin, rdir=(0, 1, 0)
            )
        if self.sender() == self.SetZView:
            self.ASDCamera.set_Camera_viewUp(
                ren=self.ren, renWin=self.renWin, rdir=(0, 0, 1)
            )
        if self.sender() == self.SetCamButton:
            self.ASDCamera.Update_Camera(Window=self, ren=self.ren, renWin=self.renWin)
        # -----------------------------------------------------------------------
        # Controlling the parallel scale
        # -----------------------------------------------------------------------
        if self.sender() == self.ParallelScaleLineEdit:
            line = True
            slider = False
            self.ASDCamera.ChangeParallelProj(
                ren=self.ren,
                renWin=self.renWin,
                line=line,
                slider=slider,
                MainWindow=self,
            )
        if self.sender() == self.ParallelScaleSlider:
            line = False
            slider = True
            self.ASDCamera.ChangeParallelProj(
                ren=self.ren,
                renWin=self.renWin,
                line=line,
                slider=slider,
                MainWindow=self,
            )
        if self.sender() == self.ParallelProjectBox:
            self.ASDCamera.toggle_projections(
                renWin=self.renWin,
                window=self,
                ren=self.ren,
                checked=self.ParallelProjectBox.isChecked(),
            )
        if self.sender() == self.CamSaveButton:
            # Get and print the current camera settings
            # camera_settings = self.ASDCamera.get_camera_settings()
            self.ASDCamera.save_camera_settings()

        if self.sender() == self.CamLoadButton:
            # Get and print the current camera settings
            self.ASDCamera.load_camera_settings()
            # camera_settings = self.ASDCamera.get_camera_settings()
            self.ASDVizOpt.update_dock_info(current_Actors=self.MomActors, Window=self)

        return

    ##########################################################################
    # Wrapper to handle the clipper actions
    ##########################################################################
    def clipperHandler(self):
        """
        Handles the clipping operation for different visualization types.

        Depending on the value of `self.viz_type`, this method selects the appropriate
        actors (MomActors, NeighActors, or EneActors) and updates the clipper using
        the `ASDGenActors.UpdateClipper` method with the selected actors and other
        relevant parameters.
        """
        if self.viz_type == "M":
            current_Actors = self.MomActors
        if self.viz_type == "N":
            current_Actors = self.NeighActors
        if self.viz_type == "E":
            current_Actors = self.EneActors
        self.ASDGenActors.UpdateClipper(
            window=self,
            current_Actors=current_Actors,
            ASDVizOpt=self.ASDVizOpt,
            renWin=self.renWin,
            viz_type=self.viz_type,
        )
        return

    ##########################################################################
    # Function that calls the taking of a Snapshot of the current rendering window
    ##########################################################################

    def Snapshot(self):
        """
        Capture and save a screenshot of the current visualization.

        This method captures a screenshot of the current visualization using the ASDVizOpt.Screenshot method.
        It increments the number_of_screenshots counter after saving the screenshot.

        Args:
            None

        Returns:
            None
        """
        self.ASDVizOpt.Screenshot(
            renWin=self.renWin,
            number_of_screenshots=self.number_of_screenshots,
            png_mode=self.actionSave_png.isChecked(),
            pov_mode=self.actionSave_pov.isChecked(),
        )
        self.number_of_screenshots = self.number_of_screenshots + 1
        return

    ##########################################################################
    # Function that calls for updating the glyph resolutions
    ##########################################################################

    def Quality_control(self, value):
        """
        Updates the glyph quality in the visualization.
        """
        self.ASDVizOpt.GlyphQualityUpdate(
            window=self,
            value=value,
            viz_type=self.viz_type,
            mode=self.mode,
            renWin=self.renWin,
        )

    ##########################################################################
    # Function that calls for toggling FXAA
    ##########################################################################

    def FXAA_control(self, check):
        """
        Toggles the FXAA (Fast Approximate Anti-Aliasing) setting.
        """
        self.ASDVizOpt.toggle_FXAA(check=check, ren=self.ren, renWin=self.renWin)

    ##########################################################################
    # Function that calls for toggling surface texture
    ##########################################################################

    def Texture_control(self, check):
        """
        Toggles the texture control in the visualization options.
        """
        self.ASDTexture.toggle_Texture(
            check=check, actor=self.MomActors, texfile=self.texturefile
        )

    ##########################################################################
    # Function that calls for toggling ORM texture
    ##########################################################################

    def ORMTexture_control(self, check):
        """
        Toggles the ORM texture visualization based on the given check state.
        """
        self.ASDTexture.toggle_ORMTexture(
            check=check, actor=self.MomActors, texfile=self.ORMtexturefile
        )

    ##########################################################################
    # Function that calls for toggling ORM texture
    ##########################################################################

    def NTexture_control(self, check):
        """
        Toggles the NTexture visualization option.
        """
        self.ASDTexture.toggle_NTexture(
            check=check, actor=self.MomActors, texfile=self.Ntexturefile
        )

    ##########################################################################
    # Function that calls for toggling ORM texture
    ##########################################################################

    def ETexture_control(self, check):
        """
        Toggles the ETexture visualization option.
        """
        self.ASDTexture.toggle_ETexture(
            check=check, actor=self.MomActors, texfile=self.Etexturefile
        )

    ##########################################################################
    # Function that calls for toggling ORM texture
    ##########################################################################

    def ATexture_control(self, check):
        """
        Toggles the ATexture visualization option.
        """
        self.ASDTexture.toggle_ATexture(
            check=check, actor=self.MomActors, texfile=self.Atexturefile
        )

    ##########################################################################
    # Function that calls for toggling SSAO
    ##########################################################################

    def SSAO_control(self, check):
        """
        Toggle the SSAO (Screen Space Ambient Occlusion) control.
        """
        self.ASDVizOpt.toggle_SSAO(check=check, ren=self.ren)

    ##########################################################################
    # Function that calls for toggling shadows
    ##########################################################################
    # def Shadow_control(self, check):
    #    self.ASDVizOpt.toggle_Shadows(check=check,ren=self.ren, renWin=self.renWin)

    ##########################################################################
    # Function that calls for toggling HDRI
    ##########################################################################
    def HDRI_control(self, check):
        """
        Toggles HDRI visualization in the ASD visualization options.
        """
        self.ASDTexture.toggle_HDRI(
            check=check, ren=self.ren, renWin=self.renWin, hdrifile=self.hdrifile
        )
        return

    ##########################################################################
    # Function that calls for toggling skybox
    ##########################################################################
    def SkyBox_control(self, check):
        """
        Toggles the SkyBox visualization option.
        """
        self.ASDTexture.toggle_SkyBox(
            check=check, actor=self.MomActors, skyboxfile=self.hdrifile
        )

        return

    ##########################################################################
    # Finding the file name for the HDR file
    ##########################################################################
    def getHDRIFile(self):
        """
        Retrieves the HDRI file name and updates the UI elements based on the file's existence.
        """
        self.hdrifile = self.ASDTexture.getHDRIFileName(window=self)
        self.hdrifile_gotten = len(self.hdrifile) > 0
        if self.hdrifile_gotten:
            self.HDRICheck.setEnabled(True)
            self.SkyBoxCheck.setEnabled(True)
        return

    ##########################################################################
    # Finding the file name for the texture image
    ##########################################################################
    def getTextureFile(self):
        """
        Retrieves the texture file name and updates the UI accordingly.
        """
        self.texturefile = self.ASDTexture.getTextureFileName(window=self)
        self.texturefile_gotten = len(self.texturefile) > 0
        if self.texturefile_gotten:
            self.TextureCheck.setEnabled(True)
        return

    ##########################################################################
    # Finding the file name for the ORM texture image
    ##########################################################################
    def getORMTextureFile(self):
        """
        Retrieves the ORM texture file name and updates the UI accordingly.
        """
        self.ORMtexturefile = self.ASDTexture.getTextureFileName(window=self)
        self.ORMtexturefile_gotten = len(self.ORMtexturefile) > 0
        if self.ORMtexturefile_gotten:
            self.ORMTextureCheck.setEnabled(True)
        return

    ##########################################################################
    # Finding the file name for the normal texture image
    ##########################################################################
    def getNTextureFile(self):
        """
        Retrieves the texture file name and updates the UI accordingly.
        """
        self.Ntexturefile = self.ASDTexture.getTextureFileName(window=self)
        self.Ntexturefile_gotten = len(self.Ntexturefile) > 0
        if self.Ntexturefile_gotten:
            self.NTextureCheck.setEnabled(True)
        return

    ##########################################################################
    # Finding the file name for the anisotropy texture image
    ##########################################################################
    def getATextureFile(self):
        """
        Retrieves a texture file name and updates the UI accordingly.
        """
        self.Atexturefile = self.ASDTexture.getTextureFileName(window=self)
        self.Atexturefile_gotten = len(self.Atexturefile) > 0
        if self.Atexturefile_gotten:
            self.ATextureCheck.setEnabled(True)
        return

    ##########################################################################
    # Finding the file name for the emissive texture image
    ##########################################################################
    def getETextureFile(self):
        """
        Retrieves the texture file name and updates the UI accordingly.
        """
        self.Etexturefile = self.ASDTexture.getTextureFileName(window=self)
        self.Etexturefile_gotten = len(self.Etexturefile) > 0
        if self.Etexturefile_gotten:
            self.ETextureCheck.setEnabled(True)
        return

    ##########################################################################
    # Function that calls for toggling specular scattering
    ##########################################################################
    def RenSpecular_control(self, value):
        """
        Updates the specular rendering option with the given value.
        """
        self.MomActors.RenSpecularUpdate(value=value, renWin=self.renWin)

    ##########################################################################
    # Function that calls for toggling specular scattering
    ##########################################################################
    def RenSpecularPower_control(self, value):
        """
        Updates the specular power in the visualization options.
        """
        self.MomActors.RenSpecularPowerUpdate(value=value, renWin=self.renWin)

    ##########################################################################
    # Function that calls for toggling ambient scattering
    ##########################################################################
    def RenAmbient_control(self, value):
        """
        Updates the ambient rendering settings.
        """
        self.MomActors.RenAmbientUpdate(value=value, renWin=self.renWin)

    ##########################################################################
    # Function that calls for toggling diffuse scattering
    ##########################################################################
    def RenDiffuse_control(self, value):
        """
        Updates the rendering window with the given diffuse value.
        """
        self.MomActors.RenDiffuseUpdate(value=value, renWin=self.renWin)

    ##########################################################################
    # Function that calls for toggling PBR Emission value
    ##########################################################################
    def PBREmission_control(self, value):
        """
        Controls the PBR emission update with the given value.
        """
        self.MomActors.PBREmissionUpdate(value=value, ren=self.ren, renWin=self.renWin)

    ##########################################################################
    # Function that calls for toggling PBR Occlusion value
    ##########################################################################
    def PBROcclusion_control(self, value):
        """
        Controls the PBROcclusion update with the given value.
        """
        self.MomActors.PBROcclusionUpdate(value=value, ren=self.ren, renWin=self.renWin)

    ##########################################################################
    # Function that calls for toggling PBR Roughness value
    ##########################################################################
    def PBRRoughness_control(self, value):
        """
        Updates the PBR roughness value in the visualization options.
        """
        self.MomActors.PBRRoughnessUpdate(value=value, renWin=self.renWin)

    ##########################################################################
    # Function that calls for toggling PBR Roughness value
    ##########################################################################
    def PBRMetallic_control(self, value):
        """
        Updates the PBR metallic value in the visualization options.
        """
        self.MomActors.PBRMetallicUpdate(value=value, renWin=self.renWin)

    # --------------------------------------------------------------------------------
    # @brief Playback control for the animation of different movies, either for moments
    # or energies.
    # @details Playback control for the animation of different movies, either for moments
    # or energies. This allow full playback, pause and individual increase or decrease
    # of the images, it ensures that the images can only be increased/decreased
    # to the boundaries permitted by the data.
    # @author Jonathan Chico
    # --------------------------------------------------------------------------------
    def Playback_control(self):
        """
        Controls playback, pause, and navigation of the movie.
        """
        # -----------------------------------------------------------------------
        # Play the movie
        # -----------------------------------------------------------------------
        if self.sender() == self.PlayButton:
            if self.PlayButton.isChecked():
                self.iren.AddObserver("TimerEvent", self.Playback)
                self.timerId = self.iren.CreateRepeatingTimer(100)
                self.iren.SetStillUpdateRate(25.0)
                self.iren.SetDesiredUpdateRate(25.0)
                self.iren.Start()
            else:
                self.iren.DestroyTimer(self.timerId)
        # -----------------------------------------------------------------------
        # Pause the rendering
        # -----------------------------------------------------------------------
        if self.sender() == self.PauseButton:
            self.PlayButton.setChecked(False)
            self.iren.DestroyTimer(self.timerId)
        # -----------------------------------------------------------------------
        # Go to the previous image
        # -----------------------------------------------------------------------
        if self.sender() == self.previousButton:
            if self.current_time > 0:
                self.current_time -= 1
                self.UpdateImage()
        # -----------------------------------------------------------------------
        # Advance the current image
        # -----------------------------------------------------------------------
        if self.sender() == self.nextButton:
            if self.current_time < self.ASDdata.number_time_steps - 1:
                self.current_time += 1
                self.UpdateImage()
        return

    ##########################################################################
    # @brief Function to control the playback of the animation, whilst taking a snapshot
    # and updating the necessary data structures.
    # @details Function to control the playback of the animation, whilst taking a snapshot
    # and updating the necessary data structures. This function updates the magnetic
    # moments, the site dependent energy, as well as the timers to ensure that the
    # visualization finishes with the last image.
    # @author Jonathan Chico
    ##########################################################################
    def Playback(self, event, obj):
        """
        Handles playback events for updating visualization and taking snapshots.
        """
        if self.viz_type == "M":
            print("UpdateImage:", self.__class__.__name__)
            # -------------------------------------------------------------------
            # If the current time is the final time of the measurement destroy the
            # timer
            # -------------------------------------------------------------------
            if self.current_time == self.ASDdata.number_time_steps - 1:
                self.iren.DestroyTimer(self.timerId)
            # -------------------------------------------------------------------
            # Update the moments
            # -------------------------------------------------------------------
            self.MomActors.UpdateMoments(
                window=self,
                ASDdata=self.ASDdata,
                ASDGenActors=self.ASDGenActors,
                renWin=self.renWin,
            )
            # -------------------------------------------------------------------
            # Take a snapshot
            # -------------------------------------------------------------------
            # self.renWin.Render()
            self.Snapshot()
            # -------------------------------------------------------------------
            # increase the current timer
            # -------------------------------------------------------------------
            self.current_time += 1
        elif self.viz_type == "E":
            # -------------------------------------------------------------------
            # If the current time is the final time of the measurement destroy the
            # timer
            # -------------------------------------------------------------------
            if self.current_time == self.ASDdata.number_time_steps - 1:
                self.iren.DestroyTimer(self.timerId)
            # -------------------------------------------------------------------
            # Update the energy
            # -------------------------------------------------------------------
            self.EneActors.UpdateEnergy(
                window=self,
                ASDdata=self.ASDdata,
                ASDGenActors=self.ASDGenActors,
                renWin=self.renWin,
            )
            # -------------------------------------------------------------------
            # Take a snapshot
            # -------------------------------------------------------------------
            self.Snapshot()
            # -------------------------------------------------------------------
            # increase the current timer
            # -------------------------------------------------------------------
            self.current_time += 1
        return

    ##########################################################################
    # Individual update of the image, either by increasing the timer count
    # by one or by minus one
    ##########################################################################
    def UpdateImage(self):
        """
        Updates the visualization based on the current visualization type.
        """
        if self.viz_type == "M":
            print("UpdateImage:", self.__class__.__name__)
            # -------------------------------------------------------------------
            # Update the moments
            # -------------------------------------------------------------------
            self.MomActors.UpdateMoments(
                window=self,
                ASDdata=self.ASDdata,
                ASDGenActors=self.ASDGenActors,
                renWin=self.renWin,
            )
        elif self.viz_type == "E":
            # -------------------------------------------------------------------
            # Update the energy
            # -------------------------------------------------------------------
            self.EneActors.UpdateEnergy(
                window=self,
                ASDdata=self.ASDdata,
                ASDGenActors=self.ASDGenActors,
                renWin=self.renWin,
            )
        return

    ##########################################################################
    # Select the energy actor
    ##########################################################################
    def toggle_EneActor(self):
        """
        Toggles the visibility of EneActors based on the sender of the signal.
        """
        if self.sender() == self.EneDensButton and self.EneDensButton.isChecked():
            self.EneActors.EneDensActor.VisibilityOn()
            self.EneActors.EneActor.VisibilityOff()
            self.renWin.Render()
        if self.sender() == self.EneSiteGlyphs and self.EneSiteGlyphs.isChecked():
            self.EneActors.EneDensActor.VisibilityOff()
            self.EneActors.EneActor.VisibilityOn()
            self.renWin.Render()
        return

    ##########################################################################
    # Update the UI
    ##########################################################################
    def UpdateRenderer(self):
        """
        Update the renderer based on the state of the SpinsBox and SpinGlyphSelectBox.
        """
        if self.sender() == self.SpinsBox:
            if self.SpinsBox.isChecked():
                self.SpinGlyphSelectBox.setEnabled(True)
            else:
                self.SpinGlyphSelectBox.setEnabled(False)
        self.renWin.Render()
        return

    ##########################################################################
    # Interactive Simulations and Dock
    ##########################################################################
    def SetSDSliderValue(self, NSimulations):
        """
        Updates the slider value display with the number of simulations.
        """
        self.IntSDSliderVal.setText(f"Simulations: {10*NSimulations}")

    def SetMCSliderValue(self, NSimulations):
        """
        Sets the text of IntMCSliderVal to display the number of simulations.
        """
        self.IntMCSliderVal.setText(f"Simulations: {10*NSimulations}")

    def IntButtons(self):
        """
        Run a number of simulations depending slider inputs from the user using
        one of two simulation modes.
        """

        if self.sender() == self.IntSStepButton:
            ASDInteractiveTab.UpdateIntInputs(self)
            for _ in range(1 * self.IntSDSlider.value()):
                self.InteractiveVtk.S_Step()
        if self.sender() == self.IntMCMSimButton:
            ASDInteractiveTab.UpdateIntInputs(self)
            for __ in range(1 * self.IntMCSlider.value()):
                self.InteractiveVtk.M_step()
        if self.sender() == self.IntMCHSimButton:
            ASDInteractiveTab.UpdateIntInputs(self)
            for __ in range(1 * self.IntMCSlider.value()):
                self.InteractiveVtk.H_step()
        if self.sender() == self.IntResetButton:
            print("Reset button pressed")
            self.InteractiveVtk.Reset()
        if self.sender() == self.IntMomentButton:
            print("Moment button pressed")
            self.InteractiveVtk.read_moments()

    def UpdateInteractiveVtk(self):
        """Update text in the interactive window."""

        ASDInteractiveTab.UpdateIntInputs(self)

        self.InteractiveVtk.UpdateTemperature()
        self.InteractiveVtk.UpdateBfield()

    def InteractiveScreenshot(self):
        """
        Captures an interactive screenshot using the InteractiveVtk instance.
        """
        self.InteractiveVtk.Screenshot()

    def InteractiveScreenshotTic(self, tic):
        """
        Toggles the screenshot mode based on the provided tic value.
        """
        if tic:
            print("Taking screenshots")
            self.InteractiveVtk.film = True
        else:
            print("Not taking screenshots")
            self.InteractiveVtk.film = False

    def CheckForInteractorFiles(self):
        """
        Check if we have any input/output files and determines if
        the interactive window is ready to launch. Shows error window
        with missing files if check fails.

        Returns:
                Check   :   Bool
        """

        restartfile, coordfile = "dummystring", "dummystring"
        Check = False

        # if len(glob.glob("restart.????????.out")) > 0:
        #     restartfile = glob.glob("restart.????????.out")[0]
        # if len(glob.glob("coord.????????.out")) > 0:
        #     coordfile = glob.glob("coord.????????.out")[0]
        if len(self.ASDInputGen.posfile) == 0 and path.exists("inpsd.dat"):
            posfile, momfile = self.ASDInputGen.GetPosMomFiles()
            self.ASDInputGen.posfile = glob.glob(posfile)[0]
            print("posfile:", self.ASDInputGen.posfile)
        if len(self.ASDInputGen.momfile) == 0 and path.exists("momfile"):
            self.ASDInputGen.momfile = glob.glob(momfile)[0]
            print("momfile:", self.ASDInputGen.momfile)

        InputChecklist = [
            path.exists("inpsd.dat"),
            path.exists(self.ASDInputGen.posfile),
            path.exists(self.ASDInputGen.momfile),
        ]

        if all(x for x in InputChecklist):
            Check = True

        # Error message
        Files = ["inpsd.dat", "posfile", "momfile"]
        MissingFiles = ", ".join(
            [file for index, file in enumerate(Files) if not InputChecklist[index]]
        )

        if not Check:
            self.InteractiveErrorWindow = ASDInputWindows.ErrorWindow()
            self.InteractiveErrorWindow.FunMsg.setText("Task failed successfully!")
            self.InteractiveErrorWindow.ErrorMsg.setText(
                f"Could not launch interactive simulation. Missing input files: {MissingFiles}"
            )
            self.InteractiveErrorWindow.show()
        return Check

    ###########################################################################

    def RunSimulation(self):
        """
        Run simulation using the uppasd module. Checks if a inpsd.dat file
        is present in current directory and asks for overwrite.

        Author: Erik Karpelin
        """

        # import uppasd as asd

        if not path.isfile("inpsd.dat"):
            print("inpsd.dat not found, creating from asd_gui")
            self.WriteInputFile()

        if self.ASDsim is not None:
            print("Running simulation from RunSimulation")
            self.ASDsim.init_simulation()
            self.ASDsim.run_simulation()

        return

    def SetStructureTemplate(self, structure):
        """Relay function to handle the structure templates."""
        self.ASDInputGen.SetStructureTemplate(self, structure)
        return

    def ResetInputs(self):
        """Relay function to handle the reset button."""
        self.ASDInputGen.ResetInputs(self)
        return

    def MagnonQuickSetup(self):
        """Relay function for magnon quick setup"""
        self.ASDInputGen.MagnonQuickSetup(self)

    def ImportSystem(self):
        """Relay function for importing system"""
        self.ASDInputGen.import_system(self)
