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

import ASD_GUI.Input_Creator.ASDInputGen as ASDInputGen
from enum import Enum
from PyQt6 import QtWidgets
from PyQt6.QtCore import Qt
from PyQt6.QtWidgets import QMainWindow

class Backend(Enum):
    UppASD_VTK = 1
    UppASD_MAT = 2
    UppASD_INP = 3
    UppASD_INT = 4
################################################################################
# @brief Class that defines the main window where all the widgets and rendering take place.
# @details It controls the actions which take place in the GUI. It defined the main window
# that allows for the following features:
# - Visual input file creation.
# - Matplotlib plotting of several key \c UppASD outputs.
# - VTK rendering of 3D \c UppASD data.
# @author Jonathan Chico
################################################################################


class UppASDVizMainWindow(QMainWindow):
    ############################################################################
    # @brief Class constructor for the main window
    # @details Class constructor for the main waindow. It initializes the inclusion
    # of several auxiliary classes that are used to setup the GUI functionality.
    # @author Jonathan Chico
    ############################################################################
    def __init__(self):
        from ASD_GUI.PLOT import ASDPlots2D
        from ASD_GUI.UI import ASDInputWindows
        from ASD_GUI.VTK_Viz import ASDVTKReading
        from ASD_GUI.PLOT import ASDPlotsReading
        from ASD_GUI.VTK_Viz import ASDVTKGenActors
        from ASD_GUI.VTK_Viz import ASDVTKVizOptions
        import ASD_GUI.UI.ASDInteractiveTab as ASDInteractiveTab
        import ASD_GUI.ASD_Interactive.interactiveASD as IntASD
        from matplotlib.figure import Figure
        from mpl_toolkits.mplot3d import Axes3D
        from vtk import vtkOpenGLRenderer, vtkInteractorStyleTrackballCamera
        from matplotlib.backends.backend_qtagg import FigureCanvas
        # from matplotlib.backends.backend_qt5agg import FigureCanvas
        from vtk.qt.QVTKRenderWindowInteractor import QVTKRenderWindowInteractor
        super(UppASDVizMainWindow, self).__init__()
        # -----------------------------------------------------------------------
        # Define the array containing the file names necessary for visualization
        # -----------------------------------------------------------------------
        self.file_names = [None]*6
        self.plotfile_names = [None]*7
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
        self.ASDVizOpt = ASDVTKVizOptions.ASDVizOptions()
        self.ASDGenActors = ASDVTKGenActors.ASDGenActors()
        self.ASDPlotData = ASDPlotsReading.ReadPlotData()
        self.ASDPlots2D = ASDPlots2D.Abstract2DPlot()
        self.ASDCorrelationPlots = ASDPlots2D.Correlation_Plots()
        self.ASDInputGen = ASDInputGen.ASDInputGen()
        self.Restart_Window = ASDInputWindows.Restart_Window()
        self.Posfile_Window = ASDInputWindows.Posfile_Window()
        self.Momfile_Window = ASDInputWindows.Momfile_Window()
        self.InitPhase_Window = ASDInputWindows.InitPhase_Window()
        self.Jfile_Window = ASDInputWindows.Jfile_Window()
        self.DMfile_Window = ASDInputWindows.DMfile_Window()
        self.Kfile_Window = ASDInputWindows.Kfile_Window()
        self.InteractiveDockWidget = ASDInteractiveTab.InteractiveDock(self)
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
        # Initialize and setup the necessary structures for the UI
        # -----------------------------------------------------------------------
        self.SetupUI()
        self.InitUI()
        # -----------------------------------------------------------------------
        # Set the Plotting Canvas for the matplotlib plots
        # -----------------------------------------------------------------------
        self.InitPlotUI()
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
        self.Plotting_Figure3D.add_subplot(111, projection='3d')
        self.Plotting_Figure3D.patch.set_alpha(0.0)
        self.Plotting_ax3D = self.Plotting_Figure3D.axes[0]
        self.Plot3DLayout.addWidget(self.Plotting_canvas3D)

        # Interactive object
        self.InteractiveVtk = IntASD.InteractiveASD(self.Intren, self.IntrenWin, self.Intiren)

        return
    ############################################################################
    # @brief Wrapper for the writing of the input file
    # @details This function first create the dictionary of key words, this is then
    # populated with the parameters from the GUI. The dictionary is then cleaned
    # and written to file.
    # @author Jonathan Chico
    ############################################################################

    def WriteInputFile(self):
        self.ASDInputGen.ASDSetDefaults()
        self.ASDInputGen.ASDInputGatherer(self)
        self.ASDInputGen.clean_var()
        self.ASDInputGen.write_inpsd()
        return

    ############################################################################
    # @brief Initialize the UI and set the relevant actions
    # @details Initialize the UI and set the relevant actions. Defines the Toolbars
    # and calls for their initialization and population, as well as the reading of the
    # .ui file defining the properties of the window.
    # Also sets up several validators to forbid erroneous data to be fed into the
    # GUI.
    # @author Jonathan Chico
    ############################################################################

    def SetupUI(self):
        import os
        from PyQt6 import uic
        from PyQt6.QtGui import QIntValidator, QDoubleValidator
        from PyQt6.QtWidgets import QToolBar, QVBoxLayout
        from ASD_GUI.UI.ASDMenuToolbar import VTK_Menu_and_Toolbar_Setup, Plot_Menu_and_Toolbar_Setup, Input_Toolbar_Setup, InteractiveDock_Setup
        self.VTKToolBar = QToolBar()
        self.MatPlotToolbar = QToolBar()
        self.InputToolbar = QToolBar()
        # -----------------------------------------------------------------------
        # Set up UI from Designer file
        # -----------------------------------------------------------------------
        path = os.path.dirname(os.path.abspath(__file__))
        uic.loadUi(os.path.join(path, 'ASD_Viz.ui'), self)
        self.chooseBackend()
        self.ModeSelector.currentChanged.connect(self.chooseBackend)
        Plot_Menu_and_Toolbar_Setup(self)
        VTK_Menu_and_Toolbar_Setup(self)
        Input_Toolbar_Setup(self)
        InteractiveDock_Setup(self)
        self.CorrOptsBox.setEnabled(True)
        self.NeighValidator = QIntValidator()
        self.IntegerValidator = QIntValidator()
        self.IntegerValidator.setRange(0, 99999999)
        self.PosDoubleValidator = QDoubleValidator()
        self.PosDoubleValidator.setRange(0, 99999999.9999)
        self.PosDoubleValidator.setDecimals(10)
        self.DoubleValidator = QDoubleValidator()
        self.DoubleValidator.setDecimals(10)
        self.ASDInputGen.ASDInputConstrainer(self)
        self.InpSqwSCStep.setValidator(self.IntegerValidator)
        self.InpPlotDt.setValidator(self.PosDoubleValidator)
        self.AMSDisplayLayout = QVBoxLayout()
        self.AMSDisplayOpts.setLayout(self.AMSDisplayLayout)
        self.ResetUI()
        self.ASDInputGen.ASDSetDefaults()
        self.Posfile_Window.InpPosDone.clicked.connect(self.update_names)
        self.Momfile_Window.InpMomDone.clicked.connect(self.update_names)
        self.Jfile_Window.InpJfileDone.clicked.connect(self.update_names)
        self.DMfile_Window.InpDMfileDone.clicked.connect(self.update_names)
        self.Jfile_Window.InJfileGenVectors.clicked.connect(lambda: 
                                        self.Jfile_Window.GenerateVectorsFromCell(self))
        self.DMfile_Window.InDMfileGenVectors.clicked.connect(lambda: 
                                        self.DMfile_Window.GenerateVectorsFromCell(self))
        self.Restart_Window.InpRestAppendButton.clicked.connect(
            self.create_restart)
        self.Restart_Window.InpRestartDone.clicked.connect(self.create_restart)
        self.Restart_Window.InpRestartDone.clicked.connect(self.update_names)
        self.InitPhase_Window.InitPhaseDoneButton.clicked.connect(
            self.getInitPhase)
        return
    ############################################################################
    # @brief Wrapper to create the restartfile
    # @author Jonathan Chico
    ############################################################################

    def create_restart(self):
        self.Restart_Window.write_restartfile(self.ASDInputGen)
        return
    ############################################################################
    # Choose which kind of backend one will use to display VTK based visualizations,
    # matplotlib based visualizations
    ############################################################################

    def chooseBackend(self):
        import ASD_GUI.UI.ASDInteractiveTab as ASDInteractive
        if self.ModeSelector.currentIndex() == 0:
            self.backend = Backend.UppASD_VTK
            if not self.VTKWidgetPresent:
                self.VTKWidget_Layout.addWidget(self.vtkWidget)
                self.VTKWidgetPresent = True
        elif self.ModeSelector.currentIndex() == 1:
            self.backend = Backend.UppASD_MAT
            self.vtkWidget.setVisible(False)
        elif self.ModeSelector.currentIndex() == 2:
            self.backend = Backend.UppASD_INP
            self.vtkWidget.setVisible(False)
        elif self.ModeSelector.currentIndex() == 3:
            self.backend = Backend.UppASD_INT
            if not self.IntWidgetPresent:
                self.InteractiveWidget_Layout.addWidget(self.IntVtkWidget)
                self.IntWidgetPresent = True
            if self.CheckForInteractorFiles() and not self.IntLaunched:
                ASDInteractive.InitializeInteractor(self, self.InteractiveVtk)
                self.IntLaunched = True
        self.ResetUI()
        return
    ############################################################################
    # Reset the UI to change between the VTK based visualization and the matplotlib
    # based visualization
    ############################################################################

    def ResetUI(self):
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
    ############################################################################
    # Initialization of some of the UI properties
    ############################################################################

    def InitUI(self):
        self.EneMainBox.setEnabled(False)
        self.CamMainBox.setEnabled(False)
        self.MagMainGroup.setEnabled(False)
        self.NeighMainBox.setEnabled(False)
        self.SceneOptMainBox.setEnabled(False)
        self.SpinGlyphSelectBox.setEnabled(True)
        self.ClippBox.setEnabled(False)
        self.ClippBox.setChecked(False)
        self.ClusBox.setVisible(False)
        self.KMCCheck.setVisible(False)
        self.SceneOptMainBox.setEnabled(True)
        self.CamMainBox.setEnabled(True)
        self.actionSave_pov.setEnabled(True)
        self.actionSave_png.setEnabled(True)
        self.ClippBox.setEnabled(True)
        self.actionDisplayMagDens.setEnabled(False)
        self.actionX_ProjMagDens.setEnabled(False)
        self.actionY_ProjMagDens.setEnabled(False)
        self.actionZ_ProjMagDens.setEnabled(False)
        self.file_names[0] = self.ASDdata.posfiles
        self.file_names[1] = self.ASDdata.magnetization
        self.file_names[2] = self.ASDdata.kmcfiles
        self.file_names[3] = self.ASDdata.structfiles
        self.file_names[4] = self.ASDdata.enefiles
        self.file_names[5] = self.ASDdata.dmdatafiles
        self.ProgressBar.setValue(0)
        self.ProgressLabel.setText(
            '   {:}%'.format(int(self.ProgressBar.value())))
        return
    ############################################################################
    # Initialization of some of the UI properties for 2D plots
    ############################################################################

    def InitPlotUI(self):
        self.plotfile_names[0] = self.ASDPlotData.yamlfile
        self.plotfile_names[1] = self.ASDPlotData.amsfile
        self.plotfile_names[2] = self.ASDPlotData.sqwfile
        self.plotfile_names[3] = self.ASDPlotData.averages
        self.plotfile_names[4] = self.ASDPlotData.trajectory
        self.plotfile_names[5] = self.ASDPlotData.totenergy
        self.plotfile_names[6] = self.ASDPlotData.qfile
        self.SqwProjBox.setEnabled(False)
        self.SqwColorMapSelect.setEnabled(False)
        self.AveOpts.setEnabled(False)
        self.EneOpts.setEnabled(False)
        self.AMSDisplayOpts.setVisible(False)
        return
    ############################################################################
    # Finding the file name for the VTK plots
    ############################################################################

    def getFile(self):
        self.ASDdata.getFileName(window=self)
        return
    ############################################################################
    # Finding the file name for the matplotlib plots
    ############################################################################

    def getPlotFile(self):
        self.ASDPlotData.getFileName(window=self)
        return
    ############################################################################
    # Finding the file name for the input file generation
    ############################################################################

    def getInpFile(self):
        self.ASDInputGen.getFileName(window=self)
        return
    ############################################################################
    ############################################################################

    def getInitPhase(self):
        if self.sender() == self.InitPhase_Window.InitPhaseDoneButton:
            self.init_phase_data = self.InitPhase_Window.init_phase_data
        return
    ############################################################################
    # @brief Open auxiliary windows for the inputfile creation GUI
    ############################################################################

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
            self.Posfile_Window.posfile_gotten = False
            # if self.InpCheckRandAlloy.isChecked():
            #     self.Posfile_Window.InPosBox.setEnabled(False)
            #     self.Posfile_Window.InPosBox.setVisible(False)
            #     self.Posfile_Window.InPosBoxRand.setEnabled(True)
            #     self.Posfile_Window.InPosBoxRand.setVisible(True)
            self.Posfile_Window.CheckForFile(self)
            self.Posfile_Window.InPosBox.setEnabled(True)
            self.Posfile_Window.InPosBox.setVisible(True)
            self.Posfile_Window.InPosBoxRand.setEnabled(False)
            self.Posfile_Window.InPosBoxRand.setVisible(False)
            self.Posfile_Window.show()
        if self.sender() == self.InpMomButtonCreate:
            self.Momfile_Window.momfile_gotten = False
            self.Momfile_Window.CheckForFile(self)
            self.Momfile_Window.show()
        if self.sender() == self.InpJfileButtonCreate:
            self.Jfile_Window.jfile_gotten = False
            self.Jfile_Window.CheckForFile(self)
            self.Jfile_Window.show()
        if self.sender() == self.InpDMButtonCreate:
            self.DMfile_Window.DMfile_gotten = False
            self.DMfile_Window.CheckForFile(self)
            self.DMfile_Window.show()
        if self.sender() == self.InpKfileButtonCreate:
            self.Kfile_Window.Kfile_gotten = False
            self.Kfile_Window.CheckForFile(self)
            self.Kfile_Window.show()
        if self.sender() == self.InpSetPhases:
            if self.InpInitLLG.isChecked():
                self.InitPhase_Window.IpNphaseBox.setEnabled(True)
                self.InitPhase_Window.IpNphaseBox.setVisible(True)
                self.InitPhase_Window.MCannealBox.setEnabled(False)
                self.InitPhase_Window.MCannealBox.setVisible(False)
                self.InitPhase_Window.IpVPOBox.setVisible(False)
                self.InitPhase_Window.IpVPOBox.setEnabled(False)
                self.InitPhase_Window.InitPhaseAddButton.setEnabled(True)
                self.InitPhase_Window.InitPhaseDelButton.setEnabled(True)
            if self.InpInitMcMet.isChecked() or self.InpInitMcHeat.isChecked():
                self.InitPhase_Window.IpNphaseBox.setEnabled(False)
                self.InitPhase_Window.IpNphaseBox.setVisible(False)
                self.InitPhase_Window.MCannealBox.setEnabled(True)
                self.InitPhase_Window.MCannealBox.setVisible(True)
                self.InitPhase_Window.IpVPOBox.setVisible(False)
                self.InitPhase_Window.IpVPOBox.setEnabled(False)
                self.InitPhase_Window.InitPhaseAddButton.setEnabled(True)
                self.InitPhase_Window.InitPhaseDelButton.setEnabled(True)
            if self.InpInitVPO.isChecked():
                self.InitPhase_Window.IpNphaseBox.setEnabled(False)
                self.InitPhase_Window.IpNphaseBox.setVisible(False)
                self.InitPhase_Window.MCannealBox.setEnabled(False)
                self.InitPhase_Window.MCannealBox.setVisible(False)
                self.InitPhase_Window.IpVPOBox.setVisible(True)
                self.InitPhase_Window.IpVPOBox.setEnabled(True)
                self.InitPhase_Window.InitPhaseAddButton.setEnabled(False)
                self.InitPhase_Window.InitPhaseDelButton.setEnabled(False)

            self.InitPhase_Window.show()
        return
    ############################################################################
    ############################################################################

    def update_names(self):
        self.ASDInputGen.update_file_name(window=self)
        return
    ############################################################################
    # @brief Function to determine if the restartfile can be created.
    # @details This function will test if the lattice vectors have been defined,
    # as well as the posfile and momfile, that is everything which is necessary to
    # generate a restartfile
    # @author Jonathan Chico
    ############################################################################

    def check_for_restart(self):
        import ASD_GUI.UI.ASDInputWindows as ASDInputWindows

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
            if (self.Momfile_Window.momfile_gotten or self.ASDInputGen.momfile_gotten) and\
                    (self.Posfile_Window.posfile_gotten or self.ASDInputGen.posfile_gotten):
                self.Restart_Window.restart_pre_generation(self.ASDInputGen)
                self.Restart_Window.show()
            else:
                self.Create_Restart_Error_Window = ASDInputWindows.Error_Window()
                self.Create_Restart_Error_Window.FunMsg.setText(
                    "I'm sorry, Dave. I'm afraid I can't do that.")
                self.Create_Restart_Error_Window.ErrorMsg.setText(
                    "Error: A posfile and a momfile must be first be defined")
                self.Create_Restart_Error_Window.show()
                print('Error: A posfile and a momfile must be first be defined')
        else:
            self.Create_Restart_Error_Window = ASDInputWindows.Error_Window()
            self.Create_Restart_Error_Window.FunMsg.setText(
                "I'm sorry, Dave. I'm afraid I can't do that.")
            self.Create_Restart_Error_Window.ErrorMsg.setText(
                "Error: The unit cell vectors need to be defined")
            self.Create_Restart_Error_Window.show()
            print('Error: The unit cell vectors need to be defined')
        return
    ############################################################################
    # Function to select the appropriate data to plot
    ############################################################################

    def PlottingSelector(self):
        from PyQt6.QtCore import QSignalBlocker
        # -----------------------------------------------------------------------
        # Plot the spin-spin correlation function
        # -----------------------------------------------------------------------
        if self.sender() == self.actionS_q_w:
            self.plotting_mode = 'correlation'
            self.MatToolBox.setCurrentIndex(0)
            self.PlotStacked.setCurrentIndex(0)
            if not self.ASDPlotData.not_read_ams:
                self.AMSDispCheckBox.setChecked(True)
                QSignalBlocker(self.AMSDispCheckBox)
            if not self.ASDPlotData.not_read_sqw:
                self.SqwDispCheckBox.setChecked(True)
                QSignalBlocker(self.SqwDispCheckBox)
        # -----------------------------------------------------------------------
        # Plot the averages
        # -----------------------------------------------------------------------
        if self.sender() == self.actionAverages:
            self.plotting_mode = 'averages'
            self.MatToolBox.setCurrentIndex(2)
            self.PlotStacked.setCurrentIndex(0)
        # -----------------------------------------------------------------------
        # Plot the energies
        # -----------------------------------------------------------------------
        if self.sender() == self.actionTotEnergy:
            self.plotting_mode = 'energy'
            self.MatToolBox.setCurrentIndex(3)
            self.PlotStacked.setCurrentIndex(0)
            # -------------------------------------------------------------------
            # Check if the 2D axis exists if it does not create it
            # -------------------------------------------------------------------
        if self.sender() == self.actionTrajectory:
            self.plotting_mode = 'trajectory'
            self.MatToolBox.setCurrentIndex(1)
            self.PlotStacked.setCurrentIndex(1)
        self.Plotting_Figure.canvas.draw()
        self.InitPlotUI()
        self.ASDPlotData.PlotReadingWrapper(self.plotfile_names, self)
        if not self.ASDPlotData.not_read_ams:
            self.ams_data_x = self.ASDPlotData.ams_data_x
            self.ams_data_y = self.ASDPlotData.ams_data_y
            self.ams_label = self.ASDPlotData.ams_label
        self.InitPlotUI()
        self.set_ams_checkboxes()
        return
    ############################################################################
    # @brief Function for the creation of checkboxes for the ams display
    # @details This should allow for the dynamical creation of checkboxes for each
    # branch in the ams. It also connects it to a function that prunes the data
    # so that it can be selectively plotted.
    # @author Jonathan Chico
    ############################################################################

    def set_ams_checkboxes(self):
        from PyQt6.QtWidgets import QCheckBox

        self.AMSCheckboxes = dict()
        for ii in reversed(range(self.AMSDisplayLayout.count())):
            self.AMSDisplayLayout.itemAt(ii).widget().setParent(None)
        if self.ASDPlotData.ams_file_present:
            # -------------------------------------------------------------------
            # Create checkboxes for the AMS branches
            # -------------------------------------------------------------------
            for ii in range(0, len(self.ASDPlotData.ams_data_y)):
                name = 'ams_branch_{}'.format(ii)
                checkbox = QCheckBox()
                checkbox.setObjectName(name)
                checkbox.setText('Display Branch {: 4d}'.format(ii+1))
                checkbox.setChecked(True)
                checkbox.toggled.connect(self.AMS_PrunePlot)
                self.AMSDisplayLayout.addWidget(checkbox)
                self.AMSCheckboxes[name] = checkbox
        return
    ############################################################################
    # Select the projection of the S(q,w)
    ############################################################################

    def SQW_Proj_Select(self):
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
    ############################################################################
    # Select the colormap over which the S(q,w) will be plotted
    ############################################################################

    def Sqw_ColorMapSelector(self):
        if self.sender() == self.SqwCoolwarm:
            self.plot2D_cmap_indx = 0
        if self.sender() == self.SqwSpectral:
            self.plot2D_cmap_indx = 1
        if self.sender() == self.SqwBlackbody:
            self.plot2D_cmap_indx = 2
        self.PlottingWrapper()
        return
    ############################################################################
    # Plotting the directions of the magnetization
    ############################################################################

    def PlotMagDirSelector(self):
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
    ############################################################################
    # Changing the marker size of the lines
    ############################################################################

    def PlotLineChanger(self, value):
        self.ASDPlots2D.linewidth = value/2.0
        self.PlottingWrapper()
    ############################################################################
    # Changing the marker size of the lines
    ############################################################################

    def PlotMarkerChanger(self, value):
        self.ASDPlots2D.markersize = value/2.0
        self.PlottingWrapper()
    ############################################################################
    # Changing the marker size of the lines
    ############################################################################

    def PlotXGridToggle(self):
        self.ASDPlots2D.xgrid = not self.ASDPlots2D.xgrid
        self.PlottingWrapper()
    ############################################################################
    # Changing the marker size of the lines
    ############################################################################

    def PlotYGridToggle(self):
        self.ASDPlots2D.ygrid = not self.ASDPlots2D.ygrid
        self.PlottingWrapper()
    ############################################################################
    #  Toggling SQW grid lines on/off
    ############################################################################

    def PlotSQWGridToggle(self):
        self.ASDCorrelationPlots.grid = not self.ASDCorrelationPlots.grid
        self.PlottingWrapper()
    ############################################################################
    #  Toggling SQW grid lines on/off
    ############################################################################

    def PlotAMSGridToggle(self):
        self.ASDPlots2D.amsgrid = not self.ASDPlots2D.amsgrid
        self.PlottingWrapper()
    ############################################################################
    # Changing the width of S(q,w) plots
    ############################################################################

    def SqwWidthChanger(self, value):
        self.ASDCorrelationPlots.sigma_w = self.ASDCorrelationPlots.w_min*value
        self.ABCorrWidthTX.setText(
            f'{self.ASDCorrelationPlots.w_min*value:.3f}')
        self.PlottingWrapper()
    ############################################################################
    # Plotting the components of the energy
    ############################################################################

    def PlotEneCompSelector(self):
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
    ############################################################################
    ############################################################################

    def ToggleInitPhase(self):
        from ASD_GUI.UI.ASDMenuToolbar import UpdateUI
        UpdateUI(self)
        return

    def ToggleHessians(self):
        from ASD_GUI.UI.ASDMenuToolbar import UpdateUI
        UpdateUI(self)
        return
    ############################################################################
    # @brief Function to selective plot the ams branches
    # @details Function to selectively plot the ams branches. It functions by
    # finding which of the checkboxes identifying each branch is selected
    # after this it creates a new data set that contains only the necessary data
    # @author Jonathan Chico
    ############################################################################

    def AMS_PrunePlot(self):
        self.ams_data_y = []
        self.ams_data_x = []
        self.ams_label = []
        for ii in range(len(self.AMSCheckboxes)):
            name = 'ams_branch_{}'.format(ii)
            if self.AMSCheckboxes[name].isChecked():
                self.ams_data_x.append(self.ASDPlotData.ams_data_x[ii])
                self.ams_data_y.append(self.ASDPlotData.ams_data_y[ii])
                self.ams_label.append(self.ASDPlotData.ams_label[ii])
        self.PlottingWrapper()
        return
    ############################################################################
    # @brief Wrapper function that takes care of plotting the selected plot
    # @details Wrapper function that takes care of plotting the selected plot, it allows
    # the user to choose between the following different types of plots
    #   - Spin-Spin correlation functions
    #       - S(q,w)
    #       - AMS
    #   - Magnetization averages
    #   - Single spin trajectories
    # @author Jonathan Chico
    ############################################################################

    def PlottingWrapper(self):
        """ Wrapper function that takes care of plotting the selected plot, it allows
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
        from ASD_GUI.PLOT import ASDPlots2D
        from ASD_GUI.UI import ASDInputWindows
        # -----------------------------------------------------------------------
        # Plotting the spin-spin correlation function
        # -----------------------------------------------------------------------
        if self.sender() == self.AMSDispCheckBox or self.sender() == self.SqwDispCheckBox:
            self.plotting_mode = 'correlation'
            if self.ASDPlotData.not_read_sqw or self.ASDPlotData.not_read_ams:
                self.ASDPlotData.PlotReadingWrapper(self.plotfile_names, self)
                if self.AMSDispCheckBox.isChecked():
                    self.ams_data_x = self.ASDPlotData.ams_data_x
                    self.ams_data_y = self.ASDPlotData.ams_data_y
                    self.ams_label = self.ASDPlotData.ams_label
                    self.set_ams_checkboxes()

        if self.plotting_mode == 'correlation':
            # -------------------------------------------------------------------
            # Perform the actual plotting
            # -------------------------------------------------------------------
            self.SqwProjBox.setEnabled(True)
            self.SqwColorMapSelect.setEnabled(True)
            self.SqwDisplayOpts.setEnabled(True)
            # -------------------------------------------------------------------
            # Plotting the S(q,w)
            # -------------------------------------------------------------------
            if self.SqwDispCheckBox.isChecked() and not self.AMSDispCheckBox.isChecked():
                if self.ASDPlotData.sqw_file_present:
                    self.ASDCorrelationPlots.Sqw_Plot(self.Plotting_ax, self.ASDPlotData.sqw_data,
                                                      self.SQW_proj_indx, self.ASDPlotData.sqw_labels, self.plot2D_cmap_indx,
                                                      self.ASDPlotData.ax_limits, self.ASDPlotData.q_labels, self.ASDPlotData.q_idx)
                    self.AMSDisplayOpts.setVisible(False)
                    self.AMSDisplayOpts.setEnabled(False)
                else:
                    self.sqw_Error_Window = ASDInputWindows.Error_Window()
                    self.sqw_Error_Window.FunMsg.setText(
                        "I'm sorry, Dave. I'm afraid I can't do that.")
                    self.sqw_Error_Window.ErrorMsg.setText(
                        "Error: No 'sqw.*.out' file.")
                    self.sqw_Error_Window.show()
                    self.SqwDispCheckBox.setChecked(False)
                    print("No 'sqw.*.out' file.")
            # -------------------------------------------------------------------
            # Plotting the AMS
            # -------------------------------------------------------------------
            elif self.AMSDispCheckBox.isChecked() and not self.SqwDispCheckBox.isChecked():
                if self.ASDPlotData.ams_file_present:
                    self.ASDPlots2D.LinePlot(self.Plotting_ax, self.ams_data_x,
                                             self.ams_data_y, self.ams_label, self.ASDPlotData.ams_ax_label,
                                             tick_labels=self.ASDPlotData.q_labels, tick_idx=self.ASDPlotData.q_idx)
                    self.AMSDisplayOpts.setVisible(True)
                    self.AMSDisplayOpts.setEnabled(True)
                else:
                    self.ams_Error_Window = ASDInputWindows.Error_Window()
                    self.ams_Error_Window.FunMsg.setText(
                        "I'm sorry, Dave. I'm afraid I can't do that.")
                    self.ams_Error_Window.ErrorMsg.setText(
                        "Error: No 'ams.*.out' file.")
                    self.ams_Error_Window.show()
                    self.AMSDispCheckBox.setChecked(False)
                    print("No 'ams.*.out' file.")
            # -------------------------------------------------------------------
            # Plotting the S(q,w) and the AMS
            # -------------------------------------------------------------------
            if self.SqwDispCheckBox.isChecked() and self.AMSDispCheckBox.isChecked():
                if self.ASDPlotData.sqw_file_present and self.ASDPlotData.ams_file_present:
                    self.ASDCorrelationPlots.AMS_Sqw_Plot(self.Plotting_ax, self.ASDPlotData.sqw_data,
                                                          self.SQW_proj_indx, self.ASDPlotData.sqw_labels, self.ams_data_x,
                                                          self.ams_data_y, self.ASDPlotData.hf_scale, self.plot2D_cmap_indx,
                                                          self.ASDPlotData.ax_limits, self.ASDPlotData.q_labels, self.ASDPlotData.q_idx)
                    self.AMSDisplayOpts.setVisible(True)
                    self.AMSDisplayOpts.setEnabled(True)
                else:
                    self.ams_sqw_Error_Window = ASDInputWindows.Error_Window()
                    self.ams_sqw_Error_Window.FunMsg.setText(
                        "I'm sorry, Dave. I'm afraid I can't do that.")
                    self.ams_sqw_Error_Window.ErrorMsg.setText(
                        "Error: No 'ams.*.out' or 'sqw.*.out' file.")
                    self.ams_sqw_Error_Window.show()
                    print("No 'ams.*.out' or 'sqw.*.out' file.")
                    self.SqwProjBox.setEnabled(False)
                    self.SqwColorMapSelect.setEnabled(False)
                    self.AMSDispCheckBox.setChecked(False)
                    self.SqwDispCheckBox.setChecked(False)
        # -----------------------------------------------------------------------
        # Plotting the average magnetization
        # -----------------------------------------------------------------------
        if self.plotting_mode == 'averages':
            self.AveOpts.setEnabled(True)
            if self.ASDPlotData.ave_file_present:
                curr_data_x = []
                curr_data_y = []
                curr_labels = []
                for ii in range(len(self.MagDirIndx)):
                    curr_data_x.append(
                        self.ASDPlotData.mitr_data[self.MagDirIndx[ii]])
                    curr_data_y.append(
                        self.ASDPlotData.mag_data[self.MagDirIndx[ii]])
                    curr_labels.append(
                        self.ASDPlotData.mag_labels[self.MagDirIndx[ii]])
                if len(self.MagDirIndx) > 0:
                    self.ASDPlots2D.LinePlot(self.Plotting_ax, curr_data_x, curr_data_y,
                                             curr_labels, self.ASDPlotData.mag_axes)
                else:
                    print('Select at least one direction to plot')
            else:
                self.ave_Error_Window = ASDInputWindows.Error_Window()
                self.ave_Error_Window.FunMsg.setText(
                    "I'm sorry, Dave. I'm afraid I can't do that.")
                self.ave_Error_Window.ErrorMsg.setText(
                    "Error: No 'averages.*.out' file.")
                self.ave_Error_Window.show()
                print("No 'averages.*.out' file.")
        # -----------------------------------------------------------------------
        # Plotting the total energy
        # -----------------------------------------------------------------------
        if self.plotting_mode == 'energy':
            self.EneOpts.setEnabled(True)
            if self.ASDPlotData.ene_file_present:
                curr_data_x = []
                curr_data_y = []
                curr_labels = []
                for ii in range(len(self.EneIndx)):
                    curr_data_x.append(
                        self.ASDPlotData.eitr_data[self.EneIndx[ii]])
                    curr_data_y.append(
                        self.ASDPlotData.ene_data[self.EneIndx[ii]])
                    curr_labels.append(
                        self.ASDPlotData.ene_labels[self.EneIndx[ii]])
                if len(self.EneIndx) > 0:
                    self.ASDPlots2D.LinePlot(self.Plotting_ax, curr_data_x, curr_data_y,
                                             curr_labels, self.ASDPlotData.ene_axes)
                else:
                    print('Select at least one energy component to plot')
            else:
                self.ene_Error_Window = ASDInputWindows.Error_Window()
                self.ene_Error_Window.FunMsg.setText(
                    "I'm sorry, Dave. I'm afraid I can't do that.")
                self.ene_Error_Window.ErrorMsg.setText(
                    "Error: No 'totenergy.*.out' file.")
                self.ene_Error_Window.show()
                print("No 'totenergy.*.out' file.")
        # -----------------------------------------------------------------------
        # Plotting the single magnetic moment trajectories
        # -----------------------------------------------------------------------
        if self.plotting_mode == 'trajectory':
            if self.ASDPlotData.trajectory_file_present:
                self.ASDPlots2D.TrajPlot(self.Plotting_ax3D, self.ASDPlotData.traj_data_x,
                                         self.ASDPlotData.traj_data_y, self.ASDPlotData.traj_data_z, self.ASDPlotData.traj_label)
            else:
                self.traj_Error_Window = ASDInputWindows.Error_Window()
                self.traj_Error_Window.FunMsg.setText(
                    "I'm sorry, Dave. I'm afraid I can't do that.")
                self.traj_Error_Window.ErrorMsg.setText(
                    "Error: No 'trajectory.*.out' file.")
                self.traj_Error_Window.show()
                print("No 'trajectory.*.out' file.")
        self.Plotting_Figure.canvas.draw()
        self.Plotting_Figure.canvas.flush_events()
        return
    ############################################################################
    # @brief Function to save the current figure to file
    # @author Jonathan Chico
    ############################################################################

    def SaveFig(self):
        from PyQt6.QtWidgets import QFileDialog
        fig_name, _ = QFileDialog.getSaveFileName(self, 'Save File')
        if len(self.InpFigDPI.text()) > 0:
            dpi = int(self.InpFigDPI.text())
        else:
            dpi = 800
        if self.plotting_mode != 'trajectory':
            fig_plot = self.Plotting_canvas
            fig_plot.print_figure(fig_name, dpi=dpi)
        else:
            fig_plot = self.Plotting_canvas3D
            fig_plot.print_figure(fig_name, dpi=dpi)
        return
    ############################################################################
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
    ############################################################################

    def AddActors(self):
        """Wrapper function that takes care of adding the necessary actors and the
        options for the different types of visualizations. It controls the visualization of:
            * Restartfiles
            * Momentsfiles
            * Energy
            * Exchange neighbours
            * DM neighbours

        Author
        ----------
        Jonathan Chico

        """
        from PyQt6.QtWidgets import QLabel
        from ASD_GUI.VTK_Viz import ASDVTKEneActors
        from ASD_GUI.VTK_Viz import ASDVTKMomActors
        from ASD_GUI.VTK_Viz import ASDVTKNeighActors
        try:
            self.ASDGenActors.scalar_bar_widget
        except:
            pass
        else:
            self.ASDGenActors.reset_GenActors()
        self.ren.RemoveAllViewProps()
        self.ren.ResetCamera()
        self.InitUI()
        # -----------------------------------------------------------------------
        # This takes care of setting up the options for the visualization of the
        # magnetic moments obtained from the restart file
        # -----------------------------------------------------------------------
        if self.sender() == self.actionMagnetization:
            # -------------------------------------------------------------------
            # Call the Moments class
            # -------------------------------------------------------------------
            self.MomActors = ASDVTKMomActors.ASDMomActors()
            self.viz_type = 'M'
            self.mode = 1
            self.current_time = 0
            self.MagMainGroup.setEnabled(True)
            self.VizToolBox.setCurrentIndex(0)
            self.menuMagnetisation_Opts.setEnabled(True)
            self.actionDisplayMagDens.setEnabled(True)
            self.actionX_ProjMagDens.setEnabled(True)
            self.actionY_ProjMagDens.setEnabled(True)
            self.actionZ_ProjMagDens.setEnabled(True)
            self.PlayButton.setEnabled(True)
            self.PauseButton.setEnabled(True)
            self.nextButton.setEnabled(True)
            self.previousButton.setEnabled(True)
            # -------------------------------------------------------------------
            # Add the data structures with regards to reading the data
            # -------------------------------------------------------------------
            self.ASDdata.ReadingWrapper(mode=self.mode, viz_type=self.viz_type,
                                        file_names=self.file_names, window=self)
            if not self.ASDdata.error_trap:
                self.MomActors.Add_MomActors(ren=self.ren, renWin=self.renWin,
                                             iren=self.iren, ASDdata=self.ASDdata, window=self)
                self.ASDVizOpt.update_dock_info(
                    current_Actors=self.MomActors, Window=self)
                # ---------------------------------------------------------------
                # Setup several global variables
                # ---------------------------------------------------------------
                self.ASDVizOpt.lut = self.MomActors.lut
                # ---------------------------------------------------------------
                # Add the general widgets such as the scalar bar and the axes
                # ---------------------------------------------------------------
                self.ASDGenActors.Add_GenActors(iren=self.iren, renWin=self.renWin,
                                                method=self.MomActors.MagDensMethod, lut=self.ASDVizOpt.lut,
                                                ren=self.ren, window=self, current_Actors=self.MomActors,
                                                flag2D=self.ASDdata.flag2D)
                # ---------------------------------------------------------------
                # Update the UI
                # ---------------------------------------------------------------
                if self.ASDdata.cluster_flag:
                    self.ClusBox.setVisible(True)
                    self.ClusBox.setChecked(True)
                    self.ASDGenActors.Add_ClusterActors(ASDdata=self.ASDdata, iren=self.iren,
                                                        renWin=self.renWin, ren=self.ren)
                if self.ASDdata.kmc_flag:
                    self.KMCCheck.setVisible(True)
                # ---------------------------------------------------------------
                # Print the visualization message
                # ---------------------------------------------------------------
                print('Visualization of magnetic moments mode chosen')
                self.current_time = self.current_time+1
        # -----------------------------------------------------------------------
        # This takes care of setting up the options for the Neighbour visualization
        # -----------------------------------------------------------------------
        if self.sender() == self.actionNeighbours:
            # -------------------------------------------------------------------
            # Call the Neighbour class
            # -------------------------------------------------------------------
            self.NeighActors = ASDVTKNeighActors.ASDNeighActors()
            self.mode = 1
            self.viz_type = 'N'
            self.NeighMainBox.setEnabled(True)
            self.VizToolBox.setCurrentIndex(2)
            self.PlayButton.setEnabled(False)
            self.PauseButton.setEnabled(False)
            self.nextButton.setEnabled(False)
            self.previousButton.setEnabled(False)
            # -------------------------------------------------------------------
            # Add the data structures with regards to reading the data
            # -------------------------------------------------------------------
            self.ASDdata.ReadingWrapper(mode=self.mode, viz_type=self.viz_type,
                                        file_names=self.file_names, window=self)
            if not self.ASDdata.error_trap:
                self.NeighActors.Add_NeighActors(ren=self.ren, renWin=self.renWin,
                                                 iren=self.iren, ASDdata=self.ASDdata, mode=self.mode)
                # ---------------------------------------------------------------
                # Set several global variables
                # ---------------------------------------------------------------
                self.ASDVizOpt.lut = self.NeighActors.lut
                # ---------------------------------------------------------------
                # Add the general widgets such as the scalar bar and the axes
                # ---------------------------------------------------------------
                self.ASDGenActors.Add_GenActors(iren=self.iren, renWin=self.renWin,
                                                method=self.NeighActors.NeighGlyph3D, lut=self.ASDVizOpt.lut,
                                                ren=self.ren, window=self, current_Actors=self.NeighActors, flag2D=True)
                # ---------------------------------------------------------------
                # Update the labels for the neighbour mode
                # ---------------------------------------------------------------
                self.NeighSelectSlider.setMaximum(self.NeighActors.SLMax)
                self.NeighSelectSlider.setMinimum(1)
                self.NeighNumberLabel.setText(
                    'Number of neighbours = {: 4d}'.format(self.NeighActors.NumNeigh))
                self.NeighValidator.setRange(1, self.ASDdata.nrAtoms)
                self.NeighSelectLineEdit.setValidator(self.NeighValidator)
                self.ASDVizOpt.update_dock_info(
                    current_Actors=self.NeighActors, Window=self)
                # ---------------------------------------------------------------
                # Update the UI
                # ---------------------------------------------------------------
                self.NeighTypesLabels = dict()
                for ii in range(0, self.ASDdata.num_types_total):
                    name = 'label_neigh_{}'.format(ii)
                    label = QLabel()
                    label.setObjectName(name)
                    label.setText(
                        'Num. Neighbours Type {: 4d} = {: 4d}'.format(ii+1, 0))
                    self.NeighInfoLayout.addWidget(label)
                    self.NeighTypesLabels[name] = label
                for ii in range(0, self.ASDdata.num_types):
                    name = 'label_neigh_{}'.format(
                        int(self.ASDdata.types[ii]-1))
                    self.NeighTypesLabels[name].setText('Num. Neighbours Type {: 4d} = {: 4d}'.format(
                        ii+1, self.ASDdata.types_counters[ii]))
                # ---------------------------------------------------------------
                # Visualize the embedded cluster into the system
                # ---------------------------------------------------------------
                if self.ASDdata.cluster_flag:
                    self.ClusBox.setVisible(True)
                    self.ClusBox.setChecked(True)
                    self.ASDGenActors.Add_ClusterActors(ASDdata=self.ASDdata, iren=self.iren,
                                                        renWin=self.renWin, ren=self.ren)
                # ---------------------------------------------------------------
                # Print the visualization message
                # ---------------------------------------------------------------
                print('Visualization of the neighbour map mode chosen')
                print('Viewing the struct file')
        # -----------------------------------------------------------------------
        # This takes care of setting up the options for the DM Neighbour visualization
        # -----------------------------------------------------------------------
        if self.sender() == self.actionDM_Neigh:
            # -------------------------------------------------------------------
            # Call the Neighbour class
            # -------------------------------------------------------------------
            self.NeighActors = ASDVTKNeighActors.ASDNeighActors()
            self.mode = 2
            self.viz_type = 'N'
            self.NeighMainBox.setEnabled(True)
            self.VizToolBox.setCurrentIndex(2)
            self.PlayButton.setEnabled(False)
            self.PauseButton.setEnabled(False)
            self.nextButton.setEnabled(False)
            self.previousButton.setEnabled(False)
            # -------------------------------------------------------------------
            # Add the data structures with regards to reading the data
            # -------------------------------------------------------------------
            self.ASDdata.ReadingWrapper(mode=self.mode, viz_type=self.viz_type,
                                        file_names=self.file_names, window=self)
            if not self.ASDdata.error_trap:
                self.NeighActors.Add_NeighActors(ren=self.ren, renWin=self.renWin,
                                                 iren=self.iren, ASDdata=self.ASDdata, mode=self.mode)
                # ---------------------------------------------------------------
                # Set several global variables
                # ---------------------------------------------------------------
                self.ASDVizOpt.lut = self.NeighActors.lut
                # ---------------------------------------------------------------
                # Add the general widgets such as the scalar bar and the axes
                # ---------------------------------------------------------------
                self.ASDGenActors.Add_GenActors(iren=self.iren, renWin=self.renWin,
                                                method=self.NeighActors.NeighGlyph3D, lut=self.ASDVizOpt.lut,
                                                ren=self.ren, window=self, current_Actors=self.NeighActors, flag2D=True)
                # ---------------------------------------------------------------
                # Update the labels for the neighbour mode
                # ---------------------------------------------------------------
                self.NeighSelectSlider.setMaximum(self.NeighActors.SLMax)
                self.NeighSelectSlider.setMinimum(1)
                self.NeighNumberLabel.setText(
                    'Number of neighbours = {: 4d}'.format(self.NeighActors.NumNeigh))
                self.ASDVizOpt.update_dock_info(
                    current_Actors=self.NeighActors, Window=self)
                # ---------------------------------------------------------------
                # Update the UI
                # ---------------------------------------------------------------
                self.NeighTypesLabels = dict()
                for ii in range(0, self.ASDdata.num_types_total):
                    name = 'label_neigh_{}'.format(ii)
                    label = QLabel()
                    label.setObjectName(name)
                    label.setText(
                        'Num. Neighbours Type {: 4d} = {: 4d}'.format(ii+1, 0))
                    self.NeighInfoLayout.addWidget(label)
                    self.NeighTypesLabels[name] = label
                for ii in range(0, self.ASDdata.num_types):
                    name = 'label_neigh_{}'.format(
                        int(self.ASDdata.types[ii]-1))
                    self.NeighTypesLabels[name].setText('Num. Neighbours Type {: 4d} = {: 4d}'.format(
                        ii+1, self.ASDdata.types_counters[ii]))
                    # -----------------------------------------------------------
                    # Visualize the embedded cluster into the system
                    # -----------------------------------------------------------
                if self.ASDdata.cluster_flag:
                    self.ClusBox.setVisible(True)
                    self.ClusBox.setChecked(True)
                    self.ASDGenActors.Add_ClusterActors(ASDdata=self.ASDdata, iren=self.iren,
                                                        renWin=self.renWin, ren=self.ren)
                # ---------------------------------------------------------------
                # Print the visualization message
                # ---------------------------------------------------------------
                print('Visualization of the neighbour map mode chosen')
                print('Viewing the struct file')
        # -----------------------------------------------------------------------
        # This takes care of setting up the options for the Energy visualization
        # -----------------------------------------------------------------------
        if self.sender() == self.actionEnergy:
            self.EneActors = ASDVTKEneActors.ASDEneActors()
            self.viz_type = 'E'
            self.mode = 1
            self.current_time = 0
            self.VizToolBox.setCurrentIndex(1)
            self.EneMainBox.setEnabled(True)
            self.PlayButton.setEnabled(True)
            self.PauseButton.setEnabled(True)
            self.nextButton.setEnabled(True)
            self.previousButton.setEnabled(True)
            # -------------------------------------------------------------------
            # Add the data structures with regards to reading the data
            # -------------------------------------------------------------------
            self.ASDdata.ReadingWrapper(mode=self.mode, viz_type=self.viz_type,
                                        file_names=self.file_names, window=self)
            if not self.ASDdata.error_trap:
                self.EneActors.Add_EneActors(ren=self.ren, renWin=self.renWin,
                                             iren=self.iren, ASDdata=self.ASDdata)
                self.ASDVizOpt.update_dock_info(
                    current_Actors=self.EneActors, Window=self)
                # ---------------------------------------------------------------
                # Setup several global variables
                # ---------------------------------------------------------------
                self.ASDVizOpt.lut = self.EneActors.lut
                # ---------------------------------------------------------------
                # Add the general widgets such as the scalar bar and the axes
                # ---------------------------------------------------------------
                self.ASDGenActors.Add_GenActors(iren=self.iren, renWin=self.renWin,
                                                method=self.EneActors.EneDensMethod, lut=self.ASDVizOpt.lut,
                                                ren=self.ren, window=self, current_Actors=self.EneActors,
                                                flag2D=self.ASDdata.flag2D)
                # ---------------------------------------------------------------
                # Update the UI
                # ---------------------------------------------------------------
                if self.ASDdata.cluster_flag:
                    self.ClusBox.setVisible(True)
                    self.ClusBox.setChecked(True)
                    self.ASDGenActors.Add_ClusterActors(ASDdata=self.ASDdata, iren=self.iren,
                                                        renWin=self.renWin, ren=self.ren)
                # ---------------------------------------------------------------
                # Print the visualization message
                # ---------------------------------------------------------------
                print('Visualization of the energy mode chosen')
                print('Viewing the localenergy file')
        return
    ############################################################################
    # @brief Enable rgb-values for single color
    # @author Anders Bergman
    ############################################################################

    def toggle_singlecolor(self, check):
        if check:
            self.RGBRedColorSlider.setEnabled(True)
            self.RGBGreenColorSlider.setEnabled(True)
            self.RGBBlueColorSlider.setEnabled(True)
        else:
            self.RGBRedColorSlider.setEnabled(False)
            self.RGBGreenColorSlider.setEnabled(False)
            self.RGBBlueColorSlider.setEnabled(False)

        return
    ############################################################################
    # @brief Toggle grayscale background on/off
    # @author Anders Bergman
    ############################################################################

    def toggle_bwSinglecolor(self, check):

        self.bwSinglecolor = check
        rgb = [self.RGBRedColorSlider.value(),
               self.RGBGreenColorSlider.value(),
               self.RGBBlueColorSlider.value()]
        bw = int(sum(rgb)/3)

        if check:
            self.RGBRedColorSlider.setValue(bw)
            self.RGBGreenColorSlider.setValue(bw)
            self.RGBRedColorSlider.setValue(bw)

        return
    ############################################################################
    # @brief Toggle depth of field focus
    # @author Anders Bergman
    ############################################################################

    def toggle_focus(self, check):
        self.ASDVizOpt.toggle_Focus(
            check=check, ren=self.ren, renWin=self.renWin)

    ############################################################################
    # @brief Toggle focal disk
    # @author Anders Bergman
    ############################################################################
    def FocalDisk_control(self, value):
        self.ASDVizOpt.setFocalDisk(
            value=value, ren=self.ren, renWin=self.renWin)

    ############################################################################
    # @brief Toggle depth of field focus
    # @author Anders Bergman
    ############################################################################

    def toggle_autofocus(self, check):
        self.ASDVizOpt.toggle_autoFocus(check=check, renWin=self.renWin)

    ############################################################################
    # @brief Toggle grayscale background on/off
    # @author Anders Bergman
    ############################################################################

    def toggle_bwBackground(self, check):

        self.bwBackground = check
        rgb = [self.RGBRedBackgroundSlider.value(),
               self.RGBGreenBackgroundSlider.value(),
               self.RGBBlueBackgroundSlider.value()]
        bw = int(sum(rgb)/3)

        if check:
            self.RGBRedBackgroundSlider.setValue(bw)
            self.RGBGreenBackgroundSlider.setValue(bw)
            self.RGBRedBackgroundSlider.setValue(bw)

        return
    ############################################################################
    # @brief Update rgb-values for single color coloring
    # @author Anders Bergman
    ############################################################################

    def set_singlecolor(self, value):

        if self.bwSinglecolor:
            self.RGBRedColorSlider.setValue(value)
            self.RGBGreenColorSlider.setValue(value)
            self.RGBBlueColorSlider.setValue(value)

        rgb = [self.RGBRedColorSlider.value(),
               self.RGBGreenColorSlider.value(),
               self.RGBBlueColorSlider.value()]

        self.ASDVizOpt.set_RGBcolor(window=self, rgb=rgb,
                                    flag2D=self.ASDdata.flag2D, viz_type=self.viz_type, renWin=self.renWin)

        return
    ############################################################################
    # @brief Update rgb-values for the background
    # @author Anders Bergman
    ############################################################################

    def set_background(self, value):

        if self.bwBackground:
            self.RGBRedBackgroundSlider.setValue(value)
            self.RGBGreenBackgroundSlider.setValue(value)
            self.RGBBlueBackgroundSlider.setValue(value)

        rgb = [self.RGBRedBackgroundSlider.value(),
               self.RGBGreenBackgroundSlider.value(),
               self.RGBBlueBackgroundSlider.value()]

        self.ASDVizOpt.set_RGBbackground(
            rgb=rgb, ren=self.ren, renWin=self.renWin)

        return
    ############################################################################
    # @brief Set the lookup table for the actors
    # @details Set the lookup table for the actors, it also allows for the change
    # of the scale type for the plotting, either linear or logarithmic scale.
    # @author Jonathan Chico
    ############################################################################

    def set_lut_db(self, mapnum):
        from vtkmodules.vtkCommonColor import vtkColorSeries, vtkNamedColors

        colorSeries = vtkColorSeries()

        if mapnum <= 3:
            self.ASDVizOpt.set_colormap_db(window=self, mapnum=mapnum,
                                           flag2D=self.ASDdata.flag2D, viz_type=self.viz_type, renWin=self.renWin)
        elif mapnum == 4:  # Spectrum
            colorSeries.SetColorScheme(vtkColorSeries.SPECTRUM)
        elif mapnum == 5:  # Warm
            colorSeries.SetColorScheme(vtkColorSeries.WARM)
        elif mapnum == 6:  # Cool
            colorSeries.SetColorScheme(vtkColorSeries.COOL)
        elif mapnum == 7:  # Blues
            colorSeries.SetColorScheme(vtkColorSeries.BLUES)
        elif mapnum == 8:  # Wildflower
            colorSeries.SetColorScheme(vtkColorSeries.WILD_FLOWER)
        elif mapnum == 9:  # Citrus
            colorSeries.SetColorScheme(vtkColorSeries.CITRUS)
        elif mapnum == 10:  # BREWER_DIVERGING_PURPLE_ORANGE_11
            colorSeries.SetColorScheme(
                vtkColorSeries.BREWER_DIVERGING_PURPLE_ORANGE_11)
        elif mapnum == 11:  # Citrus
            colorSeries.SetColorScheme(
                vtkColorSeries.BREWER_DIVERGING_BROWN_BLUE_GREEN_11)
        elif mapnum == 12:  # Citrus
            colorSeries.SetColorScheme(
                vtkColorSeries.BREWER_SEQUENTIAL_BLUE_GREEN_9)
        elif mapnum == 13:  # Citrus
            colorSeries.SetColorScheme(
                vtkColorSeries.BREWER_SEQUENTIAL_YELLOW_ORANGE_BROWN_9)
        elif mapnum == 14:  # Citrus
            colorSeries.SetColorScheme(
                vtkColorSeries.BREWER_SEQUENTIAL_BLUE_PURPLE_9)
        elif mapnum == 15:  # Citrus
            colorSeries.SetColorScheme(
                vtkColorSeries.BREWER_DIVERGING_SPECTRAL_11)
        elif mapnum == 16:  # Citrus
            colorSeries.SetColorScheme(
                vtkColorSeries.BREWER_QUALITATIVE_ACCENT)
        elif mapnum == 17:  # Citrus
            colorSeries.SetColorScheme(vtkColorSeries.BREWER_QUALITATIVE_DARK2)
        elif mapnum == 18:  # Citrus
            colorSeries.SetColorScheme(
                vtkColorSeries.BREWER_QUALITATIVE_PASTEL1)
        elif mapnum == 19:  # Citrus
            colorSeries.SetColorScheme(
                vtkColorSeries.BREWER_QUALITATIVE_PASTEL2)
        elif mapnum == 20:  # Citrus
            colorSeries.SetColorScheme(vtkColorSeries.BREWER_QUALITATIVE_SET1)
        elif mapnum == 21:  # Citrus
            colorSeries.SetColorScheme(vtkColorSeries.BREWER_QUALITATIVE_SET2)
        elif mapnum == 22:  # Citrus
            colorSeries.SetColorScheme(vtkColorSeries.BREWER_QUALITATIVE_SET3)

        if mapnum > 3:
            colorSeries.BuildLookupTable(
                self.ASDVizOpt.lut, vtkColorSeries.ORDINAL)
        self.ASDVizOpt.lut.Build()

        self.renWin.Render()
        return
    ############################################################################
    # @brief Set the lookup table for the actors
    # @details Set the lookup table for the actors, it also allows for the change
    # of the scale type for the plotting, either linear or logarithmic scale.
    # @author Jonathan Chico
    ############################################################################

    def set_lut(self):
        # self.ASDVizOpt.set_colormap(window=self,flag2D=self.ASDdata.flag2D,\
        # viz_type=self.viz_type,renWin=self.renWin)
        if self.sender() == self.LinearScale and self.LinearScale.isChecked():
            self.ASDVizOpt.lut.SetScaleToLinear()
        if self.sender() == self.LogScale and self.LogScale.isChecked():
            self.ASDVizOpt.lut.SetScaleToLog10()
        return
    ############################################################################
    # @brief Set the projection of the vectors
    # @details Set the projection of the vectors and the magnetization continuum
    # allowing one to set independent projections of the continuum visualization
    # and the spins.
    # @author Jonathan Chico
    ############################################################################

    def set_projection(self):
        """Set the projection of the vectors and the magnetization continuum
        allowing one to set independent projections of the continuum visualization
        and the spins.

        Author
        ----------
        Jonathan Chico
        """
        if self.sender() == self.DensX and self.DensX.isChecked():
            self.ASDVizOpt.set_projection(type='density', axis=0)
        if self.sender() == self.DensY and self.DensY.isChecked():
            self.ASDVizOpt.set_projection(type='density', axis=1)
        if self.sender() == self.DensZ and self.DensZ.isChecked():
            self.ASDVizOpt.set_projection(type='density', axis=2)
        if self.sender() == self.SpinX and self.SpinX.isChecked():
            self.ASDVizOpt.set_projection(type='spins', axis=0)
        if self.sender() == self.SpinY and self.SpinY.isChecked():
            self.ASDVizOpt.set_projection(type='spins', axis=1)
        if self.sender() == self.SpinZ and self.SpinZ.isChecked():
            self.ASDVizOpt.set_projection(type='spins', axis=2)
        self.renWin.Render()
        return
    ############################################################################
    # Display the different energy contributions
    ############################################################################

    def set_energy_proj(self):
        if self.sender() == self.TotEneButton and self.TotEneButton.isChecked():
            self.EneActors.src.GetPointData().SetScalars(
                self.ASDdata.energies[0])
            self.EneActors.EneMapper.SetScalarRange(
                self.EneActors.src.GetScalarRange())
            # self.EneActors.EneDensMap.SetScalarRange(self.EneActors.src.GetScalarRange())
            self.ASDGenActors.clipperMapper.SetScalarRange(
                self.EneActors.src.GetScalarRange())
        if self.sender() == self.ExcEneButton and self.ExcEneButton.isChecked():
            self.EneActors.src.GetPointData().SetScalars(
                self.ASDdata.energies[1])
            self.EneActors.EneMapper.SetScalarRange(
                self.EneActors.src.GetScalarRange())
            # self.EneActors.EneDensMap.SetScalarRange(self.EneActors.src.GetScalarRange())
            self.ASDGenActors.clipperMapper.SetScalarRange(
                self.EneActors.src.GetScalarRange())
        if self.sender() == self.DMEneButton and self.DMEneButton.isChecked():
            self.EneActors.src.GetPointData().SetScalars(
                self.ASDdata.energies[2])
            self.EneActors.EneMapper.SetScalarRange(
                self.EneActors.src.GetScalarRange())
            # self.EneActors.EneDensMap.SetScalarRange(self.EneActors.src.GetScalarRange())
            self.ASDGenActors.clipperMapper.SetScalarRange(
                self.EneActors.src.GetScalarRange())
        if self.sender() == self.AniEneButton and self.AniEneButton.isChecked():
            self.EneActors.src.GetPointData().SetScalars(
                self.ASDdata.energies[3])
            self.EneActors.EneMapper.SetScalarRange(
                self.EneActors.src.GetScalarRange())
            # self.EneActors.EneDensMap.SetScalarRange(self.EneActors.src.GetScalarRange())
            self.ASDGenActors.clipperMapper.SetScalarRange(
                self.EneActors.src.GetScalarRange())
        if self.sender() == self.BqEneButton and self.BqEneButton.isChecked():
            self.EneActors.src.GetPointData().SetScalars(
                self.ASDdata.energies[4])
            self.EneActors.EneMapper.SetScalarRange(
                self.EneActors.src.GetScalarRange())
            # self.EneActors.EneDensMap.SetScalarRange(self.EneActors.src.GetScalarRange())
            self.ASDGenActors.clipperMapper.SetScalarRange(
                self.EneActors.src.GetScalarRange())
        if self.sender() == self.BqDMEneButton and self.BqDMEneButton.isChecked():
            self.EneActors.src.GetPointData().SetScalars(
                self.ASDdata.energies[5])
            self.EneActors.EneMapper.SetScalarRange(
                self.EneActors.src.GetScalarRange())
            # self.EneActors.EneDensMap.SetScalarRange(self.EneActors.src.GetScalarRange())
            self.ASDGenActors.clipperMapper.SetScalarRange(
                self.EneActors.src.GetScalarRange())
        if self.sender() == self.PdEneButton and self.PdEneButton.isChecked():
            self.EneActors.src.GetPointData().SetScalars(
                self.ASDdata.energies[6])
            self.EneActors.EneMapper.SetScalarRange(
                self.EneActors.src.GetScalarRange())
            # self.EneActors.EneDensMap.SetScalarRange(self.EneActors.src.GetScalarRange())
            self.ASDGenActors.clipperMapper.SetScalarRange(
                self.EneActors.src.GetScalarRange())
        if self.sender() == self.BextEneButton and self.BextEneButton.isChecked():
            self.EneActors.src.GetPointData().SetScalars(
                self.ASDdata.energies[7])
            self.EneActors.EneMapper.SetScalarRange(
                self.EneActors.src.GetScalarRange())
            # self.EneActors.EneDensMap.SetScalarRange(self.EneActors.src.GetScalarRange())
            self.ASDGenActors.clipperMapper.SetScalarRange(
                self.EneActors.src.GetScalarRange())
        if self.sender() == self.DipEneButton and self.DipEneButton.isChecked():
            self.EneActors.src.GetPointData().SetScalars(
                self.ASDdata.energies[8])
            self.EneActors.EneMapper.SetScalarRange(
                self.EneActors.src.GetScalarRange())
            # self.EneActors.EneDensMap.SetScalarRange(self.EneActors.src.GetScalarRange())
            self.ASDGenActors.clipperMapper.SetScalarRange(
                self.EneActors.src.GetScalarRange())
        if self.sender() == self.ChirEneButton and self.ChirEneButton.isChecked():
            self.EneActors.src.GetPointData().SetScalars(
                self.ASDdata.energies[9])
            self.EneActors.EneMapper.SetScalarRange(
                self.EneActors.src.GetScalarRange())
            # self.EneActors.EneDensMap.SetScalarRange(self.EneActors.src.GetScalarRange())
            self.ASDGenActors.clipperMapper.SetScalarRange(
                self.EneActors.src.GetScalarRange())
        self.EneActors.EneDensMap.Modified()
        self.EneActors.EneMapper.Update()
        self.renWin.Render()
        return
    ############################################################################
    # Function to change the type of glyphs that display the individual magnetic
    # moments
    ############################################################################

    def ChangeGlyphs(self):
        """Function to change the type of glyphs that display the individual magnetic
        moments

        Author
        ----------
        Jonathan Chico
        """
        if self.sender() == self.SpinBarButton:
            if self.SpinBarButton.isChecked():
                self.ASDVizOpt.ChangeSpinGlyph(
                    renWin=self.renWin, keyword='Bars')
                self.SpinCenterCheck.setEnabled(False)
        if self.sender() == self.SpinCubeButton:
            if self.SpinCubeButton.isChecked():
                self.ASDVizOpt.ChangeSpinGlyph(
                    renWin=self.renWin, keyword='Cubes')
                self.SpinCenterCheck.setEnabled(False)
        if self.sender() == self.SpinSphereButton:
            if self.SpinSphereButton.isChecked():
                self.ASDVizOpt.ChangeSpinGlyph(
                    renWin=self.renWin, keyword='Spheres')
                self.SpinCenterCheck.setEnabled(False)
        if self.sender() == self.SpinArrowButton:
            if self.SpinArrowButton.isChecked():
                self.ASDVizOpt.ChangeSpinGlyph(
                    renWin=self.renWin, keyword='Arrows')
                self.SpinCenterCheck.setChecked(False)
                self.SpinCenterCheck.setEnabled(True)
        if self.sender() == self.SpinConeButton:
            if self.SpinConeButton.isChecked():
                self.ASDVizOpt.ChangeSpinGlyph(
                    renWin=self.renWin, keyword='Cones')
                self.SpinCenterCheck.setEnabled(False)
        if self.sender() == self.SpinCenterCheck:
            if self.SpinCenterCheck.isChecked():
                self.ASDVizOpt.ChangeSpinGlyph(
                    renWin=self.renWin, keyword='CenterOn')
            else:
                self.ASDVizOpt.ChangeSpinGlyph(
                    renWin=self.renWin, keyword='CenterOff')
        self.renWin.Render()
        return
    ############################################################################
    # Function to change the shading of the glyphs that display the individual
    # magnetic moments
    ############################################################################

    def ChangeShading(self):
        """Function to change the type of glyphs that display the individual magnetic
        moments

        Author
        ----------
        Anders Bergman, Jonathan Chico
        """
        if self.sender() == self.FlatShadeButton:
            if self.FlatShadeButton.isChecked():
                self.ASDVizOpt.ChangeSpinShade(
                    renWin=self.renWin, keyword='Flat')

        if self.sender() == self.GouraudShadeButton:
            if self.GouraudShadeButton.isChecked():
                self.ASDVizOpt.ChangeSpinShade(
                    renWin=self.renWin, keyword='Gouraud')

        if self.sender() == self.PBRShadeButton:
            if self.PBRShadeButton.isChecked():
                self.ASDVizOpt.ChangeSpinShade(
                    renWin=self.renWin, keyword='PBR')

        if self.sender() == self.PhongShadeButton:
            if self.PhongShadeButton.isChecked():
                self.ASDVizOpt.ChangeSpinShade(
                    renWin=self.renWin, keyword='Phong')

        self.renWin.Render()
        return
    ############################################################################
    # Update the neighbours
    ############################################################################

    def NeighbourControl(self):
        self.NeighActors.UpdateNeighbour(window=self, ASDdata=self.ASDdata,
                                         ASDGenActors=self.ASDGenActors, renWin=self.renWin, mode=self.mode)
        return
    ############################################################################
    # Wrapper function to handle the camera functions
    ############################################################################

    def camera_handler(self):
        # -----------------------------------------------------------------------
        # Reset the camera to the original position
        # -----------------------------------------------------------------------
        if self.sender() == self.CamResetButton:
            if self.viz_type == 'M':
                self.ASDVizOpt.reset_camera(ren=self.ren, renWin=self.renWin,
                                            current_Actors=self.MomActors)
            elif self.viz_type == 'N':
                self.ASDVizOpt.reset_camera(ren=self.ren, renWin=self.renWin,
                                            current_Actors=self.NeighActors)
            elif self.viz_type == 'E':
                self.ASDVizOpt.reset_camera(ren=self.ren, renWin=self.renWin,
                                            current_Actors=self.EneActors)
        # -----------------------------------------------------------------------
        # Controlling what is up in the camera
        # -----------------------------------------------------------------------
        if self.sender() == self.SetXView:
            self.ASDVizOpt.set_Camera_viewUp(
                ren=self.ren, renWin=self.renWin, dir=(1, 0, 0))
        if self.sender() == self.SetYView:
            self.ASDVizOpt.set_Camera_viewUp(
                ren=self.ren, renWin=self.renWin, dir=(0, 1, 0))
        if self.sender() == self.SetZView:
            self.ASDVizOpt.set_Camera_viewUp(
                ren=self.ren, renWin=self.renWin, dir=(0, 0, 1))
        if self.sender() == self.SetCamButton:
            self.ASDVizOpt.Update_Camera(
                Window=self, ren=self.ren, renWin=self.renWin)
        # -----------------------------------------------------------------------
        # Controlling the parallel scale
        # -----------------------------------------------------------------------
        if self.sender() == self.ParallelScaleLineEdit:
            line = True
            slider = False
            self.ASDVizOpt.ChangeParallelProj(ren=self.ren, renWin=self.renWin,
                                              line=line, slider=slider, MainWindow=self)
        if self.sender() == self.ParallelScaleSlider:
            line = False
            slider = True
            self.ASDVizOpt.ChangeParallelProj(ren=self.ren, renWin=self.renWin,
                                              line=line, slider=slider, MainWindow=self)
        if self.sender() == self.ParallelProjectBox:
            self.ASDVizOpt.toggle_projections(renWin=self.renWin, window=self,
                                              ren=self.ren, checked=self.ParallelProjectBox.isChecked())
        return
    ############################################################################
    # Wrapper to handle the clipper actions
    ############################################################################

    def clipperHandler(self):
        if self.viz_type == 'M':
            current_Actors = self.MomActors
        if self.viz_type == 'N':
            current_Actors = self.NeighActors
        if self.viz_type == 'E':
            current_Actors = self.EneActors
        self.ASDGenActors.UpdateClipper(window=self, current_Actors=current_Actors,
                                        ASDVizOpt=self.ASDVizOpt, renWin=self.renWin, viz_type=self.viz_type)
        return
    ############################################################################
    # Function that calls the taking of a Snapshot of the current rendering window
    ############################################################################

    def Snapshot(self):
        self.ASDVizOpt.Screenshot(renWin=self.renWin, number_of_screenshots=self.number_of_screenshots,
                                  png_mode=self.actionSave_png.isChecked(), pov_mode=self.actionSave_pov.isChecked())
        self.number_of_screenshots = self.number_of_screenshots+1
        return
    ############################################################################
    # Function that calls for updating the glyph resolutions
    ############################################################################

    def Quality_control(self, value):
        self.ASDVizOpt.GlyphQualityUpdate(
            value=value, viz_type=self.viz_type, mode=self.mode, renWin=self.renWin)
    ############################################################################
    # Function that calls for toggling FXAA
    ############################################################################

    def FXAA_control(self, check):
        self.ASDVizOpt.toggle_FXAA(
            check=check, ren=self.ren, renWin=self.renWin)
    ############################################################################
    # Function that calls for toggling surface texture
    ############################################################################

    def Texture_control(self, check):
        self.ASDVizOpt.toggle_Texture(check=check, ren=self.ren, renWin=self.renWin,
                                      texfile=self.texturefile)
    ############################################################################
    # Function that calls for toggling ORM texture
    ############################################################################

    def ORMTexture_control(self, check):
        self.ASDVizOpt.toggle_ORMTexture(check=check, ren=self.ren, renWin=self.renWin,
                                         texfile=self.ORMtexturefile)
    ############################################################################
    # Function that calls for toggling ORM texture
    ############################################################################

    def NTexture_control(self, check):
        self.ASDVizOpt.toggle_NTexture(check=check, ren=self.ren, renWin=self.renWin,
                                       texfile=self.Ntexturefile)
    ############################################################################
    # Function that calls for toggling ORM texture
    ############################################################################

    def ETexture_control(self, check):
        self.ASDVizOpt.toggle_ETexture(check=check, ren=self.ren, renWin=self.renWin,
                                       texfile=self.Etexturefile)
    ############################################################################
    # Function that calls for toggling ORM texture
    ############################################################################

    def ATexture_control(self, check):
        self.ASDVizOpt.toggle_ATexture(check=check, ren=self.ren, renWin=self.renWin,
                                       texfile=self.Atexturefile)
    ############################################################################
    # Function that calls for toggling SSAO
    ############################################################################

    def SSAO_control(self, check):
        self.ASDVizOpt.toggle_SSAO(check=check, ren=self.ren)
    ############################################################################
    # Function that calls for toggling shadows
    ############################################################################
    # def Shadow_control(self, check):
    #    self.ASDVizOpt.toggle_Shadows(check=check,ren=self.ren, renWin=self.renWin)

    ############################################################################
    # Function that calls for toggling HDRI
    ############################################################################
    def HDRI_control(self, check):
        self.ASDVizOpt.toggle_HDRI(check=check, ren=self.ren, renWin=self.renWin,
                                   hdrifile=self.hdrifile)
        return
    ############################################################################
    # Function that calls for toggling skybox
    ############################################################################

    def SkyBox_control(self, check):
        self.ASDVizOpt.toggle_SkyBox(check=check, ren=self.ren, renWin=self.renWin,
                                     skyboxfile=self.hdrifile)
        return
    ############################################################################
    # Finding the file name for the HDR file
    ############################################################################

    def getHDRIFile(self):
        self.hdrifile = self.ASDVizOpt.getHDRIFileName(window=self)
        self.hdrifile_gotten = len(self.hdrifile) > 0
        if self.hdrifile_gotten:
            self.HDRICheck.setEnabled(True)
            self.SkyBoxCheck.setEnabled(True)
        return
    ############################################################################
    # Finding the file name for the texture image
    ############################################################################

    def getTextureFile(self):
        self.texturefile = self.ASDVizOpt.getTextureFileName(window=self)
        self.texturefile_gotten = len(self.texturefile) > 0
        if self.texturefile_gotten:
            self.TextureCheck.setEnabled(True)
        return
    ############################################################################
    # Finding the file name for the ORM texture image
    ############################################################################

    def getORMTextureFile(self):
        self.ORMtexturefile = self.ASDVizOpt.getTextureFileName(window=self)
        self.ORMtexturefile_gotten = len(self.ORMtexturefile) > 0
        if self.ORMtexturefile_gotten:
            self.ORMTextureCheck.setEnabled(True)
        return
    ############################################################################
    # Finding the file name for the normal texture image
    ############################################################################

    def getNTextureFile(self):
        self.Ntexturefile = self.ASDVizOpt.getTextureFileName(window=self)
        self.Ntexturefile_gotten = len(self.Ntexturefile) > 0
        if self.Ntexturefile_gotten:
            self.NTextureCheck.setEnabled(True)
        return
    ############################################################################
    # Finding the file name for the anisotropy texture image
    ############################################################################

    def getATextureFile(self):
        self.Atexturefile = self.ASDVizOpt.getTextureFileName(window=self)
        self.Atexturefile_gotten = len(self.Atexturefile) > 0
        if self.Atexturefile_gotten:
            self.ATextureCheck.setEnabled(True)
        return
    ############################################################################
    # Finding the file name for the emissive texture image
    ############################################################################

    def getETextureFile(self):
        self.Etexturefile = self.ASDVizOpt.getTextureFileName(window=self)
        self.Etexturefile_gotten = len(self.Etexturefile) > 0
        if self.Etexturefile_gotten:
            self.ETextureCheck.setEnabled(True)
        return
    ############################################################################
    # Function that calls for toggling specular scattering
    ############################################################################

    def RenSpecular_control(self, value):
        self.ASDVizOpt.RenSpecularUpdate(value=value, renWin=self.renWin)
    ############################################################################
    # Function that calls for toggling specular scattering
    ############################################################################

    def RenSpecularPower_control(self, value):
        self.ASDVizOpt.RenSpecularPowerUpdate(value=value, renWin=self.renWin)
    ############################################################################
    # Function that calls for toggling ambient scattering
    ############################################################################

    def RenAmbient_control(self, value):
        self.ASDVizOpt.RenAmbientUpdate(value=value, renWin=self.renWin)
    ############################################################################
    # Function that calls for toggling diffuse scattering
    ############################################################################

    def RenDiffuse_control(self, value):
        self.ASDVizOpt.RenDiffuseUpdate(value=value, renWin=self.renWin)
    ############################################################################
    # Function that calls for toggling PBR Emission value
    ############################################################################

    def PBREmission_control(self, value):
        self.ASDVizOpt.PBREmissionUpdate(
            value=value, ren=self.ren, renWin=self.renWin)
    ############################################################################
    # Function that calls for toggling PBR Occlusion value
    ############################################################################

    def PBROcclusion_control(self, value):
        self.ASDVizOpt.PBROcclusionUpdate(
            value=value, ren=self.ren, renWin=self.renWin)
    ############################################################################
    # Function that calls for toggling PBR Roughness value
    ############################################################################

    def PBRRoughness_control(self, value):
        self.ASDVizOpt.PBRRoughnessUpdate(value=value, renWin=self.renWin)
    ############################################################################
    # Function that calls for toggling PBR Roughness value
    ############################################################################

    def PBRMetallic_control(self, value):
        self.ASDVizOpt.PBRMetallicUpdate(value=value, renWin=self.renWin)

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
        # -----------------------------------------------------------------------
        # Play the movie
        # -----------------------------------------------------------------------
        if self.sender() == self.PlayButton:
            if self.PlayButton.isChecked():
                self.iren.AddObserver('TimerEvent', self.Playback)
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
            if self.current_time < self.ASDdata.number_time_steps-1:
                self.current_time += 1
                self.UpdateImage()
        return
    ############################################################################
    # @brief Function to control the playback of the animation, whilst taking a snapshot
    # and updating the necessary data structures.
    # @details Function to control the playback of the animation, whilst taking a snapshot
    # and updating the necessary data structures. This function updates the magnetic
    # moments, the site dependent energy, as well as the timers to ensure that the
    # visualization finishes with the last image.
    # @author Jonathan Chico
    ############################################################################

    def Playback(self, event, obj):
        if self.viz_type == 'M':
            # -------------------------------------------------------------------
            # If the current time is the final time of the measurement destroy the
            # timer
            # -------------------------------------------------------------------
            if self.current_time == self.ASDdata.number_time_steps-1:
                self.iren.DestroyTimer(self.timerId)
            # -------------------------------------------------------------------
            # Update the moments
            # -------------------------------------------------------------------
            self.MomActors.UpdateMoments(window=self, ASDdata=self.ASDdata,
                                         ASDGenActors=self.ASDGenActors, renWin=self.renWin)
            # -------------------------------------------------------------------
            # Take a snapshot
            # -------------------------------------------------------------------
            #self.renWin.Render()
            self.Snapshot()
            # -------------------------------------------------------------------
            # increase the current timer
            # -------------------------------------------------------------------
            self.current_time += 1
        elif self.viz_type == 'E':
            # -------------------------------------------------------------------
            # If the current time is the final time of the measurement destroy the
            # timer
            # -------------------------------------------------------------------
            if self.current_time == self.ASDdata.number_time_steps-1:
                self.iren.DestroyTimer(self.timerId)
            # -------------------------------------------------------------------
            # Update the energy
            # -------------------------------------------------------------------
            self.EneActors.UpdateEnergy(window=self, ASDdata=self.ASDdata,
                                        ASDGenActors=self.ASDGenActors, renWin=self.renWin)
            # -------------------------------------------------------------------
            # Take a snapshot
            # -------------------------------------------------------------------
            self.Snapshot()
            # -------------------------------------------------------------------
            # increase the current timer
            # -------------------------------------------------------------------
            self.current_time += 1
        return
    ############################################################################
    # Individual update of the image, either by increasing the timer count
    # by one or by minus one
    ############################################################################

    def UpdateImage(self):
        if self.viz_type == 'M':
            # -------------------------------------------------------------------
            # Update the moments
            # -------------------------------------------------------------------
            self.MomActors.UpdateMoments(window=self, ASDdata=self.ASDdata,
                                         ASDGenActors=self.ASDGenActors, renWin=self.renWin)
        elif self.viz_type == 'E':
            # -------------------------------------------------------------------
            # Update the energy
            # -------------------------------------------------------------------
            self.EneActors.UpdateEnergy(window=self, ASDdata=self.ASDdata,
                                        ASDGenActors=self.ASDGenActors, renWin=self.renWin)
        return
    ############################################################################
    # Select the energy actor
    ############################################################################

    def toggle_EneActor(self):
        if self.sender() == self.EneDensButton and self.EneDensButton.isChecked():
            self.EneActors.EneDensActor.VisibilityOn()
            self.EneActors.EneActor.VisibilityOff()
            self.renWin.Render()
        if self.sender() == self.EneSiteGlyphs and self.EneSiteGlyphs.isChecked():
            self.EneActors.EneDensActor.VisibilityOff()
            self.EneActors.EneActor.VisibilityOn()
            self.renWin.Render()
        return
    ############################################################################
    # Update the UI
    ############################################################################

    def UpdateRenderer(self):
        if self.sender() == self.SpinsBox:
            if self.SpinsBox.isChecked():
                self.SpinGlyphSelectBox.setEnabled(True)
            else:
                self.SpinGlyphSelectBox.setEnabled(False)
        self.renWin.Render()
        return

    ############################################################################
    # Interactive Simulations and Dock
    ############################################################################    

    def SetSDSliderValue(self, NSimulations):
        self.IntSDSliderVal.setText(f'Simulations: {10*NSimulations}')

    def SetMCSliderValue(self, NSimulations):
        self.IntMCSliderVal.setText(f'Simulations: {10*NSimulations}')      
    
    def IntButtons(self):
        """
        Run a number of simulations depending slider inputs from the user using 
        one of two simulation modes.
        """
        import ASD_GUI.UI.ASDInteractiveTab as IntTab
        #import uppasd as asd
        if self.sender() == self.IntSStepButton:
            IntTab.UpdateIntInputs(self)
            for _ in range(10*self.IntSDSlider.value()):
                self.InteractiveVtk.S_Step()
        if self.sender() == self.IntMCSimButton:
            IntTab.UpdateIntInputs(self)
            for __ in range(10*self.IntMCSlider.value()):
                self.InteractiveVtk.M_step()
        if self.sender() == self.IntResetButton:
            self.InteractiveVtk.Reset()

    def UpdateInteractiveVtk(self):
        """Update text in the interactive window. """
        import ASD_GUI.UI.ASDInteractiveTab as IntTab

        IntTab.UpdateIntInputs(self)

        self.InteractiveVtk.UpdateTemperature()
        self.InteractiveVtk.UpdateBfield()

    def InteractiveScreenshot(self):
        self.InteractiveVtk.Screenshot()

    def CheckForInteractorFiles(self):
        """
        Check if we have any input/output files and determines if 
        the interactive window is ready to launch. Shows error window
        with missing files if check fails. 

        Returns:
                Check   :   Bool
        """

        import os.path as path
        #import uppasd as asd
        import glob
        import ASD_GUI.UI.ASDInputWindows as ASDInputWindows

        restartfile, coordfile = 'dummystring', 'dummystring'
        Check = False

        if len(glob.glob("restart.????????.out")) > 0:
            restartfile = glob.glob("restart.????????.out")[0]    
        if len(glob.glob("coord.????????.out")) > 0:
            coordfile = glob.glob("coord.????????.out")[0]
        if len(self.ASDInputGen.posfile) == 0 and path.exists('posfile'):
            self.ASDInputGen.posfile = glob.glob('posfile')[0]
        if len(self.ASDInputGen.momfile) == 0 and path.exists('momfile'):
            self.ASDInputGen.momfile = glob.glob('momfile')[0]
            
        InputChecklist = [path.exists('inpsd.dat'), path.exists(self.ASDInputGen.posfile),
                           path.exists(self.ASDInputGen.momfile)]
        OutputChecklist = [path.exists(restartfile), path.exists(coordfile)]

        if  all(x == True for x in InputChecklist) and any(x == False for x in OutputChecklist):
            print('Input found, but no output. Running uppasd...')
            try:
                import uppasd as asd
                asd.pyasd.runuppasd()
            except:
                pass
            Check = True

        if  all(x == True for x in OutputChecklist) and all(x == True for x in InputChecklist):
            Check = True

        # Error message 
        Files = ['inpsd.dat', 'posfile', 'momfile']
        MissingFiles = ', '.join([file for index, file in enumerate(Files) if not InputChecklist[index]])

        if not Check:
            self.InteractiveErrorWindow = ASDInputWindows.Error_Window()
            self.InteractiveErrorWindow.FunMsg.setText('Task failed successfully!')
            self.InteractiveErrorWindow.ErrorMsg.setText(f'Could not launch interactive simulation. Missing input files: {MissingFiles}')
            self.InteractiveErrorWindow.show()
        return Check
    
    ###########################################################################

    def RunSimulation(self):
        """
        Run simulation using the uppasd module. Checks if a inpsd.dat file
        is present in current directory and asks for overwrite. 

        Author: Erik Karpelin
        """

        #import uppasd as asd 
        import os.path as path

        if not path.isfile('inpsd.dat'):
            print('inpsd.dat not found, creating from asd_gui')
            self.WriteInputFile()
        try:
            import uppasd as asd 
            asd.pyasd.runuppasd()
        except:
            pass
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
        self.ASDInputGen.import_system(self)
