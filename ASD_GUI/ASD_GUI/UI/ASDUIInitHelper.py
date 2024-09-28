"""
ASDUIInitHelper.py

This module initializes and sets up the user interface for the ASD GUI application.




    Args:
        window (QMainWindow): The main window of the application.


    Args:
        window (QMainWindow): The main window of the application.
"""
# pylint: disable=invalid-name, no-name-in-module, no-member
import os
from enum import Enum

from PyQt6 import uic
from PyQt6.QtGui import QDoubleValidator, QIntValidator
from PyQt6.QtWidgets import QToolBar, QVBoxLayout

# from matplotlib.backends.backend_qt5agg import FigureCanvas

from ASD_GUI.UI.ASDMenuToolbar import (
    Input_Toolbar_Setup,
    InteractiveDock_Setup,
    Plot_Menu_and_Toolbar_Setup,
    VTK_Menu_and_Toolbar_Setup,
)


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
# @brief Initialize the UI and set the relevant actions
# @details Initialize the UI and set the relevant actions. Defines the Toolbars
# and calls for their initialization and population, as well as the reading of the
# .ui file defining the properties of the window.
# Also sets up several validators to forbid erroneous data to be fed into the
# GUI.
# @author Jonathan Chico
##########################################################################
def SetupUI(window):
    """
    Set up the user interface and connect signals to slots.
    """
    window.VTKToolBar = QToolBar()
    window.MatPlotToolbar = QToolBar()
    window.InputToolbar = QToolBar()
    # -----------------------------------------------------------------------
    # Set up UI from Designer file
    # -----------------------------------------------------------------------
    path = os.path.dirname(os.path.abspath(__file__))
    uic.loadUi(os.path.join(path, "ASD_Viz.ui"), window)
    window.chooseBackend()
    window.ModeSelector.currentChanged.connect(window.chooseBackend)
    Plot_Menu_and_Toolbar_Setup(window)
    VTK_Menu_and_Toolbar_Setup(window)
    Input_Toolbar_Setup(window)
    InteractiveDock_Setup(window)
    window.CorrOptsBox.setEnabled(True)
    window.NeighValidator = QIntValidator()
    window.IntegerValidator = QIntValidator()
    window.IntegerValidator.setRange(0, 99999999)
    window.PosDoubleValidator = QDoubleValidator()
    window.PosDoubleValidator.setRange(0, 99999999.9999)
    window.PosDoubleValidator.setDecimals(10)
    window.DoubleValidator = QDoubleValidator()
    window.DoubleValidator.setDecimals(10)
    window.ASDInputGen.ASDInputConstrainer(window)
    window.InpSqwSCStep.setValidator(window.IntegerValidator)
    window.InpPlotDt.setValidator(window.PosDoubleValidator)
    window.AMSDisplayLayout = QVBoxLayout()
    window.AMSDisplayOpts.setLayout(window.AMSDisplayLayout)
    window.ResetUI()
    window.ASDInputGen.ASDSetDefaults()
    window.PosfileWindow.InpPosDone.clicked.connect(window.update_names)
    window.MomfileWindow.InpMomDone.clicked.connect(window.update_names)
    window.JfileWindow.InpJfileDone.clicked.connect(window.update_names)
    window.DMfileWindow.InpDMfileDone.clicked.connect(window.update_names)
    window.JfileWindow.InJfileGenVectors.clicked.connect(
        lambda: window.JfileWindow.GenerateVectorsFromCell(window)
    )
    window.DMfileWindow.InDMfileGenVectors.clicked.connect(
        lambda: window.DMfileWindow.GenerateVectorsFromCell(window)
    )
    window.RestartWindow.InpRestAppendButton.clicked.connect(
        window.create_restart)
    window.RestartWindow.InpRestartDone.clicked.connect(window.create_restart)
    window.RestartWindow.InpRestartDone.clicked.connect(window.update_names)
    window.InitPhaseWindow.InitPhaseDoneButton.clicked.connect(
        window.getInitPhase)
    return

##########################################################################
# Initialization of some of the UI properties
##########################################################################


def InitUI(window):
    """
    Initializes the user interface components and sets their initial states.
    """
    #if True:
    if not window.initialized:
        window.EneMainBox.setEnabled(False)
        window.CamMainBox.setEnabled(False)
        window.MagMainGroup.setEnabled(False)
        window.NeighMainBox.setEnabled(False)
        window.SceneOptMainBox.setEnabled(False)
        window.SpinGlyphSelectBox.setEnabled(True)
        window.ClippBox.setEnabled(False)
        window.ClippBox.setChecked(False)
        window.ClusBox.setVisible(False)
        window.KMCCheck.setVisible(False)
        window.SceneOptMainBox.setEnabled(True)
        window.CamMainBox.setEnabled(True)
        window.actionSave_pov.setEnabled(True)
        window.actionSave_png.setEnabled(True)
        window.ClippBox.setEnabled(True)
        window.actionDisplayMagDens.setEnabled(False)
        window.actionX_ProjMagDens.setEnabled(False)
        window.actionY_ProjMagDens.setEnabled(False)
        window.actionZ_ProjMagDens.setEnabled(False)
    window.file_names[0] = window.ASDdata.posfiles
    window.file_names[1] = window.ASDdata.magnetization
    window.file_names[2] = window.ASDdata.kmcfiles
    window.file_names[3] = window.ASDdata.structfiles
    window.file_names[4] = window.ASDdata.enefiles
    window.file_names[5] = window.ASDdata.dmdatafiles
    window.ProgressBar.setValue(0)
    window.ProgressLabel.setText(f"   {int(window.ProgressBar.value())}%")
    window.initialized = True
    return
