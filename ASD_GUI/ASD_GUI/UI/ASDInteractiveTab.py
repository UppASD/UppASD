"""
This file contains the code for the InteractiveTab in the GUI.
The InteractiveTab provides options for interactive simulation and visualization.
"""
from PyQt6.QtWidgets import (
    QWidget,
    QDockWidget,
    QHBoxLayout,
    QWidget,
    QToolBox,
    QGroupBox,
    QVBoxLayout,
    QPushButton,
    QSlider,
    QLabel,
    QLineEdit,
    QGridLayout,
    QDoubleSpinBox,
    QCheckBox
)
from PyQt6.QtGui import QDoubleValidator, QIntValidator
from PyQt6.QtCore import Qt
import numpy as np


class QTextEditLogger:
    def __init__(self, text_edit):
        self.text_edit = text_edit

    def write(self, message):
        self.text_edit.append(message)
        self.text_edit.ensureCursorVisible()  # Optional: Scroll to the end

    def flush(self):
        pass  # Needed for compatibility with the `sys.stdout` interface
 

def InteractiveDock(window):
    """
    Dock creator for the Interactive page of the GUI.

    Input:
            window      :   QMainWindow

    Returns:
            self.InDock :   QDockWidget
    """
    window.IntDock = QDockWidget("Options", window)
    Docklayout = QHBoxLayout()
    window.IntContents = QWidget()
    window.IntContents.setLayout(Docklayout)
    window.IntToolBox = QToolBox()

    ############################################################################
    # General Inputs
    ############################################################################
    window.IntGeneralBox = QGroupBox()
    GeneralLayout = QVBoxLayout()
    window.TemperatureMagBox = QGroupBox()
    TemperatureMagLayout = QGridLayout()
    window.IntGeneralBox.setTitle("General")
    window.IntGeneralBox.setMaximumSize(300, 900)

    # Temp + Bfield
    window.IntTempLabel = QLabel("Temp")
    window.IntB_xLabel = QLabel("B_x")
    window.IntB_yLabel = QLabel("B_y")
    window.IntB_zLabel = QLabel("B_z")

    PosReal = QDoubleValidator()
    BfieldValidator = QDoubleValidator()
    PosInts = QIntValidator()
    PosReal.setRange(0.0, 99999.9)
    BfieldValidator.setRange(-99999.9, 99999.9)
    PosInts.setRange(0, 99999)

    window.IntTempLine = QDoubleSpinBox()
    window.IntTempLine.setRange(0.0, 99999.9)
    window.IntTempLine.setDecimals(3)
    # window.IntTempLine = QLineEdit()

    # window.IntTempLine.setValidator(PosReal)

    # window.IntB_xLine = QLineEdit()
    # window.IntB_yLine = QLineEdit()
    # window.IntB_zLine = QLineEdit()
    window.IntB_xLine = QDoubleSpinBox()
    window.IntB_yLine = QDoubleSpinBox()
    window.IntB_zLine = QDoubleSpinBox()
    window.IntB_xLine.setDecimals(3)
    window.IntB_yLine.setDecimals(3)
    window.IntB_zLine.setDecimals(3)
    window.IntB_xLine.setRange(-99999.9, 99999.9)
    window.IntB_yLine.setRange(-99999.9, 99999.9)
    window.IntB_zLine.setRange(-99999.9, 99999.9)
    # window.IntB_xLine.setValidator(BfieldValidator)
    # window.IntB_yLine.setValidator(BfieldValidator)
    # window.IntB_zLine.setValidator(BfieldValidator)

    TemperatureMagList = np.array(
        [
            window.IntTempLabel,
            window.IntTempLine,
            window.IntB_xLabel,
            window.IntB_xLine,
            window.IntB_yLabel,
            window.IntB_yLine,
            window.IntB_zLabel,
            window.IntB_zLine,
        ]
    ).reshape(4, 2)

    for row_i, row in enumerate(TemperatureMagList):
        for column_i, element in enumerate(row):
            TemperatureMagLayout.addWidget(element, row_i, column_i)

    window.TemperatureMagBox.setLayout(TemperatureMagLayout)

    # Buttons
    window.IntResetButton = QPushButton()
    window.IntResetButton.setText("Reset Simulation")

    window.IntMomentButton = QPushButton()
    window.IntMomentButton.setText("Load moments")

    window.IntScreenshotBox = QGroupBox()
    window.IntScreenshot = QPushButton()
    window.IntScreenshot.setText("Screenshot")
    
    # Screenshot checkbox
    window.IntScreenshotBox = QGroupBox()
    ScreenshotLayout = QHBoxLayout()
    window.IntScreenshotTic = QCheckBox("Continuous screenshots")
    ScreenshotLayout.addWidget(window.IntScreenshot)
    ScreenshotLayout.addWidget(window.IntScreenshotTic)
    window.IntScreenshotBox.setLayout(ScreenshotLayout)

    ############################################################################
    # Spin Dynamics inputs
    ############################################################################
    window.IntSDBox = QGroupBox()
    SDLayout = QVBoxLayout()
    SliderBox = QGroupBox()
    SliderLayout = QVBoxLayout()
    SliderBox.setMaximumSize(300, 75)
    window.IntSDBox.setTitle("Spin Dynamics")
    window.IntSDBox.setMaximumSize(300, 325)

    # Buttons + slider
    window.IntSStepButton = QPushButton()
    window.IntSStepButton.setText("Run LLG")
    # window.IntSStepButton.setText('Run Simulation')
    SDLayout.addWidget(window.IntSStepButton)

    window.IntSDSlider = QSlider(Qt.Orientation.Horizontal)
    window.IntSDSliderVal = QLabel("", window)
    window.IntSDSlider.setTickPosition(QSlider.TickPosition.TicksBelow)
    window.IntSDSlider.setRange(1, 100)
    window.IntSDSlider.setTickInterval(10)
    SliderLayout.addWidget(window.IntSDSlider)
    SliderLayout.addWidget(window.IntSDSliderVal)
    SliderBox.setLayout(SliderLayout)
    SDLayout.addWidget(SliderBox)
    window.SetSDSliderValue(1)

    # Steps and stepsizes
    window.IntSDStepBox = QGroupBox()
    StepBoxLayout = QGridLayout()
    window.IntSDStepBox.setTitle("Simulation steps")
    window.IntSDSteps = QLineEdit()
    window.IntSDStepSize = QLineEdit()
    window.IntSDSteps.setValidator(PosInts)
    window.IntSDStepSize.setValidator(PosReal)
    StepLabel = QLabel("Steps")
    TimeStepLabel = QLabel("Time step [s]")

    StepBoxLayout.addWidget(StepLabel, 0, 0)
    StepBoxLayout.addWidget(window.IntSDSteps, 0, 1)
    StepBoxLayout.addWidget(TimeStepLabel, 1, 0)
    StepBoxLayout.addWidget(window.IntSDStepSize, 1, 1)

    window.IntSDStepBox.setLayout(StepBoxLayout)
    SDLayout.addWidget(window.IntSDStepBox)
    window.IntSDStepBox.setMaximumSize(300, 100)

    window.IntSDBox.setLayout(SDLayout)
    # window.IntToolBox.addItem(window.IntSDBox, 'Spin Dynamics')

    ############################################################################
    # MonteCarlo inputs
    ############################################################################
    window.IntMCBox = QGroupBox()
    MCLayout = QVBoxLayout()
    MCSliderBox = QGroupBox()
    MCSliderLayout = QVBoxLayout()
    MCSliderBox.setMaximumSize(300, 75)
    window.IntMCBox.setTitle("Monte-Carlo")
    window.IntMCBox.setMaximumSize(300, 275)

    # Buttons and slider
    # window.IntMCSimButton = QPushButton()
    # window.IntMCSimButton.setText('Run Simulation')
    # MCLayout.addWidget(window.IntMCSimButton)
    window.IntMCMSimButton = QPushButton()
    window.IntMCMSimButton.setText("Run Metropolis")
    window.IntMCHSimButton = QPushButton()
    window.IntMCHSimButton.setText("Run Heat Bath")
    MCLayout.addWidget(window.IntMCMSimButton)
    MCLayout.addWidget(window.IntMCHSimButton)

    window.IntMCSlider = QSlider(Qt.Orientation.Horizontal)
    window.IntMCSliderVal = QLabel("", window)
    window.IntMCSlider.setTickPosition(QSlider.TickPosition.TicksBelow)
    window.IntMCSlider.setRange(1, 100)
    window.IntMCSlider.setTickInterval(10)
    MCSliderLayout.addWidget(window.IntMCSlider)
    MCSliderLayout.addWidget(window.IntMCSliderVal)
    MCSliderBox.setLayout(MCSliderLayout)
    MCLayout.addWidget(MCSliderBox)
    window.SetMCSliderValue(1)

    # Steps and stepsizes
    window.IntMCStepBox = QGroupBox()
    MCStepBoxLayout = QGridLayout()
    window.IntMCStepBox.setTitle("Simulation steps")

    window.IntMCSteps = QLineEdit()
    window.IntMCSteps.setValidator(PosInts)
    MCStepLabel = QLabel("Steps")

    MCStepBoxLayout.addWidget(MCStepLabel, 0, 0)
    MCStepBoxLayout.addWidget(window.IntMCSteps, 0, 1)

    window.IntMCStepBox.setLayout(MCStepBoxLayout)
    MCLayout.addWidget(window.IntMCStepBox)
    window.IntMCStepBox.setMaximumSize(300, 100)

    window.IntMCBox.setLayout(MCLayout)
    # window.IntToolBox.addItem(window.IntMCBox, 'Monte-Carlo')
    GeneralWidgetList = [
        window.IntResetButton,
        window.IntMomentButton,
        window.IntScreenshotBox,
        window.TemperatureMagBox,
        window.IntSDBox,
        window.IntMCBox,
    ]

    for w in GeneralWidgetList:
        GeneralLayout.addWidget(w)

    window.IntGeneralBox.setLayout(GeneralLayout)
    window.IntToolBox.addItem(window.IntGeneralBox, "")

    ############################################################################
    # Add everything to the Dock
    ############################################################################
    window.IntDock.setEnabled(True)
    window.IntDock.setVisible(True)
    window.IntDock.setMinimumSize(350, 606)
    window.IntDock.setMaximumSize(350, 1000)

    Docklayout.addWidget(window.IntToolBox)
    window.IntDock.setWidget(window.IntContents)
    window.addDockWidget(Qt.DockWidgetArea.RightDockWidgetArea, window.IntDock)
    window.IntDock.setParent(window)

    return window.IntDock


def InitializeInteractor(window):
    """
    Set up the interactive window based on the files in the
    current directory.

    Inputs:
            window          :   QMainWindow
            InteractiveVtk  :   InteractiveASD object, defined in interactiveASD.py.
    """
    if window.InteractiveVtk.asd is None:
        print("ASD object not found")
        return
    
    window.InteractiveVtk.Launch()

    CurrentSDStep = str(window.InteractiveVtk.asd.inputdata.get_nstep())
    window.IntSDSteps.setPlaceholderText(CurrentSDStep)

    CurrentMCStep = str(window.InteractiveVtk.asd.inputdata.get_mcnstep())
    window.IntMCSteps.setPlaceholderText(CurrentMCStep)

    CurrentTemp = window.InteractiveVtk.asd.inputdata.get_temp()
    # window.IntTempLine.setPlaceholderText(CurrentTemp)
    window.IntTempLine.setValue(CurrentTemp)

    CurrentTimeStep = str(window.InteractiveVtk.asd.inputdata.get_delta_t())
    window.IntSDStepSize.setPlaceholderText(CurrentTimeStep)

    CurrentMagField = window.InteractiveVtk.asd.inputdata.get_hfield()

    window.IntB_xLine.setValue(CurrentMagField[0])
    window.IntB_yLine.setValue(CurrentMagField[1])
    window.IntB_zLine.setValue(CurrentMagField[2])


def UpdateIntInputs(window):
    """
    Update uppasd inputdata based on inputs from the GUI.

    Inputs:
            window  :   QMainWindow
    """
    if window.InteractiveVtk.asd is None:
        return

    if len(window.IntTempLine.text()) > 0:
        NewTemp = float(window.IntTempLine.text())
        window.InteractiveVtk.asd.inputdata.update_temp(NewTemp)

    if len(window.IntSDSteps.text()) > 0:
        NewStep = int(window.IntSDSteps.text())
        window.InteractiveVtk.asd.inputdata.update_nstep(NewStep)

    if len(window.IntMCSteps.text()) > 0:
        NewMCStep = int(window.IntMCSteps.text())
        try:
            window.InteractiveVtk.asd.inputdata.update_mcnstep(NewMCStep)
        except ValueError:
            pass

    if len(window.IntSDStepSize.text()) > 0:
        NewTimStep = float(window.IntSDStepSize.text())
        try:
            window.InteractiveVtk.asd.inputdata.update_delta_t(NewTimStep)
        except ValueError:
            pass

    MagArray = [window.IntB_xLine, window.IntB_yLine, window.IntB_zLine]
    CurrentMagField = window.InteractiveVtk.asd.inputdata.get_hfield()
    NewMagField = []
    for component_i, component in enumerate(MagArray):
        if len(component.text()) > 0:
            NewMagField.append(float(component.text()))
        else:
            NewMagField.append(CurrentMagField[component_i])
    window.InteractiveVtk.asd.inputdata.update_hfield(np.array(NewMagField))
