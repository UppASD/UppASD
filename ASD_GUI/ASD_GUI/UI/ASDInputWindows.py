"""@package ASDInputWindows
Class to define the different auxiliary windows used for the writing of the inpsd.dat

Has the window classes for:
- Creation fo the restartfile.
- Creation of the momfile.
- Creation of the posfile.
- Creation of the initial phase.

Author
----------
Jonathan Chico
"""

import os

import numpy as np
import pandas as pd
from PyQt6 import uic
from PyQt6.QtCore import QSignalBlocker
from PyQt6.QtGui import QIntValidator
from PyQt6.QtWidgets import QDialog, QFileDialog, QLineEdit, QTableWidgetItem

import ASD_GUI.Extras.nn_list_maker.create_neighbour_list as create_neighbour_list
import ASD_GUI.Input_Creator.ASDInputGen as ASDInputgen
from ASD_GUI.Extras.nn_list_maker.read_uppasd import read_posfile
from ASD_GUI.Extras.nn_list_maker.structure import get_full_nnlist
from ASD_GUI.Input_Creator.ASDInputAux import (
    create_coord,
    create_spiral,
    write_domain_wall,
    write_skyrmion,
)


class InitPhaseWindow(QDialog):
    """Class responsible for the creation of the init phase window, this window contains
    tables that depending on the type of input that one is going to generate sets up
    the initial phase for either MC or ASD.

    Author
    ----------
    Jonathan Chico
    """

    def __init__(self, parent=None):

        super(InitPhaseWindow, self).__init__(parent)
        path = os.path.dirname(os.path.abspath(__file__))
        uic.loadUi(os.path.join(path, "InitPhase_Creator.ui"), self)
        self.InitPhaseAddButton.clicked.connect(self.table_control)
        self.InitPhaseDelButton.clicked.connect(self.table_control)
        self.InitPhaseDoneButton.clicked.connect(self.window_close)
        self.InitPhaseCancelButton.clicked.connect(self.window_close)
        return

    def table_control(self):
        """
        Function to control the addition and removal of rows in the table
        defining the initial phase.
        The user can on runtime add or delete rows until a
        minimum of one row remains. New rows are created with dummy text in them.
        TODO: Must add validators to the cells so that only the appropriate 
            data types can be fed in into the cells.
        Author
        ----------
        Jonathan Chico
        """
        # If the initial phase is MC
        if self.MCannealBox.isEnabled():
            if self.sender() == self.InitPhaseAddButton:
                rowMcanneal = self.MCannealTable.rowCount()
                self.MCannealTable.insertRow(rowMcanneal)
                text = [1, 0.0]
                for ii, value in enumerate(text):
                    item = QTableWidgetItem()
                    item.setText(str(value))
                    self.MCannealTable.setItem(rowMcanneal, ii, item)
            if self.sender() == self.InitPhaseDelButton:
                rowMcanneal = self.MCannealTable.rowCount()
                # Make sure that one cannot delete the last entry
                if rowMcanneal > 1:
                    self.MCannealTable.removeRow(rowMcanneal - 1)
        # If the initial phase is LLG
        if self.IpNphaseBox.isEnabled():
            if self.sender() == self.InitPhaseAddButton:
                rownphase = self.IpNphaseTable.rowCount()
                self.IpNphaseTable.insertRow(rownphase)
                text = [1, 0.0, 1e-16, 0.001]
                for ii, value in enumerate(text):
                    item = QTableWidgetItem()
                    item.setText(str(value))
                    self.IpNphaseTable.setItem(rownphase, ii, item)
            if self.sender() == self.InitPhaseDelButton:
                rownphase = self.IpNphaseTable.rowCount()
                # Make sure that one cannot delete the last entry
                if rownphase > 1:
                    self.IpNphaseTable.removeRow(rownphase - 1)
        return

    def window_close(self):
        """Function to close the window defining the initial phase and storing the
        needed data for passing it to the input file.

        Author
        ----------
        Jonathan Chico
        """

        # ----------------------------------------------------------------------------
        # Clean the table
        # ----------------------------------------------------------------------------
        if self.sender() == self.InitPhaseCancelButton:
            for ii in range(1, self.MCannealTable.rowCount()):
                self.MCannealTable.removeRow(ii)
            for ii in range(1, self.IpNphaseTable.rowCount()):
                self.IpNphaseTable.removeRow(ii)
            self.close()
        # ----------------------------------------------------------------------------
        # Populate the needed data
        # ----------------------------------------------------------------------------
        if self.sender() == self.InitPhaseDoneButton:
            self.init_phase_data = []
            if self.MCannealBox.isEnabled():
                self.init_phase_data.append(self.MCannealTable.rowCount())
                for row in range(0, self.MCannealTable.rowCount()):
                    tmp = []
                    for col in range(0, self.MCannealTable.columnCount()):
                        if col < 1:
                            tmp.append(int(self.MCannealTable.item(row, col).text()))
                        else:
                            tmp.append(float(self.MCannealTable.item(row, col).text()))
                    self.init_phase_data.append(tmp)
            if self.IpNphaseTable.isEnabled():
                self.init_phase_data.append(self.IpNphaseTable.rowCount())
                for row in range(0, self.IpNphaseTable.rowCount()):
                    tmp = []
                    for col in range(0, self.IpNphaseTable.columnCount()):
                        if col < 1:
                            tmp.append(int(self.IpNphaseTable.item(row, col).text()))
                        else:
                            tmp.append(float(self.IpNphaseTable.item(row, col).text()))
                    self.init_phase_data.append(tmp)
            if self.IpVPOBox.isEnabled():
                for row in range(0, self.VPOTable.rowCount()):
                    tmp = []
                    for col in range(0, self.VPOTable.columnCount()):
                        if col < 1:
                            tmp.append(int(self.VPOTable.item(row, col).text()))
                        else:
                            tmp.append(float(self.VPOTable.item(row, col).text()))
                    self.init_phase_data.append(tmp)
            self.close()
        return


################################################################################
# @brief Class responsible for the creation of error windows.
# @details Class responsible for the creation of error windows, this window contains
# a placeholder message that can be modified to display the precise error message
# @author Jonathan Chico
################################################################################
class ErrorWindow(QDialog):
    """Class responsible for the creation of error windows, this window contains
    a placeholder message that can be modified to display the precise error message.

    Author
    ----------
    Jonathan Chico
    """

    def __init__(self, parent=None):

        super(ErrorWindow, self).__init__(parent)
        path = os.path.dirname(os.path.abspath(__file__))
        uic.loadUi(os.path.join(path, "Error_Msg.ui"), self)
        self.ErrorMsgAccept.clicked.connect(self.close)
        return


################################################################################
# @brief Class responsible for the creation of information windows.
# @details Class responsible for the creation of information windows, this
# window contains a placeholder message that can be modified to display
# a more precise information message (not error)
# @author Jonathan Chico, Anders Bergman
################################################################################
class InfoWindow(QDialog):
    """Class responsible for the creation of information windows, this
    window contains a placeholder message that can be modified to display
    a more precise information message (not error)

    Author
    ----------
    Jonathan Chico, Anders Bergman
    """

    def __init__(self, parent=None):

        super(InfoWindow, self).__init__(parent)
        path = os.path.dirname(os.path.abspath(__file__))
        uic.loadUi(os.path.join(path, "Info_Msg.ui"), self)
        self.InfoMsgAccept.clicked.connect(self.close)
        return


################################################################################
## @brief Class containing the definitions necessary for the creation of the window
# responsible for the restartfile creation. As well as the definition of the
# actions associated with that window.
#
# @author Jonathan Chico
################################################################################
class RestartWindow(QDialog):
    """Class containing the definitions necessary for the creation of the window
    responsible for the restartfile creation. As well as the definition of the
    actions associated with that window.

    Author
    ----------
    Jonathan Chico
    """

    def __init__(self, parent=None):
        """Class constructor for the restart window, it handles the loading of the
        UI as well as the setting of the actions for the objects in the window.

        Author
        ----------
        Jonathan Chico
        """

        super(RestartWindow, self).__init__(parent)
        # ----------------------------------------------------------------------------
        # Load UI
        # ----------------------------------------------------------------------------
        path = os.path.dirname(os.path.abspath(__file__))
        uic.loadUi(os.path.join(path, "Restart_Creator.ui"), self)
        # ----------------------------------------------------------------------------
        # Connecting actions
        # ----------------------------------------------------------------------------
        self.DWOptBox.clicked.connect(self.restart_selector)
        self.SkxOptBox.clicked.connect(self.restart_selector)
        self.HLOptBox.clicked.connect(self.restart_selector)
        self.InpDwTypeBox.clicked.connect(self.restart_selector)
        self.InpDWVortex.clicked.connect(self.restart_selector)
        self.InpRestartCancel.clicked.connect(self.close)
        self.InpDWPosN1Slider.valueChanged.connect(self.update_sliders_DW)
        self.InpDWPosN1.editingFinished.connect(self.update_sliders_DW)
        self.InpDWPosN2Slider.valueChanged.connect(self.update_sliders_DW)
        self.InpDWPosN2.editingFinished.connect(self.update_sliders_DW)
        self.InpDWPosN3Slider.valueChanged.connect(self.update_sliders_DW)
        self.InpDWPosN3.editingFinished.connect(self.update_sliders_DW)
        self.InpRestFullEnsButton.clicked.connect(self.UpdateRestartUI)
        self.InpRestSingleEnsButton.clicked.connect(self.UpdateRestartUI)
        # ----------------------------------------------------------------------------
        # Setting validators
        # ----------------------------------------------------------------------------
        self.InpSkyrOrder.setValidator(QIntValidator())
        # ----------------------------------------------------------------------------
        # Initializing counters
        # ----------------------------------------------------------------------------
        self.curr_image = 0
        return

    def UpdateRestartUI(self):
        if self.InpRestSingleEnsButton.isChecked():
            self.InpRestAppendButton.setEnabled(True)
        else:
            self.InpRestAppendButton.setEnabled(False)
        return

    ############################################################################
    ## @brief Function that ensures that only one type of restartfile can be selected
    # at a given time in the restartfile window.
    # @author Jonathan Chico
    ############################################################################
    def restart_selector(self):
        """Function that ensures that only one type of restartfile can be selected
        at a given time in the restartfile window.
        Author
        ----------
        Jonathan Chico
        """

        if self.sender() == self.DWOptBox:
            if self.DWOptBox.isChecked():
                self.SkxOptBox.setChecked(False)
                QSignalBlocker(self.SkxOptBox)
                self.HLOptBox.setChecked(False)
                QSignalBlocker(self.HLOptBox)
        if self.sender() == self.SkxOptBox:
            if self.SkxOptBox.isChecked():
                self.DWOptBox.setChecked(False)
                QSignalBlocker(self.DWOptBox)
                self.HLOptBox.setChecked(False)
                QSignalBlocker(self.HLOptBox)
        if self.sender() == self.HLOptBox:
            if self.HLOptBox.isChecked():
                self.SkxOptBox.setChecked(False)
                QSignalBlocker(self.SkxOptBox)
                self.DWOptBox.setChecked(False)
                QSignalBlocker(self.DWOptBox)
        if self.sender() == self.InpDwTypeBox:
            if self.InpDwTypeBox.isChecked():
                self.InpDWVortex.setChecked(False)
                QSignalBlocker(self.InpDWVortex)
        if self.sender() == self.InpDWVortex:
            if self.InpDWVortex.isChecked():
                self.InpDwTypeBox.setChecked(False)
                QSignalBlocker(self.InpDwTypeBox)
        return

    ############################################################################
    # @brief Function to prepare for the generation of the restartfile from the
    # GUI.
    # @details It mainly consists in preparing data such as the magnitude of the
    # moments calling the creation of the coordinates and setting up the sliders for the
    # restarfile creator.
    # @author Jonathan Chico
    ############################################################################
    def restart_pre_generation(self, inp_data):
        """Function to prepare for the generation of the restartfile from the
        GUI. It mainly consists in preparing data such as the magnitude of the
        moments calling the creation of the coordinates and setting up the sliders for the
        restarfile creator.

        Args
        ----------
        - inp_data: Class containing the dictionary definitions needed to generate the
        restartfile

        Author
        ----------
        Jonathan Chico

        """

        # Read the posfile as defined
        self.Bas = pd.read_csv(
            inp_data.posfile, header=None, delim_whitespace=True, usecols=[2, 3, 4]
        ).values
        # Read the momfile as defined
        self.mom = pd.read_csv(
            inp_data.momfile, header=None, delim_whitespace=True, usecols=[2]
        ).values
        # Read the lattice vectors
        self.cell = np.asarray(inp_data.UppASDKeywords["geometry"]["cell"])
        # Read the number of repetitions along the different directions
        self.ncell = np.asarray(inp_data.UppASDKeywords["geometry"]["ncell"])
        # Set the number of images
        self.Mensemble = inp_data.UppASDKeywords["Mag"]["Mensemble"]
        # Set the size of the blocks for memory locality and/or macrospin approach
        self.block_size = inp_data.UppASDKeywords["Hamiltonian"]["block_size"]
        # Calculate the atomic coordinates
        (self.coord, self.mom_mag) = create_coord(
            self.cell, self.ncell, self.Bas, self.block_size, self.mom
        )
        # Set the number of atoms in the system
        self.Natom = len(self.coord)
        # Set the values for the sliders
        self.InpDWPosN1Slider.setMinimum(1)
        self.InpDWPosN1Slider.setMaximum(self.ncell[0])
        self.InpDWPosN2Slider.setMinimum(1)
        self.InpDWPosN2Slider.setMaximum(self.ncell[1])
        self.InpDWPosN3Slider.setMinimum(1)
        self.InpDWPosN3Slider.setMaximum(self.ncell[2])
        self.InpDWPosN1Validator = QIntValidator()
        self.InpDWPosN1Validator.setRange(1, self.ncell[0])
        self.InpDWPosN1.setValidator(self.InpDWPosN1Validator)
        self.InpDWPosN2Validator = QIntValidator()
        self.InpDWPosN2Validator.setRange(1, self.ncell[1])
        self.InpDWPosN2.setValidator(self.InpDWPosN2Validator)
        self.InpDWPosN3Validator = QIntValidator()
        self.InpDWPosN3Validator.setRange(1, self.ncell[2])
        self.InpDWPosN3.setValidator(self.InpDWPosN3Validator)
        # Write the coordinate file
        coord_name, _ = QFileDialog.getSaveFileName(self, "Save Coordinate File")
        coord_fmt = "{: 6d}  {: 4.8f}  {: 4.8f}  {: 4.8f}\n"
        coord_file = open(coord_name, "w")
        for iat in range(0, self.Natom):
            coord_file.write(
                coord_fmt.format(
                    iat + 1, self.coord[iat, 0], self.coord[iat, 1], self.coord[iat, 2]
                )
            )
        return

    ############################################################################
    # @brief Wrapper to write the different types of configurations via the restartfile
    # GUI.
    # @details It calls the necessary functions to ensure that a restarfile can
    # be appropriately written.
    ############################################################################
    def write_restartfile(self, inp_data):
        """Wrapper to write the different types of configurations via the restartfile
        GUI. It calls the necessary functions to ensure that a restarfile can
        be appropriately written.

        Author
        ----------
        Jonathan Chico

        """

        if self.sender() == self.InpRestAppendButton:
            self.Mensemble = 1
        # -----------------------------------------------------------------------
        # Write a restartfile with a domain wall
        # -----------------------------------------------------------------------
        if self.DWOptBox.isChecked():
            DWInfo = dict()
            # -------------------------------------------------------------------
            # Domain wall chirality
            # -------------------------------------------------------------------
            if self.DwChirLeft.isChecked():
                DWInfo["chirality"] = -1.0
            if self.DwChirRight.isChecked():
                DWInfo["chirality"] = 1.0
            # -------------------------------------------------------------------
            # Type of domain wall Neel or Bloch
            # -------------------------------------------------------------------
            if self.InpDwNeel.isChecked():
                DWInfo["DWtype"] = "Neel"
            if self.InpDwBloch.isChecked():
                DWInfo["DWtype"] = "Bloch"
            # -------------------------------------------------------------------
            # Planar domain wall
            # -------------------------------------------------------------------
            if self.InpDwTypeBox.isChecked():
                DWInfo["type"] = "planar"
                # ---------------------------------------------------------------
                # Setting the center of the domain wall
                # ---------------------------------------------------------------
                if self.DWPlane_x.isChecked():
                    DWInfo["plane"] = 0
                    Num = int(self.InpDWPosN1.text()) - 1
                    mid = (
                        max(self.Bas[:, DWInfo["plane"]])
                        + min(self.Bas[:, DWInfo["plane"]])
                    ) * 0.5
                    cent = Num * self.cell[DWInfo["plane"], DWInfo["plane"]] + mid
                    DWInfo["center"] = cent
                if self.DWPlane_y.isChecked():
                    DWInfo["plane"] = 1
                    Num = int(self.InpDWPosN2.text()) - 1
                    mid = (
                        max(self.Bas[:, DWInfo["plane"]])
                        + min(self.Bas[:, DWInfo["plane"]])
                    ) * 0.5
                    cent = Num * self.cell[DWInfo["plane"], DWInfo["plane"]] + mid
                    DWInfo["center"] = cent
                if self.DWPlane_z.isChecked():
                    DWInfo["plane"] = 2
                    Num = int(self.InpDWPosN3.text()) - 1
                    mid = (
                        max(self.Bas[:, DWInfo["plane"]])
                        + min(self.Bas[:, DWInfo["plane"]])
                    ) * 0.5
                    cent = Num * self.cell[DWInfo["plane"], DWInfo["plane"]] + mid
                # ---------------------------------------------------------------
                # Setting the magnetization easy axis
                # ---------------------------------------------------------------
                if self.InpDwMagMx.isChecked():
                    DWInfo["easy_axis"] = "x"
                if self.InpDwMagMy.isChecked():
                    DWInfo["easy_axis"] = "y"
                if self.InpDwMagMz.isChecked():
                    DWInfo["easy_axis"] = "z"
                # ---------------------------------------------------------------
                # Setting the domain wall width
                # ---------------------------------------------------------------
                if len(self.InpDWWidth.text()) > 0:
                    DWInfo["width"] = float(self.InpDWWidth.text())
                else:
                    DWInfo["width"] = 1.0
                    print("No width set a value of 1.0 assumed")
            # -------------------------------------------------------------------
            # Vortex domain wall
            # -------------------------------------------------------------------
            elif self.InpDWVortex.isChecked():
                DWInfo["type"] = "vortex"
                # ---------------------------------------------------------------
                # Setting the vortex core radius
                # ---------------------------------------------------------------
                if len(self.InpDWCoreRad.text()) > 0:
                    DWInfo["radius"] = float(self.InpDWCoreRad.text())
                else:
                    DWInfo["radius"] = 1.0
                    print("No radius set a value of 1.0 assumed")
                # ---------------------------------------------------------------
                # Setting the vortex polarity
                # ---------------------------------------------------------------
                DWInfo["polarity"] = int(self.DwCorPol.value())
                # ---------------------------------------------------------------
                # Setting the vortex core center
                # ---------------------------------------------------------------
                DWInfo["center"] = [
                    int(self.DwCorPosN1.text()),
                    int(self.DwCorPosN2.text()),
                    int(self.DwCorPosN3.text()),
                ]
            # -------------------------------------------------------------------
            # Create the domain wall
            # -------------------------------------------------------------------
            mag = write_domain_wall(self.Natom, self.Mensemble, self.coord, DWInfo)
        # -----------------------------------------------------------------------
        # Write a restartfile with an isolated skyrmion
        # -----------------------------------------------------------------------
        if self.SkxOptBox.isChecked():
            SkxInfo = dict()
            # -------------------------------------------------------------------
            # Find the input position of the skyrmion center
            # -------------------------------------------------------------------
            if len(self.SkPosN1.text()):
                N1 = int(self.SkPosN1.text())
            else:
                N1 = 0
            if len(self.SkPosN2.text()):
                N2 = int(self.SkPosN2.text())
            else:
                N2 = 0
            if len(self.SkPosN3.text()):
                N3 = int(self.SkPosN3.text())
            else:
                N3 = 0
            Num = [N1 - 1, N2 - 1, N3 - 1]
            mid = [
                (max(self.Bas[:, 0]) + min(self.Bas[:, 0])) * 0.5,
                (max(self.Bas[:, 1]) + min(self.Bas[:, 1])) * 0.5,
                (max(self.Bas[:, 2]) + min(self.Bas[:, 2])) * 0.5,
            ]
            # -------------------------------------------------------------------
            # Set the skyrmion center position
            # -------------------------------------------------------------------
            SkxInfo["center"] = (
                Num[0] * self.cell[0, :]
                + Num[1] * self.cell[1, :]
                + Num[2] * self.cell[2, :]
                + mid[:]
            )
            # -------------------------------------------------------------------
            # Set the skyrmion type, Bloch or Neel
            # -------------------------------------------------------------------
            if self.SkxTypeBloch.isChecked():
                SkxInfo["type"] = np.pi * 0.5
            if self.SkxTypeNeel.isChecked():
                SkxInfo["type"] = 0
            # -------------------------------------------------------------------
            # Set the core polarity
            # -------------------------------------------------------------------
            if int(self.InpSkyrPol.value()) > 0:
                SkxInfo["polarity"] = 0.0
            else:
                SkxInfo["polarity"] = np.pi
            # -------------------------------------------------------------------
            # Set the skyrmion order, i.e. its topological charge
            # -------------------------------------------------------------------
            if len(self.InpSkyrOrder.text()) > 0:
                SkxInfo["order"] = int(self.InpSkyrOrder.text())
            else:
                SkxInfo["order"] = 1
                print("No skyrmion order chosen, set to 1")
            # -------------------------------------------------------------------
            # Set the skyrmion handness
            # -------------------------------------------------------------------
            if self.InpSkyrHand.currentText() == "Clockwise":
                SkxInfo["handness"] = 1
            else:
                SkxInfo["handness"] = -1
            # -------------------------------------------------------------------
            # Set the skyrmion radius
            # -------------------------------------------------------------------
            if len(self.InpSkxRad.text()) > 0:
                SkxInfo["width"] = float(self.InpSkxRad.text())
            else:
                SkxInfo["width"] = 1.0
            # -------------------------------------------------------------------
            # Create the skyrmion
            # -------------------------------------------------------------------
            mag = write_skyrmion(self.Natom, self.Mensemble, self.coord, SkxInfo)
        if self.HLOptBox.isChecked():
            HLInfo = dict()
            # -------------------------------------------------------------------
            # Read the cone angle
            # -------------------------------------------------------------------
            if len(self.InpHLConeAngle.text()) > 0:
                HLInfo["cone_angle"] = float(self.InpHLConeAngle.text()) * np.pi / 180.0
            else:
                HLInfo["cone_angle"] = 90.0
            if HLInfo["cone_angle"] == 0:
                print("Cone angle set to zero, you will get a Ferromagnet")
            # -------------------------------------------------------------------
            # Read the wave vector of the spiral
            # -------------------------------------------------------------------
            if len(self.HLKX.text()) > 0:
                KX = float(self.HLKX.text())
            else:
                KX = 0.0
            if len(self.HLKY.text()) > 0:
                KY = float(self.HLKY.text())
            else:
                KY = 0.0
            if len(self.HLKZ.text()) > 0:
                KZ = float(self.HLKZ.text())
            else:
                KZ = 0.0
            HLInfo["prop_vector"] = [KX, KY, KZ]
            # -------------------------------------------------------------------
            # Read the pitch vector of the spiral
            # -------------------------------------------------------------------
            if len(self.InpHLPitchX.text()) > 0:
                PX = float(self.InpHLPitchX.text())
            else:
                PX = 0.0
            if len(self.InpHLPitchY.text()) > 0:
                PY = float(self.InpHLPitchY.text())
            else:
                PY = 0.0
            if len(self.InpHLPitchZ.text()) > 0:
                PZ = float(self.InpHLPitchZ.text())
            else:
                PZ = 0.0
            HLInfo["pitch_vector"] = [PX, PY, PZ]
            # -------------------------------------------------------------------
            # Normalize the pitch vector just in case
            # -------------------------------------------------------------------
            norm = np.sqrt(
                HLInfo["pitch_vector"][0] ** 2
                + HLInfo["pitch_vector"][1] ** 2
                + HLInfo["pitch_vector"][2] ** 2
            )
            HLInfo["pitch_vector"] = HLInfo["pitch_vector"] / norm
            # -------------------------------------------------------------------
            # Set the chirality of the spin spiral
            # -------------------------------------------------------------------
            if self.InpHLHandness.currentText() == "Clockwise":
                HLInfo["handness"] = 1
            else:
                HLInfo["handness"] = -1
            # -------------------------------------------------------------------
            # Create the spiral
            # -------------------------------------------------------------------
            mag = create_spiral(self.Natom, self.Mensemble, self.coord, HLInfo)
        # -----------------------------------------------------------------------
        # Restartfile format
        # -----------------------------------------------------------------------
        restart_fmt = "{:8d}{:8d}{:8d}  {: 16.8E}{: 16.8E}{: 16.8E}{: 16.8E}\n"
        # -----------------------------------------------------------------------
        # Set the restart file name
        # -----------------------------------------------------------------------
        if self.curr_image == 0:
            restart_name, _ = QFileDialog.getSaveFileName(self, "Save Restart File")
            # -----------------------------------------------------------------------
            # Save the restarfile name to the input
            # -----------------------------------------------------------------------
            RestartWindow.restartfile_name = restart_name
        # -----------------------------------------------------------------------
        # Actually write the restartfile
        # -----------------------------------------------------------------------
        if len(RestartWindow.restartfile_name) > 0:
            # ------------------------------------------------------------------------
            # If it is the first image generate a new file and write the header
            # ------------------------------------------------------------------------
            if self.curr_image == 0:
                restart = open(RestartWindow.restartfile_name, "w")
                restart.write(f"#" * 55 + "\n")
                restart.write(f"# File type: R\n")
                restart.write(f"# Simulation type: Init\n")
                restart.write(f"# Number of atoms: {self.Natom:8d}\n")
                restart.write(f"# Number of ensembles: {self.Mensemble:8d}\n")
                restart.write("#" * 55 + "\n")
                restart.write(f"{'#iter':>8s}{'ens':>8s}{'iatom':>8s}\
                    {'|Mom|':>16s}{'M_x':>16s}{'M_y':>16s}{'M_z':>16s}\n")
            # ------------------------------------------------------------------------
            # If the image is larger than 1 and one writes each image repone the file to append
            # ------------------------------------------------------------------------
            elif self.curr_image > 0 and self.InpRestSingleEnsButton.isChecked():
                restart = open(RestartWindow.restartfile_name, "a")
            # ------------------------------------------------------------------------
            # Write to the restartfile
            # ------------------------------------------------------------------------
            if self.sender() == self.InpRestartDone:
                if self.InpRestFullEnsButton.isChecked():
                    for ens in range(0, self.Mensemble):
                        for iat in range(0, self.Natom):
                            restart.write(
                                restart_fmt.format(
                                    -1,
                                    ens + 1,
                                    iat + 1,
                                    self.mom_mag[iat],
                                    mag[ens, iat, 0],
                                    mag[ens, iat, 1],
                                    mag[ens, iat, 2],
                                )
                            )
                    restart.close()
            elif self.sender() == self.InpRestAppendButton:
                for ens in range(0, self.Mensemble):
                    for iat in range(0, self.Natom):
                        restart.write(
                            restart_fmt.format(
                                -1,
                                ens + 1,
                                iat + 1,
                                self.mom_mag[iat],
                                mag[ens, iat, 0],
                                mag[ens, iat, 1],
                                mag[ens, iat, 2],
                            )
                        )
                restart.close()
            # ------------------------------------------------------------------------
            # If one is done close the window
            # ------------------------------------------------------------------------
            if self.sender() == self.InpRestartDone:
                self.close()
            # ------------------------------------------------------------------------
            # Increase the image count by one
            # ------------------------------------------------------------------------
            self.curr_image += 1
        else:
            print("No restartfile name given")
        return

    ############################################################################
    # @brief Update the sliders and the line edits for the domain wall position
    # @author Jonathan Chico
    ############################################################################
    def update_sliders_DW(self):
        """Update the sliders and the line edits for the domain wall position

        Author
        ----------
        Jonathan Chico
        """
        if self.sender() == self.InpDWPosN1:
            self.InpDWPosN1Slider.setValue(int(self.InpDWPosN1.text()))
        if self.sender() == self.InpDWPosN1Slider:
            self.InpDWPosN1.setText(str(self.InpDWPosN1Slider.value()))
        if self.sender() == self.InpDWPosN2:
            self.InpDWPosN2Slider.setValue(int(self.InpDWPosN2.text()))
        if self.sender() == self.InpDWPosN2Slider:
            self.InpDWPosN2.setText(str(self.InpDWPosN2Slider.value()))
        if self.sender() == self.InpDWPosN3:
            self.InpDWPosN3Slider.setValue(int(self.InpDWPosN3.text()))
        if self.sender() == self.InpDWPosN3Slider:
            self.InpDWPosN3.setText(str(self.InpDWPosN3Slider.value()))
        return


################################################################################
# @brief Class containing the defintions and actions needed for the display of the
# window handling the creation of the posfile inside the GUI.
# @author Jonathan Chico
################################################################################
class PosfileWindow(QDialog):
    """Class containing the defintions and actions needed for the display of the
    window handling the creation of the posfile inside the GUI.

    Author
    ----------
    Jonathan Chico
    """

    def __init__(self, parent=None):

        super(PosfileWindow, self).__init__(parent)
        path = os.path.dirname(os.path.abspath(__file__))
        uic.loadUi(os.path.join(path, "Posfile_Creator.ui"), self)
        self.InPosAddRow.clicked.connect(self.table_control)
        self.InPosAddRowRand.clicked.connect(self.table_control)
        self.InPosDelRow.clicked.connect(self.table_control)
        self.InPosDelRowRand.clicked.connect(self.table_control)
        self.InpPosCancel.clicked.connect(self.window_close)
        self.InpPosDone.clicked.connect(self.window_close)
        PosfileWindow.posfile_gotten = False
        PosfileWindow.posfile_name = "./posfile.dat"
        return

    ############################################################################
    ## @brief Function to control the addition and removal of rows in the table
    # defining the \c posfile.
    # @details The user can on runtime add or delete rows until a
    # minimum of one row remains. New rows are created with dummy text in them.
    # @todo Must add validators to the cells so that only the appropriate data types
    # can be fed in into the cells.
    # @author Jonathan Chico
    ############################################################################
    def table_control(self):
        """Function to control the addition and removal of rows in the table
        defining the posfile.
        The user can on runtime add or delete rows until a
        minimum of one row remains. New rows are created with dummy text in them.
        Must add validators to the cells so that only the appropriate data types
        can be fed in into the cells.

        Author
        ----------
        Jonathan Chico
        """

        if self.sender() == self.InPosAddRow:
            rowPosition = self.InPosTable.rowCount()
            self.InPosTable.insertRow(rowPosition)
            text = [rowPosition + 1, rowPosition + 1, 0.0, 0.0, 0.0]
            for ii in range(0, len(text)):
                item = QTableWidgetItem()
                item.setText(str(text[ii]))
                self.InPosTable.setItem(rowPosition, ii, item)
        if self.sender() == self.InPosAddRowRand:
            rowPosition = self.InPosTableRand.rowCount()
            self.InPosTableRand.insertRow(rowPosition)
            text = [
                rowPosition + 1,
                rowPosition + 1,
                rowPosition + 1,
                1.0,
                0.0,
                0.0,
                0.0,
            ]
            for ii, value in enumerate(text):
                item = QTableWidgetItem()
                item.setText(str(value))
                self.InPosTableRand.setItem(rowPosition, ii, item)
        if self.sender() == self.InPosDelRow:
            rowPosition = self.InPosTable.rowCount()
            # Make sure that one cannot delete the last entry
            if rowPosition > 1:
                self.InPosTable.removeRow(rowPosition - 1)
        if self.sender() == self.InPosDelRowRand:
            rowPosition = self.InPosTableRand.rowCount()
            # Make sure that one cannot delete the last entry
            if rowPosition > 1:
                self.InPosTableRand.removeRow(rowPosition - 1)
        return

    def CheckForFile(self, mainwindow):
        """If a jfile have already been selected, input it into the creator."""

        if len(mainwindow.ASDInputGen.posfile) > 0:
            posfile = np.genfromtxt(
                mainwindow.ASDInputGen.posfile.split("/")[-1], ndmin=2
            )
            posfile = [list(line) for line in posfile]
            self.posfile_gotten = True

            Table = self.InPosTable
            Table.setRowCount(0)
            for row, line in enumerate(posfile):
                Table.insertRow(row)
                for column, element in enumerate(line):
                    item = QTableWidgetItem(str(element))
                    Table.setItem(row, column, item)

    ############################################################################
    ## @brief Function handling the what the Cancel and Done buttons do in the \c posfile
    # window.
    # @details The Cancel button removes all the rows except for the first one and closes
    # the window.
    # The Done button reads the data in the cells and writes a \c posfile dubbed
    # posfile.dat in \c UppASD format.
    # @author Jonathan Chico
    ############################################################################
    def window_close(self):
        """Function handling the what the Cancel and Done buttons do in the posfile
        window.
        The Cancel button removes all the rows except for the first one and closes
        the window.
        The Done button reads the data in the cells and writes a posfile dubbed
        posfile.dat in UppASD format.

        Author
        ----------
        Jonathan Chico
        """

        if self.sender() == self.InpPosCancel:
            self.InPosTable.setRowCount(1)
            for column, value in enumerate([1, 1, 0.0, 0.0, 0.0]):
                item = QTableWidgetItem(str(value))
                self.InPosTable.setItem(0, column, item)
            self.close()
        if self.sender() == self.InpPosDone:
            posfile_name = open(PosfileWindow.posfile_name, "w")
            if self.InPosBoxRand.isEnabled():
                for row in range(0, self.InPosTableRand.rowCount()):
                    for col in range(0, self.InPosTableRand.columnCount()):
                        if col < 3:
                            entry = int(self.InPosTableRand.item(row, col).text())
                        else:
                            entry = float(self.InPosTableRand.item(row, col).text())
                        posfile_name.write(f"{entry}  ")
                    posfile_name.write(f"{entry}  ")
            if self.InPosBox.isEnabled():
                for row in range(0, self.InPosTable.rowCount()):
                    for col in range(0, self.InPosTable.columnCount()):
                        if col < 2:
                            entry = int(float(self.InPosTable.item(row, col).text()))
                        else:
                            entry = float(self.InPosTable.item(row, col).text())
                        posfile_name.write(f"{entry}  ")
                    posfile_name.write("\n")
            self.close()
            self.posfile_gotten = True
        return


################################################################################
## @brief Class containing the defintions and actions needed for the display of the
# window handling the creation of the momfile inside the GUI.
# @author Jonathan Chico
################################################################################
class MomfileWindow(QDialog):
    """Class containing the defintions and actions needed for the display of the
    window handling the creation of the momfile inside the GUI.

    Author
    ----------
    Jonathan Chico
    """

    def __init__(self, parent=None):

        super(MomfileWindow, self).__init__(parent)
        path = os.path.dirname(os.path.abspath(__file__))
        uic.loadUi(os.path.join(path, "Momfile_Creator.ui"), self)
        self.InMomAddRow.clicked.connect(self.table_control)
        self.InMomDelRow.clicked.connect(self.table_control)
        self.InpMomCancel.clicked.connect(self.window_close)
        self.InpMomDone.clicked.connect(self.window_close)
        MomfileWindow.momfile_gotten = False
        MomfileWindow.momfile_name = "./momfile.dat"
        return

    ############################################################################
    ## @brief Function to control the addition and removal of rows in the table
    # defining the \c momfile.
    # @details The user can on runtime add or delete rows until a minimum of one row remains.
    # New rows are created with dummy text in them.
    # @todo Must add validators to the cells so that only the appropriate data types
    # can be fed in into the cells.
    # @author Jonathan Chico
    ############################################################################
    def table_control(self):
        """Function to control the addition and removal of rows in the table
        defining the momfile.
        The user can on runtime add or delete rows until a minimum of one row remains.
        New rows are created with dummy text in them.

        Author
        ----------
        Jonathan Chico
        """

        if self.sender() == self.InMomAddRow:
            rowPosition = self.InMomTable.rowCount()
            self.InMomTable.insertRow(rowPosition)
            text = [rowPosition + 1, 1, 1.0, 0.0, 0.0, 1.0]
            for ii, value in enumerate(text):
                item = QTableWidgetItem()
                item.setText(str(value))
                self.InMomTable.setItem(rowPosition, ii, item)
        if self.sender() == self.InMomDelRow:
            rowPosition = self.InMomTable.rowCount()
            # Make sure that one cannot delete the last entry
            if rowPosition > 1:
                self.InMomTable.removeRow(rowPosition - 1)
        return

    def CheckForFile(self, mainwindow):
        """If a momfile have already been selected, input it into the creator."""

        if len(mainwindow.ASDInputGen.momfile) > 0:
            momfile = np.genfromtxt(
                mainwindow.ASDInputGen.momfile.split("/")[-1], ndmin=2
            )
            momfile = [list(line) for line in momfile]
            self.momfile_gotten = True

            Table = self.InMomTable
            Table.setRowCount(0)
            for row, line in enumerate(momfile):
                Table.insertRow(row)
                for column, element in enumerate(line):
                    item = QTableWidgetItem(str(element))
                    Table.setItem(row, column, item)

    ############################################################################
    ## @brief Function handling the what the Cancel and Done buttons do in the \c momfile
    # window.
    # @details The Cancel button removes all the rows except for the first one and closes
    # the window.
    # The Done button reads the data in the cells and writes a momfile dubbed
    # momfile.dat in \c UppASD format.
    # @author Jonathan Chico
    ############################################################################
    def window_close(self):
        """Function handling the what the Cancel and Done buttons do in the momfile
        window.
        The Cancel button removes all the rows except for the first one and closes
        the window.
        The Done button reads the data in the cells and writes a momfile dubbed
        momfile.dat in UppASD format.

        Author
        ----------
        Jonathan Chico
        """

        if self.sender() == self.InpMomCancel:
            self.InMomTable.setRowCount(1)
            for column, value in enumerate([1, 1, 1.0, 0.0, 0.0, 1.0]):
                item = QTableWidgetItem(str(value))
                self.InMomTable.setItem(0, column, item)
            self.close()
        if self.sender() == self.InpMomDone:
            momfile_name = open(MomfileWindow.momfile_name, "w")
            for row in range(0, self.InMomTable.rowCount()):
                for col in range(0, self.InMomTable.columnCount()):
                    if col < 2:
                        entry = int(float(self.InMomTable.item(row, col).text()))
                    else:
                        entry = float(self.InMomTable.item(row, col).text())
                    momfile_name.write(f"{entry}  ")
                momfile_name.write("\n")
            self.close()
        self.momfile_gotten = True
        return


class JfileWindow(QDialog):
    """ "
    Class containing the defintions and actions needed for the display of the
    window handling the creation of the jfile inside the GUI. Class modified from
    MomfileWindow by Erik Karpelin.

    """

    def __init__(self, parent=None):

        super(JfileWindow, self).__init__(parent)
        path = os.path.dirname(os.path.abspath(__file__))
        uic.loadUi(os.path.join(path, "Jfile_Creator.ui"), self)
        self.InJfileAddRow.clicked.connect(self.table_control)
        self.InJfileDelRow.clicked.connect(self.table_control)
        self.InpJfileCancel.clicked.connect(self.window_close)
        self.InpJfileDone.clicked.connect(self.window_close)
        JfileWindow.jfile_gotten = False
        JfileWindow.jfile_name = "./jfile"
        return

    def table_control(self):
        """
        Function to control the addition and removal of rows in the table
        defining the jfile.
        The user can on runtime add or delete rows until a minimum of one row remains.
        New rows are created with dummy text in them.

        """

        if self.sender() == self.InJfileAddRow:
            rowPosition = self.InJfileTable.rowCount()
            self.InJfileTable.insertRow(rowPosition)
            text = [1, 1, 1.0, 0.0, 0.0, 0.0]
            for ii, value in enumerate(text):
                item = QTableWidgetItem()
                item.setText(str(value))
                self.InJfileTable.setItem(rowPosition, ii, item)
        if self.sender() == self.InJfileDelRow:
            rowPosition = self.InJfileTable.rowCount()
            # Make sure that one cannot delete the last entry
            if rowPosition > 1:
                self.InJfileTable.removeRow(rowPosition - 1)
        return

    def CheckForFile(self, mainwindow):
        """If a jfile have already been selected, input it into the creator."""

        if len(mainwindow.ASDInputGen.jfile) > 0:
            jfile = np.genfromtxt(mainwindow.ASDInputGen.jfile.split("/")[-1], ndmin=2)
            jfile = [list(line) for line in jfile]
            self.jfile_gotten = True

            Table = self.InJfileTable
            Table.setRowCount(0)
            for row, line in enumerate(jfile):
                Table.insertRow(row)
                for column, element in enumerate(line):
                    item = QTableWidgetItem(str(element))
                    Table.setItem(row, column, item)

    def GenerateVectorsFromCell(self, mainwindow):
        """
        Handles the generation of neighbour vector arrays and inputs
        them into the Jfile creation window.

        Input:
                mainwindow  :   QWindow object for the main UI window

        """

        Basis = np.array(
            [
                coord.text()
                for coord in mainwindow.findChildren(QLineEdit)
                if "InpLineEditC" in coord.objectName()
            ]
        ).reshape(3, 3)

        if "" in Basis:
            print("Input-error: Empty string encountered in cell")
            return

        Table = self.InJfileTable
        Table.setRowCount(0)
        CutoffRadius = int(self.InJfileNNCutoff.value())
        Positions, numbers = read_posfile(ASDInputgen.ASDInputGen.posfile)
        Cell = (np.float64(Basis), Positions, numbers)

        for i_site, site in enumerate(Positions):
            NeighbourVectors, NeighbourTypes, _ = get_full_nnlist(
                Cell, i_site, CutoffRadius, in_cell_only=False
            )

            VectorDict = {}
            CurrentSiteVector = np.ones((len(NeighbourVectors), 1)) * int(
                numbers[i_site]
            )
            KeyVector = np.hstack(
                (CurrentSiteVector, NeighbourTypes.T.reshape(CurrentSiteVector.shape))
            )

            for index, key in enumerate(KeyVector):
                key = " ".join([str(int(i)) for i in key])
                if key not in VectorDict:
                    VectorDict[key] = []
                VectorDict[key].append(list(NeighbourVectors[index]))

            if self.InJfileSymCheck.isChecked():
                for key in VectorDict:
                    VectorDict[key] = (
                        create_neighbour_list.reduce_vectors_from_symmetry(
                            Cell, np.array(VectorDict[key])
                        )
                    )

            self.InsertVectorsInTable(VectorDict, Table)

    def InsertVectorsInTable(self, VectorDict, Table):
        """
        Helper function to GenerateVectorsFromCell which inputs vectors
        into the file creation table

        Inputs:
                VectorDict  :   dictonary with interaction numbering as keys
                                and vectors as values
                Table       :   QTableWidget item for jfile table

        """

        for key in VectorDict:
            for vector in VectorDict[key]:
                Row = [key[0], key[-1], vector[0], vector[1], vector[2], 0]
                row = Table.rowCount()
                Table.insertRow(row)
                for column, value in enumerate(Row):
                    item = QTableWidgetItem(str(value))
                    Table.setItem(row, column, item)

    def window_close(self):
        """
        Function handling the what the Cancel and Done buttons do in the momfile
        window.
        The Cancel button removes all the rows except for the first, resets
        inputs for generation of vectors and closes the window.
        The Done button reads the data in the cells and writes a jfile dubbed
        jfile in UppASD format.
        Modified from momfile_window creation.

        """

        if self.sender() == self.InpJfileCancel:
            self.InJfileTable.setRowCount(1)
            for column, value in enumerate([1, 1, 1.0, 0.0, 0.0, 0.0]):
                item = QTableWidgetItem(str(value))
                self.InJfileTable.setItem(0, column, item)
            self.InJfileSymCheck.setChecked(False)
            self.InJfileNNCutoff.setValue(0)
            self.close()
        if self.sender() == self.InpJfileDone:
            jfile_name = open(JfileWindow.jfile_name, "w", encoding="utf-8")
            for row in range(0, self.InJfileTable.rowCount()):
                if self.InJfileTable.item(row, 5).text() == "0":
                    pass
                else:
                    for col in range(0, self.InJfileTable.columnCount()):
                        if col < 2:
                            entry = int(float(self.InJfileTable.item(row, col).text()))
                        else:
                            entry = float(self.InJfileTable.item(row, col).text())
                        jfile_name.write("{entry}  ".format(**locals()))
                    jfile_name.write("\n")
                self.close()
        self.jfile_gotten = True
        return


class DMfileWindow(QDialog):
    """ "
    Class containing the defintions and actions needed for the display of the
    window handling the creation of the jfile inside the GUI. Class modified from
    MomfileWindow by Erik Karpelin.

    """

    def __init__(self, parent=None):

        super(DMfileWindow, self).__init__(parent)
        path = os.path.dirname(os.path.abspath(__file__))
        uic.loadUi(os.path.join(path, "DMfile_Creator.ui"), self)
        self.InDMfileAddRow.clicked.connect(self.table_control)
        self.InDMfileDelRow.clicked.connect(self.table_control)
        self.InpDMfileCancel.clicked.connect(self.window_close)
        self.InpDMfileDone.clicked.connect(self.window_close)
        DMfileWindow.DMfile_gotten = False
        DMfileWindow.DMfile_name = "./dmfile"
        return

    def table_control(self):
        """
        Function to control the addition and removal of rows in the table
        defining the jfile.
        The user can on runtime add or delete rows until a minimum of one row remains.
        New rows are created with dummy text in them.

        """

        if self.sender() == self.InDMfileAddRow:
            rowPosition = self.InDMfileTable.rowCount()
            self.InDMfileTable.insertRow(rowPosition)
            text = [
                1,
                1,
                1.0,
                0.0,
                0.0,
                1.0,
                0.0,
                0.0,
            ]
            for ii, value in enumerate(text):
                item = QTableWidgetItem()
                item.setText(str(value))
                self.InDMfileTable.setItem(rowPosition, ii, item)
        if self.sender() == self.InDMfileDelRow:
            rowPosition = self.InDMfileTable.rowCount()
            # Make sure that one cannot delete the last entry
            if rowPosition > 1:
                self.InDMfileTable.removeRow(rowPosition - 1)
        return

    def CheckForFile(self, mainwindow):
        """If a jfile have already been selected, input it into the creator."""

        if len(mainwindow.ASDInputGen.dmfile) > 0:
            dmfile = np.genfromtxt(
                mainwindow.ASDInputGen.dmfile.split("/")[-1], ndmin=2
            )
            self.DMfile_gotten = True

            Table = self.InDMfileTable
            Table.setRowCount(0)
            for row, line in enumerate(dmfile):
                Table.insertRow(row)
                for column, element in enumerate(line):
                    item = QTableWidgetItem(str(element))
                    Table.setItem(row, column, item)

    def GenerateVectorsFromCell(self, mainwindow):
        """
        Handles the generation of neighbour vector arrays and inputs
        them into the DMfile creation window.

        Input:
                mainwindow  :   QWindow object for the main UI window

        """
        Basis = np.array(
            [
                coord.text()
                for coord in mainwindow.findChildren(QLineEdit)
                if "InpLineEditC" in coord.objectName()
            ]
        ).reshape(3, 3)

        if "" in Basis:
            print("Input-error: Empty string encountered in cell")
            return

        Table = self.InDMfileTable
        Table.setRowCount(0)
        CutoffRadius = int(self.InDMfileNNCutoff.value())
        Positions, numbers = read_posfile(ASDInputgen.ASDInputGen.posfile)
        Cell = (np.float64(Basis), Positions, numbers)

        for i_site, site in enumerate(Positions):
            NeighbourVectors, NeighbourTypes, _ = get_full_nnlist(
                Cell, i_site, CutoffRadius, in_cell_only=False
            )

            VectorDict = {}
            CurrentSiteVector = np.ones((len(NeighbourVectors), 1)) * int(
                numbers[i_site]
            )
            KeyVector = np.hstack(
                (CurrentSiteVector, NeighbourTypes.T.reshape(CurrentSiteVector.shape))
            )

            for index, key in enumerate(KeyVector):
                key = " ".join([str(int(i)) for i in key])
                if key not in VectorDict:
                    VectorDict[key] = []
                VectorDict[key].append(list(NeighbourVectors[index]))

            self.InsertVectorsInTable(VectorDict, Table)

    def InsertVectorsInTable(self, VectorDict, Table):
        """
        Helper function to GenerateVectorsFromCell which inputs vectors
        into the file creation table

        Inputs:
                VectorDict  :   dictonary with interaction numbering as keys
                                and vectors as values
                Table       :   QTableWidget
        """

        for key in VectorDict:
            for vector in VectorDict[key]:
                Row = [key[0], key[-1], vector[0], vector[1], vector[2], 0, 0, 0]
                row = Table.rowCount()
                Table.insertRow(row)
                for column, value in enumerate(Row):
                    item = QTableWidgetItem(str(value))
                    Table.setItem(row, column, item)

    def window_close(self):
        """
        Function handling the what the Cancel and Done buttons do in the file creation
        window.
        The Cancel button removes all the rows except for the first, resets
        inputs for generation of vectors and closes the window.
        The Done button reads the data in the cells and writes a DM-file.
        Modified from momfile_window creation.

        """

        if self.sender() == self.InpDMfileCancel:
            self.InDMfileTable.setRowCount(1)
            for column, value in enumerate([1, 1, 1.0, 0.0, 0.0, 1.0, 0.0, 0.0]):
                item = QTableWidgetItem(str(value))
                self.InDMfileTable.setItem(0, column, item)
            self.InDMfileNNCutoff.setValue(0)
            self.close()
        if self.sender() == self.InpDMfileDone:
            DMfile_name = open(DMfileWindow.DMfile_name, "w", encoding="utf-8")
            for row in range(0, self.InDMfileTable.rowCount()):
                DMVector = [
                    self.InDMfileTable.item(row, 5).text(),
                    self.InDMfileTable.item(row, 6).text(),
                    self.InDMfileTable.item(row, 7).text(),
                ]
                if DMVector == ["0", "0", "0"]:
                    pass
                else:
                    for col in range(0, self.InDMfileTable.columnCount()):
                        if col < 2:
                            entry = int(float(self.InDMfileTable.item(row, col).text()))
                        else:
                            entry = float(self.InDMfileTable.item(row, col).text())
                        DMfile_name.write(f"{entry}  ")
                    DMfile_name.write("\n")
                self.close()
        self.DMfile_gotten = True
        return


class KfileWindow(QDialog):
    """ "
    Class containing the defintions and actions needed for the display of the
    window handling the creation of the jfile inside the GUI. Class modified from
    MomfileWindow by Erik Karpelin.

    """

    def __init__(self, parent=None):

        super(KfileWindow, self).__init__(parent)
        path = os.path.dirname(os.path.abspath(__file__))
        uic.loadUi(os.path.join(path, "Kfile_Creator.ui"), self)
        self.InKfileAddRow.clicked.connect(self.table_control)
        self.InKfileDelRow.clicked.connect(self.table_control)
        self.InpKfileCancel.clicked.connect(self.window_close)
        self.InpKfileDone.clicked.connect(self.window_close)
        KfileWindow.Kfile_gotten = False
        KfileWindow.Kfile_name = "./kfile"
        return

    def table_control(self):
        """
        Function to control the addition and removal of rows in the table
        defining the jfile.
        The user can on runtime add or delete rows until a minimum of one row remains.
        New rows are created with dummy text in them.

        """

        if self.sender() == self.InKfileAddRow:
            rowPosition = self.InKfileTable.rowCount()
            self.InKfileTable.insertRow(rowPosition)
            text = [
                1,
                1,
                1.0,
                0.0,
                0.0,
                1.0,
                0.0,
                0.0,
            ]
            for ii, value in enumerate(text):
                item = QTableWidgetItem()
                item.setText(str(value))
                self.InKfileTable.setItem(rowPosition, ii, item)
        if self.sender() == self.InKfileDelRow:
            rowPosition = self.InKfileTable.rowCount()
            # Make sure that one cannot delete the last entry
            if rowPosition > 1:
                self.InKfileTable.removeRow(rowPosition - 1)
        return

    def CheckForFile(self, mainwindow):
        """If a jfile have already been selected, input it into the creator."""

        if len(mainwindow.ASDInputGen.kfile) > 0:
            kfile = np.genfromtxt(mainwindow.ASDInputGen.kfile.split("/")[-1], ndmin=2)
            self.Kfile_gotten = True

            Table = self.InKfileTable
            Table.setRowCount(0)
            for row, line in enumerate(kfile):
                Table.insertRow(row)
                for column, element in enumerate(line):
                    item = QTableWidgetItem(str(element))
                    Table.setItem(row, column, item)

    def window_close(self):
        """
        Function handling the what the Cancel and Done buttons do in the file creation
        window.
        The Cancel button removes all the rows except for the first, resets
        inputs for generation of vectors and closes the window.
        The Done button reads the data in the cells and writes a K-file.
        Modified from momfile_window creation.

        """

        if self.sender() == self.InpKfileCancel:
            self.InKfileTable.setRowCount(1)
            for column, value in enumerate([1, 1, 1.0, 0.0, 0.0, 1.0, 0.0, 0.0]):
                item = QTableWidgetItem(str(value))
                self.InKfileTable.setItem(0, column, item)
            self.InKfileNNCutoff.setValue(0)
            self.close()
        if self.sender() == self.InpKfileDone:
            Kfile_name = open(KfileWindow.Kfile_name, "w", encoding="utf-8")
            for row in range(0, self.InKfileTable.rowCount()):
                KVector = [
                    self.InKfileTable.item(row, 5).text(),
                    self.InKfileTable.item(row, 6).text(),
                    self.InKfileTable.item(row, 7).text(),
                ]
                if KVector == ["0", "0", "0"]:
                    pass
                else:
                    for col in range(0, self.InKfileTable.columnCount()):
                        if col < 2:
                            entry = int(float(self.InKfileTable.item(row, col).text()))
                        else:
                            entry = float(self.InKfileTable.item(row, col).text())
                        Kfile_name.write(f"{entry}  ")
                    Kfile_name.write("\n")
                self.close()
        self.Kfile_gotten = True
        return
