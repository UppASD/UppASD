"""@package ASDInputGen
Class containing the main defintions for the inpsd.dat writer.

This class has the information needed to generate the dictionary used to print the
inpsd.dat, collect the data from the GUI and print it to file, both in inpsd.dat
format and in .yaml format.

Author
----------
Jonathan Chico
"""
################################################################################
# @brief Class containing the structures needed for the creation of the inpsd.dat
# via the GUI.
# @details It defines all the variables set in the GUI as part of a dictionary
# which is populated by certain default values. These are then purged to ensure
# that one keeps only the minimal set of input variables needed.
# @author Jonathan Chico
################################################################################

class ASDInputGen():
    """Class containing the structures needed for the creation of the inpsd.dat
    via the GUI. It defines all the variables set in the GUI as part of a dictionary
    which is populated by certain default values. These are then purged to ensure
    that one keeps only the minimal set of input variables needed.

    Author
    ----------
    Jonathan Chico
    """

    def __init__(self):
        import collections
        from PyQt6.QtGui import QIntValidator, QDoubleValidator
        ASDInputGen.UppASDKeywords = collections.OrderedDict()
        ASDInputGen.IntegerValidator = QIntValidator()
        ASDInputGen.IntegerValidator.setRange(0, 99999999)
        ASDInputGen.PosDoubleValidator = QDoubleValidator()
        ASDInputGen.PosDoubleValidator.setRange(0, 99999999.9999)
        ASDInputGen.PosDoubleValidator.setDecimals(10)
        ASDInputGen.DoubleValidator = QDoubleValidator()
        ASDInputGen.DoubleValidator.setDecimals(10)
        ASDInputGen.momfile = []
        ASDInputGen.momfile_in = []
        ASDInputGen.momfile_fi = []
        ASDInputGen.posfile = []
        ASDInputGen.jfile = []
        ASDInputGen.dmfile = []
        ASDInputGen.kfile = []
        ASDInputGen.restartfile = []
        ASDInputGen.pdfile = []
        ASDInputGen.bqfile = []
        ASDInputGen.bqdmfile = []
        ASDInputGen.qfile = []
        ASDInputGen.momfile_gotten = False
        ASDInputGen.posfile_gotten = False
        return
    ############################################################################
    # @brief Function to find the needed file names for the input file created.
    # @author Jonathan Chico
    ############################################################################

    def getFileName(self, window):
        from PyQt6 import QtWidgets
        import os

        # Check implemented to prevent window to open on unchecking the button. 
        if window.sender().__class__.__name__ != "QPushButton" and not window.sender().isChecked():
            pass

        else:
            dlg = QtWidgets.QFileDialog()
            dlg.setFileMode(QtWidgets.QFileDialog.FileMode.AnyFile)
            dlg.setDirectory('.')
            if dlg.exec():
                filename = './' + os.path.split(dlg.selectedFiles()[0])[-1]
                if window.sender() == window.InpPosButtonSelect:
                    ASDInputGen.posfile_gotten = False
                    ASDInputGen.posfile = filename
                    ASDInputGen.posfile_gotten = True
                if window.sender() == window.InpMomButtonSelect:
                    ASDInputGen.momfile_gotten = False
                    ASDInputGen.momfile = filename
                    ASDInputGen.momfile_gotten = True
                if window.sender() == window.InpInitMag4ReadButton:
                    ASDInputGen.restartfile = filename
                if window.sender() == window.InpXCCheck:
                    ASDInputGen.jfile = filename
                if window.sender() == window.InpDMCheck and window.InpDMCheck.isChecked():
                    ASDInputGen.dmfile = filename
                if window.sender() == window.InpMAECheck and window.InpMAECheck.isChecked():
                    ASDInputGen.kfile = filename
                if window.sender() == window.InpPseudoCheck and window.InpPseudoCheck.isChecked():
                    ASDInputGen.pdfile = filename
                if window.sender() == window.InpBqCheck and window.InpBqCheck.isChecked():
                    ASDInputGen.bqfile = filename
                if window.sender() == window.InpBqDMCheck and window.InpBqDMCheck.isChecked():
                    ASDInputGen.bqdmfile = filename
                if window.sender() == window.InpSetIniMomfileButton:
                    ASDInputGen.momfile_in = filename
                if window.sender() == window.InpSetFinMomfileButton:
                    ASDInputGen.momfile_fi = filename
                if window.sender() == window.InpSqQpoints:
                    ASDInputGen.qfile = filename
                if window.sender() == window.InpJfileButtonSelect: 
                    window.InpXCCheck.setChecked(True)
                    ASDInputGen.jfile = filename
                if window.sender() == window.InpKfileButtonSelect:
                    ASDInputGen.kfile = filename
                    window.InpMAECheck.setChecked(True)
                if window.sender() == window.InpDMButtonSelect: 
                    window.InpDMCheck.setChecked(True)
                    ASDInputGen.dmfile = filename
                if window.sender() == window.InpPseudDipSelect:
                    ASDInputGen.pdfile = filename
                    window.InpPseudoCheck.setChecked(True)
                if window.sender() == window.InpBiQSelect:
                    ASDInputGen.bqfile = filename
                    window.InpBqCheck.setChecked(True)
                if window.sender() == window.InpBiQDMSelect:
                    ASDInputGen.bqdmfile = filename
                    window.InpBqDMCheck.setChecked(True)

            elif window.sender().__class__.__name__  == "QCheckBox" and window.sender().isChecked():
                window.sender().setChecked(False)

        return
    ############################################################################
    # @brief Updating the file names in the input data when files are generated by the GUI.
    # @author Jonathan Chico
    ############################################################################

    def update_file_name(self, window):
        if window.sender() == window.Posfile_Window.InpPosDone:
            ASDInputGen.posfile = window.Posfile_Window.posfile_name
            self.posfile_gotten = True
        if window.sender() == window.Momfile_Window.InpMomDone:
            ASDInputGen.momfile = window.Momfile_Window.momfile_name
            self.momfile_gotten = True
        if window.sender() == window.Restart_Window.InpRestartDone:
            ASDInputGen.restartfile = window.Restart_Window.restartfile_name
        if window.sender() == window.Jfile_Window.InpJfileDone:
            ASDInputGen.jfile = window.Jfile_Window.jfile_name
            window.InpXCCheck.setChecked(True)
        if window.sender() == window.DMfile_Window.InpDMfileDone:
            ASDInputGen.dmfile = window.DMfile_Window.DMfile_name
            window.InpDMCheck.setChecked(True)
        return
    ############################################################################
    # @brief Function to get the information needed for the inpsd.dat from the
    # GUI to pass it to the dictionary.
    # @author Jonathan Chico
    ############################################################################

    def ASDInputGatherer(self, window):

        # Find the simulation name
        if len(window.InpLineEditSimid.text()) > 0:
            ASDInputGen.UppASDKeywords['general']['simid'] =\
                str(window.InpLineEditSimid.text())
           
        # Obtain the lattice constant of the system
        if len(window.InpLineEditAlat.text()) > 0:
            ASDInputGen.UppASDKeywords['geometry']['alat'] =\
                float(self.text_to_num(window.InpLineEditAlat.text()))
           
        # The number of repetitions of the unit cell
        if len(window.InpN1.text()) > 0:
            ASDInputGen.UppASDKeywords['geometry']['ncell'][0] =\
                int(self.text_to_num(window.InpN1.text()))
        if len(window.InpN2.text()) > 0:
            ASDInputGen.UppASDKeywords['geometry']['ncell'][1] =\
                int(self.text_to_num(window.InpN2.text()))
        if len(window.InpN3.text()) > 0:
            ASDInputGen.UppASDKeywords['geometry']['ncell'][2] =\
                int(self.text_to_num(window.InpN3.text()))
           
        # Set the boundary conditions
        if window.InpPBCCheckC1.isChecked():
            ASDInputGen.UppASDKeywords['geometry']['BC'][0] = 'P'
        else:
            ASDInputGen.UppASDKeywords['geometry']['BC'][0] = 0
        if window.InpPBCCheckC2.isChecked():
            ASDInputGen.UppASDKeywords['geometry']['BC'][1] = 'P'
        else:
            ASDInputGen.UppASDKeywords['geometry']['BC'][1] = 0
        if window.InpPBCCheckC3.isChecked():
            ASDInputGen.UppASDKeywords['geometry']['BC'][2] = 'P'
        else:
            ASDInputGen.UppASDKeywords['geometry']['BC'][2] = 0

        # Filling up the lattice vectors
        UiBasisVectors = [[window.InpLineEditC1_x, window.InpLineEditC1_y, window.InpLineEditC1_z],
                            [window.InpLineEditC2_x, window.InpLineEditC2_y, window.InpLineEditC2_z],
                            [window.InpLineEditC3_x, window.InpLineEditC3_y, window.InpLineEditC3_z]]
        for i, vector in enumerate(UiBasisVectors):
            for j, coord in enumerate(vector):
                if len(coord.text()) > 0:
                    ASDInputGen.UppASDKeywords['geometry']['cell'][i][j] =\
                    float(coord.text())
       
        # # Check for random alloys
        # if window.InpCheckRandAlloy.isChecked():
        #     ASDInputGen.UppASDKeywords['geometry']['do_ralloy'] = 1
        # else:
        #     ASDInputGen.UppASDKeywords['geometry']['do_ralloy'] = 0

        # Check for Hamiltonian interactions
        if window.InpXCCheck.isChecked():
            ASDInputGen.UppASDKeywords['general']['do_prnstruct'] = 1
            if len(ASDInputGen.jfile) > 0:
                ASDInputGen.UppASDKeywords['Hamiltonian']['exchange'] = ASDInputGen.jfile
            else:
                print('No exchange file name given, assuming "./jfile"')
        if window.InpDMCheck.isChecked():
            ASDInputGen.UppASDKeywords['Hamiltonian']['dm'] = ASDInputGen.dmfile
            ASDInputGen.UppASDKeywords['Hamiltonian']['do_dm'] = 1
        else:
            ASDInputGen.UppASDKeywords['Hamiltonian']['do_dm'] = 0
        if window.InpMAECheck.isChecked():
            ASDInputGen.UppASDKeywords['Hamiltonian']['anisotropy'] = ASDInputGen.kfile
            ASDInputGen.UppASDKeywords['Hamiltonian']['do_anisotropy'] = 1
        else:
            ASDInputGen.UppASDKeywords['Hamiltonian']['do_anisotropy'] = 0
        if window.InpPseudoCheck.isChecked():
            ASDInputGen.UppASDKeywords['Hamiltonian']['pd'] = ASDInputGen.pdfile
            ASDInputGen.UppASDKeywords['Hamiltonian']['do_pd'] = 1
        else:
            ASDInputGen.UppASDKeywords['Hamiltonian']['do_pd'] = 0
        if window.InpBqCheck.isChecked():
            ASDInputGen.UppASDKeywords['Hamiltonian']['bq'] = ASDInputGen.bqfile
            ASDInputGen.UppASDKeywords['Hamiltonian']['do_bq'] = 1
        else:
            ASDInputGen.UppASDKeywords['Hamiltonian']['do_bq'] = 0
        if window.InpBqCheck.isChecked():
            ASDInputGen.UppASDKeywords['Hamiltonian']['biqdm'] = ASDInputGen.bqdmfile
            ASDInputGen.UppASDKeywords['Hamiltonian']['do_biqdm'] = 1
        else:
            ASDInputGen.UppASDKeywords['Hamiltonian']['do_biqdm'] = 0

        # Filling up the magnetic field in the measurement phase
        ASDInputGen.UppASDKeywords['Hamiltonian']['hfield'] =\
            [float(self.text_to_num(window.InpBextMeasure_x.text())),
             float(self.text_to_num(window.InpBextMeasure_y.text())),
             float(self.text_to_num(window.InpBextMeasure_z.text()))]
        
        # Filling up the magnetic field in the initial phase
        ASDInputGen.UppASDKeywords['Hamiltonian']['ip_hfield'] =\
            [float(self.text_to_num(window.InpBextInit_x.text())),
             float(self.text_to_num(window.InpBextInit_y.text())),
             float(self.text_to_num(window.InpBextInit_z.text()))]
        
        # Damping for an SD measurement phase
        if len(window.InpASDLLGDamp.text()) > 0:
            ASDInputGen.UppASDKeywords['LLG_mphase']['damping'] =\
                float(self.text_to_num(window.InpASDLLGDamp.text()))
            
        # Time step for an SD measurement phase
        if len(window.InpASDLLGDT.text()) > 0:
            ASDInputGen.UppASDKeywords['LLG_mphase']['timestep'] =\
                float(self.text_to_num(window.InpASDLLGDT.text()))
            
        # Number of steps for an SD measurement phase
        if len(window.InpASDLLGSteps.text()) > 0:
            ASDInputGen.UppASDKeywords['LLG_mphase']['Nstep'] =\
                int(self.text_to_num(window.InpASDLLGSteps.text()))
        
        # Number of MC steps in the measurement phase
        if len(window.InpMCSteps.text()) > 0:
            ASDInputGen.UppASDKeywords['MC_mphase']['mcnstep'] =\
                int(self.text_to_num(window.InpMCSteps.text()))
        
        # Number of ensembles in the simulation
        if len(window.InpMensemble.text()) > 0:
            ASDInputGen.UppASDKeywords['Mag']['Mensemble'] =\
                int(self.text_to_num(window.InpMensemble.text()))
            
        # Check if the trajectories will be printed
        if window.InpGlobalTrajBox.isChecked():
            ASDInputGen.UppASDKeywords['trajectories']['do_tottraj'] = 'Y'
        else:
            ASDInputGen.UppASDKeywords['trajectories']['do_tottraj'] = 'N'

        # Information about the printing of the trajectories
        if len(window.InitTTrajStep.text()) > 0:
            ASDInputGen.UppASDKeywords['trajectories']['tottraj_step'] =\
                int(self.text_to_num(window.InitTTrajStep.text()))
        if len(window.InitTTrajBuff.text()) > 0:
            ASDInputGen.UppASDKeywords['trajectories']['tottraj_buff'] =\
                int(self.text_to_num(window.InitTTrajBuff.text()))
            
        # Check if the cumulants will be printed
        if window.InpCumuBox.isChecked():
            ASDInputGen.UppASDKeywords['cumulants']['do_cumu'] = 'Y'
        else:
            ASDInputGen.UppASDKeywords['cumulants']['do_cumu'] = 'N'

        # Options for the printing of the cumulants
        if len(window.InpCumuStep.text()) > 0:
            ASDInputGen.UppASDKeywords['cumulants']['cumu_step'] =\
                int(self.text_to_num(window.InpCumuStep.text()))
        if len(window.InpCumuBuff.text()) > 0:
            ASDInputGen.UppASDKeywords['cumulants']['cumu_buff'] =\
                int(self.text_to_num(window.InpCumuBuff.text()))
            
        # Check if the averages should be printed
        if window.InpAveBox.isChecked():
            ASDInputGen.UppASDKeywords['averages']['do_avrg'] = 'Y'
        if not window.InpAveBox.isChecked():
            ASDInputGen.UppASDKeywords['averages']['do_avrg'] = 'N'
        
        # Options for the printing of the averages
        if len(window.InpAveStep.text()) > 0:
            ASDInputGen.UppASDKeywords['averages']['avrg_step'] =\
                int(self.text_to_num(window.InpAveStep.text()))
        if len(window.InpAveBuff.text()) > 0:
            ASDInputGen.UppASDKeywords['averages']['avrg_buff'] =\
                int(self.text_to_num(window.InpAveBuff.text()))
        
        # Check if the site projected average will be printed
        if window.InpSiteAveCheck.isChecked():
            ASDInputGen.UppASDKeywords['averages']['do_proj_avrg'] = 'A'
        
        # Check if the type projected average will be printed
        elif window.InpTypeAveCheck.isChecked():
            ASDInputGen.UppASDKeywords['averages']['do_proj_avrg'] = 'Y'
        else:
            ASDInputGen.UppASDKeywords['averages']['do_proj_avrg'] = 'N'
        
        # Check if the chemically projected average will be printed
        if window.InpChemAveCheck.isChecked():
            ASDInputGen.UppASDKeywords['averages']['do_projch_avrg'] = 'Y'
        else:
            ASDInputGen.UppASDKeywords['averages']['do_projch_avrg'] = 'N'
        
        # Check if the dynamical structure factor will be calculated
        if window.ImpSqwBox.isChecked():
            ASDInputGen.UppASDKeywords['Mag_corr']['do_sc'] = 'Q'
        if not window.ImpSqwBox.isChecked():
            ASDInputGen.UppASDKeywords['Mag_corr']['do_sc'] = 'N'
        
        # Check if the static correlation is calculated
        if window.InpScCheck.isChecked():
            ASDInputGen.UppASDKeywords['Mag_corr']['do_sc'] = 'C'
        
        # Options for the structure factor
        if len(window.InpScStep.text()) > 0:
            ASDInputGen.UppASDKeywords['Mag_corr']['sc_step'] =\
                int(self.text_to_num(window.InpScStep.text()))
        if len(window.InpScNStep.text()) > 0:
            ASDInputGen.UppASDKeywords['Mag_corr']['sc_nstep'] =\
                int(self.text_to_num(window.InpScNStep.text()))
        
        # Check if the AMS will be printed
        if window.InpAMSCheck.isChecked():
            ASDInputGen.UppASDKeywords['Mag_corr']['do_ams'] = 'Y'
            if len(ASDInputGen.qfile) == 0:
                import ASD_GUI.Extras.preQ as preQ
                Positions, numbers = preQ.read_posfile(ASDInputGen.posfile)
                Cell = (ASDInputGen.UppASDKeywords['geometry']['cell'], Positions, numbers)
                preQ.Run(Cell, ASDInputGen.UppASDKeywords['general']['simid'], ASDInputGen.UppASDKeywords['geometry']['ncell'], 'C')
                ASDInputGen.UppASDKeywords['Mag_corr']['qfile'] = 'qfile.kpath'

        else:
            ASDInputGen.UppASDKeywords['Mag_corr']['do_ams'] = 'N'
        
        # Check if the STT is going to be considered
        if window.InpSTTBox.isChecked():
            if window.InpZhangLiCheck.isChecked():
                ASDInputGen.UppASDKeywords['spintorque']['stt'] = 'A'
            if window.InpSlonwceskiCheck.isChecked():
                ASDInputGen.UppASDKeywords['spintorque']['stt'] = 'Y'
            # Check if the SHE is considered
            if window.InpSHEBox.isChecked():
                ASDInputGen.UppASDKeywords['spintorque']['do_she'] = 'Y'
            else:
                ASDInputGen.UppASDKeywords['spintorque']['do_she'] = 'N'
        else:
            ASDInputGen.UppASDKeywords['spintorque']['stt'] = 'N'
       
        # Non-adiabatic parameter for the STT
        ASDInputGen.UppASDKeywords['spintorque']['adibeta'] =\
            float(self.text_to_num(window.InpBeta.text()))
       
        # Fill up the spin current vector for the STT
        ASDInputGen.UppASDKeywords['spintorque']['jvec'] =\
            [float(self.text_to_num(window.InpJvec_x.text())),
             float(self.text_to_num(window.InpJvec_y.text())),
             float(self.text_to_num(window.InpJvec_z.text()))]
       
        # Parameters for the SHE
        ASDInputGen.UppASDKeywords['spintorque']['thick_ferro'] =\
            float(self.text_to_num(window.InpSheThickness.text()))
        ASDInputGen.UppASDKeywords['spintorque']['she_angle'] =\
            float(self.text_to_num(window.InpShe.text()))
       
        # Check if SOT is going to be considered
        if window.InpSOTBox.isChecked():
            ASDInputGen.UppASDKeywords['spintorque']['do_sot'] = 'Y'
        else:
            ASDInputGen.UppASDKeywords['spintorque']['do_sot'] = 'N'
       
        # Parameters for the SOT
        ASDInputGen.UppASDKeywords['spintorque']['sot_field'] =\
            float(self.text_to_num(window.InpSOTFL.text()))
        ASDInputGen.UppASDKeywords['spintorque']['sot_damping'] =\
            float(self.text_to_num(window.InpSOTDL.text()))
       
        # Fill up the spin polarization vector for the SOT
        ASDInputGen.UppASDKeywords['spintorque']['sot_pol_vec'] =\
            [float(self.text_to_num(window.InpSOTPol_x.text())),
             float(self.text_to_num(window.InpSOTPol_y.text())),
             float(self.text_to_num(window.InpSOTPol_z.text()))]
        ASDInputGen.UppASDKeywords['LLG_mphase']['SDEAlgh'] =\
            int(window.InpASDLLGAlgh.value())
       
        # Set the symmetry of the exchange Hamiltonian
        ASDInputGen.UppASDKeywords['Hamiltonian']['Sym'] =\
            int(window.InpPairSym.value())
       
        # Set the maptype of the pairwise Hamiltonian
        ASDInputGen.UppASDKeywords['Hamiltonian']['maptype'] =\
            int(window.ImpPairMaptype.value())
       
        # Check if the dipole-dipole interaction should be considered
        if window.InpDipBruteForceCheck.isChecked():
            ASDInputGen.UppASDKeywords['Hamiltonian']['do_dip'] = 1
        if window.InpDipMacroCheck.isChecked():
            ASDInputGen.UppASDKeywords['Hamiltonian']['do_dip'] = 2
        if window.InpDipFFTCheck.isChecked():
            ASDInputGen.UppASDKeywords['Hamiltonian']['do_dip'] = 3
        if not window.InpDipBox.isChecked():
            ASDInputGen.UppASDKeywords['Hamiltonian']['do_dip'] = 0
       
        # Check which mode is going to be used for the measurement phase
        if window.InpMeasureLLG.isChecked():
            ASDInputGen.UppASDKeywords['general']['mode'] = 'S'
        if window.InpMeasureMCMet.isChecked():
            ASDInputGen.UppASDKeywords['general']['mode'] = 'M'
        if window.InpMeasureMCHeat.isChecked():
            ASDInputGen.UppASDKeywords['general']['mode'] = 'H'
        if window.InpMeasureGNEB.isChecked():
            ASDInputGen.UppASDKeywords['general']['mode'] = 'G'
            if window.InpGNEBMeasureBox.isChecked():
                ASDInputGen.UppASDKeywords['GNEB_mphase']['do_gneb'] = 'Y'
                ASDInputGen.UppASDKeywords['GNEB_mphase']['mep_ftol'] =              \
                    float(self.text_to_num(window.InpGNEBMEPTol.text()))
                ASDInputGen.UppASDKeywords['GNEB_mphase']['mep_itrmax'] =            \
                    int(self.text_to_num(window.InpGNEBMEPSteps.text()))
                ASDInputGen.UppASDKeywords['GNEB_mphase']['eig_zero'] =              \
                    float(self.text_to_num(window.InpMinEig.text()))
            else:
                ASDInputGen.UppASDKeywords['GNEB_mphase']['do_gneb'] = 'N'
            if window.InpMeasureGNEBCI.isChecked():
                ASDInputGen.UppASDKeywords['GNEB_mphase']['do_gneb_ci'] = 'Y'
                ASDInputGen.UppASDKeywords['GNEB_mphase']['mep_ftol_ci'] =           \
                    float(self.text_to_num(window.InpGNEBCITol.text()))
            else:
                ASDInputGen.UppASDKeywords['GNEB_mphase']['do_gneb_ci'] = 'N'
       
        # Check which is the initial magnetization of the system
        if window.InpInitmag1Check.isChecked():
            ASDInputGen.UppASDKeywords['Mag']['initmag'] = 1
        if window.InpInitmag2Box.isChecked():
            ASDInputGen.UppASDKeywords['Mag']['initmag'] = 2
        if window.InpInitmag3Check.isChecked():
            ASDInputGen.UppASDKeywords['Mag']['initmag'] = 3
        if window.InpInitMag4Check.isChecked():
            ASDInputGen.UppASDKeywords['Mag']['initmag'] = 4
            if len(ASDInputGen.restartfile) > 0:
                ASDInputGen.UppASDKeywords['Mag']['restartfile'] = ASDInputGen.restartfile
            else:
                print('No restartfile name given assuming "./restart.dummy.dat"')
                ASDInputGen.UppASDKeywords['Mag']['restartfile'] = './restart.dummy.dat'
        if window.InpInitmag7Check.isChecked():
            ASDInputGen.UppASDKeywords['Mag']['initmag'] = 7
            if len(ASDInputGen.restartfile) > 0:
                ASDInputGen.UppASDKeywords['Mag']['restartfile'] = ASDInputGen.restartfile
            else:
                print('No restartfile name given assuming "./restart.dummy.dat"')
                ASDInputGen.UppASDKeywords['Mag']['restartfile'] = './restart.dummy.dat'
        if window.InpInitMag6Check.isChecked():
            ASDInputGen.UppASDKeywords['Mag']['initmag'] = 6
            if len(ASDInputGen.momfile_in) > 0:
                ASDInputGen.UppASDKeywords['Mag']['momfile_i'] = ASDInputGen.momfile_in
            else:
                print('No momfile_i name given assuming "./momfile_i.dat"')
                ASDInputGen.UppASDKeywords['Mag']['momfile_i'] = './momfile_i.dat'
            if len(ASDInputGen.momfile_fi) > 0:
                ASDInputGen.UppASDKeywords['Mag']['momfile_f'] = ASDInputGen.momfile_fi
            else:
                print('No momfile_f name given assuming "./momfile_f.dat"')
                ASDInputGen.UppASDKeywords['Mag']['momfile_f'] = './momfile_f.dat'
        if window.InpIniFinCheck.isChecked():
            ASDInputGen.UppASDKeywords['Mag']['initpath'] = 1
        if window.InpFullPathChek.isChecked():
            ASDInputGen.UppASDKeywords['Mag']['initpath'] = 2
        if window.InpRelGNEBCheck.isChecked():
            ASDInputGen.UppASDKeywords['Mag']['relaxed_if'] = 'Y'
        else:
            ASDInputGen.UppASDKeywords['Mag']['relaxed_if'] = 'N'
       
        # Check which is the initial phase of the calculation
        if window.InpInitBox.isChecked():
            if window.InpInitLLG.isChecked():
                ASDInputGen.UppASDKeywords['general']['ip_mode'] = 'S'
                ASDInputGen.UppASDKeywords['LLG_iphase']['ip_nphase'] =\
                    window.init_phase_data[0]
                ASDInputGen.UppASDKeywords['LLG_iphase'][''] = window.init_phase_data[1:]
            if window.InpInitMcMet.isChecked():
                ASDInputGen.UppASDKeywords['general']['ip_mode'] = 'M'
                ASDInputGen.UppASDKeywords['MC_iphase']['ip_mcanneal'] =\
                    window.init_phase_data[0]
                ASDInputGen.UppASDKeywords['MC_iphase'][''] = window.init_phase_data[1:]
            if window.InpInitMcHeat.isChecked():
                ASDInputGen.UppASDKeywords['general']['ip_mode'] = 'H'
                ASDInputGen.UppASDKeywords['MC_iphase']['ip_mcanneal'] =\
                    window.init_phase_data[0]
                ASDInputGen.UppASDKeywords['MC_iphase'][''] = window.init_phase_data[1:]
            if window.InpInitVPO.isChecked():
                ASDInputGen.UppASDKeywords['general']['ip_mode'] = 'G'
                ASDInputGen.UppASDKeywords['VPO_iphase']['min_itrmax'] = window.init_phase_data[0][0]
                ASDInputGen.UppASDKeywords['VPO_iphase']['spring'] = window.init_phase_data[0][1]
                ASDInputGen.UppASDKeywords['VPO_iphase']['vpo_mass'] = window.init_phase_data[0][2]
                ASDInputGen.UppASDKeywords['VPO_iphase']['vpo_dt'] = window.init_phase_data[0][3]
                ASDInputGen.UppASDKeywords['VPO_iphase']['min_ftol'] = window.init_phase_data[0][4]
        else:
            ASDInputGen.UppASDKeywords['general']['ip_mode'] = 'N'
       
        # Check if the energy will be printed
        if window.InpEneTotalEneCheck.isChecked():
            ASDInputGen.UppASDKeywords['energy']['plotenergy'] = 1
        if window.InpEneSiteEneCheck.isChecked():
            ASDInputGen.UppASDKeywords['energy']['plotenergy'] = 2
        if not window.InpEneTotalEneCheck.isChecked() and not window.InpEneSiteEneCheck.isChecked():
            ASDInputGen.UppASDKeywords['energy']['plotenergy'] = 0
       
        # Check if the skyrmion number should be printed
        if window.SkxNumBox.isChecked():
            ASDInputGen.UppASDKeywords['topology']['skyrno'] = 'Y'
        else:
            ASDInputGen.UppASDKeywords['topology']['skyrno'] = 'N'
       
        # Check for the Hessian options
        if window.HessFinCheck.isChecked():
            ASDInputGen.UppASDKeywords['Hessians']['do_hess_fin'] = 'Y'
        else:
            ASDInputGen.UppASDKeywords['Hessians']['do_hess_fin'] = 'N'
        if window.HessInitCheck.isChecked():
            ASDInputGen.UppASDKeywords['Hessians']['do_hess_ini'] = 'Y'
        else:
            ASDInputGen.UppASDKeywords['Hessians']['do_hess_ini'] = 'N'
        if window.HessSPCheck.isChecked():
            ASDInputGen.UppASDKeywords['Hessians']['do_hess_sp'] = 'Y'
        else:
            ASDInputGen.UppASDKeywords['Hessians']['do_hess_sp'] = 'N'
       
        # Options for the printing of the skyrmion number
        if len(window.InpSkxStep.text()) > 0:
            ASDInputGen.UppASDKeywords['topology']['skyno_step'] =\
                int(self.text_to_num(window.InpSkxStep.text()))
        if len(window.InitSkxBuff.text()) > 0:
            ASDInputGen.UppASDKeywords['topology']['skyno_buff'] =\
                int(self.text_to_num(window.InitSkxBuff.text()))
        if len(window.InpDipBlockSizeLineEdit.text()) > 0:
            ASDInputGen.UppASDKeywords['Hamiltonian']['block_size'] =\
                int(self.text_to_num(window.InpDipBlockSizeLineEdit.text()))
        else:
            ASDInputGen.UppASDKeywords['Hamiltonian']['block_size'] = 1
      
        # Check temperature of system
        if len(window.InpASDLLGTemp.text()) > 0:
            ASDInputGen.UppASDKeywords['general']['Temp'] = \
                int(self.text_to_num(window.InpASDLLGTemp.text()))
          
        if len(window.InpMCTemp.text()) > 0:
            ASDInputGen.UppASDKeywords['general']['Temp'] = \
                int(self.text_to_num(window.InpMCTemp.text()))

        # Quick relaxation
        if window.InpQuickRelax.isChecked():
            ASDInputGen.UppASDKeywords['general']['ip_mode'] = 'SX'
            ASDInputGen.UppASDKeywords['MC_iphase']['ip_mcanneal'] = 4
            ASDInputGen.UppASDKeywords['MC_iphase'][''] =\
                  [[1000, 0.0001], [1000, 0.0001],[1000, 0.0001],[1000, 0.0001]]

        # Check mom/pos-file names
        if len(ASDInputGen.posfile) > 0:
            ASDInputGen.UppASDKeywords['geometry']['posfile'] = ASDInputGen.posfile
        else:
            print('No posfile name given, assuming "./posfile"')
        if len(ASDInputGen.momfile) > 0:
            ASDInputGen.UppASDKeywords['geometry']['momfile'] = ASDInputGen.momfile
        elif len(ASDInputGen.posfile) != 0:
            print('No momfile name given, generating "momfile.dummy"')
            ASDInputGen.UppASDKeywords['geometry']['momfile'] = ['momfile.dummy']
            vector = [1, 0, 0]
            self.GenDummyMomfile(vector, 1)
        else:
            print('No momfile name given, assuming "./momfile"')

        return
    
    ############################################################################

    def SetStructureTemplate(self, window, structure):
        """
        Function handling the structure templates by writing the disired
        structure to the UI.

        Inputs:
                window      :   QMainWindow
                structure   :   string passed from button press
        
        Author: Erik Karpelin
        """
        import numpy as np

        ASDInputGen.posfile_gotten = True
        ASDInputGen.posfile = './posfile'

        # Set symmetry of lattice
        if structure == 'sc' or 'bcc' or 'bcc2' or 'fcc':
            window.InpPairSym.setValue(1)
        if structure == 'hcp':
            window.InpPairSym.setValue(3)

        UiBasisVectors = np.array([[window.InpLineEditC1_x, window.InpLineEditC1_y, window.InpLineEditC1_z],
                                [window.InpLineEditC2_x, window.InpLineEditC2_y, window.InpLineEditC2_z],
                                [window.InpLineEditC3_x, window.InpLineEditC3_y, window.InpLineEditC3_z]])

        SimpleCubic = BodyCentric = np.array([['1.000', '0.000', '0.000'],
                                            ['0.000', '1.000', '0.000'],
                                            ['0.000', '0.000', '1.000']])
        
        FaceCentered = np.array([['0.000', '0.500', '0.500'],
                                ['0.500', '0.000', '0.500'],
                                ['0.500', '0.500', '0.000']])
        
        HCPStructure = np.array([['1.000', '0.000', '0.000'],
                                ['-0.500', '0.8660254037', '0.000'],
                                ['0.000', '0.000', '1.6329931618']])
        

        self.ClearBasisVectors(UiBasisVectors)
        
        if structure == 'sc':
            self.InsertBasisvectors(UiBasisVectors, SimpleCubic)
            Positions = '1 1    0.000   0.000   0.000'
        
        if structure == 'bcc':
            self.InsertBasisvectors(UiBasisVectors, BodyCentric)
            Positions = '1 1    0.000   0.000   0.000\n'\
                        '2 1    0.500   0.500   0.500'
            
        if structure == 'bcc2':
            self.InsertBasisvectors(UiBasisVectors, BodyCentric)
            Positions = '1 1    0.000   0.000   0.000\n'\
                        '2 2    0.500   0.500   0.500'
            
        if structure == 'fcc':
            self.InsertBasisvectors(UiBasisVectors, FaceCentered)
            Positions = '1 1    0.000   0.000   0.000'

        if structure == 'hcp':
            self.InsertBasisvectors(UiBasisVectors, HCPStructure)
            Positions = '1 1    0.00000000000   0.00000000000   0.00000000000\n'\
                        '2 1    0.00000000000   0.57735026919   0.81649658092'

        self.GenerateFile('./posfile', Positions)

    def InsertBasisvectors(self, UiBasisVectors, SpecifiedStructure): 
        """Helper function to write the specified basis to the UI. """
        for i, row in enumerate(UiBasisVectors):
            for j, QLine in enumerate(row):
                QLine.setText(SpecifiedStructure[i,j])
    
    def ClearBasisVectors(self, UiBasisVectors): 
        """Helper function which clears the basis vectors. """
        for row in UiBasisVectors:
            for QtLine in row:
                QtLine.clear()
 
    def GenerateFile(self, filename, string): 
        """Helper function which generates a file given a filename and string for input."""
        file = open(filename, 'w')
        file.write(string)
        file.close()
        return

    def ResetInputs(self, window):
        """
        Reads all QlineEdit-, QCheckBox-, QGroupBox- and 
        QRadioButton-objects from the UI and checks if these are input objects. If 
        they are, the function then clears all inputs and resets the UI to default
        using a reference Reset list or dict depending on object. 

        Inputs:
                window  :   QMainWindow

        Author: Erik Karpelin
        """
        
        import numpy as np
        import os
        from PyQt6.QtWidgets import QLineEdit, QCheckBox, QSpinBox, QGroupBox, QRadioButton

        LineEditList    = window.findChildren(QLineEdit)
        CheckBoxList    = window.findChildren(QCheckBox)
        SpinBoxList     = window.findChildren(QSpinBox)
        GroupBoxList    = window.findChildren(QGroupBox)
        RadioButtonList = window.findChildren(QRadioButton)

        CheckResetList = ['InpPBCCheckC1', 'InpPBCCheckC2', 'InpPBCCheckC3',
                              'InpInitmag3Check', 'InpLLGMeasureBox', 'InpAveBox']
        RadioResetList = ['InpMeasureLLG', 'InpMeasureLLG', 'InpZhangLiCheck',
                           'InpSiteAveCheck', 'InpInitLLG', 'InpDipBruteForceCheck']
        SpinResetDict = {'InpPairSym': 0, 'ImpPairMaptype': 1, 'InpASDLLGAlgh': 1}

        for Box in GroupBoxList:
            if 'Inp' in Box.objectName():
                Box.setChecked(False)
            if Box.objectName() in CheckResetList:
                Box.setChecked(True)
        
        for Line in LineEditList:
            if 'Inp' in Line.objectName():
                Line.clear()
            
        for Check in CheckBoxList:
            if 'Inp' in Check.objectName():
                Check.setChecked(False)
            if Check.objectName() in CheckResetList:
                Check.setChecked(True)

        for Button in RadioButtonList:
            if Button.objectName() in RadioResetList and Button.isChecked() == False:
                Button.toggle()
        
        for SpinBox in SpinBoxList:
            if SpinBox.objectName() in SpinResetDict.keys():
                SpinBox.setValue(SpinResetDict[SpinBox.objectName()])
        
        return
     
    def RemoveInputFile(self, fileobject):
        """
        Helper function to remove specified files
        
        Inputs:
                fileobject  :   list containing name of file
        """
        import os 
        if len(fileobject) > 0:
            if isinstance(fileobject, str):
                os.remove(fileobject.split('/')[-1])
            else:
                os.remove(fileobject[-1])
            fileobject = []
        return
    
    def MagnonQuickSetup(self, window):
        """Setup the GUI for linear spin wave spectra. """
        window.InpQuickRelax.setChecked(True)
        window.InpAMSCheck.setChecked(True)
        return 
    
    def GenDummyMomfile(self, vector, configuration):
        """ Generate a momfile given a vector, a configuration (Ferro/AntiFerro)
            and a posfile in the current directory. 

            Input:
                    vector          :   list
                    configuration   :   int
        """
        import numpy as np
        ASDInputGen.momfile_gotten = True
        ASDInputGen.momfile = './momfile.dummy'
        posfile = np.genfromtxt(self.posfile, ndmin = 2)
        vector = np.array(vector).astype(int)
        numbers = posfile[:,0:2].astype(int)
        magnitude = np.ones((len(numbers), 1))

        if len(numbers) > 1:
            for i in range(1, len(numbers)): 
                moment_vectors = np.vstack((vector, configuration*vector))
        else:
            moment_vectors = [vector]

        momfile_array = np.hstack((numbers, magnitude, moment_vectors))
        momfile_string = ''

        for line in momfile_array:
            for element in line:
                momfile_string += str(int(element)) + "\t"
            momfile_string += "\n"

        self.GenerateFile(ASDInputGen.momfile, momfile_string)

    def update_gui_from_import(self, window, output_files, lattice, input_type, lattice_const):
        import numpy as np
        from PyQt6.QtWidgets import QLineEdit
        """ Update GUI with files and input from file-parser. """
        if len(output_files['jfile']) > 0:
            ASDInputGen.jfile = output_files['jfile']
            window.InpXCCheck.setChecked(True)
        if len(output_files['dmfile']) > 0:
            ASDInputGen.dmfile = output_files['dmfile']
            window.InpDMCheck.setChecked(True)

        ASDInputGen.posfile = output_files['posfile']
        ASDInputGen.posfile_gotten = True
        ASDInputGen.momfile = output_files['momfile']
        ASDInputGen.momfile_gotten = True

        if input_type in ['cif', 'SPRKKR']:
            gui_lines = np.array([line for line in window.findChildren(QLineEdit)
                                if 'InpLineEditC' in line.objectName()]).reshape(3,3)
            lattice = np.array([[str(i) for i in lattice[0]],
                                [str(j) for j in lattice[1]],
                                [str(k) for k in lattice[2]]])
        
        if input_type == 'cif':
            ASDInputGen.UppASDKeywords['geometry']['posfiletype'] = 'D'
            window.InpLineEditAlat.setText(lattice_const)
            self.InsertBasisvectors(gui_lines, lattice)
        
        if input_type == 'SPRKKR': 
            window.InpLineEditAlat.setText(str('{:2.5e}'.format(lattice_const*1e-10)))
            self.InsertBasisvectors(gui_lines, lattice)

        if input_type == 'RSLMTO':
            window.InpLineEditAlat.setText(str('{:2.5e}'.format(lattice_const*1e-10)))
            self.SetStructureTemplate(window, lattice)


    def import_system(self, window):
        """Select import file."""

        from PyQt6 import QtWidgets
        import ASD_GUI.Input_Creator.System_import.ASDImportSys as ASDImportSys
        import ASD_GUI.Input_Creator.System_import.SPRKKR_parser as SPRKKR_parser 

        dlg = QtWidgets.QFileDialog()
        dlg.setFileMode(QtWidgets.QFileDialog.FileMode.AnyFile)
        dlg.setDirectory('.')

        if dlg.exec():
            if window.sender() == window.InpImportCIFButton:
                filename = dlg.selectedFiles()[0]
                output_files, lattice, alat = ASDImportSys.parse_cif(filename)
                input_type = 'cif'
            if window.sender() == window.InpImportSPRKKRButton:
                filename = dlg.selectedFiles()[0]
                output_files, lattice, alat = SPRKKR_parser.parse_sprkkr(filename)
                input_type = 'SPRKKR'
            if window.sender() == window.InpImportRSLMTOButton:
                filename = dlg.selectedFiles()[0]
                output_files, lattice, alat = ASDImportSys.parse_rs_lmto(filename)
                input_type = 'RSLMTO'

            output_files = {'jfile' : output_files[0], 'dmfile': output_files[1],
                            'posfile': output_files[2], 'momfile': output_files[3]}
            self.update_gui_from_import(window, output_files, lattice, input_type, alat)

    ############################################################################
    # @brief Function to ensure that if there is no entry it is set to zero
    # @author Jonathan Chico
    ############################################################################

    def text_to_num(self, text):
        if len(text) > 0:
            num = float(text)
        else:
            num = 0.0
        return num
    ############################################################################
    # @brief Function to constrain the data that can be put into the line-edits.
    # @details The GUI has several editable entries, to ensure that one can only
    # set the correct type of data validators are defined to ensure that only the
    # correct data types can be entered.
    # @author Jonathan Chico
    ############################################################################

    def ASDInputConstrainer(self, window):
        window.InpN1.setValidator(ASDInputGen.IntegerValidator)
        window.InpN2.setValidator(ASDInputGen.IntegerValidator)
        window.InpN3.setValidator(ASDInputGen.IntegerValidator)
        window.InpDipBlockSizeLineEdit.setValidator(
            ASDInputGen.IntegerValidator)
        window.InpMensemble.setValidator(ASDInputGen.IntegerValidator)
        window.InpASDLLGSteps.setValidator(ASDInputGen.IntegerValidator)
        window.InpMCSteps.setValidator(ASDInputGen.IntegerValidator)
        window.InpGNEBMEPSteps.setValidator(ASDInputGen.IntegerValidator)
        window.InpAveStep.setValidator(ASDInputGen.IntegerValidator)
        window.InpAveBuff.setValidator(ASDInputGen.IntegerValidator)
        window.InpCumuStep.setValidator(ASDInputGen.IntegerValidator)
        window.InpCumuBuff.setValidator(ASDInputGen.IntegerValidator)
        window.InitTTrajStep.setValidator(ASDInputGen.IntegerValidator)
        window.InitTTrajBuff.setValidator(ASDInputGen.IntegerValidator)
        window.InpLineEditC1_x.setValidator(ASDInputGen.DoubleValidator)
        window.InpLineEditC1_y.setValidator(ASDInputGen.DoubleValidator)
        window.InpLineEditC1_z.setValidator(ASDInputGen.DoubleValidator)
        window.InpLineEditC2_x.setValidator(ASDInputGen.DoubleValidator)
        window.InpLineEditC2_y.setValidator(ASDInputGen.DoubleValidator)
        window.InpLineEditC2_z.setValidator(ASDInputGen.DoubleValidator)
        window.InpLineEditC3_x.setValidator(ASDInputGen.DoubleValidator)
        window.InpLineEditC3_y.setValidator(ASDInputGen.DoubleValidator)
        window.InpLineEditC3_z.setValidator(ASDInputGen.DoubleValidator)
        window.InpLineEditAlat.setValidator(ASDInputGen.PosDoubleValidator)
        window.InpBextInit_x.setValidator(ASDInputGen.DoubleValidator)
        window.InpBextInit_y.setValidator(ASDInputGen.DoubleValidator)
        window.InpBextInit_z.setValidator(ASDInputGen.DoubleValidator)
        window.InpBextMeasure_x.setValidator(ASDInputGen.DoubleValidator)
        window.InpBextMeasure_y.setValidator(ASDInputGen.DoubleValidator)
        window.InpBextMeasure_z.setValidator(ASDInputGen.DoubleValidator)
        window.InpInitMag2Theta.setValidator(ASDInputGen.DoubleValidator)
        window.InpInitMag2Phi.setValidator(ASDInputGen.DoubleValidator)
        window.InpASDLLGTemp.setValidator(ASDInputGen.PosDoubleValidator)
        window.InpASDLLGDT.setValidator(ASDInputGen.PosDoubleValidator)
        window.InpASDLLGDamp.setValidator(ASDInputGen.PosDoubleValidator)
        window.InpMCTemp.setValidator(ASDInputGen.PosDoubleValidator)
        window.InpMinEig.setValidator(ASDInputGen.PosDoubleValidator)
        window.InpGNEBMEPTol.setValidator(ASDInputGen.PosDoubleValidator)
        window.InpGNEBCITol.setValidator(ASDInputGen.PosDoubleValidator)
        window.InpShe.setValidator(ASDInputGen.DoubleValidator)
        window.InpSheThickness.setValidator(ASDInputGen.PosDoubleValidator)
        window.InpBeta.setValidator(ASDInputGen.PosDoubleValidator)
        window.InpJvec_x.setValidator(ASDInputGen.DoubleValidator)
        window.InpJvec_y.setValidator(ASDInputGen.DoubleValidator)
        window.InpJvec_z.setValidator(ASDInputGen.DoubleValidator)
        window.InpSOTFL.setValidator(ASDInputGen.DoubleValidator)
        window.InpSOTDL.setValidator(ASDInputGen.DoubleValidator)
        window.InpSOTPol_x.setValidator(ASDInputGen.DoubleValidator)
        window.InpSOTPol_y.setValidator(ASDInputGen.DoubleValidator)
        window.InpSOTPol_z.setValidator(ASDInputGen.DoubleValidator)
        window.InpScNStep.setValidator(ASDInputGen.IntegerValidator)
        window.InpScStep.setValidator(ASDInputGen.IntegerValidator)
        return
    ############################################################################
    # @brief Creates a dictionary with the entries needed to generate the inpsd.dat
    # @details The dictionary is created and default values are setup so that a
    # minimal simulation can be run.
    # @author Jonathan Chico
    ############################################################################

    def ASDSetDefaults(self):
        import collections
        # General simulation variables
        ASDInputGen.UppASDKeywords['general'] = collections.OrderedDict()
        ASDInputGen.UppASDKeywords['general']['simid'] = '_UppASD_'
        ASDInputGen.UppASDKeywords['general']['mode'] = 'S'
        ASDInputGen.UppASDKeywords['general']['ip_mode'] = 'N'
        ASDInputGen.UppASDKeywords['general']['Temp'] = 0.001
        ASDInputGen.UppASDKeywords['general']['do_prnstruct'] = 0
        # Geometry variables
        ASDInputGen.UppASDKeywords['geometry'] = collections.OrderedDict()
        ASDInputGen.UppASDKeywords['geometry']['ncell'] = [1, 1, 1]
        ASDInputGen.UppASDKeywords['geometry']['BC'] = ['P', 'P', 'P']
        ASDInputGen.UppASDKeywords['geometry']['cell'] = [
            [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]
        ASDInputGen.UppASDKeywords['geometry']['posfile'] = './posfile'
        ASDInputGen.UppASDKeywords['geometry']['posfiletype'] = 'C'
        ASDInputGen.UppASDKeywords['geometry']['momfile'] = './momfile'
        ASDInputGen.UppASDKeywords['geometry']['do_ralloy'] = 0
        ASDInputGen.UppASDKeywords['geometry']['alat'] = 3e-10
        # Hamiltonian variables
        ASDInputGen.UppASDKeywords['Hamiltonian'] = collections.OrderedDict()
        ASDInputGen.UppASDKeywords['Hamiltonian']['maptype'] = 2
        ASDInputGen.UppASDKeywords['Hamiltonian']['Sym'] = 0
        ASDInputGen.UppASDKeywords['Hamiltonian']['exchange'] = './jfile'
        ASDInputGen.UppASDKeywords['Hamiltonian']['do_dm'] = 0
        ASDInputGen.UppASDKeywords['Hamiltonian']['dm'] = './dmfile'
        ASDInputGen.UppASDKeywords['Hamiltonian']['do_anisotropy'] = 0
        ASDInputGen.UppASDKeywords['Hamiltonian']['anisotropy'] = './kfile'
        ASDInputGen.UppASDKeywords['Hamiltonian']['do_pd'] = 0
        ASDInputGen.UppASDKeywords['Hamiltonian']['pd'] = './pdfile'
        ASDInputGen.UppASDKeywords['Hamiltonian']['do_bq'] = 0
        ASDInputGen.UppASDKeywords['Hamiltonian']['bq'] = './bqfile'
        ASDInputGen.UppASDKeywords['Hamiltonian']['do_biqdm'] = 0
        ASDInputGen.UppASDKeywords['Hamiltonian']['biqdm'] = './biqdmfile'
        ASDInputGen.UppASDKeywords['Hamiltonian']['do_dip'] = 0
        ASDInputGen.UppASDKeywords['Hamiltonian']['block_size'] = 1
        ASDInputGen.UppASDKeywords['Hamiltonian']['hfield'] = [0.0, 0.0, 0.0]
        ASDInputGen.UppASDKeywords['Hamiltonian']['ip_hfield'] = [
            0.0, 0.0, 0.0]
        # LLG measure  variables
        ASDInputGen.UppASDKeywords['LLG_mphase'] = collections.OrderedDict()
        ASDInputGen.UppASDKeywords['LLG_mphase']['SDEAlgh'] = 1
        ASDInputGen.UppASDKeywords['LLG_mphase']['damping'] = 0.01
        ASDInputGen.UppASDKeywords['LLG_mphase']['timestep'] = 1e-16
        ASDInputGen.UppASDKeywords['LLG_mphase']['Nstep'] = 1000
        # LLG measure  variables
        ASDInputGen.UppASDKeywords['GNEB_mphase'] = collections.OrderedDict()
        ASDInputGen.UppASDKeywords['GNEB_mphase']['do_gneb'] = 'N'
        ASDInputGen.UppASDKeywords['GNEB_mphase']['do_gneb_ci'] = 'N'
        ASDInputGen.UppASDKeywords['GNEB_mphase']['mep_ftol'] = 0.001
        ASDInputGen.UppASDKeywords['GNEB_mphase']['mep_ftol_ci'] = 0.00000001
        ASDInputGen.UppASDKeywords['GNEB_mphase']['mep_itrmax'] = 10000000
        ASDInputGen.UppASDKeywords['GNEB_mphase']['eig_zero'] = 0.0001
        # LLG initial phase variables
        ASDInputGen.UppASDKeywords['LLG_iphase'] = collections.OrderedDict()
        ASDInputGen.UppASDKeywords['LLG_iphase']['ip_nphase'] = 0
        # VPO initial phase variables
        ASDInputGen.UppASDKeywords['VPO_iphase'] = collections.OrderedDict()
        # MC measure variables
        ASDInputGen.UppASDKeywords['MC_mphase'] = collections.OrderedDict()
        ASDInputGen.UppASDKeywords['MC_mphase']['mcnstep'] = 10000
        # MC initial phase variables
        ASDInputGen.UppASDKeywords['MC_iphase'] = collections.OrderedDict()
        ASDInputGen.UppASDKeywords['MC_iphase']['ip_mcanneal'] = 0
        # Magnetization variables
        ASDInputGen.UppASDKeywords['Mag'] = collections.OrderedDict()
        ASDInputGen.UppASDKeywords['Mag']['Mensemble'] = 1
        ASDInputGen.UppASDKeywords['Mag']['initmag'] = 3
        ASDInputGen.UppASDKeywords['Mag']['restartfile'] = './restart.dummy.dat'
        # Correlation variables
        ASDInputGen.UppASDKeywords['Mag_corr'] = collections.OrderedDict()
        ASDInputGen.UppASDKeywords['Mag_corr']['do_sc'] = 'N'
        ASDInputGen.UppASDKeywords['Mag_corr']['sc_step'] = 10
        ASDInputGen.UppASDKeywords['Mag_corr']['sc_nstep'] = 100
        ASDInputGen.UppASDKeywords['Mag_corr']['qpoints'] = 'F'
        ASDInputGen.UppASDKeywords['Mag_corr']['qfile'] = './qfile'
        ASDInputGen.UppASDKeywords['Mag_corr']['do_ams'] = 'N'
        # Spin-Torques variables
        ASDInputGen.UppASDKeywords['spintorque'] = collections.OrderedDict()
        ASDInputGen.UppASDKeywords['spintorque']['stt'] = 'N'
        ASDInputGen.UppASDKeywords['spintorque']['adibeta'] = 0.01
        ASDInputGen.UppASDKeywords['spintorque']['jvec'] = [0.0, 0.0, 0.0]
        ASDInputGen.UppASDKeywords['spintorque']['do_sot'] = 'N'
        ASDInputGen.UppASDKeywords['spintorque']['do_she'] = 'N'
        ASDInputGen.UppASDKeywords['spintorque']['sot_field'] = 0.0
        ASDInputGen.UppASDKeywords['spintorque']['sot_damping'] = 0.0
        ASDInputGen.UppASDKeywords['spintorque']['sot_pol_vec'] = [
            0.0, 0.0, 0.0]
        ASDInputGen.UppASDKeywords['spintorque']['thick_ferro'] = 1.0
        ASDInputGen.UppASDKeywords['spintorque']['she_angle'] = 0.0
        # Prn avrg variables
        ASDInputGen.UppASDKeywords['averages'] = collections.OrderedDict()
        ASDInputGen.UppASDKeywords['averages']['do_avrg'] = 'Y'
        ASDInputGen.UppASDKeywords['averages']['avrg_step'] = 1000
        ASDInputGen.UppASDKeywords['averages']['avrg_buff'] = 100
        ASDInputGen.UppASDKeywords['averages']['do_proj_avrg'] = 'N'
        ASDInputGen.UppASDKeywords['averages']['do_projch_avrg'] = 'N'
        # Prn tottraj variables
        ASDInputGen.UppASDKeywords['trajectories'] = collections.OrderedDict()
        ASDInputGen.UppASDKeywords['trajectories']['do_tottraj'] = 'N'
        ASDInputGen.UppASDKeywords['trajectories']['tottraj_step'] = 1000
        ASDInputGen.UppASDKeywords['trajectories']['tottraj_buff'] = 100
        # Prn cumulants variables
        ASDInputGen.UppASDKeywords['cumulants'] = collections.OrderedDict()
        ASDInputGen.UppASDKeywords['cumulants']['do_cumu'] = 'Y'
        ASDInputGen.UppASDKeywords['cumulants']['cumu_step'] = 1000
        ASDInputGen.UppASDKeywords['cumulants']['cumu_buff'] = 100
        # Prn skyrmion
        ASDInputGen.UppASDKeywords['topology'] = collections.OrderedDict()
        ASDInputGen.UppASDKeywords['topology']['skyrno'] = 'N'
        ASDInputGen.UppASDKeywords['topology']['skyno_step'] = 1000
        ASDInputGen.UppASDKeywords['topology']['skyno_buff'] = 100
        # Prn energy
        ASDInputGen.UppASDKeywords['energy'] = collections.OrderedDict()
        ASDInputGen.UppASDKeywords['energy']['plotenergy'] = 1
        # Hessian
        ASDInputGen.UppASDKeywords['Hessians'] = collections.OrderedDict()
        ASDInputGen.UppASDKeywords['Hessians']['do_hess_ini'] = 'N'
        ASDInputGen.UppASDKeywords['Hessians']['do_hess_fin'] = 'N'
        ASDInputGen.UppASDKeywords['Hessians']['do_hess_sp'] = 'N'
        return
    ############################################################################
    # @brief Function to clean the input generator dictionary to remove empty entries
    # @details The function makes sute to eliminate un-needed entries to get the
    # minimal inpsd.dat
    # @author Jonathan Chico
    ############################################################################

    def clean_var(self):
        tol = 1e-10
        if ASDInputGen.UppASDKeywords['spintorque']['stt'] == 'N' and ASDInputGen.UppASDKeywords['spintorque']['do_she'] == 'N':
            del ASDInputGen.UppASDKeywords['spintorque']['jvec']
        if ASDInputGen.UppASDKeywords['spintorque']['stt'] == 'N':
            del ASDInputGen.UppASDKeywords['spintorque']['adibeta']
            del ASDInputGen.UppASDKeywords['spintorque']['stt']
        if ASDInputGen.UppASDKeywords['spintorque']['do_she'] == 'N':
            del ASDInputGen.UppASDKeywords['spintorque']['she_angle']
            del ASDInputGen.UppASDKeywords['spintorque']['thick_ferro']
            del ASDInputGen.UppASDKeywords['spintorque']['do_she']
        if ASDInputGen.UppASDKeywords['spintorque']['do_sot'] == 'N':
            del ASDInputGen.UppASDKeywords['spintorque']['sot_pol_vec']
            del ASDInputGen.UppASDKeywords['spintorque']['sot_damping']
            del ASDInputGen.UppASDKeywords['spintorque']['sot_field']
            del ASDInputGen.UppASDKeywords['spintorque']['do_sot']
        # S(q,w) and AMS flags
        if ASDInputGen.UppASDKeywords['Mag_corr']['qpoints'] != 'F':
            del ASDInputGen.UppASDKeywords['Mag_corr']['qfile']
        if ASDInputGen.UppASDKeywords['Mag_corr']['do_sc'] == 'N' and ASDInputGen.UppASDKeywords['Mag_corr']['do_ams'] == 'N':
            del ASDInputGen.UppASDKeywords['Mag_corr']['qfile']
            del ASDInputGen.UppASDKeywords['Mag_corr']['qpoints']
        if ASDInputGen.UppASDKeywords['Mag_corr']['do_sc'] == 'N':
            del ASDInputGen.UppASDKeywords['Mag_corr']['sc_step']
            del ASDInputGen.UppASDKeywords['Mag_corr']['sc_nstep']
            del ASDInputGen.UppASDKeywords['Mag_corr']['do_sc']
        if ASDInputGen.UppASDKeywords['Mag_corr']['do_ams'] == 'N':
            del ASDInputGen.UppASDKeywords['Mag_corr']['do_ams']
        # dipolar flags
        if ASDInputGen.UppASDKeywords['Hamiltonian']['do_dip'] != 2:
            del ASDInputGen.UppASDKeywords['Hamiltonian']['block_size']
        if ASDInputGen.UppASDKeywords['Hamiltonian']['do_dip'] == 0:
            del ASDInputGen.UppASDKeywords['Hamiltonian']['do_dip']
        # DMI flags
        if ASDInputGen.UppASDKeywords['Hamiltonian']['do_dm'] == 0:
            del ASDInputGen.UppASDKeywords['Hamiltonian']['do_dm']
            del ASDInputGen.UppASDKeywords['Hamiltonian']['dm']
        # Anisotropy flags
        if ASDInputGen.UppASDKeywords['Hamiltonian']['do_anisotropy'] == 0:
            del ASDInputGen.UppASDKeywords['Hamiltonian']['do_anisotropy']
            del ASDInputGen.UppASDKeywords['Hamiltonian']['anisotropy']
        # Biquadratic interaction flags
        if ASDInputGen.UppASDKeywords['Hamiltonian']['do_bq'] == 0:
            del ASDInputGen.UppASDKeywords['Hamiltonian']['do_bq']
            del ASDInputGen.UppASDKeywords['Hamiltonian']['bq']
        # Pseudo dipolar flags
        if ASDInputGen.UppASDKeywords['Hamiltonian']['do_pd'] == 0:
            del ASDInputGen.UppASDKeywords['Hamiltonian']['do_pd']
            del ASDInputGen.UppASDKeywords['Hamiltonian']['pd']
        # Biquadratic DM interaction flags
        if ASDInputGen.UppASDKeywords['Hamiltonian']['do_biqdm'] == 0:
            del ASDInputGen.UppASDKeywords['Hamiltonian']['do_biqdm']
            del ASDInputGen.UppASDKeywords['Hamiltonian']['biqdm']
        # IP Mode flags
        if ASDInputGen.UppASDKeywords['general']['ip_mode'] == 'N':
            del ASDInputGen.UppASDKeywords['general']['ip_mode']
            del ASDInputGen.UppASDKeywords['Hamiltonian']['ip_hfield']
            del ASDInputGen.UppASDKeywords['MC_iphase']
            del ASDInputGen.UppASDKeywords['LLG_iphase']
            del ASDInputGen.UppASDKeywords['VPO_iphase']
        else:
            if ASDInputGen.UppASDKeywords['general']['ip_mode'] != 'M' and ASDInputGen.UppASDKeywords['general']['ip_mode'] != 'H' and ASDInputGen.UppASDKeywords['general']['ip_mode'] != 'SX':
                del ASDInputGen.UppASDKeywords['MC_iphase']
            if ASDInputGen.UppASDKeywords['general']['ip_mode'] != 'S':
                del ASDInputGen.UppASDKeywords['LLG_iphase']
            if ASDInputGen.UppASDKeywords['general']['ip_mode'] != 'G':
                del ASDInputGen.UppASDKeywords['VPO_iphase']
        # Mode flags
        if ASDInputGen.UppASDKeywords['general']['mode'] != 'M' and ASDInputGen.UppASDKeywords['general']['mode'] != 'H':
            del ASDInputGen.UppASDKeywords['MC_mphase']
        if ASDInputGen.UppASDKeywords['general']['mode'] != 'S':
            del ASDInputGen.UppASDKeywords['LLG_mphase']
        if ASDInputGen.UppASDKeywords['general']['mode'] != 'G':
            del ASDInputGen.UppASDKeywords['GNEB_mphase']
        # Initmag flags
        if ASDInputGen.UppASDKeywords['Mag']['initmag'] != 4 and ASDInputGen.UppASDKeywords['Mag']['initmag'] != 7:
            del ASDInputGen.UppASDKeywords['Mag']['restartfile']
        # Measurement field flag
        if ASDInputGen.UppASDKeywords['Hamiltonian']['hfield'][0]**2 +\
                ASDInputGen.UppASDKeywords['Hamiltonian']['hfield'][1]**2 +\
                ASDInputGen.UppASDKeywords['Hamiltonian']['hfield'][2]**2 < tol:
            del ASDInputGen.UppASDKeywords['Hamiltonian']['hfield']
        # Random alloy flags
        if ASDInputGen.UppASDKeywords['geometry']['do_ralloy'] == 0:
            del ASDInputGen.UppASDKeywords['geometry']['do_ralloy']
        # Prn trajectories
        if ASDInputGen.UppASDKeywords['trajectories']['do_tottraj'] == 'N':
            del ASDInputGen.UppASDKeywords['trajectories']
        # Prn averages
        if ASDInputGen.UppASDKeywords['averages']['do_avrg'] == 'N':
            del ASDInputGen.UppASDKeywords['averages']
        # Prn topology
        if ASDInputGen.UppASDKeywords['topology']['skyrno'] == 'N':
            del ASDInputGen.UppASDKeywords['topology']
        # Prn cumulants
        if ASDInputGen.UppASDKeywords['cumulants']['do_cumu'] == 'N':
            del ASDInputGen.UppASDKeywords['cumulants']
        # Prn averages
        if ASDInputGen.UppASDKeywords['averages']['do_avrg'] == 'N':
            del ASDInputGen.UppASDKeywords['averages']
        else:
            if ASDInputGen.UppASDKeywords['averages']['do_proj_avrg'] == 'N':
                del ASDInputGen.UppASDKeywords['averages']['do_proj_avrg']
            if ASDInputGen.UppASDKeywords['averages']['do_projch_avrg'] == 'N':
                del ASDInputGen.UppASDKeywords['averages']['do_projch_avrg']
        # Hessians
        if ASDInputGen.UppASDKeywords['Hessians']['do_hess_ini'] == 'N':
            del ASDInputGen.UppASDKeywords['Hessians']['do_hess_ini']
        if ASDInputGen.UppASDKeywords['Hessians']['do_hess_fin'] == 'N':
            del ASDInputGen.UppASDKeywords['Hessians']['do_hess_fin']
        if ASDInputGen.UppASDKeywords['Hessians']['do_hess_sp'] == 'N':
            del ASDInputGen.UppASDKeywords['Hessians']['do_hess_sp']
        if ASDInputGen.UppASDKeywords['Mag']['relaxed_if'] == 'N':
            del ASDInputGen.UppASDKeywords['Mag']['relaxed_if']
        # for name in ASDInputGen.UppASDKeywords:
        #    if len(ASDInputGen.UppASDKeywords[name])==0:
        #        del ASDInputGen.UppASDKeywords[name]
        return
    ############################################################################
    # @brief Function to write a standard inpsd.dat and a inpsd.yaml file
    # @author Jonathan Chico
    ############################################################################
    def write_inpsd(self):
        """Function to write a standard inpsd.dat and a inpsd.yaml file

        Author
        ----------
        Jonathan Chico
        """
        import collections
        import yaml

        yaml.add_representer(collections.OrderedDict, lambda dumper, data: dumper.represent_mapping(
            'tag:yaml.org,2002:map', data.items()))
        with open('inpsd.yaml', 'w') as outfile:
            yaml.dump(ASDInputGen.UppASDKeywords,
                      outfile, default_flow_style=False)

        inpsd_file = open('inpsd.dat', 'w')
        for name in ASDInputGen.UppASDKeywords:
            for descriptor in ASDInputGen.UppASDKeywords[name]:
                current = ASDInputGen.UppASDKeywords[name][descriptor]
                if isinstance(current, list):
                    if len(descriptor) > 0:
                        inpsd_file.write('{descriptor}  '.format(**locals()))
                    for ii in range(len(current)):
                        line = current[ii]
                        if isinstance(line, list):
                            for jj in range(len(line)):
                                entry = line[jj]
                                inpsd_file.write(
                                    '{entry}  '.format(**locals()))
                            inpsd_file.write('\n')
                        else:
                            inpsd_file.write('{line}  '.format(**locals()))
                    inpsd_file.write('\n')
                elif isinstance(current, tuple):
                    inpsd_file.write('{descriptor} '.format(**locals()))
                    for ii in range(len(current)):
                        entry = current[ii]
                        inpsd_file.write('{entry}  '.format(**locals()))
                    inpsd_file.write('\n')
                else:
                    inpsd_file.write(
                        '{descriptor}  {current}\n'.format(**locals()))
            inpsd_file.write('\n')
        return
