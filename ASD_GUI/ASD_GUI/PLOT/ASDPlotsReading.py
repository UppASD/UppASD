"""  @package ASDPlotsReading
Contains the definitions needed to read the data for the matplotlib plots.
It allows one also to find the input files making use of dialog windows.

Author
----------
Jonathan Chico
"""
################################################################################
# @brief Class to read the data needed to plot via the GUI.
# @details It has a function to obtain the file names via a dialog window. It can also
# obtain the needed file names and data for the plotting making use of the uppasd.yaml
# file.
# @author Jonathan Chico
################################################################################
class ReadPlotData():
    ############################################################################
    # @brief Class constructor
    # @details Class constructor for the reading of the plotting data. Contains
    # the definitions of several containers for the file names and variables which
    # are used during plotting.
    # @author Jonathan Chico
    ############################################################################
    def __init__(self):
        ReadPlotData.sc_step=1
        ReadPlotData.timestep=1
        ReadPlotData.h_mev=4.135667662e-12
        ReadPlotData.amsfile=[]
        ReadPlotData.sqwfile=[]
        ReadPlotData.yamlfile=[]
        ReadPlotData.averages=[]
        ReadPlotData.totenergy=[]
        ReadPlotData.trajectory=[]
        ReadPlotData.not_read_sqw=True
        ReadPlotData.not_read_ams=True
        ReadPlotData.not_read_yaml=True
        ReadPlotData.not_read_averages=True
        ReadPlotData.not_read_totenergy=True
        ReadPlotData.not_read_trajectory=True

        ReadPlotData.ams_file_present=False
        ReadPlotData.sqw_file_present=False
        ReadPlotData.ave_file_present=False
        ReadPlotData.ene_file_present=False
        ReadPlotData.yaml_file_present=False
        ReadPlotData.trajectory_file_present=False
        return

    ############################################################################
    # @brief Function to get the file names for the different types of plots
    # @author Jonathan Chico
    ############################################################################
    def getFileName(self,window):
        from PyQt5 import QtWidgets
        dlg = QtWidgets.QFileDialog()
        dlg.setFileMode(QtWidgets.QFileDialog.AnyFile)
        dlg.setDirectory('.')
        if dlg.exec_():
            if window.sender()==window.actionYaml_File:
                ReadPlotData.yamlfile=dlg.selectedFiles()[0]
                ReadPlotData.not_read_yaml=True
            if window.sender()==window.actionS_q_w_File:
                ReadPlotData.sqwfile=dlg.selectedFiles()[0]
                ReadPlotData.not_read_sqw=True
            if window.sender()==window.actionAMS_File:
                ReadPlotData.amsfile=dlg.selectedFiles()[0]
                ReadPlotData.not_read_ams=True
            if window.sender()==window.actionAverages_File:
                ReadPlotData.averages=dlg.selectedFiles()[0]
                ReadPlotData.not_read_averages=True
            if window.sender()==window.actionTrajectory_File:
                ReadPlotData.trajectory=dlg.selectedFiles()
                ReadPlotData.not_read_trajectory=True
            if window.sender()==window.actionTot_Energy_File:
                ReadPlotData.totenergy=dlg.selectedFiles()
                ReadPlotData.not_read_totenergy=True
        return
    ############################################################################
    # @brief Wrapper to read the plotting data
    # @details This wrapper allows one to read the following files:
    #
    #   - Yaml file to read the measurement flags resulting from an \c UppASD simulation.
    #   - Averages file.
    #   - Single atom trajectory file.
    #   - Dynamical structure factor file \f$\mathbf{S}(\mathbf{q},\omega)\f$
    #   - Adiabatic magnon spectra file.
    # @author Jonathan Chico
    ############################################################################
    def PlotReadingWrapper(self,file_names,window):
        import UI.ASDInputWindows
        import glob
        ReadPlotData.yamlfile   = file_names[0]
        ReadPlotData.amsfile    = file_names[1]
        ReadPlotData.sqwfile    = file_names[2]
        ReadPlotData.averages   = file_names[3]
        ReadPlotData.trajectory = file_names[4]
        ReadPlotData.totenergy  = file_names[5]

        if len(ReadPlotData.yamlfile)>0:
            ReadPlotData.UppASDYamlInfo=ReadPlotData.yamlfile
            #-------------------------------------------------------------------
            # If the yaml file is found one must read it and find its data first
            #-------------------------------------------------------------------
            if ReadPlotData.not_read_yaml:
                self.Yaml_Read_Wrapper(ReadPlotData.UppASDYamlInfo,window)
                ReadPlotData.not_read_yaml=False
                ReadPlotData.yaml_file_present=True
        else:
            print("No file name selected from menu. Trying to find a 'asd.yaml' file")
            ReadPlotData.yamlfile=glob.glob("asd.yaml")
            if len(ReadPlotData.yamlfile)>0:
                ReadPlotData.yamlfile=ReadPlotData.yamlfile[0]
                ReadPlotData.UppASDYamlInfo=ReadPlotData.yamlfile
                # If the yaml file is found one must read it and find its data first
                if ReadPlotData.not_read_yaml:
                    self.Yaml_Read_Wrapper(ReadPlotData.UppASDYamlInfo,window)
                    ReadPlotData.not_read_yaml=False
                    ReadPlotData.yaml_file_present=True
            else:
                print("No yaml file found, searching for other files")
                if len(window.InpPlotDt.text())>0:
                    ReadPlotData.timestep=float(window.InpPlotDt.text())*1e9
                else:
                    ReadPlotData.timestep=1

                if len(window.InpSqwSCStep.text())>0:
                    ReadPlotData.sc_step=int(window.InpSqwSCStep.text())
                else:
                    ReadPlotData.sc_step=1

                if len(window.InpSqwSCNStep.text())>0:
                    ReadPlotData.sc_nstep=int(window.InpSqwSCNStep.text())
                else:
                    ReadPlotData.sc_nstep=1

                ReadPlotData.hf_scale=ReadPlotData.h_mev/(ReadPlotData.timestep*1e-9*ReadPlotData.sc_step*ReadPlotData.sc_nstep)/2
                #---------------------------------------------------------------
                # Trying to find an ams file
                #---------------------------------------------------------------
                if ((window.sender()==window.AMSDispCheckBox or window.sender()==window.SqwDispCheckBox) \
                and window.AMSDispCheckBox.isChecked()) \
                or (window.sender()==window.actionS_q_w and window.AMSDispCheckBox.isChecked()):
                    if len(ReadPlotData.amsfile)>0:
                        if ReadPlotData.not_read_ams:
                            (ReadPlotData.ams_data_x,ReadPlotData.ams_data_y,ReadPlotData.ams_ax_label,\
                            ReadPlotData.ams_label)=\
                            self.read_ams(ReadPlotData.amsfile)
                            ReadPlotData.not_read_ams=False
                            ReadPlotData.ams_file_present=True
                    else:
                        print("No file name selected from menu. Trying to find a 'ams.*.out' file")
                        ReadPlotData.amsfile=glob.glob("ams.*.out")
                        if len(ReadPlotData.amsfile)>0:
                            ReadPlotData.amsfile=ReadPlotData.amsfile[0]
                            if ReadPlotData.not_read_ams:
                                (ReadPlotData.ams_data_x,ReadPlotData.ams_data_y,\
                                ReadPlotData.ams_ax_label,ReadPlotData.ams_label)=\
                                self.read_ams(ReadPlotData.amsfile)
                                ReadPlotData.not_read_ams=False
                                ReadPlotData.ams_file_present=True
                #---------------------------------------------------------------
                # Trying to find the sqw file
                #---------------------------------------------------------------
                if ((window.sender()==window.SqwDispCheckBox or window.sender()==window.AMSDispCheckBox)
                and window.SqwDispCheckBox.isChecked()) \
                or (window.sender()==window.actionS_q_w and window.SqwDispCheckBox.isChecked()):
                    if ReadPlotData.sc_step==1:
                        window.sc_step_Error_Window=UI.ASDInputWindows.Error_Window()
                        window.sc_step_Error_Window.FunMsg.setText("Expect the Unexpected")
                        window.sc_step_Error_Window.ErrorMsg.setText("Error: No sampling frequency set assuming 1 the scale of the plot might be wrong.")
                        window.sc_step_Error_Window.show()
                    if ReadPlotData.sc_nstep==1:
                        window.sc_nstep_Error_Window=UI.ASDInputWindows.Error_Window()
                        window.sc_nstep_Error_Window.FunMsg.setText("Expect the Unexpected")
                        window.sc_nstep_Error_Window.ErrorMsg.setText("Error: No number frequencies set assuming 1 the scale of the plot might be wrong.")
                        window.sc_step_Error_Window.show()
                    if ReadPlotData.timestep==1:
                        window.time_step_Error_Window=UI.ASDInputWindows.Error_Window()
                        window.time_step_Error_Window.FunMsg.setText("Expect the Unexpected")
                        window.time_step_Error_Window.ErrorMsg.setText("Error: No time step set assuming 1 the scale of the plot might be wrong.")
                        window.time_step_Error_Window.show()
                    if len(ReadPlotData.sqwfile)>0:
                        if ReadPlotData.not_read_sqw:
                            (ReadPlotData.sqw_data,ReadPlotData.sqw_labels,ReadPlotData.ax_limits)=\
                            self.read_sqw(ReadPlotData.sqwfile)
                            ReadPlotData.not_read_sqw=False
                            ReadPlotData.sqw_file_present=True
                    else:
                        print("No file name selected from menu. Trying to find a 'sqw.*.out' file")
                        ReadPlotData.sqwfile=glob.glob("sqw.*.out")
                        if len(ReadPlotData.sqwfile)>0:
                            ReadPlotData.sqwfile=ReadPlotData.sqwfile[0]
                            if ReadPlotData.not_read_sqw:
                                (ReadPlotData.sqw_data,ReadPlotData.sqw_labels,ReadPlotData.ax_limits)=\
                                self.read_sqw(ReadPlotData.sqwfile)
                                ReadPlotData.not_read_sqw=False
                                ReadPlotData.sqw_file_present=True
                #---------------------------------------------------------------
                # Trying to find the averages file
                #---------------------------------------------------------------
                if window.sender()==window.actionAverages:
                    if len(ReadPlotData.averages)>0:
                        if ReadPlotData.not_read_averages:
                            (ReadPlotData.mag_data,ReadPlotData.itr_data,           \
                            ReadPlotData.mag_labels)=                               \
                            self.read_gen_plot_data(ReadPlotData.averages)
                            if ReadPlotData.timestep!=1:
                                ReadPlotData.mag_axes=['Time [ns]',r'Magnetization [$\mu_B$]']
                            else:
                                ReadPlotData.mag_axes=['# Iterations',r'Magnetization [$\mu_B$]']
                            ReadPlotData.not_read_averages=False
                            ReadPlotData.ave_file_present=True
                    else:
                        print("No file name selected from menu. Trying to find a 'averages.*.out' file")
                        ReadPlotData.averages=glob.glob("averages.*.out")
                        if len(ReadPlotData.averages)>0:
                            ReadPlotData.averages=ReadPlotData.averages[0]
                            if ReadPlotData.not_read_averages:
                                (ReadPlotData.mag_data,ReadPlotData.itr_data,       \
                                ReadPlotData.mag_labels)=                           \
                                self.read_gen_plot_data(ReadPlotData.averages)
                                if ReadPlotData.timestep!=1:
                                    ReadPlotData.mag_axes=['Time [ns]',r'Magnetization [$\mu_B$]']
                                else:
                                    ReadPlotData.mag_axes=['# Iterations',r'Magnetization [$\mu_B$]']
                                ReadPlotData.not_read_averages=False
                                ReadPlotData.ave_file_present=True
                #---------------------------------------------------------------
                # Trying to find the totalenergy file
                #---------------------------------------------------------------
                if window.sender()==window.actionTotEnergy:
                    if len(ReadPlotData.totenergy)>0:
                        if ReadPlotData.not_read_totenergy:
                            (ReadPlotData.ene_data,ReadPlotData.itr_data,           \
                            ReadPlotData.ene_labels)=                               \
                            self.read_gen_plot_data(ReadPlotData.totenergy)
                            if ReadPlotData.timestep!=1:
                                ReadPlotData.ene_axes=['Time [ns]',r'Energy/spin [mRy]']
                            else:
                                ReadPlotData.ene_axes=['# Iterations',r'Energy/spin [mRy]']
                            ReadPlotData.not_read_totenergy=False
                            ReadPlotData.ene_file_present=True
                    else:
                        print("No file name selected from menu. Trying to find a 'totenergy.*.out' file")
                        ReadPlotData.totenergy=glob.glob("totenergy.*.out")
                        if len(ReadPlotData.totenergy)>0:
                            ReadPlotData.totenergy=ReadPlotData.totenergy[0]
                            if ReadPlotData.not_read_totenergy:
                                (ReadPlotData.ene_data,ReadPlotData.itr_data,       \
                                ReadPlotData.ene_labels)=                           \
                                self.read_gen_plot_data(ReadPlotData.totenergy)
                                if ReadPlotData.timestep!=1:
                                    ReadPlotData.ene_axes=['Time [ns]',r'Energy/spin [mRy]']
                                else:
                                    ReadPlotData.ene_axes=['# Iterations',r'Energy/spin [mRy]']
                                ReadPlotData.not_read_totenergy=False
                                ReadPlotData.ene_file_present=True
                #---------------------------------------------------------------
                # Trying to find the trajectory files
                #---------------------------------------------------------------
                if window.sender()==window.actionTrajectory:
                    if len(ReadPlotData.trajectory)>0:
                        if ReadPlotData.not_read_trajectory:
                            (ReadPlotData.traj_label,ReadPlotData.traj_data_x,\
                            ReadPlotData.traj_data_y,ReadPlotData.traj_data_z)=\
                            self.read_trajectories(ReadPlotData.trajectory)
                            ReadPlotData.not_read_trajectory=False
                            ReadPlotData.trajectory_file_present=True
                    else:
                        print("No file name selected from menu. Trying to find a 'trajectory.*.out' file")
                        ReadPlotData.trajectory=glob.glob("trajectory.*.out")
                        if len(ReadPlotData.trajectory)>0:
                            (ReadPlotData.traj_label,ReadPlotData.traj_data_x,\
                            ReadPlotData.traj_data_y,ReadPlotData.traj_data_z)=\
                            self.read_trajectories(ReadPlotData.trajectory)
                            ReadPlotData.not_read_trajectory=False
                            ReadPlotData.trajectory_file_present=True
        return
    ############################################################################
    # @brief Wrapper to read all the file names and data as defined by the yaml file
    # @details This finds the file names and parameters needed for the plotting
    # @author Jonathan Chico
    ############################################################################
    def Yaml_Read_Wrapper(self,filename,window):
        import yaml
        import glob
        with open(filename, 'r') as stream:
            try:
                ReadPlotData.sim=yaml.load(stream)
            except yaml.YAMLError as exc:
                print(exc)
        if "timestep" in ReadPlotData.sim["siminfo"]:
            ReadPlotData.timestep=float(ReadPlotData.sim["siminfo"]["timestep"])*1e9
            window.InpPlotDt.setText(str(ReadPlotData.timestep*1e-9))
        if "sc_step" in ReadPlotData.sim["siminfo"]:
            ReadPlotData.sc_step=int(ReadPlotData.sim["siminfo"]["sc_step"])
            window.InpSqwSCStep.setText(str(ReadPlotData.sc_step))
        if "sc_step" in ReadPlotData.sim["siminfo"]:
            ReadPlotData.sc_nstep=int(ReadPlotData.sim["siminfo"]["sc_nstep"])
            window.InpSqwSCNStep.setText(str(ReadPlotData.sc_nstep))
        #############################################################################
        # Read the averages
        #############################################################################
        if ReadPlotData.sim["measurables"]["averages"]:
            ReadPlotData.averages="averages."+ReadPlotData.sim["simid"]+".out"
            (ReadPlotData.mag_data,ReadPlotData.itr_data,ReadPlotData.mag_labels)=  \
            self.read_gen_plot_data(ReadPlotData.averages)
            if ReadPlotData.timestep!=1:
                ReadPlotData.mag_axes=['Time [ns]',r'Magnetization [$\mu_B$]']
            else:
                ReadPlotData.mag_axes=['# Iterations',r'Magnetization [$\mu_B$]']
            ReadPlotData.not_read_averages=False
            ReadPlotData.ave_file_present=True
        #############################################################################
        # Read the energy
        #############################################################################
        if ReadPlotData.sim["measurables"]["totenergy"]:
            (ReadPlotData.ene_data,ReadPlotData.itr_data,ReadPlotData.ene_labels)=  \
            self.read_gen_plot_data(ReadPlotData.totenergy)
            if ReadPlotData.timestep!=1:
                ReadPlotData.ene_axes=['Time [ns]',r'Energy/spin [mRy]']
            else:
                ReadPlotData.ene_axes=['# Iterations',r'Energy/spin [mRy]']
            ReadPlotData.not_read_totenergy=False
            ReadPlotData.ene_file_present=True
        #############################################################################
        # Read the trajectories
        #############################################################################
        if ReadPlotData.sim["measurables"]["trajectories"]:
            ReadPlotData.trajectory=glob.glob("trajectory."+ReadPlotData.sim["simid"]+"*.out")
            (ReadPlotData.traj_label,ReadPlotData.traj_data_x,ReadPlotData.traj_data_y,\
            ReadPlotData.traj_data_z)=\
            self.read_trajectories(ReadPlotData.trajectory)
            ReadPlotData.not_read_trajectory=False
            ReadPlotData.trajectory_file_present=True
        #############################################################################
        # Read the momentfile
        #############################################################################
        if ReadPlotData.sim["measurables"]["moments"]:
            ReadPlotData.moments="moments."+ReadPlotData.sim["simid"]+".out"
        #############################################################################
        # Read the Sqw
        #############################################################################
        if ReadPlotData.sim["measurables"]["sqw"]:
            ReadPlotData.hf_scale=ReadPlotData.h_mev/(ReadPlotData.timestep*1e-9*ReadPlotData.sc_step)
            ReadPlotData.sqwfile="sqw."+ReadPlotData.sim["simid"]+".out"
            (ReadPlotData.sqw_data,ReadPlotData.sqw_labels,ReadPlotData.ax_limits)=\
            self.read_sqw(ReadPlotData.sqwfile)
            ReadPlotData.not_read_sqw=False
            ReadPlotData.sqw_file_present=True
        #############################################################################
        # Read the AMS
        #############################################################################
        if ReadPlotData.sim["measurables"]["ams"]:
            ReadPlotData.amsfile="ams."+ReadPlotData.sim["simid"]+".out"
            (ReadPlotData.ams_data_x,ReadPlotData.ams_data_y,ReadPlotData.ams_ax_label,\
            ReadPlotData.ams_label)=\
            self.read_ams(ReadPlotData.amsfile)
            ReadPlotData.not_read_ams=False
            ReadPlotData.ams_file_present=True
        return
    ############################################################################
    # @brief Read the data for the \f$\mathbf{S}(\mathbf{q},\omega)\f$ and postprocess it for plotting
    # @details Reads the data and it post-process it to ensure that the spin-spin
    # correlation function can be properly ploted.
    # @author Anders Bergman and Jonathan Chico
    ############################################################################
    def read_sqw(self,filename):
        import pandas as pd
        import numpy as np
        from scipy import signal
        sqwa = pd.read_csv(filename,header=None,delim_whitespace=True).values
        qd=int(sqwa[sqwa.shape[0]-1,0])
        ed=int(sqwa[sqwa.shape[0]-1,4])
        sqw_data=[]
        sigma=1.50
        gauss=signal.gaussian(ed,std=sigma)
        ax_limits=[min(sqwa[:,0]),max(sqwa[:,0]),min(sqwa[:,4])*ReadPlotData.hf_scale,max(sqwa[:,4])*ReadPlotData.hf_scale]
        #-----------------------------------------------------------------------
        # Perform a convolution with a windowing function for each q-point
        #-----------------------------------------------------------------------
        for iq in range(0,qd):
            indx=np.where(sqwa[:,0]==(iq+1))
            for ii in range(0,4):
                sqwa[indx[0],ii+5]=signal.convolve(sqwa[indx[0],ii+5],gauss,mode='same')
        #-----------------------------------------------------------------------
        # Find the peaks and normalize the data
        #-----------------------------------------------------------------------
        for ii in range(5,len(sqwa[0])):
            sqw=np.transpose((np.reshape(sqwa[:,ii],(qd,ed))[:,0:int(ed/2)]))
            sqw_peaks=np.argmax(sqw,axis=0)
            normMat=np.diag(1.0/np.amax(sqw,axis=0))
            sqw=np.matmul(sqw,normMat)
            sqw_data.append(sqw)
        sqw_labels=[r'$S_x(q,\omega)$ [meV]',r'$S_y(q,\omega)$ [meV]',r'$S_z(q,\omega)$ [meV]',r'$S^2(q,\omega)$ [meV]']
        return sqw_data,sqw_labels,ax_limits
    ############################################################################
    # @brief Function to read the AMS file
    # @details It sets up the data such that one can plot the different branches
    # of the AMS.
    # @author Jonathan Chico
    ############################################################################
    def read_ams(self,filename):
        import pandas as pd
        data=pd.read_csv(filename,header=None,delim_whitespace=True,skiprows=2).values
        ams_data_x=[]
        ams_data_y=[]
        ams_label=[]
        for ii in range(1,len(data[0])-1):
            ams_data_x.append(data[:,-1])
            ams_data_y.append(data[:,ii])
            if ii==1:
                ams_label.append('AMS')
            else:
                ams_label.append(None)
        ams_ax_label=[r'q',r'Energy [meV]']
        return ams_data_x,ams_data_y,ams_ax_label,ams_label
    ############################################################################
    # @brief Wrapper to read the trajectory files
    # @details Postprocess the data to create an array of arrays to ensure that
    # one can plot several trajectories at the same time
    # @author Jonathan Chico
    ############################################################################
    def read_trajectories(self,filename):
        import pandas as pd
        traj_label=[]
        traj_data_x=[]
        traj_data_y=[]
        traj_data_z=[]
        for ii in range(len(filename)):
            curr_filename=filename[ii]
            ind=curr_filename.find('.')
            ind=ind+10 # To correct for the simid
            ind_2=curr_filename[ind:].find('.')
            ntraj=curr_filename[ind:ind+ind_2]
            data=pd.read_csv(curr_filename,header=None,delim_whitespace=True,usecols=[2,3,4]).values
            traj_label.append('Atom='+ntraj)
            traj_data_x.append(data[:,0])
            traj_data_y.append(data[:,1])
            traj_data_z.append(data[:,2])
        return traj_label,traj_data_x,traj_data_y,traj_data_z

    #################################################################################
    # General function to read-in data for 2D plots that have headers, e.g. averages, energies
    #################################################################################
    def read_gen_plot_data(self,filename):
        import pandas as pd
        data_df=pd.read_csv(filename,header=0,delim_whitespace=True)
        data=data_df.values
        t_data=[]
        itr_data=[]
        data_labels=[]
        data[:,0]=data[:,0]*ReadPlotData.timestep
        for ii in range(1,len(data[0])):
            itr_data.append(data[:,0])
            t_data.append(data[:,ii])
            data_labels.append(str('$'+data_df.columns[ii]+'$'))
        return t_data,itr_data,data_labels

