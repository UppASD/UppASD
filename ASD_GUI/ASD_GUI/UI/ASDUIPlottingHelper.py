"""
ASDUIPlottingHelper.py

This module provides helper functions for initializing and managing the plotting UI in the ASD GUI.
It includes functions for setting up the UI, selecting data to plot, creating checkboxes for AMS
branches, and saving figures.

Functions:
- InitPlotUI(window): Initializes the plot UI.
- PlottingSelector(window): Handles the selection and plotting of different data types based on the
    sender action.
- set_ams_checkboxes(window): Initializes and sets up checkboxes for AMS branches in the UI.
- PlottingWrapper(window): Wrapper function that takes care of plotting the selected plot, allowing
    the user to choose between different types of plots.
- SaveFig(self): Saves the current figure to a file with specified DPI.
"""
# pylint: disable=invalid-name, no-name-in-module, no-member
from PyQt6.QtCore import QSignalBlocker
from PyQt6.QtWidgets import QCheckBox, QFileDialog

# from matplotlib.backends.backend_qt5agg import FigureCanvas

from ASD_GUI.UI import ASDInputWindows


##########################################################################
# Initialization of some of the UI properties for 2D plots
##########################################################################
def InitPlotUI(window):
    """
    Initializes the plot UI by setting file names and disabling certain UI elements.
    """
    window.plotfile_names[0] = window.ASDPlotData.yamlfile
    window.plotfile_names[1] = window.ASDPlotData.amsfile
    window.plotfile_names[2] = window.ASDPlotData.sqwfile
    window.plotfile_names[3] = window.ASDPlotData.averages
    window.plotfile_names[4] = window.ASDPlotData.trajectory
    window.plotfile_names[5] = window.ASDPlotData.totenergy
    window.plotfile_names[6] = window.ASDPlotData.qfile
    window.SqwProjBox.setEnabled(False)
    window.SqwColorMapSelect.setEnabled(False)
    window.AveOpts.setEnabled(False)
    window.EneOpts.setEnabled(False)
    window.AMSDisplayOpts.setVisible(False)
    return


##########################################################################
# Function to select the appropriate data to plot
##########################################################################
def PlottingSelector(window):
    """
    Handles the selection and plotting of different data types based on the sender action.
    """
    # -----------------------------------------------------------------------
    # Plot the spin-spin correlation function
    # -----------------------------------------------------------------------
    if window.sender() == window.actionS_q_w:
        window.plotting_mode = "correlation"
        window.MatToolBox.setCurrentIndex(0)
        window.PlotStacked.setCurrentIndex(0)
        if not window.ASDPlotData.not_read_ams:
            window.AMSDispCheckBox.setChecked(True)
            QSignalBlocker(window.AMSDispCheckBox)
        if not window.ASDPlotData.not_read_sqw:
            window.SqwDispCheckBox.setChecked(True)
            QSignalBlocker(window.SqwDispCheckBox)
    # -----------------------------------------------------------------------
    # Plot the averages
    # -----------------------------------------------------------------------
    if window.sender() == window.actionAverages:
        window.plotting_mode = "averages"
        window.MatToolBox.setCurrentIndex(2)
        window.PlotStacked.setCurrentIndex(0)
    # -----------------------------------------------------------------------
    # Plot the energies
    # -----------------------------------------------------------------------
    if window.sender() == window.actionTotEnergy:
        window.plotting_mode = "energy"
        window.MatToolBox.setCurrentIndex(3)
        window.PlotStacked.setCurrentIndex(0)
        # -------------------------------------------------------------------
        # Check if the 2D axis exists if it does not create it
        # -------------------------------------------------------------------
    if window.sender() == window.actionTrajectory:
        window.plotting_mode = "trajectory"
        window.MatToolBox.setCurrentIndex(1)
        window.PlotStacked.setCurrentIndex(1)
    window.Plotting_Figure.canvas.draw()
    InitPlotUI(window)
    window.ASDPlotData.PlotReadingWrapper(window.plotfile_names, window)
    if not window.ASDPlotData.not_read_ams:
        window.ams_data_x = window.ASDPlotData.ams_data_x
        window.ams_data_y = window.ASDPlotData.ams_data_y
        window.ams_label = window.ASDPlotData.ams_label
    InitPlotUI(window)
    set_ams_checkboxes(window)
    return


##########################################################################
# @brief Function for the creation of checkboxes for the ams display
# @details This should allow for the dynamical creation of checkboxes for each
# branch in the ams. It also connects it to a function that prunes the data
# so that it can be selectively plotted.
# @author Jonathan Chico
##########################################################################
def set_ams_checkboxes(window):
    """
    Initializes and sets up checkboxes for AMS branches in the UI.
    """
    window.AMSCheckboxes = dict()
    for ii in reversed(range(window.AMSDisplayLayout.count())):
        window.AMSDisplayLayout.itemAt(ii).widget().setParent(None)
    if window.ASDPlotData.ams_file_present:
        # -------------------------------------------------------------------
        # Create checkboxes for the AMS branches
        # -------------------------------------------------------------------
        for ii in range(0, len(window.ASDPlotData.ams_data_y)):
            name = f"ams_branch_{ii}"
            checkbox = QCheckBox()
            checkbox.setObjectName(name)
            checkbox.setText(f"Display Branch {ii + 1: 4d}")
            checkbox.setChecked(True)
            checkbox.toggled.connect(window.AMS_PrunePlot)
            window.AMSDisplayLayout.addWidget(checkbox)
            window.AMSCheckboxes[name] = checkbox
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
def PlottingWrapper(window):
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

    # -----------------------------------------------------------------------
    # Plotting the spin-spin correlation function
    # -----------------------------------------------------------------------
    if (
        window.sender() == window.AMSDispCheckBox
        or window.sender() == window.SqwDispCheckBox
    ):
        window.plotting_mode = "correlation"
        if window.ASDPlotData.not_read_sqw or window.ASDPlotData.not_read_ams:
            window.ASDPlotData.PlotReadingWrapper(window.plotfile_names, window)
            if window.AMSDispCheckBox.isChecked():
                window.ams_data_x = window.ASDPlotData.ams_data_x
                window.ams_data_y = window.ASDPlotData.ams_data_y
                window.ams_label = window.ASDPlotData.ams_label
                window.set_ams_checkboxes()

    if window.plotting_mode == "correlation":
        # -------------------------------------------------------------------
        # Perform the actual plotting
        # -------------------------------------------------------------------
        window.SqwProjBox.setEnabled(True)
        window.SqwColorMapSelect.setEnabled(True)
        window.SqwDisplayOpts.setEnabled(True)
        # -------------------------------------------------------------------
        # Plotting the S(q,w)
        # -------------------------------------------------------------------
        if (
            window.SqwDispCheckBox.isChecked()
            and not window.AMSDispCheckBox.isChecked()
        ):
            if window.ASDPlotData.sqw_file_present:
                window.ASDCorrelationPlots.Sqw_Plot(
                    window.Plotting_ax,
                    window.ASDPlotData.sqw_data,
                    window.SQW_proj_indx,
                    window.ASDPlotData.sqw_labels,
                    window.plot2D_cmap_indx,
                    window.ASDPlotData.ax_limits,
                    window.ASDPlotData.q_labels,
                    window.ASDPlotData.q_idx,
                )
                window.AMSDisplayOpts.setVisible(False)
                window.AMSDisplayOpts.setEnabled(False)
            else:
                window.sqw_Error_Window = ASDInputWindows.ErrorWindow()
                window.sqw_Error_Window.FunMsg.setText(
                    "I'm sorry, Dave. I'm afraid I can't do that."
                )
                window.sqw_Error_Window.ErrorMsg.setText(
                    "Error: No 'sqw.*.out' file."
                )
                window.sqw_Error_Window.show()
                window.SqwDispCheckBox.setChecked(False)
                print("No 'sqw.*.out' file.")
        # -------------------------------------------------------------------
        # Plotting the AMS
        # -------------------------------------------------------------------
        elif (
            window.AMSDispCheckBox.isChecked()
            and not window.SqwDispCheckBox.isChecked()
        ):
            if window.ASDPlotData.ams_file_present:
                window.ASDPlots2D.LinePlot(
                    window.Plotting_ax,
                    window.ams_data_x,
                    window.ams_data_y,
                    window.ams_label,
                    window.ASDPlotData.ams_ax_label,
                    tick_labels=window.ASDPlotData.q_labels,
                    tick_idx=window.ASDPlotData.q_idx,
                )
                window.AMSDisplayOpts.setVisible(True)
                window.AMSDisplayOpts.setEnabled(True)
            else:
                window.ams_Error_Window = ASDInputWindows.ErrorWindow()
                window.ams_Error_Window.FunMsg.setText(
                    "I'm sorry, Dave. I'm afraid I can't do that."
                )
                window.ams_Error_Window.ErrorMsg.setText(
                    "Error: No 'ams.*.out' file."
                )
                window.ams_Error_Window.show()
                window.AMSDispCheckBox.setChecked(False)
                print("No 'ams.*.out' file.")
        # -------------------------------------------------------------------
        # Plotting the S(q,w) and the AMS
        # -------------------------------------------------------------------
        if window.SqwDispCheckBox.isChecked() and window.AMSDispCheckBox.isChecked():
            if (
                window.ASDPlotData.sqw_file_present
                and window.ASDPlotData.ams_file_present
            ):
                window.ASDCorrelationPlots.AMS_Sqw_Plot(
                    window.Plotting_ax,
                    window.ASDPlotData.sqw_data,
                    window.SQW_proj_indx,
                    window.ASDPlotData.sqw_labels,
                    window.ams_data_x,
                    window.ams_data_y,
                    window.ASDPlotData.hf_scale,
                    window.plot2D_cmap_indx,
                    window.ASDPlotData.ax_limits,
                    window.ASDPlotData.q_labels,
                    window.ASDPlotData.q_idx,
                )
                window.AMSDisplayOpts.setVisible(True)
                window.AMSDisplayOpts.setEnabled(True)
            else:
                window.ams_sqw_Error_Window = ASDInputWindows.ErrorWindow()
                window.ams_sqw_Error_Window.FunMsg.setText(
                    "I'm sorry, Dave. I'm afraid I can't do that."
                )
                window.ams_sqw_Error_Window.ErrorMsg.setText(
                    "Error: No 'ams.*.out' or 'sqw.*.out' file."
                )
                window.ams_sqw_Error_Window.show()
                print("No 'ams.*.out' or 'sqw.*.out' file.")
                window.SqwProjBox.setEnabled(False)
                window.SqwColorMapSelect.setEnabled(False)
                window.AMSDispCheckBox.setChecked(False)
                window.SqwDispCheckBox.setChecked(False)
    # -----------------------------------------------------------------------
    # Plotting the average magnetization
    # -----------------------------------------------------------------------
    if window.plotting_mode == "averages":
        window.AveOpts.setEnabled(True)
        if window.ASDPlotData.ave_file_present:
            curr_data_x = []
            curr_data_y = []
            curr_labels = []
            for ii, _ in enumerate(window.MagDirIndx):
                curr_data_x.append(window.ASDPlotData.mitr_data[window.MagDirIndx[ii]])
                curr_data_y.append(window.ASDPlotData.mag_data[window.MagDirIndx[ii]])
                curr_labels.append(window.ASDPlotData.mag_labels[window.MagDirIndx[ii]])
            if len(window.MagDirIndx) > 0:
                window.ASDPlots2D.LinePlot(
                    window.Plotting_ax,
                    curr_data_x,
                    curr_data_y,
                    curr_labels,
                    window.ASDPlotData.mag_axes,
                )
            else:
                print("Select at least one direction to plot")
        else:
            window.ave_Error_Window = ASDInputWindows.ErrorWindow()
            window.ave_Error_Window.FunMsg.setText(
                "I'm sorry, Dave. I'm afraid I can't do that."
            )
            window.ave_Error_Window.ErrorMsg.setText(
                "Error: No 'averages.*.out' file."
            )
            window.ave_Error_Window.show()
            print("No 'averages.*.out' file.")
    # -----------------------------------------------------------------------
    # Plotting the total energy
    # -----------------------------------------------------------------------
    if window.plotting_mode == "energy":
        window.EneOpts.setEnabled(True)
        if window.ASDPlotData.ene_file_present:
            curr_data_x = []
            curr_data_y = []
            curr_labels = []
            for ii, index in enumerate(window.EneIndx):
                curr_data_x.append(window.ASDPlotData.eitr_data[index])
                curr_data_y.append(window.ASDPlotData.ene_data[index])
                curr_labels.append(window.ASDPlotData.ene_labels[index])
            if len(window.EneIndx) > 0:
                window.ASDPlots2D.LinePlot(
                    window.Plotting_ax,
                    curr_data_x,
                    curr_data_y,
                    curr_labels,
                    window.ASDPlotData.ene_axes,
                )
            else:
                print("Select at least one energy component to plot")
        else:
            window.ene_Error_Window = ASDInputWindows.ErrorWindow()
            window.ene_Error_Window.FunMsg.setText(
                "I'm sorry, Dave. I'm afraid I can't do that."
            )
            window.ene_Error_Window.ErrorMsg.setText(
                "Error: No 'totenergy.*.out' file."
            )
            window.ene_Error_Window.show()
            print("No 'totenergy.*.out' file.")
    # -----------------------------------------------------------------------
    # Plotting the single magnetic moment trajectories
    # -----------------------------------------------------------------------
    if window.plotting_mode == "trajectory":
        if window.ASDPlotData.trajectory_file_present:
            window.ASDPlots2D.TrajPlot(
                window.Plotting_ax3D,
                window.ASDPlotData.traj_data_x,
                window.ASDPlotData.traj_data_y,
                window.ASDPlotData.traj_data_z,
                window.ASDPlotData.traj_label,
            )
        else:
            window.traj_Error_Window = ASDInputWindows.ErrorWindow()
            window.traj_Error_Window.FunMsg.setText(
                "I'm sorry, Dave. I'm afraid I can't do that."
            )
            window.traj_Error_Window.ErrorMsg.setText(
                "Error: No 'trajectory.*.out' file."
            )
            window.traj_Error_Window.show()
            print("No 'trajectory.*.out' file.")
    window.Plotting_Figure.canvas.draw()
    window.Plotting_Figure.canvas.flush_events()
    return


##########################################################################
# @brief Function to save the current figure to file
# @author Jonathan Chico
##########################################################################
def SaveFig(self):
    """
    Save the current figure to a file with specified DPI.
    """
    fig_name, _ = QFileDialog.getSaveFileName(self, "Save File")
    if len(self.InpFigDPI.text()) > 0:
        dpi = int(self.InpFigDPI.text())
    else:
        dpi = 800
    if self.plotting_mode != "trajectory":
        fig_plot = self.Plotting_canvas
        fig_plot.print_figure(fig_name, dpi=dpi)
    else:
        fig_plot = self.Plotting_canvas3D
        fig_plot.print_figure(fig_name, dpi=dpi)
    return
