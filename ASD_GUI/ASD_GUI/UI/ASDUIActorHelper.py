"""
Wrapper function that takes care of adding the necessary actors and the
options for the different types of visualizations.

This class controls the visualization of:
    - Restartfiles
    - Momentsfiles
    - Energy
    - Exchange neighbours
    - DM neighbours

Attributes
window : object
    The main window object which contains various UI elements and data structures.

Methods
-------
AddActors(window)
    Wrapper function that takes care of adding the necessary actors and the
    options for the different types of visualizations.
"""
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
# pylint: disable=invalid-name, no-name-in-module, no-member
import io
from PyQt6.QtWidgets import QLabel

from ASD_GUI.UI import ASDUIInitHelper
from ASD_GUI.VTK_Viz import ASDVTKEneActors, ASDVTKNeighActors


def AddActors(window, viz_type=None):
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
    try:
        window.ASDGenActors.scalar_bar_widget
    except AttributeError:
        pass
    else:
        window.ASDGenActors.reset_GenActors()
    window.ren.RemoveAllViewProps()
    window.ren.ResetCamera()
    ASDUIInitHelper.InitUI(window)
    # Storing settings virtually if they exist
    if window.runtime_settings is not None:
        window.UISettings.gather_dicts(window)
        window.UISettings.write_to_yaml(window.runtime_settings)
        # window.UISettings.write_to_yaml("test.yaml")
    # -----------------------------------------------------------------------
    # This takes care of setting up the options for the visualization of the
    # magnetic moments obtained from the restart file
    # -----------------------------------------------------------------------
    if window.sender() == window.actionMagnetization or viz_type == "M":
        # -------------------------------------------------------------------
        # Call the Moments class
        # -------------------------------------------------------------------
        # self.MomActors = ASDVTKMomActors.ASDMomActors()
        print("Magnetization mode chosen")
        window.viz_type = "M"
        window.mode = 1
        window.current_time = 0
        window.MagMainGroup.setEnabled(True)
        window.VizToolBox.setCurrentIndex(0)
        window.menuMagnetisation_Opts.setEnabled(True)
        window.actionDisplayMagDens.setEnabled(True)
        window.actionX_ProjMagDens.setEnabled(True)
        window.actionY_ProjMagDens.setEnabled(True)
        window.actionZ_ProjMagDens.setEnabled(True)
        window.PlayButton.setEnabled(True)
        window.PauseButton.setEnabled(True)
        window.nextButton.setEnabled(True)
        window.previousButton.setEnabled(True)
        # -------------------------------------------------------------------
        # Add the data structures with regards to reading the data
        # -------------------------------------------------------------------
        window.ASDdata.ReadingWrapper(
            mode=window.mode,
            viz_type=window.viz_type,
            file_names=window.file_names,
            window=window,
        )
        if not window.ASDdata.error_trap:
            window.MomActors.Add_MomActors(
                ren=window.ren,
                renWin=window.renWin,
                iren=window.iren,
                ASDdata=window.ASDdata,
                window=window,
            )
            window.ASDVizOpt.update_dock_info(
                current_actors=window.MomActors, window=window
            )
            window.current_actor = window.MomActors
            # ---------------------------------------------------------------
            # Setup several global variables
            # ---------------------------------------------------------------
            window.ASDColor.lut = window.MomActors.lut
            # ---------------------------------------------------------------
            # Add the general widgets such as the scalar bar and the axes
            # ---------------------------------------------------------------
            print("Adding the general actors")
            window.ASDGenActors.Add_GenActors(
                iren=window.iren,
                renWin=window.renWin,
                method=window.MomActors.MagDensMethod,
                lut=window.ASDColor.lut,
                ren=window.ren,
                window=window,
                current_Actors=window.MomActors,
                flag2D=window.ASDdata.flag2D,
            )
            # ---------------------------------------------------------------
            # Update the UI
            # ---------------------------------------------------------------
            if window.ASDdata.cluster_flag:
                window.ClusBox.setVisible(True)
                window.ClusBox.setChecked(True)
                window.ASDGenActors.Add_ClusterActors(
                    ASDdata=window.ASDdata,
                    iren=window.iren,
                    renWin=window.renWin,
                    ren=window.ren,
                )
            if window.ASDdata.kmc_flag:
                window.KMCCheck.setVisible(True)
            # ---------------------------------------------------------------
            # Print the visualization message
            # ---------------------------------------------------------------
            print("Visualization of magnetic moments mode chosen")
            window.current_time = window.current_time + 1
            # ---------------------------------------------------------------
            # Re-initialize previous setup (or save current)
            # ---------------------------------------------------------------
            if window.runtime_settings is None:
                print("No settings file found, regenerating the settings")
                window.UISettings.gather_dicts(window)
                window.runtime_settings = io.StringIO()
                window.UISettings.write_to_yaml(window.runtime_settings)
            else:
                # window.UISettings.read_from_yaml('test.yaml')
                print("Settings file found, reading the settings")
                window.UISettings.read_from_yaml(window.runtime_settings)
                window.UISettings.restore_from_settings(window)
                window.renWin.Render()
    # -----------------------------------------------------------------------
    # This takes care of setting up the options for the Neighbour visualization
    # -----------------------------------------------------------------------
    if window.sender() == window.actionNeighbours or viz_type == "N":
        # -------------------------------------------------------------------
        # Call the Neighbour class
        # -------------------------------------------------------------------
        window.NeighActors = ASDVTKNeighActors.ASDNeighActors()
        window.mode = 1
        window.viz_type = "N"
        window.NeighMainBox.setEnabled(True)
        window.VizToolBox.setCurrentIndex(2)
        window.PlayButton.setEnabled(False)
        window.PauseButton.setEnabled(False)
        window.nextButton.setEnabled(False)
        window.previousButton.setEnabled(False)
        # -------------------------------------------------------------------
        # Add the data structures with regards to reading the data
        # -------------------------------------------------------------------
        window.ASDdata.ReadingWrapper(
            mode=window.mode,
            viz_type=window.viz_type,
            file_names=window.file_names,
            window=window,
        )
        if not window.ASDdata.error_trap:
            window.NeighActors.Add_NeighActors(
                ren=window.ren,
                renWin=window.renWin,
                iren=window.iren,
                ASDdata=window.ASDdata,
                mode=window.mode,
            )
            # ---------------------------------------------------------------
            # Set several global variables
            # ---------------------------------------------------------------
            window.ASDColor.lut = window.NeighActors.lut
            # ---------------------------------------------------------------
            # Add the general widgets such as the scalar bar and the axes
            # ---------------------------------------------------------------
            window.ASDGenActors.Add_GenActors(
                iren=window.iren,
                renWin=window.renWin,
                method=window.NeighActors.NeighGlyph3D,
                lut=window.ASDColor.lut,
                ren=window.ren,
                window=window,
                current_Actors=window.NeighActors,
                flag2D=True,
            )
            # ---------------------------------------------------------------
            # Update the labels for the neighbour mode
            # ---------------------------------------------------------------
            window.NeighSelectSlider.setMaximum(window.NeighActors.SLMax)
            window.NeighSelectSlider.setMinimum(1)
            window.NeighNumberLabel.setText(
                f"Number of neighbours = {window.NeighActors.NumNeigh: 4d}"
            )
            window.NeighValidator.setRange(1, window.ASDdata.nrAtoms)
            window.NeighSelectLineEdit.setValidator(window.NeighValidator)
            window.ASDVizOpt.update_dock_info(
                current_actors=window.NeighActors, window=window
            )
            window.current_actor = window.NeighActors
            # ---------------------------------------------------------------
            # Update the UI
            # ---------------------------------------------------------------
            window.NeighTypesLabels = dict()
            for ii in range(0, window.ASDdata.num_types_total):
                name = f"label_neigh_{ii}"
                label = QLabel()
                label.setObjectName(name)
                label.setText(f"Num. Neighbours Type {ii + 1: 4d} = {0: 4d}")
                window.NeighInfoLayout.addWidget(label)
                window.NeighTypesLabels[name] = label
            for ii in range(0, window.ASDdata.num_types):
                name = f"label_neigh_{int(window.ASDdata.types[ii] - 1)}"
                window.NeighTypesLabels[name].setText(
                    f"Num. Neighbours Type {ii + 1: 4d} = {window.ASDdata.types_counters[ii]: 4d}"
                )
            # ---------------------------------------------------------------
            # Visualize the embedded cluster into the system
            # ---------------------------------------------------------------
            if window.ASDdata.cluster_flag:
                window.ClusBox.setVisible(True)
                window.ClusBox.setChecked(True)
                window.ASDGenActors.Add_ClusterActors(
                    ASDdata=window.ASDdata,
                    iren=window.iren,
                    renWin=window.renWin,
                    ren=window.ren,
                )
            # ---------------------------------------------------------------
            # Print the visualization message
            # ---------------------------------------------------------------
            print("Visualization of the neighbour map mode chosen")
            print("Viewing the struct file")
    # -----------------------------------------------------------------------
    # This takes care of setting up the options for the DM Neighbour visualization
    # -----------------------------------------------------------------------
    if window.sender() == window.actionDM_Neigh or viz_type == "D":
        # -------------------------------------------------------------------
        # Call the Neighbour class
        # -------------------------------------------------------------------
        window.NeighActors = ASDVTKNeighActors.ASDNeighActors()
        window.mode = 2
        window.viz_type = "N"
        window.NeighMainBox.setEnabled(True)
        window.VizToolBox.setCurrentIndex(2)
        window.PlayButton.setEnabled(False)
        window.PauseButton.setEnabled(False)
        window.nextButton.setEnabled(False)
        window.previousButton.setEnabled(False)
        # -------------------------------------------------------------------
        # Add the data structures with regards to reading the data
        # -------------------------------------------------------------------
        window.ASDdata.ReadingWrapper(
            mode=window.mode,
            viz_type=window.viz_type,
            file_names=window.file_names,
            window=window,
        )
        if not window.ASDdata.error_trap:
            window.NeighActors.Add_NeighActors(
                ren=window.ren,
                renWin=window.renWin,
                iren=window.iren,
                ASDdata=window.ASDdata,
                mode=window.mode,
            )
            # ---------------------------------------------------------------
            # Set several global variables
            # ---------------------------------------------------------------
            window.ASDColor.lut = window.NeighActors.lut
            # ---------------------------------------------------------------
            # Add the general widgets such as the scalar bar and the axes
            # ---------------------------------------------------------------
            window.ASDGenActors.Add_GenActors(
                iren=window.iren,
                renWin=window.renWin,
                method=window.NeighActors.NeighGlyph3D,
                lut=window.ASDColor.lut,
                ren=window.ren,
                window=window,
                current_Actors=window.NeighActors,
                flag2D=True,
            )
            # ---------------------------------------------------------------
            # Update the labels for the neighbour mode
            # ---------------------------------------------------------------
            window.NeighSelectSlider.setMaximum(window.NeighActors.SLMax)
            window.NeighSelectSlider.setMinimum(1)
            window.NeighNumberLabel.setText(
                f"Number of neighbours = {window.NeighActors.NumNeigh: 4d}"
            )
            window.ASDVizOpt.update_dock_info(
                current_actors=window.NeighActors, window=window
            )
            window.current_actor = window.NeighActors
            # ---------------------------------------------------------------
            # Update the UI
            # ---------------------------------------------------------------
            window.NeighTypesLabels = dict()
            for ii in range(0, window.ASDdata.num_types_total):
                name = f"label_neigh_{ii}"
                label = QLabel()
                label.setObjectName(name)
                label.setText(f"Num. Neighbours Type {ii + 1: 4d} = {0: 4d}")
                window.NeighInfoLayout.addWidget(label)
                window.NeighTypesLabels[name] = label
            for ii in range(0, window.ASDdata.num_types):
                name = f"label_neigh_{int(window.ASDdata.types[ii] - 1)}"
                window.NeighTypesLabels[name].setText(
                    f"Num. Neighbours Type {ii + 1: 4d} = {window.ASDdata.types_counters[ii]: 4d}"
                )
                # -----------------------------------------------------------
                # Visualize the embedded cluster into the system
                # -----------------------------------------------------------
            if window.ASDdata.cluster_flag:
                window.ClusBox.setVisible(True)
                window.ClusBox.setChecked(True)
                window.ASDGenActors.Add_ClusterActors(
                    ASDdata=window.ASDdata,
                    iren=window.iren,
                    renWin=window.renWin,
                    ren=window.ren,
                )
            # ---------------------------------------------------------------
            # Print the visualization message
            # ---------------------------------------------------------------
            print("Visualization of the neighbour map mode chosen")
            print("Viewing the struct file")
    # -----------------------------------------------------------------------
    # This takes care of setting up the options for the Energy visualization
    # -----------------------------------------------------------------------
    if window.sender() == window.actionEnergy or viz_type == "E":
        window.EneActors = ASDVTKEneActors.ASDEneActors()
        window.viz_type = "E"
        window.mode = 1
        window.current_time = 0
        window.VizToolBox.setCurrentIndex(1)
        window.EneMainBox.setEnabled(True)
        window.PlayButton.setEnabled(True)
        window.PauseButton.setEnabled(True)
        window.nextButton.setEnabled(True)
        window.previousButton.setEnabled(True)
        # -------------------------------------------------------------------
        # Add the data structures with regards to reading the data
        # -------------------------------------------------------------------
        window.ASDdata.ReadingWrapper(
            mode=window.mode,
            viz_type=window.viz_type,
            file_names=window.file_names,
            window=window,
        )
        if not window.ASDdata.error_trap:
            window.EneActors.Add_EneActors(
                ren=window.ren,
                renWin=window.renWin,
                iren=window.iren,
                ASDdata=window.ASDdata,
            )
            window.ASDVizOpt.update_dock_info(
                current_actors=window.EneActors, window=window
            )
            window.current_actor = window.EneActors
            # ---------------------------------------------------------------
            # Setup several global variables
            # ---------------------------------------------------------------
            window.ASDColor.lut = window.EneActors.lut
            # ---------------------------------------------------------------
            # Add the general widgets such as the scalar bar and the axes
            # ---------------------------------------------------------------
            window.ASDGenActors.Add_GenActors(
                iren=window.iren,
                renWin=window.renWin,
                method=window.EneActors.EneDensMethod,
                lut=window.ASDColor.lut,
                ren=window.ren,
                window=window,
                current_Actors=window.EneActors,
                flag2D=window.ASDdata.flag2D,
            )
            # ---------------------------------------------------------------
            # Update the UI
            # ---------------------------------------------------------------
            if window.ASDdata.cluster_flag:
                window.ClusBox.setVisible(True)
                window.ClusBox.setChecked(True)
                window.ASDGenActors.Add_ClusterActors(
                    ASDdata=window.ASDdata,
                    iren=window.iren,
                    renWin=window.renWin,
                    ren=window.ren,
                )
            # ---------------------------------------------------------------
            # Print the visualization message
            # ---------------------------------------------------------------
            print("Visualization of the energy mode chosen")
            print("Viewing the localenergy file")
            
    return
