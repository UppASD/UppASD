"""@package ASDMenuToolbar
Set of auxiliary functions for the Menu and Toolbar options of the ASD_GUI.

Contains the functions that setup the toolbars and many of the connectors for the
UI actions. Also contains the wrapper function for the UI Update, thus allowing for
the synchronization of several options which yield a same outcome in the UIself.

Author
----------
Jonathan Chico
"""
# pylint: disable=invalid-name, no-name-in-module, no-member

from PyQt6.QtGui import QDoubleValidator
from PyQt6.QtWidgets import QProgressBar, QLabel, QStyle, QToolButton
from PyQt6.QtCore import QSignalBlocker


##########################################################################
# @brief Function defining the toolbar and actions for the VTK backend.
# @details This contains several buttons handling the different visualization
# modes available in the VTK API.
# @author Jonathan Chico
##########################################################################
def VTK_Menu_and_Toolbar_Setup(window):
    """Functions defining the toolbar and actions associated to it for the VTK backend.
    This contains several buttons and actions handling the different types of visualization
    modes available in the VTK API.

    Args
    ----------
        window: QMainWindow object where the toolbar is located.

    Author
    ----------
    Jonathan Chico
    """

    window.VTKToolBar.setFixedHeight(24)
    window.ASD_VTK_Layout.insertWidget(0, window.VTKToolBar)
    # ---------------------------------------------------------------------------
    # Button to exit the program
    # ---------------------------------------------------------------------------
    window.VTKExitButton = QToolButton()
    window.VTKExitButton.setText("Exit")
    window.VTKExitButton.setToolTip("Exit the UppASD visualizer")
    window.VTKExitButton.setWhatsThis("Exit the UppASD visualizer")
    window.VTKExitButton.setStatusTip("Exit the UppASD visualizer")
    window.VTKExitButton.clicked.connect(window.close)
    # ---------------------------------------------------------------------------
    # button to take a snapshot of the visualization
    # ---------------------------------------------------------------------------
    window.actionTake_Snapshot = QToolButton()
    window.actionTake_Snapshot.setText("Snapshot")
    window.actionTake_Snapshot.setToolTip("Take a snapshot of the current image")
    window.actionTake_Snapshot.setWhatsThis("Take a snapshot of the current image")
    window.actionTake_Snapshot.setStatusTip("Take a snapshot of the current image")
    # ---------------------------------------------------------------------------
    # Play button to visualize the restartfile
    # ---------------------------------------------------------------------------
    window.actionMagnetization = QToolButton()
    window.actionMagnetization.setText("Magnetization")
    window.actionMagnetization.setToolTip(
        "Visualize a magnetic configuration (restart/moment) from UppASD"
    )
    window.actionMagnetization.setWhatsThis(
        "Visualize a magnetic configuration (restart/moment) from UppASD"
    )
    window.actionMagnetization.setStatusTip(
        "Visualize a magnetic configuration (restart/moment) from UppASD"
    )
    # ---------------------------------------------------------------------------
    # Play button to visualize the momentsfile
    # ---------------------------------------------------------------------------
    window.actionMoments = QToolButton()
    window.actionMoments.setText("Moments")
    window.actionMoments.setToolTip("Visualize the momentsfile from UppASD")
    window.actionMoments.setWhatsThis("Visualize the momentsfile from UppASD")
    window.actionMoments.setStatusTip("Visualize the momentsfile from UppASD")
    # ---------------------------------------------------------------------------
    # Play button to visualize the energy
    # ---------------------------------------------------------------------------
    window.actionEnergy = QToolButton()
    window.actionEnergy.setText("Energy")
    window.actionEnergy.setToolTip("Visualize the site dependent energy from UppASD")
    window.actionEnergy.setWhatsThis("Visualize the site dependent energy from UppASD")
    window.actionEnergy.setStatusTip("Visualize the site dependent energy from UppASD")
    # ---------------------------------------------------------------------------
    # Play button to visualize the exchange neighbours
    # ---------------------------------------------------------------------------
    window.actionNeighbours = QToolButton()
    window.actionNeighbours.setText("Neighbours")
    window.actionNeighbours.setToolTip("Visualize the exchange interaction neighbours")
    window.actionNeighbours.setWhatsThis(
        "Visualize the exchange interaction neighbours"
    )
    window.actionNeighbours.setStatusTip(
        "Visualize the exchange interaction neighbours"
    )
    # ---------------------------------------------------------------------------
    # Play button to visualize the DM neighbours
    # ---------------------------------------------------------------------------
    window.actionDM_Neigh = QToolButton()
    window.actionDM_Neigh.setText("DM Neighbours")
    window.actionDM_Neigh.setToolTip("Visualize the DM interaction neighbours")
    window.actionDM_Neigh.setWhatsThis("Visualize the DM interaction neighbours")
    window.actionDM_Neigh.setStatusTip("Visualize the DM interaction neighbours")
    # ---------------------------------------------------------------------------
    # Play button to start the rendering of the scene
    # ---------------------------------------------------------------------------
    window.PlayButton = QToolButton()
    window.PlayButton.setCheckable(True)
    window.PlayButton.setChecked(False)
    window.PlayButton.setEnabled(False)
    window.PlayButton.setIcon(
        window.style().standardIcon(QStyle.StandardPixmap.SP_MediaPlay)
    )
    window.PlayButton.setToolTip("Start/Pause the visualization of animations")
    window.PlayButton.setWhatsThis("Start/Pause the visualization of animations")
    window.PlayButton.setStatusTip("Start/Pause the visualization of animations")
    window.PlayButton.clicked.connect(window.Playback_control)
    # ---------------------------------------------------------------------------
    # Play button to stop the rendering of the scene
    # ---------------------------------------------------------------------------
    window.PauseButton = QToolButton()
    window.PauseButton.setEnabled(False)
    window.PauseButton.setIcon(
        window.style().standardIcon(QStyle.StandardPixmap.SP_MediaPause)
    )
    window.PauseButton.setToolTip("Pause the visualization of animations")
    window.PauseButton.setWhatsThis("Pause the visualization of animations")
    window.PauseButton.setStatusTip("Pause the visualization of animations")
    window.PauseButton.clicked.connect(window.Playback_control)
    # ---------------------------------------------------------------------------
    # Button to go to the previous snapshot
    # ---------------------------------------------------------------------------
    window.previousButton = QToolButton()
    window.previousButton.setEnabled(False)
    window.previousButton.setIcon(
        window.style().standardIcon(QStyle.StandardPixmap.SP_MediaSkipBackward)
    )
    window.previousButton.setToolTip("Go to the previous image")
    window.previousButton.setWhatsThis("Go to the previous image")
    window.previousButton.setStatusTip("Go to the previous image")
    window.previousButton.clicked.connect(window.Playback_control)
    # ---------------------------------------------------------------------------
    # Button to go to the next snapshot
    # ---------------------------------------------------------------------------
    window.nextButton = QToolButton()
    window.nextButton.setEnabled(False)
    window.nextButton.setIcon(
        window.style().standardIcon(QStyle.StandardPixmap.SP_MediaSkipForward)
    )
    window.nextButton.setToolTip("Go to the next image")
    window.nextButton.setWhatsThis("Go to the next image")
    window.nextButton.setStatusTip("Go to the next image")
    window.nextButton.clicked.connect(window.Playback_control)
    # ---------------------------------------------------------------------------
    # Progress bar showing the progress of the rendering
    # ---------------------------------------------------------------------------
    window.ProgressBar = QProgressBar()
    window.ProgressBar.setValue(0)
    window.ProgressBar.setTextVisible(False)
    window.ProgressBar.setToolTip("Progress bar of the simulation")
    window.ProgressBar.setWhatsThis("Progress bar of the simulation")
    window.ProgressBar.setStatusTip("Progress bar of the simulation")
    window.ProgressLabel = QLabel()
    window.ProgressLabel.setText(f" {int(window.ProgressBar.value())}% ")
    # ---------------------------------------------------------------------------
    # Adding the extra buttons to the toolbar
    # ---------------------------------------------------------------------------
    window.VTKToolBar.addWidget(window.VTKExitButton)
    window.VTKToolBar.addSeparator()
    window.VTKToolBar.addWidget(window.actionTake_Snapshot)
    window.VTKToolBar.addSeparator()
    window.VTKToolBar.addWidget(window.actionMagnetization)
    window.VTKToolBar.addSeparator()
    window.VTKToolBar.addWidget(window.actionEnergy)
    window.VTKToolBar.addSeparator()
    window.VTKToolBar.addWidget(window.actionNeighbours)
    window.VTKToolBar.addWidget(window.actionDM_Neigh)
    window.VTKToolBar.addSeparator()
    window.VTKToolBar.addWidget(window.previousButton)
    window.VTKToolBar.addWidget(window.PlayButton)
    window.VTKToolBar.addWidget(window.PauseButton)
    window.VTKToolBar.addWidget(window.nextButton)
    window.VTKToolBar.addSeparator()
    window.VTKToolBar.addWidget(window.ProgressBar)
    window.VTKToolBar.addWidget(window.ProgressLabel)
    # ---------------------------------------------------------------------------
    # Adding the actions for the Menu
    # ---------------------------------------------------------------------------
    window.actionKMC_File.triggered.connect(window.getFile)
    window.actionStruct_File.triggered.connect(window.getFile)
    window.actionDM_File.triggered.connect(window.getFile)
    window.actionMagnetization_File.triggered.connect(window.getFile)
    window.actionCoordinate_File.triggered.connect(window.getFile)
    window.actionEnergy_File.triggered.connect(window.getFile)
    # ---------------------------------------------------------------------------
    # Adding the actions for the Toolbar
    # ---------------------------------------------------------------------------
    window.actionEnergy.clicked.connect(window.AddActors)
    window.actionMagnetization.clicked.connect(window.AddActors)
    window.actionNeighbours.clicked.connect(window.AddActors)
    window.actionDM_Neigh.clicked.connect(window.AddActors)
    # ---------------------------------------------------------------------------
    # Adding the actions to the colormaps
    # ---------------------------------------------------------------------------
    window.ColorMapBox.activated.connect(window.set_lut_db)
    window.SingleColorBox.toggled.connect(window.toggle_singlecolor)
    window.RGBRedColorSlider.valueChanged.connect(window.set_singlecolor)
    window.RGBRedColorSlider.valueChanged.connect(window.UpdateRenderer)
    window.RGBGreenColorSlider.valueChanged.connect(window.set_singlecolor)
    window.RGBGreenColorSlider.valueChanged.connect(window.UpdateRenderer)
    window.RGBBlueColorSlider.valueChanged.connect(window.set_singlecolor)
    window.RGBBlueColorSlider.valueChanged.connect(window.UpdateRenderer)
    window.BWSinglecolorCheck.clicked.connect(window.toggle_bwSinglecolor)
    window.LinearScale.toggled.connect(window.set_lut)
    window.LogScale.toggled.connect(window.set_lut)
    # ---------------------------------------------------------------------------
    # Adding the actions to the background
    # ---------------------------------------------------------------------------
    window.RGBRedBackgroundSlider.valueChanged.connect(window.set_background)
    window.RGBRedBackgroundSlider.valueChanged.connect(window.UpdateRenderer)
    window.RGBGreenBackgroundSlider.valueChanged.connect(window.set_background)
    window.RGBGreenBackgroundSlider.valueChanged.connect(window.UpdateRenderer)
    window.RGBBlueBackgroundSlider.valueChanged.connect(window.set_background)
    window.RGBBlueBackgroundSlider.valueChanged.connect(window.UpdateRenderer)
    window.BWBackgroundCheck.clicked.connect(window.toggle_bwBackground)
    # ---------------------------------------------------------------------------
    # Adding the actions to the moment options
    # ---------------------------------------------------------------------------
    window.ContourCheck.toggled.connect(window.ASDVizOpt.toggle_contours)
    window.SpinCheck.toggled.connect(window.ASDVizOpt.toggle_directions)
    window.ClusBox.toggled.connect(window.ASDVizOpt.toggle_cluster)
    window.ClusBox.toggled.connect(window.UpdateRenderer)
    window.KMCCheck.toggled.connect(window.ASDVizOpt.toggle_KMC)
    window.KMCCheck.toggled.connect(window.UpdateRenderer)
    # ---------------------------------------------------------------------------
    # Adding the actions to the magnetization density
    # ---------------------------------------------------------------------------
    window.DensBox.toggled.connect(window.ASDVizOpt.toggle_density)
    window.DensX.toggled.connect(window.set_projection)
    window.DensY.toggled.connect(window.set_projection)
    window.DensZ.toggled.connect(window.set_projection)
    window.SpinX.toggled.connect(window.set_projection)
    window.SpinY.toggled.connect(window.set_projection)
    window.SpinZ.toggled.connect(window.set_projection)
    window.actionDisplayMagDens.toggled.connect(window.ASDVizOpt.toggle_density)
    # ---------------------------------------------------------------------------
    # Adding the actions to the atomic spins
    # ---------------------------------------------------------------------------
    window.SpinsBox.toggled.connect(window.ASDVizOpt.toggle_spins)
    window.SpinsBox.toggled.connect(window.UpdateRenderer)
    window.SpinArrowButton.toggled.connect(window.ChangeGlyphs)
    window.SpinCubeButton.toggled.connect(window.ChangeGlyphs)
    window.SpinBarButton.toggled.connect(window.ChangeGlyphs)
    window.SpinSphereButton.toggled.connect(window.ChangeGlyphs)
    window.SpinConeButton.toggled.connect(window.ChangeGlyphs)
    window.SpinSize.valueChanged.connect(window.ASDVizOpt.ChangeSpinsSize)
    window.SpinSize.valueChanged.connect(window.UpdateRenderer)
    window.SpinCenterCheck.toggled.connect(window.ChangeGlyphs)
    # ---------------------------------------------------------------------------
    # Adding shading actions to the spins
    # ---------------------------------------------------------------------------
    window.FlatShadeButton.toggled.connect(window.ChangeShading)
    window.GouraudShadeButton.toggled.connect(window.ChangeShading)
    window.PhongShadeButton.toggled.connect(window.ChangeShading)
    window.PBRShadeButton.toggled.connect(window.ChangeShading)
    # ---------------------------------------------------------------------------
    # Adding the actions to the atoms
    # ---------------------------------------------------------------------------
    window.AtomsBox.toggled.connect(window.ASDVizOpt.toggle_atoms)
    window.AtomsBox.toggled.connect(window.UpdateRenderer)
    window.AtomSize.valueChanged.connect(window.ASDVizOpt.ChangeAtomsSize)
    window.AtomSize.valueChanged.connect(window.UpdateRenderer)
    window.AtomOpaq.valueChanged.connect(window.ASDVizOpt.ChangeAtomsOpaq)
    window.AtomOpaq.valueChanged.connect(window.UpdateRenderer)
    window.AtomQuali.valueChanged.connect(window.ASDVizOpt.ChangeAtomsQuali)
    window.AtomQuali.valueChanged.connect(window.UpdateRenderer)
    # ---------------------------------------------------------------------------
    # Adding the actions to the neighbours
    # ---------------------------------------------------------------------------
    window.NeighAtomCheck.toggled.connect(window.ASDVizOpt.toggle_NAtoms)
    window.NeighAtomCheck.toggled.connect(window.UpdateRenderer)
    window.NeighNeighsCheck.toggled.connect(window.ASDVizOpt.toggle_Neigh)
    window.NeighNeighsCheck.toggled.connect(window.UpdateRenderer)
    window.NeighOpacitySlider.valueChanged.connect(window.ASDVizOpt.NeighOpacityUpdate)
    window.NeighOpacitySlider.valueChanged.connect(window.UpdateRenderer)
    window.AtomOpacitySlider.valueChanged.connect(window.ASDVizOpt.AtomOpacityUpdate)
    window.AtomOpacitySlider.valueChanged.connect(window.UpdateRenderer)
    window.NeighSelectLineEdit.editingFinished.connect(window.NeighbourControl)
    window.NeighSelectSlider.valueChanged.connect(window.NeighbourControl)
    # ---------------------------------------------------------------------------
    # Adding the actions to the camera
    # ---------------------------------------------------------------------------
    window.SetCamButton.clicked.connect(window.camera_handler)
    window.CamResetButton.clicked.connect(window.camera_handler)
    window.SetXView.clicked.connect(window.camera_handler)
    window.SetYView.clicked.connect(window.camera_handler)
    window.SetZView.clicked.connect(window.camera_handler)
    window.actionTake_Snapshot.clicked.connect(window.Snapshot)
    window.ParallelProjectBox.toggled.connect(window.camera_handler)
    window.ParallelScaleSlider.valueChanged.connect(window.camera_handler)
    window.ParallelScaleLineEdit.editingFinished.connect(window.camera_handler)
    # ---------------------------------------------------------------------------
    # Adding the actions to the general actors
    # ---------------------------------------------------------------------------
    window.ScalarBarCheck.toggled.connect(window.ASDVizOpt.toggle_ScalarBar)
    window.ScalarBarCheck.toggled.connect(window.UpdateRenderer)
    window.AxesCheck.toggled.connect(window.ASDVizOpt.toggle_Axes)
    window.AxesCheck.toggled.connect(window.UpdateRenderer)
    window.ClippBox.toggled.connect(window.clipperHandler)
    window.ClippPlaneXCheck.toggled.connect(window.clipperHandler)
    window.ClippPlaneYCheck.toggled.connect(window.clipperHandler)
    window.ClippPlaneZCheck.toggled.connect(window.clipperHandler)

    # ---------------------------------------------------------------------------
    # Actions for advanced visualization options
    # ---------------------------------------------------------------------------
    window.ClippingPlaneSlider.valueChanged.connect(window.clipperHandler)
    window.GlyphQualitySlider.valueChanged.connect(window.Quality_control)
    window.FocusBox.toggled.connect(window.toggle_focus)
    window.AutoFocusCheck.toggled.connect(window.toggle_autofocus)
    window.FocusSlider.valueChanged.connect(window.FocalDisk_control)
    window.FXAACheck.toggled.connect(window.FXAA_control)
    window.FXAACheck.toggled.connect(window.UpdateRenderer)
    window.SSAOCheck.toggled.connect(window.SSAO_control)
    window.SSAOCheck.toggled.connect(window.UpdateRenderer)
    window.HDRICheck.toggled.connect(window.HDRI_control)
    window.HDRICheck.toggled.connect(window.UpdateRenderer)
    window.HDRIButtonSelect.clicked.connect(window.getHDRIFile)
    window.SkyBoxCheck.toggled.connect(window.SkyBox_control)
    window.SkyBoxCheck.toggled.connect(window.UpdateRenderer)
    # window.ShadowCheck.toggled.connect(window.Shadow_control)
    # window.ShadowCheck.toggled.connect(window.UpdateRenderer)
    # Texture controls
    window.TextureSelect.clicked.connect(window.getTextureFile)
    window.ORMTextureSelect.clicked.connect(window.getORMTextureFile)
    window.NTextureSelect.clicked.connect(window.getNTextureFile)
    window.ATextureSelect.clicked.connect(window.getATextureFile)
    window.ETextureSelect.clicked.connect(window.getETextureFile)
    #
    window.TextureCheck.toggled.connect(window.Texture_control)
    window.ORMTextureCheck.toggled.connect(window.ORMTexture_control)
    window.NTextureCheck.toggled.connect(window.NTexture_control)
    window.ATextureCheck.toggled.connect(window.ATexture_control)
    window.ETextureCheck.toggled.connect(window.ETexture_control)
    #
    window.RenDiffuseSlider.valueChanged.connect(window.RenDiffuse_control)
    window.RenAmbientSlider.valueChanged.connect(window.RenAmbient_control)
    window.RenSpecularSlider.valueChanged.connect(window.RenSpecular_control)
    window.RenSpecularPowerSlider.valueChanged.connect(window.RenSpecularPower_control)
    window.PBREmissionSlider.valueChanged.connect(window.PBREmission_control)
    window.PBROcclusionSlider.valueChanged.connect(window.PBROcclusion_control)
    window.PBRRoughnessSlider.valueChanged.connect(window.PBRRoughness_control)
    window.PBRMetallicSlider.valueChanged.connect(window.PBRMetallic_control)
    # ---------------------------------------------------------------------------
    # Adding the action to display the time step labels
    # ---------------------------------------------------------------------------
    window.TimeStepBox.toggled.connect(window.ASDVizOpt.toggle_time_label)
    window.TimeStepBox.toggled.connect(window.UpdateRenderer)
    # ---------------------------------------------------------------------------
    # Adding the actions to the energy contributions
    # ---------------------------------------------------------------------------
    window.TotEneButton.toggled.connect(window.set_energy_proj)
    window.ExcEneButton.toggled.connect(window.set_energy_proj)
    window.DMEneButton.toggled.connect(window.set_energy_proj)
    window.AniEneButton.toggled.connect(window.set_energy_proj)
    window.BqEneButton.toggled.connect(window.set_energy_proj)
    window.BqDMEneButton.toggled.connect(window.set_energy_proj)
    window.PdEneButton.toggled.connect(window.set_energy_proj)
    window.BextEneButton.toggled.connect(window.set_energy_proj)
    window.DipEneButton.toggled.connect(window.set_energy_proj)
    window.ChirEneButton.toggled.connect(window.set_energy_proj)
    window.EneDensButton.toggled.connect(window.toggle_EneActor)
    window.EneSiteGlyphs.toggled.connect(window.toggle_EneActor)
    # ---------------------------------------------------------------------------
    # Setting up validators
    # ---------------------------------------------------------------------------
    window.CamElevationLineEdit.setValidator(QDoubleValidator())
    window.CamAzimuthLineEdit.setValidator(QDoubleValidator())
    window.CamRollLineEdit.setValidator(QDoubleValidator())
    window.CamPitchLineEdit.setValidator(QDoubleValidator())
    window.CamYawLineEdit.setValidator(QDoubleValidator())
    window.FocalPosX.setValidator(QDoubleValidator())
    window.FocalPosY.setValidator(QDoubleValidator())
    window.FocalPosZ.setValidator(QDoubleValidator())
    window.CamPosX.setValidator(QDoubleValidator())
    window.CamPosY.setValidator(QDoubleValidator())
    window.CamPosZ.setValidator(QDoubleValidator())
    window.ParallelScaleLineEdit.setValidator(QDoubleValidator())
    return


##########################################################################
# @brief Function defining the toolbar and actions for the matplotlib backend.
# @details This contains several buttons handling the different types of plots
# that can be performed in the matplotlib API
# @author Jonathan Chico
##########################################################################
def Plot_Menu_and_Toolbar_Setup(window):
    """Functions defining the toolbar and actions associated to it for the Matplotlib backend.
    This contains several buttons and actions handling the different types of plots
    modes available in the matplotlib API.

    Args
    ----------
        window: QMainWindow object where the toolbar is located.

    Author
    ----------
    Jonathan Chico
    """

    window.MatPlotToolbar.setFixedHeight(24)
    window.ASD_PY_Layout.insertWidget(0, window.MatPlotToolbar)
    # ---------------------------------------------------------------------------
    # Button to exit the program
    # ---------------------------------------------------------------------------
    window.PlotExitButton = QToolButton()
    window.PlotExitButton.setText("Exit")
    window.PlotExitButton.setToolTip("Exit the UppASD visualizer")
    window.PlotExitButton.setWhatsThis("Exit the UppASD visualizer")
    window.PlotExitButton.setStatusTip("Exit the UppASD visualizer")
    window.PlotExitButton.clicked.connect(window.close)
    window.PlotSaveButton = QToolButton()
    window.PlotSaveButton.setText("Save")
    window.PlotSaveButton.setToolTip("Save the current image")
    window.PlotSaveButton.setWhatsThis("Save the current image")
    window.PlotSaveButton.setStatusTip("Save the current image")
    window.PlotSaveButton.clicked.connect(window.SaveFig)
    # ---------------------------------------------------------------------------
    # Button to plot the correlation function
    # ---------------------------------------------------------------------------
    window.actionS_q_w = QToolButton()
    window.actionS_q_w.setText("Correlation")
    window.actionS_q_w.setToolTip("Plot the correlation function")
    window.actionS_q_w.setWhatsThis("Plot the correlation function")
    window.actionS_q_w.setStatusTip("Plot the correlation function")
    window.actionS_q_w.clicked.connect(window.PlottingSelector)
    window.actionS_q_w.clicked.connect(window.PlottingWrapper)
    # ---------------------------------------------------------------------------
    # Button to plot the magnetization averages
    # ---------------------------------------------------------------------------
    window.actionAverages = QToolButton()
    window.actionAverages.setText("Averages")
    window.actionAverages.setToolTip("Plot the average magnetization")
    window.actionAverages.setWhatsThis("Plot the average magnetization")
    window.actionAverages.setStatusTip("Plot the average magnetization")
    window.actionAverages.clicked.connect(window.PlottingSelector)
    window.actionAverages.clicked.connect(window.PlottingWrapper)
    # ---------------------------------------------------------------------------
    # Button to plot the total energies
    # ---------------------------------------------------------------------------
    window.actionTotEnergy = QToolButton()
    window.actionTotEnergy.setText("Energies")
    window.actionTotEnergy.setToolTip("Plot the total energy")
    window.actionTotEnergy.setWhatsThis("Plot the total energy")
    window.actionTotEnergy.setStatusTip("Plot the total energy")
    window.actionTotEnergy.clicked.connect(window.PlottingSelector)
    window.actionTotEnergy.clicked.connect(window.PlottingWrapper)

    # ---------------------------------------------------------------------------
    # Button to plot the magnetic moment trajectories
    # ---------------------------------------------------------------------------
    window.actionTrajectory = QToolButton()
    window.actionTrajectory.setText("Trajectory")
    window.actionTrajectory.setToolTip("Plot the single moment trajectory")
    window.actionTrajectory.setWhatsThis("Plot the single moment trajectory")
    window.actionTrajectory.setStatusTip("Plot the single moment trajectory")
    window.actionTrajectory.clicked.connect(window.PlottingSelector)
    window.actionTrajectory.clicked.connect(window.PlottingWrapper)
    # ---------------------------------------------------------------------------
    # Adding the buttons to the toolbar
    # ---------------------------------------------------------------------------
    window.MatPlotToolbar.addWidget(window.PlotExitButton)
    window.MatPlotToolbar.addSeparator()
    window.MatPlotToolbar.addWidget(window.PlotSaveButton)
    window.MatPlotToolbar.addSeparator()
    window.MatPlotToolbar.addWidget(window.actionS_q_w)
    window.MatPlotToolbar.addWidget(window.actionAverages)
    window.MatPlotToolbar.addWidget(window.actionTotEnergy)
    window.MatPlotToolbar.addWidget(window.actionTrajectory)
    window.MatPlotToolbar.addSeparator()
    # --------------------------------------------------------------------------------
    # Setting misc actions
    # --------------------------------------------------------------------------------
    window.SqwCoolwarm.toggled.connect(window.Sqw_ColorMapSelector)
    window.SqwSpectral.toggled.connect(window.Sqw_ColorMapSelector)
    window.SqwBlackbody.toggled.connect(window.Sqw_ColorMapSelector)
    window.Sqw_x.toggled.connect(window.SQW_Proj_Select)
    window.Sqw_y.toggled.connect(window.SQW_Proj_Select)
    window.Sqw_z.toggled.connect(window.SQW_Proj_Select)
    window.Sqw_2.toggled.connect(window.SQW_Proj_Select)
    window.AMSDispCheckBox.toggled.connect(window.PlottingWrapper)
    window.SqwDispCheckBox.toggled.connect(window.PlottingWrapper)
    window.ABCorrWidth.sliderMoved.connect(window.SqwWidthChanger)
    window.Plot_M_x.toggled.connect(window.PlotMagDirSelector)
    window.Plot_M_y.toggled.connect(window.PlotMagDirSelector)
    window.Plot_M_z.toggled.connect(window.PlotMagDirSelector)
    window.Plot_M_tot.toggled.connect(window.PlotMagDirSelector)
    window.ABLineWidth.valueChanged.connect(window.PlotLineChanger)
    window.ABMarkerSize.valueChanged.connect(window.PlotMarkerChanger)
    window.ABXMajGrid.toggled.connect(window.PlotXGridToggle)
    window.ABYMajGrid.toggled.connect(window.PlotYGridToggle)
    window.ABAMSGrid.toggled.connect(window.PlotSQWGridToggle)
    window.ABAMSGrid.toggled.connect(window.PlotAMSGridToggle)
    # --------------------------------------------------------------------------------
    # Setting energy actions
    # --------------------------------------------------------------------------------
    window.EneTotCheck.toggled.connect(window.PlotEneCompSelector)
    window.EneExcCheck.toggled.connect(window.PlotEneCompSelector)
    window.EneAniCheck.toggled.connect(window.PlotEneCompSelector)
    window.EneDMCheck.toggled.connect(window.PlotEneCompSelector)
    window.EnePdCheck.toggled.connect(window.PlotEneCompSelector)
    window.EneBqDMCheck.toggled.connect(window.PlotEneCompSelector)
    window.EneBqCheck.toggled.connect(window.PlotEneCompSelector)
    window.EneDipCheck.toggled.connect(window.PlotEneCompSelector)
    window.EneExtCheck.toggled.connect(window.PlotEneCompSelector)
    window.EneLSFCheck.toggled.connect(window.PlotEneCompSelector)
    window.EneChirCheck.toggled.connect(window.PlotEneCompSelector)
    # --------------------------------------------------------------------------------
    # Setting file actions
    # --------------------------------------------------------------------------------
    window.actionYaml_File.triggered.connect(window.getPlotFile)
    window.actionS_q_w_File.triggered.connect(window.getPlotFile)
    window.actionAMS_File.triggered.connect(window.getPlotFile)
    window.actionAverages_File.triggered.connect(window.getPlotFile)
    window.actionTrajectory_File.triggered.connect(window.getPlotFile)
    return


##########################################################################
# @brief Toolbar and UI connections for the Input generator functions.
# @author Jonathan Chico
##########################################################################
def Input_Toolbar_Setup(window):
    """Functions defining the toolbar and actions associated to it for the Input generator backend.
    This contains several buttons and actions handling the different functions for the
    generation of the input from the GUI.

    Args
    ----------
        window: QMainWindow object where the toolbar is located.

    Author
    ----------
    Jonathan Chico
    """

    window.InputToolbar.setFixedHeight(24)
    window.ASDInp_Layout.insertWidget(0, window.InputToolbar)
    # --------------------------------------------------------------------------------
    # Button to exit the program
    # --------------------------------------------------------------------------------
    window.InputExitButton = QToolButton()
    window.InputExitButton.setText("Exit")
    window.InputExitButton.setToolTip("Exit the UppASD visualizer")
    window.InputExitButton.setWhatsThis("Exit the UppASD visualizer")
    window.InputExitButton.setStatusTip("Exit the UppASD visualizer")
    window.InputExitButton.clicked.connect(window.close)
    window.InputToolbar.addWidget(window.InputExitButton)
    window.InputToolbar.addSeparator()
    # --------------------------------------------------------------------------------
    # Set actions
    # --------------------------------------------------------------------------------
    window.InpDoneButton.clicked.connect(window.WriteInputFile)
    window.InpDMCheck.clicked.connect(window.getInpFile)
    window.InpXCCheck.clicked.connect(window.getInpFile)
    window.InpMAECheck.clicked.connect(window.getInpFile)
    window.InpPseudoCheck.clicked.connect(window.getInpFile)
    window.InpBqCheck.clicked.connect(window.getInpFile)
    window.InpBqDMCheck.clicked.connect(window.getInpFile)
    window.InpPseudDipSelect.clicked.connect(window.getInpFile)
    window.InpBiQSelect.clicked.connect(window.getInpFile)
    window.InpBiQDMSelect.clicked.connect(window.getInpFile)
    window.InpInitMag4ReadButton.clicked.connect(window.getInpFile)
    window.InpPosButtonSelect.clicked.connect(window.getInpFile)
    window.InpMomButtonSelect.clicked.connect(window.getInpFile)
    window.InpSetFinMomfileButton.clicked.connect(window.getInpFile)
    window.InpSetIniMomfileButton.clicked.connect(window.getInpFile)
    window.InpInitMag4CreateButton.clicked.connect(window.OpenWindow)
    window.InpPosButtonCreate.clicked.connect(window.OpenWindow)
    window.InpMomButtonCreate.clicked.connect(window.OpenWindow)
    window.InpInitBox.clicked.connect(window.ToggleInitPhase)
    window.InpMeasureLLG.clicked.connect(window.ToggleInitPhase)
    window.InpLLGMeasureBox.clicked.connect(window.ToggleInitPhase)
    window.InpMeasureMCMet.clicked.connect(window.ToggleInitPhase)
    window.InpMeasureMCHeat.clicked.connect(window.ToggleInitPhase)
    window.InpMeasureMCBox.clicked.connect(window.ToggleInitPhase)
    window.InpSetPhases.clicked.connect(window.OpenWindow)
    window.InpMeasureGNEB.clicked.connect(window.ToggleHessians)
    window.InpMeasureMCMet.clicked.connect(window.ToggleHessians)
    window.InpMeasureMCHeat.clicked.connect(window.ToggleHessians)
    window.InpMeasureLLG.clicked.connect(window.ToggleHessians)
    window.InpMeasureGnebBox.clicked.connect(window.ToggleHessians)
    window.InpConfBox.clicked.connect(window.ToggleHessians)
    window.InpInitmag2Box.clicked.connect(window.ToggleHessians)
    window.InpInitmag1Check.clicked.connect(window.ToggleHessians)
    window.InpInitmag3Check.clicked.connect(window.ToggleHessians)
    window.InpInitMag6Check.clicked.connect(window.ToggleHessians)
    window.InpInitMag4Check.clicked.connect(window.ToggleHessians)
    window.InpInitmag7Check.clicked.connect(window.ToggleHessians)
    window.InpCancelButton.clicked.connect(window.ResetInputs)
    window.InpRunSimButton.clicked.connect(window.RunSimulation)
    window.InpJfileButtonSelect.clicked.connect(window.getInpFile)
    window.InpJfileButtonCreate.clicked.connect(window.OpenWindow)
    window.InpSqQpoints.clicked.connect(window.getInpFile)
    window.InpMagnonQuickButton.clicked.connect(window.MagnonQuickSetup)
    window.InpDMButtonSelect.clicked.connect(window.getInpFile)
    window.InpDMButtonCreate.clicked.connect(window.OpenWindow)
    window.InpKfileButtonSelect.clicked.connect(window.getInpFile)
    window.InpKfileButtonCreate.clicked.connect(window.OpenWindow)
    window.InpImportCIFButton.clicked.connect(window.ImportSystem)
    window.InpImportSPRKKRButton.clicked.connect(window.ImportSystem)
    window.InpImportRSLMTOButton.clicked.connect(window.ImportSystem)

    # Structure Templates
    window.InpTemplateSCButton.clicked.connect(
        lambda: window.SetStructureTemplate("sc")
    )
    window.InpTemplateBCCButton.clicked.connect(
        lambda: window.SetStructureTemplate("bcc")
    )
    window.InpTemplateBCC2TypesButton.clicked.connect(
        lambda: window.SetStructureTemplate("bcc2")
    )
    window.InpTemplateFCCButton.clicked.connect(
        lambda: window.SetStructureTemplate("fcc")
    )
    window.InpTemplateHCPButton.clicked.connect(
        lambda: window.SetStructureTemplate("hcp")
    )

    return


def InteractiveDock_Setup(window):
    """
    Interface for buttons related to the interactive simulations.

    Inputs:
            window  :   QMainWindow
    """
    window.IntSStepButton.clicked.connect(window.IntButtons)
    window.IntResetButton.clicked.connect(window.IntButtons)
    window.IntMomentButton.clicked.connect(window.IntButtons)
    window.IntMCMSimButton.clicked.connect(window.IntButtons)
    window.IntMCHSimButton.clicked.connect(window.IntButtons)
    window.IntSDSlider.valueChanged.connect(window.SetSDSliderValue)
    window.IntMCSlider.valueChanged.connect(window.SetMCSliderValue)
    window.IntTempLine.editingFinished.connect(window.UpdateInteractiveVtk)
    window.IntB_xLine.editingFinished.connect(window.UpdateInteractiveVtk)
    window.IntB_yLine.editingFinished.connect(window.UpdateInteractiveVtk)
    window.IntB_zLine.editingFinished.connect(window.UpdateInteractiveVtk)
    window.IntScreenshot.clicked.connect(window.InteractiveScreenshot)
    window.IntScreenshotTic.toggled.connect(window.InteractiveScreenshotTic)


# ------------------------------------------------------------------------------------
# @brief Function to update the UI objects.
# @details Function to update the UI objects. Namely used to deal with the update between
# changes in an object which could have been obtained from other inputs, e.g. a menu
# that can change the same properties that the docket. This should keep the UI  consistent
# and ensure clarity.
# @author Jonathan Chico
# ------------------------------------------------------------------------------------
def UpdateUI(window):
    """Function to update the UI objects. Namely used to deal with the update between
    changes in an object which could have been obtained from other inputs, e.g. a menu
    that can change the same properties that the docket. This should keep the UI consistent
    and ensure clarity.

    Args
    ----------
        window: QMainWindow object where the toolbar is located.

    Author
    ----------
    Jonathan Chico
    """

    if window.sender() == window.InpInitBox:
        if window.InpInitBox.isChecked():
            window.InitPhaseOptionsBox.setEnabled(True)
        else:
            window.InitPhaseOptionsBox.setEnabled(False)
    elif window.sender() == window.InpMeasureLLG:
        if window.InpMeasureLLG.isChecked():
            window.InpLLGMeasureBox.setChecked(True)
            QSignalBlocker(window.InpLLGMeasureBox)
            window.InpMeasureMCBox.setChecked(False)
            QSignalBlocker(window.InpMeasureMCBox)
            window.InpMeasureGnebBox.setChecked(False)
            QSignalBlocker(window.InpMeasureGnebBox)
            window.HessFinCheck.setEnabled(False)
            window.HessInitCheck.setEnabled(False)
            window.HessSPCheck.setEnabled(False)
    elif window.sender() == window.InpLLGMeasureBox:
        if window.InpLLGMeasureBox.isChecked():
            window.InpMeasureLLG.setChecked(True)
            QSignalBlocker(window.InpMeasureLLG)
            window.InpMeasureMCBox.setChecked(False)
            QSignalBlocker(window.InpMeasureMCBox)
            window.InpMeasureGnebBox.setChecked(False)
            QSignalBlocker(window.InpMeasureGnebBox)
            window.HessFinCheck.setEnabled(False)
            window.HessInitCheck.setEnabled(False)
            window.HessSPCheck.setEnabled(False)
    elif (
        window.sender() == window.InpMeasureMCMet
        or window.sender() == window.InpMeasureMCHeat
    ):
        if window.InpMeasureMCMet.isChecked() or window.InpMeasureMCHeat.isChecked():
            window.InpLLGMeasureBox.setChecked(False)
            QSignalBlocker(window.InpLLGMeasureBox)
            window.InpMeasureMCBox.setChecked(True)
            QSignalBlocker(window.InpMeasureMCBox)
            window.InpMeasureGnebBox.setChecked(False)
            QSignalBlocker(window.InpMeasureGnebBox)
            window.HessFinCheck.setEnabled(False)
            window.HessInitCheck.setEnabled(False)
            window.HessSPCheck.setEnabled(False)
    elif window.sender() == window.InpMeasureMCBox:
        if window.InpMeasureMCBox.isChecked():
            window.InpLLGMeasureBox.setChecked(False)
            QSignalBlocker(window.InpLLGMeasureBox)
            window.InpMeasureMCMet.setChecked(True)
            QSignalBlocker(window.InpMeasureMCMet)
            window.InpMeasureGnebBox.setChecked(False)
            QSignalBlocker(window.InpMeasureGnebBox)
            window.HessFinCheck.setEnabled(False)
            window.HessInitCheck.setEnabled(False)
            window.HessSPCheck.setEnabled(False)
    elif window.sender() == window.InpMeasureGNEB:
        if window.InpMeasureGNEB.isChecked():
            window.InpMeasureGnebBox.setChecked(True)
            QSignalBlocker(window.InpMeasureGnebBox)
            window.InpMeasureMCMet.setChecked(False)
            QSignalBlocker(window.InpMeasureMCMet)
            window.InpLLGMeasureBox.setChecked(False)
            QSignalBlocker(window.InpLLGMeasureBox)
            window.HessFinCheck.setEnabled(True)
            window.HessInitCheck.setEnabled(True)
            window.HessSPCheck.setEnabled(True)
            window.InpGNEBMeasureBox.setEnabled(True)
            window.InpMeasureGNEBCI.setEnabled(True)
    elif window.sender() == window.InpMeasureGnebBox:
        if window.InpMeasureGnebBox.isChecked():
            window.InpMeasureGNEB.setChecked(True)
            QSignalBlocker(window.InpMeasureGnebBox)
            window.InpMeasureMCMet.setChecked(False)
            QSignalBlocker(window.InpMeasureMCMet)
            window.InpLLGMeasureBox.setChecked(False)
            QSignalBlocker(window.InpLLGMeasureBox)
            window.HessFinCheck.setEnabled(True)
            window.HessInitCheck.setEnabled(True)
            window.HessSPCheck.setEnabled(True)
            window.InpGNEBMeasureBox.setEnabled(True)
            window.InpMeasureGNEBCI.setEnabled(True)
    elif window.sender() == window.InpInitmag1Check:
        if window.InpInitmag1Check.isChecked():
            window.InpInitmag3Check.setChecked(False)
            window.InpInitmag2Box.setChecked(False)
            window.InpInitMag6Check.setChecked(False)
            window.InpConfBox.setChecked(False)
            window.InpInitMag4Check.setChecked(False)
            window.InpInitmag7Check.setChecked(False)
    elif window.sender() == window.InpInitmag3Check:
        if window.InpInitmag3Check.isChecked():
            window.InpInitmag1Check.setChecked(False)
            window.InpInitmag2Box.setChecked(False)
            window.InpInitMag6Check.setChecked(False)
            window.InpConfBox.setChecked(False)
            window.InpInitMag4Check.setChecked(False)
            window.InpInitmag7Check.setChecked(False)
    elif window.sender() == window.InpInitmag2Box:
        if window.InpInitmag2Box.isChecked():
            window.InpInitmag3Check.setChecked(False)
            window.InpInitmag1Check.setChecked(False)
            window.InpInitMag6Check.setChecked(False)
            window.InpConfBox.setChecked(False)
            window.InpInitMag4Check.setChecked(False)
            window.InpInitmag7Check.setChecked(False)
    elif window.sender() == window.InpInitMag6Check:
        if window.InpInitMag6Check.isChecked():
            window.InpInitmag2Box.setChecked(False)
            window.InpInitmag3Check.setChecked(False)
            window.InpInitmag1Check.setChecked(False)
            window.InpConfBox.setChecked(False)
            window.InpInitMag4Check.setChecked(False)
            window.InpInitmag7Check.setChecked(False)
    elif window.sender() == window.InpConfBox:
        if window.InpConfBox.isChecked():
            window.InpInitmag2Box.setChecked(False)
            window.InpInitmag3Check.setChecked(False)
            window.InpInitmag1Check.setChecked(False)
            window.InpInitMag6Check.setChecked(False)
            window.InpInitMag4Check.setChecked(False)
            window.InpInitmag7Check.setChecked(False)
    elif window.sender() == window.InpInitMag4Check:
        if window.InpInitMag4Check.isChecked():
            window.InpInitmag7Check.setChecked(False)
    elif window.sender() == window.InpInitmag7Check:
        if window.InpInitmag7Check.isChecked():
            window.InpInitMag4Check.setChecked(False)
    return
