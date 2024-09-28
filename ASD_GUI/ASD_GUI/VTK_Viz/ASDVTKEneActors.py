"""@package ASDVTKEneActors
Contains the class that contains the objects needed for the viualization of the site
dependent energy obtained from an UppASD calculation.
The rendering of this energy is time dependent that is one can obtain movies showing
the time evolution of the local contributions to the energy

Author
----------
Jonathan Chico
"""
# pylint: disable=invalid-name, no-name-in-module, no-member
import numpy as np
import vtk


##########################################################################
# @brief Wrapper class containing the functions for the energy visualization.
# @details Wrapper class containing the functions for the energy visualization.
# This class defines the actors for both volume rendering of the energy, as well
# as glyph enegry rendering.
# It also deals with the update of the energy when different contributions wish
# to be plotted.
# @author Jonathan Chico
##########################################################################


class ASDEneActors:
    ##########################################################################
    # @brief Function where the actual VTK actors are defined for the Energy visualization.
    # @details Function where the actual VTK actors are defined for the Energy visualization.
    # It makes sure that the camera is located in such a way as to conform to the
    # sample dimensions.
    #
    # It also tries to find the best way to visualize 2D and 3D like structures by
    # choosing on the fly between different tesellation methods.
    # @author Jonathan Chico
    ##########################################################################
    def __init__(self):
        self.active = False
    
    def Add_EneActors(self, ren, renWin, iren, ASDdata):
        """
        Add energy actors to the VTK renderer.

        This method initializes various VTK components and adds energy actors
        to the provided renderer.  It supports both 2D and 3D visualizations
        based on the ASDdata flag.

        Parameters:
        ren (vtkRenderer): The VTK renderer to which actors will be added.
        renWin (vtkRenderWindow): The VTK render window.
        iren (vtkRenderWindowInteractor): The VTK render window interactor.
        ASDdata (ASDData): The data structure containing the necessary information
        for visualization.

        Returns:
        None

        Author:
        Jonathan Chico
        """
        # -----------------------------------------------------------------------
        # Initialize variables
        # -----------------------------------------------------------------------
        ASDEneActors.timer_count = 0
        ASDEneActors.camera_pos = np.zeros(3, dtype=np.float32)
        ASDEneActors.camera_focal = np.zeros(3, dtype=np.float32)
        ASDEneActors.camera_yaw = 0.0
        ASDEneActors.camera_roll = 0.0
        ASDEneActors.camera_pitch = 0.0
        ASDEneActors.camera_azimuth = 0.0
        ASDEneActors.camera_elevation = 0.0
        ASDEneActors.cluster_disp = ASDdata.cluster_flag
        # -----------------------------------------------------------------------
        # Look up tables for colors
        # -----------------------------------------------------------------------
        # This is a diverging RWB color mapping based on the work of Kenneth
        # Moreland and with the vtk examples provided by Andrew Maclean
        if ASDdata.flag2D:
            self.lut = (
                vtk.vtkLookupTable()
            )  # Lookup table to assign color to the actors
            num_colors = 256
            self.lut.SetNumberOfTableValues(num_colors)
            # Function to map data to the colormap
            self.transfer_func = vtk.vtkColorTransferFunction()
            self.transfer_func.SetColorSpaceToDiverging()
            self.transfer_func.AddRGBPoint(0, 0.230, 0.299, 0.754)
            self.transfer_func.AddRGBPoint(1, 0.706, 0.016, 0.150)
            for ii, ss in enumerate(
                [float(xx) / float(num_colors) for xx in range(num_colors)]
            ):
                cc = self.transfer_func.GetColor(ss)
                self.lut.SetTableValue(ii, cc[0], cc[1], cc[2], 1.0)
            self.lut.Build()
        else:
            self.lut = (
                vtk.vtkLookupTable()
            )  # Lookup table to assign color to the actors
            num_colors = 256
            self.lut.SetNumberOfTableValues(num_colors)
            # Function to map data to the colormap
            self.transfer_func = vtk.vtkColorTransferFunction()
            self.transfer_func.SetColorSpaceToDiverging()
            self.transfer_func.AddRGBPoint(-0, 0.230, 0.299, 0.754)
            self.transfer_func.AddRGBPoint(1, 0.706, 0.016, 0.150)
            for ii, ss in enumerate(
                [float(xx) / float(num_colors) for xx in range(num_colors)]
            ):
                cc = self.transfer_func.GetColor(ss)
                self.lut.SetTableValue(ii, cc[0], cc[1], cc[2], 1.0)
            self.lut.Build()
        # -----------------------------------------------------------------------
        # Data structures for the generation of the smooth grid
        # -----------------------------------------------------------------------
        # Passing the data from the full system to the PolyData
        ASDEneActors.src = vtk.vtkPolyData()
        ASDEneActors.src.SetPoints(ASDdata.coord)
        ASDEneActors.src.GetPointData().SetScalars(ASDdata.energies[0])
        ASDEneActors.scalar_range = ASDEneActors.src.GetScalarRange()
        # -----------------------------------------------------------------------
        # Finding useful geometrical information of the sample
        # -----------------------------------------------------------------------
        # Finding the middle of the sample
        # Also making sure that if the sample is 2D one has no problem with boudings
        # this is mostly useful if splatters are used
        (
            ASDEneActors.xmin,
            ASDEneActors.xmax,
            ASDEneActors.ymin,
            ASDEneActors.ymax,
            ASDEneActors.zmin,
            ASDEneActors.zmax,
        ) = ASDEneActors.src.GetBounds()
        if ASDEneActors.xmin == ASDEneActors.xmax:
            ASDEneActors.xmin = 0.0
            ASDEneActors.xmax = 1.0
        if ASDEneActors.ymin == ASDEneActors.ymax:
            ASDEneActors.ymin = 0.0
            ASDEneActors.ymax = 1.0
        if ASDEneActors.zmin == ASDEneActors.zmax:
            ASDEneActors.zmin = 0.0
            ASDEneActors.zmax = 1.0
        ASDEneActors.xmid = (ASDEneActors.xmin + ASDEneActors.xmax) * 0.5
        ASDEneActors.ymid = (ASDEneActors.ymin + ASDEneActors.ymax) * 0.5
        ASDEneActors.zmid = (ASDEneActors.zmin + ASDEneActors.zmax) * 0.5
        ASDEneActors.height = (
            max(ASDEneActors.xmax, ASDEneActors.ymax, ASDEneActors.zmax) * 1.75
        )
        # Auxiliary data to find max distance in the x-direction
        self.dist_x = np.absolute(ASDEneActors.xmax - ASDEneActors.xmin)
        # Auxiliary data to find max distance in the y-direction
        self.dist_y = np.absolute(ASDEneActors.ymax - ASDEneActors.ymin)
        # Auxiliary data to find max distance in the z-direction
        self.dist_z = np.absolute(ASDEneActors.zmax - ASDEneActors.zmin)
        ASDEneActors.camera_pos[0] = ASDEneActors.xmid
        ASDEneActors.camera_pos[1] = ASDEneActors.ymid
        ASDEneActors.camera_pos[2] = ASDEneActors.height
        ASDEneActors.camera_focal[0] = ASDEneActors.xmid
        ASDEneActors.camera_focal[1] = ASDEneActors.ymid
        ASDEneActors.camera_focal[2] = ASDEneActors.zmid
        # The delaunay tesellation seems to be the best way to transform the point cloud
        # to a surface for volume rendering, the problem is that it is too slow for large
        # data sets, meaning that the best option is first to prune out the data to ensure
        # that one has a manageable number of data points over which to do the construction
        # surface reconstruction and splatter techniques also can be used to generate something
        # akin to the kind of surfaces we want. The issue is that they transform the data to a
        # regular mesh by default. And thus it is a problem for most kind of
        # systems
        if ASDdata.flag2D:
            # Passing the data to generate a triangulation of the data
            ASDEneActors.EneDensMethod = vtk.vtkDelaunay2D()
            ASDEneActors.EneDensMethod.SetInputData(ASDEneActors.src)
            ASDEneActors.EneDensMethod.BoundingTriangulationOff()
            ASDEneActors.EneDensMethod.SetTolerance(0.005)
            # Time the execution of the delaunay tessellation
            SM_timer = vtk.vtkExecutionTimer()
            SM_timer.SetFilter(ASDEneActors.EneDensMethod)
            ASDEneActors.EneDensMethod.Update()
            SM = SM_timer.GetElapsedWallClockTime()
            print("2D Delaunay:", SM)
            # Creating the mapper for the smooth surfaces
            ASDEneActors.EneDensMap = vtk.vtkDataSetMapper()
            ASDEneActors.EneDensMap.SetScalarRange(ASDEneActors.scalar_range)
            ASDEneActors.EneDensMap.SetInputConnection(
                ASDEneActors.EneDensMethod.GetOutputPort()
            )
            ASDEneActors.EneDensMap.SetLookupTable(self.lut)
            ASDEneActors.EneDensMap.SetColorModeToMapScalars()
            ASDEneActors.EneDensMap.Update()
            # Creating the actor for the smooth surfaces
            ASDEneActors.EneDensActor = vtk.vtkLODActor()
            ASDEneActors.EneDensActor.SetMapper(ASDEneActors.EneDensMap)
            ASDEneActors.EneDensActor.GetProperty().SetOpacity(0.75)
            ASDEneActors.EneDensActor.GetProperty().EdgeVisibilityOff()
        else:
            # -------------------------------------------------------------------
            # Setting the parameters for the visualization of 3D structures with
            # splatters
            # -------------------------------------------------------------------
            ASDEneActors.EneDensMethod = vtk.vtkShepardMethod()
            ASDEneActors.EneDensMethod.SetInputData(ASDEneActors.src)
            ASDEneActors.EneDensMethod.SetModelBounds(
                ASDEneActors.xmin,
                ASDEneActors.xmax,
                ASDEneActors.ymin,
                ASDEneActors.ymax,
                ASDEneActors.zmin,
                ASDEneActors.zmax,
            )
            # This should get rid of the problems when trying to map very thin
            # structures in 2D
            if self.dist_x == min(self.dist_x, self.dist_y, self.dist_z):
                ASDEneActors.EneDensMethod.SetSampleDimensions(
                    3, int(ASDEneActors.ymax), int(ASDEneActors.zmax)
                )
            elif self.dist_y == min(self.dist_x, self.dist_y, self.dist_z):
                ASDEneActors.EneDensMethod.SetSampleDimensions(
                    int(ASDEneActors.xmax), 3, int(ASDEneActors.zmax)
                )
            elif self.dist_z == min(self.dist_x, self.dist_y, self.dist_z):
                ASDEneActors.EneDensMethod.SetSampleDimensions(
                    int(ASDEneActors.xmax), int(ASDEneActors.ymax), 3
                )
            # This parameter determines how far in the sample (normalized to 1) the
            # method will look to interpolate, greatly affects performance
            ASDEneActors.EneDensMethod.SetMaximumDistance(0.1)
            # Time the execution of the checkerboard splatter
            SP_timer = vtk.vtkExecutionTimer()
            SP_timer.SetFilter(ASDEneActors.EneDensMethod)
            ASDEneActors.EneDensMethod.Update()
            SP = SP_timer.GetElapsedWallClockTime()
            print("3D Shepard Method:", SP)
            # Mapper for the image obtained from the 3D reconstruction method
            ASDEneActors.EneDensMap = vtk.vtkSmartVolumeMapper()
            ASDEneActors.EneDensMap.SetBlendModeToComposite()
            ASDEneActors.EneDensMap.SetInputConnection(
                ASDEneActors.EneDensMethod.GetOutputPort()
            )
            # Function for the opacity gradient
            volumeGradientOpacity = vtk.vtkPiecewiseFunction()
            volumeGradientOpacity.AddPoint(-1, 0.25)
            volumeGradientOpacity.AddPoint(0.5, 0.75)
            volumeGradientOpacity.AddPoint(1.0, 1.0)
            # Properties of the volume to be rendered
            ASDEneActors.volumeProperty = vtk.vtkVolumeProperty()
            ASDEneActors.volumeProperty.SetInterpolationType(1)
            ASDEneActors.volumeProperty.SetColor(self.transfer_func)
            ASDEneActors.volumeProperty.SetAmbient(0.6)
            ASDEneActors.volumeProperty.SetDiffuse(0.6)
            ASDEneActors.volumeProperty.SetSpecular(0.1)
            ASDEneActors.volumeProperty.SetGradientOpacity(
                volumeGradientOpacity)
            # Volume actor, this works in a different way than LOD actors
            ASDEneActors.EneDensActor = vtk.vtkVolume()
            ASDEneActors.EneDensActor.SetMapper(ASDEneActors.EneDensMap)
            ASDEneActors.EneDensActor.SetProperty(self.volumeProperty)
        # -----------------------------------------------------------------------
        # Energy spheres
        # -----------------------------------------------------------------------
        ASDEneActors.EneAtom = vtk.vtkSphereSource()
        ASDEneActors.EneAtom.SetRadius(0.50)
        ASDEneActors.EneAtom.SetThetaResolution(10)
        ASDEneActors.EneAtom.SetPhiResolution(10)
        # -----------------------------------------------------------------------
        # Set the mapper for the energies
        # -----------------------------------------------------------------------
        ASDEneActors.EneMapper = vtk.vtkGlyph3DMapper()
        ASDEneActors.EneMapper.SetSourceConnection(
            ASDEneActors.EneAtom.GetOutputPort())
        ASDEneActors.EneMapper.SetInputData(ASDEneActors.src)
        ASDEneActors.EneMapper.SetScalarRange(ASDEneActors.scalar_range)
        ASDEneActors.EneMapper.SetScaleFactor(1.00)
        ASDEneActors.EneMapper.SetScaleModeToNoDataScaling()
        ASDEneActors.EneMapper.SetLookupTable(self.lut)
        ASDEneActors.EneMapper.SetColorModeToMapScalars()
        ASDEneActors.EneMapper.Update()
        # -----------------------------------------------------------------------
        # Energy actors
        # -----------------------------------------------------------------------
        ASDEneActors.EneActor = vtk.vtkLODActor()
        ASDEneActors.EneActor.SetMapper(ASDEneActors.EneMapper)
        ASDEneActors.EneActor.GetProperty().SetSpecular(0.3)
        ASDEneActors.EneActor.GetProperty().SetSpecularPower(60)
        ASDEneActors.EneActor.GetProperty().SetAmbient(0.2)
        ASDEneActors.EneActor.GetProperty().SetDiffuse(0.8)
        ASDEneActors.EneActor.VisibilityOff()
        # -----------------------------------------------------------------------
        # Setting information of the renderer
        # -----------------------------------------------------------------------
        # Define the renderer
        # Add the actors to the scene
        if ASDdata.flag2D:
            ren.AddActor(ASDEneActors.EneDensActor)
        else:
            ren.AddViewProp(ASDEneActors.EneDensActor)
        ren.AddActor(ASDEneActors.EneActor)
        # Defining the camera directions
        ren.GetActiveCamera().Azimuth(ASDEneActors.camera_azimuth)
        ren.GetActiveCamera().Elevation(ASDEneActors.camera_elevation)
        ren.GetActiveCamera().Yaw(ASDEneActors.camera_yaw)
        ren.GetActiveCamera().Roll(ASDEneActors.camera_roll)
        ren.GetActiveCamera().Pitch(ASDEneActors.camera_pitch)
        ren.GetActiveCamera().SetFocalPoint(ASDEneActors.camera_focal)
        ren.GetActiveCamera().SetPosition(ASDEneActors.camera_pos)
        ren.GetActiveCamera().SetViewUp(0, 1, 0)
        # -----------------------------------------------------------------------
        # Start the renderer
        # -----------------------------------------------------------------------
        iren.Start()
        renWin.Render()
        return

    ##########################################################################
    # @brief Update the energy for visualization
    # @details Update the energy for visualization. This function takes care of
    # the correct energy component being plotted when certain signals are emitted
    # from the GUI.
    # @author Jonathan Chico
    ##########################################################################
    def UpdateEnergy(self, window, ASDdata, ASDGenActors, renWin):
        """
        Update the energy visualization based on the selected energy type.

        This method reads the energy data, updates the visualization actors
        based on the selected energy type, and refreshes the UI components
        such as progress bar and time label.
        It also triggers a render of the visualization window.

        Parameters:
        window (object): The main window object containing UI elements.
        ASDdata (object): The data object containing energy information and methods to
        read energy data.
        ASDGenActors (object): The general actors object for updating the time label.
        renWin (object): The render window object to trigger rendering.

        Author:
        Jonathan Chico
        """
        # -----------------------------------------------------------------------
        # Actually reading the data
        # -----------------------------------------------------------------------
        (ASDdata.energies, ASDdata.number_time_steps, ASDdata.time_sep) = (
            ASDdata.readEnergyData(
                ASDdata.eneFile,
                window.current_time,
                ASDdata.nrAtoms,
                ASDdata.number_time_steps,
            )
        )
        # -----------------------------------------------------------------------
        # If the total energy is clicked update it
        # -----------------------------------------------------------------------
        if window.TotEneButton.isChecked():
            ASDEneActors.src.GetPointData().SetScalars(ASDdata.energies[0])
            if window.EneDensButton.isChecked():
                ASDEneActors.EneDensMap.SetScalarRange(
                    ASDEneActors.src.GetScalarRange()
                )
            if window.EneSiteGlyphs.isChecked():
                ASDEneActors.EneMapper.SetScalarRange(
                    ASDEneActors.src.GetScalarRange())
        # -----------------------------------------------------------------------
        # If the exchange energy is clicked update it
        # -----------------------------------------------------------------------
        if window.ExcEneButton.isChecked():
            ASDEneActors.src.GetPointData().SetScalars(ASDdata.energies[1])
            if window.EneDensButton.isChecked():
                ASDEneActors.EneDensMap.SetScalarRange(
                    ASDEneActors.src.GetScalarRange()
                )
            if window.EneSiteGlyphs.isChecked():
                ASDEneActors.EneMapper.SetScalarRange(
                    ASDEneActors.src.GetScalarRange())
        # -----------------------------------------------------------------------
        # If the DMI energy is clicked update it
        # -----------------------------------------------------------------------
        if window.DMEneButton.isChecked():
            ASDEneActors.src.GetPointData().SetScalars(ASDdata.energies[2])
            if window.EneDensButton.isChecked():
                ASDEneActors.EneDensMap.SetScalarRange(
                    ASDEneActors.src.GetScalarRange()
                )
            if window.EneSiteGlyphs.isChecked():
                ASDEneActors.EneMapper.SetScalarRange(
                    ASDEneActors.src.GetScalarRange()
                )
        # -----------------------------------------------------------------------
        # If the anisotropy energy is clicked update it
        # -----------------------------------------------------------------------
        if window.AniEneButton.isChecked():
            ASDEneActors.src.GetPointData().SetScalars(ASDdata.energies[3])
            if window.EneDensButton.isChecked():
                ASDEneActors.EneDensMap.SetScalarRange(
                    ASDEneActors.src.GetScalarRange()
                )
            if window.EneSiteGlyphs.isChecked():
                ASDEneActors.EneMapper.SetScalarRange(
                    ASDEneActors.src.GetScalarRange())
        # -----------------------------------------------------------------------
        # If the biquadratic energy is clicked update it
        # -----------------------------------------------------------------------
        if window.BqEneButton.isChecked():
            ASDEneActors.src.GetPointData().SetScalars(ASDdata.energies[4])
            if window.EneDensButton.isChecked():
                ASDEneActors.EneDensMap.SetScalarRange(
                    ASDEneActors.src.GetScalarRange()
                )
            if window.EneSiteGlyphs.isChecked():
                ASDEneActors.EneMapper.SetScalarRange(
                    ASDEneActors.src.GetScalarRange())
        # -----------------------------------------------------------------------
        # If the biquadratic DMI energy is clicked update it
        # -----------------------------------------------------------------------
        if window.BqDMEneButton.isChecked():
            ASDEneActors.src.GetPointData().SetScalars(ASDdata.energies[5])
            if window.EneDensButton.isChecked():
                ASDEneActors.EneDensMap.SetScalarRange(
                    ASDEneActors.src.GetScalarRange()
                )
            if window.EneSiteGlyphs.isChecked():
                ASDEneActors.EneMapper.SetScalarRange(
                    ASDEneActors.src.GetScalarRange())
        # -----------------------------------------------------------------------
        # If the pseudodipolar energy is clicked update it
        # -----------------------------------------------------------------------
        if window.PdEneButton.isChecked():
            ASDEneActors.src.GetPointData().SetScalars(ASDdata.energies[6])
            if window.EneDensButton.isChecked():
                ASDEneActors.EneDensMap.SetScalarRange(
                    ASDEneActors.src.GetScalarRange()
                )
            if window.EneSiteGlyphs.isChecked():
                ASDEneActors.EneMapper.SetScalarRange(
                    ASDEneActors.src.GetScalarRange())
        # -----------------------------------------------------------------------
        # If the Zeeman energy is clicked update it
        # -----------------------------------------------------------------------
        if window.BextEneButton.isChecked():
            ASDEneActors.src.GetPointData().SetScalars(ASDdata.energies[7])
            if window.EneDensButton.isChecked():
                ASDEneActors.EneDensMap.SetScalarRange(
                    ASDEneActors.src.GetScalarRange()
                )
            if window.EneSiteGlyphs.isChecked():
                ASDEneActors.EneMapper.SetScalarRange(
                    ASDEneActors.src.GetScalarRange())
        # -----------------------------------------------------------------------
        # If the dipolar energy is clicked update it
        # -----------------------------------------------------------------------
        if window.DipEneButton.isChecked():
            ASDEneActors.src.GetPointData().SetScalars(ASDdata.energies[8])
            if window.EneDensButton.isChecked():
                ASDEneActors.EneDensMap.SetScalarRange(
                    ASDEneActors.src.GetScalarRange()
                )
            if window.EneSiteGlyphs.isChecked():
                ASDEneActors.EneMapper.SetScalarRange(
                    ASDEneActors.src.GetScalarRange())
        # -----------------------------------------------------------------------
        # If the CHIR  energy is clicked update it
        # -----------------------------------------------------------------------
        if window.ChirEneButton.isChecked():
            ASDEneActors.src.GetPointData().SetScalars(ASDdata.energies[9])
            if window.EneDensButton.isChecked():
                ASDEneActors.EneDensMap.SetScalarRange(
                    ASDEneActors.src.GetScalarRange()
                )
            if window.EneSiteGlyphs.isChecked():
                ASDEneActors.EneMapper.SetScalarRange(
                    ASDEneActors.src.GetScalarRange())
        # -----------------------------------------------------------------------
        # Update the UI
        # -----------------------------------------------------------------------
        window.ProgressBar.setValue(
            window.current_time * 100 / (ASDdata.number_time_steps - 1)
        )
        window.ProgressLabel.setText(f"   {int(window.ProgressBar.value())}%")
        time_label = f"{float(window.TimeStepLineEdit.text()) * ASDdata.time_sep[window.current_time] * 1e9: 4.2f} ns"
        ASDGenActors.time_label.SetInput(time_label)
        # -----------------------------------------------------------------------
        # Take a snapshot
        # -----------------------------------------------------------------------
        renWin.Render()
        return
