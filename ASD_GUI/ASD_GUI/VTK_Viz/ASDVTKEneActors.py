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
        self.timer_count = 0
        self.camera_pos = np.zeros(3, dtype=np.float32)
        self.camera_focal = np.zeros(3, dtype=np.float32)
        self.camera_yaw = 0.0
        self.camera_roll = 0.0
        self.camera_pitch = 0.0
        self.camera_azimuth = 0.0
        self.camera_elevation = 0.0
        self.cluster_disp = ASDdata.cluster_flag
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
        self.src = vtk.vtkPolyData()
        self.src.SetPoints(ASDdata.coord)
        self.src.GetPointData().SetScalars(ASDdata.energies[0])
        self.scalar_range = self.src.GetScalarRange()
        # -----------------------------------------------------------------------
        # Finding useful geometrical information of the sample
        # -----------------------------------------------------------------------
        # Finding the middle of the sample
        # Also making sure that if the sample is 2D one has no problem with boudings
        # this is mostly useful if splatters are used
        (
            self.xmin,
            self.xmax,
            self.ymin,
            self.ymax,
            self.zmin,
            self.zmax,
        ) = self.src.GetBounds()
        if self.xmin == self.xmax:
            self.xmin = 0.0
            self.xmax = 1.0
        if self.ymin == self.ymax:
            self.ymin = 0.0
            self.ymax = 1.0
        if self.zmin == self.zmax:
            self.zmin = 0.0
            self.zmax = 1.0
        self.xmid = (self.xmin + self.xmax) * 0.5
        self.ymid = (self.ymin + self.ymax) * 0.5
        self.zmid = (self.zmin + self.zmax) * 0.5
        self.height = max(self.xmax, self.ymax, self.zmax) * 1.75
        # Auxiliary data to find max distance in the x-direction
        self.dist_x = np.absolute(self.xmax - self.xmin)
        # Auxiliary data to find max distance in the y-direction
        self.dist_y = np.absolute(self.ymax - self.ymin)
        # Auxiliary data to find max distance in the z-direction
        self.dist_z = np.absolute(self.zmax - self.zmin)
        self.camera_pos[0] = self.xmid
        self.camera_pos[1] = self.ymid
        self.camera_pos[2] = self.height
        self.camera_focal[0] = self.xmid
        self.camera_focal[1] = self.ymid
        self.camera_focal[2] = self.zmid
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
            self.EneDensMethod = vtk.vtkDelaunay2D()
            self.EneDensMethod.SetInputData(self.src)
            self.EneDensMethod.BoundingTriangulationOff()
            self.EneDensMethod.SetTolerance(0.005)
            # Time the execution of the delaunay tessellation
            SM_timer = vtk.vtkExecutionTimer()
            SM_timer.SetFilter(self.EneDensMethod)
            self.EneDensMethod.Update()
            SM = SM_timer.GetElapsedWallClockTime()
            print("2D Delaunay:", SM)
            # Creating the mapper for the smooth surfaces
            self.EneDensMap = vtk.vtkDataSetMapper()
            self.EneDensMap.SetScalarRange(self.scalar_range)
            self.EneDensMap.SetInputConnection(self.EneDensMethod.GetOutputPort())
            self.EneDensMap.SetLookupTable(self.lut)
            self.EneDensMap.SetColorModeToMapScalars()
            self.EneDensMap.Update()
            # Creating the actor for the smooth surfaces
            self.EneDensActor = vtk.vtkLODActor()
            self.EneDensActor.SetMapper(self.EneDensMap)
            self.EneDensActor.GetProperty().SetOpacity(0.75)
            self.EneDensActor.GetProperty().EdgeVisibilityOff()
        else:
            # -------------------------------------------------------------------
            # Setting the parameters for the visualization of 3D structures with
            # splatters
            # -------------------------------------------------------------------
            self.EneDensMethod = vtk.vtkShepardMethod()
            self.EneDensMethod.SetInputData(self.src)
            self.EneDensMethod.SetModelBounds(
                self.xmin,
                self.xmax,
                self.ymin,
                self.ymax,
                self.zmin,
                self.zmax,
            )
            # This should get rid of the problems when trying to map very thin
            # structures in 2D
            if self.dist_x == min(self.dist_x, self.dist_y, self.dist_z):
                self.EneDensMethod.SetSampleDimensions(
                    3, int(self.ymax), int(self.zmax)
                )
            elif self.dist_y == min(self.dist_x, self.dist_y, self.dist_z):
                self.EneDensMethod.SetSampleDimensions(
                    int(self.xmax), 3, int(self.zmax)
                )
            elif self.dist_z == min(self.dist_x, self.dist_y, self.dist_z):
                self.EneDensMethod.SetSampleDimensions(
                    int(self.xmax), int(self.ymax), 3
                )
            # This parameter determines how far in the sample (normalized to 1) the
            # method will look to interpolate, greatly affects performance
            self.EneDensMethod.SetMaximumDistance(0.1)
            # Time the execution of the checkerboard splatter
            SP_timer = vtk.vtkExecutionTimer()
            SP_timer.SetFilter(self.EneDensMethod)
            self.EneDensMethod.Update()
            SP = SP_timer.GetElapsedWallClockTime()
            print("3D Shepard Method:", SP)
            # Mapper for the image obtained from the 3D reconstruction method
            self.EneDensMap = vtk.vtkSmartVolumeMapper()
            self.EneDensMap.SetBlendModeToComposite()
            self.EneDensMap.SetInputConnection(self.EneDensMethod.GetOutputPort())
            # Function for the opacity gradient
            volumeGradientOpacity = vtk.vtkPiecewiseFunction()
            volumeGradientOpacity.AddPoint(-1, 0.25)
            volumeGradientOpacity.AddPoint(0.5, 0.75)
            volumeGradientOpacity.AddPoint(1.0, 1.0)
            # Properties of the volume to be rendered
            self.volumeProperty = vtk.vtkVolumeProperty()
            self.volumeProperty.SetInterpolationType(1)
            self.volumeProperty.SetColor(self.transfer_func)
            self.volumeProperty.SetAmbient(0.6)
            self.volumeProperty.SetDiffuse(0.6)
            self.volumeProperty.SetSpecular(0.1)
            self.volumeProperty.SetGradientOpacity(volumeGradientOpacity)
            # Volume actor, this works in a different way than LOD actors
            self.EneDensActor = vtk.vtkVolume()
            self.EneDensActor.SetMapper(self.EneDensMap)
            self.EneDensActor.SetProperty(self.volumeProperty)
        # -----------------------------------------------------------------------
        # Energy spheres
        # -----------------------------------------------------------------------
        self.EneAtom = vtk.vtkSphereSource()
        self.EneAtom.SetRadius(0.50)
        self.EneAtom.SetThetaResolution(10)
        self.EneAtom.SetPhiResolution(10)
        # -----------------------------------------------------------------------
        # Set the mapper for the energies
        # -----------------------------------------------------------------------
        self.EneMapper = vtk.vtkGlyph3DMapper()
        self.EneMapper.SetSourceConnection(self.EneAtom.GetOutputPort())
        self.EneMapper.SetInputData(self.src)
        self.EneMapper.SetScalarRange(self.scalar_range)
        self.EneMapper.SetScaleFactor(1.00)
        self.EneMapper.SetScaleModeToNoDataScaling()
        self.EneMapper.SetLookupTable(self.lut)
        self.EneMapper.SetColorModeToMapScalars()
        self.EneMapper.Update()
        # -----------------------------------------------------------------------
        # Energy actors
        # -----------------------------------------------------------------------
        self.EneActor = vtk.vtkLODActor()
        self.EneActor.SetMapper(self.EneMapper)
        self.EneActor.GetProperty().SetSpecular(0.3)
        self.EneActor.GetProperty().SetSpecularPower(60)
        self.EneActor.GetProperty().SetAmbient(0.2)
        self.EneActor.GetProperty().SetDiffuse(0.8)
        self.EneActor.VisibilityOff()
        # -----------------------------------------------------------------------
        # Setting information of the renderer
        # -----------------------------------------------------------------------
        # Define the renderer
        # Add the actors to the scene
        if ASDdata.flag2D:
            ren.AddActor(self.EneDensActor)
        else:
            ren.AddViewProp(self.EneDensActor)
        ren.AddActor(self.EneActor)
        # Defining the camera directions
        ren.GetActiveCamera().Azimuth(self.camera_azimuth)
        ren.GetActiveCamera().Elevation(self.camera_elevation)
        ren.GetActiveCamera().Yaw(self.camera_yaw)
        ren.GetActiveCamera().Roll(self.camera_roll)
        ren.GetActiveCamera().Pitch(self.camera_pitch)
        ren.GetActiveCamera().SetFocalPoint(self.camera_focal)
        ren.GetActiveCamera().SetPosition(self.camera_pos)
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
            self.src.GetPointData().SetScalars(ASDdata.energies[0])
            if window.EneDensButton.isChecked():
                self.EneDensMap.SetScalarRange(self.src.GetScalarRange())
            if window.EneSiteGlyphs.isChecked():
                self.EneMapper.SetScalarRange(self.src.GetScalarRange())
        # -----------------------------------------------------------------------
        # If the exchange energy is clicked update it
        # -----------------------------------------------------------------------
        if window.ExcEneButton.isChecked():
            self.src.GetPointData().SetScalars(ASDdata.energies[1])
            if window.EneDensButton.isChecked():
                self.EneDensMap.SetScalarRange(self.src.GetScalarRange())
            if window.EneSiteGlyphs.isChecked():
                self.EneMapper.SetScalarRange(self.src.GetScalarRange())
        # -----------------------------------------------------------------------
        # If the DMI energy is clicked update it
        # -----------------------------------------------------------------------
        if window.DMEneButton.isChecked():
            self.src.GetPointData().SetScalars(ASDdata.energies[2])
            if window.EneDensButton.isChecked():
                self.EneDensMap.SetScalarRange(self.src.GetScalarRange())
            if window.EneSiteGlyphs.isChecked():
                self.EneMapper.SetScalarRange(self.src.GetScalarRange())
        # -----------------------------------------------------------------------
        # If the anisotropy energy is clicked update it
        # -----------------------------------------------------------------------
        if window.AniEneButton.isChecked():
            self.src.GetPointData().SetScalars(ASDdata.energies[3])
            if window.EneDensButton.isChecked():
                self.EneDensMap.SetScalarRange(self.src.GetScalarRange())
            if window.EneSiteGlyphs.isChecked():
                self.EneMapper.SetScalarRange(self.src.GetScalarRange())
        # -----------------------------------------------------------------------
        # If the biquadratic energy is clicked update it
        # -----------------------------------------------------------------------
        if window.BqEneButton.isChecked():
            self.src.GetPointData().SetScalars(ASDdata.energies[4])
            if window.EneDensButton.isChecked():
                self.EneDensMap.SetScalarRange(self.src.GetScalarRange())
            if window.EneSiteGlyphs.isChecked():
                self.EneMapper.SetScalarRange(self.src.GetScalarRange())
        # -----------------------------------------------------------------------
        # If the biquadratic DMI energy is clicked update it
        # -----------------------------------------------------------------------
        if window.BqDMEneButton.isChecked():
            self.src.GetPointData().SetScalars(ASDdata.energies[5])
            if window.EneDensButton.isChecked():
                self.EneDensMap.SetScalarRange(self.src.GetScalarRange())
            if window.EneSiteGlyphs.isChecked():
                self.EneMapper.SetScalarRange(self.src.GetScalarRange())
        # -----------------------------------------------------------------------
        # If the pseudodipolar energy is clicked update it
        # -----------------------------------------------------------------------
        if window.PdEneButton.isChecked():
            self.src.GetPointData().SetScalars(ASDdata.energies[6])
            if window.EneDensButton.isChecked():
                self.EneDensMap.SetScalarRange(self.src.GetScalarRange())
            if window.EneSiteGlyphs.isChecked():
                self.EneMapper.SetScalarRange(self.src.GetScalarRange())
        # -----------------------------------------------------------------------
        # If the Zeeman energy is clicked update it
        # -----------------------------------------------------------------------
        if window.BextEneButton.isChecked():
            self.src.GetPointData().SetScalars(ASDdata.energies[7])
            if window.EneDensButton.isChecked():
                self.EneDensMap.SetScalarRange(self.src.GetScalarRange())
            if window.EneSiteGlyphs.isChecked():
                self.EneMapper.SetScalarRange(self.src.GetScalarRange())
        # -----------------------------------------------------------------------
        # If the dipolar energy is clicked update it
        # -----------------------------------------------------------------------
        if window.DipEneButton.isChecked():
            self.src.GetPointData().SetScalars(ASDdata.energies[8])
            if window.EneDensButton.isChecked():
                self.EneDensMap.SetScalarRange(self.src.GetScalarRange())
            if window.EneSiteGlyphs.isChecked():
                self.EneMapper.SetScalarRange(self.src.GetScalarRange())
        # -----------------------------------------------------------------------
        # If the CHIR  energy is clicked update it
        # -----------------------------------------------------------------------
        if window.ChirEneButton.isChecked():
            self.src.GetPointData().SetScalars(ASDdata.energies[9])
            if window.EneDensButton.isChecked():
                self.EneDensMap.SetScalarRange(self.src.GetScalarRange())
            if window.EneSiteGlyphs.isChecked():
                self.EneMapper.SetScalarRange(self.src.GetScalarRange())
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
