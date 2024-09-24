"""@package ASDVTKMomActors
Wrapper class to add the VTK  actors for the visualization of UppASD data in
moment visualization mode.
It contains the needed data to add the actors, modify them, as well as some helper
functions to change them.

Author
----------
Jonathan Chico
"""
# pylint: disable=invalid-name, no-name-in-module, no-member

import vtk
import numpy as np
from vtk.util import numpy_support


class ASDMomActors:
    ##########################################################################
    # @brief Main wrapper to add the needed actors for visualization of the moments
    # @details Main wrapper to add the needed actors for visualization of the moments.
    # Class that contains the data structures for creation of the glyphs
    # for the visualization of the magnetic moments.
    # It also has the capacity to create tessellations for the visualization of
    # volume vendering.
    # @author Jonathan Chico
    ##########################################################################
    def Add_MomActors(self, ren, renWin, iren, ASDdata, window):
        """Main wrapper to add the needed actors for visualization of the moments.
        Class that contains the data structures for creation of the glyphs or the
        visualization of the magnetic moments. It also has the capacity to create
        tessellations for the visualization of volume vendering.

        Args:
            ren: current renderer.
            renWin: current rendering window.
            iren: current interactor for the renderer.
            ASDdata: class where the data read from the ASD simulations is stored.
            window: QMainWindow object where the visualization is performed.

        Author
        ----------
        Jonathan Chico
        """

        self.timer_count = 0
        self.camera_pos = np.zeros(3, dtype=np.float32)
        self.camera_focal = np.zeros(3, dtype=np.float32)
        self.camera_yaw = 0.0
        self.camera_roll = 0.0
        self.camera_pitch = 0.0
        self.camera_azimuth = 0.0
        self.camera_elevation = 0.0
        self.kmc_disp = ASDdata.kmc_flag
        self.cluster_disp = ASDdata.cluster_flag
        self.glob_color = ASDdata.colors
        # -----------------------------------------------------------------------
        # Look up tables for colors
        # -----------------------------------------------------------------------
        # This is a diverging RWB color mapping based on the work of Kenneth
        # Moreland and with the vtk examples provided by Andrew Maclean
        if ASDdata.flag2D:
            self.lut = vtk.vtkLookupTable()
            num_colors = 256
            self.lut.SetNumberOfTableValues(num_colors)
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
            self.lut = vtk.vtkLookupTable()
            num_colors = 256
            self.lut.SetNumberOfTableValues(num_colors)
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
        self.src.GetPointData().SetScalars(ASDdata.colors[2])
        self.src.GetPointData().SetVectors(ASDdata.moments)
        scalar_range = self.src.GetScalarRange()
        # -----------------------------------------------------------------------
        # Finding useful geometrical information of the sample
        # -----------------------------------------------------------------------
        # Finding the middle of the sample
        # Also making sure that if the sample is 2D one has no problem with bounds
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
        self.height = (
            max(self.xmax, self.ymax, self.zmax) * 1.75
        )
        self.dist_x = np.absolute(self.xmax - self.xmin)
        self.dist_y = np.absolute(self.ymax - self.ymin)
        self.dist_z = np.absolute(self.zmax - self.zmin)
        self.camera_pos[0] = self.xmid
        self.camera_pos[1] = self.ymid
        self.camera_pos[2] = self.height
        self.camera_focal[0] = self.xmid
        self.camera_focal[1] = self.ymid
        self.camera_focal[2] = self.zmid
        # -----------------------------------------------------------------------
        # Data structures for the spins
        # -----------------------------------------------------------------------
        # Passing the data from the full system to the PolyData
        self.src_spins = vtk.vtkPolyData()
        self.src_spins.SetPoints(ASDdata.coord)
        self.src_spins.GetPointData().SetScalars(ASDdata.colors[2])
        self.src_spins.GetPointData().SetVectors(ASDdata.moments)
        scalar_range_spins = self.src_spins.GetScalarRange()
        # -----------------------------------------------------------------------
        self.ungrid = vtk.vtkUnstructuredGrid()
        self.ungrid.SetPoints(ASDdata.coord)
        self.ungrid.GetPointData().SetScalars(
            self.src.GetPointData().GetScalars()
        )
        # writer = vtk.vtkUnstructuredGridWriter()
        # writer.SetFileName('out.vtk')
        # writer.SetInputData(ASDMomActors.ungrid)
        # writer.Write()
        # ASDMomActors.ungrid.SetPoints(ASDdata.coord)
        # ASDMomActors.ungrid.GetPointData().SetScalars(ASDdata.colors[2])
        # ASDMomActors.ungrid.GetPointData().SetVectors(ASDdata.moments)

        # -----------------------------------------------------------------------
        # The delaunay tessellation seems to be the best way to transform the point cloud
        # to a surface for volume rendering, the problem is that it is too slow for large
        # data sets, meaning that the best option is first to prune out the data to ensure
        # that one has a manageable number of data points over which to do the construction
        # surface reconstruction and splatter techniques also can be used to generate something
        # akin to the kind of surfaces we want. The issue is that they transform the data to a
        # regular mesh by default. And thus it is a problem for most kind of
        # systems
        if ASDdata.flag2D:
            # Passing the data to generate a triangulation of the data
            self.MagDensMethod = vtk.vtkDelaunay2D()
            self.MagDensMethod.SetInputData(self.src)
            self.MagDensMethod.BoundingTriangulationOff()
            self.MagDensMethod.SetTolerance(0.005)
            # Time the execution of the delaunay tessellation
            SM_timer = vtk.vtkExecutionTimer()
            SM_timer.SetFilter(self.MagDensMethod)
            self.MagDensMethod.Update()
            SM = SM_timer.GetElapsedWallClockTime()
            print("2D Delaunay:", SM)
            # Creating the mapper for the smooth surfaces
            self.MagDensMap = vtk.vtkDataSetMapper()
            self.MagDensMap.SetScalarRange(scalar_range)
            self.MagDensMap.SetInputConnection(
                self.MagDensMethod.GetOutputPort()
            )
            self.MagDensMap.SetLookupTable(self.lut)
            self.MagDensMap.SetColorModeToMapScalars()
            self.MagDensMap.Update()
            # Creating the actor for the smooth surfaces
            self.MagDensActor = vtk.vtkLODActor()
            self.MagDensActor.SetMapper(self.MagDensMap)
            self.MagDensActor.GetProperty().SetOpacity(0.75)
            self.MagDensActor.GetProperty().EdgeVisibilityOff()
            if window.DensBox.isChecked():
                self.MagDensActor.VisibilityOn()
            else:
                self.MagDensActor.VisibilityOff()
        else:
            # -------------------------------------------------------------------
            # Setting the parameters for the visualization of 3D structures with
            # splatters
            # -------------------------------------------------------------------
            self.MagDensMethod = vtk.vtkGaussianSplatter()
            self.MagDensMethod.SetInputData(self.src)
            # print( ASDMomActors.MagDensMethod.GetSampleDimensions())
            # ASDMomActors.MagDensMethod.SetSampleDimensions([10,10,2])
            # print( ASDMomActors.MagDensMethod.GetSampleDimensions())
            # -------------------------------------------------------------------
            # Options for the Gaussian splatter. These determine the quality of the
            # rendering, increasing the radius smoothens out the volume but performance
            # decreases rapidly
            # -------------------------------------------------------------------
            dist = np.asarray(self.src.GetPoint(0)) - np.asarray(
                self.src.GetPoint(1)
            )
            norm = np.sqrt(dist.dot(dist))
            if norm < 10:
                rad_fac = 0.040
            else:
                rad_fac = 0.40
            self.MagDensMethod.SetRadius(rad_fac)
            self.MagDensMethod.ScalarWarpingOn()
            # -------------------------------------------------------------------
            # The exponent factor determines how fast the gaussian splatter decay
            # they again can be used to improve quality at the sake of rendering time
            # -------------------------------------------------------------------
            self.MagDensMethod.SetExponentFactor(-10)
            self.MagDensMethod.NormalWarpingOn()
            self.MagDensMethod.SetEccentricity(10)
            # -------------------------------------------------------------------
            # The Null value can be used to try to eliminate contributions not belonging
            # to the actual sample
            # -------------------------------------------------------------------
            self.MagDensMethod.SetNullValue(-10)
            # -------------------------------------------------------------------
            # Set the actual size of the rendering model
            # -------------------------------------------------------------------
            self.MagDensMethod.SetModelBounds(
                self.xmin,
                self.xmax,
                self.ymin,
                self.ymax,
                self.zmin,
                self.zmax,
            )
            # This should get rid of the problems when trying to map very thin
            # structures in 2D
            if (
                self.dist_x == min(self.dist_x, self.dist_y, self.dist_z)
                and self.dist_x < 3
            ):
                self.MagDensMethod.SetSampleDimensions(
                    3, int(self.ymax), int(self.zmax)
                )
            elif (
                self.dist_y == min(self.dist_x, self.dist_y, self.dist_z)
                and self.dist_y < 3
            ):
                self.MagDensMethod.SetSampleDimensions(
                    int(self.xmax), 3, int(self.zmax)
                )
            elif (
                self.dist_z == min(self.dist_x, self.dist_y, self.dist_z)
                and self.dist_z < 3
            ):
                self.MagDensMethod.SetSampleDimensions(
                    int(self.xmax), int(self.ymax), 3
                )
            # Timming for the execution of the volume creation
            SP_timer = vtk.vtkExecutionTimer()
            SP_timer.SetFilter(self.MagDensMethod)
            self.MagDensMethod.Update()
            SP = SP_timer.GetElapsedWallClockTime()
            print("3D vtkGaussianSplatter Method:", SP)
            # Scalar opacities
            funcOpacityScalar = vtk.vtkPiecewiseFunction()
            funcOpacityScalar.AddPoint(-1.00, 0.00)
            funcOpacityScalar.AddPoint(0.00, 0.05)
            funcOpacityScalar.AddPoint(0.50, 0.50)
            funcOpacityScalar.AddPoint(0.75, 1.00)
            # Gradient opacities
            volumeGradientOpacity = vtk.vtkPiecewiseFunction()
            volumeGradientOpacity.AddPoint(0.000, 0.0)
            volumeGradientOpacity.AddPoint(0.001, 1.0)
            volumeGradientOpacity.AddPoint(1.000, 1.0)
            # Volume properties
            self.volumeProperty = vtk.vtkVolumeProperty()
            self.volumeProperty.SetColor(self.transfer_func)
            self.volumeProperty.SetInterpolationTypeToLinear()
            self.volumeProperty.SetAmbient(0.6)
            self.volumeProperty.SetDiffuse(0.6)
            self.volumeProperty.SetSpecular(0.1)
            self.volumeProperty.SetGradientOpacity(volumeGradientOpacity)
            self.volumeProperty.SetScalarOpacity(funcOpacityScalar)
            # Volume Mapper
            # ASDMomActors.MagDensMethod = vtk.vtkDataSetTriangleFilter()
            # ASDMomActors.MagDensMethod.SetInputData(ASDMomActors.ungrid)
            # ASDMomActors.MagDensMethod.Update()

            self.MagDensMap = vtk.vtkSmartVolumeMapper()
            # ASDMomActors.MagDensMap = vtk.vtkUnstructuredGridVolumeRayCastMapper()
            self.MagDensMap.SetInputConnection(
                self.MagDensMethod.GetOutputPort()
            )
            # ASDMomActors.MagDensMap.Update()

            # Volume Actor
            self.MagDensActor = vtk.vtkVolume()
            self.MagDensActor.SetMapper(self.MagDensMap)
            self.MagDensActor.SetProperty(self.volumeProperty)
            if window.DensBox.isChecked():
                self.MagDensActor.VisibilityOn()
            else:
                self.MagDensActor.VisibilityOff()
            self.MagDensActor.Update()
        # -----------------------------------------------------------------------
        # Data structures for the contours
        # -----------------------------------------------------------------------
        # Define the contour filters
        contours = vtk.vtkContourFilter()
        contours.SetInputConnection(self.MagDensMethod.GetOutputPort())
        # This generates the contours, it will do 5 between the -1 and 0.5
        # range
        cont_num = 5
        range_cont = (-1, 0.5)
        contours.GenerateValues(cont_num, range_cont)
        # Map the contours to graphical primitives
        contMapper = vtk.vtkPolyDataMapper()
        contMapper.SetInputConnection(contours.GetOutputPort())
        contMapper.SetScalarVisibility(False)  # colored contours
        contMapper.SetScalarRange(scalar_range)
        # Create an actor for the contours
        self.contActor = vtk.vtkLODActor()
        self.contActor.SetMapper(contMapper)
        self.contActor.GetProperty().SetColor(0, 0, 0)
        self.contActor.GetProperty().SetLineWidth(1.0)
        self.contActor.VisibilityOff()

        # -----------------------------------------------------------------------
        # Setting information of the spins
        # -----------------------------------------------------------------------
        # Create vectors
        self.spinarrow = vtk.vtkArrowSource()
        self.spinarrow.SetTipRadius(0.20)
        self.spinarrow.SetShaftRadius(0.10)
        self.spinarrow.SetTipResolution(16)
        self.spinarrow.SetShaftResolution(16)

        self.spinarrowtriangles = vtk.vtkTriangleFilter()
        self.spinarrowtriangles.SetInputConnection(
            self.spinarrow.GetOutputPort()
        )

        # Calculate normals for shading
        self.spinarrownormals = vtk.vtkPolyDataNormals()
        self.spinarrownormals.SetInputConnection(
            self.spinarrowtriangles.GetOutputPort()
        )

        self.spinarrowtcoords = vtk.vtkTextureMapToCylinder()
        self.spinarrowtcoords.SetInputConnection(
            self.spinarrownormals.GetOutputPort()
        )
        self.spinarrowtcoords.PreventSeamOn()

        self.spinarrowtangents = vtk.vtkPolyDataTangents()
        self.spinarrowtangents.SetInputConnection(
            self.spinarrowtcoords.GetOutputPort()
        )

        # Create the mapper for the spins
        self.SpinMapper = vtk.vtkGlyph3DMapper()
        self.SpinMapper.SetSourceConnection(
            self.spinarrowtangents.GetOutputPort()
        )

        self.SpinMapper.SetInputData(self.src_spins)
        self.SpinMapper.SetScalarRange(scalar_range_spins)
        self.SpinMapper.SetScaleFactor(0.50)
        self.SpinMapper.SetScaleModeToNoDataScaling()
        self.SpinMapper.SetLookupTable(self.lut)
        self.SpinMapper.SetColorModeToMapScalars()
        self.SpinMapper.Update()

        # Define the vector actor for the spins
        self.Spins = vtk.vtkActor()
        # ASDMomActors.Spins = vtk.vtkLODActor()
        self.Spins.SetMapper(self.SpinMapper)
        # ASDMomActors.Spins.GetProperty().SetInterpolationToPBR()
        self.Spins.GetProperty().SetInterpolationToGouraud()
        self.Spins.GetProperty().SetSpecular(0.4)
        self.Spins.GetProperty().SetSpecularPower(80)
        self.Spins.GetProperty().SetAmbient(0.6)
        self.Spins.GetProperty().SetDiffuse(0.4)
        self.Spins.GetProperty().SetEdgeTint(0.0, 0.0, 0.0)
        if window.SpinsBox.isChecked():
            self.Spins.VisibilityOn()
        else:
            self.Spins.VisibilityOff()

        # -----------------------------------------------------------------------
        # Setting information for the atoms
        # -----------------------------------------------------------------------
        # Create vectors
        self.AtomSphere = vtk.vtkTexturedSphereSource()
        self.AtomSphere.SetRadius(0.10)
        self.AtomSphere.SetThetaResolution(12)
        self.AtomSphere.SetPhiResolution(12)

        # Create the mapper for the atoms
        self.AtomMapper = vtk.vtkGlyph3DMapper()
        self.AtomMapper.SetSourceConnection(
            self.AtomSphere.GetOutputPort()
        )
        self.AtomMapper.SetInputData(self.src_spins)
        self.AtomMapper.SetScaleFactor(1.00)
        self.AtomMapper.SetScaleModeToNoDataScaling()
        self.AtomMapper.ScalarVisibilityOff()
        # ASDMomActors.AtomMapper.SetLookupTable(self.lut)
        # ASDMomActors.AtomMapper.SetColorModeToMapScalars()
        self.AtomMapper.Update()
        # Define the sphere actor for the atoms
        colors = vtk.vtkNamedColors()
        self.Atoms = vtk.vtkLODActor()
        self.Atoms.SetMapper(self.AtomMapper)
        self.Atoms.GetProperty().SetInterpolationToGouraud()
        # ASDMomActors.Atoms.GetProperty().SetInterpolationToPBR()
        self.Atoms.GetProperty().SetAmbient(0.8)
        self.Atoms.GetProperty().SetDiffuse(0.8)
        self.Atoms.GetProperty().SetColor(colors.GetColor3d("Silver"))
        if window.AtomsBox.isChecked():
            self.Atoms.VisibilityOn()
        else:
            self.Atoms.VisibilityOff()

        # -----------------------------------------------------------------------
        # Setting information for the skybox actor
        # -----------------------------------------------------------------------
        self.SkyBox = vtk.vtkSkybox()
        self.SkyBox.VisibilityOff()
        ren.AddActor(self.SkyBox)

        if ASDdata.kmc_flag:
            # -------------------------------------------------------------------
            # Setting data structures for the KMC particle visualization
            # -------------------------------------------------------------------
            self.KMC_src = vtk.vtkPolyData()
            self.KMC_src.SetPoints(ASDdata.coord_KMC)
            # Atom sphere
            KMC_part = vtk.vtkSphereSource()
            KMC_part.SetRadius(1.75)
            KMC_part.SetThetaResolution(12)
            KMC_part.SetPhiResolution(12)
            # Atom glyph
            KMC_part_mapper = vtk.vtkGlyph3DMapper()
            KMC_part_mapper.SetInputData(self.KMC_src)
            KMC_part_mapper.SetSourceConnection(KMC_part.GetOutputPort())
            KMC_part_mapper.SetScaleFactor(0.5)
            KMC_part_mapper.ClampingOn()
            KMC_part_mapper.SetScaleModeToNoDataScaling()
            KMC_part_mapper.SetColorModeToMapScalars()
            KMC_part_mapper.Update()
            # Atoms actors
            self.KMC_part_actor = vtk.vtkLODActor()
            self.KMC_part_actor.SetMapper(KMC_part_mapper)
            self.KMC_part_actor.GetProperty().SetOpacity(0.9)
            self.KMC_part_actor.GetProperty().SetColor(0.0, 0.0, 1.0)
            self.KMC_part_actor.GetProperty().EdgeVisibilityOn()
            self.KMC_part_actor.GetProperty().SetEdgeColor(0, 0, 0)
        # -----------------------------------------------------------------------
        # Setting information of the renderer
        # -----------------------------------------------------------------------
        # Define the renderer
        # Add the actors to the scene
        if ASDdata.flag2D:
            ren.AddActor(self.MagDensActor)
        else:
            ren.AddViewProp(self.MagDensActor)
        ren.AddActor(self.Spins)
        # ren.AddActor(ASDMomActors.vector)
        ren.AddActor(self.contActor)
        ren.AddActor(self.Atoms)
        # If the KMC particles are present add them to the renderer
        if ASDdata.kmc_flag:
            ren.AddActor(self.KMC_part_actor)
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
        print(" Done")
        return

    ##########################################################################
    # @brief Update the magnetic moments for the visualization
    # @author Jonathan Chico
    ##########################################################################
    def UpdateMoments(
        self, window, ASDdata, ASDGenActors, renWin, is_interactive=False
    ):
        """Function to update the visualization of the moments as one advances in time
        or in ensembles.
        Args:
            window: QMainWindow object where the visualizations are contained.
            ASData: class containing the data read from the UppASD simulations.
            ASDGenActors: class of general actors for the VTK visualization.
            renWin: VTK window where the visualization is performed.

        Author
        ----------
        Jonathan Chico
        """
        # -----------------------------------------------------------------------
        # Read the actual data of the magnetic moments
        # -----------------------------------------------------------------------
        print("Reading data for time step:", window.current_time)
        print(self.__class__.__name__)
        (
            ASDdata.moments,
            ASDdata.colors,
            ASDdata.number_time_steps,
            ASDdata.time_sep,
        ) = ASDdata.readVectorsData(
            ASDdata.MagFile,
            window.current_time,
            ASDdata.nrAtoms,
            ASDdata.number_time_steps,
        )
        # -----------------------------------------------------------------------
        # Update the colors
        # -----------------------------------------------------------------------
        if window.DensX.isChecked():
            self.src.GetPointData().SetScalars(ASDdata.colors[0])
        if window.SpinX.isChecked():
            self.src_spins.GetPointData().SetScalars(ASDdata.colors[0])
        if window.DensY.isChecked():
            self.src.GetPointData().SetScalars(ASDdata.colors[1])
        if window.SpinY.isChecked():
            self.src_spins.GetPointData().SetScalars(ASDdata.colors[1])
        if window.DensZ.isChecked():
            self.src.GetPointData().SetScalars(ASDdata.colors[2])
        if window.SpinZ.isChecked():
            self.src_spins.GetPointData().SetScalars(ASDdata.colors[2])
        # -----------------------------------------------------------------------
        # Update the vectors
        # -----------------------------------------------------------------------
        self.src.GetPointData().SetVectors(ASDdata.moments)
        self.src_spins.GetPointData().SetVectors(ASDdata.moments)
        # -----------------------------------------------------------------------
        # If coordinate file is animated, update coordinates
        # -----------------------------------------------------------------------
        if ASDdata.nrAtoms < ASDdata.full_coord.shape[0]:
            t_off = window.current_time * ASDdata.nrAtoms
            ASDdata.coord.SetData(
                numpy_support.numpy_to_vtk(
                    ASDdata.full_coord[t_off: t_off + ASDdata.nrAtoms]
                )
            )
            self.src.SetPoints(ASDdata.coord)
        # -----------------------------------------------------------------------
        # Update the general actors
        # -----------------------------------------------------------------------
        window.ProgressBar.setValue(
            int((window.current_time - 1) * 100 / (ASDdata.number_time_steps - 1))
        )
        window.ProgressLabel.setText(f"   {int(window.ProgressBar.value())}%")
        time_label = f"{float(window.TimeStepLineEdit.text()) * ASDdata.time_sep[window.current_time-1] * 1e9: 4.2f} ns"
        ASDGenActors.time_label.SetInput(time_label)
        # -----------------------------------------------------------------------
        # Render the window
        # -----------------------------------------------------------------------
        renWin.Render()
        return

    ##########################################################################
    # @brief Change the type of glyphs used to visualize the atomistic spins
    # @details Change the type of glyphs used to visualize the atomistic spins.
    # It allows the user to change the glyphs for the atomistic spins on the fly
    # the user can choose between the following actors:
    #
    #   - Cone
    #   - Arrows
    #   - Spheres
    #   - Cubes
    # @author Jonathan Chico
    ##########################################################################
    def ChangeSpinGlyph(self, keyword):
        """Change the type of glyphs used to visualize the atomistic spins.
        It allows the user to change the glyphs for the atomistic spins on the fly
        the user can choose between the following actors:
            * Cone
            * Arrows
            * Spheres
            * Cubes

        Args:
            renWin: current VTK rendering window.
            keyword: (str) identifier keyword to choose between the different types of glyphs.

        Author
        ----------
        Jonathan Chico
        """

        if keyword == "Cubes":
            try:
                del self.spinarrow
            except AttributeError:
                pass
            try:
                del self.spinsphere
            except AttributeError:
                pass
            try:
                del self.spincones
            except AttributeError:
                pass
            self.spincube = vtk.vtkCubeSource()
            self.spincube.SetXLength(1.0)
            self.spincube.SetYLength(1.0)
            self.spincube.SetZLength(1.0)

            # Calculate TCoords for texturing
            # self.spincubetmap = vtk.vtkTextureMapToSphere()
            # self.spincubetmap.SetInputConnection(self.spincube.GetOutputPort())
            # self.spincubetmap.PreventSeamOn()

            self.SpinMapper.SetSourceConnection(
                self.spincube.GetOutputPort()
            )
            # self.SpinMapper.SetSourceConnection(self.spincubetmap.GetOutputPort())
            self.SpinMapper.ClampingOn()
            self.SpinMapper.OrientOff()
        if keyword == "Bars":
            try:
                del self.spinarrow
            except AttributeError:
                pass
            try:
                del self.spinsphere
            except AttributeError:
                pass
            try:
                del self.spincones
            except AttributeError:
                pass
            self.spincube = vtk.vtkCubeSource()
            self.spincube.SetXLength(2.0)
            self.spincube.SetYLength(0.4)
            self.spincube.SetZLength(0.4)

            # Calculate TCoords for texturing
            self.spincubetmap = vtk.vtkTextureMapToCylinder()
            self.spincubetmap.SetInputConnection(
                self.spincube.GetOutputPort()
            )
            # self.spincubetmap.AutomaticCylinderGenerationOff()
            # self.spincubetmap.SetPoint1([ 1.0,0.0,0.0])
            # self.spincubetmap.SetPoint2([-1.0,0.0,0.0])

            # Calculate TCoords for texturing
            # self.spincubetmap = vtk.vtkTextureMapToSphere()
            # self.spincubetmap.SetInputConnection(self.spincube.GetOutputPort())
            # self.spincubetmap.AutomaticSphereGenerationOff()
            # self.spincubetmap.SetCenter([ 0.0,0.0,0.0])

            self.spincubetmap.PreventSeamOff()

            self.SpinMapper.SetSourceConnection(
                self.spincubetmap.GetOutputPort()
            )
            self.SpinMapper.ClampingOn()
            self.SpinMapper.OrientOn()

        if keyword == "Spheres":
            try:
                del self.spinarrow
            except AttributeError:
                pass
            try:
                del self.spincube
            except AttributeError:
                pass
            try:
                del self.spincones
            except AttributeError:
                pass
            self.spinsphere = vtk.vtkTexturedSphereSource()
            self.spinsphere.SetRadius(0.50)
            self.spinsphere.SetThetaResolution(12)
            self.spinsphere.SetPhiResolution(12)

            # Placeholder comment for testing tangent extraction for normal textures
            # tritri = vtk.vtkTriangleFilter()
            # tritri.SetInputConnection(self.spinsphere.GetOutputPort())
            # tritan = vtk.vtkPolyDataTangents()
            # tritan.SetInputConnection(tritri.GetOutputPort())
            # self.SpinMapper.SetSourceConnection(tritan.GetOutputPort())

            self.SpinMapper.SetSourceConnection(
                self.spinsphere.GetOutputPort()
            )
            self.SpinMapper.ClampingOn()
            self.SpinMapper.OrientOn()
            # self.SpinMapper.OrientOff()
        if keyword == "Arrows":
            try:
                del self.spinsphere
            except AttributeError:
                pass
            try:
                del self.spincube
            except AttributeError:
                pass
            try:
                del self.spincones
            except AttributeError:
                pass

            # Create vectors
            self.spinarrow = vtk.vtkArrowSource()
            self.spinarrow.SetTipRadius(0.20)
            self.spinarrow.SetShaftRadius(0.10)
            self.spinarrow.SetTipResolution(12)
            self.spinarrow.SetShaftResolution(12)

            # Calculate normals for shading
            self.spinarrownormals = vtk.vtkPolyDataNormals()
            self.spinarrownormals.SetInputConnection(
                self.spinarrow.GetOutputPort()
            )

            # Calculate TCoords for texturing
            self.spinarrownormalstmap = vtk.vtkTextureMapToCylinder()
            self.spinarrownormalstmap.SetInputConnection(
                self.spinarrownormals.GetOutputPort()
            )
            self.spinarrownormalstmap.PreventSeamOn()

            self.SpinMapper.SetSourceConnection(
                self.spinarrownormalstmap.GetOutputPort()
            )
            self.SpinMapper.OrientOn()
            self.SpinMapper.Update()

        if keyword == "CenterOn":
            self.spinarrow.SetArrowOriginToCenter()

        if keyword == "CenterOff":
            self.spinarrow.SetArrowOriginToDefault()

        if keyword == "Cones":
            try:
                del self.spinsphere
            except AttributeError:
                pass
            try:
                del self.spincube
            except AttributeError:
                pass
            try:
                del self.spinarrow
            except AttributeError:
                pass

            self.spincones = vtk.vtkConeSource()
            self.spincones.SetRadius(0.50)
            self.spincones.SetHeight(1.00)
            self.spincones.SetResolution(12)

            # Calculate normals for shading
            self.spinconenormals = vtk.vtkPolyDataNormals()
            self.spinconenormals.SetInputConnection(
                self.spincones.GetOutputPort()
            )

            # Calculate TCoords for texturing
            self.spinconeormalstmap = vtk.vtkTextureMapToCylinder()
            self.spinconeormalstmap.SetInputConnection(
                self.spinconenormals.GetOutputPort()
            )
            self.spinconeormalstmap.PreventSeamOn()

            self.SpinMapper.SetSourceConnection(
                self.spinconeormalstmap.GetOutputPort()
            )
            self.SpinMapper.OrientOn()
            self.SpinMapper.Update()
        return
    
    ##########################################################################
    # Toggle options for the contours
    ##########################################################################
    def toggle_contours(self, check):
        """
        Toggles the visibility of contour actors based on the check parameter.
        """
        if check:
            self.contActor.VisibilityOn()
        else:
            self.contActor.VisibilityOff()
        return

    ##########################################################################
    # Toggle the directions arrows
    ##########################################################################
    def toggle_directions(self, check):
        """
        Toggles the visibility of MomActors.vector based on the check parameter.
        """
        if check:
            self.vector.VisibilityOn()
        else:
            self.vector.VisibilityOff()
        return

    ##########################################################################
    # Toggle the directions arrows
    ##########################################################################
    def toggle_spins(self, check):
        """
        Toggles the visibility of the spin actors based on the check value.
        """
        if check:
            self.Spins.VisibilityOn()
        else:
            self.Spins.VisibilityOff()
        return

    ##########################################################################
    # Toggle the atomic spheres
    ##########################################################################
    def toggle_atoms(self, check):
        """
        Toggles the visibility of atoms based on the given check value.
        """
        if check:
            self.Atoms.VisibilityOn()
        else:
            self.Atoms.VisibilityOff()
        return

    ##########################################################################
    # Toggle the magnetization density
    ##########################################################################
    def toggle_density(self, check):
        """
        Toggles the visibility of the magnetic density actor based on the check value.
        """
        if check:
            self.MagDensActor.VisibilityOn()
        else:
            self.MagDensActor.VisibilityOff()
        return

    ##########################################################################
    # Set the color of the magnetization density along the x axis projection
    ##########################################################################
    def set_projection(self, atype, axis):
        """
        Set the projection type and axis for visualization.
        """
        if atype == "density":
            self.src.GetPointData().SetScalars(
                self.glob_color[axis]
            )
        elif atype == "spins":
            self.src_spins.GetPointData().SetScalars(
                self.glob_color[axis]
            )
        return

    ##########################################################################
    # Set the size of the spins via the slider
    ##########################################################################
    def ChangeSpinsSize(self, value):
        """
        Adjusts the scale factor of the SpinMapper based on the provided value.
        """
        self.SpinMapper.SetScaleFactor(0.50 * value / 10)
        return

    ##########################################################################
    # Set the interpolation of the spin glyphs
    ##########################################################################
    def ChangeSpinShade(self, keyword):
        """
        Change the shading model of the Spin and Atom actors based on the given keyword.
        """
        if keyword == "Flat":
            if hasattr(self, "Spins"):
                self.Spins.GetProperty().SetInterpolationToFlat()
            if hasattr(self, "Atoms"):
                self.Atoms.GetProperty().SetInterpolationToFlat()
        elif keyword == "Gouraud":
            if hasattr(self, "Spins"):
                self.Spins.GetProperty().SetInterpolationToGouraud()
            if hasattr(self, "Atoms"):
                self.Atoms.GetProperty().SetInterpolationToGouraud()
        elif keyword == "PBR":
            if hasattr(self, "Spins"):
                self.Spins.GetProperty().SetInterpolationToPBR()
                self.Spins.GetProperty().SetMetallic(0.5)
            if hasattr(self, "Atoms"):
                self.Atoms.GetProperty().SetInterpolationToPBR()
                self.Atoms.GetProperty().SetMetallic(0.5)
        elif keyword == "Phong":
            if hasattr(self, "Spins"):
                self.Spins.GetProperty().SetInterpolationToPhong()
            if hasattr(self, "Atoms"):
                self.Atoms.GetProperty().SetInterpolationToPhong()

        return

    ##########################################################################
    # Set the ambient scattering of the spin glyphs
    ##########################################################################
    def RenAmbientUpdate(self, value, renWin):
        """
        Update the ambient property of MomActors.Spins and render the window.
        """
        if hasattr(self, "Spins"):
            self.Spins.GetProperty().SetAmbient(float(value * 0.02))
        renWin.Render()
        return

    ##########################################################################
    # Set the diffuse scattering of the spin glyphs
    ##########################################################################
    def RenDiffuseUpdate(self, value, renWin):
        """
        Updates the diffuse property of MomActors.Spins and renders the window.
        """
        if hasattr(self, "Spins"):
            self.Spins.GetProperty().SetDiffuse(float(value * 0.01))
        renWin.Render()
        return

    ##########################################################################
    # Set the specular scattering of the spin glyphs
    ##########################################################################
    def RenSpecularUpdate(self, value, renWin):
        """
        Updates the specular property of MomActors' Spins and renders the window.
        """
        if hasattr(self, "Spins"):
            self.Spins.GetProperty().SetSpecular(float(value * 0.01))
        renWin.Render()
        return

    ##########################################################################
    # Set the specular scattering of the spin glyphs
    ##########################################################################
    def RenSpecularPowerUpdate(self, value, renWin):
        """
        Updates the specular power of the Spins actor and renders the window.
        """
        if hasattr(self, "Spins"):
            self.Spins.GetProperty().SetSpecularPower(float(value))
        renWin.Render()
        return

    ##########################################################################
    # Set the PBR emission value of the spin glyphs
    ##########################################################################
    def PBREmissionUpdate(self, value, ren, renWin):
        """
        Update the emissive factor for the PBR material and render the window.
        """
        emvec = [float(value * 0.01), float(value * 0.01), float(value * 0.01)]
        if hasattr(self, "Spins"):
            self.Spins.GetProperty().SetEmissiveFactor(emvec)
        if hasattr(self, "Atoms"):
            self.Atoms.GetProperty().SetEmissiveFactor(emvec)
        renWin.Render()
        return

    ##########################################################################
    # Set the PBR occlusion value of the spin glyphs
    ##########################################################################
    def PBROcclusionUpdate(self, value, ren, renWin):
        """
        Updates the occlusion strength for Spins and Atoms actors and renders the window.
        """
        if hasattr(self, "Spins"):
            self.Spins.GetProperty().SetOcclusionStrength(float(value * 0.01))
        if hasattr(self, "Atoms"):
            self.Atoms.GetProperty().SetOcclusionStrength(float(value * 0.01))
        renWin.Render()
        return

    ##########################################################################
    # Set the PBR roughness value of the spin glyphs
    ##########################################################################
    def PBRRoughnessUpdate(self, value, renWin):
        """
        Updates the roughness property for Spins and Atoms actors and renders the window.
        """
        if hasattr(self, "Spins"):
            self.Spins.GetProperty().SetRoughness(float(value * 0.01))
        if hasattr(self, "Atoms"):
            self.Atoms.GetProperty().SetRoughness(float(value * 0.01))
        renWin.Render()
        return

    ##########################################################################
    # Set the PBR metallic value of the spin glyphs
    ##########################################################################
    def PBRMetallicUpdate(self, value, renWin):
        """
        Updates the metallic property of MomActors' Spins and Atoms and renders the window.
        """
        if hasattr(self, "Spins"):
            self.Spins.GetProperty().SetMetallic(float(value * 0.01))
        if hasattr(self, "Atoms"):
            self.Atoms.GetProperty().SetMetallic(float(value * 0.01))
        renWin.Render()
        return

    ##########################################################################
    # Set the size of the atoms via the slider
    ##########################################################################
    def ChangeAtomsSize(self, value):
        """
        Adjusts the size of atoms in the visualization by scaling the atom mapper.
        """
        self.AtomMapper.SetScaleFactor(1.00 * value / 10.0)
        return

    ##########################################################################
    # Set the size of the atoms via the slider
    ##########################################################################
    def ChangeAtomsOpaq(self, value):
        """
        Adjusts the opacity of atom actors based on the given value.
        """
        self.Atoms.GetProperty().SetOpacity(value * 0.01)
        return

    ##########################################################################
    # Set the quality of the atoms via the slider
    ##########################################################################
    def ChangeAtomsQuali(self, value):
        """
        Adjusts the resolution of the atom sphere visualization.
        """
        self.AtomSphere.SetThetaResolution(value)
        self.AtomSphere.SetPhiResolution(value)
        return