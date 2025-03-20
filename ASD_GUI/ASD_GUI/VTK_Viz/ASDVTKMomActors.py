"""@package ASDVTKMomActors
Wrapper class to add the VTK  actors for the visualization of UppASD data in
moment visualization mode.
It contains the needed data to add the actors, modify them, as well as some helper
functions to change them.

Author
----------
Jonathan Chico
"""
import vtk
import numpy as np
from vtk import vtkPoints
from vtk.util import numpy_support

class ASDMomActors():
    ############################################################################
    # @brief Main wrapper to add the needed actors for visualization of the moments
    # @details Main wrapper to add the needed actors for visualization of the moments.
    # Class that contains the data structures for creation of the glyphs
    # for the visualization of the magnetic moments.
    # It also has the capacity to create tessellations for the visualization of
    # volume vendering.
    # @author Jonathan Chico
    ############################################################################
    def Add_MomActors(self,ren,renWin,iren,ASDdata,window):
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


        ASDMomActors.timer_count=0
        ASDMomActors.camera_pos=np.zeros(3,dtype=np.float32)
        ASDMomActors.camera_focal=np.zeros(3,dtype=np.float32)
        ASDMomActors.camera_yaw=0.0
        ASDMomActors.camera_roll=0.0
        ASDMomActors.camera_pitch=0.0
        ASDMomActors.camera_azimuth=0.0
        ASDMomActors.camera_elevation=0.0
        ASDMomActors.kmc_disp=ASDdata.kmc_flag
        ASDMomActors.cluster_disp=ASDdata.cluster_flag
        ASDMomActors.glob_color=ASDdata.colors
        #-----------------------------------------------------------------------
        # Look up tables for colors
        #-----------------------------------------------------------------------
        # This is a diverging RWB color mapping based on the work of Kenneth
        # Moreland and with the vtk examples provided by Andrew Maclean
        if ASDdata.flag2D:
            ASDMomActors.lut = vtk.vtkLookupTable()
            num_colors = 256
            ASDMomActors.lut.SetNumberOfTableValues(num_colors)
            ASDMomActors.transfer_func = vtk.vtkColorTransferFunction()
            ASDMomActors.transfer_func.SetColorSpaceToDiverging()
            ASDMomActors.transfer_func.AddRGBPoint(0, 0.230, 0.299, 0.754)
            ASDMomActors.transfer_func.AddRGBPoint(1, 0.706, 0.016, 0.150)
            for ii,ss in enumerate([float(xx)/float(num_colors) for xx in range(num_colors)]):
                cc = ASDMomActors.transfer_func.GetColor(ss)
                ASDMomActors.lut.SetTableValue(ii, cc[0], cc[1], cc[2], 1.0)
            ASDMomActors.lut.Build()
        else:
            ASDMomActors.lut = vtk.vtkLookupTable()
            num_colors = 256
            ASDMomActors.lut.SetNumberOfTableValues(num_colors)
            ASDMomActors.transfer_func = vtk.vtkColorTransferFunction()
            ASDMomActors.transfer_func.SetColorSpaceToDiverging()
            ASDMomActors.transfer_func.AddRGBPoint(-0, 0.230, 0.299, 0.754)
            ASDMomActors.transfer_func.AddRGBPoint( 1, 0.706, 0.016, 0.150)
            for ii,ss in enumerate([float(xx)/float(num_colors) for xx in range(num_colors)]):
                cc = ASDMomActors.transfer_func.GetColor(ss)
                ASDMomActors.lut.SetTableValue(ii, cc[0], cc[1], cc[2], 1.0)
            ASDMomActors.lut.Build()
        #-----------------------------------------------------------------------
        # Data structures for the generation of the smooth grid
        #-----------------------------------------------------------------------
        # Passing the data from the full system to the PolyData
        ASDMomActors.src=vtk.vtkPolyData()
        ASDMomActors.src.SetPoints(ASDdata.coord)
        ASDMomActors.src.GetPointData().SetScalars(ASDdata.colors[2])
        ASDMomActors.src.GetPointData().SetVectors(ASDdata.moments)
        scalar_range = ASDMomActors.src.GetScalarRange()
        #-----------------------------------------------------------------------
        # Finding useful geometrical information of the sample
        #-----------------------------------------------------------------------
        # Finding the middle of the sample
        # Also making sure that if the sample is 2D one has no problem with bounds
        # this is mostly useful if splatters are used
        (ASDMomActors.xmin,ASDMomActors.xmax,ASDMomActors.ymin,ASDMomActors.ymax,\
        ASDMomActors.zmin,ASDMomActors.zmax)= ASDMomActors.src.GetBounds()
        if ASDMomActors.xmin==ASDMomActors.xmax:
            ASDMomActors.xmin=0.0
            ASDMomActors.xmax=1.0
        if ASDMomActors.ymin==ASDMomActors.ymax:
            ASDMomActors.ymin=0.0
            ASDMomActors.ymax=1.0
        if ASDMomActors.zmin==ASDMomActors.zmax:
            ASDMomActors.zmin=0.0
            ASDMomActors.zmax=1.0
        ASDMomActors.xmid = (ASDMomActors.xmin+ASDMomActors.xmax)*0.5
        ASDMomActors.ymid = (ASDMomActors.ymin+ASDMomActors.ymax)*0.5
        ASDMomActors.zmid = (ASDMomActors.zmin+ASDMomActors.zmax)*0.5
        ASDMomActors.height=max(ASDMomActors.xmax,ASDMomActors.ymax,ASDMomActors.zmax)*1.75
        self.dist_x=np.absolute(ASDMomActors.xmax-ASDMomActors.xmin)
        self.dist_y=np.absolute(ASDMomActors.ymax-ASDMomActors.ymin)
        self.dist_z=np.absolute(ASDMomActors.zmax-ASDMomActors.zmin)
        ASDMomActors.camera_pos[0]=ASDMomActors.xmid
        ASDMomActors.camera_pos[1]=ASDMomActors.ymid
        ASDMomActors.camera_pos[2]=ASDMomActors.height
        ASDMomActors.camera_focal[0]=ASDMomActors.xmid
        ASDMomActors.camera_focal[1]=ASDMomActors.ymid
        ASDMomActors.camera_focal[2]=ASDMomActors.zmid
        #-----------------------------------------------------------------------
        # Data structures for the spins
        #-----------------------------------------------------------------------
        # Passing the data from the full system to the PolyData
        ASDMomActors.src_spins=vtk.vtkPolyData()
        ASDMomActors.src_spins.SetPoints(ASDdata.coord)
        ASDMomActors.src_spins.GetPointData().SetScalars(ASDdata.colors[2])
        ASDMomActors.src_spins.GetPointData().SetVectors(ASDdata.moments)
        scalar_range_spins = ASDMomActors.src_spins.GetScalarRange()
        #-----------------------------------------------------------------------
        ASDMomActors.ungrid=vtk.vtkUnstructuredGrid()
        ASDMomActors.ungrid.SetPoints(ASDdata.coord)
        ASDMomActors.ungrid.GetPointData().SetScalars(ASDMomActors.src.GetPointData().GetScalars())
        #writer = vtk.vtkUnstructuredGridWriter()
        #writer.SetFileName('out.vtk')
        #writer.SetInputData(ASDMomActors.ungrid)
        #writer.Write()
        #ASDMomActors.ungrid.SetPoints(ASDdata.coord)
        #ASDMomActors.ungrid.GetPointData().SetScalars(ASDdata.colors[2])
        #ASDMomActors.ungrid.GetPointData().SetVectors(ASDdata.moments)

        #-----------------------------------------------------------------------
        # The delaunay tessellation seems to be the best way to transform the point cloud
        # to a surface for volume rendering, the problem is that it is too slow for large
        # data sets, meaning that the best option is first to prune out the data to ensure
        # that one has a manageable number of data points over which to do the construction
        # surface reconstruction and splatter techniques also can be used to generate something
        # akin to the kind of surfaces we want. The issue is that they transform the data to a
        # regular mesh by default. And thus it is a problem for most kind of systems
        if ASDdata.flag2D:
            # Passing the data to generate a triangulation of the data
            ASDMomActors.MagDensMethod = vtk.vtkDelaunay2D()
            ASDMomActors.MagDensMethod.SetInputData(ASDMomActors.src)
            ASDMomActors.MagDensMethod.BoundingTriangulationOff()
            ASDMomActors.MagDensMethod.SetTolerance(0.005)
            # Time the execution of the delaunay tessellation
            SM_timer = vtk.vtkExecutionTimer()
            SM_timer.SetFilter(ASDMomActors.MagDensMethod)
            ASDMomActors.MagDensMethod.Update()
            SM = SM_timer.GetElapsedWallClockTime()
            print ("2D Delaunay:", SM)
            # Creating the mapper for the smooth surfaces
            ASDMomActors.MagDensMap = vtk.vtkDataSetMapper()
            ASDMomActors.MagDensMap.SetScalarRange(scalar_range)
            ASDMomActors.MagDensMap.SetInputConnection(ASDMomActors.MagDensMethod.GetOutputPort())
            ASDMomActors.MagDensMap.SetLookupTable(ASDMomActors.lut)
            ASDMomActors.MagDensMap.SetColorModeToMapScalars()
            ASDMomActors.MagDensMap.Update()
            # Creating the actor for the smooth surfaces
            ASDMomActors.MagDensActor = vtk.vtkLODActor()
            ASDMomActors.MagDensActor.SetMapper(ASDMomActors.MagDensMap)
            ASDMomActors.MagDensActor.GetProperty().SetOpacity(0.75)
            ASDMomActors.MagDensActor.GetProperty().EdgeVisibilityOff()
            if window.DensBox.isChecked():
                ASDMomActors.MagDensActor.VisibilityOn()
            else:
                ASDMomActors.MagDensActor.VisibilityOff()
        else:
            #-------------------------------------------------------------------
            # Setting the parameters for the visualization of 3D structures with
            # splatters
            #-------------------------------------------------------------------
            ASDMomActors.MagDensMethod = vtk.vtkGaussianSplatter()
            ASDMomActors.MagDensMethod.SetInputData(ASDMomActors.src)
            #print( ASDMomActors.MagDensMethod.GetSampleDimensions())
            #ASDMomActors.MagDensMethod.SetSampleDimensions([10,10,2])
            #print( ASDMomActors.MagDensMethod.GetSampleDimensions())
            #-------------------------------------------------------------------
            # Options for the Gaussian splatter. These determine the quality of the
            # rendering, increasing the radius smoothens out the volume but performance
            # decreases rapidly
            #-------------------------------------------------------------------
            dist=(np.asarray(ASDMomActors.src.GetPoint(0))-np.asarray(ASDMomActors.src.GetPoint(1)))
            norm=np.sqrt(dist.dot(dist))
            if norm<10:
                rad_fac=0.040
            else:
                rad_fac=0.40
            ASDMomActors.MagDensMethod.SetRadius(rad_fac)
            ASDMomActors.MagDensMethod.ScalarWarpingOn()
            #-------------------------------------------------------------------
            # The exponent factor determines how fast the gaussian splatter decay
            # they again can be used to improve quality at the sake of rendering time
            #-------------------------------------------------------------------
            ASDMomActors.MagDensMethod.SetExponentFactor(-10)
            ASDMomActors.MagDensMethod.NormalWarpingOn()
            ASDMomActors.MagDensMethod.SetEccentricity(10)
            #-------------------------------------------------------------------
            # The Null value can be used to try to eliminate contributions not belonging
            # to the actual sample
            #-------------------------------------------------------------------
            ASDMomActors.MagDensMethod.SetNullValue(-10)
            #-------------------------------------------------------------------
            # Set the actual size of the rendering model
            #-------------------------------------------------------------------
            ASDMomActors.MagDensMethod.SetModelBounds(ASDMomActors.xmin,ASDMomActors.xmax,\
            ASDMomActors.ymin,ASDMomActors.ymax,ASDMomActors.zmin,ASDMomActors.zmax)
            # This should get rid of the problems when trying to map very thin structures in 2D
            if self.dist_x==min(self.dist_x,self.dist_y,self.dist_z) and self.dist_x<3:
                ASDMomActors.MagDensMethod.SetSampleDimensions(3,int(ASDMomActors.ymax),int(ASDMomActors.zmax))
            elif self.dist_y==min(self.dist_x,self.dist_y,self.dist_z) and self.dist_y<3:
                ASDMomActors.MagDensMethod.SetSampleDimensions(int(ASDMomActors.xmax),3,int(ASDMomActors.zmax))
            elif self.dist_z==min(self.dist_x,self.dist_y,self.dist_z) and self.dist_z<3:
                ASDMomActors.MagDensMethod.SetSampleDimensions(int(ASDMomActors.xmax),int(ASDMomActors.ymax),3)
            # Timming for the execution of the volume creation
            SP_timer = vtk.vtkExecutionTimer()
            SP_timer.SetFilter(ASDMomActors.MagDensMethod)
            ASDMomActors.MagDensMethod.Update()
            SP = SP_timer.GetElapsedWallClockTime()
            print ("3D vtkGaussianSplatter Method:", SP)
            # Scalar opacities
            funcOpacityScalar = vtk.vtkPiecewiseFunction()
            funcOpacityScalar.AddPoint(-1.00,0.00)
            funcOpacityScalar.AddPoint( 0.00,0.05)
            funcOpacityScalar.AddPoint( 0.50,0.50)
            funcOpacityScalar.AddPoint( 0.75,1.00)
            # Gradient opacities
            volumeGradientOpacity = vtk.vtkPiecewiseFunction()
            volumeGradientOpacity.AddPoint(0.000,0.0)
            volumeGradientOpacity.AddPoint(0.001,1.0)
            volumeGradientOpacity.AddPoint(1.000,1.0)
            # Volume properties
            ASDMomActors.volumeProperty = vtk.vtkVolumeProperty()
            ASDMomActors.volumeProperty.SetColor(ASDMomActors.transfer_func)
            ASDMomActors.volumeProperty.SetInterpolationTypeToLinear()
            ASDMomActors.volumeProperty.SetAmbient(0.6)
            ASDMomActors.volumeProperty.SetDiffuse(0.6)
            ASDMomActors.volumeProperty.SetSpecular(0.1)
            ASDMomActors.volumeProperty.SetGradientOpacity(volumeGradientOpacity)
            ASDMomActors.volumeProperty.SetScalarOpacity(funcOpacityScalar)
            # Volume Mapper
            #ASDMomActors.MagDensMethod = vtk.vtkDataSetTriangleFilter()
            #ASDMomActors.MagDensMethod.SetInputData(ASDMomActors.ungrid)
            #ASDMomActors.MagDensMethod.Update()

            ASDMomActors.MagDensMap = vtk.vtkSmartVolumeMapper()
            #ASDMomActors.MagDensMap = vtk.vtkUnstructuredGridVolumeRayCastMapper()
            ASDMomActors.MagDensMap.SetInputConnection(ASDMomActors.MagDensMethod.GetOutputPort())
            #ASDMomActors.MagDensMap.Update()

            # Volume Actor
            ASDMomActors.MagDensActor = vtk.vtkVolume()
            ASDMomActors.MagDensActor.SetMapper(ASDMomActors.MagDensMap)
            ASDMomActors.MagDensActor.SetProperty(self.volumeProperty)
            if window.DensBox.isChecked():
                ASDMomActors.MagDensActor.VisibilityOn()
            else:
                ASDMomActors.MagDensActor.VisibilityOff()
            ASDMomActors.MagDensActor.Update()
        #-----------------------------------------------------------------------
        # Data structures for the contours
        #-----------------------------------------------------------------------
        # Define the contour filters
        contours = vtk.vtkContourFilter()
        contours.SetInputConnection(ASDMomActors.MagDensMethod.GetOutputPort())
        # This generates the contours, it will do 5 between the -1 and 0.5 range
        cont_num=5
        range_cont=(-1,0.5)
        contours.GenerateValues(cont_num,range_cont)
        # Map the contours to graphical primitives
        contMapper = vtk.vtkPolyDataMapper()
        contMapper.SetInputConnection(contours.GetOutputPort())
        contMapper.SetScalarVisibility(False) # colored contours
        contMapper.SetScalarRange(scalar_range)
        # Create an actor for the contours
        ASDMomActors.contActor = vtk.vtkLODActor()
        ASDMomActors.contActor.SetMapper(contMapper)
        ASDMomActors.contActor.GetProperty().SetColor(0, 0, 0)
        ASDMomActors.contActor.GetProperty().SetLineWidth(1.0)
        ASDMomActors.contActor.VisibilityOff()
        
        #-----------------------------------------------------------------------
        # Setting information of the spins
        #-----------------------------------------------------------------------
        # Create vectors
        ASDMomActors.spinarrow = vtk.vtkArrowSource()
        ASDMomActors.spinarrow.SetTipRadius(0.20)
        ASDMomActors.spinarrow.SetShaftRadius(0.10)
        ASDMomActors.spinarrow.SetTipResolution(12)
        ASDMomActors.spinarrow.SetShaftResolution(12)
        
        ASDMomActors.spinarrowtriangles = vtk.vtkTriangleFilter()
        ASDMomActors.spinarrowtriangles.SetInputConnection(ASDMomActors.spinarrow.GetOutputPort())

        # Calculate normals for shading
        ASDMomActors.spinarrownormals = vtk.vtkPolyDataNormals()
        ASDMomActors.spinarrownormals.SetInputConnection(ASDMomActors.spinarrowtriangles.GetOutputPort())

        ASDMomActors.spinarrowtcoords = vtk.vtkTextureMapToCylinder()
        ASDMomActors.spinarrowtcoords.SetInputConnection(ASDMomActors.spinarrownormals.GetOutputPort())
        ASDMomActors.spinarrowtcoords.PreventSeamOn()

        ASDMomActors.spinarrowtangents = vtk.vtkPolyDataTangents()
        ASDMomActors.spinarrowtangents.SetInputConnection(ASDMomActors.spinarrowtcoords.GetOutputPort())

        # Create the mapper for the spins
        ASDMomActors.SpinMapper = vtk.vtkGlyph3DMapper()
        ASDMomActors.SpinMapper.SetSourceConnection(ASDMomActors.spinarrowtangents.GetOutputPort())

        ASDMomActors.SpinMapper.SetInputData(ASDMomActors.src_spins)
        ASDMomActors.SpinMapper.SetScalarRange(scalar_range_spins)
        ASDMomActors.SpinMapper.SetScaleFactor(0.50)
        ASDMomActors.SpinMapper.SetScaleModeToNoDataScaling()
        ASDMomActors.SpinMapper.SetLookupTable(self.lut)
        ASDMomActors.SpinMapper.SetColorModeToMapScalars()
        ASDMomActors.SpinMapper.Update()
        
        # Define the vector actor for the spins
        ASDMomActors.Spins = vtk.vtkActor()
        #ASDMomActors.Spins = vtk.vtkLODActor()
        ASDMomActors.Spins.SetMapper(ASDMomActors.SpinMapper)
        #ASDMomActors.Spins.GetProperty().SetInterpolationToPBR()
        ASDMomActors.Spins.GetProperty().SetInterpolationToGouraud()
        ASDMomActors.Spins.GetProperty().SetSpecular(0.4)
        ASDMomActors.Spins.GetProperty().SetSpecularPower(80)
        ASDMomActors.Spins.GetProperty().SetAmbient(0.6)
        ASDMomActors.Spins.GetProperty().SetDiffuse(0.4)
        ASDMomActors.Spins.GetProperty().SetEdgeTint(0.0,0.0,0.0)
        if window.SpinsBox.isChecked():
            ASDMomActors.Spins.VisibilityOn()
        else:
            ASDMomActors.Spins.VisibilityOff()


        #-----------------------------------------------------------------------
        # Setting information for the atoms
        #-----------------------------------------------------------------------
        # Create vectors
        ASDMomActors.AtomSphere = vtk.vtkTexturedSphereSource()
        ASDMomActors.AtomSphere.SetRadius(0.10)
        ASDMomActors.AtomSphere.SetThetaResolution(12)
        ASDMomActors.AtomSphere.SetPhiResolution(12)

        
        # Create the mapper for the atoms
        ASDMomActors.AtomMapper = vtk.vtkGlyph3DMapper()
        ASDMomActors.AtomMapper.SetSourceConnection(ASDMomActors.AtomSphere.GetOutputPort())
        ASDMomActors.AtomMapper.SetInputData(ASDMomActors.src_spins)
        ASDMomActors.AtomMapper.SetScaleFactor(1.00)
        ASDMomActors.AtomMapper.SetScaleModeToNoDataScaling()
        ASDMomActors.AtomMapper.ScalarVisibilityOff()
        #ASDMomActors.AtomMapper.SetLookupTable(self.lut)
        #ASDMomActors.AtomMapper.SetColorModeToMapScalars()
        ASDMomActors.AtomMapper.Update()
        # Define the sphere actor for the atoms
        colors = vtk.vtkNamedColors()
        ASDMomActors.Atoms = vtk.vtkLODActor()
        ASDMomActors.Atoms.SetMapper(ASDMomActors.AtomMapper)
        ASDMomActors.Atoms.GetProperty().SetInterpolationToGouraud()
        #ASDMomActors.Atoms.GetProperty().SetInterpolationToPBR()
        ASDMomActors.Atoms.GetProperty().SetAmbient(0.8)
        ASDMomActors.Atoms.GetProperty().SetDiffuse(0.8)
        ASDMomActors.Atoms.GetProperty().SetColor(colors.GetColor3d("Silver"))
        if window.AtomsBox.isChecked():
            ASDMomActors.Atoms.VisibilityOn()
        else:
            ASDMomActors.Atoms.VisibilityOff()

        #-----------------------------------------------------------------------
        # Setting information for the skybox actor
        #-----------------------------------------------------------------------
        ASDMomActors.SkyBox = vtk.vtkSkybox()
        ASDMomActors.SkyBox.VisibilityOff()
        ren.AddActor(ASDMomActors.SkyBox)

        if (ASDdata.kmc_flag):
            #-------------------------------------------------------------------
            # Setting data structures for the KMC particle visualization
            #-------------------------------------------------------------------
            ASDMomActors.KMC_src=vtk.vtkPolyData()
            ASDMomActors.KMC_src.SetPoints(ASDdata.coord_KMC)
            # Atom sphere
            KMC_part = vtk.vtkSphereSource()
            KMC_part.SetRadius(1.75)
            KMC_part.SetThetaResolution(12)
            KMC_part.SetPhiResolution(12)
            # Atom glyph
            KMC_part_mapper = vtk.vtkGlyph3DMapper()
            KMC_part_mapper.SetInputData(ASDMomActors.KMC_src)
            KMC_part_mapper.SetSourceConnection(KMC_part.GetOutputPort())
            KMC_part_mapper.SetScaleFactor(0.5)
            KMC_part_mapper.ClampingOn()
            KMC_part_mapper.SetScaleModeToNoDataScaling()
            KMC_part_mapper.SetColorModeToMapScalars()
            KMC_part_mapper.Update()
            # Atoms actors
            ASDMomActors.KMC_part_actor = vtk.vtkLODActor()
            ASDMomActors.KMC_part_actor.SetMapper(KMC_part_mapper)
            ASDMomActors.KMC_part_actor.GetProperty().SetOpacity(0.9)
            ASDMomActors.KMC_part_actor.GetProperty().SetColor(0.0, 0.0, 1.0)
            ASDMomActors.KMC_part_actor.GetProperty().EdgeVisibilityOn()
            ASDMomActors.KMC_part_actor.GetProperty().SetEdgeColor(0,0,0)
        #-----------------------------------------------------------------------
        # Setting information of the renderer
        #-----------------------------------------------------------------------
        # Define the renderer
        # Add the actors to the scene
        if ASDdata.flag2D:
            ren.AddActor(ASDMomActors.MagDensActor)
        else:
            ren.AddViewProp(ASDMomActors.MagDensActor)
        ren.AddActor(ASDMomActors.Spins)
        #ren.AddActor(ASDMomActors.vector)
        ren.AddActor(ASDMomActors.contActor)
        ren.AddActor(ASDMomActors.Atoms)
        #If the KMC particles are present add them to the renderer
        if ASDdata.kmc_flag:
            ren.AddActor(ASDMomActors.KMC_part_actor)
        # Defining the camera directions
        ren.GetActiveCamera().Azimuth(ASDMomActors.camera_azimuth)
        ren.GetActiveCamera().Elevation(ASDMomActors.camera_elevation)
        ren.GetActiveCamera().Yaw(ASDMomActors.camera_yaw)
        ren.GetActiveCamera().Roll(ASDMomActors.camera_roll)
        ren.GetActiveCamera().Pitch(ASDMomActors.camera_pitch)
        ren.GetActiveCamera().SetFocalPoint(ASDMomActors.camera_focal)
        ren.GetActiveCamera().SetPosition(ASDMomActors.camera_pos)
        ren.GetActiveCamera().SetViewUp(0,1,0)
        #-----------------------------------------------------------------------
        # Start the renderer
        #-----------------------------------------------------------------------
        iren.Start()
        renWin.Render()
        print(' Done')
        return;

    ############################################################################
    # @brief Update the magnetic moments for the visualization
    # @author Jonathan Chico
    ############################################################################
    def UpdateMoments(self,window,ASDdata,ASDGenActors,renWin, is_interactive=False):
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
        from vtk import vtkPoints
        from vtk.util import numpy_support
        #-----------------------------------------------------------------------
        # Read the actual data of the magnetic moments
        #-----------------------------------------------------------------------
        print('Reading data for time step:',window.current_time)
        print(self.__class__.__name__)
        (ASDdata.moments,ASDdata.colors,ASDdata.number_time_steps,ASDdata.time_sep)=\
        ASDdata.readVectorsData(ASDdata.MagFile,window.current_time,\
        ASDdata.nrAtoms,ASDdata.number_time_steps)
        #-----------------------------------------------------------------------
        # Update the colors
        #-----------------------------------------------------------------------
        if window.DensX.isChecked():
            ASDMomActors.src.GetPointData().SetScalars(ASDdata.colors[0])
        if window.SpinX.isChecked():
            ASDMomActors.src_spins.GetPointData().SetScalars(ASDdata.colors[0])
        if window.DensY.isChecked():
            ASDMomActors.src.GetPointData().SetScalars(ASDdata.colors[1])
        if window.SpinY.isChecked():
            ASDMomActors.src_spins.GetPointData().SetScalars(ASDdata.colors[1])
        if window.DensZ.isChecked():
            ASDMomActors.src.GetPointData().SetScalars(ASDdata.colors[2])
        if window.SpinZ.isChecked():
            ASDMomActors.src_spins.GetPointData().SetScalars(ASDdata.colors[2])
        #-----------------------------------------------------------------------
        # Update the vectors
        #-----------------------------------------------------------------------
        ASDMomActors.src.GetPointData().SetVectors(ASDdata.moments)
        ASDMomActors.src_spins.GetPointData().SetVectors(ASDdata.moments)
        #-----------------------------------------------------------------------
        # If coordinate file is animated, update coordinates
        #-----------------------------------------------------------------------
        if ASDdata.nrAtoms<ASDdata.full_coord.shape[0]:
            t_off = window.current_time * ASDdata.nrAtoms
            ASDdata.coord.SetData(numpy_support.numpy_to_vtk(ASDdata.full_coord[t_off:t_off+ASDdata.nrAtoms]))
            ASDMomActors.src.SetPoints(ASDdata.coord)
        #-----------------------------------------------------------------------
        # Update the general actors
        #-----------------------------------------------------------------------
        window.ProgressBar.setValue(int((window.current_time-1)*100/(ASDdata.number_time_steps-1)))
        window.ProgressLabel.setText(f'   {int(window.ProgressBar.value())}%')
        time_label = f"{float(window.TimeStepLineEdit.text()) * ASDdata.time_sep[window.current_time-1] * 1e9: 4.2f} ns"
        ASDGenActors.time_label.SetInput(time_label)
        #-----------------------------------------------------------------------
        # Render the window
        #-----------------------------------------------------------------------
        renWin.Render()
        return
