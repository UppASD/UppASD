""" @package ASDVTKNeighActors
This is the routine where the actual VTK actors are created for the
Neighbour visualization mode. It creates a series of spheres to indicate
the atomic positions in the lattice, as well as a visual indication of the
neighbour cloud, where the color of the neighbour actors is given by the
magnitude of the exchange interaction.
It also contains a mode to add vectors for the information of the DM neighbour
interaction, with the color of the glyphs being dictated by the magnitude of the
DM interaction. The position of the vectors is located at the neighbour atom.

Author
----------
Jonathan Chico
"""
# pylint: disable=invalid-name, no-name-in-module, no-member
import vtk
import numpy as np


class ASDNeighActors:
    
    def __init__(self):
        self.active = False
    
    ##########################################################################
    # @brief This defines the actors for the visualization of neighbours from the struct file
    # @details This defines the actors for the visualization of neighbours from the struct file.
    # It allows for the visualization of the Heisenberg exchange interaction and
    # of the DM vectors.
    # @author Jonathan Chico
    ##########################################################################
    def Add_NeighActors(self, ren, renWin, iren, ASDdata, mode):
        """This defines the actors for the visualization of neighbours from the struct file.
        It allows for the visualization of the Heisenberg exchange interaction and
        of the DM vectors.
        Args:
            ren: current renderer for the VTK visualization.
            renWin: VTK window where the visualization is performed.
            iren: current window interactor for the VTK visualization
            ASData: class containing the data read from the UppASD simulations.
            mode: value indicating whether one is visualizing
            the Heisenberg exchange or DMI vectors (1/2)

        Author
        ----------
        Jonathan Chico
        """
        ASDNeighActors.timer_count = 0
        ASDNeighActors.camera_pos = np.zeros(3, dtype=np.float32)
        ASDNeighActors.camera_focal = np.zeros(3, dtype=np.float32)
        ASDNeighActors.camera_yaw = 0.0
        ASDNeighActors.camera_roll = 0.0
        ASDNeighActors.camera_pitch = 0.0
        ASDNeighActors.camera_azimuth = 0.0
        ASDNeighActors.camera_elevation = 0.0
        # -----------------------------------------------------------------------
        # Data structures for the atoms in the neighbour map
        # -----------------------------------------------------------------------
        AtomGrid = vtk.vtkPolyData()
        AtomGrid.SetPoints(ASDdata.coord)
        ASDNeighActors.SLMax = AtomGrid.GetNumberOfPoints()
        # -----------------------------------------------------------------------
        # Atom sphere
        # -----------------------------------------------------------------------
        Atom = vtk.vtkSphereSource()
        Atom.SetRadius(0.50)
        Atom.SetThetaResolution(10)
        Atom.SetPhiResolution(10)
        # -----------------------------------------------------------------------
        # Atom glyph
        # -----------------------------------------------------------------------
        Atoms = vtk.vtkGlyph3DMapper()
        Atoms.SetInputData(AtomGrid)
        Atoms.SetSourceConnection(Atom.GetOutputPort())
        Atoms.SetScaleFactor(0.5)
        Atoms.ClampingOn()
        Atoms.SetScaleModeToNoDataScaling()
        Atoms.SetColorModeToMapScalars()
        Atoms.Update()
        # -----------------------------------------------------------------------
        # Atoms actors
        # -----------------------------------------------------------------------
        ASDNeighActors.AtomsActor = vtk.vtkLODActor()
        ASDNeighActors.AtomsActor.SetMapper(Atoms)
        ASDNeighActors.AtomsActor.GetProperty().SetOpacity(0.9)
        ASDNeighActors.AtomsActor.GetProperty().SetColor(0.5, 0.5, 0.5)
        # -----------------------------------------------------------------------
        # Data structures for the neighbours in the neighbour mapper
        # -----------------------------------------------------------------------
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
        # -----------------------------------------------------------------------
        # Grid for neighbours
        # -----------------------------------------------------------------------
        ASDNeighActors.NeighGrid = vtk.vtkPolyData()
        ASDNeighActors.NeighGrid.SetPoints(ASDdata.neighs)
        ASDNeighActors.NeighGrid.GetPointData().SetScalars(ASDdata.neigh_colors)
        if mode == 2:
            ASDNeighActors.NeighGrid.GetPointData().SetVectors(ASDdata.DM_vectors)
        ASDNeighActors.NumNeigh = ASDNeighActors.NeighGrid.GetNumberOfPoints()
        scalar_range = ASDNeighActors.NeighGrid.GetScalarRange()
        # -----------------------------------------------------------------------
        # Finding useful geometrical information of the sample
        # Finding the middle of the sample
        # Also making sure that if the sample is 2D one has no problem with boudings
        # this is mostly useful if splatters are used
        # -----------------------------------------------------------------------
        (
            ASDNeighActors.xmin,
            ASDNeighActors.xmax,
            ASDNeighActors.ymin,
            ASDNeighActors.ymax,
            ASDNeighActors.zmin,
            ASDNeighActors.zmax,
        ) = ASDNeighActors.NeighGrid.GetBounds()
        if ASDNeighActors.xmin == ASDNeighActors.xmax:
            ASDNeighActors.xmin = 0.0
            ASDNeighActors.xmax = 1.0
        if ASDNeighActors.ymin == ASDNeighActors.ymax:
            ASDNeighActors.ymin = 0.0
            ASDNeighActors.ymax = 1.0
        if ASDNeighActors.zmin == ASDNeighActors.zmax:
            ASDNeighActors.zmin = 0.0
            ASDNeighActors.zmax = 1.0
        ASDNeighActors.xmid = (ASDNeighActors.xmin + ASDNeighActors.xmax) * 0.5
        ASDNeighActors.ymid = (ASDNeighActors.ymin + ASDNeighActors.ymax) * 0.5
        ASDNeighActors.zmid = (ASDNeighActors.zmin + ASDNeighActors.zmax) * 0.5
        ASDNeighActors.height = (
            max(ASDNeighActors.xmax,
                ASDNeighActors.ymax,
                ASDNeighActors.zmax) * 1.75
        )
        self.dist_x = np.absolute(ASDNeighActors.xmax - ASDNeighActors.xmin)
        self.dist_y = np.absolute(ASDNeighActors.ymax - ASDNeighActors.ymin)
        self.dist_z = np.absolute(ASDNeighActors.zmax - ASDNeighActors.zmin)
        ASDNeighActors.camera_pos[0] = ASDNeighActors.xmid
        ASDNeighActors.camera_pos[1] = ASDNeighActors.ymid
        ASDNeighActors.camera_pos[2] = ASDNeighActors.height
        ASDNeighActors.camera_focal[0] = ASDNeighActors.xmid
        ASDNeighActors.camera_focal[1] = ASDNeighActors.ymid
        ASDNeighActors.camera_focal[2] = ASDNeighActors.zmid
        # -----------------------------------------------------------------------
        # Neighbour glyphs
        # -----------------------------------------------------------------------
        if mode == 1:
            ASDNeighActors.NeighGlyphs = vtk.vtkSphereSource()
            ASDNeighActors.NeighGlyphs.SetRadius(0.50)
            ASDNeighActors.NeighGlyphs.SetThetaResolution(40)
            ASDNeighActors.NeighGlyphs.SetPhiResolution(40)
        else:
            ASDNeighActors.NeighGlyphs = vtk.vtkArrowSource()
            ASDNeighActors.NeighGlyphs.SetTipRadius(0.20)
            ASDNeighActors.NeighGlyphs.SetShaftRadius(0.10)
            ASDNeighActors.NeighGlyphs.SetTipResolution(40)
            ASDNeighActors.NeighGlyphs.SetShaftResolution(40)
        # -----------------------------------------------------------------------
        # Glyph source
        # -----------------------------------------------------------------------
        ASDNeighActors.NeighGlyph3D = vtk.vtkGlyph3D()
        ASDNeighActors.NeighGlyph3D.SetSourceConnection(
            ASDNeighActors.NeighGlyphs.GetOutputPort()
        )
        ASDNeighActors.NeighGlyph3D.SetVectorModeToUseNormal()
        ASDNeighActors.NeighGlyph3D.SetInputData(ASDNeighActors.NeighGrid)
        if mode == 1:
            ASDNeighActors.NeighGlyph3D.SetScaleFactor(1.05)
        elif mode == 2:
            ASDNeighActors.NeighGlyph3D.SetScaleFactor(1.25)
        ASDNeighActors.NeighGlyph3D.SetColorModeToColorByScalar()
        ASDNeighActors.NeighGlyph3D.SetScaleModeToDataScalingOff()
        if mode == 2:
            ASDNeighActors.NeighGlyph3D.SetVectorModeToUseVector()
        ASDNeighActors.NeighGlyph3D.Update()
        # -----------------------------------------------------------------------
        # Set up Neighbour glyphs
        # -----------------------------------------------------------------------
        ASDNeighActors.NeighMapper = vtk.vtkPolyDataMapper()
        ASDNeighActors.NeighMapper.SetInputConnection(
            ASDNeighActors.NeighGlyph3D.GetOutputPort()
        )
        ASDNeighActors.NeighMapper.SetScalarRange(scalar_range)
        ASDNeighActors.NeighMapper.SetLookupTable(self.lut)
        ASDNeighActors.NeighMapper.SetColorModeToMapScalars()
        ASDNeighActors.NeighMapper.Update()
        # -----------------------------------------------------------------------
        # Neighbour actors
        # -----------------------------------------------------------------------
        ASDNeighActors.NeighActor = vtk.vtkLODActor()
        ASDNeighActors.NeighActor.SetMapper(ASDNeighActors.NeighMapper)
        ASDNeighActors.NeighActor.GetProperty().SetSpecular(0.3)
        ASDNeighActors.NeighActor.GetProperty().SetSpecularPower(60)
        ASDNeighActors.NeighActor.GetProperty().SetAmbient(0.2)
        ASDNeighActors.NeighActor.GetProperty().SetDiffuse(0.8)
        # -----------------------------------------------------------------------
        # Grid for the center atom
        # -----------------------------------------------------------------------
        ASDNeighActors.CenterGrid = vtk.vtkPolyData()
        ASDNeighActors.CenterGrid.SetPoints(ASDdata.atomCenter)
        ASDNeighActors.CenterGrid.Modified()
        # -----------------------------------------------------------------------
        # Source for the center atom, a sphere
        # -----------------------------------------------------------------------
        ASDNeighActors.CenterSource = vtk.vtkSphereSource()
        ASDNeighActors.CenterSource.SetRadius(0.50)
        ASDNeighActors.CenterSource.SetThetaResolution(40)
        ASDNeighActors.CenterSource.SetPhiResolution(40)
        # -----------------------------------------------------------------------
        # Mapper for the center actor
        # -----------------------------------------------------------------------
        ASDNeighActors.Center = vtk.vtkGlyph3DMapper()
        ASDNeighActors.Center.SetInputData(ASDNeighActors.CenterGrid)
        ASDNeighActors.Center.SetSourceConnection(
            ASDNeighActors.CenterSource.GetOutputPort()
        )
        ASDNeighActors.Center.SetScaleFactor(1.25)
        ASDNeighActors.Center.SetScaleModeToNoDataScaling()
        ASDNeighActors.Center.Update()
        # -----------------------------------------------------------------------
        # Actor for the center atom
        # -----------------------------------------------------------------------
        ASDNeighActors.CenterActor = vtk.vtkLODActor()
        ASDNeighActors.CenterActor.SetMapper(ASDNeighActors.Center)
        ASDNeighActors.CenterActor.GetProperty().SetColor(0.4, 0.8, 0.4)
        ASDNeighActors.CenterActor.GetProperty().SetSpecular(0.3)
        ASDNeighActors.CenterActor.GetProperty().SetSpecularPower(60)
        ASDNeighActors.CenterActor.GetProperty().SetAmbient(0.2)
        ASDNeighActors.CenterActor.GetProperty().SetDiffuse(0.8)
        # -----------------------------------------------------------------------
        # Defining the camera directions
        # -----------------------------------------------------------------------
        ren.GetActiveCamera().Azimuth(0)
        ren.GetActiveCamera().Elevation(0)
        ren.GetActiveCamera().SetFocalPoint(
            ASDNeighActors.xmid, ASDNeighActors.ymid, ASDNeighActors.zmid
        )
        ren.GetActiveCamera().SetPosition(
            ASDNeighActors.xmid, ASDNeighActors.ymid, ASDNeighActors.height
        )
        ren.GetActiveCamera().Azimuth(ASDNeighActors.camera_azimuth)
        ren.GetActiveCamera().Elevation(ASDNeighActors.camera_elevation)
        ren.GetActiveCamera().Yaw(ASDNeighActors.camera_yaw)
        ren.GetActiveCamera().Roll(ASDNeighActors.camera_roll)
        ren.GetActiveCamera().Pitch(ASDNeighActors.camera_pitch)
        ren.GetActiveCamera().SetFocalPoint(ASDNeighActors.camera_focal)
        ren.GetActiveCamera().SetPosition(ASDNeighActors.camera_pos)
        ren.GetActiveCamera().SetViewUp(0, 1, 0)
        # -----------------------------------------------------------------------
        # Adding the actors for the neighbour mapping
        # -----------------------------------------------------------------------
        ren.AddActor(ASDNeighActors.NeighActor)
        ren.AddActor(ASDNeighActors.AtomsActor)
        ren.AddActor(ASDNeighActors.CenterActor)
        iren.Start()
        renWin.Render()
        return

    ##########################################################################
    # Updates the neighbours for rendering
    ##########################################################################

    def UpdateNeighbour(self, window, ASDdata, ASDGenActors, renWin, mode):
        """Function to update the visualization of the neighbours as the center atom is moved
        for both the Heisenberg exchange and the DM vectors. It also updated the label showing
        what types of atoms are contained in the neighbour cloud.
        Args:
            window: QMainWindow object where the visualizations are contained.
            ASData: class containing the data read from the UppASD simulations.
            ASDGenActors: class of general actors for the VTK visualization.
            renWin: VTK window where the visualization is performed.
            mode: value indicating whether one is visualizing
            the Heisenberg exchange or DMI vectors (1/2)

        Author
        ----------
        Jonathan Chico
        """
        # -----------------------------------------------------------------------
        # This makes sure that the neighbour is not updated twice
        # -----------------------------------------------------------------------
        slid_value = window.NeighSelectSlider.value()
        line_value = int(window.NeighSelectLineEdit.text())
        if slid_value != line_value:
            # -------------------------------------------------------------------
            # If one edits the line manually
            # -------------------------------------------------------------------
            if window.NeighSelectLineEdit.isModified():
                if mode == 1:
                    # -----------------------------------------------------------
                    # Read the data
                    # -----------------------------------------------------------
                    (
                        ASDdata.neighs,
                        ASDdata.atomCenter,
                        ASDdata.nTypes,
                        ASDdata.neigh_colors,
                        ASDdata.nNeighs,
                        ASDdata.num_types,
                        ASDdata.types_counters,
                        ASDdata.types,
                    ) = ASDdata.setNeighbours(
                        ASDdata.neighbours,
                        line_value - 1,
                        ASDdata.coord,
                        ASDdata.Neigh_strength,
                        ASDdata.curr_atom,
                        ASDdata.neigh_types,
                    )
                    # -----------------------------------------------------------
                    # Update the data
                    # -----------------------------------------------------------
                    ASDNeighActors.NeighGrid.SetPoints(ASDdata.neighs)
                    ASDNeighActors.NeighGrid.GetPointData().SetScalars(
                        ASDdata.neigh_colors
                    )
                    ASDNeighActors.NeighMapper.SetScalarRange(
                        ASDNeighActors.NeighGrid.GetScalarRange()
                    )
                    ASDNeighActors.NeighMapper.Update()
                    ASDNeighActors.CenterGrid.SetPoints(ASDdata.atomCenter)
                    ASDNeighActors.Center.Update()
                elif mode == 2:
                    # -----------------------------------------------------------
                    # Read the data
                    # -----------------------------------------------------------
                    # Calculate the neighbours to iAtom for the DM interactions
                    (
                        ASDdata.neighs,
                        ASDdata.atomCenter,
                        ASDdata.nTypes,
                        ASDdata.neigh_colors,
                        ASDdata.DM_vectors,
                        ASDdata.nNeighs,
                        ASDdata.num_types,
                        ASDdata.types_counters,
                        ASDdata.types,
                    ) = ASDdata.setDMNeighbours(
                        ASDdata.neighbours,
                        line_value - 1,
                        ASDdata.coord,
                        ASDdata.DM_vec,
                        ASDdata.DM_strength,
                        ASDdata.curr_atom,
                        ASDdata.neigh_types,
                    )
                    # -----------------------------------------------------------
                    # Update the data
                    # -----------------------------------------------------------
                    ASDNeighActors.NeighGrid.SetPoints(ASDdata.neighs)
                    ASDNeighActors.NeighGrid.GetPointData().SetScalars(
                        ASDdata.neigh_colors
                    )
                    ASDNeighActors.NeighGrid.GetPointData().SetVectors(
                        ASDdata.DM_vectors
                    )
                    ASDNeighActors.NeighMapper.SetScalarRange(
                        ASDNeighActors.NeighGrid.GetScalarRange()
                    )
                    ASDNeighActors.NeighMapper.Update()
                    ASDNeighActors.CenterGrid.SetPoints(ASDdata.atomCenter)
                    ASDNeighActors.Center.Update()
                # ---------------------------------------------------------------
                # Update the UI
                # ---------------------------------------------------------------
                window.NeighNumberLabel.setText(
                    f"Number of neighbours={ASDdata.nNeighs}"
                )
                window.NeighSelectSlider.setValue(line_value)
                window.NeighSelectLineEdit.setModified(False)
                for ii in range(0, ASDdata.num_types_total):
                    name = f"label_neigh_{ii}"
                    window.NeighTypesLabels[name].setText(
                        f"Num. Neighbours Type {ii + 1: 4d} = {0: 4d}"
                    )
                for ii in range(0, ASDdata.num_types):
                    name = f"label_neigh_{int(ASDdata.types[ii] - 1)}"
                    window.NeighTypesLabels[name].setText(
                        f"Num. Neighbours Type {ii + 1: 4d} = {ASDdata.types_counters[ii]: 4d}"
                    )
            # -------------------------------------------------------------------
            # If one updates the slider
            # -------------------------------------------------------------------
            else:
                if mode == 1:
                    # -----------------------------------------------------------
                    # Read the data
                    # -----------------------------------------------------------
                    (
                        ASDdata.neighs,
                        ASDdata.atomCenter,
                        ASDdata.nTypes,
                        ASDdata.neigh_colors,
                        ASDdata.nNeighs,
                        ASDdata.num_types,
                        ASDdata.types_counters,
                        ASDdata.types,
                    ) = ASDdata.setNeighbours(
                        ASDdata.neighbours,
                        slid_value - 1,
                        ASDdata.coord,
                        ASDdata.Neigh_strength,
                        ASDdata.curr_atom,
                        ASDdata.neigh_types,
                    )
                    # -----------------------------------------------------------
                    # Update the data
                    # -----------------------------------------------------------
                    ASDNeighActors.NeighGrid.SetPoints(ASDdata.neighs)
                    ASDNeighActors.NeighGrid.GetPointData().SetScalars(
                        ASDdata.neigh_colors
                    )
                    ASDNeighActors.NeighMapper.SetScalarRange(
                        ASDNeighActors.NeighGrid.GetScalarRange()
                    )
                    ASDNeighActors.NeighMapper.Update()
                    ASDNeighActors.CenterGrid.SetPoints(ASDdata.atomCenter)
                    ASDNeighActors.Center.Update()
                elif mode == 2:
                    # -----------------------------------------------------------
                    # Read the data
                    # -----------------------------------------------------------
                    # Calculate the neighbours to iAtom for the DM interactions
                    (
                        ASDdata.neighs,
                        ASDdata.atomCenter,
                        ASDdata.nTypes,
                        ASDdata.neigh_colors,
                        ASDdata.DM_vectors,
                        ASDdata.nNeighs,
                        ASDdata.num_types,
                        ASDdata.types_counters,
                        ASDdata.types,
                    ) = ASDdata.setDMNeighbours(
                        ASDdata.neighbours,
                        line_value - 1,
                        ASDdata.coord,
                        ASDdata.DM_vec,
                        ASDdata.DM_strength,
                        ASDdata.curr_atom,
                        ASDdata.neigh_types,
                    )
                    # -----------------------------------------------------------
                    # Update the data
                    # -----------------------------------------------------------
                    ASDNeighActors.NeighGrid.SetPoints(ASDdata.neighs)
                    ASDNeighActors.NeighGrid.GetPointData().SetScalars(
                        ASDdata.neigh_colors
                    )
                    ASDNeighActors.NeighGrid.GetPointData().SetVectors(
                        ASDdata.DM_vectors
                    )
                    ASDNeighActors.NeighMapper.SetScalarRange(
                        ASDNeighActors.NeighGrid.GetScalarRange()
                    )
                    ASDNeighActors.NeighMapper.Update()
                    ASDNeighActors.CenterGrid.SetPoints(ASDdata.atomCenter)
                    ASDNeighActors.Center.Update()
                # ---------------------------------------------------------------
                # Update the UI
                # ---------------------------------------------------------------
                window.NeighNumberLabel.setText(
                    f"Number of neighbours={ASDdata.nNeighs}"
                )
                window.NeighSelectLineEdit.setText(str(int(slid_value)))
                for ii in range(0, ASDdata.num_types_total):
                    name = f"label_neigh_{ii}"
                    window.NeighTypesLabels[name].setText(
                        f"Num. Neighbours Type {ii + 1: 4d} = {0: 4d}"
                    )
                for ii in range(0, ASDdata.num_types):
                    name = f"label_neigh_{int(ASDdata.types[ii] - 1)}"
                    window.NeighTypesLabels[name].setText(
                        f"Num. Neighbours Type {ii + 1: 4d} = {ASDdata.types_counters[ii]: 4d}"
                    )
        # -----------------------------------------------------------------------
        # Render the window
        # -----------------------------------------------------------------------
        renWin.Render()
        return
