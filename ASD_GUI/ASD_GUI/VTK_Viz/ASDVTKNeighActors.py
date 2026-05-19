"""@package ASDVTKNeighActors
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
    """
    Class for visualizing neighbors in a VTK environment.
    """
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
        self.timer_count = 0
        self.camera_pos = np.zeros(3, dtype=np.float32)
        self.camera_focal = np.zeros(3, dtype=np.float32)
        self.camera_yaw = 0.0
        self.camera_roll = 0.0
        self.camera_pitch = 0.0
        self.camera_azimuth = 0.0
        self.camera_elevation = 0.0
        # -----------------------------------------------------------------------
        # Data structures for the atoms in the neighbour map
        # -----------------------------------------------------------------------
        AtomGrid = vtk.vtkPolyData()
        AtomGrid.SetPoints(ASDdata.coord)
        self.SLMax = AtomGrid.GetNumberOfPoints()
        # -----------------------------------------------------------------------
        # Atom sphere
        # -----------------------------------------------------------------------
        Atom = vtk.vtkSphereSource()
        Atom.SetRadius(0.50)
        Atom.SetThetaResolution(16)
        Atom.SetPhiResolution(16)
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
        self.AtomsActor = vtk.vtkLODActor()
        self.AtomsActor.SetMapper(Atoms)
        self.AtomsActor.GetProperty().SetOpacity(0.9)
        self.AtomsActor.GetProperty().SetColor(0.5, 0.5, 0.5)
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
        self.NeighGrid = vtk.vtkPolyData()
        self.NeighGrid.SetPoints(ASDdata.neighs)
        self.NeighGrid.GetPointData().SetScalars(ASDdata.neigh_colors)
        if mode == 2:
            self.NeighGrid.GetPointData().SetVectors(ASDdata.DM_vectors)
        self.NumNeigh = self.NeighGrid.GetNumberOfPoints()
        scalar_range = self.NeighGrid.GetScalarRange()
        # -----------------------------------------------------------------------
        # Finding useful geometrical information of the sample
        # Finding the middle of the sample
        # Also making sure that if the sample is 2D one has no problem with boudings
        # this is mostly useful if splatters are used
        # -----------------------------------------------------------------------
        (
            self.xmin,
            self.xmax,
            self.ymin,
            self.ymax,
            self.zmin,
            self.zmax,
        ) = self.NeighGrid.GetBounds()
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
        # Neighbour glyphs
        # -----------------------------------------------------------------------
        if mode == 1:
            self.NeighGlyphs = vtk.vtkSphereSource()
            self.NeighGlyphs.SetRadius(0.50)
            self.NeighGlyphs.SetThetaResolution(16)
            self.NeighGlyphs.SetPhiResolution(16)
        else:
            self.NeighGlyphs = vtk.vtkArrowSource()
            self.NeighGlyphs.SetTipRadius(0.20)
            self.NeighGlyphs.SetShaftRadius(0.10)
            self.NeighGlyphs.SetTipResolution(16)
            self.NeighGlyphs.SetShaftResolution(16)
        # -----------------------------------------------------------------------
        # Glyph source
        # -----------------------------------------------------------------------
        self.NeighGlyph3D = vtk.vtkGlyph3D()
        self.NeighGlyph3D.SetSourceConnection(self.NeighGlyphs.GetOutputPort())
        self.NeighGlyph3D.SetVectorModeToUseNormal()
        self.NeighGlyph3D.SetInputData(self.NeighGrid)
        if mode == 1:
            self.NeighGlyph3D.SetScaleFactor(1.05)
        elif mode == 2:
            self.NeighGlyph3D.SetScaleFactor(1.25)
        self.NeighGlyph3D.SetColorModeToColorByScalar()
        self.NeighGlyph3D.SetScaleModeToDataScalingOff()
        if mode == 2:
            self.NeighGlyph3D.SetVectorModeToUseVector()
        self.NeighGlyph3D.Update()
        # -----------------------------------------------------------------------
        # Set up Neighbour glyphs
        # -----------------------------------------------------------------------
        self.NeighMapper = vtk.vtkPolyDataMapper()
        self.NeighMapper.SetInputConnection(self.NeighGlyph3D.GetOutputPort())
        self.NeighMapper.SetScalarRange(scalar_range)
        self.NeighMapper.SetLookupTable(self.lut)
        self.NeighMapper.SetColorModeToMapScalars()
        self.NeighMapper.Update()
        # -----------------------------------------------------------------------
        # Neighbour actors
        # -----------------------------------------------------------------------
        self.NeighActor = vtk.vtkLODActor()
        self.NeighActor.SetMapper(self.NeighMapper)
        self.NeighActor.GetProperty().SetSpecular(0.3)
        self.NeighActor.GetProperty().SetSpecularPower(60)
        self.NeighActor.GetProperty().SetAmbient(0.2)
        self.NeighActor.GetProperty().SetDiffuse(0.8)
        # -----------------------------------------------------------------------
        # Grid for the center atom
        # -----------------------------------------------------------------------
        self.CenterGrid = vtk.vtkPolyData()
        self.CenterGrid.SetPoints(ASDdata.atomCenter)
        self.CenterGrid.Modified()
        # -----------------------------------------------------------------------
        # Source for the center atom, a sphere
        # -----------------------------------------------------------------------
        self.CenterSource = vtk.vtkSphereSource()
        self.CenterSource.SetRadius(0.50)
        self.CenterSource.SetThetaResolution(16)
        self.CenterSource.SetPhiResolution(16)
        # -----------------------------------------------------------------------
        # Mapper for the center actor
        # -----------------------------------------------------------------------
        self.Center = vtk.vtkGlyph3DMapper()
        self.Center.SetInputData(self.CenterGrid)
        self.Center.SetSourceConnection(self.CenterSource.GetOutputPort())
        self.Center.SetScaleFactor(1.25)
        self.Center.SetScaleModeToNoDataScaling()
        self.Center.Update()
        # -----------------------------------------------------------------------
        # Actor for the center atom
        # -----------------------------------------------------------------------
        self.CenterActor = vtk.vtkLODActor()
        self.CenterActor.SetMapper(self.Center)
        self.CenterActor.GetProperty().SetColor(0.4, 0.8, 0.4)
        self.CenterActor.GetProperty().SetSpecular(0.3)
        self.CenterActor.GetProperty().SetSpecularPower(60)
        self.CenterActor.GetProperty().SetAmbient(0.2)
        self.CenterActor.GetProperty().SetDiffuse(0.8)
        # -----------------------------------------------------------------------
        # Defining the camera directions
        # -----------------------------------------------------------------------
        ren.GetActiveCamera().Azimuth(0)
        ren.GetActiveCamera().Elevation(0)
        ren.GetActiveCamera().SetFocalPoint(self.xmid, self.ymid, self.zmid)
        ren.GetActiveCamera().SetPosition(self.xmid, self.ymid, self.height)
        ren.GetActiveCamera().Azimuth(self.camera_azimuth)
        ren.GetActiveCamera().Elevation(self.camera_elevation)
        ren.GetActiveCamera().Yaw(self.camera_yaw)
        ren.GetActiveCamera().Roll(self.camera_roll)
        ren.GetActiveCamera().Pitch(self.camera_pitch)
        ren.GetActiveCamera().SetFocalPoint(self.camera_focal)
        ren.GetActiveCamera().SetPosition(self.camera_pos)
        ren.GetActiveCamera().SetViewUp(0, 1, 0)
        # -----------------------------------------------------------------------
        # Adding the actors for the neighbour mapping
        # -----------------------------------------------------------------------
        ren.AddActor(self.NeighActor)
        ren.AddActor(self.AtomsActor)
        ren.AddActor(self.CenterActor)
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
                    self.NeighGrid.SetPoints(ASDdata.neighs)
                    self.NeighGrid.GetPointData().SetScalars(ASDdata.neigh_colors)
                    self.NeighMapper.SetScalarRange(self.NeighGrid.GetScalarRange())
                    self.NeighMapper.Update()
                    self.CenterGrid.SetPoints(ASDdata.atomCenter)
                    self.Center.Update()
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
                    self.NeighGrid.SetPoints(ASDdata.neighs)
                    self.NeighGrid.GetPointData().SetScalars(ASDdata.neigh_colors)
                    self.NeighGrid.GetPointData().SetVectors(ASDdata.DM_vectors)
                    self.NeighMapper.SetScalarRange(self.NeighGrid.GetScalarRange())
                    self.NeighMapper.Update()
                    self.CenterGrid.SetPoints(ASDdata.atomCenter)
                    self.Center.Update()
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
                    self.NeighGrid.SetPoints(ASDdata.neighs)
                    self.NeighGrid.GetPointData().SetScalars(ASDdata.neigh_colors)
                    self.NeighMapper.SetScalarRange(self.NeighGrid.GetScalarRange())
                    self.NeighMapper.Update()
                    self.CenterGrid.SetPoints(ASDdata.atomCenter)
                    self.Center.Update()
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
                    self.NeighGrid.SetPoints(ASDdata.neighs)
                    self.NeighGrid.GetPointData().SetScalars(ASDdata.neigh_colors)
                    self.NeighGrid.GetPointData().SetVectors(ASDdata.DM_vectors)
                    self.NeighMapper.SetScalarRange(self.NeighGrid.GetScalarRange())
                    self.NeighMapper.Update()
                    self.CenterGrid.SetPoints(ASDdata.atomCenter)
                    self.Center.Update()
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
