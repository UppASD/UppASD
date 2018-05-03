#!/usr/bin/env vtkpython
################################################################################
# CLASS: ASDNeighActors
# @author Jonathan Chico (08/09/2017)
# @description
# Wrapper class to add the VTK  actors for the visualization of UppASD data in
# neighbour visualization mode.
################################################################################
import vtk
import ASDVTKReading

class ASDNeighActors():

    ###########################################################################
    # This defines the actors for the visualization of neighbours from the struct file
    ############################################################################
    def AddASDNeigh_actors(self,ren,renWin,mode,viz_type,iren):
        # Reading the data for the nieghbours
        ASD_data=ASDVTKReading.ASDReading()
        ASD_data.ReadingWrapper(mode,viz_type)

        # Set the maximum number of entried in the slider to be equal to the number of atoms
        ########################################################################
        # Data structures for the atoms in the neighbour map
        ########################################################################
        AtomGrid=vtk.vtkPolyData()
        AtomGrid.SetPoints(ASD_data.coord)
        ASDNeighActors.SLMax=AtomGrid.GetNumberOfPoints()
        # Atom sphere
        Atom = vtk.vtkSphereSource()
        Atom.SetRadius(0.50)
        Atom.SetThetaResolution(10)
        Atom.SetPhiResolution(10)
        # Atom glyph
        Atoms = vtk.vtkGlyph3DMapper()
        Atoms.SetInputData(AtomGrid)
        Atoms.SetSourceConnection(Atom.GetOutputPort())
        Atoms.SetScaleFactor(0.5)
        Atoms.ClampingOn()
        Atoms.SetScaleModeToNoDataScaling()
        Atoms.SetColorModeToMapScalars()
        Atoms.Update()
        # Atoms actors
        ASDNeighActors.AtomsActor = vtk.vtkLODActor()
        ASDNeighActors.AtomsActor.SetMapper(Atoms)
        ASDNeighActors.AtomsActor.GetProperty().SetOpacity(0.9)
        ASDNeighActors.AtomsActor.GetProperty().SetColor(0.5, 0.5, 0.5)
        ########################################################################
        # Data structures for the neighbours in the neighbour mapper
        ########################################################################
        #Grid for neighbours
        ASDNeighActors.NeighGrid=vtk.vtkPolyData()
        ASDNeighActors.NeighGrid.SetPoints(ASD_data.neighs)
        ASDNeighActors.NeighGrid.GetPointData().SetScalars(ASD_data.nTypes)
        ASDNeighActors.NumNeigh=ASDNeighActors.NeighGrid.GetNumberOfPoints()
        #Neighbour glyphs
        ASDNeighActors.Neigh = vtk.vtkSphereSource()
        ASDNeighActors.Neigh.SetRadius(0.50)
        ASDNeighActors.Neigh.SetThetaResolution(20)
        ASDNeighActors.Neigh.SetPhiResolution(20)
        # Set up Neighbour glyphs
        ASDNeighActors.Neighs = vtk.vtkGlyph3DMapper()
        ASDNeighActors.Neighs.SetInputData(ASDNeighActors.NeighGrid)
        ASDNeighActors.Neighs.SetSourceConnection(ASDNeighActors.Neigh.GetOutputPort())
        ASDNeighActors.Neighs.SetScaleFactor(1.25)
        ASDNeighActors.Neighs.SetColorModeToMapScalars()
        ASDNeighActors.Neighs.SetScaleModeToNoDataScaling()
        ASDNeighActors.Neighs.SetScalarRange(ASDNeighActors.NeighGrid.GetPointData().GetScalars().GetRange())
        # Color lookup table for neighbours
        lut = vtk.vtkLookupTable()
        lut.SetNumberOfColors(2)
        lut.SetTableValue(0,0.0,1.0,0.0,0.5)
        lut.SetTableValue(1,0.0,0.0,1.0,0.5)
        lut.SetTableRange(ASD_data.nTypes.GetRange())
        lut.Build()
        # Fill up Lookup Table with appropiate colors
        ASDNeighActors.Neighs.SetLookupTable(lut)
        ASDNeighActors.Neighs.Update()
        # Neighbour actors
        ASDNeighActors.NeighActor = vtk.vtkLODActor()
        ASDNeighActors.NeighActor.SetMapper(ASDNeighActors.Neighs)
        ASDNeighActors.NeighActor.GetProperty().SetSpecularColor(0.4, 0.8, 0.4)
        ASDNeighActors.NeighActor.GetProperty().SetSpecular(0.3)
        ASDNeighActors.NeighActor.GetProperty().SetSpecularPower(60)
        ASDNeighActors.NeighActor.GetProperty().SetAmbient(0.2)
        ASDNeighActors.NeighActor.GetProperty().SetDiffuse(0.8)
        # Setting up information needed for the camera
        self.xmin,self.xmax = self.AtomsActor.GetXRange()
        self.ymin,self.ymax = self.AtomsActor.GetYRange()
        self.zmin,self.zmax = self.AtomsActor.GetZRange()
        self.xmid = (self.xmin+self.xmax)/2
        self.ymid = (self.ymin+self.ymax)/2
        self.zmid = (self.zmin+self.zmax)/2
        self.height=max(self.xmax,self.ymax,self.zmax)*1.75
        # Defining the camera directions
        ren.GetActiveCamera().Azimuth(0)
        ren.GetActiveCamera().Elevation(0)
        ren.GetActiveCamera().SetFocalPoint(self.xmid,self.ymid,self.zmid)
        ren.GetActiveCamera().SetPosition(self.xmid,self.ymid,self.height)
        # Adding the actors for the neighbour mapping
        ren.AddActor(ASDNeighActors.NeighActor)
        ren.AddActor(ASDNeighActors.AtomsActor)
        iren.Start()
        renWin.Render()
