##########################################################################
# CLASS: ASDGenActors
# @author Jonathan Chico (15/05/2018)
# @description
# Wrapper class to add the VTK  actors for the visualization of auxiliary visualization data.
# It contains the needed data to add the actors, modify them, as well as some helper
# functions to change them.
##########################################################################
# pylint: disable=invalid-name, no-name-in-module, no-member
import vtk


class ASDGenActors:
    """Class that defines general actors which can be used in the different visualization types.
    It includes axes to indicate the orientation, scalar bar to indicate the range of the variables
    and cluster gaussian splatters to indicate the presence of clusters."""

    ##########################################################################
    # Tries to eliminate the general actors
    ##########################################################################
    def reset_GenActors(window):
        """
        Reset the visualization actors by deleting their references.
        """
        del ASDGenActors.scalar_bar
        del ASDGenActors.scalar_bar_widget
        del ASDGenActors.scalarBarRep
        del ASDGenActors.axes
        del ASDGenActors.OrientMarker
        return

    ##########################################################################
    # Set general actors such as the scalar bar and the axes widget
    ##########################################################################
    def Add_GenActors(
        self, iren, renWin, method, lut, ren, window, current_Actors, flag2D
    ):
        """
        Add general actors to the VTK renderer.

        This method sets up various VTK actors and widgets, including text labels,
        clippers, scalar bars, and axes, and adds them to the renderer.

        Parameters:
        - iren: vtkRenderWindowInteractor
            The interactor for the render window.
        - renWin: vtkRenderWindow
            The render window.
        - method: vtkAlgorithm
            The VTK algorithm used for generating the data.
        - lut: vtkLookupTable
            The lookup table for scalar mapping.
        - ren: vtkRenderer
            The renderer to which actors are added.
        - window: vtkRenderWindow
            The render window.
        - current_Actors: object
            An object containing current actor properties.
        - flag2D: bool
            A flag indicating whether the data is 2D or 3D.

        Returns:
        None

        Author:
        Jonathan Chico
        """

        # Create the TextActor
        ASDGenActors.time_label = vtk.vtkTextActor()
        ASDGenActors.time_label.SetInput(f"{0.00: 4.2f} ns")
        ASDGenActors.time_label.GetTextProperty().SetColor((0, 0, 0))
        # Create the text representation. Used for positioning the text_actor
        ASDGenActors.time_label_rep = vtk.vtkTextRepresentation()
        ASDGenActors.time_label_rep.GetPositionCoordinate().SetValue(0.80, 0.90)
        ASDGenActors.time_label_rep.GetPosition2Coordinate().SetValue(0.10, 0.10)
        #######################################################################
        # Creating the actual widget
        #######################################################################
        ASDGenActors.time_label_widget = vtk.vtkTextWidget()
        ASDGenActors.time_label_widget.SetRepresentation(ASDGenActors.time_label_rep)
        ASDGenActors.time_label_widget.SetInteractor(iren)
        ASDGenActors.time_label_widget.SetTextActor(ASDGenActors.time_label)
        ASDGenActors.time_label_widget.SelectableOff()
        ASDGenActors.time_label_widget.Off()
        #######################################################################
        # Creation of the data structures for the data clipping
        #######################################################################
        # Right now this only can clip polydata, which is fine for 2D structures
        # however, for the 3d Delaunay tessellation, the output is an unstructured
        # grid, which means that annother type of clipper is required
        #######################################################################
        ASDGenActors.plane = vtk.vtkPlane()
        ASDGenActors.plane.SetOrigin(
            current_Actors.xmin, current_Actors.ymid, current_Actors.zmid
        )
        ASDGenActors.plane.SetNormal(1, 0, 0)
        #######################################################################
        # Check which kind of clipper must be used, as 2D and 3D data must be
        # treated differently
        #######################################################################
        if flag2D:
            ASDGenActors.clipper = vtk.vtkClipPolyData()
        else:
            ASDGenActors.clipper = vtk.vtkClipVolume()
        # Set the common variables for the clipper mapper
        ASDGenActors.clipper.SetInputConnection(method.GetOutputPort())
        ASDGenActors.clipper.SetClipFunction(ASDGenActors.plane)
        ASDGenActors.clipper.InsideOutOn()
        # Mapper of the clipper
        ASDGenActors.clipperMapper = vtk.vtkDataSetMapper()
        ASDGenActors.clipperMapper.SetScalarRange(method.GetInput().GetScalarRange())
        ASDGenActors.clipperMapper.SetInputConnection(
            ASDGenActors.clipper.GetOutputPort()
        )
        ASDGenActors.clipperMapper.SetLookupTable(lut)
        # Creating the actor
        ASDGenActors.clipperActor = vtk.vtkLODActor()
        ASDGenActors.clipperActor.SetMapper(ASDGenActors.clipperMapper)
        ASDGenActors.clipperActor.VisibilityOff()
        # Adding the actor to the scene
        ren.AddActor(ASDGenActors.clipperActor)
        #######################################################################
        # Setting the information for the scalar bar widget
        #######################################################################
        # Create the scalar bar actor
        ASDGenActors.scalar_bar = vtk.vtkScalarBarActor()
        ASDGenActors.scalar_bar.SetLookupTable(lut)
        ASDGenActors.scalar_bar.GetLabelTextProperty().SetColor(0.0, 0.0, 0.0)
        ASDGenActors.scalar_bar.SetNumberOfLabels(5)
        ASDGenActors.scalar_bar.GetLabelTextProperty().ShadowOff()
        ASDGenActors.scalar_bar.GetLabelTextProperty().BoldOn()
        ASDGenActors.scalar_bar.GetLabelTextProperty().ItalicOff()
        ASDGenActors.scalar_bar.GetLabelTextProperty().SetFontSize(8)
        ASDGenActors.scalar_bar.SetLabelFormat("%-#6.1E")
        ASDGenActors.scalar_bar.SetBarRatio(0.5)
        ASDGenActors.scalar_bar.DrawBackgroundOn()
        ASDGenActors.scalar_bar.DrawTickLabelsOn()
        # Create the scalar bar widget
        ASDGenActors.scalar_bar_widget = vtk.vtkScalarBarWidget()
        ASDGenActors.scalar_bar_widget.SetScalarBarActor(ASDGenActors.scalar_bar)
        # Representation to actually control where the scalar bar is
        ASDGenActors.scalarBarRep = ASDGenActors.scalar_bar_widget.GetRepresentation()
        ASDGenActors.scalarBarRep.SetOrientation(0)  # 0 = Horizontal, 1 = Vertical
        ASDGenActors.scalarBarRep.GetPositionCoordinate().SetValue(0.30, 0.05)
        ASDGenActors.scalarBarRep.GetPosition2Coordinate().SetValue(0.50, 0.05)
        ASDGenActors.scalar_bar_widget.SetInteractor(iren)
        ASDGenActors.scalar_bar_widget.On()
        #######################################################################
        # Setting the information for the axes widget
        #######################################################################
        # Create the axes actor
        try:
            ASDGenActors.axes
        except BaseException:
            ASDGenActors.axes = vtk.vtkAxesActor()
            ASDGenActors.axes.SetShaftTypeToCylinder()
            ASDGenActors.axes.SetCylinderRadius(0.05)
            ASDGenActors.axes.SetNormalizedShaftLength(0.85, 0.85, 0.85)
            ASDGenActors.axes.SetNormalizedTipLength(0.40, 0.40, 0.40)
            ASDGenActors.axes.SetConeResolution(40)
            ASDGenActors.axes.SetCylinderResolution(40)
            # The properties of the text can be controlled independently
            ASDGenActors.axes.GetXAxisCaptionActor2D().GetCaptionTextProperty().SetColor(
                0.0, 0.0, 0.0
            )
            ASDGenActors.axes.GetYAxisCaptionActor2D().GetCaptionTextProperty().SetColor(
                0.0, 0.0, 0.0
            )
            ASDGenActors.axes.GetZAxisCaptionActor2D().GetCaptionTextProperty().SetColor(
                0.0, 0.0, 0.0
            )
            ASDGenActors.axes.GetXAxisCaptionActor2D().GetCaptionTextProperty().ShadowOff()
            ASDGenActors.axes.GetYAxisCaptionActor2D().GetCaptionTextProperty().ShadowOff()
            ASDGenActors.axes.GetZAxisCaptionActor2D().GetCaptionTextProperty().ShadowOff()
        else:
            pass
        #######################################################################
        # The axes actor is then used as an orientation marker widget, the advantage
        # of setting it up as a widget is that it is interactive and one can move it
        # and that it moves as the zoom changes
        # Must make sure that the widget is part of the main class so that it can
        # be actually rendered and no segfaults occur
        #######################################################################
        try:
            ASDGenActors.OrientMarker
        except BaseException:
            ASDGenActors.OrientMarker = vtk.vtkOrientationMarkerWidget()
            ASDGenActors.OrientMarker.SetOutlineColor(0.9300, 0.5700, 0.1300)
            ASDGenActors.OrientMarker.SetOrientationMarker(ASDGenActors.axes)
            ASDGenActors.OrientMarker.SetViewport(0.0, 0.0, 0.2, 0.2)
            ASDGenActors.OrientMarker.SetInteractor(iren)
            ASDGenActors.OrientMarker.EnabledOn()  # <== application freeze-crash
            ASDGenActors.OrientMarker.InteractiveOn()
        else:
            pass
        iren.Start()
        renWin.Render()
        return

    ##########################################################################
    # Setting up data for adding the impurity cluster
    ##########################################################################
    def Add_ClusterActors(self, ASDdata, iren, renWin, ren):

        #######################################################################
        # Data structures for the impurity cluster
        #######################################################################
        # Passing the data from the cluster to the PolyData
        src_clus = vtk.vtkPolyData()
        src_clus.SetPoints(ASDdata.coord_c)
        src_clus.GetPointData().SetScalars(ASDdata.colors_clus)
        # Passing the data from the selected impurities
        src_imp = vtk.vtkPolyData()
        src_imp.SetPoints(ASDdata.points_clus_imp)
        src_imp.GetPointData().SetScalars(ASDdata.colors_imp)
        src_imp.Modified()
        # Setting up the gaussian splatter for the clusters
        atomSource = vtk.vtkGaussianSplatter()
        atomSource.SetInputData(src_clus)
        atomSource.SetRadius(0.100)
        atomSource.ScalarWarpingOff()
        atomSource.SetExponentFactor(-20)
        atomSource.Update()
        bound = atomSource.GetModelBounds()
        atomSource.SetModelBounds(
            bound[0], bound[1], bound[2], bound[3], bound[4] * 0.25, bound[5] * 0.25
        )
        atomSource.Update()
        # Setting up a contour filter
        atomSurface = vtk.vtkContourFilter()
        atomSurface.SetInputConnection(atomSource.GetOutputPort())
        atomSurface.SetValue(0, 0.01)
        # Setting up the mapper
        atomMapper = vtk.vtkPolyDataMapper()
        atomMapper.SetInputConnection(atomSurface.GetOutputPort())
        atomMapper.ScalarVisibilityOff()
        # Creating the actor for the smooth surfaces
        ASDGenActors.atom = vtk.vtkActor()
        ASDGenActors.atom.SetMapper(atomMapper)
        ASDGenActors.atom.GetProperty().SetColor(0.0, 0.0, 0.0)
        ASDGenActors.atom.GetProperty().EdgeVisibilityOff()
        ASDGenActors.atom.GetProperty().SetSpecularPower(30)
        ASDGenActors.atom.GetProperty().SetAmbient(0.2)
        ASDGenActors.atom.GetProperty().SetDiffuse(0.8)
        ASDGenActors.atom.GetProperty().SetOpacity(0.25)
        # Set up imp sources
        atomSource_imp = vtk.vtkSphereSource()
        atomSource_imp.SetRadius(2.5)
        atomSource_imp.SetThetaResolution(20)
        atomSource_imp.SetPhiResolution(20)
        # Mapping the spheres to the actual points on the selected impurities
        atomMapper_imp = vtk.vtkGlyph3DMapper()
        atomMapper_imp.SetInputData(src_imp)
        atomMapper_imp.SetSourceConnection(atomSource_imp.GetOutputPort())
        atomMapper_imp.SetScaleFactor(0.2)
        atomMapper_imp.SetScaleModeToNoDataScaling()
        atomMapper_imp.Update()
        # Creating the selected impurity actors
        ASDGenActors.atom_imp = vtk.vtkLODActor()
        ASDGenActors.atom_imp.SetMapper(atomMapper_imp)
        ASDGenActors.atom_imp.GetProperty().SetSpecular(0.3)
        ASDGenActors.atom_imp.GetProperty().SetSpecularPower(30)
        ASDGenActors.atom_imp.GetProperty().SetAmbient(0.2)
        ASDGenActors.atom_imp.GetProperty().SetDiffuse(0.8)
        # If there is information about the cluster add the needed actors
        ren.AddActor(ASDGenActors.atom)
        ren.AddActor(ASDGenActors.atom_imp)
        #######################################################################
        # Start the renderer
        #######################################################################
        iren.Start()
        renWin.Render()
        return

    ##########################################################################
    # Update the clipper actor for each of the possible visualization types
    ##########################################################################

    def UpdateClipper(self, window, current_Actors, ASDVizOpt, renWin, viz_type):
        """
        Update the clipping settings based on the sender widget.

        This method updates the clipping plane or toggles the clipper based on the 
        sender widget. It supports clipping in the x, y, and z directions and updates 
        the clipping plane location based on the slider value.

        Parameters:
        window (QWidget): The window containing the sender widget.
        current_Actors (object): The current actors in the visualization.
        ASDVizOpt (object): The visualization options handler.
        renWin (object): The render window.
        viz_type (str): The type of visualization ("M", "N", or "E").

        Author:
        Jonathan Chico
        """
        if window.sender() == window.ClippBox:
            rdir = (1, 0, 0)
            origin = (current_Actors.xmin, current_Actors.ymid, current_Actors.zmid)
            vmin = current_Actors.xmin
            vmax = current_Actors.xmax
            if viz_type == "M":
                ASDVizOpt.toggle_clipper(
                    check=window.ClippBox.isChecked(),
                    current_Actors=current_Actors.MagDensActor,
                    rdir=rdir,
                    window=window,
                    origin=origin,
                    vmin=vmin,
                    vmax=vmax,
                    renWin=renWin,
                )
            if viz_type == "N":
                ASDVizOpt.toggle_clipper(
                    check=window.ClippBox.isChecked(),
                    current_Actors=current_Actors.NeighActor,
                    rdir=rdir,
                    window=window,
                    origin=origin,
                    vmin=vmin,
                    vmax=vmax,
                    renWin=renWin,
                )
            if viz_type == "E":
                ASDVizOpt.toggle_clipper(
                    check=window.ClippBox.isChecked(),
                    current_Actors=current_Actors.EneDensActor,
                    rdir=rdir,
                    window=window,
                    origin=origin,
                    vmin=vmin,
                    vmax=vmax,
                    renWin=renWin,
                )
        #######################################################################
        # set the clipping plane to be in the x-direction
        #######################################################################
        if window.sender() == window.ClippPlaneXCheck:
            rdir = (1, 0, 0)
            origin = (current_Actors.xmin, current_Actors.ymid, current_Actors.zmid)
            vmin = current_Actors.xmin
            vmax = current_Actors.xmax
            ASDVizOpt.set_clipp_plane(
                rdir=rdir, window=window, origin=origin, vmin=vmin, vmax=vmax, renWin=renWin
            )
        #######################################################################
        # set the clipping plane to be in the y-direction
        #######################################################################
        if window.sender() == window.ClippPlaneYCheck:
            rdir = (0, 1, 0)
            origin = (current_Actors.xmid, current_Actors.ymin, current_Actors.zmid)
            vmin = current_Actors.ymin
            vmax = current_Actors.ymax
            ASDVizOpt.set_clipp_plane(
                rdir=rdir, window=window, origin=origin, vmin=vmin, vmax=vmax, renWin=renWin
            )
        #######################################################################
        # set the clipping plane to be in the z-direction
        #######################################################################
        if window.sender() == window.ClippPlaneZCheck:
            rdir = (0, 0, 1)
            origin = (current_Actors.xmid, current_Actors.ymid, current_Actors.zmin)
            vmin = current_Actors.zmin
            vmax = current_Actors.zmax
            ASDVizOpt.set_clipp_plane(
                rdir=rdir, window=window, origin=origin, vmin=vmin, vmax=vmax, renWin=renWin
            )
        #######################################################################
        # Update the clipping plane location
        #######################################################################
        if window.sender() == window.ClippingPlaneSlider:
            if window.ClippPlaneXCheck.isChecked():
                origin = (
                    window.ClippingPlaneSlider.value(),
                    current_Actors.ymid,
                    current_Actors.zmid,
                )
                ASDVizOpt.ClippingUpdate(origin=origin, window=window, renWin=renWin)
            if window.ClippPlaneYCheck.isChecked():
                origin = (
                    current_Actors.xmid,
                    window.ClippingPlaneSlider.value(),
                    current_Actors.zmid,
                )
                ASDVizOpt.ClippingUpdate(origin=origin, window=window, renWin=renWin)
            if window.ClippPlaneZCheck.isChecked():
                origin = (
                    current_Actors.xmid,
                    current_Actors.ymid,
                    window.ClippingPlaneSlider.value(),
                )
                ASDVizOpt.ClippingUpdate(origin=origin, window=window, renWin=renWin)
        return
