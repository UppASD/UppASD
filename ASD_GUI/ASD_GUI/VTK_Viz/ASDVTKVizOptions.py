"""@package ASDVTKVizOptions
Contains a class with a set of functions dealing with the visualization options
in the VTK visualization mode. Dealing with visibility of objects, size of glyphs,
types of glyphs, colormap used etc.

Author
----------
Jonathan Chico
"""
# pylint: disable=invalid-name, no-name-in-module, no-member

import vtk
from PyQt6 import QtWidgets
from vtk import vtkPNGWriter, vtkPOVExporter, vtkWindowToImageFilter
from vtkmodules.vtkIOImage import vtkHDRReader
from vtkmodules.vtkRenderingCore import vtkTexture

from ASD_GUI.VTK_Viz import ASDVTKMomActors
from ASD_GUI.VTK_Viz import ASDVTKGenActors
from ASD_GUI.VTK_Viz import ASDVTKEneActors
from ASD_GUI.VTK_Viz import ASDVTKNeighActors
##########################################################################
# @brief Class containing the majority of the actions to update the visualizer
# @details Class containing the majority of the actions to update the visualizer.
# It handles the taking of snapshots, the toggling of projections, actor visibility
# and colormaps. It also controls many of the camera options so that they can be
# updated during runtime, allowing for finer control of the visualization.
# @author Jonathan Chico
##########################################################################


class ASDVizOptions:
    # lut = vtkLookupTable()
    def __init__(self):
        self.GenActors = ASDVTKGenActors.ASDGenActors()
        self.EneActors = ASDVTKEneActors.ASDEneActors()
        self.MomActors = ASDVTKMomActors.ASDMomActors()
        self.NeighActors = ASDVTKNeighActors.ASDNeighActors()
        
    ##########################################################################
    # @ brief A function that takes a renderwindow and saves its contents to a .png file
    # @author Anders Bergman
    ##########################################################################
    def Screenshot(self, renWin, number_of_screenshots, png_mode, pov_mode):
        """Function to take the rendering window and save it to file, either in .png
        or in .pov format.

        Args:
            renWin: current rendering window.
            number_of_screenshots: current number of the screenshot that is being saved.
            png_mode: (logical) variable indicated if the scene should be stored as a png
            pov_mode: (logical) variable indicated if the scene should be stored as a pov

        Author
        ----------
        Anders Bergman
        """

        win2im = vtkWindowToImageFilter()
        win2im.SetInput(renWin)
        win2im.Update()
        win2im.SetInputBufferTypeToRGBA()
        win2im.ReadFrontBufferOff()
        # -----------------------------------------------------------------------
        # Save snapshot as a '.pov'
        # -----------------------------------------------------------------------
        if pov_mode:
            povexp = vtkPOVExporter()
            povexp.SetInput(renWin)
            renWin.Render()
            povexp.SetFileName(f"snap{number_of_screenshots:05d}.pov")
            povexp.Write()
        # -----------------------------------------------------------------------
        # Save snapshot as a '.png'
        # -----------------------------------------------------------------------
        if png_mode:
            toPNG = vtkPNGWriter()
            toPNG.SetFileName(f"snap{number_of_screenshots:05d}.png")
            toPNG.SetInputConnection(win2im.GetOutputPort())
            toPNG.Write()
        return

    ##########################################################################
    # Update the dock window information when the camera is setup
    ##########################################################################
    def update_dock_info(self, current_Actors, Window):
        """
        Updates the dock information in the GUI with the current camera settings.
        """
        Window.FocalPosX.setText(str(current_Actors.camera_focal[0]))
        Window.FocalPosY.setText(str(current_Actors.camera_focal[1]))
        Window.FocalPosZ.setText(str(current_Actors.camera_focal[2]))
        Window.CamPosX.setText(str(current_Actors.camera_pos[0]))
        Window.CamPosY.setText(str(current_Actors.camera_pos[1]))
        Window.CamPosZ.setText(str(current_Actors.camera_pos[2]))
        Window.CamYawLineEdit.setText(str(current_Actors.camera_yaw))
        Window.CamRollLineEdit.setText(str(current_Actors.camera_roll))
        Window.CamPitchLineEdit.setText(str(current_Actors.camera_pitch))
        Window.CamAzimuthLineEdit.setText(str(current_Actors.camera_azimuth))
        Window.CamElevationLineEdit.setText(str(current_Actors.camera_elevation))
        return

    ##########################################################################
    # Toggle option for the axes
    ##########################################################################
    def toggle_Axes(self, check):
        """
        Toggles the visibility of the orientation marker based on the check value.
        """
        if check:
            self.GenActors.OrientMarker.SetEnabled(1)
        else:
            self.GenActors.OrientMarker.SetEnabled(0)
        return

    ##########################################################################
    # Toggle option for the scalar bar
    ##########################################################################
    def toggle_ScalarBar(self, check):
        """
        Toggles the visibility of the scalar bar widget based on the check parameter.
        """
        if check:
            self.GenActors.scalar_bar_widget.SetEnabled(1)
        else:
            self.GenActors.scalar_bar_widget.SetEnabled(0)
        return

    ##########################################################################
    # Toggle the visualization of the embedded cluster
    ##########################################################################
    def toggle_cluster(self, check):
        """
        Toggles the visibility of atom and atom_imp actors based on the check parameter.
        """
        if check:
            self.GenActors.atom.VisibilityOn()
            self.GenActors.atom_imp.VisibilityOn()
        else:
            self.GenActors.atom.VisibilityOff()
            self.GenActors.atom_imp.VisibilityOff()
        return

    ##########################################################################
    # Toggle the KMC particle visualization
    ##########################################################################
    def toggle_KMC(self, check):
        """
        Toggles the visibility of the KMC part actor based on the check value.
        """
        if check:
            self.MomActors.KMC_part_actor.VisibilityOn()
        else:
            self.MomActors.KMC_part_actor.VisibilityOff()
        return

    ##########################################################################
    # Toggle the plane clipper
    ##########################################################################
    def toggle_clipper(
        self, check, current_Actors, rdir, window, origin, vmin, vmax, renWin
    ):
        """
        Toggles the visibility of the clipper and current actors based on the check flag.
        """
        if check:
            self.GenActors.clipperActor.VisibilityOn()
            current_Actors.VisibilityOff()
            self.set_clipp_plane(rdir, window, origin, vmin, vmax, renWin)
        else:
            self.GenActors.clipperActor.VisibilityOff()
            current_Actors.VisibilityOn()
            renWin.Render()
        return

    ##########################################################################
    # Toggle the time label
    ##########################################################################
    def toggle_time_label(self, check):
        """
        Toggles the visibility of the time label widget based on the check value.
        """
        if check:
            self.GenActors.time_label_widget.On()
        else:
            self.GenActors.time_label_widget.Off()
        return

    ##########################################################################
    # Set the clipping plane such that the normal plane is 'origin'
    ##########################################################################
    def set_clipp_plane(self, rdir, window, origin, vmin, vmax, renWin):
        """
        Configures and sets the clipping plane for the visualization.
        """
        self.GenActors.plane.SetOrigin(origin)
        self.GenActors.plane.SetNormal(rdir)
        window.ClippingPlaneSlider.setMinimum(int(vmin))
        window.ClippingPlaneSlider.setMaximum(int(vmax))
        window.ClippingPlaneSlider.setValue(int(vmin))
        renWin.Render()
        return

    ##########################################################################
    # Set the position of the clipping plane via the slider
    ##########################################################################
    def ClippingUpdate(self, origin, window, renWin):
        """
        Updates the clipping plane position and refreshes the render window.
        """
        self.GenActors.plane.SetOrigin(origin)
        window.ClipPlaneLabel.setText(
            f"Clip. Plane Pos.={float(origin[0]):.1f},{float(origin[1]):.1f},{float(origin[2]):.1f}"
        )
        renWin.Render()
        return


    ##########################################################################
    # Toggle the atoms for the neighbour map
    ##########################################################################
    def toggle_NAtoms(self, check):
        """
        Toggles the visibility of the AtomsActor based on the check value.
        """
        if check:
            self.NeighActors.AtomsActor.VisibilityOn()
        else:
            self.NeighActors.AtomsActor.VisibilityOff()
        return

    ##########################################################################
    # Toggle the neighbour cloud for the neighbour map
    ##########################################################################
    def toggle_Neigh(self, check):
        """
        Toggles the visibility of the NeighActor based on the check value.
        """
        if check:
            self.NeighActors.NeighActor.VisibilityOn()
        else:
            self.NeighActors.NeighActor.VisibilityOff()
        return

    ##########################################################################
    # Set the opacity of the neighbour spheres
    ##########################################################################
    def NeighOpacityUpdate(self, value):
        """
        Update the opacity of the NeighActor based on the given value.
        """
        self.NeighActors.NeighActor.GetProperty().SetOpacity(value * 0.1)
        return

    ##########################################################################
    # Set the opacity of the atom spheres
    ##########################################################################
    def AtomOpacityUpdate(self, value):
        """
        Updates the opacity of the AtomsActor based on the given value.
        """
        self.NeighActors.AtomsActor.GetProperty().SetOpacity(value * 0.1)
        return

    ##########################################################################
    # Toggle SSAO on/off
    ##########################################################################
    def toggle_SSAO(self, check, ren):
        """
        Toggles the Screen Space Ambient Occlusion (SSAO) effect on the renderer.
        """

        if check:
            ren.UseSSAOOn()
            ren.SetSSAOKernelSize(512)
            ren.SetSSAORadius(3.0)
            ren.SetSSAOBias(0.1)
            ren.SSAOBlurOff()

            # self.toggle_HDRI(check=check,ren=ren)

        else:
            ren.UseSSAOOff()

        return

    ##########################################################################
    # Toggle automatic focal point determination on/off
    ##########################################################################
    def setFocalDisk(self, value, ren, renWin):
        """
        Adjusts the focal disk of the active camera and renders the window.
        """

        ren.GetActiveCamera().SetFocalDisk(value / 200.0)

        renWin.Render()

    ##########################################################################
    # Toggle automatic focal point determination on/off
    ##########################################################################
    def toggle_autoFocus(self, check, renWin):
        """
        Toggles the automatic focal distance for the depth of field pass and renders the window.
        """

        if check:
            self.dofPass.AutomaticFocalDistanceOn()
        else:
            self.dofPass.AutomaticFocalDistanceOff()

        renWin.Render()
        return

    ##########################################################################
    # Toggle depth of field focus on/off
    ##########################################################################
    def toggle_Focus(self, check, ren, renWin):
        """
        Toggles the Depth of Field (DOF) effect in the VTK renderer.
        """

        if check:
            # create the basic VTK render steps
            self.basicPasses = vtk.vtkRenderStepsPass()

            self.dofPass = vtk.vtkDepthOfFieldPass()
            self.dofPass.AutomaticFocalDistanceOff()
            self.dofPass.SetDelegatePass(self.basicPasses)

            # Tell the renderer to use our render pass pipeline.
            ren.GetActiveCamera().SetFocalDisk(0.5)
            ren.SetPass(self.dofPass)

            renWin.Render()

        else:
            print("DOF render pass can not be disabled.")
            ren.ReleaseGraphicsResources(renWin)

        return

    ##########################################################################
    # Toggle FXAA on/off
    ##########################################################################
    def toggle_FXAA(self, check, ren, renWin):
        """
        Toggles FXAA (Fast Approximate Anti-Aliasing) on or off for the given renderer and render window.
        """
        if check:
            ren.UseFXAAOn()
            ren.GetFXAAOptions().SetUseHighQualityEndpoints(True)
            renWin.SetMultiSamples(4)
        else:
            ren.UseFXAAOff()
        return

    ##########################################################################
    # Update glyph resolution
    ##########################################################################
    def GlyphQualityUpdate(self, window, value, viz_type, mode, renWin):
        """
        Updates the glyph quality for different visualization types and renders the window.
        """
        if viz_type == "M":
            try:
                window.MomActors.spinarrow.SetTipResolution(value)
                window.MomActors.spinarrow.SetShaftResolution(value)
            except AttributeError:
                pass
            try:
                window.MomActors.spinsphere.SetThetaResolution(value)
                window.MomActors.spinsphere.SetPhiResolution(value)
            except AttributeError:
                pass
            try:
                window.MomActors.spincones.SetResolution(value)
            except AttributeError:
                pass

        if viz_type == "N":
            if mode == 1:
                window.NeighActors.NeighGlyphs.SetThetaResolution(value)
                window.NeighActors.NeighGlyphs.SetPhiResolution(value)
            if mode == 2:
                window.NeighActors.NeighGlyphs.SetTipResolution(value)
                window.NeighActors.NeighGlyphs.SetShaftResolution(value)
        if viz_type == "E":
            window.EneActors.EneAtom.SetThetaResolution(value)
            window.EneActors.EneAtom.SetPhiResolution(value)

        renWin.Render()
        return
