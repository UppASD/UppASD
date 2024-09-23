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
from vtkmodules.vtkCommonCore import vtkLookupTable
from vtkmodules.vtkIOImage import vtkHDRReader
from vtkmodules.vtkRenderingCore import vtkTexture

from ASD_GUI.VTK_Viz import (ASDVTKEneActors, ASDVTKGenActors, ASDVTKMomActors,
                             ASDVTKNeighActors)

##########################################################################
# @brief Class containing the majority of the actions to update the visualizer
# @details Class containing the majority of the actions to update the visualizer.
# It handles the taking of snapshots, the toggling of projections, actor visibility
# and colormaps. It also controls many of the camera options so that they can be
# updated during runtime, allowing for finer control of the visualization.
# @author Jonathan Chico
##########################################################################


class ASDVizOptions:

    # GenActors = ASDVTKGenActors.ASDGenActors()
    # EneActors = ASDVTKEneActors.ASDEneActors()
    # MomActors = ASDVTKMomActors.ASDMomActors()
    # NeighActors = ASDVTKNeighActors.ASDNeighActors()

    lut = vtkLookupTable()

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
    # Toggle options for the contours
    ##########################################################################
    def toggle_contours(self, check):
        """
        Toggles the visibility of contour actors based on the check parameter.
        """
        if check:
            self.MomActors.contActor.VisibilityOn()
        else:
            self.MomActors.contActor.VisibilityOff()
        return

    ##########################################################################
    # Toggle the directions arrows
    ##########################################################################
    def toggle_directions(self, check):
        """
        Toggles the visibility of MomActors.vector based on the check parameter.
        """
        if check:
            self.MomActors.vector.VisibilityOn()
        else:
            self.MomActors.vector.VisibilityOff()
        return

    ##########################################################################
    # Toggle the directions arrows
    ##########################################################################
    def toggle_spins(self, check):
        """
        Toggles the visibility of the spin actors based on the check value.
        """
        if check:
            self.MomActors.Spins.VisibilityOn()
        else:
            self.MomActors.Spins.VisibilityOff()
        return

    ##########################################################################
    # Toggle the atomic spheres
    ##########################################################################
    def toggle_atoms(self, check):
        """
        Toggles the visibility of atoms based on the given check value.
        """
        if check:
            self.MomActors.Atoms.VisibilityOn()
        else:
            self.MomActors.Atoms.VisibilityOff()
        return

    ##########################################################################
    # Toggle the magnetization density
    ##########################################################################
    def toggle_density(self, check):
        """
        Toggles the visibility of the magnetic density actor based on the check value.
        """
        if check:
            self.MomActors.MagDensActor.VisibilityOn()
        else:
            self.MomActors.MagDensActor.VisibilityOff()
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
    # Set the color of the magnetization density along the x axis projection
    ##########################################################################
    def set_projection(self, atype, axis):
        """
        Set the projection type and axis for visualization.
        """
        if atype == "density":
            self.MomActors.src.GetPointData().SetScalars(
                self.MomActors.glob_color[axis]
            )
        elif atype == "spins":
            self.MomActors.src_spins.GetPointData().SetScalars(
                self.MomActors.glob_color[axis]
            )
        return

    ##########################################################################
    # Set the size of the spins via the slider
    ##########################################################################
    def ChangeSpinsSize(self, value):
        """
        Adjusts the scale factor of the SpinMapper based on the provided value.
        """
        self.MomActors.SpinMapper.SetScaleFactor(0.50 * value / 10)
        return

    ##########################################################################
    # Set the interpolation of the spin glyphs
    ##########################################################################
    def ChangeSpinShade(self, renWin, keyword):
        """
        Change the shading model of the Spin and Atom actors based on the given keyword.
        """
        if keyword == "Flat":
            if hasattr(self.MomActors, "Spins"):
                self.MomActors.Spins.GetProperty().SetInterpolationToFlat()
        elif keyword == "Gouraud":
            if hasattr(self.MomActors, "Spins"):
                self.MomActors.Spins.GetProperty().SetInterpolationToGouraud()
            if hasattr(self.MomActors, "Atoms"):
                self.MomActors.Atoms.GetProperty().SetInterpolationToGouraud()
        elif keyword == "PBR":
            if hasattr(self.MomActors, "Spins"):
                self.MomActors.Spins.GetProperty().SetInterpolationToPBR()
                self.MomActors.Spins.GetProperty().SetMetallic(0.5)
            if hasattr(self.MomActors, "Atoms"):
                self.MomActors.Atoms.GetProperty().SetInterpolationToPBR()
                self.MomActors.Atoms.GetProperty().SetMetallic(0.5)
        elif keyword == "Phong":
            if hasattr(self.MomActors, "Spins"):
                self.MomActors.Spins.GetProperty().SetInterpolationToPhong()
            if hasattr(self.MomActors, "Atoms"):
                self.MomActors.Atoms.GetProperty().SetInterpolationToPhong()

        renWin.Render()
        return

    ##########################################################################
    # Set the ambient scattering of the spin glyphs
    ##########################################################################
    def RenAmbientUpdate(self, value, renWin):
        """
        Update the ambient property of MomActors.Spins and render the window.
        """
        if hasattr(self.MomActors, "Spins"):
            self.MomActors.Spins.GetProperty().SetAmbient(float(value * 0.02))
        renWin.Render()
        return

    ##########################################################################
    # Set the diffuse scattering of the spin glyphs
    ##########################################################################
    def RenDiffuseUpdate(self, value, renWin):
        """
        Updates the diffuse property of MomActors.Spins and renders the window.
        """
        if hasattr(self.MomActors, "Spins"):
            self.MomActors.Spins.GetProperty().SetDiffuse(float(value * 0.01))
        renWin.Render()
        return

    ##########################################################################
    # Set the specular scattering of the spin glyphs
    ##########################################################################
    def RenSpecularUpdate(self, value, renWin):
        """
        Updates the specular property of MomActors' Spins and renders the window.
        """
        if hasattr(self.MomActors, "Spins"):
            self.MomActors.Spins.GetProperty().SetSpecular(float(value * 0.01))
        renWin.Render()
        return

    ##########################################################################
    # Set the specular scattering of the spin glyphs
    ##########################################################################
    def RenSpecularPowerUpdate(self, value, renWin):
        """
        Updates the specular power of the Spins actor and renders the window.
        """
        if hasattr(self.MomActors, "Spins"):
            self.MomActors.Spins.GetProperty().SetSpecularPower(float(value))
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
        if hasattr(self.MomActors, "Spins"):
            self.MomActors.Spins.GetProperty().SetEmissiveFactor(emvec)
        if hasattr(self.MomActors, "Atoms"):
            self.MomActors.Atoms.GetProperty().SetEmissiveFactor(emvec)
        renWin.Render()
        return

    ##########################################################################
    # Set the PBR occlusion value of the spin glyphs
    ##########################################################################
    def PBROcclusionUpdate(self, value, ren, renWin):
        """
        Updates the occlusion strength for Spins and Atoms actors and renders the window.
        """
        if hasattr(self.MomActors, "Spins"):
            self.MomActors.Spins.GetProperty().SetOcclusionStrength(
                float(value * 0.01)
            )
        if hasattr(self.MomActors, "Atoms"):
            self.MomActors.Atoms.GetProperty().SetOcclusionStrength(
                float(value * 0.01)
            )
        renWin.Render()
        return

    ##########################################################################
    # Set the PBR roughness value of the spin glyphs
    ##########################################################################
    def PBRRoughnessUpdate(self, value, renWin):
        """
        Updates the roughness property for Spins and Atoms actors and renders the window.
        """
        if hasattr(self.MomActors, "Spins"):
            self.MomActors.Spins.GetProperty().SetRoughness(
                float(value * 0.01)
            )
        if hasattr(self.MomActors, "Atoms"):
            self.MomActors.Atoms.GetProperty().SetRoughness(
                float(value * 0.01)
            )
        renWin.Render()
        return

    ##########################################################################
    # Set the PBR metallic value of the spin glyphs
    ##########################################################################
    def PBRMetallicUpdate(self, value, renWin):
        """
        Updates the metallic property of MomActors' Spins and Atoms and renders the window.
        """
        if hasattr(self.MomActors, "Spins"):
            self.MomActors.Spins.GetProperty().SetMetallic(float(value * 0.01))
        if hasattr(self.MomActors, "Atoms"):
            self.MomActors.Atoms.GetProperty().SetMetallic(float(value * 0.01))
        renWin.Render()
        return

    ##########################################################################
    # Set the size of the atoms via the slider
    ##########################################################################
    def ChangeAtomsSize(self, value):
        """
        Adjusts the size of atoms in the visualization by scaling the atom mapper.
        """
        self.MomActors.AtomMapper.SetScaleFactor(1.00 * value / 10.0)
        return

    ##########################################################################
    # Set the size of the atoms via the slider
    ##########################################################################
    def ChangeAtomsOpaq(self, value):
        """
        Adjusts the opacity of atom actors based on the given value.
        """
        self.MomActors.Atoms.GetProperty().SetOpacity(value * 0.01)
        return

    ##########################################################################
    # Set the quality of the atoms via the slider
    ##########################################################################
    def ChangeAtomsQuali(self, value):
        """
        Adjusts the resolution of the atom sphere visualization.
        """
        self.MomActors.AtomSphere.SetThetaResolution(value)
        self.MomActors.AtomSphere.SetPhiResolution(value)
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
    # @brief Function to find the needed file names for the HDR file
    ##########################################################################
    def getHDRIFileName(self, window):
        """
        Open a file dialog to select an existing HDR image file.
        """

        dlg = QtWidgets.QFileDialog()
        dlg.setFileMode(QtWidgets.QFileDialog.FileMode.ExistingFile)
        hdrifile = dlg.getOpenFileName(
            caption="Open HDR file", directory=".", filter="HDR images (*.pic *.hdr)"
        )[0]
        return hdrifile

    ##########################################################################
    # @brief Function to find the needed file names for the texture files
    ##########################################################################
    def getTextureFileName(self, window):
        """
        Opens a file dialog to select an existing texture file and returns its file path.
        """

        dlg = QtWidgets.QFileDialog()
        dlg.setFileMode(QtWidgets.QFileDialog.FileMode.ExistingFile)
        texturefile = dlg.getOpenFileName(
            caption="Open texture file", directory=".", filter="Images (*.png)"
        )[0]
        return texturefile

    ##########################################################################
    # Toggle surface texture
    ##########################################################################
    def toggle_Texture(self, check, ren, renWin, texfile):
        """
        Toggles the texture on the MomActors' Spins property based on the check flag.
        """

        if check:

            albedo_reader = vtk.vtkPNGReader()
            albedo_reader.SetFileName(texfile)

            albedo = vtk.vtkTexture()
            albedo.SetInputConnection(albedo_reader.GetOutputPort())
            albedo.UseSRGBColorSpaceOn()
            albedo.InterpolateOn()
            albedo.MipmapOn()

            self.MomActors.Spins.GetProperty().SetBaseColorTexture(albedo)
        else:
            self.MomActors.Spins.GetProperty().RemoveTexture("albedoTex")

        renWin.Render()

        return

    ##########################################################################
    # Toggle ORM texture
    ##########################################################################
    def toggle_ORMTexture(self, check, ren, renWin, texfile):
        """
        Toggles the ORM texture on or off for the Spins actor based on the check flag.
        """

        if check:

            material_reader = vtk.vtkPNGReader()
            material_reader.SetFileName(texfile)

            material = vtk.vtkTexture()
            material.SetInputConnection(material_reader.GetOutputPort())
            material.InterpolateOn()
            material.MipmapOn()

            self.MomActors.Spins.GetProperty().SetORMTexture(material)
        else:
            self.MomActors.Spins.GetProperty().RemoveTexture("materialTex")

        renWin.Render()

        return

    ##########################################################################
    # Toggle anisotropy texture
    ##########################################################################
    def toggle_ATexture(self, check, ren, renWin, texfile):
        """
        Toggles the anisotropy texture on or off based on the check parameter.
        """

        if check:

            anisotropy_reader = vtk.vtkPNGReader()
            anisotropy_reader.SetFileName(texfile)

            anisotropy = vtk.vtkTexture()
            anisotropy.SetInputConnection(anisotropy_reader.GetOutputPort())
            anisotropy.InterpolateOn()
            anisotropy.MipmapOn()

            self.MomActors.Spins.GetProperty().SetAnisotropyTexture(anisotropy)

        else:
            self.MomActors.Spins.GetProperty().RemoveTexture("anisotropyTex")

        renWin.Render()

        return

    ##########################################################################
    # Toggle normal texture
    ##########################################################################
    def toggle_NTexture(self, check, ren, renWin, texfile):
        """
        Toggles the application of a normal texture to the MomActors' Spins based on the check flag.
        """

        if check:

            normal_reader = vtk.vtkPNGReader()
            normal_reader.SetFileName(texfile)

            normal = vtk.vtkTexture()
            normal.InterpolateOn()
            normal.MipmapOn()
            normal.SetInputConnection(normal_reader.GetOutputPort())

            self.MomActors.SpinShader = self.MomActors.Spins.GetShaderProperty()

            self.MomActors.SpinShader.AddVertexShaderReplacement(
                "//VTK::Normal::Dec",  # replace the normal block
                True,  # before the standard replacements
                "//VTK::Normal::Dec\n"  # we still want the default
                "in vec3 tangentMC;\n"
                "out vec3 tangentVCVSOutput;\n",
                False,  # only do it once
            )
            self.MomActors.SpinShader.AddVertexShaderReplacement(
                "//VTK::Normal::Impl",  # replace the normal block
                True,  # before the standard replacements
                "//VTK::Normal::Impl\n"  # we still want the default
                "  tangentVCVSOutput = normalMatrix * tangentMC;\n",
                False,  # only do it once
            )
            self.MomActors.Spins.GetProperty().SetNormalTexture(normal)

        else:
            self.MomActors.SpinShader.ClearAllVertexShaderReplacements()
            self.MomActors.Spins.GetProperty().RemoveTexture("normalTex")

        renWin.Render()

        return

    ##########################################################################
    # Toggle ORM texture
    ##########################################################################
    def toggle_ETexture(self, check, ren, renWin, texfile):
        """
        Toggles the emissive texture on or off for the MomActors' Spins property.
        """

        if check:

            emissive_reader = vtk.vtkPNGReader()
            emissive_reader.SetFileName(texfile)

            emissive = vtk.vtkTexture()
            emissive.SetInputConnection(emissive_reader.GetOutputPort())
            emissive.UseSRGBColorSpaceOn()
            emissive.InterpolateOn()
            emissive.MipmapOn()

            self.MomActors.Spins.GetProperty().SetEmissiveTexture(emissive)

        else:
            self.MomActors.Spins.GetProperty().RemoveTexture("emissiveTex")

        renWin.Render()

        return

    ##########################################################################
    # Toggle Skybox on/off
    ##########################################################################
    def toggle_SkyBox(self, check, ren, renWin, skyboxfile):
        """
        Toggles the visibility of a skybox in the VTK renderer.
        """

        if check:
            reader = vtkHDRReader()
            reader.SetFileName(skyboxfile)
            reader.Update()

            texture = vtkTexture()
            texture.InterpolateOn()
            texture.MipmapOn()
            texture.SetColorModeToDirectScalars()
            texture.SetInputConnection(reader.GetOutputPort())

            self.MomActors.SkyBox.SetTexture(texture)
            self.MomActors.SkyBox.SetProjectionToSphere()
            self.MomActors.SkyBox.VisibilityOn()
        else:
            self.MomActors.SkyBox.VisibilityOff()

        renWin.Render()

        return

    ##########################################################################
    # Toggle HDRI on/off
    ##########################################################################
    def toggle_HDRI(self, check, ren, renWin, hdrifile):
        """
        Toggles HDRI (High Dynamic Range Imaging) for the given renderer and render window.
        """
        # import vtk

        if check:
            reader = vtkHDRReader()
            reader.SetFileName(hdrifile)
            reader.Update()

            texture = vtkTexture()
            texture.InterpolateOn()
            texture.MipmapOn()
            texture.SetColorModeToDirectScalars()
            texture.SetInputConnection(reader.GetOutputPort())

            # ren.RemoveAllLights()
            ren.AutomaticLightCreationOff()
            ren.UseImageBasedLightingOn()
            ren.UseSphericalHarmonicsOn()
            ren.SetEnvironmentTexture(texture)

        else:
            ren.UseSphericalHarmonicsOff()
            ren.UseImageBasedLightingOff()
            ren.AutomaticLightCreationOn()

        renWin.Render()
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
    # Toggle shadows on/off
    ##########################################################################
    # def toggle_Shadows(self,check, ren, renWin):
    # from vtkmodules.vtkRenderingOpenGL2 import (
    # vtkCameraPass,
    # vtkOpaquePass,
    # vtkRenderPassCollection,
    # vtkSequencePass,
    # vtkShadowMapPass
    # )
    # print('Toggle shadows', check)
    # if check:
    # seq = vtkSequencePass()

    # passes = vtkRenderPassCollection()

    # shadows = vtkShadowMapPass()
    # passes.AddItem(shadows.GetShadowMapBakerPass())
    # passes.AddItem(shadows)

    # opaque = vtkOpaquePass()
    # passes.AddItem(opaque)

    # seq.SetPasses(passes)

    # camera_p = vtkCameraPass()
    # camera_p.SetDelegatePass(seq)

    # Tell the renderer to use our render pass pipeline.
    # ren.SetPass(camera_p)
    # renWin.Render()

    # return
    ##########################################################################
    # Update glyph resolution
    ##########################################################################
    def GlyphQualityUpdate(self, value, viz_type, mode, renWin):
        """
        Updates the glyph quality for different visualization types and renders the window.
        """
        if viz_type == "M":
            try:
                self.MomActors.spinarrow.SetTipResolution(value)
                self.MomActors.spinarrow.SetShaftResolution(value)
            except AttributeError:
                pass
            try:
                self.MomActors.spinsphere.SetThetaResolution(value)
                self.MomActors.spinsphere.SetPhiResolution(value)
            except AttributeError:
                pass
            try:
                self.MomActors.spincones.SetResolution(value)
            except AttributeError:
                pass

        if viz_type == "N":
            if mode == 1:
                self.NeighActors.NeighGlyphs.SetThetaResolution(value)
                self.NeighActors.NeighGlyphs.SetPhiResolution(value)
            if mode == 2:
                self.NeighActors.NeighGlyphs.SetTipResolution(value)
                self.NeighActors.NeighGlyphs.SetShaftResolution(value)
        if viz_type == "E":
            self.EneActors.EneAtom.SetThetaResolution(value)
            self.EneActors.EneAtom.SetPhiResolution(value)

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

    def ChangeSpinGlyph(self, renWin, keyword):
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
                del self.MomActors.spinarrow
            except AttributeError:
                pass
            try:
                del self.MomActors.spinsphere
            except AttributeError:
                pass
            try:
                del self.MomActors.spincones
            except AttributeError:
                pass
            self.MomActors.spincube = vtk.vtkCubeSource()
            self.MomActors.spincube.SetXLength(1.0)
            self.MomActors.spincube.SetYLength(1.0)
            self.MomActors.spincube.SetZLength(1.0)

            # Calculate TCoords for texturing
            # self.MomActors.spincubetmap = vtk.vtkTextureMapToSphere()
            # self.MomActors.spincubetmap.SetInputConnection(self.MomActors.spincube.GetOutputPort())
            # self.MomActors.spincubetmap.PreventSeamOn()

            self.MomActors.SpinMapper.SetSourceConnection(
                self.MomActors.spincube.GetOutputPort()
            )
            # self.MomActors.SpinMapper.SetSourceConnection(self.MomActors.spincubetmap.GetOutputPort())
            self.MomActors.SpinMapper.ClampingOn()
            self.MomActors.SpinMapper.OrientOff()
            renWin.Render()
        if keyword == "Bars":
            try:
                del self.MomActors.spinarrow
            except AttributeError:
                pass
            try:
                del self.MomActors.spinsphere
            except AttributeError:
                pass
            try:
                del self.MomActors.spincones
            except AttributeError:
                pass
            self.MomActors.spincube = vtk.vtkCubeSource()
            self.MomActors.spincube.SetXLength(2.0)
            self.MomActors.spincube.SetYLength(0.4)
            self.MomActors.spincube.SetZLength(0.4)

            # Calculate TCoords for texturing
            self.MomActors.spincubetmap = vtk.vtkTextureMapToCylinder()
            self.MomActors.spincubetmap.SetInputConnection(
                self.MomActors.spincube.GetOutputPort()
            )
            # self.MomActors.spincubetmap.AutomaticCylinderGenerationOff()
            # self.MomActors.spincubetmap.SetPoint1([ 1.0,0.0,0.0])
            # self.MomActors.spincubetmap.SetPoint2([-1.0,0.0,0.0])

            # Calculate TCoords for texturing
            # self.MomActors.spincubetmap = vtk.vtkTextureMapToSphere()
            # self.MomActors.spincubetmap.SetInputConnection(self.MomActors.spincube.GetOutputPort())
            # self.MomActors.spincubetmap.AutomaticSphereGenerationOff()
            # self.MomActors.spincubetmap.SetCenter([ 0.0,0.0,0.0])

            self.MomActors.spincubetmap.PreventSeamOff()

            self.MomActors.SpinMapper.SetSourceConnection(
                self.MomActors.spincubetmap.GetOutputPort()
            )
            self.MomActors.SpinMapper.ClampingOn()
            self.MomActors.SpinMapper.OrientOn()
            renWin.Render()

        if keyword == "Spheres":
            try:
                del self.MomActors.spinarrow
            except AttributeError:
                pass
            try:
                del self.MomActors.spincube
            except AttributeError:
                pass
            try:
                del self.MomActors.spincones
            except AttributeError:
                pass
            self.MomActors.spinsphere = vtk.vtkTexturedSphereSource()
            self.MomActors.spinsphere.SetRadius(0.50)
            self.MomActors.spinsphere.SetThetaResolution(12)
            self.MomActors.spinsphere.SetPhiResolution(12)

            # Placeholder comment for testing tangent extraction for normal textures
            # tritri = vtk.vtkTriangleFilter()
            # tritri.SetInputConnection(self.MomActors.spinsphere.GetOutputPort())
            # tritan = vtk.vtkPolyDataTangents()
            # tritan.SetInputConnection(tritri.GetOutputPort())
            # self.MomActors.SpinMapper.SetSourceConnection(tritan.GetOutputPort())

            self.MomActors.SpinMapper.SetSourceConnection(
                self.MomActors.spinsphere.GetOutputPort()
            )
            self.MomActors.SpinMapper.ClampingOn()
            self.MomActors.SpinMapper.OrientOn()
            # self.MomActors.SpinMapper.OrientOff()
            renWin.Render()
        if keyword == "Arrows":
            try:
                del self.MomActors.spinsphere
            except AttributeError:
                pass
            try:
                del self.MomActors.spincube
            except AttributeError:
                pass
            try:
                del self.MomActors.spincones
            except AttributeError:
                pass

            # Create vectors
            self.MomActors.spinarrow = vtk.vtkArrowSource()
            self.MomActors.spinarrow.SetTipRadius(0.20)
            self.MomActors.spinarrow.SetShaftRadius(0.10)
            self.MomActors.spinarrow.SetTipResolution(12)
            self.MomActors.spinarrow.SetShaftResolution(12)

            # Calculate normals for shading
            self.MomActors.spinarrownormals = vtk.vtkPolyDataNormals()
            self.MomActors.spinarrownormals.SetInputConnection(
                self.MomActors.spinarrow.GetOutputPort()
            )

            # Calculate TCoords for texturing
            self.MomActors.spinarrownormalstmap = vtk.vtkTextureMapToCylinder()
            self.MomActors.spinarrownormalstmap.SetInputConnection(
                self.MomActors.spinarrownormals.GetOutputPort()
            )
            self.MomActors.spinarrownormalstmap.PreventSeamOn()

            self.MomActors.SpinMapper.SetSourceConnection(
                self.MomActors.spinarrownormalstmap.GetOutputPort()
            )
            self.MomActors.SpinMapper.OrientOn()
            self.MomActors.SpinMapper.Update()

            renWin.Render()
        if keyword == "CenterOn":
            self.MomActors.spinarrow.SetArrowOriginToCenter()
            renWin.Render()

        if keyword == "CenterOff":
            self.MomActors.spinarrow.SetArrowOriginToDefault()
            renWin.Render()

        if keyword == "Cones":
            try:
                del self.MomActors.spinsphere
            except AttributeError:
                pass
            try:
                del self.MomActors.spincube
            except AttributeError:
                pass
            try:
                del self.MomActors.spinarrow
            except AttributeError:
                pass

            self.MomActors.spincones = vtk.vtkConeSource()
            self.MomActors.spincones.SetRadius(0.50)
            self.MomActors.spincones.SetHeight(1.00)
            self.MomActors.spincones.SetResolution(12)

            # Calculate normals for shading
            self.MomActors.spinconenormals = vtk.vtkPolyDataNormals()
            self.MomActors.spinconenormals.SetInputConnection(
                self.MomActors.spincones.GetOutputPort()
            )

            # Calculate TCoords for texturing
            self.MomActors.spinconeormalstmap = vtk.vtkTextureMapToCylinder()
            self.MomActors.spinconeormalstmap.SetInputConnection(
                self.MomActors.spinconenormals.GetOutputPort()
            )
            self.MomActors.spinconeormalstmap.PreventSeamOn()

            self.MomActors.SpinMapper.SetSourceConnection(
                self.MomActors.spinconeormalstmap.GetOutputPort()
            )
            self.MomActors.SpinMapper.OrientOn()
            self.MomActors.SpinMapper.Update()
            renWin.Render()
        return

