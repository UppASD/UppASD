"""@package ASDVTKTexture
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
from vtkmodules.vtkIOImage import vtkHDRReader
from vtkmodules.vtkRenderingCore import vtkTexture


##########################################################################
# @brief Class containing the majority of the actions to update the visualizer
# @details Class containing the majority of the actions to update the visualizer.
# It handles the taking of snapshots, the toggling of projections, actor visibility
# and colormaps. It also controls many of the camera options so that they can be
# updated during runtime, allowing for finer control of the visualization.
# @author Jonathan Chico
##########################################################################
class ASDTexture:
    """
    ASDTexture class provides methods to manage and toggle various textures and HDRI for VTK actors.

    Methods:
        __init__(): Initializes texture attributes.
        gather_settings(): Gathers texture settings into a dictionary.
        getHDRIFileName(window): Opens a file dialog to select an HDR image file.
        getTextureFileName(window): Opens a file dialog to select a texture file.
        toggle_Texture(check, actor, texfile): Toggles the albedo texture on or off.
        toggle_ORMTexture(check, actor, texfile): Toggles the ORM texture on or off.
        toggle_ATexture(check, actor, texfile): Toggles the anisotropy texture on or off.
        toggle_NTexture(check, actor, texfile): Toggles the normal texture on or off.
        toggle_ETexture(check, actor, texfile): Toggles the emissive texture on or off.
        toggle_HDRI(check, ren, renWin, hdrifile): Toggles HDRI for the renderer and render window.
        toggle_SkyBox(check, actor, skyboxfile): Toggles the visibility of a skybox.
    """
    # lut = vtkLookupTable()
    def __init__(self):
        self.albedoTex = None
        self.materialTex = None
        self.anisotropyTex = None
        self.normalTex = None
        self.emissiveTex = None
        self.skybox = False

    def gather_settings(self):
        """
        Gathers texture settings into the provided dictionary.
        """
        settings = {}
        settings["albedoTex"] = self.albedoTex
        settings["materialTex"] = self.materialTex
        settings["anisotropyTex"] = self.anisotropyTex
        settings["normalTex"] = self.normalTex
        settings["emissiveTex"] = self.emissiveTex
        settings["skybox"] = self.skybox

        return settings

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
    def toggle_Texture(self, check, actor, texfile):
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

            actor.Spins.GetProperty().SetBaseColorTexture(albedo)
        else:
            actor.Spins.GetProperty().RemoveTexture("albedoTex")

        return

    ##########################################################################
    # Toggle ORM texture
    ##########################################################################
    def toggle_ORMTexture(self, check, actor, texfile):
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

            actor.Spins.GetProperty().SetORMTexture(material)
        else:
            actor.Spins.GetProperty().RemoveTexture("materialTex")

        return

    ##########################################################################
    # Toggle anisotropy texture
    ##########################################################################
    def toggle_ATexture(self, check, actor, texfile):
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

            actor.Spins.GetProperty().SetAnisotropyTexture(anisotropy)

        else:
            actor.Spins.GetProperty().RemoveTexture("anisotropyTex")

        return

    ##########################################################################
    # Toggle normal texture
    ##########################################################################
    def toggle_NTexture(self, check, actor, texfile):
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

            actor.SpinShader = actor.Spins.GetShaderProperty()

            actor.SpinShader.AddVertexShaderReplacement(
                "//VTK::Normal::Dec",  # replace the normal block
                True,  # before the standard replacements
                "//VTK::Normal::Dec\n"  # we still want the default
                "in vec3 tangentMC;\n"
                "out vec3 tangentVCVSOutput;\n",
                False,  # only do it once
            )
            actor.SpinShader.AddVertexShaderReplacement(
                "//VTK::Normal::Impl",  # replace the normal block
                True,  # before the standard replacements
                "//VTK::Normal::Impl\n"  # we still want the default
                "  tangentVCVSOutput = normalMatrix * tangentMC;\n",
                False,  # only do it once
            )
            actor.Spins.GetProperty().SetNormalTexture(normal)

        else:
            actor.SpinShader.ClearAllVertexShaderReplacements()
            actor.Spins.GetProperty().RemoveTexture("normalTex")

        return

    ##########################################################################
    # Toggle ORM texture
    ##########################################################################
    def toggle_ETexture(self, check, actor, texfile):
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

            actor.Spins.GetProperty().SetEmissiveTexture(emissive)

        else:
            actor.Spins.GetProperty().RemoveTexture("emissiveTex")

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

    def toggle_SkyBox(self, check, actor, skyboxfile):
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

            actor.SkyBox.SetTexture(texture)
            actor.SkyBox.SetProjectionToSphere()
            actor.SkyBox.VisibilityOn()
            self.skybox = True
        else:
            actor.SkyBox.VisibilityOff()
            self.skybox = False

        return
