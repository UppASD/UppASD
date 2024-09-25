"""
ASDUISettings module

This module provides the ASDUISettings class, which manages and stores settings
for the ASD GUI application.
"""
# pylint: disable=invalid-name, no-name-in-module, no-member

import json
import os
import yaml


# Add a constructor for tuples
def tuple_constructor(loader, node):
    """
    Constructs a tuple from a YAML sequence node.
    """
    return tuple(loader.construct_sequence(node))


yaml.SafeLoader.add_constructor("tag:yaml.org,2002:python/tuple", tuple_constructor)


class ASDUISettings:
    """
    A class to manage and store settings for the ASD GUI application.

    Methods
    -------
    __init__():
        Initializes the settings dictionary.
    gather_dicts(window):
        Gathers settings from various components of the window.
    write_to_json(file_path):
        Writes the settings to a JSON file.
    write_to_yaml(file_path):
        Writes the settings to a YAML file.
    read_from_json(file_path):
        Reads the settings from a JSON file.
    read_from_yaml(file_path):
        Reads the settings from a YAML file.
    """

    def __init__(self):
        """
        Initializes the ASDUISettings class with default settings.
        """
        self.settings = {}

    def gather_dicts(self, window):
        """
        Gathers various settings from the provided window object into a dictionary.

        Args:
            window: The window object containing settings for MomActors, ASDTexture, ASDColor, etc.

        Sets:
            self.settings (dict): A dictionary containing all gathered settings.
        """
        # The MomActors settings
        mom_settings = {
            "spin_resolution": window.MomActors.spin_resolution,
            "atom_resolution": window.MomActors.atom_resolution,
            "spin_glyph": window.MomActors.spin_glyph,
            "center_on": window.MomActors.center_on,
            "contour": window.MomActors.contour,
            "vector_directions": window.MomActors.vector_directions,
            "show_spins": window.MomActors.show_spins,
            "show_atoms": window.MomActors.show_atoms,
            "show_density": window.MomActors.show_density,
            "projection_type": window.MomActors.projection_type,
            "projection_vector": window.MomActors.projection_vector,
            "spin_size": window.MomActors.spin_size,
            "spin_shade": window.MomActors.spin_shade,
            "ambient": window.MomActors.ambient,
            "specular": window.MomActors.specular,
            "specular_power": window.MomActors.specular_power,
            "diffuse": window.MomActors.diffuse,
            "pbr_emission": window.MomActors.pbr_emission,
            "pbr_occlusion": window.MomActors.pbr_occlusion,
            "pbr_roughness": window.MomActors.pbr_roughness,
            "pbr_metallic": window.MomActors.pbr_metallic,
            "atom_size": window.MomActors.atom_size,
            "atom_opacity": window.MomActors.atom_opacity,
        }
        # The ASDVTKTexture settings
        tex_settings = {
            "albedo": window.ASDTexture.albedo,
            "albedoTex": window.ASDTexture.albedoTex,
            "material": window.ASDTexture.material,
            "materialTex": window.ASDTexture.materialTex,
            "anisotropy": window.ASDTexture.anisotropy,
            "anisotropyTex": window.ASDTexture.anisotropyTex,
            "emissive": window.ASDTexture.emissive,
            "emissiveTex": window.ASDTexture.emissiveTex,
            "normal": window.ASDTexture.normal,
            "normalTex": window.ASDTexture.normalTex,
            "skybox": window.ASDTexture.skybox,
            "skyboxfile": window.ASDTexture.skyboxfile,
            "hdri": window.ASDTexture.hdri,
            "hdrifile": window.ASDTexture.hdrifile,
        }
        # The ASDVizOptions settings
        viz_settings = {
            "ssao": window.ASDVizOpt.ssao,
            "fxxa": window.ASDVizOpt.fxxa,
            "axes": window.ASDVizOpt.axes,
            "colorbar": window.ASDVizOpt.colorbar,
            "time_label": window.ASDVizOpt.time_label,
            "focal_disc": window.ASDVizOpt.focal_disc,
            "auto_focus": window.ASDVizOpt.auto_focus,
            "focal_blur": window.ASDVizOpt.focal_blur,
        }
        # The ASDColor settings
        color_settings = {
            "num_colors": window.ASDColor.num_colors,
            "mapnum": window.ASDColor.mapnum,
            "rgb": window.ASDColor.rgb,
            "rgb_background": window.ASDColor.rgb_background,
            "lut_scale": window.ASDColor.lut_scale,
        }

        # The ASDCamera settings
        cam_settings = {
            "position": window.ASDCamera.camera.GetPosition(),
            "focal_point": window.ASDCamera.camera.GetFocalPoint(),
            "view_up": window.ASDCamera.camera.GetViewUp(),
            "clipping_range": window.ASDCamera.camera.GetClippingRange(),
            "parallel_scale": window.ASDCamera.camera.GetParallelScale(),
            "is_parallel": window.ASDCamera.camera.GetParallelProjection(),
        }

        # Gather all settings into the settings dictionary
        self.settings = {
            "MomActors": mom_settings,
            "ASDVTKTexture": tex_settings,
            "ASDVizOptions": viz_settings,
            "ASDColors": color_settings,
            "ASDCamera": cam_settings,
        }

    def write_to_json(self, file_path):
        """
        Write the settings to a JSON file.

        Args:
            file_path (str): The path to the JSON file.
        """
        with open(file_path, "w", encoding="utf-8") as json_file:
            json.dump(self.settings, json_file, indent=4)

    def write_to_yaml(self, file_path):
        """
        Write the settings to a YAML file.

        Args:
            file_path (str): The path to the YAML file to write.
        """
        with open(file_path, "w", encoding="utf-8") as yaml_file:
            yaml.dump(self.settings, yaml_file, default_flow_style=False)

    def read_from_json(self, file_path):
        """
        Reads settings from a JSON file and updates the settings attribute.

        Args:
            file_path (str): Path to the JSON file.
        """
        if os.path.exists(file_path):
            with open(file_path, "r", encoding="utf-8") as json_file:
                self.settings = json.load(json_file)
        else:
            raise FileNotFoundError(f"No such file: '{file_path}'")

    def read_from_yaml(self, file_path):
        """
        Reads settings from a YAML file and loads them into the instance.

        Args:
            file_path (str): Path to the YAML file.
        """
        if os.path.exists(file_path):
            with open(file_path, "r", encoding="utf-8") as yaml_file:
                self.settings = yaml.safe_load(yaml_file)
        else:
            raise FileNotFoundError(f"No such file: '{file_path}'")

    def restore_from_settings(self, window):
        """
        Restore display settings from a saved configuration to the provided window object.
        """
        print("Restoring settings from file")

        # Restore ASDCamera settings
        cam_settings = self.settings.get("ASDCamera", {})
        if "position" in cam_settings:
            window.ASDCamera.camera.SetPosition(cam_settings["position"])
        if "focal_point" in cam_settings:
            window.ASDCamera.camera.SetFocalPoint(cam_settings["focal_point"])
        if "view_up" in cam_settings:
            window.ASDCamera.camera.SetViewUp(cam_settings["view_up"])
        if "clipping_range" in cam_settings:
            window.ASDCamera.camera.SetClippingRange(cam_settings["clipping_range"])
        if "parallel_scale" in cam_settings:
            window.ASDCamera.camera.SetParallelScale(cam_settings["parallel_scale"])
        if "is_parallel" in cam_settings:
            window.ASDCamera.camera.SetParallelProjection(cam_settings["is_parallel"])

        # Restore ASDColor settings
        color_settings = self.settings.get("ASDColors", {})
        if (
            "num_colors" in color_settings
        ):
            window.ASDColor.set_lut_db(window, color_settings.get("mapnum"))
        if "rgb" in color_settings and window.ASDColor.num_colors == 1:
            window.ASDColor.set_RGBcolor(color_settings["rgb"], window.viz_type)
        if "rgb_background" in color_settings:
            window.ASDColor.set_RGBbackground(
                rgb=color_settings["rgb_background"], ren=window.ren
            )

        # Restore ASDVTKTexture settings
        tex_settings = self.settings.get("ASDVTKTexture", {})
        if tex_settings.get("albedo") and tex_settings.get("albedoTex") is not None:
            window.ASDTexture.toggle_Texture(
                True, window.MomActors, tex_settings["albedoTex"]
            )
        if tex_settings.get("material") and tex_settings.get("materialTex") is not None:
            window.ASDTexture.toggle_ORMTexture(
                True, window.MomActors, tex_settings["materialTex"]
            )
        if (
            tex_settings.get("anisotropy")
            and tex_settings.get("anisotropyTex") is not None
        ):
            window.ASDTexture.toggle_ATexture(
                True, window.MomActors, tex_settings["anisotropyTex"]
            )
        if tex_settings.get("normal") and tex_settings.get("normalTex") is not None:
            window.ASDTexture.toggle_NTexture(
                True, window.MomActors, tex_settings["normalTex"]
            )
        if tex_settings.get("emissive") and tex_settings.get("emissiveTex") is not None:
            window.ASDTexture.toggle_ETexture(
                True, window.MomActors, tex_settings["emissiveTex"]
            )
        if tex_settings.get("skybox") and tex_settings.get("skyboxfile") is not None:
            window.ASDTexture.toggle_Texture(
                True, window.MomActors, tex_settings["skyboxfile"]
            )
        if tex_settings.get("hdri") and tex_settings.get("hdrifile") is not None:
            window.ASDTexture.toggle_HDRI(
                True, window.ren, window.renWin, tex_settings["hdrifile"]
            )

        # Restore MomActors settings
        mom_settings = self.settings.get("MomActors", {})
        if "spin_resolution" in mom_settings:
            window.MomActors.spin_resolution = mom_settings["spin_resolution"]
        if "atom_resolution" in mom_settings:
            window.MomActors.atom_resolution = mom_settings["atom_resolution"]
        if "spin_glyph" in mom_settings:
            window.MomActors.ChangeSpinGlyph(mom_settings["spin_glyph"])
        if "center_on" in mom_settings:
            window.MomActors.ChangeSpinGlyph(mom_settings["center_on"])
        if "contour" in mom_settings:
            window.MomActors.toggle_contours(mom_settings["contour"])
        # if "vector_directions" in mom_settings:
        #     window.MomActors.toggle_directions(mom_settings["vector_directions"])
        if "show_spins" in mom_settings:
            window.MomActors.toggle_spins(mom_settings["show_spins"])
        if "show_atoms" in mom_settings:
            window.MomActors.toggle_atoms(mom_settings["show_atoms"])
        if "show_density" in mom_settings:
            window.MomActors.toggle_density(mom_settings["show_density"])
        if "projection_type" in mom_settings and "projection_vector" in mom_settings:
            window.MomActors.set_projection(
                mom_settings["projection_type"], mom_settings["projection_vector"]
            )
        if "spin_size" in mom_settings:
            window.MomActors.ChangeSpinsSize(mom_settings["spin_size"])
        if "spin_shade" in mom_settings:
            window.MomActors.ChangeSpinShade(mom_settings["spin_shade"])
        if "ambient" in mom_settings:
            window.MomActors.RenAmbientUpdate(
                mom_settings["ambient"] / 0.02, window.renWin
            )
        if "diffuse" in mom_settings:
            window.MomActors.RenDiffuseUpdate(
                mom_settings["diffuse"] / 0.01, window.renWin
            )
        if "specular" in mom_settings:
            window.MomActors.RenSpecularUpdate(
                mom_settings["specular"] / 0.01, window.renWin
            )
        if "specular_power" in mom_settings:
            window.MomActors.RenSpecularPowerUpdate(
                mom_settings["specular_power"] / 0.01, window.renWin
            )
        if "pbr_emission" in mom_settings:
            window.MomActors.PBREmissionUpdate(
                mom_settings["pbr_emission"] / 0.01, window.ren, window.renWin
            )
        if "pbr_occlusion" in mom_settings:
            window.MomActors.PBROcclusionUpdate(
                mom_settings["pbr_occlusion"] / 0.01, window.ren, window.renWin
            )
        if "pbr_roughness" in mom_settings:
            window.MomActors.PBRRoughnessUpdate(
                mom_settings["pbr_roughness"] / 0.01, window.renWin
            )
        if "pbr_metallic" in mom_settings:
            window.MomActors.PBRMetallicUpdate(
                mom_settings["pbr_metallic"] / 0.01, window.renWin
            )
        if "atom_size" in mom_settings:
            window.MomActors.ChangeAtomsSize(mom_settings["atom_size"] * 10.0)
        if "atom_opacity" in mom_settings:
            window.MomActors.ChangeAtomsOpaq(mom_settings["atom_opacity"])

        # Restore ASDVizOptions settings
        viz_settings = self.settings.get("ASDVizOptions", {})
        if "ssao" in viz_settings:
            window.ASDVizOpt.toggle_SSAO(viz_settings["ssao"], window.ren)
        if "fxxa" in viz_settings:
            window.ASDVizOpt.toggle_FXAA(
                viz_settings["fxxa"], window.ren, window.renWin
            )
        if "axes" in viz_settings:
            window.ASDVizOpt.toggle_Axes(window, viz_settings["axes"])
        if "colorbar" in viz_settings:
            window.ASDVizOpt.toggle_ScalarBar(window, viz_settings["colorbar"])
        if "time_label" in viz_settings:
            window.ASDVizOpt.toggle_time_label(window, viz_settings["time_label"])
        if "focal_blur" in viz_settings:
            window.ASDVizOpt.toggle_Focus(
                viz_settings.get("focal_blur", False), window.ren, window.renWin
            )
            window.ASDVizOpt.toggle_autoFocus(
                viz_settings.get("auto_focus", False), window.renWin
            )
            window.ASDVizOpt.setFocalDisk(
                viz_settings.get("focal_disc", 100), window.ren, window.renWin
            )

        window.renWin.Render()