"""@package ASDVTKColor
Contains a class with a set of functions dealing with the visualization options
in the VTK visualization mode. Dealing with visibility of objects, size of glyphs,
types of glyphs, colormap used etc.

Author
----------
Jonathan Chico, Anders Bergman
"""
# pylint: disable=invalid-name, no-name-in-module, no-member

import vtk
from vtkmodules.vtkCommonColor import vtkColorSeries
from vtkmodules.vtkCommonCore import vtkLookupTable
##########################################################################
# @brief Class containing the majority of the actions to update the visualizer
# @details Class containing the majority of the actions to update the visualizer.
# It handles the taking of snapshots, the toggling of projections, actor visibility
# and colormaps. It also controls many of the camera options so that they can be
# updated during runtime, allowing for finer control of the visualization.
# @author Jonathan Chico
##########################################################################


class ASDVTKColor:
    """
    ASDVTKColor class provides methods to set and manage color maps and lookup tables (LUTs)
    for visualizations using the VTK library. It allows users to choose between different color
    schemes and apply them to various actors in the visualization.

    Methods:
        set_colormap(self, window, flag2D, viz_type, renWin):
            Selects the type of colormap to be used for different actors.
        set_colormap_db(self, mapnum, viz_type):
            Selects the type of colormap to be used for different actors based on a map number.
        set_RGBcolor(self, rgb, viz_type):
            Sets the RGB color for the lookup table (LUT) and color transfer function.
        set_RGBbackground(self, rgb, ren):
        set_singlecolor(self, window, value):
            Sets the single color for the RGB sliders and updates the visualization.
        set_background(self, window, value):
        set_lut_db(self, window, mapnum):
        set_lut_scale(self, window):
    """

    def __init__(self):
        """
        Initializes the ASDVTKColor class.
        """
        self.lut = vtk.vtkLookupTable()
        self.transfer_func = vtk.vtkColorTransferFunction()
        self.NUM_LUT_COLORS = 256
        self.NUM_RGB_COLORS = 1
        self.num_colors = self.NUM_LUT_COLORS
        self.mapnum = 1
        self.rgb = [0.0, 0.0, 0.0]
        self.rgb_background = [1.0, 1.0, 1.0]
        self.lut_scale = "Linear"

    def get_color_settings(self):
        """
        Stores the current settings of num_colors, mapnum, and rgb to a dictionary.

        Returns:
        dict: A dictionary containing the current settings.
        """
        settings = {
            "num_colors": self.num_colors,
            "mapnum": self.mapnum,
            "rgb": self.rgb,
            "rgb_background": self.rgb_background,
            "lut_scale": self.lut_scale,
        }
        return settings

    def set_color_settings(self, settings, window):
        """
        Sets the color settings based on the provided settings dictionary.

        Parameters:
        settings (dict): A dictionary containing the settings with
                         keys 'num_colors', 'mapnum', and 'rgb'.
        window (object): The window object where the color settings will be applied.

        Returns:
        None
        """
        self.num_colors = settings.get("num_colors", self.NUM_LUT_COLORS)
        self.mapnum = settings.get("mapnum", 1)
        self.rgb = settings.get("rgb", [0.0, 0.0, 0.0])
        self.rgb_background = settings.get("rgb_background", [1.0, 1.0, 1.0])
        self.lut_scale = settings.get("lut_scale", "Linear")

        if self.num_colors == 256:
            self.set_lut_db(window, self.mapnum)
        elif self.num_colors == 2:
            self.set_RGBcolor(self.rgb, window.viz_type)

    def Series_to_LUT(self, colorSeries):
        """
        Convert a vtkColorSeries to a vtkLookupTable.

        Parameters:
        colorSeries (object): A vtkColorSeries object.

        Returns:
        object: A vtkLookupTable object.
        """
        lut = vtkLookupTable()
        colorSeries.BuildLookupTable(lut, vtkColorSeries.ORDINAL)
        return lut

    def LUT_to_Series(self, lut):
        """
        Convert a vtkLookupTable to a vtkColorSeries.

        Parameters:
        lut (object): A vtkLookupTable object.

        Returns:
        object: A vtkColorSeries object.
        """
        colorSeries = vtkColorSeries()
        colorSeries.SetColorScheme(vtkColorSeries.ORDINAL)
        colorSeries.BuildLookupTable(lut, vtkColorSeries.ORDINAL)
        return colorSeries

    def LUT_to_TransferFunction(self, lut):
        """
        Convert a vtkLookupTable to a vtkColorTransferFunction.

        Parameters:
        lut (object): A vtkLookupTable object.

        Returns:
        object: A vtkColorTransferFunction object.
        """
        transfer_func = vtk.vtkColorTransferFunction()
        for ii, ss in enumerate(
            [
                float(xx) / float(lut.GetNumberOfTableValues())
                for xx in range(lut.GetNumberOfTableValues())
            ]
        ):
            cc = lut.GetTableValue(ii)
            transfer_func.AddRGBPoint(ss, cc[0], cc[1], cc[2])
        return transfer_func

    def Series_to_TransferFunction(self, colorSeries):
        """
        Convert a vtkColorSeries to a vtkColorTransferFunction.

        Parameters:
        colorSeries (object): A vtkColorSeries object.

        Returns:
        object: A vtkColorTransferFunction object.
        """
        lut = vtkLookupTable()
        colorSeries.BuildLookupTable(lut, vtkColorSeries.ORDINAL)
        transfer_func = vtk.vtkColorTransferFunction()
        for ii, ss in enumerate(
            [
                float(xx) / float(lut.GetNumberOfTableValues())
                for xx in range(lut.GetNumberOfTableValues())
            ]
        ):
            cc = lut.GetTableValue(ii)
            transfer_func.AddRGBPoint(ss, cc[0], cc[1], cc[2])
        return transfer_func

    def ColorActors(self, window, flag2D, viz_type, renWin):
        """
        Apply color settings to actors based on visualization type and dimensionality.

        Parameters:
        window (object): The window containing the actors.
        flag2D (bool): Flag indicating if the visualization is 2D.
        viz_type (str): Type of visualization ('M', 'N', 'E').
        renWin (object): The render window to update.

        Returns:
        None
        """
        if viz_type == "M":
            if flag2D:
                window.MomActors.MagDensMap.SetLookupTable(self.lut)
                # window.MomActors.SpinMapper.SetLookupTable(self.lut)
                # window.ASDGenActors.scalar_bar.SetLookupTable(self.lut)
                window.ASDGenActors.clipperMapper.SetLookupTable(self.lut)
            else:
                window.MomActors.volumeProperty.SetColor(self.transfer_func)
                # window.MomActors.SpinMapper.SetLookupTable(
                #     self.transfer_func
                # )
                # window.ASDGenActors.scalar_bar.SetLookupTable(
                #     self.transfer_func
                # )
                window.ASDGenActors.clipperMapper.SetLookupTable(self.transfer_func)
        elif viz_type == "N":
            window.NeighActors.NeighMapper.SetLookupTable(self.lut)
            window.ASDGenActors.scalar_bar.SetLookupTable(self.lut)
            window.ASDGenActors.clipperMapper.SetLookupTable(self.lut)
        elif viz_type == "E":
            if flag2D:
                window.EneActors.EneDensMap.SetLookupTable(self.lut)
                window.ASDGenActors.scalar_bar.SetLookupTable(self.lut)
                window.ASDGenActors.clipperMapper.SetLookupTable(self.lut)
            else:
                window.EneActors.volumeProperty.SetColor(self.transfer_func)
                window.ASDGenActors.scalar_bar.SetLookupTable(self.transfer_func)
                window.ASDGenActors.clipperMapper.SetLookupTable(self.transfer_func)
        # -----------------------------------------------------------------------
        # Render the scene
        # -----------------------------------------------------------------------
        renWin.Render()
        return

    ##########################################################################
    # @brief Select the type of colormap that will be used for the different actors
    # @details Select the type of colormap that will be used for the different actors
    # It allows the user to choose between the following color schemes:
    #
    #   - Coolwarm
    #   - RdGy
    #   - Spectral
    #   - BlackBody
    # @author Jonathan Chico
    ##########################################################################
    def set_colormap_db(self, mapnum, viz_type):
        """Select the type of colormap that will be used for the different actors
        It allows the user to choose between the following color schemes:
            * Coolwarm (mapnum 0)
            * RdGy (mapnum 1)
            * Spectral (mapnum 2)
            * BlackBody (mapnum 3)

        Args:
            window: QMainWindow where the visualizations are being carried out.
            flag2D: (logical) identifier indicating whether the system is in 2D or 3D.
            viz_type: (str) identifier for the different types of visualization
                possible in the VTK API.
            renWin: current VTK rendering window.

        Author
        ----------
        Jonathan Chico
        """
        # self.lut = vtk.vtkLookupTable()
        self.lut.SetNumberOfTableValues(self.num_colors)
        self.transfer_func = vtk.vtkColorTransferFunction()
        # -----------------------------------------------------------------------
        # Set the color map to be given by the diverging Coolwarm scheme by Kenneth Moreland
        # -----------------------------------------------------------------------
        if mapnum == 0:
            self.transfer_func.SetColorSpaceToDiverging()
            if (viz_type == "M") or (viz_type == "E") or viz_type == "N":
                self.transfer_func.AddRGBPoint(0, 0.230, 0.299, 0.754)
                self.transfer_func.AddRGBPoint(1, 0.706, 0.016, 0.150)
            else:
                self.transfer_func.AddRGBPoint(-1, 0.230, 0.299, 0.754)
                self.transfer_func.AddRGBPoint(1, 0.706, 0.016, 0.150)
        # -----------------------------------------------------------------------
        # Set the color to be given by the black body function
        # -----------------------------------------------------------------------
        if mapnum == 1:
            self.transfer_func.SetColorSpaceToRGB()
            if (viz_type == "M") or (viz_type == "E") or viz_type == "N":
                self.transfer_func.AddRGBPoint(0.00, 0.000, 0.000, 0.000)
                self.transfer_func.AddRGBPoint(0.39, 0.698, 0.133, 0.133)
                self.transfer_func.AddRGBPoint(0.58, 0.890, 0.412, 0.020)
                self.transfer_func.AddRGBPoint(0.89, 0.902, 0.902, 0.208)
                self.transfer_func.AddRGBPoint(1.00, 1.000, 1.000, 1.000)
            else:
                self.transfer_func.AddRGBPoint(-1.0, 0.0, 0.0, 0.0)
                self.transfer_func.AddRGBPoint(-0.5, 0.9, 0.0, 0.0)
                self.transfer_func.AddRGBPoint(0.5, 0.9, 0.9, 0.0)
                self.transfer_func.AddRGBPoint(1.0, 1.0, 1.0, 1.0)
        # -----------------------------------------------------------------------
        # Set the color map to be given by the diverging RdGy
        # -----------------------------------------------------------------------
        if mapnum == 2:
            self.transfer_func.SetColorSpaceToDiverging()
            if (viz_type == "M") or (viz_type == "E") or viz_type == "N":
                self.transfer_func.AddRGBPoint(0.0, 0.79216, 0.00000, 0.12549)
                self.transfer_func.AddRGBPoint(0.5, 1.00000, 1.00000, 1.00000)
                self.transfer_func.AddRGBPoint(1.0, 0.25098, 0.25098, 0.25098)
            else:
                self.transfer_func.AddRGBPoint(-1.0, 0.79216, 0.00000, 0.12549)
                self.transfer_func.AddRGBPoint(0.0, 1.00000, 1.00000, 1.00000)
                self.transfer_func.AddRGBPoint(1.0, 0.25098, 0.25098, 0.25098)
        # -----------------------------------------------------------------------
        # Set the color map to be given by the diverging spectral clor map
        # -----------------------------------------------------------------------
        if mapnum == 3:
            self.transfer_func.SetColorSpaceToRGB()
            if (viz_type == "M") or (viz_type == "E") or viz_type == "N":
                self.transfer_func.AddRGBPoint(0.00, 0.61961, 0.00392, 0.25882)
                self.transfer_func.AddRGBPoint(0.25, 0.95686, 0.42745, 0.26275)
                self.transfer_func.AddRGBPoint(0.50, 1.00000, 1.00000, 0.74902)
                self.transfer_func.AddRGBPoint(0.75, 0.40000, 0.76078, 0.64706)
                self.transfer_func.AddRGBPoint(1.00, 0.36863, 0.30980, 0.63529)
            else:
                self.transfer_func.AddRGBPoint(-1.00, 0.61961, 0.00392, 0.25882)
                self.transfer_func.AddRGBPoint(-0.50, 0.95686, 0.42745, 0.26275)
                self.transfer_func.AddRGBPoint(0.00, 1.00000, 1.00000, 0.74902)
                self.transfer_func.AddRGBPoint(0.50, 0.40000, 0.76078, 0.64706)
                self.transfer_func.AddRGBPoint(1.00, 0.36863, 0.30980, 0.63529)
        # -----------------------------------------------------------------------
        # High-jacking this scheme for Single colors
        if mapnum == -1:
            self.transfer_func.SetColorSpaceToRGB()
            if (viz_type == "M") or (viz_type == "E") or viz_type == "N":
                self.transfer_func.AddRGBPoint(0.0, 1.00000, 1.00000, 1.00000)
                self.transfer_func.AddRGBPoint(1.0, 1.00000, 1.00000, 1.00000)
            else:
                self.transfer_func.AddRGBPoint(-1.0, 1.00000, 1.00000, 1.00000)
                self.transfer_func.AddRGBPoint(1.0, 1.00000, 1.00000, 1.00000)
        # Construct the lut with the selected colomap
        # -----------------------------------------------------------------------
        for ii, ss in enumerate(
            [float(xx) / float(self.num_colors) for xx in range(self.num_colors)]
        ):
            cc = self.transfer_func.GetColor(ss)
            self.lut.SetTableValue(ii, cc[0], cc[1], cc[2], 1.0)
        self.lut.Build()
        return

    ##########################################################################
    # @brief Setup a single color colormap
    # @details Used since glyphs are colored after colormaps and not single colors
    # @author Anders Bergman
    ##########################################################################
    def set_RGBcolor(self, rgb, viz_type):
        """
        Set the RGB color for the lookup table (LUT) and color transfer function.

        Parameters:
        rgb (list): A list of three integers representing the RGB color values.
        window (object): The window object where the color will be applied.
        flag2D (bool): A flag indicating if the visualization is 2D.
        viz_type (str): The type of visualization (e.g., "M", "E", "N").
        renWin (object): The render window object.

        Returns:
        None

        Author:
        Anders Bergman
        """

        self.num_colors = 2
        red = rgb[0] / 255.0
        green = rgb[1] / 255.0
        blue = rgb[2] / 255.0

        self.lut.SetNumberOfTableValues(self.num_colors)
        self.transfer_func.SetColorSpaceToRGB()
        if (viz_type == "M") or (viz_type == "E") or viz_type == "N":
            self.transfer_func.AddRGBPoint(0.00, red, green, blue)
            self.transfer_func.AddRGBPoint(1.00, red, green, blue)
        else:
            self.transfer_func.AddRGBPoint(-1.0, red, green, blue)
            self.transfer_func.AddRGBPoint(1.0, red, green, blue)

        for ii, ss in enumerate(
            [float(xx) / float(self.num_colors) for xx in range(self.num_colors)]
        ):
            cc = self.transfer_func.GetColor(ss)
            self.lut.SetTableValue(ii, cc[0], cc[1], cc[2], 1.0)

        self.lut.Build()
        return

    ##########################################################################
    # @brief Set the background color given RGB values
    # @author Anders Bergman
    ##########################################################################
    def set_RGBbackground(self, rgb, ren):
        """
        Sets the background color of the renderer using an RGB tuple.
        """

        nrgb = [i / 255.0 for i in rgb]
        ren.SetBackground(nrgb)

        return

    ##########################################################################
    # @brief Update rgb-values for single color coloring
    # @author Anders Bergman
    ##########################################################################
    def set_singlecolor(self, window, value):
        """
        Set the single color for the RGB sliders and update the visualization.
        """
        self.num_colors = self.NUM_RGB_COLORS
        if window.bwSinglecolor:
            window.RGBRedColorSlider.setValue(value)
            window.RGBGreenColorSlider.setValue(value)
            window.RGBBlueColorSlider.setValue(value)

        rgb = [
            window.RGBRedColorSlider.value(),
            window.RGBGreenColorSlider.value(),
            window.RGBBlueColorSlider.value(),
        ]

        self.set_RGBcolor(
            rgb=rgb,
            viz_type=window.viz_type,
        )
        self.rgb = rgb
        return

    ##########################################################################
    # @brief Update rgb-values for the background
    # @author Anders Bergman
    ##########################################################################
    def set_background(self, window, value):
        """
        Sets the background color based on the provided value.
        """
        if window.bwBackground:
            window.RGBRedBackgroundSlider.setValue(value)
            window.RGBGreenBackgroundSlider.setValue(value)
            window.RGBBlueBackgroundSlider.setValue(value)

        rgb = [
            window.RGBRedBackgroundSlider.value(),
            window.RGBGreenBackgroundSlider.value(),
            window.RGBBlueBackgroundSlider.value(),
        ]

        self.set_RGBbackground(rgb=rgb, ren=window.ren)
        self.rgb_background = rgb

        return

    ##########################################################################
    # @brief Set the lookup table for the actors
    # @details Set the lookup table for the actors, it also allows for the change
    # of the scale type for the plotting, either linear or logarithmic scale.
    # @author Jonathan Chico
    ##########################################################################
    def set_lut_db(self, window, mapnum):
        """
        Sets the lookup table (LUT) for the visualization based on the provided colormap number.
        """
        self.num_colors = self.NUM_LUT_COLORS
        self.mapnum = mapnum
        colorSeries = vtkColorSeries()

        if mapnum <= 3:
            self.set_colormap_db(mapnum=mapnum, viz_type=window.viz_type)
            self.transfer_func = self.LUT_to_TransferFunction(self.lut)
        elif mapnum == 4:  # Spectrum
            colorSeries.SetColorScheme(vtkColorSeries.SPECTRUM)
        elif mapnum == 5:  # Warm
            colorSeries.SetColorScheme(vtkColorSeries.WARM)
        elif mapnum == 6:  # Cool
            colorSeries.SetColorScheme(vtkColorSeries.COOL)
        elif mapnum == 7:  # Blues
            colorSeries.SetColorScheme(vtkColorSeries.BLUES)
        elif mapnum == 8:  # Wildflower
            colorSeries.SetColorScheme(vtkColorSeries.WILD_FLOWER)
        elif mapnum == 9:  # Citrus
            colorSeries.SetColorScheme(vtkColorSeries.CITRUS)
        elif mapnum == 10:  # BREWER_DIVERGING_PURPLE_ORANGE_11
            colorSeries.SetColorScheme(vtkColorSeries.BREWER_DIVERGING_PURPLE_ORANGE_11)
        elif mapnum == 11:  # Citrus
            colorSeries.SetColorScheme(
                vtkColorSeries.BREWER_DIVERGING_BROWN_BLUE_GREEN_11
            )
        elif mapnum == 12:  # Citrus
            colorSeries.SetColorScheme(vtkColorSeries.BREWER_SEQUENTIAL_BLUE_GREEN_9)
        elif mapnum == 13:  # Citrus
            colorSeries.SetColorScheme(
                vtkColorSeries.BREWER_SEQUENTIAL_YELLOW_ORANGE_BROWN_9
            )
        elif mapnum == 14:  # Citrus
            colorSeries.SetColorScheme(vtkColorSeries.BREWER_SEQUENTIAL_BLUE_PURPLE_9)
        elif mapnum == 15:  # Citrus
            colorSeries.SetColorScheme(vtkColorSeries.BREWER_DIVERGING_SPECTRAL_11)
        elif mapnum == 16:  # Citrus
            colorSeries.SetColorScheme(vtkColorSeries.BREWER_QUALITATIVE_ACCENT)
        elif mapnum == 17:  # Citrus
            colorSeries.SetColorScheme(vtkColorSeries.BREWER_QUALITATIVE_DARK2)
        elif mapnum == 18:  # Citrus
            colorSeries.SetColorScheme(vtkColorSeries.BREWER_QUALITATIVE_PASTEL1)
        elif mapnum == 19:  # Citrus
            colorSeries.SetColorScheme(vtkColorSeries.BREWER_QUALITATIVE_PASTEL2)
        elif mapnum == 20:  # Citrus
            colorSeries.SetColorScheme(vtkColorSeries.BREWER_QUALITATIVE_SET1)
        elif mapnum == 21:  # Citrus
            colorSeries.SetColorScheme(vtkColorSeries.BREWER_QUALITATIVE_SET2)
        elif mapnum == 22:  # Citrus
            colorSeries.SetColorScheme(vtkColorSeries.BREWER_QUALITATIVE_SET3)

        if mapnum > 3:
            colorSeries.BuildLookupTable(self.lut, vtkColorSeries.ORDINAL)
            self.transfer_func = self.Series_to_TransferFunction(colorSeries)
        self.lut.Build()

        self.ColorActors(
            window=window,
            flag2D=window.ASDdata.flag2D,
            viz_type=window.viz_type,
            renWin=window.renWin,
        )
        return

    ##########################################################################
    # @brief Set the lookup table for the actors
    # @details Set the lookup table for the actors, it also allows for the change
    # of the scale type for the plotting, either linear or logarithmic scale.
    # @author Jonathan Chico
    ##########################################################################
    def set_lut_scale(self, window):
        """
        Sets the lookup table (LUT) scale based on the sender's state.
        """
        # self.ASDVizOpt.set_colormap(window=self,flag2D=self.ASDdata.flag2D,\
        # viz_type=self.viz_type,renWin=self.renWin)
        if window.sender() == window.LinearScale and window.LinearScale.isChecked():
            self.lut.SetScaleToLinear()
            self.lut_scale = "Linear"
        if window.sender() == window.LogScale and window.LogScale.isChecked():
            self.lut.SetScaleToLog10()
            self.lut_scale = "Logarithmic"
        return
