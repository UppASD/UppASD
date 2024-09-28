"""
ASDVTKCamera Module
This module provides the CameraManager class, which manages VTK camera settings,
including saving and loading configurations, resetting the camera,
updating camera values from a GUI, toggling parallel projection,
changing parallel projection parameters, and setting the camera's view up direction.

Classes:
    CameraManager: Manages VTK camera settings.

CameraManager Methods:
    __init__(self, camera):
        Initializes the CameraManager with a VTK camera.
    get_camera_settings(self):
        Retrieves the current camera settings.
    set_camera_settings(self, settings):
        Sets the camera to the stored settings.
    save_camera_settings(self, filename='camera_settings.json'):
        Saves the current camera settings to a JSON file.
    load_camera_settings(self, filename='camera_settings.json'):
        Loads camera settings from a JSON file and applies them.
    reset_camera(self, ren, renWin, current_Actors):
        Resets the camera to the initial position.
    Update_Camera(self, Window, ren, renWin):
        Updates the camera values for the current visualization from the GUI.
    toggle_projections(self, renWin, window, ren, checked):
        Toggles the parallel projection mode of the active camera.
    ChangeParallelProj(self, ren, renWin, line, slider, MainWindow):
        Adjusts the parallel projection scale based on input from a line edit or slider.
    set_Camera_viewUp(self, ren, renWin, rdir):
        Sets the camera's view up direction and renders the window.

Author:
    Anders Bergman, Jonathan Chico
"""
# pylint: disable=invalid-name, no-name-in-module, no-member

import yaml

from vtk import vtkPNGWriter, vtkPOVExporter, vtkWindowToImageFilter


class CameraManager:
    """
    Manages VTK camera settings, including saving and loading configurations.
    """

    def __init__(self, camera):
        self.camera = camera
        self.settings = {}

    def get_camera_settings(self):
        """
        Retrieve the current camera settings.

        Returns:
            dict: A dictionary containing the current camera settings, including:
                - 'position': The camera's position.
                - 'focal_point': The camera's focal point.
                - 'view_up': The camera's view up vector.
                - 'clipping_range': The camera's clipping range.
                - 'parallel_scale': The camera's parallel scale.
                - 'is_parallel': A boolean indicating whether the camera is in
                    parallel projection mode.
        """
        self.settings = {
            "position": self.camera.GetPosition(),
            "focal_point": self.camera.GetFocalPoint(),
            "view_up": self.camera.GetViewUp(),
            "clipping_range": self.camera.GetClippingRange(),
            "parallel_scale": self.camera.GetParallelScale(),
            "is_parallel": self.camera.GetParallelProjection(),
        }
        return self.settings

    def set_camera_settings(self, settings):
        """
        Set the camera to the stored settings.

        Args:
            settings (dict): A dictionary containing the camera settings to apply, including:
                - 'position': The camera's position.
                - 'focal_point': The camera's focal point.
                - 'view_up': The camera's view up vector.
                - 'clipping_range': The camera's clipping range.
                - 'parallel_scale': The camera's parallel scale.
                - 'is_parallel': A boolean indicating whether the camera is in
                    parallel projection mode.
        """
        self.camera.SetPosition(settings["position"])
        self.camera.SetFocalPoint(settings["focal_point"])
        self.camera.SetViewUp(settings["view_up"])
        self.camera.SetClippingRange(settings["clipping_range"])
        self.camera.SetParallelScale(settings["parallel_scale"])
        self.camera.SetParallelProjection(settings["is_parallel"])

    def save_camera_settings(self, filename="camera_settings.yaml"):
        """
        Save the current camera settings to a YAML file.

        Args:
            filename (str, optional): The name of the file to save the settings to.
            Defaults to 'camera_settings.yaml'.
        """
        self.get_camera_settings()
        with open(filename, "w", encoding="utf-8") as f:
            yaml.dump(self.settings, f, default_flow_style=False)
        print(f"Camera settings saved to {filename}.")

    def load_camera_settings(self, filename="camera_settings.yaml"):
        """
        Load camera settings from a YAML file and apply them.

        Args:
            filename (str, optional): The name of the file to load the settings from.
            Defaults to 'camera_settings.yaml'.
        """
        try:
            with open(filename, "r", encoding="utf-8") as f:
                self.settings = yaml.safe_load(f)
            self.set_camera_settings(self.settings)
            print(f"Camera settings loaded from {filename}.")
        except FileNotFoundError:
            print(f"File {filename} not found. Unable to load settings.")

    ##########################################################################
    # @brief Function to reset the camera to the initial positions
    # @author Jonathan Chico
    ##########################################################################
    def reset_camera(self, ren, renWin, current_Actors=None):
        """Function to reset the camera to the initial position.

        Args:
            ren: current VTK renderer.
            renWin: current VTK rendering window.
            current_Actors: current actors which are being visualized.

        Author
        ----------
        Jonathan Chico
        """
        # Defining the camera directions
        if current_Actors is None:
            ren.GetActiveCamera().SetFocalPoint(self.settings["focal_point"])
            ren.GetActiveCamera().SetPosition(self.settings["position"])
            ren.GetActiveCamera().SetViewUp(self.settings["view_up"])
        else:
            ren.GetActiveCamera().SetFocalPoint(
                current_Actors.xmid, current_Actors.ymid, current_Actors.zmid
            )
            ren.GetActiveCamera().SetPosition(
                current_Actors.xmid, current_Actors.ymid, current_Actors.height
            )
            ren.GetActiveCamera().SetViewUp(0, 1, 0)
            # ren.GetActiveCamera().Azimuth(0)
            # ren.GetActiveCamera().Elevation(0)
            # ren.GetActiveCamera().Yaw(0)
            # ren.GetActiveCamera().Roll(0)
            # ren.GetActiveCamera().Pitch(0)

        renWin.Render()
        return

    ##########################################################################
    # @brief Update the camera values for those defined by the user in the GUI.
    # @author Jonathan Chico
    ##########################################################################
    def Update_Camera(self, Window, ren, renWin):
        """Function to update the value of the camera for the current visualization from
        the values entered by the user in the GUI.

        Args:
            Window: QMainWindow where the visualizations are being carried out.
            ren: current VTK renderer.
            renWin: current VTK rendering window.

        Author
        ----------
        Jonathan Chico
        """
        camera_focal = [0] * 3
        camera_pos = [0] * 3
        camera_pos[0] = float(Window.CamPosX.text())
        camera_pos[1] = float(Window.CamPosY.text())
        camera_pos[2] = float(Window.CamPosZ.text())
        camera_focal[0] = float(Window.FocalPosX.text())
        camera_focal[1] = float(Window.FocalPosY.text())
        camera_focal[2] = float(Window.FocalPosZ.text())
        ren.GetActiveCamera().SetFocalPoint(camera_focal)
        ren.GetActiveCamera().SetPosition(camera_pos)
        ren.GetActiveCamera().Elevation(float(Window.CamElevationLineEdit.text()))
        ren.GetActiveCamera().Azimuth(float(Window.CamAzimuthLineEdit.text()))
        ren.GetActiveCamera().Pitch(float(Window.CamPitchLineEdit.text()))
        ren.GetActiveCamera().Roll(float(Window.CamRollLineEdit.text()))
        ren.GetActiveCamera().Yaw(float(Window.CamYawLineEdit.text()))

        renWin.Render()
        return

    ##########################################################################
    # @brief Toggle the parallel projection for the camera
    # @author Jonathan Chico
    ##########################################################################
    def toggle_projections(self, renWin, window, ren, checked):
        """
        Toggles the parallel projection mode of the active camera in the given renderer.

        Parameters:
        renWin (vtkRenderWindow): The render window to be updated.
        window (QMainWindow): The main window containing UI elements.
        ren (vtkRenderer): The renderer whose camera projection mode is to be toggled.
        checked (bool): If True, enables parallel projection; otherwise, disables it.

        Author: Jonathan Chico
        """
        if checked:
            ren.GetActiveCamera().ParallelProjectionOn()
            window.ParallelScaleLineEdit.setText(
                str(window.ParallelScaleSlider.value())
            )
            ren.GetActiveCamera().SetParallelScale(
                float(window.ParallelScaleLineEdit.text())
            )
            renWin.Render()
        else:
            ren.GetActiveCamera().ParallelProjectionOff()
            renWin.Render()

        return

    ##########################################################################
    # @brief Change the parameters of the parallel projection
    # @author Jonathan Chico
    ##########################################################################
    def ChangeParallelProj(self, ren, renWin, line, slider, MainWindow):
        """
        Adjusts the parallel projection scale based on input from a line edit or slider.
        """
        if line:
            MainWindow.ParallelScaleSlider.setValue(
                float(MainWindow.ParallelScaleLineEdit.text())
            )
            ren.GetActiveCamera().SetParallelScale(
                float(MainWindow.ParallelScaleLineEdit.text())
            )
            renWin.Render()
        if slider:
            ren.GetActiveCamera().SetParallelScale(
                MainWindow.ParallelScaleSlider.value()
            )
            MainWindow.ParallelScaleLineEdit.setText(
                str(MainWindow.ParallelScaleSlider.value())
            )
            renWin.Render()

        return

    ##########################################################################
    # Set the camera view up to be defined to the (1,0,0)
    ##########################################################################
    def set_Camera_viewUp(self, ren, renWin, rdir):
        """
        Sets the camera's view up direction and renders the window.
        """
        ren.GetActiveCamera().SetViewUp(rdir)
        renWin.Render()
        return

    ##########################################################################
    # Wrapper function to handle the camera functions
    ##########################################################################
    def camera_handler(self, window):
        """
        Handles various camera operations based on the sender of the signal.

        This method performs different camera-related actions such as resetting the camera,
        setting the camera view direction, updating the camera, and controlling the parallel scale.
        The specific action is determined by the sender of the signal.

        Actions:
        - Reset the camera to the original position if the sender is CamResetButton.
        - Set the camera view direction to X, Y, or Z axis if the sender is SetXView,
             SetYView, or SetZView respectively.
        - Update the camera if the sender is SetCamButton.
        - Change the parallel projection scale based on input from ParallelScaleLineEdit
            or ParallelScaleSlider.
        - Toggle parallel projections if the sender is ParallelProjectBox.
        """
        # -----------------------------------------------------------------------
        # Reset the camera to the original position
        # -----------------------------------------------------------------------
        if window.sender() == window.CamResetButton:
            if window.viz_type == "M":
                self.reset_camera(
                    ren=window.ren,
                    renWin=window.renWin,
                    current_Actors=window.MomActors,
                )
            elif window.viz_type == "N":
                self.reset_camera(
                    ren=window.ren,
                    renWin=window.renWin,
                    current_Actors=window.NeighActors,
                )
            elif window.viz_type == "E":
                self.reset_camera(
                    ren=window.ren,
                    renWin=window.renWin,
                    current_Actors=window.EneActors,
                )
        # -----------------------------------------------------------------------
        # Controlling what is up in the camera
        # -----------------------------------------------------------------------
        if window.sender() == window.SetXView:
            self.set_Camera_viewUp(ren=window.ren, renWin=window.renWin, rdir=(1, 0, 0))
        if window.sender() == window.SetYView:
            self.set_Camera_viewUp(ren=window.ren, renWin=window.renWin, rdir=(0, 1, 0))
        if window.sender() == window.SetZView:
            self.set_Camera_viewUp(ren=window.ren, renWin=window.renWin, rdir=(0, 0, 1))
        if window.sender() == window.SetCamButton:
            self.Update_Camera(Window=window, ren=window.ren, renWin=window.renWin)
        # -----------------------------------------------------------------------
        # Controlling the parallel scale
        # -----------------------------------------------------------------------
        if window.sender() == window.ParallelScaleLineEdit:
            line = True
            slider = False
            self.ChangeParallelProj(
                ren=window.ren,
                renWin=window.renWin,
                line=line,
                slider=slider,
                MainWindow=window,
            )
        if window.sender() == window.ParallelScaleSlider:
            line = False
            slider = True
            self.ChangeParallelProj(
                ren=window.ren,
                renWin=window.renWin,
                line=line,
                slider=slider,
                MainWindow=window,
            )
        if window.sender() == window.ParallelProjectBox:
            self.toggle_projections(
                renWin=window.renWin,
                window=window,
                ren=window.ren,
                checked=window.ParallelProjectBox.isChecked(),
            )
        if window.sender() == window.CamSaveButton:
            # Get and print the current camera settings
            # camera_settings = self.ASDCamera.get_camera_settings()
            self.save_camera_settings()

        if window.sender() == window.CamLoadButton:
            # Get and print the current camera settings
            self.load_camera_settings()
            # camera_settings = self.ASDCamera.get_camera_settings()
            window.ASDVizOpt.update_dock_info(
                current_actors=window.MomActors, window=window
            )

        return

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
