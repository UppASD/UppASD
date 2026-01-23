"""ASDUICLIHelper Module
This module contains the ASDUICLIHelper class, which is a helper class for initializing and
managing the state of the ASD GUI application based on command-line interface (CLI) arguments.
Author:
    Anders Bergman
"""
# pylint: disable=invalid-name, no-name-in-module, no-member

import os
import sys
from ASD_GUI.UI import ASDUIActorHelper


class ASDUICLIHelper:
    """
    ASDUICLIHelper is a helper class for initializing and managing the state of the ASD GUI
    application based on command-line interface (CLI) arguments.

    Attributes:
        viz_type (str): The type of visualization to be used. It can be "M" for magnetization,
                        "N" for structure, or "E" for energy.
        file_list (dict): A dictionary containing file paths for various data files required
                          by the application.
        file_exist (dict): A dictionary indicating whether each file in `file_list` exists.

    Methods:
        __init__(self, window, args):
            Initializes the ASDUICLIHelper instance with the given window and CLI arguments.
        InitActorsFromCLI(self, window):
            Initializes actors in the given window based on the CLI settings.
        init_settings_from_cli(self, window):
            Initialize settings from the command line interface.
    """

    def __init__(self, window, args):
        """
        Initializes the ASDUICLIHelper class with the given window and arguments.

        Args:
            window: The window object that contains ASDdata attributes.
            args: An object containing file paths for various data files.

        Attributes:
            viz_type (str): The type of visualization to be used.
            file_list (dict): A dictionary mapping file types to their respective paths.
            file_exist (dict): A dictionary indicating the existence of each file.

        The method checks the existence of each file in the file_list and updates the
        file_exist dictionary accordingly. It also sets the corresponding attributes in
        the window's ASDdata based on the available files.
        """
        self.viz_type = None

        self.file_list = {
            "coordfile": args.coords,
            "momentfile": args.moments,
            "structfile": args.neighbours,
            "dmfile": args.dmi_neighbours,
            "enefile": args.energies,
            "settings_file": args.settings,
        }

        self.file_exist = {}
        self.has_errors = False

        for key, file_path in self.file_list.items():
            if file_path:
                if os.path.isfile(file_path):
                    self.file_exist[key] = True
                else:
                    self.file_exist[key] = False
                    self.has_errors = True
                    # Show error for missing file
                    error_msg = f"Error: {self._get_file_description(key)} '{file_path}' not found!"
                    print(error_msg)
                    sys.exit(1)
            else:
                self.file_exist[key] = False

        if self.has_errors:
            # Don't proceed with setting up files if there are errors
            return

        if self.file_exist["coordfile"]:
            window.ASDdata.posfiles = self.file_list["coordfile"]
        if self.file_exist["momentfile"]:
            window.ASDdata.magnetization = self.file_list["momentfile"]
            self.viz_type = "M"
        if self.file_exist["structfile"]:
            window.ASDdata.structfiles = self.file_list["structfile"]
            if window.viz_type is None:
                self.viz_type = "N"
        if self.file_exist["dmfile"]:
            window.ASDdata.dmdatafiles = self.file_list["dmfile"]
            if window.viz_type is None:
                self.viz_type = "N"
        if self.file_exist["enefile"]:
            window.ASDdata.enefiles = self.file_list["enefile"]
            if window.viz_type is None:
                self.viz_type = "E"
        if self.file_exist["settings_file"]:
            window.settings_file = self.file_list["settings_file"]

    def _get_file_description(self, key):
        """Get a human-readable description for a file type."""
        descriptions = {
            "coordfile": "Coordinate",
            "momentfile": "Magnetic moment",
            "structfile": "Structure",
            "dmfile": "DM data",
            "enefile": "Energy",
            "settings_file": "Settings",
        }
        return descriptions.get(key, "File")

    def InitActorsFromCLI(self, window):
        """
        Initializes actors in the given window based on the command-line interface (CLI)
        settings.

        Args:
            window: The window object where actors will be initialized. The window object
                should have a `viz_type` attribute.

        Prints:
            A message indicating that the `InitActorsFromCLI` method has been called.

        Side Effects:
            Sets the `viz_type` attribute of the window object.
            Calls `ASDUIActorHelper.AddActors` with the window and its `viz_type` if `viz_type
            is set.
        """
        if self.has_errors:
            # Don't initialize actors if there were file errors
            return
            
        window.viz_type = self.viz_type
        if self.viz_type:
            ASDUIActorHelper.AddActors(window, window.viz_type)

    def InitSettingsFromCLI(self, window):
        """
        Initialize settings from the command line interface.
        """
        if self.has_errors:
            # Don't load settings if there were file errors
            return
            
        if self.file_exist["settings_file"]:
            window.LoadSettings()
