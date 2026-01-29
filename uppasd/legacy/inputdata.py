"""
InputData class for UppASD data
"""
import logging
import numpy as np

logger = logging.getLogger(__name__)


class InputData:
    """
    InpitData class for UppASD data
    """

    def __init__(self, asd):
        """
        Initializes an instance of the InputData class.
        """
        self.driver = asd
        self.nstep = None
        self.mcnstep = None
        self.temp = None
        self.delta_t = None
        self.hfield = None
        self.iptemp = None
        self.iphfield = None
        self.ipmode = None

    def get_nstep(self):
        """
        Retrieves the value of nstep from the _uppasd module.

        Returns:
            int: The value of nstep.
        """
        if hasattr(self.driver, "get_nstep"):
            self.nstep = self.driver.get_nstep()
        return self.nstep

    def get_mcnstep(self):
        """
        Retrieves the value of mcnstep from the _uppasd module.

        Returns:
            int: The value of mcnstep.
        """
        if hasattr(self.driver, "get_mcnstep"):
            self.mcnstep = self.driver.get_mcnstep()
        return self.mcnstep

    def get_temp(self):
        """
        Retrieves the value of temp from the _uppasd module.

        Returns:
            float: The value of temp.
        """
        if hasattr(self.driver, "get_temperature"):
            self.temp = self.driver.get_temperature()
        return self.temp

    def get_delta_t(self):
        """
        Retrieves the value of delta_t from the _uppasd module.

        Returns:
            float: The value of delta_t.
        """
        if hasattr(self.driver, "get_delta_t"):
            self.delta_t = self.driver.get_delta_t()
        return self.delta_t

    def get_hfield(self):
        """
        Retrieves the value of hfield from the _uppasd module.

        Returns:
            list: The value of hfield.
        """
        if hasattr(self.driver, "get_hfield"):
            self.hfield = self.driver.get_hfield()
        return self.hfield

    def get_iptemp(self):
        """
        Retrieves the value of iptemp from the _uppasd module.

        Returns:
            float: The value of iptemp.
        """
        if hasattr(self.driver, "get_iptemperature"):
            self.iptemp = self.driver.get_iptemperature()
        return self.iptemp

    def get_iphfield(self):
        """
        Retrieves the value of iphfield from the _uppasd module.

        Returns:
            list: The value of iphfield.
        """
        if hasattr(self.driver, "get_iphfield"):
            self.iphfield = self.driver.get_iphfield()
        return self.iphfield

    def get_ipmode(self):
        """
        Retrieves the value of ipmode from the _uppasd module.

        Returns:
            str: The value of ipmode.
        """
        if hasattr(self.driver, "get_ipmode"):
            self.ipmode = self.driver.get_ipmode()
        return self.ipmode

    def update_nstep(self, value=None):
        """
        Sets the value of nstep in the _uppasd module using the class variable value.
        If value is given, updates the original value.

        Args:
            value (int, optional): The new value of nstep. Defaults to None.
        """
        if value is not None:
            self.nstep = value
        if hasattr(self.driver, "set_nstep"):
            try:
                self.driver.set_nstep(self.nstep)
            except Exception as e:
                logger.warning(f"Failed to update nstep: {e}")
        elif hasattr(self.driver, "put_nstep"):
            try:
                self.driver.put_nstep(self.nstep)
            except Exception as e:
                logger.warning(f"Failed to update nstep: {e}")

    def update_mcnstep(self, value=None):
        """
        Sets the value of mcnstep in the _uppasd module using the class variable value.
        If value is given, updates the original value.

        Args:
            value (int, optional): The new value of mcnstep. Defaults to None.
        """
        if value is not None:
            self.mcnstep = value
        if hasattr(self.driver, "set_mcnstep"):
            try:
                self.driver.set_mcnstep(self.mcnstep)
            except Exception as e:
                logger.warning(f"Failed to update mcnstep: {e}")
        elif hasattr(self.driver, "put_mcnstep"):
            try:
                self.driver.put_mcnstep(self.mcnstep)
            except Exception as e:
                logger.warning(f"Failed to update mcnstep: {e}")

    def update_temp(self, value=None):
        """
        Sets the value of temp in the _uppasd module using the class variable value.
        If value is given, updates the original value.

        Args:
            value (float, optional): The new value of temp. Defaults to None.
        """
        if value is not None:
            self.temp = value
        self.driver.put_temperature(self.temp)

    def update_delta_t(self, value=None):
        """
        Sets the value of delta_t in the _uppasd module using the class variable value.
        If value is given, updates the original value.

        Args:
            value (float, optional): The new value of delta_t. Defaults to None.
        """
        if value is not None:
            self.delta_t = value
        if hasattr(self.driver, "set_timestep"):
            try:
                self.driver.set_timestep(self.delta_t)
            except Exception as e:
                logger.warning(f"Failed to update delta_t: {e}")
        elif hasattr(self.driver, "set_delta_t"):
            try:
                self.driver.set_delta_t(self.delta_t)
            except Exception as e:
                logger.warning(f"Failed to update delta_t: {e}")
        elif hasattr(self.driver, "put_delta_t"):
            try:
                self.driver.put_delta_t(self.delta_t)
            except Exception as e:
                logger.warning(f"Failed to update delta_t: {e}")

    def update_hfield(self, value=None):
        """
        Sets the value of hfield in the _uppasd module using the class variable value.
        If value is given, updates the original value.

        Args:
            value (list, optional): The new value of hfield. Defaults to None.
        """
        if value is not None:
            self.hfield = value
        if hasattr(self.driver, "set_field"):
            try:
                # set_field expects shape (3, natom, mensemble) or (3,) for uniform field
                field = np.array(self.hfield, dtype=np.float64)
                if field.ndim == 1 and field.shape[0] == 3:
                    # Uniform field: expand to (3, natom, mensemble)
                    natom = self.driver.get_natom() if hasattr(self.driver, "get_natom") else 1
                    mensemble = self.driver.get_mensemble() if hasattr(self.driver, "get_mensemble") else 1
                    field = np.tile(field[:, np.newaxis, np.newaxis], (1, natom, mensemble))
                self.driver.set_field(field, natom, mensemble)
            except Exception as e:
                logger.warning(f"Failed to update hfield via set_field: {e}")
        elif hasattr(self.driver, "put_hfield"):
            try:
                self.driver.put_hfield(np.array(self.hfield, dtype=np.float64))
            except Exception as e:
                logger.warning(f"Failed to update hfield via put_hfield: {e}")

    def update_iptemp(self, value=None):
        """
        Sets the value of iptemp in the _uppasd module using the class variable value.
        If value is given, updates the original value.

        Args:
            value (float, optional): The new value of iptemp. Defaults to None.
        """
        if value is not None:
            self.iptemp = value
        if hasattr(self.driver, "set_iptemperature"):
            try:
                self.driver.set_iptemperature(self.iptemp)
            except Exception as e:
                logger.warning(f"Failed to update iptemp: {e}")
        elif hasattr(self.driver, "put_iptemperature"):
            try:
                self.driver.put_iptemperature(self.iptemp)
            except Exception as e:
                logger.warning(f"Failed to update iptemp: {e}")

    def update_iphfield(self, value=None):
        """
        Sets the value of iphfield in the _uppasd module using the class variable value.
        If value is given, updates the original value.

        Args:
            value (list, optional): The new value of iphfield. Defaults to None.
        """
        if value is not None:
            self.iphfield = value
        if hasattr(self.driver, "put_iphfield"):
            try:
                self.driver.put_iphfield(np.array(self.iphfield, dtype=np.float64))
            except Exception as e:
                logger.warning(f"Failed to update iphfield: {e}")

    def update_ipmode(self, value=None):
        """
        Sets the value of ipmode in the _uppasd module using the class variable value.
        If value is given, updates the original value.

        Args:
            value (str, optional): The new value of ipmode. Defaults to None.
        """
        if value is not None:
            self.ipmode = value
        if hasattr(self.driver, "put_ipmode"):
            self.driver.put_ipmode(self.ipmode)

    def get_all(self):
        """
        Calls all the get_* functions to update the class variables.
        """
        self.get_nstep()
        self.get_mcnstep()
        self.get_temp()
        self.get_delta_t()
        self.get_hfield()
        self.get_iptemp()
        self.get_iphfield()
        self.get_ipmode()
