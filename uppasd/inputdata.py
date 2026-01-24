"""
InputData class for UppASD data
"""


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
        self.nstep = self.driver.get_nstep()
        return self.nstep

    def get_mcnstep(self):
        """
        Retrieves the value of mcnstep from the _uppasd module.

        Returns:
            int: The value of mcnstep.
        """
        self.mcnstep = self.driver.get_mcnstep()
        return self.mcnstep

    def get_temp(self):
        """
        Retrieves the value of temp from the _uppasd module.

        Returns:
            float: The value of temp.
        """
        self.temp = self.driver.get_temperature()
        return self.temp

    def get_delta_t(self):
        """
        Retrieves the value of delta_t from the _uppasd module.

        Returns:
            float: The value of delta_t.
        """
        self.delta_t = self.driver.get_delta_t()
        return self.delta_t

    def get_hfield(self):
        """
        Retrieves the value of hfield from the _uppasd module.

        Returns:
            list: The value of hfield.
        """
        self.hfield = self.driver.get_hfield()
        return self.hfield

    def get_iptemp(self):
        """
        Retrieves the value of iptemp from the _uppasd module.

        Returns:
            float: The value of iptemp.
        """
        self.iptemp = self.driver.get_iptemperature()
        return self.iptemp

    def get_iphfield(self):
        """
        Retrieves the value of iphfield from the _uppasd module.

        Returns:
            list: The value of iphfield.
        """
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
        self.driver.put_nstep(self.nstep)

    def update_mcnstep(self, value=None):
        """
        Sets the value of mcnstep in the _uppasd module using the class variable value.
        If value is given, updates the original value.

        Args:
            value (int, optional): The new value of mcnstep. Defaults to None.
        """
        if value is not None:
            self.mcnstep = value
        self.driver.put_mcnstep(self.mcnstep)

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
        self.driver.put_delta_t(self.delta_t)

    def update_hfield(self, value=None):
        """
        Sets the value of hfield in the _uppasd module using the class variable value.
        If value is given, updates the original value.

        Args:
            value (list, optional): The new value of hfield. Defaults to None.
        """
        if value is not None:
            self.hfield = value
        self.driver.put_hfield(self.hfield)

    def update_iptemp(self, value=None):
        """
        Sets the value of iptemp in the _uppasd module using the class variable value.
        If value is given, updates the original value.

        Args:
            value (float, optional): The new value of iptemp. Defaults to None.
        """
        if value is not None:
            self.iptemp = value
        self.driver.put_iptemperature(self.iptemp)

    def update_iphfield(self, value=None):
        """
        Sets the value of iphfield in the _uppasd module using the class variable value.
        If value is given, updates the original value.

        Args:
            value (list, optional): The new value of iphfield. Defaults to None.
        """
        if value is not None:
            self.iphfield = value
        self.driver.put_iphfield(self.iphfield)

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
