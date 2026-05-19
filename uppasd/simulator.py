import numpy as np
from uppasd import pyasd as _asd
from uppasd import inputdata as _inputdata


class Simulator:
    """
    Class for atomistic spin dynamics simulations.
    """

    def __init__(self):
        """
        Initialize the simulation.

        Parameters:
        - moments: The initial moments of the spins.
        """
        # self.natom, self.mensemble = _asd.setupall()
        self.natom = 1
        self.mensemble = 1
        self.moments = np.zeros((3, self.natom, self.mensemble))
        self.fields = np.zeros((3, self.natom, self.mensemble))
        self.coords = np.zeros((3, self.natom))
        self.energy = np.float64(0.0)

        self.inputdata = _inputdata.InputData(_asd)

    def run_uppasd(self):
        """
        Run the UppASD simulation.
        """
        _asd.runuppasd()

    def init_simulation(self):
        """
        Initialize the atomistic spin dynamics simulation.

        This method initializes the atomistic spin dynamics simulation
        by calling the necessary functions from the `asd` module.
        """
        # Initialize the simulation
        self.natom, self.mensemble = _asd.setupall()
        self.inputdata.get_all()
        self.get_moments()
        self.get_coords()

    def print_logo(self):
        """
        Print the UppASD logo.
        """
        _asd.printlogo()

    def run_simulation(self):
        """
        Run the atomistic spin dynamics simulation.

        This method executes the steps required to run the atomistic spin dynamics simulation.
        It calls the necessary functions from the `asd` module to print the logo,
        initialize the phase, measure the simulation, and perform cleanup afterwards.
        """
        # Print the logo
        _asd.printlogo()

        # Initialize the phase
        _asd.initialphase()

        # Measure the simulation
        _asd.measure()

        # Get the moments
        self.moments = _asd.get_emom(natom=self.natom, mensemble=self.mensemble)

        # Perform cleanup
        _asd.cleanup()

        return self.moments

    def relax(self, mode="M", temperature=0.0):
        """
        Relax the system using Monte Carlo simulations.
        """
        self.moments = _asd.relax(
            natom=self.natom,
            mensemble=self.mensemble,
            imode=mode,
            itemperature=temperature,
        )

    def calculate_energy(self):
        """
        Calculate the total energy of the system.

        Returns:
        - energy: The total energy of the system.
        """
        self.energy = _asd.get_energy()

    def get_coords(self):
        """
        Update the moments of the spins based on the simulation dynamics.
        """
        self.coords = _asd.get_coords(natom=self.natom)

    def get_moments(self):
        """
        Update the moments of the spins based on the simulation dynamics.
        """
        self.moments = np.copy(
            _asd.get_emom(natom=self.natom, mensemble=self.mensemble)
        )

    def put_moments(self, moments):
        """
        Update the moments of the spins based on the simulation dynamics.
        """
        _asd.put_emom(emom=moments, natom=self.natom, mensemble=self.mensemble)
        return

    def get_fields(self):
        """
        Apply the external magnetic fields to the spins.
        """
        self.fields = _asd.get_beff(natom=self.natom, mensemble=self.mensemble)

    def put_fields(self, fields):
        """
        Apply the external magnetic fields to the spins.
        """
        _asd.put_beff(beff=fields, natom=self.natom, mensemble=self.mensemble)

    def evolve(self, evolution_type: str = "initial"):
        """
        Evolve the system according to the `inpsd.dat` file.
        """
        if evolution_type == "initial":
            _asd.initialphase()
        if evolution_type == "measure":
            _asd.measure()
