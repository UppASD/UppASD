import numpy as np
from uppasd import pyasd as _asd

class simulator:
    """
    Class for atomistic spin dynamics simulations.
    """

    def __init__(self):
        """
        Initialize the simulation.

        Parameters:
        - moments: The initial moments of the spins.
        """
        self.natom, self.mensemble = _asd.setupall()
        self.moments = np.zeros((3, self.natom, self.mensemble))
        self.fields = np.zeros((3, self.natom, self.mensemble))
        self.energy = np.float64(0.0)

    def run_simulation(self):
        """
        Run the atomistic spin dynamics simulation.

        This method executes the steps required to run the atomistic spin dynamics simulation.
        It calls the necessary functions from the `asd` module to print the logo, initialize the phase,
        measure the simulation, and perform cleanup afterwards.
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
    
    def relax(self, mode='M'):
        """
        Relax the system using Monte Carlo simulations.
        """
        self.moments = _asd.relax(natom=self.natom, mensemble=self.mensemble, imode=mode)

    def calculate_energy(self):
        """
        Calculate the total energy of the system.

        Returns:
        - energy: The total energy of the system.
        """
        self.energy = _asd.totalenergy()


    def get_moments(self):
        """
        Update the moments of the spins based on the simulation dynamics.
        """
        self.moments = _asd.get_emom(natom=self.natom, mensemble=self.mensemble)

    def put_moments(self, moments):
        """
        Update the moments of the spins based on the simulation dynamics.
        """
        _asd.put_emom(moments=moments,natom=self.natom, mensemble=self.mensemble)
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
        _asd.get_beff(fields=fields, natom=self.natom, mensemble=self.mensemble)
        
    def evolve(self, evolution_type: str = 'initial'):
        """
        Evolve the system according to the `inpsd.dat` file.
        """
        if evolution_type == 'initial':
            _asd.initialphase()
        if evolution_type == 'measure':
            _asd.measure()

    # def print_simulation_info(self):
    #     """
    #     Print information about the simulation.
    #     """
    #     # TODO: Implement the simulation info printing logic here
    #     pass
