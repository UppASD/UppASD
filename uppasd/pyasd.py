"""
Module pyasd

This module provides a Python interface to the UppASD library, which is written in Fortran.
It contains functions for running the UppASD simulation, performing relaxation methods,
and retrieving measured quantities.

The module defines the following functions:
   -  runuppasd: Run the UppASD simulation.
   -  sanitycheck: Perform a sanity check.
   -  numprocs: Get the number of processors.
   -  printlogo: Print the UppASD logo.
   -  setupall: Setup the UppASD simulation.
   -  initialphase: Perform the initial phase of the UppASD simulation.
   -  measure: Perform the measurement step of the UppASD simulation.
   -  cleanup: Clean up the UppASD simulation.
   -  relaxmontecarlo: Perform the Monte Carlo relaxation method.
   -  relaxmetropolis: Perform the Metropolis relaxation method.
   -  relaxheatbath: Perform the Heat Bath relaxation method.
   -  relaxllg: Perform the LLG relaxation method.
   -  relax: Perform the relaxation method.

"""

from __future__ import absolute_import, division, print_function

import _uppasd


# Driver routines below
def runuppasd():
    """
    Run the UppASD simulation.

    This function calls the underlying Fortran routine to run the UppASD simulation.

    """
    _uppasd.runuppasd()


def sanitycheck():
    """
    Perform a sanity check.

    This function calls the underlying Fortran routine to perform a sanity check.

    """
    _uppasd.sanitycheck()


def numprocs():
    """
    Get the number of processors.

    This function calls the underlying Fortran routine to get the number of processors.


    Returns
    -------
    nprocs : int
        The number of processors.
    """
    nprocs = _uppasd.numprocs()
    return nprocs


def printlogo():
    """
    Print the UppASD logo.

    This function calls the underlying Fortran routine to print the UppASD logo.

    """
    _uppasd.printlogo()


def setupall():
    """
    Setup the UppASD simulation.

    This function calls the underlying Fortran routine to setup the UppASD simulation.


    Returns
    -------
    natom : int
        The number of atoms.
    mensemble : int
        The number of ensembles.
    """
    natom, mensemble = _uppasd.setupall()
    return natom, mensemble


def initialphase():
    """
    Perform the initial phase of the UppASD simulation.

    This function calls the underlying Fortran routine to perform the initial phase 
    of the UppASD simulation.

    """
    _uppasd.initialphase()


def measure():
    """
    Perform the measurement step of the UppASD simulation.

    This function calls the underlying Fortran routine to perform the measurement step 
    of the UppASD simulation.

    """
    _uppasd.measure()


def cleanup():
    """
    Clean up the UppASD simulation.

    This function calls the underlying Fortran routine to clean up the UppASD simulation.

    """
    _uppasd.cleanup()


def relaxmontecarlo(natom, mensemble):
    """
    Perform the Monte Carlo relaxation method.

    This function calls the underlying Fortran routine to perform the Monte Carlo relaxation method.


    Parameters
    ----------
    natom : int
        The number of atoms.
    mensemble : int
        The number of ensembles.

    Returns
    -------
    moments : float array
        The moments of the relaxation.
    """
    moments = _uppasd.relaxmontecarlo(natom, mensemble)
    return moments


def relaxmetropolis(natom, mensemble):
    """
    Perform the Metropolis relaxation method.

    This function calls the underlying Fortran routine to perform the Metropolis relaxation method.


    Parameters
    ----------
    natom : int
        The number of atoms.
    mensemble : int
        The number of ensembles.

    Returns
    -------
    moments : float array
        The moments of the relaxation.
    """
    moments = _uppasd.relaxmetropolis(natom, mensemble)
    return moments


def relaxheatbath(natom, mensemble):
    """
    Perform the Heat Bath relaxation method.

    This function calls the underlying Fortran routine to perform the Heat Bath relaxation method.


    Parameters
    ----------
    natom : int
        The number of atoms.
    mensemble : int
        The number of ensembles.

    Returns
    -------
    moments : float array
        The moments of the relaxation.
    """
    moments = _uppasd.relaxheatbath(natom, mensemble)
    return moments


def relaxllg(natom, mensemble):
    """
    Perform the LLG relaxation method.

    This function calls the underlying Fortran routine to perform the LLG relaxation method.


    Parameters
    ----------
    natom : int
        The number of atoms.
    mensemble : int
        The number of ensembles.

    Returns
    -------
    moments : float array
        The moments of the relaxation.
    """
    moments = _uppasd.relaxllg(natom, mensemble)
    return moments


def relax(
    natom,
    mensemble,
    imode="S",
    instep=10,
    itemperature=0.0,
    itimestep=1.0e-16,
    idamping=0.5,
):
    """
    Perform the relaxation method.

    This function calls the underlying Fortran routine to perform the relaxation method.


    Parameters
    ----------
    natom : int
        The number of atoms.
    mensemble : int
        The number of ensembles.
    imode : str, optional
        The relaxation method. Default is "S" (LLG).
    instep : int, optional
        The number of relaxation steps. Default is 10.
    itemperature : float, optional
        The temperature. Default is 0.0.
    itimestep : float, optional
        The time step. Default is 1.0e-16.
    idamping : float, optional
        The damping factor. Default is 0.5.

    Returns
    -------
    moments : float array
        The moments of the relaxation.
    """
    method = {}
    method["S"] = "LLG"
    method["M"] = "Metropolis"
    method["H"] = "Heat Bath"
    print(f"Performing relaxation using the {method[imode]} method for {instep} steps.")
    print("Temperature:", itemperature)
    moments = _uppasd.relax(
        natom, mensemble, imode, instep, itemperature, itimestep, idamping
    )
    return moments

# Measurable routines below


def get_emom(natom, mensemble):
    """
    Get the effective magnetic moment.

    This function calls the underlying Fortran routine to get the effective magnetic moment.


    Parameters
    ----------
    natom : int
        The number of atoms.
    mensemble : int
        The number of ensembles.

    Returns
    -------
    emom : float array
        The effective magnetic moment.
    """
    emom = _uppasd.get_emom(natom, mensemble)
    return emom


def put_emom(emom, natom, mensemble):
    """
    Put the effective magnetic moment.

    This function calls the underlying Fortran routine to put the effective magnetic moment.


    Parameters
    ----------
    emom : float array
        The effective magnetic moment.
    natom : int
        The number of atoms.
    mensemble : int
        The number of ensembles.
    """
    _uppasd.put_emom(emom, natom, mensemble)


def get_beff(natom, mensemble):
    """
    Get the effective magnetic field.

    This function calls the underlying Fortran routine to get the effective magnetic field.


    Parameters
    ----------
    natom : int
        The number of atoms.
    mensemble : int
        The number of ensembles.

    Returns
    -------
    beff : float array
        The effective magnetic field.
    """
    beff = _uppasd.get_beff(natom, mensemble)
    return beff


def put_beff(beff, natom, mensemble):
    """
    Put the effective magnetic field.

    This function calls the underlying Fortran routine to put the effective magnetic field.


    Parameters
    ----------
    beff : float array
        The effective magnetic field.
    natom : int
        The number of atoms.
    mensemble : int
        The number of ensembles.
    """
    _uppasd.put_beff(beff, natom, mensemble)


def get_energy():
    """
    Get the total energy.

    This function calls the underlying Fortran routine to get the total energy.


    Returns
    -------
    energy : float
        The total energy.
    """
    energy = _uppasd.get_energy()
    return energy


def get_nstep():
    """
    Get the number of steps.

    This function calls the underlying Fortran routine to get the number of steps.


    Returns
    -------
    nstep : int
        The number of steps.
    """
    nstep = _uppasd.get_nstep()
    return nstep

def get_hfield():
    """
    Get the magnetic field.

    This function calls the underlying Fortran routine to get the magnetic field.


    Returns
    -------
    hfield : float
        The magnetic field.
    """
    hfield = _uppasd.get_hfield()
    return hfield

def put_hfield(hfield):
    """
    Put the magnetic field.

    This function calls the underlying Fortran routine to put the magnetic field.


    Parameters
    ----------
    hfield : float
        The magnetic field.
    """
    _uppasd.put_hfield(hfield)

def get_iphfield():
    """
    Get the initial magnetic field.

    This function calls the underlying Fortran routine to get the initial magnetic field.


    Returns
    -------
    iphfield : float
        The initial magnetic field.
    """
    iphfield = _uppasd.get_iphfield()
    return iphfield

def put_iphfield(iphfield):
    """
    Put the initial magnetic field.

    This function calls the underlying Fortran routine to put the initial magnetic field.


    Parameters
    ----------
    iphfield : float
        The initial magnetic field.
    """
    _uppasd.put_iphfield(iphfield)

def get_mcnstep():
    """
    Get the Monte Carlo number of steps.

    This function calls the underlying Fortran routine to get the Monte Carlo number of steps.


    Returns
    -------
    mcnstep : int
        The Monte Carlo number of steps.
    """
    mcnstep = _uppasd.get_mcnstep()
    return mcnstep

def get_temperature():
    """
    Get the temperature.

    This function calls the underlying Fortran routine to get the temperature.


    Returns
    -------
    temperature : float
        The temperature.
    """
    temperature = _uppasd.get_temperature()
    return temperature

def put_temperature(temperature):
    """
    Put the temperature.

    This function calls the underlying Fortran routine to put the temperature.


    Parameters
    ----------
    temperature : float
        The temperature.
    """
    _uppasd.put_temperature(temperature)

def get_iptemperature():
    """
    Get the initial temperature.

    This function calls the underlying Fortran routine to get the initial temperature.


    Returns
    -------
    iptemperature : float
        The initial temperature.
    """
    iptemperature = _uppasd.get_iptemperature()
    return iptemperature

def put_iptemperature(temperature):
    """
    Put the initial temperature.

    This function calls the underlying Fortran routine to set the initial temperature.


    Returns
    -------
    iptemperature : float
        The initial temperature.
    """
    _uppasd.put_iptemperature(temperature)

def get_delta_t():
    """
    Get the time step.

    This function calls the underlying Fortran routine to get the time step.


    Returns
    -------
    timestep : float
        The time step.
    """
    timestep = _uppasd.get_delta_t()
    return timestep
