"""
Module pyasd


Defined at /home/andersb/CrossPlatform/UppASD/source/pyasd.f90 lines 41-248

"""

from __future__ import absolute_import, division, print_function

import logging

import _uppasd
import numpy as np

_arrays = {}
_objs = {}


# Driver routines below
def runuppasd():
    """
    runuppasd()


    Defined at /home/andersb/CrossPlatform/UppASD/source/pyasd.f90 lines 53-55


    """
    _uppasd.runuppasd()


def sanitycheck():
    """
    sanitycheck()


    Defined at /home/andersb/CrossPlatform/UppASD/source/pyasd.f90 lines 60-62


    """
    _uppasd.sanitycheck()


def numprocs():
    """
    nprocs = numprocs()


    Defined at /home/andersb/CrossPlatform/UppASD/source/pyasd.f90 lines 67-71


    Returns
    -------
    nprocs : int

    """
    nprocs = _uppasd.numprocs()
    return nprocs


def printlogo():
    """
    printlogo()


    Defined at /home/andersb/CrossPlatform/UppASD/source/pyasd.f90 lines 76-78


    """
    _uppasd.printlogo()


def setupall():
    """
    setupall()


    Defined at /home/andersb/CrossPlatform/UppASD/source/pyasd.f90 lines 83-85


    """
    natom, mensemble = _uppasd.setupall()
    return natom, mensemble


def initialphase():
    """
    initialphase()


    Defined at /home/andersb/CrossPlatform/UppASD/source/pyasd.f90 lines 90-92


    """
    _uppasd.initialphase()


def measure():
    """
    measure()


    Defined at /home/andersb/CrossPlatform/UppASD/source/pyasd.f90 lines 97-99


    """
    _uppasd.measure()


def cleanup():
    """
    cleanup()


    Defined at /home/andersb/CrossPlatform/UppASD/source/pyasd.f90 lines 104-106


    """
    _uppasd.cleanup()


def relaxmontecarlo(natom, mensemble):
    """
    relaxmontecarlo()


    Defined at /home/andersb/CrossPlatform/UppASD/source/pyasd.f90 lines 111-113


    """
    moments = _uppasd.relaxmontecarlo(natom, mensemble)
    return moments


def relaxmetropolis(natom, mensemble):
    """
    relaxmetropolis()


    Defined at /home/andersb/CrossPlatform/UppASD/source/pyasd.f90 lines 115-119


    """
    moments = _uppasd.relaxmetropolis(natom, mensemble)
    return moments


def relaxheatbath(natom, mensemble):
    """
    relaxheatbath()


    Defined at /home/andersb/CrossPlatform/UppASD/source/pyasd.f90 lines 121-125


    """
    moments = _uppasd.relaxheatbath(natom, mensemble)
    return moments

def relaxllg(natom, mensemble):
    """
    relaxllg()


    Defined at /home/andersb/CrossPlatform/UppASD/source/pyasd.f90 lines 135-137


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
    relax()


    Defined at /home/andersb/CrossPlatform/UppASD/source/pyasd.f90 lines 135-137


    """
    method = {}
    method["S"] = "LLG"
    method["M"] = "Metropolis"
    method["H"] = "Heat Bath"
    print(f"Performing relaxation using the {method[imode]} method for {instep} steps.")
    moments = _uppasd.relax(
        natom, mensemble, imode, instep, itemperature, itimestep, idamping
    )
    return moments

# Measurable rotuines below


def get_emom(natom, mensemble):
    """
    emom = get_emom(natom, mensemble)


    Defined at /home/andersb/CrossPlatform/UppASD/source/pyasd.f90 lines 238-242


    Parameters
    ----------
    natom : int
    mensemble : int

    Returns
    -------
    emom : float array

    """
    emom = _uppasd.get_emom(natom, mensemble)
    return emom


def put_emom(emom, natom, mensemble):
    """
    put_emom(emom, natom, mensemble)


    Defined at /home/andersb/CrossPlatform/UppASD/source/pyasd.f90 lines 244-248


    Parameters
    ----------
    emom : float array
    natom : int
    mensemble : int

    """
    _uppasd.put_emom(emom, natom, mensemble)


def get_beff(natom, mensemble):
    """
    emom = get_emom(natom, mensemble)


    Defined at /home/andersb/CrossPlatform/UppASD/source/pyasd.f90 lines 238-242


    Parameters
    ----------
    natom : int
    mensemble : int

    Returns
    -------
    emom : float array

    """
    beff = _uppasd.get_beff(natom, mensemble)
    return beff


def put_beff(beff, natom, mensemble):
    """
    put_emom(emom, natom, mensemble)


    Defined at /home/andersb/CrossPlatform/UppASD/source/pyasd.f90 lines 244-248


    Parameters
    ----------
    emom : float array
    natom : int
    mensemble : int

    """
    _uppasd.put_beff(beff, natom, mensemble)


def get_energy():
    """
    energy = totalenergy()


    Defined at /home/andersb/CrossPlatform/UppASD/source/pyasd.f90 lines 232-236


    Returns
    -------
    energy : float

    """
    energy = _uppasd.get_energy()
    return energy
