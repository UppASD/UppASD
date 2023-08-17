"""
Module pyasd


Defined at /home/andersb/CrossPlatform/UppASD/source/pyasd.f90 lines 41-248

"""
from __future__ import print_function, absolute_import, division
import _uppasd
import logging

_arrays = {}
_objs = {}

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
    _uppasd.setupall()

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

def relaxmontecarlo():
    """
    relaxmontecarlo()
    
    
    Defined at /home/andersb/CrossPlatform/UppASD/source/pyasd.f90 lines 111-113
    
    
    """
    _uppasd.relaxmontecarlo()

def relaxmetropolis():
    """
    relaxmetropolis()
    
    
    Defined at /home/andersb/CrossPlatform/UppASD/source/pyasd.f90 lines 115-119
    
    
    """
    _uppasd.relaxmetropolis()

def relaxheatbath():
    """
    relaxheatbath()
    
    
    Defined at /home/andersb/CrossPlatform/UppASD/source/pyasd.f90 lines 121-125
    
    
    """
    _uppasd.relaxheatbath()

def relaxmultiscale():
    """
    relaxmultiscale()
    
    
    Defined at /home/andersb/CrossPlatform/UppASD/source/pyasd.f90 lines 127-129
    
    
    """
    _uppasd.relaxmultiscale()

def relaxsldmontecarlo():
    """
    relaxsldmontecarlo()
    
    
    Defined at /home/andersb/CrossPlatform/UppASD/source/pyasd.f90 lines 131-133
    
    
    """
    _uppasd.relaxsldmontecarlo()

def relaxllg():
    """
    relaxllg()
    
    
    Defined at /home/andersb/CrossPlatform/UppASD/source/pyasd.f90 lines 135-137
    
    
    """
    _uppasd.relaxllg()

def relaxmd():
    """
    relaxmd()
    
    
    Defined at /home/andersb/CrossPlatform/UppASD/source/pyasd.f90 lines 139-141
    
    
    """
    _uppasd.relaxmd()

def relaxsldllg():
    """
    relaxsldllg()
    
    
    Defined at /home/andersb/CrossPlatform/UppASD/source/pyasd.f90 lines 143-145
    
    
    """
    _uppasd.relaxsldllg()

def relaxgneb():
    """
    relaxgneb()
    
    
    Defined at /home/andersb/CrossPlatform/UppASD/source/pyasd.f90 lines 180-182
    
    
    """
    _uppasd.relaxgneb()

def runmontecarlo():
    """
    runmontecarlo()
    
    
    Defined at /home/andersb/CrossPlatform/UppASD/source/pyasd.f90 lines 184-186
    
    
    """
    _uppasd.runmontecarlo()

def runmultiscale():
    """
    runmultiscale()
    
    
    Defined at /home/andersb/CrossPlatform/UppASD/source/pyasd.f90 lines 188-190
    
    
    """
    _uppasd.runmultiscale()

def runllglite():
    """
    runllglite()
    
    
    Defined at /home/andersb/CrossPlatform/UppASD/source/pyasd.f90 lines 192-194
    
    
    """
    _uppasd.runllglite()

def runllg():
    """
    runllg()
    
    
    Defined at /home/andersb/CrossPlatform/UppASD/source/pyasd.f90 lines 196-198
    
    
    """
    _uppasd.runllg()

def runllgcuda():
    """
    runllgcuda()
    
    
    Defined at /home/andersb/CrossPlatform/UppASD/source/pyasd.f90 lines 200-202
    
    
    """
    _uppasd.runllgcuda()

def runsldmontecarlo():
    """
    runsldmontecarlo()
    
    
    Defined at /home/andersb/CrossPlatform/UppASD/source/pyasd.f90 lines 204-206
    
    
    """
    _uppasd.runsldmontecarlo()

def runld():
    """
    runld()
    
    
    Defined at /home/andersb/CrossPlatform/UppASD/source/pyasd.f90 lines 208-210
    
    
    """
    _uppasd.runld()

def runsldllg():
    """
    runsldllg()
    
    
    Defined at /home/andersb/CrossPlatform/UppASD/source/pyasd.f90 lines 212-214
    
    
    """
    _uppasd.runsldllg()

def runsldllgimplicit():
    """
    runsldllgimplicit()
    
    
    Defined at /home/andersb/CrossPlatform/UppASD/source/pyasd.f90 lines 216-218
    
    
    """
    _uppasd.runsldllgimplicit()

def rungneb():
    """
    rungneb()
    
    
    Defined at /home/andersb/CrossPlatform/UppASD/source/pyasd.f90 lines 228-230
    
    
    """
    _uppasd.rungneb()

def totalenergy():
    """
    energy = totalenergy()
    
    
    Defined at /home/andersb/CrossPlatform/UppASD/source/pyasd.f90 lines 232-236
    
    
    Returns
    -------
    energy : float
    
    """
    energy = _uppasd.totalenergy()
    return energy


_array_initialisers = []
_dt_array_initialisers = []

try:
    for func in _array_initialisers:
        func()
except ValueError:
    logging.debug('unallocated array(s) detected on import of module "pyasd".')

for func in _dt_array_initialisers:
    func()
