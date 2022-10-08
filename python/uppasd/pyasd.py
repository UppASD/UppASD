"""
Module pyasd


Defined at /home/andersb/CrossPlatform/UppASD/source/pyasd.f90 lines 41-248

"""
from __future__ import print_function, absolute_import, division
from uppasd import _uppasd
import f90wrap.runtime
import logging

_arrays = {}
_objs = {}

def runuppasd():
    """
    runuppasd()
    
    
    Defined at /home/andersb/CrossPlatform/UppASD/source/pyasd.f90 lines 53-55
    
    
    """
    _uppasd.f90wrap_runuppasd()

def sanitycheck():
    """
    sanitycheck()
    
    
    Defined at /home/andersb/CrossPlatform/UppASD/source/pyasd.f90 lines 60-62
    
    
    """
    _uppasd.f90wrap_sanitycheck()

def numprocs():
    """
    nprocs = numprocs()
    
    
    Defined at /home/andersb/CrossPlatform/UppASD/source/pyasd.f90 lines 67-71
    
    
    Returns
    -------
    nprocs : int
    
    """
    nprocs = _uppasd.f90wrap_numprocs()
    return nprocs

def printlogo():
    """
    printlogo()
    
    
    Defined at /home/andersb/CrossPlatform/UppASD/source/pyasd.f90 lines 76-78
    
    
    """
    _uppasd.f90wrap_printlogo()

def setupall():
    """
    setupall()
    
    
    Defined at /home/andersb/CrossPlatform/UppASD/source/pyasd.f90 lines 83-85
    
    
    """
    _uppasd.f90wrap_setupall()

def initialphase():
    """
    initialphase()
    
    
    Defined at /home/andersb/CrossPlatform/UppASD/source/pyasd.f90 lines 90-92
    
    
    """
    _uppasd.f90wrap_initialphase()

def measure():
    """
    measure()
    
    
    Defined at /home/andersb/CrossPlatform/UppASD/source/pyasd.f90 lines 97-99
    
    
    """
    _uppasd.f90wrap_measure()

def cleanup():
    """
    cleanup()
    
    
    Defined at /home/andersb/CrossPlatform/UppASD/source/pyasd.f90 lines 104-106
    
    
    """
    _uppasd.f90wrap_cleanup()

def relaxmontecarlo():
    """
    relaxmontecarlo()
    
    
    Defined at /home/andersb/CrossPlatform/UppASD/source/pyasd.f90 lines 111-113
    
    
    """
    _uppasd.f90wrap_relaxmontecarlo()

def relaxmetropolis():
    """
    relaxmetropolis()
    
    
    Defined at /home/andersb/CrossPlatform/UppASD/source/pyasd.f90 lines 115-119
    
    
    """
    _uppasd.f90wrap_relaxmetropolis()

def relaxheatbath():
    """
    relaxheatbath()
    
    
    Defined at /home/andersb/CrossPlatform/UppASD/source/pyasd.f90 lines 121-125
    
    
    """
    _uppasd.f90wrap_relaxheatbath()

def relaxmultiscale():
    """
    relaxmultiscale()
    
    
    Defined at /home/andersb/CrossPlatform/UppASD/source/pyasd.f90 lines 127-129
    
    
    """
    _uppasd.f90wrap_relaxmultiscale()

def relaxsldmontecarlo():
    """
    relaxsldmontecarlo()
    
    
    Defined at /home/andersb/CrossPlatform/UppASD/source/pyasd.f90 lines 131-133
    
    
    """
    _uppasd.f90wrap_relaxsldmontecarlo()

def relaxllg():
    """
    relaxllg()
    
    
    Defined at /home/andersb/CrossPlatform/UppASD/source/pyasd.f90 lines 135-137
    
    
    """
    _uppasd.f90wrap_relaxllg()

def relaxmd():
    """
    relaxmd()
    
    
    Defined at /home/andersb/CrossPlatform/UppASD/source/pyasd.f90 lines 139-141
    
    
    """
    _uppasd.f90wrap_relaxmd()

def relaxsldllg():
    """
    relaxsldllg()
    
    
    Defined at /home/andersb/CrossPlatform/UppASD/source/pyasd.f90 lines 143-145
    
    
    """
    _uppasd.f90wrap_relaxsldllg()

def relaxgneb():
    """
    relaxgneb()
    
    
    Defined at /home/andersb/CrossPlatform/UppASD/source/pyasd.f90 lines 180-182
    
    
    """
    _uppasd.f90wrap_relaxgneb()

def runmontecarlo():
    """
    runmontecarlo()
    
    
    Defined at /home/andersb/CrossPlatform/UppASD/source/pyasd.f90 lines 184-186
    
    
    """
    _uppasd.f90wrap_runmontecarlo()

def runmultiscale():
    """
    runmultiscale()
    
    
    Defined at /home/andersb/CrossPlatform/UppASD/source/pyasd.f90 lines 188-190
    
    
    """
    _uppasd.f90wrap_runmultiscale()

def runllglite():
    """
    runllglite()
    
    
    Defined at /home/andersb/CrossPlatform/UppASD/source/pyasd.f90 lines 192-194
    
    
    """
    _uppasd.f90wrap_runllglite()

def runllg():
    """
    runllg()
    
    
    Defined at /home/andersb/CrossPlatform/UppASD/source/pyasd.f90 lines 196-198
    
    
    """
    _uppasd.f90wrap_runllg()

def runllgcuda():
    """
    runllgcuda()
    
    
    Defined at /home/andersb/CrossPlatform/UppASD/source/pyasd.f90 lines 200-202
    
    
    """
    _uppasd.f90wrap_runllgcuda()

def runsldmontecarlo():
    """
    runsldmontecarlo()
    
    
    Defined at /home/andersb/CrossPlatform/UppASD/source/pyasd.f90 lines 204-206
    
    
    """
    _uppasd.f90wrap_runsldmontecarlo()

def runld():
    """
    runld()
    
    
    Defined at /home/andersb/CrossPlatform/UppASD/source/pyasd.f90 lines 208-210
    
    
    """
    _uppasd.f90wrap_runld()

def runsldllg():
    """
    runsldllg()
    
    
    Defined at /home/andersb/CrossPlatform/UppASD/source/pyasd.f90 lines 212-214
    
    
    """
    _uppasd.f90wrap_runsldllg()

def runsldllgimplicit():
    """
    runsldllgimplicit()
    
    
    Defined at /home/andersb/CrossPlatform/UppASD/source/pyasd.f90 lines 216-218
    
    
    """
    _uppasd.f90wrap_runsldllgimplicit()

def rungneb():
    """
    rungneb()
    
    
    Defined at /home/andersb/CrossPlatform/UppASD/source/pyasd.f90 lines 228-230
    
    
    """
    _uppasd.f90wrap_rungneb()

def totalenergy():
    """
    energy = totalenergy()
    
    
    Defined at /home/andersb/CrossPlatform/UppASD/source/pyasd.f90 lines 232-236
    
    
    Returns
    -------
    energy : float
    
    """
    energy = _uppasd.f90wrap_totalenergy()
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
