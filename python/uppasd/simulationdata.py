"""
Module simulationdata


Defined at /home/andersb/CrossPlatform/UppASD/source/System/simulationdata.f90 \
    lines 2-13

"""
from __future__ import print_function, absolute_import, division
from uppasd import _uppasd
import f90wrap.runtime
import logging

_arrays = {}
_objs = {}

def get_lambda1():
    """
    Element lambda1 ftype=real(dblprec) pytype=float
    
    
    Defined at /home/andersb/CrossPlatform/UppASD/source/System/simulationdata.f90 \
        line 8
    
    """
    return _uppasd.f90wrap_simulationdata__get__lambda1()

def set_lambda1(lambda1):
    _uppasd.f90wrap_simulationdata__set__lambda1(lambda1)

def get_lambda2():
    """
    Element lambda2 ftype=real(dblprec) pytype=float
    
    
    Defined at /home/andersb/CrossPlatform/UppASD/source/System/simulationdata.f90 \
        line 9
    
    """
    return _uppasd.f90wrap_simulationdata__get__lambda2()

def set_lambda2(lambda2):
    _uppasd.f90wrap_simulationdata__set__lambda2(lambda2)

def get_rstep():
    """
    Element rstep ftype=integer  pytype=int
    
    
    Defined at /home/andersb/CrossPlatform/UppASD/source/System/simulationdata.f90 \
        line 10
    
    """
    return _uppasd.f90wrap_simulationdata__get__rstep()

def set_rstep(rstep):
    _uppasd.f90wrap_simulationdata__set__rstep(rstep)

def get_mstep():
    """
    Element mstep ftype=integer  pytype=int
    
    
    Defined at /home/andersb/CrossPlatform/UppASD/source/System/simulationdata.f90 \
        line 11
    
    """
    return _uppasd.f90wrap_simulationdata__get__mstep()

def set_mstep(mstep):
    _uppasd.f90wrap_simulationdata__set__mstep(mstep)

def get_bn():
    """
    Element bn ftype=real(dblprec) pytype=float
    
    
    Defined at /home/andersb/CrossPlatform/UppASD/source/System/simulationdata.f90 \
        line 12
    
    """
    return _uppasd.f90wrap_simulationdata__get__bn()

bn = get_bn()

def get_total_energy():
    """
    Element total_energy ftype=real(dblprec) pytype=float
    
    
    Defined at /home/andersb/CrossPlatform/UppASD/source/System/simulationdata.f90 \
        line 14
    
    """
    return _uppasd.f90wrap_simulationdata__get__total_energy()

def set_total_energy(total_energy):
    _uppasd.f90wrap_simulationdata__set__total_energy(total_energy)


_array_initialisers = []
_dt_array_initialisers = []

try:
    for func in _array_initialisers:
        func()
except ValueError:
    logging.debug('unallocated array(s) detected on import of module \
        "simulationdata".')

for func in _dt_array_initialisers:
    func()
