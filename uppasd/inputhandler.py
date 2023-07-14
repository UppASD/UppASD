"""
Module inputhandler


Defined at /home/andersb/CrossPlatform/UppASD/source/Input/inputhandler.f90 \
    lines 17-1286

"""
from __future__ import print_function, absolute_import, division
from uppasd import _uppasd
import f90wrap.runtime
import logging

_arrays = {}
_objs = {}

def read_parameters(ifile):
    """
    read_parameters(ifile)
    
    
    Defined at /home/andersb/CrossPlatform/UppASD/source/Input/inputhandler.f90 \
        lines 45-1265
    
    Parameters
    ----------
    ifile : int
    
    ------------------------------------------------------------------------
     START OF MISC VARIABLES
    ------------------------------------------------------------------------
    """
    _uppasd.f90wrap_read_parameters(ifile=ifile)

def change_constants():
    """
    change_constants()
    
    
    Defined at /home/andersb/CrossPlatform/UppASD/source/Input/inputhandler.f90 \
        lines 1270-1286
    
    
    """
    _uppasd.f90wrap_change_constants()

def get_sane_input():
    """
    Element sane_input ftype=logical pytype=bool
    
    
    Defined at /home/andersb/CrossPlatform/UppASD/source/Input/inputhandler.f90 line \
        23
    
    """
    return _uppasd.f90wrap_inputhandler__get__sane_input()

def set_sane_input(sane_input):
    _uppasd.f90wrap_inputhandler__set__sane_input(sane_input)


_array_initialisers = []
_dt_array_initialisers = []

try:
    for func in _array_initialisers:
        func()
except ValueError:
    logging.debug('unallocated array(s) detected on import of module \
        "inputhandler".')

for func in _dt_array_initialisers:
    func()
