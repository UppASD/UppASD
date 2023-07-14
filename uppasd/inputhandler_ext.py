"""
Module inputhandler_ext


Defined at /home/andersb/CrossPlatform/UppASD/source/Input/inputhandler_ext.f90 \
    lines 17-2156

"""
from __future__ import print_function, absolute_import, division
from uppasd import _uppasd
import f90wrap.runtime
import logging

_arrays = {}
_objs = {}

def read_positions():
    """
    read_positions()
    
    
    Defined at /home/andersb/CrossPlatform/UppASD/source/Input/inputhandler_ext.f90 \
        lines 51-114
    
    
    """
    _uppasd.f90wrap_read_positions()

def read_positions_alloy():
    """
    read_positions_alloy()
    
    
    Defined at /home/andersb/CrossPlatform/UppASD/source/Input/inputhandler_ext.f90 \
        lines 123-193
    
    
    """
    _uppasd.f90wrap_read_positions_alloy()

def read_moments(landeg_global):
    """
    read_moments(landeg_global)
    
    
    Defined at /home/andersb/CrossPlatform/UppASD/source/Input/inputhandler_ext.f90 \
        lines 204-287
    
    Parameters
    ----------
    landeg_global : float
    
    """
    _uppasd.f90wrap_read_moments(landeg_global=landeg_global)

def read_fixed_moments(landeg_global):
    """
    read_fixed_moments(landeg_global)
    
    
    Defined at /home/andersb/CrossPlatform/UppASD/source/Input/inputhandler_ext.f90 \
        lines 298-383
    
    Parameters
    ----------
    landeg_global : float
    
    """
    _uppasd.f90wrap_read_fixed_moments(landeg_global=landeg_global)

def read_exchange(self):
    """
    read_exchange(self)
    
    
    Defined at /home/andersb/CrossPlatform/UppASD/source/Input/inputhandler_ext.f90 \
        lines 399-554
    
    Parameters
    ----------
    ham_inp : Ham_Inp_T
    
    """
    _uppasd.f90wrap_read_exchange(ham_inp=self._handle)

def read_exchange_tensor():
    """
    read_exchange_tensor()
    
    
    Defined at /home/andersb/CrossPlatform/UppASD/source/Input/inputhandler_ext.f90 \
        lines 564-566
    
    
    """
    _uppasd.f90wrap_read_exchange_tensor()

def read_exchange_build_tensor():
    """
    read_exchange_build_tensor()
    
    
    Defined at /home/andersb/CrossPlatform/UppASD/source/Input/inputhandler_ext.f90 \
        lines 576-578
    
    
    """
    _uppasd.f90wrap_read_exchange_build_tensor()

def read_exchange_getmaxnoshells(filename):
    """
    no_shells, flines = read_exchange_getmaxnoshells(filename)
    
    
    Defined at /home/andersb/CrossPlatform/UppASD/source/Input/inputhandler_ext.f90 \
        lines 812-864
    
    Parameters
    ----------
    filename : str
    
    Returns
    -------
    no_shells : int
    flines : int
    
    """
    no_shells, flines = \
        _uppasd.f90wrap_read_exchange_getmaxnoshells(filename=filename)
    return no_shells, flines

def read_exchange_getneighvec(r_red, r_tmp, isite, jsite):
    """
    read_exchange_getneighvec(r_red, r_tmp, isite, jsite)
    
    
    Defined at /home/andersb/CrossPlatform/UppASD/source/Input/inputhandler_ext.f90 \
        lines 873-903
    
    Parameters
    ----------
    r_red : float array
    r_tmp : float array
    isite : int
    jsite : int
    
    """
    _uppasd.f90wrap_read_exchange_getneighvec(r_red=r_red, r_tmp=r_tmp, isite=isite, \
        jsite=jsite)

def read_anisotropy_alloy():
    """
    read_anisotropy_alloy()
    
    
    Defined at /home/andersb/CrossPlatform/UppASD/source/Input/inputhandler_ext.f90 \
        lines 909-920
    
    
    """
    _uppasd.f90wrap_read_anisotropy_alloy()

def read_anisotropy():
    """
    read_anisotropy()
    
    
    Defined at /home/andersb/CrossPlatform/UppASD/source/Input/inputhandler_ext.f90 \
        lines 926-940
    
    
    """
    _uppasd.f90wrap_read_anisotropy()

def read_dmdata():
    """
    read_dmdata()
    
    
    Defined at /home/andersb/CrossPlatform/UppASD/source/Input/inputhandler_ext.f90 \
        lines 946-1082
    
    
    """
    _uppasd.f90wrap_read_dmdata()

def read_chirdata():
    """
    read_chirdata()
    
    
    Defined at /home/andersb/CrossPlatform/UppASD/source/Input/inputhandler_ext.f90 \
        lines 1094-1226
    
    
    """
    _uppasd.f90wrap_read_chirdata()

def read_fourxdata():
    """
    read_fourxdata()
    
    
    Defined at /home/andersb/CrossPlatform/UppASD/source/Input/inputhandler_ext.f90 \
        lines 1237-1369
    
    
    """
    _uppasd.f90wrap_read_fourxdata()

def read_pddata():
    """
    read_pddata()
    
    
    Defined at /home/andersb/CrossPlatform/UppASD/source/Input/inputhandler_ext.f90 \
        lines 1375-1497
    
    
    """
    _uppasd.f90wrap_read_pddata()

def read_biqdmdata():
    """
    read_biqdmdata()
    
    
    Defined at /home/andersb/CrossPlatform/UppASD/source/Input/inputhandler_ext.f90 \
        lines 1503-1616
    
    
    """
    _uppasd.f90wrap_read_biqdmdata()

def read_bqdata():
    """
    read_bqdata()
    
    
    Defined at /home/andersb/CrossPlatform/UppASD/source/Input/inputhandler_ext.f90 \
        lines 1622-1734
    
    
    """
    _uppasd.f90wrap_read_bqdata()

def read_ringdata():
    """
    read_ringdata()
    
    
    Defined at /home/andersb/CrossPlatform/UppASD/source/Input/inputhandler_ext.f90 \
        lines 1740-1847
    
    
    """
    _uppasd.f90wrap_read_ringdata()

def read_ip_damping():
    """
    read_ip_damping()
    
    
    Defined at /home/andersb/CrossPlatform/UppASD/source/Input/inputhandler_ext.f90 \
        lines 1890-1926
    
    
    """
    _uppasd.f90wrap_read_ip_damping()

def read_ip_damping_alloy():
    """
    read_ip_damping_alloy()
    
    
    Defined at /home/andersb/CrossPlatform/UppASD/source/Input/inputhandler_ext.f90 \
        lines 1932-1968
    
    
    """
    _uppasd.f90wrap_read_ip_damping_alloy()

def read_damping():
    """
    read_damping()
    
    
    Defined at /home/andersb/CrossPlatform/UppASD/source/Input/inputhandler_ext.f90 \
        lines 1974-2008
    
    
    """
    _uppasd.f90wrap_read_damping()

def read_damping_alloy():
    """
    read_damping_alloy()
    
    
    Defined at /home/andersb/CrossPlatform/UppASD/source/Input/inputhandler_ext.f90 \
        lines 2014-2046
    
    
    """
    _uppasd.f90wrap_read_damping_alloy()

def read_barriers():
    """
    read_barriers()
    
    
    Defined at /home/andersb/CrossPlatform/UppASD/source/Input/inputhandler_ext.f90 \
        lines 2057-2119
    
    
    """
    _uppasd.f90wrap_read_barriers()


_array_initialisers = []
_dt_array_initialisers = []

try:
    for func in _array_initialisers:
        func()
except ValueError:
    logging.debug('unallocated array(s) detected on import of module \
        "inputhandler_ext".')

for func in _dt_array_initialisers:
    func()
