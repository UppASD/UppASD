"""
Module uppasd


Defined at /home/andersb/CrossPlatform/UppASD/source/uppasd.f90 lines 41-1494

"""
from __future__ import print_function, absolute_import, division
from uppasd import _uppasd
import f90wrap.runtime
import logging

_arrays = {}
_objs = {}

def main():
    """
    main()
    
    
    Defined at /home/andersb/CrossPlatform/UppASD/source/uppasd.f90 lines 65-110
    
    
    ==============================================================
     Check if inpsd.dat exists and whether it contains anything
    --------------------------------------------------------------
    """
    _uppasd.f90wrap_main()

def number_of_active_processors():
    """
    nprocs = number_of_active_processors()
    
    
    Defined at /home/andersb/CrossPlatform/UppASD/source/uppasd.f90 lines 123-132
    
    
    Returns
    -------
    nprocs : int
    
    """
    nprocs = _uppasd.f90wrap_number_of_active_processors()
    return nprocs

def run_initial_phase():
    """
    run_initial_phase()
    
    
    Defined at /home/andersb/CrossPlatform/UppASD/source/uppasd.f90 lines 149-216
    
    
    """
    _uppasd.f90wrap_run_initial_phase()

def run_measurement_phase():
    """
    run_measurement_phase()
    
    
    Defined at /home/andersb/CrossPlatform/UppASD/source/uppasd.f90 lines 231-362
    
    
    """
    _uppasd.f90wrap_run_measurement_phase()

def cleanup_simulation():
    """
    cleanup_simulation()
    
    
    Defined at /home/andersb/CrossPlatform/UppASD/source/uppasd.f90 lines 374-475
    
    
    """
    _uppasd.f90wrap_cleanup_simulation()

def setup_simulation():
    """
    setup_simulation()
    
    
    Defined at /home/andersb/CrossPlatform/UppASD/source/uppasd.f90 lines 484-1236
    
    
    ----------------------------
     Print warning if warranted
    ----------------------------
    """
    _uppasd.f90wrap_setup_simulation()

def sd_timing(nstep):
    """
    sd_timing(nstep)
    
    
    Defined at /home/andersb/CrossPlatform/UppASD/source/uppasd.f90 lines 1249-1267
    
    Parameters
    ----------
    nstep : int
    
    """
    _uppasd.f90wrap_sd_timing(nstep=nstep)

def allocate_general(flag):
    """
    allocate_general(flag)
    
    
    Defined at /home/andersb/CrossPlatform/UppASD/source/uppasd.f90 lines 1276-1328
    
    Parameters
    ----------
    flag : int
    
    """
    _uppasd.f90wrap_allocate_general(flag=flag)

def print_logo():
    """
    print_logo()
    
    
    Defined at /home/andersb/CrossPlatform/UppASD/source/uppasd.f90 lines 1337-1389
    
    
    """
    _uppasd.f90wrap_print_logo()

def print_siminfo():
    """
    print_siminfo()
    
    
    Defined at /home/andersb/CrossPlatform/UppASD/source/uppasd.f90 lines 1398-1436
    
    
    """
    _uppasd.f90wrap_print_siminfo()

def check_format():
    """
    check_format()
    
    
    Defined at /home/andersb/CrossPlatform/UppASD/source/uppasd.f90 lines 1448-1462
    
    
    """
    _uppasd.f90wrap_check_format()

def calculate_energy(outenergy=None):
    """
    calculate_energy([outenergy])
    
    
    Defined at /home/andersb/CrossPlatform/UppASD/source/uppasd.f90 lines 1464-1494
    
    Parameters
    ----------
    outenergy : float
    
    """
    _uppasd.f90wrap_calculate_energy(outenergy=outenergy)

def get_time_a():
    """
    Element time_a ftype=real(dblprec) pytype=float
    
    
    Defined at /home/andersb/CrossPlatform/UppASD/source/uppasd.f90 line 61
    
    """
    return _uppasd.f90wrap_uppasd__get__time_a()

def set_time_a(time_a):
    _uppasd.f90wrap_uppasd__set__time_a(time_a)

def get_time_b():
    """
    Element time_b ftype=real(dblprec) pytype=float
    
    
    Defined at /home/andersb/CrossPlatform/UppASD/source/uppasd.f90 line 61
    
    """
    return _uppasd.f90wrap_uppasd__get__time_b()

def set_time_b(time_b):
    _uppasd.f90wrap_uppasd__set__time_b(time_b)

def get_time_c():
    """
    Element time_c ftype=real(dblprec) pytype=float
    
    
    Defined at /home/andersb/CrossPlatform/UppASD/source/uppasd.f90 line 61
    
    """
    return _uppasd.f90wrap_uppasd__get__time_c()

def set_time_c(time_c):
    _uppasd.f90wrap_uppasd__set__time_c(time_c)

def get_time_d():
    """
    Element time_d ftype=real(dblprec) pytype=float
    
    
    Defined at /home/andersb/CrossPlatform/UppASD/source/uppasd.f90 line 61
    
    """
    return _uppasd.f90wrap_uppasd__get__time_d()

def set_time_d(time_d):
    _uppasd.f90wrap_uppasd__set__time_d(time_d)

def get_time_e():
    """
    Element time_e ftype=real(dblprec) pytype=float
    
    
    Defined at /home/andersb/CrossPlatform/UppASD/source/uppasd.f90 line 61
    
    """
    return _uppasd.f90wrap_uppasd__get__time_e()

def set_time_e(time_e):
    _uppasd.f90wrap_uppasd__set__time_e(time_e)

def get_time_i():
    """
    Element time_i ftype=real(dblprec) pytype=float
    
    
    Defined at /home/andersb/CrossPlatform/UppASD/source/uppasd.f90 line 62
    
    """
    return _uppasd.f90wrap_uppasd__get__time_i()

def set_time_i(time_i):
    _uppasd.f90wrap_uppasd__set__time_i(time_i)

def get_time_ii():
    """
    Element time_ii ftype=real(dblprec) pytype=float
    
    
    Defined at /home/andersb/CrossPlatform/UppASD/source/uppasd.f90 line 62
    
    """
    return _uppasd.f90wrap_uppasd__get__time_ii()

def set_time_ii(time_ii):
    _uppasd.f90wrap_uppasd__set__time_ii(time_ii)

def get_time_iii():
    """
    Element time_iii ftype=real(dblprec) pytype=float
    
    
    Defined at /home/andersb/CrossPlatform/UppASD/source/uppasd.f90 line 62
    
    """
    return _uppasd.f90wrap_uppasd__get__time_iii()

def set_time_iii(time_iii):
    _uppasd.f90wrap_uppasd__set__time_iii(time_iii)

def get_time_iv():
    """
    Element time_iv ftype=real(dblprec) pytype=float
    
    
    Defined at /home/andersb/CrossPlatform/UppASD/source/uppasd.f90 line 62
    
    """
    return _uppasd.f90wrap_uppasd__get__time_iv()

def set_time_iv(time_iv):
    _uppasd.f90wrap_uppasd__set__time_iv(time_iv)

def get_time_v():
    """
    Element time_v ftype=real(dblprec) pytype=float
    
    
    Defined at /home/andersb/CrossPlatform/UppASD/source/uppasd.f90 line 62
    
    """
    return _uppasd.f90wrap_uppasd__get__time_v()

def set_time_v(time_v):
    _uppasd.f90wrap_uppasd__set__time_v(time_v)

def get_nprocs():
    """
    Element nprocs ftype=integer  pytype=int
    
    
    Defined at /home/andersb/CrossPlatform/UppASD/source/uppasd.f90 line 63
    
    """
    return _uppasd.f90wrap_uppasd__get__nprocs()

def set_nprocs(nprocs):
    _uppasd.f90wrap_uppasd__set__nprocs(nprocs)


_array_initialisers = []
_dt_array_initialisers = []

try:
    for func in _array_initialisers:
        func()
except ValueError:
    logging.debug('unallocated array(s) detected on import of module "uppasd".')

for func in _dt_array_initialisers:
    func()
