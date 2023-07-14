"""
Module momentdata


Defined at /home/andersb/CrossPlatform/UppASD/source/System/momentdata.f90 lines \
    4-118

"""
from __future__ import print_function, absolute_import, division
from uppasd import _uppasd
import f90wrap.runtime
import logging

_arrays = {}
_objs = {}

def initializemomentdata(moments, natom, mensemble):
    """
    initializemomentdata(moments, natom, mensemble)
    
    
    Defined at /home/andersb/CrossPlatform/UppASD/source/System/momentdata.f90 lines \
        22-42
    
    Parameters
    ----------
    moments : float array
    natom : int
    mensemble : int
    
    """
    _uppasd.f90wrap_initializemomentdata(moments=moments, natom=natom, \
        mensemble=mensemble)

def allocate_multiscale(flag, natom=None, mensemble=None):
    """
    allocate_multiscale(flag[, natom, mensemble])
    
    
    Defined at /home/andersb/CrossPlatform/UppASD/source/System/momentdata.f90 lines \
        45-61
    
    Parameters
    ----------
    flag : int
    natom : int
    mensemble : int
    
    """
    _uppasd.f90wrap_allocate_multiscale(flag=flag, natom=natom, mensemble=mensemble)

def allocate_mmoms(flag, natom=None, mensemble=None):
    """
    allocate_mmoms(flag[, natom, mensemble])
    
    
    Defined at /home/andersb/CrossPlatform/UppASD/source/System/momentdata.f90 lines \
        64-92
    
    Parameters
    ----------
    flag : int
    natom : int
    mensemble : int
    
    """
    _uppasd.f90wrap_allocate_mmoms(flag=flag, natom=natom, mensemble=mensemble)

def allocate_emoms(natom, mensemble, flag):
    """
    allocate_emoms(natom, mensemble, flag)
    
    
    Defined at /home/andersb/CrossPlatform/UppASD/source/System/momentdata.f90 lines \
        95-118
    
    Parameters
    ----------
    natom : int
    mensemble : int
    flag : int
    
    """
    _uppasd.f90wrap_allocate_emoms(natom=natom, mensemble=mensemble, flag=flag)

def get_array_emom():
    """
    Element emom ftype=real(dblprec) pytype=float
    
    
    Defined at /home/andersb/CrossPlatform/UppASD/source/System/momentdata.f90 line \
        10
    
    """
    global emom
    array_ndim, array_type, array_shape, array_handle = \
        _uppasd.f90wrap_momentdata__array__emom(f90wrap.runtime.empty_handle)
    if array_handle in _arrays:
        emom = _arrays[array_handle]
    else:
        emom = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                f90wrap.runtime.empty_handle,
                                _uppasd.f90wrap_momentdata__array__emom)
        _arrays[array_handle] = emom
    return emom

def set_array_emom(emom):
    emom[...] = emom

def get_array_mmom():
    """
    Element mmom ftype=real(dblprec) pytype=float
    
    
    Defined at /home/andersb/CrossPlatform/UppASD/source/System/momentdata.f90 line \
        11
    
    """
    global mmom
    array_ndim, array_type, array_shape, array_handle = \
        _uppasd.f90wrap_momentdata__array__mmom(f90wrap.runtime.empty_handle)
    if array_handle in _arrays:
        mmom = _arrays[array_handle]
    else:
        mmom = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                f90wrap.runtime.empty_handle,
                                _uppasd.f90wrap_momentdata__array__mmom)
        _arrays[array_handle] = mmom
    return mmom

def set_array_mmom(mmom):
    mmom[...] = mmom

def get_array_mmomi():
    """
    Element mmomi ftype=real(dblprec) pytype=float
    
    
    Defined at /home/andersb/CrossPlatform/UppASD/source/System/momentdata.f90 line \
        12
    
    """
    global mmomi
    array_ndim, array_type, array_shape, array_handle = \
        _uppasd.f90wrap_momentdata__array__mmomi(f90wrap.runtime.empty_handle)
    if array_handle in _arrays:
        mmomi = _arrays[array_handle]
    else:
        mmomi = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                f90wrap.runtime.empty_handle,
                                _uppasd.f90wrap_momentdata__array__mmomi)
        _arrays[array_handle] = mmomi
    return mmomi

def set_array_mmomi(mmomi):
    mmomi[...] = mmomi

def get_array_emom2():
    """
    Element emom2 ftype=real(dblprec) pytype=float
    
    
    Defined at /home/andersb/CrossPlatform/UppASD/source/System/momentdata.f90 line \
        13
    
    """
    global emom2
    array_ndim, array_type, array_shape, array_handle = \
        _uppasd.f90wrap_momentdata__array__emom2(f90wrap.runtime.empty_handle)
    if array_handle in _arrays:
        emom2 = _arrays[array_handle]
    else:
        emom2 = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                f90wrap.runtime.empty_handle,
                                _uppasd.f90wrap_momentdata__array__emom2)
        _arrays[array_handle] = emom2
    return emom2

def set_array_emom2(emom2):
    emom2[...] = emom2

def get_array_emomm():
    """
    Element emomm ftype=real(dblprec) pytype=float
    
    
    Defined at /home/andersb/CrossPlatform/UppASD/source/System/momentdata.f90 line \
        14
    
    """
    global emomm
    array_ndim, array_type, array_shape, array_handle = \
        _uppasd.f90wrap_momentdata__array__emomm(f90wrap.runtime.empty_handle)
    if array_handle in _arrays:
        emomm = _arrays[array_handle]
    else:
        emomm = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                f90wrap.runtime.empty_handle,
                                _uppasd.f90wrap_momentdata__array__emomm)
        _arrays[array_handle] = emomm
    return emomm

def set_array_emomm(emomm):
    emomm[...] = emomm

def get_array_mmom2():
    """
    Element mmom2 ftype=real(dblprec) pytype=float
    
    
    Defined at /home/andersb/CrossPlatform/UppASD/source/System/momentdata.f90 line \
        15
    
    """
    global mmom2
    array_ndim, array_type, array_shape, array_handle = \
        _uppasd.f90wrap_momentdata__array__mmom2(f90wrap.runtime.empty_handle)
    if array_handle in _arrays:
        mmom2 = _arrays[array_handle]
    else:
        mmom2 = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                f90wrap.runtime.empty_handle,
                                _uppasd.f90wrap_momentdata__array__mmom2)
        _arrays[array_handle] = mmom2
    return mmom2

def set_array_mmom2(mmom2):
    mmom2[...] = mmom2

def get_array_mmom0():
    """
    Element mmom0 ftype=real(dblprec) pytype=float
    
    
    Defined at /home/andersb/CrossPlatform/UppASD/source/System/momentdata.f90 line \
        16
    
    """
    global mmom0
    array_ndim, array_type, array_shape, array_handle = \
        _uppasd.f90wrap_momentdata__array__mmom0(f90wrap.runtime.empty_handle)
    if array_handle in _arrays:
        mmom0 = _arrays[array_handle]
    else:
        mmom0 = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                f90wrap.runtime.empty_handle,
                                _uppasd.f90wrap_momentdata__array__mmom0)
        _arrays[array_handle] = mmom0
    return mmom0

def set_array_mmom0(mmom0):
    mmom0[...] = mmom0

def get_array_multiscalebackbuffer():
    """
    Element multiscalebackbuffer ftype=real(dblprec) pytype=float
    
    
    Defined at /home/andersb/CrossPlatform/UppASD/source/System/momentdata.f90 line \
        17
    
    """
    global multiscalebackbuffer
    array_ndim, array_type, array_shape, array_handle = \
        _uppasd.f90wrap_momentdata__array__multiscalebackbuffer(f90wrap.runtime.empty_handle)
    if array_handle in _arrays:
        multiscalebackbuffer = _arrays[array_handle]
    else:
        multiscalebackbuffer = \
            f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                f90wrap.runtime.empty_handle,
                                _uppasd.f90wrap_momentdata__array__multiscalebackbuffer)
        _arrays[array_handle] = multiscalebackbuffer
    return multiscalebackbuffer

def set_array_multiscalebackbuffer(multiscalebackbuffer):
    multiscalebackbuffer[...] = multiscalebackbuffer

def get_multiscalebackbufferhead():
    """
    Element multiscalebackbufferhead ftype=integer  pytype=int
    
    
    Defined at /home/andersb/CrossPlatform/UppASD/source/System/momentdata.f90 line \
        18
    
    """
    return _uppasd.f90wrap_momentdata__get__multiscalebackbufferhead()

def set_multiscalebackbufferhead(multiscalebackbufferhead):
    _uppasd.f90wrap_momentdata__set__multiscalebackbufferhead(multiscalebackbufferhead)


_array_initialisers = [get_array_emom, get_array_mmom, get_array_mmomi, \
    get_array_emom2, get_array_emomm, get_array_mmom2, get_array_mmom0, \
    get_array_multiscalebackbuffer]
_dt_array_initialisers = []

try:
    for func in _array_initialisers:
        func()
except ValueError:
    logging.debug('unallocated array(s) detected on import of module "momentdata".')

for func in _dt_array_initialisers:
    func()
