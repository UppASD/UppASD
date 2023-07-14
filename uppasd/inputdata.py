"""
Module inputdata


Defined at /home/andersb/CrossPlatform/UppASD/source/Input/inputdata.f90 lines \
    7-792

"""
from __future__ import print_function, absolute_import, division
from uppasd import _uppasd
import f90wrap.runtime
import logging

_arrays = {}
_objs = {}

def set_input_defaults():
    """
    set_input_defaults()
    
    
    Defined at /home/andersb/CrossPlatform/UppASD/source/Input/inputdata.f90 lines \
        280-487
    
    
    """
    _uppasd.f90wrap_set_input_defaults()

def allocate_initmag(flag, na=None, nchmax=None, conf_num=None):
    """
    allocate_initmag(flag[, na, nchmax, conf_num])
    
    
    Defined at /home/andersb/CrossPlatform/UppASD/source/Input/inputdata.f90 lines \
        493-526
    
    Parameters
    ----------
    flag : int
    na : int
    nchmax : int
    conf_num : int
    
    """
    _uppasd.f90wrap_allocate_initmag(flag=flag, na=na, nchmax=nchmax, \
        conf_num=conf_num)

def allocate_chemicalinput(flag, na=None, nchmax=None):
    """
    allocate_chemicalinput(flag[, na, nchmax])
    
    
    Defined at /home/andersb/CrossPlatform/UppASD/source/Input/inputdata.f90 lines \
        532-557
    
    Parameters
    ----------
    flag : int
    na : int
    nchmax : int
    
    """
    _uppasd.f90wrap_allocate_chemicalinput(flag=flag, na=na, nchmax=nchmax)

def reshape_hamiltonianinput():
    """
    reshape_hamiltonianinput()
    
    
    Defined at /home/andersb/CrossPlatform/UppASD/source/Input/inputdata.f90 lines \
        564-606
    
    
    """
    _uppasd.f90wrap_reshape_hamiltonianinput()

def allocate_hamiltonianinput(self, flag, no_shells=None):
    """
    allocate_hamiltonianinput(self, flag[, no_shells])
    
    
    Defined at /home/andersb/CrossPlatform/UppASD/source/Input/inputdata.f90 lines \
        612-785
    
    Parameters
    ----------
    ham_inp : Ham_Inp_T
    flag : int
    no_shells : int
    
    """
    _uppasd.f90wrap_allocate_hamiltonianinput(ham_inp=self._handle, flag=flag, \
        no_shells=no_shells)

def allocate_nn(i):
    """
    allocate_nn(i)
    
    
    Defined at /home/andersb/CrossPlatform/UppASD/source/Input/inputdata.f90 lines \
        787-792
    
    Parameters
    ----------
    i : int
    
    """
    _uppasd.f90wrap_allocate_nn(i=i)

def get_n1():
    """
    Element n1 ftype=integer  pytype=int
    
    
    Defined at /home/andersb/CrossPlatform/UppASD/source/Input/inputdata.f90 line 17
    
    """
    return _uppasd.f90wrap_inputdata__get__n1()

def set_n1(n1):
    _uppasd.f90wrap_inputdata__set__n1(n1)

def get_n2():
    """
    Element n2 ftype=integer  pytype=int
    
    
    Defined at /home/andersb/CrossPlatform/UppASD/source/Input/inputdata.f90 line 18
    
    """
    return _uppasd.f90wrap_inputdata__get__n2()

def set_n2(n2):
    _uppasd.f90wrap_inputdata__set__n2(n2)

def get_n3():
    """
    Element n3 ftype=integer  pytype=int
    
    
    Defined at /home/andersb/CrossPlatform/UppASD/source/Input/inputdata.f90 line 19
    
    """
    return _uppasd.f90wrap_inputdata__get__n3()

def set_n3(n3):
    _uppasd.f90wrap_inputdata__set__n3(n3)

def get_na():
    """
    Element na ftype=integer  pytype=int
    
    
    Defined at /home/andersb/CrossPlatform/UppASD/source/Input/inputdata.f90 line 20
    
    """
    return _uppasd.f90wrap_inputdata__get__na()

def set_na(na):
    _uppasd.f90wrap_inputdata__set__na(na)

def get_nt():
    """
    Element nt ftype=integer  pytype=int
    
    
    Defined at /home/andersb/CrossPlatform/UppASD/source/Input/inputdata.f90 line 21
    
    """
    return _uppasd.f90wrap_inputdata__get__nt()

def set_nt(nt):
    _uppasd.f90wrap_inputdata__set__nt(nt)

def get_sym():
    """
    Element sym ftype=integer  pytype=int
    
    
    Defined at /home/andersb/CrossPlatform/UppASD/source/Input/inputdata.f90 line 22
    
    """
    return _uppasd.f90wrap_inputdata__get__sym()

def set_sym(sym):
    _uppasd.f90wrap_inputdata__set__sym(sym)

def get_nham():
    """
    Element nham ftype=integer  pytype=int
    
    
    Defined at /home/andersb/CrossPlatform/UppASD/source/Input/inputdata.f90 line 23
    
    """
    return _uppasd.f90wrap_inputdata__get__nham()

def set_nham(nham):
    _uppasd.f90wrap_inputdata__set__nham(nham)

def get_natom():
    """
    Element natom ftype=integer  pytype=int
    
    
    Defined at /home/andersb/CrossPlatform/UppASD/source/Input/inputdata.f90 line 24
    
    """
    return _uppasd.f90wrap_inputdata__get__natom()

def set_natom(natom):
    _uppasd.f90wrap_inputdata__set__natom(natom)

def get_natom_full():
    """
    Element natom_full ftype=integer  pytype=int
    
    
    Defined at /home/andersb/CrossPlatform/UppASD/source/Input/inputdata.f90 line 25
    
    """
    return _uppasd.f90wrap_inputdata__get__natom_full()

def set_natom_full(natom_full):
    _uppasd.f90wrap_inputdata__set__natom_full(natom_full)

def get_set_landeg():
    """
    Element set_landeg ftype=integer  pytype=int
    
    
    Defined at /home/andersb/CrossPlatform/UppASD/source/Input/inputdata.f90 line 26
    
    """
    return _uppasd.f90wrap_inputdata__get__set_landeg()

def set_set_landeg(set_landeg):
    _uppasd.f90wrap_inputdata__set__set_landeg(set_landeg)

def get_block_size():
    """
    Element block_size ftype=integer  pytype=int
    
    
    Defined at /home/andersb/CrossPlatform/UppASD/source/Input/inputdata.f90 line 27
    
    """
    return _uppasd.f90wrap_inputdata__get__block_size()

def set_block_size(block_size):
    _uppasd.f90wrap_inputdata__set__block_size(block_size)

def get_metatype():
    """
    Element metatype ftype=integer  pytype=int
    
    
    Defined at /home/andersb/CrossPlatform/UppASD/source/Input/inputdata.f90 line 28
    
    """
    return _uppasd.f90wrap_inputdata__get__metatype()

def set_metatype(metatype):
    _uppasd.f90wrap_inputdata__set__metatype(metatype)

def get_metanumb():
    """
    Element metanumb ftype=integer  pytype=int
    
    
    Defined at /home/andersb/CrossPlatform/UppASD/source/Input/inputdata.f90 line 29
    
    """
    return _uppasd.f90wrap_inputdata__get__metanumb()

def set_metanumb(metanumb):
    _uppasd.f90wrap_inputdata__set__metanumb(metanumb)

def get_array_acomp():
    """
    Element acomp ftype=integer pytype=int
    
    
    Defined at /home/andersb/CrossPlatform/UppASD/source/Input/inputdata.f90 line 30
    
    """
    global acomp
    array_ndim, array_type, array_shape, array_handle = \
        _uppasd.f90wrap_inputdata__array__acomp(f90wrap.runtime.empty_handle)
    if array_handle in _arrays:
        acomp = _arrays[array_handle]
    else:
        acomp = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                f90wrap.runtime.empty_handle,
                                _uppasd.f90wrap_inputdata__array__acomp)
        _arrays[array_handle] = acomp
    return acomp

def set_array_acomp(acomp):
    acomp[...] = acomp

def get_array_asite():
    """
    Element asite ftype=integer pytype=int
    
    
    Defined at /home/andersb/CrossPlatform/UppASD/source/Input/inputdata.f90 line 31
    
    """
    global asite
    array_ndim, array_type, array_shape, array_handle = \
        _uppasd.f90wrap_inputdata__array__asite(f90wrap.runtime.empty_handle)
    if array_handle in _arrays:
        asite = _arrays[array_handle]
    else:
        asite = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                f90wrap.runtime.empty_handle,
                                _uppasd.f90wrap_inputdata__array__asite)
        _arrays[array_handle] = asite
    return asite

def set_array_asite(asite):
    asite[...] = asite

def get_array_anumb_inp():
    """
    Element anumb_inp ftype=integer pytype=int
    
    
    Defined at /home/andersb/CrossPlatform/UppASD/source/Input/inputdata.f90 line 32
    
    """
    global anumb_inp
    array_ndim, array_type, array_shape, array_handle = \
        _uppasd.f90wrap_inputdata__array__anumb_inp(f90wrap.runtime.empty_handle)
    if array_handle in _arrays:
        anumb_inp = _arrays[array_handle]
    else:
        anumb_inp = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                f90wrap.runtime.empty_handle,
                                _uppasd.f90wrap_inputdata__array__anumb_inp)
        _arrays[array_handle] = anumb_inp
    return anumb_inp

def set_array_anumb_inp(anumb_inp):
    anumb_inp[...] = anumb_inp

def get_array_atype_inp():
    """
    Element atype_inp ftype=integer pytype=int
    
    
    Defined at /home/andersb/CrossPlatform/UppASD/source/Input/inputdata.f90 line 33
    
    """
    global atype_inp
    array_ndim, array_type, array_shape, array_handle = \
        _uppasd.f90wrap_inputdata__array__atype_inp(f90wrap.runtime.empty_handle)
    if array_handle in _arrays:
        atype_inp = _arrays[array_handle]
    else:
        atype_inp = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                f90wrap.runtime.empty_handle,
                                _uppasd.f90wrap_inputdata__array__atype_inp)
        _arrays[array_handle] = atype_inp
    return atype_inp

def set_array_atype_inp(atype_inp):
    atype_inp[...] = atype_inp

def get_posfiletype():
    """
    Element posfiletype ftype=character  pytype=str
    
    
    Defined at /home/andersb/CrossPlatform/UppASD/source/Input/inputdata.f90 line 34
    
    """
    return _uppasd.f90wrap_inputdata__get__posfiletype()

def set_posfiletype(posfiletype):
    _uppasd.f90wrap_inputdata__set__posfiletype(posfiletype)

def get_do_sortcoup():
    """
    Element do_sortcoup ftype=character  pytype=str
    
    
    Defined at /home/andersb/CrossPlatform/UppASD/source/Input/inputdata.f90 line 35
    
    """
    return _uppasd.f90wrap_inputdata__get__do_sortcoup()

def set_do_sortcoup(do_sortcoup):
    _uppasd.f90wrap_inputdata__set__do_sortcoup(do_sortcoup)

def get_bc1():
    """
    Element bc1 ftype=character(len=1) pytype=str
    
    
    Defined at /home/andersb/CrossPlatform/UppASD/source/Input/inputdata.f90 line 36
    
    """
    return _uppasd.f90wrap_inputdata__get__bc1()

def set_bc1(bc1):
    _uppasd.f90wrap_inputdata__set__bc1(bc1)

def get_bc2():
    """
    Element bc2 ftype=character(len=1) pytype=str
    
    
    Defined at /home/andersb/CrossPlatform/UppASD/source/Input/inputdata.f90 line 37
    
    """
    return _uppasd.f90wrap_inputdata__get__bc2()

def set_bc2(bc2):
    _uppasd.f90wrap_inputdata__set__bc2(bc2)

def get_bc3():
    """
    Element bc3 ftype=character(len=1) pytype=str
    
    
    Defined at /home/andersb/CrossPlatform/UppASD/source/Input/inputdata.f90 line 38
    
    """
    return _uppasd.f90wrap_inputdata__get__bc3()

def set_bc3(bc3):
    _uppasd.f90wrap_inputdata__set__bc3(bc3)

def get_posfile():
    """
    Element posfile ftype=character(len=35) pytype=str
    
    
    Defined at /home/andersb/CrossPlatform/UppASD/source/Input/inputdata.f90 line 39
    
    """
    return _uppasd.f90wrap_inputdata__get__posfile()

def set_posfile(posfile):
    _uppasd.f90wrap_inputdata__set__posfile(posfile)

def get_alat():
    """
    Element alat ftype=real(dblprec) pytype=float
    
    
    Defined at /home/andersb/CrossPlatform/UppASD/source/Input/inputdata.f90 line 40
    
    """
    return _uppasd.f90wrap_inputdata__get__alat()

def set_alat(alat):
    _uppasd.f90wrap_inputdata__set__alat(alat)

def get_scalefac():
    """
    Element scalefac ftype=real(dblprec) pytype=float
    
    
    Defined at /home/andersb/CrossPlatform/UppASD/source/Input/inputdata.f90 line 41
    
    """
    return _uppasd.f90wrap_inputdata__get__scalefac()

def set_scalefac(scalefac):
    _uppasd.f90wrap_inputdata__set__scalefac(scalefac)

def get_landeg_glob():
    """
    Element landeg_glob ftype=real(dblprec) pytype=float
    
    
    Defined at /home/andersb/CrossPlatform/UppASD/source/Input/inputdata.f90 line 42
    
    """
    return _uppasd.f90wrap_inputdata__get__landeg_glob()

def set_landeg_glob(landeg_glob):
    _uppasd.f90wrap_inputdata__set__landeg_glob(landeg_glob)

def get_array_c1():
    """
    Element c1 ftype=real(dblprec) pytype=float
    
    
    Defined at /home/andersb/CrossPlatform/UppASD/source/Input/inputdata.f90 line 43
    
    """
    global c1
    array_ndim, array_type, array_shape, array_handle = \
        _uppasd.f90wrap_inputdata__array__c1(f90wrap.runtime.empty_handle)
    if array_handle in _arrays:
        c1 = _arrays[array_handle]
    else:
        c1 = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                f90wrap.runtime.empty_handle,
                                _uppasd.f90wrap_inputdata__array__c1)
        _arrays[array_handle] = c1
    return c1

def set_array_c1(c1):
    c1[...] = c1

def get_array_c2():
    """
    Element c2 ftype=real(dblprec) pytype=float
    
    
    Defined at /home/andersb/CrossPlatform/UppASD/source/Input/inputdata.f90 line 44
    
    """
    global c2
    array_ndim, array_type, array_shape, array_handle = \
        _uppasd.f90wrap_inputdata__array__c2(f90wrap.runtime.empty_handle)
    if array_handle in _arrays:
        c2 = _arrays[array_handle]
    else:
        c2 = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                f90wrap.runtime.empty_handle,
                                _uppasd.f90wrap_inputdata__array__c2)
        _arrays[array_handle] = c2
    return c2

def set_array_c2(c2):
    c2[...] = c2

def get_array_c3():
    """
    Element c3 ftype=real(dblprec) pytype=float
    
    
    Defined at /home/andersb/CrossPlatform/UppASD/source/Input/inputdata.f90 line 45
    
    """
    global c3
    array_ndim, array_type, array_shape, array_handle = \
        _uppasd.f90wrap_inputdata__array__c3(f90wrap.runtime.empty_handle)
    if array_handle in _arrays:
        c3 = _arrays[array_handle]
    else:
        c3 = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                f90wrap.runtime.empty_handle,
                                _uppasd.f90wrap_inputdata__array__c3)
        _arrays[array_handle] = c3
    return c3

def set_array_c3(c3):
    c3[...] = c3

def get_array_bas():
    """
    Element bas ftype=real(dblprec) pytype=float
    
    
    Defined at /home/andersb/CrossPlatform/UppASD/source/Input/inputdata.f90 line 46
    
    """
    global bas
    array_ndim, array_type, array_shape, array_handle = \
        _uppasd.f90wrap_inputdata__array__bas(f90wrap.runtime.empty_handle)
    if array_handle in _arrays:
        bas = _arrays[array_handle]
    else:
        bas = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                f90wrap.runtime.empty_handle,
                                _uppasd.f90wrap_inputdata__array__bas)
        _arrays[array_handle] = bas
    return bas

def set_array_bas(bas):
    bas[...] = bas

def get_array_bas0():
    """
    Element bas0 ftype=real(dblprec) pytype=float
    
    
    Defined at /home/andersb/CrossPlatform/UppASD/source/Input/inputdata.f90 line 47
    
    """
    global bas0
    array_ndim, array_type, array_shape, array_handle = \
        _uppasd.f90wrap_inputdata__array__bas0(f90wrap.runtime.empty_handle)
    if array_handle in _arrays:
        bas0 = _arrays[array_handle]
    else:
        bas0 = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                f90wrap.runtime.empty_handle,
                                _uppasd.f90wrap_inputdata__array__bas0)
        _arrays[array_handle] = bas0
    return bas0

def set_array_bas0(bas0):
    bas0[...] = bas0

def get_array_landeg_ch():
    """
    Element landeg_ch ftype=real(dblprec) pytype=float
    
    
    Defined at /home/andersb/CrossPlatform/UppASD/source/Input/inputdata.f90 line 48
    
    """
    global landeg_ch
    array_ndim, array_type, array_shape, array_handle = \
        _uppasd.f90wrap_inputdata__array__landeg_ch(f90wrap.runtime.empty_handle)
    if array_handle in _arrays:
        landeg_ch = _arrays[array_handle]
    else:
        landeg_ch = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                f90wrap.runtime.empty_handle,
                                _uppasd.f90wrap_inputdata__array__landeg_ch)
        _arrays[array_handle] = landeg_ch
    return landeg_ch

def set_array_landeg_ch(landeg_ch):
    landeg_ch[...] = landeg_ch

def get_do_mom_legacy():
    """
    Element do_mom_legacy ftype=character(len=1) pytype=str
    
    
    Defined at /home/andersb/CrossPlatform/UppASD/source/Input/inputdata.f90 line 52
    
    """
    return _uppasd.f90wrap_inputdata__get__do_mom_legacy()

def set_do_mom_legacy(do_mom_legacy):
    _uppasd.f90wrap_inputdata__set__do_mom_legacy(do_mom_legacy)

def get_renorm_coll():
    """
    Element renorm_coll ftype=character(len=1) pytype=str
    
    
    Defined at /home/andersb/CrossPlatform/UppASD/source/Input/inputdata.f90 line 53
    
    """
    return _uppasd.f90wrap_inputdata__get__renorm_coll()

def set_renorm_coll(renorm_coll):
    _uppasd.f90wrap_inputdata__set__renorm_coll(renorm_coll)

def get_ind_mom_flag():
    """
    Element ind_mom_flag ftype=character(len=1) pytype=str
    
    
    Defined at /home/andersb/CrossPlatform/UppASD/source/Input/inputdata.f90 line 54
    
    """
    return _uppasd.f90wrap_inputdata__get__ind_mom_flag()

def set_ind_mom_flag(ind_mom_flag):
    _uppasd.f90wrap_inputdata__set__ind_mom_flag(ind_mom_flag)

def get_ind_mom_type():
    """
    Element ind_mom_type ftype=integer  pytype=int
    
    
    Defined at /home/andersb/CrossPlatform/UppASD/source/Input/inputdata.f90 line 55
    
    """
    return _uppasd.f90wrap_inputdata__get__ind_mom_type()

def set_ind_mom_type(ind_mom_type):
    _uppasd.f90wrap_inputdata__set__ind_mom_type(ind_mom_type)

def get_momfile():
    """
    Element momfile ftype=character(len=35) pytype=str
    
    
    Defined at /home/andersb/CrossPlatform/UppASD/source/Input/inputdata.f90 line 56
    
    """
    return _uppasd.f90wrap_inputdata__get__momfile()

def set_momfile(momfile):
    _uppasd.f90wrap_inputdata__set__momfile(momfile)

def get_momfile_i():
    """
    Element momfile_i ftype=character(len=35) pytype=str
    
    
    Defined at /home/andersb/CrossPlatform/UppASD/source/Input/inputdata.f90 line 57
    
    """
    return _uppasd.f90wrap_inputdata__get__momfile_i()

def set_momfile_i(momfile_i):
    _uppasd.f90wrap_inputdata__set__momfile_i(momfile_i)

def get_momfile_f():
    """
    Element momfile_f ftype=character(len=35) pytype=str
    
    
    Defined at /home/andersb/CrossPlatform/UppASD/source/Input/inputdata.f90 line 58
    
    """
    return _uppasd.f90wrap_inputdata__get__momfile_f()

def set_momfile_f(momfile_f):
    _uppasd.f90wrap_inputdata__set__momfile_f(momfile_f)

def get_ind_tol():
    """
    Element ind_tol ftype=real(dblprec) pytype=float
    
    
    Defined at /home/andersb/CrossPlatform/UppASD/source/Input/inputdata.f90 line 59
    
    """
    return _uppasd.f90wrap_inputdata__get__ind_tol()

def set_ind_tol(ind_tol):
    _uppasd.f90wrap_inputdata__set__ind_tol(ind_tol)

def get_amp_rnd():
    """
    Element amp_rnd ftype=real(dblprec) pytype=float
    
    
    Defined at /home/andersb/CrossPlatform/UppASD/source/Input/inputdata.f90 line 60
    
    """
    return _uppasd.f90wrap_inputdata__get__amp_rnd()

def set_amp_rnd(amp_rnd):
    _uppasd.f90wrap_inputdata__set__amp_rnd(amp_rnd)

def get_amp_rnd_path():
    """
    Element amp_rnd_path ftype=real(dblprec) pytype=float
    
    
    Defined at /home/andersb/CrossPlatform/UppASD/source/Input/inputdata.f90 line 61
    
    """
    return _uppasd.f90wrap_inputdata__get__amp_rnd_path()

def set_amp_rnd_path(amp_rnd_path):
    _uppasd.f90wrap_inputdata__set__amp_rnd_path(amp_rnd_path)

def get_array_ind_mom():
    """
    Element ind_mom ftype=integer pytype=int
    
    
    Defined at /home/andersb/CrossPlatform/UppASD/source/Input/inputdata.f90 line 62
    
    """
    global ind_mom
    array_ndim, array_type, array_shape, array_handle = \
        _uppasd.f90wrap_inputdata__array__ind_mom(f90wrap.runtime.empty_handle)
    if array_handle in _arrays:
        ind_mom = _arrays[array_handle]
    else:
        ind_mom = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                f90wrap.runtime.empty_handle,
                                _uppasd.f90wrap_inputdata__array__ind_mom)
        _arrays[array_handle] = ind_mom
    return ind_mom

def set_array_ind_mom(ind_mom):
    ind_mom[...] = ind_mom

def get_array_ammom_inp():
    """
    Element ammom_inp ftype=real(dblprec) pytype=float
    
    
    Defined at /home/andersb/CrossPlatform/UppASD/source/Input/inputdata.f90 line 63
    
    """
    global ammom_inp
    array_ndim, array_type, array_shape, array_handle = \
        _uppasd.f90wrap_inputdata__array__ammom_inp(f90wrap.runtime.empty_handle)
    if array_handle in _arrays:
        ammom_inp = _arrays[array_handle]
    else:
        ammom_inp = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                f90wrap.runtime.empty_handle,
                                _uppasd.f90wrap_inputdata__array__ammom_inp)
        _arrays[array_handle] = ammom_inp
    return ammom_inp

def set_array_ammom_inp(ammom_inp):
    ammom_inp[...] = ammom_inp

def get_array_aemom_inp():
    """
    Element aemom_inp ftype=real(dblprec) pytype=float
    
    
    Defined at /home/andersb/CrossPlatform/UppASD/source/Input/inputdata.f90 line 64
    
    """
    global aemom_inp
    array_ndim, array_type, array_shape, array_handle = \
        _uppasd.f90wrap_inputdata__array__aemom_inp(f90wrap.runtime.empty_handle)
    if array_handle in _arrays:
        aemom_inp = _arrays[array_handle]
    else:
        aemom_inp = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                f90wrap.runtime.empty_handle,
                                _uppasd.f90wrap_inputdata__array__aemom_inp)
        _arrays[array_handle] = aemom_inp
    return aemom_inp

def set_array_aemom_inp(aemom_inp):
    aemom_inp[...] = aemom_inp

def get_maptype():
    """
    Element maptype ftype=integer  pytype=int
    
    
    Defined at /home/andersb/CrossPlatform/UppASD/source/Input/inputdata.f90 line 65
    
    """
    return _uppasd.f90wrap_inputdata__get__maptype()

def set_maptype(maptype):
    _uppasd.f90wrap_inputdata__set__maptype(maptype)

def get_pre_jfile():
    """
    Element pre_jfile ftype=character(len=30) pytype=str
    
    
    Defined at /home/andersb/CrossPlatform/UppASD/source/Input/inputdata.f90 line 67
    
    """
    return _uppasd.f90wrap_inputdata__get__pre_jfile()

def set_pre_jfile(pre_jfile):
    _uppasd.f90wrap_inputdata__set__pre_jfile(pre_jfile)

def get_array_jfile():
    """
    Element jfile ftype=character(len=30) pytype=str
    
    
    Defined at /home/andersb/CrossPlatform/UppASD/source/Input/inputdata.f90 line 68
    
    """
    global jfile
    array_ndim, array_type, array_shape, array_handle = \
        _uppasd.f90wrap_inputdata__array__jfile(f90wrap.runtime.empty_handle)
    if array_handle in _arrays:
        jfile = _arrays[array_handle]
    else:
        jfile = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                f90wrap.runtime.empty_handle,
                                _uppasd.f90wrap_inputdata__array__jfile)
        _arrays[array_handle] = jfile
    return jfile

def set_array_jfile(jfile):
    jfile[...] = jfile

def get_array_jfiled():
    """
    Element jfiled ftype=character(len=30) pytype=str
    
    
    Defined at /home/andersb/CrossPlatform/UppASD/source/Input/inputdata.f90 line 69
    
    """
    global jfiled
    array_ndim, array_type, array_shape, array_handle = \
        _uppasd.f90wrap_inputdata__array__jfiled(f90wrap.runtime.empty_handle)
    if array_handle in _arrays:
        jfiled = _arrays[array_handle]
    else:
        jfiled = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                f90wrap.runtime.empty_handle,
                                _uppasd.f90wrap_inputdata__array__jfiled)
        _arrays[array_handle] = jfiled
    return jfiled

def set_array_jfiled(jfiled):
    jfiled[...] = jfiled

def get_minalgo():
    """
    Element minalgo ftype=integer  pytype=int
    
    
    Defined at /home/andersb/CrossPlatform/UppASD/source/Input/inputdata.f90 line 73
    
    """
    return _uppasd.f90wrap_inputdata__get__minalgo()

def set_minalgo(minalgo):
    _uppasd.f90wrap_inputdata__set__minalgo(minalgo)

def get_minitrmax():
    """
    Element minitrmax ftype=integer  pytype=int
    
    
    Defined at /home/andersb/CrossPlatform/UppASD/source/Input/inputdata.f90 line 74
    
    """
    return _uppasd.f90wrap_inputdata__get__minitrmax()

def set_minitrmax(minitrmax):
    _uppasd.f90wrap_inputdata__set__minitrmax(minitrmax)

def get_mintraj_step():
    """
    Element mintraj_step ftype=integer  pytype=int
    
    
    Defined at /home/andersb/CrossPlatform/UppASD/source/Input/inputdata.f90 line 75
    
    """
    return _uppasd.f90wrap_inputdata__get__mintraj_step()

def set_mintraj_step(mintraj_step):
    _uppasd.f90wrap_inputdata__set__mintraj_step(mintraj_step)

def get_vpodt():
    """
    Element vpodt ftype=real(dblprec) pytype=float
    
    
    Defined at /home/andersb/CrossPlatform/UppASD/source/Input/inputdata.f90 line 76
    
    """
    return _uppasd.f90wrap_inputdata__get__vpodt()

def set_vpodt(vpodt):
    _uppasd.f90wrap_inputdata__set__vpodt(vpodt)

def get_minftol():
    """
    Element minftol ftype=real(dblprec) pytype=float
    
    
    Defined at /home/andersb/CrossPlatform/UppASD/source/Input/inputdata.f90 line 77
    
    """
    return _uppasd.f90wrap_inputdata__get__minftol()

def set_minftol(minftol):
    _uppasd.f90wrap_inputdata__set__minftol(minftol)

def get_vpomass():
    """
    Element vpomass ftype=real(dblprec) pytype=float
    
    
    Defined at /home/andersb/CrossPlatform/UppASD/source/Input/inputdata.f90 line 78
    
    """
    return _uppasd.f90wrap_inputdata__get__vpomass()

def set_vpomass(vpomass):
    _uppasd.f90wrap_inputdata__set__vpomass(vpomass)

def get_initpath():
    """
    Element initpath ftype=integer  pytype=int
    
    
    Defined at /home/andersb/CrossPlatform/UppASD/source/Input/inputdata.f90 line 82
    
    """
    return _uppasd.f90wrap_inputdata__get__initpath()

def set_initpath(initpath):
    _uppasd.f90wrap_inputdata__set__initpath(initpath)

def get_mepitrmax():
    """
    Element mepitrmax ftype=integer  pytype=int
    
    
    Defined at /home/andersb/CrossPlatform/UppASD/source/Input/inputdata.f90 line 83
    
    """
    return _uppasd.f90wrap_inputdata__get__mepitrmax()

def set_mepitrmax(mepitrmax):
    _uppasd.f90wrap_inputdata__set__mepitrmax(mepitrmax)

def get_meptraj_step():
    """
    Element meptraj_step ftype=integer  pytype=int
    
    
    Defined at /home/andersb/CrossPlatform/UppASD/source/Input/inputdata.f90 line 84
    
    """
    return _uppasd.f90wrap_inputdata__get__meptraj_step()

def set_meptraj_step(meptraj_step):
    _uppasd.f90wrap_inputdata__set__meptraj_step(meptraj_step)

def get_spring():
    """
    Element spring ftype=real(dblprec) pytype=float
    
    
    Defined at /home/andersb/CrossPlatform/UppASD/source/Input/inputdata.f90 line 85
    
    """
    return _uppasd.f90wrap_inputdata__get__spring()

def set_spring(spring):
    _uppasd.f90wrap_inputdata__set__spring(spring)

def get_mepftol():
    """
    Element mepftol ftype=real(dblprec) pytype=float
    
    
    Defined at /home/andersb/CrossPlatform/UppASD/source/Input/inputdata.f90 line 86
    
    """
    return _uppasd.f90wrap_inputdata__get__mepftol()

def set_mepftol(mepftol):
    _uppasd.f90wrap_inputdata__set__mepftol(mepftol)

def get_mepftol_ci():
    """
    Element mepftol_ci ftype=real(dblprec) pytype=float
    
    
    Defined at /home/andersb/CrossPlatform/UppASD/source/Input/inputdata.f90 line 87
    
    """
    return _uppasd.f90wrap_inputdata__get__mepftol_ci()

def set_mepftol_ci(mepftol_ci):
    _uppasd.f90wrap_inputdata__set__mepftol_ci(mepftol_ci)

def get_do_gneb():
    """
    Element do_gneb ftype=character(len=1) pytype=str
    
    
    Defined at /home/andersb/CrossPlatform/UppASD/source/Input/inputdata.f90 line 88
    
    """
    return _uppasd.f90wrap_inputdata__get__do_gneb()

def set_do_gneb(do_gneb):
    _uppasd.f90wrap_inputdata__set__do_gneb(do_gneb)

def get_do_gneb_ci():
    """
    Element do_gneb_ci ftype=character(len=1) pytype=str
    
    
    Defined at /home/andersb/CrossPlatform/UppASD/source/Input/inputdata.f90 line 89
    
    """
    return _uppasd.f90wrap_inputdata__get__do_gneb_ci()

def set_do_gneb_ci(do_gneb_ci):
    _uppasd.f90wrap_inputdata__set__do_gneb_ci(do_gneb_ci)

def get_do_norm_rx():
    """
    Element do_norm_rx ftype=character(len=1) pytype=str
    
    
    Defined at /home/andersb/CrossPlatform/UppASD/source/Input/inputdata.f90 line 90
    
    """
    return _uppasd.f90wrap_inputdata__get__do_norm_rx()

def set_do_norm_rx(do_norm_rx):
    _uppasd.f90wrap_inputdata__set__do_norm_rx(do_norm_rx)

def get_en_zero():
    """
    Element en_zero ftype=character(len=1) pytype=str
    
    
    Defined at /home/andersb/CrossPlatform/UppASD/source/Input/inputdata.f90 line 91
    
    """
    return _uppasd.f90wrap_inputdata__get__en_zero()

def set_en_zero(en_zero):
    _uppasd.f90wrap_inputdata__set__en_zero(en_zero)

def get_relaxed_if():
    """
    Element relaxed_if ftype=character(len=1) pytype=str
    
    
    Defined at /home/andersb/CrossPlatform/UppASD/source/Input/inputdata.f90 line 92
    
    """
    return _uppasd.f90wrap_inputdata__get__relaxed_if()

def set_relaxed_if(relaxed_if):
    _uppasd.f90wrap_inputdata__set__relaxed_if(relaxed_if)

def get_fixed_if():
    """
    Element fixed_if ftype=character(len=1) pytype=str
    
    
    Defined at /home/andersb/CrossPlatform/UppASD/source/Input/inputdata.f90 line 93
    
    """
    return _uppasd.f90wrap_inputdata__get__fixed_if()

def set_fixed_if(fixed_if):
    _uppasd.f90wrap_inputdata__set__fixed_if(fixed_if)

def get_prn_gneb_fields():
    """
    Element prn_gneb_fields ftype=character(len=1) pytype=str
    
    
    Defined at /home/andersb/CrossPlatform/UppASD/source/Input/inputdata.f90 line 94
    
    """
    return _uppasd.f90wrap_inputdata__get__prn_gneb_fields()

def set_prn_gneb_fields(prn_gneb_fields):
    _uppasd.f90wrap_inputdata__set__prn_gneb_fields(prn_gneb_fields)

def get_do_hess_ini():
    """
    Element do_hess_ini ftype=character(len=1) pytype=str
    
    
    Defined at /home/andersb/CrossPlatform/UppASD/source/Input/inputdata.f90 line 98
    
    """
    return _uppasd.f90wrap_inputdata__get__do_hess_ini()

def set_do_hess_ini(do_hess_ini):
    _uppasd.f90wrap_inputdata__set__do_hess_ini(do_hess_ini)

def get_do_hess_fin():
    """
    Element do_hess_fin ftype=character(len=1) pytype=str
    
    
    Defined at /home/andersb/CrossPlatform/UppASD/source/Input/inputdata.f90 line 99
    
    """
    return _uppasd.f90wrap_inputdata__get__do_hess_fin()

def set_do_hess_fin(do_hess_fin):
    _uppasd.f90wrap_inputdata__set__do_hess_fin(do_hess_fin)

def get_do_hess_sp():
    """
    Element do_hess_sp ftype=character(len=1) pytype=str
    
    
    Defined at /home/andersb/CrossPlatform/UppASD/source/Input/inputdata.f90 line \
        100
    
    """
    return _uppasd.f90wrap_inputdata__get__do_hess_sp()

def set_do_hess_sp(do_hess_sp):
    _uppasd.f90wrap_inputdata__set__do_hess_sp(do_hess_sp)

def get_eig_0():
    """
    Element eig_0 ftype=real(dblprec) pytype=float
    
    
    Defined at /home/andersb/CrossPlatform/UppASD/source/Input/inputdata.f90 line \
        101
    
    """
    return _uppasd.f90wrap_inputdata__get__eig_0()

def set_eig_0(eig_0):
    _uppasd.f90wrap_inputdata__set__eig_0(eig_0)

def get_is_afm():
    """
    Element is_afm ftype=character(len=1) pytype=str
    
    
    Defined at /home/andersb/CrossPlatform/UppASD/source/Input/inputdata.f90 line \
        102
    
    """
    return _uppasd.f90wrap_inputdata__get__is_afm()

def set_is_afm(is_afm):
    _uppasd.f90wrap_inputdata__set__is_afm(is_afm)

def get_sample_num():
    """
    Element sample_num ftype=integer  pytype=int
    
    
    Defined at /home/andersb/CrossPlatform/UppASD/source/Input/inputdata.f90 line \
        106
    
    """
    return _uppasd.f90wrap_inputdata__get__sample_num()

def set_sample_num(sample_num):
    _uppasd.f90wrap_inputdata__set__sample_num(sample_num)

def get_simid():
    """
    Element simid ftype=character(len=8) pytype=str
    
    
    Defined at /home/andersb/CrossPlatform/UppASD/source/Input/inputdata.f90 line \
        110
    
    """
    return _uppasd.f90wrap_inputdata__get__simid()

def set_simid(simid):
    _uppasd.f90wrap_inputdata__set__simid(simid)

def get_mensemble():
    """
    Element mensemble ftype=integer  pytype=int
    
    
    Defined at /home/andersb/CrossPlatform/UppASD/source/Input/inputdata.f90 line \
        111
    
    """
    return _uppasd.f90wrap_inputdata__get__mensemble()

def set_mensemble(mensemble):
    _uppasd.f90wrap_inputdata__set__mensemble(mensemble)

def get_tseed():
    """
    Element tseed ftype=integer  pytype=int
    
    
    Defined at /home/andersb/CrossPlatform/UppASD/source/Input/inputdata.f90 line \
        112
    
    """
    return _uppasd.f90wrap_inputdata__get__tseed()

def set_tseed(tseed):
    _uppasd.f90wrap_inputdata__set__tseed(tseed)

def get_llg():
    """
    Element llg ftype=integer  pytype=int
    
    
    Defined at /home/andersb/CrossPlatform/UppASD/source/Input/inputdata.f90 line \
        115
    
    """
    return _uppasd.f90wrap_inputdata__get__llg()

def set_llg(llg):
    _uppasd.f90wrap_inputdata__set__llg(llg)

def get_nstep():
    """
    Element nstep ftype=integer  pytype=int
    
    
    Defined at /home/andersb/CrossPlatform/UppASD/source/Input/inputdata.f90 line \
        116
    
    """
    return _uppasd.f90wrap_inputdata__get__nstep()

def set_nstep(nstep):
    _uppasd.f90wrap_inputdata__set__nstep(nstep)

def get_sdealgh():
    """
    Element sdealgh ftype=integer  pytype=int
    
    
    Defined at /home/andersb/CrossPlatform/UppASD/source/Input/inputdata.f90 line \
        117
    
    """
    return _uppasd.f90wrap_inputdata__get__sdealgh()

def set_sdealgh(sdealgh):
    _uppasd.f90wrap_inputdata__set__sdealgh(sdealgh)

def get_ipsdealgh():
    """
    Element ipsdealgh ftype=integer  pytype=int
    
    
    Defined at /home/andersb/CrossPlatform/UppASD/source/Input/inputdata.f90 line \
        118
    
    """
    return _uppasd.f90wrap_inputdata__get__ipsdealgh()

def set_ipsdealgh(ipsdealgh):
    _uppasd.f90wrap_inputdata__set__ipsdealgh(ipsdealgh)

def get_aunits():
    """
    Element aunits ftype=character(len=1) pytype=str
    
    
    Defined at /home/andersb/CrossPlatform/UppASD/source/Input/inputdata.f90 line \
        119
    
    """
    return _uppasd.f90wrap_inputdata__get__aunits()

def set_aunits(aunits):
    _uppasd.f90wrap_inputdata__set__aunits(aunits)

def get_perp():
    """
    Element perp ftype=character(len=1) pytype=str
    
    
    Defined at /home/andersb/CrossPlatform/UppASD/source/Input/inputdata.f90 line \
        120
    
    """
    return _uppasd.f90wrap_inputdata__get__perp()

def set_perp(perp):
    _uppasd.f90wrap_inputdata__set__perp(perp)

def get_mompar():
    """
    Element mompar ftype=integer  pytype=int
    
    
    Defined at /home/andersb/CrossPlatform/UppASD/source/Input/inputdata.f90 line \
        124
    
    """
    return _uppasd.f90wrap_inputdata__get__mompar()

def set_mompar(mompar):
    _uppasd.f90wrap_inputdata__set__mompar(mompar)

def get_heisout():
    """
    Element heisout ftype=integer  pytype=int
    
    
    Defined at /home/andersb/CrossPlatform/UppASD/source/Input/inputdata.f90 line \
        125
    
    """
    return _uppasd.f90wrap_inputdata__get__heisout()

def set_heisout(heisout):
    _uppasd.f90wrap_inputdata__set__heisout(heisout)

def get_evolveout():
    """
    Element evolveout ftype=integer  pytype=int
    
    
    Defined at /home/andersb/CrossPlatform/UppASD/source/Input/inputdata.f90 line \
        126
    
    """
    return _uppasd.f90wrap_inputdata__get__evolveout()

def set_evolveout(evolveout):
    _uppasd.f90wrap_inputdata__set__evolveout(evolveout)

def get_plotenergy():
    """
    Element plotenergy ftype=integer  pytype=int
    
    
    Defined at /home/andersb/CrossPlatform/UppASD/source/Input/inputdata.f90 line \
        127
    
    """
    return _uppasd.f90wrap_inputdata__get__plotenergy()

def set_plotenergy(plotenergy):
    _uppasd.f90wrap_inputdata__set__plotenergy(plotenergy)

def get_do_hoc_debug():
    """
    Element do_hoc_debug ftype=integer  pytype=int
    
    
    Defined at /home/andersb/CrossPlatform/UppASD/source/Input/inputdata.f90 line \
        128
    
    """
    return _uppasd.f90wrap_inputdata__get__do_hoc_debug()

def set_do_hoc_debug(do_hoc_debug):
    _uppasd.f90wrap_inputdata__set__do_hoc_debug(do_hoc_debug)

def get_do_prnstruct():
    """
    Element do_prnstruct ftype=integer  pytype=int
    
    
    Defined at /home/andersb/CrossPlatform/UppASD/source/Input/inputdata.f90 line \
        129
    
    """
    return _uppasd.f90wrap_inputdata__get__do_prnstruct()

def set_do_prnstruct(do_prnstruct):
    _uppasd.f90wrap_inputdata__set__do_prnstruct(do_prnstruct)

def get_do_storeham():
    """
    Element do_storeham ftype=integer  pytype=int
    
    
    Defined at /home/andersb/CrossPlatform/UppASD/source/Input/inputdata.f90 line \
        130
    
    """
    return _uppasd.f90wrap_inputdata__get__do_storeham()

def set_do_storeham(do_storeham):
    _uppasd.f90wrap_inputdata__set__do_storeham(do_storeham)

def get_do_prn_poscar():
    """
    Element do_prn_poscar ftype=integer  pytype=int
    
    
    Defined at /home/andersb/CrossPlatform/UppASD/source/Input/inputdata.f90 line \
        131
    
    """
    return _uppasd.f90wrap_inputdata__get__do_prn_poscar()

def set_do_prn_poscar(do_prn_poscar):
    _uppasd.f90wrap_inputdata__set__do_prn_poscar(do_prn_poscar)

def get_do_prn_elk():
    """
    Element do_prn_elk ftype=integer  pytype=int
    
    
    Defined at /home/andersb/CrossPlatform/UppASD/source/Input/inputdata.f90 line \
        132
    
    """
    return _uppasd.f90wrap_inputdata__get__do_prn_elk()

def set_do_prn_elk(do_prn_elk):
    _uppasd.f90wrap_inputdata__set__do_prn_elk(do_prn_elk)

def get_do_read_elk():
    """
    Element do_read_elk ftype=integer  pytype=int
    
    
    Defined at /home/andersb/CrossPlatform/UppASD/source/Input/inputdata.f90 line \
        133
    
    """
    return _uppasd.f90wrap_inputdata__get__do_read_elk()

def set_do_read_elk(do_read_elk):
    _uppasd.f90wrap_inputdata__set__do_read_elk(do_read_elk)

def get_compensate_drift():
    """
    Element compensate_drift ftype=integer  pytype=int
    
    
    Defined at /home/andersb/CrossPlatform/UppASD/source/Input/inputdata.f90 line \
        134
    
    """
    return _uppasd.f90wrap_inputdata__get__compensate_drift()

def set_compensate_drift(compensate_drift):
    _uppasd.f90wrap_inputdata__set__compensate_drift(compensate_drift)

def get_do_sparse():
    """
    Element do_sparse ftype=character(len=1) pytype=str
    
    
    Defined at /home/andersb/CrossPlatform/UppASD/source/Input/inputdata.f90 line \
        135
    
    """
    return _uppasd.f90wrap_inputdata__get__do_sparse()

def set_do_sparse(do_sparse):
    _uppasd.f90wrap_inputdata__set__do_sparse(do_sparse)

def get_do_reduced():
    """
    Element do_reduced ftype=character(len=1) pytype=str
    
    
    Defined at /home/andersb/CrossPlatform/UppASD/source/Input/inputdata.f90 line \
        136
    
    """
    return _uppasd.f90wrap_inputdata__get__do_reduced()

def set_do_reduced(do_reduced):
    _uppasd.f90wrap_inputdata__set__do_reduced(do_reduced)

def get_temp():
    """
    Element temp ftype=real(dblprec) pytype=float
    
    
    Defined at /home/andersb/CrossPlatform/UppASD/source/Input/inputdata.f90 line \
        140
    
    """
    return _uppasd.f90wrap_inputdata__get__temp()

def set_temp(temp):
    _uppasd.f90wrap_inputdata__set__temp(temp)

def get_delta_t():
    """
    Element delta_t ftype=real(dblprec) pytype=float
    
    
    Defined at /home/andersb/CrossPlatform/UppASD/source/Input/inputdata.f90 line \
        141
    
    """
    return _uppasd.f90wrap_inputdata__get__delta_t()

def set_delta_t(delta_t):
    _uppasd.f90wrap_inputdata__set__delta_t(delta_t)

def get_relaxtime():
    """
    Element relaxtime ftype=real(dblprec) pytype=float
    
    
    Defined at /home/andersb/CrossPlatform/UppASD/source/Input/inputdata.f90 line \
        142
    
    """
    return _uppasd.f90wrap_inputdata__get__relaxtime()

def set_relaxtime(relaxtime):
    _uppasd.f90wrap_inputdata__set__relaxtime(relaxtime)

def get_mplambda1():
    """
    Element mplambda1 ftype=real(dblprec) pytype=float
    
    
    Defined at /home/andersb/CrossPlatform/UppASD/source/Input/inputdata.f90 line \
        143
    
    """
    return _uppasd.f90wrap_inputdata__get__mplambda1()

def set_mplambda1(mplambda1):
    _uppasd.f90wrap_inputdata__set__mplambda1(mplambda1)

def get_mplambda2():
    """
    Element mplambda2 ftype=real(dblprec) pytype=float
    
    
    Defined at /home/andersb/CrossPlatform/UppASD/source/Input/inputdata.f90 line \
        144
    
    """
    return _uppasd.f90wrap_inputdata__get__mplambda2()

def set_mplambda2(mplambda2):
    _uppasd.f90wrap_inputdata__set__mplambda2(mplambda2)

def get_mode():
    """
    Element mode ftype=character(len=2) pytype=str
    
    
    Defined at /home/andersb/CrossPlatform/UppASD/source/Input/inputdata.f90 line \
        145
    
    """
    return _uppasd.f90wrap_inputdata__get__mode()

def set_mode(mode):
    _uppasd.f90wrap_inputdata__set__mode(mode)

def get_array_hfield():
    """
    Element hfield ftype=real(dblprec) pytype=float
    
    
    Defined at /home/andersb/CrossPlatform/UppASD/source/Input/inputdata.f90 line \
        146
    
    """
    global hfield
    array_ndim, array_type, array_shape, array_handle = \
        _uppasd.f90wrap_inputdata__array__hfield(f90wrap.runtime.empty_handle)
    if array_handle in _arrays:
        hfield = _arrays[array_handle]
    else:
        hfield = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                f90wrap.runtime.empty_handle,
                                _uppasd.f90wrap_inputdata__array__hfield)
        _arrays[array_handle] = hfield
    return hfield

def set_array_hfield(hfield):
    hfield[...] = hfield

def get_mpnlines():
    """
    Element mpnlines ftype=integer  pytype=int
    
    
    Defined at /home/andersb/CrossPlatform/UppASD/source/Input/inputdata.f90 line \
        150
    
    """
    return _uppasd.f90wrap_inputdata__get__mpnlines()

def set_mpnlines(mpnlines):
    _uppasd.f90wrap_inputdata__set__mpnlines(mpnlines)

def get_array_mpdamping1():
    """
    Element mpdamping1 ftype=real(dblprec) pytype=float
    
    
    Defined at /home/andersb/CrossPlatform/UppASD/source/Input/inputdata.f90 line \
        151
    
    """
    global mpdamping1
    array_ndim, array_type, array_shape, array_handle = \
        _uppasd.f90wrap_inputdata__array__mpdamping1(f90wrap.runtime.empty_handle)
    if array_handle in _arrays:
        mpdamping1 = _arrays[array_handle]
    else:
        mpdamping1 = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                f90wrap.runtime.empty_handle,
                                _uppasd.f90wrap_inputdata__array__mpdamping1)
        _arrays[array_handle] = mpdamping1
    return mpdamping1

def set_array_mpdamping1(mpdamping1):
    mpdamping1[...] = mpdamping1

def get_array_mpdamping2():
    """
    Element mpdamping2 ftype=real(dblprec) pytype=float
    
    
    Defined at /home/andersb/CrossPlatform/UppASD/source/Input/inputdata.f90 line \
        152
    
    """
    global mpdamping2
    array_ndim, array_type, array_shape, array_handle = \
        _uppasd.f90wrap_inputdata__array__mpdamping2(f90wrap.runtime.empty_handle)
    if array_handle in _arrays:
        mpdamping2 = _arrays[array_handle]
    else:
        mpdamping2 = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                f90wrap.runtime.empty_handle,
                                _uppasd.f90wrap_inputdata__array__mpdamping2)
        _arrays[array_handle] = mpdamping2
    return mpdamping2

def set_array_mpdamping2(mpdamping2):
    mpdamping2[...] = mpdamping2

def get_array_mpdampingalloy1():
    """
    Element mpdampingalloy1 ftype=real(dblprec) pytype=float
    
    
    Defined at /home/andersb/CrossPlatform/UppASD/source/Input/inputdata.f90 line \
        153
    
    """
    global mpdampingalloy1
    array_ndim, array_type, array_shape, array_handle = \
        _uppasd.f90wrap_inputdata__array__mpdampingalloy1(f90wrap.runtime.empty_handle)
    if array_handle in _arrays:
        mpdampingalloy1 = _arrays[array_handle]
    else:
        mpdampingalloy1 = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                f90wrap.runtime.empty_handle,
                                _uppasd.f90wrap_inputdata__array__mpdampingalloy1)
        _arrays[array_handle] = mpdampingalloy1
    return mpdampingalloy1

def set_array_mpdampingalloy1(mpdampingalloy1):
    mpdampingalloy1[...] = mpdampingalloy1

def get_array_mpdampingalloy2():
    """
    Element mpdampingalloy2 ftype=real(dblprec) pytype=float
    
    
    Defined at /home/andersb/CrossPlatform/UppASD/source/Input/inputdata.f90 line \
        154
    
    """
    global mpdampingalloy2
    array_ndim, array_type, array_shape, array_handle = \
        _uppasd.f90wrap_inputdata__array__mpdampingalloy2(f90wrap.runtime.empty_handle)
    if array_handle in _arrays:
        mpdampingalloy2 = _arrays[array_handle]
    else:
        mpdampingalloy2 = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                f90wrap.runtime.empty_handle,
                                _uppasd.f90wrap_inputdata__array__mpdampingalloy2)
        _arrays[array_handle] = mpdampingalloy2
    return mpdampingalloy2

def set_array_mpdampingalloy2(mpdampingalloy2):
    mpdampingalloy2[...] = mpdampingalloy2

def get_do_site_damping():
    """
    Element do_site_damping ftype=character(len=1) pytype=str
    
    
    Defined at /home/andersb/CrossPlatform/UppASD/source/Input/inputdata.f90 line \
        155
    
    """
    return _uppasd.f90wrap_inputdata__get__do_site_damping()

def set_do_site_damping(do_site_damping):
    _uppasd.f90wrap_inputdata__set__do_site_damping(do_site_damping)

def get_mp_dampfile():
    """
    Element mp_dampfile ftype=character(len=35) pytype=str
    
    
    Defined at /home/andersb/CrossPlatform/UppASD/source/Input/inputdata.f90 line \
        156
    
    """
    return _uppasd.f90wrap_inputdata__get__mp_dampfile()

def set_mp_dampfile(mp_dampfile):
    _uppasd.f90wrap_inputdata__set__mp_dampfile(mp_dampfile)

def get_do_mc():
    """
    Element do_mc ftype=integer  pytype=int
    
    
    Defined at /home/andersb/CrossPlatform/UppASD/source/Input/inputdata.f90 line \
        160
    
    """
    return _uppasd.f90wrap_inputdata__get__do_mc()

def set_do_mc(do_mc):
    _uppasd.f90wrap_inputdata__set__do_mc(do_mc)

def get_mcnstep():
    """
    Element mcnstep ftype=integer  pytype=int
    
    
    Defined at /home/andersb/CrossPlatform/UppASD/source/Input/inputdata.f90 line \
        161
    
    """
    return _uppasd.f90wrap_inputdata__get__mcnstep()

def set_mcnstep(mcnstep):
    _uppasd.f90wrap_inputdata__set__mcnstep(mcnstep)

def get_mcavrg_step():
    """
    Element mcavrg_step ftype=integer  pytype=int
    
    
    Defined at /home/andersb/CrossPlatform/UppASD/source/Input/inputdata.f90 line \
        162
    
    """
    return _uppasd.f90wrap_inputdata__get__mcavrg_step()

def set_mcavrg_step(mcavrg_step):
    _uppasd.f90wrap_inputdata__set__mcavrg_step(mcavrg_step)

def get_mcavrg_buff():
    """
    Element mcavrg_buff ftype=integer  pytype=int
    
    
    Defined at /home/andersb/CrossPlatform/UppASD/source/Input/inputdata.f90 line \
        163
    
    """
    return _uppasd.f90wrap_inputdata__get__mcavrg_buff()

def set_mcavrg_buff(mcavrg_buff):
    _uppasd.f90wrap_inputdata__set__mcavrg_buff(mcavrg_buff)

def get_ipnphase():
    """
    Element ipnphase ftype=integer  pytype=int
    
    
    Defined at /home/andersb/CrossPlatform/UppASD/source/Input/inputdata.f90 line \
        167
    
    """
    return _uppasd.f90wrap_inputdata__get__ipnphase()

def set_ipnphase(ipnphase):
    _uppasd.f90wrap_inputdata__set__ipnphase(ipnphase)

def get_ipmcnphase():
    """
    Element ipmcnphase ftype=integer  pytype=int
    
    
    Defined at /home/andersb/CrossPlatform/UppASD/source/Input/inputdata.f90 line \
        168
    
    """
    return _uppasd.f90wrap_inputdata__get__ipmcnphase()

def set_ipmcnphase(ipmcnphase):
    _uppasd.f90wrap_inputdata__set__ipmcnphase(ipmcnphase)

def get_ipmode():
    """
    Element ipmode ftype=character(len=2) pytype=str
    
    
    Defined at /home/andersb/CrossPlatform/UppASD/source/Input/inputdata.f90 line \
        169
    
    """
    return _uppasd.f90wrap_inputdata__get__ipmode()

def set_ipmode(ipmode):
    _uppasd.f90wrap_inputdata__set__ipmode(ipmode)

def get_array_ipnstep():
    """
    Element ipnstep ftype=integer pytype=int
    
    
    Defined at /home/andersb/CrossPlatform/UppASD/source/Input/inputdata.f90 line \
        170
    
    """
    global ipnstep
    array_ndim, array_type, array_shape, array_handle = \
        _uppasd.f90wrap_inputdata__array__ipnstep(f90wrap.runtime.empty_handle)
    if array_handle in _arrays:
        ipnstep = _arrays[array_handle]
    else:
        ipnstep = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                f90wrap.runtime.empty_handle,
                                _uppasd.f90wrap_inputdata__array__ipnstep)
        _arrays[array_handle] = ipnstep
    return ipnstep

def set_array_ipnstep(ipnstep):
    ipnstep[...] = ipnstep

def get_array_ipmcnstep():
    """
    Element ipmcnstep ftype=integer pytype=int
    
    
    Defined at /home/andersb/CrossPlatform/UppASD/source/Input/inputdata.f90 line \
        171
    
    """
    global ipmcnstep
    array_ndim, array_type, array_shape, array_handle = \
        _uppasd.f90wrap_inputdata__array__ipmcnstep(f90wrap.runtime.empty_handle)
    if array_handle in _arrays:
        ipmcnstep = _arrays[array_handle]
    else:
        ipmcnstep = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                f90wrap.runtime.empty_handle,
                                _uppasd.f90wrap_inputdata__array__ipmcnstep)
        _arrays[array_handle] = ipmcnstep
    return ipmcnstep

def set_array_ipmcnstep(ipmcnstep):
    ipmcnstep[...] = ipmcnstep

def get_array_iphfield():
    """
    Element iphfield ftype=real(dblprec) pytype=float
    
    
    Defined at /home/andersb/CrossPlatform/UppASD/source/Input/inputdata.f90 line \
        172
    
    """
    global iphfield
    array_ndim, array_type, array_shape, array_handle = \
        _uppasd.f90wrap_inputdata__array__iphfield(f90wrap.runtime.empty_handle)
    if array_handle in _arrays:
        iphfield = _arrays[array_handle]
    else:
        iphfield = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                f90wrap.runtime.empty_handle,
                                _uppasd.f90wrap_inputdata__array__iphfield)
        _arrays[array_handle] = iphfield
    return iphfield

def set_array_iphfield(iphfield):
    iphfield[...] = iphfield

def get_array_ipdelta_t():
    """
    Element ipdelta_t ftype=real(dblprec) pytype=float
    
    
    Defined at /home/andersb/CrossPlatform/UppASD/source/Input/inputdata.f90 line \
        173
    
    """
    global ipdelta_t
    array_ndim, array_type, array_shape, array_handle = \
        _uppasd.f90wrap_inputdata__array__ipdelta_t(f90wrap.runtime.empty_handle)
    if array_handle in _arrays:
        ipdelta_t = _arrays[array_handle]
    else:
        ipdelta_t = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                f90wrap.runtime.empty_handle,
                                _uppasd.f90wrap_inputdata__array__ipdelta_t)
        _arrays[array_handle] = ipdelta_t
    return ipdelta_t

def set_array_ipdelta_t(ipdelta_t):
    ipdelta_t[...] = ipdelta_t

def get_array_iplambda1():
    """
    Element iplambda1 ftype=real(dblprec) pytype=float
    
    
    Defined at /home/andersb/CrossPlatform/UppASD/source/Input/inputdata.f90 line \
        174
    
    """
    global iplambda1
    array_ndim, array_type, array_shape, array_handle = \
        _uppasd.f90wrap_inputdata__array__iplambda1(f90wrap.runtime.empty_handle)
    if array_handle in _arrays:
        iplambda1 = _arrays[array_handle]
    else:
        iplambda1 = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                f90wrap.runtime.empty_handle,
                                _uppasd.f90wrap_inputdata__array__iplambda1)
        _arrays[array_handle] = iplambda1
    return iplambda1

def set_array_iplambda1(iplambda1):
    iplambda1[...] = iplambda1

def get_array_iplambda2():
    """
    Element iplambda2 ftype=real(dblprec) pytype=float
    
    
    Defined at /home/andersb/CrossPlatform/UppASD/source/Input/inputdata.f90 line \
        175
    
    """
    global iplambda2
    array_ndim, array_type, array_shape, array_handle = \
        _uppasd.f90wrap_inputdata__array__iplambda2(f90wrap.runtime.empty_handle)
    if array_handle in _arrays:
        iplambda2 = _arrays[array_handle]
    else:
        iplambda2 = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                f90wrap.runtime.empty_handle,
                                _uppasd.f90wrap_inputdata__array__iplambda2)
        _arrays[array_handle] = iplambda2
    return iplambda2

def set_array_iplambda2(iplambda2):
    iplambda2[...] = iplambda2

def get_array_iptemp():
    """
    Element iptemp ftype=real(dblprec) pytype=float
    
    
    Defined at /home/andersb/CrossPlatform/UppASD/source/Input/inputdata.f90 line \
        176
    
    """
    global iptemp
    array_ndim, array_type, array_shape, array_handle = \
        _uppasd.f90wrap_inputdata__array__iptemp(f90wrap.runtime.empty_handle)
    if array_handle in _arrays:
        iptemp = _arrays[array_handle]
    else:
        iptemp = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                f90wrap.runtime.empty_handle,
                                _uppasd.f90wrap_inputdata__array__iptemp)
        _arrays[array_handle] = iptemp
    return iptemp

def set_array_iptemp(iptemp):
    iptemp[...] = iptemp

def get_array_ipdamping1():
    """
    Element ipdamping1 ftype=real(dblprec) pytype=float
    
    
    Defined at /home/andersb/CrossPlatform/UppASD/source/Input/inputdata.f90 line \
        180
    
    """
    global ipdamping1
    array_ndim, array_type, array_shape, array_handle = \
        _uppasd.f90wrap_inputdata__array__ipdamping1(f90wrap.runtime.empty_handle)
    if array_handle in _arrays:
        ipdamping1 = _arrays[array_handle]
    else:
        ipdamping1 = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                f90wrap.runtime.empty_handle,
                                _uppasd.f90wrap_inputdata__array__ipdamping1)
        _arrays[array_handle] = ipdamping1
    return ipdamping1

def set_array_ipdamping1(ipdamping1):
    ipdamping1[...] = ipdamping1

def get_array_ipdamping2():
    """
    Element ipdamping2 ftype=real(dblprec) pytype=float
    
    
    Defined at /home/andersb/CrossPlatform/UppASD/source/Input/inputdata.f90 line \
        181
    
    """
    global ipdamping2
    array_ndim, array_type, array_shape, array_handle = \
        _uppasd.f90wrap_inputdata__array__ipdamping2(f90wrap.runtime.empty_handle)
    if array_handle in _arrays:
        ipdamping2 = _arrays[array_handle]
    else:
        ipdamping2 = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                f90wrap.runtime.empty_handle,
                                _uppasd.f90wrap_inputdata__array__ipdamping2)
        _arrays[array_handle] = ipdamping2
    return ipdamping2

def set_array_ipdamping2(ipdamping2):
    ipdamping2[...] = ipdamping2

def get_array_ipdampingalloy1():
    """
    Element ipdampingalloy1 ftype=real(dblprec) pytype=float
    
    
    Defined at /home/andersb/CrossPlatform/UppASD/source/Input/inputdata.f90 line \
        182
    
    """
    global ipdampingalloy1
    array_ndim, array_type, array_shape, array_handle = \
        _uppasd.f90wrap_inputdata__array__ipdampingalloy1(f90wrap.runtime.empty_handle)
    if array_handle in _arrays:
        ipdampingalloy1 = _arrays[array_handle]
    else:
        ipdampingalloy1 = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                f90wrap.runtime.empty_handle,
                                _uppasd.f90wrap_inputdata__array__ipdampingalloy1)
        _arrays[array_handle] = ipdampingalloy1
    return ipdampingalloy1

def set_array_ipdampingalloy1(ipdampingalloy1):
    ipdampingalloy1[...] = ipdampingalloy1

def get_array_ipdampingalloy2():
    """
    Element ipdampingalloy2 ftype=real(dblprec) pytype=float
    
    
    Defined at /home/andersb/CrossPlatform/UppASD/source/Input/inputdata.f90 line \
        183
    
    """
    global ipdampingalloy2
    array_ndim, array_type, array_shape, array_handle = \
        _uppasd.f90wrap_inputdata__array__ipdampingalloy2(f90wrap.runtime.empty_handle)
    if array_handle in _arrays:
        ipdampingalloy2 = _arrays[array_handle]
    else:
        ipdampingalloy2 = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                f90wrap.runtime.empty_handle,
                                _uppasd.f90wrap_inputdata__array__ipdampingalloy2)
        _arrays[array_handle] = ipdampingalloy2
    return ipdampingalloy2

def set_array_ipdampingalloy2(ipdampingalloy2):
    ipdampingalloy2[...] = ipdampingalloy2

def get_do_site_ip_damping():
    """
    Element do_site_ip_damping ftype=character(len=1) pytype=str
    
    
    Defined at /home/andersb/CrossPlatform/UppASD/source/Input/inputdata.f90 line \
        184
    
    """
    return _uppasd.f90wrap_inputdata__get__do_site_ip_damping()

def set_do_site_ip_damping(do_site_ip_damping):
    _uppasd.f90wrap_inputdata__set__do_site_ip_damping(do_site_ip_damping)

def get_ip_dampfile():
    """
    Element ip_dampfile ftype=character(len=35) pytype=str
    
    
    Defined at /home/andersb/CrossPlatform/UppASD/source/Input/inputdata.f90 line \
        185
    
    """
    return _uppasd.f90wrap_inputdata__get__ip_dampfile()

def set_ip_dampfile(ip_dampfile):
    _uppasd.f90wrap_inputdata__set__ip_dampfile(ip_dampfile)

def get_mseed():
    """
    Element mseed ftype=integer  pytype=int
    
    
    Defined at /home/andersb/CrossPlatform/UppASD/source/Input/inputdata.f90 line \
        189
    
    """
    return _uppasd.f90wrap_inputdata__get__mseed()

def set_mseed(mseed):
    _uppasd.f90wrap_inputdata__set__mseed(mseed)

def get_roteul():
    """
    Element roteul ftype=integer  pytype=int
    
    
    Defined at /home/andersb/CrossPlatform/UppASD/source/Input/inputdata.f90 line \
        190
    
    """
    return _uppasd.f90wrap_inputdata__get__roteul()

def set_roteul(roteul):
    _uppasd.f90wrap_inputdata__set__roteul(roteul)

def get_initmag():
    """
    Element initmag ftype=integer  pytype=int
    
    
    Defined at /home/andersb/CrossPlatform/UppASD/source/Input/inputdata.f90 line \
        191
    
    """
    return _uppasd.f90wrap_inputdata__get__initmag()

def set_initmag(initmag):
    _uppasd.f90wrap_inputdata__set__initmag(initmag)

def get_initneigh():
    """
    Element initneigh ftype=integer  pytype=int
    
    
    Defined at /home/andersb/CrossPlatform/UppASD/source/Input/inputdata.f90 line \
        192
    
    """
    return _uppasd.f90wrap_inputdata__get__initneigh()

def set_initneigh(initneigh):
    _uppasd.f90wrap_inputdata__set__initneigh(initneigh)

def get_phi0():
    """
    Element phi0 ftype=real(dblprec) pytype=float
    
    
    Defined at /home/andersb/CrossPlatform/UppASD/source/Input/inputdata.f90 line \
        193
    
    """
    return _uppasd.f90wrap_inputdata__get__phi0()

def set_phi0(phi0):
    _uppasd.f90wrap_inputdata__set__phi0(phi0)

def get_theta0():
    """
    Element theta0 ftype=real(dblprec) pytype=float
    
    
    Defined at /home/andersb/CrossPlatform/UppASD/source/Input/inputdata.f90 line \
        194
    
    """
    return _uppasd.f90wrap_inputdata__get__theta0()

def set_theta0(theta0):
    _uppasd.f90wrap_inputdata__set__theta0(theta0)

def get_mavg0():
    """
    Element mavg0 ftype=real(dblprec) pytype=float
    
    
    Defined at /home/andersb/CrossPlatform/UppASD/source/Input/inputdata.f90 line \
        195
    
    """
    return _uppasd.f90wrap_inputdata__get__mavg0()

def set_mavg0(mavg0):
    _uppasd.f90wrap_inputdata__set__mavg0(mavg0)

def get_initimp():
    """
    Element initimp ftype=real(dblprec) pytype=float
    
    
    Defined at /home/andersb/CrossPlatform/UppASD/source/Input/inputdata.f90 line \
        196
    
    """
    return _uppasd.f90wrap_inputdata__get__initimp()

def set_initimp(initimp):
    _uppasd.f90wrap_inputdata__set__initimp(initimp)

def get_initconc():
    """
    Element initconc ftype=real(dblprec) pytype=float
    
    
    Defined at /home/andersb/CrossPlatform/UppASD/source/Input/inputdata.f90 line \
        197
    
    """
    return _uppasd.f90wrap_inputdata__get__initconc()

def set_initconc(initconc):
    _uppasd.f90wrap_inputdata__set__initconc(initconc)

def get_initrotang():
    """
    Element initrotang ftype=real(dblprec) pytype=float
    
    
    Defined at /home/andersb/CrossPlatform/UppASD/source/Input/inputdata.f90 line \
        198
    
    """
    return _uppasd.f90wrap_inputdata__get__initrotang()

def set_initrotang(initrotang):
    _uppasd.f90wrap_inputdata__set__initrotang(initrotang)

def get_array_rotang():
    """
    Element rotang ftype=real(dblprec) pytype=float
    
    
    Defined at /home/andersb/CrossPlatform/UppASD/source/Input/inputdata.f90 line \
        199
    
    """
    global rotang
    array_ndim, array_type, array_shape, array_handle = \
        _uppasd.f90wrap_inputdata__array__rotang(f90wrap.runtime.empty_handle)
    if array_handle in _arrays:
        rotang = _arrays[array_handle]
    else:
        rotang = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                f90wrap.runtime.empty_handle,
                                _uppasd.f90wrap_inputdata__array__rotang)
        _arrays[array_handle] = rotang
    return rotang

def set_array_rotang(rotang):
    rotang[...] = rotang

def get_array_initrotvec():
    """
    Element initrotvec ftype=real(dblprec) pytype=float
    
    
    Defined at /home/andersb/CrossPlatform/UppASD/source/Input/inputdata.f90 line \
        200
    
    """
    global initrotvec
    array_ndim, array_type, array_shape, array_handle = \
        _uppasd.f90wrap_inputdata__array__initrotvec(f90wrap.runtime.empty_handle)
    if array_handle in _arrays:
        initrotvec = _arrays[array_handle]
    else:
        initrotvec = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                f90wrap.runtime.empty_handle,
                                _uppasd.f90wrap_inputdata__array__initrotvec)
        _arrays[array_handle] = initrotvec
    return initrotvec

def set_array_initrotvec(initrotvec):
    initrotvec[...] = initrotvec

def get_array_initpropvec():
    """
    Element initpropvec ftype=real(dblprec) pytype=float
    
    
    Defined at /home/andersb/CrossPlatform/UppASD/source/Input/inputdata.f90 line \
        201
    
    """
    global initpropvec
    array_ndim, array_type, array_shape, array_handle = \
        _uppasd.f90wrap_inputdata__array__initpropvec(f90wrap.runtime.empty_handle)
    if array_handle in _arrays:
        initpropvec = _arrays[array_handle]
    else:
        initpropvec = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                f90wrap.runtime.empty_handle,
                                _uppasd.f90wrap_inputdata__array__initpropvec)
        _arrays[array_handle] = initpropvec
    return initpropvec

def set_array_initpropvec(initpropvec):
    initpropvec[...] = initpropvec

def get_initexc():
    """
    Element initexc ftype=character(len=1) pytype=str
    
    
    Defined at /home/andersb/CrossPlatform/UppASD/source/Input/inputdata.f90 line \
        202
    
    """
    return _uppasd.f90wrap_inputdata__get__initexc()

def set_initexc(initexc):
    _uppasd.f90wrap_inputdata__set__initexc(initexc)

def get_restartfile():
    """
    Element restartfile ftype=character(len=35) pytype=str
    
    
    Defined at /home/andersb/CrossPlatform/UppASD/source/Input/inputdata.f90 line \
        203
    
    """
    return _uppasd.f90wrap_inputdata__get__restartfile()

def set_restartfile(restartfile):
    _uppasd.f90wrap_inputdata__set__restartfile(restartfile)

def get_demagvol():
    """
    Element demagvol ftype=real(dblprec) pytype=float
    
    
    Defined at /home/andersb/CrossPlatform/UppASD/source/Input/inputdata.f90 line \
        207
    
    """
    return _uppasd.f90wrap_inputdata__get__demagvol()

def set_demagvol(demagvol):
    _uppasd.f90wrap_inputdata__set__demagvol(demagvol)

def get_demag():
    """
    Element demag ftype=character(len=1) pytype=str
    
    
    Defined at /home/andersb/CrossPlatform/UppASD/source/Input/inputdata.f90 line \
        208
    
    """
    return _uppasd.f90wrap_inputdata__get__demag()

def set_demag(demag):
    _uppasd.f90wrap_inputdata__set__demag(demag)

def get_demag1():
    """
    Element demag1 ftype=character(len=1) pytype=str
    
    
    Defined at /home/andersb/CrossPlatform/UppASD/source/Input/inputdata.f90 line \
        209
    
    """
    return _uppasd.f90wrap_inputdata__get__demag1()

def set_demag1(demag1):
    _uppasd.f90wrap_inputdata__set__demag1(demag1)

def get_demag2():
    """
    Element demag2 ftype=character(len=1) pytype=str
    
    
    Defined at /home/andersb/CrossPlatform/UppASD/source/Input/inputdata.f90 line \
        210
    
    """
    return _uppasd.f90wrap_inputdata__get__demag2()

def set_demag2(demag2):
    _uppasd.f90wrap_inputdata__set__demag2(demag2)

def get_demag3():
    """
    Element demag3 ftype=character(len=1) pytype=str
    
    
    Defined at /home/andersb/CrossPlatform/UppASD/source/Input/inputdata.f90 line \
        211
    
    """
    return _uppasd.f90wrap_inputdata__get__demag3()

def set_demag3(demag3):
    _uppasd.f90wrap_inputdata__set__demag3(demag3)

def get_array_sitenatomfld():
    """
    Element sitenatomfld ftype=real(dblprec) pytype=float
    
    
    Defined at /home/andersb/CrossPlatform/UppASD/source/Input/inputdata.f90 line \
        215
    
    """
    global sitenatomfld
    array_ndim, array_type, array_shape, array_handle = \
        _uppasd.f90wrap_inputdata__array__sitenatomfld(f90wrap.runtime.empty_handle)
    if array_handle in _arrays:
        sitenatomfld = _arrays[array_handle]
    else:
        sitenatomfld = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                f90wrap.runtime.empty_handle,
                                _uppasd.f90wrap_inputdata__array__sitenatomfld)
        _arrays[array_handle] = sitenatomfld
    return sitenatomfld

def set_array_sitenatomfld(sitenatomfld):
    sitenatomfld[...] = sitenatomfld

def get_do_bpulse():
    """
    Element do_bpulse ftype=integer  pytype=int
    
    
    Defined at /home/andersb/CrossPlatform/UppASD/source/Input/inputdata.f90 line \
        219
    
    """
    return _uppasd.f90wrap_inputdata__get__do_bpulse()

def set_do_bpulse(do_bpulse):
    _uppasd.f90wrap_inputdata__set__do_bpulse(do_bpulse)

def get_locfield():
    """
    Element locfield ftype=character  pytype=str
    
    
    Defined at /home/andersb/CrossPlatform/UppASD/source/Input/inputdata.f90 line \
        220
    
    """
    return _uppasd.f90wrap_inputdata__get__locfield()

def set_locfield(locfield):
    _uppasd.f90wrap_inputdata__set__locfield(locfield)

def get_bpulsefile():
    """
    Element bpulsefile ftype=character(len=35) pytype=str
    
    
    Defined at /home/andersb/CrossPlatform/UppASD/source/Input/inputdata.f90 line \
        221
    
    """
    return _uppasd.f90wrap_inputdata__get__bpulsefile()

def set_bpulsefile(bpulsefile):
    _uppasd.f90wrap_inputdata__set__bpulsefile(bpulsefile)

def get_siteatomfile():
    """
    Element siteatomfile ftype=character(len=35) pytype=str
    
    
    Defined at /home/andersb/CrossPlatform/UppASD/source/Input/inputdata.f90 line \
        222
    
    """
    return _uppasd.f90wrap_inputdata__get__siteatomfile()

def set_siteatomfile(siteatomfile):
    _uppasd.f90wrap_inputdata__set__siteatomfile(siteatomfile)

def get_locfieldfile():
    """
    Element locfieldfile ftype=character(len=35) pytype=str
    
    
    Defined at /home/andersb/CrossPlatform/UppASD/source/Input/inputdata.f90 line \
        223
    
    """
    return _uppasd.f90wrap_inputdata__get__locfieldfile()

def set_locfieldfile(locfieldfile):
    _uppasd.f90wrap_inputdata__set__locfieldfile(locfieldfile)

def get_conf_num():
    """
    Element conf_num ftype=integer  pytype=int
    
    
    Defined at /home/andersb/CrossPlatform/UppASD/source/Input/inputdata.f90 line \
        227
    
    """
    return _uppasd.f90wrap_inputdata__get__conf_num()

def set_conf_num(conf_num):
    _uppasd.f90wrap_inputdata__set__conf_num(conf_num)

def get_gsconf_num():
    """
    Element gsconf_num ftype=integer  pytype=int
    
    
    Defined at /home/andersb/CrossPlatform/UppASD/source/Input/inputdata.f90 line \
        228
    
    """
    return _uppasd.f90wrap_inputdata__get__gsconf_num()

def set_gsconf_num(gsconf_num):
    _uppasd.f90wrap_inputdata__set__gsconf_num(gsconf_num)

def get_lsf_metric():
    """
    Element lsf_metric ftype=integer  pytype=int
    
    
    Defined at /home/andersb/CrossPlatform/UppASD/source/Input/inputdata.f90 line \
        229
    
    """
    return _uppasd.f90wrap_inputdata__get__lsf_metric()

def set_lsf_metric(lsf_metric):
    _uppasd.f90wrap_inputdata__set__lsf_metric(lsf_metric)

def get_lsf_window():
    """
    Element lsf_window ftype=real(dblprec) pytype=float
    
    
    Defined at /home/andersb/CrossPlatform/UppASD/source/Input/inputdata.f90 line \
        230
    
    """
    return _uppasd.f90wrap_inputdata__get__lsf_window()

def set_lsf_window(lsf_window):
    _uppasd.f90wrap_inputdata__set__lsf_window(lsf_window)

def get_do_lsf():
    """
    Element do_lsf ftype=character(len=1) pytype=str
    
    
    Defined at /home/andersb/CrossPlatform/UppASD/source/Input/inputdata.f90 line \
        231
    
    """
    return _uppasd.f90wrap_inputdata__get__do_lsf()

def set_do_lsf(do_lsf):
    _uppasd.f90wrap_inputdata__set__do_lsf(do_lsf)

def get_lsf_field():
    """
    Element lsf_field ftype=character(len=1) pytype=str
    
    
    Defined at /home/andersb/CrossPlatform/UppASD/source/Input/inputdata.f90 line \
        232
    
    """
    return _uppasd.f90wrap_inputdata__get__lsf_field()

def set_lsf_field(lsf_field):
    _uppasd.f90wrap_inputdata__set__lsf_field(lsf_field)

def get_lsf_interpolate():
    """
    Element lsf_interpolate ftype=character(len=1) pytype=str
    
    
    Defined at /home/andersb/CrossPlatform/UppASD/source/Input/inputdata.f90 line \
        233
    
    """
    return _uppasd.f90wrap_inputdata__get__lsf_interpolate()

def set_lsf_interpolate(lsf_interpolate):
    _uppasd.f90wrap_inputdata__set__lsf_interpolate(lsf_interpolate)

def get_lsffile():
    """
    Element lsffile ftype=character(len=35) pytype=str
    
    
    Defined at /home/andersb/CrossPlatform/UppASD/source/Input/inputdata.f90 line \
        234
    
    """
    return _uppasd.f90wrap_inputdata__get__lsffile()

def set_lsffile(lsffile):
    _uppasd.f90wrap_inputdata__set__lsffile(lsffile)

def get_spintemp_step():
    """
    Element spintemp_step ftype=integer  pytype=int
    
    
    Defined at /home/andersb/CrossPlatform/UppASD/source/Input/inputdata.f90 line \
        238
    
    """
    return _uppasd.f90wrap_inputdata__get__spintemp_step()

def set_spintemp_step(spintemp_step):
    _uppasd.f90wrap_inputdata__set__spintemp_step(spintemp_step)

def get_logsamp():
    """
    Element logsamp ftype=character(len=1) pytype=str
    
    
    Defined at /home/andersb/CrossPlatform/UppASD/source/Input/inputdata.f90 line \
        239
    
    """
    return _uppasd.f90wrap_inputdata__get__logsamp()

def set_logsamp(logsamp):
    _uppasd.f90wrap_inputdata__set__logsamp(logsamp)

def get_do_spintemp():
    """
    Element do_spintemp ftype=character(len=1) pytype=str
    
    
    Defined at /home/andersb/CrossPlatform/UppASD/source/Input/inputdata.f90 line \
        240
    
    """
    return _uppasd.f90wrap_inputdata__get__do_spintemp()

def set_do_spintemp(do_spintemp):
    _uppasd.f90wrap_inputdata__set__do_spintemp(do_spintemp)

def get_real_time_measure():
    """
    Element real_time_measure ftype=character(len=1) pytype=str
    
    
    Defined at /home/andersb/CrossPlatform/UppASD/source/Input/inputdata.f90 line \
        241
    
    """
    return _uppasd.f90wrap_inputdata__get__real_time_measure()

def set_real_time_measure(real_time_measure):
    _uppasd.f90wrap_inputdata__set__real_time_measure(real_time_measure)

def get_nchmax():
    """
    Element nchmax ftype=integer  pytype=int
    
    
    Defined at /home/andersb/CrossPlatform/UppASD/source/Input/inputdata.f90 line \
        245
    
    """
    return _uppasd.f90wrap_inputdata__get__nchmax()

def set_nchmax(nchmax):
    _uppasd.f90wrap_inputdata__set__nchmax(nchmax)

def get_do_ralloy():
    """
    Element do_ralloy ftype=integer  pytype=int
    
    
    Defined at /home/andersb/CrossPlatform/UppASD/source/Input/inputdata.f90 line \
        246
    
    """
    return _uppasd.f90wrap_inputdata__get__do_ralloy()

def set_do_ralloy(do_ralloy):
    _uppasd.f90wrap_inputdata__set__do_ralloy(do_ralloy)

def get_array_nch():
    """
    Element nch ftype=integer pytype=int
    
    
    Defined at /home/andersb/CrossPlatform/UppASD/source/Input/inputdata.f90 line \
        247
    
    """
    global nch
    array_ndim, array_type, array_shape, array_handle = \
        _uppasd.f90wrap_inputdata__array__nch(f90wrap.runtime.empty_handle)
    if array_handle in _arrays:
        nch = _arrays[array_handle]
    else:
        nch = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                f90wrap.runtime.empty_handle,
                                _uppasd.f90wrap_inputdata__array__nch)
        _arrays[array_handle] = nch
    return nch

def set_array_nch(nch):
    nch[...] = nch

def get_array_achtype_ch():
    """
    Element achtype_ch ftype=integer pytype=int
    
    
    Defined at /home/andersb/CrossPlatform/UppASD/source/Input/inputdata.f90 line \
        248
    
    """
    global achtype_ch
    array_ndim, array_type, array_shape, array_handle = \
        _uppasd.f90wrap_inputdata__array__achtype_ch(f90wrap.runtime.empty_handle)
    if array_handle in _arrays:
        achtype_ch = _arrays[array_handle]
    else:
        achtype_ch = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                f90wrap.runtime.empty_handle,
                                _uppasd.f90wrap_inputdata__array__achtype_ch)
        _arrays[array_handle] = achtype_ch
    return achtype_ch

def set_array_achtype_ch(achtype_ch):
    achtype_ch[...] = achtype_ch

def get_array_chconc():
    """
    Element chconc ftype=real(dblprec) pytype=float
    
    
    Defined at /home/andersb/CrossPlatform/UppASD/source/Input/inputdata.f90 line \
        249
    
    """
    global chconc
    array_ndim, array_type, array_shape, array_handle = \
        _uppasd.f90wrap_inputdata__array__chconc(f90wrap.runtime.empty_handle)
    if array_handle in _arrays:
        chconc = _arrays[array_handle]
    else:
        chconc = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                f90wrap.runtime.empty_handle,
                                _uppasd.f90wrap_inputdata__array__chconc)
        _arrays[array_handle] = chconc
    return chconc

def set_array_chconc(chconc):
    chconc[...] = chconc

def get_oldformat():
    """
    Element oldformat ftype=logical pytype=bool
    
    
    Defined at /home/andersb/CrossPlatform/UppASD/source/Input/inputdata.f90 line \
        253
    
    """
    return _uppasd.f90wrap_inputdata__get__oldformat()

def set_oldformat(oldformat):
    _uppasd.f90wrap_inputdata__set__oldformat(oldformat)

def get_rngpol():
    """
    Element rngpol ftype=character(len=1) pytype=str
    
    
    Defined at /home/andersb/CrossPlatform/UppASD/source/Input/inputdata.f90 line \
        257
    
    """
    return _uppasd.f90wrap_inputdata__get__rngpol()

def set_rngpol(rngpol):
    _uppasd.f90wrap_inputdata__set__rngpol(rngpol)

def get_ziggurat():
    """
    Element ziggurat ftype=character(len=1) pytype=str
    
    
    Defined at /home/andersb/CrossPlatform/UppASD/source/Input/inputdata.f90 line \
        258
    
    """
    return _uppasd.f90wrap_inputdata__get__ziggurat()

def set_ziggurat(ziggurat):
    _uppasd.f90wrap_inputdata__set__ziggurat(ziggurat)

def get_gpu_mode():
    """
    Element gpu_mode ftype=integer  pytype=int
    
    
    Defined at /home/andersb/CrossPlatform/UppASD/source/Input/inputdata.f90 line \
        262
    
    """
    return _uppasd.f90wrap_inputdata__get__gpu_mode()

def set_gpu_mode(gpu_mode):
    _uppasd.f90wrap_inputdata__set__gpu_mode(gpu_mode)

def get_gpu_rng():
    """
    Element gpu_rng ftype=integer  pytype=int
    
    
    Defined at /home/andersb/CrossPlatform/UppASD/source/Input/inputdata.f90 line \
        263
    
    """
    return _uppasd.f90wrap_inputdata__get__gpu_rng()

def set_gpu_rng(gpu_rng):
    _uppasd.f90wrap_inputdata__set__gpu_rng(gpu_rng)

def get_gpu_rng_seed():
    """
    Element gpu_rng_seed ftype=integer  pytype=int
    
    
    Defined at /home/andersb/CrossPlatform/UppASD/source/Input/inputdata.f90 line \
        264
    
    """
    return _uppasd.f90wrap_inputdata__get__gpu_rng_seed()

def set_gpu_rng_seed(gpu_rng_seed):
    _uppasd.f90wrap_inputdata__set__gpu_rng_seed(gpu_rng_seed)

def get_prn_ovf():
    """
    Element prn_ovf ftype=character(len=1) pytype=str
    
    
    Defined at /home/andersb/CrossPlatform/UppASD/source/Input/inputdata.f90 line \
        268
    
    """
    return _uppasd.f90wrap_inputdata__get__prn_ovf()

def set_prn_ovf(prn_ovf):
    _uppasd.f90wrap_inputdata__set__prn_ovf(prn_ovf)

def get_read_ovf():
    """
    Element read_ovf ftype=character(len=1) pytype=str
    
    
    Defined at /home/andersb/CrossPlatform/UppASD/source/Input/inputdata.f90 line \
        269
    
    """
    return _uppasd.f90wrap_inputdata__get__read_ovf()

def set_read_ovf(read_ovf):
    _uppasd.f90wrap_inputdata__set__read_ovf(read_ovf)

def get_do_multiscale():
    """
    Element do_multiscale ftype=logical pytype=bool
    
    
    Defined at /home/andersb/CrossPlatform/UppASD/source/Input/inputdata.f90 line \
        271
    
    """
    return _uppasd.f90wrap_inputdata__get__do_multiscale()

def set_do_multiscale(do_multiscale):
    _uppasd.f90wrap_inputdata__set__do_multiscale(do_multiscale)

def get_do_prnmultiscale():
    """
    Element do_prnmultiscale ftype=logical pytype=bool
    
    
    Defined at /home/andersb/CrossPlatform/UppASD/source/Input/inputdata.f90 line \
        272
    
    """
    return _uppasd.f90wrap_inputdata__get__do_prnmultiscale()

def set_do_prnmultiscale(do_prnmultiscale):
    _uppasd.f90wrap_inputdata__set__do_prnmultiscale(do_prnmultiscale)

def get_multiscale_file_name():
    """
    Element multiscale_file_name ftype=character(len=260) pytype=str
    
    
    Defined at /home/andersb/CrossPlatform/UppASD/source/Input/inputdata.f90 line \
        273
    
    """
    return _uppasd.f90wrap_inputdata__get__multiscale_file_name()

def set_multiscale_file_name(multiscale_file_name):
    _uppasd.f90wrap_inputdata__set__multiscale_file_name(multiscale_file_name)

def get_multiscale_old_format():
    """
    Element multiscale_old_format ftype=character(len=1) pytype=str
    
    
    Defined at /home/andersb/CrossPlatform/UppASD/source/Input/inputdata.f90 line \
        274
    
    """
    return _uppasd.f90wrap_inputdata__get__multiscale_old_format()

def set_multiscale_old_format(multiscale_old_format):
    _uppasd.f90wrap_inputdata__set__multiscale_old_format(multiscale_old_format)


_array_initialisers = [get_array_acomp, get_array_asite, get_array_anumb_inp, \
    get_array_atype_inp, get_array_c1, get_array_c2, get_array_c3, \
    get_array_bas, get_array_bas0, get_array_landeg_ch, get_array_ind_mom, \
    get_array_ammom_inp, get_array_aemom_inp, get_array_jfile, get_array_jfiled, \
    get_array_hfield, get_array_mpdamping1, get_array_mpdamping2, \
    get_array_mpdampingalloy1, get_array_mpdampingalloy2, get_array_ipnstep, \
    get_array_ipmcnstep, get_array_iphfield, get_array_ipdelta_t, \
    get_array_iplambda1, get_array_iplambda2, get_array_iptemp, \
    get_array_ipdamping1, get_array_ipdamping2, get_array_ipdampingalloy1, \
    get_array_ipdampingalloy2, get_array_rotang, get_array_initrotvec, \
    get_array_initpropvec, get_array_sitenatomfld, get_array_nch, \
    get_array_achtype_ch, get_array_chconc]
_dt_array_initialisers = []

try:
    for func in _array_initialisers:
        func()
except ValueError:
    logging.debug('unallocated array(s) detected on import of module "inputdata".')

for func in _dt_array_initialisers:
    func()
