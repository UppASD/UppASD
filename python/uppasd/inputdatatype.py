"""
Module inputdatatype


Defined at /home/andersb/CrossPlatform/UppASD/source/Input/inputdatatype.f90 \
    lines 7-125

"""
from __future__ import print_function, absolute_import, division
from uppasd import _uppasd
import f90wrap.runtime
import logging

_arrays = {}
_objs = {}

@f90wrap.runtime.register_class("uppasd.ham_inp_t")
class ham_inp_t(f90wrap.runtime.FortranDerivedType):
    """
    Type(name=ham_inp_t)
    
    
    Defined at /home/andersb/CrossPlatform/UppASD/source/Input/inputdatatype.f90 \
        lines 11-124
    
    """
    def __init__(self, handle=None):
        """
        self = Ham_Inp_T()
        
        
        Defined at /home/andersb/CrossPlatform/UppASD/source/Input/inputdatatype.f90 \
            lines 11-124
        
        
        Returns
        -------
        this : Ham_Inp_T
        	Object to be constructed
        
        
        Automatically generated constructor for ham_inp_t
        """
        f90wrap.runtime.FortranDerivedType.__init__(self)
        result = _uppasd.f90wrap_ham_inp_t_initialise()
        self._handle = result[0] if isinstance(result, tuple) else result
    
    def __del__(self):
        """
        Destructor for class Ham_Inp_T
        
        
        Defined at /home/andersb/CrossPlatform/UppASD/source/Input/inputdatatype.f90 \
            lines 11-124
        
        Parameters
        ----------
        this : Ham_Inp_T
        	Object to be destructed
        
        
        Automatically generated destructor for ham_inp_t
        """
        if self._alloc:
            _uppasd.f90wrap_ham_inp_t_finalise(this=self._handle)
    
    @property
    def do_jtensor(self):
        """
        Element do_jtensor ftype=integer  pytype=int
        
        
        Defined at /home/andersb/CrossPlatform/UppASD/source/Input/inputdatatype.f90 \
            line 16
        
        """
        return _uppasd.f90wrap_ham_inp_t__get__do_jtensor(self._handle)
    
    @do_jtensor.setter
    def do_jtensor(self, do_jtensor):
        _uppasd.f90wrap_ham_inp_t__set__do_jtensor(self._handle, do_jtensor)
    
    @property
    def calc_jtensor(self):
        """
        Element calc_jtensor ftype=logical pytype=bool
        
        
        Defined at /home/andersb/CrossPlatform/UppASD/source/Input/inputdatatype.f90 \
            line 17
        
        """
        return _uppasd.f90wrap_ham_inp_t__get__calc_jtensor(self._handle)
    
    @calc_jtensor.setter
    def calc_jtensor(self, calc_jtensor):
        _uppasd.f90wrap_ham_inp_t__set__calc_jtensor(self._handle, calc_jtensor)
    
    @property
    def map_multiple(self):
        """
        Element map_multiple ftype=logical pytype=bool
        
        
        Defined at /home/andersb/CrossPlatform/UppASD/source/Input/inputdatatype.f90 \
            line 18
        
        """
        return _uppasd.f90wrap_ham_inp_t__get__map_multiple(self._handle)
    
    @map_multiple.setter
    def map_multiple(self, map_multiple):
        _uppasd.f90wrap_ham_inp_t__set__map_multiple(self._handle, map_multiple)
    
    @property
    def max_no_shells(self):
        """
        Element max_no_shells ftype=integer  pytype=int
        
        
        Defined at /home/andersb/CrossPlatform/UppASD/source/Input/inputdatatype.f90 \
            line 19
        
        """
        return _uppasd.f90wrap_ham_inp_t__get__max_no_shells(self._handle)
    
    @max_no_shells.setter
    def max_no_shells(self, max_no_shells):
        _uppasd.f90wrap_ham_inp_t__set__max_no_shells(self._handle, max_no_shells)
    
    @property
    def nn(self):
        """
        Element nn ftype=integer pytype=int
        
        
        Defined at /home/andersb/CrossPlatform/UppASD/source/Input/inputdatatype.f90 \
            line 20
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _uppasd.f90wrap_ham_inp_t__array__nn(self._handle)
        if array_handle in self._arrays:
            nn = self._arrays[array_handle]
        else:
            nn = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _uppasd.f90wrap_ham_inp_t__array__nn)
            self._arrays[array_handle] = nn
        return nn
    
    @nn.setter
    def nn(self, nn):
        self.nn[...] = nn
    
    @property
    def nntype(self):
        """
        Element nntype ftype=integer pytype=int
        
        
        Defined at /home/andersb/CrossPlatform/UppASD/source/Input/inputdatatype.f90 \
            line 21
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _uppasd.f90wrap_ham_inp_t__array__nntype(self._handle)
        if array_handle in self._arrays:
            nntype = self._arrays[array_handle]
        else:
            nntype = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _uppasd.f90wrap_ham_inp_t__array__nntype)
            self._arrays[array_handle] = nntype
        return nntype
    
    @nntype.setter
    def nntype(self, nntype):
        self.nntype[...] = nntype
    
    @property
    def redcoord(self):
        """
        Element redcoord ftype=real(dblprec) pytype=float
        
        
        Defined at /home/andersb/CrossPlatform/UppASD/source/Input/inputdatatype.f90 \
            line 22
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _uppasd.f90wrap_ham_inp_t__array__redcoord(self._handle)
        if array_handle in self._arrays:
            redcoord = self._arrays[array_handle]
        else:
            redcoord = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _uppasd.f90wrap_ham_inp_t__array__redcoord)
            self._arrays[array_handle] = redcoord
        return redcoord
    
    @redcoord.setter
    def redcoord(self, redcoord):
        self.redcoord[...] = redcoord
    
    @property
    def jc(self):
        """
        Element jc ftype=real(dblprec) pytype=float
        
        
        Defined at /home/andersb/CrossPlatform/UppASD/source/Input/inputdatatype.f90 \
            line 23
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _uppasd.f90wrap_ham_inp_t__array__jc(self._handle)
        if array_handle in self._arrays:
            jc = self._arrays[array_handle]
        else:
            jc = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _uppasd.f90wrap_ham_inp_t__array__jc)
            self._arrays[array_handle] = jc
        return jc
    
    @jc.setter
    def jc(self, jc):
        self.jc[...] = jc
    
    @property
    def jcd(self):
        """
        Element jcd ftype=real(dblprec) pytype=float
        
        
        Defined at /home/andersb/CrossPlatform/UppASD/source/Input/inputdatatype.f90 \
            line 24
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _uppasd.f90wrap_ham_inp_t__array__jcd(self._handle)
        if array_handle in self._arrays:
            jcd = self._arrays[array_handle]
        else:
            jcd = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _uppasd.f90wrap_ham_inp_t__array__jcd)
            self._arrays[array_handle] = jcd
        return jcd
    
    @jcd.setter
    def jcd(self, jcd):
        self.jcd[...] = jcd
    
    @property
    def jc_tens(self):
        """
        Element jc_tens ftype=real(dblprec) pytype=float
        
        
        Defined at /home/andersb/CrossPlatform/UppASD/source/Input/inputdatatype.f90 \
            line 25
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _uppasd.f90wrap_ham_inp_t__array__jc_tens(self._handle)
        if array_handle in self._arrays:
            jc_tens = self._arrays[array_handle]
        else:
            jc_tens = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _uppasd.f90wrap_ham_inp_t__array__jc_tens)
            self._arrays[array_handle] = jc_tens
        return jc_tens
    
    @jc_tens.setter
    def jc_tens(self, jc_tens):
        self.jc_tens[...] = jc_tens
    
    @property
    def exc_inter(self):
        """
        Element exc_inter ftype=character(len=1) pytype=str
        
        
        Defined at /home/andersb/CrossPlatform/UppASD/source/Input/inputdatatype.f90 \
            line 26
        
        """
        return _uppasd.f90wrap_ham_inp_t__get__exc_inter(self._handle)
    
    @exc_inter.setter
    def exc_inter(self, exc_inter):
        _uppasd.f90wrap_ham_inp_t__set__exc_inter(self._handle, exc_inter)
    
    @property
    def jij_scale(self):
        """
        Element jij_scale ftype=real(dblprec) pytype=float
        
        
        Defined at /home/andersb/CrossPlatform/UppASD/source/Input/inputdatatype.f90 \
            line 30
        
        """
        return _uppasd.f90wrap_ham_inp_t__get__jij_scale(self._handle)
    
    @jij_scale.setter
    def jij_scale(self, jij_scale):
        _uppasd.f90wrap_ham_inp_t__set__jij_scale(self._handle, jij_scale)
    
    @property
    def ea_model(self):
        """
        Element ea_model ftype=logical pytype=bool
        
        
        Defined at /home/andersb/CrossPlatform/UppASD/source/Input/inputdatatype.f90 \
            line 31
        
        """
        return _uppasd.f90wrap_ham_inp_t__get__ea_model(self._handle)
    
    @ea_model.setter
    def ea_model(self, ea_model):
        _uppasd.f90wrap_ham_inp_t__set__ea_model(self._handle, ea_model)
    
    @property
    def ea_sigma(self):
        """
        Element ea_sigma ftype=real(dblprec) pytype=float
        
        
        Defined at /home/andersb/CrossPlatform/UppASD/source/Input/inputdatatype.f90 \
            line 32
        
        """
        return _uppasd.f90wrap_ham_inp_t__get__ea_sigma(self._handle)
    
    @ea_sigma.setter
    def ea_sigma(self, ea_sigma):
        _uppasd.f90wrap_ham_inp_t__set__ea_sigma(self._handle, ea_sigma)
    
    @property
    def do_anisotropy(self):
        """
        Element do_anisotropy ftype=integer  pytype=int
        
        
        Defined at /home/andersb/CrossPlatform/UppASD/source/Input/inputdatatype.f90 \
            line 36
        
        """
        return _uppasd.f90wrap_ham_inp_t__get__do_anisotropy(self._handle)
    
    @do_anisotropy.setter
    def do_anisotropy(self, do_anisotropy):
        _uppasd.f90wrap_ham_inp_t__set__do_anisotropy(self._handle, do_anisotropy)
    
    @property
    def random_anisotropy_density(self):
        """
        Element random_anisotropy_density ftype=real(dblprec) pytype=float
        
        
        Defined at /home/andersb/CrossPlatform/UppASD/source/Input/inputdatatype.f90 \
            line 37
        
        """
        return _uppasd.f90wrap_ham_inp_t__get__random_anisotropy_density(self._handle)
    
    @random_anisotropy_density.setter
    def random_anisotropy_density(self, random_anisotropy_density):
        _uppasd.f90wrap_ham_inp_t__set__random_anisotropy_density(self._handle, \
            random_anisotropy_density)
    
    @property
    def random_anisotropy(self):
        """
        Element random_anisotropy ftype=logical pytype=bool
        
        
        Defined at /home/andersb/CrossPlatform/UppASD/source/Input/inputdatatype.f90 \
            line 38
        
        """
        return _uppasd.f90wrap_ham_inp_t__get__random_anisotropy(self._handle)
    
    @random_anisotropy.setter
    def random_anisotropy(self, random_anisotropy):
        _uppasd.f90wrap_ham_inp_t__set__random_anisotropy(self._handle, \
            random_anisotropy)
    
    @property
    def mult_axis(self):
        """
        Element mult_axis ftype=character(len=1) pytype=str
        
        
        Defined at /home/andersb/CrossPlatform/UppASD/source/Input/inputdatatype.f90 \
            line 39
        
        """
        return _uppasd.f90wrap_ham_inp_t__get__mult_axis(self._handle)
    
    @mult_axis.setter
    def mult_axis(self, mult_axis):
        _uppasd.f90wrap_ham_inp_t__set__mult_axis(self._handle, mult_axis)
    
    @property
    def kfile(self):
        """
        Element kfile ftype=character(len=35) pytype=str
        
        
        Defined at /home/andersb/CrossPlatform/UppASD/source/Input/inputdatatype.f90 \
            line 40
        
        """
        return _uppasd.f90wrap_ham_inp_t__get__kfile(self._handle)
    
    @kfile.setter
    def kfile(self, kfile):
        _uppasd.f90wrap_ham_inp_t__set__kfile(self._handle, kfile)
    
    @property
    def anisotropytype(self):
        """
        Element anisotropytype ftype=integer pytype=int
        
        
        Defined at /home/andersb/CrossPlatform/UppASD/source/Input/inputdatatype.f90 \
            line 41
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _uppasd.f90wrap_ham_inp_t__array__anisotropytype(self._handle)
        if array_handle in self._arrays:
            anisotropytype = self._arrays[array_handle]
        else:
            anisotropytype = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _uppasd.f90wrap_ham_inp_t__array__anisotropytype)
            self._arrays[array_handle] = anisotropytype
        return anisotropytype
    
    @anisotropytype.setter
    def anisotropytype(self, anisotropytype):
        self.anisotropytype[...] = anisotropytype
    
    @property
    def anisotropytype_diff(self):
        """
        Element anisotropytype_diff ftype=integer pytype=int
        
        
        Defined at /home/andersb/CrossPlatform/UppASD/source/Input/inputdatatype.f90 \
            line 42
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _uppasd.f90wrap_ham_inp_t__array__anisotropytype_diff(self._handle)
        if array_handle in self._arrays:
            anisotropytype_diff = self._arrays[array_handle]
        else:
            anisotropytype_diff = \
                f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _uppasd.f90wrap_ham_inp_t__array__anisotropytype_diff)
            self._arrays[array_handle] = anisotropytype_diff
        return anisotropytype_diff
    
    @anisotropytype_diff.setter
    def anisotropytype_diff(self, anisotropytype_diff):
        self.anisotropytype_diff[...] = anisotropytype_diff
    
    @property
    def anisotropy(self):
        """
        Element anisotropy ftype=real(dblprec) pytype=float
        
        
        Defined at /home/andersb/CrossPlatform/UppASD/source/Input/inputdatatype.f90 \
            line 43
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _uppasd.f90wrap_ham_inp_t__array__anisotropy(self._handle)
        if array_handle in self._arrays:
            anisotropy = self._arrays[array_handle]
        else:
            anisotropy = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _uppasd.f90wrap_ham_inp_t__array__anisotropy)
            self._arrays[array_handle] = anisotropy
        return anisotropy
    
    @anisotropy.setter
    def anisotropy(self, anisotropy):
        self.anisotropy[...] = anisotropy
    
    @property
    def anisotropy_diff(self):
        """
        Element anisotropy_diff ftype=real(dblprec) pytype=float
        
        
        Defined at /home/andersb/CrossPlatform/UppASD/source/Input/inputdatatype.f90 \
            line 44
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _uppasd.f90wrap_ham_inp_t__array__anisotropy_diff(self._handle)
        if array_handle in self._arrays:
            anisotropy_diff = self._arrays[array_handle]
        else:
            anisotropy_diff = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _uppasd.f90wrap_ham_inp_t__array__anisotropy_diff)
            self._arrays[array_handle] = anisotropy_diff
        return anisotropy_diff
    
    @anisotropy_diff.setter
    def anisotropy_diff(self, anisotropy_diff):
        self.anisotropy_diff[...] = anisotropy_diff
    
    @property
    def do_dm(self):
        """
        Element do_dm ftype=integer  pytype=int
        
        
        Defined at /home/andersb/CrossPlatform/UppASD/source/Input/inputdatatype.f90 \
            line 48
        
        """
        return _uppasd.f90wrap_ham_inp_t__get__do_dm(self._handle)
    
    @do_dm.setter
    def do_dm(self, do_dm):
        _uppasd.f90wrap_ham_inp_t__set__do_dm(self._handle, do_dm)
    
    @property
    def max_no_dmshells(self):
        """
        Element max_no_dmshells ftype=integer  pytype=int
        
        
        Defined at /home/andersb/CrossPlatform/UppASD/source/Input/inputdatatype.f90 \
            line 49
        
        """
        return _uppasd.f90wrap_ham_inp_t__get__max_no_dmshells(self._handle)
    
    @max_no_dmshells.setter
    def max_no_dmshells(self, max_no_dmshells):
        _uppasd.f90wrap_ham_inp_t__set__max_no_dmshells(self._handle, max_no_dmshells)
    
    @property
    def dmfile(self):
        """
        Element dmfile ftype=character(len=35) pytype=str
        
        
        Defined at /home/andersb/CrossPlatform/UppASD/source/Input/inputdatatype.f90 \
            line 50
        
        """
        return _uppasd.f90wrap_ham_inp_t__get__dmfile(self._handle)
    
    @dmfile.setter
    def dmfile(self, dmfile):
        _uppasd.f90wrap_ham_inp_t__set__dmfile(self._handle, dmfile)
    
    @property
    def dm_scale(self):
        """
        Element dm_scale ftype=real(dblprec) pytype=float
        
        
        Defined at /home/andersb/CrossPlatform/UppASD/source/Input/inputdatatype.f90 \
            line 51
        
        """
        return _uppasd.f90wrap_ham_inp_t__get__dm_scale(self._handle)
    
    @dm_scale.setter
    def dm_scale(self, dm_scale):
        _uppasd.f90wrap_ham_inp_t__set__dm_scale(self._handle, dm_scale)
    
    @property
    def dm_nn(self):
        """
        Element dm_nn ftype=integer pytype=int
        
        
        Defined at /home/andersb/CrossPlatform/UppASD/source/Input/inputdatatype.f90 \
            line 52
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _uppasd.f90wrap_ham_inp_t__array__dm_nn(self._handle)
        if array_handle in self._arrays:
            dm_nn = self._arrays[array_handle]
        else:
            dm_nn = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _uppasd.f90wrap_ham_inp_t__array__dm_nn)
            self._arrays[array_handle] = dm_nn
        return dm_nn
    
    @dm_nn.setter
    def dm_nn(self, dm_nn):
        self.dm_nn[...] = dm_nn
    
    @property
    def dm_redcoord(self):
        """
        Element dm_redcoord ftype=real(dblprec) pytype=float
        
        
        Defined at /home/andersb/CrossPlatform/UppASD/source/Input/inputdatatype.f90 \
            line 53
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _uppasd.f90wrap_ham_inp_t__array__dm_redcoord(self._handle)
        if array_handle in self._arrays:
            dm_redcoord = self._arrays[array_handle]
        else:
            dm_redcoord = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _uppasd.f90wrap_ham_inp_t__array__dm_redcoord)
            self._arrays[array_handle] = dm_redcoord
        return dm_redcoord
    
    @dm_redcoord.setter
    def dm_redcoord(self, dm_redcoord):
        self.dm_redcoord[...] = dm_redcoord
    
    @property
    def dm_inpvect(self):
        """
        Element dm_inpvect ftype=real(dblprec) pytype=float
        
        
        Defined at /home/andersb/CrossPlatform/UppASD/source/Input/inputdatatype.f90 \
            line 54
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _uppasd.f90wrap_ham_inp_t__array__dm_inpvect(self._handle)
        if array_handle in self._arrays:
            dm_inpvect = self._arrays[array_handle]
        else:
            dm_inpvect = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _uppasd.f90wrap_ham_inp_t__array__dm_inpvect)
            self._arrays[array_handle] = dm_inpvect
        return dm_inpvect
    
    @dm_inpvect.setter
    def dm_inpvect(self, dm_inpvect):
        self.dm_inpvect[...] = dm_inpvect
    
    @property
    def do_chir(self):
        """
        Element do_chir ftype=integer  pytype=int
        
        
        Defined at /home/andersb/CrossPlatform/UppASD/source/Input/inputdatatype.f90 \
            line 58
        
        """
        return _uppasd.f90wrap_ham_inp_t__get__do_chir(self._handle)
    
    @do_chir.setter
    def do_chir(self, do_chir):
        _uppasd.f90wrap_ham_inp_t__set__do_chir(self._handle, do_chir)
    
    @property
    def max_no_chirshells(self):
        """
        Element max_no_chirshells ftype=integer  pytype=int
        
        
        Defined at /home/andersb/CrossPlatform/UppASD/source/Input/inputdatatype.f90 \
            line 59
        
        """
        return _uppasd.f90wrap_ham_inp_t__get__max_no_chirshells(self._handle)
    
    @max_no_chirshells.setter
    def max_no_chirshells(self, max_no_chirshells):
        _uppasd.f90wrap_ham_inp_t__set__max_no_chirshells(self._handle, \
            max_no_chirshells)
    
    @property
    def chirfile(self):
        """
        Element chirfile ftype=character(len=35) pytype=str
        
        
        Defined at /home/andersb/CrossPlatform/UppASD/source/Input/inputdatatype.f90 \
            line 60
        
        """
        return _uppasd.f90wrap_ham_inp_t__get__chirfile(self._handle)
    
    @chirfile.setter
    def chirfile(self, chirfile):
        _uppasd.f90wrap_ham_inp_t__set__chirfile(self._handle, chirfile)
    
    @property
    def chir_nn(self):
        """
        Element chir_nn ftype=integer pytype=int
        
        
        Defined at /home/andersb/CrossPlatform/UppASD/source/Input/inputdatatype.f90 \
            line 61
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _uppasd.f90wrap_ham_inp_t__array__chir_nn(self._handle)
        if array_handle in self._arrays:
            chir_nn = self._arrays[array_handle]
        else:
            chir_nn = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _uppasd.f90wrap_ham_inp_t__array__chir_nn)
            self._arrays[array_handle] = chir_nn
        return chir_nn
    
    @chir_nn.setter
    def chir_nn(self, chir_nn):
        self.chir_nn[...] = chir_nn
    
    @property
    def chir_inpval(self):
        """
        Element chir_inpval ftype=real(dblprec) pytype=float
        
        
        Defined at /home/andersb/CrossPlatform/UppASD/source/Input/inputdatatype.f90 \
            line 62
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _uppasd.f90wrap_ham_inp_t__array__chir_inpval(self._handle)
        if array_handle in self._arrays:
            chir_inpval = self._arrays[array_handle]
        else:
            chir_inpval = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _uppasd.f90wrap_ham_inp_t__array__chir_inpval)
            self._arrays[array_handle] = chir_inpval
        return chir_inpval
    
    @chir_inpval.setter
    def chir_inpval(self, chir_inpval):
        self.chir_inpval[...] = chir_inpval
    
    @property
    def chir_redcoord(self):
        """
        Element chir_redcoord ftype=real(dblprec) pytype=float
        
        
        Defined at /home/andersb/CrossPlatform/UppASD/source/Input/inputdatatype.f90 \
            line 63
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _uppasd.f90wrap_ham_inp_t__array__chir_redcoord(self._handle)
        if array_handle in self._arrays:
            chir_redcoord = self._arrays[array_handle]
        else:
            chir_redcoord = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _uppasd.f90wrap_ham_inp_t__array__chir_redcoord)
            self._arrays[array_handle] = chir_redcoord
        return chir_redcoord
    
    @chir_redcoord.setter
    def chir_redcoord(self, chir_redcoord):
        self.chir_redcoord[...] = chir_redcoord
    
    @property
    def do_fourx(self):
        """
        Element do_fourx ftype=integer  pytype=int
        
        
        Defined at /home/andersb/CrossPlatform/UppASD/source/Input/inputdatatype.f90 \
            line 67
        
        """
        return _uppasd.f90wrap_ham_inp_t__get__do_fourx(self._handle)
    
    @do_fourx.setter
    def do_fourx(self, do_fourx):
        _uppasd.f90wrap_ham_inp_t__set__do_fourx(self._handle, do_fourx)
    
    @property
    def max_no_fourxshells(self):
        """
        Element max_no_fourxshells ftype=integer  pytype=int
        
        
        Defined at /home/andersb/CrossPlatform/UppASD/source/Input/inputdatatype.f90 \
            line 68
        
        """
        return _uppasd.f90wrap_ham_inp_t__get__max_no_fourxshells(self._handle)
    
    @max_no_fourxshells.setter
    def max_no_fourxshells(self, max_no_fourxshells):
        _uppasd.f90wrap_ham_inp_t__set__max_no_fourxshells(self._handle, \
            max_no_fourxshells)
    
    @property
    def fourxfile(self):
        """
        Element fourxfile ftype=character(len=35) pytype=str
        
        
        Defined at /home/andersb/CrossPlatform/UppASD/source/Input/inputdatatype.f90 \
            line 69
        
        """
        return _uppasd.f90wrap_ham_inp_t__get__fourxfile(self._handle)
    
    @fourxfile.setter
    def fourxfile(self, fourxfile):
        _uppasd.f90wrap_ham_inp_t__set__fourxfile(self._handle, fourxfile)
    
    @property
    def fourx_nn(self):
        """
        Element fourx_nn ftype=integer pytype=int
        
        
        Defined at /home/andersb/CrossPlatform/UppASD/source/Input/inputdatatype.f90 \
            line 70
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _uppasd.f90wrap_ham_inp_t__array__fourx_nn(self._handle)
        if array_handle in self._arrays:
            fourx_nn = self._arrays[array_handle]
        else:
            fourx_nn = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _uppasd.f90wrap_ham_inp_t__array__fourx_nn)
            self._arrays[array_handle] = fourx_nn
        return fourx_nn
    
    @fourx_nn.setter
    def fourx_nn(self, fourx_nn):
        self.fourx_nn[...] = fourx_nn
    
    @property
    def fourx_inpval(self):
        """
        Element fourx_inpval ftype=real(dblprec) pytype=float
        
        
        Defined at /home/andersb/CrossPlatform/UppASD/source/Input/inputdatatype.f90 \
            line 71
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _uppasd.f90wrap_ham_inp_t__array__fourx_inpval(self._handle)
        if array_handle in self._arrays:
            fourx_inpval = self._arrays[array_handle]
        else:
            fourx_inpval = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _uppasd.f90wrap_ham_inp_t__array__fourx_inpval)
            self._arrays[array_handle] = fourx_inpval
        return fourx_inpval
    
    @fourx_inpval.setter
    def fourx_inpval(self, fourx_inpval):
        self.fourx_inpval[...] = fourx_inpval
    
    @property
    def fourx_redcoord(self):
        """
        Element fourx_redcoord ftype=real(dblprec) pytype=float
        
        
        Defined at /home/andersb/CrossPlatform/UppASD/source/Input/inputdatatype.f90 \
            line 72
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _uppasd.f90wrap_ham_inp_t__array__fourx_redcoord(self._handle)
        if array_handle in self._arrays:
            fourx_redcoord = self._arrays[array_handle]
        else:
            fourx_redcoord = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _uppasd.f90wrap_ham_inp_t__array__fourx_redcoord)
            self._arrays[array_handle] = fourx_redcoord
        return fourx_redcoord
    
    @fourx_redcoord.setter
    def fourx_redcoord(self, fourx_redcoord):
        self.fourx_redcoord[...] = fourx_redcoord
    
    @property
    def do_pd(self):
        """
        Element do_pd ftype=integer  pytype=int
        
        
        Defined at /home/andersb/CrossPlatform/UppASD/source/Input/inputdatatype.f90 \
            line 76
        
        """
        return _uppasd.f90wrap_ham_inp_t__get__do_pd(self._handle)
    
    @do_pd.setter
    def do_pd(self, do_pd):
        _uppasd.f90wrap_ham_inp_t__set__do_pd(self._handle, do_pd)
    
    @property
    def max_no_pdshells(self):
        """
        Element max_no_pdshells ftype=integer  pytype=int
        
        
        Defined at /home/andersb/CrossPlatform/UppASD/source/Input/inputdatatype.f90 \
            line 77
        
        """
        return _uppasd.f90wrap_ham_inp_t__get__max_no_pdshells(self._handle)
    
    @max_no_pdshells.setter
    def max_no_pdshells(self, max_no_pdshells):
        _uppasd.f90wrap_ham_inp_t__set__max_no_pdshells(self._handle, max_no_pdshells)
    
    @property
    def pdfile(self):
        """
        Element pdfile ftype=character(len=35) pytype=str
        
        
        Defined at /home/andersb/CrossPlatform/UppASD/source/Input/inputdatatype.f90 \
            line 78
        
        """
        return _uppasd.f90wrap_ham_inp_t__get__pdfile(self._handle)
    
    @pdfile.setter
    def pdfile(self, pdfile):
        _uppasd.f90wrap_ham_inp_t__set__pdfile(self._handle, pdfile)
    
    @property
    def pd_nn(self):
        """
        Element pd_nn ftype=integer pytype=int
        
        
        Defined at /home/andersb/CrossPlatform/UppASD/source/Input/inputdatatype.f90 \
            line 79
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _uppasd.f90wrap_ham_inp_t__array__pd_nn(self._handle)
        if array_handle in self._arrays:
            pd_nn = self._arrays[array_handle]
        else:
            pd_nn = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _uppasd.f90wrap_ham_inp_t__array__pd_nn)
            self._arrays[array_handle] = pd_nn
        return pd_nn
    
    @pd_nn.setter
    def pd_nn(self, pd_nn):
        self.pd_nn[...] = pd_nn
    
    @property
    def pd_redcoord(self):
        """
        Element pd_redcoord ftype=real(dblprec) pytype=float
        
        
        Defined at /home/andersb/CrossPlatform/UppASD/source/Input/inputdatatype.f90 \
            line 80
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _uppasd.f90wrap_ham_inp_t__array__pd_redcoord(self._handle)
        if array_handle in self._arrays:
            pd_redcoord = self._arrays[array_handle]
        else:
            pd_redcoord = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _uppasd.f90wrap_ham_inp_t__array__pd_redcoord)
            self._arrays[array_handle] = pd_redcoord
        return pd_redcoord
    
    @pd_redcoord.setter
    def pd_redcoord(self, pd_redcoord):
        self.pd_redcoord[...] = pd_redcoord
    
    @property
    def pd_inpvect(self):
        """
        Element pd_inpvect ftype=real(dblprec) pytype=float
        
        
        Defined at /home/andersb/CrossPlatform/UppASD/source/Input/inputdatatype.f90 \
            line 81
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _uppasd.f90wrap_ham_inp_t__array__pd_inpvect(self._handle)
        if array_handle in self._arrays:
            pd_inpvect = self._arrays[array_handle]
        else:
            pd_inpvect = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _uppasd.f90wrap_ham_inp_t__array__pd_inpvect)
            self._arrays[array_handle] = pd_inpvect
        return pd_inpvect
    
    @pd_inpvect.setter
    def pd_inpvect(self, pd_inpvect):
        self.pd_inpvect[...] = pd_inpvect
    
    @property
    def do_biqdm(self):
        """
        Element do_biqdm ftype=integer  pytype=int
        
        
        Defined at /home/andersb/CrossPlatform/UppASD/source/Input/inputdatatype.f90 \
            line 85
        
        """
        return _uppasd.f90wrap_ham_inp_t__get__do_biqdm(self._handle)
    
    @do_biqdm.setter
    def do_biqdm(self, do_biqdm):
        _uppasd.f90wrap_ham_inp_t__set__do_biqdm(self._handle, do_biqdm)
    
    @property
    def max_no_biqdmshells(self):
        """
        Element max_no_biqdmshells ftype=integer  pytype=int
        
        
        Defined at /home/andersb/CrossPlatform/UppASD/source/Input/inputdatatype.f90 \
            line 86
        
        """
        return _uppasd.f90wrap_ham_inp_t__get__max_no_biqdmshells(self._handle)
    
    @max_no_biqdmshells.setter
    def max_no_biqdmshells(self, max_no_biqdmshells):
        _uppasd.f90wrap_ham_inp_t__set__max_no_biqdmshells(self._handle, \
            max_no_biqdmshells)
    
    @property
    def biqdmfile(self):
        """
        Element biqdmfile ftype=character(len=35) pytype=str
        
        
        Defined at /home/andersb/CrossPlatform/UppASD/source/Input/inputdatatype.f90 \
            line 87
        
        """
        return _uppasd.f90wrap_ham_inp_t__get__biqdmfile(self._handle)
    
    @biqdmfile.setter
    def biqdmfile(self, biqdmfile):
        _uppasd.f90wrap_ham_inp_t__set__biqdmfile(self._handle, biqdmfile)
    
    @property
    def biqdm_nn(self):
        """
        Element biqdm_nn ftype=integer pytype=int
        
        
        Defined at /home/andersb/CrossPlatform/UppASD/source/Input/inputdatatype.f90 \
            line 88
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _uppasd.f90wrap_ham_inp_t__array__biqdm_nn(self._handle)
        if array_handle in self._arrays:
            biqdm_nn = self._arrays[array_handle]
        else:
            biqdm_nn = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _uppasd.f90wrap_ham_inp_t__array__biqdm_nn)
            self._arrays[array_handle] = biqdm_nn
        return biqdm_nn
    
    @biqdm_nn.setter
    def biqdm_nn(self, biqdm_nn):
        self.biqdm_nn[...] = biqdm_nn
    
    @property
    def biqdm_redcoord(self):
        """
        Element biqdm_redcoord ftype=real(dblprec) pytype=float
        
        
        Defined at /home/andersb/CrossPlatform/UppASD/source/Input/inputdatatype.f90 \
            line 89
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _uppasd.f90wrap_ham_inp_t__array__biqdm_redcoord(self._handle)
        if array_handle in self._arrays:
            biqdm_redcoord = self._arrays[array_handle]
        else:
            biqdm_redcoord = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _uppasd.f90wrap_ham_inp_t__array__biqdm_redcoord)
            self._arrays[array_handle] = biqdm_redcoord
        return biqdm_redcoord
    
    @biqdm_redcoord.setter
    def biqdm_redcoord(self, biqdm_redcoord):
        self.biqdm_redcoord[...] = biqdm_redcoord
    
    @property
    def biqdm_inpvect(self):
        """
        Element biqdm_inpvect ftype=real(dblprec) pytype=float
        
        
        Defined at /home/andersb/CrossPlatform/UppASD/source/Input/inputdatatype.f90 \
            line 90
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _uppasd.f90wrap_ham_inp_t__array__biqdm_inpvect(self._handle)
        if array_handle in self._arrays:
            biqdm_inpvect = self._arrays[array_handle]
        else:
            biqdm_inpvect = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _uppasd.f90wrap_ham_inp_t__array__biqdm_inpvect)
            self._arrays[array_handle] = biqdm_inpvect
        return biqdm_inpvect
    
    @biqdm_inpvect.setter
    def biqdm_inpvect(self, biqdm_inpvect):
        self.biqdm_inpvect[...] = biqdm_inpvect
    
    @property
    def do_bq(self):
        """
        Element do_bq ftype=integer  pytype=int
        
        
        Defined at /home/andersb/CrossPlatform/UppASD/source/Input/inputdatatype.f90 \
            line 94
        
        """
        return _uppasd.f90wrap_ham_inp_t__get__do_bq(self._handle)
    
    @do_bq.setter
    def do_bq(self, do_bq):
        _uppasd.f90wrap_ham_inp_t__set__do_bq(self._handle, do_bq)
    
    @property
    def max_no_bqshells(self):
        """
        Element max_no_bqshells ftype=integer  pytype=int
        
        
        Defined at /home/andersb/CrossPlatform/UppASD/source/Input/inputdatatype.f90 \
            line 95
        
        """
        return _uppasd.f90wrap_ham_inp_t__get__max_no_bqshells(self._handle)
    
    @max_no_bqshells.setter
    def max_no_bqshells(self, max_no_bqshells):
        _uppasd.f90wrap_ham_inp_t__set__max_no_bqshells(self._handle, max_no_bqshells)
    
    @property
    def bqfile(self):
        """
        Element bqfile ftype=character(len=35) pytype=str
        
        
        Defined at /home/andersb/CrossPlatform/UppASD/source/Input/inputdatatype.f90 \
            line 96
        
        """
        return _uppasd.f90wrap_ham_inp_t__get__bqfile(self._handle)
    
    @bqfile.setter
    def bqfile(self, bqfile):
        _uppasd.f90wrap_ham_inp_t__set__bqfile(self._handle, bqfile)
    
    @property
    def bq_nn(self):
        """
        Element bq_nn ftype=integer pytype=int
        
        
        Defined at /home/andersb/CrossPlatform/UppASD/source/Input/inputdatatype.f90 \
            line 97
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _uppasd.f90wrap_ham_inp_t__array__bq_nn(self._handle)
        if array_handle in self._arrays:
            bq_nn = self._arrays[array_handle]
        else:
            bq_nn = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _uppasd.f90wrap_ham_inp_t__array__bq_nn)
            self._arrays[array_handle] = bq_nn
        return bq_nn
    
    @bq_nn.setter
    def bq_nn(self, bq_nn):
        self.bq_nn[...] = bq_nn
    
    @property
    def jc_bq(self):
        """
        Element jc_bq ftype=real(dblprec) pytype=float
        
        
        Defined at /home/andersb/CrossPlatform/UppASD/source/Input/inputdatatype.f90 \
            line 98
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _uppasd.f90wrap_ham_inp_t__array__jc_bq(self._handle)
        if array_handle in self._arrays:
            jc_bq = self._arrays[array_handle]
        else:
            jc_bq = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _uppasd.f90wrap_ham_inp_t__array__jc_bq)
            self._arrays[array_handle] = jc_bq
        return jc_bq
    
    @jc_bq.setter
    def jc_bq(self, jc_bq):
        self.jc_bq[...] = jc_bq
    
    @property
    def bq_redcoord(self):
        """
        Element bq_redcoord ftype=real(dblprec) pytype=float
        
        
        Defined at /home/andersb/CrossPlatform/UppASD/source/Input/inputdatatype.f90 \
            line 99
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _uppasd.f90wrap_ham_inp_t__array__bq_redcoord(self._handle)
        if array_handle in self._arrays:
            bq_redcoord = self._arrays[array_handle]
        else:
            bq_redcoord = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _uppasd.f90wrap_ham_inp_t__array__bq_redcoord)
            self._arrays[array_handle] = bq_redcoord
        return bq_redcoord
    
    @bq_redcoord.setter
    def bq_redcoord(self, bq_redcoord):
        self.bq_redcoord[...] = bq_redcoord
    
    @property
    def do_ring(self):
        """
        Element do_ring ftype=integer  pytype=int
        
        
        Defined at /home/andersb/CrossPlatform/UppASD/source/Input/inputdatatype.f90 \
            line 103
        
        """
        return _uppasd.f90wrap_ham_inp_t__get__do_ring(self._handle)
    
    @do_ring.setter
    def do_ring(self, do_ring):
        _uppasd.f90wrap_ham_inp_t__set__do_ring(self._handle, do_ring)
    
    @property
    def max_no_ringshells(self):
        """
        Element max_no_ringshells ftype=integer  pytype=int
        
        
        Defined at /home/andersb/CrossPlatform/UppASD/source/Input/inputdatatype.f90 \
            line 104
        
        """
        return _uppasd.f90wrap_ham_inp_t__get__max_no_ringshells(self._handle)
    
    @max_no_ringshells.setter
    def max_no_ringshells(self, max_no_ringshells):
        _uppasd.f90wrap_ham_inp_t__set__max_no_ringshells(self._handle, \
            max_no_ringshells)
    
    @property
    def ringfile(self):
        """
        Element ringfile ftype=character(len=35) pytype=str
        
        
        Defined at /home/andersb/CrossPlatform/UppASD/source/Input/inputdatatype.f90 \
            line 105
        
        """
        return _uppasd.f90wrap_ham_inp_t__get__ringfile(self._handle)
    
    @ringfile.setter
    def ringfile(self, ringfile):
        _uppasd.f90wrap_ham_inp_t__set__ringfile(self._handle, ringfile)
    
    @property
    def ring_nn(self):
        """
        Element ring_nn ftype=integer pytype=int
        
        
        Defined at /home/andersb/CrossPlatform/UppASD/source/Input/inputdatatype.f90 \
            line 106
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _uppasd.f90wrap_ham_inp_t__array__ring_nn(self._handle)
        if array_handle in self._arrays:
            ring_nn = self._arrays[array_handle]
        else:
            ring_nn = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _uppasd.f90wrap_ham_inp_t__array__ring_nn)
            self._arrays[array_handle] = ring_nn
        return ring_nn
    
    @ring_nn.setter
    def ring_nn(self, ring_nn):
        self.ring_nn[...] = ring_nn
    
    @property
    def ring_redcoord_ij(self):
        """
        Element ring_redcoord_ij ftype=real(dblprec) pytype=float
        
        
        Defined at /home/andersb/CrossPlatform/UppASD/source/Input/inputdatatype.f90 \
            line 107
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _uppasd.f90wrap_ham_inp_t__array__ring_redcoord_ij(self._handle)
        if array_handle in self._arrays:
            ring_redcoord_ij = self._arrays[array_handle]
        else:
            ring_redcoord_ij = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _uppasd.f90wrap_ham_inp_t__array__ring_redcoord_ij)
            self._arrays[array_handle] = ring_redcoord_ij
        return ring_redcoord_ij
    
    @ring_redcoord_ij.setter
    def ring_redcoord_ij(self, ring_redcoord_ij):
        self.ring_redcoord_ij[...] = ring_redcoord_ij
    
    @property
    def ring_redcoord_ik(self):
        """
        Element ring_redcoord_ik ftype=real(dblprec) pytype=float
        
        
        Defined at /home/andersb/CrossPlatform/UppASD/source/Input/inputdatatype.f90 \
            line 108
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _uppasd.f90wrap_ham_inp_t__array__ring_redcoord_ik(self._handle)
        if array_handle in self._arrays:
            ring_redcoord_ik = self._arrays[array_handle]
        else:
            ring_redcoord_ik = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _uppasd.f90wrap_ham_inp_t__array__ring_redcoord_ik)
            self._arrays[array_handle] = ring_redcoord_ik
        return ring_redcoord_ik
    
    @ring_redcoord_ik.setter
    def ring_redcoord_ik(self, ring_redcoord_ik):
        self.ring_redcoord_ik[...] = ring_redcoord_ik
    
    @property
    def ring_redcoord_il(self):
        """
        Element ring_redcoord_il ftype=real(dblprec) pytype=float
        
        
        Defined at /home/andersb/CrossPlatform/UppASD/source/Input/inputdatatype.f90 \
            line 109
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _uppasd.f90wrap_ham_inp_t__array__ring_redcoord_il(self._handle)
        if array_handle in self._arrays:
            ring_redcoord_il = self._arrays[array_handle]
        else:
            ring_redcoord_il = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _uppasd.f90wrap_ham_inp_t__array__ring_redcoord_il)
            self._arrays[array_handle] = ring_redcoord_il
        return ring_redcoord_il
    
    @ring_redcoord_il.setter
    def ring_redcoord_il(self, ring_redcoord_il):
        self.ring_redcoord_il[...] = ring_redcoord_il
    
    @property
    def jc_ring(self):
        """
        Element jc_ring ftype=real(dblprec) pytype=float
        
        
        Defined at /home/andersb/CrossPlatform/UppASD/source/Input/inputdatatype.f90 \
            line 110
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _uppasd.f90wrap_ham_inp_t__array__jc_ring(self._handle)
        if array_handle in self._arrays:
            jc_ring = self._arrays[array_handle]
        else:
            jc_ring = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _uppasd.f90wrap_ham_inp_t__array__jc_ring)
            self._arrays[array_handle] = jc_ring
        return jc_ring
    
    @jc_ring.setter
    def jc_ring(self, jc_ring):
        self.jc_ring[...] = jc_ring
    
    @property
    def do_dip(self):
        """
        Element do_dip ftype=integer  pytype=int
        
        
        Defined at /home/andersb/CrossPlatform/UppASD/source/Input/inputdatatype.f90 \
            line 114
        
        """
        return _uppasd.f90wrap_ham_inp_t__get__do_dip(self._handle)
    
    @do_dip.setter
    def do_dip(self, do_dip):
        _uppasd.f90wrap_ham_inp_t__set__do_dip(self._handle, do_dip)
    
    @property
    def read_dipole(self):
        """
        Element read_dipole ftype=character(len=1) pytype=str
        
        
        Defined at /home/andersb/CrossPlatform/UppASD/source/Input/inputdatatype.f90 \
            line 115
        
        """
        return _uppasd.f90wrap_ham_inp_t__get__read_dipole(self._handle)
    
    @read_dipole.setter
    def read_dipole(self, read_dipole):
        _uppasd.f90wrap_ham_inp_t__set__read_dipole(self._handle, read_dipole)
    
    @property
    def print_dip_tensor(self):
        """
        Element print_dip_tensor ftype=character(len=1) pytype=str
        
        
        Defined at /home/andersb/CrossPlatform/UppASD/source/Input/inputdatatype.f90 \
            line 116
        
        """
        return _uppasd.f90wrap_ham_inp_t__get__print_dip_tensor(self._handle)
    
    @print_dip_tensor.setter
    def print_dip_tensor(self, print_dip_tensor):
        _uppasd.f90wrap_ham_inp_t__set__print_dip_tensor(self._handle, print_dip_tensor)
    
    @property
    def qdip_files(self):
        """
        Element qdip_files ftype=character(len=30) pytype=str
        
        
        Defined at /home/andersb/CrossPlatform/UppASD/source/Input/inputdatatype.f90 \
            line 117
        
        """
        return _uppasd.f90wrap_ham_inp_t__get__qdip_files(self._handle)
    
    @qdip_files.setter
    def qdip_files(self, qdip_files):
        _uppasd.f90wrap_ham_inp_t__set__qdip_files(self._handle, qdip_files)
    
    @property
    def do_ewald(self):
        """
        Element do_ewald ftype=character(len=1) pytype=str
        
        
        Defined at /home/andersb/CrossPlatform/UppASD/source/Input/inputdatatype.f90 \
            line 121
        
        """
        return _uppasd.f90wrap_ham_inp_t__get__do_ewald(self._handle)
    
    @do_ewald.setter
    def do_ewald(self, do_ewald):
        _uppasd.f90wrap_ham_inp_t__set__do_ewald(self._handle, do_ewald)
    
    @property
    def rmax(self):
        """
        Element rmax ftype=integer pytype=int
        
        
        Defined at /home/andersb/CrossPlatform/UppASD/source/Input/inputdatatype.f90 \
            line 122
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _uppasd.f90wrap_ham_inp_t__array__rmax(self._handle)
        if array_handle in self._arrays:
            rmax = self._arrays[array_handle]
        else:
            rmax = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _uppasd.f90wrap_ham_inp_t__array__rmax)
            self._arrays[array_handle] = rmax
        return rmax
    
    @rmax.setter
    def rmax(self, rmax):
        self.rmax[...] = rmax
    
    @property
    def kmax(self):
        """
        Element kmax ftype=integer pytype=int
        
        
        Defined at /home/andersb/CrossPlatform/UppASD/source/Input/inputdatatype.f90 \
            line 123
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _uppasd.f90wrap_ham_inp_t__array__kmax(self._handle)
        if array_handle in self._arrays:
            kmax = self._arrays[array_handle]
        else:
            kmax = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _uppasd.f90wrap_ham_inp_t__array__kmax)
            self._arrays[array_handle] = kmax
        return kmax
    
    @kmax.setter
    def kmax(self, kmax):
        self.kmax[...] = kmax
    
    @property
    def ewald_alpha(self):
        """
        Element ewald_alpha ftype=real(dblprec) pytype=float
        
        
        Defined at /home/andersb/CrossPlatform/UppASD/source/Input/inputdatatype.f90 \
            line 124
        
        """
        return _uppasd.f90wrap_ham_inp_t__get__ewald_alpha(self._handle)
    
    @ewald_alpha.setter
    def ewald_alpha(self, ewald_alpha):
        _uppasd.f90wrap_ham_inp_t__set__ewald_alpha(self._handle, ewald_alpha)
    
    def __str__(self):
        ret = ['<ham_inp_t>{\n']
        ret.append('    do_jtensor : ')
        ret.append(repr(self.do_jtensor))
        ret.append(',\n    calc_jtensor : ')
        ret.append(repr(self.calc_jtensor))
        ret.append(',\n    map_multiple : ')
        ret.append(repr(self.map_multiple))
        ret.append(',\n    max_no_shells : ')
        ret.append(repr(self.max_no_shells))
        ret.append(',\n    nn : ')
        ret.append(repr(self.nn))
        ret.append(',\n    nntype : ')
        ret.append(repr(self.nntype))
        ret.append(',\n    redcoord : ')
        ret.append(repr(self.redcoord))
        ret.append(',\n    jc : ')
        ret.append(repr(self.jc))
        ret.append(',\n    jcd : ')
        ret.append(repr(self.jcd))
        ret.append(',\n    jc_tens : ')
        ret.append(repr(self.jc_tens))
        ret.append(',\n    exc_inter : ')
        ret.append(repr(self.exc_inter))
        ret.append(',\n    jij_scale : ')
        ret.append(repr(self.jij_scale))
        ret.append(',\n    ea_model : ')
        ret.append(repr(self.ea_model))
        ret.append(',\n    ea_sigma : ')
        ret.append(repr(self.ea_sigma))
        ret.append(',\n    do_anisotropy : ')
        ret.append(repr(self.do_anisotropy))
        ret.append(',\n    random_anisotropy_density : ')
        ret.append(repr(self.random_anisotropy_density))
        ret.append(',\n    random_anisotropy : ')
        ret.append(repr(self.random_anisotropy))
        ret.append(',\n    mult_axis : ')
        ret.append(repr(self.mult_axis))
        ret.append(',\n    kfile : ')
        ret.append(repr(self.kfile))
        ret.append(',\n    anisotropytype : ')
        ret.append(repr(self.anisotropytype))
        ret.append(',\n    anisotropytype_diff : ')
        ret.append(repr(self.anisotropytype_diff))
        ret.append(',\n    anisotropy : ')
        ret.append(repr(self.anisotropy))
        ret.append(',\n    anisotropy_diff : ')
        ret.append(repr(self.anisotropy_diff))
        ret.append(',\n    do_dm : ')
        ret.append(repr(self.do_dm))
        ret.append(',\n    max_no_dmshells : ')
        ret.append(repr(self.max_no_dmshells))
        ret.append(',\n    dmfile : ')
        ret.append(repr(self.dmfile))
        ret.append(',\n    dm_scale : ')
        ret.append(repr(self.dm_scale))
        ret.append(',\n    dm_nn : ')
        ret.append(repr(self.dm_nn))
        ret.append(',\n    dm_redcoord : ')
        ret.append(repr(self.dm_redcoord))
        ret.append(',\n    dm_inpvect : ')
        ret.append(repr(self.dm_inpvect))
        ret.append(',\n    do_chir : ')
        ret.append(repr(self.do_chir))
        ret.append(',\n    max_no_chirshells : ')
        ret.append(repr(self.max_no_chirshells))
        ret.append(',\n    chirfile : ')
        ret.append(repr(self.chirfile))
        ret.append(',\n    chir_nn : ')
        ret.append(repr(self.chir_nn))
        ret.append(',\n    chir_inpval : ')
        ret.append(repr(self.chir_inpval))
        ret.append(',\n    chir_redcoord : ')
        ret.append(repr(self.chir_redcoord))
        ret.append(',\n    do_fourx : ')
        ret.append(repr(self.do_fourx))
        ret.append(',\n    max_no_fourxshells : ')
        ret.append(repr(self.max_no_fourxshells))
        ret.append(',\n    fourxfile : ')
        ret.append(repr(self.fourxfile))
        ret.append(',\n    fourx_nn : ')
        ret.append(repr(self.fourx_nn))
        ret.append(',\n    fourx_inpval : ')
        ret.append(repr(self.fourx_inpval))
        ret.append(',\n    fourx_redcoord : ')
        ret.append(repr(self.fourx_redcoord))
        ret.append(',\n    do_pd : ')
        ret.append(repr(self.do_pd))
        ret.append(',\n    max_no_pdshells : ')
        ret.append(repr(self.max_no_pdshells))
        ret.append(',\n    pdfile : ')
        ret.append(repr(self.pdfile))
        ret.append(',\n    pd_nn : ')
        ret.append(repr(self.pd_nn))
        ret.append(',\n    pd_redcoord : ')
        ret.append(repr(self.pd_redcoord))
        ret.append(',\n    pd_inpvect : ')
        ret.append(repr(self.pd_inpvect))
        ret.append(',\n    do_biqdm : ')
        ret.append(repr(self.do_biqdm))
        ret.append(',\n    max_no_biqdmshells : ')
        ret.append(repr(self.max_no_biqdmshells))
        ret.append(',\n    biqdmfile : ')
        ret.append(repr(self.biqdmfile))
        ret.append(',\n    biqdm_nn : ')
        ret.append(repr(self.biqdm_nn))
        ret.append(',\n    biqdm_redcoord : ')
        ret.append(repr(self.biqdm_redcoord))
        ret.append(',\n    biqdm_inpvect : ')
        ret.append(repr(self.biqdm_inpvect))
        ret.append(',\n    do_bq : ')
        ret.append(repr(self.do_bq))
        ret.append(',\n    max_no_bqshells : ')
        ret.append(repr(self.max_no_bqshells))
        ret.append(',\n    bqfile : ')
        ret.append(repr(self.bqfile))
        ret.append(',\n    bq_nn : ')
        ret.append(repr(self.bq_nn))
        ret.append(',\n    jc_bq : ')
        ret.append(repr(self.jc_bq))
        ret.append(',\n    bq_redcoord : ')
        ret.append(repr(self.bq_redcoord))
        ret.append(',\n    do_ring : ')
        ret.append(repr(self.do_ring))
        ret.append(',\n    max_no_ringshells : ')
        ret.append(repr(self.max_no_ringshells))
        ret.append(',\n    ringfile : ')
        ret.append(repr(self.ringfile))
        ret.append(',\n    ring_nn : ')
        ret.append(repr(self.ring_nn))
        ret.append(',\n    ring_redcoord_ij : ')
        ret.append(repr(self.ring_redcoord_ij))
        ret.append(',\n    ring_redcoord_ik : ')
        ret.append(repr(self.ring_redcoord_ik))
        ret.append(',\n    ring_redcoord_il : ')
        ret.append(repr(self.ring_redcoord_il))
        ret.append(',\n    jc_ring : ')
        ret.append(repr(self.jc_ring))
        ret.append(',\n    do_dip : ')
        ret.append(repr(self.do_dip))
        ret.append(',\n    read_dipole : ')
        ret.append(repr(self.read_dipole))
        ret.append(',\n    print_dip_tensor : ')
        ret.append(repr(self.print_dip_tensor))
        ret.append(',\n    qdip_files : ')
        ret.append(repr(self.qdip_files))
        ret.append(',\n    do_ewald : ')
        ret.append(repr(self.do_ewald))
        ret.append(',\n    rmax : ')
        ret.append(repr(self.rmax))
        ret.append(',\n    kmax : ')
        ret.append(repr(self.kmax))
        ret.append(',\n    ewald_alpha : ')
        ret.append(repr(self.ewald_alpha))
        ret.append('}')
        return ''.join(ret)
    
    _dt_array_initialisers = []
    


_array_initialisers = []
_dt_array_initialisers = []

try:
    for func in _array_initialisers:
        func()
except ValueError:
    logging.debug('unallocated array(s) detected on import of module \
        "inputdatatype".')

for func in _dt_array_initialisers:
    func()
