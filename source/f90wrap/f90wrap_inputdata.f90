! Module inputdata defined in file Input/inputdata.f90

subroutine f90wrap_set_input_defaults
    use inputdata, only: set_input_defaults
    implicit none
    
    call set_input_defaults()
end subroutine f90wrap_set_input_defaults

subroutine f90wrap_allocate_initmag(na, nchmax, conf_num, flag)
    use inputdata, only: allocate_initmag
    implicit none
    
    integer, optional, intent(in) :: na
    integer, optional, intent(in) :: nchmax
    integer, optional, intent(in) :: conf_num
    integer, intent(in) :: flag
    call allocate_initmag(NA=na, Nchmax=nchmax, conf_num=conf_num, flag=flag)
end subroutine f90wrap_allocate_initmag

subroutine f90wrap_allocate_chemicalinput(na, nchmax, flag)
    use inputdata, only: allocate_chemicalinput
    implicit none
    
    integer, optional, intent(in) :: na
    integer, optional, intent(in) :: nchmax
    integer, intent(in) :: flag
    call allocate_chemicalinput(NA=na, Nchmax=nchmax, flag=flag)
end subroutine f90wrap_allocate_chemicalinput

subroutine f90wrap_reshape_hamiltonianinput
    use inputdata, only: reshape_hamiltonianinput
    implicit none
    
    call reshape_hamiltonianinput()
end subroutine f90wrap_reshape_hamiltonianinput

subroutine f90wrap_allocate_hamiltonianinput(ham_inp, no_shells, flag)
    use inputdata, only: allocate_hamiltonianinput
    use inputdatatype, only: ham_inp_t
    implicit none
    
    type ham_inp_t_ptr_type
        type(ham_inp_t), pointer :: p => NULL()
    end type ham_inp_t_ptr_type
    type(ham_inp_t_ptr_type) :: ham_inp_ptr
    integer, intent(in), dimension(2) :: ham_inp
    integer, intent(in), optional :: no_shells
    integer, intent(in) :: flag
    ham_inp_ptr = transfer(ham_inp, ham_inp_ptr)
    call allocate_hamiltonianinput(ham_inp=ham_inp_ptr%p, no_shells=no_shells, flag=flag)
end subroutine f90wrap_allocate_hamiltonianinput

subroutine f90wrap_allocate_nn(i)
    use inputdata, only: allocate_nn
    implicit none
    
    integer, intent(in) :: i
    call allocate_nn(i=i)
end subroutine f90wrap_allocate_nn

subroutine f90wrap_inputdata__get__N1(f90wrap_N1)
    use inputdata, only: inputdata_N1 => N1
    implicit none
    integer, intent(out) :: f90wrap_N1
    
    f90wrap_N1 = inputdata_N1
end subroutine f90wrap_inputdata__get__N1

subroutine f90wrap_inputdata__set__N1(f90wrap_N1)
    use inputdata, only: inputdata_N1 => N1
    implicit none
    integer, intent(in) :: f90wrap_N1
    
    inputdata_N1 = f90wrap_N1
end subroutine f90wrap_inputdata__set__N1

subroutine f90wrap_inputdata__get__N2(f90wrap_N2)
    use inputdata, only: inputdata_N2 => N2
    implicit none
    integer, intent(out) :: f90wrap_N2
    
    f90wrap_N2 = inputdata_N2
end subroutine f90wrap_inputdata__get__N2

subroutine f90wrap_inputdata__set__N2(f90wrap_N2)
    use inputdata, only: inputdata_N2 => N2
    implicit none
    integer, intent(in) :: f90wrap_N2
    
    inputdata_N2 = f90wrap_N2
end subroutine f90wrap_inputdata__set__N2

subroutine f90wrap_inputdata__get__N3(f90wrap_N3)
    use inputdata, only: inputdata_N3 => N3
    implicit none
    integer, intent(out) :: f90wrap_N3
    
    f90wrap_N3 = inputdata_N3
end subroutine f90wrap_inputdata__get__N3

subroutine f90wrap_inputdata__set__N3(f90wrap_N3)
    use inputdata, only: inputdata_N3 => N3
    implicit none
    integer, intent(in) :: f90wrap_N3
    
    inputdata_N3 = f90wrap_N3
end subroutine f90wrap_inputdata__set__N3

subroutine f90wrap_inputdata__get__NA(f90wrap_NA)
    use inputdata, only: inputdata_NA => NA
    implicit none
    integer, intent(out) :: f90wrap_NA
    
    f90wrap_NA = inputdata_NA
end subroutine f90wrap_inputdata__get__NA

subroutine f90wrap_inputdata__set__NA(f90wrap_NA)
    use inputdata, only: inputdata_NA => NA
    implicit none
    integer, intent(in) :: f90wrap_NA
    
    inputdata_NA = f90wrap_NA
end subroutine f90wrap_inputdata__set__NA

subroutine f90wrap_inputdata__get__NT(f90wrap_NT)
    use inputdata, only: inputdata_NT => NT
    implicit none
    integer, intent(out) :: f90wrap_NT
    
    f90wrap_NT = inputdata_NT
end subroutine f90wrap_inputdata__get__NT

subroutine f90wrap_inputdata__set__NT(f90wrap_NT)
    use inputdata, only: inputdata_NT => NT
    implicit none
    integer, intent(in) :: f90wrap_NT
    
    inputdata_NT = f90wrap_NT
end subroutine f90wrap_inputdata__set__NT

subroutine f90wrap_inputdata__get__Sym(f90wrap_Sym)
    use inputdata, only: inputdata_Sym => Sym
    implicit none
    integer, intent(out) :: f90wrap_Sym
    
    f90wrap_Sym = inputdata_Sym
end subroutine f90wrap_inputdata__get__Sym

subroutine f90wrap_inputdata__set__Sym(f90wrap_Sym)
    use inputdata, only: inputdata_Sym => Sym
    implicit none
    integer, intent(in) :: f90wrap_Sym
    
    inputdata_Sym = f90wrap_Sym
end subroutine f90wrap_inputdata__set__Sym

subroutine f90wrap_inputdata__get__nHam(f90wrap_nHam)
    use inputdata, only: inputdata_nHam => nHam
    implicit none
    integer, intent(out) :: f90wrap_nHam
    
    f90wrap_nHam = inputdata_nHam
end subroutine f90wrap_inputdata__get__nHam

subroutine f90wrap_inputdata__set__nHam(f90wrap_nHam)
    use inputdata, only: inputdata_nHam => nHam
    implicit none
    integer, intent(in) :: f90wrap_nHam
    
    inputdata_nHam = f90wrap_nHam
end subroutine f90wrap_inputdata__set__nHam

subroutine f90wrap_inputdata__get__Natom(f90wrap_Natom)
    use inputdata, only: inputdata_Natom => Natom
    implicit none
    integer, intent(out) :: f90wrap_Natom
    
    f90wrap_Natom = inputdata_Natom
end subroutine f90wrap_inputdata__get__Natom

subroutine f90wrap_inputdata__set__Natom(f90wrap_Natom)
    use inputdata, only: inputdata_Natom => Natom
    implicit none
    integer, intent(in) :: f90wrap_Natom
    
    inputdata_Natom = f90wrap_Natom
end subroutine f90wrap_inputdata__set__Natom

subroutine f90wrap_inputdata__get__Natom_full(f90wrap_Natom_full)
    use inputdata, only: inputdata_Natom_full => Natom_full
    implicit none
    integer, intent(out) :: f90wrap_Natom_full
    
    f90wrap_Natom_full = inputdata_Natom_full
end subroutine f90wrap_inputdata__get__Natom_full

subroutine f90wrap_inputdata__set__Natom_full(f90wrap_Natom_full)
    use inputdata, only: inputdata_Natom_full => Natom_full
    implicit none
    integer, intent(in) :: f90wrap_Natom_full
    
    inputdata_Natom_full = f90wrap_Natom_full
end subroutine f90wrap_inputdata__set__Natom_full

subroutine f90wrap_inputdata__get__set_landeg(f90wrap_set_landeg)
    use inputdata, only: inputdata_set_landeg => set_landeg
    implicit none
    integer, intent(out) :: f90wrap_set_landeg
    
    f90wrap_set_landeg = inputdata_set_landeg
end subroutine f90wrap_inputdata__get__set_landeg

subroutine f90wrap_inputdata__set__set_landeg(f90wrap_set_landeg)
    use inputdata, only: inputdata_set_landeg => set_landeg
    implicit none
    integer, intent(in) :: f90wrap_set_landeg
    
    inputdata_set_landeg = f90wrap_set_landeg
end subroutine f90wrap_inputdata__set__set_landeg

subroutine f90wrap_inputdata__get__block_size(f90wrap_block_size)
    use inputdata, only: inputdata_block_size => block_size
    implicit none
    integer, intent(out) :: f90wrap_block_size
    
    f90wrap_block_size = inputdata_block_size
end subroutine f90wrap_inputdata__get__block_size

subroutine f90wrap_inputdata__set__block_size(f90wrap_block_size)
    use inputdata, only: inputdata_block_size => block_size
    implicit none
    integer, intent(in) :: f90wrap_block_size
    
    inputdata_block_size = f90wrap_block_size
end subroutine f90wrap_inputdata__set__block_size

subroutine f90wrap_inputdata__get__metatype(f90wrap_metatype)
    use inputdata, only: inputdata_metatype => metatype
    implicit none
    integer, intent(out) :: f90wrap_metatype
    
    f90wrap_metatype = inputdata_metatype
end subroutine f90wrap_inputdata__get__metatype

subroutine f90wrap_inputdata__set__metatype(f90wrap_metatype)
    use inputdata, only: inputdata_metatype => metatype
    implicit none
    integer, intent(in) :: f90wrap_metatype
    
    inputdata_metatype = f90wrap_metatype
end subroutine f90wrap_inputdata__set__metatype

subroutine f90wrap_inputdata__get__metanumb(f90wrap_metanumb)
    use inputdata, only: inputdata_metanumb => metanumb
    implicit none
    integer, intent(out) :: f90wrap_metanumb
    
    f90wrap_metanumb = inputdata_metanumb
end subroutine f90wrap_inputdata__get__metanumb

subroutine f90wrap_inputdata__set__metanumb(f90wrap_metanumb)
    use inputdata, only: inputdata_metanumb => metanumb
    implicit none
    integer, intent(in) :: f90wrap_metanumb
    
    inputdata_metanumb = f90wrap_metanumb
end subroutine f90wrap_inputdata__set__metanumb

subroutine f90wrap_inputdata__array__acomp(dummy_this, nd, dtype, dshape, dloc)
    use parameters
    use profiling
    use inputdatatype
    use inputdata, only: inputdata_acomp => acomp
    use, intrinsic :: iso_c_binding, only : c_int
    implicit none
    integer, intent(in) :: dummy_this(2)
    integer(c_int), intent(out) :: nd
    integer(c_int), intent(out) :: dtype
    integer(c_int), dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 1
    dtype = 5
    if (allocated(inputdata_acomp)) then
        dshape(1:1) = shape(inputdata_acomp)
        dloc = loc(inputdata_acomp)
    else
        dloc = 0
    end if
end subroutine f90wrap_inputdata__array__acomp

subroutine f90wrap_inputdata__array__asite(dummy_this, nd, dtype, dshape, dloc)
    use parameters
    use profiling
    use inputdatatype
    use inputdata, only: inputdata_asite => asite
    use, intrinsic :: iso_c_binding, only : c_int
    implicit none
    integer, intent(in) :: dummy_this(2)
    integer(c_int), intent(out) :: nd
    integer(c_int), intent(out) :: dtype
    integer(c_int), dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 1
    dtype = 5
    if (allocated(inputdata_asite)) then
        dshape(1:1) = shape(inputdata_asite)
        dloc = loc(inputdata_asite)
    else
        dloc = 0
    end if
end subroutine f90wrap_inputdata__array__asite

subroutine f90wrap_inputdata__array__anumb_inp(dummy_this, nd, dtype, dshape, dloc)
    use parameters
    use profiling
    use inputdatatype
    use inputdata, only: inputdata_anumb_inp => anumb_inp
    use, intrinsic :: iso_c_binding, only : c_int
    implicit none
    integer, intent(in) :: dummy_this(2)
    integer(c_int), intent(out) :: nd
    integer(c_int), intent(out) :: dtype
    integer(c_int), dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 1
    dtype = 5
    if (allocated(inputdata_anumb_inp)) then
        dshape(1:1) = shape(inputdata_anumb_inp)
        dloc = loc(inputdata_anumb_inp)
    else
        dloc = 0
    end if
end subroutine f90wrap_inputdata__array__anumb_inp

subroutine f90wrap_inputdata__array__atype_inp(dummy_this, nd, dtype, dshape, dloc)
    use parameters
    use profiling
    use inputdatatype
    use inputdata, only: inputdata_atype_inp => atype_inp
    use, intrinsic :: iso_c_binding, only : c_int
    implicit none
    integer, intent(in) :: dummy_this(2)
    integer(c_int), intent(out) :: nd
    integer(c_int), intent(out) :: dtype
    integer(c_int), dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 1
    dtype = 5
    if (allocated(inputdata_atype_inp)) then
        dshape(1:1) = shape(inputdata_atype_inp)
        dloc = loc(inputdata_atype_inp)
    else
        dloc = 0
    end if
end subroutine f90wrap_inputdata__array__atype_inp

subroutine f90wrap_inputdata__get__posfiletype(f90wrap_posfiletype)
    use inputdata, only: inputdata_posfiletype => posfiletype
    implicit none
    character(1), intent(out) :: f90wrap_posfiletype
    
    f90wrap_posfiletype = inputdata_posfiletype
end subroutine f90wrap_inputdata__get__posfiletype

subroutine f90wrap_inputdata__set__posfiletype(f90wrap_posfiletype)
    use inputdata, only: inputdata_posfiletype => posfiletype
    implicit none
    character(1), intent(in) :: f90wrap_posfiletype
    
    inputdata_posfiletype = f90wrap_posfiletype
end subroutine f90wrap_inputdata__set__posfiletype

subroutine f90wrap_inputdata__get__do_sortcoup(f90wrap_do_sortcoup)
    use inputdata, only: inputdata_do_sortcoup => do_sortcoup
    implicit none
    character(1), intent(out) :: f90wrap_do_sortcoup
    
    f90wrap_do_sortcoup = inputdata_do_sortcoup
end subroutine f90wrap_inputdata__get__do_sortcoup

subroutine f90wrap_inputdata__set__do_sortcoup(f90wrap_do_sortcoup)
    use inputdata, only: inputdata_do_sortcoup => do_sortcoup
    implicit none
    character(1), intent(in) :: f90wrap_do_sortcoup
    
    inputdata_do_sortcoup = f90wrap_do_sortcoup
end subroutine f90wrap_inputdata__set__do_sortcoup

subroutine f90wrap_inputdata__get__BC1(f90wrap_BC1)
    use inputdata, only: inputdata_BC1 => BC1
    implicit none
    character(1), intent(out) :: f90wrap_BC1
    
    f90wrap_BC1 = inputdata_BC1
end subroutine f90wrap_inputdata__get__BC1

subroutine f90wrap_inputdata__set__BC1(f90wrap_BC1)
    use inputdata, only: inputdata_BC1 => BC1
    implicit none
    character(1), intent(in) :: f90wrap_BC1
    
    inputdata_BC1 = f90wrap_BC1
end subroutine f90wrap_inputdata__set__BC1

subroutine f90wrap_inputdata__get__BC2(f90wrap_BC2)
    use inputdata, only: inputdata_BC2 => BC2
    implicit none
    character(1), intent(out) :: f90wrap_BC2
    
    f90wrap_BC2 = inputdata_BC2
end subroutine f90wrap_inputdata__get__BC2

subroutine f90wrap_inputdata__set__BC2(f90wrap_BC2)
    use inputdata, only: inputdata_BC2 => BC2
    implicit none
    character(1), intent(in) :: f90wrap_BC2
    
    inputdata_BC2 = f90wrap_BC2
end subroutine f90wrap_inputdata__set__BC2

subroutine f90wrap_inputdata__get__BC3(f90wrap_BC3)
    use inputdata, only: inputdata_BC3 => BC3
    implicit none
    character(1), intent(out) :: f90wrap_BC3
    
    f90wrap_BC3 = inputdata_BC3
end subroutine f90wrap_inputdata__get__BC3

subroutine f90wrap_inputdata__set__BC3(f90wrap_BC3)
    use inputdata, only: inputdata_BC3 => BC3
    implicit none
    character(1), intent(in) :: f90wrap_BC3
    
    inputdata_BC3 = f90wrap_BC3
end subroutine f90wrap_inputdata__set__BC3

subroutine f90wrap_inputdata__get__posfile(f90wrap_posfile)
    use inputdata, only: inputdata_posfile => posfile
    implicit none
    character(35), intent(out) :: f90wrap_posfile
    
    f90wrap_posfile = inputdata_posfile
end subroutine f90wrap_inputdata__get__posfile

subroutine f90wrap_inputdata__set__posfile(f90wrap_posfile)
    use inputdata, only: inputdata_posfile => posfile
    implicit none
    character(35), intent(in) :: f90wrap_posfile
    
    inputdata_posfile = f90wrap_posfile
end subroutine f90wrap_inputdata__set__posfile

subroutine f90wrap_inputdata__get__alat(f90wrap_alat)
    use inputdata, only: inputdata_alat => alat
    implicit none
    real(8), intent(out) :: f90wrap_alat
    
    f90wrap_alat = inputdata_alat
end subroutine f90wrap_inputdata__get__alat

subroutine f90wrap_inputdata__set__alat(f90wrap_alat)
    use inputdata, only: inputdata_alat => alat
    implicit none
    real(8), intent(in) :: f90wrap_alat
    
    inputdata_alat = f90wrap_alat
end subroutine f90wrap_inputdata__set__alat

subroutine f90wrap_inputdata__get__scalefac(f90wrap_scalefac)
    use inputdata, only: inputdata_scalefac => scalefac
    implicit none
    real(8), intent(out) :: f90wrap_scalefac
    
    f90wrap_scalefac = inputdata_scalefac
end subroutine f90wrap_inputdata__get__scalefac

subroutine f90wrap_inputdata__set__scalefac(f90wrap_scalefac)
    use inputdata, only: inputdata_scalefac => scalefac
    implicit none
    real(8), intent(in) :: f90wrap_scalefac
    
    inputdata_scalefac = f90wrap_scalefac
end subroutine f90wrap_inputdata__set__scalefac

subroutine f90wrap_inputdata__get__Landeg_glob(f90wrap_Landeg_glob)
    use inputdata, only: inputdata_Landeg_glob => Landeg_glob
    implicit none
    real(8), intent(out) :: f90wrap_Landeg_glob
    
    f90wrap_Landeg_glob = inputdata_Landeg_glob
end subroutine f90wrap_inputdata__get__Landeg_glob

subroutine f90wrap_inputdata__set__Landeg_glob(f90wrap_Landeg_glob)
    use inputdata, only: inputdata_Landeg_glob => Landeg_glob
    implicit none
    real(8), intent(in) :: f90wrap_Landeg_glob
    
    inputdata_Landeg_glob = f90wrap_Landeg_glob
end subroutine f90wrap_inputdata__set__Landeg_glob

subroutine f90wrap_inputdata__array__C1(dummy_this, nd, dtype, dshape, dloc)
    use parameters
    use profiling
    use inputdatatype
    use inputdata, only: inputdata_c1 => c1
    use, intrinsic :: iso_c_binding, only : c_int
    implicit none
    integer, intent(in) :: dummy_this(2)
    integer(c_int), intent(out) :: nd
    integer(c_int), intent(out) :: dtype
    integer(c_int), dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 1
    dtype = 12
    dshape(1:1) = shape(inputdata_C1)
    dloc = loc(inputdata_C1)
end subroutine f90wrap_inputdata__array__C1

subroutine f90wrap_inputdata__array__C2(dummy_this, nd, dtype, dshape, dloc)
    use parameters
    use profiling
    use inputdatatype
    use inputdata, only: inputdata_c2 => c2
    use, intrinsic :: iso_c_binding, only : c_int
    implicit none
    integer, intent(in) :: dummy_this(2)
    integer(c_int), intent(out) :: nd
    integer(c_int), intent(out) :: dtype
    integer(c_int), dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 1
    dtype = 12
    dshape(1:1) = shape(inputdata_C2)
    dloc = loc(inputdata_C2)
end subroutine f90wrap_inputdata__array__C2

subroutine f90wrap_inputdata__array__C3(dummy_this, nd, dtype, dshape, dloc)
    use parameters
    use profiling
    use inputdatatype
    use inputdata, only: inputdata_c3 => c3
    use, intrinsic :: iso_c_binding, only : c_int
    implicit none
    integer, intent(in) :: dummy_this(2)
    integer(c_int), intent(out) :: nd
    integer(c_int), intent(out) :: dtype
    integer(c_int), dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 1
    dtype = 12
    dshape(1:1) = shape(inputdata_C3)
    dloc = loc(inputdata_C3)
end subroutine f90wrap_inputdata__array__C3

subroutine f90wrap_inputdata__array__Bas(dummy_this, nd, dtype, dshape, dloc)
    use parameters
    use profiling
    use inputdatatype
    use inputdata, only: inputdata_bas => bas
    use, intrinsic :: iso_c_binding, only : c_int
    implicit none
    integer, intent(in) :: dummy_this(2)
    integer(c_int), intent(out) :: nd
    integer(c_int), intent(out) :: dtype
    integer(c_int), dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 2
    dtype = 12
    if (allocated(inputdata_Bas)) then
        dshape(1:2) = shape(inputdata_Bas)
        dloc = loc(inputdata_Bas)
    else
        dloc = 0
    end if
end subroutine f90wrap_inputdata__array__Bas

subroutine f90wrap_inputdata__array__Bas0(dummy_this, nd, dtype, dshape, dloc)
    use parameters
    use profiling
    use inputdatatype
    use inputdata, only: inputdata_bas0 => bas0
    use, intrinsic :: iso_c_binding, only : c_int
    implicit none
    integer, intent(in) :: dummy_this(2)
    integer(c_int), intent(out) :: nd
    integer(c_int), intent(out) :: dtype
    integer(c_int), dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 2
    dtype = 12
    if (allocated(inputdata_Bas0)) then
        dshape(1:2) = shape(inputdata_Bas0)
        dloc = loc(inputdata_Bas0)
    else
        dloc = 0
    end if
end subroutine f90wrap_inputdata__array__Bas0

subroutine f90wrap_inputdata__array__Landeg_ch(dummy_this, nd, dtype, dshape, dloc)
    use parameters
    use profiling
    use inputdatatype
    use inputdata, only: inputdata_landeg_ch => landeg_ch
    use, intrinsic :: iso_c_binding, only : c_int
    implicit none
    integer, intent(in) :: dummy_this(2)
    integer(c_int), intent(out) :: nd
    integer(c_int), intent(out) :: dtype
    integer(c_int), dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 3
    dtype = 12
    if (allocated(inputdata_Landeg_ch)) then
        dshape(1:3) = shape(inputdata_Landeg_ch)
        dloc = loc(inputdata_Landeg_ch)
    else
        dloc = 0
    end if
end subroutine f90wrap_inputdata__array__Landeg_ch

subroutine f90wrap_inputdata__get__do_mom_legacy(f90wrap_do_mom_legacy)
    use inputdata, only: inputdata_do_mom_legacy => do_mom_legacy
    implicit none
    character(1), intent(out) :: f90wrap_do_mom_legacy
    
    f90wrap_do_mom_legacy = inputdata_do_mom_legacy
end subroutine f90wrap_inputdata__get__do_mom_legacy

subroutine f90wrap_inputdata__set__do_mom_legacy(f90wrap_do_mom_legacy)
    use inputdata, only: inputdata_do_mom_legacy => do_mom_legacy
    implicit none
    character(1), intent(in) :: f90wrap_do_mom_legacy
    
    inputdata_do_mom_legacy = f90wrap_do_mom_legacy
end subroutine f90wrap_inputdata__set__do_mom_legacy

subroutine f90wrap_inputdata__get__renorm_coll(f90wrap_renorm_coll)
    use inputdata, only: inputdata_renorm_coll => renorm_coll
    implicit none
    character(1), intent(out) :: f90wrap_renorm_coll
    
    f90wrap_renorm_coll = inputdata_renorm_coll
end subroutine f90wrap_inputdata__get__renorm_coll

subroutine f90wrap_inputdata__set__renorm_coll(f90wrap_renorm_coll)
    use inputdata, only: inputdata_renorm_coll => renorm_coll
    implicit none
    character(1), intent(in) :: f90wrap_renorm_coll
    
    inputdata_renorm_coll = f90wrap_renorm_coll
end subroutine f90wrap_inputdata__set__renorm_coll

subroutine f90wrap_inputdata__get__ind_mom_flag(f90wrap_ind_mom_flag)
    use inputdata, only: inputdata_ind_mom_flag => ind_mom_flag
    implicit none
    character(1), intent(out) :: f90wrap_ind_mom_flag
    
    f90wrap_ind_mom_flag = inputdata_ind_mom_flag
end subroutine f90wrap_inputdata__get__ind_mom_flag

subroutine f90wrap_inputdata__set__ind_mom_flag(f90wrap_ind_mom_flag)
    use inputdata, only: inputdata_ind_mom_flag => ind_mom_flag
    implicit none
    character(1), intent(in) :: f90wrap_ind_mom_flag
    
    inputdata_ind_mom_flag = f90wrap_ind_mom_flag
end subroutine f90wrap_inputdata__set__ind_mom_flag

subroutine f90wrap_inputdata__get__ind_mom_type(f90wrap_ind_mom_type)
    use inputdata, only: inputdata_ind_mom_type => ind_mom_type
    implicit none
    integer, intent(out) :: f90wrap_ind_mom_type
    
    f90wrap_ind_mom_type = inputdata_ind_mom_type
end subroutine f90wrap_inputdata__get__ind_mom_type

subroutine f90wrap_inputdata__set__ind_mom_type(f90wrap_ind_mom_type)
    use inputdata, only: inputdata_ind_mom_type => ind_mom_type
    implicit none
    integer, intent(in) :: f90wrap_ind_mom_type
    
    inputdata_ind_mom_type = f90wrap_ind_mom_type
end subroutine f90wrap_inputdata__set__ind_mom_type

subroutine f90wrap_inputdata__get__momfile(f90wrap_momfile)
    use inputdata, only: inputdata_momfile => momfile
    implicit none
    character(35), intent(out) :: f90wrap_momfile
    
    f90wrap_momfile = inputdata_momfile
end subroutine f90wrap_inputdata__get__momfile

subroutine f90wrap_inputdata__set__momfile(f90wrap_momfile)
    use inputdata, only: inputdata_momfile => momfile
    implicit none
    character(35), intent(in) :: f90wrap_momfile
    
    inputdata_momfile = f90wrap_momfile
end subroutine f90wrap_inputdata__set__momfile

subroutine f90wrap_inputdata__get__momfile_i(f90wrap_momfile_i)
    use inputdata, only: inputdata_momfile_i => momfile_i
    implicit none
    character(35), intent(out) :: f90wrap_momfile_i
    
    f90wrap_momfile_i = inputdata_momfile_i
end subroutine f90wrap_inputdata__get__momfile_i

subroutine f90wrap_inputdata__set__momfile_i(f90wrap_momfile_i)
    use inputdata, only: inputdata_momfile_i => momfile_i
    implicit none
    character(35), intent(in) :: f90wrap_momfile_i
    
    inputdata_momfile_i = f90wrap_momfile_i
end subroutine f90wrap_inputdata__set__momfile_i

subroutine f90wrap_inputdata__get__momfile_f(f90wrap_momfile_f)
    use inputdata, only: inputdata_momfile_f => momfile_f
    implicit none
    character(35), intent(out) :: f90wrap_momfile_f
    
    f90wrap_momfile_f = inputdata_momfile_f
end subroutine f90wrap_inputdata__get__momfile_f

subroutine f90wrap_inputdata__set__momfile_f(f90wrap_momfile_f)
    use inputdata, only: inputdata_momfile_f => momfile_f
    implicit none
    character(35), intent(in) :: f90wrap_momfile_f
    
    inputdata_momfile_f = f90wrap_momfile_f
end subroutine f90wrap_inputdata__set__momfile_f

subroutine f90wrap_inputdata__get__ind_tol(f90wrap_ind_tol)
    use inputdata, only: inputdata_ind_tol => ind_tol
    implicit none
    real(8), intent(out) :: f90wrap_ind_tol
    
    f90wrap_ind_tol = inputdata_ind_tol
end subroutine f90wrap_inputdata__get__ind_tol

subroutine f90wrap_inputdata__set__ind_tol(f90wrap_ind_tol)
    use inputdata, only: inputdata_ind_tol => ind_tol
    implicit none
    real(8), intent(in) :: f90wrap_ind_tol
    
    inputdata_ind_tol = f90wrap_ind_tol
end subroutine f90wrap_inputdata__set__ind_tol

subroutine f90wrap_inputdata__get__amp_rnd(f90wrap_amp_rnd)
    use inputdata, only: inputdata_amp_rnd => amp_rnd
    implicit none
    real(8), intent(out) :: f90wrap_amp_rnd
    
    f90wrap_amp_rnd = inputdata_amp_rnd
end subroutine f90wrap_inputdata__get__amp_rnd

subroutine f90wrap_inputdata__set__amp_rnd(f90wrap_amp_rnd)
    use inputdata, only: inputdata_amp_rnd => amp_rnd
    implicit none
    real(8), intent(in) :: f90wrap_amp_rnd
    
    inputdata_amp_rnd = f90wrap_amp_rnd
end subroutine f90wrap_inputdata__set__amp_rnd

subroutine f90wrap_inputdata__get__amp_rnd_path(f90wrap_amp_rnd_path)
    use inputdata, only: inputdata_amp_rnd_path => amp_rnd_path
    implicit none
    real(8), intent(out) :: f90wrap_amp_rnd_path
    
    f90wrap_amp_rnd_path = inputdata_amp_rnd_path
end subroutine f90wrap_inputdata__get__amp_rnd_path

subroutine f90wrap_inputdata__set__amp_rnd_path(f90wrap_amp_rnd_path)
    use inputdata, only: inputdata_amp_rnd_path => amp_rnd_path
    implicit none
    real(8), intent(in) :: f90wrap_amp_rnd_path
    
    inputdata_amp_rnd_path = f90wrap_amp_rnd_path
end subroutine f90wrap_inputdata__set__amp_rnd_path

subroutine f90wrap_inputdata__array__ind_mom(dummy_this, nd, dtype, dshape, dloc)
    use parameters
    use profiling
    use inputdatatype
    use inputdata, only: inputdata_ind_mom => ind_mom
    use, intrinsic :: iso_c_binding, only : c_int
    implicit none
    integer, intent(in) :: dummy_this(2)
    integer(c_int), intent(out) :: nd
    integer(c_int), intent(out) :: dtype
    integer(c_int), dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 2
    dtype = 5
    if (allocated(inputdata_ind_mom)) then
        dshape(1:2) = shape(inputdata_ind_mom)
        dloc = loc(inputdata_ind_mom)
    else
        dloc = 0
    end if
end subroutine f90wrap_inputdata__array__ind_mom

subroutine f90wrap_inputdata__array__ammom_inp(dummy_this, nd, dtype, dshape, dloc)
    use parameters
    use profiling
    use inputdatatype
    use inputdata, only: inputdata_ammom_inp => ammom_inp
    use, intrinsic :: iso_c_binding, only : c_int
    implicit none
    integer, intent(in) :: dummy_this(2)
    integer(c_int), intent(out) :: nd
    integer(c_int), intent(out) :: dtype
    integer(c_int), dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 3
    dtype = 12
    if (allocated(inputdata_ammom_inp)) then
        dshape(1:3) = shape(inputdata_ammom_inp)
        dloc = loc(inputdata_ammom_inp)
    else
        dloc = 0
    end if
end subroutine f90wrap_inputdata__array__ammom_inp

subroutine f90wrap_inputdata__array__aemom_inp(dummy_this, nd, dtype, dshape, dloc)
    use parameters
    use profiling
    use inputdatatype
    use inputdata, only: inputdata_aemom_inp => aemom_inp
    use, intrinsic :: iso_c_binding, only : c_int
    implicit none
    integer, intent(in) :: dummy_this(2)
    integer(c_int), intent(out) :: nd
    integer(c_int), intent(out) :: dtype
    integer(c_int), dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 4
    dtype = 12
    if (allocated(inputdata_aemom_inp)) then
        dshape(1:4) = shape(inputdata_aemom_inp)
        dloc = loc(inputdata_aemom_inp)
    else
        dloc = 0
    end if
end subroutine f90wrap_inputdata__array__aemom_inp

subroutine f90wrap_inputdata__get__maptype(f90wrap_maptype)
    use inputdata, only: inputdata_maptype => maptype
    implicit none
    integer, intent(out) :: f90wrap_maptype
    
    f90wrap_maptype = inputdata_maptype
end subroutine f90wrap_inputdata__get__maptype

subroutine f90wrap_inputdata__set__maptype(f90wrap_maptype)
    use inputdata, only: inputdata_maptype => maptype
    implicit none
    integer, intent(in) :: f90wrap_maptype
    
    inputdata_maptype = f90wrap_maptype
end subroutine f90wrap_inputdata__set__maptype

subroutine f90wrap_inputdata__get__pre_jfile(f90wrap_pre_jfile)
    use inputdata, only: inputdata_pre_jfile => pre_jfile
    implicit none
    character(30), intent(out) :: f90wrap_pre_jfile
    
    f90wrap_pre_jfile = inputdata_pre_jfile
end subroutine f90wrap_inputdata__get__pre_jfile

subroutine f90wrap_inputdata__set__pre_jfile(f90wrap_pre_jfile)
    use inputdata, only: inputdata_pre_jfile => pre_jfile
    implicit none
    character(30), intent(in) :: f90wrap_pre_jfile
    
    inputdata_pre_jfile = f90wrap_pre_jfile
end subroutine f90wrap_inputdata__set__pre_jfile

subroutine f90wrap_inputdata__array__jfile(dummy_this, nd, dtype, dshape, dloc)
    use parameters
    use profiling
    use inputdatatype
    use inputdata, only: inputdata_jfile => jfile
    use, intrinsic :: iso_c_binding, only : c_int
    implicit none
    integer, intent(in) :: dummy_this(2)
    integer(c_int), intent(out) :: nd
    integer(c_int), intent(out) :: dtype
    integer(c_int), dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 2
    dtype = 2
    if (allocated(inputdata_jfile)) then
        dshape(1:2) = (/len(inputdata_jfile(1)), shape(inputdata_jfile)/)
        dloc = loc(inputdata_jfile)
    else
        dloc = 0
    end if
end subroutine f90wrap_inputdata__array__jfile

subroutine f90wrap_inputdata__array__jfileD(dummy_this, nd, dtype, dshape, dloc)
    use parameters
    use profiling
    use inputdatatype
    use inputdata, only: inputdata_jfiled => jfiled
    use, intrinsic :: iso_c_binding, only : c_int
    implicit none
    integer, intent(in) :: dummy_this(2)
    integer(c_int), intent(out) :: nd
    integer(c_int), intent(out) :: dtype
    integer(c_int), dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 2
    dtype = 2
    if (allocated(inputdata_jfileD)) then
        dshape(1:2) = (/len(inputdata_jfileD(1)), shape(inputdata_jfileD)/)
        dloc = loc(inputdata_jfileD)
    else
        dloc = 0
    end if
end subroutine f90wrap_inputdata__array__jfileD

subroutine f90wrap_inputdata__get__minalgo(f90wrap_minalgo)
    use inputdata, only: inputdata_minalgo => minalgo
    implicit none
    integer, intent(out) :: f90wrap_minalgo
    
    f90wrap_minalgo = inputdata_minalgo
end subroutine f90wrap_inputdata__get__minalgo

subroutine f90wrap_inputdata__set__minalgo(f90wrap_minalgo)
    use inputdata, only: inputdata_minalgo => minalgo
    implicit none
    integer, intent(in) :: f90wrap_minalgo
    
    inputdata_minalgo = f90wrap_minalgo
end subroutine f90wrap_inputdata__set__minalgo

subroutine f90wrap_inputdata__get__minitrmax(f90wrap_minitrmax)
    use inputdata, only: inputdata_minitrmax => minitrmax
    implicit none
    integer, intent(out) :: f90wrap_minitrmax
    
    f90wrap_minitrmax = inputdata_minitrmax
end subroutine f90wrap_inputdata__get__minitrmax

subroutine f90wrap_inputdata__set__minitrmax(f90wrap_minitrmax)
    use inputdata, only: inputdata_minitrmax => minitrmax
    implicit none
    integer, intent(in) :: f90wrap_minitrmax
    
    inputdata_minitrmax = f90wrap_minitrmax
end subroutine f90wrap_inputdata__set__minitrmax

subroutine f90wrap_inputdata__get__mintraj_step(f90wrap_mintraj_step)
    use inputdata, only: inputdata_mintraj_step => mintraj_step
    implicit none
    integer, intent(out) :: f90wrap_mintraj_step
    
    f90wrap_mintraj_step = inputdata_mintraj_step
end subroutine f90wrap_inputdata__get__mintraj_step

subroutine f90wrap_inputdata__set__mintraj_step(f90wrap_mintraj_step)
    use inputdata, only: inputdata_mintraj_step => mintraj_step
    implicit none
    integer, intent(in) :: f90wrap_mintraj_step
    
    inputdata_mintraj_step = f90wrap_mintraj_step
end subroutine f90wrap_inputdata__set__mintraj_step

subroutine f90wrap_inputdata__get__vpodt(f90wrap_vpodt)
    use inputdata, only: inputdata_vpodt => vpodt
    implicit none
    real(8), intent(out) :: f90wrap_vpodt
    
    f90wrap_vpodt = inputdata_vpodt
end subroutine f90wrap_inputdata__get__vpodt

subroutine f90wrap_inputdata__set__vpodt(f90wrap_vpodt)
    use inputdata, only: inputdata_vpodt => vpodt
    implicit none
    real(8), intent(in) :: f90wrap_vpodt
    
    inputdata_vpodt = f90wrap_vpodt
end subroutine f90wrap_inputdata__set__vpodt

subroutine f90wrap_inputdata__get__minftol(f90wrap_minftol)
    use inputdata, only: inputdata_minftol => minftol
    implicit none
    real(8), intent(out) :: f90wrap_minftol
    
    f90wrap_minftol = inputdata_minftol
end subroutine f90wrap_inputdata__get__minftol

subroutine f90wrap_inputdata__set__minftol(f90wrap_minftol)
    use inputdata, only: inputdata_minftol => minftol
    implicit none
    real(8), intent(in) :: f90wrap_minftol
    
    inputdata_minftol = f90wrap_minftol
end subroutine f90wrap_inputdata__set__minftol

subroutine f90wrap_inputdata__get__vpomass(f90wrap_vpomass)
    use inputdata, only: inputdata_vpomass => vpomass
    implicit none
    real(8), intent(out) :: f90wrap_vpomass
    
    f90wrap_vpomass = inputdata_vpomass
end subroutine f90wrap_inputdata__get__vpomass

subroutine f90wrap_inputdata__set__vpomass(f90wrap_vpomass)
    use inputdata, only: inputdata_vpomass => vpomass
    implicit none
    real(8), intent(in) :: f90wrap_vpomass
    
    inputdata_vpomass = f90wrap_vpomass
end subroutine f90wrap_inputdata__set__vpomass

subroutine f90wrap_inputdata__get__initpath(f90wrap_initpath)
    use inputdata, only: inputdata_initpath => initpath
    implicit none
    integer, intent(out) :: f90wrap_initpath
    
    f90wrap_initpath = inputdata_initpath
end subroutine f90wrap_inputdata__get__initpath

subroutine f90wrap_inputdata__set__initpath(f90wrap_initpath)
    use inputdata, only: inputdata_initpath => initpath
    implicit none
    integer, intent(in) :: f90wrap_initpath
    
    inputdata_initpath = f90wrap_initpath
end subroutine f90wrap_inputdata__set__initpath

subroutine f90wrap_inputdata__get__mepitrmax(f90wrap_mepitrmax)
    use inputdata, only: inputdata_mepitrmax => mepitrmax
    implicit none
    integer, intent(out) :: f90wrap_mepitrmax
    
    f90wrap_mepitrmax = inputdata_mepitrmax
end subroutine f90wrap_inputdata__get__mepitrmax

subroutine f90wrap_inputdata__set__mepitrmax(f90wrap_mepitrmax)
    use inputdata, only: inputdata_mepitrmax => mepitrmax
    implicit none
    integer, intent(in) :: f90wrap_mepitrmax
    
    inputdata_mepitrmax = f90wrap_mepitrmax
end subroutine f90wrap_inputdata__set__mepitrmax

subroutine f90wrap_inputdata__get__meptraj_step(f90wrap_meptraj_step)
    use inputdata, only: inputdata_meptraj_step => meptraj_step
    implicit none
    integer, intent(out) :: f90wrap_meptraj_step
    
    f90wrap_meptraj_step = inputdata_meptraj_step
end subroutine f90wrap_inputdata__get__meptraj_step

subroutine f90wrap_inputdata__set__meptraj_step(f90wrap_meptraj_step)
    use inputdata, only: inputdata_meptraj_step => meptraj_step
    implicit none
    integer, intent(in) :: f90wrap_meptraj_step
    
    inputdata_meptraj_step = f90wrap_meptraj_step
end subroutine f90wrap_inputdata__set__meptraj_step

subroutine f90wrap_inputdata__get__spring(f90wrap_spring)
    use inputdata, only: inputdata_spring => spring
    implicit none
    real(8), intent(out) :: f90wrap_spring
    
    f90wrap_spring = inputdata_spring
end subroutine f90wrap_inputdata__get__spring

subroutine f90wrap_inputdata__set__spring(f90wrap_spring)
    use inputdata, only: inputdata_spring => spring
    implicit none
    real(8), intent(in) :: f90wrap_spring
    
    inputdata_spring = f90wrap_spring
end subroutine f90wrap_inputdata__set__spring

subroutine f90wrap_inputdata__get__mepftol(f90wrap_mepftol)
    use inputdata, only: inputdata_mepftol => mepftol
    implicit none
    real(8), intent(out) :: f90wrap_mepftol
    
    f90wrap_mepftol = inputdata_mepftol
end subroutine f90wrap_inputdata__get__mepftol

subroutine f90wrap_inputdata__set__mepftol(f90wrap_mepftol)
    use inputdata, only: inputdata_mepftol => mepftol
    implicit none
    real(8), intent(in) :: f90wrap_mepftol
    
    inputdata_mepftol = f90wrap_mepftol
end subroutine f90wrap_inputdata__set__mepftol

subroutine f90wrap_inputdata__get__mepftol_ci(f90wrap_mepftol_ci)
    use inputdata, only: inputdata_mepftol_ci => mepftol_ci
    implicit none
    real(8), intent(out) :: f90wrap_mepftol_ci
    
    f90wrap_mepftol_ci = inputdata_mepftol_ci
end subroutine f90wrap_inputdata__get__mepftol_ci

subroutine f90wrap_inputdata__set__mepftol_ci(f90wrap_mepftol_ci)
    use inputdata, only: inputdata_mepftol_ci => mepftol_ci
    implicit none
    real(8), intent(in) :: f90wrap_mepftol_ci
    
    inputdata_mepftol_ci = f90wrap_mepftol_ci
end subroutine f90wrap_inputdata__set__mepftol_ci

subroutine f90wrap_inputdata__get__do_gneb(f90wrap_do_gneb)
    use inputdata, only: inputdata_do_gneb => do_gneb
    implicit none
    character(1), intent(out) :: f90wrap_do_gneb
    
    f90wrap_do_gneb = inputdata_do_gneb
end subroutine f90wrap_inputdata__get__do_gneb

subroutine f90wrap_inputdata__set__do_gneb(f90wrap_do_gneb)
    use inputdata, only: inputdata_do_gneb => do_gneb
    implicit none
    character(1), intent(in) :: f90wrap_do_gneb
    
    inputdata_do_gneb = f90wrap_do_gneb
end subroutine f90wrap_inputdata__set__do_gneb

subroutine f90wrap_inputdata__get__do_gneb_ci(f90wrap_do_gneb_ci)
    use inputdata, only: inputdata_do_gneb_ci => do_gneb_ci
    implicit none
    character(1), intent(out) :: f90wrap_do_gneb_ci
    
    f90wrap_do_gneb_ci = inputdata_do_gneb_ci
end subroutine f90wrap_inputdata__get__do_gneb_ci

subroutine f90wrap_inputdata__set__do_gneb_ci(f90wrap_do_gneb_ci)
    use inputdata, only: inputdata_do_gneb_ci => do_gneb_ci
    implicit none
    character(1), intent(in) :: f90wrap_do_gneb_ci
    
    inputdata_do_gneb_ci = f90wrap_do_gneb_ci
end subroutine f90wrap_inputdata__set__do_gneb_ci

subroutine f90wrap_inputdata__get__do_norm_rx(f90wrap_do_norm_rx)
    use inputdata, only: inputdata_do_norm_rx => do_norm_rx
    implicit none
    character(1), intent(out) :: f90wrap_do_norm_rx
    
    f90wrap_do_norm_rx = inputdata_do_norm_rx
end subroutine f90wrap_inputdata__get__do_norm_rx

subroutine f90wrap_inputdata__set__do_norm_rx(f90wrap_do_norm_rx)
    use inputdata, only: inputdata_do_norm_rx => do_norm_rx
    implicit none
    character(1), intent(in) :: f90wrap_do_norm_rx
    
    inputdata_do_norm_rx = f90wrap_do_norm_rx
end subroutine f90wrap_inputdata__set__do_norm_rx

subroutine f90wrap_inputdata__get__en_zero(f90wrap_en_zero)
    use inputdata, only: inputdata_en_zero => en_zero
    implicit none
    character(1), intent(out) :: f90wrap_en_zero
    
    f90wrap_en_zero = inputdata_en_zero
end subroutine f90wrap_inputdata__get__en_zero

subroutine f90wrap_inputdata__set__en_zero(f90wrap_en_zero)
    use inputdata, only: inputdata_en_zero => en_zero
    implicit none
    character(1), intent(in) :: f90wrap_en_zero
    
    inputdata_en_zero = f90wrap_en_zero
end subroutine f90wrap_inputdata__set__en_zero

subroutine f90wrap_inputdata__get__relaxed_if(f90wrap_relaxed_if)
    use inputdata, only: inputdata_relaxed_if => relaxed_if
    implicit none
    character(1), intent(out) :: f90wrap_relaxed_if
    
    f90wrap_relaxed_if = inputdata_relaxed_if
end subroutine f90wrap_inputdata__get__relaxed_if

subroutine f90wrap_inputdata__set__relaxed_if(f90wrap_relaxed_if)
    use inputdata, only: inputdata_relaxed_if => relaxed_if
    implicit none
    character(1), intent(in) :: f90wrap_relaxed_if
    
    inputdata_relaxed_if = f90wrap_relaxed_if
end subroutine f90wrap_inputdata__set__relaxed_if

subroutine f90wrap_inputdata__get__fixed_if(f90wrap_fixed_if)
    use inputdata, only: inputdata_fixed_if => fixed_if
    implicit none
    character(1), intent(out) :: f90wrap_fixed_if
    
    f90wrap_fixed_if = inputdata_fixed_if
end subroutine f90wrap_inputdata__get__fixed_if

subroutine f90wrap_inputdata__set__fixed_if(f90wrap_fixed_if)
    use inputdata, only: inputdata_fixed_if => fixed_if
    implicit none
    character(1), intent(in) :: f90wrap_fixed_if
    
    inputdata_fixed_if = f90wrap_fixed_if
end subroutine f90wrap_inputdata__set__fixed_if

subroutine f90wrap_inputdata__get__prn_gneb_fields(f90wrap_prn_gneb_fields)
    use inputdata, only: inputdata_prn_gneb_fields => prn_gneb_fields
    implicit none
    character(1), intent(out) :: f90wrap_prn_gneb_fields
    
    f90wrap_prn_gneb_fields = inputdata_prn_gneb_fields
end subroutine f90wrap_inputdata__get__prn_gneb_fields

subroutine f90wrap_inputdata__set__prn_gneb_fields(f90wrap_prn_gneb_fields)
    use inputdata, only: inputdata_prn_gneb_fields => prn_gneb_fields
    implicit none
    character(1), intent(in) :: f90wrap_prn_gneb_fields
    
    inputdata_prn_gneb_fields = f90wrap_prn_gneb_fields
end subroutine f90wrap_inputdata__set__prn_gneb_fields

subroutine f90wrap_inputdata__get__do_hess_ini(f90wrap_do_hess_ini)
    use inputdata, only: inputdata_do_hess_ini => do_hess_ini
    implicit none
    character(1), intent(out) :: f90wrap_do_hess_ini
    
    f90wrap_do_hess_ini = inputdata_do_hess_ini
end subroutine f90wrap_inputdata__get__do_hess_ini

subroutine f90wrap_inputdata__set__do_hess_ini(f90wrap_do_hess_ini)
    use inputdata, only: inputdata_do_hess_ini => do_hess_ini
    implicit none
    character(1), intent(in) :: f90wrap_do_hess_ini
    
    inputdata_do_hess_ini = f90wrap_do_hess_ini
end subroutine f90wrap_inputdata__set__do_hess_ini

subroutine f90wrap_inputdata__get__do_hess_fin(f90wrap_do_hess_fin)
    use inputdata, only: inputdata_do_hess_fin => do_hess_fin
    implicit none
    character(1), intent(out) :: f90wrap_do_hess_fin
    
    f90wrap_do_hess_fin = inputdata_do_hess_fin
end subroutine f90wrap_inputdata__get__do_hess_fin

subroutine f90wrap_inputdata__set__do_hess_fin(f90wrap_do_hess_fin)
    use inputdata, only: inputdata_do_hess_fin => do_hess_fin
    implicit none
    character(1), intent(in) :: f90wrap_do_hess_fin
    
    inputdata_do_hess_fin = f90wrap_do_hess_fin
end subroutine f90wrap_inputdata__set__do_hess_fin

subroutine f90wrap_inputdata__get__do_hess_sp(f90wrap_do_hess_sp)
    use inputdata, only: inputdata_do_hess_sp => do_hess_sp
    implicit none
    character(1), intent(out) :: f90wrap_do_hess_sp
    
    f90wrap_do_hess_sp = inputdata_do_hess_sp
end subroutine f90wrap_inputdata__get__do_hess_sp

subroutine f90wrap_inputdata__set__do_hess_sp(f90wrap_do_hess_sp)
    use inputdata, only: inputdata_do_hess_sp => do_hess_sp
    implicit none
    character(1), intent(in) :: f90wrap_do_hess_sp
    
    inputdata_do_hess_sp = f90wrap_do_hess_sp
end subroutine f90wrap_inputdata__set__do_hess_sp

subroutine f90wrap_inputdata__get__eig_0(f90wrap_eig_0)
    use inputdata, only: inputdata_eig_0 => eig_0
    implicit none
    real(8), intent(out) :: f90wrap_eig_0
    
    f90wrap_eig_0 = inputdata_eig_0
end subroutine f90wrap_inputdata__get__eig_0

subroutine f90wrap_inputdata__set__eig_0(f90wrap_eig_0)
    use inputdata, only: inputdata_eig_0 => eig_0
    implicit none
    real(8), intent(in) :: f90wrap_eig_0
    
    inputdata_eig_0 = f90wrap_eig_0
end subroutine f90wrap_inputdata__set__eig_0

subroutine f90wrap_inputdata__get__is_afm(f90wrap_is_afm)
    use inputdata, only: inputdata_is_afm => is_afm
    implicit none
    character(1), intent(out) :: f90wrap_is_afm
    
    f90wrap_is_afm = inputdata_is_afm
end subroutine f90wrap_inputdata__get__is_afm

subroutine f90wrap_inputdata__set__is_afm(f90wrap_is_afm)
    use inputdata, only: inputdata_is_afm => is_afm
    implicit none
    character(1), intent(in) :: f90wrap_is_afm
    
    inputdata_is_afm = f90wrap_is_afm
end subroutine f90wrap_inputdata__set__is_afm

subroutine f90wrap_inputdata__get__sample_num(f90wrap_sample_num)
    use inputdata, only: inputdata_sample_num => sample_num
    implicit none
    integer, intent(out) :: f90wrap_sample_num
    
    f90wrap_sample_num = inputdata_sample_num
end subroutine f90wrap_inputdata__get__sample_num

subroutine f90wrap_inputdata__set__sample_num(f90wrap_sample_num)
    use inputdata, only: inputdata_sample_num => sample_num
    implicit none
    integer, intent(in) :: f90wrap_sample_num
    
    inputdata_sample_num = f90wrap_sample_num
end subroutine f90wrap_inputdata__set__sample_num

subroutine f90wrap_inputdata__get__simid(f90wrap_simid)
    use inputdata, only: inputdata_simid => simid
    implicit none
    character(8), intent(out) :: f90wrap_simid
    
    f90wrap_simid = inputdata_simid
end subroutine f90wrap_inputdata__get__simid

subroutine f90wrap_inputdata__set__simid(f90wrap_simid)
    use inputdata, only: inputdata_simid => simid
    implicit none
    character(8), intent(in) :: f90wrap_simid
    
    inputdata_simid = f90wrap_simid
end subroutine f90wrap_inputdata__set__simid

subroutine f90wrap_inputdata__get__Mensemble(f90wrap_Mensemble)
    use inputdata, only: inputdata_Mensemble => Mensemble
    implicit none
    integer, intent(out) :: f90wrap_Mensemble
    
    f90wrap_Mensemble = inputdata_Mensemble
end subroutine f90wrap_inputdata__get__Mensemble

subroutine f90wrap_inputdata__set__Mensemble(f90wrap_Mensemble)
    use inputdata, only: inputdata_Mensemble => Mensemble
    implicit none
    integer, intent(in) :: f90wrap_Mensemble
    
    inputdata_Mensemble = f90wrap_Mensemble
end subroutine f90wrap_inputdata__set__Mensemble

subroutine f90wrap_inputdata__get__tseed(f90wrap_tseed)
    use inputdata, only: inputdata_tseed => tseed
    implicit none
    integer, intent(out) :: f90wrap_tseed
    
    f90wrap_tseed = inputdata_tseed
end subroutine f90wrap_inputdata__get__tseed

subroutine f90wrap_inputdata__set__tseed(f90wrap_tseed)
    use inputdata, only: inputdata_tseed => tseed
    implicit none
    integer, intent(in) :: f90wrap_tseed
    
    inputdata_tseed = f90wrap_tseed
end subroutine f90wrap_inputdata__set__tseed

subroutine f90wrap_inputdata__get__llg(f90wrap_llg)
    use inputdata, only: inputdata_llg => llg
    implicit none
    integer, intent(out) :: f90wrap_llg
    
    f90wrap_llg = inputdata_llg
end subroutine f90wrap_inputdata__get__llg

subroutine f90wrap_inputdata__set__llg(f90wrap_llg)
    use inputdata, only: inputdata_llg => llg
    implicit none
    integer, intent(in) :: f90wrap_llg
    
    inputdata_llg = f90wrap_llg
end subroutine f90wrap_inputdata__set__llg

subroutine f90wrap_inputdata__get__nstep(f90wrap_nstep)
    use inputdata, only: inputdata_nstep => nstep
    implicit none
    integer, intent(out) :: f90wrap_nstep
    
    f90wrap_nstep = inputdata_nstep
end subroutine f90wrap_inputdata__get__nstep

subroutine f90wrap_inputdata__set__nstep(f90wrap_nstep)
    use inputdata, only: inputdata_nstep => nstep
    implicit none
    integer, intent(in) :: f90wrap_nstep
    
    inputdata_nstep = f90wrap_nstep
end subroutine f90wrap_inputdata__set__nstep

subroutine f90wrap_inputdata__get__SDEalgh(f90wrap_SDEalgh)
    use inputdata, only: inputdata_SDEalgh => SDEalgh
    implicit none
    integer, intent(out) :: f90wrap_SDEalgh
    
    f90wrap_SDEalgh = inputdata_SDEalgh
end subroutine f90wrap_inputdata__get__SDEalgh

subroutine f90wrap_inputdata__set__SDEalgh(f90wrap_SDEalgh)
    use inputdata, only: inputdata_SDEalgh => SDEalgh
    implicit none
    integer, intent(in) :: f90wrap_SDEalgh
    
    inputdata_SDEalgh = f90wrap_SDEalgh
end subroutine f90wrap_inputdata__set__SDEalgh

subroutine f90wrap_inputdata__get__ipSDEalgh(f90wrap_ipSDEalgh)
    use inputdata, only: inputdata_ipSDEalgh => ipSDEalgh
    implicit none
    integer, intent(out) :: f90wrap_ipSDEalgh
    
    f90wrap_ipSDEalgh = inputdata_ipSDEalgh
end subroutine f90wrap_inputdata__get__ipSDEalgh

subroutine f90wrap_inputdata__set__ipSDEalgh(f90wrap_ipSDEalgh)
    use inputdata, only: inputdata_ipSDEalgh => ipSDEalgh
    implicit none
    integer, intent(in) :: f90wrap_ipSDEalgh
    
    inputdata_ipSDEalgh = f90wrap_ipSDEalgh
end subroutine f90wrap_inputdata__set__ipSDEalgh

subroutine f90wrap_inputdata__get__aunits(f90wrap_aunits)
    use inputdata, only: inputdata_aunits => aunits
    implicit none
    character(1), intent(out) :: f90wrap_aunits
    
    f90wrap_aunits = inputdata_aunits
end subroutine f90wrap_inputdata__get__aunits

subroutine f90wrap_inputdata__set__aunits(f90wrap_aunits)
    use inputdata, only: inputdata_aunits => aunits
    implicit none
    character(1), intent(in) :: f90wrap_aunits
    
    inputdata_aunits = f90wrap_aunits
end subroutine f90wrap_inputdata__set__aunits

subroutine f90wrap_inputdata__get__perp(f90wrap_perp)
    use inputdata, only: inputdata_perp => perp
    implicit none
    character(1), intent(out) :: f90wrap_perp
    
    f90wrap_perp = inputdata_perp
end subroutine f90wrap_inputdata__get__perp

subroutine f90wrap_inputdata__set__perp(f90wrap_perp)
    use inputdata, only: inputdata_perp => perp
    implicit none
    character(1), intent(in) :: f90wrap_perp
    
    inputdata_perp = f90wrap_perp
end subroutine f90wrap_inputdata__set__perp

subroutine f90wrap_inputdata__get__mompar(f90wrap_mompar)
    use inputdata, only: inputdata_mompar => mompar
    implicit none
    integer, intent(out) :: f90wrap_mompar
    
    f90wrap_mompar = inputdata_mompar
end subroutine f90wrap_inputdata__get__mompar

subroutine f90wrap_inputdata__set__mompar(f90wrap_mompar)
    use inputdata, only: inputdata_mompar => mompar
    implicit none
    integer, intent(in) :: f90wrap_mompar
    
    inputdata_mompar = f90wrap_mompar
end subroutine f90wrap_inputdata__set__mompar

subroutine f90wrap_inputdata__get__heisout(f90wrap_heisout)
    use inputdata, only: inputdata_heisout => heisout
    implicit none
    integer, intent(out) :: f90wrap_heisout
    
    f90wrap_heisout = inputdata_heisout
end subroutine f90wrap_inputdata__get__heisout

subroutine f90wrap_inputdata__set__heisout(f90wrap_heisout)
    use inputdata, only: inputdata_heisout => heisout
    implicit none
    integer, intent(in) :: f90wrap_heisout
    
    inputdata_heisout = f90wrap_heisout
end subroutine f90wrap_inputdata__set__heisout

subroutine f90wrap_inputdata__get__evolveout(f90wrap_evolveout)
    use inputdata, only: inputdata_evolveout => evolveout
    implicit none
    integer, intent(out) :: f90wrap_evolveout
    
    f90wrap_evolveout = inputdata_evolveout
end subroutine f90wrap_inputdata__get__evolveout

subroutine f90wrap_inputdata__set__evolveout(f90wrap_evolveout)
    use inputdata, only: inputdata_evolveout => evolveout
    implicit none
    integer, intent(in) :: f90wrap_evolveout
    
    inputdata_evolveout = f90wrap_evolveout
end subroutine f90wrap_inputdata__set__evolveout

subroutine f90wrap_inputdata__get__plotenergy(f90wrap_plotenergy)
    use inputdata, only: inputdata_plotenergy => plotenergy
    implicit none
    integer, intent(out) :: f90wrap_plotenergy
    
    f90wrap_plotenergy = inputdata_plotenergy
end subroutine f90wrap_inputdata__get__plotenergy

subroutine f90wrap_inputdata__set__plotenergy(f90wrap_plotenergy)
    use inputdata, only: inputdata_plotenergy => plotenergy
    implicit none
    integer, intent(in) :: f90wrap_plotenergy
    
    inputdata_plotenergy = f90wrap_plotenergy
end subroutine f90wrap_inputdata__set__plotenergy

subroutine f90wrap_inputdata__get__do_hoc_debug(f90wrap_do_hoc_debug)
    use inputdata, only: inputdata_do_hoc_debug => do_hoc_debug
    implicit none
    integer, intent(out) :: f90wrap_do_hoc_debug
    
    f90wrap_do_hoc_debug = inputdata_do_hoc_debug
end subroutine f90wrap_inputdata__get__do_hoc_debug

subroutine f90wrap_inputdata__set__do_hoc_debug(f90wrap_do_hoc_debug)
    use inputdata, only: inputdata_do_hoc_debug => do_hoc_debug
    implicit none
    integer, intent(in) :: f90wrap_do_hoc_debug
    
    inputdata_do_hoc_debug = f90wrap_do_hoc_debug
end subroutine f90wrap_inputdata__set__do_hoc_debug

subroutine f90wrap_inputdata__get__do_prnstruct(f90wrap_do_prnstruct)
    use inputdata, only: inputdata_do_prnstruct => do_prnstruct
    implicit none
    integer, intent(out) :: f90wrap_do_prnstruct
    
    f90wrap_do_prnstruct = inputdata_do_prnstruct
end subroutine f90wrap_inputdata__get__do_prnstruct

subroutine f90wrap_inputdata__set__do_prnstruct(f90wrap_do_prnstruct)
    use inputdata, only: inputdata_do_prnstruct => do_prnstruct
    implicit none
    integer, intent(in) :: f90wrap_do_prnstruct
    
    inputdata_do_prnstruct = f90wrap_do_prnstruct
end subroutine f90wrap_inputdata__set__do_prnstruct

subroutine f90wrap_inputdata__get__do_storeham(f90wrap_do_storeham)
    use inputdata, only: inputdata_do_storeham => do_storeham
    implicit none
    integer, intent(out) :: f90wrap_do_storeham
    
    f90wrap_do_storeham = inputdata_do_storeham
end subroutine f90wrap_inputdata__get__do_storeham

subroutine f90wrap_inputdata__set__do_storeham(f90wrap_do_storeham)
    use inputdata, only: inputdata_do_storeham => do_storeham
    implicit none
    integer, intent(in) :: f90wrap_do_storeham
    
    inputdata_do_storeham = f90wrap_do_storeham
end subroutine f90wrap_inputdata__set__do_storeham

subroutine f90wrap_inputdata__get__do_prn_poscar(f90wrap_do_prn_poscar)
    use inputdata, only: inputdata_do_prn_poscar => do_prn_poscar
    implicit none
    integer, intent(out) :: f90wrap_do_prn_poscar
    
    f90wrap_do_prn_poscar = inputdata_do_prn_poscar
end subroutine f90wrap_inputdata__get__do_prn_poscar

subroutine f90wrap_inputdata__set__do_prn_poscar(f90wrap_do_prn_poscar)
    use inputdata, only: inputdata_do_prn_poscar => do_prn_poscar
    implicit none
    integer, intent(in) :: f90wrap_do_prn_poscar
    
    inputdata_do_prn_poscar = f90wrap_do_prn_poscar
end subroutine f90wrap_inputdata__set__do_prn_poscar

subroutine f90wrap_inputdata__get__do_prn_elk(f90wrap_do_prn_elk)
    use inputdata, only: inputdata_do_prn_elk => do_prn_elk
    implicit none
    integer, intent(out) :: f90wrap_do_prn_elk
    
    f90wrap_do_prn_elk = inputdata_do_prn_elk
end subroutine f90wrap_inputdata__get__do_prn_elk

subroutine f90wrap_inputdata__set__do_prn_elk(f90wrap_do_prn_elk)
    use inputdata, only: inputdata_do_prn_elk => do_prn_elk
    implicit none
    integer, intent(in) :: f90wrap_do_prn_elk
    
    inputdata_do_prn_elk = f90wrap_do_prn_elk
end subroutine f90wrap_inputdata__set__do_prn_elk

subroutine f90wrap_inputdata__get__do_read_elk(f90wrap_do_read_elk)
    use inputdata, only: inputdata_do_read_elk => do_read_elk
    implicit none
    integer, intent(out) :: f90wrap_do_read_elk
    
    f90wrap_do_read_elk = inputdata_do_read_elk
end subroutine f90wrap_inputdata__get__do_read_elk

subroutine f90wrap_inputdata__set__do_read_elk(f90wrap_do_read_elk)
    use inputdata, only: inputdata_do_read_elk => do_read_elk
    implicit none
    integer, intent(in) :: f90wrap_do_read_elk
    
    inputdata_do_read_elk = f90wrap_do_read_elk
end subroutine f90wrap_inputdata__set__do_read_elk

subroutine f90wrap_inputdata__get__compensate_drift(f90wrap_compensate_drift)
    use inputdata, only: inputdata_compensate_drift => compensate_drift
    implicit none
    integer, intent(out) :: f90wrap_compensate_drift
    
    f90wrap_compensate_drift = inputdata_compensate_drift
end subroutine f90wrap_inputdata__get__compensate_drift

subroutine f90wrap_inputdata__set__compensate_drift(f90wrap_compensate_drift)
    use inputdata, only: inputdata_compensate_drift => compensate_drift
    implicit none
    integer, intent(in) :: f90wrap_compensate_drift
    
    inputdata_compensate_drift = f90wrap_compensate_drift
end subroutine f90wrap_inputdata__set__compensate_drift

subroutine f90wrap_inputdata__get__do_sparse(f90wrap_do_sparse)
    use inputdata, only: inputdata_do_sparse => do_sparse
    implicit none
    character(1), intent(out) :: f90wrap_do_sparse
    
    f90wrap_do_sparse = inputdata_do_sparse
end subroutine f90wrap_inputdata__get__do_sparse

subroutine f90wrap_inputdata__set__do_sparse(f90wrap_do_sparse)
    use inputdata, only: inputdata_do_sparse => do_sparse
    implicit none
    character(1), intent(in) :: f90wrap_do_sparse
    
    inputdata_do_sparse = f90wrap_do_sparse
end subroutine f90wrap_inputdata__set__do_sparse

subroutine f90wrap_inputdata__get__do_reduced(f90wrap_do_reduced)
    use inputdata, only: inputdata_do_reduced => do_reduced
    implicit none
    character(1), intent(out) :: f90wrap_do_reduced
    
    f90wrap_do_reduced = inputdata_do_reduced
end subroutine f90wrap_inputdata__get__do_reduced

subroutine f90wrap_inputdata__set__do_reduced(f90wrap_do_reduced)
    use inputdata, only: inputdata_do_reduced => do_reduced
    implicit none
    character(1), intent(in) :: f90wrap_do_reduced
    
    inputdata_do_reduced = f90wrap_do_reduced
end subroutine f90wrap_inputdata__set__do_reduced

subroutine f90wrap_inputdata__get__Temp(f90wrap_Temp)
    use inputdata, only: inputdata_Temp => Temp
    implicit none
    real(8), intent(out) :: f90wrap_Temp
    
    f90wrap_Temp = inputdata_Temp
end subroutine f90wrap_inputdata__get__Temp

subroutine f90wrap_inputdata__set__Temp(f90wrap_Temp)
    use inputdata, only: inputdata_Temp => Temp
    implicit none
    real(8), intent(in) :: f90wrap_Temp
    
    inputdata_Temp = f90wrap_Temp
end subroutine f90wrap_inputdata__set__Temp

subroutine f90wrap_inputdata__get__delta_t(f90wrap_delta_t)
    use inputdata, only: inputdata_delta_t => delta_t
    implicit none
    real(8), intent(out) :: f90wrap_delta_t
    
    f90wrap_delta_t = inputdata_delta_t
end subroutine f90wrap_inputdata__get__delta_t

subroutine f90wrap_inputdata__set__delta_t(f90wrap_delta_t)
    use inputdata, only: inputdata_delta_t => delta_t
    implicit none
    real(8), intent(in) :: f90wrap_delta_t
    
    inputdata_delta_t = f90wrap_delta_t
end subroutine f90wrap_inputdata__set__delta_t

subroutine f90wrap_inputdata__get__relaxtime(f90wrap_relaxtime)
    use inputdata, only: inputdata_relaxtime => relaxtime
    implicit none
    real(8), intent(out) :: f90wrap_relaxtime
    
    f90wrap_relaxtime = inputdata_relaxtime
end subroutine f90wrap_inputdata__get__relaxtime

subroutine f90wrap_inputdata__set__relaxtime(f90wrap_relaxtime)
    use inputdata, only: inputdata_relaxtime => relaxtime
    implicit none
    real(8), intent(in) :: f90wrap_relaxtime
    
    inputdata_relaxtime = f90wrap_relaxtime
end subroutine f90wrap_inputdata__set__relaxtime

subroutine f90wrap_inputdata__get__mplambda1(f90wrap_mplambda1)
    use inputdata, only: inputdata_mplambda1 => mplambda1
    implicit none
    real(8), intent(out) :: f90wrap_mplambda1
    
    f90wrap_mplambda1 = inputdata_mplambda1
end subroutine f90wrap_inputdata__get__mplambda1

subroutine f90wrap_inputdata__set__mplambda1(f90wrap_mplambda1)
    use inputdata, only: inputdata_mplambda1 => mplambda1
    implicit none
    real(8), intent(in) :: f90wrap_mplambda1
    
    inputdata_mplambda1 = f90wrap_mplambda1
end subroutine f90wrap_inputdata__set__mplambda1

subroutine f90wrap_inputdata__get__mplambda2(f90wrap_mplambda2)
    use inputdata, only: inputdata_mplambda2 => mplambda2
    implicit none
    real(8), intent(out) :: f90wrap_mplambda2
    
    f90wrap_mplambda2 = inputdata_mplambda2
end subroutine f90wrap_inputdata__get__mplambda2

subroutine f90wrap_inputdata__set__mplambda2(f90wrap_mplambda2)
    use inputdata, only: inputdata_mplambda2 => mplambda2
    implicit none
    real(8), intent(in) :: f90wrap_mplambda2
    
    inputdata_mplambda2 = f90wrap_mplambda2
end subroutine f90wrap_inputdata__set__mplambda2

subroutine f90wrap_inputdata__get__mode(f90wrap_mode)
    use inputdata, only: inputdata_mode => mode
    implicit none
    character(2), intent(out) :: f90wrap_mode
    
    f90wrap_mode = inputdata_mode
end subroutine f90wrap_inputdata__get__mode

subroutine f90wrap_inputdata__set__mode(f90wrap_mode)
    use inputdata, only: inputdata_mode => mode
    implicit none
    character(2), intent(in) :: f90wrap_mode
    
    inputdata_mode = f90wrap_mode
end subroutine f90wrap_inputdata__set__mode

subroutine f90wrap_inputdata__array__hfield(dummy_this, nd, dtype, dshape, dloc)
    use parameters
    use profiling
    use inputdatatype
    use inputdata, only: inputdata_hfield => hfield
    use, intrinsic :: iso_c_binding, only : c_int
    implicit none
    integer, intent(in) :: dummy_this(2)
    integer(c_int), intent(out) :: nd
    integer(c_int), intent(out) :: dtype
    integer(c_int), dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 1
    dtype = 12
    dshape(1:1) = shape(inputdata_hfield)
    dloc = loc(inputdata_hfield)
end subroutine f90wrap_inputdata__array__hfield

subroutine f90wrap_inputdata__get__mpNlines(f90wrap_mpNlines)
    use inputdata, only: inputdata_mpNlines => mpNlines
    implicit none
    integer, intent(out) :: f90wrap_mpNlines
    
    f90wrap_mpNlines = inputdata_mpNlines
end subroutine f90wrap_inputdata__get__mpNlines

subroutine f90wrap_inputdata__set__mpNlines(f90wrap_mpNlines)
    use inputdata, only: inputdata_mpNlines => mpNlines
    implicit none
    integer, intent(in) :: f90wrap_mpNlines
    
    inputdata_mpNlines = f90wrap_mpNlines
end subroutine f90wrap_inputdata__set__mpNlines

subroutine f90wrap_inputdata__array__mpdamping1(dummy_this, nd, dtype, dshape, dloc)
    use parameters
    use profiling
    use inputdatatype
    use inputdata, only: inputdata_mpdamping1 => mpdamping1
    use, intrinsic :: iso_c_binding, only : c_int
    implicit none
    integer, intent(in) :: dummy_this(2)
    integer(c_int), intent(out) :: nd
    integer(c_int), intent(out) :: dtype
    integer(c_int), dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 1
    dtype = 12
    if (allocated(inputdata_mpdamping1)) then
        dshape(1:1) = shape(inputdata_mpdamping1)
        dloc = loc(inputdata_mpdamping1)
    else
        dloc = 0
    end if
end subroutine f90wrap_inputdata__array__mpdamping1

subroutine f90wrap_inputdata__array__mpdamping2(dummy_this, nd, dtype, dshape, dloc)
    use parameters
    use profiling
    use inputdatatype
    use inputdata, only: inputdata_mpdamping2 => mpdamping2
    use, intrinsic :: iso_c_binding, only : c_int
    implicit none
    integer, intent(in) :: dummy_this(2)
    integer(c_int), intent(out) :: nd
    integer(c_int), intent(out) :: dtype
    integer(c_int), dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 1
    dtype = 12
    if (allocated(inputdata_mpdamping2)) then
        dshape(1:1) = shape(inputdata_mpdamping2)
        dloc = loc(inputdata_mpdamping2)
    else
        dloc = 0
    end if
end subroutine f90wrap_inputdata__array__mpdamping2

subroutine f90wrap_inputdata__array__mpdampingalloy1(dummy_this, nd, dtype, dshape, dloc)
    use parameters
    use profiling
    use inputdatatype
    use inputdata, only: inputdata_mpdampingalloy1 => mpdampingalloy1
    use, intrinsic :: iso_c_binding, only : c_int
    implicit none
    integer, intent(in) :: dummy_this(2)
    integer(c_int), intent(out) :: nd
    integer(c_int), intent(out) :: dtype
    integer(c_int), dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 2
    dtype = 12
    if (allocated(inputdata_mpdampingalloy1)) then
        dshape(1:2) = shape(inputdata_mpdampingalloy1)
        dloc = loc(inputdata_mpdampingalloy1)
    else
        dloc = 0
    end if
end subroutine f90wrap_inputdata__array__mpdampingalloy1

subroutine f90wrap_inputdata__array__mpdampingalloy2(dummy_this, nd, dtype, dshape, dloc)
    use parameters
    use profiling
    use inputdatatype
    use inputdata, only: inputdata_mpdampingalloy2 => mpdampingalloy2
    use, intrinsic :: iso_c_binding, only : c_int
    implicit none
    integer, intent(in) :: dummy_this(2)
    integer(c_int), intent(out) :: nd
    integer(c_int), intent(out) :: dtype
    integer(c_int), dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 2
    dtype = 12
    if (allocated(inputdata_mpdampingalloy2)) then
        dshape(1:2) = shape(inputdata_mpdampingalloy2)
        dloc = loc(inputdata_mpdampingalloy2)
    else
        dloc = 0
    end if
end subroutine f90wrap_inputdata__array__mpdampingalloy2

subroutine f90wrap_inputdata__get__do_site_damping(f90wrap_do_site_damping)
    use inputdata, only: inputdata_do_site_damping => do_site_damping
    implicit none
    character(1), intent(out) :: f90wrap_do_site_damping
    
    f90wrap_do_site_damping = inputdata_do_site_damping
end subroutine f90wrap_inputdata__get__do_site_damping

subroutine f90wrap_inputdata__set__do_site_damping(f90wrap_do_site_damping)
    use inputdata, only: inputdata_do_site_damping => do_site_damping
    implicit none
    character(1), intent(in) :: f90wrap_do_site_damping
    
    inputdata_do_site_damping = f90wrap_do_site_damping
end subroutine f90wrap_inputdata__set__do_site_damping

subroutine f90wrap_inputdata__get__mp_dampfile(f90wrap_mp_dampfile)
    use inputdata, only: inputdata_mp_dampfile => mp_dampfile
    implicit none
    character(35), intent(out) :: f90wrap_mp_dampfile
    
    f90wrap_mp_dampfile = inputdata_mp_dampfile
end subroutine f90wrap_inputdata__get__mp_dampfile

subroutine f90wrap_inputdata__set__mp_dampfile(f90wrap_mp_dampfile)
    use inputdata, only: inputdata_mp_dampfile => mp_dampfile
    implicit none
    character(35), intent(in) :: f90wrap_mp_dampfile
    
    inputdata_mp_dampfile = f90wrap_mp_dampfile
end subroutine f90wrap_inputdata__set__mp_dampfile

subroutine f90wrap_inputdata__get__do_mc(f90wrap_do_mc)
    use inputdata, only: inputdata_do_mc => do_mc
    implicit none
    integer, intent(out) :: f90wrap_do_mc
    
    f90wrap_do_mc = inputdata_do_mc
end subroutine f90wrap_inputdata__get__do_mc

subroutine f90wrap_inputdata__set__do_mc(f90wrap_do_mc)
    use inputdata, only: inputdata_do_mc => do_mc
    implicit none
    integer, intent(in) :: f90wrap_do_mc
    
    inputdata_do_mc = f90wrap_do_mc
end subroutine f90wrap_inputdata__set__do_mc

subroutine f90wrap_inputdata__get__mcnstep(f90wrap_mcnstep)
    use inputdata, only: inputdata_mcnstep => mcnstep
    implicit none
    integer, intent(out) :: f90wrap_mcnstep
    
    f90wrap_mcnstep = inputdata_mcnstep
end subroutine f90wrap_inputdata__get__mcnstep

subroutine f90wrap_inputdata__set__mcnstep(f90wrap_mcnstep)
    use inputdata, only: inputdata_mcnstep => mcnstep
    implicit none
    integer, intent(in) :: f90wrap_mcnstep
    
    inputdata_mcnstep = f90wrap_mcnstep
end subroutine f90wrap_inputdata__set__mcnstep

subroutine f90wrap_inputdata__get__mcavrg_step(f90wrap_mcavrg_step)
    use inputdata, only: inputdata_mcavrg_step => mcavrg_step
    implicit none
    integer, intent(out) :: f90wrap_mcavrg_step
    
    f90wrap_mcavrg_step = inputdata_mcavrg_step
end subroutine f90wrap_inputdata__get__mcavrg_step

subroutine f90wrap_inputdata__set__mcavrg_step(f90wrap_mcavrg_step)
    use inputdata, only: inputdata_mcavrg_step => mcavrg_step
    implicit none
    integer, intent(in) :: f90wrap_mcavrg_step
    
    inputdata_mcavrg_step = f90wrap_mcavrg_step
end subroutine f90wrap_inputdata__set__mcavrg_step

subroutine f90wrap_inputdata__get__mcavrg_buff(f90wrap_mcavrg_buff)
    use inputdata, only: inputdata_mcavrg_buff => mcavrg_buff
    implicit none
    integer, intent(out) :: f90wrap_mcavrg_buff
    
    f90wrap_mcavrg_buff = inputdata_mcavrg_buff
end subroutine f90wrap_inputdata__get__mcavrg_buff

subroutine f90wrap_inputdata__set__mcavrg_buff(f90wrap_mcavrg_buff)
    use inputdata, only: inputdata_mcavrg_buff => mcavrg_buff
    implicit none
    integer, intent(in) :: f90wrap_mcavrg_buff
    
    inputdata_mcavrg_buff = f90wrap_mcavrg_buff
end subroutine f90wrap_inputdata__set__mcavrg_buff

subroutine f90wrap_inputdata__get__ipnphase(f90wrap_ipnphase)
    use inputdata, only: inputdata_ipnphase => ipnphase
    implicit none
    integer, intent(out) :: f90wrap_ipnphase
    
    f90wrap_ipnphase = inputdata_ipnphase
end subroutine f90wrap_inputdata__get__ipnphase

subroutine f90wrap_inputdata__set__ipnphase(f90wrap_ipnphase)
    use inputdata, only: inputdata_ipnphase => ipnphase
    implicit none
    integer, intent(in) :: f90wrap_ipnphase
    
    inputdata_ipnphase = f90wrap_ipnphase
end subroutine f90wrap_inputdata__set__ipnphase

subroutine f90wrap_inputdata__get__ipmcnphase(f90wrap_ipmcnphase)
    use inputdata, only: inputdata_ipmcnphase => ipmcnphase
    implicit none
    integer, intent(out) :: f90wrap_ipmcnphase
    
    f90wrap_ipmcnphase = inputdata_ipmcnphase
end subroutine f90wrap_inputdata__get__ipmcnphase

subroutine f90wrap_inputdata__set__ipmcnphase(f90wrap_ipmcnphase)
    use inputdata, only: inputdata_ipmcnphase => ipmcnphase
    implicit none
    integer, intent(in) :: f90wrap_ipmcnphase
    
    inputdata_ipmcnphase = f90wrap_ipmcnphase
end subroutine f90wrap_inputdata__set__ipmcnphase

subroutine f90wrap_inputdata__get__ipmode(f90wrap_ipmode)
    use inputdata, only: inputdata_ipmode => ipmode
    implicit none
    character(2), intent(out) :: f90wrap_ipmode
    
    f90wrap_ipmode = inputdata_ipmode
end subroutine f90wrap_inputdata__get__ipmode

subroutine f90wrap_inputdata__set__ipmode(f90wrap_ipmode)
    use inputdata, only: inputdata_ipmode => ipmode
    implicit none
    character(2), intent(in) :: f90wrap_ipmode
    
    inputdata_ipmode = f90wrap_ipmode
end subroutine f90wrap_inputdata__set__ipmode

subroutine f90wrap_inputdata__array__ipnstep(dummy_this, nd, dtype, dshape, dloc)
    use parameters
    use profiling
    use inputdatatype
    use inputdata, only: inputdata_ipnstep => ipnstep
    use, intrinsic :: iso_c_binding, only : c_int
    implicit none
    integer, intent(in) :: dummy_this(2)
    integer(c_int), intent(out) :: nd
    integer(c_int), intent(out) :: dtype
    integer(c_int), dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 1
    dtype = 5
    if (allocated(inputdata_ipnstep)) then
        dshape(1:1) = shape(inputdata_ipnstep)
        dloc = loc(inputdata_ipnstep)
    else
        dloc = 0
    end if
end subroutine f90wrap_inputdata__array__ipnstep

subroutine f90wrap_inputdata__array__ipmcnstep(dummy_this, nd, dtype, dshape, dloc)
    use parameters
    use profiling
    use inputdatatype
    use inputdata, only: inputdata_ipmcnstep => ipmcnstep
    use, intrinsic :: iso_c_binding, only : c_int
    implicit none
    integer, intent(in) :: dummy_this(2)
    integer(c_int), intent(out) :: nd
    integer(c_int), intent(out) :: dtype
    integer(c_int), dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 1
    dtype = 5
    if (allocated(inputdata_ipmcnstep)) then
        dshape(1:1) = shape(inputdata_ipmcnstep)
        dloc = loc(inputdata_ipmcnstep)
    else
        dloc = 0
    end if
end subroutine f90wrap_inputdata__array__ipmcnstep

subroutine f90wrap_inputdata__array__iphfield(dummy_this, nd, dtype, dshape, dloc)
    use parameters
    use profiling
    use inputdatatype
    use inputdata, only: inputdata_iphfield => iphfield
    use, intrinsic :: iso_c_binding, only : c_int
    implicit none
    integer, intent(in) :: dummy_this(2)
    integer(c_int), intent(out) :: nd
    integer(c_int), intent(out) :: dtype
    integer(c_int), dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 1
    dtype = 12
    dshape(1:1) = shape(inputdata_iphfield)
    dloc = loc(inputdata_iphfield)
end subroutine f90wrap_inputdata__array__iphfield

subroutine f90wrap_inputdata__array__ipdelta_t(dummy_this, nd, dtype, dshape, dloc)
    use parameters
    use profiling
    use inputdatatype
    use inputdata, only: inputdata_ipdelta_t => ipdelta_t
    use, intrinsic :: iso_c_binding, only : c_int
    implicit none
    integer, intent(in) :: dummy_this(2)
    integer(c_int), intent(out) :: nd
    integer(c_int), intent(out) :: dtype
    integer(c_int), dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 1
    dtype = 12
    if (allocated(inputdata_ipdelta_t)) then
        dshape(1:1) = shape(inputdata_ipdelta_t)
        dloc = loc(inputdata_ipdelta_t)
    else
        dloc = 0
    end if
end subroutine f90wrap_inputdata__array__ipdelta_t

subroutine f90wrap_inputdata__array__iplambda1(dummy_this, nd, dtype, dshape, dloc)
    use parameters
    use profiling
    use inputdatatype
    use inputdata, only: inputdata_iplambda1 => iplambda1
    use, intrinsic :: iso_c_binding, only : c_int
    implicit none
    integer, intent(in) :: dummy_this(2)
    integer(c_int), intent(out) :: nd
    integer(c_int), intent(out) :: dtype
    integer(c_int), dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 1
    dtype = 12
    if (allocated(inputdata_iplambda1)) then
        dshape(1:1) = shape(inputdata_iplambda1)
        dloc = loc(inputdata_iplambda1)
    else
        dloc = 0
    end if
end subroutine f90wrap_inputdata__array__iplambda1

subroutine f90wrap_inputdata__array__iplambda2(dummy_this, nd, dtype, dshape, dloc)
    use parameters
    use profiling
    use inputdatatype
    use inputdata, only: inputdata_iplambda2 => iplambda2
    use, intrinsic :: iso_c_binding, only : c_int
    implicit none
    integer, intent(in) :: dummy_this(2)
    integer(c_int), intent(out) :: nd
    integer(c_int), intent(out) :: dtype
    integer(c_int), dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 1
    dtype = 12
    if (allocated(inputdata_iplambda2)) then
        dshape(1:1) = shape(inputdata_iplambda2)
        dloc = loc(inputdata_iplambda2)
    else
        dloc = 0
    end if
end subroutine f90wrap_inputdata__array__iplambda2

subroutine f90wrap_inputdata__array__ipTemp(dummy_this, nd, dtype, dshape, dloc)
    use parameters
    use profiling
    use inputdatatype
    use inputdata, only: inputdata_iptemp => iptemp
    use, intrinsic :: iso_c_binding, only : c_int
    implicit none
    integer, intent(in) :: dummy_this(2)
    integer(c_int), intent(out) :: nd
    integer(c_int), intent(out) :: dtype
    integer(c_int), dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 1
    dtype = 12
    if (allocated(inputdata_ipTemp)) then
        dshape(1:1) = shape(inputdata_ipTemp)
        dloc = loc(inputdata_ipTemp)
    else
        dloc = 0
    end if
end subroutine f90wrap_inputdata__array__ipTemp

subroutine f90wrap_inputdata__array__ipdamping1(dummy_this, nd, dtype, dshape, dloc)
    use parameters
    use profiling
    use inputdatatype
    use inputdata, only: inputdata_ipdamping1 => ipdamping1
    use, intrinsic :: iso_c_binding, only : c_int
    implicit none
    integer, intent(in) :: dummy_this(2)
    integer(c_int), intent(out) :: nd
    integer(c_int), intent(out) :: dtype
    integer(c_int), dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 2
    dtype = 12
    if (allocated(inputdata_ipdamping1)) then
        dshape(1:2) = shape(inputdata_ipdamping1)
        dloc = loc(inputdata_ipdamping1)
    else
        dloc = 0
    end if
end subroutine f90wrap_inputdata__array__ipdamping1

subroutine f90wrap_inputdata__array__ipdamping2(dummy_this, nd, dtype, dshape, dloc)
    use parameters
    use profiling
    use inputdatatype
    use inputdata, only: inputdata_ipdamping2 => ipdamping2
    use, intrinsic :: iso_c_binding, only : c_int
    implicit none
    integer, intent(in) :: dummy_this(2)
    integer(c_int), intent(out) :: nd
    integer(c_int), intent(out) :: dtype
    integer(c_int), dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 2
    dtype = 12
    if (allocated(inputdata_ipdamping2)) then
        dshape(1:2) = shape(inputdata_ipdamping2)
        dloc = loc(inputdata_ipdamping2)
    else
        dloc = 0
    end if
end subroutine f90wrap_inputdata__array__ipdamping2

subroutine f90wrap_inputdata__array__ipdampingalloy1(dummy_this, nd, dtype, dshape, dloc)
    use parameters
    use profiling
    use inputdatatype
    use inputdata, only: inputdata_ipdampingalloy1 => ipdampingalloy1
    use, intrinsic :: iso_c_binding, only : c_int
    implicit none
    integer, intent(in) :: dummy_this(2)
    integer(c_int), intent(out) :: nd
    integer(c_int), intent(out) :: dtype
    integer(c_int), dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 3
    dtype = 12
    if (allocated(inputdata_ipdampingalloy1)) then
        dshape(1:3) = shape(inputdata_ipdampingalloy1)
        dloc = loc(inputdata_ipdampingalloy1)
    else
        dloc = 0
    end if
end subroutine f90wrap_inputdata__array__ipdampingalloy1

subroutine f90wrap_inputdata__array__ipdampingalloy2(dummy_this, nd, dtype, dshape, dloc)
    use parameters
    use profiling
    use inputdatatype
    use inputdata, only: inputdata_ipdampingalloy2 => ipdampingalloy2
    use, intrinsic :: iso_c_binding, only : c_int
    implicit none
    integer, intent(in) :: dummy_this(2)
    integer(c_int), intent(out) :: nd
    integer(c_int), intent(out) :: dtype
    integer(c_int), dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 3
    dtype = 12
    if (allocated(inputdata_ipdampingalloy2)) then
        dshape(1:3) = shape(inputdata_ipdampingalloy2)
        dloc = loc(inputdata_ipdampingalloy2)
    else
        dloc = 0
    end if
end subroutine f90wrap_inputdata__array__ipdampingalloy2

subroutine f90wrap_inputdata__get__do_site_ip_damping(f90wrap_do_site_ip_damping)
    use inputdata, only: inputdata_do_site_ip_damping => do_site_ip_damping
    implicit none
    character(1), intent(out) :: f90wrap_do_site_ip_damping
    
    f90wrap_do_site_ip_damping = inputdata_do_site_ip_damping
end subroutine f90wrap_inputdata__get__do_site_ip_damping

subroutine f90wrap_inputdata__set__do_site_ip_damping(f90wrap_do_site_ip_damping)
    use inputdata, only: inputdata_do_site_ip_damping => do_site_ip_damping
    implicit none
    character(1), intent(in) :: f90wrap_do_site_ip_damping
    
    inputdata_do_site_ip_damping = f90wrap_do_site_ip_damping
end subroutine f90wrap_inputdata__set__do_site_ip_damping

subroutine f90wrap_inputdata__get__ip_dampfile(f90wrap_ip_dampfile)
    use inputdata, only: inputdata_ip_dampfile => ip_dampfile
    implicit none
    character(35), intent(out) :: f90wrap_ip_dampfile
    
    f90wrap_ip_dampfile = inputdata_ip_dampfile
end subroutine f90wrap_inputdata__get__ip_dampfile

subroutine f90wrap_inputdata__set__ip_dampfile(f90wrap_ip_dampfile)
    use inputdata, only: inputdata_ip_dampfile => ip_dampfile
    implicit none
    character(35), intent(in) :: f90wrap_ip_dampfile
    
    inputdata_ip_dampfile = f90wrap_ip_dampfile
end subroutine f90wrap_inputdata__set__ip_dampfile

subroutine f90wrap_inputdata__get__mseed(f90wrap_mseed)
    use inputdata, only: inputdata_mseed => mseed
    implicit none
    integer, intent(out) :: f90wrap_mseed
    
    f90wrap_mseed = inputdata_mseed
end subroutine f90wrap_inputdata__get__mseed

subroutine f90wrap_inputdata__set__mseed(f90wrap_mseed)
    use inputdata, only: inputdata_mseed => mseed
    implicit none
    integer, intent(in) :: f90wrap_mseed
    
    inputdata_mseed = f90wrap_mseed
end subroutine f90wrap_inputdata__set__mseed

subroutine f90wrap_inputdata__get__roteul(f90wrap_roteul)
    use inputdata, only: inputdata_roteul => roteul
    implicit none
    integer, intent(out) :: f90wrap_roteul
    
    f90wrap_roteul = inputdata_roteul
end subroutine f90wrap_inputdata__get__roteul

subroutine f90wrap_inputdata__set__roteul(f90wrap_roteul)
    use inputdata, only: inputdata_roteul => roteul
    implicit none
    integer, intent(in) :: f90wrap_roteul
    
    inputdata_roteul = f90wrap_roteul
end subroutine f90wrap_inputdata__set__roteul

subroutine f90wrap_inputdata__get__initmag(f90wrap_initmag)
    use inputdata, only: inputdata_initmag => initmag
    implicit none
    integer, intent(out) :: f90wrap_initmag
    
    f90wrap_initmag = inputdata_initmag
end subroutine f90wrap_inputdata__get__initmag

subroutine f90wrap_inputdata__set__initmag(f90wrap_initmag)
    use inputdata, only: inputdata_initmag => initmag
    implicit none
    integer, intent(in) :: f90wrap_initmag
    
    inputdata_initmag = f90wrap_initmag
end subroutine f90wrap_inputdata__set__initmag

subroutine f90wrap_inputdata__get__initneigh(f90wrap_initneigh)
    use inputdata, only: inputdata_initneigh => initneigh
    implicit none
    integer, intent(out) :: f90wrap_initneigh
    
    f90wrap_initneigh = inputdata_initneigh
end subroutine f90wrap_inputdata__get__initneigh

subroutine f90wrap_inputdata__set__initneigh(f90wrap_initneigh)
    use inputdata, only: inputdata_initneigh => initneigh
    implicit none
    integer, intent(in) :: f90wrap_initneigh
    
    inputdata_initneigh = f90wrap_initneigh
end subroutine f90wrap_inputdata__set__initneigh

subroutine f90wrap_inputdata__get__phi0(f90wrap_phi0)
    use inputdata, only: inputdata_phi0 => phi0
    implicit none
    real(8), intent(out) :: f90wrap_phi0
    
    f90wrap_phi0 = inputdata_phi0
end subroutine f90wrap_inputdata__get__phi0

subroutine f90wrap_inputdata__set__phi0(f90wrap_phi0)
    use inputdata, only: inputdata_phi0 => phi0
    implicit none
    real(8), intent(in) :: f90wrap_phi0
    
    inputdata_phi0 = f90wrap_phi0
end subroutine f90wrap_inputdata__set__phi0

subroutine f90wrap_inputdata__get__theta0(f90wrap_theta0)
    use inputdata, only: inputdata_theta0 => theta0
    implicit none
    real(8), intent(out) :: f90wrap_theta0
    
    f90wrap_theta0 = inputdata_theta0
end subroutine f90wrap_inputdata__get__theta0

subroutine f90wrap_inputdata__set__theta0(f90wrap_theta0)
    use inputdata, only: inputdata_theta0 => theta0
    implicit none
    real(8), intent(in) :: f90wrap_theta0
    
    inputdata_theta0 = f90wrap_theta0
end subroutine f90wrap_inputdata__set__theta0

subroutine f90wrap_inputdata__get__mavg0(f90wrap_mavg0)
    use inputdata, only: inputdata_mavg0 => mavg0
    implicit none
    real(8), intent(out) :: f90wrap_mavg0
    
    f90wrap_mavg0 = inputdata_mavg0
end subroutine f90wrap_inputdata__get__mavg0

subroutine f90wrap_inputdata__set__mavg0(f90wrap_mavg0)
    use inputdata, only: inputdata_mavg0 => mavg0
    implicit none
    real(8), intent(in) :: f90wrap_mavg0
    
    inputdata_mavg0 = f90wrap_mavg0
end subroutine f90wrap_inputdata__set__mavg0

subroutine f90wrap_inputdata__get__initimp(f90wrap_initimp)
    use inputdata, only: inputdata_initimp => initimp
    implicit none
    real(8), intent(out) :: f90wrap_initimp
    
    f90wrap_initimp = inputdata_initimp
end subroutine f90wrap_inputdata__get__initimp

subroutine f90wrap_inputdata__set__initimp(f90wrap_initimp)
    use inputdata, only: inputdata_initimp => initimp
    implicit none
    real(8), intent(in) :: f90wrap_initimp
    
    inputdata_initimp = f90wrap_initimp
end subroutine f90wrap_inputdata__set__initimp

subroutine f90wrap_inputdata__get__initconc(f90wrap_initconc)
    use inputdata, only: inputdata_initconc => initconc
    implicit none
    real(8), intent(out) :: f90wrap_initconc
    
    f90wrap_initconc = inputdata_initconc
end subroutine f90wrap_inputdata__get__initconc

subroutine f90wrap_inputdata__set__initconc(f90wrap_initconc)
    use inputdata, only: inputdata_initconc => initconc
    implicit none
    real(8), intent(in) :: f90wrap_initconc
    
    inputdata_initconc = f90wrap_initconc
end subroutine f90wrap_inputdata__set__initconc

subroutine f90wrap_inputdata__get__initrotang(f90wrap_initrotang)
    use inputdata, only: inputdata_initrotang => initrotang
    implicit none
    real(8), intent(out) :: f90wrap_initrotang
    
    f90wrap_initrotang = inputdata_initrotang
end subroutine f90wrap_inputdata__get__initrotang

subroutine f90wrap_inputdata__set__initrotang(f90wrap_initrotang)
    use inputdata, only: inputdata_initrotang => initrotang
    implicit none
    real(8), intent(in) :: f90wrap_initrotang
    
    inputdata_initrotang = f90wrap_initrotang
end subroutine f90wrap_inputdata__set__initrotang

subroutine f90wrap_inputdata__array__rotang(dummy_this, nd, dtype, dshape, dloc)
    use parameters
    use profiling
    use inputdatatype
    use inputdata, only: inputdata_rotang => rotang
    use, intrinsic :: iso_c_binding, only : c_int
    implicit none
    integer, intent(in) :: dummy_this(2)
    integer(c_int), intent(out) :: nd
    integer(c_int), intent(out) :: dtype
    integer(c_int), dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 1
    dtype = 12
    dshape(1:1) = shape(inputdata_rotang)
    dloc = loc(inputdata_rotang)
end subroutine f90wrap_inputdata__array__rotang

subroutine f90wrap_inputdata__array__initrotvec(dummy_this, nd, dtype, dshape, dloc)
    use parameters
    use profiling
    use inputdatatype
    use inputdata, only: inputdata_initrotvec => initrotvec
    use, intrinsic :: iso_c_binding, only : c_int
    implicit none
    integer, intent(in) :: dummy_this(2)
    integer(c_int), intent(out) :: nd
    integer(c_int), intent(out) :: dtype
    integer(c_int), dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 1
    dtype = 12
    dshape(1:1) = shape(inputdata_initrotvec)
    dloc = loc(inputdata_initrotvec)
end subroutine f90wrap_inputdata__array__initrotvec

subroutine f90wrap_inputdata__array__initpropvec(dummy_this, nd, dtype, dshape, dloc)
    use parameters
    use profiling
    use inputdatatype
    use inputdata, only: inputdata_initpropvec => initpropvec
    use, intrinsic :: iso_c_binding, only : c_int
    implicit none
    integer, intent(in) :: dummy_this(2)
    integer(c_int), intent(out) :: nd
    integer(c_int), intent(out) :: dtype
    integer(c_int), dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 1
    dtype = 12
    dshape(1:1) = shape(inputdata_initpropvec)
    dloc = loc(inputdata_initpropvec)
end subroutine f90wrap_inputdata__array__initpropvec

subroutine f90wrap_inputdata__get__initexc(f90wrap_initexc)
    use inputdata, only: inputdata_initexc => initexc
    implicit none
    character(1), intent(out) :: f90wrap_initexc
    
    f90wrap_initexc = inputdata_initexc
end subroutine f90wrap_inputdata__get__initexc

subroutine f90wrap_inputdata__set__initexc(f90wrap_initexc)
    use inputdata, only: inputdata_initexc => initexc
    implicit none
    character(1), intent(in) :: f90wrap_initexc
    
    inputdata_initexc = f90wrap_initexc
end subroutine f90wrap_inputdata__set__initexc

subroutine f90wrap_inputdata__get__restartfile(f90wrap_restartfile)
    use inputdata, only: inputdata_restartfile => restartfile
    implicit none
    character(35), intent(out) :: f90wrap_restartfile
    
    f90wrap_restartfile = inputdata_restartfile
end subroutine f90wrap_inputdata__get__restartfile

subroutine f90wrap_inputdata__set__restartfile(f90wrap_restartfile)
    use inputdata, only: inputdata_restartfile => restartfile
    implicit none
    character(35), intent(in) :: f90wrap_restartfile
    
    inputdata_restartfile = f90wrap_restartfile
end subroutine f90wrap_inputdata__set__restartfile

subroutine f90wrap_inputdata__get__demagvol(f90wrap_demagvol)
    use inputdata, only: inputdata_demagvol => demagvol
    implicit none
    real(8), intent(out) :: f90wrap_demagvol
    
    f90wrap_demagvol = inputdata_demagvol
end subroutine f90wrap_inputdata__get__demagvol

subroutine f90wrap_inputdata__set__demagvol(f90wrap_demagvol)
    use inputdata, only: inputdata_demagvol => demagvol
    implicit none
    real(8), intent(in) :: f90wrap_demagvol
    
    inputdata_demagvol = f90wrap_demagvol
end subroutine f90wrap_inputdata__set__demagvol

subroutine f90wrap_inputdata__get__demag(f90wrap_demag)
    use inputdata, only: inputdata_demag => demag
    implicit none
    character(1), intent(out) :: f90wrap_demag
    
    f90wrap_demag = inputdata_demag
end subroutine f90wrap_inputdata__get__demag

subroutine f90wrap_inputdata__set__demag(f90wrap_demag)
    use inputdata, only: inputdata_demag => demag
    implicit none
    character(1), intent(in) :: f90wrap_demag
    
    inputdata_demag = f90wrap_demag
end subroutine f90wrap_inputdata__set__demag

subroutine f90wrap_inputdata__get__demag1(f90wrap_demag1)
    use inputdata, only: inputdata_demag1 => demag1
    implicit none
    character(1), intent(out) :: f90wrap_demag1
    
    f90wrap_demag1 = inputdata_demag1
end subroutine f90wrap_inputdata__get__demag1

subroutine f90wrap_inputdata__set__demag1(f90wrap_demag1)
    use inputdata, only: inputdata_demag1 => demag1
    implicit none
    character(1), intent(in) :: f90wrap_demag1
    
    inputdata_demag1 = f90wrap_demag1
end subroutine f90wrap_inputdata__set__demag1

subroutine f90wrap_inputdata__get__demag2(f90wrap_demag2)
    use inputdata, only: inputdata_demag2 => demag2
    implicit none
    character(1), intent(out) :: f90wrap_demag2
    
    f90wrap_demag2 = inputdata_demag2
end subroutine f90wrap_inputdata__get__demag2

subroutine f90wrap_inputdata__set__demag2(f90wrap_demag2)
    use inputdata, only: inputdata_demag2 => demag2
    implicit none
    character(1), intent(in) :: f90wrap_demag2
    
    inputdata_demag2 = f90wrap_demag2
end subroutine f90wrap_inputdata__set__demag2

subroutine f90wrap_inputdata__get__demag3(f90wrap_demag3)
    use inputdata, only: inputdata_demag3 => demag3
    implicit none
    character(1), intent(out) :: f90wrap_demag3
    
    f90wrap_demag3 = inputdata_demag3
end subroutine f90wrap_inputdata__get__demag3

subroutine f90wrap_inputdata__set__demag3(f90wrap_demag3)
    use inputdata, only: inputdata_demag3 => demag3
    implicit none
    character(1), intent(in) :: f90wrap_demag3
    
    inputdata_demag3 = f90wrap_demag3
end subroutine f90wrap_inputdata__set__demag3

subroutine f90wrap_inputdata__array__sitenatomfld(dummy_this, nd, dtype, dshape, dloc)
    use parameters
    use profiling
    use inputdatatype
    use inputdata, only: inputdata_sitenatomfld => sitenatomfld
    use, intrinsic :: iso_c_binding, only : c_int
    implicit none
    integer, intent(in) :: dummy_this(2)
    integer(c_int), intent(out) :: nd
    integer(c_int), intent(out) :: dtype
    integer(c_int), dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 2
    dtype = 12
    if (allocated(inputdata_sitenatomfld)) then
        dshape(1:2) = shape(inputdata_sitenatomfld)
        dloc = loc(inputdata_sitenatomfld)
    else
        dloc = 0
    end if
end subroutine f90wrap_inputdata__array__sitenatomfld

subroutine f90wrap_inputdata__get__do_bpulse(f90wrap_do_bpulse)
    use inputdata, only: inputdata_do_bpulse => do_bpulse
    implicit none
    integer, intent(out) :: f90wrap_do_bpulse
    
    f90wrap_do_bpulse = inputdata_do_bpulse
end subroutine f90wrap_inputdata__get__do_bpulse

subroutine f90wrap_inputdata__set__do_bpulse(f90wrap_do_bpulse)
    use inputdata, only: inputdata_do_bpulse => do_bpulse
    implicit none
    integer, intent(in) :: f90wrap_do_bpulse
    
    inputdata_do_bpulse = f90wrap_do_bpulse
end subroutine f90wrap_inputdata__set__do_bpulse

subroutine f90wrap_inputdata__get__locfield(f90wrap_locfield)
    use inputdata, only: inputdata_locfield => locfield
    implicit none
    character(1), intent(out) :: f90wrap_locfield
    
    f90wrap_locfield = inputdata_locfield
end subroutine f90wrap_inputdata__get__locfield

subroutine f90wrap_inputdata__set__locfield(f90wrap_locfield)
    use inputdata, only: inputdata_locfield => locfield
    implicit none
    character(1), intent(in) :: f90wrap_locfield
    
    inputdata_locfield = f90wrap_locfield
end subroutine f90wrap_inputdata__set__locfield

subroutine f90wrap_inputdata__get__bpulsefile(f90wrap_bpulsefile)
    use inputdata, only: inputdata_bpulsefile => bpulsefile
    implicit none
    character(35), intent(out) :: f90wrap_bpulsefile
    
    f90wrap_bpulsefile = inputdata_bpulsefile
end subroutine f90wrap_inputdata__get__bpulsefile

subroutine f90wrap_inputdata__set__bpulsefile(f90wrap_bpulsefile)
    use inputdata, only: inputdata_bpulsefile => bpulsefile
    implicit none
    character(35), intent(in) :: f90wrap_bpulsefile
    
    inputdata_bpulsefile = f90wrap_bpulsefile
end subroutine f90wrap_inputdata__set__bpulsefile

subroutine f90wrap_inputdata__get__siteatomfile(f90wrap_siteatomfile)
    use inputdata, only: inputdata_siteatomfile => siteatomfile
    implicit none
    character(35), intent(out) :: f90wrap_siteatomfile
    
    f90wrap_siteatomfile = inputdata_siteatomfile
end subroutine f90wrap_inputdata__get__siteatomfile

subroutine f90wrap_inputdata__set__siteatomfile(f90wrap_siteatomfile)
    use inputdata, only: inputdata_siteatomfile => siteatomfile
    implicit none
    character(35), intent(in) :: f90wrap_siteatomfile
    
    inputdata_siteatomfile = f90wrap_siteatomfile
end subroutine f90wrap_inputdata__set__siteatomfile

subroutine f90wrap_inputdata__get__locfieldfile(f90wrap_locfieldfile)
    use inputdata, only: inputdata_locfieldfile => locfieldfile
    implicit none
    character(35), intent(out) :: f90wrap_locfieldfile
    
    f90wrap_locfieldfile = inputdata_locfieldfile
end subroutine f90wrap_inputdata__get__locfieldfile

subroutine f90wrap_inputdata__set__locfieldfile(f90wrap_locfieldfile)
    use inputdata, only: inputdata_locfieldfile => locfieldfile
    implicit none
    character(35), intent(in) :: f90wrap_locfieldfile
    
    inputdata_locfieldfile = f90wrap_locfieldfile
end subroutine f90wrap_inputdata__set__locfieldfile

subroutine f90wrap_inputdata__get__conf_num(f90wrap_conf_num)
    use inputdata, only: inputdata_conf_num => conf_num
    implicit none
    integer, intent(out) :: f90wrap_conf_num
    
    f90wrap_conf_num = inputdata_conf_num
end subroutine f90wrap_inputdata__get__conf_num

subroutine f90wrap_inputdata__set__conf_num(f90wrap_conf_num)
    use inputdata, only: inputdata_conf_num => conf_num
    implicit none
    integer, intent(in) :: f90wrap_conf_num
    
    inputdata_conf_num = f90wrap_conf_num
end subroutine f90wrap_inputdata__set__conf_num

subroutine f90wrap_inputdata__get__gsconf_num(f90wrap_gsconf_num)
    use inputdata, only: inputdata_gsconf_num => gsconf_num
    implicit none
    integer, intent(out) :: f90wrap_gsconf_num
    
    f90wrap_gsconf_num = inputdata_gsconf_num
end subroutine f90wrap_inputdata__get__gsconf_num

subroutine f90wrap_inputdata__set__gsconf_num(f90wrap_gsconf_num)
    use inputdata, only: inputdata_gsconf_num => gsconf_num
    implicit none
    integer, intent(in) :: f90wrap_gsconf_num
    
    inputdata_gsconf_num = f90wrap_gsconf_num
end subroutine f90wrap_inputdata__set__gsconf_num

subroutine f90wrap_inputdata__get__lsf_metric(f90wrap_lsf_metric)
    use inputdata, only: inputdata_lsf_metric => lsf_metric
    implicit none
    integer, intent(out) :: f90wrap_lsf_metric
    
    f90wrap_lsf_metric = inputdata_lsf_metric
end subroutine f90wrap_inputdata__get__lsf_metric

subroutine f90wrap_inputdata__set__lsf_metric(f90wrap_lsf_metric)
    use inputdata, only: inputdata_lsf_metric => lsf_metric
    implicit none
    integer, intent(in) :: f90wrap_lsf_metric
    
    inputdata_lsf_metric = f90wrap_lsf_metric
end subroutine f90wrap_inputdata__set__lsf_metric

subroutine f90wrap_inputdata__get__lsf_window(f90wrap_lsf_window)
    use inputdata, only: inputdata_lsf_window => lsf_window
    implicit none
    real(8), intent(out) :: f90wrap_lsf_window
    
    f90wrap_lsf_window = inputdata_lsf_window
end subroutine f90wrap_inputdata__get__lsf_window

subroutine f90wrap_inputdata__set__lsf_window(f90wrap_lsf_window)
    use inputdata, only: inputdata_lsf_window => lsf_window
    implicit none
    real(8), intent(in) :: f90wrap_lsf_window
    
    inputdata_lsf_window = f90wrap_lsf_window
end subroutine f90wrap_inputdata__set__lsf_window

subroutine f90wrap_inputdata__get__do_lsf(f90wrap_do_lsf)
    use inputdata, only: inputdata_do_lsf => do_lsf
    implicit none
    character(1), intent(out) :: f90wrap_do_lsf
    
    f90wrap_do_lsf = inputdata_do_lsf
end subroutine f90wrap_inputdata__get__do_lsf

subroutine f90wrap_inputdata__set__do_lsf(f90wrap_do_lsf)
    use inputdata, only: inputdata_do_lsf => do_lsf
    implicit none
    character(1), intent(in) :: f90wrap_do_lsf
    
    inputdata_do_lsf = f90wrap_do_lsf
end subroutine f90wrap_inputdata__set__do_lsf

subroutine f90wrap_inputdata__get__lsf_field(f90wrap_lsf_field)
    use inputdata, only: inputdata_lsf_field => lsf_field
    implicit none
    character(1), intent(out) :: f90wrap_lsf_field
    
    f90wrap_lsf_field = inputdata_lsf_field
end subroutine f90wrap_inputdata__get__lsf_field

subroutine f90wrap_inputdata__set__lsf_field(f90wrap_lsf_field)
    use inputdata, only: inputdata_lsf_field => lsf_field
    implicit none
    character(1), intent(in) :: f90wrap_lsf_field
    
    inputdata_lsf_field = f90wrap_lsf_field
end subroutine f90wrap_inputdata__set__lsf_field

subroutine f90wrap_inputdata__get__lsf_interpolate(f90wrap_lsf_interpolate)
    use inputdata, only: inputdata_lsf_interpolate => lsf_interpolate
    implicit none
    character(1), intent(out) :: f90wrap_lsf_interpolate
    
    f90wrap_lsf_interpolate = inputdata_lsf_interpolate
end subroutine f90wrap_inputdata__get__lsf_interpolate

subroutine f90wrap_inputdata__set__lsf_interpolate(f90wrap_lsf_interpolate)
    use inputdata, only: inputdata_lsf_interpolate => lsf_interpolate
    implicit none
    character(1), intent(in) :: f90wrap_lsf_interpolate
    
    inputdata_lsf_interpolate = f90wrap_lsf_interpolate
end subroutine f90wrap_inputdata__set__lsf_interpolate

subroutine f90wrap_inputdata__get__lsffile(f90wrap_lsffile)
    use inputdata, only: inputdata_lsffile => lsffile
    implicit none
    character(35), intent(out) :: f90wrap_lsffile
    
    f90wrap_lsffile = inputdata_lsffile
end subroutine f90wrap_inputdata__get__lsffile

subroutine f90wrap_inputdata__set__lsffile(f90wrap_lsffile)
    use inputdata, only: inputdata_lsffile => lsffile
    implicit none
    character(35), intent(in) :: f90wrap_lsffile
    
    inputdata_lsffile = f90wrap_lsffile
end subroutine f90wrap_inputdata__set__lsffile

subroutine f90wrap_inputdata__get__spintemp_step(f90wrap_spintemp_step)
    use inputdata, only: inputdata_spintemp_step => spintemp_step
    implicit none
    integer, intent(out) :: f90wrap_spintemp_step
    
    f90wrap_spintemp_step = inputdata_spintemp_step
end subroutine f90wrap_inputdata__get__spintemp_step

subroutine f90wrap_inputdata__set__spintemp_step(f90wrap_spintemp_step)
    use inputdata, only: inputdata_spintemp_step => spintemp_step
    implicit none
    integer, intent(in) :: f90wrap_spintemp_step
    
    inputdata_spintemp_step = f90wrap_spintemp_step
end subroutine f90wrap_inputdata__set__spintemp_step

subroutine f90wrap_inputdata__get__logsamp(f90wrap_logsamp)
    use inputdata, only: inputdata_logsamp => logsamp
    implicit none
    character(1), intent(out) :: f90wrap_logsamp
    
    f90wrap_logsamp = inputdata_logsamp
end subroutine f90wrap_inputdata__get__logsamp

subroutine f90wrap_inputdata__set__logsamp(f90wrap_logsamp)
    use inputdata, only: inputdata_logsamp => logsamp
    implicit none
    character(1), intent(in) :: f90wrap_logsamp
    
    inputdata_logsamp = f90wrap_logsamp
end subroutine f90wrap_inputdata__set__logsamp

subroutine f90wrap_inputdata__get__do_spintemp(f90wrap_do_spintemp)
    use inputdata, only: inputdata_do_spintemp => do_spintemp
    implicit none
    character(1), intent(out) :: f90wrap_do_spintemp
    
    f90wrap_do_spintemp = inputdata_do_spintemp
end subroutine f90wrap_inputdata__get__do_spintemp

subroutine f90wrap_inputdata__set__do_spintemp(f90wrap_do_spintemp)
    use inputdata, only: inputdata_do_spintemp => do_spintemp
    implicit none
    character(1), intent(in) :: f90wrap_do_spintemp
    
    inputdata_do_spintemp = f90wrap_do_spintemp
end subroutine f90wrap_inputdata__set__do_spintemp

subroutine f90wrap_inputdata__get__real_time_measure(f90wrap_real_time_measure)
    use inputdata, only: inputdata_real_time_measure => real_time_measure
    implicit none
    character(1), intent(out) :: f90wrap_real_time_measure
    
    f90wrap_real_time_measure = inputdata_real_time_measure
end subroutine f90wrap_inputdata__get__real_time_measure

subroutine f90wrap_inputdata__set__real_time_measure(f90wrap_real_time_measure)
    use inputdata, only: inputdata_real_time_measure => real_time_measure
    implicit none
    character(1), intent(in) :: f90wrap_real_time_measure
    
    inputdata_real_time_measure = f90wrap_real_time_measure
end subroutine f90wrap_inputdata__set__real_time_measure

subroutine f90wrap_inputdata__get__Nchmax(f90wrap_Nchmax)
    use inputdata, only: inputdata_Nchmax => Nchmax
    implicit none
    integer, intent(out) :: f90wrap_Nchmax
    
    f90wrap_Nchmax = inputdata_Nchmax
end subroutine f90wrap_inputdata__get__Nchmax

subroutine f90wrap_inputdata__set__Nchmax(f90wrap_Nchmax)
    use inputdata, only: inputdata_Nchmax => Nchmax
    implicit none
    integer, intent(in) :: f90wrap_Nchmax
    
    inputdata_Nchmax = f90wrap_Nchmax
end subroutine f90wrap_inputdata__set__Nchmax

subroutine f90wrap_inputdata__get__do_ralloy(f90wrap_do_ralloy)
    use inputdata, only: inputdata_do_ralloy => do_ralloy
    implicit none
    integer, intent(out) :: f90wrap_do_ralloy
    
    f90wrap_do_ralloy = inputdata_do_ralloy
end subroutine f90wrap_inputdata__get__do_ralloy

subroutine f90wrap_inputdata__set__do_ralloy(f90wrap_do_ralloy)
    use inputdata, only: inputdata_do_ralloy => do_ralloy
    implicit none
    integer, intent(in) :: f90wrap_do_ralloy
    
    inputdata_do_ralloy = f90wrap_do_ralloy
end subroutine f90wrap_inputdata__set__do_ralloy

subroutine f90wrap_inputdata__array__Nch(dummy_this, nd, dtype, dshape, dloc)
    use parameters
    use profiling
    use inputdatatype
    use inputdata, only: inputdata_nch => nch
    use, intrinsic :: iso_c_binding, only : c_int
    implicit none
    integer, intent(in) :: dummy_this(2)
    integer(c_int), intent(out) :: nd
    integer(c_int), intent(out) :: dtype
    integer(c_int), dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 1
    dtype = 5
    if (allocated(inputdata_Nch)) then
        dshape(1:1) = shape(inputdata_Nch)
        dloc = loc(inputdata_Nch)
    else
        dloc = 0
    end if
end subroutine f90wrap_inputdata__array__Nch

subroutine f90wrap_inputdata__array__achtype_ch(dummy_this, nd, dtype, dshape, dloc)
    use parameters
    use profiling
    use inputdatatype
    use inputdata, only: inputdata_achtype_ch => achtype_ch
    use, intrinsic :: iso_c_binding, only : c_int
    implicit none
    integer, intent(in) :: dummy_this(2)
    integer(c_int), intent(out) :: nd
    integer(c_int), intent(out) :: dtype
    integer(c_int), dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 2
    dtype = 5
    if (allocated(inputdata_achtype_ch)) then
        dshape(1:2) = shape(inputdata_achtype_ch)
        dloc = loc(inputdata_achtype_ch)
    else
        dloc = 0
    end if
end subroutine f90wrap_inputdata__array__achtype_ch

subroutine f90wrap_inputdata__array__chconc(dummy_this, nd, dtype, dshape, dloc)
    use parameters
    use profiling
    use inputdatatype
    use inputdata, only: inputdata_chconc => chconc
    use, intrinsic :: iso_c_binding, only : c_int
    implicit none
    integer, intent(in) :: dummy_this(2)
    integer(c_int), intent(out) :: nd
    integer(c_int), intent(out) :: dtype
    integer(c_int), dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 2
    dtype = 12
    if (allocated(inputdata_chconc)) then
        dshape(1:2) = shape(inputdata_chconc)
        dloc = loc(inputdata_chconc)
    else
        dloc = 0
    end if
end subroutine f90wrap_inputdata__array__chconc

subroutine f90wrap_inputdata__get__oldformat(f90wrap_oldformat)
    use inputdata, only: inputdata_oldformat => oldformat
    implicit none
    logical, intent(out) :: f90wrap_oldformat
    
    f90wrap_oldformat = inputdata_oldformat
end subroutine f90wrap_inputdata__get__oldformat

subroutine f90wrap_inputdata__set__oldformat(f90wrap_oldformat)
    use inputdata, only: inputdata_oldformat => oldformat
    implicit none
    logical, intent(in) :: f90wrap_oldformat
    
    inputdata_oldformat = f90wrap_oldformat
end subroutine f90wrap_inputdata__set__oldformat

subroutine f90wrap_inputdata__get__rngpol(f90wrap_rngpol)
    use inputdata, only: inputdata_rngpol => rngpol
    implicit none
    character(1), intent(out) :: f90wrap_rngpol
    
    f90wrap_rngpol = inputdata_rngpol
end subroutine f90wrap_inputdata__get__rngpol

subroutine f90wrap_inputdata__set__rngpol(f90wrap_rngpol)
    use inputdata, only: inputdata_rngpol => rngpol
    implicit none
    character(1), intent(in) :: f90wrap_rngpol
    
    inputdata_rngpol = f90wrap_rngpol
end subroutine f90wrap_inputdata__set__rngpol

subroutine f90wrap_inputdata__get__ziggurat(f90wrap_ziggurat)
    use inputdata, only: inputdata_ziggurat => ziggurat
    implicit none
    character(1), intent(out) :: f90wrap_ziggurat
    
    f90wrap_ziggurat = inputdata_ziggurat
end subroutine f90wrap_inputdata__get__ziggurat

subroutine f90wrap_inputdata__set__ziggurat(f90wrap_ziggurat)
    use inputdata, only: inputdata_ziggurat => ziggurat
    implicit none
    character(1), intent(in) :: f90wrap_ziggurat
    
    inputdata_ziggurat = f90wrap_ziggurat
end subroutine f90wrap_inputdata__set__ziggurat

subroutine f90wrap_inputdata__get__gpu_mode(f90wrap_gpu_mode)
    use inputdata, only: inputdata_gpu_mode => gpu_mode
    implicit none
    integer, intent(out) :: f90wrap_gpu_mode
    
    f90wrap_gpu_mode = inputdata_gpu_mode
end subroutine f90wrap_inputdata__get__gpu_mode

subroutine f90wrap_inputdata__set__gpu_mode(f90wrap_gpu_mode)
    use inputdata, only: inputdata_gpu_mode => gpu_mode
    implicit none
    integer, intent(in) :: f90wrap_gpu_mode
    
    inputdata_gpu_mode = f90wrap_gpu_mode
end subroutine f90wrap_inputdata__set__gpu_mode

subroutine f90wrap_inputdata__get__gpu_rng(f90wrap_gpu_rng)
    use inputdata, only: inputdata_gpu_rng => gpu_rng
    implicit none
    integer, intent(out) :: f90wrap_gpu_rng
    
    f90wrap_gpu_rng = inputdata_gpu_rng
end subroutine f90wrap_inputdata__get__gpu_rng

subroutine f90wrap_inputdata__set__gpu_rng(f90wrap_gpu_rng)
    use inputdata, only: inputdata_gpu_rng => gpu_rng
    implicit none
    integer, intent(in) :: f90wrap_gpu_rng
    
    inputdata_gpu_rng = f90wrap_gpu_rng
end subroutine f90wrap_inputdata__set__gpu_rng

subroutine f90wrap_inputdata__get__gpu_rng_seed(f90wrap_gpu_rng_seed)
    use inputdata, only: inputdata_gpu_rng_seed => gpu_rng_seed
    implicit none
    integer, intent(out) :: f90wrap_gpu_rng_seed
    
    f90wrap_gpu_rng_seed = inputdata_gpu_rng_seed
end subroutine f90wrap_inputdata__get__gpu_rng_seed

subroutine f90wrap_inputdata__set__gpu_rng_seed(f90wrap_gpu_rng_seed)
    use inputdata, only: inputdata_gpu_rng_seed => gpu_rng_seed
    implicit none
    integer, intent(in) :: f90wrap_gpu_rng_seed
    
    inputdata_gpu_rng_seed = f90wrap_gpu_rng_seed
end subroutine f90wrap_inputdata__set__gpu_rng_seed

subroutine f90wrap_inputdata__get__prn_ovf(f90wrap_prn_ovf)
    use inputdata, only: inputdata_prn_ovf => prn_ovf
    implicit none
    character(1), intent(out) :: f90wrap_prn_ovf
    
    f90wrap_prn_ovf = inputdata_prn_ovf
end subroutine f90wrap_inputdata__get__prn_ovf

subroutine f90wrap_inputdata__set__prn_ovf(f90wrap_prn_ovf)
    use inputdata, only: inputdata_prn_ovf => prn_ovf
    implicit none
    character(1), intent(in) :: f90wrap_prn_ovf
    
    inputdata_prn_ovf = f90wrap_prn_ovf
end subroutine f90wrap_inputdata__set__prn_ovf

subroutine f90wrap_inputdata__get__read_ovf(f90wrap_read_ovf)
    use inputdata, only: inputdata_read_ovf => read_ovf
    implicit none
    character(1), intent(out) :: f90wrap_read_ovf
    
    f90wrap_read_ovf = inputdata_read_ovf
end subroutine f90wrap_inputdata__get__read_ovf

subroutine f90wrap_inputdata__set__read_ovf(f90wrap_read_ovf)
    use inputdata, only: inputdata_read_ovf => read_ovf
    implicit none
    character(1), intent(in) :: f90wrap_read_ovf
    
    inputdata_read_ovf = f90wrap_read_ovf
end subroutine f90wrap_inputdata__set__read_ovf

subroutine f90wrap_inputdata__get__do_multiscale(f90wrap_do_multiscale)
    use inputdata, only: inputdata_do_multiscale => do_multiscale
    implicit none
    logical, intent(out) :: f90wrap_do_multiscale
    
    f90wrap_do_multiscale = inputdata_do_multiscale
end subroutine f90wrap_inputdata__get__do_multiscale

subroutine f90wrap_inputdata__set__do_multiscale(f90wrap_do_multiscale)
    use inputdata, only: inputdata_do_multiscale => do_multiscale
    implicit none
    logical, intent(in) :: f90wrap_do_multiscale
    
    inputdata_do_multiscale = f90wrap_do_multiscale
end subroutine f90wrap_inputdata__set__do_multiscale

subroutine f90wrap_inputdata__get__do_prnmultiscale(f90wrap_do_prnmultiscale)
    use inputdata, only: inputdata_do_prnmultiscale => do_prnmultiscale
    implicit none
    logical, intent(out) :: f90wrap_do_prnmultiscale
    
    f90wrap_do_prnmultiscale = inputdata_do_prnmultiscale
end subroutine f90wrap_inputdata__get__do_prnmultiscale

subroutine f90wrap_inputdata__set__do_prnmultiscale(f90wrap_do_prnmultiscale)
    use inputdata, only: inputdata_do_prnmultiscale => do_prnmultiscale
    implicit none
    logical, intent(in) :: f90wrap_do_prnmultiscale
    
    inputdata_do_prnmultiscale = f90wrap_do_prnmultiscale
end subroutine f90wrap_inputdata__set__do_prnmultiscale

subroutine f90wrap_inputdata__get__multiscale_file_name(f90wrap_multiscale_file_name)
    use inputdata, only: inputdata_multiscale_file_name => multiscale_file_name
    implicit none
    character(260), intent(out) :: f90wrap_multiscale_file_name
    
    f90wrap_multiscale_file_name = inputdata_multiscale_file_name
end subroutine f90wrap_inputdata__get__multiscale_file_name

subroutine f90wrap_inputdata__set__multiscale_file_name(f90wrap_multiscale_file_name)
    use inputdata, only: inputdata_multiscale_file_name => multiscale_file_name
    implicit none
    character(260), intent(in) :: f90wrap_multiscale_file_name
    
    inputdata_multiscale_file_name = f90wrap_multiscale_file_name
end subroutine f90wrap_inputdata__set__multiscale_file_name

subroutine f90wrap_inputdata__get__multiscale_old_format(f90wrap_multiscale_old_format)
    use inputdata, only: inputdata_multiscale_old_format => multiscale_old_format
    implicit none
    character(1), intent(out) :: f90wrap_multiscale_old_format
    
    f90wrap_multiscale_old_format = inputdata_multiscale_old_format
end subroutine f90wrap_inputdata__get__multiscale_old_format

subroutine f90wrap_inputdata__set__multiscale_old_format(f90wrap_multiscale_old_format)
    use inputdata, only: inputdata_multiscale_old_format => multiscale_old_format
    implicit none
    character(1), intent(in) :: f90wrap_multiscale_old_format
    
    inputdata_multiscale_old_format = f90wrap_multiscale_old_format
end subroutine f90wrap_inputdata__set__multiscale_old_format

! End of module inputdata defined in file Input/inputdata.f90

