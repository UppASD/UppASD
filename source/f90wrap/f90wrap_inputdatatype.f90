! Module inputdatatype defined in file /Users/andersb/Jobb/UppASD_release/UppASD_release/source/Input/inputdatatype.f90

subroutine f90wrap_ham_inp_t__get__do_jtensor(this, f90wrap_do_jtensor)
    use inputdatatype, only: ham_inp_t
    implicit none
    type ham_inp_t_ptr_type
        type(ham_inp_t), pointer :: p => NULL()
    end type ham_inp_t_ptr_type
    integer, intent(in)   :: this(2)
    type(ham_inp_t_ptr_type) :: this_ptr
    integer, intent(out) :: f90wrap_do_jtensor
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_do_jtensor = this_ptr%p%do_jtensor
end subroutine f90wrap_ham_inp_t__get__do_jtensor

subroutine f90wrap_ham_inp_t__set__do_jtensor(this, f90wrap_do_jtensor)
    use inputdatatype, only: ham_inp_t
    implicit none
    type ham_inp_t_ptr_type
        type(ham_inp_t), pointer :: p => NULL()
    end type ham_inp_t_ptr_type
    integer, intent(in)   :: this(2)
    type(ham_inp_t_ptr_type) :: this_ptr
    integer, intent(in) :: f90wrap_do_jtensor
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%do_jtensor = f90wrap_do_jtensor
end subroutine f90wrap_ham_inp_t__set__do_jtensor

subroutine f90wrap_ham_inp_t__get__calc_jtensor(this, f90wrap_calc_jtensor)
    use inputdatatype, only: ham_inp_t
    implicit none
    type ham_inp_t_ptr_type
        type(ham_inp_t), pointer :: p => NULL()
    end type ham_inp_t_ptr_type
    integer, intent(in)   :: this(2)
    type(ham_inp_t_ptr_type) :: this_ptr
    logical, intent(out) :: f90wrap_calc_jtensor
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_calc_jtensor = this_ptr%p%calc_jtensor
end subroutine f90wrap_ham_inp_t__get__calc_jtensor

subroutine f90wrap_ham_inp_t__set__calc_jtensor(this, f90wrap_calc_jtensor)
    use inputdatatype, only: ham_inp_t
    implicit none
    type ham_inp_t_ptr_type
        type(ham_inp_t), pointer :: p => NULL()
    end type ham_inp_t_ptr_type
    integer, intent(in)   :: this(2)
    type(ham_inp_t_ptr_type) :: this_ptr
    logical, intent(in) :: f90wrap_calc_jtensor
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%calc_jtensor = f90wrap_calc_jtensor
end subroutine f90wrap_ham_inp_t__set__calc_jtensor

subroutine f90wrap_ham_inp_t__get__map_multiple(this, f90wrap_map_multiple)
    use inputdatatype, only: ham_inp_t
    implicit none
    type ham_inp_t_ptr_type
        type(ham_inp_t), pointer :: p => NULL()
    end type ham_inp_t_ptr_type
    integer, intent(in)   :: this(2)
    type(ham_inp_t_ptr_type) :: this_ptr
    logical, intent(out) :: f90wrap_map_multiple
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_map_multiple = this_ptr%p%map_multiple
end subroutine f90wrap_ham_inp_t__get__map_multiple

subroutine f90wrap_ham_inp_t__set__map_multiple(this, f90wrap_map_multiple)
    use inputdatatype, only: ham_inp_t
    implicit none
    type ham_inp_t_ptr_type
        type(ham_inp_t), pointer :: p => NULL()
    end type ham_inp_t_ptr_type
    integer, intent(in)   :: this(2)
    type(ham_inp_t_ptr_type) :: this_ptr
    logical, intent(in) :: f90wrap_map_multiple
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%map_multiple = f90wrap_map_multiple
end subroutine f90wrap_ham_inp_t__set__map_multiple

subroutine f90wrap_ham_inp_t__get__max_no_shells(this, f90wrap_max_no_shells)
    use inputdatatype, only: ham_inp_t
    implicit none
    type ham_inp_t_ptr_type
        type(ham_inp_t), pointer :: p => NULL()
    end type ham_inp_t_ptr_type
    integer, intent(in)   :: this(2)
    type(ham_inp_t_ptr_type) :: this_ptr
    integer, intent(out) :: f90wrap_max_no_shells
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_max_no_shells = this_ptr%p%max_no_shells
end subroutine f90wrap_ham_inp_t__get__max_no_shells

subroutine f90wrap_ham_inp_t__set__max_no_shells(this, f90wrap_max_no_shells)
    use inputdatatype, only: ham_inp_t
    implicit none
    type ham_inp_t_ptr_type
        type(ham_inp_t), pointer :: p => NULL()
    end type ham_inp_t_ptr_type
    integer, intent(in)   :: this(2)
    type(ham_inp_t_ptr_type) :: this_ptr
    integer, intent(in) :: f90wrap_max_no_shells
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%max_no_shells = f90wrap_max_no_shells
end subroutine f90wrap_ham_inp_t__set__max_no_shells

subroutine f90wrap_ham_inp_t__array__nn(this, nd, dtype, dshape, dloc)
    use inputdatatype, only: ham_inp_t
    use, intrinsic :: iso_c_binding, only : c_int
    implicit none
    type ham_inp_t_ptr_type
        type(ham_inp_t), pointer :: p => NULL()
    end type ham_inp_t_ptr_type
    integer(c_int), intent(in) :: this(2)
    type(ham_inp_t_ptr_type) :: this_ptr
    integer(c_int), intent(out) :: nd
    integer(c_int), intent(out) :: dtype
    integer(c_int), dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 1
    dtype = 5
    this_ptr = transfer(this, this_ptr)
    if (allocated(this_ptr%p%nn)) then
        dshape(1:1) = shape(this_ptr%p%nn)
        dloc = loc(this_ptr%p%nn)
    else
        dloc = 0
    end if
end subroutine f90wrap_ham_inp_t__array__nn

subroutine f90wrap_ham_inp_t__array__nntype(this, nd, dtype, dshape, dloc)
    use inputdatatype, only: ham_inp_t
    use, intrinsic :: iso_c_binding, only : c_int
    implicit none
    type ham_inp_t_ptr_type
        type(ham_inp_t), pointer :: p => NULL()
    end type ham_inp_t_ptr_type
    integer(c_int), intent(in) :: this(2)
    type(ham_inp_t_ptr_type) :: this_ptr
    integer(c_int), intent(out) :: nd
    integer(c_int), intent(out) :: dtype
    integer(c_int), dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 2
    dtype = 5
    this_ptr = transfer(this, this_ptr)
    if (allocated(this_ptr%p%nntype)) then
        dshape(1:2) = shape(this_ptr%p%nntype)
        dloc = loc(this_ptr%p%nntype)
    else
        dloc = 0
    end if
end subroutine f90wrap_ham_inp_t__array__nntype

subroutine f90wrap_ham_inp_t__array__redcoord(this, nd, dtype, dshape, dloc)
    use inputdatatype, only: ham_inp_t
    use, intrinsic :: iso_c_binding, only : c_int
    implicit none
    type ham_inp_t_ptr_type
        type(ham_inp_t), pointer :: p => NULL()
    end type ham_inp_t_ptr_type
    integer(c_int), intent(in) :: this(2)
    type(ham_inp_t_ptr_type) :: this_ptr
    integer(c_int), intent(out) :: nd
    integer(c_int), intent(out) :: dtype
    integer(c_int), dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 3
    dtype = 12
    this_ptr = transfer(this, this_ptr)
    if (allocated(this_ptr%p%redcoord)) then
        dshape(1:3) = shape(this_ptr%p%redcoord)
        dloc = loc(this_ptr%p%redcoord)
    else
        dloc = 0
    end if
end subroutine f90wrap_ham_inp_t__array__redcoord

subroutine f90wrap_ham_inp_t__array__jc(this, nd, dtype, dshape, dloc)
    use inputdatatype, only: ham_inp_t
    use, intrinsic :: iso_c_binding, only : c_int
    implicit none
    type ham_inp_t_ptr_type
        type(ham_inp_t), pointer :: p => NULL()
    end type ham_inp_t_ptr_type
    integer(c_int), intent(in) :: this(2)
    type(ham_inp_t_ptr_type) :: this_ptr
    integer(c_int), intent(out) :: nd
    integer(c_int), intent(out) :: dtype
    integer(c_int), dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 5
    dtype = 12
    this_ptr = transfer(this, this_ptr)
    if (allocated(this_ptr%p%jc)) then
        dshape(1:5) = shape(this_ptr%p%jc)
        dloc = loc(this_ptr%p%jc)
    else
        dloc = 0
    end if
end subroutine f90wrap_ham_inp_t__array__jc

subroutine f90wrap_ham_inp_t__array__jcD(this, nd, dtype, dshape, dloc)
    use inputdatatype, only: ham_inp_t
    use, intrinsic :: iso_c_binding, only : c_int
    implicit none
    type ham_inp_t_ptr_type
        type(ham_inp_t), pointer :: p => NULL()
    end type ham_inp_t_ptr_type
    integer(c_int), intent(in) :: this(2)
    type(ham_inp_t_ptr_type) :: this_ptr
    integer(c_int), intent(out) :: nd
    integer(c_int), intent(out) :: dtype
    integer(c_int), dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 5
    dtype = 12
    this_ptr = transfer(this, this_ptr)
    if (allocated(this_ptr%p%jcD)) then
        dshape(1:5) = shape(this_ptr%p%jcD)
        dloc = loc(this_ptr%p%jcD)
    else
        dloc = 0
    end if
end subroutine f90wrap_ham_inp_t__array__jcD

subroutine f90wrap_ham_inp_t__array__jc_tens(this, nd, dtype, dshape, dloc)
    use inputdatatype, only: ham_inp_t
    use, intrinsic :: iso_c_binding, only : c_int
    implicit none
    type ham_inp_t_ptr_type
        type(ham_inp_t), pointer :: p => NULL()
    end type ham_inp_t_ptr_type
    integer(c_int), intent(in) :: this(2)
    type(ham_inp_t_ptr_type) :: this_ptr
    integer(c_int), intent(out) :: nd
    integer(c_int), intent(out) :: dtype
    integer(c_int), dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 6
    dtype = 12
    this_ptr = transfer(this, this_ptr)
    if (allocated(this_ptr%p%jc_tens)) then
        dshape(1:6) = shape(this_ptr%p%jc_tens)
        dloc = loc(this_ptr%p%jc_tens)
    else
        dloc = 0
    end if
end subroutine f90wrap_ham_inp_t__array__jc_tens

subroutine f90wrap_ham_inp_t__get__exc_inter(this, f90wrap_exc_inter)
    use inputdatatype, only: ham_inp_t
    implicit none
    type ham_inp_t_ptr_type
        type(ham_inp_t), pointer :: p => NULL()
    end type ham_inp_t_ptr_type
    integer, intent(in)   :: this(2)
    type(ham_inp_t_ptr_type) :: this_ptr
    character(1), intent(out) :: f90wrap_exc_inter
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_exc_inter = this_ptr%p%exc_inter
end subroutine f90wrap_ham_inp_t__get__exc_inter

subroutine f90wrap_ham_inp_t__set__exc_inter(this, f90wrap_exc_inter)
    use inputdatatype, only: ham_inp_t
    implicit none
    type ham_inp_t_ptr_type
        type(ham_inp_t), pointer :: p => NULL()
    end type ham_inp_t_ptr_type
    integer, intent(in)   :: this(2)
    type(ham_inp_t_ptr_type) :: this_ptr
    character(1), intent(in) :: f90wrap_exc_inter
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%exc_inter = f90wrap_exc_inter
end subroutine f90wrap_ham_inp_t__set__exc_inter

subroutine f90wrap_ham_inp_t__get__jij_scale(this, f90wrap_jij_scale)
    use inputdatatype, only: ham_inp_t
    implicit none
    type ham_inp_t_ptr_type
        type(ham_inp_t), pointer :: p => NULL()
    end type ham_inp_t_ptr_type
    integer, intent(in)   :: this(2)
    type(ham_inp_t_ptr_type) :: this_ptr
    real(8), intent(out) :: f90wrap_jij_scale
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_jij_scale = this_ptr%p%jij_scale
end subroutine f90wrap_ham_inp_t__get__jij_scale

subroutine f90wrap_ham_inp_t__set__jij_scale(this, f90wrap_jij_scale)
    use inputdatatype, only: ham_inp_t
    implicit none
    type ham_inp_t_ptr_type
        type(ham_inp_t), pointer :: p => NULL()
    end type ham_inp_t_ptr_type
    integer, intent(in)   :: this(2)
    type(ham_inp_t_ptr_type) :: this_ptr
    real(8), intent(in) :: f90wrap_jij_scale
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%jij_scale = f90wrap_jij_scale
end subroutine f90wrap_ham_inp_t__set__jij_scale

subroutine f90wrap_ham_inp_t__get__ea_model(this, f90wrap_ea_model)
    use inputdatatype, only: ham_inp_t
    implicit none
    type ham_inp_t_ptr_type
        type(ham_inp_t), pointer :: p => NULL()
    end type ham_inp_t_ptr_type
    integer, intent(in)   :: this(2)
    type(ham_inp_t_ptr_type) :: this_ptr
    logical, intent(out) :: f90wrap_ea_model
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_ea_model = this_ptr%p%ea_model
end subroutine f90wrap_ham_inp_t__get__ea_model

subroutine f90wrap_ham_inp_t__set__ea_model(this, f90wrap_ea_model)
    use inputdatatype, only: ham_inp_t
    implicit none
    type ham_inp_t_ptr_type
        type(ham_inp_t), pointer :: p => NULL()
    end type ham_inp_t_ptr_type
    integer, intent(in)   :: this(2)
    type(ham_inp_t_ptr_type) :: this_ptr
    logical, intent(in) :: f90wrap_ea_model
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%ea_model = f90wrap_ea_model
end subroutine f90wrap_ham_inp_t__set__ea_model

subroutine f90wrap_ham_inp_t__get__ea_sigma(this, f90wrap_ea_sigma)
    use inputdatatype, only: ham_inp_t
    implicit none
    type ham_inp_t_ptr_type
        type(ham_inp_t), pointer :: p => NULL()
    end type ham_inp_t_ptr_type
    integer, intent(in)   :: this(2)
    type(ham_inp_t_ptr_type) :: this_ptr
    real(8), intent(out) :: f90wrap_ea_sigma
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_ea_sigma = this_ptr%p%ea_sigma
end subroutine f90wrap_ham_inp_t__get__ea_sigma

subroutine f90wrap_ham_inp_t__set__ea_sigma(this, f90wrap_ea_sigma)
    use inputdatatype, only: ham_inp_t
    implicit none
    type ham_inp_t_ptr_type
        type(ham_inp_t), pointer :: p => NULL()
    end type ham_inp_t_ptr_type
    integer, intent(in)   :: this(2)
    type(ham_inp_t_ptr_type) :: this_ptr
    real(8), intent(in) :: f90wrap_ea_sigma
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%ea_sigma = f90wrap_ea_sigma
end subroutine f90wrap_ham_inp_t__set__ea_sigma

subroutine f90wrap_ham_inp_t__get__do_anisotropy(this, f90wrap_do_anisotropy)
    use inputdatatype, only: ham_inp_t
    implicit none
    type ham_inp_t_ptr_type
        type(ham_inp_t), pointer :: p => NULL()
    end type ham_inp_t_ptr_type
    integer, intent(in)   :: this(2)
    type(ham_inp_t_ptr_type) :: this_ptr
    integer, intent(out) :: f90wrap_do_anisotropy
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_do_anisotropy = this_ptr%p%do_anisotropy
end subroutine f90wrap_ham_inp_t__get__do_anisotropy

subroutine f90wrap_ham_inp_t__set__do_anisotropy(this, f90wrap_do_anisotropy)
    use inputdatatype, only: ham_inp_t
    implicit none
    type ham_inp_t_ptr_type
        type(ham_inp_t), pointer :: p => NULL()
    end type ham_inp_t_ptr_type
    integer, intent(in)   :: this(2)
    type(ham_inp_t_ptr_type) :: this_ptr
    integer, intent(in) :: f90wrap_do_anisotropy
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%do_anisotropy = f90wrap_do_anisotropy
end subroutine f90wrap_ham_inp_t__set__do_anisotropy

subroutine f90wrap_ham_inp_t__get__random_anisotropy_density(this, f90wrap_random_anisotropy_density)
    use inputdatatype, only: ham_inp_t
    implicit none
    type ham_inp_t_ptr_type
        type(ham_inp_t), pointer :: p => NULL()
    end type ham_inp_t_ptr_type
    integer, intent(in)   :: this(2)
    type(ham_inp_t_ptr_type) :: this_ptr
    real(8), intent(out) :: f90wrap_random_anisotropy_density
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_random_anisotropy_density = this_ptr%p%random_anisotropy_density
end subroutine f90wrap_ham_inp_t__get__random_anisotropy_density

subroutine f90wrap_ham_inp_t__set__random_anisotropy_density(this, f90wrap_random_anisotropy_density)
    use inputdatatype, only: ham_inp_t
    implicit none
    type ham_inp_t_ptr_type
        type(ham_inp_t), pointer :: p => NULL()
    end type ham_inp_t_ptr_type
    integer, intent(in)   :: this(2)
    type(ham_inp_t_ptr_type) :: this_ptr
    real(8), intent(in) :: f90wrap_random_anisotropy_density
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%random_anisotropy_density = f90wrap_random_anisotropy_density
end subroutine f90wrap_ham_inp_t__set__random_anisotropy_density

subroutine f90wrap_ham_inp_t__get__random_anisotropy(this, f90wrap_random_anisotropy)
    use inputdatatype, only: ham_inp_t
    implicit none
    type ham_inp_t_ptr_type
        type(ham_inp_t), pointer :: p => NULL()
    end type ham_inp_t_ptr_type
    integer, intent(in)   :: this(2)
    type(ham_inp_t_ptr_type) :: this_ptr
    logical, intent(out) :: f90wrap_random_anisotropy
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_random_anisotropy = this_ptr%p%random_anisotropy
end subroutine f90wrap_ham_inp_t__get__random_anisotropy

subroutine f90wrap_ham_inp_t__set__random_anisotropy(this, f90wrap_random_anisotropy)
    use inputdatatype, only: ham_inp_t
    implicit none
    type ham_inp_t_ptr_type
        type(ham_inp_t), pointer :: p => NULL()
    end type ham_inp_t_ptr_type
    integer, intent(in)   :: this(2)
    type(ham_inp_t_ptr_type) :: this_ptr
    logical, intent(in) :: f90wrap_random_anisotropy
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%random_anisotropy = f90wrap_random_anisotropy
end subroutine f90wrap_ham_inp_t__set__random_anisotropy

subroutine f90wrap_ham_inp_t__get__mult_axis(this, f90wrap_mult_axis)
    use inputdatatype, only: ham_inp_t
    implicit none
    type ham_inp_t_ptr_type
        type(ham_inp_t), pointer :: p => NULL()
    end type ham_inp_t_ptr_type
    integer, intent(in)   :: this(2)
    type(ham_inp_t_ptr_type) :: this_ptr
    character(1), intent(out) :: f90wrap_mult_axis
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_mult_axis = this_ptr%p%mult_axis
end subroutine f90wrap_ham_inp_t__get__mult_axis

subroutine f90wrap_ham_inp_t__set__mult_axis(this, f90wrap_mult_axis)
    use inputdatatype, only: ham_inp_t
    implicit none
    type ham_inp_t_ptr_type
        type(ham_inp_t), pointer :: p => NULL()
    end type ham_inp_t_ptr_type
    integer, intent(in)   :: this(2)
    type(ham_inp_t_ptr_type) :: this_ptr
    character(1), intent(in) :: f90wrap_mult_axis
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%mult_axis = f90wrap_mult_axis
end subroutine f90wrap_ham_inp_t__set__mult_axis

subroutine f90wrap_ham_inp_t__get__kfile(this, f90wrap_kfile)
    use inputdatatype, only: ham_inp_t
    implicit none
    type ham_inp_t_ptr_type
        type(ham_inp_t), pointer :: p => NULL()
    end type ham_inp_t_ptr_type
    integer, intent(in)   :: this(2)
    type(ham_inp_t_ptr_type) :: this_ptr
    character(35), intent(out) :: f90wrap_kfile
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_kfile = this_ptr%p%kfile
end subroutine f90wrap_ham_inp_t__get__kfile

subroutine f90wrap_ham_inp_t__set__kfile(this, f90wrap_kfile)
    use inputdatatype, only: ham_inp_t
    implicit none
    type ham_inp_t_ptr_type
        type(ham_inp_t), pointer :: p => NULL()
    end type ham_inp_t_ptr_type
    integer, intent(in)   :: this(2)
    type(ham_inp_t_ptr_type) :: this_ptr
    character(35), intent(in) :: f90wrap_kfile
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%kfile = f90wrap_kfile
end subroutine f90wrap_ham_inp_t__set__kfile

subroutine f90wrap_ham_inp_t__array__anisotropytype(this, nd, dtype, dshape, dloc)
    use inputdatatype, only: ham_inp_t
    use, intrinsic :: iso_c_binding, only : c_int
    implicit none
    type ham_inp_t_ptr_type
        type(ham_inp_t), pointer :: p => NULL()
    end type ham_inp_t_ptr_type
    integer(c_int), intent(in) :: this(2)
    type(ham_inp_t_ptr_type) :: this_ptr
    integer(c_int), intent(out) :: nd
    integer(c_int), intent(out) :: dtype
    integer(c_int), dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 2
    dtype = 5
    this_ptr = transfer(this, this_ptr)
    if (allocated(this_ptr%p%anisotropytype)) then
        dshape(1:2) = shape(this_ptr%p%anisotropytype)
        dloc = loc(this_ptr%p%anisotropytype)
    else
        dloc = 0
    end if
end subroutine f90wrap_ham_inp_t__array__anisotropytype

subroutine f90wrap_ham_inp_t__array__anisotropytype_diff(this, nd, dtype, dshape, dloc)
    use inputdatatype, only: ham_inp_t
    use, intrinsic :: iso_c_binding, only : c_int
    implicit none
    type ham_inp_t_ptr_type
        type(ham_inp_t), pointer :: p => NULL()
    end type ham_inp_t_ptr_type
    integer(c_int), intent(in) :: this(2)
    type(ham_inp_t_ptr_type) :: this_ptr
    integer(c_int), intent(out) :: nd
    integer(c_int), intent(out) :: dtype
    integer(c_int), dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 2
    dtype = 5
    this_ptr = transfer(this, this_ptr)
    if (allocated(this_ptr%p%anisotropytype_diff)) then
        dshape(1:2) = shape(this_ptr%p%anisotropytype_diff)
        dloc = loc(this_ptr%p%anisotropytype_diff)
    else
        dloc = 0
    end if
end subroutine f90wrap_ham_inp_t__array__anisotropytype_diff

subroutine f90wrap_ham_inp_t__array__anisotropy(this, nd, dtype, dshape, dloc)
    use inputdatatype, only: ham_inp_t
    use, intrinsic :: iso_c_binding, only : c_int
    implicit none
    type ham_inp_t_ptr_type
        type(ham_inp_t), pointer :: p => NULL()
    end type ham_inp_t_ptr_type
    integer(c_int), intent(in) :: this(2)
    type(ham_inp_t_ptr_type) :: this_ptr
    integer(c_int), intent(out) :: nd
    integer(c_int), intent(out) :: dtype
    integer(c_int), dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 3
    dtype = 12
    this_ptr = transfer(this, this_ptr)
    if (allocated(this_ptr%p%anisotropy)) then
        dshape(1:3) = shape(this_ptr%p%anisotropy)
        dloc = loc(this_ptr%p%anisotropy)
    else
        dloc = 0
    end if
end subroutine f90wrap_ham_inp_t__array__anisotropy

subroutine f90wrap_ham_inp_t__array__anisotropy_diff(this, nd, dtype, dshape, dloc)
    use inputdatatype, only: ham_inp_t
    use, intrinsic :: iso_c_binding, only : c_int
    implicit none
    type ham_inp_t_ptr_type
        type(ham_inp_t), pointer :: p => NULL()
    end type ham_inp_t_ptr_type
    integer(c_int), intent(in) :: this(2)
    type(ham_inp_t_ptr_type) :: this_ptr
    integer(c_int), intent(out) :: nd
    integer(c_int), intent(out) :: dtype
    integer(c_int), dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 3
    dtype = 12
    this_ptr = transfer(this, this_ptr)
    if (allocated(this_ptr%p%anisotropy_diff)) then
        dshape(1:3) = shape(this_ptr%p%anisotropy_diff)
        dloc = loc(this_ptr%p%anisotropy_diff)
    else
        dloc = 0
    end if
end subroutine f90wrap_ham_inp_t__array__anisotropy_diff

subroutine f90wrap_ham_inp_t__get__do_dm(this, f90wrap_do_dm)
    use inputdatatype, only: ham_inp_t
    implicit none
    type ham_inp_t_ptr_type
        type(ham_inp_t), pointer :: p => NULL()
    end type ham_inp_t_ptr_type
    integer, intent(in)   :: this(2)
    type(ham_inp_t_ptr_type) :: this_ptr
    integer, intent(out) :: f90wrap_do_dm
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_do_dm = this_ptr%p%do_dm
end subroutine f90wrap_ham_inp_t__get__do_dm

subroutine f90wrap_ham_inp_t__set__do_dm(this, f90wrap_do_dm)
    use inputdatatype, only: ham_inp_t
    implicit none
    type ham_inp_t_ptr_type
        type(ham_inp_t), pointer :: p => NULL()
    end type ham_inp_t_ptr_type
    integer, intent(in)   :: this(2)
    type(ham_inp_t_ptr_type) :: this_ptr
    integer, intent(in) :: f90wrap_do_dm
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%do_dm = f90wrap_do_dm
end subroutine f90wrap_ham_inp_t__set__do_dm

subroutine f90wrap_ham_inp_t__get__max_no_dmshells(this, f90wrap_max_no_dmshells)
    use inputdatatype, only: ham_inp_t
    implicit none
    type ham_inp_t_ptr_type
        type(ham_inp_t), pointer :: p => NULL()
    end type ham_inp_t_ptr_type
    integer, intent(in)   :: this(2)
    type(ham_inp_t_ptr_type) :: this_ptr
    integer, intent(out) :: f90wrap_max_no_dmshells
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_max_no_dmshells = this_ptr%p%max_no_dmshells
end subroutine f90wrap_ham_inp_t__get__max_no_dmshells

subroutine f90wrap_ham_inp_t__set__max_no_dmshells(this, f90wrap_max_no_dmshells)
    use inputdatatype, only: ham_inp_t
    implicit none
    type ham_inp_t_ptr_type
        type(ham_inp_t), pointer :: p => NULL()
    end type ham_inp_t_ptr_type
    integer, intent(in)   :: this(2)
    type(ham_inp_t_ptr_type) :: this_ptr
    integer, intent(in) :: f90wrap_max_no_dmshells
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%max_no_dmshells = f90wrap_max_no_dmshells
end subroutine f90wrap_ham_inp_t__set__max_no_dmshells

subroutine f90wrap_ham_inp_t__get__dmfile(this, f90wrap_dmfile)
    use inputdatatype, only: ham_inp_t
    implicit none
    type ham_inp_t_ptr_type
        type(ham_inp_t), pointer :: p => NULL()
    end type ham_inp_t_ptr_type
    integer, intent(in)   :: this(2)
    type(ham_inp_t_ptr_type) :: this_ptr
    character(35), intent(out) :: f90wrap_dmfile
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_dmfile = this_ptr%p%dmfile
end subroutine f90wrap_ham_inp_t__get__dmfile

subroutine f90wrap_ham_inp_t__set__dmfile(this, f90wrap_dmfile)
    use inputdatatype, only: ham_inp_t
    implicit none
    type ham_inp_t_ptr_type
        type(ham_inp_t), pointer :: p => NULL()
    end type ham_inp_t_ptr_type
    integer, intent(in)   :: this(2)
    type(ham_inp_t_ptr_type) :: this_ptr
    character(35), intent(in) :: f90wrap_dmfile
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%dmfile = f90wrap_dmfile
end subroutine f90wrap_ham_inp_t__set__dmfile

subroutine f90wrap_ham_inp_t__get__dm_scale(this, f90wrap_dm_scale)
    use inputdatatype, only: ham_inp_t
    implicit none
    type ham_inp_t_ptr_type
        type(ham_inp_t), pointer :: p => NULL()
    end type ham_inp_t_ptr_type
    integer, intent(in)   :: this(2)
    type(ham_inp_t_ptr_type) :: this_ptr
    real(8), intent(out) :: f90wrap_dm_scale
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_dm_scale = this_ptr%p%dm_scale
end subroutine f90wrap_ham_inp_t__get__dm_scale

subroutine f90wrap_ham_inp_t__set__dm_scale(this, f90wrap_dm_scale)
    use inputdatatype, only: ham_inp_t
    implicit none
    type ham_inp_t_ptr_type
        type(ham_inp_t), pointer :: p => NULL()
    end type ham_inp_t_ptr_type
    integer, intent(in)   :: this(2)
    type(ham_inp_t_ptr_type) :: this_ptr
    real(8), intent(in) :: f90wrap_dm_scale
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%dm_scale = f90wrap_dm_scale
end subroutine f90wrap_ham_inp_t__set__dm_scale

subroutine f90wrap_ham_inp_t__array__dm_nn(this, nd, dtype, dshape, dloc)
    use inputdatatype, only: ham_inp_t
    use, intrinsic :: iso_c_binding, only : c_int
    implicit none
    type ham_inp_t_ptr_type
        type(ham_inp_t), pointer :: p => NULL()
    end type ham_inp_t_ptr_type
    integer(c_int), intent(in) :: this(2)
    type(ham_inp_t_ptr_type) :: this_ptr
    integer(c_int), intent(out) :: nd
    integer(c_int), intent(out) :: dtype
    integer(c_int), dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 1
    dtype = 5
    this_ptr = transfer(this, this_ptr)
    if (allocated(this_ptr%p%dm_nn)) then
        dshape(1:1) = shape(this_ptr%p%dm_nn)
        dloc = loc(this_ptr%p%dm_nn)
    else
        dloc = 0
    end if
end subroutine f90wrap_ham_inp_t__array__dm_nn

subroutine f90wrap_ham_inp_t__array__dm_redcoord(this, nd, dtype, dshape, dloc)
    use inputdatatype, only: ham_inp_t
    use, intrinsic :: iso_c_binding, only : c_int
    implicit none
    type ham_inp_t_ptr_type
        type(ham_inp_t), pointer :: p => NULL()
    end type ham_inp_t_ptr_type
    integer(c_int), intent(in) :: this(2)
    type(ham_inp_t_ptr_type) :: this_ptr
    integer(c_int), intent(out) :: nd
    integer(c_int), intent(out) :: dtype
    integer(c_int), dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 3
    dtype = 12
    this_ptr = transfer(this, this_ptr)
    if (allocated(this_ptr%p%dm_redcoord)) then
        dshape(1:3) = shape(this_ptr%p%dm_redcoord)
        dloc = loc(this_ptr%p%dm_redcoord)
    else
        dloc = 0
    end if
end subroutine f90wrap_ham_inp_t__array__dm_redcoord

subroutine f90wrap_ham_inp_t__array__dm_inpvect(this, nd, dtype, dshape, dloc)
    use inputdatatype, only: ham_inp_t
    use, intrinsic :: iso_c_binding, only : c_int
    implicit none
    type ham_inp_t_ptr_type
        type(ham_inp_t), pointer :: p => NULL()
    end type ham_inp_t_ptr_type
    integer(c_int), intent(in) :: this(2)
    type(ham_inp_t_ptr_type) :: this_ptr
    integer(c_int), intent(out) :: nd
    integer(c_int), intent(out) :: dtype
    integer(c_int), dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 5
    dtype = 12
    this_ptr = transfer(this, this_ptr)
    if (allocated(this_ptr%p%dm_inpvect)) then
        dshape(1:5) = shape(this_ptr%p%dm_inpvect)
        dloc = loc(this_ptr%p%dm_inpvect)
    else
        dloc = 0
    end if
end subroutine f90wrap_ham_inp_t__array__dm_inpvect

subroutine f90wrap_ham_inp_t__get__do_sa(this, f90wrap_do_sa)
    use inputdatatype, only: ham_inp_t
    implicit none
    type ham_inp_t_ptr_type
        type(ham_inp_t), pointer :: p => NULL()
    end type ham_inp_t_ptr_type
    integer, intent(in)   :: this(2)
    type(ham_inp_t_ptr_type) :: this_ptr
    integer, intent(out) :: f90wrap_do_sa
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_do_sa = this_ptr%p%do_sa
end subroutine f90wrap_ham_inp_t__get__do_sa

subroutine f90wrap_ham_inp_t__set__do_sa(this, f90wrap_do_sa)
    use inputdatatype, only: ham_inp_t
    implicit none
    type ham_inp_t_ptr_type
        type(ham_inp_t), pointer :: p => NULL()
    end type ham_inp_t_ptr_type
    integer, intent(in)   :: this(2)
    type(ham_inp_t_ptr_type) :: this_ptr
    integer, intent(in) :: f90wrap_do_sa
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%do_sa = f90wrap_do_sa
end subroutine f90wrap_ham_inp_t__set__do_sa

subroutine f90wrap_ham_inp_t__get__max_no_sashells(this, f90wrap_max_no_sashells)
    use inputdatatype, only: ham_inp_t
    implicit none
    type ham_inp_t_ptr_type
        type(ham_inp_t), pointer :: p => NULL()
    end type ham_inp_t_ptr_type
    integer, intent(in)   :: this(2)
    type(ham_inp_t_ptr_type) :: this_ptr
    integer, intent(out) :: f90wrap_max_no_sashells
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_max_no_sashells = this_ptr%p%max_no_sashells
end subroutine f90wrap_ham_inp_t__get__max_no_sashells

subroutine f90wrap_ham_inp_t__set__max_no_sashells(this, f90wrap_max_no_sashells)
    use inputdatatype, only: ham_inp_t
    implicit none
    type ham_inp_t_ptr_type
        type(ham_inp_t), pointer :: p => NULL()
    end type ham_inp_t_ptr_type
    integer, intent(in)   :: this(2)
    type(ham_inp_t_ptr_type) :: this_ptr
    integer, intent(in) :: f90wrap_max_no_sashells
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%max_no_sashells = f90wrap_max_no_sashells
end subroutine f90wrap_ham_inp_t__set__max_no_sashells

subroutine f90wrap_ham_inp_t__get__safile(this, f90wrap_safile)
    use inputdatatype, only: ham_inp_t
    implicit none
    type ham_inp_t_ptr_type
        type(ham_inp_t), pointer :: p => NULL()
    end type ham_inp_t_ptr_type
    integer, intent(in)   :: this(2)
    type(ham_inp_t_ptr_type) :: this_ptr
    character(35), intent(out) :: f90wrap_safile
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_safile = this_ptr%p%safile
end subroutine f90wrap_ham_inp_t__get__safile

subroutine f90wrap_ham_inp_t__set__safile(this, f90wrap_safile)
    use inputdatatype, only: ham_inp_t
    implicit none
    type ham_inp_t_ptr_type
        type(ham_inp_t), pointer :: p => NULL()
    end type ham_inp_t_ptr_type
    integer, intent(in)   :: this(2)
    type(ham_inp_t_ptr_type) :: this_ptr
    character(35), intent(in) :: f90wrap_safile
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%safile = f90wrap_safile
end subroutine f90wrap_ham_inp_t__set__safile

subroutine f90wrap_ham_inp_t__get__sa_scale(this, f90wrap_sa_scale)
    use inputdatatype, only: ham_inp_t
    implicit none
    type ham_inp_t_ptr_type
        type(ham_inp_t), pointer :: p => NULL()
    end type ham_inp_t_ptr_type
    integer, intent(in)   :: this(2)
    type(ham_inp_t_ptr_type) :: this_ptr
    real(8), intent(out) :: f90wrap_sa_scale
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_sa_scale = this_ptr%p%sa_scale
end subroutine f90wrap_ham_inp_t__get__sa_scale

subroutine f90wrap_ham_inp_t__set__sa_scale(this, f90wrap_sa_scale)
    use inputdatatype, only: ham_inp_t
    implicit none
    type ham_inp_t_ptr_type
        type(ham_inp_t), pointer :: p => NULL()
    end type ham_inp_t_ptr_type
    integer, intent(in)   :: this(2)
    type(ham_inp_t_ptr_type) :: this_ptr
    real(8), intent(in) :: f90wrap_sa_scale
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%sa_scale = f90wrap_sa_scale
end subroutine f90wrap_ham_inp_t__set__sa_scale

subroutine f90wrap_ham_inp_t__array__sa_nn(this, nd, dtype, dshape, dloc)
    use inputdatatype, only: ham_inp_t
    use, intrinsic :: iso_c_binding, only : c_int
    implicit none
    type ham_inp_t_ptr_type
        type(ham_inp_t), pointer :: p => NULL()
    end type ham_inp_t_ptr_type
    integer(c_int), intent(in) :: this(2)
    type(ham_inp_t_ptr_type) :: this_ptr
    integer(c_int), intent(out) :: nd
    integer(c_int), intent(out) :: dtype
    integer(c_int), dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 1
    dtype = 5
    this_ptr = transfer(this, this_ptr)
    if (allocated(this_ptr%p%sa_nn)) then
        dshape(1:1) = shape(this_ptr%p%sa_nn)
        dloc = loc(this_ptr%p%sa_nn)
    else
        dloc = 0
    end if
end subroutine f90wrap_ham_inp_t__array__sa_nn

subroutine f90wrap_ham_inp_t__array__sa_redcoord(this, nd, dtype, dshape, dloc)
    use inputdatatype, only: ham_inp_t
    use, intrinsic :: iso_c_binding, only : c_int
    implicit none
    type ham_inp_t_ptr_type
        type(ham_inp_t), pointer :: p => NULL()
    end type ham_inp_t_ptr_type
    integer(c_int), intent(in) :: this(2)
    type(ham_inp_t_ptr_type) :: this_ptr
    integer(c_int), intent(out) :: nd
    integer(c_int), intent(out) :: dtype
    integer(c_int), dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 3
    dtype = 12
    this_ptr = transfer(this, this_ptr)
    if (allocated(this_ptr%p%sa_redcoord)) then
        dshape(1:3) = shape(this_ptr%p%sa_redcoord)
        dloc = loc(this_ptr%p%sa_redcoord)
    else
        dloc = 0
    end if
end subroutine f90wrap_ham_inp_t__array__sa_redcoord

subroutine f90wrap_ham_inp_t__array__sa_inpvect(this, nd, dtype, dshape, dloc)
    use inputdatatype, only: ham_inp_t
    use, intrinsic :: iso_c_binding, only : c_int
    implicit none
    type ham_inp_t_ptr_type
        type(ham_inp_t), pointer :: p => NULL()
    end type ham_inp_t_ptr_type
    integer(c_int), intent(in) :: this(2)
    type(ham_inp_t_ptr_type) :: this_ptr
    integer(c_int), intent(out) :: nd
    integer(c_int), intent(out) :: dtype
    integer(c_int), dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 5
    dtype = 12
    this_ptr = transfer(this, this_ptr)
    if (allocated(this_ptr%p%sa_inpvect)) then
        dshape(1:5) = shape(this_ptr%p%sa_inpvect)
        dloc = loc(this_ptr%p%sa_inpvect)
    else
        dloc = 0
    end if
end subroutine f90wrap_ham_inp_t__array__sa_inpvect

subroutine f90wrap_ham_inp_t__get__do_chir(this, f90wrap_do_chir)
    use inputdatatype, only: ham_inp_t
    implicit none
    type ham_inp_t_ptr_type
        type(ham_inp_t), pointer :: p => NULL()
    end type ham_inp_t_ptr_type
    integer, intent(in)   :: this(2)
    type(ham_inp_t_ptr_type) :: this_ptr
    integer, intent(out) :: f90wrap_do_chir
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_do_chir = this_ptr%p%do_chir
end subroutine f90wrap_ham_inp_t__get__do_chir

subroutine f90wrap_ham_inp_t__set__do_chir(this, f90wrap_do_chir)
    use inputdatatype, only: ham_inp_t
    implicit none
    type ham_inp_t_ptr_type
        type(ham_inp_t), pointer :: p => NULL()
    end type ham_inp_t_ptr_type
    integer, intent(in)   :: this(2)
    type(ham_inp_t_ptr_type) :: this_ptr
    integer, intent(in) :: f90wrap_do_chir
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%do_chir = f90wrap_do_chir
end subroutine f90wrap_ham_inp_t__set__do_chir

subroutine f90wrap_ham_inp_t__get__max_no_chirshells(this, f90wrap_max_no_chirshells)
    use inputdatatype, only: ham_inp_t
    implicit none
    type ham_inp_t_ptr_type
        type(ham_inp_t), pointer :: p => NULL()
    end type ham_inp_t_ptr_type
    integer, intent(in)   :: this(2)
    type(ham_inp_t_ptr_type) :: this_ptr
    integer, intent(out) :: f90wrap_max_no_chirshells
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_max_no_chirshells = this_ptr%p%max_no_chirshells
end subroutine f90wrap_ham_inp_t__get__max_no_chirshells

subroutine f90wrap_ham_inp_t__set__max_no_chirshells(this, f90wrap_max_no_chirshells)
    use inputdatatype, only: ham_inp_t
    implicit none
    type ham_inp_t_ptr_type
        type(ham_inp_t), pointer :: p => NULL()
    end type ham_inp_t_ptr_type
    integer, intent(in)   :: this(2)
    type(ham_inp_t_ptr_type) :: this_ptr
    integer, intent(in) :: f90wrap_max_no_chirshells
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%max_no_chirshells = f90wrap_max_no_chirshells
end subroutine f90wrap_ham_inp_t__set__max_no_chirshells

subroutine f90wrap_ham_inp_t__get__chirfile(this, f90wrap_chirfile)
    use inputdatatype, only: ham_inp_t
    implicit none
    type ham_inp_t_ptr_type
        type(ham_inp_t), pointer :: p => NULL()
    end type ham_inp_t_ptr_type
    integer, intent(in)   :: this(2)
    type(ham_inp_t_ptr_type) :: this_ptr
    character(35), intent(out) :: f90wrap_chirfile
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_chirfile = this_ptr%p%chirfile
end subroutine f90wrap_ham_inp_t__get__chirfile

subroutine f90wrap_ham_inp_t__set__chirfile(this, f90wrap_chirfile)
    use inputdatatype, only: ham_inp_t
    implicit none
    type ham_inp_t_ptr_type
        type(ham_inp_t), pointer :: p => NULL()
    end type ham_inp_t_ptr_type
    integer, intent(in)   :: this(2)
    type(ham_inp_t_ptr_type) :: this_ptr
    character(35), intent(in) :: f90wrap_chirfile
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%chirfile = f90wrap_chirfile
end subroutine f90wrap_ham_inp_t__set__chirfile

subroutine f90wrap_ham_inp_t__array__chir_nn(this, nd, dtype, dshape, dloc)
    use inputdatatype, only: ham_inp_t
    use, intrinsic :: iso_c_binding, only : c_int
    implicit none
    type ham_inp_t_ptr_type
        type(ham_inp_t), pointer :: p => NULL()
    end type ham_inp_t_ptr_type
    integer(c_int), intent(in) :: this(2)
    type(ham_inp_t_ptr_type) :: this_ptr
    integer(c_int), intent(out) :: nd
    integer(c_int), intent(out) :: dtype
    integer(c_int), dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 1
    dtype = 5
    this_ptr = transfer(this, this_ptr)
    if (allocated(this_ptr%p%chir_nn)) then
        dshape(1:1) = shape(this_ptr%p%chir_nn)
        dloc = loc(this_ptr%p%chir_nn)
    else
        dloc = 0
    end if
end subroutine f90wrap_ham_inp_t__array__chir_nn

subroutine f90wrap_ham_inp_t__array__chir_inpval(this, nd, dtype, dshape, dloc)
    use inputdatatype, only: ham_inp_t
    use, intrinsic :: iso_c_binding, only : c_int
    implicit none
    type ham_inp_t_ptr_type
        type(ham_inp_t), pointer :: p => NULL()
    end type ham_inp_t_ptr_type
    integer(c_int), intent(in) :: this(2)
    type(ham_inp_t_ptr_type) :: this_ptr
    integer(c_int), intent(out) :: nd
    integer(c_int), intent(out) :: dtype
    integer(c_int), dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 4
    dtype = 12
    this_ptr = transfer(this, this_ptr)
    if (allocated(this_ptr%p%chir_inpval)) then
        dshape(1:4) = shape(this_ptr%p%chir_inpval)
        dloc = loc(this_ptr%p%chir_inpval)
    else
        dloc = 0
    end if
end subroutine f90wrap_ham_inp_t__array__chir_inpval

subroutine f90wrap_ham_inp_t__array__chir_redcoord(this, nd, dtype, dshape, dloc)
    use inputdatatype, only: ham_inp_t
    use, intrinsic :: iso_c_binding, only : c_int
    implicit none
    type ham_inp_t_ptr_type
        type(ham_inp_t), pointer :: p => NULL()
    end type ham_inp_t_ptr_type
    integer(c_int), intent(in) :: this(2)
    type(ham_inp_t_ptr_type) :: this_ptr
    integer(c_int), intent(out) :: nd
    integer(c_int), intent(out) :: dtype
    integer(c_int), dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 4
    dtype = 12
    this_ptr = transfer(this, this_ptr)
    if (allocated(this_ptr%p%chir_redcoord)) then
        dshape(1:4) = shape(this_ptr%p%chir_redcoord)
        dloc = loc(this_ptr%p%chir_redcoord)
    else
        dloc = 0
    end if
end subroutine f90wrap_ham_inp_t__array__chir_redcoord

subroutine f90wrap_ham_inp_t__get__do_fourx(this, f90wrap_do_fourx)
    use inputdatatype, only: ham_inp_t
    implicit none
    type ham_inp_t_ptr_type
        type(ham_inp_t), pointer :: p => NULL()
    end type ham_inp_t_ptr_type
    integer, intent(in)   :: this(2)
    type(ham_inp_t_ptr_type) :: this_ptr
    integer, intent(out) :: f90wrap_do_fourx
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_do_fourx = this_ptr%p%do_fourx
end subroutine f90wrap_ham_inp_t__get__do_fourx

subroutine f90wrap_ham_inp_t__set__do_fourx(this, f90wrap_do_fourx)
    use inputdatatype, only: ham_inp_t
    implicit none
    type ham_inp_t_ptr_type
        type(ham_inp_t), pointer :: p => NULL()
    end type ham_inp_t_ptr_type
    integer, intent(in)   :: this(2)
    type(ham_inp_t_ptr_type) :: this_ptr
    integer, intent(in) :: f90wrap_do_fourx
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%do_fourx = f90wrap_do_fourx
end subroutine f90wrap_ham_inp_t__set__do_fourx

subroutine f90wrap_ham_inp_t__get__max_no_fourxshells(this, f90wrap_max_no_fourxshells)
    use inputdatatype, only: ham_inp_t
    implicit none
    type ham_inp_t_ptr_type
        type(ham_inp_t), pointer :: p => NULL()
    end type ham_inp_t_ptr_type
    integer, intent(in)   :: this(2)
    type(ham_inp_t_ptr_type) :: this_ptr
    integer, intent(out) :: f90wrap_max_no_fourxshells
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_max_no_fourxshells = this_ptr%p%max_no_fourxshells
end subroutine f90wrap_ham_inp_t__get__max_no_fourxshells

subroutine f90wrap_ham_inp_t__set__max_no_fourxshells(this, f90wrap_max_no_fourxshells)
    use inputdatatype, only: ham_inp_t
    implicit none
    type ham_inp_t_ptr_type
        type(ham_inp_t), pointer :: p => NULL()
    end type ham_inp_t_ptr_type
    integer, intent(in)   :: this(2)
    type(ham_inp_t_ptr_type) :: this_ptr
    integer, intent(in) :: f90wrap_max_no_fourxshells
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%max_no_fourxshells = f90wrap_max_no_fourxshells
end subroutine f90wrap_ham_inp_t__set__max_no_fourxshells

subroutine f90wrap_ham_inp_t__get__fourxfile(this, f90wrap_fourxfile)
    use inputdatatype, only: ham_inp_t
    implicit none
    type ham_inp_t_ptr_type
        type(ham_inp_t), pointer :: p => NULL()
    end type ham_inp_t_ptr_type
    integer, intent(in)   :: this(2)
    type(ham_inp_t_ptr_type) :: this_ptr
    character(35), intent(out) :: f90wrap_fourxfile
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_fourxfile = this_ptr%p%fourxfile
end subroutine f90wrap_ham_inp_t__get__fourxfile

subroutine f90wrap_ham_inp_t__set__fourxfile(this, f90wrap_fourxfile)
    use inputdatatype, only: ham_inp_t
    implicit none
    type ham_inp_t_ptr_type
        type(ham_inp_t), pointer :: p => NULL()
    end type ham_inp_t_ptr_type
    integer, intent(in)   :: this(2)
    type(ham_inp_t_ptr_type) :: this_ptr
    character(35), intent(in) :: f90wrap_fourxfile
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%fourxfile = f90wrap_fourxfile
end subroutine f90wrap_ham_inp_t__set__fourxfile

subroutine f90wrap_ham_inp_t__array__fourx_nn(this, nd, dtype, dshape, dloc)
    use inputdatatype, only: ham_inp_t
    use, intrinsic :: iso_c_binding, only : c_int
    implicit none
    type ham_inp_t_ptr_type
        type(ham_inp_t), pointer :: p => NULL()
    end type ham_inp_t_ptr_type
    integer(c_int), intent(in) :: this(2)
    type(ham_inp_t_ptr_type) :: this_ptr
    integer(c_int), intent(out) :: nd
    integer(c_int), intent(out) :: dtype
    integer(c_int), dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 1
    dtype = 5
    this_ptr = transfer(this, this_ptr)
    if (allocated(this_ptr%p%fourx_nn)) then
        dshape(1:1) = shape(this_ptr%p%fourx_nn)
        dloc = loc(this_ptr%p%fourx_nn)
    else
        dloc = 0
    end if
end subroutine f90wrap_ham_inp_t__array__fourx_nn

subroutine f90wrap_ham_inp_t__array__fourx_inpval(this, nd, dtype, dshape, dloc)
    use inputdatatype, only: ham_inp_t
    use, intrinsic :: iso_c_binding, only : c_int
    implicit none
    type ham_inp_t_ptr_type
        type(ham_inp_t), pointer :: p => NULL()
    end type ham_inp_t_ptr_type
    integer(c_int), intent(in) :: this(2)
    type(ham_inp_t_ptr_type) :: this_ptr
    integer(c_int), intent(out) :: nd
    integer(c_int), intent(out) :: dtype
    integer(c_int), dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 4
    dtype = 12
    this_ptr = transfer(this, this_ptr)
    if (allocated(this_ptr%p%fourx_inpval)) then
        dshape(1:4) = shape(this_ptr%p%fourx_inpval)
        dloc = loc(this_ptr%p%fourx_inpval)
    else
        dloc = 0
    end if
end subroutine f90wrap_ham_inp_t__array__fourx_inpval

subroutine f90wrap_ham_inp_t__array__fourx_redcoord(this, nd, dtype, dshape, dloc)
    use inputdatatype, only: ham_inp_t
    use, intrinsic :: iso_c_binding, only : c_int
    implicit none
    type ham_inp_t_ptr_type
        type(ham_inp_t), pointer :: p => NULL()
    end type ham_inp_t_ptr_type
    integer(c_int), intent(in) :: this(2)
    type(ham_inp_t_ptr_type) :: this_ptr
    integer(c_int), intent(out) :: nd
    integer(c_int), intent(out) :: dtype
    integer(c_int), dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 4
    dtype = 12
    this_ptr = transfer(this, this_ptr)
    if (allocated(this_ptr%p%fourx_redcoord)) then
        dshape(1:4) = shape(this_ptr%p%fourx_redcoord)
        dloc = loc(this_ptr%p%fourx_redcoord)
    else
        dloc = 0
    end if
end subroutine f90wrap_ham_inp_t__array__fourx_redcoord

subroutine f90wrap_ham_inp_t__get__do_pd(this, f90wrap_do_pd)
    use inputdatatype, only: ham_inp_t
    implicit none
    type ham_inp_t_ptr_type
        type(ham_inp_t), pointer :: p => NULL()
    end type ham_inp_t_ptr_type
    integer, intent(in)   :: this(2)
    type(ham_inp_t_ptr_type) :: this_ptr
    integer, intent(out) :: f90wrap_do_pd
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_do_pd = this_ptr%p%do_pd
end subroutine f90wrap_ham_inp_t__get__do_pd

subroutine f90wrap_ham_inp_t__set__do_pd(this, f90wrap_do_pd)
    use inputdatatype, only: ham_inp_t
    implicit none
    type ham_inp_t_ptr_type
        type(ham_inp_t), pointer :: p => NULL()
    end type ham_inp_t_ptr_type
    integer, intent(in)   :: this(2)
    type(ham_inp_t_ptr_type) :: this_ptr
    integer, intent(in) :: f90wrap_do_pd
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%do_pd = f90wrap_do_pd
end subroutine f90wrap_ham_inp_t__set__do_pd

subroutine f90wrap_ham_inp_t__get__max_no_pdshells(this, f90wrap_max_no_pdshells)
    use inputdatatype, only: ham_inp_t
    implicit none
    type ham_inp_t_ptr_type
        type(ham_inp_t), pointer :: p => NULL()
    end type ham_inp_t_ptr_type
    integer, intent(in)   :: this(2)
    type(ham_inp_t_ptr_type) :: this_ptr
    integer, intent(out) :: f90wrap_max_no_pdshells
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_max_no_pdshells = this_ptr%p%max_no_pdshells
end subroutine f90wrap_ham_inp_t__get__max_no_pdshells

subroutine f90wrap_ham_inp_t__set__max_no_pdshells(this, f90wrap_max_no_pdshells)
    use inputdatatype, only: ham_inp_t
    implicit none
    type ham_inp_t_ptr_type
        type(ham_inp_t), pointer :: p => NULL()
    end type ham_inp_t_ptr_type
    integer, intent(in)   :: this(2)
    type(ham_inp_t_ptr_type) :: this_ptr
    integer, intent(in) :: f90wrap_max_no_pdshells
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%max_no_pdshells = f90wrap_max_no_pdshells
end subroutine f90wrap_ham_inp_t__set__max_no_pdshells

subroutine f90wrap_ham_inp_t__get__pdfile(this, f90wrap_pdfile)
    use inputdatatype, only: ham_inp_t
    implicit none
    type ham_inp_t_ptr_type
        type(ham_inp_t), pointer :: p => NULL()
    end type ham_inp_t_ptr_type
    integer, intent(in)   :: this(2)
    type(ham_inp_t_ptr_type) :: this_ptr
    character(35), intent(out) :: f90wrap_pdfile
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_pdfile = this_ptr%p%pdfile
end subroutine f90wrap_ham_inp_t__get__pdfile

subroutine f90wrap_ham_inp_t__set__pdfile(this, f90wrap_pdfile)
    use inputdatatype, only: ham_inp_t
    implicit none
    type ham_inp_t_ptr_type
        type(ham_inp_t), pointer :: p => NULL()
    end type ham_inp_t_ptr_type
    integer, intent(in)   :: this(2)
    type(ham_inp_t_ptr_type) :: this_ptr
    character(35), intent(in) :: f90wrap_pdfile
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%pdfile = f90wrap_pdfile
end subroutine f90wrap_ham_inp_t__set__pdfile

subroutine f90wrap_ham_inp_t__array__pd_nn(this, nd, dtype, dshape, dloc)
    use inputdatatype, only: ham_inp_t
    use, intrinsic :: iso_c_binding, only : c_int
    implicit none
    type ham_inp_t_ptr_type
        type(ham_inp_t), pointer :: p => NULL()
    end type ham_inp_t_ptr_type
    integer(c_int), intent(in) :: this(2)
    type(ham_inp_t_ptr_type) :: this_ptr
    integer(c_int), intent(out) :: nd
    integer(c_int), intent(out) :: dtype
    integer(c_int), dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 1
    dtype = 5
    this_ptr = transfer(this, this_ptr)
    if (allocated(this_ptr%p%pd_nn)) then
        dshape(1:1) = shape(this_ptr%p%pd_nn)
        dloc = loc(this_ptr%p%pd_nn)
    else
        dloc = 0
    end if
end subroutine f90wrap_ham_inp_t__array__pd_nn

subroutine f90wrap_ham_inp_t__array__pd_redcoord(this, nd, dtype, dshape, dloc)
    use inputdatatype, only: ham_inp_t
    use, intrinsic :: iso_c_binding, only : c_int
    implicit none
    type ham_inp_t_ptr_type
        type(ham_inp_t), pointer :: p => NULL()
    end type ham_inp_t_ptr_type
    integer(c_int), intent(in) :: this(2)
    type(ham_inp_t_ptr_type) :: this_ptr
    integer(c_int), intent(out) :: nd
    integer(c_int), intent(out) :: dtype
    integer(c_int), dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 3
    dtype = 12
    this_ptr = transfer(this, this_ptr)
    if (allocated(this_ptr%p%pd_redcoord)) then
        dshape(1:3) = shape(this_ptr%p%pd_redcoord)
        dloc = loc(this_ptr%p%pd_redcoord)
    else
        dloc = 0
    end if
end subroutine f90wrap_ham_inp_t__array__pd_redcoord

subroutine f90wrap_ham_inp_t__array__pd_inpvect(this, nd, dtype, dshape, dloc)
    use inputdatatype, only: ham_inp_t
    use, intrinsic :: iso_c_binding, only : c_int
    implicit none
    type ham_inp_t_ptr_type
        type(ham_inp_t), pointer :: p => NULL()
    end type ham_inp_t_ptr_type
    integer(c_int), intent(in) :: this(2)
    type(ham_inp_t_ptr_type) :: this_ptr
    integer(c_int), intent(out) :: nd
    integer(c_int), intent(out) :: dtype
    integer(c_int), dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 5
    dtype = 12
    this_ptr = transfer(this, this_ptr)
    if (allocated(this_ptr%p%pd_inpvect)) then
        dshape(1:5) = shape(this_ptr%p%pd_inpvect)
        dloc = loc(this_ptr%p%pd_inpvect)
    else
        dloc = 0
    end if
end subroutine f90wrap_ham_inp_t__array__pd_inpvect

subroutine f90wrap_ham_inp_t__get__do_biqdm(this, f90wrap_do_biqdm)
    use inputdatatype, only: ham_inp_t
    implicit none
    type ham_inp_t_ptr_type
        type(ham_inp_t), pointer :: p => NULL()
    end type ham_inp_t_ptr_type
    integer, intent(in)   :: this(2)
    type(ham_inp_t_ptr_type) :: this_ptr
    integer, intent(out) :: f90wrap_do_biqdm
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_do_biqdm = this_ptr%p%do_biqdm
end subroutine f90wrap_ham_inp_t__get__do_biqdm

subroutine f90wrap_ham_inp_t__set__do_biqdm(this, f90wrap_do_biqdm)
    use inputdatatype, only: ham_inp_t
    implicit none
    type ham_inp_t_ptr_type
        type(ham_inp_t), pointer :: p => NULL()
    end type ham_inp_t_ptr_type
    integer, intent(in)   :: this(2)
    type(ham_inp_t_ptr_type) :: this_ptr
    integer, intent(in) :: f90wrap_do_biqdm
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%do_biqdm = f90wrap_do_biqdm
end subroutine f90wrap_ham_inp_t__set__do_biqdm

subroutine f90wrap_ham_inp_t__get__max_no_biqdmshells(this, f90wrap_max_no_biqdmshells)
    use inputdatatype, only: ham_inp_t
    implicit none
    type ham_inp_t_ptr_type
        type(ham_inp_t), pointer :: p => NULL()
    end type ham_inp_t_ptr_type
    integer, intent(in)   :: this(2)
    type(ham_inp_t_ptr_type) :: this_ptr
    integer, intent(out) :: f90wrap_max_no_biqdmshells
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_max_no_biqdmshells = this_ptr%p%max_no_biqdmshells
end subroutine f90wrap_ham_inp_t__get__max_no_biqdmshells

subroutine f90wrap_ham_inp_t__set__max_no_biqdmshells(this, f90wrap_max_no_biqdmshells)
    use inputdatatype, only: ham_inp_t
    implicit none
    type ham_inp_t_ptr_type
        type(ham_inp_t), pointer :: p => NULL()
    end type ham_inp_t_ptr_type
    integer, intent(in)   :: this(2)
    type(ham_inp_t_ptr_type) :: this_ptr
    integer, intent(in) :: f90wrap_max_no_biqdmshells
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%max_no_biqdmshells = f90wrap_max_no_biqdmshells
end subroutine f90wrap_ham_inp_t__set__max_no_biqdmshells

subroutine f90wrap_ham_inp_t__get__biqdmfile(this, f90wrap_biqdmfile)
    use inputdatatype, only: ham_inp_t
    implicit none
    type ham_inp_t_ptr_type
        type(ham_inp_t), pointer :: p => NULL()
    end type ham_inp_t_ptr_type
    integer, intent(in)   :: this(2)
    type(ham_inp_t_ptr_type) :: this_ptr
    character(35), intent(out) :: f90wrap_biqdmfile
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_biqdmfile = this_ptr%p%biqdmfile
end subroutine f90wrap_ham_inp_t__get__biqdmfile

subroutine f90wrap_ham_inp_t__set__biqdmfile(this, f90wrap_biqdmfile)
    use inputdatatype, only: ham_inp_t
    implicit none
    type ham_inp_t_ptr_type
        type(ham_inp_t), pointer :: p => NULL()
    end type ham_inp_t_ptr_type
    integer, intent(in)   :: this(2)
    type(ham_inp_t_ptr_type) :: this_ptr
    character(35), intent(in) :: f90wrap_biqdmfile
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%biqdmfile = f90wrap_biqdmfile
end subroutine f90wrap_ham_inp_t__set__biqdmfile

subroutine f90wrap_ham_inp_t__array__biqdm_nn(this, nd, dtype, dshape, dloc)
    use inputdatatype, only: ham_inp_t
    use, intrinsic :: iso_c_binding, only : c_int
    implicit none
    type ham_inp_t_ptr_type
        type(ham_inp_t), pointer :: p => NULL()
    end type ham_inp_t_ptr_type
    integer(c_int), intent(in) :: this(2)
    type(ham_inp_t_ptr_type) :: this_ptr
    integer(c_int), intent(out) :: nd
    integer(c_int), intent(out) :: dtype
    integer(c_int), dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 1
    dtype = 5
    this_ptr = transfer(this, this_ptr)
    if (allocated(this_ptr%p%biqdm_nn)) then
        dshape(1:1) = shape(this_ptr%p%biqdm_nn)
        dloc = loc(this_ptr%p%biqdm_nn)
    else
        dloc = 0
    end if
end subroutine f90wrap_ham_inp_t__array__biqdm_nn

subroutine f90wrap_ham_inp_t__array__biqdm_redcoord(this, nd, dtype, dshape, dloc)
    use inputdatatype, only: ham_inp_t
    use, intrinsic :: iso_c_binding, only : c_int
    implicit none
    type ham_inp_t_ptr_type
        type(ham_inp_t), pointer :: p => NULL()
    end type ham_inp_t_ptr_type
    integer(c_int), intent(in) :: this(2)
    type(ham_inp_t_ptr_type) :: this_ptr
    integer(c_int), intent(out) :: nd
    integer(c_int), intent(out) :: dtype
    integer(c_int), dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 3
    dtype = 12
    this_ptr = transfer(this, this_ptr)
    if (allocated(this_ptr%p%biqdm_redcoord)) then
        dshape(1:3) = shape(this_ptr%p%biqdm_redcoord)
        dloc = loc(this_ptr%p%biqdm_redcoord)
    else
        dloc = 0
    end if
end subroutine f90wrap_ham_inp_t__array__biqdm_redcoord

subroutine f90wrap_ham_inp_t__array__biqdm_inpvect(this, nd, dtype, dshape, dloc)
    use inputdatatype, only: ham_inp_t
    use, intrinsic :: iso_c_binding, only : c_int
    implicit none
    type ham_inp_t_ptr_type
        type(ham_inp_t), pointer :: p => NULL()
    end type ham_inp_t_ptr_type
    integer(c_int), intent(in) :: this(2)
    type(ham_inp_t_ptr_type) :: this_ptr
    integer(c_int), intent(out) :: nd
    integer(c_int), intent(out) :: dtype
    integer(c_int), dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 5
    dtype = 12
    this_ptr = transfer(this, this_ptr)
    if (allocated(this_ptr%p%biqdm_inpvect)) then
        dshape(1:5) = shape(this_ptr%p%biqdm_inpvect)
        dloc = loc(this_ptr%p%biqdm_inpvect)
    else
        dloc = 0
    end if
end subroutine f90wrap_ham_inp_t__array__biqdm_inpvect

subroutine f90wrap_ham_inp_t__get__do_bq(this, f90wrap_do_bq)
    use inputdatatype, only: ham_inp_t
    implicit none
    type ham_inp_t_ptr_type
        type(ham_inp_t), pointer :: p => NULL()
    end type ham_inp_t_ptr_type
    integer, intent(in)   :: this(2)
    type(ham_inp_t_ptr_type) :: this_ptr
    integer, intent(out) :: f90wrap_do_bq
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_do_bq = this_ptr%p%do_bq
end subroutine f90wrap_ham_inp_t__get__do_bq

subroutine f90wrap_ham_inp_t__set__do_bq(this, f90wrap_do_bq)
    use inputdatatype, only: ham_inp_t
    implicit none
    type ham_inp_t_ptr_type
        type(ham_inp_t), pointer :: p => NULL()
    end type ham_inp_t_ptr_type
    integer, intent(in)   :: this(2)
    type(ham_inp_t_ptr_type) :: this_ptr
    integer, intent(in) :: f90wrap_do_bq
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%do_bq = f90wrap_do_bq
end subroutine f90wrap_ham_inp_t__set__do_bq

subroutine f90wrap_ham_inp_t__get__max_no_bqshells(this, f90wrap_max_no_bqshells)
    use inputdatatype, only: ham_inp_t
    implicit none
    type ham_inp_t_ptr_type
        type(ham_inp_t), pointer :: p => NULL()
    end type ham_inp_t_ptr_type
    integer, intent(in)   :: this(2)
    type(ham_inp_t_ptr_type) :: this_ptr
    integer, intent(out) :: f90wrap_max_no_bqshells
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_max_no_bqshells = this_ptr%p%max_no_bqshells
end subroutine f90wrap_ham_inp_t__get__max_no_bqshells

subroutine f90wrap_ham_inp_t__set__max_no_bqshells(this, f90wrap_max_no_bqshells)
    use inputdatatype, only: ham_inp_t
    implicit none
    type ham_inp_t_ptr_type
        type(ham_inp_t), pointer :: p => NULL()
    end type ham_inp_t_ptr_type
    integer, intent(in)   :: this(2)
    type(ham_inp_t_ptr_type) :: this_ptr
    integer, intent(in) :: f90wrap_max_no_bqshells
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%max_no_bqshells = f90wrap_max_no_bqshells
end subroutine f90wrap_ham_inp_t__set__max_no_bqshells

subroutine f90wrap_ham_inp_t__get__bqfile(this, f90wrap_bqfile)
    use inputdatatype, only: ham_inp_t
    implicit none
    type ham_inp_t_ptr_type
        type(ham_inp_t), pointer :: p => NULL()
    end type ham_inp_t_ptr_type
    integer, intent(in)   :: this(2)
    type(ham_inp_t_ptr_type) :: this_ptr
    character(35), intent(out) :: f90wrap_bqfile
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_bqfile = this_ptr%p%bqfile
end subroutine f90wrap_ham_inp_t__get__bqfile

subroutine f90wrap_ham_inp_t__set__bqfile(this, f90wrap_bqfile)
    use inputdatatype, only: ham_inp_t
    implicit none
    type ham_inp_t_ptr_type
        type(ham_inp_t), pointer :: p => NULL()
    end type ham_inp_t_ptr_type
    integer, intent(in)   :: this(2)
    type(ham_inp_t_ptr_type) :: this_ptr
    character(35), intent(in) :: f90wrap_bqfile
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%bqfile = f90wrap_bqfile
end subroutine f90wrap_ham_inp_t__set__bqfile

subroutine f90wrap_ham_inp_t__array__bq_nn(this, nd, dtype, dshape, dloc)
    use inputdatatype, only: ham_inp_t
    use, intrinsic :: iso_c_binding, only : c_int
    implicit none
    type ham_inp_t_ptr_type
        type(ham_inp_t), pointer :: p => NULL()
    end type ham_inp_t_ptr_type
    integer(c_int), intent(in) :: this(2)
    type(ham_inp_t_ptr_type) :: this_ptr
    integer(c_int), intent(out) :: nd
    integer(c_int), intent(out) :: dtype
    integer(c_int), dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 1
    dtype = 5
    this_ptr = transfer(this, this_ptr)
    if (allocated(this_ptr%p%bq_nn)) then
        dshape(1:1) = shape(this_ptr%p%bq_nn)
        dloc = loc(this_ptr%p%bq_nn)
    else
        dloc = 0
    end if
end subroutine f90wrap_ham_inp_t__array__bq_nn

subroutine f90wrap_ham_inp_t__array__jc_bq(this, nd, dtype, dshape, dloc)
    use inputdatatype, only: ham_inp_t
    use, intrinsic :: iso_c_binding, only : c_int
    implicit none
    type ham_inp_t_ptr_type
        type(ham_inp_t), pointer :: p => NULL()
    end type ham_inp_t_ptr_type
    integer(c_int), intent(in) :: this(2)
    type(ham_inp_t_ptr_type) :: this_ptr
    integer(c_int), intent(out) :: nd
    integer(c_int), intent(out) :: dtype
    integer(c_int), dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 4
    dtype = 12
    this_ptr = transfer(this, this_ptr)
    if (allocated(this_ptr%p%jc_bq)) then
        dshape(1:4) = shape(this_ptr%p%jc_bq)
        dloc = loc(this_ptr%p%jc_bq)
    else
        dloc = 0
    end if
end subroutine f90wrap_ham_inp_t__array__jc_bq

subroutine f90wrap_ham_inp_t__array__bq_redcoord(this, nd, dtype, dshape, dloc)
    use inputdatatype, only: ham_inp_t
    use, intrinsic :: iso_c_binding, only : c_int
    implicit none
    type ham_inp_t_ptr_type
        type(ham_inp_t), pointer :: p => NULL()
    end type ham_inp_t_ptr_type
    integer(c_int), intent(in) :: this(2)
    type(ham_inp_t_ptr_type) :: this_ptr
    integer(c_int), intent(out) :: nd
    integer(c_int), intent(out) :: dtype
    integer(c_int), dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 3
    dtype = 12
    this_ptr = transfer(this, this_ptr)
    if (allocated(this_ptr%p%bq_redcoord)) then
        dshape(1:3) = shape(this_ptr%p%bq_redcoord)
        dloc = loc(this_ptr%p%bq_redcoord)
    else
        dloc = 0
    end if
end subroutine f90wrap_ham_inp_t__array__bq_redcoord

subroutine f90wrap_ham_inp_t__get__do_ring(this, f90wrap_do_ring)
    use inputdatatype, only: ham_inp_t
    implicit none
    type ham_inp_t_ptr_type
        type(ham_inp_t), pointer :: p => NULL()
    end type ham_inp_t_ptr_type
    integer, intent(in)   :: this(2)
    type(ham_inp_t_ptr_type) :: this_ptr
    integer, intent(out) :: f90wrap_do_ring
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_do_ring = this_ptr%p%do_ring
end subroutine f90wrap_ham_inp_t__get__do_ring

subroutine f90wrap_ham_inp_t__set__do_ring(this, f90wrap_do_ring)
    use inputdatatype, only: ham_inp_t
    implicit none
    type ham_inp_t_ptr_type
        type(ham_inp_t), pointer :: p => NULL()
    end type ham_inp_t_ptr_type
    integer, intent(in)   :: this(2)
    type(ham_inp_t_ptr_type) :: this_ptr
    integer, intent(in) :: f90wrap_do_ring
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%do_ring = f90wrap_do_ring
end subroutine f90wrap_ham_inp_t__set__do_ring

subroutine f90wrap_ham_inp_t__get__max_no_ringshells(this, f90wrap_max_no_ringshells)
    use inputdatatype, only: ham_inp_t
    implicit none
    type ham_inp_t_ptr_type
        type(ham_inp_t), pointer :: p => NULL()
    end type ham_inp_t_ptr_type
    integer, intent(in)   :: this(2)
    type(ham_inp_t_ptr_type) :: this_ptr
    integer, intent(out) :: f90wrap_max_no_ringshells
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_max_no_ringshells = this_ptr%p%max_no_ringshells
end subroutine f90wrap_ham_inp_t__get__max_no_ringshells

subroutine f90wrap_ham_inp_t__set__max_no_ringshells(this, f90wrap_max_no_ringshells)
    use inputdatatype, only: ham_inp_t
    implicit none
    type ham_inp_t_ptr_type
        type(ham_inp_t), pointer :: p => NULL()
    end type ham_inp_t_ptr_type
    integer, intent(in)   :: this(2)
    type(ham_inp_t_ptr_type) :: this_ptr
    integer, intent(in) :: f90wrap_max_no_ringshells
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%max_no_ringshells = f90wrap_max_no_ringshells
end subroutine f90wrap_ham_inp_t__set__max_no_ringshells

subroutine f90wrap_ham_inp_t__get__ringfile(this, f90wrap_ringfile)
    use inputdatatype, only: ham_inp_t
    implicit none
    type ham_inp_t_ptr_type
        type(ham_inp_t), pointer :: p => NULL()
    end type ham_inp_t_ptr_type
    integer, intent(in)   :: this(2)
    type(ham_inp_t_ptr_type) :: this_ptr
    character(35), intent(out) :: f90wrap_ringfile
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_ringfile = this_ptr%p%ringfile
end subroutine f90wrap_ham_inp_t__get__ringfile

subroutine f90wrap_ham_inp_t__set__ringfile(this, f90wrap_ringfile)
    use inputdatatype, only: ham_inp_t
    implicit none
    type ham_inp_t_ptr_type
        type(ham_inp_t), pointer :: p => NULL()
    end type ham_inp_t_ptr_type
    integer, intent(in)   :: this(2)
    type(ham_inp_t_ptr_type) :: this_ptr
    character(35), intent(in) :: f90wrap_ringfile
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%ringfile = f90wrap_ringfile
end subroutine f90wrap_ham_inp_t__set__ringfile

subroutine f90wrap_ham_inp_t__array__ring_nn(this, nd, dtype, dshape, dloc)
    use inputdatatype, only: ham_inp_t
    use, intrinsic :: iso_c_binding, only : c_int
    implicit none
    type ham_inp_t_ptr_type
        type(ham_inp_t), pointer :: p => NULL()
    end type ham_inp_t_ptr_type
    integer(c_int), intent(in) :: this(2)
    type(ham_inp_t_ptr_type) :: this_ptr
    integer(c_int), intent(out) :: nd
    integer(c_int), intent(out) :: dtype
    integer(c_int), dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 1
    dtype = 5
    this_ptr = transfer(this, this_ptr)
    if (allocated(this_ptr%p%ring_nn)) then
        dshape(1:1) = shape(this_ptr%p%ring_nn)
        dloc = loc(this_ptr%p%ring_nn)
    else
        dloc = 0
    end if
end subroutine f90wrap_ham_inp_t__array__ring_nn

subroutine f90wrap_ham_inp_t__array__ring_redcoord_ij(this, nd, dtype, dshape, dloc)
    use inputdatatype, only: ham_inp_t
    use, intrinsic :: iso_c_binding, only : c_int
    implicit none
    type ham_inp_t_ptr_type
        type(ham_inp_t), pointer :: p => NULL()
    end type ham_inp_t_ptr_type
    integer(c_int), intent(in) :: this(2)
    type(ham_inp_t_ptr_type) :: this_ptr
    integer(c_int), intent(out) :: nd
    integer(c_int), intent(out) :: dtype
    integer(c_int), dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 3
    dtype = 12
    this_ptr = transfer(this, this_ptr)
    if (allocated(this_ptr%p%ring_redcoord_ij)) then
        dshape(1:3) = shape(this_ptr%p%ring_redcoord_ij)
        dloc = loc(this_ptr%p%ring_redcoord_ij)
    else
        dloc = 0
    end if
end subroutine f90wrap_ham_inp_t__array__ring_redcoord_ij

subroutine f90wrap_ham_inp_t__array__ring_redcoord_ik(this, nd, dtype, dshape, dloc)
    use inputdatatype, only: ham_inp_t
    use, intrinsic :: iso_c_binding, only : c_int
    implicit none
    type ham_inp_t_ptr_type
        type(ham_inp_t), pointer :: p => NULL()
    end type ham_inp_t_ptr_type
    integer(c_int), intent(in) :: this(2)
    type(ham_inp_t_ptr_type) :: this_ptr
    integer(c_int), intent(out) :: nd
    integer(c_int), intent(out) :: dtype
    integer(c_int), dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 3
    dtype = 12
    this_ptr = transfer(this, this_ptr)
    if (allocated(this_ptr%p%ring_redcoord_ik)) then
        dshape(1:3) = shape(this_ptr%p%ring_redcoord_ik)
        dloc = loc(this_ptr%p%ring_redcoord_ik)
    else
        dloc = 0
    end if
end subroutine f90wrap_ham_inp_t__array__ring_redcoord_ik

subroutine f90wrap_ham_inp_t__array__ring_redcoord_il(this, nd, dtype, dshape, dloc)
    use inputdatatype, only: ham_inp_t
    use, intrinsic :: iso_c_binding, only : c_int
    implicit none
    type ham_inp_t_ptr_type
        type(ham_inp_t), pointer :: p => NULL()
    end type ham_inp_t_ptr_type
    integer(c_int), intent(in) :: this(2)
    type(ham_inp_t_ptr_type) :: this_ptr
    integer(c_int), intent(out) :: nd
    integer(c_int), intent(out) :: dtype
    integer(c_int), dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 3
    dtype = 12
    this_ptr = transfer(this, this_ptr)
    if (allocated(this_ptr%p%ring_redcoord_il)) then
        dshape(1:3) = shape(this_ptr%p%ring_redcoord_il)
        dloc = loc(this_ptr%p%ring_redcoord_il)
    else
        dloc = 0
    end if
end subroutine f90wrap_ham_inp_t__array__ring_redcoord_il

subroutine f90wrap_ham_inp_t__array__jc_ring(this, nd, dtype, dshape, dloc)
    use inputdatatype, only: ham_inp_t
    use, intrinsic :: iso_c_binding, only : c_int
    implicit none
    type ham_inp_t_ptr_type
        type(ham_inp_t), pointer :: p => NULL()
    end type ham_inp_t_ptr_type
    integer(c_int), intent(in) :: this(2)
    type(ham_inp_t_ptr_type) :: this_ptr
    integer(c_int), intent(out) :: nd
    integer(c_int), intent(out) :: dtype
    integer(c_int), dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 1
    dtype = 12
    this_ptr = transfer(this, this_ptr)
    if (allocated(this_ptr%p%jc_ring)) then
        dshape(1:1) = shape(this_ptr%p%jc_ring)
        dloc = loc(this_ptr%p%jc_ring)
    else
        dloc = 0
    end if
end subroutine f90wrap_ham_inp_t__array__jc_ring

subroutine f90wrap_ham_inp_t__get__do_dip(this, f90wrap_do_dip)
    use inputdatatype, only: ham_inp_t
    implicit none
    type ham_inp_t_ptr_type
        type(ham_inp_t), pointer :: p => NULL()
    end type ham_inp_t_ptr_type
    integer, intent(in)   :: this(2)
    type(ham_inp_t_ptr_type) :: this_ptr
    integer, intent(out) :: f90wrap_do_dip
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_do_dip = this_ptr%p%do_dip
end subroutine f90wrap_ham_inp_t__get__do_dip

subroutine f90wrap_ham_inp_t__set__do_dip(this, f90wrap_do_dip)
    use inputdatatype, only: ham_inp_t
    implicit none
    type ham_inp_t_ptr_type
        type(ham_inp_t), pointer :: p => NULL()
    end type ham_inp_t_ptr_type
    integer, intent(in)   :: this(2)
    type(ham_inp_t_ptr_type) :: this_ptr
    integer, intent(in) :: f90wrap_do_dip
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%do_dip = f90wrap_do_dip
end subroutine f90wrap_ham_inp_t__set__do_dip

subroutine f90wrap_ham_inp_t__get__read_dipole(this, f90wrap_read_dipole)
    use inputdatatype, only: ham_inp_t
    implicit none
    type ham_inp_t_ptr_type
        type(ham_inp_t), pointer :: p => NULL()
    end type ham_inp_t_ptr_type
    integer, intent(in)   :: this(2)
    type(ham_inp_t_ptr_type) :: this_ptr
    character(1), intent(out) :: f90wrap_read_dipole
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_read_dipole = this_ptr%p%read_dipole
end subroutine f90wrap_ham_inp_t__get__read_dipole

subroutine f90wrap_ham_inp_t__set__read_dipole(this, f90wrap_read_dipole)
    use inputdatatype, only: ham_inp_t
    implicit none
    type ham_inp_t_ptr_type
        type(ham_inp_t), pointer :: p => NULL()
    end type ham_inp_t_ptr_type
    integer, intent(in)   :: this(2)
    type(ham_inp_t_ptr_type) :: this_ptr
    character(1), intent(in) :: f90wrap_read_dipole
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%read_dipole = f90wrap_read_dipole
end subroutine f90wrap_ham_inp_t__set__read_dipole

subroutine f90wrap_ham_inp_t__get__print_dip_tensor(this, f90wrap_print_dip_tensor)
    use inputdatatype, only: ham_inp_t
    implicit none
    type ham_inp_t_ptr_type
        type(ham_inp_t), pointer :: p => NULL()
    end type ham_inp_t_ptr_type
    integer, intent(in)   :: this(2)
    type(ham_inp_t_ptr_type) :: this_ptr
    character(1), intent(out) :: f90wrap_print_dip_tensor
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_print_dip_tensor = this_ptr%p%print_dip_tensor
end subroutine f90wrap_ham_inp_t__get__print_dip_tensor

subroutine f90wrap_ham_inp_t__set__print_dip_tensor(this, f90wrap_print_dip_tensor)
    use inputdatatype, only: ham_inp_t
    implicit none
    type ham_inp_t_ptr_type
        type(ham_inp_t), pointer :: p => NULL()
    end type ham_inp_t_ptr_type
    integer, intent(in)   :: this(2)
    type(ham_inp_t_ptr_type) :: this_ptr
    character(1), intent(in) :: f90wrap_print_dip_tensor
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%print_dip_tensor = f90wrap_print_dip_tensor
end subroutine f90wrap_ham_inp_t__set__print_dip_tensor

subroutine f90wrap_ham_inp_t__get__qdip_files(this, f90wrap_qdip_files)
    use inputdatatype, only: ham_inp_t
    implicit none
    type ham_inp_t_ptr_type
        type(ham_inp_t), pointer :: p => NULL()
    end type ham_inp_t_ptr_type
    integer, intent(in)   :: this(2)
    type(ham_inp_t_ptr_type) :: this_ptr
    character(30), intent(out) :: f90wrap_qdip_files
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_qdip_files = this_ptr%p%qdip_files
end subroutine f90wrap_ham_inp_t__get__qdip_files

subroutine f90wrap_ham_inp_t__set__qdip_files(this, f90wrap_qdip_files)
    use inputdatatype, only: ham_inp_t
    implicit none
    type ham_inp_t_ptr_type
        type(ham_inp_t), pointer :: p => NULL()
    end type ham_inp_t_ptr_type
    integer, intent(in)   :: this(2)
    type(ham_inp_t_ptr_type) :: this_ptr
    character(30), intent(in) :: f90wrap_qdip_files
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%qdip_files = f90wrap_qdip_files
end subroutine f90wrap_ham_inp_t__set__qdip_files

subroutine f90wrap_ham_inp_t__get__do_ewald(this, f90wrap_do_ewald)
    use inputdatatype, only: ham_inp_t
    implicit none
    type ham_inp_t_ptr_type
        type(ham_inp_t), pointer :: p => NULL()
    end type ham_inp_t_ptr_type
    integer, intent(in)   :: this(2)
    type(ham_inp_t_ptr_type) :: this_ptr
    character(1), intent(out) :: f90wrap_do_ewald
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_do_ewald = this_ptr%p%do_ewald
end subroutine f90wrap_ham_inp_t__get__do_ewald

subroutine f90wrap_ham_inp_t__set__do_ewald(this, f90wrap_do_ewald)
    use inputdatatype, only: ham_inp_t
    implicit none
    type ham_inp_t_ptr_type
        type(ham_inp_t), pointer :: p => NULL()
    end type ham_inp_t_ptr_type
    integer, intent(in)   :: this(2)
    type(ham_inp_t_ptr_type) :: this_ptr
    character(1), intent(in) :: f90wrap_do_ewald
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%do_ewald = f90wrap_do_ewald
end subroutine f90wrap_ham_inp_t__set__do_ewald

subroutine f90wrap_ham_inp_t__array__RMAX(this, nd, dtype, dshape, dloc)
    use inputdatatype, only: ham_inp_t
    use, intrinsic :: iso_c_binding, only : c_int
    implicit none
    type ham_inp_t_ptr_type
        type(ham_inp_t), pointer :: p => NULL()
    end type ham_inp_t_ptr_type
    integer(c_int), intent(in) :: this(2)
    type(ham_inp_t_ptr_type) :: this_ptr
    integer(c_int), intent(out) :: nd
    integer(c_int), intent(out) :: dtype
    integer(c_int), dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 1
    dtype = 5
    this_ptr = transfer(this, this_ptr)
    dshape(1:1) = shape(this_ptr%p%RMAX)
    dloc = loc(this_ptr%p%RMAX)
end subroutine f90wrap_ham_inp_t__array__RMAX

subroutine f90wrap_ham_inp_t__array__KMAX(this, nd, dtype, dshape, dloc)
    use inputdatatype, only: ham_inp_t
    use, intrinsic :: iso_c_binding, only : c_int
    implicit none
    type ham_inp_t_ptr_type
        type(ham_inp_t), pointer :: p => NULL()
    end type ham_inp_t_ptr_type
    integer(c_int), intent(in) :: this(2)
    type(ham_inp_t_ptr_type) :: this_ptr
    integer(c_int), intent(out) :: nd
    integer(c_int), intent(out) :: dtype
    integer(c_int), dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 1
    dtype = 5
    this_ptr = transfer(this, this_ptr)
    dshape(1:1) = shape(this_ptr%p%KMAX)
    dloc = loc(this_ptr%p%KMAX)
end subroutine f90wrap_ham_inp_t__array__KMAX

subroutine f90wrap_ham_inp_t__get__Ewald_alpha(this, f90wrap_Ewald_alpha)
    use inputdatatype, only: ham_inp_t
    implicit none
    type ham_inp_t_ptr_type
        type(ham_inp_t), pointer :: p => NULL()
    end type ham_inp_t_ptr_type
    integer, intent(in)   :: this(2)
    type(ham_inp_t_ptr_type) :: this_ptr
    real(8), intent(out) :: f90wrap_Ewald_alpha
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_Ewald_alpha = this_ptr%p%Ewald_alpha
end subroutine f90wrap_ham_inp_t__get__Ewald_alpha

subroutine f90wrap_ham_inp_t__set__Ewald_alpha(this, f90wrap_Ewald_alpha)
    use inputdatatype, only: ham_inp_t
    implicit none
    type ham_inp_t_ptr_type
        type(ham_inp_t), pointer :: p => NULL()
    end type ham_inp_t_ptr_type
    integer, intent(in)   :: this(2)
    type(ham_inp_t_ptr_type) :: this_ptr
    real(8), intent(in) :: f90wrap_Ewald_alpha
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%Ewald_alpha = f90wrap_Ewald_alpha
end subroutine f90wrap_ham_inp_t__set__Ewald_alpha

subroutine f90wrap_ham_inp_t_initialise(this)
    use inputdatatype, only: ham_inp_t
    implicit none
    
    type ham_inp_t_ptr_type
        type(ham_inp_t), pointer :: p => NULL()
    end type ham_inp_t_ptr_type
    type(ham_inp_t_ptr_type) :: this_ptr
    integer, intent(out), dimension(2) :: this
    allocate(this_ptr%p)
    this = transfer(this_ptr, this)
end subroutine f90wrap_ham_inp_t_initialise

subroutine f90wrap_ham_inp_t_finalise(this)
    use inputdatatype, only: ham_inp_t
    implicit none
    
    type ham_inp_t_ptr_type
        type(ham_inp_t), pointer :: p => NULL()
    end type ham_inp_t_ptr_type
    type(ham_inp_t_ptr_type) :: this_ptr
    integer, intent(in), dimension(2) :: this
    this_ptr = transfer(this, this_ptr)
    deallocate(this_ptr%p)
end subroutine f90wrap_ham_inp_t_finalise

! End of module inputdatatype defined in file /Users/andersb/Jobb/UppASD_release/UppASD_release/source/Input/inputdatatype.f90

