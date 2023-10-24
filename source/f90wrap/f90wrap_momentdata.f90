! Module momentdata defined in file System/momentdata.f90

subroutine f90wrap_allocate_mmoms(natom, mensemble, flag)
    use momentdata, only: allocate_mmoms
    implicit none
    
    integer, intent(in), optional :: natom
    integer, intent(in), optional :: mensemble
    integer, intent(in) :: flag
    call allocate_mmoms(Natom=natom, Mensemble=mensemble, flag=flag)
end subroutine f90wrap_allocate_mmoms

subroutine f90wrap_allocate_emoms(natom, mensemble, flag)
    use momentdata, only: allocate_emoms
    implicit none
    
    integer, intent(in) :: natom
    integer, intent(in) :: mensemble
    integer, intent(in) :: flag
    call allocate_emoms(Natom=natom, Mensemble=mensemble, flag=flag)
end subroutine f90wrap_allocate_emoms

subroutine f90wrap_momentdata__array__emom(dummy_this, nd, dtype, dshape, dloc)
    use parameters
    use profiling
    use momentdata, only: momentdata_emom => emom
    use, intrinsic :: iso_c_binding, only : c_int
    implicit none
    integer, intent(in) :: dummy_this(2)
    integer(c_int), intent(out) :: nd
    integer(c_int), intent(out) :: dtype
    integer(c_int), dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 3
    dtype = 12
    if (allocated(momentdata_emom)) then
        dshape(1:3) = shape(momentdata_emom)
        dloc = loc(momentdata_emom)
    else
        dloc = 0
    end if
end subroutine f90wrap_momentdata__array__emom

subroutine f90wrap_momentdata__array__mmom(dummy_this, nd, dtype, dshape, dloc)
    use parameters
    use profiling
    use momentdata, only: momentdata_mmom => mmom
    use, intrinsic :: iso_c_binding, only : c_int
    implicit none
    integer, intent(in) :: dummy_this(2)
    integer(c_int), intent(out) :: nd
    integer(c_int), intent(out) :: dtype
    integer(c_int), dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 2
    dtype = 12
    if (allocated(momentdata_mmom)) then
        dshape(1:2) = shape(momentdata_mmom)
        dloc = loc(momentdata_mmom)
    else
        dloc = 0
    end if
end subroutine f90wrap_momentdata__array__mmom

subroutine f90wrap_momentdata__array__mmomi(dummy_this, nd, dtype, dshape, dloc)
    use parameters
    use profiling
    use momentdata, only: momentdata_mmomi => mmomi
    use, intrinsic :: iso_c_binding, only : c_int
    implicit none
    integer, intent(in) :: dummy_this(2)
    integer(c_int), intent(out) :: nd
    integer(c_int), intent(out) :: dtype
    integer(c_int), dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 2
    dtype = 12
    if (allocated(momentdata_mmomi)) then
        dshape(1:2) = shape(momentdata_mmomi)
        dloc = loc(momentdata_mmomi)
    else
        dloc = 0
    end if
end subroutine f90wrap_momentdata__array__mmomi

subroutine f90wrap_momentdata__array__emom2(dummy_this, nd, dtype, dshape, dloc)
    use parameters
    use profiling
    use momentdata, only: momentdata_emom2 => emom2
    use, intrinsic :: iso_c_binding, only : c_int
    implicit none
    integer, intent(in) :: dummy_this(2)
    integer(c_int), intent(out) :: nd
    integer(c_int), intent(out) :: dtype
    integer(c_int), dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 3
    dtype = 12
    if (allocated(momentdata_emom2)) then
        dshape(1:3) = shape(momentdata_emom2)
        dloc = loc(momentdata_emom2)
    else
        dloc = 0
    end if
end subroutine f90wrap_momentdata__array__emom2

subroutine f90wrap_momentdata__array__emomM(dummy_this, nd, dtype, dshape, dloc)
    use parameters
    use profiling
    use momentdata, only: momentdata_emomm => emomm
    use, intrinsic :: iso_c_binding, only : c_int
    implicit none
    integer, intent(in) :: dummy_this(2)
    integer(c_int), intent(out) :: nd
    integer(c_int), intent(out) :: dtype
    integer(c_int), dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 3
    dtype = 12
    if (allocated(momentdata_emomM)) then
        dshape(1:3) = shape(momentdata_emomM)
        dloc = loc(momentdata_emomM)
    else
        dloc = 0
    end if
end subroutine f90wrap_momentdata__array__emomM

subroutine f90wrap_momentdata__array__mmom2(dummy_this, nd, dtype, dshape, dloc)
    use parameters
    use profiling
    use momentdata, only: momentdata_mmom2 => mmom2
    use, intrinsic :: iso_c_binding, only : c_int
    implicit none
    integer, intent(in) :: dummy_this(2)
    integer(c_int), intent(out) :: nd
    integer(c_int), intent(out) :: dtype
    integer(c_int), dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 2
    dtype = 12
    if (allocated(momentdata_mmom2)) then
        dshape(1:2) = shape(momentdata_mmom2)
        dloc = loc(momentdata_mmom2)
    else
        dloc = 0
    end if
end subroutine f90wrap_momentdata__array__mmom2

subroutine f90wrap_momentdata__array__mmom0(dummy_this, nd, dtype, dshape, dloc)
    use parameters
    use profiling
    use momentdata, only: momentdata_mmom0 => mmom0
    use, intrinsic :: iso_c_binding, only : c_int
    implicit none
    integer, intent(in) :: dummy_this(2)
    integer(c_int), intent(out) :: nd
    integer(c_int), intent(out) :: dtype
    integer(c_int), dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 2
    dtype = 12
    if (allocated(momentdata_mmom0)) then
        dshape(1:2) = shape(momentdata_mmom0)
        dloc = loc(momentdata_mmom0)
    else
        dloc = 0
    end if
end subroutine f90wrap_momentdata__array__mmom0

! End of module momentdata defined in file System/momentdata.f90

