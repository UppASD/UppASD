! Module inputhandler_ext defined in file Input/inputhandler_ext.f90

subroutine f90wrap_read_positions
    use inputhandler_ext, only: read_positions
    implicit none
    
    call read_positions()
end subroutine f90wrap_read_positions

subroutine f90wrap_read_positions_alloy
    use inputhandler_ext, only: read_positions_alloy
    implicit none
    
    call read_positions_alloy()
end subroutine f90wrap_read_positions_alloy

subroutine f90wrap_read_moments(landeg_global)
    use inputhandler_ext, only: read_moments
    implicit none
    
    real(8), intent(in) :: landeg_global
    call read_moments(Landeg_global=landeg_global)
end subroutine f90wrap_read_moments

subroutine f90wrap_read_fixed_moments(landeg_global)
    use inputhandler_ext, only: read_fixed_moments
    implicit none
    
    real(8), intent(in) :: landeg_global
    call read_fixed_moments(Landeg_global=landeg_global)
end subroutine f90wrap_read_fixed_moments

subroutine f90wrap_read_exchange(ham_inp)
    use inputhandler_ext, only: read_exchange
    use inputdatatype, only: ham_inp_t
    implicit none
    
    type ham_inp_t_ptr_type
        type(ham_inp_t), pointer :: p => NULL()
    end type ham_inp_t_ptr_type
    type(ham_inp_t_ptr_type) :: ham_inp_ptr
    integer, intent(in), dimension(2) :: ham_inp
    ham_inp_ptr = transfer(ham_inp, ham_inp_ptr)
    call read_exchange(ham_inp=ham_inp_ptr%p)
end subroutine f90wrap_read_exchange

subroutine f90wrap_read_exchange_tensor
    use inputhandler_ext, only: read_exchange_tensor
    implicit none
    
    call read_exchange_tensor()
end subroutine f90wrap_read_exchange_tensor

subroutine f90wrap_read_exchange_build_tensor
    use inputhandler_ext, only: read_exchange_build_tensor
    implicit none
    
    call read_exchange_build_tensor()
end subroutine f90wrap_read_exchange_build_tensor

subroutine f90wrap_read_exchange_getmaxnoshells(no_shells, flines)
    use inputhandler_ext, only: read_exchange_getmaxnoshells
    implicit none
    
    integer, intent(out) :: no_shells
    integer, intent(out) :: flines
    call read_exchange_getmaxnoshells(no_shells=no_shells, flines=flines)
end subroutine f90wrap_read_exchange_getmaxnoshells

subroutine f90wrap_read_exchange_getneighvec(r_red, r_tmp, isite, jsite)
    use inputhandler_ext, only: read_exchange_getneighvec
    implicit none
    
    real(8), dimension(3), intent(inout) :: r_red
    real(8), dimension(3), intent(in) :: r_tmp
    integer, intent(in) :: isite
    integer, intent(in) :: jsite
    call read_exchange_getneighvec(r_red=r_red, r_tmp=r_tmp, isite=isite, jsite=jsite)
end subroutine f90wrap_read_exchange_getneighvec

subroutine f90wrap_read_anisotropy_alloy
    use inputhandler_ext, only: read_anisotropy_alloy
    implicit none
    
    call read_anisotropy_alloy()
end subroutine f90wrap_read_anisotropy_alloy

subroutine f90wrap_read_anisotropy
    use inputhandler_ext, only: read_anisotropy
    implicit none
    
    call read_anisotropy()
end subroutine f90wrap_read_anisotropy

subroutine f90wrap_read_dmdata
    use inputhandler_ext, only: read_dmdata
    implicit none
    
    call read_dmdata()
end subroutine f90wrap_read_dmdata

subroutine f90wrap_read_sadata
    use inputhandler_ext, only: read_sadata
    implicit none
    
    call read_sadata()
end subroutine f90wrap_read_sadata

subroutine f90wrap_read_chirdata
    use inputhandler_ext, only: read_chirdata
    implicit none
    
    call read_chirdata()
end subroutine f90wrap_read_chirdata

subroutine f90wrap_read_fourxdata
    use inputhandler_ext, only: read_fourxdata
    implicit none
    
    call read_fourxdata()
end subroutine f90wrap_read_fourxdata

subroutine f90wrap_read_pddata
    use inputhandler_ext, only: read_pddata
    implicit none
    
    call read_pddata()
end subroutine f90wrap_read_pddata

subroutine f90wrap_read_biqdmdata
    use inputhandler_ext, only: read_biqdmdata
    implicit none
    
    call read_biqdmdata()
end subroutine f90wrap_read_biqdmdata

subroutine f90wrap_read_bqdata
    use inputhandler_ext, only: read_bqdata
    implicit none
    
    call read_bqdata()
end subroutine f90wrap_read_bqdata

subroutine f90wrap_read_ringdata
    use inputhandler_ext, only: read_ringdata
    implicit none
    
    call read_ringdata()
end subroutine f90wrap_read_ringdata

subroutine f90wrap_read_ip_damping
    use inputhandler_ext, only: read_ip_damping
    implicit none
    
    call read_ip_damping()
end subroutine f90wrap_read_ip_damping

subroutine f90wrap_read_ip_damping_alloy
    use inputhandler_ext, only: read_ip_damping_alloy
    implicit none
    
    call read_ip_damping_alloy()
end subroutine f90wrap_read_ip_damping_alloy

subroutine f90wrap_read_damping
    use inputhandler_ext, only: read_damping
    implicit none
    
    call read_damping()
end subroutine f90wrap_read_damping

subroutine f90wrap_read_damping_alloy
    use inputhandler_ext, only: read_damping_alloy
    implicit none
    
    call read_damping_alloy()
end subroutine f90wrap_read_damping_alloy

subroutine f90wrap_read_barriers
    use inputhandler_ext, only: read_barriers
    implicit none
    
    call read_barriers()
end subroutine f90wrap_read_barriers

! End of module inputhandler_ext defined in file Input/inputhandler_ext.f90

