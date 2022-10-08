! Module uppasd defined in file /Users/andersb/Jobb/UppASD_release/UppASD_release/source/uppasd.f90

subroutine f90wrap_main
    use uppasd, only: main
    implicit none
    
    call main()
end subroutine f90wrap_main

subroutine f90wrap_number_of_active_processors(ret_nprocs)
    use uppasd, only: number_of_active_processors
    implicit none
    
    integer, intent(out) :: ret_nprocs
    ret_nprocs = number_of_active_processors()
end subroutine f90wrap_number_of_active_processors

subroutine f90wrap_run_initial_phase
    use uppasd, only: run_initial_phase
    implicit none
    
    call run_initial_phase()
end subroutine f90wrap_run_initial_phase

subroutine f90wrap_run_measurement_phase
    use uppasd, only: run_measurement_phase
    implicit none
    
    call run_measurement_phase()
end subroutine f90wrap_run_measurement_phase

subroutine f90wrap_cleanup_simulation
    use uppasd, only: cleanup_simulation
    implicit none
    
    call cleanup_simulation()
end subroutine f90wrap_cleanup_simulation

subroutine f90wrap_setup_simulation
    use uppasd, only: setup_simulation
    implicit none
    
    call setup_simulation()
end subroutine f90wrap_setup_simulation

subroutine f90wrap_sd_timing(nstep)
    use uppasd, only: sd_timing
    implicit none
    
    integer, intent(in) :: nstep
    call sd_timing(Nstep=nstep)
end subroutine f90wrap_sd_timing

subroutine f90wrap_allocate_general(flag)
    use uppasd, only: allocate_general
    implicit none
    
    integer, intent(in) :: flag
    call allocate_general(flag=flag)
end subroutine f90wrap_allocate_general

subroutine f90wrap_print_logo
    use uppasd, only: print_logo
    implicit none
    
    call print_logo()
end subroutine f90wrap_print_logo

subroutine f90wrap_print_siminfo
    use uppasd, only: print_siminfo
    implicit none
    
    call print_siminfo()
end subroutine f90wrap_print_siminfo

subroutine f90wrap_check_format
    use uppasd, only: check_format
    implicit none
    
    call check_format()
end subroutine f90wrap_check_format

subroutine f90wrap_calculate_energy(outenergy)
    use uppasd, only: calculate_energy
    implicit none
    
    real(8), optional :: outenergy
    call calculate_energy(outenergy=outenergy)
end subroutine f90wrap_calculate_energy

subroutine f90wrap_uppasd__get__time_a(f90wrap_time_a)
    use uppasd, only: uppasd_time_a => time_a
    implicit none
    real(8), intent(out) :: f90wrap_time_a
    
    f90wrap_time_a = uppasd_time_a
end subroutine f90wrap_uppasd__get__time_a

subroutine f90wrap_uppasd__set__time_a(f90wrap_time_a)
    use uppasd, only: uppasd_time_a => time_a
    implicit none
    real(8), intent(in) :: f90wrap_time_a
    
    uppasd_time_a = f90wrap_time_a
end subroutine f90wrap_uppasd__set__time_a

subroutine f90wrap_uppasd__get__time_b(f90wrap_time_b)
    use uppasd, only: uppasd_time_b => time_b
    implicit none
    real(8), intent(out) :: f90wrap_time_b
    
    f90wrap_time_b = uppasd_time_b
end subroutine f90wrap_uppasd__get__time_b

subroutine f90wrap_uppasd__set__time_b(f90wrap_time_b)
    use uppasd, only: uppasd_time_b => time_b
    implicit none
    real(8), intent(in) :: f90wrap_time_b
    
    uppasd_time_b = f90wrap_time_b
end subroutine f90wrap_uppasd__set__time_b

subroutine f90wrap_uppasd__get__time_c(f90wrap_time_c)
    use uppasd, only: uppasd_time_c => time_c
    implicit none
    real(8), intent(out) :: f90wrap_time_c
    
    f90wrap_time_c = uppasd_time_c
end subroutine f90wrap_uppasd__get__time_c

subroutine f90wrap_uppasd__set__time_c(f90wrap_time_c)
    use uppasd, only: uppasd_time_c => time_c
    implicit none
    real(8), intent(in) :: f90wrap_time_c
    
    uppasd_time_c = f90wrap_time_c
end subroutine f90wrap_uppasd__set__time_c

subroutine f90wrap_uppasd__get__time_d(f90wrap_time_d)
    use uppasd, only: uppasd_time_d => time_d
    implicit none
    real(8), intent(out) :: f90wrap_time_d
    
    f90wrap_time_d = uppasd_time_d
end subroutine f90wrap_uppasd__get__time_d

subroutine f90wrap_uppasd__set__time_d(f90wrap_time_d)
    use uppasd, only: uppasd_time_d => time_d
    implicit none
    real(8), intent(in) :: f90wrap_time_d
    
    uppasd_time_d = f90wrap_time_d
end subroutine f90wrap_uppasd__set__time_d

subroutine f90wrap_uppasd__get__time_e(f90wrap_time_e)
    use uppasd, only: uppasd_time_e => time_e
    implicit none
    real(8), intent(out) :: f90wrap_time_e
    
    f90wrap_time_e = uppasd_time_e
end subroutine f90wrap_uppasd__get__time_e

subroutine f90wrap_uppasd__set__time_e(f90wrap_time_e)
    use uppasd, only: uppasd_time_e => time_e
    implicit none
    real(8), intent(in) :: f90wrap_time_e
    
    uppasd_time_e = f90wrap_time_e
end subroutine f90wrap_uppasd__set__time_e

subroutine f90wrap_uppasd__get__time_I(f90wrap_time_I)
    use uppasd, only: uppasd_time_I => time_I
    implicit none
    real(8), intent(out) :: f90wrap_time_I
    
    f90wrap_time_I = uppasd_time_I
end subroutine f90wrap_uppasd__get__time_I

subroutine f90wrap_uppasd__set__time_I(f90wrap_time_I)
    use uppasd, only: uppasd_time_I => time_I
    implicit none
    real(8), intent(in) :: f90wrap_time_I
    
    uppasd_time_I = f90wrap_time_I
end subroutine f90wrap_uppasd__set__time_I

subroutine f90wrap_uppasd__get__time_II(f90wrap_time_II)
    use uppasd, only: uppasd_time_II => time_II
    implicit none
    real(8), intent(out) :: f90wrap_time_II
    
    f90wrap_time_II = uppasd_time_II
end subroutine f90wrap_uppasd__get__time_II

subroutine f90wrap_uppasd__set__time_II(f90wrap_time_II)
    use uppasd, only: uppasd_time_II => time_II
    implicit none
    real(8), intent(in) :: f90wrap_time_II
    
    uppasd_time_II = f90wrap_time_II
end subroutine f90wrap_uppasd__set__time_II

subroutine f90wrap_uppasd__get__time_III(f90wrap_time_III)
    use uppasd, only: uppasd_time_III => time_III
    implicit none
    real(8), intent(out) :: f90wrap_time_III
    
    f90wrap_time_III = uppasd_time_III
end subroutine f90wrap_uppasd__get__time_III

subroutine f90wrap_uppasd__set__time_III(f90wrap_time_III)
    use uppasd, only: uppasd_time_III => time_III
    implicit none
    real(8), intent(in) :: f90wrap_time_III
    
    uppasd_time_III = f90wrap_time_III
end subroutine f90wrap_uppasd__set__time_III

subroutine f90wrap_uppasd__get__time_IV(f90wrap_time_IV)
    use uppasd, only: uppasd_time_IV => time_IV
    implicit none
    real(8), intent(out) :: f90wrap_time_IV
    
    f90wrap_time_IV = uppasd_time_IV
end subroutine f90wrap_uppasd__get__time_IV

subroutine f90wrap_uppasd__set__time_IV(f90wrap_time_IV)
    use uppasd, only: uppasd_time_IV => time_IV
    implicit none
    real(8), intent(in) :: f90wrap_time_IV
    
    uppasd_time_IV = f90wrap_time_IV
end subroutine f90wrap_uppasd__set__time_IV

subroutine f90wrap_uppasd__get__time_V(f90wrap_time_V)
    use uppasd, only: uppasd_time_V => time_V
    implicit none
    real(8), intent(out) :: f90wrap_time_V
    
    f90wrap_time_V = uppasd_time_V
end subroutine f90wrap_uppasd__get__time_V

subroutine f90wrap_uppasd__set__time_V(f90wrap_time_V)
    use uppasd, only: uppasd_time_V => time_V
    implicit none
    real(8), intent(in) :: f90wrap_time_V
    
    uppasd_time_V = f90wrap_time_V
end subroutine f90wrap_uppasd__set__time_V

subroutine f90wrap_uppasd__get__nprocs(f90wrap_nprocs)
    use uppasd, only: uppasd_nprocs => nprocs
    implicit none
    integer, intent(out) :: f90wrap_nprocs
    
    f90wrap_nprocs = uppasd_nprocs
end subroutine f90wrap_uppasd__get__nprocs

subroutine f90wrap_uppasd__set__nprocs(f90wrap_nprocs)
    use uppasd, only: uppasd_nprocs => nprocs
    implicit none
    integer, intent(in) :: f90wrap_nprocs
    
    uppasd_nprocs = f90wrap_nprocs
end subroutine f90wrap_uppasd__set__nprocs

! End of module uppasd defined in file /Users/andersb/Jobb/UppASD_release/UppASD_release/source/uppasd.f90

