! Module pyasd defined in file /Users/andersb/Jobb/UppASD_release/UppASD_release/source/pyasd.f90

subroutine f90wrap_runuppasd
    use pyasd, only: runuppasd
    implicit none
    
    call runuppasd()
end subroutine f90wrap_runuppasd

subroutine f90wrap_sanitycheck
    use pyasd, only: sanitycheck
    implicit none
    
    call sanitycheck()
end subroutine f90wrap_sanitycheck

subroutine f90wrap_numprocs(ret_nprocs)
    use pyasd, only: numprocs
    implicit none
    
    integer, intent(out) :: ret_nprocs
    ret_nprocs = numprocs()
end subroutine f90wrap_numprocs

subroutine f90wrap_printlogo
    use pyasd, only: printlogo
    implicit none
    
    call printlogo()
end subroutine f90wrap_printlogo

subroutine f90wrap_setupall
    use pyasd, only: setupall
    implicit none
    
    call setupall()
end subroutine f90wrap_setupall

subroutine f90wrap_initialphase
    use pyasd, only: initialphase
    implicit none
    
    call initialphase()
end subroutine f90wrap_initialphase

subroutine f90wrap_measure
    use pyasd, only: measure
    implicit none
    
    call measure()
end subroutine f90wrap_measure

subroutine f90wrap_cleanup
    use pyasd, only: cleanup
    implicit none
    
    call cleanup()
end subroutine f90wrap_cleanup

subroutine f90wrap_relaxmontecarlo
    use pyasd, only: relaxmontecarlo
    implicit none
    
    call relaxmontecarlo()
end subroutine f90wrap_relaxmontecarlo

subroutine f90wrap_relaxmetropolis
    use pyasd, only: relaxmetropolis
    implicit none
    
    call relaxmetropolis()
end subroutine f90wrap_relaxmetropolis

subroutine f90wrap_relaxheatbath
    use pyasd, only: relaxheatbath
    implicit none
    
    call relaxheatbath()
end subroutine f90wrap_relaxheatbath

subroutine f90wrap_relaxmultiscale
    use pyasd, only: relaxmultiscale
    implicit none
    
    call relaxmultiscale()
end subroutine f90wrap_relaxmultiscale

subroutine f90wrap_relaxllg
    use pyasd, only: relaxllg
    implicit none
    
    call relaxllg()
end subroutine f90wrap_relaxllg

subroutine f90wrap_relaxmd
    use pyasd, only: relaxmd
    implicit none
    
    call relaxmd()
end subroutine f90wrap_relaxmd

subroutine f90wrap_relaxsldllg
    use pyasd, only: relaxsldllg
    implicit none
    
    call relaxsldllg()
end subroutine f90wrap_relaxsldllg

subroutine f90wrap_relaxgneb
    use pyasd, only: relaxgneb
    implicit none
    
    call relaxgneb()
end subroutine f90wrap_relaxgneb

subroutine f90wrap_runmontecarlo
    use pyasd, only: runmontecarlo
    implicit none
    
    call runmontecarlo()
end subroutine f90wrap_runmontecarlo

subroutine f90wrap_runmultiscale
    use pyasd, only: runmultiscale
    implicit none
    
    call runmultiscale()
end subroutine f90wrap_runmultiscale

subroutine f90wrap_runllglite
    use pyasd, only: runllglite
    implicit none
    
    call runllglite()
end subroutine f90wrap_runllglite

subroutine f90wrap_runllg
    use pyasd, only: runllg
    implicit none
    
    call runllg()
end subroutine f90wrap_runllg

subroutine f90wrap_runllgcuda
    use pyasd, only: runllgcuda
    implicit none
    
    call runllgcuda()
end subroutine f90wrap_runllgcuda

subroutine f90wrap_runld
    use pyasd, only: runld
    implicit none
    
    call runld()
end subroutine f90wrap_runld

subroutine f90wrap_runsldllg
    use pyasd, only: runsldllg
    implicit none
    
    call runsldllg()
end subroutine f90wrap_runsldllg

subroutine f90wrap_runsldllgimplicit
    use pyasd, only: runsldllgimplicit
    implicit none
    
    call runsldllgimplicit()
end subroutine f90wrap_runsldllgimplicit

subroutine f90wrap_rungneb
    use pyasd, only: rungneb
    implicit none
    
    call rungneb()
end subroutine f90wrap_rungneb

subroutine f90wrap_totalenergy(ret_energy)
    use pyasd, only: totalenergy
    implicit none
    
    real(8), intent(out) :: ret_energy
    ret_energy = totalenergy()
end subroutine f90wrap_totalenergy

! End of module pyasd defined in file /Users/andersb/Jobb/UppASD_release/UppASD_release/source/pyasd.f90

