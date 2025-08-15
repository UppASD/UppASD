! Module simulationdata defined in file System/simulationdata.f90

subroutine f90wrap_simulationdata__get__lambda1(f90wrap_lambda1)
    use simulationdata, only: simulationdata_lambda1 => lambda1
    implicit none
    real(8), intent(out) :: f90wrap_lambda1
    
    f90wrap_lambda1 = simulationdata_lambda1
end subroutine f90wrap_simulationdata__get__lambda1

subroutine f90wrap_simulationdata__set__lambda1(f90wrap_lambda1)
    use simulationdata, only: simulationdata_lambda1 => lambda1
    implicit none
    real(8), intent(in) :: f90wrap_lambda1
    
    simulationdata_lambda1 = f90wrap_lambda1
end subroutine f90wrap_simulationdata__set__lambda1

subroutine f90wrap_simulationdata__get__lambda2(f90wrap_lambda2)
    use simulationdata, only: simulationdata_lambda2 => lambda2
    implicit none
    real(8), intent(out) :: f90wrap_lambda2
    
    f90wrap_lambda2 = simulationdata_lambda2
end subroutine f90wrap_simulationdata__get__lambda2

subroutine f90wrap_simulationdata__set__lambda2(f90wrap_lambda2)
    use simulationdata, only: simulationdata_lambda2 => lambda2
    implicit none
    real(8), intent(in) :: f90wrap_lambda2
    
    simulationdata_lambda2 = f90wrap_lambda2
end subroutine f90wrap_simulationdata__set__lambda2

subroutine f90wrap_simulationdata__get__rstep(f90wrap_rstep)
    use simulationdata, only: simulationdata_rstep => rstep
    implicit none
    integer, intent(out) :: f90wrap_rstep
    
    f90wrap_rstep = simulationdata_rstep
end subroutine f90wrap_simulationdata__get__rstep

subroutine f90wrap_simulationdata__set__rstep(f90wrap_rstep)
    use simulationdata, only: simulationdata_rstep => rstep
    implicit none
    integer, intent(in) :: f90wrap_rstep
    
    simulationdata_rstep = f90wrap_rstep
end subroutine f90wrap_simulationdata__set__rstep

subroutine f90wrap_simulationdata__get__mstep(f90wrap_mstep)
    use simulationdata, only: simulationdata_mstep => mstep
    implicit none
    integer, intent(out) :: f90wrap_mstep
    
    f90wrap_mstep = simulationdata_mstep
end subroutine f90wrap_simulationdata__get__mstep

subroutine f90wrap_simulationdata__set__mstep(f90wrap_mstep)
    use simulationdata, only: simulationdata_mstep => mstep
    implicit none
    integer, intent(in) :: f90wrap_mstep
    
    simulationdata_mstep = f90wrap_mstep
end subroutine f90wrap_simulationdata__set__mstep

subroutine f90wrap_simulationdata__get__bn(f90wrap_bn)
    use simulationdata, only: simulationdata_bn => bn
    implicit none
    real(8), intent(out) :: f90wrap_bn
    
    f90wrap_bn = simulationdata_bn
end subroutine f90wrap_simulationdata__get__bn

subroutine f90wrap_simulationdata__get__total_energy(f90wrap_total_energy)
    use simulationdata, only: simulationdata_total_energy => total_energy
    implicit none
    real(8), intent(out) :: f90wrap_total_energy
    
    f90wrap_total_energy = simulationdata_total_energy
end subroutine f90wrap_simulationdata__get__total_energy

subroutine f90wrap_simulationdata__set__total_energy(f90wrap_total_energy)
    use simulationdata, only: simulationdata_total_energy => total_energy
    implicit none
    real(8), intent(in) :: f90wrap_total_energy
    
    simulationdata_total_energy = f90wrap_total_energy
end subroutine f90wrap_simulationdata__set__total_energy

! End of module simulationdata defined in file System/simulationdata.f90

