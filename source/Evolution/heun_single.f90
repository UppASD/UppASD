!> A single step Heun solver for the LLG-equations
!> @copyright
!! Copyright (C) 2008-2018 UppASD group
!! This file is distributed under the terms of the
!! GNU General Public License. 
!! See http://www.gnu.org/copyleft/gpl.txt
module Heun_single
  use Parameters
  use Profiling
   use ErrorHandling

  implicit none

  private
  public :: evolve2, evolve3


contains 


  !> The evolution routine for single step Heun
  subroutine evolve2(Natom, Mensemble, Landeg, llg, bn, lambda1_array, lambda2_array, beff, beff2, &
             field1, field2, mmomi, emom2, compensate_drift, deltat, Temp_array,temprescale,thermal_field)
    !
    use Constants
    use RandomNumbers, only : ranv, rng_gaussian
    !
    implicit none

    integer, intent(in) :: Natom !< Number of atoms in system
    integer, intent(in) :: Mensemble !< Number of ensembles 
    real(dblprec), dimension(Natom), intent(in) :: Landeg  !< Gyromagnetic ratio
    integer, intent(in) :: llg !< Type of equation of motion (1=LLG)
    real(dblprec), intent(in) :: bn !< Scaling factor for LLG equation (parameter=1.0d0)
    real(dblprec), dimension(Natom), intent(in) :: lambda1_array !< Damping parameter
    real(dblprec), dimension(Natom), intent(in) :: lambda2_array !< Additional damping parameter (not used for llg=1)
    real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: beff !< Total effective field from application of Hamiltonian
    real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: beff2 !< External field from application of Hamiltonian
    real(dblprec), dimension(3,Natom,Mensemble), intent(inout), optional :: thermal_field
    real(dblprec), dimension(3,Mensemble), intent(in) :: field1 !< Average internal effective field
    real(dblprec), dimension(3,Mensemble), intent(in) :: field2 !< Average external effective field
    real(dblprec), dimension(Natom,Mensemble), intent(in) :: mmomi !< Inverse of magnitude of magnetic moments
    real(dblprec), dimension(3,Natom,Mensemble), intent(inout) :: emom2  !< Final (or temporary) unit moment vector
    integer, intent(in) :: compensate_drift !< Correct for drift in RNG
    real(dblprec), intent(in) :: deltat !< Time step
    real(dblprec), dimension(Natom), intent(in) :: Temp_array  !< Temperature
    real(dblprec), intent(in) :: temprescale  !< Temperature rescaling from QHB

    real(dblprec) :: D
    real(dblprec) :: etot
    real(dblprec) :: dt, mu, sigma, avf1, avf2
    real(dblprec) :: tk1,tk2
    real(dblprec) :: rx,ry,rz
    real(dblprec), dimension(Natom) :: Dk

    integer :: i

            call ErrorHandling_missing('Euler solver')

  end subroutine evolve2


  !> An alternative evolution routine for single step Heun
  subroutine evolve3(Natom, Mensemble, Landeg, llg, lambda1_array, beff, field1, field2, &
             mmomi, emom2, compensate_drift, deltat,Temp_array, temprescale, thermal_field) 
    use Constants
    use RandomNumbers, only : rng_norm, ranv

    implicit none

    integer, intent(in) :: Natom !< Number of atoms in system
    integer, intent(in) :: Mensemble !< Number of ensembles 
    real(dblprec), dimension(Natom), intent(in) :: Landeg  !< Gyromagnetic ratio
    integer, intent(in) :: llg !< Type of equation of motion (1=LLG)
    real(dblprec), dimension(Natom), intent(in) :: lambda1_array !< Damping parameter
    real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: beff !< Total effective field from application of Hamiltonian
    real(dblprec), dimension(3,Natom,Mensemble), intent(inout), optional :: thermal_field
    real(dblprec), dimension(3,Mensemble), intent(in) :: field1 !< Average internal effective field
    real(dblprec), dimension(3,Mensemble), intent(in) :: field2 !< Average external effective field
    real(dblprec), dimension(Natom,Mensemble), intent(in) :: mmomi !< Inverse of magnitude of magnetic moments
    real(dblprec), dimension(3,Natom,Mensemble), intent(inout) :: emom2  !< Final (or temporary) unit moment vector
    integer, intent(in) :: compensate_drift !< Correct for drift in RNG
    real(dblprec), intent(in) :: deltat !< Time step
    real(dblprec), dimension(Natom), intent(in) :: Temp_array  !< Temperature
    real(dblprec), intent(in) :: temprescale  !< Temperature rescaling from QHB

    real(dblprec) :: D
    real(dblprec) :: etot
    real(dblprec) :: dt, mu, sigma, avf1, avf2
    real(dblprec) :: tk1
    real(dblprec) :: rx,ry,rz
    real(dblprec), dimension(Natom) :: Dk

    integer :: i

            call ErrorHandling_missing('Euler solver')

  end subroutine evolve3


  !> Perform a single Heun step
  subroutine heun2(Natom, Mensemble, Landeg,lambda1_array, lambda2_array, llg,beff, beff2,emom2,dt,b2h1,b2h2)
    !
    use RandomNumbers, only : ranv
    !
    implicit none

    integer, intent(in) :: Natom !< Number of atoms in system
    integer, intent(in) :: Mensemble !< Number of ensembles 
    real(dblprec), dimension(Natom), intent(in) :: Landeg  !< Gyromagnetic ratio
    real(dblprec), dimension(Natom), intent(in) :: lambda1_array !< Damping parameter
    real(dblprec), dimension(Natom), intent(in) :: lambda2_array !< Additional damping parameter (not used for llg=1)
    integer, intent(in) :: llg !< Type of equation of motion (1=LLG)
    real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: beff !< Total effective field from application of Hamiltonian
    real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: beff2 !< External field from application of Hamiltonian
    real(dblprec), dimension(3,Natom,Mensemble), intent(inout) :: emom2  !< Final (or temporary) unit moment vector
    real(dblprec), intent(in) :: dt !< Time step
    real(dblprec), intent(in) :: b2h1 !< Scale factor for effective field
    real(dblprec), intent(in) :: b2h2 !< Scale factor for applied field

    real(dblprec) :: b1xx, b1xy, b1xz
    real(dblprec) :: b1yx, b1yy, b1yz
    real(dblprec) :: b1zx, b1zy, b1zz
    real(dblprec) :: a1x, a1y, a1z
    real(dblprec) :: etx, ety, etz
    real(dblprec) :: dtg
    !
    real(dblprec) :: dwx, dwy, dwz
    real(dblprec) :: e1x, e1y, e1z
    real(dblprec) :: hlpx,hlpy,hlpz
    real(dblprec) :: hesx,hesy,hesz,prot
    real(dblprec), dimension(Natom) :: clambda1
 
    !
    ! ... Local variables ...
    integer :: i
    !
    !
            call ErrorHandling_missing('Euler solver')

  end subroutine heun2


  !> Perform a single Heun step
  subroutine heun3(Natom, Mensemble, Landeg, lambda1_array, beff, emom2, dt)

    use RandomNumbers, only : ranv

    implicit none

    integer, intent(in) :: Natom !< Number of atoms in system
    integer, intent(in) :: Mensemble !< Number of ensembles 
    real(dblprec), dimension(Natom), intent(in) :: Landeg  !< Gyromagnetic ratio
    real(dblprec), dimension(Natom), intent(in) :: lambda1_array !< Damping parameter
    real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: beff !< Total effective field from application of Hamiltonian
    real(dblprec), dimension(3,Natom,Mensemble), intent(inout) :: emom2  !< Final (or temporary) unit moment vector

    real(dblprec) :: b1xx, b1xy, b1xz
    real(dblprec) :: b1yx, b1yy, b1yz
    real(dblprec) :: b1zx, b1zy, b1zz
    real(dblprec) :: a1x, a1y, a1z
    real(dblprec) :: etx, ety, etz
    real(dblprec), intent(in) :: dt !< Time step
    real(dblprec) :: dtg
    !
    real(dblprec) :: dwx, dwy, dwz
    real(dblprec) :: e1x, e1y, e1z
    real(dblprec) :: hlpx,hlpy,hlpz
    real(dblprec) :: hesx,hesy,hesz,prot
    !
    ! ... Local variables ...
    integer :: i
    !
    !

            call ErrorHandling_missing('Euler solver')

  end subroutine heun3

end module heun_single
