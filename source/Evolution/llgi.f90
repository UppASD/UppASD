!-------------------------------------------------------------------------------
! MODULE: LLGI
!> \brief Heun like solver for the stochastic LLGI-equation
!! \details Heun type solver for the generalized LLGI equation
!! with magnetic accleration terms. If relaxtime is put to zero it correponds to standard LLG eq.
!! \todo Replace unit length moment vectors emom with full lenght vector emomM
!> @author
!> Lars Bergvist
!> @copyright
!! Copyright (C) 2008-2018 UppASD group
!! This file is distributed under the terms of the
!! GNU General Public License.
!! See http://www.gnu.org/copyleft/gpl.txt
!-------------------------------------------------------------------------------
module LLGI
   use Parameters
   use Profiling
   use ErrorHandling

   implicit none
   !
   real(dblprec), dimension(:,:,:), allocatable :: btherm !< Thermal stochastic field
   real(dblprec), dimension(:,:,:), allocatable :: mtmp  !< Temporary new solution m(t+dt)
   real(dblprec), dimension(:,:,:), allocatable :: mpast  !< Saving past solution m(t-dt)
   real(dblprec), dimension(:,:,:), allocatable :: bdup !< Resulting effective field

   private

   public :: llgi_evolve_first, llgi_evolve_second, allocate_llgifields


contains


   !-----------------------------------------------------------------------------
   ! SUBROUTINE: llgi_evolve_first
   !> @brief First step (predictor) of the Heun solver for the LLG-I eq, calculates the stochastic field preform update
   !-----------------------------------------------------------------------------
   subroutine llgi_evolve_first(Natom, Mensemble, lambda1_array, beff, b2eff, emom, emom2, emomM, mmom, delta_t, relaxtime,Temp_array,temprescale, thermal_field)
      !
      implicit none
      !
      integer, intent(in) :: Natom !< Number of atoms in system
      integer, intent(in) :: Mensemble !< Number of ensembles
      real(dblprec), dimension(Natom), intent(in) :: lambda1_array !< Damping parameter
      real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: beff !< Total effective field from application of Hamiltonian
      real(dblprec), dimension(3,Natom,Mensemble), intent(out) :: b2eff !< Temporary storage of magnetic field
      real(dblprec), dimension(3,Natom,Mensemble), intent(inout) :: thermal_field
      real(dblprec), dimension(3,Natom,Mensemble), intent(out) :: emom   !< Current unit moment vector m(t),mp(t+1)
      real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: emom2  !< Final (or temporary) unit moment vector m(t)
      real(dblprec), dimension(3,Natom,Mensemble), intent(out) :: emomM  !< Current magnetic moment vector
      real(dblprec), dimension(Natom,Mensemble), intent(in) :: mmom !< Magnitude of magnetic moments
      real(dblprec), intent(in) :: delta_t !< Time step
      real(dblprec), intent(in) :: relaxtime !< Relaxation time for inertial regime
      real(dblprec), dimension(Natom), intent(in) :: Temp_array !< Temperature (array)
      real(dblprec), intent(inout) :: temprescale  !< Temperature rescaling from QHB

      integer :: i,k

      call ErrorHandling_missing('Inertia')
      b2eff=0.0d0;emom=0.0d0;emomM=0.0d0;

   end subroutine llgi_evolve_first


   !-----------------------------------------------------------------------------
   ! SUBROUTINE: llgi_evolve_second
   !> @brief Second step of solver, calculates the corrected effective field from
   !> the predicted effective fields.
   !-----------------------------------------------------------------------------
   subroutine llgi_evolve_second(Natom, Mensemble, lambda1_array, beff, b2eff, emom, emom2, delta_t,relaxtime)
      !
      implicit none
      !
      integer, intent(in) :: Natom !< Number of atoms in system
      integer, intent(in) :: Mensemble !< Number of ensembles
      real(dblprec), dimension(Natom), intent(in) :: lambda1_array !< Damping parameter
      real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: beff !< Total effective field from application of Hamiltonian
      real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: b2eff !< Temporary storage of magnetic field
      real(dblprec), dimension(3,Natom,Mensemble), intent(out) :: emom   !< Current unit moment vector
      real(dblprec), dimension(3,Natom,Mensemble), intent(inout) :: emom2  !< Final (or temporary) unit moment vector
      real(dblprec), intent(in) :: delta_t !< Time step
      real(dblprec), intent(in) :: relaxtime !< Relaxation time for inertial regime
      !
      integer :: i,k
      real(dblprec) :: etot

      call ErrorHandling_missing('Inertia')
      emom=0.0d0

   end subroutine llgi_evolve_second

   !
   !-----------------------------------------------------------------------------
   ! SUBROUTINE: thermfield
   !> @brief Calculates stochastic field
   !> @todo Confirm the magnitude of the stochastic field as a function
   !> of temperature
   !-----------------------------------------------------------------------------
   subroutine thermfield(Natom, Mensemble, lambda1_array, mmom, deltat,Temp_array,temprescale)
      !
      use Constants, only : k_bolt, gama, mub
      use RandomNumbers, only : rng_gaussian

      implicit none

      integer, intent(in) :: Natom !< Number of atoms in system
      integer, intent(in) :: Mensemble !< Number of ensembles
      real(dblprec), dimension(Natom), intent(in) :: lambda1_array !< Damping parameter
      real(dblprec), dimension(Natom,Mensemble), intent(in) :: mmom !< Magnitude of magnetic moments
      real(dblprec), intent(in) :: deltat !< Time step
      real(dblprec), dimension(Natom) ,intent(in) :: Temp_array  !< Temperature (array)
      real(dblprec), intent(in) :: temprescale  !< Temperature rescaling from QHB

      real(dblprec), dimension(Natom) :: Dp
      real(dblprec) :: mu, sigma

      integer :: i,k

      call ErrorHandling_missing('Inertia')

   end subroutine thermfield


   !-----------------------------------------------------------------------------
   ! SUBROUTINE: llgiheun
   !> @brief Heun solver for LLG-I equation
   !-----------------------------------------------------------------------------
   subroutine llgiheun(Natom, Mensemble,lambda1_array,delta_t,relaxtime,emom,emom2)
      !
      use Constants, only : gama
      implicit none

      integer, intent(in) :: Natom !< Number of atoms in system
      integer, intent(in) :: Mensemble !< Number of ensembles
      real(dblprec), dimension(Natom), intent(in) :: lambda1_array !< Damping parameter
      real(dblprec), intent(in) :: delta_t !< Time step
      real(dblprec), intent(in) :: relaxtime !< Relaxation time for inertial regime
      real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: emom   !< Current unit moment vector, past timestep
      real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: emom2  !< Final (or temporary) unit moment vector
      !
      !   Local variables
      !
      real(dblprec), dimension(3,3) :: A,Ainv
      real(dblprec),dimension(3) :: B   ! b vector
      real(dblprec) :: detmatrix,dtg
      real(dblprec) :: e1x, e1y, e1z
      real(dblprec) :: etx, ety, etz
      real(dblprec) :: hx, hy, hz
      real(dblprec), dimension(Natom) :: alfa, beta, zeta

      integer :: i,k

      call ErrorHandling_missing('Inertia')

   end subroutine llgiheun


   !-----------------------------------------------------------------------------
   ! SUBROUTINE: allocate_llgifields
   !> @brief Allocates work arrays for the Heun-I solver
   !-----------------------------------------------------------------------------
   subroutine allocate_llgifields(Natom,Mensemble,flag)

      implicit none

      integer, intent(in) :: Natom !< Number of atoms in system
      integer, intent(in) :: Mensemble !< Number of ensembles
      integer, intent(in) :: flag !< Allocate or deallocate (1/-1)

      integer :: i_all, i_stat

      if(flag>0) then
         allocate(mtmp(3,Natom,Mensemble),stat=i_stat)
         call memocc(i_stat,product(shape(mtmp))*kind(mtmp),'mtmp','allocate_llgifields')
         allocate(mpast(3,Natom,Mensemble),stat=i_stat)
         call memocc(i_stat,product(shape(mpast))*kind(mpast),'mpast','allocate_depondtfields')
         mpast=0.d0
         allocate(btherm(3,Natom,Mensemble),stat=i_stat)
         call memocc(i_stat,product(shape(btherm))*kind(btherm),'btherm','allocate_llgifields')
         allocate(bdup(3,Natom,Mensemble),stat=i_stat)
         call memocc(i_stat,product(shape(bdup))*kind(bdup),'bdup','allocate_llgifields')
      else
         i_all=-product(shape(mtmp))*kind(mtmp)
         deallocate(mtmp,stat=i_stat)
         call memocc(i_stat,i_all,'mtmp','allocate_systemdata')
         i_all=-product(shape(mpast))*kind(mpast)
         deallocate(mpast,stat=i_stat)
         call memocc(i_stat,i_all,'mpast','allocate_systemdata')
         i_all=-product(shape(btherm))*kind(btherm)
         deallocate(btherm,stat=i_stat)
         call memocc(i_stat,i_all,'btherm','allocate_systemdata')
         i_all=-product(shape(bdup))*kind(bdup)
         deallocate(bdup,stat=i_stat)
         call memocc(i_stat,i_all,'bdup','allocate_systemdata')
      end if

   end subroutine allocate_llgifields


end module LLGI
