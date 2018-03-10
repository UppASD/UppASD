!-------------------------------------------------------------------------------
! MODULE: ADAPTIVETIMESTEPPING
!>  Data and routines for adaptive time stepping
!> @author
!> A. Bergman, T. Nystrand etc
!> @copyright
!! Copyright (C) 2008-2018 UppASD group
!! This file is distributed under the terms of the
!! GNU General Public License.
!! See http://www.gnu.org/copyleft/gpl.txt
!-------------------------------------------------------------------------------
module adaptivetimestepping

   use Parameters
   use Profiling
   use Constants
   use Correlation, only : calc_corr_w, w, do_sc, sc_nstep
   use ErrorHandling

   implicit none

   real(dblprec) :: omega_max

   public :: calculate_timestep, calculate_omegainit

contains

   !--------------------------------------------------------------------------
   !
   ! DESCRIPTION
   !> @brief
   !> Main driver for adaptive time stepping based on .....
   !---------------------------------------------------------------------------------
   subroutine adapt_time_step(Natom,Mensemble,beff,omega_max,larmor_numrev,larmor_thr,rstep,mstep,nstep,totalsimtime,&
         therm_fields,do_sc,sc_step,sc_tidx,sd_phaseflag,adapt_time_interval,&
         adapt_step,adaptive_time_flag,deltat_correction_flag,delta_t)

      implicit none

      integer, intent(in) :: Natom                                            !< Number of atoms in the system
      integer, intent(in) :: Mensemble                                        !< Number of ensembles
      integer, intent(in) :: larmor_numrev                                    !< Number of time steps per revolution needed to suppress numerical damping
      integer, intent(in) :: rstep, mstep                                     !< Simulation parameters from 0sd from the MP
      integer, intent(in) :: sc_tidx                                          !< Current sampling time indicator
      integer, intent(in) :: adapt_step
      integer, intent(in) :: adapt_time_interval
      integer, intent(inout) :: sc_step                                       !< Sampling period between spin correlation measurements
      integer, intent(inout) :: nstep                                   !< Number of steps in measurement phase
      character(len=1), intent(in) :: do_sc

      logical, intent(in) :: sd_phaseflag                                     !< Spin dynamics phase indicator
      logical, intent(in) :: adaptive_time_flag
      logical, intent(inout) :: deltat_correction_flag

      real(dblprec), intent(in) :: larmor_thr                                 !< Larmor frequency threshold
      real(dblprec), intent(inout) :: delta_t                                 !< Time step of the simulation
      real(dblprec), intent(inout) :: omega_max                               !< Previously saved maximum angular frequency
      real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: beff         !< Total effective field from application of Hamiltonian
      real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: therm_fields !< Thermal stochastic field
      real(dblprec) :: totalsimtime                                           !< Total simulation time

      call ErrorHandling_missing('Adaptive time')

   end subroutine adapt_time_step

   !> Calculate new timestep
   subroutine calculate_timestep(Natom, Mensemble, beff, omega_max, larmor_numrev, larmor_thr, mstep, nstep, totalsimtime, &
         therm_fields, sc_step, sc_tidx, sd_phaseflag,delta_t)

      implicit none

      integer, intent(in) :: Natom                                            !< Number of atoms in the system
      integer, intent(in) :: Mensemble                                        !< Number of ensembles
      integer, intent(in) :: larmor_numrev                                    !< Number of time steps per revolution needed to suppress numerical damping
      integer, intent(in) :: mstep                                     !< Simulation parameters from 0sd from the MP
      integer, intent(in) :: sc_tidx                                          !< Current sampling time indicator
      integer, intent(inout) :: sc_step                                       !< Sampling period between spin correlation measurements
      integer, intent(inout) :: nstep                                   !< Number of steps in measurement phase

      logical, intent(in) :: sd_phaseflag                                     !< Spin dynamics phase indicator

      real(dblprec), intent(in) :: larmor_thr                                 !< Larmor frequency threshold
      real(dblprec), intent(inout) :: delta_t                                 !< Time step of the simulation
      real(dblprec), intent(inout) :: omega_max                               !< Previously saved maximum angular frequency
      real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: beff         !< Total effective field from application of Hamiltonian
      real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: therm_fields !< Thermal stochastic field
      real(dblprec) :: totalsimtime                                           !< Total simulation time

      real(dblprec) :: pi, t_frac, frac, delta_t_old, sc_fact
      real(dblprec) :: old_omega_max, relprec
      real(dblprec), dimension(Natom,Mensemble) :: beff_abs
      real(dblprec), dimension(Natom,Mensemble) :: bthm_abs                   !< Magnitude of magnetic moments

      integer :: ik, i, k, timestep_frac                    !< Elapsed time of the simuation
      integer :: int_frac

      call ErrorHandling_missing('Adaptive time')

   end subroutine calculate_timestep

   !> Maximum angular frequency
   subroutine calculate_omegainit(omega_max, larmor_numrev, delta_t)
      implicit none

      real(dblprec), intent(out) :: omega_max                          !< Previously saved maximum angular frequency
      integer, intent(in) :: larmor_numrev                             !< Number of time steps per revolution needed to suppress numerical damping
      real(dblprec), intent(in) :: delta_t                             !< Time step of the simulation

      real(dblprec) :: pi

      omega_max=0.0d0
      call ErrorHandling_missing('Adaptive time')

   end subroutine calculate_omegainit

end module adaptivetimestepping
