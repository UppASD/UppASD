!>  Data and routines for adaptive time stepping
!> @author
!> A. Bergman, T. Nystrand etc
!> @copyright
!> GNU Public License.
module adaptivetimestepping

   use Parameters
   use Profiling
   use Constants
   use Omegas, only : calc_corr_w
!  use Correlation, only : sc, do_sc, uc
!   use CorrelationType

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
         therm_fields,do_sc,sc_step,sc_nstep,sc_tidx,sd_phaseflag,adapt_time_interval,&
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
      integer, intent(inout) :: sc_nstep                                       !< Sampling period between spin correlation measurements
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

      ! Standard case
      if(adaptive_time_flag) then
         if (mod(mstep-rstep-1,adapt_time_interval)==0) then
            call calculate_timestep(Natom, Mensemble, beff, omega_max, larmor_numrev, larmor_thr, mstep, nstep, totalsimtime,&
               therm_fields, do_sc, sc_step, sc_nstep, sc_tidx, sd_phaseflag, delta_t)
         end if
      end if

   end subroutine adapt_time_step

   !> Calculate new timestep
   subroutine calculate_timestep(Natom, Mensemble, beff, omega_max, larmor_numrev, larmor_thr, mstep, nstep, totalsimtime, &
         therm_fields, do_sc, sc_step, sc_nstep, sc_tidx, sd_phaseflag,delta_t)

      implicit none

      integer, intent(in) :: Natom                                            !< Number of atoms in the system
      integer, intent(in) :: Mensemble                                        !< Number of ensembles
      integer, intent(in) :: larmor_numrev                                    !< Number of time steps per revolution needed to suppress numerical damping
      integer, intent(in) :: mstep                                     !< Simulation parameters from 0sd from the MP
      character(len=1), intent(in) :: do_sc
      integer, intent(in) :: sc_tidx                                          !< Current sampling time indicator
      integer, intent(inout) :: sc_step                                       !< Sampling period between spin correlation measurements
      integer, intent(inout) :: sc_nstep                                      !< Sampling period between spin correlation measurements
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

      pi = 4.*atan(1.0_dblprec)

      old_omega_max=omega_max
      !$omp parallel do default(shared) schedule(static) private(ik,i,k)
      do ik = 1,Natom*Mensemble
         i=mod(ik-1,Natom)+1
         k=int((ik-1)/Natom)+1
         beff_abs(i,k) = norm2(beff(:,i,k))
         bthm_abs(i,k) = norm2(therm_fields(:,i,k))
      end do
      !$omp end parallel do

      ! Compute the maximum Larmor frequency and the relative change
      omega_max = maxval(beff_abs) + maxval(bthm_abs)
      relprec = abs((old_omega_max-omega_max)/omega_max)

      ! If the relative change in Larmor precession exceeds the threshold value, then
      ! update the time step with respect to the Larmor frequency and the number of time steps per rev. specified in the input handler
      ! Moreover, one needs to adjust the total number of measurement steps that the simulation is running over, i.e., nstep
      if(relprec > larmor_thr) then
         delta_t_old = delta_t
         delta_t = 2*pi/(larmor_numrev*omega_max)
         frac = delta_t/delta_t_old
         int_frac = int(log(frac)/log(2.0_dblprec))   ! Set the fraction as a power of two
         timestep_frac = (int_frac-1)*100
         delta_t = int_frac*delta_t_old

         t_frac = totalsimtime/delta_t
         if(t_frac.GT.huge(nstep)) then
            stop
         end if

         nstep = nint(t_frac)

         ! Print output to terminal
         write(*,'(2x,a)') 'Performs correction of the time step.'
         write(*,'(2x,a,f10.2)') 'Now using a ', timestep_frac, '% larger time step than before.'
         if(nstep.EQ.0) then
            write(*,'(2x,a)') 'Attention: The new time step fell below the total simulation time for the phase being run. This will terminate the simulation immediately.'
         end if
         write(*,'(2x,a,i8,a,i8,a)') 'Current progress: ', mstep,'out of', nstep, 'steps completed.'

         !! Perform correction of sampling correlation frequencies, if used, and while in the sd measurement phase only
         if(do_sc=='C'.and.(.not.sd_phaseflag)) then
            sc_fact = int_frac ! Assuming the same fractional change between samples as the fractional change in delta t
            !sc%sc_step = ceiling(sc%sc_step/sc_fact)+1 ! Round off to closest upper integer and add 1, hence sc_step >= 2 (see 0sd.f90, line~1194)
            write(*,'(2x,a,i8)') 'New spin correlation sampling period: ', sc_step
            !call calc_corr_w(delta_t,sc) ! Perform frequency correction
         end if
      end if

   end subroutine calculate_timestep

   !> Maximum angular frequency
   subroutine calculate_omegainit(omega_max, larmor_numrev, delta_t)
      implicit none

      real(dblprec), intent(out) :: omega_max                          !< Previously saved maximum angular frequency
      integer, intent(in) :: larmor_numrev                             !< Number of time steps per revolution needed to suppress numerical damping
      real(dblprec), intent(in) :: delta_t                             !< Time step of the simulation

      real(dblprec) :: pi

      pi = 4.*atan(1.0_dblprec)

      ! The initial Larmor precession depends on the input specified time step
      omega_max = 2*pi/(larmor_numrev*delta_t)

   end subroutine calculate_omegainit

end module adaptivetimestepping
