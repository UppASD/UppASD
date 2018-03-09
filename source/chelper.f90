!------------CHELPER-------------!
!> @author
!> Niklas Fejes
!>
!> @brief
!> Helper for c code
!--------------------------------!
module Chelper
   use InputData
   use SimulationData,   only : rstep, lambda1
   use MomentData,       only : emomM, emom, mmom, mmom0, mmom2, emom2, mmomi
   use ChemicalData,     only : asite_ch, achem_ch,atype_ch
   use AutoCorrelation,  only : nspinwait, spinwait
   use MicroWaveField,   only : mwffield
   use Constants,        only : gama, mub, k_bolt
   use HamiltonianData,  only : ncoup, nlist, nlistsize, max_no_neigh, &
      dmlist, dm_vect, dmlistsize, max_no_dmneigh,ind_nlist,ind_nlistsize

   use prn_averages,     only : calc_and_print_cumulant, do_avrg, mavg, binderc,  avrg_step, do_cumu, cumu_step, cumu_buff
   use prn_trajectories, only : do_tottraj, ntraj, tottraj_buff,tottraj_step, traj_step
   use Temperature,      only : temp_array
   use Spinicedata,      only : vert_ice_coord
   use Fielddata,        only : thermal_field, beff, beff1, beff3,  b2eff, external_field
   use Systemdata,       only : coord, atype

   use Measurements,     only : measure, do_measurements, flush_measurements, calc_mavrg
   use UpdateMoments,    only : moment_update


   implicit none


   private

   public :: fortran_do_measurements, fortran_measure, fortran_measure_moment, fortran_moment_update, &
      fortran_flush_measurements, FortranData_Initiate, fortran_calc_simulation_status_variables

contains

   subroutine array_test(A,B,arr)
      implicit none
      integer, intent(in) :: A
      integer, intent(in) :: B
      real(dblprec), dimension(3,A,B), intent(out) :: arr
      integer :: i,j,c

      c=0
      do i=1,A
         do j=1,B
            arr(1,i,j) = c
            c=c+1
         end do
      end do

   end subroutine array_test

   !---------------------------------------------------------------------
   !> @brief
   !> Calculates and returns the magnetic avergage to C++ simulation.
   !> Needed so that the right mavrg and bindc is printed during simulation.
   !> Binderc is used as a pointer and is only updated if not already calculated
   !> elsewhere.
   !
   !> @author
   !> Thomas Nystrand
   !---------------------------------------------------------------------
   subroutine fortran_calc_simulation_status_variables(mavrg)
      implicit none
      real(dblprec), intent(inout) :: mavrg
      call calc_mavrg(Natom, Mensemble, emomM, mavrg)
      if(do_cumu=='N') then
         call calc_and_print_cumulant(Natom, Mensemble, emomM, simid, Temp, plotenergy, cumu_buff, .false.)
      endif
   end subroutine fortran_calc_simulation_status_variables

   ! Measurements with pre-set parameters
   subroutine fortran_measure(cmstep)
      implicit none
      integer, intent(in) :: cmstep !< Current simulation step

      call measure(Natom, Mensemble, NT, NA, N1, N2, N3, simid, &
         cmstep, &
         emom, emomM, mmom, & ! in MomentData
         Nchmax, do_ralloy, Natom_full,&
         asite_ch, achem_ch, atype, & ! in ChemicalData, SystemData
         plotenergy, Temp, &
         real_time_measure, delta_t, logsamp, max_no_neigh, nlist ,ncoup ,nlistsize ,&
         thermal_field, beff,beff1,beff3,coord,ind_mom,ind_nlistsize,ind_nlist,atype_ch)
   end subroutine fortran_measure


   ! Measurements with pre-set parameters
   subroutine fortran_measure_moment(ext_emomM, ext_emom, ext_mmom, ext_mstep)
      implicit none
      real(dblprec), dimension(3,Natom, Mensemble), intent(in) :: ext_emom
      real(dblprec), dimension(3,Natom, Mensemble), intent(in) :: ext_emomM
      real(dblprec), dimension(Natom, Mensemble), intent(in)   :: ext_mmom
      integer, intent(in) :: ext_mstep

      call measure(Natom, Mensemble, NT, NA, N1, N2, N3, simid, &
         ext_mstep, &
         ext_emom, ext_emomM, ext_mmom, &
         Nchmax, do_ralloy, Natom_full,&
         asite_ch, achem_ch, atype, &
         plotenergy, Temp, &
         real_time_measure, delta_t, logsamp, max_no_neigh, nlist ,ncoup ,nlistsize ,&
         thermal_field, beff,beff1,beff3,coord,ind_mom,ind_nlistsize,ind_nlist,atype_ch)
   end subroutine fortran_measure_moment


   ! Do measurements with pre-set parameters
   subroutine fortran_do_measurements(cmstep, do_copy)
      implicit none
      integer, intent(in) :: cmstep !< Current simulation step
      integer, intent(out) :: do_copy !< Flag if copy or not

      call do_measurements(cmstep, do_avrg, do_tottraj, avrg_step, ntraj, tottraj_step, &
         traj_step, do_cumu, cumu_step, logsamp, do_copy)
   end subroutine fortran_do_measurements



   ! Moment update with pre-set parameters
   subroutine fortran_moment_update()
      implicit none
      call moment_update(Natom, Mensemble, mmom, mmom0, mmom2, emom, emom2, emomM, mmomi, mompar, initexc)
   end subroutine fortran_moment_update



   ! Flush measurements with pre-set parameters
   subroutine fortran_flush_measurements(cmstep)
      implicit none
      integer, intent(in) :: cmstep !< Current simulation step
      call flush_measurements(Natom, Mensemble, NT, NA, N1, N2, N3, simid, cmstep,emom, mmom, &
         Nchmax,atype,real_time_measure,mcnstep,do_ralloy,Natom_full,atype_ch,ind_mom)
   end subroutine fortran_flush_measurements



   ! Initiate pointers for C/C++ implementation
   !> Calls functions in fortrandata.cpp
   subroutine FortranData_Initiate(stt,btorque)
      implicit none
      character(len=1), intent(in) :: STT       !< Treat spin transfer torque? (Y/N)
      real(dblprec), dimension(3,Natom, Mensemble), intent(inout) :: btorque !< Field from (m x dm/dr)

      call FortranData_setConstants(stt,SDEalgh,rstep,nstep,Natom,Mensemble,max_no_neigh, &
         delta_t,gama,k_bolt,mub,mplambda1, &
         binderc,mavg,mompar,initexc,do_dm,max_no_dmneigh)

      call FortranData_setMatrices(ncoup(1,1,1),nlist(1,1),nlistsize(1),beff(1,1,1),b2eff(1,1,1), &
         emomM(1,1,1),emom(1,1,1),emom2(1,1,1),external_field(1,1,1), &
         mmom(1,1),btorque(1,1,1),Temp_array(1),mmom0(1,1),mmom2(1,1),mmomi(1,1), &
         dm_vect(1,1,1),dmlist(1,1),dmlistsize(1))

      call FortranData_setInputData(gpu_mode, gpu_rng, gpu_rng_seed)

   end subroutine FortranData_Initiate



end module Chelper
