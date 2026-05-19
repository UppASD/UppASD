!------------CHELPER-------------!
!> @author
!> Niklas Fejes
!>
!> @brief
!> Helper for c code
!--------------------------------!
module Chelper
   use iso_c_binding, only : c_double
   use InputData
   use SimulationData,   only : rstep, lambda1
   use MomentData,       only : emomM, emom, mmom, mmom0, mmom2, emom2, mmomi
   use ChemicalData,     only : asite_ch, achem_ch,atype_ch
   use AutoCorrelation,  only : nspinwait, spinwait
   use MicroWaveField,   only : mwffield
   use Constants,        only : gama, mub, k_bolt
   use HamiltonianData,  only : ham

   use prn_averages,     only : calc_and_print_cumulant, do_avrg, mavg, binderc, avrg_step, do_cumu, cumu_step, cumu_buff
   use prn_trajectories, only : do_tottraj, ntraj, tottraj_buff,tottraj_step, traj_step
   use Temperature,      only : temp_array, ipTemp_array
   use Spinicedata,      only : vert_ice_coord
   use Fielddata,        only : thermal_field, beff, beff1, beff3,  b2eff, external_field
   use Gradients,        only : dxyz_vec, dxyz_atom, dxyz_list
   use Systemdata,       only : coord, atype

   use Measurements,     only : measure, do_measurements, flush_measurements, calc_mavrg
   use UpdateMoments,    only : moment_update

   use Correlation
   use Correlation_core
   use Correlation_print, only : print_gk, print_gkw, print_gkt
   use AutoCorrelation,       only : autocorr_sample, do_autocorr
   use prn_cudameasurements,  only : print_observable
   use ChemicalData, only : achtype
   use MetaTypes
   use Omegas

   implicit none


   private

   public :: fortran_do_measurements,fortran_measure,fortran_measure_moment,        &
      fortran_moment_update,fortran_flush_measurements,FortranData_Initiate,        &
      fortran_calc_simulation_status_variables, fortran_print_correlations,          &
      fortran_measure_correlations, fortran_measure_rest, fortran_print_measurables

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
   !> Calculates and returns the magnetic average to C++ simulation.
   !> Needed so that the right mavrg and bindc is printed during simulation.
   !> Binderc is used as a pointer and is only updated if not already calculated
   !> elsewhere.
   !
   !> @author
   !> Thomas Nystrand
   !---------------------------------------------------------------------
   subroutine fortran_calc_simulation_status_variables(mavrg) bind(C,name='fortran_calc_simulation_status_variables')
      implicit none
      real(c_double), intent(inout) :: mavrg
      real(dblprec) :: mavrg_loc

      mavrg_loc = real(mavrg, kind=dblprec)
      call calc_mavrg(Natom, Mensemble, emomM, mavrg_loc)
      if(do_cumu=='N') then
         call calc_and_print_cumulant(Natom,Mensemble,emomM,simid,Temp,1.0_dblprec, &
            0.0_dblprec,plotenergy,cumu_buff,.false.)
      endif
      mavrg = real(mavrg_loc, kind=c_double)
   end subroutine fortran_calc_simulation_status_variables

   ! Measurements with pre-set parameters
   subroutine fortran_measure(cmstep)
      implicit none
      integer, intent(in) :: cmstep !< Current simulation step

      integer :: cgk_flag
      cgk_flag=0

      call measure(Natom,Mensemble,NT,NA,nHam,N1,N2,N3,simid,cmstep,emom,emomM,mmom,&
         Nchmax,do_ralloy,Natom_full,asite_ch,achem_ch,atype,plotenergy,Temp,       &
         1.0_dblprec,0.0_dblprec,real_time_measure,delta_t,logsamp,ham%max_no_neigh,ham%nlist,  &
         ham%ncoup,ham%nlistsize,ham%aham,thermal_field,beff,beff1,beff3,coord,     &
         ham%ind_list_full,ham%ind_nlistsize,ham%ind_nlist,ham%max_no_neigh_ind,    &
         ham%sus_ind,do_mom_legacy,mode)

      ! Spin correlation
      ! Sample magnetic moments for correlation functions
      call correlation_wrapper(Natom,Mensemble,coord,simid,emomM,cmstep,delta_t,  &
      NT_meta,atype_meta,Nchmax,achtype,sc,do_sc,do_sr,cgk_flag)

   end subroutine fortran_measure


   ! Measurements with pre-set parameters
   subroutine fortran_measure_moment(ext_emomM, ext_emom, ext_mmom, ext_mstep)
      implicit none
      real(dblprec), dimension(3,Natom, Mensemble), intent(in) :: ext_emom
      real(dblprec), dimension(3,Natom, Mensemble), intent(in) :: ext_emomM
      real(dblprec), dimension(Natom, Mensemble), intent(in)   :: ext_mmom
      integer, intent(in) :: ext_mstep

      integer :: cgk_flag
      cgk_flag=0

      call measure(Natom,Mensemble,NT,NA,nHam,N1,N2,N3,simid,ext_mstep,ext_emom,    &
         ext_emomM,ext_mmom,Nchmax,do_ralloy,Natom_full,asite_ch,achem_ch,atype,    &
         plotenergy,Temp,1.0_dblprec,0.0_dblprec,real_time_measure,delta_t,logsamp,             &
         ham%max_no_neigh,ham%nlist,ham%ncoup,ham%nlistsize,ham%aham,thermal_field, &
         beff,beff1,beff3,coord,ham%ind_list_full,ham%ind_nlistsize,ham%ind_nlist,  &
         ham%max_no_neigh_ind,ham%sus_ind,do_mom_legacy,mode)


      ! Spin correlation
      ! Sample magnetic moments for correlation functions
         call correlation_wrapper(Natom,Mensemble,coord,simid,emomM,ext_mstep,delta_t,  &
         NT_meta,atype_meta,Nchmax,achtype,sc,do_sc,do_sr,cgk_flag)

   end subroutine fortran_measure_moment

   ! Correlation-only sampling entry used by newer GPU measurement pipeline.
   subroutine fortran_measure_correlations(ext_emomM, ext_emom, ext_mmom, ext_mstep)
      implicit none
      real(dblprec), dimension(3,Natom, Mensemble), intent(in) :: ext_emom
      real(dblprec), dimension(3,Natom, Mensemble), intent(in) :: ext_emomM
      real(dblprec), dimension(Natom, Mensemble), intent(in)   :: ext_mmom
      integer, intent(in) :: ext_mstep
      integer :: cgk_flag

      cgk_flag=0
      call correlation_wrapper(Natom,Mensemble,coord,simid,ext_emomM,ext_mstep,delta_t,  &
         NT_meta,atype_meta,Nchmax,achtype,sc,do_sc,do_sr,cgk_flag)
   end subroutine fortran_measure_correlations

   ! Measurement entry for quantities that remain on the Fortran side in GPU runs.
   subroutine fortran_measure_rest(ext_emomM, ext_emom, ext_mmom, ext_beff, ext_mstep)
      implicit none
      real(dblprec), dimension(3,Natom, Mensemble), intent(in) :: ext_emom
      real(dblprec), dimension(3,Natom, Mensemble), intent(in) :: ext_emomM
      real(dblprec), dimension(Natom, Mensemble), intent(in)   :: ext_mmom
      real(dblprec), dimension(3,Natom, Mensemble), intent(inout) :: ext_beff
      integer, intent(in) :: ext_mstep

      call measure(Natom,Mensemble,NT,NA,nHam,N1,N2,N3,simid,ext_mstep,ext_emom,    &
         ext_emomM,ext_mmom,Nchmax,do_ralloy,Natom_full,asite_ch,achem_ch,atype,    &
         plotenergy,Temp,1.0_dblprec,0.0_dblprec,real_time_measure,delta_t,logsamp, &
         ham%max_no_neigh,ham%nlist,ham%ncoup,ham%nlistsize,ham%aham,thermal_field, &
         ext_beff,beff1,beff3,coord,ham%ind_list_full,ham%ind_nlistsize,             &
         ham%ind_nlist,ham%max_no_neigh_ind,ham%sus_ind,do_mom_legacy,mode)
   end subroutine fortran_measure_rest

   ! Print GPU-accumulated correlations through existing Fortran printers.
   subroutine fortran_print_correlations()
      implicit none
      if(do_sc=='C'.or.do_sc=='Y') then
         call print_gk(NT_meta, Nchmax, sc, sc, simid, sc%label)
      endif
      if(do_sc=='Q'.or.do_sc=='Y') then
         call print_gkw(NT_meta, Nchmax, sc, sc, simid, sc%label)
      endif
      if(do_sc=='T'.or.do_sc=='Y') then
         call print_gkt(NT_meta, Nchmax, sc, sc, simid, sc%label)
      endif
   end subroutine fortran_print_correlations

   ! Print measurables calculated in GPU code via Fortran-side printer.
   subroutine fortran_print_measurables(obs_step, obs_buff, indxb_obs, obs_name, &
        obs_label, obs_dim, obs_buffer, mstep)
      implicit none
      integer, intent(in) :: obs_step, obs_buff, obs_dim
      real(dblprec), dimension(:), allocatable, intent(in) :: indxb_obs
      real(dblprec), dimension(:,:,:), allocatable, intent(in) :: obs_buffer
      character(len=16), intent(in) :: obs_name
      character(len=16), dimension(:), allocatable, intent(in) :: obs_label
      integer, intent(in) :: mstep

      call print_observable(simid, Mensemble, obs_name, obs_step, obs_buff, &
         obs_dim, indxb_obs, obs_buffer, obs_label, real_time_measure, delta_t, mstep)
   end subroutine fortran_print_measurables


   ! Do measurements with pre-set parameters
   subroutine fortran_do_measurements(cmstep, do_copy)
      implicit none
      integer, intent(in) :: cmstep !< Current simulation step
      integer, intent(out) :: do_copy !< Flag if copy or not

      call do_measurements(cmstep,do_avrg,do_tottraj,avrg_step,ntraj,tottraj_step,  &
         traj_step,do_cumu,cumu_step,logsamp,do_copy,do_gpu_measurements)
   end subroutine fortran_do_measurements



   ! Moment update with pre-set parameters
   subroutine fortran_moment_update()
      implicit none
      call moment_update(Natom,Mensemble,mmom,mmom0,mmom2,emom,emom2,emomM,mmomi,   &
         mompar,initexc)
   end subroutine fortran_moment_update



   ! Flush measurements with pre-set parameters
   subroutine fortran_flush_measurements(cmstep)
      implicit none
      integer, intent(in) :: cmstep !< Current simulation step
      call flush_measurements(Natom,Mensemble,NT,NA,N1,N2,N3,simid,cmstep,emom,mmom,&
         Nchmax,atype,real_time_measure,mcnstep,ham%ind_list_full,do_mom_legacy,mode)
   end subroutine fortran_flush_measurements



   ! Initiate pointers for C/C++ implementation
   !> Calls functions in fortrandata.cpp
   subroutine FortranData_Initiate(stt,btorque)
      implicit none
      character(len=1), intent(in) :: STT       !< Treat spin transfer torque? (Y/N)
      real(dblprec), dimension(3,Natom, Mensemble), intent(inout) :: btorque !< Field from (m x dm/dr)

      call FortranData_setConstants(stt,SDEalgh,rstep,nstep,Natom,Mensemble,        &
         ham%max_no_neigh,delta_t,gama,k_bolt,mub,mplambda1,binderc,mavg,mompar,    &
         initexc,ham_inp%do_dm,ham%max_no_dmneigh, ham_inp%do_jtensor, ham_inp%do_anisotropy, nHam)

      !call FortranData_setMatrices(ham%ncoup(1,1,1),ham%nlist(1,1),ham%nlistsize(1),&
      !   beff(1,1,1),b2eff(1,1,1),emomM(1,1,1),emom(1,1,1),emom2(1,1,1),            &
      !   external_field(1,1,1),mmom(1,1),btorque(1,1,1),Temp_array(1),mmom0(1,1),   &
      !   mmom2(1,1),mmomi(1,1),ham%dm_vect(1,1,1),ham%dmlist(1,1),ham%dmlistsize(1))


      call FortranData_setHamiltonian(ham%ncoup,ham%nlist,ham%nlistsize, &
         ham%dm_vect,ham%dmlist,ham%dmlistsize, &
         ham%kaniso, ham%eaniso, ham%taniso, ham%sb, &
         ham%j_tens, ham%aHam, &
         external_field, btorque,Temp_array, &
         ipTemp, ipmcnstep, ipTemp_array, ipnstep, ipdelta_t, iplambda1)

      call FortranData_setLattice(beff, b2eff, emomM, emom, emom2, mmom, mmom0, mmom2, mmomi, &
         dxyz_vec, dxyz_atom, dxyz_list)

      call FortranData_setInputData(gpu_mode, gpu_rng, gpu_rng_seed)

   end subroutine FortranData_Initiate

end module Chelper
