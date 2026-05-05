!------------CHELPER-------------!
!> @author
!> Niklas Fejes
!>
!> @brief
!> Helper for c code
!--------------------------------!
module Chelper
   use iso_c_binding
   use InputData
   use SimulationData,   only : rstep, lambda1
   use MomentData,       only : emomM, emom, mmom, mmom0, mmom2, emom2, mmomi
   use ChemicalData,     only : asite_ch, achem_ch,atype_ch
   !use AutoCorrelation
   use AutoCorrelation
   use MicroWaveField,   only : mwffield
   use Constants,        only : gama, mub, k_bolt, mry
   use HamiltonianData,  only : ham

   use prn_averages,  
   use prn_topology,     only : skyno, skyno_step, skyno_buff
   use Gradients,        only : dxyz_vec, dxyz_atom, dxyz_list
   use Energy,           only : eavg_buff, eavg2_buff, eavg4_buff, eavrg_step, eavrg_buff, calc_energy
   use prn_trajectories, only : do_tottraj, ntraj, tottraj_buff, tottraj_step, &
        traj_step, traj_buff, traj_atom, mmomb, mmomb_traj, emomb, emomb_traj
   !use Temperature,      only : temp_array, iptemp_array
   use Temperature
   use Spinicedata,      only : vert_ice_coord
   use Fielddata,        only : thermal_field, beff, beff1, beff3,  b2eff, external_field, time_external_field
   use Systemdata,       only : coord, atype

   use Measurements,     only : measure, do_measurements, flush_measurements, calc_mavrg
   use UpdateMoments,    only : moment_update

   use Correlation
   use Correlation_core
   use Correlation_Print
   use Correlation_type
   use Correlation_utils, only: find_rmid
   use Omegas
   use Qvectors
   use Correlation_utils
   use AutoCorrelation,  only : autocorr_sample, do_autocorr, spinwait, autocorr_buff, indxb_ac
   use ChemicalData, only : achtype
   use MetaTypes
   use macrocells
   use OptimizationRoutines


   use prn_cudameasurements,   only :  print_observable, print_trajectory

   implicit none

   interface
         subroutine FortranData_setCorrelations(q, r_mid, coord, w, m_k, m_kw, m_kt, deltat_corr, scstep_arr, sc_nsamp, sc_tidx,&
               atype, achtype, m_k_proj, m_k_projch, m_kt_proj, m_kt_projch, m_kw_proj, m_kw_projch) &
                bind(C, name="fortrandata_setcorrelations_")
         import :: c_double, c_int
         real(c_double)    :: q(*)
         real(c_double)    :: r_mid(*)
         real(c_double)    :: coord(*)
         real(c_double)    :: w(*)
         complex(c_double) :: m_k(*)
         complex(c_double) :: m_kw(*)
         complex(c_double) :: m_kt(*)
         complex(c_double) :: m_k_proj(*)
         complex(c_double) :: m_kw_proj(*)
         complex(c_double) :: m_kt_proj(*)
         complex(c_double) :: m_k_projch(*)
         complex(c_double) :: m_kw_projch(*)
         complex(c_double) :: m_kt_projch(*)
         real(c_double)    :: deltat_corr(*)
         real(c_double)    :: scstep_arr(*)
         integer(c_int), intent(inout) :: sc_nsamp
         integer(c_int), intent(inout) :: sc_tidx
         integer(c_int), intent(inout) :: atype(*)
         integer(c_int), intent(inout) :: achtype(*)
      end subroutine FortranData_setCorrelations
   end interface


   private

   public :: fortran_do_measurements,fortran_measure,fortran_measure_moment,        &
      fortran_moment_update,fortran_flush_measurements,FortranData_Initiate,        &
      fortran_calc_simulation_status_variables, fortran_print_measurables,          &
      fortran_print_correlations, fortran_measure_correlations,                     &
      fortran_measure_rest

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
      real(dblprec), intent(inout) :: mavrg
      call calc_mavrg(Natom, Mensemble, emomM, mavrg)
      if(do_cumu=='N') then
         call calc_and_print_cumulant(Natom,Mensemble,emomM,simid,Temp,1.0_dblprec, &
            0.0_dblprec,plotenergy,cumu_buff,.false.)
      endif
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
      !call correlation_wrapper(Natom,Mensemble,coord,simid,emomM,cmstep,delta_t,  &
      !NT_meta,atype_meta,Nchmax,achtype,sc,do_sc,do_sr,cgk_flag)

   end subroutine fortran_measure

      ! Correlations on Fortran side
   subroutine fortran_measure_correlations(ext_emomM, ext_emom, ext_mmom, ext_mstep)
      implicit none
      real(dblprec), dimension(3,Natom, Mensemble), intent(in) :: ext_emom
      real(dblprec), dimension(3,Natom, Mensemble), intent(in) :: ext_emomM
      real(dblprec), dimension(Natom, Mensemble), intent(in)   :: ext_mmom
      integer, intent(in) :: ext_mstep

      integer :: cgk_flag
      cgk_flag=0

      ! Spin correlation
      ! Sample magnetic moments for correlation functions
         call correlation_wrapper(Natom,Mensemble,coord,simid,emomM,ext_mstep,delta_t,  &
         NT_meta,atype_meta,Nchmax,achtype,sc,do_sc,do_sr,cgk_flag)

   end subroutine fortran_measure_correlations


   ! Measurements with pre-set parameters
   subroutine fortran_measure_moment(ext_emomM, ext_emom, ext_mmom, ext_mstep)
      implicit none
      real(dblprec), dimension(3,Natom, Mensemble), intent(in) :: ext_emom
      real(dblprec), dimension(3,Natom, Mensemble), intent(in) :: ext_emomM
      real(dblprec), dimension(Natom, Mensemble), intent(in)   :: ext_mmom
      integer, intent(in) :: ext_mstep
      real(dblprec) ::  totene, totenergy 

      integer :: cgk_flag
      cgk_flag=0
      
      call measure(Natom,Mensemble,NT,NA,nHam,N1,N2,N3,simid,ext_mstep,ext_emom,    &
         ext_emomM,ext_mmom,Nchmax,do_ralloy,Natom_full,asite_ch,achem_ch,atype,    &
         plotenergy,Temp,1.0_dblprec,0.0_dblprec,real_time_measure,delta_t,logsamp,             &
         ham%max_no_neigh,ham%nlist,ham%ncoup,ham%nlistsize,ham%aham,thermal_field, &
         beff,beff1,beff3,coord,ham%ind_list_full,ham%ind_nlistsize,ham%ind_nlist,  &
         ham%max_no_neigh_ind,ham%sus_ind,do_mom_legacy,mode)

      ! Calculate total and term resolved energies
         if(plotenergy>0) then
            if (mod(ext_mstep-1,ene_step)==0) then
               call calc_energy(nHam,ext_mstep,Natom,Nchmax, &
                  conf_num,Mensemble,Natom,Num_macro,1,         &
                  plotenergy,Temp,delta_t,do_lsf,        &
                  lsf_field,lsf_interpolate,real_time_measure,simid,cell_index,            &
                  macro_nlistsize,ext_mmom,ext_emom,ext_emomM,emomM_macro,external_field,              &
                  time_external_field,max_no_constellations,maxNoConstl,                   &
                  unitCellType,constlNCoup,constellations,OPT_flag,                        &
                  constellationsNeighType,totene,NA,N1,N2,N3)
            end if
         endif


      ! Spin correlation
      ! Sample magnetic moments for correlation functions
      !   call correlation_wrapper(Natom,Mensemble,coord,simid,emomM,ext_mstep,delta_t,  &
      !   NT_meta,atype_meta,Nchmax,achtype,sc,do_sc,do_sr,cgk_flag)

   end subroutine fortran_measure_moment

      ! Measurements not implemeted or not planned to be implementedon GPU
   subroutine fortran_measure_rest(ext_emomM, ext_emom, ext_mmom, ext_beff, ext_mstep)
      use Math_functions, only : f_logstep
      use prn_trajectories, only : print_trajectories
      use Polarization,     only : print_pol
      use prn_fields,       only : print_fields
      use Temperature

      implicit none
      real(dblprec), dimension(3,Natom, Mensemble), intent(in) :: ext_emom
      real(dblprec), dimension(Natom, Mensemble), intent(in) :: ext_mmom
      real(dblprec), dimension(3,Natom, Mensemble), intent(in) :: ext_emomM
      real(dblprec), dimension(3,Natom, Mensemble), intent(in) :: ext_beff
      integer, intent(in) :: ext_mstep

      integer :: sstep

      sstep = f_logstep(ext_mstep,logsamp)

      ! Sample spin temperature
      if (do_spintemp=='Y') then
         if(mod(ext_mstep,spintemp_step)==0) then
            call spintemperature(Natom,Mensemble,ext_mstep,1,simid,ext_emomM,ext_beff,1)
         end if
      endif

      call print_trajectories(Natom,sstep,ext_mstep,Mensemble,ext_emom,ext_mmom,delta_t,        &
         real_time_measure,simid,do_mom_legacy,mode)

      call print_pol(sstep,ext_mstep,Natom,Mensemble,ham%max_no_neigh,ham%nlist,ham%nlistsize,ext_emomM,&
         delta_t,simid,real_time_measure)

      !if (do_autocorr=='Y') then
      !   if ( mod(sstep-1,ac_step)==0.or.any(spinwaitt==sstep-1)) then
      !      call buffer_autocorr(Natom,Mensemble,do_autocorr,ext_mstep,nspinwait,spinwait,ext_emom,&
      !         ext_emomM,bcount_ac,delta_t,real_time_measure)
      !      if (bcount_ac==ac_buff) then
      !         ! Write the autocorrelation buffer to file
      !         call prn_autocorr(Natom,simid,do_autocorr,nspinwait,real_time_measure)
      !         ! Reset statistics buffer
      !         bcount_ac=1
      !      else
      !         bcount_ac=bcount_ac+1
      !      endif
      !   endif
      !end if

      !call print_fields(ext_mstep,sstep,Natom,Mensemble,simid,real_time_measure,delta_t,&
      !   ext_beff,thermal_field,beff1,beff3,ext_emom)

   end subroutine fortran_measure_rest


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
      integer, intent(in) :: cmstep !< Current simulation stepfind_rmid(rmid,coord,Natom)
      call flush_measurements(Natom,Mensemble,NT,NA,N1,N2,N3,simid,cmstep,emom,mmom,&
         Nchmax,atype,real_time_measure,mcnstep,ham%ind_list_full,do_mom_legacy,mode)
   end subroutine fortran_flush_measurements

      ! print GPU calculated correlations
   subroutine fortran_print_correlations()
      implicit none
      integer ::cgk_flag
      cgk_flag = 2
      !type(corr_t), intent(inout) :: cc !< Derived type for correlation data
      if(do_sc=='C'.or.do_sc=='Y') then
         call print_gk(NT, Nchmax, sc, sc, simid, sc%label)
      endif
     if(do_sc=='Q'.or.do_sc=='Y') then
         call print_gkw(NT, Nchmax, sc, sc, simid, sc%label)
      endif
     if(do_sc=='T'.or.do_sc=='Y') then
         call print_gkt(NT, Nchmax, sc, sc, simid, sc%label)

      endif
   end subroutine fortran_print_correlations



   ! Initiate pointers for C/C++ implementation
   !> Calls functions in fortrandata.cpp
   subroutine FortranData_Initiate(stt,btorque,cc, phase)
      
      implicit none
      character(len=1), intent(in) :: STT !< Treat spi p_sc_max_nstn transfer torque? (Y/N)
      character(len=1), intent(in) :: phase !< initial or measurement (I/M)
      type(corr_t), intent(inout) :: cc !< Derived type for correlation data
      real(dblprec), dimension(3,Natom, Mensemble), intent(inout) :: btorque !< Field from (m x dm/dr)
      integer :: zeroflag = 0
      
      !!!TODO: replace those with actual variables 
      integer :: ene_step = 10
      integer :: ene_buff = 100
      !character(len=1)::do_projch_avrg = 'Y'
      !character(len=1)::do_cumu_proj = 'Y'


      if(cc%do_proj=='C'.or.cc%do_proj=='Y'.or.cc%do_proj=='T'.or.cc%do_proj=='Q'.or.cc%do_projch=='C'.or.cc%do_projch=='Y'.or.cc%do_projch=='Q'.or.cc%do_projch=='T') then
         print *, "Projections are not available in GPU correlations yet, please use do_gpu_correlations 0"
         return  
      end if

      if(do_gpu_correlations=='Y'.and.phase=='M') then
         ! Calculate r_mid for GPU correlations, as it is needed for the correlation calculations and not calculated on the Fortran side otherwise

            print *,'AB first test, do_sc = ', do_sc
         if (do_sc=='Y' .or. do_sc=='C') then
            ! Initializing Fortran-side correlation arrays. Here for static correlation function
            cc%gk_flag=2
            allocate(cc%m_k(3,nq))
            cc%m_k=0.0_dblprec
            ! call find_rmid(r_mid,coord,Natom)
            ! Possible alternative: Call the Fortran routine with init flag = 0
            ! zeroflag = 0
            ! call calc_gk2(Natom, Mensemble, NT,atype,Nchmax,achtype, cc, coord, simid, emomM, zeroflag)
         end if
         if (do_sc=='Y' .or. do_sc=='Q') then
            ! Initializing Fortran-side correlation arrays. Here for dynamic correlation function
            cc%gkt_flag=2
            allocate(cc%m_kt(3,nq,cc%sc_max_nstep))
            cc%m_kt=0.0_dblprec
            call allocate_deltatcorr(.true.,cc)
            allocate(cc%m_kw(3,nq,cc%nw))
            cc%m_kw=0.0_dblprec
            ! Possible alternative: Call the Fortran routine with init flag = 0
            ! zeroflag = 0
            ! cc%gkt_flag=0
            ! call calc_gkt(Natom, Mensemble, NT,atype,Nchmax,achtype, cc, coord, emomM, zeroflag)
         end if
         call find_rmid(r_mid,coord,Natom)

      end if
     
      print *,' AB shape of m_k', shape(cc%m_k)
      print *,' AB allocated?', allocated(cc%m_k)

      !print *, 'FORTRAN EMOOOOM', size(emomM), size(emomM,1), size(emomM,2)



      call FortranData_setFlags(ham_inp%do_dm, ham_inp%do_jtensor, ham_inp%do_anisotropy, &
           do_avrg, do_proj_avrg, do_projch_avrg, do_cumu, do_cumu_proj, plotenergy, do_autocorr, do_tottraj, ntraj, &
           do_gpu_measurements, skyno, do_sc, do_gpu_correlations, real_time_measure, &
           cc%do_proj, cc%do_projch, do_ralloy)

      call FortranData_setConstants(stt,SDEalgh,rstep,nstep,Natom,Mensemble, &
         ham%max_no_neigh,delta_t,gama,k_bolt,mub,mplambda1,binderc,mavg,mompar, &
         initexc,ham%max_no_dmneigh,nHam, Temp, ipmcnphase, mcnstep, ipnphase, &
         avrg_step, avrg_buff, cumu_step, cumu_buff, ene_step, ene_buff, &
         tottraj_step, tottraj_buff, skyno_step, skyno_buff, nq, sc_window_fun, &
         cc%nw, cc%sc_sep, cc%sc_step, cc%sc_max_nstep, nspinwait, ac_step, ac_buff, &
         NT, Nchmax, mry, NA, Natom_full)

      call FortranData_setHamiltonian(ham%ncoup,ham%nlist,ham%nlistsize, &
         ham%dm_vect,ham%dmlist,ham%dmlistsize, &
         ham%kaniso, ham%eaniso, ham%taniso, ham%sb, &
         ham%j_tens, ham%aHam, &
         external_field, btorque,Temp_array, &
         ipTemp, ipmcnstep, ipTemp_array, ipnstep)

      call FortranData_setLattice(beff, b2eff, emomM, emom, emom2, mmom, mmom0, mmom2, mmomi, &
         dxyz_vec, dxyz_atom, dxyz_list)


      call FortranData_setMeasurables( &
           mavg_buff, mavg2_buff, mavg4_buff, &
           mavg_buff_proj, mavg2_buff_proj, mavg4_buff_proj, &
           binderc, avrgmcum, avrgm2cum, avrgm4cum, &
           eavg_buff, eavg2_buff, &
           traj_step, traj_buff, traj_atom, &
           mmomb, mmomb_traj, emomb, emomb_traj, &
           spinwaitt, spinwait, &
           indxb_ac, autocorr_buff, &
           achem_ch, asite_ch) !TODO: get rd of those once printed on CPU

      call FortranData_setCorrelations(q, r_mid, coord, cc%w, cc%m_k, cc%m_kw, cc%m_kt, cc%deltat_corr, &
          cc%scstep_arr, cc%sc_nsamp, cc%sc_tidx, atype, achtype, cc%m_k_proj, cc%m_k_projch, &
          cc%m_kt_proj, cc%m_kt_projch, cc%m_kw_proj, cc%m_kw_projch)


      call FortranData_setInputData(gpu_mode, gpu_rng, gpu_rng_seed)

   end subroutine FortranData_Initiate

       ! Print measurables calculated in CUDA
   subroutine fortran_print_measurables(obs_step, obs_buff, indxb_obs, obs_name, &
        obs_label, obs_dim, obs_buffer, mstep)
      implicit none
      integer, intent(in) :: obs_step, obs_buff, obs_dim
      real(dblprec), dimension(:), allocatable, intent(in) :: indxb_obs
      real(dblprec), dimension(:,:,:), allocatable, intent(in) :: obs_buffer
      character(len=16), intent(in) :: obs_name !< Observable name
      character(len=16), dimension(:), allocatable, intent(in) :: obs_label
      integer, intent(in) :: mstep !< Current simulation step
      call print_observable(simid, Mensemble, obs_name, obs_step, obs_buff, &
      obs_dim, indxb_obs, obs_buffer, obs_label, real_time_measure, delta_t, mstep)
   end subroutine fortran_print_measurables


end module Chelper


!   subroutine fortran_print_measurables(obs_step, obs_buff, indxb_obs, obs_name, &
!        obs_label, obs_dim, obs_buffer, mstep)
!      implicit none
!      integer, intent(in) :: obs_step, obs_buff, obs_dim
!      real(dblprec), dimension(obs_buff), intent(in) :: indxb_obs
!      real(dblprec), dimension(obs_dim, Natom, Mensemble), intent(in) :: obs_buffer
!      character(len=16), intent(in) :: obs_name !< Observable name
!      character(len=16), dimension(obs_dim) :: obs_label
!      integer, intent(in) :: mstep !< Current simulation step
!      call print_observable(simid, Mensemble, obs_name, obs_step, obs_buff, &
!      obs_dim, indxb_obs, obs_buffer, obs_label, real_time_measure, delta_t, mstep)
!   end subroutine fortran_print_measurables
