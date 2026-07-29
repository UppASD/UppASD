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
   use SimulationData,   only : rstep, mstep, lambda1
   use MomentData,       only : emomM, emom, mmom, mmom0, mmom2, emom2, mmomi
   use ChemicalData,     only : asite_ch, achem_ch,atype_ch
   !use AutoCorrelation
   use AutoCorrelation
   use MicroWaveField,   only : mwffield
   use Constants,        only : gama, mub, k_bolt, mry
   use HamiltonianData,  only : ham
   use AdaptiveCGProduction, only : adaptive_cg_state, adaptive_cg_is_enabled

   use prn_averages,  only : avrg_buff, avrg_step, avrgm2cum, avrgm4cum, avrgmcum, &
        binderc, calc_and_print_cumulant, cumu_buff, cumu_step, do_avrg, do_cumu, do_cumu_proj, do_proj_avrg, do_projch_avrg, &
        mavg, mavg_buff, mavg2_buff, mavg4_buff, mavg_buff_proj, mavg2_buff_proj, &
        mavg4_buff_proj
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
   use Omegas
   use Qvectors
   use Correlation_utils
   use ChemicalData, only : achtype
   use MetaTypes
   use macrocells


   use prn_cudameasurements,   only :  print_observable, print_trajectory

   implicit none

   interface
      subroutine FortranData_setFlags(do_dm, do_jtensor, do_anisotropy, do_avrg, do_proj_avrg, do_projch_avrg, &
            do_cumu, do_cumu_proj, plotenergy, do_autocorr, do_tottraj, ntraj, do_gpu_measurements, skyno, &
            do_sc, do_gpu_correlations, do_gpu_timings, real_time_measure, do_sc_proj, do_sc_projch, do_ralloy) &
            bind(C, name="fortrandata_setflags_")
         import :: c_int, c_char
         integer(c_int), intent(inout) :: do_dm, do_jtensor, do_anisotropy
         character(c_char), intent(inout) :: do_avrg, do_proj_avrg, do_projch_avrg
         character(c_char), intent(inout) :: do_cumu, do_cumu_proj
         integer(c_int), intent(inout) :: plotenergy
         character(c_char), intent(inout) :: do_autocorr, do_tottraj
         integer(c_int), intent(inout) :: ntraj
         character(c_char), intent(inout) :: do_gpu_measurements, skyno, do_sc, do_gpu_correlations, do_gpu_timings
         character(c_char), intent(inout) :: real_time_measure, do_sc_proj, do_sc_projch
         integer(c_int), intent(inout) :: do_ralloy
      end subroutine FortranData_setFlags

      subroutine FortranData_setConstants(stt, SDEalgh, rstep, nstep, Natom, Mensemble, max_no_neigh, &
            delta_t, gamma, k_bolt, mub, mplambda1, binderc, mavg, mompar, initexc, max_no_dmneigh, nHam, &
            Temp, ipmcnphase, mcnstep, ipnphase, avrg_step, avrg_buff, cumu_step, cumu_buff, ene_step, ene_buff, &
            tottraj_step, tottraj_buff, skyno_step, skyno_buff, nq, sc_window_fun, nw, sc_sep, sc_step, &
            sc_max_nstep, nspinwait, ac_step, ac_buff, nt, Nchmax, mry, NA, Natom_full) &
            bind(C, name="fortrandata_setconstants_")
         import :: c_int, c_char, c_double
         character(c_char), intent(in) :: stt
         character(c_char), intent(inout) :: initexc
         integer(c_int), intent(inout) :: SDEalgh, rstep, nstep, Natom, Mensemble, max_no_neigh
         real(c_double), intent(inout) :: delta_t, gamma, k_bolt, mub, mplambda1, binderc, mavg
         integer(c_int), intent(inout) :: mompar, max_no_dmneigh, nHam
         real(c_double), intent(inout) :: Temp
         integer(c_int), intent(inout) :: ipmcnphase, mcnstep, ipnphase, avrg_step, avrg_buff, cumu_step, cumu_buff
         integer(c_int), intent(inout) :: ene_step, ene_buff, tottraj_step, tottraj_buff, skyno_step, skyno_buff
         integer(c_int), intent(inout) :: nq, sc_window_fun, nw, sc_sep, sc_step, sc_max_nstep, nspinwait
         integer(c_int), intent(inout) :: ac_step, ac_buff, nt, Nchmax, NA, Natom_full
         real(c_double), intent(inout) :: mry
      end subroutine FortranData_setConstants

      subroutine FortranData_setHamiltonian(ncoup, nlist, nlistsize, dmvect, dmlist, dmlistsize, kaniso, &
            eaniso, taniso, sb, j_tensor, aHam, external_field, btorque, Temp_array, ipTemp, ipmcnstep, &
            ipTemp_array, ipnstep, ipdelta_t, iplambda1) bind(C, name="fortrandata_sethamiltonian_")
         import :: c_int, c_double
         real(c_double), intent(inout) :: ncoup(*), dmvect(*), kaniso(*), eaniso(*), sb(*), j_tensor(*)
         integer(c_int), intent(inout) :: nlist(*), nlistsize(*), dmlist(*), dmlistsize(*), taniso(*), aHam(*)
         real(c_double), intent(inout) :: external_field(*), btorque(*), Temp_array(*), ipTemp(*)
         integer(c_int), intent(inout) :: ipmcnstep(*)
         real(c_double), intent(inout) :: ipTemp_array(*), ipdelta_t(*), iplambda1(*)
         integer(c_int), intent(inout) :: ipnstep(*)
      end subroutine FortranData_setHamiltonian

      subroutine FortranData_setLattice(beff, b2eff, emomM, emom, emom2, mmom, mmom0, mmom2, mmomi, &
            dxyz_vec, dxyz_atom, dxyz_list, atype) bind(C, name="fortrandata_setlattice_")
         import :: c_int, c_double
         real(c_double), intent(inout) :: beff(*), b2eff(*), emomM(*), emom(*), emom2(*)
         real(c_double), intent(inout) :: mmom(*), mmom0(*), mmom2(*), mmomi(*), dxyz_vec(*)
         integer(c_int), intent(inout) :: dxyz_atom(*), dxyz_list(*), atype(*)
      end subroutine FortranData_setLattice

      subroutine FortranData_setMeasurables(mavg_buff, mavg2_buff, mavg4_buff, mavg_buff_proj, mavg2_buff_proj, &
            mavg4_buff_proj, binderc, avrgmcum, avrgm2cum, avrgm4cum, eavg_buff, eavg2_buff, traj_step, &
            traj_buff, traj_atom, mmomb, mmomb_traj, emomb, emomb_traj, spinwaitt, spinwait, indxb_ac, &
            autocorr_buff, achem_ch, asite_ch) bind(C, name="fortrandata_setmeasurables_")
         import :: c_int, c_double
         real(c_double), intent(inout) :: mavg_buff(*), mavg2_buff(*), mavg4_buff(*), mavg_buff_proj(*)
         real(c_double), intent(inout) :: mavg2_buff_proj(*), mavg4_buff_proj(*)
         real(c_double), intent(inout) :: binderc, avrgmcum, avrgm2cum, avrgm4cum
         real(c_double), intent(inout) :: eavg_buff(*), eavg2_buff(*), mmomb(*), mmomb_traj(*), emomb(*)
         real(c_double), intent(inout) :: emomb_traj(*), spinwait(*), indxb_ac(*), autocorr_buff(*)
         integer(c_int), intent(inout) :: traj_step(*), traj_buff(*), traj_atom(*), spinwaitt(*), achem_ch(*), asite_ch(*)
      end subroutine FortranData_setMeasurables

      subroutine FortranData_setInputData(gpu_mode, gpu_rng, gpu_rng_seed) &
            bind(C, name="fortrandata_setinputdata_")
         import :: c_int
         integer(c_int), intent(inout) :: gpu_mode, gpu_rng, gpu_rng_seed
      end subroutine FortranData_setInputData

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

      subroutine FortranData_setGpuGeometry(N1, N2, N3, NA, BC1, BC2, BC3, C1, C2, C3, Bas, alat, do_gpu_convolution) &
               bind(C, name="fortrandata_setgpugeometry_")
         import :: c_int, c_char, c_double
         integer(c_int), intent(inout) :: N1
         integer(c_int), intent(inout) :: N2
         integer(c_int), intent(inout) :: N3
         integer(c_int), intent(inout) :: NA
         character(c_char), intent(inout) :: BC1
         character(c_char), intent(inout) :: BC2
         character(c_char), intent(inout) :: BC3
         real(c_double), intent(inout) :: C1(*)
         real(c_double), intent(inout) :: C2(*)
         real(c_double), intent(inout) :: C3(*)
         real(c_double), intent(inout) :: Bas(*)
         real(c_double), intent(inout) :: alat
         character(c_char), intent(inout) :: do_gpu_convolution
      end subroutine FortranData_setGpuGeometry

      subroutine FortranData_setGpuDipole(mode, surface, alpha, rcut, tol, mesh) &
               bind(C, name="fortrandata_setgpudipole_")
         import :: c_int, c_double
         integer(c_int), intent(inout) :: mode, surface
         real(c_double), intent(inout) :: alpha, rcut, tol
         integer(c_int), intent(inout) :: mesh(*)
      end subroutine FortranData_setGpuDipole

      subroutine FortranData_setMacrocell(do_dip, Num_macro, block_x, block_y, block_z, cell_index, macro_nlistsize, &
            macro_center, macro_min_coord, macro_max_coord) &
               bind(C, name="fortrandata_setmacrocell_")
         import :: c_int, c_double
         integer(c_int), intent(inout) :: do_dip, Num_macro, block_x, block_y, block_z
         integer(c_int), intent(inout) :: cell_index(*), macro_nlistsize(*)
         real(c_double), intent(inout) :: macro_center(*), macro_min_coord(*), macro_max_coord(*)
      end subroutine FortranData_setMacrocell

      subroutine FortranData_clearMacrocell() bind(C, name="fortrandata_clearmacell_")
      end subroutine FortranData_clearMacrocell

      subroutine FortranData_setPmeMacrocell(Num_macro, macro_grid, cell_index, macro_nlistsize, &
            macro_center, macro_min_coord, macro_max_coord) bind(C, name="fortrandata_setpmemacrocell_")
         import :: c_int, c_double
         integer(c_int), intent(inout) :: Num_macro, macro_grid(*)
         integer(c_int), intent(inout) :: cell_index(*), macro_nlistsize(*)
         real(c_double), intent(inout) :: macro_center(*), macro_min_coord(*), macro_max_coord(*)
      end subroutine FortranData_setPmeMacrocell

      subroutine FortranData_clearPmeMacrocell() bind(C, name="fortrandata_clearpmemacrocell_")
      end subroutine FortranData_clearPmeMacrocell

      ! CG-09 optional staging seam.  The canonical topology arrays retain
      ! their Fortran ids; the GPU owner validates and copies them during the
      ! normal GpuSimulation allocation lifecycle.
      subroutine FortranData_setAdaptiveTopology(geometry_mode, atoms, blocks, basis, fft_channels, &
            fft_grid_channels, dynamic_channels, ensembles, selector_criteria, repetition_shape, &
            block_shape, block_grid, cell_vectors, block_vectors, atom_to_block, atom_to_basis, &
            atom_to_dynamic_channel, atom_to_fft_channel, atom_to_fft_grid_index, &
            basis_to_dynamic_channel, basis_to_fft_channel, block_atom_count, block_atom_offset, &
            block_atoms, block_grid_coordinate, block_basis_population, block_fft_population, &
            block_dynamic_population, block_center, block_volume, block_state, pending_state, &
            state_age, transition_epoch, selector_scores, coarse_moment, coarse_direction, &
            coarse_field, channel_moment_sum) bind(C, name="fortrandata_setadaptivetopology_")
         import :: c_int, c_double
         integer(c_int), intent(inout) :: geometry_mode, atoms, blocks, basis, fft_channels, fft_grid_channels
         integer(c_int), intent(inout) :: dynamic_channels, ensembles, selector_criteria
         integer(c_int), intent(inout) :: repetition_shape(*), block_shape(*), block_grid(*)
         real(c_double), intent(inout) :: cell_vectors(*), block_vectors(*)
         integer(c_int), intent(inout) :: atom_to_block(*), atom_to_basis(*)
         integer(c_int), intent(inout) :: atom_to_dynamic_channel(*), atom_to_fft_channel(*)
         integer(c_int), intent(inout) :: atom_to_fft_grid_index(*), basis_to_dynamic_channel(*)
         integer(c_int), intent(inout) :: basis_to_fft_channel(*), block_atom_count(*)
         integer(c_int), intent(inout) :: block_atom_offset(*), block_atoms(*)
         integer(c_int), intent(inout) :: block_grid_coordinate(*), block_basis_population(*)
         integer(c_int), intent(inout) :: block_fft_population(*), block_dynamic_population(*)
         real(c_double), intent(inout) :: block_center(*), block_volume(*)
         integer(c_int), intent(inout) :: block_state(*), pending_state(*)
         integer(c_int), intent(inout) :: state_age(*), transition_epoch(*)
         real(c_double), intent(inout) :: selector_scores(*), coarse_moment(*)
         real(c_double), intent(inout) :: coarse_direction(*), coarse_field(*)
         real(c_double), intent(inout) :: channel_moment_sum(*)
      end subroutine FortranData_setAdaptiveTopology

      subroutine FortranData_clearAdaptiveTopology() bind(C, name="fortrandata_clearadaptivetopology_")
      end subroutine FortranData_clearAdaptiveTopology

      subroutine FortranData_setAdaptiveKernels(atom_moment, projection_block, projection_weight, &
            bonds, bond_atom, bond_matrix, selector_edges, selector_edge, inverse_block_transpose, &
            exchange_stiffness, spiralization, anisotropy_axis_count, anisotropy_axis, &
            anisotropy_k1, anisotropy_k2, normalization_floor, magnetic_moment_si, gamma_per_ts, &
            damping, adaptive_mask, update_interval, refine_threshold, coarsen_threshold, &
            minimum_dwell, buffer_dilation, reconstruction_scheme, cone_angle_rad) &
            bind(C, name="fortrandata_setadaptivekernels_")
         import :: c_int, c_double
         real(c_double), intent(inout) :: atom_moment(*), projection_weight(*), bond_matrix(*)
         real(c_double), intent(inout) :: inverse_block_transpose(*), exchange_stiffness(*)
         real(c_double), intent(inout) :: spiralization(*), anisotropy_axis(*)
         real(c_double), intent(inout) :: anisotropy_k1(*), anisotropy_k2(*)
         real(c_double), intent(inout) :: normalization_floor, magnetic_moment_si, gamma_per_ts, damping
         real(c_double), intent(inout) :: refine_threshold, coarsen_threshold, cone_angle_rad
         integer(c_int), intent(inout) :: projection_block(*), bonds, bond_atom(*)
         integer(c_int), intent(inout) :: selector_edges, selector_edge(*), anisotropy_axis_count(*)
         integer(c_int), intent(inout) :: adaptive_mask, update_interval, minimum_dwell
         integer(c_int), intent(inout) :: buffer_dilation, reconstruction_scheme
      end subroutine FortranData_setAdaptiveKernels
   end interface


   private

   public :: fortran_do_measurements,fortran_measure,fortran_measure_moment,        &
      fortran_moment_update,fortran_flush_measurements,FortranData_Initiate,        &
      fortran_calc_simulation_status_variables, fortran_print_measurables,          &
      fortran_print_correlations, fortran_measure_correlations,                     &
      fortran_measure_rest, fortran_do_correlations
   public :: FortranData_setAdaptiveTopology, FortranData_clearAdaptiveTopology

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
   subroutine fortran_measure(cmstep) bind(C, name="fortran_measure")
      implicit none
      integer(c_int), intent(in) :: cmstep !< Current simulation step

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
   subroutine fortran_measure_correlations(ext_emomM, ext_emom, ext_mmom, ext_mstep) &
         bind(C, name="fortran_measure_correlations")
      implicit none
      real(c_double), intent(in) :: ext_emom(*)
      real(c_double), intent(in) :: ext_emomM(*)
      real(c_double), intent(in) :: ext_mmom(*)
      integer(c_int), intent(in) :: ext_mstep

      integer :: cgk_flag
      cgk_flag=0

      ! Spin correlation
      ! Sample magnetic moments for correlation functions
         call correlation_wrapper(Natom,Mensemble,coord,simid,emomM,ext_mstep,delta_t,  &
         NT_meta,atype_meta,Nchmax,achtype,sc,do_sc,do_sr,cgk_flag)

   end subroutine fortran_measure_correlations


   ! Measurements with pre-set parameters
   subroutine fortran_measure_moment(ext_emomM, ext_emom, ext_mmom, ext_mstep) &
         bind(C, name="fortran_measure_moment")
      implicit none
      type(c_ptr), value :: ext_emomM, ext_emom, ext_mmom
      integer(c_int), intent(in) :: ext_mstep
      real(c_double), pointer :: f_emomM(:,:,:), f_emom(:,:,:), f_mmom(:,:)
      real(dblprec) ::  totene, totenergy 

      integer :: cgk_flag
      cgk_flag=0
      call c_f_pointer(ext_emomM, f_emomM, [3, Natom, Mensemble])
      call c_f_pointer(ext_emom, f_emom, [3, Natom, Mensemble])
      call c_f_pointer(ext_mmom, f_mmom, [Natom, Mensemble])
      
      call measure(Natom,Mensemble,NT,NA,nHam,N1,N2,N3,simid,ext_mstep,f_emom,    &
         f_emomM,f_mmom,Nchmax,do_ralloy,Natom_full,asite_ch,achem_ch,atype,    &
         plotenergy,Temp,1.0_dblprec,0.0_dblprec,real_time_measure,delta_t,logsamp,             &
         ham%max_no_neigh,ham%nlist,ham%ncoup,ham%nlistsize,ham%aham,thermal_field, &
         beff,beff1,beff3,coord,ham%ind_list_full,ham%ind_nlistsize,ham%ind_nlist,  &
         ham%max_no_neigh_ind,ham%sus_ind,do_mom_legacy,mode)

      ! GPU measurements own the energy reduction and totenergy writer.  In that
      ! mode do_dip intentionally remains zero on the Fortran side, so invoking
      ! the legacy calculator here would append a competing all-zero Dip/total
      ! row after the GPU writer has produced the physical result.
         if(plotenergy>0 .and. do_gpu_measurements /= 'Y') then
            if (mod(ext_mstep-1,ene_step)==0) then
               call calc_energy(nHam,ext_mstep,Natom,Nchmax, &
                  conf_num,Mensemble,Natom,Num_macro,1,         &
                  plotenergy,Temp,delta_t,do_lsf,        &
                  lsf_field,lsf_interpolate,real_time_measure,simid,cell_index,            &
                  macro_nlistsize,f_mmom,f_emom,f_emomM,emomM_macro,external_field,              &
                  time_external_field,totene,NA,N1,N2,N3)
            end if
         endif


      ! Spin correlation
      ! Sample magnetic moments for correlation functions
      !   call correlation_wrapper(Natom,Mensemble,coord,simid,emomM,ext_mstep,delta_t,  &
      !   NT_meta,atype_meta,Nchmax,achtype,sc,do_sc,do_sr,cgk_flag)

   end subroutine fortran_measure_moment

      ! Measurements not implemeted or not planned to be implementedon GPU
   subroutine fortran_measure_rest(ext_emomM, ext_emom, ext_mmom, ext_beff, ext_mstep) &
         bind(C, name="fortran_measure_rest")
      use Math_functions, only : f_logstep
      use prn_trajectories, only : print_trajectories
      use Polarization,     only : print_pol
      use prn_fields,       only : print_fields
      use Temperature

      implicit none
      type(c_ptr), value :: ext_emomM, ext_emom, ext_mmom, ext_beff
      integer(c_int), intent(in) :: ext_mstep
      real(c_double), pointer :: f_emomM(:,:,:), f_emom(:,:,:), f_mmom(:,:), f_beff(:,:,:)

      integer :: sstep
      call c_f_pointer(ext_emomM, f_emomM, [3, Natom, Mensemble])
      call c_f_pointer(ext_emom, f_emom, [3, Natom, Mensemble])
      call c_f_pointer(ext_mmom, f_mmom, [Natom, Mensemble])
      call c_f_pointer(ext_beff, f_beff, [3, Natom, Mensemble])

      sstep = f_logstep(ext_mstep,logsamp)

      ! Sample spin temperature
      if (do_spintemp=='Y') then
         if(mod(ext_mstep,spintemp_step)==0) then
            call spintemperature(Natom,Mensemble,ext_mstep,1,simid,f_emomM,f_beff,1)
         end if
      endif

      call print_trajectories(Natom,sstep,ext_mstep,Mensemble,f_emom,f_mmom,delta_t,        &
         real_time_measure,simid,do_mom_legacy,mode)

      call print_pol(sstep,ext_mstep,Natom,Mensemble,ham%max_no_neigh,ham%nlist,ham%nlistsize,f_emomM,&
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
   subroutine fortran_do_measurements(cmstep, do_copy) bind(C, name="fortran_do_measurements")
      implicit none
      integer(c_int), intent(in) :: cmstep !< Current simulation step
      integer(c_int), intent(out) :: do_copy !< Flag if copy or not

      call do_measurements(cmstep,do_avrg,do_tottraj,avrg_step,ntraj,tottraj_step,  &
           traj_step,do_cumu,cumu_step,logsamp,do_copy,do_gpu_measurements)
   end subroutine fortran_do_measurements

   !---------------------------------------------------------------------
   !> @brief
   !> Report whether the current step is a spin-correlation sampling step.
   !> Mirrors the internal gating of correlation_wrapper (calc_gk2 for the
   !> static S(q), calc_gkt for the dynamic S(q,t)) so that FortranCorrelation
   !> hands emomM to the CPU sampler on exactly the steps it consumes it. The
   !> measurement cadence (fortran_do_measurements: avrg_step/cumu_step/...) does
   !> not in general coincide with the correlation cadence (sc_step/sc_sep), so
   !> the two paths must use separate gates.
   !---------------------------------------------------------------------
   subroutine fortran_do_correlations(cmstep, do_copy) bind(C, name="fortran_do_correlations")
      implicit none
      integer(c_int), intent(in)  :: cmstep  !< Current simulation step
      integer(c_int), intent(out) :: do_copy !< Flag if moment must be copied

      do_copy = 0
      ! Static S(q): sampled every sc_sep steps
      if (do_sc=='C'.or.do_sc=='Y') then
         if (sc%sc_sep>0) then
            if (mod(cmstep-1,sc%sc_sep)==0) do_copy = 1
         end if
      end if
      ! Dynamic S(q,t): sampled every sc_step steps (offset by 1, as in
      ! correlation_wrapper). correlation_wrapper self-limits on sc_max_nstep,
      ! so no sc_tidx guard is needed here (and it would race the async queue).
      if (do_sc=='Q'.or.do_sc=='T'.or.do_sc=='Y') then
         if (sc%sc_step>0) then
            if (mod(cmstep-1,sc%sc_step)==1) do_copy = 1
         end if
      end if
   end subroutine fortran_do_correlations



   ! Moment update with pre-set parameters
   subroutine fortran_moment_update() bind(C, name="fortran_moment_update")
      implicit none
      call moment_update(Natom,Mensemble,mmom,mmom0,mmom2,emom,emom2,emomM,mmomi,   &
         mompar,initexc)
   end subroutine fortran_moment_update



   ! Flush measurements with pre-set parameters
   subroutine fortran_flush_measurements(cmstep) bind(C, name="fortran_flush_measurements")
      implicit none
      integer(c_int), intent(in) :: cmstep !< Current simulation stepfind_rmid(rmid,coord,Natom)
      call flush_measurements(Natom,Mensemble,NT,NA,N1,N2,N3,simid,cmstep,emom,mmom,&
         Nchmax,atype,real_time_measure,mcnstep,ham%ind_list_full,do_mom_legacy,mode)
   end subroutine fortran_flush_measurements

      ! print GPU calculated correlations
   subroutine fortran_print_correlations() bind(C, name="fortran_print_correlations")
      implicit none
      integer ::cgk_flag
      cgk_flag = 2
      !type(corr_t), intent(inout) :: cc !< Derived type for correlation data
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

   subroutine fortran_set_mstep(value) bind(C, name="fortran_set_mstep")
      integer(c_int), intent(in) :: value
      mstep = value
   end subroutine fortran_set_mstep

   subroutine fortran_get_mstep(value) bind(C, name="fortran_get_mstep")
      integer(c_int), intent(out) :: value
      value = mstep
   end subroutine fortran_get_mstep



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
      integer :: ene_step = 100
      integer :: ene_buff = 100
      !character(len=1)::do_projch_avrg = 'Y'
      !character(len=1)::do_cumu_proj = 'Y'


      !if(cc%do_proj=='C'.or.cc%do_proj=='Y'.or.cc%do_proj=='T'.or.cc%do_proj=='Q'.or.cc%do_projch=='C'.or.cc%do_projch=='Y'.or.cc%do_projch=='Q'.or.cc%do_projch=='T') then
      !   print *, "Projections are not available in GPU correlations yet, please use do_gpu_correlations 0"
      !   return  
      !end if

      if(do_gpu_correlations=='Y'.and.phase=='M') then
         ! Calculate r_mid for GPU correlations, as it is needed for the correlation calculations and not calculated on the Fortran side otherwise

         if (do_sc=='Y' .or. do_sc=='C') then
            ! Initializing Fortran-side correlation arrays. Here for static correlation function
            cc%gk_flag=2
            allocate(cc%m_k(3,nq))
            cc%m_k=0.0_dblprec
            if(cc%do_proj=='C'.or.cc%do_proj=='Y') then
               allocate(cc%m_k_proj(3,NT_meta,nq))
               cc%m_k_proj=0.0_dblprec
            end if

            if(cc%do_projch=='C'.or.cc%do_projch=='Y') then
               allocate(cc%m_k_projch(3,Nchmax,nq))
               cc%m_k_projch=0.0_dblprec
            end if

            ! call find_rmid(r_mid,coord,Natom)
            ! Possible alternative: Call the Fortran routine with init flag = 0
            ! zeroflag = 0
            ! call calc_gk2(Natom, Mensemble, NT,atype,Nchmax,achtype, cc, coord, simid, emomM, zeroflag)
         end if

         ! Projected static correlations are needed whenever projected mode includes static terms,
         ! even if base do_sc='Q' and m_k itself is not allocated.
         if ((cc%do_proj=='C'.or.cc%do_proj=='Y') .and. (.not.allocated(cc%m_k_proj))) then
            allocate(cc%m_k_proj(3,NT_meta,nq))
            cc%m_k_proj=0.0_dblprec
         end if

         if ((cc%do_projch=='C'.or.cc%do_projch=='Y') .and. (.not.allocated(cc%m_k_projch))) then
            allocate(cc%m_k_projch(3,Nchmax,nq))
            cc%m_k_projch=0.0_dblprec
         end if
         if (do_sc=='Y' .or. do_sc=='Q') then
            ! Initializing Fortran-side correlation arrays. Here for dynamic correlation function
            cc%gkt_flag=2
            allocate(cc%m_kt(3,nq,cc%sc_max_nstep))
            cc%m_kt=0.0_dblprec
            call allocate_deltatcorr(.true.,cc)
            allocate(cc%m_kw(3,nq,cc%nw))
            cc%m_kw=0.0_dblprec
            if(cc%do_proj=='Q'.or.cc%do_proj=='T'.or.cc%do_proj=='Y') then
               allocate(cc%m_kt_proj(3,NT_meta,nq,cc%sc_max_nstep))
               cc%m_kt_proj=0.0_dblprec
               allocate(cc%m_kw_proj(3,NT_meta,nq,cc%nw))
               cc%m_kw_proj=0.0_dblprec
            end if

            if(cc%do_projch=='Q'.or.cc%do_projch=='T'.or.cc%do_projch=='Y') then
               allocate(cc%m_kt_projch(3,Nchmax,nq,cc%sc_max_nstep))
               cc%m_kt_projch=0.0_dblprec
               allocate(cc%m_kw_projch(3,Nchmax,nq,cc%nw))
               cc%m_kw_projch=0.0_dblprec
            end if

            ! Possible alternative: Call the Fortran routine with init flag = 0
            ! zeroflag = 0
            ! cc%gkt_flag=0
            ! call calc_gkt(Natom, Mensemble, NT,atype,Nchmax,achtype, cc, coord, emomM, zeroflag)
         end if
         call find_rmid(r_mid,coord,Natom)

      end if
     

      call FortranData_setFlags(ham_inp%do_dm, ham_inp%do_jtensor, ham_inp%do_anisotropy, &
           do_avrg, do_proj_avrg, do_projch_avrg, do_cumu, do_cumu_proj, plotenergy, do_autocorr, do_tottraj, ntraj, &
           do_gpu_measurements, skyno, do_sc, do_gpu_correlations, do_gpu_timings, real_time_measure, &
           cc%do_proj, cc%do_projch, do_ralloy)

      call FortranData_setConstants(stt,SDEalgh,rstep,nstep,Natom,Mensemble, &
         ham%max_no_neigh,delta_t,gama,k_bolt,mub,mplambda1,binderc,mavg,mompar, &
         initexc,ham%max_no_dmneigh,nHam, Temp, ipmcnphase, mcnstep, ipnphase, &
         avrg_step, avrg_buff, cumu_step, cumu_buff, ene_step, ene_buff, &
         tottraj_step, tottraj_buff, skyno_step, skyno_buff, nq, sc_window_fun, &
         cc%nw, cc%sc_sep, cc%sc_step, cc%sc_max_nstep, nspinwait, ac_step, ac_buff, &
         NT_meta, Nchmax, mry, NA, Natom_full)

      call FortranData_setGpuGeometry(N1, N2, N3, NA, BC1, BC2, BC3, C1, C2, C3, Bas, alat, do_gpu_convolution)

      select case(trim(gpu_dipole_mode))
      case('OFF')
         gpu_dipole_mode_id=0
      case('EWALD3D_FFT')
         gpu_dipole_mode_id=1
      case('OPEN_FFT')
         gpu_dipole_mode_id=2
      case('PME3D')
         error stop 'gpu_dipole_mode PME3D was renamed; use EWALD3D_FFT'
      case default
         error stop 'Invalid gpu_dipole_mode (OFF, EWALD3D_FFT, OPEN_FFT)'
      end select
      select case(trim(gpu_dipole_surface))
      case('TINFOIL')
         gpu_dipole_surface_id=0
      case('VACUUM_SPHERE')
         gpu_dipole_surface_id=1
      case default
         error stop 'Invalid gpu_dipole_surface (TINFOIL, VACUUM_SPHERE)'
      end select
      if(gpu_dipole_alpha < 0.0_dblprec .or. gpu_dipole_rcut < 0.0_dblprec .or. any(gpu_dipole_mesh < 0)) then
         error stop 'GPU dipole alpha, cutoff, and mesh must be non-negative (zero selects auto)'
      endif
      call FortranData_setGpuDipole(gpu_dipole_mode_id, gpu_dipole_surface_id, gpu_dipole_alpha, &
         gpu_dipole_rcut, gpu_dipole_tol, gpu_dipole_mesh)

      ! Macrocell arrays are allocated only when macro cells are requested.
      ! Do not pass an unallocated Fortran allocatable through the C ABI.
      if (Num_macro > 0) then
         call FortranData_setMacrocell(ham_inp%do_dip, Num_macro, block_size_x, block_size_y, block_size_z, &
            cell_index, macro_nlistsize, gpu_macro_center, gpu_macro_min_coord, gpu_macro_max_coord)
      else
         call FortranData_clearMacrocell()
      endif
      if (pme_Num_macro > 0) then
         call FortranData_setPmeMacrocell(pme_Num_macro, pme_macro_grid, pme_cell_index, pme_macro_nlistsize, &
            pme_macro_center, pme_macro_min_coord, pme_macro_max_coord)
      else
         call FortranData_clearPmeMacrocell()
      endif

      call FortranData_setHamiltonian(ham%ncoup,ham%nlist,ham%nlistsize, &
         ham%dm_vect,ham%dmlist,ham%dmlistsize, &
         ham%kaniso, ham%eaniso, ham%taniso, ham%sb, &
         ham%j_tens, ham%aHam, &
         external_field, btorque,Temp_array, &
          ipTemp, ipmcnstep, ipTemp_array, ipnstep, ipdelta_t, iplambda1)

      call FortranData_setLattice(beff, b2eff, emomM, emom, emom2, mmom, mmom0, mmom2, mmomi, &
         dxyz_vec, dxyz_atom, dxyz_list, atype)


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
          cc%scstep_arr, cc%sc_nsamp, cc%sc_tidx, atype_meta, achtype, cc%m_k_proj, cc%m_k_projch, &
          cc%m_kt_proj, cc%m_kt_projch, cc%m_kw_proj, cc%m_kw_projch)

      ! Feature-off is an explicit null sentinel.  Enabled GPU runs stage the
      ! complete production topology/runtime/kernel inventory here, before
      ! GpuSimulation performs memory preflight and copies it to owned storage.
      call FortranData_clearAdaptiveTopology()
      if (adaptive_cg_is_enabled() .and. adaptive_cg_state%gpu_requested) then
         call FortranData_setAdaptiveTopology( &
            adaptive_cg_state%topology%geometry_mode,adaptive_cg_state%topology%n_atoms, &
            adaptive_cg_state%topology%n_spatial_blocks,adaptive_cg_state%topology%n_basis, &
            adaptive_cg_state%topology%n_fft_channels_per_block, &
            adaptive_cg_state%topology%n_fft_grid_channels, &
            adaptive_cg_state%topology%n_dynamic_channels,Mensemble, &
            adaptive_cg_state%gpu_selector_criteria,adaptive_cg_state%topology%repetition_shape, &
            adaptive_cg_state%topology%block_shape,adaptive_cg_state%topology%block_grid, &
            adaptive_cg_state%topology%cell_vectors,adaptive_cg_state%topology%block_vectors, &
            adaptive_cg_state%topology%atom_to_block,adaptive_cg_state%topology%atom_to_basis, &
            adaptive_cg_state%topology%atom_to_dynamic_channel, &
            adaptive_cg_state%topology%atom_to_fft_channel, &
            adaptive_cg_state%topology%atom_to_fft_grid_index, &
            adaptive_cg_state%topology%basis_to_dynamic_channel, &
            adaptive_cg_state%topology%basis_to_fft_channel, &
            adaptive_cg_state%topology%block_atom_count, &
            adaptive_cg_state%topology%block_atom_offset,adaptive_cg_state%topology%block_atoms, &
            adaptive_cg_state%topology%block_grid_coordinate, &
            adaptive_cg_state%topology%block_basis_population, &
            adaptive_cg_state%topology%block_fft_channel_population, &
            adaptive_cg_state%topology%block_dynamic_channel_population, &
            adaptive_cg_state%topology%block_center,adaptive_cg_state%topology%block_volume, &
            adaptive_cg_state%runtime%hybrid%block_state, &
            adaptive_cg_state%gpu_pending_state,adaptive_cg_state%gpu_state_age, &
            adaptive_cg_state%gpu_transition_epoch,adaptive_cg_state%gpu_selector_scores, &
            adaptive_cg_state%runtime%coarse_resultant_mub,adaptive_cg_state%coarse_direction, &
            adaptive_cg_state%gpu_coarse_field,adaptive_cg_state%runtime%channel_moment_sum_mub)
         call FortranData_setAdaptiveKernels( &
            adaptive_cg_state%atom_moment_mub,adaptive_cg_state%projection%stencil_block, &
            adaptive_cg_state%projection%shape_weight,adaptive_cg_state%gpu_bonds, &
            adaptive_cg_state%bond_atom,adaptive_cg_state%bond_matrix_j, &
            adaptive_cg_state%gpu_selector_edges, &
            adaptive_cg_state%bond_atom,adaptive_cg_state%tensor%inverse_block_transpose_m1, &
            adaptive_cg_state%tensor%exchange_stiffness_j_per_m, &
            adaptive_cg_state%tensor%spiralization_j_per_m2, &
            adaptive_cg_state%gpu_anisotropy_axis_count,adaptive_cg_state%gpu_anisotropy_axis, &
            adaptive_cg_state%gpu_anisotropy_k1,adaptive_cg_state%gpu_anisotropy_k2, &
            adaptive_cg_state%projection%normalization_floor, &
            adaptive_cg_state%gpu_magnetic_moment_si, &
            adaptive_cg_state%tensor%channel_gamma_per_t_s, &
            adaptive_cg_state%tensor%channel_damping,adaptive_cg_state%gpu_adaptive_mask, &
            adaptive_cg%update_interval,adaptive_cg%refine_threshold, &
            adaptive_cg%coarsen_threshold,adaptive_cg%minimum_dwell_updates, &
            adaptive_cg%buffer_blocks,adaptive_cg_state%gpu_reconstruction_scheme, &
            adaptive_cg_state%reconstruction%cone_angle_rad)
      endif

      call FortranData_setInputData(gpu_mode, gpu_rng, gpu_rng_seed)

   end subroutine FortranData_Initiate

       ! Print measurables calculated in CUDA
   subroutine fortran_print_measurables(obs_step, obs_buff, indxb_obs, obs_name, &
        obs_label, obs_dim, obs_buffer, mstep) bind(C, name="fortran_print_measurables")
      implicit none
      integer(c_int), intent(in) :: obs_step, obs_buff, obs_dim, mstep
      type(c_ptr), value :: indxb_obs
      character(c_char), intent(in) :: obs_name(*), obs_label(*)
      type(c_ptr), value :: obs_buffer
      real(c_double), pointer :: f_indxb_obs(:)
      real(c_double), pointer :: f_obs_buffer(:,:,:)
      character(len=16) :: f_obs_name
      character(len=16), allocatable :: f_obs_label(:)
      integer :: i, j

      f_obs_name = ''
      do i = 1, len(f_obs_name)
         f_obs_name(i:i) = obs_name(i)
      end do

      allocate(f_obs_label(obs_dim))
      do i = 1, obs_dim
         f_obs_label(i) = ''
         do j = 1, len(f_obs_label(i))
            f_obs_label(i)(j:j) = obs_label((i - 1) * len(f_obs_label(i)) + j)
         end do
      end do

      call c_f_pointer(indxb_obs, f_indxb_obs, [obs_buff])
      call c_f_pointer(obs_buffer, f_obs_buffer, [obs_dim, obs_buff, Mensemble])
      call print_observable(simid, Mensemble, f_obs_name, obs_step, obs_buff, &
         obs_dim, f_indxb_obs, f_obs_buffer, f_obs_label, real_time_measure, delta_t, mstep)
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
