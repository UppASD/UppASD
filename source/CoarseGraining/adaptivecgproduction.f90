!-------------------------------------------------------------------------------
! MODULE: AdaptiveCGProduction
!> Production ownership, validation, setup, CPU dispatch, diagnostics, and
!> cleanup for the first supported adaptive coarse-graining capability.
!-------------------------------------------------------------------------------
module AdaptiveCGProduction

   use, intrinsic :: iso_c_binding, only : c_int, c_double
   use, intrinsic :: iso_fortran_env, only : int64
   use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
   use Parameters, only : dblprec
   use Constants, only : gama, mub, pi
   use InputData
   use SystemData, only : coord, anumb, Landeg
   use MomentData, only : emom, mmom, emomM
   use HamiltonianData, only : ham
   use Damping, only : lambda1_array
   use macrocells, only : block_size_x, block_size_y, block_size_z
   use stiffness, only : coarse_material_type, coarse_lattice_input_type, eta_min, eta_max, &
      extract_coarse_material_from_uppasd, validate_coarse_material_small_q, &
      COARSE_MATERIAL_OK
   use QHB, only : do_qhb
   use Temperature, only : do_3tm
   use AutoCorrelation, only : do_autocorr
   use KMCData, only : do_kmc
   use FixedMom, only : do_fixed_mom
   use BlockTopology, only : block_topology_type, build_block_topology, &
      REGULAR_REPLICATED_CELL, BLOCK_TOPOLOGY_OK
   use CoarseTensorOperator, only : coarse_tensor_operator_type, &
      coarse_operator_options_type, setup_coarse_tensor_operator, coarse_llg_rhs, &
      COARSE_TENSOR_OK
   use SmoothProjectedOperator, only : smooth_projected_operator_type, &
      setup_smooth_projected_operator, SMOOTH_PROJECTED_OK
   use StaticHybridOperator, only : static_hybrid_operator_type, &
      static_hybrid_energy_type, setup_static_hybrid_operator, &
      evaluate_static_hybrid_operator, STATIC_HYBRID_OK
   use BlockSelector, only : selector_configuration_type, selector_registry_type, &
      selector_evaluation_type, selector_requests_type, register_max_neighbour_misalignment, &
      evaluate_selector_registry, combine_selector_requests, SELECTOR_OK
   use AdaptiveHybridSolver, only : adaptive_hybrid_runtime_type, &
      adaptive_reconstruction_configuration_type, setup_adaptive_hybrid_runtime, &
      apply_adaptive_transitions, reconstruct_block_aligned, &
      reconstruct_block_constrained_cone, ADAPTIVE_HYBRID_OK, &
      ADAPTIVE_STAGE_COMPLETE_STEP, RECONSTRUCTION_ALIGNED, &
      RECONSTRUCTION_CONSTRAINED_CONE

   implicit none
   private

   integer, parameter, public :: ADAPTIVE_CG_PRODUCTION_OK = 0
   integer, parameter, public :: ADAPTIVE_CG_PRODUCTION_REJECTED = 1
   integer, parameter, public :: ADAPTIVE_CG_PRODUCTION_SETUP_FAILED = 2

   type, public :: adaptive_cg_production_state_type
      logical :: ready = .false.
      logical :: gpu_requested = .false.
      logical :: adaptive_mask = .false.
      integer :: completed_steps = 0
      integer(int64) :: field_evaluations = 0_int64
      integer(int64) :: active_atom_updates = 0_int64
      integer(int64) :: active_block_updates = 0_int64
      integer(int64) :: accepted_transitions = 0_int64
      integer(int64) :: rejected_transitions = 0_int64
      real(dblprec) :: field_seconds = 0.0_dblprec
      real(dblprec) :: integration_seconds = 0.0_dblprec
      real(dblprec) :: reconstruction_seconds = 0.0_dblprec
      real(dblprec) :: selector_seconds = 0.0_dblprec
      real(dblprec) :: last_atom_field_sum_t = 0.0_dblprec
      real(dblprec) :: last_atom_field_norm2_t2 = 0.0_dblprec
      real(dblprec) :: last_coarse_field_sum_t = 0.0_dblprec
      real(dblprec) :: last_coarse_field_norm2_t2 = 0.0_dblprec
      type(static_hybrid_energy_type) :: last_energy
      type(block_topology_type) :: topology
      type(coarse_material_type) :: material
      type(coarse_tensor_operator_type) :: tensor
      type(smooth_projected_operator_type) :: projection
      type(static_hybrid_operator_type) :: hybrid
      type(adaptive_hybrid_runtime_type) :: runtime
      type(selector_configuration_type) :: selector_configuration
      type(selector_registry_type) :: selector_registry
      type(adaptive_reconstruction_configuration_type) :: reconstruction
      logical, allocatable :: initial_fine_mask(:)
      logical, allocatable :: hard_fine_mask(:)
      real(dblprec), allocatable :: atom_moment_mub(:)
      real(dblprec), allocatable :: coarse_direction(:,:,:,:)
      integer, allocatable :: bond_atom(:,:)
      real(dblprec), allocatable :: bond_matrix_j(:,:,:)
      real(dblprec), allocatable :: bond_displacement_m(:,:)
      ! Persistent host staging for the normal Fortran/C GPU setup seam.
      integer(c_int) :: gpu_selector_criteria = 1_c_int
      integer(c_int) :: gpu_bonds = 0_c_int
      integer(c_int) :: gpu_selector_edges = 0_c_int
      integer(c_int) :: gpu_adaptive_mask = 0_c_int
      integer(c_int) :: gpu_reconstruction_scheme = 1_c_int
      integer(c_int) :: gpu_diagnostics = 0_c_int
      integer(c_int) :: gpu_buffer_dilation = 0_c_int
      real(c_double) :: gpu_magnetic_moment_si = 9.274009994d-24
      integer(c_int), allocatable :: gpu_pending_state(:)
      integer(c_int), allocatable :: gpu_state_age(:)
      integer(c_int), allocatable :: gpu_transition_epoch(:)
      integer(c_int), allocatable :: gpu_anisotropy_axis_count(:)
      real(c_double), allocatable :: gpu_selector_scores(:,:)
      real(c_double), allocatable :: gpu_coarse_field(:,:,:,:)
      real(c_double), allocatable :: gpu_anisotropy_axis(:,:,:)
      real(c_double), allocatable :: gpu_anisotropy_k1(:,:)
      real(c_double), allocatable :: gpu_anisotropy_k2(:,:)
   end type adaptive_cg_production_state_type

   type(adaptive_cg_production_state_type), save, public :: adaptive_cg_state

   public :: setup_adaptive_cg_production
   public :: cleanup_adaptive_cg_production
   public :: adaptive_cg_is_enabled
   public :: adaptive_cg_cpu_step
   public :: print_adaptive_cg_summary

contains

   logical function adaptive_cg_is_enabled()
      adaptive_cg_is_enabled = adaptive_cg_state%ready
   end function adaptive_cg_is_enabled

   subroutine setup_adaptive_cg_production(status, diagnostic)
      integer, intent(out) :: status
      character(len=*), intent(out) :: diagnostic

      type(coarse_lattice_input_type) :: validation_input
      type(coarse_operator_options_type) :: options
      integer(c_int) :: topology_status, selector_status
      integer :: local_status, basis, atom, ensemble, block
      integer :: basis_channel(NA)
      real(dblprec) :: cell_vectors(3,3), inverse_block(3,3)
      real(dblprec) :: fractional_block_coordinate(3,Natom)
      real(dblprec) :: channel_gamma(1), channel_damping(1)
      real(dblprec) :: q(3,3), normal(3,3), phase(1,3)
      real(dblprec), allocatable :: production_ncoup(:,:,:,:)
      real(dblprec), allocatable :: production_dm_vect(:,:,:)
      integer, allocatable :: production_dmlist(:,:), production_dmlistsize(:)
      character(len=512) :: message

      status = ADAPTIVE_CG_PRODUCTION_OK
      diagnostic = ''
      if (same_word(adaptive_cg%enabled,'N')) return

      call validate_configuration(status,diagnostic)
      if (status /= ADAPTIVE_CG_PRODUCTION_OK) return
      adaptive_cg_state%gpu_requested = affirmative(do_gpu) .and. affirmative(do_gpu_llg)
      adaptive_cg_state%adaptive_mask = same_word(adaptive_cg%mask_mode,'ADAPTIVE')

      basis_channel = -1
      do basis = 1, NA
         if (abs(ammom_inp(basis,1,1)) > tiny(1.0_dblprec)) basis_channel(basis) = 1
      end do
      if (same_word(adaptive_cg%channel_mode,'FILE')) then
         call read_channel_file(trim(adaptive_cg%channel_file),basis_channel,status,diagnostic)
         if (status /= ADAPTIVE_CG_PRODUCTION_OK) return
      end if
      if (maxval(basis_channel) /= 1 .or. any(basis_channel > 1)) then
         call reject('cg_channel_mode supports exactly one dynamical FM channel',status,diagnostic)
         return
      end if

      cell_vectors(:,1) = alat*C1
      cell_vectors(:,2) = alat*C2
      cell_vectors(:,3) = alat*C3
      call build_block_topology(adaptive_cg_state%topology,REGULAR_REPLICATED_CELL,NA, &
         (/N1,N2,N3/),Natom,(/block_size_x,block_size_y,block_size_z/),cell_vectors, &
         basis_channel,topology_status,message)
      if (topology_status /= BLOCK_TOPOLOGY_OK) then
         call setup_failed('BlockTopology: '//trim(message),status,diagnostic)
         return
      end if

      allocate(adaptive_cg_state%initial_fine_mask(adaptive_cg_state%topology%n_spatial_blocks), &
         adaptive_cg_state%hard_fine_mask(adaptive_cg_state%topology%n_spatial_blocks))
      adaptive_cg_state%initial_fine_mask = .true.
      adaptive_cg_state%hard_fine_mask = .false.
      if (len_trim(adaptive_cg%static_mask_file) > 0) then
         adaptive_cg_state%initial_fine_mask = .false.
         call read_mask_file(trim(adaptive_cg%static_mask_file), &
            adaptive_cg_state%initial_fine_mask,status,diagnostic)
         if (status /= ADAPTIVE_CG_PRODUCTION_OK) then
            call cleanup_adaptive_cg_production()
            return
         end if
         if (adaptive_cg_state%adaptive_mask) &
            adaptive_cg_state%hard_fine_mask = adaptive_cg_state%initial_fine_mask
      end if

      channel_gamma(1) = Landeg(1)
      channel_damping(1) = lambda1_array(1)
      allocate(production_ncoup(1,size(ham%ncoup,1),size(ham%ncoup,2),1))
      production_ncoup(1,:,:,1)=ham%ncoup(:,:,1)
      allocate(production_dmlist(0,Natom),production_dmlistsize(nHam), &
         production_dm_vect(3,0,nHam))
      production_dmlistsize = 0
      call extract_coarse_material_from_uppasd(NA,N1,N2,N3,Natom,nHam,1,Nchmax, &
         ham%max_no_neigh,0,eta_min,eta_max,anumb,ham%aham,ham%nlistsize, &
         ham%nlist,production_dmlistsize,production_dmlist,alat,C1,C2,C3,BC1,BC2,BC3, &
         coord,ammom_inp,production_ncoup,production_dm_vect,basis_channel, &
         adaptive_cg_state%material, &
         local_status,message,channel_gamma,channel_damping,validation_input)
      if (local_status /= COARSE_MATERIAL_OK) then
         call setup_failed('coarse material extraction: '//trim(message),status,diagnostic)
         call cleanup_adaptive_cg_production()
         return
      end if
      q = 0.0_dblprec
      q(1,1) = 1.0d-6/alat
      q(2,2) = 1.0d-6/alat
      q(3,3) = 1.0d-6/alat
      normal = 0.0_dblprec
      normal(3,:) = 1.0_dblprec
      phase = 0.0_dblprec
      call validate_coarse_material_small_q(validation_input,adaptive_cg_state%material, &
         q,normal,phase,1.0d-6,1.0d-3,local_status,message)
      if (local_status /= COARSE_MATERIAL_OK) then
         call setup_failed('material convention/small-q validation: '//trim(message),status,diagnostic)
         call cleanup_adaptive_cg_production()
         return
      end if

      options%temperature_k = 0.0_dblprec
      call setup_coarse_tensor_operator(adaptive_cg_state%tensor,adaptive_cg_state%topology, &
         adaptive_cg_state%material,options,local_status,message)
      if (local_status /= COARSE_TENSOR_OK) then
         call setup_failed('coarse tensor setup: '//trim(message),status,diagnostic)
         call cleanup_adaptive_cg_production()
         return
      end if

      allocate(adaptive_cg_state%atom_moment_mub(Natom))
      adaptive_cg_state%atom_moment_mub = mmom(:,1)
      inverse_block = inverse3(adaptive_cg_state%topology%block_vectors)
      do atom = 1, Natom
         fractional_block_coordinate(:,atom) = matmul(inverse_block,alat*coord(:,atom))
      end do
      call setup_smooth_projected_operator(adaptive_cg_state%projection, &
         adaptive_cg_state%topology,fractional_block_coordinate, &
         adaptive_cg_state%atom_moment_mub,local_status,message)
      if (local_status /= SMOOTH_PROJECTED_OK) then
         call setup_failed('projection setup: '//trim(message),status,diagnostic)
         call cleanup_adaptive_cg_production()
         return
      end if

      call build_unique_bonds(status,diagnostic)
      if (status /= ADAPTIVE_CG_PRODUCTION_OK) then
         call cleanup_adaptive_cg_production()
         return
      end if
      call setup_static_hybrid_operator(adaptive_cg_state%hybrid, &
         adaptive_cg_state%topology,adaptive_cg_state%tensor,adaptive_cg_state%projection, &
         adaptive_cg_state%initial_fine_mask,adaptive_cg_state%bond_atom, &
         adaptive_cg_state%bond_displacement_m,adaptive_cg%buffer_blocks,local_status,message)
      if (local_status /= STATIC_HYBRID_OK) then
         call setup_failed('hybrid ownership setup: '//trim(message),status,diagnostic)
         call cleanup_adaptive_cg_production()
         return
      end if

      allocate(adaptive_cg_state%coarse_direction(3,1, &
         adaptive_cg_state%topology%n_spatial_blocks,Mensemble))
      adaptive_cg_state%coarse_direction = 0.0_dblprec
      do ensemble = 1, Mensemble
         do block = 1, adaptive_cg_state%topology%n_spatial_blocks
            adaptive_cg_state%coarse_direction(:,1,block,ensemble) = &
               block_mean_direction(block,ensemble)
         end do
      end do
      call setup_adaptive_hybrid_runtime(adaptive_cg_state%runtime, &
         adaptive_cg_state%topology,adaptive_cg_state%hybrid, &
         adaptive_cg_state%atom_moment_mub,emom,adaptive_cg_state%coarse_direction, &
         local_status,message)
      if (local_status /= ADAPTIVE_HYBRID_OK) then
         call setup_failed('adaptive runtime setup: '//trim(message),status,diagnostic)
         call cleanup_adaptive_cg_production()
         return
      end if

      adaptive_cg_state%selector_configuration%refine_threshold = adaptive_cg%refine_threshold
      adaptive_cg_state%selector_configuration%coarsen_threshold = adaptive_cg%coarsen_threshold
      adaptive_cg_state%selector_configuration%update_interval = adaptive_cg%update_interval
      adaptive_cg_state%selector_configuration%minimum_dwell_updates = &
         adaptive_cg%minimum_dwell_updates
      adaptive_cg_state%selector_configuration%buffer_dilation_blocks = adaptive_cg%buffer_blocks
      adaptive_cg_state%reconstruction%scheme = merge(RECONSTRUCTION_CONSTRAINED_CONE, &
         RECONSTRUCTION_ALIGNED,same_word(adaptive_cg%reconstruction,'CONE'))
      adaptive_cg_state%reconstruction%cone_angle_rad = adaptive_cg%cone_angle_deg*pi/180.0_dblprec
      adaptive_cg_state%reconstruction%energy_jump_limit_j = adaptive_cg%energy_jump_limit_j
      if (adaptive_cg_state%adaptive_mask) then
         call register_max_neighbour_misalignment(adaptive_cg_state%selector_registry, &
            selector_status,message)
         if (selector_status /= SELECTOR_OK) then
            call setup_failed('cg_selector: '//trim(message),status,diagnostic)
            call cleanup_adaptive_cg_production()
            return
         end if
      end if

      if (adaptive_cg_state%gpu_requested) call allocate_gpu_staging()

      adaptive_cg_state%ready = .true.
      call print_resolved_configuration()
   end subroutine setup_adaptive_cg_production

   subroutine allocate_gpu_staging()
      integer :: blocks
      blocks = adaptive_cg_state%topology%n_spatial_blocks
      allocate(adaptive_cg_state%gpu_pending_state(blocks), &
         adaptive_cg_state%gpu_state_age(blocks), &
         adaptive_cg_state%gpu_transition_epoch(blocks), &
         adaptive_cg_state%gpu_anisotropy_axis_count(blocks), &
         adaptive_cg_state%gpu_selector_scores(1,blocks), &
         adaptive_cg_state%gpu_coarse_field(3,1,blocks,Mensemble), &
         adaptive_cg_state%gpu_anisotropy_axis(3,2,blocks), &
         adaptive_cg_state%gpu_anisotropy_k1(2,blocks), &
         adaptive_cg_state%gpu_anisotropy_k2(2,blocks))
      adaptive_cg_state%gpu_pending_state = int( &
         adaptive_cg_state%runtime%hybrid%block_state,c_int)
      adaptive_cg_state%gpu_bonds = size(adaptive_cg_state%bond_atom,2)
      adaptive_cg_state%gpu_selector_edges = adaptive_cg_state%gpu_bonds
      adaptive_cg_state%gpu_adaptive_mask = &
         merge(1_c_int,0_c_int,adaptive_cg_state%adaptive_mask)
      adaptive_cg_state%gpu_reconstruction_scheme = &
         adaptive_cg_state%reconstruction%scheme
      adaptive_cg_state%gpu_diagnostics = int(adaptive_cg%diagnostics,c_int)
      adaptive_cg_state%gpu_buffer_dilation = int(maxval( &
         adaptive_cg_state%runtime%hybrid%buffer_width_blocks),c_int)
      adaptive_cg_state%gpu_state_age = adaptive_cg_state%runtime%selector%state_age
      adaptive_cg_state%gpu_transition_epoch = &
         adaptive_cg_state%runtime%selector%transition_epoch
      adaptive_cg_state%gpu_anisotropy_axis_count = 0_c_int
      adaptive_cg_state%gpu_selector_scores = 0.0_c_double
      adaptive_cg_state%gpu_coarse_field = 0.0_c_double
      adaptive_cg_state%gpu_anisotropy_axis = 0.0_c_double
      adaptive_cg_state%gpu_anisotropy_k1 = 0.0_c_double
      adaptive_cg_state%gpu_anisotropy_k2 = 0.0_c_double
   end subroutine allocate_gpu_staging

   subroutine validate_configuration(status,diagnostic)
      integer, intent(out) :: status
      character(len=*), intent(out) :: diagnostic

      status = ADAPTIVE_CG_PRODUCTION_OK
      diagnostic = ''
      if (.not. same_word(adaptive_cg%enabled,'Y')) then
         call reject('do_adaptive_cg must be Y or N',status,diagnostic); return
      end if
      if (.not. same_word(adaptive_cg%operator,'TENSOR') .and. &
          .not. same_word(adaptive_cg%operator,'PROJECTED')) then
         call reject('cg_operator must be TENSOR or PROJECTED',status,diagnostic); return
      end if
      if (.not. same_word(adaptive_cg%mask_mode,'STATIC') .and. &
          .not. same_word(adaptive_cg%mask_mode,'ADAPTIVE')) then
         call reject('cg_mask_mode must be STATIC or ADAPTIVE',status,diagnostic); return
      end if
      if (.not. same_word(adaptive_cg%selector,'MAX_ANGLE')) then
         call reject('cg_selector supports MAX_ANGLE only',status,diagnostic); return
      end if
      if (.not. ieee_is_finite(adaptive_cg%refine_threshold) .or. &
          .not. ieee_is_finite(adaptive_cg%coarsen_threshold) .or. &
          adaptive_cg%coarsen_threshold < 0.0_dblprec .or. &
          adaptive_cg%refine_threshold > 2.0_dblprec .or. &
          adaptive_cg%coarsen_threshold > adaptive_cg%refine_threshold) then
         call reject('cg_coarsen_threshold/refine_threshold require 0 <= coarsen <= refine <= 2', &
            status,diagnostic); return
      end if
      if (adaptive_cg%update_interval <= 0) then
         call reject('cg_update_interval must be positive',status,diagnostic); return
      end if
      if (adaptive_cg%minimum_dwell_updates < 0) then
         call reject('cg_minimum_dwell_updates must be nonnegative',status,diagnostic); return
      end if
      if (adaptive_cg%buffer_blocks < 0) then
         call reject('cg_buffer_blocks must be nonnegative',status,diagnostic); return
      end if
      if (.not. same_word(adaptive_cg%channel_mode,'BASIS') .and. &
          .not. same_word(adaptive_cg%channel_mode,'FILE')) then
         call reject('cg_channel_mode must be BASIS or FILE',status,diagnostic); return
      end if
      if (same_word(adaptive_cg%channel_mode,'FILE') .and. &
          len_trim(adaptive_cg%channel_file) == 0) then
         call reject('cg_channel_file is required when cg_channel_mode=FILE',status,diagnostic); return
      end if
      if (.not. same_word(adaptive_cg%reconstruction,'ALIGNED') .and. &
          .not. same_word(adaptive_cg%reconstruction,'CONE')) then
         call reject('cg_reconstruction must be ALIGNED or CONE',status,diagnostic); return
      end if
      if (.not. ieee_is_finite(adaptive_cg%cone_angle_deg) .or. &
          adaptive_cg%cone_angle_deg < 0.0_dblprec .or. adaptive_cg%cone_angle_deg > 90.0_dblprec) then
         call reject('cg_cone_angle is in degrees and must lie in [0,90]',status,diagnostic); return
      end if
      if (.not. ieee_is_finite(adaptive_cg%energy_jump_limit_j) .or. &
          adaptive_cg%energy_jump_limit_j < 0.0_dblprec) then
         call reject('cg_energy_jump_limit is in joules and must be nonnegative',status,diagnostic); return
      end if
      if (adaptive_cg%diagnostics < 0 .or. adaptive_cg%diagnostics > 3) then
         call reject('cg_diagnostics must be 0, 1, 2, or 3',status,diagnostic); return
      end if
      if (mode /= 'S') then
         call reject('mode: adaptive coarse graining supports spin dynamics (S) only',status,diagnostic); return
      end if
      if (ipmode /= 'N') then
         call reject('ip_mode: adaptive coarse graining requires no initial phase; restart/MC/SD initial phases are unsupported', &
            status,diagnostic); return
      end if
      if (initmag == 4) then
         call reject('initmag=4 restart is unsupported until adaptive state serialization is implemented',status,diagnostic); return
      end if
      if (SDEalgh /= 1 .or. llg /= 1) then
         call reject('SDEalgh/llg: adaptive coarse graining supports deterministic fixed-length Heun LLG (1) only', &
            status,diagnostic); return
      end if
      if (Temp /= 0.0_dblprec .or. index('QRPT', do_qhb) > 0 .or. &
          index('YE', do_3tm) > 0) then
         call reject('Temp/do_qhb/do_3tm: adaptive coarse graining requires deterministic T=0 dynamics', &
            status,diagnostic); return
      end if
      if (any((/BC1,BC2,BC3/) /= 'P')) then
         call reject('BC: adaptive coarse graining currently requires P P P',status,diagnostic); return
      end if
      if (NA <= 0 .or. Natom /= NA*N1*N2*N3 .or. (Natom > 1 .and. NA == Natom)) then
         call reject('geometry: adaptive coarse graining requires a regular replicated cell, not an explicit device', &
            status,diagnostic); return
      end if
      if (any((/block_size_x,block_size_y,block_size_z/) <= 0) .or. &
          any(mod((/N1,N2,N3/),(/block_size_x,block_size_y,block_size_z/)) /= 0)) then
         call reject('block_size_x/y/z must be positive and divide ncell exactly; partial edge blocks are unsupported', &
            status,diagnostic); return
      end if
      if (do_ralloy /= 0 .or. Nchmax /= 1 .or. do_lsf /= 'N' .or. ind_mom_flag /= 'N') then
         call reject('do_ralloy/Nchmax/do_lsf/ind_mom_flag: only one fixed-length FM channel is supported', &
            status,diagnostic); return
      end if
      if (ham_inp%do_jtensor /= 0 .or. ham_inp%do_dm /= 0 .or. ham_inp%do_anisotropy /= 0 .or. &
          ham_inp%do_sa /= 0 .or. ham_inp%do_pd /= 0 .or. ham_inp%do_bq /= 0 .or. &
          ham_inp%do_biqdm /= 0 .or. ham_inp%do_ring /= 0 .or. ham_inp%do_chir /= 0 .or. &
          ham_inp%do_fourx /= 0) then
         call reject('Hamiltonian capability: production adaptive CG currently supports scalar Heisenberg exchange only', &
            status,diagnostic); return
      end if
      if (ham_inp%do_dip /= 0) then
         call reject('do_dip: legacy atomistic/macrocell dipole is unsupported; use GPU EWALD3D_FFT', &
            status,diagnostic); return
      end if
      if (trim(gpu_dipole_mode) /= 'OFF') then
         if (.not. (affirmative(do_gpu) .and. affirmative(do_gpu_llg))) then
            call reject('gpu_dipole_mode: adaptive FFT dipole requires do_gpu=Y and do_gpu_llg=Y', &
               status,diagnostic); return
         end if
         if (trim(gpu_dipole_mode) /= 'EWALD3D_FFT') then
            call reject('gpu_dipole_mode: adaptive periodic CG supports EWALD3D_FFT only', &
               status,diagnostic); return
         end if
      end if
      if (any(abs(hfield) > 0.0_dblprec) .or. do_bpulse /= 0 .or. demag == 'Y') then
         call reject('hfield/do_bpulse/demag: external or time-dependent fields are not supported by the first production CG path', &
            status,diagnostic); return
      end if
      if (do_sparse == 'Y' .or. do_reduced == 'Y' .or. mompar /= 0 .or. do_fixed_mom == 'Y') then
         call reject('do_sparse/do_reduced/mompar/do_fixed_mom is unsupported with adaptive coarse graining', &
            status,diagnostic); return
      end if
      if (plotenergy /= 0) then
         call reject('plotenergy: legacy atomistic energy output is unsupported; use cg_diagnostics',status,diagnostic); return
      end if
      if (do_autocorr == 'Y' .or. do_spintemp == 'Y' .or. do_kmc == 'Y') then
         call reject('do_autocorr/do_spintemp/do_kmc: post-step legacy samplers are unsupported with adaptive CG', &
            status,diagnostic); return
      end if
      if (.not. allocated(Landeg) .or. .not. allocated(lambda1_array) .or. &
          maxval(abs(Landeg-Landeg(1))) > 1.0d-12*max(1.0_dblprec,abs(Landeg(1))) .or. &
          maxval(abs(lambda1_array-lambda1_array(1))) > 1.0d-12) then
         call reject('Landeg/do_site_damping: the single-FM path requires uniform gamma and damping',status,diagnostic); return
      end if
   end subroutine validate_configuration

   subroutine adaptive_cg_cpu_step(step,atom_direction,atom_magnitude,atom_moment,status,diagnostic)
      integer, intent(in) :: step
      real(dblprec), intent(inout) :: atom_direction(:,:,:)
      real(dblprec), intent(in) :: atom_magnitude(:,:)
      real(dblprec), intent(out) :: atom_moment(:,:,:)
      integer, intent(out) :: status
      character(len=*), intent(out) :: diagnostic

      real(dblprec), allocatable :: atom0(:,:,:), coarse0(:,:,:,:), atom_predictor(:,:,:)
      real(dblprec), allocatable :: coarse_predictor(:,:,:,:), atom_field0(:,:,:), atom_field1(:,:,:)
      real(dblprec), allocatable :: coarse_field0(:,:,:,:), coarse_field1(:,:,:,:)
      real(dblprec), allocatable :: coarse_rhs0(:,:), coarse_rhs1(:,:)
      real(dblprec) :: rhs0(3), rhs1(3), candidate(3), norm, time0, time1
      integer :: ensemble, atom, block, local_status

      status = ADAPTIVE_CG_PRODUCTION_SETUP_FAILED
      diagnostic = ''
      if (.not. adaptive_cg_state%ready) then
         diagnostic = 'adaptive CPU step called without a ready production runtime'; return
      end if
      allocate(atom0,source=atom_direction)
      allocate(coarse0,source=adaptive_cg_state%coarse_direction)
      allocate(atom_predictor,source=atom_direction)
      allocate(coarse_predictor,source=adaptive_cg_state%coarse_direction)
      allocate(atom_field0(3,Natom,Mensemble),atom_field1(3,Natom,Mensemble))
      allocate(coarse_field0(3,1,adaptive_cg_state%topology%n_spatial_blocks,Mensemble))
      allocate(coarse_field1(3,1,adaptive_cg_state%topology%n_spatial_blocks,Mensemble))
      allocate(coarse_rhs0(3,adaptive_cg_state%topology%n_spatial_blocks), &
         coarse_rhs1(3,adaptive_cg_state%topology%n_spatial_blocks))

      call evaluate_all_ensembles(atom0,coarse0,atom_field0,coarse_field0,local_status,diagnostic)
      if (local_status /= STATIC_HYBRID_OK) return
      call cpu_time(time0)
      do ensemble = 1, Mensemble
         call coarse_llg_rhs(adaptive_cg_state%tensor,coarse0(:,1,:,ensemble), &
            coarse_field0(:,1,:,ensemble),coarse_rhs0,local_status,diagnostic)
         if (local_status /= COARSE_TENSOR_OK) return
         do atom = 1, Natom
            if (.not. adaptive_cg_state%runtime%hybrid%atomistic_atom(atom)) cycle
            call llg_rhs(atom0(:,atom,ensemble),atom_field0(:,atom,ensemble),Landeg(1), &
               lambda1_array(1),rhs0)
            candidate = atom0(:,atom,ensemble)+delta_t*rhs0
            atom_predictor(:,atom,ensemble) = candidate/sqrt(sum(candidate*candidate))
         end do
         do block = 1, adaptive_cg_state%topology%n_spatial_blocks
            if (.not. adaptive_cg_state%runtime%hybrid%coarse_block(block)) cycle
            candidate = coarse0(:,1,block,ensemble)+delta_t*coarse_rhs0(:,block)
            coarse_predictor(:,1,block,ensemble) = candidate/sqrt(sum(candidate*candidate))
         end do
      end do
      call cpu_time(time1)
      adaptive_cg_state%integration_seconds = adaptive_cg_state%integration_seconds + &
         max(0.0_dblprec,time1-time0)
      call evaluate_all_ensembles(atom_predictor,coarse_predictor,atom_field1,coarse_field1, &
         local_status,diagnostic)
      if (local_status /= STATIC_HYBRID_OK) return
      call cpu_time(time0)
      do ensemble = 1, Mensemble
         call coarse_llg_rhs(adaptive_cg_state%tensor,coarse0(:,1,:,ensemble), &
            coarse_field0(:,1,:,ensemble),coarse_rhs0,local_status,diagnostic)
         if (local_status /= COARSE_TENSOR_OK) return
         call coarse_llg_rhs(adaptive_cg_state%tensor,coarse_predictor(:,1,:,ensemble), &
            coarse_field1(:,1,:,ensemble),coarse_rhs1,local_status,diagnostic)
         if (local_status /= COARSE_TENSOR_OK) return
         do atom = 1, Natom
            if (.not. adaptive_cg_state%runtime%hybrid%atomistic_atom(atom)) cycle
            call llg_rhs(atom0(:,atom,ensemble),atom_field0(:,atom,ensemble),Landeg(1), &
               lambda1_array(1),rhs0)
            call llg_rhs(atom_predictor(:,atom,ensemble),atom_field1(:,atom,ensemble), &
               Landeg(1),lambda1_array(1),rhs1)
            candidate = atom0(:,atom,ensemble)+0.5_dblprec*delta_t*(rhs0+rhs1)
            atom_direction(:,atom,ensemble) = candidate/sqrt(sum(candidate*candidate))
         end do
         do block = 1, adaptive_cg_state%topology%n_spatial_blocks
            if (.not. adaptive_cg_state%runtime%hybrid%coarse_block(block)) cycle
            candidate = coarse0(:,1,block,ensemble)+0.5_dblprec*delta_t* &
               (coarse_rhs0(:,block)+coarse_rhs1(:,block))
            adaptive_cg_state%coarse_direction(:,1,block,ensemble) = &
               candidate/sqrt(sum(candidate*candidate))
         end do
      end do
      call cpu_time(time1)
      adaptive_cg_state%integration_seconds = adaptive_cg_state%integration_seconds + &
         max(0.0_dblprec,time1-time0)
      call cpu_time(time0)
      call reconstruct_coarse_atoms(atom_direction,status,diagnostic)
      if (status /= ADAPTIVE_CG_PRODUCTION_OK) return
      call cpu_time(time1)
      adaptive_cg_state%reconstruction_seconds = adaptive_cg_state%reconstruction_seconds + &
         max(0.0_dblprec,time1-time0)
      if (adaptive_cg_state%adaptive_mask) then
         call cpu_time(time0)
         call update_adaptive_mask(step,atom_direction,status,diagnostic)
         if (status /= ADAPTIVE_CG_PRODUCTION_OK) return
         call cpu_time(time1)
         adaptive_cg_state%selector_seconds = adaptive_cg_state%selector_seconds + &
            max(0.0_dblprec,time1-time0)
      end if
      do ensemble = 1, Mensemble
         do atom = 1, Natom
            atom_moment(:,atom,ensemble) = atom_magnitude(atom,ensemble)*atom_direction(:,atom,ensemble)
         end do
      end do
      adaptive_cg_state%completed_steps = adaptive_cg_state%completed_steps+1
      adaptive_cg_state%active_atom_updates = adaptive_cg_state%active_atom_updates + &
         int(size(adaptive_cg_state%runtime%active_atom_list),int64)
      adaptive_cg_state%active_block_updates = adaptive_cg_state%active_block_updates + &
         int(size(adaptive_cg_state%runtime%active_coarse_list),int64)
      status = ADAPTIVE_CG_PRODUCTION_OK
   end subroutine adaptive_cg_cpu_step

   subroutine evaluate_all_ensembles(atom_direction,coarse_direction,atom_field,coarse_field,status,diagnostic)
      real(dblprec), intent(in) :: atom_direction(:,:,:), coarse_direction(:,:,:,:)
      real(dblprec), intent(out) :: atom_field(:,:,:), coarse_field(:,:,:,:)
      integer, intent(out) :: status
      character(len=*), intent(out) :: diagnostic
      type(static_hybrid_energy_type) :: energy
      integer :: ensemble
      real(dblprec) :: time0, time1

      call cpu_time(time0)
      adaptive_cg_state%last_energy = static_hybrid_energy_type()
      do ensemble = 1, Mensemble
         call evaluate_static_hybrid_operator(adaptive_cg_state%runtime%hybrid, &
            atom_direction(:,:,ensemble),coarse_direction(:,:,:,ensemble), &
            adaptive_cg_state%bond_matrix_j,atom_field(:,:,ensemble), &
            coarse_field(:,:,:,ensemble),energy,status,diagnostic)
         if (status /= STATIC_HYBRID_OK) return
         adaptive_cg_state%field_evaluations = adaptive_cg_state%field_evaluations+1_int64
         adaptive_cg_state%last_energy%atomistic_bilinear_j = &
            adaptive_cg_state%last_energy%atomistic_bilinear_j + energy%atomistic_bilinear_j
         adaptive_cg_state%last_energy%atomistic_onsite_j = &
            adaptive_cg_state%last_energy%atomistic_onsite_j + energy%atomistic_onsite_j
         adaptive_cg_state%last_energy%coarse%exchange_j = &
            adaptive_cg_state%last_energy%coarse%exchange_j + energy%coarse%exchange_j
         adaptive_cg_state%last_energy%coarse%spiralization_j = &
            adaptive_cg_state%last_energy%coarse%spiralization_j + energy%coarse%spiralization_j
         adaptive_cg_state%last_energy%coarse%anisotropy_j = &
            adaptive_cg_state%last_energy%coarse%anisotropy_j + energy%coarse%anisotropy_j
         adaptive_cg_state%last_energy%coarse%external_j = &
            adaptive_cg_state%last_energy%coarse%external_j + energy%coarse%external_j
         adaptive_cg_state%last_energy%coarse%dipole_j = &
            adaptive_cg_state%last_energy%coarse%dipole_j + energy%coarse%dipole_j
         adaptive_cg_state%last_energy%total_j = &
            adaptive_cg_state%last_energy%total_j + energy%total_j
      end do
      adaptive_cg_state%last_atom_field_sum_t = sum(atom_field)
      adaptive_cg_state%last_atom_field_norm2_t2 = sum(atom_field*atom_field)
      adaptive_cg_state%last_coarse_field_sum_t = sum(coarse_field)
      adaptive_cg_state%last_coarse_field_norm2_t2 = sum(coarse_field*coarse_field)
      call cpu_time(time1)
      adaptive_cg_state%field_seconds = adaptive_cg_state%field_seconds + &
         max(0.0_dblprec,time1-time0)
   end subroutine evaluate_all_ensembles

   subroutine update_adaptive_mask(step,atom_direction,status,diagnostic)
      integer, intent(in) :: step
      real(dblprec), intent(inout) :: atom_direction(:,:,:)
      integer, intent(out) :: status
      character(len=*), intent(out) :: diagnostic
      type(selector_evaluation_type) :: evaluation
      type(selector_requests_type) :: requests
      integer(c_int) :: selector_status
      integer :: before, after

      if (mod(step,adaptive_cg%update_interval) /= 0) then
         status = ADAPTIVE_CG_PRODUCTION_OK; diagnostic=''; return
      end if
      call evaluate_selector_registry(adaptive_cg_state%selector_registry, &
         adaptive_cg_state%runtime%selector,atom_direction(:,:,1), &
         adaptive_cg_state%topology%atom_to_block,adaptive_cg_state%bond_atom(1,:), &
         adaptive_cg_state%bond_atom(2,:),evaluation,selector_status,diagnostic)
      if (selector_status /= SELECTOR_OK) then
         status = ADAPTIVE_CG_PRODUCTION_SETUP_FAILED; return
      end if
      call combine_selector_requests(evaluation,adaptive_cg_state%selector_configuration, &
         requests,selector_status,diagnostic)
      if (selector_status /= SELECTOR_OK) then
         status = ADAPTIVE_CG_PRODUCTION_SETUP_FAILED; return
      end if
      before = 0
      if (allocated(adaptive_cg_state%runtime%transition_log%event)) &
         before = size(adaptive_cg_state%runtime%transition_log%event)
      call apply_adaptive_transitions(adaptive_cg_state%runtime,adaptive_cg_state%topology, &
         requests,adaptive_cg_state%hard_fine_mask,evaluation, &
         adaptive_cg_state%selector_configuration,adaptive_cg_state%reconstruction, &
         int(step,c_int),ADAPTIVE_STAGE_COMPLETE_STEP,adaptive_cg_state%atom_moment_mub, &
         atom_direction,adaptive_cg_state%coarse_direction,production_energy_evaluator, &
         status,diagnostic)
      if (status /= ADAPTIVE_HYBRID_OK) then
         status = ADAPTIVE_CG_PRODUCTION_SETUP_FAILED; return
      end if
      after = before
      if (allocated(adaptive_cg_state%runtime%transition_log%event)) &
         after = size(adaptive_cg_state%runtime%transition_log%event)
      if (after > before) then
         adaptive_cg_state%accepted_transitions = adaptive_cg_state%accepted_transitions + &
            count(adaptive_cg_state%runtime%transition_log%event(before+1:after)%accepted)
         adaptive_cg_state%rejected_transitions = adaptive_cg_state%rejected_transitions + &
            count(.not. adaptive_cg_state%runtime%transition_log%event(before+1:after)%accepted)
         if (adaptive_cg%diagnostics >= 2) call print_transition_events(before+1,after)
      end if
      if (adaptive_cg%diagnostics >= 2) call print_resolution_state('step',step)
      status = ADAPTIVE_CG_PRODUCTION_OK
   end subroutine update_adaptive_mask

   subroutine production_energy_evaluator(operator,atom_direction,coarse_direction,energy_j,status,diagnostic)
      type(static_hybrid_operator_type), intent(in) :: operator
      real(dblprec), intent(in) :: atom_direction(:,:,:), coarse_direction(:,:,:,:)
      real(dblprec), intent(out) :: energy_j
      integer, intent(out) :: status
      character(len=*), intent(out) :: diagnostic
      real(dblprec) :: atom_field(3,Natom), coarse_field(3,1,operator%n_blocks)
      type(static_hybrid_energy_type) :: energy
      integer :: ensemble

      energy_j = 0.0_dblprec
      do ensemble = 1, size(atom_direction,3)
         call evaluate_static_hybrid_operator(operator,atom_direction(:,:,ensemble), &
            coarse_direction(:,:,:,ensemble),adaptive_cg_state%bond_matrix_j, &
            atom_field,coarse_field,energy,status,diagnostic)
         if (status /= STATIC_HYBRID_OK) return
         energy_j = energy_j+energy%total_j
      end do
   end subroutine production_energy_evaluator

   subroutine reconstruct_coarse_atoms(atom_direction,status,diagnostic)
      real(dblprec), intent(inout) :: atom_direction(:,:,:)
      integer, intent(out) :: status
      character(len=*), intent(out) :: diagnostic
      real(dblprec), allocatable :: requested(:,:,:)
      integer(int64), allocatable :: seeds(:,:)
      integer :: block, local_status, ensemble

      allocate(requested(3,1,Mensemble),seeds(1,Mensemble))
      do block = 1, adaptive_cg_state%topology%n_spatial_blocks
         if (.not. adaptive_cg_state%runtime%hybrid%coarse_block(block)) cycle
         do ensemble = 1, Mensemble
            requested(:,1,ensemble) = adaptive_cg_state%runtime%channel_moment_sum_mub(1,block,ensemble) * &
               adaptive_cg_state%coarse_direction(:,1,block,ensemble)
         end do
         if (adaptive_cg_state%reconstruction%scheme == RECONSTRUCTION_ALIGNED) then
            call reconstruct_block_aligned(adaptive_cg_state%topology,block, &
               adaptive_cg_state%atom_moment_mub,requested,atom_direction,1.0d-12, &
               local_status,diagnostic)
         else
            call reconstruct_block_constrained_cone(adaptive_cg_state%topology,block, &
               adaptive_cg_state%atom_moment_mub,requested,atom_direction, &
               adaptive_cg_state%reconstruction%cone_angle_rad,1.0d-12, &
               adaptive_cg_state%reconstruction%global_seed, &
               adaptive_cg_state%runtime%selector%transition_epoch(block),seeds, &
               local_status,diagnostic)
         end if
         if (local_status /= ADAPTIVE_HYBRID_OK) then
            status = ADAPTIVE_CG_PRODUCTION_SETUP_FAILED; return
         end if
      end do
      status = ADAPTIVE_CG_PRODUCTION_OK
      diagnostic = ''
   end subroutine reconstruct_coarse_atoms

   subroutine build_unique_bonds(status,diagnostic)
      integer, intent(out) :: status
      character(len=*), intent(out) :: diagnostic
      integer :: directed, atom, neighbour, target, iham, bond, found
      integer, allocatable :: pair_i(:), pair_j(:)
      real(dblprec), allocatable :: pair_matrix(:,:,:), pair_disp(:,:)
      real(dblprec) :: displacement(3), coefficient

      directed = sum(ham%nlistsize(ham%aham))
      allocate(pair_i(max(1,directed)),pair_j(max(1,directed)), &
         pair_matrix(3,3,max(1,directed)),pair_disp(3,max(1,directed)))
      pair_matrix = 0.0_dblprec
      pair_disp = 0.0_dblprec
      bond = 0
      do atom = 1, Natom
         iham = ham%aham(atom)
         do neighbour = 1, ham%nlistsize(iham)
            target = ham%nlist(neighbour,atom)
            if (target < 1 .or. target > Natom .or. target == atom) cycle
            found = 0
            if (bond > 0) then
               do found = 1, bond
                  if (pair_i(found) == min(atom,target) .and. pair_j(found) == max(atom,target)) exit
               end do
               if (found > bond) found = 0
            end if
            if (found == 0) then
               bond = bond+1
               found = bond
               pair_i(found) = min(atom,target)
               pair_j(found) = max(atom,target)
               call wrapped_displacement(pair_i(found),pair_j(found),displacement)
               pair_disp(:,found) = displacement
            end if
            coefficient = 0.5_dblprec*mub*mmom(atom,1)*mmom(target,1) * &
               ham%ncoup(neighbour,iham,1)
            pair_matrix(1,1,found) = pair_matrix(1,1,found)+coefficient
            pair_matrix(2,2,found) = pair_matrix(2,2,found)+coefficient
            pair_matrix(3,3,found) = pair_matrix(3,3,found)+coefficient
         end do
      end do
      if (bond == 0) then
         call setup_failed('Hamiltonian contains no usable scalar exchange bonds',status,diagnostic)
         return
      end if
      allocate(adaptive_cg_state%bond_atom(2,bond), &
         adaptive_cg_state%bond_matrix_j(3,3,bond), &
         adaptive_cg_state%bond_displacement_m(3,bond))
      adaptive_cg_state%bond_atom(1,:) = pair_i(1:bond)
      adaptive_cg_state%bond_atom(2,:) = pair_j(1:bond)
      adaptive_cg_state%bond_matrix_j = pair_matrix(:,:,1:bond)
      adaptive_cg_state%bond_displacement_m = pair_disp(:,1:bond)
      status = ADAPTIVE_CG_PRODUCTION_OK
      diagnostic = ''
   end subroutine build_unique_bonds

   subroutine read_mask_file(filename,mask,status,diagnostic)
      character(len=*), intent(in) :: filename
      logical, intent(inout) :: mask(:)
      integer, intent(out) :: status
      character(len=*), intent(out) :: diagnostic
      logical :: seen(size(mask))
      character(len=512) :: line
      character(len=32) :: state
      integer :: unit, ios, block, marker

      seen = .false.
      open(newunit=unit,file=filename,status='old',action='read',iostat=ios)
      if (ios /= 0) then
         call setup_failed('cg_static_mask_file cannot open '//trim(filename),status,diagnostic); return
      end if
      do
         read(unit,'(a)',iostat=ios) line
         if (ios < 0) exit
         if (ios > 0) then
            close(unit); call setup_failed('cg_static_mask_file read error',status,diagnostic); return
         end if
         call strip_comment(line)
         if (len_trim(line) == 0) cycle
         read(line,*,iostat=ios) block,state
         if (ios /= 0 .or. block < 1 .or. block > size(mask)) then
            close(unit); call setup_failed('cg_static_mask_file expects: one_based_block_id FINE|COARSE', &
               status,diagnostic); return
         end if
         if (seen(block)) then
            close(unit); call setup_failed('cg_static_mask_file contains duplicate block id',status,diagnostic); return
         end if
         if (same_word(state,'FINE')) then
            mask(block)=.true.
         else if (same_word(state,'COARSE')) then
            mask(block)=.false.
         else
            close(unit); call setup_failed('cg_static_mask_file state must be FINE or COARSE',status,diagnostic); return
         end if
         seen(block)=.true.
      end do
      close(unit)
      status = ADAPTIVE_CG_PRODUCTION_OK
      diagnostic = ''
   end subroutine read_mask_file

   subroutine read_channel_file(filename,map,status,diagnostic)
      character(len=*), intent(in) :: filename
      integer, intent(inout) :: map(:)
      integer, intent(out) :: status
      character(len=*), intent(out) :: diagnostic
      logical :: seen(size(map))
      character(len=512) :: line
      integer :: unit, ios, basis, channel

      seen=.false.
      open(newunit=unit,file=filename,status='old',action='read',iostat=ios)
      if (ios /= 0) then
         call setup_failed('cg_channel_file cannot open '//trim(filename),status,diagnostic); return
      end if
      do
         read(unit,'(a)',iostat=ios) line
         if (ios < 0) exit
         if (ios > 0) then
            close(unit); call setup_failed('cg_channel_file read error',status,diagnostic); return
         end if
         call strip_comment(line)
         if (len_trim(line)==0) cycle
         read(line,*,iostat=ios) basis,channel
         if (ios /= 0 .or. basis < 1 .or. basis > size(map) .or. &
             (channel /= 1 .and. channel /= -1)) then
            close(unit); call setup_failed('cg_channel_file expects: one_based_basis_id 1|-1',status,diagnostic); return
         end if
         if (seen(basis)) then
            close(unit); call setup_failed('cg_channel_file contains duplicate basis id',status,diagnostic); return
         end if
         map(basis)=channel
         seen(basis)=.true.
      end do
      close(unit)
      status=ADAPTIVE_CG_PRODUCTION_OK
      diagnostic=''
   end subroutine read_channel_file

   subroutine cleanup_adaptive_cg_production()
      adaptive_cg_state = adaptive_cg_production_state_type()
   end subroutine cleanup_adaptive_cg_production

   subroutine print_resolved_configuration()
      integer :: fine, coarse
      if (adaptive_cg%diagnostics == 0) return
      fine=count(adaptive_cg_state%runtime%hybrid%fine_seed)
      coarse=count(adaptive_cg_state%runtime%hybrid%coarse_block)
      write(*,'(a)') 'AdaptiveCG: capability accepted: regular periodic single-FM deterministic Heun'
      write(*,'(a,a,a,a)') 'AdaptiveCG: operator=',trim(adaptive_cg%operator), &
         ' mask_mode=',trim(adaptive_cg%mask_mode)
      write(*,'(a,a,a,a,a,i0)') 'AdaptiveCG: selector=',trim(adaptive_cg%selector), &
         ' reconstruction=',trim(adaptive_cg%reconstruction), &
         ' diagnostics=',adaptive_cg%diagnostics
      write(*,'(a,es24.16,a,es24.16,a,i0,a,i0,a,i0)') &
         'AdaptiveCG: refine_threshold=',adaptive_cg%refine_threshold, &
         ' coarsen_threshold=',adaptive_cg%coarsen_threshold, &
         ' update_interval=',adaptive_cg%update_interval, &
         ' minimum_dwell=',adaptive_cg%minimum_dwell_updates, &
         ' buffer_blocks=',adaptive_cg%buffer_blocks
      write(*,'(a,a,a,a)') 'AdaptiveCG: channel_mode=',trim(adaptive_cg%channel_mode), &
         ' fft_source_mapping=basis-resolved-to-single-dynamical-channel'
      write(*,'(a,a)') 'AdaptiveCG: gpu_dipole_mode=',trim(gpu_dipole_mode)
      write(*,'(a,3(i0,1x),a,3(i0,1x))') 'AdaptiveCG: block_size=',block_size_x,block_size_y,block_size_z, &
         ' block_grid=',adaptive_cg_state%topology%block_grid
      write(*,'(a,i0,a,i0,a,i0,a,i0)') 'AdaptiveCG: atoms=',Natom,' blocks=', &
         adaptive_cg_state%topology%n_spatial_blocks,' channels=1 initial_fine=',fine
      write(*,'(a,i0,a,i0,a,i0)') 'AdaptiveCG: initial_coarse=',coarse, &
         ' active_atoms=',size(adaptive_cg_state%runtime%active_atom_list), &
         ' active_blocks=',size(adaptive_cg_state%runtime%active_coarse_list)
      write(*,'(a,i0,a,i0)') 'AdaptiveCG: interface_atoms=', &
         size(adaptive_cg_state%runtime%interface_atom_list), &
         ' owned_cpu_bytes=',adaptive_owned_cpu_bytes()
      if (adaptive_cg%diagnostics >= 2) call print_resolution_state('initial',0)
   end subroutine print_resolved_configuration

   subroutine print_adaptive_cg_summary()
      integer(int64) :: baseline
      if (.not. adaptive_cg_state%ready .or. adaptive_cg%diagnostics == 0) return
      if (adaptive_cg_state%gpu_requested) then
         write(*,'(a)') 'AdaptiveCG: GPU lifecycle complete; device work and timing summary reported above'
         return
      end if
      baseline=int(adaptive_cg_state%completed_steps,int64)*int(Natom,int64)
      write(*,'(a,i0,a,i0,a,i0)') 'AdaptiveCG: completed_steps=',adaptive_cg_state%completed_steps, &
         ' field_evaluations=',adaptive_cg_state%field_evaluations, &
         ' active_atom_updates=',adaptive_cg_state%active_atom_updates
      write(*,'(a,i0,a,i0,a,i0)') 'AdaptiveCG: active_block_updates=', &
         adaptive_cg_state%active_block_updates,' accepted_transitions=', &
         adaptive_cg_state%accepted_transitions,' rejected_transitions=', &
         adaptive_cg_state%rejected_transitions
      write(*,'(a,i0,a,i0)') 'AdaptiveCG: baseline_atom_updates=',baseline, &
         ' reduced_atom_updates=',max(0_int64,baseline-adaptive_cg_state%active_atom_updates)
      write(*,'(a,es24.16,a,es24.16,a,es24.16,a,es24.16)') &
         'AdaptiveCG: last_energy_j atomistic_bilinear=', &
         adaptive_cg_state%last_energy%atomistic_bilinear_j, &
         ' atomistic_onsite=',adaptive_cg_state%last_energy%atomistic_onsite_j, &
         ' coarse_exchange=',adaptive_cg_state%last_energy%coarse%exchange_j, &
         ' coarse_spiralization=',adaptive_cg_state%last_energy%coarse%spiralization_j
      write(*,'(a,es24.16,a,es24.16,a,es24.16)') &
         'AdaptiveCG: last_energy_j coarse_anisotropy=', &
         adaptive_cg_state%last_energy%coarse%anisotropy_j, &
         ' coarse_external=',adaptive_cg_state%last_energy%coarse%external_j, &
         ' coarse_dipole=',adaptive_cg_state%last_energy%coarse%dipole_j
      write(*,'(a,es24.16)') 'AdaptiveCG: last_total_energy_j=', &
         adaptive_cg_state%last_energy%total_j
      write(*,'(a,es24.16,a,es24.16,a,es24.16,a,es24.16)') &
         'AdaptiveCG: last_field_checksums_t atom_sum=',adaptive_cg_state%last_atom_field_sum_t, &
         ' atom_norm2=',adaptive_cg_state%last_atom_field_norm2_t2, &
         ' coarse_sum=',adaptive_cg_state%last_coarse_field_sum_t, &
         ' coarse_norm2=',adaptive_cg_state%last_coarse_field_norm2_t2
      write(*,'(a,es24.16,a,es24.16,a,es24.16,a,es24.16)') &
         'AdaptiveCG: phase_seconds field=',adaptive_cg_state%field_seconds, &
         ' integration=',adaptive_cg_state%integration_seconds, &
         ' reconstruction=',adaptive_cg_state%reconstruction_seconds, &
         ' selector=',adaptive_cg_state%selector_seconds
      write(*,'(a,i0)') 'AdaptiveCG: owned_cpu_bytes=',adaptive_owned_cpu_bytes()
      call print_resolution_state('final',adaptive_cg_state%completed_steps)
      write(*,'(a,es24.16,a,es24.16)') 'AdaptiveCG: trajectory_checksums direction_sum=', &
         sum(emom),' direction_norm2=',sum(emom*emom)
   end subroutine print_adaptive_cg_summary

   subroutine print_transition_events(first,last)
      integer, intent(in) :: first, last
      integer :: event

      do event = first, last
         write(*,'(a,i0,a,i0,a,i0,a,i0,a,l1,a,a,a,a,a,3(es24.16,1x))') &
            'AdaptiveCG: transition step=', &
            adaptive_cg_state%runtime%transition_log%event(event)%synchronization_step, &
            ' block=',adaptive_cg_state%runtime%transition_log%event(event)%block, &
            ' old_state=',adaptive_cg_state%runtime%transition_log%event(event)%old_state, &
            ' new_state=',adaptive_cg_state%runtime%transition_log%event(event)%new_state, &
            ' accepted=',adaptive_cg_state%runtime%transition_log%event(event)%accepted, &
            ' reason=',trim(adaptive_cg_state%runtime%transition_log%event(event)%selector_reason), &
            ' outcome=',trim(adaptive_cg_state%runtime%transition_log%event(event)%outcome), &
            ' energies_j=',adaptive_cg_state%runtime%transition_log%event(event)%energy_before_j, &
            adaptive_cg_state%runtime%transition_log%event(event)%energy_after_j, &
            adaptive_cg_state%runtime%transition_log%event(event)%energy_jump_j
      end do
   end subroutine print_transition_events

   subroutine print_resolution_state(label,step)
      character(len=*), intent(in) :: label
      integer, intent(in) :: step
      integer :: block

      write(*,'(a,a,a,i0,a)',advance='no') 'AdaptiveCG: resolution_state label=', &
         trim(label),' step=',step,' values='
      do block = 1, adaptive_cg_state%topology%n_spatial_blocks
         if (block > 1) write(*,'(a)',advance='no') ','
         write(*,'(i0)',advance='no') adaptive_cg_state%runtime%hybrid%block_state(block)
      end do
      write(*,*)
      write(*,'(a,i0,a,i0,a,i0)') 'AdaptiveCG: resolution_counts fine=', &
         count(adaptive_cg_state%runtime%hybrid%block_state == 2), &
         ' interface=',count(adaptive_cg_state%runtime%hybrid%block_state == 1), &
         ' coarse=',count(adaptive_cg_state%runtime%hybrid%block_state == 0)
   end subroutine print_resolution_state

   integer(int64) function adaptive_owned_cpu_bytes()
      integer(int64) :: bytes

      bytes = 0_int64
      if (allocated(adaptive_cg_state%initial_fine_mask)) bytes=bytes+int(size( &
         adaptive_cg_state%initial_fine_mask),int64)*storage_size(.true.)/8
      if (allocated(adaptive_cg_state%hard_fine_mask)) bytes=bytes+int(size( &
         adaptive_cg_state%hard_fine_mask),int64)*storage_size(.true.)/8
      if (allocated(adaptive_cg_state%atom_moment_mub)) bytes=bytes+int(size( &
         adaptive_cg_state%atom_moment_mub),int64)*storage_size(1.0_dblprec)/8
      if (allocated(adaptive_cg_state%coarse_direction)) bytes=bytes+int(size( &
         adaptive_cg_state%coarse_direction),int64)*storage_size(1.0_dblprec)/8
      if (allocated(adaptive_cg_state%bond_atom)) bytes=bytes+int(size( &
         adaptive_cg_state%bond_atom),int64)*storage_size(1)/8
      if (allocated(adaptive_cg_state%bond_matrix_j)) bytes=bytes+int(size( &
         adaptive_cg_state%bond_matrix_j),int64)*storage_size(1.0_dblprec)/8
      if (allocated(adaptive_cg_state%bond_displacement_m)) bytes=bytes+int(size( &
         adaptive_cg_state%bond_displacement_m),int64)*storage_size(1.0_dblprec)/8
      if (allocated(adaptive_cg_state%runtime%active_atom_list)) bytes=bytes+int(size( &
         adaptive_cg_state%runtime%active_atom_list),int64)*storage_size(1)/8
      if (allocated(adaptive_cg_state%runtime%active_coarse_list)) bytes=bytes+int(size( &
         adaptive_cg_state%runtime%active_coarse_list),int64)*storage_size(1)/8
      if (allocated(adaptive_cg_state%runtime%interface_atom_list)) bytes=bytes+int(size( &
         adaptive_cg_state%runtime%interface_atom_list),int64)*storage_size(1)/8
      adaptive_owned_cpu_bytes = bytes
   end function adaptive_owned_cpu_bytes

   function block_mean_direction(block,ensemble) result(direction)
      integer, intent(in) :: block, ensemble
      real(dblprec) :: direction(3), norm
      integer :: position, atom
      direction=0.0_dblprec
      do position=adaptive_cg_state%topology%block_atom_offset(block)+1, &
            adaptive_cg_state%topology%block_atom_offset(block+1)
         atom=adaptive_cg_state%topology%block_atoms(position)
         if (adaptive_cg_state%topology%atom_to_dynamic_channel(atom)==1) &
            direction=direction+mmom(atom,ensemble)*emom(:,atom,ensemble)
      end do
      norm=sqrt(sum(direction*direction))
      if (norm <= tiny(norm)) direction=(/0.0_dblprec,0.0_dblprec,1.0_dblprec/)
      if (norm > tiny(norm)) direction=direction/norm
   end function block_mean_direction

   pure subroutine llg_rhs(direction,field,gamma,damping,rhs)
      real(dblprec), intent(in) :: direction(3),field(3),gamma,damping
      real(dblprec), intent(out) :: rhs(3)
      real(dblprec) :: cross1(3),cross2(3)
      cross1=cross(direction,field)
      cross2=cross(direction,cross1)
      rhs=-gamma*(cross1+damping*cross2)/(1.0_dblprec+damping*damping)
   end subroutine llg_rhs

   pure function cross(a,b) result(c)
      real(dblprec), intent(in) :: a(3),b(3)
      real(dblprec) :: c(3)
      c=(/a(2)*b(3)-a(3)*b(2),a(3)*b(1)-a(1)*b(3),a(1)*b(2)-a(2)*b(1)/)
   end function cross

   subroutine wrapped_displacement(i,j,value)
      integer, intent(in) :: i,j
      real(dblprec), intent(out) :: value(3)
      integer :: x,y,z
      real(dblprec) :: candidate(3)
      value=alat*(coord(:,j)-coord(:,i))
      do z=-1,1; do y=-1,1; do x=-1,1
         candidate=alat*(coord(:,j)-coord(:,i)+real(x*N1,dblprec)*C1+ &
            real(y*N2,dblprec)*C2+real(z*N3,dblprec)*C3)
         if (sum(candidate*candidate)<sum(value*value)) value=candidate
      end do; end do; end do
   end subroutine wrapped_displacement

   pure function inverse3(matrix) result(inverse)
      real(dblprec), intent(in) :: matrix(3,3)
      real(dblprec) :: inverse(3,3), determinant
      determinant=matrix(1,1)*(matrix(2,2)*matrix(3,3)-matrix(2,3)*matrix(3,2))- &
         matrix(1,2)*(matrix(2,1)*matrix(3,3)-matrix(2,3)*matrix(3,1))+ &
         matrix(1,3)*(matrix(2,1)*matrix(3,2)-matrix(2,2)*matrix(3,1))
      inverse(1,1)=(matrix(2,2)*matrix(3,3)-matrix(2,3)*matrix(3,2))/determinant
      inverse(1,2)=(matrix(1,3)*matrix(3,2)-matrix(1,2)*matrix(3,3))/determinant
      inverse(1,3)=(matrix(1,2)*matrix(2,3)-matrix(1,3)*matrix(2,2))/determinant
      inverse(2,1)=(matrix(2,3)*matrix(3,1)-matrix(2,1)*matrix(3,3))/determinant
      inverse(2,2)=(matrix(1,1)*matrix(3,3)-matrix(1,3)*matrix(3,1))/determinant
      inverse(2,3)=(matrix(1,3)*matrix(2,1)-matrix(1,1)*matrix(2,3))/determinant
      inverse(3,1)=(matrix(2,1)*matrix(3,2)-matrix(2,2)*matrix(3,1))/determinant
      inverse(3,2)=(matrix(1,2)*matrix(3,1)-matrix(1,1)*matrix(3,2))/determinant
      inverse(3,3)=(matrix(1,1)*matrix(2,2)-matrix(1,2)*matrix(2,1))/determinant
   end function inverse3

   pure logical function affirmative(value)
      character(len=*), intent(in) :: value
      affirmative=same_word(value,'Y')
   end function affirmative

   pure logical function same_word(left,right)
      character(len=*), intent(in) :: left,right
      character(len=len(left)) :: a
      character(len=len(right)) :: b
      integer :: i,code
      a=adjustl(left); b=adjustl(right)
      do i=1,len(a)
         code=iachar(a(i:i)); if(code>=iachar('a').and.code<=iachar('z')) a(i:i)=achar(code-32)
      end do
      do i=1,len(b)
         code=iachar(b(i:i)); if(code>=iachar('a').and.code<=iachar('z')) b(i:i)=achar(code-32)
      end do
      same_word=trim(a)==trim(b)
   end function same_word

   subroutine strip_comment(line)
      character(len=*), intent(inout) :: line
      integer :: i,p
      p=0
      do i=1,len_trim(line)
         if (index('#%!;',line(i:i))>0) then; p=i; exit; end if
      end do
      if (p>0) line(p:)=' '
   end subroutine strip_comment

   subroutine reject(message,status,diagnostic)
      character(len=*), intent(in) :: message
      integer, intent(out) :: status
      character(len=*), intent(out) :: diagnostic
      status=ADAPTIVE_CG_PRODUCTION_REJECTED
      diagnostic='AdaptiveCG setup rejected: '//trim(message)
   end subroutine reject

   subroutine setup_failed(message,status,diagnostic)
      character(len=*), intent(in) :: message
      integer, intent(out) :: status
      character(len=*), intent(out) :: diagnostic
      status=ADAPTIVE_CG_PRODUCTION_SETUP_FAILED
      diagnostic='AdaptiveCG setup failed: '//trim(message)
   end subroutine setup_failed

end module AdaptiveCGProduction
