!-------------------------------------------------------------------------------
! MODULE: AdaptiveHybridSolver
!> @brief Transactional adaptive-resolution updates for the CPU hybrid solver.
!>
!> Selector decisions are first evaluated on a private runtime copy.  Each
!> proposed block transition is then restricted/reconstructed, has its static
!> hybrid ownership rebuilt, and is passed to the caller's complete hybrid
!> energy evaluator.  Only a transition satisfying the configured energy-jump
!> limit is published.  Rejected candidates leave directions, selector epochs,
!> active lists, and interaction ownership unchanged.
!>
!> This is deterministic fixed-length reconstruction only.  It makes no
!> finite-temperature or local-thermalization claim.
!-------------------------------------------------------------------------------
module AdaptiveHybridSolver

   use, intrinsic :: iso_c_binding, only : c_int
   use, intrinsic :: iso_fortran_env, only : int64
   use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
   use Parameters, only : dblprec
   use BlockTopology, only : block_topology_type, REGULAR_REPLICATED_CELL
   use BlockSelector, only : selector_configuration_type, selector_requests_type, &
      selector_evaluation_type, block_selector_runtime_type, selector_transition_log_type, &
      selector_transition_event_type, &
      initialize_selector_runtime, advance_selector_state, SELECTOR_OK, &
      RESOLUTION_ATOMISTIC, RESOLUTION_COARSE
   use StaticHybridOperator, only : static_hybrid_operator_type, &
      rebuild_static_hybrid_ownership, STATIC_HYBRID_OK, STATIC_HYBRID_BUFFER

   implicit none
   private

   integer, parameter, public :: ADAPTIVE_HYBRID_OK = 0
   integer, parameter, public :: ADAPTIVE_HYBRID_INVALID_TOPOLOGY = 1
   integer, parameter, public :: ADAPTIVE_HYBRID_INVALID_STATE = 2
   integer, parameter, public :: ADAPTIVE_HYBRID_INVALID_CONFIGURATION = 3
   integer, parameter, public :: ADAPTIVE_HYBRID_INVALID_STAGE = 4
   integer, parameter, public :: ADAPTIVE_HYBRID_RECONSTRUCTION_FAILED = 5
   integer, parameter, public :: ADAPTIVE_HYBRID_ENERGY_FAILED = 6
   integer, parameter, public :: ADAPTIVE_HYBRID_ALLOCATION_FAILED = 7

   integer, parameter, public :: ADAPTIVE_STAGE_PREDICTOR = 1
   integer, parameter, public :: ADAPTIVE_STAGE_CORRECTOR = 2
   integer, parameter, public :: ADAPTIVE_STAGE_COMPLETE_STEP = 3

   integer, parameter, public :: RECONSTRUCTION_ALIGNED = 1
   integer, parameter, public :: RECONSTRUCTION_CONSTRAINED_CONE = 2

   type, public :: adaptive_reconstruction_configuration_type
      integer :: scheme = RECONSTRUCTION_ALIGNED
      real(dblprec) :: cone_angle_rad = 0.0_dblprec
      real(dblprec) :: resultant_tolerance = 1.0d-12
      real(dblprec) :: energy_jump_limit_j = huge(1.0_dblprec)
      integer(int64) :: global_seed = 1_int64
   end type adaptive_reconstruction_configuration_type

   type, public :: adaptive_transition_event_type
      integer(c_int) :: synchronization_step = 0_c_int
      integer(c_int) :: block = 0_c_int
      integer(c_int) :: old_state = RESOLUTION_COARSE
      integer(c_int) :: new_state = RESOLUTION_COARSE
      integer(c_int) :: previous_dwell_age = 0_c_int
      integer(c_int) :: transition_epoch = 0_c_int
      logical :: accepted = .false.
      character(len=48) :: selector_reason = ''
      character(len=48) :: outcome = ''
      integer :: reconstruction_scheme = RECONSTRUCTION_ALIGNED
      real(dblprec) :: energy_before_j = 0.0_dblprec
      real(dblprec) :: energy_after_j = 0.0_dblprec
      real(dblprec) :: energy_jump_j = 0.0_dblprec
      real(dblprec), allocatable :: criterion_score(:)
      integer(int64), allocatable :: reconstruction_seed(:,:) ! (channel,ensemble)
   end type adaptive_transition_event_type

   type, public :: adaptive_transition_log_type
      type(adaptive_transition_event_type), allocatable :: event(:)
   end type adaptive_transition_log_type

   type, public :: adaptive_hybrid_runtime_type
      logical :: ready = .false.
      type(block_selector_runtime_type) :: selector
      type(static_hybrid_operator_type) :: hybrid
      real(dblprec), allocatable :: coarse_resultant_mub(:,:,:,:) ! (3,channel,block,ensemble)
      real(dblprec), allocatable :: channel_moment_sum_mub(:,:,:) ! (channel,block,ensemble)
      integer, allocatable :: active_atom_list(:)
      integer, allocatable :: active_coarse_list(:)
      integer, allocatable :: interface_atom_list(:)
      type(adaptive_transition_log_type) :: transition_log
   end type adaptive_hybrid_runtime_type

   abstract interface
      subroutine adaptive_energy_evaluator(operator,atom_direction,coarse_direction, &
            energy_j,status,diagnostic)
         import :: static_hybrid_operator_type, dblprec
         type(static_hybrid_operator_type), intent(in) :: operator
         real(dblprec), intent(in) :: atom_direction(:,:,:)
         real(dblprec), intent(in) :: coarse_direction(:,:,:,:)
         real(dblprec), intent(out) :: energy_j
         integer, intent(out) :: status
         character(len=*), intent(out) :: diagnostic
      end subroutine adaptive_energy_evaluator
   end interface

   public :: restrict_channel_moments
   public :: reconstruct_block_aligned
   public :: reconstruct_block_constrained_cone
   public :: deterministic_reconstruction_seed
   public :: setup_adaptive_hybrid_runtime
   public :: apply_adaptive_transitions
   public :: clear_adaptive_transition_log

contains

   !> Channel-resolved R = sum_i mu_i s_i.  A zero resultant is explicitly
   !> marked undefined and is never normalized, so compensated channels remain
   !> independent and a compensated net vector cannot trigger NaNs.
   subroutine restrict_channel_moments(topology,atom_moment_mub,atom_direction, &
         resultant_mub,moment_sum_mub,direction,direction_defined,status,diagnostic)
      type(block_topology_type), intent(in) :: topology
      real(dblprec), intent(in) :: atom_moment_mub(:)
      real(dblprec), intent(in) :: atom_direction(:,:,:)
      real(dblprec), intent(out) :: resultant_mub(:,:,:,:)
      real(dblprec), intent(out) :: moment_sum_mub(:,:,:)
      real(dblprec), intent(out) :: direction(:,:,:,:)
      logical, intent(out) :: direction_defined(:,:,:)
      integer, intent(out) :: status
      character(len=*), intent(out) :: diagnostic

      integer :: atom, block, channel, ensemble
      real(dblprec) :: norm, scale

      status = ADAPTIVE_HYBRID_INVALID_STATE
      diagnostic = ''
      if (.not. topology%ready .or. topology%geometry_mode /= REGULAR_REPLICATED_CELL) then
         status = ADAPTIVE_HYBRID_INVALID_TOPOLOGY
         diagnostic = 'Moment restriction requires a ready regular block topology'
         return
      end if
      if (size(atom_moment_mub) /= topology%n_atoms .or. &
          any(shape(atom_direction) /= (/3,int(topology%n_atoms),size(atom_direction,3)/)) .or. &
          any(shape(resultant_mub) /= (/3,int(topology%n_dynamic_channels), &
             int(topology%n_spatial_blocks),size(atom_direction,3)/)) .or. &
          any(shape(moment_sum_mub) /= (/int(topology%n_dynamic_channels), &
             int(topology%n_spatial_blocks),size(atom_direction,3)/)) .or. &
          any(shape(direction) /= shape(resultant_mub)) .or. &
          any(shape(direction_defined) /= shape(moment_sum_mub))) then
         diagnostic = 'Restriction arrays do not match topology/channel/ensemble dimensions'
         return
      end if
      if (.not. all(ieee_is_finite(atom_moment_mub)) .or. &
          .not. all(ieee_is_finite(atom_direction)) .or. any(atom_moment_mub < 0.0_dblprec)) then
         diagnostic = 'Restriction requires finite directions and nonnegative moments'
         return
      end if

      resultant_mub = 0.0_dblprec
      moment_sum_mub = 0.0_dblprec
      direction = 0.0_dblprec
      direction_defined = .false.
      do atom = 1, topology%n_atoms
         channel = topology%atom_to_dynamic_channel(atom)
         if (channel <= 0) cycle
         if (atom_moment_mub(atom) <= 0.0_dblprec) then
            diagnostic = 'Every magnetic atom requires a positive moment for restriction'
            return
         end if
         block = topology%atom_to_block(atom)
         do ensemble = 1, size(atom_direction,3)
            resultant_mub(:,channel,block,ensemble) = &
               resultant_mub(:,channel,block,ensemble) + &
               atom_moment_mub(atom)*atom_direction(:,atom,ensemble)
            moment_sum_mub(channel,block,ensemble) = &
               moment_sum_mub(channel,block,ensemble) + atom_moment_mub(atom)
         end do
      end do
      do ensemble = 1, size(atom_direction,3)
         do block = 1, topology%n_spatial_blocks
            do channel = 1, topology%n_dynamic_channels
               norm = norm3(resultant_mub(:,channel,block,ensemble))
               scale = max(1.0_dblprec,moment_sum_mub(channel,block,ensemble))
               if (norm > 64.0_dblprec*epsilon(1.0_dblprec)*scale) then
                  direction(:,channel,block,ensemble) = &
                     resultant_mub(:,channel,block,ensemble)/norm
                  direction_defined(channel,block,ensemble) = .true.
               end if
            end do
         end do
      end do
      status = ADAPTIVE_HYBRID_OK
   end subroutine restrict_channel_moments

   subroutine reconstruct_block_aligned(topology,block,atom_moment_mub, &
         requested_resultant_mub,atom_direction,tolerance,status,diagnostic)
      type(block_topology_type), intent(in) :: topology
      integer, intent(in) :: block
      real(dblprec), intent(in) :: atom_moment_mub(:)
      real(dblprec), intent(in) :: requested_resultant_mub(:,:,:)
      real(dblprec), intent(inout) :: atom_direction(:,:,:)
      real(dblprec), intent(in) :: tolerance
      integer, intent(out) :: status
      character(len=*), intent(out) :: diagnostic

      integer :: atom, channel, ensemble
      real(dblprec) :: total, requested_norm, mean_direction(3)

      call validate_reconstruction_inputs(topology,block,atom_moment_mub, &
         requested_resultant_mub,atom_direction,tolerance,status,diagnostic)
      if (status /= ADAPTIVE_HYBRID_OK) return
      do ensemble = 1, size(atom_direction,3)
         do channel = 1, topology%n_dynamic_channels
            total = channel_block_moment(topology,block,channel,atom_moment_mub)
            requested_norm = norm3(requested_resultant_mub(:,channel,ensemble))
            if (requested_norm <= tolerance*max(1.0_dblprec,total) .or. &
                abs(requested_norm-total) > tolerance*max(1.0_dblprec,total)) then
               status = ADAPTIVE_HYBRID_RECONSTRUCTION_FAILED
               diagnostic = 'Aligned reconstruction requires the full nonzero channel moment'
               return
            end if
            mean_direction = requested_resultant_mub(:,channel,ensemble)/requested_norm
            do atom = 1, topology%n_atoms
               if (topology%atom_to_block(atom) == block .and. &
                   topology%atom_to_dynamic_channel(atom) == channel) then
                  atom_direction(:,atom,ensemble) = mean_direction
               end if
            end do
         end do
      end do
      status = ADAPTIVE_HYBRID_OK
   end subroutine reconstruct_block_aligned

   !> Deterministic constrained cone sampling.  Random transverse vectors are
   !> moment-weight de-biased first.  Scaling all of them by one bisection
   !> parameter preserves exactly zero transverse resultant while monotonically
   !> matching the requested longitudinal resultant.
   subroutine reconstruct_block_constrained_cone(topology,block,atom_moment_mub, &
         requested_resultant_mub,atom_direction,cone_angle_rad,tolerance, &
         global_seed,transition_epoch,seeds,status,diagnostic)
      type(block_topology_type), intent(in) :: topology
      integer, intent(in) :: block
      real(dblprec), intent(in) :: atom_moment_mub(:)
      real(dblprec), intent(in) :: requested_resultant_mub(:,:,:)
      real(dblprec), intent(inout) :: atom_direction(:,:,:)
      real(dblprec), intent(in) :: cone_angle_rad, tolerance
      integer(int64), intent(in) :: global_seed
      integer(c_int), intent(in) :: transition_epoch
      integer(int64), intent(out) :: seeds(:,:)
      integer, intent(out) :: status
      character(len=*), intent(out) :: diagnostic

      integer :: atom, channel, ensemble, i, nmember, iteration
      integer, allocatable :: member(:)
      integer(int64) :: rng
      real(dblprec), allocatable :: tx(:), ty(:), mu(:)
      real(dblprec) :: total, requested_norm, mean_direction(3), e1(3), e2(3)
      real(dblprec) :: angle, radius, mean_x, mean_y, maximum_radius
      real(dblprec) :: lower, upper, amplitude, longitudinal, residual(3)

      call validate_reconstruction_inputs(topology,block,atom_moment_mub, &
         requested_resultant_mub,atom_direction,tolerance,status,diagnostic)
      if (status /= ADAPTIVE_HYBRID_OK) return
      if (.not. ieee_is_finite(cone_angle_rad) .or. cone_angle_rad < 0.0_dblprec .or. &
          cone_angle_rad > 0.5_dblprec*acos(-1.0_dblprec) .or. &
          any(shape(seeds) /= (/int(topology%n_dynamic_channels),size(atom_direction,3)/))) then
         status = ADAPTIVE_HYBRID_INVALID_CONFIGURATION
         diagnostic = 'Cone angle must be in [0,pi/2] and seed storage must be (channel,ensemble)'
         return
      end if

      do channel = 1, topology%n_dynamic_channels
         nmember = count(topology%atom_to_block == block .and. &
            topology%atom_to_dynamic_channel == channel)
         allocate(member(nmember),tx(nmember),ty(nmember),mu(nmember))
         i = 0
         do atom = 1, topology%n_atoms
            if (topology%atom_to_block(atom) == block .and. &
                topology%atom_to_dynamic_channel(atom) == channel) then
               i = i+1
               member(i) = atom
               mu(i) = atom_moment_mub(atom)
            end if
         end do
         total = sum(mu)
         do ensemble = 1, size(atom_direction,3)
            requested_norm = norm3(requested_resultant_mub(:,channel,ensemble))
            if (requested_norm <= tolerance*max(1.0_dblprec,total) .or. &
                requested_norm > total+tolerance*max(1.0_dblprec,total)) then
               status = ADAPTIVE_HYBRID_RECONSTRUCTION_FAILED
               diagnostic = 'Requested cone resultant is zero or exceeds the channel moment sum'
               return
            end if
            mean_direction = requested_resultant_mub(:,channel,ensemble)/requested_norm
            call transverse_basis(mean_direction,e1,e2)
            seeds(channel,ensemble) = deterministic_reconstruction_seed( &
               global_seed,block,channel,ensemble,int(transition_epoch))
            rng = seeds(channel,ensemble)
            do i = 1, nmember
               angle = 2.0_dblprec*acos(-1.0_dblprec)*next_uniform(rng)
               radius = 0.5_dblprec+0.5_dblprec*next_uniform(rng)
               tx(i) = radius*cos(angle)
               ty(i) = radius*sin(angle)
            end do
            mean_x = sum(mu*tx)/total
            mean_y = sum(mu*ty)/total
            tx = tx-mean_x
            ty = ty-mean_y
            maximum_radius = maxval(sqrt(tx*tx+ty*ty))
            if (maximum_radius <= tiny(1.0_dblprec)) then
               if (abs(requested_norm-total) > tolerance*max(1.0_dblprec,total)) then
                  status = ADAPTIVE_HYBRID_RECONSTRUCTION_FAILED
                  diagnostic = 'A one-spin channel cannot represent a reduced cone resultant'
                  return
               end if
               tx = 0.0_dblprec
               ty = 0.0_dblprec
            else
               tx = tx/maximum_radius
               ty = ty/maximum_radius
            end if

            lower = 0.0_dblprec
            upper = sin(cone_angle_rad)
            longitudinal = cone_longitudinal(upper,tx,ty,mu)
            if (requested_norm < longitudinal-tolerance*max(1.0_dblprec,total)) then
               status = ADAPTIVE_HYBRID_RECONSTRUCTION_FAILED
               diagnostic = 'Configured cone angle cannot represent the requested resultant'
               return
            end if
            do iteration = 1, 100
               amplitude = 0.5_dblprec*(lower+upper)
               longitudinal = cone_longitudinal(amplitude,tx,ty,mu)
               if (longitudinal > requested_norm) then
                  lower = amplitude
               else
                  upper = amplitude
               end if
               if (abs(longitudinal-requested_norm) <= &
                   0.1_dblprec*tolerance*max(1.0_dblprec,total)) exit
            end do
            amplitude = 0.5_dblprec*(lower+upper)
            do i = 1, nmember
               radius = amplitude*sqrt(tx(i)*tx(i)+ty(i)*ty(i))
               atom_direction(:,member(i),ensemble) = &
                  sqrt(max(0.0_dblprec,1.0_dblprec-radius*radius))*mean_direction + &
                  amplitude*tx(i)*e1+amplitude*ty(i)*e2
            end do
            residual = -requested_resultant_mub(:,channel,ensemble)
            do i = 1, nmember
               residual = residual+mu(i)*atom_direction(:,member(i),ensemble)
            end do
            if (norm3(residual) > tolerance*max(1.0_dblprec,total) .or. &
                maxval(abs([(norm3(atom_direction(:,member(i),ensemble))-1.0_dblprec, &
                   i=1,nmember)])) > 16.0_dblprec*tolerance) then
               status = ADAPTIVE_HYBRID_RECONSTRUCTION_FAILED
               diagnostic = 'Constrained cone iteration did not match the requested resultant'
               return
            end if
         end do
         deallocate(member,tx,ty,mu)
      end do
      status = ADAPTIVE_HYBRID_OK
   end subroutine reconstruct_block_constrained_cone

   pure integer(int64) function deterministic_reconstruction_seed(global_seed,block, &
         channel,ensemble,epoch) result(seed)
      integer(int64), intent(in) :: global_seed
      integer, intent(in) :: block, channel, ensemble, epoch
      integer(int64), parameter :: modulus = 2147483647_int64

      seed = modulo(global_seed,modulus)
      seed = modulo(seed+104729_int64*modulo(int(block,int64),modulus),modulus)
      seed = modulo(seed+130363_int64*modulo(int(channel,int64),modulus),modulus)
      seed = modulo(seed+155921_int64*modulo(int(ensemble,int64),modulus),modulus)
      seed = modulo(seed+196613_int64*modulo(int(epoch,int64),modulus),modulus)
      if (seed <= 0_int64) seed = seed+modulus-1_int64
   end function deterministic_reconstruction_seed

   subroutine setup_adaptive_hybrid_runtime(runtime,topology,hybrid,atom_moment_mub, &
         atom_direction,coarse_direction,status,diagnostic)
      type(adaptive_hybrid_runtime_type), intent(out) :: runtime
      type(block_topology_type), intent(in) :: topology
      type(static_hybrid_operator_type), intent(in) :: hybrid
      real(dblprec), intent(in) :: atom_moment_mub(:)
      real(dblprec), intent(in) :: atom_direction(:,:,:)
      real(dblprec), intent(in) :: coarse_direction(:,:,:,:)
      integer, intent(out) :: status
      character(len=*), intent(out) :: diagnostic

      integer(c_int) :: selector_status
      integer :: block, channel, ensemble
      character(len=512) :: selector_diagnostic
      real(dblprec), allocatable :: restricted_direction(:,:,:,:)
      logical, allocatable :: defined(:,:,:)

      runtime%ready = .false.
      status = ADAPTIVE_HYBRID_INVALID_STATE
      diagnostic = ''
      if (.not. topology%ready .or. .not. hybrid%ready .or. &
          topology%n_atoms /= hybrid%n_atoms .or. &
          topology%n_spatial_blocks /= hybrid%n_blocks .or. &
          topology%n_dynamic_channels /= size(coarse_direction,2) .or. &
          any(shape(atom_direction) /= (/3,int(topology%n_atoms),size(atom_direction,3)/)) .or. &
          any(shape(coarse_direction) /= (/3,int(topology%n_dynamic_channels), &
             int(topology%n_spatial_blocks),size(atom_direction,3)/)) .or. &
          .not. all(ieee_is_finite(coarse_direction))) then
         diagnostic = 'Adaptive setup requires matching topology, hybrid, atom, and coarse states'
         return
      end if
      runtime%hybrid = hybrid
      call initialize_selector_runtime(runtime%selector,int(topology%n_spatial_blocks), &
         status=selector_status,diagnostic=selector_diagnostic)
      if (selector_status /= SELECTOR_OK) then
         diagnostic = trim(selector_diagnostic)
         return
      end if
      where (hybrid%fine_seed)
         runtime%selector%resolution_state = RESOLUTION_ATOMISTIC
      elsewhere
         runtime%selector%resolution_state = RESOLUTION_COARSE
      end where
      allocate(runtime%coarse_resultant_mub(3,topology%n_dynamic_channels, &
         topology%n_spatial_blocks,size(atom_direction,3)), &
         runtime%channel_moment_sum_mub(topology%n_dynamic_channels, &
         topology%n_spatial_blocks,size(atom_direction,3)), &
         restricted_direction(3,topology%n_dynamic_channels,topology%n_spatial_blocks, &
         size(atom_direction,3)),defined(topology%n_dynamic_channels, &
         topology%n_spatial_blocks,size(atom_direction,3)))
      call restrict_channel_moments(topology,atom_moment_mub,atom_direction, &
         runtime%coarse_resultant_mub,runtime%channel_moment_sum_mub, &
         restricted_direction,defined,status,diagnostic)
      if (status /= ADAPTIVE_HYBRID_OK) return
      do ensemble = 1, size(atom_direction,3)
         do block = 1, topology%n_spatial_blocks
            if (runtime%selector%resolution_state(block) /= RESOLUTION_COARSE) cycle
            do channel = 1, topology%n_dynamic_channels
               runtime%coarse_resultant_mub(:,channel,block,ensemble) = &
                  runtime%channel_moment_sum_mub(channel,block,ensemble)* &
                  coarse_direction(:,channel,block,ensemble)
            end do
         end do
      end do
      call rebuild_active_lists(runtime,topology,status,diagnostic)
      if (status /= ADAPTIVE_HYBRID_OK) return
      runtime%ready = .true.
   end subroutine setup_adaptive_hybrid_runtime

   !> Apply selector transitions only at a complete-step synchronization point.
   !> The energy evaluator must call the complete static hybrid energy path for
   !> every ensemble and return their total.
   subroutine apply_adaptive_transitions(runtime,topology,requests,hard_atomistic_mask, &
         evaluation,selector_configuration,reconstruction_configuration, &
         synchronization_step,integration_stage,atom_moment_mub,atom_direction, &
         coarse_direction,energy_evaluator,status,diagnostic)
      type(adaptive_hybrid_runtime_type), intent(inout) :: runtime
      type(block_topology_type), intent(in) :: topology
      type(selector_requests_type), intent(in) :: requests
      logical, intent(in) :: hard_atomistic_mask(:)
      type(selector_evaluation_type), intent(in) :: evaluation
      type(selector_configuration_type), intent(in) :: selector_configuration
      type(adaptive_reconstruction_configuration_type), intent(in) :: reconstruction_configuration
      integer(c_int), intent(in) :: synchronization_step
      integer, intent(in) :: integration_stage
      real(dblprec), intent(in) :: atom_moment_mub(:)
      real(dblprec), intent(inout) :: atom_direction(:,:,:)
      real(dblprec), intent(inout) :: coarse_direction(:,:,:,:)
      procedure(adaptive_energy_evaluator) :: energy_evaluator
      integer, intent(out) :: status
      character(len=*), intent(out) :: diagnostic

      type(block_selector_runtime_type) :: decision, working
      type(selector_transition_log_type) :: decisions
      type(static_hybrid_operator_type) :: candidate_hybrid
      real(dblprec), allocatable :: candidate_atom(:,:,:), candidate_coarse(:,:,:,:)
      real(dblprec), allocatable :: resultants(:,:,:,:), sums(:,:,:), directions(:,:,:,:)
      logical, allocatable :: defined(:,:,:), has_decision(:)
      integer, allocatable :: candidate_active_atom(:), candidate_active_coarse(:)
      integer, allocatable :: candidate_interface_atom(:)
      integer(int64), allocatable :: seeds(:,:)
      integer(c_int), allocatable :: old_age(:), old_epoch(:), old_state(:)
      integer(c_int) :: selector_status
      integer :: event, block, channel, ensemble, dependency_status
      real(dblprec) :: energy_before, energy_after, energy_jump
      logical :: accepted
      character(len=512) :: dependency_diagnostic, outcome

      status = ADAPTIVE_HYBRID_INVALID_STATE
      diagnostic = ''
      if (.not. runtime%ready) then
         diagnostic = 'Adaptive hybrid runtime is not initialized'
         return
      end if
      if (.not. topology%ready .or. topology%n_atoms /= runtime%hybrid%n_atoms .or. &
          topology%n_spatial_blocks /= runtime%hybrid%n_blocks .or. &
          size(atom_moment_mub) /= topology%n_atoms .or. &
          any(shape(atom_direction) /= (/3,int(topology%n_atoms),size(atom_direction,3)/)) .or. &
          any(shape(coarse_direction) /= (/3,int(topology%n_dynamic_channels), &
             int(topology%n_spatial_blocks),size(atom_direction,3)/)) .or. &
          any(shape(runtime%coarse_resultant_mub) /= shape(coarse_direction)) .or. &
          any(shape(runtime%channel_moment_sum_mub) /= &
             (/int(topology%n_dynamic_channels),int(topology%n_spatial_blocks), &
             size(atom_direction,3)/)) .or. .not. all(ieee_is_finite(atom_moment_mub)) .or. &
          .not. all(ieee_is_finite(atom_direction)) .or. &
          .not. all(ieee_is_finite(coarse_direction))) then
         diagnostic = 'Adaptive transition state does not match the initialized topology and ensembles'
         return
      end if
      if (integration_stage /= ADAPTIVE_STAGE_COMPLETE_STEP) then
         status = ADAPTIVE_HYBRID_INVALID_STAGE
         diagnostic = 'Adaptive transitions are allowed only between complete integration steps'
         return
      end if
      if (.not. valid_reconstruction_configuration(reconstruction_configuration)) then
         status = ADAPTIVE_HYBRID_INVALID_CONFIGURATION
         diagnostic = 'Adaptive reconstruction scheme, tolerance, cone angle, or energy limit is invalid'
         return
      end if

      decision = runtime%selector
      old_state = runtime%selector%resolution_state
      old_age = runtime%selector%state_age
      old_epoch = runtime%selector%transition_epoch
      call advance_selector_state(decision,requests,hard_atomistic_mask,evaluation, &
         selector_configuration,synchronization_step,decisions,selector_status, &
         dependency_diagnostic)
      if (selector_status /= SELECTOR_OK) then
         diagnostic = 'Selector transition proposal failed: '//trim(dependency_diagnostic)
         return
      end if
      if (.not. allocated(decisions%event)) then
         runtime%selector = decision
         status = ADAPTIVE_HYBRID_OK
         return
      end if

      call energy_evaluator(runtime%hybrid,atom_direction,coarse_direction, &
         energy_before,dependency_status,dependency_diagnostic)
      if (dependency_status /= 0 .or. .not. ieee_is_finite(energy_before)) then
         status = ADAPTIVE_HYBRID_ENERGY_FAILED
         diagnostic = 'Current hybrid energy evaluation failed: '//trim(dependency_diagnostic)
         return
      end if

      working = runtime%selector
      allocate(has_decision(topology%n_spatial_blocks),seeds(topology%n_dynamic_channels, &
         size(atom_direction,3)),resultants(3,topology%n_dynamic_channels, &
         topology%n_spatial_blocks,size(atom_direction,3)), &
         sums(topology%n_dynamic_channels,topology%n_spatial_blocks,size(atom_direction,3)), &
         directions(3,topology%n_dynamic_channels,topology%n_spatial_blocks, &
         size(atom_direction,3)),defined(topology%n_dynamic_channels, &
         topology%n_spatial_blocks,size(atom_direction,3)))
      has_decision = .false.

      do event = 1, size(decisions%event)
         status = ADAPTIVE_HYBRID_OK
         block = decisions%event(event)%block
         has_decision(block) = .true.
         candidate_atom = atom_direction
         candidate_coarse = coarse_direction
         candidate_hybrid = runtime%hybrid
         accepted = .false.
         outcome = ''
         seeds = 0_int64

         if (decisions%event(event)%new_state == RESOLUTION_COARSE) then
            call restrict_channel_moments(topology,atom_moment_mub,candidate_atom, &
               resultants,sums,directions,defined,status,dependency_diagnostic)
            if (status == ADAPTIVE_HYBRID_OK .and. all(defined(:,block,:))) then
               do ensemble = 1, size(atom_direction,3)
                  do channel = 1, topology%n_dynamic_channels
                     candidate_coarse(:,channel,block,ensemble) = &
                        directions(:,channel,block,ensemble)
                  end do
               end do
            else
               status = ADAPTIVE_HYBRID_RECONSTRUCTION_FAILED
               outcome = 'undefined-restricted-resultant'
            end if
         else
            select case (reconstruction_configuration%scheme)
            case (RECONSTRUCTION_ALIGNED)
               call reconstruct_block_aligned(topology,block,atom_moment_mub, &
                  runtime%coarse_resultant_mub(:,:,block,:),candidate_atom, &
                  reconstruction_configuration%resultant_tolerance,status, &
                  dependency_diagnostic)
            case (RECONSTRUCTION_CONSTRAINED_CONE)
               call reconstruct_block_constrained_cone(topology,block,atom_moment_mub, &
                  runtime%coarse_resultant_mub(:,:,block,:),candidate_atom, &
                  reconstruction_configuration%cone_angle_rad, &
                  reconstruction_configuration%resultant_tolerance, &
                  reconstruction_configuration%global_seed,decision%transition_epoch(block), &
                  seeds,status,dependency_diagnostic)
            end select
            if (status /= ADAPTIVE_HYBRID_OK) outcome = 'reconstruction-rejected'
         end if

         if (status == ADAPTIVE_HYBRID_OK) then
            working%resolution_state(block) = decisions%event(event)%new_state
            call rebuild_static_hybrid_ownership(candidate_hybrid,topology, &
               working%resolution_state == RESOLUTION_ATOMISTIC,dependency_status, &
               dependency_diagnostic)
            if (dependency_status /= STATIC_HYBRID_OK) then
               working%resolution_state(block) = decisions%event(event)%old_state
               status = ADAPTIVE_HYBRID_INVALID_STATE
               outcome = 'ownership-rebuild-failed'
            end if
         end if

         energy_after = energy_before
         energy_jump = 0.0_dblprec
         if (status == ADAPTIVE_HYBRID_OK) then
            call energy_evaluator(candidate_hybrid,candidate_atom,candidate_coarse, &
               energy_after,dependency_status,dependency_diagnostic)
            if (dependency_status /= 0 .or. .not. ieee_is_finite(energy_after)) then
               working%resolution_state(block) = decisions%event(event)%old_state
               status = ADAPTIVE_HYBRID_ENERGY_FAILED
               outcome = 'energy-evaluation-failed'
            else
               energy_jump = energy_after-energy_before
               accepted = abs(energy_jump) <= reconstruction_configuration%energy_jump_limit_j
               if (accepted) then
                  outcome = 'accepted'
               else
                  working%resolution_state(block) = decisions%event(event)%old_state
                  outcome = 'energy-jump-rejected'
               end if
            end if
         end if
         if (accepted) then
            call make_active_lists(candidate_hybrid,topology,candidate_active_atom, &
               candidate_active_coarse,candidate_interface_atom,dependency_status, &
               dependency_diagnostic)
            if (dependency_status /= ADAPTIVE_HYBRID_OK) then
               accepted = .false.
               status = dependency_status
               outcome = 'active-list-rebuild-failed'
            end if
         end if

         call append_adaptive_event(runtime%transition_log,decisions%event(event), &
            old_age(block),decision%transition_epoch(block),accepted,trim(outcome), &
            reconstruction_configuration%scheme,energy_before,energy_after,energy_jump,seeds)

         if (accepted) then
            atom_direction = candidate_atom
            coarse_direction = candidate_coarse
            runtime%hybrid = candidate_hybrid
            call move_alloc(candidate_active_atom,runtime%active_atom_list)
            call move_alloc(candidate_active_coarse,runtime%active_coarse_list)
            call move_alloc(candidate_interface_atom,runtime%interface_atom_list)
            working%state_age(block) = decision%state_age(block)
            working%transition_epoch(block) = decision%transition_epoch(block)
            if (decisions%event(event)%new_state == RESOLUTION_COARSE) then
               runtime%coarse_resultant_mub(:,:,block,:) = resultants(:,:,block,:)
               runtime%channel_moment_sum_mub(:,block,:) = sums(:,block,:)
            end if
            energy_before = energy_after
         else
            working%resolution_state(block) = old_state(block)
            working%state_age(block) = old_age(block)
            working%transition_epoch(block) = old_epoch(block)
         end if
      end do

      do block = 1, topology%n_spatial_blocks
         if (.not. has_decision(block)) working%state_age(block) = decision%state_age(block)
      end do
      runtime%selector = working
      status = ADAPTIVE_HYBRID_OK
   end subroutine apply_adaptive_transitions

   subroutine clear_adaptive_transition_log(log)
      type(adaptive_transition_log_type), intent(inout) :: log
      if (allocated(log%event)) deallocate(log%event)
   end subroutine clear_adaptive_transition_log

   subroutine rebuild_active_lists(runtime,topology,status,diagnostic)
      type(adaptive_hybrid_runtime_type), intent(inout) :: runtime
      type(block_topology_type), intent(in) :: topology
      integer, intent(out) :: status
      character(len=*), intent(out) :: diagnostic
      integer, allocatable :: active_atom(:), active_coarse(:), interface_atom(:)

      call make_active_lists(runtime%hybrid,topology,active_atom,active_coarse, &
         interface_atom,status,diagnostic)
      if (status /= ADAPTIVE_HYBRID_OK) return
      call move_alloc(active_atom,runtime%active_atom_list)
      call move_alloc(active_coarse,runtime%active_coarse_list)
      call move_alloc(interface_atom,runtime%interface_atom_list)
   end subroutine rebuild_active_lists

   subroutine make_active_lists(operator,topology,active_atom,active_coarse, &
         interface_atom,status,diagnostic)
      type(static_hybrid_operator_type), intent(in) :: operator
      type(block_topology_type), intent(in) :: topology
      integer, allocatable, intent(out) :: active_atom(:), active_coarse(:)
      integer, allocatable, intent(out) :: interface_atom(:)
      integer, intent(out) :: status
      character(len=*), intent(out) :: diagnostic
      integer :: atom, block, na, nc, ni, allocation_status

      status = ADAPTIVE_HYBRID_INVALID_STATE
      diagnostic = ''
      if (.not. operator%ready .or. &
          size(operator%atomistic_atom) /= topology%n_atoms) then
         diagnostic = 'Cannot build active lists from mismatched hybrid ownership'
         return
      end if
      na = count(operator%atomistic_atom)
      nc = count(operator%coarse_block)
      ni = 0
      do atom = 1, topology%n_atoms
         block = topology%atom_to_block(atom)
         if (operator%block_state(block) == STATIC_HYBRID_BUFFER) ni = ni+1
      end do
      allocate(active_atom(na),active_coarse(nc),interface_atom(ni),stat=allocation_status)
      if (allocation_status /= 0) then
         status = ADAPTIVE_HYBRID_ALLOCATION_FAILED
         diagnostic = 'Unable to allocate rebuilt adaptive active lists'
         return
      end if
      na = 0
      nc = 0
      ni = 0
      do atom = 1, topology%n_atoms
         block = topology%atom_to_block(atom)
         if (operator%atomistic_atom(atom)) then
            na = na+1
            active_atom(na) = atom
         end if
         if (operator%block_state(block) == STATIC_HYBRID_BUFFER) then
            ni = ni+1
            interface_atom(ni) = atom
         end if
      end do
      do block = 1, topology%n_spatial_blocks
         if (operator%coarse_block(block)) then
            nc = nc+1
            active_coarse(nc) = block
         end if
      end do
      status = ADAPTIVE_HYBRID_OK
   end subroutine make_active_lists

   subroutine validate_reconstruction_inputs(topology,block,atom_moment_mub, &
         requested_resultant_mub,atom_direction,tolerance,status,diagnostic)
      type(block_topology_type), intent(in) :: topology
      integer, intent(in) :: block
      real(dblprec), intent(in) :: atom_moment_mub(:), requested_resultant_mub(:,:,:)
      real(dblprec), intent(in) :: atom_direction(:,:,:), tolerance
      integer, intent(out) :: status
      character(len=*), intent(out) :: diagnostic

      status = ADAPTIVE_HYBRID_INVALID_STATE
      diagnostic = ''
      if (.not. topology%ready .or. block < 1 .or. block > topology%n_spatial_blocks .or. &
          size(atom_moment_mub) /= topology%n_atoms .or. &
          any(shape(atom_direction) /= (/3,int(topology%n_atoms),size(atom_direction,3)/)) .or. &
          any(shape(requested_resultant_mub) /= (/3,int(topology%n_dynamic_channels), &
             size(atom_direction,3)/)) .or. .not. ieee_is_finite(tolerance) .or. &
          tolerance <= 0.0_dblprec .or. .not. all(ieee_is_finite(atom_moment_mub)) .or. &
          .not. all(ieee_is_finite(requested_resultant_mub)) .or. &
          any(atom_moment_mub < 0.0_dblprec)) then
         diagnostic = 'Reconstruction inputs do not match the topology or contain invalid values'
         return
      end if
      status = ADAPTIVE_HYBRID_OK
   end subroutine validate_reconstruction_inputs

   pure real(dblprec) function channel_block_moment(topology,block,channel,atom_moment_mub)
      type(block_topology_type), intent(in) :: topology
      integer, intent(in) :: block, channel
      real(dblprec), intent(in) :: atom_moment_mub(:)
      integer :: atom

      channel_block_moment = 0.0_dblprec
      do atom = 1, topology%n_atoms
         if (topology%atom_to_block(atom) == block .and. &
             topology%atom_to_dynamic_channel(atom) == channel) then
            channel_block_moment = channel_block_moment+atom_moment_mub(atom)
         end if
      end do
   end function channel_block_moment

   pure real(dblprec) function cone_longitudinal(amplitude,tx,ty,mu)
      real(dblprec), intent(in) :: amplitude, tx(:), ty(:), mu(:)
      cone_longitudinal = sum(mu*sqrt(max(0.0_dblprec, &
         1.0_dblprec-amplitude*amplitude*(tx*tx+ty*ty))))
   end function cone_longitudinal

   subroutine transverse_basis(direction,e1,e2)
      real(dblprec), intent(in) :: direction(3)
      real(dblprec), intent(out) :: e1(3), e2(3)
      real(dblprec) :: reference(3)

      if (abs(direction(3)) < 0.8_dblprec) then
         reference = (/0.0_dblprec,0.0_dblprec,1.0_dblprec/)
      else
         reference = (/0.0_dblprec,1.0_dblprec,0.0_dblprec/)
      end if
      e1 = cross(reference,direction)
      e1 = e1/norm3(e1)
      e2 = cross(direction,e1)
   end subroutine transverse_basis

   real(dblprec) function next_uniform(state)
      integer(int64), intent(inout) :: state
      integer(int64), parameter :: modulus = 2147483647_int64
      state = modulo(48271_int64*state,modulus)
      if (state <= 0_int64) state = 1_int64
      next_uniform = real(state,dblprec)/real(modulus,dblprec)
   end function next_uniform

   pure function cross(a,b) result(c)
      real(dblprec), intent(in) :: a(3), b(3)
      real(dblprec) :: c(3)
      c = (/a(2)*b(3)-a(3)*b(2),a(3)*b(1)-a(1)*b(3), &
         a(1)*b(2)-a(2)*b(1)/)
   end function cross

   pure real(dblprec) function norm3(vector)
      real(dblprec), intent(in) :: vector(3)
      norm3 = sqrt(sum(vector*vector))
   end function norm3

   logical function valid_reconstruction_configuration(configuration)
      type(adaptive_reconstruction_configuration_type), intent(in) :: configuration
      valid_reconstruction_configuration = &
         (configuration%scheme == RECONSTRUCTION_ALIGNED .or. &
          configuration%scheme == RECONSTRUCTION_CONSTRAINED_CONE) .and. &
         ieee_is_finite(configuration%cone_angle_rad) .and. &
         configuration%cone_angle_rad >= 0.0_dblprec .and. &
         configuration%cone_angle_rad <= 0.5_dblprec*acos(-1.0_dblprec) .and. &
         ieee_is_finite(configuration%resultant_tolerance) .and. &
         configuration%resultant_tolerance > 0.0_dblprec .and. &
         ieee_is_finite(configuration%energy_jump_limit_j) .and. &
         configuration%energy_jump_limit_j >= 0.0_dblprec
   end function valid_reconstruction_configuration

   subroutine append_adaptive_event(log,decision,previous_age,epoch,accepted,outcome, &
         scheme,energy_before,energy_after,energy_jump,seeds)
      type(adaptive_transition_log_type), intent(inout) :: log
      type(selector_transition_event_type), intent(in) :: decision
      integer(c_int), intent(in) :: previous_age, epoch
      logical, intent(in) :: accepted
      character(len=*), intent(in) :: outcome
      integer, intent(in) :: scheme
      real(dblprec), intent(in) :: energy_before, energy_after, energy_jump
      integer(int64), intent(in) :: seeds(:,:)
      type(adaptive_transition_event_type), allocatable :: expanded(:)
      integer :: old_count

      old_count = 0
      if (allocated(log%event)) old_count = size(log%event)
      allocate(expanded(old_count+1))
      if (old_count > 0) expanded(:old_count) = log%event
      expanded(old_count+1)%synchronization_step = decision%synchronization_step
      expanded(old_count+1)%block = decision%block
      expanded(old_count+1)%old_state = decision%old_state
      expanded(old_count+1)%new_state = decision%new_state
      expanded(old_count+1)%previous_dwell_age = previous_age
      expanded(old_count+1)%transition_epoch = epoch
      expanded(old_count+1)%accepted = accepted
      expanded(old_count+1)%selector_reason = decision%reason
      expanded(old_count+1)%outcome = outcome
      expanded(old_count+1)%reconstruction_scheme = scheme
      expanded(old_count+1)%energy_before_j = energy_before
      expanded(old_count+1)%energy_after_j = energy_after
      expanded(old_count+1)%energy_jump_j = energy_jump
      allocate(expanded(old_count+1)%criterion_score(size(decision%score)), &
         expanded(old_count+1)%reconstruction_seed(size(seeds,1),size(seeds,2)))
      expanded(old_count+1)%criterion_score = decision%score
      expanded(old_count+1)%reconstruction_seed = seeds
      call move_alloc(expanded,log%event)
   end subroutine append_adaptive_event

end module AdaptiveHybridSolver
