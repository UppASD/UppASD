program test_adaptive_hybrid_solver

   use, intrinsic :: iso_c_binding, only : c_double, c_int
   use, intrinsic :: iso_fortran_env, only : int64
   use Parameters, only : dblprec
   use BlockTopology
   use BlockSelector
   use stiffness, only : coarse_material_type
   use CoarseTensorOperator
   use SmoothProjectedOperator
   use StaticHybridOperator
   use AdaptiveHybridSolver

   implicit none

   integer :: failures
   real(dblprec), parameter :: pi = acos(-1.0_dblprec)
   real(dblprec), allocatable :: evaluator_matrix(:,:,:)
   real(dblprec) :: evaluator_mask_penalty

   failures = 0
   evaluator_mask_penalty = 0.0_dblprec
   call test_channel_restriction_and_reconstruction()
   call test_transactional_chain_transitions()
   call test_adaptive_skyrmion_transition()
   call test_polarization_forces_refine_of_coarse_block()
   call test_transition_workspace_lifecycle()

   if (failures /= 0) then
      write(*,'(a,i0)') 'adaptive hybrid solver tests failed: ',failures
      stop 1
   end if
   write(*,'(a)') 'adaptive hybrid solver tests passed'

contains

   subroutine test_channel_restriction_and_reconstruction()
      type(block_topology_type) :: topology
      integer :: atom, channel, ensemble, status
      integer(c_int) :: topology_status
      character(len=512) :: message
      real(dblprec) :: cell(3,3), moment(8), spin(3,8,2)
      real(dblprec) :: resultant(3,2,1,2), moment_sum(2,1,2)
      real(dblprec) :: direction(3,2,1,2), requested(3,2,2)
      real(dblprec) :: aligned(3,8,2), cone_a(3,8,2), cone_b(3,8,2)
      logical :: defined(2,1,2)
      integer(int64) :: seeds_a(2,2), seeds_b(2,2)

      cell = 0.0_dblprec
      cell(1,1) = 2.0d-10
      cell(2,2) = 2.0d-10
      cell(3,3) = 2.0d-10
      call build_block_topology(topology,REGULAR_REPLICATED_CELL,2,(/4,1,1/),8, &
         (/4,1,1/),cell,(/1,2/),topology_status,message)
      call check(topology_status == BLOCK_TOPOLOGY_OK,'multi-channel topology builds')
      moment = (/1.0_dblprec,1.4_dblprec,0.8_dblprec,1.2_dblprec, &
         1.1_dblprec,0.9_dblprec,1.3_dblprec,0.7_dblprec/)
      spin = 0.0_dblprec
      do atom = 1, 8
         channel = topology%atom_to_dynamic_channel(atom)
         if (channel == 1) then
            spin(:,atom,1) = (/1.0_dblprec,0.0_dblprec,0.0_dblprec/)
         else
            spin(:,atom,1) = (/-1.0_dblprec,0.0_dblprec,0.0_dblprec/)
         end if
         spin(:,atom,2) = spin(:,atom,1)
      end do
      call restrict_channel_moments(topology,moment,spin,resultant,moment_sum, &
         direction,defined,status,message)
      call check(status == ADAPTIVE_HYBRID_OK .and. all(defined), &
         'opposite compensated channels remain independently normalized')
      call check(norm3(resultant(:,1,1,1)+resultant(:,2,1,1)) < 0.7_dblprec .and. &
         direction(1,1,1,1) > 0.99_dblprec .and. direction(1,2,1,1) < -0.99_dblprec, &
         'multi-channel restriction does not normalize the compensated net vector')

      ! Force a true zero resultant inside channel one; it must be reported as
      ! undefined rather than divided by zero.
      spin(:,topology%block_atoms(1),1) = (/1.0_dblprec,0.0_dblprec,0.0_dblprec/)
      spin(:,topology%block_atoms(3),1) = (/-1.0_dblprec,0.0_dblprec,0.0_dblprec/)
      spin(:,topology%block_atoms(5),1) = (/1.0_dblprec,0.0_dblprec,0.0_dblprec/)
      spin(:,topology%block_atoms(7),1) = (/-1.0_dblprec,0.0_dblprec,0.0_dblprec/)
      moment = 1.0_dblprec
      call restrict_channel_moments(topology,moment,spin,resultant,moment_sum, &
         direction,defined,status,message)
      call check(.not. defined(1,1,1) .and. &
         maxval(abs(direction(:,1,1,1))) <= tiny(1.0_dblprec), &
         'a zero channel resultant is left undefined without normalization')

      requested = 0.0_dblprec
      do ensemble = 1, 2
         requested(:,1,ensemble) = 4.0_dblprec*unit((/0.7_dblprec,0.2_dblprec,0.4_dblprec/))
         requested(:,2,ensemble) = 4.0_dblprec*unit((/-0.3_dblprec,0.8_dblprec,0.2_dblprec/))
      end do
      aligned = spin
      call reconstruct_block_aligned(topology,1,moment,requested,aligned,1.0d-12, &
         status,message)
      call restrict_channel_moments(topology,moment,aligned,resultant,moment_sum, &
         direction,defined,status,message)
      call check(status == ADAPTIVE_HYBRID_OK .and. &
         maxval(abs(resultant(:,:,1,:)-requested)) < 2.0d-14, &
         'aligned reconstruction exactly restores the full channel resultants')

      requested = 0.96_dblprec*requested
      cone_a = aligned
      cone_b = aligned
      call reconstruct_block_constrained_cone(topology,1,moment,requested,cone_a, &
         0.8_dblprec,1.0d-12,99173_int64,7_c_int,seeds_a,status,message)
      call check(status == ADAPTIVE_HYBRID_OK,'constrained cone reconstruction succeeds: '//trim(message))
      call reconstruct_block_constrained_cone(topology,1,moment,requested,cone_b, &
         0.8_dblprec,1.0d-12,99173_int64,7_c_int,seeds_b,status,message)
      call restrict_channel_moments(topology,moment,cone_a,resultant,moment_sum, &
         direction,defined,status,message)
      call check(maxval(abs(resultant(:,:,1,:)-requested)) < 5.0d-12 .and. &
         maxval(abs(cone_a-cone_b)) <= tiny(1.0_dblprec) .and. all(seeds_a == seeds_b), &
         'cone replay is deterministic and matches the requested resultant')
      call check(seeds_a(1,1) /= seeds_a(2,1) .and. seeds_a(1,1) /= seeds_a(1,2) .and. &
         deterministic_reconstruction_seed(99173_int64,1,1,1,7) /= &
         deterministic_reconstruction_seed(99173_int64,2,1,1,7) .and. &
         deterministic_reconstruction_seed(99173_int64,1,1,1,7) /= &
         deterministic_reconstruction_seed(99173_int64,1,1,1,8), &
         'seeds resolve block/channel/ensemble/epoch identities')
   end subroutine test_channel_restriction_and_reconstruction

   subroutine test_transactional_chain_transitions()
      type(block_topology_type) :: topology
      type(static_hybrid_operator_type) :: hybrid, hybrid_before
      type(adaptive_hybrid_runtime_type) :: runtime
      type(selector_requests_type) :: requests
      type(selector_evaluation_type) :: evaluation
      type(selector_configuration_type) :: selector_configuration
      type(adaptive_reconstruction_configuration_type) :: reconstruction
      integer, parameter :: n = 32
      integer :: atom, block, status, cycle
      integer(c_int) :: state_before(n)
      character(len=512) :: message
      integer :: bonds(2,n)
      real(dblprec) :: displacement(3,n), matrix(3,3,n), moment(n)
      real(dblprec) :: atom_direction(3,n,1), atom_before(3,n,1)
      real(dblprec) :: coarse_direction(3,1,n,1), x, profile
      logical :: hard(n)

      call make_chain_fixture(n,topology,hybrid,bonds,displacement,matrix,status,message)
      call check(status == ADAPTIVE_HYBRID_OK,'adaptive chain fixture builds: '//trim(message))
      evaluator_matrix = matrix
      moment = 1.7_dblprec
      do atom = 1, n
         x = real(atom-1,dblprec)
         profile = tanh((x-10.5_dblprec)/3.5_dblprec) - &
            tanh((x-25.5_dblprec)/3.5_dblprec) - 1.0_dblprec
         atom_direction(:,atom,1) = unit((/sqrt(max(0.0_dblprec,1.0_dblprec-profile*profile)), &
            0.0_dblprec,profile/))
         coarse_direction(:,1,atom,1) = atom_direction(:,atom,1)
      end do
      call setup_adaptive_hybrid_runtime(runtime,topology,hybrid,moment,atom_direction, &
         coarse_direction,status,message)
      call check(status == ADAPTIVE_HYBRID_OK,'adaptive chain runtime initializes')
      call synthetic_selector(n,evaluation,requests)
      hard = .false.
      block = 12
      requests%refine(block) = .true.
      reconstruction%energy_jump_limit_j = huge(1.0_dblprec)

      state_before = runtime%selector%resolution_state
      call apply_adaptive_transitions(runtime,topology,requests,hard,evaluation, &
         selector_configuration,reconstruction,0_c_int,ADAPTIVE_STAGE_PREDICTOR, &
         moment,atom_direction,coarse_direction,hybrid_energy,status,message)
      call check(status == ADAPTIVE_HYBRID_INVALID_STAGE .and. &
         all(runtime%selector%resolution_state == state_before), &
         'selector cannot transition inside the predictor stage')

      ! A deterministic synthetic energy jump forces rejection.  Every live
      ! state component and ownership list must retain its prior value.
      evaluator_mask_penalty = 1.0d-18
      reconstruction%energy_jump_limit_j = 0.5d-18
      atom_before = atom_direction
      hybrid_before = runtime%hybrid
      call apply_adaptive_transitions(runtime,topology,requests,hard,evaluation, &
         selector_configuration,reconstruction,0_c_int,ADAPTIVE_STAGE_COMPLETE_STEP, &
         moment,atom_direction,coarse_direction,hybrid_energy,status,message)
      call check(status == ADAPTIVE_HYBRID_OK .and. &
         runtime%selector%resolution_state(block) == RESOLUTION_COARSE .and. &
         runtime%selector%transition_epoch(block) == 0 .and. &
         maxval(abs(atom_direction-atom_before)) <= tiny(1.0_dblprec) .and. &
         all(runtime%hybrid%atomistic_bond_owner .eqv. hybrid_before%atomistic_bond_owner), &
         'energy-jump rejection rolls back state, epoch, directions, and ownership')
      call check(allocated(runtime%transition_log%event) .and. &
         .not. runtime%transition_log%event(1)%accepted .and. &
         trim(runtime%transition_log%event(1)%outcome) == 'energy-jump-rejected' .and. &
         abs(runtime%transition_log%event(1)%energy_jump_j) > &
         reconstruction%energy_jump_limit_j, &
         'rejected transition logs its measured energy jump')

      evaluator_mask_penalty = 0.0_dblprec
      reconstruction%energy_jump_limit_j = huge(1.0_dblprec)
      do cycle = 1, 3
         requests%refine = .false.
         requests%coarsen = .false.
         requests%refine(block) = .true.
         call apply_adaptive_transitions(runtime,topology,requests,hard,evaluation, &
            selector_configuration,reconstruction,int(2*cycle-1,c_int), &
            ADAPTIVE_STAGE_COMPLETE_STEP,moment,atom_direction,coarse_direction, &
            hybrid_energy,status,message)
         call check(status == ADAPTIVE_HYBRID_OK .and. &
            runtime%selector%resolution_state(block) == RESOLUTION_ATOMISTIC, &
            'domain-wall block refines at a complete step')
         call check_runtime_invariants(runtime,topology,bonds)

         requests%refine = .false.
         requests%coarsen = .false.
         requests%coarsen(block) = .true.
         call apply_adaptive_transitions(runtime,topology,requests,hard,evaluation, &
            selector_configuration,reconstruction,int(2*cycle,c_int), &
            ADAPTIVE_STAGE_COMPLETE_STEP,moment,atom_direction,coarse_direction, &
            hybrid_energy,status,message)
         call check(status == ADAPTIVE_HYBRID_OK .and. &
            runtime%selector%resolution_state(block) == RESOLUTION_COARSE, &
            'domain-wall block coarsens after restriction')
         call check_runtime_invariants(runtime,topology,bonds)
      end do
      call check(count([(runtime%transition_log%event(cycle)%accepted, &
         cycle=1,size(runtime%transition_log%event))]) == 6, &
         'repeated domain-wall transitions are accepted and logged')
   end subroutine test_transactional_chain_transitions

   !> RCG-06A (F-17) allocation lifecycle negative control: the persistent
   !> per-candidate-transition workspace must be allocated exactly once by
   !> setup_adaptive_hybrid_runtime, reused (not reallocated) on every later
   !> preflight call with matching geometry, rejected with a diagnostic
   !> rather than silently corrupting memory if the geometry has drifted,
   !> and fully deallocated on cleanup.
   subroutine test_transition_workspace_lifecycle()
      type(block_topology_type) :: topology
      type(static_hybrid_operator_type) :: hybrid
      type(adaptive_hybrid_runtime_type) :: runtime
      integer, parameter :: n = 16
      integer :: status
      character(len=512) :: message
      integer :: bonds(2,n)
      real(dblprec) :: displacement(3,n), matrix(3,3,n), moment(n)
      real(dblprec) :: atom_direction(3,n,1), coarse_direction(3,1,n,1)
      integer :: atom

      call make_chain_fixture(n,topology,hybrid,bonds,displacement,matrix,status,message)
      call check(status == ADAPTIVE_HYBRID_OK,'lifecycle fixture builds: '//trim(message))
      moment = 1.7_dblprec
      do atom = 1, n
         atom_direction(:,atom,1) = (/0.0_dblprec,0.0_dblprec,1.0_dblprec/)
         coarse_direction(:,1,atom,1) = atom_direction(:,atom,1)
      end do
      call setup_adaptive_hybrid_runtime(runtime,topology,hybrid,moment,atom_direction, &
         coarse_direction,status,message)
      call check(status == ADAPTIVE_HYBRID_OK,'lifecycle runtime initializes')
      call check(runtime%transition_workspace_ready .and. &
         allocated(runtime%candidate_atom) .and. &
         all(shape(runtime%candidate_atom) == (/3,n,1/)) .and. &
         allocated(runtime%candidate_seeds), &
         'transition workspace is allocated once at setup with matching shape')

      ! Reuse: a second preflight call at the same geometry must succeed
      ! without altering the already-correct allocation.
      call ensure_transition_workspace(runtime,topology,1,status,message)
      call check(status == ADAPTIVE_HYBRID_OK .and. &
         all(shape(runtime%candidate_atom) == (/3,n,1/)), &
         'transition workspace preflight reuses an already-correct allocation')

      ! Drift: simulate a geometry change after allocation (this should
      ! never happen through the public API, but the preflight must still
      ! reject it explicitly rather than silently reading/writing past the
      ! stale allocation on the next apply_adaptive_transitions call).
      deallocate(runtime%candidate_atom)
      allocate(runtime%candidate_atom(3,n+1,1))
      call ensure_transition_workspace(runtime,topology,1,status,message)
      call check(status == ADAPTIVE_HYBRID_ALLOCATION_FAILED .and. len_trim(message) > 0, &
         'transition workspace preflight rejects a post-allocation geometry drift')

      ! Cleanup: reassigning the runtime to its default structure
      ! constructor (as cleanup_adaptive_cg_production does for the whole
      ! adaptive_cg_state) must deallocate every workspace field.
      runtime = adaptive_hybrid_runtime_type()
      call check(.not. runtime%transition_workspace_ready .and. &
         .not. allocated(runtime%candidate_atom) .and. &
         .not. allocated(runtime%candidate_seeds), &
         'cleanup deallocates the transition workspace')

      ! Resize: a fresh setup at a different topology size must allocate at
      ! the new size, not reuse the (already-deallocated) prior geometry.
      call setup_adaptive_hybrid_runtime(runtime,topology,hybrid,moment,atom_direction, &
         coarse_direction,status,message)
      call check(status == ADAPTIVE_HYBRID_OK .and. runtime%transition_workspace_ready .and. &
         all(shape(runtime%candidate_atom) == (/3,n,1/)), &
         'transition workspace reallocates correctly after cleanup and re-setup')
   end subroutine test_transition_workspace_lifecycle

   subroutine test_adaptive_skyrmion_transition()
      type(block_topology_type) :: topology
      type(static_hybrid_operator_type) :: hybrid
      type(adaptive_hybrid_runtime_type) :: runtime
      type(selector_requests_type) :: requests
      type(selector_evaluation_type) :: evaluation
      type(selector_configuration_type) :: selector_configuration
      type(adaptive_reconstruction_configuration_type) :: reconstruction
      integer, parameter :: nx = 8, ny = 8, n = nx*ny
      integer :: atom, x, y, block, status
      character(len=512) :: message
      integer :: bonds(2,2*n)
      real(dblprec) :: displacement(3,2*n), matrix(3,3,2*n), moment(n)
      real(dblprec) :: atom_direction(3,n,1), coarse_direction(3,1,n,1)
      real(dblprec) :: dx, dy, radius, theta, phi
      logical :: hard(n)

      call make_square_fixture(nx,ny,topology,hybrid,bonds,displacement,matrix,status,message)
      call check(status == ADAPTIVE_HYBRID_OK,'adaptive square fixture builds')
      evaluator_matrix = matrix
      moment = 1.7_dblprec
      do y = 0, ny-1
         do x = 0, nx-1
            atom = 1+x+nx*y
            dx = real(x,dblprec)-3.5_dblprec
            dy = real(y,dblprec)-3.5_dblprec
            radius = sqrt(dx*dx+dy*dy)
            theta = pi*exp(-(radius/2.2_dblprec)**2)
            phi = atan2(dy,dx)+0.5_dblprec*pi
            atom_direction(:,atom,1) = (/sin(theta)*cos(phi),sin(theta)*sin(phi),cos(theta)/)
            coarse_direction(:,1,atom,1) = atom_direction(:,atom,1)
         end do
      end do
      call setup_adaptive_hybrid_runtime(runtime,topology,hybrid,moment,atom_direction, &
         coarse_direction,status,message)
      call synthetic_selector(n,evaluation,requests)
      hard = .false.
      block = 1+3+nx*3
      requests%refine(block) = .true.
      reconstruction%energy_jump_limit_j = huge(1.0_dblprec)
      call apply_adaptive_transitions(runtime,topology,requests,hard,evaluation, &
         selector_configuration,reconstruction,1_c_int,ADAPTIVE_STAGE_COMPLETE_STEP, &
         moment,atom_direction,coarse_direction,hybrid_energy,status,message)
      call check(status == ADAPTIVE_HYBRID_OK .and. &
         runtime%selector%resolution_state(block) == RESOLUTION_ATOMISTIC .and. &
         runtime%transition_log%event(1)%accepted .and. &
         ieee_is_finite_scalar(runtime%transition_log%event(1)%energy_jump_j), &
         'skyrmion adaptive transition is accepted with a finite measured jump')
      call check_runtime_invariants(runtime,topology,bonds)
   end subroutine test_adaptive_skyrmion_transition

   !> RCG-03: an already-coarse block whose dormant atomistic state has since
   !> become polarization-unsafe must be forced back to atomistic only at an
   !> accepted synchronization point, with a reason and ratio distinct from
   !> the pre-existing static-mask hard exclusion.  Production reconstruction
   !> (reconstruct_coarse_atoms) always rebuilds a coarse block's dormant
   !> atoms to match its coarse channel exactly, so a genuinely reduced
   !> resultant can only be observed at the operator layer exercised here
   !> (see docs/RCG-03_POLARIZATION_ANISOTROPY_EVIDENCE.md for why a
   !> production-executable fixture cannot reach this path).
   subroutine test_polarization_forces_refine_of_coarse_block()
      type(block_topology_type) :: topology
      type(static_hybrid_operator_type) :: hybrid
      type(adaptive_hybrid_runtime_type) :: runtime
      type(selector_requests_type) :: requests
      type(selector_evaluation_type) :: evaluation
      type(selector_configuration_type) :: selector_configuration
      type(adaptive_reconstruction_configuration_type) :: reconstruction
      integer, parameter :: n = 4
      integer :: status, hybrid_status
      integer(c_int) :: state_before
      character(len=512) :: message
      integer :: bonds(2,n)
      real(dblprec) :: displacement(3,n), matrix(3,3,n), moment(n)
      real(dblprec) :: atom_direction(3,n,1), coarse_direction(3,1,1,1)
      real(dblprec) :: resultant(3,1,1,1), moment_sum(1,1,1), direction(3,1,1,1)
      real(dblprec) :: ratio(1)
      logical :: defined(1,1,1), unsafe(1)

      call make_single_block_fixture(n,topology,hybrid,bonds,displacement,matrix,status,message)
      call check(status == ADAPTIVE_HYBRID_OK,'single-block polarization fixture builds: '//trim(message))
      evaluator_matrix = matrix
      evaluator_mask_penalty = 0.0_dblprec
      moment = 1.7_dblprec
      atom_direction = 0.0_dblprec
      atom_direction(3,:,1) = 1.0_dblprec
      coarse_direction(:,1,1,1) = (/0.0_dblprec,0.0_dblprec,1.0_dblprec/)
      call setup_adaptive_hybrid_runtime(runtime,topology,hybrid,moment,atom_direction, &
         coarse_direction,status,message)
      call check(status == ADAPTIVE_HYBRID_OK,'single-block polarization runtime initializes')
      call check(runtime%selector%resolution_state(1) == RESOLUTION_COARSE, &
         'the block starts already coarse, before any evolution is applied')

      ! Evolve the block's dormant atomistic state into genuine, honest low
      ! polarization: three of its four members stay aligned and one flips
      ! fully, giving an exact resultant/moment-sum ratio of 0.5.
      atom_direction(3,4,1) = -1.0_dblprec
      call restrict_channel_moments(topology,moment,atom_direction,resultant,moment_sum, &
         direction,defined,hybrid_status,message)
      call check(hybrid_status == ADAPTIVE_HYBRID_OK .and. defined(1,1,1), &
         'the evolved block state has a well-defined but reduced resultant')
      call evaluate_polarization_gate(topology,resultant,moment_sum,defined,0.9_dblprec, &
         unsafe,hybrid_status,message,ratio)
      call check(hybrid_status == ADAPTIVE_HYBRID_OK .and. unsafe(1) .and. &
         abs(ratio(1)-0.5_dblprec) < 1.0d-12, &
         'the polarization gate measures the exact 0.5 ratio and flags the block unsafe')

      allocate(evaluation%score(1,1),evaluation%refine(1,1),evaluation%coarsen(1,1), &
         requests%refine(1),requests%coarsen(1))
      evaluation%score = 0.0_dblprec
      evaluation%refine = .false.
      evaluation%coarsen = .false.
      requests%refine = .false.
      requests%coarsen = .false.
      reconstruction%energy_jump_limit_j = huge(1.0_dblprec)

      state_before = runtime%selector%resolution_state(1)
      call apply_adaptive_transitions(runtime,topology,requests,unsafe,evaluation, &
         selector_configuration,reconstruction,0_c_int,ADAPTIVE_STAGE_PREDICTOR, &
         moment,atom_direction,coarse_direction,hybrid_energy,status,message, &
         unsafe,ratio)
      call check(status == ADAPTIVE_HYBRID_INVALID_STAGE .and. &
         runtime%selector%resolution_state(1) == state_before, &
         'an already-coarse polarization-unsafe block cannot be forced to refine mid-integrator-stage')

      call apply_adaptive_transitions(runtime,topology,requests,unsafe,evaluation, &
         selector_configuration,reconstruction,0_c_int,ADAPTIVE_STAGE_COMPLETE_STEP, &
         moment,atom_direction,coarse_direction,hybrid_energy,status,message, &
         unsafe,ratio)
      call check(status == ADAPTIVE_HYBRID_OK .and. &
         runtime%selector%resolution_state(1) == RESOLUTION_ATOMISTIC, &
         'the already-coarse block is forced atomistic at the next accepted synchronization point')
      call check(allocated(runtime%transition_log%event) .and. &
         size(runtime%transition_log%event) == 1 .and. &
         runtime%transition_log%event(1)%accepted .and. &
         trim(runtime%transition_log%event(1)%selector_reason) == 'polarization-unsafe' .and. &
         abs(runtime%transition_log%event(1)%polarization_ratio-0.5_dblprec) < 1.0d-12, &
         'the forced transition is logged with its own reason and the measured polarization ratio')
   end subroutine test_polarization_forces_refine_of_coarse_block

   subroutine make_single_block_fixture(n,topology,hybrid,bonds,displacement,matrix,status,message)
      integer, intent(in) :: n
      type(block_topology_type), intent(out) :: topology
      type(static_hybrid_operator_type), intent(out) :: hybrid
      integer, intent(out) :: bonds(2,n)
      real(dblprec), intent(out) :: displacement(3,n), matrix(3,3,n)
      integer, intent(out) :: status
      character(len=*), intent(out) :: message
      type(coarse_tensor_operator_type) :: tensor
      type(smooth_projected_operator_type) :: projection
      real(dblprec) :: coordinate(3,n)
      logical :: mask(1)

      ! width=n puts every atom in the single spatial block, unlike
      ! make_chain_fixture's width=1 (one atom per block) -- needed so the
      ! polarization ratio of a block can vary at all.  mask is per-block
      ! (one block here), not per-atom.
      call make_operator_fixture(n,n,topology,tensor,projection,bonds,displacement, &
         matrix,coordinate,status,message)
      if (status /= ADAPTIVE_HYBRID_OK) return
      mask = .false.
      call setup_static_hybrid_operator(hybrid,topology,tensor,projection,mask,bonds, &
         displacement,0,status,message)
      if (status == STATIC_HYBRID_OK) status = ADAPTIVE_HYBRID_OK
   end subroutine make_single_block_fixture

   subroutine hybrid_energy(operator,atom_direction,coarse_direction,energy_j,status,diagnostic)
      type(static_hybrid_operator_type), intent(inout) :: operator
      real(dblprec), intent(in) :: atom_direction(:,:,:)
      real(dblprec), intent(in) :: coarse_direction(:,:,:,:)
      real(dblprec), intent(out) :: energy_j
      integer, intent(out) :: status
      character(len=*), intent(out) :: diagnostic
      type(static_hybrid_energy_type) :: energy
      real(dblprec) :: fine_field(3,operator%n_atoms)
      real(dblprec) :: coarse_field(3,1,operator%n_blocks)
      integer :: ensemble

      energy_j = evaluator_mask_penalty*real(count(operator%fine_seed),dblprec)
      status = 0
      diagnostic = ''
      do ensemble = 1, size(atom_direction,3)
         call evaluate_static_hybrid_operator(operator,atom_direction(:,:,ensemble), &
            coarse_direction(:,:,:,ensemble),evaluator_matrix,fine_field,coarse_field, &
            energy,status,diagnostic)
         if (status /= STATIC_HYBRID_OK) return
         energy_j = energy_j+energy%total_j
      end do
   end subroutine hybrid_energy

   subroutine make_chain_fixture(n,topology,hybrid,bonds,displacement,matrix,status,message)
      integer, intent(in) :: n
      type(block_topology_type), intent(out) :: topology
      type(static_hybrid_operator_type), intent(out) :: hybrid
      integer, intent(out) :: bonds(2,n)
      real(dblprec), intent(out) :: displacement(3,n), matrix(3,3,n)
      integer, intent(out) :: status
      character(len=*), intent(out) :: message
      type(coarse_tensor_operator_type) :: tensor
      type(smooth_projected_operator_type) :: projection
      real(dblprec) :: coordinate(3,n)
      logical :: mask(n)

      call make_operator_fixture(n,1,topology,tensor,projection,bonds,displacement, &
         matrix,coordinate,status,message)
      if (status /= ADAPTIVE_HYBRID_OK) return
      mask = .false.
      call setup_static_hybrid_operator(hybrid,topology,tensor,projection,mask,bonds, &
         displacement,0,status,message)
      if (status == STATIC_HYBRID_OK) status = ADAPTIVE_HYBRID_OK
   end subroutine make_chain_fixture

   subroutine make_square_fixture(nx,ny,topology,hybrid,bonds,displacement,matrix,status,message)
      integer, intent(in) :: nx, ny
      type(block_topology_type), intent(out) :: topology
      type(static_hybrid_operator_type), intent(out) :: hybrid
      integer, intent(out) :: bonds(:,:)
      real(dblprec), intent(out) :: displacement(:,:), matrix(:,:,:)
      integer, intent(out) :: status
      character(len=*), intent(out) :: message
      type(coarse_tensor_operator_type) :: tensor
      type(smooth_projected_operator_type) :: projection
      type(coarse_material_type) :: material
      type(coarse_operator_options_type) :: options
      integer :: atom, x, y, right, up, bond
      real(dblprec), parameter :: a = 2.0d-10, pair_j = 3.0d-21
      real(dblprec) :: cell(3,3), exchange(3,3), dmi(3,3)
      real(dblprec) :: coordinate(3,nx*ny), moment(nx*ny)
      logical :: mask(nx*ny)

      cell = 0.0_dblprec
      cell(1,1) = a
      cell(2,2) = a
      cell(3,3) = a
      call build_block_topology(topology,REGULAR_REPLICATED_CELL,1,(/nx,ny,1/), &
         nx*ny,(/1,1,1/),cell,(/1/),status,message)
      exchange = 0.0_dblprec
      exchange(1,1) = pair_j/(2.0_dblprec*a)
      exchange(2,2) = exchange(1,1)
      dmi = 0.0_dblprec
      call make_material(cell,exchange,dmi,material)
      call setup_coarse_tensor_operator(tensor,topology,material,options,status,message)
      coordinate = 0.0_dblprec
      moment = 1.7_dblprec
      bond = 0
      do y = 0, ny-1
         do x = 0, nx-1
            atom = 1+x+nx*y
            coordinate(:,atom) = (/real(x,dblprec),real(y,dblprec),0.0_dblprec/)
            right = 1+modulo(x+1,nx)+nx*y
            up = 1+x+nx*modulo(y+1,ny)
            bond = bond+1
            bonds(:,bond) = (/atom,right/)
            displacement(:,bond) = (/a,0.0_dblprec,0.0_dblprec/)
            matrix(:,:,bond) = 0.0_dblprec
            matrix(1,1,bond) = pair_j
            matrix(2,2,bond) = pair_j
            matrix(3,3,bond) = pair_j
            bond = bond+1
            bonds(:,bond) = (/atom,up/)
            displacement(:,bond) = (/0.0_dblprec,a,0.0_dblprec/)
            matrix(:,:,bond) = matrix(:,:,bond-1)
         end do
      end do
      call setup_smooth_projected_operator(projection,topology,coordinate,moment,status,message)
      mask = .false.
      call setup_static_hybrid_operator(hybrid,topology,tensor,projection,mask,bonds, &
         displacement,0,status,message)
      if (status == STATIC_HYBRID_OK) status = ADAPTIVE_HYBRID_OK
   end subroutine make_square_fixture

   subroutine make_operator_fixture(n,width,topology,tensor,projection,bonds, &
         displacement,matrix,coordinate,status,message)
      integer, intent(in) :: n, width
      type(block_topology_type), intent(out) :: topology
      type(coarse_tensor_operator_type), intent(out) :: tensor
      type(smooth_projected_operator_type), intent(out) :: projection
      integer, intent(out) :: bonds(2,n)
      real(dblprec), intent(out) :: displacement(3,n), matrix(3,3,n), coordinate(3,n)
      integer, intent(out) :: status
      character(len=*), intent(out) :: message
      type(coarse_material_type) :: material
      type(coarse_operator_options_type) :: options
      integer :: atom
      real(dblprec), parameter :: a = 2.0d-10, pair_j = 3.0d-21
      real(dblprec) :: cell(3,3), exchange(3,3), dmi(3,3), moment(n)

      cell = 0.0_dblprec
      cell(1,1) = a
      cell(2,2) = a
      cell(3,3) = a
      call build_block_topology(topology,REGULAR_REPLICATED_CELL,1,(/n,1,1/),n, &
         (/width,1,1/),cell,(/1/),status,message)
      exchange = 0.0_dblprec
      exchange(1,1) = pair_j/(2.0_dblprec*a)
      dmi = 0.0_dblprec
      call make_material(cell,exchange,dmi,material)
      call setup_coarse_tensor_operator(tensor,topology,material,options,status,message)
      coordinate = 0.0_dblprec
      moment = 1.7_dblprec
      do atom = 1, n
         coordinate(1,atom) = (real(atom,dblprec)-0.5_dblprec)/real(width,dblprec)-0.5_dblprec
         bonds(:,atom) = (/atom,modulo(atom,n)+1/)
         displacement(:,atom) = (/a,0.0_dblprec,0.0_dblprec/)
         matrix(:,:,atom) = 0.0_dblprec
         matrix(1,1,atom) = pair_j
         matrix(2,2,atom) = pair_j
         matrix(3,3,atom) = pair_j
      end do
      call setup_smooth_projected_operator(projection,topology,coordinate,moment,status,message)
      if (status == SMOOTH_PROJECTED_OK) status = ADAPTIVE_HYBRID_OK
   end subroutine make_operator_fixture

   subroutine make_material(cell,exchange,dmi,material)
      real(dblprec), intent(in) :: cell(3,3), exchange(3,3), dmi(3,3)
      type(coarse_material_type), intent(out) :: material

      material%ready = .true.
      material%n_basis = 1
      material%n_channels = 1
      material%cell_volume_m3 = cell(1,1)*cell(2,2)*cell(3,3)
      allocate(material%basis_to_channel(1),material%basis_moment_mub(1), &
         material%channel_moment_mub(1),material%channel_gamma(1), &
         material%channel_damping(1),material%exchange_stiffness(3,3,1,1), &
         material%spiralization(3,3,1,1))
      material%basis_to_channel = 1
      material%basis_moment_mub = 1.7_dblprec
      material%channel_moment_mub = 1.7_dblprec
      material%channel_gamma = 1.76d11
      material%channel_damping = 0.01_dblprec
      material%exchange_stiffness(:,:,1,1) = exchange
      material%spiralization(:,:,1,1) = dmi
      material%diagnostics%channel_dynamics_parameters_validated = .true.
      material%diagnostics%small_q_energy_validated = .true.
   end subroutine make_material

   subroutine synthetic_selector(n,evaluation,requests)
      integer, intent(in) :: n
      type(selector_evaluation_type), intent(out) :: evaluation
      type(selector_requests_type), intent(out) :: requests
      allocate(evaluation%score(1,n),evaluation%refine(1,n),evaluation%coarsen(1,n), &
         requests%refine(n),requests%coarsen(n))
      evaluation%score = 0.0_dblprec
      evaluation%refine = .false.
      evaluation%coarsen = .true.
      requests%refine = .false.
      requests%coarsen = .false.
   end subroutine synthetic_selector

   subroutine check_runtime_invariants(runtime,topology,bonds)
      type(adaptive_hybrid_runtime_type), intent(in) :: runtime
      type(block_topology_type), intent(in) :: topology
      integer, intent(in) :: bonds(:,:)
      integer :: bond
      logical :: expected

      call check(size(runtime%active_atom_list) == count(runtime%hybrid%atomistic_atom) .and. &
         size(runtime%active_coarse_list) == count(runtime%hybrid%coarse_block) .and. &
         size(runtime%interface_atom_list) == &
         sum(topology%block_atom_count,mask=runtime%hybrid%block_state == STATIC_HYBRID_BUFFER), &
         'accepted transition rebuilds compact active/interface lists')
      call check(count(runtime%hybrid%atomistic_block .and. runtime%hybrid%coarse_block) == 0 .and. &
         count(runtime%hybrid%atomistic_block .or. runtime%hybrid%coarse_block) == &
         topology%n_spatial_blocks,'accepted transition preserves block-state partition')
      do bond = 1, size(bonds,2)
         expected = runtime%hybrid%active_bond(bond) .and. &
            (runtime%hybrid%atomistic_atom(bonds(1,bond)) .or. &
             runtime%hybrid%atomistic_atom(bonds(2,bond)))
         call check(runtime%hybrid%atomistic_bond_owner(bond) .eqv. expected, &
            'accepted transition rebuilds interface bond ownership')
      end do
   end subroutine check_runtime_invariants

   pure function unit(vector) result(normalized)
      real(dblprec), intent(in) :: vector(3)
      real(dblprec) :: normalized(3)
      normalized = vector/norm3(vector)
   end function unit

   pure real(dblprec) function norm3(vector)
      real(dblprec), intent(in) :: vector(3)
      norm3 = sqrt(sum(vector*vector))
   end function norm3

   logical function ieee_is_finite_scalar(value)
      use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
      real(dblprec), intent(in) :: value
      ieee_is_finite_scalar = ieee_is_finite(value)
   end function ieee_is_finite_scalar

   subroutine check(condition,message)
      logical, intent(in) :: condition
      character(len=*), intent(in) :: message
      if (.not. condition) then
         failures = failures+1
         write(*,'(a)') 'FAIL: '//trim(message)
      end if
   end subroutine check

end program test_adaptive_hybrid_solver
