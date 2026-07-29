program test_block_selector

   use, intrinsic :: iso_c_binding, only : c_double, c_int
   use BlockSelector

   implicit none

   integer :: failures

   failures = 0
   call test_misalignment_registry_is_read_only()
   call test_hysteresis_dwell_interval_and_log()
   call test_static_mask_and_combination()
   call test_hard_exclusion_and_buffer_dilation()
   call test_deterministic_replay()

   if (failures /= 0) then
      write(*,'(a,i0)') 'block selector tests failed: ', failures
      stop 1
   end if
   write(*,'(a)') 'block selector tests passed'

contains

   subroutine test_misalignment_registry_is_read_only()
      type(selector_registry_type) :: registry
      type(block_selector_runtime_type) :: runtime
      type(selector_evaluation_type) :: evaluation
      integer(c_int) :: status
      character(len=256) :: diagnostic
      real(c_double) :: spins(3,4)
      integer(c_int) :: atom_block(4), neighbour_i(3), neighbour_j(3), before(2)

      call initialize_selector_runtime(runtime, 2, status=status, diagnostic=diagnostic)
      call register_max_neighbour_misalignment(registry, status, diagnostic)
      spins = 0.0_c_double
      spins(3,:) = (/1.0_c_double, 1.0_c_double, -1.0_c_double, 1.0_c_double/)
      atom_block = (/1_c_int,1_c_int,2_c_int,2_c_int/)
      neighbour_i = (/1_c_int,2_c_int,3_c_int/)
      neighbour_j = (/2_c_int,3_c_int,4_c_int/)
      before = runtime%resolution_state
      call evaluate_selector_registry(registry, runtime, spins, atom_block, neighbour_i, neighbour_j, &
         evaluation, status, diagnostic)
      call check(status == SELECTOR_OK, 'misalignment registry evaluation succeeds')
      call check(all(abs(evaluation%score(1,:) - (/2.0_c_double,2.0_c_double/)) < 1.0e-14_c_double), &
         'maximum neighbour misalignment contributes to both endpoint blocks')
      call check(all(runtime%resolution_state == before), &
         'criteria do not mutate the block runtime')
      call check(all(.not. evaluation%refine) .and. all(evaluation%coarsen), &
         'score criterion returns data rather than state requests')
   end subroutine test_misalignment_registry_is_read_only

   subroutine test_hysteresis_dwell_interval_and_log()
      type(block_selector_runtime_type) :: runtime
      type(selector_evaluation_type) :: evaluation
      type(selector_requests_type) :: requests
      type(selector_transition_log_type) :: transition_log
      type(selector_configuration_type) :: configuration
      integer(c_int) :: status
      character(len=256) :: diagnostic
      logical :: hard(1)

      configuration%refine_threshold = 0.60_c_double
      configuration%coarsen_threshold = 0.20_c_double
      configuration%update_interval = 2_c_int
      configuration%minimum_dwell_updates = 2_c_int
      hard = .false.
      call initialize_selector_runtime(runtime, 1, status=status, diagnostic=diagnostic)
      call synthetic_evaluation((/0.80_c_double/), evaluation)
      call combine_selector_requests(evaluation, configuration, requests, status, diagnostic)
      call advance_selector_state(runtime, requests, hard, evaluation, configuration, 0_c_int, &
         transition_log, status, diagnostic)
      call check(runtime%resolution_state(1) == RESOLUTION_ATOMISTIC, 'high score refines at synchronization')
      call check(allocated(transition_log%event) .and. size(transition_log%event) == 1, &
         'refinement appends one deterministic transition')
      call check(trim(transition_log%event(1)%reason) == 'refine-request' .and. &
         abs(transition_log%event(1)%score(1)-0.80_c_double) < 1.0e-14_c_double, &
         'transition log includes reason and score data')

      call synthetic_evaluation((/0.10_c_double/), evaluation)
      call combine_selector_requests(evaluation, configuration, requests, status, diagnostic)
      call advance_selector_state(runtime, requests, hard, evaluation, configuration, 1_c_int, &
         transition_log, status, diagnostic)
      call check(runtime%resolution_state(1) == RESOLUTION_ATOMISTIC .and. size(transition_log%event) == 1, &
         'off-interval selector update is a no-op')
      call advance_selector_state(runtime, requests, hard, evaluation, configuration, 2_c_int, &
         transition_log, status, diagnostic)
      call check(runtime%resolution_state(1) == RESOLUTION_ATOMISTIC, &
         'dwell age suppresses the first coarsen request')
      call advance_selector_state(runtime, requests, hard, evaluation, configuration, 4_c_int, &
         transition_log, status, diagnostic)
      call check(runtime%resolution_state(1) == RESOLUTION_COARSE .and. size(transition_log%event) == 2, &
         'dwell age permits coarsening only after two selector updates')

      ! A synthetic high/low chatter sequence cannot flip at every update.
      call synthetic_evaluation((/0.80_c_double/), evaluation)
      call combine_selector_requests(evaluation, configuration, requests, status, diagnostic)
      call advance_selector_state(runtime, requests, hard, evaluation, configuration, 6_c_int, &
         transition_log, status, diagnostic)
      call synthetic_evaluation((/0.10_c_double/), evaluation)
      call combine_selector_requests(evaluation, configuration, requests, status, diagnostic)
      call advance_selector_state(runtime, requests, hard, evaluation, configuration, 8_c_int, &
         transition_log, status, diagnostic)
      call check(runtime%resolution_state(1) == RESOLUTION_COARSE .and. size(transition_log%event) == 2, &
         'synthetic chatter sequence does not chatter')
   end subroutine test_hysteresis_dwell_interval_and_log

   subroutine test_static_mask_and_combination()
      type(selector_registry_type) :: registry
      type(block_selector_runtime_type) :: runtime
      type(selector_evaluation_type) :: evaluation
      type(selector_requests_type) :: requests
      type(selector_transition_log_type) :: transition_log
      type(selector_configuration_type) :: configuration
      integer(c_int) :: status
      character(len=256) :: diagnostic
      logical :: static_mask(2), hard(2)

      static_mask = (/.true.,.false./)
      hard = .false.
      configuration%refine_threshold = 0.70_c_double
      configuration%coarsen_threshold = 0.20_c_double
      call initialize_selector_runtime(runtime, 2, status=status, diagnostic=diagnostic)
      call register_static_atomistic_mask(registry, static_mask, status, diagnostic)
      ! No spin/neighbour inputs: a static mask is itself a usable criterion.
      call evaluate_selector_registry(registry, runtime, evaluation=evaluation, status=status, diagnostic=diagnostic)
      call combine_selector_requests(evaluation, configuration, requests, status, diagnostic)
      call advance_selector_state(runtime, requests, hard, evaluation, configuration, 0_c_int, &
         transition_log, status, diagnostic)
      call check(status == SELECTOR_OK .and. all(runtime%resolution_state == &
         (/RESOLUTION_ATOMISTIC,RESOLUTION_COARSE/)), 'static mask works without score evaluation')
      call check(all(requests%refine .eqv. (/.true.,.false./)) .and. &
         all(requests%coarsen .eqv. (/.false.,.true./)), &
         'static criterion participates in OR-refine and AND-coarsen policy')

      ! An additional high score can request refinement, but static exclusion
      ! still vetoes coarsening only where its mask is set.
      call append_synthetic_row(evaluation, (/0.10_c_double,0.90_c_double/))
      call combine_selector_requests(evaluation, configuration, requests, status, diagnostic)
      call check(all(requests%refine .eqv. (/.true.,.true./)) .and. &
         all(requests%coarsen .eqv. (/.false.,.false./)), 'multiple criteria combine as specified')
   end subroutine test_static_mask_and_combination

   subroutine test_hard_exclusion_and_buffer_dilation()
      type(block_selector_runtime_type) :: runtime
      type(selector_evaluation_type) :: evaluation
      type(selector_requests_type) :: requests
      type(selector_transition_log_type) :: transition_log
      type(selector_configuration_type) :: configuration
      integer(c_int) :: status, effective(4)
      character(len=256) :: diagnostic
      logical :: hard(4), periodic(3)

      hard = (/.true.,.false.,.false.,.false./)
      periodic = (/.true.,.false.,.false./)
      call initialize_selector_runtime(runtime, 4, status=status, diagnostic=diagnostic)
      call synthetic_evaluation((/0.0_c_double,0.0_c_double,0.0_c_double,0.0_c_double/), evaluation)
      call combine_selector_requests(evaluation, configuration, requests, status, diagnostic)
      call advance_selector_state(runtime, requests, hard, evaluation, configuration, 0_c_int, &
         transition_log, status, diagnostic)
      call check(runtime%resolution_state(1) == RESOLUTION_ATOMISTIC, 'hard unsupported block remains atomistic')
      call dilate_atomistic_blocks(runtime, (/4,1,1/), 1, periodic, effective, status, diagnostic)
      call check(all(effective == (/RESOLUTION_ATOMISTIC,RESOLUTION_BUFFER,RESOLUTION_COARSE,RESOLUTION_BUFFER/)), &
         'periodic buffer dilation wraps at the boundary')
      periodic = .false.
      call dilate_atomistic_blocks(runtime, (/4,1,1/), 1, periodic, effective, status, diagnostic)
      call check(all(effective == (/RESOLUTION_ATOMISTIC,RESOLUTION_BUFFER,RESOLUTION_COARSE,RESOLUTION_COARSE/)), &
         'nonperiodic dilation clips at the boundary')
   end subroutine test_hard_exclusion_and_buffer_dilation

   subroutine test_deterministic_replay()
      type(block_selector_runtime_type) :: first_runtime, second_runtime
      type(selector_evaluation_type) :: evaluation
      type(selector_requests_type) :: requests
      type(selector_transition_log_type) :: first_log, second_log
      type(selector_configuration_type) :: configuration
      integer(c_int) :: status
      character(len=256) :: diagnostic
      logical :: hard(2)
      integer :: step, event
      real(c_double) :: score_sequence(4,2)

      configuration%refine_threshold = 0.70_c_double
      configuration%coarsen_threshold = 0.10_c_double
      configuration%minimum_dwell_updates = 1_c_int
      hard = .false.
      score_sequence = reshape((/0.8_c_double,0.0_c_double, 0.8_c_double,0.0_c_double, &
         0.0_c_double,0.0_c_double, 0.0_c_double,0.0_c_double/), shape(score_sequence))
      call initialize_selector_runtime(first_runtime, 2, status=status, diagnostic=diagnostic)
      call initialize_selector_runtime(second_runtime, 2, status=status, diagnostic=diagnostic)
      do step = 1, size(score_sequence,1)
         call synthetic_evaluation(score_sequence(step,:), evaluation)
         call combine_selector_requests(evaluation, configuration, requests, status, diagnostic)
         call advance_selector_state(first_runtime, requests, hard, evaluation, configuration, int(step-1,c_int), &
            first_log, status, diagnostic)
         call advance_selector_state(second_runtime, requests, hard, evaluation, configuration, int(step-1,c_int), &
            second_log, status, diagnostic)
      end do
      call check(all(first_runtime%resolution_state == second_runtime%resolution_state) .and. &
         size(first_log%event) == size(second_log%event), 'identical scores produce reproducible transition decisions')
      do event = 1, size(first_log%event)
         call check(first_log%event(event)%block == second_log%event(event)%block .and. &
            first_log%event(event)%synchronization_step == second_log%event(event)%synchronization_step .and. &
            trim(first_log%event(event)%reason) == trim(second_log%event(event)%reason) .and. &
            all(first_log%event(event)%score == second_log%event(event)%score), 'transition log order is deterministic')
      end do
   end subroutine test_deterministic_replay

   subroutine synthetic_evaluation(score, evaluation)
      real(c_double), intent(in) :: score(:)
      type(selector_evaluation_type), intent(out) :: evaluation

      allocate(evaluation%score(1,size(score)), evaluation%refine(1,size(score)), evaluation%coarsen(1,size(score)))
      evaluation%score(1,:) = score
      evaluation%refine = .false.
      evaluation%coarsen = .true.
   end subroutine synthetic_evaluation

   subroutine append_synthetic_row(evaluation, score)
      type(selector_evaluation_type), intent(inout) :: evaluation
      real(c_double), intent(in) :: score(:)
      real(c_double), allocatable :: expanded_score(:,:)
      logical, allocatable :: expanded_refine(:,:), expanded_coarsen(:,:)
      integer :: old_rows

      old_rows = size(evaluation%score,1)
      allocate(expanded_score(old_rows+1,size(score)), expanded_refine(old_rows+1,size(score)), &
         expanded_coarsen(old_rows+1,size(score)))
      expanded_score(:old_rows,:) = evaluation%score
      expanded_refine(:old_rows,:) = evaluation%refine
      expanded_coarsen(:old_rows,:) = evaluation%coarsen
      expanded_score(old_rows+1,:) = score
      expanded_refine(old_rows+1,:) = .false.
      expanded_coarsen(old_rows+1,:) = .true.
      call move_alloc(expanded_score,evaluation%score)
      call move_alloc(expanded_refine,evaluation%refine)
      call move_alloc(expanded_coarsen,evaluation%coarsen)
   end subroutine append_synthetic_row

   subroutine check(condition, label)
      logical, intent(in) :: condition
      character(len=*), intent(in) :: label

      if (.not. condition) then
         failures = failures + 1
         write(*,'(a)') 'FAIL: '//trim(label)
      end if
   end subroutine check

end program test_block_selector
