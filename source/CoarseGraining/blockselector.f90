!-------------------------------------------------------------------------------
! MODULE: BlockSelector
!> @brief Pure selection criteria plus an explicit adaptive-resolution sync.
!>
!> This module deliberately knows nothing about Hamiltonians, integrators,
!> fields, restriction, or reconstruction.  Criterion evaluation is read-only;
!> advance_selector_state is the only routine that changes block state.
!-------------------------------------------------------------------------------
module BlockSelector

   use, intrinsic :: iso_c_binding, only : c_double, c_int

   implicit none
   private

   integer(c_int), parameter, public :: SELECTOR_OK = 0_c_int
   integer(c_int), parameter, public :: SELECTOR_INVALID_ARGUMENT = 1_c_int
   integer(c_int), parameter, public :: SELECTOR_ALLOCATION_FAILED = 2_c_int
   integer(c_int), parameter, public :: SELECTOR_UNKNOWN_CRITERION = 3_c_int

   integer(c_int), parameter, public :: RESOLUTION_ATOMISTIC = 1_c_int
   integer(c_int), parameter, public :: RESOLUTION_BUFFER = 2_c_int
   integer(c_int), parameter, public :: RESOLUTION_COARSE = 3_c_int

   integer(c_int), parameter, public :: CRITERION_MAX_NEIGHBOUR_MISALIGNMENT = 1_c_int
   integer(c_int), parameter, public :: CRITERION_STATIC_MASK = 2_c_int

   type, public :: selector_configuration_type
      real(c_double) :: refine_threshold = 0.25_c_double
      real(c_double) :: coarsen_threshold = 0.10_c_double
      integer(c_int) :: update_interval = 1_c_int
      integer(c_int) :: minimum_dwell_updates = 0_c_int
      integer(c_int) :: buffer_dilation_blocks = 0_c_int
   end type selector_configuration_type

   ! A criterion descriptor contains immutable setup only.  In particular, it
   ! contains no pointer to a runtime and cannot mutate resolution_state.
   type, public :: selector_criterion_type
      integer(c_int) :: kind = 0_c_int
      character(len=48) :: name = ''
      logical, allocatable :: static_atomistic_mask(:)
   end type selector_criterion_type

   type, public :: selector_registry_type
      type(selector_criterion_type), allocatable :: criterion(:)
   end type selector_registry_type

   type, public :: selector_evaluation_type
      ! score(criterion, block), refine(criterion, block), and
      ! coarsen(criterion, block) are criterion outputs, never state storage.
      real(c_double), allocatable :: score(:,:)
      logical, allocatable :: refine(:,:)
      logical, allocatable :: coarsen(:,:)
   end type selector_evaluation_type

   type, public :: selector_requests_type
      logical, allocatable :: refine(:)
      logical, allocatable :: coarsen(:)
   end type selector_requests_type

   type, public :: block_selector_runtime_type
      integer(c_int), allocatable :: resolution_state(:)
      integer(c_int), allocatable :: state_age(:)
      integer(c_int), allocatable :: transition_epoch(:)
   end type block_selector_runtime_type

   type, public :: selector_transition_event_type
      integer(c_int) :: synchronization_step = 0_c_int
      integer(c_int) :: block = 0_c_int
      integer(c_int) :: old_state = RESOLUTION_COARSE
      integer(c_int) :: new_state = RESOLUTION_COARSE
      character(len=48) :: reason = ''
      real(c_double), allocatable :: score(:)
      logical, allocatable :: refine_request(:)
      logical, allocatable :: coarsen_request(:)
   end type selector_transition_event_type

   type, public :: selector_transition_log_type
      type(selector_transition_event_type), allocatable :: event(:)
   end type selector_transition_log_type

   public :: register_max_neighbour_misalignment
   public :: register_static_atomistic_mask
   public :: evaluate_selector_registry
   public :: combine_selector_requests
   public :: initialize_selector_runtime
   public :: advance_selector_state
   public :: dilate_atomistic_blocks
   public :: clear_transition_log

contains

   subroutine register_max_neighbour_misalignment(registry, status, diagnostic)
      type(selector_registry_type), intent(inout) :: registry
      integer(c_int), intent(out) :: status
      character(len=*), intent(out) :: diagnostic
      type(selector_criterion_type) :: descriptor

      descriptor%kind = CRITERION_MAX_NEIGHBOUR_MISALIGNMENT
      descriptor%name = 'maximum-neighbour-misalignment'
      call append_criterion(registry, descriptor, status, diagnostic)
   end subroutine register_max_neighbour_misalignment

   subroutine register_static_atomistic_mask(registry, atomistic_mask, status, diagnostic)
      type(selector_registry_type), intent(inout) :: registry
      logical, intent(in) :: atomistic_mask(:)
      integer(c_int), intent(out) :: status
      character(len=*), intent(out) :: diagnostic
      type(selector_criterion_type) :: descriptor
      integer :: allocation_status

      status = SELECTOR_OK
      diagnostic = ''
      if (size(atomistic_mask) == 0) then
         status = SELECTOR_INVALID_ARGUMENT
         diagnostic = 'A static atomistic mask must contain at least one block'
         return
      end if
      descriptor%kind = CRITERION_STATIC_MASK
      descriptor%name = 'static-atomistic-mask'
      allocate(descriptor%static_atomistic_mask(size(atomistic_mask)), stat=allocation_status)
      if (allocation_status /= 0) then
         status = SELECTOR_ALLOCATION_FAILED
         diagnostic = 'Unable to allocate the static atomistic mask criterion'
         return
      end if
      descriptor%static_atomistic_mask = atomistic_mask
      call append_criterion(registry, descriptor, status, diagnostic)
   end subroutine register_static_atomistic_mask

   !> Evaluate all registered criteria.  runtime is intentionally intent(in):
   !> a criterion may inspect existing state in a future extension, but cannot
   !> change it.  neighbour_i/j are one-based atom ids and each edge is used
   !> symmetrically, so a cross-block edge contributes to both endpoint blocks.
   subroutine evaluate_selector_registry(registry, runtime, spins, atom_to_block, &
         neighbour_i, neighbour_j, evaluation, status, diagnostic)
      type(selector_registry_type), intent(in) :: registry
      type(block_selector_runtime_type), intent(in) :: runtime
      real(c_double), intent(in), optional :: spins(:,:)
      integer(c_int), intent(in), optional :: atom_to_block(:), neighbour_i(:), neighbour_j(:)
      type(selector_evaluation_type), intent(out) :: evaluation
      integer(c_int), intent(out) :: status
      character(len=*), intent(out) :: diagnostic

      integer :: nblocks, ncriteria, criterion, allocation_status

      status = SELECTOR_OK
      diagnostic = ''
      if (.not. allocated(runtime%resolution_state)) then
         call fail_evaluation('Selector runtime has not been initialized')
         return
      end if
      nblocks = size(runtime%resolution_state)
      ncriteria = 0
      if (allocated(registry%criterion)) ncriteria = size(registry%criterion)
      if (ncriteria == 0) then
         call fail_evaluation('Selector registry contains no criteria')
         return
      end if
      allocate(evaluation%score(ncriteria,nblocks), evaluation%refine(ncriteria,nblocks), &
         evaluation%coarsen(ncriteria,nblocks), stat=allocation_status)
      if (allocation_status /= 0) then
         status = SELECTOR_ALLOCATION_FAILED
         diagnostic = 'Unable to allocate selector criterion outputs'
         return
      end if
      evaluation%score = 0.0_c_double
      evaluation%refine = .false.
      evaluation%coarsen = .true.

      do criterion = 1, ncriteria
         select case (registry%criterion(criterion)%kind)
         case (CRITERION_MAX_NEIGHBOUR_MISALIGNMENT)
            if (.not. present(spins) .or. .not. present(atom_to_block) .or. &
                .not. present(neighbour_i) .or. .not. present(neighbour_j)) then
               call fail_evaluation('Misalignment criterion requires spins, block ids, and neighbour pairs')
               return
            end if
            call evaluate_misalignment(spins, atom_to_block, neighbour_i, neighbour_j, &
               evaluation%score(criterion,:), status, diagnostic)
            if (status /= SELECTOR_OK) return
            ! Threshold comparisons happen in the combiner/state policy; this
            ! criterion returns score data only.
         case (CRITERION_STATIC_MASK)
            if (.not. allocated(registry%criterion(criterion)%static_atomistic_mask) .or. &
                size(registry%criterion(criterion)%static_atomistic_mask) /= nblocks) then
               call fail_evaluation('Static atomistic mask size does not match selector runtime')
               return
            end if
            evaluation%refine(criterion,:) = registry%criterion(criterion)%static_atomistic_mask
            evaluation%coarsen(criterion,:) = .not. registry%criterion(criterion)%static_atomistic_mask
         case default
            status = SELECTOR_UNKNOWN_CRITERION
            diagnostic = 'Selector registry contains an unknown criterion kind'
            return
         end select
      end do

   contains

      subroutine fail_evaluation(message)
         character(len=*), intent(in) :: message

         status = SELECTOR_INVALID_ARGUMENT
         diagnostic = message
      end subroutine fail_evaluation
   end subroutine evaluate_selector_registry

   !> The configured policy is OR for refinement and AND for coarsening.
   !> Score criteria are converted to requests here, rather than inside the
   !> criteria, preserving the criterion/state mutation boundary.
   subroutine combine_selector_requests(evaluation, configuration, requests, status, diagnostic)
      type(selector_evaluation_type), intent(in) :: evaluation
      type(selector_configuration_type), intent(in) :: configuration
      type(selector_requests_type), intent(out) :: requests
      integer(c_int), intent(out) :: status
      character(len=*), intent(out) :: diagnostic

      integer :: ncriteria, nblocks, criterion, allocation_status

      status = SELECTOR_OK
      diagnostic = ''
      if (.not. allocated(evaluation%score) .or. .not. allocated(evaluation%refine) .or. &
          .not. allocated(evaluation%coarsen)) then
         call fail_combination('Selector evaluation is not allocated')
         return
      end if
      ncriteria = size(evaluation%score,1)
      nblocks = size(evaluation%score,2)
      if (size(evaluation%refine,1) /= ncriteria .or. size(evaluation%refine,2) /= nblocks .or. &
          size(evaluation%coarsen,1) /= ncriteria .or. size(evaluation%coarsen,2) /= nblocks) then
         call fail_combination('Selector evaluation arrays have incompatible shapes')
         return
      end if
      if (configuration%coarsen_threshold > configuration%refine_threshold .or. &
          configuration%update_interval <= 0 .or. configuration%minimum_dwell_updates < 0) then
         call fail_combination('Selector thresholds or update cadence are invalid')
         return
      end if
      allocate(requests%refine(nblocks), requests%coarsen(nblocks), stat=allocation_status)
      if (allocation_status /= 0) then
         status = SELECTOR_ALLOCATION_FAILED
         diagnostic = 'Unable to allocate combined selector requests'
         return
      end if
      requests%refine = .false.
      requests%coarsen = .true.
      do criterion = 1, ncriteria
         requests%refine = requests%refine .or. evaluation%refine(criterion,:) .or. &
            (evaluation%score(criterion,:) >= configuration%refine_threshold)
         requests%coarsen = requests%coarsen .and. evaluation%coarsen(criterion,:) .and. &
            (evaluation%score(criterion,:) <= configuration%coarsen_threshold)
      end do

   contains

      subroutine fail_combination(message)
         character(len=*), intent(in) :: message

         status = SELECTOR_INVALID_ARGUMENT
         diagnostic = message
      end subroutine fail_combination
   end subroutine combine_selector_requests

   subroutine initialize_selector_runtime(runtime, nblocks, initial_state, status, diagnostic)
      type(block_selector_runtime_type), intent(out) :: runtime
      integer, intent(in) :: nblocks
      integer(c_int), intent(in), optional :: initial_state
      integer(c_int), intent(out) :: status
      character(len=*), intent(out) :: diagnostic

      integer :: allocation_status
      integer(c_int) :: state

      status = SELECTOR_OK
      diagnostic = ''
      state = RESOLUTION_COARSE
      if (present(initial_state)) state = initial_state
      if (nblocks <= 0 .or. .not. valid_resolution_state(state)) then
         status = SELECTOR_INVALID_ARGUMENT
         diagnostic = 'Selector runtime requires positive blocks and a valid initial state'
         return
      end if
      allocate(runtime%resolution_state(nblocks), runtime%state_age(nblocks), &
         runtime%transition_epoch(nblocks), stat=allocation_status)
      if (allocation_status /= 0) then
         status = SELECTOR_ALLOCATION_FAILED
         diagnostic = 'Unable to allocate selector runtime state'
         return
      end if
      runtime%resolution_state = state
      ! A freshly initialized state may respond at the first synchronization.
      runtime%state_age = huge(0_c_int)
      runtime%transition_epoch = 0_c_int
   end subroutine initialize_selector_runtime

   !> The sole block-state mutation API.  Calls not on update_interval are
   !> intentional no-ops, so the caller may invoke it at every safe sync point.
   subroutine advance_selector_state(runtime, requests, hard_atomistic_mask, evaluation, &
         configuration, synchronization_step, transition_log, status, diagnostic)
      type(block_selector_runtime_type), intent(inout) :: runtime
      type(selector_requests_type), intent(in) :: requests
      logical, intent(in) :: hard_atomistic_mask(:)
      type(selector_evaluation_type), intent(in) :: evaluation
      type(selector_configuration_type), intent(in) :: configuration
      integer(c_int), intent(in) :: synchronization_step
      type(selector_transition_log_type), intent(inout) :: transition_log
      integer(c_int), intent(out) :: status
      character(len=*), intent(out) :: diagnostic

      integer :: block, nblocks
      integer(c_int) :: old_state, new_state
      character(len=48) :: reason

      status = SELECTOR_OK
      diagnostic = ''
      if (.not. runtime_is_valid(runtime) .or. .not. allocated(requests%refine) .or. &
          .not. allocated(requests%coarsen) .or. .not. allocated(evaluation%score)) then
         call fail_advance('Selector runtime, requests, or scores are not initialized')
         return
      end if
      nblocks = size(runtime%resolution_state)
      if (size(requests%refine) /= nblocks .or. size(requests%coarsen) /= nblocks .or. &
          size(hard_atomistic_mask) /= nblocks .or. size(evaluation%score,2) /= nblocks) then
         call fail_advance('Selector state, request, exclusion, and score sizes must match')
         return
      end if
      if (configuration%update_interval <= 0 .or. configuration%minimum_dwell_updates < 0) then
         call fail_advance('Selector update interval and dwell age must be nonnegative')
         return
      end if
      if (mod(synchronization_step, configuration%update_interval) /= 0_c_int) return

      do block = 1, nblocks
         old_state = runtime%resolution_state(block)
         new_state = old_state
         reason = ''
         if (runtime%state_age(block) < huge(0_c_int)) then
            runtime%state_age(block) = runtime%state_age(block) + 1_c_int
         end if
         if (hard_atomistic_mask(block)) then
            new_state = RESOLUTION_ATOMISTIC
            reason = 'hard-atomistic-exclusion'
         else if (old_state == RESOLUTION_COARSE .and. requests%refine(block) .and. &
                  runtime%state_age(block) >= configuration%minimum_dwell_updates) then
            new_state = RESOLUTION_ATOMISTIC
            reason = 'refine-request'
         else if (old_state == RESOLUTION_ATOMISTIC .and. requests%coarsen(block) .and. &
                  runtime%state_age(block) >= configuration%minimum_dwell_updates) then
            new_state = RESOLUTION_COARSE
            reason = 'coarsen-request'
         end if
         if (new_state /= old_state) then
            call append_transition(transition_log, synchronization_step, block, old_state, new_state, &
               reason, evaluation%score(:,block), evaluation%refine(:,block), evaluation%coarsen(:,block), &
               status, diagnostic)
            if (status /= SELECTOR_OK) return
            runtime%resolution_state(block) = new_state
            runtime%state_age(block) = 0_c_int
            runtime%transition_epoch(block) = runtime%transition_epoch(block) + 1_c_int
         end if
      end do

   contains

      subroutine fail_advance(message)
         character(len=*), intent(in) :: message

         status = SELECTOR_INVALID_ARGUMENT
         diagnostic = message
      end subroutine fail_advance
   end subroutine advance_selector_state

   !> Derive a buffer mask without changing runtime resolution state.  Dilation
   !> uses Chebyshev distance on the regular block grid and can wrap per axis.
   subroutine dilate_atomistic_blocks(runtime, block_grid, dilation_blocks, periodic, effective_state, &
         status, diagnostic)
      type(block_selector_runtime_type), intent(in) :: runtime
      integer, intent(in) :: block_grid(3)
      integer, intent(in) :: dilation_blocks
      logical, intent(in) :: periodic(3)
      integer(c_int), intent(out) :: effective_state(:)
      integer(c_int), intent(out) :: status
      character(len=*), intent(out) :: diagnostic

      integer :: nblocks, block, source, dx, dy, dz, source_coordinate(3), target_coordinate(3)

      status = SELECTOR_OK
      diagnostic = ''
      if (.not. runtime_is_valid(runtime)) then
         call fail_dilation('Selector runtime is not initialized')
         return
      end if
      nblocks = size(runtime%resolution_state)
      if (any(block_grid <= 0) .or. product(block_grid) /= nblocks .or. dilation_blocks < 0 .or. &
          size(effective_state) /= nblocks) then
         call fail_dilation('Invalid block grid, dilation width, or effective-state size')
         return
      end if
      effective_state = RESOLUTION_COARSE
      where (runtime%resolution_state == RESOLUTION_ATOMISTIC)
         effective_state = RESOLUTION_ATOMISTIC
      end where
      do source = 1, nblocks
         if (runtime%resolution_state(source) /= RESOLUTION_ATOMISTIC) cycle
         source_coordinate = block_coordinate(source, block_grid)
         do dz = -dilation_blocks, dilation_blocks
            do dy = -dilation_blocks, dilation_blocks
               do dx = -dilation_blocks, dilation_blocks
                  target_coordinate = source_coordinate + (/dx,dy,dz/)
                  if (.not. normalize_coordinate(target_coordinate, block_grid, periodic)) cycle
                  block = block_number(target_coordinate, block_grid)
                  if (effective_state(block) /= RESOLUTION_ATOMISTIC) effective_state(block) = RESOLUTION_BUFFER
               end do
            end do
         end do
      end do

   contains

      subroutine fail_dilation(message)
         character(len=*), intent(in) :: message

         status = SELECTOR_INVALID_ARGUMENT
         diagnostic = message
      end subroutine fail_dilation
   end subroutine dilate_atomistic_blocks

   subroutine clear_transition_log(transition_log)
      type(selector_transition_log_type), intent(inout) :: transition_log

      if (allocated(transition_log%event)) deallocate(transition_log%event)
   end subroutine clear_transition_log

   subroutine evaluate_misalignment(spins, atom_to_block, neighbour_i, neighbour_j, score, status, diagnostic)
      real(c_double), intent(in) :: spins(:,:)
      integer(c_int), intent(in) :: atom_to_block(:), neighbour_i(:), neighbour_j(:)
      real(c_double), intent(out) :: score(:)
      integer(c_int), intent(out) :: status
      character(len=*), intent(out) :: diagnostic

      integer :: edge, atom_i, atom_j, block_i, block_j
      real(c_double) :: q

      status = SELECTOR_OK
      diagnostic = ''
      if (size(spins,1) /= 3 .or. size(spins,2) /= size(atom_to_block) .or. &
          size(neighbour_i) /= size(neighbour_j) .or. any(atom_to_block < 1) .or. &
          any(atom_to_block > size(score))) then
         status = SELECTOR_INVALID_ARGUMENT
         diagnostic = 'Misalignment inputs have incompatible dimensions or block ids'
         return
      end if
      score = 0.0_c_double
      do edge = 1, size(neighbour_i)
         atom_i = neighbour_i(edge)
         atom_j = neighbour_j(edge)
         if (atom_i < 1 .or. atom_i > size(atom_to_block) .or. atom_j < 1 .or. atom_j > size(atom_to_block)) then
            status = SELECTOR_INVALID_ARGUMENT
            diagnostic = 'Misalignment neighbour pair contains an invalid atom id'
            return
         end if
         q = max(0.0_c_double, 1.0_c_double - dot_product(spins(:,atom_i),spins(:,atom_j)))
         block_i = atom_to_block(atom_i)
         block_j = atom_to_block(atom_j)
         score(block_i) = max(score(block_i),q)
         score(block_j) = max(score(block_j),q)
      end do
   end subroutine evaluate_misalignment

   subroutine append_criterion(registry, descriptor, status, diagnostic)
      type(selector_registry_type), intent(inout) :: registry
      type(selector_criterion_type), intent(in) :: descriptor
      integer(c_int), intent(out) :: status
      character(len=*), intent(out) :: diagnostic

      type(selector_criterion_type), allocatable :: expanded(:)
      integer :: old_count, allocation_status

      status = SELECTOR_OK
      diagnostic = ''
      old_count = 0
      if (allocated(registry%criterion)) old_count = size(registry%criterion)
      allocate(expanded(old_count+1), stat=allocation_status)
      if (allocation_status /= 0) then
         status = SELECTOR_ALLOCATION_FAILED
         diagnostic = 'Unable to grow selector criterion registry'
         return
      end if
      if (old_count > 0) expanded(:old_count) = registry%criterion
      expanded(old_count+1) = descriptor
      call move_alloc(expanded, registry%criterion)
   end subroutine append_criterion

   subroutine append_transition(transition_log, synchronization_step, block, old_state, new_state, reason, &
         score, refine_request, coarsen_request, status, diagnostic)
      type(selector_transition_log_type), intent(inout) :: transition_log
      integer(c_int), intent(in) :: synchronization_step, block, old_state, new_state
      character(len=*), intent(in) :: reason
      real(c_double), intent(in) :: score(:)
      logical, intent(in) :: refine_request(:), coarsen_request(:)
      integer(c_int), intent(out) :: status
      character(len=*), intent(out) :: diagnostic

      type(selector_transition_event_type), allocatable :: expanded(:)
      integer :: old_count, allocation_status

      status = SELECTOR_OK
      diagnostic = ''
      old_count = 0
      if (allocated(transition_log%event)) old_count = size(transition_log%event)
      allocate(expanded(old_count+1), stat=allocation_status)
      if (allocation_status /= 0) then
         status = SELECTOR_ALLOCATION_FAILED
         diagnostic = 'Unable to append a selector transition log event'
         return
      end if
      if (old_count > 0) expanded(:old_count) = transition_log%event
      expanded(old_count+1)%synchronization_step = synchronization_step
      expanded(old_count+1)%block = int(block,c_int)
      expanded(old_count+1)%old_state = old_state
      expanded(old_count+1)%new_state = new_state
      expanded(old_count+1)%reason = reason
      allocate(expanded(old_count+1)%score(size(score)), expanded(old_count+1)%refine_request(size(refine_request)), &
         expanded(old_count+1)%coarsen_request(size(coarsen_request)), stat=allocation_status)
      if (allocation_status /= 0) then
         status = SELECTOR_ALLOCATION_FAILED
         diagnostic = 'Unable to allocate selector transition score data'
         return
      end if
      expanded(old_count+1)%score = score
      expanded(old_count+1)%refine_request = refine_request
      expanded(old_count+1)%coarsen_request = coarsen_request
      call move_alloc(expanded, transition_log%event)
   end subroutine append_transition

   logical function runtime_is_valid(runtime)
      type(block_selector_runtime_type), intent(in) :: runtime

      runtime_is_valid = allocated(runtime%resolution_state) .and. allocated(runtime%state_age) .and. &
         allocated(runtime%transition_epoch)
      if (runtime_is_valid) runtime_is_valid = size(runtime%state_age) == size(runtime%resolution_state) .and. &
         size(runtime%transition_epoch) == size(runtime%resolution_state) .and. &
         all(runtime%resolution_state == RESOLUTION_ATOMISTIC .or. &
             runtime%resolution_state == RESOLUTION_COARSE)
   end function runtime_is_valid

   logical function valid_resolution_state(state)
      integer(c_int), intent(in) :: state

      valid_resolution_state = state == RESOLUTION_ATOMISTIC .or. state == RESOLUTION_COARSE
   end function valid_resolution_state

   pure function block_coordinate(block, block_grid) result(coordinate)
      integer, intent(in) :: block, block_grid(3)
      integer :: coordinate(3), zero_based

      zero_based = block - 1
      coordinate(1) = mod(zero_based,block_grid(1))
      zero_based = zero_based / block_grid(1)
      coordinate(2) = mod(zero_based,block_grid(2))
      coordinate(3) = zero_based / block_grid(2)
   end function block_coordinate

   pure integer function block_number(coordinate, block_grid)
      integer, intent(in) :: coordinate(3), block_grid(3)

      block_number = 1 + coordinate(1) + block_grid(1) * (coordinate(2) + block_grid(2) * coordinate(3))
   end function block_number

   logical function normalize_coordinate(coordinate, block_grid, periodic)
      integer, intent(inout) :: coordinate(3)
      integer, intent(in) :: block_grid(3)
      logical, intent(in) :: periodic(3)
      integer :: axis

      normalize_coordinate = .true.
      do axis = 1, 3
         if (coordinate(axis) < 0 .or. coordinate(axis) >= block_grid(axis)) then
            if (.not. periodic(axis)) then
               normalize_coordinate = .false.
               return
            end if
            coordinate(axis) = modulo(coordinate(axis),block_grid(axis))
         end if
      end do
   end function normalize_coordinate

end module BlockSelector
