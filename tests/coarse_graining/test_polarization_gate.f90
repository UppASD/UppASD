program test_polarization_gate

   use, intrinsic :: iso_c_binding, only : c_double, c_int
   use Parameters, only : dblprec
   use BlockTopology, only : block_topology_type, build_block_topology, &
      REGULAR_REPLICATED_CELL, BLOCK_TOPOLOGY_OK
   use AdaptiveHybridSolver, only : evaluate_polarization_gate, &
      ADAPTIVE_HYBRID_OK, ADAPTIVE_HYBRID_INVALID_TOPOLOGY

   implicit none

   integer :: failures

   failures = 0
   call test_threshold_boundaries()
   call test_zero_resultant_never_normalizes()
   call test_invalid_arguments()

   if (failures /= 0) then
      write(*,'(a,i0)') 'polarization gate tests failed: ', failures
      stop 1
   end if
   write(*,'(a)') 'polarization gate tests passed'

contains

   !> POL-THRESHOLD: below, exactly at, above, and roundoff-scale-below the
   !> accepted threshold must each select the documented side of the "< "
   !> boundary (safe iff resultant/moment-sum ratio is not strictly below
   !> the threshold).
   subroutine test_threshold_boundaries()
      type(block_topology_type) :: topology
      integer(c_int) :: topology_status
      character(len=512) :: message
      real(c_double) :: cell(3,3)
      real(dblprec) :: resultant_mub(3,1,5,1), moment_sum_mub(1,5,1)
      logical :: direction_defined(1,5,1), unsafe_block(5)
      integer :: status
      character(len=256) :: diagnostic

      cell = identity_cell()
      call build_block_topology(topology, REGULAR_REPLICATED_CELL, 1, (/5,1,1/), 5, &
         (/1,1,1/), cell, (/1/), topology_status, message)
      call check(topology_status == BLOCK_TOPOLOGY_OK, 'threshold-boundary topology builds')
      if (.not. topology%ready) return

      moment_sum_mub = 1.0_dblprec
      direction_defined = .true.
      resultant_mub = 0.0_dblprec
      resultant_mub(1,1,1,1) = 0.85_dblprec                     ! strictly below
      resultant_mub(1,1,2,1) = 0.9_dblprec                      ! exactly at
      resultant_mub(1,1,3,1) = 0.95_dblprec                     ! strictly above
      resultant_mub(1,1,4,1) = 0.0_dblprec                      ! cancelled (see below)
      direction_defined(1,4,1) = .false.
      resultant_mub(1,1,5,1) = 0.9_dblprec - 1.0e-12_dblprec    ! roundoff-scale below

      call evaluate_polarization_gate(topology, resultant_mub, moment_sum_mub, &
         direction_defined, 0.9_dblprec, unsafe_block, status, diagnostic)
      call check(status == ADAPTIVE_HYBRID_OK, 'threshold-boundary gate evaluates: '//trim(diagnostic))
      call check(unsafe_block(1), 'ratio strictly below threshold is unsafe')
      call check(.not. unsafe_block(2), 'ratio exactly at threshold is safe')
      call check(.not. unsafe_block(3), 'ratio strictly above threshold is safe')
      call check(unsafe_block(4), 'undefined direction is unsafe')
      call check(unsafe_block(5), 'roundoff-scale below-threshold ratio is unsafe')
   end subroutine test_threshold_boundaries

   !> POL-CANCELLATION negative control: a block whose channel resultant is
   !> exactly zero (compensated single-channel moments) must never be judged
   !> safe to coarsen, and the gate must not attempt to normalize it. This is
   !> the defect the pre-fix code (no gate at all) would have let through.
   subroutine test_zero_resultant_never_normalizes()
      type(block_topology_type) :: topology
      integer(c_int) :: topology_status
      character(len=512) :: message
      real(c_double) :: cell(3,3)
      real(dblprec) :: resultant_mub(3,1,1,1), moment_sum_mub(1,1,1)
      logical :: direction_defined(1,1,1), unsafe_block(1)
      integer :: status
      character(len=256) :: diagnostic

      cell = identity_cell()
      call build_block_topology(topology, REGULAR_REPLICATED_CELL, 1, (/1,1,1/), 1, &
         (/1,1,1/), cell, (/1/), topology_status, message)
      call check(topology_status == BLOCK_TOPOLOGY_OK, 'cancellation topology builds')
      if (.not. topology%ready) return

      resultant_mub = 0.0_dblprec
      moment_sum_mub = 2.0_dblprec
      direction_defined = .false.

      call evaluate_polarization_gate(topology, resultant_mub, moment_sum_mub, &
         direction_defined, 0.9_dblprec, unsafe_block, status, diagnostic)
      call check(status == ADAPTIVE_HYBRID_OK, 'cancellation gate evaluates: '//trim(diagnostic))
      call check(unsafe_block(1), 'exactly-cancelled channel resultant remains unsafe')
   end subroutine test_zero_resultant_never_normalizes

   subroutine test_invalid_arguments()
      type(block_topology_type) :: topology, not_ready_topology
      integer(c_int) :: topology_status
      character(len=512) :: message
      real(c_double) :: cell(3,3)
      real(dblprec) :: resultant_mub(3,1,1,1), moment_sum_mub(1,1,1)
      real(dblprec) :: mismatched_resultant(3,1,2,1)
      logical :: direction_defined(1,1,1), unsafe_block(1), mismatched_unsafe(2)
      integer :: status
      character(len=256) :: diagnostic

      cell = identity_cell()
      call build_block_topology(topology, REGULAR_REPLICATED_CELL, 1, (/1,1,1/), 1, &
         (/1,1,1/), cell, (/1/), topology_status, message)
      call check(topology_status == BLOCK_TOPOLOGY_OK, 'invalid-argument topology builds')
      if (.not. topology%ready) return

      resultant_mub = 0.0_dblprec
      moment_sum_mub = 1.0_dblprec
      direction_defined = .true.

      call evaluate_polarization_gate(not_ready_topology, resultant_mub, moment_sum_mub, &
         direction_defined, 0.9_dblprec, unsafe_block, status, diagnostic)
      call check(status == ADAPTIVE_HYBRID_INVALID_TOPOLOGY, &
         'a not-ready topology is rejected: '//trim(diagnostic))

      call evaluate_polarization_gate(topology, mismatched_resultant, moment_sum_mub, &
         direction_defined, 0.9_dblprec, mismatched_unsafe, status, diagnostic)
      call check(status /= ADAPTIVE_HYBRID_OK, &
         'mismatched array shapes are rejected: '//trim(diagnostic))

      call evaluate_polarization_gate(topology, resultant_mub, moment_sum_mub, &
         direction_defined, 0.0_dblprec, unsafe_block, status, diagnostic)
      call check(status /= ADAPTIVE_HYBRID_OK, &
         'a zero threshold is rejected: '//trim(diagnostic))

      call evaluate_polarization_gate(topology, resultant_mub, moment_sum_mub, &
         direction_defined, 1.5_dblprec, unsafe_block, status, diagnostic)
      call check(status /= ADAPTIVE_HYBRID_OK, &
         'a threshold above one is rejected: '//trim(diagnostic))
   end subroutine test_invalid_arguments

   function identity_cell() result(cell)
      real(c_double) :: cell(3,3)

      cell = 0.0_c_double
      cell(1,1) = 1.0_c_double
      cell(2,2) = 1.0_c_double
      cell(3,3) = 1.0_c_double
   end function identity_cell

   subroutine check(condition, label)
      logical, intent(in) :: condition
      character(len=*), intent(in) :: label

      if (.not. condition) then
         failures = failures + 1
         write(*,'(a)') 'FAIL: '//trim(label)
      end if
   end subroutine check

end program test_polarization_gate
