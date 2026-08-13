!> RCG-07 (F-08): dilation-complexity fixtures for StaticHybridOperator's
!> `rebuild_static_hybrid_ownership`.
!>
!> Two independent checks:
!>
!> 1. `test_local_stencil_matches_brute_force_reference` -- the local
!>    directional-stencil rebuild must produce a block-ownership map
!>    bitwise identical to a brute-force all-block-by-all-seed reference
!>    transcribed here (the *algorithm* the parent blueprint accepts:
!>    every block within the periodic Chebyshev-wrapped interaction-derived
!>    buffer width of any fine seed is atomistic), independent of
!>    `StaticHybridOperator`'s own implementation.
!> 2. `test_dilation_scales_linearly_not_quadratically` -- a negative
!>    control on wall-clock growth. With a fixed fine-block fraction and a
!>    fixed physical buffer width, the previous all-block-by-all-seed scan
!>    performed O(n_seeds*n_blocks) = O(n_blocks^2) work; the local stencil
!>    performs O(n_seeds*stencil_volume) = O(n_blocks) work. Measured on the
!>    pre-fix commit, the assertion below fails (observed ratio far above
!>    the quadratic-vs-linear discriminating threshold); on the fix it
!>    passes (see docs/RCG-07_CPU_ALGORITHMIC_SCALING_EVIDENCE.md).
program test_static_hybrid_dilation_scaling

   use, intrinsic :: iso_fortran_env, only : int64
   use Parameters, only : dblprec
   use BlockTopology
   use stiffness, only : coarse_material_type
   use CoarseTensorOperator
   use SmoothProjectedOperator
   use StaticHybridOperator

   implicit none

   integer :: failures

   failures = 0
   call test_local_stencil_matches_brute_force_reference()
   call test_dilation_scales_linearly_not_quadratically()

   if (failures /= 0) then
      write(*,'(a,i0)') 'static hybrid dilation scaling tests failed: ',failures
      stop 1
   end if
   write(*,'(a)') 'static hybrid dilation scaling tests passed'

contains

   subroutine test_local_stencil_matches_brute_force_reference()
      type(block_topology_type) :: topology
      type(static_hybrid_operator_type) :: hybrid
      logical, allocatable :: fine_mask(:), expected_atomistic(:)
      integer, parameter :: nblocks = 97
      integer :: status, seed, block
      character(len=512) :: message
      integer :: delta(3), periodic_delta(3)

      call build_chain_hybrid(nblocks,fine_stride=11,topology=topology,hybrid=hybrid, &
         fine_mask=fine_mask,status=status,message=message)
      call check(status == STATIC_HYBRID_OK,'brute-force reference fixture builds: '//trim(message))

      allocate(expected_atomistic(nblocks))
      expected_atomistic = fine_mask
      do seed = 1, nblocks
         if (.not. fine_mask(seed)) cycle
         do block = 1, nblocks
            delta = abs(topology%block_grid_coordinate(:,block) - &
               topology%block_grid_coordinate(:,seed))
            periodic_delta = min(delta,topology%block_grid-delta)
            if (all(periodic_delta <= hybrid%buffer_width_blocks)) then
               expected_atomistic(block) = .true.
            end if
         end do
      end do

      call check(count(fine_mask) > 1 .and. count(fine_mask) < nblocks, &
         'brute-force reference exercises more than one non-trivial seed')
      call check(all(hybrid%atomistic_block .eqv. expected_atomistic), &
         'local-stencil rebuild matches brute-force all-block-by-all-seed reference')
      call check(all(hybrid%coarse_block .eqv. .not. expected_atomistic), &
         'coarse_block remains the exact complement of atomistic_block')
   end subroutine test_local_stencil_matches_brute_force_reference

   !> Negative control (RCG remediation blueprint SS2.3): asserts a scaling
   !> bound that a quadratic implementation cannot meet. Quadratic
   !> prediction for an 8x block-count increase at fixed fine fraction is
   !> ~64x wall time; linear prediction is ~8x. The threshold below (20x)
   !> sits comfortably between the two, with margin for measurement noise,
   !> so it discriminates the two algorithms without being timing-flaky.
   subroutine test_dilation_scales_linearly_not_quadratically()
      integer, parameter :: n_sizes = 3
      integer, parameter :: sizes(n_sizes) = (/2048,8192,16384/)
      integer, parameter :: fine_stride = 8
      integer, parameter :: repetitions = 20
      integer, parameter :: trials = 5
      real(dblprec) :: elapsed(n_sizes), trial_seconds
      integer :: level, trial, status
      character(len=512) :: message

      ! Noise-robust timing: an OS scheduling stall or CPU contention from
      ! an unrelated process can only ever inflate a measured duration, so
      ! the minimum across several independent trials is a better estimate
      ! of the algorithm's true cost than a single measurement or a mean
      ! (which a single stalled trial would skew upward). Matches the
      ! repeated-sample idiom already used by this codebase's GPU adaptive
      ! benchmark (median/MAD over repetitions).
      do level = 1, n_sizes
         elapsed(level) = huge(1.0_dblprec)
         do trial = 1, trials
            call time_rebuild(sizes(level),fine_stride,repetitions,trial_seconds,status,message)
            call check(status == STATIC_HYBRID_OK,'scaling fixture builds: '//trim(message))
            elapsed(level) = min(elapsed(level),trial_seconds)
         end do
      end do

      write(*,'(a)') 'RCG-07 dilation scaling (blocks, min-of-trials rebuild seconds):'
      do level = 1, n_sizes
         write(*,'(a,i0,a,es12.5)') '  blocks=',sizes(level),' seconds=',elapsed(level)
      end do

      call check(elapsed(1) >= 0.0_dblprec .and. elapsed(n_sizes) >= 0.0_dblprec, &
         'measured wall-clock timings are non-negative')
      ! Guard the ratio against a degenerate near-zero denominator (a
      ! system_clock resolution floor, not an algorithmic property) by
      ! flooring it at the smallest resolvable duration used elsewhere in
      ! this module's timing.
      call check(elapsed(n_sizes) < 20.0_dblprec*max(elapsed(1),1.0d-7), &
         'dilation wall time grows sub-quadratically (linear-consistent) with block count')
   end subroutine test_dilation_scales_linearly_not_quadratically

   subroutine time_rebuild(nblocks,fine_stride,repetitions,mean_seconds,status,message)
      integer, intent(in) :: nblocks, fine_stride, repetitions
      real(dblprec), intent(out) :: mean_seconds
      integer, intent(out) :: status
      character(len=*), intent(out) :: message
      type(block_topology_type) :: topology
      type(static_hybrid_operator_type) :: hybrid
      logical, allocatable :: fine_mask(:)
      integer(int64) :: tick0, tick1, tick_rate
      integer :: rep, rebuild_status
      character(len=512) :: rebuild_message

      call build_chain_hybrid(nblocks,fine_stride,topology,hybrid,fine_mask,status,message)
      if (status /= STATIC_HYBRID_OK) return

      call system_clock(count=tick0,count_rate=tick_rate)
      do rep = 1, repetitions
         call rebuild_static_hybrid_ownership(hybrid,topology,fine_mask,rebuild_status,rebuild_message)
      end do
      call system_clock(count=tick1)
      if (rebuild_status /= STATIC_HYBRID_OK) then
         status = rebuild_status
         message = 'repeated rebuild failed: '//trim(rebuild_message)
         return
      end if
      if (tick_rate > 0_int64) then
         mean_seconds = real(tick1-tick0,dblprec)/real(tick_rate,dblprec)/real(repetitions,dblprec)
      else
         mean_seconds = 0.0_dblprec
      end if
      status = STATIC_HYBRID_OK
   end subroutine time_rebuild

   !> Builds a periodic 1D chain of `nblocks` single-atom blocks (block
   !> size 1, so buffer_width_blocks==(1,1,1) from the fixed nearest-
   !> neighbour bond distance regardless of `nblocks`) with every
   !> `fine_stride`-th block marked as a fine seed. Uses allocatable arrays
   !> throughout (RCG-06A: no O(nblocks) automatic/stack arrays), so
   !> `nblocks` can scale into the tens of thousands for the scaling
   !> fixture above.
   subroutine build_chain_hybrid(nblocks,fine_stride,topology,hybrid,fine_mask,status,message)
      integer, intent(in) :: nblocks, fine_stride
      type(block_topology_type), intent(out) :: topology
      type(static_hybrid_operator_type), intent(out) :: hybrid
      logical, allocatable, intent(out) :: fine_mask(:)
      integer, intent(out) :: status
      character(len=*), intent(out) :: message
      type(coarse_tensor_operator_type) :: tensor
      type(smooth_projected_operator_type) :: projection
      type(coarse_material_type) :: material
      type(coarse_operator_options_type) :: options
      integer, allocatable :: bonds(:,:)
      real(dblprec), allocatable :: displacement(:,:), matrix(:,:,:), coordinate(:,:), moment(:)
      real(dblprec), parameter :: a = 2.0d-10, pair_j = 3.0d-21
      real(dblprec) :: cell(3,3), exchange(3,3), dmi(3,3)
      integer :: atom

      cell = 0.0_dblprec
      cell(1,1) = a; cell(2,2) = a; cell(3,3) = a
      call build_block_topology(topology,REGULAR_REPLICATED_CELL,1,(/nblocks,1,1/), &
         nblocks,(/1,1,1/),cell,(/1/),status,message)
      if (status /= BLOCK_TOPOLOGY_OK) then
         status = STATIC_HYBRID_INVALID_TOPOLOGY
         return
      end if

      exchange = 0.0_dblprec
      exchange(1,1) = pair_j/(2.0_dblprec*a)
      dmi = 0.0_dblprec
      material%ready = .true.
      material%n_basis = 1
      material%n_channels = 1
      material%cell_volume_m3 = a**3
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
      call setup_coarse_tensor_operator(tensor,topology,material,options,status,message)
      if (status /= COARSE_TENSOR_OK) return

      allocate(bonds(2,nblocks),displacement(3,nblocks),matrix(3,3,nblocks), &
         coordinate(3,nblocks),moment(nblocks))
      moment = 1.7_dblprec
      do atom = 1, nblocks
         coordinate(:,atom) = 0.0_dblprec
         bonds(:,atom) = (/atom,modulo(atom,nblocks)+1/)
         displacement(:,atom) = (/a,0.0_dblprec,0.0_dblprec/)
         matrix(:,:,atom) = 0.0_dblprec
         matrix(1,1,atom) = pair_j; matrix(2,2,atom) = pair_j; matrix(3,3,atom) = pair_j
      end do
      call setup_smooth_projected_operator(projection,topology,coordinate,moment,status,message)
      if (status /= SMOOTH_PROJECTED_OK) return

      allocate(fine_mask(nblocks))
      fine_mask = .false.
      do atom = 1, nblocks, fine_stride
         fine_mask(atom) = .true.
      end do

      call setup_static_hybrid_operator(hybrid,topology,tensor,projection,fine_mask,bonds, &
         displacement,0,status,message)
   end subroutine build_chain_hybrid

   subroutine check(condition,message)
      logical, intent(in) :: condition
      character(len=*), intent(in) :: message

      if (condition) then
         write(*,'(a,a)') 'PASS: ',trim(message)
      else
         write(*,'(a,a)') 'FAIL: ',trim(message)
         failures = failures+1
      end if
   end subroutine check

end program test_static_hybrid_dilation_scaling
