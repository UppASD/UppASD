program test_block_topology

   use, intrinsic :: iso_c_binding, only : c_double, c_int
   use BlockTopology

   implicit none

   integer :: failures

   failures = 0
   call check(.not. canonical_block_topology%ready, &
      'feature-off canonical topology is not ready')
   call check(.not. block_topology_has_allocations(canonical_block_topology), &
      'feature-off canonical topology has no allocations')

   call test_na_one()
   call test_na_two_and_distinct_maps()
   call test_na_greater_than_two()
   call test_non_cubic_blocks()
   call test_block_size_one()
   call test_skew_cell()
   call test_invalid_preflight()
   call test_explicit_device_rejection()

   call check(.not. block_topology_has_allocations(canonical_block_topology), &
      'local builders do not allocate the feature-off canonical topology')

   if (failures /= 0) then
      write(*,'(a,i0)') 'block topology tests failed: ', failures
      stop 1
   end if
   write(*,'(a)') 'block topology tests passed'

contains

   subroutine test_na_one()
      type(block_topology_type) :: topology
      integer(c_int) :: status
      character(len=512) :: message
      real(c_double) :: cell(3,3)

      cell = identity_cell()
      call build_block_topology(topology, REGULAR_REPLICATED_CELL, 1, (/1,1,1/), 1, &
         (/1,1,1/), cell, (/1/), status, message)
      call check(status == BLOCK_TOPOLOGY_OK .and. topology%ready, &
         'NA=1, 1x1x1 topology builds')
      if (.not. topology%ready) return
      call check(topology%n_spatial_blocks == 1 .and. topology%n_dynamic_channels == 1, &
         'NA=1 counts are exact')
      call check(all(topology%block_atom_offset == (/0,1/)) .and. &
         all(topology%block_atoms == (/1/)), 'NA=1 CSR membership is exact')
   end subroutine test_na_one

   subroutine test_na_two_and_distinct_maps()
      type(block_topology_type) :: topology
      integer(c_int) :: status
      character(len=512) :: message
      real(c_double) :: cell(3,3)

      cell = identity_cell()
      call build_block_topology(topology, REGULAR_REPLICATED_CELL, 2, (/2,1,1/), 4, &
         (/1,1,1/), cell, (/1,1/), status, message)
      call check(status == BLOCK_TOPOLOGY_OK, 'NA=2 topology builds')
      if (.not. topology%ready) return
      call check(all(topology%atom_to_basis == (/1,2,1,2/)), &
         'atom-to-basis map is basis resolved')
      call check(all(topology%atom_to_dynamic_channel == (/1,1,1,1/)), &
         'atom-to-dynamical-channel map can merge basis sites')
      call check(all(topology%atom_to_fft_channel == (/1,2,1,2/)), &
         'atom-to-FFT-channel map remains basis resolved')
      call check(all(topology%atom_to_fft_grid_index == (/1,2,3,4/)), &
         'global FFT-grid channels use the PME identity')
      call check(topology%n_basis == 2 .and. topology%n_fft_channels_per_block == 2 .and. &
         topology%n_dynamic_channels == 1, 'basis, FFT, and dynamical counts are separate')
   end subroutine test_na_two_and_distinct_maps

   subroutine test_na_greater_than_two()
      type(block_topology_type) :: topology
      integer(c_int) :: status
      character(len=512) :: message
      real(c_double) :: cell(3,3)

      cell = identity_cell()
      call build_block_topology(topology, REGULAR_REPLICATED_CELL, 4, (/2,1,1/), 8, &
         (/2,1,1/), cell, (/1,2,1,-1/), status, message)
      call check(status == BLOCK_TOPOLOGY_OK, 'NA>2 topology builds')
      if (.not. topology%ready) return
      call check(topology%n_dynamic_channels == 2, 'NA>2 channel count excludes nonmagnetic basis')
      call check(all(topology%block_dynamic_channel_population(:,1) == (/4,2/)), &
         'NA>2 magnetic channel populations are uniform and nonempty')
      call check(count(topology%atom_to_dynamic_channel == -1) == 2, &
         'nonmagnetic basis atoms retain channel id -1')
   end subroutine test_na_greater_than_two

   subroutine test_non_cubic_blocks()
      type(block_topology_type) :: topology
      integer(c_int) :: status
      character(len=512) :: message
      real(c_double) :: cell(3,3)
      integer :: atom, block, first_atom
      integer :: membership_count(48)

      cell = identity_cell()
      call build_block_topology(topology, REGULAR_REPLICATED_CELL, 1, (/4,6,2/), 48, &
         (/2,3,1/), cell, (/1/), status, message)
      call check(status == BLOCK_TOPOLOGY_OK, 'non-cubic block topology builds')
      if (.not. topology%ready) return
      call check(all(topology%block_grid == (/2,2,2/)) .and. &
         topology%n_spatial_blocks == 8, 'non-cubic block grid is exact')
      call check(all(topology%block_atom_count == 6), &
         'non-cubic blocks have uniform atom populations')
      membership_count = 0
      do atom = 1, size(topology%block_atoms)
         membership_count(topology%block_atoms(atom)) = &
            membership_count(topology%block_atoms(atom)) + 1
      end do
      call check(all(membership_count == 1), &
         'non-cubic CSR membership covers every atom exactly once')
      do block = 1, 8
         first_atom = 1 + 6 * (block - 1)
         call check(all(topology%atom_to_block(first_atom:first_atom+5) == block), &
            'spatial block ids follow the FFT x-fastest block grid')
      end do
      call check(regular_spatial_block_id((/1,1,1/),(/2,2,2/)) == 8, &
         'FFT spatial block id helper is x-fastest')
   end subroutine test_non_cubic_blocks

   subroutine test_block_size_one()
      type(block_topology_type) :: topology
      integer(c_int) :: status
      character(len=512) :: message
      real(c_double) :: cell(3,3)
      integer(c_int) :: expected_atoms(12), expected_blocks(12)
      integer :: atom

      cell = identity_cell()
      call build_block_topology(topology, REGULAR_REPLICATED_CELL, 3, (/2,2,1/), 12, &
         (/1,1,1/), cell, (/1,2,3/), status, message)
      call check(status == BLOCK_TOPOLOGY_OK, 'block-size-one topology builds')
      if (.not. topology%ready) return
      expected_atoms = [(int(atom,c_int), atom=1,12)]
      do atom = 1, 12
         expected_blocks(atom) = int(1 + (atom-1)/3, c_int)
      end do
      call check(all(topology%block_atoms == expected_atoms), &
         'block-size-one CSR preserves exact atom order')
      call check(all(topology%atom_to_fft_grid_index == expected_atoms), &
         'block-size-one FFT mapping is exact')
      call check(all(topology%atom_to_block == expected_blocks), &
         'block-size-one spatial mapping is exact')
   end subroutine test_block_size_one

   subroutine test_skew_cell()
      type(block_topology_type) :: topology
      integer(c_int) :: status
      character(len=512) :: message
      real(c_double) :: cell(3,3), expected_vectors(3,3), expected_center(3)

      cell(:,1) = (/2.0_c_double, 0.0_c_double, 0.0_c_double/)
      cell(:,2) = (/0.5_c_double, 3.0_c_double, 0.0_c_double/)
      cell(:,3) = (/0.25_c_double, 0.75_c_double, 4.0_c_double/)
      expected_vectors = cell
      expected_vectors(:,1) = 2.0_c_double * cell(:,1)
      expected_center = 0.5_c_double * sum(expected_vectors, dim=2)
      call build_block_topology(topology, REGULAR_REPLICATED_CELL, 1, (/2,1,1/), 2, &
         (/2,1,1/), cell, (/1/), status, message)
      call check(status == BLOCK_TOPOLOGY_OK, 'skew-cell topology builds')
      if (.not. topology%ready) return
      call check(maxval(abs(topology%cell_vectors-cell)) < 1.0e-14_c_double, &
         'skew replicated-cell vectors are stored')
      call check(maxval(abs(topology%block_vectors-expected_vectors)) < 1.0e-14_c_double, &
         'skew physical block vectors are stored')
      call check(abs(topology%block_volume(1)-48.0_c_double) < 1.0e-14_c_double, &
         'skew block volume is exact')
      call check(maxval(abs(topology%block_center(:,1)-expected_center)) < 1.0e-14_c_double, &
         'skew block center is exact')
   end subroutine test_skew_cell

   subroutine test_invalid_preflight()
      type(block_topology_type) :: topology
      integer(c_int) :: status
      character(len=512) :: message
      real(c_double) :: cell(3,3)

      cell = identity_cell()
      call build_block_topology(topology, REGULAR_REPLICATED_CELL, 1, (/2,2,1/), 3, &
         (/1,1,1/), cell, (/1/), status, message)
      call check(status == BLOCK_TOPOLOGY_INVALID_ATOM_COUNT .and. &
         .not. block_topology_has_allocations(topology), &
         'invalid atom count fails before allocation')

      call build_block_topology(topology, REGULAR_REPLICATED_CELL, 1, (/3,2,1/), 6, &
         (/2,1,1/), cell, (/1/), status, message)
      call check(status == BLOCK_TOPOLOGY_INVALID_BLOCK_SHAPE .and. &
         .not. block_topology_has_allocations(topology), &
         'nondivisible block shape fails before allocation')

      call build_block_topology(topology, REGULAR_REPLICATED_CELL, 1, (/2,2,1/), 4, &
         (/0,1,1/), cell, (/1/), status, message)
      call check(status == BLOCK_TOPOLOGY_INVALID_BLOCK_SHAPE .and. &
         .not. block_topology_has_allocations(topology), &
         'nonpositive block shape fails before allocation')

      call build_block_topology(topology, REGULAR_REPLICATED_CELL, 2, (/2,1,1/), 4, &
         (/1,1,1/), cell, (/-1,-1/), status, message)
      call check(status == BLOCK_TOPOLOGY_INVALID_CHANNEL_MAP .and. &
         .not. block_topology_has_allocations(topology), &
         'empty magnetic channel map fails before allocation')

      call build_block_topology(topology, REGULAR_REPLICATED_CELL, 2, (/2,1,1/), 4, &
         (/1,1,1/), cell, (/1,3/), status, message)
      call check(status == BLOCK_TOPOLOGY_INVALID_CHANNEL_MAP .and. &
         .not. block_topology_has_allocations(topology), &
         'empty channel id gap fails before allocation')
   end subroutine test_invalid_preflight

   subroutine test_explicit_device_rejection()
      type(block_topology_type) :: topology
      integer(c_int) :: status
      character(len=512) :: message
      real(c_double) :: cell(3,3)

      cell = identity_cell()
      call build_block_topology(topology, EXPLICIT_DEVICE, 6, (/1,1,1/), 6, &
         (/1,1,1/), cell, (/1,2,3,4,5,6/), status, message)
      call check(status == BLOCK_TOPOLOGY_INVALID_GEOMETRY, &
         'explicit-device geometry is rejected')
      call check(index(message,'EXPLICIT_DEVICE') > 0 .and. index(message,'NA=Natom') > 0, &
         'explicit-device diagnostic names geometry and NA=Natom')
      call check(index(message,'no NA-sized channel arrays were allocated') > 0 .and. &
         .not. block_topology_has_allocations(topology), &
         'explicit-device rejection allocates no topology channel arrays')
   end subroutine test_explicit_device_rejection

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

end program test_block_topology
