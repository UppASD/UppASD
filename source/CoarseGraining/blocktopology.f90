!-------------------------------------------------------------------------------
! MODULE: BlockTopology
!> @brief Immutable regular spatial-block metadata for adaptive coarse graining.
!>
!> Construction is explicit: declaring a topology does not allocate anything.
!> A successfully built topology contains only scalar descriptors and plain
!> contiguous arrays.  Runtime masks, moments, fields, and state transitions
!> deliberately do not belong here.
!-------------------------------------------------------------------------------
module BlockTopology

   use, intrinsic :: iso_c_binding, only : c_double, c_int
   use, intrinsic :: iso_fortran_env, only : int64
   use, intrinsic :: ieee_arithmetic, only : ieee_is_finite

   implicit none
   private

   ! Geometry capability enum.  Only the regular replicated cell is currently
   ! constructible; the other values reserve explicit capability identities.
   integer(c_int), parameter, public :: GEOMETRY_UNSET = 0_c_int
   integer(c_int), parameter, public :: REGULAR_REPLICATED_CELL = 1_c_int
   integer(c_int), parameter, public :: EXPLICIT_DEVICE = 2_c_int
   integer(c_int), parameter, public :: IRREGULAR_SPATIAL_BINS = 3_c_int

   integer(c_int), parameter, public :: BLOCK_TOPOLOGY_OK = 0_c_int
   integer(c_int), parameter, public :: BLOCK_TOPOLOGY_INVALID_GEOMETRY = 1_c_int
   integer(c_int), parameter, public :: BLOCK_TOPOLOGY_INVALID_EXTENT = 2_c_int
   integer(c_int), parameter, public :: BLOCK_TOPOLOGY_INVALID_ATOM_COUNT = 3_c_int
   integer(c_int), parameter, public :: BLOCK_TOPOLOGY_INVALID_BLOCK_SHAPE = 4_c_int
   integer(c_int), parameter, public :: BLOCK_TOPOLOGY_INVALID_CHANNEL_MAP = 5_c_int
   integer(c_int), parameter, public :: BLOCK_TOPOLOGY_INVALID_CELL = 6_c_int
   integer(c_int), parameter, public :: BLOCK_TOPOLOGY_ALLOCATION_FAILED = 7_c_int
   integer(c_int), parameter, public :: BLOCK_TOPOLOGY_INTERNAL_ERROR = 8_c_int

   type, public :: block_topology_type
      logical :: ready = .false.
      integer(c_int) :: geometry_mode = GEOMETRY_UNSET
      integer(c_int) :: n_atoms = 0_c_int
      integer(c_int) :: n_spatial_blocks = 0_c_int
      integer(c_int) :: n_basis = 0_c_int
      integer(c_int) :: n_fft_channels_per_block = 0_c_int
      integer(c_int) :: n_fft_grid_channels = 0_c_int
      integer(c_int) :: n_dynamic_channels = 0_c_int
      integer(c_int) :: repetition_shape(3) = 0_c_int
      integer(c_int) :: block_shape(3) = 0_c_int
      integer(c_int) :: block_grid(3) = 0_c_int

      ! Physical vectors are columns: cell_vectors(:,r) is one replicated-cell
      ! step and block_vectors(:,r) is one spatial-block step.
      real(c_double) :: cell_vectors(3,3) = 0.0_c_double
      real(c_double) :: block_vectors(3,3) = 0.0_c_double

      ! Atom maps preserve the input atom numbering.  FFT channel is local to a
      ! spatial block; fft_grid_index is the existing global PME identity.
      integer(c_int), allocatable :: atom_to_block(:)
      integer(c_int), allocatable :: atom_to_basis(:)
      integer(c_int), allocatable :: atom_to_dynamic_channel(:)
      integer(c_int), allocatable :: atom_to_fft_channel(:)
      integer(c_int), allocatable :: atom_to_fft_grid_index(:)

      integer(c_int), allocatable :: basis_to_dynamic_channel(:)
      integer(c_int), allocatable :: basis_to_fft_channel(:)

      ! Zero-based CSR offsets make the arrays directly usable from C/C++.
      ! Fortran members of block b occupy offset(b)+1:offset(b+1).
      integer(c_int), allocatable :: block_atom_count(:)
      integer(c_int), allocatable :: block_atom_offset(:)
      integer(c_int), allocatable :: block_atoms(:)
      integer(c_int), allocatable :: block_grid_coordinate(:,:)

      integer(c_int), allocatable :: block_basis_population(:,:)
      integer(c_int), allocatable :: block_fft_channel_population(:,:)
      integer(c_int), allocatable :: block_dynamic_channel_population(:,:)

      real(c_double), allocatable :: block_center(:,:)
      real(c_double), allocatable :: block_volume(:)
   end type block_topology_type

   ! This object remains allocation-free until an owning feature explicitly
   ! passes it to build_block_topology.
   type(block_topology_type), public, save :: canonical_block_topology

   public :: build_block_topology
   public :: block_topology_has_allocations
   public :: regular_spatial_block_id
   public :: regular_fft_grid_channel

contains

   !> Existing FFT spatial-grid ordering: x is fastest, followed by y and z.
   pure integer(c_int) function regular_spatial_block_id(block_coordinate, block_grid) result(block_id)
      integer, intent(in) :: block_coordinate(3), block_grid(3)

      block_id = int(1 + block_coordinate(1) + block_grid(1) * &
         (block_coordinate(2) + block_grid(2) * block_coordinate(3)), c_int)
   end function regular_spatial_block_id

   !> Existing basis-resolved PME channel identity.
   pure integer(c_int) function regular_fft_grid_channel(basis_id, block_id, n_basis) result(channel_id)
      integer, intent(in) :: basis_id, block_id, n_basis

      channel_id = int(basis_id + n_basis * (block_id - 1), c_int)
   end function regular_fft_grid_channel

   logical function block_topology_has_allocations(topology)
      type(block_topology_type), intent(in) :: topology

      block_topology_has_allocations = allocated(topology%atom_to_block) .or. &
         allocated(topology%atom_to_basis) .or. &
         allocated(topology%atom_to_dynamic_channel) .or. &
         allocated(topology%atom_to_fft_channel) .or. &
         allocated(topology%atom_to_fft_grid_index) .or. &
         allocated(topology%basis_to_dynamic_channel) .or. &
         allocated(topology%basis_to_fft_channel) .or. &
         allocated(topology%block_atom_count) .or. &
         allocated(topology%block_atom_offset) .or. &
         allocated(topology%block_atoms) .or. &
         allocated(topology%block_grid_coordinate) .or. &
         allocated(topology%block_basis_population) .or. &
         allocated(topology%block_fft_channel_population) .or. &
         allocated(topology%block_dynamic_channel_population) .or. &
         allocated(topology%block_center) .or. allocated(topology%block_volume)
   end function block_topology_has_allocations

   !> Build a complete regular topology or return an allocation-free failure.
   !>
   !> atom_to_block/atom_to_basis/etc. are indexed by the canonical global
   !> atom index (geometry.f90's I0+I1*NA+I2*N1*NA+I3*N2*N1*NA, matching
   !> magnetizationinit.f90's loop nesting), the same numbering emom/ham%nlist
   !> use -- NOT the block-major traversal create_pme_macrocell_layout uses
   !> internally for its own, separate pme_cell_index array. The builder
   !> never renumbers an atom or modifies an atomistic array.
   subroutine build_block_topology(topology, geometry_mode, NA, repetitions, Natom, &
         block_shape, cell_vectors, basis_dynamic_channel, status, diagnostic)
      type(block_topology_type), intent(out) :: topology
      integer, intent(in) :: geometry_mode, NA, repetitions(3), Natom, block_shape(3)
      real(c_double), intent(in) :: cell_vectors(3,3)
      integer, intent(in) :: basis_dynamic_channel(:)
      integer(c_int), intent(out) :: status
      character(len=*), intent(out) :: diagnostic

      integer :: allocation_status
      integer :: atom, basis, block, channel, position
      integer :: ib1, ib2, ib3, i1, i2, i3
      integer :: grid(3), block_coordinate(3), cursor_value
      integer :: nblocks, n_dynamic_channels, expected_basis_population
      integer(int64) :: expected_atoms, nblocks_wide, fft_channels_wide
      integer(c_int), allocatable :: cursor(:)
      logical, allocatable :: seen_atom(:)
      real(c_double) :: volume

      status = BLOCK_TOPOLOGY_OK
      diagnostic = ''

      ! Capability rejection deliberately precedes every allocation, especially
      ! any basis/channel-shaped allocation in an NA=Natom device input.
      if (geometry_mode == EXPLICIT_DEVICE) then
         if (Natom > 0 .and. NA == Natom) then
            call fail(BLOCK_TOPOLOGY_INVALID_GEOMETRY, &
               'Adaptive coarse graining rejects EXPLICIT_DEVICE geometry with NA=Natom; '// &
               'coordinate-binned magnetic channels are not implemented and no NA-sized channel arrays were allocated')
         else
            call fail(BLOCK_TOPOLOGY_INVALID_GEOMETRY, &
               'Adaptive coarse graining rejects EXPLICIT_DEVICE geometry; coordinate-binned magnetic channels '// &
               'are not implemented and no topology channel arrays were allocated')
         end if
         return
      end if
      if (Natom > 1 .and. NA == Natom) then
         call fail(BLOCK_TOPOLOGY_INVALID_GEOMETRY, &
            'Adaptive coarse graining identifies NA=Natom as EXPLICIT_DEVICE geometry and rejects it; '// &
            'coordinate-binned magnetic channels are not implemented and no NA-sized channel arrays were allocated')
         return
      end if
      if (geometry_mode /= REGULAR_REPLICATED_CELL) then
         call fail(BLOCK_TOPOLOGY_INVALID_GEOMETRY, &
            'Adaptive coarse graining supports REGULAR_REPLICATED_CELL geometry only')
         return
      end if
      if (NA <= 0 .or. Natom <= 0 .or. any(repetitions <= 0)) then
         call fail(BLOCK_TOPOLOGY_INVALID_EXTENT, &
            'REGULAR_REPLICATED_CELL requires positive NA, Natom, N1, N2, and N3')
         return
      end if

      expected_atoms = int(NA, int64) * int(repetitions(1), int64) * &
         int(repetitions(2), int64) * int(repetitions(3), int64)
      if (expected_atoms /= int(Natom, int64)) then
         call fail(BLOCK_TOPOLOGY_INVALID_ATOM_COUNT, &
            'REGULAR_REPLICATED_CELL requires Natom=NA*N1*N2*N3')
         return
      end if
      if (any(block_shape <= 0)) then
         call fail(BLOCK_TOPOLOGY_INVALID_BLOCK_SHAPE, &
            'Adaptive coarse graining requires positive block dimensions')
         return
      end if
      if (any(mod(repetitions, block_shape) /= 0)) then
         call fail(BLOCK_TOPOLOGY_INVALID_BLOCK_SHAPE, &
            'Adaptive coarse graining requires block dimensions that divide N1, N2, and N3 exactly')
         return
      end if
      if (size(basis_dynamic_channel) /= NA) then
         call fail(BLOCK_TOPOLOGY_INVALID_CHANNEL_MAP, &
            'The basis-to-dynamical-channel map must contain exactly NA entries')
         return
      end if
      if (any(basis_dynamic_channel == 0) .or. any(basis_dynamic_channel < -1)) then
         call fail(BLOCK_TOPOLOGY_INVALID_CHANNEL_MAP, &
            'Dynamical-channel ids must be positive, or -1 for a nonmagnetic basis site')
         return
      end if
      n_dynamic_channels = maxval(basis_dynamic_channel)
      if (n_dynamic_channels <= 0) then
         call fail(BLOCK_TOPOLOGY_INVALID_CHANNEL_MAP, &
            'Adaptive coarse graining requires at least one nonempty magnetic dynamical channel')
         return
      end if
      do channel = 1, n_dynamic_channels
         if (count(basis_dynamic_channel == channel) == 0) then
            call fail(BLOCK_TOPOLOGY_INVALID_CHANNEL_MAP, &
               'Dynamical-channel ids must be contiguous and every magnetic channel must be nonempty')
            return
         end if
      end do
      if (.not. all(ieee_is_finite(cell_vectors))) then
         call fail(BLOCK_TOPOLOGY_INVALID_CELL, &
            'REGULAR_REPLICATED_CELL vectors must contain only finite values')
         return
      end if

      grid = repetitions / block_shape
      nblocks_wide = int(grid(1), int64) * int(grid(2), int64) * int(grid(3), int64)
      fft_channels_wide = int(NA, int64) * nblocks_wide
      if (nblocks_wide > int(huge(nblocks), int64) .or. &
          fft_channels_wide > int(huge(nblocks), int64)) then
         call fail(BLOCK_TOPOLOGY_INVALID_EXTENT, &
            'REGULAR_REPLICATED_CELL topology counts exceed the supported integer range')
         return
      end if
      nblocks = int(nblocks_wide)

      topology%block_vectors = cell_vectors
      topology%block_vectors(:,1) = real(block_shape(1), c_double) * cell_vectors(:,1)
      topology%block_vectors(:,2) = real(block_shape(2), c_double) * cell_vectors(:,2)
      topology%block_vectors(:,3) = real(block_shape(3), c_double) * cell_vectors(:,3)
      volume = abs(determinant3(topology%block_vectors))
      if (.not. ieee_is_finite(volume) .or. volume <= tiny(volume)) then
         call fail(BLOCK_TOPOLOGY_INVALID_CELL, &
            'REGULAR_REPLICATED_CELL block vectors must span a finite nonzero volume')
         return
      end if

      allocate(topology%atom_to_block(Natom), topology%atom_to_basis(Natom), &
         topology%atom_to_dynamic_channel(Natom), topology%atom_to_fft_channel(Natom), &
         topology%atom_to_fft_grid_index(Natom), topology%basis_to_dynamic_channel(NA), &
         topology%basis_to_fft_channel(NA), topology%block_atom_count(nblocks), &
         topology%block_atom_offset(nblocks + 1), topology%block_atoms(Natom), &
         topology%block_grid_coordinate(3,nblocks), &
         topology%block_basis_population(NA,nblocks), &
         topology%block_fft_channel_population(NA,nblocks), &
         topology%block_dynamic_channel_population(n_dynamic_channels,nblocks), &
         topology%block_center(3,nblocks), topology%block_volume(nblocks), &
         cursor(nblocks), seen_atom(Natom), stat=allocation_status)
      if (allocation_status /= 0) then
         call clear_topology(topology)
         call fail(BLOCK_TOPOLOGY_ALLOCATION_FAILED, &
            'Unable to allocate the validated regular block topology')
         return
      end if

      topology%atom_to_block = 0_c_int
      topology%atom_to_basis = 0_c_int
      topology%atom_to_dynamic_channel = -1_c_int
      topology%atom_to_fft_channel = 0_c_int
      topology%atom_to_fft_grid_index = 0_c_int
      topology%basis_to_dynamic_channel = int(basis_dynamic_channel, c_int)
      topology%basis_to_fft_channel = [(int(basis, c_int), basis=1,NA)]
      topology%block_atom_count = 0_c_int
      topology%block_atom_offset = 0_c_int
      topology%block_atoms = 0_c_int
      topology%block_grid_coordinate = 0_c_int
      topology%block_basis_population = 0_c_int
      topology%block_fft_channel_population = 0_c_int
      topology%block_dynamic_channel_population = 0_c_int
      topology%block_center = 0.0_c_double
      topology%block_volume = volume

      ! Atom numbering must match the canonical global atom index used
      ! everywhere else in the codebase (geometry.f90's I0+I1*NA+I2*N1*NA+
      ! I3*N2*N1*NA, basis fastest then I1, I2, I3 slowest -- see
      ! magnetizationinit.f90's identical loop nesting), NOT a block-major
      ! traversal counter: atom_to_block is queried downstream (selector
      ! misalignment, polarization-gate channel restriction, static/adaptive
      ! ownership) using real atom indices taken from emom/bond lists, which
      ! are always in this canonical order, with no translation layer. A
      ! block-major counter here silently mislabels which physical atom
      ! belongs to which spatial block whenever a block spans more than one
      ! cell along any axis (RCG-04G finding: previously masked because every
      ! prior fixture's state was uniform, random, or checked only
      ! atom-count aggregates, none of which are sensitive to *which* atom
      ! carries which block label).
      atom = 0
      do ib3 = 0, grid(3) - 1
         do ib2 = 0, grid(2) - 1
            do ib1 = 0, grid(1) - 1
               block_coordinate = (/ib1, ib2, ib3/)
               block = int(regular_spatial_block_id(block_coordinate, grid))
               topology%block_grid_coordinate(:,block) = int(block_coordinate, c_int)
               topology%block_center(:,block) = &
                  real(ib1 * block_shape(1), c_double) * cell_vectors(:,1) + &
                  real(ib2 * block_shape(2), c_double) * cell_vectors(:,2) + &
                  real(ib3 * block_shape(3), c_double) * cell_vectors(:,3) + &
                  0.5_c_double * sum(topology%block_vectors, dim=2)
               do i3 = ib3 * block_shape(3), (ib3 + 1) * block_shape(3) - 1
                  do i2 = ib2 * block_shape(2), (ib2 + 1) * block_shape(2) - 1
                     do i1 = ib1 * block_shape(1), (ib1 + 1) * block_shape(1) - 1
                        do basis = 1, NA
                           atom = basis + i1 * NA + i2 * repetitions(1) * NA + &
                              i3 * repetitions(2) * repetitions(1) * NA
                           topology%atom_to_block(atom) = int(block, c_int)
                           topology%atom_to_basis(atom) = int(basis, c_int)
                           topology%atom_to_dynamic_channel(atom) = &
                              int(basis_dynamic_channel(basis), c_int)
                           topology%atom_to_fft_channel(atom) = int(basis, c_int)
                           topology%atom_to_fft_grid_index(atom) = &
                              regular_fft_grid_channel(basis, block, NA)
                           topology%block_atom_count(block) = topology%block_atom_count(block) + 1_c_int
                           topology%block_basis_population(basis,block) = &
                              topology%block_basis_population(basis,block) + 1_c_int
                           topology%block_fft_channel_population(basis,block) = &
                              topology%block_fft_channel_population(basis,block) + 1_c_int
                           channel = basis_dynamic_channel(basis)
                           if (channel > 0) then
                              topology%block_dynamic_channel_population(channel,block) = &
                                 topology%block_dynamic_channel_population(channel,block) + 1_c_int
                           end if
                        end do
                     end do
                  end do
               end do
            end do
         end do
      end do
      if (any(topology%atom_to_block < 1) .or. any(topology%atom_to_block > nblocks)) then
         call clear_topology(topology)
         call fail(BLOCK_TOPOLOGY_INTERNAL_ERROR, &
            'Internal error: regular topology did not assign every atom')
         return
      end if

      topology%block_atom_offset(1) = 0_c_int
      do block = 1, nblocks
         topology%block_atom_offset(block + 1) = topology%block_atom_offset(block) + &
            topology%block_atom_count(block)
         cursor(block) = topology%block_atom_offset(block) + 1_c_int
      end do
      seen_atom = .false.
      do atom = 1, Natom
         block = topology%atom_to_block(atom)
         if (block < 1 .or. block > nblocks) then
            call clear_topology(topology)
            call fail(BLOCK_TOPOLOGY_INTERNAL_ERROR, &
               'Internal error: atom has an invalid spatial-block membership')
            return
         end if
         position = cursor(block)
         if (position < topology%block_atom_offset(block) + 1 .or. &
             position > topology%block_atom_offset(block + 1)) then
            call clear_topology(topology)
            call fail(BLOCK_TOPOLOGY_INTERNAL_ERROR, &
               'Internal error: spatial-block CSR membership overflow')
            return
         end if
         topology%block_atoms(position) = int(atom, c_int)
         cursor(block) = cursor(block) + 1_c_int
      end do
      do position = 1, Natom
         atom = topology%block_atoms(position)
         if (atom < 1 .or. atom > Natom .or. seen_atom(atom)) then
            call clear_topology(topology)
            call fail(BLOCK_TOPOLOGY_INTERNAL_ERROR, &
               'Internal error: CSR membership does not cover each atom exactly once')
            return
         end if
         seen_atom(atom) = .true.
      end do
      if (.not. all(seen_atom)) then
         call clear_topology(topology)
         call fail(BLOCK_TOPOLOGY_INTERNAL_ERROR, &
            'Internal error: CSR membership omits at least one atom')
         return
      end if

      expected_basis_population = product(block_shape)
      if (any(topology%block_basis_population /= expected_basis_population)) then
         call clear_topology(topology)
         call fail(BLOCK_TOPOLOGY_INTERNAL_ERROR, &
            'REGULAR_REPLICATED_CELL requires a uniform population of every basis site in every block')
         return
      end if
      do channel = 1, n_dynamic_channels
         cursor_value = count(basis_dynamic_channel == channel) * expected_basis_population
         if (cursor_value <= 0 .or. &
             any(topology%block_dynamic_channel_population(channel,:) /= cursor_value)) then
            call clear_topology(topology)
            call fail(BLOCK_TOPOLOGY_INTERNAL_ERROR, &
               'REGULAR_REPLICATED_CELL requires every magnetic channel to be nonempty and uniform')
            return
         end if
      end do

      topology%ready = .true.
      topology%geometry_mode = REGULAR_REPLICATED_CELL
      topology%n_atoms = int(Natom, c_int)
      topology%n_spatial_blocks = int(nblocks, c_int)
      topology%n_basis = int(NA, c_int)
      topology%n_fft_channels_per_block = int(NA, c_int)
      topology%n_fft_grid_channels = int(fft_channels_wide, c_int)
      topology%n_dynamic_channels = int(n_dynamic_channels, c_int)
      topology%repetition_shape = int(repetitions, c_int)
      topology%block_shape = int(block_shape, c_int)
      topology%block_grid = int(grid, c_int)
      topology%cell_vectors = cell_vectors

   contains

      subroutine fail(error_status, text)
         integer(c_int), intent(in) :: error_status
         character(len=*), intent(in) :: text

         status = error_status
         diagnostic = text
      end subroutine fail

   end subroutine build_block_topology

   pure real(c_double) function determinant3(matrix) result(determinant)
      real(c_double), intent(in) :: matrix(3,3)

      determinant = matrix(1,1) * (matrix(2,2) * matrix(3,3) - matrix(2,3) * matrix(3,2)) - &
         matrix(1,2) * (matrix(2,1) * matrix(3,3) - matrix(2,3) * matrix(3,1)) + &
         matrix(1,3) * (matrix(2,1) * matrix(3,2) - matrix(2,2) * matrix(3,1))
   end function determinant3

   subroutine clear_topology(topology)
      type(block_topology_type), intent(inout) :: topology

      if (allocated(topology%atom_to_block)) deallocate(topology%atom_to_block)
      if (allocated(topology%atom_to_basis)) deallocate(topology%atom_to_basis)
      if (allocated(topology%atom_to_dynamic_channel)) deallocate(topology%atom_to_dynamic_channel)
      if (allocated(topology%atom_to_fft_channel)) deallocate(topology%atom_to_fft_channel)
      if (allocated(topology%atom_to_fft_grid_index)) deallocate(topology%atom_to_fft_grid_index)
      if (allocated(topology%basis_to_dynamic_channel)) deallocate(topology%basis_to_dynamic_channel)
      if (allocated(topology%basis_to_fft_channel)) deallocate(topology%basis_to_fft_channel)
      if (allocated(topology%block_atom_count)) deallocate(topology%block_atom_count)
      if (allocated(topology%block_atom_offset)) deallocate(topology%block_atom_offset)
      if (allocated(topology%block_atoms)) deallocate(topology%block_atoms)
      if (allocated(topology%block_grid_coordinate)) deallocate(topology%block_grid_coordinate)
      if (allocated(topology%block_basis_population)) deallocate(topology%block_basis_population)
      if (allocated(topology%block_fft_channel_population)) deallocate(topology%block_fft_channel_population)
      if (allocated(topology%block_dynamic_channel_population)) deallocate(topology%block_dynamic_channel_population)
      if (allocated(topology%block_center)) deallocate(topology%block_center)
      if (allocated(topology%block_volume)) deallocate(topology%block_volume)

      topology%ready = .false.
      topology%geometry_mode = GEOMETRY_UNSET
      topology%n_atoms = 0_c_int
      topology%n_spatial_blocks = 0_c_int
      topology%n_basis = 0_c_int
      topology%n_fft_channels_per_block = 0_c_int
      topology%n_fft_grid_channels = 0_c_int
      topology%n_dynamic_channels = 0_c_int
      topology%repetition_shape = 0_c_int
      topology%block_shape = 0_c_int
      topology%block_grid = 0_c_int
      topology%cell_vectors = 0.0_c_double
      topology%block_vectors = 0.0_c_double
   end subroutine clear_topology

end module BlockTopology
