!-------------------------------------------------------------------------------
! MODULE: SmoothProjectedOperator
!> @brief Double-precision CPU reference for smooth atomistic energy projection.
!>
!> Coarse variables live at the integer nodes of the regular block grid.  The
!> caller supplies each atom position in those global fractional block
!> coordinates (a block centre is an integer node).  Eight periodic trilinear
!> weights form v_a = sum_I w_aI M_I, followed by s_a = v_a/|v_a|.
!>
!> The exact normalization Jacobian is
!>
!>   d s_a / d M_I = w_aI (I - s_a s_a^T) / |v_a|.
!>
!> Consequently the field restriction is its moment-weighted adjoint:
!>
!>   B_I = sum_a (mu_a/mu_I) w_aI/|v_a|
!>                    (I - s_a s_a^T) B_a.
!>
!> Bilinear bonds use E_ij = -s_i^T K_ij s_j, with each physical bond supplied
!> once.  A general K represents scalar exchange, antisymmetric DMI, and tensor
!> exchange.  Bond indices and matrices are read-only; no atom or interaction
!> array is reordered, folded, or modified.  This is deliberately a clear CPU
!> reference implementation, not a production-performance or GPU path.
!-------------------------------------------------------------------------------
module SmoothProjectedOperator

   use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
   use Parameters, only : dblprec
   use BlockTopology, only : block_topology_type, REGULAR_REPLICATED_CELL, &
      regular_spatial_block_id
   use stiffness, only : COARSE_MUB_SI

   implicit none
   private

   integer, parameter, public :: SMOOTH_PROJECTED_OK = 0
   integer, parameter, public :: SMOOTH_PROJECTED_INVALID_TOPOLOGY = 1
   integer, parameter, public :: SMOOTH_PROJECTED_INVALID_COORDINATES = 2
   integer, parameter, public :: SMOOTH_PROJECTED_INVALID_MOMENTS = 3
   integer, parameter, public :: SMOOTH_PROJECTED_INVALID_STATE = 4
   integer, parameter, public :: SMOOTH_PROJECTED_ZERO_INTERPOLANT = 5
   integer, parameter, public :: SMOOTH_PROJECTED_INVALID_BONDS = 6
   integer, parameter, public :: SMOOTH_PROJECTED_PERIODIC = 1

   type, public :: smooth_projected_operator_type
      logical :: ready = .false.
      integer :: boundary_mode(3) = SMOOTH_PROJECTED_PERIODIC
      integer :: n_atoms = 0
      integer :: n_blocks = 0
      integer :: n_channels = 0
      integer :: block_grid(3) = 0
      real(dblprec) :: normalization_floor = 1.0d-12
      integer, allocatable :: atom_channel(:)
      integer, allocatable :: stencil_block(:,:)       ! (8,atom)
      real(dblprec), allocatable :: shape_weight(:,:)  ! (8,atom)
      real(dblprec), allocatable :: atom_moment_mub(:)
      real(dblprec), allocatable :: coarse_moment_mub(:,:) ! (channel,block)
      ! RCG-06A (F-13/F-17): persistent per-operator scratch, sized once in
      ! setup_smooth_projected_operator, replacing the natoms-sized automatic
      ! (stack) arrays previously declared inside prolongate_smooth_directions,
      ! restrict_projected_atomistic_field, and evaluate_projected_bilinear.
      ! Written and consumed within a single call only. See
      ! coarse_tensor_operator_type's matching comment for the OpenMP caveat.
      real(dblprec), allocatable :: scratch_local_norm(:)
      real(dblprec), allocatable :: scratch_restrict_direction(:,:)
      real(dblprec), allocatable :: scratch_restrict_raw_norm(:)
      real(dblprec), allocatable :: scratch_bilinear_direction(:,:)
      real(dblprec), allocatable :: scratch_bilinear_raw_norm(:)
      real(dblprec), allocatable :: scratch_bilinear_fine_field(:,:)
   end type smooth_projected_operator_type

   public :: setup_smooth_projected_operator
   public :: prolongate_smooth_directions
   public :: restrict_projected_atomistic_field
   public :: evaluate_projected_bilinear

contains

   !> Precompute periodic trilinear stencils without changing atom ordering.
   subroutine setup_smooth_projected_operator(operator, topology, &
         atom_fractional_block_coordinate, atom_moment_mub, status, diagnostic, &
         boundary_mode, normalization_floor)
      type(smooth_projected_operator_type), intent(out) :: operator
      type(block_topology_type), intent(in) :: topology
      real(dblprec), intent(in) :: atom_fractional_block_coordinate(:,:)
      real(dblprec), intent(in) :: atom_moment_mub(:)
      integer, intent(out) :: status
      character(len=*), intent(out) :: diagnostic
      integer, intent(in), optional :: boundary_mode(3)
      real(dblprec), intent(in), optional :: normalization_floor

      integer :: atom, channel, block, corner, dx, dy, dz
      integer :: base(3), coordinate(3)
      real(dblprec) :: fraction(3), wx, wy, wz

      operator%ready = .false.
      status = SMOOTH_PROJECTED_INVALID_TOPOLOGY
      diagnostic = ''

      if (.not. topology%ready .or. &
          topology%geometry_mode /= REGULAR_REPLICATED_CELL) then
         diagnostic = 'Smooth projection requires a ready REGULAR_REPLICATED_CELL topology'
         return
      end if
      if (present(boundary_mode)) then
         if (any(boundary_mode /= SMOOTH_PROJECTED_PERIODIC)) then
            diagnostic = 'Smooth projection currently supports explicit periodic boundaries only'
            return
         end if
      end if
      if (size(atom_fractional_block_coordinate,1) /= 3 .or. &
          size(atom_fractional_block_coordinate,2) /= topology%n_atoms .or. &
          .not. all(ieee_is_finite(atom_fractional_block_coordinate)) .or. &
          any(abs(atom_fractional_block_coordinate) > &
              0.5_dblprec*real(huge(0),dblprec))) then
         status = SMOOTH_PROJECTED_INVALID_COORDINATES
         diagnostic = 'Atom fractional block coordinates must be finite, representable, and have shape (3,natoms)'
         return
      end if
      if (size(atom_moment_mub) /= topology%n_atoms .or. &
          .not. all(ieee_is_finite(atom_moment_mub)) .or. &
          any(atom_moment_mub < 0.0_dblprec)) then
         status = SMOOTH_PROJECTED_INVALID_MOMENTS
         diagnostic = 'Atom moments must be finite, nonnegative, and have length natoms'
         return
      end if
      do atom = 1, topology%n_atoms
         channel = topology%atom_to_dynamic_channel(atom)
         if (channel > 0 .and. atom_moment_mub(atom) <= 0.0_dblprec) then
            status = SMOOTH_PROJECTED_INVALID_MOMENTS
            diagnostic = 'Every magnetic atom must have a strictly positive moment'
            return
         end if
      end do
      if (present(normalization_floor)) then
         if (.not. ieee_is_finite(normalization_floor) .or. &
             normalization_floor <= 0.0_dblprec) then
            status = SMOOTH_PROJECTED_INVALID_STATE
            diagnostic = 'The interpolation normalization floor must be finite and positive'
            return
         end if
         operator%normalization_floor = normalization_floor
      end if

      operator%n_atoms = topology%n_atoms
      operator%n_blocks = topology%n_spatial_blocks
      operator%n_channels = topology%n_dynamic_channels
      operator%block_grid = topology%block_grid
      allocate(operator%atom_channel(operator%n_atoms), &
         operator%stencil_block(8,operator%n_atoms), &
         operator%shape_weight(8,operator%n_atoms), &
         operator%atom_moment_mub(operator%n_atoms), &
         operator%coarse_moment_mub(operator%n_channels,operator%n_blocks))
      operator%atom_channel = topology%atom_to_dynamic_channel
      operator%atom_moment_mub = atom_moment_mub
      operator%coarse_moment_mub = 0.0_dblprec
      allocate(operator%scratch_local_norm(operator%n_atoms), &
         operator%scratch_restrict_direction(3,operator%n_atoms), &
         operator%scratch_restrict_raw_norm(operator%n_atoms), &
         operator%scratch_bilinear_direction(3,operator%n_atoms), &
         operator%scratch_bilinear_raw_norm(operator%n_atoms), &
         operator%scratch_bilinear_fine_field(3,operator%n_atoms))

      do atom = 1, operator%n_atoms
         channel = operator%atom_channel(atom)
         if (channel > 0) then
            block = topology%atom_to_block(atom)
            operator%coarse_moment_mub(channel,block) = &
               operator%coarse_moment_mub(channel,block) + atom_moment_mub(atom)
         end if

         base = floor(atom_fractional_block_coordinate(:,atom))
         fraction = atom_fractional_block_coordinate(:,atom) - real(base,dblprec)
         corner = 0
         do dz = 0, 1
            wz = merge(1.0_dblprec-fraction(3),fraction(3),dz == 0)
            do dy = 0, 1
               wy = merge(1.0_dblprec-fraction(2),fraction(2),dy == 0)
               do dx = 0, 1
                  wx = merge(1.0_dblprec-fraction(1),fraction(1),dx == 0)
                  corner = corner + 1
                  coordinate = modulo(base + (/dx,dy,dz/),operator%block_grid)
                  operator%stencil_block(corner,atom) = &
                     regular_spatial_block_id(coordinate,operator%block_grid)
                  operator%shape_weight(corner,atom) = wx*wy*wz
               end do
            end do
         end do
      end do

      if (any(operator%coarse_moment_mub <= 0.0_dblprec)) then
         status = SMOOTH_PROJECTED_INVALID_MOMENTS
         diagnostic = 'Every block/channel variable requires a positive coarse moment'
         return
      end if
      if (maxval(abs(sum(operator%shape_weight,dim=1)-1.0_dblprec)) > &
          64.0_dblprec*epsilon(1.0_dblprec)) then
         status = SMOOTH_PROJECTED_INVALID_COORDINATES
         diagnostic = 'Internal trilinear stencil failed the partition-of-unity check'
         return
      end if

      operator%ready = .true.
      status = SMOOTH_PROJECTED_OK
   end subroutine setup_smooth_projected_operator

   !> Apply trilinear interpolation and exact pointwise spin normalization.
   subroutine prolongate_smooth_directions(operator, coarse_variable, &
         atom_direction, status, diagnostic, raw_norm)
      type(smooth_projected_operator_type), intent(inout) :: operator
      real(dblprec), intent(in) :: coarse_variable(:,:,:)
      real(dblprec), intent(out) :: atom_direction(:,:)
      integer, intent(out) :: status
      character(len=*), intent(out) :: diagnostic
      real(dblprec), intent(out), optional :: raw_norm(:)

      call prolongate_with_norm(operator%ready,operator%n_atoms,operator%n_channels, &
         operator%n_blocks,operator%atom_channel,operator%stencil_block, &
         operator%shape_weight,operator%normalization_floor,coarse_variable, &
         atom_direction,operator%scratch_local_norm,status,diagnostic)
      if (present(raw_norm)) then
         if (size(raw_norm) == operator%n_atoms) then
            raw_norm = operator%scratch_local_norm
         else if (status == SMOOTH_PROJECTED_OK) then
            status = SMOOTH_PROJECTED_INVALID_STATE
            diagnostic = 'Raw interpolation norms must have length natoms'
         end if
      end if
   end subroutine prolongate_smooth_directions

   !> Apply the moment-weighted adjoint of the normalized prolongation.
   subroutine restrict_projected_atomistic_field(operator, coarse_variable, &
         atom_field_t, coarse_field_t, status, diagnostic)
      type(smooth_projected_operator_type), intent(inout) :: operator
      real(dblprec), intent(in) :: coarse_variable(:,:,:)
      real(dblprec), intent(in) :: atom_field_t(:,:)
      real(dblprec), intent(out) :: coarse_field_t(:,:,:)
      integer, intent(out) :: status
      character(len=*), intent(out) :: diagnostic

      if (size(atom_field_t,1) /= 3 .or. &
          size(atom_field_t,2) /= operator%n_atoms .or. &
          .not. all(ieee_is_finite(atom_field_t))) then
         status = SMOOTH_PROJECTED_INVALID_STATE
         diagnostic = 'Atomistic field must be finite with shape (3,natoms)'
         return
      end if
      call prolongate_with_norm(operator%ready,operator%n_atoms,operator%n_channels, &
         operator%n_blocks,operator%atom_channel,operator%stencil_block, &
         operator%shape_weight,operator%normalization_floor,coarse_variable, &
         operator%scratch_restrict_direction,operator%scratch_restrict_raw_norm, &
         status,diagnostic)
      if (status /= SMOOTH_PROJECTED_OK) return
      call restrict_with_state(operator,operator%scratch_restrict_direction, &
         operator%scratch_restrict_raw_norm,atom_field_t,coarse_field_t,status,diagnostic)
   end subroutine restrict_projected_atomistic_field

   !> Evaluate read-only unique-pair bilinear bonds and their projected field.
   subroutine evaluate_projected_bilinear(operator, coarse_variable, bond_atom, &
         bond_matrix_j, energy_j, coarse_field_t, status, diagnostic, &
         atom_direction, atom_field_t)
      type(smooth_projected_operator_type), intent(inout) :: operator
      real(dblprec), intent(in) :: coarse_variable(:,:,:)
      integer, intent(in) :: bond_atom(:,:)
      real(dblprec), intent(in) :: bond_matrix_j(:,:,:)
      real(dblprec), intent(out) :: energy_j
      real(dblprec), intent(out) :: coarse_field_t(:,:,:)
      integer, intent(out) :: status
      character(len=*), intent(out) :: diagnostic
      real(dblprec), intent(out), optional :: atom_direction(:,:)
      real(dblprec), intent(out), optional :: atom_field_t(:,:)

      integer :: atom_i, atom_j, bond, nbonds

      energy_j = 0.0_dblprec
      if (.not. operator%ready) then
         status = SMOOTH_PROJECTED_INVALID_STATE
         diagnostic = 'Cannot evaluate an uninitialized smooth projected operator'
         return
      end if
      nbonds = size(bond_atom,2)
      if (size(bond_atom,1) /= 2 .or. size(bond_matrix_j,1) /= 3 .or. &
          size(bond_matrix_j,2) /= 3 .or. size(bond_matrix_j,3) /= nbonds .or. &
          .not. all(ieee_is_finite(bond_matrix_j))) then
         status = SMOOTH_PROJECTED_INVALID_BONDS
         diagnostic = 'Bilinear input requires bond_atom(2,nbond) and finite K(3,3,nbond)'
         return
      end if
      if (nbonds > 0) then
         if (any(bond_atom < 1) .or. any(bond_atom > operator%n_atoms) .or. &
             any(bond_atom(1,:) == bond_atom(2,:))) then
            status = SMOOTH_PROJECTED_INVALID_BONDS
            diagnostic = 'Every bilinear bond must join two distinct magnetic atom indices'
            return
         end if
         do bond = 1, nbonds
            if (operator%atom_channel(bond_atom(1,bond)) <= 0 .or. &
                operator%atom_channel(bond_atom(2,bond)) <= 0) then
               status = SMOOTH_PROJECTED_INVALID_BONDS
               diagnostic = 'Bilinear bonds may reference magnetic atoms only'
               return
            end if
         end do
      end if

      associate (direction => operator%scratch_bilinear_direction, &
            raw_norm => operator%scratch_bilinear_raw_norm, &
            fine_field => operator%scratch_bilinear_fine_field)
      call prolongate_with_norm(operator%ready,operator%n_atoms,operator%n_channels, &
         operator%n_blocks,operator%atom_channel,operator%stencil_block, &
         operator%shape_weight,operator%normalization_floor,coarse_variable, &
         direction,raw_norm,status,diagnostic)
      if (status /= SMOOTH_PROJECTED_OK) return
      fine_field = 0.0_dblprec
      do bond = 1, nbonds
         atom_i = bond_atom(1,bond)
         atom_j = bond_atom(2,bond)
         energy_j = energy_j - dot_product(direction(:,atom_i), &
            matmul(bond_matrix_j(:,:,bond),direction(:,atom_j)))
         fine_field(:,atom_i) = fine_field(:,atom_i) + &
            matmul(bond_matrix_j(:,:,bond),direction(:,atom_j)) / &
            (COARSE_MUB_SI*operator%atom_moment_mub(atom_i))
         fine_field(:,atom_j) = fine_field(:,atom_j) + &
            matmul(transpose(bond_matrix_j(:,:,bond)),direction(:,atom_i)) / &
            (COARSE_MUB_SI*operator%atom_moment_mub(atom_j))
      end do

      call restrict_with_state(operator,direction,raw_norm,fine_field, &
         coarse_field_t,status,diagnostic)
      if (status /= SMOOTH_PROJECTED_OK) return
      if (present(atom_direction)) then
         if (any(shape(atom_direction) /= shape(direction))) then
            status = SMOOTH_PROJECTED_INVALID_STATE
            diagnostic = 'Optional atom directions must have shape (3,natoms)'
            return
         end if
         atom_direction = direction
      end if
      if (present(atom_field_t)) then
         if (any(shape(atom_field_t) /= shape(fine_field))) then
            status = SMOOTH_PROJECTED_INVALID_STATE
            diagnostic = 'Optional atom fields must have shape (3,natoms)'
            return
         end if
         atom_field_t = fine_field
      end if
      end associate
   end subroutine evaluate_projected_bilinear

   !> RCG-06A: takes explicit operator fields rather than the whole operator
   !> so that a caller may pass one of operator's own scratch arrays as
   !> `atom_direction`/`raw_norm` without violating Fortran's actual-argument
   !> aliasing rules -- see coarsetensoroperator.f90's matching comment on
   !> physical_forward_gradient for the underlying hazard this avoids.
   subroutine prolongate_with_norm(ready,n_atoms,n_channels,n_blocks,atom_channel, &
         stencil_block,shape_weight,normalization_floor,coarse_variable,atom_direction, &
         raw_norm,status,diagnostic)
      logical, intent(in) :: ready
      integer, intent(in) :: n_atoms, n_channels, n_blocks
      integer, intent(in) :: atom_channel(:)
      integer, intent(in) :: stencil_block(:,:)
      real(dblprec), intent(in) :: shape_weight(:,:)
      real(dblprec), intent(in) :: normalization_floor
      real(dblprec), intent(in) :: coarse_variable(:,:,:)
      real(dblprec), intent(out) :: atom_direction(:,:)
      real(dblprec), intent(out) :: raw_norm(:)
      integer, intent(out) :: status
      character(len=*), intent(out) :: diagnostic

      integer :: atom, channel, corner, block
      real(dblprec) :: interpolated(3)

      status = SMOOTH_PROJECTED_INVALID_STATE
      diagnostic = ''
      if (.not. ready) then
         diagnostic = 'Cannot prolong with an uninitialized smooth projected operator'
         return
      end if
      if (size(coarse_variable,1) /= 3 .or. &
          size(coarse_variable,2) /= n_channels .or. &
          size(coarse_variable,3) /= n_blocks .or. &
          .not. all(ieee_is_finite(coarse_variable))) then
         diagnostic = 'Coarse variables must be finite with shape (3,nchannels,nblocks)'
         return
      end if
      if (size(atom_direction,1) /= 3 .or. &
          size(atom_direction,2) /= n_atoms .or. &
          size(raw_norm) /= n_atoms) then
         diagnostic = 'Prolongation outputs must have atom-compatible shapes'
         return
      end if

      atom_direction = 0.0_dblprec
      raw_norm = 0.0_dblprec
      do atom = 1, n_atoms
         channel = atom_channel(atom)
         if (channel <= 0) cycle
         interpolated = 0.0_dblprec
         do corner = 1, 8
            block = stencil_block(corner,atom)
            interpolated = interpolated + shape_weight(corner,atom) * &
               coarse_variable(:,channel,block)
         end do
         raw_norm(atom) = sqrt(sum(interpolated*interpolated))
         if (.not. ieee_is_finite(raw_norm(atom)) .or. &
             raw_norm(atom) <= normalization_floor) then
            status = SMOOTH_PROJECTED_ZERO_INTERPOLANT
            diagnostic = 'Normalized prolongation encountered a zero or cancelling interpolant'
            return
         end if
         atom_direction(:,atom) = interpolated/raw_norm(atom)
      end do
      status = SMOOTH_PROJECTED_OK
   end subroutine prolongate_with_norm

   subroutine restrict_with_state(operator,atom_direction,raw_norm,atom_field_t, &
         coarse_field_t,status,diagnostic)
      type(smooth_projected_operator_type), intent(in) :: operator
      real(dblprec), intent(in) :: atom_direction(:,:), raw_norm(:)
      real(dblprec), intent(in) :: atom_field_t(:,:)
      real(dblprec), intent(out) :: coarse_field_t(:,:,:)
      integer, intent(out) :: status
      character(len=*), intent(out) :: diagnostic

      integer :: atom, channel, corner, block
      real(dblprec) :: tangent_field(3), factor

      if (size(coarse_field_t,1) /= 3 .or. &
          size(coarse_field_t,2) /= operator%n_channels .or. &
          size(coarse_field_t,3) /= operator%n_blocks) then
         status = SMOOTH_PROJECTED_INVALID_STATE
         diagnostic = 'Coarse field must have shape (3,nchannels,nblocks)'
         return
      end if
      coarse_field_t = 0.0_dblprec
      do atom = 1, operator%n_atoms
         channel = operator%atom_channel(atom)
         if (channel <= 0) cycle
         tangent_field = atom_field_t(:,atom) - atom_direction(:,atom) * &
            dot_product(atom_direction(:,atom),atom_field_t(:,atom))
         do corner = 1, 8
            block = operator%stencil_block(corner,atom)
            factor = operator%atom_moment_mub(atom) * &
               operator%shape_weight(corner,atom) / &
               (operator%coarse_moment_mub(channel,block)*raw_norm(atom))
            coarse_field_t(:,channel,block) = &
               coarse_field_t(:,channel,block) + factor*tangent_field
         end do
      end do
      status = SMOOTH_PROJECTED_OK
      diagnostic = ''
   end subroutine restrict_with_state

end module SmoothProjectedOperator
