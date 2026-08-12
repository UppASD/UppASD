!-------------------------------------------------------------------------------
! MODULE: StaticHybridOperator
!> @brief Conservative static fine/buffer/coarse CPU reference dispatch.
!>
!> Fine seeds are dilated periodically by the largest active short-range
!> interaction radius plus a caller-selected safety layer.  Fine and buffer
!> atoms are real atomistic degrees of freedom.  Atoms in coarse blocks are
!> smooth-prolongated ghosts whenever an owned atomistic bond crosses the
!> interface.  The ghost reaction is returned with the exact normalized
!> prolongation adjoint.
!>
!> Short-range ownership is deliberately non-overlapping:
!>
!> * an atomistic bond is owned iff at least one endpoint is fine/buffer;
!> * a tensor gradient density is owned iff its complete stencil is coarse;
!> * an on-site term is owned by the active representation of its block.
!>
!> The regular-grid FFT dipole remains an independent all-grid owner.  A
!> validated mask may be rebuilt between complete integration steps by the
!> adaptive CPU orchestration layer.
!-------------------------------------------------------------------------------
module StaticHybridOperator

   use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
   use Parameters, only : dblprec
   use BlockTopology, only : block_topology_type, REGULAR_REPLICATED_CELL
   use stiffness, only : COARSE_MUB_SI
   use CoarseTensorOperator, only : coarse_tensor_operator_type, &
      coarse_energy_terms_type, evaluate_coarse_tensor_operator, COARSE_TENSOR_OK
   use SmoothProjectedOperator, only : smooth_projected_operator_type, &
      prolongate_smooth_directions, restrict_projected_atomistic_field, &
      SMOOTH_PROJECTED_OK

   implicit none
   private

   integer, parameter, public :: STATIC_HYBRID_OK = 0
   integer, parameter, public :: STATIC_HYBRID_INVALID_TOPOLOGY = 1
   integer, parameter, public :: STATIC_HYBRID_INVALID_MASK = 2
   integer, parameter, public :: STATIC_HYBRID_INVALID_INTERACTIONS = 3
   integer, parameter, public :: STATIC_HYBRID_INVALID_STATE = 4
   integer, parameter, public :: STATIC_HYBRID_UNSUPPORTED_MODE = 5
   integer, parameter, public :: STATIC_HYBRID_DEPENDENCY_ERROR = 6

   integer, parameter, public :: STATIC_HYBRID_COARSE = 0
   integer, parameter, public :: STATIC_HYBRID_BUFFER = 1
   integer, parameter, public :: STATIC_HYBRID_FINE = 2

   type, public :: static_hybrid_energy_type
      real(dblprec) :: atomistic_bilinear_j = 0.0_dblprec
      real(dblprec) :: atomistic_onsite_j = 0.0_dblprec
      type(coarse_energy_terms_type) :: coarse
      real(dblprec) :: total_j = 0.0_dblprec
   end type static_hybrid_energy_type

   type, public :: static_hybrid_operator_type
      logical :: ready = .false.
      integer :: n_atoms = 0
      integer :: n_blocks = 0
      integer :: n_bonds = 0
      integer :: safety_dilation_blocks = 0
      integer :: buffer_width_blocks(3) = 0
      real(dblprec) :: maximum_interaction_radius_m = 0.0_dblprec
      integer, allocatable :: block_state(:)
      logical, allocatable :: fine_seed(:)
      logical, allocatable :: atomistic_block(:)
      logical, allocatable :: coarse_block(:)
      logical, allocatable :: atomistic_atom(:)
      integer, allocatable :: bond_atom(:,:)
      logical, allocatable :: active_bond(:)
      logical, allocatable :: atomistic_bond_owner(:)
      type(coarse_tensor_operator_type) :: tensor
      type(smooth_projected_operator_type) :: projection
      ! RCG-06A (F-13/F-17): persistent per-operator scratch, sized once in
      ! setup_static_hybrid_operator, replacing the natoms/nblocks-sized
      ! automatic (stack) arrays previously declared inside
      ! evaluate_static_hybrid_operator. Written and consumed within a single
      ! call only; never passed to a subroutine call alongside `operator`
      ! itself (each is passed as its own actual argument), so no Fortran
      ! argument-aliasing hazard applies here -- contrast
      ! coarse_tensor_operator_type's scratch, which required decoupling
      ! physical_forward_gradient/add_physical_gradient_transpose from
      ! `operator` for exactly that reason.
      real(dblprec), allocatable :: scratch_effective_direction(:,:)
      real(dblprec), allocatable :: scratch_ghost_direction(:,:)
      real(dblprec), allocatable :: scratch_atomistic_field(:,:)
      real(dblprec), allocatable :: scratch_ghost_reaction_field(:,:)
      real(dblprec), allocatable :: scratch_reaction_field(:,:,:)
      real(dblprec), allocatable :: scratch_tensor_field(:,:)
   end type static_hybrid_operator_type

   public :: setup_static_hybrid_operator
   public :: rebuild_static_hybrid_ownership
   public :: evaluate_static_hybrid_operator

contains

   subroutine setup_static_hybrid_operator(operator,topology,tensor,projection, &
         fine_mask,bond_atom,bond_displacement_m,safety_dilation_blocks,status, &
         diagnostic,active_bond)
      type(static_hybrid_operator_type), intent(out) :: operator
      type(block_topology_type), intent(in) :: topology
      type(coarse_tensor_operator_type), intent(in) :: tensor
      type(smooth_projected_operator_type), intent(in) :: projection
      logical, intent(in) :: fine_mask(:)
      integer, intent(in) :: bond_atom(:,:)
      real(dblprec), intent(in) :: bond_displacement_m(:,:)
      integer, intent(in) :: safety_dilation_blocks
      integer, intent(out) :: status
      character(len=*), intent(out) :: diagnostic
      logical, intent(in), optional :: active_bond(:)

      integer :: block, bond
      real(dblprec) :: inverse_block_m1(3,3), fractional_radius

      operator%ready = .false.
      status = STATIC_HYBRID_INVALID_TOPOLOGY
      diagnostic = ''
      if (.not. topology%ready .or. topology%geometry_mode /= REGULAR_REPLICATED_CELL) then
         diagnostic = 'Static hybrid setup requires a ready REGULAR_REPLICATED_CELL topology'
         return
      end if
      if (.not. tensor%ready .or. .not. projection%ready .or. &
          tensor%nblocks /= topology%n_spatial_blocks .or. &
          projection%n_blocks /= topology%n_spatial_blocks .or. &
          projection%n_atoms /= topology%n_atoms .or. projection%n_channels /= 1) then
         status = STATIC_HYBRID_DEPENDENCY_ERROR
         diagnostic = 'Static hybrid setup requires matching ready tensor and one-channel projection operators'
         return
      end if
      if (tensor%use_uniform_coarse_dipole) then
         status = STATIC_HYBRID_UNSUPPORTED_MODE
         diagnostic = 'FFT dipole is an independent all-grid owner and cannot be embedded in short-range hybrid dispatch'
         return
      end if
      if (size(fine_mask) /= topology%n_spatial_blocks) then
         status = STATIC_HYBRID_INVALID_MASK
         diagnostic = 'Static fine mask must have length nblocks'
         return
      end if
      if (safety_dilation_blocks < 0) then
         status = STATIC_HYBRID_INVALID_MASK
         diagnostic = 'Safety dilation must be a nonnegative number of blocks'
         return
      end if
      if (size(bond_atom,1) /= 2 .or. &
          size(bond_displacement_m,1) /= 3 .or. &
          size(bond_displacement_m,2) /= size(bond_atom,2) .or. &
          .not. all(ieee_is_finite(bond_displacement_m))) then
         status = STATIC_HYBRID_INVALID_INTERACTIONS
         diagnostic = 'Interactions require bond_atom(2,nbond) and finite physical displacements(3,nbond)'
         return
      end if
      if (size(bond_atom,2) > 0) then
         if (any(bond_atom < 1) .or. any(bond_atom > topology%n_atoms) .or. &
             any(bond_atom(1,:) == bond_atom(2,:))) then
            status = STATIC_HYBRID_INVALID_INTERACTIONS
            diagnostic = 'Every active short-range bond must join two distinct magnetic atoms'
            return
         end if
      end if
      if (present(active_bond)) then
         if (size(active_bond) /= size(bond_atom,2)) then
            status = STATIC_HYBRID_INVALID_INTERACTIONS
            diagnostic = 'Active-interaction mask must have length nbonds'
            return
         end if
      end if

      operator%n_atoms = topology%n_atoms
      operator%n_blocks = topology%n_spatial_blocks
      operator%n_bonds = size(bond_atom,2)
      operator%safety_dilation_blocks = safety_dilation_blocks
      operator%tensor = tensor
      operator%projection = projection
      allocate(operator%block_state(operator%n_blocks), &
         operator%fine_seed(operator%n_blocks), &
         operator%atomistic_block(operator%n_blocks), &
         operator%coarse_block(operator%n_blocks), &
         operator%atomistic_atom(operator%n_atoms), &
         operator%bond_atom(2,operator%n_bonds), &
         operator%active_bond(operator%n_bonds), &
         operator%atomistic_bond_owner(operator%n_bonds))
      operator%fine_seed = fine_mask
      operator%bond_atom = bond_atom
      operator%active_bond = .true.
      if (present(active_bond)) operator%active_bond = active_bond

      allocate(operator%scratch_effective_direction(3,operator%n_atoms), &
         operator%scratch_ghost_direction(3,operator%n_atoms), &
         operator%scratch_atomistic_field(3,operator%n_atoms), &
         operator%scratch_ghost_reaction_field(3,operator%n_atoms), &
         operator%scratch_reaction_field(3,1,operator%n_blocks), &
         operator%scratch_tensor_field(3,operator%n_blocks))

      operator%maximum_interaction_radius_m = 0.0_dblprec
      do bond = 1, operator%n_bonds
         if (.not. operator%active_bond(bond)) cycle
         operator%maximum_interaction_radius_m = max( &
            operator%maximum_interaction_radius_m, &
            sqrt(sum(bond_displacement_m(:,bond)**2)))
      end do
      inverse_block_m1 = inverse3(topology%block_vectors)
      do block = 1, 3
         fractional_radius = operator%maximum_interaction_radius_m * &
            sqrt(sum(inverse_block_m1(block,:)**2))
         operator%buffer_width_blocks(block) = &
            ceiling(max(0.0_dblprec,fractional_radius-64.0_dblprec*epsilon(1.0_dblprec))) + &
            safety_dilation_blocks
      end do

      operator%ready = .true.
      call rebuild_static_hybrid_ownership(operator,topology,fine_mask,status,diagnostic)
      if (status /= STATIC_HYBRID_OK) operator%ready = .false.
   end subroutine setup_static_hybrid_operator

   !> Rebuild all mask-derived work ownership without changing immutable
   !> interaction/operator data.  The caller applies this only to a temporary
   !> operator and publishes that copy after a transition has passed its energy
   !> check, so a rejected transition cannot corrupt the active dispatch.
   subroutine rebuild_static_hybrid_ownership(operator,topology,fine_mask,status,diagnostic)
      type(static_hybrid_operator_type), intent(inout) :: operator
      type(block_topology_type), intent(in) :: topology
      logical, intent(in) :: fine_mask(:)
      integer, intent(out) :: status
      character(len=*), intent(out) :: diagnostic

      integer :: atom, block, bond, seed
      integer :: delta(3), periodic_delta(3)

      status = STATIC_HYBRID_INVALID_STATE
      diagnostic = ''
      if (.not. operator%ready .or. .not. topology%ready .or. &
          topology%geometry_mode /= REGULAR_REPLICATED_CELL .or. &
          topology%n_atoms /= operator%n_atoms .or. &
          topology%n_spatial_blocks /= operator%n_blocks) then
         diagnostic = 'Hybrid ownership rebuild requires matching ready operator and topology'
         return
      end if
      if (size(fine_mask) /= operator%n_blocks) then
         status = STATIC_HYBRID_INVALID_MASK
         diagnostic = 'Adaptive fine mask must have length nblocks'
         return
      end if
      if (.not. allocated(operator%fine_seed) .or. &
          .not. allocated(operator%atomistic_block) .or. &
          .not. allocated(operator%coarse_block) .or. &
          .not. allocated(operator%block_state) .or. &
          .not. allocated(operator%atomistic_atom) .or. &
          .not. allocated(operator%atomistic_bond_owner)) then
         diagnostic = 'Hybrid ownership arrays are not initialized'
         return
      end if

      operator%fine_seed = fine_mask
      operator%atomistic_block = fine_mask
      do seed = 1, operator%n_blocks
         if (.not. fine_mask(seed)) cycle
         do block = 1, operator%n_blocks
            delta = abs(topology%block_grid_coordinate(:,block) - &
               topology%block_grid_coordinate(:,seed))
            periodic_delta = min(delta,topology%block_grid-delta)
            if (all(periodic_delta <= operator%buffer_width_blocks)) then
               operator%atomistic_block(block) = .true.
            end if
         end do
      end do
      operator%coarse_block = .not. operator%atomistic_block
      operator%block_state = STATIC_HYBRID_COARSE
      where (operator%atomistic_block) operator%block_state = STATIC_HYBRID_BUFFER
      where (fine_mask) operator%block_state = STATIC_HYBRID_FINE

      do atom = 1, operator%n_atoms
         operator%atomistic_atom(atom) = operator%atomistic_block( &
            topology%atom_to_block(atom))
      end do
      do bond = 1, operator%n_bonds
         operator%atomistic_bond_owner(bond) = operator%active_bond(bond) .and. &
            (operator%atomistic_atom(operator%bond_atom(1,bond)) .or. &
             operator%atomistic_atom(operator%bond_atom(2,bond)))
      end do
      status = STATIC_HYBRID_OK
   end subroutine rebuild_static_hybrid_ownership

   subroutine evaluate_static_hybrid_operator(operator,fine_direction,coarse_direction, &
         bond_matrix_j,fine_field_t,coarse_field_t,energies,status,diagnostic, &
         atomistic_onsite_energy_j,atomistic_onsite_field_t,coarse_external_field_t)
      type(static_hybrid_operator_type), intent(inout) :: operator
      real(dblprec), intent(in) :: fine_direction(:,:)
      real(dblprec), intent(in) :: coarse_direction(:,:,:)
      real(dblprec), intent(in) :: bond_matrix_j(:,:,:)
      real(dblprec), intent(out) :: fine_field_t(:,:)
      real(dblprec), intent(out) :: coarse_field_t(:,:,:)
      type(static_hybrid_energy_type), intent(out) :: energies
      integer, intent(out) :: status
      character(len=*), intent(out) :: diagnostic
      real(dblprec), intent(in), optional :: atomistic_onsite_energy_j(:)
      real(dblprec), intent(in), optional :: atomistic_onsite_field_t(:,:)
      real(dblprec), intent(in), optional :: coarse_external_field_t(:,:)

      integer :: atom, atom_i, atom_j, block, bond, dependency_status
      character(len=512) :: dependency_diagnostic

      energies = static_hybrid_energy_type()
      status = STATIC_HYBRID_INVALID_STATE
      diagnostic = ''
      if (.not. operator%ready) then
         diagnostic = 'Cannot evaluate an uninitialized static hybrid operator'
         return
      end if
      if (size(fine_direction,1) /= 3 .or. &
          size(fine_direction,2) /= operator%n_atoms .or. &
          any(shape(coarse_direction) /= (/3,1,operator%n_blocks/)) .or. &
          any(shape(fine_field_t) /= (/3,operator%n_atoms/)) .or. &
          any(shape(coarse_field_t) /= (/3,1,operator%n_blocks/)) .or. &
          .not. all(ieee_is_finite(fine_direction)) .or. &
          .not. all(ieee_is_finite(coarse_direction))) then
         diagnostic = 'Hybrid directions and fields must be finite and match atom/block shapes'
         return
      end if
      if (size(bond_matrix_j,1) /= 3 .or. size(bond_matrix_j,2) /= 3 .or. &
          size(bond_matrix_j,3) /= operator%n_bonds .or. &
          .not. all(ieee_is_finite(bond_matrix_j))) then
         status = STATIC_HYBRID_INVALID_INTERACTIONS
         diagnostic = 'Bilinear matrices must be finite with shape (3,3,nbonds)'
         return
      end if
      if (present(atomistic_onsite_energy_j) .neqv. present(atomistic_onsite_field_t)) then
         diagnostic = 'Atomistic on-site energy and field inputs must be supplied together'
         return
      end if
      if (present(atomistic_onsite_energy_j)) then
         if (size(atomistic_onsite_energy_j) /= operator%n_atoms .or. &
             any(shape(atomistic_onsite_field_t) /= (/3,operator%n_atoms/)) .or. &
             .not. all(ieee_is_finite(atomistic_onsite_energy_j)) .or. &
             .not. all(ieee_is_finite(atomistic_onsite_field_t))) then
            diagnostic = 'Atomistic on-site arrays must be finite and atom-sized'
            return
         end if
      end if
      do atom = 1, operator%n_atoms
         if (operator%atomistic_atom(atom)) then
            if (abs(sum(fine_direction(:,atom)**2)-1.0_dblprec) > 1.0d-8) then
               diagnostic = 'Every active fine/buffer direction must have unit norm'
               return
            end if
         end if
      end do

      associate (effective_direction => operator%scratch_effective_direction, &
            ghost_direction => operator%scratch_ghost_direction, &
            atomistic_field => operator%scratch_atomistic_field, &
            ghost_reaction_field => operator%scratch_ghost_reaction_field, &
            reaction_field => operator%scratch_reaction_field, &
            tensor_field => operator%scratch_tensor_field)
      call prolongate_smooth_directions(operator%projection,coarse_direction, &
         ghost_direction,dependency_status,dependency_diagnostic)
      if (dependency_status /= SMOOTH_PROJECTED_OK) then
         status = STATIC_HYBRID_DEPENDENCY_ERROR
         diagnostic = 'Ghost prolongation failed: '//trim(dependency_diagnostic)
         return
      end if
      effective_direction = ghost_direction
      do atom = 1, operator%n_atoms
         if (operator%atomistic_atom(atom)) effective_direction(:,atom) = fine_direction(:,atom)
      end do

      atomistic_field = 0.0_dblprec
      do bond = 1, operator%n_bonds
         if (.not. operator%atomistic_bond_owner(bond)) cycle
         atom_i = operator%bond_atom(1,bond)
         atom_j = operator%bond_atom(2,bond)
         energies%atomistic_bilinear_j = energies%atomistic_bilinear_j - &
            dot_product(effective_direction(:,atom_i), &
            matmul(bond_matrix_j(:,:,bond),effective_direction(:,atom_j)))
         atomistic_field(:,atom_i) = atomistic_field(:,atom_i) + &
            matmul(bond_matrix_j(:,:,bond),effective_direction(:,atom_j)) / &
            (COARSE_MUB_SI*operator%projection%atom_moment_mub(atom_i))
         atomistic_field(:,atom_j) = atomistic_field(:,atom_j) + &
            matmul(transpose(bond_matrix_j(:,:,bond)),effective_direction(:,atom_i)) / &
            (COARSE_MUB_SI*operator%projection%atom_moment_mub(atom_j))
      end do

      if (present(atomistic_onsite_energy_j)) then
         do atom = 1, operator%n_atoms
            if (.not. operator%atomistic_atom(atom)) cycle
            energies%atomistic_onsite_j = energies%atomistic_onsite_j + &
               atomistic_onsite_energy_j(atom)
            atomistic_field(:,atom) = atomistic_field(:,atom) + &
               atomistic_onsite_field_t(:,atom)
         end do
      end if

      fine_field_t = 0.0_dblprec
      ghost_reaction_field = 0.0_dblprec
      do atom = 1, operator%n_atoms
         if (operator%atomistic_atom(atom)) then
            fine_field_t(:,atom) = atomistic_field(:,atom)
         else
            ghost_reaction_field(:,atom) = atomistic_field(:,atom)
         end if
      end do
      call restrict_projected_atomistic_field(operator%projection,coarse_direction, &
         ghost_reaction_field,reaction_field,dependency_status,dependency_diagnostic)
      if (dependency_status /= SMOOTH_PROJECTED_OK) then
         status = STATIC_HYBRID_DEPENDENCY_ERROR
         diagnostic = 'Ghost reaction restriction failed: '//trim(dependency_diagnostic)
         return
      end if

      if (present(coarse_external_field_t)) then
         call evaluate_coarse_tensor_operator(operator%tensor,coarse_direction(:,1,:), &
            tensor_field,energies%coarse,dependency_status,dependency_diagnostic, &
            external_field_t=coarse_external_field_t, &
            interaction_owner=operator%coarse_block,onsite_owner=operator%coarse_block)
      else
         call evaluate_coarse_tensor_operator(operator%tensor,coarse_direction(:,1,:), &
            tensor_field,energies%coarse,dependency_status,dependency_diagnostic, &
            interaction_owner=operator%coarse_block,onsite_owner=operator%coarse_block)
      end if
      if (dependency_status /= COARSE_TENSOR_OK) then
         status = STATIC_HYBRID_DEPENDENCY_ERROR
         diagnostic = 'Masked tensor evaluation failed: '//trim(dependency_diagnostic)
         return
      end if
      coarse_field_t = 0.0_dblprec
      do block = 1, operator%n_blocks
         if (operator%coarse_block(block)) then
            coarse_field_t(:,1,block) = tensor_field(:,block) + reaction_field(:,1,block)
         end if
      end do
      energies%total_j = energies%atomistic_bilinear_j + &
         energies%atomistic_onsite_j + energies%coarse%total_j
      status = STATIC_HYBRID_OK
      end associate
   end subroutine evaluate_static_hybrid_operator

   pure function inverse3(matrix) result(inverse)
      real(dblprec), intent(in) :: matrix(3,3)
      real(dblprec) :: inverse(3,3), determinant

      determinant = matrix(1,1)*(matrix(2,2)*matrix(3,3)-matrix(2,3)*matrix(3,2)) - &
         matrix(1,2)*(matrix(2,1)*matrix(3,3)-matrix(2,3)*matrix(3,1)) + &
         matrix(1,3)*(matrix(2,1)*matrix(3,2)-matrix(2,2)*matrix(3,1))
      inverse(1,1) = (matrix(2,2)*matrix(3,3)-matrix(2,3)*matrix(3,2))/determinant
      inverse(1,2) = (matrix(1,3)*matrix(3,2)-matrix(1,2)*matrix(3,3))/determinant
      inverse(1,3) = (matrix(1,2)*matrix(2,3)-matrix(1,3)*matrix(2,2))/determinant
      inverse(2,1) = (matrix(2,3)*matrix(3,1)-matrix(2,1)*matrix(3,3))/determinant
      inverse(2,2) = (matrix(1,1)*matrix(3,3)-matrix(1,3)*matrix(3,1))/determinant
      inverse(2,3) = (matrix(1,3)*matrix(2,1)-matrix(1,1)*matrix(2,3))/determinant
      inverse(3,1) = (matrix(2,1)*matrix(3,2)-matrix(2,2)*matrix(3,1))/determinant
      inverse(3,2) = (matrix(1,2)*matrix(3,1)-matrix(1,1)*matrix(3,2))/determinant
      inverse(3,3) = (matrix(1,1)*matrix(2,2)-matrix(1,2)*matrix(2,1))/determinant
   end function inverse3

end module StaticHybridOperator
