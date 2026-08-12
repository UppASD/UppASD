!-------------------------------------------------------------------------------
! MODULE: CoarseTensorOperator
!> @brief Double-precision CPU reference operator for an all-coarse FM grid.
!>
!> This module owns no atomistic state.  It consumes the immutable canonical
!> topology and the accepted SI coarse-material descriptor, then evaluates the
!> CG-01 discrete energy with transparent periodic reference loops.  The
!> uniformly coarse dipole solver remains the owner of its field; that field is
!> accepted here so that it participates in the common field and energy report.
!-------------------------------------------------------------------------------
module CoarseTensorOperator

   use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
   use Parameters, only : dblprec
   use BlockTopology, only : block_topology_type, REGULAR_REPLICATED_CELL, &
      regular_spatial_block_id
   use stiffness, only : coarse_material_type, coarse_material_runtime_status, &
      COARSE_RUNTIME_SINGLE_FM, COARSE_MATERIAL_OK, COARSE_DMI_ENERGY_PLUS, &
      COARSE_FIELD_DERIVATIVE_MINUS, COARSE_MUB_SI

   implicit none
   private

   integer, parameter, public :: COARSE_TENSOR_OK = 0
   integer, parameter, public :: COARSE_TENSOR_INVALID_TOPOLOGY = 1
   integer, parameter, public :: COARSE_TENSOR_INVALID_MATERIAL = 2
   integer, parameter, public :: COARSE_TENSOR_UNSUPPORTED_GEOMETRY = 3
   integer, parameter, public :: COARSE_TENSOR_UNSUPPORTED_SOLVER = 4
   integer, parameter, public :: COARSE_TENSOR_UNSUPPORTED_TEMPERATURE = 5
   integer, parameter, public :: COARSE_TENSOR_UNSUPPORTED_MODE = 6
   integer, parameter, public :: COARSE_TENSOR_INVALID_ANISOTROPY = 7
   integer, parameter, public :: COARSE_TENSOR_INVALID_STATE = 8

   integer, parameter, public :: COARSE_BOUNDARY_PERIODIC = 1
   integer, parameter, public :: COARSE_SOLVER_DETERMINISTIC_HEUN = 1
   integer, parameter, public :: COARSE_ANISOTROPY_NONE = 0
   integer, parameter, public :: COARSE_ANISOTROPY_UNIAXIAL = 1

   type, public :: coarse_operator_options_type
      integer :: boundary_mode(3) = COARSE_BOUNDARY_PERIODIC
      integer :: solver_mode = COARSE_SOLVER_DETERMINISTIC_HEUN
      real(dblprec) :: temperature_k = 0.0_dblprec
      logical :: stochastic_field = .false.
      logical :: adaptive_switching = .false.
      logical :: interface_coupling = .false.
      logical :: restart_state = .false.
      logical :: time_dependent_external_field = .false.
      logical :: use_uniform_coarse_dipole = .false.
   end type coarse_operator_options_type

   !> Per-block one- or two-axis UppASD-compatible anisotropy, in J/m^3.
   type, public :: coarse_anisotropy_type
      integer, allocatable :: kind(:)
      integer, allocatable :: axis_count(:)
      real(dblprec), allocatable :: axis(:,:,:) ! (spin, axis, block)
      real(dblprec), allocatable :: k1(:,:)      ! (axis, block)
      real(dblprec), allocatable :: k2(:,:)      ! (axis, block)
   end type coarse_anisotropy_type

   type, public :: coarse_energy_terms_type
      real(dblprec) :: exchange_j = 0.0_dblprec
      real(dblprec) :: spiralization_j = 0.0_dblprec
      real(dblprec) :: anisotropy_j = 0.0_dblprec
      real(dblprec) :: external_j = 0.0_dblprec
      real(dblprec) :: dipole_j = 0.0_dblprec
      real(dblprec) :: total_j = 0.0_dblprec
   end type coarse_energy_terms_type

   type, public :: coarse_field_terms_type
      real(dblprec), allocatable :: exchange_t(:,:)
      real(dblprec), allocatable :: spiralization_t(:,:)
      real(dblprec), allocatable :: anisotropy_t(:,:)
      real(dblprec), allocatable :: external_t(:,:)
      real(dblprec), allocatable :: dipole_t(:,:)
   end type coarse_field_terms_type

   type, public :: coarse_tensor_operator_type
      logical :: ready = .false.
      integer :: nblocks = 0
      integer :: block_grid(3) = 0
      integer, allocatable :: block_coordinate(:,:)
      real(dblprec) :: block_vectors_m(3,3) = 0.0_dblprec
      real(dblprec) :: inverse_block_transpose_m1(3,3) = 0.0_dblprec
      real(dblprec) :: block_volume_m3 = 0.0_dblprec
      real(dblprec) :: block_moment_mub = 0.0_dblprec
      real(dblprec) :: exchange_stiffness_j_per_m(3,3) = 0.0_dblprec
      real(dblprec) :: spiralization_j_per_m2(3,3) = 0.0_dblprec
      real(dblprec) :: channel_gamma_per_t_s = 0.0_dblprec
      real(dblprec) :: channel_damping = 0.0_dblprec
      logical :: use_uniform_coarse_dipole = .false.
      type(coarse_anisotropy_type) :: anisotropy
      ! RCG-06A (F-17): persistent per-operator scratch, sized once in
      ! setup_coarse_tensor_operator and reused by every
      ! evaluate_coarse_tensor_operator call instead of the ten nblocks-sized
      ! temporaries previously heap-allocated on every call (this operator is
      ! evaluated twice per Heun stage per ensemble, every adaptive step).
      ! Written and consumed within a single call only; not durable state
      ! between calls, and (like every other RCG-06A operator scratch field)
      ! not safe to share across a future OpenMP-parallel evaluation of the
      ! SAME operator instance from different ensembles -- see RCG-07.
      real(dblprec), allocatable :: scratch_gradient(:,:,:)
      real(dblprec), allocatable :: scratch_exchange_derivative(:,:)
      real(dblprec), allocatable :: scratch_dmi_derivative(:,:)
      real(dblprec), allocatable :: scratch_anisotropy_derivative(:,:)
      real(dblprec), allocatable :: scratch_work(:,:)
      real(dblprec), allocatable :: scratch_exchange_field(:,:)
      real(dblprec), allocatable :: scratch_dmi_field(:,:)
      real(dblprec), allocatable :: scratch_anisotropy_field(:,:)
      real(dblprec), allocatable :: scratch_external_field(:,:)
      real(dblprec), allocatable :: scratch_dipole_field(:,:)
   end type coarse_tensor_operator_type

   public :: setup_coarse_tensor_operator
   public :: evaluate_coarse_tensor_operator
   public :: coarse_llg_rhs

contains

   !> Validate the complete initial capability matrix before retaining state.
   subroutine setup_coarse_tensor_operator(operator, topology, material, options, &
         status, diagnostic, anisotropy)
      type(coarse_tensor_operator_type), intent(out) :: operator
      type(block_topology_type), intent(in) :: topology
      type(coarse_material_type), intent(in) :: material
      type(coarse_operator_options_type), intent(in) :: options
      integer, intent(out) :: status
      character(len=*), intent(out) :: diagnostic
      type(coarse_anisotropy_type), intent(in), optional :: anisotropy

      integer :: runtime_status, cells_per_block, block, axis_index
      real(dblprec) :: cell_volume, determinant, axis_norm, scale
      character(len=512) :: runtime_diagnostic

      operator%ready = .false.
      status = COARSE_TENSOR_INVALID_TOPOLOGY
      diagnostic = ''

      if (.not. topology%ready) then
         diagnostic = 'All-coarse tensor setup requires a ready canonical topology'
         return
      end if
      if (topology%geometry_mode /= REGULAR_REPLICATED_CELL) then
         status = COARSE_TENSOR_UNSUPPORTED_GEOMETRY
         diagnostic = 'All-coarse tensor setup supports REGULAR_REPLICATED_CELL geometry only'
         return
      end if
      if (topology%n_dynamic_channels /= 1 .or. &
          size(topology%block_dynamic_channel_population,1) /= 1) then
         status = COARSE_TENSOR_UNSUPPORTED_MODE
         diagnostic = 'All-coarse tensor setup supports exactly one ferromagnetic dynamical channel'
         return
      end if
      if (any(topology%block_dynamic_channel_population(1,:) <= 0)) then
         diagnostic = 'Every all-coarse block must contain the ferromagnetic channel'
         return
      end if
      if (any(options%boundary_mode /= COARSE_BOUNDARY_PERIODIC)) then
         status = COARSE_TENSOR_UNSUPPORTED_GEOMETRY
         diagnostic = 'All-coarse tensor setup supports periodic block derivatives only'
         return
      end if
      if (options%solver_mode /= COARSE_SOLVER_DETERMINISTIC_HEUN) then
         status = COARSE_TENSOR_UNSUPPORTED_SOLVER
         diagnostic = 'All-coarse tensor setup supports the deterministic Heun LLG algebra only'
         return
      end if
      if (.not. ieee_is_finite(options%temperature_k) .or. &
          options%temperature_k /= 0.0_dblprec .or. options%stochastic_field) then
         status = COARSE_TENSOR_UNSUPPORTED_TEMPERATURE
         diagnostic = 'All-coarse tensor setup requires T=0 and no stochastic field'
         return
      end if
      if (options%adaptive_switching .or. options%interface_coupling .or. &
          options%restart_state .or. options%time_dependent_external_field) then
         status = COARSE_TENSOR_UNSUPPORTED_MODE
         diagnostic = 'CG-04 is static all-coarse only: adaptive, interface, restart, and time-dependent modes are rejected'
         return
      end if

      call coarse_material_runtime_status(material, COARSE_RUNTIME_SINGLE_FM, &
         runtime_status, runtime_diagnostic)
      if (runtime_status /= COARSE_MATERIAL_OK) then
         status = COARSE_TENSOR_INVALID_MATERIAL
         diagnostic = 'All-coarse tensor material rejected: '//trim(runtime_diagnostic)
         return
      end if
      if (material%metadata%dmi_energy_sign /= COARSE_DMI_ENERGY_PLUS .or. &
          material%metadata%field_derivative_sign /= COARSE_FIELD_DERIVATIVE_MINUS .or. &
          trim(material%metadata%convention_version) /= 'CG-01-approved-v1') then
         status = COARSE_TENSOR_INVALID_MATERIAL
         diagnostic = 'All-coarse tensor material does not use the approved CG-01 sign and derivative convention'
         return
      end if
      if (.not. allocated(material%channel_moment_mub) .or. &
          .not. allocated(material%channel_gamma) .or. &
          .not. allocated(material%channel_damping) .or. &
          .not. allocated(material%basis_to_channel) .or. &
          .not. allocated(material%basis_moment_mub) .or. &
          .not. allocated(material%exchange_stiffness) .or. &
          .not. allocated(material%spiralization)) then
         status = COARSE_TENSOR_INVALID_MATERIAL
         diagnostic = 'All-coarse tensor material is missing mandatory channel or tensor storage'
         return
      end if
      if (size(material%channel_moment_mub) /= 1 .or. &
          size(material%channel_gamma) /= 1 .or. &
          size(material%channel_damping) /= 1 .or. &
          any(shape(material%exchange_stiffness) /= (/3,3,1,1/)) .or. &
          any(shape(material%spiralization) /= (/3,3,1,1/))) then
         status = COARSE_TENSOR_INVALID_MATERIAL
         diagnostic = 'All-coarse tensor material arrays must retain one explicit channel dimension'
         return
      end if
      if (material%n_basis /= topology%n_basis .or. &
          size(material%basis_to_channel) /= topology%n_basis .or. &
          size(material%basis_moment_mub) /= topology%n_basis .or. &
          any(material%basis_to_channel /= topology%basis_to_dynamic_channel)) then
         status = COARSE_TENSOR_INVALID_MATERIAL
         diagnostic = 'Material basis/channel map must match the canonical topology exactly'
         return
      end if
      if (.not. all(ieee_is_finite(material%exchange_stiffness(:,:,1,1))) .or. &
          .not. all(ieee_is_finite(material%spiralization(:,:,1,1))) .or. &
          maxval(abs(material%exchange_stiffness(:,:,1,1) - &
                     transpose(material%exchange_stiffness(:,:,1,1)))) > &
          1.0d-12 * max(tiny(1.0_dblprec), &
                         maxval(abs(material%exchange_stiffness(:,:,1,1))))) then
         status = COARSE_TENSOR_INVALID_MATERIAL
         diagnostic = 'Single-channel exchange stiffness must be finite and Cartesian-symmetric'
         return
      end if
      if (.not. ieee_is_finite(material%cell_volume_m3) .or. &
          .not. all(ieee_is_finite(material%basis_moment_mub)) .or. &
          .not. all(ieee_is_finite(material%channel_moment_mub)) .or. &
          .not. all(ieee_is_finite(material%channel_gamma)) .or. &
          .not. all(ieee_is_finite(material%channel_damping)) .or. &
          material%cell_volume_m3 <= 0.0_dblprec .or. &
          material%channel_moment_mub(1) <= 0.0_dblprec .or. &
          material%channel_gamma(1) <= 0.0_dblprec .or. &
          material%channel_damping(1) < 0.0_dblprec) then
         status = COARSE_TENSOR_INVALID_MATERIAL
         diagnostic = 'All-coarse tensor dynamics require positive moment/gamma and nonnegative damping'
         return
      end if

      determinant = determinant3(topology%block_vectors)
      cell_volume = abs(determinant3(topology%cell_vectors))
      if (.not. ieee_is_finite(determinant) .or. abs(determinant) <= tiny(determinant)) then
         status = COARSE_TENSOR_UNSUPPORTED_GEOMETRY
         diagnostic = 'Physical block vectors must span a finite nonzero volume'
         return
      end if
      scale = max(cell_volume, material%cell_volume_m3)
      if (abs(cell_volume-material%cell_volume_m3) > 1.0d-10*scale) then
         status = COARSE_TENSOR_INVALID_MATERIAL
         diagnostic = 'Material SI cell volume does not match the canonical physical cell metric'
         return
      end if
      if (maxval(abs(topology%block_volume-abs(determinant))) > &
          1.0d-12*max(abs(determinant),tiny(determinant))) then
         status = COARSE_TENSOR_UNSUPPORTED_GEOMETRY
         diagnostic = 'CG-04 requires a uniform physical block volume'
         return
      end if

      operator%nblocks = topology%n_spatial_blocks
      operator%block_grid = topology%block_grid
      operator%block_vectors_m = topology%block_vectors
      operator%block_volume_m3 = abs(determinant)
      operator%inverse_block_transpose_m1 = transpose(inverse3(topology%block_vectors))
      cells_per_block = product(topology%block_shape)
      if (any(topology%block_dynamic_channel_population(1,:) /= &
          count(material%basis_to_channel == 1)*cells_per_block)) then
         status = COARSE_TENSOR_INVALID_MATERIAL
         diagnostic = 'Material channel population does not match every canonical block'
         return
      end if
      operator%block_moment_mub = material%channel_moment_mub(1) * real(cells_per_block,dblprec)
      operator%exchange_stiffness_j_per_m = material%exchange_stiffness(:,:,1,1)
      operator%spiralization_j_per_m2 = material%spiralization(:,:,1,1)
      operator%channel_gamma_per_t_s = material%channel_gamma(1)
      operator%channel_damping = material%channel_damping(1)
      operator%use_uniform_coarse_dipole = options%use_uniform_coarse_dipole
      allocate(operator%block_coordinate(3,operator%nblocks))
      operator%block_coordinate = topology%block_grid_coordinate

      allocate(operator%anisotropy%kind(operator%nblocks), &
         operator%anisotropy%axis_count(operator%nblocks), &
         operator%anisotropy%axis(3,2,operator%nblocks), &
         operator%anisotropy%k1(2,operator%nblocks), &
         operator%anisotropy%k2(2,operator%nblocks))
      operator%anisotropy%kind = COARSE_ANISOTROPY_NONE
      operator%anisotropy%axis_count = 0
      operator%anisotropy%axis = 0.0_dblprec
      operator%anisotropy%k1 = 0.0_dblprec
      operator%anisotropy%k2 = 0.0_dblprec

      if (present(anisotropy)) then
         if (.not. allocated(anisotropy%kind) .or. &
             .not. allocated(anisotropy%axis_count) .or. &
             .not. allocated(anisotropy%axis) .or. &
             .not. allocated(anisotropy%k1) .or. .not. allocated(anisotropy%k2)) then
            status = COARSE_TENSOR_INVALID_ANISOTROPY
            diagnostic = 'Anisotropy descriptor must allocate kind, count, axes, K1, and K2 together'
            return
         end if
         if (size(anisotropy%kind) /= operator%nblocks .or. &
             size(anisotropy%axis_count) /= operator%nblocks .or. &
             any(shape(anisotropy%axis) /= (/3,2,operator%nblocks/)) .or. &
             any(shape(anisotropy%k1) /= (/2,operator%nblocks/)) .or. &
             any(shape(anisotropy%k2) /= (/2,operator%nblocks/))) then
            status = COARSE_TENSOR_INVALID_ANISOTROPY
            diagnostic = 'Anisotropy descriptor dimensions must match (spin,axis,block)=(3,2,nblocks)'
            return
         end if
         if (.not. all(ieee_is_finite(anisotropy%axis)) .or. &
             .not. all(ieee_is_finite(anisotropy%k1)) .or. &
             .not. all(ieee_is_finite(anisotropy%k2))) then
            status = COARSE_TENSOR_INVALID_ANISOTROPY
            diagnostic = 'Anisotropy axes and coefficients must be finite SI values'
            return
         end if
         do block = 1, operator%nblocks
            if (anisotropy%kind(block) /= COARSE_ANISOTROPY_NONE .and. &
                anisotropy%kind(block) /= COARSE_ANISOTROPY_UNIAXIAL) then
               status = COARSE_TENSOR_INVALID_ANISOTROPY
               diagnostic = 'CG-04 supports only none or UppASD-compatible uniaxial anisotropy'
               return
            end if
            if (anisotropy%axis_count(block) < 0 .or. anisotropy%axis_count(block) > 2 .or. &
                (anisotropy%kind(block) == COARSE_ANISOTROPY_NONE .and. &
                 anisotropy%axis_count(block) /= 0) .or. &
                (anisotropy%kind(block) == COARSE_ANISOTROPY_UNIAXIAL .and. &
                 anisotropy%axis_count(block) < 1)) then
               status = COARSE_TENSOR_INVALID_ANISOTROPY
               diagnostic = 'Uniaxial anisotropy requires one or two axes; none requires zero axes'
               return
            end if
            do axis_index = 1, anisotropy%axis_count(block)
               axis_norm = sqrt(sum(anisotropy%axis(:,axis_index,block)**2))
               if (abs(axis_norm-1.0_dblprec) > 1.0d-12) then
                  status = COARSE_TENSOR_INVALID_ANISOTROPY
                  diagnostic = 'Every supported anisotropy axis must be a normalized Cartesian vector'
                  return
               end if
            end do
         end do
         operator%anisotropy%kind = anisotropy%kind
         operator%anisotropy%axis_count = anisotropy%axis_count
         operator%anisotropy%axis = anisotropy%axis
         operator%anisotropy%k1 = anisotropy%k1
         operator%anisotropy%k2 = anisotropy%k2
      end if

      allocate(operator%scratch_gradient(3,3,operator%nblocks), &
         operator%scratch_exchange_derivative(3,operator%nblocks), &
         operator%scratch_dmi_derivative(3,operator%nblocks), &
         operator%scratch_anisotropy_derivative(3,operator%nblocks), &
         operator%scratch_work(3,operator%nblocks), &
         operator%scratch_exchange_field(3,operator%nblocks), &
         operator%scratch_dmi_field(3,operator%nblocks), &
         operator%scratch_anisotropy_field(3,operator%nblocks), &
         operator%scratch_external_field(3,operator%nblocks), &
         operator%scratch_dipole_field(3,operator%nblocks))

      operator%ready = .true.
      status = COARSE_TENSOR_OK
      diagnostic = ''
   end subroutine setup_coarse_tensor_operator

   !> Evaluate explicit per-term energies and their unconstrained SI fields.
   !>
   !> Optional ownership masks are used by the static hybrid reference.  A
   !> gradient-density term is retained only when its complete forward
   !> difference stencil belongs to interaction_owner.  Local anisotropy and
   !> external-field terms use onsite_owner.  The uniformly coarse FFT dipole
   !> remains a separate all-grid owner and is deliberately not masked here.
   subroutine evaluate_coarse_tensor_operator(operator, direction, field_t, energies, &
         status, diagnostic, external_field_t, uniform_coarse_dipole_field_t, term_fields, &
         interaction_owner, onsite_owner)
      type(coarse_tensor_operator_type), intent(inout) :: operator
      real(dblprec), intent(in) :: direction(:,:)
      real(dblprec), intent(out) :: field_t(:,:)
      type(coarse_energy_terms_type), intent(out) :: energies
      integer, intent(out) :: status
      character(len=*), intent(out) :: diagnostic
      real(dblprec), intent(in), optional :: external_field_t(:,:)
      real(dblprec), intent(in), optional :: uniform_coarse_dipole_field_t(:,:)
      type(coarse_field_terms_type), intent(out), optional :: term_fields
      logical, intent(in), optional :: interaction_owner(:)
      logical, intent(in), optional :: onsite_owner(:)

      integer :: block, p, q, k, axis_index
      real(dblprec) :: c, density, conversion
      real(dblprec) :: basis(3), local_derivative(3)
      logical :: owned

      energies = coarse_energy_terms_type()
      status = COARSE_TENSOR_INVALID_STATE
      diagnostic = ''
      if (.not. operator%ready) then
         diagnostic = 'Cannot evaluate an uninitialized all-coarse tensor operator'
         return
      end if
      if (size(direction,1) /= 3 .or. size(direction,2) /= operator%nblocks .or. &
          size(field_t,1) /= 3 .or. size(field_t,2) /= operator%nblocks) then
         diagnostic = 'Coarse directions and fields must have shape (3,nblocks)'
         return
      end if
      if (.not. all(ieee_is_finite(direction))) then
         diagnostic = 'Coarse directions must be finite'
         return
      end if
      do block = 1, operator%nblocks
         if (abs(sum(direction(:,block)**2)-1.0_dblprec) > 1.0d-8) then
            diagnostic = 'Every coarse direction must have unit norm'
            return
         end if
      end do
      if (present(external_field_t)) then
         if (size(external_field_t,1) /= 3 .or. &
             size(external_field_t,2) /= operator%nblocks .or. &
             .not. all(ieee_is_finite(external_field_t))) then
            diagnostic = 'Static external field must be finite with shape (3,nblocks)'
            return
         end if
      end if
      if (present(interaction_owner)) then
         if (size(interaction_owner) /= operator%nblocks) then
            diagnostic = 'Interaction ownership mask must have length nblocks'
            return
         end if
      end if
      if (present(onsite_owner)) then
         if (size(onsite_owner) /= operator%nblocks) then
            diagnostic = 'On-site ownership mask must have length nblocks'
            return
         end if
      end if
      if (operator%use_uniform_coarse_dipole .neqv. present(uniform_coarse_dipole_field_t)) then
         diagnostic = 'Uniform coarse dipole field presence must match the setup capability'
         return
      end if
      if (present(uniform_coarse_dipole_field_t)) then
         if (size(uniform_coarse_dipole_field_t,1) /= 3 .or. &
             size(uniform_coarse_dipole_field_t,2) /= operator%nblocks .or. &
             .not. all(ieee_is_finite(uniform_coarse_dipole_field_t))) then
            diagnostic = 'Uniform coarse dipole field must be finite with shape (3,nblocks)'
            return
         end if
      end if

      associate (gradient => operator%scratch_gradient, &
            exchange_derivative => operator%scratch_exchange_derivative, &
            dmi_derivative => operator%scratch_dmi_derivative, &
            anisotropy_derivative => operator%scratch_anisotropy_derivative, &
            work => operator%scratch_work, &
            exchange_field => operator%scratch_exchange_field, &
            dmi_field => operator%scratch_dmi_field, &
            anisotropy_field => operator%scratch_anisotropy_field, &
            external_field => operator%scratch_external_field, &
            dipole_field => operator%scratch_dipole_field)
      exchange_derivative = 0.0_dblprec
      dmi_derivative = 0.0_dblprec
      anisotropy_derivative = 0.0_dblprec
      external_field = 0.0_dblprec
      dipole_field = 0.0_dblprec

      call physical_forward_gradient(operator%nblocks,operator%block_grid, &
         operator%block_coordinate,operator%inverse_block_transpose_m1,direction,gradient)

      do p = 1, 3
         do q = 1, 3
            work = 0.0_dblprec
            do block = 1, operator%nblocks
               owned = gradient_term_owned(operator,interaction_owner,block,p,q)
               if (.not. owned) cycle
               energies%exchange_j = energies%exchange_j + operator%block_volume_m3 * &
                  operator%exchange_stiffness_j_per_m(p,q) * &
                  dot_product(gradient(:,p,block),gradient(:,q,block))
               work(:,block) = operator%block_volume_m3 * &
                  operator%exchange_stiffness_j_per_m(p,q) * gradient(:,q,block)
            end do
            call add_physical_gradient_transpose(operator%nblocks,operator%block_grid, &
               operator%block_coordinate,operator%inverse_block_transpose_m1,p,work, &
               exchange_derivative,2.0_dblprec)
         end do
      end do

      do k = 1, 3
         basis = 0.0_dblprec
         basis(k) = 1.0_dblprec
         do p = 1, 3
            work = 0.0_dblprec
            do block = 1, operator%nblocks
               owned = gradient_term_owned(operator,interaction_owner,block,p)
               if (.not. owned) cycle
               energies%spiralization_j = energies%spiralization_j + &
                  operator%block_volume_m3 * operator%spiralization_j_per_m2(k,p) * &
                  dot_product(basis,cross3(direction(:,block),gradient(:,p,block)))
               dmi_derivative(:,block) = dmi_derivative(:,block) - &
                  operator%block_volume_m3 * operator%spiralization_j_per_m2(k,p) * &
                  cross3(basis,gradient(:,p,block))
               work(:,block) = operator%block_volume_m3 * &
                  operator%spiralization_j_per_m2(k,p) * cross3(basis,direction(:,block))
            end do
            call add_physical_gradient_transpose(operator%nblocks,operator%block_grid, &
               operator%block_coordinate,operator%inverse_block_transpose_m1,p,work, &
               dmi_derivative,1.0_dblprec)
         end do
      end do

      do block = 1, operator%nblocks
         if (present(onsite_owner)) then
            if (.not. onsite_owner(block)) cycle
         end if
         local_derivative = 0.0_dblprec
         do axis_index = 1, operator%anisotropy%axis_count(block)
            c = dot_product(direction(:,block),operator%anisotropy%axis(:,axis_index,block))
            density = operator%anisotropy%k1(axis_index,block)*c*c + &
               2.0_dblprec*operator%anisotropy%k2(axis_index,block)*c*c - &
               operator%anisotropy%k2(axis_index,block)*c**4
            energies%anisotropy_j = energies%anisotropy_j + &
               operator%block_volume_m3*density
            local_derivative = local_derivative + 2.0_dblprec*c * &
               (operator%anisotropy%k1(axis_index,block) + &
                2.0_dblprec*operator%anisotropy%k2(axis_index,block)*(1.0_dblprec-c*c)) * &
               operator%anisotropy%axis(:,axis_index,block)
         end do
         anisotropy_derivative(:,block) = operator%block_volume_m3*local_derivative
      end do

      if (present(external_field_t)) then
         do block = 1, operator%nblocks
            if (present(onsite_owner)) then
               if (.not. onsite_owner(block)) cycle
            end if
            external_field(:,block) = external_field_t(:,block)
            energies%external_j = energies%external_j - COARSE_MUB_SI * &
               operator%block_moment_mub * dot_product(external_field(:,block),direction(:,block))
         end do
      end if
      if (present(uniform_coarse_dipole_field_t)) then
         dipole_field = uniform_coarse_dipole_field_t
         do block = 1, operator%nblocks
            energies%dipole_j = energies%dipole_j - 0.5_dblprec*COARSE_MUB_SI * &
               operator%block_moment_mub * dot_product(dipole_field(:,block),direction(:,block))
         end do
      end if

      conversion = -1.0_dblprec/(COARSE_MUB_SI*operator%block_moment_mub)
      exchange_field = conversion*exchange_derivative
      dmi_field = conversion*dmi_derivative
      anisotropy_field = conversion*anisotropy_derivative
      field_t = exchange_field + dmi_field + anisotropy_field + external_field + dipole_field
      energies%total_j = energies%exchange_j + energies%spiralization_j + &
         energies%anisotropy_j + energies%external_j + energies%dipole_j

      if (present(term_fields)) then
         allocate(term_fields%exchange_t(3,operator%nblocks), &
            term_fields%spiralization_t(3,operator%nblocks), &
            term_fields%anisotropy_t(3,operator%nblocks), &
            term_fields%external_t(3,operator%nblocks), &
            term_fields%dipole_t(3,operator%nblocks))
         term_fields%exchange_t = exchange_field
         term_fields%spiralization_t = dmi_field
         term_fields%anisotropy_t = anisotropy_field
         term_fields%external_t = external_field
         term_fields%dipole_t = dipole_field
      end if
      end associate

      status = COARSE_TENSOR_OK
      diagnostic = ''
   end subroutine evaluate_coarse_tensor_operator

   logical function gradient_term_owned(operator,owner,block,p,q) result(owned)
      type(coarse_tensor_operator_type), intent(in) :: operator
      logical, intent(in), optional :: owner(:)
      integer, intent(in) :: block, p
      integer, intent(in), optional :: q

      integer :: direction_index, plus_block
      integer :: coordinate(3)

      owned = .true.
      if (.not. present(owner)) return
      if (.not. owner(block)) then
         owned = .false.
         return
      end if
      do direction_index = 1, 3
         if (operator%inverse_block_transpose_m1(p,direction_index) == 0.0_dblprec) then
            if (.not. present(q)) cycle
            if (operator%inverse_block_transpose_m1(q,direction_index) == 0.0_dblprec) cycle
         end if
         coordinate = operator%block_coordinate(:,block)
         coordinate(direction_index) = modulo(coordinate(direction_index)+1, &
            operator%block_grid(direction_index))
         plus_block = regular_spatial_block_id(coordinate,operator%block_grid)
         if (.not. owner(plus_block)) then
            owned = .false.
            return
         end if
      end do
   end function gradient_term_owned

   !> Deterministic Gilbert LLG right-hand side, matching the existing algebra.
   subroutine coarse_llg_rhs(operator, direction, field_t, derivative_per_s, status, diagnostic)
      type(coarse_tensor_operator_type), intent(in) :: operator
      real(dblprec), intent(in) :: direction(:,:), field_t(:,:)
      real(dblprec), intent(out) :: derivative_per_s(:,:)
      integer, intent(out) :: status
      character(len=*), intent(out) :: diagnostic

      integer :: block
      real(dblprec) :: prefactor

      status = COARSE_TENSOR_INVALID_STATE
      diagnostic = ''
      if (.not. operator%ready .or. size(direction,1) /= 3 .or. &
          size(direction,2) /= operator%nblocks .or. &
          any(shape(field_t) /= shape(direction)) .or. &
          any(shape(derivative_per_s) /= shape(direction))) then
         diagnostic = 'Coarse LLG arrays require a ready operator and shape (3,nblocks)'
         return
      end if
      if (.not. all(ieee_is_finite(direction)) .or. .not. all(ieee_is_finite(field_t))) then
         diagnostic = 'Coarse LLG direction and field must be finite'
         return
      end if
      prefactor = -operator%channel_gamma_per_t_s / &
         (1.0_dblprec+operator%channel_damping**2)
      do block = 1, operator%nblocks
         derivative_per_s(:,block) = prefactor * &
            (cross3(direction(:,block),field_t(:,block)) + operator%channel_damping * &
             cross3(direction(:,block),cross3(direction(:,block),field_t(:,block))))
      end do
      status = COARSE_TENSOR_OK
   end subroutine coarse_llg_rhs

   !> RCG-06A: takes explicit geometry rather than the whole operator so that
   !> a caller may pass one of operator's own scratch arrays as `gradient`
   !> without violating Fortran's actual-argument aliasing rules (an
   !> intent(out)/intent(inout) actual must not overlap any other actual in
   !> the same reference; `operator` and `operator%scratch_gradient` would
   !> overlap if both were passed together). Geometry-only arguments carry no
   !> such risk since they share no storage with the caller's scratch arrays.
   subroutine physical_forward_gradient(nblocks, block_grid, block_coordinate, &
         inverse_block_transpose_m1, values, gradient)
      integer, intent(in) :: nblocks
      integer, intent(in) :: block_grid(3)
      integer, intent(in) :: block_coordinate(:,:)
      real(dblprec), intent(in) :: inverse_block_transpose_m1(3,3)
      real(dblprec), intent(in) :: values(:,:)
      real(dblprec), intent(out) :: gradient(:,:,:)

      integer :: block, direction_index, physical_index, plus_block
      integer :: coordinate(3)

      gradient = 0.0_dblprec
      do block = 1, nblocks
         do direction_index = 1, 3
            coordinate = block_coordinate(:,block)
            coordinate(direction_index) = modulo(coordinate(direction_index)+1, &
               block_grid(direction_index))
            plus_block = regular_spatial_block_id(coordinate,block_grid)
            do physical_index = 1, 3
               gradient(:,physical_index,block) = gradient(:,physical_index,block) + &
                  inverse_block_transpose_m1(physical_index,direction_index) * &
                  (values(:,plus_block)-values(:,block))
            end do
         end do
      end do
   end subroutine physical_forward_gradient

   !> RCG-06A: see physical_forward_gradient's header comment -- the same
   !> operator/scratch aliasing hazard applies here (`result` may be one of
   !> operator's own scratch arrays), so geometry is passed explicitly.
   subroutine add_physical_gradient_transpose(nblocks, block_grid, block_coordinate, &
         inverse_block_transpose_m1, physical_index, values, result, scale)
      integer, intent(in) :: nblocks
      integer, intent(in) :: block_grid(3)
      integer, intent(in) :: block_coordinate(:,:)
      real(dblprec), intent(in) :: inverse_block_transpose_m1(3,3)
      integer, intent(in) :: physical_index
      real(dblprec), intent(in) :: values(:,:)
      real(dblprec), intent(inout) :: result(:,:)
      real(dblprec), intent(in) :: scale

      integer :: block, direction_index, plus_block
      integer :: coordinate(3)
      real(dblprec) :: coefficient

      do block = 1, nblocks
         do direction_index = 1, 3
            coefficient = scale * inverse_block_transpose_m1(physical_index,direction_index)
            coordinate = block_coordinate(:,block)
            coordinate(direction_index) = modulo(coordinate(direction_index)+1, &
               block_grid(direction_index))
            plus_block = regular_spatial_block_id(coordinate,block_grid)
            result(:,plus_block) = result(:,plus_block) + coefficient*values(:,block)
            result(:,block) = result(:,block) - coefficient*values(:,block)
         end do
      end do
   end subroutine add_physical_gradient_transpose

   pure function cross3(left, right) result(cross)
      real(dblprec), intent(in) :: left(3), right(3)
      real(dblprec) :: cross(3)

      cross(1) = left(2)*right(3)-left(3)*right(2)
      cross(2) = left(3)*right(1)-left(1)*right(3)
      cross(3) = left(1)*right(2)-left(2)*right(1)
   end function cross3

   pure real(dblprec) function determinant3(matrix) result(determinant)
      real(dblprec), intent(in) :: matrix(3,3)

      determinant = matrix(1,1)*(matrix(2,2)*matrix(3,3)-matrix(2,3)*matrix(3,2)) - &
         matrix(1,2)*(matrix(2,1)*matrix(3,3)-matrix(2,3)*matrix(3,1)) + &
         matrix(1,3)*(matrix(2,1)*matrix(3,2)-matrix(2,2)*matrix(3,1))
   end function determinant3

   pure function inverse3(matrix) result(inverse)
      real(dblprec), intent(in) :: matrix(3,3)
      real(dblprec) :: inverse(3,3), determinant

      determinant = determinant3(matrix)
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

end module CoarseTensorOperator
