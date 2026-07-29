!-------------------------------------------------------------------------------
! MODULE: stiffnes
!> @brief
!> Intended to use for calculation of the stiffness constant for FM materials
!
!> @details
!> Routines to calculate the ferromagnetic exchange stiffness for ordered systems and random alloys
!> The formalism used here is the one used by Pajda et al. PRB 64, 174402, where the exchange interactions are summed
!> and weighted by the square of the distance to obtain the spin wave stiffness D. This formalism is then exted to multi-sublattice
!> cases and random alloys to be able to explore more complex situations
!> @author
!> Jonathan Chico
!> M. Pereiro ---> Improvment of the regresion method by using Pade aproximants
!> L. Bergqvist -> Added Tc-MFA
!> @copyright
!> GNU Public License
!-------------------------------------------------------------------------------
module stiffness

   use Constants
   use Profiling
   use, intrinsic :: ieee_arithmetic, only : ieee_is_finite

   implicit none

   integer :: eta_max !< Number of convergence parameters for the stiffness
   integer :: eta_min !< Minimum convergence parameters for the stiffness
   character(len=1) :: do_stiffness !< Calculate the spin wave stiffness of a ferromagnet
   character(len=1) :: do_dm_stiffness !< Calculate the DMI spiralization
   character(len=1) :: prn_J0_matrix   !< Print the full site dependent J0 matrix
   real(dblprec) :: A_xc     !< Exchange stiffness constant from rational fit J/m
   real(dblprec) :: M_sat    !< Saturation magnetization muB/m^3
   real(dblprec) :: Dxc_fit  !< Spin wave stiffnes from rational fit
   real(dblprec) :: cell_vol !< Unit cell volume m^3
   real(dblprec) :: A_xc_lsq !< Exchange stiffness constant from LSQ fit
   real(dblprec) :: D_err_fit   !< Error form spin wave stiffnes rational fit
   real(dblprec) :: Dxc_fit_lsq !< Spin wave stiffness LSQ fit
   real(dblprec), dimension(:,:), allocatable :: J0_matrix !< Matrix being used to calculate the Tc-MFA
   real(dblprec), dimension(:,:), allocatable :: D_xc_stiffness_matrix !< Spin wave stiffness tensor rational fit
   real(dblprec), dimension(:,:), allocatable :: A_xc_stiffness_matrix !< Exchange stiffness constant rational fit
   real(dblprec), dimension(:,:), allocatable :: D_xc_stiffness_matrix_lsq !< Spin wave stiffness tensor LSQ fit
   real(dblprec), dimension(:,:), allocatable :: A_xc_stiffness_matrix_lsq !< Exchange stiffness constant LSQ fit
   real(dblprec), dimension(:,:), allocatable :: DM0_mat !< DMI spiralization matrix [meVA]
   real(dblprec), dimension(:,:), allocatable :: DM0_mat_lsq !< DMI spiralization matrix (LSQ fit) [meVA]

   real(dblprec), dimension(:,:,:), allocatable :: J0_matrix_alloy !< Exchange matrix alloy 
   real(dblprec), dimension(:), allocatable :: Axc_fit_alloy !< Exchange stiffness matrix alloy
   real(dblprec), dimension(:), allocatable :: Dxc_fit_alloy !< Spin wave  stiffness matrix alloy
   real(dblprec), dimension(:), allocatable :: Tc_alloy !< Tc-MFA alloy

   ! Coarse-material extraction status and convention enums.  These values are
   ! deliberately independent of the legacy stiffness output flags above.
   integer, parameter :: COARSE_MATERIAL_OK = 0
   integer, parameter :: COARSE_MATERIAL_INVALID_INPUT = 1
   integer, parameter :: COARSE_MATERIAL_FIT_FAILED = 2
   integer, parameter :: COARSE_MATERIAL_VALIDATION_FAILED = 3
   integer, parameter :: COARSE_MATERIAL_RUNTIME_GATED = 4

   integer, parameter :: COARSE_RUNTIME_SINGLE_FM = 1
   integer, parameter :: COARSE_RUNTIME_FERRI_AFM = 2

   integer, parameter :: COARSE_DMI_ENERGY_PLUS = 1
   integer, parameter :: COARSE_FIELD_DERIVATIVE_MINUS = -1
   integer, parameter :: COARSE_MODE_EXTRACTION_UNVALIDATED = 0
   integer, parameter :: COARSE_MODE_EXTRACTION_ACOUSTIC_OPTICAL = 1
   real(dblprec), parameter :: COARSE_MUB_SI = 9.274009994d-24

   !> Directed pair coefficients and convergence controls used by the
   !> side-effect-light calculation service.  Pair energies are coefficients
   !> of the CG-01 Hamiltonian before its directed-pair factor 1/2:
   !>
   !>   E_J = -1/2 sum K_ij m_i.m_j
   !>   E_D = +1/2 sum L_ij.(m_i x m_j)
   !>
   !> K and L are in joules and displacement_m is Cartesian, in metres.
   type :: coarse_lattice_input_type
      integer :: n_basis = 0
      integer :: n_channels = 0
      real(dblprec) :: cell_volume_m3 = 0.0_dblprec
      integer :: fit_first = 1
      integer :: fit_last = 0
      integer, allocatable :: basis_to_channel(:)
      real(dblprec), allocatable :: basis_moment_mub(:)
      real(dblprec), allocatable :: channel_gamma(:)
      real(dblprec), allocatable :: channel_damping(:)
      real(dblprec), allocatable :: eta_inverse_m(:)
      integer, allocatable :: exchange_source_basis(:)
      integer, allocatable :: exchange_target_basis(:)
      real(dblprec), allocatable :: exchange_displacement_m(:,:)
      real(dblprec), allocatable :: exchange_pair_energy_j(:)
      integer, allocatable :: dmi_source_basis(:)
      integer, allocatable :: dmi_target_basis(:)
      real(dblprec), allocatable :: dmi_displacement_m(:,:)
      real(dblprec), allocatable :: dmi_pair_energy_j(:,:)
   end type coarse_lattice_input_type

   !> Raw real-space sums.  No eigenvalue or sublattice-mode reduction is
   !> applied.  The final index of A and D is the regularization sample.
   type :: coarse_lattice_sums_type
      integer :: n_basis = 0
      integer :: n_channels = 0
      integer :: n_eta = 0
      integer :: fit_first = 1
      integer :: fit_last = 0
      real(dblprec), allocatable :: eta_inverse_m(:)
      real(dblprec), allocatable :: basis_local_exchange(:,:)
      real(dblprec), allocatable :: channel_local_exchange(:,:)
      real(dblprec), allocatable :: basis_exchange_stiffness(:,:,:,:,:)
      real(dblprec), allocatable :: channel_exchange_stiffness(:,:,:,:,:)
      real(dblprec), allocatable :: basis_spiralization(:,:,:,:,:)
      real(dblprec), allocatable :: channel_spiralization(:,:,:,:,:)
   end type coarse_lattice_sums_type

   type :: coarse_material_metadata_type
      character(len=16) :: energy_unit = 'J'
      character(len=16) :: length_unit = 'm'
      character(len=16) :: moment_unit = 'muB'
      character(len=16) :: field_unit = 'T'
      character(len=24) :: channel_gamma_unit = 's^-1 T^-1'
      character(len=24) :: channel_damping_unit = 'dimensionless'
      character(len=24) :: local_exchange_unit = 'J m^-3'
      character(len=24) :: exchange_stiffness_unit = 'J m^-1'
      character(len=24) :: spiralization_unit = 'J m^-2'
      character(len=96) :: exchange_order = &
         '(space_p,space_q,left_channel,right_channel)'
      character(len=96) :: spiralization_order = &
         '(cross_spin_k,space_p,left_channel,right_channel)'
      character(len=256) :: hamiltonian_convention = &
         'E_J=-1/2 sum K_ij m_i.m_j; E_D=+1/2 sum L_ij.(m_i x m_j)'
      character(len=128) :: pair_convention = &
         'central-cell outgoing directed pairs; reciprocal partners retained'
      character(len=64) :: regularization = 'exp(-eta*Cartesian_pair_distance)'
      character(len=32) :: convention_version = 'CG-01-approved-v1'
      integer :: dmi_energy_sign = COARSE_DMI_ENERGY_PLUS
      integer :: field_derivative_sign = COARSE_FIELD_DERIVATIVE_MINUS
      logical :: automatic_moment_sign_inference = .false.
   end type coarse_material_metadata_type

   type :: coarse_material_diagnostics_type
      character(len=48) :: fit_method = 'uninitialized'
      logical :: fit_performed = .false.
      integer :: fit_first = 0
      integer :: fit_last = 0
      integer :: fit_sample_count = 0
      real(dblprec) :: eta_min_inverse_m = 0.0_dblprec
      real(dblprec) :: eta_max_inverse_m = 0.0_dblprec
      real(dblprec) :: exchange_reciprocity_max_j_per_m = 0.0_dblprec
      real(dblprec) :: regularization_exchange_delta_j_per_m = 0.0_dblprec
      real(dblprec) :: regularization_dmi_delta_j_per_m2 = 0.0_dblprec
      real(dblprec), allocatable :: exchange_fit_coeff(:,:,:,:,:)
      real(dblprec), allocatable :: spiralization_fit_coeff(:,:,:,:,:)
      real(dblprec), allocatable :: exchange_fit_rms(:,:,:,:)
      real(dblprec), allocatable :: spiralization_fit_rms(:,:,:,:)
      logical :: small_q_energy_validated = .false.
      integer :: small_q_sample_count = 0
      real(dblprec) :: small_q_max_inverse_m = 0.0_dblprec
      real(dblprec) :: small_q_max_abs_error_j_per_m3 = 0.0_dblprec
      real(dblprec) :: small_q_max_scaled_error = 0.0_dblprec
      real(dblprec), allocatable :: small_q_atomistic_energy_j_per_m3(:)
      real(dblprec), allocatable :: small_q_continuum_energy_j_per_m3(:)
      real(dblprec), allocatable :: small_q_abs_error_j_per_m3(:)
      logical :: channel_dynamics_parameters_validated = .false.
      integer :: mode_extraction = COARSE_MODE_EXTRACTION_UNVALIDATED
      logical :: acoustic_mode_extraction_validated = .false.
      logical :: optical_mode_extraction_validated = .false.
      character(len=256) :: runtime_gate_reason = &
         'small-q validation has not been performed'
   end type coarse_material_diagnostics_type

   !> Typed material result consumed by future coarse operators.  The raw
   !> member retains every basis/channel regularized sum before fitting.
   type :: coarse_material_type
      logical :: ready = .false.
      integer :: n_basis = 0
      integer :: n_channels = 0
      real(dblprec) :: cell_volume_m3 = 0.0_dblprec
      integer, allocatable :: basis_to_channel(:)
      real(dblprec), allocatable :: basis_moment_mub(:)
      real(dblprec), allocatable :: channel_moment_mub(:)
      real(dblprec), allocatable :: channel_gamma(:)
      real(dblprec), allocatable :: channel_damping(:)
      real(dblprec), allocatable :: local_exchange(:,:)
      real(dblprec), allocatable :: exchange_stiffness(:,:,:,:)
      real(dblprec), allocatable :: spiralization(:,:,:,:)
      type(coarse_lattice_sums_type) :: raw
      type(coarse_material_metadata_type) :: metadata
      type(coarse_material_diagnostics_type) :: diagnostics
   end type coarse_material_type

   public :: coarse_lattice_input_type
   public :: coarse_lattice_sums_type
   public :: coarse_material_metadata_type
   public :: coarse_material_diagnostics_type
   public :: coarse_material_type
   public :: calculate_coarse_lattice_sums
   public :: fit_coarse_material
   public :: extract_coarse_material
   public :: extract_coarse_material_from_uppasd
   public :: validate_coarse_material_small_q
   public :: validate_coarse_material_two_sublattice_modes
   public :: coarse_material_runtime_status
   public


contains

   !---------------------------------------------------------------------------
   !> Validate the allocation, indexing, units, and convergence controls of a
   !> typed lattice input.  This routine performs no I/O and changes no module
   !> state.
   !---------------------------------------------------------------------------
   subroutine validate_coarse_lattice_input(input, status, diagnostic)
      type(coarse_lattice_input_type), intent(in) :: input
      integer, intent(out) :: status
      character(len=*), intent(out) :: diagnostic

      integer :: channel, n_exchange, n_dmi, n_eta

      status = COARSE_MATERIAL_INVALID_INPUT
      diagnostic = ''

      if (input%n_basis <= 0 .or. input%n_channels <= 0) then
         diagnostic = 'Coarse-material extraction requires positive basis and channel counts'
         return
      end if
      if (.not. ieee_is_finite(input%cell_volume_m3) .or. &
          input%cell_volume_m3 <= tiny(input%cell_volume_m3)) then
         diagnostic = 'Coarse-material extraction requires a finite positive SI cell volume'
         return
      end if
      if (.not. allocated(input%basis_to_channel) .or. &
          .not. allocated(input%basis_moment_mub) .or. &
          .not. allocated(input%eta_inverse_m)) then
         diagnostic = 'Basis map, basis moments, and convergence parameters are mandatory'
         return
      end if
      if (size(input%basis_to_channel) /= input%n_basis .or. &
          size(input%basis_moment_mub) /= input%n_basis) then
         diagnostic = 'Basis map and moment arrays must contain exactly n_basis entries'
         return
      end if
      if (any(input%basis_to_channel == 0) .or. &
          any(input%basis_to_channel < -1) .or. &
          any(input%basis_to_channel > input%n_channels)) then
         diagnostic = 'Basis sites must map to a positive dynamical channel, or -1 if nonmagnetic'
         return
      end if
      do channel = 1, input%n_channels
         if (count(input%basis_to_channel == channel) == 0) then
            diagnostic = 'Dynamical channel ids must be contiguous and nonempty'
            return
         end if
      end do
      if (.not. all(ieee_is_finite(input%basis_moment_mub)) .or. &
          any(input%basis_moment_mub < 0.0_dblprec)) then
         diagnostic = 'Basis moment magnitudes must be finite and nonnegative in muB'
         return
      end if
      if (any(input%basis_moment_mub <= 0.0_dblprec .and. &
              input%basis_to_channel > 0)) then
         diagnostic = 'Every magnetic dynamical-channel basis site requires a positive moment magnitude'
         return
      end if
      if (allocated(input%channel_gamma)) then
         if (size(input%channel_gamma) /= input%n_channels .or. &
             .not. all(ieee_is_finite(input%channel_gamma))) then
            diagnostic = 'Channel gyromagnetic factors must be finite and have n_channels entries'
            return
         end if
      end if
      if (allocated(input%channel_damping)) then
         if (size(input%channel_damping) /= input%n_channels .or. &
             .not. all(ieee_is_finite(input%channel_damping))) then
            diagnostic = 'Channel damping must be finite and have n_channels entries'
            return
         end if
      end if

      n_eta = size(input%eta_inverse_m)
      if (n_eta <= 0 .or. .not. all(ieee_is_finite(input%eta_inverse_m)) .or. &
          any(input%eta_inverse_m < 0.0_dblprec)) then
         diagnostic = 'At least one finite nonnegative regularization eta is required'
         return
      end if
      if (input%fit_first < 1 .or. input%fit_last < input%fit_first .or. &
          input%fit_last > n_eta) then
         diagnostic = 'The convergence fit range must lie inside eta_inverse_m'
         return
      end if

      if (.not. allocated(input%exchange_source_basis) .or. &
          .not. allocated(input%exchange_target_basis) .or. &
          .not. allocated(input%exchange_displacement_m) .or. &
          .not. allocated(input%exchange_pair_energy_j)) then
         diagnostic = 'The exchange directed-pair table must be allocated, including for zero pairs'
         return
      end if
      n_exchange = size(input%exchange_source_basis)
      if (size(input%exchange_target_basis) /= n_exchange .or. &
          size(input%exchange_displacement_m,1) /= 3 .or. &
          size(input%exchange_displacement_m,2) /= n_exchange .or. &
          size(input%exchange_pair_energy_j) /= n_exchange) then
         diagnostic = 'Exchange directed-pair table extents are inconsistent'
         return
      end if
      if (n_exchange > 0) then
         if (any(input%exchange_source_basis < 1) .or. &
             any(input%exchange_source_basis > input%n_basis) .or. &
             any(input%exchange_target_basis < 1) .or. &
             any(input%exchange_target_basis > input%n_basis)) then
            diagnostic = 'Exchange pair basis indices lie outside 1:n_basis'
            return
         end if
         if (.not. all(ieee_is_finite(input%exchange_displacement_m)) .or. &
             .not. all(ieee_is_finite(input%exchange_pair_energy_j))) then
            diagnostic = 'Exchange pair displacements and energies must be finite SI values'
            return
         end if
      end if

      if (.not. allocated(input%dmi_source_basis) .or. &
          .not. allocated(input%dmi_target_basis) .or. &
          .not. allocated(input%dmi_displacement_m) .or. &
          .not. allocated(input%dmi_pair_energy_j)) then
         diagnostic = 'The DMI directed-pair table must be allocated, including for zero pairs'
         return
      end if
      n_dmi = size(input%dmi_source_basis)
      if (size(input%dmi_target_basis) /= n_dmi .or. &
          size(input%dmi_displacement_m,1) /= 3 .or. &
          size(input%dmi_displacement_m,2) /= n_dmi .or. &
          size(input%dmi_pair_energy_j,1) /= 3 .or. &
          size(input%dmi_pair_energy_j,2) /= n_dmi) then
         diagnostic = 'DMI directed-pair table extents are inconsistent'
         return
      end if
      if (n_dmi > 0) then
         if (any(input%dmi_source_basis < 1) .or. &
             any(input%dmi_source_basis > input%n_basis) .or. &
             any(input%dmi_target_basis < 1) .or. &
             any(input%dmi_target_basis > input%n_basis)) then
            diagnostic = 'DMI pair basis indices lie outside 1:n_basis'
            return
         end if
         if (.not. all(ieee_is_finite(input%dmi_displacement_m)) .or. &
             .not. all(ieee_is_finite(input%dmi_pair_energy_j))) then
            diagnostic = 'DMI pair displacements and energies must be finite SI values'
            return
         end if
      end if

      status = COARSE_MATERIAL_OK
      diagnostic = ''
   end subroutine validate_coarse_lattice_input

   !---------------------------------------------------------------------------
   !> Calculate raw, regularized real-space lattice sums in Cartesian SI units.
   !>
   !> For directed exchange coefficient K_ij and DMI coefficient L_ij,
   !>
   !> A_pq^(ab) = (1/4V) sum_(i in a,j in b) K_ij r_p r_q
   !> D_kp^(ab) = (1/2V) sum_(i in a,j in b) L_ij,k r_p .
   !>
   !> The factors follow directly from the CG-01 directed-pair Hamiltonian.
   !> No file unit, eigenproblem, fitting routine, or module-global result is
   !> touched.
   !---------------------------------------------------------------------------
   subroutine calculate_coarse_lattice_sums(input, sums, status, diagnostic)
      type(coarse_lattice_input_type), intent(in) :: input
      type(coarse_lattice_sums_type), intent(out) :: sums
      integer, intent(out) :: status
      character(len=*), intent(out) :: diagnostic

      integer :: pair, eta_index, p, q, k
      integer :: left_basis, right_basis, left_channel, right_channel
      real(dblprec) :: distance_m, weight, volume
      real(dblprec), allocatable :: basis_local_directed(:,:)
      real(dblprec), allocatable :: channel_local_directed(:,:)

      call validate_coarse_lattice_input(input, status, diagnostic)
      if (status /= COARSE_MATERIAL_OK) return

      sums%n_basis = input%n_basis
      sums%n_channels = input%n_channels
      sums%n_eta = size(input%eta_inverse_m)
      sums%fit_first = input%fit_first
      sums%fit_last = input%fit_last
      allocate(sums%eta_inverse_m(sums%n_eta))
      allocate(sums%basis_local_exchange(input%n_basis,input%n_basis))
      allocate(sums%channel_local_exchange(input%n_channels,input%n_channels))
      allocate(sums%basis_exchange_stiffness(3,3,input%n_basis,input%n_basis,sums%n_eta))
      allocate(sums%channel_exchange_stiffness(3,3,input%n_channels,input%n_channels,sums%n_eta))
      allocate(sums%basis_spiralization(3,3,input%n_basis,input%n_basis,sums%n_eta))
      allocate(sums%channel_spiralization(3,3,input%n_channels,input%n_channels,sums%n_eta))
      allocate(basis_local_directed(input%n_basis,input%n_basis))
      allocate(channel_local_directed(input%n_channels,input%n_channels))

      sums%eta_inverse_m = input%eta_inverse_m
      sums%basis_local_exchange = 0.0_dblprec
      sums%channel_local_exchange = 0.0_dblprec
      sums%basis_exchange_stiffness = 0.0_dblprec
      sums%channel_exchange_stiffness = 0.0_dblprec
      sums%basis_spiralization = 0.0_dblprec
      sums%channel_spiralization = 0.0_dblprec
      basis_local_directed = 0.0_dblprec
      channel_local_directed = 0.0_dblprec
      volume = input%cell_volume_m3

      do pair = 1, size(input%exchange_source_basis)
         left_basis = input%exchange_source_basis(pair)
         right_basis = input%exchange_target_basis(pair)
         left_channel = input%basis_to_channel(left_basis)
         right_channel = input%basis_to_channel(right_basis)
         if (left_basis /= right_basis) then
            basis_local_directed(left_basis,right_basis) = &
               basis_local_directed(left_basis,right_basis) - &
               0.5_dblprec * input%exchange_pair_energy_j(pair) / volume
         end if
         if (left_channel > 0 .and. right_channel > 0 .and. &
             left_channel /= right_channel) then
            channel_local_directed(left_channel,right_channel) = &
               channel_local_directed(left_channel,right_channel) - &
               0.5_dblprec * input%exchange_pair_energy_j(pair) / volume
         end if
         distance_m = sqrt(sum(input%exchange_displacement_m(:,pair)**2))
         do eta_index = 1, sums%n_eta
            weight = exp(-input%eta_inverse_m(eta_index) * distance_m)
            do p = 1, 3
               do q = 1, 3
                  sums%basis_exchange_stiffness(p,q,left_basis,right_basis,eta_index) = &
                     sums%basis_exchange_stiffness(p,q,left_basis,right_basis,eta_index) + &
                     0.25_dblprec * input%exchange_pair_energy_j(pair) * &
                     input%exchange_displacement_m(p,pair) * &
                     input%exchange_displacement_m(q,pair) * weight / volume
                  if (left_channel > 0 .and. right_channel > 0) then
                     sums%channel_exchange_stiffness(p,q,left_channel,right_channel,eta_index) = &
                        sums%channel_exchange_stiffness(p,q,left_channel,right_channel,eta_index) + &
                        0.25_dblprec * input%exchange_pair_energy_j(pair) * &
                        input%exchange_displacement_m(p,pair) * &
                        input%exchange_displacement_m(q,pair) * weight / volume
                  end if
               end do
            end do
         end do
      end do

      do left_basis = 1, input%n_basis
         do right_basis = left_basis + 1, input%n_basis
            sums%basis_local_exchange(left_basis,right_basis) = &
               basis_local_directed(left_basis,right_basis) + &
               basis_local_directed(right_basis,left_basis)
            sums%basis_local_exchange(right_basis,left_basis) = &
               sums%basis_local_exchange(left_basis,right_basis)
         end do
      end do
      do left_channel = 1, input%n_channels
         do right_channel = left_channel + 1, input%n_channels
            sums%channel_local_exchange(left_channel,right_channel) = &
               channel_local_directed(left_channel,right_channel) + &
               channel_local_directed(right_channel,left_channel)
            sums%channel_local_exchange(right_channel,left_channel) = &
               sums%channel_local_exchange(left_channel,right_channel)
         end do
      end do

      do pair = 1, size(input%dmi_source_basis)
         left_basis = input%dmi_source_basis(pair)
         right_basis = input%dmi_target_basis(pair)
         left_channel = input%basis_to_channel(left_basis)
         right_channel = input%basis_to_channel(right_basis)
         distance_m = sqrt(sum(input%dmi_displacement_m(:,pair)**2))
         do eta_index = 1, sums%n_eta
            weight = exp(-input%eta_inverse_m(eta_index) * distance_m)
            do k = 1, 3
               do p = 1, 3
                  sums%basis_spiralization(k,p,left_basis,right_basis,eta_index) = &
                     sums%basis_spiralization(k,p,left_basis,right_basis,eta_index) + &
                     0.5_dblprec * input%dmi_pair_energy_j(k,pair) * &
                     input%dmi_displacement_m(p,pair) * weight / volume
                  if (left_channel > 0 .and. right_channel > 0) then
                     sums%channel_spiralization(k,p,left_channel,right_channel,eta_index) = &
                        sums%channel_spiralization(k,p,left_channel,right_channel,eta_index) + &
                        0.5_dblprec * input%dmi_pair_energy_j(k,pair) * &
                        input%dmi_displacement_m(p,pair) * weight / volume
                  end if
               end do
            end do
         end do
      end do

      status = COARSE_MATERIAL_OK
      diagnostic = ''
   end subroutine calculate_coarse_lattice_sums

   !---------------------------------------------------------------------------
   !> Fit one regularized lattice-sum series with a quadratic least-squares
   !> extrapolation to eta=0.  Scaling eta before forming the normal equations
   !> avoids powers of SI inverse metres in the solve.
   !---------------------------------------------------------------------------
   subroutine fit_regularization_series(x, y, coefficient, rms, status)
      real(dblprec), intent(in) :: x(:), y(:)
      real(dblprec), intent(out) :: coefficient(3)
      real(dblprec), intent(out) :: rms
      integer, intent(out) :: status

      integer :: i, row, col, solve_status
      real(dblprec) :: xscale, z, prediction
      real(dblprec) :: normal(3,3), rhs(3), scaled_coefficient(3)

      coefficient = 0.0_dblprec
      rms = 0.0_dblprec
      status = COARSE_MATERIAL_FIT_FAILED
      if (size(x) /= size(y) .or. size(x) < 3) return

      xscale = maxval(abs(x))
      if (xscale <= tiny(xscale)) return
      normal = 0.0_dblprec
      rhs = 0.0_dblprec
      do i = 1, size(x)
         z = x(i) / xscale
         do row = 1, 3
            rhs(row) = rhs(row) + y(i) * z**(row-1)
            do col = 1, 3
               normal(row,col) = normal(row,col) + z**(row+col-2)
            end do
         end do
      end do
      call solve_dense_3x3(normal, rhs, scaled_coefficient, solve_status)
      if (solve_status /= COARSE_MATERIAL_OK) return

      coefficient(1) = scaled_coefficient(1)
      coefficient(2) = scaled_coefficient(2) / xscale
      coefficient(3) = scaled_coefficient(3) / (xscale*xscale)
      do i = 1, size(x)
         prediction = coefficient(1) + coefficient(2)*x(i) + coefficient(3)*x(i)*x(i)
         rms = rms + (prediction-y(i))**2
      end do
      rms = sqrt(rms / real(size(x),dblprec))
      status = COARSE_MATERIAL_OK
   end subroutine fit_regularization_series

   subroutine solve_dense_3x3(matrix, rhs, solution, status)
      real(dblprec), intent(in) :: matrix(3,3), rhs(3)
      real(dblprec), intent(out) :: solution(3)
      integer, intent(out) :: status

      integer :: column, row, pivot
      real(dblprec) :: augmented(3,4), row_copy(4), scale, tolerance

      augmented(:,1:3) = matrix
      augmented(:,4) = rhs
      solution = 0.0_dblprec
      scale = maxval(abs(matrix))
      tolerance = 128.0_dblprec * epsilon(scale) * max(1.0_dblprec,scale)
      do column = 1, 3
         pivot = column - 1 + maxloc(abs(augmented(column:3,column)),dim=1)
         if (abs(augmented(pivot,column)) <= tolerance) then
            status = COARSE_MATERIAL_FIT_FAILED
            return
         end if
         if (pivot /= column) then
            row_copy = augmented(column,:)
            augmented(column,:) = augmented(pivot,:)
            augmented(pivot,:) = row_copy
         end if
         augmented(column,:) = augmented(column,:) / augmented(column,column)
         do row = 1, 3
            if (row /= column) then
               augmented(row,:) = augmented(row,:) - &
                  augmented(row,column) * augmented(column,:)
            end if
         end do
      end do
      solution = augmented(:,4)
      status = COARSE_MATERIAL_OK
   end subroutine solve_dense_3x3

   !---------------------------------------------------------------------------
   !> Convert raw uncollapsed lattice sums into a typed coarse material.
   !> Fitting is deliberately separate from calculation so callers can inspect
   !> every eta sample or apply another reviewed extrapolation.
   !---------------------------------------------------------------------------
   subroutine fit_coarse_material(input, sums, material, status, diagnostic)
      type(coarse_lattice_input_type), intent(in) :: input
      type(coarse_lattice_sums_type), intent(in) :: sums
      type(coarse_material_type), intent(out) :: material
      integer, intent(out) :: status
      character(len=*), intent(out) :: diagnostic

      integer :: p, q, left_channel, right_channel, basis
      integer :: fit_status, input_status, nfit, zero_location(1), zero_index
      real(dblprec) :: coefficient(3), rms, zero_tolerance
      real(dblprec) :: reciprocity_error
      logical :: use_fit

      status = COARSE_MATERIAL_INVALID_INPUT
      diagnostic = ''
      call validate_coarse_lattice_input(input, input_status, diagnostic)
      if (input_status /= COARSE_MATERIAL_OK) return
      if (sums%n_basis /= input%n_basis .or. sums%n_channels /= input%n_channels .or. &
          sums%n_eta /= size(input%eta_inverse_m)) then
         diagnostic = 'Raw lattice sums do not match the typed extraction input'
         return
      end if
      if (.not. allocated(sums%eta_inverse_m) .or. &
          .not. allocated(sums%basis_local_exchange) .or. &
          .not. allocated(sums%channel_local_exchange) .or. &
          .not. allocated(sums%basis_exchange_stiffness) .or. &
          .not. allocated(sums%channel_exchange_stiffness) .or. &
          .not. allocated(sums%basis_spiralization) .or. &
          .not. allocated(sums%channel_spiralization)) then
         diagnostic = 'Raw lattice sums are incomplete'
         return
      end if
      if (size(sums%eta_inverse_m) /= sums%n_eta .or. &
          any(shape(sums%basis_local_exchange) /= (/input%n_basis,input%n_basis/)) .or. &
          any(shape(sums%channel_local_exchange) /= (/input%n_channels,input%n_channels/)) .or. &
          any(shape(sums%basis_exchange_stiffness) /= &
             (/3,3,input%n_basis,input%n_basis,sums%n_eta/)) .or. &
          any(shape(sums%channel_exchange_stiffness) /= &
             (/3,3,input%n_channels,input%n_channels,sums%n_eta/)) .or. &
          any(shape(sums%basis_spiralization) /= &
             (/3,3,input%n_basis,input%n_basis,sums%n_eta/)) .or. &
          any(shape(sums%channel_spiralization) /= &
             (/3,3,input%n_channels,input%n_channels,sums%n_eta/))) then
         diagnostic = 'Raw lattice-sum array extents are inconsistent'
         return
      end if
      if (sums%fit_first < 1 .or. sums%fit_last > sums%n_eta .or. &
          sums%fit_last < sums%fit_first) then
         diagnostic = 'Raw lattice sums contain an invalid convergence fit range'
         return
      end if

      nfit = sums%fit_last - sums%fit_first + 1
      use_fit = nfit >= 3
      zero_location = minloc(abs(sums%eta_inverse_m))
      zero_index = zero_location(1)
      zero_tolerance = 64.0_dblprec * epsilon(1.0_dblprec) * &
         max(1.0_dblprec,maxval(abs(sums%eta_inverse_m)))
      if (.not. use_fit .and. abs(sums%eta_inverse_m(zero_index)) > zero_tolerance) then
         diagnostic = 'Fewer than three fit samples require an explicit eta=0 lattice sum'
         status = COARSE_MATERIAL_FIT_FAILED
         return
      end if

      material%n_basis = input%n_basis
      material%n_channels = input%n_channels
      material%cell_volume_m3 = input%cell_volume_m3
      allocate(material%basis_to_channel(input%n_basis))
      allocate(material%basis_moment_mub(input%n_basis))
      allocate(material%channel_moment_mub(input%n_channels))
      allocate(material%channel_gamma(input%n_channels))
      allocate(material%channel_damping(input%n_channels))
      allocate(material%local_exchange(input%n_channels,input%n_channels))
      allocate(material%exchange_stiffness(3,3,input%n_channels,input%n_channels))
      allocate(material%spiralization(3,3,input%n_channels,input%n_channels))
      material%basis_to_channel = input%basis_to_channel
      material%basis_moment_mub = input%basis_moment_mub
      material%channel_moment_mub = 0.0_dblprec
      do basis = 1, input%n_basis
         if (input%basis_to_channel(basis) > 0) then
            material%channel_moment_mub(input%basis_to_channel(basis)) = &
               material%channel_moment_mub(input%basis_to_channel(basis)) + &
               input%basis_moment_mub(basis)
         end if
      end do
      material%channel_gamma = 0.0_dblprec
      material%channel_damping = 0.0_dblprec
      if (allocated(input%channel_gamma)) material%channel_gamma = input%channel_gamma
      if (allocated(input%channel_damping)) material%channel_damping = input%channel_damping
      material%diagnostics%channel_dynamics_parameters_validated = &
         allocated(input%channel_gamma) .and. allocated(input%channel_damping)
      if (material%diagnostics%channel_dynamics_parameters_validated) then
         material%diagnostics%channel_dynamics_parameters_validated = &
            all(input%channel_gamma > 0.0_dblprec) .and. &
            all(input%channel_damping >= 0.0_dblprec)
      end if
      material%local_exchange = sums%channel_local_exchange
      material%exchange_stiffness = 0.0_dblprec
      material%spiralization = 0.0_dblprec
      material%raw = sums

      allocate(material%diagnostics%exchange_fit_coeff(3,3,3,input%n_channels,input%n_channels))
      allocate(material%diagnostics%spiralization_fit_coeff(3,3,3,input%n_channels,input%n_channels))
      allocate(material%diagnostics%exchange_fit_rms(3,3,input%n_channels,input%n_channels))
      allocate(material%diagnostics%spiralization_fit_rms(3,3,input%n_channels,input%n_channels))
      material%diagnostics%exchange_fit_coeff = 0.0_dblprec
      material%diagnostics%spiralization_fit_coeff = 0.0_dblprec
      material%diagnostics%exchange_fit_rms = 0.0_dblprec
      material%diagnostics%spiralization_fit_rms = 0.0_dblprec
      material%diagnostics%fit_first = sums%fit_first
      material%diagnostics%fit_last = sums%fit_last
      material%diagnostics%fit_sample_count = nfit
      material%diagnostics%eta_min_inverse_m = minval( &
         sums%eta_inverse_m(sums%fit_first:sums%fit_last))
      material%diagnostics%eta_max_inverse_m = maxval( &
         sums%eta_inverse_m(sums%fit_first:sums%fit_last))

      do left_channel = 1, input%n_channels
         do right_channel = 1, input%n_channels
            do p = 1, 3
               do q = 1, 3
                  if (use_fit) then
                     call fit_regularization_series( &
                        sums%eta_inverse_m(sums%fit_first:sums%fit_last), &
                        sums%channel_exchange_stiffness(p,q,left_channel,right_channel, &
                           sums%fit_first:sums%fit_last), coefficient, rms, fit_status)
                     if (fit_status /= COARSE_MATERIAL_OK) then
                        diagnostic = 'Quadratic exchange regularization fit is singular'
                        status = COARSE_MATERIAL_FIT_FAILED
                        return
                     end if
                  else
                     coefficient = 0.0_dblprec
                     coefficient(1) = sums%channel_exchange_stiffness( &
                        p,q,left_channel,right_channel,zero_index)
                     rms = 0.0_dblprec
                  end if
                  material%exchange_stiffness(p,q,left_channel,right_channel) = coefficient(1)
                  material%diagnostics%exchange_fit_coeff(:,p,q,left_channel,right_channel) = coefficient
                  material%diagnostics%exchange_fit_rms(p,q,left_channel,right_channel) = rms

                  if (use_fit) then
                     call fit_regularization_series( &
                        sums%eta_inverse_m(sums%fit_first:sums%fit_last), &
                        sums%channel_spiralization(p,q,left_channel,right_channel, &
                           sums%fit_first:sums%fit_last), coefficient, rms, fit_status)
                     if (fit_status /= COARSE_MATERIAL_OK) then
                        diagnostic = 'Quadratic DMI regularization fit is singular'
                        status = COARSE_MATERIAL_FIT_FAILED
                        return
                     end if
                  else
                     coefficient = 0.0_dblprec
                     coefficient(1) = sums%channel_spiralization( &
                        p,q,left_channel,right_channel,zero_index)
                     rms = 0.0_dblprec
                  end if
                  material%spiralization(p,q,left_channel,right_channel) = coefficient(1)
                  material%diagnostics%spiralization_fit_coeff(:,p,q,left_channel,right_channel) = coefficient
                  material%diagnostics%spiralization_fit_rms(p,q,left_channel,right_channel) = rms
               end do
            end do
         end do
      end do

      material%diagnostics%fit_performed = use_fit
      if (use_fit) then
         material%diagnostics%fit_method = 'quadratic least squares in eta'
      else
         material%diagnostics%fit_method = 'direct eta=0 lattice sum'
      end if
      material%diagnostics%regularization_exchange_delta_j_per_m = &
         maxval(abs(material%exchange_stiffness - &
         sums%channel_exchange_stiffness(:,:,:,:,zero_index)))
      material%diagnostics%regularization_dmi_delta_j_per_m2 = &
         maxval(abs(material%spiralization - &
         sums%channel_spiralization(:,:,:,:,zero_index)))

      reciprocity_error = 0.0_dblprec
      do left_channel = 1, input%n_channels
         do right_channel = 1, input%n_channels
            do p = 1, 3
               do q = 1, 3
                  reciprocity_error = max(reciprocity_error, abs( &
                     material%exchange_stiffness(p,q,left_channel,right_channel) - &
                     material%exchange_stiffness(q,p,right_channel,left_channel)))
               end do
            end do
         end do
      end do
      material%diagnostics%exchange_reciprocity_max_j_per_m = reciprocity_error
      if (input%n_channels > 1) then
         material%diagnostics%runtime_gate_reason = &
            'ferri/AFM runtime is gated: acoustic and optical mode extraction is unvalidated'
      else
         material%diagnostics%runtime_gate_reason = &
            'single-channel runtime requires a passing direct small-q energy validation'
      end if
      material%ready = .true.
      status = COARSE_MATERIAL_OK
      diagnostic = ''
   end subroutine fit_coarse_material

   !> Convenience composition of the independent calculation and fitting APIs.
   subroutine extract_coarse_material(input, material, status, diagnostic)
      type(coarse_lattice_input_type), intent(in) :: input
      type(coarse_material_type), intent(out) :: material
      integer, intent(out) :: status
      character(len=*), intent(out) :: diagnostic

      type(coarse_lattice_sums_type) :: sums

      call calculate_coarse_lattice_sums(input, sums, status, diagnostic)
      if (status /= COARSE_MATERIAL_OK) return
      call fit_coarse_material(input, sums, material, status, diagnostic)
   end subroutine extract_coarse_material

   !---------------------------------------------------------------------------
   !> Side-effect-light adapter from the ordered UppASD Hamiltonian arrays to
   !> the typed SI extraction service.  It deliberately does not infer a
   !> ferrimagnetic/AFM channel map from signed input moments.
   !---------------------------------------------------------------------------
   subroutine extract_coarse_material_from_uppasd(NA,N1,N2,N3,Natom,nHam,mconf, &
         Nchmax,max_no_neigh,max_no_dmneigh,eta_min,eta_max,anumb,aham,nlistsize, &
         nlist,dmlistsize,dmlist,alat,C1,C2,C3,BC1,BC2,BC3,coord,ammom_inp,ncoup,dm_vect, &
         basis_to_channel,material,status,diagnostic,channel_gamma,channel_damping)

      integer, intent(in) :: NA, N1, N2, N3, Natom, nHam, mconf, Nchmax
      integer, intent(in) :: max_no_neigh, max_no_dmneigh, eta_min, eta_max
      integer, intent(in) :: anumb(Natom), aham(Natom), nlistsize(nHam)
      integer, intent(in) :: nlist(max_no_neigh,Natom), dmlistsize(nHam)
      integer, intent(in) :: dmlist(max_no_dmneigh,Natom)
      integer, intent(in) :: basis_to_channel(NA)
      real(dblprec), intent(in) :: alat, C1(3), C2(3), C3(3), coord(3,Natom)
      real(dblprec), intent(in) :: ammom_inp(:,:,:), ncoup(:,:,:,:), dm_vect(:,:,:)
      character(len=1), intent(in) :: BC1, BC2, BC3
      type(coarse_material_type), intent(out) :: material
      integer, intent(out) :: status
      character(len=*), intent(out) :: diagnostic
      real(dblprec), intent(in), optional :: channel_gamma(:), channel_damping(:)

      type(coarse_lattice_input_type) :: input
      integer :: i, pair, neighbor, iatom, jatom, iham, countstart
      integer :: n_exchange, n_dmi, n_channels, eta_index
      real(dblprec) :: cell_matrix(3,3), displacement(3)
      real(dblprec) :: moment_i, moment_j

      status = COARSE_MATERIAL_INVALID_INPUT
      diagnostic = ''
      if (NA <= 0 .or. Nchmax <= 0 .or. Natom /= NA*N1*N2*N3 .or. nHam <= 0) then
         diagnostic = 'UppASD coarse extraction requires a regular replicated cell'
         return
      end if
      if (Natom > 1 .and. NA == Natom) then
         diagnostic = 'Typed coarse extraction rejects EXPLICIT_DEVICE geometry with NA=Natom'
         return
      end if
      if (.not. ieee_is_finite(alat) .or. alat <= 0.0_dblprec .or. &
          .not. all(ieee_is_finite(C1)) .or. .not. all(ieee_is_finite(C2)) .or. &
          .not. all(ieee_is_finite(C3)) .or. .not. all(ieee_is_finite(coord))) then
         diagnostic = 'UppASD coarse extraction requires finite SI geometry and positive alat'
         return
      end if
      if (mconf < 1 .or. mconf > size(ammom_inp,3) .or. &
          size(ammom_inp,1) < NA .or. size(ammom_inp,2) < 1) then
         diagnostic = 'UppASD moment configuration is outside ammom_inp'
         return
      end if
      if (size(ncoup,1) < 1 .or. size(ncoup,2) < max_no_neigh .or. &
          size(ncoup,3) < nHam .or. size(ncoup,4) < mconf .or. &
          size(dm_vect,1) < 3 .or. size(dm_vect,2) < max_no_dmneigh .or. &
          size(dm_vect,3) < nHam) then
         diagnostic = 'UppASD Hamiltonian array extents do not match the neighbour tables'
         return
      end if
      if (any(nlistsize < 0) .or. any(nlistsize > max_no_neigh) .or. &
          any(dmlistsize < 0) .or. any(dmlistsize > max_no_dmneigh)) then
         diagnostic = 'UppASD neighbour-list sizes lie outside their allocated extents'
         return
      end if
      if (eta_max < 0 .or. eta_min < 0 .or. eta_min > eta_max) then
         diagnostic = 'UppASD stiffness eta_min:eta_max is invalid'
         return
      end if
      if (abs(mub-COARSE_MUB_SI) > 1.0d-8*COARSE_MUB_SI) then
         diagnostic = 'Typed coarse extraction requires physical SI constants, not reduced units'
         return
      end if
      n_channels = maxval(basis_to_channel)
      if (n_channels <= 0 .or. any(basis_to_channel == 0) .or. &
          any(basis_to_channel < -1)) then
         diagnostic = 'An explicit basis-to-channel map is required; use -1 only for nonmagnetic sites'
         return
      end if
      if (present(channel_gamma)) then
         if (size(channel_gamma) /= n_channels) then
            diagnostic = 'channel_gamma must have max(basis_to_channel) entries'
            return
         end if
      end if
      if (present(channel_damping)) then
         if (size(channel_damping) /= n_channels) then
            diagnostic = 'channel_damping must have max(basis_to_channel) entries'
            return
         end if
      end if

      countstart = (N1/2)*NA + (N2/2)*N1*NA + (N3/2)*N2*N1*NA
      do i = 1, NA
         iatom = countstart + i
         if (iatom < 1 .or. iatom > Natom .or. anumb(iatom) /= i) then
            diagnostic = 'Central-cell atom ordering is not the regular basis-fast UppASD layout'
            return
         end if
         if (aham(iatom) < 1 .or. aham(iatom) > nHam) then
            diagnostic = 'Central-cell Hamiltonian lookup lies outside 1:nHam'
            return
         end if
      end do

      n_exchange = 0
      n_dmi = 0
      do i = 1, NA
         iham = aham(countstart+i)
         n_exchange = n_exchange + nlistsize(iham)
         n_dmi = n_dmi + dmlistsize(iham)
      end do

      input%n_basis = NA
      input%n_channels = n_channels
      cell_matrix(:,1) = alat*C1
      cell_matrix(:,2) = alat*C2
      cell_matrix(:,3) = alat*C3
      input%cell_volume_m3 = abs(determinant_3x3(cell_matrix))
      input%fit_first = eta_min + 1
      input%fit_last = eta_max + 1
      allocate(input%basis_to_channel(NA), input%basis_moment_mub(NA))
      allocate(input%channel_gamma(n_channels), input%channel_damping(n_channels))
      allocate(input%eta_inverse_m(eta_max+1))
      allocate(input%exchange_source_basis(n_exchange))
      allocate(input%exchange_target_basis(n_exchange))
      allocate(input%exchange_displacement_m(3,n_exchange))
      allocate(input%exchange_pair_energy_j(n_exchange))
      allocate(input%dmi_source_basis(n_dmi))
      allocate(input%dmi_target_basis(n_dmi))
      allocate(input%dmi_displacement_m(3,n_dmi))
      allocate(input%dmi_pair_energy_j(3,n_dmi))
      input%basis_to_channel = basis_to_channel
      input%basis_moment_mub = abs(ammom_inp(1:NA,1,mconf))
      input%channel_gamma = 0.0_dblprec
      input%channel_damping = 0.0_dblprec
      if (present(channel_gamma)) input%channel_gamma = channel_gamma
      if (present(channel_damping)) input%channel_damping = channel_damping
      do eta_index = 1, eta_max + 1
         input%eta_inverse_m(eta_index) = 0.10_dblprec * real(eta_index-1,dblprec) / alat
      end do

      pair = 0
      do i = 1, NA
         iatom = countstart + i
         iham = aham(iatom)
         moment_i = abs(ammom_inp(i,1,mconf))
         do neighbor = 1, nlistsize(iham)
            jatom = nlist(neighbor,iatom)
            if (jatom < 1 .or. jatom > Natom) then
               diagnostic = 'Exchange neighbor index lies outside 1:Natom'
               return
            end if
            if (anumb(jatom) < 1 .or. anumb(jatom) > NA) then
               diagnostic = 'Exchange neighbor basis index lies outside 1:NA'
               return
            end if
            call coarse_wrapped_coordinate_difference(N1,N2,N3,C1,C2,C3,BC1,BC2,BC3, &
               coord(:,jatom),coord(:,iatom),displacement)
            moment_j = abs(ammom_inp(anumb(jatom),1,mconf))
            pair = pair + 1
            input%exchange_source_basis(pair) = i
            input%exchange_target_basis(pair) = anumb(jatom)
            input%exchange_displacement_m(:,pair) = alat*displacement
            input%exchange_pair_energy_j(pair) = COARSE_MUB_SI * &
               ncoup(1,neighbor,iham,mconf) * moment_i * moment_j
         end do
      end do

      pair = 0
      do i = 1, NA
         iatom = countstart + i
         iham = aham(iatom)
         moment_i = abs(ammom_inp(i,1,mconf))
         do neighbor = 1, dmlistsize(iham)
            jatom = dmlist(neighbor,iatom)
            if (jatom < 1 .or. jatom > Natom) then
               diagnostic = 'DMI neighbor index lies outside 1:Natom'
               return
            end if
            if (anumb(jatom) < 1 .or. anumb(jatom) > NA) then
               diagnostic = 'DMI neighbor basis index lies outside 1:NA'
               return
            end if
            call coarse_wrapped_coordinate_difference(N1,N2,N3,C1,C2,C3,BC1,BC2,BC3, &
               coord(:,jatom),coord(:,iatom),displacement)
            moment_j = abs(ammom_inp(anumb(jatom),1,mconf))
            pair = pair + 1
            input%dmi_source_basis(pair) = i
            input%dmi_target_basis(pair) = anumb(jatom)
            input%dmi_displacement_m(:,pair) = alat*displacement
            input%dmi_pair_energy_j(:,pair) = COARSE_MUB_SI * &
               dm_vect(:,neighbor,iham) * moment_i * moment_j
         end do
      end do

      call extract_coarse_material(input, material, status, diagnostic)
   end subroutine extract_coarse_material_from_uppasd

   pure subroutine coarse_wrapped_coordinate_difference(N1,N2,N3,C1,C2,C3, &
         BC1,BC2,BC3,target_coordinate,source_coordinate,difference)
      integer, intent(in) :: N1, N2, N3
      real(dblprec), intent(in) :: C1(3), C2(3), C3(3)
      character(len=1), intent(in) :: BC1, BC2, BC3
      real(dblprec), intent(in) :: target_coordinate(3), source_coordinate(3)
      real(dblprec), intent(out) :: difference(3)

      integer :: x, y, z, xmin, xmax, ymin, ymax, zmin, zmax
      real(dblprec) :: original(3), candidate(3)

      original = target_coordinate-source_coordinate
      difference = original
      xmin = 0
      xmax = 0
      ymin = 0
      ymax = 0
      zmin = 0
      zmax = 0
      if (BC1 == 'P') then
         xmin = -1
         xmax = 1
      end if
      if (BC2 == 'P') then
         ymin = -1
         ymax = 1
      end if
      if (BC3 == 'P') then
         zmin = -1
         zmax = 1
      end if
      do z = zmin, zmax
         do y = ymin, ymax
            do x = xmin, xmax
               candidate = original + real(x*N1,dblprec)*C1 + &
                  real(y*N2,dblprec)*C2 + real(z*N3,dblprec)*C3
               if (sum(candidate*candidate) < sum(difference*difference)) &
                  difference = candidate
            end do
         end do
      end do
   end subroutine coarse_wrapped_coordinate_difference

   !---------------------------------------------------------------------------
   !> Compare the fitted tensor model with direct atomistic planar-spiral
   !> energies for caller-supplied small Cartesian q vectors.
   !>
   !> channel_phase fixes the q=0 phase of each channel for each sample;
   !> phases 0 and pi therefore cover acoustic and controlled optical fixtures.
   !> Passing these energy checks does not claim that dynamical acoustic and
   !> optical eigenmodes have been extracted.
   !---------------------------------------------------------------------------
   subroutine validate_coarse_material_small_q(input, material, q_inverse_m, &
         spiral_normal, channel_phase, absolute_tolerance_j_per_m3, &
         relative_tolerance, status, diagnostic)
      type(coarse_lattice_input_type), intent(in) :: input
      type(coarse_material_type), intent(inout) :: material
      real(dblprec), intent(in) :: q_inverse_m(:,:)
      real(dblprec), intent(in) :: spiral_normal(:,:)
      real(dblprec), intent(in) :: channel_phase(:,:)
      real(dblprec), intent(in) :: absolute_tolerance_j_per_m3
      real(dblprec), intent(in) :: relative_tolerance
      integer, intent(out) :: status
      character(len=*), intent(out) :: diagnostic

      integer :: sample, pair, p, q, k, left_channel, right_channel
      integer :: n_samples, input_status
      real(dblprec) :: phase_zero, phase, exact_energy, continuum_energy
      real(dblprec) :: normal_norm, error, scale, allowed, scaled_error
      logical :: all_pass

      status = COARSE_MATERIAL_INVALID_INPUT
      diagnostic = ''
      call validate_coarse_lattice_input(input, input_status, diagnostic)
      if (input_status /= COARSE_MATERIAL_OK) return
      if (.not. material%ready .or. material%n_basis /= input%n_basis .or. &
          material%n_channels /= input%n_channels) then
         diagnostic = 'Small-q validation requires a matching ready coarse material'
         return
      end if
      if (size(q_inverse_m,1) /= 3) then
         diagnostic = 'Small-q vectors must have shape (3,n_samples)'
         return
      end if
      n_samples = size(q_inverse_m,2)
      if (n_samples <= 0 .or. size(spiral_normal,1) /= 3 .or. &
          size(spiral_normal,2) /= n_samples .or. &
          size(channel_phase,1) /= input%n_channels .or. &
          size(channel_phase,2) /= n_samples) then
         diagnostic = 'Small-q normals and channel phases do not match the sample count'
         return
      end if
      if (.not. all(ieee_is_finite(q_inverse_m)) .or. &
          .not. all(ieee_is_finite(spiral_normal)) .or. &
          .not. all(ieee_is_finite(channel_phase)) .or. &
          absolute_tolerance_j_per_m3 < 0.0_dblprec .or. relative_tolerance < 0.0_dblprec) then
         diagnostic = 'Small-q validation inputs and tolerances must be finite and nonnegative'
         return
      end if

      if (allocated(material%diagnostics%small_q_atomistic_energy_j_per_m3)) &
         deallocate(material%diagnostics%small_q_atomistic_energy_j_per_m3)
      if (allocated(material%diagnostics%small_q_continuum_energy_j_per_m3)) &
         deallocate(material%diagnostics%small_q_continuum_energy_j_per_m3)
      if (allocated(material%diagnostics%small_q_abs_error_j_per_m3)) &
         deallocate(material%diagnostics%small_q_abs_error_j_per_m3)
      allocate(material%diagnostics%small_q_atomistic_energy_j_per_m3(n_samples))
      allocate(material%diagnostics%small_q_continuum_energy_j_per_m3(n_samples))
      allocate(material%diagnostics%small_q_abs_error_j_per_m3(n_samples))
      material%diagnostics%small_q_atomistic_energy_j_per_m3 = 0.0_dblprec
      material%diagnostics%small_q_continuum_energy_j_per_m3 = 0.0_dblprec
      material%diagnostics%small_q_abs_error_j_per_m3 = 0.0_dblprec
      material%diagnostics%small_q_sample_count = n_samples
      material%diagnostics%small_q_max_inverse_m = 0.0_dblprec
      material%diagnostics%small_q_max_abs_error_j_per_m3 = 0.0_dblprec
      material%diagnostics%small_q_max_scaled_error = 0.0_dblprec
      all_pass = .true.

      do sample = 1, n_samples
         normal_norm = sqrt(sum(spiral_normal(:,sample)**2))
         if (abs(normal_norm-1.0_dblprec) > 1.0d-10) then
            diagnostic = 'Every small-q spiral normal must be a Cartesian unit vector'
            material%diagnostics%small_q_energy_validated = .false.
            status = COARSE_MATERIAL_INVALID_INPUT
            return
         end if
         exact_energy = 0.0_dblprec
         do pair = 1, size(input%exchange_source_basis)
            left_channel = input%basis_to_channel(input%exchange_source_basis(pair))
            right_channel = input%basis_to_channel(input%exchange_target_basis(pair))
            if (left_channel < 1 .or. right_channel < 1) cycle
            phase_zero = channel_phase(right_channel,sample) - &
               channel_phase(left_channel,sample)
            phase = phase_zero + dot_product(q_inverse_m(:,sample), &
               input%exchange_displacement_m(:,pair))
            exact_energy = exact_energy + 0.5_dblprec * &
               input%exchange_pair_energy_j(pair) * &
               (cos(phase_zero)-cos(phase)) / input%cell_volume_m3
         end do
         do pair = 1, size(input%dmi_source_basis)
            left_channel = input%basis_to_channel(input%dmi_source_basis(pair))
            right_channel = input%basis_to_channel(input%dmi_target_basis(pair))
            if (left_channel < 1 .or. right_channel < 1) cycle
            phase_zero = channel_phase(right_channel,sample) - &
               channel_phase(left_channel,sample)
            phase = phase_zero + dot_product(q_inverse_m(:,sample), &
               input%dmi_displacement_m(:,pair))
            exact_energy = exact_energy + 0.5_dblprec * &
               dot_product(input%dmi_pair_energy_j(:,pair),spiral_normal(:,sample)) * &
               (sin(phase)-sin(phase_zero)) / input%cell_volume_m3
         end do

         continuum_energy = 0.0_dblprec
         do left_channel = 1, input%n_channels
            do right_channel = 1, input%n_channels
               phase_zero = channel_phase(right_channel,sample) - &
                  channel_phase(left_channel,sample)
               do p = 1, 3
                  do q = 1, 3
                     continuum_energy = continuum_energy + &
                        material%exchange_stiffness(p,q,left_channel,right_channel) * &
                        q_inverse_m(p,sample) * q_inverse_m(q,sample) * cos(phase_zero)
                  end do
                  do k = 1, 3
                     continuum_energy = continuum_energy + &
                        material%spiralization(k,p,left_channel,right_channel) * &
                        spiral_normal(k,sample) * q_inverse_m(p,sample) * cos(phase_zero)
                  end do
               end do
            end do
         end do

         error = abs(exact_energy-continuum_energy)
         scale = max(abs(exact_energy),abs(continuum_energy))
         allowed = absolute_tolerance_j_per_m3 + relative_tolerance*scale
         scaled_error = error / max(allowed,tiny(allowed))
         all_pass = all_pass .and. error <= allowed
         material%diagnostics%small_q_atomistic_energy_j_per_m3(sample) = exact_energy
         material%diagnostics%small_q_continuum_energy_j_per_m3(sample) = continuum_energy
         material%diagnostics%small_q_abs_error_j_per_m3(sample) = error
         material%diagnostics%small_q_max_inverse_m = max( &
            material%diagnostics%small_q_max_inverse_m, &
            sqrt(sum(q_inverse_m(:,sample)**2)))
         material%diagnostics%small_q_max_abs_error_j_per_m3 = max( &
            material%diagnostics%small_q_max_abs_error_j_per_m3,error)
         material%diagnostics%small_q_max_scaled_error = max( &
            material%diagnostics%small_q_max_scaled_error,scaled_error)
      end do

      material%diagnostics%small_q_energy_validated = all_pass
      if (all_pass) then
         if (material%n_channels == 1) then
            material%diagnostics%runtime_gate_reason = ''
         else
            material%diagnostics%runtime_gate_reason = &
               'small-q energies pass, but ferri/AFM acoustic and optical mode extraction remains unvalidated'
         end if
         status = COARSE_MATERIAL_OK
         diagnostic = ''
      else
         material%diagnostics%runtime_gate_reason = &
            'direct small-q atomistic and continuum energies disagree'
         status = COARSE_MATERIAL_VALIDATION_FAILED
         diagnostic = material%diagnostics%runtime_gate_reason
      end if
   end subroutine validate_coarse_material_small_q

   !---------------------------------------------------------------------------
   !> Validate the q=0 acoustic and optical frequencies of the deliberately
   !> limited two-sublattice reference model.  The caller supplies frequencies
   !> independently obtained from an atomistic Hamiltonian or an accepted
   !> small-q dynamical matrix.  For E_loc=L_12 m_1.m_2 the reference has a
   !> zero acoustic frequency and
   !>
   !> |L_12| V_cell/mu_B (gamma_1/M_1 + gamma_2/M_2)
   !>
   !> for the relative (optical) rotation.  This is intentionally not a
   !> general alloy or finite-temperature reduction.
   !---------------------------------------------------------------------------
   subroutine validate_coarse_material_two_sublattice_modes(material, &
         atomistic_acoustic_frequency_per_s, atomistic_optical_frequency_per_s, &
         relative_tolerance, status, diagnostic)
      type(coarse_material_type), intent(inout) :: material
      real(dblprec), intent(in) :: atomistic_acoustic_frequency_per_s
      real(dblprec), intent(in) :: atomistic_optical_frequency_per_s
      real(dblprec), intent(in) :: relative_tolerance
      integer, intent(out) :: status
      character(len=*), intent(out) :: diagnostic

      real(dblprec) :: expected_optical, acoustic_scale, optical_scale

      status = COARSE_MATERIAL_INVALID_INPUT
      diagnostic = ''
      if (.not. material%ready .or. material%n_channels /= 2 .or. &
          .not. allocated(material%local_exchange) .or. &
          .not. allocated(material%channel_moment_mub) .or. &
          .not. allocated(material%channel_gamma) .or. &
          .not. allocated(material%channel_damping)) then
         diagnostic = 'Two-sublattice mode validation requires a ready two-channel material'
         return
      end if
      if (.not. material%diagnostics%channel_dynamics_parameters_validated .or. &
          .not. all(ieee_is_finite(material%local_exchange)) .or. &
          .not. ieee_is_finite(atomistic_acoustic_frequency_per_s) .or. &
          .not. ieee_is_finite(atomistic_optical_frequency_per_s) .or. &
          .not. ieee_is_finite(relative_tolerance) .or. relative_tolerance < 0.0_dblprec) then
         diagnostic = 'Two-sublattice mode inputs must be finite and dynamics parameters accepted'
         return
      end if
      if (any(material%channel_moment_mub <= 0.0_dblprec) .or. &
          any(material%channel_gamma <= 0.0_dblprec) .or. &
          any(material%channel_damping < 0.0_dblprec)) then
         diagnostic = 'Two-sublattice mode validation requires positive moments/gamma and nonnegative damping'
         return
      end if

      expected_optical = abs(material%local_exchange(1,2)) * material%cell_volume_m3 / &
         COARSE_MUB_SI * (material%channel_gamma(1)/material%channel_moment_mub(1) + &
                          material%channel_gamma(2)/material%channel_moment_mub(2))
      acoustic_scale = max(1.0_dblprec,abs(atomistic_optical_frequency_per_s))
      optical_scale = max(1.0_dblprec,abs(expected_optical), &
         abs(atomistic_optical_frequency_per_s))
      material%diagnostics%acoustic_mode_extraction_validated = &
         abs(atomistic_acoustic_frequency_per_s) <= relative_tolerance*acoustic_scale
      material%diagnostics%optical_mode_extraction_validated = &
         abs(atomistic_optical_frequency_per_s-expected_optical) <= &
         relative_tolerance*optical_scale
      material%diagnostics%mode_extraction = COARSE_MODE_EXTRACTION_UNVALIDATED
      if (material%diagnostics%acoustic_mode_extraction_validated .and. &
          material%diagnostics%optical_mode_extraction_validated) then
         material%diagnostics%mode_extraction = COARSE_MODE_EXTRACTION_ACOUSTIC_OPTICAL
         material%diagnostics%runtime_gate_reason = ''
         status = COARSE_MATERIAL_OK
      else
         material%diagnostics%runtime_gate_reason = &
            'two-sublattice atomistic acoustic/optical frequencies disagree with local-exchange dynamics'
         status = COARSE_MATERIAL_VALIDATION_FAILED
         diagnostic = material%diagnostics%runtime_gate_reason
      end if
   end subroutine validate_coarse_material_two_sublattice_modes

   !---------------------------------------------------------------------------
   !> Explicit runtime capability gate.  CG-03 intentionally has no operation
   !> that marks acoustic/optical extraction as validated; that belongs to the
   !> later reviewed multi-channel task.
   !---------------------------------------------------------------------------
   subroutine coarse_material_runtime_status(material, requested_runtime, status, diagnostic)
      type(coarse_material_type), intent(in) :: material
      integer, intent(in) :: requested_runtime
      integer, intent(out) :: status
      character(len=*), intent(out) :: diagnostic

      status = COARSE_MATERIAL_RUNTIME_GATED
      diagnostic = ''
      if (.not. material%ready) then
         diagnostic = 'Coarse material is not ready'
         return
      end if
      select case (requested_runtime)
      case (COARSE_RUNTIME_SINGLE_FM)
         if (material%n_channels /= 1) then
            diagnostic = 'Single-FM runtime requires exactly one dynamical channel'
         else if (.not. material%diagnostics%channel_dynamics_parameters_validated) then
            diagnostic = 'Single-FM runtime requires positive channel gamma and nonnegative damping'
         else if (.not. material%diagnostics%small_q_energy_validated) then
            diagnostic = 'Single-FM runtime requires direct small-q atomistic energy validation'
         else
            status = COARSE_MATERIAL_OK
            diagnostic = ''
         end if
      case (COARSE_RUNTIME_FERRI_AFM)
         if (material%n_channels < 2) then
            diagnostic = 'Ferri/AFM runtime requires at least two explicit dynamical channels'
         else if (.not. material%diagnostics%channel_dynamics_parameters_validated) then
            diagnostic = 'Ferri/AFM runtime requires positive channel gamma and nonnegative damping'
         else if (material%diagnostics%mode_extraction /= &
                  COARSE_MODE_EXTRACTION_ACOUSTIC_OPTICAL .or. &
                  .not. material%diagnostics%acoustic_mode_extraction_validated .or. &
                  .not. material%diagnostics%optical_mode_extraction_validated) then
            diagnostic = 'Ferri/AFM runtime is gated until acoustic and optical mode extraction is validated'
         else
            status = COARSE_MATERIAL_OK
            diagnostic = ''
         end if
      case default
         diagnostic = 'Unknown coarse-material runtime request'
      end select
   end subroutine coarse_material_runtime_status

   pure real(dblprec) function determinant_3x3(matrix) result(determinant)
      real(dblprec), intent(in) :: matrix(3,3)

      determinant = matrix(1,1) * (matrix(2,2)*matrix(3,3)-matrix(2,3)*matrix(3,2)) - &
         matrix(1,2) * (matrix(2,1)*matrix(3,3)-matrix(2,3)*matrix(3,1)) + &
         matrix(1,3) * (matrix(2,1)*matrix(3,2)-matrix(2,2)*matrix(3,1))
   end function determinant_3x3

   !---------------------------------------------------------------------------
   !> @brief
   !> Calculate the exchange stiffness and the spin wave stiffness for a ferromagnet
   !
   !> @author
   !> Jonathan Chico
   !>
   !> @date 10/02/2017 - Jonathan Chico
   !> - Added the calculation of the stiffness tensor, that is now also a matrix
   !>   will be written which can be useful for non cubic systems
   !---------------------------------------------------------------------------
   subroutine ferro_stiffness(NT,NA,N1,N2,N3,hdim,mconf,Natom,nHam,Nchmax,eta_max,eta_min,&
         conf_num,do_ralloy,Natom_full,max_no_neigh,Nch,anumb,atype,nlistsize,&
         atype_ch,asite_ch,achem_ch,nlist,alat,C1,C2,C3,coord,chconc,ammom_inp,ncoup,aham, &
         A_xc,M_sat,Dxc_fit,cell_vol,A_xc_lsq,D_err_fit,Dxc_fit_lsq,J0_matrix,&
         D_xc_stiffness_matrix,A_xc_stiffness_matrix,D_xc_stiffness_matrix_lsq,&
         A_xc_stiffness_matrix_lsq)

      use Constants
      use Math_functions, only : f_wrap_coord_diff

      implicit none

      ! .. Input variables
      integer, intent(in) :: NT  !< Number of types of atoms
      integer, intent(in) :: NA  !< Number of atoms in one cell
      integer, intent(in) :: N1  !< Number of cell repetitions in x direction
      integer, intent(in) :: N2  !< Number of cell repetitions in y direction
      integer, intent(in) :: N3  !< Number of cell repetitions in z direction
      integer, intent(in) :: hdim   !< Number of elements in Hamiltonian element (scalar or vector)
      integer, intent(in) :: mconf  !< LSF ground state configuration
      integer, intent(in) :: Natom  !< Number of atoms in system
      integer, intent(in) :: nHam   !< Number of atoms in Hamiltonian
      integer, intent(in) :: Nchmax !< Max number of chemical components on each site in cell
      integer, intent(in) :: eta_max  !< Number of convergence parameters for the stiffness
      integer, intent(in) :: eta_min  !< Minimum  convergence parameters for the stiffness
      integer, intent(in) :: conf_num !< Number of LSF configurations
      integer, intent(in) :: do_ralloy    !< Random alloy simulation (0/1)
      integer, intent(in) :: Natom_full   !< Number of atoms for full system (=Natom if not dilute)
      integer, intent(in) :: max_no_neigh !< Calculated maximum of neighbours for exchange
      integer, dimension(NA), intent(in) :: Nch !< Number of chemical components on each site in cell
      integer, dimension(Natom), intent(in) :: anumb !< Atom number in cell
      integer, dimension(Natom), intent(in) :: atype !< Type of atom
      integer, dimension(nHam), intent(in) :: nlistsize !< Size of neighbour list for Heisenberg exchange couplings
      integer, dimension(Natom_full), intent(in) :: atype_ch !< Actual site of atom for dilute system
      integer, dimension(Natom_full), intent(in) :: asite_ch !< Actual site of atom for dilute system
      integer, dimension(Natom_full), intent(in) :: achem_ch !< Chemical type of atoms (reduced list)
      integer, dimension(max_no_neigh,Natom), intent(in) :: nlist  !< Neighbour list for Heisenberg exchange couplings

      real(dblprec), intent(in) :: alat !< Lattice parameter
      real(dblprec), dimension(3), intent(in) :: C1 !< First lattice vector
      real(dblprec), dimension(3), intent(in) :: C2 !< Second lattice vector
      real(dblprec), dimension(3), intent(in) :: C3 !< Third lattice vector
      real(dblprec), dimension(3,Natom), intent(in) :: coord !< Coordinates of atoms
      real(dblprec), dimension(NA,Nchmax), intent(in) :: chconc !< Chemical concentration on sites
      real(dblprec), dimension(NA,Nchmax,conf_num), intent(in) :: ammom_inp !< Magnetic moment directions from input (for alloys)
      real(dblprec), dimension(hdim,max_no_neigh,nHam,conf_num), intent(in) :: ncoup !< Heisenberg exchange couplings
      integer, dimension(Natom), intent(in) :: aham !< Hamiltonian look-up table

      ! .. Output variables
      real(dblprec), intent(out) :: A_xc     !< Exchange stiffness constant from rational fit J/m
      real(dblprec), intent(out) :: M_sat    !< Saturation magnetization muB/m^3
      real(dblprec), intent(out) :: Dxc_fit  !< Spin wave stiffnes from rational fit
      real(dblprec), intent(out) :: cell_vol !< Unit cell volume m^3
      real(dblprec), intent(out) :: A_xc_lsq !< Exchange stiffness constant from LSQ fit
      real(dblprec), intent(out) :: D_err_fit   !< Error form spin wave stiffnes rational fit
      real(dblprec), intent(out) :: Dxc_fit_lsq !< Spin wave stiffness LSQ fit
      real(dblprec), dimension(NA,NA), intent(out) :: J0_matrix !< Matrix being used to calculate the Tc-MFA
      real(dblprec), dimension(3,3), intent(out) :: D_xc_stiffness_matrix !< Spin wave stiffness tensor rational fit
      real(dblprec), dimension(3,3), intent(out) :: A_xc_stiffness_matrix !< Exchange stiffness constant rational fit
      real(dblprec), dimension(3,3), intent(out) :: D_xc_stiffness_matrix_lsq !< Spin wave stiffness tensor LSQ fit
      real(dblprec), dimension(3,3), intent(out) :: A_xc_stiffness_matrix_lsq !< Exchange stiffness constant LSQ fit

      ! .. Local variables
      integer :: ii, jj,i_stat,i_all
      integer :: lwork, info, alwork
      integer :: iatom, jatom, iham
      integer :: katom, I1, I2, I3, countstart
      integer :: i, k, eta, ich, eta_redu

      real(dblprec) :: total_mom
      real(dblprec) :: fcinv, jij, jijsign
      real(dblprec) :: rij2, rij, rfit, drfit, stiff_par

      real(dblprec), dimension(3) :: rcoord
      real(dblprec), dimension(0:eta_max) :: temp_x, eig_val, eig_val_temp
      real(dblprec), dimension(eta_max-(eta_min-1),3) :: lmatrix
      real(dblprec), dimension(eta_max-(eta_min-1)) :: dvector
      real(dblprec), dimension(0:eta_max,3,3) :: eig_val_mat

      real(dblprec), dimension(:), allocatable :: wres
      real(dblprec), dimension(:), allocatable :: awork
      real(dblprec), dimension(:,:,:), allocatable :: etemp
      real(dblprec), dimension(:,:,:), allocatable :: stiff_matrix !< Matrix being used to calculate the ferromagnetic exchange stiffness
      real(dblprec), dimension(:,:,:,:,:), allocatable :: D_matrix !< Matrix used to calculate the FM exchange stiffness tensor

      real(dblprec), dimension(:), allocatable :: work
      real(dblprec), dimension(:), allocatable :: rwork
      real(dblprec), dimension(:), allocatable :: iwres
      real(dblprec), dimension(:), allocatable :: iwres_mat
      real(dblprec), dimension(:), allocatable :: rwres
      real(dblprec), dimension(:), allocatable :: rwres_mat
      real(dblprec), dimension(:), allocatable :: cwres
      real(dblprec), dimension(:), allocatable :: cwres_mat
      real(dblprec), dimension(:,:),allocatable :: ctemp
      real(dblprec), dimension(:,:),allocatable :: A_inplace
      real(dblprec), dimension(:,:),allocatable :: A_mat_inplace

      lwork=16*na
      alwork=((eta_max-(eta_min-1))*16)

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Allocate work arrays
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      allocate(cwres(NA),stat=i_stat)
      call memocc(i_stat,product(shape(cwres))*kind(cwres),'cwres','ferro_stiffness')
      allocate(rwres(NA),stat=i_stat)
      call memocc(i_stat,product(shape(rwres))*kind(rwres),'rwres','ferro_stiffness')
      allocate(iwres(NA),stat=i_stat)
      call memocc(i_stat,product(shape(iwres))*kind(iwres),'iwres','ferro_stiffness')
      allocate(work(lwork),stat=i_stat)
      call memocc(i_stat,product(shape(work))*kind(work),'work','ferro_stiffness')
      allocate(cwres_mat(NA),stat=i_stat)
      call memocc(i_stat,product(shape(cwres_mat))*kind(cwres_mat),'cwres_mat','ferro_stiffness')
      allocate(rwres_mat(NA),stat=i_stat)
      call memocc(i_stat,product(shape(rwres_mat))*kind(rwres_mat),'rwres_mat','ferro_stiffness')
      allocate(iwres_mat(NA),stat=i_stat)
      call memocc(i_stat,product(shape(iwres_mat))*kind(iwres_mat),'iwres_mat','ferro_stiffness')
      allocate(ctemp(NA,NA),stat=i_stat)
      call memocc(i_stat,product(shape(ctemp))*kind(ctemp),'ctemp','ferro_stiffness')
      allocate(rwork(lwork),stat=i_stat)
      call memocc(i_stat,product(shape(rwork))*kind(rwork),'rwork','ferro_stiffness')
      allocate(wres(eta_max),stat=i_stat)
      call memocc(i_stat,product(shape(wres))*kind(wres),'wres','ferro_stiffness')
      allocate(awork(alwork),stat=i_stat)
      call memocc(i_stat,product(shape(awork))*kind(awork),'awork','ferro_stiffness')
      allocate(A_inplace(NA,NA),stat=i_stat)
      call memocc(i_stat,product(shape(A_inplace))*kind(A_inplace),'A_inplace','ferro_stiffness')
      allocate(A_mat_inplace(NA,NA),stat=i_stat)
      call memocc(i_stat,product(shape(A_mat_inplace))*kind(A_mat_inplace),'A_mat_inplace','ferro_stiffness')
      allocate(etemp(NA,NA,0:eta_max),stat=i_stat)
      call memocc(i_stat,product(shape(etemp))*kind(etemp),'etemp','ferro_stiffness')
      allocate(stiff_matrix(NA,NA,0:eta_max),stat=i_stat)
      call memocc(i_stat,product(shape(stiff_matrix))*kind(stiff_matrix),'stiff_matrix','ferro_stiffness')
      allocate(D_matrix(3,3,NA,NA,0:eta_max),stat=i_stat)
      call memocc(i_stat,product(shape(D_matrix))*kind(D_matrix),'D_matrix','ferro_stiffness')

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Initialization of variables
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      work=0.0_dblprec
      wres=0.0_dblprec
      etemp=0.0_dblprec
      ctemp=0.0_dblprec
      rwork=0.0_dblprec
      cwres=0.0_dblprec
      iwres=0.0_dblprec
      rwres=0.0_dblprec
      temp_x=0.0_dblprec
      eig_val=0.0_dblprec
      lmatrix=0.0_dblprec
      dvector=0.0_dblprec
      cwres_mat=0.0_dblprec
      iwres_mat=0.0_dblprec
      rwres_mat=0.0_dblprec
      total_mom=0.0_dblprec
      J0_matrix=0.0_dblprec
      A_inplace=0.0_dblprec
      eig_val_mat=0.0_dblprec
      eig_val_temp=0.0_dblprec
      stiff_matrix=0.0_dblprec
      A_mat_inplace=0.0_dblprec
      D_matrix=0.0_dblprec
      D_xc_stiffness_matrix=0.0_dblprec
      A_xc_stiffness_matrix=0.0_dblprec
      D_xc_stiffness_matrix_lsq=0.0_dblprec
      A_xc_stiffness_matrix_lsq=0.0_dblprec
      D_err_fit=0.0_dblprec

      fcinv=mub/mry
      I1 = N1/2
      I2 = N2/2
      I3 = N3/2
      countstart = 0+I1*NA+I2*N1*NA+I3*N2*N1*NA
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Calculate the volume of the cell in meters
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      cell_vol=(C1(1)*C2(2)*C3(3)-C1(1)*C2(3)*C3(2)+&
         C1(2)*C2(3)*C3(1)-C1(2)*C2(1)*C3(3)+&
         C1(3)*C2(1)*C3(2)-C1(3)*C2(2)*C3(1))*alat**3
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Calculating the total moment of the cell
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      if (do_ralloy==0) then
         do i=1,NA
            total_mom=total_mom+ammom_inp(i,1,mconf)
         enddo
         ! In case of random alloy calculate the weighted average
      else if (do_ralloy==1) then
         do i=1,NA
            do ich=1,Nch(i)
               total_mom=total_mom+ammom_inp(i,ich,mconf)*chconc(i,ich)
            enddo
         enddo
      endif

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Calculating the saturation magnetization
      M_sat=total_mom/cell_vol
      write(420,'(f14.6,g20.8)') total_mom,M_sat
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Calculate the stiffness for checmically pure systems
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      if (do_ralloy==0) then

         ! Need to create a matrix which includes inter and intra sublattice interactions
         ! Now must loop over the convergency factor to make sure that the sum is well defined
         do eta=0, eta_max

            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            ! Loop over atoms in the unit cell and sum up the exchange interactions
            do i=1, NA
               iatom=i+countstart
               iham=aham(iatom)
               ! Loop over the neighbors for each atom in the unit cell
               do jatom=1, nlistsize(iham)
                  ! Neighbouring atom
                  katom=nlist(jatom,iatom)

                  ! Distace vector between the atoms
                  call f_wrap_coord_diff(Natom,coord,katom,iatom,rcoord)

                  ! Distance vector between neighbouring atoms
                  rij2=rcoord(1)**2+rcoord(2)**2+rcoord(3)**2
                  ! Distance between neighbouring atoms
                  rij=sqrt(rij2)

                  ! Calculating the "real" Jijs in mRyd
                  jij=ncoup(1,jatom,iham,mconf)*fcinv*ammom_inp(anumb(katom),1,mconf)*ammom_inp(anumb(iatom),1,mconf)*0.5_dblprec
                  ! Which is the sign between the magnetic moments (To take care of AFM interactions)
                  jijsign=sign(ammom_inp(anumb(katom),1,mconf),ammom_inp(anumb(iatom),1,mconf))/abs(ammom_inp(anumb(katom),1,mconf))
                  ! Calculate the proportionality parameter
                  stiff_par=(2.0_dblprec/3.0_dblprec)*jij*jijsign*rij2*(alat**2)/(sqrt(abs(ammom_inp(anumb(katom),1,mconf))*abs(ammom_inp(anumb(iatom),1,mconf))))
                  ! The actual stiffness matrix
                  stiff_matrix(i,anumb(katom),eta)=stiff_matrix(i,anumb(katom),eta)+stiff_par*exp(-0.10_dblprec*eta*rij)

                  ! Calculate the stiffness matrix
                  do ii=1,3
                     do jj=1,3
                        ! Create a parameter for the calculation of the stiffness
                        stiff_par=2.0_dblprec*jij*jijsign*rcoord(ii)*rcoord(jj)*(alat**2)/&
                           (sqrt(abs(ammom_inp(anumb(katom),1,mconf))*abs(ammom_inp(anumb(iatom),1,mconf))))
                        ! Save the parameter for the stiffness matrix
                        D_matrix(ii,jj,i,anumb(katom),eta)=D_matrix(ii,jj,i,anumb(katom),eta)+&
                           stiff_par*exp(-0.10_dblprec*eta*sqrt(abs(rcoord(ii)*rcoord(jj))))
                     enddo
                  enddo
                  J0_matrix(i,anumb(katom))=J0_matrix(i,anumb(katom))+(2.0_dblprec/3.0_dblprec)*jij*jijsign
               enddo

            enddo

            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            ! Calculating the eigenvalues for the scalar stiffness
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            ! The matrix is transformed to the more familiar units of meVÅ^2
            A_inplace(1:NA,1:NA)=stiff_matrix(1:NA,1:NA,eta)*1d20*ry_ev

            ! The eigenvalues for the spin wave stiffness are calculated using LAPACK
            call dgeev('N','N',NA, A_inplace, NA, rwres, iwres, ctemp, NA, etemp(1,1,1), NA, WORK, LWORK, INFO)
            if(info.ne.0) then
               print '(2x,a,i4)', 'Problem in zgeev 1:',info
            end if
            if(eta==0) then
               write(2420,'(3f14.6)') A_inplace/ry_ev
               write(1420,'(3f14.6)') maxval((rwres))/ry_ev
               write(420,'(3f14.6)') maxval((rwres))
            end if

            ! Temporal x-axis for the fitting to a polynomial
            temp_x(eta)=0.10_dblprec*eta
            ! Finding the maximum eigenvalue of the exchange stiffness matrix
            eig_val(eta)=maxval((rwres))
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            ! End of calculation of eignevalues for the scalar stiffness
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            ! Calculating the eigenvalues for the stiffness tensor
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            do ii=1,3
               do jj=1,3
                  ! Transform to meVA^2
                  A_mat_inplace(1:NA,1:NA)=D_matrix(ii,jj,1:NA,1:NA,eta)*1d20*ry_ev
                  ! Calculate eigenvalues for each component of the matrix
                  call dgeev('N','N',NA, A_mat_inplace, NA, rwres_mat, iwres_mat, ctemp, NA, etemp(1,1,1), NA, WORK, LWORK, INFO)
                  if(info.ne.0) then
                     print '(2x,a,i4)', 'Problem in zgeev 2:',info
                  end if
                  eig_val_mat(eta,ii,jj)=maxval((rwres_mat))
               enddo
            enddo
            !  AB hack
            if(eta==0) then
               write(2420,'(3f14.6)') D_matrix(:,:,1,1,eta)*1d20
               write(1420,'(3f14.6)') eig_val_mat(eta,:,:)/ry_ev
               write(420,'(3f14.6)') eig_val_mat(eta,:,:)
            end if
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            ! End of calculation of eignevalues for the stiffness tensor
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         enddo

         ! Defining the reduced eta, common for scalar and tensor
         eta_redu=eta_max-(eta_min-1)
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         ! Calculating the scalar stiffness with rational polynomials
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

         ! A polynomial fit is performed to obtain the spin wave stiffness
         call ratint(temp_x(eta_min:eta_max),eig_val(eta_min:eta_max),eta_redu,0.0_dblprec,rfit,drfit)
         ! Calculate the exchange stiffness from micromagnetics with rational
         ! polynomials
         A_xc=rfit*M_sat*1e-20/(1000*Joule_ev*4.0_dblprec)
         ! Storing the spin wave stiffness scalar
         Dxc_fit=rfit
         D_err_fit=drfit

         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         ! End of calculating the scalar stiffness with rational polynomials
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         ! Calculating the stiffness tensor with rational polynomials
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         do ii=1,3
            do jj=1,3
               eig_val_temp(1:eta_max)=eig_val_mat(1:eta_max,ii,jj)
               ! A polynomial fit is performed to obtain the spin wave stiffness
               call ratint(temp_x(eta_min:eta_max),eig_val_temp(eta_min:eta_max),eta_redu,0.0_dblprec,rfit,drfit)
               ! Calculate the exchange stiffness from micromagnetics with rational
               ! polynomials
               A_xc_stiffness_matrix(ii,jj)=rfit*1e-20*M_sat/(1000*Joule_ev*4.0_dblprec)
               ! Saving the spin wave stiffness matrix
               D_xc_stiffness_matrix(ii,jj)=rfit
            enddo
         enddo
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         ! End of calculating the stiffness tensor with rational polynomials
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         ! Calculating the scalar exchange stiffness with LSQ
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         ! Creating a vector for the eta parameter
         do eta=eta_min,eta_max
            k=eta-(eta_min-1)
            lmatrix(k,1)=1._dblprec ; lmatrix(k,2)=0.1_dblprec*eta ; lmatrix(k,3)=(0.1_dblprec*eta)**2; dvector(k)=eig_val(eta)
         enddo

         call dgels('N',eta_max-(eta_min-1),3,1,lmatrix,eta_max-(eta_min-1),dvector,eta_max-(eta_min-1),awork,alwork,info)
         if(info.ne.0) then
            print '(2x,a,i4)', 'Problem in dgels:',info
         end if
         ! Calculate the Exchange stiffness from micromagnetics
         A_xc_lsq=dvector(1)*1e-20*M_sat/(1000*Joule_ev*4.0_dblprec)

         Dxc_fit_lsq=dvector(1)

         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         ! End of calculation of the scalar exchange stiffness with LSQ
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         ! Calculating the exchange stiffness tensor with LSQ
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

         do ii=1, 3
            do jj=1,3
               ! Creating a vector for the eta parameter
               do eta=eta_min,eta_max
                  k=eta-(eta_min-1)
                  lmatrix(k,1)=1._dblprec ; lmatrix(k,2)=0.1_dblprec*eta ; lmatrix(k,3)=(0.1_dblprec*eta)**2; dvector(k)=eig_val_mat(eta,ii,jj)
               enddo

               call dgels('N',eta_max-(eta_min-1),3,1,lmatrix,eta_max-(eta_min-1),dvector,eta_max-(eta_min-1),awork,alwork,info)
               if(info.ne.0) then
                  print '(2x,a,i4)', 'Problem in dgels:',info
               end if
               ! Calculate the exchange stiffness from micromagnetics
               A_xc_stiffness_matrix_lsq(ii,jj)=dvector(1)*1e-20*M_sat/(1000*Joule_ev*4.0_dblprec)

               D_xc_stiffness_matrix_lsq(ii,jj)=dvector(1)

            enddo
         enddo

         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         ! End of calculation of the exchange stiffness tensor with LSQ
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         ! Calculation of the MF Tc
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

         J0_matrix=J0_matrix/eta_max
         A_inplace(1:NA,1:NA)=(J0_matrix*ry_ev*1e-3)/k_bolt_ev
         ! The eigenvalues for the spin wave stiffness are calculated using LAPACK
         call dgeev('N','N',NA, A_inplace, NA, rwres, iwres, ctemp, NA, etemp(1,1,1), NA, WORK, LWORK, INFO)
         if(info.ne.0) then
            print '(2x,a,i4)', 'Problem in zgeev 3:',info
         end if
         eig_val(1)=maxval((rwres))
         write(*,1009) 'Tc-MFA from stiffness :',eig_val(1),'K'

         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         ! End of calculation of the MF Tc
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         ! Calculate the stiffness for the random alloy
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      else if (do_ralloy==1) then
         ! Loop over the convergency factor to make sure that the sum is well defined

         do eta=1, eta_max

            ! Loop over atoms the atoms in the unit cell
            do i=1,NA
               iatom=i+countstart
               iham=aham(iatom)
               ! Loop over the neighbors for each atom in the unit cell
               do jatom=1, nlistsize(iham)

                  ! Neighbouring atom
                  katom=nlist(jatom,iatom)

                  ! Distace vector between the atoms
                  call f_wrap_coord_diff(Natom,coord,katom,iatom,rcoord)

                  ! Distance vector between neighbouring atoms
                  rij2=rcoord(1)**2+rcoord(2)**2+rcoord(3)**2
                  ! Distance between neighbouring atoms
                  rij=sqrt(rij2)

                  ! Calculating the "real" Jijs in mRyd
                  jij=ncoup(1,jatom,iham,mconf)*fcinv*ammom_inp(asite_ch(katom),achem_ch(katom),mconf)*&
                     ammom_inp(asite_ch(iatom),achem_ch(iatom),mconf)*0.5_dblprec
                  ! Which is the sign between the magnetic moments (To take care of AFM interactions)
                  jijsign=sign(ammom_inp(asite_ch(katom),achem_ch(katom),mconf),ammom_inp(asite_ch(iatom),achem_ch(iatom),mconf))/&
                     abs(ammom_inp(asite_ch(katom),achem_ch(katom),mconf))

                  ! Parameter for the stiffness
                  stiff_par=(2.0_dblprec/3.0_dblprec)*jij*jijsign*rij2*(alat**2)/(sqrt(abs(ammom_inp(asite_ch(katom),achem_ch(katom),mconf))*&
                     abs(ammom_inp(asite_ch(iatom),achem_ch(iatom),mconf))))
                  ! The actual stiffness matrix
                  stiff_matrix(i,asite_ch(katom),eta)=stiff_matrix(i,asite_ch(katom),eta)+&
                     stiff_par*exp(-0.10_dblprec*eta*rij)

                  ! Calculate the stiffness matrix
                  do ii=1,3
                     do jj=1,3
                        D_matrix(ii,jj,i,anumb(katom),eta)=D_matrix(ii,jj,i,asite_ch(katom),eta)+&
                           2.0_dblprec*jij*jijsign*rcoord(ii)*rcoord(jj)*(alat**2)*exp(-0.10_dblprec*eta*rij)/&
                           (sqrt(abs(ammom_inp(asite_ch(katom),achem_ch(katom),mconf))*abs(ammom_inp(asite_ch(iatom),achem_ch(iatom),mconf))))
                     enddo
                  enddo

                  ! Filling the J0 matrix
                  J0_matrix(i,asite_ch(katom))=J0_matrix(i,asite_ch(katom))+(2.0_dblprec/3.0_dblprec)*jij*jijsign
               enddo

            enddo

            ! Calculation of the eigenvalues for the stiffness matrix

            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            ! Calculating the eigenvalues for the scalar stiffness
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            ! The matrix is transformed to the more familiar units of meVÅ^2
            A_inplace(1:NA,1:NA)=stiff_matrix(1:NA,1:NA,eta)*1d20*ry_ev

            ! The eigenvalues for the spin wave stiffness are calculated using LAPACK
            call dgeev('N','N',NA, A_inplace, NA, rwres, iwres, ctemp, NA, etemp(1,1,1), NA, WORK, LWORK, INFO)
            if(info.ne.0) then
               print '(2x,a,i4)', 'Problem in zgeev 4:',info
            end if

            ! Temporal x-axis for the fitting to a polynomial
            temp_x(eta)=0.10_dblprec*eta
            ! Finding the maximum eigenvalue of the exchange stiffness matrix
            eig_val(eta)=maxval((rwres))
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            ! End of calculation of eignevalues for the scalar stiffness
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            ! Calculating the eigenvalues for the stiffness tensor
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            do ii=1,3
               do jj=1,3
                  ! Transform to meVA^2
                  A_mat_inplace(1:NA,1:NA)=D_matrix(ii,jj,1:NA,1:NA,eta)*1d20*ry_ev
                  ! Calculate eigenvalues for each component of the matrix
                  call dgeev('N','N',NA, A_mat_inplace, NA, rwres_mat, iwres_mat, ctemp, NA, etemp(1,1,1), NA, WORK, LWORK, INFO)
                  eig_val_mat(eta,ii,jj)=maxval((rwres_mat))
               enddo
            enddo
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            ! End of calculation of eignevalues for the stiffness tensor
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         enddo

         ! Defining the reduced eta, common for scalar and tensor
         eta_redu=eta_max-(eta_min-1)
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         ! Calculating the scalar stiffness with rational polynomials
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         ! A polynomial fit is performed to obtain the spin wave stiffness
         call ratint(temp_x(eta_min:eta_max),eig_val(eta_min:eta_max),eta_redu,0.0_dblprec,rfit,drfit)
         ! Calculate the exchange stiffness from micromagnetics with rational
         ! polynomials
         A_xc=rfit*M_sat*1e-20/(1000*Joule_ev*4.0_dblprec)
         ! Storing the spin wave stiffness scalar
         Dxc_fit=rfit
         D_err_fit=drfit

         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         ! End of calculating the scalar stiffness with rational polynomials
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         ! Calculating the stiffness tensor with rational polynomials
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         do ii=1,3
            do jj=1,3
               eig_val_temp(1:eta_max)=eig_val_mat(1:eta_max,ii,jj)
               ! A polynomial fit is performed to obtain the spin wave stiffness
               call ratint(temp_x(eta_min:eta_max),eig_val_temp(eta_min:eta_max),eta_redu,0.0_dblprec,rfit,drfit)
               ! Calculate the exchange stiffness from micromagnetics with rational
               ! polynomials
               A_xc_stiffness_matrix(ii,jj)=rfit*1e-20*M_sat/(1000*Joule_ev*4.0_dblprec)
               ! Saving the spin wave stiffness matrix
               D_xc_stiffness_matrix(ii,jj)=rfit
            enddo
         enddo
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         ! End of calculating the stiffness tensor with rational polynomials
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         ! Calculating the scalar exchange stiffness with LSQ
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         ! Creating a vector for the eta parameter
         do eta=eta_min,eta_max
            k=eta-(eta_min-1)
            lmatrix(k,1)=1._dblprec ; lmatrix(k,2)=0.1_dblprec*eta ; lmatrix(k,3)=(0.1_dblprec*eta)**2; dvector(k)=eig_val(eta)
         enddo

         call dgels('N',eta_max-(eta_min-1),3,1,lmatrix,eta_max-(eta_min-1),dvector,eta_max-(eta_min-1),awork,alwork,info)
         if(info.ne.0) then
            print '(2x,a,i4)', 'Problem in dgels:',info
         end if
         ! Calculate the Exchange stiffness from micromagnetics
         A_xc_lsq=dvector(1)*1e-20*M_sat/(1000*Joule_ev*4)

         Dxc_fit_lsq=dvector(1)

         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         ! End of calculation of the scalar exchange stiffness with LSQ
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         ! Calculating the exchange stiffness tensor with LSQ
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

         do ii=1, 3
            do jj=1,3
               ! Creating a vector for the eta parameter
               do eta=eta_min,eta_max
                  k=eta-(eta_min-1)
                  lmatrix(k,1)=1._dblprec ; lmatrix(k,2)=0.1_dblprec*eta ; lmatrix(k,3)=(0.1_dblprec*eta)**2; dvector(k)=eig_val_mat(eta,ii,jj)
               enddo

               call dgels('N',eta_max-(eta_min-1),3,1,lmatrix,eta_max-(eta_min-1),dvector,eta_max-(eta_min-1),awork,alwork,info)
               if(info.ne.0) then
                  print '(2x,a,i4)', 'Problem in dgels:',info
               end if
               ! Calculate the Exchange stiffness from micromagnetics
               A_xc_stiffness_matrix_lsq(ii,jj)=dvector(1)*1e-20*M_sat/(1000*Joule_ev*4)

               D_xc_stiffness_matrix_lsq(ii,jj)=dvector(1)

            enddo
         enddo

         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         ! End of calculation of the exchange stiffness tensor with LSQ
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         ! Calculation of the MF Tc
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

         J0_matrix=J0_matrix/eta_max
         A_inplace(1:NA,1:NA)=(J0_matrix*ry_ev*1e-3)/k_bolt_ev
         ! The eigenvalues for the spin wave stiffness are calculated using LAPACK
         call dgeev('N','N',NA, A_inplace, NA, rwres, iwres, ctemp, NA, etemp(1,1,1), NA, WORK, LWORK, INFO)
         if(info.ne.0) then
            print '(2x,a,i4)', 'Problem in zgeev 5:',info
         end if
         eig_val(1)=maxval((rwres))
!         write(*,1009) 'Tc-MFA from stiffness :',eig_val(1),'K'

         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         ! End of calculation of the MF Tc
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      endif

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! End of stiffness calculations
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      ! The issue for the chemical dissorder is that there might not be different chemical types in the same cell
      i_all=-product(shape(work))*kind(work)
      deallocate(work,stat=i_stat)
      call memocc(i_stat,i_all,'work','ferro_stiffness')
      i_all=-product(shape(wres))*kind(wres)
      deallocate(wres,stat=i_stat)
      call memocc(i_stat,i_all,'wres','ferro_stiffness')
      i_all=-product(shape(cwres))*kind(cwres)
      deallocate(cwres,stat=i_stat)
      call memocc(i_stat,i_all,'cwres','ferro_stiffness')
      i_all=-product(shape(ctemp))*kind(ctemp)
      deallocate(ctemp,stat=i_stat)
      call memocc(i_stat,i_all,'ctemp','ferro_stiffness')
      i_all=-product(shape(rwork))*kind(rwork)
      deallocate(rwork,stat=i_stat)
      call memocc(i_stat,i_all,'rwork','ferro_stiffness')
      i_all=-product(shape(awork))*kind(awork)
      deallocate(awork,stat=i_stat)
      call memocc(i_stat,i_all,'awork','ferro_stiffness')
      i_all=-product(shape(etemp))*kind(etemp)
      deallocate(etemp,stat=i_stat)
      call memocc(i_stat,i_all,'etemp','ferro_stiffness')
      i_all=-product(shape(D_matrix))*kind(D_matrix)
      deallocate(D_matrix,stat=i_stat)
      call memocc(i_stat,i_all,'D_matrix','ferro_stiffness')
      i_all=-product(shape(A_inplace))*kind(A_inplace)
      deallocate(A_inplace,stat=i_stat)
      call memocc(i_stat,i_all,'A_inplace','ferro_stiffness')
      i_all=-product(shape(cwres_mat))*kind(cwres_mat)
      deallocate(cwres_mat,stat=i_stat)
      call memocc(i_stat,i_all,'cwres_mat','ferro_stiffness')
      i_all=-product(shape(stiff_matrix))*kind(stiff_matrix)
      deallocate(stiff_matrix,stat=i_stat)
      call memocc(i_stat,i_all,'stiff_matrix','ferro_stiffness')
      i_all=-product(shape(A_mat_inplace))*kind(A_mat_inplace)
      deallocate(A_mat_inplace,stat=i_stat)
      call memocc(i_stat,i_all,'A_mat_inplace','ferro_stiffness')

      1009 format (2x,a,f10.1,2x,a)
   end subroutine ferro_stiffness

   !---------------------------------------------------------------------------
   !> @brief
   !> Calculation of the tensorial DMI stiffness
   !
   !> @author
   !> Anders Bergman
   !> @date 10/02/2017 - Jonathan Chico
   !> - Routine originally written by Anders Bergman, implemented in ASD by
   !> Jonathan Chico
   !> - Modification to make compatible with new printing routine
   !---------------------------------------------------------------------------
   subroutine DMI_stiffness(NT,NA,N1,N2,N3,Natom,nHam, Nchmax,eta_max,eta_min,max_no_dmneigh,anumb,&
         dmlistsize,dmlist,alat,coord,ammom_inp,dm_vect,DM0_mat,DM0_mat_lsq,aham)

      use Constants
      use Math_functions, only : f_wrap_coord_diff

      implicit none

      ! .. Input variables
      integer, intent(in) :: NT  !< Number of types of atoms
      integer, intent(in) :: NA  !< Number of atoms in one cell
      integer, intent(in) :: N1  !< Number of cell repetitions in x direction
      integer, intent(in) :: N2  !< Number of cell repetitions in y direction
      integer, intent(in) :: N3  !< Number of cell repetitions in z direction
      integer, intent(in) :: Natom  !< Number of atoms in system
      integer, intent(in) :: nHam   !< Number of atoms in Hamiltonian
      integer, intent(in) :: Nchmax !< Max number of chemical components on each site in cell
      integer, intent(in) :: eta_max  !< Number of convergence parameters for the stiffness
      integer, intent(in) :: eta_min  !< Minimum  convergence parameters for the stiffness
      integer, intent(in) :: max_no_dmneigh !< Calculated maximum of neighbours for exchange
      integer, dimension(Natom), intent(in) :: anumb !< Atom number in cell
      integer, dimension(nHam), intent(in) :: dmlistsize !< Size of neighbour list for DM
      integer, dimension(max_no_dmneigh,Natom), intent(in) :: dmlist !< Neighbour list for Heisenberg exchange couplings

      real(dblprec), intent(in) :: alat !< Lattice parameter
      real(dblprec), dimension(3,Natom), intent(in) :: coord !< Coordinates of atoms
      real(dblprec), dimension(NA,Nchmax), intent(in) :: ammom_inp !< Magnetic moment directions from input (for alloys)
      real(dblprec), dimension(3,max_no_dmneigh,nHam), intent(in) :: dm_vect !< Heisenberg exchange couplings
      integer, dimension(Natom), intent(in) :: aham  !< Hamiltonian look-up table

      ! .. Output variables
      real(dblprec), dimension(3,3), intent(out) :: DM0_mat !< DMI spiralization matrix [meVA]
      real(dblprec), dimension(3,3), intent(out) :: DM0_mat_lsq !< DMI spiralization matrix [meVA]

      ! .. Local variables
      integer :: info,ii,jj,i_stat,i_all
      integer :: i, j, k, alwork, eta
      integer :: iatom, jatom, iham
      integer :: katom, I1, I2, I3, countstart
      real(dblprec) :: rij2, rij
      real(dblprec) :: fcinv, dmsign
      real(dblprec) :: dm_mag_par
      real(dblprec), dimension(3) :: rcoord, DM_xc, dm_stiff,dij
      real(dblprec), dimension(:),allocatable :: eigenvals
      real(dblprec), dimension(:,:,:,:,:), allocatable :: dm_mat !< Matrix being used to calculate the DM stiffness
      complex(dblprec), dimension(:,:),allocatable :: eigenvecs
      real(dblprec), dimension(eta_max-(eta_min-1),3) :: lmatrix
      real(dblprec), dimension(eta_max-(eta_min-1)) :: dvector
      real(dblprec), dimension(:), allocatable :: awork
      real(dblprec), dimension(:,:,:), allocatable :: DM0_mat_eta

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Allocate working arrays
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      alwork=((eta_max-(eta_min-1))*16)
      allocate(eigenvals(NA),stat=i_stat)
      call memocc(i_stat,product(shape(eigenvals))*kind(eigenvals),'eigenvals','DMI_stiffness')
      allocate(eigenvecs(NA,NA),stat=i_stat)
      call memocc(i_stat,product(shape(eigenvecs))*kind(eigenvecs),'eigenvecs','DMI_stiffness')
      allocate(dm_mat(3,3,NA,NA,0:eta_max),stat=i_stat)
      call memocc(i_stat,product(shape(dm_mat))*kind(eigenvecs),'dm_mat','DMI_stiffness')
      allocate(DM0_mat_eta(3,3,0:eta_max),stat=i_stat)
      call memocc(i_stat,product(shape(DM0_mat_eta))*kind(DM0_mat_eta),'DM0_mat_eta','DMI_stiffness')
      allocate(awork(alwork),stat=i_stat)
      call memocc(i_stat,product(shape(awork))*kind(awork),'awork','DMI_stiffness')
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Initialize variables
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      DM_xc=0.0_dblprec
      rcoord=0.0_dblprec
      dm_mat=0.0_dblprec
      DM0_mat=0.0_dblprec
      DM0_mat_eta=0.0_dblprec
      dm_stiff=0.0_dblprec
      eigenvals=0.0_dblprec
      lmatrix=0.0_dblprec
      dvector=0.0_dblprec

      I1 = N1/2
      I2 = N2/2
      I3 = N3/2
      fcinv=mub/mry
      countstart = 0+I1*NA+I2*N1*NA+I3*N2*N1*NA
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      ! This routine might be redundant
      if (anumb(countstart+1)/=1) then ! could be mod(countstart,NA) ==0 aswell
         do i = 1,NA
            if (anumb(countstart+i+1)==1) then
               countstart = countstart+i
               exit
            else if (anumb(countstart-i+1)==1) then
               countstart = countstart-i
               exit
            end if
         end do
      end if


      do eta=0,eta_max
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         ! Loop over atoms in the unit cell and sum up the DMI vectors
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         do i=1, NA
            iatom=i+countstart
            iham=aham(iatom)
            ! Loop over the neighbors for each atom in the unit cell
            do jatom=1, dmlistsize(iham)

               ! Neighbouring atom
               katom=dmlist(jatom,iatom)
               ! Distance VECTOR betwwen neighbouring atoms
               call f_wrap_coord_diff(Natom,coord,katom,iatom,rcoord)

               ! Distance between the atoms
               rij2=rcoord(1)**2+rcoord(2)**2+rcoord(3)**2
               rij=sqrt(rij2)

               ! Calculating the "real" DMI in mRyd
               dij(1)=dm_vect(1,jatom,iham)*fcinv*ammom_inp(anumb(katom),1)*ammom_inp(anumb(iatom),1)*0.5_dblprec
               dij(2)=dm_vect(2,jatom,iham)*fcinv*ammom_inp(anumb(katom),1)*ammom_inp(anumb(iatom),1)*0.5_dblprec
               dij(3)=dm_vect(3,jatom,iham)*fcinv*ammom_inp(anumb(katom),1)*ammom_inp(anumb(iatom),1)*0.5_dblprec
               ! Which is the sign between the magnetic moments (To take care of AFM interactions)
               dmsign=sign(ammom_inp(anumb(katom),1),ammom_inp(anumb(iatom),1))/abs(ammom_inp(anumb(katom),1))

               ! Variable to weight in the moments and lattice parameter
               dm_mag_par=alat/sqrt(abs(ammom_inp(anumb(katom),1))*abs(ammom_inp(anumb(iatom),1)))

               ! The actual DM stiffness matrix in here it will be looped over spin
               ! and lattice space
               do ii=1,3
                  do jj=1,3
                     ! Notice that this matrix mixes both lattice and spin degrees of
                     ! freedom as the dij acts on the spins and rij is a vector in
                     ! real space
                     !dm_mat(ii,jj,i,anumb(katom))=dm_mat(ii,jj,i,anumb(katom))+dmsign*rcoord(ii)*dij(jj)*&
                     !                             dm_mag_par
                     dm_mat(ii,jj,i,anumb(katom),eta)=dm_mat(ii,jj,i,anumb(katom),eta)+dmsign*rcoord(ii)*dij(jj)*&
                        dm_mag_par*exp(-0.10_dblprec*eta*rij)
                  enddo
               enddo

            enddo

         enddo
      end do

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! End of calculation of DM spiralization matrix
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      !!! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!! ! Loop to calculate the eigenvalues of the DM spiralization matrix
      !!! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      DM0_mat_eta=0.0_dblprec
      do eta=0,eta_max
         do i=1,3
            do j=1,3
               do ii=1,na
                  do jj=1,na
                     DM0_mat_eta(i,j,eta)=DM0_mat_eta(i,j,eta)+dm_mat(i,j,ii,jj,eta)*1d10*ry_ev
                  end do
               end do
            enddo
         enddo
      end do
      DM0_mat=DM0_mat_eta(:,:,0)
      !!! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!! ! End of calculation of the DM spiralzation eigenvalues
      !!! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! LSQ over eta
      do i=1,3
         do j=1,3
            ! Creating a vector for the eta parameter
            do eta=eta_min,eta_max
               k=eta-(eta_min-1)
               lmatrix(k,1)=1._dblprec ; lmatrix(k,2)=0.1_dblprec*eta ; lmatrix(k,3)=(0.1_dblprec*eta)**2; dvector(k)=DM0_mat_eta(i,j,eta)
            enddo

            call dgels('N',eta_max-(eta_min-1),3,1,lmatrix,eta_max-(eta_min-1),dvector,eta_max-(eta_min-1),awork,alwork,info)
            if(info.ne.0) then
               print '(2x,a,i4)', 'Problem in dgels:',info
            end if
            DM0_mat_lsq(i,j)=dvector(1)
         enddo
      enddo

      i_all=-product(shape(eigenvals))*kind(eigenvals)
      deallocate(eigenvals,stat=i_stat)
      call memocc(i_stat,i_all,'eigenvals','DMI_stiffness')
      i_all=-product(shape(eigenvecs))*kind(eigenvecs)
      deallocate(eigenvecs,stat=i_stat)
      call memocc(i_stat,i_all,'eigenvecs','DMI_stiffness')
      i_all=-product(shape(dm_mat))*kind(dm_mat)
      deallocate(dm_mat,stat=i_stat)
      call memocc(i_stat,i_all,'dm_mat','DMI_stiffness')
      i_all=-product(shape(DM0_mat_eta))*kind(DM0_mat_eta)
      deallocate(DM0_mat_eta,stat=i_stat)
      call memocc(i_stat,i_all,'DM0_mat_eta','DMI_stiffness')
      i_all=-product(shape(awork))*kind(awork)
      deallocate(awork,stat=i_stat)
      call memocc(i_stat,i_all,'awork','DMI_stiffness')

   end subroutine DMI_stiffness

   !---------------------------------------------------------------------------
   !> @brief
   !> Fitting via a rational polynomial approach (from Numerical Recipes)
   !
   !> @author
   !> Manuel Pereiro 
   !---------------------------------------------------------------------------
   subroutine ratint(xa,ya,n,x,y,dy)
      ! Largest expected value of n, and a small number.
      ! Given arrays xa and ya, each of length n, and given a value of x,
      ! this routine returns a value of y and an accuracy estimate dy. The
      ! value returned is that of the diagonal rational function, evaluated at x,
      ! which passes through the n points (xai, yai), i = 1...n.
      implicit none

      ! .. Input variables
      integer :: n
      real(dblprec), intent(in) :: xa(:),ya(:),x

      ! .. Output variables
      real(dblprec), intent(out) :: y,dy

      ! .. Local variables
      integer :: i,m,ns
      real :: dd,h,hh,t,w
      real(dblprec) :: TINY
      parameter (TINY=1.e-25)
      real(dblprec), dimension(n) :: c,d

      ns=1
      hh=abs(x-xa(1))

      do i=1,n
         h=abs(x-xa(i))
         if (h.eq.0.) then
            y=ya(i)
            dy=0.0
            return
         else if (h.lt.hh) then
            ns=i
            hh=h
         endif
         c(i)=ya(i)
         d(i)=ya(i)+TINY !The TINY part is needed to prevent a rare zero over zero condition.
      end do

      y=ya(ns)
      ns=ns-1

      do m=1,n-1
         do i=1,n-m
            w=c(i+1)-d(i)
            h=xa(i+m)-x !h will never be zero, since this was tested in the initialization loop
            w=c(i+1)-d(i)
            h=xa(i+m)-x
            t=(xa(i)-x)*d(i)/h
            dd=t-c(i+1)
            !if(dd.eq.0.) pause 'failure in ratint' !This error condition indicates that the interpolating
            if(dd.eq.0.) stop 'failure in ratint' !This error condition indicates that the interpolating
            dd=w/dd                                !function has a pole at the requested value of x.
            d(i)=c(i+1)*dd
            c(i)=t*dd
         end do

         if (2*ns.lt.n-m) then
            dy=c(ns+1)
         else
            dy=d(ns)
            ns=ns-1
         endif

         y=y+dy
      end do

   end subroutine ratint

   !---------------------------------------------------------------------------
   !> @brief
   !> Calculate the exchange stiffness and the spin wave stifness for random alloy
   !
   !> @author
   !> Lars Bergqvist
   !>
   !> @date 08/02/2018 - Lars Bergqvist
   !---------------------------------------------------------------------------
   subroutine ferro_random_stiffness(NT,NA,N1,N2,N3,hdim,mconf,Natom,nHam,Nchmax,eta_max,eta_min,&
         conf_num,do_ralloy,Natom_full,max_no_neigh,Nch,anumb,atype,nlistsize,&
         atype_ch,asite_ch,achem_ch,nlist,alat,C1,C2,C3,coord,chconc,ammom_inp,ncoup,aham, &
         Axc_fit_alloy,Dxc_fit_alloy,Tc_alloy,J0_matrix_alloy)
         

      use Constants
      use Math_functions, only : f_wrap_coord_diff

      implicit none

      ! .. Input variables
      integer, intent(in) :: NT  !< Number of types of atoms
      integer, intent(in) :: NA  !< Number of atoms in one cell
      integer, intent(in) :: N1  !< Number of cell repetitions in x direction
      integer, intent(in) :: N2  !< Number of cell repetitions in y direction
      integer, intent(in) :: N3  !< Number of cell repetitions in z direction
      integer, intent(in) :: hdim   !< Number of elements in Hamiltonian element (scalar or vector)
      integer, intent(in) :: mconf  !< LSF ground state configuration
      integer, intent(in) :: Natom  !< Number of atoms in system
      integer, intent(in) :: nHam   !< Number of atoms in Hamiltonian
      integer, intent(in) :: Nchmax !< Max number of chemical components on each site in cell
      integer, intent(in) :: eta_max  !< Number of convergence parameters for the stiffness
      integer, intent(in) :: eta_min  !< Minimum  convergence parameters for the stiffness
      integer, intent(in) :: conf_num !< Number of LSF configurations
      integer, intent(in) :: do_ralloy    !< Random alloy simulation (0/1)
      integer, intent(in) :: Natom_full   !< Number of atoms for full system (=Natom if not dilute)
      integer, intent(in) :: max_no_neigh !< Calculated maximum of neighbours for exchange
      integer, dimension(NA), intent(in) :: Nch !< Number of chemical components on each site in cell
      integer, dimension(Natom), intent(in) :: anumb !< Atom number in cell
      integer, dimension(Natom), intent(in) :: atype !< Type of atom
      integer, dimension(nHam), intent(in) :: nlistsize !< Size of neighbour list for Heisenberg exchange couplings
      integer, dimension(Natom_full), intent(in) :: atype_ch !< Actual site of atom for dilute system
      integer, dimension(Natom_full), intent(in) :: asite_ch !< Actual site of atom for dilute system
      integer, dimension(Natom_full), intent(in) :: achem_ch !< Chemical type of atoms (reduced list)
      integer, dimension(max_no_neigh,Natom), intent(in) :: nlist  !< Neighbour list for Heisenberg exchange couplings

      real(dblprec), intent(in) :: alat !< Lattice parameter
      real(dblprec), dimension(3), intent(in) :: C1 !< First lattice vector
      real(dblprec), dimension(3), intent(in) :: C2 !< Second lattice vector
      real(dblprec), dimension(3), intent(in) :: C3 !< Third lattice vector
      real(dblprec), dimension(3,Natom), intent(in) :: coord !< Coordinates of atoms
      real(dblprec), dimension(NA,Nchmax), intent(in) :: chconc !< Chemical concentration on sites
      real(dblprec), dimension(NA,Nchmax,conf_num), intent(in) :: ammom_inp !< Magnetic moment directions from input (for alloys)
      real(dblprec), dimension(hdim,max_no_neigh,nHam,conf_num), intent(in) :: ncoup !< Heisenberg exchange couplings
      integer, dimension(Natom), intent(in) :: aham !< Hamiltonian look-up table

      ! .. Output variables
      real(dblprec),dimension(natom),intent(out) :: Axc_fit_alloy !< Exchange stiffness alloy
      real(dblprec),dimension(natom),intent(out) :: Dxc_fit_alloy !< Spin wave stiffness alloy
      real(dblprec),dimension(natom),intent(out) :: Tc_alloy  !< Tc-MFA alloy
      real(dblprec),dimension(na,na,natom),intent(out) :: J0_matrix_alloy !< Exchange matrix alloy

      ! .. Local variables
      integer :: i_stat, i_all, ia
      integer :: lwork, info, alwork
      integer :: iatom, jatom, iham
      integer :: katom, I1, I2, I3, countstart
      integer :: i, k, eta, ich

      real(dblprec) :: M_sat    !< Saturation magnetization muB/m^3
      real(dblprec) :: cell_vol !< Unit cell volume m^3
      real(dblprec) :: total_mom
      real(dblprec) :: fcinv, jij, jijsign
      real(dblprec) :: rij2, rij, stiff_par

      real(dblprec), dimension(3) :: rcoord
      real(dblprec), dimension(eta_max) :: eig_val, eig_val_temp
      real(dblprec), dimension(eta_max-(eta_min-1),3) :: lmatrix
      real(dblprec), dimension(eta_max-(eta_min-1)) :: dvector

      real(dblprec), dimension(:), allocatable :: wres
      real(dblprec), dimension(:), allocatable :: awork
      real(dblprec), dimension(:,:,:), allocatable :: etemp
      real(dblprec), dimension(:,:,:), allocatable :: stiff_matrix !< Matrix being used to calculate the ferromagnetic exchange stiffness

      real(dblprec), dimension(:), allocatable :: work
      real(dblprec), dimension(:), allocatable :: rwork
      real(dblprec), dimension(:), allocatable :: iwres
      real(dblprec), dimension(:), allocatable :: rwres
      real(dblprec), dimension(:), allocatable :: cwres
      real(dblprec), dimension(:,:),allocatable :: ctemp
      real(dblprec), dimension(:,:),allocatable :: A_inplace


      lwork=16*na
      alwork=((eta_max-(eta_min-1))*16)

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Allocate work arrays
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      allocate(cwres(NA),stat=i_stat)
      call memocc(i_stat,product(shape(cwres))*kind(cwres),'cwres','ferro_stiffness')
      allocate(rwres(NA),stat=i_stat)
      call memocc(i_stat,product(shape(rwres))*kind(rwres),'rwres','ferro_stiffness')
      allocate(iwres(NA),stat=i_stat)
      call memocc(i_stat,product(shape(iwres))*kind(iwres),'iwres','ferro_stiffness')
      allocate(work(lwork),stat=i_stat)
      call memocc(i_stat,product(shape(work))*kind(work),'work','ferro_stiffness')
      allocate(ctemp(NA,NA),stat=i_stat)
      call memocc(i_stat,product(shape(ctemp))*kind(ctemp),'ctemp','ferro_stiffness')
      allocate(rwork(lwork),stat=i_stat)
      call memocc(i_stat,product(shape(rwork))*kind(rwork),'rwork','ferro_stiffness')
      allocate(wres(eta_max),stat=i_stat)
      call memocc(i_stat,product(shape(wres))*kind(wres),'wres','ferro_stiffness')
      allocate(awork(alwork),stat=i_stat)
      call memocc(i_stat,product(shape(awork))*kind(awork),'awork','ferro_stiffness')
      allocate(A_inplace(NA,NA),stat=i_stat)
      call memocc(i_stat,product(shape(A_inplace))*kind(A_inplace),'A_inplace','ferro_stiffness')
      allocate(etemp(NA,NA,eta_max),stat=i_stat)
      call memocc(i_stat,product(shape(etemp))*kind(etemp),'etemp','ferro_stiffness')
      allocate(stiff_matrix(NA,NA,eta_max),stat=i_stat)
      call memocc(i_stat,product(shape(stiff_matrix))*kind(stiff_matrix),'stiff_matrix','ferro_stiffness')

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Initialization of variables
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      work=0.0_dblprec
      wres=0.0_dblprec
      etemp=0.0_dblprec
      ctemp=0.0_dblprec
      rwork=0.0_dblprec
      cwres=0.0_dblprec
      iwres=0.0_dblprec
      rwres=0.0_dblprec
      eig_val=0.0_dblprec
      lmatrix=0.0_dblprec
      dvector=0.0_dblprec
      total_mom=0.0_dblprec
      A_inplace=0.0_dblprec
      eig_val_temp=0.0_dblprec
      stiff_matrix=0.0_dblprec
      D_err_fit=0.0_dblprec

      Axc_fit_alloy=0.0_dblprec
      Dxc_fit_alloy=0.0_dblprec
      Tc_alloy=0.0_dblprec
      J0_matrix_alloy=0.0_dblprec

      fcinv=mub/mry
      I1 = N1/2
      I2 = N2/2
      I3 = N3/2
      countstart = 0+I1*NA+I2*N1*NA+I3*N2*N1*NA
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Calculate the volume of the cell in meters
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      cell_vol=(C1(1)*C2(2)*C3(3)-C1(1)*C2(3)*C3(2)+&
         C1(2)*C2(3)*C3(1)-C1(2)*C2(1)*C3(3)+&
         C1(3)*C2(1)*C3(2)-C1(3)*C2(2)*C3(1))*alat**3
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Calculating the total moment of the cell
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      do i=1,NA
         do ich=1,Nch(i)
            total_mom=total_mom+ammom_inp(i,ich,mconf)*chconc(i,ich)
         enddo
      enddo

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Calculating the saturation magnetization
      M_sat=total_mom/cell_vol
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!$omp parallel do default(shared),private(ia,k,eta,i,iatom,iham,jatom,katom,rcoord,rij2,rij,jij,jijsign,stiff_par,A_inplace,rwres,iwres,ctemp,work,info,lmatrix,dvector,awork,eig_val,etemp),reduction(+:stiff_matrix,J0_matrix_alloy)
      do ia=1,natom
         stiff_matrix=0.0_dblprec
         do eta=1, eta_max

            ! Loop over atoms the atoms in the unit cell
            do i=1,NA
               iatom=i+(ia-1)
               iham=aham(iatom)
               ! Loop over the neighbors for each atom in the unit cell
               do jatom=1, nlistsize(iham)

                  ! Neighbouring atom
                  katom=nlist(jatom,iatom)

                  ! Distace vector between the atoms
                  call f_wrap_coord_diff(Natom,coord,katom,iatom,rcoord)

                  ! Distance vector between neighbouring atoms
                  rij2=rcoord(1)**2+rcoord(2)**2+rcoord(3)**2
                  ! Distance between neighbouring atoms
                  rij=sqrt(rij2)

                  ! Calculating the "real" Jijs in mRyd
                  jij=ncoup(1,jatom,iham,mconf)*fcinv*ammom_inp(asite_ch(katom),achem_ch(katom),mconf)*&
                     ammom_inp(asite_ch(iatom),achem_ch(iatom),mconf)*0.5_dblprec
                  ! Which is the sign between the magnetic moments (To take care of AFM interactions)
                  jijsign=sign(ammom_inp(asite_ch(katom),achem_ch(katom),mconf),ammom_inp(asite_ch(iatom),achem_ch(iatom),mconf))/&
                     abs(ammom_inp(asite_ch(katom),achem_ch(katom),mconf))

                  ! Parameter for the stiffness
                  stiff_par=(2.0_dblprec/3.0_dblprec)*jij*jijsign*rij2*(alat**2)/(sqrt(abs(ammom_inp(asite_ch(katom),achem_ch(katom),mconf))*&
                     abs(ammom_inp(asite_ch(iatom),achem_ch(iatom),mconf))))
                  ! The actual stiffness matrix
                  stiff_matrix(i,asite_ch(katom),eta)=stiff_matrix(i,asite_ch(katom),eta)+&
                     stiff_par*exp(-0.10_dblprec*eta*rij)

                  ! Filling the J0 matrix
                  J0_matrix_alloy(i,asite_ch(katom),ia)=J0_matrix_alloy(i,asite_ch(katom),ia)+(2.0_dblprec/3.0_dblprec)*jij*jijsign
               enddo

            enddo

            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            ! Calculating the eigenvalues for the scalar stiffness
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            ! The matrix is transformed to the more familiar units of meVÅ^2
            A_inplace(1:NA,1:NA)=stiff_matrix(1:NA,1:NA,eta)*1d20*ry_ev

            ! The eigenvalues for the spin wave stiffness are calculated using LAPACK
            call dgeev('N','N',NA, A_inplace, NA, rwres, iwres, ctemp, NA, etemp(1,1,1), NA, WORK, LWORK, INFO)
            if(info.ne.0) then
               print '(2x,a,i4)', 'Problem in zgeev 4:',info
            end if

            ! Finding the maximum eigenvalue of the exchange stiffness matrix
            eig_val(eta)=maxval((rwres))
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            ! End of calculation of eignevalues for the scalar stiffness
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

         enddo
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         ! Calculating the scalar exchange stiffness with LSQ
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         ! Creating a vector for the eta parameter
         do eta=eta_min,eta_max
            k=eta-(eta_min-1)
            lmatrix(k,1)=1._dblprec ; lmatrix(k,2)=0.1_dblprec*eta ; lmatrix(k,3)=(0.1_dblprec*eta)**2; dvector(k)=eig_val(eta)
         enddo

         call dgels('N',eta_max-(eta_min-1),3,1,lmatrix,eta_max-(eta_min-1),dvector,eta_max-(eta_min-1),awork,alwork,info)
         if(info.ne.0) then
            print '(2x,a,i4)', 'Problem in dgels:',info
         end if
         ! Calculate the Exchange stiffness from micromagnetics
         Axc_fit_alloy(ia)=1e12*dvector(1)*1e-20*M_sat/(1000*Joule_ev*4)

         Dxc_fit_alloy(ia)=dvector(1)

         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         ! End of calculation of the scalar exchange stiffness with LSQ
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         ! Calculation of the MF Tc
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

         J0_matrix_alloy(:,:,ia)=J0_matrix_alloy(:,:,ia)/eta_max
         A_inplace(1:NA,1:NA)=(J0_matrix_alloy(:,:,ia)*ry_ev*1e-3)/k_bolt_ev
         ! The eigenvalues for the spin wave stiffness are calculated using LAPACK
         call dgeev('N','N',NA, A_inplace, NA, rwres, iwres, ctemp, NA, etemp(1,1,1), NA, WORK, LWORK, INFO)
         if(info.ne.0) then
            print '(2x,a,i4)', 'Problem in zgeev 5:',info
         end if
         Tc_alloy(ia)=maxval((rwres))
!         write(*,1009) 'Tc-MFA from stiffness :',Tc_alloy(ia),'K'

         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         ! End of calculation of the MF Tc
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      enddo
!$omp end parallel do

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! End of stiffness calculations
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      ! The issue for the chemical dissorder is that there might not be different chemical types in the same cell
      i_all=-product(shape(work))*kind(work)
      deallocate(work,stat=i_stat)
      call memocc(i_stat,i_all,'work','ferro_stiffness')
      i_all=-product(shape(wres))*kind(wres)
      deallocate(wres,stat=i_stat)
      call memocc(i_stat,i_all,'wres','ferro_stiffness')
      i_all=-product(shape(cwres))*kind(cwres)
      deallocate(cwres,stat=i_stat)
      call memocc(i_stat,i_all,'cwres','ferro_stiffness')
      i_all=-product(shape(ctemp))*kind(ctemp)
      deallocate(ctemp,stat=i_stat)
      call memocc(i_stat,i_all,'ctemp','ferro_stiffness')
      i_all=-product(shape(rwork))*kind(rwork)
      deallocate(rwork,stat=i_stat)
      call memocc(i_stat,i_all,'rwork','ferro_stiffness')
      i_all=-product(shape(awork))*kind(awork)
      deallocate(awork,stat=i_stat)
      call memocc(i_stat,i_all,'awork','ferro_stiffness')
      i_all=-product(shape(etemp))*kind(etemp)
      deallocate(etemp,stat=i_stat)
      call memocc(i_stat,i_all,'etemp','ferro_stiffness')
      i_all=-product(shape(A_inplace))*kind(A_inplace)
      deallocate(A_inplace,stat=i_stat)
      call memocc(i_stat,i_all,'A_inplace','ferro_stiffness')
      i_all=-product(shape(stiff_matrix))*kind(stiff_matrix)
      deallocate(stiff_matrix,stat=i_stat)
      call memocc(i_stat,i_all,'stiff_matrix','ferro_stiffness')

      1009 format (2x,a,f10.1,2x,a)
   end subroutine ferro_random_stiffness

end module stiffness
