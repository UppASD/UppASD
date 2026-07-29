program test_stiffness_material

   use Parameters, only : dblprec
   use stiffness

   implicit none

   integer :: failures

   failures = 0
   call test_uppasd_adapter()
   call test_nearest_neighbor_fm_and_small_q()
   call test_anisotropic_exchange()
   call test_skew_cell()
   call test_regularization_diagnostics()
   call test_two_basis_raw_data_and_runtime_gate()

   if (failures /= 0) then
      write(*,'(a,i0)') 'coarse stiffness material tests failed: ', failures
      stop 1
   end if
   write(*,'(a)') 'coarse stiffness material tests passed'

contains

   subroutine test_uppasd_adapter()
      type(coarse_material_type) :: material
      integer :: status
      character(len=512) :: message
      integer :: anumb(3), aham(3), nlistsize(1), dmlistsize(1)
      integer :: nlist(2,3), dmlist(2,3), basis_to_channel(1)
      real(dblprec), parameter :: a = 2.0d-10
      real(dblprec), parameter :: exchange_k = 4.0d-21
      real(dblprec), parameter :: dmi_l = 1.0d-24
      real(dblprec) :: C1(3), C2(3), C3(3), coord(3,3)
      real(dblprec) :: ammom_inp(1,1,1), ncoup(1,2,1,1), dm_vect(3,2,1)
      real(dblprec) :: channel_gamma(1), channel_damping(1)

      anumb = 1
      aham = 1
      nlistsize = 2
      dmlistsize = 2
      nlist = 1
      dmlist = 1
      nlist(:,2) = (/1,3/)
      dmlist(:,2) = (/1,3/)
      basis_to_channel = 1
      C1 = (/1.0d0,0.0d0,0.0d0/)
      C2 = (/0.0d0,1.0d0,0.0d0/)
      C3 = (/0.0d0,0.0d0,1.0d0/)
      coord = 0.0_dblprec
      coord(1,:) = (/0.0d0,1.0d0,2.0d0/)
      ammom_inp = 1.0_dblprec
      ncoup = exchange_k/COARSE_MUB_SI
      dm_vect = 0.0_dblprec
      dm_vect(3,1,1) = -dmi_l/COARSE_MUB_SI
      dm_vect(3,2,1) = dmi_l/COARSE_MUB_SI
      channel_gamma = 3.5d11
      channel_damping = 0.05d0

      call extract_coarse_material_from_uppasd(1,3,1,1,3,1,1,1,2,2,0,0, &
         anumb,aham,nlistsize,nlist,dmlistsize,dmlist,a,C1,C2,C3,'P','O','O', &
         coord,ammom_inp,ncoup,dm_vect,basis_to_channel,material,status,message, &
         channel_gamma,channel_damping)
      call check(status == COARSE_MATERIAL_OK .and. material%ready, &
         'ordered UppASD Hamiltonian adapter supplies a typed material')
      if (status /= COARSE_MATERIAL_OK) return
      call check_close(material%exchange_stiffness(1,1,1,1),exchange_k/(2*a),1.0d-13, &
         'UppASD field coefficient converts to SI exchange stiffness')
      call check_close(material%spiralization(3,1,1,1), &
         dmi_l*a/(a**3),1.0d-13, &
         'UppASD DMI field coefficient converts with the approved sign')
   end subroutine test_uppasd_adapter

   subroutine test_nearest_neighbor_fm_and_small_q()
      type(coarse_lattice_input_type) :: input
      type(coarse_material_type) :: material
      integer :: status, pair
      character(len=512) :: message
      real(dblprec), parameter :: a = 2.0d-10
      real(dblprec), parameter :: exchange_k = 4.0d-21
      real(dblprec), parameter :: dmi_l = 1.0d-24
      real(dblprec) :: expected_a, expected_d
      real(dblprec) :: qvec(3,4), normal(3,4), phase(1,4)

      call initialize_input(input,1,1,a**3,6,2,1)
      do pair = 1, 6
         input%exchange_source_basis(pair) = 1
         input%exchange_target_basis(pair) = 1
         input%exchange_pair_energy_j(pair) = exchange_k
      end do
      input%exchange_displacement_m = 0.0_dblprec
      input%exchange_displacement_m(1,1:2) = (/a,-a/)
      input%exchange_displacement_m(2,3:4) = (/a,-a/)
      input%exchange_displacement_m(3,5:6) = (/a,-a/)

      input%dmi_source_basis = 1
      input%dmi_target_basis = 1
      input%dmi_displacement_m = 0.0_dblprec
      input%dmi_displacement_m(1,:) = (/a,-a/)
      input%dmi_pair_energy_j = 0.0_dblprec
      input%dmi_pair_energy_j(3,:) = (/dmi_l,-dmi_l/)

      call extract_coarse_material(input,material,status,message)
      call check(status == COARSE_MATERIAL_OK .and. material%ready, &
         'nearest-neighbour FM material extracts')
      if (.not. material%ready) return

      expected_a = exchange_k/(2.0_dblprec*a)
      expected_d = dmi_l*a/input%cell_volume_m3
      call check_close(material%exchange_stiffness(1,1,1,1),expected_a,1.0d-13, &
         'simple-cubic nearest-neighbour A_xx is analytic')
      call check_close(material%exchange_stiffness(2,2,1,1),expected_a,1.0d-13, &
         'simple-cubic nearest-neighbour A_yy is analytic')
      call check_close(material%exchange_stiffness(3,3,1,1),expected_a,1.0d-13, &
         'simple-cubic nearest-neighbour A_zz is analytic')
      call check_close(material%spiralization(3,1,1,1),expected_d,1.0d-13, &
         'DMI handedness gives positive D_zx for D(+x)=+z')
      call check(material%metadata%dmi_energy_sign == COARSE_DMI_ENERGY_PLUS .and. &
         material%metadata%field_derivative_sign == COARSE_FIELD_DERIVATIVE_MINUS, &
         'CG-01 DMI and field signs are machine-readable')
      call check(index(material%metadata%exchange_stiffness_unit,'J m^-1') > 0 .and. &
         index(material%metadata%spiralization_order,'cross_spin_k') > 0, &
         'tensor units and ordering are self-describing')
      call check(size(material%raw%basis_exchange_stiffness,3) == 1 .and. &
         size(material%raw%channel_exchange_stiffness,4) == 1, &
         'raw basis and channel tensors remain exposed')
      call check(-expected_d/(2.0_dblprec*expected_a) < 0.0_dblprec, &
         'positive D_zx selects the CG-01 negative-q handedness')

      qvec = 0.0_dblprec
      qvec(1,:) = (/1.0d5,-1.0d5,2.0d5,-2.0d5/)
      normal = 0.0_dblprec
      normal(3,:) = 1.0_dblprec
      phase = 0.0_dblprec
      call validate_coarse_material_small_q(input,material,qvec,normal,phase, &
         1.0d-7,1.0d-7,status,message)
      call check(status == COARSE_MATERIAL_OK .and. &
         material%diagnostics%small_q_energy_validated, &
         'direct small-q atomistic FM energies match A q^2 + D q')
      call check(material%diagnostics%small_q_atomistic_energy_j_per_m3(1) > 0.0_dblprec .and. &
         material%diagnostics%small_q_atomistic_energy_j_per_m3(2) < 0.0_dblprec, &
         'direct DMI spiral energies retain the approved handedness')
      call coarse_material_runtime_status(material,COARSE_RUNTIME_SINGLE_FM,status,message)
      call check(status == COARSE_MATERIAL_OK, &
         'validated one-channel FM descriptor passes the extraction runtime gate')
   end subroutine test_nearest_neighbor_fm_and_small_q

   subroutine test_anisotropic_exchange()
      type(coarse_lattice_input_type) :: input
      type(coarse_material_type) :: material
      integer :: status, pair
      character(len=512) :: message
      real(dblprec), parameter :: a = 2.5d-10
      real(dblprec), parameter :: kx = 6.0d-21
      real(dblprec), parameter :: ky = 3.0d-21
      real(dblprec), parameter :: kz = 1.5d-21

      call initialize_input(input,1,1,a**3,6,0,1)
      input%exchange_displacement_m = 0.0_dblprec
      input%exchange_displacement_m(1,1:2) = (/a,-a/)
      input%exchange_displacement_m(2,3:4) = (/a,-a/)
      input%exchange_displacement_m(3,5:6) = (/a,-a/)
      do pair = 1, 6
         input%exchange_source_basis(pair) = 1
         input%exchange_target_basis(pair) = 1
      end do
      input%exchange_pair_energy_j = (/kx,kx,ky,ky,kz,kz/)

      call extract_coarse_material(input,material,status,message)
      call check(status == COARSE_MATERIAL_OK, 'anisotropic exchange material extracts')
      if (status /= COARSE_MATERIAL_OK) return
      call check_close(material%exchange_stiffness(1,1,1,1),kx/(2*a),1.0d-13, &
         'anisotropic A_xx follows the x-bond exchange')
      call check_close(material%exchange_stiffness(2,2,1,1),ky/(2*a),1.0d-13, &
         'anisotropic A_yy follows the y-bond exchange')
      call check_close(material%exchange_stiffness(3,3,1,1),kz/(2*a),1.0d-13, &
         'anisotropic A_zz follows the z-bond exchange')
      call check(maxval(abs(material%exchange_stiffness(:,:,1,1) - &
         diagonal3((/kx/(2*a),ky/(2*a),kz/(2*a)/)))) < 1.0d-24, &
         'anisotropic nearest-neighbour tensor has no spurious off-diagonal terms')
   end subroutine test_anisotropic_exchange

   subroutine test_skew_cell()
      type(coarse_lattice_input_type) :: input
      type(coarse_material_type) :: material
      integer :: status, pair, direction
      character(len=512) :: message
      real(dblprec), parameter :: exchange_k = 3.0d-21
      real(dblprec) :: cell(3,3), expected(3,3), volume

      cell(:,1) = (/2.0d-10,0.0d0,0.0d0/)
      cell(:,2) = (/0.6d-10,1.5d-10,0.0d0/)
      cell(:,3) = (/0.2d-10,0.3d-10,1.2d-10/)
      volume = abs(det3(cell))
      call initialize_input(input,1,1,volume,6,0,1)
      do pair = 1, 6
         input%exchange_source_basis(pair) = 1
         input%exchange_target_basis(pair) = 1
         input%exchange_pair_energy_j(pair) = exchange_k
      end do
      do direction = 1, 3
         input%exchange_displacement_m(:,2*direction-1) = cell(:,direction)
         input%exchange_displacement_m(:,2*direction) = -cell(:,direction)
      end do

      call extract_coarse_material(input,material,status,message)
      call check(status == COARSE_MATERIAL_OK, 'skew-cell material extracts')
      if (status /= COARSE_MATERIAL_OK) return
      expected = 0.0_dblprec
      do direction = 1, 3
         expected = expected + exchange_k/(2.0_dblprec*volume) * &
            outer_product(cell(:,direction),cell(:,direction))
      end do
      call check(maxval(abs(material%exchange_stiffness(:,:,1,1)-expected)) < 1.0d-24, &
         'skew-cell sum uses physical Cartesian displacements and volume')
      call check(abs(material%exchange_stiffness(1,2,1,1)) > 1.0d-14, &
         'skew-cell tensor retains its analytic Cartesian cross term')
   end subroutine test_skew_cell

   subroutine test_regularization_diagnostics()
      type(coarse_lattice_input_type) :: input
      type(coarse_material_type) :: material
      integer :: status
      character(len=512) :: message
      real(dblprec), parameter :: a = 2.0d-10
      real(dblprec), parameter :: exchange_k = 4.0d-21
      real(dblprec) :: expected

      call initialize_input(input,1,1,a**3,2,0,4)
      input%eta_inverse_m = (/1.0d7,2.0d7,3.0d7,4.0d7/)
      input%fit_first = 1
      input%fit_last = 4
      input%exchange_source_basis = 1
      input%exchange_target_basis = 1
      input%exchange_pair_energy_j = exchange_k
      input%exchange_displacement_m = 0.0_dblprec
      input%exchange_displacement_m(1,:) = (/a,-a/)

      call extract_coarse_material(input,material,status,message)
      call check(status == COARSE_MATERIAL_OK .and. &
         material%diagnostics%fit_performed, &
         'regularized sums use a separately reported extrapolation fit')
      if (status /= COARSE_MATERIAL_OK) return
      expected = exchange_k/(2.0_dblprec*a)
      call check_close(material%exchange_stiffness(1,1,1,1),expected,5.0d-8, &
         'quadratic eta extrapolation recovers the nearest-neighbour intercept')
      call check(material%diagnostics%fit_sample_count == 4 .and. &
         abs(material%diagnostics%eta_min_inverse_m-1.0d7) < 1.0d-12 .and. &
         abs(material%diagnostics%eta_max_inverse_m-4.0d7) < 1.0d-12, &
         'fit range and convergence parameters are retained')
      call check(allocated(material%diagnostics%exchange_fit_coeff) .and. &
         allocated(material%diagnostics%exchange_fit_rms) .and. &
         material%diagnostics%regularization_exchange_delta_j_per_m > 0.0_dblprec, &
         'fit coefficients, residuals, and regularization delta are exposed')
      call check(maxval(abs(material%raw%eta_inverse_m-input%eta_inverse_m)) < 1.0d-12, &
         'all raw convergence samples remain available before fitting')
   end subroutine test_regularization_diagnostics

   subroutine test_two_basis_raw_data_and_runtime_gate()
      type(coarse_lattice_input_type) :: input
      type(coarse_material_type) :: material
      integer :: status
      character(len=512) :: message
      real(dblprec), parameter :: a = 2.0d-10
      real(dblprec), parameter :: exchange_k = 2.0d-21
      real(dblprec) :: expected_cross_a, expected_local
      real(dblprec) :: qvec(3,4), normal(3,4), phase(2,4)

      call initialize_input(input,2,2,a**3,4,0,1)
      input%basis_to_channel = (/1,2/)
      input%basis_moment_mub = (/1.0d0,1.5d0/)
      input%channel_gamma = (/2.0d0,2.0d0/)
      input%channel_damping = (/0.05d0,0.05d0/)
      input%exchange_source_basis = (/1,1,2,2/)
      input%exchange_target_basis = (/2,2,1,1/)
      input%exchange_pair_energy_j = exchange_k
      input%exchange_displacement_m = 0.0_dblprec
      input%exchange_displacement_m(1,:) = (/0.5d0*a,-0.5d0*a,0.5d0*a,-0.5d0*a/)

      call extract_coarse_material(input,material,status,message)
      call check(status == COARSE_MATERIAL_OK .and. material%n_channels == 2, &
         'two-basis/two-channel material extracts without eigenvalue collapse')
      if (status /= COARSE_MATERIAL_OK) return
      expected_cross_a = exchange_k*a*a/(8.0_dblprec*input%cell_volume_m3)
      expected_local = -2.0_dblprec*exchange_k/input%cell_volume_m3
      call check_close(material%raw%basis_exchange_stiffness(1,1,1,2,1), &
         expected_cross_a,1.0d-13,'basis-resolved A_12 remains exposed')
      call check_close(material%raw%basis_exchange_stiffness(1,1,2,1,1), &
         expected_cross_a,1.0d-13,'basis-resolved A_21 remains exposed')
      call check_close(material%exchange_stiffness(1,1,1,2),expected_cross_a,1.0d-13, &
         'channel A_12 is retained rather than replaced by a maximum eigenvalue')
      call check_close(material%local_exchange(1,2),expected_local,1.0d-13, &
         'two-channel local interbasis exchange is retained in J m^-3')
      call check(all(abs(material%channel_moment_mub-(/1.0d0,1.5d0/)) < 1.0d-14), &
         'channel moments remain separate for a two-basis material')

      qvec = 0.0_dblprec
      qvec(1,:) = (/1.0d5,-1.0d5,1.5d5,-1.5d5/)
      normal = 0.0_dblprec
      normal(3,:) = 1.0_dblprec
      phase = 0.0_dblprec
      phase(2,3:4) = acos(-1.0_dblprec)
      call validate_coarse_material_small_q(input,material,qvec,normal,phase, &
         1.0d-8,1.0d-7,status,message)
      call check(status == COARSE_MATERIAL_OK, &
         'two-basis acoustic and optical phase-sector small-q energies agree')
      call check(.not. material%diagnostics%acoustic_mode_extraction_validated .and. &
         .not. material%diagnostics%optical_mode_extraction_validated, &
         'static spiral checks do not masquerade as dynamical mode extraction')
      call coarse_material_runtime_status(material,COARSE_RUNTIME_FERRI_AFM,status,message)
      call check(status == COARSE_MATERIAL_RUNTIME_GATED .and. &
         index(message,'acoustic and optical') > 0, &
         'ferri/AFM runtime fails explicitly until acoustic and optical modes are validated')
   end subroutine test_two_basis_raw_data_and_runtime_gate

   subroutine initialize_input(input,n_basis,n_channels,volume,n_exchange,n_dmi,n_eta)
      type(coarse_lattice_input_type), intent(out) :: input
      integer, intent(in) :: n_basis, n_channels, n_exchange, n_dmi, n_eta
      real(dblprec), intent(in) :: volume

      input%n_basis = n_basis
      input%n_channels = n_channels
      input%cell_volume_m3 = volume
      input%fit_first = 1
      input%fit_last = n_eta
      allocate(input%basis_to_channel(n_basis), input%basis_moment_mub(n_basis))
      allocate(input%channel_gamma(n_channels), input%channel_damping(n_channels))
      allocate(input%eta_inverse_m(n_eta))
      allocate(input%exchange_source_basis(n_exchange))
      allocate(input%exchange_target_basis(n_exchange))
      allocate(input%exchange_displacement_m(3,n_exchange))
      allocate(input%exchange_pair_energy_j(n_exchange))
      allocate(input%dmi_source_basis(n_dmi))
      allocate(input%dmi_target_basis(n_dmi))
      allocate(input%dmi_displacement_m(3,n_dmi))
      allocate(input%dmi_pair_energy_j(3,n_dmi))
      input%basis_to_channel = 1
      input%basis_moment_mub = 1.0_dblprec
      input%channel_gamma = 2.0_dblprec
      input%channel_damping = 0.05_dblprec
      input%eta_inverse_m = 0.0_dblprec
      input%exchange_source_basis = 0
      input%exchange_target_basis = 0
      input%exchange_displacement_m = 0.0_dblprec
      input%exchange_pair_energy_j = 0.0_dblprec
      input%dmi_source_basis = 0
      input%dmi_target_basis = 0
      input%dmi_displacement_m = 0.0_dblprec
      input%dmi_pair_energy_j = 0.0_dblprec
   end subroutine initialize_input

   function diagonal3(values) result(matrix)
      real(dblprec), intent(in) :: values(3)
      real(dblprec) :: matrix(3,3)

      matrix = 0.0_dblprec
      matrix(1,1) = values(1)
      matrix(2,2) = values(2)
      matrix(3,3) = values(3)
   end function diagonal3

   function outer_product(left,right) result(matrix)
      real(dblprec), intent(in) :: left(3), right(3)
      real(dblprec) :: matrix(3,3)
      integer :: i, j

      do i = 1, 3
         do j = 1, 3
            matrix(i,j) = left(i)*right(j)
         end do
      end do
   end function outer_product

   pure real(dblprec) function det3(matrix) result(determinant)
      real(dblprec), intent(in) :: matrix(3,3)

      determinant = matrix(1,1)*(matrix(2,2)*matrix(3,3)-matrix(2,3)*matrix(3,2)) - &
         matrix(1,2)*(matrix(2,1)*matrix(3,3)-matrix(2,3)*matrix(3,1)) + &
         matrix(1,3)*(matrix(2,1)*matrix(3,2)-matrix(2,2)*matrix(3,1))
   end function det3

   subroutine check_close(actual,expected,relative_tolerance,label)
      real(dblprec), intent(in) :: actual, expected, relative_tolerance
      character(len=*), intent(in) :: label
      real(dblprec) :: scale

      scale = max(abs(expected),1.0d-30)
      call check(abs(actual-expected) <= relative_tolerance*scale,label)
   end subroutine check_close

   subroutine check(condition,label)
      logical, intent(in) :: condition
      character(len=*), intent(in) :: label

      if (.not. condition) then
         failures = failures + 1
         write(*,'(a)') 'FAIL: '//trim(label)
      end if
   end subroutine check

end program test_stiffness_material
