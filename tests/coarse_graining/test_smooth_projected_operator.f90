program test_smooth_projected_operator

   use, intrinsic :: iso_fortran_env, only : int64
   use Parameters, only : dblprec
   use BlockTopology
   use stiffness
   use CoarseTensorOperator
   use SmoothProjectedOperator

   implicit none

   integer :: failures

   failures = 0
   call test_shape_weights_boundaries_and_order()
   call test_channel_separation()
   call test_normalized_adjoint_and_bilinear_derivative()
   call test_block_size_one()
   call test_spiral_mesh_scaling_and_tensor_limit()

   if (failures /= 0) then
      write(*,'(a,i0)') 'smooth projected operator tests failed: ', failures
      stop 1
   end if
   write(*,'(a)') 'smooth projected operator tests passed'

contains

   subroutine test_shape_weights_boundaries_and_order()
      type(block_topology_type) :: topology
      type(smooth_projected_operator_type) :: operator
      integer :: atom, block, status
      integer, parameter :: n = 8
      integer, parameter :: permutation(n) = (/4,1,7,2,8,3,6,5/)
      character(len=512) :: message
      real(dblprec) :: cell(3,3), coordinate(3,n), moment(n)
      real(dblprec) :: coarse(3,1,n), atom_direction(3,n), expected(3)

      call cubic_topology(n,1,topology,cell)
      moment = 1.0_dblprec
      coordinate = 0.0_dblprec
      do atom = 1, n
         coordinate(1,atom) = real(permutation(atom)-1,dblprec)
      end do
      call setup_smooth_projected_operator(operator,topology,coordinate,moment,status,message)
      call check(status == SMOOTH_PROJECTED_OK,'trilinear stencil setup: '//trim(message))
      call check(maxval(abs(sum(operator%shape_weight,dim=1)-1.0_dblprec)) < 2.0d-15, &
         'trilinear shape weights form a partition of unity')
      call check(minval(operator%shape_weight) >= 0.0_dblprec .and. &
         maxval(operator%stencil_block) <= n .and. minval(operator%stencil_block) >= 1, &
         'periodic shape stencils contain valid nonnegative weights')

      do block = 1, n
         expected = (/1.0_dblprec,0.17_dblprec*real(block,dblprec), &
            -0.11_dblprec*real(block,dblprec)/)
         coarse(:,1,block) = expected/sqrt(sum(expected*expected))
      end do
      call prolongate_smooth_directions(operator,coarse,atom_direction,status,message)
      call check(status == SMOOTH_PROJECTED_OK,'permuted atom ordering prolongates')
      do atom = 1, n
         call check(maxval(abs(atom_direction(:,atom)-coarse(:,1,permutation(atom)))) < 2.0d-15, &
            'prolongation preserves caller atom index '//integer_text(atom))
      end do

      coordinate(:,1) = (/-0.25_dblprec,0.0_dblprec,0.0_dblprec/)
      coordinate(:,2) = (/real(n,dblprec)-0.25_dblprec,0.0_dblprec,0.0_dblprec/)
      call setup_smooth_projected_operator(operator,topology,coordinate,moment,status,message)
      call check(status == SMOOTH_PROJECTED_OK,'periodic boundary coordinates setup')
      call check(all(operator%stencil_block(:,1) == operator%stencil_block(:,2)) .and. &
         maxval(abs(operator%shape_weight(:,1)-operator%shape_weight(:,2))) < 2.0d-15, &
         'coordinates differing by one grid period have identical stencils')

      call setup_smooth_projected_operator(operator,topology,coordinate,moment,status,message, &
         boundary_mode=(/SMOOTH_PROJECTED_PERIODIC,99,SMOOTH_PROJECTED_PERIODIC/))
      call check(status == SMOOTH_PROJECTED_INVALID_TOPOLOGY .and. &
         index(message,'periodic') > 0,'nonperiodic boundary modes reject explicitly')
   end subroutine test_shape_weights_boundaries_and_order

   subroutine test_channel_separation()
      type(block_topology_type) :: topology
      type(smooth_projected_operator_type) :: operator
      integer, parameter :: nblock = 4, natom = 2*nblock
      integer :: atom, block, channel, status
      character(len=512) :: message
      real(dblprec) :: cell(3,3), coordinate(3,natom), moment(natom)
      real(dblprec) :: coarse(3,2,nblock), atom_direction(3,natom)

      cell = 0.0_dblprec
      cell(1,1) = 1.0d-10
      cell(2,2) = 1.0d-10
      cell(3,3) = 1.0d-10
      call build_block_topology(topology,REGULAR_REPLICATED_CELL,2,(/nblock,1,1/), &
         natom,(/1,1,1/),cell,(/1,2/),status,message)
      call check(status == BLOCK_TOPOLOGY_OK,'two-channel test topology builds')
      coordinate = 0.0_dblprec
      moment = 1.0_dblprec
      do atom = 1, natom
         block = topology%atom_to_block(atom)
         coordinate(1,atom) = real(block-1,dblprec)
      end do
      do block = 1, nblock
         coarse(:,1,block) = (/1.0_dblprec,0.1_dblprec*real(block,dblprec),0.0_dblprec/)
         coarse(:,2,block) = (/-0.2_dblprec,1.0_dblprec,0.1_dblprec*real(block,dblprec)/)
         coarse(:,1,block) = coarse(:,1,block)/sqrt(sum(coarse(:,1,block)**2))
         coarse(:,2,block) = coarse(:,2,block)/sqrt(sum(coarse(:,2,block)**2))
      end do
      call setup_smooth_projected_operator(operator,topology,coordinate,moment,status,message)
      call prolongate_smooth_directions(operator,coarse,atom_direction,status,message)
      call check(status == SMOOTH_PROJECTED_OK,'two-channel prolongation evaluates')
      do atom = 1, natom
         block = topology%atom_to_block(atom)
         channel = topology%atom_to_dynamic_channel(atom)
         call check(maxval(abs(atom_direction(:,atom)-coarse(:,channel,block))) < 2.0d-15, &
            'prolongation keeps block/channel identities distinct')
      end do
   end subroutine test_channel_separation

   subroutine test_normalized_adjoint_and_bilinear_derivative()
      type(block_topology_type) :: topology
      type(smooth_projected_operator_type) :: operator
      integer, parameter :: ncell = 8, nblock = 4
      integer :: atom, block, component, status
      integer :: bonds(2,ncell), bonds_before(2,ncell)
      integer(int64) :: matrix_before(9*ncell), coordinate_before(3*ncell)
      character(len=512) :: message
      real(dblprec), parameter :: epsilon_fd = 1.0d-6
      real(dblprec) :: cell(3,3), coordinate(3,ncell)
      real(dblprec) :: moment(ncell), matrix(3,3,ncell)
      real(dblprec) :: coarse(3,1,nblock), plus(3,1,nblock), minus(3,1,nblock)
      real(dblprec) :: coarse_field(3,1,nblock), unused_field(3,1,nblock)
      real(dblprec) :: atom_direction(3,ncell), atom_field(3,ncell)
      real(dblprec) :: plus_direction(3,ncell), minus_direction(3,ncell)
      real(dblprec) :: energy, plus_energy, minus_energy, numerical, analytic
      real(dblprec) :: atom_functional_plus, atom_functional_minus

      call cubic_topology(ncell,2,topology,cell)
      moment = [(1.1_dblprec+0.03_dblprec*real(atom,dblprec),atom=1,ncell)]
      coordinate = 0.0_dblprec
      do atom = 1, ncell
         coordinate(1,atom) = (real(atom,dblprec)-0.5_dblprec)/2.0_dblprec-0.5_dblprec
      end do
      call setup_smooth_projected_operator(operator,topology,coordinate,moment,status,message)
      call check(status == SMOOTH_PROJECTED_OK,'normalized-adjoint fixture setup: '//trim(message))

      do block = 1, nblock
         coarse(:,1,block) = (/0.8_dblprec+0.04_dblprec*real(block,dblprec), &
            0.21_dblprec*real(block-2,dblprec), &
            0.35_dblprec-0.03_dblprec*real(block,dblprec)/)
      end do
      do atom = 1, ncell
         bonds(:,atom) = (/atom,modulo(atom,ncell)+1/)
         matrix(:,:,atom) = reshape((/ &
            1.7d-21, 0.2d-21,-0.1d-21, &
           -0.3d-21, 1.2d-21, 0.4d-21, &
            0.1d-21,-0.2d-21, 0.9d-21/),(/3,3/))
         matrix(:,:,atom) = matrix(:,:,atom)*(1.0_dblprec+0.02_dblprec*real(atom,dblprec))
         atom_field(:,atom) = (/0.03_dblprec*real(atom,dblprec), &
            -0.07_dblprec+0.01_dblprec*real(atom,dblprec),0.04_dblprec/)
      end do
      bonds_before = bonds
      matrix_before = transfer(matrix,matrix_before)
      coordinate_before = transfer(coordinate,coordinate_before)

      call evaluate_projected_bilinear(operator,coarse,bonds,matrix,energy,coarse_field, &
         status,message,atom_direction)
      call check(status == SMOOTH_PROJECTED_OK,'general bilinear projected evaluation')
      call check(all(bonds == bonds_before) .and. &
         all(transfer(matrix,matrix_before) == matrix_before) .and. &
         all(transfer(coordinate,coordinate_before) == coordinate_before), &
         'evaluation leaves source interactions, coordinates, and atom ordering unchanged')

      block = 2
      component = 3
      plus = coarse
      minus = coarse
      plus(component,1,block) = plus(component,1,block)+epsilon_fd
      minus(component,1,block) = minus(component,1,block)-epsilon_fd
      call evaluate_projected_bilinear(operator,plus,bonds,matrix,plus_energy,unused_field, &
         status,message)
      call evaluate_projected_bilinear(operator,minus,bonds,matrix,minus_energy,unused_field, &
         status,message)
      numerical = (plus_energy-minus_energy)/(2.0_dblprec*epsilon_fd)
      analytic = -COARSE_MUB_SI*operator%coarse_moment_mub(1,block) * &
         coarse_field(component,1,block)
      call check_close(numerical,analytic,3.0d-8, &
         'projected bilinear field is the finite-difference energy derivative')

      call restrict_projected_atomistic_field(operator,coarse,atom_field,coarse_field, &
         status,message)
      call check(status == SMOOTH_PROJECTED_OK,'standalone adjoint restriction evaluates')
      plus = coarse
      minus = coarse
      plus(component,1,block) = plus(component,1,block)+epsilon_fd
      minus(component,1,block) = minus(component,1,block)-epsilon_fd
      call prolongate_smooth_directions(operator,plus,plus_direction,status,message)
      call prolongate_smooth_directions(operator,minus,minus_direction,status,message)
      atom_functional_plus = -COARSE_MUB_SI * &
         sum(spread(moment,1,3)*atom_field*plus_direction)
      atom_functional_minus = -COARSE_MUB_SI * &
         sum(spread(moment,1,3)*atom_field*minus_direction)
      numerical = (atom_functional_plus-atom_functional_minus)/(2.0_dblprec*epsilon_fd)
      analytic = -COARSE_MUB_SI*operator%coarse_moment_mub(1,block) * &
         coarse_field(component,1,block)
      call check_close(numerical,analytic,3.0d-8, &
         'moment-weighted restriction is the exact normalized-prolongation adjoint')

      coarse = 0.0_dblprec
      call prolongate_smooth_directions(operator,coarse,atom_direction,status,message)
      call check(status == SMOOTH_PROJECTED_ZERO_INTERPOLANT .and. &
         index(message,'cancelling') > 0,'zero normalized interpolants reject diagnostically')
   end subroutine test_normalized_adjoint_and_bilinear_derivative

   subroutine test_block_size_one()
      type(block_topology_type) :: topology
      type(smooth_projected_operator_type) :: operator
      integer, parameter :: n = 6
      integer :: atom, status
      integer :: bonds(2,n)
      integer :: dmi_bond(2,1)
      character(len=512) :: message
      real(dblprec), parameter :: pair_j = 2.3d-21
      real(dblprec), parameter :: dmi_j = 0.4d-21
      real(dblprec) :: cell(3,3), coordinate(3,n), moment(n), matrix(3,3,n)
      real(dblprec) :: dmi_matrix(3,3,1), atom_field(3,n)
      real(dblprec) :: coarse(3,1,n), atom_direction(3,n), field(3,1,n)
      real(dblprec) :: angle, energy, expected

      call cubic_topology(n,1,topology,cell)
      coordinate = 0.0_dblprec
      moment = 1.7_dblprec
      matrix = 0.0_dblprec
      do atom = 1, n
         coordinate(1,atom) = real(atom-1,dblprec)
         angle = 2.0_dblprec*pi*real(atom-1,dblprec)/real(n,dblprec)
         coarse(:,1,atom) = (/cos(angle),sin(angle),0.0_dblprec/)
         bonds(:,atom) = (/atom,modulo(atom,n)+1/)
         matrix(1,1,atom) = pair_j
         matrix(2,2,atom) = pair_j
         matrix(3,3,atom) = pair_j
      end do
      call setup_smooth_projected_operator(operator,topology,coordinate,moment,status,message)
      call evaluate_projected_bilinear(operator,coarse,bonds,matrix,energy,field,status, &
         message,atom_direction)
      expected = -real(n,dblprec)*pair_j*cos(2.0_dblprec*pi/real(n,dblprec))
      call check(status == SMOOTH_PROJECTED_OK,'block-size-one projected evaluation')
      call check(maxval(abs(atom_direction-coarse(:,1,:))) < 2.0d-15, &
         'block-size-one prolongation is the identity')
      call check_close(energy,expected,2.0d-15, &
         'block-size-one projected energy equals the atomistic bond energy')

      coarse = 0.0_dblprec
      coarse(1,1,:) = 1.0_dblprec
      coarse(:,1,2) = (/0.0_dblprec,1.0_dblprec,0.0_dblprec/)
      dmi_bond(:,1) = (/1,2/)
      dmi_matrix = 0.0_dblprec
      dmi_matrix(1,2,1) = -dmi_j
      dmi_matrix(2,1,1) = dmi_j
      call evaluate_projected_bilinear(operator,coarse,dmi_bond,dmi_matrix,energy, &
         field,status,message,atom_field_t=atom_field)
      call check_close(energy,dmi_j,2.0d-15, &
         'antisymmetric bond matrix has the accepted positive DMI energy sign')
      call check_close(atom_field(1,1),-dmi_j/(COARSE_MUB_SI*moment(1)),2.0d-15, &
         'antisymmetric bond matrix has the accepted DMI field handedness')
   end subroutine test_block_size_one

   subroutine test_spiral_mesh_scaling_and_tensor_limit()
      integer, parameter :: ncell = 256
      integer, parameter :: widths(3) = (/2,4,8/)
      real(dblprec), parameter :: a = 2.5d-10
      real(dblprec), parameter :: pair_j = 3.0d-21
      real(dblprec) :: projected_excess(3), tensor_energy(3), rigid_excess(3)
      real(dblprec) :: field_relative_error(3)
      integer :: item

      do item = 1, size(widths)
         call spiral_case(ncell,widths(item),a,pair_j,projected_excess(item), &
            tensor_energy(item),rigid_excess(item),field_relative_error(item))
         call check_close(projected_excess(item),tensor_energy(item),1.2d-2, &
            'smooth and tensor long-wave spiral energies agree for block width '// &
            integer_text(widths(item)))
         call check(field_relative_error(item) < 2.0d-2, &
            'smooth and tensor long-wave spiral fields agree for block width '// &
            integer_text(widths(item)))
      end do
      call check(maxval(projected_excess)/minval(projected_excess) < 1.015_dblprec, &
         'smooth spin-spiral exchange has mesh-independent leading scaling')
      call check(rigid_excess(3)/rigid_excess(1) > 3.8_dblprec, &
         'test-only piecewise-constant projection exposes wrong block-size scaling')
      call check(abs(projected_excess(3)-tensor_energy(3)) > &
         abs(projected_excess(1)-tensor_energy(1)), &
         'smooth/tensor discrepancy grows with qh as a discretization error')
   end subroutine test_spiral_mesh_scaling_and_tensor_limit

   subroutine spiral_case(ncell,width,a,pair_j,projected_excess,tensor_energy, &
         rigid_excess,field_relative_error)
      integer, intent(in) :: ncell, width
      real(dblprec), intent(in) :: a, pair_j
      real(dblprec), intent(out) :: projected_excess, tensor_energy
      real(dblprec), intent(out) :: rigid_excess, field_relative_error
      type(block_topology_type) :: topology
      type(smooth_projected_operator_type) :: projected
      type(coarse_tensor_operator_type) :: tensor
      type(coarse_operator_options_type) :: options
      type(coarse_material_type) :: material
      type(coarse_energy_terms_type) :: tensor_terms
      integer :: atom, block, nblock, status
      integer, allocatable :: bonds(:,:)
      character(len=512) :: message
      real(dblprec) :: cell(3,3), exchange(3,3), dmi(3,3), q, angle
      real(dblprec), allocatable :: coordinate(:,:), moment(:), matrix(:,:,:)
      real(dblprec), allocatable :: coarse3(:,:,:), coarse2(:,:), atom_direction(:,:)
      real(dblprec), allocatable :: projected_field(:,:,:), tensor_field(:,:), zero(:,:)
      real(dblprec), allocatable :: projected_tangent(:,:), tensor_tangent(:,:)
      real(dblprec) :: projected_energy, field_scale

      nblock = ncell/width
      call cubic_topology(ncell,width,topology,cell,a)
      allocate(coordinate(3,ncell),moment(ncell),bonds(2,ncell),matrix(3,3,ncell))
      allocate(coarse3(3,1,nblock),coarse2(3,nblock),atom_direction(3,ncell))
      allocate(projected_field(3,1,nblock),tensor_field(3,nblock),zero(3,nblock))
      allocate(projected_tangent(3,nblock),tensor_tangent(3,nblock))
      coordinate = 0.0_dblprec
      moment = 1.7_dblprec
      matrix = 0.0_dblprec
      q = 2.0_dblprec*pi/(real(ncell,dblprec)*a)
      do atom = 1, ncell
         coordinate(1,atom) = (real(atom,dblprec)-0.5_dblprec)/real(width,dblprec)-0.5_dblprec
         bonds(:,atom) = (/atom,modulo(atom,ncell)+1/)
         matrix(1,1,atom) = pair_j
         matrix(2,2,atom) = pair_j
         matrix(3,3,atom) = pair_j
      end do
      do block = 1, nblock
         angle = q*(real(block,dblprec)-0.5_dblprec)*real(width,dblprec)*a
         coarse3(:,1,block) = (/cos(angle),sin(angle),0.0_dblprec/)
      end do
      coarse2 = coarse3(:,1,:)

      call setup_smooth_projected_operator(projected,topology,coordinate,moment,status,message)
      call check(status == SMOOTH_PROJECTED_OK,'spiral projected setup: '//trim(message))
      call evaluate_projected_bilinear(projected,coarse3,bonds,matrix,projected_energy, &
         projected_field,status,message,atom_direction)
      call check(status == SMOOTH_PROJECTED_OK,'spiral projected evaluation: '//trim(message))
      projected_excess = projected_energy+real(ncell,dblprec)*pair_j

      exchange = 0.0_dblprec
      exchange(1,1) = pair_j/(2.0_dblprec*a)
      exchange(2,2) = exchange(1,1)
      exchange(3,3) = exchange(1,1)
      dmi = 0.0_dblprec
      call make_material(cell,exchange,dmi,material)
      call setup_coarse_tensor_operator(tensor,topology,material,options,status,message)
      call check(status == COARSE_TENSOR_OK,'spiral tensor setup: '//trim(message))
      zero = 0.0_dblprec
      call evaluate_coarse_tensor_operator(tensor,coarse2,tensor_field,tensor_terms, &
         status,message,external_field_t=zero)
      call check(status == COARSE_TENSOR_OK,'spiral tensor evaluation: '//trim(message))
      tensor_energy = tensor_terms%exchange_j

      rigid_excess = 0.0_dblprec
      do atom = 1, ncell
         block = (atom-1)/width+1
         rigid_excess = rigid_excess + pair_j * &
            (1.0_dblprec-dot_product(coarse2(:,block), &
             coarse2(:,modulo(atom,ncell)/width+1)))
      end do

      ! Normalized prolongation removes the radial field gauge.  Compare the
      ! physical tangent fields on a weakly phase-modulated long-wave spiral.
      do block = 1, nblock
         angle = q*(real(block,dblprec)-0.5_dblprec)*real(width,dblprec)*a
         angle = angle+0.05_dblprec*sin(2.0_dblprec*angle)
         coarse3(:,1,block) = (/cos(angle),sin(angle),0.0_dblprec/)
      end do
      coarse2 = coarse3(:,1,:)
      call evaluate_projected_bilinear(projected,coarse3,bonds,matrix,projected_energy, &
         projected_field,status,message)
      call evaluate_coarse_tensor_operator(tensor,coarse2,tensor_field,tensor_terms, &
         status,message,external_field_t=zero)
      do block = 1, nblock
         projected_tangent(:,block) = projected_field(:,1,block) - coarse2(:,block) * &
            dot_product(coarse2(:,block),projected_field(:,1,block))
         tensor_tangent(:,block) = tensor_field(:,block) - coarse2(:,block) * &
            dot_product(coarse2(:,block),tensor_field(:,block))
      end do
      field_scale = maxval(abs(tensor_tangent))
      field_relative_error = maxval(abs(projected_tangent-tensor_tangent))/field_scale
   end subroutine spiral_case

   subroutine cubic_topology(ncell,width,topology,cell,a_in)
      integer, intent(in) :: ncell, width
      type(block_topology_type), intent(out) :: topology
      real(dblprec), intent(out) :: cell(3,3)
      real(dblprec), intent(in), optional :: a_in
      integer :: status
      character(len=512) :: message
      real(dblprec) :: a

      a = 1.0d-10
      if (present(a_in)) a = a_in
      cell = 0.0_dblprec
      cell(1,1) = a
      cell(2,2) = a
      cell(3,3) = a
      call build_block_topology(topology,REGULAR_REPLICATED_CELL,1,(/ncell,1,1/), &
         ncell,(/width,1,1/),cell,(/1/),status,message)
      call check(status == BLOCK_TOPOLOGY_OK,'test topology builds: '//trim(message))
   end subroutine cubic_topology

   subroutine make_material(cell,exchange,dmi,material)
      real(dblprec), intent(in) :: cell(3,3), exchange(3,3), dmi(3,3)
      type(coarse_material_type), intent(out) :: material

      material%ready = .true.
      material%n_basis = 1
      material%n_channels = 1
      material%cell_volume_m3 = abs(cell(1,1)*cell(2,2)*cell(3,3))
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
   end subroutine make_material

   character(len=16) function integer_text(value)
      integer, intent(in) :: value
      write(integer_text,'(i0)') value
   end function integer_text

   subroutine check(condition,message)
      logical, intent(in) :: condition
      character(len=*), intent(in) :: message

      if (.not. condition) then
         failures = failures+1
         write(*,'(a)') 'FAIL: '//trim(message)
      end if
   end subroutine check

   subroutine check_close(actual,expected,relative_tolerance,message)
      real(dblprec), intent(in) :: actual, expected, relative_tolerance
      character(len=*), intent(in) :: message
      real(dblprec) :: scale

      scale = max(abs(expected),abs(actual),tiny(1.0_dblprec))
      call check(abs(actual-expected) <= relative_tolerance*scale,message)
      if (abs(actual-expected) > relative_tolerance*scale) then
         write(*,'(a,2es24.15)') '  actual/expected: ',actual,expected
      end if
   end subroutine check_close

end program test_smooth_projected_operator
