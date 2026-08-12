program test_coarse_tensor_operator

   use Parameters, only : dblprec
   use BlockTopology
   use stiffness
   use CoarseTensorOperator

   implicit none

   integer :: failures

   failures = 0
   call test_uniform_rotation_and_metric()
   call test_energy_derivatives_and_reporting()
   call test_small_q_chirality_and_block_one()
   call test_small_q_chiral_minimum()
   call test_domain_wall_refinement()
   call test_llg_and_setup_rejections()
   call test_dipole_unmasked_and_exactly_once()

   if (failures /= 0) then
      write(*,'(a,i0)') 'coarse tensor operator tests failed: ', failures
      stop 1
   end if
   write(*,'(a)') 'coarse tensor operator tests passed'

contains

   subroutine test_uniform_rotation_and_metric()
      type(block_topology_type) :: topology
      type(coarse_material_type) :: material
      type(coarse_tensor_operator_type) :: operator
      type(coarse_operator_options_type) :: options
      type(coarse_energy_terms_type) :: energy, rotated_energy
      type(coarse_field_terms_type) :: terms
      integer :: status, block, n
      character(len=512) :: message
      real(dblprec), parameter :: a = 2.0d-10
      real(dblprec), parameter :: exchange_a = 8.0d-12
      real(dblprec) :: cell(3,3), exchange(3,3), dmi(3,3), angle, expected
      real(dblprec) :: rotation(3,3), skew(3,3), inverse_transpose(3,3)
      real(dblprec) :: delta(3), discrete_gradient(3,3), phase_step(3)
      real(dblprec), allocatable :: direction(:,:), rotated(:,:), field(:,:), zero_field(:,:)

      cell = 0.0_dblprec
      cell(1,1) = a
      cell(2,2) = a
      cell(3,3) = a
      exchange = 0.0_dblprec
      exchange(1,1) = exchange_a
      exchange(2,2) = 0.7_dblprec*exchange_a
      exchange(3,3) = 1.3_dblprec*exchange_a
      dmi = 0.0_dblprec
      dmi(3,1) = 1.0d-3
      call make_fixture((/8,3,1/),(/1,1,1/),cell,exchange,dmi,topology,material,operator,options)

      allocate(direction(3,operator%nblocks),rotated(3,operator%nblocks), &
         field(3,operator%nblocks),zero_field(3,operator%nblocks))
      direction = 0.0_dblprec
      direction(1,:) = 1.0_dblprec
      zero_field = 0.0_dblprec
      call evaluate_coarse_tensor_operator(operator,direction,field,energy,status,message, &
         external_field_t=zero_field,term_fields=terms)
      call check(status == COARSE_TENSOR_OK,'uniform state evaluates')
      call check(maxval(abs(terms%exchange_t))+maxval(abs(terms%spiralization_t)) < 1.0d-13, &
         'uniform state has zero exchange and DMI fields')
      call check(abs(energy%exchange_j)+abs(energy%spiralization_j) < 1.0d-40, &
         'uniform state has zero exchange and DMI energies')

      do block = 1, operator%nblocks
         angle = 2.0_dblprec*pi*real(mod(block-1,8),dblprec)/8.0_dblprec + &
            0.17_dblprec*real((block-1)/8,dblprec)
         direction(:,block) = (/cos(angle),0.3_dblprec*sin(angle),sin(angle)/)
         direction(:,block) = direction(:,block)/vector_norm(direction(:,block))
      end do
      rotation = reshape((/0.0_dblprec,1.0_dblprec,0.0_dblprec, &
                          0.0_dblprec,0.0_dblprec,1.0_dblprec, &
                          1.0_dblprec,0.0_dblprec,0.0_dblprec/),(/3,3/))
      rotated = matmul(rotation,direction)
      call evaluate_coarse_tensor_operator(operator,direction,field,energy,status,message, &
         external_field_t=zero_field)
      call evaluate_coarse_tensor_operator(operator,rotated,field,rotated_energy,status,message, &
         external_field_t=zero_field)
      call check_close(rotated_energy%exchange_j,energy%exchange_j,2.0d-13, &
         'global spin rotation leaves tensor exchange energy invariant')

      skew(:,1) = (/a,0.0_dblprec,0.0_dblprec/)
      skew(:,2) = (/0.35_dblprec*a,0.9_dblprec*a,0.0_dblprec/)
      skew(:,3) = (/0.15_dblprec*a,0.2_dblprec*a,1.1_dblprec*a/)
      exchange = 0.0_dblprec
      exchange(1,1) = exchange_a
      exchange(2,2) = exchange_a
      exchange(3,3) = exchange_a
      dmi = 0.0_dblprec
      call make_fixture((/8,3,2/),(/1,1,1/),skew,exchange,dmi,topology,material,operator,options)
      n = operator%nblocks
      call allocate_if_needed(direction,n,rotated,field,zero_field)
      phase_step = (/2.0_dblprec*pi/8.0_dblprec,2.0_dblprec*pi/3.0_dblprec,0.0_dblprec/)
      do block = 1, n
         angle = dot_product(phase_step,real(operator%block_coordinate(:,block),dblprec))
         direction(:,block) = (/cos(angle),sin(angle),0.0_dblprec/)
      end do
      zero_field = 0.0_dblprec
      call evaluate_coarse_tensor_operator(operator,direction,field,energy,status,message, &
         external_field_t=zero_field)
      inverse_transpose = transpose(inverse_matrix(skew))
      discrete_gradient = 0.0_dblprec
      do block = 1, 3
         delta = (/cos(phase_step(block))-1.0_dblprec,sin(phase_step(block)),0.0_dblprec/)
         discrete_gradient = discrete_gradient + spread(delta,2,3) * &
            spread(inverse_transpose(:,block),1,3)
      end do
      expected = real(n,dblprec)*abs(det3(skew))*exchange_a*sum(discrete_gradient**2)
      call check_close(energy%exchange_j,expected,3.0d-13, &
         'skew-cell stencil applies the physical inverse-transpose metric')
   end subroutine test_uniform_rotation_and_metric

   subroutine test_energy_derivatives_and_reporting()
      type(block_topology_type) :: topology
      type(coarse_material_type) :: material
      type(coarse_tensor_operator_type) :: operator
      type(coarse_operator_options_type) :: options
      type(coarse_anisotropy_type) :: anisotropy
      type(coarse_energy_terms_type) :: energy
      type(coarse_field_terms_type) :: terms
      integer :: status, block, term
      character(len=512) :: message
      real(dblprec), parameter :: a = 2.5d-10
      real(dblprec) :: cell(3,3), exchange(3,3), dmi(3,3)
      real(dblprec) :: direction(3,12), field(3,12), external(3,12), dipole(3,12)
      real(dblprec) :: tangent(3), numerical, analytic

      cell = 0.0_dblprec
      cell(1,1) = a
      cell(2,2) = a
      cell(3,3) = a
      exchange = reshape((/9.0d-12,1.0d-12,0.0d0, &
                           1.0d-12,7.0d-12,0.5d-12, &
                           0.0d0,0.5d-12,5.0d-12/),(/3,3/))
      dmi = reshape((/0.2d-3,-0.1d-3,0.3d-3, &
                      0.4d-3,0.2d-3,-0.2d-3, &
                      0.8d-3,-0.3d-3,0.1d-3/),(/3,3/))
      options%use_uniform_coarse_dipole = .true.
      call make_topology_and_material((/6,2,1/),(/1,1,1/),cell,exchange,dmi,topology,material)
      call make_anisotropy(topology%n_spatial_blocks,anisotropy,-3.0d5,0.7d5)
      anisotropy%axis_count = 2
      anisotropy%axis(1,2,:) = 1.0_dblprec
      anisotropy%k1(2,:) = 0.8d5
      anisotropy%k2(2,:) = -0.2d5
      call setup_coarse_tensor_operator(operator,topology,material,options,status,message,anisotropy)
      call check(status == COARSE_TENSOR_OK,'full supported term fixture sets up')
      if (status /= COARSE_TENSOR_OK) return

      do block = 1, operator%nblocks
         direction(:,block) = (/0.4_dblprec+0.03_dblprec*block, &
            sin(0.37_dblprec*block),0.7_dblprec-0.02_dblprec*block/)
         direction(:,block) = direction(:,block)/vector_norm(direction(:,block))
         external(:,block) = (/0.1_dblprec,-0.03_dblprec,0.02_dblprec/) * &
            (1.0_dblprec+0.01_dblprec*block)
      end do
      call local_symmetric_dipole(direction,dipole)
      call evaluate_coarse_tensor_operator(operator,direction,field,energy,status,message, &
         external,dipole,terms)
      call check(status == COARSE_TENSOR_OK,'all supported terms evaluate together')
      call check_close(energy%total_j,energy%exchange_j+energy%spiralization_j+ &
         energy%anisotropy_j+energy%external_j+energy%dipole_j,1.0d-14, &
         'per-term energy report sums exactly to total')
      call check(maxval(abs(field-(terms%exchange_t+terms%spiralization_t+ &
         terms%anisotropy_t+terms%external_t+terms%dipole_t))) < &
         1.0d-13*max(1.0_dblprec,maxval(abs(field))), &
         'per-term fields sum exactly to total field')

      tangent = (/0.3_dblprec,-0.7_dblprec,0.2_dblprec/)
      tangent = tangent-dot_product(tangent,direction(:,3))*direction(:,3)
      tangent = tangent/vector_norm(tangent)
      do term = 1, 5
         call finite_energy_derivative(operator,direction,external,3,tangent,term,numerical)
         select case (term)
         case (1)
            analytic = -COARSE_MUB_SI*operator%block_moment_mub*dot_product(terms%exchange_t(:,3),tangent)
         case (2)
            analytic = -COARSE_MUB_SI*operator%block_moment_mub*dot_product(terms%spiralization_t(:,3),tangent)
         case (3)
            analytic = -COARSE_MUB_SI*operator%block_moment_mub*dot_product(terms%anisotropy_t(:,3),tangent)
         case (4)
            analytic = -COARSE_MUB_SI*operator%block_moment_mub*dot_product(terms%external_t(:,3),tangent)
         case default
            analytic = -COARSE_MUB_SI*operator%block_moment_mub*dot_product(terms%dipole_t(:,3),tangent)
         end select
         call check_close(numerical,analytic,2.0d-6,'energy derivative matches field for term '//digit(term))
      end do
   end subroutine test_energy_derivatives_and_reporting

   subroutine test_small_q_chirality_and_block_one()
      type(block_topology_type) :: topology
      type(coarse_material_type) :: material
      type(coarse_tensor_operator_type) :: operator
      type(coarse_operator_options_type) :: options
      type(coarse_energy_terms_type) :: positive, negative
      type(coarse_field_terms_type) :: terms
      integer :: status, block, block_width
      character(len=512) :: message
      real(dblprec), parameter :: a = 2.0d-10
      real(dblprec), parameter :: pair_j = 4.0d-21
      real(dblprec), parameter :: gamma = 3.2d11
      real(dblprec) :: cell(3,3), exchange(3,3), dmi(3,3), q, angle
      real(dblprec), allocatable :: direction(:,:), field(:,:), zero(:,:)
      real(dblprec) :: atomistic_frequency, coarse_frequency, torque(3), relative_error

      cell = 0.0_dblprec
      cell(1,1) = a
      cell(2,2) = a
      cell(3,3) = a
      exchange = 0.0_dblprec
      exchange(1,1) = pair_j/(2.0_dblprec*a)
      dmi = 0.0_dblprec
      dmi(3,1) = 1.5d-3
      q = 2.0_dblprec*pi/(64.0_dblprec*a)

      do block_width = 1, 4, 3
         call make_fixture((/64,1,1/),(/block_width,1,1/),cell,exchange,dmi, &
            topology,material,operator,options)
         material%channel_gamma(1) = gamma
         call setup_coarse_tensor_operator(operator,topology,material,options,status,message)
         allocate(direction(3,operator%nblocks),field(3,operator%nblocks),zero(3,operator%nblocks))
         zero = 0.0_dblprec
         direction = 0.0_dblprec
         do block = 1, operator%nblocks
            angle = q*operator%block_vectors_m(1,1)*real(operator%block_coordinate(1,block),dblprec)
            direction(1,block) = 1.0d-6*cos(angle)
            direction(3,block) = sqrt(1.0_dblprec-direction(1,block)**2)
         end do
         call evaluate_coarse_tensor_operator(operator,direction,field,positive,status,message, &
            external_field_t=zero,term_fields=terms)
         torque = cross_product(direction(:,1),terms%exchange_t(:,1))
         coarse_frequency = gamma*abs(torque(2))/1.0d-6
         atomistic_frequency = gamma*2.0_dblprec*pair_j*(1.0_dblprec-cos(q*a))/ &
            (COARSE_MUB_SI*material%channel_moment_mub(1))
         relative_error = abs(coarse_frequency-atomistic_frequency)/atomistic_frequency
         if (block_width == 1) then
            call check(relative_error < 2.0d-8, &
               'block-size-one exchange dispersion equals nearest-neighbour atomistic result')
         else
            call check(relative_error < 2.0d-2, &
               'small-q coarse exchange dispersion matches atomistic reference')
         end if
         deallocate(direction,field,zero)
      end do

      call make_fixture((/64,1,1/),(/2,1,1/),cell,exchange,dmi,topology,material,operator,options)
      allocate(direction(3,operator%nblocks),field(3,operator%nblocks),zero(3,operator%nblocks))
      zero = 0.0_dblprec
      do block = 1, operator%nblocks
         angle = q*operator%block_vectors_m(1,1)*real(operator%block_coordinate(1,block),dblprec)
         direction(:,block) = (/cos(angle),sin(angle),0.0_dblprec/)
      end do
      call evaluate_coarse_tensor_operator(operator,direction,field,positive,status,message, &
         external_field_t=zero)
      direction(2,:) = -direction(2,:)
      call evaluate_coarse_tensor_operator(operator,direction,field,negative,status,message, &
         external_field_t=zero)
      call check(positive%spiralization_j > 0.0_dblprec .and. &
         negative%spiralization_j < 0.0_dblprec, &
         'positive D_zx gives E(+q)>E(-q) and the approved negative-q chirality')
   end subroutine test_small_q_chirality_and_block_one

   subroutine test_small_q_chiral_minimum()
      ! DMI-SPIRAL-Q: for E/V=A q^2+D_zx q, the independently prescribed
      ! continuum minimum is q_min=-D_zx/(2A).  The same physical chain is
      ! sampled at two block resolutions; both must select the left-handed
      ! (negative-q) mode and approach the analytic magnitude.
      type(block_topology_type) :: topology
      type(coarse_material_type) :: material
      type(coarse_tensor_operator_type) :: operator
      type(coarse_operator_options_type) :: options
      type(coarse_energy_terms_type) :: energy, minimum
      integer, parameter :: ncell = 256
      integer :: width, mode, best_mode, block, status
      character(len=512) :: message
      real(dblprec), parameter :: a = 2.0d-10, exchange_a = 1.0d-11
      real(dblprec) :: cell(3,3), exchange(3,3), dmi(3,3), dmi_zx
      real(dblprec) :: length, q, q_minimum, q_analytic, angle
      real(dblprec), allocatable :: direction(:,:), field(:,:), zero(:,:)

      cell = 0.0_dblprec
      cell(1,1) = a
      cell(2,2) = a
      cell(3,3) = a
      exchange = 0.0_dblprec
      exchange(1,1) = exchange_a
      length = real(ncell,dblprec)*a
      q_analytic = -2.0_dblprec*acos(-1.0_dblprec)/length
      dmi_zx = -2.0_dblprec*exchange_a*q_analytic
      dmi = 0.0_dblprec
      dmi(3,1) = dmi_zx

      do width = 1, 4, 3
         call make_fixture((/ncell,1,1/),(/width,1,1/),cell,exchange,dmi, &
            topology,material,operator,options)
         allocate(direction(3,operator%nblocks),field(3,operator%nblocks),zero(3,operator%nblocks))
         zero = 0.0_dblprec
         minimum%total_j = huge(1.0_dblprec)
         best_mode = 0
         do mode = -3,3
            q = 2.0_dblprec*acos(-1.0_dblprec)*real(mode,dblprec)/length
            do block = 1,operator%nblocks
               angle = q*operator%block_vectors_m(1,1)*real(operator%block_coordinate(1,block),dblprec)
               direction(:,block) = (/cos(angle),sin(angle),0.0_dblprec/)
            end do
            call evaluate_coarse_tensor_operator(operator,direction,field,energy,status,message, &
               external_field_t=zero)
            call check(status == COARSE_TENSOR_OK,'small-q chiral chain evaluates')
            if (energy%total_j < minimum%total_j) then
               minimum = energy
               best_mode = mode
            end if
         end do
         q_minimum = 2.0_dblprec*acos(-1.0_dblprec)*real(best_mode,dblprec)/length
         call check(best_mode < 0,'positive D_zx selects the analytic left-handed (negative-q) minimum')
         call check_close(q_minimum,q_analytic,2.0d-1, &
            'small-q chiral-chain minimum magnitude matches -D_zx/(2A)')
         deallocate(direction,field,zero)
      end do
   end subroutine test_small_q_chiral_minimum

   subroutine test_domain_wall_refinement()
      type(block_topology_type) :: topology
      type(coarse_material_type) :: material
      type(coarse_tensor_operator_type) :: operator
      type(coarse_operator_options_type) :: options
      type(coarse_anisotropy_type) :: anisotropy
      type(coarse_energy_terms_type) :: energy
      integer :: status, refinement, block, width_cells
      character(len=512) :: message
      real(dblprec), parameter :: a = 0.5d-9
      real(dblprec), parameter :: exchange_a = 1.0d-11
      real(dblprec), parameter :: easy_k = 4.0d5
      real(dblprec) :: cell(3,3), exchange(3,3), dmi(3,3), wall_delta, length
      real(dblprec) :: x, x1, x2, theta, expected, baseline, numerical_width
      real(dblprec) :: energy_error(4), width_error(4)
      real(dblprec), allocatable :: direction(:,:), field(:,:), zero(:,:), sampled_theta(:)

      cell = 0.0_dblprec
      cell(1,1) = a
      cell(2,2) = a
      cell(3,3) = a
      exchange = 0.0_dblprec
      exchange(1,1) = exchange_a
      dmi = 0.0_dblprec
      wall_delta = sqrt(exchange_a/easy_k)
      length = 256.0_dblprec*a
      x1 = 0.25_dblprec*length
      x2 = 0.75_dblprec*length
      expected = 8.0_dblprec*sqrt(exchange_a*easy_k)*a*a

      do refinement = 1, 4
         width_cells = 2**(4-refinement)
         call make_topology_and_material((/256,1,1/),(/width_cells,1,1/),cell, &
            exchange,dmi,topology,material)
         call make_anisotropy(topology%n_spatial_blocks,anisotropy,-easy_k,0.0_dblprec)
         call setup_coarse_tensor_operator(operator,topology,material,options,status,message,anisotropy)
         allocate(direction(3,operator%nblocks),field(3,operator%nblocks), &
            zero(3,operator%nblocks),sampled_theta(operator%nblocks))
         zero = 0.0_dblprec
         do block = 1, operator%nblocks
            x = (real(operator%block_coordinate(1,block),dblprec)+0.5_dblprec) * &
               operator%block_vectors_m(1,1)
            theta = wall_angle(x,x1,wall_delta)-wall_angle(x,x2,wall_delta)
            sampled_theta(block) = theta
            direction(:,block) = (/sin(theta),0.0_dblprec,cos(theta)/)
         end do
         call evaluate_coarse_tensor_operator(operator,direction,field,energy,status,message, &
            external_field_t=zero)
         baseline = -easy_k*operator%block_volume_m3*real(operator%nblocks,dblprec)
         energy_error(refinement) = abs((energy%exchange_j+energy%anisotropy_j-baseline)-expected)/expected
         block = nint(x1/operator%block_vectors_m(1,1))
         block = max(1,min(operator%nblocks-1,block))
         numerical_width = operator%block_vectors_m(1,1) / &
            abs(sampled_theta(block+1)-sampled_theta(block))
         width_error(refinement) = abs(numerical_width-wall_delta)/wall_delta
         deallocate(direction,field,zero,sampled_theta)
      end do
      call check(all(energy_error(2:4) < energy_error(1:3)) .and. energy_error(4) < 2.0d-3, &
         'two-wall energy converges to 8 sqrt(A K) under block refinement')
      call check(all(width_error(2:4) < width_error(1:3)) .and. width_error(4) < 2.0d-3, &
         'sampled domain-wall width converges to sqrt(A/K) under block refinement')
   end subroutine test_domain_wall_refinement

   subroutine test_llg_and_setup_rejections()
      type(block_topology_type) :: topology, bad_topology
      type(coarse_material_type) :: material, bad_material
      type(coarse_tensor_operator_type) :: operator
      type(coarse_operator_options_type) :: options
      type(coarse_anisotropy_type) :: anisotropy
      integer :: status
      character(len=512) :: message
      real(dblprec) :: cell(3,3), exchange(3,3), dmi(3,3)
      real(dblprec) :: direction(3,2), field(3,2), rhs(3,2), expected(3)

      cell = 0.0_dblprec
      cell(1,1) = 2.0d-10
      cell(2,2) = 2.0d-10
      cell(3,3) = 2.0d-10
      exchange = 0.0_dblprec
      dmi = 0.0_dblprec
      call make_topology_and_material((/2,1,1/),(/1,1,1/),cell,exchange,dmi,topology,material)
      call setup_coarse_tensor_operator(operator,topology,material,options,status,message)
      direction = 0.0_dblprec
      direction(1,:) = 1.0_dblprec
      field = 0.0_dblprec
      field(3,:) = 2.0_dblprec
      call coarse_llg_rhs(operator,direction,field,rhs,status,message)
      expected = operator%channel_gamma_per_t_s*2.0_dblprec / &
         (1.0_dblprec+operator%channel_damping**2) * &
         (/0.0_dblprec,1.0_dblprec,operator%channel_damping/)
      call check(maxval(abs(rhs(:,1)-expected)) < 1.0d-14*maxval(abs(expected)), &
         'coarse deterministic LLG algebra matches the accepted Gilbert sign and damping')

      options%temperature_k = 1.0_dblprec
      call setup_coarse_tensor_operator(operator,topology,material,options,status,message)
      call check(status == COARSE_TENSOR_UNSUPPORTED_TEMPERATURE .and. index(message,'T=0') > 0, &
         'finite temperature rejects cleanly during setup')
      options = coarse_operator_options_type()
      options%solver_mode = 99
      call setup_coarse_tensor_operator(operator,topology,material,options,status,message)
      call check(status == COARSE_TENSOR_UNSUPPORTED_SOLVER .and. index(message,'Heun') > 0, &
         'unsupported solver rejects cleanly during setup')
      options = coarse_operator_options_type()
      options%boundary_mode(1) = 99
      call setup_coarse_tensor_operator(operator,topology,material,options,status,message)
      call check(status == COARSE_TENSOR_UNSUPPORTED_GEOMETRY .and. index(message,'periodic') > 0, &
         'unsupported boundary geometry rejects cleanly during setup')
      options = coarse_operator_options_type()
      options%adaptive_switching = .true.
      call setup_coarse_tensor_operator(operator,topology,material,options,status,message)
      call check(status == COARSE_TENSOR_UNSUPPORTED_MODE .and. index(message,'static all-coarse') > 0, &
         'adaptive behavior is explicitly rejected rather than hidden in CG-04')
      options = coarse_operator_options_type()
      options%interface_coupling = .true.
      call setup_coarse_tensor_operator(operator,topology,material,options,status,message)
      call check(status == COARSE_TENSOR_UNSUPPORTED_MODE, &
         'interface behavior is explicitly rejected rather than hidden in CG-04')

      bad_material = material
      bad_material%n_channels = 2
      call setup_coarse_tensor_operator(operator,topology,bad_material, &
         coarse_operator_options_type(),status,message)
      call check(status == COARSE_TENSOR_INVALID_MATERIAL .and. index(message,'exactly one') > 0, &
         'unsupported multichannel material rejects through the accepted runtime gate')
      bad_topology = topology
      bad_topology%geometry_mode = EXPLICIT_DEVICE
      call setup_coarse_tensor_operator(operator,bad_topology,material, &
         coarse_operator_options_type(),status,message)
      call check(status == COARSE_TENSOR_UNSUPPORTED_GEOMETRY, &
         'unsupported explicit-device geometry rejects during setup')

      call make_anisotropy(topology%n_spatial_blocks,anisotropy,-1.0d5,0.0_dblprec)
      anisotropy%kind(1) = 99
      call setup_coarse_tensor_operator(operator,topology,material, &
         coarse_operator_options_type(),status,message,anisotropy)
      call check(status == COARSE_TENSOR_INVALID_ANISOTROPY .and. index(message,'uniaxial') > 0, &
         'unsupported anisotropy kind rejects cleanly during setup')
   end subroutine test_llg_and_setup_rejections

   !> RCG-05F: the CPU-side counterpart to
   !> tests/coarse_graining/test_gpu_adaptive_runtime.cpp's GPU-only "FFT
   !> dipole included exactly once" assertions (lines 367,407,414,421 there).
   !> Confirms, on a genuinely anisotropic-block-shape, non-orthogonal (skew)
   !> cell fixture (not the cubic ones the module docstring at
   !> coarsetensoroperator.f90:345-347 was presumably last checked against),
   !> that the uniform coarse FFT dipole (coarsetensoroperator.f90:507-513)
   !> is (1) never gated by interaction_owner/onsite_owner -- the "deliberately
   !> not masked here" claim -- (2) included for every block, and (3) equals
   !> the analytic all-grid sum exactly, i.e. neither doubled nor dropped
   !> anywhere, and does not perturb the masked exchange/DMI/anisotropy/
   !> external terms it is summed alongside.
   subroutine test_dipole_unmasked_and_exactly_once()
      type(block_topology_type) :: topology
      type(coarse_material_type) :: material
      type(coarse_tensor_operator_type) :: operator
      type(coarse_operator_options_type) :: options
      type(coarse_anisotropy_type) :: anisotropy
      type(coarse_energy_terms_type) :: energy_all, energy_none, energy_mixed, energy_mixed_nodip
      type(coarse_field_terms_type) :: terms_all, terms_none, terms_mixed
      integer :: status, block, n
      character(len=512) :: message
      real(dblprec), parameter :: a = 2.0d-10
      real(dblprec) :: skew(3,3), exchange(3,3), dmi(3,3), expected_dipole_j
      real(dblprec), allocatable :: direction(:,:), field(:,:), external(:,:), dipole(:,:)
      real(dblprec), allocatable :: zero_dipole(:,:), field_mixed(:,:)
      logical, allocatable :: owner_all(:), owner_none(:), owner_mixed(:)

      ! Genuinely non-orthogonal (skew) atom cell, stretched anisotropically
      ! into blocks (block_shape=(1,2,3), like RCG-05B's unequal-width
      ! fixtures): a1 is left unskewed, a2/a3 carry off-diagonal components,
      ! matching the same isolate-one-variable-at-a-time convention
      ! tests/coarse_graining/e2e/ownership_dipole_skew/README.md documents.
      skew(:,1) = (/a,0.0_dblprec,0.0_dblprec/)
      skew(:,2) = (/0.3_dblprec*a,0.85_dblprec*a,0.0_dblprec/)
      skew(:,3) = (/0.1_dblprec*a,0.15_dblprec*a,1.2_dblprec*a/)
      exchange = 0.0_dblprec
      exchange(1,1) = 8.0d-12
      exchange(2,2) = 6.0d-12
      exchange(3,3) = 4.0d-12
      dmi = 0.0_dblprec
      dmi(3,1) = 0.5d-3

      options%use_uniform_coarse_dipole = .true.
      call make_topology_and_material((/6,4,6/),(/1,2,3/),skew,exchange,dmi,topology,material)
      call make_anisotropy(topology%n_spatial_blocks,anisotropy,-2.0d5,0.5d5)
      call setup_coarse_tensor_operator(operator,topology,material,options,status,message,anisotropy)
      call check(status == COARSE_TENSOR_OK,'dipole exactly-once fixture sets up: '//trim(message))
      if (status /= COARSE_TENSOR_OK) return

      n = operator%nblocks
      allocate(direction(3,n),field(3,n),external(3,n),dipole(3,n),zero_dipole(3,n),field_mixed(3,n))
      zero_dipole = 0.0_dblprec
      allocate(owner_all(n),owner_none(n),owner_mixed(n))
      do block = 1, n
         direction(:,block) = (/0.5_dblprec+0.02_dblprec*block,sin(0.29_dblprec*block), &
            0.6_dblprec-0.01_dblprec*block/)
         direction(:,block) = direction(:,block)/vector_norm(direction(:,block))
         external(:,block) = (/0.05_dblprec,-0.02_dblprec,0.01_dblprec/)
         dipole(:,block) = (/0.03_dblprec,0.07_dblprec,-0.04_dblprec/) * &
            (1.0_dblprec+0.005_dblprec*block)
      end do
      owner_all = .true.
      owner_none = .false.
      owner_mixed = .true.
      do block = 1, n, 3
         owner_mixed(block) = .false.
      end do

      call evaluate_coarse_tensor_operator(operator,direction,field,energy_all,status,message, &
         external_field_t=external,uniform_coarse_dipole_field_t=dipole,term_fields=terms_all, &
         interaction_owner=owner_all,onsite_owner=owner_all)
      call check(status == COARSE_TENSOR_OK,'all-owned + dipole evaluates: '//trim(message))
      call evaluate_coarse_tensor_operator(operator,direction,field,energy_none,status,message, &
         external_field_t=external,uniform_coarse_dipole_field_t=dipole,term_fields=terms_none, &
         interaction_owner=owner_none,onsite_owner=owner_none)
      call check(status == COARSE_TENSOR_OK,'none-owned + dipole evaluates: '//trim(message))
      call evaluate_coarse_tensor_operator(operator,direction,field_mixed,energy_mixed,status,message, &
         external_field_t=external,uniform_coarse_dipole_field_t=dipole,term_fields=terms_mixed, &
         interaction_owner=owner_mixed,onsite_owner=owner_mixed)
      call check(status == COARSE_TENSOR_OK,'mixed-owned + dipole evaluates: '//trim(message))

      ! (1) Never gated by the mask: dipole energy/field identical whether
      ! every block, no block, or a mixed subset is "owned".
      call check_close(energy_none%dipole_j,energy_all%dipole_j,1.0d-13, &
         'dipole energy is independent of interaction_owner/onsite_owner (none vs all)')
      call check_close(energy_mixed%dipole_j,energy_all%dipole_j,1.0d-13, &
         'dipole energy is independent of interaction_owner/onsite_owner (mixed vs all)')
      call check(maxval(abs(terms_none%dipole_t-terms_all%dipole_t)) < &
         1.0d-13*max(1.0_dblprec,maxval(abs(terms_all%dipole_t))), &
         'dipole field is independent of interaction_owner/onsite_owner (none vs all)')
      call check(maxval(abs(terms_mixed%dipole_t-terms_all%dipole_t)) < &
         1.0d-13*max(1.0_dblprec,maxval(abs(terms_all%dipole_t))), &
         'dipole field is independent of interaction_owner/onsite_owner (mixed vs all)')
      call check(maxval(abs(terms_none%dipole_t-dipole)) < &
         1.0d-13*max(1.0_dblprec,maxval(abs(dipole))), &
         'dipole field term equals the supplied all-grid field exactly for every block')

      ! (2)+(3) Exactly once, all-grid: equals the analytic sum over EVERY
      ! block (not just owned ones), not zero, not doubled.
      expected_dipole_j = 0.0_dblprec
      do block = 1, n
         expected_dipole_j = expected_dipole_j - 0.5_dblprec*COARSE_MUB_SI* &
            operator%block_moment_mub*dot_product(dipole(:,block),direction(:,block))
      end do
      call check_close(energy_all%dipole_j,expected_dipole_j,1.0d-12, &
         'dipole energy equals the analytic all-grid sum exactly once (all-owned)')
      call check_close(energy_none%dipole_j,expected_dipole_j,1.0d-12, &
         'dipole energy equals the analytic all-grid sum exactly once (none-owned)')

      ! Does not perturb the masked terms it is summed alongside. Since this
      ! operator was set up with use_uniform_coarse_dipole=.true., the
      ! dipole argument cannot be omitted (coarsetensoroperator.f90:414
      ! rejects a presence mismatch against setup capability) -- so "no
      ! dipole" is exercised as a genuinely zero-valued dipole field through
      ! the exact same code path, isolating the VALUE's effect rather than
      ! a different code path.
      call evaluate_coarse_tensor_operator(operator,direction,field,energy_mixed_nodip,status,message, &
         external_field_t=external,uniform_coarse_dipole_field_t=zero_dipole, &
         interaction_owner=owner_mixed,onsite_owner=owner_mixed)
      call check(status == COARSE_TENSOR_OK,'mixed-owned with zero dipole evaluates: '//trim(message))
      call check_close(energy_mixed%exchange_j,energy_mixed_nodip%exchange_j,1.0d-13, &
         'adding a nonzero dipole does not perturb the masked exchange energy')
      call check_close(energy_mixed%spiralization_j,energy_mixed_nodip%spiralization_j,1.0d-13, &
         'adding a nonzero dipole does not perturb the masked DMI energy')
      call check_close(energy_mixed%anisotropy_j,energy_mixed_nodip%anisotropy_j,1.0d-13, &
         'adding a nonzero dipole does not perturb the masked anisotropy energy')
      call check_close(energy_mixed%external_j,energy_mixed_nodip%external_j,1.0d-13, &
         'adding a nonzero dipole does not perturb the masked external-field energy')
      call check(energy_mixed_nodip%dipole_j == 0.0_dblprec, &
         'a zero-valued dipole field contributes exactly zero dipole energy')

      ! Per-term energies/fields (including dipole) still sum exactly to the
      ! total -- the same additivity contract test_energy_derivatives_and_
      ! reporting already checks without a mask, now confirmed under one.
      call check_close(energy_mixed%total_j, &
         energy_mixed%exchange_j+energy_mixed%spiralization_j+energy_mixed%anisotropy_j+ &
         energy_mixed%external_j+energy_mixed%dipole_j,1.0d-13, &
         'per-term energies (including dipole) sum exactly to the masked total')
      call check(maxval(abs(field_mixed-(terms_mixed%exchange_t+terms_mixed%spiralization_t+ &
         terms_mixed%anisotropy_t+terms_mixed%external_t+terms_mixed%dipole_t))) < &
         1.0d-13*max(1.0_dblprec,maxval(abs(field_mixed))), &
         'per-term fields (including dipole) sum exactly to the masked total field')
   end subroutine test_dipole_unmasked_and_exactly_once

   subroutine make_fixture(repetitions,block_shape,cell,exchange,dmi,topology,material,operator,options)
      integer, intent(in) :: repetitions(3), block_shape(3)
      real(dblprec), intent(in) :: cell(3,3), exchange(3,3), dmi(3,3)
      type(block_topology_type), intent(out) :: topology
      type(coarse_material_type), intent(out) :: material
      type(coarse_tensor_operator_type), intent(out) :: operator
      type(coarse_operator_options_type), intent(in) :: options
      integer :: status
      character(len=512) :: message

      call make_topology_and_material(repetitions,block_shape,cell,exchange,dmi,topology,material)
      call setup_coarse_tensor_operator(operator,topology,material,options,status,message)
      call check(status == COARSE_TENSOR_OK,'reference operator setup: '//trim(message))
   end subroutine make_fixture

   subroutine make_topology_and_material(repetitions,block_shape,cell,exchange,dmi,topology,material)
      integer, intent(in) :: repetitions(3), block_shape(3)
      real(dblprec), intent(in) :: cell(3,3), exchange(3,3), dmi(3,3)
      type(block_topology_type), intent(out) :: topology
      type(coarse_material_type), intent(out) :: material
      integer :: status, natom
      character(len=512) :: message

      natom = product(repetitions)
      call build_block_topology(topology,REGULAR_REPLICATED_CELL,1,repetitions,natom, &
         block_shape,cell,(/1/),status,message)
      call check(status == BLOCK_TOPOLOGY_OK,'fixture topology builds: '//trim(message))

      material%ready = .true.
      material%n_basis = 1
      material%n_channels = 1
      material%cell_volume_m3 = abs(det3(cell))
      allocate(material%basis_to_channel(1),material%basis_moment_mub(1), &
         material%channel_moment_mub(1),material%channel_gamma(1), &
         material%channel_damping(1),material%exchange_stiffness(3,3,1,1), &
         material%spiralization(3,3,1,1))
      material%basis_to_channel = 1
      material%basis_moment_mub = 1.7_dblprec
      material%channel_moment_mub = 1.7_dblprec
      material%channel_gamma = 3.2d11
      material%channel_damping = 0.07_dblprec
      material%exchange_stiffness(:,:,1,1) = exchange
      material%spiralization(:,:,1,1) = dmi
      material%diagnostics%channel_dynamics_parameters_validated = .true.
      material%diagnostics%small_q_energy_validated = .true.
   end subroutine make_topology_and_material

   subroutine make_anisotropy(nblocks,anisotropy,k1,k2)
      integer, intent(in) :: nblocks
      real(dblprec), intent(in) :: k1, k2
      type(coarse_anisotropy_type), intent(out) :: anisotropy

      allocate(anisotropy%kind(nblocks),anisotropy%axis_count(nblocks), &
         anisotropy%axis(3,2,nblocks),anisotropy%k1(2,nblocks),anisotropy%k2(2,nblocks))
      anisotropy%kind = COARSE_ANISOTROPY_UNIAXIAL
      anisotropy%axis_count = 1
      anisotropy%axis = 0.0_dblprec
      anisotropy%axis(3,1,:) = 1.0_dblprec
      anisotropy%k1 = 0.0_dblprec
      anisotropy%k2 = 0.0_dblprec
      anisotropy%k1(1,:) = k1
      anisotropy%k2(1,:) = k2
   end subroutine make_anisotropy

   subroutine finite_energy_derivative(operator,direction,external,selected,tangent,term,numerical)
      type(coarse_tensor_operator_type), intent(in) :: operator
      real(dblprec), intent(in) :: direction(:,:), external(:,:), tangent(3)
      integer, intent(in) :: selected, term
      real(dblprec), intent(out) :: numerical
      type(coarse_energy_terms_type) :: plus_energy, minus_energy
      integer :: status
      character(len=512) :: message
      real(dblprec), parameter :: epsilon = 2.0d-7
      real(dblprec) :: plus(3,size(direction,2)), minus(3,size(direction,2))
      real(dblprec) :: field(3,size(direction,2)), dipole(3,size(direction,2))

      plus = direction
      minus = direction
      plus(:,selected) = plus(:,selected)+epsilon*tangent
      minus(:,selected) = minus(:,selected)-epsilon*tangent
      plus(:,selected) = plus(:,selected)/vector_norm(plus(:,selected))
      minus(:,selected) = minus(:,selected)/vector_norm(minus(:,selected))
      call local_symmetric_dipole(plus,dipole)
      call evaluate_coarse_tensor_operator(operator,plus,field,plus_energy,status,message,external,dipole)
      call local_symmetric_dipole(minus,dipole)
      call evaluate_coarse_tensor_operator(operator,minus,field,minus_energy,status,message,external,dipole)
      select case (term)
      case (1)
         numerical = (plus_energy%exchange_j-minus_energy%exchange_j)/(2.0_dblprec*epsilon)
      case (2)
         numerical = (plus_energy%spiralization_j-minus_energy%spiralization_j)/(2.0_dblprec*epsilon)
      case (3)
         numerical = (plus_energy%anisotropy_j-minus_energy%anisotropy_j)/(2.0_dblprec*epsilon)
      case (4)
         numerical = (plus_energy%external_j-minus_energy%external_j)/(2.0_dblprec*epsilon)
      case default
         numerical = (plus_energy%dipole_j-minus_energy%dipole_j)/(2.0_dblprec*epsilon)
      end select
   end subroutine finite_energy_derivative

   subroutine local_symmetric_dipole(direction,field)
      real(dblprec), intent(in) :: direction(:,:)
      real(dblprec), intent(out) :: field(:,:)
      real(dblprec), parameter :: tensor(3,3) = reshape((/ &
         -0.12_dblprec,0.03_dblprec,0.01_dblprec, &
          0.03_dblprec,0.08_dblprec,-0.02_dblprec, &
          0.01_dblprec,-0.02_dblprec,0.04_dblprec/),(/3,3/))
      integer :: block

      do block = 1, size(direction,2)
         field(:,block) = matmul(tensor,direction(:,block))
      end do
   end subroutine local_symmetric_dipole

   subroutine allocate_if_needed(first,n,second,third,fourth)
      real(dblprec), allocatable, intent(inout) :: first(:,:), second(:,:), third(:,:), fourth(:,:)
      integer, intent(in) :: n

      if (allocated(first)) deallocate(first)
      if (allocated(second)) deallocate(second)
      if (allocated(third)) deallocate(third)
      if (allocated(fourth)) deallocate(fourth)
      allocate(first(3,n),second(3,n),third(3,n),fourth(3,n))
   end subroutine allocate_if_needed

   pure real(dblprec) function wall_angle(x,center,width)
      real(dblprec), intent(in) :: x, center, width
      real(dblprec) :: argument

      argument = max(-50.0_dblprec,min(50.0_dblprec,(x-center)/width))
      wall_angle = 2.0_dblprec*atan(exp(argument))
   end function wall_angle

   pure function cross_product(left,right) result(cross)
      real(dblprec), intent(in) :: left(3), right(3)
      real(dblprec) :: cross(3)

      cross = (/left(2)*right(3)-left(3)*right(2), &
         left(3)*right(1)-left(1)*right(3), &
         left(1)*right(2)-left(2)*right(1)/)
   end function cross_product

   pure real(dblprec) function det3(matrix)
      real(dblprec), intent(in) :: matrix(3,3)

      det3 = matrix(1,1)*(matrix(2,2)*matrix(3,3)-matrix(2,3)*matrix(3,2)) - &
         matrix(1,2)*(matrix(2,1)*matrix(3,3)-matrix(2,3)*matrix(3,1)) + &
         matrix(1,3)*(matrix(2,1)*matrix(3,2)-matrix(2,2)*matrix(3,1))
   end function det3

   pure function inverse_matrix(matrix) result(inverse)
      real(dblprec), intent(in) :: matrix(3,3)
      real(dblprec) :: inverse(3,3), determinant

      determinant = det3(matrix)
      inverse(1,1) = (matrix(2,2)*matrix(3,3)-matrix(2,3)*matrix(3,2))/determinant
      inverse(1,2) = (matrix(1,3)*matrix(3,2)-matrix(1,2)*matrix(3,3))/determinant
      inverse(1,3) = (matrix(1,2)*matrix(2,3)-matrix(1,3)*matrix(2,2))/determinant
      inverse(2,1) = (matrix(2,3)*matrix(3,1)-matrix(2,1)*matrix(3,3))/determinant
      inverse(2,2) = (matrix(1,1)*matrix(3,3)-matrix(1,3)*matrix(3,1))/determinant
      inverse(2,3) = (matrix(1,3)*matrix(2,1)-matrix(1,1)*matrix(2,3))/determinant
      inverse(3,1) = (matrix(2,1)*matrix(3,2)-matrix(2,2)*matrix(3,1))/determinant
      inverse(3,2) = (matrix(1,2)*matrix(3,1)-matrix(1,1)*matrix(3,2))/determinant
      inverse(3,3) = (matrix(1,1)*matrix(2,2)-matrix(1,2)*matrix(2,1))/determinant
   end function inverse_matrix

   pure real(dblprec) function vector_norm(vector)
      real(dblprec), intent(in) :: vector(3)
      vector_norm = sqrt(sum(vector*vector))
   end function vector_norm

   character(len=1) function digit(value)
      integer, intent(in) :: value
      write(digit,'(i1)') value
   end function digit

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

end program test_coarse_tensor_operator
