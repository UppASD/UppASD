program test_coarse_dispersion

   ! RCG-04-FU4: measure the production coarse tensor operator on periodic
   ! Fourier modes.  The benchmark deliberately calls the same production
   ! setup/evaluation entry points used by the all-coarse runtime; the
   ! formulas below are only the independent symbol oracle.
   use Parameters, only : dblprec
   use BlockTopology
   use stiffness, only : coarse_material_type, COARSE_MUB_SI
   use CoarseTensorOperator

   implicit none

   integer :: failures
   real(dblprec), parameter :: pi = acos(-1.0_dblprec)

   failures = 0
   call test_periodic_symbols()

   if (failures /= 0) then
      write(*,'(a,i0)') 'coarse dispersion tests failed: ', failures
      stop 1
   end if
   write(*,'(a)') 'coarse dispersion tests passed'

contains

   subroutine test_periodic_symbols()
      integer, parameter :: n = 64
      integer, parameter :: n_scales = 4, n_modes = 6
      integer, parameter :: modes(n_modes) = (/1,2,4,8,16,24/)
      real(dblprec), parameter :: h_values(n_scales) = &
         (/1.0d-10,2.0d-10,4.0d-10,8.0d-10/)
      real(dblprec), parameter :: exchange_a = 1.0d-11
      real(dblprec), parameter :: dmi_d = 2.0d-3
      real(dblprec), parameter :: epsilon = 1.0d-3
      real(dblprec) :: cell(3,3), exchange(3,3), dmi(3,3)
      real(dblprec) :: qh, theta, q, h, expected_exchange_symbol
      real(dblprec) :: expected_dmi_symbol, expected_exchange_energy
      real(dblprec) :: expected_dmi_energy, expected_exchange_curvature
      real(dblprec) :: expected_dmi_curvature, exchange_ratio, dmi_ratio
      real(dblprec) :: exchange_continuum_ratio, dmi_continuum_ratio
      real(dblprec) :: exchange_curvature, dmi_curvature
      real(dblprec) :: exchange_energy, dmi_energy, transverse_dot
      real(dblprec) :: max_exchange_residual, max_dmi_residual
      real(dblprec) :: sum_exchange_ratio, sum_dmi_ratio
      real(dblprec) :: direction(3,n), field(3,n)
      type(block_topology_type) :: topology
      type(coarse_material_type) :: exchange_material, dmi_material
      type(coarse_tensor_operator_type) :: exchange_operator, dmi_operator
      type(coarse_operator_options_type) :: options
      type(coarse_energy_terms_type) :: energy
      integer :: scale_index, mode_index, block, status
      character(len=512) :: message

      max_exchange_residual = 0.0_dblprec
      max_dmi_residual = 0.0_dblprec
      sum_exchange_ratio = 0.0_dblprec
      sum_dmi_ratio = 0.0_dblprec

      write(*,'(a)') 'RCG-04-FU4 periodic production-operator dispersion sweep'
      write(*,'(a)') 'columns: h_nm mode qh exchange_continuum_ratio exchange_energy_over_discrete exchange_field_over_discrete dmi_continuum_ratio dmi_energy_over_discrete dmi_field_over_discrete'

      do scale_index = 1, n_scales
         h = h_values(scale_index)
         cell = 0.0_dblprec
         cell(1,1) = h
         cell(2,2) = h
         cell(3,3) = h

         exchange = 0.0_dblprec
         exchange = 0.0_dblprec
         exchange(1,1) = exchange_a
         dmi = 0.0_dblprec
         call make_material_and_topology(n,cell,exchange,dmi,topology, &
            exchange_material,status,message)
         call check(status == BLOCK_TOPOLOGY_OK,'exchange topology builds: '//trim(message))
         if (status /= BLOCK_TOPOLOGY_OK) cycle
         call setup_coarse_tensor_operator(exchange_operator,topology, &
            exchange_material,options,status,message)
         call check(status == COARSE_TENSOR_OK,'exchange operator sets up: '//trim(message))
         if (status /= COARSE_TENSOR_OK) cycle

         dmi = 0.0_dblprec
         dmi(3,1) = dmi_d
         call make_material_and_topology(n,cell,0.0_dblprec*exchange,dmi, &
            topology,dmi_material,status,message)
         call check(status == BLOCK_TOPOLOGY_OK,'DMI topology builds: '//trim(message))
         if (status /= BLOCK_TOPOLOGY_OK) cycle
         call setup_coarse_tensor_operator(dmi_operator,topology,dmi_material, &
            options,status,message)
         call check(status == COARSE_TENSOR_OK,'DMI operator sets up: '//trim(message))
         if (status /= COARSE_TENSOR_OK) cycle

         do mode_index = 1, n_modes
            qh = 2.0_dblprec*pi*real(modes(mode_index),dblprec)/real(n,dblprec)
            theta = qh
            q = qh/h
            direction(3,:) = sqrt(1.0_dblprec-epsilon**2)
            do block = 1, n
               direction(1,block) = epsilon*cos(theta*real(block-1,dblprec))
               direction(2,block) = epsilon*sin(theta*real(block-1,dblprec))
            end do

            call evaluate_coarse_tensor_operator(exchange_operator,direction,field, &
               energy,status,message)
            call check(status == COARSE_TENSOR_OK,'exchange Fourier mode evaluates: '//trim(message))
            exchange_energy = energy%exchange_j
            transverse_dot = sum(field(1,:)*direction(1,:) + &
               field(2,:)*direction(2,:))
            exchange_curvature = -transverse_dot/(real(n,dblprec)*epsilon**2)

            call evaluate_coarse_tensor_operator(dmi_operator,direction,field, &
               energy,status,message)
            call check(status == COARSE_TENSOR_OK,'DMI Fourier mode evaluates: '//trim(message))
            dmi_energy = energy%spiralization_j
            transverse_dot = sum(field(1,:)*direction(1,:) + &
               field(2,:)*direction(2,:))
            dmi_curvature = -transverse_dot/(real(n,dblprec)*epsilon**2)

            expected_exchange_symbol = 4.0_dblprec*sin(0.5_dblprec*qh)**2/h**2
            expected_dmi_symbol = sin(qh)/h
            exchange_continuum_ratio = expected_exchange_symbol/q**2
            dmi_continuum_ratio = expected_dmi_symbol/q
            expected_exchange_energy = real(n,dblprec)*h**3*exchange_a* &
               epsilon**2*expected_exchange_symbol
            expected_dmi_energy = real(n,dblprec)*h**3*dmi_d*epsilon**2* &
               expected_dmi_symbol
            expected_exchange_curvature = 2.0_dblprec*h**3*exchange_a* &
               expected_exchange_symbol/(COARSE_MUB_SI*exchange_operator%block_moment_mub)
            expected_dmi_curvature = 2.0_dblprec*h**3*dmi_d*expected_dmi_symbol/ &
               (COARSE_MUB_SI*dmi_operator%block_moment_mub)

            exchange_ratio = exchange_energy/expected_exchange_energy
            dmi_ratio = dmi_energy/expected_dmi_energy
            write(*,'(f8.3,1x,i4,1x,f8.4,1x,6(es14.6,1x))') h*1.0d9, &
               modes(mode_index), qh, exchange_continuum_ratio, exchange_ratio, &
               exchange_curvature/expected_exchange_curvature, dmi_continuum_ratio, &
               dmi_ratio, dmi_curvature/expected_dmi_curvature
            call check_close(exchange_energy,expected_exchange_energy,2.0d-12, &
               'exchange energy follows production-stencil symbol')
            call check_close(exchange_curvature,expected_exchange_curvature,2.0d-12, &
               'exchange field curvature follows production-stencil symbol')
            call check_close(dmi_energy,expected_dmi_energy,2.0d-12, &
               'DMI energy follows production-stencil symbol')
            call check_close(dmi_curvature,expected_dmi_curvature,2.0d-12, &
               'DMI field curvature follows production-stencil symbol')
            max_exchange_residual = max(max_exchange_residual, &
               abs(exchange_ratio-1.0_dblprec), &
               abs(exchange_curvature/expected_exchange_curvature-1.0_dblprec))
            max_dmi_residual = max(max_dmi_residual,abs(dmi_ratio-1.0_dblprec), &
               abs(dmi_curvature/expected_dmi_curvature-1.0_dblprec))
            sum_exchange_ratio = sum_exchange_ratio + exchange_ratio
            sum_dmi_ratio = sum_dmi_ratio + dmi_ratio
         end do
      end do

      write(*,'(a,es14.6)') 'max exchange relative residual: ',max_exchange_residual
      write(*,'(a,es14.6)') 'max DMI relative residual: ',max_dmi_residual
      write(*,'(a,es14.6)') 'exchange energy multiplicative fit: ', &
         sum_exchange_ratio/real(n_scales*n_modes,dblprec)
      write(*,'(a,es14.6)') 'DMI energy multiplicative fit: ', &
         sum_dmi_ratio/real(n_scales*n_modes,dblprec)
   end subroutine test_periodic_symbols

   subroutine make_material_and_topology(n,cell,exchange,dmi,topology,material,status,message)
      integer, intent(in) :: n
      real(dblprec), intent(in) :: cell(3,3), exchange(3,3), dmi(3,3)
      type(block_topology_type), intent(out) :: topology
      type(coarse_material_type), intent(out) :: material
      integer, intent(out) :: status
      character(len=*), intent(out) :: message

      call build_block_topology(topology,REGULAR_REPLICATED_CELL,1,(/n,1,1/),n, &
         (/1,1,1/),cell,(/1/),status,message)
      if (status /= BLOCK_TOPOLOGY_OK) return
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
      material%channel_gamma = 1.76d11
      material%channel_damping = 0.0_dblprec
      material%exchange_stiffness(:,:,1,1) = exchange
      material%spiralization(:,:,1,1) = dmi
      material%diagnostics%channel_dynamics_parameters_validated = .true.
      material%diagnostics%small_q_energy_validated = .true.
   end subroutine make_material_and_topology

   pure real(dblprec) function det3(matrix)
      real(dblprec), intent(in) :: matrix(3,3)
      det3 = matrix(1,1)*(matrix(2,2)*matrix(3,3)-matrix(2,3)*matrix(3,2)) - &
         matrix(1,2)*(matrix(2,1)*matrix(3,3)-matrix(2,3)*matrix(3,1)) + &
         matrix(1,3)*(matrix(2,1)*matrix(3,2)-matrix(2,2)*matrix(3,1))
   end function det3

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
      scale = max(abs(actual),abs(expected),tiny(1.0_dblprec))
      call check(abs(actual-expected) <= relative_tolerance*scale,message)
   end subroutine check_close

end program test_coarse_dispersion
