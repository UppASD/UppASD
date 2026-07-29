program test_multichannel_coarse_tensor_operator

   use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
   use Parameters, only : dblprec
   use BlockTopology
   use stiffness
   use CoarseTensorOperator, only : COARSE_TENSOR_OK, coarse_operator_options_type
   use MultiChannelCoarseTensorOperator

   implicit none

   integer :: failures

   failures = 0
   call test_two_sublattice_material_and_modes()
   call test_channel_resolved_operator_and_texture()
   if (failures /= 0) then
      write(*,'(a,i0)') 'multichannel coarse tensor operator tests failed: ',failures
      stop 1
   end if
   write(*,'(a)') 'multichannel coarse tensor operator tests passed'

contains

   subroutine test_two_sublattice_material_and_modes()
      type(coarse_lattice_input_type) :: input
      type(coarse_material_type) :: material
      integer :: status
      character(len=512) :: message
      real(dblprec), parameter :: a=2.0d-10, pair=3.0d-21
      real(dblprec) :: omega

      call make_two_channel_input(input,a,pair)
      call extract_coarse_material(input,material,status,message)
      call check(status == COARSE_MATERIAL_OK .and. material%n_channels == 2, &
         'two-channel extraction preserves explicit channels')
      call check_close(material%local_exchange(1,2),-2.0_dblprec*pair/(a**3),1.0d-13, &
         'intersublattice exchange retains the atomistic sign')
      call check(all(material%channel_moment_mub == (/1.0_dblprec,1.5_dblprec/)), &
         'separate channel moments are not collapsed')
      omega=abs(material%local_exchange(1,2))*material%cell_volume_m3/COARSE_MUB_SI * &
         (material%channel_gamma(1)/material%channel_moment_mub(1)+ &
          material%channel_gamma(2)/material%channel_moment_mub(2))
      call validate_coarse_material_two_sublattice_modes(material,0.0_dblprec,omega,1.0d-12,status,message)
      call check(status == COARSE_MATERIAL_OK .and. material%diagnostics%acoustic_mode_extraction_validated .and. &
         material%diagnostics%optical_mode_extraction_validated, &
         'controlled atomistic acoustic and optical modes accept the material')
      call coarse_material_runtime_status(material,COARSE_RUNTIME_FERRI_AFM,status,message)
      call check(status == COARSE_MATERIAL_OK, 'runtime accepts only after the two-mode diagnostic passes')
   end subroutine test_two_sublattice_material_and_modes

   subroutine test_channel_resolved_operator_and_texture()
      type(coarse_lattice_input_type) :: input
      type(coarse_material_type) :: material
      type(block_topology_type) :: topology
      type(multichannel_coarse_tensor_operator_type) :: operator
      type(multichannel_coarse_energy_terms_type) :: energy
      type(multichannel_coarse_field_terms_type) :: terms
      type(coarse_operator_options_type) :: options
      integer :: status, block
      integer :: channel_map(2), repetitions(3)
      character(len=512) :: message
      real(dblprec), parameter :: a=2.0d-10, pair=3.0d-21
      real(dblprec) :: cell(3,3), omega, expected_field(3)
      real(dblprec) :: direction(3,2,8), field(3,2,8), rhs(3,2,8), angle

      call make_two_channel_input(input,a,pair)
      call extract_coarse_material(input,material,status,message)
      omega=abs(material%local_exchange(1,2))*material%cell_volume_m3/COARSE_MUB_SI * &
         (material%channel_gamma(1)/material%channel_moment_mub(1)+ &
          material%channel_gamma(2)/material%channel_moment_mub(2))
      call validate_coarse_material_two_sublattice_modes(material,0.0_dblprec,omega,1.0d-12,status,message)
      cell=0.0_dblprec; cell(1,1)=a; cell(2,2)=a; cell(3,3)=a
      channel_map=(/1,2/); repetitions=(/8,1,1/)
      call build_block_topology(topology,REGULAR_REPLICATED_CELL,2,repetitions,16,(/1,1,1/),cell,channel_map,status,message)
      call check(status == BLOCK_TOPOLOGY_OK, 'two-channel topology builds')
      call setup_multichannel_coarse_tensor_operator(operator,topology,material,options,status,message)
      call check(status == COARSE_TENSOR_OK, 'accepted two-channel material enables CPU reference only')
      if (status /= COARSE_TENSOR_OK) return

      direction=0.0_dblprec; direction(3,1,:)=1.0_dblprec; direction(3,2,:)=-1.0_dblprec
      call evaluate_multichannel_coarse_tensor_operator(operator,direction,field,energy,status,message,terms)
      call check(status == COARSE_TENSOR_OK .and. all(ieee_is_finite(field)), &
         'compensated AFM blocks retain finite independent channel fields')
      call multichannel_coarse_llg_rhs(operator,direction,field,rhs,status,message)
      call check(status == COARSE_TENSOR_OK .and. maxval(abs(rhs)) < 1.0d-5, &
         'collinear AFM ground state has zero sublattice torque without net-moment normalization')

      direction(:,1,:)=spread((/0.0_dblprec,0.0_dblprec,1.0_dblprec/),2,8)
      direction(:,2,:)=spread((/1.0_dblprec,0.0_dblprec,0.0_dblprec/),2,8)
      call evaluate_multichannel_coarse_tensor_operator(operator,direction,field,energy,status,message,terms)
      expected_field=-material%local_exchange(1,2)*operator%block_volume_m3 / &
         (COARSE_MUB_SI*operator%block_moment_mub(1,1))*direction(:,2,1)
      call check_close(maxval(abs(terms%local_exchange_t(:,1,1)-expected_field)),0.0_dblprec,1.0d-11, &
         'local intersublattice field has the negative energy derivative sign')

      do block=1,operator%nblocks
         angle=2.0_dblprec*atan(exp((real(block,dblprec)-4.5_dblprec)/1.5_dblprec))
         direction(:,1,block)=(/sin(angle),0.0_dblprec,cos(angle)/)
         direction(:,2,block)=-direction(:,1,block)
      end do
      call evaluate_multichannel_coarse_tensor_operator(operator,direction,field,energy,status,message,terms)
      call check(status == COARSE_TENSOR_OK .and. abs(energy%exchange_j) > 0.0_dblprec .and. &
         all(ieee_is_finite(field)), 'two-sublattice AFM domain-wall texture evaluates with channel-resolved stiffness')
   end subroutine test_channel_resolved_operator_and_texture

   subroutine make_two_channel_input(input,a,pair)
      type(coarse_lattice_input_type), intent(out) :: input
      real(dblprec), intent(in) :: a,pair
      input%n_basis=2; input%n_channels=2; input%cell_volume_m3=a**3; input%fit_first=1; input%fit_last=1
      allocate(input%basis_to_channel(2),input%basis_moment_mub(2),input%channel_gamma(2),input%channel_damping(2),input%eta_inverse_m(1))
      allocate(input%exchange_source_basis(4),input%exchange_target_basis(4),input%exchange_displacement_m(3,4),input%exchange_pair_energy_j(4))
      allocate(input%dmi_source_basis(0),input%dmi_target_basis(0),input%dmi_displacement_m(3,0),input%dmi_pair_energy_j(3,0))
      input%basis_to_channel=(/1,2/); input%basis_moment_mub=(/1.0_dblprec,1.5_dblprec/)
      input%channel_gamma=(/2.0d11,3.0d11/); input%channel_damping=(/0.01_dblprec,0.03_dblprec/); input%eta_inverse_m=0.0_dblprec
      input%exchange_source_basis=(/1,1,2,2/); input%exchange_target_basis=(/2,2,1,1/); input%exchange_pair_energy_j=pair
      input%exchange_displacement_m=0.0_dblprec; input%exchange_displacement_m(1,:)=(/0.5_dblprec*a,-0.5_dblprec*a,0.5_dblprec*a,-0.5_dblprec*a/)
   end subroutine make_two_channel_input

   subroutine check(condition,message)
      logical,intent(in)::condition
      character(len=*),intent(in)::message
      if (.not. condition) then
         failures=failures+1; write(*,'(a)') 'FAIL: '//trim(message)
      end if
   end subroutine check

   subroutine check_close(actual,expected,tolerance,message)
      real(dblprec),intent(in)::actual,expected,tolerance
      character(len=*),intent(in)::message
      call check(abs(actual-expected) <= tolerance*max(1.0_dblprec,abs(actual),abs(expected)),message)
   end subroutine check_close

end program test_multichannel_coarse_tensor_operator
