!-------------------------------------------------------------------------------
! MODULE: MultiChannelCoarseTensorOperator
!> @brief Controlled two-sublattice CPU reference for CG-11.
!>
!> This is deliberately separate from the established single-channel operator.
!> It consumes only a material whose acoustic/optical extraction diagnostics
!> have passed, keeps every channel direction and moment independently, and
!> does not deposit a net compensated macrospin.
!-------------------------------------------------------------------------------
module MultiChannelCoarseTensorOperator

   use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
   use Parameters, only : dblprec
   use BlockTopology, only : block_topology_type, REGULAR_REPLICATED_CELL, &
      regular_spatial_block_id
   use stiffness, only : coarse_material_type, coarse_material_runtime_status, &
      COARSE_RUNTIME_FERRI_AFM, COARSE_MATERIAL_OK, COARSE_MUB_SI
   use CoarseTensorOperator, only : COARSE_TENSOR_OK, COARSE_TENSOR_INVALID_TOPOLOGY, &
      COARSE_TENSOR_INVALID_MATERIAL, COARSE_TENSOR_UNSUPPORTED_GEOMETRY, &
      COARSE_TENSOR_UNSUPPORTED_SOLVER, COARSE_TENSOR_UNSUPPORTED_TEMPERATURE, &
      COARSE_TENSOR_UNSUPPORTED_MODE, COARSE_TENSOR_INVALID_STATE, &
      COARSE_BOUNDARY_PERIODIC, COARSE_SOLVER_DETERMINISTIC_HEUN, &
      coarse_operator_options_type

   implicit none
   private

   type, public :: multichannel_coarse_tensor_operator_type
      logical :: ready = .false.
      integer :: nblocks = 0
      integer :: nchannels = 0
      integer :: block_grid(3) = 0
      integer, allocatable :: block_coordinate(:,:)
      real(dblprec) :: block_vectors_m(3,3) = 0.0_dblprec
      real(dblprec) :: inverse_block_transpose_m1(3,3) = 0.0_dblprec
      real(dblprec) :: block_volume_m3 = 0.0_dblprec
      real(dblprec), allocatable :: block_moment_mub(:,:) ! (channel,block)
      real(dblprec), allocatable :: local_exchange_j_per_m3(:,:)
      real(dblprec), allocatable :: exchange_stiffness_j_per_m(:,:,:,:)
      real(dblprec), allocatable :: spiralization_j_per_m2(:,:,:,:)
      real(dblprec), allocatable :: channel_gamma_per_t_s(:)
      real(dblprec), allocatable :: channel_damping(:)
   end type multichannel_coarse_tensor_operator_type

   type, public :: multichannel_coarse_energy_terms_type
      real(dblprec) :: local_exchange_j = 0.0_dblprec
      real(dblprec) :: exchange_j = 0.0_dblprec
      real(dblprec) :: spiralization_j = 0.0_dblprec
      real(dblprec) :: total_j = 0.0_dblprec
   end type multichannel_coarse_energy_terms_type

   type, public :: multichannel_coarse_field_terms_type
      real(dblprec), allocatable :: local_exchange_t(:,:,:)
      real(dblprec), allocatable :: exchange_t(:,:,:)
      real(dblprec), allocatable :: spiralization_t(:,:,:)
   end type multichannel_coarse_field_terms_type

   public :: setup_multichannel_coarse_tensor_operator
   public :: evaluate_multichannel_coarse_tensor_operator
   public :: multichannel_coarse_llg_rhs

contains

   subroutine setup_multichannel_coarse_tensor_operator(operator,topology,material,options, &
         status,diagnostic)
      type(multichannel_coarse_tensor_operator_type), intent(out) :: operator
      type(block_topology_type), intent(in) :: topology
      type(coarse_material_type), intent(in) :: material
      type(coarse_operator_options_type), intent(in) :: options
      integer, intent(out) :: status
      character(len=*), intent(out) :: diagnostic

      integer :: runtime_status, channel, block, cells_per_block, c1, c2
      real(dblprec) :: determinant, cell_volume, scale
      character(len=512) :: runtime_diagnostic

      status = COARSE_TENSOR_INVALID_TOPOLOGY
      diagnostic = ''
      if (.not. topology%ready .or. topology%geometry_mode /= REGULAR_REPLICATED_CELL) then
         diagnostic = 'Two-sublattice tensor setup requires a ready regular topology'
         return
      end if
      if (topology%n_dynamic_channels /= 2 .or. material%n_channels /= 2) then
         status = COARSE_TENSOR_UNSUPPORTED_MODE
         diagnostic = 'CG-11 CPU reference supports exactly two explicit dynamical channels'
         return
      end if
      if (any(topology%block_dynamic_channel_population <= 0)) then
         diagnostic = 'Every controlled two-sublattice block must contain both channels'
         return
      end if
      if (any(options%boundary_mode /= COARSE_BOUNDARY_PERIODIC)) then
         status = COARSE_TENSOR_UNSUPPORTED_GEOMETRY
         diagnostic = 'Two-sublattice tensor setup supports periodic derivatives only'
         return
      end if
      if (options%solver_mode /= COARSE_SOLVER_DETERMINISTIC_HEUN) then
         status = COARSE_TENSOR_UNSUPPORTED_SOLVER
         diagnostic = 'Two-sublattice tensor setup supports deterministic Heun LLG only'
         return
      end if
      if (.not. ieee_is_finite(options%temperature_k) .or. options%temperature_k /= 0.0_dblprec .or. &
          options%stochastic_field) then
         status = COARSE_TENSOR_UNSUPPORTED_TEMPERATURE
         diagnostic = 'Two-sublattice tensor setup requires T=0 and no stochastic field'
         return
      end if
      if (options%adaptive_switching .or. options%interface_coupling .or. options%restart_state .or. &
          options%time_dependent_external_field .or. options%use_uniform_coarse_dipole) then
         status = COARSE_TENSOR_UNSUPPORTED_MODE
         diagnostic = 'CG-11 is a static CPU two-sublattice reference; adaptive, interface, restart, dipole, and time-dependent modes are rejected'
         return
      end if
      call coarse_material_runtime_status(material,COARSE_RUNTIME_FERRI_AFM,runtime_status,runtime_diagnostic)
      if (runtime_status /= COARSE_MATERIAL_OK) then
         status = COARSE_TENSOR_INVALID_MATERIAL
         diagnostic = 'Two-sublattice material rejected: '//trim(runtime_diagnostic)
         return
      end if
      if (.not. allocated(material%basis_to_channel) .or. .not. allocated(material%basis_moment_mub) .or. &
          .not. allocated(material%channel_moment_mub) .or. .not. allocated(material%channel_gamma) .or. &
          .not. allocated(material%channel_damping) .or. .not. allocated(material%local_exchange) .or. &
          .not. allocated(material%exchange_stiffness) .or. .not. allocated(material%spiralization) .or. &
          material%n_basis /= topology%n_basis .or. size(material%basis_to_channel) /= topology%n_basis .or. &
          any(material%basis_to_channel /= topology%basis_to_dynamic_channel)) then
         status = COARSE_TENSOR_INVALID_MATERIAL
         diagnostic = 'Two-sublattice material arrays and basis/channel map must match topology'
         return
      end if
      if (any(shape(material%local_exchange) /= (/2,2/)) .or. &
          any(shape(material%exchange_stiffness) /= (/3,3,2,2/)) .or. &
          any(shape(material%spiralization) /= (/3,3,2,2/)) .or. &
          .not. all(ieee_is_finite(material%local_exchange)) .or. &
          .not. all(ieee_is_finite(material%exchange_stiffness)) .or. &
          .not. all(ieee_is_finite(material%spiralization)) .or. &
          .not. all(ieee_is_finite(material%channel_moment_mub)) .or. &
          .not. all(ieee_is_finite(material%channel_gamma)) .or. &
          .not. all(ieee_is_finite(material%channel_damping)) .or. &
          any(material%channel_moment_mub <= 0.0_dblprec) .or. any(material%channel_gamma <= 0.0_dblprec) .or. &
          any(material%channel_damping < 0.0_dblprec)) then
         status = COARSE_TENSOR_INVALID_MATERIAL
         diagnostic = 'Two-sublattice material must contain finite tensors and valid per-channel dynamics'
         return
      end if
      if (maxval(abs(material%local_exchange-transpose(material%local_exchange))) > &
          1.0d-12*max(1.0_dblprec,maxval(abs(material%local_exchange)))) then
         status = COARSE_TENSOR_INVALID_MATERIAL
         diagnostic = 'Two-sublattice local exchange must be channel-symmetric'
         return
      end if
      ! RCG-03: mixed-derivative accuracy of the gradient discretization
      ! relies on Cartesian symmetry of the per-channel-pair spatial
      ! exchange stiffness tensor, mirroring the single-channel assertion in
      ! CoarseTensorOperator; local_exchange's channel symmetry above does
      ! not cover this separate spatial tensor.
      do c2 = 1, 2
         do c1 = 1, 2
            if (maxval(abs(material%exchange_stiffness(:,:,c1,c2) - &
                       transpose(material%exchange_stiffness(:,:,c1,c2)))) > &
                1.0d-12*max(tiny(1.0_dblprec), &
                   maxval(abs(material%exchange_stiffness(:,:,c1,c2))))) then
               status = COARSE_TENSOR_INVALID_MATERIAL
               diagnostic = 'Two-sublattice exchange stiffness must be Cartesian-symmetric per channel pair'
               return
            end if
         end do
      end do
      determinant = determinant3(topology%block_vectors)
      cell_volume = abs(determinant3(topology%cell_vectors))
      scale = max(cell_volume,material%cell_volume_m3)
      if (abs(determinant) <= tiny(determinant) .or. material%cell_volume_m3 <= 0.0_dblprec .or. &
          abs(cell_volume-material%cell_volume_m3) > 1.0d-10*scale) then
         status = COARSE_TENSOR_INVALID_MATERIAL
         diagnostic = 'Two-sublattice material volume must match the physical topology cell metric'
         return
      end if

      operator%nblocks = topology%n_spatial_blocks
      operator%nchannels = 2
      operator%block_grid = topology%block_grid
      operator%block_vectors_m = topology%block_vectors
      operator%block_volume_m3 = abs(determinant)
      operator%inverse_block_transpose_m1 = transpose(inverse3(topology%block_vectors))
      allocate(operator%block_coordinate(3,operator%nblocks),operator%block_moment_mub(2,operator%nblocks), &
         operator%local_exchange_j_per_m3(2,2),operator%exchange_stiffness_j_per_m(3,3,2,2), &
         operator%spiralization_j_per_m2(3,3,2,2),operator%channel_gamma_per_t_s(2),operator%channel_damping(2))
      operator%block_coordinate = topology%block_grid_coordinate
      cells_per_block = product(topology%block_shape)
      do channel = 1, 2
         if (any(topology%block_dynamic_channel_population(channel,:) /= &
             count(material%basis_to_channel == channel)*cells_per_block)) then
            status = COARSE_TENSOR_INVALID_MATERIAL
            diagnostic = 'Material channel population does not match every canonical block'
            return
         end if
         do block = 1, operator%nblocks
            operator%block_moment_mub(channel,block) = material%channel_moment_mub(channel)*cells_per_block
         end do
      end do
      operator%local_exchange_j_per_m3 = material%local_exchange
      operator%exchange_stiffness_j_per_m = material%exchange_stiffness
      operator%spiralization_j_per_m2 = material%spiralization
      operator%channel_gamma_per_t_s = material%channel_gamma
      operator%channel_damping = material%channel_damping
      operator%ready = .true.
      status = COARSE_TENSOR_OK
   end subroutine setup_multichannel_coarse_tensor_operator

   subroutine evaluate_multichannel_coarse_tensor_operator(operator,direction,field_t,energies, &
         status,diagnostic,term_fields)
      type(multichannel_coarse_tensor_operator_type), intent(in) :: operator
      real(dblprec), intent(in) :: direction(:,:,:)
      real(dblprec), intent(out) :: field_t(:,:,:)
      type(multichannel_coarse_energy_terms_type), intent(out) :: energies
      integer, intent(out) :: status
      character(len=*), intent(out) :: diagnostic
      type(multichannel_coarse_field_terms_type), intent(out), optional :: term_fields

      integer :: block, channel, other, p, q, k
      real(dblprec) :: e(3), conversion
      real(dblprec), allocatable :: gradient(:,:,:,:), local_derivative(:,:,:), exchange_derivative(:,:,:), dmi_derivative(:,:,:)
      real(dblprec), allocatable :: work(:,:), local_field(:,:,:), exchange_field(:,:,:), dmi_field(:,:,:)

      energies = multichannel_coarse_energy_terms_type()
      status = COARSE_TENSOR_INVALID_STATE
      diagnostic = ''
      if (.not. operator%ready .or. any(shape(direction) /= (/3,operator%nchannels,operator%nblocks/)) .or. &
          any(shape(field_t) /= shape(direction))) then
         diagnostic = 'Two-sublattice directions and fields must have shape (3,channel,block)'
         return
      end if
      if (.not. all(ieee_is_finite(direction))) then
         diagnostic = 'Two-sublattice directions must be finite'
         return
      end if
      do block=1,operator%nblocks
         do channel=1,operator%nchannels
            if (abs(sum(direction(:,channel,block)**2)-1.0_dblprec) > 1.0d-8) then
               diagnostic = 'Every two-sublattice channel direction must have unit norm'
               return
            end if
         end do
      end do
      allocate(gradient(3,3,operator%nchannels,operator%nblocks),local_derivative(3,operator%nchannels,operator%nblocks), &
         exchange_derivative(3,operator%nchannels,operator%nblocks),dmi_derivative(3,operator%nchannels,operator%nblocks), &
         local_field(3,operator%nchannels,operator%nblocks),exchange_field(3,operator%nchannels,operator%nblocks), &
         dmi_field(3,operator%nchannels,operator%nblocks),work(3,operator%nblocks))
      local_derivative=0.0_dblprec; exchange_derivative=0.0_dblprec; dmi_derivative=0.0_dblprec
      call multichannel_gradient(operator,direction,gradient)

      do block=1,operator%nblocks
         do channel=1,operator%nchannels-1
            do other=channel+1,operator%nchannels
               energies%local_exchange_j = energies%local_exchange_j + operator%block_volume_m3 * &
                  operator%local_exchange_j_per_m3(channel,other)*dot_product(direction(:,channel,block),direction(:,other,block))
               local_derivative(:,channel,block) = local_derivative(:,channel,block) + operator%block_volume_m3 * &
                  operator%local_exchange_j_per_m3(channel,other)*direction(:,other,block)
               local_derivative(:,other,block) = local_derivative(:,other,block) + operator%block_volume_m3 * &
                  operator%local_exchange_j_per_m3(channel,other)*direction(:,channel,block)
            end do
         end do
      end do
      do channel=1,operator%nchannels
         do other=1,operator%nchannels
            do p=1,3
               do q=1,3
                  do block=1,operator%nblocks
                     energies%exchange_j = energies%exchange_j + operator%block_volume_m3* &
                        operator%exchange_stiffness_j_per_m(p,q,channel,other)* &
                        dot_product(gradient(:,p,channel,block),gradient(:,q,other,block))
                     work(:,block) = operator%block_volume_m3*operator%exchange_stiffness_j_per_m(p,q,channel,other)*gradient(:,q,other,block)
                  end do
                  call add_gradient_transpose(operator,p,work,exchange_derivative(:,channel,:),1.0_dblprec)
               end do
            end do
         end do
      end do
      do channel=1,operator%nchannels
         do other=1,operator%nchannels
            do k=1,3
               e=0.0_dblprec; e(k)=1.0_dblprec
               do p=1,3
                  work=0.0_dblprec
                  do block=1,operator%nblocks
                     energies%spiralization_j = energies%spiralization_j + operator%block_volume_m3* &
                        operator%spiralization_j_per_m2(k,p,channel,other)* &
                        dot_product(e,cross3(direction(:,channel,block),gradient(:,p,other,block)))
                     dmi_derivative(:,channel,block) = dmi_derivative(:,channel,block) - operator%block_volume_m3* &
                        operator%spiralization_j_per_m2(k,p,channel,other)*cross3(e,gradient(:,p,other,block))
                     work(:,block) = operator%block_volume_m3*operator%spiralization_j_per_m2(k,p,channel,other)*cross3(e,direction(:,channel,block))
                  end do
                  call add_gradient_transpose(operator,p,work,dmi_derivative(:,other,:),1.0_dblprec)
               end do
            end do
         end do
      end do
      do block=1,operator%nblocks
         do channel=1,operator%nchannels
            conversion=-1.0_dblprec/(COARSE_MUB_SI*operator%block_moment_mub(channel,block))
            local_field(:,channel,block)=conversion*local_derivative(:,channel,block)
            exchange_field(:,channel,block)=conversion*exchange_derivative(:,channel,block)
            dmi_field(:,channel,block)=conversion*dmi_derivative(:,channel,block)
         end do
      end do
      field_t=local_field+exchange_field+dmi_field
      energies%total_j=energies%local_exchange_j+energies%exchange_j+energies%spiralization_j
      if (present(term_fields)) then
         allocate(term_fields%local_exchange_t(3,operator%nchannels,operator%nblocks), &
            term_fields%exchange_t(3,operator%nchannels,operator%nblocks),term_fields%spiralization_t(3,operator%nchannels,operator%nblocks))
         term_fields%local_exchange_t=local_field; term_fields%exchange_t=exchange_field; term_fields%spiralization_t=dmi_field
      end if
      status=COARSE_TENSOR_OK
   end subroutine evaluate_multichannel_coarse_tensor_operator

   subroutine multichannel_coarse_llg_rhs(operator,direction,field_t,derivative_per_s,status,diagnostic)
      type(multichannel_coarse_tensor_operator_type), intent(in) :: operator
      real(dblprec), intent(in) :: direction(:,:,:),field_t(:,:,:)
      real(dblprec), intent(out) :: derivative_per_s(:,:,:)
      integer, intent(out) :: status
      character(len=*), intent(out) :: diagnostic
      integer :: block,channel
      real(dblprec) :: prefactor
      status=COARSE_TENSOR_INVALID_STATE; diagnostic=''
      if (.not. operator%ready .or. any(shape(direction) /= (/3,operator%nchannels,operator%nblocks/)) .or. &
          any(shape(field_t) /= shape(direction)) .or. any(shape(derivative_per_s) /= shape(direction)) .or. &
          .not. all(ieee_is_finite(direction)) .or. .not. all(ieee_is_finite(field_t))) then
         diagnostic='Two-sublattice LLG arrays require a ready operator and finite (3,channel,block) arrays'
         return
      end if
      do block=1,operator%nblocks
         do channel=1,operator%nchannels
            prefactor=-operator%channel_gamma_per_t_s(channel)/(1.0_dblprec+operator%channel_damping(channel)**2)
            derivative_per_s(:,channel,block)=prefactor*(cross3(direction(:,channel,block),field_t(:,channel,block)) + &
               operator%channel_damping(channel)*cross3(direction(:,channel,block),cross3(direction(:,channel,block),field_t(:,channel,block))))
         end do
      end do
      status=COARSE_TENSOR_OK
   end subroutine multichannel_coarse_llg_rhs

   subroutine multichannel_gradient(operator,values,gradient)
      type(multichannel_coarse_tensor_operator_type), intent(in) :: operator
      real(dblprec), intent(in) :: values(:,:,:)
      real(dblprec), intent(out) :: gradient(:,:,:,:)
      integer :: block,channel,r,p,plus_block,coordinate(3)
      gradient=0.0_dblprec
      do block=1,operator%nblocks
         do r=1,3
            coordinate=operator%block_coordinate(:,block); coordinate(r)=modulo(coordinate(r)+1,operator%block_grid(r))
            plus_block=regular_spatial_block_id(coordinate,operator%block_grid)
            do p=1,3
               do channel=1,operator%nchannels
                  gradient(:,p,channel,block)=gradient(:,p,channel,block)+operator%inverse_block_transpose_m1(p,r)*(values(:,channel,plus_block)-values(:,channel,block))
               end do
            end do
         end do
      end do
   end subroutine multichannel_gradient

   subroutine add_gradient_transpose(operator,p,values,result,scale)
      type(multichannel_coarse_tensor_operator_type), intent(in) :: operator
      integer, intent(in) :: p
      real(dblprec), intent(in) :: values(:,:)
      real(dblprec), intent(inout) :: result(:,:)
      real(dblprec), intent(in) :: scale
      integer :: block,r,plus_block,coordinate(3)
      real(dblprec) :: coefficient
      do block=1,operator%nblocks
         do r=1,3
            coefficient=scale*operator%inverse_block_transpose_m1(p,r)
            coordinate=operator%block_coordinate(:,block); coordinate(r)=modulo(coordinate(r)+1,operator%block_grid(r))
            plus_block=regular_spatial_block_id(coordinate,operator%block_grid)
            result(:,plus_block)=result(:,plus_block)+coefficient*values(:,block)
            result(:,block)=result(:,block)-coefficient*values(:,block)
         end do
      end do
   end subroutine add_gradient_transpose

   pure function cross3(left,right) result(cross)
      real(dblprec), intent(in) :: left(3),right(3)
      real(dblprec) :: cross(3)
      cross=(/left(2)*right(3)-left(3)*right(2),left(3)*right(1)-left(1)*right(3),left(1)*right(2)-left(2)*right(1)/)
   end function cross3

   pure real(dblprec) function determinant3(matrix) result(determinant)
      real(dblprec), intent(in) :: matrix(3,3)
      determinant=matrix(1,1)*(matrix(2,2)*matrix(3,3)-matrix(2,3)*matrix(3,2))-matrix(1,2)*(matrix(2,1)*matrix(3,3)-matrix(2,3)*matrix(3,1))+matrix(1,3)*(matrix(2,1)*matrix(3,2)-matrix(2,2)*matrix(3,1))
   end function determinant3

   pure function inverse3(matrix) result(inverse)
      real(dblprec), intent(in) :: matrix(3,3)
      real(dblprec) :: inverse(3,3),determinant
      determinant=determinant3(matrix)
      inverse(1,1)=(matrix(2,2)*matrix(3,3)-matrix(2,3)*matrix(3,2))/determinant; inverse(1,2)=(matrix(1,3)*matrix(3,2)-matrix(1,2)*matrix(3,3))/determinant; inverse(1,3)=(matrix(1,2)*matrix(2,3)-matrix(1,3)*matrix(2,2))/determinant
      inverse(2,1)=(matrix(2,3)*matrix(3,1)-matrix(2,1)*matrix(3,3))/determinant; inverse(2,2)=(matrix(1,1)*matrix(3,3)-matrix(1,3)*matrix(3,1))/determinant; inverse(2,3)=(matrix(1,3)*matrix(2,1)-matrix(1,1)*matrix(2,3))/determinant
      inverse(3,1)=(matrix(2,1)*matrix(3,2)-matrix(2,2)*matrix(3,1))/determinant; inverse(3,2)=(matrix(1,2)*matrix(3,1)-matrix(1,1)*matrix(3,2))/determinant; inverse(3,3)=(matrix(1,1)*matrix(2,2)-matrix(1,2)*matrix(2,1))/determinant
   end function inverse3

end module MultiChannelCoarseTensorOperator
