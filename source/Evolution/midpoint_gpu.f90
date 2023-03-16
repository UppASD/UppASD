!-------------------------------------------------------------------------------
!  MODULE: Midpoint
!> @brief
!> The semi-implicit midpoint solver for the LLG-equations
!> @details Ref: J.H. Mentink et al, J. Phys.: Condens. Matter, 22, 176001 (2010)
!> @authors
!> Johan Mentink
!> Jonathan Chico
!> Anders Bergman
!> @copyright
!> GNU Public License.
!-------------------------------------------------------------------------------
module Midpoint_gpu

   use Profiling
   use Parameters

   implicit none

   real(dblprec), dimension(:,:,:), allocatable :: btorque_full !< Resulting effective field

   private :: btorque_full

contains

   !----------------------------------------------------------------------------
   !  SUBROUTINE: smodeulermpt
   !> @brief
   !> Semi-implicit midpoint variants. All consist of 2 steps (t and f)
   !> Names according to note Michael Tretyakov
   !> First step of midpoint solver
   !> @note Jonathan Chico: Added Zhang-Li spin transfer torque term into the field
   !----------------------------------------------------------------------------
   subroutine smodeulermpt_gpu(Natom,Mensemble,Landeg,bn,lambda1_array,beff,emom,emom2, &
      emomM,mmom,deltat,thermal_field,STT,do_she,do_sot,btorque,she_btorque,        &
      sot_btorque,Nred,red_atom_list)

      use Constants
      use RandomNumbers, only : ranv

      implicit none

      integer, intent(in) :: Nred            !< Number of moments that can be updated
      integer, intent(in) :: Natom           !< Number of atoms in system
      integer, intent(in) :: Mensemble       !< Number of ensembles
      real(dblprec), intent(in) :: bn        !< Scaling factor for LLG equation (parameter=1.0_dblprec)
      real(dblprec), intent(in) :: deltat    !< Time step
      character(len=1), intent(in) :: STT    !< Treat spin transfer torque?
      character(len=1), intent(in) :: do_she !< Treat the SHE spin transfer torque
      character(len=1), intent(in) :: do_sot !< Treat the general SOT model
      integer, dimension(Nred), intent(in) :: red_atom_list !< Reduced list containing atoms allowed to evolve in a fixed moment calculation
      real(dblprec), dimension(Natom), intent(in) :: Landeg  !< Gyromagnetic ratio
      real(dblprec), dimension(Natom), intent(in) :: lambda1_array   !< Damping parameter
      real(dblprec), dimension(Natom,Mensemble), intent(in) :: mmom  !< Magnitude of magnetic moments
      real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: emom         !< Current unit moment vector
      real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: btorque      !< Spin transfer torque
      real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: she_btorque  !< SHE spin transfer torque
      real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: sot_btorque  !< Spin orbit torque
      ! .. In/out variables
      real(dblprec), dimension(3,Natom,Mensemble), intent(inout) :: beff   !< Total effective field from application of Hamiltonian
      real(dblprec), dimension(3,Natom,Mensemble), intent(inout) :: emom2  !< Final (or temporary) unit moment vector
      real(dblprec), dimension(3,Natom,Mensemble), intent(inout) :: emomM  !< Current magnetic moment vector
      real(dblprec), dimension(3,Natom,Mensemble), intent(inout) :: thermal_field   !< Thermal field

      ! ... Local variables ...
      integer :: i, j,ired
      real(dblprec) :: lldamp

      ! deterministic variables
      real(dblprec),dimension(3) :: a1
      !
      ! stochastic variables
      real(dblprec),dimension(3) :: s1
      !
      ! auxilary variables
      real(dblprec),dimension(3) :: A
      real(dblprec) :: detAi
      real(dblprec),dimension(3) :: a2
      !
      ! time steps
      real(dblprec) :: dt,sqrtdt
      real(dblprec) :: dtg, sqrtdtg
      !
      ! spins
      real(dblprec),dimension(3)  :: et

      btorque_full=0.0_dblprec

      ! Call parameters
      !$omp target enter data map(to:Nred,Natom,Mensemble,bn,deltat,STT,do_she,do_sot)
      !$omp target enter data map(to:red_atom_list,Landeg,lambda1_array,mmom,emom,btorque,she_btorque,sot_btorque)
      !$omp target enter data map(to:beff,emom2,emomM,thermal_field)
      ! Local variables
      !!!$omp target enter data map(to:i,j,ired,lldamp,a1,s1,A,a2)
      !!!$omp target enter data map(to:dt,sqrtdt,dtg,sqrtdtg,et)
      !$omp target teams distribute parallel do simd collapse(2) private(ired,i,j,et,s1,a1,A,detAi,a2,dt,dtg,sqrtdt,sqrtdtg,lldamp)
      do ired=1,Nred
         do j=1,Mensemble
            i=red_atom_list(ired)
            lldamp=1.0_dblprec/(1.0_dblprec+lambda1_array(i)*lambda1_array(i))
            dt=deltat*bn*gama*lldamp !dimm. less time
            sqrtdt=sqrt(dt)
            dtg=dt*Landeg(i)
            sqrtdtg=sqrtdt*Landeg(i)
            ! a1 = -b1 - lambda*(e1 cross b1)
            a1(1)=-btorque_full(1,i,j)-beff(1,i,j)-lambda1_array(i)*(emom(2,i,j)*beff(3,i,j)-emom(3,i,j)*beff(2,i,j))
            a1(2)=-btorque_full(2,i,j)-beff(2,i,j)-lambda1_array(i)*(emom(3,i,j)*beff(1,i,j)-emom(1,i,j)*beff(3,i,j))
            a1(3)=-btorque_full(3,i,j)-beff(3,i,j)-lambda1_array(i)*(emom(1,i,j)*beff(2,i,j)-emom(2,i,j)*beff(1,i,j))
            !
            ! s1 is stochastic counterpart of a1
            s1(1)=-ranv(1,i,j)-lambda1_array(i)*(emom(2,i,j)*ranv(3,i,j)-emom(3,i,j)*ranv(2,i,j))
            s1(2)=-ranv(2,i,j)-lambda1_array(i)*(emom(3,i,j)*ranv(1,i,j)-emom(1,i,j)*ranv(3,i,j))
            s1(3)=-ranv(3,i,j)-lambda1_array(i)*(emom(1,i,j)*ranv(2,i,j)-emom(2,i,j)*ranv(1,i,j))

            thermal_field(:,i,j)=s1(:)
            !
            !
            ! semi-implicitness midpoint requires solution of linear system:
            ! A*e2 = At*e1, At=transpose(A) => e2=inv(A)*At*e1
            ! A = I + skew(dt*a1/2 + sqrt(dt)*s1/2)
            ! write A*e2=a2, a2=At*e1 => e2=inv(A)*a2
            ! Ax,Ay,Az off-diagonal components of A
            ! solve with Cramers' rule => define detAi=1/determinant(A)
            !
            A(:)=0.5_dblprec*dtg*a1(:)+0.5_dblprec*sqrtdtg*s1(:)

            detAi=1.0_dblprec/(1.0_dblprec+sum(A(:)*A(:)))
            !

            a2(1)=emom(1,i,j)+emom(2,i,j)*A(3)-emom(3,i,j)*A(2)
            a2(2)=emom(2,i,j)+emom(3,i,j)*A(1)-emom(1,i,j)*A(3)
            a2(3)=emom(3,i,j)+emom(1,i,j)*A(2)-emom(2,i,j)*A(1)
            !
            et(1)=a2(1)*(1+A(1)*A(1))+a2(2)*(A(1)*A(2)+A(3))+a2(3)*(A(1)*A(3)-A(2))
            et(2)=a2(1)*(A(2)*A(1)-A(3))+a2(2)*(1+A(2)*A(2))+a2(3)*(A(2)*A(3)+A(1))
            et(3)=a2(1)*(A(3)*A(1)+A(2))+a2(2)*(A(3)*A(2)-A(1))+a2(3)*(1+A(3)*A(3))

            ! now use et for writing et'=(e1+et)/2 in emom2
            et(:)=0.5_dblprec*(emom(:,i,j)+et(:)*detAi)
            !
            ! write et'=(e1+et)/2 in emom2
            emom2(:,i,j)=et(:)

            ! write new emomM
            ! effective field for step f approximated with et
            ! no new variables (e.g emomt) to reduce memory requirements
            emomM(:,i,j)=et(:)*mmom(i,j)
         enddo
      enddo
      !$omp end target teams distribute parallel do simd
      return

   end subroutine smodeulermpt_gpu

   !-----------------------------------------------------------------------------
   !  SUBROUTINE: modeulermpf
   !> @brief
   !> Second step of midpoint solver
   !-----------------------------------------------------------------------------
   subroutine modeulermpf_gpu(Natom,Mensemble,Landeg,bn,lambda1_array,beff,emom,emom2,  &
      deltat,STT,do_she,do_sot,btorque,she_btorque,sot_btorque,Nred,red_atom_list)

      use Constants
      use RandomNumbers, only : ranv

      implicit none

      integer, intent(in) :: Nred   !< Number of moments that can be updated
      integer, intent(in) :: Natom  !< Number of atoms in system
      integer, intent(in) :: Mensemble !< Number of ensembles
      real(dblprec), intent(in) :: bn !< Scaling factor for LLG equation (parameter=1.0_dblprec)
      real(dblprec), intent(in) :: deltat !< Time step
      character(len=1), intent(in) :: STT    !< Treat spin transfer torque?
      character(len=1), intent(in) :: do_she !< Treat the SHE spin transfer torque
      character(len=1), intent(in) :: do_sot !< Treat the general SOT model
      integer, dimension(Nred), intent(in) :: red_atom_list !< Reduced list containing atoms allowed to evolve in a fixed moment calculation
      real(dblprec), dimension(Natom), intent(in) :: Landeg  !< Gyromagnetic ratio
      real(dblprec), dimension(Natom), intent(in) :: lambda1_array !< Damping parameter
      real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: beff !< Total effective field from application of Hamiltonian
      real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: btorque      !< Spin transfer torque
      real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: she_btorque  !< SHE spin transfer torque
      real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: sot_btorque  !< Spin orbit torque
      ! .. In/out variables
      real(dblprec), dimension(3,Natom,Mensemble), intent(inout) :: emom   !< Current unit moment vector
      real(dblprec), dimension(3,Natom,Mensemble), intent(inout) :: emom2  !< Final (or temporary) unit moment vector

      ! deterministic variables
      real(dblprec),dimension(3) :: a1
      !
      ! stochastic variables
      real(dblprec),dimension(3) :: s1
      !
      ! auxilary variables
      real(dblprec),dimension(3) :: A
      real(dblprec) :: detAi
      real(dblprec),dimension(3) :: a2
      !
      ! time steps
      real(dblprec) :: dt,sqrtdt
      real(dblprec) :: dtg, sqrtdtg
      !
      ! spins
      real(dblprec),dimension(3) :: etp
      !
      ! ... Local variables ...
      integer :: i, j, ired
      real(dblprec) :: lldamp

      btorque_full=0.0_dblprec

      ! Call parameters
      !$omp target enter data map(to:Nred,Natom,Mensemble,bn,deltat,STT,do_she,do_sot)
      !!!$omp target enter data map(to:red_atom_list,Landeg,lambda1_array,mmom,emom,btorque,she_btorque,sot_btorque)
      !!!$omp target enter data map(to:beff,emom2,emomM,thermal_field)
      ! Local variables
      !!!$omp target enter data map(to:i,j,ired,lldamp,a1,s1,A,a2)
      !!!$omp target enter data map(to:dt,sqrtdt,dtg,sqrtdtg,et)
      !$omp target teams distribute parallel do simd collapse(2) private(ired,i,j,etp,s1,a1,A,detAi,a2,dt,dtg,sqrtdt,sqrtdtg,lldamp)
      do ired=1,Nred
         do j=1,Mensemble
            i=red_atom_list(ired)
            lldamp=1.0_dblprec/(1.0_dblprec+lambda1_array(i)*lambda1_array(i))
            dt=deltat*bn*gama*lldamp !dimm. less time
            sqrtdt=sqrt(dt)
            dtg=dt*Landeg(i)
            sqrtdtg=sqrtdt*Landeg(i)
            etp(:)=emom2(:,i,j) ! load etp=et' back from emom2
            !
            ! a1 = -b1 - lambda*(et cross b1)
            a1(1)=-btorque_full(1,i,j)-beff(1,i,j)-lambda1_array(i)*(etp(2)*beff(3,i,j)-etp(3)*beff(2,i,j))
            a1(2)=-btorque_full(2,i,j)-beff(2,i,j)-lambda1_array(i)*(etp(3)*beff(1,i,j)-etp(1)*beff(3,i,j))
            a1(3)=-btorque_full(3,i,j)-beff(3,i,j)-lambda1_array(i)*(etp(1)*beff(2,i,j)-etp(2)*beff(1,i,j))
            !
            ! s1 is stochastic counterpart of a1
            s1(1)=-ranv(1,i,j)-lambda1_array(i)*(etp(2)*ranv(3,i,j)-etp(3)*ranv(2,i,j))
            s1(2)=-ranv(2,i,j)-lambda1_array(i)*(etp(3)*ranv(1,i,j)-etp(1)*ranv(3,i,j))
            s1(3)=-ranv(3,i,j)-lambda1_array(i)*(etp(1)*ranv(2,i,j)-etp(2)*ranv(1,i,j))
            !
            ! semi-implicitness midpoint requires solution of linear system:
            ! A*e2 = At*e1, At=transpose(A) => e2=inv(A)*At*e1
            ! A = I + skew(dt*a1/2 + sqrt(dt)*s1/2)
            ! write A*e2=a2, a2=At*e1 => e2=inv(A)*a2
            ! Ax,Ay,Az off-diagonal components of A
            ! solve with Cramers' rule => define detAi=1/determinant(A)
            !
            A(:)=0.5_dblprec*dtg*a1(:)+0.5_dblprec*sqrtdtg*s1(:)
            detAi=1.0_dblprec/(1.0_dblprec+sum(A(:)*A(:)))
            !
            a2(1)=emom(1,i,j)+emom(2,i,j)*A(3)-emom(3,i,j)*A(2)
            a2(2)=emom(2,i,j)+emom(3,i,j)*A(1)-emom(1,i,j)*A(3)
            a2(3)=emom(3,i,j)+emom(1,i,j)*A(2)-emom(2,i,j)*A(1)

            ! now use etp simply to write emom2
            etp(1)=a2(1)*(1+A(1)*A(1))+a2(2)*(A(1)*A(2)+A(3))+a2(3)*(A(1)*A(3)-A(2))
            etp(2)=a2(1)*(A(2)*A(1)-A(3))+a2(2)*(1+A(2)*A(2))+a2(3)*(A(2)*A(3)+A(1))
            etp(3)=a2(1)*(A(3)*A(1)+A(2))+a2(2)*(A(3)*A(2)-A(1))+a2(3)*(1+A(3)*A(3))
            !
            emom2(:,i,j)=etp(:)*detAi

         end do
      end do
      !$omp end target teams distribute parallel do simd

      return

   end subroutine modeulermpf_gpu


   !----------------------------------------------------------------------------
   ! SUBROUTINE: allocate_aux_midpoint_fields
   !> @brief Allocation of auxilary fields for the treatment of STT and SOT based torques
   !> @author Jonathan Chico
   !----------------------------------------------------------------------------
   subroutine allocate_aux_midpoint_fields_gpu(flag,Natom,Mensemble)

      implicit none

      integer, intent(in) :: flag   !< Allocate or deallocate (1/-1)
      integer, intent(in), optional :: Natom !< Number of atoms in system
      integer, intent(in), optional :: Mensemble   !< Number of ensembles

      integer :: i_stat,i_all

      if (flag>0) then
         allocate(btorque_full(3,Natom,Mensemble),stat=i_stat)
         call memocc(i_stat,product(shape(btorque_full))*kind(btorque_full),'btorque_full','allocate_aux_midpoint_fields')
         btorque_full=0.0_dblprec
      else
         i_all=-product(shape(btorque_full))*kind(btorque_full)
         deallocate(btorque_full,stat=i_stat)
         call memocc(i_stat,i_all,'btorque_full','allocate_aux_midpoint_fields')
      endif

   end subroutine allocate_aux_midpoint_fields_gpu

end module midpoint_gpu
