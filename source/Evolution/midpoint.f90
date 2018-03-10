!-------------------------------------------------------------------------------
!  MODULE: Midpoint
!> @brief
!> The semi-implicit midpoint solver for the LLG-equations
!> @details Ref: J.H. Mentink et al, J. Phys.: Condens. Matter, 22, 176001 (2010)
!> @author
!> Johan Mentink
!-------------------------------------------------------------------------------
module Midpoint
   use Parameters
   use Profiling

contains


   !-----------------------------------------------------------------------------
   !> SUBROUTINE: smodeulermpt
   !> @brief
   !> Semi-implicit midpoint variants
   !> All consist of 2 steps (t and f)
   !> Names according to note Michael Tretyakov
   !> First step of midpoint solver
   !-----------------------------------------------------------------------------
   subroutine smodeulermpt(Natom, Mensemble, Landeg,bn, lambda1_array, beff, emom, &
         emom2, emomM, mmom, deltat,thermal_field)

      use Constants
      use RandomNumbers, only : ranv

      implicit none

      integer, intent(in) :: Natom !< Number of atoms in system
      integer, intent(in) :: Mensemble !< Number of ensembles
      real(dblprec), dimension(Natom), intent(in) :: Landeg  !< Gyromagnetic ratio
      real(dblprec), intent(in) :: bn !< Scaling factor for LLG equation (parameter=1.0d0)
      real(dblprec), dimension(Natom), intent(in) :: lambda1_array !< Damping parameter
      real(dblprec), dimension(3,Natom,Mensemble), intent(inout) :: beff !< Total effective field from application of Hamiltonian
      real(dblprec), dimension(3,Natom,Mensemble), intent(inout) :: thermal_field
      real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: emom   !< Current unit moment vector
      real(dblprec), dimension(3,Natom,Mensemble), intent(inout) :: emom2  !< Final (or temporary) unit moment vector
      real(dblprec), dimension(3,Natom,Mensemble), intent(inout) :: emomM  !< Current magnetic moment vector
      real(dblprec), dimension(Natom,Mensemble), intent(in) :: mmom !< Magnitude of magnetic moments
      real(dblprec), intent(in) :: deltat !< Time step

      ! ... Local variables ...
      integer :: i, j, ij
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

      !$omp parallel do default(shared) schedule(static) &
      !$omp private(i,j,et,s1,a1,A,detAi,a2,dt,dtg,sqrtdt,sqrtdtg,lldamp),collapse(2)
      do i=1,Natom
         do j=1,Mensemble

            lldamp=1.0D0/(1.0D0+lambda1_array(i)*lambda1_array(i))
            dt=deltat*bn*gama*lldamp !dimm. less time
            sqrtdt=sqrt(dt)
            dtg=dt*Landeg(i)
            sqrtdtg=sqrtdt*Landeg(i)

            ! a1 = -b1 - lambda*(e1 cross b1)
            a1(1)=-beff(1,i,j)-lambda1_array(i)*(emom(2,i,j)*beff(3,i,j)-emom(3,i,j)*beff(2,i,j))
            a1(2)=-beff(2,i,j)-lambda1_array(i)*(emom(3,i,j)*beff(1,i,j)-emom(1,i,j)*beff(3,i,j))
            a1(3)=-beff(3,i,j)-lambda1_array(i)*(emom(1,i,j)*beff(2,i,j)-emom(2,i,j)*beff(1,i,j))
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
            A(:)=0.5d0*dtg*a1(:)+0.5d0*sqrtdtg*s1(:)

            detAi=1.0d0/(1.0d0+sum(A(:)*A(:)))
            !
            a2(1)=emom(1,i,j)+emom(2,i,j)*A(3)-emom(3,i,j)*A(2)
            a2(2)=emom(2,i,j)+emom(3,i,j)*A(1)-emom(1,i,j)*A(3)
            a2(3)=emom(3,i,j)+emom(1,i,j)*A(2)-emom(2,i,j)*A(1)
            !
            et(1)=a2(1)*(1+A(1)*A(1))+a2(2)*(A(1)*A(2)+A(3))+a2(3)*(A(1)*A(3)-A(2))
            et(2)=a2(1)*(A(2)*A(1)-A(3))+a2(2)*(1+A(2)*A(2))+a2(3)*(A(2)*A(3)+A(1))
            et(3)=a2(1)*(A(3)*A(1)+A(2))+a2(2)*(A(3)*A(2)-A(1))+a2(3)*(1+A(3)*A(3))

            ! now use et for writing et'=(e1+et)/2 in emom2
            et(:)=0.5d0*(emom(:,i,j)+et(:)*detAi)
            !
            ! write et'=(e1+et)/2 in emom2
            emom2(:,i,j)=et(:)

            ! write new emomM
            ! effective field for step f approximated with et
            ! no new variables (e.g emomt) to reduce memory requirements
            emomM(:,i,j)=et(:)*mmom(i,j)
         enddo
      enddo
      !$omp end parallel do
      return

   end subroutine smodeulermpt


   !-----------------------------------------------------------------------------
   !> SUBROUTINE: modeulermpf
   !> @brief
   !> Second step of midpoint solver
   !-----------------------------------------------------------------------------
   subroutine modeulermpf(Natom, Mensemble, Landeg, bn, lambda1_array, beff, emom, emom2, deltat)
      use Constants
      use RandomNumbers, only : ranv

      implicit none

      integer, intent(in) :: Natom !< Number of atoms in system
      integer, intent(in) :: Mensemble !< Number of ensembles
      real(dblprec), dimension(Natom), intent(in) :: Landeg  !< Gyromagnetic ratio
      real(dblprec), intent(in) :: bn !< Scaling factor for LLG equation (parameter=1.0d0)
      real(dblprec), dimension(Natom), intent(in) :: lambda1_array !< Damping parameter
      real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: beff !< Total effective field from application of Hamiltonian
      real(dblprec), dimension(3,Natom,Mensemble), intent(inout) :: emom   !< Current unit moment vector
      real(dblprec), dimension(3,Natom,Mensemble), intent(inout) :: emom2  !< Final (or temporary) unit moment vector
      real(dblprec), intent(in) :: deltat !< Time step

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
      integer :: i, j, ij
      real(dblprec) :: lldamp
      !
      !$omp parallel do default(shared) schedule(static) &
      !$omp private(i,j,etp,s1,a1,A,detAi,a2,dt,dtg,sqrtdt,sqrtdtg,lldamp),collapse(2)

      do i=1,Natom
         do j=1,Mensemble
            lldamp=1.0D0/(1.0D0+lambda1_array(i)*lambda1_array(i))
            dt=deltat*bn*gama*lldamp !dimm. less time
            sqrtdt=sqrt(dt)
            dtg=dt*Landeg(i)
            sqrtdtg=sqrtdt*Landeg(i)
            etp(:)=emom2(:,i,j) ! load etp=et' back from emom2
            !
            ! a1 = -b1 - lambda*(et cross b1)
            a1(1)=-beff(1,i,j)-lambda1_array(i)*(etp(2)*beff(3,i,j)-etp(3)*beff(2,i,j))
            a1(2)=-beff(2,i,j)-lambda1_array(i)*(etp(3)*beff(1,i,j)-etp(1)*beff(3,i,j))
            a1(3)=-beff(3,i,j)-lambda1_array(i)*(etp(1)*beff(2,i,j)-etp(2)*beff(1,i,j))
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
            A(:)=0.5d0*dtg*a1(:)+0.5d0*sqrtdtg*s1(:)
            detAi=1.0d0/(1.0d0+sum(A(:)*A(:)))
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
      !$omp end parallel do

      return

   end subroutine modeulermpf


end module midpoint
