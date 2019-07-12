!------------------------------------------------------------------------------------
! MODULE: LLGI
!> @brief Heun like solver for the stochastic LLGI-equation
!> @details Heun type solver for the generalized LLGI equation
!> with magnetic accleration terms. If relaxtime is put to zero it correponds to standard LLG eq.
!> @todo Replace unit length moment vectors emom with full lenght vector emomM
!> @author
!> Lars Bergqvist
!> @copyright
!> GNU Public License.
!------------------------------------------------------------------------------------
module LLGI
   use Parameters
   use Profiling

   implicit none
   !
   real(dblprec), dimension(:,:,:), allocatable :: btherm !< Thermal stochastic field
   real(dblprec), dimension(:,:,:), allocatable :: mtmp  !< Temporary new solution m(t+dt)
   real(dblprec), dimension(:,:,:), allocatable :: mpast  !< Saving past solution m(t-dt)
   real(dblprec), dimension(:,:,:), allocatable :: bdup !< Resulting effective field

   private

   public :: llgi_evolve_first, llgi_evolve_second, allocate_llgifields

contains

   !---------------------------------------------------------------------------------
   ! SUBROUTINE: llgi_evolve_first
   !> @brief First step (predictor) of the Heun solver for the LLG-I eq, calculates the stochastic field preform update
   !---------------------------------------------------------------------------------
   subroutine llgi_evolve_first(Natom,Mensemble,lambda1_array,beff,b2eff,emom,emom2,&
      emomM,mmom,delta_t,relaxtime,Temp_array,temprescale,thermal_field,Nred,       &
      red_atom_list)
      !
      implicit none
      !
      integer, intent(in) :: Nred   !< Number of moments that can be updated
      integer, intent(in) :: Natom  !< Number of atoms in system
      integer, intent(in) :: Mensemble !< Number of ensembles
      real(dblprec), intent(in) :: delta_t !< Time step
      real(dblprec), intent(in) :: relaxtime !< Relaxation time for inertial regime
      integer, dimension(Nred), intent(in) :: red_atom_list !< Reduced list containing atoms allowed to evolve in a fixed moment calculation
      real(dblprec), dimension(Natom), intent(in) :: Temp_array      !< Temperature (array)
      real(dblprec), dimension(Natom), intent(in) :: lambda1_array   !< Damping parameter
      real(dblprec), dimension(Natom,Mensemble), intent(in) :: mmom     !< Magnitude of magnetic moments
      real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: beff   !< Total effective field from application of Hamiltonian
      real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: emom2  !< Final (or temporary) unit moment vector m(t)
      real(dblprec), dimension(3,Natom,Mensemble), intent(out) :: b2eff !< Temporary storage of magnetic field
      real(dblprec), dimension(3,Natom,Mensemble), intent(out) :: emom  !< Current unit moment vector m(t),mp(t+1)
      real(dblprec), dimension(3,Natom,Mensemble), intent(out) :: emomM !< Current magnetic moment vector
      real(dblprec), intent(inout) :: temprescale  !< Temperature rescaling from QHB
      real(dblprec), dimension(3,Natom,Mensemble), intent(inout) :: thermal_field

      integer :: i,k

      ! LLGI equation, Heun like predictor-corrector solver
      ! Calculate stochastic field
      call thermfield(Natom, Mensemble,lambda1_array,mmom,delta_t,Temp_array,temprescale)

      thermal_field=btherm
      ! Construct field
      bdup=beff+btherm

      ! Load back solution m(t-dt) from previous timestep to emom
      emom=mpast

      ! Predictor step of Heun solver
      call llgiheun(Natom,Mensemble,lambda1_array,delta_t,relaxtime,emom,emom2,Nred,&
         red_atom_list)

      ! Copy predicted m(t+dt) to emom for heisge and b(t) to b2eff
      ! At this stage, m(t-dt)=mpast,m(t)=emom2,mp(t+1)=emom=mtmp
      emom=mtmp
      !$omp parallel do schedule(static) default(shared) private(i,k)
      do i=1,Natom
         do k=1,Mensemble
            emomM(:,i,k)=mtmp(:,i,k)*mmom(i,k)
         end do
      end do
      !$omp end parallel do

      b2eff=bdup
   end subroutine llgi_evolve_first

   !---------------------------------------------------------------------------------
   ! SUBROUTINE: llgi_evolve_second
   !> @brief Second step of solver, calculates the corrected effective field from
   !> the predicted effective fields.
   !---------------------------------------------------------------------------------
   subroutine llgi_evolve_second(Natom,Mensemble,lambda1_array,beff,b2eff,emom,     &  
      emom2,delta_t,relaxtime,Nred,red_atom_list)
      !
      implicit none
      !
      integer, intent(in) :: Nred	!< Number of moments that can be updated
      integer, intent(in) :: Natom !< Number of atoms in system
      integer, intent(in) :: Mensemble !< Number of ensembles
      integer, dimension(Nred), intent(in) :: red_atom_list !< Reduced list containing atoms allowed to evolve in a fixed moment calculation
      real(dblprec), dimension(Natom), intent(in) :: lambda1_array !< Damping parameter
      real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: beff !< Total effective field from application of Hamiltonian
      real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: b2eff !< Temporary storage of magnetic field
      real(dblprec), dimension(3,Natom,Mensemble), intent(out) :: emom   !< Current unit moment vector
      real(dblprec), dimension(3,Natom,Mensemble), intent(inout) :: emom2  !< Final (or temporary) unit moment vector
      real(dblprec), intent(in) :: delta_t !< Time step
      real(dblprec), intent(in) :: relaxtime !< Relaxation time for inertial regime
      !
      integer :: i,k
      real(dblprec) :: etot

      ! Construct corrected field
      bdup=0.5_dblprec*(beff+btherm)+0.5_dblprec*b2eff

      ! restore m(t-1) to emom from mpast and update m(t) to mpast
      emom=mpast
      mpast=emom2

      ! Corrector step of Heun solver
      call llgiheun(Natom, Mensemble,lambda1_array,delta_t,relaxtime,emom,emom2,    &
         Nred,red_atom_list)

      ! Save m(t) to emom (can be removed)
      emom=emom2

      ! Copy corrected solution m(t+dt) to emom2 and normalize solution
      emom2=mtmp

      !$omp parallel do schedule(static) default(shared) private(i,k,etot)
      do i=1,Natom
         do k=1,Mensemble
            etot=1.0_dblprec/norm2(emom2(:,i,k))
            emom2(1,i,k)=emom2(1,i,k)*etot
            emom2(2,i,k)=emom2(2,i,k)*etot
            emom2(3,i,k)=emom2(3,i,k)*etot
         end do
      end do
      !$omp end parallel do

   end subroutine llgi_evolve_second

   !-----------------------------------------------------------------------------
   ! SUBROUTINE: thermfield
   !> @brief Calculates stochastic field
   !> @todo Confirm the magnitude of the stochastic field as a function
   !> of temperature
   !-----------------------------------------------------------------------------
   subroutine thermfield(Natom, Mensemble, lambda1_array, mmom, deltat,Temp_array,temprescale)
      !
      use Constants, only : k_bolt, gama, mub
      use RandomNumbers, only : rng_gaussian

      implicit none

      integer, intent(in) :: Natom !< Number of atoms in system
      integer, intent(in) :: Mensemble !< Number of ensembles
      real(dblprec), dimension(Natom), intent(in) :: lambda1_array !< Damping parameter
      real(dblprec), dimension(Natom,Mensemble), intent(in) :: mmom !< Magnitude of magnetic moments
      real(dblprec), intent(in) :: deltat !< Time step
      real(dblprec), dimension(Natom) ,intent(in) :: Temp_array  !< Temperature (array)
      real(dblprec), intent(in) :: temprescale  !< Temperature rescaling from QHB

      real(dblprec), dimension(Natom) :: Dp
      real(dblprec) :: mu, sigma

      integer :: i,k

      ! LL equations ONE universal damping
      Dp=(2._dblprec*lambda1_array*k_bolt)/(deltat*gama*mub)    !Ok (LLG)  (still needs to be checked with Fokker-Planck)

      call rng_gaussian(btherm,3*Natom*Mensemble,1.0_dblprec)

      ! LLG equations ONE universal damping
      mu=0_dblprec
      do k=1, Mensemble
         do i=1, Natom
            sigma=sqrt(Dp(i)*temprescale*Temp_array(i)/mmom(i,k))
            btherm(1,i,k)=btherm(1,i,k)*sigma
            btherm(2,i,k)=btherm(2,i,k)*sigma
            btherm(3,i,k)=btherm(3,i,k)*sigma
         end do
      end do

   end subroutine thermfield

   !---------------------------------------------------------------------------------
   ! SUBROUTINE: llgiheun
   !> @brief Heun solver for LLG-I equation
   !---------------------------------------------------------------------------------
   subroutine llgiheun(Natom, Mensemble,lambda1_array,delta_t,relaxtime,emom,emom2, &
      Nred,red_atom_list)
      !
      use Constants, only : gama
      implicit none

      integer, intent(in) :: Nred   !< Number of moments that can be updated
      integer, intent(in) :: Natom !< Number of atoms in system
      integer, intent(in) :: Mensemble !< Number of ensembles
      real(dblprec), intent(in) :: delta_t !< Time step
      real(dblprec), intent(in) :: relaxtime !< Relaxation time for inertial regime
      integer, dimension(Nred), intent(in) :: red_atom_list !< Reduced list containing atoms allowed to evolve in a fixed moment calculation
      real(dblprec), dimension(Natom), intent(in) :: lambda1_array !< Damping parameter
      real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: emom   !< Current unit moment vector, past timestep
      real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: emom2  !< Final (or temporary) unit moment vector
      ! ..Local variables
      real(dblprec), dimension(3,3) :: A,Ainv
      real(dblprec),dimension(3) :: B   ! b vector
      real(dblprec) :: detmatrix,dtg
      real(dblprec) :: e1x, e1y, e1z
      real(dblprec) :: etx, ety, etz
      real(dblprec) :: hx, hy, hz
      real(dblprec), dimension(Natom) :: alfa, beta, zeta

      integer :: i,k,ired

      alfa=lambda1_array
      beta=-alfa*(1._dblprec+relaxtime/delta_t)
      zeta=-alfa*(relaxtime/delta_t)
      dtg=-delta_t*gama
      !$omp parallel do default(shared) &
      !$omp private(i,ired,k,etx,ety,etz,e1x,e1y,e1z,hx,hy,hz,A &
      !$omp ,B,detmatrix,Ainv)

      do ired=1, Nred
         do k=1, Mensemble
            i=red_atom_list(ired)
            e1x=emom2(1,i,k)   !m(i)
            e1y=emom2(2,i,k)
            e1z=emom2(3,i,k)

            etx=emom(1,i,k)    !m(i-1)
            ety=emom(2,i,k)
            etz=emom(3,i,k)

            hx=bdup(1,i,k)   !heff(i)
            hy=bdup(2,i,k)
            hz=bdup(3,i,k)

            ! Construct A matrix
            A(1,1)=1.0_dblprec
            A(1,2)=-beta(i)*e1z
            A(1,3)=beta(i)*e1y
            A(2,1)=beta(i)*e1z
            A(2,2)=1.0_dblprec
            A(2,3)=-beta(i)*e1x
            A(3,1)=-beta(i)*e1y
            A(3,2)=beta(i)*e1x
            A(3,3)=1.0_dblprec

            ! Construct inversion of A
            detmatrix=1.0_dblprec-A(2,3)*A(3,2)-A(1,2)*A(2,1)-A(1,3)*A(3,1)
            detmatrix=1.0_dblprec/detmatrix
            Ainv(1,1)=(1.0_dblprec-A(3,2)*A(2,3))*detmatrix
            Ainv(1,2)=(A(1,3)*A(3,2)-A(1,2))*detmatrix
            Ainv(1,3)=(A(1,2)*A(2,3)-A(1,3))*detmatrix
            Ainv(2,1)=(A(2,3)*A(3,1)-A(2,1))*detmatrix
            Ainv(2,2)=(1.0_dblprec-A(3,1)*A(1,3))*detmatrix
            Ainv(2,3)=(A(1,3)*A(2,1)-A(2,3))*detmatrix
            Ainv(3,1)=(A(2,1)*A(3,2)-A(3,1))*detmatrix
            Ainv(3,2)=(A(1,2)*A(3,1)-A(3,2))*detmatrix
            Ainv(3,3)=(1.0_dblprec-A(2,1)*A(1,2))*detmatrix

            ! Construct B vector
            B(1)=e1x+dtg*(e1y*hz-e1z*hy)+zeta(i)*(e1z*ety-e1y*etz)
            B(2)=e1y+dtg*(e1z*hx-e1x*hz)+zeta(i)*(e1x*etz-e1z*etx)
            B(3)=e1z+dtg*(e1x*hy-e1y*hx)+zeta(i)*(e1y*etx-e1x*ety)

            mtmp(1,i,k)=Ainv(1,1)*B(1)+Ainv(1,2)*B(2)+Ainv(1,3)*B(3)
            mtmp(2,i,k)=Ainv(2,1)*B(1)+Ainv(2,2)*B(2)+Ainv(2,3)*B(3)
            mtmp(3,i,k)=Ainv(3,1)*B(1)+Ainv(3,2)*B(2)+Ainv(3,3)*B(3)
         end do
      end do
      !$omp end parallel do
   end subroutine llgiheun

   !-----------------------------------------------------------------------------
   ! SUBROUTINE: allocate_llgifields
   !> @brief Allocates work arrays for the Heun-I solver
   !-----------------------------------------------------------------------------
   subroutine allocate_llgifields(Natom,Mensemble,flag)

      implicit none

      integer, intent(in) :: Natom !< Number of atoms in system
      integer, intent(in) :: Mensemble !< Number of ensembles
      integer, intent(in) :: flag !< Allocate or deallocate (1/-1)

      integer :: i_all, i_stat

      if(flag>0) then
         allocate(mtmp(3,Natom,Mensemble),stat=i_stat)
         call memocc(i_stat,product(shape(mtmp))*kind(mtmp),'mtmp','allocate_llgifields')
         allocate(mpast(3,Natom,Mensemble),stat=i_stat)
         call memocc(i_stat,product(shape(mpast))*kind(mpast),'mpast','allocate_depondtfields')
         mpast=0.0_dblprec
         allocate(btherm(3,Natom,Mensemble),stat=i_stat)
         call memocc(i_stat,product(shape(btherm))*kind(btherm),'btherm','allocate_llgifields')
         allocate(bdup(3,Natom,Mensemble),stat=i_stat)
         call memocc(i_stat,product(shape(bdup))*kind(bdup),'bdup','allocate_llgifields')
      else
         i_all=-product(shape(mtmp))*kind(mtmp)
         deallocate(mtmp,stat=i_stat)
         call memocc(i_stat,i_all,'mtmp','allocate_systemdata')
         i_all=-product(shape(mpast))*kind(mpast)
         deallocate(mpast,stat=i_stat)
         call memocc(i_stat,i_all,'mpast','allocate_systemdata')
         i_all=-product(shape(btherm))*kind(btherm)
         deallocate(btherm,stat=i_stat)
         call memocc(i_stat,i_all,'btherm','allocate_systemdata')
         i_all=-product(shape(bdup))*kind(bdup)
         deallocate(bdup,stat=i_stat)
         call memocc(i_stat,i_all,'bdup','allocate_systemdata')
      end if

   end subroutine allocate_llgifields

end module LLGI
