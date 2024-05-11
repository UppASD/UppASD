!-------------------------------------------------------------------------------
!  MODULE: Depondt_ms
!> @brief
!> The Depondt solver for the stochastic LLG-equation for Multiscale
!> @authors
!> Edgar Mendez
!> Nikos Ntallis
!> Anders Bergman
!> Manuel Pereiro 
!> @copyright
!> GNU Public License.
!> @details In principle the solver is of Heun type but uses rotations to
!> keep the magnitudes of the moments.
!> Ref: Ph. Depondt and F.G. Mertens, J. Phys.: Condens. Matter 21, 336005 (2009)
!> @todo Replace unit length moment vectors emom with full lenght vector emomM
!-------------------------------------------------------------------------------
module Depondt_ms
   use MultiscaleDampingBand
   use Multiscale, only : multiscaleBackbuffer, multiscaleBackbufferHead 
!  use Depondt
   use Profiling
   use Parameters
   use Profiling
   !
   implicit none

   real(dblprec), dimension(:,:,:), allocatable :: mrod !< Rotated magnetic moments
   real(dblprec), dimension(:,:,:), allocatable :: btherm !< Thermal stochastic field
   real(dblprec), dimension(:,:,:), allocatable :: bloc  !< Local effective field
   real(dblprec), dimension(:,:,:), allocatable :: bdup !< Resulting effective field
   real(dblprec), dimension(:,:,:), allocatable :: dedt


!   abstract interface
!     function Dmdt(atom,ensemble) result(d)
!       import :: dblprec
!       integer, intent(in) :: atom,ensemble
!       real(dblprec), dimension(3) :: d
!     end function Dmdt
!   end interface

   private

   public :: depondt_evolve_first_ms, depondt_evolve_second_ms

contains

   subroutine depondt_evolve_first_ms(Natom,Nred,Mensemble,lambda1_array,beff,b2eff,   &
         btorque, emom, emom2, emomM, mmom, delta_t, Temp_array, temprescale,stt,      &
         thermal_field,do_she,she_btorque,do_sot,sot_btorque,red_atom_list,dband)
      !
      use Constants, only : k_bolt, gama, mub
      use RandomNumbers, only : rng_gaussian, rng_gaussianP

      implicit none
      !
      integer, intent(in) :: Nred            !< Number of atoms that evolve
      integer, intent(in) :: Natom           !< Number of atoms in system
      integer, intent(in) :: Mensemble       !< Number of ensembles
      real(dblprec), intent(in) :: delta_t   !< Time step
      character(len=1), intent(in) :: STT    !< Treat spin transfer torque?
      character(len=1), intent(in) :: do_she !< Treat the spin hall effect transfer torque
      character(len=1), intent(in) :: do_sot !< Treat the general SOT model
      integer, dimension(Nred), intent(in) :: red_atom_list !< List of indices of atoms that evolve
      real(dblprec), dimension(Natom), intent(in) :: Temp_array !< Temperature (array)
      real(dblprec), dimension(Natom), intent(in) :: lambda1_array !< Damping parameter
      real(dblprec), dimension(Natom,Mensemble), intent(in) :: mmom !< Magnitude of magnetic moments
      real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: beff !< Total effective field from application of Hamiltonian
      real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: btorque     !< Spin transfer torque
      real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: she_btorque !< Spin Hall effect spin transfer torque
      real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: sot_btorque !< Spin orbit torque

      real(dblprec), dimension(3,Natom,Mensemble), intent(out) :: b2eff !< Temporary storage of magnetic field
      real(dblprec), dimension(3,Natom,Mensemble), intent(out) :: emom2 !< Final (or temporary) unit moment vector
      real(dblprec), dimension(3,Natom,Mensemble), intent(out) :: emomM !< Current magnetic moment vector
      real(dblprec), intent(inout) :: temprescale  !< Temperature rescaling from QHB
      real(dblprec), dimension(3,Natom,Mensemble), intent(inout) :: emom     !< Current unit moment vector
      real(dblprec), dimension(3,Natom,Mensemble), intent(inout) :: thermal_field
      type(DampingBandData), intent(in) :: dband !< Multiscale damping band parameters real(dblprec), dimension(3,Natom,Mensemble), intent(inout) :: thermal_field

      integer :: i,k,ired
      real(dblprec) :: v,Bnorm,hx,hy,hz
      real(dblprec) :: u,cosv,sinv,lldamp
      real(dblprec) :: sigma, Dp, she_fac, stt_fac,sot_fac
      real(dblprec), dimension(3,Natom,Mensemble) :: bdup !< Resulting effective field
      bdup=0.0_dblprec

      if(stt/='N') then
         stt_fac=1.0_dblprec
         !!!$omp parallel do default(shared) private(i,ired,k)  schedule(static) collapse(2)
         do k=1,Mensemble
            do i=1,Natom
               !i=red_atom_list(ired)
               ! Adding STT and SHE torques if present (prefactor instead of if-statement)
               bdup(:,i,k)=bdup(:,i,k)+stt_fac*btorque(:,i,k)
            end do
         end do
         !!!$omp end parallel do
      else
         stt_fac=0.0_dblprec
      end if

      if(do_she/='N') then
         she_fac=1.0_dblprec
         !!!$omp parallel do default(shared) private(i,ired,k)  schedule(static) collapse(2)
         do k=1,Mensemble
            do i=1,Natom
               !i=red_atom_list(ired)
               ! Adding STT and SHE torques if present (prefactor instead of if-statement)
               bdup(:,i,k)= bdup(:,i,k)+she_fac*she_btorque(:,i,k)
            end do
         end do
         !!!$omp end parallel do
      else
         she_fac=0.0_dblprec
      end if
      if(do_sot/='N') then
         sot_fac=1.0_dblprec
         !!!$omp parallel do default(shared) private(i,ired,k)  schedule(static) collapse(2)
         do k=1,Mensemble
            do i=1,Natom
               !i=red_atom_list(ired)
               ! Adding STT, SHE and SOT torques if present (prefactor instead of if-statement)
               bdup(:,i,k)= bdup(:,i,k)+sot_fac*sot_btorque(:,i,k)
            end do
         end do
         !!!$omp end parallel do
      else
         sot_fac=0.0_dblprec
      end if

      call rng_gaussianP(btherm,3*Natom*Mensemble,1.0_dblprec)

      if(.not.allocated(dedt)) then
         allocate(dedt(3,Natom,Mensemble))
      end if

      ! Dupont recipe J. Phys.: Condens. Matter 21 (2009) 336005
      !!!$omp parallel do default(shared) private(ired,i,k,Dp,sigma,Bnorm,hx,hy,hz,v,lldamp,cosv,sinv,u)  schedule(static) collapse(2)
      do k=1,Mensemble
         do i=1,Natom

            !i=red_atom_list(ired)
            ! Thermal field
            !   LL equations ONE universal damping
            Dp=(2.0_dblprec*lambda1_array(i)*k_bolt)/(delta_t*gama*mub)   !LLG
            sigma=sqrt(Dp*temprescale*Temp_array(i)/mmom(i,k))
            btherm(:,i,k)=btherm(:,i,k)*sigma

            ! Construct local field
            bloc(:,i,k)=beff(:,i,k)+btherm(:,i,k)
            thermal_field(:,i,k)=btherm(:,i,k)

            ! Construct effective field (including damping term)
            bdup(1,i,k)=bdup(1,i,k)+bloc(1,i,k)+lambda1_array(i)*emom(2,i,k)*bloc(3,i,k)-lambda1_array(i)*emom(3,i,k)*bloc(2,i,k)
            bdup(2,i,k)=bdup(2,i,k)+bloc(2,i,k)+lambda1_array(i)*emom(3,i,k)*bloc(1,i,k)-lambda1_array(i)*emom(1,i,k)*bloc(3,i,k)
            bdup(3,i,k)=bdup(3,i,k)+bloc(3,i,k)+lambda1_array(i)*emom(1,i,k)*bloc(2,i,k)-lambda1_array(i)*emom(2,i,k)*bloc(1,i,k)

            dedt(:,i,k) =  diff_first(i,k)
         enddo
      enddo

      call buildbeff_ms(Natom, Mensemble,lambda1_array,emom,btorque,stt,do_she,&
         she_btorque,dband,dedt)
         !she_btorque,dband,diff_first)

      !!!$omp parallel do default(shared) private(ired,i,k,Dp,sigma,Bnorm,hx,hy,hz,v,lldamp,cosv,sinv,u)  schedule(static) collapse(2)
      do k=1,Mensemble
         do i=1,Natom

            ! Set up rotation matrices and perform rotations
            lldamp=1.0_dblprec/(1.0_dblprec+lambda1_array(i)**2)
            Bnorm=bdup(1,i,k)**2+bdup(2,i,k)**2+bdup(3,i,k)**2
            Bnorm=sqrt(Bnorm)+1.0d-15
            hx=bdup(1,i,k)/Bnorm
            hy=bdup(2,i,k)/Bnorm
            hz=bdup(3,i,k)/Bnorm
            ! Euler
            !v=0.0_dblprec
            ! Heun
            v=Bnorm*delta_t*gama*lldamp
            ! Ralston
            !v=Bnorm*delta_t*gama*lldamp*2.0_dblprec/3.0_dblprec
            ! Midpoint
            !v=Bnorm*delta_t*gama*lldamp*0.5_dblprec
            cosv=cos(v)
            sinv=sin(v)
            u=1.0_dblprec-cosv
            mrod(1,i,k)=hx*hx*u*emom(1,i,k)+cosv*emom(1,i,k)+hx*hy*u*emom(2,i,k)-hz*sinv*emom(2,i,k)+ &
               hx*hz*u*emom(3,i,k)+hy*sinv*emom(3,i,k)
            mrod(2,i,k)=hy*hx*u*emom(1,i,k)+hz*sinv*emom(1,i,k)+hy*hy*u*emom(2,i,k)+cosv*emom(2,i,k)+ &
               hy*hz*u*emom(3,i,k)-hx*sinv*emom(3,i,k)
            mrod(3,i,k)=hx*hz*u*emom(1,i,k)-hy*sinv*emom(1,i,k)+hz*hy*u*emom(2,i,k)+hx*sinv*emom(2,i,k)+ &
               hz*hz*u*emom(3,i,k)+cosv*emom(3,i,k)

            ! copy m(t) to emom2 and m(t+dt) to emom for heisge, save b(t)
            emom2(:,i,k)=emom(:,i,k)
            emomM(:,i,k)=mrod(:,i,k)*mmom(i,k)

            emom(:,i,k)=mrod(:,i,k)

            b2eff(:,i,k)=bdup(:,i,k)

         end do
      end do
      !!$omp end parallel do

   contains

      function diff_first(atom,ensemble) result (d)
         use Constants, only : gama      
         implicit none
         integer, intent(in) :: atom, ensemble
         real(dblprec), dimension(3) :: d

         real(dblprec), parameter :: A = 3.0_dblprec / 2.0_dblprec
         real(dblprec), parameter :: B = -2.0_dblprec
         real(dblprec), parameter :: C = 1.0_dblprec / 2.0_dblprec

         real(dblprec),dimension(3) :: numerator
         ! current previous and second previous
         real(dblprec),dimension(3) :: emom_0, emom_1, emom_2 
         integer :: current, prev, prev2
         real(dblprec) :: dt

         dt = delta_t * gama

         current = multiscaleBackbufferHead
         prev = multiscaleBackbufferHead-1
         prev2 = multiscaleBackbufferHead-2

         if (prev2 < 1) then
            prev2 = prev2 + ubound(multiscaleBackbuffer,4)
            if (prev < 1) then
               prev = prev + ubound(multiscaleBackbuffer,4)
            end if      
         end if

         emom_0 = multiscaleBackbuffer(:,atom,ensemble,current)
         emom_1 = multiscaleBackbuffer(:,atom,ensemble,prev)
         emom_2 = multiscaleBackbuffer(:,atom,ensemble,prev2)

         numerator = A*emom_0 + B*emom_1 + C*emom_2
         d = numerator / dt
         return

      end function diff_first
   end subroutine depondt_evolve_first_ms

   !-----------------------------------------------------------------------------
   !  SUBROUTINE: depondt_evolve_second
   !> @brief
   !> Second step of Depond solver, calculates the corrected effective field from
   !> the predicted effective fields. Rotates the moments in the corrected field
   !
   !> @author Anders Bergman
   !-----------------------------------------------------------------------------
   subroutine depondt_evolve_second_ms(Natom,Nred,Mensemble,lambda1_array,beff,b2eff,  &
         btorque, emom, emom2, delta_t, stt,do_she,she_btorque,do_sot,sot_btorque,  &
         red_atom_list,dband)

      use Constants, only : gama
      !
      implicit none
      !
      integer, intent(in) :: Nred   !< Number of atoms that evolve
      integer, intent(in) :: Natom  !< Number of atoms in system
      integer, intent(in) :: Mensemble !< Number of ensembles
      real(dblprec), intent(in) :: delta_t !< Time step
      character(len=1), intent(in) :: STT    !< Treat spin transfer torque?
      character(len=1), intent(in) :: do_she !< Treat the SHE spin transfer torque
      character(len=1), intent(in) :: do_sot !< Treat the general SOT model
      integer, dimension(Nred), intent(in) :: red_atom_list !< List of indices of atoms that evolve
      real(dblprec), dimension(Natom), intent(in) :: lambda1_array !< Damping parameter
      real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: beff !< Total effective field from application of Hamiltonian
      real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: b2eff !< Temporary storage of magnetic field
      real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: btorque !< Spin transfer torque
      real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: she_btorque !< SHE spin transfer torque
      real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: sot_btorque !< Spin orbit torque
      real(dblprec), dimension(3,Natom,Mensemble), intent(out) :: emom   !< Current unit moment vector
      real(dblprec), dimension(3,Natom,Mensemble), intent(out) :: emom2  !< Final (or temporary) unit moment vector
      type(DampingBandData), intent(in) :: dband !< Multiscale damping band parameters
      integer :: i,k,ired
      real(dblprec) :: v,Bnorm,hx,hy,hz
      real(dblprec) :: u,cosv,sinv,lldamp, she_fac, stt_fac,sot_fac

      if(stt/='N') then
         stt_fac=1.0_dblprec
      else
         stt_fac=0.0_dblprec
      end if

      if(do_she/='N') then
         she_fac=1.0_dblprec
      else
         she_fac=0.0_dblprec
      end if
      if(do_sot/='N') then
         sot_fac=1.0_dblprec
      else
         sot_fac=0.0_dblprec
      end if

      bdup(:,:,:)=0.0_dblprec
      if(stt=='Y') then
         !!$omp parallel do default(shared) private(i,ired,k)  schedule(static) collapse(2)
         do k=1,Mensemble
            do i=1,Natom
               !i=red_atom_list(ired)
               ! Adding STT and SHE torques if present (prefactor instead of if-statement)
               bdup(:,i,k)=bdup(:,i,k)+stt_fac*btorque(:,i,k)+she_fac*she_btorque(:,i,k)
            end do
         end do
         !!$omp end parallel do
      endif
      if (do_sot=='Y') then
         !!$omp parallel do default(shared) private(i,ired,k)  schedule(static) collapse(2)
         do k=1,Mensemble
            do i=1,Natom
               !i=red_atom_list(ired)
               ! Adding STT and SHE torques if present (prefactor instead of if-statement)
               bdup(:,i,k)=bdup(:,i,k)+sot_fac*sot_btorque(:,i,k)
            end do
         end do
         !!$omp end parallel do
      endif
      if (do_she=='Y') then
         !!$omp parallel do default(shared) private(i,ired,k)  schedule(static) collapse(2)
         do k=1,Mensemble
            do i=1,Natom
               !i=red_atom_list(ired)
               ! Adding STT and SHE torques if present (prefactor instead of if-statement)
               bdup(:,i,k)=bdup(:,i,k)+she_fac*she_btorque(:,i,k)
            end do
         end do
         !!$omp end parallel do
      endif

      !!$omp parallel do default(shared) private(ired,i,k,Bnorm,hx,hy,hz,v,lldamp,cosv,sinv,u) schedule(static) collapse(2)
      do k=1,Mensemble
         do i=1,Natom

            !i=red_atom_list(ired)
            ! Construct local field
            bloc(:,i,k)=beff(:,i,k)+btherm(:,i,k)

            ! Construct effective field (including damping term)
            bdup(1,i,k)=bdup(1,i,k)+bloc(1,i,k)+lambda1_array(i)*emom(2,i,k)*bloc(3,i,k)-lambda1_array(i)*emom(3,i,k)*bloc(2,i,k)
            bdup(2,i,k)=bdup(2,i,k)+bloc(2,i,k)+lambda1_array(i)*emom(3,i,k)*bloc(1,i,k)-lambda1_array(i)*emom(1,i,k)*bloc(3,i,k)
            bdup(3,i,k)=bdup(3,i,k)+bloc(3,i,k)+lambda1_array(i)*emom(1,i,k)*bloc(2,i,k)-lambda1_array(i)*emom(2,i,k)*bloc(1,i,k)

            dedt(:,i,k)=diff_second(i,k)
         end do
      end do

      call buildbeff_ms(Natom, Mensemble,lambda1_array,emom, btorque,stt,do_she,&
         she_btorque,dband,dedt)
         !she_btorque,dband,diff_second)

      !!$omp parallel do default(shared) private(ired,i,k,Bnorm,hx,hy,hz,v,lldamp,cosv,sinv,u) schedule(static) collapse(2)
      do k=1,Mensemble
         do i=1,Natom


            ! Corrected field
            ! Euler
            !bdup(:,i,k)=0.00_dblprec*bdup(:,i,k)+1.00_dblprec*b2eff(:,i,k)
            ! Heun
            bdup(:,i,k)=0.50_dblprec*bdup(:,i,k)+0.50_dblprec*b2eff(:,i,k)
            ! Ralston
            !bdup(:,i,k)=0.75_dblprec*bdup(:,i,k)+0.25_dblprec*b2eff(:,i,k)
            ! Midpoint
            !bdup(:,i,k)=1.00_dblprec*bdup(:,i,k)+0.00_dblprec*b2eff(:,i,k)
            !
            emom(:,i,k)=emom2(:,i,k)

            ! Set up rotation matrices and perform rotations
            lldamp=1.0_dblprec/(1.0_dblprec+lambda1_array(i)**2)
            Bnorm=bdup(1,i,k)**2+bdup(2,i,k)**2+bdup(3,i,k)**2
            Bnorm=sqrt(Bnorm)+1.0d-15
            hx=bdup(1,i,k)/Bnorm
            hy=bdup(2,i,k)/Bnorm
            hz=bdup(3,i,k)/Bnorm
            v=Bnorm*delta_t*gama*lldamp
            cosv=cos(v)
            sinv=sin(v)
            u=1.0_dblprec-cosv
            mrod(1,i,k)=hx*hx*u*emom(1,i,k)+cosv*emom(1,i,k)+hx*hy*u*emom(2,i,k)-hz*sinv*emom(2,i,k)+ &
               hx*hz*u*emom(3,i,k)+hy*sinv*emom(3,i,k)
            mrod(2,i,k)=hy*hx*u*emom(1,i,k)+hz*sinv*emom(1,i,k)+hy*hy*u*emom(2,i,k)+cosv*emom(2,i,k)+ &
               hy*hz*u*emom(3,i,k)-hx*sinv*emom(3,i,k)
            mrod(3,i,k)=hx*hz*u*emom(1,i,k)-hy*sinv*emom(1,i,k)+hz*hy*u*emom(2,i,k)+hx*sinv*emom(2,i,k)+ &
               hz*hz*u*emom(3,i,k)+cosv*emom(3,i,k)

            ! Final update
            emom2(:,i,k)=mrod(:,i,k)
         end do
      end do
      !!$omp end parallel do


   contains

      function diff_second(atom,ensemble) result (d)
         use Constants, only : gama
         implicit none
         integer, intent(in) :: atom, ensemble
         real(dblprec), dimension(3) :: d

         real(dblprec), parameter :: A = 3.0_dblprec / 2.0_dblprec
         real(dblprec), parameter :: B = -2.0_dblprec
         real(dblprec), parameter :: C = 1.0_dblprec / 2.0_dblprec

         real(dblprec),dimension(3) :: numerator
         ! current previous and second previous
         real(dblprec),dimension(3) :: emom_0, emom_1, emom_2 
         integer :: current, prev
         real(dblprec) :: dt

         dt = delta_t * gama

         current = multiscaleBackbufferHead
         prev = multiscaleBackbufferHead-1

         if (prev < 1) then
            prev = prev + ubound(multiscaleBackbuffer,4)
         end if

         emom_0 = emom(:,atom,ensemble)
         emom_1 = multiscaleBackbuffer(:,atom,ensemble,current)
         emom_2 = multiscaleBackbuffer(:,atom,ensemble,prev)

         numerator = A*emom_0 + B*emom_1 + C*emom_2
         d = numerator / dt
         return      
      end function diff_second
   end subroutine depondt_evolve_second_ms


   subroutine buildbeff_ms(Natom, Mensemble,lambda1_array,emom, btorque, stt,do_she,&
         she_btorque,dband,dedt)

      implicit none
      !
      integer, intent(in) :: Natom !< Number of atoms in system
      integer, intent(in) :: Mensemble !< Number of ensembles 
      real(dblprec), dimension(Natom), intent(in) :: lambda1_array !< Damping parameter
      real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: emom   !< Current unit moment vector
      real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: btorque !< Spin transfer torque
      real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: she_btorque !< SHE spin transfer torque
      character(len=1), intent(in) :: STT !< Treat spin transfer torque? 
      character(len=1), intent(in) :: do_she !< Treat SHE spin transfer torque
      type(DampingBandData), intent(in) :: dband !< Multiscale damping band parameters
      real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: dedt  !<  

      !
      integer :: i, k, index
      real(dblprec) :: ma1,ma2,ma3, gamma, sdnorm

      !
      if(dband%enable) then
         !!$omp parallel do default(shared) private(i,k) schedule(dynamic)
         do i=1,Natom
            index = dband%interpolation%indices(i)
            do k=1,Mensemble
               if(index .ne. 0) then
                  sdnorm = sum(dedt(:,i,k)**2) ** 0.25_dblprec         
                  gamma = dband%coefficients(index) * sdnorm
                  ma1 = (1-lambda1_array(i)**2)*gamma*dband%preinterpolation(1,index,k)
                  ma2 = (1-lambda1_array(i)**2)*gamma*dband%preinterpolation(2,index,k) 
                  ma3 = (1-lambda1_array(i)**2)*gamma*dband%preinterpolation(3,index,k)

                  bdup(1,i,k) = bdup(1,i,k) + (emom(2,i,k)*ma3 - emom(3,i,k)*ma2)
                  bdup(2,i,k) = bdup(2,i,k) + (emom(3,i,k)*ma1 - emom(1,i,k)*ma3)
                  bdup(3,i,k) = bdup(3,i,k) + (emom(1,i,k)*ma2 - emom(2,i,k)*ma1)
               end if
            end do
         end do
         !!$omp end parallel do
      end if

   end subroutine buildbeff_ms

end module Depondt_ms
