!-------------------------------------------------------------------------------
!> MODULE: Heun_proper
!> @brief
!> A predictor-corrector Heun solver for the LLG-equation
!> @author
!! Johan Hellsvik
!> @copyright
!! Copyright (C) 2008-2018 UppASD group
!! This file is distributed under the terms of the
!! GNU General Public License. 
!! See http://www.gnu.org/copyleft/gpl.txt
!-------------------------------------------------------------------------------
module Heun_proper
   use Parameters
   use Profiling
   !
   implicit none

   private

   public :: evolve4p, evolve4f

   !
contains


   !-----------------------------------------------------------------------------
   !> SUBROUTINE: evolve4p
   !> @brief
   !> First part of the integrator
   !-----------------------------------------------------------------------------
   subroutine evolve4p(Natom, Mensemble, Landeg, llg, lambda1_array, compensate_drift, &
         beff, b2eff, mmomi, emom2, emom, emomM, mmom, deltat,Temp_array,temprescale,thermal_field)
      !
      use Constants
      use RandomNumbers, only : ranv, rng_gaussian
      !
      implicit none

      integer, intent(in) :: Natom !< Number of atoms in system
      integer, intent(in) :: Mensemble !< Number of ensembles
      real(dblprec), dimension(Natom), intent(in) :: Landeg  !< Gyromagnetic ratio
      integer, intent(in) :: llg !< Type of equation of motion (1=LLG)
      integer, intent(in) :: compensate_drift !< Correct for drift in RNG
      real(dblprec), dimension(Natom), intent(in) :: lambda1_array !< Damping parameter
      real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: beff !< Total effective field from application of Hamiltonian
      real(dblprec), dimension(3,Natom,Mensemble), intent(out) :: b2eff !< Temporary storage of magnetic field
      real(dblprec), dimension(3,Natom,Mensemble), intent(inout) :: thermal_field
      real(dblprec), dimension(Natom,Mensemble), intent(in) :: mmomi !< Inverse of magnitude of magnetic moments
      real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: emom2  !< Final (or temporary) unit moment vector
      real(dblprec), dimension(3,Natom,Mensemble), intent(out) :: emom  !< Current unit moment vector
      real(dblprec), dimension(3,Natom,Mensemble), intent(out) :: emomM !< Current magnetic moment vector
      real(dblprec), dimension(Natom,Mensemble), intent(in) :: mmom !< Magnitude of magnetic moments
      real(dblprec), intent(in) :: deltat !< Time step
      real(dblprec), dimension(Natom), intent(in) :: Temp_array  !< Temperature
      real(dblprec), intent(in) :: temprescale  !< Temperature rescaling from QHB

      real(dblprec) :: D
      real(dblprec) :: dt, mu, sigma
      real(dblprec) :: rx,ry,rz
      real(dblprec), dimension(Natom) :: Dk

      integer :: i,k

      !JohanH-Dec23: Temporarily evolve3 and heun3
      !work only for Mensemble=1
      !   LLG equations ONE universal damping
      Dk=0.0d0
      if (llg==1) then
         Dk=lambda1_array/(1+lambda1_array**2)*k_bolt/mub
         dt=deltat*gama !rescaled time
      else
         stop 'Only llg=1 is currently supported'
      endif
      !
      call rng_gaussian(ranv,3*Natom*Mensemble,1.d0)
      mu=0d0

      !$omp parallel do default(shared) private(k,i,D,sigma) collapse(2)
      do k=1, Mensemble
         do i=1, Natom
            D=Dk(i)*mmomi(i,k)*Temp_array(i)*temprescale
            sigma=sqrt(dt*2*D)
            ranv(1,i,k)=ranv(1,i,k)*sigma
            ranv(2,i,k)=ranv(2,i,k)*sigma
            ranv(3,i,k)=ranv(3,i,k)*sigma
         end do
      end do
      !$omp end parallel do

      if(compensate_drift==1) then
         do k=1, Mensemble
            rx=0.0d0;ry=0.0d0;rz=0.0d0
            do i=1, Natom
               rx=rx+ranv(1,i,k)
               ry=ry+ranv(2,i,k)
               rz=rz+ranv(3,i,k)
            end do
            rx=rx/Natom
            ry=ry/Natom
            rz=rz/Natom
            do i=1, Natom
               ranv(1,i,k)=ranv(1,i,k)-rx
               ranv(2,i,k)=ranv(2,i,k)-ry
               ranv(3,i,k)=ranv(3,i,k)-rz
            end do
         end do
      end if
      call heun4p(Natom, Mensemble, Landeg, lambda1_array,  beff,  b2eff, emom, emomM, emom2, mmom, dt)

      thermal_field=ranv

   end subroutine evolve4p


   !-----------------------------------------------------------------------------
   !> SUBROUTINE: evolve4f
   !> @brief
   !> Second part of the integrator
   !-----------------------------------------------------------------------------
   subroutine evolve4f(Natom, Mensemble, Landeg, lambda1_array, llg, beff,  b2eff, emom2, emom, deltat)
      !
      use Constants
      !
      implicit none

      integer, intent(in) :: Natom !< Number of atoms in system
      integer, intent(in) :: Mensemble !< Number of independent simulations
      real(dblprec), dimension(Natom), intent(in) :: Landeg  !< Gyromagnetic ratio
      integer, intent(in) :: llg !< Type of equation of motion (1=LLG)
      real(dblprec), dimension(Natom), intent(in) :: lambda1_array !< Damping parameter
      real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: beff !< Total effective field from application of Hamiltonian
      real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: b2eff !< Temporary storage of magnetic field
      real(dblprec), dimension(3,Natom,Mensemble), intent(out) :: emom2  !< Final (or temporary) unit moment vector
      real(dblprec), dimension(3,Natom,Mensemble), intent(out) :: emom   !< Current unit moment vector
      real(dblprec), intent(in) :: deltat !< Time step

      real(dblprec) :: etot, dt

      integer :: i,k

      if (llg==1) then
         dt=deltat*gama !rescaled time
      else
         stop 'Only llg=1 is currently supported'
      endif

      call heun4f(Natom, Mensemble, Landeg, lambda1_array, beff, b2eff, emom, emom2, dt)

      !$omp parallel do schedule(static) default(shared) private(i,k,etot) collapse(2)
      do i=1,Natom
         do k=1,Mensemble
            etot=1/sqrt(emom2(1,i,k)**2+emom2(2,i,k)**2+emom2(3,i,k)**2)
            emom2(1,i,k)=emom2(1,i,k)*etot
            emom2(2,i,k)=emom2(2,i,k)*etot
            emom2(3,i,k)=emom2(3,i,k)*etot
         end do
      end do
      !$omp end parallel do

   end subroutine evolve4f


   !-----------------------------------------------------------------------------
   !> SUBROUTINE: heun4p
   !> @brief
   !> Performs predictor step
   !-----------------------------------------------------------------------------
   subroutine heun4p(Natom, Mensemble, Landeg, lambda1_array,  beff,  b2eff, emom, emomM, emom2, mmom, dt)
      !
      use RandomNumbers, only : ranv
      !
      implicit none

      integer, intent(in) :: Natom !< Number of atoms in system
      integer, intent(in) :: Mensemble !< Number of ensembles
      real(dblprec), dimension(Natom), intent(in) :: Landeg  !< Gyromagnetic ratio
      real(dblprec), dimension(Natom), intent(in) :: lambda1_array !< Damping parameter
      real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: beff !< Total effective field from application of Hamiltonian
      real(dblprec), dimension(3,Natom,Mensemble), intent(out) :: b2eff !< Temporary storage of magnetic field
      real(dblprec), dimension(3,Natom,Mensemble), intent(out) :: emom   !< Current unit moment vector
      real(dblprec), dimension(3,Natom,Mensemble), intent(out) :: emomM  !< Current magnetic moment vector
      real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: emom2  !< Final (or temporary) unit moment vector
      real(dblprec), dimension(Natom,Mensemble), intent(in) :: mmom !< Magnitude of magnetic moments
      real(dblprec), intent(in) :: dt !< Time step

      real(dblprec) :: b1xx, b1xy, b1xz
      real(dblprec) :: b1yx, b1yy, b1yz
      real(dblprec) :: b1zx, b1zy, b1zz
      real(dblprec) :: a1x, a1y, a1z
      real(dblprec) :: etx, ety, etz
      real(dblprec) :: dtg
      real(dblprec) :: dwx, dwy, dwz
      real(dblprec) :: e1x, e1y, e1z
      real(dblprec) :: hlpx,hlpy,hlpz
      real(dblprec) :: hesx,hesy,hesz,prot
      !
      ! ... Local variables ...
      integer :: i,k
      !
      !$omp parallel do default(shared) &

      !$omp private(i,e1x,e1y,e1z,dwx,dwy,dwz,hesx,hesy,hesz &
      !$omp ,hlpx,hlpy,hlpz,etx,ety,etz,a1x,a1y,a1z,prot,b1xx &
      !$omp ,b1xy,b1xz,b1yx,b1yy,b1yz,b1zx,b1zy,b1zz,dtg) collapse(2)
      do i=1,Natom
         do k=1,Mensemble

            !Store away the initial effective field that is used
            !for the support value
            b2eff(1,i,k)=beff(1,i,k)
            b2eff(2,i,k)=beff(2,i,k)
            b2eff(3,i,k)=beff(3,i,k)
            e1x=emom2(1,i,k)
            e1y=emom2(2,i,k)
            e1z=emom2(3,i,k)
            dtg=dt*Landeg(i)
            dwx=ranv(1,i,k)*Landeg(i)
            dwy=ranv(2,i,k)*Landeg(i)
            dwz=ranv(3,i,k)*Landeg(i)
            hesx=beff(1,i,k)
            hlpx=lambda1_array(i)*beff(1,i,k)
            hesy=beff(2,i,k)
            hlpy=lambda1_array(i)*beff(2,i,k)
            hesz=beff(3,i,k)
            hlpz=lambda1_array(i)*beff(3,i,k)

            ! Calculate a1 and b1
            prot=e1x*hlpx+e1y*hlpy+e1z*hlpz
            a1x=hlpx+e1z*hesy-e1y*hesz-e1x*prot
            a1y=hlpy-e1z*hesx+e1x*hesz-e1y*prot
            a1z=hlpz+e1y*hesx-e1x*hesy-e1z*prot
            b1xx=          lambda1_array(i)-e1x*e1x*lambda1_array(i)
            b1xy=e1z-e1x*e1y*lambda1_array(i)
            b1xz=-e1y-e1x*e1z*lambda1_array(i)
            b1yx=-e1z-e1y*e1x*lambda1_array(i)
            b1yy=          lambda1_array(i)-e1y*e1y*lambda1_array(i)
            b1yz=e1x-e1y*e1z*lambda1_array(i)
            b1zx=e1y-e1z*e1x*lambda1_array(i)
            b1zy=-e1x-e1z*e1y*lambda1_array(i)
            b1zz=          lambda1_array(i)-e1z*e1z*lambda1_array(i)

            etx=e1x+a1x*dtg+b1xx*dwx+b1xy*dwy+b1xz*dwz
            ety=e1y+a1y*dtg+b1yx*dwx+b1yy*dwy+b1yz*dwz
            etz=e1z+a1z*dtg+b1zx*dwx+b1zy*dwy+b1zz*dwz

            !Store the euler-type support values in emomM
            !so that the effective field can be recalculated
            !in heisge.
            emom(1,i,k)=etx
            emom(2,i,k)=ety
            emom(3,i,k)=etz
            emomM(1,i,k)=emom(1,i,k)*mmom(i,k)
            emomM(2,i,k)=emom(2,i,k)*mmom(i,k)
            emomM(3,i,k)=emom(3,i,k)*mmom(i,k)
         end do
      end do
      !$omp end parallel do

      return

   end subroutine heun4p

   !-----------------------------------------------------------------------------
   !> SUBROUTINE: heun4f
   !> @brief
   !> Performs corrector step
   !-----------------------------------------------------------------------------
   subroutine heun4f(Natom, Mensemble, Landeg, lambda1_array, beff, b2eff, emom, emom2, dt)

      use RandomNumbers, only : ranv
      !

      implicit none

      integer, intent(in) :: Natom !< Number of atoms in system
      integer, intent(in) :: Mensemble !< Number of ensembles
      real(dblprec), dimension(Natom), intent(in) :: Landeg  !< Gyromagnetic ratio
      real(dblprec), dimension(Natom), intent(in) :: lambda1_array !< Damping parameter
      real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: beff !< Total effective field from application of Hamiltonian
      real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: b2eff !< Temporary storage of magnetic field
      real(dblprec), dimension(3,Natom,Mensemble), intent(inout) :: emom   !< Current unit moment vector
      real(dblprec), dimension(3,Natom,Mensemble), intent(inout) :: emom2  !< Final (or temporary) unit moment vector
      real(dblprec), intent(in) :: dt !< Time step

      ! ... Local variables ...
      integer :: i,k

      real(dblprec) :: dtg
      real(dblprec) :: dwx, dwy, dwz

      real(dblprec) :: bxx, bxy, bxz
      real(dblprec) :: byx, byy, byz
      real(dblprec) :: bzx, bzy, bzz
      real(dblprec) :: ax, ay, az

      real(dblprec) :: b1xx, b1xy, b1xz
      real(dblprec) :: b1yx, b1yy, b1yz
      real(dblprec) :: b1zx, b1zy, b1zz
      real(dblprec) :: a1x, a1y, a1z

      real(dblprec) :: b2xx, b2xy, b2xz
      real(dblprec) :: b2yx, b2yy, b2yz
      real(dblprec) :: b2zx, b2zy, b2zz
      real(dblprec) :: a2x, a2y, a2z

      real(dblprec) :: h1px,h1py,h1pz
      real(dblprec) :: h1esx,h1esy,h1esz
      real(dblprec) :: h2px,h2py,h2pz
      real(dblprec) :: h2esx,h2esy,h2esz

      real(dblprec) :: etx, ety, etz
      real(dblprec) :: e1x, e1y, e1z
      real(dblprec) :: e2x, e2y, e2z

      real(dblprec) :: prot

      !Executable statements

      !$omp parallel do default(shared) &
      !$omp private(i,dtg, dwx,dwy,dwz, &
      !$omp bxx, bxy, bxz, byx, byy, byz, bzx, bzy, bzz, ax, ay, az, &
      !$omp b1xx, b1xy, b1xz, b1yx, b1yy, b1yz, b1zx, b1zy, b1zz, a1x, a1y, a1z, &
      !$omp b2xx, b2xy, b2xz, b2yx, b2yy, b2yz, b2zx, b2zy, b2zz, a2x, a2y, a2z, &
      !$omp h1px, h1py, h1pz, h1esx, h1esy, h1esz, &
      !$omp h2px, h2py, h2pz, h2esx, h2esy, h2esz, &
      !$omp etx, ety, etz, e1x, e1y, e1z, e2x, e2y, e2z ,prot) collapse(2)
      do i=1,Natom
         do k=1,Mensemble
            dtg=dt*Landeg(i)
            dwx=ranv(1,i,k)*Landeg(i)
            dwy=ranv(2,i,k)*Landeg(i)
            dwz=ranv(3,i,k)*Landeg(i)

            e1x=emom(1,i,k)
            e1y=emom(2,i,k)
            e1z=emom(3,i,k)
            e2x=emom2(1,i,k)
            e2y=emom2(2,i,k)
            e2z=emom2(3,i,k)

            h1esx=beff(1,i,k)
            h1px=lambda1_array(i)*beff(1,i,k)
            h1esy=beff(2,i,k)
            h1py=lambda1_array(i)*beff(2,i,k)
            h1esz=beff(3,i,k)
            h1pz=lambda1_array(i)*beff(3,i,k)

            h2esx=b2eff(1,i,k)
            h2px=lambda1_array(i)*b2eff(1,i,k)
            h2esy=b2eff(2,i,k)
            h2py=lambda1_array(i)*b2eff(2,i,k)
            h2esz=b2eff(3,i,k)
            h2pz=lambda1_array(i)*b2eff(3,i,k)

            !Calculate Ai_ytilde and Bik_ytilde
            prot=e1x*h1px+e1y*h1py+e1z*h1pz
            a1x=h1px+e1z*h1esy-e1y*h1esz-e1x*prot
            a1y=h1py-e1z*h1esx+e1x*h1esz-e1y*prot
            a1z=h1pz+e1y*h1esx-e1x*h1esy-e1z*prot
            b1xx=          lambda1_array(i)-e1x*e1x*lambda1_array(i)
            b1xy=e1z-e1x*e1y*lambda1_array(i)
            b1xz=-e1y-e1x*e1z*lambda1_array(i)
            b1yx=-e1z-e1y*e1x*lambda1_array(i)
            b1yy=          lambda1_array(i)-e1y*e1y*lambda1_array(i)
            b1yz=e1x-e1y*e1z*lambda1_array(i)
            b1zx=e1y-e1z*e1x*lambda1_array(i)
            b1zy=-e1x-e1z*e1y*lambda1_array(i)
            b1zz=          lambda1_array(i)-e1z*e1z-lambda1_array(i)

            !Calculate Ai and Bik
            prot=e2x*h2px+e2y*h2py+e2z*h2pz
            a2x=h2px+e2z*h2esy-e2y*h2esz-e2x*prot
            a2y=h2py-e2z*h2esx+e2x*h2esz-e2y*prot
            a2z=h2pz+e2y*h2esx-e2x*h2esy-e2z*prot
            b2xx=          lambda1_array(i)-e2x*e2x*lambda1_array(i)
            b2xy= e2z-e2x*e2y*lambda1_array(i)
            b2xz=-e2y-e2x*e2z*lambda1_array(i)
            b2yx=-e2z-e2y*e2x*lambda1_array(i)
            b2yy=          lambda1_array(i)-e2y*e2y*lambda1_array(i)
            b2yz= e2x-e2y*e2z*lambda1_array(i)
            b2zx= e2y-e2z*e2x*lambda1_array(i)
            b2zy=-e2x-e2z*e2y*lambda1_array(i)
            b2zz=          lambda1_array(i)-e2z*e2z*lambda1_array(i)

            !Calculate Aitot_tot=Ai_ytilde+Ai and
            !Bik_tot=Bik_ytilde+Bik
            ax=a1x+a2x
            ay=a1y+a2y
            az=a1z+a2z
            bxx=b1xx+b2xx
            bxy=b1xy+b2xy
            bxz=b1xz+b2xz
            byx=b1yx+b2yx
            byy=b1yy+b2yy
            byz=b1yz+b2yz
            bzx=b1zx+b2zx
            bzy=b1zy+b2zy
            bzz=b1zz+b2zz

            etx=ax*dtg+bxx*dwx+bxy*dwy+bxz*dwz
            ety=ay*dtg+byx*dwx+byy*dwy+byz*dwz
            etz=az*dtg+bzx*dwx+bzy*dwy+bzz*dwz
            emom2(1,i,k)=e2x+0.5d0*etx
            emom2(2,i,k)=e2y+0.5d0*ety
            emom2(3,i,k)=e2z+0.5d0*etz
         end do
      end do
      !$omp end parallel do

      return

   end subroutine heun4f


end module heun_proper
