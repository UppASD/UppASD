!------------------------------------------------------------------------------------
!  MODULE: Heun_proper
!> @brief
!> A predictor-corrector Heun solver for the LLG-equation
!> @copyright
!> GNU Public License.
!------------------------------------------------------------------------------------
module Heun_proper
   use Parameters
   use Profiling
   !
   implicit none

   real(dblprec), dimension(:,:,:), allocatable :: btorque_full !< Resulting effective field

   private

   public :: evolve4p, evolve4f, allocate_aux_heun_fields

   !
contains

   !---------------------------------------------------------------------------------
   !  SUBROUTINE: evolve4p
   !> @brief
   !> First part of the integrator
   !---------------------------------------------------------------------------------
   subroutine evolve4p(Natom,Mensemble,Landeg,llg,lambda1_array,compensate_drift,   &
      beff,b2eff,mmomi,emom2,emom,emomM,mmom,deltat,Temp_array,temprescale,         &
      thermal_field,STT,do_she,do_sot,btorque,she_btorque,sot_btorque,Nred,         &
      red_atom_list)

      use Constants
      use RandomNumbers, only : ranv, rng_gaussian
      !
      implicit none

      integer, intent(in) :: llg    !< Type of equation of motion (1=LLG)
      integer, intent(in) :: Nred   !< Number of moments that can be updated
      integer, intent(in) :: Natom  !< Number of atoms in system
      integer, intent(in) :: Mensemble !< Number of ensembles
      integer, intent(in) :: compensate_drift !< Correct for drift in RNG
      real(dblprec), intent(in) :: deltat !< Time step
      real(dblprec), intent(in) :: temprescale  !< Temperature rescaling from QHB
      character(len=1), intent(in) :: STT    !< Treat spin transfer torque?
      character(len=1), intent(in) :: do_she !< Treat the SHE spin transfer torque
      character(len=1), intent(in) :: do_sot !< Treat the general SOT model
      integer, dimension(Nred), intent(in) :: red_atom_list !< Reduced list containing atoms allowed to evolve in a fixed moment calculation
      real(dblprec), dimension(Natom), intent(in) :: Landeg !< Gyromagnetic ratio
      real(dblprec), dimension(Natom), intent(in) :: lambda1_array !< Damping parameter
      real(dblprec), dimension(Natom), intent(in) :: Temp_array  !< Temperature
      real(dblprec), dimension(Natom,Mensemble), intent(in) :: mmom !< Magnitude of magnetic moments
      real(dblprec), dimension(Natom,Mensemble), intent(in) :: mmomi !< Inverse of magnitude of magnetic moments
      real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: beff !< Total effective field from application of Hamiltonian
      real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: emom2  !< Final (or temporary) unit moment vector
      real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: btorque      !< Spin transfer torque
      real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: she_btorque  !< SHE spin transfer torque
      real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: sot_btorque  !< Spin orbit torque
      !.. Output variables
      real(dblprec), dimension(3,Natom,Mensemble), intent(out) :: emom  !< Current unit moment vector
      real(dblprec), dimension(3,Natom,Mensemble), intent(out) :: emomM !< Current magnetic moment vector
      real(dblprec), dimension(3,Natom,Mensemble), intent(out) :: b2eff !< Temporary storage of magnetic field
      ! .. In/out variables
      real(dblprec), dimension(3,Natom,Mensemble), intent(inout) :: thermal_field   !< Thermal field

      real(dblprec) :: D
      real(dblprec) :: dt, mu, sigma
      real(dblprec) :: rx,ry,rz
      real(dblprec), dimension(Natom) :: Dk

      integer :: i,k,ired

      !JohanH-Dec23: Temporarily evolve3 and heun3
      !work only for Mensemble=1
      !   LLG equations ONE universal damping
      Dk=0.0_dblprec
      if (llg==1) then
         Dk=lambda1_array/(1+lambda1_array**2)*k_bolt/mub
         dt=deltat*gama !rescaled time
      else
         stop 'Only llg=1 is currently supported'
      endif
      !
      call rng_gaussian(ranv,3*Natom*Mensemble,1.0_dblprec)
      mu=0_dblprec

      !$omp parallel do default(shared) private(ired,k,i,D,sigma) collapse(2)
      do k=1, Mensemble
         do ired=1, Nred
            i=red_atom_list(ired)
            D=Dk(i)*mmomi(i,k)*Temp_array(i)*temprescale
            sigma=sqrt(dt*2*D)
            ranv(1,i,k)=ranv(1,i,k)*sigma
            ranv(2,i,k)=ranv(2,i,k)*sigma
            ranv(3,i,k)=ranv(3,i,k)*sigma
         end do
      end do
      !$omp end parallel do

      btorque_full=0.0_dblprec
      if(stt/='N') then
         !$omp parallel do default(shared) private(ired,i,k)  schedule(static) collapse(2)
         do k=1,Mensemble
            do ired=1,Nred
               i=red_atom_list(ired)
               ! Adding STT and SHE torques if present (prefactor instead of if-statement)
               btorque_full(:,i,k)=btorque_full(:,i,k)+btorque(:,i,k)
            end do
         end do
         !$omp end parallel do
      end if

      if(do_she/='N') then
         !$omp parallel do default(shared) private(ired,i,k)  schedule(static) collapse(2)
         do k=1,Mensemble
            do ired=1,Nred
               i=red_atom_list(ired)
               ! Adding STT and SHE torques if present (prefactor instead of if-statement)
               btorque_full(:,i,k)= btorque_full(:,i,k)+she_btorque(:,i,k)
            end do
         end do
         !$omp end parallel do
      end if
      if(do_sot/='N') then
         !$omp parallel do default(shared) private(ired,i,k)  schedule(static) collapse(2)
         do k=1,Mensemble
            do ired=1,Nred
               i=red_atom_list(ired)
               ! Adding STT, SHE and SOT torques if present (prefactor instead of if-statement)
               btorque_full(:,i,k)= btorque_full(:,i,k)+sot_btorque(:,i,k)
            end do
         end do
         !$omp end parallel do
      end if

      if(compensate_drift==1) then
         do k=1, Mensemble
            rx=0.0_dblprec;ry=0.0_dblprec;rz=0.0_dblprec
            do ired=1, Nred
               i=red_atom_list(ired)
               rx=rx+ranv(1,i,k)
               ry=ry+ranv(2,i,k)
               rz=rz+ranv(3,i,k)
            end do
            rx=rx/Natom
            ry=ry/Natom
            rz=rz/Natom
            do ired=1, Nred
               i=red_atom_list(ired)
               ranv(1,i,k)=ranv(1,i,k)-rx
               ranv(2,i,k)=ranv(2,i,k)-ry
               ranv(3,i,k)=ranv(3,i,k)-rz
            end do
         end do
      end if
      call heun4p(Natom,Mensemble,Landeg,lambda1_array,beff,b2eff,emom,emomM,emom2, &
         mmom,dt,Nred,red_atom_list)

      thermal_field=ranv

   end subroutine evolve4p

   !---------------------------------------------------------------------------------
   !  SUBROUTINE: evolve4f
   !> @brief
   !> Second part of the integrator
   !---------------------------------------------------------------------------------
   subroutine evolve4f(Natom,Mensemble,Landeg,lambda1_array,llg,beff,b2eff,emom2,   &
      emom,deltat,STT,do_she,do_sot,btorque,she_btorque,sot_btorque,Nred,red_atom_list)
      !
      use Constants
      !
      implicit none

      integer, intent(in) :: llg    !< Type of equation of motion (1=LLG)
      integer, intent(in) :: Nred   !< Number of moments that can be updated
      integer, intent(in) :: Natom  !< Number of atoms in system
      integer, intent(in) :: Mensemble !< Number of independent simulations
      real(dblprec), intent(in) :: deltat !< Time step
      character(len=1), intent(in) :: STT    !< Treat spin transfer torque?
      character(len=1), intent(in) :: do_she !< Treat the SHE spin transfer torque
      character(len=1), intent(in) :: do_sot !< Treat the general SOT model
      integer, dimension(Nred), intent(in) :: red_atom_list !< Reduced list containing atoms allowed to evolve in a fixed moment calculation
      real(dblprec), dimension(Natom), intent(in) :: Landeg !< Gyromagnetic ratio
      real(dblprec), dimension(Natom), intent(in) :: lambda1_array !< Damping parameter
      real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: beff !< Total effective field from application of Hamiltonian
      real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: b2eff !< Temporary storage of magnetic field
      real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: btorque      !< Spin transfer torque
      real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: she_btorque  !< SHE spin transfer torque
      real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: sot_btorque  !< Spin orbit torque
      ! .. Output variables
      real(dblprec), dimension(3,Natom,Mensemble), intent(out) :: emom2  !< Final (or temporary) unit moment vector
      real(dblprec), dimension(3,Natom,Mensemble), intent(out) :: emom   !< Current unit moment vector

      real(dblprec) :: etot, dt

      integer :: i,k,ired

      if (llg==1) then
         dt=deltat*gama !rescaled time
      else
         stop 'Only llg=1 is currently supported'
      endif

      btorque_full=0.0_dblprec
      if(stt/='N') then
         !$omp parallel do default(shared) private(ired,i,k)  schedule(static) collapse(2)
         do k=1,Mensemble
            do ired=1,Nred
               i=red_atom_list(ired)
               ! Adding STT and SHE torques if present (prefactor instead of if-statement)
               btorque_full(:,i,k)=btorque_full(:,i,k)+btorque(:,i,k)
            end do
         end do
         !$omp end parallel do
      end if

      if(do_she/='N') then
         !$omp parallel do default(shared) private(ired,i,k)  schedule(static) collapse(2)
         do k=1,Mensemble
            do ired=1,Nred
               i=red_atom_list(ired)
               ! Adding STT and SHE torques if present (prefactor instead of if-statement)
               btorque_full(:,i,k)= btorque_full(:,i,k)+she_btorque(:,i,k)
            end do
         end do
         !$omp end parallel do
      end if
      if(do_sot/='N') then
         !$omp parallel do default(shared) private(ired,i,k)  schedule(static) collapse(2)
         do k=1,Mensemble
            do ired=1,Nred
               i=red_atom_list(ired)
               ! Adding STT, SHE and SOT torques if present (prefactor instead of if-statement)
               btorque_full(:,i,k)= btorque_full(:,i,k)+sot_btorque(:,i,k)
            end do
         end do
         !$omp end parallel do
      end if
      call heun4f(Natom,Mensemble,Landeg,lambda1_array,beff,b2eff,emom,emom2,dt,    &
         Nred,red_atom_list)

      !$omp parallel do schedule(static) default(shared) private(i,k,etot) collapse(2)
      do i=1,Natom
         do k=1,Mensemble
            etot=1.0_dblprec/norm2(emom2(:,i,k))
            emom2(1,i,k)=emom2(1,i,k)*etot
            emom2(2,i,k)=emom2(2,i,k)*etot
            emom2(3,i,k)=emom2(3,i,k)*etot
         end do
      end do
      !$omp end parallel do

   end subroutine evolve4f

   !---------------------------------------------------------------------------------
   !  SUBROUTINE: heun4p
   !> @brief
   !> Performs predictor step
   !---------------------------------------------------------------------------------
   subroutine heun4p(Natom,Mensemble,Landeg,lambda1_array,beff,b2eff,emom,emomM,    &
      emom2,mmom,dt,Nred,red_atom_list)
      !
      use RandomNumbers, only : ranv
      !
      implicit none

      integer, intent(in) :: Nred      !< Number of moments that can be updated
      integer, intent(in) :: Natom     !< Number of atoms in system
      integer, intent(in) :: Mensemble !< Number of ensembles
      real(dblprec), intent(in) :: dt  !< Time step
      integer, dimension(Nred), intent(in) :: red_atom_list !< Reduced list containing atoms allowed to evolve in a fixed moment calculation
      real(dblprec), dimension(Natom), intent(in) :: Landeg !< Gyromagnetic ratio
      real(dblprec), dimension(Natom), intent(in) :: lambda1_array      !< Damping parameter
      real(dblprec), dimension(Natom,Mensemble), intent(in) :: mmom     !< Magnitude of magnetic moments
      real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: beff   !< Total effective field from application of Hamiltonian
      real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: emom2  !< Final (or temporary) unit moment vector
      ! .. Output variables
      real(dblprec), dimension(3,Natom,Mensemble), intent(out) :: b2eff !< Temporary storage of magnetic field
      real(dblprec), dimension(3,Natom,Mensemble), intent(out) :: emom  !< Current unit moment vector
      real(dblprec), dimension(3,Natom,Mensemble), intent(out) :: emomM !< Current magnetic moment vector

      real(dblprec) :: b1xx, b1xy, b1xz
      real(dblprec) :: bsttx, bstty, bsttz
      real(dblprec) :: b1yx, b1yy, b1yz
      real(dblprec) :: b1zx, b1zy, b1zz
      real(dblprec) :: a1x, a1y, a1z
      real(dblprec) :: etx, ety, etz
      real(dblprec) :: dtg
      real(dblprec) :: dwx, dwy, dwz
      real(dblprec) :: e1x, e1y, e1z
      real(dblprec) :: hlpx,hlpy,hlpz
      real(dblprec) :: hesx,hesy,hesz,prot
      real(dblprec) :: lldamp
      !
      ! ... Local variables ...
      integer :: i,k,ired
      !
      !$omp parallel do default(shared) &

      !$omp private(i,ired,e1x,e1y,e1z,dwx,dwy,dwz,hesx,hesy,hesz &
      !$omp ,hlpx,hlpy,hlpz,etx,ety,etz,a1x,a1y,a1z,prot,b1xx &
      !$omp ,b1xy,b1xz,b1yx,b1yy,b1yz,b1zx,b1zy,b1zz,dtg,lldamp &
      !$omp ,bsttx,bstty,bsttz) collapse(2)
      do ired=1,Nred
         do k=1,Mensemble
            i=red_atom_list(ired)
            !Store away the initial effective field that is used
            !for the support value
            b2eff(1,i,k)=beff(1,i,k)
            b2eff(2,i,k)=beff(2,i,k)
            b2eff(3,i,k)=beff(3,i,k)
            bsttx=btorque_full(1,i,k)
            bstty=btorque_full(2,i,k)
            bsttz=btorque_full(3,i,k)
            e1x=emom2(1,i,k)
            e1y=emom2(2,i,k)
            e1z=emom2(3,i,k)
            lldamp=1.0_dblprec/(1.0_dblprec+lambda1_array(i)**2)
            dtg=dt*Landeg(i)*lldamp
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
            ! a1_\mu=\alpha B_\mu + (m\times B)_\mu -m_\mu (\alpha m \cdot B)
            a1x=hlpx+e1z*(hesy+bstty)-e1y*(hesz+bsttz)-e1x*prot
            a1y=hlpy-e1z*(hesx+bsttx)+e1x*(hesz+bsttz)-e1y*prot
            a1z=hlpz+e1y*(hesx+bsttx)-e1x*(hesy+bstty)-e1z*prot
            ! b1 is a matrix with diagonal entries given by b_\mu\mu=(\alpha-\alpha m_\mu m_\mu)
            ! and off diagonal entries b_\mu\nu are given by
            ! b_\mu\nu= m_\mu\times m_\nu -\alpha m_\mu m_\nu
            b1xx=          lambda1_array(i)-e1x*e1x*lambda1_array(i)
            b1xy= e1z-e1x*e1y*lambda1_array(i)
            b1xz=-e1y-e1x*e1z*lambda1_array(i)
            b1yx=-e1z-e1y*e1x*lambda1_array(i)
            b1yy=          lambda1_array(i)-e1y*e1y*lambda1_array(i)
            b1yz= e1x-e1y*e1z*lambda1_array(i)
            b1zx= e1y-e1z*e1x*lambda1_array(i)
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

   !---------------------------------------------------------------------------------
   !  SUBROUTINE: heun4f
   !> @brief
   !> Performs corrector step
   !---------------------------------------------------------------------------------
   subroutine heun4f(Natom,Mensemble,Landeg,lambda1_array,beff,b2eff,emom,emom2,dt, &
      Nred,red_atom_list)

      use RandomNumbers, only : ranv

      implicit none

      integer, intent(in) :: Nred   !< Number of moments that can be updated
      integer, intent(in) :: Natom  !< Number of atoms in system
      integer, intent(in) :: Mensemble !< Number of ensembles
      real(dblprec), intent(in) :: dt  !< Time step
      integer, dimension(Nred), intent(in) :: red_atom_list    !< Reduced list containing atoms allowed to evolve in a fixed moment calculation
      real(dblprec), dimension(Natom), intent(in) :: Landeg    !< Gyromagnetic ratio
      real(dblprec), dimension(Natom), intent(in) :: lambda1_array !< Damping parameter
      real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: beff !< Total effective field from application of Hamiltonian
      real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: b2eff !< Temporary storage of magnetic field
      real(dblprec), dimension(3,Natom,Mensemble), intent(inout) :: emom   !< Current unit moment vector
      real(dblprec), dimension(3,Natom,Mensemble), intent(inout) :: emom2  !< Final (or temporary) unit moment vector

      ! ... Local variables ...
      integer :: i,k, ired

      real(dblprec) :: dtg
      real(dblprec) :: dwx, dwy, dwz

      real(dblprec) :: bxx, bxy, bxz
      real(dblprec) :: bsttx, bstty, bsttz
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
      real(dblprec) :: lldamp
      !Executable statements

      !$omp parallel do default(shared) &
      !$omp private(i,ired,dtg, dwx,dwy,dwz, &
      !$omp bxx, bxy, bxz, byx, byy, byz, bzx, bzy, bzz, ax, ay, az, &
      !$omp b1xx, b1xy, b1xz, b1yx, b1yy, b1yz, b1zx, b1zy, b1zz, a1x, a1y, a1z, &
      !$omp b2xx, b2xy, b2xz, b2yx, b2yy, b2yz, b2zx, b2zy, b2zz, a2x, a2y, a2z, &
      !$omp h1px, h1py, h1pz, h1esx, h1esy, h1esz, &
      !$omp h2px, h2py, h2pz, h2esx, h2esy, h2esz, &
      !$omp etx, ety, etz, e1x, e1y, e1z, e2x, e2y, e2z ,prot,lldamp,&
      !$omp bsttx,bstty,bsttz) collapse(2)
      do ired=1,Nred
         do k=1,Mensemble
            i=red_atom_list(ired)
            lldamp=1.0_dblprec/(1.0_dblprec+lambda1_array(i)**2)
            dtg=dt*Landeg(i)*lldamp
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

            bsttx=btorque_full(1,i,k)
            bstty=btorque_full(2,i,k)
            bsttz=btorque_full(3,i,k)

            !Calculate Ai_ytilde and Bik_ytilde
            prot=e1x*h1px+e1y*h1py+e1z*h1pz
            a1x=h1px+e1z*(h1esy+bstty)-e1y*(h1esz+bsttz)-e1x*prot
            a1y=h1py-e1z*(h1esx+bsttx)+e1x*(h1esz+bsttz)-e1y*prot
            a1z=h1pz+e1y*(h1esx+bsttx)-e1x*(h1esy+bstty)-e1z*prot
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
            emom2(1,i,k)=e2x+0.5_dblprec*etx
            emom2(2,i,k)=e2y+0.5_dblprec*ety
            emom2(3,i,k)=e2z+0.5_dblprec*etz
         end do
      end do
      !$omp end parallel do

      return

   end subroutine heun4f

   !---------------------------------------------------------------------------------
   ! SUBROUTINE: allocate_aux_heun_fields
   !> @brief Allocation of auxilary fields for the treatment of STT and SOT based torques
   !> @author Jonathan Chico
   !---------------------------------------------------------------------------------
   subroutine allocate_aux_heun_fields(flag,Natom,Mensemble)

      implicit none

      integer, intent(in) :: flag   !< Allocate or deallocate (1/-1)
      integer, intent(in), optional :: Natom !< Number of atoms in system
      integer, intent(in), optional :: Mensemble   !< Number of ensembles

      integer :: i_stat,i_all

      if (flag > 0) then
         if (.not. allocated(btorque_full)) then
            allocate(btorque_full(3,Natom,Mensemble), stat=i_stat)
            call memocc(i_stat, product(shape(btorque_full))*kind(btorque_full), 'btorque_full', 'allocate_aux_heun_fields')
         end if
         btorque_full = 0.0_dblprec
      else
         if (allocated(btorque_full)) then
            i_all = -product(shape(btorque_full))*kind(btorque_full)
            deallocate(btorque_full, stat=i_stat)
            call memocc(i_stat, i_all, 'btorque_full', 'allocate_aux_heun_fields')
         end if
      endif

   end subroutine allocate_aux_heun_fields

end module heun_proper
