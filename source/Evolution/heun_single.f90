!-------------------------------------------------------------------------------
! MODULE: Heun_single
!> @brief A single step Heun solver for the LLG-equations
!> @copyright
!> GNU Public License.
!-------------------------------------------------------------------------------
module Heun_single
   use Parameters
   use Profiling

   implicit none

   private
   public :: evolve2, evolve3


contains

   !-----------------------------------------------------------------------------
   ! SUBROUTINE: evolve2
   !> @brief The evolution routine for single step Heun
   !-----------------------------------------------------------------------------
   subroutine evolve2(Natom,Mensemble,Landeg,llg,bn,lambda1_array,lambda2_array,    &
      beff,beff2,field1,field2,mmomi,emom2,compensate_drift,deltat,Temp_array,      &
      temprescale,thermal_field,Nred,red_atom_list)
      !
      use Constants
      use RandomNumbers, only : ranv, rng_gaussian
      !
      implicit none

      integer, intent(in) :: llg                !< Type of equation of motion (1=LLG)
      integer, intent(in) :: Nred
      integer, intent(in) :: Natom              !< Number of atoms in system
      integer, intent(in) :: Mensemble          !< Number of ensembles
      integer, intent(in) :: compensate_drift   !< Correct for drift in RNG
      real(dblprec), intent(in) :: bn           !< Scaling factor for LLG equation (parameter=1.0_dblprec)
      real(dblprec), intent(in) :: deltat       !< Time step
      real(dblprec), intent(in) :: temprescale  !< Temperature rescaling from QHB
      integer, dimension(Nred), intent(in) :: red_atom_list
      real(dblprec), dimension(Natom), intent(in) :: Landeg          !< Gyromagnetic ratio
      real(dblprec), dimension(Natom), intent(in) :: Temp_array      !< Temperature
      real(dblprec), dimension(Natom), intent(in) :: lambda1_array   !< Damping parameter
      real(dblprec), dimension(Natom), intent(in) :: lambda2_array   !< Additional damping parameter (not used for llg=1)
      real(dblprec), dimension(3,Mensemble), intent(in) :: field1    !< Average internal effective field
      real(dblprec), dimension(3,Mensemble), intent(in) :: field2    !< Average external effective field
      real(dblprec), dimension(Natom,Mensemble), intent(in) :: mmomi !< Inverse of magnitude of magnetic moments
      real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: beff   !< Total effective field from application of Hamiltonian
      real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: beff2  !< External field from application of Hamiltonian
      real(dblprec), dimension(3,Natom,Mensemble), intent(inout) :: emom2  !< Final (or temporary) unit moment vector
      real(dblprec), dimension(3,Natom,Mensemble), intent(inout), optional :: thermal_field

      real(dblprec) :: D
      real(dblprec) :: etot
      real(dblprec) :: dt, mu, sigma, avf1, avf2
      real(dblprec) :: tk1,tk2
      real(dblprec) :: rx,ry,rz
      real(dblprec), dimension(Natom) :: Dk

      integer :: i, ired

      !JohanH-Dec23: Temporarily evolve2 and heun2
      !work only for Mensemble=1
      avf1=norm2(field1(:,1))
      avf2=norm2(field2(:,1))

      ! LL equations ONE universal damping
      Dk=0.0_dblprec
      if (llg==0) then
         Dk=(lambda1_array*k_bolt/gama/(mub))*(gama/bn)  !last factor for dim. less.
         dt=deltat*bn*gama !dim. less time
         tk1=1.0_dblprec/bn
         tk2=0.0_dblprec

         ! LLG equations ONE universal damping
      else if (llg==1) then
         Dk=(lambda1_array/(1+lambda1_array**2)*k_bolt/gama/(mub))*(gama/bn)  !last factor for dim. less.
         dt=deltat*bn*gama !dim. less time
         tk1=1.0_dblprec/bn
         tk2=0.0_dblprec

         ! LL equations TWO damping parameters
      else if(llg==2) then
         Dk=(lambda1_array*avf1+lambda2_array*avf2)/(avf1+avf2)*(k_bolt/gama/(mub))*(gama/bn)
         !last factor for dim. less.
         dt=deltat*bn*gama !dim. less time
         tk1=1.0_dblprec/bn
         tk2=1.0_dblprec/bn

         ! LLG equations TWO damping parameters, fluctuations included in damping1
      else if(llg==3) then
         Dk=(lambda1_array*avf1+lambda2_array*avf2)/(avf1+avf2)/(1+lambda1_array**2)*(k_bolt/gama/(mub))*(gama/bn)
         !last factor for dim. less.
         dt=deltat*bn*gama !dim. less time
         tk1=1.0_dblprec/bn
         tk2=1.0_dblprec/bn

         ! LL equations TWO damping parameters, but use thermal fluctuations corresponding to one universal damping
      else if(llg==4) then
         Dk=(lambda1_array*k_bolt/gama/(mub))*(gama/bn)  !last factor for dim. less.
         dt=deltat*bn*gama !dim. less time
         tk1=1.0_dblprec/bn
         tk2=1.0_dblprec/bn

         ! LLG equations TWO damping parameters, but use thermal fluctuations corresponding to one universl damping
      else if(llg==5) then
         Dk=(lambda1_array/(1+lambda1_array**2)*k_bolt/gama/(mub))*(gama/bn)  !last factor for dim. less.
         dt=deltat*bn*gama !dim. less time
         tk1=1.0_dblprec/bn
         tk2=0.0_dblprec

         ! LL equations ONE universal damping write energies
      else if(llg==6) then
         Dk=(lambda1_array*k_bolt/gama/(mub))*(gama/bn)  !last factor for dim. less.
         dt=deltat*bn*gama !dim. less time

         ! LLG equations ONE universal damping write energies
      else if (llg==7) then
         Dk=(lambda1_array/(1.0_dblprec+lambda1_array**2)*k_bolt/gama/(mub))*(gama/bn) !last factor for dim. less.
         dt=deltat*bn*gama !dim. less time
         tk1=1.0_dblprec/bn
         tk2=0.0_dblprec
      endif

      call rng_gaussian(ranv,3*Natom*Mensemble,1.0_dblprec)

      do ired=1, Nred
         i=red_atom_list(ired)
         D=Dk(i)*mmomi(i,1)*Temp_array(i)*temprescale
         sigma=sqrt(dt*2*D)
         mu=0.0_dblprec
         ranv(1,i,1)=ranv(1,i,1)*sigma
         ranv(2,i,1)=ranv(2,i,1)*sigma
         ranv(3,i,1)=ranv(3,i,1)*sigma
      end do

      if(compensate_drift==1) then
         rx=0.0_dblprec;ry=0.0_dblprec;rz=0.0_dblprec
         do ired=1, Nred
            i=red_atom_list(ired)
            rx=rx+ranv(1,i,1)
            ry=ry+ranv(2,i,1)
            rz=rz+ranv(3,i,1)
         end do
         rx=rx/Natom
         ry=ry/Natom
         rz=rz/Natom
         do ired=1, Nred
            i=red_atom_list(ired)
            ranv(1,i,1)=ranv(1,i,1)-rx
            ranv(2,i,1)=ranv(2,i,1)-ry
            ranv(3,i,1)=ranv(3,i,1)-rz
         end do
      end if

      call heun2(Natom,Mensemble,Landeg,lambda1_array,lambda2_array,llg,beff,beff2, &
         emom2,dt,tk1,tk2,Nred,red_atom_list)

      !$omp parallel do schedule(static) default(shared) private(i,etot)
      do i=1,Natom
         etot=1/norm2(emom2(:,i,1))
         emom2(1,i,1)=emom2(1,i,1)*etot
         emom2(2,i,1)=emom2(2,i,1)*etot
         emom2(3,i,1)=emom2(3,i,1)*etot
      end do
      !$omp end parallel do

      thermal_field=ranv

   end subroutine evolve2

   !-----------------------------------------------------------------------------
   ! SUBROUTINE: evolve3
   !> @brief An alternative evolution routine for single step Heun
   !-----------------------------------------------------------------------------
   subroutine evolve3(Natom,Mensemble,Landeg,llg,lambda1_array,beff,field1,field2,  &
      mmomi,emom2,compensate_drift,deltat,Temp_array,temprescale,thermal_field,Nred,&
      red_atom_list)

      use Constants
      use RandomNumbers, only : rng_norm, ranv

      implicit none

      integer, intent(in) :: llg 	!< Type of equation of motion (1=LLG)
      integer, intent(in) :: Nred	!< Number of moments that can be updated
      integer, intent(in) :: Natom 	!< Number of atoms in system
      integer, intent(in) :: Mensemble !< Number of ensembles
      integer, intent(in) :: compensate_drift 	!< Correct for drift in RNG
      real(dblprec), intent(in) :: deltat 		!< Time step
      real(dblprec), intent(in) :: temprescale  !< Temperature rescaling from QHB
      integer, dimension(Nred), intent(in) :: red_atom_list	!< Reduced list containing atoms allowed to evolve in a fixed moment calculation
      real(dblprec), dimension(Natom), intent(in) :: Landeg  !< Gyromagnetic ratio
      real(dblprec), dimension(Natom), intent(in) :: Temp_array  !< Temperature
      real(dblprec), dimension(Natom), intent(in) :: lambda1_array !< Damping parameter
      real(dblprec), dimension(3,Mensemble), intent(in) :: field1 !< Average internal effective field
      real(dblprec), dimension(3,Mensemble), intent(in) :: field2 !< Average external effective field
      real(dblprec), dimension(Natom,Mensemble), intent(in) :: mmomi !< Inverse of magnitude of magnetic moments
      real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: beff !< Total effective field from application of Hamiltonian
      !.. In/out variables
      real(dblprec), dimension(3,Natom,Mensemble), intent(inout) :: emom2  !< Final (or temporary) unit moment vector
      real(dblprec), dimension(3,Natom,Mensemble), intent(inout), optional :: thermal_field	!< Thermal field

      real(dblprec) :: D
      real(dblprec) :: etot
      real(dblprec) :: dt, mu, sigma, avf1, avf2
      real(dblprec) :: tk1
      real(dblprec) :: rx,ry,rz
      real(dblprec), dimension(Natom) :: Dk

      integer :: i,ired

      !JohanH-Dec23: Temporarily evolve3 and heun3
      !work only for Mensemble=1
      avf1=norm2(field1(:,1))
      avf2=norm2(field2(:,1))

      ! LLG equations ONE universal damping
      Dk=0.0_dblprec
      if (llg==1) then
         Dk=lambda1_array/(1+lambda1_array**2)*k_bolt/mub
         dt=deltat*gama !rescaled time
         tk1=1.0_dblprec
      else
         stop 'Only llg=1 is currently supported'
      endif

		do ired=1, Nred
			i=red_atom_list(ired)
         D=Dk(i)*mmomi(i,1)*Temp_array(i)*temprescale
         sigma=sqrt(dt*2*D)
         mu=0_dblprec
         ranv(1,i,1)=rng_norm(mu,sigma)
         ranv(2,i,1)=rng_norm(mu,sigma)
         ranv(3,i,1)=rng_norm(mu,sigma)
      end do

      if(compensate_drift==1) then
         rx=0.0_dblprec;ry=0.0_dblprec;rz=0.0_dblprec
			do ired=1, Nred
				i=red_atom_list(ired)
            rx=rx+ranv(1,i,1)
            ry=ry+ranv(2,i,1)
            rz=rz+ranv(3,i,1)
         end do
         rx=rx/Natom
         ry=ry/Natom
         rz=rz/Natom
			do ired=1, Nred
				i=red_atom_list(ired)
            ranv(1,i,1)=ranv(1,i,1)-rx
            ranv(2,i,1)=ranv(2,i,1)-ry
            ranv(3,i,1)=ranv(3,i,1)-rz
         end do
      end if

		call heun3(Natom, Mensemble, Landeg, lambda1_array, beff, emom2, dt,Nred,		&
			red_atom_list)

      !$omp parallel do schedule(static) default(shared) private(i,etot)
      do i=1,Natom
         etot=1.0_dblprec/norm2(emom2(:,i,1))
         emom2(1,i,1)=emom2(1,i,1)*etot
         emom2(2,i,1)=emom2(2,i,1)*etot
         emom2(3,i,1)=emom2(3,i,1)*etot
      end do
      !$omp end parallel do

      thermal_field=ranv

   end subroutine evolve3

   !-----------------------------------------------------------------------------
   ! SUBROUTINE: heun2
   !> @brief Perform a single Heun step
   !-----------------------------------------------------------------------------
   subroutine heun2(Natom,Mensemble,Landeg,lambda1_array,lambda2_array,llg,beff,    &
      beff2,emom2,dt,b2h1,b2h2,Nred,red_atom_list)
      !
      use RandomNumbers, only : ranv
      !
      implicit none

		integer, intent(in) :: llg 	!< Type of equation of motion (1=LLG)
		integer, intent(in) :: Nred	!< Number of moments that can be updated
		integer, intent(in) :: Natom 	!< Number of atoms in system
      integer, intent(in) :: Mensemble !< Number of ensembles
      real(dblprec), intent(in) :: dt !< Time step
      real(dblprec), intent(in) :: b2h1 !< Scale factor for effective field
		real(dblprec), intent(in) :: b2h2 !< Scale factor for applied field
		integer, dimension(Nred), intent(in) :: red_atom_list	!< Reduced list containing atoms allowed to evolve in a fixed moment calculation
      real(dblprec), dimension(Natom), intent(in) :: Landeg !< Gyromagnetic ratio
      real(dblprec), dimension(Natom), intent(in) :: lambda1_array !< Damping parameter
      real(dblprec), dimension(Natom), intent(in) :: lambda2_array !< Additional damping parameter (not used for llg=1)
      real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: beff !< Total effective field from application of Hamiltonian
      real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: beff2 !< External field from application of Hamiltonian
      real(dblprec), dimension(3,Natom,Mensemble), intent(inout) :: emom2  !< Final (or temporary) unit moment vector

      real(dblprec) :: b1xx, b1xy, b1xz
      real(dblprec) :: b1yx, b1yy, b1yz
      real(dblprec) :: b1zx, b1zy, b1zz
      real(dblprec) :: a1x, a1y, a1z
      real(dblprec) :: etx, ety, etz
      real(dblprec) :: dtg
      !
      real(dblprec) :: dwx, dwy, dwz
      real(dblprec) :: e1x, e1y, e1z
      real(dblprec) :: hlpx,hlpy,hlpz
      real(dblprec) :: hesx,hesy,hesz,prot
      real(dblprec), dimension(Natom) :: clambda1

      ! ... Local variables ...
      integer :: i, ired
      !
      if(llg==1.or.llg==3.or.llg==5.or.llg==7) then
         clambda1=lambda1_array
      else
         clambda1=0
      end if

      !$omp parallel do default(shared) &
      !$omp private(i,ired,e1x,e1y,e1z,dwx,dwy,dwz,hesx,hesy,hesz &
      !$omp ,hlpx,hlpy,hlpz,etx,ety,etz,a1x,a1y,a1z,prot,b1xx &
      !$omp ,b1xy,b1xz,b1yx,b1yy,b1yz,b1zx,b1zy,b1zz,dtg)
		do ired=1,Nred
			i=red_atom_list(ired)
         e1x=emom2(1,i,1)
         e1y=emom2(2,i,1)
         e1z=emom2(3,i,1)
         dtg=dt*Landeg(i)
         dwx=ranv(1,i,1)*Landeg(i)
         dwy=ranv(2,i,1)*Landeg(i)
         dwz=ranv(3,i,1)*Landeg(i)
         hesx=beff(1,i,1)*b2h1+beff2(1,i,1)*b2h2
         hlpx=lambda1_array(i)*beff(1,i,1)*b2h1+lambda2_array(i)*beff2(1,i,1)*b2h2
         hesy=beff(2,i,1)*b2h1+beff2(2,i,1)*b2h2
         hlpy=lambda1_array(i)*beff(2,i,1)*b2h1+lambda2_array(i)*beff2(2,i,1)*b2h2
         hesz=beff(3,i,1)*b2h1+beff2(3,i,1)*b2h2
         hlpz=lambda1_array(i)*beff(3,i,1)*b2h1+lambda2_array(i)*beff2(3,i,1)*b2h2
         !
         prot=e1x*hlpx+e1y*hlpy+e1z*hlpz
         a1x=hlpx+e1z*hesy-e1y*hesz-e1x*prot
         a1y=hlpy-e1z*hesx+e1x*hesz-e1y*prot
         a1z=hlpz+e1y*hesx-e1x*hesy-e1z*prot
         !
         b1xx=          clambda1(i)-e1x*e1x*clambda1(i)
         b1xy=e1z-e1x*e1y*clambda1(i)
         b1xz=-e1y-e1x*e1z*clambda1(i)
         b1yx=-e1z-e1y*e1x*clambda1(i)
         b1yy=          clambda1(i)-e1y*e1y*clambda1(i)
         b1yz=e1x-e1y*e1z*clambda1(i)
         b1zx=e1y-e1z*e1x*clambda1(i)
         b1zy=-e1x-e1z*e1y*clambda1(i)
         b1zz=          clambda1(i)-e1z*e1z*clambda1(i)
         !
         etx=e1x+a1x*dtg+b1xx*dwx+b1xy*dwy+b1xz*dwz
         ety=e1y+a1y*dtg+b1yx*dwx+b1yy*dwy+b1yz*dwz
         etz=e1z+a1z*dtg+b1zx*dwx+b1zy*dwy+b1zz*dwz
         !
         prot=etx*hlpx+ety*hlpy+etz*hlpz
         a1x=a1x+hlpx+etz*hesy-ety*hesz-etx*prot
         a1y=a1y+hlpy-etz*hesx+etx*hesz-ety*prot
         a1z=a1z+hlpz+ety*hesx-etx*hesy-etz*prot
         !
         b1xx=b1xx+          clambda1(i)-etx*etx*clambda1(i)
         b1xy=b1xy+etz-etx*ety*clambda1(i)
         b1xz=b1xz-ety-etx*etz*clambda1(i)
         b1yx=b1yx-etz-ety*etx*clambda1(i)
         b1yy=b1yy+          clambda1(i)-ety*ety*clambda1(i)
         b1yz=b1yz+etx-ety*etz*clambda1(i)
         b1zx=b1zx+ety-etz*etx*clambda1(i)
         b1zy=b1zy-etx-etz*ety*clambda1(i)
         b1zz=b1zz+          clambda1(i)-etz*etz*clambda1(i)
         !
         etx=a1x*dtg+b1xx*dwx+b1xy*dwy+b1xz*dwz
         ety=a1y*dtg+b1yx*dwx+b1yy*dwy+b1yz*dwz
         etz=a1z*dtg+b1zx*dwx+b1zy*dwy+b1zz*dwz

         emom2(1,i,1)=e1x+0.5_dblprec*etx
         emom2(2,i,1)=e1y+0.5_dblprec*ety
         emom2(3,i,1)=e1z+0.5_dblprec*etz
      end do
      !$omp end parallel do

      return

   end subroutine heun2

   !-----------------------------------------------------------------------------
   ! SUBROUTINE: heun3
   !> @brief Perform a single Heun step
   !-----------------------------------------------------------------------------
   subroutine heun3(Natom,Mensemble,Landeg,lambda1_array,beff,emom2,dt,Nred,red_atom_list)

      use RandomNumbers, only : ranv

      implicit none

		integer, intent(in) :: Nred	!< Number of moments that can be updated
      integer, intent(in) :: Natom 	!< Number of atoms in system
      integer, intent(in) :: Mensemble !< Number of ensembles
      real(dblprec), intent(in) :: dt 	!< Time step
		integer, dimension(Nred), intent(in) :: red_atom_list		!< Reduced list containing atoms allowed to evolve in a fixed moment calculation
		real(dblprec), dimension(Natom), intent(in) :: Landeg  	!< Gyromagnetic ratio
      real(dblprec), dimension(Natom), intent(in) :: lambda1_array 		!< Damping parameter
      real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: beff 	!< Total effective field from application of Hamiltonian
      real(dblprec), dimension(3,Natom,Mensemble), intent(inout) :: emom2  !< Final (or temporary) unit moment vector

      real(dblprec) :: b1xx, b1xy, b1xz
      real(dblprec) :: b1yx, b1yy, b1yz
      real(dblprec) :: b1zx, b1zy, b1zz
      real(dblprec) :: a1x, a1y, a1z
      real(dblprec) :: etx, ety, etz
      real(dblprec) :: dtg
      !
      real(dblprec) :: dwx, dwy, dwz
      real(dblprec) :: e1x, e1y, e1z
      real(dblprec) :: hlpx,hlpy,hlpz
      real(dblprec) :: hesx,hesy,hesz,prot
      !
      ! ... Local variables ...
      integer :: i,ired

      !$omp parallel do default(shared) &
      !$omp private(i,ired,e1x,e1y,e1z,dwx,dwy,dwz,hesx,hesy,hesz &
      !$omp ,hlpx,hlpy,hlpz,etx,ety,etz,a1x,a1y,a1z,prot,b1xx &
      !$omp ,b1xy,b1xz,b1yx,b1yy,b1yz,b1zx,b1zy,b1zz,dtg)
      do ired=1,Nred
        	i=red_atom_list(ired)
         e1x=emom2(1,i,1)
         e1y=emom2(2,i,1)
         e1z=emom2(3,i,1)
         dtg=dt*Landeg(i)
         dwx=ranv(1,i,1)*Landeg(i)
         dwy=ranv(2,i,1)*Landeg(i)
         dwz=ranv(3,i,1)*Landeg(i)

         !eventually variables hesx, hesy, hesz
         !can be removed
         hesx=beff(1,i,1)
         hlpx=lambda1_array(i)*beff(1,i,1)
         hesy=beff(2,i,1)
         hlpy=lambda1_array(i)*beff(2,i,1)
         hesz=beff(3,i,1)
         hlpz=lambda1_array(i)*beff(3,i,1)

         prot=e1x*hlpx+e1y*hlpy+e1z*hlpz
         a1x=hlpx+e1z*hesy-e1y*hesz-e1x*prot
         a1y=hlpy-e1z*hesx+e1x*hesz-e1y*prot
         a1z=hlpz+e1y*hesx-e1x*hesy-e1z*prot
         !
         b1xx=          lambda1_array(i)-e1x*e1x*lambda1_array(i)
         b1xy= e1z-e1x*e1y*lambda1_array(i)
         b1xz=-e1y-e1x*e1z*lambda1_array(i)
         b1yx=-e1z-e1y*e1x*lambda1_array(i)
         b1yy=          lambda1_array(i)-e1y*e1y*lambda1_array(i)
         b1yz= e1x-e1y*e1z*lambda1_array(i)
         b1zx= e1y-e1z*e1x*lambda1_array(i)
         b1zy=-e1x-e1z*e1y*lambda1_array(i)
         b1zz=          lambda1_array(i)-e1z*e1z*lambda1_array(i)
         !
         etx=e1x+a1x*dtg+b1xx*dwx+b1xy*dwy+b1xz*dwz
         ety=e1y+a1y*dtg+b1yx*dwx+b1yy*dwy+b1yz*dwz
         etz=e1z+a1z*dtg+b1zx*dwx+b1zy*dwy+b1zz*dwz
         !
         prot=etx*hlpx+ety*hlpy+etz*hlpz
         a1x=a1x+hlpx+etz*hesy-ety*hesz-etx*prot
         a1y=a1y+hlpy-etz*hesx+etx*hesz-ety*prot
         a1z=a1z+hlpz+ety*hesx-etx*hesy-etz*prot
         !
         b1xx=b1xx+          lambda1_array(i)-etx*etx*lambda1_array(i)
         b1xy=b1xy+etz-etx*ety*lambda1_array(i)
         b1xz=b1xz-ety-etx*etz*lambda1_array(i)
         b1yx=b1yx-etz-ety*etx*lambda1_array(i)
         b1yy=b1yy+          lambda1_array(i)-ety*ety*lambda1_array(i)
         b1yz=b1yz+etx-ety*etz*lambda1_array(i)
         b1zx=b1zx+ety-etz*etx*lambda1_array(i)
         b1zy=b1zy-etx-etz*ety*lambda1_array(i)
         b1zz=b1zz+          lambda1_array(i)-etz*etz*lambda1_array(i)
         !
         etx=a1x*dtg+b1xx*dwx+b1xy*dwy+b1xz*dwz
         ety=a1y*dtg+b1yx*dwx+b1yy*dwy+b1yz*dwz
         etz=a1z*dtg+b1zx*dwx+b1zy*dwy+b1zz*dwz
         emom2(1,i,1)=e1x+0.5_dblprec*etx
         emom2(2,i,1)=e1y+0.5_dblprec*ety
         emom2(3,i,1)=e1z+0.5_dblprec*etz
      end do
      !$omp end parallel do

      return

   end subroutine heun3

end module heun_single
