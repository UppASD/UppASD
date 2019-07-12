!-------------------------------------------------------------------------------
! MODULE: VPO
!> @brief Routines needed to implement the velocity projection optimization (VPO) scheme
!> @author Pavel Bessarab
!-------------------------------------------------------------------------------
module VPO
   use Parameters
   use HamiltonianActions
   use math_functions
   implicit none

   public

contains

   !---------------------------------------------------------------------------------
   ! SUBROUTINE: normalize_vec
   !> @brief Normalize vector
   !> @note Jonathan Chico: Redundant with definition in math_functions, to be removed
   !> @author Pavel Bessarab
   !---------------------------------------------------------------------------------
   subroutine normalize_vec(vec_size,vec)

      implicit none

      integer, intent(in) :: vec_size  !! Size of the vector
      real(dblprec), dimension(vec_size), intent(inout) :: vec !! Input vector
      ! .. Local variables
      real(dblprec) :: tmp

      tmp = norm2(vec)
      vec(:) = vec(:)/tmp

   end subroutine normalize_vec

   !---------------------------------------------------------------------------------
   ! SUBROUTINE: convert_force
   !> @brief Convert effective fields to the format used in VPO
   !> @author Pavel Bessarab
   !---------------------------------------------------------------------------------
   subroutine convert_force(Natom,nim,mmom,emom,beff)

      implicit none

      integer, intent(in) :: nim !! Number of images
      integer, intent(in) :: Natom !! Number of atoms in system
      real(dblprec), dimension(Natom,nim), intent(in) :: mmom !! Magnitude of the magnetic moments
      real(dblprec), dimension(3,Natom,nim), intent(in) :: emom !! Current magnetic moment vector
      real(dblprec), dimension(3,Natom,nim), intent(inout) :: beff !! Current magnetic field
      ! .. Local variables
      integer :: ii,jj
      real(dblprec) :: tmp

      do jj=1,nim
         do ii=1,Natom
            tmp = 0.0_dblprec
            beff(:,ii,jj) = mmom(ii,jj)*beff(:,ii,jj)
            tmp = sum(beff(:,ii,jj)*emom(:,ii,jj))
            beff(:,ii,jj) = beff(:,ii,jj)-tmp*emom(:,ii,jj)
         end do
      end do
   end subroutine convert_force

   !-----------------------------------------------------------------------------
   ! SUBROUTINE: convert_mom
   !> @brief Convert normalized moments stored in the format used in VPO to emomM
   !> @author Pavel Bessarab
   !-----------------------------------------------------------------------------
   subroutine convert_mom(Natom,coo,mmom,emomM)
      implicit none
      integer, intent(in) :: Natom !< Number of atoms in system
      real(dblprec), dimension(Natom), intent(in) :: mmom !< Magnitude of the magnetic moments
      real(dblprec), dimension(3*Natom), intent(in) :: coo
      real(dblprec), dimension(3,Natom), intent(out) :: emomM !< Current magnetic moment vector
      ! .. Local variables
      integer :: ii

      !$omp parallel do default(shared), private(ii)
      do ii=1,Natom
         emomM(1,ii) = coo(3*ii-2)*mmom(ii)
         emomM(2,ii) = coo(3*ii-1)*mmom(ii)
         emomM(3,ii) = coo(3*ii)*mmom(ii)
      end do
      !$omp end parallel do

   end subroutine convert_mom

   !-----------------------------------------------------------------------------
   ! SUBROUTINE: calc_axis_all
   !> @brief Calculate axes of rotation based on magnetic configuration and effective fields
   !> @author Pavel Bessarab
   !-----------------------------------------------------------------------------
   subroutine calc_axis_all(Natom,m_in,f_in,ax)
      implicit none
      integer, intent(in) :: Natom !< Number of atoms in system
      real(dblprec), dimension(3,Natom) ,intent(in) :: m_in
      real(dblprec), dimension(3,Natom) ,intent(in) :: f_in
      real(dblprec), dimension(3,Natom), intent(out) :: ax
      ! .. Local variables
      integer :: ii

      !$omp parallel do default(shared), private(ii)
      do ii=1,Natom
         ax(:,ii) = calc_axis(m_in(:,ii),f_in(:,ii))
      end do
      !$omp end parallel do
   end subroutine calc_axis_all

   !---------------------------------------------------------------------------------
   ! FUNCTION: calc_axis
   !> @brief Calculate axis of rotation based on two 3-vectors
   !> @author Pavel Bessarab
   !---------------------------------------------------------------------------------
   function calc_axis(m_in,f_in)

      implicit none

      real(dblprec), dimension(3), intent(in) :: m_in
      real(dblprec), dimension(3), intent(in) :: f_in
      real(dblprec), dimension(3) :: calc_axis
      real(dblprec) :: tmp,a,b,c,eps=epsilon(a)
      real(dblprec), dimension(3) :: xx,yy

      tmp = norm2(f_in)

      if (tmp<eps) then
         a = 1.0_dblprec
         b = 1.0_dblprec
         c = 1.0_dblprec
      else
         a = f_in(1)
         b = f_in(2)
         c = f_in(3)
      end if

      yy(1) = sign(a,f_in(1))
      yy(2) = sign(b,f_in(2))
      yy(3) = sign(c,f_in(3))

      xx = f_cross_product(m_in,yy)
      calc_axis(:) = f_normalize_vec(xx,3)

   end function calc_axis

   !---------------------------------------------------------------------------------
   ! SUBROUTINE: quick_min_coo
   !> @brief Update coordinates
   !> @author Pavel Bessarab
   !---------------------------------------------------------------------------------
   subroutine quick_min_coo(Natom,poi_in,vel,f,poi_out,mass,dt)

      implicit none

      integer, intent(in) :: Natom !< Number of atoms in system
      real(dblprec), intent(in) :: dt
      real(dblprec), intent(in) :: mass
      real(dblprec), dimension(3,Natom), intent(in) :: f
      real(dblprec), dimension(3,Natom), intent(in) :: vel
      real(dblprec), dimension(3,Natom), intent(in) :: poi_in
      ! .. Output variables
      real(dblprec), dimension(3,Natom),intent(out) :: poi_out
      ! .. Local variables
      integer :: ii,jj
      real(dblprec) :: mmom

      !$omp parallel do default(shared), private(ii,jj,mmom)
      do ii=1,Natom
         mmom = norm2(poi_in(:,ii))
         do jj=1,3
            poi_out(jj,ii) = poi_in(jj,ii)+(vel(jj,ii)+0.50_dblprec*f(jj,ii)/mass*dt)*dt*mmom
         end do
         poi_out(:,ii)=f_normalize_vec(poi_out(:,ii),3)
         do jj=1,3
            poi_out(jj,ii) = poi_out(jj,ii)*mmom
         end do
      end do
      !$omp end parallel do

   end subroutine quick_min_coo

   !---------------------------------------------------------------------------------
   ! SUBROUTINE: quick_min_vel
   !> @brief Update velocities
   !> @author Pavel Bessarab
   !---------------------------------------------------------------------------------
   subroutine quick_min_vel(Natom,vel,f1,f2,ax,ang,mass,dt)

      implicit none

      integer, intent(in) :: Natom !< Number of atoms in system
      real(dblprec), intent(in) :: dt
      real(dblprec), intent(in) :: mass
      real(dblprec), dimension(Natom), intent(in) :: ang
      real(dblprec), dimension(3,Natom), intent(in) :: f1
      real(dblprec), dimension(3,Natom), intent(in) :: f2
      real(dblprec), dimension(3,Natom), intent(in) :: ax
      real(dblprec), dimension(3,Natom), intent(inout) :: vel
      ! .. Local variables
      integer :: ii,jj
      real(dblprec), dimension(:,:), allocatable :: rv,rf

      allocate(rv(3,Natom),rf(3,Natom))

      call rotate_vec_all(Natom,vel,ax,ang,rv)
      call rotate_vec_all(Natom,f1,ax,ang,rf)

      !$omp parallel do default(shared), private(ii,jj), collapse(2)
      do ii=1,Natom
         do jj=1,3
            vel(jj,ii) = rv(jj,ii)+0.50_dblprec*(rf(jj,ii)+f2(jj,ii))/mass*dt
         end do
      end do
      !$omp end parallel do

      deallocate(rv,rf)

   end subroutine quick_min_vel

   !---------------------------------------------------------------------------------
   ! SUBROUTINE: rotate_vec_all
   !> @brief Rotate all 3-vectors stored in 3*Natom-array around axes ax by angle ang
   !> @author Pavel Bessarab
   !---------------------------------------------------------------------------------
   subroutine rotate_vec_all(Natom,v_in,ax,ang,v_out)
      implicit none
      integer, intent(in) :: Natom !< Number of atoms in system
      real(dblprec), dimension(Natom), intent(in) :: ang
      real(dblprec), dimension(3,Natom), intent(in) :: ax
      real(dblprec), dimension(3,Natom), intent(in) :: v_in
      real(dblprec), dimension(3,Natom), intent(out) :: v_out
      ! .. Local variables
      integer :: ii

      !$omp parallel do default(shared), private(ii)
      do ii=1,Natom
         call rotate_vec(v_in(:,ii),ax(:,ii),ang(ii),v_out(:,ii))
      end do
      !$omp end parallel do
   end subroutine rotate_vec_all

   !---------------------------------------------------------------------------------
   ! SUBROUTINE: rotate_vec
   !> @brief Rotate 3-vector v_in by angle ang around axis ax
   !> @author Pavel Bessarab
   !---------------------------------------------------------------------------------
   subroutine rotate_vec(v_in,ax,ang,v_out)

      implicit none

      real(dblprec), intent(in) :: ang
      real(dblprec), dimension(3), intent(in) :: ax
      real(dblprec), dimension(3), intent(in) :: v_in
      real(dblprec), dimension(3), intent(out) :: v_out
      ! .. Local variables
      real(dblprec) :: tmp,sinang,cosang

      tmp = 0.0_dblprec
      tmp = sum(ax(:)*v_in(:))

      sinang = sin(ang)
      cosang = cos(ang)

      v_out(1) = v_in(1)*cosang + sinang*(ax(2)*v_in(3)-ax(3)*v_in(2))+(1.0_dblprec-cosang)*tmp*ax(1)
      v_out(2) = v_in(2)*cosang - sinang*(ax(1)*v_in(3)-ax(3)*v_in(1))+(1.0_dblprec-cosang)*tmp*ax(2)
      v_out(3) = v_in(3)*cosang + sinang*(ax(1)*v_in(2)-ax(2)*v_in(1))+(1.0_dblprec-cosang)*tmp*ax(3)
   end subroutine rotate_vec

   !---------------------------------------------------------------------------------
   ! SUBROUTINE: calc_ang_all
   !> @brief Calculate angle between each pair of 3-vectors stored in 3*Natom-arrays
   !> @author Pavel Bessarab
   !---------------------------------------------------------------------------------
   subroutine calc_ang_all(Natom,n1,n2,ang)

      implicit none

      integer, intent(in) :: Natom !< Number of atoms in system
      real(dblprec), dimension(3,Natom), intent(in) :: n1
      real(dblprec), dimension(3,Natom), intent(in) :: n2
      real(dblprec), dimension(Natom),intent(out) :: ang
      ! .. Local variables
      integer :: ii

      !$omp parallel do default(shared), private(ii)
      do ii=1,Natom
         ang(ii) = calc_ang(n1(:,ii),n2(:,ii))
      end do
      !$omp end parallel do
   end subroutine calc_ang_all

   !---------------------------------------------------------------------------------
   ! FUNCTION: calc_ang
   !> @brief Calculate angle between two 3-vectors
   !> @author Pavel Bessarab
   !---------------------------------------------------------------------------------
   real(dblprec) function calc_ang(n1,n2)
      implicit none
      real(dblprec), dimension(3), intent(in) :: n1 !n1 and n2 have to be normalized
      real(dblprec), dimension(3), intent(in) :: n2 !n1 and n2 have to be normalized
      ! .. Local variables
      real(dblprec) :: prod,nn
      real(dblprec), dimension(3) :: normal

      normal = f_cross_product(n1,n2)

      nn = norm2(normal)
      prod = sum(n1(:)*n2(:))

      calc_ang = atan2(nn,prod)
      return

   end function calc_ang

   !---------------------------------------------------------------------------------
   ! SUBROUTINE: project_vel
   !> @brief Project velocity on the direction of force
   !> @author Pavel Bessarab
   !---------------------------------------------------------------------------------
   subroutine project_vel(Natom,vel,f)

      implicit none

      integer, intent(in) :: Natom !< Number of atoms in system
      real(dblprec), dimension(Natom), intent(in) :: f
      real(dblprec), dimension(Natom), intent(inout) :: vel
      ! .. Local variables
      real(dblprec) :: fv,fd
      integer :: ii

      fv = 0.0_dblprec
      fd = 0.0_dblprec

      fv = sum(vel(:)*f(:))
      fd = sum(f(:)*f(:))
      if (fv<0.0_dblprec) then
         vel(:)=0.0_dblprec
      else
         vel(:)=f(:)*fv/fd
      end if
   end subroutine project_vel

   !-----------------------------------------------------------------------------
   ! SUBROUTINE: project_vel
   !> @brief Project velocity on the direction of force (used with GNEB)
   !> @author Pavel Bessarab
   !-----------------------------------------------------------------------------
   subroutine project_vel_gneb(Natom,nim,vel,f)
      implicit none
      integer, intent(in) :: nim !< Number of images 
      integer, intent(in) :: Natom !< Number of atoms in system
      real(dblprec), dimension(3,Natom,nim), intent(in) :: f
      real(dblprec), dimension(3,Natom,nim), intent(inout) :: vel
      ! .. Local variables
      real(dblprec) :: fv,fd
      integer :: ii,jj,kk

      fv = 0.0_dblprec
      fd = 0.0_dblprec

      !$omp parallel do default(shared), private(ii,jj,kk), reduction(+:fv,fd), collapse(3)
      do ii=1,nim
         do jj=1,Natom
            do kk=1,3
               fv = fv + vel(kk,jj,ii)*f(kk,jj,ii)
               fd = fd+f(kk,jj,ii)*f(kk,jj,ii)
            end do
         end do
      end do
      !$omp end parallel do
      if (fv<0.0_dblprec) then
         vel(:,:,:) = 0.0_dblprec
      else
         vel(:,:,:) = f(:,:,:)*fv/fd
      end if
   end subroutine project_vel_gneb

   !-----------------------------------------------------------------------------
   ! SUBROUTINE: normalize_vec_all
   !> @brief Normalize each 3-vector stored in the 3*Natom-array
   !> @author Pavel Bessarab
   !-----------------------------------------------------------------------------
   subroutine normalize_vec_all(Natom,vec)
      implicit none
      integer, intent(in) :: Natom !< Number of atoms in system
      real(dblprec), dimension(3*Natom), intent(inout) :: vec !< Vector to be normalized
      ! .. Local variables
      real(dblprec) :: tmp
      integer :: ii

      !$omp parallel do default(shared), private(ii,tmp)
      do ii=1,Natom
         tmp = sqrt(vec(3*ii-2)*vec(3*ii-2)+vec(3*ii-1)*vec(3*ii-1)+vec(3*ii)*vec(3*ii))
         vec(3*ii-2) = vec(3*ii-2)/tmp
         vec(3*ii-1) = vec(3*ii-1)/tmp
         vec(3*ii)   = vec(3*ii)/tmp
      end do
      !$omp end parallel do
   end subroutine normalize_vec_all

end module VPO
