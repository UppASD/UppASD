!-------------------------------------------------------------------------------
! MODULE: PathInit
!> @brief Routines needed to initialize path for geodesic nudged elastic band calculation
!> @author Pavel Bessarab
!-------------------------------------------------------------------------------
module PathInit
   use Parameters
   use Constants
   implicit none

   private
   public :: save_path, save_sp, geodesic_path

contains

   !-------------------------------------------------------------------------------
   ! SUBROUTINE: geodesic_path_one
   !> @brief Calculation of the geodesic path from the intial and final states for a given atom
   !> @author Pavel Bessarab
   !-------------------------------------------------------------------------------
   subroutine geodesic_path_one(nim,ni,nf,path)

      use RandomNumbers, only : rng_uniform
      use math_functions
      use VPO, only: calc_ang,calc_axis

      implicit none

      integer, intent(in) :: nim !< Number of images
      real(dblprec), dimension(3), intent(in) :: ni !< Initial configuration for the i-th atom
      real(dblprec), dimension(3), intent(in) :: nf !< Final configuration for the i-th atom
      real(dblprec), dimension(3,nim), intent(out) :: path !< Path for for the i-th atom
      ! .. Local variables
      integer :: ii,jj
      real(dblprec) :: dtheta, theta, angle, pr, tmp, eps=epsilon(angle)
      real(dblprec), dimension(3) :: ax, vec
      real(dblprec), dimension(2) :: rn

      angle = calc_ang(ni,nf)
      if (angle<eps) then
         do ii = 1,3
            vec(ii) = (nf(ii) - ni(ii))/(nim-1)
         end do
         path(:,1) = ni(:)
         path(:,nim) = nf(:)
         do ii=2,nim-1
            tmp = 0.0_dblprec
            do jj=1,3
               path(jj,ii) = ni(jj) + real((ii-1),dblprec)*vec(jj)
               tmp = tmp + path(jj,ii)*path(jj,ii)
            end do
            tmp = sqrt(tmp)
            path(:,ii) = path(:,ii)/tmp
         end do

      elseif (abs(angle-pi)<eps) then
         tmp = 0.0_dblprec
         do while (tmp<eps)
            pr = 0.0_dblprec
            do jj=1,3
               call rng_uniform(rn,2)
               ax(jj) = sign(1.0_dblprec,2.0_dblprec*rn(1)-1.0_dblprec)*(rn(2)+1.0_dblprec)
               pr = pr + ax(jj)*ni(jj)
            end do
            ax(:) = ax(:) - pr*ni(:)
            tmp = norm2(ax)
         end do
         tmp = norm2(ax)
         ax = ax/tmp
         dtheta = pi/(nim-1)
         path(:,1) = ni(:)
         path(:,nim) = nf(:)

         do ii=2,nim-1
            theta = (ii-1)*dtheta
            path(1,ii) = ni(1)*cos(theta) + sin(theta)*(ax(2)*ni(3)-ax(3)*ni(2))
            path(2,ii) = ni(2)*cos(theta) - sin(theta)*(ax(1)*ni(3)-ax(3)*ni(1))
            path(3,ii) = ni(3)*cos(theta) + sin(theta)*(ax(1)*ni(2)-ax(2)*ni(1))
            path(:,ii) = f_normalize_vec(path(:,ii),3)
         end do

      else
         ax = calc_axis(ni,nf)
         dtheta = angle/(nim-1)

         path(:,1) = ni(:)
         path(:,nim) = nf(:)

         do ii=2,nim-1
            theta = (ii-1)*dtheta
            path(1,ii) = ni(1)*cos(theta) + sin(theta)*(ax(2)*ni(3)-ax(3)*ni(2))
            path(2,ii) = ni(2)*cos(theta) - sin(theta)*(ax(1)*ni(3)-ax(3)*ni(1))
            path(3,ii) = ni(3)*cos(theta) + sin(theta)*(ax(1)*ni(2)-ax(2)*ni(1))
            path(:,ii) = f_normalize_vec(path(:,ii),3)
         end do

      end if

   end subroutine geodesic_path_one

   !-------------------------------------------------------------------------------
   ! SUBROUTINE: geodesic_path
   !> @brief Calculation of geodesic path from the initial and fianl states
   !> @author Pavel Bessarab
   !-------------------------------------------------------------------------------
   subroutine geodesic_path(Natom,nim,ni,nf,amp_rnd,emom)
      use RandomNumbers, only : rng_uniform
      implicit none
      integer, intent(in) :: nim !< Number of images
      integer, intent(in) :: Natom  !< Number of atoms in system
      real(dblprec), intent(in) :: amp_rnd   !< Amplitude of the random noise
      real(dblprec), dimension(3,Natom), intent(in) :: nf !< Final configuration
      real(dblprec), dimension(3,Natom), intent(in) :: ni !< Initial configuration
      real(dblprec), dimension(3,Natom,nim), intent(out) :: emom  !< Current unit moment vector
      !.. Local variables
      integer :: ii,jj
      real(dblprec) :: mmom_tmp
      real(dblprec), dimension(3) :: u, v,rn

      do ii=1,Natom
         call geodesic_path_one(nim,ni(:,ii),nf(:,ii),emom(:,ii,:))
      end do

      do ii=1,nim
         do jj=1, Natom
            if ((ii==1).or.(ii==nim)) then
               u(:) = 0.0_dblprec
            else
               call rng_uniform(rn,3)
               u(:)=2.0_dblprec*(rn(:)-0.50_dblprec)
            end if

            v(1:3) = emom(1:3,jj,ii)+amp_rnd*u(1:3)
            mmom_tmp = norm2(v)
            emom(1:3,jj,ii) = v(1:3)/mmom_tmp
         end do
      end do

   end subroutine geodesic_path

   !-------------------------------------------------------------------------------
   ! SUBROUTINE: save_path
   !> @brief Print the path to file. Ensembles correspond to images in GNEB method
   !> @author Pavel Bessarab
   !-------------------------------------------------------------------------------
   subroutine save_path(Natom,Mensemble,simid,prn_mode,emom,mmom,mode,do_mom_legacy)
      use ErrorHandling,        only : ErrorHandling_ERROR
      use Restart, only : prn_mag_conf
      !.. Implicit declarations
      implicit none

      integer, intent(in) :: Natom     !< Number of atoms in system
      integer, intent(in) :: Mensemble !< Number of ensembles
      integer, intent(in) :: prn_mode  !< 0 - states before relaxation, 1 - states after relaxation
      character(len=1), intent(in) :: mode            !< Type of simulation
      character(len=1), intent(in) :: do_mom_legacy   !< Flag to print/read moments in legacy output
      character(len=8), intent(in) :: simid           !< Name of simulation
      real(dblprec), dimension(Natom,Mensemble), intent(in) :: mmom   !< Magnitude of the magnetic moment
      real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: emom !< Current unit moment vector

      integer :: ii, jj
      character(len=35) :: filn

      if (do_mom_legacy.ne.'Y') then
         if (prn_mode==0) then
            call prn_mag_conf(Natom,0,Mensemble,'M',simid,mmom,emom,'_path_in',mode)
         elseif (prn_mode==1) then
            call prn_mag_conf(Natom,0,Mensemble,'M',simid,mmom,emom,'_path',mode)
         elseif (prn_mode==2) then
            call prn_mag_conf(Natom,0,Mensemble,'M',simid,mmom,emom,'_if_restart',  &
               mode)
         elseif (prn_mode==3) then
            call prn_mag_conf(Natom,0,Mensemble,'M',simid,mmom,emom,'_path_restart',&
               mode)
         elseif (prn_mode==4) then
            call prn_mag_conf(Natom,0,Mensemble,'M',simid,mmom,emom,'_sp',mode)
         else
            write (*,*) "Error writing initial and final states to file"
            stop
         end if
      else
         if (prn_mode==0) then
            filn = 'moment_path_in.'//trim(adjustl(simid))//'.out'
         elseif (prn_mode==1) then
            filn = 'moment_path.'//trim(adjustl(simid))//'.out'
         elseif (prn_mode==2) then
            filn = 'moment_if_restart.'//trim(adjustl(simid))//'.out'
         elseif (prn_mode==3) then
            filn = 'moment_path_restart.'//trim(adjustl(simid))//'.out'
         elseif (prn_mode==4) then
            filn = 'moment_sp.'//trim(adjustl(simid))//'.out'
         else
            call ErrorHandling_ERROR("Error writing path to file!")
         end if
         open(4, file=filn, access = 'sequential',action = 'write', status = 'replace')
         do ii=1, Mensemble
            do jj=1, Natom
               write (4,10002) ii,jj,emom(1,jj,ii),emom(2,jj,ii),emom(3,jj,ii),mmom(jj,ii)
            end do
         end do
         close(4)
      endif
      return
      call ErrorHandling_ERROR("Error writing path to file!")
10002 format (i8,2x,i8,2x,2x, es16.8E3,2x,es16.8E3,2x,es16.8E3,2x,es16.8E3)

   end subroutine save_path

   !-------------------------------------------------------------------------------
   ! SUBROUTINE: save_sp
   !> @brief Print the SP configuration to file. Ensembles correspond to images in GNEB method
   !> @author Pavel Bessarab
   !-------------------------------------------------------------------------------
   subroutine save_sp(Natom, simid,mmom,emom)
      use ErrorHandling,        only : ErrorHandling_ERROR
      !.. Implicit declarations
      implicit none

      integer, intent(in) :: Natom           !< Number of atoms in system
      character(len=8), intent(in) :: simid  !< Name of simulation
      real(dblprec), dimension(Natom), intent(in) :: mmom   !< Magnitude of magnetic moments
      real(dblprec), dimension(3,Natom), intent(in) :: emom !< Current unit moment vector

      integer :: ii, jj
      character(len=35) :: filn

      filn = 'moment_sp.'//trim(adjustl(simid))//'.out'

      open(ofileno, file=filn, access = 'sequential',action = 'write', status = 'replace')
      ii=1
      do jj=1, Natom
         write (ofileno,10002) ii, jj, emom(1,jj), emom(2,jj), emom(3,jj), mmom(jj)
      end do
      close(ofileno)
      return
      call ErrorHandling_ERROR("Error writing path to file!")
      10002 format (i8,2x,i8,2x,2x, es16.8E3,2x,es16.8E3,2x,es16.8E3,2x,es16.8E3)

   end subroutine save_sp

end module PathInit
