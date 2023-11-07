!------------------------------------------------------------------------------------
!> @brief
!> Routines used to calculate topological properties of the magnetic system
!> Mostly related to the calculation of the skyrmion number
!
!> @author
!> Anders Bergman
!> Jonathan Chico ---> Reorganized in different modules, added site and type dependance
!> @copyright
!> GNU Public License
!------------------------------------------------------------------------------------
module Topology

   use Parameters
   use Profiling

   ! Parameters for the printing
   integer :: skyno_step !< Interval for sampling the skyrmion number
   integer :: skyno_buff !< Buffer size for the sampling of the skyrmion number
   character(len=1) :: skyno      !< Perform skyrmion number measurement
   character(len=1) :: do_proj_skyno !< Perform type dependent skyrmion number measurement
   character(len=1) :: do_skyno_den  !< Perform site dependent skyrmion number measurement
   character(len=1) :: do_skyno_cmass  !< Perform center-of-mass skyrmion number measurement

   integer :: nsimp !< Number of simplices
   integer, dimension(:,:), allocatable :: simp !< Array for storing Delaunay simplices

   public

contains

   !---------------------------------------------------------------------------------
   !> @brief
   !> Calculates the total skyrmion number of the system
   !
   !> @author
   !> Anders Bergman
   !---------------------------------------------------------------------------------
   real(dblprec) function pontryagin_no(Natom,Mensemble,emomM,grad_mom)
      use Constants

      implicit none

      integer, intent(in) :: Natom !< Number of atoms in system
      integer, intent(in) :: Mensemble !< Number of ensembles
      real(dblprec), dimension(3,Natom, Mensemble), intent(in) :: emomM  !< Current magnetic moment vector
      real(dblprec), dimension(3,3,Natom, Mensemble), intent(in) :: grad_mom  !< Gradient of magnetic moment vector

      integer :: iatom, k
      real(dblprec) :: thesum,cvec_x,cvec_y,cvec_z

      thesum=0.0_dblprec

      !$omp parallel do default(shared) private(iatom,k,cvec_x,cvec_y,cvec_z) reduction(+:thesum)
      do iatom=1, Natom
         do k=1, Mensemble
            cvec_x=grad_mom(2,1,iatom,k)*grad_mom(3,2,iatom,k)-grad_mom(3,1,iatom,k)*grad_mom(2,2,iatom,k)
            cvec_y=grad_mom(3,1,iatom,k)*grad_mom(1,2,iatom,k)-grad_mom(1,1,iatom,k)*grad_mom(3,2,iatom,k)
            cvec_z=grad_mom(1,1,iatom,k)*grad_mom(2,2,iatom,k)-grad_mom(2,1,iatom,k)*grad_mom(1,2,iatom,k)
            thesum=thesum+emomM(1,iatom,k)*cvec_x+emomM(2,iatom,k)*cvec_y+emomM(3,iatom,k)*cvec_z
         end do
      end do
      !$omp end parallel do

      pontryagin_no=thesum/pi/Mensemble
      !
      return
      !
   end function pontryagin_no

   !---------------------------------------------------------------------------------
   !> @brief
   !> Calculates the total skyrmion number of the system using triangulation
   !
   !> @author
   !> Anders Bergman
   !---------------------------------------------------------------------------------
   real(dblprec) function pontryagin_tri(Natom,Mensemble,emom)
      use Constants
      use math_functions

      implicit none

      integer, intent(in) :: Natom !< Number of atoms in system
      integer, intent(in) :: Mensemble !< Number of ensembles
      real(dblprec), dimension(3,Natom, Mensemble), intent(in) :: emom  !< Current magnetic moment vector


      integer :: k,i1,i2,i3, isimp
      real(dblprec), dimension(3) :: m1, m2, m3
      real(dblprec) :: thesum,q,qq, m1m2m3, m1m2, m1m3, m2m3

      thesum=0.0_dblprec

      !$omp parallel do default(shared) private(isimp,k,q,qq,m1,m2,m3,m1m2m3,m1m2,m1m3,m2m3) reduction(+:thesum)
      do isimp=1,nsimp
         do k=1, Mensemble
            m1 = emom(:,simp(1,isimp),k)
            m2 = emom(:,simp(2,isimp),k)
            m3 = emom(:,simp(3,isimp),k)
            m1m2m3=f_volume(m1,m2,m3)
            m1m2=dot_product(m1,m2)
            m1m3=dot_product(m1,m3)
            m2m3=dot_product(m2,m3)
            qq=m1m2m3/(1.0_dblprec+m1m2+m1m3+m2m3)
            q=2.0_dblprec*atan(qq)
            thesum=thesum+q
         end do
      end do
      !$omp end parallel do

      pontryagin_tri=thesum/(4.0_dblprec*pi)/Mensemble
      !
      return
      !
   end function pontryagin_tri

   !---------------------------------------------------------------------------------
   !> @brief
   !> Calculates projected skyrmion numbers of the system using triangulation
   !
   !> @author
   !> Anders Bergman
   !---------------------------------------------------------------------------------
function pontryagin_tri_proj(NA, Natom,Mensemble,emom)
      use Constants
      use math_functions

      implicit none

      integer, intent(in) :: NA    !< Number of atoms in unit cell
      integer, intent(in) :: Natom !< Number of atoms in system
      integer, intent(in) :: Mensemble !< Number of ensembles
      real(dblprec), dimension(3,Natom, Mensemble), intent(in) :: emom  !< Current magnetic moment vector


      integer :: k,i1,i2,i3, isimp, isite
      real(dblprec), dimension(3) :: m1, m2, m3
      real(dblprec) :: thesum,q,qq, m1m2m3, m1m2, m1m3, m2m3

      real(dblprec), dimension(NA) :: pontryagin_tri_proj 
      real(dblprec), dimension(NA) :: thesum_proj

      thesum_proj=0.0_dblprec

      !!$omp parallel do default(shared) private(isimp,k,q,qq,m1,m2,m3,m1m2m3,m1m2,m1m3,m2m3) reduction(+:thesum)
      do k=1, Mensemble
         do isite=1,NA
            do isimp=isite,nsimp+isite-1, NA
               m1 = emom(:,simp(1,isimp),k)
               m2 = emom(:,simp(2,isimp),k)
               m3 = emom(:,simp(3,isimp),k)
               m1m2m3=f_volume(m1,m2,m3)
               m1m2=dot_product(m1,m2)
               m1m3=dot_product(m1,m3)
               m2m3=dot_product(m2,m3)
               qq=m1m2m3/(1.0_dblprec+m1m2+m1m3+m2m3)
               q=2.0_dblprec*atan(qq)
               thesum_proj(isite)=thesum_proj(isite)+q
            end do
         end do
      end do
      !!$omp end parallel do

      pontryagin_tri_proj=thesum_proj/(4.0_dblprec*pi)/Mensemble
      !
      return
      !
   end function pontryagin_tri_proj
   !---------------------------------------------------------------------------------
   !> @brief
   !> Calculates the local skyrmion number density of the system using triangulation
   !
   !> @author
   !> Anders Bergman
   !---------------------------------------------------------------------------------
   real(dblprec) function pontryagin_tri_dens(iatom,Natom,Mensemble,emom)
      use Constants
      use math_functions

      implicit none

      integer, intent(in) :: Natom !< Number of atoms in system
      integer, intent(in) :: Mensemble !< Number of ensembles
      real(dblprec), dimension(3,Natom, Mensemble), intent(in) :: emom  !< Current magnetic moment vector
      integer, intent(in) :: iatom !< Current atom


      integer :: k,i1,i2,i3, isimp
      real(dblprec), dimension(3) :: m1, m2, m3
      real(dblprec) :: thesum,q,qq, m1m2m3, m1m2, m1m3, m2m3

      thesum=0.0_dblprec

      !isimp=2*iatom.

      do k=1, Mensemble
         do isimp=2*iatom-1,2*iatom
            m1 = emom(:,simp(1,isimp),k)
            m2 = emom(:,simp(2,isimp),k)
            m3 = emom(:,simp(3,isimp),k)
            m1m2m3=f_volume(m1,m2,m3)
            m1m2=dot_product(m1,m2)
            m1m3=dot_product(m1,m3)
            m2m3=dot_product(m2,m3)
            qq=m1m2m3/(1.0_dblprec+m1m2+m1m3+m2m3)
            q=2.0_dblprec*atan(qq)
            thesum=thesum+q
         end do
      end do

      pontryagin_tri_dens=thesum/(4.0_dblprec*pi)/Mensemble
      !
      return
      !
   end function pontryagin_tri_dens

   !---------------------------------------------------------------------------------
   !> @brief
   !> Calculation of the site dependent skyrmion number
   !
   !> @author
   !> Jonathan Chico
   !---------------------------------------------------------------------------------
   real(dblprec) function pontryagin_no_density(iatom,Natom,Mensemble,emomM,grad_mom)

      use Constants

      implicit none

      !.. Input variables
      integer, intent(in) :: iatom !< Current atomic position
      integer, intent(in) :: Natom !< Number of atoms in system
      integer, intent(in) :: Mensemble !< Number of ensembles
      real(dblprec), dimension(3,Natom, Mensemble), intent(in) :: emomM  !< Current magnetic moment vector
      real(dblprec), dimension(3,3,Natom, Mensemble), intent(in) :: grad_mom  !< Gradient of magnetic moment vector

      ! .. Local variables
      integer :: k
      real(dblprec) :: thesum,cvec_x,cvec_y,cvec_z

      thesum=0.0_dblprec

      do k=1, Mensemble
         cvec_x=grad_mom(2,1,iatom,k)*grad_mom(3,2,iatom,k)-grad_mom(3,1,iatom,k)*grad_mom(2,2,iatom,k)
         cvec_y=grad_mom(3,1,iatom,k)*grad_mom(1,2,iatom,k)-grad_mom(1,1,iatom,k)*grad_mom(3,2,iatom,k)
         cvec_z=grad_mom(1,1,iatom,k)*grad_mom(2,2,iatom,k)-grad_mom(2,1,iatom,k)*grad_mom(1,2,iatom,k)
         thesum=thesum+emomM(1,iatom,k)*cvec_x+emomM(2,iatom,k)*cvec_y+emomM(3,iatom,k)*cvec_z
      end do

      pontryagin_no_density=thesum/pi/Mensemble

   end function pontryagin_no_density

   !---------------------------------------------------------------------------------
   !> @brief
   !> Calculation of the type dependent skyrmion number
   !
   !> @author
   !> Jonathan Chico
   !---------------------------------------------------------------------------------
   function proj_pontryagin_no(NT,Natom,Mensemble,atype,emomM,proj_grad_mom)
      use Constants

      implicit none

      integer, intent(in) :: NT    !< Number of types of atoms
      integer, intent(in) :: Natom !< Number of atoms in system
      integer, intent(in) :: Mensemble !< Number of ensembles
      integer, dimension(Natom), intent(in) :: atype !< Type of atom
      real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: emomM  !< Current magnetic moment vector
      real(dblprec), dimension(3,3,Natom,Mensemble,NT), intent(in) :: proj_grad_mom  !< Gradient of magnetic moment vector

      real(dblprec), dimension(NT) :: proj_pontryagin_no

      integer :: iatom, k,ii
      real(dblprec), dimension(NT) :: thesum,cvec_x,cvec_y,cvec_z

      thesum=0.0_dblprec

      !$omp parallel do default(shared) private(iatom,k,cvec_x,cvec_y,cvec_z,ii) reduction(+:thesum)
      do iatom=1, Natom
         do k=1, Mensemble
            ii=atype(iatom)
            cvec_x(ii)=proj_grad_mom(2,1,iatom,k,ii)*proj_grad_mom(3,2,iatom,k,ii)-proj_grad_mom(3,1,iatom,k,ii)*proj_grad_mom(2,2,iatom,k,ii)
            cvec_y(ii)=proj_grad_mom(3,1,iatom,k,ii)*proj_grad_mom(1,2,iatom,k,ii)-proj_grad_mom(1,1,iatom,k,ii)*proj_grad_mom(3,2,iatom,k,ii)
            cvec_z(ii)=proj_grad_mom(1,1,iatom,k,ii)*proj_grad_mom(2,2,iatom,k,ii)-proj_grad_mom(2,1,iatom,k,ii)*proj_grad_mom(1,2,iatom,k,ii)
            thesum(ii)=thesum(ii)+emomM(1,iatom,k)*cvec_x(ii)+emomM(2,iatom,k)*cvec_y(ii)+emomM(3,iatom,k)*cvec_z(ii)
         end do
      end do
      !$omp end parallel do

      proj_pontryagin_no=thesum/pi/Mensemble
      !
      return
      !
   end function proj_pontryagin_no

   !---------------------------------------------------------------------------------
   !> @brief
   !> Constructing a Delaunay triangulation from an a priori known triangular
   !lattice
   !
   !> @author
   !> Anders Bergman
   !---------------------------------------------------------------------------------
   subroutine delaunay_tri_tri(nx,ny,nz, NT)
      use Constants

      implicit none

      integer, intent(in) :: nx    !< Unit cell repetitions in X
      integer, intent(in) :: ny    !< Unit cell repetitions in Y
      integer, intent(in) :: nz    !< Unit cell repetitions in Z
      integer, intent(in) :: NT    !< Number of types of atoms


      integer :: x,y,z, ix,iy,iz, it

      !!! if (NT>1) then
      !!!    write(*,'(1x,a)') 'Hard-coded triangulation only works for one atom per unit cell.' 
      !!!    stop
      !!! endif

      nsimp = 2*nx*ny*nz*nt
      allocate(simp(3,nsimp))
      simp=1

      nsimp=0
      do z=1,nz
         do y=1,ny
            do x=1,nx
               do it=1,NT
                  ! Triangle 1  [0 -> +x -> +y -> 0]
                  nsimp=nsimp+1
                  simp(1,nsimp)=NT*wrap_idx(x,y,z,nx,ny,nz)+it-NT
                  simp(2,nsimp)=NT*wrap_idx(modulo(x,nx)+1,y,z,nx,ny,nz)+it-NT
                  simp(3,nsimp)=NT*wrap_idx(x,modulo(y,ny)+1,z,nx,ny,nz)+it-NT
               end do
            end do
         end do
      end do
      do z=1,nz
         do y=1,ny
            do x=1,nx
               do it=1,NT
                  ! Triangle  2  [ 0 -> +y -> +y - x -> 0]
                  nsimp=nsimp+1
                  simp(1,nsimp)=NT*wrap_idx(x,y,z,nx,ny,nz)+it-NT
                  simp(2,nsimp)=NT*wrap_idx(x,modulo(y,ny)+1,z,nx,ny,nz)+it-NT
                  simp(3,nsimp)=NT*wrap_idx(modulo(x-2,nx)+1,modulo(y,ny)+1,z,nx,ny,nz)+it-NT
               end do
            end do
         end do
      end do


      if (nsimp .ne. 2*(nx)*(ny)*nz*nt) then
         write(*,'(1x,a,2i8)') 'Warning: Delaunay failure', nsimp, (2*nx-1)*(ny-1)*nz
         print *,shape(simp)
      end if


      !
      return

      contains 

         integer function wrap_idx(x,y,z,nx,ny,nz)
            !
            implicit none
            !
            integer, intent(in) :: x,y,z
            integer, intent(in) :: nx,ny,nz

               wrap_idx=nx*ny*(z-1)+nx*(y-1)+x

            end function wrap_idx
      !
   end subroutine delaunay_tri_tri

end module Topology
