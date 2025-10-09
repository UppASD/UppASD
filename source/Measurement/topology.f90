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
   use Systemdata, only : coord

   ! Parameters for the printing
   integer :: skyno_step !< Interval for sampling the skyrmion number
   integer :: skyno_buff !< Buffer size for the sampling of the skyrmion number
   character(len=1) :: skyno      !< Perform skyrmion number measurement
   character(len=1) :: do_proj_skyno !< Perform type dependent skyrmion number measurement
   character(len=1) :: do_skyno_den  !< Perform site dependent skyrmion number measurement
   character(len=1) :: do_skyno_cmass  !< Perform center-of-mass skyrmion number measurement

   integer :: nsimp !< Number of simplices
   integer, dimension(:,:), allocatable :: simp !< Array for storing Delaunay simplices

   real(dblprec) :: chi_avg !< Average scalar chirality (instantaneous)
   real(dblprec) :: chi_cavg  = 0.0_dblprec !< Average scalar chirality (cumulative)
   real(dblprec), dimension(3) :: kappa_cavg = 0.0_dblprec !< Average scalar chirality (cumulative)
   real(dblprec), dimension(3) :: kappa_csum = 0.0_dblprec !< Cumulative sum of the vector chirality
   integer :: n_chi_cavg = 0 !< Number of times the average scalar chirality has been calculated

   ! Variables for OAM 
   character(len=1) :: do_oam = 'N' !< Perform OAM measurement
   integer :: oam_step !< Interval for sampling OAM
   integer :: oam_buff !< Buffer size for the sampling of OAM
   character(len=1) :: print_mesh = 'N' !< Print triangulation mesh to file
   real(dblprec), dimension(:, :, :), allocatable :: mu_arr !< Array for dynamic magnetization
   real(dblprec), dimension(:, :), allocatable :: S0_arr !< Array for static magnetization
   real(dblprec), dimension(:), allocatable:: Lz_t !< Array for OAM
   real(dblprec), dimension(:, :), allocatable:: Lz_i !< Array for local OAM
   real(dblprec) :: Lz_avg !< Average OAM (instantaneous)
   real(dblprec) :: Lz_tot !< Total OAM (instantaneous)
   real(dblprec) :: Lz_csum  = 0.0_dblprec !< Cumulative sum of the OAM
   integer :: n_Lz_cavg = 0 !< Number of times the average OAM has been calculated

   integer :: step_counter = 0

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
subroutine delaunay_tri_tri(nx,ny,nz,NT,coords)
   use Constants
   implicit none

   integer, intent(in) :: nx,ny,nz,NT
   real(dblprec), intent(in) :: coords(3,nx*ny*nz*NT)

   integer :: x,y,z,it, nsimp_max
   integer :: i00,i10,i01,i11
   real(dblprec) :: d2_diag1,d2_diag2
   real(dblprec), dimension(3) :: r00,r10,r01,r11

   nsimp_max = 2*nx*ny*nz*NT
   allocate(simp(3,nsimp_max))
   nsimp = 0

   do z=1,nz
      do y=1,ny
         do x=1,nx
            do it=1,NT
               ! FIXED: Proper indexing calculation
               ! Each cell has NT atoms, wrap_idx gives 1-based cell index
               i00 = NT*(wrap_idx(x,       y,       z,nx,ny,nz)-1) + it
               i10 = NT*(wrap_idx(modulo(x,nx)+1, y,       z,nx,ny,nz)-1) + it
               i01 = NT*(wrap_idx(x,       modulo(y,ny)+1,z,nx,ny,nz)-1) + it
               i11 = NT*(wrap_idx(modulo(x,nx)+1, modulo(y,ny)+1,z,nx,ny,nz)-1) + it

               r00 = coords(:,i00)
               r10 = coords(:,i10)
               r01 = coords(:,i01)
               r11 = coords(:,i11)

               ! Compare diagonals to choose triangulation
               d2_diag1 = sum((r00-r11)**2)
               d2_diag2 = sum((r10-r01)**2)

               if (d2_diag1 <= d2_diag2) then
                  nsimp = nsimp+1 ; simp(:,nsimp) = [i00,i10,i11]
                  nsimp = nsimp+1 ; simp(:,nsimp) = [i00,i11,i01]
               else
                  nsimp = nsimp+1 ; simp(:,nsimp) = [i00,i10,i01]
                  nsimp = nsimp+1 ; simp(:,nsimp) = [i10,i11,i01]
               end if

            end do
         end do
      end do
   end do

   write(*,'(1x,a,i8,a)') 'Triangulation created with ', nsimp, ' simplices'

   contains
      integer function wrap_idx(x,y,z,nx,ny,nz)
         implicit none
         integer, intent(in) :: x,y,z,nx,ny,nz
         wrap_idx = nx*ny*(z-1) + nx*(y-1) + x
      end function wrap_idx
end subroutine delaunay_tri_tri


   !---------------------------------------------------------------------------------
   !> @brief
   !> Calculates the scalar chirality using triangulation
   !
   !> @author
   !> Anders Bergman
   !---------------------------------------------------------------------------------
   function chirality_tri(Natom,Mensemble,emom) result(kappa_avg)
      use constants
      use math_functions, only : f_cross_product
      use InputData, only: N1, N2, N3, NA
      implicit none
      integer, intent(in) :: Natom, Mensemble
      real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: emom

      real(dblprec), dimension(3) :: kappa_tot   ! global vector chirality
      real(dblprec), dimension(3) :: kappa_avg   ! average vector chirality
      real(dblprec) :: chi_tot                   ! optional scalar version
      real(dblprec), dimension(3) :: m1,m2,m3
      real(dblprec), dimension(3) :: c12,c23,c31
      integer :: k,isimp

      ! Ensure triangulation is set up
      if (nsimp == 0) then
         write(*,'(1x, a)') "Setting up triangulation for chirality calculation"
         call delaunay_tri_tri(N1, N2, N3, NA, coord)
      end if

      kappa_tot = 0.0_dblprec
      chi_tot   = 0.0_dblprec

      !$omp parallel do default(shared) private(isimp,k,m1,m2,m3,c12,c23,c31)  &
      !$omp& reduction(+:kappa_tot,chi_tot)
      do isimp = 1, nsimp
         do k = 1, Mensemble
            m1 = emom(:,simp(1,isimp),k)
            m2 = emom(:,simp(2,isimp),k)
            m3 = emom(:,simp(3,isimp),k)

            ! pairwise cross-products
            c12 = f_cross_product(m1,m2)   ! c12 = m1 × m2
            c23 = f_cross_product(m2,m3)
            c31 = f_cross_product(m3,m1)

            ! accumulate vector chirality for this triangle
            kappa_tot = kappa_tot + c12 + c23 + c31

            ! optional: scalar chirality
            chi_tot = chi_tot + dot_product(m1,c23)   ! m1·(m2×m3)
         end do
      end do
      !$omp end parallel do

      kappa_avg = kappa_tot / real(nsimp*Mensemble, dblprec)
      ! Chi avg stored in module data
      chi_avg = chi_tot / real(nsimp*Mensemble, dblprec)
   end function chirality_tri

!===============================================================
!   Main driver routine  (flag = 0/1/2)
!===============================================================
   subroutine calculate_oam(Natom,Mensemble,emom,mstep,flag)
      use math_functions, only : f_cross_product
      use SystemData, only : coord
      use InputData, only : simid, N1, N2, N3, NA, C1, C2, C3
      implicit none
      integer,          intent(in) :: Natom, Mensemble, mstep, flag
      real(dblprec),    intent(in) :: emom(3,Natom,Mensemble)
   
      integer :: i,k,t,v1,v2,v3
      real(dblprec):: n_hat(3), cross_mu_n(3)
      real(dblprec):: dmu_dx(3), dmu_dy(3)
      !
      character(len=30) :: filn
   
   !----------------------------------------------------------------
   select case(flag)
   
   !------------------ 0 :  allocate & initialise ------------------
   case (0)
   
      ! Ensure triangulation is set up for solid-angle OAM calculation
      if (nsimp == 0) then
         write(*,'(1x, a)') "Setting up triangulation for OAM calculation"
         call delaunay_tri_tri(N1, N2, N3, NA, coord)
      end if

      ! Allocate arrays (keeping some for compatibility, but not all are needed for triangulation approach)
      allocate(S0_arr(3,Natom), mu_arr(3,Natom,Mensemble))
      allocate(Lz_i(Natom,Mensemble))
      allocate(Lz_t(0:mstep))
   
      S0_arr = emom(:,:,1)               ! ground state from first frame
      step_counter = 0
   
   !------------------ 1 :  sample current step --------------------
   case (1)
      if ( mod(mstep-1,oam_step)==0) then
      step_counter = step_counter + 1

      ! Use solid-angle triangulation approach for OAM calculation
      ! Lz_tot = oam_tri_phase(Natom, Mensemble, emom, coord)
      Lz_tot = oam_tri_improved(Natom, Mensemble, emom, coord)
      ! Lz_tot = oam_tri(Natom, Mensemble, emom, coord)
      Lz_avg = Lz_tot / real(Natom, dblprec)

      n_Lz_cavg = n_Lz_cavg + 1
      Lz_csum = Lz_csum + Lz_avg

      ! write the OAM data
      write(filn,'(''oam.'',a,''.out'')') trim(simid)
      if (mstep <= 1) then
         open(ofileno,file=filn, status="replace")
         write (ofileno,*) "#step      Lz_tot          Lz_avg         Lz_cavg"
         close(ofileno)
      end if
      open(ofileno,file=filn, position="append")
      write (ofileno,240) mstep, Lz_tot, Lz_avg, Lz_csum/n_Lz_cavg
      close(ofileno)

   end if
   
   !------------------ 2 :  FFT / output / free --------------------
   case (2)
      !call fft_and_print(step_counter,Lz_t)   ! (replace with real FFT)
   
      deallocate(S0_arr, mu_arr, Lz_i, Lz_t)
      step_counter = 0
   end select
   !----------------------------------------------------------------

   240 format(i8,2x,3f16.8)
   end subroutine calculate_oam
  
!===============================================================
!> Compute orbital angular momentum Lz from triangulated phase windings
!> - Uses transverse fluctuations relative to ground state S0_arr
!> - Works for arbitrary noncollinear reference states
!===============================================================
real(dblprec) function oam_tri(Natom, Mensemble, emom, coords) result(Lz)
   use Constants
   use InputData, only: C1,C2,C3,N1,N2,N3
   implicit none
   integer, intent(in) :: Natom, Mensemble
   real(dblprec), intent(in) :: emom(3,Natom,Mensemble)
   real(dblprec), intent(in) :: coords(3,Natom)

   integer :: isimp, k, i1,i2,i3
   real(dblprec) :: Lx,Ly,Lz_box
   real(dblprec) :: Lsum, rho_t, gamma_t, area_t
   real(dblprec) :: x1,y1,x2,y2,x3,y3, cross
   real(dblprec) :: th1,th2,th3,dth12,dth23,dth31
   real(dblprec) :: mu1(2),mu2(2),mu3(2)
   real(dblprec) :: ex(3),ey(3),ez(3)
   real(dblprec) :: Mavg(3),norm

   !--- box lengths
   Lx = N1*sqrt(dot_product(C1,C1))
   Ly = N2*sqrt(dot_product(C2,C2))
   Lz_box = N3*sqrt(dot_product(C3,C3))

   !--- choose a global frame from average magnetization of ensemble[1]
   Mavg = 0.0_dblprec
   do i1=1,Natom
      Mavg = Mavg + emom(:,i1,1)
   end do
   Mavg = Mavg/sqrt(dot_product(Mavg,Mavg))

   call make_orthonormal_basis(Mavg,ex,ey,ez)  ! ez=Mavg, ex,ey span transverse plane

   Lsum = 0.0_dblprec

   !$omp parallel do default(shared) private(isimp,k,i1,i2,i3,mu1,mu2,mu3,th1,th2,th3, &
   !$omp dth12,dth23,dth31,gamma_t,rho_t,x1,y1,x2,y2,x3,y3,cross,area_t,norm) reduction(+:Lsum)
   do isimp=1,nsimp
      i1 = simp(1,isimp); i2 = simp(2,isimp); i3 = simp(3,isimp)

      do k=1,Mensemble
         ! project onto fixed (ex,ey)
         mu1 = [ dot_product(emom(:,i1,k),ex), dot_product(emom(:,i1,k),ey) ]
         mu2 = [ dot_product(emom(:,i2,k),ex), dot_product(emom(:,i2,k),ey) ]
         mu3 = [ dot_product(emom(:,i3,k),ex), dot_product(emom(:,i3,k),ey) ]

         th1 = atan2(mu1(2),mu1(1))
         th2 = atan2(mu2(2),mu2(1))
         th3 = atan2(mu3(2),mu3(1))

         ! unwrap phase differences
         dth12 = modulo(th2-th1+pi,2*pi)-pi
         dth23 = modulo(th3-th2+pi,2*pi)-pi
         dth31 = modulo(th1-th3+pi,2*pi)-pi
         gamma_t = dth12 + dth23 + dth31

         ! amplitude (mean squared norm)
         rho_t = ( dot_product(mu1,mu1) + dot_product(mu2,mu2) + dot_product(mu3,mu3) ) / 3.0_dblprec

         ! area (2D, using coords projected to xy plane)
         x1=coords(1,i1); y1=coords(2,i1)
         x2=coords(1,i2); y2=coords(2,i2)
         x3=coords(1,i3); y3=coords(2,i3)
         cross = (x2-x1)*(y3-y1) - (y2-y1)*(x3-x1)
         area_t = 0.5_dblprec*abs(cross)

         Lsum = Lsum + rho_t * gamma_t * area_t
      end do
   end do
   !$omp end parallel do

   Lz = Lsum/Mensemble
end function oam_tri

!===============================================================
!> @brief
!> IMPROVED: Compute OAM using solid-angle formulation
!> This is more robust for noncollinear spin configurations
!===============================================================
real(dblprec) function oam_tri_improved(Natom, Mensemble, emom, coords) result(Lz)
   use Constants
   use InputData, only: C1,C2,C3,N1,N2,N3
   use math_functions, only: f_volume
   implicit none
   integer, intent(in) :: Natom, Mensemble
   real(dblprec), intent(in) :: emom(3,Natom,Mensemble)
   real(dblprec), intent(in) :: coords(3,Natom)

   integer :: isimp, k, i1,i2,i3
   real(dblprec) :: Lsum, Omega_t, area_t, rho_t
   real(dblprec) :: x1,y1,x2,y2,x3,y3, cross
   real(dblprec) :: m1(3),m2(3),m3(3), dm1(3),dm2(3),dm3(3)
   real(dblprec) :: m1m2m3, m1m2, m1m3, m2m3, qq
   real(dblprec) :: Lx, Ly

   ! System dimensions
   Lx = N1*sqrt(dot_product(C1,C1))
   Ly = N2*sqrt(dot_product(C2,C2))

   Lsum = 0.0_dblprec

   !$omp parallel do default(shared) private(isimp,k,i1,i2,i3,m1,m2,m3,dm1,dm2,dm3, &
   !$omp m1m2m3,m1m2,m1m3,m2m3,qq,Omega_t,rho_t,x1,y1,x2,y2,x3,y3,cross,area_t) &
   !$omp reduction(+:Lsum)
   do isimp=1,nsimp
      i1 = simp(1,isimp); i2 = simp(2,isimp); i3 = simp(3,isimp)

      ! Get triangle coordinates (xy projection)
      x1=coords(1,i1); y1=coords(2,i1)
      x2=coords(1,i2); y2=coords(2,i2)
      x3=coords(1,i3); y3=coords(2,i3)

      ! Triangle area (2D cross product)
      cross = (x2-x1)*(y3-y1) - (y2-y1)*(x3-x1)
      area_t = 0.5_dblprec*abs(cross)

      do k=1,Mensemble
         ! Current spin configuration
         m1 = emom(:,i1,k)
         m2 = emom(:,i2,k)
         m3 = emom(:,i3,k)

         ! Deviations from ground state (if tracking dynamics)
         dm1 = m1 - S0_arr(:,i1)
         dm2 = m2 - S0_arr(:,i2)
         dm3 = m3 - S0_arr(:,i3)

         ! Average transverse fluctuation amplitude
         rho_t = (dot_product(dm1,dm1) + dot_product(dm2,dm2) + dot_product(dm3,dm3))/3.0_dblprec

         ! Solid angle (skyrmion number formula)
         m1m2m3 = f_volume(m1,m2,m3)
         m1m2 = dot_product(m1,m2)
         m1m3 = dot_product(m1,m3)
         m2m3 = dot_product(m2,m3)
         qq = m1m2m3/(1.0_dblprec+m1m2+m1m3+m2m3)
         Omega_t = 2.0_dblprec*atan(qq)

         ! OAM contribution: ρ_t * Ω_t * A_t
         Lsum = Lsum + rho_t * Omega_t * area_t
      end do
   end do
   !$omp end parallel do

   Lz = Lsum / Mensemble
end function oam_tri_improved

!===============================================================
!> @brief  
!> ALTERNATIVE: Phase-winding OAM with proper local frames
!> Uses local coordinate system at each vertex for better accuracy
!===============================================================
real(dblprec) function oam_tri_phase(Natom, Mensemble, emom, coords) result(Lz)
   use Constants
   use InputData, only: C1,C2,C3,N1,N2,N3
   implicit none
   integer, intent(in) :: Natom, Mensemble
   real(dblprec), intent(in) :: emom(3,Natom,Mensemble)
   real(dblprec), intent(in) :: coords(3,Natom)

   integer :: isimp, k, i1,i2,i3
   real(dblprec) :: Lsum, rho_t, gamma_t, area_t
   real(dblprec) :: x1,y1,x2,y2,x3,y3, cross
   real(dblprec) :: th1,th2,th3,dth12,dth23,dth31
   real(dblprec) :: mu1(2),mu2(2),mu3(2)
   real(dblprec) :: ex(3),ey(3),ez(3)
   real(dblprec) :: dm1(3),dm2(3),dm3(3)

   Lsum = 0.0_dblprec

   !$omp parallel do default(shared) private(isimp,k,i1,i2,i3,dm1,dm2,dm3, &
   !$omp ex,ey,ez,mu1,mu2,mu3,th1,th2,th3,dth12,dth23,dth31,gamma_t,rho_t, &
   !$omp x1,y1,x2,y2,x3,y3,cross,area_t) reduction(+:Lsum)
   do isimp=1,nsimp
      i1 = simp(1,isimp); i2 = simp(2,isimp); i3 = simp(3,isimp)

      ! Triangle area
      x1=coords(1,i1); y1=coords(2,i1)
      x2=coords(1,i2); y2=coords(2,i2)
      x3=coords(1,i3); y3=coords(2,i3)
      cross = (x2-x1)*(y3-y1) - (y2-y1)*(x3-x1)
      area_t = 0.5_dblprec*abs(cross)

      do k=1,Mensemble
         ! Deviations from ground state
         dm1 = emom(:,i1,k) - S0_arr(:,i1)
         dm2 = emom(:,i2,k) - S0_arr(:,i2)
         dm3 = emom(:,i3,k) - S0_arr(:,i3)

         ! Use LOCAL frame at each vertex for projection
         call local_frame(S0_arr(:,i1), ex,ey,ez)
         mu1 = [dot_product(dm1,ex), dot_product(dm1,ey)]
         th1 = atan2(mu1(2),mu1(1))

         call local_frame(S0_arr(:,i2), ex,ey,ez)
         mu2 = [dot_product(dm2,ex), dot_product(dm2,ey)]
         th2 = atan2(mu2(2),mu2(1))

         call local_frame(S0_arr(:,i3), ex,ey,ez)
         mu3 = [dot_product(dm3,ex), dot_product(dm3,ey)]
         th3 = atan2(mu3(2),mu3(1))

         ! Average amplitude
         rho_t = (dot_product(mu1,mu1) + dot_product(mu2,mu2) + dot_product(mu3,mu3))/3.0_dblprec

         ! Unwrap phase differences
         dth12 = modulo(th2-th1+pi,2*pi)-pi
         dth23 = modulo(th3-th2+pi,2*pi)-pi
         dth31 = modulo(th1-th3+pi,2*pi)-pi
         gamma_t = dth12 + dth23 + dth31

         Lsum = Lsum + rho_t * gamma_t * area_t
      end do
   end do
   !$omp end parallel do

   Lz = Lsum / Mensemble
end function oam_tri_phase

subroutine make_orthonormal_basis(nz,ex,ey,ez)
   real(dblprec), intent(in)  :: nz(3)
   real(dblprec), intent(out) :: ex(3),ey(3),ez(3)
   real(dblprec) :: tmp(3)
   ez = nz
   if (abs(ez(1))<0.9_dblprec) then
      tmp = [1.0_dblprec,0.0_dblprec,0.0_dblprec]
   else
      tmp = [0.0_dblprec,1.0_dblprec,0.0_dblprec]
   end if
   ex = tmp - dot_product(tmp,ez)*ez
   ex = ex/sqrt(dot_product(ex,ex))
   ey = [ ez(2)*ex(3)-ez(3)*ex(2), ez(3)*ex(1)-ez(1)*ex(3), ez(1)*ex(2)-ez(2)*ex(1) ]
end subroutine

! real(dblprec) function oam_tri(Natom, Mensemble, emom, coords)
!    use Constants
!    use InputData, only: C1, C2, C3, N1, N2, N3
!    implicit none
!    integer, intent(in) :: Natom
!    integer, intent(in) :: Mensemble
!    real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: emom
!    real(dblprec), dimension(3,Natom), intent(in) :: coords   ! atom positions (x,y,z)
!    ! real(dblprec), dimension(3,Natom), intent(in) :: S0_arr   ! reference spins at t=0
! 
!    integer :: k, isimp, i1,i2,i3
!    real(dblprec) :: Lsum, rho_t, gamma_t, area_t
!    real(dblprec) :: x1,y1,x2,y2,x3,y3, cross
!    real(dblprec) :: z1,z2,z3
!    real(dblprec) :: th1, th2, th3, dth12, dth23, dth31
!    real(dblprec) :: mu(2), norm
! 
!    ! Box dimensions for minimal image
!    real(dblprec) :: Lx, Ly, Lz
! 
!    ! Variables for triangle validity check
!    real(dblprec) :: dist2, dist3, threshold
! 
!    ! local frame basis
!    real(dblprec) :: ex(3), ey(3), ez(3), psi(3)
! 
!    ! Compute box dimensions
!    Lx = N1 * sqrt(dot_product(C1,C1))
!    Ly = N2 * sqrt(dot_product(C2,C2))
!    Lz = N3 * sqrt(dot_product(C3,C3))
! 
!    ! print *, 'OAM: System dimensions (Lx, Ly, Lz): ', Lx, Ly, Lz
! 
!    Lsum = 0.0_dblprec
! 
!    ex = 0.0_dblprec; ex(1) = 1.0_dblprec
!    ey = 0.0_dblprec; ey(2) = 1.0_dblprec
!    ez = 0.0_dblprec; ez(3) = 1.0_dblprec
!    !$omp parallel do default(shared) private(isimp,k,i1,i2,i3,th1,th2,th3,dth12,dth23,dth31,gamma_t,rho_t, &
!    !$omp   x1,y1,z1,x2,y2,z2,x3,y3,z3,cross,area_t,mu,norm,ex,ey,ez) reduction(+:Lsum)
!    do isimp=1,nsimp
!       i1 = simp(1,isimp); i2 = simp(2,isimp); i3 = simp(3,isimp)
!       do k=1,Mensemble
!          ! coordinates
!          x1=coords(1,i1); y1=coords(2,i1); z1=coords(3,i1)
!          x2=coords(1,i2); y2=coords(2,i2); z2=coords(3,i2)
!          x3=coords(1,i3); y3=coords(2,i3); z3=coords(3,i3)
! 
!          ! Apply minimal image correction relative to vertex 1 using fractional coordinates
!          ! call minimal_image_correction(x2, y2, z2, x1, y1, z1, C1, C2, C3, N1, N2, N3)
!          ! call minimal_image_correction(x3, y3, z3, x1, y1, z1, C1, C2, C3, N1, N2, N3)
! 
!          ! Check if triangle is valid (not spanning across periodic boundaries)
!          dist2 = sqrt((x2-x1)**2 + (y2-y1)**2 + (z2-z1)**2)
!          dist3 = sqrt((x3-x1)**2 + (y3-y1)**2 + (z3-z1)**2)
!          threshold = min(Lx, Ly) / 2.0_dblprec
!          if (dist2 > threshold .or. dist3 > threshold) cycle
! 
!          rho_t = 0.0_dblprec
! 
!          ! --- vertex 1 ---
!          call local_frame(S0_arr(:,i1), ex,ey,ez)
!          !psi = emom(:,i1,k) - S0_arr(:,i1)
!          mu = project_transverse(emom(:,i1,k)-S0_arr(:,i1), ex,ey,ez)
!          !print *, ' 1 mu = ', isimp, mu, atan2(mu(2),mu(1))
!          th1 = atan2(mu(2),mu(1))
!          norm = mu(1)*mu(1)+mu(2)*mu(2)
!          rho_t = rho_t + norm
! 
!          ! --- vertex 2 ---
!          call local_frame(S0_arr(:,i2), ex,ey,ez)
!          !psi = emom(:,i2,k) - S0_arr(:,i2)
!          mu = project_transverse(emom(:,i2,k)-S0_arr(:,i2), ex,ey,ez)
!          !print *, ' 2 mu = ', isimp, mu, atan2(mu(2),mu(1))
!          th2 = atan2(mu(2),mu(1))
!          norm = mu(1)*mu(1)+mu(2)*mu(2)
!          rho_t = rho_t + norm
! 
!          ! --- vertex 3 ---
!          call local_frame(S0_arr(:,i3), ex,ey,ez)
!          !psi = emom(:,i3,k) - S0_arr(:,i3)
!          mu = project_transverse(emom(:,i3,k)-S0_arr(:,i3), ex,ey,ez)
!          !print *, ' 3 mu = ', isimp, mu, atan2(mu(2),mu(1))
!          th3 = atan2(mu(2),mu(1))
!          norm = mu(1)*mu(1)+mu(2)*mu(2)
!          rho_t = rho_t + norm
! 
!          rho_t = rho_t/3.0_dblprec
! 
!          ! unwrap phase differences
!          dth12 = modulo(th2-th1+pi,2*pi)-pi
!          dth23 = modulo(th3-th2+pi,2*pi)-pi
!          dth31 = modulo(th1-th3+pi,2*pi)-pi
!          gamma_t = dth12 + dth23 + dth31
! 
!          ! triangle area (2D cross product)
!          cross = (x2-x1)*(y3-y1) - (y2-y1)*(x3-x1)
!          area_t = 0.5_dblprec*abs(cross)
! 
!          ! OAM contribution
!          Lsum = Lsum + rho_t * gamma_t * area_t
!       end do
!    end do
!    !$omp end parallel do
! 
!    oam_tri = Lsum / Mensemble
!    return
! end function oam_tri
! 
!===============================================================
!> Apply minimal image correction to coordinates relative to a reference point
!> using fractional coordinates for general lattices
!===============================================================
subroutine minimal_image_correction(x, y, z, x_ref, y_ref, z_ref, C1, C2, C3, N1, N2, N3)
   use Constants
   implicit none
   real(dblprec), intent(inout) :: x, y, z
   real(dblprec), intent(in) :: x_ref, y_ref, z_ref
   real(dblprec), dimension(3), intent(in) :: C1, C2, C3
   integer, intent(in) :: N1, N2, N3
   
   real(dblprec) :: dx, dy, dz, s1, s2, s3
   real(dblprec), dimension(3) :: cross23, cross31, cross12, r_diff
   real(dblprec) :: vol
   
   ! Vector difference
   dx = x - x_ref
   dy = y - y_ref  
   dz = z - z_ref
   r_diff = [dx, dy, dz]
   
   ! Calculate cross products and volume for fractional coordinate conversion
   cross23(1) = C2(2)*C3(3) - C2(3)*C3(2)
   cross23(2) = C2(3)*C3(1) - C2(1)*C3(3)
   cross23(3) = C2(1)*C3(2) - C2(2)*C3(1)
   vol = dot_product(C1, cross23)  ! C1 · (C2 × C3)
   
   cross31(1) = C3(2)*C1(3) - C3(3)*C1(2)
   cross31(2) = C3(3)*C1(1) - C3(1)*C1(3)
   cross31(3) = C3(1)*C1(2) - C3(2)*C1(1)
   
   cross12(1) = C1(2)*C2(3) - C1(3)*C2(2)
   cross12(2) = C1(3)*C2(1) - C1(1)*C2(3)
   cross12(3) = C1(1)*C2(2) - C1(2)*C2(1)
   
   ! Convert to fractional coordinates
   s1 = dot_product(r_diff, cross23) / vol
   s2 = dot_product(r_diff, cross31) / vol
   s3 = dot_product(r_diff, cross12) / vol
   
   ! Apply periodic boundary conditions for supercell
   s1 = s1 - real(N1, dblprec) * nint(s1 / real(N1, dblprec))
   s2 = s2 - real(N2, dblprec) * nint(s2 / real(N2, dblprec))
   s3 = s3 - real(N3, dblprec) * nint(s3 / real(N3, dblprec))
   
   ! Convert back to Cartesian coordinates relative to reference
   x = x_ref + s1*C1(1) + s2*C2(1) + s3*C3(1)
   y = y_ref + s1*C1(2) + s2*C2(2) + s3*C3(2)
   z = z_ref + s1*C1(3) + s2*C2(3) + s3*C3(3)
end subroutine minimal_image_correction

subroutine local_frame(n0, ex, ey, ez)
   use Constants
   implicit none
   real(dblprec), intent(in)  :: n0(3)
   real(dblprec), intent(out) :: ex(3), ey(3), ez(3)
   real(dblprec) :: tmp(3)

   ez = n0 / sqrt(dot_product(n0,n0))
   if (abs(ez(1)) < 0.9_dblprec) then
      tmp = [1.0_dblprec, 0.0_dblprec, 0.0_dblprec]
   else
      tmp = [0.0_dblprec, 1.0_dblprec, 0.0_dblprec]
   end if
   ex = tmp - dot_product(tmp,ez)*ez
   ex = ex / sqrt(dot_product(ex,ex))
   ey(1) = ez(2)*ex(3) - ez(3)*ex(2)
   ey(2) = ez(3)*ex(1) - ez(1)*ex(3)
   ey(3) = ez(1)*ex(2) - ez(2)*ex(1)
end subroutine local_frame

!===============================================================
!> Project a fluctuation vector mu_vec into the transverse plane
!> defined by (ex,ey)
!===============================================================
function project_transverse(mu_vec, ex, ey, ez) result(mu)
   use Constants
   implicit none
   real(dblprec), intent(in) :: mu_vec(3), ex(3), ey(3), ez(3)
   real(dblprec) :: mu(2)
   mu(1) = dot_product(mu_vec, ex)
   mu(2) = dot_product(mu_vec, ey)
end function project_transverse

!======================================================================
!> @brief
!> Write the triangulation mesh to a file for visualization/debugging
!> Outputs both vertex coordinates and triangle connectivity
!> Applies minimal image correction to coordinates relative to origin for compact visualization
!======================================================================
subroutine print_triangulation_mesh(filename, coords, Natom, simid, C1, C2, C3, N1, N2, N3)
   use Constants
   implicit none

   character(len=*), intent(in) :: filename
   real(dblprec), intent(in) :: coords(3,*)  ! atom coordinates
   integer, intent(in) :: Natom              ! number of atoms
   character(len=*), intent(in) :: simid     ! simulation ID for filename
   real(dblprec), dimension(3), intent(in) :: C1, C2, C3
   integer, intent(in) :: N1, N2, N3

   integer :: i, j, ios, valid_count
   character(len=100) :: filn
   real(dblprec) :: x, y, z, x1, y1, z1, x2, y2, z2, x3, y3, z3, dist2, dist3, threshold, Lx, Ly, Lz

   ! Only print if flag is set
   if (print_mesh /= 'Y') return

   ! Create filename
   filn = trim(filename) // '.' // trim(simid) // '.mesh'

   open(unit=1001, file=trim(filn), status='replace', action='write', iostat=ios)
   if (ios /= 0) then
      write(*,*) 'Error opening mesh file: ', trim(filn)
      return
   end if

   ! Compute system dimensions
   Lx = real(N1, dblprec) * sqrt(dot_product(C1, C1))
   Ly = real(N2, dblprec) * sqrt(dot_product(C2, C2))
   Lz = real(N3, dblprec) * sqrt(dot_product(C3, C3))
   print *, 'System dimensions (Lx, Ly, Lz): ', Lx, Ly, Lz
   threshold = min(Lx, Ly) / 2.0_dblprec

   ! Write header
   write(1001,'(a)') '# Triangulation mesh file'
   write(1001,'(a,i8)') '# Total triangles in triangulation: ', nsimp
   write(1001,'(a)') '# Format: triangle_index vertex1_index vertex2_index vertex3_index'
   write(1001,'(a)') '# Followed by vertex coordinates: vertex_index x y z (minimal image corrected relative to origin)'
   write(1001,'(a,f12.6)') '# Triangle validity threshold: ', threshold

   ! Write triangle connectivity (only valid triangles)
   write(1001,'(a)') '# Triangle connectivity (only valid triangles):'
   valid_count = 0
   do i = 1, nsimp
      ! Get original coordinates
      x1 = coords(1,simp(1,i)); y1 = coords(2,simp(1,i)); z1 = coords(3,simp(1,i))
      x2 = coords(1,simp(2,i)); y2 = coords(2,simp(2,i)); z2 = coords(3,simp(2,i))
      x3 = coords(1,simp(3,i)); y3 = coords(2,simp(3,i)); z3 = coords(3,simp(3,i))
      
      ! Apply minimal image correction relative to vertex 1
      ! call minimal_image_correction(x2, y2, z2, x1, y1, z1, C1, C2, C3, N1, N2, N3)
      ! call minimal_image_correction(x3, y3, z3, x1, y1, z1, C1, C2, C3, N1, N2, N3)
      
      ! Check validity
      dist2 = sqrt((x2-x1)**2 + (y2-y1)**2 + (z2-z1)**2)
      dist3 = sqrt((x3-x1)**2 + (y3-y1)**2 + (z3-z1)**2)
      
      if (dist2 <= threshold .and. dist3 <= threshold) then
         valid_count = valid_count + 1
         write(1001,'(i8,3i12)') valid_count, simp(1,i), simp(2,i), simp(3,i)
      end if
   end do
   write(1001,'(a,i8)') '# Valid triangles written: ', valid_count

   ! Write vertex coordinates (with minimal image correction relative to origin)
   write(1001,'(a)') '# Vertex coordinates (minimal image corrected relative to origin):'
   do i = 1, Natom
      x = coords(1,i)
      y = coords(2,i)
      z = coords(3,i)
      ! Apply minimal image correction relative to origin
      ! call minimal_image_correction(x, y, z, 0.0_dblprec, 0.0_dblprec, 0.0_dblprec, C1, C2, C3, N1, N2, N3)
      write(1001,'(i8,3f16.8)') i, x, y, z
   end do

   close(1001)
   write(*,'(1x,a,a)') 'Triangulation mesh written to: ', trim(filn)

end subroutine print_triangulation_mesh

!===============================================================
!> @brief
!> Validation routine to check triangulation quality
!===============================================================
subroutine validate_triangulation(coords, Natom, C1, C2, C3, N1, N2, N3)
   use Constants
   implicit none
   real(dblprec), intent(in) :: coords(3,Natom)
   integer, intent(in) :: Natom
   real(dblprec), dimension(3), intent(in) :: C1, C2, C3
   integer, intent(in) :: N1, N2, N3

   integer :: isimp, i1,i2,i3, nerr
   real(dblprec) :: x1,y1,z1,x2,y2,z2,x3,y3,z3
   real(dblprec) :: d12,d23,d31, dmax, Lx,Ly,Lz, threshold
   real(dblprec) :: area, cross

   Lx = N1*sqrt(dot_product(C1,C1))
   Ly = N2*sqrt(dot_product(C2,C2))
   Lz = N3*sqrt(dot_product(C3,C3))
   threshold = 0.5_dblprec*min(Lx,Ly)

   nerr = 0
   dmax = 0.0_dblprec

   do isimp=1,nsimp
      i1 = simp(1,isimp); i2 = simp(2,isimp); i3 = simp(3,isimp)
      
      ! Check indices
      if (i1<1 .or. i1>Natom .or. i2<1 .or. i2>Natom .or. i3<1 .or. i3>Natom) then
         write(*,*) 'ERROR: Triangle ',isimp,' has invalid indices:',i1,i2,i3
         nerr = nerr + 1
         cycle
      endif

      x1=coords(1,i1); y1=coords(2,i1); z1=coords(3,i1)
      x2=coords(1,i2); y2=coords(2,i2); z2=coords(3,i2)
      x3=coords(1,i3); y3=coords(2,i3); z3=coords(3,i3)

      ! Check edge lengths
      d12 = sqrt((x2-x1)**2 + (y2-y1)**2 + (z2-z1)**2)
      d23 = sqrt((x3-x2)**2 + (y3-y2)**2 + (z3-z2)**2)
      d31 = sqrt((x1-x3)**2 + (y1-y3)**2 + (z1-z3)**2)
      dmax = max(dmax, d12, d23, d31)

      if (d12>threshold .or. d23>threshold .or. d31>threshold) then
         write(*,'(a,i6,a,3f8.3)') 'WARNING: Triangle ',isimp, &
                ' has long edge (PBC artifact?): ',d12,d23,d31
      endif

      ! Check area (degenerate triangles)
      cross = (x2-x1)*(y3-y1) - (y2-y1)*(x3-x1)
      area = 0.5_dblprec*abs(cross)
      if (area < 1.0e-10_dblprec) then
         write(*,*) 'WARNING: Triangle ',isimp,' has zero area'
      endif
   end do

   write(*,'(1x,a,i8)') 'Triangulation validation complete. Errors: ', nerr
   write(*,'(1x,a,f12.6)') 'Maximum edge length: ', dmax
   write(*,'(1x,a,f12.6)') 'PBC threshold: ', threshold
end subroutine validate_triangulation

end module Topology
