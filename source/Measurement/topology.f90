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

   real(dblprec) :: chi_avg !< Average scalar chirality (instantaneous)
   real(dblprec) :: chi_cavg  = 0.0_dblprec !< Average scalar chirality (cumulative)
   real(dblprec), dimension(3) :: kappa_cavg = 0.0_dblprec !< Average scalar chirality (cumulative)
   real(dblprec), dimension(3) :: kappa_csum = 0.0_dblprec !< Cumulative sum of the vector chirality
   integer :: n_chi_cavg = 0 !< Number of times the average scalar chirality has been calculated

   integer, dimension(:,:), allocatable :: triangle_list !< List of triangles for each layer
   integer :: n_triangles !< Number of triangles in the triangulation
   integer :: n_layers !< Number of layers in the triangulation

   ! Variables for OAM 
   character(len=1) :: do_oam !< Perform OAM measurement
   integer :: oam_step !< Interval for sampling OAM
   integer :: oam_buff !< Buffer size for the sampling of OAM
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
         !do y=1,ny
         !   do x=1,nx
         do y=1,ny
            do x=1,nx
               do it=1,NT
                  ! Triangle 1  [0 -> +x -> +y -> 0]
                  nsimp=nsimp+1
                  simp(1,nsimp)=NT*wrap_idx(x,y,z,nx,ny,nz)+it-NT
                  ! simp(2,nsimp)=NT*wrap_idx(x+1,y,z,nx,ny,nz)+it-NT
                  ! simp(3,nsimp)=NT*wrap_idx(x,y+1,z,nx,ny,nz)+it-NT
                  simp(2,nsimp)=NT*wrap_idx(modulo(x,nx)+1,y,z,nx,ny,nz)+it-NT
                  simp(3,nsimp)=NT*wrap_idx(x,modulo(y,ny)+1,z,nx,ny,nz)+it-NT
               end do
            end do
         end do
      end do
      do z=1,nz
         !do y=1,ny
         !   do x=1,nx
         do y=1,ny
            do x=1,nx
               do it=1,NT
                  ! Triangle  2  [ 0 -> +y -> +y - x -> 0]
                  nsimp=nsimp+1
                  simp(1,nsimp)=NT*wrap_idx(x,y,z,nx,ny,nz)+it-NT
                  ! simp(2,nsimp)=NT*wrap_idx(x,y+1,z,nx,ny,nz)+it-NT
                  ! simp(3,nsimp)=NT*wrap_idx(x-1,y+1,z,nx,ny,nz)+it-NT
                  simp(2,nsimp)=NT*wrap_idx(x,modulo(y,ny)+1,z,nx,ny,nz)+it-NT
                  simp(3,nsimp)=NT*wrap_idx(modulo(x-2,nx)+1,modulo(y,ny)+1,z,nx,ny,nz)+it-NT
               end do
            end do
         end do
      end do


      ! if (nsimp .ne. 2*(nx)*(ny)*nz*nt) then
      !    write(*,'(1x,a,2i8)') 'Warning: Delaunay failure', nsimp, (2*nx-1)*(ny-1)*nz
      !    print *,shape(simp)
      ! end if


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

!======================================================================
   subroutine triangulate_layers ( nAtoms, coords )
      !----------------------------------------------------------------------
      ! Build triangle_list(3,n_triangles) for all atomic layers (z–slices)
      ! using a centroid-fan triangulation that rejects zero-area triangles.
      !----------------------------------------------------------------------
      
         implicit none
         integer,          intent(in)  :: nAtoms
         real(dblprec),    intent(in)  :: coords(3,nAtoms)
      
      
         ! ---- local -------------------------------------------------------
         integer                 :: i,j,k, count, out_idx, nLayers
         real(dblprec), allocatable :: zvals(:), uniqueZ(:)
         integer,       allocatable :: layerIdx(:)
         real(dblprec), allocatable :: angles(:), rho2(:)
         real(dblprec) :: cx, cy, dx, dy, tempAng, tempR
         integer       :: tempIdx
         real(dblprec) :: A, tol
         integer, allocatable :: itmp(:,:)
         real(dblprec), allocatable :: dtmp(:)
      
         tol = 1.0d-8
      
      !----------------------------------------------------------------------
      ! 1.  Collect unique z-coordinates  (layer detection)
      !----------------------------------------------------------------------
         allocate(zvals(nAtoms))
         zvals = coords(3,:)
      
         allocate(uniqueZ(nAtoms))      ! temporary oversize
         nLayers = 0
         do i = 1, nAtoms
            if (nLayers == 0) then
               nLayers = 1
               uniqueZ(1) = zvals(i)
            else
               do j = 1, nLayers
                  if (abs(zvals(i)-uniqueZ(j)) < tol) exit
               end do
               if (j > nLayers) then
                  nLayers = nLayers + 1
                  uniqueZ(nLayers) = zvals(i)
               end if
            end if
         end do
          allocate(dtmp(nLayers))
          dtmp = uniqueZ(1:nLayers)
          call move_alloc(dtmp, uniqueZ)   ! reshape uniqueZ to minimal size and transfer data
      
      !----------------------------------------------------------------------
      ! 2. First pass:  count triangles  (fan gives  count-2  per layer)
      !----------------------------------------------------------------------
         n_triangles = 0
         do k = 1, nLayers
            count = count_atoms_in_layer(coords, uniqueZ(k), tol)
            if (count >= 3) n_triangles = n_triangles + (count - 2)
         end do
         allocate(triangle_list(3, n_triangles))
      
         out_idx = 0         ! running index into triangle_list
      !----------------------------------------------------------------------
      ! 3. Process each layer
      !----------------------------------------------------------------------
         do k = 1, nLayers
            call collect_layer_indices(coords, uniqueZ(k), tol, layerIdx, count)
            if (count < 3) cycle
      
            ! ---- centroid --------------------------------------------------
            cx = sum(coords(1,layerIdx(1:count))) / real(count, dblprec)
            cy = sum(coords(2,layerIdx(1:count))) / real(count, dblprec)
      
            ! ---- polar coordinates (angle & radius²) -----------------------
            allocate(angles(0:count), rho2(0:count))
            do i = 1, count
               dx = coords(1,layerIdx(i)) - cx
               dy = coords(2,layerIdx(i)) - cy
               angles(i) = atan2(dy, dx)
               rho2(i)   = dx*dx + dy*dy
            end do
      
            ! ---- insertion sort by (angle , rho²) --------------------------
            do i = 2, count
               tempAng = angles(i) ; tempR = rho2(i) ; tempIdx = layerIdx(i)
               j = i - 1
               do while ( j>=1 .and. &
                 ( angles(j) > tempAng  .or. &
                   ( abs(angles(j)-tempAng) < 1.0e-12_dblprec .and. rho2(j) > tempR ) ) )
                  angles(j+1)   = angles(j)
                  rho2(j+1)     = rho2(j)
                  layerIdx(j+1) = layerIdx(j)
                  j = j - 1
               end do
               angles(j+1)   = tempAng
               rho2(j+1)     = tempR
               layerIdx(j+1) = tempIdx
            end do
      
            ! ---- fan triangulation, skip zero-area triples -----------------
            do i = 2, count-1
               A = 0.5d0 * ( (coords(1,layerIdx(i  ))-coords(1,layerIdx(1)))* &
                             (coords(2,layerIdx(i+1))-coords(2,layerIdx(1))) - &
                             (coords(1,layerIdx(i+1))-coords(1,layerIdx(1)))* &
                             (coords(2,layerIdx(i  ))-coords(2,layerIdx(1))) )
               if (abs(A) < 1.0e-10_dblprec) cycle   ! collinear -> skip
               out_idx = out_idx + 1
               triangle_list(1,out_idx) = layerIdx(1)
               triangle_list(2,out_idx) = layerIdx(i)
               triangle_list(3,out_idx) = layerIdx(i+1)
            end do
      
            deallocate(angles, rho2, layerIdx)
         end do
      
         ! update globals
         n_triangles = out_idx
         allocate(itmp(3, n_triangles))
         itmp = triangle_list(:,1:n_triangles)
         call move_alloc(itmp, triangle_list)   ! reshape triangle_list to minimal size and transfer data
         n_layers    = nLayers

         ! Pretty print the triangle list
         ! if (allocated(triangle_list)) then
         !    print *, "Triangle list: Total number of triangles =", size(triangle_list, dim=2)
         !    do i = 1, size(triangle_list, dim=2)
         !       write(*,'(A,I5,A,3(I5))') "Triangle ", i, ": ", triangle_list(1,i), triangle_list(2,i), triangle_list(3,i)
         !    end do
         ! else
         !    print *, "Triangle list is not allocated."
         ! end if
         ! print *, "Number of layers: ", n_layers
         ! print *, "Number of triangles: ", n_triangles
         ! print *,'__________________________________________________________'

      
      contains
      !----------------------------------------------------------------------
      function count_atoms_in_layer(c,z0,tol) result(cnt)
         real(dblprec), intent(in) :: c(3,nAtoms), z0, tol
         integer :: cnt, ii
         cnt = 0
         do ii=1,nAtoms
            if (abs(c(3,ii)-z0) < tol) cnt = cnt+1
         end do
      end function count_atoms_in_layer
      !----------------------------------------------------------------------
      subroutine collect_layer_indices(c,z0,tol,idx,ct)
         real(dblprec), intent(in) :: c(3,nAtoms), z0, tol
         integer, allocatable, intent(out) :: idx(:)
         integer,               intent(out) :: ct
         integer :: ii
         integer, allocatable :: idx_tmp(:)
      
         ct=0
         allocate(idx(nAtoms))     ! temporary oversize
         do ii=1,nAtoms
            if (abs(c(3,ii)-z0) < tol) then
               ct = ct+1
               idx(ct) = ii
            end if
         end do
         allocate(idx_tmp(ct))
         idx_tmp = idx(1:ct)
         call move_alloc(idx_tmp, idx)   ! reshape idx to minimal size and transfer data
      end subroutine collect_layer_indices
      !----------------------------------------------------------------------
      end subroutine triangulate_layers
      !======================================================================


!!!    subroutine triangulate_layers(nAtoms, coords)
!!!       implicit none
!!!       integer, intent(in) :: nAtoms         ! Total number of atoms
!!!       real(dblprec), dimension(3, nAtoms), intent(in) :: coords  ! Atom coordinates (x,y,z)
!!! 
!!!       integer :: i, j, k, nLayers, totalTri, idx, count, out_idx
!!!       real(dblprec), allocatable, dimension(:) :: zvals, uniqueZ
!!!       integer, allocatable, dimension(:) :: layerIdx
!!!       real(dblprec), allocatable, dimension(:) :: angles
!!!       integer :: tempIdx
!!!       real(dblprec) :: tempAngle
!!!       real(dblprec), allocatable :: tmp(:)
!!! 
!!!       real(dblprec) :: cx, cy
!!!       real(dblprec) :: tol
!!!       tol = 1.0e-8_dblprec
!!! 
!!!       !---------------------------------------------------------------------
!!!       ! Determine unique z values (layers)
!!!       allocate(zvals(nAtoms))
!!!       do i = 1, nAtoms
!!!           zvals(i) = coords(3,i)
!!!       end do
!!! 
!!!       allocate(uniqueZ(nAtoms))
!!!       nLayers = 0
!!!       do i = 1, nAtoms
!!!           if (nLayers == 0) then
!!!                nLayers = 1
!!!                uniqueZ(1) = zvals(i)
!!!           else
!!!                do j = 1, nLayers
!!!                    if (abs(zvals(i)-uniqueZ(j)) < tol) exit
!!!                end do
!!!                if (j > nLayers) then
!!!                    nLayers = nLayers + 1
!!!                    uniqueZ(nLayers) = zvals(i)
!!!                end if
!!!           end if
!!!       end do
!!!       if (nLayers < size(uniqueZ)) then
!!!          allocate(tmp(nLayers))
!!!          tmp = uniqueZ(1:nLayers)
!!!          deallocate(uniqueZ)
!!!          allocate(uniqueZ(nLayers))
!!!          uniqueZ = tmp
!!!          deallocate(tmp)
!!!       end if
!!!       !---------------------------------------------------------------------
!!!       ! First pass: count total number of triangles (for convex, non‐degenerate layers)
!!!       totalTri = 0
!!!       do k = 1, nLayers
!!!           count = 0
!!!           do i = 1, nAtoms
!!!                if (abs(coords(3,i)-uniqueZ(k)) < tol) count = count + 1
!!!           end do
!!!           if (count >= 3) totalTri = totalTri + (count - 2)
!!!       end do
!!! 
!!!       if (totalTri <= 0) then
!!!           allocate(triangle_list(3, 0))
!!!           return
!!!       end if
!!! 
!!!       allocate(triangle_list(3, totalTri))
!!!       out_idx = 0
!!! 
!!!       !---------------------------------------------------------------------
!!!       ! Process each layer separately: triangulate using a fan method.
!!!       ! For each layer, sort the points in order around the centroid.
!!!       do k = 1, nLayers
!!!           ! Get indices for atoms in the current layer
!!!           count = 0
!!!           allocate(layerIdx(nAtoms))
!!!           do i = 1, nAtoms
!!!                if (abs(coords(3,i)-uniqueZ(k)) < tol) then
!!!                    count = count + 1
!!!                    layerIdx(count) = i
!!!                end if
!!!           end do
!!! 
!!!           if (count < 3) then
!!!                deallocate(layerIdx)
!!!                cycle
!!!           end if
!!! 
!!!           ! Compute centroid in the xy-plane
!!!           cx = 0.0_dblprec
!!!           cy = 0.0_dblprec
!!!           do i = 1, count
!!!                cx = cx + coords(1, layerIdx(i))
!!!                cy = cy + coords(2, layerIdx(i))
!!!           end do
!!!           cx = cx / real(count, dblprec)
!!!           cy = cy / real(count, dblprec)
!!! 
!!!           ! Allocate temporary array to hold angle of each point relative to centroid
!!!           allocate(angles(0:count))
!!!           angles = 0.0_dblprec
!!!           do i = 1, count
!!!                angles(i) = atan2(coords(2, layerIdx(i)) - cy, coords(1, layerIdx(i)) - cx)
!!!           end do
!!! 
!!!           ! Simple insertion sort of layerIdx and angles based on the angle values
!!!           do i = 2, count
!!!                tempAngle = angles(i)
!!!                tempIdx = layerIdx(i)
!!!                j = i - 1
!!!                do while (j >= 1 .and. angles(j) > tempAngle)
!!!                    angles(j+1) = angles(j)
!!!                    layerIdx(j+1) = layerIdx(j)
!!!                    j = j - 1
!!!                end do
!!!                angles(j+1) = tempAngle
!!!                layerIdx(j+1) = tempIdx
!!!           end do
!!! 
!!!           ! Fan triangulation: use first point as anchor
!!!           do i = 2, count-1
!!!                out_idx = out_idx + 1
!!!                triangle_list(1, out_idx) = layerIdx(1)
!!!                triangle_list(2, out_idx) = layerIdx(i)
!!!                triangle_list(3, out_idx) = layerIdx(i+1)
!!!           end do
!!! 
!!!           deallocate(angles, layerIdx)
!!!       end do
!!! 
!!!       ! Pretty print the triangle list
!!!       ! if (allocated(triangle_list)) then
!!!       !    print *, "Triangle list: Total number of triangles =", size(triangle_list, dim=2)
!!!       !    do i = 1, size(triangle_list, dim=2)
!!!       !       write(*,'(A,I5,A,3(I5))') "Triangle ", i, ": ", triangle_list(1,i), triangle_list(2,i), triangle_list(3,i)
!!!       !    end do
!!!       ! else
!!!       !    print *, "Triangle list is not allocated."
!!!       ! end if
!!! 
!!!       n_triangles = totalTri
!!!       n_layers = nLayers
!!!    end subroutine triangulate_layers
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
      implicit none
      integer, intent(in) :: Natom, Mensemble
      real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: emom

      real(dblprec), dimension(3) :: kappa_tot   ! global vector chirality
      real(dblprec), dimension(3) :: kappa_avg   ! average vector chirality
      real(dblprec) :: chi_tot                   ! optional scalar version
      real(dblprec), dimension(3) :: m1,m2,m3
      real(dblprec), dimension(3) :: c12,c23,c31
      integer :: k,i_triangle

      kappa_tot = 0.0_dblprec
      chi_tot   = 0.0_dblprec

      !$omp parallel do default(shared) private(i_triangle,k,m1,m2,m3,c12,c23,c31)  &
      !$omp& reduction(+:kappa_tot,chi_tot)
      do i_triangle = 1, n_triangles
         do k = 1, Mensemble
            m1 = emom(:,triangle_list(1,i_triangle),k)
            m2 = emom(:,triangle_list(2,i_triangle),k)
            m3 = emom(:,triangle_list(3,i_triangle),k)

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

      kappa_avg = kappa_tot / real(n_triangles*Mensemble, dblprec)
      ! Chi avg stored in module data
      chi_avg = chi_tot / real(n_triangles*Mensemble, dblprec)
   end function chirality_tri

!===============================================================
!   Main driver routine  (flag = 0/1/2)
!===============================================================
   subroutine calculate_oam(Natom,Mensemble,emom,mstep,flag)
      use math_functions, only : f_cross_product
      use SystemData, only : coord
      use InputData, only : simid, N1, N2, N3, NA
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
         call delaunay_tri_tri(N1, N2, N3, NA)
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
      Lz_tot = oam_tri(Natom, Mensemble, emom, coord)
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
real(dblprec) function oam_tri(Natom, Mensemble, emom, coords) !, S0_arr)
   use Constants
   implicit none
   integer, intent(in) :: Natom
   integer, intent(in) :: Mensemble
   real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: emom
   real(dblprec), dimension(3,Natom), intent(in) :: coords   ! atom positions (x,y,z)
   ! real(dblprec), dimension(3,Natom), intent(in) :: S0_arr   ! reference spins at t=0

   integer :: k, isimp, i1,i2,i3
   real(dblprec) :: Lsum, rho_t, gamma_t, area_t
   real(dblprec) :: x1,y1,x2,y2,x3,y3, cross
   real(dblprec) :: th1, th2, th3, dth12, dth23, dth31
   real(dblprec) :: mu(2), norm

   ! local frame basis
   real(dblprec) :: ex(3), ey(3), ez(3)

   Lsum = 0.0_dblprec

   !$omp parallel do default(shared) private(isimp,k,i1,i2,i3,th1,th2,th3,dth12,dth23,dth31,gamma_t,rho_t, &
   !$omp   x1,y1,x2,y2,x3,y3,cross,area_t,mu,norm,ex,ey,ez) reduction(+:Lsum)
   do isimp=1,nsimp
      i1 = simp(1,isimp); i2 = simp(2,isimp); i3 = simp(3,isimp)
      do k=1,Mensemble
         ! coordinates
         x1=coords(1,i1); y1=coords(2,i1)
         x2=coords(1,i2); y2=coords(2,i2)
         x3=coords(1,i3); y3=coords(2,i3)

         rho_t = 0.0_dblprec

         ! --- vertex 1 ---
         call local_frame(S0_arr(:,i1), ex,ey,ez)
         mu = project_transverse(emom(:,i1,k)-S0_arr(:,i1), ex,ey,ez)
         th1 = atan2(mu(2),mu(1))
         norm = mu(1)*mu(1)+mu(2)*mu(2)
         rho_t = rho_t + norm

         ! --- vertex 2 ---
         call local_frame(S0_arr(:,i2), ex,ey,ez)
         mu = project_transverse(emom(:,i2,k)-S0_arr(:,i2), ex,ey,ez)
         th2 = atan2(mu(2),mu(1))
         norm = mu(1)*mu(1)+mu(2)*mu(2)
         rho_t = rho_t + norm

         ! --- vertex 3 ---
         call local_frame(S0_arr(:,i3), ex,ey,ez)
         mu = project_transverse(emom(:,i3,k)-S0_arr(:,i3), ex,ey,ez)
         th3 = atan2(mu(2),mu(1))
         norm = mu(1)*mu(1)+mu(2)*mu(2)
         rho_t = rho_t + norm

         rho_t = rho_t/3.0_dblprec

         ! unwrap phase differences
         dth12 = modulo(th2-th1+pi,2*pi)-pi
         dth23 = modulo(th3-th2+pi,2*pi)-pi
         dth31 = modulo(th1-th3+pi,2*pi)-pi
         gamma_t = dth12 + dth23 + dth31

         ! triangle area (2D cross product)
         cross = (x2-x1)*(y3-y1) - (y2-y1)*(x3-x1)
         area_t = 0.5_dblprec*abs(cross)

         ! OAM contribution
         Lsum = Lsum + rho_t * gamma_t * area_t
      end do
   end do
   !$omp end parallel do

   oam_tri = Lsum / Mensemble
   return
end function oam_tri

!===============================================================
!> Construct a local orthonormal frame (ex,ey,ez)
!> from a reference spin n0 (ez = n0/|n0|)
!===============================================================
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

!!! real(dblprec) function oam_tri(Natom, Mensemble, emom, coords)
!!!    use Constants
!!!    implicit none
!!!    integer, intent(in) :: Natom
!!!    integer, intent(in) :: Mensemble
!!!    real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: emom
!!!    real(dblprec), dimension(3,Natom), intent(in) :: coords   ! atom positions (x,y,z)
!!! 
!!!    integer :: k, isimp, i1,i2,i3
!!!    real(dblprec) :: Lsum, rho_t, gamma_t, area_t
!!!    real(dblprec) :: x1,y1,x2,y2,x3,y3, cross
!!!    real(dblprec) :: th1, th2, th3, dth12, dth23, dth31
!!!    real(dblprec) :: mx,my, norm
!!!    real(dblprec) :: n_hat(3), transverse(3), mu_vec(3)
!!! 
!!!    Lsum = 0.0_dblprec
!!! 
!!!    ! !$omp parallel do default(shared) private(isimp,k,i1,i2,i3,th1,th2,th3,dth12,dth23,dth31,gamma_t,rho_t, &
!!!    ! !     x1,y1,x2,y2,x3,y3,cross,area_t,mx,my,norm,n_hat,transverse,mu_vec) reduction(+:Lsum)
!!!    do isimp=1,nsimp
!!!       i1 = simp(1,isimp); i2 = simp(2,isimp); i3 = simp(3,isimp)
!!!       do k=1,Mensemble
!!!          ! coordinates
!!!          x1=coords(1,i1); y1=coords(2,i1)
!!!          x2=coords(1,i2); y2=coords(2,i2)
!!!          x3=coords(1,i3); y3=coords(2,i3)
!!! 
!!!          ! transverse components and phases for vertex 1
!!!          n_hat = S0_arr(:,i1)/norm2(S0_arr(:,i1))
!!!          mu_vec = emom(:,i1,k) - S0_arr(:,i1)  ! subtract ground state
!!!          transverse = mu_vec - dot_product(mu_vec, n_hat)*n_hat  ! project out longitudinal component
!!!          mx = transverse(1); my = transverse(2)
!!!          th1 = atan2(my,mx)
!!!          norm = mx*mx+my*my
!!!          rho_t = norm
!!! 
!!!          ! transverse components and phases for vertex 2
!!!          n_hat = S0_arr(:,i2)/norm2(S0_arr(:,i2))
!!!          mu_vec = emom(:,i2,k) - S0_arr(:,i2)  ! subtract ground state
!!!          transverse = mu_vec - dot_product(mu_vec, n_hat)*n_hat  ! project out longitudinal component
!!!          mx = transverse(1); my = transverse(2)
!!!          th2 = atan2(my,mx)
!!!          rho_t = rho_t + mx*mx+my*my
!!! 
!!!          ! transverse components and phases for vertex 3
!!!          n_hat = S0_arr(:,i3)/norm2(S0_arr(:,i3))
!!!          mu_vec = emom(:,i3,k) - S0_arr(:,i3)  ! subtract ground state
!!!          transverse = mu_vec - dot_product(mu_vec, n_hat)*n_hat  ! project out longitudinal component
!!!          mx = transverse(1); my = transverse(2)
!!!          th3 = atan2(my,mx)
!!!          rho_t = rho_t + mx*mx+my*my
!!! 
!!!          rho_t = rho_t/3.0_dblprec
!!! 
!!!          ! unwrap phase differences
!!!          dth12 = modulo(th2-th1+pi,2*pi)-pi
!!!          dth23 = modulo(th3-th2+pi,2*pi)-pi
!!!          dth31 = modulo(th1-th3+pi,2*pi)-pi
!!! 
!!!          gamma_t = dth12 + dth23 + dth31   ! total winding on triangle
!!! 
!!!          ! triangle area (2D cross product)
!!!          cross = (x2-x1)*(y3-y1) - (y2-y1)*(x3-x1)
!!!          area_t = 0.5_dblprec*abs(cross)
!!! 
!!!          ! OAM contribution from this triangle
!!!          Lsum = Lsum + rho_t * gamma_t * area_t
!!!       end do
!!!    end do
!!!    !!$omp end parallel do
!!! 
!!!    oam_tri = Lsum / Mensemble
!!!    ! print *, 'OAM (triangulation) = ', oam_tri / real(Natom*Mensemble, dblprec)
!!!    return
!!! end function oam_tri


end module Topology
