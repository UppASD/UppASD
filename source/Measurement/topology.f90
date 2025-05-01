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

   subroutine triangulate_layers(nAtoms, coords)
      implicit none
      integer, intent(in) :: nAtoms         ! Total number of atoms
      real(dblprec), dimension(3, nAtoms), intent(in) :: coords  ! Atom coordinates (x,y,z)

      integer :: i, j, k, nLayers, totalTri, idx, count, out_idx
      real(dblprec), allocatable, dimension(:) :: zvals, uniqueZ
      integer, allocatable, dimension(:) :: layerIdx
      real(dblprec), allocatable, dimension(:) :: angles
      integer :: tempIdx
      real(dblprec) :: tempAngle
      real(dblprec), allocatable :: tmp(:)

      real(dblprec) :: cx, cy
      real(dblprec) :: tol
      tol = 1.0e-8_dblprec

      !---------------------------------------------------------------------
      ! Determine unique z values (layers)
      allocate(zvals(nAtoms))
      do i = 1, nAtoms
          zvals(i) = coords(3,i)
      end do

      allocate(uniqueZ(nAtoms))
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
      if (nLayers < size(uniqueZ)) then
         allocate(tmp(nLayers))
         tmp = uniqueZ(1:nLayers)
         deallocate(uniqueZ)
         allocate(uniqueZ(nLayers))
         uniqueZ = tmp
         deallocate(tmp)
      end if
      !---------------------------------------------------------------------
      ! First pass: count total number of triangles (for convex, non‐degenerate layers)
      totalTri = 0
      do k = 1, nLayers
          count = 0
          do i = 1, nAtoms
               if (abs(coords(3,i)-uniqueZ(k)) < tol) count = count + 1
          end do
          if (count >= 3) totalTri = totalTri + (count - 2)
      end do

      if (totalTri <= 0) then
          allocate(triangle_list(3, 0))
          return
      end if

      allocate(triangle_list(3, totalTri))
      out_idx = 0

      !---------------------------------------------------------------------
      ! Process each layer separately: triangulate using a fan method.
      ! For each layer, sort the points in order around the centroid.
      do k = 1, nLayers
          ! Get indices for atoms in the current layer
          count = 0
          allocate(layerIdx(nAtoms))
          do i = 1, nAtoms
               if (abs(coords(3,i)-uniqueZ(k)) < tol) then
                   count = count + 1
                   layerIdx(count) = i
               end if
          end do

          if (count < 3) then
               deallocate(layerIdx)
               cycle
          end if

          ! Compute centroid in the xy-plane
          cx = 0.0_dblprec
          cy = 0.0_dblprec
          do i = 1, count
               cx = cx + coords(1, layerIdx(i))
               cy = cy + coords(2, layerIdx(i))
          end do
          cx = cx / real(count, dblprec)
          cy = cy / real(count, dblprec)

          ! Allocate temporary array to hold angle of each point relative to centroid
          allocate(angles(count))
          do i = 1, count
               angles(i) = atan2(coords(2, layerIdx(i)) - cy, coords(1, layerIdx(i)) - cx)
          end do

          ! Simple insertion sort of layerIdx and angles based on the angle values
          do i = 2, count
               tempAngle = angles(i)
               tempIdx = layerIdx(i)
               j = i - 1
               do while (j >= 1 .and. angles(j) > tempAngle)
                   angles(j+1) = angles(j)
                   layerIdx(j+1) = layerIdx(j)
                   j = j - 1
               end do
               angles(j+1) = tempAngle
               layerIdx(j+1) = tempIdx
          end do

          ! Fan triangulation: use first point as anchor
          do i = 2, count-1
               out_idx = out_idx + 1
               triangle_list(1, out_idx) = layerIdx(1)
               triangle_list(2, out_idx) = layerIdx(i)
               triangle_list(3, out_idx) = layerIdx(i+1)
          end do

          deallocate(angles, layerIdx)
      end do

      ! Pretty print the triangle list
      ! if (allocated(triangle_list)) then
      !    print *, "Triangle list: Total number of triangles =", size(triangle_list, dim=2)
      !    do i = 1, size(triangle_list, dim=2)
      !       write(*,'(A,I5,A,3(I5))') "Triangle ", i, ": ", triangle_list(1,i), triangle_list(2,i), triangle_list(3,i)
      !    end do
      ! else
      !    print *, "Triangle list is not allocated."
      ! end if

      n_triangles = totalTri
      n_layers = nLayers
   end subroutine triangulate_layers
   !---------------------------------------------------------------------------------
   !> @brief
   !> Calculates the scalar chirality using triangulation
   !
   !> @author
   !> Anders Bergman
   !---------------------------------------------------------------------------------
   function chirality_tri(Natom,Mensemble,emom) result(kappa_avg)
      use constants
      use math_functions
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

      ! print *,'-------------------------------'
      ! print '(3f10.4)', emom(:,:,1)
      ! print *,'============================='
      !$omp parallel do default(shared) private(i_triangle,k,m1,m2,m3,c12,c23,c31)  &
      !$omp& reduction(+:kappa_tot,chi_tot)
      do i_triangle = 1, n_triangles
         ! print *, 'Triangle', i_triangle, ':'
         ! print *, "Corners = ", triangle_list(:,i_triangle)
         do k = 1, Mensemble
            m1 = emom(:,triangle_list(1,i_triangle),k)
            m2 = emom(:,triangle_list(2,i_triangle),k)
            m3 = emom(:,triangle_list(3,i_triangle),k)
            ! print *, 'm1 = ', m1
            ! print *, 'm2 = ', m2
            ! print *, 'm3 = ', m3

            ! pairwise cross-products
            c12 = f_cross_product(m1,m2)   ! c12 = m1 × m2
            c23 = f_cross_product(m2,m3)
            c31 = f_cross_product(m3,m1)
            !   print *, 'c12 = ', c12
            !   print *, 'c23 = ', c23
            !   print *, 'c31 = ', c31

            ! accumulate vector chirality for this triangle
            kappa_tot = kappa_tot + c12 + c23 + c31

            ! optional: scalar chirality
            chi_tot = chi_tot + dot_product(m1,c23)   ! m1·(m2×m3)
         end do
      end do
      !$omp end parallel do

      ! print *, 'kappa_tot = ', kappa_tot
      kappa_avg = kappa_tot / real(n_triangles*Mensemble, dblprec)
      ! print *, 'kappa_avg = ', kappa_avg
      ! print *, 'kappa_avg = ', kappa_avg
      ! Chi avg stored in module data
      chi_avg = chi_tot / real(n_triangles*Mensemble, dblprec)
      ! print *, 'chi_avg = ', chi_avg
   end function chirality_tri

   !!! function knb_pol_tri (Natom, Mensemble, emom)  result (P_avg)
   !!!    !----------------------------------------------------------------------
   !!!    !   Triangle-based Katsura–Nagaosa–Balatsky polarisation
   !!!    !   Each triangle is taken once with a fixed orientation.
   !!!    !   The three bond contributions are summed inside that triangle.
   !!!    !   Because every bond is shared by two triangles, we divide by 2
   !!!    !   at the end to eliminate the double counting.
   !!!    !----------------------------------------------------------------------
   !!!    use constants                                    ! includes dblprec, pi, C_KNB
   !!!    use math_functions                               ! cross_product()
   !!!    implicit none

   !!!    integer,  intent(in) :: Natom, Mensemble
   !!!    real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: emom

   !!!    ! -- return value
   !!!    real(dblprec), dimension(3) :: P_avg            ! averaged over all triangles & ensembles

   !!!    ! -- local variables
   !!!    real(dblprec), dimension(3) :: P_tot            ! running total
   !!!    real(dblprec), dimension(3) :: m1, m2, m3
   !!!    real(dblprec), dimension(3) :: sxs, sxt, sxu    ! cross-products for each bond
   !!!    integer :: isimp, k
   !!!    integer :: i1,i2,i3                             ! vertex indices for the triangle


   !!!    ! Initialisation
   !!!    P_tot = 0.0_dblprec

   !!! !$omp parallel do default(shared) private(isimp,k,m1,m2,m3,sxs,sxt,sxu,i1,i2,i3) reduction(+:P_tot)
   !!!    do isimp = 1, nsimp
   !!!       i1 = simp(1,isimp)
   !!!       i2 = simp(2,isimp)
   !!!       i3 = simp(3,isimp)

   !!!       do k = 1, Mensemble
   !!!          ! -------- get the three spins on this triangle
   !!!          m1 = emom(:, i1, k)
   !!!          m2 = emom(:, i2, k)
   !!!          m3 = emom(:, i3, k)

   !!!          ! -------- e_12 × (S1×S2)
   !!!          call cross_product (m1, m2, sxs)                 ! S1×S2
   !!!          call cross_product (evec(:,i1,i2), sxs, sxs)     ! e_12 × ...

   !!!          ! -------- e_23 × (S2×S3)
   !!!          call cross_product (m2, m3, sxt)
   !!!          call cross_product (evec(:,i2,i3), sxt, sxt)

   !!!          ! -------- e_31 × (S3×S1)
   !!!          call cross_product (m3, m1, sxu)
   !!!          call cross_product (evec(:,i3,i1), sxu, sxu)

   !!!          ! -------- add to global sum
   !!!          P_tot = P_tot + sxs + sxt + sxu
   !!!       end do
   !!!    end do
   !!! !$omp end parallel do

   !!!    ! -- normalise:  divide by 2 for double counting, by nsimp*Mens for ensemble,
   !!!    P_avg =P_tot / real(2*nsimp*Mensemble, dblprec)

   !!! end function knb_pol_tri

end module Topology
