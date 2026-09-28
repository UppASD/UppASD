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
   use Mesh2D, only : nsimp, simp, site_tri_ptr, site_tri_idx

   ! Parameters for the printing
   integer :: skyno_step !< Interval for sampling the skyrmion number
   integer :: skyno_buff !< Buffer size for the sampling of the skyrmion number
   character(len=1) :: skyno      !< Perform skyrmion number measurement
   character(len=1) :: do_proj_skyno !< Perform type dependent skyrmion number measurement
   character(len=1) :: do_skyno_den  !< Perform site dependent skyrmion number measurement
   character(len=1) :: do_skyno_cmass  !< Perform center-of-mass skyrmion number measurement

   real(dblprec) :: chi_avg !< Average scalar chirality (instantaneous)
   real(dblprec) :: chi_cavg  = 0.0_dblprec !< Average scalar chirality (cumulative)
   real(dblprec), dimension(3) :: kappa_cavg = 0.0_dblprec !< Average scalar chirality (cumulative)
   real(dblprec), dimension(3) :: kappa_csum = 0.0_dblprec !< Cumulative sum of the vector chirality
   integer :: n_chi_cavg = 0 !< Number of times the average scalar chirality has been calculated

   character(len=1) :: print_mesh = 'N' !< Print triangulation mesh to file

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


      integer :: k, isimp
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


      integer :: k, isimp, isite
      real(dblprec), dimension(3) :: m1, m2, m3
      real(dblprec) :: q,qq, m1m2m3, m1m2, m1m3, m2m3

      real(dblprec), dimension(NA) :: pontryagin_tri_proj 
      real(dblprec), dimension(NA) :: thesum_proj

      thesum_proj=0.0_dblprec

      !!$omp parallel do default(shared) private(isimp,k,q,qq,m1,m2,m3,m1m2m3,m1m2,m1m3,m2m3,isite) reduction(+:thesum_proj)
      do k=1, Mensemble
         do isimp=1,nsimp
            isite=mod(simp(1,isimp)-1,NA)+1
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


      integer :: k, isimp, itri
      real(dblprec), dimension(3) :: m1, m2, m3
      real(dblprec) :: thesum,q,qq, m1m2m3, m1m2, m1m3, m2m3

      thesum=0.0_dblprec

      !isimp=2*iatom.

      if (.not.allocated(site_tri_ptr)) then
         pontryagin_tri_dens=0.0_dblprec
         return
      end if

      do k=1, Mensemble
         do itri=site_tri_ptr(iatom),site_tri_ptr(iatom+1)-1
            isimp=site_tri_idx(itri)
            m1 = emom(:,simp(1,isimp),k)
            m2 = emom(:,simp(2,isimp),k)
            m3 = emom(:,simp(3,isimp),k)
            m1m2m3=f_volume(m1,m2,m3)
            m1m2=dot_product(m1,m2)
            m1m3=dot_product(m1,m3)
            m2m3=dot_product(m2,m3)
            qq=m1m2m3/(1.0_dblprec+m1m2+m1m3+m2m3)
            q=2.0_dblprec*atan(qq)
            thesum=thesum+q/3.0_dblprec
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
   !> Calculates the scalar chirality using triangulation
   !
   !> @author
   !> Anders Bergman
   !---------------------------------------------------------------------------------
   function chirality_tri(Natom,Mensemble,emom) result(kappa_avg)
      use math_functions, only : f_cross_product
      implicit none
      integer, intent(in) :: Natom, Mensemble
      real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: emom

      real(dblprec), dimension(3) :: kappa_tot   ! global vector chirality
      real(dblprec), dimension(3) :: kappa_avg   ! average vector chirality
      real(dblprec) :: chi_tot                   ! optional scalar version
      real(dblprec), dimension(3) :: m1,m2,m3
      real(dblprec), dimension(3) :: c12,c23,c31
      integer :: k,isimp
      logical, save :: empty_mesh_warning = .false.

      if (nsimp==0) then
         if (.not.empty_mesh_warning) then
            write(*,'(1x,a)') 'WARNING: chirality_tri called with an empty mesh; returning zeros.'
            empty_mesh_warning=.true.
         end if
         kappa_avg=0.0_dblprec
         chi_avg=0.0_dblprec
         return
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
!======================================================================
!> @brief
!> Write the triangulation mesh to a file for visualization/debugging
!> Outputs both vertex coordinates and triangle connectivity
!> Writes the coordinates and connectivity currently held by the mesh
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

   integer :: i, ios, valid_count
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
   write(1001,'(a)') '# Followed by vertex coordinates: vertex_index x y z'
   write(1001,'(a,f12.6)') '# Triangle validity threshold: ', threshold

   ! Write triangle connectivity (only valid triangles)
   write(1001,'(a)') '# Triangle connectivity (only valid triangles):'
   valid_count = 0
   do i = 1, nsimp
      ! Get original coordinates
      x1 = coords(1,simp(1,i)); y1 = coords(2,simp(1,i)); z1 = coords(3,simp(1,i))
      x2 = coords(1,simp(2,i)); y2 = coords(2,simp(2,i)); z2 = coords(3,simp(2,i))
      x3 = coords(1,simp(3,i)); y3 = coords(2,simp(3,i)); z3 = coords(3,simp(3,i))
      
      ! Check validity
      dist2 = sqrt((x2-x1)**2 + (y2-y1)**2 + (z2-z1)**2)
      dist3 = sqrt((x3-x1)**2 + (y3-y1)**2 + (z3-z1)**2)
      
      if (dist2 <= threshold .and. dist3 <= threshold) then
         valid_count = valid_count + 1
         write(1001,'(i8,3i12)') valid_count, simp(1,i), simp(2,i), simp(3,i)
      end if
   end do
   write(1001,'(a,i8)') '# Valid triangles written: ', valid_count

   ! Write vertex coordinates
   write(1001,'(a)') '# Vertex coordinates:'
   do i = 1, Natom
      x = coords(1,i)
      y = coords(2,i)
      z = coords(3,i)
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
