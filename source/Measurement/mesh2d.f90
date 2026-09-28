!------------------------------------------------------------------------------------
!> @brief Shared two-dimensional finite-element mesh for topology measurements.
!>
!> The mesh stores triangle connectivity together with local, minimum-image
!> geometry.  Consumers can therefore use the same orientation, areas and FEM
!> coefficients without modifying the system coordinates.
!------------------------------------------------------------------------------------
module Mesh2D

   use Parameters
   use Profiling

   implicit none

   private

   integer, public :: nsimp = 0
   integer, allocatable, public :: simp(:,:)          ! (3,nsimp), counter-clockwise
   real(dblprec), allocatable, public :: tri_area(:)  ! (nsimp) > 0
   real(dblprec), allocatable, public :: site_area(:) ! A_i = sum_D A_D / 3
   real(dblprec), allocatable, public :: site_wsum(:) ! sum_D A_D
   real(dblprec), allocatable, public :: grad_b(:,:)  ! (3,nsimp), FEM d/dx
   real(dblprec), allocatable, public :: grad_c(:,:)  ! (3,nsimp), FEM d/dy
   integer, allocatable, public :: site_tri_ptr(:)    ! CSR row pointers
   integer, allocatable, public :: site_tri_idx(:)    ! CSR triangle indices
   integer, public :: ndegenerate = 0
   real(dblprec), public :: mesh_cell_area = 0.0_dblprec
   logical, save :: non_xy_mesh_warning = .false.

   public :: mesh2d_build, mesh2d_report, mesh2d_release

contains

   !---------------------------------------------------------------------------------
   !> @brief Build the shared two-dimensional mesh.
   !>
   !> Periodic vertices are placed by minimum image relative to the first vertex
   !> of each cell.  Open directions omit their wrap cells.  Every z layer gets
   !> its own independent xy triangulation; triangles never connect layers.
   !---------------------------------------------------------------------------------
   subroutine mesh2d_build(N1,N2,N3,NA,coord,C1,C2,C3,BC1,BC2,BC3)

      implicit none

      integer, intent(in) :: N1, N2, N3, NA
      real(dblprec), intent(in) :: coord(3,*)
      real(dblprec), intent(in) :: C1(3), C2(3), C3(3)
      character(len=1), intent(in) :: BC1, BC2, BC3

      integer :: i_stat, i_all, max_tri, tri_count, total_inc
      integer :: x, y, z, ixp, iyp, it, i00, i10, i01, i11, layer_offset
      integer :: ia, iv, isite, tri_index, natom
      integer, allocatable :: simp_work(:,:), site_count(:), cursor(:)
      real(dblprec), allocatable :: area_work(:), b_work(:,:), c_work(:,:)
      real(dblprec) :: cell1(2), cell2(2), inv_cell(2,2), det_cell
      real(dblprec) :: r00(2), rel10(2), rel01(2), rel11(2)
      real(dblprec) :: p(3,2), d1(2), d2(2)
      real(dblprec) :: d2_diag1, d2_diag2
      integer :: vertices(3)
      logical :: periodic_x, periodic_y

      ! C3 is intentionally part of the public interface although this is an xy mesh.
      if (size(C3) < 3) error stop 'Mesh2D: invalid third lattice vector'
      if (N1 < 1 .or. N2 < 1 .or. N3 < 1 .or. NA < 1) then
         error stop 'Mesh2D: invalid mesh dimensions'
      end if

      call mesh2d_release()

      natom = N1*N2*N3*NA
      max_tri = max(1,2*N1*N2*N3*NA)
      periodic_x = BC1=='P'
      periodic_y = BC2=='P'

      if (N2==1 .and. N3>1 .and. .not.non_xy_mesh_warning) then
         write(*,'(1x,a)') 'WARNING: Mesh2D received N2=1 and N3>1; the system is not in an xy plane.'
         non_xy_mesh_warning=.true.
      end if

      cell1 = real(N1,dblprec)*C1(1:2)
      cell2 = real(N2,dblprec)*C2(1:2)
      det_cell = cell1(1)*cell2(2)-cell1(2)*cell2(1)
      mesh_cell_area=abs(det_cell)/real(N1*N2,dblprec)
      inv_cell = 0.0_dblprec
      if (abs(det_cell)>1.0e-14_dblprec) then
         inv_cell(1,1)= cell2(2)/det_cell
         inv_cell(1,2)=-cell2(1)/det_cell
         inv_cell(2,1)=-cell1(2)/det_cell
         inv_cell(2,2)= cell1(1)/det_cell
      end if

      allocate(simp_work(3,max_tri),stat=i_stat)
      call memocc(i_stat,product(shape(simp_work))*kind(simp_work),'simp_work','mesh2d_build')
      allocate(area_work(max_tri),stat=i_stat)
      call memocc(i_stat,product(shape(area_work))*kind(area_work),'area_work','mesh2d_build')
      allocate(b_work(3,max_tri),stat=i_stat)
      call memocc(i_stat,product(shape(b_work))*kind(b_work),'b_work','mesh2d_build')
      allocate(c_work(3,max_tri),stat=i_stat)
      call memocc(i_stat,product(shape(c_work))*kind(c_work),'c_work','mesh2d_build')

      tri_count = 0
      ndegenerate = 0

      do z=1,N3
         layer_offset=NA*N1*N2*(z-1)
         do y=1,N2
            if (y==N2 .and. .not.periodic_y) cycle
            iyp = modulo(y,N2)+1
            do x=1,N1
               if (x==N1 .and. .not.periodic_x) cycle
               ixp = modulo(x,N1)+1
               do it=1,NA
                  i00 = layer_offset+NA*((y-1)*N1+(x-1))+it
                  i10 = layer_offset+NA*((y-1)*N1+(ixp-1))+it
                  i01 = layer_offset+NA*((iyp-1)*N1+(x-1))+it
                  i11 = layer_offset+NA*((iyp-1)*N1+(ixp-1))+it

                  r00=coord(1:2,i00)
                  call minimum_image_2d(coord(1:2,i10)-r00,cell1,cell2,inv_cell, &
                     det_cell,periodic_x,periodic_y,rel10)
                  call minimum_image_2d(coord(1:2,i01)-r00,cell1,cell2,inv_cell, &
                     det_cell,periodic_x,periodic_y,rel01)
                  call minimum_image_2d(coord(1:2,i11)-r00,cell1,cell2,inv_cell, &
                     det_cell,periodic_x,periodic_y,rel11)

                  d2_diag1=dot_product(rel11,rel11)
                  d1=rel10-rel01
                  d2_diag2=dot_product(d1,d1)

                  if (d2_diag1<=d2_diag2) then
                     vertices=(/i00,i10,i11/)
                     p=0.0_dblprec
                     p(2,:)=rel10
                     p(3,:)=rel11
                     call append_triangle(vertices,p,tri_count,ndegenerate,simp_work, &
                        area_work,b_work,c_work)

                     vertices=(/i00,i11,i01/)
                     p=0.0_dblprec
                     p(2,:)=rel11
                     p(3,:)=rel01
                     call append_triangle(vertices,p,tri_count,ndegenerate,simp_work, &
                        area_work,b_work,c_work)
                  else
                     vertices=(/i00,i10,i01/)
                     p=0.0_dblprec
                     p(2,:)=rel10
                     p(3,:)=rel01
                     call append_triangle(vertices,p,tri_count,ndegenerate,simp_work, &
                        area_work,b_work,c_work)

                     vertices=(/i10,i11,i01/)
                     p(1,:)=rel10
                     p(2,:)=rel11
                     p(3,:)=rel01
                     call append_triangle(vertices,p,tri_count,ndegenerate,simp_work, &
                        area_work,b_work,c_work)
                  end if
               end do
            end do
         end do
      end do

      nsimp=tri_count
      allocate(simp(3,nsimp),stat=i_stat)
      call memocc(i_stat,product(shape(simp))*kind(simp),'simp','mesh2d_build')
      allocate(tri_area(nsimp),stat=i_stat)
      call memocc(i_stat,product(shape(tri_area))*kind(tri_area),'tri_area','mesh2d_build')
      allocate(grad_b(3,nsimp),stat=i_stat)
      call memocc(i_stat,product(shape(grad_b))*kind(grad_b),'grad_b','mesh2d_build')
      allocate(grad_c(3,nsimp),stat=i_stat)
      call memocc(i_stat,product(shape(grad_c))*kind(grad_c),'grad_c','mesh2d_build')
      if (nsimp>0) then
         simp=simp_work(:,1:nsimp)
         tri_area=area_work(1:nsimp)
         grad_b=b_work(:,1:nsimp)
         grad_c=c_work(:,1:nsimp)
      end if

      allocate(site_area(natom),stat=i_stat)
      call memocc(i_stat,product(shape(site_area))*kind(site_area),'site_area','mesh2d_build')
      allocate(site_wsum(natom),stat=i_stat)
      call memocc(i_stat,product(shape(site_wsum))*kind(site_wsum),'site_wsum','mesh2d_build')
      allocate(site_count(natom),stat=i_stat)
      call memocc(i_stat,product(shape(site_count))*kind(site_count),'site_count','mesh2d_build')
      site_area=0.0_dblprec
      site_wsum=0.0_dblprec
      site_count=0

      do tri_index=1,nsimp
         do iv=1,3
            isite=simp(iv,tri_index)
            site_area(isite)=site_area(isite)+tri_area(tri_index)/3.0_dblprec
            site_wsum(isite)=site_wsum(isite)+tri_area(tri_index)
            site_count(isite)=site_count(isite)+1
         end do
      end do

      allocate(site_tri_ptr(natom+1),stat=i_stat)
      call memocc(i_stat,product(shape(site_tri_ptr))*kind(site_tri_ptr),'site_tri_ptr','mesh2d_build')
      site_tri_ptr(1)=1
      do ia=1,natom
         site_tri_ptr(ia+1)=site_tri_ptr(ia)+site_count(ia)
      end do
      total_inc=site_tri_ptr(natom+1)-1
      allocate(site_tri_idx(total_inc),stat=i_stat)
      call memocc(i_stat,product(shape(site_tri_idx))*kind(site_tri_idx),'site_tri_idx','mesh2d_build')
      allocate(cursor(natom),stat=i_stat)
      call memocc(i_stat,product(shape(cursor))*kind(cursor),'cursor','mesh2d_build')
      cursor=site_tri_ptr(1:natom)
      do tri_index=1,nsimp
         do iv=1,3
            isite=simp(iv,tri_index)
            site_tri_idx(cursor(isite))=tri_index
            cursor(isite)=cursor(isite)+1
         end do
      end do

      i_all=-product(shape(cursor))*kind(cursor)
      deallocate(cursor,stat=i_stat)
      call memocc(i_stat,i_all,'cursor','mesh2d_build')
      i_all=-product(shape(site_count))*kind(site_count)
      deallocate(site_count,stat=i_stat)
      call memocc(i_stat,i_all,'site_count','mesh2d_build')
      i_all=-product(shape(simp_work))*kind(simp_work)
      deallocate(simp_work,stat=i_stat)
      call memocc(i_stat,i_all,'simp_work','mesh2d_build')
      i_all=-product(shape(area_work))*kind(area_work)
      deallocate(area_work,stat=i_stat)
      call memocc(i_stat,i_all,'area_work','mesh2d_build')
      i_all=-product(shape(b_work))*kind(b_work)
      deallocate(b_work,stat=i_stat)
      call memocc(i_stat,i_all,'b_work','mesh2d_build')
      i_all=-product(shape(c_work))*kind(c_work)
      deallocate(c_work,stat=i_stat)
      call memocc(i_stat,i_all,'c_work','mesh2d_build')

   end subroutine mesh2d_build

   !---------------------------------------------------------------------------------
   !> @brief Print the mesh diagnostic required by the trajectory OAM contract.
   !---------------------------------------------------------------------------------
   subroutine mesh2d_report()
      implicit none

      real(dblprec) :: total_area, cell_area

      total_area=0.0_dblprec
      if (allocated(tri_area)) total_area=sum(tri_area)
      cell_area=mesh_cell_area
      write(*,'(a,i0,a,es24.16,a,es24.16,a,i0)') 'Mesh2D: ntri= ',nsimp, &
         ' total_area= ',total_area,' cell_area= ',cell_area,' degenerate= ',ndegenerate
   end subroutine mesh2d_report

   !---------------------------------------------------------------------------------
   !> @brief Release all mesh-owned storage.
   !---------------------------------------------------------------------------------
   subroutine mesh2d_release()
      implicit none
      integer :: i_stat, i_all

      if (allocated(simp)) then
         i_all=-product(shape(simp))*kind(simp)
         deallocate(simp,stat=i_stat)
         call memocc(i_stat,i_all,'simp','mesh2d_release')
      end if
      if (allocated(tri_area)) then
         i_all=-product(shape(tri_area))*kind(tri_area)
         deallocate(tri_area,stat=i_stat)
         call memocc(i_stat,i_all,'tri_area','mesh2d_release')
      end if
      if (allocated(site_area)) then
         i_all=-product(shape(site_area))*kind(site_area)
         deallocate(site_area,stat=i_stat)
         call memocc(i_stat,i_all,'site_area','mesh2d_release')
      end if
      if (allocated(site_wsum)) then
         i_all=-product(shape(site_wsum))*kind(site_wsum)
         deallocate(site_wsum,stat=i_stat)
         call memocc(i_stat,i_all,'site_wsum','mesh2d_release')
      end if
      if (allocated(grad_b)) then
         i_all=-product(shape(grad_b))*kind(grad_b)
         deallocate(grad_b,stat=i_stat)
         call memocc(i_stat,i_all,'grad_b','mesh2d_release')
      end if
      if (allocated(grad_c)) then
         i_all=-product(shape(grad_c))*kind(grad_c)
         deallocate(grad_c,stat=i_stat)
         call memocc(i_stat,i_all,'grad_c','mesh2d_release')
      end if
      if (allocated(site_tri_ptr)) then
         i_all=-product(shape(site_tri_ptr))*kind(site_tri_ptr)
         deallocate(site_tri_ptr,stat=i_stat)
         call memocc(i_stat,i_all,'site_tri_ptr','mesh2d_release')
      end if
      if (allocated(site_tri_idx)) then
         i_all=-product(shape(site_tri_idx))*kind(site_tri_idx)
         deallocate(site_tri_idx,stat=i_stat)
         call memocc(i_stat,i_all,'site_tri_idx','mesh2d_release')
      end if
      nsimp=0
      ndegenerate=0
      mesh_cell_area=0.0_dblprec
   end subroutine mesh2d_release

   subroutine minimum_image_2d(delta,cell1,cell2,inv_cell,det_cell,periodic_x,periodic_y,local)
      implicit none
      real(dblprec), intent(in) :: delta(2), cell1(2), cell2(2), inv_cell(2,2), det_cell
      logical, intent(in) :: periodic_x, periodic_y
      real(dblprec), intent(out) :: local(2)
      real(dblprec) :: reduced(2), matrix(2,2)

      if (abs(det_cell)<=1.0e-14_dblprec) then
         local=delta
         return
      end if
      matrix(:,1)=cell1
      matrix(:,2)=cell2
      reduced=matmul(inv_cell,delta)
      if (periodic_x) reduced(1)=reduced(1)-real(nint(reduced(1)),dblprec)
      if (periodic_y) reduced(2)=reduced(2)-real(nint(reduced(2)),dblprec)
      local=matmul(matrix,reduced)
   end subroutine minimum_image_2d

   subroutine append_triangle(vertices,position,tri_count,degenerate_count,simp_work,area_work,b_work,c_work)
      implicit none
      integer, intent(in) :: vertices(3)
      real(dblprec), intent(inout) :: position(3,2)
      integer, intent(inout) :: tri_count, degenerate_count
      integer, intent(inout) :: simp_work(:,:)
      real(dblprec), intent(inout) :: area_work(:), b_work(:,:), c_work(:,:)

      integer :: tmp_index
      real(dblprec) :: cross, area, tmp_pos(2)
      integer :: ordered_vertices(3)

      cross=(position(2,1)-position(1,1))*(position(3,2)-position(1,2))- &
         (position(2,2)-position(1,2))*(position(3,1)-position(1,1))
      area=0.5_dblprec*abs(cross)
      if (area<1.0e-12_dblprec) then
         degenerate_count=degenerate_count+1
         return
      end if

      ordered_vertices=vertices
      if (cross<0.0_dblprec) then
         tmp_index=ordered_vertices(2)
         ordered_vertices(2)=ordered_vertices(3)
         ordered_vertices(3)=tmp_index
         tmp_pos=position(2,:)
         position(2,:)=position(3,:)
         position(3,:)=tmp_pos
         cross=-cross
      end if

      tri_count=tri_count+1
      simp_work(:,tri_count)=ordered_vertices
      area_work(tri_count)=0.5_dblprec*cross
      b_work(:,tri_count)=(/position(2,2)-position(3,2),position(3,2)-position(1,2), &
         position(1,2)-position(2,2)/)/cross
      c_work(:,tri_count)=(/position(3,1)-position(2,1),position(1,1)-position(3,1), &
         position(2,1)-position(1,1)/)/cross
   end subroutine append_triangle

end module Mesh2D
