!-------------------------------------------------------------------------------
! MODULE: HamiltonianTargetOrder
!> @brief Persistent, non-invasive CPU target traversal order.
!-------------------------------------------------------------------------------
module HamiltonianTargetOrder

   use Parameters

   implicit none
   private

   public :: setup_target_order
   public :: target_order_range

contains

   ! Build either the identity order or a deterministic 3-D Morton order.
   ! The values in target_order are the original physical atom IDs.
   subroutine setup_target_order(coord,nlistsize,aHam,target_order,target_work_prefix, &
         target_total_work,use_sfc,target_order_weighted)
      implicit none

      real(dblprec), intent(in) :: coord(:,:)
      integer, intent(in) :: nlistsize(:), aHam(:)
      integer, allocatable, intent(inout) :: target_order(:)
      integer(kind=8), allocatable, intent(inout) :: target_work_prefix(:)
      integer(kind=8), intent(out) :: target_total_work
      logical, intent(in) :: use_sfc
      logical, intent(out) :: target_order_weighted

      integer :: natom, i, ih, work_i, min_work, max_work
      integer(kind=8), allocatable :: keys(:), tmp_keys(:)
      integer, allocatable :: tmp_order(:)

      natom=size(aHam)
      if (size(coord,1) /= 3 .or. size(coord,2) /= natom) then
         error stop 'setup_target_order: coordinate dimensions do not match atom count'
      endif
      if (natom < 1 .or. size(nlistsize) < 1) then
         error stop 'setup_target_order: invalid Hamiltonian dimensions'
      endif
      if (any(aHam < 1) .or. any(aHam > size(nlistsize))) then
         error stop 'setup_target_order: invalid Hamiltonian lookup index'
      endif

      if (allocated(target_order)) deallocate(target_order)
      if (allocated(target_work_prefix)) deallocate(target_work_prefix)
      allocate(target_order(natom),target_work_prefix(natom))
      do i=1,natom
         target_order(i)=i
      end do

      if (use_sfc .and. natom > 1) then
         allocate(keys(natom),tmp_keys(natom),tmp_order(natom))
         call morton_keys(coord,keys)
         call sort_by_key(target_order,keys,tmp_order,tmp_keys)
         deallocate(keys,tmp_keys,tmp_order)
      endif

      target_work_prefix=0_8
      min_work=huge(min_work)
      max_work=0
      target_total_work=0_8
      do i=1,natom
         ih=aHam(target_order(i))
         work_i=max(1,nlistsize(ih))
         min_work=min(min_work,work_i)
         max_work=max(max_work,work_i)
         target_total_work=target_total_work+int(work_i,kind=8)
         target_work_prefix(i)=target_total_work
      end do
      target_order_weighted=(max_work /= min_work)

      if (use_sfc) then
         write(*,'(2x,a)') 'CPU Hamiltonian target order: MORTON'
      else
         write(*,'(2x,a)') 'CPU Hamiltonian target order: NATURAL'
      endif
      if (target_order_weighted) then
         write(*,'(2x,a)') 'CPU Hamiltonian partition: weighted static'
      else
         write(*,'(2x,a)') 'CPU Hamiltonian partition: static'
      endif
   end subroutine setup_target_order


   ! Return the inclusive sequence-index interval assigned to one thread.
   ! The cuts are based on the cumulative target work, so the atom order is
   ! contiguous and no dynamic scheduling is required.
   subroutine target_order_range(prefix,total_work,nthreads,thread_id,q_start,q_stop)
      implicit none

      integer(kind=8), intent(in) :: prefix(:), total_work
      integer, intent(in) :: nthreads, thread_id
      integer, intent(out) :: q_start, q_stop
      integer :: left_cut, right_cut

      if (nthreads < 1 .or. thread_id < 0 .or. thread_id >= nthreads) then
         error stop 'target_order_range: invalid thread partition'
      endif
      left_cut=target_cut(prefix,total_work,nthreads,thread_id)
      right_cut=target_cut(prefix,total_work,nthreads,thread_id+1)
      q_start=left_cut+1
      q_stop=right_cut
   end subroutine target_order_range


   integer function target_cut(prefix,total_work,nthreads,cut_number)
      implicit none

      integer(kind=8), intent(in) :: prefix(:), total_work
      integer, intent(in) :: nthreads, cut_number
      integer :: lo, hi, mid
      integer(kind=8) :: target

      if (cut_number <= 0) then
         target_cut=0
         return
      endif
      if (cut_number >= nthreads) then
         target_cut=size(prefix)
         return
      endif

      target=(total_work*int(cut_number,kind=8))/int(nthreads,kind=8)
      lo=1
      hi=size(prefix)
      do while (lo < hi)
         mid=(lo+hi)/2
         if (prefix(mid) >= target) then
            hi=mid
         else
            lo=mid+1
         endif
      end do
      target_cut=lo
   end function target_cut


   subroutine morton_keys(coord,keys)
      implicit none

      real(dblprec), intent(in) :: coord(:,:)
      integer(kind=8), intent(out) :: keys(:)
      real(dblprec) :: cmin(3), crange(3), scaled
      integer(kind=8) :: q(3)
      integer :: i, d, bit
      integer(kind=8), parameter :: max_quant=1048575_8

      do d=1,3
         cmin(d)=minval(coord(d,:))
         crange(d)=maxval(coord(d,:))-cmin(d)
      end do
      do i=1,size(coord,2)
         q=0_8
         do d=1,3
            if (crange(d) > 0.0_dblprec) then
               scaled=(coord(d,i)-cmin(d))/crange(d)*real(max_quant,dblprec)
               q(d)=min(max_quant,max(0_8,int(scaled,kind=8)))
            endif
         end do
         keys(i)=0_8
         do bit=0,19
            keys(i)=ior(keys(i),ishft(iand(q(1),ishft(1_8,bit)),2*bit))
            keys(i)=ior(keys(i),ishft(iand(q(2),ishft(1_8,bit)),2*bit+1))
            keys(i)=ior(keys(i),ishft(iand(q(3),ishft(1_8,bit)),2*bit+2))
         end do
      end do
   end subroutine morton_keys


   subroutine sort_by_key(order,keys,tmp_order,tmp_keys)
      implicit none

      integer, intent(inout) :: order(:)
      integer(kind=8), intent(inout) :: keys(:)
      integer, intent(out) :: tmp_order(:)
      integer(kind=8), intent(out) :: tmp_keys(:)
      integer :: width, left, mid, right, i, j, k

      width=1
      do while (width < size(order))
         left=1
         do while (left <= size(order))
            mid=min(left+width,size(order)+1)
            right=min(left+2*width,size(order)+1)
            i=left
            j=mid
            do k=left,right-1
               if (i >= mid) then
                  tmp_keys(k)=keys(j)
                  tmp_order(k)=order(j)
                  j=j+1
               else if (j >= right) then
                  tmp_keys(k)=keys(i)
                  tmp_order(k)=order(i)
                  i=i+1
               else if (less_pair(keys(i),order(i),keys(j),order(j))) then
                  tmp_keys(k)=keys(i)
                  tmp_order(k)=order(i)
                  i=i+1
               else
                  tmp_keys(k)=keys(j)
                  tmp_order(k)=order(j)
                  j=j+1
               endif
            end do
            left=right
         end do
         keys=tmp_keys
         order=tmp_order
         width=2*width
      end do
   end subroutine sort_by_key


   logical function less_pair(key_a,id_a,key_b,id_b)
      implicit none
      integer(kind=8), intent(in) :: key_a,key_b
      integer, intent(in) :: id_a,id_b

      less_pair=(key_a < key_b) .or. (key_a == key_b .and. id_a < id_b)
   end function less_pair

end module HamiltonianTargetOrder
