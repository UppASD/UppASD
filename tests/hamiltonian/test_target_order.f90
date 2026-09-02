program test_target_order

   use Parameters
   use HamiltonianTargetOrder, only : setup_target_order, target_order_range

   implicit none

   integer, parameter :: natom=8, nthreads=3
   real(dblprec) :: coord(3,natom)
   integer :: aHam(natom), nlistsize(2), order(natom)
   integer, allocatable :: target_order(:)
   integer(kind=8), allocatable :: prefix(:)
   integer(kind=8) :: total
   logical :: weighted
   integer :: q_start, q_stop, thread_id, previous_stop

   coord=0.0_dblprec
   coord(:,1)=(/1.0_dblprec,1.0_dblprec,0.0_dblprec/)
   coord(:,2)=(/0.0_dblprec,0.0_dblprec,0.0_dblprec/)
   coord(:,3)=(/1.0_dblprec,0.0_dblprec,0.0_dblprec/)
   coord(:,4)=(/0.0_dblprec,1.0_dblprec,0.0_dblprec/)
   coord(:,5)=(/0.0_dblprec,0.0_dblprec,1.0_dblprec/)
   coord(:,6)=(/1.0_dblprec,0.0_dblprec,1.0_dblprec/)
   coord(:,7)=(/0.0_dblprec,1.0_dblprec,1.0_dblprec/)
   coord(:,8)=(/1.0_dblprec,1.0_dblprec,1.0_dblprec/)
   aHam=(/1,2,1,2,1,2,1,2/)
   nlistsize=(/2,2/)

   call setup_target_order(coord,nlistsize,aHam,target_order,prefix,total,.false.,weighted)
   order=(/(target_order(thread_id),thread_id=1,natom)/)
   call check(all(order==(/1,2,3,4,5,6,7,8/)), 'natural target order is identity')
   call check(.not.weighted, 'homogeneous order reports ordinary static partition')
   call check(total==16_8, 'natural cumulative work total is exact')

   nlistsize=(/1,3/)
   call setup_target_order(coord,nlistsize,aHam,target_order,prefix,total,.true.,weighted)
   order=(/(target_order(thread_id),thread_id=1,natom)/)
   call check(all(order==(/2,3,4,1,5,6,7,8/)), 'Morton order is deterministic')
   call check(weighted, 'heterogeneous order enables weighted partitioning')
   call check(is_permutation(order), 'Morton order preserves physical atom IDs')

   previous_stop=0
   do thread_id=0,nthreads-1
      call target_order_range(prefix,total,nthreads,thread_id,q_start,q_stop)
      call check(q_start==previous_stop+1, 'weighted ranges are contiguous')
      previous_stop=q_stop
   end do
   call check(previous_stop==natom, 'weighted ranges cover every target')

   deallocate(target_order,prefix)
   write(*,'(a)') 'CPU-HAM-02A target-order tests passed'

contains

   logical function is_permutation(values)
      integer, intent(in) :: values(:)
      integer :: i
      is_permutation=.true.
      do i=1,size(values)
         if (count(values==i) /= 1) is_permutation=.false.
      end do
   end function is_permutation

   subroutine check(condition,label)
      logical, intent(in) :: condition
      character(len=*), intent(in) :: label
      if (.not.condition) error stop 'FAILED: '//trim(label)
   end subroutine check

end program test_target_order
