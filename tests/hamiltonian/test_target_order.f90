program test_target_order

   use Parameters
   use HamiltonianData, only : ham
   use HamiltonianActions, only : effective_field
   use InputData, only : ham_inp
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
   call benchmark_natural_target_traversal()
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

   subroutine benchmark_natural_target_traversal()
      integer, parameter :: benchmark_natom=4096, benchmark_ensembles=4
      integer, parameter :: benchmark_neighbours=6, repetitions=40
      integer :: atom,slot,ensemble,repeat,count_rate,count_start,count_stop
      real(dblprec) :: spin(3,benchmark_natom,benchmark_ensembles)
      real(dblprec) :: mmom(benchmark_natom,benchmark_ensembles)
      real(dblprec) :: external_field(3,benchmark_natom,benchmark_ensembles)
      real(dblprec) :: time_field(3,benchmark_natom,benchmark_ensembles)
      real(dblprec) :: beff(3,benchmark_natom,benchmark_ensembles)
      real(dblprec) :: beff1(3,benchmark_natom,benchmark_ensembles)
      real(dblprec) :: beff2(3,benchmark_natom,benchmark_ensembles)
      real(dblprec) :: emomM_macro(3,1,benchmark_ensembles)
      real(dblprec) :: energy,native_seconds,ordered_seconds
      integer :: cell_index(benchmark_natom),macro_nlistsize(1)

      ham_inp%do_jtensor=0
      ham_inp%exc_inter='N'
      ham_inp%do_dm=0
      ham_inp%do_sa=0
      ham_inp%do_pd=0
      ham_inp%do_biqdm=0
      ham_inp%do_bq=0
      ham_inp%do_ring=0
      ham_inp%do_chir=0
      ham_inp%do_anisotropy=0
      ham_inp%do_dip=0
      ham_inp%mult_axis='N'

      allocate(ham%aHam(benchmark_natom),ham%nlistsize(benchmark_natom), &
         ham%nlist(benchmark_neighbours,benchmark_natom),ham%ncoup(benchmark_neighbours,1,1))
      ham%aHam=1
      ham%nlistsize=benchmark_neighbours
      ham%ncoup=0.0_dblprec
      do slot=1,benchmark_neighbours
         ham%ncoup(slot,1,1)=0.1_dblprec*real(slot,dblprec)
      end do
      do atom=1,benchmark_natom
         do slot=1,benchmark_neighbours
            ham%nlist(slot,atom)=modulo(atom+slot-1,benchmark_natom)+1
         end do
      end do
      do ensemble=1,benchmark_ensembles
         do atom=1,benchmark_natom
            mmom(atom,ensemble)=1.0_dblprec
            spin(:,atom,ensemble)=(/sin(real(atom+ensemble,dblprec)), &
               cos(real(2*atom+ensemble,dblprec)),0.1_dblprec*real(atom,dblprec)/)
         end do
      end do
      external_field=0.0_dblprec
      time_field=0.0_dblprec
      emomM_macro=0.0_dblprec
      cell_index=1
      macro_nlistsize=benchmark_natom
      allocate(ham%target_order(benchmark_natom))
      ham%target_order=(/(atom,atom=1,benchmark_natom)/)
      ham%target_order_weighted=.false.
      ham%target_order_sfc=.false.

      call system_clock(count_start,count_rate=count_rate)
      do repeat=1,repetitions
         call effective_field(benchmark_natom,benchmark_ensembles,1,benchmark_natom,spin,mmom, &
            external_field,time_field,beff,beff1,beff2,energy,1,cell_index,emomM_macro, &
            macro_nlistsize,1,16,16,16,measure_energy=.false.)
      end do
      call system_clock(count_stop)
      native_seconds=real(count_stop-count_start,dblprec)/real(count_rate,dblprec)

      ham%target_order_sfc=.true.
      call system_clock(count_start,count_rate=count_rate)
      do repeat=1,repetitions
         call effective_field(benchmark_natom,benchmark_ensembles,1,benchmark_natom,spin,mmom, &
            external_field,time_field,beff,beff1,beff2,energy,1,cell_index,emomM_macro, &
            macro_nlistsize,1,16,16,16,measure_energy=.false.)
      end do
      call system_clock(count_stop)
      ordered_seconds=real(count_stop-count_start,dblprec)/real(count_rate,dblprec)
      call check(native_seconds > 0.0_dblprec .and. ordered_seconds > 0.0_dblprec, &
         'natural target traversal benchmark records both paths')
      write(*,'(a,3(es12.4,1x),a,i0)') 'CPU-HAM-08 natural target traversal native/ordered/ratio=', &
         native_seconds,ordered_seconds,ordered_seconds/native_seconds,' ensembles=',benchmark_ensembles

      deallocate(ham%target_order,ham%aHam,ham%nlistsize,ham%nlist,ham%ncoup)
   end subroutine benchmark_natural_target_traversal

end program test_target_order
