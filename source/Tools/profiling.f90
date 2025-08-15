!> Routines for monitoing time and memory usage
!> @author
!> Luigi Genovese, Anders Bergman
!> @copyright
!> GNU Public License.
module Profiling
   use Parameters

   implicit none

   integer :: do_meminfo = 0     !< Print memory allocation info to `meminfo` file (0/1)

   public

contains

   !control the memory occupation by calculating the overall size in bytes of the allocated arrays
   !usage:
   ! when allocating allocating an array "stuff" of dimension n in the routine "dosome"
   ! allocate(stuff(n),stat=i_stat)
   ! call memocc(i_stat,product(shape(stuff))*kind(stuff),'stuff','dosome')
   ! when deallocating
   ! i_all=-product(shape(stuff))*kind(stuff)
   ! deallocate(stuff,stat=i_stat)
   ! call memocc(i_stat,i_all,'stuff','dosome')
   ! the counters are initialized with
   ! call memocc(0,iproc,'count','start') (iproc = mpi rank, nproc=mpi size)
   ! and stopped with
   ! call memocc(0,0,'count','stop')
   ! at the end of the calculation a short report is printed on the screen
   ! some information can be also written on disk following the needs
   ! This file is distributed under the terms of the
   ! GNU General Public License, see http://www.gnu.org/copyleft/gpl.txt .
   ! Copyright (C) Luigi Genovese, CEA Grenoble, France, 2007
   !> Memory profiling routine
   subroutine memocc(istat,isize,array,routine)
      implicit none
      character(len=*), intent(in) :: array
      character(len=*), intent(in) :: routine
      integer, intent(in) :: istat
      integer, intent(in) :: isize

      ! Local variables
      character(len=36) :: maxroutine,locroutine
      character(len=36) :: maxarray,locarray
      integer :: nalloc,ndealloc,locpeak,locmemory,iproc
      integer :: dblsize
      integer(kind=8) :: memory,maxmemory
      character(len=1) :: allocationflag

      save :: memory,nalloc,ndealloc,maxroutine,maxarray,maxmemory
      save :: locroutine,locarray,locpeak,locmemory,iproc


      dblsize=1
      !
      select case(array)
      case('count')
         if (routine=='start') then
            memory=0
            maxmemory=0
            nalloc=0
            ndealloc=0
            locroutine='routine'
            locarray='array'
            locmemory=0
            locpeak=0
            iproc=isize
            !open the writing file for the root process
            if (iproc == 0) then
               open(unit=mfileno,file='meminfo',status='unknown')
               if (do_meminfo>0) then
                  write(mfileno,'(a32,1x,a20,3(1x,a12))')&
                     '(Data in kB)             Routine','    Peak Array',&
                     'Routine Mem','Total Mem','Action'
               end if
            end if
         else if (routine=='stop' .and. iproc==0) then
            write(*,'(1x,a)')&
               '-------------------MEMORY CONSUMPTION REPORT-----------------------'
            write(*,'(1x,2(i0,a),i0)')&
               nalloc,' allocations and ',ndealloc,' deallocations, remaining memory(B):',memory
            write(*,'(1x,a,i0,a)') 'Memory occupation peak: ',maxmemory/1024,' kB'
            write(*,'(1x,5(a))') 'For the array: "',trim(maxarray),'" in routine "',trim(maxroutine),'"'
            write(*,'(1x,a)')&
               '-----------------END MEMORY CONSUMPTION REPORT---------------------'
            if (do_meminfo > 0 ) then
               close(mfileno)
            else
               close(mfileno,status='delete')
            end if
         end if

      case default
         !control of the allocation/deallocation status
         if (istat/=0) then
            if (isize>=0) then
               write(*,*)' subroutine ',routine,': problem of allocation of array ',array
               stop
            else if (isize<0) then
               write(*,*)' subroutine ',routine,': problem of deallocation of array ',array
               stop
            end if
         end if
         select case(iproc)
         case (0)
            ! To be used for inspecting an array which is not deallocated
            if (isize>0) then
               allocationflag = 'A'
            elseif (isize<0) then
               allocationflag = 'D'
            else
               allocationflag = '?'
            end if

            if(do_meminfo>0) write(mfileno,'(a32,1x,a20,3(1x,i12),1x,a12)') &
               trim(routine),trim(array),isize*dblsize,memory,maxmemory,allocationflag
            if (trim(locroutine) /= routine) then
               locroutine=routine
               locmemory=isize*dblsize
               locpeak=isize*dblsize
            else
               locmemory=locmemory+isize*dblsize
               if (locmemory > locpeak) then
                  locpeak=locmemory
               end if
            end if
            locarray=array
            memory=memory+isize*dblsize
            if (memory > maxmemory) then
               maxmemory=memory
               maxroutine=routine
               maxarray=array
            end if
            if (isize>0) then
               nalloc=nalloc+1
            else if (isize<0) then
               ndealloc=ndealloc+1
            end if
         case default
            return
         end select

      end select

   end subroutine memocc

   !> Time profining routine
   subroutine timing(iproc,category,action)
      implicit none
      !Variables
      integer, intent(in) :: iproc !< Current processor
      character(len=14), intent(in) :: category
      character(len=2), intent(in) :: action ! possibilities: INitialize, ON, OFf, REsults

      !Local variables
      logical :: parallel,init
      integer, parameter :: ncat=14 ! define timing categories
      integer :: i,ii,nprocs
      integer, external :: omp_get_num_threads,OMP_GET_NUM_PROCS

      !cputime routine gives a real
      real(kind=8) :: total,total0,time,time0
      real(kind=8) :: pc,total_pc
      real(kind=8) :: flops(ncat),timesum(ncat+1),timemax(ncat+1),timemin(ncat+1)
      save :: time0,init,timesum,total0,parallel

      character(len=14), dimension(ncat), parameter :: cats = (/ &
         'Startup      ' , &
         'Initial      ' , &
         'Measurement  ' , &
         'Hamiltonian  ' , &
         'Evolution    ' , &
         'RNG          ' , &
         'MonteCarlo   ' , &
         'Energy       ' , &
         'Moments      ' , &
         'PrintRestart ' , &
         'LattCorr     ' , &
         'SpinCorr     ' , &
         'ChernNumber  ' , &
         'Dipolar Int. '/)

      !$omp parallel
      !$omp master
      nprocs=omp_get_num_threads()
      !$omp end master
      !$omp end parallel

      if (action.eq.'IN') then ! INIT
         call cpu_time(total0)
         do i=1,ncat
            flops(i)=0.0_dblprec
            timesum(i)=0.0_dblprec
         enddo
         parallel=trim(category).eq.'parallel'
         init=.false.

      else if (action.eq.'RE') then ! RESULT
         if (init.neqv..false.) then
            print *, 'TIMING INITIALIZED BEFORE RESULTS'
            stop
         endif
         ! sum results over all processor
         call cpu_time(total)
         total=total-total0
         timesum(ncat+1)=total
         if (parallel) then
         else
            do i=1,ncat+1
               timemax(i)=timesum(i)
               timemin(i)=timesum(i)
            enddo
         endif
         total=timemax(ncat+1)
         if (iproc.eq.0) then
            write(*,'(1x,a)')&
               '----------------TIME CONSUMPTION REPORT----------------------'
            write(*,*) 'CATEGORY                 TIME(sec)           PERCENT'
            total_pc=0.0_dblprec
            do i=1,ncat
               pc=100.0_dblprec*timemax(i)/real(total,kind=8)
               if(timemax(i)>0.0_dblprec) write(*,'(2x,a14,10x,1pe9.2,10x,0pf8.1 )') &
                  cats(i),timemax(i)/nprocs,pc
               total_pc=total_pc+pc
            enddo
            write(*,'(1x,61("-"))')
            write(*,'(1x,a,10x,1pe9.2,6x,a,0pf5.1)') 'Total wall time',&
               total/nprocs,'Sampled percent ',total_pc
         endif

      else

         ii=100000
         do i=1,ncat
            if (trim(category).eq.trim(cats(i))) then
               ii=i
               exit
            endif
         enddo
         if (ii.eq.100000) then
            print *, 'ACTION  ',action
            stop 'TIMING CATEGORY NOT DEFINED'
         end if
         if (action.eq.'ON') then ! ON
            if (init.neqv..false.) then
               print *, cats(ii),': TIMING INITIALIZED BEFORE READ'
               stop
            endif
            call cpu_time(time0)
            init=.true.
         else if (action.eq.'OF') then ! OFF
            if (init.neqv..true.) then
               print *, cats(ii), 'not initialized'
               stop
            endif
            call cpu_time(time)
            timesum(ii)=timesum(ii)+time-time0
            init=.false.
         else
            stop 'TIMING ACTION UNDEFINED'
         endif

      endif
   end subroutine timing

end module
