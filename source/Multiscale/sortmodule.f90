!> authors
!> Edgar Mendez
!> Nikos  Ntallis
!> Manuel Pereiro

module SortModule
implicit none

  integer, pointer, save :: sortIntArrayCtx(:) => null()
  procedure(compareFunctionIf), pointer, save :: heapsortCompareCtx => null()
  procedure(swapFunctionIf), pointer, save :: heapsortSwapCtx => null()
  integer, save :: heapsortStartIndexCtx = 1

    interface sort
        procedure sortGeneric
        procedure sortIntArray
    end interface
    
    interface
       integer function compareFunctionIf(i, j)
         implicit none
         integer, intent(in) :: i, j
       end function compareFunctionIf
       subroutine swapFunctionIf(i, j)
         implicit none
         integer, intent(in) :: i, j
       end subroutine swapFunctionIf
    end interface

  public searchSortedArray, sort
  public quicksort, heapsort, introsort
private
contains

!> Returns an index of an entry in sortedarray that has the value val
!! If no such value exist the function returns -1
pure integer function searchSortedArray(val, sortedarray)
implicit none
    integer, dimension(:), intent(in) :: sortedarray
    integer, intent(in) :: val
    
    integer :: upperLimit, lowerLimit, searchIndex

    upperLimit = ubound(sortedarray, 1)
    lowerLimit = 1
    searchIndex = (upperLimit + lowerLimit)/2
    if (upperLimit < lowerLimit) then
        searchSortedArray = -1
        return
    endif
    if (val < sortedarray(1) .or. val > sortedarray(upperLimit)) then
        searchSortedArray = -1
        return
    endif
    do while (upperLimit >= lowerLimit)
        if (sortedarray(searchIndex) == val) then
            searchSortedArray = searchIndex
            return
        endif
        if (sortedarray(searchIndex) < val) then
            lowerLimit = searchIndex + 1
        else
            upperLimit = searchIndex - 1
        endif
        searchIndex = (upperLimit + lowerLimit)/2
    enddo
    searchSortedArray = -1
end function searchSortedArray

subroutine sortIntArray(intArray)
implicit none
  integer, dimension(:), target, intent(inout) :: intArray

  sortIntArrayCtx => intArray
  call sortGeneric(sortIntArrayCompare, sortIntArraySwap, 1, ubound(intArray, 1))
  nullify(sortIntArrayCtx)
end subroutine sortIntArray

integer function sortIntArrayCompare(i, j)
implicit none
  integer, intent(in) :: i, j

  sortIntArrayCompare = sortIntArrayCtx(j) - sortIntArrayCtx(i)
end function sortIntArrayCompare

subroutine sortIntArraySwap(i, j)
implicit none
  integer, intent(in) :: i, j
  integer :: tempVal

  tempVal = sortIntArrayCtx(i)
  sortIntArrayCtx(i) = sortIntArrayCtx(j)
  sortIntArrayCtx(j) = tempVal
end subroutine sortIntArraySwap

subroutine sortGeneric(compareFunction, swapFunction, startIndex, endIndex)
implicit none
    interface
        integer function compareFunction(i, j)
        implicit none
            integer, intent(in) :: i, j
        end function
        subroutine swapFunction(i, j)
        implicit none
            integer, intent(in) :: i, j
        end subroutine swapFunction
    end interface
    integer, intent(in) :: startIndex
    integer, intent(in) :: endIndex
    
    call introsort(compareFunction, swapFunction, startIndex, endIndex)
end subroutine sortGeneric

recursive subroutine quicksort(compareFunction, swapFunction, startIndex, endIndex)
implicit none
    procedure(compareFunctionIf) :: compareFunction
    procedure(swapFunctionIf) :: swapFunction
    integer, intent(in) :: startIndex
    integer, intent(in) :: endIndex

    integer :: i, j
    integer :: partitionIndex
    
    if (endIndex - startIndex <= 20) then
        do i = startIndex + 1, endIndex
            do j = i, startIndex + 1, -1
                if (compareFunction(j - 1, j) < 0) then
                    call swapFunction(j - 1, j)
                else
                    exit
                end if
            end do
        end do
    elseif (startIndex < endIndex) then
        call quicksort_partition(compareFunction, swapFunction, startIndex, endIndex, partitionIndex)
        call quicksort(compareFunction, swapFunction, startIndex, partitionIndex - 1)
        call quicksort(compareFunction, swapFunction, partitionIndex + 1, endIndex)
    endif
    
end subroutine quicksort

subroutine quicksort_partition(compare, swap, startIndex, endIndex, p)
  implicit none
  procedure(compareFunctionIf) :: compare
  procedure(swapFunctionIf) :: swap
  integer, intent(in) :: startIndex
  integer, intent(in) :: endIndex
  integer, intent(out) :: p
  
  integer :: pivotIndex
  integer :: i, j
  !call swap(startIndex, startIndex + (endIndex - startIndex)/2)
  pivotIndex = startIndex
  i = startIndex - 1
  j = endIndex + 1
  do while (.true.)
     do while(.true.)
        i = i + 1
        if (compare(i, pivotIndex) <= 0) exit
     enddo
     do while (.true.)
        j = j - 1
        if (compare(j, pivotIndex) >= 0) exit
     enddo
     
     if (i >= j) then
        p = j
        call swap(pivotIndex, p)
        return
     endif
     
     call swap(i, j)
  enddo
end subroutine quicksort_partition
    
    
subroutine heapsort(compareFunction, swapFunction, startIndex, endIndex)
  implicit none
  procedure(compareFunctionIf) :: compareFunction
  procedure(swapFunctionIf) :: swapFunction
  integer, intent(in) :: startIndex
  integer, intent(in) :: endIndex

  integer :: last
    heapsortCompareCtx => compareFunction
    heapsortSwapCtx => swapFunction
    heapsortStartIndexCtx = startIndex

  if(startIndex < endIndex) then
     last = endIndex
     !! Build heap, first is largest (and root)
      call heapsortBuildHeap(1,1+endIndex-startIndex)
     
     ! pull and swap, decrease from the right
     do while(last > startIndex)
        call swapFunction(startIndex,last)
        last = last-1
        call heapsortDownHeap(1,1+last-startIndex)
     end do
  end if

    nullify(heapsortCompareCtx)
    nullify(heapsortSwapCtx)
end subroutine heapsort

  integer function heapsortOffsetCompare(i,j)
  implicit none
    integer, intent(in) :: i,j

    heapsortOffsetCompare = heapsortCompareCtx(i+heapsortStartIndexCtx-1, j+heapsortStartIndexCtx-1)
  end function heapsortOffsetCompare

  subroutine heapsortOffsetSwap(i,j)
  implicit none
    integer, intent(in) :: i,j

    call heapsortSwapCtx(i+heapsortStartIndexCtx-1, j+heapsortStartIndexCtx-1)
  end subroutine heapsortOffsetSwap

  subroutine heapsortBuildHeap(startIndex, endIndex)
  implicit none
    integer, intent(in) :: startIndex
    integer, intent(in) :: endIndex
    integer :: i

    do i = endIndex/2,startIndex,-1
      call heapsortDownHeap(i,endIndex)
    end do
  end subroutine heapsortBuildHeap

  recursive subroutine heapsortDownHeap(startIndex, endIndex)
  implicit none
    integer, intent(in) :: startIndex
    integer, intent(in) :: endIndex

    integer :: max_child
    integer :: child1, child2
    max_child = 0
    child1 = startIndex*2
    if (child1 <= endIndex) then
      child2 = child1+1
      if(child2 <= endIndex) then
        if(heapsortOffsetCompare(child1,child2) > 0) then
          max_child = child2
        else
          max_child = child1
        end if
      else
        max_child = child1
      end if
    end if

    if(max_child /= 0) then
      if(heapsortOffsetCompare(max_child,startIndex)<0) then
        call heapsortOffsetSwap(max_child,startIndex)
        call heapsortDownHeap(max_child,endIndex)
      end if
    end if
  end subroutine heapsortDownHeap


subroutine introsort(compareFunction, swapFunction, startIndex, endIndex)
  implicit none
    procedure(compareFunctionIf) :: compareFunction
    procedure(swapFunctionIf) :: swapFunction
    integer, intent(in) :: startIndex
    integer, intent(in) :: endIndex

    real,parameter :: inv_log_2 = 1.0/log(2.0)
    
    integer :: depthLimit

    if(endIndex /= startIndex) then
       depthLimit = floor(log(max(1.0,real(endIndex-startIndex)))*inv_log_2)
       call introsort_aux(compareFunction, swapFunction, startIndex,endIndex, 0)
    end if
    
  contains

    recursive subroutine introsort_aux(compareFunction, swapFunction, &
         startIndex, endIndex, depth)
      implicit none
      procedure(compareFunctionIf) :: compareFunction
      procedure(swapFunctionIf) :: swapFunction
      integer, intent(in) :: startIndex
      integer, intent(in) :: endIndex
      integer, intent(in) :: depth
      
      integer :: partitionIndex

      if (depth > depthLimit) then
         call heapsort(compareFunction,swapFunction,startIndex,endIndex)
      else
         if (startIndex < endIndex) then
            call quicksort_partition(compareFunction, swapFunction, startIndex, endIndex, partitionIndex)
            call introsort_aux(compareFunction, swapFunction, startIndex, partitionIndex - 1, depth+1)
            call introsort_aux(compareFunction, swapFunction, partitionIndex + 1, endIndex, depth+1)
         endif
      end if
    end subroutine introsort_aux    
  end subroutine introsort
end module SortModule
