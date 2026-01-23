!> authors
!> Edgar Mendez
!> Nikos  Ntallis
!> Manuel Pereiro

module DynamicArray
    use Parameters
implicit none
    integer, parameter :: initialSize = 50
    real, parameter :: growCoeff = 1.25
    type DynArrayReal
        real(dblprec), dimension(:), pointer :: values
        integer :: length
    endtype
    type DynArrayInt
        integer, dimension(:), pointer :: values
        integer :: length
    endtype
    
    !> Initializes a dynamic array.
    !! @param array (out) the array to initialize
    !! @param size  (opt in) the number of elements to preallocate
    interface newArray
        procedure newIntArray
        procedure newRealArray
    end interface newArray

    !> Releases a dynamic array´s memory, effectively leaving it empty.
    !! The array can still be used.
    !! @param array (inout) the array to release
    interface deallocArray
        procedure deallocIntArray
        procedure deallocRealArray
    end interface

    !> Appends one element to the tail of a given array.
    !! If the array is too small to hold the element,
    !! it will grow to fit growCoeff*length elements.
    !! @param array Where to put the element
    !! @param value Value to append
    interface addEntry
        procedure addIntEntry
        procedure addRealEntry
    end interface
    
    !> Erases all elements from the array but does not release
    !! the memory used to hold them.
    !! @param array  Array
    interface clearArray
        procedure clearIntArray
        procedure clearRealArray
    end interface clearArray

    
    !> Ensures that the specified number of elements fits in the array.
    !! The array is expanded as needed to fit all the elements.
    !! @param array  Array
    !! @param length Number of elements that must fit in the array.
    interface ensureAllocLength
        procedure ensureIntAllocLength
        procedure ensureRealAllocLength
    end interface ensureAllocLength
    
public newArray, addEntry, clearArray, ensureAllocLength, DynArrayReal, &
       DynArrayInt, deallocArray

private
contains

subroutine newIntArray(array, size)
implicit none
    type(DynArrayInt), intent(out) :: array
    Integer, optional, intent(in) :: size
    
    nullify(array%values)
    array%length = 0
    
    if(present(size)) then
        if (size > 0) then
            call ensureAllocLength(array,size)
        endif
    end if
end subroutine newIntArray

subroutine newRealArray(array, size)
implicit none
    type(DynArrayReal), intent(out) :: array
    integer, optional, intent(in) :: size
    
    nullify(array%values)
    array%length = 0
    
    if(present(size)) then
        if (size > 0) then
            call ensureAllocLength(array,size)
        end if
    end if
end subroutine newRealArray

subroutine deallocIntArray(array)
implicit none
    type(DynArrayInt), intent(inout) :: array
    
    array%length = 0
    if (associated(array%values)) then
        deallocate(array%values)
        nullify(array%values)
    endif
end subroutine

subroutine deallocRealArray(array)
implicit none
    type(DynArrayReal), intent(inout) :: array
    
    array%length = 0
    if (associated(array%values)) then
        deallocate(array%values)
        nullify(array%values)
    endif
end subroutine

subroutine addIntEntry(array, val)
implicit none
    type(DynArrayInt), intent(inout) :: array
    integer, intent(in) :: val
    
    integer :: allocated

    if (.not. associated(array%values)) then
       allocated = -1
       array%length = 0
    else
       allocated = ubound(array%values,1)
    endif

    if (allocated < array%length + 1) then 
       call ensureIntAllocLength(array, &
            max(initialSize, ceiling(array%length*growCoeff)))
    endif
    array%length = array%length + 1
    array%values(array%length) = val
    
end subroutine addIntEntry

subroutine addRealEntry(array, val)
implicit none
    type(DynArrayReal), intent(inout) :: array
    real(dblprec), intent(in) :: val

    integer :: allocated
    
    if (.not. associated(array%values)) then
       allocated = -1
       array%length = 0
    else
       allocated = ubound(array%values,1)
    endif

    if (allocated < array%length + 1) then 
       call ensureRealAllocLength(array, &
            max(initialSize, ceiling(array%length*growCoeff)))
    endif
    array%length = array%length + 1
    array%values(array%length) = val
end subroutine addRealEntry
  
pure subroutine clearRealArray(array)
implicit none
    type(DynArrayReal), intent(inout) :: array
    array%length = 0
end subroutine clearRealArray

pure subroutine clearIntArray(array)
implicit none
    type(DynArrayInt), intent(inout) :: array
    array%length = 0
end subroutine clearIntArray

subroutine ensureRealAllocLength(array, length)
implicit none
    type(DynArrayReal), intent(inout) :: array
    integer, intent(in) :: length
    
    integer :: allocatedSize
    real(dblprec), dimension(:), pointer :: tmpValues

    if (.not. associated(array%values)) then
        allocate(array%values(length))
        array%length = 0
    else
       allocatedSize = ubound(array%values, 1)
       if (allocatedSize < length) then
          tmpValues => array%values
          allocate(array%values(length))
          array%values(1:allocatedSize) = tmpValues
          deallocate(tmpValues)
       endif
    end if
end subroutine

subroutine ensureIntAllocLength(array, length)
implicit none
    type(DynArrayInt), intent(inout) :: array
    integer, intent(in) :: length
    
    integer :: allocatedSize
    integer, dimension(:), pointer :: tmpValues

    if (.not. associated(array%values)) then
        allocate(array%values(length))
        array%length = 0
    else
       allocatedSize = ubound(array%values, 1)
       if (allocatedSize < length) then
          tmpValues => array%values
          allocate(array%values(length))
          array%values(1:allocatedSize) = tmpValues
          deallocate(tmpValues)
       endif
    endif
end subroutine

end module
