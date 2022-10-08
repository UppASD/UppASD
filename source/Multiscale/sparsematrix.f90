!> authors
!> Edgar Mendez
!> Nikos  Ntallis
!> Manuel Pereiro

module SparseMatrix
    use Parameters
    use DynamicArray
    use SortModule
implicit none
    type SpMatrix
        type(DynArrayReal) :: entries
        type(DynArrayInt) :: row
        type(DynArrayInt) :: col
    endtype

public SpMatrix, addMatrixVectorProduct, addTransposedMatrixVectorProduct, &
     allocSpMatrix, deallocSpMatrix, sortMatrixByRow, sortMatrixByRowAndColumn, &
     addMatrixEntry, removeZeros

private
contains

!> Computes res = res + A*vect
subroutine addMatrixVectorProduct(mat, vect, res)
implicit none
    type(SpMatrix), intent(in) :: mat
    real(dblprec), dimension(:), intent(in) :: vect
    real(dblprec), dimension(:), intent(inout) :: res
    
    integer :: i
    do i = 1, mat%entries%length
        res(mat%row%values(i)) = res(mat%row%values(i)) + mat%entries%values(i) * vect(mat%col%values(i))
    end do
end subroutine

!> Computes res = res + A^t * vect
subroutine addTransposedMatrixVectorProduct(mat, vect, res)
implicit none
    type(SpMatrix), intent(in) :: mat
    real(dblprec), dimension(:), intent(in) :: vect
    real(dblprec), dimension(:), intent(inout) :: res
    
    integer :: i
    do i = 1, mat%entries%length
        res(mat%col%values(i)) = res(mat%col%values(i)) + mat%entries%values(i) * vect(mat%row%values(i))
    end do
end subroutine

subroutine addMatrixEntry(matrix, row, col, val)
implicit none
    type(SpMatrix), intent(inout) :: matrix
    integer, intent(in) :: row
    integer, intent(in) :: col
    real(dblprec), intent(in) :: val
    
    call addEntry(matrix%entries, val)
    call addEntry(matrix%row, row)
    call addEntry(matrix%col, col)
end subroutine addMatrixEntry

subroutine removeZeros(matrix)
implicit none
    type(SpMatrix), intent(inout) :: matrix
    
    integer :: currentIndex
    integer :: lastIndex
    
    lastIndex = matrix%entries%length
    
    currentIndex = 1
    do while (currentIndex <= lastIndex)
        if (abs(matrix%entries%values(currentIndex)) < 1d-12) then
            matrix%entries%values(currentIndex) = matrix%entries%values(lastIndex)
            matrix%col%values(currentIndex) = matrix%col%values(lastIndex)
            matrix%row%values(currentIndex) = matrix%row%values(lastIndex)
            lastIndex = lastIndex - 1
        else
            currentIndex = currentIndex + 1
        end if
    end do
    matrix%entries%length = currentIndex - 1
    matrix%row%length = currentIndex - 1
    matrix%col%length = currentIndex - 1
end subroutine removeZeros

subroutine allocSpMatrix(matrix)
implicit none
    type(SpMatrix), intent(inout) :: matrix
    
    call newArray(matrix%entries)
    call newArray(matrix%row)
    call newArray(matrix%col)
end subroutine allocSpMatrix

subroutine deallocSpMatrix(matrix)
implicit none
    type(SpMatrix), intent(inout) :: matrix
    call deallocArray(matrix%entries)
    call deallocArray(matrix%col)
    call deallocArray(matrix%row)
end subroutine deallocSpMatrix

subroutine sortMatrixByRowAndColumn(matrix)
implicit none
    type(SpMatrix), intent(inout) :: matrix

    call sort(compareEntries, swapEntries, 1, matrix%entries%length)
contains
  pure integer function compareEntries(entry1, entry2)
    implicit none
    integer, intent(in) :: entry1
    integer, intent(in) :: entry2
    
    compareEntries = matrix%row%values(entry2) - matrix%row%values(entry1)
    if (compareEntries == 0) then
       compareEntries = matrix%col%values(entry2) - matrix%col%values(entry1)
    end if
  end function compareEntries
    
  subroutine swapEntries(entry1, entry2)
    implicit none
    integer, intent(in) :: entry1
    integer, intent(in) :: entry2
        
    real(dblprec) :: tmpEntry
    integer :: tmpRow
    integer :: tmpCol
    
    tmpEntry = matrix%entries%values(entry1)
        tmpRow = matrix%row%values(entry1)
        tmpCol = matrix%col%values(entry1)
        
        matrix%entries%values(entry1) = matrix%entries%values(entry2)
        matrix%row%values(entry1) = matrix%row%values(entry2)
        matrix%col%values(entry1) = matrix%col%values(entry2)
        
        matrix%entries%values(entry2) = tmpEntry
        matrix%row%values(entry2) = tmpRow
        matrix%col%values(entry2) = tmpCol
      end subroutine swapEntries
  end subroutine sortMatrixByRowAndColumn


subroutine sortMatrixByRow(matrix)
implicit none
    type(SpMatrix), intent(inout) :: matrix

    call sort(compareEntries, swapEntries, 1, matrix%entries%length)
contains
    pure integer function compareEntries(entry1, entry2)
    implicit none
        integer, intent(in) :: entry1
        integer, intent(in) :: entry2
    
        compareEntries = matrix%row%values(entry2) - matrix%row%values(entry1)
    end function compareEntries
    
    subroutine swapEntries(entry1, entry2)
    implicit none
        integer, intent(in) :: entry1
        integer, intent(in) :: entry2
        
        real(dblprec) :: tmpEntry
        integer :: tmpRow
        integer :: tmpCol
        
        tmpEntry = matrix%entries%values(entry1)
        tmpRow = matrix%row%values(entry1)
        tmpCol = matrix%col%values(entry1)
        
        matrix%entries%values(entry1) = matrix%entries%values(entry2)
        matrix%row%values(entry1) = matrix%row%values(entry2)
        matrix%col%values(entry1) = matrix%col%values(entry2)
        
        matrix%entries%values(entry2) = tmpEntry
        matrix%row%values(entry2) = tmpRow
        matrix%col%values(entry2) = tmpCol
    end subroutine swapEntries
end subroutine sortMatrixByRow

end module SparseMatrix
