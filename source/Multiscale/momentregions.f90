!> This module parses the moments fields and allows to
!!  retrieve the magnetic moment affecting a region.
!> authors
!> Edgar Mendez
!> Nikos  Ntallis
!> Manuel Pereiro
module MomentRegions
  use Parameters
  use MultiscaleFileParser
  use FileInput
  use ShapeModule

  implicit none

  type MomentParameters
     real(dblprec) :: magnitude
     real(dblprec), dimension(3) :: direction     
  end type MomentParameters

  !> Single-linked null-ended list of moment descriptiors.
  !! Empty lists are represented by a nullified pointer
  type MomentList
     type(ShapeList)  :: shapes
     type(MomentParameters) :: parameters     
     type(MomentList), pointer :: next
  end type MomentList
    
  private

  public &
       MomentList, parseMomentList,  &
       deallocateMomentList, momentFromPoint
  
  
contains
  
  !> Creates a new moment list with one element.
  !! Use concatMomentLists with this function to make arbitrary lists.
  subroutine newMomentList(list)
    implicit none
    type(MomentList), pointer, intent(out) :: list
    allocate(list)
    nullify(list%next)
    call allocateShapeList(list % shapes)
    list%parameters%magnitude = 1d0 
    list%parameters%direction = (/1,0,0/)
  end subroutine newMomentList

  !> Concatenates two moment lists.
  !! Never deallocate l2, it will be deallocated with the result.
  !! The result is left in l1
  subroutine concatMomentLists(l1,l2)
    implicit none
    type(MomentList), pointer, intent(in) :: l1
    type(MomentList), pointer, intent(in) :: l2

    type(MomentList), pointer :: last_l1
    
    last_l1 => l1
    do while (associated(last_l1%next))
       last_l1=>last_l1%next
    end do

    last_l1%next => l2
  end subroutine concatMomentLists
  
  !> Deallocate a moment list (the whole chain from the given element).
  !! Shapes inside each element are also freed
  !! The pointer is nullified.
  subroutine deallocateMomentList(list)
    implicit none
    type(MomentList), pointer, intent(inout) :: list
    type(MomentList), pointer :: cur
    type(MomentList), pointer :: next

    cur => list 
    do while (associated(cur))       
       next => cur%next
       call deallocateShapeList(cur%shapes)
       deallocate(cur)
       cur => next
    end do
    nullify(list)

  end subroutine deallocateMomentList

  !> Given a point and a moment list, returns a pointer to the
  !! first element in the list that covers the point.
  !! If no element in the list covers the point, a null pointer is
  !! returned instead.
  function momentFromPoint(point,momlist) result (p)
    implicit none
    real(dblprec), dimension(3), intent(in)   :: point
    type(MomentList), intent(in), pointer :: momlist
    type(MomentList), pointer :: p

    p => momlist
    do while (associated(p))
       if(isPointInsideShapeList(p%shapes, point)) return
       p => p%next       
    end do
  end function momentFromPoint


  !> Reads moments parameters from a fileData
  subroutine parseMomentParameters(fData, parameters)
    use FileInput
  implicit none
    type(FileData), intent(inout) :: fData
    type(MomentParameters), intent(out) :: parameters

    real(dblprec) :: norm
    integer :: i
        
    ! Parse Magnitude
    call parseReal(fData,parameters%magnitude)
    if (fData%ierr /= 0) return
    call readNextWord(fData)

    call skipOptionalDelimiter(fData)
    ! Parse vector and normalize
    do i=1,3
       call parseReal(fData,parameters%direction(i))
       if (fData%ierr /= 0) then
          call createErrorMsg(fData,1, 'Unexpected symbol or end of line '//&
               'when reading moment direction vector.')
          return
       end if
       call readNextWord(fData)
    end do
    norm = sqrt(sum(parameters%direction**2))
    parameters%direction = parameters%direction / norm

  end subroutine parseMomentParameters


     
  !> Parses a moment list from a FileData.
  subroutine parseMomentList(fData, momlst)
    use FileInput
  implicit none
    type(FileData), intent(inout) :: fData
    type(MomentList), intent(out), pointer :: momlst
    

    type(MomentList), pointer :: new_momlst
    call newMomentList(new_momlst)

    if (.not. isAtEndOfLine(fData) .and. fData%word /= "file") then
       call parseMomentParameters(fData, new_momlst%parameters)
       call skipOptionalDelimiter(fData)
    else
       call createErrorMsg(fData,1, 'Expected moment parameters')
       return
    end if
    
    if (fData%ierr /= 0) then
       call deallocateMomentList(new_momlst)
       return
    end if

    call parseShapeList(fData, new_momlst%shapes)
    if (associated(momlst)) then
       call concatMomentLists(momlst, new_momlst)
    else
       momlst => new_momlst
    end if
    
  end subroutine parseMomentList

  
end module MomentRegions
