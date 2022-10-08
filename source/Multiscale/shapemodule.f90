!> authors
!> Edgar Mendez
!> Nikos  Ntallis
!> Manuel Pereiro

module ShapeModule
    use Parameters
implicit none
    type SphereShape
        real(dblprec), dimension(3) :: center
        real(dblprec) :: radius
    end type
    type BoxShape
        real(dblprec), dimension(3) :: corner
        real(dblprec), dimension(3) :: sizes
    end type
    type GenericShape
        type(BoxShape), pointer :: box
        type(SphereShape), pointer :: sphere
        type(GenericShape), pointer :: nextShape !In order for them to be in a linked list
    end type
    type ShapeList
        type(GenericShape), pointer :: dummyBeforeFirst
        type(GenericShape), pointer :: lastElement
    end type
    
    !> Checks if a point is inside a shape
    interface isPointInsideShape
        procedure isPointInsideSphere
        procedure isPointInsideBox
        procedure isPointInsideGenericShape
        procedure isPointInsideShapeList
    end interface
    
private
        
public deallocateShape, isPointInsideShape, deallocateShapeList, &
       parseShapeList, allocateShapeList, BoxShape, SphereShape,&
       ShapeList, printShapeList, isPointInsideShapeList
contains

subroutine parseShapeList(fData, list)
    use MultiscaleFileParser
    use FileInput
implicit none
    type(FileData), intent(inout) :: fData
    type(ShapeList), intent(inout) :: list
    
    if (.not. isAtEndOfLine(fData)) then
        call tryToParseFromExternalFile(fData, parseList)
        if (fData%ierr /= 0) return
    else
        call readNextLine(fData)
        if (fData%ierr /= 0) return
        call parseList(fData)
    end if
  contains
    subroutine parseList(fData)
    implicit none
        type(FileData), intent(inout) :: fData
    
        do while (.true.)
            call toLowerCase(fData%word)
            select case (trim(fData%word))
            case('sphere')
                call appendNewShape(list)
                call readNextWord(fData)
                call parseSphere(fData, list%lastElement)
                if (fData%ierr /= 0) return
            case('box')
                call appendNewShape(list)
                call readNextWord(fData)
                call parseBox(fData, list%lastElement)
                if (fData%ierr /= 0) return
            case default
                exit
            end select
            if (fData%ierr /= 0) exit
        enddo
    end subroutine
end subroutine parseShapeList

subroutine parseSphere(fData, shape)
    use MultiscaleFileParser
implicit none
    type(FileData), intent(inout) :: fData
    type(GenericShape), intent(out) :: shape
    
    integer :: i
    allocate(shape%sphere)
    do i = 1, 3
        call parseReal(fData, shape%sphere%center(i))
        if (fData%ierr /= 0) return
        call readNextWord(fData)
    enddo
    call skipOptionalDelimiter(fData) 
    call parseReal(fData, shape%sphere%radius)
    if (fData%ierr /= 0) return
    call readNextLine(fData)
end subroutine parseSphere

subroutine parseBox(fData, shape)
    use MultiscaleFileParser
implicit none
    type(FileData), intent(inout) :: fData
    type(GenericShape), intent(out) :: shape
    
    integer :: i

    allocate(shape%box)
    do i = 1, 3
        call parseReal(fData, shape%box%corner(i))
        if (fData%ierr /= 0) return
        call readNextWord(fData)
    enddo
    call skipOptionalDelimiter(fData) 
    do i = 1, 3
        call parseReal(fData, shape%box%sizes(i))
        if (fData%ierr /= 0) return
        call readNextWord(fData)
    enddo
    call readNextLine(fData)
end subroutine parseBox

subroutine printShapeList(list)
implicit none
    type(ShapeList), intent(in) :: list
    
    type(GenericShape), pointer :: currentShape
    
    currentShape => list%dummyBeforeFirst%nextShape
    do while (associated(currentShape))
        call printShape(currentShape)
        currentShape => currentShape%nextShape
    end do
end subroutine

subroutine printShape(genShape)
implicit none
    type(GenericShape), intent(in) :: genShape
    
    if (associated(genShape%sphere)) then
        print *, 'Sphere Radius: ', genShape%sphere%radius, ' Center: ', genShape%sphere%center
    elseif (associated(genShape%box)) then
        print *, 'Box Corner: ', genShape%box%corner, ' Size: ', genShape%box%sizes
    else
        print *, 'The shape is not associated with anything, this is probably a bug.'
    end if

end subroutine

subroutine allocateGenericShape(genShape)
implicit none
    type(GenericShape), intent(inout) :: genShape
    
    nullify(genShape%box)
    nullify(genShape%sphere)
    nullify(genShape%nextShape)
end subroutine

subroutine deallocateShape(genShape)
implicit none
    type(GenericShape), intent(inout) :: genShape
    if (associated(genShape%box)) then
        deallocate(genShape%box)
    elseif (associated(genShape%sphere)) then
        deallocate(genShape%sphere)
    endif
    nullify(genShape%nextShape)
end subroutine deallocateShape

subroutine allocateShapeList(list)
implicit none
    type(ShapeList), intent(inout) :: list
    
    allocate(list%dummyBeforeFirst)
    call allocateGenericShape(list%dummyBeforeFirst)
    list%lastElement => list%dummyBeforeFirst
end subroutine

subroutine deallocateShapeList(list)
implicit none
    type(ShapeList), intent(inout) :: list
    
    type(GenericShape), pointer :: currentShape
    type(GenericShape), pointer :: nextShape
    
    currentShape => list%dummyBeforeFirst
    do while (associated(currentShape))
        nextShape => currentShape%nextShape
        call deallocateShape(currentShape)
        deallocate(currentShape)
        currentShape => nextShape
    end do
    nullify(list%dummyBeforeFirst)
    nullify(list%lastElement)
end subroutine

!> Appends a new shape to the end of list.
!! @param[in,out] list The list that will be appended to.
subroutine appendNewShape(list)
implicit none
    type(ShapeList), intent(inout) :: list
    
    allocate(list%lastElement%nextShape)
    list%lastElement => list%lastElement%nextShape
    call allocateGenericShape(list%lastElement)
end subroutine

!> Check if a point is inside a generic shape
!! @param[in] s The generic shape
!! @param[in] point The point.
logical function isPointInsideGenericShape(s, point)
implicit none
    type(GenericShape), intent(in) :: s
    real(dblprec), dimension(3), intent(in) :: point

    if (associated(s%box)) then
        isPointInsideGenericShape = isPointInsideBox(s%box, point)
    elseif (associated(s%sphere)) then
        isPointInsideGenericShape = isPointInsideSphere(s%sphere, point)
    else
        print *, 'ERROR:'
        print *, 'The generic shape is not associated with anything'
        stop
    endif
end function isPointInsideGenericShape

!> Check if a point is inside one of the shapes in list
!! @param[in] list
!! @param[in] point
logical function isPointInsideShapeList(list, point) result(inside)
implicit none
    type(ShapeList), intent(in) :: list
    real(dblprec), dimension(3), intent(in) :: point
    
    type(GenericShape), pointer :: currentShape
    
    currentShape => list%dummyBeforeFirst%nextShape
    do while (associated(currentShape))
        if (isPointInsideShape(currentShape, point)) then
            inside = .true.
            return
        end if
        currentShape => currentShape%nextShape
    end do
    inside = .false.
end function

logical function isPointInsideSphere(s, point)
implicit none
    type(SphereShape), intent(in) :: s
    real(dblprec), dimension(3), intent(in) :: point
    
    isPointInsideSphere = sum((s%center - point)**2) < s%radius**2
end function isPointInsideSphere

logical function isPointInsideBox(b, point)
implicit none
    type(BoxShape), intent(in) :: b
    real(dblprec), dimension(3), intent(in) :: point
    
    isPointInsideBox = (b%corner(1) <= point(1) .and. b%corner(1) + b%sizes(1) >= point(1)) &
                    .and. (b%corner(2) <= point(2) .and. b%corner(2) + b%sizes(2) >= point(2)) &
                    .and. (b%corner(3) <= point(3) .and. b%corner(3) + b%sizes(3) >= point(3))
end function isPointInsideBox

end module ShapeModule
