!> authors
!> Edgar Mendez
!> Nikos  Ntallis
!> Manuel Pereiro

module InteractionInput
  use Parameters
  use MultiscaleFileParser
  use FileInput
  use FiniteDifference
implicit none

  type InteractionInfo
     real(dblprec), dimension(:, :), allocatable :: interactionValue
     real(dblprec), dimension(:, :), allocatable :: dirVect
     integer      , dimension(:, :), allocatable :: unitCellDirection
     integer, dimension(:), allocatable :: atomType1
     integer, dimension(:), allocatable :: atomType2
  end type InteractionInfo
  
  type AtomSetupInfo
     integer, dimension(:)   , allocatable :: atomTypes
     integer, dimension(:, :), allocatable :: fromUnitcellLocation
     integer :: nrOfAtoms
  end type AtomSetupInfo

  private
  public &
       InteractionInfo, getInteractionValue, parseInteractionList, &
       initialiseInteractionInfo, deallocateInteractionInfo, printInteractionInfo, &
       getMaxInteractionRadius, &
       AtomSetupInfo, allocateAtomSetupInfo, deallocateAtomSetupInfo

contains

  subroutine initialiseInteractionInfo(info)
    implicit none
    type(InteractionInfo), intent(inout) :: info
  end subroutine initialiseInteractionInfo

  !> Print interaction details
  subroutine printInteractionInfo(info)
    implicit none
    type(InteractionInfo), intent(in) :: info

    integer :: i
    
    if (allocated(info%dirVect)) then
       print *, 'Directional data'
       do i = 1, ubound(info%dirVect, 2)
          print *,&
               'Dir vect: ', info%dirVect(:, i), &
               'Value: ', info%interactionValue(:,i)
       end do
    elseif (allocated(info%unitCellDirection)) then
       print *, 'Unitcell data'
       do i = 1, ubound(info%unitCellDirection, 2)
          print *,&
               'Unitcell vect: ', info%unitCellDirection(:, i),&
               'Value: ', info%interactionValue(:,i)
       end do
    else
       print *, 'No interactions specified.'
    end if
  end subroutine printInteractionInfo


  !> Deallocate an interaction info structure
  subroutine deallocateInteractionInfo(info)
    implicit none
    type(InteractionInfo), intent(inout) :: info
    
    if (allocated(info%dirVect)) deallocate(info%dirVect)
    if (allocated(info%unitCellDirection)) deallocate(info%unitCellDirection)
    if (allocated(info%atomType1)) deallocate(info%atomType1)
    if (allocated(info%atomType2)) deallocate(info%atomType2)
    if (allocated(info%interactionValue)) deallocate(info%interactionValue) 
  end subroutine deallocateInteractionInfo

  !> Parse an atom-to-atom interaction.
  !! This is valid for dm as well as for the exchange, or any other interaction
  !! defined by a scalar or a vector.  
  subroutine parseInteractionList(fData, components, interactionList)
    use FileInput
  implicit none
    type(FileData), intent(inout) :: fData
    integer, intent(in) :: components
    type(InteractionInfo), intent(out) :: interactionList

    procedure(parserFunction), pointer :: listParser
    
    if (trim(fData%word) == 'unitcell') then
       listParser => parseUnitcellList
       call readNextWord(fData)
    else
       listParser => parseList
    end if
    if (isAtEndOfLine(fData)) then
       call readNextLine(fData)
       if (fData%ierr /= 0) return
       call listParser(fData)
    else
       call tryToParseFromExternalFile(fData, listParser)
       if (fData%ierr /= 0) return
    end if
  contains

    subroutine parseList(fData)
      implicit none
      type(FileData), intent(inout) :: fData

      integer :: i, j
      integer :: nrOfLinks

      nrOfLinks = countLinesStartWithInteger(fData)
      allocate(interactionList%atomType1(nrOfLinks))
      allocate(interactionList%atomType2(nrOfLinks))
      allocate(interactionList%dirVect(3, nrOfLinks))
      allocate(interactionList%interactionValue(components, nrOfLinks))
      do i = 1, nrOfLinks
         call parseInt(fData, interactionList%atomType1(i))
         if (fData%ierr /= 0) return
         call readNextWord(fData)
         
         call skipOptionalDelimiter(fData)
         
         call parseInt(fData, interactionList%atomType2(i))
         if (fData%ierr /= 0) return
         call readNextWord(fData)
         
         call skipOptionalDelimiter(fData)
         
         do j = 1, 3
            call parseReal(fData, interactionList%dirVect(j, i))
            if (fData%ierr /= 0) return
            call readNextWord(fData)
         enddo

         call skipOptionalDelimiter(fData)

         do j = 1, components
            call parseReal(fData, interactionList%interactionValue(j,i))
            if (fData%ierr /= 0) return
            if (j /= components) then
               call readNextWord(fData)
            end if
         end do
         
         call readNextLine(fData)
         if (fData%ierr /= 0) return
      enddo
    end subroutine parseList

    subroutine parseUnitcellList(fData)
      implicit none
      type(FileData), intent(inout) :: fData

      integer :: i, j
      integer :: nrOfLinks

      nrOfLinks = countLinesStartWithInteger(fData)
      allocate(interactionList%atomType1(nrOfLinks))
      allocate(interactionList%atomType2(nrOfLinks))
      allocate(interactionList%unitcellDirection(3, nrOfLinks))
      allocate(interactionList%interactionValue(components,nrOfLinks))
      do i = 1, nrOfLinks
         call parseInt(fData, interactionList%atomType1(i))
         if (fData%ierr /= 0) return
         call readNextWord(fData)

         call skipOptionalDelimiter(fData)

         call parseInt(fData, interactionList%atomType2(i))
         if (fData%ierr /= 0) return
         call readNextWord(fData)

         call skipOptionalDelimiter(fData)

         do j = 1, 3
            call parseInt(fData, interactionList%unitcellDirection(j, i))
            if (fData%ierr /= 0) return
            call readNextWord(fData)
         enddo

         call skipOptionalDelimiter(fData)

         do j = 1, components
            call parseReal(fData, interactionList%interactionValue(j,i))
            if (fData%ierr /= 0) return
            if (j /= components) then
               call readNextWord(fData)
            end if
         end do

         call readNextLine(fData)
         if (fData%ierr /= 0) return
      enddo
    end subroutine parseUnitcellList
  end subroutine parseInteractionList





  pure real(dblprec) function getMaxInteractionRadius(excInfo, unitcellSize) result(maxDistance)
    implicit none
    type(InteractionInfo), intent(in) :: excInfo
    real(dblprec), dimension(3), intent(in) :: unitcellSize

    integer :: i
    real(dblprec) :: currentDistance
    real(dblprec) :: unitcellDiagonal

    unitcellDiagonal = sqrt(sum(unitcellSize**2))

    maxDistance = 0.0_dblprec
    if (allocated(excInfo%dirVect)) then
       do i = 1, ubound(excInfo%dirVect, 2)
          currentDistance = sqrt(sum(excInfo%dirVect(:, i)**2))
          maxDistance = max(maxDistance, currentDistance)
       end do
    elseif (allocated(excInfo%unitCellDirection)) then
       do i = 1, ubound(excInfo%unitCellDirection, 2)
          currentDistance = sqrt(real(sum((excInfo%unitcellDirection(:, i) + 1)**2)))
          maxDistance = max(maxDistance, currentDistance)
       end do
    end if
  end function getMaxInteractionRadius

  pure function getInteractionValue ( &
    excInfo, space, maxNrOfUnitcells, &
    atomPositions, atomInfo, atomI, atomJ, tolerance)&
       result(coeff)
    use GeometryModule, only : SpaceStruct, getDirectionalVector
  implicit none
    type(InteractionInfo), intent(in) :: excInfo
    type(SpaceStruct), intent(in) :: space
    integer, dimension(3), intent(in) :: maxNrOfUnitcells
    real(dblprec), dimension(:, :), intent(in) :: atomPositions
    type(AtomSetupInfo), intent(in) :: atomInfo
    real(dblprec),dimension(ubound(excInfo%interactionValue,1)) :: coeff
    
    integer, intent(in) :: atomI
    integer, intent(in) :: atomJ
    real(dblprec), intent(in) :: tolerance

    integer :: i

    coeff = 0d0
    if (allocated(excInfo%dirVect)) then
       do i = 1, ubound(excInfo%dirVect, 2)
          if (atomTypesMatch(i, atomI, atomJ) &
               .and. directionalVectorMatch(i, atomI, atomJ,tolerance)) then
             coeff = excInfo%interactionValue(:,i)
             return
          endif
       end do
    elseif (allocated(excInfo%unitcellDirection)) then
       do i = 1, ubound(excInfo%unitcellDirection, 2)
          if (atomTypesMatch(i, atomI, atomJ) &
             .and. unitcellDirectionMatch(i, atomI, atomJ)) then             
             coeff = excInfo%interactionValue(:,i)
             return             
          end if
       end do
    end if
  contains
    pure logical function atomTypesMatch(i, atomI, atomJ)
      implicit none
      integer, intent(in) :: i, atomI, atomJ

      atomTypesMatch = excInfo%atomType1(i) == atomInfo%atomTypes(atomI) .and. &
           excInfo%atomType2(i) == atomInfo%atomTypes(atomJ)
    end function atomTypesMatch

    pure logical function directionalVectorMatch(i, atomI, atomJ, tolerance)
      implicit none
      integer, intent(in) :: i, atomI, atomJ
      real(dblprec), intent(in) :: tolerance
      
      real(dblprec), dimension(3) :: dirVect

      call getDirectionalVector(space, &
           atomPositions(:, atomI), atomPositions(:, atomJ), dirVect)        
      directionalVectorMatch = (maxval(abs(dirVect - excInfo%dirVect(:, i))) < tolerance)
      
    end function directionalVectorMatch

    pure logical function unitcellDirectionMatch(i, atomI, atomJ)
      implicit none
      integer, intent(in) :: i, atomI, atomJ

      integer, dimension(3) :: location
      logical :: inside
      call modularGrid(maxNrOfUnitcells, space%periodicBoundary, &
           atomInfo%fromUnitcellLocation(:, atomI) + excInfo%unitcellDirection(:, i),&
           location, inside)

      if (.not. inside) then
         unitcellDirectionMatch = .false.
         return
      end if

      unitcellDirectionMatch = all(location == atomInfo%fromUnitcellLocation(:, atomJ))
    end function unitcellDirectionMatch
  end function getInteractionValue



  subroutine allocateAtomSetupInfo(atomSetup, length)
    implicit none
    type(AtomSetupInfo), intent(inout) :: atomSetup
    integer, intent(in) :: length

    allocate(atomSetup%atomTypes(length))
    allocate(atomSetup%fromUnitcellLocation(3, length))
    atomSetup%nrOfAtoms = length
  end subroutine allocateAtomSetupInfo

  subroutine deallocateAtomSetupInfo(atomSetup)
    implicit none
    type(AtomSetupInfo), intent(inout) :: atomSetup

    deallocate(atomSetup%atomTypes)
    deallocate(atomSetup%fromUnitcellLocation)
    atomSetup%nrOfAtoms = 0
  end subroutine deallocateAtomSetupInfo


end module InteractionInput
