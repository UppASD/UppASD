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

   type(InteractionInfo), pointer, save :: parseInteractionListCtx => null()
   integer, save :: parseInteractionComponentsCtx = 0

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
      type(InteractionInfo), target, intent(out) :: interactionList

    procedure(parserFunction), pointer :: listParser

    parseInteractionListCtx => interactionList
    parseInteractionComponentsCtx = components

    if (trim(fData%word) == 'unitcell') then
       listParser => parseUnitcellListCallback
       call readNextWord(fData)
    else
       listParser => parseListCallback
    end if
    if (isAtEndOfLine(fData)) then
       call readNextLine(fData)
       if (fData%ierr /= 0) return
       call listParser(fData)
    else
       call tryToParseFromExternalFile(fData, listParser)
       if (fData%ierr /= 0) return
    end if

    nullify(parseInteractionListCtx)
    parseInteractionComponentsCtx = 0
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
    integer, intent(in) :: atomI
    integer, intent(in) :: atomJ
    real(dblprec), intent(in) :: tolerance
    real(dblprec),dimension(ubound(excInfo%interactionValue,1)) :: coeff

    integer :: i
    real(dblprec), dimension(3) :: dirVect
    integer, dimension(3) :: location
    logical :: inside

    coeff = 0d0
    if (allocated(excInfo%dirVect)) then
       do i = 1, ubound(excInfo%dirVect, 2)
          if (excInfo%atomType1(i) == atomInfo%atomTypes(atomI) .and. &
               excInfo%atomType2(i) == atomInfo%atomTypes(atomJ)) then
             call getDirectionalVector(space, atomPositions(:, atomI), atomPositions(:, atomJ), dirVect)
             if (maxval(abs(dirVect - excInfo%dirVect(:, i))) < tolerance) then
                coeff = excInfo%interactionValue(:,i)
                return
             end if
          endif
       end do
    elseif (allocated(excInfo%unitCellDirection)) then
       do i = 1, ubound(excInfo%unitCellDirection, 2)
          if (excInfo%atomType1(i) == atomInfo%atomTypes(atomI) .and. &
              excInfo%atomType2(i) == atomInfo%atomTypes(atomJ)) then
             call modularGrid(maxNrOfUnitcells, space%periodicBoundary, &
                  atomInfo%fromUnitcellLocation(:, atomI) + excInfo%unitcellDirection(:, i), &
                  location, inside)

             if (.not. inside) cycle
             if (all(location == atomInfo%fromUnitcellLocation(:, atomJ))) then
                coeff = excInfo%interactionValue(:,i)
                return
             end if
          end if
       end do
    end if
  end function getInteractionValue

  subroutine parseListCallback(fData)
    implicit none
    type(FileData), intent(inout) :: fData

    integer :: i, j
    integer :: nrOfLinks

    nrOfLinks = countLinesStartWithInteger(fData)
    allocate(parseInteractionListCtx%atomType1(nrOfLinks))
    allocate(parseInteractionListCtx%atomType2(nrOfLinks))
    allocate(parseInteractionListCtx%dirVect(3, nrOfLinks))
    allocate(parseInteractionListCtx%interactionValue(parseInteractionComponentsCtx, nrOfLinks))
    do i = 1, nrOfLinks
       call parseInt(fData, parseInteractionListCtx%atomType1(i))
       if (fData%ierr /= 0) return
       call readNextWord(fData)

       call skipOptionalDelimiter(fData)

       call parseInt(fData, parseInteractionListCtx%atomType2(i))
       if (fData%ierr /= 0) return
       call readNextWord(fData)

       call skipOptionalDelimiter(fData)

       do j = 1, 3
          call parseReal(fData, parseInteractionListCtx%dirVect(j, i))
          if (fData%ierr /= 0) return
          call readNextWord(fData)
       enddo

       call skipOptionalDelimiter(fData)

       do j = 1, parseInteractionComponentsCtx
          call parseReal(fData, parseInteractionListCtx%interactionValue(j,i))
          if (fData%ierr /= 0) return
          if (j /= parseInteractionComponentsCtx) then
             call readNextWord(fData)
          end if
       end do

       call readNextLine(fData)
       if (fData%ierr /= 0) return
    enddo
  end subroutine parseListCallback

  subroutine parseUnitcellListCallback(fData)
    implicit none
    type(FileData), intent(inout) :: fData

    integer :: i, j
    integer :: nrOfLinks

    nrOfLinks = countLinesStartWithInteger(fData)
    allocate(parseInteractionListCtx%atomType1(nrOfLinks))
    allocate(parseInteractionListCtx%atomType2(nrOfLinks))
    allocate(parseInteractionListCtx%unitcellDirection(3, nrOfLinks))
    allocate(parseInteractionListCtx%interactionValue(parseInteractionComponentsCtx,nrOfLinks))
    do i = 1, nrOfLinks
       call parseInt(fData, parseInteractionListCtx%atomType1(i))
       if (fData%ierr /= 0) return
       call readNextWord(fData)

       call skipOptionalDelimiter(fData)

       call parseInt(fData, parseInteractionListCtx%atomType2(i))
       if (fData%ierr /= 0) return
       call readNextWord(fData)

       call skipOptionalDelimiter(fData)

       do j = 1, 3
          call parseInt(fData, parseInteractionListCtx%unitcellDirection(j, i))
          if (fData%ierr /= 0) return
          call readNextWord(fData)
       enddo

       call skipOptionalDelimiter(fData)

       do j = 1, parseInteractionComponentsCtx
          call parseReal(fData, parseInteractionListCtx%interactionValue(j,i))
          if (fData%ierr /= 0) return
          if (j /= parseInteractionComponentsCtx) then
             call readNextWord(fData)
          end if
       end do

       call readNextLine(fData)
       if (fData%ierr /= 0) return
    enddo
  end subroutine parseUnitcellListCallback



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
