!> authors
!> Edgar Mendez
!> Nikos  Ntallis
!> Manuel Pereiro

module AtomGenerator
    use Parameters
    use GeometryModule
    use ShapeModule
    use FiniteDifference
    use KdTreeModule
implicit none

    type AtomCell
        real(dblprec), dimension(:, :), allocatable :: positions
        real(dblprec), dimension(:, :), allocatable :: moments
        real(dblprec), dimension(:), allocatable :: momentMagnitudes
        real(dblprec), dimension(3) :: size
        integer :: nrOfAtoms
        integer, dimension(:), allocatable :: chemicalTypes
        integer, dimension(:), allocatable :: zones
        integer, dimension(:, :), allocatable :: unitcellLocation
    end type

    type AtomLinkedList
        type(AtomLinkedListNode), pointer :: dummyElementBeforeFirst
        type(AtomLinkedListNode), pointer :: lastElement
    end type

    type AtomLinkedListNode
        type(AtomLinkedListNode), pointer :: nextAtom
        real(dblprec), dimension(3) :: positions
        real(dblprec), dimension(3) :: moment
        integer, dimension(3) :: unitcellLocation
        integer :: atomType
        integer :: zone
        integer :: fromAtomInUnitCell
    end type

    type AtomZoneList
       type(AtomZoneList), pointer :: next
       type(ShapeList) :: shape
       integer :: zone_id
    end type AtomZoneList
    
public AtomZoneList, deallocateZoneList, parseZoneList, &         
       generateAtoms, generatePaddingAtoms, &
       allocateAtomLinkedList, extractPositions, extractTypes, &
       extractFromUnitcellLocation, generateFiniteDifferenceAtoms, &
       extractMoments, deallocateAtomList, &
       allocateAtomCell, deallocateAtomCell, &
       AtomCell, AtomLinkedList, parseAtomCell

private 

contains

!< Releases an AtomZoneList
recursive subroutine deallocateZoneList(l)
implicit none
  type(AtomZoneList), pointer, intent(inout) :: l
  type(AtomZoneList), pointer :: next

  if(associated(l)) then
     next => l%next
     call deallocateZoneList(next)
     deallocate(l)     
  end if
  nullify(l)
end subroutine deallocateZoneList

!< inserts a new element at the beginning of the given AtomZoneList
subroutine insertZone(shape, zone_id, list)
  implicit none
  type(ShapeList), intent(in) :: shape
  integer, intent(in) :: zone_id
  type(AtomZoneList), pointer, intent(inout) :: list
  type(AtomZoneList), pointer :: next

  next => list
  allocate(list)
  list%next => next
  list%zone_id = zone_id
  list%shape = shape
  
end subroutine insertZone

  !> Reads anisotropy parameters from a fileData
  subroutine parseZoneList(fData, atomZones)
    use MultiscaleFileParser
    use FileInput
  implicit none
    type(FileData), intent(inout) :: fData
    type(AtomZoneList), pointer, intent(out) :: atomZones

    integer :: zone
    type(ShapeList) :: shapes
    call allocateShapeList(shapes)

    if (.not. isAtEndOfLine(fData) .and. fData%word /= "file") then
       ! Parse atom zone id
       call parseInt(fData,zone)
       if (fData%ierr /= 0) then
          call createErrorMsg(fData, 1, 'Expected zone identifier (an integer)')
          return
       end if
       call readNextWord(fData)
       if (zone < 0) then
          call createErrorMsg(fData, 1, 'Zone id must be positive.')
          return
       end if
    end if

    if (fData%ierr /= 0) then
       call deallocateShapeList(shapes)
       return
    end if

    call parseShapeList(fData, shapes)
    call insertZone(shapes, zone, atomZones)
    
  end subroutine parseZoneList

!< Returns the id of the first given zone that contains the point.
!! Returns 0 if no zone is found
recursive function atomZoneFromPoint(zones, point) result (zone)
  use ShapeModule, only: ShapeList, isPointInsideShapeList
implicit none
  type(AtomZoneList), pointer, intent(in) :: zones
  real(dblprec), dimension(3), intent(in) :: point
  integer :: zone

  zone = 0
  if (associated(zones)) then
     if (isPointInsideShapeList(zones%shape, point)) then
        zone = zones%zone_id
     else
        zone = atomZoneFromPoint(zones%next, point)
     end if
  endif
end function atomZoneFromPoint


function countAtomsInCell(fData) result(natoms)
  use MultiscaleFileParser
  use FileInput
implicit none
  type(FileData), intent(inout) :: fData
  integer :: natoms
  natoms = 0
  call peekFile(fData, countAtoms)
contains
  subroutine countAtoms(fData)
    use MultiscaleFileParser
  implicit none
    type(FileData), intent(inout) :: fData
    integer :: testAtomType
    integer :: errVal
    
    do while (.true.)
       read(fData%word, *, iostat=errVal) testAtomType
       if (errVal /= 0) then
          if (trim(fData%word) .ne. 'zone') then
             exit
          end if
       else
          natoms = natoms + 1
       end if
       
       call skipToNextLine(fData)
       if (fData%ierr /= 0) then
          fData%ierr = 0
          exit
       endif

    enddo
  end subroutine countAtoms
end function countAtomsInCell
  


!< Reads an atomic unit cell specification
subroutine parseAtomCell(fData, cell)
    use MultiscaleFileParser
    use FileInput
implicit none
    type(FileData), intent(inout) :: fData
    type(AtomCell), intent(inout) :: cell
    
    integer :: i, j
    real(dblprec) :: tmpMag
    integer :: nrOfAtoms

    do i = 1, 3
        call parseReal(fData, cell%size(i))
        if (fData%ierr /= 0) return
        call readNextWord(fData)
    end do
    if (isAtEndOfLine(fData)) then
        call readNextLine(fData)
        if (fData%ierr /= 0) return
        call atomListParser(fData)
    else
        call tryToParseFromExternalFile(fData, atomListParser)
    end if

  contains
    subroutine atomListParser(fData)
    implicit none
        type(FileData), intent(inout) :: fData
        integer :: zone, i
        zone = 0
        nrOfAtoms = countAtomsInCell(fData)
        
        call allocateAtomCell(cell, nrOfAtoms)
        
        i = 1
        do while(i <= nrOfAtoms)
           if (trim(fData%word) == 'zone') then
              call readNextWord(fData)
              call parseInt(fData, zone)
              if (zone < 0) then
                 call createErrorMsg(fData, 1, 'The zone identifier must be positive')
                 return 
              end if
           else
              call parseInt(fData, cell%chemicalTypes(i))
              if (fData%ierr /= 0) return
              call readNextWord(fData)

              call skipOptionalDelimiter(fData)
            
              do j = 1, 3
                 call parseReal(fData, cell%positions(j, i))
                 if (fData%ierr /= 0) return
                 call readNextWord(fData)
              enddo
              
              call skipOptionalDelimiter(fData)
              
              call parseReal(fData, cell%momentMagnitudes(i))
              if (fData%ierr /= 0) return
              call readNextWord(fData)
              
              call skipOptionalDelimiter(fData)
              
              do j = 1, 3
                 call parseReal(fData, cell%moments(j, i))
                 if (fData%ierr /= 0) return
                 call readNextWord(fData)
              enddo
              tmpMag = sqrt(sum(cell%moments(:, i)**2))
              if (tmpMag < 1d-3) then
                 call createErrorMsg(fData, 1, 'The directional moment have too small magnitude.')
                 return
              endif
              cell%moments(:, i) = cell%moments(:, i)/tmpMag
              cell%zones(i) = zone
              i = i + 1
           end if

           call readNextLine(fData)
           if (fData%ierr /= 0) return              

        enddo
      end subroutine atomListParser

        
        
      
end subroutine parseAtomCell

subroutine allocateAtomCell(atoms, nrOfAtoms)
implicit none
    type(AtomCell), intent(inout) :: atoms
    integer, intent(in) :: nrOfAtoms
    
    allocate(atoms%chemicalTypes(nrOfAtoms))
    allocate(atoms%positions(3, nrOfAtoms))
    allocate(atoms%moments(3, nrOfAtoms))
    allocate(atoms%momentMagnitudes(nrOfAtoms))
    allocate(atoms%unitcellLocation(3, nrOfAtoms))
    allocate(atoms%zones(nrOfAtoms))
    atoms%nrOfAtoms = nrOfAtoms
end subroutine allocateAtomCell

subroutine deallocateAtomCell(atoms)
implicit none
    type(AtomCell), intent(inout) :: atoms
    
    deallocate(atoms%chemicalTypes)
    deallocate(atoms%positions)
    deallocate(atoms%momentMagnitudes)
    deallocate(atoms%moments)
    deallocate(atoms%unitcellLocation)
    deallocate(atoms%zones)
    atoms%nrOfAtoms = 0
end subroutine deallocateAtomCell

!> Generates all the atoms from the unitCell. The unitCell is
!! repeated over the axis aligned box that have one corner in (0, 0, 0)
!! and the other in mesh%universeSize. It will only add atoms
!! in the regions which are atomistic cubes in the finite difference mesh.
!! It will not add atoms if the atom will be placed inside a holeShape.
!! This subroutine handles 1, 2 and 3 dimensions. Which dimension is choosen is
!! from mesh%spatDimension.
!! @param[in] mesh The finite difference mesh.
!! @param[in] unitCell The unit cell that will be repeated.
!! @param[in] holeShapes Shapes that defines the holes.
!! @param[out] atoms The atoms that are generated.
!! @param[out] nrOfAtoms The number of atoms that were generated.
subroutine generateAtoms(mesh, unitCell, holeShapes, atomZones, atoms, nrOfAtoms)
implicit none
    type(FiniteDiffMesh), intent(in) :: mesh
    type(AtomCell), intent(in) :: unitCell
    type(ShapeList), intent(in) :: holeShapes
    type(AtomZoneList), pointer, intent(in) :: atomZones
    type(AtomLinkedList), intent(inout) :: atoms
    integer, intent(out) :: nrOfAtoms
    
    call generatePositions(mesh%space, unitCell, isPointInAtomZone, nrOfAtoms, atoms)
    call initializeMoments(atoms, unitCell)
    call initializeAtomType(atoms, unitCell)
contains
    logical function isPointInAtomZone(point, zone)
    implicit none
        real(dblprec), dimension(3), intent(in) :: point
        integer, intent(in) :: zone
        isPointInAtomZone = &
             isPointInsideAtomisticBox(mesh, point) &
             .and. .not. isPointInsideShape(holeShapes, point) &
             .and. atomZoneFromPoint(atomZones, point) == zone
        
    end function isPointInAtomZone
end subroutine generateAtoms

!> Generates the padding atoms. The unit cell is repeated over the universe.
!! Those atoms that then are outside the atomisitc region but close enough to it, is added as a padding atom.
!! The distance to the atomistic region is definied as the distance to the closest atom.
!! @param[in] mesh A finite difference mesh that tells where the atomistic/non atomistic part of the space is.
!! @param[in] unitCell The unitcell that will be repeated.
!! @param[in] atomPositions The positions of the atoms
!! @param[in] atomTree A tree that contains all the real atoms
!! @param[in] paddingWidth The width of the padding band.
!! @param[out] paddingAtoms The generated padding atoms.
!! @param[out] nrOfAtoms The number of atoms that were generated.
subroutine generatePaddingAtoms(mesh, unitCell, atomPositions, atomTree, paddingWidth, paddingAtoms, nrOfAtoms)
implicit none
    type(FiniteDiffMesh), intent(in) :: mesh
    type(AtomCell), intent(in) :: unitCell
    real(dblprec), dimension(:, :), intent(in) :: atomPositions
    type(KdTree), intent(in) :: atomTree
    real(dblprec), intent(in) :: paddingWidth
    type(AtomLinkedList), intent(inout) :: paddingAtoms
    integer, intent(out) :: nrOfAtoms
    
    call generatePositions(mesh%space, unitCell, isInsidePaddingRegion, nrOfAtoms, paddingAtoms)
    call initializeMoments(paddingAtoms, unitcell)
    call initializeAtomType(paddingAtoms, unitcell)

  contains
    
    logical function isInsidePaddingRegion(point, zone)
    implicit none
        real(dblprec), dimension(3), intent(in) :: point
        integer, intent(in) :: zone
        
        real(dblprec) :: dist

        if (isPointInsideAtomisticBox(mesh, point) .or. getBoundaryDistanceAt(mesh, point) > paddingWidth) then
            isInsidePaddingRegion = .false.
        else
            call distToClosest(mesh%space, atomPositions, atomTree, point, distance = dist)
            isInsidePaddingRegion = dist < paddingWidth .and. zone == 0
        end if
    end function isInsidePaddingRegion
end subroutine generatePaddingAtoms

!< Generates positions for the atoms by repeating the unit cell.
!! @param[in] space
!! @param[in] unitCell
!! @param[in] filter A function that filters the positions, if the filter return true then the atom will be added, otherwise it will not.
!! @param[out] nrOfAtomsGenerated The number of positions generated
!! @param[out] generatedAtoms The generated atoms
subroutine generatePositions(space, unitCell, filter, nrOfAtomsGenerated, generatedAtoms)
implicit none
    type(SpaceStruct), intent(in) :: space
    type(AtomCell), intent(in) :: unitCell
    interface
       logical function filter(point, zone)
         use Parameters
         implicit none
         real(dblprec), dimension(3), intent(in) :: point
         integer, intent(in) :: zone
       end function filter
    end interface
    integer, intent(out) :: nrOfAtomsGenerated
    type(AtomLinkedList), intent(inout) :: generatedAtoms
    
    integer :: atomInUnitCell
    integer, dimension(3) :: unitCellLocation
    real(dblprec), dimension(3) :: atomPosition
    integer :: zone
    logical :: done
    
    type(AtomLinkedListNode), pointer :: currentAtom
    
    nrOfAtomsGenerated = 0
    do atomInUnitCell = 1, unitcell%nrOfAtoms
        done = .false.
        unitCellLocation = (/ 0, 1, 1 /)
        do while (.true.)
            call nextUnitCellLocation(space, unitCell%positions(:, atomInUnitCell), unitCell%size, unitCellLocation, &
                    atomPosition, done)
            if (done) exit
            if (all(atomPosition(1:space%spatDimension) < space%universeSize(1:space%spatDimension))) then
                zone = unitCell%zones(atomInUnitCell)
                if (filter(atomPosition,zone)) then
                    nrOfAtomsGenerated = nrOfAtomsGenerated + 1
                    allocate(currentAtom)
                    nullify(currentAtom%nextAtom)
                    currentAtom%positions(1:space%spatDimension) = atomPosition(1:space%spatDimension)
                    currentAtom%positions((space%spatDimension + 1):3) = 0
                    currentAtom%fromAtomInUnitCell = atomInUnitCell
                    currentAtom%unitcellLocation = unitCellLocation
                    call insertElement(generatedAtoms, currentAtom)
                    nullify(currentAtom)
                end if
            endif
        end do
    end do
  contains
    
    subroutine nextUnitCellLocation(space, startPoint, stepSize, unitCellLocation, pos, done)
    implicit none
        type(SpaceStruct), intent(in) :: space
        real(dblprec), dimension(3), intent(in) :: startPoint
        real(dblprec), dimension(3), intent(in) :: stepSize
        integer, dimension(3), intent(inout) :: unitCellLocation
        real(dblprec), dimension(3), intent(out) :: pos
        logical, intent(out) :: done
    
        integer :: i
        
        pos = 0
        done = .false.
        do i = 1, space%spatDimension
            unitCellLocation(i) = unitCellLocation(i) + 1
            pos(i) = (unitCellLocation(i) - 1) * stepSize(i) + startPoint(i)
            if (pos(i) >= space%universeSize(i) - 1d-6) then
                if (i >= space%spatDimension) then
                    done = .true.
                    exit
                else
                    pos(i) = startPoint(i)
                    unitCellLocation(i) = 1
                endif
            else
                exit
            end if
        enddo
        pos = (unitcellLocation - 1) * stepSize + startPoint
    end subroutine nextUnitCellLocation
end subroutine generatePositions

!< Initializes the moments from the unit cell.
!! @param[in] atomList The atoms that the moments will be initialized for
!! @param[in] atomcell The unitcell that the moments will be initialized from.
subroutine initializeMoments(atomList, unitcell)
implicit none
    type(AtomLinkedList), intent(in) :: atomList
    type(AtomCell), intent(in) :: unitcell
    
    call traverseAtomLinkedList(atomList, setMoment)
  contains
    subroutine setMoment(atom, i)
    implicit none
        type(AtomLinkedListNode), pointer, intent(in) :: atom
        integer, intent(in) :: i
        
        if (.false.) print *, i !Hacky solution to prevent warning about unused
        atom%moment = unitcell%moments(:, atom%fromAtomInUnitCell)
        atom%moment = atom%moment / sqrt(sum(atom%moment**2))
        atom%moment = atom%moment * unitcell%momentMagnitudes(atom%fromAtomInUnitCell)
    end subroutine setMoment
end subroutine initializeMoments

!< Initializes the atom types
!! @param[in] atomList The atoms that the types should be initialized for
!! @param[in] unitcell The unitcell that the initialization will be using
subroutine initializeAtomType(atomList, unitcell)
implicit none
    type(AtomLinkedList), intent(in) :: atomList
    type(AtomCell), intent(in) :: unitcell
    
    call traverseAtomLinkedList(atomList, setAtomType)
  contains
    subroutine setAtomType(atom, i)
    implicit none
        type(AtomLinkedListNode), pointer, intent(in) :: atom
        integer, intent(in) :: i
        
        if (.false.) print *, i !Hacky solution to prevent warning about unused
        atom%atomType = unitcell%chemicalTypes(atom%fromAtomInUnitCell)
    end subroutine
end subroutine initializeAtomType

!> Generates the atoms in the finite differences mesh. The indices for the atoms will also be generated, a negative index indicates that the
!! node should be interpolated.
!! @param[in] mesh The finite difference mesh
!! @param[in] nrOfAtomsBefore The number of atoms that is already inside generatedAtoms, this is used so the inidices for the fintie difference becomes correct.
!! @param[in] moment_magnitude The magnitude of the moments.
!! @param[out] atomIndices An array that will match the grid points in mesh, each entry says which atom in atoms that is stored there, 0 indicates no atom, a negative index indicates that it is a interpolation point.
!! @param[out] atoms The atoms that are generated.
!! @param[out] nrOfAtoms The number of atoms that were generated.
subroutine generateFiniteDifferenceAtoms(mesh, nrOfAtomsBefore, moment_magnitude, atomIndices, generatedAtoms, nrOfAtoms)
implicit none
    type(FiniteDiffMesh), intent(in) :: mesh
    integer, intent(in) :: nrOfAtomsBefore
    real(dblprec), intent(in) :: moment_magnitude
    integer, dimension(:, :, :), allocatable, intent(out) :: atomIndices
    type(AtomLinkedList), intent(inout) :: generatedAtoms
    integer, intent(out) :: nrOfAtoms
    
    integer :: i, j, k
    type(AtomLinkedListNode), pointer :: currentAtom
    
    allocate(atomIndices(mesh%nrOfGridPoints(1), mesh%nrOfGridPoints(2), mesh%nrOfGridPoints(3)))
    
    nrOfAtoms = 0
    do k = 1, mesh%nrOfGridPoints(3)
        do j = 1, mesh%nrOfGridPoints(2)
            do i = 1, mesh%nrOfGridPoints(1)
                if (shouldPlaceFiniteDiffAtom(mesh, (/i, j, k/))) then
                    allocate(currentAtom)
                    nullify(currentAtom%nextAtom)

                    nrOfAtoms = nrOfAtoms + 1
                    currentAtom%positions = ((/ i, j, k /) - 1)*mesh%boxSize
                    currentAtom%moment = (/ moment_magnitude, 0.0_dblprec, 0.0_dblprec /)
                    currentAtom%fromAtomInUnitCell = -1
                    currentAtom%unitcellLocation = -1
                    call insertElement(generatedAtoms, currentAtom)
                    nullify(currentAtom)
                    
                    if (isFiniteDiffGridPoint(mesh, (/i, j, k/))) then
                        atomIndices(i, j, k) = nrOfAtoms + nrOfAtomsBefore
                    else
                        atomIndices(i, j, k) = -(nrOfAtoms + nrOfAtomsBefore)
                    endif
                else
                    atomIndices(i, j, k) = 0
                endif
            enddo
        enddo
    enddo
    
  contains
    
    logical function shouldPlaceFiniteDiffAtom(mesh, ind)
    implicit none
        type(FiniteDiffMesh), intent(in) :: mesh
        integer, dimension(3) :: ind
        
        shouldPlaceFiniteDiffAtom = isFiniteDiffGridPoint(mesh, ind) .or. &
                                    isInterpolationGridPoint(mesh, ind)
    end function shouldPlaceFiniteDiffAtom
end subroutine generateFiniteDifferenceAtoms

!< Extract all the positions from atomList
!! @param[in] atomList
!! @param[out] positions The positions
subroutine extractPositions(atomList, positions)
implicit none
    type(AtomLinkedList), intent(in) :: atomList
    real(dblprec), dimension(:, :), intent(inout) :: positions
    
    call traverseAtomLinkedList(atomList, storePosition)
  contains
    subroutine storePosition(atom, i)
    implicit none
        type(AtomLinkedListNode), pointer, intent(in) :: atom
        integer, intent(in) :: i
        positions(:, i) = atom%positions
    end subroutine
end subroutine extractPositions

!< Extracts the moments from the atom list.
!! @param[in] atomList
!! @param[in] moments
subroutine extractMoments(atomList, moments)
implicit none
    type(AtomLinkedList), intent(in) :: atomList
    real(dblprec), dimension(:, :), intent(inout) :: moments
    
    call traverseAtomLinkedList(atomList, storeMoments)
  contains
    subroutine storeMoments(atom, i)
    implicit none
        type(AtomLinkedListNode), pointer, intent(in) :: atom
        integer, intent(in) :: i
        moments(:, i) = atom%moment
    end subroutine
end subroutine extractMoments

!< Extracts the types of the in the atom list
!! @param[in] atomList
!! @param[in,out] types
subroutine extractTypes(atomList, types)
implicit none
    type(AtomLinkedList), intent(in) :: atomList
    integer, dimension(:), intent(inout) :: types
    
    call traverseAtomLinkedList(atomList, storeTypes)
  contains
    subroutine storeTypes(atom, i)
    implicit none
        type(AtomLinkedListNode), pointer, intent(in) :: atom
        integer, intent(in) :: i
        
        types(i) = atom%atomType
    end subroutine
end subroutine extractTypes

subroutine extractFromUnitcellLocation(atomList, unitcellLocation)
implicit none
    type(AtomLinkedList), intent(in) :: atomList
    integer, dimension(:, :), intent(inout) :: unitcellLocation
    
    call traverseAtomLinkedList(atomList, storeUnitcellLocation)
  contains
    subroutine storeUnitcellLocation(atom, i)
    implicit none
        type(AtomLinkedListNode), pointer, intent(in) :: atom
        integer, intent(in) :: i
        
        unitcellLocation(:, i) = atom%unitcellLocation
    end subroutine
end subroutine extractFromUnitcellLocation

subroutine traverseAtomLinkedList(atomList, func)
implicit none
    type(AtomLinkedList), intent(in) :: atomList
    interface
        subroutine func(atom, atomCount)
        import :: AtomLinkedListNode
        implicit none
            type(AtomLinkedListNode), pointer, intent(in) :: atom
            integer, intent(in) :: atomCount
        end subroutine func
    end interface

    type(AtomLinkedListNode), pointer :: currentAtom
    integer :: i
    currentAtom => atomList%dummyElementBeforeFirst%nextAtom
    i = 1
    do while (associated(currentAtom))
        call func(currentAtom, i)
        i = i + 1
        currentAtom => currentAtom%nextAtom
    end do
end subroutine

subroutine allocateAtomLinkedList(atomList)
    type(AtomLinkedList), intent(inout) :: atomList
    
    allocate(atomList%dummyElementBeforeFirst)
    nullify(atomList%dummyElementBeforeFirst%nextAtom)
    atomList%lastElement => atomList%dummyElementBeforeFirst
end subroutine

subroutine deallocateAtomList(atomList)
implicit none
    type(AtomLinkedList), intent(inout) :: atomList
   
    type(AtomLinkedListNode), pointer :: nextAtom
    type(AtomLinkedListNode), pointer :: currentAtom
    
    currentAtom => atomList%dummyElementBeforeFirst
    do while (associated(currentAtom))
        nextAtom => currentAtom%nextAtom
        deallocate(currentAtom)
        currentAtom => nextAtom
    end do
    
    nullify(atomList%dummyElementBeforeFirst)
    nullify(atomList%lastElement)
end subroutine deallocateAtomList

!< Inserts an element into a list
!! @param[in,out] atomList
!! @param[in,out] atom
subroutine insertElement(atomList, atom)
implicit none
    type(AtomLinkedList), intent(inout) :: atomList
    type(AtomLinkedListNode), target, intent(inout) :: atom
    
    atomList%lastElement%nextAtom => atom
    atomList%lastElement => atomList%lastElement%nextAtom
    nullify(atom%nextAtom)
end subroutine

end module AtomGenerator
