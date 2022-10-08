module FiniteDifference
    use Parameters
    use ShapeModule
    use GeometryModule
    use DynamicArray
implicit none
    !> The Finite difference mesh consist of a number of grid points
    !! the points are located from (0, 0, 0) to universeSize.
    !! The spacing between the points are boxSize.
    !! A cube is the cube that 8 mesh points are corners to (in 3 dimensions)
    !! if we are in 2 dimensions the cube is a square, and in 1 dimension it is an interval
    !> authors
    !> Edgar Mendez
    !> Nikos  Ntallis
    !> Manuel Pereiro

    type FiniteDiffMesh
        type(SpaceStruct) :: space
        logical*1, dimension(:, :, :), allocatable :: isAtomisticCube
        real(dblprec), dimension(:, :, :), allocatable :: boundaryDistance
        integer, dimension(3) :: nrOfGridPoints
        integer, dimension(3) :: nrOfBoxes
        real(dblprec), dimension(3) :: boxSize
    end type

public createFiniteDiffMesh, deallocateFiniteDiffMesh, isPointInsideAtomisticBox, &
       getBoundaryDistanceAt, isFiniteDiffGridPoint, isInterpolationGridPoint, &
       getCornersOfCubeAtPoint, createFiniteDiffLinks, FiniteDiffMesh, getNonInterpolationIndices, &
       modularGrid, getInterpolationIndices

private
contains

subroutine getCornersOfCubeAtPoint(mesh, finiteDiffIndices, point, corners)
implicit none
    type(FiniteDiffMesh), intent(in) :: mesh
    integer, dimension(:, :, :), intent(in) :: finiteDiffIndices
    real(dblprec), dimension(3), intent(in) :: point
    integer, dimension(2**mesh%space%spatDimension), intent(out) :: corners
    
    integer, dimension(3) :: boxIndex
    integer, dimension(3) :: neighIndex
    integer :: i, j, k
    integer :: cornerIndex
    logical :: inside
    integer :: kMax, jMax
    
    jMax = 0
    kMax = 0
    if (mesh%space%spatDimension >= 2) then
        jMax = 1
    end if
    if (mesh%space%spatDimension == 3) then
        kMax = 1
    end if
    
    boxIndex = floor(point / mesh%boxSize) + 1
    
    cornerIndex = 1
    do k = 0, kMax
        do j = 0, jMax
            do i = 0, 1
                call modularGrid(mesh%nrOfGridPoints, mesh%space%periodicBoundary, &
                             boxIndex + (/ i, j, k /), neighIndex, inside)
                if (.not. inside) then
                    print *, 'Something went wrong in getCornersOfCubeAtPoint'
                    print *, 'CubeIndex: ', boxIndex
                    print *, 'Point: ', point
                    print *, 'NrOfGridPoints: ', mesh%nrOfGridPoints
                    stop
                endif
                corners(cornerIndex) = finiteDiffIndices(neighIndex(1), neighIndex(2), neighIndex(3))
                cornerIndex = cornerIndex + 1
            enddo
        enddo
    enddo
end subroutine getCornersOfCubeAtPoint

!< Returns all the indices of the atoms that are not interpolation nodes
!! @param[in] meshIndices All the indices of the atoms
!! @param[in, out] indices The indices that are not interpolation indices
subroutine getNonInterpolationIndices(meshIndices, indices)
    use DynamicArray
implicit none
    integer, dimension(:, :, :), intent(in) :: meshIndices
    type(DynArrayInt), intent(inout) :: indices
    
    integer :: i, j, k
    
    do k = 1, ubound(meshIndices, 3)
        do j = 1, ubound(meshIndices, 2)
            do i = 1, ubound(meshIndices, 1)
                if (meshIndices(i, j, k) > 0) then
                    call addEntry(indices, meshIndices(i, j, k))
                end if
            end do
        end do
    end do
end subroutine

!< Returns all the indices of the atoms that are interpolation nodes
!! @param[in] meshIndices All the indices of the atoms
!! @param[in, out] indices The indices that are not interpolation indices
subroutine getInterpolationIndices(meshIndices, indices)
    use DynamicArray
implicit none
    integer, dimension(:, :, :), intent(in) :: meshIndices
    type(DynArrayInt), intent(inout) :: indices
    
    integer :: i, j, k
    
    do k = 1, ubound(meshIndices, 3)
        do j = 1, ubound(meshIndices, 2)
            do i = 1, ubound(meshIndices, 1)
                if (meshIndices(i, j, k) < 0) then
                    call addEntry(indices, -meshIndices(i, j, k))
                end if
            end do
        end do
    end do
end subroutine

!> Creates the finite difference mesh. The mesh will definie which cubes in it are atomistic and which are not.
!! Those mesh cubes whose corner is inside an atomistic shape will be determined to be atomistic.
!! @param[in] space
!! @param[in] nrOfBoxes How many boxes there should be in the mesh
!! @param[in] atomisticShapes The shapes that defines the atomistic region
!! @param[out] mesh The created mesh
subroutine createFiniteDiffMesh(space, nrOfBoxes, atomisticShapes, mesh)
implicit none
    type(SpaceStruct), intent(in) :: space
    integer, dimension(3), intent(in) :: nrOfBoxes
    type(ShapeList), intent(in) :: atomisticShapes
    type(FiniteDiffMesh), intent(out) :: mesh

    integer :: i, j, k
    real(dblprec), dimension(3) :: point
    real(dblprec), dimension(3) :: boxSize

    mesh%space = space
    do i = 1, 3
        if (i > space%spatDimension) then
            boxSize(i) = 0
            mesh%nrOfBoxes(i) = 1
            mesh%nrOfGridPoints(i) = 1
            mesh%space%universeSize(i) = 1
        else
            boxSize(i) = space%universeSize(i)/nrOfBoxes(i)
            mesh%nrOfBoxes(i) = nrOfBoxes(i)
            if (space%periodicBoundary(i)) then
                mesh%nrOfGridPoints(i) = nrOfBoxes(i)
            else
                mesh%nrOfGridPoints(i) = nrOfBoxes(i) + 1
            endif
        endif
    enddo
    
    allocate(mesh%isAtomisticCube(mesh%nrOfBoxes(1), mesh%nrOfBoxes(2), mesh%nrOfBoxes(3)))
    do i = 1, mesh%nrOfBoxes(3)
        point(3) = boxSize(3) * (i - 1)
        do j = 1, mesh%nrOfBoxes(2)
            point(2) = boxSize(2) * (j - 1)
            do k = 1, mesh%nrOfBoxes(1)
                point(1) = boxSize(1) * (k - 1)
                mesh%isAtomisticCube(k, j, i) = isBoxAtomistic(atomisticShapes, space%spatDimension, point, boxSize)
            enddo
        enddo
    enddo
    
    mesh%boxSize = mesh%space%universeSize / mesh%nrOfBoxes
    call generateBoundaryDistance(mesh)
end subroutine createFiniteDiffMesh

!< Links the finite difference nodes
!! @param[in] mesh
!! @param[in] continuous_exchange_coef
!! @param[in] continuum_dm
!! @param[in] meshIndices The indices of the atoms in the mesh
!! @param[in,out] exchange The links between the nodes will be added in this matrix
subroutine createFiniteDiffLinks(mesh, continuous_exchange_coef, continuum_dm, &
     meshIndices, exchange, dmMatrix)
    use SparseMatrix
implicit none
    type(FiniteDiffMesh), intent(in) :: mesh
    real(dblprec), intent(in), dimension(3)   :: continuous_exchange_coef
    real(dblprec), intent(in), dimension(3,3) :: continuum_dm 
    integer, dimension(:, :, :), intent(in) :: meshIndices
    type(SpMatrix), intent(inout) :: exchange
    type(SpMatrix), intent(inout) :: dmMatrix
    
    integer :: dimIterator
    integer, dimension(3) :: currentIndex
    integer, dimension(3) :: leftIndex, rightIndex
    integer, dimension(3) :: lookupCoordLeft
    integer, dimension(3) :: lookupCoordRight
    logical, dimension(3) :: wrapped_l,wrapped_r
    logical :: leftInside, rightInside
    integer :: i, j, k
    
    do k = 1, ubound(meshIndices, 3)
       do j = 1, ubound(meshIndices, 2)
          do i = 1, ubound(meshIndices, 1)
             if (meshIndices(i, j, k) <= 0) cycle
             currentIndex = (/ i, j, k /)
             
             do dimIterator = 1, mesh%space%spatDimension
                leftIndex = currentIndex
                rightIndex = currentIndex
                leftIndex(dimIterator) = currentIndex(dimIterator) - 1
                rightIndex(dimIterator) = currentIndex(dimIterator) + 1
                call modularGrid(ubound(meshIndices), mesh%space%periodicBoundary,&
                     leftIndex, lookupCoordLeft, leftInside, wrapped_l)
                call modularGrid(ubound(meshIndices), mesh%space%periodicBoundary,&
                     rightIndex, lookupCoordRight, rightInside, wrapped_r)

                
                if (.not. leftInside .and. rightInside) then
                   call addVacuumBoundary(exchange,dmMatrix, currentIndex, &
                        lookupCoordRight, mesh, meshIndices, &
                        dimIterator, continuous_exchange_coef, continuum_dm)
                elseif (.not. rightInside .and. leftInside) then
                   call addVacuumBoundary(exchange,dmMatrix, currentIndex, &
                        lookupCoordLeft, mesh, meshIndices, &
                        dimIterator, continuous_exchange_coef, continuum_dm)
                elseif (leftInside .and. rightInside) then                   
                   call linkToNeighbour(exchange,dmMatrix, currentIndex, &
                        lookupCoordLeft, &
                        mesh, meshIndices, dimIterator, &
                        continuous_exchange_coef, continuum_dm, &
                        wrapped_l(dimIterator))
                   call linkToNeighbour(exchange,dmMatrix, currentIndex, &
                        lookupCoordRight, &
                        mesh, meshIndices, dimIterator, &
                        continuous_exchange_coef, continuum_dm, &
                        wrapped_r(dimIterator))
                else
                   print *, 'Something is wrong in the linking ' &
                        // 'of the finite difference'
                   stop 1
                endif
                
             enddo
          enddo
       enddo
    enddo

 contains
   subroutine addVacuumBoundary(exchangeMatrix, dmMatrix, fromIndex, toIndex,&
        mesh, meshIndices, inDimension, continuous_exchange_coef, continuum_dm)
     use dm, only : addDmSp
    implicit none
        type(SpMatrix), intent(inout) :: exchangeMatrix
        type(SpMatrix), intent(inout) :: dmMatrix
        integer, dimension(3), intent(in) :: fromIndex
        integer, dimension(3), intent(in) :: toIndex
        type(FiniteDiffMesh), intent(in) :: mesh
        integer, dimension(:, :, :), intent(in) :: meshIndices
        integer, intent(in) :: inDimension
        real(dblprec), intent(in), dimension(3)   :: continuous_exchange_coef
        real(dblprec), intent(in), dimension(3,3) :: continuum_dm

        integer :: atomIndexFrom, atomIndexTo
        real(dblprec) :: exchangeCoef
        real(dblprec), dimension(3) :: link_continuum_dm

        atomIndexFrom = abs(meshIndices(fromIndex(1), fromIndex(2), fromIndex(3)))
        atomIndexTo = abs(meshIndices(toIndex(1), toIndex(2), toIndex(3)))
        exchangeCoef = 2 * continuous_exchange_coef(inDimension) / mesh%boxSize(inDimension)**2
        link_continuum_dm = 2*continuum_dm(:,inDimension) / mesh%boxSize(inDimension)
        
        if (fromIndex(inDimension) < toIndex(inDimension)) then
           link_continuum_dm = -link_continuum_dm
        end if
        
        call addMatrixEntry(exchangeMatrix, &
             atomIndexFrom, atomIndexTo,&
             exchangeCoef)        
        
        call addDmSp(dmMatrix, atomIndexFrom, atomIndexTo, &
             link_continuum_dm)
        
    end subroutine addVacuumBoundary

    subroutine linkToNeighbour(exchangeMatrix, dmMatrix, &
         fromIndex, neighIndex, mesh, meshIndices, &
         inDimension, continuous_exchange_coef, continuum_dm, &
         wrapped)
      use dm, only : addDmSp
    implicit none
        type(SpMatrix), intent(inout) :: exchangeMatrix
        type(SpMatrix), intent(inout) :: dmMatrix
        integer, dimension(3), intent(in) :: fromIndex
        integer, dimension(3), intent(in) :: neighIndex
        type(FiniteDiffMesh), intent(in) :: mesh
        integer, dimension(:, :, :), intent(in) :: meshIndices
        integer, intent(in) :: inDimension
        real(dblprec), intent(in), dimension(3)   :: continuous_exchange_coef 
        real(dblprec), intent(in), dimension(3,3) :: continuum_dm
        logical, intent(in) :: wrapped
        
        integer :: atomIndexFrom, atomNeighIndex
        real(dblprec) :: exchangeCoef
        real(dblprec), dimension(3) :: link_continuum_dm
        real(dblprec), dimension(3) :: link_continuum_dm_signed
        
        atomIndexFrom = abs(meshIndices(fromIndex(1), fromIndex(2), fromIndex(3)))
        atomNeighIndex = abs(meshIndices(neighIndex(1), neighIndex(2), neighIndex(3)))
        exchangeCoef = continuous_exchange_coef(inDimension) / mesh%boxSize(inDimension)**2
        if (wrapped) then
           link_continuum_dm = &
                -continuum_dm(:,inDimension) / mesh%boxSize(inDimension)
        else
           link_continuum_dm = &
                continuum_dm(:,inDimension) / mesh%boxSize(inDimension)
        end if
        
        call addMatrixEntry(exchangeMatrix, atomIndexFrom, atomNeighIndex, &
             exchangeCoef)
        
        if (fromIndex(inDimension) < neighIndex(inDimension)) then
           link_continuum_dm_signed = -link_continuum_dm
        else
           link_continuum_dm_signed = link_continuum_dm
        end if

        call addDmSp(dmMatrix, atomIndexFrom, atomneighIndex, &
             link_continuum_dm_signed)

      end subroutine linkToNeighbour
    
    
end subroutine createFiniteDiffLinks


!> Computes where an index is stored in a grid that might have
!! periodic boundary conditions. If the index is outside the grid
!! at a periodic boundary, the index is moved inside the grid where it belongs.
!! @param maxPoints The number of grid points in the mesh
!! @param periodicBoundary A .true. value indicates that it the boundary is periodic
!! @param coord The grid points coordinate
!! @param lookupCoord Where the coordinate is in the mesh
!! @param inside Tells if the point is inside the mesh or not, if it is not the lookupCoord is meaningless
!! @param wrapped True if the periodic boundary is was used to make the connection
pure subroutine modularGrid(maxPoints, periodicBoundary, coord, lookupCoord, inside, wrapped)
implicit none
    integer, dimension(3), intent(in) :: maxPoints
    logical, dimension(3), intent(in) :: periodicBoundary
    integer, dimension(3), intent(in) :: coord
    integer, dimension(3), intent(out) :: lookupCoord
    logical, intent(out) :: inside
    logical, dimension(3), intent(out), optional :: wrapped

    integer :: i
    logical, dimension(3) :: wrap
    
    inside = .true.
    wrap = .false.
    lookupCoord = coord    

    do i = 1, 3
       if (periodicBoundary(i)) then
          if (coord(i) <= 0 .or. coord(i) > maxPoints(i)) then
             lookupCoord(i) = modulo(coord(i) - 1, maxPoints(i)) + 1
             wrap(i) = .true.
          end if
        else
            if (coord(i) <= 0 .or. coord(i) > maxPoints(i)) then
                inside = .false.
                lookupCoord(i) = 0
            endif
        endif
     enddo
     if(present(wrapped)) then
        wrapped = wrap
     end if
end subroutine modularGrid

!< Generates the distance to the atomistic-continuous boundary
!! this is calculated for each box in the mesh.
!! @param[inout] mesh
subroutine generateBoundaryDistance(mesh)
implicit none
    type(FiniteDiffMesh), intent(inout) :: mesh
    
    real(dblprec) :: newDistance
    type LinkedListIndices
        integer, dimension(3) :: cubeInd
        integer, dimension(3) :: origin
        type(LinkedListIndices), pointer :: nextIndex
    end type
    type(LinkedListIndices), pointer :: startElement
    type(LinkedListIndices), pointer :: endElement
    nullify(startElement)
    nullify(endElement)

    allocate(mesh%boundaryDistance(mesh%nrOfBoxes(1), mesh%nrOfBoxes(2), mesh%nrOfBoxes(3)))

    mesh%boundaryDistance = huge(mesh%boundaryDistance(1, 1, 1))
    
    call setDistanceToZeroAtBoundary(mesh)
    call computeSmallestDistanceForAllBoxes();

    if (associated(endElement)) deallocate(endElement)
    
 contains
 
    subroutine setDistanceToZeroAtBoundary(mesh)
    implicit none
        type(FiniteDiffMesh), intent(inout) :: mesh
        integer :: i, j, k
        integer, dimension(3) :: boxIndex
        do k = 1, mesh%nrOfBoxes(3)
            do j = 1, mesh%nrOfBoxes(2)
                do i = 1, mesh%nrOfBoxes(1)
                    boxIndex = (/ i, j, k /)
                    if (isAtBoundary(mesh, boxIndex)) then
                        mesh%boundaryDistance(i, j, k) = 0
                        call insertElement(boxIndex, boxIndex)
                    endif
                enddo
            enddo
        enddo
    end subroutine
    
    subroutine computeSmallestDistanceForAllBoxes()
    implicit none
        integer :: i, j, k, ii
        integer, dimension(3) :: nghbrIndex
        logical :: inside
        type(LinkedListIndices), pointer :: currentElement

        do while (associated(startElement) .and. associated(startElement%nextIndex))
            currentElement => startElement
            startElement => startElement%nextIndex
            
            do ii = 0, 3**mesh%space%spatDimension - 1
                i = modulo(modulo(ii, 3) + 1, 3) - 1
                j = modulo(modulo(ii/3, 3) + 1, 3) - 1
                k = modulo(modulo(ii/9, 3) + 1, 3) - 1
                call modularGrid(mesh%nrOfBoxes, mesh%space%periodicBoundary, currentElement%cubeInd + (/ i, j, k /), &
                        nghbrIndex, inside)
                if (inside) then
                    newDistance = getDistance(mesh%space, (nghbrIndex - 1)*mesh%boxSize, &
                                                (currentElement%origin - 1) * mesh%boxSize)
                    if (mesh%boundaryDistance(nghbrIndex(1), nghbrIndex(2), nghbrIndex(3)) > newDistance) then
                        mesh%boundaryDistance(nghbrIndex(1), nghbrIndex(2), nghbrIndex(3)) = newDistance
                        call insertElement(nghbrIndex, currentElement%origin)
                    endif
                endif
            enddo
            deallocate(currentElement)
            nullify(currentElement)
        enddo
    end subroutine
    
    subroutine insertElement(boxIndex, originIndex)
    implicit none
        integer, dimension(3) :: boxIndex
        integer, dimension(3) :: originIndex
        
        if (.not. associated(startElement)) then
            allocate(startElement)
            nullify(startElement%nextIndex)
            endElement => startElement
        endif
        endElement%cubeInd = boxIndex
        endElement%origin = originIndex
        allocate(endElement%nextIndex)
        endElement => endElement%nextIndex
        nullify(endElement%nextIndex)
    end subroutine
    
    logical function isAtBoundary(mesh, cubeIndex)
    implicit none
        type(FiniteDiffMesh), intent(in) :: mesh
        integer, dimension(3) :: cubeIndex
        
        logical :: middleType
        
        
        integer :: ii
        integer :: i, j, k
        integer, dimension(3) :: nghbrIndex
        logical :: nghbrInsideMesh
        
        middleType = mesh%isAtomisticCube(cubeIndex(1), cubeIndex(2), cubeIndex(3))
        
        do ii = 0, 3**mesh%space%spatDimension - 1
            i = modulo(modulo(ii, 3) + 1, 3) - 1
            j = modulo(modulo(ii/3, 3) + 1, 3) - 1
            k = modulo(modulo(ii/9, 3) + 1, 3) - 1
            if (i == 0 .and. j == 0 .and. k == 0) cycle
            
            call modularGrid(mesh%nrOfBoxes, mesh%space%periodicBoundary, cubeIndex + (/ i, j, k /), &
                        nghbrIndex, nghbrInsideMesh)
            if (nghbrInsideMesh) then
                if (mesh%isAtomisticCube(nghbrIndex(1), nghbrIndex(2), nghbrIndex(3)) .neqv. middleType) then
                    isAtBoundary = .true.
                    return
                endif
            endif
        enddo
        isAtBoundary = .false.
    end function isAtBoundary
end subroutine generateBoundaryDistance

subroutine deallocateFiniteDiffMesh(mesh)
implicit none
    type(FiniteDiffMesh), intent(inout) :: mesh
    
    deallocate(mesh%isAtomisticCube)
    deallocate(mesh%boundaryDistance)
end subroutine

!> Returns the computed distance to the atomistic-continuous boundary 
!! @param mesh The mesh that it conserns.
!! @param point The point that the distance is taken from.
real(dblprec) function getBoundaryDistanceAt(mesh, point) result(boundaryDistance)
implicit none
    type(FiniteDiffMesh), intent(in) :: mesh
    real(dblprec), dimension(3), intent(in) :: point
    
    integer, dimension(3) :: meshIndex
    
    meshIndex = floor(point / (mesh%space%universeSize + 1d-8) * mesh%nrOfBoxes) + 1
    boundaryDistance = mesh%boundaryDistance(meshIndex(1), meshIndex(2), meshIndex(3))
end function getBoundaryDistanceAt

!> Returns true if the point is inside an atomistic box. Note that there is
!! a small extension of the continuous region at the atomistic-continuous boundary. I.e
!! there is a small region at the boundary where the point is inside an atomistic box, but
!! this function will say that it is not. This is to prevent atoms from beeing placed
!! exactly at the boundary of the continuous region.
!! @param[in] mesh
!! @param[out] point The point that is look if it is inside an atomistic box
logical function isPointInsideAtomisticBox(mesh, point) result(inside)
implicit none
    type(FiniteDiffMesh), intent(in) :: mesh
    real(dblprec), dimension(3), intent(in) :: point

    integer, dimension(3) :: meshIndex
    integer, dimension(3) :: neighIndex
    integer, dimension(3) :: lookupIndex
    logical :: insideNeigh
    integer :: i, j, k
    integer, dimension(3) :: maxNeigh
    real(dblprec), dimension(3) :: coordInsideBox
    real(dblprec) :: distToClosestAtomisticBox
    real(dblprec) :: currDist
    real(dblprec), dimension(3, -1:1) :: distances
    
    meshIndex = floor((point / mesh%space%universeSize) * mesh%nrOfBoxes) + 1
    !print *,meshIndex
    coordInsideBox = point - (meshIndex - 1) * mesh%boxSize
    
    distances(1, :) = (/ coordInsideBox(1), 0.0_dblprec, mesh%boxSize(1) - coordInsideBox(1) /)
    distances(2, :) = (/ coordInsideBox(2), 0.0_dblprec, mesh%boxSize(2) - coordInsideBox(2) /)
    distances(3, :) = (/ coordInsideBox(3), 0.0_dblprec, mesh%boxSize(3) - coordInsideBox(3) /)

    maxNeigh = merge((/ 1, 1, 1 /), (/ 0, 0, 0 /), (/ mesh%space%spatDimension >= 1, &
                                                      mesh%space%spatDimension >= 2, &
                                                      mesh%space%spatDimension == 3 /))

    distToClosestAtomisticBox = huge(distToClosestAtomisticBox)
    do k = -maxNeigh(3), maxNeigh(3)
        do j = -maxNeigh(2), maxNeigh(2)
            do i = -maxNeigh(1), maxNeigh(1)
                if (i == 0 .and. j == 0 .and. k == 0) cycle
                
                neighIndex = meshIndex + (/ i, j, k /)
                call modularGrid(mesh%nrOfBoxes, mesh%space%periodicBoundary, neighIndex, lookupIndex, insideNeigh)
                if (.not. insideNeigh) cycle
                if (mesh%isAtomisticCube(lookupIndex(1), lookupIndex(2), lookupIndex(3))) cycle
                currDist = sqrt(distances(1, i)**2 + distances(2, j)**2 + distances(3, k)**2)
                distToClosestAtomisticBox = min(distToClosestAtomisticBox, currDist)
            enddo
        end do
    end do
    inside = mesh%isAtomisticCube(meshIndex(1), meshIndex(2), meshIndex(3)) .and. &
            distToClosestAtomisticBox > 1d-6
end function isPointInsideAtomisticBox

!> Tells wheter a grid point is a finite difference grid point or not. It is a finite difference grid point
!! if any mesh cube surrounding it is a finite difference cube. This means that it is either on the boundary
!! between the atomistic region or completely inside the continuous region.
!! @param mesh
!! @param gridPoint The indices of the grid point that is checked if it is a finite difference grid point.
logical function isFiniteDiffGridPoint(mesh, gridPoint) result(nghbr)
implicit none
    type(FiniteDiffMesh), intent(in) :: mesh
    integer, dimension(3), intent(in) :: gridPoint
    
    integer :: i, j, k
    integer :: minI, minJ, minK
    integer :: maxI, maxJ, maxK
    integer, dimension(3) :: ind
    logical :: inside
    
    minI = max(1, gridPoint(1) - 1)
    minJ = max(1, gridPoint(2) - 1)
    minK = max(1, gridPoint(3) - 1)
    
    maxI = min(mesh%nrOfBoxes(1), gridPoint(1))
    maxJ = min(mesh%nrOfBoxes(2), gridPoint(2))
    maxK = min(mesh%nrOfBoxes(3), gridPoint(3))
    
    do k = gridPoint(3) - 1, gridPoint(3)
        do j = gridPoint(2) - 1, gridPoint(2)
            do i = gridPoint(1) - 1, gridPoint(1)
                call modularGrid(mesh%nrOfBoxes, mesh%space%periodicBoundary, (/i, j, k/), ind, inside)
                if (inside) then 
                    if (.not. mesh%isAtomisticCube(ind(1), ind(2), ind(3))) then
                        nghbr = .true.
                        return
                    endif
                endif
            enddo
        enddo
    enddo
    nghbr = .false.
end function isFiniteDiffGridPoint

logical function isInterpolationGridPoint(mesh, gridPoint) result(interpolation)
implicit none
    type(FiniteDiffMesh), intent(in) :: mesh
    integer, dimension(3), intent(in) :: gridPoint
    
    integer :: i, j, k
    integer, dimension(3) :: ind
    logical :: inside
    
    interpolation = .false.
    ! This iterates over a shape that looks like
    !  **
    ! ****
    ! ****
    !  **
    ! in two dimension. In three dimensions the shape looks similar,
    ! the edges of the cube is being removed.
    ! The corner closest to (0, 0, 0) in this shape is located at
    ! gridPoint - 2
    do k = gridPoint(3) - 2, gridPoint(3) + 1
        do j = gridPoint(2) - 2, gridPoint(2) + 1
            do i = gridPoint(1) - 2, gridPoint(1) + 1
                if (count((/ i == gridPoint(1) - 2, j == gridPoint(2) - 2, k == gridPoint(3) - 2/)) &
                 + count((/i == gridPoint(1) + 1, j == gridPoint(2) + 1, k == gridPoint(3) + 1/)) >= 2) then
                    cycle
                endif
                call modularGrid(mesh%nrOfBoxes, mesh%space%periodicBoundary, (/ i, j, k /), ind, inside)
                if (inside) then
                    if (.not. mesh%isAtomisticCube(ind(1), ind(2), ind(3))) then
                        interpolation = .true.
                        return
                    end if
                end if
            enddo
        enddo
    enddo
end function isInterpolationGridPoint

!< If a corner of a box is inside an atomistic shape, then the box is considered to
!! be atomistic.
!! @param[in] atomisticShapes The atomistic shapes.
!! @param[in] spatDimension The dimension of the space.
!! @param[in] boxCorner The corner of the box that is closest to (0, 0, 0).
!! @param[in] boxSize The size of the box.
logical function isBoxAtomistic(atomisticShapes, spatDimension, boxCorner, boxSize)
implicit none
    type(ShapeList), intent(in) :: atomisticShapes
    integer, intent(in) :: spatDimension
    real(dblprec), dimension(3) :: boxCorner
    real(dblprec), dimension(3) :: boxSize
    
    integer :: i
    real(dblprec), dimension(3, 2**spatDimension) :: cornerPoints
    do i = 0, 2**spatDimension -1
        cornerPoints(:, i + 1) = boxCorner + (/ boxSize(1) * ibits(i, 0, 1), &
                        boxSize(2) * ibits(i, 1, 1), boxSize(3) * ibits(i, 2, 1) /)
    end do
    
    isBoxAtomistic = .false.
    do i = 1, 2**spatDimension
        if (isPointInsideShape(atomisticShapes, cornerPoints(:, i))) then
            isBoxAtomistic = .true.
            return
        end if
    end do
end function isBoxAtomistic

end module FiniteDifference
