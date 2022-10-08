!> authors
!> Edgar Mendez
!> Nikos  Ntallis
!> Manuel Pereiro

module AtomLinker
    use Parameters
    use KdTreeModule
    use GeometryModule
    use SparseMatrix
    use DynamicArray
implicit none
    
    interface
        pure real(dblprec) function exchangeLaw(atomI, atomJ)
            use Parameters
        implicit none
            integer, intent(in) :: atomI
            integer, intent(in) :: atomJ
        end function
    end interface

    
public exchangeLaw, createPaddingInterpolationWeights, createExchangeMatrix, &
        createFiniteDiffInterpolationLinks
private
contains


subroutine createPaddingInterpolationWeights(atomPositions, paddingIndices, mesh, finiteDiffIndices, weights)
    use FiniteDifference
implicit none
    real(dblprec), dimension(:, :), intent(in) :: atomPositions
    integer, dimension(:), intent(in) :: paddingIndices
    type(FiniteDiffMesh), intent(in) :: mesh
    integer, dimension(:, :, :) :: finiteDiffIndices
    type(SpMatrix), intent(inout) :: weights
    
    integer, dimension(2**mesh%space%spatDimension) :: corners
    real(dblprec), dimension(2**mesh%space%spatDimension) :: tmpWeights
    real(dblprec), dimension(3) :: atomPosition,cornerPosition,direction
    integer :: i, j
    integer :: atomIndex
    
    do i = 1, ubound(paddingIndices, 1)
       atomIndex = paddingIndices(i)
       atomPosition = atomPositions(:,atomIndex)
       call getCornersOfCubeAtPoint(mesh, finiteDiffIndices, atomPosition, corners)
       do j = 1, 2**mesh%space%spatDimension
          if (corners(j) /= 0) then
             cornerPosition = atomPositions(:, abs(corners(j)))
             call getDirectionalVector(mesh%space, atomPosition, cornerPosition, &
                  direction)
             tmpWeights(j) = product(1 - abs(direction/mesh%boxSize(:)))
          end if
        enddo
        tmpWeights = tmpWeights / sum(tmpWeights)
        do j = 1, 2**mesh%space%spatDimension
            call addMatrixEntry(weights, atomIndex, abs(corners(j)), tmpWeights(j))
        enddo
    enddo
end subroutine createPaddingInterpolationWeights

!> Creates the weights for the finite difference interpolation nodes
!! The output will be a weight matrix w_ij where i is the index of the finite difference node
!! and j is the index of the atom in atomGeom
!! @param[in] space
!! @param[in] atomPositions Array holding positions for all atoms.
!! @param[in] atomIndices Tree indexing finite difference elements.
!! @param[in] boxSize The size of the cubes in the finite difference mesh
!! @param[in] atomLatSp Maximal distance between two atoms with a defined interaction.
!! @param[in] meshIndices The indices of the finite difference atoms ordered in a mesh
!! @param[in,out] weights The weight matrix
subroutine createFiniteDiffInterpolationLinks(mesh, space, atomPositions, atomIndices, &
     boxSize, atomLatSp, &
     meshIndices, weights)
    use AreaCoefficients
    use ShapeModule, only: BoxShape
    use SparseMatrix
    use KdTreeModule, only: KdTree
    use FiniteDifference
  implicit none
    type(FiniteDiffMesh), intent(in) :: mesh
    type(SpaceStruct), intent(in) :: space
    real(dblprec), dimension(:, :), intent(in) :: atomPositions
    type(KdTree), intent(in) :: atomIndices
    real(dblprec), dimension(3), intent(in) :: boxSize
    real(dblprec), intent(in) :: atomLatSp
    integer, dimension(:, :, :), intent(in) :: meshIndices
    type(SpMatrix), intent(inout) :: weights
    
    integer :: i, j, k
    type(DynArrayReal) :: interpolationWeights
    type(DynArrayInt) :: interpolationAtoms
    type(BoxShape) :: interpolationArea
    real(dblprec), dimension(3) :: expansion
    integer, parameter :: ITERATIVE_DEPTH = 6
    
    expansion = (/ 1d0, 1d0, 1d0 /) * ( atomLatSp * 1.01d0 )
    
    call newArray(interpolationAtoms)
    call newArray(interpolationWeights)
  
    do k = 1, ubound(meshIndices, 3)
       do j = 1, ubound(meshIndices, 2)
            do i = 1, ubound(meshIndices, 1)
                if (meshIndices(i, j, k) >= 0) cycle
                interpolationArea%corner = atomPositions(:, -meshIndices(i, j, k)) - boxSize*0.5_dblprec
                interpolationArea%sizes = boxSize
                call cutBoxToUniverse(space, interpolationArea)
                
                 
                !call iterativeAreaCoefficients(interpolationArea, expansion, atomPositions,&
                !     space, atomIndices, ITERATIVE_DEPTH, interpolationAtoms, interpolationWeights)

                call boxesAreaCoefficients(mesh, interpolationArea, atomLatSp, atomPositions, atomIndices, &
                     interpolationAtoms, interpolationWeights)
                
                call addWeightsToMatrix(-meshIndices(i, j, k), interpolationAtoms, interpolationWeights)
                
            enddo
        enddo
    enddo

    call deallocArray(interpolationAtoms)
    call deallocArray(interpolationWeights)
  contains
    subroutine cutBoxToUniverse(space, box)
    implicit none
        type(SpaceStruct), intent(in) :: space
        type(BoxShape), intent(inout) :: box
        
        integer :: i
        
        do i = 1, space%spatDimension
            if (space%periodicBoundary(i)) cycle

            if (box%corner(i) < 0) then
                box%sizes(i) = box%sizes(i) + box%corner(i)
                box%corner(i) = 0
            end if
            
            if (box%corner(i) + box%sizes(i) > space%universeSize(i)) then
                box%sizes(i) = space%universeSize(i) - box%corner(i)
            end if
        end do
    end subroutine 
    subroutine addWeightsToMatrix(currentAtom, interpolationAtoms, interpolationWeights)
    implicit none
        integer, intent(in) :: currentAtom
        type(DynArrayInt), intent(in) :: interpolationAtoms
        type(DynArrayReal), intent(in) :: interpolationWeights
        
        integer :: i
        do i = 1, interpolationAtoms%length
            call addMatrixEntry(weights, currentAtom, interpolationAtoms%values(i), &
                    interpolationWeights%values(i))
        end do
    end subroutine addWeightsToMatrix
end subroutine  createFiniteDiffInterpolationLinks


!> Creates the exchange matrix where entry ij is the exchange from atom i to j. This only links
!! the real, partially coarse grained and fully coarse grained atoms.
!! @param[in] space
!! @param[in] atomPositions The positions of the atoms
!! @param[in] totalTree A tree that contains all the atoms that should be linked to.
!! @param[in] realTree The real atoms
!! @param[in] partTree The partially coarse grained atoms
!! @param[in] coarseTree The coarse grained atoms
!! @param[in] realExchangeLaw The exchange law that descibes the exchange between the real atoms
!! @param[in] maxRealRadius The maximum exchange radius for the real atoms
!! @param[in] coarseExchangeLaw The exchange law that descibes the exchange between the coarse grained atoms
!! @param[in] maxCoarseRadius The maximum exchange radius for the coarse graiend exchange
!! @param[inout] exchange The resulting exchange matrix will be appended to this preallocated matrix.
subroutine createExchangeMatrix(&
     space, atomPositions, errorTolerance, &
     totalTree, realTree, partTree, coarseTree, &
     realExchangeLaw, maxRealRadius, coarseExchangeLaw, maxCoarseRadius, exchange)
implicit none
    type(SpaceStruct), intent(in) :: space
    real(dblprec), dimension(:, :), intent(in) :: atomPositions
    real(dblprec), intent(in) :: errorTolerance
    type(KdTree), intent(in) :: totalTree, realTree, partTree, coarseTree
    procedure(exchangeLaw) :: realExchangeLaw, coarseExchangeLaw
    real(dblprec), intent(in) :: maxRealRadius
    real(dblprec), intent(in) :: maxCoarseRadius
    type(SpMatrix), intent(inout) :: exchange

    integer :: i
    type(DynArrayInt) :: iIndex, jIndex
    type(DynArrayReal) :: realExchange, coarseExchange
    type(DynArrayReal) :: distancesSquared
    type(SpMatrix) :: knownExchanges
    type(SpMatrix) :: constraintMatrix
    real(dblprec), dimension(:), allocatable :: constraintVector
    real(dblprec), dimension(:), allocatable :: preferedExchange
    real(dblprec), dimension(:), allocatable :: fittedExchanges
    
    call linkAtomsToNeighbours(space, atomPositions, totalTree, realTree%indices, realExchangeLaw, maxRealRadius, exchange)
    call linkAtomsToNeighbours(space, atomPositions, totalTree, coarseTree%indices, coarseExchangeLaw, maxCoarseRadius, exchange)
    
    call newArray(iIndex)
    call newArray(jIndex)
    call newArray(realExchange)
    call newArray(coarseExchange)
    call newArray(distancesSquared)
    call allocSpMatrix(knownExchanges)

    ! Find J's from fully coarse-grained and real atoms to pcg
    call linkAtomsToNeighbours(space, atomPositions, coarseTree, partTree%indices, &
         coarseExchangeLaw, maxCoarseRadius, knownExchanges)
    call linkAtomsToNeighbours(space, atomPositions, realTree, partTree%indices, &
         realExchangeLaw, maxRealRadius, knownExchanges)
    ! Find links inside pcg
    call getUnknownExchanges(space, atomPositions, partTree, partTree%indices, &
         realExchangeLaw, maxRealRadius, coarseExchangeLaw, maxCoarseRadius, &
         errorTolerance, &
         iIndex, jIndex, realExchange, coarseExchange, distancesSquared)
        
    if (iIndex%length /= 0) then 
       call createConstraints(space, atomPositions, partTree%indices, knownExchanges,&
            iIndex, jIndex, constraintMatrix, constraintVector)
       call getPreferedExchange(space, atomPositions, realTree, coarseTree, &
            iIndex, jIndex, realExchange, coarseExchange, preferedExchange)
       call weightedMinimizationWithConstraints(preferedExchange,&
            distancesSquared%values(1:distancesSquared%length),&
            constraintMatrix, constraintVector, fittedExchanges)
    
       do i = 1, iIndex%length
          call addMatrixEntry(exchange, iIndex%values(i), jIndex%values(i), fittedExchanges(i))
          call addMatrixEntry(exchange, jIndex%values(i), iIndex%values(i), fittedExchanges(i))
       end do
       deallocate(preferedExchange)
       deallocate(constraintVector)
       deallocate(fittedExchanges)
       call deallocArray(distancesSquared)
       call deallocArray(iIndex)
       call deallocArray(jIndex)
       call deallocArray(realExchange)
       call deallocArray(coarseExchange)
       call deallocSpMatrix(constraintMatrix)
    end if
    do i = 1, knownExchanges%entries%length
        call addMatrixEntry(exchange, knownExchanges%row%values(i), knownExchanges%col%values(i), knownExchanges%entries%values(i))
    end do
    

    call deallocSpMatrix(knownExchanges)
end subroutine createExchangeMatrix

!> Links the atoms to its neighbours given a exchange law
!! @param[in] space
!! @param[in] atomPositions The positions of the atoms
!! @param[in] neighboursIndices The indices of the potential neighbours that should be linked to
!! @param[in] atomIndices The indices of the atoms that should be linked from.
!! @param[in] law The exchange law that will be linked with
!! @param[in] maxExchangeRadius The maximum exchange radius, outside this radius the exchange is assumed to be zero.
!! @param[in] exchange The result will be added to this preallocated matrix.
subroutine linkAtomsToNeighbours(space, atomPositions, neighboursIndices, atomIndices, law, maxExchangeRadius, exchange)
implicit none
    type(SpaceStruct), intent(in) :: space
    real(dblprec), dimension(:, :), intent(in) :: atomPositions
    type(KdTree), intent(in) :: neighboursIndices
    integer, dimension(:), intent(in) :: atomIndices
    procedure(exchangeLaw) :: law
    real(dblprec), intent(in) :: maxExchangeRadius
    type(SpMatrix), intent(inout) :: exchange
    
    type(DynArrayInt) :: neighIndices
    real(dblprec) :: exchangeCoef
    integer :: atomIndex
    integer :: i, j
    
    call newArray(neighIndices)
    
    do i = 1, ubound(atomIndices, 1)
        atomIndex = atomIndices(i)
        
        call getNeighbours(space, atomPositions, neighboursIndices, atomPositions(:, atomIndex), &
                        maxExchangeRadius, neighIndices)
        do j = 1, neighIndices%length
           if (atomIndex /= neighIndices%values(j)) then
              exchangeCoef = law(atomIndex, neighIndices%values(j))
              if (abs(exchangeCoef) > 1d-12) then
                 call addMatrixEntry(exchange, atomIndex, neighIndices%values(j), exchangeCoef)
              end if
           end if
        end do
    end do
    
    call deallocArray(neighIndices)
end subroutine linkAtomsToNeighbours

subroutine getUnknownExchanges(space, atomPositions, tree, indices, &
     realExchangeLaw, maxRealExchangeRadius,&
     coarseExchangeLaw, maxCoarseExchangeRadius, &
     errorTolerance, &
     iIndex, jIndex, realExchanges, coarseExchanges, linkDistancesSquared)
implicit none
    type(SpaceStruct), intent(in) :: space
    real(dblprec), dimension(:, :), intent(in) :: atomPositions
    type(KdTree), intent(in) :: tree
    integer, dimension(:), intent(in) :: indices
    procedure(exchangeLaw) :: realExchangeLaw
    real(dblprec) :: maxRealExchangeRadius
    procedure(exchangeLaw) :: coarseExchangeLaw
    real(dblprec) :: maxCoarseExchangeRadius
    real(dblprec) :: errorTolerance
    type(DynArrayInt), intent(inout) :: iIndex
    type(DynArrayInt), intent(inout) :: jIndex
    type(DynArrayReal), intent(inout) :: realExchanges
    type(DynArrayReal), intent(inout) :: coarseExchanges
    type(DynArrayReal), intent(inout) :: linkDistancesSquared

    real(dblprec), dimension(3) :: dirVect    
    integer :: atomIndex,neighbourIndex
    integer :: i, j
    type(DynArrayInt) :: neighIndices
    real(dblprec) :: squaredDistBetweenAtoms
    real(dblprec) :: coarseGrainedExchange, realExchange

    call newArray(neighIndices)
    call clearArray(iIndex)
    call clearArray(jIndex)
    do i = 1, ubound(indices, 1)
       atomIndex = indices(i)
       call getNeighbours(space, atomPositions, tree, atomPositions(:, atomIndex),&
            1.1 * max(maxRealExchangeRadius, maxCoarseExchangeRadius), neighIndices)
       do j = 1, neighIndices%length
          if (atomIndex < neighIndices%values(j)) then
             neighbourIndex = neighIndices%values(j)
             
             coarseGrainedExchange =  coarseExchangeLaw(atomIndex, neighbourIndex)
             realExchange = realExchangeLaw(atomIndex, neighbourIndex) 
             if(realExchange > 1d-12 .or. coarseGrainedExchange > 1d-12) then
                call getDirectionalVector(space, &
                     atomPositions(:, atomIndex),&
                     atomPositions(:, neighbourIndex), dirVect)

                squaredDistBetweenAtoms = sum(dirVect(1:space%spatDimension)**2)
                if (squaredDistBetweenAtoms > (maxRealExchangeRadius+errorTolerance)**2) &
                     realExchange = 0
                if (squaredDistBetweenAtoms > (maxCoarseExchangeRadius+errorTolerance)**2) &
                     coarseGrainedExchange = 0

                call addEntry(realExchanges, realExchange)
                call addEntry(coarseExchanges, coarseGrainedExchange)
                
                call addEntry(iIndex, atomIndex)
                call addEntry(jIndex, neighbourIndex)
                call addEntry(linkDistancesSquared,squaredDistBetweenAtoms)
             end if
          end if
       end do
    end do
    call deallocArray(neighIndices)
end subroutine getUnknownExchanges

subroutine createConstraints(space, atomPositions, partiallyCoarseGrained, knownExchanges, iIndex, jIndex, M, b)
    use SortModule
implicit none
    type(SpaceStruct), intent(in) :: space
    real(dblprec), dimension(:, :), intent(in) :: atomPositions
    integer, dimension(:), intent(in) :: partiallyCoarseGrained
    type(SpMatrix), intent(in) :: knownExchanges
    type(DynArrayInt), intent(in) :: iIndex
    type(DynArrayInt), intent(in) :: jIndex
    type(SpMatrix), intent(out) :: M
    real(dblprec), dimension(:), allocatable, intent(out) :: b
    
    integer :: linkId, knownEntry
    integer :: iRow, jRow
    integer :: coordIndex
    integer, dimension(:), allocatable :: sortedPartiallyCoarseGrained
    real(dblprec), dimension(3) :: dirVect

    allocate(sortedPartiallyCoarseGrained(size(partiallyCoarseGrained)))
    sortedPartiallyCoarseGrained = partiallyCoarseGrained
    call sort(sortedPartiallyCoarseGrained)
    call allocSpMatrix(M)
    do linkId = 1, iIndex%length
        call getDirectionalVector(space, atomPositions(:, iIndex%values(linkId)), atomPositions(:, jIndex%values(linkId)), dirVect)
        iRow = searchSortedArray(iIndex%values(linkId), sortedPartiallyCoarseGrained)
        jRow = searchSortedArray(jIndex%values(linkId), sortedPartiallyCoarseGrained)
        do coordIndex = 1, space%spatDimension
            call addMatrixEntry(M, (iRow - 1)*space%spatDimension + 1 + coordIndex - 1, linkId, dirVect(coordIndex))
            call addMatrixEntry(M, (jRow - 1)*space%spatDimension + 1 + coordIndex - 1, linkId, -dirVect(coordIndex))
        end do
    end do
    
    allocate(b(space%spatDimension * ubound(partiallyCoarseGrained, 1)))
    b = 0.0_dblprec
    do knownEntry = 1, knownExchanges%entries%length
        call getDirectionalVector(space, atomPositions(:, knownExchanges%row%values(knownEntry)), &
                                atomPositions(:, knownExchanges%col%values(knownEntry)), dirVect)
        iRow = searchSortedArray(knownExchanges%row%values(knownEntry), sortedPartiallyCoarseGrained)
        iRow = (iRow - 1)*space%spatDimension + 1
        do coordIndex = 1, space%spatDimension
            b(iRow + coordIndex - 1) = b(iRow + coordIndex - 1) - knownExchanges%entries%values(knownEntry)*dirVect(coordIndex)
        end do
     end do
    deallocate(sortedPartiallyCoarseGrained)
end subroutine createConstraints

subroutine getPreferedExchange(space, atomPositions, realTree, coarseGrainedTree, &
     iIndex, jIndex, realExchanges, coarseExchanges, preferedExchange)
implicit none
    type(SpaceStruct), intent(in) :: space
    real(dblprec), dimension(:, :), intent(in) :: atomPositions
    type(KdTree), intent(in) :: realTree
    type(KdTree), intent(in) :: coarseGrainedTree
    type(DynArrayInt), intent(in) :: iIndex
    type(DynArrayInt), intent(in) :: jIndex
    type(DynArrayReal), intent(in) :: realExchanges
    type(DynArrayReal), intent(in) :: coarseExchanges
    real(dblprec), dimension(:), allocatable, intent(out) :: preferedExchange

    integer :: i
    real(dblprec), dimension(3) :: dirVect
    real(dblprec), dimension(3) :: middle
    real(dblprec) :: coarseGrainedExchange, realExchange
    real(dblprec) :: distToReal, distToCoarse
    
    allocate(preferedExchange(iIndex%length))
    do i = 1, iIndex%length
        call getDirectionalVector(space, atomPositions(:, iIndex%values(i)), &
            atomPositions(:, jIndex%values(i)), dirVect)
        middle = atomPositions(:, iIndex%values(i)) + dirVect * 0.5_dblprec

        coarseGrainedExchange = coarseExchanges%values(i)
        realExchange = realExchanges%values(i)
        
        call distToClosest(space, atomPositions, realTree, middle, distance=distToReal)
        call distToClosest(space, atomPositions, coarseGrainedTree, middle, distance=distToCoarse)
        
        preferedExchange(i) = (coarseGrainedExchange * distToReal + realExchange * distToCoarse) / (distToReal + distToCoarse)        
    end do
end subroutine getPreferedExchange

!> Minimizes the norm of ||W(x - p)||, where W is diagonal, subject to Mx = b
!! NOTE: The constraint matrix is being changed in this subroutine, so it will
!! be unusable after a call to this subroutine
!! @param preferedResult The p vector
!! @param weights The square root of the diagonal of the matrix W
!! @param constraintMatrix The M matrix
!! @param constraintVector The b vector
!! @param result The x vector that minimizes the norm
subroutine weightedMinimizationWithConstraints(preferedResult, weights, constraintMatrix, constraintVector, result)
    use lsqrModule
implicit none
    real(dblprec), dimension(:), intent(in) :: preferedResult
    real(dblprec), dimension(:), intent(in) :: weights
    type(SpMatrix), intent(inout) :: constraintMatrix
    real(dblprec), dimension(:), intent(in) :: constraintVector
    real(dblprec), dimension(:), allocatable, intent(out) :: result
    
    real(dblprec), dimension(:), allocatable :: tempVector
    
    !Variables that are used by lsqr
    integer :: istop, itn
    real(dblprec) :: Anorm, Acond, rnorm, Arnorm, xnorm
    real(dblprec), dimension(:), allocatable :: se
        
    allocate(tempVector(ubound(constraintVector, 1)))
    allocate(result(ubound(weights, 1)))
    tempVector = 0.0_dblprec
    call addMatrixVectorProduct(constraintMatrix, preferedResult, tempVector)
    tempVector = constraintVector - tempVector
    call multMatrixWithInvertedDiagonalMatrix(constraintMatrix, weights) ! results in that constraintMatrix is M_w = M*W^(-1)
    
    result = preferedResult
    call lsqr(ubound(constraintVector, 1), ubound(result, 1), addWeightedConstraintVectorProduct, &
            addTransposedWeightedConstraintVectorProduct, tempVector, 0.0_dblprec, .false., result, se, &
            1d-12, 1d-12, 10d12, 40*ubound(constraintVector, 1), -1, istop, itn, Anorm, Acond, &
            rnorm, Arnorm, xnorm)

    if(istop >= 4) then
       print *, "WARNING: The least-squares problem of partial coarse-grained atoms"
       print *, " is poorly conditioned, a solution within bound could not be found."
       stop ""
    end if
    
    result = result/weights + preferedResult
    deallocate(tempVector)
  contains
    
    !> Computes M*W^(-1) where W is view as the diagonal of a diagonal matrix.
    subroutine multMatrixWithInvertedDiagonalMatrix(matrix, diagonalMatrix)
    implicit none
        type(SpMatrix), intent(inout) :: matrix
        real(dblprec), dimension(:), intent(in) :: diagonalMatrix
        
        integer :: i
        do i = 1, matrix%entries%length
            matrix%entries%values(i) = matrix%entries%values(i) / diagonalMatrix(matrix%col%values(i))
        end do
    end subroutine multMatrixWithInvertedDiagonalMatrix
    
    !> Computes accum = accum + M_w * x
    subroutine addWeightedConstraintVectorProduct(rows, columns, x, accum)
    implicit none
        integer, intent(in) :: rows
        integer, intent(in) :: columns
        real(dblprec), dimension(columns), intent(in) :: x
        real(dblprec), dimension(rows), intent(inout) :: accum
        
        call addMatrixVectorProduct(constraintMatrix, x, accum)
    end subroutine
    
    !> Computes result = result + M_w^t * x
    subroutine addTransposedWeightedConstraintVectorProduct(rows, columns, result, x)
    implicit none
        integer, intent(in) :: rows
        integer, intent(in) :: columns
        real(dblprec), dimension(columns), intent(inout) :: result
        real(dblprec), dimension(rows), intent(in) :: x
        
        call addTransposedMatrixVectorProduct(constraintMatrix, x, result)
    end subroutine 
end subroutine weightedMinimizationWithConstraints

end module
