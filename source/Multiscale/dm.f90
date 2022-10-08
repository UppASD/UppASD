!> authors
!> Edgar Mendez
!> Nikos  Ntallis
!> Manuel Pereiro

module DM
  use Parameters
  use KdTreeModule
  use GeometryModule
  use SparseMatrix
  use DynamicArray
  implicit none
  
  interface
     pure function dmLaw(atomI, atomJ) result(dm)
       use Parameters
     implicit none
       integer, intent(in) :: atomI
       integer, intent(in) :: atomJ
       real(dblprec),dimension(3) :: dm
     end function dmLaw
  end interface

  private    
  public &
       dmLaw, createDmMatrix, addDmSp

contains


  !> Add a dm vector to a dynamic array.
  !! Uses three indices in a dynamic array to hold one dm vector.
  subroutine addDmDyn(dm, val)
    implicit none
    type(DynArrayReal), intent(inout) :: dm
    real(dblprec), dimension(3), intent(in) :: val

    integer :: k
    
    do k=1,ubound(val,1)
       call addEntry(dm,val(k))
    end do
    
  end subroutine addDmDyn
  
  !> Fetch a dm vector from a dynamic array.
  subroutine getDmDyn(dm, index, val)
    implicit none
    type(DynArrayReal), intent(in) :: dm
    integer, intent(in) :: index
    real(dblprec), dimension(3), intent(inout) :: val
    
    integer :: begin,end
    begin = (index-1)*3+1
    end = begin+ubound(val,1)-1

    if (end > dm%length) then
       stop "Assert: Read DD past end"
    end if
    
    val = dm%values(begin:end)
    
  end subroutine getDmDyn

  
  !> Add a DM vector to a sparse matrix.
  !! It uses three consecutive actual columns to store a dm vector
  !! the given row and column are however increasing unit-by-unit, instead
  !! of in steps of three.
  subroutine addDmSp(dm, row,col, val)
    implicit none
    type(SpMatrix), intent(inout) :: dm
    integer, intent(in) :: row, col
    real(dblprec), dimension(3), intent(in) :: val

    integer :: k
    
    do k=1,ubound(val,1)
       call addMatrixEntry(dm,row,(col-1)*3+k,val(k))
    end do
    
  end subroutine addDmSp
  
  !> Fetch a dm vector from a sparse matrix.
  subroutine getDmSp(dm, index, row,col,val)
    implicit none
    type(SpMatrix), intent(in) :: dm
    integer, intent(in) :: index
    integer, intent(out) :: row,col
    real(dblprec), dimension(3), intent(inout) :: val
    
    integer :: begin,end
    begin = (index-1)*ubound(val,1)+1    
    end= begin+ubound(val,1)-1    
     
    row = dm%row%values(begin)
    col =(dm%col%values(begin)-1)/3+1
    val = dm%entries%values(begin:end)
    
  end subroutine getDmSp

  
  
  !> Creates the Dzyaloshinsky-Moriya matrix where entry [3*i,j 3i+1,j 3i+2,j]
  !! is the dm interaction vector from atom i to j. This only links
  !! the real, partially coarse grained and fully coarse grained atoms.
  !! @param[in] space
  !! @param[in] atomPositions The positions of the atoms
  !! @param[in] totalTree A tree that contains all the atoms that should be linked to.
  !! @param[in] realTree The real atoms
  !! @param[in] partTree The partially coarse grained atoms
  !! @param[in] coarseTree The coarse grained appptoms
  !! @param[in] realDmLaw The exchange law that descibes the exchange between the real atoms
  !! @param[in] maxRealRadius The maximum exchange radius for the real atoms
  !! @param[in] coarseDmLaw The exchange law that descibes the exchange between the coarse grained atoms
  !! @param[in] maxCoarseRadius The maximum exchange radius for the coarse graiend exchange
  !! @param[inout] dm The resulting dm matrix will be appended to this preallocated matrix.
  subroutine createDmMatrix(&
       space, atomPositions, errorTolerance,&       
       totalTree, realTree, partTree, coarseTree, &
       realDmLaw, maxRealRadius, coarseDmLaw, maxCoarseRadius, dm)
    implicit none
    type(SpaceStruct), intent(in) :: space
    real(dblprec), dimension(:, :), intent(in) :: atomPositions
    real(dblprec), intent(in) :: errorTolerance
    type(KdTree), intent(in) :: totalTree, realTree, partTree, coarseTree
    procedure(dmLaw) :: realDmLaw, coarseDmLaw
    real(dblprec), intent(in) :: maxRealRadius
    real(dblprec), intent(in) :: maxCoarseRadius
    type(SpMatrix), intent(inout) :: dm

    integer :: i, begin,end ,r,c
    real(dblprec),dimension(3) :: v
    type(DynArrayInt) :: iIndex, jIndex
    type(DynArrayReal) :: realDm, coarseDm
    type(DynArrayReal) :: distancesSquared
    type(DynArrayReal) :: De
    type(SpMatrix) :: knownDm
    type(SpMatrix) :: constraintMatrix
    real(dblprec), dimension(:), allocatable :: constraintVector
    real(dblprec), dimension(:), allocatable :: preferredDm
    real(dblprec), dimension(:), allocatable :: fittedDms

    call linkAtomsToNeighbours(space, atomPositions, totalTree, &
         realTree%indices, realDmLaw, maxRealRadius, dm)
    call linkAtomsToNeighbours(space, atomPositions, totalTree, &
         coarseTree%indices, coarseDmLaw, maxCoarseRadius, dm)

    call newArray(iIndex)
    call newArray(jIndex)
    call newArray(realDm)
    call newArray(coarseDm)
    call newArray(distancesSquared)
    call newArray(De)
    call allocSpMatrix(knownDm)

    call linkAtomsToNeighbours(space, atomPositions, coarseTree, partTree%indices,&
         coarseDmLaw, maxCoarseRadius, knownDm)
    call linkAtomsToNeighbours(space, atomPositions, realTree, partTree%indices,&
         realDmLaw, maxRealRadius, knownDm)
    
    call getUnknownDm(space, atomPositions, partTree, partTree%indices, &
         realDmLaw, maxRealRadius, coarseDmLaw, maxCoarseRadius, &
         errorTolerance, De, iIndex, jIndex, realDm, coarseDm, distancesSquared)
    
    if (iIndex%length /= 0) then
       
       call createConstraints(space, atomPositions, partTree%indices, knownDm,&
            iIndex, jIndex, De, constraintMatrix, constraintVector)
       call getPreferredDm(space, atomPositions, realTree, coarseTree, &
            iIndex, jIndex, realDm, coarseDm, preferredDm)
       
       print *, "Sizes: ", &
            " links: ", iIndex%length , &
            " cm:",maxval(constraintMatrix%row%values(1:constraintMatrix%row%length)),&
            "x",   maxval(constraintMatrix%col%values(1:constraintMatrix%col%length)),&
            " cv:",ubound(constraintVector,1), &
            " De:",De%length, &
            " PD:",preferredDm, &
            " cv:",constraintVector

       do i = 1, constraintMatrix%row%length / 3
          call getDmSp(constraintMatrix, i, r,c,v)
          print *,"CM : ",r,c,v
       end do
       
       call weightedMinimizationWithConstraints(preferredDm,&
            distancesSquared%values(1:distancesSquared%length),&
            constraintMatrix, constraintVector, fittedDms)
       
       begin = 1
       do i = 1, iIndex%length
          end   = begin + 3 - 1
          call addDmSp(dm, iIndex%values(i), &
               jIndex%values(i), fittedDms(begin:end))
          !call addDmSp(dm, jIndex%values(i), &
          !     iIndex%values(i), fittedDms(begin:end))
          begin = begin + 3
       end do
       
       deallocate(preferredDm)
       deallocate(constraintVector)
       deallocate(fittedDms)
       call deallocArray(distancesSquared)
       call deallocArray(iIndex)
       call deallocArray(jIndex)
       call deallocArray(realDm)
       call deallocArray(coarseDm)
       call deallocSpMatrix(constraintMatrix)
    end if
    
    do i = 1, knownDm%entries%length
       call addMatrixEntry(dm, knownDm%row%values(i), knownDm%col%values(i), knownDm%entries%values(i))
    end do

    call deallocSpMatrix(knownDm)
  end subroutine createDmMatrix

  
  !> Links the atoms to its neighbours given a exchange law
  !! @param[in] space
  !! @param[in] atomPositions The positions of the atoms
  !! @param[in] neighboursIndices The indices of the potential neighbours that should be linked to
  !! @param[in] atomIndices The indices of the atoms that should be linked from.
  !! @param[in] law The exchange law that will be linked with
  !! @param[in] maxDmRadius The maximum exchange radius, outside this radius the exchange is assumed to be zero.
  !! @param[in,out] dm The result will be added to this preallocated matrix.
  subroutine linkAtomsToNeighbours(space, atomPositions, neighboursIndices, atomIndices, law, maxDmRadius, dm)
    implicit none
    type(SpaceStruct), intent(in) :: space
    real(dblprec), dimension(:, :), intent(in) :: atomPositions
    type(KdTree), intent(in) :: neighboursIndices
    integer, dimension(:), intent(in) :: atomIndices
    procedure(dmLaw) :: law
    real(dblprec), intent(in) :: maxDmRadius
    type(SpMatrix), intent(inout) :: dm

    type(DynArrayInt) :: nearbyAtoms
    real(dblprec), dimension(3) :: dmCoef
    integer :: atomIndex
    integer :: i, j

    call newArray(nearbyAtoms)

    do i = 1, ubound(atomIndices, 1)
       atomIndex = atomIndices(i)
       
       call getNeighbours(space, atomPositions, neighboursIndices,&
            atomPositions(:, atomIndex), maxDmRadius, nearbyAtoms)
       
       do j = 1, nearbyAtoms%length
          if (atomIndex /= nearbyAtoms%values(j)) then
             dmCoef = law(atomIndex, nearbyAtoms%values(j))
             if (norm2(dmCoef) > 1d-12) then
                call addDmSp(dm, atomIndex, nearbyAtoms%values(j), dmCoef)
             end if
          end if
       end do
    end do

    call deallocArray(nearbyAtoms)
  end subroutine linkAtomsToNeighbours
  
  
  subroutine getUnknownDm(space, atomPositions, tree, indices, &
       realDmLaw, maxRealDmRadius, &
       coarseDmLaw, maxCoarseDmRadius, &
       errorTolerance, &
       De, iIndex, jIndex, realDms, coarseDms, linkDistancesSquared)
    implicit none
    type(SpaceStruct), intent(in) :: space
    real(dblprec), dimension(:, :), intent(in) :: atomPositions
    type(KdTree), intent(in) :: tree
    integer, dimension(:), intent(in) :: indices
    procedure(dmLaw) :: realDmLaw
    real(dblprec) :: maxRealDmRadius
    procedure(dmLaw) :: coarseDmLaw
    real(dblprec) :: maxCoarseDmRadius
    real(dblprec) :: errorTolerance
    type(DynArrayReal), intent(inout) :: De
    type(DynArrayInt),  intent(inout) :: iIndex
    type(DynArrayInt),  intent(inout) :: jIndex
    type(DynArrayReal), intent(inout) :: realDms
    type(DynArrayReal), intent(inout) :: coarseDms
    type(DynArrayReal), intent(inout) :: linkDistancesSquared
    
    real(dblprec), dimension(3) :: dirVect    
    integer :: atomIndex,neighbourIndex
    integer :: i, j
    type(DynArrayInt) :: neighIndices
    real(dblprec) :: squaredDistBetweenAtoms
    real(dblprec),dimension(3) :: coarseGrainedDm, realDm
    real(dblprec) :: DeAccum
    
    call newArray(neighIndices)
    call clearArray(iIndex)
    call clearArray(jIndex)
    do i = 1, ubound(indices, 1)
       atomIndex = indices(i)
       call getNeighbours(space, &
            atomPositions, tree, atomPositions(:, atomIndex), &
            1.1 * max(maxRealDmRadius, maxCoarseDmRadius), neighIndices)
       DeAccum = 0
       do j = 1, neighIndices%length
          if (atomIndex < neighIndices%values(j)) then
             neighbourIndex = neighIndices%values(j)

             coarseGrainedDm =  coarseDmLaw(atomIndex, neighbourIndex)
             realDm = realDmLaw(atomIndex, neighbourIndex) 
             if(norm2(realDm) > 1d-12 .or. norm2(coarseGrainedDm) > 1d-12) then
                call getDirectionalVector(space, &
                     atomPositions(:, atomIndex),&
                     atomPositions(:, neighbourIndex), dirVect)

                squaredDistBetweenAtoms = sum(dirVect**2)
                if(squaredDistBetweenAtoms > (maxRealDmRadius + errorTolerance)**2)&
                     realDm = 0
                if(squaredDistBetweenAtoms > (maxCoarseDmRadius + errorTolerance)**2)&
                     coarseGrainedDm = 0
                
                call addDmDyn(realDms, realDm)
                call addDmDyn(coarseDms, coarseGrainedDm)

                call addEntry(iIndex, atomIndex)
                call addEntry(jIndex, neighbourIndex)
                call addEntry(linkDistancesSquared,squaredDistBetweenAtoms)
                DeAccum = DeAccum + sum(dirVect*(realDm + coarseGrainedDm))
             end if
          end if
       end do
       
       call addEntry(De, DeAccum * 5.0d-1)
       
    end do
    call deallocArray(neighIndices)
  end subroutine getUnknownDm

  subroutine createConstraints(space, atomPositions, partiallyCoarseGrained, knownDms, iIndex, jIndex, De, M, b)
    use SortModule
    implicit none
    type(SpaceStruct), intent(in) :: space
    real(dblprec), dimension(:, :), intent(in) :: atomPositions
    integer, dimension(:), intent(in) :: partiallyCoarseGrained
    type(SpMatrix), intent(in) :: knownDms
    type(DynArrayInt), intent(in) :: iIndex
    type(DynArrayInt), intent(in) :: jIndex
    type(DynArrayReal), intent(inout) :: De
    type(SpMatrix), intent(out) :: M   
    real(dblprec), dimension(:), allocatable, intent(out) :: b

    integer :: linkId, knownEntry, de_id, last_i
    integer :: iRow, jRow, row,col
    integer :: coordIndex
    integer, dimension(:), allocatable :: sortedPartiallyCoarseGrained
    real(dblprec), dimension(3) :: dirVect, dmVect

    sortedPartiallyCoarseGrained = partiallyCoarseGrained
    call sort(sortedPartiallyCoarseGrained)
    allocate(b(3 * ubound(partiallyCoarseGrained, 1)))

    b = 0.0_dblprec
    
    call allocSpMatrix(M)
    de_id = 0
    last_i = 0
    do linkId = 1, iIndex%length

       call getDirectionalVector(space, &
            atomPositions(:, iIndex%values(linkId)), &
            atomPositions(:, jIndex%values(linkId)), &
            dirVect)

       iRow = searchSortedArray(&
            iIndex%values(linkId),&
            sortedPartiallyCoarseGrained)
       jRow = searchSortedArray(&
            jIndex%values(linkId),&
            sortedPartiallyCoarseGrained)

       
       call addDmSp(M, linkId, iRow, dirVect)
       call addDmSp(M, linkId, jRow, dirVect)  
       if (iIndex%values(linkId) /= last_i) then
          de_id = de_id + 1
          last_i = iIndex%values(linkId)
       end if
       b(linkId) = De%values(de_id)       
       
    end do

    do knownEntry = 1, knownDms%row%length/3
       call getDmSp(knownDms, knownEntry, row,col,dmVect)
       call getDirectionalVector(space, &
            atomPositions(:, row), &
            atomPositions(:, col), &
            dirVect)
       
       iRow = searchSortedArray(row, sortedPartiallyCoarseGrained)
       !if(iRow > 0) then ! Not sure about this check
          iRow = (iRow - 1)*3
          do coordIndex = 1, 3
             b(iRow + coordIndex) = b(iRow + coordIndex) - sum(dmVect*dirVect(coordIndex))
          end do
       !end if
    end do
    deallocate(sortedPartiallyCoarseGrained)
  end subroutine createConstraints

  subroutine getPreferredDm(space, atomPositions, realTree, coarseGrainedTree, &
       iIndex, jIndex, realDms, coarseDms, preferredDm)
    implicit none
    type(SpaceStruct), intent(in) :: space
    real(dblprec), dimension(:, :), intent(in) :: atomPositions
    type(KdTree), intent(in) :: realTree
    type(KdTree), intent(in) :: coarseGrainedTree
    type(DynArrayInt), intent(in) :: iIndex
    type(DynArrayInt), intent(in) :: jIndex
    type(DynArrayReal), intent(in) :: realDms
    type(DynArrayReal), intent(in) :: coarseDms
    real(dblprec), dimension(:), allocatable, intent(out) :: preferredDm

    integer :: i
    real(dblprec), dimension(3) :: dirVect
    real(dblprec), dimension(3) :: middle
    real(dblprec), dimension(3) :: coarseGrainedDm, realDm
    real(dblprec) :: distToReal, distToCoarse

    allocate(preferredDm(iIndex%length*3))
    do i = 1, iIndex%length
       call getDirectionalVector(space, atomPositions(:, iIndex%values(i)), &
            atomPositions(:, jIndex%values(i)), dirVect)
       middle = atomPositions(:, iIndex%values(i)) + dirVect * 0.5_dblprec


       call getDmDyn(coarseDms,i,coarseGrainedDm)
       call getDmDyn(realDms,i,realDm)       

       call distToClosest(space, atomPositions, realTree, middle,&
            distance=distToReal)
       call distToClosest(space, atomPositions, coarseGrainedTree, middle,&
            distance=distToCoarse)

       preferredDm((i-1)*3+1:i*3) = &
            (coarseGrainedDm * distToReal + realDm * distToCoarse) &
            / (distToReal + distToCoarse)        
    end do
  end subroutine getPreferredDm

  !> Minimizes the norm of ||W(x - p)||, where W is diagonal, subject to Mx = b
  !! NOTE: The constraint matrix is being changed in this subroutine, so it will
  !! be unusable after a call to this subroutine
  !! @param preferredResult The p vector
  !! @param weights The square root of the diagonal of the matrix W
  !! @param constraintMatrix The M matrix
  !! @param constraintVector The b vector
  !! @param result The x vector that minimizes the norm
  subroutine weightedMinimizationWithConstraints(preferredResult, weights, constraintMatrix, constraintVector, result)
    use lsqrModule
    implicit none
    real(dblprec), dimension(:), intent(in) :: preferredResult
    real(dblprec), dimension(:), intent(in) :: weights
    type(SpMatrix), intent(inout) :: constraintMatrix
    real(dblprec), dimension(:), intent(in) :: constraintVector
    real(dblprec), dimension(:), allocatable, intent(inout) :: result

    real(dblprec), dimension(:), allocatable :: tempVector

    !Variables that are used by lsqr
    integer :: istop, itn, i,j
    real(dblprec) :: Anorm, Acond, rnorm, Arnorm, xnorm
    real(dblprec), dimension(:), allocatable :: se
    
    allocate(tempVector(ubound(constraintVector, 1)))
    allocate(result(ubound(weights, 1)))
    tempVector = 0.0_dblprec
    call addMatrixVectorProduct(constraintMatrix, preferredResult, tempVector)
    tempVector = constraintVector - tempVector
    call multMatrixWithInvertedDiagonalMatrix(constraintMatrix, weights) ! results in that constraintMatrix is M_w = M*W^(-1)
    
    result = preferredResult
    call lsqr(ubound(constraintVector, 1), ubound(result, 1), addWeightedConstraintVectorProduct, &
         addTransposedWeightedConstraintVectorProduct, tempVector, 0.0_dblprec, .false., result, se, &
         1d-12, 1d-12, 10d12, 40*ubound(constraintVector, 1), -1, istop, itn, Anorm, Acond, &
         rnorm, Arnorm, xnorm)
    
    print *,"W: " ,weights
    do i = 1, ubound(weights,1)
       do j = 1, 3
          result((i-1)*3+j) = result((i-1)*3+j) / weights(i)
       end do
    end do
    result = result + preferredResult
    deallocate(tempVector)
    
  contains

    !> Computes M*W^(-1) where W is the diagonal of a diagonal matrix.
    subroutine multMatrixWithInvertedDiagonalMatrix(matrix, diagonalMatrix)
      implicit none
      type(SpMatrix), intent(inout) :: matrix
      real(dblprec), dimension(:), intent(in) :: diagonalMatrix

      integer :: i
      do i = 1, matrix%entries%length
         matrix%entries%values(i) = matrix%entries%values(i) / diagonalMatrix((matrix%col%values(i)-1)/3+1)
      end do
    end subroutine multMatrixWithInvertedDiagonalMatrix

    !> Computes result = result + M_w * x
    subroutine addWeightedConstraintVectorProduct(rows, columns, x, result)
      implicit none
      integer, intent(in) :: rows
      integer, intent(in) :: columns
      real(dblprec), dimension(columns), intent(in) :: x
      real(dblprec), dimension(rows), intent(inout) :: result

      call addMatrixVectorProduct(constraintMatrix, x, result)
    end subroutine addWeightedConstraintVectorProduct

    !> Computes result = result + M_w^t * x
    subroutine addTransposedWeightedConstraintVectorProduct(rows, columns, result, x)
      implicit none
      integer, intent(in) :: rows
      integer, intent(in) :: columns
      real(dblprec), dimension(columns), intent(inout) :: result
      real(dblprec), dimension(rows), intent(in) :: x

      call addTransposedMatrixVectorProduct(constraintMatrix, x, result)
    end subroutine addTransposedWeightedConstraintVectorProduct
  end subroutine weightedMinimizationWithConstraints


    
end module DM
