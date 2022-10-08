!> authors
!> Edgar Mendez
!> Nikos  Ntallis
!> Manuel Pereiro

module GradientIndices
  use parameters

  private
  public &
       generateDirectionalDerivativeLinks

contains

  !> Fills a matrix with weights to calculate gradient in the continuum.
  subroutine generateDirectionalDerivativeLinks(&
       realAtomTree, paddingAtomTree, positions, atomLatSp, mesh, finiteDiffIndices,&
       windowSize, vector, gradientCoeffs)
    use KdTreeModule
    use GeometryModule
    use FiniteDifference
    use DynamicArray
    use SparseMatrix
    use ShapeModule, only: BoxShape
    use SortModule
  implicit none
    type(KdTree), intent(in) :: realAtomTree
    type(KdTree), intent(in) :: paddingAtomTree
    real(dblprec), dimension(:, :), intent(in) :: positions
    real(dblprec), intent(in) :: atomLatSp    
    type(FiniteDiffMesh), intent(in) :: mesh
    integer, dimension(:, :, :), intent(in) :: finiteDiffIndices
    real(dblprec), dimension(3) :: windowSize, vector
    type(SpMatrix), intent(inout) :: gradientCoeffs

    integer :: i,j,k
    type(DynArrayInt) :: paddingNodes 


    call newArray(paddingNodes)
        
    !! Calculate links between atoms
    call AddAtomLinks(mesh,realAtomTree%indices,realAtomTree, paddingAtomTree, &
         windowSize, vector, positions, atomLatSp, gradientCoeffs)

    !! Links between nodes
    ! ToDo: Putting this triple loop here is opening a can of worms.
    !       This should be accessed somehow by a member of finitedifference
    do k = 1, ubound(finiteDiffIndices, 3)
       do j = 1, ubound(finiteDiffIndices, 2)
          do i = 1, ubound(finiteDiffIndices, 1)
             if (finiteDiffIndices(i, j, k) > 0) then ! regular node
                call addMeshLink((/i,j,k/), mesh, finiteDiffIndices, &
                     vector, gradientCoeffs)
             end if
          end do
       end do
    end do

        
    
    !! Sorting the matrix by row, the algorithm to convert an SpMatrix to
    !!  the UppASD interpolation structure is much simpler for sorted matrices.
    !! Sorting by column too isn't more expensive and could (hopefully)
    !!  bring a small cache/branch pred. hit rate improvement.
    call sortMatrixByRowAndColumn(gradientCoeffs)


    !if(0 .eq. 1) pause "Eraseme"
    !call allocSpMatrix(reference)
    !call MikhaAddAtomLinks(mesh,(/3378/),realAtomTree, paddingAtomTree, &
    !     windowSize, vector, positions, atomLatSp, reference)
    !call sortMatrixByRowAndColumn(gradientCoeffs)
!
!    print *,"Atom"
!    do i = 1,gradientCoeffs%row%length
!       if (gradientCoeffs%row%values(i).eq.3378) then
!          j = i
!          do while (gradientCoeffs%row%values(j) .eq. 3378)
!             do k=1,reference%col%length
!                if(gradientCoeffs%col%values(j) == reference%col%values(k)) then
!                   print *,"Ref ", gradientCoeffs%col%values(j), &
!                        " Mik ", gradientCoeffs%entries%values(j), &
!                        "  Me ", reference%entries%values(k), &
!                        " rt: ", gradientCoeffs%entries%values(j)/reference%entries%values(k)                
!                end if
!             end do
!             j=j+1
!          end do
!          exit
!       end if
!    end do
!
  end subroutine generateDirectionalDerivativeLinks


  
  subroutine addMeshLink(index, mesh, finiteDiffIndices, &
       vector, gradientCoeffs)
    use FiniteDifference
    use SparseMatrix
    implicit none
    integer, dimension(3) :: index
    type(FiniteDiffMesh), intent(in) :: mesh
    integer, dimension(:, :, :), intent(in) :: finiteDiffIndices
    real(dblprec),dimension(3),intent(in) :: vector
    type(SpMatrix), intent(inout) :: gradientCoeffs

    integer :: dimIterator, offset, dim
    logical :: inside
    integer, dimension(3) :: neighIndex, lookup
    real(dblprec) :: weight, selfWeight ! Narcissist coefficient
    integer :: neighbour, self

    weight = 0
    selfWeight = 0
    self = abs(finiteDiffIndices(index(1),index(2),index(3)))
    dim = mesh%space%spatDimension
    ! Add the bounding nodes of pre_pos to neighbours and directions.
    do dimIterator = 1, dim
       do offset = -1,1,2
          neighIndex = index
          neighIndex(dimIterator) = neighIndex(dimIterator) + offset
          call modularGrid(ubound(finiteDiffIndices), mesh%space%periodicBoundary,&
               neighIndex, lookup, inside)
          
          weight = vector(dimIterator) / (2.0*mesh%boxSize(dimIterator))
          weight = sign(weight, vector(dimIterator) * offset) 
          
          if(.not. inside) then
             selfWeight = selfWeight + weight
          else
             neighbour = abs(finiteDiffIndices(lookup(1),lookup(2),lookup(3)))
             
             if(abs(weight) > 1.0d-10) then
                call addMatrixEntry(gradientCoeffs,&
                     self,neighbour, weight)
             end if
          end if
       end do
    end do           
    if(abs(selfWeight) > 1.0d-10) then
       call addMatrixEntry(gradientCoeffs,&
            self,self, selfWeight)
    end if
    
  end subroutine addMeshLink

  subroutine addAtomLinks(&
       mesh, indices, &
       realAtomTree, paddingAtomTree, windowSize, &
       vector, positions, atomLatSp, gradientCoeffs)
    use AreaCoefficients, only: BoxesAreaCoefficients
    use KdTreeModule
    use GeometryModule
    use FiniteDifference
    use DynamicArray
    use SparseMatrix
    use ShapeModule, only: BoxShape
    use MomentInterpolant
    implicit none
    type(FiniteDiffMesh), intent(in) :: mesh
    integer, dimension(:), intent(in) :: indices
    type(KdTree), intent(in) :: realAtomTree
    type(KdTree), intent(in) :: paddingAtomTree
    real(dblprec), dimension(3), intent(in) :: windowSize
    real(dblprec), dimension(3), intent(in) :: vector
    real(dblprec), dimension(:,:), intent(in) :: positions
    real(dblprec), intent(in) :: atomLatSp
    type(SpMatrix), intent(inout) :: gradientCoeffs

    integer :: c, i, j, nRealAtoms, nPaddingAtoms, iatom
    real(dblprec), dimension(3) :: position
    integer, allocatable, dimension(:) :: allIndices

    type(BoxShape) :: box
    type(DynArrayInt)  :: localAtoms
    type(DynArrayReal) :: localAreas
    type(DynArrayReal) :: gradient_weights
    type(DynArrayReal) :: gradient_pre_weights

    real(dblprec), dimension(:,:), allocatable :: K
    real(dblprec), dimension(:,:), allocatable :: M       
    real(dblprec), dimension(:,:), allocatable :: A,Inv

    !! Number of columns of the matrix K
    integer,dimension(3),parameter :: kCoefficientsPerDim =(/3, 6, 10/)
    integer :: nKcoefficients, nNeigh

    real(dblprec) :: delta_norm, strength
    real(dblprec),dimension(3) :: gradient_delta

    type(KdTree) :: totalTree
    real(dblprec), dimension(3,2) :: eval_positions
    real(dblprec) :: a_scale,b_scale

    nRealAtoms = ubound(realAtomTree%indices,1)
    nPaddingAtoms = ubound(paddingAtomTree%indices,1)
    allocate(allIndices(nRealAtoms+nPaddingAtoms))
    allIndices(1:nRealAtoms) = realAtomTree%indices
    allIndices(nRealAtoms+1:nRealAtoms+nPaddingAtoms) = paddingAtomTree%indices
    call buildKdTree(totalTree, positions, allIndices)

    delta_norm =  atomLatSp * 0.01
    strength = sqrt(sum(vector**2))
    gradient_delta = (vector / strength) * delta_norm

    !! Allocate A and Inv
    nKcoefficients = kCoefficientsPerDim(mesh%space%spatDimension)
    allocate(A(nKcoefficients,nKcoefficients))
    allocate(Inv(nKcoefficients,nKcoefficients))
    
    ! Alloc locals (20 is an orientative size, may grow larger)
    call newArray(gradient_weights,20)
    call newArray(gradient_pre_weights,20)
    call newArray(localAtoms,20)
    call newArray(localAreas,20)

    ! Build K
    nKcoefficients  = kCoefficientsPerDim(mesh%space%spatDimension)

    ! Build box
    box%sizes = windowSize

    do i=1,ubound(indices,1)
       iatom = indices(i)
       position = positions(:,iatom)
       box%corner = positions(1:3,iatom) - (windowSize/2)

       call boxesAreaCoefficients(mesh,box, atomLatSp, positions, &
            totalTree, localAtoms, localAreas)
       
       ! Remove zero elemets
       c=localAreas%length
       j = 1
       do while (j <= c)
          if(abs(localAreas%values(j)) .lt. 1d-15) then
             localAreas%values(j) = localAreas%values(c)
             localAtoms%values(j) = localAtoms%values(c)
             c = c-1
          else
             j = j+1
          end if
       end do
       localAreas%length = c
       nNeigh = localAreas%length

       !localAreas%values(1:nNeigh) = -1*localAreas%values(1:nNeigh)
       
       !! Ensure K and M are N-by-nCoeff and nCoeff-by-N respectively,
       !!  where N is the number of atoms in scope
       call ensureSize(nNeigh, nKcoefficients, K)
       call ensureSize(nKcoefficients, nNeigh, M)

       !! Calculate K
       call Kcoeffs(mesh%space, positions,box,localAtoms%values(1:nNeigh), K)

       !! Compute A = K^t·diag(w)·K
       call kt_diagw_k(K,localAreas%values(1:nNeigh), A) 
       !! Compute Inv = A^-1
       call invert(A,Inv)
       !! Compute M = Inv·K^t·diag(w)
       call a_kt_diagw(Inv,K,localAreas%values(1:nNeigh), M)

       !! Compute the gradient coefficients vector
       eval_positions(:,1) = position + gradient_delta;
       eval_positions(:,2) = position - gradient_delta;

       if(iatom == 3005) then
          print *,"iatom",iatom
          print *,"neighs", localAtoms%values(1:nNeigh)          
          print *,"Kcoefs",K
          print *,"M",M
       end if
       
       !! Calculate K for solution positions
       call Kcoeffs(mesh%space, eval_positions,box,(/1,2/), K)
       !! Compute the weights at both eval points
       ! Note: gfortran will warn here that k is possibly uninitialized.
       ! It is allocated in ensureSize up there,
       ! I don't know how to supress the warning here
       call rowMatrixProd(K, M, ubound(K,2), nNeigh, 1, gradient_weights)
       call rowMatrixProd(K, M, ubound(K,2), nNeigh, 2, gradient_pre_weights)
       
       !! Normalize both coefficient vectors
       !! So that M·v gives a moment of adequate magnitude
       a_scale = sum(gradient_weights%values(1:nNeigh))
       b_scale = sum(gradient_pre_weights%values(1:nNeigh))

       gradient_weights%values(1:nNeigh) = &
            gradient_weights%values(1:nNeigh) * b_scale
       gradient_pre_weights%values(1:nNeigh) = &
            gradient_pre_weights%values(1:nNeigh) * a_scale
       !! The gradient is (M·v2 - M·v1) / (2|d|)
       !! which is M · (v2-v1)/(2|d|)
       !! Uppasd provides M for each timestep, we calculate here the rest
       gradient_weights%values(1:nNeigh) = &
            ( &
            gradient_weights%values(1:nNeigh)       &
            - gradient_pre_weights%values(1:nNeigh) &
            )/ (delta_norm*2.0*a_scale*b_scale);

       !! Insert coefs into the matrix
       do j=1,nNeigh
          if(abs(gradient_weights%values(j)) > 1.0d-10) then
             call addMatrixEntry(gradientCoeffs,iatom,&
                  localAtoms%values(j), gradient_weights%values(j) * strength)
          end if
       end do

    end do
  end subroutine addAtomLinks
  
  ! Let @ denote the outer prod.
  ! Let V denote the differential operator (nabla).
  ! Let ~= mean 'approximately equal'.
  ! Let s_ij be the direction vector from atom i to atom j.
  ! Let d be the vector director of our derivative.
  ! Assume the lattice is simmetric. That is
  !  FORALL s_ij EXISTS atom k s.t: s_ij = s_ik AND k nearest neighbour of i
  ! Our calculation follows from:
  ! Sum_j (s_ij@sij) V m_i ~= (sum_j s_ij@m_j)
  ! Let B = (Sum_j s_ij@sij), then
  ! V m_i ~= B^-1 (sum_j s_ij@m_j)
  ! d·(V m_i) ~= d·B^-1 (sum_j s_ij@m_j)
  subroutine MikhaAddAtomLinks(&
       mesh, indices, &
       realAtomTree, paddingAtomTree, windowSize, &
       gradient_delta, positions, atomLatSp, gradientCoeffs)
    use AreaCoefficients, only: BoxesAreaCoefficients
    use KdTreeModule
    use GeometryModule
    use FiniteDifference
    use DynamicArray
    use SparseMatrix
    use ShapeModule, only: BoxShape
    use MomentInterpolant
  implicit none
    type(FiniteDiffMesh), intent(in) :: mesh
    integer, dimension(:), intent(in) :: indices
    type(KdTree), intent(in) :: realAtomTree
    type(KdTree), intent(in) :: paddingAtomTree
    real(dblprec), dimension(3), intent(in) :: windowSize
    real(dblprec), dimension(3), intent(in) :: gradient_delta
    real(dblprec), dimension(:,:), intent(in) :: positions
    real(dblprec), intent(in) :: atomLatSp
    type(SpMatrix), intent(inout) :: gradientCoeffs

    integer :: dim
    integer :: nRealAtoms, nPaddingAtoms, i,iatom,j,k,l
    integer, allocatable, dimension(:) :: allIndices
    type(DynArrayInt) :: neighbours
    real(dblprec), dimension(mesh%space%spatDimension) :: v,vBi 
    real(dblprec), dimension(3,3) :: B,Bi ! Need not to be 3x3, may be smaller
    real(dblprec),dimension(3) :: position,s
    real(dblprec) :: weight
    type(KdTree) :: totalTree

    dim = mesh%space%spatDimension
    
    v = gradient_delta(1:dim)    
    
    call newArray(neighbours,10)
    
    nRealAtoms = ubound(realAtomTree%indices,1)
    nPaddingAtoms = ubound(paddingAtomTree%indices,1)
    allocate(allIndices(nRealAtoms+nPaddingAtoms))
    allIndices(1:nRealAtoms) = realAtomTree%indices
    allIndices(nRealAtoms+1:nRealAtoms+nPaddingAtoms) = paddingAtomTree%indices
    call buildKdTree(totalTree, positions, allIndices)
        
    do i = 1,ubound(indices,1)
       iatom = indices(i)
       position = positions(:,iatom)
       
       call getNeighbours(mesh%space, positions, totalTree, position,&
            atomLatSp * 1.1_dblprec, neighbours)

       ! B_i = sum_j s_ij@s_ij 
       B = 0.0_dblprec
       do j = 1,neighbours%length          
          ! s is our s_ij
          call getDirectionalVector(mesh%space, &
               position, positions(:,neighbours%values(j)), s)
          do k = 1,dim
             do l = 1,dim
                B(k,l) = B(k,l) + (s(k)*s(l))
             end do
          end do
       end do
       ! Bi = B^-1
       call invert(B(1:dim,1:dim),Bi(1:dim,1:dim))
       ! vBi = v·B^-1
       do k=1,dim
          vBi(k) = dot(v,Bi(1:dim,k))
       end do
       
       do j = 1,neighbours%length          
          ! s is our s_ij
          call getDirectionalVector(mesh%space, &
               position, positions(:,neighbours%values(j)), s)          
          weight = dot(vBi,s(1:dim))
          if(abs(weight) > 1.0d-17) then
             call addMatrixEntry(gradientCoeffs,iatom,&
                  neighbours%values(j), weight)
          end if           
       end do
       
    end do
    
  end subroutine MikhaAddAtomLinks
  
end module GradientIndices
