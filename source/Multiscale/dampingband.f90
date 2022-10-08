!> The damping band is a low-pass filter added in the atomistic layer
!! when multiscale is in use. This band deals with reflections and high-frequency
!! waves that arise from the interface.
!! The dampig band is effectively a region in the atomistic layer where the equation
!! that governs the spin needs to incorporate an additional term.
!> authors
!> Edgar Mendez
!> Nikos  Ntallis
!> Manuel Pereiro
module DampingBand
  use DynamicArray
  use Parameters
  use GeometryModule
  use ShapeModule
  use KdTreeModule
  
public dampingPositionalCoefficients, dampingAreaCoeff
private
contains

  
  !> Calculates the coefficient that controls the strength of the effect of
  !! the damping band depending on the distance of each atom to the boundary.
  !! @param[in] dampingIndices Integer array of indices of atoms inside the band.
  !! @param[in] nonDampingIndices Tree index of atoms outside the band.
  !! @param[in] space Space topology
  !! @param[in] positions Array of atom positions
  !! @param[in] paddingIndices Tree index of padding atom positions.
  !! @param[in] atomLatSp Maximal distance between two atoms with a defined interaction. 
  !! @param[in] dampingWidth Width of the damping band.
  !! @param[in,out] coeffs Calculated coefficients
  subroutine dampingPositionalCoefficients(dampingIndices, nonDampingIndices, &
       space, positions, paddingIndices, atomLatSp, &
       strength, dampingWidth, coeffs)
        use KdTreeModule
    implicit none
    integer, dimension(:),  intent(in)         :: dampingIndices
    type(KdTree), intent(in)                   :: nonDampingIndices
    type(SpaceStruct), intent(in)              :: space
    real(dblprec), dimension(:, :), intent(in) :: positions
    type(KdTree), intent(in)                   :: paddingIndices
    real(dblprec), intent(in)                  :: atomLatSp
    real(dblprec), intent(in)                  :: strength
    real(dblprec), intent(in)                  :: dampingWidth
    real(dblprec), dimension(:), allocatable, intent(inout) :: coeffs 

    real(dblprec) :: dist2inside, dist2outside
    integer :: i,atom, closest_in, nrDampingAtoms
    
    nrDampingAtoms = ubound(dampingIndices,1)
    
    if(.not. allocated(coeffs)) then
       allocate(coeffs(nrDampingAtoms))
    end if
    
    do i=1,nrDampingAtoms
       atom = dampingIndices(i)
       call distToClosest(space, positions, paddingIndices, positions(:, atom), distance=dist2outside)
       closest_in = 0
       call distToClosest(space, positions, nonDampingIndices, positions(:,atom), distance=dist2inside, closest=closest_in)
       if (closest_in .eq. 0) then
          dist2inside = dampingWidth - dist2outside
       end if
       coeffs(i) = strength * (dist2inside / &
            (dist2inside+dist2outside-atomLatSp) )**2
    end do
    
  end subroutine dampingPositionalCoefficients

  !> Calculates the coefficients used for local spin averaging in the damping band.
  !! To retrieve the list of neighbors of the Nth atom in the damping zone use
  !! atomNeighbors(indices(n):indices(n+1)-1)
  !! @param mesh Finite differences mesh 
  !! @param space Space topology
  !! @param positions Atom positions
  !! @param atomIndices K-d tree of all atoms
  !! @param dampingIndices Integer array of indices of atoms inside the band.
  !! @param atomLatSp Maximal distance between two atoms with a defined interaction.
  !! @param windowSize The sides of the window where areas are calculated.
  !! @param dampingAreaCoeffs (inout) A sparse matrix g zone the index of the first neighboring atom in atomNeighbors (and thus, the index of its coefficient). An additional index pointing to the end of the array is added for convenience.
  subroutine dampingAreaCoeff(mesh, space, positions, atomIndices, dampingIndices,&
       atomLatSp, windowSize, dampingAreaCoeffs)
    use AreaCoefficients
    use SparseMatrix
    use FiniteDifference
    use MomentInterpolant
  implicit none
    type(FiniteDiffMesh), intent(in) :: mesh
    type(SpaceStruct), intent(in) :: space
    real(dblprec), dimension(:, :), intent(in) :: positions
    type(KdTree), intent(in) :: atomIndices
    integer, dimension(:) :: dampingIndices
    real(dblprec), intent(in) :: atomLatSp
    real(dblprec), dimension(:), intent(in) :: windowSize
    
    type(SpMatrix), intent(inout) :: dampingAreaCoeffs

    
    type(BoxShape) :: box
    integer :: i,j,atom,c
    
    type(DynArrayInt)  :: localAtoms
    type(DynArrayReal) :: localAreas
    
    real(dblprec), dimension(:,:), allocatable :: K
    real(dblprec), dimension(:,:), allocatable :: M

    real(dblprec), dimension(:,:), allocatable :: A,Inv
   
    !! Number of columns of the matrix K
    integer,dimension(3),parameter :: kCoefficientsPerDim =(/3, 6, 10/)
    integer :: nKcoefficients, nDampingAtoms, nLocalAtoms

    call newArray(localAtoms)
    call newArray(localAreas)
    
    nKcoefficients  = kCoefficientsPerDim(space%spatDimension)
    box%sizes = windowSize
    nDampingAtoms = ubound(dampingIndices,1)
    
    !! Allocate A and Inv
    allocate(A(nKcoefficients,nKcoefficients))
    allocate(Inv(nKcoefficients,nKcoefficients))

    do i=1,nDampingAtoms
       atom = dampingIndices(i)
       box%corner = positions(1:3,atom) - (windowSize/2)
       
       call boxesAreaCoefficients(mesh,box, atomLatSp, positions, &
            atomIndices, localAtoms, localAreas)

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
       nLocalAtoms = localAreas%length

       !! Ensure K and M are N-by-nCoeff and nCoeff-by-N respectively,
       !!  where N is the number of atoms in scope
       call ensureSize(nLocalAtoms, nKcoefficients, K)
       call ensureSize(nKcoefficients, nLocalAtoms, M)
       
       !! Calculate K
       call Kcoeffs(space, positions,box,localAtoms%values(1:nLocalAtoms), K)

       !! Compute A = K^t·diag(w)·K
       call kt_diagw_k(K,localAreas%values(1:nLocalAtoms), A) 
       !! Compute Inv = A^-1
       call invert(A,Inv)
       !! Compute M = Inv·K^t·diag(w)
       call a_kt_diagw(Inv,K,localAreas%values(1:nLocalAtoms), M)
       
       do j=1,nLocalAtoms
          call addMatrixEntry(dampingAreaCoeffs,atom,&
                   localAtoms%values(j),M(ubound(M,1),j))
       end do
    end do
    
    call deallocArray(localAtoms)
    call deallocArray(localAreas)
    deallocate(A)
    deallocate(Inv)
    
  end subroutine dampingAreaCoeff

end module dampingband
