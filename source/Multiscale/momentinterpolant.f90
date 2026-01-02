!>This module holds utility functions to solve one specific least-squares problem
!! that is used to find interpolation coefficients for each atom within a region.
!! The problem we solve is fitting a set of quadratic functions
!! f_j = a_1j·x^2+a_2j·y^2+a_3j·z^2 + a_4j·xy + a_5j·xz + a_6j·yz + a_7j·x + a_8j·y + a_9j·z + a_10j
!! Each f_j corresponds to one atom in the region.
!! We constrain f_j(pos.atom_j) = 1, f_j(pos.atom_i) = 0 if i/=j
!! And obtain the coefficients a_ij in a matrix.
!! We use weights when solving the problem.
!! The problem is written in matrix form as
!!  (Kt W K) M = Kt W I
!! solved as:
!!  M = (Kt W K)^-1 Kt W I
!! Where Kt is the transpose of K, which is the matrix of dependent variables of
!! the polynomials, that is x^2 y^2 z^2 xy...
!! W is a diagonal matrix containing the weights
!! M is the matrix of coeficients a_ij, solution of the problem
!! I is the Identity, which happens to be the matrix expressing the aforementioned constraints.
!! The inverse is solved by gaussian elimination (see below)
!! Many of the operations are performed inplace to save memory and reduce allocations.
!! The problem will be ill-conditioned when positions are too close to each other
!! and unsolvable if too few or repeated positions are given.
!! In such cases NaN´s will appear, beware!
!> authors
!> Edgar Mendez
!> Nikos  Ntallis
!> Manuel Pereiro

module MomentInterpolant
  use Parameters
implicit none
  private
  public  a_kt_diagw, invert, kt_diagw_k, Kcoeffs, ensureSize, rowMatrixProd, dot
contains

  
  !> Vector-vector scalar product
  function dot(row,col) result(r)
  implicit none
    real(dblprec), dimension(:), intent(in) :: row, col
    real(dblprec) :: r

    r = sum(row(:) * col(:))
  end function dot
  
  
  !> Multiplies one and only one row from A to the whole matrix B
  subroutine rowMatrixProd(A, B, nACols, nBCols, row, product)
    use DynamicArray
    implicit none
    real(dblprec), dimension(:,:), intent(in) :: A       
    real(dblprec), dimension(:,:), intent(in) :: B
    integer, intent(in) :: nACols
    integer, intent(in) :: nBCols
    integer, intent(in) :: row
    type(DynArrayReal), intent(inout) :: product

    integer :: i
    real(dblprec) :: x

    call clearArray(product)
    do i = 1,nBCols
       x = sum(A(row,1:nACols)*B(1:nACols,i))
       call addEntry(product,x)
    end do

  end subroutine rowMatrixProd
  
  !! Ensures that A shaped as n-by-m matrix.
  !! The values could be discarded, the size is never reduced, only grown.
  subroutine ensureSize(m, n, A)
    integer, intent(in) :: m
    integer, intent(in) :: n
    real(dblprec), dimension(:,:), allocatable, intent(inout) :: A

    if (.not. allocated(A)) then
       allocate(A(m,n))    
    elseif(ubound(A,1) < m .or. ubound(A,2) < n) then
       deallocate(A)
       allocate(A(m,n))
    endif

  end subroutine ensureSize

  !! Calculate K given positions
  subroutine Kcoeffs(space, positions,box,atomIndices, K)
    use GeometryModule
    use ShapeModule, only: BoxShape
    implicit none
    type(SpaceStruct), intent(in) :: space
    real(dblprec), dimension(:, :), intent(in) :: positions
    type(BoxShape), intent(in)   :: box
    integer, dimension(:), intent(in) :: atomIndices
    real(dblprec), dimension(:,:), intent(inout) :: K

    Integer :: i,atom
    real(dblprec), dimension(3) :: centre
    real(dblprec), dimension(3,ubound(atomIndices,1)) :: normpos

    centre = (box%corner + box%sizes/2)
    do i=1,ubound(atomIndices,1)
       atom = atomIndices(i)
       !! Care, gets the direction considering periodic boundaries!
       !! Notice: When the domain is too thin, an atom that is specified twice
       !! (e.g: due to PBC) will be assigned the same vector twice.
       !! This shouldn´t be happening right now only because areaCoefficients will
       !! list each atom once, but could generally cause tricky situations.
       call getDirectionalVector(space, &
            positions(:,atom), centre, normpos(:,i))
       normpos(:,i) = normpos(:,i) / (box%sizes/2)
    end do
    if(space%spatDimension == 1) then
       do i=1,ubound(normpos,2)
          K(i,1:3) = (/ normpos(1,i)**2, normpos(1,i), 1.0_dblprec /)
       enddo
    elseif(space%spatDimension == 2) then
       do i=1,ubound(normpos,2)
          K(i,1:6) = (/ normpos(1,i)**2, normpos(2,i)**2, normpos(1,i)*normpos(2,i), &
               normpos(1,i),normpos(2,i), 1.0_dblprec /)          
       enddo
    else
       do i=1,ubound(normpos,2)
          K(i,1:10) = &
               (/ normpos(1,i)**2, normpos(2,i)**2, normpos(3,i)**2, &
               normpos(1,i)*normpos(2,i), normpos(1,i)*normpos(3,i), normpos(2,i)*normpos(3,i), &
               normpos(1,i), normpos(2,i), normpos(3,i), 1.0_dblprec /)
       end do
    end if
  end subroutine Kcoeffs

  !! Calculates Kt W K, see module header for details.
  !! Out must be allocated beforehand
  pure subroutine kt_diagw_k(A,W, Out)
    real(dblprec), dimension(:,:), intent(in)  :: A
    real(dblprec), dimension(:), intent(in)    :: W
    real(dblprec), dimension(:,:), intent(out) :: Out

    integer :: i,j,k
    real(dblprec) :: accum
    !! E = (At diag(W) A), 
    !! thus E_ij = a_ik w_i a_kj
    !! It´s symmetric, we calculate each element once.
    do j=1,ubound(A,2)
       do i=j,ubound(A,2)
          accum = 0
          do k=1,ubound(W,1)
             accum = accum + (a(k,i) * w(k) * a(k,j))
          end do
          Out(i,j) = accum
          Out(j,i) = accum
       end do
    end do

  end subroutine kt_diagw_k

  !> Calculates A´s inverse in Inv using gaussian elimination
  !! Destroys A in the process
  subroutine invert(A,Inv)
    real(dblprec), dimension(:,:), intent(inout)  :: A
    real(dblprec), dimension(:,:), intent(out) :: Inv

    real(dblprec) :: lambda
    integer :: i,j

    do i=1,ubound(Inv,1)
       do j=1,ubound(Inv,2)
          Inv(i,j) = 0
          if(i.eq.j) Inv(i,j) = 1
       end do
    end do

    !! Zeros below diagonal:
    do i=1,ubound(Inv,1)-1
       do j=(i+1),ubound(Inv,1)
          lambda = A(j,i) / A(i,i)
          A(j,:) = A(j,:) - A(i,:)*lambda
          Inv(j,:) = Inv(j,:) - Inv(i,:)*lambda
       end do
    end do

    !! Zeros above
    do i=1,ubound(Inv,1)-1
       do j=(i+1),ubound(Inv,1)
          lambda = A(i,j) / A(j,j)
          A(i,:) = A(i,:) - A(j,:)*lambda
          Inv(i,:) = Inv(i,:) - Inv(j,:)*lambda
       end do
    end do

    !! Ones in diag (not really caring about A)
    do i=1,ubound(Inv,1)
       Inv(i,:) = Inv(i,:) / A(i,i)
       !A(i,:) = A(i,:) / A(i,i)
    end do
  end subroutine invert
  
  !! A Kt W
  subroutine a_kt_diagw(A,Ks,w, Out)
    implicit none
    real(dblprec), dimension(:,:), intent(inout)  :: A
    real(dblprec), dimension(:,:), intent(inout) :: Ks
    real(dblprec), dimension(:), intent(inout)  :: w
    real(dblprec), dimension(:,:), intent(out)  :: Out

    integer :: i,j,k
    real(dblprec) :: accum
    !! E = A K^t diag(W)
    !! thus E_ij = (Σ(k) A_ik Ks_jk ) W_j
    do j=1,ubound(W,1)
       do i=1,ubound(A,1)
          accum = 0
          do k=1,ubound(A,2)
             accum = accum + A(i,k)*Ks(j,k)
          end do
          Out(i,j) = accum * w(j)
       end do
    end do
  end subroutine a_kt_diagw

  
end module MomentInterpolant
