module neighbor_sfc_mod
  use iso_fortran_env, only: wp => real64, i32 => int32
  implicit none
  private
  public :: reorder_atoms_morton

contains

  !**********************************************************************
  ! Function: morton3D
  ! Purpose : Compute 3D Morton (Z-order) code by bit-interleaving
  ! Inputs  : ix, iy, iz - integer coordinates (0..2^B-1)
  !           B         - bits per dimension
  ! Returns : interleaved 3B-bit code
  !**********************************************************************
  function morton3D(ix, iy, iz, B) result(code)
    integer(kind=i32), intent(in) :: ix, iy, iz, B
    integer(kind=i32)             :: code
    integer                        :: i
    code = 0_i32
    do i = 0, B-1
      code = ieor(code, ishft(ibits(ix, i, 1),      3*i    ))
      code = ieor(code, ishft(ibits(iy, i, 1),      3*i + 1))
      code = ieor(code, ishft(ibits(iz, i, 1),      3*i + 2))
    end do
  end function morton3D


  !**********************************************************************
  ! Subroutine: quicksort_codes
  ! Purpose   : Sort codes(:) in-place and carry along perm(:)
  !**********************************************************************
  recursive subroutine quicksort_codes(codes, perm, left, right)
    integer(kind=i32), intent(inout) :: codes(:), perm(:)
    integer, intent(in)              :: left, right
    integer(kind=i32)                :: pivot, temp32
    integer                           :: i, j

    i = left; j = right; pivot = codes((left + right) / 2)
    do while (i <= j)
      do while (codes(i) < pivot)
        i = i + 1
      end do
      do while (codes(j) > pivot)
        j = j - 1
      end do
      if (i <= j) then
        temp32 = codes(i); codes(i) = codes(j); codes(j) = temp32
        temp32 = perm(i);  perm(i)  = perm(j);  perm(j)  = temp32
        i = i + 1; j = j - 1
      end if
    end do
    if (left < j)  call quicksort_codes(codes, perm, left, j)
    if (i < right) call quicksort_codes(codes, perm, i, right)
  end subroutine quicksort_codes


  !**********************************************************************
  ! Subroutine: reorder_atoms_morton
  ! Purpose   : Reorder atoms by Morton code for better cache locality.
  ! Inputs    : coords(:,1:N) - atom positions
  !             N              - number of atoms
  !             B              - bits per dimension (e.g. 10)
  ! Outputs   : coords reordered in-place
  !**********************************************************************
  subroutine reorder_atoms_morton(coords, N, B)
    real(wp), intent(inout)        :: coords(3, N)
    integer,   intent(in)          :: N, B

    integer(kind=i32), allocatable :: codes(:), perm(:)
    real(wp), allocatable          :: coords_tmp(:,:)
    real(wp)                       :: xmin, xmax, ymin, ymax, zmin, zmax
    real(wp)                       :: x_norm, y_norm, z_norm
    integer(kind=i32)             :: ix, iy, iz, max_q
    integer                        :: i

    ! Compute domain bounds from coords
    xmin = minval(coords(1,1:N)); xmax = maxval(coords(1,1:N))
    ymin = minval(coords(2,1:N)); ymax = maxval(coords(2,1:N))
    zmin = minval(coords(3,1:N)); zmax = maxval(coords(3,1:N))

    ! Allocate temporary arrays
    allocate(codes(N), perm(N), coords_tmp(3, N))

    ! Initialize permutation locally
    do i = 1, N
      perm(i) = i
    end do

    ! Precompute quantization max integer
    max_q = ishft(1_i32, B) - 1_i32

    ! Compute Morton codes
    do i = 1, N
      x_norm = (coords(1,i) - xmin) / max(xmax - xmin, 1.0_wp)
      y_norm = (coords(2,i) - ymin) / max(ymax - ymin, 1.0_wp)
      z_norm = (coords(3,i) - zmin) / max(zmax - zmin, 1.0_wp)
      ix = int(min(max(x_norm,0.0_wp),1.0_wp) * real(max_q, wp))
      iy = int(min(max(y_norm,0.0_wp),1.0_wp) * real(max_q, wp))
      iz = int(min(max(z_norm,0.0_wp),1.0_wp) * real(max_q, wp))
      codes(i) = morton3D(ix, iy, iz, B)
    end do

    ! Sort codes and permute indices
    call quicksort_codes(codes, perm, 1, N)

    ! Reorder coordinates in-place
    do i = 1, N
      coords_tmp(:, i) = coords(:, perm(i))
    end do
    coords = coords_tmp

    ! Clean up temporaries
    deallocate(codes, perm, coords_tmp)
  end subroutine reorder_atoms_morton

end module neighbor_sfc_mod
