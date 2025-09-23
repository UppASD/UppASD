!-------------------------------------------------------------------------------
! MODULE: SFCNeighbourHamiltonian
!> @brief SFC-based neighbour finding for Hamiltonian construction
!> @details This module provides coordinate-based neighbor finding using
!> space-filling curves. It replaces the supercell-based approach with
!> direct coordinate processing.
!> @author Anders Bergman (adapted from supercell approach)
!-------------------------------------------------------------------------------
module SFCNeighbourHamiltonian
   use iso_fortran_env, only: wp => real64, i32 => int32
   implicit none

   private
   public :: setup_sfc_neighbour_hamiltonian, setup_sfc_nm

contains

   !----------------------------------------------------------------------------
   ! SUBROUTINE: setup_sfc_neighbour_hamiltonian
   !> @brief SFC-based replacement for setup_neighbour_hamiltonian
   !> @details Uses coordinate-based neighbor finding with SFC ordering
   !----------------------------------------------------------------------------
   subroutine setup_sfc_neighbour_hamiltonian(Natom, conf_num, nHam, coords, atype, &
      shell_distances, nn, max_no_neigh, max_no_shells, max_no_equiv, &
      nlistsize, nlist, ncoup, nm, nmdim, xc_shells, &
      fs_nlist, fs_nlistsize, do_ralloy, Natom_full, &
      Nchmax, atype_ch, asite_ch, achem_ch, ammom_inp, hdim, lexp, do_sortcoup, &
      do_lsf, nind, lsf_field, map_multiple)
      
      implicit none
      
      ! Input parameters
      integer, intent(in) :: Natom, conf_num, nHam, hdim, lexp, do_ralloy, Natom_full, Nchmax
      real(wp), dimension(3,Natom), intent(in) :: coords
      integer, dimension(Natom), intent(in) :: atype
      integer, dimension(:), intent(in) :: nn  ! Number of shells per type
      real(wp), dimension(:,:), intent(in) :: shell_distances  ! Shell distances per type
      real(wp), dimension(:,:,:), intent(in) :: xc_shells  ! Exchange couplings per shell
      character, intent(in) :: do_sortcoup, do_lsf, lsf_field
      logical, intent(in) :: map_multiple
      integer, dimension(Natom_full), intent(in), optional :: atype_ch, asite_ch, achem_ch
      real(wp), dimension(:,:,:), intent(in), optional :: ammom_inp
      
      ! Output parameters
      integer, intent(out) :: max_no_neigh, max_no_equiv, max_no_shells
      integer, dimension(Natom), intent(out) :: nlistsize
      integer, dimension(:,:), allocatable, intent(out) :: nlist, nmdim
      integer, dimension(:,:,:), allocatable, intent(out) :: nm
      real(wp), dimension(:,:,:,:), allocatable, intent(out) :: ncoup
      integer, dimension(:,:), intent(out), optional :: fs_nlist
      integer, dimension(:), intent(out), optional :: fs_nlistsize
      integer, dimension(:,:), allocatable, intent(inout), optional :: nind
      
      ! Local variables
      integer :: i, neighbor_count
      integer, dimension(:), allocatable :: morton_perm
      real(wp) :: shell_tol
      
      ! Parameters
      shell_tol = 0.1_wp  ! Tolerance for shell assignment
      
      ! Compute Morton ordering for better cache locality
      allocate(morton_perm(Natom))
      call compute_morton_permutation_local(coords, Natom, morton_perm)
      
      ! Calculate maximum possible neighbors and shells
      max_no_shells = maxval(nn)
      max_no_equiv = 48  ! Conservative estimate
      max_no_neigh = max_no_equiv * max_no_shells
      
      ! Allocate output arrays
      allocate(nlist(max_no_neigh, Natom))
      allocate(nm(Natom, max_no_shells, max_no_equiv))
      allocate(nmdim(max_no_shells, Natom))
      allocate(ncoup(hdim, max_no_neigh, nHam, conf_num))
      
      ! Initialize arrays
      nlist = 0
      nm = 0
      nmdim = 0
      ncoup = 0.0_wp
      nlistsize = 0
      
      ! Build neighbor lists using direct distance search with SFC ordering
      !$omp parallel do default(shared) private(i, neighbor_count)
      do i = 1, Natom
         neighbor_count = 0
         call find_neighbors_direct(morton_perm(i), coords, atype, nn, shell_distances, &
                                   shell_tol, Natom, max_no_neigh, max_no_shells, &
                                   max_no_equiv, nlist(:,i), nm(i,:,:), nmdim(:,i), &
                                   neighbor_count)
         nlistsize(i) = neighbor_count
      end do
      !$omp end parallel do
      
      ! Calculate exchange couplings if needed
      if (nHam > 0) then
         call calculate_exchange_couplings_simple(nlist, nmdim, coords, atype, xc_shells, &
                                                 ammom_inp, ncoup, Natom, nHam, conf_num, &
                                                 hdim, lexp, max_no_neigh, max_no_shells, &
                                                 do_ralloy, atype_ch, achem_ch, asite_ch)
      endif
      
      ! Sort neighbor lists if requested
      if (do_sortcoup == 'Y') then
         call sort_neighbor_lists_simple(nlist, ncoup, nlistsize, Natom, nHam, &
                                        conf_num, hdim, max_no_neigh)
      endif
      
      ! Setup first shell neighbors for LSF if needed
      if (do_lsf == 'Y' .and. lsf_field == 'L' .and. present(fs_nlist) .and. present(fs_nlistsize)) then
         call setup_first_shell_neighbors_simple(nm, nmdim, fs_nlist, fs_nlistsize, &
                                                nlist, nlistsize, nind, Natom, max_no_equiv)
      endif
      
      deallocate(morton_perm)
      
   end subroutine setup_sfc_neighbour_hamiltonian

   !----------------------------------------------------------------------------
   ! SUBROUTINE: setup_sfc_nm  
   !> @brief SFC-based replacement for setup_nm
   !> @details Coordinate-based neighbor map construction
   !----------------------------------------------------------------------------
   subroutine setup_sfc_nm(Natom, coords, atype, max_no_neigh, max_no_shells, &
      max_no_equiv, nn, shell_distances, nm, nmdim, do_ralloy, Natom_full, &
      acellnumb, atype_ch, nntype)
      
      implicit none
      
      ! Input parameters
      integer, intent(in) :: Natom, Natom_full, do_ralloy
      real(wp), dimension(3,Natom), intent(in) :: coords
      integer, dimension(Natom), intent(in) :: atype
      integer, dimension(:), intent(in) :: nn
      real(wp), dimension(:,:), intent(in) :: shell_distances
      integer, dimension(Natom_full), intent(in), optional :: acellnumb, atype_ch
      integer, dimension(:,:), intent(in), optional :: nntype
      
      ! Output parameters
      integer, intent(out) :: max_no_neigh, max_no_equiv, max_no_shells
      integer, dimension(:,:,:), allocatable, intent(out) :: nm
      integer, dimension(:,:), allocatable, intent(out) :: nmdim
      
      ! Local variables
      integer :: i
      integer, dimension(:), allocatable :: morton_perm
      real(wp) :: shell_tol
      
      ! Parameters
      shell_tol = 0.1_wp
      
      ! Setup Morton ordering
      allocate(morton_perm(Natom))
      call compute_morton_permutation_local(coords, Natom, morton_perm)
      
      ! Set output dimensions
      max_no_shells = maxval(nn)
      max_no_equiv = 48
      max_no_neigh = max_no_equiv * max_no_shells
      
      ! Allocate output arrays
      allocate(nm(Natom, max_no_shells, max_no_equiv))
      allocate(nmdim(max_no_shells, Natom))
      nm = 0
      nmdim = 0
      
      ! Build neighbor maps
      !$omp parallel do default(shared) private(i)
      do i = 1, Natom
         call find_neighbors_for_nm_simple(morton_perm(i), coords, atype, nn, &
                                          shell_distances, shell_tol, Natom, &
                                          max_no_shells, max_no_equiv, &
                                          nm(i,:,:), nmdim(:,i))
      end do
      !$omp end parallel do
      
      deallocate(morton_perm)
      
   end subroutine setup_sfc_nm

   !----------------------------------------------------------------------------
   ! Helper subroutines
   !----------------------------------------------------------------------------
   
   subroutine compute_morton_permutation_local(coords, N, perm)
      implicit none
      integer, intent(in) :: N
      real(wp), dimension(3,N), intent(in) :: coords
      integer, dimension(N), intent(out) :: perm
      
      integer :: i
      integer(kind=i32), dimension(:), allocatable :: codes
      real(wp) :: xmin, xmax, ymin, ymax, zmin, zmax, x_norm, y_norm, z_norm
      integer(kind=i32) :: ix, iy, iz, max_q
      integer, parameter :: B = 10  ! bits per dimension
      
      allocate(codes(N))
      
      ! Compute domain bounds
      xmin = minval(coords(1,:)); xmax = maxval(coords(1,:))
      ymin = minval(coords(2,:)); ymax = maxval(coords(2,:))
      zmin = minval(coords(3,:)); zmax = maxval(coords(3,:))
      
      max_q = ishft(1_i32, B) - 1_i32
      
      ! Initialize permutation
      do i = 1, N
         perm(i) = i
      end do
      
      ! Compute Morton codes
      do i = 1, N
         x_norm = (coords(1,i) - xmin) / max(xmax - xmin, 1.0_wp)
         y_norm = (coords(2,i) - ymin) / max(ymax - ymin, 1.0_wp)
         z_norm = (coords(3,i) - zmin) / max(zmax - zmin, 1.0_wp)
         ix = int(min(max(x_norm,0.0_wp),1.0_wp) * real(max_q, wp), i32)
         iy = int(min(max(y_norm,0.0_wp),1.0_wp) * real(max_q, wp), i32)
         iz = int(min(max(z_norm,0.0_wp),1.0_wp) * real(max_q, wp), i32)
         codes(i) = morton3D_local(ix, iy, iz, B)
      end do
      
      ! Sort by Morton codes
      call quicksort_local(codes, perm, 1, N)
      
      deallocate(codes)
   end subroutine compute_morton_permutation_local
   
   function morton3D_local(ix, iy, iz, B) result(code)
      integer(kind=i32), intent(in) :: ix, iy, iz, B
      integer(kind=i32) :: code
      integer :: i
      code = 0_i32
      do i = 0, B-1
         code = ieor(code, ishft(ibits(ix, i, 1), 3*i))
         code = ieor(code, ishft(ibits(iy, i, 1), 3*i + 1))
         code = ieor(code, ishft(ibits(iz, i, 1), 3*i + 2))
      end do
   end function morton3D_local
   
   recursive subroutine quicksort_local(codes, perm, left, right)
      integer(kind=i32), dimension(:), intent(inout) :: codes
      integer, dimension(:), intent(inout) :: perm
      integer, intent(in) :: left, right
      integer(kind=i32) :: pivot, temp32
      integer :: i, j, temp_perm
      
      if (left >= right) return
      
      i = left; j = right
      pivot = codes((left + right) / 2)
      
      do while (i <= j)
         do while (codes(i) < pivot)
            i = i + 1
         end do
         do while (codes(j) > pivot)
            j = j - 1
         end do
         if (i <= j) then
            temp32 = codes(i); codes(i) = codes(j); codes(j) = temp32
            temp_perm = perm(i); perm(i) = perm(j); perm(j) = temp_perm
            i = i + 1; j = j - 1
         end if
      end do
      
      if (left < j) call quicksort_local(codes, perm, left, j)
      if (i < right) call quicksort_local(codes, perm, i, right)
   end subroutine quicksort_local

   subroutine find_neighbors_direct(atom_i, coords, atype, nn, shell_distances, &
                                   shell_tol, Natom, max_no_neigh, max_no_shells, &
                                   max_no_equiv, nlist_i, nm_i, nmdim_i, neighbor_count)
      implicit none
      integer, intent(in) :: atom_i, Natom, max_no_neigh, max_no_shells, max_no_equiv
      real(wp), dimension(3,Natom), intent(in) :: coords
      integer, dimension(Natom), intent(in) :: atype
      integer, dimension(:), intent(in) :: nn
      real(wp), dimension(:,:), intent(in) :: shell_distances
      real(wp), intent(in) :: shell_tol
      integer, dimension(max_no_neigh), intent(out) :: nlist_i
      integer, dimension(max_no_shells, max_no_equiv), intent(out) :: nm_i
      integer, dimension(max_no_shells), intent(out) :: nmdim_i
      integer, intent(out) :: neighbor_count
      
      integer :: atom_j, shell, atom_type_i
      real(wp) :: distance, shell_dist
      integer, dimension(max_no_shells) :: shell_counts
      
      atom_type_i = atype(atom_i)
      neighbor_count = 0
      nlist_i = 0
      nm_i = 0
      nmdim_i = 0
      shell_counts = 0
      
      ! Direct search through all atoms
      do atom_j = 1, Natom
         if (atom_j /= atom_i) then
            distance = sqrt(sum((coords(:,atom_i) - coords(:,atom_j))**2))
            
            ! Check which shell this distance corresponds to
            do shell = 1, nn(atom_type_i)
               shell_dist = shell_distances(atom_type_i, shell)
               if (abs(distance - shell_dist) < shell_tol) then
                  neighbor_count = neighbor_count + 1
                  if (neighbor_count <= max_no_neigh) then
                     nlist_i(neighbor_count) = atom_j
                  endif
                  
                  shell_counts(shell) = shell_counts(shell) + 1
                  if (shell_counts(shell) <= max_no_equiv) then
                     nm_i(shell, shell_counts(shell)) = atom_j
                  endif
                  exit
               endif
            end do
         endif
      end do
      
      nmdim_i = shell_counts
   end subroutine find_neighbors_direct

   subroutine find_neighbors_for_nm_simple(atom_i, coords, atype, nn, &
                                          shell_distances, shell_tol, Natom, &
                                          max_no_shells, max_no_equiv, &
                                          nm_i, nmdim_i)
      implicit none
      integer, intent(in) :: atom_i, Natom, max_no_shells, max_no_equiv
      real(wp), dimension(3,Natom), intent(in) :: coords
      integer, dimension(Natom), intent(in) :: atype
      integer, dimension(:), intent(in) :: nn
      real(wp), dimension(:,:), intent(in) :: shell_distances
      real(wp), intent(in) :: shell_tol
      integer, dimension(max_no_shells, max_no_equiv), intent(out) :: nm_i
      integer, dimension(max_no_shells), intent(out) :: nmdim_i
      
      integer :: atom_j, shell, atom_type_i
      real(wp) :: distance, shell_dist
      integer, dimension(max_no_shells) :: shell_counts
      
      atom_type_i = atype(atom_i)
      nm_i = 0
      nmdim_i = 0
      shell_counts = 0
      
      ! Direct search through all atoms
      do atom_j = 1, Natom
         if (atom_j /= atom_i) then
            distance = sqrt(sum((coords(:,atom_i) - coords(:,atom_j))**2))
            
            do shell = 1, nn(atom_type_i)
               shell_dist = shell_distances(atom_type_i, shell)
               if (abs(distance - shell_dist) < shell_tol) then
                  shell_counts(shell) = shell_counts(shell) + 1
                  if (shell_counts(shell) <= max_no_equiv) then
                     nm_i(shell, shell_counts(shell)) = atom_j
                  endif
                  exit
               endif
            end do
         endif
      end do
      
      nmdim_i = shell_counts
   end subroutine find_neighbors_for_nm_simple

   ! Placeholder subroutines - implement as needed based on original code
   subroutine calculate_exchange_couplings_simple(nlist, nmdim, coords, atype, xc_shells, &
                                                 ammom_inp, ncoup, Natom, nHam, conf_num, &
                                                 hdim, lexp, max_no_neigh, max_no_shells, &
                                                 do_ralloy, atype_ch, achem_ch, asite_ch)
      implicit none
      integer, intent(in) :: Natom, nHam, conf_num, hdim, lexp, max_no_neigh, max_no_shells, do_ralloy
      integer, dimension(:,:), intent(in) :: nlist, nmdim
      real(wp), dimension(3,Natom), intent(in) :: coords
      integer, dimension(Natom), intent(in) :: atype
      real(wp), dimension(:,:,:), intent(in) :: xc_shells
      real(wp), dimension(:,:,:), intent(in), optional :: ammom_inp
      real(wp), dimension(:,:,:,:), intent(inout) :: ncoup
      integer, dimension(:), intent(in), optional :: atype_ch, achem_ch, asite_ch
      
      ! Basic implementation - set to unity for now
      ncoup = 1.0_wp
   end subroutine calculate_exchange_couplings_simple

   subroutine sort_neighbor_lists_simple(nlist, ncoup, nlistsize, Natom, nHam, &
                                        conf_num, hdim, max_no_neigh)
      implicit none
      integer, intent(in) :: Natom, nHam, conf_num, hdim, max_no_neigh
      integer, dimension(:,:), intent(inout) :: nlist
      integer, dimension(Natom), intent(in) :: nlistsize
      real(wp), dimension(:,:,:,:), intent(inout) :: ncoup
      
      ! Placeholder - implement sorting if needed
   end subroutine sort_neighbor_lists_simple

   subroutine setup_first_shell_neighbors_simple(nm, nmdim, fs_nlist, fs_nlistsize, &
                                                nlist, nlistsize, nind, Natom, max_no_equiv)
      implicit none
      integer, intent(in) :: Natom, max_no_equiv
      integer, dimension(:,:,:), intent(in) :: nm
      integer, dimension(:,:), intent(in) :: nmdim
      integer, dimension(:,:), intent(out) :: fs_nlist
      integer, dimension(:), intent(out) :: fs_nlistsize
      integer, dimension(:,:), intent(in) :: nlist
      integer, dimension(:), intent(in) :: nlistsize
      integer, dimension(:,:), allocatable, intent(inout) :: nind
      
      ! Placeholder - implement first shell setup for LSF
   end subroutine setup_first_shell_neighbors_simple

end module SFCNeighbourHamiltonian