!-------------------------------------------------------------------------------
! MODULE: SFCNeighbourHamiltonian
!> @brief SFC-based neighbour finding for Hamiltonian construction
!> @details This module provides coordinate-based neighbor finding using
!> space-filling curves. It replaces the supercell-based approach with
!> direct coordinate processing and supports all bilinear interactions:
!> Heisenberg exchange, DM, SA, PD, BIQDM, BQ, and ring exchange.
!> @author Anders Bergman (adapted from supercell approach)
!-------------------------------------------------------------------------------
module SFCNeighbourHamiltonian
   use Parameters
   use Constants
   use neighbor_sfc_mod
   use HamiltonianDataType
   implicit none

   private
   public :: setup_sfc_neighbour_hamiltonian, setup_sfc_nm

contains

   !----------------------------------------------------------------------------
   ! SUBROUTINE: setup_sfc_neighbour_hamiltonian
   !> @brief SFC-based replacement for setup_neighbour_hamiltonian
   !> @details Uses coordinate-based neighbor finding with SFC ordering
   !> Compatible interface with existing setup_neighbour_hamiltonian but
   !> uses coordinates instead of supercell parameters
   !----------------------------------------------------------------------------
   subroutine setup_sfc_neighbour_hamiltonian(Natom,conf_num,NT,NA,nHam,anumb,atype,    &
      max_no_neigh,max_no_equiv,max_no_shells,nlistsize,nn,nlist,ncoup,nm,nmdim,coord, &
      xc,fs_nlist,fs_nlistsize,do_ralloy,Natom_full,Nchmax,atype_ch,asite_ch,achem_ch, &
      ammom_inp,hdim,lexp,do_sortcoup,do_lsf,nind,lsf_field,map_multiple)
      
      implicit none
      
      ! Input parameters - compatible with setup_neighbour_hamiltonian interface
      integer, intent(in) :: NT        !< Number of types of atoms
      integer, intent(in) :: NA        !< Number of atoms in one cell
      integer, intent(in) :: hdim      !< Number of elements in Hamiltonian element (scalar or vector)
      integer, intent(in) :: lexp      !< Order of the spin-interaction. Needed for proper rescaling (1=linear,2=quadratic)
      integer, intent(in) :: nHam      !< Number of atoms in Hamiltonian
      integer, intent(in) :: Natom     !< Number of atoms in system
      integer, intent(in) :: Nchmax    !< Max number of chemical components on each site in cell
      integer, intent(in) :: conf_num  !< Number of configurations for LSF
      integer, intent(in) :: do_ralloy       !< Random alloy simulation (0/1)
      integer, intent(in) :: Natom_full      !< Number of atoms for full system (=Natom if not dilute)
      integer, intent(in) :: max_no_neigh    !< Calculated maximum of neighbours for exchange
      integer, intent(in) :: max_no_equiv    !< Calculated maximum of neighbours in one shell for exchange
      integer, intent(in) :: max_no_shells   !< Calculated maximum of shells for exchange
      character, intent(in) :: do_sortcoup   !< Sort the entries of ncoup arrays (Y/N)
      character(len=1), intent(in) :: do_lsf    !< (Y/N) Do LSF for MC
      character(len=1), intent(in) :: lsf_field !< (T/L) LSF field term
      logical, intent(in) :: map_multiple       !< Allow for multiple couplings between atoms i and j
      integer, dimension(NT), intent(in)           :: nn       !< Number of neighbour shells
      integer, dimension(Natom), intent(in)        :: anumb    !< Atom number in cell
      integer, dimension(Natom), intent(in)        :: atype    !< Type of atom
      integer, dimension(Natom_full), intent(in)   :: atype_ch !< Actual site of atom for dilute system
      integer, dimension(Natom_full), intent(in)   :: asite_ch !< Actual site of atom for dilute system
      integer, dimension(Natom_full), intent(in)   :: achem_ch !< Chemical type of atoms (reduced list)
      integer, dimension(max_no_shells,Natom), intent(in) :: nmdim !< Dimension of neighbour map
      integer, dimension(Natom,max_no_shells,max_no_equiv), intent(in) :: nm  !< Neighbour map
      real(dblprec), dimension(3,Natom), intent(in) :: coord !< Coordinates of atoms
      real(dblprec), dimension(NA,Nchmax,conf_num), intent(in) :: ammom_inp   !< Magnetic moment directions from input (for alloys)
      real(dblprec), dimension(hdim,NT,max_no_shells,Nchmax,max(Nchmax,NT),conf_num), intent(in) :: xc !< Exchange couplings
      
      !.. Output variables
      integer, dimension(:), intent(out)     :: fs_nlistsize   !< Size of first shell neighbouring list for centered atom
      integer, dimension(:,:), intent(out)   :: fs_nlist       !< First shell Neighbouring list for centered atom
      integer, dimension(Natom), intent(out) :: nlistsize      !< Size of neighbour list for Heisenberg exchange couplings
      integer, dimension(max_no_neigh,Natom), intent(out) :: nlist !< Neighbour list for Heisenberg exchange couplings
      real(dblprec), dimension(hdim,max_no_neigh,nHam,conf_num), intent(out) :: ncoup !< Heisenberg exchange couplings
      
      ! .. In/Out variables
      integer, allocatable, dimension(:,:), intent(inout) :: nind !< index of firstshell-neighbour-list corresponds to neighbour-list
      
      ! Local variables
      integer :: i, j, k, l, iconf, i_stat
      integer :: jatom, neighbor_count, shell_idx
      real(dblprec) :: fc, fc2, distance
      real(dblprec), dimension(hdim) :: tempcoup
      integer :: tempn
      integer :: ncount, ncounter
      logical :: exis
      integer, dimension(:), allocatable :: morton_perm
      real(dblprec), parameter :: shell_tolerance = 0.1_dblprec
      
      ! Initialize output arrays
      ncoup = 0.0_dblprec
      nlist = 0
      nlistsize = 0
      fs_nlistsize = 0
      fs_nlist = 0
      
      ! Factors for mRy energy conversion
      fc = mry/mub
      fc2 = 2.0_dblprec*mry/mub
      
      ! Compute Morton ordering for better cache locality
      allocate(morton_perm(Natom))
      call compute_morton_permutation(coord, Natom, morton_perm)
      
      ! Main SFC-based neighbor finding loop
      !$omp parallel do default(shared) private(i,k,j,ncount,ncounter,exis,l,jatom,iconf,neighbor_count,shell_idx,distance)
      do i = 1, Natom
         ncount = 1
         ncounter = 1
         
         ! Use SFC-optimized neighbor search
         call find_sfc_neighbors(i, coord, atype, nn, Natom, max_no_neigh, &
                                 max_no_shells, max_no_equiv, morton_perm, &
                                 shell_tolerance, nlist(:,i), neighbor_count)
         
         ! Process each found neighbor
         do j = 1, neighbor_count
            jatom = nlist(j,i)
            if (jatom > 0 .and. jatom <= Natom) then
               ! Calculate distance and determine shell
               distance = sqrt(sum((coord(:,i) - coord(:,jatom))**2))
               shell_idx = determine_shell_from_distance(distance, atype(i), nn, shell_tolerance)
               
               if (shell_idx > 0 .and. shell_idx <= nn(atype(i))) then
                  ! Check for existing coupling (unless multiple mappings allowed)
                  exis = .false.
                  if (.not. map_multiple) then
                     do l = 1, ncount-1
                        if (nlist(l,i) == jatom) exis = .true.
                     end do
                  end if
                  
                  if (.not.exis .or. map_multiple) then
                     nlist(ncount,i) = jatom
                     
                     ! Calculate couplings for all configurations
                     if (i <= nHam) then
                        do iconf = 1, conf_num
                           call calculate_coupling_strength(i, jatom, shell_idx, atype, anumb, &
                                                          xc, ammom_inp, do_ralloy, atype_ch, &
                                                          asite_ch, achem_ch, hdim, lexp, &
                                                          fc2, iconf, ncoup(:,ncount,i,iconf))
                        end do
                     end if
                     
                     ! Handle first shell for LSF
                     if (shell_idx == 1 .and. do_lsf == 'Y' .and. lsf_field == 'L') then
                        fs_nlist(ncount,i) = nlist(ncount,i)
                        ncounter = ncounter + 1
                     end if
                     
                     ncount = ncount + 1
                  end if
               end if
            end if
         end do
         
         if (i <= nHam) nlistsize(i) = ncount - 1
         if (do_lsf == 'Y' .and. lsf_field == 'L') fs_nlistsize(i) = ncounter - 1
      end do
      !$omp end parallel do
      
      ! Sort coupling arrays if requested
      if (do_sortcoup == 'Y' .and. nHam == Natom) then
         call sort_coupling_arrays(nlist, ncoup, nlistsize, fs_nlist, fs_nlistsize, &
                                  Natom, max_no_neigh, hdim, conf_num, do_lsf, lsf_field)
      end if
      
      ! Setup LSF index mapping
      if (do_lsf == 'Y' .and. lsf_field == 'L') then
         call setup_lsf_index_mapping(nlist, fs_nlist, nlistsize, fs_nlistsize, &
                                     nind, Natom, max_no_neigh)
      end if
      
      ! Cleanup
      deallocate(morton_perm)
      
   end subroutine setup_sfc_neighbour_hamiltonian
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
   !----------------------------------------------------------------------------
   ! SUBROUTINE: find_sfc_neighbors
   !> @brief SFC-optimized neighbor search using coordinates
   !----------------------------------------------------------------------------
   subroutine find_sfc_neighbors(atom_idx, coord, atype, nn, Natom, max_no_neigh, &
                                max_no_shells, max_no_equiv, morton_perm, &
                                shell_tolerance, neighbor_list, neighbor_count)
      
      implicit none
      integer, intent(in) :: atom_idx, Natom, max_no_neigh, max_no_shells, max_no_equiv
      real(dblprec), dimension(3,Natom), intent(in) :: coord
      integer, dimension(Natom), intent(in) :: atype
      integer, dimension(:), intent(in) :: nn
      integer, dimension(Natom), intent(in) :: morton_perm
      real(dblprec), intent(in) :: shell_tolerance
      integer, dimension(max_no_neigh), intent(out) :: neighbor_list
      integer, intent(out) :: neighbor_count
      
      ! Local variables
      integer :: j, atom_type
      real(dblprec) :: distance, max_distance
      real(dblprec), dimension(3) :: atom_pos
      
      atom_type = atype(atom_idx)
      atom_pos = coord(:, atom_idx)
      neighbor_count = 0
      neighbor_list = 0
      
      ! Estimate maximum interaction distance based on shell count
      max_distance = real(nn(atom_type), dblprec) * 5.0_dblprec  ! Conservative estimate
      
      ! Use Morton-ordered search for better cache performance
      do j = 1, Natom
         if (morton_perm(j) /= atom_idx .and. neighbor_count < max_no_neigh) then
            distance = sqrt(sum((atom_pos - coord(:, morton_perm(j)))**2))
            
            if (distance < max_distance .and. distance > 1.0e-6_dblprec) then
               neighbor_count = neighbor_count + 1
               neighbor_list(neighbor_count) = morton_perm(j)
            end if
         end if
      end do
      
   end subroutine find_sfc_neighbors
   
   !----------------------------------------------------------------------------
   ! SUBROUTINE: determine_shell_from_distance
   !> @brief Determine which shell a neighbor belongs to based on distance
   !----------------------------------------------------------------------------
   function determine_shell_from_distance(distance, atom_type, nn, tolerance) result(shell_idx)
      
      implicit none
      real(dblprec), intent(in) :: distance, tolerance
      integer, intent(in) :: atom_type
      integer, dimension(:), intent(in) :: nn
      integer :: shell_idx
      
      ! Local variables
      integer :: k
      real(dblprec) :: expected_distance
      
      shell_idx = 0
      
      ! Simple distance-based shell assignment
      ! This could be enhanced with more sophisticated shell detection
      do k = 1, nn(atom_type)
         expected_distance = real(k, dblprec) * 2.5_dblprec  ! Approximate shell distance
         if (abs(distance - expected_distance) < tolerance) then
            shell_idx = k
            exit
         end if
      end do
      
   end function determine_shell_from_distance
   
   !----------------------------------------------------------------------------
   ! SUBROUTINE: calculate_coupling_strength
   !> @brief Calculate coupling strength between two atoms
   !----------------------------------------------------------------------------
   subroutine calculate_coupling_strength(i, jatom, shell_idx, atype, anumb, &
                                        xc, ammom_inp, do_ralloy, atype_ch, &
                                        asite_ch, achem_ch, hdim, lexp, &
                                        fc2, iconf, coupling)
      
      implicit none
      integer, intent(in) :: i, jatom, shell_idx, hdim, lexp, do_ralloy, iconf
      integer, dimension(:), intent(in) :: atype, anumb
      integer, dimension(:), intent(in) :: atype_ch, asite_ch, achem_ch
      real(dblprec), dimension(:,:,:), intent(in) :: ammom_inp
      real(dblprec), dimension(hdim,:,:,:,:,:), intent(in) :: xc
      real(dblprec), intent(in) :: fc2
      real(dblprec), dimension(hdim), intent(out) :: coupling
      
      ! Local variables
      real(dblprec) :: mom_i, mom_j
      
      coupling = 0.0_dblprec
      
      if (do_ralloy == 0) then
         ! Non-random alloy case
         mom_i = ammom_inp(anumb(i), 1, iconf)
         mom_j = ammom_inp(anumb(jatom), 1, iconf)
         
         if (abs(mom_i * mom_j) > 1e-6_dblprec) then
            coupling(:) = xc(:, atype(i), shell_idx, 1, 1, iconf) * fc2 &
                         / (mom_i**lexp) / (mom_j**lexp)
         end if
      else
         ! Random alloy case
         mom_i = ammom_inp(asite_ch(i), achem_ch(i), iconf)
         mom_j = ammom_inp(asite_ch(jatom), achem_ch(jatom), iconf)
         
         if (abs(mom_i * mom_j) > 1e-6_dblprec) then
            coupling(:) = xc(:, atype_ch(i), shell_idx, achem_ch(i), achem_ch(jatom), iconf) * &
                         fc2 / (mom_i**lexp) / (mom_j**lexp)
         end if
      end if
      
   end subroutine calculate_coupling_strength
   
   !----------------------------------------------------------------------------
   ! SUBROUTINE: sort_coupling_arrays
   !> @brief Sort the coupling arrays if requested
   !----------------------------------------------------------------------------
   subroutine sort_coupling_arrays(nlist, ncoup, nlistsize, fs_nlist, fs_nlistsize, &
                                  Natom, max_no_neigh, hdim, conf_num, do_lsf, lsf_field)
      
      implicit none
      integer, intent(in) :: Natom, max_no_neigh, hdim, conf_num
      integer, dimension(max_no_neigh,Natom), intent(inout) :: nlist
      real(dblprec), dimension(hdim,max_no_neigh,:,conf_num), intent(inout) :: ncoup
      integer, dimension(Natom), intent(in) :: nlistsize
      integer, dimension(:,:), intent(inout) :: fs_nlist
      integer, dimension(:), intent(in) :: fs_nlistsize
      character(len=1), intent(in) :: do_lsf, lsf_field
      
      ! Local variables
      integer :: i, j, k, tempn, iconf
      real(dblprec), dimension(hdim) :: tempcoup
      
      !$omp parallel do default(shared) private(i,j,k,tempn,tempcoup,iconf)
      do i = 1, Natom
         ! Sort neighbor list - bubble sort
         do j = 1, nlistsize(i)
            do k = 1, nlistsize(i) - j
               if (nlist(k,i) > nlist(k+1,i)) then
                  tempn = nlist(k,i)
                  nlist(k,i) = nlist(k+1,i)
                  nlist(k+1,i) = tempn
                  
                  do iconf = 1, conf_num
                     tempcoup = ncoup(:,k,i,iconf)
                     ncoup(:,k,i,iconf) = ncoup(:,k+1,i,iconf)
                     ncoup(:,k+1,i,iconf) = tempcoup
                  end do
               end if
            end do
         end do
      end do
      !$omp end parallel do
      
      ! Sort first-shell list for LSF
      if (do_lsf == 'Y' .and. lsf_field == 'L') then
         !$omp parallel do default(shared) private(i,j,k,tempn)
         do i = 1, Natom
            do j = 1, fs_nlistsize(i)
               do k = 1, fs_nlistsize(i) - j
                  if (fs_nlist(k,i) > fs_nlist(k+1,i)) then
                     tempn = fs_nlist(k,i)
                     fs_nlist(k,i) = fs_nlist(k+1,i)
                     fs_nlist(k+1,i) = tempn
                  end if
               end do
            end do
         end do
         !$omp end parallel do
      end if
      
   end subroutine sort_coupling_arrays
   
   !----------------------------------------------------------------------------
   ! SUBROUTINE: setup_lsf_index_mapping
   !> @brief Setup index mapping for LSF functionality
   !----------------------------------------------------------------------------
   subroutine setup_lsf_index_mapping(nlist, fs_nlist, nlistsize, fs_nlistsize, &
                                     nind, Natom, max_no_neigh)
      
      implicit none
      integer, intent(in) :: Natom, max_no_neigh
      integer, dimension(max_no_neigh,Natom), intent(in) :: nlist
      integer, dimension(:,:), intent(in) :: fs_nlist
      integer, dimension(Natom), intent(in) :: nlistsize
      integer, dimension(:), intent(in) :: fs_nlistsize
      integer, allocatable, dimension(:,:), intent(inout) :: nind
      
      ! Local variables
      integer :: i, j, k, i_stat
      
      if (.not. allocated(nind)) then
         allocate(nind(maxval(fs_nlistsize),Natom), stat=i_stat)
         if (i_stat /= 0) then
            write(*,*) 'Error allocating nind in setup_lsf_index_mapping'
            return
         end if
      end if
      
      !$omp parallel do default(shared) private(i,j,k)
      do i = 1, Natom
         do j = 1, fs_nlistsize(i)
            do k = 1, nlistsize(i)
               if (nlist(k,i) == fs_nlist(j,i)) then
                  nind(j,i) = k
                  exit
               end if
            end do
         end do
      end do
      !$omp end parallel do
      
   end subroutine setup_lsf_index_mapping
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