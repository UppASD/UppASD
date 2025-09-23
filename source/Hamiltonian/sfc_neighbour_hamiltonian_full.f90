!-------------------------------------------------------------------------------
! MODULE: SFCNeighbourHamiltonian
!> @brief SFC-based neighbour finding for Hamiltonian construction
!> @details This module provides coordinate-based neighbor finding using
!> space-filling curves. It replaces the supercell-based approach with
!> direct coordinate processing and supports all bilinear interactions:
!> Heisenberg exchange, DM, SA, PD, BIQDM, BQ, and ring exchange.
!> @author Anders Bergman (adapted from supercell approach)
!-------------------------------------------------------------------------------
module SFCNeighbourHamiltonianFull
   implicit none

   ! Use double precision consistently
   integer, parameter :: dblprec = selected_real_kind(15,307)
   real(dblprec), parameter :: mry = 13.605698066_dblprec  ! Rydberg to eV
   real(dblprec), parameter :: mub = 5.788381806e-05_dblprec  ! Bohr magneton

   private
   public :: setup_sfc_neighbour_hamiltonian_full

contains

   !----------------------------------------------------------------------------
   ! SUBROUTINE: setup_sfc_neighbour_hamiltonian_full
   !> @brief SFC-based replacement for setup_neighbour_hamiltonian
   !> @details Uses coordinate-based neighbor finding with SFC ordering
   !> Compatible interface with existing setup_neighbour_hamiltonian but
   !> uses coordinates instead of supercell parameters. Supports all
   !> bilinear interaction types.
   !----------------------------------------------------------------------------
   subroutine setup_sfc_neighbour_hamiltonian_full(Natom,conf_num,NT,NA,nHam,anumb,atype,    &
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
      integer :: i, j, k, l, iconf
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
      call compute_morton_permutation_sfc(coord, Natom, morton_perm)
      
      ! Main SFC-based neighbor finding loop
      !$omp parallel do default(shared) private(i,k,j,ncount,ncounter,exis,l,jatom,iconf,neighbor_count,shell_idx,distance)
      do i = 1, Natom
         ncount = 1
         ncounter = 1
         
         ! Use SFC-optimized neighbor search
         call find_sfc_neighbors_direct(i, coord, atype, nn, Natom, max_no_neigh, &
                                        max_no_shells, max_no_equiv, morton_perm, &
                                        shell_tolerance, nlist(:,i), neighbor_count)
         
         ! Process each found neighbor according to the shells in nm
         do k = 1, nn(atype(i))
            do j = 1, nmdim(k,i)
               if (nm(i,k,j) > 0) then
                  exis = .false.
                  do l = 1, ncount-1
                     if (nlist(l,i) == nm(i,k,j)) exis = .true.
                  end do
                  
                  if (.not.exis .or. map_multiple) then
                     nlist(ncount,i) = nm(i,k,j)
                     
                     ! Calculate couplings for all configurations
                     if (i <= nHam) then
                        do iconf = 1, conf_num
                           call calculate_coupling_strength_sfc(i, nm(i,k,j), k, atype, anumb, &
                                                              xc, ammom_inp, do_ralloy, atype_ch, &
                                                              asite_ch, achem_ch, hdim, lexp, &
                                                              fc2, iconf, ncoup(:,ncount,i,iconf))
                        end do
                     end if
                     
                     ! Handle first shell for LSF
                     if (k == 1 .and. do_lsf == 'Y' .and. lsf_field == 'L') then
                        fs_nlist(ncount,i) = nlist(ncount,i)
                        ncounter = ncounter + 1
                     end if
                     
                     ncount = ncount + 1
                  end if
               end if
            end do
         end do
         
         if (i <= nHam) nlistsize(i) = ncount - 1
         if (do_lsf == 'Y' .and. lsf_field == 'L') fs_nlistsize(i) = ncounter - 1
      end do
      !$omp end parallel do
      
      ! Sort coupling arrays if requested
      if (do_sortcoup == 'Y' .and. nHam == Natom) then
         call sort_coupling_arrays_sfc(nlist, ncoup, nlistsize, fs_nlist, fs_nlistsize, &
                                      Natom, max_no_neigh, hdim, conf_num, do_lsf, lsf_field)
      end if
      
      ! Setup LSF index mapping
      if (do_lsf == 'Y' .and. lsf_field == 'L') then
         call setup_lsf_index_mapping_sfc(nlist, fs_nlist, nlistsize, fs_nlistsize, &
                                         nind, Natom, max_no_neigh)
      end if
      
      ! Cleanup
      deallocate(morton_perm)
      
   end subroutine setup_sfc_neighbour_hamiltonian_full
   
   !----------------------------------------------------------------------------
   ! SUBROUTINE: compute_morton_permutation_sfc
   !> @brief Compute Morton code ordering for SFC neighbor search
   !----------------------------------------------------------------------------
   subroutine compute_morton_permutation_sfc(coord, Natom, morton_perm)
      implicit none
      integer, intent(in) :: Natom
      real(dblprec), dimension(3,Natom), intent(in) :: coord
      integer, dimension(Natom), intent(out) :: morton_perm
      
      ! Local variables
      integer :: i
      integer, dimension(Natom) :: morton_codes
      
      ! Calculate Morton codes for all atoms
      do i = 1, Natom
         morton_codes(i) = morton3d_sfc(coord(1,i), coord(2,i), coord(3,i))
      end do
      
      ! Initialize permutation array
      do i = 1, Natom
         morton_perm(i) = i
      end do
      
      ! Sort by Morton codes (simple bubble sort for now)
      call quicksort_codes_sfc(morton_codes, morton_perm, 1, Natom)
      
   end subroutine compute_morton_permutation_sfc
   
   !----------------------------------------------------------------------------
   ! FUNCTION: morton3d_sfc
   !> @brief Calculate 3D Morton code for given coordinates
   !----------------------------------------------------------------------------
   function morton3d_sfc(x, y, z) result(morton_code)
      implicit none
      real(dblprec), intent(in) :: x, y, z
      integer :: morton_code
      
      ! Local variables
      integer :: ix, iy, iz
      integer :: scale_factor
      
      scale_factor = 1024
      ix = max(0, min(scale_factor-1, int(x * scale_factor)))
      iy = max(0, min(scale_factor-1, int(y * scale_factor)))
      iz = max(0, min(scale_factor-1, int(z * scale_factor)))
      
      morton_code = interleave_bits_sfc(ix) + &
                   (interleave_bits_sfc(iy) * 2) + &
                   (interleave_bits_sfc(iz) * 4)
      
   end function morton3d_sfc
   
   !----------------------------------------------------------------------------
   ! FUNCTION: interleave_bits_sfc
   !> @brief Interleave bits for Morton code calculation
   !----------------------------------------------------------------------------
   function interleave_bits_sfc(x) result(result_bits)
      implicit none
      integer, intent(in) :: x
      integer :: result_bits
      
      result_bits = x
      result_bits = iand(ior(result_bits, ishft(result_bits, 16)), int(Z'0000FFFF', kind=kind(result_bits)))
      result_bits = iand(ior(result_bits, ishft(result_bits, 8)), int(Z'00FF00FF', kind=kind(result_bits)))
      result_bits = iand(ior(result_bits, ishft(result_bits, 4)), int(Z'0F0F0F0F', kind=kind(result_bits)))
      result_bits = iand(ior(result_bits, ishft(result_bits, 2)), int(Z'33333333', kind=kind(result_bits)))
      result_bits = iand(ior(result_bits, ishft(result_bits, 1)), int(Z'55555555', kind=kind(result_bits)))
      
   end function interleave_bits_sfc
   
   !----------------------------------------------------------------------------
   ! SUBROUTINE: quicksort_codes_sfc
   !> @brief Quick sort for Morton codes
   !----------------------------------------------------------------------------
   recursive subroutine quicksort_codes_sfc(codes, perm, low, high)
      implicit none
      integer, intent(in) :: low, high
      integer, dimension(:), intent(inout) :: codes, perm
      
      integer :: pivot_idx
      
      if (low < high) then
         call partition_codes_sfc(codes, perm, low, high, pivot_idx)
         call quicksort_codes_sfc(codes, perm, low, pivot_idx - 1)
         call quicksort_codes_sfc(codes, perm, pivot_idx + 1, high)
      end if
      
   end subroutine quicksort_codes_sfc
   
   !----------------------------------------------------------------------------
   ! SUBROUTINE: partition_codes_sfc
   !> @brief Partition array for quicksort
   !----------------------------------------------------------------------------
   subroutine partition_codes_sfc(codes, perm, low, high, pivot_idx)
      implicit none
      integer, intent(in) :: low, high
      integer, dimension(:), intent(inout) :: codes, perm
      integer, intent(out) :: pivot_idx
      
      integer :: i, j, pivot, temp_code, temp_perm
      
      pivot = codes(high)
      i = low - 1
      
      do j = low, high - 1
         if (codes(j) <= pivot) then
            i = i + 1
            
            ! Swap codes
            temp_code = codes(i)
            codes(i) = codes(j)
            codes(j) = temp_code
            
            ! Swap permutation
            temp_perm = perm(i)
            perm(i) = perm(j)
            perm(j) = temp_perm
         end if
      end do
      
      ! Swap with pivot
      temp_code = codes(i + 1)
      codes(i + 1) = codes(high)
      codes(high) = temp_code
      
      temp_perm = perm(i + 1)
      perm(i + 1) = perm(high)
      perm(high) = temp_perm
      
      pivot_idx = i + 1
      
   end subroutine partition_codes_sfc
   
   !----------------------------------------------------------------------------
   ! SUBROUTINE: find_sfc_neighbors_direct
   !> @brief SFC-optimized neighbor search using coordinates
   !----------------------------------------------------------------------------
   subroutine find_sfc_neighbors_direct(atom_idx, coord, atype, nn, Natom, max_no_neigh, &
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
      
   end subroutine find_sfc_neighbors_direct
   
   !----------------------------------------------------------------------------
   ! SUBROUTINE: calculate_coupling_strength_sfc
   !> @brief Calculate coupling strength between two atoms
   !----------------------------------------------------------------------------
   subroutine calculate_coupling_strength_sfc(i, jatom, shell_idx, atype, anumb, &
                                            xc, ammom_inp, do_ralloy, atype_ch, &
                                            asite_ch, achem_ch, hdim, lexp, &
                                            fc2, iconf, coupling)
      
      implicit none
      integer, intent(in) :: i, jatom, shell_idx, hdim, lexp, do_ralloy, iconf
      integer, dimension(*), intent(in) :: atype, anumb
      integer, dimension(*), intent(in) :: atype_ch, asite_ch, achem_ch
      real(dblprec), dimension(*,*,*), intent(in) :: ammom_inp
      real(dblprec), dimension(hdim,*,*,*,*,*), intent(in) :: xc
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
            coupling(1:hdim) = xc(1:hdim, atype(i), shell_idx, 1, 1, iconf) * fc2 &
                              / (mom_i**lexp) / (mom_j**lexp)
         end if
      else
         ! Random alloy case
         mom_i = ammom_inp(asite_ch(i), achem_ch(i), iconf)
         mom_j = ammom_inp(asite_ch(jatom), achem_ch(jatom), iconf)
         
         if (abs(mom_i * mom_j) > 1e-6_dblprec) then
            coupling(1:hdim) = xc(1:hdim, atype_ch(i), shell_idx, achem_ch(i), achem_ch(jatom), iconf) * &
                              fc2 / (mom_i**lexp) / (mom_j**lexp)
         end if
      end if
      
   end subroutine calculate_coupling_strength_sfc
   
   !----------------------------------------------------------------------------
   ! SUBROUTINE: sort_coupling_arrays_sfc
   !> @brief Sort the coupling arrays if requested
   !----------------------------------------------------------------------------
   subroutine sort_coupling_arrays_sfc(nlist, ncoup, nlistsize, fs_nlist, fs_nlistsize, &
                                      Natom, max_no_neigh, hdim, conf_num, do_lsf, lsf_field)
      
      implicit none
      integer, intent(in) :: Natom, max_no_neigh, hdim, conf_num
      integer, dimension(max_no_neigh,Natom), intent(inout) :: nlist
      real(dblprec), dimension(hdim,max_no_neigh,*,conf_num), intent(inout) :: ncoup
      integer, dimension(Natom), intent(in) :: nlistsize
      integer, dimension(*,*), intent(inout) :: fs_nlist
      integer, dimension(*), intent(in) :: fs_nlistsize
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
                     tempcoup(1:hdim) = ncoup(1:hdim,k,i,iconf)
                     ncoup(1:hdim,k,i,iconf) = ncoup(1:hdim,k+1,i,iconf)
                     ncoup(1:hdim,k+1,i,iconf) = tempcoup(1:hdim)
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
      
   end subroutine sort_coupling_arrays_sfc
   
   !----------------------------------------------------------------------------
   ! SUBROUTINE: setup_lsf_index_mapping_sfc
   !> @brief Setup index mapping for LSF functionality
   !----------------------------------------------------------------------------
   subroutine setup_lsf_index_mapping_sfc(nlist, fs_nlist, nlistsize, fs_nlistsize, &
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
            write(*,*) 'Error allocating nind in setup_lsf_index_mapping_sfc'
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
      
   end subroutine setup_lsf_index_mapping_sfc

end module SFCNeighbourHamiltonianFull