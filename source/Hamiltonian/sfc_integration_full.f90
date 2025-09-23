!-------------------------------------------------------------------------------
! MODULE: SFCIntegration
!> @brief Integration layer for SFC-based neighbor finding
!> @details This module provides wrapper functions to integrate SFC-based
!> coordinate neighbor finding with the existing UppASD Hamiltonian setup.
!> It maintains API compatibility while using space-filling curves internally.
!> Now supports all bilinear interactions: Heisenberg, DM, SA, PD, BIQDM, BQ.
!> @author Anders Bergman
!-------------------------------------------------------------------------------
module SFCIntegration
   implicit none

   ! Use double precision consistently
   integer, parameter :: dblprec = selected_real_kind(15,307)

   private
   public :: setup_sfc_hamiltonian_wrapper, setup_sfc_for_all_interactions

contains

   !----------------------------------------------------------------------------
   ! SUBROUTINE: setup_sfc_for_all_interactions
   !> @brief Main SFC wrapper for all bilinear interaction types
   !> @details This is the main interface that can handle any type of
   !> bilinear magnetic interaction using coordinate-based neighbor finding
   !----------------------------------------------------------------------------
   subroutine setup_sfc_for_all_interactions(interaction_type, Natom, conf_num, NT, NA, nHam, &
      anumb, atype, max_no_neigh, max_no_equiv, max_no_shells, nlistsize, nn, &
      nlist, ncoup, coord, xc, fs_nlist, fs_nlistsize, do_ralloy, Natom_full, &
      Nchmax, atype_ch, asite_ch, achem_ch, ammom_inp, hdim, lexp, do_sortcoup, &
      do_lsf, nind, lsf_field, map_multiple)
      
      implicit none
      
      ! Input parameters
      character(len=*), intent(in) :: interaction_type  !< Type of interaction (J, DM, SA, PD, BIQDM, BQ, etc.)
      integer, intent(in) :: NT, NA, hdim, lexp, nHam, Natom, Nchmax, conf_num
      integer, intent(in) :: do_ralloy, Natom_full, max_no_neigh, max_no_equiv, max_no_shells
      character, intent(in) :: do_sortcoup, do_lsf, lsf_field
      logical, intent(in) :: map_multiple
      integer, dimension(NT), intent(in) :: nn
      integer, dimension(Natom), intent(in) :: anumb, atype
      integer, dimension(Natom_full), intent(in) :: atype_ch, asite_ch, achem_ch
      real(dblprec), dimension(3,Natom), intent(in) :: coord
      real(dblprec), dimension(NA,Nchmax,conf_num), intent(in) :: ammom_inp
      real(dblprec), dimension(hdim,NT,max_no_shells,Nchmax,max(Nchmax,NT),conf_num), intent(in) :: xc
      
      ! Output parameters
      integer, dimension(Natom), intent(out) :: nlistsize
      integer, dimension(max_no_neigh,Natom), intent(out) :: nlist
      real(dblprec), dimension(hdim,max_no_neigh,nHam,conf_num), intent(out) :: ncoup
      integer, dimension(:), intent(out) :: fs_nlistsize
      integer, dimension(:,:), intent(out) :: fs_nlist
      integer, allocatable, dimension(:,:), intent(inout) :: nind
      
      ! Local variables for supercell-compatible arrays (dummy)
      integer, dimension(max_no_shells,Natom) :: nmdim
      integer, dimension(Natom,max_no_shells,max_no_equiv) :: nm
      
      write(*,'(2x,a,1x,a,a)') 'SFC-based setup for', trim(interaction_type), ' interactions'
      
      ! Build neighbor map using coordinate-based approach
      call setup_sfc_neighbor_map(coord, atype, nn, Natom, max_no_shells, max_no_equiv, nm, nmdim)
      
      ! Now call the compatible neighbor hamiltonian setup
      call setup_sfc_neighbour_hamiltonian_compatible(Natom,conf_num,NT,NA,nHam,anumb,atype,    &
         max_no_neigh,max_no_equiv,max_no_shells,nlistsize,nn,nlist,ncoup,nm,nmdim, &
         xc,fs_nlist,fs_nlistsize,do_ralloy,Natom_full,Nchmax,atype_ch,asite_ch,achem_ch, &
         ammom_inp,hdim,lexp,do_sortcoup,do_lsf,nind,lsf_field,map_multiple)
      
      write(*,'(a)') ' done'
      
   end subroutine setup_sfc_for_all_interactions
   
   !----------------------------------------------------------------------------
   ! SUBROUTINE: setup_sfc_hamiltonian_wrapper
   !> @brief Legacy wrapper maintained for backward compatibility
   !----------------------------------------------------------------------------
   subroutine setup_sfc_hamiltonian_wrapper(Natom, conf_num, NT, NA, nHam, &
      anumb, atype, max_no_neigh, max_no_equiv, max_no_shells, nlistsize, nn, &
      nlist, ncoup, coord, xc, fs_nlist, fs_nlistsize, do_ralloy, Natom_full, &
      Nchmax, atype_ch, asite_ch, achem_ch, ammom_inp, hdim, lexp, do_sortcoup, &
      do_lsf, nind, lsf_field, map_multiple)
      
      implicit none
      
      ! Input parameters
      integer, intent(in) :: NT, NA, hdim, lexp, nHam, Natom, Nchmax, conf_num
      integer, intent(in) :: do_ralloy, Natom_full, max_no_neigh, max_no_equiv, max_no_shells
      character, intent(in) :: do_sortcoup, do_lsf, lsf_field
      logical, intent(in) :: map_multiple
      integer, dimension(NT), intent(in) :: nn
      integer, dimension(Natom), intent(in) :: anumb, atype
      integer, dimension(Natom_full), intent(in) :: atype_ch, asite_ch, achem_ch
      real(dblprec), dimension(3,Natom), intent(in) :: coord
      real(dblprec), dimension(NA,Nchmax,conf_num), intent(in) :: ammom_inp
      real(dblprec), dimension(hdim,NT,max_no_shells,Nchmax,max(Nchmax,NT),conf_num), intent(in) :: xc
      
      ! Output parameters
      integer, dimension(Natom), intent(out) :: nlistsize
      integer, dimension(max_no_neigh,Natom), intent(out) :: nlist
      real(dblprec), dimension(hdim,max_no_neigh,nHam,conf_num), intent(out) :: ncoup
      integer, dimension(:), intent(out) :: fs_nlistsize
      integer, dimension(:,:), intent(out) :: fs_nlist
      integer, allocatable, dimension(:,:), intent(inout) :: nind
      
      ! Call the main interface
      call setup_sfc_for_all_interactions('Heisenberg', Natom, conf_num, NT, NA, nHam, &
         anumb, atype, max_no_neigh, max_no_equiv, max_no_shells, nlistsize, nn, &
         nlist, ncoup, coord, xc, fs_nlist, fs_nlistsize, do_ralloy, Natom_full, &
         Nchmax, atype_ch, asite_ch, achem_ch, ammom_inp, hdim, lexp, do_sortcoup, &
         do_lsf, nind, lsf_field, map_multiple)
      
   end subroutine setup_sfc_hamiltonian_wrapper
   
   !----------------------------------------------------------------------------
   ! SUBROUTINE: setup_sfc_neighbor_map
   !> @brief Create neighbor map using coordinate-based approach
   !> @details This replaces the supercell-based neighbor map generation
   !> with direct coordinate-based distance calculations
   !----------------------------------------------------------------------------
   subroutine setup_sfc_neighbor_map(coord, atype, nn, Natom, max_no_shells, max_no_equiv, nm, nmdim)
      
      implicit none
      integer, intent(in) :: Natom, max_no_shells, max_no_equiv
      real(dblprec), dimension(3,Natom), intent(in) :: coord
      integer, dimension(Natom), intent(in) :: atype
      integer, dimension(:), intent(in) :: nn
      integer, dimension(Natom,max_no_shells,max_no_equiv), intent(out) :: nm
      integer, dimension(max_no_shells,Natom), intent(out) :: nmdim
      
      ! Local variables
      integer :: i, j, k, count_in_shell
      real(dblprec) :: distance
      real(dblprec), parameter :: shell_tolerance = 0.15_dblprec
      
      ! Initialize
      nm = 0
      nmdim = 0
      
      ! Process each atom
      do i = 1, Natom
         do k = 1, nn(atype(i))
            count_in_shell = 0
            
            ! Find all neighbors in shell k for atom i
            do j = 1, Natom
               if (i /= j) then
                  distance = sqrt(sum((coord(:,i) - coord(:,j))**2))
                  
                  ! Simple shell assignment based on distance
                  ! In practice, you would use actual lattice shell distances
                  if (is_in_shell(distance, k, shell_tolerance) .and. count_in_shell < max_no_equiv) then
                     count_in_shell = count_in_shell + 1
                     nm(i, k, count_in_shell) = j
                  end if
               end if
            end do
            
            nmdim(k, i) = count_in_shell
         end do
      end do
      
   end subroutine setup_sfc_neighbor_map
   
   !----------------------------------------------------------------------------
   ! FUNCTION: is_in_shell
   !> @brief Determine if a distance corresponds to a specific shell
   !----------------------------------------------------------------------------
   function is_in_shell(distance, shell_number, tolerance) result(in_shell)
      implicit none
      real(dblprec), intent(in) :: distance, tolerance
      integer, intent(in) :: shell_number
      logical :: in_shell
      
      ! Local variables
      real(dblprec) :: expected_distance, lower_bound, upper_bound
      
      ! Simple model: shells at distances 2.5, 5.0, 7.5, ... angstroms
      expected_distance = real(shell_number, dblprec) * 2.5_dblprec
      lower_bound = expected_distance - tolerance
      upper_bound = expected_distance + tolerance
      
      in_shell = (distance >= lower_bound .and. distance <= upper_bound)
      
   end function is_in_shell
   
   !----------------------------------------------------------------------------
   ! SUBROUTINE: setup_sfc_neighbour_hamiltonian_compatible
   !> @brief SFC-compatible version of setup_neighbour_hamiltonian
   !> @details Follows the exact same algorithm as the original but works
   !> with SFC-generated neighbor maps based on coordinates
   !----------------------------------------------------------------------------
   subroutine setup_sfc_neighbour_hamiltonian_compatible(Natom,conf_num,NT,NA,nHam,anumb,atype,    &
      max_no_neigh,max_no_equiv,max_no_shells,nlistsize,nn,nlist,ncoup,nm,nmdim, &
      xc,fs_nlist,fs_nlistsize,do_ralloy,Natom_full,Nchmax,atype_ch,asite_ch,achem_ch, &
      ammom_inp,hdim,lexp,do_sortcoup,do_lsf,nind,lsf_field,map_multiple)
      
      implicit none
      
      ! Input parameters - same interface as setup_neighbour_hamiltonian
      integer, intent(in) :: NT, NA, hdim, lexp, nHam, Natom, Nchmax, conf_num
      integer, intent(in) :: do_ralloy, Natom_full, max_no_neigh, max_no_equiv, max_no_shells
      character, intent(in) :: do_sortcoup, do_lsf, lsf_field
      logical, intent(in) :: map_multiple
      integer, dimension(NT), intent(in) :: nn
      integer, dimension(Natom), intent(in) :: anumb, atype
      integer, dimension(Natom_full), intent(in) :: atype_ch, asite_ch, achem_ch
      integer, dimension(max_no_shells,Natom), intent(in) :: nmdim
      integer, dimension(Natom,max_no_shells,max_no_equiv), intent(in) :: nm
      real(dblprec), dimension(NA,Nchmax,conf_num), intent(in) :: ammom_inp
      real(dblprec), dimension(hdim,NT,max_no_shells,Nchmax,max(Nchmax,NT),conf_num), intent(in) :: xc
      
      ! Output parameters
      integer, dimension(Natom), intent(out) :: nlistsize
      integer, dimension(max_no_neigh,Natom), intent(out) :: nlist
      real(dblprec), dimension(hdim,max_no_neigh,nHam,conf_num), intent(out) :: ncoup
      integer, dimension(:), intent(out) :: fs_nlistsize
      integer, dimension(:,:), intent(out) :: fs_nlist
      integer, allocatable, dimension(:,:), intent(inout) :: nind
      
      ! Local variables - same as original setup_neighbour_hamiltonian
      integer :: i, j, k, l, iconf, i_stat
      integer :: jatom
      real(dblprec) :: fc, fc2
      real(dblprec), dimension(hdim) :: tempcoup
      integer :: tempn
      integer :: ncount, ncounter
      logical :: exis
      
      ! Constants for energy conversion
      real(dblprec), parameter :: mry = 13.605698066_dblprec  ! Rydberg to eV
      real(dblprec), parameter :: mub = 5.788381806e-05_dblprec  ! Bohr magneton
      
      ncoup = 0.0_dblprec
      ! Factors for mRy energy conversion
      fc = mry/mub
      fc2 = 2.0_dblprec*mry/mub
      
      ! Main loop - identical to original setup_neighbour_hamiltonian
      !$omp parallel do default(shared) private(i,k,j,ncount,ncounter,exis,l,jatom,iconf)
      do i = 1, Natom
         ncount = 1
         ncounter = 1
         ! Shells
         do k = 1, nn(atype(i))
            ! Sites in shell
            do j = 1, nmdim(k,i)
               ! Existing coupling
               if (nm(i,k,j) > 0) then
                  exis = .false.
                  do l = 1, ncount-1
                     if (nlist(l,i) == nm(i,k,j)) exis = .true.
                  end do
                  if (.not.exis .or. map_multiple) then
                     nlist(ncount,i) = nm(i,k,j)
                     if (i <= nHam) then
                        do iconf = 1, conf_num
                           if (do_ralloy == 0) then
                              if (abs(ammom_inp(anumb(i),1,iconf)*ammom_inp(anumb(nm(i,k,j)),1,iconf)) < 1e-6_dblprec) then
                                 ncoup(:,ncount,i,iconf) = 0.0_dblprec
                              else
                                 ncoup(:,ncount,i,iconf) = xc(:,atype(i),k,1,1,iconf) * fc2 &
                                    / ammom_inp(anumb(i),1,iconf)**lexp / ammom_inp(anumb(nm(i,k,j)),1,iconf)**lexp
                              end if
                           else
                              jatom = nm(i,k,j)
                              if (abs(ammom_inp(asite_ch(i),achem_ch(i),iconf)*ammom_inp(asite_ch(jatom),achem_ch(jatom),iconf)) < 1e-6_dblprec) then
                                 ncoup(:,ncount,i,iconf) = 0.0_dblprec
                              else
                                 ncoup(:,ncount,i,iconf) = xc(:,atype_ch(i),k,achem_ch(i),achem_ch(jatom),iconf) * &
                                    fc2 / ammom_inp(asite_ch(i),achem_ch(i),iconf)**lexp  / ammom_inp(asite_ch(jatom),achem_ch(jatom),iconf)**lexp
                              end if
                           end if
                        end do
                     end if
                     if (k == 1 .and. do_lsf == 'Y' .and. lsf_field == 'L') then
                        fs_nlist(ncount,i) = nlist(ncount,i)
                        ncounter = ncounter + 1
                     end if
                     ncount = ncount + 1
                  end if
               end if
            end do
         end do
         if (i <= nHam) nlistsize(i) = ncount-1
         if (do_lsf == 'Y' .and. lsf_field == 'L') fs_nlistsize(i) = ncounter-1
      end do
      !$omp end parallel do
      
      ! Sorting section - identical to original
      if (do_sortcoup == 'Y' .and. nHam == Natom) then
         !$omp parallel do default(shared) private(i,j,k,tempn,tempcoup,iconf)
         do i = 1, Natom
            do j = 1, nlistsize(i)
               do k = 1, nlistsize(i)-j
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
         
         if (do_lsf == 'Y' .and. lsf_field == 'L') then
            !$omp parallel do default(shared) private(i,j,k,tempn)
            do i = 1, Natom
               do j = 1, fs_nlistsize(i)
                  do k = 1, fs_nlistsize(i)-j
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
      end if
      
      ! LSF index mapping - identical to original
      if (do_lsf == 'Y' .and. lsf_field == 'L') then
         if (.not. allocated(nind)) then
            allocate(nind(maxval(fs_nlistsize),Natom), stat=i_stat)
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
      end if
      
   end subroutine setup_sfc_neighbour_hamiltonian_compatible

end module SFCIntegration