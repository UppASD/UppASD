!-------------------------------------------------------------------------------
! MODULE: SFCIntegration
!> @brief Integration helpers for SFC-based neighbor finding
!> @details This module provides utilities to integrate SFC-based neighbor
!> finding with the existing Hamiltonian setup infrastructure.
!> @author Anders Bergman
!-------------------------------------------------------------------------------
module SFCIntegration
   use iso_fortran_env, only: wp => real64
   use SFCNeighbourHamiltonian, only: setup_sfc_neighbour_hamiltonian, setup_sfc_nm
   implicit none

   private
   public :: setup_sfc_hamiltonian_wrapper, convert_to_sfc_coords

contains

   !----------------------------------------------------------------------------
   ! SUBROUTINE: setup_sfc_hamiltonian_wrapper
   !> @brief Wrapper to interface SFC neighbor finding with existing code
   !> @details Converts from supercell-based inputs to coordinate-based SFC approach
   !----------------------------------------------------------------------------
   subroutine setup_sfc_hamiltonian_wrapper(Natom, NT, NA, N1, N2, N3, C1, C2, C3, &
      Bas, atype, anumb, coord, nn, redcoord, conf_num, nHam, &
      max_no_neigh, max_no_shells, max_no_equiv, nlistsize, nlist, ncoup, nm, nmdim, &
      xc, fs_nlist, fs_nlistsize, do_ralloy, Natom_full, Nchmax, &
      atype_ch, asite_ch, achem_ch, ammom_inp, hdim, lexp, do_sortcoup, &
      do_lsf, nind, lsf_field, map_multiple, use_sfc_method)
      
      implicit none
      
      ! System geometry inputs (supercell-based)
      integer, intent(in) :: Natom, NT, NA, N1, N2, N3, conf_num, nHam
      integer, intent(in) :: hdim, lexp, do_ralloy, Natom_full, Nchmax
      real(wp), dimension(3), intent(in) :: C1, C2, C3
      real(wp), dimension(3,NA), intent(in) :: Bas
      real(wp), dimension(3,Natom), intent(in) :: coord
      integer, dimension(Natom), intent(in) :: anumb, atype
      integer, dimension(NT), intent(in) :: nn
      real(wp), dimension(NT,maxval(nn),3), intent(in) :: redcoord
      real(wp), dimension(:,:,:,:,:,:), intent(in) :: xc
      character, intent(in) :: do_sortcoup, do_lsf, lsf_field
      logical, intent(in) :: map_multiple, use_sfc_method
      integer, dimension(Natom_full), intent(in), optional :: atype_ch, asite_ch, achem_ch
      real(wp), dimension(:,:,:), intent(in), optional :: ammom_inp
      
      ! Outputs
      integer, intent(out) :: max_no_neigh, max_no_shells, max_no_equiv
      integer, dimension(Natom), intent(out) :: nlistsize
      integer, dimension(:,:), allocatable, intent(out) :: nlist, nmdim
      integer, dimension(:,:,:), allocatable, intent(out) :: nm
      real(wp), dimension(:,:,:,:), allocatable, intent(out) :: ncoup
      integer, dimension(:,:), intent(out), optional :: fs_nlist
      integer, dimension(:), intent(out), optional :: fs_nlistsize
      integer, dimension(:,:), allocatable, intent(inout), optional :: nind
      
      ! Local variables
      real(wp), dimension(:,:), allocatable :: shell_distances
      real(wp), dimension(:,:,:), allocatable :: xc_shells_dummy
      integer :: itype, ishell
      
      if (use_sfc_method) then
         ! Convert redcoord to shell distances
         allocate(shell_distances(NT, maxval(nn)))
         allocate(xc_shells_dummy(NT, maxval(nn), hdim))
         
         do itype = 1, NT
            do ishell = 1, nn(itype)
               ! Convert reduced coordinates to distances
               shell_distances(itype, ishell) = sqrt(sum(redcoord(itype, ishell, :)**2))
               ! Extract exchange couplings for this shell
               xc_shells_dummy(itype, ishell, :) = xc(:, itype, ishell, 1, 1, 1)
            end do
         end do
         
         ! Call SFC-based neighbor finding
         call setup_sfc_neighbour_hamiltonian(Natom, conf_num, nHam, coord, atype, &
            shell_distances, nn, max_no_neigh, max_no_shells, max_no_equiv, &
            nlistsize, nlist, ncoup, nm, nmdim, xc_shells_dummy, &
            fs_nlist, fs_nlistsize, do_ralloy, Natom_full, &
            Nchmax, atype_ch, asite_ch, achem_ch, ammom_inp, hdim, lexp, do_sortcoup, &
            do_lsf, nind, lsf_field, map_multiple)
         
         deallocate(shell_distances, xc_shells_dummy)
         
      else
         ! Fall back to original supercell-based method
         ! This would call the original setup_neighbour_hamiltonian
         ! For now, just indicate this path was taken
         write(*,*) 'SFC Integration: Falling back to original supercell method'
         stop 'Original method not implemented in this wrapper'
      endif
      
   end subroutine setup_sfc_hamiltonian_wrapper

   !----------------------------------------------------------------------------
   ! SUBROUTINE: convert_to_sfc_coords
   !> @brief Convert supercell geometry to coordinate array
   !> @details Helper function to build coordinate array from supercell parameters
   !----------------------------------------------------------------------------
   subroutine convert_to_sfc_coords(N1, N2, N3, NA, C1, C2, C3, Bas, coords)
      implicit none
      integer, intent(in) :: N1, N2, N3, NA
      real(wp), dimension(3), intent(in) :: C1, C2, C3
      real(wp), dimension(3,NA), intent(in) :: Bas
      real(wp), dimension(3,N1*N2*N3*NA), intent(out) :: coords
      
      integer :: ix, iy, iz, ia, idx
      
      idx = 0
      do iz = 0, N3-1
         do iy = 0, N2-1
            do ix = 0, N1-1
               do ia = 1, NA
                  idx = idx + 1
                  coords(1,idx) = Bas(1,ia) + real(ix,wp)*C1(1) + real(iy,wp)*C2(1) + real(iz,wp)*C3(1)
                  coords(2,idx) = Bas(2,ia) + real(ix,wp)*C1(2) + real(iy,wp)*C2(2) + real(iz,wp)*C3(2)
                  coords(3,idx) = Bas(3,ia) + real(ix,wp)*C1(3) + real(iy,wp)*C2(3) + real(iz,wp)*C3(3)
               end do
            end do
         end do
      end do
   end subroutine convert_to_sfc_coords

end module SFCIntegration