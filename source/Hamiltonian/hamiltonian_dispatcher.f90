!-------------------------------------------------------------------------------
! MODULE: HamiltonianSetupDispatcher
!> @brief Dispatcher module for choosing between traditional and SFC Hamiltonian setup
!> @details This module provides a unified interface that can dispatch between
!> the traditional supercell-based Hamiltonian setup and the new SFC-based
!> coordinate approach. The choice is controlled by input parameters.
!> @author Anders Bergman
!-------------------------------------------------------------------------------
module HamiltonianSetupDispatcher
   use Parameters
   use Constants
   use SFCIntegration
   use SFCHamiltonianInterface
   ! NOTE: HamiltonianInit would be imported here for traditional setup
   ! use HamiltonianInit, only: setup_neighbour_hamiltonian
   implicit none

   private
   public :: setup_hamiltonian_dispatch, setup_neighbour_hamiltonian_dispatch
   public :: use_sfc_setup

   ! Module variables for configuration
   logical :: use_sfc_setup = .false.  !< Flag to use SFC-based setup instead of traditional

contains

   !----------------------------------------------------------------------------
   ! SUBROUTINE: setup_hamiltonian_dispatch
   !> @brief Main dispatcher for Hamiltonian setup
   !> @details Chooses between traditional supercell and SFC coordinate-based
   !> Hamiltonian setup based on configuration
   !----------------------------------------------------------------------------
   subroutine setup_hamiltonian_dispatch(interaction_type, Natom, conf_num, NT, NA, nHam, &
      anumb, atype, max_no_neigh, max_no_equiv, max_no_shells, nlistsize, nn, &
      nlist, ncoup, coord, xc, fs_nlist, fs_nlistsize, do_ralloy, Natom_full, &
      Nchmax, atype_ch, asite_ch, achem_ch, ammom_inp, hdim, lexp, do_sortcoup, &
      do_lsf, nind, lsf_field, map_multiple, use_sfc)
      
      implicit none
      
      ! Input parameters
      character(len=*), intent(in) :: interaction_type
      integer, intent(in) :: NT, NA, hdim, lexp, nHam, Natom, Nchmax, conf_num
      integer, intent(in) :: do_ralloy, Natom_full, max_no_neigh, max_no_equiv, max_no_shells
      character, intent(in) :: do_sortcoup, do_lsf, lsf_field
      logical, intent(in) :: map_multiple
      logical, intent(in), optional :: use_sfc  !< Force SFC usage if present
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
      
      ! Local variables
      logical :: use_sfc_local
      
      ! Determine which setup method to use
      if (present(use_sfc)) then
         use_sfc_local = use_sfc
      else
         use_sfc_local = use_sfc_setup
      end if
      
      if (use_sfc_local) then
         write(*,'(a)') 'Using SFC-based Hamiltonian setup'
         call setup_sfc_all_interactions(interaction_type, Natom, conf_num, NT, NA, nHam, &
            anumb, atype, max_no_neigh, max_no_equiv, max_no_shells, nlistsize, nn, &
            nlist, ncoup, coord, xc, fs_nlist, fs_nlistsize, do_ralloy, Natom_full, &
            Nchmax, atype_ch, asite_ch, achem_ch, ammom_inp, hdim, lexp, do_sortcoup, &
            do_lsf, nind, lsf_field, map_multiple)
      else
         write(*,'(a)') 'Using traditional supercell-based Hamiltonian setup'
         
         ! NOTE: Traditional setup requires supercell neighbor map (nm, nmdim)
         ! which must be computed from supercell lattice parameters (N1,N2,N3)
         ! The traditional path would be:
         ! 1. call setup_neighbour_hamiltonian(Natom,conf_num,NT,NA,nHam,anumb,atype, &
         !       max_no_neigh,max_no_equiv,max_no_shells,nlistsize,nn,nlist,ncoup,nm,nmdim,xc, &
         !       fs_nlist,fs_nlistsize,do_ralloy,Natom_full,Nchmax,atype_ch,asite_ch,achem_ch, &
         !       ammom_inp,hdim,lexp,do_sortcoup,do_lsf,nind,lsf_field,map_multiple)
         ! But this requires nm,nmdim from prior supercell setup_hamiltonian call
         
         write(*,'(a)') 'WARNING: Traditional setup requires supercell neighbor map!'
         write(*,'(a)') 'Falling back to SFC setup for now.'
         call setup_sfc_all_interactions(interaction_type, Natom, conf_num, NT, NA, nHam, &
            anumb, atype, max_no_neigh, max_no_equiv, max_no_shells, nlistsize, nn, &
            nlist, ncoup, coord, xc, fs_nlist, fs_nlistsize, do_ralloy, Natom_full, &
            Nchmax, atype_ch, asite_ch, achem_ch, ammom_inp, hdim, lexp, do_sortcoup, &
            do_lsf, nind, lsf_field, map_multiple)
      end if
      
   end subroutine setup_hamiltonian_dispatch

   !----------------------------------------------------------------------------
   ! SUBROUTINE: setup_neighbour_hamiltonian_dispatch
   !> @brief Dispatcher specifically for neighbor-based Hamiltonian setup
   !> @details Legacy interface maintained for backward compatibility
   !----------------------------------------------------------------------------
   subroutine setup_neighbour_hamiltonian_dispatch(Natom, conf_num, NT, NA, nHam, &
      anumb, atype, max_no_neigh, max_no_equiv, max_no_shells, nlistsize, nn, &
      nlist, ncoup, coord, xc, fs_nlist, fs_nlistsize, do_ralloy, Natom_full, &
      Nchmax, atype_ch, asite_ch, achem_ch, ammom_inp, hdim, lexp, do_sortcoup, &
      do_lsf, nind, lsf_field, map_multiple)
      
      implicit none
      
      ! Input parameters (same as setup_hamiltonian_dispatch but without interaction_type)
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
      
      ! Determine interaction type based on hdim
      character(len=20) :: interaction_type
      
      select case (hdim)
      case (1)
         if (lexp == 1) then
            interaction_type = 'Heisenberg'
         else
            interaction_type = 'Biquadratic'
         end if
      case (3)
         interaction_type = 'DM'  ! Could be DM or SA - assume DM for now
      case (9)
         interaction_type = 'PD'
      case default
         interaction_type = 'Heisenberg'
      end select
      
      ! Call the main dispatcher
      call setup_hamiltonian_dispatch(interaction_type, Natom, conf_num, NT, NA, nHam, &
         anumb, atype, max_no_neigh, max_no_equiv, max_no_shells, nlistsize, nn, &
         nlist, ncoup, coord, xc, fs_nlist, fs_nlistsize, do_ralloy, Natom_full, &
         Nchmax, atype_ch, asite_ch, achem_ch, ammom_inp, hdim, lexp, do_sortcoup, &
         do_lsf, nind, lsf_field, map_multiple)
      
   end subroutine setup_neighbour_hamiltonian_dispatch

   !----------------------------------------------------------------------------
   ! SUBROUTINE: set_sfc_mode
   !> @brief Configure whether to use SFC-based setup
   !----------------------------------------------------------------------------
   subroutine set_sfc_mode(use_sfc)
      implicit none
      logical, intent(in) :: use_sfc
      
      use_sfc_setup = use_sfc
      if (use_sfc) then
         write(*,'(a)') 'SFC-based Hamiltonian setup enabled'
      else
         write(*,'(a)') 'Traditional supercell-based Hamiltonian setup enabled'
      end if
      
   end subroutine set_sfc_mode

   !----------------------------------------------------------------------------
   ! FUNCTION: get_sfc_mode
   !> @brief Check if SFC-based setup is enabled
   !----------------------------------------------------------------------------
   function get_sfc_mode() result(is_enabled)
      implicit none
      logical :: is_enabled
      
      is_enabled = use_sfc_setup
      
   end function get_sfc_mode

end module HamiltonianSetupDispatcher