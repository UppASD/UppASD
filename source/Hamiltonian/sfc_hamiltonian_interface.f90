!-------------------------------------------------------------------------------
! MODULE: SFCHamiltonianInterface
!> @brief Complete interface for SFC-based Hamiltonian setup supporting all interaction types
!> @details This module provides the main entry points for setting up all types of
!> bilinear magnetic interactions using SFC-based coordinate neighbor finding.
!> It supports: Heisenberg (J), Dzyaloshinskii-Moriya (DM), Symmetric Anisotropic (SA),
!> Pseudo-Dipolar (PD), BIQDM, Biquadratic (BQ), and Ring exchange interactions.
!> @author Anders Bergman
!-------------------------------------------------------------------------------
module SFCHamiltonianInterface
   use SFCIntegration
   implicit none

   ! Use double precision consistently
   integer, parameter :: dblprec = selected_real_kind(15,307)

   private
   public :: setup_sfc_heisenberg, setup_sfc_dm_interaction, setup_sfc_sa_interaction
   public :: setup_sfc_pd_interaction, setup_sfc_biqdm_interaction, setup_sfc_bq_interaction
   public :: setup_sfc_ring_exchange, setup_sfc_all_interactions

contains

   !----------------------------------------------------------------------------
   ! SUBROUTINE: setup_sfc_heisenberg
   !> @brief Setup Heisenberg exchange interactions using SFC (hdim=1)
   !----------------------------------------------------------------------------
   subroutine setup_sfc_heisenberg(Natom, conf_num, NT, NA, nHam, &
      anumb, atype, max_no_neigh, max_no_equiv, max_no_shells, nlistsize, nn, &
      nlist, ncoup, coord, xc, fs_nlist, fs_nlistsize, do_ralloy, Natom_full, &
      Nchmax, atype_ch, asite_ch, achem_ch, ammom_inp, do_sortcoup, &
      do_lsf, nind, lsf_field, map_multiple)
      
      implicit none
      
      ! Input parameters
      integer, intent(in) :: NT, NA, nHam, Natom, Nchmax, conf_num
      integer, intent(in) :: do_ralloy, Natom_full, max_no_neigh, max_no_equiv, max_no_shells
      character, intent(in) :: do_sortcoup, do_lsf, lsf_field
      logical, intent(in) :: map_multiple
      integer, dimension(NT), intent(in) :: nn
      integer, dimension(Natom), intent(in) :: anumb, atype
      integer, dimension(Natom_full), intent(in) :: atype_ch, asite_ch, achem_ch
      real(dblprec), dimension(3,Natom), intent(in) :: coord
      real(dblprec), dimension(NA,Nchmax,conf_num), intent(in) :: ammom_inp
      real(dblprec), dimension(1,NT,max_no_shells,Nchmax,max(Nchmax,NT),conf_num), intent(in) :: xc
      
      ! Output parameters
      integer, dimension(Natom), intent(out) :: nlistsize
      integer, dimension(max_no_neigh,Natom), intent(out) :: nlist
      real(dblprec), dimension(1,max_no_neigh,nHam,conf_num), intent(out) :: ncoup
      integer, dimension(:), intent(out) :: fs_nlistsize
      integer, dimension(:,:), intent(out) :: fs_nlist
      integer, allocatable, dimension(:,:), intent(inout) :: nind
      
      write(*,'(a)') 'Setting up Heisenberg exchange with SFC neighbor finding...'
      
      call setup_sfc_for_all_interactions('Heisenberg', Natom, conf_num, NT, NA, nHam, &
         anumb, atype, max_no_neigh, max_no_equiv, max_no_shells, nlistsize, nn, &
         nlist, ncoup, coord, xc, fs_nlist, fs_nlistsize, do_ralloy, Natom_full, &
         Nchmax, atype_ch, asite_ch, achem_ch, ammom_inp, 1, 1, do_sortcoup, &
         do_lsf, nind, lsf_field, map_multiple)
      
   end subroutine setup_sfc_heisenberg

   !----------------------------------------------------------------------------
   ! SUBROUTINE: setup_sfc_dm_interaction
   !> @brief Setup Dzyaloshinskii-Moriya interactions using SFC (hdim=3)
   !----------------------------------------------------------------------------
   subroutine setup_sfc_dm_interaction(Natom, conf_num, NT, NA, nHam, &
      anumb, atype, max_no_neigh, max_no_equiv, max_no_shells, nlistsize, nn, &
      nlist, ncoup, coord, xc, fs_nlist, fs_nlistsize, do_ralloy, Natom_full, &
      Nchmax, atype_ch, asite_ch, achem_ch, ammom_inp, do_sortcoup, &
      do_lsf, nind, lsf_field, map_multiple)
      
      implicit none
      
      ! Input parameters
      integer, intent(in) :: NT, NA, nHam, Natom, Nchmax, conf_num
      integer, intent(in) :: do_ralloy, Natom_full, max_no_neigh, max_no_equiv, max_no_shells
      character, intent(in) :: do_sortcoup, do_lsf, lsf_field
      logical, intent(in) :: map_multiple
      integer, dimension(NT), intent(in) :: nn
      integer, dimension(Natom), intent(in) :: anumb, atype
      integer, dimension(Natom_full), intent(in) :: atype_ch, asite_ch, achem_ch
      real(dblprec), dimension(3,Natom), intent(in) :: coord
      real(dblprec), dimension(NA,Nchmax,conf_num), intent(in) :: ammom_inp
      real(dblprec), dimension(3,NT,max_no_shells,Nchmax,max(Nchmax,NT),conf_num), intent(in) :: xc
      
      ! Output parameters
      integer, dimension(Natom), intent(out) :: nlistsize
      integer, dimension(max_no_neigh,Natom), intent(out) :: nlist
      real(dblprec), dimension(3,max_no_neigh,nHam,conf_num), intent(out) :: ncoup
      integer, dimension(:), intent(out) :: fs_nlistsize
      integer, dimension(:,:), intent(out) :: fs_nlist
      integer, allocatable, dimension(:,:), intent(inout) :: nind
      
      write(*,'(a)') 'Setting up Dzyaloshinskii-Moriya interactions with SFC neighbor finding...'
      
      call setup_sfc_for_all_interactions('DM', Natom, conf_num, NT, NA, nHam, &
         anumb, atype, max_no_neigh, max_no_equiv, max_no_shells, nlistsize, nn, &
         nlist, ncoup, coord, xc, fs_nlist, fs_nlistsize, do_ralloy, Natom_full, &
         Nchmax, atype_ch, asite_ch, achem_ch, ammom_inp, 3, 1, do_sortcoup, &
         do_lsf, nind, lsf_field, map_multiple)
      
   end subroutine setup_sfc_dm_interaction

   !----------------------------------------------------------------------------
   ! SUBROUTINE: setup_sfc_sa_interaction
   !> @brief Setup Symmetric Anisotropic interactions using SFC (hdim=3)
   !----------------------------------------------------------------------------
   subroutine setup_sfc_sa_interaction(Natom, conf_num, NT, NA, nHam, &
      anumb, atype, max_no_neigh, max_no_equiv, max_no_shells, nlistsize, nn, &
      nlist, ncoup, coord, xc, fs_nlist, fs_nlistsize, do_ralloy, Natom_full, &
      Nchmax, atype_ch, asite_ch, achem_ch, ammom_inp, do_sortcoup, &
      do_lsf, nind, lsf_field, map_multiple)
      
      implicit none
      
      ! Input parameters
      integer, intent(in) :: NT, NA, nHam, Natom, Nchmax, conf_num
      integer, intent(in) :: do_ralloy, Natom_full, max_no_neigh, max_no_equiv, max_no_shells
      character, intent(in) :: do_sortcoup, do_lsf, lsf_field
      logical, intent(in) :: map_multiple
      integer, dimension(NT), intent(in) :: nn
      integer, dimension(Natom), intent(in) :: anumb, atype
      integer, dimension(Natom_full), intent(in) :: atype_ch, asite_ch, achem_ch
      real(dblprec), dimension(3,Natom), intent(in) :: coord
      real(dblprec), dimension(NA,Nchmax,conf_num), intent(in) :: ammom_inp
      real(dblprec), dimension(3,NT,max_no_shells,Nchmax,max(Nchmax,NT),conf_num), intent(in) :: xc
      
      ! Output parameters
      integer, dimension(Natom), intent(out) :: nlistsize
      integer, dimension(max_no_neigh,Natom), intent(out) :: nlist
      real(dblprec), dimension(3,max_no_neigh,nHam,conf_num), intent(out) :: ncoup
      integer, dimension(:), intent(out) :: fs_nlistsize
      integer, dimension(:,:), intent(out) :: fs_nlist
      integer, allocatable, dimension(:,:), intent(inout) :: nind
      
      write(*,'(a)') 'Setting up Symmetric Anisotropic interactions with SFC neighbor finding...'
      
      call setup_sfc_for_all_interactions('SA', Natom, conf_num, NT, NA, nHam, &
         anumb, atype, max_no_neigh, max_no_equiv, max_no_shells, nlistsize, nn, &
         nlist, ncoup, coord, xc, fs_nlist, fs_nlistsize, do_ralloy, Natom_full, &
         Nchmax, atype_ch, asite_ch, achem_ch, ammom_inp, 3, 1, do_sortcoup, &
         do_lsf, nind, lsf_field, map_multiple)
      
   end subroutine setup_sfc_sa_interaction

   !----------------------------------------------------------------------------
   ! SUBROUTINE: setup_sfc_pd_interaction
   !> @brief Setup Pseudo-Dipolar interactions using SFC (hdim=9)
   !----------------------------------------------------------------------------
   subroutine setup_sfc_pd_interaction(Natom, conf_num, NT, NA, nHam, &
      anumb, atype, max_no_neigh, max_no_equiv, max_no_shells, nlistsize, nn, &
      nlist, ncoup, coord, xc, fs_nlist, fs_nlistsize, do_ralloy, Natom_full, &
      Nchmax, atype_ch, asite_ch, achem_ch, ammom_inp, do_sortcoup, &
      do_lsf, nind, lsf_field, map_multiple)
      
      implicit none
      
      ! Input parameters
      integer, intent(in) :: NT, NA, nHam, Natom, Nchmax, conf_num
      integer, intent(in) :: do_ralloy, Natom_full, max_no_neigh, max_no_equiv, max_no_shells
      character, intent(in) :: do_sortcoup, do_lsf, lsf_field
      logical, intent(in) :: map_multiple
      integer, dimension(NT), intent(in) :: nn
      integer, dimension(Natom), intent(in) :: anumb, atype
      integer, dimension(Natom_full), intent(in) :: atype_ch, asite_ch, achem_ch
      real(dblprec), dimension(3,Natom), intent(in) :: coord
      real(dblprec), dimension(NA,Nchmax,conf_num), intent(in) :: ammom_inp
      real(dblprec), dimension(9,NT,max_no_shells,Nchmax,max(Nchmax,NT),conf_num), intent(in) :: xc
      
      ! Output parameters
      integer, dimension(Natom), intent(out) :: nlistsize
      integer, dimension(max_no_neigh,Natom), intent(out) :: nlist
      real(dblprec), dimension(9,max_no_neigh,nHam,conf_num), intent(out) :: ncoup
      integer, dimension(:), intent(out) :: fs_nlistsize
      integer, dimension(:,:), intent(out) :: fs_nlist
      integer, allocatable, dimension(:,:), intent(inout) :: nind
      
      write(*,'(a)') 'Setting up Pseudo-Dipolar interactions with SFC neighbor finding...'
      
      call setup_sfc_for_all_interactions('PD', Natom, conf_num, NT, NA, nHam, &
         anumb, atype, max_no_neigh, max_no_equiv, max_no_shells, nlistsize, nn, &
         nlist, ncoup, coord, xc, fs_nlist, fs_nlistsize, do_ralloy, Natom_full, &
         Nchmax, atype_ch, asite_ch, achem_ch, ammom_inp, 9, 1, do_sortcoup, &
         do_lsf, nind, lsf_field, map_multiple)
      
   end subroutine setup_sfc_pd_interaction

   !----------------------------------------------------------------------------
   ! SUBROUTINE: setup_sfc_biqdm_interaction
   !> @brief Setup BIQDM interactions using SFC (hdim=1)
   !----------------------------------------------------------------------------
   subroutine setup_sfc_biqdm_interaction(Natom, conf_num, NT, NA, nHam, &
      anumb, atype, max_no_neigh, max_no_equiv, max_no_shells, nlistsize, nn, &
      nlist, ncoup, coord, xc, fs_nlist, fs_nlistsize, do_ralloy, Natom_full, &
      Nchmax, atype_ch, asite_ch, achem_ch, ammom_inp, do_sortcoup, &
      do_lsf, nind, lsf_field, map_multiple)
      
      implicit none
      
      ! Input parameters
      integer, intent(in) :: NT, NA, nHam, Natom, Nchmax, conf_num
      integer, intent(in) :: do_ralloy, Natom_full, max_no_neigh, max_no_equiv, max_no_shells
      character, intent(in) :: do_sortcoup, do_lsf, lsf_field
      logical, intent(in) :: map_multiple
      integer, dimension(NT), intent(in) :: nn
      integer, dimension(Natom), intent(in) :: anumb, atype
      integer, dimension(Natom_full), intent(in) :: atype_ch, asite_ch, achem_ch
      real(dblprec), dimension(3,Natom), intent(in) :: coord
      real(dblprec), dimension(NA,Nchmax,conf_num), intent(in) :: ammom_inp
      real(dblprec), dimension(1,NT,max_no_shells,Nchmax,max(Nchmax,NT),conf_num), intent(in) :: xc
      
      ! Output parameters
      integer, dimension(Natom), intent(out) :: nlistsize
      integer, dimension(max_no_neigh,Natom), intent(out) :: nlist
      real(dblprec), dimension(1,max_no_neigh,nHam,conf_num), intent(out) :: ncoup
      integer, dimension(:), intent(out) :: fs_nlistsize
      integer, dimension(:,:), intent(out) :: fs_nlist
      integer, allocatable, dimension(:,:), intent(inout) :: nind
      
      write(*,'(a)') 'Setting up BIQDM interactions with SFC neighbor finding...'
      
      call setup_sfc_for_all_interactions('BIQDM', Natom, conf_num, NT, NA, nHam, &
         anumb, atype, max_no_neigh, max_no_equiv, max_no_shells, nlistsize, nn, &
         nlist, ncoup, coord, xc, fs_nlist, fs_nlistsize, do_ralloy, Natom_full, &
         Nchmax, atype_ch, asite_ch, achem_ch, ammom_inp, 1, 2, do_sortcoup, &
         do_lsf, nind, lsf_field, map_multiple)
      
   end subroutine setup_sfc_biqdm_interaction

   !----------------------------------------------------------------------------
   ! SUBROUTINE: setup_sfc_bq_interaction
   !> @brief Setup Biquadratic interactions using SFC (hdim=1)
   !----------------------------------------------------------------------------
   subroutine setup_sfc_bq_interaction(Natom, conf_num, NT, NA, nHam, &
      anumb, atype, max_no_neigh, max_no_equiv, max_no_shells, nlistsize, nn, &
      nlist, ncoup, coord, xc, fs_nlist, fs_nlistsize, do_ralloy, Natom_full, &
      Nchmax, atype_ch, asite_ch, achem_ch, ammom_inp, do_sortcoup, &
      do_lsf, nind, lsf_field, map_multiple)
      
      implicit none
      
      ! Input parameters
      integer, intent(in) :: NT, NA, nHam, Natom, Nchmax, conf_num
      integer, intent(in) :: do_ralloy, Natom_full, max_no_neigh, max_no_equiv, max_no_shells
      character, intent(in) :: do_sortcoup, do_lsf, lsf_field
      logical, intent(in) :: map_multiple
      integer, dimension(NT), intent(in) :: nn
      integer, dimension(Natom), intent(in) :: anumb, atype
      integer, dimension(Natom_full), intent(in) :: atype_ch, asite_ch, achem_ch
      real(dblprec), dimension(3,Natom), intent(in) :: coord
      real(dblprec), dimension(NA,Nchmax,conf_num), intent(in) :: ammom_inp
      real(dblprec), dimension(1,NT,max_no_shells,Nchmax,max(Nchmax,NT),conf_num), intent(in) :: xc
      
      ! Output parameters
      integer, dimension(Natom), intent(out) :: nlistsize
      integer, dimension(max_no_neigh,Natom), intent(out) :: nlist
      real(dblprec), dimension(1,max_no_neigh,nHam,conf_num), intent(out) :: ncoup
      integer, dimension(:), intent(out) :: fs_nlistsize
      integer, dimension(:,:), intent(out) :: fs_nlist
      integer, allocatable, dimension(:,:), intent(inout) :: nind
      
      write(*,'(a)') 'Setting up Biquadratic interactions with SFC neighbor finding...'
      
      call setup_sfc_for_all_interactions('BQ', Natom, conf_num, NT, NA, nHam, &
         anumb, atype, max_no_neigh, max_no_equiv, max_no_shells, nlistsize, nn, &
         nlist, ncoup, coord, xc, fs_nlist, fs_nlistsize, do_ralloy, Natom_full, &
         Nchmax, atype_ch, asite_ch, achem_ch, ammom_inp, 1, 2, do_sortcoup, &
         do_lsf, nind, lsf_field, map_multiple)
      
   end subroutine setup_sfc_bq_interaction

   !----------------------------------------------------------------------------
   ! SUBROUTINE: setup_sfc_ring_exchange
   !> @brief Setup Ring exchange interactions using SFC (hdim=1)
   !----------------------------------------------------------------------------
   subroutine setup_sfc_ring_exchange(Natom, conf_num, NT, NA, nHam, &
      anumb, atype, max_no_neigh, max_no_equiv, max_no_shells, nlistsize, nn, &
      nlist, ncoup, coord, xc, fs_nlist, fs_nlistsize, do_ralloy, Natom_full, &
      Nchmax, atype_ch, asite_ch, achem_ch, ammom_inp, do_sortcoup, &
      do_lsf, nind, lsf_field, map_multiple)
      
      implicit none
      
      ! Input parameters
      integer, intent(in) :: NT, NA, nHam, Natom, Nchmax, conf_num
      integer, intent(in) :: do_ralloy, Natom_full, max_no_neigh, max_no_equiv, max_no_shells
      character, intent(in) :: do_sortcoup, do_lsf, lsf_field
      logical, intent(in) :: map_multiple
      integer, dimension(NT), intent(in) :: nn
      integer, dimension(Natom), intent(in) :: anumb, atype
      integer, dimension(Natom_full), intent(in) :: atype_ch, asite_ch, achem_ch
      real(dblprec), dimension(3,Natom), intent(in) :: coord
      real(dblprec), dimension(NA,Nchmax,conf_num), intent(in) :: ammom_inp
      real(dblprec), dimension(1,NT,max_no_shells,Nchmax,max(Nchmax,NT),conf_num), intent(in) :: xc
      
      ! Output parameters
      integer, dimension(Natom), intent(out) :: nlistsize
      integer, dimension(max_no_neigh,Natom), intent(out) :: nlist
      real(dblprec), dimension(1,max_no_neigh,nHam,conf_num), intent(out) :: ncoup
      integer, dimension(:), intent(out) :: fs_nlistsize
      integer, dimension(:,:), intent(out) :: fs_nlist
      integer, allocatable, dimension(:,:), intent(inout) :: nind
      
      write(*,'(a)') 'Setting up Ring exchange interactions with SFC neighbor finding...'
      
      call setup_sfc_for_all_interactions('RingExchange', Natom, conf_num, NT, NA, nHam, &
         anumb, atype, max_no_neigh, max_no_equiv, max_no_shells, nlistsize, nn, &
         nlist, ncoup, coord, xc, fs_nlist, fs_nlistsize, do_ralloy, Natom_full, &
         Nchmax, atype_ch, asite_ch, achem_ch, ammom_inp, 1, 1, do_sortcoup, &
         do_lsf, nind, lsf_field, map_multiple)
      
   end subroutine setup_sfc_ring_exchange

   !----------------------------------------------------------------------------
   ! SUBROUTINE: setup_sfc_all_interactions
   !> @brief Master routine to setup all interaction types with proper dispatch
   !> @details This routine can be called with any interaction type and will
   !> automatically dispatch to the correct specific setup routine
   !----------------------------------------------------------------------------
   subroutine setup_sfc_all_interactions(interaction_type, Natom, conf_num, NT, NA, nHam, &
      anumb, atype, max_no_neigh, max_no_equiv, max_no_shells, nlistsize, nn, &
      nlist, ncoup, coord, xc, fs_nlist, fs_nlistsize, do_ralloy, Natom_full, &
      Nchmax, atype_ch, asite_ch, achem_ch, ammom_inp, hdim, lexp, do_sortcoup, &
      do_lsf, nind, lsf_field, map_multiple)
      
      implicit none
      
      ! Input parameters
      character(len=*), intent(in) :: interaction_type
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
      
      write(*,'(a,1x,a,1x,a,i0,a,i0,a)') 'SFC setup for', trim(interaction_type), &
         'interactions (hdim=', hdim, ', lexp=', lexp, ')'
      
      ! Direct call to the universal SFC integration
      call setup_sfc_for_all_interactions(interaction_type, Natom, conf_num, NT, NA, nHam, &
         anumb, atype, max_no_neigh, max_no_equiv, max_no_shells, nlistsize, nn, &
         nlist, ncoup, coord, xc, fs_nlist, fs_nlistsize, do_ralloy, Natom_full, &
         Nchmax, atype_ch, asite_ch, achem_ch, ammom_inp, hdim, lexp, do_sortcoup, &
         do_lsf, nind, lsf_field, map_multiple)
      
   end subroutine setup_sfc_all_interactions

end module SFCHamiltonianInterface