!-------------------------------------------------------------------------------
! MODULE: HamiltonianData
!> @brief Data and allocation routines for the Hamiltonian
!> @copyright
!! Copyright (C) 2008-2018 UppASD group
!! This file is distributed under the terms of the
!! GNU General Public License.
!! See http://www.gnu.org/copyleft/gpl.txt
!-------------------------------------------------------------------------------
module HamiltonianData
   use Parameters
   use Profiling
   !
   implicit none
   !
   ! From setup
   ! Variables for Heisenberg exchange
   integer ::  max_no_neigh      !< Calculated maximum of neighbours for exchange
   integer, dimension(:), allocatable :: nlistsize       !< Size of neighbour list for Heisenberg exchange couplings
   integer, dimension(:,:), allocatable :: nlist         !< Neighbour list for Heisenberg exchange couplings
   real(dblprec), dimension(:,:,:), allocatable :: ncoup   !< Heisenberg exchange couplings
   real(dblprec), dimension(:,:,:), allocatable :: ncoupD   !< Heisenberg exchange couplings (DLM)
   real(dblprec), dimension(:,:,:,:), allocatable :: j_tens   !< Exchange tensor (SKKR style)
   ! Variables for DMI
   integer ::  max_no_dmneigh    !< Calculated maximum of neighbours for DM exchange
   integer, dimension(:), allocatable :: dmlistsize      !< Size of neighbour list for DM
   integer, dimension(:,:), allocatable :: dmlist        !< List of neighbours for DM
   real(dblprec), dimension(:,:,:), allocatable :: dm_vect    !< Dzyaloshinskii-Moriya exchange vector
   ! Variables for PD exchange
   integer ::  max_no_pdneigh    !< Calculated maximum of neighbours for PD exchange
   integer, dimension(:), allocatable :: pdlistsize      !< Size of neighbour list for PD
   integer, dimension(:,:), allocatable :: pdlist        !< List of neighbours for PD
   real(dblprec), dimension(:,:,:), allocatable :: pd_vect    !< Pseudo-Dipolar exchange vector
   ! Variables for BIQDM interactions
   integer ::  max_no_biqdmneigh !< Calculated maximum of neighbours for BIQDM exchange
   integer, dimension(:), allocatable :: biqdmlistsize   !< Size of neighbour list for BIQDM
   integer, dimension(:,:), allocatable :: biqdmlist     !< List of neighbours for BIQDM
   real(dblprec), dimension(:,:,:), allocatable :: biqdm_vect !< BIQDM exchange vector
   ! Variables for BQ interactions
   integer, dimension(:), allocatable :: bqlistsize      !< Size of neighbour list for BQ
   integer, dimension(:,:), allocatable :: bqlist        !< List of neighbours for BQ
   real(dblprec), dimension(:,:), allocatable :: j_bq         !< Biquadratic exchange couplings
   ! Variables for anisotropy
   integer, dimension(:), allocatable :: taniso          !< Type of anisotropy (0-2)
   integer, dimension(:), allocatable :: taniso_diff     !< Type of anisotropy (0-2)
   real(dblprec), dimension(:), allocatable :: sb             !< Ratio between Cubic and Uniaxial anisotropy
   real(dblprec), dimension(:), allocatable :: sb_diff        !< Ratio between Cubic and Uniaxial anisotropy
   real(dblprec), dimension(:,:), allocatable :: kaniso       !< Anisotropy constant
   real(dblprec), dimension(:,:), allocatable :: kaniso_diff  !< Anisotropy constant
   real(dblprec), dimension(:,:), allocatable :: eaniso       !< Unit anisotropy vector
   real(dblprec), dimension(:,:), allocatable :: eaniso_diff  !< Unit anisotropy vector
   ! Variables for induced moments
   integer, dimension(:), allocatable :: ind_nlistsize !< Size of the list for the induced moments
   integer, dimension(:,:), allocatable :: ind_nlist   !< Neighbour list between induced moments and their first permanent moments
   real(dblprec), dimension(:), allocatable :: sus_ind !< Scaling factor for the magneitc moment of the induced moments
   ! Variables for LSF
   integer, dimension(:), allocatable :: fs_nlistsize    !< Size of first shell neighbouring list for centered atom
   integer, dimension(:,:), allocatable :: fs_nlist      !< First shell Neighbouring list for centered atom
   integer, dimension(:,:), allocatable :: nind          !< Index of firstshell-neighbour-list corresponds to neighbour-list
   ! Variables for dipolar
   real(dblprec), dimension(:,:,:,:), allocatable :: Qdip     !< Matrix for dipole-dipole interaction

   public

contains


   !> Allocate arrays for anisotropy
   subroutine allocate_anisotropies(Natom,mult_axis,flag)
      implicit none

      integer, intent(in),optional :: Natom !< Number of atoms in system
      integer, intent(in) :: flag  !< Allocate or deallocate (1/-1)
      character(len=1), intent(in) :: mult_axis
      integer :: i_all, i_stat

      ! Allocate arrays for anisotropy
      if(flag>0) then
         allocate(taniso(Natom),stat=i_stat)
         call memocc(i_stat,product(shape(taniso))*kind(taniso),'taniso','setup_anisotropies')
         allocate(eaniso(3,Natom),stat=i_stat)
         call memocc(i_stat,product(shape(eaniso))*kind(eaniso),'eaniso','setup_anisotropies')
         allocate(kaniso(2,Natom),stat=i_stat)
         call memocc(i_stat,product(shape(kaniso))*kind(kaniso),'kaniso','setup_anisotropies')
         allocate(sb(Natom),stat=i_stat)
         call memocc(i_stat,product(shape(sb))*kind(sb),'sb','setup_anisotropies')

         if (mult_axis=='Y') then
            allocate(taniso_diff(Natom),stat=i_stat)
            call memocc(i_stat,product(shape(taniso_diff))*kind(taniso_diff),'taniso_diff','setup_anisotropies')
            allocate(eaniso_diff(3,Natom),stat=i_stat)
            call memocc(i_stat,product(shape(eaniso_diff))*kind(eaniso_diff),'eaniso_diff','setup_anisotropies')
            allocate(kaniso_diff(2,Natom),stat=i_stat)
            call memocc(i_stat,product(shape(kaniso_diff))*kind(kaniso_diff),'kaniso_diff','setup_anisotropies')
            allocate(sb_diff(Natom),stat=i_stat)
            call memocc(i_stat,product(shape(sb_diff))*kind(sb_diff),'sb_diff','setup_anisotropies')
         endif

      else
         i_all=-product(shape(taniso))*kind(taniso)
         deallocate(taniso,stat=i_stat)
         call memocc(i_stat,i_all,'taniso','setup_anisotropies')
         i_all=-product(shape(eaniso))*kind(eaniso)
         deallocate(eaniso,stat=i_stat)
         call memocc(i_stat,i_all,'eaniso','setup_anisotropies')
         i_all=-product(shape(kaniso))*kind(kaniso)
         deallocate(kaniso,stat=i_stat)
         call memocc(i_stat,i_all,'kaniso','setup_anisotropies')
         i_all=-product(shape(sb))*kind(sb)
         deallocate(sb,stat=i_stat)
         call memocc(i_stat,i_all,'sb','setup_anisotropies')

         if (mult_axis=='Y') then
            i_all=-product(shape(taniso_diff))*kind(taniso_diff)
            deallocate(taniso_diff,stat=i_stat)
            call memocc(i_stat,i_all,'taniso_diff','setup_anisotropies')
            i_all=-product(shape(eaniso_diff))*kind(eaniso_diff)
            deallocate(eaniso_diff,stat=i_stat)
            call memocc(i_stat,i_all,'eaniso_diff','setup_anisotropies')
            i_all=-product(shape(kaniso_diff))*kind(kaniso_diff)
            deallocate(kaniso_diff,stat=i_stat)
            call memocc(i_stat,i_all,'kaniso_diff','setup_anisotropies')
            i_all=-product(shape(sb_diff))*kind(sb_diff)
            deallocate(sb_diff,stat=i_stat)
            call memocc(i_stat,i_all,'sb_diff','setup_anisotropies')
         endif
      end if

   end subroutine allocate_anisotropies


   !> Allocate arrays for Heisenberg Hamiltonian
   subroutine allocate_hamiltoniandata(Natom,conf_num,max_no_neigh,do_jtensor,do_lsf,&
         ind_mom_flag,flag,lsf_field,exc_inter)
      implicit none

      integer, intent(in) ::  do_jtensor           !<  Use SKKR style exchange tensor (0=off, 1=on)
      integer, intent(in),optional :: Natom        !< Number of atoms in system
      integer, intent(in),optional :: conf_num !< Number of configurations for LSF
      integer, intent(in),optional :: max_no_neigh !< Calculated maximum of neighbours for exchange
      character(len=1), intent(in) :: do_lsf       !< Including LSF energy
      character(len=1), intent(in) :: lsf_field       !< LSF field term
      character(len=1), intent(in) :: ind_mom_flag !< Flag for the induced moments
      character(len=1), intent(in) :: exc_inter !< Flag for interpolations of exchange
      integer, intent(in) :: flag  !< Allocate or deallocate (1/-1)

      integer :: i_all, i_stat
      ! Exchange
      if(flag>0) then
         allocate(nlistsize(Natom),stat=i_stat)
         call memocc(i_stat,product(shape(nlistsize))*kind(nlistsize),'nlistsize','allocate_hamiltoniandata')
         allocate(nlist(max_no_neigh,Natom),stat=i_stat)
         call memocc(i_stat,product(shape(nlist))*kind(nlist),'nlist','allocate_hamiltoniandata')

         if(do_lsf=='Y' .and. lsf_field=='L') then
            allocate(fs_nlistsize(Natom),stat=i_stat)
            call memocc(i_stat,product(shape(fs_nlistsize))*kind(fs_nlistsize),'fs_nlistsize','allocate_hamiltoniandata')
            allocate(fs_nlist(max_no_neigh,Natom),stat=i_stat)
            call memocc(i_stat,product(shape(fs_nlist))*kind(fs_nlist),'fs_nlist','allocate_hamiltoniandata')
         endif

         if (ind_mom_flag=='Y') then
            allocate(ind_nlistsize(Natom),stat=i_stat)
            call memocc(i_stat,product(shape(ind_nlistsize))*kind(ind_nlistsize),'ind_nlistsize','allocate_hamiltoniandata')
            allocate(ind_nlist(max_no_neigh,Natom),stat=i_stat)
            call memocc(i_stat,product(shape(ind_nlist))*kind(ind_nlist),'ind_nlist','allocate_hamiltoniandata')
            allocate(sus_ind(Natom),stat=i_stat)
            call memocc(i_stat,product(shape(sus_ind))*kind(sus_ind), 'sus_ind','allocate_hamiltoniandata')
         endif

         if (do_jtensor/=1) then
            allocate(ncoup(max_no_neigh,Natom,conf_num),stat=i_stat)
            call memocc(i_stat,product(shape(ncoup))*kind(ncoup),'ncoup','allocate_hamiltoniandata')
            if (exc_inter=='Y') then
               allocate(ncoupD(max_no_neigh,Natom,conf_num),stat=i_stat)
               call memocc(i_stat,product(shape(ncoupD))*kind(ncoupD),'ncoupD','allocate_hamiltoniandata')
            endif
         else
            allocate(j_tens(3,3,max_no_neigh,Natom),stat=i_stat)
            call memocc(i_stat,product(shape(j_tens))*kind(j_tens),'j_tens','allocate_hamiltoniandata')
         end if
      else
         i_all=-product(shape(nlistsize))*kind(nlistsize)
         deallocate(nlistsize,stat=i_stat)
         call memocc(i_stat,i_all,'nlistsize','allocate_hamiltoniandata')
         i_all=-product(shape(nlist))*kind(nlist)
         deallocate(nlist,stat=i_stat)
         call memocc(i_stat,i_all,'nlist','allocate_hamiltoniandata')

         if(do_lsf=='Y' .and. lsf_field=='L') then
            i_all=-product(shape(fs_nlistsize))*kind(fs_nlistsize)
            deallocate(fs_nlistsize,stat=i_stat)
            call memocc(i_stat,i_all,'fs_nlistsize','allocate_hamiltoniandata')
            i_all=-product(shape(fs_nlist))*kind(fs_nlist)
            deallocate(fs_nlist,stat=i_stat)
            call memocc(i_stat,i_all,'fs_nlist','allocate_hamiltoniandata')
            i_all=-product(shape(nind))*kind(nind)
            deallocate(nind,stat=i_stat)
            call memocc(i_stat,i_all,'nind','allocate_hamiltoniandata')
         end if

         if (ind_mom_flag=='Y') then
            i_all=-product(shape(ind_nlistsize))*kind(ind_nlistsize)
            deallocate(ind_nlistsize,stat=i_stat)
            call memocc(i_stat,i_all,'ind_nlistsize','allocate_hamiltoniandata')
            i_all=-product(shape(ind_nlist))*kind(ind_nlist)
            deallocate(ind_nlist,stat=i_stat)
            call memocc(i_stat,i_all,'ind_nlist','allocate_hamiltoniandata')
            i_all=-product(shape(sus_ind))*kind(sus_ind)
            deallocate(sus_ind,stat=i_stat)
            call memocc(i_stat,i_all,'sus_ind','allocate_hamiltoniandata')

         endif

         if (do_jtensor/=1) then
            i_all=-product(shape(ncoup))*kind(ncoup)
            deallocate(ncoup,stat=i_stat)
            call memocc(i_stat,i_all,'ncoup','allocate_hamiltoniandata')
            if (exc_inter=='Y') then
               i_all=-product(shape(ncoupD))*kind(ncoupD)
               deallocate(ncoupD,stat=i_stat)
               call memocc(i_stat,i_all,'ncoupD','allocate_hamiltoniandata')
            endif
         else
            i_all=-product(shape(j_tens))*kind(j_tens)
            deallocate(j_tens,stat=i_stat)
            call memocc(i_stat,i_all,'j_tens','allocate_hamiltoniandata')
         end if
      end if

   end subroutine allocate_hamiltoniandata


   !> Allocate arrays for Dzyaloshinskii-Moriya Hamiltonian
   subroutine allocate_dmhamiltoniandata(Natom,max_no_dmneigh,flag)
      implicit none

      integer, optional, intent(in) :: Natom !< Number of atoms in system
      integer, optional, intent(in) :: max_no_dmneigh !< Calculated number of neighbours with DM interactions
      integer, intent(in) :: flag !< Allocate or deallocate (1/-1)

      integer :: i_all, i_stat

      if(flag>0) then
         allocate(dmlistsize(Natom),stat=i_stat)
         call memocc(i_stat,product(shape(dmlistsize))*kind(dmlistsize),'dmlistsize','allocate_dmhamiltoniandata')
         allocate(dmlist(max_no_dmneigh,Natom),stat=i_stat)
         call memocc(i_stat,product(shape(dmlist))*kind(dmlist),'dmlist','allocate_dmhamiltoniandata')
         allocate(dm_vect(3,max_no_dmneigh,Natom),stat=i_stat)
         call memocc(i_stat,product(shape(dm_vect))*kind(dm_vect),'dm_vect','allocate_dmhamiltoniandata')
      else
         i_all=-product(shape(dmlistsize))*kind(dmlistsize)
         deallocate(dmlistsize,stat=i_stat)
         call memocc(i_stat,i_all,'dmlistsize','allocate_dmhamiltoniandata')
         i_all=-product(shape(dmlist))*kind(dmlist)
         deallocate(dmlist,stat=i_stat)
         call memocc(i_stat,i_all,'dmlist','allocate_dmhamiltoniandata')
         i_all=-product(shape(dm_vect))*kind(dm_vect)
         deallocate(dm_vect,stat=i_stat)
         call memocc(i_stat,i_all,'dm_vect','allocate_dmhamiltoniandata')
      end if

   end subroutine allocate_dmhamiltoniandata


   !> Allocate arrays for Pseudo-Dipolar Hamiltonian
   subroutine allocate_pdhamiltoniandata(Natom,nn_pd_tot,flag)
      implicit none

      integer, optional, intent(in) :: Natom !< Number of atoms in system
      integer, optional, intent(in) :: nn_pd_tot !< Calculated number of neighbours with PD interactions
      integer, intent(in) :: flag !< Allocate or deallocate (1/-1)

      integer :: i_all, i_stat

      if(flag>0) then
         allocate(pdlistsize(Natom),stat=i_stat)
         call memocc(i_stat,product(shape(pdlistsize))*kind(pdlistsize),'pdlistsize','allocate_pdhamiltoniandata')
         allocate(pdlist(nn_pd_tot,Natom),stat=i_stat)
         call memocc(i_stat,product(shape(pdlist))*kind(pdlist),'pdlist','allocate_pdhamiltoniandata')
         allocate(pd_vect(6,nn_pd_tot,Natom),stat=i_stat)
         call memocc(i_stat,product(shape(pd_vect))*kind(pd_vect),'pd_vect','allocate_pdhamiltoniandata')
      else
         i_all=-product(shape(pdlistsize))*kind(pdlistsize)
         deallocate(pdlistsize,stat=i_stat)
         call memocc(i_stat,i_all,'pdlistsize','allocate_pdhamiltoniandata')
         i_all=-product(shape(pdlist))*kind(pdlist)
         deallocate(pdlist,stat=i_stat)
         call memocc(i_stat,i_all,'pdlist','allocate_pdhamiltoniandata')
         i_all=-product(shape(pd_vect))*kind(pd_vect)
         deallocate(pd_vect,stat=i_stat)
         call memocc(i_stat,i_all,'pd_vect','allocate_pdhamiltoniandata')
      end if

   end subroutine allocate_pdhamiltoniandata


   !> Allocate arrays for BIQDM Hamiltonian
   subroutine allocate_biqdmhamiltoniandata(Natom,nn_biqdm_tot,flag)
      implicit none

      integer, optional, intent(in) :: Natom !< Number of atoms in system
      integer, optional, intent(in) :: nn_biqdm_tot !< Calculated number of neighbours with BIQDM interactions
      integer, intent(in) :: flag !< Allocate or deallocate (1/-1)

      integer :: i_all, i_stat

      if(flag>0) then
         allocate(biqdmlistsize(Natom),stat=i_stat)
         call memocc(i_stat,product(shape(biqdmlistsize))*kind(biqdmlistsize),'biqdmlistsize','allocate_biqdmhamiltoniandata')
         allocate(biqdmlist(nn_biqdm_tot,Natom),stat=i_stat)
         call memocc(i_stat,product(shape(biqdmlist))*kind(biqdmlist),'biqdmlist','allocate_biqdmhamiltoniandata')
         allocate(biqdm_vect(1,nn_biqdm_tot,Natom),stat=i_stat)
         call memocc(i_stat,product(shape(biqdm_vect))*kind(biqdm_vect),'biqdm_vect','allocate_biqdmhamiltoniandata')
      else
         i_all=-product(shape(biqdmlistsize))*kind(biqdmlistsize)
         deallocate(biqdmlistsize,stat=i_stat)
         call memocc(i_stat,i_all,'biqdmlistsize','allocate_biqdmhamiltoniandata')
         i_all=-product(shape(biqdmlist))*kind(biqdmlist)
         deallocate(biqdmlist,stat=i_stat)
         call memocc(i_stat,i_all,'biqdmlist','allocate_biqdmhamiltoniandata')
         i_all=-product(shape(biqdm_vect))*kind(biqdm_vect)
         deallocate(biqdm_vect,stat=i_stat)
         call memocc(i_stat,i_all,'biqdm_vect','allocate_biqdmhamiltoniandata')
      end if

   end subroutine allocate_biqdmhamiltoniandata


   !> Allocate arrays for biquadratic exchange Hamiltonian
   subroutine allocate_bqhamiltoniandata(Natom,nn_bq_tot,flag)
      implicit none

      integer, optional, intent(in) :: Natom !< Number of atoms in system
      integer, optional, intent(in) :: nn_bq_tot !< Calculated number of neighbours with BQ interactions
      integer, intent(in) :: flag !< Allocate or deallocate (1/-1)

      integer :: i_all, i_stat

      ! Exchange
      if(flag>0) then
         allocate(bqlistsize(Natom),stat=i_stat)
         call memocc(i_stat,product(shape(bqlistsize))*kind(bqlistsize),'bqlistsize','allocate_bqhamiltoniandata')
         allocate(bqlist(nn_bq_tot,Natom),stat=i_stat)
         call memocc(i_stat,product(shape(bqlist))*kind(bqlist),'bqlist','allocate_bqhamiltoniandata')
         allocate(j_bq(nn_bq_tot,Natom),stat=i_stat)
         call memocc(i_stat,product(shape(j_bq))*kind(j_bq),'j_bq','allocate_bqhamiltoniandata')
      else
         i_all=-product(shape(bqlistsize))*kind(bqlistsize)
         deallocate(bqlistsize,stat=i_stat)
         call memocc(i_stat,i_all,'bqlistsize','allocate_bqhamiltoniandata')
         i_all=-product(shape(bqlist))*kind(bqlist)
         deallocate(bqlist,stat=i_stat)
         call memocc(i_stat,i_all,'bqlist','allocate_bqhamiltoniandata')
         i_all=-product(shape(j_bq))*kind(j_bq)
         deallocate(j_bq,stat=i_stat)
         call memocc(i_stat,i_all,'j_bq','allocate_bqhamiltoniandata')
      end if

   end subroutine allocate_bqhamiltoniandata


   !> Allocate arrays for dipole matrix
   subroutine allocate_dipole(Natom,flag)
      implicit none

      integer, intent(in) :: Natom !< Number of atoms in system
      integer, intent(in) :: flag !< Allocate or deallocate (1/-1)

      integer :: i_all, i_stat

      !  Allocate Q matrix
      if(flag>0) then
         allocate(Qdip(3,3,natom,natom),stat=i_stat)
         call memocc(i_stat,product(shape(Qdip))*kind(Qdip),'Qdip','allocate_dipole')
      else
         i_all=-product(shape(Qdip))*kind(Qdip)
         deallocate(Qdip,stat=i_stat)
         call memocc(i_stat,i_all,'Qdip','allocate_dipole')
      end if

   end subroutine allocate_dipole
   !

   subroutine scalar_to_tensor(Natom,do_dm)
      !
      implicit none
      !
      !
      integer, intent(in) :: Natom !< Number of atoms in system
      integer, intent(in) :: do_dm   !< Add Dzyaloshinskii-Moriya (DM) term to Hamiltonian (0/1)
      !
      integer :: i,j, dn, jn, k
      !
      do i=1,Natom
         j_tens(:,:,:,i)=0.0d0
         !Exchange term
         do j=1,nlistsize(i)
            j_tens(1,1,j,i) = ncoup(j,i,1)
            j_tens(2,2,j,i) = ncoup(j,i,1)
            j_tens(3,3,j,i) = ncoup(j,i,1)
         end do

         ! Dzyaloshinskii-Moriya term
         if(do_dm==1) then
            do j=1,dmlistsize(i)
               dn  = dmlist(j,i)
               do k=1,nlistsize(i)
                  jn  = nlist(k,i)
                  if(jn==dn) then
                     j_tens(2,3,j,i) = dm_vect(1,j,i)
                     j_tens(1,3,j,i) = dm_vect(2,j,i)
                     j_tens(1,2,j,i) = dm_vect(3,j,i)
                  end if
               end do
            end do
         end if

      end do

   end subroutine scalar_to_tensor

end module HamiltonianData
