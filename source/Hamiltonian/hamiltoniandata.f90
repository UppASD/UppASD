!-------------------------------------------------------------------------------
! MODULE: HamiltonianData
!> @brief Data and allocation routines for the Hamiltonian
!> @copyright
!> GNU Public License.
!-------------------------------------------------------------------------------
module HamiltonianData
   use Profiling
   use Parameters
   !
   implicit none
   !
   ! From setup
   type ham_t
      integer, dimension(:), allocatable :: aHam !< Lookup-table for Hamiltonian
      ! Variables for Heisenberg exchange
      integer ::  max_no_neigh                                 !< Calculated maximum of neighbours for exchange
      integer, dimension(:), allocatable :: nlistsize          !< Size of neighbour list for Heisenberg exchange couplings
      integer, dimension(:,:), allocatable :: nlist            !< Neighbour list for Heisenberg exchange couplings
      real(dblprec), dimension(:,:,:), allocatable :: ncoup    !< Heisenberg exchange couplings
      real(dblprec), dimension(:,:,:), allocatable :: ncoupD   !< Heisenberg exchange couplings (DLM)
      real(dblprec), dimension(:,:,:,:), allocatable :: j_tens !< Exchange tensor (SKKR style)
      !! Reduced variables for Heisenberg exchange
      !integer :: NA_red                                          !< Number of reduced atoms
      !integer, dimension(:), allocatable :: nlistsize_red        !< Size of neighbour list for Heisenberg exchange couplings (reduced)
      !integer, dimension(:,:), allocatable :: nlist_red          !< Neighbour list for Heisenberg exchange couplings (reduced)
      !real(dblprec), dimension(:,:,:), allocatable :: ncoup_red  !< Heisenberg exchange couplings (reduced)
      !real(dblprec), dimension(:,:,:), allocatable :: ncoupD_red !< Heisenberg exchange couplings (DLM) (reduced)
      ! Variables for DMI
      integer ::  max_no_dmneigh                               !< Calculated maximum of neighbours for DM exchange
      integer, dimension(:), allocatable :: dmlistsize         !< Size of neighbour list for DM
      integer, dimension(:,:), allocatable :: dmlist           !< List of neighbours for DM
      real(dblprec), dimension(:,:,:), allocatable :: dm_vect  !< Dzyaloshinskii-Moriya exchange vector
      ! Variables for CHIR exchange
      integer :: nn_chir_tot                                     !< Calculated number of neighbours with CHIR interactions
      integer ::  max_no_chirneigh                               !< Calculated maximum of neighbours for CHIR exchange
      integer, dimension(:), allocatable :: chirlistsize         !< Size of neighbour list for CHIR
      integer, dimension(:,:,:), allocatable :: chirlist           !< List of neighbours for CHIR
      real(dblprec), dimension(:,:), allocatable :: chir_coup  !< scalar chirality exchange coupling
      ! Variables for PD exchange
      integer :: nn_pd_tot                                     !< Calculated number of neighbours with PD interactions
      integer ::  max_no_pdneigh                               !< Calculated maximum of neighbours for PD exchange
      integer, dimension(:), allocatable :: pdlistsize         !< Size of neighbour list for PD
      integer, dimension(:,:), allocatable :: pdlist           !< List of neighbours for PD
      real(dblprec), dimension(:,:,:), allocatable :: pd_vect  !< Pseudo-Dipolar exchange vector
      ! Variables for BIQDM interactions
      integer :: nn_biqdm_tot                                     !< Calculated number of neighbours with BIQDM interactions
      integer ::  max_no_biqdmneigh                               !< Calculated maximum of neighbours for BIQDM exchange
      integer, dimension(:), allocatable :: biqdmlistsize         !< Size of neighbour list for BIQDM
      integer, dimension(:,:), allocatable :: biqdmlist           !< List of neighbours for BIQDM
      real(dblprec), dimension(:,:,:), allocatable :: biqdm_vect  !< BIQDM exchange vector
      ! Variables for BQ interactions
      integer :: nn_bq_tot                               !< Calculated number of neighbours with BQ interactions
      integer, dimension(:), allocatable :: bqlistsize   !< Size of neighbour list for BQ
      integer, dimension(:,:), allocatable :: bqlist     !< List of neighbours for BQ
      real(dblprec), dimension(:,:), allocatable :: j_bq !< Biquadratic exchange couplings
      ! Variables for anisotropy
      integer, dimension(:), allocatable :: taniso                !< Type of anisotropy (0-2)
      integer, dimension(:), allocatable :: taniso_diff           !< Type of anisotropy (0-2)
      real(dblprec), dimension(:), allocatable :: sb              !< Ratio between Cubic and Uniaxial anisotropy
      real(dblprec), dimension(:), allocatable :: sb_diff         !< Ratio between Cubic and Uniaxial anisotropy
      real(dblprec), dimension(:,:), allocatable :: kaniso        !< Anisotropy constant
      real(dblprec), dimension(:,:), allocatable :: kaniso_diff   !< Anisotropy constant
      real(dblprec), dimension(:,:), allocatable :: eaniso        !< Unit anisotropy vector
      real(dblprec), dimension(:,:), allocatable :: eaniso_diff   !< Unit anisotropy vector
      ! Variables for induced moments
      integer :: fix_num                                    !< Number of "fixed" moments
      integer :: max_no_neigh_ind                           !< Number of nearest neighbours for the induced treatment
      integer, dimension(:), allocatable :: fix_list        !< List containing the "fixed" moments
      integer, dimension(:), allocatable :: ind_list_full   !< Indication of whether a given moment is induced/fixed 1/0
      integer, dimension(:), allocatable :: ind_nlistsize   !< Size of the list for the induced moments
      integer, dimension(:,:), allocatable :: ind_nlist     !< Neighbour list between induced moments and their first permanent moments
      integer, dimension(:), allocatable :: fix_nlistsize   !< Size of the list for the permanent moments
      integer, dimension(:,:), allocatable :: fix_nlist     !< Neighbour list between permanent moments and their first induced moments
      real(dblprec), dimension(:), allocatable :: sus_ind   !< Scaling factor for the magneitc moment of the induced moments
      ! Variables for LSF
      integer, dimension(:), allocatable :: fs_nlistsize    !< Size of first shell neighbouring list for centered atom
      integer, dimension(:,:), allocatable :: nind          !< Index of firstshell-neighbour-list corresponds to neighbour-list
      integer, dimension(:,:), allocatable :: fs_nlist      !< First shell Neighbouring list for centered atom
      ! Variables for dipolar
      real(dblprec), dimension(:,:,:,:), allocatable :: Qdip         !< Matrix for dipole-dipole interaction
      real(dblprec), dimension(:,:,:,:), allocatable :: Qdip_macro   !< Matrix for macro spin dipole-dipole interaction
   end type ham_t

   type(ham_t) :: ham
   public

contains

   !> Allocate arrays for anisotropy
   subroutine allocate_anisotropies(Natom,mult_axis,flag)
      implicit none

      integer, intent(in),optional :: Natom !< Number of atoms in system
      integer, intent(in) :: flag  !< Allocate or deallocate (1/-1)
      character(len=1), intent(in) :: mult_axis !< Flag to treat more than one anisotropy axis at the same time
      integer :: i_all, i_stat

      ! Allocate arrays for anisotropy
      if(flag>0) then
         allocate(ham%taniso(Natom),stat=i_stat)
         call memocc(i_stat,product(shape(ham%taniso))*kind(ham%taniso),'taniso','allocate_anisotropies')
         allocate(ham%eaniso(3,Natom),stat=i_stat)
         call memocc(i_stat,product(shape(ham%eaniso))*kind(ham%eaniso),'eaniso','allocate_anisotropies')
         allocate(ham%kaniso(2,Natom),stat=i_stat)
         call memocc(i_stat,product(shape(ham%kaniso))*kind(ham%kaniso),'kaniso','allocate_anisotropies')
         allocate(ham%sb(Natom),stat=i_stat)
         call memocc(i_stat,product(shape(ham%sb))*kind(ham%sb),'sb','allocate_anisotropies')

         if (mult_axis=='Y') then
            allocate(ham%taniso_diff(Natom),stat=i_stat)
            call memocc(i_stat,product(shape(ham%taniso_diff))*kind(ham%taniso_diff),'taniso_diff','allocate_anisotropies')
            allocate(ham%eaniso_diff(3,Natom),stat=i_stat)
            call memocc(i_stat,product(shape(ham%eaniso_diff))*kind(ham%eaniso_diff),'eaniso_diff','allocate_anisotropies')
            allocate(ham%kaniso_diff(2,Natom),stat=i_stat)
            call memocc(i_stat,product(shape(ham%kaniso_diff))*kind(ham%kaniso_diff),'kaniso_diff','allocate_anisotropies')
            allocate(ham%sb_diff(Natom),stat=i_stat)
            call memocc(i_stat,product(shape(ham%sb_diff))*kind(ham%sb_diff),'sb_diff','allocate_anisotropies')
         endif

      else
         i_all=-product(shape(ham%taniso))*kind(ham%taniso)
         deallocate(ham%taniso,stat=i_stat)
         call memocc(i_stat,i_all,'taniso','allocate_anisotropies')
         i_all=-product(shape(ham%eaniso))*kind(ham%eaniso)
         deallocate(ham%eaniso,stat=i_stat)
         call memocc(i_stat,i_all,'eaniso','allocate_anisotropies')
         i_all=-product(shape(ham%kaniso))*kind(ham%kaniso)
         deallocate(ham%kaniso,stat=i_stat)
         call memocc(i_stat,i_all,'kaniso','allocate_anisotropies')
         i_all=-product(shape(ham%sb))*kind(ham%sb)
         deallocate(ham%sb,stat=i_stat)
         call memocc(i_stat,i_all,'sb','allocate_anisotropies')

         if (mult_axis=='Y') then
            i_all=-product(shape(ham%taniso_diff))*kind(ham%taniso_diff)
            deallocate(ham%taniso_diff,stat=i_stat)
            call memocc(i_stat,i_all,'taniso_diff','allocate_anisotropies')
            i_all=-product(shape(ham%eaniso_diff))*kind(ham%eaniso_diff)
            deallocate(ham%eaniso_diff,stat=i_stat)
            call memocc(i_stat,i_all,'eaniso_diff','allocate_anisotropies')
            i_all=-product(shape(ham%kaniso_diff))*kind(ham%kaniso_diff)
            deallocate(ham%kaniso_diff,stat=i_stat)
            call memocc(i_stat,i_all,'kaniso_diff','allocate_anisotropies')
            i_all=-product(shape(ham%sb_diff))*kind(ham%sb_diff)
            deallocate(ham%sb_diff,stat=i_stat)
            call memocc(i_stat,i_all,'sb_diff','allocate_anisotropies')
         endif
      end if

   end subroutine allocate_anisotropies


   !> Allocate arrays for Heisenberg Hamiltonian
   subroutine allocate_hamiltoniandata(Natom,NA,nHam, conf_num,max_no_neigh,do_jtensor,do_lsf,&
         flag,lsf_field,exc_inter)
      implicit none

      integer, intent(in) ::  do_jtensor           !< Use SKKR style exchange tensor (0=off, 1=on)
      integer, intent(in),optional :: NA           !< Number of atoms in unit cell
      integer, intent(in),optional :: nHam         !< Number of atoms in Hamiltonian
      integer, intent(in),optional :: Natom        !< Number of atoms in system
      integer, intent(in),optional :: conf_num     !< Number of configurations for LSF
      integer, intent(in),optional :: max_no_neigh !< Calculated maximum of neighbours for exchange
      character(len=1), intent(in) :: do_lsf       !< Including LSF energy
      character(len=1), intent(in) :: lsf_field    !< LSF field term
      character(len=1), intent(in) :: exc_inter    !< Flag for interpolations of exchange
      integer, intent(in) :: flag  !< Allocate or deallocate (1/-1)

      integer :: i_all, i_stat
      ! Exchange
      if(flag>0) then
         allocate(ham%nlistsize(nHam),stat=i_stat)
         call memocc(i_stat,product(shape(ham%nlistsize))*kind(ham%nlistsize),'nlistsize','allocate_hamiltoniandata')
         ham%nlistsize=0
         !ham%NA_red=NA
         !allocate(ham%nlistsize_red(ham%NA_red),stat=i_stat)
         !call memocc(i_stat,product(shape(ham%nlistsize_red))*kind(ham%nlistsize_red),'nlistsize_red','allocate_hamiltoniandata')
         allocate(ham%nlist(max_no_neigh,Natom),stat=i_stat)
         call memocc(i_stat,product(shape(ham%nlist))*kind(ham%nlist),'nlist','allocate_hamiltoniandata')
         ham%nlist=0
         allocate(ham%aHam(Natom),stat=i_stat)
         call memocc(i_stat,product(shape(ham%aHam))*kind(ham%aHam),'aHam','allocate_hamiltoniandata')
         ham%aHam=0
         if(do_lsf=='Y' .and. lsf_field=='L') then
            allocate(ham%fs_nlistsize(Natom),stat=i_stat)
            call memocc(i_stat,product(shape(ham%fs_nlistsize))*kind(ham%fs_nlistsize),'fs_nlistsize','allocate_hamiltoniandata')
            ham%fs_nlistsize=0
            allocate(ham%fs_nlist(max_no_neigh,Natom),stat=i_stat)
            call memocc(i_stat,product(shape(ham%fs_nlist))*kind(ham%fs_nlist),'fs_nlist','allocate_hamiltoniandata')
            ham%fs_nlist=0
         endif

         if (do_jtensor/=1) then
            allocate(ham%ncoup(max_no_neigh,nHam,conf_num),stat=i_stat)
            call memocc(i_stat,product(shape(ham%ncoup))*kind(ham%ncoup),'ncoup','allocate_hamiltoniandata')
            ham%ncoup=0.0_dblprec
            !allocate(ham%ncoup_red(max_no_neigh,ham%NA_red,conf_num),stat=i_stat)
            !call memocc(i_stat,product(shape(ham%ncoup_red))*kind(ham%ncoup_red),'ncoup_red','allocate_hamiltoniandata')
            if (exc_inter=='Y') then
               allocate(ham%ncoupD(max_no_neigh,nHam,conf_num),stat=i_stat)
               call memocc(i_stat,product(shape(ham%ncoupD))*kind(ham%ncoupD),'ncoupD','allocate_hamiltoniandata')
               ham%ncoupD=0.0_dblprec
               !allocate(ham%ncoupD_red(max_no_neigh,NA,conf_num),stat=i_stat)
               !call memocc(i_stat,product(shape(ham%ncoupD_red))*kind(ham%ncoupD_red),'ncoupD_red','allocate_hamiltoniandata')
            endif
         else
            allocate(ham%j_tens(3,3,max_no_neigh,nHam),stat=i_stat)
            call memocc(i_stat,product(shape(ham%j_tens))*kind(ham%j_tens),'j_tens','allocate_hamiltoniandata')
            ham%j_tens=0.0_dblprec
         end if
      else
         i_all=-product(shape(ham%nlistsize))*kind(ham%nlistsize)
         deallocate(ham%nlistsize,stat=i_stat)
         call memocc(i_stat,i_all,'nlistsize','allocate_hamiltoniandata')
         i_all=-product(shape(ham%nlist))*kind(ham%nlist)
         deallocate(ham%nlist,stat=i_stat)
         call memocc(i_stat,i_all,'nlist','allocate_hamiltoniandata')
         i_all=-product(shape(ham%aHam))*kind(ham%aHam)
         deallocate(ham%aHam,stat=i_stat)
         call memocc(i_stat,i_all,'aHam','allocate_hamiltoniandata')
         if(do_lsf=='Y' .and. lsf_field=='L') then
            i_all=-product(shape(ham%fs_nlistsize))*kind(ham%fs_nlistsize)
            deallocate(ham%fs_nlistsize,stat=i_stat)
            call memocc(i_stat,i_all,'fs_nlistsize','allocate_hamiltoniandata')
            i_all=-product(shape(ham%fs_nlist))*kind(ham%fs_nlist)
            deallocate(ham%fs_nlist,stat=i_stat)
            call memocc(i_stat,i_all,'fs_nlist','allocate_hamiltoniandata')
            i_all=-product(shape(ham%nind))*kind(ham%nind)
            deallocate(ham%nind,stat=i_stat)
            call memocc(i_stat,i_all,'nind','allocate_hamiltoniandata')
         end if
         if (do_jtensor/=1) then
            i_all=-product(shape(ham%ncoup))*kind(ham%ncoup)
            deallocate(ham%ncoup,stat=i_stat)
            call memocc(i_stat,i_all,'ncoup','allocate_hamiltoniandata')
            if (exc_inter=='Y') then
               i_all=-product(shape(ham%ncoupD))*kind(ham%ncoupD)
               deallocate(ham%ncoupD,stat=i_stat)
               call memocc(i_stat,i_all,'ncoupD','allocate_hamiltoniandata')
            endif
         else
            i_all=-product(shape(ham%j_tens))*kind(ham%j_tens)
            deallocate(ham%j_tens,stat=i_stat)
            call memocc(i_stat,i_all,'j_tens','allocate_hamiltoniandata')
         end if
      end if

   end subroutine allocate_hamiltoniandata

   subroutine allocate_hamiltoniandata_ind(flag,Natom,max_no_neigh_ind)

      implicit none

      integer, intent(in) :: flag
      integer, intent(in), optional :: Natom
      integer, intent(in), optional :: max_no_neigh_ind

      integer :: i_stat, i_all

      if (flag>0) then
         allocate(ham%ind_nlistsize(Natom),stat=i_stat)
         call memocc(i_stat,product(shape(ham%ind_nlistsize))*kind(ham%ind_nlistsize),'ind_nlistsize','allocate_hamiltoniandata_ind')
         ham%ind_nlistsize=0
         allocate(ham%ind_nlist(max_no_neigh_ind,Natom),stat=i_stat)
         call memocc(i_stat,product(shape(ham%ind_nlist))*kind(ham%ind_nlist),'ind_nlist','allocate_hamiltoniandata_ind')
         ham%ind_nlist=0
         allocate(ham%fix_nlistsize(Natom),stat=i_stat)
         call memocc(i_stat,product(shape(ham%fix_nlistsize))*kind(ham%fix_nlistsize),'fix_nlistsize','allocate_hamiltoniandata_ind')
         ham%fix_nlistsize=0
         allocate(ham%fix_nlist(max_no_neigh_ind,Natom),stat=i_stat)
         call memocc(i_stat,product(shape(ham%fix_nlist))*kind(ham%fix_nlist),'fix_nlist','allocate_hamiltoniandata_ind')
         ham%fix_nlist=0
         allocate(ham%sus_ind(Natom),stat=i_stat)
         call memocc(i_stat,product(shape(ham%sus_ind))*kind(ham%sus_ind), 'sus_ind','allocate_hamiltoniandata_ind')
         ham%sus_ind=1.0_dblprec
      else
         i_all=-product(shape(ham%ind_nlistsize))*kind(ham%ind_nlistsize)
         deallocate(ham%ind_nlistsize,stat=i_stat)
         call memocc(i_stat,i_all,'ind_nlistsize','allocate_hamiltoniandata_ind')
         i_all=-product(shape(ham%ind_nlist))*kind(ham%ind_nlist)
         deallocate(ham%ind_nlist,stat=i_stat)
         call memocc(i_stat,i_all,'ind_nlist','allocate_hamiltoniandata_ind')
         i_all=-product(shape(ham%fix_nlistsize))*kind(ham%fix_nlistsize)
         deallocate(ham%fix_nlistsize,stat=i_stat)
         call memocc(i_stat,i_all,'fix_nlistsize','allocate_hamiltoniandata_ind')
         i_all=-product(shape(ham%fix_nlist))*kind(ham%fix_nlist)
         deallocate(ham%fix_nlist,stat=i_stat)
         call memocc(i_stat,i_all,'fix_nlist','allocate_hamiltoniandata_ind')
         i_all=-product(shape(ham%sus_ind))*kind(ham%sus_ind)
         deallocate(ham%sus_ind,stat=i_stat)
         call memocc(i_stat,i_all,'sus_ind','allocate_hamiltoniandata_ind')
         i_all=-product(shape(ham%ind_list_full))*kind(ham%ind_list_full)
         deallocate(ham%ind_list_full,stat=i_stat)
         call memocc(i_stat,i_all,'ind_list_full','allocate_hamiltoniandata_ind')
         if (allocated(ham%fix_list)) then
            i_all=-product(shape(ham%fix_list))*kind(ham%fix_list)
            deallocate(ham%fix_list,stat=i_stat)
            call memocc(i_stat,i_all,'fix_list','allocate_hamiltoniandata_ind')
         endif

      endif

   end subroutine allocate_hamiltoniandata_ind


   !> Allocate arrays for Dzyaloshinskii-Moriya Hamiltonian
   subroutine allocate_dmhamiltoniandata(Natom,nHam, max_no_dmneigh,flag)
      implicit none

      integer, optional, intent(in) :: Natom !< Number of atoms in system
      integer, optional, intent(in) :: nHam !< Number of atoms in Hamiltonian
      integer, optional, intent(in) :: max_no_dmneigh !< Calculated number of neighbours with DM interactions
      integer, intent(in) :: flag !< Allocate or deallocate (1/-1)

      integer :: i_all, i_stat

      if(flag>0) then
         allocate(ham%dmlistsize(nHam),stat=i_stat)
         call memocc(i_stat,product(shape(ham%dmlistsize))*kind(ham%dmlistsize),'dmlistsize','allocate_dmhamiltoniandata')
         allocate(ham%dmlist(max_no_dmneigh,Natom),stat=i_stat)
         call memocc(i_stat,product(shape(ham%dmlist))*kind(ham%dmlist),'dmlist','allocate_dmhamiltoniandata')
         allocate(ham%dm_vect(3,max_no_dmneigh,nHam),stat=i_stat)
         call memocc(i_stat,product(shape(ham%dm_vect))*kind(ham%dm_vect),'dm_vect','allocate_dmhamiltoniandata')
      else
         i_all=-product(shape(ham%dmlistsize))*kind(ham%dmlistsize)
         deallocate(ham%dmlistsize,stat=i_stat)
         call memocc(i_stat,i_all,'dmlistsize','allocate_dmhamiltoniandata')
         i_all=-product(shape(ham%dmlist))*kind(ham%dmlist)
         deallocate(ham%dmlist,stat=i_stat)
         call memocc(i_stat,i_all,'dmlist','allocate_dmhamiltoniandata')
         i_all=-product(shape(ham%dm_vect))*kind(ham%dm_vect)
         deallocate(ham%dm_vect,stat=i_stat)
         call memocc(i_stat,i_all,'dm_vect','allocate_dmhamiltoniandata')
      end if

   end subroutine allocate_dmhamiltoniandata


   !> Allocate arrays for scalar chirality  Hamiltonian
   subroutine allocate_chirhamiltoniandata(Natom,nHam,nn_chir_tot,flag)
      implicit none

      integer, optional, intent(in) :: Natom !< Number of atoms in system
      integer, optional, intent(in) :: nHam  !< Number of atoms in Hamiltonian
      integer, optional, intent(in) :: nn_chir_tot !< Calculated number of neighbours with chir interactions
      integer, intent(in) :: flag !< Allocate or deallocate (1/-1)

      integer :: i_all, i_stat

      if(flag>0) then
         allocate(ham%chirlistsize(nHam),stat=i_stat)
         call memocc(i_stat,product(shape(ham%chirlistsize))*kind(ham%chirlistsize),'chirlistsize','allocate_chirhamiltoniandata')
         allocate(ham%chirlist(2,nn_chir_tot,Natom),stat=i_stat)
         call memocc(i_stat,product(shape(ham%chirlist))*kind(ham%chirlist),'chirlist','allocate_chirhamiltoniandata')
         allocate(ham%chir_coup(nn_chir_tot,nHam),stat=i_stat)
         call memocc(i_stat,product(shape(ham%chir_coup))*kind(ham%chir_coup),'chir_coup','allocate_chirhamiltoniandata')
      else
         i_all=-product(shape(ham%chirlistsize))*kind(ham%chirlistsize)
         deallocate(ham%chirlistsize,stat=i_stat)
         call memocc(i_stat,i_all,'chirlistsize','allocate_chirhamiltoniandata')
         i_all=-product(shape(ham%chirlist))*kind(ham%chirlist)
         deallocate(ham%chirlist,stat=i_stat)
         call memocc(i_stat,i_all,'chirlist','allocate_chirhamiltoniandata')
         i_all=-product(shape(ham%chir_coup))*kind(ham%chir_coup)
         deallocate(ham%chir_coup,stat=i_stat)
         call memocc(i_stat,i_all,'chir_coup','allocate_chirhamiltoniandata')
      end if

   end subroutine allocate_chirhamiltoniandata


   !> Allocate arrays for Pseudo-Dipolar Hamiltonian
   subroutine allocate_pdhamiltoniandata(Natom,nHam,nn_pd_tot,flag)
      implicit none

      integer, optional, intent(in) :: Natom !< Number of atoms in system
      integer, optional, intent(in) :: nHam  !< Number of atoms in Hamiltonian
      integer, optional, intent(in) :: nn_pd_tot !< Calculated number of neighbours with PD interactions
      integer, intent(in) :: flag !< Allocate or deallocate (1/-1)

      integer :: i_all, i_stat

      if(flag>0) then
         allocate(ham%pdlistsize(nHam),stat=i_stat)
         call memocc(i_stat,product(shape(ham%pdlistsize))*kind(ham%pdlistsize),'pdlistsize','allocate_pdhamiltoniandata')
         allocate(ham%pdlist(nn_pd_tot,Natom),stat=i_stat)
         call memocc(i_stat,product(shape(ham%pdlist))*kind(ham%pdlist),'pdlist','allocate_pdhamiltoniandata')
         allocate(ham%pd_vect(6,nn_pd_tot,nHam),stat=i_stat)
         call memocc(i_stat,product(shape(ham%pd_vect))*kind(ham%pd_vect),'pd_vect','allocate_pdhamiltoniandata')
      else
         i_all=-product(shape(ham%pdlistsize))*kind(ham%pdlistsize)
         deallocate(ham%pdlistsize,stat=i_stat)
         call memocc(i_stat,i_all,'pdlistsize','allocate_pdhamiltoniandata')
         i_all=-product(shape(ham%pdlist))*kind(ham%pdlist)
         deallocate(ham%pdlist,stat=i_stat)
         call memocc(i_stat,i_all,'pdlist','allocate_pdhamiltoniandata')
         i_all=-product(shape(ham%pd_vect))*kind(ham%pd_vect)
         deallocate(ham%pd_vect,stat=i_stat)
         call memocc(i_stat,i_all,'pd_vect','allocate_pdhamiltoniandata')
      end if

   end subroutine allocate_pdhamiltoniandata


   !> Allocate arrays for BIQDM Hamiltonian
   subroutine allocate_biqdmhamiltoniandata(Natom,nHam, nn_biqdm_tot,flag)
      implicit none

      integer, optional, intent(in) :: Natom !< Number of atoms in system
      integer, optional, intent(in) :: nHam  !< Number of atoms in Hamiltonian
      integer, optional, intent(in) :: nn_biqdm_tot !< Calculated number of neighbours with BIQDM interactions
      integer, intent(in) :: flag !< Allocate or deallocate (1/-1)

      integer :: i_all, i_stat

      if(flag>0) then
         allocate(ham%biqdmlistsize(nHam),stat=i_stat)
         call memocc(i_stat,product(shape(ham%biqdmlistsize))*kind(ham%biqdmlistsize),'biqdmlistsize','allocate_biqdmhamiltoniandata')
         allocate(ham%biqdmlist(nn_biqdm_tot,Natom),stat=i_stat)
         call memocc(i_stat,product(shape(ham%biqdmlist))*kind(ham%biqdmlist),'biqdmlist','allocate_biqdmhamiltoniandata')
         allocate(ham%biqdm_vect(1,nn_biqdm_tot,nHam),stat=i_stat)
         call memocc(i_stat,product(shape(ham%biqdm_vect))*kind(ham%biqdm_vect),'biqdm_vect','allocate_biqdmhamiltoniandata')
      else
         i_all=-product(shape(ham%biqdmlistsize))*kind(ham%biqdmlistsize)
         deallocate(ham%biqdmlistsize,stat=i_stat)
         call memocc(i_stat,i_all,'biqdmlistsize','allocate_biqdmhamiltoniandata')
         i_all=-product(shape(ham%biqdmlist))*kind(ham%biqdmlist)
         deallocate(ham%biqdmlist,stat=i_stat)
         call memocc(i_stat,i_all,'biqdmlist','allocate_biqdmhamiltoniandata')
         i_all=-product(shape(ham%biqdm_vect))*kind(ham%biqdm_vect)
         deallocate(ham%biqdm_vect,stat=i_stat)
         call memocc(i_stat,i_all,'biqdm_vect','allocate_biqdmhamiltoniandata')
      end if

   end subroutine allocate_biqdmhamiltoniandata

   !----------------------------------------------------------------------------
   !> @brief Allocate arrays for biquadratic exchange Hamiltonian
   !----------------------------------------------------------------------------
   subroutine allocate_bqhamiltoniandata(Natom,nHam, nn_bq_tot,flag)
      implicit none

      integer, optional, intent(in) :: Natom !< Number of atoms in system
      integer, optional, intent(in) :: nHam  !< Number of atoms in Hamiltonian
      integer, optional, intent(in) :: nn_bq_tot !< Calculated number of neighbours with BQ interactions
      integer, intent(in) :: flag !< Allocate or deallocate (1/-1)

      integer :: i_all, i_stat

      ! Exchange
      if(flag>0) then
         allocate(ham%bqlistsize(Natom),stat=i_stat)
         call memocc(i_stat,product(shape(ham%bqlistsize))*kind(ham%bqlistsize),'bqlistsize','allocate_bqhamiltoniandata')
         allocate(ham%bqlist(nn_bq_tot,Natom),stat=i_stat)
         call memocc(i_stat,product(shape(ham%bqlist))*kind(ham%bqlist),'bqlist','allocate_bqhamiltoniandata')
         allocate(ham%j_bq(nn_bq_tot,nHam),stat=i_stat)
         call memocc(i_stat,product(shape(ham%j_bq))*kind(ham%j_bq),'j_bq','allocate_bqhamiltoniandata')
      else
         i_all=-product(shape(ham%bqlistsize))*kind(ham%bqlistsize)
         deallocate(ham%bqlistsize,stat=i_stat)
         call memocc(i_stat,i_all,'bqlistsize','allocate_bqhamiltoniandata')
         i_all=-product(shape(ham%bqlist))*kind(ham%bqlist)
         deallocate(ham%bqlist,stat=i_stat)
         call memocc(i_stat,i_all,'bqlist','allocate_bqhamiltoniandata')
         i_all=-product(shape(ham%j_bq))*kind(ham%j_bq)
         deallocate(ham%j_bq,stat=i_stat)
         call memocc(i_stat,i_all,'j_bq','allocate_bqhamiltoniandata')
      end if

   end subroutine allocate_bqhamiltoniandata

   !----------------------------------------------------------------------------
   !> @brief Allocate arrays for dipole matrix
   !----------------------------------------------------------------------------
   subroutine allocate_dipole(Natom,flag)
      implicit none

      integer, optional, intent(in) :: Natom !< Number of atoms in system
      integer, intent(in) :: flag !< Allocate or deallocate (1/-1)

      integer :: i_all, i_stat

      !  Allocate Q matrix
      if(flag>0) then
         allocate(ham%Qdip(3,3,natom,natom),stat=i_stat)
         call memocc(i_stat,product(shape(ham%Qdip))*kind(ham%Qdip),'Qdip','allocate_dipole')
      else
         i_all=-product(shape(ham%Qdip))*kind(ham%Qdip)
         deallocate(ham%Qdip,stat=i_stat)
         call memocc(i_stat,i_all,'Qdip','allocate_dipole')
      end if

   end subroutine allocate_dipole
   !
   !----------------------------------------------------------------------------
   !> @brief Allocate arrays for macro spin dipole matrix
   !> @author Jonathan Chico
   !----------------------------------------------------------------------------
   subroutine allocate_macro_dipole(Num_macro,flag)

      implicit none

      integer, intent(in) :: Num_macro !< Number of macrocells in the system
      integer, intent(in) :: flag !< Allocate or deallocate (1/-1)

      integer :: i_all, i_stat

      if (flag>0) then
         allocate(ham%Qdip_macro(3,3,Num_macro,Num_macro),stat=i_stat)
         call memocc(i_stat,product(shape(ham%Qdip_macro))*kind(ham%Qdip_macro),'Qdip_macro','allocate_macro_dipole')
         ham%Qdip_macro=0.0_dblprec
      else
         i_all=-product(shape(ham%Qdip_macro))*kind(ham%Qdip_macro)
         deallocate(ham%Qdip_macro,stat=i_stat)
         call memocc(i_stat,i_all,'Qdip_macro','allocate_macro_dipole')
      endif

   end subroutine allocate_macro_dipole


   subroutine scalar_to_tensor(nHam,do_dm)
      !
      implicit none
      !
      !
      integer, intent(in) :: nHam  !< Number of atoms in Hamiltonian
      integer, intent(in) :: do_dm   !< Add Dzyaloshinskii-Moriya (DM) term to Hamiltonian (0/1)
      !
      integer :: i,j, dn, jn, k
      !
      do i=1,nHam
         ham%j_tens(:,:,:,i)=0.0_dblprec
         !Exchange term
         do j=1,ham%nlistsize(i)
            ham%j_tens(1,1,j,i) = ham%ncoup(j,i,1)
            ham%j_tens(2,2,j,i) = ham%ncoup(j,i,1)
            ham%j_tens(3,3,j,i) = ham%ncoup(j,i,1)
         end do

         ! Dzyaloshinskii-Moriya term
         if(do_dm==1) then
            do j=1,ham%dmlistsize(i)
               dn  = ham%dmlist(j,i)
               do k=1,ham%nlistsize(i)
                  jn  = ham%nlist(k,i)
                  if(jn==dn) then
                     ham%j_tens(2,3,j,i) = ham%dm_vect(1,j,i)
                     ham%j_tens(1,3,j,i) = ham%dm_vect(2,j,i)
                     ham%j_tens(1,2,j,i) = ham%dm_vect(3,j,i)
                  end if
               end do
            end do
         end if

      end do

   end subroutine scalar_to_tensor

end module HamiltonianData
