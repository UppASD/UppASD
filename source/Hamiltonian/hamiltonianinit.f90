!-------------------------------------------------------------------------------
!  MODULE: HamiltonianInit
!> @brief
!> Routines for mounting the Hamiltonian, including exchange, DM, anisotropies
!> @author
!> Anders Bergman, Lars Bergqvist, Johan Hellsvik
!> @copyright
!> GNU Public License.
!-------------------------------------------------------------------------------
module HamiltonianInit

   use Profiling
   use Parameters

   implicit none

   integer, dimension(:,:,:,:), allocatable :: nme !< Neighbour map, elements
   integer, dimension(:,:,:), allocatable :: nm    !< Neighbour map
   integer, dimension(:,:), allocatable :: nmdim   !< Dimension of neighbour map

   private
   public :: setup_hamiltonian!, setup_reduced_hamiltonian

contains

   !----------------------------------------------------------------------------
   !  SUBROUTINE: setup_hamiltonian
   !> @brief
   !>  Main routine for mounting the Hamiltonian
   !> @todo Change mounting of DM Hamiltonian to follow the Heisenberg approach
   !> @todo Introduce additional flags for symmetry in order to enable use of symops for other couplings than Heisenberg
   !----------------------------------------------------------------------------
   subroutine setup_hamiltonian(NT,NA,N1,N2,N3,Nchmax,do_ralloy,Natom_full,         &
      Mensemble,nHam,Natom,achtype,atype_ch,asite_ch,achem_ch,acellnumb,            &
      acellnumbrev,atype,anumb,alat,C1,C2,C3,Bas,ammom_inp,coord,BC1,BC2,BC3,sym,   &
      do_jtensor,max_no_shells,nn,redcoord,jc,jcD,jc_tens,do_dm,max_no_dmshells,    &
      dm_nn,dm_redcoord,dm_inpvect,do_anisotropy,random_anisotropy_density,         &
      anisotropytype,anisotropytype_diff,anisotropy,anisotropy_diff,                &
      random_anisotropy,mult_axis,mconf,conf_num,map_multiple,do_lsf,lsf_field,     &
      exc_inter,do_bq,max_no_bqshells,bq_nn,bq_redcoord,jc_bq,do_biqdm,             &
      max_no_biqdmshells,biqdm_nn,biqdm_redcoord,biqdm_inpvect,do_pd,               &
      max_no_pdshells,pd_nn,pd_redcoord,pd_inpvect,do_dip,do_chir,max_no_chirshells,&
      chir_nn,chir_redcoord,chir_inpval,ind_mom,ind_tol,ind_mom_flag,NA_clus,       &
      NT_clus,N1_clus,N2_clus,N3_clus,Natom_clus,Nchmax_clus,clus_expand,           &
      Natom_full_clus,max_no_shells_clus,max_no_dmshells_clus,NN_clus,dm_nn_clus,   &
      achtype_clus,atype_ch_clus,asite_ch_clus,achem_ch_clus,acellnumb_clus,        &
      acellnumbrev_clus,atype_clus,anumb_clus,atype_inp_clus,anumb_inp_clus,C1_clus,&
      C2_clus,C3_clus,Bas_clus,ammom_inp_clus,redcoord_clus,jc_clus,                &
      dm_redcoord_clus,dm_inpvect_clus,do_cluster,do_anisotropy_clus,               &
      random_anisotropy_density_clus,anisotropytype_clus,anisotropytype_diff_clus,  &
      anisotropy_clus,anisotropy_diff_clus,random_anisotropy_clus,mult_axis_clus,   &
      Num_macro,max_num_atom_macro_cell,macro_nlistsize,macro_atom_nlist,block_size,&
      do_reduced,do_prnstruct,do_sortcoup,simid,print_dip_tensor,read_dipole,       &
      qdip_files)

      use LSF,             only : LSF_datareshape
      use clusters,        only : allocate_cluster_hamiltoniandata,                 &
         allocate_cluster_dmhamiltoniandata, allocate_cluster_anisotropies, ham_clus
      use InputData,       only : jij_scale, dm_scale
      use NeighbourMap,    only : setup_nm, setup_nm_nelem
      use InducedMoments
      use HamiltonianData, only : allocate_hamiltoniandata, allocate_anisotropies,  &
         allocate_dmhamiltoniandata, allocate_pdhamiltoniandata,                    &
         allocate_biqdmhamiltoniandata, allocate_bqhamiltoniandata,                 &
         allocate_chirhamiltoniandata, ham
      use PrintHamiltonian
      use DipoleManager, only : dipole_setup
      use ErrorHandling

      !.. Implicit declarations
      implicit none
      ! System/geometry variables
      integer, intent(in) :: NT         !< Number of types of atoms
      integer, intent(in) :: NA         !< Number of atoms in one cell
      integer, intent(in) :: N1         !< Number of cell repetitions in x direction
      integer, intent(in) :: N2         !< Number of cell repetitions in y direction
      integer, intent(in) :: N3         !< Number of cell repetitions in z direction
      integer, intent(in) :: Nchmax     !< Max number of chemical components on each site in cell
      integer, intent(in) :: do_ralloy  !< Random alloy simulation (0/1)
      integer, intent(in) :: Natom_full !< Number of atoms for full system (=Natom if not dilute)
      integer, intent(in) :: Mensemble
      integer, intent(inout) :: nHam    !< Number of atoms in Hamiltonian
      integer, intent(inout) :: Natom   !< Number of atoms in system
      integer, dimension(:), intent(inout) :: achtype       !< Chemical type of atoms (full list)
      integer, dimension(:), intent(inout) :: atype_ch      !< Actual type of atom for dilute system
      integer, dimension(:), intent(inout) :: asite_ch      !< Actual site of atom for dilute system
      integer, dimension(:), intent(inout) :: achem_ch      !< Chemical type of atoms (reduced list)
      integer, dimension(:), intent(inout) :: acellnumb     !< List for translating atom no. in full cell to actual cell
      integer, dimension(:), intent(inout) :: acellnumbrev  !< List for translating atom no. in actual cell to full cell
      integer, dimension(Natom), intent(inout) :: atype     !< Type of atom
      integer, dimension(Natom), intent(inout) :: anumb     !< Atom number in cell
      real(dblprec), intent(in) :: alat                     !< Lattice parameter
      real(dblprec), dimension(3), intent(in) :: C1         !< First lattice vector
      real(dblprec), dimension(3), intent(in) :: C2         !< Second lattice vector
      real(dblprec), dimension(3), intent(in) :: C3         !< Third lattice vector
      real(dblprec), dimension(3,NA), intent(inout) :: Bas !< Coordinates for basis atoms
      real(dblprec), dimension(:,:,:), intent(inout) :: ammom_inp !< Magnetic moment directions from input (for alloys)
      real(dblprec), dimension(:,:), allocatable, intent(inout) :: coord !< Coordinates of atoms
      character(len=1), intent(in) :: BC1  !< Boundary conditions in x-direction
      character(len=1), intent(in) :: BC2  !< Boundary conditions in y-direction
      character(len=1), intent(in) :: BC3  !< Boundary conditions in z-direction
      ! Heisenberg exchange variables
      integer, intent(in) :: sym              !< Symmetry of system (0-3)
      integer, intent(in) :: do_jtensor       !< Use SKKR style exchange tensor (0=off, 1=on)
      integer, intent(inout) :: max_no_shells !< Calculated maximum of shells for exchange
      integer, dimension(:), intent(in) :: nn !< Number of neighbour shells for exchange
      real(dblprec), dimension(:,:,:), intent(in) :: redcoord        !< Coordinates for Heisenberg exchange couplings
      real(dblprec), dimension(:,:,:,:,:), intent(in) :: jc          !< Exchange couplings (input)
      real(dblprec), dimension(:,:,:,:,:), intent(in) :: jcD         !< Exchange couplings (input) (DLM)
      real(dblprec), dimension(:,:,:,:,:,:), intent(in) :: jc_tens   !< Tensorial exchange (SKKR) couplings (input)
      ! DMI variables
      integer, intent(in) :: do_dm               !< Add Dzyaloshinskii-Moriya (DM) term to Hamiltonian (0/1)
      integer, intent(in) :: max_no_dmshells     !< Limit of number of shells for DM interactions
      integer, dimension(:), intent(in) :: dm_nn !< Number of neighbour shells for DM
      real(dblprec), dimension(:,:,:), intent(in) :: dm_redcoord    !< Coordinates for DM exchange couplings
      real(dblprec), dimension(:,:,:,:,:), intent(in) :: dm_inpvect !< DM couplings
      ! Anisotropy variables
      integer, intent(in) :: do_anisotropy                           !< Read anisotropy data (1/0)
      real(dblprec), intent(in) :: random_anisotropy_density         !< Density of random anisotropy
      integer, dimension(:,:), intent(in) :: anisotropytype          !< Type of anisotropies (0-2)
      integer, dimension(:,:), intent(in) :: anisotropytype_diff     !< Type of anisotropies (0-2)
      real(dblprec), dimension(:,:,:), intent(in) :: anisotropy      !< Input data for anisotropies
      real(dblprec), dimension(:,:,:), intent(in) :: anisotropy_diff !< Input data for anisotropies
      logical, intent(in) :: random_anisotropy  !< Distribute anisotropy randomly
      character(len=1), intent(in) :: mult_axis !< Flag to treat more than one anisotropy axis at the same time
      ! LSF variables
      integer, intent(in) :: mconf                 !< LSF moment ground state conf
      integer, intent(in) :: conf_num              !< Number of configurations for LSF
      logical, intent(in) :: map_multiple          !< Allow for multiple couplings between atoms i and j
      character(len=1), intent(in) ::  do_lsf      !< (Y/N) Do LSF for MC
      character(len=1), intent(in) ::  lsf_field   !< (T/L) LSF field term
      character(len=1), intent(in) ::  exc_inter   !< Interpolation of Jij between FM/DLM (Y/N)
      ! BQ interactions variables
      integer, intent(in) :: do_bq               !< Add biquadratic exchange (BQ) term to Hamiltonian (0/1)
      integer, intent(in) :: max_no_bqshells     !< Limit of number of shells for BQ interactions
      integer, dimension(:), intent(in) :: bq_nn !< Number of neighbour shells for BQ
      real(dblprec), dimension(:,:,:), intent(in) :: bq_redcoord  !< Coordinates for BQ exchange couplings
      real(dblprec), dimension(:,:,:,:), intent(in) :: jc_bq      !< Biquadratic exchange coupling (input)
      ! BQDM interactions variables
      integer, intent(in) :: do_biqdm               !< Add biquadratic DM (BIQDM) term to Hamiltonian (0/1)
      integer, intent(in) :: max_no_biqdmshells     !< Limit of number of shells for BIQDM interactions
      integer, dimension(:), intent(in) :: biqdm_nn !< Number of neighbour shells for BIQDM
      real(dblprec), dimension(:,:,:), intent(in) :: biqdm_redcoord      !< Coordinates for DM exchange couplings
      real(dblprec), dimension(:,:,:,:,:), intent(in) :: biqdm_inpvect   !< BIQDM couplings
      ! Pseudo-Dipolar interactions variables
      integer, intent(in) :: do_pd               !< Add Pseudo-Dipolar (PD) term to Hamiltonian (0/1)
      integer, intent(in) :: max_no_pdshells     !< Limit of number of shells for PD interactions
      integer, dimension(:), intent(in) :: pd_nn !< Number of neighbour shells for PD
      real(dblprec), dimension(:,:,:), intent(in) :: pd_redcoord      !< Coordinates for PD exchange couplings
      real(dblprec), dimension(:,:,:,:,:), intent(in) :: pd_inpvect   !< PD couplings
      ! Scalar chirality interactions variables
      integer, intent(in) :: do_chir               !< Add Scalar chirality (CHIR) term to Hamiltonian (0/1)
      integer, intent(in) :: max_no_chirshells     !< Limit of number of shells for CHIR interactions
      integer, dimension(:), intent(in) :: chir_nn !< Number of neighbour shells for CHIR
      real(dblprec), dimension(:,:,:,:), intent(in) :: chir_redcoord      !< Coordinates for CHIR exchange couplings
      real(dblprec), dimension(:,:,:,:), intent(in) :: chir_inpval    !< CHIR couplings
      ! Dipolar interactions variables
      integer, intent(in) :: do_dip !<  Calculate dipole-dipole contribution (0=Off, 1=Brute Force, 2=macrocell)
      ! Induced moments variables
      integer, dimension(NA,Nchmax), intent(in) :: ind_mom  !< Indication of whether a given moment is induced (1) or fixed (0)
      real(dblprec), intent(in) :: ind_tol                  !< Value for the tolerance between neighbouring shells
      character(len=1), intent(in) :: ind_mom_flag          !< Flag to indicate that there are induced moments being considered
      ! Cluster variables
      integer, intent(in) :: NA_clus         !< Number of atoms in the cluster unit cell
      integer, intent(in) :: NT_clus         !< Number of types of atoms in the cluster unit cell
      integer, intent(in) :: N1_clus         !< Number of cell repetitions in x direction for the cluster
      integer, intent(in) :: N2_clus         !< Number of cell repetitions in y direction for the cluster
      integer, intent(in) :: N3_clus         !< Number of cell repetitions in z direction for the cluster
      integer, intent(in) :: Natom_clus      !< Number of atoms in the cluster
      integer, intent(in) :: Nchmax_clus     !< Max number of chemical components on each site in cell for the cluster
      integer, intent(in) :: clus_expand
      integer, intent(in) :: Natom_full_clus !< Number of atoms for full system (=Natom if not dilute)
      integer, intent(inout) :: max_no_shells_clus    !< Calculated maximum of shells for exchange
      integer, intent(inout) :: max_no_dmshells_clus  !< Limit of number of shells for DM interactions
      integer, dimension(:), intent(in) :: NN_clus    !< Number of neighbour shells for exchange
      integer, dimension(:), intent(in) :: dm_nn_clus !< Number of neighbour shells for DM
      integer, dimension(:), intent(inout) :: achtype_clus       !< Chemical type of atoms (full list)
      integer, dimension(:), intent(inout) :: atype_ch_clus      !< Actual type of atom for dilute system
      integer, dimension(:), intent(inout) :: asite_ch_clus      !< Actual site of atom for dilute system
      integer, dimension(:), intent(inout) :: achem_ch_clus      !< Chemical type of atoms (reduced list)
      integer, dimension(:), intent(inout) :: acellnumb_clus     !< List for translating atom no. in full cell to actual cell
      integer, dimension(:), intent(inout) :: acellnumbrev_clus  !< List for translating atom no. in actual cell to full cell
      integer, dimension(Natom_clus), intent(inout) :: atype_clus    !< Type of atom for the cluster
      integer, dimension(Natom_clus), intent(inout) :: anumb_clus    !< Atom number in cell for the cluster
      integer, dimension(NA_clus), intent(inout) :: atype_inp_clus   !< Type of atom for the cluster (input)
      integer, dimension(NA_clus), intent(inout) :: anumb_inp_clus   !< Atom number in cell for the cluster (input)
      real(dblprec), dimension(3), intent(in) :: C1_clus    !< First lattice vector for the cluster
      real(dblprec), dimension(3), intent(in) :: C2_clus    !< Second lattice vector for the cluster
      real(dblprec), dimension(3), intent(in) :: C3_clus    !< Third lattice vector for the cluster
      real(dblprec), dimension(3,NA_clus), intent(inout) :: Bas_clus       !< Coordinates for basis atoms for the cluster
      real(dblprec), dimension(NA_clus,nchmax,conf_num), intent(inout) :: ammom_inp_clus !< Magnetic moment directions from input (for alloys)
      real(dblprec), dimension(:,:,:), intent(in)     :: redcoord_clus     !< Coordinates for Heisenberg exchange couplings for the cluster
      real(dblprec), dimension(:,:,:,:,:), intent(in) :: jc_clus           !< Exchange couplings (input) for the cluster
      real(dblprec), dimension(:,:,:), intent(in)     :: dm_redcoord_clus  !< Coordinates for DM exchange couplings
      real(dblprec), dimension(:,:,:,:,:), intent(in) :: dm_inpvect_clus   !< DM couplings for the cluster
      character(len=1), intent(in) :: do_cluster   !< Perform cluster simulation (Y/N)
      ! Anisotropy cluster variables
      integer, intent(in) :: do_anisotropy_clus                       !< Read anisotropy data for the cluster (1/0)
      real(dblprec), intent(in) :: random_anisotropy_density_clus     !< Density of random anisotropy for the cluster
      integer, dimension(:,:), intent(in) :: anisotropytype_clus      !< Type of anisotropies for cluster (0-2)
      integer, dimension(:,:), intent(in) :: anisotropytype_diff_clus !< Type of anisotropies for cluster (0-2)
      real(dblprec), dimension(:,:,:), intent(in) :: anisotropy_clus  !< Input data for anisotropies for cluster
      real(dblprec), dimension(:,:,:), intent(in) :: anisotropy_diff_clus !< Input data for anisotropies for cluster
      logical, intent(in) :: random_anisotropy_clus !< Distribute anisotropy randomly in the cluster
      character(len=1), intent(in) :: mult_axis_clus !< Flag to treat more than one anisotropy axis at the same time for the cluster
      ! Macrocell variables
      integer, intent(in) :: Num_macro !< Number of macrocells in the system
      integer, intent(in) :: max_num_atom_macro_cell  !< Maximum number of atoms in  a macrocell
      integer, dimension(Num_macro), intent(in) :: macro_nlistsize   !< Number of atoms per macrocell
      integer, dimension(Num_macro,max_num_atom_macro_cell), intent(in) :: macro_atom_nlist  !< List containing the information of which atoms are in a given macrocell
      integer, intent(in) :: block_size   !< Size of the blocking parameter for the macro cell creation
      character(len=1), intent(in) :: do_reduced   !< Use reduced formulation of Hamiltonian (Y/N)

      ! Misc variables
      integer, intent(in) :: do_prnstruct  !< Print Hamiltonian information (0/1)
      character, intent(in) :: do_sortcoup !< Sort the entries of ncoup arrays (Y/N)
      character(len=8), intent(in) :: simid            !< Name of simulation
      character(len=1), intent(in) :: print_dip_tensor !< Flag to print the dipole tensor
      character(len=1), intent(in) :: read_dipole !< Flag to read the dipole-dipole tensor
      character(len=30), intent(in) :: qdip_files !< Input file that contains the dipole-dipole tensor

      ! .. Local variables
      integer :: max_no_equiv       !< Calculated maximum of neighbours in one shell for exchange
      integer :: max_no_equiv_clus  !< Calculated maximum of neighbours in one shell for exchange in the cluster
      integer :: tot_max_no_neigh   !< Find which is the largest maximum number of neighbours between the cluster and the host
      integer :: i_all, i_stat

      ! Variable currently used only to match with the extended call arguments in setup_nm_nelem
      real(dblprec), dimension(27,NT,max_no_chirshells,NT,NT,48) :: chir_symtens

      integer, dimension(48,max_no_chirshells,na) :: nm_cell_symind  !< Indices for elements of the symmetry degenerate coupling tensor
      
      ! Set the Hamiltonian dimension
      if (do_reduced=='Y') then
         nHam=NA
      else
         nHam=Natom
      end if

      if (do_jtensor/=1) then

         ! Setup or import map
         if (.not.prev_map(simid)) then
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            ! Data for the cluster method
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            if (do_cluster=='Y') then
               ! The issue is that right now in the way it is defined, NA_clus=Natom_clus
               ! which is fine without repetitions, this should probably be fixed in the future
               ! the issue now is the C1,C2,C3, parameters, right now they are hard coded, the same
               ! than for the N1,N2, N3 belonging to the cluster
               write (*,'(2x,a)',advance='no') 'Set up neighbour map for exchange for the cluster'
               call setup_nm(Natom_clus,NT_clus,NA_clus,N1_clus,N2_clus,N3_clus,    &
                  C1_clus,C2_clus,C3_clus,'0','0','0',block_size,atype_clus,        &
                  Bas_clus,ham_clus%max_no_neigh_clus,max_no_shells_clus,           &
                  max_no_equiv_clus,sym,NN_clus,redcoord_clus,nm,nmdim,do_ralloy,   &
                  Natom_full_clus)
               write(*,'(a)') ' done'

               if (do_prnstruct==1.or.do_prnstruct==4) then
                  !  Print neighbor map
                  write (*,'(2x,a)',advance='no') 'Print neighbour map for exchange'
                  call clus_prnge(NA_clus,NT_clus,N1_clus,N2_clus,N3_clus,do_ralloy,&
                     Natom_clus,Nchmax_clus,Natom_full_clus,max_no_equiv_clus,      &
                     max_no_shells_clus,NN_clus,atype_clus,achtype_clus,            &
                     acellnumb_clus,acellnumbrev_clus,nmdim,nm,redcoord_clus,simid)
                  write(*,'(a)') ' done'
               end if
               ! Need to check if all the other variables such as ham%fs_ham%nlist and so on are consistent
               ! specially the atype_ch since right now the chermistry does not enter the cluster at all

               call allocate_cluster_hamiltoniandata(do_jtensor,Natom_clus,conf_num,&
                  ham_clus%max_no_neigh_clus,1)

               write (*,'(2x,a)',advance='no') 'Mount Heisenberg Hamiltonian for the cluster'
               call setup_neighbour_hamiltonian(Natom_clus,conf_num,NT_clus,NA_clus,&
                  Natom_clus,anumb_clus,atype_clus,ham_clus%max_no_neigh_clus,      &
                  max_no_equiv_clus,max_no_shells_clus,ham_clus%nlistsize_clus,     &
                  NN_clus,ham_clus%nlist_clus,ham_clus%ncoup_clus,nm,nmdim,jc_clus, &
                  ham%fs_nlist,ham%fs_nlistsize,do_ralloy,Natom_full_clus,          &
                  Nchmax_clus,atype_ch_clus,asite_ch_clus,achem_ch_clus,            &
                  ammom_inp_clus,1,1,do_sortcoup,do_lsf,ham%nind,lsf_field,         &
                  map_multiple)
               write(*,'(a)') ' done'

               call deallocate_nm()

               if (do_dm==1) then
                  ! Allocate and mount DM Hamiltonian
                  write (*,'(2x,a)',advance='no') 'Set up neighbour map for Dzyaloshinskii-Moriya exchange for the cluster'
                  call setup_nm(Natom_clus,NT_clus,NA_clus,N1_clus,N2_clus,N3_clus, &
                     C1_clus,C2_clus,C3_clus,BC1,BC2,BC3,block_size,atype_clus,     &
                     Bas_clus,ham_clus%max_no_dmneigh_clus,max_no_dmshells_clus,    &
                     max_no_equiv_clus,0,dm_nn_clus,dm_redcoord_clus,nm,nmdim,      &
                     do_ralloy,Natom_full_clus)
                  write(*,'(a)') ' done'

                  call allocate_cluster_dmhamiltoniandata(Natom_clus,               &
                     ham_clus%max_no_dmneigh_clus,1)
                  ! Transform data to general structure
                  write (*,'(2x,a)',advance='no') 'Mount Dzyaloshinskii-Moriya Hamiltonian for the cluster'
                  call setup_neighbour_hamiltonian(Natom_clus,conf_num,NT_clus,     &
                     NA_clus,Natom_clus,anumb_clus,atype_clus,                      &
                     ham_clus%max_no_dmneigh_clus,max_no_equiv_clus,                &
                     max_no_dmshells_clus,ham_clus%dmlistsize_clus,dm_nn_clus,      &
                     ham_clus%dmlist_clus,ham_clus%dm_vect_clus,nm,nmdim,           &
                     dm_inpvect_clus,ham%fs_nlist,ham%fs_nlistsize,do_ralloy,       &
                     Natom_full,Nchmax,atype_ch,asite_ch,achem_ch,ammom_inp_clus,3, &
                     1,do_sortcoup,do_lsf,ham%nind,lsf_field,map_multiple)

                     !Re-scale DMI if needed
                     if(dm_scale.ne.1.0_dblprec) ham%dm_vect=jij_scale*ham%dm_vect
                  write(*,'(a)') ' done'
                  call deallocate_nm()
               endif
               if (do_anisotropy_clus.eq.1) then
                  ! Setting up the anisotopies for the cluster
                  write(*,'(2x,a)',advance='no') "Set up anisotropies for the cluster"
                  call allocate_cluster_anisotropies(Natom_clus,mult_axis_clus,1)
                  call setup_anisotropies(Natom_clus,NA_clus,anumb_clus,            &
                     anisotropytype_clus,ham_clus%taniso_clus,ham_clus%eaniso_clus, &
                     ham_clus%kaniso_clus,ham_clus%sb_clus,anisotropy_clus,         &
                     mult_axis_clus,anisotropytype_diff_clus,                       &
                     ham_clus%taniso_diff_clus,ham_clus%eaniso_diff_clus,           &
                     ham_clus%kaniso_diff_clus,ham_clus%sb_diff_clus,               &
                     anisotropy_diff_clus,random_anisotropy_clus,                   &
                     random_anisotropy_density_clus,do_ralloy,Natom_full_clus,      &
                     Nchmax_clus,achem_ch_clus,ammom_inp_clus,mconf)
                  write(*,'(a)') ' done'
               endif
            endif
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            ! End of cluster method
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            !  Setup neighbor map
            write (*,'(2x,a)',advance='no') 'Set up neighbour map for exchange'
            call setup_nm(Natom,NT,NA,N1,N2,N3,C1,C2,C3,BC1,BC2,BC3,block_size,     &
               atype,Bas,ham%max_no_neigh,max_no_shells,max_no_equiv,sym,nn,        &
               redcoord,nm,nmdim,do_ralloy,Natom_full,acellnumb,atype_ch)
            write(*,'(a)') ' done'

            ! If one is doing the cluster method one could have that the impurity system
            ! has more neighbours than the host system, and then one should enforce
            ! that the number of neighbours corresponds to the largest system
            if (do_cluster=='Y') then
               tot_max_no_neigh=max(ham%max_no_neigh,ham_clus%max_no_neigh_clus,Natom_clus-1)
               ham%max_no_neigh=tot_max_no_neigh+clus_expand
            endif

            ! Transform data to general structure
            write (*,'(2x,a)',advance='no') 'Mount Heisenberg Hamiltonian'
            ! Actual call to allocate the hamiltonian data
            call allocate_hamiltoniandata(Natom,NA,nHam,conf_num,ham%max_no_neigh,  &
               do_jtensor,do_lsf,1,lsf_field,exc_inter)

            ! Setup the Hamiltonian look-up table
            call setup_aHam(Natom,anumb,do_reduced)
            !
            call setup_neighbour_hamiltonian(Natom,conf_num,NT,NA,nHam,anumb,atype, &
               ham%max_no_neigh,max_no_equiv,max_no_shells,ham%nlistsize,nn,        &
               ham%nlist,ham%ncoup,nm,nmdim,jc,ham%fs_nlist,ham%fs_nlistsize,       &
               do_ralloy,Natom_full,Nchmax,atype_ch,asite_ch,achem_ch,ammom_inp,1,1,&
               do_sortcoup,do_lsf,ham%nind,lsf_field,map_multiple)

            if (exc_inter=='Y') then
               call setup_neighbour_hamiltonian(Natom,conf_num,NT,NA,nHam,anumb,    &
                  atype,ham%max_no_neigh,max_no_equiv,max_no_shells,ham%nlistsize,  &
                  nn,ham%nlist,ham%ncoupD,nm,nmdim,jcD,ham%fs_nlist,                &
                  ham%fs_nlistsize,do_ralloy, Natom_full, Nchmax,atype_ch,asite_ch, &
                  achem_ch,ammom_inp,1,1,do_sortcoup,do_lsf,ham%nind,lsf_field,     &
                  map_multiple)
            endif

            write(*,'(a)') ' done'

            ! Deallocate the large neighbour map.
            call deallocate_nm()

            ! Re-scale Jij if needed
            if(jij_scale.ne.1.0_dblprec) ham%ncoup=jij_scale*ham%ncoup

            if(do_prnstruct==1.or.do_prnstruct==4) then
               write(*,'(2x,a)',advance='no') "Print exchange interaction strengths"
               call prn_exchange(NA,1,Natom,Nchmax,do_ralloy,Natom_full,            &
                  ham%max_no_neigh,simid,anumb,atype,ham%nlistsize,asite_ch,        &
                  achem_ch,ham%nlist,coord,ammom_inp,ham%ncoup,ham%aham)
               write(*,'(a)') ' done'
            end if

            if (ind_mom_flag=='Y') then
               write(*,'(2x,a)',advance='no') 'Set up neighbour map for induced moments'
               call induced_mapping(Natom,NT,NA,N1,N2,N3,sym,max_no_shells,nn,atype,&
                  ham%ind_nlistsize,ham%ind_nlist,ham%fix_nlistsize,ham%fix_nlist,  &
                  do_sortcoup,Nchmax,do_ralloy,     &
                  Natom_full,atype_ch,acellnumb,C1,C2,C3,Bas,BC1,BC2,BC3,ind_tol,   &
                  redcoord,ind_mom,block_size,ham%ind_list_full,                    &
                  ham%max_no_neigh_ind,ham%fix_num,ham%fix_list)

               write(*,'(a)') ' done'
               if (do_prnstruct==1) then
                  write(*,'(2x,a)',advance='no') 'Print neighbour map for induced moments'
                  call prn_ind_exchange(NA,Natom,Nchmax,do_ralloy,Natom_full,       &
                     ham%max_no_neigh_ind,anumb,achtype,ham%ind_nlistsize,ind_mom,  &
                     ham%ind_nlist,simid)
                  write(*,'(a)') ' done'
               endif
            endif
            if(do_prnstruct==5) then
               write(*,'(2x,a)',advance='no') "Print exchange interaction strengths (sparse format)"
               call prn_exchange_sparse(Natom,ham%max_no_neigh,ham%nlistsize,       &
                  ham%nlist,ham%ncoup,simid,1)
               write(*,'(a)') ' done'
            end if
         else ! If .out file exists
            write (*,'(2x,a)',advance='no') 'Importing exchange mapping'
            call read_exchange_getdim(Natom, ham%max_no_neigh, simid)
            call allocate_hamiltoniandata(Natom,NA,nHam,conf_num,ham%max_no_neigh,  &
            do_jtensor,do_lsf,1,lsf_field,exc_inter)
            ! Setup the Hamiltonian look-up table
            call setup_aHam(Natom,anumb,do_reduced)
            call read_exchange(NA,nHam,Natom,Nchmax,conf_num,do_ralloy,Natom_full, &
               ham%max_no_neigh,anumb,atype,asite_ch,achem_ch,ammom_inp,            &
               ham%nlistsize,ham%nlist,ham%ncoup,simid)
            write(*,'(a)') ' done'
         end if
      else ! If tensor
         !  Setup neighbor map and exchange tensor
         write (*,'(2x,a)',advance='no') 'Set up neighbour map for exchange'
         call setup_nm(Natom,NT,NA,N1,N2,N3,C1,C2,C3,BC1,BC2,BC3,block_size,atype,  &
            Bas,ham%max_no_neigh,max_no_shells,max_no_equiv,sym,nn,redcoord,nm,     &
            nmdim,do_ralloy,Natom_full,acellnumb,atype_ch)
         write(*,'(a)') ' done'

         ! Transform data to general structure
         write (*,'(2x,a)',advance='no') 'Mount Heisenberg Hamiltonian in tensorial (SKKR) form'
         call allocate_hamiltoniandata(Natom,NA,nHam,conf_num,ham%max_no_neigh,     &
            do_jtensor,do_lsf,1,lsf_field,exc_inter)

         ! Setup the Hamiltonian look-up table
         call setup_aHam(Natom,anumb,do_reduced)

         call setup_neighbour_hamiltonian(Natom,conf_num,NT,NA,nHam, anumb,atype,   &
            ham%max_no_neigh,max_no_equiv,max_no_shells,ham%nlistsize,nn,ham%nlist, &
            ham%j_tens,nm,nmdim,jc_tens,ham%fs_nlist,ham%fs_nlistsize,do_ralloy,    &
            Natom_full,Nchmax,atype_ch,asite_ch,achem_ch,ammom_inp,9,1,do_sortcoup, &
            do_lsf,ham%nind,lsf_field,map_multiple)
         write(*,'(a)') ' done'

         ! Deallocate the large neighbour map.
         call deallocate_nm()

         if(do_prnstruct==1.or.do_prnstruct==4) then
            write(*,'(2x,a)',advance='no') "Print exchange interaction strenghts"
            call prn_exchange(NA,9,Natom,Nchmax,do_ralloy,Natom_full,               &
               ham%max_no_neigh,simid,anumb,atype,ham%nlistsize,asite_ch,achem_ch,  &
               ham%nlist,coord,ammom_inp,ham%j_tens,ham%aham)
            write(*,'(a)') ' done'
         end if
      end if
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! End of tensor if
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Anisotropies
      if (do_anisotropy==1) then
         write(*,'(2x,a)',advance='no') "Set up anisotropies"
         call allocate_anisotropies(Natom,mult_axis,1)
         call setup_anisotropies(Natom,NA,anumb,anisotropytype,ham%taniso,          &
            ham%eaniso,ham%kaniso,ham%sb,anisotropy,mult_axis,anisotropytype_diff,  &
            ham%taniso_diff,ham%eaniso_diff,ham%kaniso_diff,ham%sb_diff,            &
            anisotropy_diff,random_anisotropy,random_anisotropy_density,do_ralloy,  &
            Natom_full,Nchmax,achem_ch,ammom_inp,mconf)
         write(*,'(a)') ' done'
         ! Print anisotropies
         if(do_prnstruct==1 .and. do_anisotropy == 1) then
            write(*,'(2x,a)',advance='no') "Print anisotropies"
            call prn_anisotropy(Natom,NA,anisotropytype,ham%taniso,ham%eaniso,      &
            ham%kaniso,simid,Nchmax)
               write(*,'(a)') ' done'
            end if
      endif
      ! Setup the dipole-dipole interactions
      if (do_dip>0) then
         call dipole_setup(NA,N1,N2,N3,Natom,do_dip,Num_macro,Mensemble,            &
            max_num_atom_macro_cell,macro_nlistsize,macro_atom_nlist,alat,C1,C2,C3, &
            Bas,coord,simid,read_dipole,print_dip_tensor,qdip_files,ham%Qdip,ham%Qdip_macro)
      endif

      ! Calculate mean field estimate of Tc
      if(do_jtensor/=1.and.int(N1*N2*N3).ne.1) then
         write(*,'(2x,a)',advance='no') "Calculate mean-field estimate of Tc: "
         call estimate_tc_mfa(Natom,nHam,Nchmax,conf_num,NA,atype,anumb,Natom_full, &
            atype_ch,asite_ch,achem_ch,ham%nlistsize,ham%nlist,ham%max_no_neigh,    &
            ham%ncoup,ammom_inp,mconf,do_ralloy)
      end if

      if(do_dm==1) then
         ! Allocate and mount DM Hamiltonian
         write (*,'(2x,a)',advance='no') 'Set up neighbour map for Dzyaloshinskii-Moriya exchange'
         call setup_nm(Natom,NT,NA,N1,N2,N3,C1,C2,C3,BC1,BC2,BC3,block_size,atype,  &
            Bas,ham%max_no_dmneigh,max_no_dmshells,max_no_equiv,0,dm_nn,dm_redcoord,&
            nm,nmdim,do_ralloy,Natom_full,acellnumb,atype_ch)
         write(*,'(a)') ' done'

         ! If one is doing the cluster method one could have that the impurity system
         ! has more neighbours than the host system, and then one should enforce
         ! that the number of neighbours corresponds to the largest system
         if (do_cluster=='Y') then
            tot_max_no_neigh=max(ham%max_no_dmneigh,ham_clus%max_no_dmneigh_clus,Natom_clus-1)
            ham%max_no_dmneigh=tot_max_no_neigh+clus_expand
         endif

         call allocate_dmhamiltoniandata(Natom,nHam,ham%max_no_dmneigh,1)

         ! Transform data to general structure
         write (*,'(2x,a)',advance='no') 'Mount Dzyaloshinskii-Moriya Hamiltonian'
         call setup_neighbour_hamiltonian(Natom,conf_num,NT,NA,nHam,anumb,atype,    &
            ham%max_no_dmneigh,max_no_equiv,max_no_dmshells,ham%dmlistsize,dm_nn,   &
            ham%dmlist,ham%dm_vect,nm,nmdim,dm_inpvect,ham%fs_nlist,                &
            ham%fs_nlistsize,do_ralloy,Natom_full,Nchmax,atype_ch,asite_ch,achem_ch,&
            ammom_inp,3,1,do_sortcoup,do_lsf,ham%nind,lsf_field,map_multiple)
         write(*,'(a)') ' done'

         ! Deallocate the large neighbour map.
         call deallocate_nm()

         ! Print DM interactions
         if((do_prnstruct==1.or.do_prnstruct==4).and.do_cluster/='Y') then
            write(*,'(2x,a)',advance='no') "Print Dzyaloshinskii-Moriya interactions"
            call prn_dmcoup(NA,Natom,Nchmax,do_ralloy,Natom_full,ham%max_no_dmneigh,&
               anumb,atype,ham%dmlistsize,asite_ch,achem_ch,ham%dmlist,coord,       &
               ammom_inp,ham%dm_vect,simid)
            write(*,'(a)') ' done'
         end if
      end if

      ! If LSF
      if(do_lsf=='Y') then
         write (*,'(2x,a)',advance='no') 'Set up moments map for Longitudial Fluctuation'
         call LSF_datareshape(NA, Nchmax, conf_num)
         write(*,'(a)') ' done'
      endif

      if(do_pd==1) then
         ! Allocate and mount PD Hamiltonian
         write (*,'(2x,a)',advance='no') 'Set up neighbour map for Pseudo-Dipolar exchange'
         call setup_nm(Natom,NT,NA,N1,N2,N3,C1,C2,C3,BC1,BC2,BC3,block_size,atype,  &
            Bas,ham%nn_pd_tot,max_no_pdshells,max_no_equiv,sym,pd_nn,pd_redcoord,nm,&
            nmdim,do_ralloy,Natom_full,acellnumb,atype_ch)
         write(*,'(a)') ' done'

         call allocate_pdhamiltoniandata(Natom,nHam,ham%nn_pd_tot,1)

         ! Transform data to general structure
         write (*,'(2x,a)',advance='no') 'Mount Pseudo-Dipolar Hamiltonian'
         call setup_neighbour_hamiltonian(Natom,conf_num,NT,NA,nHam,anumb,atype,    &
            ham%nn_pd_tot,max_no_equiv,max_no_pdshells,ham%pdlistsize,pd_nn,        &
            ham%pdlist,ham%pd_vect,nm,nmdim,pd_inpvect,ham%fs_nlist,                &
            ham%fs_nlistsize,do_ralloy,Natom_full,Nchmax,atype_ch,asite_ch,achem_ch,&
            ammom_inp,6,1,do_sortcoup,do_lsf,ham%nind,lsf_field,map_multiple)
         write(*,'(a)') ' done'

         ! Deallocate the large neighbour map.
         call deallocate_nm()

         ! Print PD interactions
         if(do_prnstruct==1.or.do_prnstruct==4) then
            write(*,'(2x,a)',advance='no') "Print Pseudo-Dipolar interactions"
            call prn_pdcoup(Natom,ham%nn_pd_tot,ham%pdlistsize,ham%pdlist,          &
               ham%pd_vect,simid)
            write(*,'(a)') ' done'
         end if

      end if

      if(do_chir==1) then

         call ErrorHandling_missing('Chiral interactions')
         !!! ! Allocate and mount chir Hamiltonian
         !!! write (*,'(2x,a)',advance='no') 'Set up neighbour map for MML interaction'
         !!! call setup_nm_nelem(Natom,NT,NA,N1,N2,N3,C1,C2,C3,BC1,BC2,BC3,atype,Bas,   &
         !!!    ham%nn_chir_tot,max_no_chirshells,max_no_equiv,0,chir_nn,chir_redcoord, &
         !!!    nme,nmdim,2,do_ralloy,Natom_full,acellnumb,atype_ch)

         !!! write(*,'(a)') ' done'

         !!! call allocate_chirhamiltoniandata(Natom,nHam,ham%nn_chir_tot,1)

         !!! ! Transform data to general structure
         !!! write (*,'(2x,a)',advance='no') 'Mount chir Hamiltonian'
         !!! call setup_neighbour_latticehamiltonian(Natom,NT,NA,anumb,atype,           &
         !!!    ham%nn_chir_tot,max_no_equiv,max_no_chirshells,2,1,1,ham%chirlistsize,  &
         !!!    chir_nn,ham%chirlist,ham%chir_coup,nme,nmdim,chir_inpval,do_sortcoup,   &
         !!!    'N',do_ralloy,Natom_full,Nchmax,atype_ch,asite_ch,achem_ch,ammom_inp,1,1)
         !!! write(*,'(a)') ' done'

         !!! ham%max_no_chirneigh=ham%nn_chir_tot
         !!! ! Deallocate the large neighbour map.
         !!! !call deallocate_nme()
         !!! i_all=-product(shape(nme))*kind(nme)
         !!! deallocate(nme,stat=i_stat)
         !!! call memocc(i_stat,i_all,'nme','setup_hamiltonian')

         !!! call prn_chircoup(Natom,ham%max_no_chirneigh,ham%chirlistsize,ham%chirlist,&
         !!!    ham%chir_coup,simid)

         !!! ham%chir_coup=ham%chir_coup*2.0_dblprec*mry/mub

         !print *,'ham%nn_chir_tot', ham%nn_chir_tot
         !print *,'ham%chirlistsize',ham%chirlistsize
         !print *,'ham%chirlist',ham%chirlist
         !print *,'ham%chir_coup',ham%chir_coup
     !!!    ! Rescale chir couplings if needed
     !!!    if(chir_scale.ne.1.0_dblprec) chir_tens=chir_scale*chir_tens

     !!!    ! Print chir interactions
     !!!    if(do_prnstruct==1.or.do_prnstruct==4) then
     !!!       write(*,'(2x,a)',advance='no') "Print chir interactions"
     !!!       call prn_chircoup(Natom, nn_chir_tot, chir_listsize, chir_list, chir_tens, simid)
     !!!       write(*,'(a)') ' done'
     !!!     call chk_chircoup(Natom, nn_chir_tot, chir_listsize, chir_list, chir_tens, simid)
     !!!    end if

      end if

      if(do_biqdm==1) then
         ! Allocate and mount BIQDM Hamiltonian
         write (*,'(2x,a)',advance='no') 'Set up neighbour map for BIQDM exchange'
         call setup_nm(Natom,NT,NA,N1,N2,N3,C1,C2,C3,BC1,BC2,BC3,block_size,atype,  &
            Bas,ham%nn_biqdm_tot,max_no_biqdmshells,max_no_equiv,0,biqdm_nn,        &
            biqdm_redcoord,nm,nmdim,do_ralloy,Natom_full,acellnumb,atype_ch)
         write(*,'(a)') ' done'

         call allocate_biqdmhamiltoniandata(Natom,nHam,ham%nn_biqdm_tot,1)
         ! Transform data to general structure
         write (*,'(2x,a)',advance='no') 'Mount BIQDM Hamiltonian'

         call setup_neighbour_hamiltonian(Natom,conf_num,NT,NA,nHam,anumb,atype,    &
            ham%nn_biqdm_tot,max_no_equiv,max_no_biqdmshells,ham%biqdmlistsize,     &
            biqdm_nn,ham%biqdmlist,ham%biqdm_vect,nm,nmdim,biqdm_inpvect,           &
            ham%fs_nlist,ham%fs_nlistsize,do_ralloy,Natom_full,Nchmax,atype_ch,     &
            asite_ch,achem_ch,ammom_inp,1,2,do_sortcoup,do_lsf,ham%nind,lsf_field,  &
            map_multiple)
         write(*,'(a)') ' done'

         ! Deallocate the large neighbour map.
         call deallocate_nm()

         ! Print BIQDM interactions
         if(do_prnstruct==1.or.do_prnstruct==4) then
            write(*,'(2x,a)',advance='no') "Print BIQDM interactions"
            call prn_biqdmcoup(Natom,ham%nn_biqdm_tot,ham%biqdmlistsize,            &
               ham%biqdmlist,ham%biqdm_vect, simid)
            write(*,'(a)') ' done'
         end if

      end if

      if(do_bq==1) then

         ! Allocate and mount BQ Hamiltonian
         write (*,'(2x,a)',advance='no') 'Set up neighbour map for biquadratic exchange'
         call setup_nm(Natom,NT,NA,N1,N2,N3,C1,C2,C3,BC1,BC2,BC3,block_size,atype,  &
            Bas,ham%nn_bq_tot,max_no_bqshells,max_no_equiv,sym,bq_nn,bq_redcoord,nm,&
            nmdim,do_ralloy,Natom_full,acellnumb,atype_ch)
         write(*,'(a)') ' done'

         call allocate_bqhamiltoniandata(Natom,nHam,ham%nn_bq_tot,1)

         ! Transform data to general structure
         write (*,'(2x,a)',advance='no') 'Mount biquadratic exchange Hamiltonian'
         call setup_neighbour_hamiltonian(Natom,conf_num,NT,NA,nHam,anumb,atype,    &
            ham%nn_bq_tot,max_no_equiv,max_no_bqshells,ham%bqlistsize,bq_nn,        &
            ham%bqlist,ham%j_bq,nm,nmdim,jc_bq,ham%fs_nlist,ham%fs_nlistsize,       &
            do_ralloy,Natom_full,Nchmax,atype_ch,asite_ch,achem_ch,ammom_inp,1,2,   &
            do_sortcoup,do_lsf,ham%nind,lsf_field,map_multiple)
         write(*,'(a)') ' done'

         ! Deallocate the large neighbour map.
         call deallocate_nm()

         ! Print BQ interactions
         if(do_prnstruct==1.or.do_prnstruct==4) then
            write(*,'(2x,a)',advance='no') "Print biquadratic exchange interactions"
            call prn_bqcoup(Natom,ham%nn_bq_tot,ham%bqlistsize,ham%bqlist,ham%j_bq, &
               simid)
            write(*,'(a)') ' done'
         end if
      endif

   contains

      ! Setup the Hamiltonian look-up table
      subroutine setup_aHam(Natom,anumb,do_reduced)
         !
         implicit none
         !
         integer, intent(in) :: Natom !< Number of atoms in system
         integer, dimension(Natom), intent(in) :: anumb !< Number of atom in cell
         character(len=1), intent(in) :: do_reduced       !< Use reduced formulation of Hamiltonian (T/F)
         !
         integer :: ia
         !
         if (do_reduced=='Y') then
            do ia=1,Natom
               ham%aHam(ia)=anumb(ia)
            end do
         else
            do ia=1,Natom
               ham%aHam(ia)=ia
            end do
         end if
         !
      end subroutine setup_aHam

      !> Evaluates if exchange map already exists
      logical function prev_map(simid)
         !
         implicit none
         !
         character(len=8),intent(in) :: simid !< Name of simulation
         character(len=20) :: filn

         !.. Executable statements
         write (filn,'(''struct.'',a8,''.out'')') simid
         inquire(file=filn,exist=prev_map)
         return
      end function prev_map

   end subroutine setup_hamiltonian

   !> Sets up single ion anisotropies
   subroutine setup_anisotropies(Natom,NA,anumb,anisotropytype,taniso,eaniso,kaniso,&
      sb,anisotropy,mult_axis,anisotropytype_diff,taniso_diff,eaniso_diff,          &
      kaniso_diff,sb_diff,anisotropy_diff,random_anisotropy,                        &
      random_anisotropy_density,do_ralloy,Natom_full,Nchmax,achem_ch,ammom_inp,mconf)

      use Constants
      use RandomNumbers, only: rng_uniform
      !
      implicit none
      !
      integer, intent(in) :: Natom                                      !< Number of atoms in system
      integer, intent(in) :: NA                                         !< Number of atoms in one cell
      integer, dimension(Natom), intent(in) :: anumb                    !< Atom number in cell
      integer, intent(in) :: Nchmax                                     !< Max number of chemical components on each site in cell
      integer, dimension(:,:), intent(in) :: anisotropytype             !< Type of anisotropies (0-2)
      integer, dimension(:,:), intent(in) :: anisotropytype_diff        !< Type of anisotropies (0-2)
      real(dblprec), dimension(:,:,:), intent(in) :: anisotropy         !< Input data for anisotropies
      real(dblprec), dimension(:,:,:), intent(in) :: anisotropy_diff    !< Input data for anisotropies
      logical, intent(in) :: random_anisotropy                          !< Distribute anisotropy randomly
      real(dblprec), intent(in) :: random_anisotropy_density            !< Density of random anisotropy
      integer, dimension(:), intent(out) :: taniso                      !< Type of anisotropy (0-2)
      integer, dimension(:), intent(out) :: taniso_diff                 !< Type of anisotropy (0-2)
      real(dblprec), dimension(:,:), intent(out) :: eaniso              !< Unit anisotropy vector
      real(dblprec), dimension(:,:), intent(out) :: eaniso_diff         !< Unit anisotropy vector
      real(dblprec), dimension(:,:), intent(out) :: kaniso              !< Anisotropy constant
      real(dblprec), dimension(:,:), intent(out) :: kaniso_diff         !< Anisotropy constant
      real(dblprec), dimension(:), intent(out) :: sb                    !< Ratio between the anisotropies
      real(dblprec), dimension(:), intent(out) :: sb_diff               !< Ratio between the anisotropies
      integer, intent(in) :: Natom_full                                 !< Number of atoms for full system (=Natom if not dilute)
      integer, intent(in) :: do_ralloy                                  !< Random alloy simulation (0/1)
      integer, dimension(:), intent(in) :: achem_ch                     !< Chemical type of atoms (reduced list)
      real(dblprec), dimension(:,:,:), intent(in) :: ammom_inp          !< Magnetic moment directions from input (for alloys)
      character(len=1), intent(in) :: mult_axis                         !< Flag to treat more than one anisotropy axis at the same time
      integer, intent(in) :: mconf                                      !< LSF ground state configuration

      integer :: i,j
      real(dblprec) :: fc, anorm,anorm_diff,rn(1)

      fc = mry/mub

      ! Anisotropy
      if(do_ralloy==0.and..not.random_anisotropy) then
         do i=1, Natom
            taniso(i)=anisotropytype(anumb(i),1)

            anorm=anisotropy(anumb(i),3,1)**2+anisotropy(anumb(i),4,1)**2+anisotropy(anumb(i),5,1)**2
            eaniso(1:3,i)=anisotropy(anumb(i),3:5,1)/sqrt(anorm+1.0d-15)
            if(taniso(i)==2) then
               kaniso(1,i)=fc*anisotropy(anumb(i),1,1)/ammom_inp(anumb(i),1,1)**4
               kaniso(2,i)=fc*anisotropy(anumb(i),2,1)/ammom_inp(anumb(i),1,1)**6
            else
               kaniso(1,i)=fc*anisotropy(anumb(i),1,1)/ammom_inp(anumb(i),1,1)**2
               kaniso(2,i)=fc*anisotropy(anumb(i),2,1)/ammom_inp(anumb(i),1,1)**4
            end if

            sb(i)=anisotropy(anumb(i),6,1)
         end do

         if (mult_axis=='Y') then
            do i=1, Natom
               taniso_diff(i)=anisotropytype_diff(anumb(i),1)
               anorm_diff=anisotropy_diff(anumb(i),3,1)**2+anisotropy_diff(anumb(i),4,1)**2+anisotropy_diff(anumb(i),5,1)**2
               eaniso_diff(1:3,i)=anisotropy_diff(anumb(i),3:5,1)/sqrt(anorm_diff)
               kaniso_diff(1:2,i)=fc*anisotropy_diff(anumb(i),1:2,1)/ammom_inp(anumb(i),1,1)**2
               sb_diff(i)=anisotropy_diff(anumb(i),6,1)
            end do
         endif

      else if(do_ralloy==0.and.random_anisotropy) then
         taniso=0
         do j=1, Natom*int(random_anisotropy_density)
            call rng_uniform(rn,1)
            i = int(Natom*rn(1))+1
            taniso(i)=anisotropytype(anumb(i),1)

            anorm=anisotropy(anumb(i),3,1)**2+anisotropy(anumb(i),4,1)**2+anisotropy(anumb(i),5,1)**2
            eaniso(1:3,i)=anisotropy(anumb(i),3:5,1)/sqrt(anorm)
            if(taniso(i)==2) then
               kaniso(1,i)=fc*anisotropy(anumb(i),1,1)/ammom_inp(anumb(i),1,1)**4
               kaniso(2,i)=fc*anisotropy(anumb(i),2,1)/ammom_inp(anumb(i),1,1)**6
            else
               kaniso(1,i)=fc*anisotropy(anumb(i),1,1)/ammom_inp(anumb(i),1,1)**2
               kaniso(2,i)=fc*anisotropy(anumb(i),2,1)/ammom_inp(anumb(i),1,1)**4
            end if

            sb(i)=anisotropy(anumb(i),6,1)
         end do
      else
         do i=1, Natom
            taniso(i)=anisotropytype(anumb(i),achem_ch(i))

            anorm=anisotropy(anumb(i),3,achem_ch(i))**2+anisotropy(anumb(i),4,achem_ch(i)) &
               **2+anisotropy(anumb(i),5,achem_ch(i))**2
            eaniso(1:3,i)=anisotropy(anumb(i),3:5,achem_ch(i))/sqrt(anorm)
            if(taniso(i)==2) then
               kaniso(1,i)=fc*anisotropy(anumb(i),1,achem_ch(i))/ammom_inp(anumb(i),achem_ch(i),mconf)**4
               kaniso(2,i)=fc*anisotropy(anumb(i),2,achem_ch(i))/ammom_inp(anumb(i),achem_ch(i),mconf)**6
            else
               kaniso(1,i)=fc*anisotropy(anumb(i),1,achem_ch(i))/ammom_inp(anumb(i),achem_ch(i),mconf)**2
               kaniso(2,i)=fc*anisotropy(anumb(i),2,achem_ch(i))/ammom_inp(anumb(i),achem_ch(i),mconf)**4

            end if

            sb(i)=anisotropy(anumb(i),6,achem_ch(i))
         end do

         if (mult_axis=='Y') then
            taniso_diff(i)=anisotropytype_diff(anumb(i),achem_ch(i))

            anorm_diff=anisotropy_diff(anumb(i),3,achem_ch(i))**2+anisotropy_diff(anumb(i),4,achem_ch(i)) &
               **2+anisotropy_diff(anumb(i),5,achem_ch(i))**2
            eaniso_diff(1:3,i)=anisotropy_diff(anumb(i),3:5,achem_ch(i))/sqrt(anorm_diff)
            kaniso_diff(1:2,i)=fc*anisotropy_diff(anumb(i),1:2,achem_ch(i))/ammom_inp(anumb(i),achem_ch(i),mconf)**2

            sb_diff(i)=anisotropy_diff(anumb(i),6,achem_ch(i))

         endif
      end if

   end subroutine setup_anisotropies

   !----------------------------------------------------------------------------
   ! SUBROUTINE: setup_neighbour_hamiltonian
   !> @brief Mounts the actual Heisenberg Hamiltonian (or DM Hamiltonian)
   !> @details This routine actually setups any kind of pairwise interaction between
   !> magnetic atoms, i.e. Heisenberg, DMI, BQ, etc.
   !> Specifically it setups the needed reduced lists and couplings for the system
   !----------------------------------------------------------------------------
   subroutine setup_neighbour_hamiltonian(Natom,conf_num,NT,NA,nHam,anumb,atype,    &
      max_no_neigh,max_no_equiv,max_no_shells,nlistsize,nn,nlist,ncoup,nm,nmdim,xc, &
      fs_nlist,fs_nlistsize,do_ralloy,Natom_full,Nchmax,atype_ch,asite_ch,achem_ch, &
      ammom_inp,hdim,lexp,do_sortcoup,do_lsf,nind,lsf_field,map_multiple)
      !
      use Constants
      use Sorting, only : MergeSortIR
      !
      implicit none
      !
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
      !
      integer :: i, j, k, l,iconf, i_stat
      integer :: jatom
      real(dblprec) :: fc,fc2
      real(dblprec), dimension(hdim) :: tempcoup
      integer :: tempn
      integer :: ncount,ncounter
      logical :: exis

      ncoup = 0.0_dblprec
      ! Factors for mRy energy conversion
      fc = mry/mub
      fc2 = 2.0_dblprec*mry/mub

      ! Loop over atoms
      !$omp parallel do default(shared) private(i,k,j,ncount,ncounter,exis,l,jatom,iconf)
      do i=1, Natom
         ncount = 1
         ncounter = 1
         ! Shells
         do k=1, Nn(atype(i))
            ! Sites in shell
            do j=1, nmdim(k,i)
               ! Existing coupling
               if (nm(i,k,j)>0) then
                  exis=.false.
                  do l=1,ncount-1
                     if(nlist(l,i)==nm(i,k,j)) exis=.true.
                  end do
                  if(.not.exis.or.map_multiple) then
                     nlist(ncount,i) = nm(i,k,j)
                     if (i<=nHam) then
                        do iconf=1,conf_num
                           if (do_ralloy==0) then
                              if(abs(ammom_inp(anumb(i),1,iconf)*ammom_inp(anumb(nm(i,k,j)),1,iconf))<1e-6) then
                                 ncoup(:,ncount,i,iconf)=0.0_dblprec
                              else
                                 ncoup(:,ncount,i,iconf) = xc(:,atype(i),k,1,1,iconf) * fc2 &
                                    / ammom_inp(anumb(i),1,iconf)**lexp / ammom_inp(anumb(nm(i,k,j)),1,iconf)**lexp
                              endif
                           else
                              !
                              jatom = nm(i,k,j)
                              if(abs(ammom_inp(asite_ch(i),achem_ch(i),iconf)*ammom_inp(asite_ch(jatom),achem_ch(jatom),iconf))<1e-6) then
                                 ncoup(:,ncount,i,iconf)=0.0_dblprec
                              else
                                 ncoup(:,ncount,i,iconf) = xc(:,atype_ch(i),k,achem_ch(i),achem_ch(jatom),iconf) * &
                                    fc2 / ammom_inp(asite_ch(i),achem_ch(i),iconf)**lexp  / ammom_inp(asite_ch(jatom),achem_ch(jatom),iconf)**lexp
                              endif
                           end if
                        end do
                     end if
                     if (k==1 .and. do_lsf=='Y' .and. lsf_field=='L') then  !< conditional adding fist-shell nlist to fs_nlist
                        fs_nlist(ncount,i) = nlist(ncount,i)
                        ncounter = ncounter + 1
                     end if
                     ncount = ncount + 1
                  end if
               end if
            end do
         end do
         if (i<=nHam) nlistsize(i) = ncount-1
         if (do_lsf=='Y' .and. lsf_field=='L') fs_nlistsize(i) = ncounter-1
      end do
      !omp end parallel do

      if(do_sortcoup == 'Y'.and. nHam==Natom) then
         !$omp parallel do default(shared) private(i,j,k,tempn,tempcoup,iconf)
         do i=1,Natom
            ! sort neighbour list - bubble sort...
            do j=1,nlistsize(i)
               do k=1,nlistsize(i)-j
                  if ( nlist(k,i) .gt. nlist(k+1,i) ) then
                     tempn        = nlist(k,i)
                     nlist(k,i)   = nlist(k+1,i)
                     nlist(k+1,i) = tempn
                     do iconf=1,conf_num
                        tempcoup     = ncoup(:,k,i,iconf)
                        ncoup(:,k,i,iconf)   = ncoup(:,k+1,i,iconf)
                        ncoup(:,k+1,i,iconf) = tempcoup
                     enddo
                  end if
               end do
            end do
         end do
         !$omp end parallel do

         if(do_lsf=='Y' .and. lsf_field=='L') then
            !$omp parallel do default(shared) private(i,j,k,tempn)
            do i=1,Natom
               ! In the same way sort fs_nlistsize
               do j=1,fs_nlistsize(i)
                  do k=1,fs_nlistsize(i)-j
                     if ( fs_nlist(k,i) .gt. fs_nlist(k+1,i) ) then
                        tempn        = fs_nlist(k,i)
                        fs_nlist(k,i)   = fs_nlist(k+1,i)
                        fs_nlist(k+1,i) = tempn
                     end if
                  end do
               end do
            end do
            !$omp end parallel do
         endif
      end if

      !> link atom-index from fs_nlist to nlist, in order to use ncoup

      if (do_lsf=='Y' .and. lsf_field=='L') then
         if (.not. allocated(nind)) then
            allocate(nind(maxval(fs_nlistsize),Natom),stat=i_stat)
            call memocc(i_stat,product(shape(nind))*kind(nind),'nind','setup_neighbour_hamiltonian')
         endif
         !$omp parallel do default(shared) private(i)
         do i = 1,Natom
            do j=1,fs_nlistsize(i)
               do k = 1, nlistsize(i)
                  if (nlist(k,i) == fs_nlist(j,i)) then
                     nind(j,i) = k
                     exit
                  endif
               enddo
            end do
            nind(:,i) = nind(1:fs_nlistsize(i),i)
         end do
         !$omp end parallel do
      end if

   end subroutine setup_neighbour_hamiltonian

   !> Finds the dimension of the Heisenberg Hamiltonian if read from file
   subroutine read_exchange_getdim(Natom,max_no_neigh,simid)
      !
      !.. Implicit declarations
      implicit none

      integer, intent(in) :: Natom !< Number of atoms in system
      integer, intent(out) :: max_no_neigh !< Calculated maximum of neighbours for exchange
      character(len=8),intent(in) :: simid !< Name of simulation

      !.. Local variables
      integer :: ntest
      character(len=20) :: filn
      character(len=18) :: dummy
      character(len=27) :: dummy2
      !.. Executable statements
      write (filn,'(''struct.'',a8,''.out'')') simid
      open(ifileno, file=filn)

      ! Find max no. neighbours
      max_no_neigh=1
      read (ifileno,*)
      read (ifileno,'(a18,1x,i20)') dummy,ntest
      if (ntest.ne.Natom) stop "Number of atom of previous system not the same"
      read(ifileno,'(a27,1x,i20)') dummy2,max_no_neigh
      close(ifileno)
      10001 format (5x,i8,13x,i7)

   end subroutine read_exchange_getdim

   !----------------------------------------------------------------------------
   !> Reads the Heisenberg Hamiltonian from file
   !> @note Jonathan Chico: Modified such that it corresponds to the new units
   !> of the struct1file
   !----------------------------------------------------------------------------
   subroutine read_exchange(NA,nHam,Natom,Nchmax,conf_num,do_ralloy,Natom_full,     &
      max_no_neigh,anumb,atype,asite_ch,achem_ch,ammom_inp,nlistsize,nlist,ncoup,   &
      simid)
      !
      use Constants, only : mub,mry
      !.. Implicit declarations
      implicit none

      integer, intent(in) :: NA           !< Number of atoms in one cell
      integer, intent(in) :: nHam         !< Number of atoms in Hamiltonian
      integer, intent(in) :: Natom        !< Number of atoms in system
      integer, intent(in) :: Nchmax       !< Number of chemical components in system
      integer, intent(in) :: conf_num     !< Number of configurations for LSF
      integer, intent(in) :: do_ralloy    !< Random alloy simulation (0/1)
      integer, intent(in) :: Natom_full   !< Number of atoms for full system (=Natom if not dilute)
      integer, intent(in) :: max_no_neigh !< Calculated maximum of neighbours for exchange
      integer, dimension(Natom), intent(in)        :: anumb    !< Atom number in cell
      integer, dimension(Natom), intent(in)        :: atype    !< Type of atom
      integer, dimension(Natom_full), intent(in)   :: asite_ch !< Actual site of atom for dilute system
      integer, dimension(Natom_full), intent(in)   :: achem_ch !< Chemical type of atoms (reduced list)
      real(dblprec), dimension(NA,Nchmax,conf_num), intent(in) :: ammom_inp   !< Magnetic moment directions from input (for alloys)
      integer, dimension(nHam), intent(out) :: nlistsize !< Size of neighbour list for Heisenberg exchange couplings
      integer, dimension(max_no_neigh, Natom), intent(out) :: nlist !< Neighbour list for Heisenberg exchange couplings
      real(dblprec), dimension(max_no_neigh, nHam, conf_num), intent(out) :: ncoup !< Heisenberg exchange couplings
      character(len=8),intent(in) :: simid !< Name of simulation

      !.. Local variables
      integer :: i,j,jatom,iatom,i_err
      real(dblprec) :: fc2
      real(dblprec) :: dum,ncoup_tmp
      character(len=20) :: filn
      i_err=0
      fc2 = 2.0_dblprec*mry/mub
      !.. Executable statements
      if (conf_num>1) then
         write(*,*) 'Restart does not work with LSF right now'
         stop
      endif

      write (filn,'(''struct.'',a8,''.out'')') simid
      open(ifileno, file=filn)

      ncoup=0.0_dblprec
      read (ifileno,*)
      read (ifileno,*)
      read (ifileno,*)
      read (ifileno,*)
      read (ifileno,*)
      do while(i_err.eq.0)
         ! Exchange
         read (ifileno,*,iostat=i_err) iatom,jatom,dum,dum,dum,dum,dum,ncoup_tmp
         nlistsize(iatom)=nlistsize(iatom)+1
         nlist(nlistsize(iatom),iatom)=jatom
         ncoup(nlistsize(iatom),iatom,1)=ncoup_tmp
      end do
      close(ifileno)

      if (do_ralloy==0) then
         do i=1,nHam
            do j=1, nlistsize(i)
               jatom=nlist(j,i)
               if(abs(ammom_inp(anumb(i),1,1)*ammom_inp(anumb(jatom),1,1))<1e-6) then
                  ncoup(j,i,1)=0.0_dblprec
               else
                  ncoup(j,i,1)=ncoup(j,i,1)*&
                  fc2/ammom_inp(anumb(i),1,1)/ammom_inp(anumb(jatom),1,1)
               endif
            enddo
         enddo
      else
         do i=1,nHam
            do j=1, nlistsize(i)
               jatom = nlist(j,i)
               if(abs(ammom_inp(asite_ch(i),achem_ch(i),1)*ammom_inp(asite_ch(jatom),achem_ch(jatom),1))<1e-6) then
                  ncoup(j,i,1)=0.0_dblprec
               else
                  ncoup(j,i,1)=ncoup(j,i,1)* &
                  fc2/ammom_inp(asite_ch(i),achem_ch(i),1)/ammom_inp(asite_ch(jatom),achem_ch(jatom),1)
               endif
            enddo
         enddo
      endif

      10001 format (5x,i8,13x,i7)
      10002 format (13x,5i6)
      10003 format (5es16.8)

   end subroutine read_exchange

   !----------------------------------------------------------------------------
   !> Calculate a mean field estimate of the critical temperature
   !----------------------------------------------------------------------------
   subroutine estimate_tc_mfa(Natom,nHam,Nchmax,conf_num,NA,atype,anumb,Natom_full, &
         atype_ch,asite_ch,achem_ch,nlistsize, nlist,max_no_neigh,ncoup,ammom_inp,  &
         mconf,do_ralloy)
      !
      use Constants
      !use WangLandau, only : wl_emin, wl_estep, wl_nhist, wl_maxT
      !
      implicit none
      !
      integer, intent(in) :: NA           !< Number of atoms in one cell
      integer, intent(in) :: nHam         !< Number of atoms in Hamiltonian
      integer, intent(in) :: Natom        !< Number of atoms in system
      integer, intent(in) :: Nchmax       !< Number of chemical components in system
      integer, intent(in) :: conf_num     !< Number of configurations for LSF
      integer, intent(in) :: Natom_full   !< Number of atoms for full system (=Natom if not dilute)
      integer, intent(in) :: max_no_neigh !< Calculated maximum of neighbours for exchange
      integer, dimension(Natom), intent(in) :: atype !< Type of atom
      integer, dimension(Natom), intent(in) :: anumb !< Atom number in cell
      integer, dimension(Natom_full), intent(in) :: atype_ch !< Actual site of atom for dilute system
      integer, dimension(Natom_full), intent(in) :: asite_ch !< Actual site of atom for dilute system
      integer, dimension(Natom_full), intent(in) :: achem_ch !< Chemical type of atoms (reduced list)
      integer, dimension(nHam), intent(in) :: nlistsize !< Size of neighbour list for Heisenberg exchange couplings
      integer, dimension(max_no_neigh,Natom), intent(in) :: nlist !< Neighbour list for Heisenberg exchange couplings
      real(dblprec), dimension(max_no_neigh,nHam,conf_num), intent(in) :: ncoup !< Heisenberg exchange couplings
      real(dblprec), dimension(NA,Nchmax,conf_num), intent(in) :: ammom_inp !< Magnetic moment magnitudes from input
      integer, intent(in) :: mconf !< LSF ground state configuration
      integer, intent(in) :: do_ralloy                                  !< Random alloy simulation (0/1)
      !
      real(dblprec), dimension(:,:), allocatable :: jmat
      real(dblprec), dimension(:), allocatable :: rwork
      real(dblprec), dimension(:), allocatable :: work
      real(dblprec), dimension(:), allocatable :: cwres
      real(dblprec), dimension(:),allocatable :: ctemp
      real(dblprec), dimension(:,:),allocatable :: etemp
      real(dblprec) :: em,emin,fc2
      integer :: i,j,k, i_stat, i_all,ia,iatom
      integer :: lwork, info
      real(dblprec),dimension(nchmax) :: xconc
      !
      ! Allocate J0 matrix
      allocate(jmat(na,na),stat=i_stat)
      call memocc(i_stat,product(shape(jmat))*kind(jmat),'jmat','estimate_tc_mfa')
      jmat=0.0_dblprec

      lwork=3*na
      allocate(work(lwork))
      allocate(ctemp(na))
      allocate(etemp(na,na))
      allocate(rwork(lwork))
      allocate(cwres(NA))

      fc2 = 2.0_dblprec*mry/mub
      ! Fill J0 matrix

      if (do_ralloy==0) then
         do i=1,na
            do j=1,nlistsize(i)
               k=atype(nlist(j,i))
               jmat(i,k)=jmat(i,k)+ncoup(j,i,mconf)*ammom_inp(anumb(i),1,mconf)*ammom_inp(anumb(nlist(j,i)),1,mconf)
            end do
         end do
      else
         do ia=1,natom
            i=asite_ch(ia)
            do j=1,nlistsize(ia)
               k=asite_ch(nlist(j,ia))
               jmat(i,k)=jmat(i,k)+ncoup(j,ia,mconf)*ammom_inp(asite_ch(ia),achem_ch(ia),mconf)*ammom_inp(asite_ch(nlist(j,ia)),achem_ch(nlist(j,ia)),mconf)
            end do
         end do
         jmat=jmat/natom
      endif

      ! Scale matrix to eV energy
      jmat=jmat/fc2*ry_ev

      ! Find maximum eigen value
      call dgeev('N','N',NA, jmat, NA, cwres, ctemp, etemp, NA, etemp,NA,WORK, LWORK, INFO)
      if(info.ne.0) then
         print '(2x,a,i4)', 'Problem in dgeev:',info
      end if
      em=maxval(cwres)
      emin=minval(cwres)

      ! Wang-Landau entity extraction
      !wl_emin=-1.05_dblprec*maxval(abs(cwres))*natom/ry_ev
      !wl_nhist=2*int(abs(wl_emin))
      !wl_estep=abs(2.0_dblprec*wl_emin)/wl_nhist
      !wl_maxT=3.0_dblprec*(2.0_dblprec*em/3.0_dblprec/k_bolt_ev/1000)

      ! Write result to stdout
      write(*,'(f7.1,a)') (2.0_dblprec*em/3.0_dblprec/k_bolt_ev/1000),' K.'
      write(*,'(2x,a,f10.3,a)') 'Max eigenvalue:',em,' meV'
      if(na>1) write(*,'(2x,a,f10.3,a)') 'Min eigenvalue:',emin,' meV'
!      print *,'WL debug:',wl_emin,wl_estep

      ! Deallocate array
      i_all=-product(shape(jmat))*kind(jmat)
      deallocate(jmat,stat=i_stat)
      call memocc(i_stat,i_all,'jmat','estimate_tc_mfa')


      deallocate(work)
      deallocate(ctemp)
      deallocate(rwork)
      deallocate(cwres)
      deallocate(etemp)

! Random alloy chemical supermatrix
      if (do_ralloy==1) then

! Calculate actual concentration of the impurities in the supercell
         xconc=0.0_dblprec
         do ia=1,natom
           xconc(achem_ch(ia))=xconc(achem_ch(ia))+1.0_dblprec/natom
         enddo

         ! Allocate J0 matrix
         allocate(jmat(nchmax,nchmax),stat=i_stat)
         call memocc(i_stat,product(shape(jmat))*kind(jmat),'jmat','estimate_tc_mfa')
         jmat=0.0_dblprec

         lwork=3*nchmax
         allocate(work(lwork))
         allocate(ctemp(nchmax))
         allocate(etemp(nchmax,nchmax))
         allocate(rwork(lwork))
         allocate(cwres(nchmax))

         do ia=1,natom
            i=achem_ch(ia)
            do j=1,nlistsize(ia)
               k=achem_ch(nlist(j,ia))
               jmat(i,k)=jmat(i,k)+ncoup(j,ia,mconf)*ammom_inp(asite_ch(ia),achem_ch(ia),mconf)*ammom_inp(asite_ch(nlist(j,ia)),k,mconf)/xconc(i)
            end do
         end do
         jmat=jmat/natom
         ! Scale matrix to eV energy
         jmat=jmat/fc2*ry_ev

         ! Find maximum eigen value
         call dgeev('N','N',Nchmax,jmat,Nchmax,cwres,ctemp,etemp,Nchmax,etemp,Nchmax,WORK,LWORK,INFO)
         if(info.ne.0) then
            print '(2x,a,i4)', 'Problem in dgeev:',info
         end if
         em=maxval(cwres)
         emin=minval(cwres)
         ! Write result to stdout
         write(*,'(2x,a)',advance='no') "Calculate R-mean-field estimate of Tc: "
         write(*,'(f7.1,a)') (2.0_dblprec*em/3.0_dblprec/k_bolt_ev/1000),' K.'
         write(*,'(2x,a,f10.3,a)') 'Max eigenvalue:',em,' meV'
         if(na>1) write(*,'(2x,a,f10.3,a)') 'Min eigenvalue:',emin,' meV'
         ! Deallocate array
         i_all=-product(shape(jmat))*kind(jmat)
         deallocate(jmat,stat=i_stat)
         call memocc(i_stat,i_all,'jmat','estimate_tc_mfa')
         deallocate(work)
         deallocate(ctemp)
         deallocate(rwork)
         deallocate(cwres)
         deallocate(etemp)
      endif

   end subroutine estimate_tc_mfa

   !> Deallocate neighbour map array
   subroutine deallocate_nm()

      implicit none

      integer :: i_stat, i_all

      i_all=-product(shape(nm))*kind(nm)
      deallocate(nm,stat=i_stat)
      call memocc(i_stat,i_all,'nm','deallocat_nm')
      i_all=-product(shape(nmdim))*kind(nmdim)
      deallocate(nmdim,stat=i_stat)
      call memocc(i_stat,i_all,'nmdim','deallocat_nm')

   end subroutine deallocate_nm

   !!!!> Copies ncoup data to reduced matrix
   !!!subroutine setup_reduced_hamiltonian(Natom,NA,conf_num)

   !!!   use HamiltonianData, only : ham
   !!!   !.. Implicit declarations
   !!!   implicit none

   !!!   integer, intent(in) :: Natom   !< Number of atoms in system
   !!!   integer, intent(in) :: NA         !< Number of atoms in one cell
   !!!   integer, intent(in) :: conf_num !< Number of configurations for LSF

   !!!   integer :: i,j

   !!!   do i=1,NA
   !!!      ham%nlistsize_red(i)=ham%nlistsize(i)
   !!!      do j=1,ham%max_no_neigh
   !!!         ham%ncoup_red(j,i,:)=ham%ncoup(j,i,:)
   !!!      end do
   !!!   end do

   !!!end subroutine setup_reduced_hamiltonian


end module HamiltonianInit
