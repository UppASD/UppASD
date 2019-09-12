!-------------------------------------------------------------------------------
!  MODULE: clusters
!> @brief
!> Intended to be able to embed clusters of different atoms inside a
!> host system.
!> @details The idea behind this is to give the user more flexibility
!> when treating systems that lack symmetry and that can be described
!> as an inpurity cluster embeded inside a host
!
!> @author
!> Jonathan Chico
!> @copyright
!> GNU Public License.
!> @note
!> Since this routine breaks all symmetries by construction, The functionality of
!> `do_reduced` is not implemented at all here.
!-------------------------------------------------------------------------------
module clusters
   use Parameters
   use Profiling

   implicit none

   integer :: N1_clus   !< Number of cell repetitions in x direction for the cluster
   integer :: N2_clus   !< Number of cell repetitions in y direction for the cluster
   integer :: N3_clus   !< Number of cell repetitions in z direction for the cluster
   integer :: NA_clus   !< Number of atoms in the cluster unit cell
   integer :: NT_clus   !< Number of types of atoms in the cluster unit cell
   integer :: Natom_clus      !< Number of atoms in the cluser
   integer :: Nchmax_clus     !< Max number of chemical components on each site in cell for the cluster
   integer :: Natom_full_clus !< Number of atoms for the cluster (=Natom if not dilute)
   integer :: do_anisotropy_clus     !< Read anisotropy data for the cluster (1/0)
   integer :: max_no_shells_clus     !< Actual maximum number of shells for the cluster
   integer :: max_no_dmshells_clus   !< Actual maximum number of shells for DM interactions for the cluster
   integer, dimension(:), allocatable :: NN_clus                        !< Number of neighbour shells for the cluster
   integer, dimension(:), allocatable :: Nch_clus                       !< Number of chemical components on each site in cell
   integer, dimension(:), allocatable :: dm_nn_clus                     !< No. shells of neighbours for DM for the cluster
   integer, dimension(:), allocatable :: anumb_inp_clus                 !< Atom number in cluster
   integer, dimension(:), allocatable :: atype_inp_clus                 !< Type of atom in cluster
   integer, dimension(:,:), allocatable :: anisotropytype_clus          !< Type of anisotropies for the cluster (0-2)
   integer, dimension(:,:), allocatable :: anisotropytype_diff_clus     !< Type of anisotropies when one is considering more than one anisotropy axis for the cluster
   real(dblprec) :: random_anisotropy_density_clus                      !< Density for randomness in anisotropy for the cluster
   real(dblprec), dimension(3) :: C1_clus                               !< First lattice vector for the cluster
   real(dblprec), dimension(3) :: C2_clus                               !< Second lattice vector for the cluster
   real(dblprec), dimension(3) :: C3_clus                               !< Third lattice vector for the cluster
   real(dblprec), dimension(:,:), allocatable :: Bas_clus               !< Coordinates for basis atoms for the cluster
   real(dblprec), dimension(:,:), allocatable :: chconc_clus            !< Chemical concentration on sites
   real(dblprec), dimension(:,:,:), allocatable :: redcoord_clus        !< Coordinates for Heisenberg exchange couplings for the cluster
   real(dblprec), dimension(:,:,:), allocatable :: ammom_inp_clus       !< Magnetic moment magnitudes from input for the cluster
   real(dblprec), dimension(:,:,:), allocatable :: Landeg_ch_clus       !< Gyromagnetic ratio for the cluster
   real(dblprec), dimension(:,:,:), allocatable :: anisotropy_clus      !< Input data for anisotropies for the cluster
   real(dblprec), dimension(:,:,:), allocatable :: dm_redcoord_clus     !< Neighbour vectors for DM for the cluster
   real(dblprec), dimension(:,:,:), allocatable :: anisotropy_diff_clus !< Input data for the second anisotropy axis for the cluster
   real(dblprec), dimension(:,:,:,:), allocatable :: aemom_inp_clus     !< Magnetic moment directions from input for the cluster
   real(dblprec), dimension(:,:,:,:,:), allocatable :: jc_clus          !< Exchange couplings for the cluster
   real(dblprec), dimension(:,:,:,:,:), allocatable :: dm_inpvect_clus  !< Neighbour vectors for DM for the cluster
   logical :: random_anisotropy_clus   !< Put random anisotropy in the sample for the cluster (T/F)
   character(len=1) :: do_cluster      !< Perform cluster embedding procedure
   character(len=1) :: mult_axis_clus  !< Flag to treat more than one anisotropy axis at the same time for the cluster
   character(len=35) :: kfile_clus     !< File name for the anisotropies for the cluster
   character(len=35) :: dmfile_clus    !< File name for Dzyaloshinskii-Moriya data of the cluster
   character(len=35) :: posfile_clus   !< File name for coordinates of the cluster
   character(len=35) :: momfile_clus   !< File name for magnetic moments for the cluster
   character(len=30), dimension(:), allocatable ::jfile_clus !< File name for exchange couplings for the cluster

   type ham_clus_t
      integer :: max_no_neigh_clus      !< Calculated maximum of neighbours for exchange for the cluster
      integer :: max_no_dmneigh_clus    !< Calculated maximum of neighbours for DM exchange for the cluster
      integer, dimension(:), allocatable :: taniso_clus        !< Type of anisotropy for cluster (0-2)
      integer, dimension(:), allocatable :: nlistsize_clus     !< Size of neighbour list for Heisenberg exchange couplings for the cluster
      integer, dimension(:), allocatable :: dmlistsize_clus    !< Size of neighbour list for DM for the cluster
      integer, dimension(:), allocatable :: taniso_diff_clus   !< Type of anisotropy for cluster (0-2)
      integer, dimension(:,:), allocatable :: nlist_clus       !< Neighbour list for Heisenberg exchange couplings for the cluster
      integer, dimension(:,:), allocatable :: dmlist_clus      !< List of neighbours for DM for the cluster
      integer, dimension(:,:), allocatable :: ind_mom_clus     !< Flag to decide whether a moment is induced or not
      real(dblprec), dimension(:), allocatable :: sb_clus            !< Ratio between Cubic and Uniaxial anisotropy for the cluster
      real(dblprec), dimension(:), allocatable :: sb_diff_clus       !< Ratio between Cubic and Uniaxial anisotropy for the cluster
      real(dblprec), dimension(:,:), allocatable :: kaniso_clus      !< Anisotropy constant for the cluster
      real(dblprec), dimension(:,:), allocatable :: eaniso_clus      !< Unit anisotropy vector for the cluster
      real(dblprec), dimension(:,:), allocatable :: kaniso_diff_clus !< Anisotropy constant for the cluster
      real(dblprec), dimension(:,:), allocatable :: eaniso_diff_clus !< Unit anisotropy vector for the cluster
      real(dblprec), dimension(:,:,:), allocatable :: ncoup_clus     !< Heisenberg exchange couplings for the cluster
      real(dblprec), dimension(:,:,:), allocatable :: dm_vect_clus   !< Dzyaloshinskii-Moriya exchange vector
   end type ham_clus_t

   type(ham_clus_t) :: ham_clus

contains

   !----------------------------------------------------------------------------
   ! SUBROUTINE: reading_wrapper_clus
   !> @brief Subroutine to read some of the data needed for the cluster method
   !> @author Jonathan Chico
   !----------------------------------------------------------------------------
   subroutine reading_wrapper_clus(conf_num,do_ralloy,set_landeg,do_anisotropy_clus,&
      Landeg_global,do_lsf,posfiletype,ind_mom_flag)

      implicit none

      integer, intent(in) :: conf_num     !< Number of configurations for LSF
      integer, intent(in) :: do_ralloy    !< Random alloy simulation (0/1)
      integer, intent(in) :: set_landeg   !< Set 'custom' value for gyromagnetic ration (1=yes,0=no)
      integer, intent(in) :: do_anisotropy_clus    !< Read anisotropy data for the cluster (1/0)
      real(dblprec), intent(in) :: Landeg_global   !< Default gyromagnetic ratio
      character(len=1), intent(in) :: do_lsf       !< (Y/N) Do LSF for MC
      character(len=1), intent(in) :: posfiletype  !< posfile type (D)irect or (C)arteisian coordinates in posfile
      character(len=1), intent(in) :: ind_mom_flag !< Flag to indicate that there are induced moments being considered

      if (do_ralloy==0) then
         call read_positions_clus(posfiletype)
         if (do_anisotropy_clus==1) then
            call allocate_anisotropy_input_clus(1)
            call read_anisotropy_clus()
         endif
      else
         call read_positions_clus_alloy(posfiletype)
         if (do_anisotropy_clus==1) then
            call allocate_anisotropy_input_clus(1)
            call read_anisotropy_clus()
         endif
      endif
      Natom_full_clus=NA_clus*N1_clus*N2_clus*N3_clus
      call read_moments_clus(Landeg_global,conf_num,ind_mom_flag,do_lsf,set_landeg)

   end subroutine reading_wrapper_clus

   !-----------------------------------------------------------------------------
   !  SUBROUTINE: cluster_creation
   !> @brief
   !> Wrapper routine that contains the majority of the calls needed to embed an
   !> impurity cluster inside a host in ASD.
   !> @details Of this way one can embed systems with a
   !> large number of different types inside a host with a given exchange.
   !> Currently one can take into account, different moments, Heisenberg exchange,
   !> DMI and anisotropies in the cluster.
   !
   !> @author
   !> Jonathan Chico
   !-----------------------------------------------------------------------------
   subroutine cluster_creation(NT,do_dm,Natom,Nchmax,initmag,conf_num,Mensemble,    &
      do_ralloy,Natom_full,do_jtensor,do_prnstruct,do_anisotropy_clus,index_clus,   &
      atype_clus,anumb_clus,coord,coord_clus,simid,mult_axis_clus,atype,anumb,      &
      atype_ch,asite_ch,achem_ch,achtype,acellnumb,mmom,emom,emomM,Landeg)

      use Profiling
      use Parameters
      use HamiltonianData,    only : ham
      use PrintHamiltonian,   only : prn_exchange, prn_dmcoup

      implicit none

      integer, intent(in) :: NT              !< Number of types of atoms
      integer, intent(in) :: do_dm           !< Add Dzyaloshinskii-Moriya (DM) term to Hamiltonian (0/1)
      integer, intent(in) :: Natom           !< Number of atoms in system
      integer, intent(in) :: Nchmax          !< Max number of chemical components on each site in cell
      integer, intent(in) :: initmag         !< Mode of initialization of magnetic moments (1-4)
      integer, intent(in) :: conf_num        !< Number of configurations for LSF
      integer, intent(in) :: Mensemble       !< Number of ensembles
      integer, intent(in) :: do_ralloy       !< Random alloy simulation (0/1)
      integer, intent(in) :: Natom_full      !< Number of atoms for full system (=Natom if not dilute)
      integer, intent(in) :: do_jtensor      !< Use SKKR style exchange tensor (0=off, 1=on)
      integer, intent(in) :: do_prnstruct    !< Print Hamiltonian information (0/1)
      integer, intent(in) :: do_anisotropy_clus !< Read anisotropy data for the cluster (1/0)
      integer, dimension(Natom_clus), intent(in)         :: index_clus  !< Mapping of cluster indices to host indices
      integer, dimension(Natom_clus), intent(in)         :: atype_clus  !< Type of atom for the cluster
      integer, dimension(Natom_clus), intent(in)         :: anumb_clus  !< Atom number in cell for the cluster
      real(dblprec), dimension(3,Natom), intent(in)      :: coord
      real(dblprec), dimension(3,Natom_clus), intent(in) :: coord_clus  !< Coordinates of all atoms belonging to the cluster
      character(len=8), intent(in) :: simid           !< Name of simulation
      character(len=1), intent(in) :: mult_axis_clus  !< Flag to treat more than one anisotropy axis at the same time
      ! In/out variables
      integer, dimension(Natom), intent(inout)        :: atype       !< Type of atom
      integer, dimension(Natom), intent(inout)        :: anumb       !< Atom number in cell
      integer, dimension(Natom), intent(inout)        :: achtype     !< Chemical type of atoms (full list)
      integer, dimension(Natom_full), intent(inout)   :: atype_ch    !< Actual type of atom for dilute system
      integer, dimension(Natom_full), intent(inout)   :: asite_ch    !< Actual site of atom for dilute system
      integer, dimension(Natom_full), intent(inout)   :: achem_ch    !< Chemical type of atoms (reduced list)
      integer, dimension(Natom), intent(inout)        :: acellnumb   !< List for translating atom no. in full cell to actual cell
      real(dblprec), dimension(Natom), intent(inout) :: Landeg  !< Gyromagnetic ratio
      real(dblprec), dimension(Natom,Mensemble), intent(inout) :: mmom     !< Magnitude of magnetic moments
      real(dblprec), dimension(3,Natom,Mensemble), intent(inout) :: emom   !< Current unit moment vector
      real(dblprec), dimension(3,Natom,Mensemble), intent(inout) :: emomM  !< Current magnetic moment vector
      ! .. Local variables
      integer :: iatom,iclus,jclus,jatom,neigh,ineigh_c,ens
      integer :: list_size
      integer, dimension(ham%max_no_neigh) :: temp_nlist
      character(len=30) :: filn

      ! In this loop one would identify which atoms of the host, would belong to the cluster
      !$omp parallel do default(shared) private(iatom,iclus,ens)
      do iclus=1, Natom_clus
         iatom=index_clus(iclus)
         ! Now that one knows which atoms belong to the cluster one can change their types
         atype(iatom)=atype_clus(iclus)+NT
         anumb(iatom)=iatom
         ! One can also change the magnitude and direction of the moments
         ! If one chooses to have the directions as given from the momfile then one can
         ! set them by the values goten from the momfile of the cluster
         Landeg(iatom)=Landeg_ch_clus(iclus,1,1)
         if (initmag.eq.3) then
            do ens=1, Mensemble
               emom(1,iatom,ens)=aemom_inp_clus(1,anumb_clus(iclus),1,1)
               emom(2,iatom,ens)=aemom_inp_clus(2,anumb_clus(iclus),1,1)
               emom(3,iatom,ens)=aemom_inp_clus(3,anumb_clus(iclus),1,1)
               mmom(iatom,ens)=ammom_inp_clus(anumb_clus(iclus),1,1)
               emomM(1,iatom,ens)=emom(1,iatom,ens)*mmom(iatom,ens)
               emomM(2,iatom,ens)=emom(2,iatom,ens)*mmom(iatom,ens)
               emomM(3,iatom,ens)=emom(3,iatom,ens)*mmom(iatom,ens)
            enddo
            ! Otherwise the directions are unchanged
         else
            do ens=1, Mensemble
               mmom(iatom,ens)=ammom_inp_clus(atype_clus(iclus),1,1)
               emomM(1,iatom,ens)=emom(1,iatom,ens)*mmom(iatom,ens)
               emomM(2,iatom,ens)=emom(2,iatom,ens)*mmom(iatom,ens)
               emomM(3,iatom,ens)=emom(3,iatom,ens)*mmom(iatom,ens)
            enddo
         endif
      enddo
      !$omp end parallel do
      ! Now one can overwrite the exchange interactions for the atoms belonging to the cluster
      do iatom=1,Natom
         ! Check if the current atom is on the cluster list
         if (ANY(index_clus==iatom)) then
            ! Find the index for which this happens
            iclus = minloc(abs(index_clus - iatom), 1)
            list_size=max(ham%nlistsize(iatom),ham_clus%nlistsize_clus(iclus))
            ham%nlistsize(iatom)=list_size
            ! Loop over the neighbor atoms in the cluster
            do ineigh_c=1, ham_clus%nlistsize_clus(iclus)
               jclus=ham_clus%nlist_clus(ineigh_c,iclus)
               jatom=index_clus(jclus)
               ! Find the index where the neighbour list is equal to the present item
               temp_nlist=ham%nlist(:,iatom)
               neigh=0
               if (ANY(ham%nlist(:,iatom)==jatom)) then
                  neigh=minloc(abs(temp_nlist-jatom),1)
                  ! Overwrite the interaction of the host with the one of the cluster
                  ham%nlist(neigh,iatom)=jatom
                  ham%ncoup(neigh,iatom,1)=ham_clus%ncoup_clus(ineigh_c,iclus,1)
               else
                  if (minval(temp_nlist(:ham%nlistsize(iatom)))==0) then
                     neigh=minloc(temp_nlist,1)
                     ham%nlist(neigh,iatom)=jatom
                     ham%ncoup(neigh,iatom,1)=ham_clus%ncoup_clus(ineigh_c,iclus,1)
                  else
                     ham%nlistsize(iatom)=ham%nlistsize(iatom)+1
                     ham%nlist(ham%nlistsize(iatom),iatom)=jatom
                     ham%ncoup(ham%nlistsize(iatom),iatom,1)=ham_clus%ncoup_clus(ineigh_c,iclus,1)
                  endif
               endif
            enddo
         endif
      enddo

      if (do_dm.eq.1) then
         ! Now after the exchange has been mapped one can map the DMI of the cluster
         do iatom=1,Natom
            ! Check if the current atom is on the cluster list
            if (ANY(index_clus==iatom)) then
               ! Find the index for which this happens
               iclus = minloc(abs(index_clus - iatom), 1)
               list_size=max(ham%dmlistsize(iatom),ham_clus%dmlistsize_clus(iclus))
               ham%dmlistsize(iatom)=list_size
               ! Loop over the neighbor atoms in the cluster
               do ineigh_c=1, ham_clus%dmlistsize_clus(iclus)
                  jclus=ham_clus%dmlist_clus(ineigh_c,iclus)
                  jatom=index_clus(jclus)
                  ! Find the index where the neighbour list is equal to the present item
                  temp_nlist=ham%dmlist(:,iatom)
                  neigh=0
                  if (ANY(ham%dmlist(:,iatom)==jatom)) then
                     neigh=minloc(abs(temp_nlist-jatom),1)
                     ! Overwrite the interaction of the host with the one of the cluster
                     ham%dmlist(neigh,iatom)=jatom
                     ham%dm_vect(:,neigh,iatom)=ham_clus%dm_vect_clus(:,ineigh_c,iclus)
                  else
                     if (minval(temp_nlist(:ham%dmlistsize(iatom)))==0) then
                        neigh=minloc(temp_nlist,1)
                        ham%dmlist(neigh,iatom)=jatom
                        ham%dm_vect(:,neigh,iatom)=ham_clus%dm_vect_clus(:,ineigh_c,iclus)
                     else
                        ham%dmlistsize(iatom)=ham%dmlistsize(iatom)+1
                        ham%dmlist(ham%dmlistsize(iatom),iatom)=jatom
                        ham%dm_vect(:,ham%dmlistsize(iatom),iatom)=ham_clus%dm_vect_clus(:,ineigh_c,iclus)
                     endif
                  endif
               enddo
            endif
         enddo
      endif

      if (do_anisotropy_clus.eq.1) then
         !If one has an anisotopy then one overwrites the anisotopy from the bulk
         ! with the one of the cluster
         !$omp parallel do default(shared) private(iclus,iatom)
         do iclus=1,Natom_clus
            iatom=index_clus(iclus)
            ham%sb(iatom)           = ham_clus%sb_clus(iclus)
            ham%taniso(iatom)       = ham_clus%taniso_clus(iclus)
            ham%eaniso(1:3,iatom)   = ham_clus%eaniso_clus(1:3,iclus)
            ham%kaniso(1:2,iatom)   = ham_clus%kaniso_clus(1:2,iclus)
         enddo
         !$omp end parallel do

         if (mult_axis_clus=='Y') then
            !$omp parallel do default(shared) private(iclus,iatom)
            do iclus=1,Natom_clus
               iatom=index_clus(iclus)
               ham%sb_diff(iatom)         = ham_clus%sb_diff_clus(iclus)
               ham%taniso_diff(iatom)     = ham_clus%taniso_diff_clus(iclus)
               ham%eaniso_diff(1:3,iatom) = ham_clus%eaniso_diff_clus(1:3,iclus)
               ham%kaniso_diff(1:2,iatom) = ham_clus%kaniso_diff_clus(1:2,iclus)
            end do
            !$omp end parallel do
         endif
      endif
      ! Impose symmetrization of the interactions to avoid probles at the boundaries
      call symmetrize_interactions(Natom,do_dm)
      ! Write extra information for the cluster
      write (filn,'(''clus_info.'',a8,''.out'')') simid
      open(ofileno, file=filn)
      write(ofileno,'(a)') "Clus. ind.  iatom  coord_x  coord_y  coord_z  type"
      do iclus=1,Natom_clus
         write(ofileno,'(i6,2x,i6,2x,3f12.6,i6)') iclus,index_clus(iclus),coord_clus(1:3,iclus),atype_clus(iclus)
      enddo
      close(ofileno)

      if (do_prnstruct==1) then
         write(*,'(2x,a)',advance='no') "Print exchange interaction strenghts after embedding"
         call prn_exchange(Natom,1,Natom,1,do_ralloy,Natom_full,ham%max_no_neigh,   &
            simid,anumb,atype,ham%nlistsize,asite_ch,achem_ch,ham%nlist,coord,mmom, &
            ham%ncoup,ham%aham)
         write(*,'(a)') ' done'

         if (do_dm==1) then
            write(*,'(2x,a)',advance='no') "Print Dzyaloshinskii-Moriya interactions after embedding"
            call prn_dmcoup(Natom,Natom,1,do_ralloy,Natom_full,ham%max_no_dmneigh,  &
               anumb,atype,ham%dmlistsize,asite_ch,achem_ch,ham%dmlist,coord,       &
               mmom,ham%dm_vect,simid)
            write(*,'(a)') ' done'
         endif
      endif
      ! Deallocate unnecessary arrays
      call allocate_cluster_hamiltoniandata(do_jtensor,NA_clus,conf_num,            &
         ham_clus%max_no_neigh_clus,-1)
      call allocate_hamiltonianinput_clus(Nchmax_clus,conf_num,flag=-1)
      if (do_dm==1) then
         call allocate_cluster_dmhamiltoniandata(NA_clus,                           &
            ham_clus%max_no_dmneigh_clus,-1)
      endif

      if (do_anisotropy_clus==1) then
         call allocate_cluster_anisotropies(NA_clus,mult_axis_clus,-1)
         call allocate_anisotropy_input_clus(-1)
      endif

   end subroutine cluster_creation

   !----------------------------------------------------------------------------
   ! SUBROUTINE: symmetrize_interactions
   !> @brief Makes sure that \f$ J_{ij}=J_{ji}\f$ and that \f$|D_{ij}|=|D_{ji}|f$
   !> @details Due to the way the embedding is performed there might be problems
   !> at the boundaries, this routine tries to ensure that the behaviour of the
   !> interactions at the boundaries corresponds to what one expects from the Heisenberg
   !> model.
   !> @author Jonathan Chico
   !----------------------------------------------------------------------------
   subroutine symmetrize_interactions(Natom,do_dm)

      use HamiltonianData,    only : ham

      implicit none

      integer, intent(in) :: Natom
      integer, intent(in) :: do_dm

      integer :: iatom, ineigh,jatom,jneigh
      real(dblprec) :: temp_Jij,temp_Jji, temp_Dij,temp_Dji,tol

      tol=1e-5

      do iatom=1,Natom
         ! Loop over the neighbours to iatom
         do ineigh=1,ham%nlistsize(iatom)
            ! Find the neighbouring atoms to iatom
            jatom=ham%nlist(ineigh,iatom)
            ! Loop over the neighbours of the jatom
            do jneigh=1,ham%nlistsize(jatom)
               ! If the neighbor to jatom is iatom one must symmetrize their interactions
               if (ham%nlist(jneigh,jatom)==iatom) then
                  temp_Jij=ham%ncoup(ineigh,iatom,1)
                  temp_Jji=ham%ncoup(jneigh,jatom,1)
                  ! Average over the interactions
                  ham%ncoup(ineigh,iatom,1)=(temp_Jij+temp_Jji)*0.5_dblprec
                  ham%ncoup(jneigh,jatom,1)=(temp_Jij+temp_Jji)*0.5_dblprec
               endif
            enddo
         enddo
         ! If DMI is present symmetrize those interactions as well
         if (do_dm==1) then
            ! Loop over the neighbours to iatom
            do ineigh=1,ham%dmlistsize(iatom)
               ! Find the neighbouring atoms to iatom
               jatom=ham%dmlist(ineigh,iatom)
               ! Loop over the neighbours of the jatom
               do jneigh=1,ham%dmlistsize(jatom)
                  ! If the neighbor to jatom is iatom one must symmetrize their interactions
                  if (ham%dmlist(jneigh,jatom)==iatom) then
                     temp_Dij=norm2(ham%dm_vect(:,ineigh,iatom))
                     temp_Dji=norm2(ham%dm_vect(:,jneigh,jatom))
                     if (temp_Dij<tol) then
                        temp_Dij=1.0_dblprec
                     endif
                     if (temp_Dji<tol) then
                        temp_Dji=1.0_dblprec
                     endif
                     ! Average over the interactions
                     ham%dm_vect(:,ineigh,iatom)=ham%dm_vect(:,ineigh,iatom)*(temp_Dij+temp_Dji)*0.5_dblprec/(temp_Dij)
                     ham%dm_vect(:,jneigh,jatom)=ham%dm_vect(:,jneigh,jatom)*(temp_Dij+temp_Dji)*0.5_dblprec/(temp_Dji)
                  endif
               enddo
            enddo
         endif
      enddo

   end subroutine symmetrize_interactions

   !----------------------------------------------------------------------------
   ! SUBROUTINE: set_input_defaults_clus
   !> @brief Set defaults for the input variables for the clusters application
   !> @author Jonathan Chico
   !----------------------------------------------------------------------------
   subroutine set_input_defaults_clus()

      implicit none

      real(dblprec) :: one = 1.0_dblprec
      real(dblprec) :: zero= 0.0_dblprec

      NA_clus              = 0
      NT_clus              = 0
      N1_clus              = 0
      N2_clus              = 0
      N3_clus              = 0
      Natom_clus           = 0
      Nchmax_clus          = 1
      do_anisotropy_clus   = 0
      do_cluster           = 'N'
      kfile_clus           = 'kfile_clus'
      dmfile_clus          = 'dmfile_clus'
      momfile_clus         = 'momfile_clus'
      posfile_clus         = 'posfile_clus'
      mult_axis_clus       = 'N'
      C1_clus              = (/one,zero,zero/)
      C2_clus              = (/zero,one,zero/)
      C3_clus              = (/zero,zero,one/)

   end subroutine set_input_defaults_clus


   !---------------------------------------------------------------------------
   ! SUBROUTINE read_positions_clus
   !> @brief
   !> Read Positions for the cluster atoms
   !
   !> @author
   !> Jonathan Chico, based on the routine by Anders Bergman
   !---------------------------------------------------------------------------
   subroutine read_positions_clus(posfiletype)
      !
      implicit none

      character(len=1), intent(in) :: posfiletype  !< posfile type (D)irect or (C)arteisian coordinates in posfile
      !
      integer :: flines,itype,mtype,iat,isite,msite,i_stat
      real(dblprec),dimension(3) :: tmp

      ! Open input file
      open(ifileno, file=posfile_clus)
      ! Check if input file is for random alloy
      flines=0
      mtype=0
      msite=0
      ! Pre-read file to get max no. sites, types and chemical species
      do
         read(ifileno,*,end=200)  isite,itype
         flines=flines+1
         msite=max(msite,isite)
         mtype=max(mtype,itype)
      end do
      200 continue
      rewind(ifileno)

      ! Set no. sites, types and chemical species
      if(msite/=flines) write(*,*) "WARNING: Check input file ", posfile_clus,      &
         " for inconsistent information."
      NA_clus=msite
      NT_clus=mtype
      ! Allocate input arrays
      allocate(Bas_clus(3,NA_clus),stat=i_stat)
      call memocc(i_stat,product(shape(Bas_clus))*kind(Bas_clus),'Bas_clus','read_positions_clus')
      allocate(atype_inp_clus(NA_clus),stat=i_stat)
      call memocc(i_stat,product(shape(atype_inp_clus))*kind(atype_inp_clus),'atype_inp_clus','read_positions_clus')
      allocate(anumb_inp_clus(NA_clus),stat=i_stat)
      call memocc(i_stat,product(shape(anumb_inp_clus))*kind(anumb_inp_clus),'anumb_inp_clus','read_positions_clus')

      ! Read basis atoms and setup type array
      ! Site, Type, Rx, Ry, Rz
      if (posfiletype=='C') then
         do iat=1, NA_clus
            read (ifileno,*) isite,itype,Bas_clus(1,isite),Bas_clus(2,isite),       &
               Bas_clus(3,isite)
            atype_inp_clus(isite)=itype

            ! Redundant but kept for the time beeing
            anumb_inp_clus(isite)=isite
         enddo
      elseif (posfiletype=='D') then
         do iat=1, NA_clus
            read (ifileno,*) isite, itype,  tmp(1),tmp(2),tmp(3)
            atype_inp_clus(isite)=itype
            Bas_clus(1,isite)=tmp(1)*C1_clus(1)+tmp(2)*C2_clus(1)+tmp(3)*C3_clus(1)
            Bas_clus(2,isite)=tmp(1)*C1_clus(2)+tmp(2)*C2_clus(2)+tmp(3)*C3_clus(2)
            Bas_clus(3,isite)=tmp(1)*C1_clus(3)+tmp(2)*C2_clus(3)+tmp(3)*C3_clus(3)
            ! Redundant but kept for the time beeing
            anumb_inp_clus(isite)=isite
         enddo
      endif
      close (ifileno)
   end subroutine read_positions_clus

   !----------------------------------------------------------------------------
   !> @brief
   !> Read Positions for random alloys for the cluster
   !
   !> @author
   !> Jonathan Chico, routine based in the read_positions_alloy routine by Anders Bergman
   !----------------------------------------------------------------------------
   subroutine read_positions_clus_alloy(posfiletype)
      !
      !
      implicit none
      !
      character(len=1), intent(in) :: posfiletype  !< posfile type (D)irect or (C)arteisian coordinates in posfile
      !
      integer :: flines,isite,itype,ichem
      integer :: msite,mtype,mchem,iat, i_stat
      real(dblprec) :: rconc
      real(dblprec),dimension(3) :: tmp

      ! Open input file
      open(ifileno, file=posfile_clus)

      ! Check if input file is for random alloy
      flines=0
      msite=0
      mtype=0
      mchem=0

      ! Pre-read file to get max no. sites, types and chemical species
      do
         read(ifileno,*,end=200)  isite,itype,ichem
         flines=flines+1
         msite=max(msite,isite)
         mtype=max(mtype,itype)
         mchem=max(mchem,ichem)
      end do
      200 continue
      rewind(ifileno)

      ! Set no. sites, types and chemical species
      na_clus=msite
      nt_clus=mtype
      nchmax_clus=mchem

      ! Allocate input arrays
      allocate(bas_clus(3,msite),stat=i_stat)
      call memocc(i_stat,product(shape(bas_clus))*kind(bas_clus),'bas_clus','read_positions_clus_alloy')
      allocate(atype_inp_clus(NA_clus),stat=i_stat)
      call memocc(i_stat,product(shape(atype_inp_clus))*kind(atype_inp_clus),'atype_inp_clus','read_positions_clus_alloy')
      allocate(anumb_inp_clus(NA_clus),stat=i_stat)
      call memocc(i_stat,product(shape(anumb_inp_clus))*kind(anumb_inp_clus),'anumb_inp_clus','read_positions_clus_alloy')

      ! Chemical data
      allocate(nch_clus(na_clus),stat=i_stat)
      call memocc(i_stat,product(shape(nch_clus))*kind(nch_clus),'nch_clus','read_positions_clus_alloy')
      nch_clus=0
      allocate(chconc_clus(na_clus,nchmax_clus),stat=i_stat)
      call memocc(i_stat,product(shape(chconc_clus))*kind(chconc_clus),'chconc_clus','read_positions_clus_alloy')
      chconc_clus=0.0_dblprec

      ! Read data
      ! Site  Type   Chem_type Conc, Rx, Ry, Rz
      if (posfiletype=='C') then
         do iat=1, flines
            read (ifileno,*) isite, itype, ichem, rconc, bas_clus(1,isite), bas_clus(2,isite), bas_clus(3,isite)
            nch_clus(isite)=max(nch_clus(isite),ichem)
            atype_inp_clus(isite)=itype
            anumb_inp_clus(isite)=isite
            chconc_clus(isite,ichem)=rconc
         enddo
      elseif(posfiletype=='D') then
         do iat=1, flines
            read (ifileno,*) isite, itype,ichem,rconc,  tmp(1),tmp(2),tmp(3)
            bas_clus(1,isite)=tmp(1)*C1_clus(1)+tmp(2)*C2_clus(1)+tmp(3)*C3_clus(1)
            bas_clus(2,isite)=tmp(1)*C1_clus(2)+tmp(2)*C2_clus(2)+tmp(3)*C3_clus(2)
            bas_clus(3,isite)=tmp(1)*C1_clus(3)+tmp(2)*C2_clus(3)+tmp(3)*C3_clus(3)
            nch_clus(isite)=max(nch_clus(isite),ichem)
            atype_inp_clus(isite)=itype
            anumb_inp_clus(isite)=isite
            chconc_clus(isite,ichem)=rconc
         enddo
      endif
      close (ifileno)
   end subroutine read_positions_clus_alloy

   !-----------------------------------------------------------------------------
   ! SUBROUTINE: deallocate_cluster_info
   !> @brief deallocate cluster information
   !> @author Jonathan Chico
   !-----------------------------------------------------------------------------
   subroutine deallocate_cluster_info()

      use Profiling

      implicit none

      integer :: i_all, i_stat

      i_all=-product(shape(Bas_clus))*kind(Bas_clus)
      deallocate(Bas_clus,stat=i_stat)
      call memocc(i_stat,i_all,'Bas_clus','deallocate_cluster_info')

      i_all=-product(shape(ammom_inp_clus))*kind(ammom_inp_clus)
      deallocate(ammom_inp_clus,stat=i_stat)
      call memocc(i_stat,i_all,'ammom_inp_clus','deallocate_cluster_info')

      i_all=-product(shape(aemom_inp_clus))*kind(aemom_inp_clus)
      deallocate(aemom_inp_clus,stat=i_stat)
      call memocc(i_stat,i_all,'aemom_inp_clus','deallocate_cluster_info')

      i_all=-product(shape(Landeg_ch_clus))*kind(Landeg_ch_clus)
      deallocate(Landeg_ch_clus,stat=i_stat)
      call memocc(i_stat,i_all,'Landeg_ch_clus','deallocate_cluster_info')

      i_all=-product(shape(ham_clus%ind_mom_clus))*kind(ham_clus%ind_mom_clus)
      deallocate(ham_clus%ind_mom_clus,stat=i_stat)
      call memocc(i_stat,i_all,'ham_clus%ind_mom_clus','deallocate_cluster_info')


   end subroutine deallocate_cluster_info

   !---------------------------------------------------------------------------
   ! SUBROUTINE read_moments_clus
   !> @brief
   !> Read Magnetic moments for the cluster
   !
   !> @author
   !> Jonathan Chico based on the routines by Anders Bergman
   !---------------------------------------------------------------------------
   subroutine read_moments_clus(Landeg_global,conf_num,ind_mom_flag,do_lsf,set_landeg)

      use LSF
      use Profiling
      use Parameters

      implicit none

      integer, intent(in) :: set_landeg            !< Set 'custom' value for gyromagnetic ration (1=yes,0=no)
      integer, intent(in) :: conf_num              !< Number of configurations for LSF
      character(len=1), intent(in) :: ind_mom_flag !< Flag to indicate that there are induced moments being considered
      character(len=1), intent(in) :: do_lsf       !< (Y/N) Do LSF for MC
      real(dblprec), intent(in) :: Landeg_global   !< Default gyromagnetic ratio
      !
      integer :: i_err,isite,ichem,i_stat,iconf
      real(dblprec)  :: aemom_tmp

      iconf = 1

      open(ifileno,file=trim(momfile_clus))

      !Allocate arrays according to data from position input
      allocate(ammom_inp_clus(NA_clus,nchmax_clus,conf_num),stat=i_stat)
      call memocc(i_stat,product(shape(ammom_inp_clus))*kind(ammom_inp_clus),'ammom_inp_clus','read_moments_clus')
      ammom_inp_clus=0.0_dblprec

      allocate(aemom_inp_clus(3,NA_clus,nchmax_clus,conf_num),stat=i_stat)
      call memocc(i_stat,product(shape(aemom_inp_clus))*kind(aemom_inp_clus),'aemom_inp_clus','read_moments_clus')
      aemom_inp_clus=0.0_dblprec

      allocate(Landeg_ch_clus(NA_clus,nchmax_clus,conf_num),stat=i_stat)
      call memocc(i_stat,product(shape(Landeg_ch_clus))*kind(Landeg_ch_clus),'Landeg_ch_clus','read_moments_clus')
      Landeg_ch_clus=0.0_dblprec

      allocate(ham_clus%ind_mom_clus(NA_clus,nchmax_clus),stat=i_stat)
      call memocc(i_stat,product(shape(ham_clus%ind_mom_clus))*kind(ham_clus%ind_mom_clus),'ham_clus%ind_mom_clus','read_moments_clus')
      ham_clus%ind_mom_clus=0

      i_err=0

      if(set_landeg==1) then
         ! If the induced magnetic moments flag is on one must read whether a certain moment is induced or not
         if (do_lsf=='N') then
            if (ind_mom_flag=='Y') then
               do while(i_err==0)
                  read(ifileno,*,iostat=i_err) isite,ichem,ammom_inp_clus(isite,ichem,1), &
                     aemom_inp_clus(1:3,isite,ichem,1),Landeg_ch_clus(isite,ichem,1),     &
                     ham_clus%ind_mom_clus(isite,ichem)
                  aemom_tmp=norm2(aemom_inp_clus(:,isite,ichem,1))
                  aemom_inp_clus(:,isite,ichem,1)=aemom_inp_clus(:,isite,ichem,1)/aemom_tmp
               end do
            else
               do while(i_err==0)
                  read(ifileno,*,iostat=i_err) isite,ichem,ammom_inp_clus(isite,ichem,1), &
                     aemom_inp_clus(1:3,isite,ichem,1),Landeg_ch_clus(isite,ichem,1)
                  aemom_tmp=norm2(aemom_inp_clus(:,isite,ichem,1))
                  aemom_inp_clus(:,isite,ichem,1)=aemom_inp_clus(:,isite,ichem,1)/aemom_tmp
               end do
            endif
         else ! LSF, not yet with induced moments
            ! For LSF modified momfile requires configuration number as first column
            do while(i_err==0)
               read(ifileno,*,iostat=i_err) iconf,isite,ichem,ammom_inp_clus(isite,ichem,iconf), &
                  aemom_inp_clus(1:3,isite,ichem,iconf),Landeg_ch_clus(isite,ichem,iconf)
               aemom_tmp=norm2(aemom_inp_clus(:,isite,ichem,iconf))
               aemom_inp_clus(:,isite,ichem,iconf)=aemom_inp_clus(:,isite,ichem,iconf)/aemom_tmp
            end do
         endif
      else
         Landeg_ch_clus=Landeg_global
         if (do_lsf=='N') then
            if (ind_mom_flag=='Y') then
               ! If the induced magnetic moments flag is on one must read whether a certain moment is induced or not
               do while(i_err==0)
                  read(ifileno,*,iostat=i_err) isite,ichem,ammom_inp_clus(isite,ichem,1), &
                     aemom_inp_clus(1:3,isite,ichem,1),ham_clus%ind_mom_clus(isite,ichem)
                  aemom_tmp=norm2(aemom_inp_clus(:,isite,ichem,1))
                  aemom_inp_clus(:,isite,ichem,1)=aemom_inp_clus(:,isite,ichem,1)/aemom_tmp
               end do
            else
               do while(i_err==0)
                  read(ifileno,*,iostat=i_err) isite,ichem,ammom_inp_clus(isite,ichem,1), &
                     aemom_inp_clus(1:3,isite,ichem,1)
                  aemom_tmp=norm2(aemom_inp_clus(:,isite,ichem,1))
                  aemom_inp_clus(1:3,isite,ichem,1)=aemom_inp_clus(1:3,isite,ichem,1)/aemom_tmp
               end do
            endif
         else   ! LSF
            do while(i_err==0)
               read(ifileno,*,iostat=i_err) iconf,isite,ichem,ammom_inp_clus(isite,ichem,iconf), &
                  aemom_inp_clus(1:3,isite,ichem,iconf)
               aemom_tmp=norm2(aemom_inp_clus(:,isite,ichem,iconf))
               aemom_inp_clus(:,isite,ichem,iconf)=aemom_inp_clus(:,isite,ichem,iconf)/aemom_tmp
            end do
         endif
      end if
      close(ifileno)
      !
   end subroutine read_moments_clus

   !---------------------------------------------------------------------------
   !  SUBROUTINE: read_exchange_getMaxNoShells_clus
   !> @brief
   !> Get the max no of exchange shells and lines for the cluster
   !> Helper for read exchange
   !>
   !> @author
   !> Jonathan Chico based on the routines by Anders Bergman
   !---------------------------------------------------------------------------
   subroutine read_exchange_getMaxNoShells_clus(no_shells_clus,flines_clus,ifileno2,&
      do_ralloy)

      implicit none

      integer, intent(in)                    :: ifileno2
      integer, intent(in)                    :: do_ralloy   !< Random alloy simulation (0/1)
      integer, intent(out)                   :: no_shells_clus,flines_clus
      integer                                :: mtype
      integer                                :: itype,jtype,isite,jsite,ichem,jchem
      integer                                :: i_stat,i_all
      integer, dimension(:,:,:), allocatable :: nn_tmp_clus

      ! Open input file
      ! Check if input file is for random alloy
      flines_clus=0
      mtype=0

      allocate(nn_tmp_clus(NT_clus,max(NT_clus,nchmax_clus),nchmax_clus),stat=i_stat)
      call memocc(i_stat,product(shape(nn_tmp_clus))*kind(nn_tmp_clus),'nn_tmp_clus','read_exchange')
      nn_tmp_clus=0

      ! Pre-read file to get max no. exchange shells and no. lines
      do
         if(do_ralloy==0) then
            read(ifileno2,*,end=200)  isite,jsite
            ichem=1
            jchem=1
         else
            read(ifileno2,*,end=200)  isite, jsite, ichem, jchem
         end if
         flines_clus=flines_clus+1
         itype=atype_inp_clus(isite)
         jtype=atype_inp_clus(jsite)
         if(do_ralloy==0) then
            nn_tmp_clus(itype,jtype,jchem)=nn_tmp_clus(itype,jtype,jchem)+1
         else
            nn_tmp_clus(itype,ichem,jchem)=nn_tmp_clus(itype,ichem,jchem)+1
         end if
      end do
      200 continue

      rewind(ifileno2)
      no_shells_clus=0
      if (do_ralloy==0) then
         do itype=1,NT_clus
            no_shells_clus=max(sum(nn_tmp_clus(itype,:,:)),no_shells_clus)
         end do
      else
         do itype=1,NT_clus
            do ichem=1,Nchmax_clus
               do jchem=1,Nchmax_clus
                  no_shells_clus=max(sum(nn_tmp_clus(itype,ichem,:)),no_shells_clus)
               end do
            end do
         end do
      endif
      i_all=-product(shape(nn_tmp_clus))*kind(nn_tmp_clus)
      deallocate(nn_tmp_clus,stat=i_stat)
      call memocc(i_stat,i_all,'nn_tmp_clus','read_exchange')

   end subroutine read_exchange_getMaxNoShells_clus

   !---------------------------------------------------------------------------
   !  SUBROUTINE: read_exchange_getNeighVec_clus
   !> @brief
   !> Obtaining the neighbour vector for the custer
   !>
   !> @author
   !> Jonathan Chico based of Anders Bergman routines
   !---------------------------------------------------------------------------
   subroutine read_exchange_getNeighVec_clus(r_red_clus,r_tmp_clus,isite_c,jsite_c, &
      maptype,posfiletype)

      implicit none

      integer, intent(in)                       :: isite_c,jsite_c
      integer, intent(in)                       :: maptype        !< Format for input data (1=direct format,2=bgfm style)
      character(len=1), intent(in)              :: posfiletype    !< posfile type (D)irect or (C)arteisian coordinates in posfile
      real(dblprec), dimension(3), intent(in)   :: r_tmp_clus
      real(dblprec), dimension(3), intent(out)  :: r_red_clus

      if(maptype==2) then
         ! Calculate proper neighbour vector (from "bgfm")
         r_red_clus(1)=Bas_clus(1,jsite_c)-Bas_clus(1,isite_c)+C1_clus(1)*r_tmp_clus(1)+C2_clus(1)*r_tmp_clus(2)+C3_clus(1)*r_tmp_clus(3)
         r_red_clus(2)=Bas_clus(2,jsite_c)-Bas_clus(2,isite_c)+C1_clus(2)*r_tmp_clus(1)+C2_clus(2)*r_tmp_clus(2)+C3_clus(2)*r_tmp_clus(3)
         r_red_clus(3)=Bas_clus(3,jsite_c)-Bas_clus(3,isite_c)+C1_clus(3)*r_tmp_clus(1)+C2_clus(3)*r_tmp_clus(2)+C3_clus(3)*r_tmp_clus(3)
      else
         ! Calculates neighbour vectors from direct coordinates or Cartesian
         ! coordinates, corresponding to how the atomic positions are entered
         if (posfiletype=='C') then
            r_red_clus=r_tmp_clus
         elseif (posfiletype=='D') then
            r_red_clus(1)=r_tmp_clus(1)*C1_clus(1)+r_tmp_clus(2)*C2_clus(1)+r_tmp_clus(3)*C3_clus(1)
            r_red_clus(2)=r_tmp_clus(1)*C1_clus(2)+r_tmp_clus(2)*C2_clus(2)+r_tmp_clus(3)*C3_clus(2)
            r_red_clus(3)=r_tmp_clus(1)*C1_clus(3)+r_tmp_clus(2)*C2_clus(3)+r_tmp_clus(3)*C3_clus(3)
         else
            stop 'Only posfiletype = C or D is currently supported'
         endif
      end if
   end subroutine read_exchange_getNeighVec_clus

   !---------------------------------------------------------------------------
   !  SUBROUTINE: read_anisotropy_clus
   !> @brief
   !> Reading the anisotopies for the cluster atoms
   !>
   !> @author
   !> Jonathan Chico based of Anders Bergman routines
   !---------------------------------------------------------------------------
   subroutine read_anisotropy_clus()
      !
      implicit none
      !
      integer :: m,iat

      open(ifileno, file=adjustl(kfile_clus))

      do m=1, NA_clus
         read (ifileno,*) iat,anisotropytype_clus(iat,1),anisotropy_clus(iat,1,1),  &
            anisotropy_clus(iat,2,1),anisotropy_clus(iat,3,1),                      &
            anisotropy_clus(iat,4,1),anisotropy_clus(iat,5,1),                      &
            anisotropy_clus(iat,6,1)
      enddo

      if (mult_axis_clus=='Y') then
         do m=1,NA_clus
            read(ifileno,*) iat,anisotropytype_diff_clus(iat,1),                    &
            anisotropy_diff_clus(iat,1,1),anisotropy_diff_clus(iat,2,1),            &
            anisotropy_diff_clus(iat,3,1),anisotropy_diff_clus(iat,4,1),            &
            anisotropy_diff_clus(iat,5,1),anisotropy_diff_clus(iat,6,1)
         enddo

      endif

      close (ifileno)

   end subroutine read_anisotropy_clus

   !---------------------------------------------------------------------------
   !  SUBROUTINE: allocate_hamiltonianinput_clus
   !> @brief
   !> Routine to allocate the needed data for the input of the cluster Hamiltonian
   !>
   !> @author
   !> Jonathan Chico based of Anders Bergman routines
   !---------------------------------------------------------------------------
   subroutine allocate_hamiltonianinput_clus(conf_num,no_shells_clus,flag)

      implicit none

      integer, intent(in) :: flag      !< Allocate or deallocate (1/-1)
      integer, intent(in) :: conf_num  !< Number of configurations for LSF
      integer, intent(in), optional :: no_shells_clus !< Parameter limiting number of exchange coupling shells for the cluster

      integer :: i_all, i_stat

      if (flag>0) then
         allocate(jc_clus(NT_clus,no_shells_clus,Nchmax_clus,Nchmax_clus,conf_num),stat=i_stat)
         call memocc(i_stat,product(shape(jc_clus))*kind(jc_clus),'jc_clus','allocate_hamiltonianinput_clus')
         jc_clus=0.0_dblprec
         allocate(redcoord_clus(NT_clus,no_shells_clus,3),stat=i_stat)
         call memocc(i_stat,product(shape(redcoord_clus))*kind(redcoord_clus),'redcoord_clus','allocate_hamiltonianinput_clus')
         redcoord_clus=0.0_dblprec
         allocate(NN_clus(NT_clus),stat=i_stat)
         call memocc(i_stat,product(shape(NN_clus))*kind(NN_clus),'NN_clus','allocate_hamiltonianinput_clus')
         NN_clus=0
         allocate(dm_nn_clus(NT_clus),stat=i_stat)
         call memocc(i_stat,product(shape(dm_nn_clus))*kind(dm_nn_clus),'dm_nn_clus','allocate_hamiltonianinput_clus')
         dm_nn_clus=0
      else
         i_all=-product(shape(NN_clus))*kind(NN_clus)
         deallocate(NN_clus,stat=i_stat)
         call memocc(i_stat,i_all,'NN_clus','allocate_hamiltonianinput_clus')
         i_all=-product(shape(jfile_clus))*kind(jfile_clus)
         deallocate(jfile_clus,stat=i_stat)
         call memocc(i_stat,i_all,'jfile_clus','allocate_hamiltonianinput_clus')
         i_all=-product(shape(dm_nn_clus))*kind(dm_nn_clus)
         deallocate(dm_nn_clus,stat=i_stat)
         call memocc(i_stat,i_all,'dm_nn_clus','allocate_hamiltonianinput_clus')
         i_all=-product(shape(jc_clus))*kind(jc_clus)
         deallocate(jc_clus,stat=i_stat)
         call memocc(i_stat,i_all,'jc_clus','allocate_hamiltonianinput_clus')
         i_all=-product(shape(atype_inp_clus))*kind(atype_inp_clus)
         deallocate(atype_inp_clus,stat=i_stat)
         call memocc(i_stat,i_all,'atype_inp_clus','allocate_hamiltonianinput_clust')
         i_all=-product(shape(anumb_inp_clus))*kind(anumb_inp_clus)
         deallocate(anumb_inp_clus,stat=i_stat)
         call memocc(i_stat,i_all,'anumb_inp_clus','allocate_hamiltonianinput_clus')
         i_all=-product(shape(redcoord_clus))*kind(redcoord_clus)
         deallocate(redcoord_clus,stat=i_stat)
         call memocc(i_stat,i_all,'redcoord_clus','allocate_hamiltonianinput_clus')
         if (allocated(dm_redcoord_clus)) then
            i_all=-product(shape(dm_redcoord_clus))*kind(dm_redcoord_clus)
            deallocate(dm_redcoord_clus,stat=i_stat)
            call memocc(i_stat,i_all,'dm_redcoord_clus','allocate_hamiltonianinput_clus')
         endif
         if (allocated(dm_inpvect_clus)) then
            i_all=-product(shape(dm_inpvect_clus))*kind(dm_inpvect_clus)
            deallocate(dm_inpvect_clus,stat=i_stat)
            call memocc(i_stat,i_all,'dm_inpvect_clus','allocate_hamiltonianinput_clus')
         endif

      endif

   end subroutine  allocate_hamiltonianinput_clus

   !---------------------------------------------------------------------------
   !  SUBROUTINE: allocate_anisotropy_input_clus
   !> @brief
   !> Routine to allocate the needed data for the anisotropy in the impurity cluster
   !>
   !> @author
   !> Jonathan Chico based of Anders Bergman routines
   !---------------------------------------------------------------------------
   subroutine allocate_anisotropy_input_clus(flag)

      implicit none

      integer, intent(in) :: flag  !< Allocate or deallocate (1/-1)

      integer :: i_all, i_stat

      if (flag>0) then
         allocate(anisotropytype_clus(NA_clus,Nchmax_clus),stat=i_stat)
         call memocc(i_stat,product(shape(anisotropytype_clus))*kind(anisotropytype_clus),'anisotropytype_clus','allocate_anisotropy_input_clus')
         anisotropytype_clus=0
         allocate(anisotropy_clus(NA_clus,6,Nchmax_clus),stat=i_stat)
         call memocc(i_stat,product(shape(anisotropy_clus))*kind(anisotropy_clus),'anisotropy_clus','allocate_anisotropy_input_clus')
         anisotropy_clus=0.0_dblprec

         if (mult_axis_clus=='Y') then
            allocate(anisotropytype_diff_clus(NA_clus,Nchmax_clus),stat=i_stat)
            call memocc(i_stat,product(shape(anisotropytype_diff_clus))*kind(anisotropytype_diff_clus),'anisotropytype_diff_clus','allocate_anisotropy_input_clus')
            allocate(anisotropy_diff_clus(NA_clus,6,Nchmax_clus),stat=i_stat)
            call memocc(i_stat,product(shape(anisotropy_diff_clus))*kind(anisotropy_diff_clus),'anisotropy_diff_clus','allocate_anisotropy_input_clus')
            anisotropy_diff_clus=0.0_dblprec
         endif
      else
         i_all=-product(shape(anisotropy_clus))*kind(anisotropy_clus)
         deallocate(anisotropy_clus,stat=i_stat)
         call memocc(i_stat,i_all,'anisotropy_clus','allocate_anisotropy_input_clus')
         i_all=-product(shape(anisotropytype_clus))*kind(anisotropytype_clus)
         deallocate(anisotropytype_clus,stat=i_stat)
         call memocc(i_stat,i_all,'anisotropytype_clus','allocate_anisotropy_input_clus')

         if (mult_axis_clus=='Y') then
            i_all=-product(shape(anisotropy_diff_clus))*kind(anisotropy_diff_clus)
            deallocate(anisotropy_diff_clus,stat=i_stat)
            call memocc(i_stat,i_all,'anisotropy_diff_clus','allocate_anisotropy_input_clus')
            i_all=-product(shape(anisotropytype_diff_clus))*kind(anisotropytype_diff_clus)
            deallocate(anisotropytype_diff_clus,stat=i_stat)
            call memocc(i_stat,i_all,'anisotropytype_diff_clus','allocate_anisotropy_input_clus')
         endif

      endif

   end subroutine allocate_anisotropy_input_clus

   !---------------------------------------------------------------------------
   !  SUBROUTINE: allocate_cluster_hamiltoniandata
   !> @brief
   !> Routine to allocate the needed data for the actual cluster Hamiltonian
   !>
   !> @author
   !> Jonathan Chico based of Anders Bergman routines
   !---------------------------------------------------------------------------
   subroutine allocate_cluster_hamiltoniandata(do_jtensor,NA_clus,conf_num,max_no_neigh_clus,flag)

      implicit none

      integer, intent(in) :: flag  !< Allocate or deallocate (1/-1)
      integer, intent(in) ::  do_jtensor        !< Use SKKR style exchange tensor (0=off, 1=on)
      integer, intent(in), optional :: NA_clus   !< Number of atoms in the cluster unit cell
      integer, intent(in), optional :: conf_num  !< Number of configurations for LSF
      integer, intent(in), optional :: max_no_neigh_clus !< Calculated maximum of neighbours for exchange

      integer :: i_all, i_stat
      ! Exchange
      if(flag>0) then
         allocate(ham_clus%nlistsize_clus(NA_clus),stat=i_stat)
         call memocc(i_stat,product(shape(ham_clus%nlistsize_clus))*kind(ham_clus%nlistsize_clus),'ham_clus%nlistsize_clus','allocate_cluster_hamiltoniandata')
         ham_clus%nlistsize_clus=0
         allocate(ham_clus%nlist_clus(max_no_neigh_clus,NA_clus),stat=i_stat)
         call memocc(i_stat,product(shape(ham_clus%nlist_clus))*kind(ham_clus%nlist_clus),'ham_clus%nlist_clus','allocate_cluster_hamiltoniandata')
         ham_clus%nlist_clus=0
         allocate(ham_clus%ncoup_clus(max_no_neigh_clus,NA_clus,conf_num),stat=i_stat)
         call memocc(i_stat,product(shape(ham_clus%ncoup_clus))*kind(ham_clus%ncoup_clus),'ham_clus%ncoup_clus','allocate_cluster_hamiltoniandata')
         ham_clus%ncoup_clus=0.0_dblprec
      else
         i_all=-product(shape(ham_clus%nlistsize_clus))*kind(ham_clus%nlistsize_clus)
         deallocate(ham_clus%nlistsize_clus,stat=i_stat)
         call memocc(i_stat,i_all,'ham_clus%nlistsize_clus','allocate_cluster_hamiltoniandata')
         i_all=-product(shape(ham_clus%nlist_clus))*kind(ham_clus%nlist_clus)
         deallocate(ham_clus%nlist_clus,stat=i_stat)
         call memocc(i_stat,i_all,'ham_clus%nlist_clus','allocate_cluster_hamiltoniandata')

         if (do_jtensor/=1) then
            i_all=-product(shape(ham_clus%ncoup_clus))*kind(ham_clus%ncoup_clus)
            deallocate(ham_clus%ncoup_clus,stat=i_stat)
            call memocc(i_stat,i_all,'ham_clus%ncoup_clus','allocate_cluster_hamiltoniandata')
         end if
      end if

   end subroutine allocate_cluster_hamiltoniandata

   !-----------------------------------------------------------------------------
   ! SUBROUTINE: allocate_cluster_dmhamiltoniandata
   !> @brief Allocate arrays for Dzyaloshinskii-Moriya Hamiltonian
   !> @author
   !> Jonathan Chico based of Anders Bergman routines
   !-----------------------------------------------------------------------------
   subroutine allocate_cluster_dmhamiltoniandata(NA_clus,max_no_dmneigh_clus,flag)

      implicit none

      integer, intent(in) :: flag !< Allocate or deallocate (1/-1)
      integer, optional, intent(in) :: NA_clus  !< Number of atoms in the cluster unit cell
      integer, optional, intent(in) :: max_no_dmneigh_clus !< Calculated number of neighbours with DM interactions in the cluster

      integer :: i_all, i_stat

      if(flag>0) then
         allocate(ham_clus%dmlistsize_clus(NA_clus),stat=i_stat)
         call memocc(i_stat,product(shape(ham_clus%dmlistsize_clus))*kind(ham_clus%dmlistsize_clus),'ham_clus%dmlistsize_clus','allocate_cluster_dmhamiltoniandata')
         ham_clus%dmlistsize_clus=0
         allocate(ham_clus%dmlist_clus(max_no_dmneigh_clus,NA_clus),stat=i_stat)
         call memocc(i_stat,product(shape(ham_clus%dmlist_clus))*kind(ham_clus%dmlist_clus),'ham_clus%dmlist_clus','allocate_cluster_dmhamiltoniandata')
         ham_clus%dmlist_clus=0
         allocate(ham_clus%dm_vect_clus(3,max_no_dmneigh_clus,NA_clus),stat=i_stat)
         call memocc(i_stat,product(shape(ham_clus%dm_vect_clus))*kind(ham_clus%dm_vect_clus),'ham_clus%dm_vect_clus','allocate_cluster_dmhamiltoniandata')
         ham_clus%dm_vect_clus=0.0_dblprec
      else
         i_all=-product(shape(ham_clus%dmlistsize_clus))*kind(ham_clus%dmlistsize_clus)
         deallocate(ham_clus%dmlistsize_clus,stat=i_stat)
         call memocc(i_stat,i_all,'ham_clus%dmlistsize_clus','allocate_cluster_dmhamiltoniandata')
         i_all=-product(shape(ham_clus%dmlist_clus))*kind(ham_clus%dmlist_clus)
         deallocate(ham_clus%dmlist_clus,stat=i_stat)
         call memocc(i_stat,i_all,'ham_clus%dmlist_clus','allocate_cluster_dmhamiltoniandata')
         i_all=-product(shape(ham_clus%dm_vect_clus))*kind(ham_clus%dm_vect_clus)
         deallocate(ham_clus%dm_vect_clus,stat=i_stat)
         call memocc(i_stat,i_all,'ham_clus%dm_vect_clus','allocate_cluster_dmhamiltoniandata')
      end if

   end subroutine allocate_cluster_dmhamiltoniandata

   !-----------------------------------------------------------------------------
   ! SUBROUTINE: allocate_cluster_anisotropies
   !> @brief Allocate arrays for anisotropy for the cluster
   !> @author
   !> Jonathan Chico based of Anders Bergman routines
   !-----------------------------------------------------------------------------
   subroutine allocate_cluster_anisotropies(NA_clus,mult_axis_clus,flag)

      implicit none

      integer, intent(in) :: flag  !< Allocate or deallocate (1/-1)
      character(len=1), intent(in) :: mult_axis_clus
      integer, intent(in), optional :: NA_clus  !< Number of atoms in the cluster unit cell
      ! .. Local variables
      integer :: i_all, i_stat

      ! Allocate arrays for anisotropy
      if(flag>0) then
         allocate(ham_clus%taniso_clus(NA_clus),stat=i_stat)
         call memocc(i_stat,product(shape(ham_clus%taniso_clus))*kind(ham_clus%taniso_clus),'ham_clus%taniso_clus','setup_cluster_anisotropies')
         allocate(ham_clus%eaniso_clus(3,NA_clus),stat=i_stat)
         call memocc(i_stat,product(shape(ham_clus%eaniso_clus))*kind(ham_clus%eaniso_clus),'ham_clus%eaniso_clus','setup_cluster_anisotropies')
         allocate(ham_clus%kaniso_clus(2,NA_clus),stat=i_stat)
         call memocc(i_stat,product(shape(ham_clus%kaniso_clus))*kind(ham_clus%kaniso_clus),'ham_clus%kaniso_clus','setup_cluster_anisotropies')
         allocate(ham_clus%sb_clus(NA_clus),stat=i_stat)
         call memocc(i_stat,product(shape(ham_clus%sb_clus))*kind(ham_clus%sb_clus),'ham_clus%sb_clus','setup_cluster_anisotropies')

         if (mult_axis_clus=='Y') then
            allocate(ham_clus%taniso_diff_clus(NA_clus),stat=i_stat)
            call memocc(i_stat,product(shape(ham_clus%taniso_diff_clus))*kind(ham_clus%taniso_diff_clus),'ham_clus%taniso_diff_clus','setup_cluster_anisotropies')
            allocate(ham_clus%eaniso_diff_clus(3,NA_clus),stat=i_stat)
            call memocc(i_stat,product(shape(ham_clus%eaniso_diff_clus))*kind(ham_clus%eaniso_diff_clus),'ham_clus%eaniso_diff_clus','setup_cluster_anisotropies')
            allocate(ham_clus%kaniso_diff_clus(2,NA_clus),stat=i_stat)
            call memocc(i_stat,product(shape(ham_clus%kaniso_diff_clus))*kind(ham_clus%kaniso_diff_clus),'ham_clus%kaniso_diff_clus','setup_cluster_anisotropies')
            allocate(ham_clus%sb_diff_clus(NA_clus),stat=i_stat)
            call memocc(i_stat,product(shape(ham_clus%sb_diff_clus))*kind(ham_clus%sb_diff_clus),'ham_clus%sb_diff_clus','setup_cluster_anisotropies')
         endif

      else
         i_all=-product(shape(ham_clus%taniso_clus))*kind(ham_clus%taniso_clus)
         deallocate(ham_clus%taniso_clus,stat=i_stat)
         call memocc(i_stat,i_all,'ham_clus%taniso_clus','setup_cluster_anisotropies')
         i_all=-product(shape(ham_clus%eaniso_clus))*kind(ham_clus%eaniso_clus)
         deallocate(ham_clus%eaniso_clus,stat=i_stat)
         call memocc(i_stat,i_all,'ham_clus%eaniso_clus','setup_cluster_anisotropies')
         i_all=-product(shape(ham_clus%kaniso_clus))*kind(ham_clus%kaniso_clus)
         deallocate(ham_clus%kaniso_clus,stat=i_stat)
         call memocc(i_stat,i_all,'ham_clus%kaniso_clus','setup_cluster_anisotropies')
         i_all=-product(shape(ham_clus%sb_clus))*kind(ham_clus%sb_clus)
         deallocate(ham_clus%sb_clus,stat=i_stat)
         call memocc(i_stat,i_all,'ham_clus%sb_clus','setup_cluster_anisotropies')

         if (mult_axis_clus=='Y') then
            i_all=-product(shape(ham_clus%taniso_diff_clus))*kind(ham_clus%taniso_diff_clus)
            deallocate(ham_clus%taniso_diff_clus,stat=i_stat)
            call memocc(i_stat,i_all,'ham_clus%taniso_diff_clus','setup_cluster_anisotropies')
            i_all=-product(shape(ham_clus%eaniso_diff_clus))*kind(ham_clus%eaniso_diff_clus)
            deallocate(ham_clus%eaniso_diff_clus,stat=i_stat)
            call memocc(i_stat,i_all,'ham_clus%eaniso_diff_clus','setup_cluster_anisotropies')
            i_all=-product(shape(ham_clus%kaniso_diff_clus))*kind(ham_clus%kaniso_diff_clus)
            deallocate(ham_clus%kaniso_diff_clus,stat=i_stat)
            call memocc(i_stat,i_all,'ham_clus%kaniso_diff_clus','setup_cluster_anisotropies')
            i_all=-product(shape(ham_clus%sb_diff_clus))*kind(ham_clus%sb_diff_clus)
            deallocate(ham_clus%sb_diff_clus,stat=i_stat)
            call memocc(i_stat,i_all,'ham_clus%sb_diff_clus','setup_cluster_anisotropies')
         endif
      end if

   end subroutine allocate_cluster_anisotropies

   !---------------------------------------------------------------------------
   !  SUBROUTINE: read_parameters_cluster
   !> @brief
   !> Read input parameters for the cluster.
   !
   !> @author
   !> Jonathan Chico based on the routines by Anders Bergman
   !---------------------------------------------------------------------------
   subroutine read_parameters_cluster(ifile,conf_num)
      use FileParser
      use ErrorHandling
      use InputData, only: do_anisotropy

      implicit none

      ! ... Formal Arguments ...
      integer, intent(in) :: ifile     !< File to read from
      integer, intent(in) :: conf_num  !< Number of configurations for LSF
      !
      ! ... Local Variables ...
      character(len=50) :: keyword, cache, string
      character(len=30) ::pre_jfile                           !< File name for exchange couplings
      integer :: rd_len, i_err, i_errb,i_stat,ii,i
      logical :: comment

      do
         10     continue
         ! Read file character for character until first whitespace
         keyword=""
         call bytereader(keyword,rd_len,ifile,i_errb)

         ! converting Capital letters
         call caps2small(keyword)

         ! check for comment markers (currently % and #)
         comment=(scan(trim(keyword),'%')==1).or.(scan(trim(keyword),'#')==1).or.&
            (scan(trim(keyword),'*')==1).or.(scan(trim(keyword),'=')==1.or.&
            (scan(trim(keyword),'!')==1))

         if (comment) then
            read(ifile,*)
         else
            ! Parse keyword
            keyword=trim(keyword)
            select case(keyword)

            !!!!!!!!!!!!!!!!!!!!! START OF VARIABLES FOR THE CLUSTER !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

         case('cell_clus')
            read(ifile,*,iostat=i_err) C1_clus, C2_clus, C3_clus
            if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err

         case('ncell_clus')
            read(ifile,*,iostat=i_err) N1_clus, N2_clus, N3_clus
            if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err

         case('do_cluster')
            read(ifile,*,iostat=i_err) do_cluster
            if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err

         case('momfile_clus')
            read(ifile,'(a)',iostat=i_err) cache
            if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err
            momfile_clus=adjustl(trim(cache))

         case('posfile_clus')
            read(ifile,'(a)',iostat=i_err) cache
            if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err
            posfile_clus=adjustl(trim(cache))

         case('anisotropy_clus')
            do_anisotropy_clus=1
            do_anisotropy=1
            read(ifile,'(a)',iostat=i_err) cache
            if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err
            kfile_clus=adjustl(trim(cache))

            !> - exchange
            !! Name of exchange file
         case('exchange_clus')
            read(ifile,'(a)',iostat=i_err) cache
            if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err
            allocate(jfile_clus(conf_num),stat=i_stat)
            call memocc(i_stat,product(shape(jfile_clus))*kind(jfile_clus),'jfile_clus','inputhandler')
            if (conf_num==1) then
               jfile_clus=adjustl(trim(cache))
               ! For the LSF method there are several jfiles that must be read, depending on the number of configurations
            else
               pre_jfile=adjustl(trim(cache))
               i=len_trim(pre_jfile)-1
               do ii=1,conf_num
                  if (ii<10) then
                     write(string,'("(a",i0,",i1)")') i
                     write(jfile_clus(ii),string) adjustl(trim(pre_jfile)),ii
                  else if (ii>9.and.ii<100) then
                     write(string,'("(a",i0,",i2)")') i
                     write(jfile_clus(ii),string) adjustl(trim(pre_jfile)),ii
                  else if (ii>99.and.ii<1000) then
                     write(string,'("(a",i0,",i3)")') i
                     write(jfile_clus(ii),string) adjustl(trim(pre_jfile)),ii
                  else if (ii>999.and.ii<10000) then
                     write(string,'("(a",i0,",i4)")') i
                     write(jfile_clus(ii),string) adjustl(trim(pre_jfile)),ii
                  endif
               enddo
            endif

         case('dm_clus')
            read(ifile,'(a)',iostat=i_err) cache
            if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err
            dmfile_clus=trim(adjustl(cache))
            call ErrorHandling_check_file_exists(dmfile_clus, &
               'Please specify dm <dmfile> where <dmfile> is a valid dm interaction file')
            !!!!!!!!!!!!!!!!!!!!! END OF VARIABLES FOR THE CLUSTER !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         end select
      endif
      ! End of file
      if (i_errb==20) goto 20
      ! End of row
      if (i_errb==10) goto 10
   end do

   20  continue

   rewind(ifile)
   return

end subroutine read_parameters_cluster

end module clusters
