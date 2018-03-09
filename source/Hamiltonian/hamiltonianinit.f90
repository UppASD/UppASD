!-------------------------------------------------------------------------------
!> MODULE: HamiltonianInit
!> @brief
!> Routines for mounting the Hamiltonian, including exchange, DM, anisotropies
!> @copyright
!! Copyright (C) 2008-2018 UppASD group
!! This file is distributed under the terms of the
!! GNU General Public License. 
!! See http://www.gnu.org/copyleft/gpl.txt
!-------------------------------------------------------------------------------
module HamiltonianInit
  use Parameters
  use Profiling
  implicit none

  integer, dimension(:,:,:), allocatable :: nm !< Neighbour map
  integer, dimension(:,:), allocatable :: nmdim !< Dimension of neighbour map

  private
  public :: setup_hamiltonian


contains


  !-------------------------------------------------------------------------------
  !> SUBROUTINE: setup_hamiltonian
  !> @brief
  !>  Main routine for mounting the Hamiltonian
  !! @todo Change mounting of DM Hamiltonian to follow the Heisenberg approach
  !! @todo Introduce additional flags for symmetry in order to enable use of symops for other couplings than Heisenberg
  !-------------------------------------------------------------------------------
  subroutine setup_hamiltonian(NT,NA,N1,N2,N3,Nchmax,do_ralloy,Natom_full,Natom,acellnumb,acellnumbrev,&
             achtype,atype_ch,asite_ch,achem_ch,atype,anumb,alat,C1,C2,C3,Bas,ammom_inp,coord,BC1,BC2,BC3,&
             sym,do_jtensor,max_no_neigh,max_no_shells,nn,nlistsize,nlist,redcoord,jc,jcD,jc_tens,ncoup,&
             ncoupD,j_tens,do_dm,max_no_dmshells,max_no_dmneigh,dm_nn,dmlistsize,dmlist,dm_vect,dm_redcoord,&
             dm_inpvect,do_anisotropy,taniso,taniso_diff,random_anisotropy_density,anisotropytype,&
             anisotropytype_diff,anisotropy,anisotropy_diff,sb,sb_diff,eaniso,eaniso_diff,kaniso,&
             kaniso_diff,random_anisotropy,mult_axis,mconf,conf_num,fs_nlistsize,fs_nlist,nind,&
             map_multiple,do_lsf,lsf_field,exc_inter,do_bq,nn_bq_tot,max_no_bqshells,bq_nn,bqlistsize,&
             bqlist,bq_redcoord,jc_bq,j_bq,do_biqdm,nn_biqdm_tot,max_no_biqdmshells,biqdm_nn,&
             biqdmlistsize,biqdmlist,biqdm_redcoord,biqdm_inpvect,biqdm_vect,do_pd,nn_pd_tot,&
             max_no_pdshells,pd_nn,pdlistsize,pdlist,pd_redcoord,pd_inpvect,pd_vect,do_dip,qdip,&
             ind_mom,ind_nlistsize,ind_nlist,ind_tol,ind_mom_flag,&
             do_prnstruct,do_sortcoup,simid)

    use NeighbourMap, only : setup_nm
    use InputHandler , only : allocate_hamiltonianinput
    use HamiltonianData, only : allocate_hamiltoniandata, allocate_anisotropies, allocate_dipole, &
         allocate_dmhamiltoniandata, allocate_pdhamiltoniandata, allocate_biqdmhamiltoniandata, &
         allocate_bqhamiltoniandata
    use PrintHamiltonian
    use LSF, only: LSF_datareshape
    use InducedMoments

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
    integer, intent(inout) :: Natom   !< Number of atoms in system
    integer, dimension(:), intent(inout) :: acellnumb    !< List for translating atom no. in full cell to actual cell
    integer, dimension(:), intent(inout) :: acellnumbrev !< List for translating atom no. in actual cell to full cell
    integer, dimension(:), intent(inout) :: achtype   !< Chemical type of atoms (full list)
    integer, dimension(:), intent(inout) :: atype_ch  !< Actual type of atom for dilute system
    integer, dimension(:), intent(inout) :: asite_ch  !< Actual site of atom for dilute system
    integer, dimension(:), intent(inout) :: achem_ch  !< Chemical type of atoms (reduced list)
    integer, dimension(Natom), intent(inout) :: atype !< Type of atom
    integer, dimension(Natom), intent(inout) :: anumb !< Atom number in cell
    real(dblprec) :: alat                         !< Lattice parameter
    real(dblprec), dimension(3), intent(in) :: C1 !< First lattice vector
    real(dblprec), dimension(3), intent(in) :: C2 !< Second lattice vector
    real(dblprec), dimension(3), intent(in) :: C3 !< Third lattice vector
    real(dblprec), dimension(3,NA), intent(inout) :: Bas !< Coordinates for basis atoms
    real(dblprec), dimension(:,:,:), intent(inout) :: ammom_inp !< Magnetic moment directions from input (for alloys)
    real(dblprec), dimension(:,:), allocatable, intent(inout) :: coord !< Coordinates of atoms
    character(len=1), intent(in) :: BC1  !< Boundary conditions in x-direction
    character(len=1), intent(in) :: BC2  !< Boundary conditions in y-direction
    character(len=1), intent(in) :: BC3  !< Boundary conditions in z-direction
    ! Heisenberg exchange variables
    integer, intent(in) :: sym              !< Symmetry of system (0-3)
    integer, intent(in) :: do_jtensor       !< Use SKKR style exchange tensor (0=off, 1=on)
    integer, intent(inout) :: max_no_neigh  !< Calculated maximum of neighbours for exchange
    integer, intent(inout) :: max_no_shells !< Calculated maximum of shells for exchange
    integer, dimension(:), intent(in) :: nn !< Number of neighbour shells for exchange
    integer, dimension(:), allocatable, intent(inout) :: nlistsize !< Size of neighbour list for Heisenberg exchange couplings
    integer, dimension(:,:), allocatable, intent(inout) :: nlist   !< Neighbour list for Heisenberg exchange couplings
    real(dblprec), dimension(:,:,:), intent(in) :: redcoord        !< Coordinates for Heisenberg exchange couplings
    real(dblprec), dimension(:,:,:,:,:), intent(in) :: jc          !< Exchange couplings (input)
    real(dblprec), dimension(:,:,:,:,:), intent(in) :: jcD         !< Exchange couplings (input) (DLM)
    real(dblprec), dimension(:,:,:,:,:,:), intent(in) :: jc_tens   !< Tensorial exchange (SKKR) couplings (input)
    real(dblprec), dimension(:,:,:), allocatable, intent(inout) :: ncoup      !< Heisenberg exchange couplings (mounted)
    real(dblprec), dimension(:,:,:), allocatable, intent(inout) :: ncoupD     !< Heisenberg exchange couplings (mounted) (DLM)
    real(dblprec), dimension(:,:,:,:), allocatable, intent(inout) :: j_tens   !< Tensorial exchange couplings (mounted)
    ! DMI variables
    integer, intent(in) :: do_dm               !< Add Dzyaloshinskii-Moriya (DM) term to Hamiltonian (0/1)
    integer, intent(in) :: max_no_dmshells     !< Limit of number of shells for DM interactions
    integer, intent(inout) :: max_no_dmneigh   !< Calculated number of neighbours with DM interactions
    integer, dimension(:), intent(in) :: dm_nn !< Number of neighbour shells for DM
    integer, dimension(:), allocatable, intent(inout) :: dmlistsize        !< Size of neighbour list for DM
    integer, dimension(:,:), allocatable, intent(inout) :: dmlist          !< List of neighbours for DM
    real(dblprec), dimension(:,:,:), allocatable, intent(inout) :: dm_vect !< Dzyaloshinskii-Moriya exchange vector
    real(dblprec), dimension(:,:,:), intent(in) :: dm_redcoord    !< Coordinates for DM exchange couplings
    real(dblprec), dimension(:,:,:,:,:), intent(in) :: dm_inpvect !< DM couplings
    ! Anisotropy variables
    integer, intent(in) :: do_anisotropy                               !< Read anisotropy data (1/0)
    integer, dimension(:), allocatable, intent(inout) :: taniso        !< Type of anisotropy (0-2)
    integer, dimension(:), allocatable, intent(inout) :: taniso_diff   !< Type of anisotropy (0-2)
    real(dblprec) :: random_anisotropy_density                         !< Densitu of random anisotropy
    integer, dimension(:,:), intent(in) :: anisotropytype        !< Type of anisotropies (0-2)
    integer, dimension(:,:), intent(in) :: anisotropytype_diff   !< Type of anisotropies (0-2)
    real(dblprec), dimension(:,:,:), intent(in) :: anisotropy          !< Input data for anisotropies
    real(dblprec), dimension(:,:,:), intent(in) :: anisotropy_diff     !< Input data for anisotropies
    real(dblprec), dimension(:), allocatable, intent(inout) :: sb             !< Ratio between Cubic and Uniaxial anisotropy
    real(dblprec), dimension(:), allocatable, intent(inout) :: sb_diff        !< Ratio between Cubic and Uniaxial anisotropy
    real(dblprec), dimension(:,:), allocatable, intent(inout) :: eaniso       !< Unit anisotropy vector
    real(dblprec), dimension(:,:), allocatable, intent(inout) :: eaniso_diff  !< Unit anisotropy vector
    real(dblprec), dimension(:,:), allocatable, intent(inout) :: kaniso       !< Anisotropy constant
    real(dblprec), dimension(:,:), allocatable, intent(inout) :: kaniso_diff  !< Anisotropy constant
    logical :: random_anisotropy !< Distribute anisotropy randomly
    character(len=1), intent(in) :: mult_axis !< Flag to treat more than one anisotropy axis at the same time
    ! LSF variables
    integer, intent(in) :: mconf           !< LSF moment ground state conf
    integer, intent(in) :: conf_num        !< Number of configurations for LSF
    integer, dimension(:), allocatable, intent(inout) :: fs_nlistsize !< Size of first shell neighbouring list for centered atom
    integer, dimension(:,:), allocatable,intent(inout) :: fs_nlist    !< First shell Neighbouring list for centered atom
    integer, dimension(:,:), allocatable,intent(inout) :: nind        !< index of firstshell-neighbour-list corresponds to neighbour-list
    logical, intent(in) :: map_multiple        !< Allow for multiple couplings between atoms i and j
    character(len=1), intent(in) ::  do_lsf    !< (Y/N) Do LSF for MC
    character(len=1), intent(in) ::  lsf_field !< (T/L) LSF field term
    character(len=1), intent(in) ::  exc_inter !< Interpolation of Jij between FM/DLM (Y/N)
    ! BQ interactions variables
    integer, intent(in) :: do_bq               !< Add biquadratic exchange (BQ) term to Hamiltonian (0/1)
    integer, intent(inout) :: nn_bq_tot        !< Calculated number of neighbours with BQ interactions
    integer, intent(in) :: max_no_bqshells     !< Limit of number of shells for BQ interactions
    integer, dimension(:), intent(in) :: bq_nn !< Number of neighbour shells for BQ
    integer, dimension(:), allocatable, intent(inout) :: bqlistsize    !< Size of neighbour list for BQ
    integer, dimension(:,:), allocatable, intent(inout) :: bqlist      !< List of neighbours for BQ
    real(dblprec), dimension(:,:,:), intent(in) :: bq_redcoord         !< Coordinates for BQ exchange couplings
    real(dblprec), dimension(:,:,:,:), intent(in) :: jc_bq             !< Biquadratic exchange coupling (input)
    real(dblprec), dimension(:,:), allocatable, intent(inout) :: j_bq  !< Biquadratic exchange couplings (mounted)
    ! BQDM interactions variables
    integer, intent(in) :: do_biqdm               !< Add biquadratic DM (BIQDM) term to Hamiltonian (0/1)
    integer, intent(inout) :: nn_biqdm_tot        !< Calculated number of neighbours with BIQDM interactions
    integer, intent(in) :: max_no_biqdmshells     !< Limit of number of shells for BIQDM interactions
    integer, dimension(:), intent(in) :: biqdm_nn !< Number of neighbour shells for BIQDM
    integer, dimension(:), allocatable, intent(inout) :: biqdmlistsize !< Size of neighbour list for BIQDM
    integer, dimension(:,:), allocatable, intent(inout) :: biqdmlist   !< List of neighbours for BIQDM
    real(dblprec), dimension(:,:,:), intent(in) :: biqdm_redcoord      !< Coordinates for DM exchange couplings
    real(dblprec), dimension(:,:,:,:,:), intent(in) :: biqdm_inpvect   !< BIQDM couplings
    real(dblprec), dimension(:,:,:), allocatable, intent(inout) :: biqdm_vect !< BIQDM exchange vector
    ! Pseudo-Dipolar interactions variables
    integer, intent(in) :: do_pd               !< Add Pseudo-Dipolar (PD) term to Hamiltonian (0/1)
    integer, intent(inout) :: nn_pd_tot        !< Calculated number of neighbours with PD interactions
    integer, intent(in) :: max_no_pdshells     !< Limit of number of shells for PD interactions
    integer, dimension(:), intent(in) :: pd_nn !< Number of neighbour shells for PD
    integer, dimension(:), allocatable, intent(inout) :: pdlistsize !< Size of neighbour list for PD
    integer, dimension(:,:), allocatable, intent(inout) :: pdlist   !< List of neighbours for PD
    real(dblprec), dimension(:,:,:), intent(in) :: pd_redcoord      !< Coordinates for PD exchange couplings
    real(dblprec), dimension(:,:,:,:,:), intent(in) :: pd_inpvect   !< PD couplings
    real(dblprec), dimension(:,:,:), allocatable, intent(inout) :: pd_vect    !< Pseudo-Dipolar exchange vector
    ! Dipolar interactions variables
    integer, intent(in) :: do_dip !< Calculate dipole-dipole contribution (0/1)
    real(dblprec), dimension(:,:,:,:), allocatable, intent(inout) :: qdip  !< Matrix for dipole-dipole interaction
    ! Induced moments variables
    integer, dimension(NA,Nchmax), intent(in) :: ind_mom !< Indication of whether a given moment is induced (1) or fixed (0)
    integer, dimension(:), allocatable, intent(inout) :: ind_nlistsize   !< Size of the list for the induced moments
    integer, dimension(:,:), allocatable, intent(inout) :: ind_nlist !< Neighbour list between induced moments and their first permanent moments
    real(dblprec), intent(in) :: ind_tol !< Value for the tolerance between neighbouring shells
    character(len=1), intent(in) :: ind_mom_flag !< Flag to indicate that there are induced moments being considered

    ! Misc variables
    integer, intent(in) :: do_prnstruct  !< Print Hamiltonian information (0/1)
    character, intent(in) :: do_sortcoup !< Sort the entries of ncoup arrays (Y/N)
    character(len=8) :: simid            !< Name of simulation

    ! .. Local variables
    integer :: max_no_equiv !< Calculated maximum of neighbours in one shell for exchange
    real(dblprec) :: zero=0.0_dblprec


    if (do_jtensor/=1) then

       ! Setup or import map
       if (.not.prev_map(simid)) then

          !  Setup neighbor map
          write (*,'(2x,a)',advance='no') 'Set up neighbour map for exchange'
          call setup_nm(Natom, NT, NA, N1, N2, N3, C1, C2, C3, BC1, BC2, BC3, atype, Bas, &
               max_no_neigh, max_no_shells, max_no_equiv, sym, &
               nn, redcoord, nm, nmdim, &
               do_ralloy, Natom_full, acellnumb, atype_ch)
          write(*,'(a)') ' done'

          if((do_prnstruct==1.or.do_prnstruct==4)) then
             !  Print neighbor map
             write (*,'(2x,a)',advance='no') 'Print neighbour map for exchange'
             call prnge(simid, Natom, NT, NA, N1, N2, N3, atype, &
                  max_no_shells, max_no_equiv, redcoord, nn, nm, nmdim, &
                  do_ralloy, Natom_full, Nchmax, acellnumb, acellnumbrev, achtype)

             write(*,'(a)') ' done'
          end if

          ! Transform data to general structure
          write (*,'(2x,a)',advance='no') 'Mount Heisenberg Hamiltonian'
          ! Actual call to allocate the hamiltonian data
          call allocate_hamiltoniandata(Natom,conf_num,max_no_neigh,do_jtensor,do_lsf, &
               ind_mom_flag,1,lsf_field,exc_inter)
          call setup_neighbour_hamiltonian(Natom, conf_num,NT, NA, anumb, atype, &
               max_no_neigh, max_no_equiv, max_no_shells, &
               nlistsize, nn, nlist, ncoup, nm, nmdim, jc, fs_nlist, fs_nlistsize, &
               do_ralloy, Natom_full, Nchmax, &
               atype_ch, asite_ch, achem_ch, ammom_inp, 1, 1,do_sortcoup, do_lsf,nind,lsf_field,map_multiple)

          if (exc_inter=='Y') then
             call setup_neighbour_hamiltonian(Natom, conf_num,NT, NA, anumb, atype, &
                  max_no_neigh, max_no_equiv, max_no_shells, &
                  nlistsize, nn, nlist, ncoupD, nm, nmdim, jcD, fs_nlist, fs_nlistsize, &
                  do_ralloy, Natom_full, Nchmax, &
                  atype_ch, asite_ch, achem_ch, ammom_inp, 1, 1, do_sortcoup, do_lsf,nind,lsf_field,map_multiple)
          endif

          write(*,'(a)') ' done'

          ! Deallocate the large neighbour map.
          call deallocate_nm()

          if((do_prnstruct==1.or.do_prnstruct==4)) then
             write(*,'(2x,a)',advance='no') "Print exchange interaction strengths"
             call prn_exchange(Natom, max_no_neigh, nlistsize, nlist, ncoup, simid, 1)
             write(*,'(a)') ' done'
          end if

          if (ind_mom_flag=='Y') then
             write(*,'(2x,a)',advance='no') 'Set up neighbour map for induced moments'
             call induced_mapping(Natom,NT,NA,N1,N2,N3,sym,max_no_shells,nn,atype,&
                  ind_nlistsize,ind_nlist,do_sortcoup,Nchmax,do_ralloy,Natom_full,&
                  atype_ch,acellnumb,C1,C2,C3,Bas,BC1,BC2,BC3,&
                  ind_tol,redcoord,ind_mom)

             write(*,'(a)') ' done'
             if (do_prnstruct==1) then
               write(*,'(2x,a)',advance='no') 'Print neighbour map for induced moments'
               call prn_ind_exchange(Natom,max_no_neigh,ind_nlistsize,ind_nlist,simid)
               write(*,'(a)') ' done'
             endif
          endif
          if(do_prnstruct==5) then
             write(*,'(2x,a)',advance='no') "Print exchange interaction strengths (sparse format)"
             call prn_exchange_sparse(Natom, max_no_neigh, nlistsize, nlist, ncoup, simid, 1)
             write(*,'(a)') ' done'
          end if

       else ! If .out file exists

          write (*,'(2x,a)',advance='no') 'Importing exchange mapping'
          call read_exchange_getdim(Natom, max_no_neigh, simid)
          call allocate_hamiltoniandata(Natom,conf_num, max_no_neigh,do_jtensor,do_lsf,ind_mom_flag,&
          1,lsf_field,exc_inter)
          call read_exchange(Natom,conf_num, max_no_neigh, nlistsize, nlist, ncoup, simid)
          write(*,'(a)') ' done'

       end if

    else ! If tensor

       !  Setup neighbor map and exchange tensor
       write (*,'(2x,a)',advance='no') 'Set up neighbour map for exchange'
       call setup_nm(Natom, NT, NA, N1, N2, N3, C1, C2, C3, BC1, BC2, BC3, atype, Bas, &
            max_no_neigh, max_no_shells, max_no_equiv, sym, &
            nn, redcoord, nm, nmdim, &
            do_ralloy, Natom_full, acellnumb, atype_ch)
       write(*,'(a)') ' done'

       if (do_prnstruct==1.or.do_prnstruct==4) then
          !  Print neighbor map
          write (*,'(2x,a)',advance='no') 'Print neighbour map for exchange'
          call prnge(simid, Natom, NT, NA, N1, N2, N3, atype, &
               max_no_shells, max_no_equiv, redcoord, nn, nm, nmdim, &
               do_ralloy, Natom_full, Nchmax, acellnumb, acellnumbrev, achtype)

          write(*,'(a)') ' done'
       end if

       ! Transform data to general structure
       write (*,'(2x,a)',advance='no') 'Mount Heisenberg Hamiltonian in tensorial (SKKR) form'
       call allocate_hamiltoniandata(Natom,conf_num,max_no_neigh,do_jtensor,do_lsf, ind_mom_flag,&
       1,lsf_field,exc_inter)
       call setup_neighbour_hamiltonian(Natom, conf_num,NT, NA, anumb, atype, &
            max_no_neigh, max_no_equiv, max_no_shells, &
            nlistsize, nn, nlist, j_tens,nm, nmdim, jc_tens, fs_nlist, fs_nlistsize, &
            do_ralloy, Natom_full, Nchmax, &
            atype_ch, asite_ch, achem_ch, ammom_inp, 9,1,do_sortcoup,do_lsf, nind,lsf_field,map_multiple)
       write(*,'(a)') ' done'

       ! Deallocate the large neighbour map.
       call deallocate_nm()

       if(do_prnstruct==1.or.do_prnstruct==4) then
          write(*,'(2x,a)',advance='no') "Print exchange interaction strenghts"
          call prn_exchange(Natom, max_no_neigh, nlistsize, nlist, j_tens, simid, 9)
          write(*,'(a)') ' done'
       end if

    end if

    ! Anisotropies
    write(*,'(2x,a)',advance='no') "Set up anisotropies"
    call allocate_anisotropies(Natom,mult_axis,1)
    call setup_anisotropies(Natom, NA ,anumb, anisotropytype, &
          taniso, eaniso, kaniso, sb, anisotropy, &
          mult_axis,anisotropytype_diff,taniso_diff, eaniso_diff, &
          kaniso_diff, sb_diff, anisotropy_diff, &
          random_anisotropy, random_anisotropy_density, &
          do_ralloy, Natom_full, Nchmax, achem_ch, ammom_inp,mconf)
    write(*,'(a)') ' done'

    ! Print anisotropies
    if(do_prnstruct==1 .and. do_anisotropy == 1) then
       write(*,'(2x,a)',advance='no') "Print anisotropies"
       call prn_anisotropy(Natom, NA, anisotropytype, taniso, eaniso, kaniso, simid, Nchmax)
       write(*,'(a)') ' done'
    end if

    !  Set up dipole-tensor
    if(do_dip==1) then
       write(*,'(2x,a)',advance='no') "Set up dipole-dipole matrix"
       call allocate_dipole(Natom,1)
       call setup_qdip(Natom, coord, alat, Qdip)
       write(*,'(a)') '  done'
    end if

    ! Calculate mean field estimate of Tc
    if(do_jtensor/=1.and.int(N1*N2*N3).ne.1) then
       write(*,'(2x,a)',advance='no') "Calculate mean-field estimate of Tc: "
       call estimate_tc_mfa(Natom,Nchmax,conf_num,NA,atype, anumb,Natom_full,atype_ch,&
       asite_ch,achem_ch,nlistsize, nlist, max_no_neigh,ncoup,ammom_inp,mconf,do_ralloy)
    end if

    if(do_dm==1) then
       ! Allocate and mount DM Hamiltonian
       write (*,'(2x,a)',advance='no') 'Set up neighbour map for Dzyaloshinskii-Moriya exchange'
       call setup_nm(Natom, NT, NA, N1, N2, N3, C1, C2, C3, BC1, BC2, BC3, atype, Bas, &
            max_no_dmneigh, max_no_dmshells, max_no_equiv, 0, dm_nn, dm_redcoord, nm, nmdim, &
            do_ralloy, Natom_full, acellnumb, atype_ch)
       write(*,'(a)') ' done'

      if (do_prnstruct==1) then
         write (*,'(2x,a)',advance='no') 'Print neighbour map for DMI'
         call dmprnge(simid, Natom, NT, NA, N1, N2, N3, atype, &
              max_no_dmshells, max_no_equiv, dm_redcoord, dm_nn, nm, nmdim, &
              do_ralloy, Natom_full, Nchmax, acellnumb, acellnumbrev, achtype)

         write(*,'(a)') ' done'
      endif

       call allocate_dmhamiltoniandata(Natom,max_no_dmneigh,1)

       ! Transform data to general structure
       write (*,'(2x,a)',advance='no') 'Mount Dzyaloshinskii-Moriya Hamiltonian'
       call setup_neighbour_hamiltonian(Natom, conf_num,NT, NA, anumb, atype, &
            max_no_dmneigh, max_no_equiv, max_no_dmshells, &
            dmlistsize, dm_nn, dmlist, dm_vect, nm, nmdim, dm_inpvect, fs_nlist, fs_nlistsize, &
            do_ralloy, Natom_full, Nchmax, &
            atype_ch, asite_ch, achem_ch, ammom_inp, 3, 1,do_sortcoup,do_lsf, nind,lsf_field,map_multiple)
       write(*,'(a)') ' done'

       ! Deallocate the large neighbour map.
       call deallocate_nm()

       ! Print DM interactions
       if((do_prnstruct==1.or.do_prnstruct==4)) then 
          write(*,'(2x,a)',advance='no') "Print Dzyaloshinskii-Moriya interactions"
          call prn_dmcoup(Natom, max_no_dmneigh, dmlistsize, dmlist, dm_vect, simid)
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
       call setup_nm(Natom, NT, NA, N1, N2, N3, C1, C2, C3, BC1, BC2, BC3, atype, Bas, &
            nn_pd_tot, max_no_pdshells, max_no_equiv, sym, &
            pd_nn, pd_redcoord, nm, nmdim, &
            do_ralloy, Natom_full, acellnumb, atype_ch)
       write(*,'(a)') ' done'

       call allocate_pdhamiltoniandata(Natom,nn_pd_tot,1)

       ! Transform data to general structure
       write (*,'(2x,a)',advance='no') 'Mount Pseudo-Dipolar Hamiltonian'
       call setup_neighbour_hamiltonian(Natom, conf_num,NT, NA, anumb, atype, &
            nn_pd_tot, max_no_equiv, max_no_pdshells, &
            pdlistsize, pd_nn, pdlist, pd_vect, nm, nmdim, pd_inpvect, fs_nlist, fs_nlistsize, &
            do_ralloy, Natom_full, Nchmax, &
            atype_ch, asite_ch, achem_ch, ammom_inp, 6, 1,do_sortcoup, do_lsf,nind,lsf_field,map_multiple)
       write(*,'(a)') ' done'

       ! Deallocate the large neighbour map.
       call deallocate_nm()

       ! Print PD interactions
       if(do_prnstruct==1.or.do_prnstruct==4) then
          write(*,'(2x,a)',advance='no') "Print Pseudo-Dipolar interactions"
          call prn_pdcoup(Natom, nn_pd_tot, pdlistsize, pdlist, pd_vect, simid)
          write(*,'(a)') ' done'
       end if

    end if

    if(do_biqdm==1) then

      ! Allocate and mount BIQDM Hamiltonian
      write (*,'(2x,a)',advance='no') 'Set up neighbour map for BIQDM exchange'
      call setup_nm(Natom, NT, NA, N1, N2, N3, C1, C2, C3, BC1, BC2, BC3, atype, Bas, &
           nn_biqdm_tot, max_no_biqdmshells, max_no_equiv, 0, &
           biqdm_nn, biqdm_redcoord, nm, nmdim, &
           do_ralloy, Natom_full, acellnumb, atype_ch)
      write(*,'(a)') ' done'

      call allocate_biqdmhamiltoniandata(Natom,nn_biqdm_tot,1)
      ! Transform data to general structure
      write (*,'(2x,a)',advance='no') 'Mount BIQDM Hamiltonian'

      call setup_neighbour_hamiltonian(Natom, conf_num,NT, NA, anumb, atype, &
           nn_biqdm_tot, max_no_equiv, max_no_biqdmshells, &
           biqdmlistsize, biqdm_nn, biqdmlist, biqdm_vect,nm, nmdim, biqdm_inpvect, fs_nlist, fs_nlistsize, &
           do_ralloy, Natom_full, Nchmax, &
           atype_ch, asite_ch, achem_ch, ammom_inp, 1, 2,do_sortcoup, do_lsf,nind,lsf_field,map_multiple)
      write(*,'(a)') ' done'

      ! Deallocate the large neighbour map.
      call deallocate_nm()

      ! Print BIQDM interactions
      if(do_prnstruct==1.or.do_prnstruct==4) then
        write(*,'(2x,a)',advance='no') "Print BIQDM interactions"
        call prn_biqdmcoup(Natom, nn_biqdm_tot, biqdmlistsize, biqdmlist, biqdm_vect, simid)
        write(*,'(a)') ' done'
      end if

    end if

    if(do_bq==1) then

       ! Allocate and mount BQ Hamiltonian
       write (*,'(2x,a)',advance='no') 'Set up neighbour map for biquadratic exchange'
       call setup_nm(Natom, NT, NA, N1, N2, N3, C1, C2, C3, BC1, BC2, BC3, atype, Bas, &
            nn_bq_tot, max_no_bqshells, max_no_equiv, sym, &
            bq_nn, bq_redcoord, nm, nmdim, &
            do_ralloy, Natom_full, acellnumb, atype_ch)
       write(*,'(a)') ' done'

       call allocate_bqhamiltoniandata(Natom,nn_bq_tot,1)

       ! Transform data to general structure
       write (*,'(2x,a)',advance='no') 'Mount biquadratic exchange Hamiltonian'
       call setup_neighbour_hamiltonian(Natom, conf_num,NT, NA, anumb, atype, &
            nn_bq_tot, max_no_equiv, max_no_bqshells, &
            bqlistsize, bq_nn, bqlist, j_bq, nm, nmdim, jc_bq, fs_nlist, fs_nlistsize, &
            do_ralloy, Natom_full, Nchmax, &
            atype_ch, asite_ch, achem_ch, ammom_inp, 1, 2,do_sortcoup, do_lsf,nind,lsf_field,map_multiple)
       write(*,'(a)') ' done'

       ! Deallocate the large neighbour map.
       call deallocate_nm()

       ! Print BQ interactions
       if(do_prnstruct==1.or.do_prnstruct==4) then
          write(*,'(2x,a)',advance='no') "Print biquadratic exchange interactions"
          call prn_bqcoup(Natom, nn_bq_tot, bqlistsize, bqlist, j_bq, simid)
          write(*,'(a)') ' done'
       end if
     endif

       ! Deallocate auxiliarly data of the Ewald sum

  contains


    !> Evaluates if exchange map already exists
    logical function prev_map(simid)
      !
      implicit none
      !
      character(len=8),intent(in) :: simid !< Name of simulation
      character(len=20) :: filn

      !.. Executable statements
      write (filn,'(''struct1.'',a8,''.out'')') simid
      inquire(file=filn,exist=prev_map)
      return
    end function prev_map

  end subroutine setup_hamiltonian


  !> Sets up single ion anisotropies
  subroutine setup_anisotropies(Natom, NA ,anumb, anisotropytype, &
             taniso, eaniso, kaniso, sb, anisotropy, &
             mult_axis,anisotropytype_diff,taniso_diff, eaniso_diff, &
             kaniso_diff, sb_diff, anisotropy_diff, &
             random_anisotropy, random_anisotropy_density, &
             do_ralloy, Natom_full, Nchmax, achem_ch, ammom_inp,mconf)

    use Constants
    use RandomNumbers, only: rng_uniform
    !
    implicit none
    !
    integer, intent(in) :: Natom                                      !< Number of atoms in system
    integer, intent(in) :: NA                                         !< Number of atoms in one cell
    integer, dimension(Natom), intent(in) :: anumb                    !< Atom number in cell
    integer, intent(in) :: Nchmax                                     !< Max number of chemical components on each site in cell
    integer, dimension(:,:), intent(in) :: anisotropytype       !< Type of anisotropies (0-2)
    integer, dimension(:,:), intent(in) :: anisotropytype_diff  !< Type of anisotropies (0-2)
    real(dblprec), dimension(:,:,:), intent(in) :: anisotropy         !< Input data for anisotropies
    real(dblprec), dimension(:,:,:), intent(in) :: anisotropy_diff    !< Input data for anisotropies
    logical :: random_anisotropy                                      !< Distribute anisotropy randomly
    real(dblprec) :: random_anisotropy_density                        !< Densitu of random anisotropy
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
    real(dblprec), dimension(:,:,:), intent(in) :: ammom_inp            !< Magnetic moment directions from input (for alloys)
    character(len=1), intent(in) :: mult_axis
    integer, intent(in) :: mconf                                       !< LSF ground state configuration

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

  !> Mounts the actual Heisenberg Hamiltonian (or DM Hamiltonian)
  subroutine setup_neighbour_hamiltonian(Natom,conf_num, NT, NA, anumb, atype, &
             max_no_neigh, max_no_equiv, max_no_shells, &
             nlistsize, nn, nlist, ncoup,nm, nmdim, xc, fs_nlist, fs_nlistsize, &
             do_ralloy, Natom_full, Nchmax, &
             atype_ch, asite_ch, achem_ch, ammom_inp, hdim, lexp,do_sortcoup,do_lsf, &
             nind,lsf_field, map_multiple)
    !
    use Constants
    use Sorting, only : MergeSortIR
    !
    implicit none
    !
    integer, intent(in) :: Natom !< Number of atoms in system
    integer, intent(in) :: conf_num !< Number of configurations for LSF
    integer, intent(in) :: NT  !< Number of types of atoms
    integer, intent(in) :: NA  !< Number of atoms in one cell
    integer, dimension(Natom), intent(in) :: anumb !< Atom number in cell
    integer, dimension(Natom), intent(in) :: atype !< Type of atom
    integer, dimension(Natom), intent(out):: nlistsize !< Size of neighbour list for Heisenberg exchange couplings
    integer, dimension(NT), intent(in):: nn !< Number of neighbour shells
    integer, intent(in) :: max_no_neigh !< Calculated maximum of neighbours for exchange
    integer, intent(in) :: max_no_equiv !< Calculated maximum of neighbours in one shell for exchange
    integer, intent(in) :: max_no_shells !< Calculated maximum of shells for exchange
    integer, dimension(max_no_neigh,Natom), intent(out) :: nlist !< Neighbour list for Heisenberg exchange couplings
    integer, intent(in) :: hdim  !< Number of elements in Hamiltonian element (scalar or vector)
    integer, intent(in) :: lexp  !< Order of the spin-interaction. Needed for proper rescaling (1=linear,2=quadratic)
    character :: do_sortcoup !< Sort the entries of ncoup arrays (Y/N)
    real(dblprec), dimension(hdim,max_no_neigh,Natom,conf_num), intent(out) :: ncoup !< Heisenberg exchange couplings
    integer, dimension(Natom,max_no_shells,max_no_equiv), intent(in) :: nm !< Neighbour map
    integer, dimension(max_no_shells,Natom), intent(in) :: nmdim !< Dimension of neighbour map
    integer, dimension(:,:), intent(out) :: fs_nlist !< First shell Neighbouring list for centered atom
    integer, dimension(:), intent(out) :: fs_nlistsize  !< Size of first shell neighbouring list for centered atom
    integer, intent(in) :: Nchmax !< Max number of chemical components on each site in cell
    real(dblprec), dimension(hdim,NT,max_no_shells,Nchmax,max(Nchmax,NT),conf_num), intent(in) :: xc !< Exchange couplings
    integer, intent(in) :: do_ralloy  !< Random alloy simulation (0/1)
    integer, intent(in) :: Natom_full !< Number of atoms for full system (=Natom if not dilute)
    integer, dimension(Natom_full), intent(in) :: atype_ch !< Actual site of atom for dilute system
    integer, dimension(Natom_full), intent(in) :: asite_ch !< Actual site of atom for dilute system
    integer, dimension(Natom_full), intent(in) :: achem_ch !< Chemical type of atoms (reduced list)
    real(dblprec), dimension(NA,Nchmax,conf_num), intent(in) :: ammom_inp !< Magnetic moment directions from input (for alloys)
    character(len=1), intent(in) :: do_lsf !< (Y/N) Do LSF for MC
    integer, allocatable, dimension(:,:)  :: nind !< index of firstshell-neighbour-list corresponds to neighbour-list
    character(len=1), intent(in) :: lsf_field !< (T/L) LSF field term
    logical, intent(in) :: map_multiple    !< Allow for multiple couplings between atoms i and j

    !
    integer :: i, j, k, l,iconf, i_stat
    integer :: jatom, jchtype
    real(dblprec) :: fc,fc2
    real(dblprec), dimension(hdim) :: tempcoup
    integer :: tempn
    integer :: ncount,ncounter
    logical :: exis

    ncoup = 0.0d0
    ! Factors for mRy energy conversion
    fc = mry/mub
    fc2 = 2.0d0*mry/mub

    ! Loop over atoms
    !$omp parallel do default(shared) private(i,k,j,ncount,ncounter,exis,l,jatom,jchtype,iconf)
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
                   do iconf=1,conf_num
                      if (do_ralloy==0) then
                         if(abs(ammom_inp(anumb(i),1,iconf)*ammom_inp(anumb(nm(i,k,j)),1,iconf))<1e-6) then
                                ncoup(:,ncount,i,iconf)=0.d0
                         else
                                ncoup(:,ncount,i,iconf) = xc(:,atype(i),k,1,1,iconf) * fc2 &
                              / ammom_inp(anumb(i),1,iconf)**lexp / ammom_inp(anumb(nm(i,k,j)),1,iconf)**lexp
                         endif
                      else
                         !
                         jatom = nm(i,k,j)
                         jchtype = atype_ch(jatom)
                         if(abs(ammom_inp(asite_ch(i),achem_ch(i),iconf)*ammom_inp(asite_ch(jatom),achem_ch(jatom),iconf))<1e-6) then
                             ncoup(:,ncount,i,iconf)=0.d0
                         else
                             ncoup(:,ncount,i,iconf) = xc(:,atype_ch(i),k,achem_ch(i),achem_ch(jatom),iconf) * &
                              fc2 / ammom_inp(asite_ch(i),achem_ch(i),iconf)**lexp  / ammom_inp(asite_ch(jatom),achem_ch(jatom),iconf)**lexp
                         endif
                      end if
                   end do
                   if (k==1 .and. do_lsf=='Y' .and. lsf_field=='L') then  !< conditional adding fist-shell nlist to fs_nlist
                       fs_nlist(ncount,i) = nlist(ncount,i)
                       ncounter = ncounter + 1
                   end if
                   ncount = ncount + 1
                end if
             end if
          end do
       end do
       nlistsize(i) = ncount-1
       if (do_lsf=='Y' .and. lsf_field=='L') fs_nlistsize(i) = ncounter-1
    end do
    !omp end parallel do

    if(do_sortcoup == 'Y') then
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
  subroutine read_exchange_getdim(Natom, max_no_neigh,simid)
    !

    !.. Implicit declarations
    implicit none

    integer, intent(in) :: Natom !< Number of atoms in system
    integer, intent(out) :: max_no_neigh !< Calculated maximum of neighbours for exchange
    character(len=8),intent(in) :: simid !< Name of simulation

    !.. Local variables
    integer :: i,j
    integer :: ntest
    character(len=20) :: filn

    !.. Executable statements
    write (filn,'(''struct1.'',a8,''.out'')') simid
    open(ifileno, file=filn)

    ! Find max no. neighbours
    max_no_neigh=1
    read (ifileno,*)
    do i=1,Natom
       read (ifileno,*)
       read (ifileno,10001) j,ntest
       do j=1,2*ceiling(ntest/5.0d0)
          read (ifileno,*)
       end do
       max_no_neigh=max(max_no_neigh,ntest)
    end do
    close(ifileno)

10001 format (5x,i8,13x,i7)

  end subroutine read_exchange_getdim


  !> Reads the Heisenberg Hamiltonian from file
  subroutine read_exchange(Natom, conf_num,max_no_neigh, nlistsize, nlist, ncoup, simid)
    !
    use Constants, only : mub,mry
    !.. Implicit declarations
    implicit none

    integer, intent(in) :: Natom !< Number of atoms in system
    integer, intent(in) :: conf_num !< Number of configurations for LSF
    integer, intent(in) :: max_no_neigh !< Calculated maximum of neighbours for exchange
    integer, dimension(Natom), intent(out) :: nlistsize !< Size of neighbour list for Heisenberg exchange couplings
    integer, dimension(max_no_neigh, Natom), intent(out) :: nlist !< Neighbour list for Heisenberg exchange couplings
    real(dblprec), dimension(max_no_neigh, Natom,conf_num), intent(out) :: ncoup !< Heisenberg exchange couplings
    character(len=8),intent(in) :: simid !< Name of simulation

    !.. Local variables
    integer :: i,j
    character(len=20) :: filn

    !.. Executable statements
    if (conf_num>1) then
       write(*,*) 'Restart does not work with LSF right now'
       stop
    endif


    write (filn,'(''struct1.'',a8,''.out'')') simid
    open(ifileno, file=filn)

    ncoup=0.0d0
    read (ifileno,*)
    do i=1,Natom
       ! Exchange
       read (ifileno,*)
       read (ifileno,10001) j,nlistsize(i)
       read (ifileno,10002) nlist(1:nlistsize(i),i)
       read (ifileno,10003) ncoup(1:nlistsize(i),i,1)
    end do
    close(ifileno)
    ncoup=ncoup*mry/mub
10001 format (5x,i8,13x,i7)
10002 format (13x,5i6)
10003 format (5es16.8)

  end subroutine read_exchange


  !> Set up the dipole-dipole interaction matrix
  !! @todo Change to correct scale, needs length information
  !! @todo Consider Ewald summation
  subroutine setup_qdip(Natom, coord, alat, Qdip)
    !
    use Constants
    !
    implicit none
    !
    integer, intent(in) :: Natom !< Number of atoms in system
    real(dblprec),dimension(3,Natom), intent(in) :: coord !< Coordinates of atoms
    real(dblprec) :: alat !< Lattice parameter
    real(dblprec),dimension(3,3,Natom,Natom), intent(out):: Qdip !< Matrix for dipole-dipole interaction

    ! ...Local Variables...
    integer :: i,j,mu,nu
    real(dblprec) :: fac,R2,Rm5,Rm53,R,Rm3
    real(dblprec),dimension(3) :: Rij
    !
    ! fac determines the strength of the dipole interaction
    fac= (1.0d-7/alat**3)*mub
    print *,alat,mub
    print *,fac
    !
    !$omp parallel do private(i,j,Rij,R2,R,Rm3,Rm5,Rm53,mu,nu)
    do i=1,Natom
       do j=1,Natom
          Rij=(coord(:,i)-coord(:,j)) !*alat
          R2=Rij(1)**2+Rij(2)**2+Rij(3)**2
          R=sqrt(R2)
          Rm5=R2**(-2.5)*fac
          Rm53=3.0d0*Rm5
          do mu=1,3
             do nu=1,3
                Qdip(nu,mu,j,i)=Rm53*Rij(nu)*Rij(mu)
             end do
             Qdip(mu,mu,j,i)=Qdip(mu,mu,j,i)-Rm5*R2
          end do
       end do
       Qdip(:,:,i,i)=0.0d0
    end do
    !$omp end parallel do
    !
    return
    !
  end subroutine setup_qdip


  !> Calculate a mean field estimate of the critical temperature
  subroutine estimate_tc_mfa(Natom,Nchmax,conf_num,NA, atype, anumb, Natom_full,&
             atype_ch,asite_ch,achem_ch,nlistsize, nlist, max_no_neigh, ncoup,ammom_inp,&
             mconf,do_ralloy)
    !
    use Constants
    !
    implicit none
    !
    integer, intent(in) :: Natom !< Number of atoms in system
    integer, intent(in) :: Nchmax !< Number of chemical components in system
    integer, intent(in) :: conf_num !< Number of configurations for LSF
    integer, intent(in) :: NA  !< Number of atoms in one cell
    integer, dimension(Natom), intent(in) :: atype !< Type of atom
    integer, dimension(Natom), intent(in) :: anumb !< Atom number in cell
    integer, intent(in) :: Natom_full !< Number of atoms for full system (=Natom if not dilute)
    integer, dimension(Natom_full), intent(in) :: atype_ch !< Actual site of atom for dilute system
    integer, dimension(Natom_full), intent(in) :: asite_ch !< Actual site of atom for dilute system
    integer, dimension(Natom_full), intent(in) :: achem_ch !< Chemical type of atoms (reduced list)
    integer, dimension(Natom), intent(in) :: nlistsize !< Size of neighbour list for Heisenberg exchange couplings
    integer, intent(in) :: max_no_neigh !< Calculated maximum of neighbours for exchange
    integer, dimension(max_no_neigh,Natom), intent(in) :: nlist !< Neighbour list for Heisenberg exchange couplings
    real(dblprec), dimension(max_no_neigh,Natom,conf_num), intent(in) :: ncoup !< Heisenberg exchange couplings
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
    integer :: i,j,k, i_stat, i_all
    integer :: lwork, info
    !
    ! Allocate J0 matrix
    allocate(jmat(na,na),stat=i_stat)
    call memocc(i_stat,product(shape(jmat))*kind(jmat),'jmat','estimate_tc_mfa')
    jmat=0.0d0

    lwork=3*na
    allocate(work(lwork))
    allocate(ctemp(na))
    allocate(etemp(na,na))
    allocate(rwork(lwork))
    allocate(cwres(NA))

    fc2 = 2.0d0*mry/mub
    ! Fill J0 matrix

    if (do_ralloy==0) then
       do i=1,na
          do j=1,nlistsize(i)
             k=atype(nlist(j,i))
             jmat(i,k)=jmat(i,k)+ncoup(j,i,mconf)*ammom_inp(anumb(i),1,mconf)*ammom_inp(anumb(nlist(j,i)),1,mconf)
          end do
       end do
    else
       do i=1,na
          do j=1,nlistsize(i)
             k=atype_ch(nlist(j,i))
             jmat(i,k)=jmat(i,k)+ncoup(j,i,mconf)*ammom_inp(asite_ch(i),achem_ch(i),mconf)*ammom_inp(asite_ch(nlist(j,i)),achem_ch(nlist(j,i)),mconf)
          end do
       end do

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

    ! Write result to stdout
    write(*,'(f7.1,a)') (2.0d0*em/3.0d0/k_bolt_ev/1000),' K.'
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


end module HamiltonianInit
