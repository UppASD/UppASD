
!> @copyright
!! Copyright (C) 2008-2018 UppASD group
!! This file is distributed under the terms of the
!! GNU General Public License. 
!! See http://www.gnu.org/copyleft/gpl.txt

module SpinIce

  use Parameters
  use Constants
  use SpinIceData
  use RandomNumbers, only: rng_uniform,rng_gaussian, use_vsl
   use ErrorHandling

  implicit none

  integer, dimension(:), allocatable :: nmdimt_ver_ver ! Temporary storage of dimensions for vertex-vertex neighbour map
  integer, dimension(:), allocatable :: nmdimt_ver_atom ! Temporary storage of dimension for atom-vertex neighbour map
  real(dblprec), dimension(:,:,:), allocatable :: sym_mats !< Symmetry operation matrices
  integer :: nsym !< Number of symmetry operations

  private :: nmdimt_ver_ver,nmdimt_ver_atom

  public :: ice_rule, select_cluster, setup_neighbour_list
  public :: setup_ver_nm, setup_ver_atom_nm, mc_update_spinice

contains

!> Driver for spin ice Monte Carlo

 subroutine mc_update_spinice(Natom, Mensemble, max_no_neigh, conf_num, ncoup, ncoupD, nlist, nlistsize, &
            do_dm, max_no_dmneigh, dm_vect, dmlist, dmlistsize, do_pd, nn_pd_tot, pd_vect, pdlist, pdlistsize, &
            do_biqdm, nn_biqdm_tot, biqdm_vect, biqdmlist, biqdmlistsize, do_bq, nn_bq_tot, j_bq, bqlist, bqlistsize, &
            taniso, eaniso, kaniso,sb,emomM, emom, mmom, iflip_a, extfield, &
            mult_axis,taniso_diff, eaniso_diff, kaniso_diff,sb_diff, &
            do_dip, Qdip,exc_inter, temperature,temprescale,ind_nlistsize,ind_nlist,sus_ind,ind_mom_flag)

    use montecarlo_common
    !
    implicit none

    integer, intent(in) :: Natom !< Number of atoms in system
    integer, intent(in) :: Mensemble !< Number of ensembles
    integer, intent(in) :: max_no_neigh !< Calculated maximum of neighbours for exchange
    integer, intent(in) :: conf_num  !< number of configurations for LSF
    real(dblprec), dimension(max_no_neigh,Natom,conf_num), intent(in) :: ncoup !< Heisenberg exchange couplings
    real(dblprec), dimension(max_no_neigh,Natom,conf_num), intent(in) :: ncoupD !< Heisenberg exchange couplings (DLM)
    integer, dimension(max_no_neigh,Natom), intent(in) :: nlist !< Neighbour list for Heisenberg exchange couplings
    integer, dimension(Natom),intent(in) :: nlistsize !< Size of neighbour list for Heisenberg exchange couplings
    integer, intent(in) :: do_dm   !< Add Dzyaloshinskii-Moriya (DM) term to Hamiltonian (0/1)
    integer, intent(in) :: max_no_dmneigh !< Calculated number of neighbours with DM interactions
    real(dblprec), dimension(3,max_no_dmneigh,Natom), intent(in) :: dm_vect !< Dzyaloshinskii-Moriya exchange vector
    integer, dimension(max_no_dmneigh,Natom), intent(in) :: dmlist   !< List of neighbours for DM
    integer, dimension(Natom),intent(in) :: dmlistsize !< Size of neighbour list for DM
    integer, intent(in) :: do_pd   !< Add Pseudo-Dipolar (PD) term to Hamiltonian (0/1)
    integer, intent(in) :: nn_pd_tot !< Calculated number of neighbours with PD interactions
    real(dblprec), dimension(6,nn_pd_tot,Natom), intent(in) :: pd_vect !< Pseudo-Dipolar exchange vector
    integer, dimension(nn_pd_tot,Natom), intent(in) :: pdlist   !< List of neighbours for PD
    integer, dimension(Natom),intent(in) :: pdlistsize !< Size of neighbour list for PD
    integer, intent(in) :: do_biqdm   !< Add biquadratic DM (BIQDM) term to Hamiltonian (0/1)
    integer, intent(in) :: nn_biqdm_tot !< Calculated number of neighbours with BIQDM interactions
    real(dblprec), dimension(1,nn_biqdm_tot,Natom), intent(in) :: biqdm_vect !< BIQDM exchange coupling
    integer, dimension(nn_biqdm_tot,Natom), intent(in) :: biqdmlist   !< List of neighbours for BIQDM
    integer, dimension(Natom),intent(in) :: biqdmlistsize !< Size of neighbour list for BIQDM
    integer, intent(in) :: do_bq   !< Add biquadratic exchange (BQ) term to Hamiltonian (0/1)
    integer, intent(in) :: nn_bq_tot !< Calculated number of neighbours with BQ interactions
    real(dblprec), dimension(nn_bq_tot,Natom), intent(in) :: j_bq !< Biquadratic exchange couplings
    integer, dimension(nn_bq_tot,Natom), intent(in) :: bqlist   !< List of neighbours for BQ
    integer, dimension(Natom),intent(in) :: bqlistsize !< Size of neighbour list for BQ
    integer, dimension(Natom),intent(in) :: taniso !< Type of anisotropy (0-2)
    integer, dimension(Natom),intent(in) :: taniso_diff !< Type of anisotropy (0-2)
    real(dblprec), dimension(3,Natom), intent(in) :: eaniso !< Unit anisotropy vector
    real(dblprec), dimension(3,Natom), intent(in) :: eaniso_diff !< Unit anisotropy vector
    real(dblprec), dimension(2,Natom), intent(in) :: kaniso !< Anisotropy constant
    real(dblprec), dimension(2,Natom), intent(in) :: kaniso_diff !< Anisotropy constant
    real(dblprec), dimension(Natom), intent(in) :: sb !< Ratio between anisotropy constants
    real(dblprec), dimension(Natom), intent(in) :: sb_diff !< Ratio between anisotropy constants
    real(dblprec), dimension(3,Natom,Mensemble), intent(inout) :: emomM  !< Current magnetic moment vector
    real(dblprec), dimension(3,Natom,Mensemble), intent(inout) :: emom   !< Current unit moment vector
    real(dblprec), dimension(Natom,Mensemble), intent(inout) :: mmom !< Magnitude of magnetic moments
    real(dblprec), intent(in) :: temperature !< Temperature
    real(dblprec), intent(in) :: temprescale !< Temperature rescaling according to QHB
    real(dblprec), dimension(3), intent(in) :: extfield !< External magnetic field
    integer, intent(in) :: do_dip  !<  Calculate dipole-dipole contribution (0/1)
    real(dblprec), dimension(3,3,Natom,Natom), intent(in), optional :: Qdip !< Matrix for dipole-dipole interaction
    character(len=1), intent(in) :: mult_axis
    character(len=1), intent(in)  ::  exc_inter           !< Exchange interpolation between FM/DLM (Y/N)
    integer,dimension(natom),intent(in) :: iflip_a !< Flipping pattern
    integer, dimension(Natom),intent(in) :: ind_nlistsize !< Size of neighbour list for Heisenberg exchange couplings
    integer, dimension(max_no_neigh,Natom), intent(in) :: ind_nlist !< Neighbour list for Heisenberg exchange couplings
    real(dblprec), dimension(Natom), intent(in) :: sus_ind
    character(len=1), intent(in) :: ind_mom_flag
 !.. Local scalars
    integer :: i, k, loop_len, ii, ip
    integer :: iatom
    real(dblprec) :: de, delt  !< Energy difference

    real(dblprec) :: ve, ice_count,ice_count_up,ice_count_down !< Trial vertex energy

 !.. Local arrays
    real(dblprec), dimension(3) :: newmom !< New trial moment

    real(dblprec),dimension(natom) :: flipprob_a ,newmmom_a
    real(dblprec),dimension(3,natom) :: newmom_a

    real(dblprec), dimension(3,Natom,Mensemble) :: emomM_copy  !< Copy of current magnetic moment vector

    integer, dimension(Nvertex) :: Ncount_ver
    integer, dimension(Natom) :: Ncount_atom

    call ErrorHandling_missing('Spin-ice')

 end subroutine mc_update_spinice


  ! This routine creates a list of spins that will be flipped in acordance with the spin-ice rules
  ! The list must be designed in a smart way such that one does not have to go trough all the spins everytime one selects a random one if not the process is of the order of N2 which is not efficient
  subroutine select_cluster(Natom, Mensemble,emom, Ncount_atom, loop_len, k)

    implicit none

    integer, intent(in) :: k ! Current ensemble
    integer, intent(in) :: Natom
    integer, intent(in) :: Mensemble
    integer, dimension(Natom), intent(inout) :: Ncount_atom

    real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: emom

    integer, intent(out) :: loop_len ! length of spin ice loop

    !AB
    logical :: loop_closed, atom_selected, failure, total_failure

    integer :: start_vertex, current_vertex, next_vertex, prev_vertex, flip_atom
    integer :: flip_atom_list_num, prev_flip_atom,prev_flip_atom_num_list
    !
    integer, dimension(:), allocatable :: visited_vertices, snake, atomsnake ! This are the vertices that are visited and form the cluster
    integer :: snaketail, atomlen, i, snakehead  ! This are the indexes used for the cluster

    loop_len=0
    call ErrorHandling_missing('Spin-ice')

  end subroutine select_cluster

  ! This subroutine selects the first vertex by checking if it follows the ice rules and if it has neighbours that follows the ice rules
  subroutine select_first_atom(Natom, Nvertex, Mensemble, max_no_equiv_ver, max_no_equiv_ver_atom, nlist_ver, &
             nlist_ver_atom, nlistsize_ver, nlistsize_ver_atom, start_vertex, current_vertex, flip_atom, flip_atom_num_list,atom_selected,&
             failure,en,emom,vertex_hit,vert_ice_coord)

    implicit none

    integer, intent(in) :: Natom
    integer, intent(in) :: Nvertex
    integer, intent(in) :: Mensemble
    integer, intent(in) :: en ! Current ensemble
    integer, intent(in) :: max_no_equiv_ver
    integer, intent(in) :: max_no_equiv_ver_atom
    integer, dimension(max_no_equiv_ver,Nvertex), intent(in) :: nlist_ver ! Neighbour list for vertex vertex or vertex atoms
    integer, dimension(Nvertex), intent(in):: nlistsize_ver ! Size of neighbour list for vertex vertex or vertex atoms
    integer, dimension(max_no_equiv_ver_atom,Nvertex) , intent(in) :: nlist_ver_atom
    integer, dimension(Nvertex), intent(in) :: nlistsize_ver_atom

    integer, intent(out) :: start_vertex
    integer, intent(out) :: current_vertex
    integer, intent(out) :: flip_atom
    integer, intent(out) :: flip_atom_num_list

    logical, intent(inout) :: atom_selected ! Logical variable to see if there is an atom selected
    logical, intent(out) :: failure         ! Logical variable to see if there was a failure in the algorithm
    real(dblprec), dimension(3,Nvertex,max_no_equiv_ver_atom), intent(in) :: vert_ice_coord
    real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: emom
    integer, dimension(Nvertex), intent(inout) :: vertex_hit   !List to store visited vertices

    integer, dimension(:), allocatable :: nvertex_hit   ! List to store visited neighbour vertices

    real(dblprec) :: dflip(1), ve, ice_count,ice_count_up,ice_count_down ! ve is the vertex total moment

    integer :: i, l, k
    integer :: trial_vert, trial_nvert

    start_vertex=0;current_vertex=0;flip_atom=0;flip_atom_num_list=0;failure=.true.
    call ErrorHandling_missing('Spin-ice')

  end subroutine select_first_atom

  ! This subroutine selects the first vertex by checking if it follows the ice rules and if it has neighbours that follows the ice rules
  subroutine select_new_vertex(Natom, Nvertex, Mensemble, max_no_equiv_ver, max_no_equiv_ver_atom, nlist_ver, &
             nlist_ver_atom, nlistsize_ver, nlistsize_ver_atom, prev_vertex, current_vertex, next_vertex, &
             flip_atom, prev_flip_atom, prev_flip_atom_list_num, flip_atom_list_num,&
             atom_selected,failure,en,emom, vert_ice_coord)

    implicit none

    integer, intent(in) :: Natom
    integer, intent(in) :: Nvertex
    integer, intent(in) :: Mensemble
    integer, intent(in) :: en ! Current ensemble
    integer, intent(in) :: max_no_equiv_ver
    integer, intent(in) :: max_no_equiv_ver_atom
    integer, intent(in) :: prev_flip_atom ! previous atom marked to be flipped
    integer, intent(in) :: prev_flip_atom_list_num
    integer, dimension(max_no_equiv_ver,Nvertex), intent(in) :: nlist_ver ! Neighbour list for vertex vertex or vertex atoms
    integer, dimension(Nvertex), intent(in):: nlistsize_ver ! Size of neighbour list for vertex vertex or vertex atoms
    integer, dimension(max_no_equiv_ver_atom,Nvertex) , intent(in) :: nlist_ver_atom
    integer, dimension(Nvertex), intent(in) :: nlistsize_ver_atom
    real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: emom
    real(dblprec), dimension(3,Nvertex,max_no_equiv_ver_atom), intent(in) :: vert_ice_coord
    integer, intent(in) :: prev_vertex
    integer, intent(in) :: current_vertex
    integer, intent(inout) :: next_vertex
    integer, intent(out) :: flip_atom
    integer, intent(out) :: flip_atom_list_num

    logical, intent(inout) :: atom_selected
    logical, intent(out) :: failure

    integer, dimension(:), allocatable :: nvertex_hit   !List to store visited neighbour vertices

    real(dblprec) :: dflip(1), ve, ice_count,ice_count_up,ice_count_down
    real(dblprec) :: temp_flip_sign,temp_prev_flip_sign

    integer :: i, l, k

    flip_atom=0;flip_atom_list_num=0;failure=.true.
    call ErrorHandling_missing('Spin-ice')

  end subroutine select_new_vertex

  ! Set up vertex-vertex meighbour maps
  subroutine setup_ver_nm(Nvertex, NT_ver, NA_ver, N1, N2, N3, C1_ver, C2_ver, C3_ver, BC1, BC2, BC3, atype, Bas_ver, &
             max_no_equiv_ver, sym, ver_dist_coord, nm_ver_ver, nmdim_ver_ver)
    !
    !
    implicit none
    !
    ! Quatities from the repetitions of the cell that create the square lattice, they should be the same for the atom and vertex lattice
    !
    integer, intent(in) :: Nvertex ! The number of vertices in the system
    integer, intent(in) :: NT_ver ! Number of types of vertices
    integer, intent(in) :: NA_ver  ! Number of vertices in one cell
    integer, intent(in) :: N1  ! Number of cell repetitions in x direction
    integer, intent(in) :: N2  ! Number of cell repetitions in y direction
    integer, intent(in) :: N3  ! Number of cell repetitions in z direction
    real(dblprec), dimension(3), intent(in) :: C1_ver ! First lattice vector of vertex cell
    real(dblprec), dimension(3), intent(in) :: C2_ver ! Second lattice vector of vertex cell
    real(dblprec), dimension(3), intent(in) :: C3_ver ! Third lattice vector of vertex cell
    character(len=1), intent(in) :: BC1 !< Boundary conditions in x-direction
    character(len=1), intent(in) :: BC2 !< Boundary conditions in y-direction
    character(len=1), intent(in) :: BC3 !< Boundary conditions in z-direction
    integer, dimension(Nvertex), intent(in) :: atype !< Type of atom
    integer, intent(in) :: sym !< Symmetry of system (0-3)
    !
    ! Quantites that refer to the vertices only
    !
    real(dblprec), dimension(3,NA_ver), intent(in) :: Bas_ver ! Coordinates for basis vertex elements
    real(dblprec), dimension(NT_ver,3), intent(in) :: ver_dist_coord
    ! This are the ouput variables
    integer, intent(out) :: max_no_equiv_ver ! Calculated maximum of neighbours in one shell
    integer, dimension(:,:), allocatable, intent(out) :: nm_ver_ver ! Vertex-vertex neighbour map
    integer, dimension(:), allocatable, intent(out) :: nmdim_ver_ver ! Dimension of the vertex-vertex neighbour map
    !
    integer :: i0
    integer :: i_stat,i_all
    real(dblprec) :: tol
    real(dblprec) :: detmatrix
    real(dblprec), dimension(3,3) :: invmatrix
    !
    integer :: ix,iy,iz,xc_hop,yc_hop,zc_hop,ia,counter,itype,inei
    integer :: iat,jat,j0,jx,jy,jz
    real(dblprec), dimension(3) :: bsf,cvec,icvec,rvec
    real(dblprec), dimension(:,:,:), allocatable:: nncoord_ver !< Full list of neighbours for each type
    integer, dimension(:,:), allocatable:: nm_cell
    integer, dimension(:,:,:), allocatable :: nm_trunk
    integer, dimension(:), allocatable :: nnm_cell
    logical :: is_periodic,is_dilute
    !

    call ErrorHandling_missing('Spin-ice')

  end subroutine setup_ver_nm

  ! Set up vertex-vertex meighbour maps
  subroutine setup_ver_atom_nm(Nvertex, NT_ver, NA_ver, NA_atom, N1, N2, N3, C1_ver, C2_ver, C3_ver, BC1, BC2, BC3, atype, Bas_ver, Bas_atom, &
             max_no_equiv_ver_atom, sym, ver_atom_dist_coord, nm_ver_atom, nmdim_ver_atom)
    !
    !
    implicit none
    !
    ! Quatities from the repetitions of the cell that create the square lattice, they should be the same for the atom and vertex lattice
    !
    integer, intent(in) :: NA_atom  ! Number of atoms in one cell
    integer, intent(in) :: N1  !< Number of cell repetitions in x direction
    integer, intent(in) :: N2  !< Number of cell repetitions in y direction
    integer, intent(in) :: N3  !< Number of cell repetitions in z direction
    real(dblprec), dimension(3,NA_atom), intent(in) :: Bas_atom ! Coordinates of the atoms in the unit cell
    character(len=1), intent(in) :: BC1 !< Boundary conditions in x-direction
    character(len=1), intent(in) :: BC2 !< Boundary conditions in y-direction
    character(len=1), intent(in) :: BC3 !< Boundary conditions in z-direction
    integer, dimension(Nvertex), intent(in) :: atype !< Type of atom
    integer, intent(in) :: sym !< Symmetry of system (0-3)
    !
    ! Quantites that refer to the vertices only
    !
    integer, intent(in) :: Nvertex ! Number of vertices in the system
    integer, intent(in) :: NA_ver ! Number of vertices in the unit cell
    integer, intent(in) :: NT_ver ! Number of vertex types
    real(dblprec), dimension(3), intent(in) :: C1_ver ! First lattice vector of the vertex lattice
    real(dblprec), dimension(3), intent(in) :: C2_ver ! Second lattice vector of the vertex lattice
    real(dblprec), dimension(3), intent(in) :: C3_ver ! Third lattice vector of the vertex lattice
    real(dblprec), dimension(3,NA_ver), intent(in) :: Bas_ver ! Coordinates for basis vertex elements
    real(dblprec), dimension(NT_ver,3), intent(in) :: ver_atom_dist_coord ! Distance between neighbouring atoms and vertex
    ! This are the ouput variables
    integer, intent(out) :: max_no_equiv_ver_atom ! Calculated maximum of neighbours in one shell
    integer, dimension(:,:), allocatable, intent(out) :: nm_ver_atom ! Vertex-vertex neighbour map
    integer, dimension(:), allocatable, intent(out) :: nmdim_ver_atom ! Dimension of the vertex-vertex neighbour map
    !
    integer :: I0
    integer :: i_stat,i_all
    real(dblprec) :: tol
    real(dblprec) :: detmatrix
    real(dblprec), dimension(3,3) :: invmatrix
    !
    integer :: ix,iy,iz,xc_hop,yc_hop,zc_hop,ia,counter,itype,inei
    integer :: iat,jat,j0,jx,jy,jz
    real(dblprec), dimension(3) :: bsf,cvec,icvec,rvec
    real(dblprec), dimension(:,:,:), allocatable:: nncoord_ver_atom ! Full list of atoms neighbours for each vertex
    integer, dimension(:,:), allocatable:: nm_cell
    integer, dimension(:,:,:), allocatable :: nm_trunk
    integer, dimension(:), allocatable :: nnm_cell
    logical :: is_periodic,is_dilute
    !

    call ErrorHandling_missing('Spin-ice')

  end subroutine setup_ver_atom_nm

  !> Find possible symmetry operations depending on assumed symmetry
  subroutine get_symops(isym)
    !
    implicit none
    !
    integer, intent(in) :: isym !< Type of assumed symmetry (0-3)
    !
    integer :: i,j,x,y,z,j_s,x1,x2,y1,y2
    integer :: i_stat
    integer :: sym_count
    real(dblprec) :: half,roothalf
    !

    call ErrorHandling_missing('Spin-ice')

  end subroutine get_symops

  !> Create full neighbour vertex vertex list according to symmetry
  subroutine get_full_ver_ver_list(NT_ver, ver_dist_coord, max_no_equiv_ver, nncoord_ver)
    !
    !
    implicit none
    !
    integer, intent(in) :: NT_ver ! Number of types of vertices
    real(dblprec), dimension(NT_ver,3), intent(in) :: ver_dist_coord
    integer, intent(in) :: max_no_equiv_ver ! Calculated maximum of neighbours in one shell
    real(dblprec), dimension(3,max_no_equiv_ver,NT_ver), intent(out) :: nncoord_ver

    real(dblprec) :: tol
    real(dblprec), dimension(3) :: tvect
    integer :: counter
    logical :: unique
    integer :: i, j, k ,itype, isym

    nncoord_ver=0.0d0
    call ErrorHandling_missing('Spin-ice')

  end subroutine get_full_ver_ver_list

  ! Create full neighbour vertex atom list according to symmetry
  subroutine get_full_ver_atom_list(NT_ver, ver_atom_dist_coord, max_no_equiv_ver_atom, nncoord_ver_atom)
    !
    !
    implicit none
    !
    integer, intent(in) :: NT_ver ! Number of types of vertices
    real(dblprec), dimension(NT_ver,3), intent(in) :: ver_atom_dist_coord
    integer, intent(in) :: max_no_equiv_ver_atom ! Calculated maximum of neighbours in one shell for exchange
    real(dblprec), dimension(3,max_no_equiv_ver_atom,NT_ver), intent(out) :: nncoord_ver_atom


    real(dblprec) :: tol
    real(dblprec), dimension(3) :: tvect
    integer :: counter
    logical :: unique
    integer :: i, j, k ,itype, isym

    nncoord_ver_atom=0.0d0
    call ErrorHandling_missing('Spin-ice')

  end subroutine get_full_ver_atom_list

  ! Mounts the actual Neighbour list
  subroutine setup_neighbour_list(Nvertex, max_no_equiv, nlistsize, nlist, nm, nmdim)
    !
    implicit none
    !
    integer, intent(in) :: Nvertex ! Number of vertices in system
    integer, dimension(Nvertex), intent(out):: nlistsize ! Size of neighbour list for vertex vertex or vertex atoms
    integer, intent(in) :: max_no_equiv ! Calculated maximum of neighbours in one shell for vertex vertex of vertex atoms
    integer, dimension(max_no_equiv,Nvertex), intent(out) :: nlist ! Neighbour list for vertex vertex or vertex atoms
    integer, dimension(Nvertex,max_no_equiv), intent(in) :: nm ! Vertex-Vertex or Vertex-atom neighbour map
    integer, dimension(Nvertex), intent(in) :: nmdim ! Dimension of the vertex-vertex (vertex-atom) neighbour map
    !
    integer :: i, j, l
    integer :: iver, jver
    integer :: ncount
    logical :: exis

    nlistsize=0;nlist=0
    call ErrorHandling_missing('Spin-ice')

  end subroutine setup_neighbour_list

! logical function ice_rule(Natom,Mensemble, k,emom,ii,ve,ice_count,ice_count_up,ice_count_down)
  logical function ice_rule(Natom,Mensemble, en,emom,trial_vert,total_trial,ice_count,ice_count_up, ice_count_down)

    implicit none

    integer, intent(in) :: Natom
    integer, intent(in) :: Mensemble
    integer, intent(in) :: en ! Current ensemble
    real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: emom
    integer, intent(in) :: trial_vert
    real(dblprec), intent(out) :: total_trial

    real(dblprec), intent(out) :: ice_count_up, ice_count_down
    real(dblprec), intent(out) :: ice_count

    integer :: max_ver_size
    real(dblprec) :: trial_ene
    real(dblprec) :: temp_trial
    real(dblprec), dimension(3) :: trial_ene_vec

    integer :: atom

    ice_rule=.false.;total_trial=0.0d0;ice_count_up=0.0d0;ice_count_down=0.0d0;ice_count=0.0d0
    call ErrorHandling_missing('Spin-ice')

  end function ice_rule

!> Calculate total energy for a loop in the spin Ice
subroutine calculate_spinice_energy(Natom, Mensemble, max_no_neigh, conf_num, ncoup, nlist, nlistsize, &
           do_dm, max_no_dmneigh, dm_vect, dmlist, dmlistsize, &
           do_pd, nn_pd_tot, pd_vect, pdlist, pdlistsize, &
           do_biqdm, nn_biqdm_tot, biqdm_vect, biqdmlist, biqdmlistsize, &
           do_bq, nn_bq_tot, j_bq, bqlist, bqlistsize, &
           taniso, eaniso, kaniso,sb,emomM, emom, &
           Ncount_atom, loop_len, extfield, de, k, &
           do_dip, Qdip)

  use Constants, only : mub

  !.. Implicit declarations
  implicit none

  integer, intent(in) :: Natom !< Number of atoms in system
  integer, intent(in) :: Mensemble !< Number of ensembles
  integer, intent(in) :: max_no_neigh !< Calculated maximum of neighbours for exchange
  integer, intent(in) :: conf_num   !< Number of configurations for LSF
  real(dblprec), dimension(max_no_neigh,Natom,conf_num), intent(in) :: ncoup !< Heisenberg exchange couplings
  integer, dimension(max_no_neigh,Natom), intent(in) :: nlist !< Neighbour list for Heisenberg exchange couplings
  integer, dimension(Natom),intent(in) :: nlistsize !< Size of neighbour list for Heisenberg exchange couplings
  integer, intent(in) :: do_dm   !< Add Dzyaloshinskii-Moriya (DM) term to Hamiltonian (0/1)
  integer, intent(in) :: max_no_dmneigh !< Calculated number of neighbours with DM interactions
  real(dblprec), dimension(3,max_no_dmneigh,Natom), intent(in) :: dm_vect !< Dzyaloshinskii-Moriya exchange vector
  integer, dimension(max_no_dmneigh,Natom), intent(in) :: dmlist   !< List of neighbours for DM
  integer, dimension(Natom),intent(in) :: dmlistsize !< Size of neighbour list for DM
  integer, intent(in) :: do_pd   !< Add Pseudo-Dipolar (PD) term to Hamiltonian (0/1)
  integer, intent(in) :: nn_pd_tot !< Calculated number of neighbours with PD interactions
  real(dblprec), dimension(6,nn_pd_tot,Natom), intent(in) :: pd_vect !< Pseudo-Dipolar exchange vector
  integer, dimension(nn_pd_tot,Natom), intent(in) :: pdlist   !< List of neighbours for PD
  integer, dimension(Natom),intent(in) :: pdlistsize !< Size of neighbour list for PD
  integer, intent(in) :: do_biqdm   !< Add biquadratic DM (BIQDM) term to Hamiltonian (0/1)
  integer, intent(in) :: nn_biqdm_tot !< Calculated number of neighbours with BIQDM interactions
  real(dblprec), dimension(1,nn_biqdm_tot,Natom), intent(in) :: biqdm_vect !< BIQDM exchange coupling
  integer, dimension(nn_biqdm_tot,Natom), intent(in) :: biqdmlist   !< List of neighbours for BIQDM
  integer, dimension(Natom),intent(in) :: biqdmlistsize !< Size of neighbour list for BIQDM
  integer, intent(in) :: do_bq   !< Add biquadratic exchange (BQ) term to Hamiltonian (0/1)
  integer, intent(in) :: nn_bq_tot !< Calculated number of neighbours with BQ interactions
  real(dblprec), dimension(nn_bq_tot,Natom), intent(in) :: j_bq !< Biquadratic exchange couplings
  integer, dimension(nn_bq_tot,Natom), intent(in) :: bqlist   !< List of neighbours for BQ
  integer, dimension(Natom),intent(in) :: bqlistsize !< Size of neighbour list for BQ
  integer, dimension(Natom),intent(in) :: taniso !< Type of anisotropy (0-2)
  real(dblprec), dimension(3,Natom), intent(in) :: eaniso !< Unit anisotropy vector
  real(dblprec), dimension(2,Natom), intent(in) :: kaniso !< Anisotropy constant
  real(dblprec), dimension(Natom), intent(in) :: sb!< Ratio between the Anisotropy constants
  real(dblprec), dimension(3,Natom,Mensemble), intent(inout) :: emomM  !< Current magnetic moment vector
  real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: emom   !< Current unit moment vector
  integer, intent(in) :: loop_len !< Lenght of cluster loop
  integer, dimension(Natom), intent(in) :: Ncount_atom !< List of atoms in cluster loop
  real(dblprec), dimension(3), intent(in) :: extfield !< External magnetic field
  real(dblprec), intent(out):: de  !< Energy difference
  integer, intent(in) :: k !< Current ensemble
  real(dblprec) :: aw1,aw2
  integer, intent(in) :: do_dip  !<  Calculate dipole-dipole contribution (0/1)
  real(dblprec), dimension(3,3,Natom,Natom), intent(in), optional :: Qdip !< Matrix for dipole-dipole interaction

  !.. Local scalars
  integer :: j,iloop,cloop, catom
  real(dblprec) :: tta,ttb,e_c,e_t
  real(dblprec) :: bqmdot

  !.. Local arrays

  real(dblprec), dimension(2) :: tt_loop
  real(dblprec), dimension(:,:), allocatable :: befftemp_loop

  !.. Executable statements

    de=0.0d0
    call ErrorHandling_missing('Spin-ice')

end subroutine calculate_spinice_energy


!> Flip spin in spin-ice
subroutine Spin_Ice_flip(emom, emomM, mmom, Ncount_atom, Natom, Mensemble,Temperature, de, k,loop_len)

  use Constants
  use RandomNumbers, only : rng_uniform

  !.. Implicit declarations
  implicit none

  integer, intent(in) :: Natom !< Number of atoms in system
  integer, intent(in) :: loop_len
  integer, intent(in) :: Mensemble !< Number of ensembles
  integer, dimension(Natom), intent(in) :: Ncount_atom
  real(dblprec), dimension(3,Natom,Mensemble) :: emom   !< Current unit moment vector
  real(dblprec), dimension(3,Natom,Mensemble) :: emomM  !< Current magnetic moment vector
  real(dblprec), dimension(Natom,Mensemble), intent(in) :: mmom !< Magnitude of magnetic moments
  real(dblprec), intent(in):: de  !< Energy difference
  real(dblprec), intent(in):: Temperature !< Temperature

  integer, intent(in) :: k !< Current ensemble
  integer :: i
  real(dblprec) :: beta , flipprob(1)

    call ErrorHandling_missing('Spin-ice')


end subroutine Spin_Ice_flip

end module SpinIce
