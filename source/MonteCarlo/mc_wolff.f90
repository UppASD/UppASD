!> Wolff cluster Monte Carlo
!> @copyright
!! Copyright (C) 2008-2018 UppASD group
!! This file is distributed under the terms of the
!! GNU General Public License. 
!! See http://www.gnu.org/copyleft/gpl.txt
module MC_WOLFF

  use Parameters
  use Constants
  use Profiling
   use ErrorHandling

  implicit none

  integer, dimension(:), allocatable :: wolff_list_neigh_size
  integer, dimension(:,:), allocatable :: wolff_list_neigh

  public :: wolff_list_neigh_size, wolff_list_neigh

contains

  ! This subroutine creates the neighbour list for the wolff algorithm in which the nearest neighbours are needed
  subroutine create_wolff_neigh_list(Natom, N1, N2, N3, C1, C2, C3, BC1, BC2, BC3, &
             max_no_neigh, nlistsize, nlist, coord,wolff_list_neigh_size,wolff_list_neigh)

    implicit none

    integer, intent(in) :: Natom !< Number of atoms in system
    integer, intent(in) :: N1    !< Number of cell repetitions in x direction
    integer, intent(in) :: N2    !< Number of cell repetitions in y direction
    integer, intent(in) :: N3    !< Number of cell repetitions in z direction
    real(dblprec), dimension(3), intent(in) :: C1 !< First lattice vector
    real(dblprec), dimension(3), intent(in) :: C2 !< Second lattice vector
    real(dblprec), dimension(3), intent(in) :: C3 !< Third lattice vector
    character(len=1), intent(in) :: BC1 !< Boundary conditions in x-direction
    character(len=1), intent(in) :: BC2 !< Boundary conditions in y-direction
    character(len=1), intent(in) :: BC3 !< Boundary conditions in z-direction
    integer, intent(in) :: max_no_neigh !< Calculated maximum of neighbours for exchange
    integer, dimension(Natom), intent(in) :: nlistsize !< Size of neighbour list for Heisenberg exchange couplings
    integer, dimension(max_no_neigh, Natom), intent(in) :: nlist !< Neighbour list for Heisenberg exchange couplings
    real(dblprec), dimension(3,Natom), intent(in) :: coord !< Coordinates of atoms

    integer, dimension(:), allocatable :: wolff_list_neigh_size
    integer, dimension(:,:), allocatable :: wolff_list_neigh

    real(dblprec), parameter :: rcutoff=1.05d0

    integer :: iatom, jneigh, jatom, ix, iy, iz, icount
    integer :: i_stat
    integer :: signx, signy, signz
    integer :: perx, pery, perz
    real(dblprec) :: jnorm, jpnorm, jtnorm
    real(dblprec), dimension(3) :: icoord, jcoord, jpcoord, jtcoord
    real(dblprec), dimension(3,-1:1,-1:1,-1:1) :: minmat_coord
    real(dblprec), dimension(-1:1,-1:1,-1:1) :: minmat_norm
    integer, dimension(-1:1,-1:1,-1:1) :: minmat_atom

            call ErrorHandling_missing('Wolff algorithm')

  end subroutine create_wolff_neigh_list

  ! This createds the cluster for the Wolff algorithm
  subroutine create_wolff_cluster(k_ensem,Natom,Mensemble,emom,nlistsize,nlist,&
             Temperature,conf_num,ncoup,max_no_neigh,visited_atoms)

    use RandomNumbers
    implicit none

    integer, intent(in) :: k_ensem
    integer, intent(in) :: Natom
    integer, intent(in) :: Mensemble
    integer, intent(in) :: max_no_neigh

    integer, dimension(Natom), intent(inout) :: visited_atoms
    integer, dimension(:), allocatable :: wolffsnake
    integer, dimension(:), allocatable :: border_atoms
    real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: emom

    integer, dimension(Natom), intent(in) :: nlistsize
    integer, dimension(max_no_neigh,Natom), intent(in) :: nlist
    integer, intent(in) :: conf_num !> number of configurations for LSF
    real(dblprec), intent(in) :: Temperature

    real(dblprec), dimension(max_no_neigh,Natom,conf_num), intent(in) :: ncoup

    real(dblprec) :: dflip(1)
    integer :: border_index, i_border, snake_head, snake_tail
    integer :: atom_wolff_len, selected_atom,curr_len
    logical :: loop_closed, failure, total_failure,new_border

            call ErrorHandling_missing('Wolff algorithm')

  end subroutine create_wolff_cluster

  ! Check the neighbours to see if they can be added into the cluster
  subroutine neighbour_check(k_ensem,Natom,Mensemble,selected_atom,emom,nlist,nlistsize,max_no_neigh,&
             border_atoms,atom_wolff_len,wolffsnake,border_index,curr_len,new_border,conf_num,ncoup,&
             Temperature,visited_atoms)

    use RandomNumbers

    implicit none

    integer, intent(in) :: k_ensem
    integer, intent(in) :: Natom
    integer, intent(in) :: Mensemble
    integer, intent(in) :: selected_atom
    integer, intent(in) :: max_no_neigh
    integer, intent(inout) :: atom_wolff_len
    integer, intent(inout) :: border_index
    integer, dimension(Natom), intent(inout) :: visited_atoms
    integer, dimension(Natom), intent(inout) :: border_atoms
    integer, dimension(Natom), intent(inout) :: wolffsnake
    real(dblprec), intent(in) :: Temperature
    integer, intent(in) :: conf_num   !< Number of configurations for LSF
    real(dblprec), dimension(max_no_neigh,Natom,conf_num), intent(in) :: ncoup
    real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: emom
    integer, dimension(Natom), intent(in) :: nlistsize
    integer, dimension(max_no_neigh,Natom), intent(in) :: nlist

    integer, intent(inout) :: curr_len
    logical, intent(out) :: new_border

    integer :: neigh_loop
    real(dblprec) :: flipprob(1)

            new_border=.false.
            call ErrorHandling_missing('Wolff algorithm')

  end subroutine neighbour_check

  ! This function decides if a link between two atoms is added with a certain probability
  logical function prob_add(Natom,max_no_neigh,i_atom,j_atom,Temperature,conf_num,ncoup)

     use RandomNumbers
     use Constants

     implicit none

     integer, intent(in) :: Natom
     integer, intent(in) :: max_no_neigh
     integer, intent(in) :: i_atom
     integer, intent(in) :: j_atom
     real(dblprec), intent(in) :: Temperature
     integer, intent(in) :: conf_num   !< Number of configurations for LSF
     real(dblprec), dimension(max_no_neigh,Natom,conf_num), intent(in) :: ncoup
     real(dblprec) :: beta
     real(dblprec) :: flipprob(1), prob

     prob_add=.false.
            call ErrorHandling_missing('Wolff algorithm')

  end function prob_add

end module MC_WOLFF
