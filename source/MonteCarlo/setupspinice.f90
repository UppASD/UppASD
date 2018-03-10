!> Spin-ice simulation routines
!> @copyright
!! Copyright (C) 2008-2018 UppASD group
!! This file is distributed under the terms of the
!! GNU General Public License. 
!! See http://www.gnu.org/copyleft/gpl.txt
module SetupSpinIce

   use Parameters
   use Profiling
   use Constants
   !use Geometry
   use SpinIce
   use SpinIceData
   use prn_SpinIce
   use ErrorHandling

   implicit none

   integer, dimension(:,:), allocatable :: nm_ver_ver
   integer, dimension(:,:), allocatable :: nm_ver_atom
   integer, dimension(:), allocatable :: nmdim_ver_ver ! Dimension of neighbour map
   integer, dimension(:), allocatable :: nmdim_ver_atom

   private

   public :: setup_ice_neighbours,setup_vertex_geometry
   public :: read_parameters_spinice, spin_ice_init, read_vertices

contains

   !> setup system
   subroutine setup_ice_neighbours(Natom, Mensemble,NT,NA, N1, N2, N3, BC1, BC2, BC3, atype, Bas, sym, simid,coord)

      ! use Geometry, only : coord_vertices

      !.. Implicit declarations
      implicit none

      integer, intent(inout) :: Natom ! Number of atoms in system
      integer, intent(in) :: NT ! Number of types of atoms
      integer, intent(in) :: NA  ! Number of atoms in one cell
      integer, intent(in) :: N1  ! Number of cell repetitions in x direction
      integer, intent(in) :: N2  ! Number of cell repetitions in y direction
      integer, intent(in) :: N3  ! Number of cell repetitions in z direction
      integer, intent(in) :: Mensemble
      character(len=1), intent(in) :: BC1 ! Boundary conditions in x-direction
      character(len=1), intent(in) :: BC2 ! Boundary conditions in y-direction
      character(len=1), intent(in) :: BC3 ! Boundary conditions in z-direction
      integer, dimension(Natom), intent(inout) :: atype
      real(dblprec), dimension(3,NA), intent(inout) :: Bas !< Coordinates for basis atoms
      integer, intent(in) :: sym ! Symmetry of system (0-3)
      real(dblprec), dimension(3,Natom), intent(in) :: coord
      character(len=8) :: simid ! Name of simulation


      character(len=1) :: neigh_type

      call ErrorHandling_missing('Spin-ice')


   end subroutine setup_ice_neighbours

   !> deallocation of arrays
   subroutine deallocate_nm_spin_ice()

      implicit none

      integer :: i_stat, i_all

      call ErrorHandling_missing('Spin-ice')

   end subroutine deallocate_nm_spin_ice

   !> Initialization of the variables needed for the Spin-Ice measurements
   subroutine spin_ice_init()

      implicit none

      real(dblprec) :: one=1.0_dblprec
      real(dblprec) :: zero=0.0_dblprec

      !Vertices
      NA_ver  = 0
      NT_ver  = 0
      ver_sym = 0
      C1_ver  = (/one,zero,zero/)
      C2_ver  = (/zero,one,zero/)
      C3_ver  = (/zero,zero,one/)
      vertex  = 'vertexfile'

      ver_no          = 0
      mchits          = 0
      vertex_step     = 1000
      vertex_buff     = 10
      prn_vertices    = 'N'
      loop_ave_len    = 0.0d0
      mchits_spin_ice = 0

   end subroutine spin_ice_init

   !> Setup vertices to first supercell.
   subroutine setup_globvertices(Nvertex, NA_ver, Bas_ver, C1_ver, C2_ver, C3_ver, N1, N2, N3, atype_ver, anumb_ver, do_prnstruct, simid,coord_vertices)
      !
      !
      implicit none
      !
      integer, intent(inout) :: Nvertex !< Number of atoms in system
      integer, intent(in) :: NA_ver  !< Number of atoms in one cell
      real(dblprec), dimension(3,NA_ver) , intent(inout) :: Bas_ver !< Coordinates for basis atoms
      real(dblprec), dimension(3) , intent(in) :: C1_ver !< First lattice vector
      real(dblprec), dimension(3) , intent(in) :: C2_ver !< Second lattice vector
      real(dblprec), dimension(3) , intent(in) :: C3_ver !< Third lattice vector
      real(dblprec), dimension(:,:),allocatable,intent(inout) :: coord_vertices !< Coordinates of atoms
      integer, intent(in) :: N1  !< Number of cell repetitions in x direction
      integer, intent(in) :: N2  !< Number of cell repetitions in y direction
      integer, intent(in) :: N3  !< Number of cell repetitions in z direction
      integer, dimension(Nvertex), intent(in) :: atype_ver !< Type of atom
      integer, dimension(Nvertex), intent(in) :: anumb_ver !< Atom number in cell
      integer, intent(in) :: do_prnstruct !< Print Hamiltonian information (0/1)
      character(len=8) :: simid !< Name of simulation

      integer :: i,i_stat
      integer :: I0, i1, i2, i3
      integer :: iatom
      character(len=20) :: filn
      real(dblprec) :: detmatrix
      real(dblprec), dimension(3) :: icvec, bsf
      real(dblprec), dimension(3,3) :: invmatrix
      !

      call ErrorHandling_missing('Spin-ice')

   end subroutine setup_globvertices

   !> Sets up the type of vertices in the system
   subroutine setup_vertex_type_and_numb(Nvertex, NA_ver, N1, N2, N3, atype_ver, anumb_ver, atype_inp_ver, anumb_inp_ver)
      !
      !
      implicit none

      integer, intent(in) :: Nvertex ! Number of atoms in system
      integer, intent(in) :: NA_ver  ! Number of atoms in one cell
      integer, intent(in) :: N1  ! Number of cell repetitions in x direction
      integer, intent(in) :: N2  ! Number of cell repetitions in y direction
      integer, intent(in) :: N3  ! Number of cell repetitions in z direction
      integer, dimension(Nvertex), intent(inout) :: atype_ver ! Type of atom
      integer, dimension(Nvertex), intent(inout) :: anumb_ver ! Atom number in cell
      integer, dimension(NA_ver), intent(inout) :: atype_inp_ver  ! Type of atom from input
      integer, dimension(NA_ver), intent(inout) :: anumb_inp_ver ! Atom number in cell from input
      !
      integer :: A1
      integer :: i0, i1, i2, i3
      !
      !
      call ErrorHandling_missing('Spin-ice')

   end subroutine setup_vertex_type_and_numb

   !> Wrapper routine for setting up structural and chemical information about the system
   subroutine setup_vertex_geometry( N1, N2, N3, do_prnstruct, simid)


      !
      implicit none

      integer, intent(in) :: N1  ! Number of cell repetitions in x direction
      integer, intent(in) :: N2  ! Number of cell repetitions in y direction
      integer, intent(in) :: N3  ! Number of cell repetitions in z direction
      integer, intent(in) :: do_prnstruct !< Print Hamiltonian information (0/1)
      character(len=8), intent(in) :: simid !< Name of simulation
      !
      integer :: i_stat
      !
      call ErrorHandling_missing('Spin-ice')

   end subroutine setup_vertex_geometry


   subroutine GEOMETRICAL_ICE(Natom,Nvertex,nlistsize_ver_atom,nlist_ver_atom,max_no_equiv_ver_atom,coord,coord_vertices,vert_ice_coord,simid, &
         BC1,BC2,BC3,ver_atom_dist_coord,NT_ver,atype_ver)

      implicit none

      integer, intent(in) :: NT_ver
      integer, intent(in) :: Natom
      integer, intent(in) :: Nvertex
      integer, intent(in) :: max_no_equiv_ver_atom
      integer, dimension(Nvertex), intent(in) :: atype_ver ! Type of vertex
      integer, dimension(Nvertex), intent(in) :: nlistsize_ver_atom
      integer, dimension(max_no_equiv_ver_atom,Nvertex) , intent(in) :: nlist_ver_atom

      real(dblprec), dimension(NT_ver,3), intent(in) :: ver_atom_dist_coord
      real(dblprec), dimension(3,Natom), intent(in) :: coord
      real(dblprec), dimension(3,Nvertex), intent(in) :: coord_vertices

      character(len=8), intent(in) :: simid
      character(len=1), intent(in) :: BC1 ! Boundary conditions in x-direction
      character(len=1), intent(in) :: BC2 ! Boundary conditions in y-direction
      character(len=1), intent(in) :: BC3 ! Boundary conditions in z-direction

      real(dblprec), dimension(3,Nvertex,max_no_equiv_ver_atom), intent(out) :: vert_ice_coord
      real(dblprec) :: norm,fac, new_norm
      real(dblprec), dimension(3) :: temp

      integer :: vertex, atom

      character(len=30) :: filn

      vert_ice_coord=0.0d0
      call ErrorHandling_missing('Spin-ice')

   end subroutine GEOMETRICAL_ICE


   !---------------------------------------------------------------------------
   !> @brief
   !> Read input parameters.
   !
   !> @author
   !> Anders Bergman
   !---------------------------------------------------------------------------
   subroutine read_parameters_spinice(ifile)
      use FileParser

      implicit none

      ! ... Formal Arguments ...
      integer, intent(in) :: ifile   !< File to read from
      !
      ! ... Local Variables ...
      character(len=50) :: keyword, cache
      integer :: rd_len, i_err, i_errb
      logical :: comment


      call ErrorHandling_missing('Spin-ice')

   end subroutine read_parameters_spinice

   subroutine read_vertices()
      use FileParser

      implicit none

      integer :: i_err,rd_len,i_errb, iat,i_stat
      character(len=50) :: keyword
      logical :: comment

      call ErrorHandling_missing('Spin-ice')

   end subroutine read_vertices


end module SetupSpinIce
