!> @copyright
!! Copyright (C) 2008-2018 UppASD group
!! This file is distributed under the terms of the
!! GNU General Public License. 
!! See http://www.gnu.org/copyleft/gpl.txt
module SpinIceData

   use Parameters
   use Profiling
   !
   implicit none
   !
   !From setup
   integer, dimension(:), allocatable :: anumb_ver              !< Vertex number in cell
   integer, dimension(:), allocatable :: anumb_inp_ver          !< Vertex number in cell from input
   integer, dimension(:), allocatable :: atype_ver              !< Type of vertices
   integer, dimension(:), allocatable :: atype_inp_ver          !< Vertex type in cell from input
   real(dblprec), dimension(:,:), allocatable :: coord_vertices !< Coordinates of vertices 
   real(dblprec), dimension(:,:,:), allocatable :: vert_ice_coord

   integer ::  max_no_equiv_ver       !< Calculated maximum of neighbours for vertex-vertex geometry
   integer ::  max_no_equiv_ver_atom  !< Calculated maximum of neighbours for vertex-atom geometry
   integer, dimension(:), allocatable :: nlistsize_ver      !< Size of neighbour list for vertex-vertex
   integer, dimension(:), allocatable :: nlistsize_ver_atom !< Size of neighbour list for vertex-atom
   integer, dimension(:,:), allocatable :: nlist_ver        !< Neighbour list for vertex-vertex
   integer, dimension(:,:), allocatable :: nlist_ver_atom   !< List of neighbours for vertex-atom

   ! Geometry of the vertices in the square lattice spin ice system
   integer :: NA_ver                                                 !< Number of vertices in the unit cell
   integer :: NT_ver                                                 !< Number of types of vertices
   integer :: Nvertex                                                !< Total number of vertices
   integer :: ver_sym                                                !< Symmetry of the vertex system
   character(len=35) :: vertex                                       !< Vertex file name
   real(dblprec), dimension(3) :: C1_ver                             !< First vertex lattice vector
   real(dblprec), dimension(3) :: C2_ver                             !< Second vertex lattice vector
   real(dblprec), dimension(3) :: C3_ver                             !< Third vertex lattice vector
   real(dblprec), dimension(:,:), allocatable :: Bas_ver             !< Coordinates for basis vertices
   real(dblprec), dimension(:,:), allocatable :: ver_dist_coord      !< Coordinates for vertex-vertex distance
   real(dblprec), dimension(:,:), allocatable :: ver_atom_dist_coord !< Coordinates for vertex-atom distance

   ! Vertex printing variables  
   integer :: mchits                !< Total number of spin flips
   integer :: ver_no                !< Number of vertices which follow ICE RULE
   integer :: vertex_step           !< Interval for sampling vertices
   integer :: vertex_buff           !< Buffer size for vertex printing (the user should really set it equal to tottraj buff)
   integer :: bcount_vertex         !< Counter of buffer for vertices
   integer :: mchits_spin_ice       !< Total number of cluster flips
   character(len=1) :: prn_vertices !< Flag for vertex printing

   real(dblprec) :: loop_ave_len !< Average length of a loop
   real(dblprec), dimension(:), allocatable :: indxb_vertex          !< Step counter for vertex printing
   real(dblprec), dimension(:,:,:), allocatable :: ice_count_buffer    
   real(dblprec), dimension(:,:,:), allocatable :: total_trial_buffer  !< Buffer for the total "energy" in the vertices  
   logical, dimension(:,:,:), allocatable :: ice_rule_buffer

   public


contains

   !> Allocate arrays for loop algorithm, working and orinting
   subroutine allocate_spinice_data(Natom,Mensemble,flag)

      implicit none

      integer, intent(in) :: flag  !< Allocate or deallocate (1/-1)
      integer, intent(in) :: Natom
      integer, intent(in) :: Mensemble

      integer :: i_all, i_stat

      ! Allocation of arrays and initialization
      if(flag>0) then

         allocate(nlistsize_ver(Nvertex),stat=i_stat)
         call memocc(i_stat,product(shape(nlistsize_ver))*kind(nlistsize_ver),'nlistsize_ver','allocate_spinicedata')

         allocate(nlistsize_ver_atom(Nvertex),stat=i_stat)
         call memocc(i_stat,product(shape(nlistsize_ver_atom))*kind(nlistsize_ver_atom),'nlistsize_ver_atom','allocate_spinicedata')

         allocate(nlist_ver(max_no_equiv_ver,Nvertex),stat=i_stat)
         call memocc(i_stat,product(shape(nlist_ver))*kind(nlist_ver),'nlist_ver','allocate_spinicedata')

         allocate(nlist_ver_atom(max_no_equiv_ver_atom,Nvertex),stat=i_stat)
         call memocc(i_stat,product(shape(nlist_ver_atom))*kind(nlist_ver_atom),'nlist_ver_atom','allocate_spinice_data')

         allocate(vert_ice_coord(3,Nvertex,max_no_equiv_ver_atom),stat=i_stat)
         call memocc(i_stat,product(shape(vert_ice_coord))*kind(vert_ice_coord),'vert_ice_coord','allocate_spinice_data')

         bcount_vertex=1

         if (prn_vertices=='Y') then
            allocate(indxb_vertex(vertex_buff),stat=i_stat)
            call memocc(i_stat,product(shape(indxb_vertex))*kind(indxb_vertex),'indxb_vertex','allocate_spinicedata')
            allocate(ice_count_buffer(Nvertex,Mensemble,vertex_buff),stat=i_stat)
            call memocc(i_stat,product(shape(ice_count_buffer))*kind(ice_count_buffer),'ice_count_buffer','allocate_spinicedata')
            allocate(total_trial_buffer(Nvertex,Mensemble,vertex_buff),stat=i_stat)
            call memocc(i_stat,product(shape(total_trial_buffer))*kind(total_trial_buffer),'total_trial_buffer','allocate_spinicedata')
            allocate(ice_rule_buffer(Nvertex,Mensemble,vertex_buff),stat=i_stat)
            call memocc(i_stat,product(shape(ice_rule_buffer))*kind(ice_rule_buffer),'ice_rule_buffer','allocate_spinicedata')
         endif

      else
         i_all=-product(shape(nlistsize_ver))*kind(nlistsize_ver)
         deallocate(nlistsize_ver,stat=i_stat)
         call memocc(i_stat,i_all,'nlistsize_ver','allocate_spinicedata')

         i_all=-product(shape(nlistsize_ver_atom))*kind(nlistsize_ver_atom)
         deallocate(nlistsize_ver_atom,stat=i_stat)
         call memocc(i_stat,i_all,'nlistsize_ver_atom','allocate_spinicedata')

         i_all=-product(shape(nlist_ver))*kind(nlist_ver)
         deallocate(nlist_ver,stat=i_stat)
         call memocc(i_stat,i_all,'nlist_ver','allocate_spinicedata')

         i_all=-product(shape(nlist_ver_atom))*kind(nlist_ver_atom)
         deallocate(nlist_ver_atom,stat=i_stat)
         call memocc(i_stat,i_all,'nlist_ver_atom','allocate_spinicedata')

         i_all=-product(shape(vert_ice_coord))*kind(vert_ice_coord)
         deallocate(vert_ice_coord,stat=i_stat)
         call memocc(i_stat,i_all,'vert_ice_coord','allocate_spinicedata')


         if (prn_vertices=='Y') then
            i_all=-product(shape(indxb_vertex))*kind(indxb_vertex)
            deallocate(indxb_vertex,stat=i_stat)
            call memocc(i_stat,i_all,'indxb_vertex','allocate_spinicedata')

            i_all=-product(shape(ice_rule_buffer))*kind(ice_rule_buffer)
            deallocate(ice_rule_buffer,stat=i_stat)
            call memocc(i_stat,i_all,'ice_rule_buffer','allocate_spinicedata')

            i_all=-product(shape(ice_count_buffer))*kind(ice_count_buffer)
            deallocate(ice_count_buffer,stat=i_stat)
            call memocc(i_stat,i_all,'ice_count_buffer','allocate_spinicedata')

            i_all=-product(shape(total_trial_buffer))*kind(total_trial_buffer)
            deallocate(total_trial_buffer,stat=i_stat)
            call memocc(i_stat,i_all,'total_trial_buffer','allocate_spinicedata')

         endif

      end if

   end subroutine allocate_spinice_data

   !> Allocate necessary data for vertex creation
   subroutine allocate_vertex_data(flag)

      integer, intent(in) :: flag

      integer :: i_all, i_stat

      if (flag>0) then

         allocate(anumb_ver(Nvertex),stat=i_stat)
         call memocc(i_stat,product(shape(anumb_ver))*kind(anumb_ver),'anumb_ver','allocate_spinicedata')

         allocate(atype_ver(Nvertex),stat=i_stat)
         call memocc(i_stat,product(shape(atype_ver))*kind(atype_ver),'atype_ver','allocate_spinicedata')

      else

         i_all=-product(shape(anumb_ver))*kind(anumb_ver)
         deallocate(anumb_ver,stat=i_stat)
         call memocc(i_stat,i_all,'anumb_ver','allocate_spinicedata')

         i_all=-product(shape(atype_ver))*kind(atype_ver)
         deallocate(atype_ver,stat=i_stat)
         call memocc(i_stat,i_all,'atype_ver','allocate_spinicedata')

      end if

   end subroutine allocate_vertex_data

end module SpinIceData
