!-------------------------------------------------------------------------------
! MODULE: HamiltonianMacroBlocks
!> @brief Build coarse macroblock metadata for Hamiltonian execution backends.
!> @details Reuses the existing macrocell decomposition and groups finalized
!>          Hamiltonian neighbour entries by interacting macroblock pair.
!-------------------------------------------------------------------------------
module HamiltonianMacroBlocks

   use iso_c_binding, only : c_int
   use Parameters
   use macrocells, only : do_macro_cells, Num_macro, cell_index, macro_nlistsize,  &
      macro_atom_nlist

   implicit none

   private

   type, public :: macroblock_layout_type
      logical :: ready = .false.
      integer(c_int) :: nblocks = 0
      integer(c_int) :: n_pair_groups = 0
      integer(c_int) :: n_entries = 0
      integer(c_int), allocatable :: atom_to_block(:)
      integer(c_int), allocatable :: block_atom_count(:)
      integer(c_int), allocatable :: block_atom_offset(:)
      integer(c_int), allocatable :: block_atoms(:)
      integer(c_int), allocatable :: block_neigh_count(:)
      integer(c_int), allocatable :: block_neigh_offset(:)
      integer(c_int), allocatable :: block_neigh_list(:)
      integer(c_int), allocatable :: pair_group_src(:)
      integer(c_int), allocatable :: pair_group_dst(:)
      integer(c_int), allocatable :: pair_group_entry_offset(:)
      integer(c_int), allocatable :: entry_atom_i(:)
      integer(c_int), allocatable :: entry_atom_j(:)
      integer(c_int), allocatable :: entry_ih(:)
      integer(c_int), allocatable :: entry_jslot(:)
   end type macroblock_layout_type

   type(macroblock_layout_type), public, save :: ham_macroblock_layout

   public :: build_macroblock_layout
   public :: destroy_macroblock_layout
   public :: macroblock_layout_is_ready

contains

   subroutine destroy_macroblock_layout(layout)
      type(macroblock_layout_type), intent(inout) :: layout

      if (allocated(layout%atom_to_block)) deallocate(layout%atom_to_block)
      if (allocated(layout%block_atom_count)) deallocate(layout%block_atom_count)
      if (allocated(layout%block_atom_offset)) deallocate(layout%block_atom_offset)
      if (allocated(layout%block_atoms)) deallocate(layout%block_atoms)
      if (allocated(layout%block_neigh_count)) deallocate(layout%block_neigh_count)
      if (allocated(layout%block_neigh_offset)) deallocate(layout%block_neigh_offset)
      if (allocated(layout%block_neigh_list)) deallocate(layout%block_neigh_list)
      if (allocated(layout%pair_group_src)) deallocate(layout%pair_group_src)
      if (allocated(layout%pair_group_dst)) deallocate(layout%pair_group_dst)
      if (allocated(layout%pair_group_entry_offset)) deallocate(layout%pair_group_entry_offset)
      if (allocated(layout%entry_atom_i)) deallocate(layout%entry_atom_i)
      if (allocated(layout%entry_atom_j)) deallocate(layout%entry_atom_j)
      if (allocated(layout%entry_ih)) deallocate(layout%entry_ih)
      if (allocated(layout%entry_jslot)) deallocate(layout%entry_jslot)

      layout%ready = .false.
      layout%nblocks = 0
      layout%n_pair_groups = 0
      layout%n_entries = 0
   end subroutine destroy_macroblock_layout

   logical function macroblock_layout_is_ready(layout)
      type(macroblock_layout_type), intent(in), optional :: layout

      if (present(layout)) then
         macroblock_layout_is_ready = layout%ready
      else
         macroblock_layout_is_ready = ham_macroblock_layout%ready
      end if
   end function macroblock_layout_is_ready

   subroutine build_macroblock_layout(layout, Natom, nlist, nlistsize, aHam)
      type(macroblock_layout_type), intent(inout) :: layout
      integer, intent(in) :: Natom
      integer, dimension(:,:), intent(in) :: nlist
      integer, dimension(:), intent(in) :: nlistsize
      integer, dimension(:), intent(in) :: aHam

      integer :: bi, bj, iat, jat, ih, slot
      integer :: nblocks, n_pair_groups, n_entries, n_block_neigh
      integer :: group_id, atom_pos, neigh_pos, entry_pos
      integer, allocatable :: pair_counts(:,:), pair_group_index(:,:), cursor(:)

      call destroy_macroblock_layout(layout)

      if (do_macro_cells /= 'Y') return
      if (Num_macro <= 0) return
      if (.not.allocated(cell_index)) return
      if (.not.allocated(macro_nlistsize)) return
      if (.not.allocated(macro_atom_nlist)) return

      nblocks = Num_macro

      allocate(layout%atom_to_block(Natom))
      allocate(layout%block_atom_count(nblocks))
      allocate(layout%block_atom_offset(nblocks + 1))

      layout%atom_to_block = int(cell_index(1:Natom), c_int)
      layout%block_atom_count = int(macro_nlistsize(1:nblocks), c_int)
      layout%block_atom_offset(1) = 0_c_int
      do bi = 1, nblocks
         layout%block_atom_offset(bi + 1) = layout%block_atom_offset(bi) + layout%block_atom_count(bi)
      end do

      allocate(layout%block_atoms(Natom))
      atom_pos = 1
      do bi = 1, nblocks
         do slot = 1, macro_nlistsize(bi)
            layout%block_atoms(atom_pos) = int(macro_atom_nlist(bi, slot), c_int)
            atom_pos = atom_pos + 1
         end do
      end do

      allocate(pair_counts(nblocks, nblocks))
      allocate(pair_group_index(nblocks, nblocks))
      pair_counts = 0
      pair_group_index = 0

      do iat = 1, Natom
         ih = aHam(iat)
         if (ih <= 0 .or. ih > size(nlistsize)) cycle
         bi = cell_index(iat)
         if (bi <= 0 .or. bi > nblocks) cycle
         do slot = 1, nlistsize(ih)
            jat = nlist(slot, iat)
            if (jat <= 0 .or. jat > Natom) cycle
            bj = cell_index(jat)
            if (bj <= 0 .or. bj > nblocks) cycle
            pair_counts(bi, bj) = pair_counts(bi, bj) + 1
         end do
      end do

      n_entries = sum(pair_counts)
      n_pair_groups = count(pair_counts > 0)
      n_block_neigh = 0
      do bi = 1, nblocks
         n_block_neigh = n_block_neigh + count(pair_counts(bi, :) > 0)
      end do

      layout%nblocks = int(nblocks, c_int)
      layout%n_pair_groups = int(n_pair_groups, c_int)
      layout%n_entries = int(n_entries, c_int)

      if (n_pair_groups == 0 .or. n_entries == 0) then
         deallocate(pair_counts)
         deallocate(pair_group_index)
         return
      end if

      allocate(layout%block_neigh_count(nblocks))
      allocate(layout%block_neigh_offset(nblocks + 1))
      allocate(layout%block_neigh_list(n_block_neigh))
      allocate(layout%pair_group_src(n_pair_groups))
      allocate(layout%pair_group_dst(n_pair_groups))
      allocate(layout%pair_group_entry_offset(n_pair_groups + 1))
      allocate(layout%entry_atom_i(n_entries))
      allocate(layout%entry_atom_j(n_entries))
      allocate(layout%entry_ih(n_entries))
      allocate(layout%entry_jslot(n_entries))

      layout%block_neigh_offset(1) = 0_c_int
      do bi = 1, nblocks
         layout%block_neigh_count(bi) = int(count(pair_counts(bi, :) > 0), c_int)
         layout%block_neigh_offset(bi + 1) = layout%block_neigh_offset(bi) + layout%block_neigh_count(bi)
      end do

      group_id = 0
      neigh_pos = 1
      entry_pos = 0
      do bi = 1, nblocks
         do bj = 1, nblocks
            if (pair_counts(bi, bj) <= 0) cycle
            group_id = group_id + 1
            pair_group_index(bi, bj) = group_id
            layout%pair_group_src(group_id) = int(bi, c_int)
            layout%pair_group_dst(group_id) = int(bj, c_int)
            layout%pair_group_entry_offset(group_id) = int(entry_pos, c_int)
            entry_pos = entry_pos + pair_counts(bi, bj)
            layout%block_neigh_list(neigh_pos) = int(bj, c_int)
            neigh_pos = neigh_pos + 1
         end do
      end do
      layout%pair_group_entry_offset(n_pair_groups + 1) = int(entry_pos, c_int)

      allocate(cursor(n_pair_groups))
      do group_id = 1, n_pair_groups
         cursor(group_id) = int(layout%pair_group_entry_offset(group_id)) + 1
      end do

      do iat = 1, Natom
         ih = aHam(iat)
         if (ih <= 0 .or. ih > size(nlistsize)) cycle
         bi = cell_index(iat)
         if (bi <= 0 .or. bi > nblocks) cycle
         do slot = 1, nlistsize(ih)
            jat = nlist(slot, iat)
            if (jat <= 0 .or. jat > Natom) cycle
            bj = cell_index(jat)
            if (bj <= 0 .or. bj > nblocks) cycle
            group_id = pair_group_index(bi, bj)
            if (group_id <= 0) cycle
            entry_pos = cursor(group_id)
            layout%entry_atom_i(entry_pos) = int(iat, c_int)
            layout%entry_atom_j(entry_pos) = int(jat, c_int)
            layout%entry_ih(entry_pos) = int(ih, c_int)
            layout%entry_jslot(entry_pos) = int(slot, c_int)
            cursor(group_id) = cursor(group_id) + 1
         end do
      end do

      layout%ready = .true.

      deallocate(cursor)
      deallocate(pair_counts)
      deallocate(pair_group_index)
   end subroutine build_macroblock_layout

end module HamiltonianMacroBlocks
