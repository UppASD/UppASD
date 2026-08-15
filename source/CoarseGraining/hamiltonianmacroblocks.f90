!-------------------------------------------------------------------------------
! MODULE: HamiltonianMacroBlocks
!> @brief Backend-neutral spatial layout derived from the macrocell decomposition.
!>
!> The layout deliberately does not renumber atoms or modify Hamiltonian arrays.
!> It records macrocell membership and groups the finalized exchange-neighbour
!> entries by (source block, destination block), with an additional span for
!> every source atom.  CPU, HIP/CUDA, FFT, or future reordered-storage backends
!> can therefore use the same metadata without changing the physical ordering
!> selected by the regular or SFC neighbour mapper.  The optional canonical
!> topology argument adapts the same interaction metadata to BlockTopology;
!> omitting it preserves the active legacy macrocell path.
!-------------------------------------------------------------------------------
module HamiltonianMacroBlocks

   use Parameters
   use macrocells, only : Num_macro, cell_index, macro_nlistsize, macro_atom_nlist
   use BlockTopology, only : block_topology_type

   implicit none
   private

   type, public :: macroblock_layout_type
      logical :: ready = .false.
      integer :: nblocks = 0
      integer :: n_block_pairs = 0
      integer :: n_hamiltonian_entries = 0
      integer :: max_block_atoms = 0
      ! Atom membership and block-local position, indexed by global atom.
      integer, allocatable :: atom_to_block(:)
      integer, allocatable :: atom_to_local(:)
      ! CSR block-to-atom membership: atoms for block b are in
      ! block_atoms(block_atom_offset(b)+1:block_atom_offset(b+1)).
      integer, allocatable :: block_atom_count(:)
      integer, allocatable :: block_atom_offset(:)
      integer, allocatable :: block_atoms(:)
      ! CSR source-block-to-block-pair adjacency.  A pair p describes the
      ! exchange-active directed relation pair_source_block(p) -> pair_destination_block(p).
      integer, allocatable :: source_block_pair_offset(:)
      integer, allocatable :: destination_block_for_pair(:)
      integer, allocatable :: pair_source_block(:)
      integer, allocatable :: pair_destination_block(:)
      ! CSR block-pair-to-Hamiltonian-entry and block-pair-to-source-atom maps.
      integer, allocatable :: block_pair_hamiltonian_entry_offset(:)
      integer, allocatable :: block_pair_source_atom_offset(:)
      ! CSR source atom to Hamiltonian entries, repeated once for each block pair.
      integer, allocatable :: block_pair_source_atom_entry_offset(:)
      ! Atom-level exchange-neighbour entries.  The reduced index and neighbour
      ! slot address ham%nlistsize and ham%ncoup without duplicating couplings.
      integer, allocatable :: hamiltonian_entry_source_atom(:)
      integer, allocatable :: hamiltonian_entry_destination_atom(:)
      integer, allocatable :: hamiltonian_entry_reduced_index(:)
      integer, allocatable :: hamiltonian_entry_neighbour_slot(:)
      integer, allocatable :: hamiltonian_entry_source_local_atom(:)
      integer, allocatable :: hamiltonian_entry_destination_local_atom(:)
   end type macroblock_layout_type

   type(macroblock_layout_type), public, save :: ham_macroblock_layout

   public :: build_macroblock_layout, destroy_macroblock_layout, macroblock_layout_is_ready

contains

   subroutine destroy_macroblock_layout(layout)
      type(macroblock_layout_type), intent(inout) :: layout

      if (allocated(layout%atom_to_block)) deallocate(layout%atom_to_block)
      if (allocated(layout%atom_to_local)) deallocate(layout%atom_to_local)
      if (allocated(layout%block_atom_count)) deallocate(layout%block_atom_count)
      if (allocated(layout%block_atom_offset)) deallocate(layout%block_atom_offset)
      if (allocated(layout%block_atoms)) deallocate(layout%block_atoms)
      if (allocated(layout%source_block_pair_offset)) deallocate(layout%source_block_pair_offset)
      if (allocated(layout%destination_block_for_pair)) deallocate(layout%destination_block_for_pair)
      if (allocated(layout%pair_source_block)) deallocate(layout%pair_source_block)
      if (allocated(layout%pair_destination_block)) deallocate(layout%pair_destination_block)
      if (allocated(layout%block_pair_hamiltonian_entry_offset)) deallocate(layout%block_pair_hamiltonian_entry_offset)
      if (allocated(layout%block_pair_source_atom_offset)) deallocate(layout%block_pair_source_atom_offset)
      if (allocated(layout%block_pair_source_atom_entry_offset)) deallocate(layout%block_pair_source_atom_entry_offset)
      if (allocated(layout%hamiltonian_entry_source_atom)) deallocate(layout%hamiltonian_entry_source_atom)
      if (allocated(layout%hamiltonian_entry_destination_atom)) deallocate(layout%hamiltonian_entry_destination_atom)
      if (allocated(layout%hamiltonian_entry_reduced_index)) deallocate(layout%hamiltonian_entry_reduced_index)
      if (allocated(layout%hamiltonian_entry_neighbour_slot)) deallocate(layout%hamiltonian_entry_neighbour_slot)
      if (allocated(layout%hamiltonian_entry_source_local_atom)) deallocate(layout%hamiltonian_entry_source_local_atom)
      if (allocated(layout%hamiltonian_entry_destination_local_atom)) deallocate(layout%hamiltonian_entry_destination_local_atom)
      layout%ready = .false.
      layout%nblocks = 0
      layout%n_block_pairs = 0
      layout%n_hamiltonian_entries = 0
      layout%max_block_atoms = 0
   end subroutine destroy_macroblock_layout

   logical function macroblock_layout_is_ready(layout)
      type(macroblock_layout_type), intent(in), optional :: layout
      if (present(layout)) then
         macroblock_layout_is_ready = layout%ready
      else
         macroblock_layout_is_ready = ham_macroblock_layout%ready
      end if
   end function macroblock_layout_is_ready

   subroutine build_macroblock_layout(layout, Natom, nlist, nlistsize, aHam, topology)
      type(macroblock_layout_type), intent(inout) :: layout
      integer, intent(in) :: Natom
      integer, dimension(:,:), intent(in) :: nlist
      integer, dimension(:), intent(in) :: nlistsize, aHam
      type(block_topology_type), intent(in), optional :: topology

      integer :: bi, bj, iat, jat, ih, slot, group, entry, atom_pos, local
      integer :: nblocks, n_pairs, n_entries, n_pair_src_atoms
      integer :: generation, key, running, moved, cursor
      integer, allocatable :: pair_entry_count(:)
      integer, allocatable :: local_entry_count(:), local_entry_cursor(:), target_group(:)
      integer, allocatable :: target_stamp(:), pair_cursor(:)
      integer, allocatable :: discovered_source(:), discovered_target(:), order(:), bucket(:)

      call destroy_macroblock_layout(layout)
      if (Natom <= 0 .or. size(aHam) < Natom) return
      if (present(topology)) then
         if (.not. topology%ready .or. topology%n_atoms /= Natom) return
         nblocks = topology%n_spatial_blocks
      else
         if (Num_macro <= 0) return
         if (.not. allocated(cell_index) .or. .not. allocated(macro_nlistsize) .or. &
             .not. allocated(macro_atom_nlist)) return
         if (size(cell_index) < Natom) return
         nblocks = Num_macro
      end if
      allocate(layout%atom_to_block(Natom), layout%atom_to_local(Natom))
      allocate(layout%block_atom_count(nblocks), layout%block_atom_offset(nblocks + 1))
      layout%atom_to_local = -1
      if (present(topology)) then
         layout%atom_to_block = topology%atom_to_block
         layout%block_atom_count = topology%block_atom_count
         layout%block_atom_offset = topology%block_atom_offset
      else
         layout%atom_to_block = cell_index(1:Natom)
         layout%block_atom_count = macro_nlistsize(1:nblocks)
         layout%block_atom_offset(1) = 0
         do bi = 1, nblocks
            layout%block_atom_offset(bi + 1) = layout%block_atom_offset(bi) + layout%block_atom_count(bi)
         end do
      end if
      if (layout%block_atom_offset(nblocks + 1) /= Natom) then
         call destroy_macroblock_layout(layout)
         return
      end if

      allocate(layout%block_atoms(Natom))
      atom_pos = 1
      do bi = 1, nblocks
         do local = 1, layout%block_atom_count(bi)
            if (present(topology)) then
               iat = topology%block_atoms(topology%block_atom_offset(bi) + local)
            else
               iat = macro_atom_nlist(bi, local)
            end if
            if (iat < 1 .or. iat > Natom) then
               call destroy_macroblock_layout(layout)
               return
            end if
            layout%block_atoms(atom_pos) = iat
            layout%atom_to_block(iat) = bi
            layout%atom_to_local(iat) = local - 1
            atom_pos = atom_pos + 1
         end do
      end do

      ! RCG-09D: block-pair discovery is driven by a generation-stamped marker
      ! array.  A source block's cost is proportional to the number of directed
      ! neighbour entries it owns plus the number of target blocks it actually
      ! touches; no pass clears or scans an nblocks-length array per source
      ! block, so sparse block connectivity no longer costs O(nblocks^2).
      allocate(layout%source_block_pair_offset(nblocks + 1), target_stamp(nblocks), &
         target_group(nblocks))
      target_stamp = 0
      generation = 0
      layout%source_block_pair_offset(1) = 0
      n_pairs = 0
      do bi = 1, nblocks
         generation = generation + 1
         do atom_pos = layout%block_atom_offset(bi) + 1, layout%block_atom_offset(bi + 1)
            iat = layout%block_atoms(atom_pos)
            ih = aHam(iat)
            if (ih < 1 .or. ih > size(nlistsize)) cycle
            do slot = 1, nlistsize(ih)
               jat = nlist(slot, iat)
               if (jat < 1 .or. jat > Natom) cycle
               bj = layout%atom_to_block(jat)
               if (bj < 1 .or. bj > nblocks) cycle
               if (target_stamp(bj) == generation) cycle
               target_stamp(bj) = generation
               n_pairs = n_pairs + 1
            end do
         end do
         layout%source_block_pair_offset(bi + 1) = n_pairs
      end do

      allocate(layout%destination_block_for_pair(n_pairs), layout%pair_source_block(n_pairs), &
         layout%pair_destination_block(n_pairs), layout%block_pair_hamiltonian_entry_offset(n_pairs + 1), &
         layout%block_pair_source_atom_offset(n_pairs + 1), pair_entry_count(n_pairs))
      allocate(discovered_source(max(1, n_pairs)), discovered_target(max(1, n_pairs)), &
         order(max(1, n_pairs)), bucket(nblocks), pair_cursor(nblocks))
      pair_entry_count = 0
      cursor = 0
      do bi = 1, nblocks
         generation = generation + 1
         do atom_pos = layout%block_atom_offset(bi) + 1, layout%block_atom_offset(bi + 1)
            iat = layout%block_atoms(atom_pos)
            ih = aHam(iat)
            if (ih < 1 .or. ih > size(nlistsize)) cycle
            do slot = 1, nlistsize(ih)
               jat = nlist(slot, iat)
               if (jat < 1 .or. jat > Natom) cycle
               bj = layout%atom_to_block(jat)
               if (bj < 1 .or. bj > nblocks) cycle
               if (target_stamp(bj) == generation) cycle
               target_stamp(bj) = generation
               cursor = cursor + 1
               discovered_source(cursor) = bi
               discovered_target(cursor) = bj
            end do
         end do
      end do

      ! Pairs are discovered in neighbour-entry order, but the CSR contract
      ! orders every source block's destinations ascending.  One counting sort
      ! by destination block, replayed into the per-source CSR cursors, restores
      ! that order in O(n_pairs + nblocks) instead of an nblocks scan per source.
      bucket = 0
      do cursor = 1, n_pairs
         bucket(discovered_target(cursor)) = bucket(discovered_target(cursor)) + 1
      end do
      running = 1
      do key = 1, nblocks
         moved = bucket(key)
         bucket(key) = running
         running = running + moved
      end do
      do cursor = 1, n_pairs
         key = discovered_target(cursor)
         order(bucket(key)) = cursor
         bucket(key) = bucket(key) + 1
      end do
      pair_cursor(1:nblocks) = layout%source_block_pair_offset(1:nblocks)
      do cursor = 1, n_pairs
         entry = order(cursor)
         bi = discovered_source(entry)
         pair_cursor(bi) = pair_cursor(bi) + 1
         group = pair_cursor(bi)
         layout%pair_source_block(group) = bi
         layout%pair_destination_block(group) = discovered_target(entry)
         layout%destination_block_for_pair(group) = discovered_target(entry)
      end do

      do bi = 1, nblocks
         generation = generation + 1
         do group = layout%source_block_pair_offset(bi) + 1, layout%source_block_pair_offset(bi + 1)
            bj = layout%pair_destination_block(group)
            target_stamp(bj) = generation
            target_group(bj) = group
         end do
         do atom_pos = layout%block_atom_offset(bi) + 1, layout%block_atom_offset(bi + 1)
            iat = layout%block_atoms(atom_pos)
            ih = aHam(iat)
            if (ih < 1 .or. ih > size(nlistsize)) cycle
            do slot = 1, nlistsize(ih)
               jat = nlist(slot, iat)
               if (jat < 1 .or. jat > Natom) cycle
               bj = layout%atom_to_block(jat)
               if (bj < 1 .or. bj > nblocks) cycle
               if (target_stamp(bj) /= generation) cycle
               pair_entry_count(target_group(bj)) = pair_entry_count(target_group(bj)) + 1
            end do
         end do
      end do

      layout%block_pair_hamiltonian_entry_offset(1) = 0
      do group = 1, n_pairs
         layout%block_pair_hamiltonian_entry_offset(group + 1) = &
            layout%block_pair_hamiltonian_entry_offset(group) + pair_entry_count(group)
      end do
      n_entries = layout%block_pair_hamiltonian_entry_offset(n_pairs + 1)
      n_pair_src_atoms = 0
      layout%block_pair_source_atom_offset(1) = 0
      do group = 1, n_pairs
         n_pair_src_atoms = n_pair_src_atoms + layout%block_atom_count(layout%pair_source_block(group))
         layout%block_pair_source_atom_offset(group + 1) = n_pair_src_atoms
      end do

      allocate(layout%block_pair_source_atom_entry_offset(n_pair_src_atoms + 1), &
         local_entry_count(n_pair_src_atoms), local_entry_cursor(n_pair_src_atoms))
      local_entry_count = 0
      do bi = 1, nblocks
         generation = generation + 1
         do group = layout%source_block_pair_offset(bi) + 1, layout%source_block_pair_offset(bi + 1)
            bj = layout%pair_destination_block(group)
            target_stamp(bj) = generation
            target_group(bj) = group
         end do
         do atom_pos = layout%block_atom_offset(bi) + 1, layout%block_atom_offset(bi + 1)
            iat = layout%block_atoms(atom_pos)
            ih = aHam(iat)
            if (ih < 1 .or. ih > size(nlistsize)) cycle
            do slot = 1, nlistsize(ih)
               jat = nlist(slot, iat)
               if (jat < 1 .or. jat > Natom) cycle
               bj = layout%atom_to_block(jat)
               if (bj < 1 .or. bj > nblocks) cycle
               if (target_stamp(bj) /= generation) cycle
               group = target_group(bj)
               local = layout%block_pair_source_atom_offset(group) + layout%atom_to_local(iat) + 1
               local_entry_count(local) = local_entry_count(local) + 1
            end do
         end do
      end do
      layout%block_pair_source_atom_entry_offset(1) = 0
      do group = 1, n_pairs
         entry = layout%block_pair_hamiltonian_entry_offset(group)
         do local = 1, layout%block_atom_count(layout%pair_source_block(group))
            atom_pos = layout%block_pair_source_atom_offset(group) + local
            layout%block_pair_source_atom_entry_offset(atom_pos) = entry
            entry = entry + local_entry_count(atom_pos)
         end do
      end do
      layout%block_pair_source_atom_entry_offset(n_pair_src_atoms + 1) = n_entries
      local_entry_cursor = layout%block_pair_source_atom_entry_offset(1:n_pair_src_atoms) + 1

      allocate(layout%hamiltonian_entry_source_atom(n_entries), layout%hamiltonian_entry_destination_atom(n_entries), &
         layout%hamiltonian_entry_reduced_index(n_entries), layout%hamiltonian_entry_neighbour_slot(n_entries), &
         layout%hamiltonian_entry_source_local_atom(n_entries), layout%hamiltonian_entry_destination_local_atom(n_entries))
      do bi = 1, nblocks
         generation = generation + 1
         do group = layout%source_block_pair_offset(bi) + 1, layout%source_block_pair_offset(bi + 1)
            bj = layout%pair_destination_block(group)
            target_stamp(bj) = generation
            target_group(bj) = group
         end do
         do atom_pos = layout%block_atom_offset(bi) + 1, layout%block_atom_offset(bi + 1)
            iat = layout%block_atoms(atom_pos)
            ih = aHam(iat)
            if (ih < 1 .or. ih > size(nlistsize)) cycle
            do slot = 1, nlistsize(ih)
               jat = nlist(slot, iat)
               if (jat < 1 .or. jat > Natom) cycle
               bj = layout%atom_to_block(jat)
               if (bj < 1 .or. bj > nblocks) cycle
               if (target_stamp(bj) /= generation) cycle
               group = target_group(bj)
               local = layout%block_pair_source_atom_offset(group) + layout%atom_to_local(iat) + 1
               entry = local_entry_cursor(local)
               layout%hamiltonian_entry_source_atom(entry) = iat
               layout%hamiltonian_entry_destination_atom(entry) = jat
               layout%hamiltonian_entry_reduced_index(entry) = ih
               layout%hamiltonian_entry_neighbour_slot(entry) = slot
               layout%hamiltonian_entry_source_local_atom(entry) = layout%atom_to_local(iat)
               layout%hamiltonian_entry_destination_local_atom(entry) = layout%atom_to_local(jat)
               local_entry_cursor(local) = entry + 1
            end do
         end do
      end do

      layout%nblocks = nblocks
      layout%n_block_pairs = n_pairs
      layout%n_hamiltonian_entries = n_entries
      layout%max_block_atoms = maxval(layout%block_atom_count)
      layout%ready = .true.

      deallocate(target_stamp, target_group, pair_entry_count, local_entry_count, local_entry_cursor)
      deallocate(discovered_source, discovered_target, order, bucket, pair_cursor)
   end subroutine build_macroblock_layout

end module HamiltonianMacroBlocks
