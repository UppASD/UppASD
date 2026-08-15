!> RCG-09D: equivalence and scaling fixtures for the two adaptive
!> coarse-graining setup constructions that were sparse in storage but
!> quadratic in work.
!>
!> Both fixtures follow the RCG-07 pattern: the *previous* algorithm is
!> transcribed here as an independent brute-force reference, the production
!> construction is required to reproduce it, and a wall-clock bound is
!> asserted that the previous algorithm cannot meet.
!>
!> 1. `test_unique_bonds_match_discovered_list_reference` -- the record
!>    emission + stable-counting-sort fold in `fold_directed_bonds` must
!>    reproduce the discovered-list linear-scan reference: same canonical
!>    bond count, same bond numbering (first appearance in the emission
!>    sequence), same canonical atom pair per bond, all exactly, and the same
!>    folded 3x3 matrices to roundoff.
!> 2. `test_unique_bonds_reject_periodic_image_alias` -- the RCG-09A.1
!>    capability boundary still routes through the new grouping: a canonical
!>    pair carrying two directed entries in one direction has no unambiguous
!>    periodic image and must be rejected, not merged.
!> 3. `test_unique_bond_fold_scales_subquadratically` -- negative control on
!>    growth in directed bond count.
!> 4. `test_macroblock_layout_matches_dense_reference` -- the
!>    generation-stamped block-adjacency discovery must reproduce the dense
!>    per-source-block clear/scan reference in every CSR array.
!> 5. `test_macroblock_layout_scales_subquadratically` -- negative control on
!>    growth in block count at fixed atoms per block.
program test_adaptive_setup_scaling

   use, intrinsic :: iso_fortran_env, only : int64
   use Parameters, only : dblprec
   use Constants, only : mub
   use AdaptiveCGProduction, only : fold_directed_bonds, &
      ADAPTIVE_CG_PRODUCTION_OK, ADAPTIVE_CG_PRODUCTION_REJECTED
   use macrocells, only : Num_macro, cell_index, macro_nlistsize, macro_atom_nlist
   use HamiltonianMacroBlocks, only : macroblock_layout_type, build_macroblock_layout, &
      destroy_macroblock_layout

   implicit none

   ! A four-neighbour ring is used for both constructions.  The exchange and
   ! DMI offsets deliberately differ in their second shell so that the DMI
   ! pass discovers canonical bonds the exchange pass never saw; bond
   ! numbering by first appearance is then a real constraint, not a
   ! coincidence of a single shared neighbour list.
   integer, parameter :: exchange_offset(4) = (/1,-1,2,-2/)
   integer, parameter :: dmi_offset(4) = (/1,-1,3,-3/)

   integer :: failures

   failures = 0
   call test_unique_bonds_match_discovered_list_reference()
   call test_unique_bonds_reject_periodic_image_alias()
   call test_unique_bond_fold_scales_subquadratically()
   call test_macroblock_layout_matches_dense_reference()
   call test_macroblock_layout_scales_subquadratically()

   if (failures /= 0) then
      write(*,'(a,i0)') 'adaptive setup scaling tests failed: ',failures
      stop 1
   end if
   write(*,'(a)') 'adaptive setup scaling tests passed'

contains

   !---------------------------------------------------------------------------
   ! Unique bond construction
   !---------------------------------------------------------------------------

   subroutine test_unique_bonds_match_discovered_list_reference()
      integer, parameter :: sizes(3) = (/11,64,257/)
      integer, allocatable :: aham(:), nlistsize(:), nlist(:,:), dmlistsize(:), dmlist(:,:)
      real(dblprec), allocatable :: ncoup(:,:), dm_vect(:,:,:), moment(:)
      integer, allocatable :: bond_atom(:,:), reference_atom(:,:)
      real(dblprec), allocatable :: bond_matrix(:,:,:), reference_matrix(:,:,:)
      integer :: level, natom, bond_count, reference_count, status
      character(len=512) :: diagnostic

      do level = 1, size(sizes)
         natom = sizes(level)
         call build_ring_lists(natom,aham,nlistsize,nlist,ncoup,dmlistsize,dmlist, &
            dm_vect,moment)
         call reference_fold(natom,.true.,aham,nlistsize,nlist,ncoup,moment, &
            dmlistsize,dmlist,dm_vect,reference_count,reference_atom,reference_matrix)
         call fold_directed_bonds(natom,aham,nlistsize,nlist,ncoup,moment,bond_count, &
            bond_atom,bond_matrix,status,diagnostic,dmlistsize=dmlistsize, &
            dmlist=dmlist,dm_vect=dm_vect)
         call check(status == ADAPTIVE_CG_PRODUCTION_OK, &
            'exchange+DMI ring fold is accepted: '//trim(diagnostic))
         call compare_fold(natom,'exchange+DMI',bond_count,bond_atom,bond_matrix, &
            reference_count,reference_atom,reference_matrix)

         ! Same geometry without DMI: the exchange pass alone must still
         ! reproduce the reference, so the DMI-only bonds are the only
         ! difference between the two runs.
         call reference_fold(natom,.false.,aham,nlistsize,nlist,ncoup,moment, &
            dmlistsize,dmlist,dm_vect,reference_count,reference_atom,reference_matrix)
         call fold_directed_bonds(natom,aham,nlistsize,nlist,ncoup,moment,bond_count, &
            bond_atom,bond_matrix,status,diagnostic)
         call check(status == ADAPTIVE_CG_PRODUCTION_OK, &
            'exchange-only ring fold is accepted: '//trim(diagnostic))
         call compare_fold(natom,'exchange-only',bond_count,bond_atom,bond_matrix, &
            reference_count,reference_atom,reference_matrix)
      end do
   end subroutine test_unique_bonds_match_discovered_list_reference

   subroutine compare_fold(natom,label,bond_count,bond_atom,bond_matrix, &
      reference_count,reference_atom,reference_matrix)
      integer, intent(in) :: natom, bond_count, reference_count
      character(len=*), intent(in) :: label
      integer, intent(in) :: bond_atom(:,:), reference_atom(:,:)
      real(dblprec), intent(in) :: bond_matrix(:,:,:), reference_matrix(:,:,:)
      real(dblprec) :: scale, difference
      character(len=128) :: prefix

      write(prefix,'(a,a,i0,a)') trim(label),' ring(',natom,'): '
      call check(bond_count == reference_count, &
         trim(prefix)//'canonical bond count matches the discovered-list reference')
      if (bond_count /= reference_count) return
      call check(all(bond_atom(1,1:bond_count) == reference_atom(1,1:bond_count)) .and. &
         all(bond_atom(2,1:bond_count) == reference_atom(2,1:bond_count)), &
         trim(prefix)//'bond numbering and canonical atom pairs match the reference')
      ! The reference is a separate translation unit, and this project builds
      ! Release with -Ofast, so the two compilations of the same multiply chain
      ! are entitled to differ by an ulp.  Bitwise equality of the production
      ! fold against the previous production implementation is established
      ! binary-against-binary in docs/RCG-09D_SETUP_SCALING_EVIDENCE.md; here
      ! the numerical requirement is roundoff, which any misgrouped,
      ! mis-oriented, or reordered fold would exceed by orders of magnitude.
      scale = maxval(abs(reference_matrix(:,:,1:bond_count)))
      difference = maxval(abs(bond_matrix(:,:,1:bond_count)-reference_matrix(:,:,1:bond_count)))
      write(*,'(a,a,es12.5)') trim(prefix),'max relative matrix difference = ', &
         difference/max(scale,tiny(1.0_dblprec))
      call check(difference <= 16.0_dblprec*epsilon(1.0_dblprec)*scale, &
         trim(prefix)//'folded bond matrices match the reference to roundoff')
      call check((maxval(abs(bond_matrix(1,2,1:bond_count))) > 0.0_dblprec) .eqv. &
         (label == 'exchange+DMI'), &
         trim(prefix)//'off-diagonal DMI content is present exactly when DMI is folded')
   end subroutine compare_fold

   subroutine test_unique_bonds_reject_periodic_image_alias()
      integer, parameter :: natom = 24
      integer, allocatable :: aham(:), nlistsize(:), nlist(:,:), dmlistsize(:), dmlist(:,:)
      real(dblprec), allocatable :: ncoup(:,:), dm_vect(:,:,:), moment(:)
      integer, allocatable :: bond_atom(:,:)
      real(dblprec), allocatable :: bond_matrix(:,:,:)
      integer :: bond_count, status
      character(len=512) :: diagnostic

      call build_ring_lists(natom,aham,nlistsize,nlist,ncoup,dmlistsize,dmlist,dm_vect,moment)
      ! Atom 1 now reaches atom 2 through two directed exchange entries, which
      ! the sparse production list cannot disambiguate into distinct periodic
      ! images.  Slot 3 replaces the +2 neighbour, so the reverse direction
      ! still carries exactly one entry.
      nlist(3,1) = nlist(1,1)
      ncoup(3,1) = ncoup(1,1)
      call fold_directed_bonds(natom,aham,nlistsize,nlist,ncoup,moment,bond_count, &
         bond_atom,bond_matrix,status,diagnostic)
      call check(status == ADAPTIVE_CG_PRODUCTION_REJECTED .and. &
         index(diagnostic,'periodic-image alias') > 0, &
         'a doubled directed exchange entry is rejected as a periodic-image alias')

      call build_ring_lists(natom,aham,nlistsize,nlist,ncoup,dmlistsize,dmlist,dm_vect,moment)
      ncoup(1,1) = ncoup(1,1)+0.01_dblprec
      call fold_directed_bonds(natom,aham,nlistsize,nlist,ncoup,moment,bond_count, &
         bond_atom,bond_matrix,status,diagnostic)
      call check(status == ADAPTIVE_CG_PRODUCTION_REJECTED .and. &
         index(diagnostic,'reciprocal exchange contract') > 0, &
         'a non-reciprocal J_ij still fails the folded-pair contract after grouping')
   end subroutine test_unique_bonds_reject_periodic_image_alias

   !> Negative control. The discovered-list scan compared every accepted
   !> directed contribution against every canonical pair found so far, i.e.
   !> O(E^2) for E directed entries; the record/sort/fold construction is
   !> O(E + natom). A 16x increase in E predicts ~16x for the fixed algorithm
   !> and ~256x for the previous one. The 48x threshold sits between the two
   !> with margin for measurement noise.
   subroutine test_unique_bond_fold_scales_subquadratically()
      integer, parameter :: n_sizes = 3
      integer, parameter :: sizes(n_sizes) = (/4096,16384,65536/)
      integer, parameter :: trials = 3
      real(dblprec) :: elapsed(n_sizes), ratio
      integer :: level, trial

      do level = 1, n_sizes
         elapsed(level) = huge(1.0_dblprec)
         do trial = 1, trials
            elapsed(level) = min(elapsed(level),time_fold(sizes(level)))
         end do
      end do

      write(*,'(a)') 'RCG-09D unique-bond fold scaling (atoms, directed entries, min-of-trials seconds):'
      do level = 1, n_sizes
         write(*,'(a,i0,a,i0,a,es12.5)') '  atoms=',sizes(level), &
            ' directed=',8*sizes(level),' seconds=',elapsed(level)
      end do

      ratio = elapsed(n_sizes)/max(elapsed(1),tiny(1.0_dblprec))
      write(*,'(a,es12.5)') '  growth over a 16x directed-entry increase: ',ratio
      call check(ratio < 48.0_dblprec, &
         'unique-bond fold growth over 16x directed entries stays far below the quadratic prediction')
   end subroutine test_unique_bond_fold_scales_subquadratically

   real(dblprec) function time_fold(natom) result(seconds)
      integer, intent(in) :: natom
      integer, allocatable :: aham(:), nlistsize(:), nlist(:,:), dmlistsize(:), dmlist(:,:)
      real(dblprec), allocatable :: ncoup(:,:), dm_vect(:,:,:), moment(:)
      integer, allocatable :: bond_atom(:,:)
      real(dblprec), allocatable :: bond_matrix(:,:,:)
      integer :: bond_count, status
      integer(int64) :: started, stopped, rate
      character(len=512) :: diagnostic

      call build_ring_lists(natom,aham,nlistsize,nlist,ncoup,dmlistsize,dmlist,dm_vect,moment)
      call system_clock(started,rate)
      call fold_directed_bonds(natom,aham,nlistsize,nlist,ncoup,moment,bond_count, &
         bond_atom,bond_matrix,status,diagnostic,dmlistsize=dmlistsize,dmlist=dmlist, &
         dm_vect=dm_vect)
      call system_clock(stopped)
      call check(status == ADAPTIVE_CG_PRODUCTION_OK, &
         'scaling fixture fold is accepted: '//trim(diagnostic))
      seconds = real(stopped-started,dblprec)/real(rate,dblprec)
   end function time_fold

   !> Deterministic four-neighbour periodic ring with reciprocal couplings:
   !> J_ji = J_ij and D_ji = -D_ij, both keyed on the canonical atom pair.
   subroutine build_ring_lists(natom,aham,nlistsize,nlist,ncoup,dmlistsize,dmlist, &
      dm_vect,moment)
      integer, intent(in) :: natom
      integer, allocatable, intent(out) :: aham(:), nlistsize(:), nlist(:,:), &
         dmlistsize(:), dmlist(:,:)
      real(dblprec), allocatable, intent(out) :: ncoup(:,:), dm_vect(:,:,:), moment(:)
      integer :: atom, slot, target

      if (allocated(aham)) deallocate(aham,nlistsize,nlist,ncoup,dmlistsize,dmlist,dm_vect,moment)
      allocate(aham(natom),nlistsize(natom),nlist(4,natom),ncoup(4,natom), &
         dmlistsize(natom),dmlist(4,natom),dm_vect(3,4,natom),moment(natom))
      do atom = 1, natom
         aham(atom) = atom
         nlistsize(atom) = 4
         dmlistsize(atom) = 4
         moment(atom) = 1.0_dblprec+0.25_dblprec*real(mod(atom,5),dblprec)
         do slot = 1, 4
            target = wrapped_atom(atom+exchange_offset(slot),natom)
            nlist(slot,atom) = target
            ncoup(slot,atom) = pair_exchange(min(atom,target),max(atom,target))
            target = wrapped_atom(atom+dmi_offset(slot),natom)
            dmlist(slot,atom) = target
            dm_vect(:,slot,atom) = merge(1.0_dblprec,-1.0_dblprec,atom < target) * &
               pair_dmi(min(atom,target),max(atom,target))
         end do
      end do
   end subroutine build_ring_lists

   pure integer function wrapped_atom(index,natom)
      integer, intent(in) :: index, natom
      wrapped_atom = modulo(index-1,natom)+1
   end function wrapped_atom

   pure real(dblprec) function pair_exchange(low,high)
      integer, intent(in) :: low, high
      pair_exchange = 5.0d-2+1.0d-3*real(low,dblprec)-7.0d-4*real(high,dblprec)
   end function pair_exchange

   pure function pair_dmi(low,high) result(vector)
      integer, intent(in) :: low, high
      real(dblprec) :: vector(3)
      vector(1) = 3.7d-2-1.1d-3*real(low,dblprec)
      vector(2) = -2.9d-2+9.0d-4*real(high,dblprec)
      vector(3) = 6.1d-2+5.0d-4*real(low+high,dblprec)
   end function pair_dmi

   !> Brute-force reference: the discovered-list linear scan this task
   !> replaced, transcribed verbatim in its accepted folding convention
   !> (RCG-09A.1) and independent of the production implementation.
   subroutine reference_fold(natom,use_dm,aham,nlistsize,nlist,ncoup,moment, &
      dmlistsize,dmlist,dm_vect,bond,pair_atom,pair_matrix)
      integer, intent(in) :: natom, aham(:), nlistsize(:), nlist(:,:), dmlistsize(:), dmlist(:,:)
      logical, intent(in) :: use_dm
      real(dblprec), intent(in) :: ncoup(:,:), moment(:), dm_vect(:,:,:)
      integer, intent(out) :: bond
      integer, allocatable, intent(out) :: pair_atom(:,:)
      real(dblprec), allocatable, intent(out) :: pair_matrix(:,:,:)
      integer :: directed, atom, neighbour, target, iham, found
      real(dblprec) :: coefficient, orientation, dmi_energy_j(3), dmi_matrix(3,3)

      directed = sum(nlistsize(aham(1:natom)))
      if (use_dm) directed = directed+sum(dmlistsize(aham(1:natom)))
      if (allocated(pair_atom)) deallocate(pair_atom,pair_matrix)
      allocate(pair_atom(2,max(1,directed)),pair_matrix(3,3,max(1,directed)))
      pair_matrix = 0.0_dblprec
      bond = 0
      do atom = 1, natom
         iham = aham(atom)
         do neighbour = 1, nlistsize(iham)
            target = nlist(neighbour,atom)
            if (target < 1 .or. target > natom .or. target == atom) cycle
            found = 0
            if (bond > 0) then
               do found = 1, bond
                  if (pair_atom(1,found) == min(atom,target) .and. &
                      pair_atom(2,found) == max(atom,target)) exit
               end do
               if (found > bond) found = 0
            end if
            if (found == 0) then
               bond = bond+1
               found = bond
               pair_atom(1,found) = min(atom,target)
               pair_atom(2,found) = max(atom,target)
            end if
            coefficient = 0.5_dblprec*mub*moment(atom)*moment(target)*ncoup(neighbour,iham)
            pair_matrix(1,1,found) = pair_matrix(1,1,found)+coefficient
            pair_matrix(2,2,found) = pair_matrix(2,2,found)+coefficient
            pair_matrix(3,3,found) = pair_matrix(3,3,found)+coefficient
         end do
      end do
      if (use_dm) then
         do atom = 1, natom
            iham = aham(atom)
            do neighbour = 1, dmlistsize(iham)
               target = dmlist(neighbour,atom)
               if (target < 1 .or. target > natom .or. target == atom) cycle
               found = 0
               if (bond > 0) then
                  do found = 1, bond
                     if (pair_atom(1,found) == min(atom,target) .and. &
                         pair_atom(2,found) == max(atom,target)) exit
                  end do
                  if (found > bond) found = 0
               end if
               if (found == 0) then
                  bond = bond+1
                  found = bond
                  pair_atom(1,found) = min(atom,target)
                  pair_atom(2,found) = max(atom,target)
               end if
               dmi_energy_j = mub*moment(atom)*moment(target)*dm_vect(:,neighbour,iham)
               dmi_matrix = reshape((/0.0_dblprec,dmi_energy_j(3),-dmi_energy_j(2), &
                  -dmi_energy_j(3),0.0_dblprec,dmi_energy_j(1), &
                  dmi_energy_j(2),-dmi_energy_j(1),0.0_dblprec/),(/3,3/))
               orientation = merge(1.0_dblprec,-1.0_dblprec,atom < target)
               pair_matrix(:,:,found) = pair_matrix(:,:,found)+0.5_dblprec*orientation*dmi_matrix
            end do
         end do
      end if
   end subroutine reference_fold

   !---------------------------------------------------------------------------
   ! Hamiltonian macro-block connectivity
   !---------------------------------------------------------------------------

   subroutine test_macroblock_layout_matches_dense_reference()
      integer, parameter :: block_counts(3) = (/7,33,129/)
      integer, parameter :: block_atoms(2) = (/1,3/)
      type(macroblock_layout_type) :: layout
      integer, allocatable :: aham(:), nlistsize(:), nlist(:,:)
      integer, allocatable :: reference_pair_offset(:), reference_source(:), reference_target(:)
      integer, allocatable :: reference_entry_offset(:), reference_atom_offset(:), &
         reference_atom_entry_offset(:), reference_entry_source(:), reference_entry_target(:), &
         reference_entry_slot(:)
      integer :: level, variant, nblocks, per_block, natom
      character(len=128) :: prefix

      do level = 1, size(block_counts)
         do variant = 1, size(block_atoms)
            nblocks = block_counts(level)
            per_block = block_atoms(variant)
            natom = nblocks*per_block
            write(prefix,'(a,i0,a,i0,a)') 'macroblocks(',nblocks,'x',per_block,'): '
            call build_block_chain(nblocks,per_block,aham,nlistsize,nlist)
            call build_macroblock_layout(layout,natom,nlist,nlistsize,aham)
            call check(layout%ready,trim(prefix)//'production layout is built')
            if (.not. layout%ready) cycle
            call reference_macroblock_layout(natom,nblocks,per_block,aham,nlistsize,nlist, &
               reference_pair_offset,reference_source,reference_target,reference_entry_offset, &
               reference_atom_offset,reference_atom_entry_offset,reference_entry_source, &
               reference_entry_target,reference_entry_slot)

            call check(layout%n_block_pairs == size(reference_source), &
               trim(prefix)//'block-pair count matches the dense reference')
            call check(all(layout%source_block_pair_offset == reference_pair_offset), &
               trim(prefix)//'source-block CSR offsets match the dense reference')
            if (layout%n_block_pairs == size(reference_source)) then
               call check(all(layout%pair_source_block == reference_source), &
                  trim(prefix)//'pair source blocks match the dense reference')
               call check(all(layout%pair_destination_block == reference_target), &
                  trim(prefix)//'pair destination blocks match the dense reference, in ascending order')
               call check(all(layout%destination_block_for_pair == reference_target), &
                  trim(prefix)//'destination_block_for_pair matches the dense reference')
               call check(all(layout%block_pair_hamiltonian_entry_offset == reference_entry_offset), &
                  trim(prefix)//'block-pair entry offsets match the dense reference')
               call check(all(layout%block_pair_source_atom_offset == reference_atom_offset), &
                  trim(prefix)//'block-pair source-atom offsets match the dense reference')
            end if
            call check(layout%n_hamiltonian_entries == size(reference_entry_source), &
               trim(prefix)//'Hamiltonian entry count matches the dense reference')
            if (layout%n_hamiltonian_entries == size(reference_entry_source)) then
               call check(all(layout%block_pair_source_atom_entry_offset == reference_atom_entry_offset), &
                  trim(prefix)//'per-source-atom entry offsets match the dense reference')
               call check(all(layout%hamiltonian_entry_source_atom == reference_entry_source), &
                  trim(prefix)//'entry source atoms match the dense reference')
               call check(all(layout%hamiltonian_entry_destination_atom == reference_entry_target), &
                  trim(prefix)//'entry destination atoms match the dense reference')
               call check(all(layout%hamiltonian_entry_neighbour_slot == reference_entry_slot), &
                  trim(prefix)//'entry neighbour slots match the dense reference')
            end if
            call check(count(layout%pair_source_block /= layout%pair_destination_block) > 0, &
               trim(prefix)//'fixture exercises off-diagonal block connectivity')
            call destroy_macroblock_layout(layout)
         end do
      end do
      call release_macrocell_fixture()
   end subroutine test_macroblock_layout_matches_dense_reference

   !> Negative control. Every pass of the previous construction cleared or
   !> scanned an nblocks-length array once per source block, i.e. O(nblocks^2)
   !> at fixed atoms per block; the stamped construction is O(entries +
   !> nblocks). A 16x block-count increase predicts ~16x for the fixed
   !> algorithm and ~256x for the previous one.
   subroutine test_macroblock_layout_scales_subquadratically()
      integer, parameter :: n_sizes = 3
      integer, parameter :: sizes(n_sizes) = (/1024,4096,16384/)
      integer, parameter :: trials = 3
      real(dblprec) :: elapsed(n_sizes), ratio
      integer :: level, trial

      do level = 1, n_sizes
         elapsed(level) = huge(1.0_dblprec)
         do trial = 1, trials
            elapsed(level) = min(elapsed(level),time_macroblock_layout(sizes(level)))
         end do
      end do

      write(*,'(a)') 'RCG-09D macroblock layout scaling (blocks, atoms, min-of-trials seconds):'
      do level = 1, n_sizes
         write(*,'(a,i0,a,i0,a,es12.5)') '  blocks=',sizes(level), &
            ' atoms=',2*sizes(level),' seconds=',elapsed(level)
      end do

      ratio = elapsed(n_sizes)/max(elapsed(1),tiny(1.0_dblprec))
      write(*,'(a,es12.5)') '  growth over a 16x block-count increase: ',ratio
      call check(ratio < 48.0_dblprec, &
         'macroblock layout growth over 16x blocks stays far below the quadratic prediction')
      call release_macrocell_fixture()
   end subroutine test_macroblock_layout_scales_subquadratically

   real(dblprec) function time_macroblock_layout(nblocks) result(seconds)
      integer, intent(in) :: nblocks
      integer, parameter :: per_block = 2
      type(macroblock_layout_type) :: layout
      integer, allocatable :: aham(:), nlistsize(:), nlist(:,:)
      integer(int64) :: started, stopped, rate

      call build_block_chain(nblocks,per_block,aham,nlistsize,nlist)
      call system_clock(started,rate)
      call build_macroblock_layout(layout,nblocks*per_block,nlist,nlistsize,aham)
      call system_clock(stopped)
      call check(layout%ready,'macroblock scaling fixture is built')
      seconds = real(stopped-started,dblprec)/real(rate,dblprec)
      call destroy_macroblock_layout(layout)
   end function time_macroblock_layout

   !> A periodic chain of blocks with `per_block` atoms each, laid out
   !> contiguously, whose four-neighbour atomistic ring reaches at most the
   !> immediately adjacent blocks.  Block connectivity is therefore sparse:
   !> the pair count grows linearly in nblocks even though the previous
   !> construction's per-source-block work did not.
   subroutine build_block_chain(nblocks,per_block,aham,nlistsize,nlist)
      integer, intent(in) :: nblocks, per_block
      integer, allocatable, intent(out) :: aham(:), nlistsize(:), nlist(:,:)
      integer :: natom, atom, slot, block, local

      natom = nblocks*per_block
      if (allocated(aham)) deallocate(aham,nlistsize,nlist)
      allocate(aham(natom),nlistsize(natom),nlist(4,natom))
      do atom = 1, natom
         aham(atom) = atom
         nlistsize(atom) = 4
         do slot = 1, 4
            nlist(slot,atom) = wrapped_atom(atom+exchange_offset(slot),natom)
         end do
      end do
      call release_macrocell_fixture()
      Num_macro = nblocks
      allocate(cell_index(natom),macro_nlistsize(nblocks),macro_atom_nlist(nblocks,per_block))
      macro_nlistsize = per_block
      do block = 1, nblocks
         do local = 1, per_block
            atom = (block-1)*per_block+local
            cell_index(atom) = block
            macro_atom_nlist(block,local) = atom
         end do
      end do
   end subroutine build_block_chain

   subroutine release_macrocell_fixture()
      if (allocated(cell_index)) deallocate(cell_index)
      if (allocated(macro_nlistsize)) deallocate(macro_nlistsize)
      if (allocated(macro_atom_nlist)) deallocate(macro_atom_nlist)
      Num_macro = 0
   end subroutine release_macrocell_fixture

   !> Brute-force reference: the dense per-source-block clear/scan this task
   !> replaced, transcribed independently of the production implementation.
   subroutine reference_macroblock_layout(natom,nblocks,per_block,aham,nlistsize,nlist, &
      pair_offset,pair_source,pair_target,entry_offset,atom_offset,atom_entry_offset, &
      entry_source,entry_target,entry_slot)
      integer, intent(in) :: natom, nblocks, per_block, aham(:), nlistsize(:), nlist(:,:)
      integer, allocatable, intent(out) :: pair_offset(:), pair_source(:), pair_target(:), &
         entry_offset(:), atom_offset(:), atom_entry_offset(:), entry_source(:), &
         entry_target(:), entry_slot(:)
      integer, allocatable :: atom_to_block(:), atom_to_local(:), target_group(:), &
         pair_entry_count(:), local_entry_count(:), local_entry_cursor(:)
      logical, allocatable :: seen_target(:)
      integer :: bi, bj, iat, jat, ih, slot, group, entry, local, n_pairs, n_entries, n_pair_src_atoms

      if (allocated(pair_offset)) deallocate(pair_offset,pair_source,pair_target,entry_offset, &
         atom_offset,atom_entry_offset,entry_source,entry_target,entry_slot)
      allocate(atom_to_block(natom),atom_to_local(natom),seen_target(nblocks),target_group(nblocks))
      do iat = 1, natom
         atom_to_block(iat) = (iat-1)/per_block+1
         atom_to_local(iat) = mod(iat-1,per_block)
      end do

      allocate(pair_offset(nblocks+1))
      pair_offset(1) = 0
      n_pairs = 0
      do bi = 1, nblocks
         seen_target = .false.
         do iat = (bi-1)*per_block+1, bi*per_block
            ih = aham(iat)
            do slot = 1, nlistsize(ih)
               jat = nlist(slot,iat)
               if (jat < 1 .or. jat > natom) cycle
               seen_target(atom_to_block(jat)) = .true.
            end do
         end do
         n_pairs = n_pairs+count(seen_target)
         pair_offset(bi+1) = n_pairs
      end do

      allocate(pair_source(n_pairs),pair_target(n_pairs),entry_offset(n_pairs+1), &
         atom_offset(n_pairs+1),pair_entry_count(n_pairs))
      pair_entry_count = 0
      group = 0
      do bi = 1, nblocks
         seen_target = .false.
         do iat = (bi-1)*per_block+1, bi*per_block
            ih = aham(iat)
            do slot = 1, nlistsize(ih)
               jat = nlist(slot,iat)
               if (jat < 1 .or. jat > natom) cycle
               seen_target(atom_to_block(jat)) = .true.
            end do
         end do
         do bj = 1, nblocks
            if (.not. seen_target(bj)) cycle
            group = group+1
            pair_source(group) = bi
            pair_target(group) = bj
         end do
      end do

      do bi = 1, nblocks
         target_group = 0
         do group = pair_offset(bi)+1, pair_offset(bi+1)
            target_group(pair_target(group)) = group
         end do
         do iat = (bi-1)*per_block+1, bi*per_block
            ih = aham(iat)
            do slot = 1, nlistsize(ih)
               jat = nlist(slot,iat)
               if (jat < 1 .or. jat > natom) cycle
               bj = atom_to_block(jat)
               pair_entry_count(target_group(bj)) = pair_entry_count(target_group(bj))+1
            end do
         end do
      end do

      entry_offset(1) = 0
      do group = 1, n_pairs
         entry_offset(group+1) = entry_offset(group)+pair_entry_count(group)
      end do
      n_entries = entry_offset(n_pairs+1)
      n_pair_src_atoms = 0
      atom_offset(1) = 0
      do group = 1, n_pairs
         n_pair_src_atoms = n_pair_src_atoms+per_block
         atom_offset(group+1) = n_pair_src_atoms
      end do

      allocate(atom_entry_offset(n_pair_src_atoms+1),local_entry_count(n_pair_src_atoms), &
         local_entry_cursor(n_pair_src_atoms))
      local_entry_count = 0
      do bi = 1, nblocks
         target_group = 0
         do group = pair_offset(bi)+1, pair_offset(bi+1)
            target_group(pair_target(group)) = group
         end do
         do iat = (bi-1)*per_block+1, bi*per_block
            ih = aham(iat)
            do slot = 1, nlistsize(ih)
               jat = nlist(slot,iat)
               if (jat < 1 .or. jat > natom) cycle
               group = target_group(atom_to_block(jat))
               local = atom_offset(group)+atom_to_local(iat)+1
               local_entry_count(local) = local_entry_count(local)+1
            end do
         end do
      end do
      atom_entry_offset(1) = 0
      do group = 1, n_pairs
         entry = entry_offset(group)
         do local = 1, per_block
            atom_entry_offset(atom_offset(group)+local) = entry
            entry = entry+local_entry_count(atom_offset(group)+local)
         end do
      end do
      atom_entry_offset(n_pair_src_atoms+1) = n_entries
      local_entry_cursor = atom_entry_offset(1:n_pair_src_atoms)+1

      allocate(entry_source(n_entries),entry_target(n_entries),entry_slot(n_entries))
      do bi = 1, nblocks
         target_group = 0
         do group = pair_offset(bi)+1, pair_offset(bi+1)
            target_group(pair_target(group)) = group
         end do
         do iat = (bi-1)*per_block+1, bi*per_block
            ih = aham(iat)
            do slot = 1, nlistsize(ih)
               jat = nlist(slot,iat)
               if (jat < 1 .or. jat > natom) cycle
               group = target_group(atom_to_block(jat))
               local = atom_offset(group)+atom_to_local(iat)+1
               entry = local_entry_cursor(local)
               entry_source(entry) = iat
               entry_target(entry) = jat
               entry_slot(entry) = slot
               local_entry_cursor(local) = entry+1
            end do
         end do
      end do
      deallocate(atom_to_block,atom_to_local,seen_target,target_group,pair_entry_count, &
         local_entry_count,local_entry_cursor)
   end subroutine reference_macroblock_layout

   subroutine check(condition,label)
      logical, intent(in) :: condition
      character(len=*), intent(in) :: label
      if (condition) then
         write(*,'(a,a)') 'PASS: ',trim(label)
      else
         write(*,'(a,a)') 'FAIL: ',trim(label)
         failures = failures+1
      end if
   end subroutine check

end program test_adaptive_setup_scaling
