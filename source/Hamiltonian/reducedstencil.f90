!-------------------------------------------------------------------------------
! MODULE: ReducedStencil
!> @brief Validated translational stencil representation for reduced scalar J.
!-------------------------------------------------------------------------------
module ReducedStencil

   use Parameters, only : dblprec

   implicit none
   private

   type, public :: reduced_stencil_record
      integer :: output_basis = 0
      integer :: input_basis = 0
      integer :: delta_cell(3) = 0
      integer :: cell_offset = 0
      real(dblprec) :: j = 0.0_dblprec
   end type reduced_stencil_record

   type, public :: reduced_stencil_t
      integer :: na = 0
      integer :: n1 = 0
      integer :: n2 = 0
      integer :: n3 = 0
      integer, allocatable :: record_start(:)
      type(reduced_stencil_record), allocatable :: record(:)
   end type reduced_stencil_t

   public :: reduced_stencil_eligible
   public :: build_reduced_stencil
   public :: apply_reduced_stencil
   public :: apply_reduced_stencil_target
   public :: clear_reduced_stencil
   public :: cell_basis_to_atom
   public :: atom_to_cell_basis
   public :: wrap_reduced_cell

contains

   ! Eligibility is deliberately narrower than the general production
   ! Hamiltonian.  The compact representation is exact only when every target
   ! is a translated copy of a basis target and the source map is periodic in
   ! all three directions.
   logical function reduced_stencil_eligible(natom,na,n1,n2,n3,bc1,bc2,bc3, &
         do_reduced,do_ralloy,nham,ham_index,diagnostic)
      implicit none

      integer, intent(in) :: natom, na, n1, n2, n3, do_ralloy, nham
      character(len=1), intent(in) :: bc1, bc2, bc3, do_reduced
      integer, intent(in) :: ham_index(:)
      character(len=*), intent(out), optional :: diagnostic

      reduced_stencil_eligible=.false.
      if (present(diagnostic)) diagnostic=''

      if (do_reduced /= 'Y') then
         call reject('do_reduced is not Y',diagnostic)
         return
      endif
      if (do_ralloy /= 0) then
         call reject('random-alloy Hamiltonians are not translation invariant',diagnostic)
         return
      endif
      if (bc1 /= 'P' .or. bc2 /= 'P' .or. bc3 /= 'P') then
         call reject('all three boundary conditions must be periodic',diagnostic)
         return
      endif
      if (na < 1 .or. n1 < 1 .or. n2 < 1 .or. n3 < 1) then
         call reject('lattice and basis dimensions must be positive',diagnostic)
         return
      endif
      if (natom /= na*n1*n2*n3) then
         call reject('atom layout is not a regular replicated cell',diagnostic)
         return
      endif
      if (nham /= na) then
         call reject('reduced Hamiltonian size is not the basis size',diagnostic)
         return
      endif
      if (size(ham_index) < natom) then
         call reject('Hamiltonian lookup does not cover all physical atoms',diagnostic)
         return
      endif
      if (any(ham_index(1:natom) < 1) .or. any(ham_index(1:natom) > nham)) then
         call reject('Hamiltonian lookup leaves the reduced basis',diagnostic)
         return
      endif

      reduced_stencil_eligible=.true.
      if (present(diagnostic)) diagnostic='eligible'

   contains

      subroutine reject(message,result)
         character(len=*), intent(in) :: message
         character(len=*), intent(out), optional :: result
         if (present(result)) result=message
      end subroutine reject

   end function reduced_stencil_eligible


   ! Derive the stencil from the production reduced arrays.  The first cell is
   ! the source of the compact records; every other cell is then checked
   ! against it before the representation is accepted.
   subroutine build_reduced_stencil(natom,na,n1,n2,n3,bc1,bc2,bc3,do_reduced, &
         do_ralloy,nham,ham_index,nlistsize,nlist,ncoup,stencil,ok,diagnostic)
      implicit none

      integer, intent(in) :: natom, na, n1, n2, n3, do_ralloy, nham
      character(len=1), intent(in) :: bc1, bc2, bc3, do_reduced
      integer, intent(in) :: ham_index(:), nlistsize(:), nlist(:,:)
      real(dblprec), intent(in) :: ncoup(:,:,:)
      type(reduced_stencil_t), intent(out) :: stencil
      logical, intent(out), optional :: ok
      character(len=*), intent(out), optional :: diagnostic

      integer :: b, atom, source, slot, nrec, z, start
      integer :: cell(3), source_cell(3), source_basis
      integer :: wrapped(3), expected
      logical :: success
      character(len=256) :: reason

      call clear_reduced_stencil(stencil)
      success=.false.
      reason=''

      if (.not.reduced_stencil_eligible(natom,na,n1,n2,n3,bc1,bc2,bc3, &
            do_reduced,do_ralloy,nham,ham_index,reason)) then
         call finish(success,reason,ok,diagnostic)
         return
      endif
      if (size(nlistsize) < nham .or. size(nlist,2) < natom) then
         reason='production reduced neighbour arrays have insufficient dimensions'
         call finish(success,reason,ok,diagnostic)
         return
      endif
      if (size(ncoup,1) < 1 .or. size(ncoup,2) < nham .or. size(ncoup,3) /= 1) then
         reason='production scalar-J coupling array has insufficient dimensions'
         call finish(success,reason,ok,diagnostic)
         return
      endif

      nrec=0
      do b=1,na
         if (ham_index(b) < 1 .or. ham_index(b) > nham) then
            reason='basis target does not have a valid reduced Hamiltonian index'
            call finish(success,reason,ok,diagnostic)
            return
         endif
         z=nlistsize(ham_index(b))
         if (z < 0 .or. z > size(nlist,1) .or. z > size(ncoup,1)) then
            reason='basis neighbour count is outside production array bounds'
            call finish(success,reason,ok,diagnostic)
            return
         endif
         nrec=nrec+z
      end do

      stencil%na=na
      stencil%n1=n1
      stencil%n2=n2
      stencil%n3=n3
      allocate(stencil%record_start(na+1),stencil%record(nrec))
      nrec=0
      do b=1,na
         stencil%record_start(b)=nrec+1
         atom=cell_basis_to_atom(0,0,0,b,na,n1,n2,n3)
         z=nlistsize(ham_index(atom))
         do slot=1,z
            source=nlist(slot,atom)
            if (source < 1 .or. source > natom) then
               reason='production neighbour index is outside physical atom range'
               call clear_reduced_stencil(stencil)
               call finish(success,reason,ok,diagnostic)
               return
            endif
            call atom_to_cell_basis(source,na,n1,n2,n3,source_cell,source_basis)
            nrec=nrec+1
            stencil%record(nrec)%output_basis=b
            stencil%record(nrec)%input_basis=source_basis
            stencil%record(nrec)%delta_cell=source_cell
            stencil%record(nrec)%cell_offset=source_cell(1)+n1*(source_cell(2)+n2*source_cell(3))
            stencil%record(nrec)%j=ncoup(slot,ham_index(atom),1)
         end do
      end do
      stencil%record_start(na+1)=nrec+1

      ! The reduced production path retains full nlist records for every target
      ! but obtains the count/couplings from the compact basis entry.  Check
      ! both that lookup and every translated source index are invariant.
      do atom=1,natom
         call atom_to_cell_basis(atom,na,n1,n2,n3,cell,b)
         if (ham_index(atom) /= ham_index(b)) then
            reason='Hamiltonian lookup is not translation invariant'
            call clear_reduced_stencil(stencil)
            call finish(success,reason,ok,diagnostic)
            return
         endif
         start=stencil%record_start(b)
         z=stencil%record_start(b+1)-start
         if (nlistsize(ham_index(atom)) /= z .and. cell(1)==0 .and. &
               cell(2)==0 .and. cell(3)==0) then
            reason='basis neighbour count is inconsistent with the stencil'
            call clear_reduced_stencil(stencil)
            call finish(success,reason,ok,diagnostic)
            return
         endif
         do slot=1,z
            source_basis=stencil%record(start+slot-1)%input_basis
            call wrap_reduced_cell(cell,stencil%record(start+slot-1)%delta_cell, &
               n1,n2,n3,wrapped)
            expected=cell_basis_to_atom(wrapped(1),wrapped(2),wrapped(3),source_basis, &
               na,n1,n2,n3)
            if (nlist(slot,atom) /= expected) then
               reason='translation-invariance validation failed for nlist'
               call clear_reduced_stencil(stencil)
               call finish(success,reason,ok,diagnostic)
               return
            endif
            if (ham_index(nlist(slot,atom)) < 1 .or. ham_index(nlist(slot,atom)) > nham) then
               reason='translated source has an invalid reduced Hamiltonian index'
               call clear_reduced_stencil(stencil)
               call finish(success,reason,ok,diagnostic)
               return
            endif
         end do
      end do

      success=.true.
      reason='stencil built and translation invariance validated'
      call finish(success,reason,ok,diagnostic)

   contains

      subroutine finish(success,message,ok_out,diagnostic_out)
         logical, intent(in) :: success
         character(len=*), intent(in) :: message
         logical, intent(out), optional :: ok_out
         character(len=*), intent(out), optional :: diagnostic_out
         if (present(ok_out)) ok_out=success
         if (present(diagnostic_out)) diagnostic_out=message
      end subroutine finish

   end subroutine build_reduced_stencil


   subroutine apply_reduced_stencil(stencil,spin,field)
      implicit none

      type(reduced_stencil_t), intent(in) :: stencil
      real(dblprec), intent(in) :: spin(:,:)
      real(dblprec), intent(out) :: field(:,:)

      integer :: atom, b, slot, start, stop
      integer :: cell(3), wrapped(3), source
      real(dblprec) :: coupling

      if (size(spin,1) /= 3 .or. size(field,1) /= 3 .or. &
            size(spin,2) /= size(field,2) .or. &
            size(spin,2) /= stencil%na*stencil%n1*stencil%n2*stencil%n3) then
         error stop 'apply_reduced_stencil: incompatible spin/field dimensions'
      endif
      if (.not.allocated(stencil%record_start) .or. .not.allocated(stencil%record)) then
         error stop 'apply_reduced_stencil: stencil is not allocated'
      endif

      field=0.0_dblprec
      do atom=1,size(spin,2)
         call atom_to_cell_basis(atom,stencil%na,stencil%n1,stencil%n2,stencil%n3,cell,b)
         start=stencil%record_start(b)
         stop=stencil%record_start(b+1)-1
         do slot=start,stop
            call wrap_reduced_cell(cell,stencil%record(slot)%delta_cell, &
               stencil%n1,stencil%n2,stencil%n3,wrapped)
            source=cell_basis_to_atom(wrapped(1),wrapped(2),wrapped(3), &
               stencil%record(slot)%input_basis,stencil%na,stencil%n1,stencil%n2,stencil%n3)
            coupling=stencil%record(slot)%j
            field(:,atom)=field(:,atom)+coupling*spin(:,source)
         end do
      end do
   end subroutine apply_reduced_stencil


   ! Production target kernel.  The stencil is built once during Hamiltonian
   ! setup; only the target cell and compact records are touched here.  The
   ! source atom is reconstructed from basis/cell coordinates instead of
   ! loading an arbitrary physical neighbour index.
   subroutine apply_reduced_stencil_target(stencil,atom,ensemble,spin,field)
      implicit none

      type(reduced_stencil_t), intent(in) :: stencil
      integer, intent(in) :: atom, ensemble
      real(dblprec), intent(in) :: spin(:,:,:)
      real(dblprec), intent(inout) :: field(:)

      integer :: zero_based, cell_number, cell1, cell2, cell3, basis
      integer :: slot, start, stop, source, source_cell1, source_cell2, source_cell3
      integer :: input_basis, delta1, delta2, delta3
      real(dblprec) :: bx, by, bz, coupling

      ! These checks are outside the stencil loop.  Production callers have
      ! already validated the dimensions during setup.
      if (size(spin,1) /= 3 .or. size(field) /= 3 .or. &
            atom < 1 .or. atom > stencil%na*stencil%n1*stencil%n2*stencil%n3 .or. &
            ensemble < 1 .or. ensemble > size(spin,3)) then
         error stop 'apply_reduced_stencil_target: incompatible dimensions'
      endif

      zero_based=atom-1
      basis=modulo(zero_based,stencil%na)+1
      cell_number=zero_based/stencil%na
      cell1=modulo(cell_number,stencil%n1)
      cell_number=cell_number/stencil%n1
      cell2=modulo(cell_number,stencil%n2)
      cell3=cell_number/stencil%n2

      start=stencil%record_start(basis)
      stop=stencil%record_start(basis+1)-1
      bx=field(1)
      by=field(2)
      bz=field(3)
      do slot=start,stop
         input_basis=stencil%record(slot)%input_basis
         delta1=stencil%record(slot)%delta_cell(1)
         delta2=stencil%record(slot)%delta_cell(2)
         delta3=stencil%record(slot)%delta_cell(3)
         source_cell1=cell1+delta1
         source_cell2=cell2+delta2
         source_cell3=cell3+delta3
         source=input_basis+stencil%na*(cell1+stencil%n1*(cell2+stencil%n2*cell3)+ &
            stencil%record(slot)%cell_offset)
         if (source_cell1 >= stencil%n1) source=source-stencil%na*stencil%n1
         if (source_cell2 >= stencil%n2) source=source-stencil%na*stencil%n1*stencil%n2
         if (source_cell3 >= stencil%n3) source=source-stencil%na*stencil%n1*stencil%n2*stencil%n3
         coupling=stencil%record(slot)%j
         bx=bx+coupling*spin(1,source,ensemble)
         by=by+coupling*spin(2,source,ensemble)
         bz=bz+coupling*spin(3,source,ensemble)
      end do
      field(1)=bx
      field(2)=by
      field(3)=bz
   end subroutine apply_reduced_stencil_target


   pure integer function cell_basis_to_atom(cell1,cell2,cell3,basis,na,n1,n2,n3)
      implicit none
      integer, intent(in) :: cell1,cell2,cell3,basis,na,n1,n2,n3

      if (basis < 1 .or. basis > na .or. cell1 < 0 .or. cell1 >= n1 .or. &
            cell2 < 0 .or. cell2 >= n2 .or. cell3 < 0 .or. cell3 >= n3) then
         error stop 'cell_basis_to_atom: index outside replicated cell'
      endif
      cell_basis_to_atom=basis+na*(cell1+n1*(cell2+n2*cell3))
   end function cell_basis_to_atom


   pure subroutine atom_to_cell_basis(atom,na,n1,n2,n3,cell,basis)
      implicit none
      integer, intent(in) :: atom,na,n1,n2,n3
      integer, intent(out) :: cell(3),basis
      integer :: zero_based, cell_number

      if (atom < 1 .or. atom > na*n1*n2*n3) then
         error stop 'atom_to_cell_basis: atom outside replicated cell'
      endif
      zero_based=atom-1
      basis=modulo(zero_based,na)+1
      cell_number=zero_based/na
      cell(1)=modulo(cell_number,n1)
      cell_number=cell_number/n1
      cell(2)=modulo(cell_number,n2)
      cell(3)=cell_number/n2
   end subroutine atom_to_cell_basis


   pure subroutine wrap_reduced_cell(cell,delta,n1,n2,n3,wrapped)
      implicit none
      integer, intent(in) :: cell(3),delta(3),n1,n2,n3
      integer, intent(out) :: wrapped(3)

      wrapped(1)=modulo(cell(1)+delta(1),n1)
      wrapped(2)=modulo(cell(2)+delta(2),n2)
      wrapped(3)=modulo(cell(3)+delta(3),n3)
   end subroutine wrap_reduced_cell


   subroutine clear_reduced_stencil(stencil)
      implicit none
      type(reduced_stencil_t), intent(out) :: stencil
      stencil%na=0
      stencil%n1=0
      stencil%n2=0
      stencil%n3=0
      if (allocated(stencil%record_start)) deallocate(stencil%record_start)
      if (allocated(stencil%record)) deallocate(stencil%record)
   end subroutine clear_reduced_stencil

end module ReducedStencil
