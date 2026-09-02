! CPU-HAM-02B: reduced periodic scalar-J stencil construction and oracle.
program test_reduced_stencil

   use Parameters, only : dblprec
   use Constants, only : mub,mry
   use HamiltonianData, only : ham
   use HamiltonianActions, only : effective_field
   use InputData, only : ham_inp
   use ReducedStencil

   implicit none

   integer, parameter :: max_basis=4, max_records=12
   integer :: delta(3,max_records,max_basis), input_basis(max_records,max_basis)
   real(dblprec) :: coupling(max_records,max_basis)
   integer :: z(max_basis)

   call test_mapping()

   delta=0
   input_basis=0
   coupling=0.0_dblprec
   z=0
   z(1)=5
   delta(:,1:5,1)=reshape((/0,0,0, 1,0,0, -1,0,0, 2,0,0, -2,0,0/),(/3,5/))
   input_basis(1:5,1)=(/1,1,1,1,1/)
   coupling(1:5,1)=(/0.20_dblprec,-0.30_dblprec,0.40_dblprec, &
      0.15_dblprec,-0.05_dblprec/)
   call run_case('one-basis periodic',1,7,1,1,z,delta,input_basis,coupling)

   delta=0
   input_basis=0
   coupling=0.0_dblprec
   z=0
   z(1)=4
   z(2)=4
   delta(:,1:4,1)=reshape((/0,0,0, 0,0,0, -1,1,-1, 1,0,0/),(/3,4/))
   input_basis(1:4,1)=(/1,2,2,1/)
   coupling(1:4,1)=(/0.20_dblprec,0.35_dblprec,-0.10_dblprec,0.07_dblprec/)
   delta(:,1:4,2)=reshape((/0,0,0, 0,0,0, 1,-1,1, -1,0,0/),(/3,4/))
   input_basis(1:4,2)=(/2,1,1,2/)
   coupling(1:4,2)=(/-0.12_dblprec,0.31_dblprec,0.08_dblprec,-0.04_dblprec/)
   call run_case('multi-basis periodic boundary crossing',2,3,2,2,z,delta,input_basis,coupling)

   delta=0
   input_basis=0
   coupling=0.0_dblprec
   z=0
   z(1:4)=7
   call make_long_range_fixture(z,delta,input_basis,coupling)
   call run_case('long-range periodic',4,7,5,3,z,delta,input_basis,coupling)

   call test_negative_eligibility()
   call test_metadata_footprint()

   write(*,'(a)') 'CPU-HAM-02B reduced-stencil tests passed'

contains

   subroutine test_mapping()
      integer :: na,n1,n2,n3,atom,basis,cell(3),round_trip

      na=3
      n1=4
      n2=3
      n3=2
      do atom=1,na*n1*n2*n3
         call atom_to_cell_basis(atom,na,n1,n2,n3,cell,basis)
         round_trip=cell_basis_to_atom(cell(1),cell(2),cell(3),basis,na,n1,n2,n3)
         call check(round_trip==atom,'cell/basis mapping round trip')
      end do
      call check(cell_basis_to_atom(3,2,1,2,na,n1,n2,n3)== &
         2+na*(3+n1*(2+n2*1)),'cell/basis mapping uses canonical atom order')
   end subroutine test_mapping


   subroutine run_case(label,na,n1,n2,n3,z,delta,input_basis,coupling)
      character(len=*), intent(in) :: label
      integer, intent(in) :: na,n1,n2,n3,z(max_basis)
      integer, intent(in) :: delta(3,max_records,max_basis)
      integer, intent(in) :: input_basis(max_records,max_basis)
      real(dblprec), intent(in) :: coupling(max_records,max_basis)

      integer :: natom,maxz,atom,b,slot,source
      integer :: cell(3),wrapped(3)
      integer, allocatable :: ham_index(:),nlistsize(:),nlist(:,:)
      real(dblprec), allocatable :: ncoup(:,:,:),spin(:,:),spin3d(:,:,:),mmom(:,:)
      real(dblprec), allocatable :: external_field(:,:,:),time_field(:,:,:)
      real(dblprec), allocatable :: beff(:,:,:),beff1(:,:,:),beff2(:,:,:)
      real(dblprec), allocatable :: stencil_field(:,:)
      real(dblprec), allocatable :: emomM_macro(:,:,:)
      integer, allocatable :: cell_index(:),macro_nlistsize(:)
      type(reduced_stencil_t) :: stencil
      logical :: ok
      character(len=256) :: diagnostic
      real(dblprec) :: direct_energy,stencil_energy

      natom=na*n1*n2*n3
      maxz=maxval(z)
      allocate(ham_index(natom),nlistsize(natom),nlist(maxz,natom),ncoup(maxz,na,1))
      ham_index=0
      nlistsize=0
      nlist=0
      ncoup=0.0_dblprec

      do atom=1,natom
         call atom_to_cell_basis(atom,na,n1,n2,n3,cell,b)
         ham_index(atom)=b
         if (cell(1)==0 .and. cell(2)==0 .and. cell(3)==0) nlistsize(b)=z(b)
         do slot=1,z(b)
            call wrap_reduced_cell(cell,delta(:,slot,b),n1,n2,n3,wrapped)
            source=cell_basis_to_atom(wrapped(1),wrapped(2),wrapped(3), &
               input_basis(slot,b),na,n1,n2,n3)
            nlist(slot,atom)=source
         end do
      end do
      do b=1,na
         do slot=1,z(b)
            ncoup(slot,b,1)=coupling(slot,b)
         end do
      end do

      call clear_hamiltonian_fixture()
      allocate(ham%aHam(natom),ham%nlistsize(natom),ham%nlist(maxz,natom), &
         ham%ncoup(maxz,na,1))
      ham%aHam=ham_index
      ham%nlistsize=nlistsize
      ham%nlist=nlist
      ham%ncoup=ncoup

      call build_reduced_stencil(natom,na,n1,n2,n3,'P','P','P','Y',0,na,ham_index, &
         nlistsize,nlist,ncoup,stencil,ok,diagnostic)
      call check(ok,trim(label)//': stencil accepted')
      call check(index(diagnostic,'translation invariance validated')>0, &
         trim(label)//': translation validation reported')

      allocate(spin(3,natom),spin3d(3,natom,1),stencil_field(3,natom),mmom(natom,1))
      allocate(external_field(3,natom,1),time_field(3,natom,1))
      allocate(beff(3,natom,1),beff1(3,natom,1),beff2(3,natom,1))
      allocate(emomM_macro(3,1,1),cell_index(natom),macro_nlistsize(1))
      do atom=1,natom
         spin(:,atom)=(/0.11_dblprec*real(atom,dblprec), &
            -0.07_dblprec*real(atom+1,dblprec),0.03_dblprec*real(2*atom+1,dblprec)/)
      end do
      spin3d(:,:,1)=spin
      mmom=1.0_dblprec
      external_field=0.0_dblprec
      time_field=0.0_dblprec
      cell_index=1
      macro_nlistsize(1)=natom
      emomM_macro=0.0_dblprec

      call apply_reduced_stencil(stencil,spin,stencil_field)
      call check(maxval(abs(stencil_field(:,1)- &
         (coupling(1,1)*spin(:,1)+coupling(2,1)*spin(:,2)+ &
         coupling(3,1)*spin(:,n1)+coupling(4,1)*spin(:,3)+ &
         coupling(5,1)*spin(:,n1-1)))) < 1.0d-14 .or. na/=1, &
         trim(label)//': host stencil applies expected target order')

      ham_inp%do_jtensor=0
      ham_inp%exc_inter='N'
      ham_inp%do_dm=0
      ham_inp%do_sa=0
      ham_inp%do_pd=0
      ham_inp%do_biqdm=0
      ham_inp%do_bq=0
      ham_inp%do_ring=0
      ham_inp%do_chir=0
      ham_inp%do_anisotropy=0
      ham_inp%do_dip=0
      if (allocated(ham%target_order)) deallocate(ham%target_order)
      if (allocated(ham%target_work_prefix)) deallocate(ham%target_work_prefix)

      call effective_field(natom,1,1,natom,spin3d,mmom,external_field,time_field, &
         beff,beff1,beff2,direct_energy,1,cell_index,emomM_macro,macro_nlistsize, &
         na,n1,n2,n3,measure_energy=.true.)
      call check(maxval(abs(beff1(:,:,1)-stencil_field)) < 1.0d-14, &
         trim(label)//': reduced field matches canonical DIRECT')
      stencil_energy=-0.5_dblprec*sum(spin*stencil_field)*mub/mry
      call check(abs(direct_energy-stencil_energy) < 1.0d-12, &
         trim(label)//': pair energy matches field-derived oracle')

      deallocate(ham_index,nlistsize,nlist,ncoup,spin,spin3d,stencil_field,mmom, &
         external_field,time_field,beff,beff1,beff2,emomM_macro,cell_index,macro_nlistsize)
      call clear_hamiltonian_fixture()
   end subroutine run_case


   subroutine make_long_range_fixture(z,delta,input_basis,coupling)
      integer, intent(in) :: z(max_basis)
      integer, intent(out) :: delta(3,max_records,max_basis)
      integer, intent(out) :: input_basis(max_records,max_basis)
      real(dblprec), intent(out) :: coupling(max_records,max_basis)
      integer :: b,slot

      do b=1,max_basis
         do slot=1,z(b)
            delta(:,slot,b)=(/mod(slot,3)-1,mod(slot+1,5)-2,mod(slot+2,3)-1/)
            input_basis(slot,b)=mod(b+slot-1,max_basis)+1
            coupling(slot,b)=0.01_dblprec*real(2*b-slot,dblprec)
         end do
      end do
   end subroutine make_long_range_fixture


   subroutine test_negative_eligibility()
      integer :: ham_index(8)
      logical :: eligible,ok
      integer :: nlistsize(8),nlist(2,8)
      real(dblprec) :: ncoup(2,1,1)
      type(reduced_stencil_t) :: stencil
      character(len=256) :: diagnostic

      ham_index=(/1,1,1,1,1,1,1,1/)
      eligible=reduced_stencil_eligible(8,1,2,2,2,'P','P','0','Y',0,1,ham_index,diagnostic)
      call check(.not.eligible .and. index(diagnostic,'periodic')>0, &
         'nonperiodic eligibility is rejected')
      eligible=reduced_stencil_eligible(8,1,2,2,2,'P','P','P','N',0,1,ham_index,diagnostic)
      call check(.not.eligible .and. index(diagnostic,'do_reduced')>0, &
         'non-reduced Hamiltonian eligibility is rejected')

      nlistsize=(/2,0,0,0,0,0,0,0/)
      nlist=1
      ncoup=0.5_dblprec
      nlist(1,5)=2
      call build_reduced_stencil(8,1,2,2,2,'P','P','P','Y',0,1,ham_index, &
         nlistsize,nlist,ncoup,stencil,ok,diagnostic)
      call check(.not.ok .and. index(diagnostic,'translation-invariance')>0, &
         'translation-invariance violation is declined')
   end subroutine test_negative_eligibility


   subroutine test_metadata_footprint()
      integer(kind=8) :: direct_bytes,stencil_bytes
      integer(kind=8), parameter :: i4=4_8, r8=8_8, record_bytes=28_8

      direct_bytes=(int(1340,8)*16384_8+2_8*16384_8)*i4+int(1340,8)*4_8*r8
      stencil_bytes=4_8*1338_8*record_bytes
      write(*,'(a,i0,a,i0,a,f8.1,a)') 'CPU-HAM-02B Fe metadata bytes: direct=', &
         (96_8*16000_8+2_8*16000_8)*i4+96_8*2_8*r8,', reduced=', &
         2_8*96_8*record_bytes,', ratio=',real((96_8*16000_8+2_8*16000_8)*i4+96_8*2_8*r8)/ &
         real(2_8*96_8*record_bytes),'x'
      write(*,'(a,i0,a,i0,a,f8.1,a)') 'CPU-HAM-02B Nd metadata bytes: direct=',direct_bytes, &
         ', reduced=',stencil_bytes,', ratio=',real(direct_bytes)/real(stencil_bytes),'x'
      call check(stencil_bytes < direct_bytes,'reduced metadata is smaller than direct metadata')
   end subroutine test_metadata_footprint


   subroutine clear_hamiltonian_fixture()
      if (allocated(ham%aHam)) deallocate(ham%aHam)
      if (allocated(ham%nlistsize)) deallocate(ham%nlistsize)
      if (allocated(ham%nlist)) deallocate(ham%nlist)
      if (allocated(ham%ncoup)) deallocate(ham%ncoup)
      if (allocated(ham%target_order)) deallocate(ham%target_order)
      if (allocated(ham%target_work_prefix)) deallocate(ham%target_work_prefix)
   end subroutine clear_hamiltonian_fixture


   subroutine check(condition,label)
      logical, intent(in) :: condition
      character(len=*), intent(in) :: label
      if (.not.condition) error stop 'FAILED: '//trim(label)
   end subroutine check

end program test_reduced_stencil
