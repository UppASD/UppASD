! CPU-HAM-04B: persistent scalar-J periodic CPU convolution backend parity.
program test_cpu_convolution

   use Parameters, only : dblprec
   use HamiltonianData, only : ham
   use HamiltonianActions, only : effective_field, setup_convolution_backend, &
      cleanup_convolution_backend, convolution_backend_can_apply, convolution_backend_get_stats
   use InputData, only : ham_inp, do_convolution, do_sparse
   use ReducedStencil
   use CPUConvolution

   implicit none

   call test_delta_and_uniform()
   call test_multibasis_and_ensembles()
   call test_long_range_nd_fixture()
   call test_production_backend()
   call test_ineligible_cases()
   write(*,'(a)') 'CPU-HAM-04B CPU convolution backend tests passed'

contains

   subroutine test_delta_and_uniform()
      type(reduced_stencil_t) :: stencil
      type(cpu_convolution_t) :: convolution
      real(dblprec), allocatable :: spin(:,:,:),field(:,:,:),direct(:,:,:)
      character(len=256) :: diagnostic
      logical :: ok
      integer :: atom,axis

      call make_stencil(stencil,1,5,2,1,1)
      call put_record(stencil,1,1,1,(/1,0,0/),2.5_dblprec)
      ok=cpu_convolution_init(convolution,stencil,1,diagnostic)
      call check(ok,trim(diagnostic))
      ok=cpu_convolution_build_kernel(convolution,stencil,diagnostic)
      call check(ok,trim(diagnostic))
      allocate(spin(3,10,1),field(3,10,1),direct(3,10,1))

      do atom=1,10
         do axis=1,3
            spin(axis,atom,1)=0.07_dblprec*real(axis+2*atom,dblprec)-0.11_dblprec*real(atom,dblprec)
         end do
      end do
      ok=cpu_convolution_apply(convolution,spin,field,diagnostic)
      call check(ok,trim(diagnostic))
      call apply_direct(stencil,spin,direct)
      call check(maxval(abs(field-direct)) < 2.0e-12_dblprec, &
         'delta convolution matches circular DIRECT')
      call check(abs(field(1,1,1)-2.5_dblprec*spin(1,2,1)) < 2.0e-12_dblprec, &
         'delta kernel has the target-to-source sign')

      spin=1.0_dblprec
      ok=cpu_convolution_apply(convolution,spin,field,diagnostic)
      call check(ok,trim(diagnostic))
      call check(maxval(abs(field-2.5_dblprec)) < 2.0e-12_dblprec, &
         'uniform state preserves periodic normalization')
      call check(abs(cpu_convolution_pair_energy(spin,field)+37.5_dblprec) < 2.0e-12_dblprec, &
         'pair energy uses -1/2 field relation')

      deallocate(spin,field,direct)
      call cpu_convolution_clear(convolution)
      call clear_reduced_stencil(stencil)
   end subroutine test_delta_and_uniform


   subroutine test_multibasis_and_ensembles()
      type(reduced_stencil_t) :: stencil
      type(cpu_convolution_t) :: convolution
      real(dblprec), allocatable :: spin(:,:,:),field(:,:,:),direct(:,:,:)
      character(len=256) :: diagnostic
      logical :: ok
      integer :: atom,axis,ensemble

      call make_stencil(stencil,2,3,2,2,4)
      call put_record(stencil,1,1,1,(/0,0,0/),0.50_dblprec)
      call put_record(stencil,2,1,2,(/1,0,0/),1.20_dblprec)
      call put_record(stencil,3,2,1,(/-1,0,0/),-0.70_dblprec)
      call put_record(stencil,4,2,2,(/0,0,0/),0.25_dblprec)
      ok=cpu_convolution_init(convolution,stencil,3,diagnostic)
      call check(ok,trim(diagnostic))
      ok=cpu_convolution_build_kernel(convolution,stencil,diagnostic)
      call check(ok,trim(diagnostic))
      call check(size(convolution%kernel_spectral)==convolution%spectral_cells*4, &
         'multi-basis stores exactly N_A^2 spectral kernels')
      allocate(spin(3,24,3),field(3,24,3),direct(3,24,3))
      do ensemble=1,3
         do atom=1,24
            do axis=1,3
               spin(axis,atom,ensemble)=0.013_dblprec*real(axis+atom+7*ensemble,dblprec) &
                  -0.021_dblprec*real(modulo(atom+axis+ensemble,5),dblprec)
            end do
         end do
      end do
      ok=cpu_convolution_apply(convolution,spin,field,diagnostic)
      call check(ok,trim(diagnostic))
      call apply_direct(stencil,spin,direct)
      call check(maxval(abs(field-direct)) < 3.0e-12_dblprec, &
         'multi-basis multi-ensemble convolution matches DIRECT')

      deallocate(spin,field,direct)
      call cpu_convolution_clear(convolution)
      call clear_reduced_stencil(stencil)
   end subroutine test_multibasis_and_ensembles


   subroutine test_long_range_nd_fixture()
      integer, parameter :: na=4,n1=7,n2=5,n3=3,nrecords=7
      type(reduced_stencil_t) :: stencil
      type(cpu_convolution_t) :: convolution
      real(dblprec), allocatable :: spin(:,:,:),field(:,:,:),direct(:,:,:)
      character(len=256) :: diagnostic
      logical :: ok
      integer :: a,slot,atom,axis

      call make_stencil(stencil,na,n1,n2,n3,na*nrecords)
      do a=1,na
         do slot=1,nrecords
            call put_record(stencil,(a-1)*nrecords+slot,a,modulo(a+slot-1,na)+1, &
               (/modulo(slot,3)-1,modulo(slot+1,5)-2,modulo(slot+2,3)-1/), &
               0.01_dblprec*real(2*a-slot,dblprec))
         end do
      end do
      ok=cpu_convolution_init(convolution,stencil,2,diagnostic)
      call check(ok,trim(diagnostic))
      ok=cpu_convolution_build_kernel(convolution,stencil,diagnostic)
      call check(ok,trim(diagnostic))
      allocate(spin(3,na*n1*n2*n3,2),field(3,na*n1*n2*n3,2), &
         direct(3,na*n1*n2*n3,2))
      do atom=1,size(spin,2)
         do axis=1,3
            spin(axis,atom,1)=sin(0.17_dblprec*real(axis+atom,dblprec))
            spin(axis,atom,2)=cos(0.11_dblprec*real(2*axis+atom,dblprec))
         end do
      end do
      ok=cpu_convolution_apply(convolution,spin,field,diagnostic)
      call check(ok,trim(diagnostic))
      call apply_direct(stencil,spin,direct)
      call check(maxval(abs(field-direct)) < 5.0e-12_dblprec, &
         'Nd-derived long-range convolution matches DIRECT')

      deallocate(spin,field,direct)
      call cpu_convolution_clear(convolution)
      call clear_reduced_stencil(stencil)
   end subroutine test_long_range_nd_fixture


   subroutine test_production_backend()
      integer, parameter :: na=1,n1=5,n2=2,n3=1,natom=10,maxz=3,ensembles=2
      integer :: ham_index(natom),nlistsize(natom),nlist(maxz,natom)
      real(dblprec) :: ncoup(maxz,na,1)
      real(dblprec), allocatable :: spin(:,:,:),mmom(:,:),external_field(:,:,:)
      real(dblprec), allocatable :: time_field(:,:,:),beff(:,:,:),beff1(:,:,:)
      real(dblprec), allocatable :: beff2(:,:,:),direct(:,:,:)
      real(dblprec), allocatable :: emomM_macro(:,:,:)
      integer, allocatable :: cell_index(:),macro_nlistsize(:)
      type(reduced_stencil_t) :: stencil
      character(len=256) :: diagnostic
      real(dblprec) :: direct_energy,convolution_energy,setup_seconds,pack_seconds
      real(dblprec) :: forward_seconds,spectral_seconds,inverse_seconds,unpack_seconds,apply_seconds
      integer(kind=8) :: apply_count
      integer :: atom,axis,ensemble,slot,cell(3),wrapped(3)
      logical :: ok

      ham_index=1
      nlistsize=maxz
      ncoup(:,1,1)=(/0.4_dblprec,-0.2_dblprec,0.1_dblprec/)
      do atom=1,natom
         call atom_to_cell_basis(atom,na,n1,n2,n3,cell,slot)
         call wrap_reduced_cell(cell,(/0,0,0/),n1,n2,n3,wrapped)
         nlist(1,atom)=cell_basis_to_atom(wrapped(1),wrapped(2),wrapped(3),1,na,n1,n2,n3)
         call wrap_reduced_cell(cell,(/1,0,0/),n1,n2,n3,wrapped)
         nlist(2,atom)=cell_basis_to_atom(wrapped(1),wrapped(2),wrapped(3),1,na,n1,n2,n3)
         call wrap_reduced_cell(cell,(/-1,0,0/),n1,n2,n3,wrapped)
         nlist(3,atom)=cell_basis_to_atom(wrapped(1),wrapped(2),wrapped(3),1,na,n1,n2,n3)
      end do

      call clear_production_fixture()
      allocate(ham%aHam(natom),ham%nlistsize(natom),ham%nlist(maxz,natom), &
         ham%ncoup(maxz,na,1))
      ham%aHam=ham_index
      ham%nlistsize=nlistsize
      ham%nlist=nlist
      ham%ncoup=ncoup
      call build_reduced_stencil(natom,na,n1,n2,n3,'P','P','P','Y',0,na,ham_index, &
         nlistsize,nlist,ncoup,stencil,ok,diagnostic)
      call check(ok,trim(diagnostic))
      ham%reduced_stencil=stencil

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
      do_sparse='N'
      do_convolution='N'
      allocate(spin(3,natom,ensembles),mmom(natom,ensembles),external_field(3,natom,ensembles), &
         time_field(3,natom,ensembles),beff(3,natom,ensembles),beff1(3,natom,ensembles), &
         beff2(3,natom,ensembles),direct(3,natom,ensembles),emomM_macro(3,1,ensembles), &
         cell_index(natom),macro_nlistsize(1))
      do ensemble=1,ensembles
         do atom=1,natom
            do axis=1,3
               spin(axis,atom,ensemble)=0.03_dblprec*real(axis+2*atom+5*ensemble,dblprec)- &
                  0.01_dblprec*real(modulo(axis+atom,4),dblprec)
            end do
         end do
      end do
      mmom=1.0_dblprec
      external_field=0.0_dblprec
      time_field=0.0_dblprec
      emomM_macro=0.0_dblprec
      cell_index=1
      macro_nlistsize(1)=natom

      call clear_reduced_stencil(ham%reduced_stencil)
      call effective_field(natom,ensembles,1,natom,spin,mmom,external_field,time_field, &
         beff,beff1,beff2,direct_energy,1,cell_index,emomM_macro,macro_nlistsize, &
         na,n1,n2,n3,measure_energy=.true.)
      direct=beff1

      do_convolution='Y'
      ham%reduced_stencil=stencil
      call setup_convolution_backend(natom,ensembles,0,'N',na,n1,n2,n3,'P','P','P','Y',na)
      call check(convolution_backend_can_apply(natom,ensembles,1,natom), &
         'production convolution backend is active')
      call effective_field(natom,ensembles,1,natom,spin,mmom,external_field,time_field, &
         beff,beff1,beff2,convolution_energy,1,cell_index,emomM_macro,macro_nlistsize, &
         na,n1,n2,n3,measure_energy=.true.)
      call check(maxval(abs(beff1-direct)) < 2.0e-12_dblprec, &
         'production convolution field matches DIRECT')
      call check(abs(convolution_energy-direct_energy) < 2.0e-12_dblprec, &
         'production convolution energy matches DIRECT')
      call effective_field(natom,ensembles,1,natom,spin,mmom,external_field,time_field, &
         beff,beff1,beff2,convolution_energy,1,cell_index,emomM_macro,macro_nlistsize, &
         na,n1,n2,n3,measure_energy=.false.)
      call check(maxval(abs(beff1-direct)) < 2.0e-12_dblprec, &
         'field-only production convolution remains equal to DIRECT')
      call convolution_backend_get_stats(setup_seconds,pack_seconds,forward_seconds, &
         spectral_seconds,inverse_seconds,unpack_seconds,apply_seconds,apply_count)
      write(*,'(a,es12.4,a,5(es12.4,1x),a,es12.4)') 'CPU-HAM-04B setup=',setup_seconds, &
         ' pack/forward/multiply/inverse/unpack=',pack_seconds,forward_seconds, &
         spectral_seconds,inverse_seconds,unpack_seconds,' pair-field=',apply_seconds
      call check(apply_count==2_8 .and. setup_seconds >= 0.0_dblprec .and. &
         pack_seconds >= 0.0_dblprec .and. forward_seconds >= 0.0_dblprec .and. &
         spectral_seconds >= 0.0_dblprec .and. inverse_seconds >= 0.0_dblprec .and. &
         unpack_seconds >= 0.0_dblprec .and. apply_seconds >= 0.0_dblprec, &
         'production convolution reports persistent stage statistics')

      call cleanup_convolution_backend()
      do_convolution='N'
      deallocate(spin,mmom,external_field,time_field,beff,beff1,beff2,direct,emomM_macro, &
         cell_index,macro_nlistsize)
      call clear_production_fixture()
      call clear_reduced_stencil(stencil)
   end subroutine test_production_backend


   subroutine clear_production_fixture()
      call clear_reduced_stencil(ham%reduced_stencil)
      if (allocated(ham%aHam)) deallocate(ham%aHam)
      if (allocated(ham%nlistsize)) deallocate(ham%nlistsize)
      if (allocated(ham%nlist)) deallocate(ham%nlist)
      if (allocated(ham%ncoup)) deallocate(ham%ncoup)
   end subroutine clear_production_fixture


   subroutine test_ineligible_cases()
      integer, parameter :: natom=8
      integer :: ham_index(natom),nlistsize(natom),nlist(1,natom)
      real(dblprec) :: ncoup(1,1,1)
      logical :: ok
      character(len=256) :: diagnostic
      integer :: atom

      ham_index=1
      nlistsize=1
      ncoup=1.0_dblprec
      do atom=1,natom
         nlist(1,atom)=modulo(atom,natom)+1
      end do
      ok=cpu_convolution_eligible(natom,1,2,2,2,'P','P','0','Y',0,1, &
         ham_index,nlistsize,nlist,ncoup,diagnostic)
      call check(.not.ok .and. index(diagnostic,'periodic')>0, &
         'nonperiodic geometry is rejected')
      ok=cpu_convolution_eligible(natom,1,2,2,2,'P','P','P','N',0,1, &
         ham_index,nlistsize,nlist,ncoup,diagnostic)
      call check(.not.ok .and. index(diagnostic,'do_reduced')>0, &
         'non-reduced Hamiltonian is rejected')
      nlist(1,5)=1
      ok=cpu_convolution_eligible(natom,1,2,2,2,'P','P','P','Y',0,1, &
         ham_index,nlistsize,nlist,ncoup,diagnostic)
      call check(.not.ok .and. index(diagnostic,'translation-invariance')>0, &
         'translation-invariance violation is rejected')
   end subroutine test_ineligible_cases


   subroutine make_stencil(stencil,na,n1,n2,n3,nrecords)
      type(reduced_stencil_t), intent(out) :: stencil
      integer, intent(in) :: na,n1,n2,n3,nrecords
      integer :: a

      stencil%na=na
      stencil%n1=n1
      stencil%n2=n2
      stencil%n3=n3
      allocate(stencil%record_start(na+1),stencil%record(nrecords))
      do a=1,na+1
         stencil%record_start(a)=(a-1)*(nrecords/na)+1
      end do
   end subroutine make_stencil


   subroutine put_record(stencil,slot,output_basis,input_basis,delta,j)
      type(reduced_stencil_t), intent(inout) :: stencil
      integer, intent(in) :: slot,output_basis,input_basis,delta(3)
      real(dblprec), intent(in) :: j

      stencil%record(slot)%output_basis=output_basis
      stencil%record(slot)%input_basis=input_basis
      stencil%record(slot)%delta_cell=delta
      stencil%record(slot)%j=j
   end subroutine put_record


   subroutine apply_direct(stencil,spin,field)
      type(reduced_stencil_t), intent(in) :: stencil
      real(dblprec), intent(in) :: spin(:,:,:)
      real(dblprec), intent(out) :: field(:,:,:)
      integer :: ensemble,atom,b,slot,start,stop,source
      integer :: cell(3),wrapped(3)

      field=0.0_dblprec
      do ensemble=1,size(spin,3)
         do atom=1,size(spin,2)
            call atom_to_cell_basis(atom,stencil%na,stencil%n1,stencil%n2,stencil%n3,cell,b)
            start=stencil%record_start(b)
            stop=stencil%record_start(b+1)-1
            do slot=start,stop
               call wrap_reduced_cell(cell,stencil%record(slot)%delta_cell, &
                  stencil%n1,stencil%n2,stencil%n3,wrapped)
               source=cell_basis_to_atom(wrapped(1),wrapped(2),wrapped(3), &
                  stencil%record(slot)%input_basis,stencil%na,stencil%n1,stencil%n2,stencil%n3)
               field(:,atom,ensemble)=field(:,atom,ensemble)+stencil%record(slot)%j*spin(:,source,ensemble)
            end do
         end do
      end do
   end subroutine apply_direct


   subroutine check(condition,label)
      logical, intent(in) :: condition
      character(len=*), intent(in) :: label
      if (.not.condition) error stop 'FAILED: '//trim(label)
   end subroutine check

end program test_cpu_convolution
