! CPU-HAM-04A: scalar-J periodic FFT convolution reference and oracle.
program test_cpu_convolution

   use Parameters, only : dblprec
   use ReducedStencil
   use CPUConvolution

   implicit none

   call test_delta_and_uniform()
   call test_multibasis_and_ensembles()
   call test_long_range_nd_fixture()
   call test_ineligible_cases()
   write(*,'(a)') 'CPU-HAM-04A CPU convolution tests passed'

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
