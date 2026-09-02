! CPU-HAM-09: repeated stage timing for the persistent CPU convolution path.
program benchmark_cpu_convolution

   use Parameters, only : dblprec
   use CPUFFTProvider, only : cpu_fft_provider_name, cpu_fft_provider_threads
   use CPUConvolution
   use ReducedStencil

   implicit none

   write(*,'(a,a,a,i0)') 'CPU-HAM-09 provider=',trim(cpu_fft_provider_name()), &
      ' provider_threads=',cpu_fft_provider_threads()
   call run_case('Nd scalar-J long-range',4,25,25,25,1338,.false.,1,5)
   call run_case('Fe scalar-J medium-range',2,20,20,20,96,.false.,2,5)
   call run_case('short-range scalar-J control',1,32,32,32,6,.false.,1,5)
   call run_case('2D J+D',1,32,32,1,4,.true.,2,5)
   call run_case('3D J+D multi-basis',2,16,16,16,6,.true.,2,5)

contains

   subroutine run_case(label,na,n1,n2,n3,z,with_dmi,ensembles,repetitions)
      character(len=*), intent(in) :: label
      integer, intent(in) :: na,n1,n2,n3,z,ensembles,repetitions
      logical, intent(in) :: with_dmi
      type(reduced_stencil_t) :: stencil
      type(cpu_convolution_t) :: convolution
      real(dblprec), allocatable :: spin(:,:,:),field(:,:,:),direct(:,:,:)
      real(dblprec) :: setup_start,setup_seconds,apply_seconds
      real(dblprec) :: previous(6),current(6),delta(6),max_error
      real(dblprec) :: stage_start,stage_stop
      character(len=256) :: diagnostic
      logical :: ok
      integer :: natom,rep,axis,atom,ensemble
      integer(kind=8) :: apply_count

      natom=na*n1*n2*n3
      call make_stencil(stencil,na,n1,n2,n3,z,with_dmi)
      allocate(spin(3,natom,ensembles),field(3,natom,ensembles),direct(3,natom,ensembles))
      do ensemble=1,ensembles
         do atom=1,natom
            do axis=1,3
               spin(axis,atom,ensemble)=sin(0.017_dblprec*real(axis+atom+11*ensemble,dblprec))+ &
                  0.13_dblprec*cos(0.031_dblprec*real(2*axis+atom+ensemble,dblprec))
            end do
         end do
      end do

      setup_start=wall_seconds()
      ok=cpu_convolution_init(convolution,stencil,ensembles,diagnostic)
      if (ok) ok=cpu_convolution_build_kernel(convolution,stencil,diagnostic)
      setup_seconds=wall_seconds()-setup_start
      if (.not.ok) error stop 'CPU-HAM-09 benchmark setup failed: '//trim(diagnostic)

      ! Warm up FFTW and the OpenMP runtime before collecting samples.
      ok=cpu_convolution_apply(convolution,spin,field,diagnostic)
      if (.not.ok) error stop 'CPU-HAM-09 benchmark warmup failed: '//trim(diagnostic)
      call apply_direct(stencil,spin,direct)
      max_error=maxval(abs(field-direct))
      ! Re-warm after the serial DIRECT oracle so the collected samples start
      ! with the same FFT/cache state as a production steady-state apply.
      ok=cpu_convolution_apply(convolution,spin,field,diagnostic)
      if (.not.ok) error stop 'CPU-HAM-09 benchmark second warmup failed: '//trim(diagnostic)
      call cpu_convolution_get_stats(convolution,previous(1),previous(2),previous(3), &
         previous(4),previous(5),previous(6),apply_count)

      write(*,'(a)') ''
      write(*,'(a,a,a,i0,a,i0,a,i0,a,i0,a,i0,a,l1,a,es12.4,a,es12.4)') &
         'case=',trim(label),' natom=',natom,' grid=',n1,'x',n2,'x',n3, &
         ' basis=',na,' dmi=',with_dmi,' setup_s=',setup_seconds, &
         ' parity_max=',max_error
      write(*,'(a)') 'sample wall_s pack_s forward_s spectral_s inverse_s unpack_s apply_s'
      do rep=1,repetitions
         stage_start=wall_seconds()
         ok=cpu_convolution_apply(convolution,spin,field,diagnostic)
         stage_stop=wall_seconds()
         if (.not.ok) error stop 'CPU-HAM-09 benchmark apply failed: '//trim(diagnostic)
         apply_seconds=stage_stop-stage_start
         call cpu_convolution_get_stats(convolution,current(1),current(2),current(3), &
            current(4),current(5),current(6),apply_count)
         delta=current-previous
         write(*,'(i0,7(1x,es12.4))') rep,apply_seconds,delta
         previous=current
      end do

      call cpu_convolution_clear(convolution)
      call clear_reduced_stencil(stencil)
      deallocate(spin,field,direct)
   end subroutine run_case


   subroutine make_stencil(stencil,na,n1,n2,n3,z,with_dmi)
      type(reduced_stencil_t), intent(out) :: stencil
      integer, intent(in) :: na,n1,n2,n3,z
      logical, intent(in) :: with_dmi
      integer :: a,slot,index,source_basis
      integer :: delta(3)

      stencil%na=na
      stencil%n1=n1
      stencil%n2=n2
      stencil%n3=n3
      allocate(stencil%record_start(na+1),stencil%record(na*z))
      if (with_dmi) allocate(stencil%dmi_record_start(na+1),stencil%dmi_record(na*z))
      do a=1,na
         stencil%record_start(a)=(a-1)*z+1
         if (with_dmi) stencil%dmi_record_start(a)=(a-1)*z+1
         do slot=1,z
            index=(a-1)*z+slot
            source_basis=modulo(a+slot-2,na)+1
            call benchmark_delta(slot,n1,n2,n3,delta)
            stencil%record(index)%output_basis=a
            stencil%record(index)%input_basis=source_basis
            stencil%record(index)%delta_cell=delta
            stencil%record(index)%j=0.01_dblprec*real(2*a-slot,dblprec)
            if (with_dmi) then
               stencil%dmi_record(index)%output_basis=a
               stencil%dmi_record(index)%input_basis=source_basis
               stencil%dmi_record(index)%delta_cell=delta
               stencil%dmi_record(index)%d=(/0.003_dblprec*real(slot,dblprec), &
                  -0.002_dblprec*real(a+slot,dblprec),0.001_dblprec*real(a,dblprec)/)
            endif
         end do
      end do
      stencil%record_start(na+1)=na*z+1
      if (with_dmi) stencil%dmi_record_start(na+1)=na*z+1
   end subroutine make_stencil


   subroutine benchmark_delta(slot,n1,n2,n3,delta)
      integer, intent(in) :: slot,n1,n2,n3
      integer, intent(out) :: delta(3)

      delta=0
      select case(modulo(slot-1,6))
      case(0)
         delta(1)=1
      case(1)
         delta(1)=-1
      case(2)
         delta(2)=1
      case(3)
         delta(2)=-1
      case(4)
         delta(3)=1
      case(5)
         delta(3)=-1
      end select
      if (n1 == 1) delta(1)=0
      if (n2 == 1) delta(2)=0
      if (n3 == 1) delta(3)=0
   end subroutine benchmark_delta


   subroutine apply_direct(stencil,spin,field)
      type(reduced_stencil_t), intent(in) :: stencil
      real(dblprec), intent(in) :: spin(:,:,:)
      real(dblprec), intent(out) :: field(:,:,:)
      integer :: ensemble,atom,b,slot,start,stop,source
      integer :: cell(3),wrapped(3)
      real(dblprec) :: d(3),m(3)

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
            if (allocated(stencil%dmi_record_start)) then
               start=stencil%dmi_record_start(b)
               stop=stencil%dmi_record_start(b+1)-1
               do slot=start,stop
                  call wrap_reduced_cell(cell,stencil%dmi_record(slot)%delta_cell, &
                     stencil%n1,stencil%n2,stencil%n3,wrapped)
                  source=cell_basis_to_atom(wrapped(1),wrapped(2),wrapped(3), &
                     stencil%dmi_record(slot)%input_basis,stencil%na,stencil%n1,stencil%n2,stencil%n3)
                  d=stencil%dmi_record(slot)%d
                  m=spin(:,source,ensemble)
                  field(1,atom,ensemble)=field(1,atom,ensemble)+d(2)*m(3)-d(3)*m(2)
                  field(2,atom,ensemble)=field(2,atom,ensemble)+d(3)*m(1)-d(1)*m(3)
                  field(3,atom,ensemble)=field(3,atom,ensemble)+d(1)*m(2)-d(2)*m(1)
               end do
            endif
         end do
      end do
   end subroutine apply_direct


   real(dblprec) function wall_seconds()
      integer(kind=8) :: count,rate

      call system_clock(count,rate)
      wall_seconds=real(count,dblprec)/real(rate,dblprec)
   end function wall_seconds

end program benchmark_cpu_convolution
