!-------------------------------------------------------------------------------
! MODULE: CPUConvolution
!> @brief Persistent periodic lattice convolution backend for scalar J plus DMI.
!>
!> The object is provider-neutral at the Hamiltonian boundary and owns the
!> persistent mapping, spectral kernel, FFT plans and timestep work buffers
!> used by the explicit CPU-HAM-04B production path.  Its field result is
!> validated against the canonical DIRECT implementation by the test target.
!-------------------------------------------------------------------------------
module CPUConvolution

   use, intrinsic :: iso_c_binding, only : C_DOUBLE, C_DOUBLE_COMPLEX
   use Parameters, only : dblprec
   use ReducedStencil, only : reduced_stencil_t, build_reduced_stencil, &
      clear_reduced_stencil, reduced_stencil_eligible, atom_to_cell_basis, &
      cell_basis_to_atom
   use CPUFFTProvider, only : cpu_fft_plan_t, cpu_fft_plan_many_r2c, &
      cpu_fft_plan_many_c2r, cpu_fft_execute_r2c, cpu_fft_execute_c2r, &
      cpu_fft_plan_destroy

   implicit none
   private

   type, public :: cpu_convolution_t
      integer :: na=0
      integer :: n1=0
      integer :: n2=0
      integer :: n3=0
      integer :: ncells=0
      integer :: natom=0
      integer :: ensembles=0
      integer :: spectral_cells=0
      integer :: kernel_batches=0
      logical :: ready=.false.
      logical :: kernel_ready=.false.
      type(cpu_fft_plan_t) :: forward_plan
      type(cpu_fft_plan_t) :: backward_plan
      type(cpu_fft_plan_t) :: kernel_plan
      integer, allocatable :: atom_cell(:,:)
      integer, allocatable :: atom_basis(:)
      integer, allocatable :: atom_index(:,:)
      real(C_DOUBLE), allocatable :: real_work(:)
      complex(C_DOUBLE_COMPLEX), allocatable :: spin_spectral(:)
      complex(C_DOUBLE_COMPLEX), allocatable :: field_spectral(:)
      real(C_DOUBLE), allocatable :: kernel_real(:)
      complex(C_DOUBLE_COMPLEX), allocatable :: kernel_spectral(:)
      real(dblprec) :: pack_seconds=0.0_dblprec
      real(dblprec) :: forward_seconds=0.0_dblprec
      real(dblprec) :: spectral_seconds=0.0_dblprec
      real(dblprec) :: inverse_seconds=0.0_dblprec
      real(dblprec) :: unpack_seconds=0.0_dblprec
      real(dblprec) :: apply_seconds=0.0_dblprec
      integer(kind=8) :: apply_count=0_8
   end type cpu_convolution_t

   public :: cpu_convolution_eligible
   public :: cpu_convolution_init
   public :: cpu_convolution_build_kernel
   public :: cpu_convolution_apply
   public :: cpu_convolution_pair_energy
   public :: cpu_convolution_get_stats
   public :: cpu_convolution_clear

contains

   ! Validate production data and translation invariance before a convolution
   ! object is constructed.  Building a temporary stencil here is intentional:
   ! it reuses the accepted ReducedStencil validation rather than inventing a
   ! second connectivity oracle.
   logical function cpu_convolution_eligible(natom,na,n1,n2,n3,bc1,bc2,bc3, &
         do_reduced,do_ralloy,nham,ham_index,nlistsize,nlist,ncoup,diagnostic, &
         dmlistsize,dmlist,dm_vect)
      integer, intent(in) :: natom,na,n1,n2,n3,do_ralloy,nham
      character(len=1), intent(in) :: bc1,bc2,bc3,do_reduced
      integer, intent(in) :: ham_index(:),nlistsize(:),nlist(:,:)
      real(dblprec), intent(in) :: ncoup(:,:,:)
      character(len=*), intent(out), optional :: diagnostic
      integer, intent(in), optional :: dmlistsize(:), dmlist(:,:)
      real(dblprec), intent(in), optional :: dm_vect(:,:,:)
      type(reduced_stencil_t) :: stencil
      logical :: ok
      character(len=256) :: reason

      cpu_convolution_eligible=.false.
      reason=''
      if (.not.reduced_stencil_eligible(natom,na,n1,n2,n3,bc1,bc2,bc3, &
            do_reduced,do_ralloy,nham,ham_index,reason)) then
         call set_diagnostic(reason,diagnostic)
         return
      endif
      call build_reduced_stencil(natom,na,n1,n2,n3,bc1,bc2,bc3,do_reduced, &
         do_ralloy,nham,ham_index,nlistsize,nlist,ncoup,stencil,ok,reason, &
         dmlistsize,dmlist,dm_vect)
      call clear_reduced_stencil(stencil)
      cpu_convolution_eligible=ok
      call set_diagnostic(reason,diagnostic)
   end function cpu_convolution_eligible


   logical function cpu_convolution_init(convolution,stencil,ensembles,diagnostic)
      type(cpu_convolution_t), intent(inout) :: convolution
      type(reduced_stencil_t), intent(in) :: stencil
      integer, intent(in) :: ensembles
      character(len=*), intent(out), optional :: diagnostic
      integer :: atom,natom,batch_count,kernel_batches
      logical :: ok
      character(len=256) :: reason

      call cpu_convolution_clear(convolution)
      cpu_convolution_init=.false.
      reason=''
      if (.not.allocated(stencil%record_start) .or. .not.allocated(stencil%record)) then
         reason='reduced stencil is not allocated'
         call set_diagnostic(reason,diagnostic)
         return
      endif
      if (stencil%na < 1 .or. stencil%n1 < 1 .or. stencil%n2 < 1 .or. &
            stencil%n3 < 1 .or. ensembles < 1) then
         reason='convolution dimensions and ensemble count must be positive'
         call set_diagnostic(reason,diagnostic)
         return
      endif

      convolution%na=stencil%na
      convolution%n1=stencil%n1
      convolution%n2=stencil%n2
      convolution%n3=stencil%n3
      convolution%ncells=stencil%n1*stencil%n2*stencil%n3
      convolution%natom=stencil%na*convolution%ncells
      convolution%ensembles=ensembles
      convolution%spectral_cells=(stencil%n1/2+1)*stencil%n2*stencil%n3
      natom=convolution%natom
      batch_count=3*stencil%na*ensembles
      kernel_batches=stencil%na*stencil%na
      if (allocated(stencil%dmi_record_start)) kernel_batches=4*kernel_batches
      convolution%kernel_batches=kernel_batches

      allocate(convolution%atom_cell(3,natom),convolution%atom_basis(natom), &
         convolution%atom_index(stencil%na,convolution%ncells), &
         convolution%real_work(convolution%ncells*batch_count), &
         convolution%spin_spectral(convolution%spectral_cells*batch_count), &
         convolution%field_spectral(convolution%spectral_cells*batch_count), &
         convolution%kernel_real(convolution%ncells*kernel_batches), &
         convolution%kernel_spectral(convolution%spectral_cells*kernel_batches))

      do atom=1,natom
         call atom_to_cell_basis(atom,stencil%na,stencil%n1,stencil%n2,stencil%n3, &
            convolution%atom_cell(:,atom),convolution%atom_basis(atom))
         convolution%atom_index(convolution%atom_basis(atom), &
            1+convolution%atom_cell(1,atom)+stencil%n1*(convolution%atom_cell(2,atom)+ &
            stencil%n2*convolution%atom_cell(3,atom)))=atom
      end do

      ok=cpu_fft_plan_many_r2c(convolution%forward_plan,stencil%n1,stencil%n2, &
         stencil%n3,batch_count,convolution%real_work,convolution%spin_spectral)
      if (.not.ok) then
         reason='FFT provider could not create the forward plan'
         call cpu_convolution_clear(convolution)
         call set_diagnostic(reason,diagnostic)
         return
      endif
      ok=cpu_fft_plan_many_c2r(convolution%backward_plan,stencil%n1,stencil%n2, &
         stencil%n3,batch_count,convolution%field_spectral,convolution%real_work)
      if (.not.ok) then
         reason='FFT provider could not create the backward plan'
         call cpu_convolution_clear(convolution)
         call set_diagnostic(reason,diagnostic)
         return
      endif
      ok=cpu_fft_plan_many_r2c(convolution%kernel_plan,stencil%n1,stencil%n2, &
         stencil%n3,kernel_batches,convolution%kernel_real,convolution%kernel_spectral)
      if (.not.ok) then
         reason='FFT provider could not create the kernel plan'
         call cpu_convolution_clear(convolution)
         call set_diagnostic(reason,diagnostic)
         return
      endif

      convolution%ready=.true.
      reason='persistent CPU convolution state initialized'
      call set_diagnostic(reason,diagnostic)
      cpu_convolution_init=.true.
   end function cpu_convolution_init


   logical function cpu_convolution_build_kernel(convolution,stencil,diagnostic)
      type(cpu_convolution_t), intent(inout) :: convolution
      type(reduced_stencil_t), intent(in) :: stencil
      character(len=*), intent(out), optional :: diagnostic
      integer :: a,b,slot,start,stop,delta(3),x,y,z,cell,pair,component
      character(len=256) :: reason

      cpu_convolution_build_kernel=.false.
      reason=''
      if (.not.convolution%ready) then
         reason='convolution object is not initialized'
         call set_diagnostic(reason,diagnostic)
         return
      endif
      if (stencil%na /= convolution%na .or. stencil%n1 /= convolution%n1 .or. &
            stencil%n2 /= convolution%n2 .or. stencil%n3 /= convolution%n3 .or. &
            .not.allocated(stencil%record_start) .or. &
            .not.allocated(stencil%record)) then
         reason='stencil does not match convolution object'
         call set_diagnostic(reason,diagnostic)
         return
      endif

      convolution%kernel_real=0.0_C_DOUBLE
      do a=1,convolution%na
         start=stencil%record_start(a)
         stop=stencil%record_start(a+1)-1
         do slot=start,stop
            b=stencil%record(slot)%input_basis
            delta=stencil%record(slot)%delta_cell
            ! FFT(K)*FFT(M) produces sum_d K(d) M(R-d).  The DIRECT
            ! convention is M(R+delta), hence store the kernel at -delta.
            x=modulo(-delta(1),convolution%n1)
            y=modulo(-delta(2),convolution%n2)
            z=modulo(-delta(3),convolution%n3)
            cell=x+convolution%n1*(y+convolution%n2*z)
            pair=a+convolution%na*(b-1)
            convolution%kernel_real(cell+1+convolution%ncells*(pair-1))=&
               convolution%kernel_real(cell+1+convolution%ncells*(pair-1))+ &
               stencil%record(slot)%j
         end do
         if (allocated(stencil%dmi_record_start)) then
            start=stencil%dmi_record_start(a)
            stop=stencil%dmi_record_start(a+1)-1
            do slot=start,stop
               b=stencil%dmi_record(slot)%input_basis
               delta=stencil%dmi_record(slot)%delta_cell
               x=modulo(-delta(1),convolution%n1)
               y=modulo(-delta(2),convolution%n2)
               z=modulo(-delta(3),convolution%n3)
               cell=x+convolution%n1*(y+convolution%n2*z)
               pair=a+convolution%na*(b-1)
               do component=1,3
                  convolution%kernel_real(cell+1+convolution%ncells*(pair-1+ &
                     convolution%na*convolution%na*component))=&
                     convolution%kernel_real(cell+1+convolution%ncells*(pair-1+ &
                     convolution%na*convolution%na*component))+ &
                     stencil%dmi_record(slot)%d(component)
               end do
            end do
         endif
      end do
      call cpu_fft_execute_r2c(convolution%kernel_plan,convolution%kernel_real, &
         convolution%kernel_spectral)
      convolution%kernel_ready=.true.
      if (allocated(stencil%dmi_record_start)) then
         reason='J+D spectral kernels built with 4*N_A^2 basis kernels'
      else
         reason='scalar-J spectral kernel built with N_A^2 basis kernels'
      endif
      call set_diagnostic(reason,diagnostic)
      cpu_convolution_build_kernel=.true.
   end function cpu_convolution_build_kernel


   logical function cpu_convolution_apply(convolution,spin,field,diagnostic)
      type(cpu_convolution_t), intent(inout) :: convolution
      real(dblprec), intent(in) :: spin(:,:,:)
      real(dblprec), intent(out) :: field(:,:,:)
      character(len=*), intent(out), optional :: diagnostic
      integer :: ensemble,a,b,axis,cell,atom,q,batch_in,batch_out,pair,base_in,pair_count
      real(dblprec) :: stage_start
      real(dblprec) :: apply_start
      real(C_DOUBLE) :: scale
      complex(C_DOUBLE_COMPLEX) :: value,j_value,dx,dy,dz,mx,my,mz
      character(len=256) :: reason

      cpu_convolution_apply=.false.
      reason=''
      if (.not.convolution%ready .or. .not.convolution%kernel_ready) then
         reason='convolution object and spectral kernel must be ready'
         call set_diagnostic(reason,diagnostic)
         return
      endif
      if (size(spin,1) /= 3 .or. size(spin,2) /= convolution%natom .or. &
            size(spin,3) /= convolution%ensembles .or. &
            size(field,1) /= 3 .or. size(field,2) /= convolution%natom .or. &
            size(field,3) /= convolution%ensembles) then
         reason='spin and field dimensions do not match convolution object'
         call set_diagnostic(reason,diagnostic)
         return
      endif

      call cpu_time(apply_start)
      call cpu_time(stage_start)
      convolution%real_work=0.0_C_DOUBLE
      do ensemble=1,convolution%ensembles
         do cell=0,convolution%ncells-1
            do a=1,convolution%na
               atom=convolution%atom_index(a,cell+1)+(ensemble-1)*convolution%natom
               do axis=1,3
                  batch_in=axis+3*(a-1+convolution%na*(ensemble-1))
                  convolution%real_work(cell+1+convolution%ncells*(batch_in-1))=&
                     spin(axis,atom-(ensemble-1)*convolution%natom,ensemble)
               end do
            end do
         end do
      end do
      call accumulate_seconds(stage_start,convolution%pack_seconds)

      call cpu_time(stage_start)
      call cpu_fft_execute_r2c(convolution%forward_plan,convolution%real_work, &
         convolution%spin_spectral)
      call accumulate_seconds(stage_start,convolution%forward_seconds)

      call cpu_time(stage_start)
      convolution%field_spectral=cmplx(0.0_C_DOUBLE,0.0_C_DOUBLE,kind=C_DOUBLE)
      pair_count=convolution%na*convolution%na
      do ensemble=1,convolution%ensembles
         do a=1,convolution%na
            do axis=1,3
               batch_out=axis+3*(a-1+convolution%na*(ensemble-1))
               do q=0,convolution%spectral_cells-1
                  value=cmplx(0.0_C_DOUBLE,0.0_C_DOUBLE,kind=C_DOUBLE)
                  do b=1,convolution%na
                     base_in=3*(b-1+convolution%na*(ensemble-1))
                     batch_in=axis+base_in
                     pair=a+convolution%na*(b-1)
                     j_value=convolution%kernel_spectral(q+1+ &
                        convolution%spectral_cells*(pair-1))
                     mx=convolution%spin_spectral(q+1+convolution%spectral_cells*base_in)
                     my=convolution%spin_spectral(q+1+convolution%spectral_cells*(base_in+1))
                     mz=convolution%spin_spectral(q+1+convolution%spectral_cells*(base_in+2))
                     if (axis == 1) then
                        value=value+j_value*mx
                     elseif (axis == 2) then
                        value=value+j_value*my
                     else
                        value=value+j_value*mz
                     endif
                     if (convolution%kernel_batches > pair_count) then
                        dx=convolution%kernel_spectral(q+1+convolution%spectral_cells*(pair-1+pair_count))
                        dy=convolution%kernel_spectral(q+1+convolution%spectral_cells*(pair-1+2*pair_count))
                        dz=convolution%kernel_spectral(q+1+convolution%spectral_cells*(pair-1+3*pair_count))
                        if (axis == 1) then
                           value=value+dy*mz-dz*my
                        elseif (axis == 2) then
                           value=value+dz*mx-dx*mz
                        else
                           value=value+dx*my-dy*mx
                        endif
                     endif
                  end do
                  convolution%field_spectral(q+1+convolution%spectral_cells*(batch_out-1))=value
               end do
            end do
         end do
      end do
      call accumulate_seconds(stage_start,convolution%spectral_seconds)

      call cpu_time(stage_start)
      call cpu_fft_execute_c2r(convolution%backward_plan,convolution%field_spectral, &
         convolution%real_work)
      call accumulate_seconds(stage_start,convolution%inverse_seconds)

      call cpu_time(stage_start)
      scale=1.0_C_DOUBLE/real(convolution%ncells,C_DOUBLE)
      do ensemble=1,convolution%ensembles
         do cell=0,convolution%ncells-1
            do a=1,convolution%na
               atom=convolution%atom_index(a,cell+1)
               do axis=1,3
                  batch_out=axis+3*(a-1+convolution%na*(ensemble-1))
                  field(axis,atom,ensemble)=scale*convolution%real_work(cell+1+ &
                     convolution%ncells*(batch_out-1))
               end do
            end do
         end do
      end do
      call accumulate_seconds(stage_start,convolution%unpack_seconds)
      call accumulate_seconds(apply_start,convolution%apply_seconds)
      convolution%apply_count=convolution%apply_count+1_8
      reason='FFT convolution field applied'
      call set_diagnostic(reason,diagnostic)
      cpu_convolution_apply=.true.
   end function cpu_convolution_apply


   subroutine cpu_convolution_get_stats(convolution,pack_seconds,forward_seconds, &
         spectral_seconds,inverse_seconds,unpack_seconds,apply_seconds,apply_count)
      type(cpu_convolution_t), intent(in) :: convolution
      real(dblprec), intent(out) :: pack_seconds,forward_seconds,spectral_seconds
      real(dblprec), intent(out) :: inverse_seconds,unpack_seconds,apply_seconds
      integer(kind=8), intent(out) :: apply_count

      pack_seconds=convolution%pack_seconds
      forward_seconds=convolution%forward_seconds
      spectral_seconds=convolution%spectral_seconds
      inverse_seconds=convolution%inverse_seconds
      unpack_seconds=convolution%unpack_seconds
      apply_seconds=convolution%apply_seconds
      apply_count=convolution%apply_count
   end subroutine cpu_convolution_get_stats


   pure real(dblprec) function cpu_convolution_pair_energy(spin,field)
      real(dblprec), intent(in) :: spin(:,:,:),field(:,:,:)

      if (any(shape(spin) /= shape(field))) then
         cpu_convolution_pair_energy=0.0_dblprec
      else
         cpu_convolution_pair_energy=-0.5_dblprec*sum(spin*field)
      endif
   end function cpu_convolution_pair_energy


   subroutine cpu_convolution_clear(convolution)
      type(cpu_convolution_t), intent(inout) :: convolution

      call cpu_fft_plan_destroy(convolution%forward_plan)
      call cpu_fft_plan_destroy(convolution%backward_plan)
      call cpu_fft_plan_destroy(convolution%kernel_plan)
      if (allocated(convolution%atom_cell)) deallocate(convolution%atom_cell)
      if (allocated(convolution%atom_basis)) deallocate(convolution%atom_basis)
      if (allocated(convolution%atom_index)) deallocate(convolution%atom_index)
      if (allocated(convolution%real_work)) deallocate(convolution%real_work)
      if (allocated(convolution%spin_spectral)) deallocate(convolution%spin_spectral)
      if (allocated(convolution%field_spectral)) deallocate(convolution%field_spectral)
      if (allocated(convolution%kernel_real)) deallocate(convolution%kernel_real)
      if (allocated(convolution%kernel_spectral)) deallocate(convolution%kernel_spectral)
      convolution%na=0
      convolution%n1=0
      convolution%n2=0
      convolution%n3=0
      convolution%ncells=0
      convolution%natom=0
      convolution%ensembles=0
      convolution%spectral_cells=0
      convolution%kernel_batches=0
      convolution%ready=.false.
      convolution%kernel_ready=.false.
      convolution%pack_seconds=0.0_dblprec
      convolution%forward_seconds=0.0_dblprec
      convolution%spectral_seconds=0.0_dblprec
      convolution%inverse_seconds=0.0_dblprec
      convolution%unpack_seconds=0.0_dblprec
      convolution%apply_seconds=0.0_dblprec
      convolution%apply_count=0_8
   end subroutine cpu_convolution_clear


   subroutine accumulate_seconds(time_start,total)
      real(dblprec), intent(in) :: time_start
      real(dblprec), intent(inout) :: total
      real(dblprec) :: time_stop

      call cpu_time(time_stop)
      total=total+max(0.0_dblprec,time_stop-time_start)
   end subroutine accumulate_seconds


   subroutine set_diagnostic(reason,diagnostic)
      character(len=*), intent(in) :: reason
      character(len=*), intent(out), optional :: diagnostic

      if (present(diagnostic)) diagnostic=reason
   end subroutine set_diagnostic

end module CPUConvolution
