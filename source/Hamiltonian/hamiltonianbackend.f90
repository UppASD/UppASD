!-------------------------------------------------------------------------------
! MODULE: HamiltonianBackend
!> @brief
!> CPU Hamiltonian backend selection, validation, setup, application and cleanup.
!>
!> This module owns all persistent backend state.  HamiltonianActions only
!> consumes the prepared pair field and assembles the physical field kernels.
!-------------------------------------------------------------------------------
module HamiltonianBackend

   use Parameters
   use HamiltonianData
   use InputData, only : ham_inp, cpu_ham_backend, do_sparse, do_convolution
#ifdef USE_FFTW
   use CPUConvolution, only : cpu_convolution_t, cpu_convolution_eligible, &
      cpu_convolution_init, cpu_convolution_build_kernel, cpu_convolution_apply, &
      cpu_convolution_get_stats, cpu_convolution_clear
   use CPUFFTProvider, only : cpu_fft_provider_name, cpu_fft_provider_threads
#endif

   implicit none
   private

   public :: set_reduced_direct_testing, reduced_direct_testing_enabled
   public :: resolve_cpu_ham_backend, resolve_requested_cpu_ham_backend
   public :: setup_cpu_hamiltonian_backend, cleanup_cpu_hamiltonian_backend
   public :: setup_sparse_backend, cleanup_sparse_backend
   public :: sparse_backend_can_apply, apply_sparse_exchange
   public :: sparse_backend_get_stats
   public :: setup_convolution_backend, cleanup_convolution_backend
   public :: convolution_backend_can_apply, apply_convolution_exchange
   public :: convolution_backend_get_stats
   public :: prepare_cpu_hamiltonian_backend, get_backend_pair_field

   character(len=16) :: cpu_ham_backend_resolved = 'direct'
   character(len=256) :: cpu_ham_backend_reason = 'not configured'
   logical :: cpu_ham_backend_initialized = .false.
   logical :: reduced_direct_testing = .false.
   logical :: partial_backend_fallback_reported = .false.

   ! Persistent CPU sparse-backend state. The matrix is a directed CSR
   ! representation of the canonical neighbour-list operator; physical atom
   ! indices and canonical Hamiltonian arrays are never reordered.
   integer, dimension(:), allocatable :: sparse_row_ptr
   integer, dimension(:), allocatable :: sparse_columns
   real(dblprec), dimension(:), allocatable :: sparse_values
   real(dblprec), dimension(:,:,:), allocatable :: sparse_rhs
   real(dblprec), dimension(:,:,:), allocatable :: sparse_field
   logical :: sparse_backend_ready = .false.
   integer :: sparse_natom = 0
   integer :: sparse_mensemble = 0
   integer(kind=8) :: sparse_nnz = 0_8
   real(dblprec) :: sparse_setup_seconds = 0.0_dblprec
   real(dblprec) :: sparse_pack_seconds = 0.0_dblprec
   real(dblprec) :: sparse_apply_seconds = 0.0_dblprec

#ifdef USE_FFTW
   type(cpu_convolution_t) :: convolution_backend
   real(dblprec), dimension(:,:,:), allocatable :: convolution_field
   logical :: convolution_backend_ready = .false.
   real(dblprec) :: convolution_setup_seconds = 0.0_dblprec
#endif

contains

   subroutine set_reduced_direct_testing(enabled)
      logical, intent(in) :: enabled

      reduced_direct_testing=enabled
   end subroutine set_reduced_direct_testing

   logical function reduced_direct_testing_enabled()
      reduced_direct_testing_enabled=reduced_direct_testing
   end function reduced_direct_testing_enabled

   !----------------------------------------------------------------------------
   !> Normalize a backend name without making input spelling significant.
   !----------------------------------------------------------------------------
   pure function lowercase_backend_name(value) result(lower)
      character(len=*), intent(in) :: value
      character(len=len(value)) :: lower
      integer :: i,code

      do i=1,len(value)
         code=iachar(value(i:i))
         if (code >= iachar('A') .and. code <= iachar('Z')) then
            lower(i:i)=achar(code+iachar('a')-iachar('A'))
         else
            lower(i:i)=value(i:i)
         endif
      enddo
   end function lowercase_backend_name


   !----------------------------------------------------------------------------
   !> Resolve the public CPU pair-backend vocabulary.
   !>
   !> AUTO is deliberately rejected: CPU-HAM-05 found no portable crossover
   !> rule, and a machine-dependent heuristic must not be inferred here.
   !----------------------------------------------------------------------------
   subroutine resolve_cpu_ham_backend(requested,resolved,ok,diagnostic)
      character(len=*), intent(in) :: requested
      character(len=*), intent(out) :: resolved
      logical, intent(out) :: ok
      character(len=*), intent(out) :: diagnostic
      character(len=16) :: normalized

      normalized=trim(lowercase_backend_name(requested))
      resolved=''
      diagnostic=''
      ok=.true.
      select case(normalized)
      case('direct')
         resolved='direct'
      case('sparse')
         resolved='sparse'
      case('convolution')
         resolved='convolution'
      case('auto')
         ok=.false.
         diagnostic='AUTO is disabled: CPU-HAM-05 did not establish a portable crossover rule'
      case default
         ok=.false.
         diagnostic='unsupported CPU Hamiltonian backend "'//trim(requested)// &
            '"; expected direct, sparse or convolution'
      end select
   end subroutine resolve_cpu_ham_backend


   logical function legacy_backend_flag_enabled(flag)
      character(len=*), intent(in) :: flag

      select case(trim(lowercase_backend_name(flag)))
      case('y','t','1')
         legacy_backend_flag_enabled=.true.
      case default
         legacy_backend_flag_enabled=.false.
      end select
   end function legacy_backend_flag_enabled


   !----------------------------------------------------------------------------
   !> Resolve the canonical request plus the pre-CPU-HAM-07 compatibility flags.
   !----------------------------------------------------------------------------
   subroutine resolve_requested_cpu_ham_backend(resolved,ok,diagnostic)
      character(len=*), intent(out) :: resolved
      logical, intent(out) :: ok
      character(len=*), intent(out) :: diagnostic
      character(len=16) :: requested
      logical :: legacy_sparse,legacy_convolution

      call resolve_cpu_ham_backend(cpu_ham_backend,resolved,ok,diagnostic)
      if (.not.ok) return

      legacy_sparse=legacy_backend_flag_enabled(do_sparse)
      legacy_convolution=legacy_backend_flag_enabled(do_convolution)
      if (legacy_sparse .and. legacy_convolution) then
         ok=.false.
         diagnostic='do_sparse and do_convolution cannot select two CPU Hamiltonian backends'
         return
      endif

      requested=resolved
      if (trim(resolved) == 'direct') then
         if (legacy_sparse) then
            resolved='sparse'
            diagnostic='legacy do_sparse compatibility alias'
         elseif (legacy_convolution) then
            resolved='convolution'
            diagnostic='legacy do_convolution compatibility alias'
         else
            diagnostic='explicit DIRECT request; safe general path'
         endif
      elseif (trim(resolved) == 'sparse' .and. legacy_convolution) then
         ok=.false.
         diagnostic='cpu_ham_backend=sparse conflicts with do_convolution'
      elseif (trim(resolved) == 'convolution' .and. legacy_sparse) then
         ok=.false.
         diagnostic='cpu_ham_backend=convolution conflicts with do_sparse'
      else
         diagnostic='explicit '//trim(requested)//' request'
      endif
   end subroutine resolve_requested_cpu_ham_backend


   subroutine report_cpu_ham_backend(requested,resolved,eligible,reason)
      character(len=*), intent(in) :: requested,resolved,reason
      logical, intent(in) :: eligible

      write(*,'(2x,a,a,a,a,a,l1,a,a)') 'CPU Hamiltonian backend: requested=',trim(requested), &
         ' resolved=',trim(resolved),' eligible=',eligible,' reason=',trim(reason)
   end subroutine report_cpu_ham_backend


   subroutine reject_cpu_ham_backend(reason,resolved)
      character(len=*), intent(in) :: reason
      character(len=*), intent(in), optional :: resolved
      character(len=16) :: resolved_name

      resolved_name='unresolved'
      if (present(resolved)) resolved_name=trim(resolved)
      write(*,'(/,1x,a)') 'ERROR: CPU Hamiltonian backend request rejected.'
      write(*,'(3x,a,a)') 'requested = ',trim(cpu_ham_backend)
      write(*,'(3x,a,a)') 'resolved  = ',trim(resolved_name)
      write(*,'(3x,a,a)') 'reason    = ',trim(reason)
      ! This is a user-facing configuration failure.  Keep the non-zero exit
      ! status, but avoid turning it into a compiler/runtime backtrace that
      ! obscures the actionable diagnostic above.
      stop 1
   end subroutine reject_cpu_ham_backend


   !----------------------------------------------------------------------------
   !> Construct the explicitly requested production CPU pair backend.
   !>
   !> DIRECT performs no persistent setup.  SPARSE and CONVOLUTION must prove
   !> eligibility and become ready; an explicit ineligible request is fatal so
   !> production never silently changes backend or physics.
   !----------------------------------------------------------------------------
   subroutine setup_cpu_hamiltonian_backend(Natom,Mensemble,do_ralloy,do_lsf,NA,N1,N2,N3, &
         BC1,BC2,BC3,do_reduced,nHam)
      implicit none

      integer, intent(in) :: Natom,Mensemble,do_ralloy,NA,N1,N2,N3,nHam
      character(len=1), intent(in) :: do_lsf,BC1,BC2,BC3,do_reduced
      character(len=16) :: resolved
      logical :: ok
      character(len=256) :: diagnostic

      call cleanup_cpu_hamiltonian_backend()
      call resolve_requested_cpu_ham_backend(resolved,ok,diagnostic)
      if (.not.ok) call reject_cpu_ham_backend(trim(diagnostic),resolved)

      select case(trim(resolved))
      case('direct')
         cpu_ham_backend_resolved='direct'
         cpu_ham_backend_reason='DIRECT is the safe general path; canonical HamiltonianActions terms remain active'
         cpu_ham_backend_initialized=.true.
         call report_cpu_ham_backend(cpu_ham_backend,resolved,.true.,cpu_ham_backend_reason)
      case('sparse')
         call setup_sparse_backend(Natom,Mensemble,do_ralloy,do_lsf)
         if (.not.sparse_backend_ready) then
            if (len_trim(cpu_ham_backend_reason) == 0) cpu_ham_backend_reason= &
               'scalar-J sparse backend is ineligible for this Hamiltonian'
            call reject_cpu_ham_backend(trim(cpu_ham_backend_reason),'sparse')
         endif
         cpu_ham_backend_resolved='sparse'
         cpu_ham_backend_reason='eligible scalar-J; persistent CSR backend initialized'
         cpu_ham_backend_initialized=.true.
         call report_cpu_ham_backend(cpu_ham_backend,resolved,.true.,cpu_ham_backend_reason)
      case('convolution')
         call setup_convolution_backend(Natom,Mensemble,do_ralloy,do_lsf,NA,N1,N2,N3, &
            BC1,BC2,BC3,do_reduced,nHam)
         if (.not.convolution_backend_can_apply(Natom,Mensemble,1,Natom)) then
            if (len_trim(cpu_ham_backend_reason) == 0) cpu_ham_backend_reason= &
               'periodic reduced scalar-J/DMI convolution backend is unavailable or ineligible'
            call reject_cpu_ham_backend(trim(cpu_ham_backend_reason),'convolution')
         endif
         cpu_ham_backend_resolved='convolution'
         cpu_ham_backend_reason='eligible periodic reduced scalar-J/DMI; persistent convolution backend initialized'
         cpu_ham_backend_initialized=.true.
         call report_cpu_ham_backend(cpu_ham_backend,resolved,.true.,cpu_ham_backend_reason)
      case default
         call reject_cpu_ham_backend('internal backend resolution failure',resolved)
      end select
   end subroutine setup_cpu_hamiltonian_backend


   subroutine cleanup_cpu_hamiltonian_backend()
      implicit none

      call cleanup_sparse_backend()
      call cleanup_convolution_backend()
      cpu_ham_backend_resolved='direct'
      cpu_ham_backend_reason='not configured'
      cpu_ham_backend_initialized=.false.
      reduced_direct_testing=.false.
      partial_backend_fallback_reported=.false.
   end subroutine cleanup_cpu_hamiltonian_backend

   !----------------------------------------------------------------------------
   !> Prepare the selected backend for one effective-field request.  The
   !> backend owns range eligibility, partial-range fallback reporting and
   !> application of any persistent work buffers.
   !----------------------------------------------------------------------------
   subroutine prepare_cpu_hamiltonian_backend(Natom,Mensemble,start_atom,stop_atom,emomM, &
         sparse_active,convolution_active)
      implicit none

      integer, intent(in) :: Natom,Mensemble,start_atom,stop_atom
      real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: emomM
      logical, intent(out) :: sparse_active,convolution_active

      sparse_active=sparse_backend_can_apply(Natom,Mensemble,start_atom,stop_atom)
      if (sparse_active) call apply_sparse_exchange(Mensemble,emomM)
      convolution_active=convolution_backend_can_apply(Natom,Mensemble,start_atom,stop_atom)
      if (convolution_active) call apply_convolution_exchange(Mensemble,emomM)

      if ((start_atom /= 1 .or. stop_atom /= Natom) .and. cpu_ham_backend_initialized .and. &
            (trim(cpu_ham_backend_resolved) == 'sparse' .or. trim(cpu_ham_backend_resolved) == 'convolution')) then
         if (.not.partial_backend_fallback_reported) then
            write(*,'(2x,a,a,a,i0,a,i0,a)') 'CPU Hamiltonian backend: partial range ', &
               trim(cpu_ham_backend_resolved),' request [',start_atom,',',stop_atom, &
               '] uses intentional DIRECT fallback'
            partial_backend_fallback_reported=.true.
         endif
      endif
   end subroutine prepare_cpu_hamiltonian_backend

   !----------------------------------------------------------------------------
   !> Return the pair field produced by the selected persistent backend.
   !----------------------------------------------------------------------------
   subroutine get_backend_pair_field(i,k,use_sparse_pair,use_convolution_pair,pair_field, &
         pair_includes_dmi)
      implicit none

      integer, intent(in) :: i,k
      logical, intent(in) :: use_sparse_pair,use_convolution_pair
      real(dblprec), dimension(3), intent(out) :: pair_field
      logical, intent(out) :: pair_includes_dmi

      pair_field=0.0_dblprec
      pair_includes_dmi=.false.
      if (use_convolution_pair) then
#ifdef USE_FFTW
         pair_field=convolution_field(:,i,k)
         pair_includes_dmi=(ham_inp%do_dm == 1)
#else
         error stop 'CPU convolution backend is unavailable'
#endif
      elseif (use_sparse_pair) then
         pair_field=sparse_field(i,:,k)
      endif
   end subroutine get_backend_pair_field

   !----------------------------------------------------------------------------
   !> Construct the persistent portable scalar-J CSR backend.
   !>
   !> This is an explicit scalar-J setup helper.  Production selection is
   !> owned by setup_cpu_hamiltonian_backend; legacy do_sparse callers remain
   !> supported for the focused backend tests and compatibility path.
   !----------------------------------------------------------------------------
   subroutine setup_sparse_backend(Natom,Mensemble,do_ralloy,do_lsf)
      use omp_lib, only : omp_get_wtime
      implicit none

      integer, intent(in) :: Natom, Mensemble, do_ralloy
      character(len=1), intent(in) :: do_lsf
      integer :: i, j, ih, nnz, pos, istat
      real(dblprec) :: t_start
      character(len=16) :: requested_backend
      character(len=256) :: backend_diagnostic
      logical :: backend_ok

      call cleanup_sparse_backend()
      cpu_ham_backend_reason=''
      call resolve_requested_cpu_ham_backend(requested_backend,backend_ok,backend_diagnostic)
      if (.not.backend_ok) then
         cpu_ham_backend_reason=trim(backend_diagnostic)
         write(*,'(2x,a,a)') 'Scalar-J sparse backend declined: ',trim(backend_diagnostic)
         return
      endif
      if (trim(requested_backend) /= 'sparse') return

      if (ham_inp%do_jtensor == 1 .or. ham_inp%do_dm == 1 .or. &
            ham_inp%exc_inter /= 'N' .or. do_ralloy /= 0 .or. do_lsf == 'Y') then
         cpu_ham_backend_reason='sparse supports only static scalar-J Hamiltonians without DMI, rescaling, disorder or LSF'
         write(*,'(2x,a)') 'Scalar-J sparse backend declined: unsupported Hamiltonian variant'
         return
      endif
      if (.not. allocated(ham%aHam) .or. .not. allocated(ham%nlistsize) .or. &
            .not. allocated(ham%nlist) .or. .not. allocated(ham%ncoup)) then
         cpu_ham_backend_reason='sparse exchange data unavailable'
         write(*,'(2x,a)') 'Scalar-J sparse backend declined: exchange data unavailable'
         return
      endif
      if (size(ham%aHam) /= Natom .or. size(ham%nlist,2) < Natom) then
         cpu_ham_backend_reason='sparse exchange dimensions do not cover all physical atoms'
         write(*,'(2x,a)') 'Scalar-J sparse backend declined: physical atom dimensions unavailable'
         return
      endif

      if (any(ham%aHam(1:Natom) < 1) .or. any(ham%aHam(1:Natom) > size(ham%nlistsize))) then
         cpu_ham_backend_reason='sparse Hamiltonian map leaves the available scalar-J rows'
         write(*,'(2x,a)') 'Scalar-J sparse backend declined: invalid reduced/non-reduced row map'
         return
      endif
      if (size(ham%ncoup,1) < maxval(ham%nlistsize(ham%aHam(1:Natom))) .or. &
            size(ham%ncoup,2) < maxval(ham%aHam(1:Natom)) .or. size(ham%ncoup,3) < 1) then
         cpu_ham_backend_reason='sparse coupling dimensions do not cover the selected Hamiltonian rows'
         write(*,'(2x,a)') 'Scalar-J sparse backend declined: coupling dimensions unavailable'
         return
      endif

      t_start = omp_get_wtime()
      nnz = 0
      do i=1,Natom
         ih=ham%aHam(i)
         if (ih < 1 .or. ih > size(ham%nlistsize)) then
            cpu_ham_backend_reason='sparse Hamiltonian map contains an invalid exchange index'
            write(*,'(2x,a)') 'Scalar-J sparse backend declined: invalid Hamiltonian map'
            return
         endif
         if (ham%nlistsize(ih) < 0 .or. ham%nlistsize(ih) > size(ham%nlist,1) .or. &
               ham%nlistsize(ih) > size(ham%ncoup,1)) then
            cpu_ham_backend_reason='sparse Hamiltonian row has an invalid neighbour count'
            write(*,'(2x,a)') 'Scalar-J sparse backend declined: invalid neighbour count'
            return
         endif
         nnz=nnz+ham%nlistsize(ih)
      end do

      allocate(sparse_row_ptr(Natom+1),stat=istat)
      if (istat /= 0) error stop 'Unable to allocate sparse CSR row pointers'
      allocate(sparse_columns(nnz),sparse_values(nnz),stat=istat)
      if (istat /= 0) error stop 'Unable to allocate sparse CSR entries'
      allocate(sparse_rhs(Natom,3,Mensemble),sparse_field(Natom,3,Mensemble),stat=istat)
      if (istat /= 0) error stop 'Unable to allocate sparse RHS work buffers'

      sparse_row_ptr(1)=0
      pos=0
      do i=1,Natom
         ih=ham%aHam(i)
         do j=1,ham%nlistsize(ih)
            if (ham%nlist(j,i) < 1 .or. ham%nlist(j,i) > Natom) then
               cpu_ham_backend_reason='sparse Hamiltonian map contains an invalid physical neighbour index'
               write(*,'(2x,a)') 'Scalar-J sparse backend declined: invalid physical neighbour index'
               return
            endif
            pos=pos+1
            sparse_columns(pos)=ham%nlist(j,i)
            sparse_values(pos)=ham%ncoup(j,ih,1)
         end do
         sparse_row_ptr(i+1)=pos
      end do

      sparse_natom=Natom
      sparse_mensemble=Mensemble
      sparse_nnz=int(nnz,kind=8)
      sparse_backend_ready=.true.
      sparse_setup_seconds=omp_get_wtime()-t_start
      cpu_ham_backend_resolved='sparse'
      cpu_ham_backend_initialized=.true.
      cpu_ham_backend_reason='eligible scalar-J; persistent CSR backend initialized'
      write(*,'(2x,a,i0,a,i0,a)') 'Persistent scalar-J sparse backend ready: ', &
         Natom, ' atoms, ', nnz, ' directed entries'
   end subroutine setup_sparse_backend

   !----------------------------------------------------------------------------
   !> Release all persistent sparse state. Safe for repeated setup/cleanup.
   !----------------------------------------------------------------------------
   subroutine cleanup_sparse_backend()
      implicit none

      if (allocated(sparse_row_ptr)) deallocate(sparse_row_ptr)
      if (allocated(sparse_columns)) deallocate(sparse_columns)
      if (allocated(sparse_values)) deallocate(sparse_values)
      if (allocated(sparse_rhs)) deallocate(sparse_rhs)
      if (allocated(sparse_field)) deallocate(sparse_field)
      sparse_backend_ready=.false.
      sparse_natom=0
      sparse_mensemble=0
      sparse_nnz=0_8
      sparse_setup_seconds=0.0_dblprec
      sparse_pack_seconds=0.0_dblprec
      sparse_apply_seconds=0.0_dblprec
      if (trim(cpu_ham_backend_resolved) == 'sparse') cpu_ham_backend_resolved='direct'
      cpu_ham_backend_initialized=.false.
   end subroutine cleanup_sparse_backend

   !----------------------------------------------------------------------------
   !> Return whether the explicit sparse backend can serve this field request.
   !----------------------------------------------------------------------------
   logical function sparse_backend_can_apply(Natom,Mensemble,start_atom,stop_atom)
      implicit none
      integer, intent(in) :: Natom,Mensemble,start_atom,stop_atom

      sparse_backend_can_apply = sparse_backend_ready .and. trim(cpu_ham_backend_resolved) == 'sparse' .and. &
         Natom == sparse_natom .and. Mensemble == sparse_mensemble .and. &
         start_atom == 1 .and. stop_atom == Natom
   end function sparse_backend_can_apply

   !----------------------------------------------------------------------------
   !> Apply persistent CSR J to all three Cartesian RHS columns.
   !>
   !> The sparse backend owns its OpenMP work. Packing and application are
   !> separate measured regions, and no allocation or CSR construction occurs.
   !----------------------------------------------------------------------------
   subroutine apply_sparse_exchange(Mensemble,emomM)
      use omp_lib, only : omp_get_wtime
      implicit none

      integer, intent(in) :: Mensemble
      real(dblprec), dimension(3,sparse_natom,Mensemble), intent(in) :: emomM
      integer :: i,j,k,p
      real(dblprec) :: t_start

      t_start=omp_get_wtime()
      !$omp parallel do default(shared) private(i,k) collapse(2) schedule(static)
      do k=1,Mensemble
         do i=1,sparse_natom
            sparse_rhs(i,1,k)=emomM(1,i,k)
            sparse_rhs(i,2,k)=emomM(2,i,k)
            sparse_rhs(i,3,k)=emomM(3,i,k)
         end do
      end do
      !$omp end parallel do
      sparse_pack_seconds=sparse_pack_seconds+(omp_get_wtime()-t_start)

      t_start=omp_get_wtime()
      !$omp parallel do default(shared) private(i,j,k,p) collapse(2) schedule(static)
      do k=1,Mensemble
         do i=1,sparse_natom
            sparse_field(i,1,k)=0.0_dblprec
            sparse_field(i,2,k)=0.0_dblprec
            sparse_field(i,3,k)=0.0_dblprec
            do p=sparse_row_ptr(i)+1,sparse_row_ptr(i+1)
               j=sparse_columns(p)
               sparse_field(i,1,k)=sparse_field(i,1,k)+sparse_values(p)*sparse_rhs(j,1,k)
               sparse_field(i,2,k)=sparse_field(i,2,k)+sparse_values(p)*sparse_rhs(j,2,k)
               sparse_field(i,3,k)=sparse_field(i,3,k)+sparse_values(p)*sparse_rhs(j,3,k)
            end do
         end do
      end do
      !$omp end parallel do
      sparse_apply_seconds=sparse_apply_seconds+(omp_get_wtime()-t_start)
   end subroutine apply_sparse_exchange

   !----------------------------------------------------------------------------
   !> Expose setup/apply timing for backend measurements without exposing state.
   !----------------------------------------------------------------------------
   subroutine sparse_backend_get_stats(setup_seconds,pack_seconds,apply_seconds,nnz)
      implicit none
      real(dblprec), intent(out) :: setup_seconds,pack_seconds,apply_seconds
      integer(kind=8), intent(out) :: nnz

      setup_seconds=sparse_setup_seconds
      pack_seconds=sparse_pack_seconds
      apply_seconds=sparse_apply_seconds
      nnz=sparse_nnz
   end subroutine sparse_backend_get_stats

   !----------------------------------------------------------------------------
   !> Construct the persistent periodic scalar-J/DMI convolution backend.
   !>
   !> The validated reduced stencil is the sole source for the translational
   !> kernel. FFT plans, mappings, spectra and work buffers all persist until
   !> cleanup; no setup or allocation occurs in effective_field.
   !----------------------------------------------------------------------------
   subroutine setup_convolution_backend(Natom,Mensemble,do_ralloy,do_lsf,NA,N1,N2,N3, &
         BC1,BC2,BC3,do_reduced,nHam)
      implicit none

      integer, intent(in) :: Natom,Mensemble,do_ralloy,NA,N1,N2,N3,nHam
      character(len=1), intent(in) :: do_lsf,BC1,BC2,BC3,do_reduced
      character(len=16) :: requested_backend
      character(len=256) :: backend_diagnostic
      logical :: backend_ok
#ifdef USE_FFTW
      logical :: ok
      character(len=256) :: diagnostic
      real(dblprec) :: t_start,t_end
#endif

      call cleanup_convolution_backend()
      cpu_ham_backend_reason=''
      call resolve_requested_cpu_ham_backend(requested_backend,backend_ok,backend_diagnostic)
      if (.not.backend_ok) then
         cpu_ham_backend_reason=trim(backend_diagnostic)
         return
      endif
      if (trim(requested_backend) /= 'convolution') return

#ifdef USE_FFTW
      if (ham_inp%do_jtensor == 1 .or. &
            ham_inp%do_sa == 1 .or. ham_inp%do_pd == 1 .or. &
            ham_inp%do_biqdm == 1 .or. ham_inp%do_bq == 1 .or. &
            ham_inp%do_ring == 1 .or. ham_inp%do_chir == 1 .or. &
            ham_inp%exc_inter /= 'N' .or. do_ralloy /= 0 .or. do_lsf == 'Y') then
         cpu_ham_backend_reason='convolution supports only periodic reduced scalar-J/DMI without tensor, onsite pair extensions, disorder or LSF'
         return
      endif
      if (.not.allocated(ham%reduced_stencil%record_start) .or. &
            .not.allocated(ham%reduced_stencil%record)) then
         cpu_ham_backend_reason='convolution requires an eligible reduced translational stencil'
         return
      endif
      if (ham_inp%do_dm == 1) then
         if (.not.allocated(ham%dmlistsize) .or. .not.allocated(ham%dmlist) .or. &
               .not.allocated(ham%dm_vect) .or. &
               .not.allocated(ham%reduced_stencil%dmi_record_start)) then
            cpu_ham_backend_reason='convolution DMI data is unavailable for the requested J+D Hamiltonian'
            return
         endif
         ok=cpu_convolution_eligible(Natom,NA,N1,N2,N3,BC1,BC2,BC3,do_reduced, &
            do_ralloy,nHam,ham%aHam,ham%nlistsize,ham%nlist,ham%ncoup,diagnostic, &
            dmlistsize=ham%dmlistsize,dmlist=ham%dmlist,dm_vect=ham%dm_vect)
      else
         ok=cpu_convolution_eligible(Natom,NA,N1,N2,N3,BC1,BC2,BC3,do_reduced, &
            do_ralloy,nHam,ham%aHam,ham%nlistsize,ham%nlist,ham%ncoup,diagnostic)
      endif
      if (.not.ok) then
         cpu_ham_backend_reason='convolution eligibility check failed: '//trim(diagnostic)
         return
      endif

      call cpu_time(t_start)
      ok=cpu_convolution_init(convolution_backend,ham%reduced_stencil,Mensemble,diagnostic)
      if (ok) ok=cpu_convolution_build_kernel(convolution_backend,ham%reduced_stencil,diagnostic)
      if (.not.ok) then
         cpu_ham_backend_reason='convolution setup failed: '//trim(diagnostic)
         call cleanup_convolution_backend()
         return
      endif
      allocate(convolution_field(3,Natom,Mensemble))
      convolution_field=0.0_dblprec
      call cpu_time(t_end)
      convolution_setup_seconds=t_end-t_start
      convolution_backend_ready=.true.
      cpu_ham_backend_resolved='convolution'
      cpu_ham_backend_initialized=.true.
      cpu_ham_backend_reason='eligible periodic reduced scalar-J/DMI; persistent convolution backend initialized'
      write(*,'(2x,a,i0,a,i0,a,a,a,i0)') &
         'Persistent scalar-J CPU convolution ready: ',Natom,' atoms, ',NA, &
         ' basis, provider ',trim(cpu_fft_provider_name()),' threads=', &
         cpu_fft_provider_threads()
#else
      cpu_ham_backend_reason='CPU convolution support is not compiled into this executable; rebuild with FFTW CPU support (-DUSE_FFTW=ON), or select cpu_ham_backend direct'
#endif
   end subroutine setup_convolution_backend


   subroutine cleanup_convolution_backend()
      implicit none

#ifdef USE_FFTW
      call cpu_convolution_clear(convolution_backend)
      if (allocated(convolution_field)) deallocate(convolution_field)
      convolution_backend_ready=.false.
      convolution_setup_seconds=0.0_dblprec
      if (trim(cpu_ham_backend_resolved) == 'convolution') cpu_ham_backend_resolved='direct'
      cpu_ham_backend_initialized=.false.
#endif
   end subroutine cleanup_convolution_backend


   logical function convolution_backend_can_apply(Natom,Mensemble,start_atom,stop_atom)
      implicit none
      integer, intent(in) :: Natom,Mensemble,start_atom,stop_atom

      convolution_backend_can_apply=.false.
#ifdef USE_FFTW
      convolution_backend_can_apply=convolution_backend_ready .and. trim(cpu_ham_backend_resolved) == 'convolution' .and. &
         Natom == convolution_backend%natom .and. Mensemble == convolution_backend%ensembles .and. &
         start_atom == 1 .and. stop_atom == Natom
#endif
   end function convolution_backend_can_apply


   subroutine apply_convolution_exchange(Mensemble,emomM)
      implicit none
      integer, intent(in) :: Mensemble
#ifdef USE_FFTW
      real(dblprec), dimension(3,convolution_backend%natom,Mensemble), intent(in) :: emomM
      logical :: ok
      character(len=256) :: diagnostic

      ok=cpu_convolution_apply(convolution_backend,emomM,convolution_field,diagnostic)
      if (.not.ok) error stop 'CPU convolution apply failed: '//trim(diagnostic)
#else
      real(dblprec), intent(in) :: emomM(:,:,:)
#endif
   end subroutine apply_convolution_exchange


   subroutine convolution_backend_get_stats(setup_seconds,pack_seconds,forward_seconds, &
         spectral_seconds,inverse_seconds,unpack_seconds,apply_seconds,apply_count)
      implicit none
      real(dblprec), intent(out) :: setup_seconds,pack_seconds,forward_seconds
      real(dblprec), intent(out) :: spectral_seconds,inverse_seconds,unpack_seconds,apply_seconds
      integer(kind=8), intent(out) :: apply_count

      setup_seconds=0.0_dblprec
      pack_seconds=0.0_dblprec
      forward_seconds=0.0_dblprec
      spectral_seconds=0.0_dblprec
      inverse_seconds=0.0_dblprec
      unpack_seconds=0.0_dblprec
      apply_seconds=0.0_dblprec
      apply_count=0_8
#ifdef USE_FFTW
      setup_seconds=convolution_setup_seconds
      call cpu_convolution_get_stats(convolution_backend,pack_seconds,forward_seconds, &
         spectral_seconds,inverse_seconds,unpack_seconds,apply_seconds,apply_count)
#endif
   end subroutine convolution_backend_get_stats

end module HamiltonianBackend
