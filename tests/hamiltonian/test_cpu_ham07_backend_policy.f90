! CPU-HAM-07: explicit CPU pair-backend policy and compatibility aliases.
program test_cpu_ham07_backend_policy

   use HamiltonianActions, only : resolve_cpu_ham_backend, &
      resolve_requested_cpu_ham_backend
   use InputData, only : cpu_ham_backend, do_sparse, do_convolution

   implicit none

   integer :: failures

   failures=0
   call check_public_name('direct','direct')
   call check_public_name('DIRECT','direct')
   call check_public_name('sparse','sparse')
   call check_public_name('convolution','convolution')
   call check_rejected_name('auto','AUTO is disabled')
   call check_rejected_name('not-a-backend','unsupported CPU Hamiltonian backend')

   cpu_ham_backend='direct'
   do_sparse='N'
   do_convolution='N'
   call check_requested('direct','direct','explicit DIRECT request')

   ! Existing benchmark/input files continue to resolve their legacy flags,
   ! but the production setup still applies the same capability checks.
   do_sparse='Y'
   call check_requested('direct','sparse','legacy do_sparse compatibility alias')
   do_sparse='N'
   do_convolution='Y'
   call check_requested('direct','convolution','legacy do_convolution compatibility alias')

   do_sparse='Y'
   call check_requested_rejected('direct','two CPU Hamiltonian backends')
   do_sparse='N'
   cpu_ham_backend='sparse'
   do_convolution='Y'
   call check_requested_rejected('sparse','conflicts with do_convolution')
   do_convolution='N'
   cpu_ham_backend='convolution'
   do_sparse='Y'
   call check_requested_rejected('convolution','conflicts with do_sparse')

   cpu_ham_backend='direct'
   do_sparse='N'
   do_convolution='N'
   if (failures /= 0) error stop 'CPU-HAM-07 backend policy tests failed'
   write(*,'(a)') 'CPU-HAM-07 backend policy tests passed'

contains

   subroutine check_public_name(requested,expected)
      character(len=*), intent(in) :: requested,expected
      character(len=32) :: resolved
      character(len=256) :: diagnostic
      logical :: ok

      call resolve_cpu_ham_backend(requested,resolved,ok,diagnostic)
      call check(ok .and. trim(resolved)==expected, &
         'public backend '//trim(requested)//' resolves to '//trim(expected))
   end subroutine check_public_name


   subroutine check_rejected_name(requested,expected_text)
      character(len=*), intent(in) :: requested,expected_text
      character(len=32) :: resolved
      character(len=256) :: diagnostic
      logical :: ok

      call resolve_cpu_ham_backend(requested,resolved,ok,diagnostic)
      call check(.not.ok .and. index(diagnostic,expected_text)>0, &
         'invalid public backend '//trim(requested)//' is rejected clearly')
   end subroutine check_rejected_name


   subroutine check_requested(requested,expected,expected_text)
      character(len=*), intent(in) :: requested,expected,expected_text
      character(len=32) :: resolved
      character(len=256) :: diagnostic
      logical :: ok

      cpu_ham_backend=requested
      call resolve_requested_cpu_ham_backend(resolved,ok,diagnostic)
      call check(ok .and. trim(resolved)==expected .and. index(diagnostic,expected_text)>0, &
         'requested backend resolves to '//trim(expected))
   end subroutine check_requested


   subroutine check_requested_rejected(requested,expected_text)
      character(len=*), intent(in) :: requested,expected_text
      character(len=32) :: resolved
      character(len=256) :: diagnostic
      logical :: ok

      cpu_ham_backend=requested
      call resolve_requested_cpu_ham_backend(resolved,ok,diagnostic)
      call check(.not.ok .and. index(diagnostic,expected_text)>0, &
         'conflicting backend request is rejected clearly')
   end subroutine check_requested_rejected


   subroutine check(condition,label)
      logical, intent(in) :: condition
      character(len=*), intent(in) :: label

      if (.not.condition) then
         write(*,'(a)') 'FAILED: '//trim(label)
         failures=failures+1
      endif
   end subroutine check

end program test_cpu_ham07_backend_policy
