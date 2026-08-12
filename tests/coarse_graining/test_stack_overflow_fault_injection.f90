!-------------------------------------------------------------------------------
!> RCG-06A (F-13) negative control: reproduces, in isolation, exactly the
!> automatic-(stack)-array shape and size that
!> StaticHybridOperator's evaluate_static_hybrid_operator declared before
!> RCG-06A hoisted it into persistent operator-owned scratch. n_atoms is a
!> runtime dummy argument, not a compile-time constant, so the compiler
!> cannot promote the arrays to static/global storage -- this must allocate
!> on the stack exactly as the pre-fix production code did.
!>
!> Not linked against asdlib; this is a minimal standalone reproduction, not
!> a call into production source, so it stays valid evidence even as
!> production code around it changes.
!>
!> Driven by tests/coarse_graining/run_mem_large_host.py under a constrained
!> `ulimit -s`, where it is expected to be killed (SIGSEGV) -- see
!> tests/coarse_graining/e2e/mem_large_host/README.md.
!-------------------------------------------------------------------------------
program test_stack_overflow_fault_injection
   implicit none
   integer, parameter :: dblprec = kind(1.0d0)
   integer :: n_atoms

   n_atoms = 8000
   call touch_automatic_arrays(n_atoms)
   write(*,'(a)') 'test_stack_overflow_fault_injection: did not crash'
contains

   subroutine touch_automatic_arrays(n_atoms)
      integer, intent(in) :: n_atoms
      real(dblprec) :: effective_direction(3,n_atoms)
      real(dblprec) :: ghost_direction(3,n_atoms)
      real(dblprec) :: atomistic_field(3,n_atoms)
      real(dblprec) :: ghost_reaction_field(3,n_atoms)
      integer :: atom

      do atom = 1, n_atoms
         effective_direction(:,atom) = real(atom,dblprec)
         ghost_direction(:,atom) = real(atom,dblprec)
         atomistic_field(:,atom) = real(atom,dblprec)
         ghost_reaction_field(:,atom) = real(atom,dblprec)
      end do
      write(*,'(a,es24.16)') 'checksum=', &
         sum(effective_direction)+sum(ghost_direction)+ &
         sum(atomistic_field)+sum(ghost_reaction_field)
   end subroutine touch_automatic_arrays

end program test_stack_overflow_fault_injection
