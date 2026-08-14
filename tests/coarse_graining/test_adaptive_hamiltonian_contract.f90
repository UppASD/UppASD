! RCG-09A.1: setup-contract regression for the unique-pair production fold.
! This calls the predicates used by adaptive setup; it does not duplicate the
! folding rule in a separate analytical oracle.
program test_adaptive_hamiltonian_contract
   use Parameters, only : dblprec
   use AdaptiveCGProduction, only : ADAPTIVE_CG_PRODUCTION_OK, &
      ADAPTIVE_CG_PRODUCTION_REJECTED, validate_adaptive_ensemble_moments, &
      validate_adaptive_folded_pair_contract
   implicit none

   integer, parameter :: pairs = 1
   integer :: status, failures
   integer :: pair_i(pairs), pair_j(pairs), exchange_forward_count(pairs), &
      exchange_reverse_count(pairs), dmi_forward_count(pairs), dmi_reverse_count(pairs)
   real(dblprec) :: exchange_forward(pairs), exchange_reverse(pairs), &
      dmi_forward(3,pairs), dmi_reverse(3,pairs), moments(2,2)
   character(len=512) :: diagnostic

   failures = 0
   pair_i = 1; pair_j = 2
   exchange_forward_count = 1; exchange_reverse_count = 1
   dmi_forward_count = 1; dmi_reverse_count = 1
   exchange_forward = 0.17_dblprec; exchange_reverse = 0.17_dblprec
   dmi_forward(:,1) = (/0.37_dblprec,-0.29_dblprec,0.61_dblprec/)
   ! `dmi_reverse` is canonicalized by the production setup: -D_21.
   dmi_reverse = dmi_forward

   call validate_adaptive_folded_pair_contract(pairs,pair_i,pair_j, &
      exchange_forward_count,exchange_reverse_count,dmi_forward_count, &
      dmi_reverse_count,exchange_forward,exchange_reverse,dmi_forward, &
      dmi_reverse,status,diagnostic)
   call check(status == ADAPTIVE_CG_PRODUCTION_OK, &
      'reciprocal J and D directed pair is accepted')

   dmi_reverse(3,1) = dmi_reverse(3,1) + 0.01_dblprec
   call validate_adaptive_folded_pair_contract(pairs,pair_i,pair_j, &
      exchange_forward_count,exchange_reverse_count,dmi_forward_count, &
      dmi_reverse_count,exchange_forward,exchange_reverse,dmi_forward, &
      dmi_reverse,status,diagnostic)
   call check(status == ADAPTIVE_CG_PRODUCTION_REJECTED .and. &
      index(diagnostic,'reciprocal DMI contract') > 0, &
      'D_ji /= -D_ij is rejected by the folded-pair contract')
   dmi_reverse = dmi_forward

   dmi_forward_count = 2
   call validate_adaptive_folded_pair_contract(pairs,pair_i,pair_j, &
      exchange_forward_count,exchange_reverse_count,dmi_forward_count, &
      dmi_reverse_count,exchange_forward,exchange_reverse,dmi_forward, &
      dmi_reverse,status,diagnostic)
   call check(status == ADAPTIVE_CG_PRODUCTION_REJECTED .and. &
      index(diagnostic,'periodic-image alias') > 0, &
      'multiple physical DMI entries for one atom pair are rejected')
   dmi_forward_count = 1

   exchange_reverse = 0.18_dblprec
   call validate_adaptive_folded_pair_contract(pairs,pair_i,pair_j, &
      exchange_forward_count,exchange_reverse_count,dmi_forward_count, &
      dmi_reverse_count,exchange_forward,exchange_reverse,dmi_forward, &
      dmi_reverse,status,diagnostic)
   call check(status == ADAPTIVE_CG_PRODUCTION_REJECTED .and. &
      index(diagnostic,'reciprocal exchange contract') > 0, &
      'J_ji /= J_ij is rejected by the folded-pair contract')

   moments = 1.0_dblprec
   moments(:,2) = (/1.0_dblprec,1.7_dblprec/)
   moments(:,1) = moments(:,2)
   call validate_adaptive_ensemble_moments(moments,status,diagnostic)
   call check(status == ADAPTIVE_CG_PRODUCTION_OK, &
      'identical per-atom ensemble moments are accepted')
   moments(2,2) = 1.71_dblprec
   call validate_adaptive_ensemble_moments(moments,status,diagnostic)
   call check(status == ADAPTIVE_CG_PRODUCTION_REJECTED .and. &
      index(diagnostic,'ensemble moments') > 0, &
      'different ensemble moment magnitude is rejected explicitly')

   if (failures /= 0) error stop 'adaptive Hamiltonian contract failed'
   write(*,'(a)') 'adaptive Hamiltonian contract passed'

contains

   subroutine check(condition,label)
      logical, intent(in) :: condition
      character(len=*), intent(in) :: label
      if (.not. condition) then
         write(*,'(a)') 'FAILED: '//trim(label)
         failures = failures + 1
      end if
   end subroutine check

end program test_adaptive_hamiltonian_contract
