! DMI-DIMER-ENERGY: an analytic convention oracle, deliberately independent
! of the production implementations under test.
program test_dmi_dimer_energy
   use Parameters, only : dblprec
   use Constants, only : mub
   use HamiltonianData, only : ham
   use HamiltonianActions, only : dzyaloshinskii_moriya_field
   use ApplyHamiltonian, only : heisge
   implicit none

   integer, parameter :: natom = 2, ensembles = 1, max_neigh = 1
   real(dblprec), parameter :: d = 2.5_dblprec
   real(dblprec), parameter :: tolerance = 1.0d-13
   integer :: failures
   integer :: nlistsize(natom), nlist(max_neigh,natom), dmlistsize(natom), dmlist(max_neigh,natom)
   integer :: taniso(natom), pdlist(max_neigh,natom), pdlistsize(natom)
   integer :: biqdmlist(max_neigh,natom), biqdmlistsize(natom), bqlist(max_neigh,natom), bqlistsize(natom)
   real(dblprec) :: emom_m(3,natom,ensembles), emom(3,natom,ensembles)
   real(dblprec) :: dm_vect(3,max_neigh,natom), ncoup(max_neigh,natom)
   real(dblprec) :: pd_vect(6,max_neigh,natom), biqdm_vect(1,max_neigh,natom), j_bq(max_neigh,natom)
   real(dblprec) :: eaniso(3,natom), kaniso(2,natom), sb(natom), external_field(3,natom,ensembles)
   real(dblprec) :: action_field(3), apply_field(3,natom,ensembles), apply_internal(3,natom,ensembles)
   real(dblprec) :: apply_external(3,natom,ensembles), expected_field(3,natom), energy_j

   failures = 0
   emom_m = 0.0_dblprec
   emom_m(:,1,1) = (/1.0_dblprec,0.0_dblprec,0.0_dblprec/)
   emom_m(:,2,1) = (/0.0_dblprec,1.0_dblprec,0.0_dblprec/)
   emom = emom_m

   ! The directed dimer has D_12=+D z and D_21=-D z.  Hence the exact
   ! energy is mu_B D and B_1=-D x, B_2=-D y.
   dm_vect = 0.0_dblprec
   dm_vect(3,1,1) = d
   dm_vect(3,1,2) = -d
   dmlist(1,:) = (/2,1/)
   dmlistsize = 1
   expected_field = 0.0_dblprec
   expected_field(:,1) = (/-d,0.0_dblprec,0.0_dblprec/)
   expected_field(:,2) = (/0.0_dblprec,-d,0.0_dblprec/)
   energy_j = mub*d

   call check_close(0.5_dblprec*mub*(dot_product(dm_vect(:,1,1),cross(emom_m(:,1,1),emom_m(:,2,1))) + &
      dot_product(dm_vect(:,1,2),cross(emom_m(:,2,1),emom_m(:,1,1)))),energy_j,tolerance*mub, &
      'hand-derived directed-pair DMI energy',failures)

   call initialise_action_hamiltonian(dm_vect,dmlist,dmlistsize)
   action_field = 0.0_dblprec
   call dzyaloshinskii_moriya_field(1,1,action_field,natom,ensembles,emom_m)
   call check_vector(action_field,expected_field(:,1),tolerance, &
      'HamiltonianActions DMI field uses D_ij cross M_j',failures)
   call clear_action_hamiltonian()

   nlistsize = 0
   nlist = 1
   ncoup = 0.0_dblprec
   taniso = 0
   pdlist = 1
   pdlistsize = 0
   pd_vect = 0.0_dblprec
   biqdmlist = 1
   biqdmlistsize = 0
   biqdm_vect = 0.0_dblprec
   bqlist = 1
   bqlistsize = 0
   j_bq = 0.0_dblprec
   eaniso = 0.0_dblprec
   kaniso = 0.0_dblprec
   sb = 0.0_dblprec
   external_field = 0.0_dblprec
   call heisge(natom,ensembles,max_neigh,ncoup,nlist,nlistsize,1,max_neigh,dm_vect,dmlist,dmlistsize, &
      0,max_neigh,pd_vect,pdlist,pdlistsize,0,max_neigh,biqdm_vect,biqdmlist,biqdmlistsize, &
      0,max_neigh,j_bq,bqlist,bqlistsize,taniso,eaniso,kaniso,sb,apply_field,apply_internal,apply_external, &
      emom_m,emom,external_field,0)
   call check_vector(apply_internal(:,1,1),expected_field(:,1),tolerance, &
      'ApplyHamiltonian DMI field uses D_ij cross M_j',failures)
   call check_vector(apply_internal(:,2,1),expected_field(:,2),tolerance, &
      'ApplyHamiltonian second dimer field uses D_ji cross M_i',failures)

   if (failures /= 0) error stop 'DMI-DIMER-ENERGY failed'
   write(*,'(a)') 'DMI-DIMER-ENERGY passed'

contains

   pure function cross(a,b) result(value)
      real(dblprec), intent(in) :: a(3), b(3)
      real(dblprec) :: value(3)
      value = (/a(2)*b(3)-a(3)*b(2), a(3)*b(1)-a(1)*b(3), a(1)*b(2)-a(2)*b(1)/)
   end function cross

   subroutine initialise_action_hamiltonian(vectors,neighbours,sizes)
      real(dblprec), intent(in) :: vectors(3,max_neigh,natom)
      integer, intent(in) :: neighbours(max_neigh,natom), sizes(natom)
      allocate(ham%aHam(natom),ham%dmlistsize(natom),ham%dmlist(max_neigh,natom),ham%dm_vect(3,max_neigh,natom))
      ham%aHam = (/1,2/)
      ham%dmlistsize = sizes
      ham%dmlist = neighbours
      ham%dm_vect = vectors
   end subroutine initialise_action_hamiltonian

   subroutine clear_action_hamiltonian()
      deallocate(ham%aHam,ham%dmlistsize,ham%dmlist,ham%dm_vect)
   end subroutine clear_action_hamiltonian

   subroutine check_close(actual,expected,epsilon,label,count)
      real(dblprec), intent(in) :: actual, expected, epsilon
      character(len=*), intent(in) :: label
      integer, intent(inout) :: count
      if (abs(actual-expected) > epsilon) then
         write(*,'(a,2(1x,es24.16))') trim(label)//' expected/actual:',expected,actual
         count = count + 1
      end if
   end subroutine check_close

   subroutine check_vector(actual,expected,epsilon,label,count)
      real(dblprec), intent(in) :: actual(3), expected(3), epsilon
      character(len=*), intent(in) :: label
      integer, intent(inout) :: count
      if (maxval(abs(actual-expected)) > epsilon) then
         write(*,'(a,3(1x,es16.8))') trim(label)//' expected:',expected
         write(*,'(a,3(1x,es16.8))') 'actual:',actual
         count = count + 1
      end if
   end subroutine check_vector

end program test_dmi_dimer_energy
