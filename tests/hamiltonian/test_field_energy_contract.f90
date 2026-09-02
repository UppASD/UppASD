! CPU-HAM-01A: canonical field/energy execution contract.
program test_field_energy_contract
   use Parameters, only : dblprec
   use Constants, only : mub,mry
   use HamiltonianData, only : ham
   use HamiltonianActions, only : effective_field
   use InputData, only : ham_inp
   implicit none

   integer, parameter :: natom=3, nensemble=2, nmacro=1, max_neigh=1
   real(dblprec), parameter :: tolerance=0.0_dblprec
   integer :: i, j, jj, failures
   integer :: cell_index(natom), macro_nlistsize(nmacro)
   real(dblprec) :: emomM(3,natom,nensemble), mmom(natom,nensemble)
   real(dblprec) :: external_field(3,natom,nensemble), time_external_field(3,natom,nensemble)
   real(dblprec) :: emomM_macro(3,nmacro,nensemble)
   real(dblprec) :: beff_on(3,natom,nensemble), beff_off(3,natom,nensemble)
   real(dblprec) :: beff_j(3,natom,nensemble), beff1_j(3,natom,nensemble)
   real(dblprec) :: beff2_j(3,natom,nensemble), j_reference(3,natom,nensemble)
   real(dblprec) :: beff_default(3,natom,nensemble)
   real(dblprec) :: beff1_on(3,natom,nensemble), beff1_off(3,natom,nensemble)
   real(dblprec) :: beff1_default(3,natom,nensemble)
   real(dblprec) :: beff2_on(3,natom,nensemble), beff2_off(3,natom,nensemble)
   real(dblprec) :: beff2_default(3,natom,nensemble)
   real(dblprec) :: energy_on, energy_off, energy_default, field_energy

   failures=0
   call setup_hamiltonian()

   mmom=1.0_dblprec
   emomM(:,1,1)=(/0.8_dblprec,0.3_dblprec,0.5_dblprec/)
   emomM(:,2,1)=(/-0.2_dblprec,0.9_dblprec,0.4_dblprec/)
   emomM(:,3,1)=(/0.6_dblprec,-0.7_dblprec,0.1_dblprec/)
   emomM(:,:,2)=-emomM(:,:,1)
   external_field=0.0_dblprec
   external_field(1,:,:)=0.11_dblprec
   time_external_field=0.0_dblprec
   time_external_field(2,:,:)=0.07_dblprec
   emomM_macro=0.0_dblprec

   ! Exercise the production scalar-J path directly with non-collinear
   ! moments and multiple ensembles before enabling the DMI part below.
   ham_inp%do_dm=0
   call effective_field(natom,nensemble,1,natom,emomM,mmom,external_field, &
      time_external_field,beff_j,beff1_j,beff2_j,energy_off,nmacro,cell_index, &
      emomM_macro,macro_nlistsize,1,1,1,1,measure_energy=.false.)
   j_reference=0.0_dblprec
   do j=1,nensemble
      do i=1,natom
         j_reference(:,i,j)=0.0_dblprec
         do jj=1,ham%nlistsize(ham%aHam(i))
            j_reference(1,i,j)=j_reference(1,i,j)+ham%ncoup(jj,ham%aHam(i),1)* &
               emomM(1,ham%nlist(jj,i),j)
            j_reference(2,i,j)=j_reference(2,i,j)+ham%ncoup(jj,ham%aHam(i),1)* &
               emomM(2,ham%nlist(jj,i),j)
            j_reference(3,i,j)=j_reference(3,i,j)+ham%ncoup(jj,ham%aHam(i),1)* &
               emomM(3,ham%nlist(jj,i),j)
         end do
      end do
   end do
   call check(maxval(abs(beff1_j-j_reference)) <= tolerance, &
      'scalar-J field matches the direct neighbour-list oracle')
   ham_inp%do_dm=1

   call effective_field(natom,nensemble,1,natom,emomM,mmom,external_field, &
      time_external_field,beff_on,beff1_on,beff2_on,energy_on,nmacro,cell_index, &
      emomM_macro,macro_nlistsize,1,1,1,1,measure_energy=.true.)
   call effective_field(natom,nensemble,1,natom,emomM,mmom,external_field, &
      time_external_field,beff_off,beff1_off,beff2_off,energy_off,nmacro,cell_index, &
      emomM_macro,macro_nlistsize,1,1,1,1,measure_energy=.false.)
   call effective_field(natom,nensemble,1,natom,emomM,mmom,external_field, &
      time_external_field,beff_default,beff1_default,beff2_default,energy_default,nmacro, &
      cell_index,emomM_macro,macro_nlistsize,1,1,1,1)

   call check(maxval(abs(beff_on-beff_off)) <= tolerance, &
      'total field is unchanged when energy is disabled')
   call check(maxval(abs(beff1_on-beff1_off)) <= tolerance, &
      'internal field is unchanged when energy is disabled')
   call check(maxval(abs(beff2_on-beff2_off)) <= tolerance, &
      'external field is unchanged when energy is disabled')
   call check(energy_off == 0.0_dblprec, 'field-only call returns zero energy')
   call check(energy_on /= 0.0_dblprec, 'energy-enabled call produces energy')
   field_energy=0.0_dblprec
   do i=1,natom
      field_energy=field_energy-0.5_dblprec*sum(emomM(:,i,:)*beff1_on(:,i,:))
      field_energy=field_energy-sum(emomM(:,i,:)*(external_field(:,i,:)+time_external_field(:,i,:)))
   end do
   field_energy=field_energy*mub/mry
   call check(abs(energy_on-field_energy) <= 1.0d-12, &
      'bilinear pair energy follows -1/2 moment dot canonical pair field')
   call check(maxval(abs(beff_on-beff_default)) <= tolerance, &
      'legacy default field behavior is preserved')
   call check(abs(energy_on-energy_default) <= 1.0d-12, &
      'legacy default energy behavior is preserved')

   ! Also exercise the dipole manager's energy gate, using a small dense oracle.
   ham_inp%do_dip=1
   ham%Qdip=0.0_dblprec
   do i=1,natom
      ham%Qdip(1,1,i,i)=0.13_dblprec
      ham%Qdip(2,2,i,i)=0.09_dblprec
      ham%Qdip(3,3,i,i)=0.17_dblprec
   end do
   call effective_field(natom,nensemble,1,natom,emomM,mmom,external_field, &
      time_external_field,beff_on,beff1_on,beff2_on,energy_on,nmacro,cell_index, &
      emomM_macro,macro_nlistsize,1,1,1,1,measure_energy=.true.)
   call effective_field(natom,nensemble,1,natom,emomM,mmom,external_field, &
      time_external_field,beff_off,beff1_off,beff2_off,energy_off,nmacro,cell_index, &
      emomM_macro,macro_nlistsize,1,1,1,1,measure_energy=.false.)
   call check(maxval(abs(beff_on-beff_off)) <= tolerance, &
      'dipole field is unchanged when energy is disabled')
   call check(energy_off == 0.0_dblprec, 'dipole field-only call returns zero energy')
   call check(energy_on /= 0.0_dblprec, 'dipole energy-enabled call produces energy')
   field_energy=0.0_dblprec
   do i=1,natom
      field_energy=field_energy-0.5_dblprec*sum(emomM(:,i,:)*(beff_on(:,i,:)-beff1_on(:,i,:)-beff2_on(:,i,:)))
      field_energy=field_energy-0.5_dblprec*sum(emomM(:,i,:)*beff1_on(:,i,:))
      field_energy=field_energy-sum(emomM(:,i,:)*(external_field(:,i,:)+time_external_field(:,i,:)))
   end do
   field_energy=field_energy*mub/mry
   call check(abs(energy_on-field_energy) <= 1.0d-12, &
      'dipolar pair energy follows -1/2 moment dot canonical pair field')

   call clear_hamiltonian()
   if (failures /= 0) error stop 'CPU-HAM-01A field/energy contract failed'
   write(*,'(a)') 'CPU-HAM-01A field/energy contract passed'

contains

   subroutine setup_hamiltonian()
      integer :: j

      ham_inp%do_jtensor=0
      ham_inp%exc_inter='N'
      ham_inp%do_dm=1
      ham_inp%do_sa=0
      ham_inp%do_pd=0
      ham_inp%do_biqdm=0
      ham_inp%do_bq=0
      ham_inp%do_ring=0
      ham_inp%do_chir=0
      ham_inp%do_anisotropy=0
      ham_inp%do_dip=0
      ham_inp%mult_axis='N'

      allocate(ham%aHam(natom),ham%nlistsize(natom),ham%nlist(max_neigh,natom), &
         ham%ncoup(max_neigh,natom,1),ham%dmlistsize(natom),ham%dmlist(max_neigh,natom), &
         ham%dm_vect(3,max_neigh,natom),ham%Qdip(3,3,natom,natom), &
         ham%Qdip_macro(3,3,nmacro,nmacro))
      ham%aHam=(/1,2,3/)
      ham%nlistsize=1
      ham%nlist(:,1)=(/2/)
      ham%nlist(:,2)=(/3/)
      ham%nlist(:,3)=(/1/)
      ham%ncoup(:,:,1)=0.23_dblprec
      ham%dmlistsize=1
      ham%dmlist(:,1)=(/2/)
      ham%dmlist(:,2)=(/3/)
      ham%dmlist(:,3)=(/1/)
      ham%dm_vect=0.0_dblprec
      ham%dm_vect(:,1,1)=(/0.17_dblprec,-0.08_dblprec,0.11_dblprec/)
      ham%dm_vect(:,1,2)=(/-0.04_dblprec,0.13_dblprec,0.06_dblprec/)
      ham%dm_vect(:,1,3)=(/0.09_dblprec,0.05_dblprec,-0.12_dblprec/)
      ham%Qdip=0.0_dblprec
      ham%Qdip_macro=0.0_dblprec
      cell_index=1
      macro_nlistsize=natom
      do j=1,natom
         if (ham%nlist(1,j) < 1) error stop 'invalid test neighbour'
      end do
   end subroutine setup_hamiltonian

   subroutine clear_hamiltonian()
      deallocate(ham%aHam,ham%nlistsize,ham%nlist,ham%ncoup,ham%dmlistsize, &
         ham%dmlist,ham%dm_vect,ham%Qdip,ham%Qdip_macro)
   end subroutine clear_hamiltonian

   subroutine check(condition,label)
      logical, intent(in) :: condition
      character(len=*), intent(in) :: label
      if (.not. condition) then
         write(*,'(a)') 'FAILED: '//trim(label)
         failures=failures+1
      end if
   end subroutine check

end program test_field_energy_contract
