! CPU-HAM-11: canonical global term decomposition versus local MC DeltaE.
program test_cpu_ham11_energy
   use Parameters, only : dblprec
   use Constants, only : mub,mry
   use HamiltonianData, only : ham
   use HamiltonianActions, only : effective_field, canonical_field_energy, &
      canonical_onsite_energy, HAM_TERM_EXCHANGE, HAM_TERM_DMI, HAM_TERM_ANISOTROPY, &
      HAM_TERM_COUNT
   use HamiltonianBackend, only : set_reduced_direct_testing
   use InputData, only : ham_inp, cpu_ham_backend, do_sparse, do_convolution
   use montecarlo_common, only : calculate_energy => calculate_energy
   implicit none

   integer, parameter :: natom=4, nensemble=1, nmacro=1, max_neigh=natom-1, max_macro=natom
   real(dblprec), parameter :: tolerance=1.0d-12
   integer :: i, j, x, failures
   integer :: cell_index(natom), macro_nlistsize(nmacro), macro_atom_nlist(nmacro,max_macro)
   integer :: iflip, icell
   real(dblprec) :: emomM(3,natom,nensemble), trialM(3,natom,nensemble)
   real(dblprec) :: emom(3,natom,nensemble), mmom(natom,nensemble)
   real(dblprec) :: extfield(3), external_field(3,natom,nensemble)
   real(dblprec) :: time_external_field(3,natom,nensemble)
   real(dblprec) :: emomM_macro(3,nmacro,nensemble), macro_trial(3), macro_mag_trial
   real(dblprec) :: beff(3,natom,nensemble), beff1(3,natom,nensemble)
   real(dblprec) :: beff2(3,natom,nensemble), effective_energy
   real(dblprec) :: terms(3,HAM_TERM_COUNT,natom,nensemble)
   real(dblprec) :: before_raw, after_raw, delta_global_raw, delta_mc_raw
   real(dblprec) :: de_mc
   real(dblprec) :: newmom(3), norm
   integer :: mode_dm, mode_anis

   failures=0
   call setup_hamiltonian()

   mmom=1.0_dblprec
   emomM(:,1,1)=(/0.817d0,0.231d0,-0.529d0/)
   emomM(:,2,1)=(/-0.347d0,0.861d0,0.371d0/)
   emomM(:,3,1)=(/0.592d0,-0.671d0,0.445d0/)
   emomM(:,4,1)=(/-0.713d0,-0.192d0,0.674d0/)
   do i=1,natom
      norm=sqrt(sum(emomM(:,i,1)**2))
      emomM(:,i,1)=emomM(:,i,1)/norm
   enddo
   emom=emomM
   external_field=0.0_dblprec
   time_external_field=0.0_dblprec
   emomM_macro=0.0_dblprec
   extfield=0.0_dblprec
   newmom=(/0.291d0,-0.816d0,0.499d0/)
   newmom=newmom/sqrt(sum(newmom**2))
   iflip=3

   ! First validate the pure bilinear exchange and DMI paths independently.
   do mode_anis=0,1
      if (mode_anis==1) then
         ham_inp%do_anisotropy=1
         ham%taniso=1
         ham%eaniso=0.0_dblprec
         ham%eaniso(3,:)=1.0_dblprec
         ham%kaniso=0.0_dblprec
         ham%kaniso(1,:)=0.19_dblprec
         ham%kaniso(2,:)=0.031_dblprec
      else
         ham_inp%do_anisotropy=0
      endif
      do mode_dm=0,1
         ham_inp%do_dm=mode_dm
         call global_raw_energy(emomM,before_raw)
         trialM=emomM
         trialM(:,iflip,1)=newmom
         call global_raw_energy(trialM,after_raw)
         delta_global_raw=after_raw-before_raw
         call mc_delta(mode_dm,mode_anis,de_mc)
         delta_mc_raw=de_mc/mub
         if (mode_dm==0 .and. mode_anis==0) then
            call check(abs(delta_global_raw-delta_mc_raw)<=tolerance, &
               'exchange MC DeltaE matches canonical global before/after energy')
            call check(abs((-delta_mc_raw)-delta_global_raw)>1.0d-6, &
               'exchange sign mutation is discriminating')
         elseif (mode_dm==1 .and. mode_anis==0) then
            call check(abs(delta_global_raw-delta_mc_raw)<=tolerance, &
               'DMI MC DeltaE matches canonical global before/after energy')
            call check(abs(2.0_dblprec*delta_mc_raw-delta_global_raw)>1.0d-6, &
               'DMI factor mutation is discriminating')
         elseif (mode_dm==0 .and. mode_anis==1) then
            call check(abs(delta_global_raw-delta_mc_raw)<=tolerance, &
               'onsite anisotropy MC DeltaE matches canonical global energy')
         endif
      enddo
   enddo

   call clear_hamiltonian()
   if (failures/=0) error stop 'CPU-HAM-11 energy convention parity failed'
   write(*,'(a)') 'CPU-HAM-11 energy convention parity passed'

contains

   subroutine setup_hamiltonian()
      integer :: ih
      real(dblprec) :: d(3)

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
      ham_inp%mult_axis='N'
      cpu_ham_backend='direct'
      do_sparse='N'
      do_convolution='N'
      call set_reduced_direct_testing(.false.)

      allocate(ham%aHam(natom),ham%nlistsize(natom),ham%nlist(max_neigh,natom), &
         ham%ncoup(max_neigh,natom,1),ham%dmlistsize(natom),ham%dmlist(max_neigh,natom), &
         ham%dm_vect(3,max_neigh,natom),ham%taniso(natom),ham%taniso_diff(natom), &
         ham%sb(natom),ham%sb_diff(natom),ham%kaniso(2,natom),ham%kaniso_diff(2,natom), &
         ham%eaniso(3,natom),ham%eaniso_diff(3,natom))
      ham%aHam=(/(ih,ih=1,natom)/)
      ham%nlistsize=max_neigh
      ham%dmlistsize=max_neigh
      ham%taniso=0
      ham%taniso_diff=0
      ham%sb=0.0_dblprec
      ham%sb_diff=0.0_dblprec
      ham%kaniso=0.0_dblprec
      ham%kaniso_diff=0.0_dblprec
      ham%eaniso=0.0_dblprec
      ham%eaniso_diff=0.0_dblprec

      do ih=1,natom
         j=0
         do x=1,natom
            if (x==ih) cycle
            j=j+1
            ham%nlist(j,ih)=x
            ham%dmlist(j,ih)=x
            ham%ncoup(j,ih,1)=0.27_dblprec+0.013_dblprec*real(min(ih,x),dblprec)
            if (ih<x) then
               d=(/0.071_dblprec+0.009_dblprec*real(ih,dblprec), &
                  -0.043_dblprec+0.005_dblprec*real(x,dblprec), &
                  0.058_dblprec-0.004_dblprec*real(ih+x,dblprec)/)
            else
               d=-(/0.071_dblprec+0.009_dblprec*real(x,dblprec), &
                  -0.043_dblprec+0.005_dblprec*real(ih,dblprec), &
                  0.058_dblprec-0.004_dblprec*real(ih+x,dblprec)/)
            endif
            ham%dm_vect(:,j,ih)=d
         enddo
      enddo
      cell_index=1
      macro_nlistsize=natom
      macro_atom_nlist=0
      macro_atom_nlist(1,1:natom)=(/(ih,ih=1,natom)/)
   end subroutine setup_hamiltonian

   subroutine global_raw_energy(state,raw)
      real(dblprec), intent(in) :: state(3,natom,nensemble)
      real(dblprec), intent(out) :: raw
      integer :: ii,kk

      call effective_field(natom,nensemble,1,natom,state,mmom,external_field,time_external_field, &
         beff,beff1,beff2,effective_energy,nmacro,cell_index,emomM_macro,macro_nlistsize, &
         1,1,1,1,measure_energy=.true.,term_fields=terms)
      raw=0.0_dblprec
      do kk=1,nensemble
         do ii=1,natom
            raw=raw+canonical_field_energy(HAM_TERM_EXCHANGE,state(:,ii,kk),terms(:,HAM_TERM_EXCHANGE,ii,kk))
            raw=raw+canonical_field_energy(HAM_TERM_DMI,state(:,ii,kk),terms(:,HAM_TERM_DMI,ii,kk))
            if (ham_inp%do_anisotropy==1) raw=raw+canonical_onsite_energy(ii,state(:,ii,kk),ham_inp%mult_axis)
         enddo
      enddo
      call check(abs(effective_energy*mry/mub-raw)<=tolerance, &
         'effective_field returned energy matches canonical term total')
   end subroutine global_raw_energy

   subroutine mc_delta(do_dm_value,do_anis_value,de)
      integer, intent(in) :: do_dm_value,do_anis_value
      real(dblprec), intent(out) :: de

      call calculate_energy(natom,nensemble,natom,1,do_dm_value,0,0,0,0,0,0, &
         emomM,emom,mmom,iflip,newmom,extfield,de,1,'N',0,nmacro,max_macro, &
         cell_index,macro_nlistsize,macro_atom_nlist,emomM_macro,icell, &
         macro_mag_trial,macro_trial,'N',do_anis_value,0)
   end subroutine mc_delta

   subroutine clear_hamiltonian()
      deallocate(ham%aHam,ham%nlistsize,ham%nlist,ham%ncoup,ham%dmlistsize,ham%dmlist, &
         ham%dm_vect,ham%taniso,ham%taniso_diff,ham%sb,ham%sb_diff,ham%kaniso, &
         ham%kaniso_diff,ham%eaniso,ham%eaniso_diff)
   end subroutine clear_hamiltonian

   subroutine check(condition,label)
      logical, intent(in) :: condition
      character(len=*), intent(in) :: label
      if (.not.condition) then
         write(*,'(a)') 'FAILED: '//trim(label)
         failures=failures+1
      endif
   end subroutine check

end program test_cpu_ham11_energy
