! CPU-HAM-03B: portable persistent scalar-J sparse backend parity.
program test_sparse_backend
   use Parameters, only : dblprec
   use HamiltonianData, only : ham
   use HamiltonianActions, only : effective_field, setup_sparse_backend, &
      cleanup_sparse_backend, sparse_backend_can_apply, sparse_backend_get_stats
   use InputData, only : ham_inp, do_sparse
   implicit none

   integer :: failures

   failures=0
   call run_case(7,3,'Nd')
   call run_case(10,2,'Fe')
   call check_fallback()
   call cleanup_sparse_backend()

   if (failures /= 0) error stop 'CPU-HAM-03B sparse backend parity failed'
   write(*,'(a)') 'CPU-HAM-03B sparse backend parity passed'

contains

   subroutine run_case(natom,nensemble,label)
      integer, intent(in) :: natom,nensemble
      character(len=*), intent(in) :: label
      integer, parameter :: max_neigh=4, nmacro=1
      integer :: i,j,k,nnz
      integer :: cell_index(natom),macro_nlistsize(nmacro)
      real(dblprec) :: emomM(3,natom,nensemble),mmom(natom,nensemble)
      real(dblprec) :: external_field(3,natom,nensemble),time_external_field(3,natom,nensemble)
      real(dblprec) :: emomM_macro(3,nmacro,nensemble)
      real(dblprec) :: beff_sparse(3,natom,nensemble),beff_direct(3,natom,nensemble)
      real(dblprec) :: beff1_sparse(3,natom,nensemble),beff1_direct(3,natom,nensemble)
      real(dblprec) :: beff2_sparse(3,natom,nensemble),beff2_direct(3,natom,nensemble)
      real(dblprec) :: energy_sparse,energy_direct
      real(dblprec) :: setup_seconds,pack_seconds,apply_seconds
      integer(kind=8) :: sparse_nnz

      call setup_case(natom,max_neigh)
      do k=1,nensemble
         do i=1,natom
            mmom(i,k)=1.0_dblprec+0.01_dblprec*real(mod(i+k,5),dblprec)
            do j=1,3
               emomM(j,i,k)=sin(real(17*i+5*j+3*k,dblprec))
            end do
         end do
      end do
      external_field=0.0_dblprec
      time_external_field=0.0_dblprec
      emomM_macro=0.0_dblprec
      cell_index=1
      macro_nlistsize=natom
      nnz=sum(ham%nlistsize)

      do_sparse='Y'
      call setup_sparse_backend(natom,nensemble,0,'N')
      call check(sparse_backend_can_apply(natom,nensemble,1,natom), &
         trim(label)//' sparse backend is active')
      call effective_field(natom,nensemble,1,natom,emomM,mmom,external_field, &
         time_external_field,beff_sparse,beff1_sparse,beff2_sparse,energy_sparse, &
         nmacro,cell_index,emomM_macro,macro_nlistsize,1,1,1,1,measure_energy=.true.)

      do_sparse='N'
      call effective_field(natom,nensemble,1,natom,emomM,mmom,external_field, &
         time_external_field,beff_direct,beff1_direct,beff2_direct,energy_direct, &
         nmacro,cell_index,emomM_macro,macro_nlistsize,1,1,1,1,measure_energy=.true.)

      call check(maxval(abs(beff_sparse-beff_direct)) <= 1.0d-13, &
         trim(label)//' total field matches DIRECT')
      call check(maxval(abs(beff1_sparse-beff1_direct)) <= 1.0d-13, &
         trim(label)//' internal field matches DIRECT')
      call check(maxval(abs(beff2_sparse-beff2_direct)) <= 1.0d-13, &
         trim(label)//' external field matches DIRECT')
      call check(abs(energy_sparse-energy_direct) <= 1.0d-13, &
         trim(label)//' field-derived energy matches DIRECT')

      call sparse_backend_get_stats(setup_seconds,pack_seconds,apply_seconds,sparse_nnz)
      write(*,'(a,a,a,3(es12.4,1x),a,i0)') 'CPU-HAM-03B ',trim(label), &
         ' setup/pack/apply seconds=',setup_seconds,pack_seconds,apply_seconds, &
         ' nnz=',sparse_nnz
      call check(sparse_nnz == int(nnz,kind=8),trim(label)//' CSR has exact directed nnz')
      call check(setup_seconds >= 0.0_dblprec,trim(label)//' setup cost is recorded')
      call check(pack_seconds >= 0.0_dblprec,trim(label)//' packing cost is recorded')
      call check(apply_seconds >= 0.0_dblprec,trim(label)//' sparse apply cost is recorded')

      ! Repeated setup must replace, not append to, the persistent structure.
      do_sparse='Y'
      call setup_sparse_backend(natom,nensemble,0,'N')
      call check(sparse_backend_can_apply(natom,nensemble,1,natom), &
         trim(label)//' repeated setup remains active')
      call cleanup_sparse_backend()
      call check(.not. sparse_backend_can_apply(natom,nensemble,1,natom), &
         trim(label)//' cleanup releases sparse state')
      call clear_case()
   end subroutine run_case

   subroutine setup_case(natom,max_neigh)
      integer, intent(in) :: natom,max_neigh
      integer :: i,j,n

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

      allocate(ham%aHam(natom),ham%nlistsize(natom),ham%nlist(max_neigh,natom), &
         ham%ncoup(max_neigh,natom,1),ham%ncoupD(max_neigh,natom,1))
      ham%aHam=(/(i,i=1,natom)/)
      ham%nlistsize=0
      ham%nlist=0
      ham%ncoup=0.0_dblprec
      ham%ncoupD=0.0_dblprec
      do i=1,natom
         n=mod(i,4)
         ham%nlistsize(i)=n
         do j=1,n
            ham%nlist(j,i)=mod(i+2*j-2,natom)+1
            ham%ncoup(j,i,1)=0.07_dblprec*real(i-j,dblprec)
            ham%ncoupD(j,i,1)=0.03_dblprec*real(i+j,dblprec)
         end do
      end do
   end subroutine setup_case

   subroutine check_fallback()
      integer, parameter :: natom=5,nensemble=2,max_neigh=4,nmacro=1
      integer :: cell_index(natom),macro_nlistsize(nmacro)
      real(dblprec) :: emomM(3,natom,nensemble),mmom(natom,nensemble)
      real(dblprec) :: ext(3,natom,nensemble),text(3,natom,nensemble)
      real(dblprec) :: emomM_macro(3,nmacro,nensemble)
      real(dblprec) :: beff_sparse(3,natom,nensemble),beff_direct(3,natom,nensemble)
      real(dblprec) :: beff1(3,natom,nensemble),beff2(3,natom,nensemble),energy
      integer :: i,j,k

      call setup_case(natom,max_neigh)
      do k=1,nensemble
         do i=1,natom
            mmom(i,k)=1.0_dblprec
            emomM(:,i,k)=(/0.2_dblprec*i, -0.1_dblprec*k, 0.03_dblprec*(i+k)/)
         end do
      end do
      ext=0.0_dblprec
      text=0.0_dblprec
      emomM_macro=0.0_dblprec
      cell_index=1
      macro_nlistsize=natom

      ham_inp%exc_inter='Y'
      do_sparse='Y'
      call setup_sparse_backend(natom,nensemble,0,'N')
      call check(.not. sparse_backend_can_apply(natom,nensemble,1,natom), &
         'rescaled exchange explicitly falls back from sparse')
      call effective_field(natom,nensemble,1,natom,emomM,mmom,ext,text, &
         beff_sparse,beff1,beff2,energy,nmacro,cell_index,emomM_macro, &
         macro_nlistsize,1,1,1,1,measure_energy=.false.)
      do_sparse='N'
      call effective_field(natom,nensemble,1,natom,emomM,mmom,ext,text, &
         beff_direct,beff1,beff2,energy,nmacro,cell_index,emomM_macro, &
         macro_nlistsize,1,1,1,1,measure_energy=.false.)
      call check(maxval(abs(beff_sparse-beff_direct)) <= 1.0d-13, &
         'rescaled exchange fallback matches DIRECT')
      call cleanup_sparse_backend()
      call clear_case()
   end subroutine check_fallback

   subroutine clear_case()
      if (allocated(ham%aHam)) deallocate(ham%aHam)
      if (allocated(ham%nlistsize)) deallocate(ham%nlistsize)
      if (allocated(ham%nlist)) deallocate(ham%nlist)
      if (allocated(ham%ncoup)) deallocate(ham%ncoup)
      if (allocated(ham%ncoupD)) deallocate(ham%ncoupD)
   end subroutine clear_case

   subroutine check(condition,label)
      logical, intent(in) :: condition
      character(len=*), intent(in) :: label
      if (.not. condition) then
         write(*,'(a)') 'FAILED: '//trim(label)
         failures=failures+1
      endif
   end subroutine check

end program test_sparse_backend
