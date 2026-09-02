! CPU-HAM-06: periodic J+D parity and DMI negative controls.
program test_cpu_ham06_j_plus_d

   use Parameters, only : dblprec
   use Constants, only : mub,mry
   use HamiltonianData, only : ham
   use HamiltonianActions, only : effective_field, setup_convolution_backend, &
      setup_cpu_hamiltonian_backend, cleanup_convolution_backend, convolution_backend_can_apply
   use InputData, only : ham_inp, cpu_ham_backend, do_sparse, do_convolution
   use ReducedStencil
   use CPUConvolution

   implicit none

   integer :: failures

   failures=0
   call run_case('2D skyrmion-shaped J+D',1,5,4,1,4,failures)
   call run_case('3D chiral J+D',2,4,3,3,6,failures)
   if (failures /= 0) error stop 'CPU-HAM-06 J+D tests failed'
   write(*,'(a)') 'CPU-HAM-06 J+D tests passed'

contains

   subroutine run_case(label,na,n1,n2,n3,z,failures)
      character(len=*), intent(in) :: label
      integer, intent(in) :: na,n1,n2,n3,z
      integer, intent(inout) :: failures

      integer :: natom,ensembles,maxz,maxd,atom,b,slot,axis,ensemble
      integer :: cell(3),wrapped(3),source_cell(3),source
      integer, allocatable :: ham_index(:),nlistsize(:),nlist(:,:),dmlistsize(:),dmlist(:,:)
      real(dblprec), allocatable :: ncoup(:,:,:),dm_vect(:,:,:),bad_dm(:,:,:)
      real(dblprec), allocatable :: spin(:,:,:),mmom(:,:),ext(:,:,:),text(:,:,:)
      real(dblprec), allocatable :: beff(:,:,:),beff1(:,:,:),beff2(:,:,:)
      real(dblprec), allocatable :: direct_field(:,:,:),reduced_field(:,:,:),wrong_field(:,:,:)
      real(dblprec), allocatable :: convolution_field(:,:,:)
      real(dblprec), allocatable :: emomM_macro(:,:,:)
      integer, allocatable :: cell_index(:),macro_nlistsize(:)
      type(reduced_stencil_t) :: stencil,wrong_stencil,transposed_stencil
      type(cpu_convolution_t) :: convolution
      real(dblprec) :: direct_energy,reduced_energy,convolution_energy
      real(dblprec) :: pair_energy
      real(dblprec) :: x,y,r,theta,phi,radius,pi
      character(len=256) :: diagnostic
      logical :: ok

      natom=na*n1*n2*n3
      ensembles=2
      maxz=z
      maxd=z
      allocate(ham_index(natom),nlistsize(natom),nlist(maxz,natom), &
         dmlistsize(natom),dmlist(maxd,natom),ncoup(maxz,na,1),dm_vect(3,maxd,na))
      ham_index=0
      nlistsize=0
      nlist=0
      dmlistsize=0
      dmlist=0
      ncoup=0.0_dblprec
      dm_vect=0.0_dblprec

      do atom=1,natom
         call atom_to_cell_basis(atom,na,n1,n2,n3,cell,b)
         ham_index(atom)=b
         nlistsize(atom)=z
         dmlistsize(atom)=z
         do slot=1,z
            call fixture_delta(slot,n1,n2,n3,wrapped)
            call wrap_reduced_cell(cell,wrapped,n1,n2,n3,source_cell)
            source=cell_basis_to_atom(source_cell(1),source_cell(2),source_cell(3),b,na,n1,n2,n3)
            nlist(slot,atom)=source
            dmlist(slot,atom)=source
         end do
      end do
      do b=1,na
         do slot=1,z
            ncoup(slot,b,1)=0.35_dblprec+0.04_dblprec*real(slot+b,dblprec)
            call fixture_dmi(slot,dm_vect(:,slot,b))
         end do
      end do

      call clear_fixture()
      allocate(ham%aHam(natom),ham%nlistsize(natom),ham%nlist(maxz,natom), &
         ham%ncoup(maxz,na,1),ham%dmlistsize(natom),ham%dmlist(maxd,natom), &
         ham%dm_vect(3,maxd,na))
      ham%aHam=ham_index
      ham%nlistsize=nlistsize
      ham%nlist=nlist
      ham%ncoup=ncoup
      ham%dmlistsize=dmlistsize
      ham%dmlist=dmlist
      ham%dm_vect=dm_vect

      call build_reduced_stencil(natom,na,n1,n2,n3,'P','P','P','Y',0,na,ham_index, &
         nlistsize,nlist,ncoup,stencil,ok,diagnostic,dmlistsize=dmlistsize, &
         dmlist=dmlist,dm_vect=dm_vect)
      call check(ok,trim(label)//': combined stencil accepted',failures)
      call check(allocated(stencil%dmi_record_start),trim(label)//': DMI records retained',failures)
      if (.not.ok) then
         call clear_fixture()
         deallocate(ham_index,nlistsize,nlist,dmlistsize,dmlist,ncoup,dm_vect)
         return
      endif

      allocate(spin(3,natom,ensembles),mmom(natom,ensembles),ext(3,natom,ensembles), &
         text(3,natom,ensembles),beff(3,natom,ensembles),beff1(3,natom,ensembles), &
         beff2(3,natom,ensembles),direct_field(3,natom,ensembles), &
         reduced_field(3,natom,ensembles),wrong_field(3,natom,ensembles), &
         convolution_field(3,natom,ensembles), &
         emomM_macro(3,1,ensembles),cell_index(natom),macro_nlistsize(1))
      do ensemble=1,ensembles
         do atom=1,natom
            call atom_to_cell_basis(atom,na,n1,n2,n3,cell,b)
            do axis=1,3
               spin(axis,atom,ensemble)=0.31_dblprec*sin(real(axis+3*cell(1)+ &
                  5*cell(2)+7*cell(3)+2*b+11*ensemble,dblprec)) + &
                  0.07_dblprec*cos(real(2*axis+cell(1)+3*cell(2)+b,dblprec))
            end do
         end do
      end do
      mmom=1.0_dblprec
      ext=0.0_dblprec
      text=0.0_dblprec
      emomM_macro=0.0_dblprec
      cell_index=1
      macro_nlistsize(1)=natom
      pi=acos(-1.0_dblprec)
      radius=0.5_dblprec*real(min(n1,n2),dblprec)

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
      do_sparse='N'
      do_convolution='N'

      if (n3 == 1 .and. na == 1) then
         ! Radial Néel skyrmion-like texture for the 2D parity fixture.
         do ensemble=1,ensembles
            do atom=1,natom
               call atom_to_cell_basis(atom,na,n1,n2,n3,cell,b)
               x=real(cell(1),dblprec)-0.5_dblprec*real(n1-1,dblprec)
               y=real(cell(2),dblprec)-0.5_dblprec*real(n2-1,dblprec)
               r=sqrt(x*x+y*y)
               theta=pi*max(0.0_dblprec,1.0_dblprec-r/radius)
               phi=atan2(y,x)+0.13_dblprec*real(ensemble-1,dblprec)
               spin(:,atom,ensemble)=(/sin(theta)*cos(phi),sin(theta)*sin(phi),cos(theta)/)
            end do
         end do
      endif

      ! Canonical production DIRECT is the physics oracle.
      call clear_reduced_stencil(ham%reduced_stencil)
      call effective_field(natom,ensembles,1,natom,spin,mmom,ext,text,beff,beff1,beff2, &
         direct_energy,1,cell_index,emomM_macro,macro_nlistsize,na,n1,n2,n3, &
         measure_energy=.true.)
      direct_field=beff1
      pair_energy=-0.5_dblprec*sum(spin*direct_field)
      call check(abs(direct_energy-pair_energy*mub/mry) < 2.0d-12, &
         trim(label)//': DIRECT uses field-derived pair energy',failures)

      ! The production reduced target path applies J and D separately but from
      ! one validated periodic representation.
      ham%reduced_stencil=stencil
      call effective_field(natom,ensembles,1,natom,spin,mmom,ext,text,beff,beff1,beff2, &
         reduced_energy,1,cell_index,emomM_macro,macro_nlistsize,na,n1,n2,n3, &
         measure_energy=.true.)
      reduced_field=beff1
      call compare_fields(direct_field,reduced_field,trim(label)//': reduced field parity',failures)
      call check(abs(reduced_energy-direct_energy) < 3.0d-12, &
         trim(label)//': reduced energy parity',failures)

      ! Reverse the DMI sign in the candidate representation.  It remains a
      ! valid reciprocal stencil, but must not agree with the canonical oracle.
      wrong_stencil=stencil
      do slot=1,size(wrong_stencil%dmi_record)
         wrong_stencil%dmi_record(slot)%d=-wrong_stencil%dmi_record(slot)%d
      end do
      ham%reduced_stencil=wrong_stencil
      call effective_field(natom,ensembles,1,natom,spin,mmom,ext,text,beff,beff1,beff2, &
         reduced_energy,1,cell_index,emomM_macro,macro_nlistsize,na,n1,n2,n3, &
         measure_energy=.false.)
      wrong_field=beff1
      call check(maxval(abs(wrong_field-direct_field)) > 1.0d-8, &
         trim(label)//': reversed DMI sign negative control fails parity',failures)
      call clear_reduced_stencil(wrong_stencil)

      ! A component-transposed DMI operator is also a valid reciprocal list,
      ! but it must not reproduce the canonical spin/spatial contraction.
      transposed_stencil=stencil
      do slot=1,size(transposed_stencil%dmi_record)
         transposed_stencil%dmi_record(slot)%d=(/transposed_stencil%dmi_record(slot)%d(2), &
            transposed_stencil%dmi_record(slot)%d(1),transposed_stencil%dmi_record(slot)%d(3)/)
      end do
      ham%reduced_stencil=transposed_stencil
      call effective_field(natom,ensembles,1,natom,spin,mmom,ext,text,beff,beff1,beff2, &
         reduced_energy,1,cell_index,emomM_macro,macro_nlistsize,na,n1,n2,n3, &
         measure_energy=.false.)
      wrong_field=beff1
      call check(maxval(abs(wrong_field-direct_field)) > 1.0d-8, &
         trim(label)//': transposed DMI operator negative control fails parity',failures)
      call clear_reduced_stencil(transposed_stencil)

      ! Restore the accepted stencil before exercising convolution.
      ham%reduced_stencil=stencil
      ok=cpu_convolution_init(convolution,stencil,ensembles,diagnostic)
      call check(ok,trim(label)//': J+D convolution object initialized',failures)
      if (ok) ok=cpu_convolution_build_kernel(convolution,stencil,diagnostic)
      call check(ok,trim(label)//': J+D spectral kernels built',failures)
      if (ok) then
         call check(convolution%kernel_batches == 4*na*na, &
            trim(label)//': convolution stores J,Dx,Dy,Dz basis kernels',failures)
         ok=cpu_convolution_apply(convolution,spin,convolution_field,diagnostic)
         call check(ok,trim(label)//': direct convolution apply succeeds',failures)
         if (ok) then
            call compare_fields(direct_field,convolution_field, &
               trim(label)//': direct convolution field parity',failures)
            pair_energy=-0.5_dblprec*sum(spin*convolution_field)
            call check(abs(pair_energy-(-0.5_dblprec*sum(spin*direct_field))) < 5.0d-11, &
               trim(label)//': convolution field-derived pair energy',failures)
         endif
      endif
      call cpu_convolution_clear(convolution)
      do_convolution='Y'
      cpu_ham_backend='convolution'
      call setup_cpu_hamiltonian_backend(natom,ensembles,0,'N',na,n1,n2,n3,'P','P','P','Y',na)
      call check(convolution_backend_can_apply(natom,ensembles,1,natom), &
         trim(label)//': J+D convolution backend is active',failures)
      if (convolution_backend_can_apply(natom,ensembles,1,natom)) then
         call effective_field(natom,ensembles,1,natom,spin,mmom,ext,text,beff,beff1,beff2, &
            convolution_energy,1,cell_index,emomM_macro,macro_nlistsize,na,n1,n2,n3, &
            measure_energy=.true.)
         call compare_fields(direct_field,beff1,trim(label)//': convolution field parity',failures)
         call check(abs(convolution_energy-direct_energy) < 4.0d-12, &
            trim(label)//': convolution energy parity',failures)
      endif
      call cleanup_convolution_backend()
      do_convolution='N'
      cpu_ham_backend='direct'

      deallocate(spin,mmom,ext,text,beff,beff1,beff2,direct_field,reduced_field,wrong_field, &
         convolution_field, &
         emomM_macro,cell_index,macro_nlistsize)
      call clear_reduced_stencil(stencil)
      call clear_fixture()

      ! Changing one reverse DMI endpoint to the transpose/wrong-antisymmetry
      ! convention must be rejected before any optimized apply is enabled.
      bad_dm=dm_vect
      bad_dm(:,2,:)=bad_dm(:,1,:)
      call build_reduced_stencil(natom,na,n1,n2,n3,'P','P','P','Y',0,na,ham_index, &
         nlistsize,nlist,ncoup,stencil,ok,diagnostic,dmlistsize=dmlistsize, &
         dmlist=dmlist,dm_vect=bad_dm)
      call check(.not.ok .and. index(diagnostic,'reciprocal/antisymmetric')>0, &
         trim(label)//': wrong antisymmetric operator is rejected',failures)
      if (allocated(bad_dm)) deallocate(bad_dm)
      call clear_reduced_stencil(stencil)
      deallocate(ham_index,nlistsize,nlist,dmlistsize,dmlist,ncoup,dm_vect)

   end subroutine run_case

   subroutine fixture_delta(slot,n1,n2,n3,delta)
      integer, intent(in) :: slot,n1,n2,n3
      integer, intent(out) :: delta(3)
      delta=0
      if (slot == 1) delta(1)=1
      if (slot == 2) delta(1)=n1-1
      if (slot == 3) delta(2)=1
      if (slot == 4) delta(2)=n2-1
      if (slot == 5) delta(3)=1
      if (slot == 6) delta(3)=n3-1
   end subroutine fixture_delta

   subroutine fixture_dmi(slot,vector)
      integer, intent(in) :: slot
      real(dblprec), intent(out) :: vector(3)
      vector=0.0_dblprec
      select case(slot)
      case(1); vector(3)=0.23_dblprec
      case(2); vector(3)=-0.23_dblprec
      case(3); vector(1)=0.17_dblprec
      case(4); vector(1)=-0.17_dblprec
      case(5); vector(2)=0.11_dblprec
      case(6); vector(2)=-0.11_dblprec
      end select
   end subroutine fixture_dmi

   subroutine compare_fields(expected,actual,label,count)
      real(dblprec), intent(in) :: expected(:,:,:),actual(:,:,:)
      character(len=*), intent(in) :: label
      integer, intent(inout) :: count
      real(dblprec) :: max_abs,scale
      max_abs=maxval(abs(actual-expected))
      scale=max(1.0_dblprec,maxval(abs(expected)))
      write(*,'(2x,a,2(es12.4,1x))') trim(label)//' max_abs/rel=',max_abs,max_abs/scale
      call check(max_abs < 5.0d-11 .and. max_abs/scale < 5.0d-11,label,count)
   end subroutine compare_fields

   subroutine check(condition,label,count)
      logical, intent(in) :: condition
      character(len=*), intent(in) :: label
      integer, intent(inout) :: count
      if (.not.condition) then
         write(*,'(a)') 'FAILED: '//trim(label)
         count=count+1
      endif
   end subroutine check

   subroutine clear_fixture()
      call clear_reduced_stencil(ham%reduced_stencil)
      if (allocated(ham%aHam)) deallocate(ham%aHam)
      if (allocated(ham%nlistsize)) deallocate(ham%nlistsize)
      if (allocated(ham%nlist)) deallocate(ham%nlist)
      if (allocated(ham%ncoup)) deallocate(ham%ncoup)
      if (allocated(ham%dmlistsize)) deallocate(ham%dmlistsize)
      if (allocated(ham%dmlist)) deallocate(ham%dmlist)
      if (allocated(ham%dm_vect)) deallocate(ham%dm_vect)
   end subroutine clear_fixture

end program test_cpu_ham06_j_plus_d
