#ifdef __INTEL_MKL__
   include "mkl_spblas.f90"
#endif
module Sparse
   !
   use Parameters
   use ErrorHandling
   !
#ifdef __INTEL_MKL__
   use mkl_spblas
#endif
   implicit none
   !
   real(dblprec), dimension(:), pointer :: Hscalar
   real(dblprec), dimension(:,:,:), pointer :: Hblock
   integer, dimension(:), allocatable :: Hcolumns
   integer, dimension(:), allocatable :: HpointerB
   integer, dimension(:), allocatable :: HpointerE
   !
   !integer, dimension(:,:), allocatable :: max_cut_list
   !
   integer :: maxSdim, totSdim, i_stat
   integer :: Hblocksize = 3
   !
#ifdef __INTEL_MKL__
    type(SPARSE_MATRIX_T) :: H_sparse
    type(MATRIX_DESCR) :: H_desc
#endif

   private

   public ::  setupSparseScalar,setupSparseBlock, setupSparseTensor
   public ::  effective_field_SparseScalar, effective_field_SparseBlock

contains

   subroutine setupSparseScalar(Natom,nHam,conf_num,max_no_neigh,nlist,nlistsize,ncoup, aham)
      !
      implicit none
      !
      integer, intent(in) :: Natom  !< Number of atoms in system
      integer, intent(in) :: nHam  !< Number of atoms in Hamiltonian
      integer, intent(in) :: conf_num !< Number of LSF configurations
      integer, intent(in) :: max_no_neigh !< Calculated maximum of neighbours for exchange
      integer, dimension(max_no_neigh,Natom), intent(in) :: nlist  !< Neighbour list for Heisenberg exchange couplings
      integer, dimension(nHam), intent(in) :: nlistsize !< Size of neighbour list for Heisenberg exchange couplings
      real(dblprec), dimension(max_no_neigh,nHam,conf_num), intent(in) :: ncoup !< Heisenberg exchange couplings
      integer, dimension(nHam), intent(in) :: aham !< Hamiltonian look-up table
      !
      integer :: i_atom, j_atom, j, nelem, ielem, i_ham
      !
      nelem=max_no_neigh*Natom
      call allocateSparse(Natom,nelem,1)
      !
      Hscalar=0.0_dblprec;Hcolumns=0;HpointerB=0;HpointerE=0
      !
      ielem=0
      do i_atom=1,Natom
         i_ham=aham(i_atom)
         HpointerB(i_atom)=ielem+1-1
         do j=1,nlistsize(i_ham)
            j_atom=nlist(j,i_atom)
            ielem=ielem+1
            !Hscalar(ielem)=ncoup(j_atom,i_atom,1)
            Hscalar(ielem)=ncoup(j,i_ham,1)
            !Hcolumns(ielem)=j_atom
            Hcolumns(ielem)=j_atom-1
         end do
         HpointerE(i_atom)=ielem+1-1
      end do
      !
#ifdef __INTEL_MKL__
      i_stat=MKL_SPARSE_D_CREATE_CSR(H_sparse,SPARSE_INDEX_BASE_ZERO,Natom,Natom,HpointerB,HpointerE,Hcolumns,Hscalar)
      H_desc%type=SPARSE_MATRIX_TYPE_SYMMETRIC
      !H_desc%type=SPARSE_MATRIX_TYPE_GENERAL
      H_desc%mode=SPARSE_FILL_MODE_LOWER
      H_desc%diag=SPARSE_DIAG_NON_UNIT
      print *,'Sparse matrix creation',i_stat
#endif
      !print *,Hscalar(1:ielem-1)
   end subroutine setupSparseScalar

   subroutine setupSparseTensor(Natom,nHam,conf_num,max_no_neigh,nlist,nlistsize,j_tens, aham)
      !
      implicit none
      !
      integer, intent(in) :: Natom  !< Number of atoms in system
      integer, intent(in) :: nHam   !< Number of atoms in Hamiltonian
      !integer, intent(in) :: Mensemble !< Number of ensembles
      integer, intent(in) :: conf_num !< Number of LSF configurations
      integer, intent(in) :: max_no_neigh !< Calculated maximum of neighbours for exchange
      integer, dimension(max_no_neigh,Natom), intent(in) :: nlist  !< Neighbour list for Heisenberg exchange couplings
      integer, dimension(nHam), intent(in) :: nlistsize !< Size of neighbour list for Heisenberg exchange couplings
      real(dblprec), dimension(3,3,max_no_neigh,nHam), intent(in) :: j_tens !< Heisenberg exchange couplings
      integer, dimension(nHam), intent(in) :: aham !< Hamiltonian look-up table
      !real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: external_field !< External magnetic field
      !
      integer :: i_atom,j_atom,j,nelem, ielem, i_ham
      !
      nelem=(max_no_neigh)*Natom
      call allocateSparse(Natom,nelem,1)
      !
      HBlock=0.0_dblprec;Hcolumns=0;HpointerB=0;HpointerE=0
      !
      ielem=0
      !print *,'In setupBlock'
      ielem=0
      do i_atom=1,Natom
         i_ham=aham(i_atom)
         HpointerB(i_atom)=ielem
         do j=1,nlistsize(i_ham)
            j_atom=nlist(j,i_atom)
            ielem=ielem+1
            Hblock(:,:,ielem)=transpose(j_tens(:,:,j,i_ham))
            Hcolumns(ielem)=j_atom-1
         end do
         HpointerE(i_atom)=ielem
      end do
      !!
   end subroutine setupSparseTensor

   subroutine setupSparseBlock(Natom,nHam,conf_num,max_no_neigh,nlist,nlistsize,ncoup,max_no_dmneigh,dmlistsize,dmlist,dm_vect, aham)
      !
      implicit none
      !
      integer, intent(in) :: Natom  !< Number of atoms in system
      integer, intent(in) :: nHam   !< Number of atoms in Hamiltonian
      !integer, intent(in) :: Mensemble !< Number of ensembles
      integer, intent(in) :: conf_num !< Number of LSF configurations
      integer, intent(in) :: max_no_neigh !< Calculated maximum of neighbours for exchange
      integer, dimension(max_no_neigh,Natom), intent(in) :: nlist  !< Neighbour list for Heisenberg exchange couplings
      integer, dimension(nHam), intent(in) :: nlistsize !< Size of neighbour list for Heisenberg exchange couplings
      real(dblprec), dimension(max_no_neigh,nHam,conf_num), intent(in) :: ncoup !< Heisenberg exchange couplings
      integer, intent(in) :: max_no_dmneigh !< Calculated number of neighbours with DM interactions
      integer, dimension(Natom),intent(in) :: dmlistsize !< Size of neighbour list for DM
      integer, dimension(max_no_dmneigh,Natom), intent(in) :: dmlist   !< List of neighbours for DM
      real(dblprec), dimension(3,max_no_dmneigh,Natom), intent(in) :: dm_vect !< Dzyaloshinskii-Moriya exchange vector
      integer, dimension(nHam), intent(in) :: aham !< Hamiltonian look-up table
      !real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: external_field !< External magnetic field
      !
      integer :: i_atom,j_atom,j,nelem, ielem, i_ham
      !
      nelem=(max_no_neigh)*Natom
      call allocateSparse(Natom,nelem,1)
      !
      HBlock=0.0_dblprec;Hcolumns=0;HpointerB=0;HpointerE=0
      !
      ielem=0
      if(max_no_neigh.ne.max_no_dmneigh) then
         write(*,'(1x,a)') 'Warning: J and DMI lists do not agree. Tensor setup shaky.'
      end if
      !print *,'In setupBlock'
      ielem=0
      do i_atom=1,Natom
         i_ham=aham(i_atom)
         HpointerB(i_atom)=ielem
         do j=1,nlistsize(i_ham)
            j_atom=nlist(j,i_atom)
            ielem=ielem+1
            !!Hscalar(ielem)=ncoup(j_atom,i_atom,1)
            !Hscalar(ielem)=ncoup(j,i_ham,1)
            !!Hcolumns(ielem)=j_atom
            ! Scalar (diagonal) Heisenberg
            Hblock(1,1,ielem)=ncoup(j,i_ham,1)
            Hblock(2,2,ielem)=ncoup(j,i_ham,1)
            Hblock(3,3,ielem)=ncoup(j,i_ham,1)
            ! Anti-symmetric (off-diagonal) DMI
            Hblock(1,2,ielem)=-dm_vect(3,j,i_ham)
            Hblock(1,3,ielem)= dm_vect(2,j,i_ham)
            Hblock(2,3,ielem)=-dm_vect(1,j,i_ham)
            Hblock(2,1,ielem)= dm_vect(3,j,i_ham)
            Hblock(3,1,ielem)=-dm_vect(2,j,i_ham)
            Hblock(3,2,ielem)= dm_vect(1,j,i_ham)
            Hcolumns(ielem)=j_atom-1
         end do
         HpointerE(i_atom)=ielem
      end do
      !do i_atom=1,Natom
      !   i_ham=aham(i_atom)
      !   HpointerB(i_atom)=ielem
      !   !! on-site
      !   !ielem=ielem+1
      !   !Hblock(1,1,ielem)=0.0_dblprec !external_field(1,i_atom,1)
      !   !Hblock(2,2,ielem)=0.0_dblprec !external_field(2,i_atom,1)
      !   !Hblock(3,3,ielem)=0.0_dblprec !external_field(3,i_atom,1)
      !   !Hcolumns(ielem)=i_atom-1

      !   do j=1,nlistsize(i_ham)
      !      ielem=ielem+1
      !      ! Scalar (diagonal) Heisenberg
      !      Hblock(1,1,ielem)=ncoup(j,i_ham,1)
      !      Hblock(2,2,ielem)=ncoup(j,i_ham,1)
      !      Hblock(3,3,ielem)=ncoup(j,i_ham,1)
      !      ! Anti-symmetric (off-diagonal) DMI
      !      Hblock(1,2,ielem)=-dm_vect(3,j,i_ham)
      !      Hblock(1,3,ielem)= dm_vect(2,j,i_ham)
      !      Hblock(2,3,ielem)=-dm_vect(1,j,i_ham)
      !      Hblock(2,1,ielem)= dm_vect(3,j,i_ham)
      !      Hblock(3,1,ielem)=-dm_vect(2,j,i_ham)
      !      Hblock(3,2,ielem)= dm_vect(1,j,i_ham)
      !      !
      !      j_atom=nlist(j,i_atom)
      !      Hcolumns(ielem)=j_atom-1
      !   end do
      !   HpointerE(i_atom)=ielem
      !end do
      !!Hcolumns=Hcolumns+1
      !!HpointerE=HpointerE+1
      !!HpointerB=HpointerB+1
      !! For zero-based indexing
      !!matdescra=(/'G','U','N','C','0','0'/)
      !! For one-based indexing
      !!matdescra=(/'G','U','N','F','0','0'/)
      !Hcolumns=Hcolumns+1
      !HpointerB=HpointerB+1
      !HpointerE=HpointerE+1
      !!
   end subroutine setupSparseBlock

   subroutine effective_field_SparseBlock(Natom,Mensemble,emomM,beff)
      !
      implicit none
      !
      !include 'mkl.fi'
      !
      integer, intent(in) :: Natom  !< Number of atoms in system
      integer, intent(in) :: Mensemble !< Number of ensembles 
      real(dblprec), dimension(3,Natom,Mensemble),intent(in) :: emomM  !< Current magnetic moment vector
      real(dblprec), dimension(3,Natom,Mensemble), intent(out) :: beff !< Total effective field from application of Hamiltonian
      !
      integer :: i, j, k, m, n
      character, dimension(6) :: matdescra

      real(dblprec) :: mone = 1.0_dblprec
      real(dblprec) :: zero = 0.0_dblprec

      !
      ! zero based
      matdescra=(/'G','U','N','C','0','0'/)
      !
      do k=1,Mensemble
         !call mkl_dbsrmv('N',Natom,Natom,3, mone, matdescra, HBlock, Hcolumns, HpointerB, HpointerE, emomM_t(:,:), zero, beff_t(:,:))
#ifdef __INTEL_MKL__
        call mkl_dbsrmv('N',Natom,Natom,3, mone, matdescra, HBlock, Hcolumns, HpointerB, HpointerE, emomM(:,:,k), zero, beff(:,:,k))
#else
         !call ErrorHandling_ERROR('Sparse operations only implemented for mkl')
         !$omp parallel do default(shared) private(i,j)
         do i=1,Natom
            !do j=HpointerB(i)+1,HpointerB(i+1)
            !do j=HpointerB(i),HpointerE(i)-1
            beff(:,i,k)=0.0_dblprec
            do j=HpointerB(i)+1,HpointerE(i)
               !print '(3i4,9f12.4)',i,j,Hcolumns(j)+1,Hblock(:,:,j)
               do m=1,3
                  do n=1,3
                    beff(m,i,k)=beff(m,i,k)+Hblock(n,m,j)*emomM(n,Hcolumns(j)+1,k)
                  end do
               end do
            end do
         end do
         !$omp end parallel do
         !do i=1,Natom
         !   print '(i4,4x,3f10.4,4x,3f10.4)', i, emomM(:,i,k), beff(:,i,k)
         !end do
#endif
      end do
      !
   end subroutine effective_field_SparseBlock

   subroutine effective_field_SparseScalar(Natom,Mensemble,emomM,beff)
      !
      implicit none
      !
      !include 'mkl.fi'
      !
      integer, intent(in) :: Natom  !< Number of atoms in system
      integer, intent(in) :: Mensemble !< Number of ensembles 
      real(dblprec), dimension(3,Natom,Mensemble),intent(in) :: emomM  !< Current magnetic moment vector
      real(dblprec), dimension(3,Natom,Mensemble), intent(out) :: beff !< Total effective field from application of Hamiltonian
      !
      !
      integer :: i, j, k, idx
      character, dimension(6) :: matdescra

      real(dblprec) :: mone = 1.0_dblprec
      real(dblprec) :: zero = 0.0_dblprec

      !
      ! zero based
      matdescra=(/'S','L','N','C','0','0'/)
      idx=0
      beff=0.0_dblprec
      do k=1,Mensemble
#ifdef __INTEL_MKL__
         call mkl_dcsrmm('N',Natom,3,Natom, mone, matdescra, Hscalar, Hcolumns, HpointerB, HpointerE, emomM(:,:,k), 3, zero, beff(:,:,k),3)
         !i_stat=mkl_sparse_d_mm(SPARSE_OPERATION_TRANSPOSE,mone,H_sparse,H_desc,SPARSE_LAYOUT_COLUMN_MAJOR,emomM(:,:,k),3,Natom,zero,beff(:,:,k),Natom)
!        MKL_SPARSE_D_MM(operation,alpha,A,descr,layout,x,columns,ldx,beta,y,ldy)
#else
         !call ErrorHandling_ERROR('Sparse operations only implemented for mkl')
         !$omp parallel do default(shared) private(i,j)
         do i=1,Natom
            beff(:,i,k)=0.0_dblprec
            do j=HpointerB(i)+1,HpointerE(i)
               !print *,i,j,Hcolumns(j)+1,Hscalar(j)
               beff(1:3,i,k)=beff(1:3,i,k)+Hscalar(j)*emomM(1:3,Hcolumns(j)+1,k)
            end do
         end do
         !$omp end parallel do
         !do i=1,Natom
         !print '(i4,4x,3f10.4,4x,3f10.4)', i, emomM(:,i,k), beff(:,i,k)
         !end do
#endif
      end do
      !
   end subroutine effective_field_SparseScalar


   subroutine allocateSparse(Natom,nelem,flag)
      !
      use Profiling
      !
      implicit none
      !
      integer, intent(in) :: Natom
      integer, intent(in) :: nelem
      integer, intent(in) :: flag
      !
      integer :: i_stat,i_all
      !
      if(flag>0) then
         allocate(Hblock(Hblocksize,Hblocksize,nelem),stat=i_stat)
         !allocate(Hblock(nelem,Hblocksize,Hblocksize),stat=i_stat)
         call memocc(i_stat,product(shape(Hblock))*kind(Hblock),'Hblock','allocateSparse')
         allocate(Hscalar(nelem),stat=i_stat)
         call memocc(i_stat,product(shape(Hscalar))*kind(Hscalar),'Hscalar','allocateSparse')
         allocate(Hcolumns(nelem),stat=i_stat)
         call memocc(i_stat,product(shape(Hcolumns))*kind(Hcolumns),'Hcolumns','allocSparse')
         allocate(HpointerB(Natom+1),stat=i_stat)
         call memocc(i_stat,product(shape(HpointerB))*kind(HpointerB),'HpointerB','allocSparse')
         allocate(HpointerE(Natom+1),stat=i_stat)
         call memocc(i_stat,product(shape(HpointerE))*kind(HpointerE),'HpointerE','allocSparse')
      else if(flag<0) then
         i_all=-product(shape(Hblock))*kind(Hblock)
         deallocate(Hblock,stat=i_stat)
         call memocc(i_stat,i_all,'Hblock','allocateSparse')
         i_all=-product(shape(Hscalar))*kind(Hscalar)
         deallocate(Hscalar,stat=i_stat)
         call memocc(i_stat,i_all,'Hscalar','allocateSparse')
         deallocate(Hcolumns,stat=i_stat)
         call memocc(i_stat,i_all,'Hcolumns','allocSparse')
         i_all=-product(shape(HpointerB))*kind(HpointerB)
         deallocate(HpointerB,stat=i_stat)
         call memocc(i_stat,i_all,'HpointerB','allocSparse')
         i_all=-product(shape(HpointerE))*kind(HpointerE)
         deallocate(HpointerE,stat=i_stat)
         call memocc(i_stat,i_all,'HpointerE','allocSparse')
         !i_all=-product(shape())*kind()
         !deallocate(,stat=i_stat)
         !call memocc(i_stat,i_all,'','allocateSparse')
      end if
      !
   end subroutine allocateSparse

end module Sparse
