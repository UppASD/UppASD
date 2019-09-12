!-------------------------------------------------------------------------------
! MODULE: NeighbourMap
!> @brief Routines for creating neighbour maps for non-local interactions
!> @author
!> Anders Bergman
!> @copyright
!> GNU Public License.
!-------------------------------------------------------------------------------
module NeighbourMap
   use Parameters
   use Profiling
   ! temporary
   use inputdata, only : do_hoc_debug
   use Math_functions
   
   implicit none

   integer, dimension(:,:), allocatable :: nmdimt !< Temporary storage of dimensions for neighbour map
   real(dblprec), dimension(:,:,:), allocatable :: sym_mats !< Symmetry operation matrices
   integer :: nsym !< Number of symmetry operations

   private
   public :: setup_nm, setup_nm_nelem

contains

   !----------------------------------------------------------------------------
   ! SUBROUTINE: setup_nm
   !> Set up the neighbour maps for pair wise interactions taking into consideration
   !> boundary conditions
   !----------------------------------------------------------------------------
   subroutine setup_nm(Natom,NT,NA,N1,N2,N3,C1,C2,C3,BC1,BC2,BC3,block_size,atype,  &
      Bas,max_no_neigh,max_no_shells,max_no_equiv,sym,nn,redcoord,nm,nmdim,         &
      do_ralloy,Natom_full,acellnumb,atype_ch)
      !
      implicit none
      !
      integer, intent(in) :: NT              !< Number of types of atoms
      integer, intent(in) :: NA              !< Number of atoms in one cell
      integer, intent(in) :: N1              !< Number of cell repetitions in x direction
      integer, intent(in) :: N2              !< Number of cell repetitions in y direction
      integer, intent(in) :: N3              !< Number of cell repetitions in z direction
      integer, intent(in) :: sym             !< Symmetry of system (0-3)
      integer, intent(in) :: Natom           !< Number of atoms in system
      integer, intent(in) :: Natom_full      !< Number of atoms for full system (=Natom if not dilute)
      integer, intent(in) :: block_size      !< Size of the blocking parameter for the macro cell creation
      integer, intent(in) :: max_no_shells   !< Calculated maximum of shells for exchange
      integer, dimension(NT), intent(in)     :: nn    !< Number of neighbour shells
      integer, dimension(Natom), intent(in)  :: atype !< Type of atom
      character(len=1), intent(in) :: BC1    !< Boundary conditions in x-direction
      character(len=1), intent(in) :: BC2    !< Boundary conditions in y-direction
      character(len=1), intent(in) :: BC3    !< Boundary conditions in z-direction
      real(dblprec), dimension(3), intent(in) :: C1 !< First lattice vector
      real(dblprec), dimension(3), intent(in) :: C2 !< Second lattice vector
      real(dblprec), dimension(3), intent(in) :: C3 !< Third lattice vector
      real(dblprec), dimension(3,NA), intent(in) :: Bas !< Coordinates for basis atoms
      real(dblprec), dimension(NT,max_no_shells,3), intent(in) :: redcoord !< Coordinates for Heisenberg exchange couplings
      ! .. Output variables
      integer, intent(out) :: max_no_neigh !< Calculated maximum of neighbours for exchange
      integer, intent(out) :: max_no_equiv !< Calculated maximum of neighbours in one shell for exchange
      integer, dimension(:,:,:), allocatable, intent(out) :: nm !< Neighbour map
      integer, dimension(:,:), allocatable, intent(out) :: nmdim !< Dimension of neighbour map
      ! .. Optional variables
      integer, intent(in), optional :: do_ralloy  !< Random alloy simulation (0/1)
      integer, dimension(Natom_full), optional ,intent(in) :: acellnumb !< List for translating atom no. in full cell to actual cell
      integer, dimension(Natom_full), optional ,intent(in) :: atype_ch !< Actual type of atom for dilute system
      !
      integer :: i0
      integer :: i, nelem=1
      integer :: i_stat,i_all
      integer :: nndim
      real(dblprec) :: tol
      real(dblprec) :: detmatrix
      real(dblprec), dimension(3,3) :: invmatrix
      integer :: ix,iy,iz,xc_hop,yc_hop,zc_hop,ia,counter,ishell,itype,inei
      integer :: iat,jat,j0,jx,jy,jz, iix, iiy, iiz
      real(dblprec), dimension(3) :: bsf,cvec,icvec,rvec
      real(dblprec), dimension(:,:,:,:), allocatable:: nncoord !< Full list of neighbours for each type
      integer, dimension(:,:,:), allocatable:: nm_cell
      integer, dimension(:,:,:,:), allocatable :: nm_trunk
      integer, dimension(:,:), allocatable :: nnm_cell
      logical :: is_periodic,is_dilute

      ! Set tolerance
      tol=0.00050_dblprec**2

      ! calculate max.no. of shells of neighbours
      nndim=0
      do i=1,NA
         nndim=max(nndim,NN(atype(i)))
      enddo

      ! Calculate inverse of basis matrix
      detmatrix=  C1(1)*C2(2)*C3(3)-C1(1)*C2(3)*C3(2)+&
                  C1(2)*C2(3)*C3(1)-C1(2)*C2(1)*C3(3)+&
                  C1(3)*C2(1)*C3(2)-C1(3)*C2(2)*C3(1)
      invmatrix=0.0_dblprec

      if(abs(detmatrix)>dbl_tolerance) then
         invmatrix(1,1)=(C2(2)*C3(3)-C3(2)*C2(3))/detmatrix
         invmatrix(1,2)=(C1(3)*C3(2)-C3(3)*C1(2))/detmatrix
         invmatrix(1,3)=(C1(2)*C2(3)-C2(2)*C1(3))/detmatrix
         invmatrix(2,1)=(C2(3)*C3(1)-C3(3)*C2(1))/detmatrix
         invmatrix(2,2)=(C1(1)*C3(3)-C3(1)*C1(3))/detmatrix
         invmatrix(2,3)=(C1(3)*C2(1)-C2(3)*C1(1))/detmatrix
         invmatrix(3,1)=(C2(1)*C3(2)-C3(1)*C2(2))/detmatrix
         invmatrix(3,2)=(C1(2)*C3(1)-C3(2)*C1(1))/detmatrix
         invmatrix(3,3)=(C1(1)*C2(2)-C2(1)*C1(2))/detmatrix
      end if
      max_no_neigh=0

      ! create all symmetry matrices wrt symmetry type
      call get_symops(sym)

      ! Allocate arrays
      ! Maximum no. of neighbours in each shell is determined by the symmetry (max=48)
      max_no_equiv=nsym
      ! neighbour coordinates
      allocate(nncoord(3,max_no_equiv,max_no_shells,nt),stat=i_stat)
      call memocc(i_stat,product(shape(nncoord))*kind(nncoord),'nncoord','setup_nm')
      nncoord=0.0_dblprec
      ! neighbour map information 1: atom no in first supercell
      allocate(nm_cell(max_no_equiv,max_no_shells,na),stat=i_stat)
      call memocc(i_stat,product(shape(nm_cell))*kind(nm_cell),'nm_cell','setup_nm')
      nm_cell=0
      ! neighbour map information 2: translation in basis vectors
      ! so that full neighbour vector is nm_cell+nm_trunk
      is_periodic=(BC1=='P'.or.BC2=='P'.or.BC3=='P')
      is_dilute=(N1*N2*N3>1.and.is_periodic)
      if(N1*N2*N3>1) then
         allocate(nm_trunk(3,max_no_equiv,max_no_shells,na),stat=i_stat)
         call memocc(i_stat,product(shape(nm_trunk))*kind(nm_trunk),'nm_trunk','setup_nm')
         nm_trunk=0
      end if
      ! max no. neighbours in each shell
      allocate(nnm_cell(max_no_shells,na),stat=i_stat)
      call memocc(i_stat,product(shape(nnm_cell))*kind(nnm_cell),'nnm_cell','setup_nm')
      nm_cell=0
      ! Allocate arrays for neighbour map dimensions
      allocate(nmdim(maxval(NN),Natom),stat=i_stat)
      call memocc(i_stat,product(shape(nmdim))*kind(nmdim),'nmdim','setup_nm')
      allocate(nmdimt(maxval(NN),na),stat=i_stat)
      call memocc(i_stat,product(shape(nmdimt))*kind(nmdimt),'nmdimt','setup_nm')
      nmdim=0
      nmdimt=0
      ! Create full neighbour list according to symmetry
      call get_fullnnlist(NT,NN,nelem,max_no_shells,redcoord,max_no_equiv,nncoord)
      ! Start looking in "first cell"
      if(do_hoc_debug==1) then
         write(*,*) 'Constructs nm_cell and nm_trunk'
      end if
      !$omp parallel do default(shared) private(i0,itype,ishell,counter,inei,cvec,icvec,bsf,rvec,ia)
      do I0=1,NA
         itype=atype(I0)
         if(do_hoc_debug==1) then
            write(*,'(a,i4,a)') '------- i0 ', i0, ' -------'
         end if

         ! Shell
         do ishell=1,NN(itype)
            counter=0

            ! Symmetry equivalent sites in shell
            do inei=1,nmdimt(ishell,itype)

               if(do_hoc_debug==1) then
                  write(*,'(a,i4)') 'counter ', counter
                  write(*,'(a,i4,a,3f10.6)') ' inei ', inei, '   nncoord ', nncoord(1:3,inei,ishell,itype)
                  !write(*,'(a,i4,a,i4,a,3f10.6)') ' ielem ', ielem, ' inei ', inei, '   nncoord ', nncoord(1:3,ielem,inei,ishell,itype)
               end if

               ! Coordinate vector in cartesian coordinates
               cvec(1)=nncoord(1,inei,ishell,itype)+Bas(1,i0)
               cvec(2)=nncoord(2,inei,ishell,itype)+Bas(2,i0)
               cvec(3)=nncoord(3,inei,ishell,itype)+Bas(3,i0)
               if(do_hoc_debug==1) write(*,'(a,3f10.6)') 'cvec  ',cvec(:)

               ! Find coordinate vector in basis coordinates
               icvec(1)=cvec(1)*invmatrix(1,1)+cvec(2)*invmatrix(2,1)+cvec(3)*invmatrix(3,1)
               icvec(2)=cvec(1)*invmatrix(1,2)+cvec(2)*invmatrix(2,2)+cvec(3)*invmatrix(3,2)
               icvec(3)=cvec(1)*invmatrix(1,3)+cvec(2)*invmatrix(2,3)+cvec(3)*invmatrix(3,3)
               if(do_hoc_debug==1) write(*,'(a,3f10.6)')'icvec ',icvec(:)

               ! Fold back to original cell
               bsf(1)=floor(icvec(1)+1d-6)
               bsf(2)=floor(icvec(2)+1d-6)
               bsf(3)=floor(icvec(3)+1d-6)
               if(do_hoc_debug==1) write(*,'(a,3f10.6)') 'bsf   ',bsf(:)

               ! Corresponding position of atom in cell
               rvec(1)=cvec(1)-bsf(1)*C1(1)-bsf(2)*C2(1)-bsf(3)*C3(1)
               rvec(2)=cvec(2)-bsf(1)*C1(2)-bsf(2)*C2(2)-bsf(3)*C3(2)
               rvec(3)=cvec(3)-bsf(1)*C1(3)-bsf(2)*C2(3)-bsf(3)*C3(3)
               if(do_hoc_debug==1) write(*,'(a,3f10.6)')  'rvec  ',rvec(:)

               ! loop through atoms in cell to find match
               do ia=1,NA
                  if((rvec(1)-Bas(1,ia))**2+ &
                     (rvec(2)-Bas(2,ia))**2+ &
                     (rvec(3)-Bas(3,ia))**2<tol) then

                     if(do_hoc_debug==1) write(*,'(a,i4,a,i4)') 'hit!  i0 ', i0, ' ia ', ia

                     counter=counter+1
                     if(max_no_equiv>=counter) then
                        nm_cell(counter,ishell,i0)=ia
                        if(N1*N2*N3>1) then
                           nm_trunk(1,counter,ishell,i0)=nint(bsf(1))
                           nm_trunk(2,counter,ishell,i0)=nint(bsf(2))
                           nm_trunk(3,counter,ishell,i0)=nint(bsf(3))
                        end if
                     end if
                  end if
               end do
            end do
            nnm_cell(ishell,i0)=counter
         end do
      end do
      !$omp end parallel do
      i_all=-product(shape(nncoord))*kind(nncoord)
      deallocate(nncoord,stat=i_stat)
      call memocc(i_stat,i_all,'nncoord','setup_nm')
      ! Allocate nm : neighbour map
      allocate(nm(nAtom,maxval(Nn),max_no_equiv),stat=i_stat)
      call memocc(i_stat,product(shape(nm))*kind(nm),'nm','setup_nm')
      nm=0
      ! With info in nm_cell and nm_trunk for first NA atoms, make full nm list
      do iiz=0, N3-1, block_size
         do iiy=0, N2-1, block_size
            do iix=0, N1-1, block_size
               do iz=iiz,min(iiz+block_size-1,N3-1)
                  do iy=iiy,min(iiy+block_size-1,N2-1)
                     do ix=iix,min(iix+block_size-1,N1-1)
                        do i0=1, NA
                           itype=atype(i0)
                           iat=i0+ix*NA+iy*N1*NA+iz*N2*N1*NA
                           if (do_ralloy==1) then
                              iat = acellnumb(iat)
                              if (iat /= 0) then
                                 itype = atype_ch(iat)
                              end if
                           end if
                           if (iat/=0) then
                              ! Shell
                              do ishell=1,NN(itype)
                                 if (do_ralloy==1) then
                                    nmdim(ishell,iat)=nmdimt(ishell,atype_ch(iat))
                                 else
                                    nmdim(ishell,iat)=nmdimt(ishell,atype(iat))
                                 end if
                                 ! Site in shell
                                 do inei=1,nnm_cell(ishell,i0)
                                    ! Designation of cell
                                    if(N1*N2*N3>1) then
                                       xc_hop=nm_trunk(1,inei,ishell,i0)
                                       yc_hop=nm_trunk(2,inei,ishell,i0)
                                       zc_hop=nm_trunk(3,inei,ishell,i0)
                                    else
                                       xc_hop=0
                                       yc_hop=0
                                       zc_hop=0
                                    end if
                                    ! Position in cell
                                    j0=nm_cell(inei,ishell,i0)
                                    ! Wrap around if periodic boundaries
                                    jx=xc_hop+ix
                                    if(BC1=='P') then
                                       jx=mod(jx+1000*N1,N1)
                                    else if (N1*N2*N3<=1) then
                                       jx=0
                                    end if
                                    jy=yc_hop+iy
                                    if(BC2=='P') then
                                       jy=mod(jy+1000*N2,N2)
                                    else if (N1*N2*N3<=1) then
                                       jy=0
                                    end if
                                    jz=zc_hop+iz
                                    if(BC3=='P') then
                                       jz=mod(jz+1000*N3,N3)
                                    else if (N1*N2*N3<=1) then
                                       jz=0
                                    end if
                                    ! See if atom exists, then add entry to neighbour map
                                    if(jx>=0.and.jx<N1.and.jy>=0.and.jy<N2.and.jz>=0.and.jz<N3) then
                                       jat=j0+jx*NA+jy*N1*NA+jz*N2*N1*NA
                                       if (do_ralloy==1) jat = acellnumb(jat)
                                       nm(iat,ishell,inei)=jat
                                    end if
                                    !
                                 end do
                              end do
                           end if
                        end do
                     end do
                  end do
               end do
            end do
         end do
      end do

      ! Deallocate
      if(N1*N2*N3>1) then
         i_all=-product(shape(nm_trunk))*kind(nm_trunk)
         deallocate(nm_trunk,stat=i_stat)
         call memocc(i_stat,i_all,'nm_trunk','setup_nm')
      end if

      i_all=-product(shape(nm_cell))*kind(nm_cell)
      deallocate(nm_cell,stat=i_stat)
      call memocc(i_stat,i_all,'nm_cell','setup_nm')
      !
      i_all=-product(shape(sym_mats))*kind(sym_mats)
      deallocate(sym_mats,stat=i_stat)
      call memocc(i_stat,i_all,'sym_mats','setup_nm')
      !
      ! Calculate maximum number of neighbours
      max_no_neigh=1
      do i0=1, NA
         itype=atype(i0)
         counter=0
         do ishell=1,NN(itype)
            counter = counter +nnm_cell(ishell,i0)
         end do
         max_no_neigh=max(max_no_neigh,counter)
      end do

      i_all=-product(shape(nnm_cell))*kind(nnm_cell)
      deallocate(nnm_cell,stat=i_stat)
      call memocc(i_stat,i_all,'nnm_cell','setup_nm')
      i_all=-product(shape(nmdimt))*kind(nmdimt)
      deallocate(nmdimt,stat=i_stat)
      call memocc(i_stat,i_all,'nmdimt','setup_nm')
      !
   end subroutine setup_nm

   !> Find possible symmetry operations depending on assumed symmetry
   subroutine get_symops(isym)
      !
      implicit none
      !
      integer, intent(in) :: isym !< Type of assumed symmetry (0-3)
      !
      integer :: i,j,x,y,z,j_s,x1,x2,y1,y2
      integer :: i_stat
      integer :: sym_count
      real(dblprec) :: half,roothalf
      !

      sym_count=0

      ! No symmetry
      if (isym==0) then

         sym_count=1
         allocate(sym_mats(3,3,1),stat=i_stat)
         call memocc(i_stat,product(shape(sym_mats))*kind(sym_mats),'sym_mats','get_symops')
         sym_mats=0.0_dblprec
         do i=1,3
            sym_mats(i,i,1)=1.0_dblprec
         end do

         ! Cubic symmetry
      else if(isym==1) then

         allocate(sym_mats(3,3,48),stat=i_stat)
         call memocc(i_stat,product(shape(sym_mats))*kind(sym_mats),'sym_mats','get_symops')
         sym_mats=0.0_dblprec
         sym_count=0
         do i=1,3
            do j=0,1
               j_s=(-1)**j
               do x=0,1
                  do y=0,1
                     do z=0,1
                        sym_count=sym_count+1
                        sym_mats(1,mod(i-j_s,3)+1,sym_count)=(-1.0_dblprec)**x
                        sym_mats(2,mod(i,3)+1,sym_count)=(-1.0_dblprec)**y
                        sym_mats(3,mod(i+j_s,3)+1,sym_count)=(-1.0_dblprec)**z
                     end do
                  end do
               end do
            end do
         end do

         ! Cubic symmetry in xy-plane
      else if(isym==2) then
         allocate(sym_mats(3,3,12),stat=i_stat)
         call memocc(i_stat,product(shape(sym_mats))*kind(sym_mats),'sym_mats','get_symops')
         sym_mats=0.0_dblprec
         sym_count=0
         do j=0,1
            do x=0,1
               do y=0,1
                  sym_count=sym_count+1
                  sym_mats(1,mod(j,2)+1,sym_count)=(-1.0_dblprec)**x
                  sym_mats(2,mod(j+1,2)+1,sym_count)=(-1.0_dblprec)**y
                  sym_mats(3,3,sym_count)=1.0_dblprec
               end do
            end do
         end do

         ! Hexagonal symmetry
      else if(isym==3) then

         allocate(sym_mats(3,3,24),stat=i_stat)
         call memocc(i_stat,product(shape(sym_mats))*kind(sym_mats),'sym_mats','get_symops')
         sym_mats=0.0_dblprec
         sym_count=0
         half=0.50_dblprec
         roothalf=sqrt(3.0_dblprec)*0.50_dblprec
         ! 8 ops due to 'cartesian' inversion
         do x=0,1
            do y=0,1
               do z=0,1
                  sym_count=sym_count+1
                  sym_mats(1,1,sym_count)=(-1.0_dblprec)**x
                  sym_mats(2,2,sym_count)=(-1.0_dblprec)**y
                  sym_mats(3,3,sym_count)=(-1.0_dblprec)**z
               end do
            end do
         end do
         ! 16 ops due to 'cartesian' inversion
         do x1=0,1
            do x2=0,1
               do y1=0,1
                  do y2=0,1
                     if((-1.0_dblprec)**x1*(-1.0_dblprec)**x2*(-1.0_dblprec)**y1*(-1.0_dblprec)**y2<0.0_dblprec) then
                        do z=0,1
                           sym_count=sym_count+1
                           sym_mats(1,1,sym_count)=(-1.0_dblprec)**x1*half
                           sym_mats(2,1,sym_count)=(-1.0_dblprec)**x2*roothalf
                           sym_mats(1,2,sym_count)=(-1.0_dblprec)**y1*roothalf
                           sym_mats(2,2,sym_count)=(-1.0_dblprec)**y2*half
                           sym_mats(3,3,sym_count)=(-1.0_dblprec)**z
                        end do
                     end if
                  end do
               end do
            end do
         end do

         ! Reads symmetry operations from file
      else if(isym==4) then
         open(ifileno,file='sym.mat')
         read(ifileno,*) sym_count
         allocate(sym_mats(3,3,sym_count),stat=i_stat)
         call memocc(i_stat,product(shape(sym_mats))*kind(sym_mats),'sym_mats','get_symops')
         do j=1,sym_count
            do x=1,3
               read(ifileno,*) (sym_mats(y,x,j),y=1,3)
            end do
         end do
         close(ifileno)
      end if
      nsym=sym_count
      !
   end subroutine get_symops

   !----------------------------------------------------------------------------
   ! SUBROUTINE: get_fullnnlist
   !> @todo Consider the full space group of the structure. Make use of space group libraries. The ELK space group modules?
   !> @todo The if statements to see if a coupling ij is present would then be unnecessary
   !> @todo A possible intermediate construction is to enable use of a distinct point group symmetry for each coupling
   !> @todo This intermediate construction is halfway complete but for now inactivated.
   !> Create full neighbour list according to symmetry
   !----------------------------------------------------------------------------
   subroutine get_fullnnlist(NT,NN,Nelem,max_no_shells,redcoord,max_no_equiv,nncoord)
      !
      !
      implicit none
      !
      integer, intent(in) :: nt !< Number of types of atoms
      integer, intent(in) :: nelem  !< Number of elements in each coupling (=1 for Heisenberg or other two-site couplings,=3 for 4-spin ring exchange)
      integer, intent(in) :: max_no_equiv !< Calculated maximum of neighbours in one shell for exchange
      integer, intent(in) :: max_no_shells !< Calculated maximum of shells for exchange
      integer, dimension(NT), intent(in) :: nn !< Number of neighbour shells
      real(dblprec), dimension(NT,max_no_shells,3,nelem), intent(in) :: redcoord !< Coordinates for Heisenberg exchange couplings
      real(dblprec), dimension(3,nelem,max_no_equiv,max_no_shells,nt), intent(out) :: nncoord !< Full list of neighbours for each type

      real(dblprec) :: tol
      real(dblprec), dimension(3) :: tvect
      real(dblprec), dimension(3,3) :: tvectelem
      integer :: counter
      logical :: unique
      integer :: i, j, k ,itype, ishell, isym
      integer :: ielem

      !0:th order. B = A (transformation of a scalar)
      !1:st order. Bi = Rip Ap (transformation of a vector). Use active rotation instead? Passive and active rotations the same for vectors apart from an overall sign?
      !2:nd order. Bij = Rpi Rqj Apq (transformation of a rank-2 tensor). Use active rotation instead?
      !3:th order. Bijk = Rpi Rqj Rrk Apqr (transformation of a rank-3 tensor). Use active rotation instead?
      !4:th order. Bijkl = Rpi Rqj Rrk Rsl Apqrs (transformation of a rank-4 tensor). Use active rotation instead?

      ! Tolerance
      ! To be tuned!
      tol=0.00050_dblprec

      nncoord = 0.0_dblprec
      if(do_hoc_debug==1) then
         write(*,*) 'Symmetry equivalent bond vectors'
      end if
      do itype=1,nt
         do ishell=1,NN(itype)
            if (nsym==1) then
               do k=1,3
                  nncoord(k,1:nelem,1,ishell,itype)=redcoord(itype,ishell,k,1:nelem)
               end do
               if(do_hoc_debug==1) then
                  write(*,*) '----'
                  write(*,'(a,i4,a,i4)') 'itype ', itype, ' ishell ', ishell
                  do ielem=1,nelem
                     write(*,'(a,3f10.6)') 'redcoord ', redcoord(itype,ishell,1:3,ielem)
                  end do
                  do ielem=1,nelem
                     write(*,'(a,3f10.6)') 'nncoord  ', nncoord(1:3,ielem,1,ishell,itype)
                  end do
               end if
               nmdimt(ishell,itype)=1
            else
               counter=0
               ! Loop over symmetries
               do isym=1,nsym
                  tvect=0.0_dblprec
                  tvectelem=0.0_dblprec
                  unique=.true.
                  do i=1,3
                     do j=1,3
                        tvectelem(i,1:nelem)=tvectelem(i,1:nelem) +&
                           redcoord(itype,ishell,j,1:nelem)*sym_mats(i,j,isym)
                     end do
                  end do
                  do k=1,counter
                     if(nelem == 1) then
                        if( (tvectelem(1,1)-nncoord(1,1,k,ishell,itype))**2 + &
                           (tvectelem(2,1)-nncoord(2,1,k,ishell,itype))**2 + &
                           (tvectelem(3,1)-nncoord(3,1,k,ishell,itype))**2 < tol) unique = .false.
                     end if
                     if(nelem == 2) then
                        if( (tvectelem(1,1)-nncoord(1,1,k,ishell,itype))**2 + &
                           (tvectelem(2,1)-nncoord(2,1,k,ishell,itype))**2 + &
                           (tvectelem(3,1)-nncoord(3,1,k,ishell,itype))**2 + &
                           (tvectelem(1,2)-nncoord(1,2,k,ishell,itype))**2 + &
                           (tvectelem(2,2)-nncoord(2,2,k,ishell,itype))**2 + &
                           (tvectelem(3,2)-nncoord(3,2,k,ishell,itype))**2 < tol) unique = .false.
                     end if
                     if(nelem == 3) then
                        if( (tvectelem(1,1)-nncoord(1,1,k,ishell,itype))**2 + &
                           (tvectelem(2,1)-nncoord(2,1,k,ishell,itype))**2 + &
                           (tvectelem(3,1)-nncoord(3,1,k,ishell,itype))**2 + &
                           (tvectelem(1,2)-nncoord(1,2,k,ishell,itype))**2 + &
                           (tvectelem(2,2)-nncoord(2,2,k,ishell,itype))**2 + &
                           (tvectelem(3,2)-nncoord(3,2,k,ishell,itype))**2 + &
                           (tvectelem(1,3)-nncoord(1,3,k,ishell,itype))**2 + &
                           (tvectelem(2,3)-nncoord(2,3,k,ishell,itype))**2 + &
                           (tvectelem(3,3)-nncoord(3,3,k,ishell,itype))**2 < tol) unique = .false.
                     end if
                  end do
                  if (unique) then
                     counter=counter+1
                     do i=1,3
                        nncoord(i,1:nelem,counter,ishell,itype)=tvectelem(i,1:nelem)
                     end do
                     if(do_hoc_debug==1) then
                        write(*,*) '----'
                        write(*,'(a,i4,a,i4,a,i4)') 'itype ', itype, ' ishell ', ishell, ' counter ', counter
                        do ielem=1,nelem
                           write(*,'(a,3f10.6)') 'redcoord ', redcoord(itype,ishell,1:3,ielem)
                        end do
                        do ielem=1,nelem
                           write(*,'(a,3f10.6)') 'nncoord  ', nncoord(1:3,ielem,counter,ishell,itype)
                        end do
                     end if
                  end if
               end do
               nmdimt(ishell,itype)=counter
            end if
         end do
      end do
      !
   end subroutine get_fullnnlist

   !----------------------------------------------------------------------------
   ! SUBROUTINE: setup_nm_nelem
   !> @brief Set up neighbour maps
   !----------------------------------------------------------------------------
   subroutine setup_nm_nelem(Natom, NT, NA, N1, N2, N3, C1, C2, C3, BC1, BC2, BC3, atype, Bas, &
         max_no_neigh, max_no_shells, max_no_equiv, sym, &
         nn, redcoord, nm, nmdim, nelem, &
         do_ralloy, Natom_full, acellnumb, atype_ch, &
         Nchmax, hdim, do_tens_sym_in, couptensrank, invsym, timesym, couptens, fullcouptens, nm_cell_symind)
         !Nchmax, hdim, do_tens_sym_in, couptensrank, invsym, timesym, couptens, fullcouptens)

      !
      implicit none
      !
      integer, intent(in) :: NT !< Number of types of atoms
      integer, intent(in) :: NA  !< Number of atoms in one cell
      integer, intent(in) :: N1  !< Number of cell repetitions in x direction
      integer, intent(in) :: N2  !< Number of cell repetitions in y direction
      integer, intent(in) :: N3  !< Number of cell repetitions in z direction
      integer, intent(in) :: Natom !< Number of atoms in system
      real(dblprec), dimension(3), intent(in) :: C1 !< First lattice vector
      real(dblprec), dimension(3), intent(in) :: C2 !< Second lattice vector
      real(dblprec), dimension(3), intent(in) :: C3 !< Third lattice vector
      character(len=1), intent(in) :: BC1 !< Boundary conditions in x-direction
      character(len=1), intent(in) :: BC2 !< Boundary conditions in y-direction
      character(len=1), intent(in) :: BC3 !< Boundary conditions in z-direction
      integer, dimension(Natom), intent(in) :: atype !< Type of atom
      real(dblprec), dimension(3,NA), intent(in) :: Bas !< Coordinates for basis atoms
      integer, dimension(NT), intent(in) :: nn !< Number of neighbour shells
      integer, intent(out) :: max_no_neigh !< Calculated maximum of neighbours for exchange
      integer, intent(in) :: max_no_shells !< Calculated maximum of shells for exchange
      integer, intent(in) :: nelem  !< Number of elements in each coupling (=1 for Heisenberg or other two-site couplings,=3 for 4-spin ring exchange)
      real(dblprec), dimension(NT,max_no_shells,3,nelem), intent(in) :: redcoord !< Coordinates for Heisenberg exchange couplings
      integer, intent(out) :: max_no_equiv !< Calculated maximum of neighbours in one shell for exchange
      integer, intent(in) :: sym !< Symmetry of system (0-3)
      integer, dimension(:,:,:,:), allocatable, intent(out) :: nm !< Neighbour map
      integer, dimension(:,:), allocatable, intent(out) :: nmdim !< Dimension of neighbour map
      integer, intent(in), optional :: do_ralloy  !< Random alloy simulation (0/1)
      integer, intent(in) :: Natom_full !< Number of atoms for full system (=Natom if not dilute)
      integer, dimension(Natom_full), optional ,intent(in) :: acellnumb !< List for translating atom no. in full cell to actual cell
      integer, dimension(Natom_full), optional ,intent(in) :: atype_ch !< Actual type of atom for dilute system
      
      integer, intent(in) :: Nchmax !< Number of chemical types of atoms
      integer, intent(in) :: hdim  !< Number of elements in Hamiltonian element (scalar or vector)
      logical, intent(in) :: do_tens_sym_in
      integer, intent(in) :: couptensrank !< Rank of coupling tensor
      integer, intent(in) :: invsym !< Flag for inversion symmetry of coupling tensor
      integer, intent(in) :: timesym !< Flag for time reversal symmetry of coupling tensor
      real(dblprec), dimension(hdim,NT,max_no_shells,Nchmax,Nchmax), intent(in) :: couptens !< Coupling tensor
      real(dblprec), dimension(hdim,NT,max_no_shells,Nchmax,Nchmax,48), intent(out) :: fullcouptens !< Symmetry degenerate coupling tensor
      integer, dimension(48,max_no_shells,na), intent(out) :: nm_cell_symind  !< Indices for elements of the symmetry degenerate coupling tensor
      !
      integer :: i0
      integer :: i, ielem
      integer :: i_stat,i_all
      integer :: nndim
      real(dblprec) :: tol
      real(dblprec) :: detmatrix
      real(dblprec), dimension(3,3) :: invmatrix
      integer :: ix,iy,iz,xc_hop,yc_hop,zc_hop,ia,counter,ishell,itype,inei
      integer :: iat,jat,j0,jx,jy,jz
      real(dblprec), dimension(3) :: bsf,cvec,icvec,rvec
      real(dblprec), dimension(:,:,:,:,:), allocatable:: nncoord !< Full list of neighbours for each type
      integer, dimension(:,:,:,:), allocatable:: nm_cell
      integer, dimension(:,:,:,:,:), allocatable :: nm_trunk
      integer, dimension(:,:), allocatable :: nnm_cell
      logical :: is_periodic,is_dilute
      logical :: hit
      logical, dimension(3) :: elemhit

      ! Set tolerance
      tol=0.0005_dblprec**2

      ! calculate max.no. of shells of neighbours
      nndim=0
      do i=1,NA
         nndim=max(nndim,NN(atype(i)))
      enddo

      ! Calculate inverse of basis matrix
      detmatrix=C1(1)*C2(2)*C3(3)-C1(1)*C2(3)*C3(2)+&
         C1(2)*C2(3)*C3(1)-C1(2)*C2(1)*C3(3)+&
         C1(3)*C2(1)*C3(2)-C1(3)*C2(2)*C3(1)
      invmatrix=0.0_dblprec

      if(abs(detmatrix)>dbl_tolerance) then
         invmatrix(1,1)=(C2(2)*C3(3)-C3(2)*C2(3))/detmatrix
         invmatrix(1,2)=(C1(3)*C3(2)-C3(3)*C1(2))/detmatrix
         invmatrix(1,3)=(C1(2)*C2(3)-C2(2)*C1(3))/detmatrix
         invmatrix(2,1)=(C2(3)*C3(1)-C3(3)*C2(1))/detmatrix
         invmatrix(2,2)=(C1(1)*C3(3)-C3(1)*C1(3))/detmatrix
         invmatrix(2,3)=(C1(3)*C2(1)-C2(3)*C1(1))/detmatrix
         invmatrix(3,1)=(C2(1)*C3(2)-C3(1)*C2(2))/detmatrix
         invmatrix(3,2)=(C1(2)*C3(1)-C3(2)*C1(1))/detmatrix
         invmatrix(3,3)=(C1(1)*C2(2)-C2(1)*C1(2))/detmatrix
      end if
      max_no_neigh=0

      ! create all symmetry matrices wrt symmetry type
      call get_symops(sym)

      ! Allocate arrays
      ! Maximum no. of neighbours in each shell is determined by the symmetry (max=48)
      max_no_equiv=nsym

      ! neighbour coordinates
      allocate(nncoord(3,nelem,max_no_equiv,max_no_shells,nt),stat=i_stat)
      call memocc(i_stat,product(shape(nncoord))*kind(nncoord),'nncoord','setup_nm_nelem')

      ! neighbour map information 1: atom no in first supercell
      allocate(nm_cell(nelem,max_no_equiv,max_no_shells,na),stat=i_stat)
      call memocc(i_stat,product(shape(nm_cell))*kind(nm_cell),'nm_cell','setup_nm_nelem')
      nm_cell=0

      ! Indicies of the symmetry equivalent coupling tensors
      !allocate(nm_cell_symind(max_no_equiv,max_no_shells,na),stat=i_stat)
      !call memocc(i_stat,product(shape(nm_cell_symind))*kind(nm_cell_symind),'nm_cell_symind','setup_nm_nelem')
      ! For case of no symmetry, set all elements to 1.
      if(nsym == 1) then
         nm_cell_symind=1
      else
         nm_cell_symind=0
      end if

      ! neighbour map information 2: translation in basis vectors
      ! so that full neighbour vector is nm_cell+nm_trunk
      is_periodic=(BC1=='P'.or.BC2=='P'.or.BC3=='P')
      is_dilute=(N1*N2*N3>1.and.is_periodic)
      if(N1*N2*N3>1) then
         allocate(nm_trunk(3,nelem,max_no_equiv,max_no_shells,na),stat=i_stat)
         call memocc(i_stat,product(shape(nm_trunk))*kind(nm_trunk),'nm_trunk','setup_nm_nelem')
         nm_trunk=0
      end if

      ! max no. neighbours in each shell
      allocate(nnm_cell(max_no_shells,na),stat=i_stat)
      call memocc(i_stat,product(shape(nnm_cell))*kind(nnm_cell),'nnm_cell','setup_nm_nelem')

      ! Allocate arrays for neighbour map dimensions
      allocate(nmdim(maxval(NN),natom),stat=i_stat)
      call memocc(i_stat,product(shape(nmdim))*kind(nmdim),'nmdim','setup_nm_nelem')
      allocate(nmdimt(maxval(NN),na),stat=i_stat)
      call memocc(i_stat,product(shape(nmdimt))*kind(nmdimt),'nmdimt','setup_nm_nelem')
      nmdim=0
      nmdimt=0
      itype=1

      !!! Create full neighbour list according to symmetry
      !!! call get_fullnnlist(NT, NN, Nelem, max_no_shells, redcoord, max_no_equiv, nncoord)

      ! Create full neighbour list according to symmetry
      call get_fullnnlist_nelem(NT, NN, Nelem, max_no_shells, redcoord, max_no_equiv, nncoord, &
      Nchmax, hdim, do_tens_sym_in, couptensrank, timesym, invsym, couptens, fullcouptens)
      
      ! Start looking in "first cell"
      ! Unclear if thread-safe. Deactive omp for caution. Reactivate after testing
      !!!$omp parallel do default(shared) private(i0,itype,ishell,counter,inei,cvec,icvec,bsf,rvec,ia)
      if(do_hoc_debug==1) then
         write(*,*) 'Constructs nm_cell and nm_trunk'
      end if
      do I0=1,NA
         itype=atype(I0)
         if(do_hoc_debug==1) then
            write(*,'(a,i4,a)') '------- i0 ', i0, ' -------'
         end if

         ! Shell
         do ishell=1,NN(itype)
            if(do_hoc_debug==1) then
               write(*,'(a,i4,a,i4)') 'itype ', itype, ' ishell ', ishell
            end if
            counter=0
            !original set of hit and elemhit

            ! Symmetry equivalent sites in shell
            do inei=1,nmdimt(ishell,itype)
               counter=counter+1
               !original set of hit and elemhit
               hit = .false.
               elemhit = .false.

               ! Check each of the nelem bond vectors separately
               do ielem=1,nelem
                  if(do_hoc_debug==1) then
                     write(*,'(a,i4)') 'counter ', counter
                     write(*,'(a,i4,a,i4,a,3f10.6)') ' ielem ', ielem, ' inei ', inei, '   nncoord ', nncoord(1:3,ielem,inei,ishell,itype)
                  end if
                  ! Coordinate vector in cartesian coordinates
                  cvec(1)=nncoord(1,ielem,inei,ishell,itype)+Bas(1,i0)
                  cvec(2)=nncoord(2,ielem,inei,ishell,itype)+Bas(2,i0)
                  cvec(3)=nncoord(3,ielem,inei,ishell,itype)+Bas(3,i0)
                  if(do_hoc_debug==1) write(*,'(a,3f10.6)') 'cvec  ',cvec(:)

                  ! Find coordinate vector in basis coordinates
                  icvec(1)=cvec(1)*invmatrix(1,1)+cvec(2)*invmatrix(2,1)+cvec(3)*invmatrix(3,1)
                  icvec(2)=cvec(1)*invmatrix(1,2)+cvec(2)*invmatrix(2,2)+cvec(3)*invmatrix(3,2)
                  icvec(3)=cvec(1)*invmatrix(1,3)+cvec(2)*invmatrix(2,3)+cvec(3)*invmatrix(3,3)
                  if(do_hoc_debug==1) write(*,'(a,3f10.6)')'icvec ',icvec(:)

                  ! Fold back to original cell
                  bsf(1)=floor(icvec(1)+1d-6)
                  bsf(2)=floor(icvec(2)+1d-6)
                  bsf(3)=floor(icvec(3)+1d-6)
                  if(do_hoc_debug==1) write(*,'(a,3f10.6)') 'bsf   ',bsf(:)

                  ! Corresponding position of atom in cell
                  rvec(1)=cvec(1)-bsf(1)*C1(1)-bsf(2)*C2(1)-bsf(3)*C3(1)
                  rvec(2)=cvec(2)-bsf(1)*C1(2)-bsf(2)*C2(2)-bsf(3)*C3(2)
                  rvec(3)=cvec(3)-bsf(1)*C1(3)-bsf(2)*C2(3)-bsf(3)*C3(3)
                  if(do_hoc_debug==1) write(*,'(a,3f10.6)')  'rvec  ',rvec(:)

                  ! loop through atoms in cell to find match
                  do ia=1,NA
                     if(  (rvec(1)-Bas(1,ia))**2+ &
                        (rvec(2)-Bas(2,ia))**2+ &
                        (rvec(3)-Bas(3,ia))**2<tol) then
                        elemhit(ielem) = .true.
                        if(do_hoc_debug==1) then
                           write(*,'(a,i4,a,i4)') 'hit!  i0 ', i0, ' ia ', ia
                        end if
                        if(max_no_equiv>=counter) then
                           nm_cell(ielem,counter,ishell,i0)=ia
                           nm_cell_symind(counter,ishell,i0)=inei
                           if(do_hoc_debug==1) then
                              write(*,'(a,i4,a,i4a,i4,a,i4)') 'counter ', counter,  ' ishell ', ishell, ' i0 ',  i0,  ' inei', inei
                           end if
                           if(N1*N2*N3>1) then
                              nm_trunk(1,ielem,counter,ishell,i0)=nint(bsf(1))
                              nm_trunk(2,ielem,counter,ishell,i0)=nint(bsf(2))
                              nm_trunk(3,ielem,counter,ishell,i0)=nint(bsf(3))
                           end if
                        end if
                     end if
                  end do

               end do

               ! Check that each of the 1, 2 or 3 neighbours have been found
               if(nelem == 1) hit = elemhit(1)
               if(nelem == 2) then
                  hit = elemhit(1) .and. elemhit(2)
               end if
               if(nelem == 3) then
                  hit = elemhit(1) .and. elemhit(2) .and. elemhit(3)
               end if
               if(hit .eqv. .false.) then
                  if(N1*N2*N3>1) then
                     nm_trunk(1:3,1:nelem,counter,ishell,i0)=0
                  end if
                  counter = counter -1
               end if
               nnm_cell(ishell,i0)=counter

            end do
         end do
      end do

      if(do_hoc_debug==1) then
         write(*,*) '------ nm_cell --------'
         do i0=1,NA
            do ishell=1, NN(itype)
               do inei=1, nnm_cell(ishell,i0)
                  write(*,'(a,i4,a,i2,a,i2,a,4i4)') 'i0 ', i0, ' ishell ', ishell, ' inei ', inei, &
                     ' : ', i0, nm_cell(1:nelem,inei,ishell,i0)
               end do
            end do
         end do
      end if

      i_all=-product(shape(nncoord))*kind(nncoord)
      deallocate(nncoord,stat=i_stat)
      call memocc(i_stat,i_all,'nncoord','setup_nm_nelem')

      ! Allocate nm : neighbour map
      allocate(nm(Natom,maxval(Nn),max_no_equiv,nelem),stat=i_stat)
      call memocc(i_stat,product(shape(nm))*kind(nm),'nm','setup_nm_nelem')
      nm=0

      if(do_hoc_debug==1) then
         write(*,*) 'Constructs nm'
      end if
      ! With info in nm_cell and nm_trunk for first NA atoms, make full nm list
      do iz=0, N3-1
         do iy=0, N2-1
            do ix=0, N1-1
               do i0=1, NA
                  itype=atype(i0)
                  iat=i0+ix*NA+iy*N1*NA+iz*N2*N1*NA
                  if (do_ralloy==1) then
                     iat = acellnumb(iat)
                     if (iat /= 0) then
                        itype = atype_ch(iat)
                     end if
                  end if

                  if (iat/=0) then

                     ! Shell
                     do ishell=1,NN(itype)
                        if (do_ralloy==1) then
                           nmdim(ishell,iat)=nmdimt(ishell,atype_ch(iat))
                        else
                           nmdim(ishell,iat)=nmdimt(ishell,atype(iat))
                        end if
                        !counter=1

                        ! Site in shell
                        do inei=1,nnm_cell(ishell,i0)

                           hit = .false.
                           elemhit = .false.
                           do ielem=1,nelem
                              !
                              ! Designation of cell
                              if(N1*N2*N3>1) then
                                 xc_hop=nm_trunk(1,ielem,inei,ishell,i0)
                                 yc_hop=nm_trunk(2,ielem,inei,ishell,i0)
                                 zc_hop=nm_trunk(3,ielem,inei,ishell,i0)
                              else
                                 xc_hop=0
                                 yc_hop=0
                                 zc_hop=0
                              end if

                              ! Position in cell
                              j0=nm_cell(ielem,inei,ishell,i0)

                              ! Wrap around if periodic boundaries
                              jx=xc_hop+ix
                              if(BC1=='P') then
                                 jx=mod(jx+1000*N1,N1)
                              else if (N1*N2*N3<=1) then
                                 jx=0
                              end if
                              jy=yc_hop+iy
                              if(BC2=='P') then
                                 jy=mod(jy+1000*N2,N2)
                              else if (N1*N2*N3<=1) then
                                 jy=0
                              end if
                              jz=zc_hop+iz
                              if(BC3=='P') then
                                 jz=mod(jz+1000*N3,N3)
                              else if (N1*N2*N3<=1) then
                                 jz=0
                              end if

                              ! See if atom exists, then add entry to neighbour map
                              if(jx>=0.and.jx<N1.and.jy>=0.and.jy<N2.and.jz>=0.and.jz<N3) then
                                 jat=j0+jx*NA+jy*N1*NA+jz*N2*N1*NA
                                 if (do_ralloy==1) jat = acellnumb(jat)
                                 nm(iat,ishell,inei,ielem)=jat
                                 !nm(iat,ishell,counter,ielem)=jat
                                 elemhit(ielem) = .true.
                              end if
                              !
                           end do

                           ! Check that each of the 1, 2 or 3 neighbours have been found
                           if(nelem == 1) hit = elemhit(1)
                           if(nelem == 2) then
                              hit = elemhit(1) .and. elemhit(2)
                           end if
                           if(nelem == 3) then
                              hit = elemhit(1) .and. elemhit(2) .and. elemhit(3)
                           end if
                           if(hit) then
                              write(*,*) 'All neighbours have been found'
                              !counter = counter + 1
                           else
                              nm(iat,ishell,inei,1:nelem)=0
                           end if

                        end do

                        !counter = counter -1
                        !nmdim(ishell,iat) = counter

                        !if(do_hoc_debug==1) write(*,*) 'nm hit counter ', counter

                        if(do_hoc_debug==1) then
                           write(*,'(a,i4,a,i4,a,i4)') 'iat ', iat, ' ishell ', ishell, ' nmdim ', nmdim(ishell,iat)
                        end if

                     end do

                  end if

               end do
            end do
         end do
      end do

      ! Deallocate
      if(N1*N2*N3>1) then
         i_all=-product(shape(nm_trunk))*kind(nm_trunk)
         deallocate(nm_trunk,stat=i_stat)
         call memocc(i_stat,i_all,'nm_trunk','setup_nm_nelem')
      end if

      i_all=-product(shape(nm_cell))*kind(nm_cell)
      deallocate(nm_cell,stat=i_stat)
      call memocc(i_stat,i_all,'nm_cell','setup_nm_nelem')
      !
      !i_all=-product(shape(nm_cell_symind))*kind(nm_cell_symind)
      !deallocate(nm_cell_symind,stat=i_stat)
      !call memocc(i_stat,i_all,'nm_cell_symind','setup_nm_nelem')
      !
      i_all=-product(shape(sym_mats))*kind(sym_mats)
      deallocate(sym_mats,stat=i_stat)
      call memocc(i_stat,i_all,'sym_mats','setup_nm_nelem')
      !
      ! Calculate maximum number of neighbours
      max_no_neigh=1
      do i0=1, NA
         itype=atype(i0)
         counter=0
         do ishell=1,NN(itype)
            counter = counter +nnm_cell(ishell,i0)
         end do
         max_no_neigh=max(max_no_neigh,counter)
      end do

      i_all=-product(shape(nnm_cell))*kind(nnm_cell)
      deallocate(nnm_cell,stat=i_stat)
      call memocc(i_stat,i_all,'nnm_cell','setup_nm_nelem')
      i_all=-product(shape(nmdimt))*kind(nmdimt)
      deallocate(nmdimt,stat=i_stat)
      call memocc(i_stat,i_all,'nmdimt','setup_nm_nelem')
      !

   end subroutine setup_nm_nelem


      !----------------------------------------------------------------------------
   ! SUBROUTINE: get_fullnnlist
   !> @todo Consider the full space group of the structure. Make use of space group libraries. The ELK space group modules?
   !> @todo The if statements to see if a coupling ij is present would then be unnecessary
   !> @todo A possible intermediate construction is to enable use of a distinct point group symmetry for each coupling
   !> @todo This intermediate construction is halfway complete but for now inactivated.
   !> Create full neighbour list according to symmetry
   !----------------------------------------------------------------------------
   subroutine get_fullnnlist_nelem(NT,NN,Nelem,max_no_shells,redcoord,max_no_equiv,nncoord, &
      Nchmax, hdim, do_tens_sym_in, couptensrank, timesym, invsym, couptens, fullcouptens)
      !subroutine get_fullnnlist(NT, NN, Nelem, max_no_shells, redcoord, max_no_equiv, nncoord)
      !
      !
      implicit none
      !
      integer, intent(in) :: nt !< Number of types of atoms
      integer, intent(in) :: nelem  !< Number of elements in each coupling (=1 for Heisenberg or other two-site couplings,=3 for 4-spin ring exchange)
      integer, intent(in) :: max_no_equiv !< Calculated maximum of neighbours in one shell for exchange
      integer, intent(in) :: max_no_shells !< Calculated maximum of shells for exchange
      integer, dimension(NT), intent(in) :: nn !< Number of neighbour shells
      real(dblprec), dimension(NT,max_no_shells,3,nelem), intent(in) :: redcoord !< Coordinates for Heisenberg exchange couplings
      real(dblprec), dimension(3,nelem,max_no_equiv,max_no_shells,nt), intent(out) :: nncoord !< Full list of neighbours for each type
      integer, intent(in) :: Nchmax !< Number of chemical types of atoms
      integer, intent(in) :: hdim  !< Number of elements in Hamiltonian element (scalar or vector)
      logical, intent(in) :: do_tens_sym_in
      integer, intent(in) :: couptensrank !< Rank of coupling tensor
      integer, intent(in) :: invsym !< Inversion symmetry of coupling tensor
      integer, intent(in) :: timesym !< Time reversal symmetry of coupling tensor
      real(dblprec), dimension(hdim,NT,max_no_shells,Nchmax,Nchmax), intent(in) :: couptens !< Coupling tensor
      real(dblprec), dimension(hdim,NT,max_no_shells,Nchmax,Nchmax,max_no_equiv), intent(out) :: fullcouptens !< Symmetry degenerate coupling tensor

      ! This flag is not needed. Keep for the time being.
      logical :: do_tens_sym
      
      real(dblprec) :: tol
      real(dblprec), dimension(3) :: tvect
      real(dblprec), dimension(3,3) :: tvectelem
      integer :: counter
      logical :: unique
      integer :: i, j, k ,itype, ishell, isym
      integer :: ielem
      integer :: ich, jch
      
      !0:th order. B = A (transformation of a scalar)
      !1:st order. Bi = Rip Ap (transformation of a vector). Use active rotation instead? Passive and active rotations the same for vectors apart from an overall sign?
      !2:nd order. Bij = Rpi Rqj Apq (transformation of a rank-2 tensor). Use active rotation instead?
      !3:th order. Bijk = Rpi Rqj Rrk Apqr (transformation of a rank-3 tensor). Use active rotation instead?
      !4:th order. Bijkl = Rpi Rqj Rrk Rsl Apqrs (transformation of a rank-4 tensor). Use active rotation instead?

      ! Tolerance
      ! To be tuned!
      tol=0.00050_dblprec

      do_tens_sym = .true.

      nncoord = 0.0_dblprec
      fullcouptens = 0.0_dblprec

      if(do_hoc_debug==1) then
         write(*,*) 'Symmetry equivalent bond vectors'
      end if
      do itype=1,nt
         do ishell=1,NN(itype)
            if (nsym==1) then
               do k=1,3
                  nncoord(k,1:nelem,1,ishell,itype)=redcoord(itype,ishell,k,1:nelem)
               end do
               if(do_hoc_debug==1) then
                  write(*,*) '----'
                  write(*,'(a,i4,a,i4)') 'itype ', itype, ' ishell ', ishell
                  do ielem=1,nelem
                     write(*,'(a,3f10.6)') 'redcoord ', redcoord(itype,ishell,1:3,ielem)
                  end do
                  do ielem=1,nelem
                     write(*,'(a,3f10.6)') 'nncoord  ', nncoord(1:3,ielem,1,ishell,itype)
                  end do
               end if
               do counter=1,1
               !do counter=1,max_no_equiv
                  do ich=1,Nchmax
                     do jch=1,Nchmax
                        fullcouptens(1:hdim,itype,ishell,ich,jch,counter) = &
                           couptens(1:hdim,itype,ishell,ich,jch)
                     end do
                  end do
               end do
               nmdimt(ishell,itype)=1
            else
               counter=0
               ! Loop over symmetries
               do isym=1,nsym
                  tvect=0.0_dblprec
                  tvectelem=0.0_dblprec
                  unique=.true.
                  do i=1,3
                     do j=1,3
                        tvectelem(i,1:nelem)=tvectelem(i,1:nelem) +&
                           redcoord(itype,ishell,j,1:nelem)*sym_mats(i,j,isym)
                     end do
                  end do
                  do k=1,counter
                     if(nelem == 1) then
                        if( (tvectelem(1,1)-nncoord(1,1,k,ishell,itype))**2 + &
                           (tvectelem(2,1)-nncoord(2,1,k,ishell,itype))**2 + &
                           (tvectelem(3,1)-nncoord(3,1,k,ishell,itype))**2 < tol) unique = .false.
                     end if
                     if(nelem == 2) then
                        if( (tvectelem(1,1)-nncoord(1,1,k,ishell,itype))**2 + &
                           (tvectelem(2,1)-nncoord(2,1,k,ishell,itype))**2 + &
                           (tvectelem(3,1)-nncoord(3,1,k,ishell,itype))**2 + &
                           (tvectelem(1,2)-nncoord(1,2,k,ishell,itype))**2 + &
                           (tvectelem(2,2)-nncoord(2,2,k,ishell,itype))**2 + &
                           (tvectelem(3,2)-nncoord(3,2,k,ishell,itype))**2 < tol) unique = .false.
                     end if
                     if(nelem == 3) then
                        if( (tvectelem(1,1)-nncoord(1,1,k,ishell,itype))**2 + &
                           (tvectelem(2,1)-nncoord(2,1,k,ishell,itype))**2 + &
                           (tvectelem(3,1)-nncoord(3,1,k,ishell,itype))**2 + &
                           (tvectelem(1,2)-nncoord(1,2,k,ishell,itype))**2 + &
                           (tvectelem(2,2)-nncoord(2,2,k,ishell,itype))**2 + &
                           (tvectelem(3,2)-nncoord(3,2,k,ishell,itype))**2 + &
                           (tvectelem(1,3)-nncoord(1,3,k,ishell,itype))**2 + &
                           (tvectelem(2,3)-nncoord(2,3,k,ishell,itype))**2 + &
                           (tvectelem(3,3)-nncoord(3,3,k,ishell,itype))**2 < tol) unique = .false.
                     end if
                  end do
                  if (unique) then
                     counter=counter+1
                     do i=1,3
                        nncoord(i,1:nelem,counter,ishell,itype)=tvectelem(i,1:nelem)
                     end do
                     if(do_hoc_debug==1) then
                        write(*,*) '----'
                        write(*,'(a,i4,a,i4,a,i4)') 'itype ', itype, ' ishell ', ishell, ' counter ', counter
                        do ielem=1,nelem
                           write(*,'(a,3f10.6)') 'redcoord ', redcoord(itype,ishell,1:3,ielem)
                        end do
                        do ielem=1,nelem
                           write(*,'(a,3f10.6)') 'nncoord  ', nncoord(1:3,ielem,counter,ishell,itype)
                        end do
                     end if
                     if(do_tens_sym) then
                        do ich=1,Nchmax
                           do jch=1,Nchmax
                              write(*,*) 'isym ', isym, 'couptensrank ', couptensrank
                              write(*,'(a)') 'couptens(1:hdim,itype,ishell,ich,jch) '
                              write(*,'(9es14.6)') couptens(1:hdim,itype,ishell,ich,jch)
                              !write(*,'(9es14.6)') 'couptens(1:hdim,itype,ishell,ich,jch) ', couptens(1:hdim,itype,ishell,ich,jch)
                              ! For now the timesym and invsym factors are inactive and set to 1.
                              ! They can be used if permutation operations are added so that e.g. the input of 
                              ! a coupling ijk specifies also the coupling jik. Here two examples for MMU coupling
                              ! (m_i . m_j) u_k --> (m_j . m_i) u_k      invsym=1 , timesym=1
                              ! (m_i . m_j) u_i --> (m_i . m_j) u_j      invsym=-1, timesym=1 
                              ! Additional if-statements
                              if(couptensrank == 0) then
                                 fullcouptens(1:hdim,itype,ishell,ich,jch,counter) = &
                                    couptens(1:hdim,itype,ishell,ich,jch) * timesym * invsym
                              else if(couptensrank == 1) then
                                 fullcouptens(1:hdim,itype,ishell,ich,jch,counter) = &
                                    transt1(couptens(1:hdim,itype,ishell,ich,jch),sym_mats(1:3,1:3,isym)) * timesym * invsym
                              else if (couptensrank == 2) then
                                 fullcouptens(1:hdim,itype,ishell,ich,jch,counter) = &
                                    transt2(couptens(1:hdim,itype,ishell,ich,jch),sym_mats(1:3,1:3,isym)) * timesym * invsym
                              else if(couptensrank == 3) then
                                 fullcouptens(1:hdim,itype,ishell,ich,jch,counter) = &
                                    transt3(couptens(1:hdim,itype,ishell,ich,jch),sym_mats(1:3,1:3,isym)) * timesym * invsym
                              else if(couptensrank == 4) then
                                 fullcouptens(1:hdim,itype,ishell,ich,jch,counter) = &
                                    transt4(couptens(1:hdim,itype,ishell,ich,jch),sym_mats(1:3,1:3,isym)) * timesym * invsym
                              end if
                              write(*,'(a)') 'fullcouptens(1:hdim,itype,ishell,ich,jch,counter) '
                              write(*,'(9es14.6)') fullcouptens(1:hdim,itype,ishell,ich,jch,counter)
                              !write(*,'(9es14.6)') 'fullcouptens(1:hdim,itype,ishell,ich,jch,counter) ', fullcouptens(1:hdim,itype,ishell,ich,jch,counter)
                           end do
                        end do
                     end if
                  end if
               end do
               nmdimt(ishell,itype)=counter
            end if
         end do
      end do
      !
   end subroutine get_fullnnlist_nelem

   
end module NeighbourMap
