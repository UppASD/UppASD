!> Module for findng unique sets of atoms ("metatypes") from coordination numbers
module MetaTypes
   use Parameters
   use Profiling
   !
   implicit none
   !
   integer :: NT_meta
   integer :: NA_meta
   integer, dimension(:), allocatable :: atype_meta
   integer, dimension(:), allocatable :: anumb_meta
   integer, dimension(:), allocatable :: NA_metalist


contains

   subroutine allocate_metatype(Natom,flag)

      implicit none

      integer, intent(in), optional :: Natom !< Number of atoms in system
      integer, intent(in) :: flag  !< Allocate or deallocate (1/-1)

      integer :: i_all, i_stat

      ! If flag > 0 allocate arrays
      if(flag>0) then
         allocate(atype_meta(Natom),stat=i_stat)
         call memocc(i_stat,product(shape(atype_meta))*kind(atype_meta),'atype_meta','allocate_metatype')
      else
         ! Otherwise deallocate arrays
         i_all=-product(shape(atype_meta))*kind(atype_meta)
         deallocate(atype_meta,stat=i_stat)
         call memocc(i_stat,i_all,'atype_meta','allocate_metatype')
      end if

   end subroutine allocate_metatype

   subroutine allocate_metanumb(Natom,flag)

      implicit none

      integer, intent(in), optional :: Natom !< Number of atoms in system
      integer, intent(in) :: flag  !< Allocate or deallocate (1/-1)

      integer :: i_all, i_stat

      ! If flag > 0 allocate arrays
      if(flag>0) then
         allocate(anumb_meta(Natom),stat=i_stat)
         call memocc(i_stat,product(shape(anumb_meta))*kind(anumb_meta),'anumb_meta','allocate_metanumb')
         allocate(NA_metalist(Natom),stat=i_stat)
         call memocc(i_stat,product(shape(NA_metalist))*kind(NA_metalist),'NA_metalist','allocate_metanumb')
      else
         ! Otherwise deallocate arrays
         i_all=-product(shape(anumb_meta))*kind(anumb_meta)
         deallocate(anumb_meta,stat=i_stat)
         call memocc(i_stat,i_all,'anumb_meta','allocate_metanumb')
         i_all=-product(shape(NA_metalist))*kind(NA_metalist)
         deallocate(NA_metalist,stat=i_stat)
         call memocc(i_stat,i_all,'NA_metalist','allocate_metanumb')
      end if

   end subroutine allocate_metanumb

   subroutine find_metatype_coordination(Natom,atype,ham,metatype)
      use HamiltonianDataType

      integer, intent(in) :: Natom    !< Number of atoms
      integer, dimension(Natom) , intent(in) :: atype !< Array of proper atom types
      type(ham_t), intent(in) :: ham  !< Hamiltonian
      integer, intent(in) :: metatype  !< Flag how to decide metatype

      integer :: iatom, nuniq, itype, NT
      integer, dimension(:), allocatable :: nn_max_type
      integer :: i_all, i_stat


      call allocate_metatype(Natom,1)

      NT=maxval(atype)
      !allocate(nn_max_type(nuniq),stat=i_stat)
      allocate(nn_max_type(2*NT),stat=i_stat)
      call memocc(i_stat,product(shape(nn_max_type))*kind(nn_max_type),'nn_max_type','find_metatype_coordination')
      nn_max_type=0

      ! Identify atoms with less then maximum coordination
      if (metatype==1) then
         do itype=1,NT
            do iatom=1,Natom
               if (atype(iatom)==itype) then
                  nn_max_type(itype)=max(nn_max_type(itype),ham%nlistsize(iatom))
               end if
            end do
            do iatom=1,Natom
               if (atype(iatom)==itype) then
                  if (ham%nlistsize(iatom)==nn_max_type(itype)) then
                     atype_meta(iatom)=itype
                  else
                     atype_meta(iatom)=itype+NT
                  end if
               end if
            end do
         end do
         NT_meta=maxval(atype_meta)
      else if (metatype==2) then
         do itype=1,NT
            do iatom=1,Natom
               if (atype(iatom)==itype) then
                  nn_max_type(itype)=max(nn_max_type(itype),ham%nlistsize(iatom))
               end if
            end do
            do iatom=1,Natom
               if (atype(iatom)==itype) then
                  if (ham%nlistsize(iatom)==nn_max_type(itype)) then
                     atype_meta(iatom)=1
                  else
                     atype_meta(iatom)=2
                  end if
               end if
            end do
         end do
         NT_meta=2
      else if (metatype==0) then
         NT_meta=NT
         atype_meta=atype
      end if

      i_all=-product(shape(nn_max_type))*kind(nn_max_type)
      deallocate(nn_max_type,stat=i_stat)
      call memocc(i_stat,i_all,'nn_max_type','find_metatype_coordination')

   end subroutine find_metatype_coordination


   subroutine find_metanumb(Natom,NA,N1,N2,N3,Bas,C1,C2,C3,BC1,BC2,BC3,&
      coord,atype,anumb,do_ralloy,metanumb)

      use HamiltonianDataType

      integer, intent(in) :: Natom    !< Number of atoms
      integer, intent(in) :: NA  !< Number of atoms in one cell
      integer, intent(in) :: N1  !< Number of cell repetitions in x direction
      integer, intent(in) :: N2  !< Number of cell repetitions in y direction
      integer, intent(in) :: N3  !< Number of cell repetitions in z direction
      integer, intent(in) :: do_ralloy    !< Random alloy simulation (0/1)
      real(dblprec), dimension(3), intent(in) :: C1 !< First lattice vector
      real(dblprec), dimension(3), intent(in) :: C2 !< Second lattice vector
      real(dblprec), dimension(3), intent(in) :: C3 !< Third lattice vector
      character(len=1), intent(in) :: BC1    !< Boundary conditions in x-direction
      character(len=1), intent(in) :: BC2    !< Boundary conditions in y-direction
      character(len=1), intent(in) :: BC3    !< Boundary conditions in z-direction
      ! .. In/Out variables
      integer, dimension(Natom), intent(inout)   :: atype       !< Type of atom
      integer, dimension(Natom), intent(inout)   :: anumb       !< Atom number in cell
      real(dblprec), dimension(3,NA) , intent(inout) :: Bas !< Coordinates for basis atoms
      real(dblprec), dimension(:,:), allocatable, intent(inout) :: coord !< Coordinates of atoms

      ! .. Output variables
      integer, intent(inout) :: metanumb  !<  Flag how to decide metanumb

      integer :: iatom, nuniq, ia, ia_meta, ia_count
      integer :: sx,sy,sz, repX, repY, repZ
      integer :: ibc1,ibc2,ibc3
      integer, dimension(:), allocatable :: nn_max_type
      real(dblprec) :: detmatrix
      real(dblprec), dimension(3,3) :: invmatrix
      real(dblprec), dimension(3) :: coord_shifted
      real(dblprec), dimension(3) :: coord_folded


      call allocate_metanumb(Natom,1)

      if (metanumb==0) then
         NA_meta=NA
         anumb_meta=anumb
         do ia=1,NA
            NA_metalist(ia)=ia
         end do
         return
      else if (do_ralloy.ne.0) then
         write(*,'(1x,a)') "Metanumber functionality not available for random alloys"
         return
      end if

      if(BC1=='P') then
         ibc1=0
      else
         ibc1=1
      end if
      if(BC2=='P') then
         ibc2=0
      else
         ibc2=1
      end if
      if(BC3=='P') then
         ibc3=0
      else
         ibc3=1
      end if

      if(ibc1+ibc2+ibc3==0) then
         NA_meta=NA
         anumb_meta=anumb
         return
      else
         repX=(N1-1)*ibc1+1
         repY=(N2-1)*ibc2+1
         repZ=(N3-1)*ibc3+1
         print *,repX,repY,repZ,repX*repY*repZ
         NA_meta=NA*repX*repY*repZ
      end if

      print *,'NA:', NA,' NA_meta:',NA_meta

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


      ia_count  = 0
      NA_metalist=0
      do iatom=1,natom
         ia = mod(iatom-1,NA)+1
         coord_shifted=coord(:,iatom)-coord(:,anumb(iatom))
         coord_folded=matmul(coord_shifted,invmatrix)
         sx=nint(coord_folded(1))*ibc1
         sy=nint(coord_folded(2))*ibc2
         sz=nint(coord_folded(3))*ibc3
         !print '(2x,i8,3f12.6,5x,3f12.6)', iatom, coord_shifted, coord_folded
         ia_meta = na * (sx + sy*repX + sz*repX*repY) + ia
         !NA_meta=((N1-1)*ibc1+(N2-1)*ibc2+(N3-1)*ibc3+1)*NA
         !print '(38x,5i6)', ia,ia_meta,sx,sy,sz
         !write(999,*) iatom, ia_meta
         if (ia_meta>maxval(anumb_meta)) then
            ia_count=ia_count+1
            NA_metalist(ia_count)=iatom
         end if
         anumb_meta(iatom) = ia_meta
      end do


   end subroutine find_metanumb


end module MetaTypes
