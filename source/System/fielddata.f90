!> Data for storing external and effective fields
!> @copyright
!! Copyright (C) 2008-2018 UppASD group
!! This file is distributed under the terms of the
!! GNU General Public License. 
!! See http://www.gnu.org/copyleft/gpl.txt
module FieldData
   use Parameters
   use Profiling
   !
   implicit none
   !
   real(dblprec), dimension(:,:), allocatable :: sitefld !< Site dependent applied field
   real(dblprec), dimension(:,:,:), allocatable :: beff !< Total effective field from application of Hamiltonian
   real(dblprec), dimension(:,:,:), allocatable :: beff1 !< Internal effective field from application of Hamiltonian
   real(dblprec), dimension(:,:,:), allocatable :: beff2 !< External field from application of Hamiltonian
   real(dblprec), dimension(:,:,:), allocatable :: beff3 !< Internal effective bfield from mixed spin-lattice Hamiltonians
   real(dblprec), dimension(:,:,:), allocatable :: b2eff !< Temporary storage of magnetic field
   real(dblprec), dimension(:,:,:), allocatable :: external_field !< External magnetic field
   real(dblprec), dimension(:,:,:), allocatable :: time_external_field
   real(dblprec), dimension(:,:,:), allocatable :: thermal_field
   !
   real(dblprec), dimension(:,:), allocatable :: field1 !< Average internal effective field
   real(dblprec), dimension(:,:), allocatable :: field2 !< Average external effective field


contains
   !> Allocate/deallocate arrays for external and effective fields
   subroutine allocate_fields(Natom,Mensemble,flag)

      implicit none

      integer, intent(in), optional :: Natom !< Number of atoms in system
      integer, intent(in), optional :: Mensemble !< Number of ensembles
      integer, intent(in) :: flag  !< Allocate or deallocate (1/-1)

      integer :: i_all, i_stat

      ! If flag > 0 allocate arrays
      if(flag>0) then
         allocate(beff(3,Natom,Mensemble),stat=i_stat)
         call memocc(i_stat,product(shape(beff))*kind(beff),'beff','allocate_fields')
         beff=0.0D0
         allocate(beff1(3,Natom,Mensemble),stat=i_stat)
         call memocc(i_stat,product(shape(beff1))*kind(beff1),'beff1','allocate_fields')
         beff1=0.0D0
         allocate(beff2(3,Natom,Mensemble),stat=i_stat)
         call memocc(i_stat,product(shape(beff2))*kind(beff2),'beff2','allocate_fields')
         beff2=0.0D0
         allocate(beff3(3,Natom,Mensemble),stat=i_stat)
         call memocc(i_stat,product(shape(beff3))*kind(beff3),'beff3','allocate_fields')
         beff3=0.0D0
         allocate(b2eff(3,Natom,Mensemble),stat=i_stat)
         call memocc(i_stat,product(shape(b2eff))*kind(b2eff),'b2eff','allocate_fields')
         b2eff=0.0D0
         allocate(external_field(3,Natom,Mensemble),stat=i_stat)
         call memocc(i_stat,product(shape(external_field))*kind(external_field),'external_field','allocate_fields')
         external_field=0.0D0
         allocate(field1(3,Mensemble),stat=i_stat)
         call memocc(i_stat,product(shape(field1))*kind(field1),'field1','allocate_fields')
         field1=0.0D0
         allocate(field2(3,Mensemble),stat=i_stat)
         call memocc(i_stat,product(shape(field2))*kind(field2),'field2','allocate_fields')
         field2=0.0D0
         allocate(time_external_field(3,Natom,Mensemble),stat=i_stat)
         call memocc(i_stat,product(shape(time_external_field))*kind(time_external_field),'time_external_field','allocate_fields')
         time_external_field=0.0D0
         allocate(thermal_field(3,Natom,Mensemble),stat=i_stat)
         call memocc(i_stat,product(shape(thermal_field))*kind(thermal_field),'thermal_field','allocate_fields')
         thermal_field=0.0D0
      else
         ! Otherwise deallocate arrays
         i_all=-product(shape(beff))*kind(beff)
         deallocate(beff,stat=i_stat)
         call memocc(i_stat,i_all,'beff','allocate_fields')
         i_all=-product(shape(beff1))*kind(beff1)
         deallocate(beff1,stat=i_stat)
         call memocc(i_stat,i_all,'beff1','allocate_fields')
         i_all=-product(shape(beff2))*kind(beff2)
         deallocate(beff2,stat=i_stat)
         call memocc(i_stat,i_all,'beff2','allocate_fields')
         i_all=-product(shape(beff3))*kind(beff3)
         deallocate(beff3,stat=i_stat)
         call memocc(i_stat,i_all,'beff3','allocate_fields')
         i_all=-product(shape(b2eff))*kind(b2eff)
         deallocate(b2eff,stat=i_stat)
         call memocc(i_stat,i_all,'b2eff','allocate_fields')
         i_all=-product(shape(external_field))*kind(external_field)
         deallocate(external_field,stat=i_stat)
         call memocc(i_stat,i_all,'external_field','allocate_fields')
         i_all=-product(shape(field1))*kind(field1)
         deallocate(field1,stat=i_stat)
         call memocc(i_stat,i_all,'field1','allocate_fields')
         i_all=-product(shape(field2))*kind(field2)
         deallocate(field2,stat=i_stat)
         call memocc(i_stat,i_all,'field2','allocate_fields')
         !
         if (allocated(time_external_field)) then
            i_all=-product(shape(time_external_field))*kind(time_external_field)
            deallocate(time_external_field,stat=i_stat)
            call memocc(i_stat,i_all,'time_external_field','allocate_fields')
         endif
         !
         if (allocated(thermal_field)) then
            i_all=-product(shape(thermal_field))*kind(thermal_field)
            deallocate(thermal_field,stat=i_stat)
            call memocc(i_stat,i_all,'thermal_field','allocate_fields')
         endif

      end if
   end subroutine allocate_fields

   !> Read in local fields from file
   subroutine read_local_field(NA,locfieldfile)
      !
      implicit none
      !
      integer, intent(in) :: NA  !< Number of atoms in one cell
      character(LEN=35), intent(in) :: locfieldfile !< File name for local applied field info
      !
      integer :: i, j, i_stat, i_all
      !
      print *,'-->', locfieldfile
      open(ifileno,file=locfieldfile)
      if(NA>0) then
         allocate(sitefld(3,NA),stat=i_stat)
         call memocc(i_stat,product(shape(sitefld))*kind(sitefld),'sitefld','read_local_field')
      else
         i_all=-product(shape(sitefld))*kind(sitefld)
         deallocate(sitefld,stat=i_stat)
         call memocc(i_stat,i_all,'sitefld','allocate_fields')
         return
      end if
      do i=1,NA
         read(ifileno,*) (sitefld(j,i),j=1,3)
      end do
      close(ifileno)
      return
   end subroutine read_local_field

   !> Determine stuff relating to field pulse
   logical function allocation_field_time(mwf,do_gauss,mwf_gauss,mov_gauss,mwf_mov_gauss,mwf_gauss_spatial,&
         mov_circle,mwf_mov_circle,mov_square,mwf_mov_square)

      character(len=1), intent(in) :: mwf
      character(len=1), intent(in) :: do_gauss
      character(len=1), intent(in) :: mwf_gauss
      character(len=1), intent(in) :: mov_gauss
      character(len=1), intent(in) :: mov_circle
      character(len=1), intent(in) :: mov_square
      character(len=1), intent(in) :: mwf_mov_gauss
      character(len=1), intent(in) :: mwf_mov_circle
      character(len=1), intent(in) :: mwf_mov_square
      character(len=1), intent(in) :: mwf_gauss_spatial

      if (mwf=='Y'.or.mwf=='P'.or.mwf=='S'.or.mwf=='W'.or.mwf=='I') then
         allocation_field_time=.true.

      else if (do_gauss=='P'.or.do_gauss=='Y') then
         allocation_field_time=.true.

      else if (mwf_gauss=='Y'.or.mwf_gauss=='P'.or.mwf_gauss=='S'.or.mwf_gauss=='W') then
         allocation_field_time=.true.

      else if (mov_gauss=='Y'.or.mov_gauss=='P') then
         allocation_field_time=.true.

      else if (mwf_mov_gauss=='Y'.or.mwf_mov_gauss=='P') then
         allocation_field_time=.true.

      else if (mov_circle=='Y'.or.mov_circle=='P') then
         allocation_field_time=.true.

      else if (mwf_mov_circle=='Y'.or.mwf_mov_circle=='P') then
         allocation_field_time=.true.

      else if (mov_square=='Y'.or.mov_square=='P') then
         allocation_field_time=.true.

      else if (mwf_mov_square=='Y'.or.mwf_mov_square=='P') then
         allocation_field_time=.true.

      else if (mwf_gauss_spatial=='Y'.or.mwf_gauss_spatial=='P') then
         allocation_field_time=.true.

      else
         allocation_field_time=.false.

      endif

   end function  allocation_field_time

end module FieldData
