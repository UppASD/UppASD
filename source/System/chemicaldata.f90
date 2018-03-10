!> Data for treating random alloy systems
!! @todo Consider moving allocation routine
!> @copyright
!! Copyright (C) 2008-2018 UppASD group
!! This file is distributed under the terms of the
!! GNU General Public License. 
!! See http://www.gnu.org/copyleft/gpl.txt
module ChemicalData
   use Parameters
   use Profiling
   !
   implicit none
   !
   integer, dimension(:), allocatable :: achtype !< Chemical type of atoms (full list)
   integer, dimension(:), allocatable :: acellnumb !< List for translating atom no. in full cell to actual cell
   integer, dimension(:), allocatable :: acellnumbrev !< List for translating atom no. in actual cell to full cell
   integer, dimension(:), allocatable :: achem_ch !< Chemical type of atoms (reduced list) (achem_ch(i)=achtype(acellnumbrev(i)))
   integer, dimension(:), allocatable :: atype_ch !< Actual site of atom for dilute system 
   integer, dimension(:), allocatable :: asite_ch !< Actual site of atom for dilute system 
   real(dblprec) :: chconceff !< Effective chemical concentration

   !
contains


   !> Allocate arrays for random alloys
   subroutine allocate_chemicaldata(Natom,flag)
      implicit none

      integer, optional, intent(in) :: Natom !< Number of atoms in system
      integer, intent(in) :: flag !< Allocate or deallocate (1/-1)

      integer :: i_all, i_stat

      if(flag>0) then
         allocate(achtype(Natom),stat=i_stat)
         call memocc(i_stat,product(shape(achtype))*kind(achtype),'achtype','allocate_chemicaldata')
         achtype=0
         allocate(acellnumb(Natom),stat=i_stat)
         call memocc(i_stat,product(shape(acellnumb))*kind(acellnumb),'acellnumb','allocate_chemicaldata')
         acellnumb=0
         allocate(acellnumbrev(Natom),stat=i_stat)
         call memocc(i_stat,product(shape(acellnumbrev))*kind(acellnumbrev),'acellnumbrev','allocate_chemicaldata')
         acellnumbrev=0
         allocate(achem_ch(Natom),stat=i_stat)
         call memocc(i_stat,product(shape(achem_ch))*kind(achem_ch),'achem_ch','allocate_chemicaldata')
         achem_ch=0
         allocate(asite_ch(Natom),stat=i_stat)
         call memocc(i_stat,product(shape(asite_ch))*kind(asite_ch),'asite_ch','allocate_chemicaldata')
         asite_ch=0
         allocate(atype_ch(Natom),stat=i_stat)
         call memocc(i_stat,product(shape(atype_ch))*kind(atype_ch),'atype_ch','allocate_chemicaldata')
         atype_ch=0
      else
         i_all=-product(shape(achtype))*kind(achtype)
         deallocate(achtype,stat=i_stat)
         call memocc(i_stat,i_all,'achtype','allocate_chemicaldata')
         i_all=-product(shape(acellnumb))*kind(acellnumb)
         deallocate(acellnumb,stat=i_stat)
         call memocc(i_stat,i_all,'acellnumb','allocate_chemicaldata')
         i_all=-product(shape(acellnumbrev))*kind(acellnumbrev)
         deallocate(acellnumbrev,stat=i_stat)
         call memocc(i_stat,i_all,'acellnumbrev','allocate_chemicaldata')
         i_all=-product(shape(achem_ch))*kind(achem_ch)
         deallocate(achem_ch,stat=i_stat)
         call memocc(i_stat,i_all,'achem_ch','allocate_chemicaldata')
         i_all=-product(shape(asite_ch))*kind(asite_ch)
         deallocate(asite_ch,stat=i_stat)
         call memocc(i_stat,i_all,'asite_ch','allocate_chemicaldata')
         i_all=-product(shape(atype_ch))*kind(atype_ch)
         deallocate(atype_ch,stat=i_stat)
         call memocc(i_stat,i_all,'atype_ch','allocate_chemicaldata')
      end if

   end subroutine allocate_chemicaldata


end module ChemicalData

