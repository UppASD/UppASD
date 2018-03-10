!> Data describing the simulated system
!> @copyright
!! Copyright (C) 2008-2018 UppASD group
!! This file is distributed under the terms of the
!! GNU General Public License. 
!! See http://www.gnu.org/copyleft/gpl.txt
module SystemData
   use Parameters
   use Profiling
   !
   implicit none
   !
   !From setup
   real(dblprec), dimension(:), allocatable :: Landeg  !< Gyromagnetic ratio
   integer, dimension(:), allocatable :: anumb !< Atom number in cell
   integer, dimension(:), allocatable :: atype !< Type of atom 
   real(dblprec), dimension(:,:), allocatable :: coord !< Coordinates of atoms
   !
   integer, parameter ::  limit_no_shells=350 !< Limit of number of shells for exchange interactions
   integer :: nn_dm_max=54 !< Limit of number of shells for DM interactions

   public


contains


   subroutine allocate_systemdata(Natom, flag)

      implicit none

      integer, intent(in),optional :: Natom !< Number of atoms in system
      integer, intent(in) :: flag  !< Allocate or deallocate (1/-1)

      integer :: i_stat, i_all

      if(flag>0) then
         allocate(anumb(Natom),stat=i_stat)
         call memocc(i_stat,product(shape(anumb))*kind(anumb),'anumb','read_input1')
         allocate(atype(Natom),stat=i_stat)
         call memocc(i_stat,product(shape(atype))*kind(atype),'atype','read_input1')
         allocate(Landeg(Natom),stat=i_stat)
         call memocc(i_stat,product(shape(Landeg))*kind(Landeg),'Landeg','read_input1')
      else
         i_all=-product(shape(anumb))*kind(anumb)
         deallocate(anumb,stat=i_stat)
         call memocc(i_stat,i_all,'anumb','allocate_systemdata')
         i_all=-product(shape(atype))*kind(atype)
         deallocate(atype,stat=i_stat)
         call memocc(i_stat,i_all,'atype','allocate_systemdata')
         i_all=-product(shape(Landeg))*kind(Landeg)
         deallocate(Landeg,stat=i_stat)
         call memocc(i_stat,i_all,'Landeg','allocate_systemdata')
      end if

   end subroutine allocate_systemdata


end module SystemData


