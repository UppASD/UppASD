!> Routine for saving moment snapshot to be used as restart file
!> @copyright
!! Copyright (C) 2008-2018 UppASD group
!! This file is distributed under the terms of the
!! GNU General Public License. 
!! See http://www.gnu.org/copyleft/gpl.txt
module Restart
   use Parameters
   use Profiling

   implicit none
   public


contains


   !> Prints magnetic moment to restart file
   subroutine prnrestart(Natom, Mensemble, simid, mstep, emom, mmom)
      !
      !.. Implicit declarations
      implicit none

      integer, intent(in) :: Natom !< Number of atoms in system
      integer, intent(in) :: Mensemble !< Number of ensembles 
      character(len=8), intent(in) :: simid !< Name of simulation 
      integer, intent(in) :: mstep !< Current simulation step
      real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: emom   !< Current unit moment vector
      real(dblprec), dimension(Natom,Mensemble), intent(in) :: mmom !< Magnitude of magnetic moments

      integer :: i, j
      character(len=30) :: filn

      !.. Executable statements
      write (filn,'(''restart.'',a8,''.out'')') simid
      open(ofileno, file=filn)
      write (ofileno,10002) mstep

      do i=1, Mensemble
         do j=1, Natom
            write (ofileno,10003) i, j, mmom(j,i), emom(1,j,i), emom(2,j,i), emom(3,j,i)
         end do
      end do

      close(ofileno)

      10002 format (i8)
      10003 format (i8,i8,2x,es16.8,es16.8,es16.8,es16.8)

   end subroutine prnrestart


end module Restart
