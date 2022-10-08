!> Routine for saving displacments and velocities snapshot to be used as lattice restart file
module LatticeRestart
   use Parameters
   use Profiling

   implicit none
   public


contains


   !> Prints displacements and velocities to restart file
   subroutine prnlattrestart(simid, Natom, Mensemble, mstep, uvec, vvec)
      !
      !.. Implicit declarations
      implicit none

      character(len=8), intent(in) :: simid !< Name of simulation 
      integer, intent(in) :: Natom !< Number of atoms in system
      integer, intent(in) :: Mensemble !< Number of ensembles 
      integer, intent(in) :: mstep !< Current simulation step
      real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: uvec   !< Current ionic displacement
      real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: vvec   !< Current ionic velocity

      !.. Local variables
      integer :: i, j
      character(len=30) :: filn

      !.. Executable statements
      write (filn,'(''lattrestart.'',a,''.out'')') trim(simid)
      open(9, file=filn)
      write (9,10002) mstep

      do i=1, Mensemble
         do j=1, Natom
            write (9,10003) i, j, uvec(1,j,i), uvec(2,j,i), uvec(3,j,i), vvec(1,j,i), vvec(2,j,i), vvec(3,j,i)
         end do
      end do

      close(9)

      10002 format (i8)
      10003 format (i8,i8,2x,8es16.8)

   end subroutine prnlattrestart


end module LatticeRestart
