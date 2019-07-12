!-------------------------------------------------------------------------------
! MODULE DemagField
!> @brief Data and routines for calculating demagnetizing field
!> @author Johan Hellsvik
!> @copyright
!> GNU Public License.
!-------------------------------------------------------------------------------
module DemagField
   use Parameters
   use Profiling
   !
   implicit none
   !
   real(dblprec), dimension(:,:), allocatable :: demagfld !< Demagnetization field
   !
   public


contains

   !-----------------------------------------------------------------------------
   ! SUBROUTINE: calc_demag
   !> @brief Calculate demagnetization field
   !-----------------------------------------------------------------------------
   subroutine calc_demag(Natom, Mensemble, demag1, demag2, demag3, demagvol, emomM)
      use Constants, only : mub, mu0

      !.. Implicit declatrations
      implicit none

      !.. Formal arguments
      integer, intent(in) :: Natom !< Number of atoms in system
      integer, intent(in) :: Mensemble !< Number of ensembles
      real(dblprec), dimension(3,Natom,Mensemble),intent(in) ::  emomM !< Current magnetic moment vector
      character(len=1),intent(in) :: demag1 !< Model demagnetization field in x-direction (Y/N)
      character(len=1),intent(in) :: demag2 !< Model demagnetization field in y-direction (Y/N)
      character(len=1),intent(in) :: demag3 !< Model demagnetization field in z-direction (Y/N)
      real(dblprec),intent(in) :: demagvol !< Volume for modelling demagnetization field

      !.. Scalar variables
      integer :: i, k
      real(dblprec), dimension(Mensemble) ::  mx, my, mz

      !.. Executable statements
      mx=0.0_dblprec
      my=0.0_dblprec
      mz=0.0_dblprec
      !.. Sum over moments
      do k=1, Mensemble
         do i=1, Natom
            mx(k) = mx(k) + emomM(1,i,k)
            my(k) = my(k) + emomM(2,i,k)
            mz(k) = mz(k) + emomM(3,i,k)
         end do
      end do
      !.. Demagnetization field
      do k=1, Mensemble
         if (demag1=='Y') then
            demagfld(1,k)=-mu0*mx(k)/Natom*mub/demagvol
         endif
         if (demag2=='Y') then
            demagfld(2,k)=-mu0*my(k)/Natom*mub/demagvol
         endif
         if (demag3=='Y') then
            demagfld(3,k)=-mu0*mz(k)/Natom*mub/demagvol
         endif
      end do
   end subroutine calc_demag


   !-----------------------------------------------------------------------------
   ! SUBROUTINE: allocate_demag
   !> @brief Allocate array for demagnetization field
   !-----------------------------------------------------------------------------
   subroutine allocate_demag(Mensemble)
      !
      implicit none
      !
      integer, intent(in) :: Mensemble !< Number of ensembles
      !
      integer :: i_stat
      !
      allocate(demagfld(3,Mensemble),stat=i_stat)
      call memocc(i_stat,product(shape(demagfld))*kind(demagfld),'demagfld','read_demag')
      demagfld=0.0_dblprec
      !
   end subroutine allocate_demag

end module DemagField
