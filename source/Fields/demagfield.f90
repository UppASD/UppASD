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
   real(dblprec), dimension(3) :: demag_mask
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
      !.. Demagnetization mask
      demag_mask = 0.0_dblprec
      if(demag1=='Y') demag_mask(1) = -mu0*mub/demagvol/Natom
      if(demag2=='Y') demag_mask(2) = -mu0*mub/demagvol/Natom
      if(demag3=='Y') demag_mask(3) = -mu0*mub/demagvol/Natom

      !.. Demagnetization field
      do k=1, Mensemble
         demagfld(1,k) = mx(k) * demag_mask(1)
         demagfld(2,k) = my(k) * demag_mask(2)
         demagfld(3,k) = mz(k) * demag_mask(3)
         !!! if (demag1=='Y') then
         !!!    demagfld(1,k)=-mu0*mx(k)/Natom*mub/demagvol
         !!! endif
         !!! if (demag2=='Y') then
         !!!    demagfld(2,k)=-mu0*my(k)/Natom*mub/demagvol
         !!! endif
         !!! if (demag3=='Y') then
         !!!    demagfld(3,k)=-mu0*mz(k)/Natom*mub/demagvol
         !!!    !print *,'DEMAG:', demagfld(3,k)
         !!! endif
      end do
      !print *,' AB Mavg', mx(1), my(1), mz(1)
      !print *,' AB Dema', demagfld(:,1)

   end subroutine calc_demag

   function demag_update_diff(mom_pre, mom_aft) result(db_demag)
      use Constants, only : mub, mu0

      implicit none
      real(dblprec), dimension(3), intent(in) :: mom_pre
      real(dblprec), dimension(3), intent(in) :: mom_aft

      real(dblprec), dimension(3) :: db_demag

      db_demag = dot_product(mom_pre-mom_aft,demag_mask)

   end function demag_update_diff

   function demag_update(mom_in) result(db_demag)
      use Constants, only : mub, mu0

      implicit none
      real(dblprec), dimension(3), intent(in) :: mom_in

      real(dblprec), dimension(3) :: db_demag

      db_demag = -mom_in*demag_mask

   end function demag_update

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
