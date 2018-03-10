!-------------------------------------------------------------------------------
! MODULE: UpdateMoments
!> Routines for updating magnetic moment after time evolution
!> @copyright
!! Copyright (C) 2008-2018 UppASD group
!! This file is distributed under the terms of the
!! GNU General Public License.
!! See http://www.gnu.org/copyleft/gpl.txt
!-------------------------------------------------------------------------------
module UpdateMoments
   use Parameters
   use Profiling

contains

   !-----------------------------------------------------------------------------
   ! SUBROUTINE: moment_update
   !> @brief Wrapper routine for updating the magnetic moments
   !-----------------------------------------------------------------------------
   subroutine moment_update(Natom, Mensemble, mmom, mmom0, mmom2, emom, emom2, emomM, mmomi, mompar, initexc)

      implicit none

      integer, intent(in) :: Natom        !< Number of atoms in system
      integer, intent(in) :: Mensemble    !< Number of ensembles
      real(dblprec), dimension(3,Natom,Mensemble), intent(out) :: emom    !< Current unit moment vector
      real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: emom2    !< Final (or temporary) unit moment vector
      real(dblprec), dimension(3,Natom,Mensemble), intent(out) :: emomM   !< Current magnetic moment vector
      real(dblprec), dimension(Natom,Mensemble), intent(inout) :: mmom    !< Magnitude of magnetic moments
      real(dblprec), dimension(Natom,Mensemble), intent(in), optional :: mmom0 !< Starting magnitude of magnetic moments
      real(dblprec), dimension(Natom,Mensemble), intent(out) :: mmom2 !< Temporary value of magnitude of magnetic moments
      real(dblprec), dimension(Natom,Mensemble), intent(out) :: mmomi !< Inverse of magnitude of magnetic moments
      integer, intent(in), optional :: mompar !< Parametrization of magnetic moment magnitudes (0=no)
      character(len=1), intent(in) :: initexc !< Mode of excitation of initial magnetic moments (I=vacancies, R=two magnon Raman, F=no)

      if (present(mmom0).and.present(mompar)) then
         call calcm(Natom, Mensemble,mmom2, mmom, emom2, mmom0, mompar)
      endif

      call copym(Natom, Mensemble,emom, emom2, emomM, mmom2, mmom, mmomi, initexc)

   end subroutine moment_update

   !> Transfer moments from emom2 to emom and emomM
   subroutine copym(Natom, Mensemble,emom, emom2, emomM, mmom2, mmom, mmomi, initexc)
      !
      implicit none

      integer, intent(in) :: Natom !< Number of atoms in system
      integer, intent(in) :: Mensemble !< Number of ensembles
      real(dblprec), dimension(3,Natom,Mensemble), intent(out) :: emom   !< Current unit moment vector
      real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: emom2  !< Final (or temporary) unit moment vector
      real(dblprec), dimension(3,Natom,Mensemble), intent(out) :: emomM  !< Current magnetic moment vector
      real(dblprec), dimension(Natom,Mensemble), intent(in) :: mmom2 !< Temporary value of magnitude of magnetic moments
      real(dblprec), dimension(Natom,Mensemble), intent(in) :: mmom !< Magnitude of magnetic moments
      real(dblprec), dimension(Natom,Mensemble), intent(out) :: mmomi !< Inverse of magnitude of magnetic moments
      character(len=1), intent(in) :: initexc !< Mode of excitation of initial magnetic moments (I=vacancies, R=two magnon Raman, F=no)

      real(dblprec) :: momtol = 0.000001_dblprec
      integer :: i, j

      ! Transfer moments
      if(initexc /= 'I') then

         !$omp parallel do default(shared) private(i,j) collapse(2)
         do j=1, Mensemble
            do i=1, Natom
               emom(:,i,j)  = emom2(:,i,j)
               emomM(:,i,j) = emom2(:,i,j)*mmom2(i,j)
               mmomi(i,j)   = 1.0d0/mmom(i,j)

            end do
         end do
         !$omp end parallel do

      else

         !$omp parallel do default(shared) private(i,j) collapse(2)
         do i=1, Natom
            do j=1, Mensemble
               emom(1,i,j)=emom2(1,i,j)
               emom(2,i,j)=emom2(2,i,j)
               emom(3,i,j)=emom2(3,i,j)
               emomM(1,i,j)=emom2(1,i,j)*mmom(i,j)
               emomM(2,i,j)=emom2(2,i,j)*mmom(i,j)
               emomM(3,i,j)=emom2(3,i,j)*mmom(i,j)
               if(mmom(i,j) .lt. momtol) then
                  mmomi(i,j)=1.0d0
               else
                  mmomi(i,j)=1.0d0/mmom(i,j)
               end if
            end do
         end do
         !$omp end parallel do

      end if

   end subroutine copym

   !-----------------------------------------------------------------------------
   ! SUBROUTINE: copym
   !> @brief Transfer moments from emom2 to emom and emomM
   !-----------------------------------------------------------------------------
   subroutine calcm(Natom, Mensemble,mmom2, mmom, emom2, mmom0, mompar)
      !

      implicit none

      integer, intent(in) :: Natom !< Number of atoms in system
      integer, intent(in) :: Mensemble !< Number of ensembles
      real(dblprec), dimension(Natom,Mensemble), intent(out) :: mmom2 !< Temporary value of magnitude of magnetic moments
      real(dblprec), dimension(Natom,Mensemble), intent(inout) :: mmom !< Magnitude of magnetic moments
      real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: emom2  !< Final (or temporary) unit moment vector
      real(dblprec), dimension(Natom,Mensemble), optional, intent(in) :: mmom0 !< Starting magnitude of magnetic moments
      integer, intent(in), optional :: mompar !< Parametrization of magnetic moment magnitudes (0=no)


      integer :: i, j

      if(mompar==1) then
         do j=1, Mensemble
            do i=1, Natom

               !M=M0*mz
               mmom2(i,j)=max(mmom0(i,j)*abs(emom2(3,i,j)),1d-4)
            end do
         end do
         mmom=mmom2
      else if(mompar==2) then
         do j=1, Mensemble
            do i=1, Natom

               !M=M0*mz^2
               mmom2(i,j)=max(mmom0(i,j)*(emom2(3,i,j)*emom2(3,i,j)),1.0d-4)
            end do
         end do
         mmom=mmom2
      else
         do j=1, Mensemble
            do i=1, Natom
               !Normal case:
               mmom2(i,j)=mmom(i,j)
            end do
         end do
      end if
   end subroutine calcm

end module UpdateMoments
