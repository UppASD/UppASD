!-------------------------------------------------------------------------------
! MODULE: CalculateFields
!> @brief Routines for calculate external and average fields
!> @copyright
!> GNU Public License.
!-------------------------------------------------------------------------------
module CalculateFields

   use Profiling
   use Parameters

   implicit none

   public

contains

   !----------------------------------------------------------------------------
   ! SUBROUTINE: calc_external_fields
   !> @brief Calculate the total static external field, including global and site
   !> dependent contributions
   !----------------------------------------------------------------------------
   subroutine calc_external_fields(Natom,Mensemble,NA,hfield,anumb,external_fields, &
      do_bpulse,sitefld,sitenatomfld)
      !
      implicit none

      integer, intent(in) :: NA  !< Number of atoms in one cell
      integer, intent(in) :: Natom !< Number of atoms in system
      integer, intent(in) :: Mensemble !< Number of ensembles
      integer, intent(in) :: do_bpulse  !< Add magnetic field pulse (0=no, 1-4 for different shapes)
      integer, dimension(Natom), intent(in) :: anumb !< Atom number in cell
      real(dblprec), dimension(3), intent(in) :: hfield !< Applied magnetic field
      real(dblprec), dimension(3,Natom,Mensemble), intent(inout) :: external_fields !< External magnetic field
      real(dblprec), dimension(:,:), intent(in), optional :: sitefld !< Site dependent applied field
      real(dblprec), dimension(:,:), intent(in), optional :: sitenatomfld !< Site dependent applied field for Natom

      !
      integer :: i,k
      real(dblprec), dimension(3,Mensemble) :: global_field
      !
      ! Initialize global field
      global_field=0.0_dblprec
      external_fields=0.0_dblprec
      !
      ! First find global contribution to external field
      !
      ! Static field
      do k=1,Mensemble
         global_field(1,k) = global_field(1,k)+hfield(1)
         global_field(2,k) = global_field(2,k)+hfield(2)
         global_field(3,k) = global_field(3,k)+hfield(3)
      end do
      !
      ! Then add site-dependent field to global field (This only works for atom types)
      if (do_bpulse==5) then
         do k=1,Mensemble
            do i=1,Natom
               external_fields(1,i,k) = external_fields(1,i,k) + global_field(1,k)+sitefld(1,anumb(i))
               external_fields(2,i,k) = external_fields(2,i,k) + global_field(2,k)+sitefld(2,anumb(i))
               external_fields(3,i,k) = external_fields(3,i,k) + global_field(3,k)+sitefld(3,anumb(i))
            end do
         end do
         ! Then add site-dependent field (Site dependent that is a Natom field)
      else if (do_bpulse==6) then
         do k=1,Mensemble
            do i=1,Natom
               external_fields(1,i,k) = external_fields(1,i,k)+ global_field(1,k)+sitenatomfld(1,i)
               external_fields(2,i,k) = external_fields(2,i,k)+ global_field(2,k)+sitenatomfld(2,i)
               external_fields(3,i,k) = external_fields(3,i,k)+ global_field(3,k)+sitenatomfld(3,i)
            end do
         end do
      else
         do k=1,Mensemble
            do i=1,Natom
               external_fields(1,i,k) = external_fields(1,i,k) + global_field(1,k)
               external_fields(2,i,k) = external_fields(2,i,k) + global_field(2,k)
               external_fields(3,i,k) = external_fields(3,i,k) + global_field(3,k)
            end do
         end do
      end if
      !
   end subroutine calc_external_fields

   !-----------------------------------------------------------------------------
   ! SUBROUTINE: calc_external_time_fields
   !> @brief Performs the calculation of the time-dependent fields
   !> @details This wrapper includes all time dependent fields acting over the sample
   !> it deals with both frequency dependent fields and with fields that traverse
   !> the sample.
   !> @author Jonathan Chico
   !-----------------------------------------------------------------------------
   subroutine calc_external_time_fields(Natom,Mensemble,time_external_fields,       &
      do_bpulse,demag,mwf,mwf_gauss_spatial,do_gauss,mwf_gauss,mov_gauss,           &
      mwf_mov_gauss,bpulsefield,demagfld,mwffield,gauss_mwffield,site_mwffield,     &
      gauss_spatial_site_mwffield,gauss_spatial,gauss_site_mwffield,                &
      mov_gauss_spatial,mwf_mov_gauss_spatial,mov_circle,mwf_mov_circle,            &
      mov_circle_spatial,mwf_mov_circle_spatial,mov_square,mwf_mov_square,          &
      mov_square_spatial,mwf_mov_square_spatial)

      implicit none

      integer, intent(in) :: Natom !< Number of atoms in system
      integer, intent(in) :: Mensemble !< Number of ensembles
      integer, intent(in) :: do_bpulse  !< Add magnetic field pulse (0=no, 1-4 for different shapes)

      character(len=1), intent(in) :: mwf       !< Add monochromatic microwave field (Y/P/S/W/N)
      character(len=1), intent(in) :: demag     !< Model demagnetization field (Y/N)
      character(len=1), intent(in) :: do_gauss  !< Add the Static gaussian shaped field
      character(len=1), intent(in) :: mwf_gauss !< Add frequency broadened microwave field (Y/P/S/W/N)
      character(len=1), intent(in) :: mov_gauss !< Add the moving static gaussian shaped field (Y/P/N)
      character(len=1), intent(in) :: mov_circle   !< Add the moving static circular shaped field (Y/P/N)
      character(len=1), intent(in) :: mov_square   !< Add the moving static cubic shaped field (Y/P/N)
      character(len=1), intent(in) :: mwf_mov_gauss      !< Add the moving microwave gaussian shaped field (Y/P/N)
      character(len=1), intent(in) :: mwf_mov_circle     !< Add the moving microwave circular shaped field (Y/P/N)
      character(len=1), intent(in) :: mwf_mov_square     !< Add the moving microwave cubic shaped field (Y/P/N)
      character(len=1), intent(in) :: mwf_gauss_spatial  !< Add the frequency broadened gaussian shaped microwave field (Y/P/N)

      real(dblprec), dimension(3), intent(in) :: mwffield         !< Current microwave field
      real(dblprec), dimension(3), intent(in) :: bpulsefield      !< Current applied magnetic field from pulse
      real(dblprec), dimension(3), intent(in) :: gauss_mwffield   !< Current gaussian broadened microwave field
      real(dblprec), dimension(3,Mensemble), intent(in)  ::  demagfld      !< Demagnetization field
      real(dblprec), dimension(3,Natom), intent(in)      :: site_mwffield  !< Site and time dependent applied field for Natom
      real(dblprec), dimension(3,Natom), intent(in)      :: gauss_spatial  !< Array with the static gaussian shaped field
      real(dblprec), dimension(3,Natom), intent(in)      :: mov_gauss_spatial    !< Array with the moving gaussian shaped static field
      real(dblprec), dimension(3,Natom), intent(in)      :: mov_square_spatial   !< Array with the moving cubic shaped static field
      real(dblprec), dimension(3,Natom), intent(in)      :: mov_circle_spatial   !< Array with the moving circular shaped static field
      real(dblprec), dimension(3,Natom), intent(in)      :: gauss_site_mwffield  !< Array with the gaussian broadened microwave field
      real(dblprec), dimension(3,Natom), intent(in)      :: mwf_mov_gauss_spatial   !< Array with the moving microwave gaussian shaped field
      real(dblprec), dimension(3,Natom), intent(in)      :: mwf_mov_square_spatial  !< Array with the moving microwave cubic shaped field
      real(dblprec), dimension(3,Natom), intent(in)      :: mwf_mov_circle_spatial  !< Array with the moving microwave circular shaped field
      real(dblprec), dimension(3,Natom), intent(in)      :: gauss_spatial_site_mwffield   !< Array with the gaussian broadened spatially gaussian microwave field

      real(dblprec), dimension(3,Natom,Mensemble), intent(out) :: time_external_fields !< External magnetic field

      integer :: i,k
      real(dblprec), dimension(3,Mensemble) :: global_field
      !
      ! Initialize global field
      global_field=0.0_dblprec
      time_external_fields=0.0_dblprec

      ! Pulsed field
      if (do_bpulse>0.and.do_bpulse<5) then
         do k=1,Mensemble
            global_field(1,k) = global_field(1,k)+bpulsefield(1)
            global_field(2,k) = global_field(2,k)+bpulsefield(2)
            global_field(3,k) = global_field(3,k)+bpulsefield(3)
         end do
      end if
      !
      ! Demagnetization field
      if (demag=='Y') then
         do k=1,Mensemble
            global_field(1,k) = global_field(1,k)+demagfld(1,k)
            global_field(2,k) = global_field(2,k)+demagfld(2,k)
            global_field(3,k) = global_field(3,k)+demagfld(3,k)
         end do
      endif
      !

      ! Microwave field
      if (mwf=='Y'.or.mwf=='P'.or.mwf=='I') then
         do k=1,Mensemble
            global_field(1,k) = global_field(1,k)+mwffield(1)
            global_field(2,k) = global_field(2,k)+mwffield(2)
            global_field(3,k) = global_field(3,k)+mwffield(3)
         end do
      endif
      !
      ! Global gaussian microwave field
      if (mwf_gauss=='Y'.or.mwf_gauss=='P') then
         do k=1,Mensemble
            global_field(1,k) = global_field(1,k)+gauss_mwffield(1)
            global_field(2,k) = global_field(2,k)+gauss_mwffield(2)
            global_field(3,k) = global_field(3,k)+gauss_mwffield(3)
         end do
      endif

      ! Site dependent microwave pulse field
      if (mwf=='S'.or.mwf=='W') then
         do k=1,Mensemble
            do i=1,Natom
               time_external_fields(1,i,k) = time_external_fields(1,i,k) + global_field(1,k)+site_mwffield(1,i)
               time_external_fields(2,i,k) = time_external_fields(2,i,k) + global_field(2,k)+site_mwffield(2,i)
               time_external_fields(3,i,k) = time_external_fields(3,i,k) + global_field(3,k)+site_mwffield(3,i)
            end do
         end do
      end if
      ! Site dependent microwave gaussian field
      if(mwf_gauss=='S'.or.mwf_gauss=='W') then
         do k=1,Mensemble
            do i=1,Natom
               time_external_fields(1,i,k) = time_external_fields(1,i,k) + global_field(1,k)+gauss_site_mwffield(1,i)
               time_external_fields(2,i,k) = time_external_fields(2,i,k) + global_field(2,k)+gauss_site_mwffield(2,i)
               time_external_fields(3,i,k) = time_external_fields(3,i,k) + global_field(3,k)+gauss_site_mwffield(3,i)
            end do
         end do
      endif

      if (do_gauss=='Y'.or.do_gauss=='P') then
         do k=1, Mensemble
            do i=1, Natom
               time_external_fields(1,i,k) = time_external_fields(1,i,k) + global_field(1,k)+gauss_spatial(1,i)
               time_external_fields(2,i,k) = time_external_fields(2,i,k) + global_field(2,k)+gauss_spatial(2,i)
               time_external_fields(3,i,k) = time_external_fields(3,i,k) + global_field(3,k)+gauss_spatial(3,i)
            enddo
         enddo
      endif

      if (mwf_gauss_spatial=='Y'.or.mwf_gauss_spatial=='P') then
         do k=1,Mensemble
            do i=1,Natom
               time_external_fields(1,i,k) = time_external_fields(1,i,k) + global_field(1,k)+gauss_spatial_site_mwffield(1,i)
               time_external_fields(2,i,k) = time_external_fields(2,i,k) + global_field(2,k)+gauss_spatial_site_mwffield(2,i)
               time_external_fields(3,i,k) = time_external_fields(3,i,k) + global_field(3,k)+gauss_spatial_site_mwffield(3,i)
            end do
         end do
      endif

      if (mov_gauss=='Y'.or.mov_gauss=='P') then
         do k=1,Mensemble
            do i=1,Natom
               time_external_fields(1,i,k) = time_external_fields(1,i,k) + global_field(1,k)+mov_gauss_spatial(1,i)
               time_external_fields(2,i,k) = time_external_fields(2,i,k) + global_field(2,k)+mov_gauss_spatial(2,i)
               time_external_fields(3,i,k) = time_external_fields(3,i,k) + global_field(3,k)+mov_gauss_spatial(3,i)
            end do
         end do
      endif

      if (mwf_mov_gauss=='Y'.or.mwf_mov_gauss=='P') then
         do k=1,Mensemble
            do i=1,Natom
               time_external_fields(1,i,k) = time_external_fields(1,i,k) + global_field(1,k)+mwf_mov_gauss_spatial(1,i)
               time_external_fields(2,i,k) = time_external_fields(2,i,k) + global_field(2,k)+mwf_mov_gauss_spatial(2,i)
               time_external_fields(3,i,k) = time_external_fields(3,i,k) + global_field(3,k)+mwf_mov_gauss_spatial(3,i)
            end do
         end do
      endif

      ! Calculate circular shaped moving field
      if (mov_circle=='Y'.or.mov_circle=='P') then
         do k=1,Mensemble
            do i=1,Natom
               time_external_fields(1,i,k) = time_external_fields(1,i,k) + global_field(1,k)+mov_circle_spatial(1,i)
               time_external_fields(2,i,k) = time_external_fields(2,i,k) + global_field(2,k)+mov_circle_spatial(2,i)
               time_external_fields(3,i,k) = time_external_fields(3,i,k) + global_field(3,k)+mov_circle_spatial(3,i)
            end do
         end do
      endif
      ! Calculate circular shaped microwave field
      if (mwf_mov_circle=='Y'.or.mwf_mov_circle=='P') then
         do k=1,Mensemble
            do i=1,Natom
               time_external_fields(1,i,k) = time_external_fields(1,i,k) + global_field(1,k)+mwf_mov_circle_spatial(1,i)
               time_external_fields(2,i,k) = time_external_fields(2,i,k) + global_field(2,k)+mwf_mov_circle_spatial(2,i)
               time_external_fields(3,i,k) = time_external_fields(3,i,k) + global_field(3,k)+mwf_mov_circle_spatial(3,i)
            end do
         end do
      endif

      ! Calculate cubic shaped moving field
      if (mov_square=='Y'.or.mov_square=='P') then
         do k=1,Mensemble
            do i=1,Natom
               time_external_fields(1,i,k) = time_external_fields(1,i,k) + global_field(1,k)+mov_square_spatial(1,i)
               time_external_fields(2,i,k) = time_external_fields(2,i,k) + global_field(2,k)+mov_square_spatial(2,i)
               time_external_fields(3,i,k) = time_external_fields(3,i,k) + global_field(3,k)+mov_square_spatial(3,i)
            end do
         end do
      endif
      ! Calculate cubic shaped microwave field
      if (mwf_mov_square=='Y'.or.mwf_mov_square=='P') then
         do k=1,Mensemble
            do i=1,Natom
               time_external_fields(1,i,k) = time_external_fields(1,i,k) + global_field(1,k)+mwf_mov_square_spatial(1,i)
               time_external_fields(2,i,k) = time_external_fields(2,i,k) + global_field(2,k)+mwf_mov_square_spatial(2,i)
               time_external_fields(3,i,k) = time_external_fields(3,i,k) + global_field(3,k)+mwf_mov_square_spatial(3,i)
            end do
         end do
      endif

      ! Add the fields
      do k=1,Mensemble
         do i=1,Natom
            time_external_fields(1,i,k) = time_external_fields(1,i,k) + global_field(1,k)
            time_external_fields(2,i,k) = time_external_fields(2,i,k) + global_field(2,k)
            time_external_fields(3,i,k) = time_external_fields(3,i,k) + global_field(3,k)
         end do
      end do

   end subroutine calc_external_time_fields

   !-----------------------------------------------------------------------------
   ! SUBROUTINE: calc_averagefield
   !> @brief Calculate average field for use with multiple heat baths
   !-----------------------------------------------------------------------------
   subroutine calc_averagefield(Natom,Mensemble,beff1,beff2,field1,field2)
      !
      !.. Implicit declarations
      implicit none

      integer, intent(in) :: Natom !< Number of atoms in system
      integer, intent(in) :: Mensemble !< Number of ensembles
      real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: beff1 !< Internal effective field from application of Hamiltonian
      real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: beff2 !< External field from application of Hamiltonian
      real(dblprec), dimension(3,Mensemble), intent(out) :: field1 !< Average internal effective field
      real(dblprec), dimension(3,Mensemble), intent(out) :: field2 !< Average external effective field

      !.. Local scalars
      integer :: i, j

      !.. Executable statements

      !$omp parallel do default(shared) &
      !$omp  private(i,j)
      do j=1, Mensemble
         field1(1,j)=0.0_dblprec
         field1(2,j)=0.0_dblprec
         field1(3,j)=0.0_dblprec
         field2(1,j)=0.0_dblprec
         field2(2,j)=0.0_dblprec
         field2(3,j)=0.0_dblprec
         do i=1, Natom
            field1(1:3,j)=field1(1:3,j)+beff1(1:3,i,j)
            field2(1:3,j)=field2(1:3,j)+beff2(1:3,i,j)
         end do
         field1(1:3,j)=field1(1:3,j)/Natom
         field2(1:3,j)=field2(1:3,j)/Natom
      end do
      !$omp end parallel do
   end subroutine calc_averagefield

end module calculatefields
