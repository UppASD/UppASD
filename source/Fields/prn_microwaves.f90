!> Data and routines necessary for printing all the possible microwave fields
!> @copyright
!! Copyright (C) 2008-2018 UppASD group
!! This file is distributed under the terms of the
!! GNU General Public License.
!! See http://www.gnu.org/copyleft/gpl.txt
module prn_microwaves
 
   use Parameters
   use Profiling
   use MicroWaveField

   implicit none

   real(dblprec), dimension(:), allocatable :: indxb_mwf                !< Step counter for the monochromatic microwave field
   real(dblprec), dimension(:), allocatable :: indxb_gauss              !< Step counter for the static gaussian shaped field
   real(dblprec), dimension(:), allocatable :: indxb_mwf_gauss          !< Step counter for the frequency broadened microwave field
   real(dblprec), dimension(:), allocatable :: indxb_mov_gauss          !< Step counter for the moving static gaussian shaped field
   real(dblprec), dimension(:), allocatable :: indxb_mov_circle         !< Step counter for the moving static circular shaped field
   real(dblprec), dimension(:), allocatable :: indxb_mov_square         !< Step counter for the moving static cubic shaped field
   real(dblprec), dimension(:), allocatable :: indxb_mwf_mov_gauss      !< Step counter for the moving gaussian shaped microwave field
   real(dblprec), dimension(:), allocatable :: indxb_mwf_mov_circle     !< Step counter for the moving circular shaped microwave field
   real(dblprec), dimension(:), allocatable :: indxb_mwf_mov_square     !< Step counter for the moving cubic shaped microwave field
   real(dblprec), dimension(:), allocatable :: indxb_mwf_gauss_spatial  !< Step counter for the gaussian shaped frequency broadened microwave field
   real(dblprec), dimension(:,:,:), allocatable :: mwfb               !< Buffer for the monochromatic microwave field
   real(dblprec), dimension(:,:,:), allocatable :: gaussb             !< Buffer the static gaussian shaped field
   real(dblprec), dimension(:,:,:), allocatable :: mwf_gaussb         !< Buffer for the frequency broadened microwave field
   real(dblprec), dimension(:,:,:), allocatable :: mov_gaussb         !< Buffer the moving static gaussian shaped field
   real(dblprec), dimension(:,:,:), allocatable :: mov_circleb        !< Buffer the moving static circular shaped field
   real(dblprec), dimension(:,:,:), allocatable :: mov_squareb        !< Buffer the moving static cubic shaped field
   real(dblprec), dimension(:,:,:), allocatable :: mwf_mov_gaussb     !< Buffer the moving gaussian shaped microwave field
   real(dblprec), dimension(:,:,:), allocatable :: mwf_mov_circleb    !< Buffer the moving circular shaped microwave field
   real(dblprec), dimension(:,:,:), allocatable :: mwf_mov_squareb    !< Buffer the moving cubic shaped microwave field
   real(dblprec), dimension(:,:,:), allocatable :: mwf_gauss_spatialb !< Buffer the gaussian shaped frequency broadened microwave field

   ! Counters for the measurements
   integer :: bcount_mwf               !< Counter of buffer for the monochromatic microwave field
   integer :: bcount_gauss             !< Counter of buffer for static gaussian shaped field
   integer :: bcount_mwf_gauss         !< Counter of buffer for the frequency broadened microwave field
   integer :: bcount_mov_gauss         !< Counter of buffer for the moving static gaussian shaped field
   integer :: bcount_mov_circle        !< Counter of buffer for the moving static circular shaped field
   integer :: bcount_mov_square        !< Counter of buffer for the moving static cubic shaped field
   integer :: bcount_mwf_mov_gauss     !< Counter of buffer for the moving gaussian shaped field microwave field
   integer :: bcount_mwf_mov_circle    !< Counter of buffer for the moving circular shaped field microwave field
   integer :: bcount_mwf_mov_square    !< Counter of buffer for the moving cubic shaped field microwave field
   integer :: bcount_mwf_gauss_spatial !< Counter of buffer for the gaussian shaped frequency broadened microwave field


   public
   !private

   !.. Private integers
   private :: bcount_mwf,bcount_gauss,bcount_mwf_gauss,bcount_mov_gauss,bcount_mwf_mov_gauss,&
      bcount_mwf_gauss_spatial,bcount_mov_circle,bcount_mwf_mov_circle,bcount_mov_square,bcount_mwf_mov_square
   !.. Private working arrays
   private :: indxb_mwf, indxb_gauss, indxb_mwf_gauss, indxb_mov_gauss,indxb_mov_circle,indxb_mwf_mov_gauss,&
      indxb_mwf_mov_circle,indxb_mov_square,indxb_mwf_mov_square,indxb_mwf_gauss_spatial
   private :: mwfb,gaussb,mwf_gaussb,mov_gaussb,mov_circleb,mwf_mov_gaussb,mwf_mov_circleb,mov_squareb,mwf_mov_squareb,mwf_gauss_spatialb

contains

   !> Printing wrapper for the microwave fields
   subroutine print_mwf_fields(sstep,mstep,Natom,delta_t,real_time_measure,simid)

      implicit none

      integer, intent(in) :: sstep ! Step in Logarithmic scale
      integer, intent(in) :: mstep !< Current simulation step
      integer, intent(in) :: Natom !< Number ofr atoms in the system

      real(dblprec), intent(in) :: delta_t              !< Current time step
      character(len=1), intent(in) :: real_time_measure !< Measurements displayed in real time
      character(len=8), intent(in) :: simid             !< Simulation name

      ! Monochromatic microwave field
      if (prn_mwf=='Y') then
         if (mod(sstep-1,mwf_step)==0) then
            call buffer_mwf(Natom, mstep-1, mwf,bcount_mwf,delta_t,real_time_measure,mwffield,site_mwffield)
            if (bcount_mwf==mwf_buff) then
               call print_mwf(Natom, mwf,simid,real_time_measure)
               bcount_mwf=1
            else
               bcount_mwf=bcount_mwf+1
            endif
         endif
      endif
      ! Static gaussian shaped field
      if (prn_gauss=='Y') then
         if (mod(sstep-1,gauss_step)==0) then
            call buffer_gauss(Natom, mstep-1, bcount_gauss,delta_t,real_time_measure,gauss_spatial,do_gauss)
            if (bcount_gauss==gauss_buff) then
               call print_gauss(Natom,simid,real_time_measure)
               bcount_gauss=1
            else
               bcount_gauss=bcount_gauss+1
            endif
         endif
      endif
      ! Frequency broadened microwave field
      if (prn_mwf_gauss=='Y') then
         if (mod(sstep-1,mwf_gauss_step)==0) then
            call buffer_mwf_gauss(Natom, mstep-1, mwf_gauss,bcount_mwf_gauss,delta_t,&
               real_time_measure,gauss_mwffield,gauss_site_mwffield)
            if (bcount_mwf_gauss==mwf_gauss_buff) then
               call print_mwf_gauss(Natom, mwf_gauss,simid,real_time_measure)
               bcount_mwf_gauss=1
            else
               bcount_mwf_gauss=bcount_mwf_gauss+1
            endif
         endif
      endif
      ! Moving static gaussian shaped field
      if (prn_mov_gauss=='Y') then
         if (mod(sstep-1,mov_gauss_pstep)==0) then
            call buffer_mov_gauss(Natom, mstep-1, bcount_mov_gauss,delta_t,&
               real_time_measure,mov_gauss_spatial,mov_gauss)
            if (bcount_mov_gauss==mov_gauss_buff) then
               call print_mov_gauss(Natom,simid,real_time_measure)
               bcount_mov_gauss=1
            else
               bcount_mov_gauss=bcount_mov_gauss+1
            endif
         endif
      endif
      ! Moving gaussian shaped microwave field
      if (prn_mwf_mov_gauss=='Y') then
         if (mod(sstep-1,mwf_mov_gauss_pstep)==0) then
            call buffer_mwf_mov_gauss(Natom, mstep-1, bcount_mwf_mov_gauss,delta_t,&
               real_time_measure,mwf_mov_gauss_spatial,mwf_mov_gauss)
            if (bcount_mwf_mov_gauss==mwf_mov_gauss_buff) then
               call print_mwf_mov_gauss(Natom,simid,real_time_measure)
               bcount_mwf_mov_gauss=1
            else
               bcount_mwf_mov_gauss=bcount_mwf_mov_gauss+1
            endif
         endif
      endif
      ! Moving static circular shaped field
      if (prn_mov_circle=='Y') then
         if (mod(sstep-1,mov_circle_pstep)==0) then
            call buffer_mov_circle(Natom, mstep-1, bcount_mov_circle,delta_t,&
               real_time_measure,mov_circle_spatial,mov_circle)
            if (bcount_mov_circle==mov_circle_buff) then
               call print_mov_circle(Natom,simid,real_time_measure)
               bcount_mov_circle=1
            else
               bcount_mov_circle=bcount_mov_circle+1
            endif
         endif
      endif
      ! Moving circular shaped microwave field
      if (prn_mwf_mov_circle=='Y') then
         if (mod(sstep-1,mwf_mov_circle_pstep)==0) then
            call buffer_mwf_mov_circle(Natom, mstep-1, bcount_mwf_mov_circle,delta_t,&
               real_time_measure,mwf_mov_circle_spatial,mwf_mov_circle)
            if (bcount_mwf_mov_circle==mwf_mov_circle_buff) then
               call print_mwf_mov_circle(Natom,simid,real_time_measure)
               bcount_mwf_mov_circle=1
            else
               bcount_mwf_mov_circle=bcount_mwf_mov_circle+1
            endif
         endif
      endif

      ! Moving static cubic shaped field
      if (prn_mov_square=='Y') then
         if (mod(sstep-1,mov_square_pstep)==0) then
            call buffer_mov_square(Natom, mstep-1, bcount_mov_square,delta_t,&
               real_time_measure,mov_square_spatial,mov_square)
            if (bcount_mov_square==mov_square_buff) then
               call print_mov_square(Natom,simid,real_time_measure)
               bcount_mov_square=1
            else
               bcount_mov_square=bcount_mov_square+1
            endif
         endif
      endif
      ! Moving cubic shaped microwave field
      if (prn_mwf_mov_square=='Y') then
         if (mod(sstep-1,mwf_mov_square_pstep)==0) then
            call buffer_mwf_mov_square(Natom, mstep-1, bcount_mwf_mov_square,delta_t,&
               real_time_measure,mwf_mov_square_spatial,mwf_mov_square)
            if (bcount_mwf_mov_square==mwf_mov_square_buff) then
               call print_mwf_mov_square(Natom,simid,real_time_measure)
               bcount_mwf_mov_square=1
            else
               bcount_mwf_mov_square=bcount_mwf_mov_square+1
            endif
         endif
      endif

      ! Gaussian shaped frequency broadened microwave field
      if (prn_mwf_gauss_spatial=='Y') then
         if (mod(sstep-1,mwf_gauss_spatial_step)==0) then
            call buffer_mwf_gauss_spatial(Natom, mstep-1, mwf_gauss_spatial,bcount_mwf_gauss_spatial,delta_t,&
               real_time_measure,gauss_spatial_site_mwffield)
            if (bcount_mwf_gauss_spatial==mwf_gauss_spatial_buff) then
               call print_mwf_gauss_spatial(Natom,simid,real_time_measure)
               bcount_mwf_gauss_spatial=1
            else
               bcount_mwf_gauss_spatial=bcount_mwf_gauss_spatial+1
            endif
         endif
      endif

   end subroutine print_mwf_fields

   !> Flush the microwave fields
   subroutine flush_mwf_fields(Natom,real_time_measure,simid)

      implicit none

      integer, intent(in) :: Natom !< Number ofr atoms in the system
      character(len=1), intent(in) :: real_time_measure !< Measurements displayed in real time
      character(len=8), intent(in) :: simid             !< Simulation name

      ! Flush the monochromatic microwave field
      if (prn_mwf=='Y') then
         bcount_mwf=bcount_mwf-1
         call print_mwf(Natom, mwf,simid,real_time_measure)
      endif
      ! Flush the static gaussian shaped field
      if (prn_gauss=='Y') then
         bcount_gauss=bcount_gauss-1
         call print_gauss(Natom,simid,real_time_measure)
      endif
      ! Flush the frequency broadened microwave field
      if (prn_mwf_gauss=='Y') then
         bcount_mwf_gauss=bcount_mwf_gauss-1
         call print_mwf_gauss(Natom, mwf_gauss,simid,real_time_measure)
      endif
      ! Flush the moving gaussian shaped static field
      if (prn_mov_gauss=='Y') then
         bcount_mov_gauss=bcount_mov_gauss-1
         call print_mov_gauss(Natom,simid,real_time_measure)
      endif
      ! Flush the moving gaussian shaped microwave field
      if (prn_mwf_mov_gauss=='Y') then
         bcount_mwf_mov_gauss=bcount_mwf_mov_gauss-1
         call print_mwf_mov_gauss(Natom,simid,real_time_measure)
      endif
      ! Flush the moving circular shaped static field
      if (prn_mov_circle=='Y') then
         bcount_mov_circle=bcount_mov_circle-1
         call print_mov_circle(Natom,simid,real_time_measure)
      endif
      ! Flush the moving circular shaped microwave field
      if (prn_mwf_mov_circle=='Y') then
         bcount_mwf_mov_circle=bcount_mwf_mov_circle-1
         call print_mwf_mov_circle(Natom,simid,real_time_measure)
      endif
      ! Flush the moving cubic shaped static field
      if (prn_mov_square=='Y') then
         bcount_mov_square=bcount_mov_square-1
         call print_mov_square(Natom,simid,real_time_measure)
      endif
      ! Flush the moving cubic shaped microwave field
      if (prn_mwf_mov_square=='Y') then
         bcount_mwf_mov_square=bcount_mwf_mov_square-1
         call print_mwf_mov_square(Natom,simid,real_time_measure)
      endif

      ! Flush the gaussian shaped frequency broadened microwave field
      if (prn_mwf_gauss_spatial=='Y') then
         bcount_mwf_gauss_spatial=bcount_mwf_gauss_spatial-1
         call print_mwf_gauss_spatial(Natom,simid,real_time_measure)
      endif

   end subroutine flush_mwf_fields

   !> Initialize the arrays and scalars for printing the the microwave fields
   subroutine mwf_initializations(Natom,flag)

      implicit none

      integer, intent(in) :: flag  !< Allocate or deallocate (1/-1)
      integer, intent(in) :: Natom !< Total number of atoms in the system

      !.. Local variables
      integer :: i_stat, i_all

      if (flag>0) then
         !.. Initiate microwave field measurement buffer counters
         bcount_mwf               = 1
         bcount_gauss             = 1
         bcount_mwf_gauss         = 1
         bcount_mov_gauss         = 1
         bcount_mov_circle        = 1
         bcount_mov_square        = 1
         bcount_mwf_mov_gauss     = 1
         bcount_mwf_mov_circle    = 1
         bcount_mwf_mov_square    = 1
         bcount_mwf_gauss_spatial = 1

         ! Printing the monochromatic microwave field
         if (prn_mwf=='Y') then
            allocate(indxb_mwf(mwf_buff),stat=i_stat)
            call memocc(i_stat,product(shape(indxb_mwf))*kind(indxb_mwf),'indxb_mwf','mwf_initializations')
            allocate(mwfb(3,Natom,mwf_buff),stat=i_stat)
            call memocc(i_stat,product(shape(mwfb))*kind(mwfb),'mwfb','mwf_initializations')
         endif

         ! Printing the static gaussian shaped field
         if (prn_gauss=='Y') then
            allocate(indxb_gauss(gauss_buff),stat=i_stat)
            call memocc(i_stat,product(shape(indxb_gauss))*kind(indxb_gauss),'indxb_gauss','mwf_initializations')
            allocate(gaussb(3,Natom,gauss_buff),stat=i_stat)
            call memocc(i_stat,product(shape(gaussb))*kind(gaussb),'gaussb','mwf_initializations')
         endif

         ! Printing the frequency broadened microwave field
         if (prn_mwf_gauss=='Y') then
            allocate(indxb_mwf_gauss(mwf_gauss_buff),stat=i_stat)
            call memocc(i_stat,product(shape(indxb_mwf_gauss))*kind(indxb_mwf_gauss),'indxb_mwf_gauss','mwf_initializations')
            allocate(mwf_gaussb(3,Natom,mwf_gauss_buff),stat=i_stat)
            call memocc(i_stat,product(shape(mwf_gaussb))*kind(mwf_gaussb),'mwf_gaussb','mwf_initializations')
         endif

         ! Printing the moving static gaussian shaped field
         if (prn_mov_gauss=='Y') then
            allocate(indxb_mov_gauss(mov_gauss_buff),stat=i_stat)
            call memocc(i_stat,product(shape(indxb_mov_gauss))*kind(indxb_mov_gauss),'indxb_mov_gauss','mwf_initializations')
            allocate(mov_gaussb(3,Natom,mov_gauss_buff),stat=i_stat)
            call memocc(i_stat,product(shape(mov_gaussb))*kind(mov_gaussb),'mov_gaussb','mwf_initializations')
         endif

         ! Printing the moving gaussian shaped microwave field
         if (prn_mwf_mov_gauss=='Y') then
            allocate(indxb_mwf_mov_gauss(mwf_mov_gauss_buff),stat=i_stat)
            call memocc(i_stat,product(shape(indxb_mwf_mov_gauss))*kind(indxb_mwf_mov_gauss),'indxb_mwf_mov_gauss','mwf_initializations')
            allocate(mov_gaussb(3,Natom,mov_gauss_buff),stat=i_stat)
            call memocc(i_stat,product(shape(mov_gaussb))*kind(mov_gaussb),'mov_gaussb','mwf_initializations')
         endif

         ! Printing the moving static circular shaped field
         if (prn_mov_circle=='Y') then
            allocate(indxb_mov_circle(mov_circle_buff),stat=i_stat)
            call memocc(i_stat,product(shape(indxb_mov_circle))*kind(indxb_mov_circle),'indxb_mov_circle','mwf_initializations')
            allocate(mov_circleb(3,Natom,mov_circle_buff),stat=i_stat)
            call memocc(i_stat,product(shape(mov_circleb))*kind(mov_circleb),'mov_circleb','mwf_initializations')
         endif

         ! Printing the moving circular shaped microwave field
         if (prn_mwf_mov_circle=='Y') then
            allocate(indxb_mwf_mov_circle(mwf_mov_circle_buff),stat=i_stat)
            call memocc(i_stat,product(shape(indxb_mwf_mov_circle))*kind(indxb_mwf_mov_circle),'indxb_mwf_mov_circle','mwf_initializations')
            allocate(mov_circleb(3,Natom,mov_circle_buff),stat=i_stat)
            call memocc(i_stat,product(shape(mov_circleb))*kind(mov_circleb),'mov_circleb','mwf_initializations')
         endif

         ! Printing the moving static cubic shaped field
         if (prn_mov_square=='Y') then
            allocate(indxb_mov_square(mov_square_buff),stat=i_stat)
            call memocc(i_stat,product(shape(indxb_mov_square))*kind(indxb_mov_square),'indxb_mov_square','mwf_initializations')
            allocate(mov_squareb(3,Natom,mov_square_buff),stat=i_stat)
            call memocc(i_stat,product(shape(mov_squareb))*kind(mov_squareb),'mov_squareb','mwf_initializations')
         endif

         ! Printing the moving cubic shaped microwave field
         if (prn_mwf_mov_square=='Y') then
            allocate(indxb_mwf_mov_square(mwf_mov_square_buff),stat=i_stat)
            call memocc(i_stat,product(shape(indxb_mwf_mov_square))*kind(indxb_mwf_mov_square),'indxb_mwf_mov_square','mwf_initializations')
            allocate(mov_squareb(3,Natom,mov_square_buff),stat=i_stat)
            call memocc(i_stat,product(shape(mov_squareb))*kind(mov_squareb),'mov_squareb','mwf_initializations')
         endif

         ! Printing the gaussian shaped frequency broadened microwave field
         if (prn_mwf_gauss_spatial=='Y') then
            allocate(indxb_mwf_gauss_spatial(mwf_gauss_spatial_buff),stat=i_stat)
            call memocc(i_stat,product(shape(indxb_mwf_gauss_spatial))*kind(indxb_mwf_gauss_spatial),'indxb_mwf_gauss_spatial','mwf_initializations')
            allocate(mwf_gauss_spatialb(3,Natom,mwf_gauss_spatial_buff),stat=i_stat)
            call memocc(i_stat,product(shape(mwf_gauss_spatialb))*kind(mwf_gauss_spatialb),'mwf_gauss_spatialb','mwf_initializations')
         endif

      else

         if (prn_mwf=='Y') then
            i_all=-product(shape(indxb_mwf))*kind(indxb_mwf)
            deallocate(indxb_mwf,stat=i_stat)
            call memocc(i_stat,i_all,'indxb_mwf','mwf_initializations')
            i_all=-product(shape(mwfb))*kind(mwfb)
            deallocate(mwfb,stat=i_stat)
            call memocc(i_stat,i_all,'mwfb','mwf_initializations')
         endif

         if (prn_gauss=='Y') then
            i_all=-product(shape(indxb_gauss))*kind(indxb_gauss)
            deallocate(indxb_gauss,stat=i_stat)
            call memocc(i_stat,i_all,'indxb_gauss','mwf_initializations')
            i_all=-product(shape(gaussb))*kind(gaussb)
            deallocate(gaussb,stat=i_stat)
            call memocc(i_stat,i_all,'gaussb','mwf_initializations')
         endif

         if (prn_mwf_gauss=='Y') then
            i_all=-product(shape(indxb_mwf_gauss))*kind(indxb_mwf_gauss)
            deallocate(indxb_mwf_gauss,stat=i_stat)
            call memocc(i_stat,i_all,'indxb_mwf_gauss','mwf_initializations')
            i_all=-product(shape(mwf_gaussb))*kind(mwf_gaussb)
            deallocate(mwf_gaussb,stat=i_stat)
            call memocc(i_stat,i_all,'mwf_gaussb','mwf_initializations')
         endif

         if (prn_mov_gauss=='Y') then
            i_all=-product(shape(indxb_mov_gauss))*kind(indxb_mov_gauss)
            deallocate(indxb_mov_gauss,stat=i_stat)
            call memocc(i_stat,i_all,'indxb_mov_gauss','mwf_initializations')
            i_all=-product(shape(mov_gaussb))*kind(mov_gaussb)
            deallocate(mov_gaussb,stat=i_stat)
            call memocc(i_stat,i_all,'mov_gaussb','mwf_initializations')
         endif

         if (prn_mwf_mov_gauss=='Y') then
            i_all=-product(shape(indxb_mwf_mov_gauss))*kind(indxb_mwf_mov_gauss)
            deallocate(indxb_mwf_mov_gauss,stat=i_stat)
            call memocc(i_stat,i_all,'indxb_mwf_mov_gauss','mwf_initializations')
            i_all=-product(shape(mwf_mov_gaussb))*kind(mwf_mov_gaussb)
            deallocate(mwf_mov_gaussb,stat=i_stat)
            call memocc(i_stat,i_all,'mwf_mov_gaussb','mwf_initializations')
         endif

         if (prn_mov_circle=='Y') then
            i_all=-product(shape(indxb_mov_circle))*kind(indxb_mov_circle)
            deallocate(indxb_mov_circle,stat=i_stat)
            call memocc(i_stat,i_all,'indxb_mov_circle','mwf_initializations')
            i_all=-product(shape(mov_circleb))*kind(mov_circleb)
            deallocate(mov_circleb,stat=i_stat)
            call memocc(i_stat,i_all,'mov_circleb','mwf_initializations')
         endif

         if (prn_mwf_mov_circle=='Y') then
            i_all=-product(shape(indxb_mwf_mov_circle))*kind(indxb_mwf_mov_circle)
            deallocate(indxb_mwf_mov_circle,stat=i_stat)
            call memocc(i_stat,i_all,'indxb_mwf_mov_circle','mwf_initializations')
            i_all=-product(shape(mwf_mov_circleb))*kind(mwf_mov_circleb)
            deallocate(mwf_mov_circleb,stat=i_stat)
            call memocc(i_stat,i_all,'mwf_mov_circleb','mwf_initializations')
         endif

         if (prn_mov_square=='Y') then
            i_all=-product(shape(indxb_mov_square))*kind(indxb_mov_square)
            deallocate(indxb_mov_square,stat=i_stat)
            call memocc(i_stat,i_all,'indxb_mov_square','mwf_initializations')
            i_all=-product(shape(mov_squareb))*kind(mov_squareb)
            deallocate(mov_squareb,stat=i_stat)
            call memocc(i_stat,i_all,'mov_squareb','mwf_initializations')
         endif

         if (prn_mwf_mov_square=='Y') then
            i_all=-product(shape(indxb_mwf_mov_square))*kind(indxb_mwf_mov_square)
            deallocate(indxb_mwf_mov_square,stat=i_stat)
            call memocc(i_stat,i_all,'indxb_mwf_mov_square','mwf_initializations')
            i_all=-product(shape(mwf_mov_squareb))*kind(mwf_mov_squareb)
            deallocate(mwf_mov_squareb,stat=i_stat)
            call memocc(i_stat,i_all,'mwf_mov_squareb','mwf_initializations')
         endif

         if (prn_mwf_gauss_spatial=='Y') then
            i_all=-product(shape(indxb_mwf_gauss_spatial))*kind(indxb_mwf_gauss_spatial)
            deallocate(indxb_mwf_gauss_spatial,stat=i_stat)
            call memocc(i_stat,i_all,'indxb_mwf_gauss_spatial','mwf_initializations')
            i_all=-product(shape(mwf_gauss_spatialb))*kind(mwf_gauss_spatialb)
            deallocate(mwf_gauss_spatialb,stat=i_stat)
            call memocc(i_stat,i_all,'mwf_gauss_spatialb','mwf_initializations')
         endif

      endif

   end subroutine mwf_initializations

   !> Subroutine to initialize the necesary parameters for the microwave field printing
   subroutine initialize_prn_microwaves()

      implicit none

      real(dblprec) :: zero=0.0_dblprec

      ! Monochromatic microwave field
      prn_mwf  = 'N'
      mwf_step = 1000
      mwf_buff = 100

      ! Static gaussian shaped field
      prn_gauss  = 'N'
      gauss_step = 1000
      gauss_buff = 100

      ! Frequency broadened microwave field
      prn_mwf_gauss  = 'N'
      mwf_gauss_step = 1000
      mwf_gauss_buff = 100

      ! Moving static gaussian shaped field
      prn_mov_gauss   = 'N'
      mov_gauss_buff  = 100
      mov_gauss_pstep = 1000

      ! Moving gaussian shaped microwave field
      prn_mwf_mov_gauss   = 'N'
      mwf_mov_gauss_buff  = 100
      mwf_mov_gauss_pstep = 1000

      ! Moving static circular shaped field
      prn_mov_circle   = 'N'
      mov_circle_buff  = 100
      mov_circle_pstep = 1000

      ! Moving circular shaped microwave field
      prn_mwf_mov_circle   = 'N'
      mwf_mov_circle_buff  = 100
      mwf_mov_circle_pstep = 1000

      ! Gaussian shaped frequency broadened microwave field
      prn_mwf_gauss_spatial  = 'N'
      mwf_gauss_spatial_step = 1000
      mwf_gauss_spatial_buff = 100

      !Monochromatic Microwave field
      mwf             = 'N'
      site_phase      = 'N'
      mwf_site_file   = 'mwf_site_file'
      site_phase_file = 'phase_site_file'
      mwf_pulse_time  = 0
      mwfampl         = 0.0D0
      mwffreq         = 0.0D0
      mwfdir          = (/zero,zero,zero/)

      !Gaussian frequency broadened Microwave fields
      mwf_gauss            = 'N'
      mwf_gauss_site_file  = 'mwf_gauss_site_file'
      mwf_gauss_pulse_time = 0
      mwf_gauss_ampl       = 0.0D0
      mwf_gauss_freq       = 0.0D0
      mwf_gauss_time_sigma = 0.0D0
      mwf_gauss_dir        = (/zero,zero,zero/)

      !Spatial Gaussian shaped frequency broadened microwave fields
      mwf_gauss_spatial             = 'N'
      mwf_gauss_spatial_site_file   = 'mwf_gauss_spatial_site_file'
      mwf_gauss_spatial_pulse_time  = 0
      mwf_gauss_spatial_freq        = 0.0D0
      mwf_gauss_spatial_ampl        = 0.0D0
      mwf_gauss_spatial_time_sigma  = 0.0D0
      mwf_gauss_spatial_space_sigma = (/zero,zero,zero/)

      ! Static Gaussian shaped microwave field
      do_gauss            = 'N'
      gauss_site_file     = 'gauss_site_file'
      gauss_pulse_time    = 0
      gauss_spatial_ampl  = 0.0D0
      gauss_spatial_sigma = (/zero,zero,zero/)

      !Moving gaussian shaped static field
      mov_gauss             = 'N'
      mov_gauss_file        = 'mov_gauss_file'
      mov_gauss_step        = 0
      mov_gauss_pulse_time  = 0
      mov_gauss_ampl        = 0.0D0
      mov_gauss_space_sigma = (/zero,zero,zero/)

      !Moving gaussian shaped microwave field
      mwf_mov_gauss             = 'N'
      mwf_mov_gauss_file        = 'mwf_mov_gauss_file'
      mwf_mov_gauss_step        = 0
      mwf_mov_gauss_pulse_time  = 0
      mwf_mov_gauss_ampl        = 0.0D0
      mwf_mov_gauss_freq        = 0.0D0
      mwf_mov_gauss_time_sigma  = 0.0D0
      mwf_mov_gauss_space_sigma = (/zero,zero,zero/)

      !Moving circular shaped static field
      mov_circle             = 'N'
      mov_circle_file        = 'mov_circle_file'
      mov_circle_step        = 0
      mov_circle_pulse_time  = 0
      mov_circle_ampl        = 0.0D0
      mov_circle_radius      = 0.0D0

      !Moving circular shaped microwave field
      mwf_mov_circle             = 'N'
      mwf_mov_circle_file        = 'mwf_mov_circle_file'
      mwf_mov_circle_step        = 0
      mwf_mov_circle_pulse_time  = 0
      mwf_mov_circle_ampl        = 0.0D0
      mwf_mov_circle_freq        = 0.0D0
      mwf_mov_circle_time_sigma  = 0.0D0
      mwf_mov_circle_radius      = 0.0D0

      !Moving cubic shaped static field
      mov_square             = 'N'
      mov_square_file        = 'mov_square_file'
      mov_square_step        = 0
      mov_square_pulse_time  = 0
      mov_square_ampl        = 0.0D0
      mov_square_dimensions  = (/zero,zero,zero/)

      !Moving cubic shaped microwave field
      mwf_mov_square             = 'N'
      mwf_mov_square_file        = 'mwf_mov_square_file'
      mwf_mov_square_step        = 0
      mwf_mov_square_pulse_time  = 0
      mwf_mov_square_ampl        = 0.0D0
      mwf_mov_square_freq        = 0.0D0
      mwf_mov_square_time_sigma  = 0.0D0
      mwf_mov_square_dimensions  = (/zero,zero,zero/)

   end subroutine initialize_prn_microwaves

   !> Buffer the monochromatic microwave field
   subroutine buffer_mwf(Natom, mstep, mwf,bcount_mwf,delta_t,real_time_measure,mwffield,site_mwffield)

      implicit none

      integer, intent(in) :: mstep          !< Current simulation step
      integer, intent(in) :: Natom          !< Number of atoms in system
      integer, intent(in) :: bcount_mwf     !< Counter of buffer for monochromatic microwave field
      real(dblprec), intent(in) :: delta_t  !< Current time step (used for real time measurements)
      real(dblprec), dimension(3), intent(in) :: mwffield            !< Global monochromatic microwave field
      real(dblprec), dimension(3,Natom), intent(in) :: site_mwffield !< Site dependent monochromatic microwave field
      character(len=1), intent(in) :: mwf               !< Add monochromatic microwave field (Y/P/S/W/N)
      character(len=1), intent(in) :: real_time_measure !< Real time measurement flag

      ! .. Local variables
      integer :: i

      if (mwf=='Y'.or.mwf=='P') then
         do i=1, Natom
            ! If the field is not site dependent mwffield is the one stored in the buffer
            mwfb(1:3,i,bcount_mwf) = mwffield(1:3)
         enddo
      else if (mwf=='S'.or.mwf=='W') then
         ! If the field is site dependent one saves the field by sites
         do i=1, Natom
            mwfb(1:3,i,bcount_mwf)=site_mwffield(1:3,i)
         end do
      endif

      if (real_time_measure=='Y') then
         indxb_mwf(bcount_mwf)=mstep*delta_t
      else
         indxb_mwf(bcount_mwf)=mstep
      endif

   end subroutine buffer_mwf

   !> Print the monochromatic microwave field
   subroutine print_mwf(Natom, mwf,simid,real_time_measure)
      !
      !.. Implicit declarations
      implicit none

      integer, intent(in) :: Natom        !< Number of atoms in system
      character(len=1), intent(in) :: mwf !< Add monochromatic microwave field
      character(len=1), intent(in) :: real_time_measure !< Measurements displayed in real time
      character(len=8), intent(in) :: simid !< Name of simulation

      integer :: i, j
      character(len=30) :: filn

      ! This is the microwave field which acts uniformly over the sample
      if(mwf=='Y'.or.mwf=='P') then

         write (filn,'(''mwf.'',a8,''.out'')') simid
         open(ofileno, file=filn, position="append")
         do i=1, bcount_mwf
            if (real_time_measure=='Y') then
               write (ofileno,10005) indxb_mwf(i), mwfb(1,1,i), mwfb(2,1,i), mwfb(3,1,i)
            else
               write (ofileno,10004) int(indxb_mwf(i)), mwfb(1,1,i), mwfb(2,1,i), mwfb(3,1,i)
            endif
         end do
         close(ofileno)

         ! This is a microwave field which acts over selected atoms
      else if (mwf=='S'.or.mwf=='P') then

         write (filn,'(''mwf.'',a8,''.out'')') simid
         open(ofileno, file=filn, position="append")
         do i=1, bcount_mwf
            do j=1, Natom
               if (real_time_measure=='Y') then
                  write (ofileno,10007) indxb_mwf(i),j, mwfb(1,j,i), mwfb(2,j,i), mwfb(3,j,i)
               else
                  write (ofileno,10006) int(indxb_mwf(i)),j, mwfb(1,j,i), mwfb(2,j,i), mwfb(3,j,i)
               endif
            enddo
         end do
         close(ofileno)
      endif

      return

      write (*,*) "Error writing the monochromatic microwave field file"

      10004 format (i8,3es16.8)
      10005 format (es16.8,3es16.8)
      10006 format (i8,i8,3es18.9E3)
      10007 format (es16.8,i8,3es18.9E3)

   end subroutine print_mwf

   !> Buffer the gaussian broadened microwave field
   subroutine buffer_mwf_gauss(Natom, mstep, mwf_gauss,bcount_mwf_gauss,delta_t,&
         real_time_measure,gauss_mwffield,gauss_site_mwffield)

      implicit none

      integer, intent(in) :: mstep            !< Current simulation step
      integer, intent(in) :: Natom            !< Number of atoms in system
      integer, intent(in) :: bcount_mwf_gauss !< Counter of buffer for monochromatic microwave field
      real(dblprec), intent(in) :: delta_t    !< Current time step (used for real time measurements)
      real(dblprec), dimension(3), intent(in) :: gauss_mwffield            !< Global frequency broadenend microwave field
      real(dblprec), dimension(3,Natom), intent(in) :: gauss_site_mwffield !< Site dependent frequency broadenend microwave field
      character(len=1), intent(in) :: mwf_gauss         !< Add frequency broadened microwave field
      character(len=1), intent(in) :: real_time_measure !< Real time measurement flag

      ! .. Local variables
      integer :: i

      if (mwf_gauss=='Y'.or.mwf_gauss=='P') then
         do i=1, Natom
            ! If the field is not site dependent gauss_mwffield is the one stored in the buffer
            mwf_gaussb(1:3,i,bcount_mwf_gauss) = gauss_mwffield(1:3)
         enddo
      else if (mwf_gauss=='S'.or.mwf_gauss=='W') then
         ! If the field is site dependent one saves the field by sites
         do i=1, Natom
            mwf_gaussb(1:3,i,bcount_mwf_gauss)=gauss_site_mwffield(1:3,i)
         end do
      endif

      if (real_time_measure=='Y') then
         indxb_mwf_gauss(bcount_mwf_gauss)=mstep*delta_t
      else
         indxb_mwf_gauss(bcount_mwf_gauss)=mstep
      endif

   end subroutine buffer_mwf_gauss

   !> Print the frequency broadened microwave field
   subroutine print_mwf_gauss(Natom, mwf_gauss,simid,real_time_measure)
      !
      !.. Implicit declarations
      implicit none

      integer, intent(in) :: Natom     !< Number of atoms in system
      character(len=1), intent(in) :: mwf_gauss         !< Add frequency broadened microwave field
      character(len=1), intent(in) :: real_time_measure !< Measurements displayed in real time
      character(len=8), intent(in) :: simid !< Name of simulation

      integer :: i, j
      character(len=30) :: filn

      ! This is the microwave field which acts uniformly over the sample
      if(mwf_gauss=='Y'.or.mwf_gauss=='P') then

         write (filn,'(''mwf_gauss.'',a8,''.out'')') simid
         open(ofileno, file=filn, position="append")
         do i=1, bcount_mwf_gauss
            if (real_time_measure=='Y') then
               write (ofileno,10005) indxb_mwf_gauss(i), mwf_gaussb(1,1,i), mwf_gaussb(2,1,i), mwf_gaussb(3,1,i)
            else
               write (ofileno,10004) int(indxb_mwf_gauss(i)), mwf_gaussb(1,1,i), mwf_gaussb(2,1,i), mwf_gaussb(3,1,i)
            endif
         end do
         close(ofileno)

         ! This is a microwave field which acts over selected atoms
      else if (mwf_gauss=='S'.or.mwf_gauss=='P') then

         write (filn,'(''mwf_gauss.'',a8,''.out'')') simid
         open(ofileno, file=filn, position="append")
         do i=1, bcount_mwf_gauss
            do j=1, Natom
               if (real_time_measure=='Y') then
                  write (ofileno,10007) indxb_mwf_gauss(i),j, mwf_gaussb(1,j,i), mwf_gaussb(2,j,i), mwf_gaussb(3,j,i)
               else
                  write (ofileno,10006) int(indxb_mwf_gauss(i)),j, mwf_gaussb(1,j,i), mwf_gaussb(2,j,i), mwf_gaussb(3,j,i)
               endif
            enddo
         end do
         close(ofileno)
      endif

      return

      write (*,*) "Error writing the gaussian broadened microwave field file"

      10004 format (i8,3es16.8)
      10005 format (es16.8,3es16.8)
      10006 format (i8,i8,3es18.9E3)
      10007 format (es16.8,i8,3es18.9E3)

   end subroutine print_mwf_gauss

   !> Buffer the frequency broadened gaussian shaped microwave field
   subroutine buffer_mwf_gauss_spatial(Natom, mstep, mwf_gauss_spatial,bcount_mwf_gauss_spatial,delta_t,&
         real_time_measure,gauss_spatial_site_mwffield)

      implicit none

      integer, intent(in) :: mstep            !< Current simulation step
      integer, intent(in) :: Natom            !< Number of atoms in system
      integer, intent(in) :: bcount_mwf_gauss_spatial !< Counter of buffer for monochromatic microwave field
      real(dblprec), intent(in) :: delta_t            !< Current time step (used for real time measurements)
      real(dblprec), dimension(3,Natom), intent(in) :: gauss_spatial_site_mwffield !< Array with the gaussian broadened spatially gaussian microwave field
      character(len=1), intent(in) :: mwf_gauss_spatial !< Add the frequency broadened gaussian shaped microwave field
      character(len=1), intent(in) :: real_time_measure !< Real time measurement flag

      ! .. Local variables
      integer :: i

      if (mwf_gauss_spatial=='Y'.or.mwf_gauss_spatial=='P') then
         do i=1, Natom
            mwf_gauss_spatialb(1:3,i,bcount_mwf_gauss_spatial)=gauss_spatial_site_mwffield(1:3,i)
         end do
      endif

      if (real_time_measure=='Y') then
         indxb_mwf_gauss_spatial(bcount_mwf_gauss_spatial)=mstep*delta_t
      else
         indxb_mwf_gauss_spatial(bcount_mwf_gauss_spatial)=mstep
      endif

   end subroutine buffer_mwf_gauss_spatial

   !> Print the frequency broadened gaussian shaped microwave field
   subroutine print_mwf_gauss_spatial(Natom,simid,real_time_measure)
      !
      !.. Implicit declarations
      implicit none

      integer, intent(in) :: Natom     !< Number of atoms in system
      character(len=8), intent(in) :: simid !< Name of simulation
      character(len=1), intent(in) :: real_time_measure !< Measurements displayed in real time

      !.. Local variables
      integer :: i, j
      character(len=30) :: filn

      ! Print to file the frequency broadened gaussian shaped microwave field
      write (filn,'(''mwf_gauss_s.'',a8,''.out'')') simid
      open(ofileno, file=filn, position="append")
      do i=1, bcount_mwf_gauss_spatial
         do j=1, Natom
            if (real_time_measure=='Y') then
               write (ofileno,10007) indxb_mwf_gauss_spatial(i),j, mwf_gauss_spatialb(1,j,i), mwf_gauss_spatialb(2,j,i), mwf_gauss_spatialb(3,j,i)
            else
               write (ofileno,10006) int(indxb_mwf_gauss_spatial(i)),j, mwf_gauss_spatialb(1,j,i), mwf_gauss_spatialb(2,j,i), mwf_gauss_spatialb(3,j,i)
            endif
         enddo
      end do
      close(ofileno)

      return

      write (*,*) "Error writing the gaussian shaped frequency broadened microwave field file"

      10006 format (i8,i8,3es18.9E3)
      10007 format (es16.8,i8,3es18.9E3)

   end subroutine print_mwf_gauss_spatial

   !> Buffer the static gaussian shaped field
   subroutine buffer_gauss(Natom, mstep, bcount_gauss,delta_t,&
         real_time_measure,gauss_spatial,do_gauss)

      implicit none

      integer, intent(in) :: mstep         !< Current simulation step
      integer, intent(in) :: Natom         !< Number of atoms in system
      integer, intent(in) :: bcount_gauss  !< Counter of buffer for the static gaussian shaped field
      real(dblprec), intent(in) :: delta_t !< Current time step (used for real time measurements)
      real(dblprec), dimension(3,Natom), intent(in) :: gauss_spatial !< Array with the static gaussian shaped field
      character(len=1), intent(in) :: do_gauss          !< Add the static gaussian shaped field
      character(len=1), intent(in) :: real_time_measure !< Real time measurement flag

      ! .. Local variables
      integer :: i

      if (do_gauss=='Y'.or.do_gauss=='P') then
         do i=1, Natom
            gaussb(1:3,i,bcount_gauss)=gauss_spatial(1:3,i)
         end do
      endif

      if (real_time_measure=='Y') then
         indxb_gauss(bcount_gauss)=mstep*delta_t
      else
         indxb_gauss(bcount_gauss)=mstep
      endif

   end subroutine buffer_gauss

   !> Print the static gaussian shaped field
   subroutine print_gauss(Natom,simid,real_time_measure)
      !
      !.. Implicit declarations
      implicit none

      integer, intent(in) :: Natom          !< Number of atoms in system
      character(len=8), intent(in) :: simid !< Name of simulation
      character(len=1), intent(in) :: real_time_measure !< Measurements displayed in real time

      !.. Local variables
      integer :: i, j
      character(len=30) :: filn

      ! Print ot file the static gaussian shaped field
      write (filn,'(''gauss.'',a8,''.out'')') simid
      open(ofileno, file=filn, position="append")
      do i=1, bcount_mwf_gauss_spatial
         do j=1, Natom
            if (real_time_measure=='Y') then
               write (ofileno,10007) indxb_gauss(i),j, gaussb(1,j,i), gaussb(2,j,i), gaussb(3,j,i)
            else
               write (ofileno,10006) int(indxb_gauss(i)),j, gaussb(1,j,i), gaussb(2,j,i), gaussb(3,j,i)
            endif
         enddo
      end do
      close(ofileno)

      return

      write (*,*) "Error writing the static gaussian shaped field file"

      10006 format (i8,i8,3es18.9E3)
      10007 format (es16.8,i8,3es18.9E3)

   end subroutine print_gauss

   !> Buffer the moving static gaussian shaped field
   subroutine buffer_mov_gauss(Natom, mstep, bcount_mov_gauss,delta_t,&
         real_time_measure,mov_gauss_spatial,mov_gauss)

      implicit none

      integer, intent(in) :: mstep         !< Current simulation step
      integer, intent(in) :: Natom         !< Number of atoms in system
      integer, intent(in) :: bcount_mov_gauss  !< Counter for the buffer of the moving static gaussian shaped field
      real(dblprec), intent(in) :: delta_t !< Current time step (used for real time measurements)
      real(dblprec), dimension(3,Natom), intent(in) :: mov_gauss_spatial !< Array with the moving gaussian shaped static field
      character(len=1), intent(in) :: mov_gauss         !< Add the moving static gaussian shaped field (Y/P/N)
      character(len=1), intent(in) :: real_time_measure !< Real time measurement flag

      ! .. Local variables
      integer :: i

      if (mov_gauss=='Y'.or.mov_gauss=='P') then
         do i=1, Natom
            mov_gaussb(1:3,i,bcount_mov_gauss)=mov_gauss_spatial(1:3,i)
         end do
      endif

      if (real_time_measure=='Y') then
         indxb_mov_gauss(bcount_mov_gauss)=mstep*delta_t
      else
         indxb_mov_gauss(bcount_mov_gauss)=mstep
      endif

   end subroutine buffer_mov_gauss

   !> Print the moving static gaussian shaped field
   subroutine print_mov_gauss(Natom,simid,real_time_measure)
      !
      !.. Implicit declarations
      implicit none

      integer, intent(in) :: Natom          !< Number of atoms in system
      character(len=8), intent(in) :: simid !< Name of simulation
      character(len=1), intent(in) :: real_time_measure !< Measurements displayed in real time

      !.. Local variables
      integer :: i, j
      character(len=30) :: filn

      ! Print to file the moving static gaussian shaped field
      write (filn,'(''mov_gauss.'',a8,''.out'')') simid
      open(ofileno, file=filn, position="append")
      do i=1, bcount_mov_gauss
         do j=1, Natom
            if (real_time_measure=='Y') then
               write (ofileno,10007) indxb_mov_gauss(i),j, mov_gaussb(1,j,i), mov_gaussb(2,j,i), mov_gaussb(3,j,i)
            else
               write (ofileno,10006) int(indxb_mov_gauss(i)),j, mov_gaussb(1,j,i), mov_gaussb(2,j,i), mov_gaussb(3,j,i)
            endif
         enddo
      end do
      close(ofileno)

      return

      write (*,*) "Error writing the gaussian shaped static moving field file"

      10006 format (i8,i8,3es18.9E3)
      10007 format (es16.8,i8,3es18.9E3)

   end subroutine print_mov_gauss

   !> Buffer the moving static circular shaped field
   subroutine buffer_mov_circle(Natom, mstep, bcount_mov_circle,delta_t,&
         real_time_measure,mov_circle_spatial,mov_circle)

      implicit none

      integer, intent(in) :: mstep         !< Current simulation step
      integer, intent(in) :: Natom         !< Number of atoms in system
      integer, intent(in) :: bcount_mov_circle  !< Counter for the buffer of the moving static circular shaped field
      real(dblprec), intent(in) :: delta_t !< Current time step (used for real time measurements)
      real(dblprec), dimension(3,Natom), intent(in) :: mov_circle_spatial !< Array with the moving circular shaped static field
      character(len=1), intent(in) :: mov_circle        !< Add the moving static circular shaped field (Y/P/N)
      character(len=1), intent(in) :: real_time_measure !< Real time measurement flag

      ! .. Local variables
      integer :: i

      if (mov_circle=='Y'.or.mov_circle=='P') then
         do i=1, Natom
            mov_circleb(1:3,i,bcount_mov_circle)=mov_circle_spatial(1:3,i)
         end do
      endif

      if (real_time_measure=='Y') then
         indxb_mov_circle(bcount_mov_circle)=mstep*delta_t
      else
         indxb_mov_circle(bcount_mov_circle)=mstep
      endif

   end subroutine buffer_mov_circle

   !> Print the moving static circular shaped field
   subroutine print_mov_circle(Natom,simid,real_time_measure)
      !
      !.. Implicit declarations
      implicit none

      integer, intent(in) :: Natom          !< Number of atoms in system
      character(len=8), intent(in) :: simid !< Name of simulation
      character(len=1), intent(in) :: real_time_measure !< Measurements displayed in real time

      !.. Local variables
      integer :: i, j
      character(len=30) :: filn

      ! Print to file the moving static circular shaped field
      write (filn,'(''mov_circle.'',a8,''.out'')') simid
      open(ofileno, file=filn, position="append")
      do i=1, bcount_mov_circle
         do j=1, Natom
            if (real_time_measure=='Y') then
               write (ofileno,10007) indxb_mov_circle(i),j, mov_circleb(1,j,i), mov_circleb(2,j,i), mov_circleb(3,j,i)
            else
               write (ofileno,10006) int(indxb_mov_circle(i)),j, mov_circleb(1,j,i), mov_circleb(2,j,i), mov_circleb(3,j,i)
            endif
         enddo
      end do
      close(ofileno)

      return

      write (*,*) "Error writing the circular shaped static moving field file"

      10006 format (i8,i8,3es18.9E3)
      10007 format (es16.8,i8,3es18.9E3)

   end subroutine print_mov_circle

   !> Buffer the moving microwave gaussian shaped field
   subroutine buffer_mwf_mov_gauss(Natom, mstep, bcount_mwf_mov_gauss,delta_t,&
         real_time_measure,mwf_mov_gauss_spatial,mwf_mov_gauss)

      implicit none

      integer, intent(in) :: mstep         !< Current simulation step
      integer, intent(in) :: Natom         !< Number of atoms in system
      integer, intent(in) :: bcount_mwf_mov_gauss  !< Counter of buffer for the moving microwave gaussian shaped field
      real(dblprec), intent(in) :: delta_t !< Current time step (used for real time measurements)
      real(dblprec), dimension(3,Natom), intent(in) :: mwf_mov_gauss_spatial !< Array with the moving gaussian shaped static field
      character(len=1), intent(in) :: mwf_mov_gauss         !< Add the moving static gaussian shaped field (Y/P/N)
      character(len=1), intent(in) :: real_time_measure !< Real time measurement flag

      !.. Local variables
      integer :: i

      if (mwf_mov_gauss=='Y'.or.mwf_mov_gauss=='P') then
         do i=1, Natom
            mwf_mov_gaussb(1:3,i,bcount_mwf_mov_gauss)=mwf_mov_gauss_spatial(1:3,i)
         end do
      endif

      if (real_time_measure=='Y') then
         indxb_mwf_mov_gauss(bcount_mwf_mov_gauss)=mstep*delta_t
      else
         indxb_mwf_mov_gauss(bcount_mwf_mov_gauss)=mstep
      endif

   end subroutine buffer_mwf_mov_gauss

   !> Print the moving microwave gaussian shaped field
   subroutine print_mwf_mov_gauss(Natom, simid,real_time_measure)
      !
      !.. Implicit declarations
      implicit none

      integer, intent(in) :: Natom          !< Number of atoms in system
      character(len=8), intent(in) :: simid !< Name of simulation
      character(len=1), intent(in) :: real_time_measure !< Measurements displayed in real time

      !.. Local variables
      integer :: i, j
      character(len=30) :: filn

      ! Print to file the moving microwave gaussian shaped field
      write (filn,'(''mwf_mov_gauss.'',a8,''.out'')') simid
      open(ofileno, file=filn, position="append")
      do i=1, bcount_mwf_mov_gauss
         do j=1, Natom
            if (real_time_measure=='Y') then
               write (ofileno,10007) indxb_mwf_mov_gauss(i),j, mwf_mov_gaussb(1,j,i), mwf_mov_gaussb(2,j,i), mwf_mov_gaussb(3,j,i)
            else
               write (ofileno,10006) int(indxb_mwf_mov_gauss(i)),j, mwf_mov_gaussb(1,j,i), mwf_mov_gaussb(2,j,i), mwf_mov_gaussb(3,j,i)
            endif
         enddo
      end do
      close(ofileno)

      return

      write (*,*) "Error writing the gaussian shaped static moving field file"

      10006 format (i8,i8,3es18.9E3)
      10007 format (es16.8,i8,3es18.9E3)

   end subroutine print_mwf_mov_gauss

   !> Buffer the moving microwave circular shaped field
   subroutine buffer_mwf_mov_circle(Natom, mstep, bcount_mwf_mov_circle,delta_t,&
         real_time_measure,mwf_mov_circle_spatial,mwf_mov_circle)

      implicit none

      integer, intent(in) :: mstep         !< Current simulation step
      integer, intent(in) :: Natom         !< Number of atoms in system
      integer, intent(in) :: bcount_mwf_mov_circle  !< Counter of buffer for the moving microwave circular shaped field
      real(dblprec), intent(in) :: delta_t !< Current time step (used for real time measurements)
      real(dblprec), dimension(3,Natom), intent(in) :: mwf_mov_circle_spatial !< Array with the moving circular shaped static field
      character(len=1), intent(in) :: mwf_mov_circle    !< Add the moving static circular shaped field (Y/P/N)
      character(len=1), intent(in) :: real_time_measure !< Real time measurement flag

      !.. Local variables
      integer :: i

      if (mwf_mov_circle=='Y'.or.mwf_mov_circle=='P') then
         do i=1, Natom
            mwf_mov_circleb(1:3,i,bcount_mwf_mov_circle)=mwf_mov_circle_spatial(1:3,i)
         end do
      endif

      if (real_time_measure=='Y') then
         indxb_mwf_mov_circle(bcount_mwf_mov_circle)=mstep*delta_t
      else
         indxb_mwf_mov_circle(bcount_mwf_mov_circle)=mstep
      endif

   end subroutine buffer_mwf_mov_circle

   !> Print the moving microwave circular shaped field
   subroutine print_mwf_mov_circle(Natom, simid,real_time_measure)
      !
      !.. Implicit declarations
      implicit none

      integer, intent(in) :: Natom          !< Number of atoms in system
      character(len=8), intent(in) :: simid !< Name of simulation
      character(len=1), intent(in) :: real_time_measure !< Measurements displayed in real time

      !.. Local variables
      integer :: i, j
      character(len=30) :: filn

      ! Print to file the moving microwave circular shaped field
      write (filn,'(''mwf_mov_circle.'',a8,''.out'')') simid
      open(ofileno, file=filn, position="append")
      do i=1, bcount_mwf_mov_circle
         do j=1, Natom
            if (real_time_measure=='Y') then
               write (ofileno,10007) indxb_mwf_mov_circle(i),j, mwf_mov_circleb(1,j,i), mwf_mov_circleb(2,j,i), mwf_mov_circleb(3,j,i)
            else
               write (ofileno,10006) int(indxb_mwf_mov_circle(i)),j, mwf_mov_circleb(1,j,i), mwf_mov_circleb(2,j,i), mwf_mov_circleb(3,j,i)
            endif
         enddo
      end do
      close(ofileno)

      return

      write (*,*) "Error writing the circular shaped static moving field file"

      10006 format (i8,i8,3es18.9E3)
      10007 format (es16.8,i8,3es18.9E3)

   end subroutine print_mwf_mov_circle

   !> Buffer the moving static cubic shaped field
   subroutine buffer_mov_square(Natom, mstep, bcount_mov_square,delta_t,&
         real_time_measure,mov_square_spatial,mov_square)

      implicit none

      integer, intent(in) :: mstep         !< Current simulation step
      integer, intent(in) :: Natom         !< Number of atoms in system
      integer, intent(in) :: bcount_mov_square  !< Counter for the buffer of the moving static cubic shaped field
      real(dblprec), intent(in) :: delta_t !< Current time step (used for real time measurements)
      real(dblprec), dimension(3,Natom), intent(in) :: mov_square_spatial !< Array with the moving cubic shaped static field
      character(len=1), intent(in) :: mov_square        !< Add the moving static cubic shaped field (Y/P/N)
      character(len=1), intent(in) :: real_time_measure !< Real time measurement flag

      ! .. Local variables
      integer :: i

      if (mov_square=='Y'.or.mov_square=='P') then
         do i=1, Natom
            mov_squareb(1:3,i,bcount_mov_square)=mov_square_spatial(1:3,i)
         end do
      endif

      if (real_time_measure=='Y') then
         indxb_mov_square(bcount_mov_square)=mstep*delta_t
      else
         indxb_mov_square(bcount_mov_square)=mstep
      endif

   end subroutine buffer_mov_square

   !> Print the moving static cubic shaped field
   subroutine print_mov_square(Natom,simid,real_time_measure)
      !
      !.. Implicit declarations
      implicit none

      integer, intent(in) :: Natom          !< Number of atoms in system
      character(len=8), intent(in) :: simid !< Name of simulation
      character(len=1), intent(in) :: real_time_measure !< Measurements displayed in real time

      !.. Local variables
      integer :: i, j
      character(len=30) :: filn

      ! Print to file the moving static cubic shaped field
      write (filn,'(''mov_square.'',a8,''.out'')') simid
      open(ofileno, file=filn, position="append")
      do i=1, bcount_mov_square
         do j=1, Natom
            if (real_time_measure=='Y') then
               write (ofileno,10007) indxb_mov_square(i),j, mov_squareb(1,j,i), mov_squareb(2,j,i), mov_squareb(3,j,i)
            else
               write (ofileno,10006) int(indxb_mov_square(i)),j, mov_squareb(1,j,i), mov_squareb(2,j,i), mov_squareb(3,j,i)
            endif
         enddo
      end do
      close(ofileno)

      return

      write (*,*) "Error writing the cubic shaped static moving field file"

      10006 format (i8,i8,3es18.9E3)
      10007 format (es16.8,i8,3es18.9E3)

   end subroutine print_mov_square

   !> Buffer the moving microwave cubic shaped field
   subroutine buffer_mwf_mov_square(Natom, mstep, bcount_mwf_mov_square,delta_t,&
         real_time_measure,mwf_mov_square_spatial,mwf_mov_square)

      implicit none

      integer, intent(in) :: mstep         !< Current simulation step
      integer, intent(in) :: Natom         !< Number of atoms in system
      integer, intent(in) :: bcount_mwf_mov_square  !< Counter of buffer for the moving microwave cubic shaped field
      real(dblprec), intent(in) :: delta_t !< Current time step (used for real time measurements)
      real(dblprec), dimension(3,Natom), intent(in) :: mwf_mov_square_spatial !< Array with the moving cubic shaped static field
      character(len=1), intent(in) :: mwf_mov_square    !< Add the moving static cubic shaped field (Y/P/N)
      character(len=1), intent(in) :: real_time_measure !< Real time measurement flag

      !.. Local variables
      integer :: i

      if (mwf_mov_square=='Y'.or.mwf_mov_square=='P') then
         do i=1, Natom
            mwf_mov_squareb(1:3,i,bcount_mwf_mov_square)=mwf_mov_square_spatial(1:3,i)
         end do
      endif

      if (real_time_measure=='Y') then
         indxb_mwf_mov_square(bcount_mwf_mov_square)=mstep*delta_t
      else
         indxb_mwf_mov_square(bcount_mwf_mov_square)=mstep
      endif

   end subroutine buffer_mwf_mov_square

   !> Print the moving microwave cubic shaped field
   subroutine print_mwf_mov_square(Natom, simid,real_time_measure)
      !
      !.. Implicit declarations
      implicit none

      integer, intent(in) :: Natom          !< Number of atoms in system
      character(len=8), intent(in) :: simid !< Name of simulation
      character(len=1), intent(in) :: real_time_measure !< Measurements displayed in real time

      !.. Local variables
      integer :: i, j
      character(len=30) :: filn

      ! Print to file the moving microwave cubic shaped field
      write (filn,'(''mwf_mov_square.'',a8,''.out'')') simid
      open(ofileno, file=filn, position="append")
      do i=1, bcount_mwf_mov_square
         do j=1, Natom
            if (real_time_measure=='Y') then
               write (ofileno,10007) indxb_mwf_mov_square(i),j, mwf_mov_squareb(1,j,i), mwf_mov_squareb(2,j,i), mwf_mov_squareb(3,j,i)
            else
               write (ofileno,10006) int(indxb_mwf_mov_square(i)),j, mwf_mov_squareb(1,j,i), mwf_mov_squareb(2,j,i), mwf_mov_squareb(3,j,i)
            endif
         enddo
      end do
      close(ofileno)

      return

      write (*,*) "Error writing the cubic shaped static moving field file"

      10006 format (i8,i8,3es18.9E3)
      10007 format (es16.8,i8,3es18.9E3)

   end subroutine print_mwf_mov_square

end module prn_microwaves
