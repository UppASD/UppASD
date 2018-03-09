!> Data and routines for calculating applied microwave field
!> @copyright
!! Copyright (C) 2008-2018 UppASD group
!! This file is distributed under the terms of the
!! GNU General Public License. 
!! See http://www.gnu.org/copyleft/gpl.txt
module MicroWaveField
   use Parameters
   use Profiling
   use Constants
   !
   implicit none
   !
   !Micro-wave field
   !From input
   !From simulation
   integer :: centering
   integer :: mwf_centering
   integer :: Num_mov_centers       !< Trajectory for the moving gaussian field
   integer :: Num_gauss_centers     !< Number of centers for the gaussian shaped field
   integer :: Num_mwf_mov_centers   !< Trajectory for the microwave gaussian shaped moving field
   integer :: Num_mwf_gauss_centers !< Number of centers for the gaussian shaped microwave field
   real(dblprec) :: mwftime !< Current time for microwave field evaluation
   real(dblprec), dimension(3) :: mwffield       !< Current time dependent microwave field
   real(dblprec), dimension(3) :: gauss_mwffield !< Current gaussian broadened microwave field
   integer, dimension(:), allocatable :: mov_site_center        !< Array with the trajectory of the static gaussian field
   integer, dimension(:), allocatable :: gauss_site_center      !< Array with the centers of the gaussian static fields
   integer, dimension(:), allocatable :: mwf_mov_site_center    !< Array with the trajectory of the time dependent gaussian field
   integer, dimension(:), allocatable :: mwf_gauss_site_center  !< Array with the centers of the gaussian shaped microwave field
   real(dblprec), dimension(:), allocatable :: field_site_phase       !< Array containing the site dependent phase
   real(dblprec), dimension(:,:), allocatable :: gauss_sites          !< Array containing the centers of the gaussian static field
   real(dblprec), dimension(:,:), allocatable :: gauss_spatial        !< Array with the static gaussian shaped field
   real(dblprec), dimension(:,:), allocatable :: site_mwffield        !< Array with the monochromatic microwave field
   real(dblprec), dimension(:,:), allocatable :: mwffield_sites       !< Array containing the atoms in which the monochromatic microwave field is "on"
   real(dblprec), dimension(:,:), allocatable :: mov_gauss_sites      !< Array containing the trajectory of the static gaussian shaped field
   real(dblprec), dimension(:,:), allocatable :: mov_square_sites     !< Array containing the trajectory of the static cubic shaped field
   real(dblprec), dimension(:,:), allocatable :: mov_circle_sites     !< Array containing the trajectory of the static circular shaped field
   real(dblprec), dimension(:,:), allocatable :: mov_gauss_spatial    !< Array with the moving gaussian shaped static field
   real(dblprec), dimension(:,:), allocatable :: mov_square_spatial   !< Array with the moving cubic shaped static field
   real(dblprec), dimension(:,:), allocatable :: mov_circle_spatial   !< Array with the moving circular shaped static field
   real(dblprec), dimension(:,:), allocatable :: gauss_site_mwffield  !< Array with the gaussian broadened microwave field
   real(dblprec), dimension(:,:), allocatable :: gauss_mwffield_sites !< Array containing the atoms in which the gaussian broadened microwave field is "on"
   real(dblprec), dimension(:,:), allocatable :: mwf_mov_gauss_spatial     !< Array with the moving microwave gaussian shaped field
   real(dblprec), dimension(:,:), allocatable :: mwf_mov_square_spatial    !< Array with the moving microwave cubic shaped field
   real(dblprec), dimension(:,:), allocatable :: mwf_mov_circle_spatial    !< Array with the moving microwave circular shaped field
   real(dblprec), dimension(:,:), allocatable :: mov_gauss_mwffield_sites  !< Array containing the trajectory of the gaussian shaped microwave field
   real(dblprec), dimension(:,:), allocatable :: mov_square_mwffield_sites !< Array containing the trajectory of the cubic shaped microwave field
   real(dblprec), dimension(:,:), allocatable :: mov_circle_mwffield_sites !< Array containing the trajectory of the circular shaped microwave field
   real(dblprec), dimension(:,:), allocatable :: gauss_spatial_site_mwffield  !< Array with the gaussian broadened spatially gaussian microwave field
   real(dblprec), dimension(:,:), allocatable :: gauss_spatial_mwffield_sites !< Array containing the atoms where gaussian broadened spatially gaussian microwave field is "on"
   real(dblprec), dimension(:,:), allocatable :: time_external_fields
   ! Buffers and step sizes
   integer :: mwf_buff               !< Buffer size for the monochromatic microwave field
   integer :: mwf_step               !< Interval for sampling the monochromatic microwave field
   integer :: gauss_buff             !< Buffer size for the static gaussian shaped field
   integer :: gauss_step             !< Interval for sampling the static gaussian shaped field
   integer :: mwf_gauss_buff         !< Buffer size for the frequency broadened microwave field
   integer :: mwf_gauss_step         !< Interval for sampling the frequency broadened microwave field
   integer :: mov_gauss_buff         !< Buffer size for the moving static gaussian shaped field
   integer :: mov_circle_buff        !< Buffer size for the moving static circular shaped field
   integer :: mov_square_buff        !< Buffer size for the moving static cubic shaped field
   integer :: mov_gauss_pstep        !< Interval for sampling the moving static gaussian shaped field
   integer :: mov_circle_pstep       !< Interval for sampling the moving static circular shaped field
   integer :: mov_square_pstep       !< Interval for sampling the moving static cubic shaped field
   integer :: mwf_mov_gauss_buff     !< Buffer size for the moving gaussian shaped microwave field
   integer :: mwf_mov_circle_buff    !< Buffer size for the moving circular shaped microwave field
   integer :: mwf_mov_square_buff    !< Buffer size for the moving cubic shaped microwave field
   integer :: mwf_mov_gauss_pstep    !< Interval for sampling the moving gaussian shaped microwave field
   integer :: mwf_mov_circle_pstep   !< Interval for sampling the moving circular shaped microwave field
   integer :: mwf_mov_square_pstep   !< Interval for sampling the moving cubic shaped microwave field
   integer :: mwf_gauss_spatial_buff !< Buffer size for the gaussian shaped frequency broadened microwave field
   integer :: mwf_gauss_spatial_step !< Interval for sampling the gaussian shaped frequency broadened microwave field

   ! Measurement flags
   character(len=1) :: prn_mwf               !< Flag for printing the monochromatic microwave field
   character(len=1) :: prn_gauss             !< Flag for printing the static gaussian shaped field
   character(len=1) :: prn_mwf_gauss         !< Flag for printing the frequency broadened microwave field
   character(len=1) :: prn_mov_gauss         !< Flag for printing the moving static gaussian shaped field
   character(len=1) :: prn_mov_circle        !< Flag for printing the moving static circular shaped field
   character(len=1) :: prn_mov_square        !< Flag for printing the moving static cubic shaped field
   character(len=1) :: prn_mwf_mov_gauss     !< Flag for printing the moving gaussian shaped microwave field
   character(len=1) :: prn_mwf_mov_circle    !< Flag for printing the moving circular shaped microwave field
   character(len=1) :: prn_mwf_mov_square    !< Flag for printing the moving cubic shaped microwave field
   character(len=1) :: prn_mwf_gauss_spatial !< Flag for printing the gaussian shaped frequency broadened microwave field

   !Monochromatic Microwave field
   character(len=1) :: mwf              !< Add monochromatic microwave field (Y/P/S/W/N)
   character(len=1) :: site_phase       !< Add site dependent phase to the microwavefield
   character(len=35) :: mwf_site_file   !< File name for the site dependent monochromatic micorwave field
   character(len=35) :: site_phase_file !< File name for the site dependent phase
   integer :: mwf_pulse_time !< Number of time steps in which the monochormatic microwave field is on
   real(dblprec) :: mwfampl  !< Amplitude of the monochormatic microwave field
   real(dblprec) :: mwffreq  !< Frequency of the monochormatic microwave field
   real(dblprec), dimension(3) :: mwfdir  !< Direction of the monochromatic microwave field

   !Gaussian frequency broadened Microwave fields
   character(len=1) :: mwf_gauss            !< Add frequency broadened microwave field (Y/P/S/W/N)
   character(len=35) :: mwf_gauss_site_file !< File name for the site dependent frequency broadened microwave field
   integer :: mwf_gauss_pulse_time !< Number of time steps in which the frequency broadened microwave field is on
   real(dblprec) :: mwf_gauss_ampl !< Amplitude of the frequency broadened microwave field
   real(dblprec) :: mwf_gauss_freq !< Frequency of the frequency broadened microwave field
   real(dblprec) :: mwf_gauss_time_sigma        !< Frequency sigma parameter for the frequency broadened microwave field
   real(dblprec), dimension(3) :: mwf_gauss_dir !< Direction of the frequency broadened microwave field

   !Spatial Gaussian shaped frequency broadened microwave fields
   character(len=1) :: mwf_gauss_spatial            !< Add the frequency broadened gaussian shaped microwave field (Y/P/N)
   character(len=35) :: mwf_gauss_spatial_site_file !< File name for the site dependent gaussian shaped gaussian broadened microwave field
   integer :: mwf_gauss_spatial_pulse_time !< Number of time steps in which the frequency broadened gaussian shaped microwave field is on
   real(dblprec) :: mwf_gauss_spatial_freq !< Frequency for the frequency broadened gaussian shaped microwave field
   real(dblprec) :: mwf_gauss_spatial_ampl !< Amplitude for the frequency broadened gaussian shaped microwave field
   real(dblprec) :: mwf_gauss_spatial_time_sigma !< Frequency sigma parameter for the frequency broadened gaussian shaped microwave field
   real(dblprec), dimension(3) :: mwf_gauss_spatial_space_sigma !< Spatial sigma parameter for the frequency broadened gaussian shaped microwave field

   ! Static Gaussian shaped microwave field
   character(len=1) :: do_gauss          !< Add the Static gaussian shaped field
   character(len=35) :: gauss_site_file  !< File name for the site dependent static gaussian shaped field
   integer :: gauss_pulse_time           !< Number of time steps in which the static gaussian shaped pulse is on
   real(dblprec) :: gauss_spatial_ampl   !< Amplitude of the gaussian shaped static field
   real(dblprec), dimension(3) :: gauss_spatial_sigma !< Sigma parameter for the gaussian shaped static field

   !Moving gaussian shaped static field
   character(len=1) :: mov_gauss        !< Add the moving static gaussian shaped field (Y/P/N)
   character(len=35) :: mov_gauss_file  !< Moving static gaussian shaped field trajectory file
   integer :: mov_gauss_step            !< Number of time steps in which the positions of the moving static gaussian shaped pulse is updated
   integer :: mov_gauss_pulse_time      !< Number of time steps in which the moving static gaussian shaped pulse is on
   real(dblprec) :: mov_gauss_ampl      !< Amplitude of the moving static gaussian shaped pulse
   real(dblprec), dimension(3) :: mov_gauss_space_sigma !< Sigma parameter for the moving gaussian shaped static field

   !Moving gaussian shaped microwave field
   character(len=1) :: mwf_mov_gauss       !< Add the moving microwave gaussian shaped field (Y/P/N)
   character(len=35) :: mwf_mov_gauss_file !< Moving frequency broadened gaussian shaped microwave filed trajectory file
   integer :: mwf_mov_gauss_step           !< Number of time steps in which the positions of the moving microwave gaussian shaped pulse is updated
   integer :: mwf_mov_gauss_pulse_time     !< Number of time steps in which the moving microwave gaussian shaped pulse is on
   real(dblprec) :: mwf_mov_gauss_ampl     !< Amplitude of the moving microwave gaussian shaped pulse
   real(dblprec) :: mwf_mov_gauss_freq     !< Frequency of the moving microwave gaussian shaped pulse
   real(dblprec) :: mwf_mov_gauss_time_sigma !< Sigma parameter for frequency gaussian in the moving gaussian shaped microwave field
   real(dblprec), dimension(3) :: mwf_mov_gauss_space_sigma !< Sigma parameter for the moving gaussian shaped microwave field

   !Moving circular shaped static field
   character(len=1) :: mov_circle        !< Add the moving static circular shaped field (Y/P/N)
   character(len=35) :: mov_circle_file  !< Moving static circular shaped field trajectory file
   integer :: mov_circle_step            !< Number of time steps in which the positions of the moving static circular shaped pulse is updated
   integer :: mov_circle_pulse_time      !< Number of time steps in which the moving static circular shaped pulse is on
   real(dblprec) :: mov_circle_ampl      !< Amplitude of the moving static circular shaped pulse
   real(dblprec) :: mov_circle_radius    !< Radius for the moving circular shaped static field

   !Moving circular shaped microwave field
   character(len=1) :: mwf_mov_circle         !< Add the moving microwave circular shaped field (Y/P/N)
   character(len=35) :: mwf_mov_circle_file   !< Moving frequency broadened circular shaped microwave filed trajectory file
   integer :: mwf_mov_circle_step             !< Number of time steps in which the positions of the moving microwave circular shaped pulse is updated
   integer :: mwf_mov_circle_pulse_time       !< Number of time steps in which the moving microwave circular shaped pulse is on
   real(dblprec) :: mwf_mov_circle_ampl       !< Amplitude of the moving microwave circular shaped pulse
   real(dblprec) :: mwf_mov_circle_freq       !< Frequency of the moving microwave circular shaped pulse
   real(dblprec) :: mwf_mov_circle_time_sigma !< Sigma parameter for frequency gaussian in the moving circular shaped microwave field
   real(dblprec) :: mwf_mov_circle_radius !< Sigma parameter for the moving circular shaped microwave field

   !Moving cubic shaped static field
   character(len=1) :: mov_square        !< Add the moving static cubic shaped field (Y/P/N)
   character(len=35) :: mov_square_file  !< Moving static cubic shaped field trajectory file
   integer :: mov_square_step            !< Number of time steps in which the positions of the moving static cubic shaped pulse is updated
   integer :: mov_square_pulse_time      !< Number of time steps in which the moving static cubic shaped pulse is on
   real(dblprec) :: mov_square_ampl      !< Amplitude of the moving static cubic shaped pulse
   real(dblprec), dimension(3) :: mov_square_dimensions    !< Dimensions for the moving cubic shaped static field

   !Moving cubic shaped microwave field
   character(len=1) :: mwf_mov_square         !< Add the moving microwave cubic shaped field (Y/P/N)
   character(len=35) :: mwf_mov_square_file   !< Moving frequency broadened cubic shaped microwave filed trajectory file
   integer :: mwf_mov_square_step             !< Number of time steps in which the positions of the moving microwave cubic shaped pulse is updated
   integer :: mwf_mov_square_pulse_time       !< Number of time steps in which the moving microwave cubic shaped pulse is on
   real(dblprec) :: mwf_mov_square_ampl       !< Amplitude of the moving microwave cubic shaped pulse
   real(dblprec) :: mwf_mov_square_freq       !< Frequency of the moving microwave cubic shaped pulse
   real(dblprec) :: mwf_mov_square_time_sigma !< Sigma parameter for frequency gaussian in the moving cubic shaped microwave field
   real(dblprec), dimension(3) :: mwf_mov_square_dimensions !< Dimensions for the moving cubic shaped microwave field
   !
   public

contains

   !> Wrapper routine to calculate all the microwave related fields (and the gaussian shaped pulsed field)
   ! All the fields present in the routine can be used at the same time, independently from each other
   subroutine calculate_mwf_fields(Natom,time,maxtime,delta_t,coord,flag)

      implicit none

      ! Simulation variables
      integer, intent(in) :: time       !< Current simulation step
      integer, intent(in) :: Natom      !< Number of atoms in the system
      integer, intent(in) :: maxtime    !< Total number of simulation steps
      real(dblprec), intent(in) :: delta_t       !< Current time step
      real(dblprec), dimension(3,Natom) :: coord !< Coordinates for all the atoms

      integer, intent(in) :: flag

      !.. Executable part

      ! Calculate microwave field
      if(mwf=='Y'.or.mwf=='P'.or.mwf=='I') then
         ! Global monochromatic microwave field
         call calc_mwf(mwfampl, mwfdir, mwffreq, delta_t, time, mwf_pulse_time, mwf)
      else if (mwf=='S'.or.mwf=='W') then
         ! Site dependent monochromatic microwave field
         call calc_site_mwf(Natom,mwfampl, mwffreq, delta_t, time, mwf_pulse_time, mwf)
      endif

      ! Caculate a gaussian microwave field
      if (mwf_gauss=='Y'.or.mwf_gauss=='P') then
         ! Global frequency broadened microwave field
         call calc_gaussian_mwf(mwf_gauss_ampl,mwf_gauss_freq,delta_t,&
            mwf_gauss_time_sigma,mwf_gauss_dir,time,maxtime,mwf_gauss_pulse_time,mwf_gauss)
      else if (mwf_gauss=='S'.or.mwf_gauss=='W') then
         ! Site dependent frequency broadened microwave field
         call calc_site_mwf_gauss(Natom,mwf_gauss_ampl, mwf_gauss_freq, delta_t, time, maxtime,&
            mwf_gauss_pulse_time,mwf_gauss_time_sigma,mwf_gauss)
      endif

      ! Calculate a microwave gaussian field and a gaussian spatial distribution
      if (mwf_gauss_spatial=='Y'.or.mwf_gauss_spatial=='P') then
         ! Frequency broadened gaussian shaped microwave field
         call calc_spatial_mwf_site_gauss(Natom,coord,mwf_gauss_spatial,mwf_gauss_spatial_freq,&
            delta_t,time,maxtime,mwf_gauss_spatial_pulse_time,mwf_gauss_spatial_time_sigma,&
            mwf_gauss_spatial_space_sigma,mwf_gauss_spatial_ampl)
      endif

      ! Calculate a moving static spatially resolved gaussian magnetic pulse
      if (mov_gauss=='Y'.or.mov_gauss=='P') then
         ! See if the current center of the field must be updated
         if (centering.lt.Num_mov_centers) then
            if (mod(time-1,mov_gauss_step)==0.and.flag==0) then
               centering=centering+1
            endif
         endif
         ! Moving static gaussian shaped field
         call calc_moving_gauss(Natom,coord,mov_gauss_ampl,time,mov_gauss_pulse_time,centering,&
            mov_gauss_space_sigma,mov_gauss)
      endif

      ! Calculate a moving time dependent spatially resolved gaussian magnetic pulse
      if (mwf_mov_gauss=='Y'.or.mwf_mov_gauss=='P') then
         ! See if the current center of the field must be updated
         if (centering.lt.Num_mov_centers) then
            if (mod(time-1,mwf_mov_gauss_step)==0.and.flag==0)  then
               mwf_centering=mwf_centering+1
            endif
         endif
         ! Moving gaussian shaped microwave field
         call calc_moving_gauss_mwf(Natom,coord,mwf_mov_gauss_ampl,mwf_mov_gauss_freq,mwf_mov_gauss_time_sigma,&
            time,maxtime,delta_t,mwf_mov_gauss_pulse_time,mwf_centering,mwf_mov_gauss_space_sigma,mwf_mov_gauss)
      endif

      ! Calculate a moving static spatially resolved circular magnetic pulse
      if (mov_circle=='Y'.or.mov_circle=='P') then
         ! See if the current center of the field must be updated
         if (centering.lt.Num_mov_centers) then
            if (mod(time-1,mov_circle_step)==0.and.flag==0) then
               centering=centering+1
            endif
         endif
         ! Moving static gaussian shaped field
         call calc_moving_circle(Natom,coord,mov_circle_ampl,time,mov_circle_pulse_time,centering,&
            mov_circle_radius,mov_circle)

      endif

      ! Calculate a moving time dependent spatially resolved ciruclar magnetic pulse
      if (mwf_mov_circle=='Y'.or.mwf_mov_circle=='P') then
         ! See if the current center of the field must be updated
         if (centering.lt.Num_mov_centers) then
            if (mod(time-1,mwf_mov_circle_step)==0.and.flag==0)  then
               mwf_centering=mwf_centering+1
            endif
         endif
         ! Moving gaussian shaped microwave field
         call calc_moving_circle_mwf(Natom,coord,mwf_mov_circle_ampl,mwf_mov_circle_freq,mwf_mov_circle_time_sigma,&
            time,maxtime,delta_t,mwf_mov_circle_pulse_time,mwf_mov_circle_radius,mwf_mov_circle)
      endif

      ! Calculate a moving static spatially resolved cubic magnetic pulse
      if (mov_square=='Y'.or.mov_square=='P') then
         ! See if the current center of the field must be updated
         if (centering.lt.Num_mov_centers) then
            if (mod(time-1,mov_square_step)==0.and.flag==0) then
               centering=centering+1
            endif
         endif
         ! Moving static gaussian shaped field
         call calc_moving_square(Natom,coord,mov_square_ampl,time,mov_square_pulse_time,centering,&
            mov_square_dimensions,mov_square)

      endif

      ! Calculate a moving time dependent spatially resolved cubic magnetic pulse
      if (mwf_mov_square=='Y'.or.mwf_mov_square=='P') then
         ! See if the current center of the field must be updated
         if (centering.lt.Num_mov_centers) then
            if (mod(time-1,mwf_mov_square_step)==0.and.flag==0)  then
               mwf_centering=mwf_centering+1
            endif
         endif
         ! Moving gaussian shaped microwave field
         call calc_moving_square_mwf(Natom,coord,mwf_mov_square_ampl,mwf_mov_square_freq,mwf_mov_square_time_sigma,&
            time,maxtime,delta_t,mwf_mov_square_pulse_time,mwf_mov_square_dimensions,mwf_mov_square)
      endif

      ! Calculate the static gaussian shaped field
      if (do_gauss=='Y'.or.do_gauss=='P') then
         ! Gaussian shaped static field
         call calc_spatial_gauss(Natom,time,gauss_spatial_sigma,coord,gauss_spatial_ampl,delta_t,do_gauss,gauss_pulse_time)
      endif

   end subroutine calculate_mwf_fields

   !>  Calculate monochormatic global microwave field
   subroutine calc_mwf(mwfampl, mwfdir, mwffreq, delta_t,time,mwf_pulse_time,mwf)
      !

      !.. Implicit declatrations
      implicit none

      integer, intent(in) :: time           !< Current simulation step
      integer, intent(in) :: mwf_pulse_time !< Number of time steps in which the monochormatic microwave field is on
      real(dblprec), intent(in) :: mwfampl  !< Amplitude of the monochormatic microwave field
      real(dblprec), intent(in) :: mwffreq  !< Frequency of the monochormatic microwave field
      real(dblprec), intent(in) :: delta_t  !< Current time step
      real(dblprec), intent(in), dimension(3) :: mwfdir  !< Direction of the monochromatic microwave field
      character(len=1), intent(in) :: mwf   !< Add monochromatic microwave field (Y/P/S/W/N)

      mwftime=delta_t*(time-1)
      if(mwf=='Y') then
         mwffield(1)=mwfampl*mwfdir(1)*sin(mwffreq*mwftime*2*pi)
         mwffield(2)=mwfampl*mwfdir(2)*sin(mwffreq*mwftime*2*pi)
         mwffield(3)=mwfampl*mwfdir(3)*sin(mwffreq*mwftime*2*pi)
      else if (mwf=='P') then

         if(time<=mwf_pulse_time) then
            mwffield(1)=mwfampl*mwfdir(1)*sin(mwffreq*mwftime*2*pi)
            mwffield(2)=mwfampl*mwfdir(2)*sin(mwffreq*mwftime*2*pi)
            mwffield(3)=mwfampl*mwfdir(3)*sin(mwffreq*mwftime*2*pi)
         else
            mwffield=0.d0
         endif

      else if (mwf=='I') then

         if(time<=mwf_pulse_time) then
            mwffield(1)=mwfampl*mwfdir(1)*sin(mwffreq*mwftime*2*pi)**2
            mwffield(2)=mwfampl*mwfdir(2)*sin(mwffreq*mwftime*2*pi)**2
            mwffield(3)=mwfampl*mwfdir(3)*sin(mwffreq*mwftime*2*pi)**2
         else
            mwffield=0.d0
         endif
      end if

   end subroutine calc_mwf


   !> Calculate frequency broadened global microwave field
   subroutine calc_gaussian_mwf(mwf_gauss_ampl,mwf_gauss_freq,delta_t,mwf_gauss_time_sigma,&
         mwf_gauss_dir,time,maxtime,mwf_gauss_pulse_time,mwf_gauss)

      implicit none

      real(dblprec), intent(in) :: delta_t        !< Current time step
      real(dblprec), intent(in) :: mwf_gauss_ampl !< Amplitude of the frequency broadened microwave field
      real(dblprec), intent(in) :: mwf_gauss_freq !< Frequency of the frequency broadened microwave field
      real(dblprec), intent(in) :: mwf_gauss_time_sigma        !< Frequency sigma parameter for the frequency broadened microwave field
      real(dblprec), dimension(3), intent(in) :: mwf_gauss_dir !< Direction of the frequency broadened microwave field
      integer, intent(in) :: time    !< Current simulation step
      integer, intent(in) :: maxtime !< Total number of simulation steps
      integer, intent(in) :: mwf_gauss_pulse_time !< Number of time steps in which the frequency broadened microwave field is on
      character(len=1), intent(in) :: mwf_gauss   !< Add frequency broadened microwave field (Y/P/S/W/N)

      mwftime=delta_t*(time-1)

      if (mwf_gauss=='Y') then
         gauss_mwffield(1)=mwf_gauss_ampl*mwf_gauss_dir(1)*sin(mwf_gauss_freq*mwftime*2*pi)&
            *exp((-mwftime**2)/(2*(maxtime*delta_t*mwf_gauss_time_sigma)**2))
         gauss_mwffield(2)=mwf_gauss_ampl*mwf_gauss_dir(2)*sin(mwf_gauss_freq*mwftime*2*pi)&
            *exp((-mwftime**2)/(2*(maxtime*delta_t*mwf_gauss_time_sigma)**2))
         gauss_mwffield(3)=mwf_gauss_ampl*mwf_gauss_dir(3)*sin(mwf_gauss_freq*mwftime*2*pi)&
            *exp((-mwftime**2)/(2*(maxtime*delta_t*mwf_gauss_time_sigma)**2))

      else if (mwf_gauss=='P') then

         if(time<=mwf_gauss_pulse_time) then
            gauss_mwffield(1)=mwf_gauss_ampl*mwf_gauss_dir(1)*sin(mwf_gauss_freq*mwftime*2*pi)&
               *exp((-mwftime**2)/(2*(maxtime*delta_t*mwf_gauss_time_sigma)**2))
            gauss_mwffield(2)=mwf_gauss_ampl*mwf_gauss_dir(2)*sin(mwf_gauss_freq*mwftime*2*pi)&
               *exp((-mwftime**2)/(2*(maxtime*delta_t*mwf_gauss_time_sigma)**2))
            gauss_mwffield(3)=mwf_gauss_ampl*mwf_gauss_dir(3)*sin(mwf_gauss_freq*mwftime*2*pi)&
               *exp((-mwftime**2)/(2*(maxtime*delta_t*mwf_gauss_time_sigma)**2))
         else
            gauss_mwffield=0.d0
         endif

      endif

   end subroutine calc_gaussian_mwf

   !> Calculate site dependent monochromatic microwave field
   subroutine calc_site_mwf(Natom, mwfampl, mwffreq, delta_t, time, mwf_pulse_time, mwf)
      !

      !.. Implicit declatrations
      implicit none

      integer, intent(in) :: time  !< Current simulation step
      integer, intent(in) :: Natom !< Number of atoms in the system
      integer, intent(in) :: mwf_pulse_time !< Number of time steps in which the monochormatic microwave field is on
      real(dblprec), intent(in) :: mwfampl  !< Amplitude of the monochormatic microwave field
      real(dblprec), intent(in) :: mwffreq  !< Frequency of the monochormatic microwave field
      real(dblprec), intent(in) :: delta_t  !< Current time step
      character(len=1), intent(in) :: mwf   !< Add monochromatic microwave field (Y/P/S/W/N)

      !.. Local variables
      integer :: i

      mwftime=delta_t*(time-1)

      if(mwf=='S') then
         do i=1,Natom
            site_mwffield(1,i)=mwfampl*mwffield_sites(1,i)*sin(mwffreq*mwftime*2*pi+field_site_phase(i))
            site_mwffield(2,i)=mwfampl*mwffield_sites(2,i)*sin(mwffreq*mwftime*2*pi+field_site_phase(i))
            site_mwffield(3,i)=mwfampl*mwffield_sites(3,i)*sin(mwffreq*mwftime*2*pi+field_site_phase(i))
         enddo

      else if (mwf=='W') then
         if(time<=mwf_pulse_time) then
            do i=1,Natom
               site_mwffield(1,i)=mwfampl*mwffield_sites(1,i)*sin(mwffreq*mwftime*2*pi+field_site_phase(i))
               site_mwffield(2,i)=mwfampl*mwffield_sites(2,i)*sin(mwffreq*mwftime*2*pi+field_site_phase(i))
               site_mwffield(3,i)=mwfampl*mwffield_sites(3,i)*sin(mwffreq*mwftime*2*pi+field_site_phase(i))
            enddo
         else
            site_mwffield=0.d0
         endif

      end if

   end subroutine calc_site_mwf

   ! Site and time dependent gaussian magnetic field
   ! This function is a time dependent gaussian that only acts in a certain part of the sample, the shape of the volume where the field is acting is detemined by the input file
   subroutine calc_site_mwf_gauss(Natom,mwf_gauss_ampl, mwf_gauss_freq, delta_t, time, maxtime,mwf_gauss_pulse_time,&
         mwf_gauss_time_sigma,mwf_gauss)
      !

      !.. Implicit declatrations
      implicit none

      integer, intent(in) :: time    !< Current simulation step
      integer, intent(in) :: Natom   !< Number of atoms in the system
      integer, intent(in) :: maxtime !< Total number of simulation steps
      integer, intent(in) :: mwf_gauss_pulse_time !< Time of the oscillating pulse
      real(dblprec), intent(in) :: delta_t        !< Current time step
      real(dblprec), intent(in) :: mwf_gauss_ampl !< Amplitude of microwave field
      real(dblprec), intent(in) :: mwf_gauss_freq !< Frequency of microwave field
      real(dblprec), intent(in) :: mwf_gauss_time_sigma !< Frequency sigma parameter for the frequency broadened microwave field
      character(len=1), intent(in) :: mwf_gauss         !< Add frequency broadened microwave field (Y/P/S/W/N)

      !.. Local variables
      integer :: i

      mwftime=delta_t*(time-1)

      if (mwf_gauss=='S') then
         do i=1, Natom
            gauss_site_mwffield(1,i)=mwf_gauss_ampl*gauss_mwffield_sites(1,i)*sin(mwf_gauss_freq*mwftime*2*pi)&
               *exp((-mwftime**2)/(2*(maxtime*delta_t*mwf_gauss_time_sigma)**2))
            gauss_site_mwffield(2,i)=mwf_gauss_ampl*gauss_mwffield_sites(2,i)*sin(mwf_gauss_freq*mwftime*2*pi)&
               *exp((-mwftime**2)/(2*(maxtime*delta_t*mwf_gauss_time_sigma)**2))
            gauss_site_mwffield(3,i)=mwf_gauss_ampl*gauss_mwffield_sites(3,i)*sin(mwf_gauss_freq*mwftime*2*pi)&
               *exp((-mwftime**2)/(2*(maxtime*delta_t*mwf_gauss_time_sigma)**2))
         enddo

      else if (mwf_gauss=='W') then
         if(time<=mwf_gauss_pulse_time) then
            do i=1,Natom
               gauss_site_mwffield(1,i)=mwf_gauss_ampl*gauss_mwffield_sites(1,i)*sin(mwf_gauss_freq*mwftime*2*pi)&
                  *exp((-mwftime**2)/(2*(maxtime*delta_t*mwf_gauss_time_sigma)**2))
               gauss_site_mwffield(2,i)=mwf_gauss_ampl*gauss_mwffield_sites(2,i)*sin(mwf_gauss_freq*mwftime*2*pi)&
                  *exp((-mwftime**2)/(2*(maxtime*delta_t*mwf_gauss_time_sigma)**2))
               gauss_site_mwffield(3,i)=mwf_gauss_ampl*gauss_mwffield_sites(3,i)*sin(mwf_gauss_freq*mwftime*2*pi)&
                  *exp((-mwftime**2)/(2*(maxtime*delta_t*mwf_gauss_time_sigma)**2))
            enddo
         else
            gauss_site_mwffield=0.d0
         endif

      end if

   end subroutine calc_site_mwf_gauss


   ! The one can calculate a time dependent, position dependent gaussian magnetic pulse
   ! In the input file the atom number for the center of the guassian is given, as well as the direction of the magnetic field
   ! The gaussian distribution only affects the intensity it does not affect the direction of the magnetization
   subroutine calc_spatial_mwf_site_gauss(Natom,coord,mwf_gauss_spatial,mwf_gauss_spatial_freq,delta_t, &
         time,maxtime,mwf_gauss_spatial_pulse_time,mwf_gauss_spatial_time_sigma,mwf_gauss_spatial_space_sigma,mwf_gauss_spatial_ampl)

      implicit none

      integer, intent(in) :: time    !< Current simulation step
      integer, intent(in) :: Natom   !< Number of atoms in the system
      integer, intent(in) :: maxtime !< Total number of simulation steps
      integer, intent(in) :: mwf_gauss_spatial_pulse_time !< Number of time steps in which the frequency broadened gaussian shaped microwave field is on
      real(dblprec), intent(in) :: delta_t                !< Current time step
      real(dblprec), intent(in) :: mwf_gauss_spatial_freq !< Frequency for the frequency broadened gaussian shaped microwave field
      real(dblprec), intent(in) :: mwf_gauss_spatial_ampl !< Amplitude for the frequency broadened gaussian shaped microwave field
      real(dblprec), intent(in) :: mwf_gauss_spatial_time_sigma !< Frequency sigma parameter for the frequency broadened gaussian shaped microwave field
      real(dblprec), dimension(3), intent(in) :: mwf_gauss_spatial_space_sigma !< Spatial sigma parameter for the frequency broadened gaussian shaped microwave field
      real(dblprec), dimension(3,Natom), intent(in) :: coord !< Coordinates of the atoms in the system
      character(len=1), intent(in) :: mwf_gauss_spatial      !< Add the frequency broadened gaussian shaped microwave field (Y/P/N)

      ! .. Local variables
      real(dblprec), dimension(3) :: Gauss_factor_R
      real(dblprec), dimension(3) :: R_center
      real(dblprec) :: xmax,ymax, zmax
      integer :: centers, j

      gauss_spatial_site_mwffield=0.0d0

      xmax=maxval(coord(1,:))
      ymax=maxval(coord(2,:))
      zmax=maxval(coord(3,:))

      mwftime=delta_t*(time-1)

      if (mwf_gauss_spatial=='Y') then
         do centers=1, Num_mwf_gauss_centers
            R_center(1:3)=coord(1:3,mwf_gauss_site_center(centers)) ! Coordinates of the centers of the distribution of the spatial magnetic field
            do j=1, Natom
               ! The following definitions allows us to have a gaussian function with different spreading in the three dimensions
               ! The definitions take into account the fact that one might have a 2D system
               Gauss_factor_R(1)=((coord(1,j)-R_center(1))**2)/(2*xmax*mwf_gauss_spatial_space_sigma(1)**2)
               if (abs(xmax)<dbl_tolerance) Gauss_factor_R(1)=0.0d0
               Gauss_factor_R(2)=((coord(2,j)-R_center(2))**2)/(2*ymax*mwf_gauss_spatial_space_sigma(2)**2)
               if (abs(ymax)<dbl_tolerance) Gauss_factor_R(2)=0.0d0
               Gauss_factor_R(3)=((coord(3,j)-R_center(3))**2)/(2*zmax*mwf_gauss_spatial_space_sigma(3)**2)
               if (abs(zmax)<dbl_tolerance) Gauss_factor_R(3)=0.0d0

               ! The gaussian magnetic field equation is A*sin(wt)*exp(-(t**2)/2*sigma_t**2)*exp(-[(Rx-x)**2/2*sigma_x**2]-[(Ry-y)**2/2*sigma_y**2]-[(Rz-z)**2/2*sigma_z**2])
               gauss_spatial_site_mwffield(1,j)=gauss_spatial_site_mwffield(1,j)+mwf_gauss_spatial_ampl*gauss_spatial_mwffield_sites(1,centers)*&
                  sin(mwf_gauss_spatial_freq*mwftime*2*pi)*exp((-(mwftime)**2)/(2*(maxtime*delta_t*mwf_gauss_spatial_time_sigma)**2))*&
                  exp(-(Gauss_factor_R(1)+Gauss_factor_R(2)+Gauss_factor_R(3)))
               gauss_spatial_site_mwffield(2,j)=gauss_spatial_site_mwffield(2,j)+mwf_gauss_spatial_ampl*gauss_spatial_mwffield_sites(2,centers)*&
                  sin(mwf_gauss_spatial_freq*mwftime*2*pi)*exp((-(mwftime)**2)/(2*(maxtime*delta_t*mwf_gauss_spatial_time_sigma)**2))*&
                  exp(-(Gauss_factor_R(1)+Gauss_factor_R(2)+Gauss_factor_R(3)))
               gauss_spatial_site_mwffield(3,j)=gauss_spatial_site_mwffield(3,j)+mwf_gauss_spatial_ampl*gauss_spatial_mwffield_sites(3,centers)*&
                  sin(mwf_gauss_spatial_freq*mwftime*2*pi)*exp((-(mwftime)**2)/(2*(maxtime*delta_t*mwf_gauss_spatial_time_sigma)**2))*&
                  exp(-(Gauss_factor_R(1)+Gauss_factor_R(2)+Gauss_factor_R(3)))
            enddo
         enddo
      else if (mwf_gauss_spatial=='P') then
         if(time<=mwf_gauss_spatial_pulse_time) then
            do centers=1, Num_mwf_gauss_centers
               R_center(1:3)=coord(1:3,mwf_gauss_site_center(centers)) ! Coordinates of the centers of the distribution of the spatial magnetic field
               do j=1, Natom
                  ! The following definitions allows us to have a gaussian function with different spreading in the three dimensions
                  ! The definitions take into account the fact that one might have a 2D system
                  Gauss_factor_R(1)=((coord(1,j)-R_center(1))**2)/(2*xmax*mwf_gauss_spatial_space_sigma(1)**2)
                  if (abs(xmax)<dbl_tolerance) Gauss_factor_R(1)=0.0d0
                  Gauss_factor_R(2)=((coord(2,j)-R_center(2))**2)/(2*ymax*mwf_gauss_spatial_space_sigma(2)**2)
                  if (abs(ymax)<dbl_tolerance) Gauss_factor_R(2)=0.0d0
                  Gauss_factor_R(3)=((coord(3,j)-R_center(3))**2)/(2*zmax*mwf_gauss_spatial_space_sigma(3)**2)
                  if (abs(zmax)<dbl_tolerance) Gauss_factor_R(3)=0.0d0

                  ! The gaussian magnetic field equation is A*sin(wt)*exp(-(t**2)/2*sigma_t**2)*exp(-[(Rx-x)**2/2*sigma_x**2]-[(Ry-y)**2/2*sigma_y**2]-[(Rz-z)**2/2*sigma_z**2])
                  gauss_spatial_site_mwffield(1,j)=gauss_spatial_site_mwffield(1,j)+mwf_gauss_spatial_ampl*gauss_spatial_mwffield_sites(1,centers)*&
                     sin(mwf_gauss_spatial_freq*mwftime*2*pi)*exp((-mwftime**2)/(2*(maxtime*delta_t*mwf_gauss_spatial_time_sigma)**2))*&
                     exp(-(Gauss_factor_R(1)+Gauss_factor_R(2)+Gauss_factor_R(3)))
                  gauss_spatial_site_mwffield(2,j)=gauss_spatial_site_mwffield(2,j)+mwf_gauss_spatial_ampl*gauss_spatial_mwffield_sites(2,centers)*&
                     sin(mwf_gauss_spatial_freq*mwftime*2*pi)*exp((-mwftime**2)/(2*(maxtime*delta_t*mwf_gauss_spatial_time_sigma)**2))*&
                     exp(-(Gauss_factor_R(1)+Gauss_factor_R(2)+Gauss_factor_R(3)))
                  gauss_spatial_site_mwffield(3,j)=gauss_spatial_site_mwffield(3,j)+mwf_gauss_spatial_ampl*gauss_spatial_mwffield_sites(3,centers)*&
                     sin(mwf_gauss_spatial_freq*mwftime*2*pi)*exp((-mwftime**2)/(2*(maxtime*delta_t*mwf_gauss_spatial_time_sigma)**2))*&
                     exp(-(Gauss_factor_R(1)+Gauss_factor_R(2)+Gauss_factor_R(3)))
               enddo
            enddo
         else
            gauss_spatial_site_mwffield=0.0d0
         endif
      endif

   end subroutine calc_spatial_mwf_site_gauss

   !> Calculate a static gaussian shaped field
   subroutine calc_spatial_gauss(Natom,time,gauss_spatial_sigma,coord,gauss_spatial_ampl,delta_t,do_gauss,gauss_pulse_time)

      implicit none

      integer, intent(in) :: time   !< Current simulation step
      integer, intent(in) :: Natom  !< Number of atoms in the system
      integer, intent(in) :: gauss_pulse_time !< Number of time steps in which the static gaussian shaped pulse is on
      real(dblprec), intent(in) :: delta_t    !< Current time step
      real(dblprec), intent(in) :: gauss_spatial_ampl !< Amplitude of the gaussian shaped static field
      real(dblprec), dimension(3), intent(in) :: gauss_spatial_sigma !< Sigma parameter for the gaussian shaped static field
      real(dblprec), dimension(3,Natom) , intent(in) :: coord        !< Coordinates od the atoms in the system
      character(len=1), intent(in) :: do_gauss                       !< Add the Static gaussian shaped field

      ! .. Local variables
      real(dblprec), dimension(3) :: Gauss_factor_R
      real(dblprec), dimension(3) :: R_center
      real(dblprec) :: xmax, ymax, zmax
      integer :: centers, j

      mwftime=(time-1)*delta_t

      gauss_spatial=0.0d0

      xmax=maxval(coord(1,:))
      ymax=maxval(coord(2,:))
      zmax=maxval(coord(3,:))

      if(do_gauss=='Y') then
         ! It is a contant spatially resolved localized gaussian magnetic field
         do centers=1, Num_gauss_centers
            R_center(1:3)=coord(1:3,gauss_site_center(centers))
            do j=1,Natom
               ! The following definitions allows us to have a gaussian function with different spreading in the three dimensions
               ! The definitions take into account the fact that one might have a 2D system
               Gauss_factor_R(1)=((coord(1,j)-R_center(1))**2)/(2*xmax*gauss_spatial_sigma(1)**2)
               if (abs(xmax)<dbl_tolerance) Gauss_factor_R(1)=0.0d0
               Gauss_factor_R(2)=((coord(2,j)-R_center(2))**2)/(2*ymax*gauss_spatial_sigma(2)**2)
               if (abs(ymax)<dbl_tolerance) Gauss_factor_R(2)=0.0d0
               Gauss_factor_R(3)=((coord(3,j)-R_center(3))**2)/(2*zmax*gauss_spatial_sigma(3)**2)
               if (abs(zmax)<dbl_tolerance) Gauss_factor_R(3)=0.0d0

               gauss_spatial(1,j)=gauss_spatial(1,j)+gauss_sites(1,centers)*gauss_spatial_ampl&
                  *exp(-(Gauss_factor_R(1)+Gauss_factor_R(2)+Gauss_factor_R(3)))
               gauss_spatial(2,j)=gauss_spatial(2,j)+gauss_sites(2,centers)*gauss_spatial_ampl&
                  *exp(-(Gauss_factor_R(1)+Gauss_factor_R(2)+Gauss_factor_R(3)))
               gauss_spatial(3,j)=gauss_spatial(3,j)+gauss_sites(3,centers)*gauss_spatial_ampl&
                  *exp(-(Gauss_factor_R(1)+Gauss_factor_R(2)+Gauss_factor_R(3)))
            enddo
         enddo
      else if (do_gauss=='P') then
         ! If the time is least than the time which it is supposed to be on the field is zero
         if(time<=gauss_pulse_time) then
            do centers=1, Num_gauss_centers
               R_center(1:3)=coord(1:3,gauss_site_center(centers))
               do j=1,Natom
                  ! The following definitions allows us to have a gaussian function with different spreading in the three dimensions
                  ! The definitions take into account the fact that one might have a 2D system
                  Gauss_factor_R(1)=((coord(1,j)-R_center(1))**2)/(2*xmax*gauss_spatial_sigma(1)**2)
                  if (abs(xmax)<dbl_tolerance) Gauss_factor_R(1)=0.0d0
                  Gauss_factor_R(2)=((coord(2,j)-R_center(2))**2)/(2*ymax*gauss_spatial_sigma(2)**2)
                  if (abs(ymax)<dbl_tolerance) Gauss_factor_R(2)=0.0d0
                  Gauss_factor_R(3)=((coord(3,j)-R_center(3))**2)/(2*zmax*gauss_spatial_sigma(3)**2)
                  if (abs(zmax)<dbl_tolerance) Gauss_factor_R(3)=0.0d0

                  gauss_spatial(1,j)=gauss_spatial(1,j)+gauss_sites(1,centers)*gauss_spatial_ampl&
                     *exp(-(Gauss_factor_R(1)+Gauss_factor_R(2)+Gauss_factor_R(3)))
                  gauss_spatial(2,j)=gauss_spatial(2,j)+gauss_sites(2,centers)*gauss_spatial_ampl&
                     *exp(-(Gauss_factor_R(1)+Gauss_factor_R(2)+Gauss_factor_R(3)))
                  gauss_spatial(3,j)=gauss_spatial(3,j)+gauss_sites(3,centers)*gauss_spatial_ampl&
                     *exp(-(Gauss_factor_R(1)+Gauss_factor_R(2)+Gauss_factor_R(3)))
               enddo
            enddo
         else
            gauss_spatial=0.0d0
         endif
      endif

   end subroutine calc_spatial_gauss

   !> Calculate a moving gaussian shaped static field which follows a predetermined trajectory
   subroutine calc_moving_gauss(Natom,coord,mov_gauss_ampl,time,mov_gauss_pulse_time,&
         centering,mov_gauss_space_sigma,mov_gauss)
      !
      !.. Implicit declatrations
      implicit none

      integer, intent(in) :: time       !< Current simulation step
      integer, intent(in) :: Natom      !< Number of atoms in the system
      integer, intent(in) :: centering  !< Current position of the center of the field
      integer, intent(in) :: mov_gauss_pulse_time  !< Number of time steps in which the moving static gaussian shaped pulse is on
      real(dblprec), intent(in) :: mov_gauss_ampl  !< Amplitude of microwave field
      real(dblprec), dimension(3), intent(in) :: mov_gauss_space_sigma !< Sigma parameter for the moving gaussian shaped static field
      real(dblprec), dimension(3,Natom), intent(in) :: coord           !< Coordinates of the atoms in the system
      character(len=1), intent(in) :: mov_gauss                        !< Add the moving static gaussian shaped field (Y/P/N)

      !.. Local variables
      real(dblprec) :: xmax,ymax,zmax
      real(dblprec), dimension(3) :: R_center
      real(dblprec), dimension(3) :: Gauss_factor_R

      integer :: j

      xmax=maxval(coord(1,:))
      ymax=maxval(coord(2,:))
      zmax=maxval(coord(3,:))

      ! This is a "constant field"  that is is always on
      if (mov_gauss=='Y') then
         R_center(1:3)=coord(1:3,mov_site_center(centering))
         do j=1, Natom
            Gauss_factor_R(1)=((coord(1,j)-R_center(1))**2)/(2*xmax*mov_gauss_space_sigma(1)**2)
            if (abs(xmax)<dbl_tolerance) Gauss_factor_R(1)=0.0d0
            Gauss_factor_R(2)=((coord(2,j)-R_center(2))**2)/(2*ymax*mov_gauss_space_sigma(2)**2)
            if (abs(ymax)<dbl_tolerance) Gauss_factor_R(2)=0.0d0
            Gauss_factor_R(3)=((coord(3,j)-R_center(3))**2)/(2*zmax*mov_gauss_space_sigma(3)**2)
            if (abs(zmax)<dbl_tolerance) Gauss_factor_R(3)=0.0d0

            mov_gauss_spatial(1,j)=mov_gauss_sites(1,centering)*mov_gauss_ampl&
               *exp(-(Gauss_factor_R(1)+Gauss_factor_R(2)+Gauss_factor_R(3)))
            mov_gauss_spatial(2,j)=mov_gauss_sites(2,centering)*mov_gauss_ampl&
               *exp(-(Gauss_factor_R(1)+Gauss_factor_R(2)+Gauss_factor_R(3)))
            mov_gauss_spatial(3,j)=mov_gauss_sites(3,centering)*mov_gauss_ampl&
               *exp(-(Gauss_factor_R(1)+Gauss_factor_R(2)+Gauss_factor_R(3)))

         enddo
         ! This is a "pulsed field" that is after a certain time it is set to zero
      else if (mov_gauss=='P') then
         if(time<=mov_gauss_pulse_time) then
            R_center(1:3)=coord(1:3,mov_site_center(centering))
            do j=1, Natom
               Gauss_factor_R(1)=((coord(1,j)-R_center(1))**2)/(2*xmax*mov_gauss_space_sigma(1)**2)
               if (abs(xmax)<dbl_tolerance) Gauss_factor_R(1)=0.0d0
               Gauss_factor_R(2)=((coord(2,j)-R_center(2))**2)/(2*ymax*mov_gauss_space_sigma(2)**2)
               if (abs(ymax)<dbl_tolerance) Gauss_factor_R(2)=0.0d0
               Gauss_factor_R(3)=((coord(3,j)-R_center(3))**2)/(2*zmax*mov_gauss_space_sigma(3)**2)
               if (abs(zmax)<dbl_tolerance) Gauss_factor_R(3)=0.0d0

               mov_gauss_spatial(1,j)=mov_gauss_sites(1,centering)*mov_gauss_ampl&
                  *exp(-(Gauss_factor_R(1)+Gauss_factor_R(2)+Gauss_factor_R(3)))
               mov_gauss_spatial(2,j)=mov_gauss_sites(2,centering)*mov_gauss_ampl&
                  *exp(-(Gauss_factor_R(1)+Gauss_factor_R(2)+Gauss_factor_R(3)))
               mov_gauss_spatial(3,j)=mov_gauss_sites(3,centering)*mov_gauss_ampl&
                  *exp(-(Gauss_factor_R(1)+Gauss_factor_R(2)+Gauss_factor_R(3)))
            enddo
         else
            mov_gauss_spatial=0.0d0
         endif
      end if

   end subroutine calc_moving_gauss

   !> Calculate moving gaussian shaped gaussian frequency broadened microwave field through a predetermined trajectory
   subroutine calc_moving_gauss_mwf(Natom,coord,mwf_mov_gauss_ampl,mwf_mov_gauss_freq,mwf_mov_gauss_time_sigma,&
         time,maxtime,delta_t,mwf_mov_gauss_pulse_time,mwf_centering,mwf_mov_gauss_space_sigma,mwf_mov_gauss)

      !.. Implicit declatrations
      implicit none

      integer, intent(in) :: time           !< Current simulation step
      integer, intent(in) :: Natom          !< Number of atoms in the system
      integer, intent(in) :: maxtime        !< Total number of simulation steps
      integer, intent(in) :: mwf_centering  !< Current position of the center of the field
      integer, intent(in) :: mwf_mov_gauss_pulse_time  !< Number of time steps in which the moving microwave gaussian shaped pulse is on
      real(dblprec), intent(in) :: delta_t             !< Current time step
      real(dblprec), intent(in) :: mwf_mov_gauss_ampl       !< Amplitude of the moving microwave gaussian shaped pulse
      real(dblprec), intent(in) :: mwf_mov_gauss_freq       !< Frequency of the moving microwave gaussian shaped pulse
      real(dblprec), intent(in) :: mwf_mov_gauss_time_sigma !< Sigma parameter for frequency gaussian in the moving gaussian shaped microwave field
      real(dblprec), dimension(3), intent(in) :: mwf_mov_gauss_space_sigma !< Sigma parameter for the moving gaussian shaped microwave field
      real(dblprec), dimension(3,Natom), intent(in) :: coord !< Coordinates of the atoms in the system
      character(len=1), intent(in) :: mwf_mov_gauss          !< Add the moving microwave gaussian shaped field (Y/P/N)

      !.. Local variables
      real(dblprec) :: xmax,ymax,zmax
      real(dblprec) :: mwftime
      real(dblprec), dimension(3) :: R_center
      real(dblprec), dimension(3) :: Gauss_factor_R

      integer :: j

      xmax=maxval(coord(1,:))
      ymax=maxval(coord(2,:))
      zmax=maxval(coord(3,:))

      mwftime=(time-1)*delta_t

      if (mwf_mov_gauss=='Y') then
         R_center(1:3)=coord(1:3,mwf_mov_site_center(mwf_centering))
         do j=1, Natom
            Gauss_factor_R(1)=((coord(1,j)-R_center(1))**2)/(2*xmax*mwf_mov_gauss_space_sigma(1)**2)
            if (abs(xmax)<dbl_tolerance) Gauss_factor_R(1)=0.0d0
            Gauss_factor_R(2)=((coord(2,j)-R_center(2))**2)/(2*ymax*mwf_mov_gauss_space_sigma(2)**2)
            if (abs(ymax)<dbl_tolerance) Gauss_factor_R(2)=0.0d0
            Gauss_factor_R(3)=((coord(3,j)-R_center(3))**2)/(2*zmax*mwf_mov_gauss_space_sigma(3)**2)
            if (abs(zmax)<dbl_tolerance)Gauss_factor_R(3)=0.0d0

            mwf_mov_gauss_spatial(1,j)=mov_gauss_mwffield_sites(1,mwf_centering)*mwf_mov_gauss_ampl*sin(mwf_mov_gauss_freq*mwftime*2*pi)*&
               exp(-(Gauss_factor_R(1)+Gauss_factor_R(2)+Gauss_factor_R(3)))*&
               exp((-mwftime**2)/(2*(maxtime*delta_t*mwf_mov_gauss_time_sigma)**2))
            mwf_mov_gauss_spatial(2,j)=mov_gauss_mwffield_sites(2,mwf_centering)*mwf_mov_gauss_ampl*sin(mwf_mov_gauss_freq*mwftime*2*pi)*&
               exp(-(Gauss_factor_R(1)+Gauss_factor_R(2)+Gauss_factor_R(3)))*&
               exp((-mwftime**2)/(2*(maxtime*delta_t*mwf_mov_gauss_time_sigma)**2))
            mwf_mov_gauss_spatial(3,j)=mov_gauss_mwffield_sites(3,mwf_centering)*mwf_mov_gauss_ampl*sin(mwf_mov_gauss_freq*mwftime*2*pi)*&
               exp(-(Gauss_factor_R(1)+Gauss_factor_R(2)+Gauss_factor_R(3)))*&
               exp((-mwftime**2)/(2*(maxtime*delta_t*mwf_mov_gauss_time_sigma)**2))

         enddo
      else if (mwf_mov_gauss=='P') then
         if(time<=mwf_mov_gauss_pulse_time) then
            R_center(1:3)=coord(1:3,mwf_mov_site_center(mwf_centering))
            do j=1, Natom
               Gauss_factor_R(1)=((coord(1,j)-R_center(1))**2)/(2*xmax*mwf_mov_gauss_space_sigma(1)**2)
               if (abs(xmax)<dbl_tolerance) Gauss_factor_R(1)=0.0d0
               Gauss_factor_R(2)=((coord(2,j)-R_center(2))**2)/(2*ymax*mwf_mov_gauss_space_sigma(2)**2)
               if (abs(ymax)<dbl_tolerance) Gauss_factor_R(2)=0.0d0
               Gauss_factor_R(3)=((coord(3,j)-R_center(3))**2)/(2*zmax*mwf_mov_gauss_space_sigma(3)**2)
               if (abs(zmax)<dbl_tolerance)Gauss_factor_R(3)=0.0d0

               mwf_mov_gauss_spatial(1,j)=mov_gauss_mwffield_sites(1,mwf_centering)*mwf_mov_gauss_ampl*sin(mwf_mov_gauss_freq*mwftime*2*pi)*&
                  exp(-(Gauss_factor_R(1)+Gauss_factor_R(2)+Gauss_factor_R(3)))*&
                  exp((-mwftime**2)/(2*(maxtime*delta_t*mwf_mov_gauss_time_sigma)**2))
               mwf_mov_gauss_spatial(2,j)=mov_gauss_mwffield_sites(2,mwf_centering)*mwf_mov_gauss_ampl*sin(mwf_mov_gauss_freq*mwftime*2*pi)*&
                  exp(-(Gauss_factor_R(1)+Gauss_factor_R(2)+Gauss_factor_R(3)))*&
                  exp((-mwftime**2)/(2*(maxtime*delta_t*mwf_mov_gauss_time_sigma)**2))
               mwf_mov_gauss_spatial(3,j)=mov_gauss_mwffield_sites(3,mwf_centering)*mwf_mov_gauss_ampl*sin(mwf_mov_gauss_freq*mwftime*2*pi)*&
                  exp(-(Gauss_factor_R(1)+Gauss_factor_R(2)+Gauss_factor_R(3)))*&
                  exp((-mwftime**2)/(2*(maxtime*delta_t*mwf_mov_gauss_time_sigma)**2))
            enddo
         else
            mwf_mov_gauss_spatial=0.0d0
         endif
      end if


   end subroutine calc_moving_gauss_mwf

   !> Calculate a spherical shaped static field which follows a predetermined trajectory
   subroutine calc_moving_circle(Natom,coord,mov_circle_ampl,time,mov_circle_pulse_time,&
         centering,mov_circle_radius,mov_circle)
      !
      !.. Implicit declatrations
      implicit none

      integer, intent(in) :: time       !< Current simulation step
      integer, intent(in) :: Natom      !< Number of atoms in the system
      integer, intent(in) :: centering  !< Current position of the center of the field
      integer, intent(in) :: mov_circle_pulse_time   !< Number of time steps in which the moving static gaussian shaped pulse is on
      real(dblprec), intent(in) :: mov_circle_ampl   !< Amplitude of microwave field
      real(dblprec), intent(in) :: mov_circle_radius !< Radius for the spherical field
      real(dblprec), dimension(3,Natom), intent(in) :: coord !< Coordinates of the atoms in the system
      character(len=1), intent(in) :: mov_circle             !< Add the moving static circular shaped field (Y/P/N)

      !.. Local variables
      real(dblprec) :: R_fac2
      real(dblprec), dimension(3) :: R_center

      integer :: j

      ! This is a "constant field"  that is is always on
      if (mov_circle=='Y') then
         R_center(1:3)=coord(1:3,mov_site_center(centering))
         do j=1, Natom
            R_fac2=(coord(1,j)-R_center(1))**2+(coord(2,j)-R_center(2))**2+(coord(3,j)-R_center(3))**2
            if (R_fac2<=mov_circle_radius**2) then
               mov_circle_spatial(1,j)=mov_circle_sites(1,centering)*mov_circle_ampl
               mov_circle_spatial(2,j)=mov_circle_sites(2,centering)*mov_circle_ampl
               mov_circle_spatial(3,j)=mov_circle_sites(3,centering)*mov_circle_ampl
            else
               mov_circle_spatial(1,j)=0.0D0
               mov_circle_spatial(2,j)=0.0D0
               mov_circle_spatial(3,j)=0.0D0
            endif
         enddo
         ! This is a "pulsed field" that is after a certain time it is set to zero
      else if (mov_circle=='P') then
         if(time<=mov_circle_pulse_time) then
            R_center(1:3)=coord(1:3,mov_site_center(centering))
            do j=1, Natom
               R_fac2=(coord(1,j)-R_center(1))**2+(coord(2,j)-R_center(2))**2+(coord(3,j)-R_center(3))**2
               if (R_fac2<=mov_circle_radius**2) then
                  mov_circle_spatial(1,j)=mov_circle_sites(1,centering)*mov_circle_ampl
                  mov_circle_spatial(2,j)=mov_circle_sites(2,centering)*mov_circle_ampl
                  mov_circle_spatial(3,j)=mov_circle_sites(3,centering)*mov_circle_ampl
               else
                  mov_circle_spatial(1,j)=0.0D0
                  mov_circle_spatial(2,j)=0.0D0
                  mov_circle_spatial(3,j)=0.0D0
               endif
            enddo
         else
            mov_circle_spatial=0.0d0
         endif
      end if

   end subroutine calc_moving_circle

   !> Calculate a spherical shaped microwavefield field which follows a predetermined trajectory
   subroutine calc_moving_circle_mwf(Natom,coord,mwf_mov_circle_ampl,mwf_mov_circle_freq,mwf_mov_circle_time_sigma,&
         time,maxtime,delta_t,mwf_mov_circle_pulse_time,mwf_mov_circle_radius,mwf_mov_circle)
      !
      !.. Implicit declatrations
      implicit none

      integer, intent(in) :: time           !< Current simulation step
      integer, intent(in) :: Natom          !< Number of atoms in the system
      integer, intent(in) :: maxtime        !< Total number of simulation steps
      integer, intent(in) :: mwf_mov_circle_pulse_time !< Number of time steps in which the moving microwave circular shaped pulse is on
      real(dblprec), intent(in) :: delta_t             !< Current time step
      real(dblprec), intent(in) :: mwf_mov_circle_ampl       !< Amplitude of the moving microwave circular shaped pulse
      real(dblprec), intent(in) :: mwf_mov_circle_freq       !< Frequency of the moving microwave circular shaped pulse
      real(dblprec), intent(in) :: mwf_mov_circle_time_sigma !< Sigma parameter for frequency circular shaped moving microwave field
      real(dblprec), intent(in) :: mwf_mov_circle_radius     !< Radius for the moving circular shaped microwave field
      real(dblprec), dimension(3,Natom), intent(in) :: coord !< Coordinates of the atoms in the system
      character(len=1), intent(in) :: mwf_mov_circle

      !.. Local variables
      real(dblprec) :: R_fac2
      real(dblprec), dimension(3) :: R_center

      integer :: j

      ! This is a "constant field"  that is is always on
      if (mwf_mov_circle=='Y') then
         R_center(1:3)=coord(1:3,mov_site_center(centering))
         do j=1, Natom
            R_fac2=(coord(1,j)-R_center(1))**2+(coord(2,j)-R_center(2))**2+(coord(3,j)-R_center(3))**2
            if (R_fac2<=mwf_mov_circle_radius**2) then
               mwf_mov_circle_spatial(1,j)=mov_circle_mwffield_sites(1,centering)*mwf_mov_circle_ampl*&
                  sin(mwf_mov_circle_freq*mwftime*2*pi)*exp((-mwftime**2)/(2*(maxtime*delta_t*mwf_mov_circle_time_sigma)**2))
               mwf_mov_circle_spatial(2,j)=mov_circle_mwffield_sites(2,centering)*mwf_mov_circle_ampl*&
                  sin(mwf_mov_circle_freq*mwftime*2*pi)*exp((-mwftime**2)/(2*(maxtime*delta_t*mwf_mov_circle_time_sigma)**2))
               mwf_mov_circle_spatial(3,j)=mov_circle_mwffield_sites(3,centering)*mwf_mov_circle_ampl*&
                  sin(mwf_mov_circle_freq*mwftime*2*pi)*exp((-mwftime**2)/(2*(maxtime*delta_t*mwf_mov_circle_time_sigma)**2))
            else
               mwf_mov_circle_spatial(1,j)=0.0D0
               mwf_mov_circle_spatial(2,j)=0.0D0
               mwf_mov_circle_spatial(3,j)=0.0D0
            endif
         enddo
         ! This is a "pulsed field" that is after a certain time it is set to zero
      else if (mwf_mov_circle=='P') then
         if(time<=mwf_mov_circle_pulse_time) then
            R_center(1:3)=coord(1:3,mov_site_center(centering))
            do j=1, Natom
               R_fac2=(coord(1,j)-R_center(1))**2+(coord(2,j)-R_center(2))**2+(coord(3,j)-R_center(3))**2
               if (R_fac2<=mwf_mov_circle_radius**2) then
                  mwf_mov_circle_spatial(1,j)=mov_circle_mwffield_sites(1,centering)*mwf_mov_circle_ampl*&
                     sin(mwf_mov_circle_freq*mwftime*2*pi)*exp((-mwftime**2)/(2*(maxtime*delta_t*mwf_mov_circle_time_sigma)**2))
                  mwf_mov_circle_spatial(2,j)=mov_circle_mwffield_sites(2,centering)*mwf_mov_circle_ampl*&
                     sin(mwf_mov_circle_freq*mwftime*2*pi)*exp((-mwftime**2)/(2*(maxtime*delta_t*mwf_mov_circle_time_sigma)**2))
                  mwf_mov_circle_spatial(3,j)=mov_circle_mwffield_sites(3,centering)*mwf_mov_circle_ampl*&
                     sin(mwf_mov_circle_freq*mwftime*2*pi)*exp((-mwftime**2)/(2*(maxtime*delta_t*mwf_mov_circle_time_sigma)**2))
               else
                  mwf_mov_circle_spatial(1,j)=0.0D0
                  mwf_mov_circle_spatial(2,j)=0.0D0
                  mwf_mov_circle_spatial(3,j)=0.0D0
               endif
            enddo
         else
            mwf_mov_circle_spatial=0.0d0
         endif
      end if

   end subroutine calc_moving_circle_mwf

   !> Calculate a cubic shaped static field which follows a predetermined trajectory
   subroutine calc_moving_square(Natom,coord,mov_square_ampl,time,mov_square_pulse_time,&
         centering,mov_square_dimensions,mov_square)
      !
      !.. Implicit declatrations
      implicit none

      integer, intent(in) :: time       !< Current simulation step
      integer, intent(in) :: Natom      !< Number of atoms in the system
      integer, intent(in) :: centering  !< Current position of the center of the field
      integer, intent(in) :: mov_square_pulse_time   !< Number of time steps in which the moving static cubic shaped pulse is on
      real(dblprec), intent(in) :: mov_square_ampl   !< Amplitude of microwave field
      real(dblprec), dimension(3), intent(in) :: mov_square_dimensions !< Dimensions for the cubic field
      real(dblprec), dimension(3,Natom), intent(in) :: coord !< Coordinates of the atoms in the system
      character(len=1), intent(in) :: mov_square             !< Add the moving static cubic shaped field (Y/P/N)

      !.. Local variables
      real(dblprec) :: xtemp,ytemp,ztemp
      real(dblprec), dimension(3) :: R_center

      integer :: j

      ! This is a "constant field"  that is is always on
      if (mov_square=='Y') then
         R_center(1:3)=coord(1:3,mov_site_center(centering))
         do j=1, Natom
            xtemp=abs(coord(1,j)-R_center(1))
            ytemp=abs(coord(2,j)-R_center(2))
            ztemp=abs(coord(3,j)-R_center(3))
            if ((xtemp<=mov_square_dimensions(1)).and.(ytemp<=mov_square_dimensions(2)).and.(ztemp<=mov_square_dimensions(3)))then
               mov_square_spatial(1,j)=mov_square_sites(1,centering)*mov_square_ampl
               mov_square_spatial(2,j)=mov_square_sites(2,centering)*mov_square_ampl
               mov_square_spatial(3,j)=mov_square_sites(3,centering)*mov_square_ampl
            else
               mov_square_spatial(1,j)=0.0D0
               mov_square_spatial(2,j)=0.0D0
               mov_square_spatial(3,j)=0.0D0
            endif
         enddo
         ! This is a "pulsed field" that is after a certain time it is set to zero
      else if (mov_square=='P') then
         if(time<=mov_square_pulse_time) then
            R_center(1:3)=coord(1:3,mov_site_center(centering))
            do j=1, Natom
               xtemp=abs(coord(1,j)-R_center(1))
               ytemp=abs(coord(2,j)-R_center(2))
               ztemp=abs(coord(3,j)-R_center(3))
               if ( (xtemp<=mov_square_dimensions(1)).and.(ytemp<=mov_square_dimensions(2)).and.(ztemp<=mov_square_dimensions(3)))then
                  mov_square_spatial(1,j)=mov_square_sites(1,centering)*mov_square_ampl
                  mov_square_spatial(2,j)=mov_square_sites(2,centering)*mov_square_ampl
                  mov_square_spatial(3,j)=mov_square_sites(3,centering)*mov_square_ampl
               else
                  mov_square_spatial(1,j)=0.0D0
                  mov_square_spatial(2,j)=0.0D0
                  mov_square_spatial(3,j)=0.0D0
               endif
            enddo
         else
            mov_square_spatial=0.0d0
         endif
      end if

   end subroutine calc_moving_square

   !> Calculate a cubic shaped microwavefield field which follows a predetermined trajectory
   subroutine calc_moving_square_mwf(Natom,coord,mwf_mov_square_ampl,mwf_mov_square_freq,mwf_mov_square_time_sigma,&
         time,maxtime,delta_t,mwf_mov_square_pulse_time,mwf_mov_square_dimensions,mwf_mov_square)
      !
      !.. Implicit declatrations
      implicit none

      integer, intent(in) :: time           !< Current simulation step
      integer, intent(in) :: Natom          !< Number of atoms in the system
      integer, intent(in) :: maxtime        !< Total number of simulation steps
      integer, intent(in) :: mwf_mov_square_pulse_time !< Number of time steps in which the moving microwave circular shaped pulse is on
      real(dblprec), intent(in) :: delta_t             !< Current time step
      real(dblprec), intent(in) :: mwf_mov_square_ampl       !< Amplitude of the moving microwave circular shaped pulse
      real(dblprec), intent(in) :: mwf_mov_square_freq       !< Frequency of the moving microwave circular shaped pulse
      real(dblprec), intent(in) :: mwf_mov_square_time_sigma !< Sigma parameter for frequency circular shaped moving microwave field
      real(dblprec), dimension(3), intent(in) :: mwf_mov_square_dimensions !< Dimensions for the moving circular shaped microwave field
      real(dblprec), dimension(3,Natom), intent(in) :: coord !< Coordinates of the atoms in the system
      character(len=1), intent(in) :: mwf_mov_square

      !.. Local variables
      real(dblprec) :: xtemp,ytemp,ztemp
      real(dblprec), dimension(3) :: R_center

      integer :: j

      ! This is a "constant field"  that is is always on
      if (mwf_mov_square=='Y') then
         R_center(1:3)=coord(1:3,mov_site_center(centering))
         do j=1, Natom
            xtemp=abs(coord(1,j)-R_center(1))
            ytemp=abs(coord(2,j)-R_center(2))
            ztemp=abs(coord(3,j)-R_center(3))
            if ( (xtemp<=mwf_mov_square_dimensions(1)).and.(ytemp<=mwf_mov_square_dimensions(2)).and.(ztemp<=mwf_mov_square_dimensions(3)))then
               mwf_mov_square_spatial(1,j)=mov_square_mwffield_sites(1,centering)*mwf_mov_square_ampl*&
                  sin(mwf_mov_square_freq*mwftime*2*pi)*exp((-mwftime**2)/(2*(maxtime*delta_t*mwf_mov_square_time_sigma)**2))
               mwf_mov_square_spatial(2,j)=mov_square_mwffield_sites(2,centering)*mwf_mov_square_ampl*&
                  sin(mwf_mov_square_freq*mwftime*2*pi)*exp((-mwftime**2)/(2*(maxtime*delta_t*mwf_mov_square_time_sigma)**2))
               mwf_mov_square_spatial(3,j)=mov_square_mwffield_sites(3,centering)*mwf_mov_square_ampl*&
                  sin(mwf_mov_square_freq*mwftime*2*pi)*exp((-mwftime**2)/(2*(maxtime*delta_t*mwf_mov_square_time_sigma)**2))
            else
               mwf_mov_square_spatial(1,j)=0.0D0
               mwf_mov_square_spatial(2,j)=0.0D0
               mwf_mov_square_spatial(3,j)=0.0D0
            endif
         enddo
         ! This is a "pulsed field" that is after a certain time it is set to zero
      else if (mwf_mov_square=='P') then
         if(time<=mwf_mov_square_pulse_time) then
            R_center(1:3)=coord(1:3,mov_site_center(centering))
            do j=1, Natom
               xtemp=abs(coord(1,j)-R_center(1))
               ytemp=abs(coord(2,j)-R_center(2))
               ztemp=abs(coord(3,j)-R_center(3))
               if ( (xtemp<=mwf_mov_square_dimensions(1)).and.(ytemp<=mwf_mov_square_dimensions(2)).and.(ztemp<=mwf_mov_square_dimensions(3)))then
                  mwf_mov_square_spatial(1,j)=mov_square_mwffield_sites(1,centering)*mwf_mov_square_ampl*&
                     sin(mwf_mov_square_freq*mwftime*2*pi)*exp((-mwftime**2)/(2*(maxtime*delta_t*mwf_mov_square_time_sigma)**2))
                  mwf_mov_square_spatial(2,j)=mov_square_mwffield_sites(2,centering)*mwf_mov_square_ampl*&
                     sin(mwf_mov_square_freq*mwftime*2*pi)*exp((-mwftime**2)/(2*(maxtime*delta_t*mwf_mov_square_time_sigma)**2))
                  mwf_mov_square_spatial(3,j)=mov_square_mwffield_sites(3,centering)*mwf_mov_square_ampl*&
                     sin(mwf_mov_square_freq*mwftime*2*pi)*exp((-mwftime**2)/(2*(maxtime*delta_t*mwf_mov_square_time_sigma)**2))
               else
                  mwf_mov_square_spatial(1,j)=0.0D0
                  mwf_mov_square_spatial(2,j)=0.0D0
                  mwf_mov_square_spatial(3,j)=0.0D0
               endif
            enddo
         else
            mwf_mov_square_spatial=0.0d0
         endif
      end if

   end subroutine calc_moving_square_mwf

   !> Routine to setup all the microwave fields, allocate arrays and read the relevant data
   subroutine setup_mwf_fields(Natom,status_flag)

      implicit none

      integer, intent(in) :: Natom         !< Number of atoms in the system
      integer, intent(in) :: status_flag   !< Flag to see wether the arrays must be allocated or deallocated

      !.. Local variables
      integer :: i_stat,i_all

      ! Flag to allocate the fields
      if (status_flag==1) then

         ! Reading the information for all the possible time dependent fields which are "on"
         call read_local_mwf_fields(Natom,mwf,mwf_gauss,mwf_gauss_spatial,mwf_mov_gauss,mov_gauss,&
            do_gauss,mwf_site_file,mwf_gauss_site_file,mwf_gauss_spatial_site_file,mov_gauss_file,&
            mwf_mov_gauss_file,gauss_site_file,mwf_mov_circle_file,mov_circle_file,mwf_mov_circle,mov_circle,&
            mwf_mov_square,mov_square,mov_square_file,mwf_mov_square_file,site_phase,site_phase_file)

         ! Allocating the site dependent single frequency magnetic field
         if(mwf=='S'.or.mwf=='W') then
            allocate(site_mwffield(3,Natom),stat=i_stat)
            call memocc(i_stat,product(shape(site_mwffield))*kind(site_mwffield),'site_mwffield','setup_mwf_fields')
            site_mwffield=0.0D0
         endif

         ! Allocating the site dependent gaussian magnetic field
         if(mwf_gauss=='S'.or.mwf_gauss=='W') then
            allocate(gauss_site_mwffield(3,Natom),stat=i_stat)
            call memocc(i_stat,product(shape(gauss_site_mwffield))*kind(gauss_site_mwffield),'gauss_site_mwffield','setup_mwf_fields')
            gauss_site_mwffield=0.0D0
         endif

         ! Allocate the site depending static gaussian magnetic field
         if(do_gauss=='Y'.or.do_gauss=='P') then
            allocate(gauss_spatial(3,Natom),stat=i_stat)
            call memocc(i_stat,product(shape(gauss_spatial))*kind(gauss_spatial),'gauss_spatial','setup_mwf_fields')
            gauss_spatial=0.0D0
         endif

         ! Allocate the site depending time dependent gaussian magnetic field
         if(mwf_gauss_spatial=='Y'.or.mwf_gauss_spatial=='P') then
            allocate(gauss_spatial_site_mwffield(3,Natom),stat=i_stat)
            call memocc(i_stat,product(shape(gauss_spatial_site_mwffield))*kind(gauss_spatial_site_mwffield),'gauss_spatial_site_mwffield','setup_mwf_fields')
            gauss_spatial_site_mwffield=0.0D0
         endif

         ! Allocate the static spatially resolved moving gaussian field
         if(mov_gauss=='Y'.or.mov_gauss=='P') then
            allocate(mov_gauss_spatial(3,Natom),stat=i_stat)
            call memocc(i_stat,product(shape(mov_gauss_spatial))*kind(mov_gauss_spatial),'mov_gauss_spatial','setup_mwf_fields')
            mov_gauss_spatial=0.0D0
            centering=1
         endif

         ! Allocate the time dependent spatially resolved  moving gaussian field
         if(mwf_mov_gauss=='Y'.or.mwf_mov_gauss=='P') then
            allocate(mwf_mov_gauss_spatial(3,Natom),stat=i_stat)
            call memocc(i_stat,product(shape(mwf_mov_gauss_spatial))*kind(mwf_mov_gauss_spatial),'mwf_mov_gauss_spatial','setup_mwf_fields')
            mwf_mov_gauss_spatial=0.0D0
            mwf_centering=1
         endif

         ! Allocate the static spatially resolved moving circular field
         if(mov_circle=='Y'.or.mov_circle=='P') then
            allocate(mov_circle_spatial(3,Natom),stat=i_stat)
            call memocc(i_stat,product(shape(mov_circle_spatial))*kind(mov_circle_spatial),'mov_circle_spatial','setup_mwf_fields')
            mov_circle_spatial=0.0D0
            centering=1
         endif

         ! Allocate the time dependent spatially resolved  moving circle field
         if(mwf_mov_circle=='Y'.or.mwf_mov_circle=='P') then
            allocate(mwf_mov_circle_spatial(3,Natom),stat=i_stat)
            call memocc(i_stat,product(shape(mwf_mov_circle_spatial))*kind(mwf_mov_circle_spatial),'mwf_mov_circle_spatial','setup_mwf_fields')
            mwf_mov_circle_spatial=0.0D0
            mwf_centering=1
         endif

         ! Allocate the static spatially resolved moving cubic field
         if(mov_square=='Y'.or.mov_square=='P') then
            allocate(mov_square_spatial(3,Natom),stat=i_stat)
            call memocc(i_stat,product(shape(mov_square_spatial))*kind(mov_square_spatial),'mov_square_spatial','setup_mwf_fields')
            mov_square_spatial=0.0D0
            centering=1
         endif

         ! Allocate the time dependent spatially resolved  moving cubic field
         if(mwf_mov_square=='Y'.or.mwf_mov_square=='P') then
            allocate(mwf_mov_square_spatial(3,Natom),stat=i_stat)
            call memocc(i_stat,product(shape(mwf_mov_square_spatial))*kind(mwf_mov_square_spatial),'mwf_mov_square_spatial','setup_mwf_fields')
            mwf_mov_square_spatial=0.0D0
            mwf_centering=1
         endif

         ! Flag to deallocate the present fields
      else if (status_flag==-1) then

         ! Deallocate the site dependent single frequency microwave field
         if(mwf=='S'.or.mwf=='W') then
            i_all=-product(shape(site_mwffield))*kind(site_mwffield)
            deallocate(site_mwffield,stat=i_stat)
            call memocc(i_stat,i_all,'site_mwffield','setup_mwf_fields')

            i_all=-product(shape(mwffield_sites))*kind(mwffield_sites)
            deallocate(mwffield_sites,stat=i_stat)
            call memocc(i_stat,i_all,'mwffield_sites','setup_mwf_fields')

            i_all=-product(shape(field_site_phase))*kind(field_site_phase)
            deallocate(field_site_phase,stat=i_stat)
            call memocc(i_stat,i_all,'field_site_phase','setup_mwf_fields')

         endif

         ! Deallocate the site dependent gaussian magnetic field
         if(mwf_gauss=='S'.or.mwf_gauss=='W') then
            i_all=-product(shape(gauss_site_mwffield))*kind(gauss_site_mwffield)
            deallocate(gauss_site_mwffield,stat=i_stat)
            call memocc(i_stat,i_all,'gauss_site_mwffield','setup_mwf_fields')

            i_all=-product(shape(gauss_mwffield_sites))*kind(gauss_mwffield_sites)
            deallocate(gauss_mwffield_sites,stat=i_stat)
            call memocc(i_stat,i_all,'gauss_mwffield_sites','setup_mwf_fields')

         endif

         ! Deallocate the static spatially resolved magnetic field
         if(do_gauss=='Y'.or.do_gauss=='P') then
            i_all=-product(shape(gauss_spatial))*kind(gauss_spatial)
            deallocate(gauss_spatial,stat=i_stat)
            call memocc(i_stat,i_all,'gauss_spatial','setup_mwf_fields')

            i_all=-product(shape(gauss_sites))*kind(gauss_sites)
            deallocate(gauss_sites,stat=i_stat)
            call memocc(i_stat,i_all,'gauss_sites','setup_mwf_fields')

            i_all=-product(shape(gauss_site_center))*kind(gauss_site_center)
            deallocate(gauss_site_center,stat=i_stat)
            call memocc(i_stat,i_all,'gauss_site_center','setup_mwf_fields')

         endif

         ! Deallocate the time dependent spatially resolved magnetic field
         if(mwf_gauss_spatial=='Y'.or.mwf_gauss_spatial=='P') then
            i_all=-product(shape(gauss_spatial_site_mwffield))*kind(gauss_spatial_site_mwffield)
            deallocate(gauss_spatial_site_mwffield,stat=i_stat)
            call memocc(i_stat,i_all,'gauss_spatial_site_mwffield','setup_mwf_fields')

            i_all=-product(shape(gauss_spatial_mwffield_sites))*kind(gauss_spatial_mwffield_sites)
            deallocate(gauss_spatial_mwffield_sites,stat=i_stat)
            call memocc(i_stat,i_all,'gauss_spatial_site_mwffield','setup_mwf_fields')

            i_all=-product(shape(mwf_gauss_site_center))*kind(mwf_gauss_site_center)
            deallocate(mwf_gauss_site_center,stat=i_stat)
            call memocc(i_stat,i_all,'mwf_gauss_site_center','setup_mwf_fields')

         endif

         ! Deallocate the static spatially resolved moving magnetic field
         if(mov_gauss=='Y'.or.mov_gauss=='P') then
            i_all=-product(shape(mov_gauss_spatial))*kind(mov_gauss_spatial)
            deallocate(mov_gauss_spatial,stat=i_stat)
            call memocc(i_stat,i_all,'mov_gauss_spatial','setup_mwf_fields')

            i_all=-product(shape(mov_gauss_sites))*kind(mov_gauss_sites)
            deallocate(mov_gauss_sites,stat=i_stat)
            call memocc(i_stat,i_all,'mov_gauss_sites','setup_mwf_fields')

            i_all=-product(shape(mov_site_center))*kind(mov_site_center)
            deallocate(mov_site_center,stat=i_stat)
            call memocc(i_stat,i_all,'mov_site_center','setup_mwf_fields')

         endif

         ! Deallocate the static spatially resolved moving magnetic field
         if(mwf_mov_gauss=='Y'.or.mwf_mov_gauss=='P') then
            i_all=-product(shape(mwf_mov_gauss_spatial))*kind(mwf_mov_gauss_spatial)
            deallocate(mwf_mov_gauss_spatial,stat=i_stat)
            call memocc(i_stat,i_all,'mwf_mov_gauss_spatial','setup_mwf_fields')

            i_all=-product(shape(mov_gauss_mwffield_sites))*kind(mov_gauss_mwffield_sites)
            deallocate(mov_gauss_mwffield_sites,stat=i_stat)
            call memocc(i_stat,i_all,'mov_gauss_mwffield_sites','setup_mwf_fields')

            i_all=-product(shape(mwf_mov_site_center))*kind(mwf_mov_site_center)
            deallocate(mwf_mov_site_center,stat=i_stat)
            call memocc(i_stat,i_all,'mwf_mov_site_center','setup_mwf_fields')

         endif

         ! Deallocate the static circular spatially resolved moving magnetic field
         if(mov_circle=='Y'.or.mov_circle=='P') then
            i_all=-product(shape(mov_circle_spatial))*kind(mov_circle_spatial)
            deallocate(mov_circle_spatial,stat=i_stat)
            call memocc(i_stat,i_all,'mov_circle_spatial','setup_mwf_fields')

            i_all=-product(shape(mov_circle_sites))*kind(mov_circle_sites)
            deallocate(mov_circle_sites,stat=i_stat)
            call memocc(i_stat,i_all,'mov_circle_sites','setup_mwf_fields')

            i_all=-product(shape(mov_site_center))*kind(mov_site_center)
            deallocate(mov_site_center,stat=i_stat)
            call memocc(i_stat,i_all,'mov_site_center','setup_mwf_fields')

         endif

         ! Deallocate the static spatially resolved moving magnetic field
         if(mwf_mov_circle=='Y'.or.mwf_mov_circle=='P') then
            i_all=-product(shape(mwf_mov_circle_spatial))*kind(mwf_mov_circle_spatial)
            deallocate(mwf_mov_circle_spatial,stat=i_stat)
            call memocc(i_stat,i_all,'mov_circle_spatial','setup_mwf_fields')

            i_all=-product(shape(mov_circle_mwffield_sites))*kind(mov_circle_mwffield_sites)
            deallocate(mov_circle_mwffield_sites,stat=i_stat)
            call memocc(i_stat,i_all,'mov_circle_mwffield_sites','setup_mwf_fields')

            i_all=-product(shape(mwf_mov_site_center))*kind(mwf_mov_site_center)
            deallocate(mwf_mov_site_center,stat=i_stat)
            call memocc(i_stat,i_all,'mwf_mov_site_center','setup_mwf_fields')

         endif

         ! Deallocate the static cubic spatially resolved moving magnetic field
         if(mov_square=='Y'.or.mov_square=='P') then
            i_all=-product(shape(mov_square_spatial))*kind(mov_square_spatial)
            deallocate(mov_square_spatial,stat=i_stat)
            call memocc(i_stat,i_all,'mov_square_spatial','setup_mwf_fields')

            i_all=-product(shape(mov_square_sites))*kind(mov_square_sites)
            deallocate(mov_square_sites,stat=i_stat)
            call memocc(i_stat,i_all,'mov_square_sites','setup_mwf_fields')

            i_all=-product(shape(mov_site_center))*kind(mov_site_center)
            deallocate(mov_site_center,stat=i_stat)
            call memocc(i_stat,i_all,'mov_site_center','setup_mwf_fields')

         endif

         ! Deallocate the static spatially resolved moving magnetic field
         if(mwf_mov_square=='Y'.or.mwf_mov_square=='P') then
            i_all=-product(shape(mwf_mov_square_spatial))*kind(mwf_mov_square_spatial)
            deallocate(mwf_mov_square_spatial,stat=i_stat)
            call memocc(i_stat,i_all,'mov_square_spatial','setup_mwf_fields')

            i_all=-product(shape(mov_square_mwffield_sites))*kind(mov_square_mwffield_sites)
            deallocate(mov_square_mwffield_sites,stat=i_stat)
            call memocc(i_stat,i_all,'mov_square_mwffield_sites','setup_mwf_fields')

            i_all=-product(shape(mwf_mov_site_center))*kind(mwf_mov_site_center)
            deallocate(mwf_mov_site_center,stat=i_stat)
            call memocc(i_stat,i_all,'mwf_mov_site_center','setup_mwf_fields')

         endif

      end if

   end subroutine setup_mwf_fields

   !> Subroutine to read the site dependent microwave fields
   subroutine read_local_mwf_fields(Natom,mwf,mwf_gauss,mwf_gauss_spatial,mwf_mov_gauss,mov_gauss,&
         do_gauss,mwf_site_file,mwf_gauss_site_file,mwf_gauss_spatial_site_file,mov_gauss_file,&
         mwf_mov_gauss_file,gauss_site_file,mwf_mov_circle_file,mov_circle_file,mwf_mov_circle,&
         mov_circle,mwf_mov_square,mov_square,mov_square_file,mwf_mov_square_file,site_phase, &
         site_phase_file)

      implicit none

      integer, intent(in) :: Natom !< Number of atoms in the system
      character(len=1), intent(in) :: mwf                !< Add monochromatic microwave field (Y/P/S/W/N)
      character(len=1), intent(in) :: do_gauss           !< Add the Static gaussian shaped field
      character(len=1), intent(in) :: mwf_gauss          !< Add frequency broadened microwave field (Y/P/S/W/N)
      character(len=1), intent(in) :: mov_gauss          !< Add the moving static gaussian shaped field (Y/P/N)
      character(len=1), intent(in) :: site_phase         !< Add the site dependent phase (Y/N)
      character(len=1), intent(in) :: mov_circle         !< Add the moving static circular shaped field (Y/P/N)
      character(len=1), intent(in) :: mov_square         !< Add the moving static cubic shaped field (Y/P/N)
      character(len=1), intent(in) :: mwf_mov_gauss      !< Add the moving microwave gaussian shaped field (Y/P/N)
      character(len=1), intent(in) :: mwf_mov_circle     !< Add the moving microwave circular shaped field (Y/P/N)
      character(len=1), intent(in) :: mwf_mov_square     !< Add the moving microwave cubic shaped field (Y/P/N)
      character(len=1), intent(in) :: mwf_gauss_spatial  !< Add the frequency broadened gaussian shaped microwave field (Y/P/N)
      character(len=35), intent(in) :: mwf_site_file    !< File name for the site dependent monochromatic micorwave field
      character(len=35), intent(in) :: mov_gauss_file   !< Moving static gaussian shaped field trajectory file
      character(len=35), intent(in) :: site_phase_file  !< File for site dependent phase for the microwavefield
      character(len=35), intent(in) :: mov_circle_file  !< Moving static circular shaped field trajectory file
      character(len=35), intent(in) :: mov_square_file  !< Moving static cubic shaped field trajectory file
      character(len=35), intent(in) :: gauss_site_file  !< File name for the site dependent static gaussian shaped field
      character(len=35), intent(in) :: mwf_mov_gauss_file !< Moving frequency broadened gaussian shaped microwave filed trajectory file
      character(len=35), intent(in) :: mwf_mov_circle_file !< Moving frequency broadened circular shaped microwave filed trajectory file
      character(len=35), intent(in) :: mwf_mov_square_file !< Moving frequency broadened cubic shaped microwave filed trajectory file
      character(len=35), intent(in) :: mwf_gauss_site_file !< File name for the site dependent frequency broadened microwave field
      character(len=35), intent(in) :: mwf_gauss_spatial_site_file !< File name for the site dependent gaussian shaped gaussian broadened microwave field


      !.. Local variables
      integer :: i, flines, isite, i_stat

      ! Read the file for the monochormatic site dependent microwave field
      if (mwf=='S'.or.mwf=='W') then

         open  (ofileno, file=trim(mwf_site_file))

         flines=0
         ! Pre-read file to get number of lines
         do
            read(222,*,end=200)  isite
            flines=flines+1
         end do
         200 continue
         rewind(222)

         ! If the size of the file is NATOM then there is no problem
         if (flines.ne.Natom) write(*,'(2x,a)') 'WARNING: Size of the MWFFIELD_SITES is not NATOM'
         ! Allocate the site-dependent field
         allocate(mwffield_sites(3,flines),stat=i_stat)
         call memocc(i_stat,product(shape(mwffield_sites))*kind(mwffield_sites),'mwffield_sites','read_local_mwf_fields')
         mwffield_sites=0.0D0

         allocate(field_site_phase(flines),stat=i_stat)
         call memocc(i_stat,product(shape(field_site_phase))*kind(field_site_phase),'field_site_phase','read_local_mwf_fields')
         field_site_phase=0.0D0

         do i=1, flines
            read(222,*) isite, mwffield_sites(1,isite), mwffield_sites(2,isite), mwffield_sites(3,isite)
         end do

         close (ofileno)



         if (site_phase=='Y') then
            open  (ofileno,file=trim(site_phase_file))
            flines=0
            do
               read(458,*,end=250) isite
               flines=flines+1
            end do
            250 continue
            rewind(458)
            if (flines.ne.Natom) write(*,'(2x,a)') 'WARNING: Size of the MWFFIELD_SITES is not NATOM'

            do i=1, flines
               read(458,*) isite, field_site_phase(isite)
            end do
         endif
      endif

      ! Read the file for the frequency broadened site dependent microwave field
      if (mwf_gauss=='S'.or.mwf_gauss=='W') then

         open  (ofileno, file=trim(mwf_gauss_site_file))

         flines=0
         ! Pre-read file to get number of lines
         do
            read(223,*,end=400)  isite
            flines=flines+1
         end do
         400 continue
         rewind(223)

         ! If the size of the file is NATOM then there is no problem
         if (flines.ne.Natom) write(*,'(2x,a)') 'WARNING: Size of the GAUSS_MWFFIELD_SITES is not NATOM'
         ! Allocate the site-dependent field
         allocate(gauss_mwffield_sites(3,flines),stat=i_stat)
         call memocc(i_stat,product(shape(gauss_mwffield_sites))*kind(gauss_mwffield_sites),'gauss_mwffield_sites','read_local_mwf_fields')
         gauss_mwffield_sites=0.0D0

         do i=1, flines
            read(223,*) isite, gauss_mwffield_sites(1,isite), gauss_mwffield_sites(2,isite), gauss_mwffield_sites(3,isite)
         end do

         close (ofileno)

      endif

      ! Read the file for the frequency broadened, gaussian shaped site dependent microwave field
      if (mwf_gauss_spatial=='Y'.or.mwf_gauss_spatial=='P') then

         open  (ofileno, file=trim(mwf_gauss_spatial_site_file))

         flines=0
         ! Pre-read file to get number of lines
         do
            read(224,*,end=600)  isite
            flines=flines+1
         end do
         600 continue
         rewind(224)

         ! Allocate the site-dependent field
         allocate(gauss_spatial_mwffield_sites(3,flines),stat=i_stat)
         call memocc(i_stat,product(shape(gauss_spatial_mwffield_sites))*kind(gauss_spatial_mwffield_sites),'gauss_spatial_mwffield_sites','read_local_mwf_fields')
         gauss_spatial_mwffield_sites=0.0D0

         allocate(mwf_gauss_site_center(flines),stat=i_stat)
         call memocc(i_stat,product(shape(mwf_gauss_site_center))*kind(mwf_gauss_site_center),'mwf_gauss_site_center','read_local_mwf_fields')
         mwf_gauss_site_center=0

         Num_mwf_gauss_centers=flines
         do i=1, flines
            read(224,*) mwf_gauss_site_center(i), gauss_spatial_mwffield_sites(1,i), gauss_spatial_mwffield_sites(2,i), gauss_spatial_mwffield_sites(3,i)
         end do

         close (ofileno)

      endif

      ! Read the file for the moving gaussian shape microwave field
      if (mwf_mov_gauss=='Y'.or.mwf_mov_gauss=='P') then
         open  (ofileno, file=trim(mwf_mov_gauss_file))

         flines=0
         ! Pre-read file to get number of lines
         do
            read(225,*,end=800)  isite
            flines=flines+1
         end do
         800 continue
         rewind(225)

         ! Allocate the site-dependent field
         allocate(mov_gauss_mwffield_sites(3,flines),stat=i_stat)
         call memocc(i_stat,product(shape(mov_gauss_mwffield_sites))*kind(mov_gauss_mwffield_sites),'mov_gauss_mwffield_sites','read_local_mwf_fields')
         mov_gauss_mwffield_sites=0.0D0

         allocate(mwf_mov_site_center(flines),stat=i_stat)
         call memocc(i_stat,product(shape(mwf_mov_site_center))*kind(mwf_mov_site_center),'mwf_mov_site_center','read_local_mwf_fields')
         mwf_mov_site_center=0
         ! If the size of the field is not NATOM there will be a problem
         Num_mwf_mov_centers=flines
         do i=1, flines
            read(225,*) mwf_mov_site_center(i), mov_gauss_mwffield_sites(1,i), mov_gauss_mwffield_sites(2,i), mov_gauss_mwffield_sites(3,i)
         end do

         close (ofileno)

      endif

      ! Read the file for the moving circular shape microwave field
      if (mwf_mov_circle=='Y'.or.mwf_mov_circle=='P') then
         open  (ofileno, file=trim(mwf_mov_circle_file))

         flines=0
         ! Pre-read file to get number of lines
         do
            read(265,*,end=880)  isite
            flines=flines+1
         end do
         880 continue
         rewind(265)

         ! Allocate the site-dependent field
         allocate(mov_circle_mwffield_sites(3,flines),stat=i_stat)
         call memocc(i_stat,product(shape(mov_circle_mwffield_sites))*kind(mov_circle_mwffield_sites),'mov_circle_mwffield_sites','read_local_mwf_fields')
         mov_circle_mwffield_sites=0.0D0

         allocate(mwf_mov_site_center(flines),stat=i_stat)
         call memocc(i_stat,product(shape(mwf_mov_site_center))*kind(mwf_mov_site_center),'mwf_mov_site_center','read_local_mwf_fields')
         mwf_mov_site_center=0
         ! If the size of the field is not NATOM there will be a problem
         Num_mwf_mov_centers=flines
         do i=1, flines
            read(265,*) mwf_mov_site_center(i), mov_circle_mwffield_sites(1,i), mov_circle_mwffield_sites(2,i), mov_circle_mwffield_sites(3,i)
         end do

         close (ofileno)

      endif

      ! Read the file for the moving cubic shape microwave field
      if (mwf_mov_square=='Y'.or.mwf_mov_square=='P') then
         open  (ofileno, file=trim(mwf_mov_square_file))

         flines=0
         ! Pre-read file to get number of lines
         do
            read(265,*,end=890)  isite
            flines=flines+1
         end do
         890 continue
         rewind(265)

         ! Allocate the site-dependent field
         allocate(mov_square_mwffield_sites(3,flines),stat=i_stat)
         call memocc(i_stat,product(shape(mov_square_mwffield_sites))*kind(mov_square_mwffield_sites),'mov_square_mwffield_sites','read_local_mwf_fields')
         mov_square_mwffield_sites=0.0D0

         allocate(mwf_mov_site_center(flines),stat=i_stat)
         call memocc(i_stat,product(shape(mwf_mov_site_center))*kind(mwf_mov_site_center),'mwf_mov_site_center','read_local_mwf_fields')
         mwf_mov_site_center=0
         ! If the size of the field is not NATOM there will be a problem
         Num_mwf_mov_centers=flines
         do i=1, flines
            read(265,*) mwf_mov_site_center(i), mov_square_mwffield_sites(1,i), mov_square_mwffield_sites(2,i), mov_square_mwffield_sites(3,i)
         end do

         close (ofileno)

      endif

      ! Read the file for the moving gaussian shape static field
      if (mov_gauss=='Y'.or.mov_gauss=='P') then
         open  (ofileno, file=trim(mov_gauss_file))

         flines=0
         ! Pre-read file to get number of lines
         do
            read(226,*,end=500)  isite
            flines=flines+1
         end do
         500 continue
         rewind(226)

         ! Allocate the site-dependent field
         allocate(mov_gauss_sites(3,flines),stat=i_stat)
         call memocc(i_stat,product(shape(mov_gauss_sites))*kind(mov_gauss_sites),'mov_gauss_sites','read_local_mwf_fields')
         mov_gauss_sites=0.0D0

         allocate(mov_site_center(flines),stat=i_stat)
         call memocc(i_stat,product(shape(mov_site_center))*kind(mov_site_center),'mov_site_center','read_local_mwf_fields')
         mov_site_center=0
         ! If the size of the field is not NATOM there will be a problem
         Num_mov_centers=flines
         do i=1, flines
            read(226,*) mov_site_center(i), mov_gauss_sites(1,i), mov_gauss_sites(2,i), mov_gauss_sites(3,i)
         end do

         close (ofileno)

      endif

      ! Read the file for the moving circular shape static field
      if (mov_circle=='Y'.or.mov_circle=='P') then
         open  (ofileno, file=trim(mov_circle_file))
         flines=0
         ! Pre-read file to get number of lines
         do
            read(266,*,end=700)  isite
            flines=flines+1
         end do
         700 continue
         rewind(266)

         ! Allocate the site-dependent field
         allocate(mov_circle_sites(3,flines),stat=i_stat)
         call memocc(i_stat,product(shape(mov_circle_sites))*kind(mov_circle_sites),'mov_circle_sites','read_local_mwf_fields')
         mov_circle_sites=0.0D0

         allocate(mov_site_center(flines),stat=i_stat)
         call memocc(i_stat,product(shape(mov_site_center))*kind(mov_site_center),'mov_site_center','read_local_mwf_fields')
         mov_site_center=0
         ! If the size of the field is not NATOM there will be a problem
         Num_mov_centers=flines
         do i=1, flines
            read(266,*) mov_site_center(i), mov_circle_sites(1,i), mov_circle_sites(2,i), mov_circle_sites(3,i)
         end do

         close (ofileno)

      endif

      ! Read the file for the moving cubic shape static field
      if (mov_square=='Y'.or.mov_square=='P') then
         open  (ofileno, file=trim(mov_square_file))
         flines=0
         ! Pre-read file to get number of lines
         do
            read(266,*,end=750)  isite
            flines=flines+1
         end do
         750 continue
         rewind(266)

         ! Allocate the site-dependent field
         allocate(mov_square_sites(3,flines),stat=i_stat)
         call memocc(i_stat,product(shape(mov_square_sites))*kind(mov_square_sites),'mov_square_sites','read_local_mwf_fields')
         mov_square_sites=0.0D0

         allocate(mov_site_center(flines),stat=i_stat)
         call memocc(i_stat,product(shape(mov_site_center))*kind(mov_site_center),'mov_site_center','read_local_mwf_fields')
         mov_site_center=0
         ! If the size of the field is not NATOM there will be a problem
         Num_mov_centers=flines
         do i=1, flines
            read(266,*) mov_site_center(i), mov_square_sites(1,i), mov_square_sites(2,i), mov_square_sites(3,i)
         end do

         close (ofileno)

      endif

      ! Read the file for the static gaussian shaped field
      if (do_gauss=='Y'.or.do_gauss=='P') then
         open  (ofileno, file=trim(gauss_site_file))

         flines=0
         ! Pre-read file to get number of lines
         do
            read(227,*,end=100)  isite
            flines=flines+1
         end do
         100 continue
         rewind(227)

         ! Allocate the site-dependent field
         allocate(gauss_sites(3,flines),stat=i_stat)
         call memocc(i_stat,product(shape(gauss_sites))*kind(gauss_sites),'gauss_sites','read_local_mwf_fields')
         gauss_sites=0.0D0

         allocate(gauss_site_center(flines),stat=i_stat)
         call memocc(i_stat,product(shape(gauss_site_center))*kind(gauss_site_center),'gauss_site_center','read_local_mwf_fields')
         gauss_site_center=0
         ! If the size of the field is not NATOM there will be a problem
         Num_gauss_centers=flines
         do i=1, flines
            read(227,*) gauss_site_center(i), gauss_sites(1,i), gauss_sites(2,i), gauss_sites(3,i)
         end do

         close (ofileno)

      endif

   end subroutine read_local_mwf_fields

   !---------------------------------------------------------------------------
   !> @brief
   !> Read input parameters.
   !
   !> @author
   !> Anders Bergman
   !---------------------------------------------------------------------------
   subroutine read_parameters_mwf(ifile)
      use FileParser

      implicit none

      ! ... Formal Arguments ...
      integer, intent(in) :: ifile   !< File to read from
      !
      ! ... Local Variables ...
      character(len=50) :: keyword, cache
      integer :: rd_len, i_err, i_errb
      logical :: comment

      do
         10     continue
         ! Read file character for character until first whitespace
         keyword=""
         call bytereader(keyword,rd_len,ifile,i_errb)

         ! converting Capital letters
         call caps2small(keyword)

         ! check for comment markers (currently % and #)
         comment=(scan(trim(keyword),'%')==1).or.(scan(trim(keyword),'#')==1).or.&
            (scan(trim(keyword),'*')==1).or.(scan(trim(keyword),'=')==1.or.&
            (scan(trim(keyword),'!')==1))

         if (comment) then
            read(ifile,*)
         else
            ! Parse keyword
            keyword=trim(keyword)
            select case(keyword)
         case('mwf')  ! Flag for the monochormatic microwave field
            read(ifile,*,iostat=i_err) mwf
            if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err

         case('mwfdir')  ! Direction of the monochormatic microwave field
            read(ifile,*,iostat=i_err) mwfdir
            if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err

         case('mwfampl')  ! Amplitude of the monochormatic microwave field
            read(ifile,*,iostat=i_err) mwfampl
            if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err

         case('mwffreq') ! Frequency of the monochormatic microwave field
            read(ifile,*,iostat=i_err) mwffreq
            if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err

         case('mwf_pulse_time') ! Number of time steps in which the monochormatic microwave field is on
            read(ifile,*,iostat=i_err) mwf_pulse_time
            if(i_err/=0) write(*,*) 'ERROR:  Reading ', trim(keyword),' data',i_err

         case('mwf_site_file') ! Site dependent monochromatic microwave field
            read(ifile,'(a)',iostat=i_err) cache
            if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err
            mwf_site_file=adjustl(trim(cache))

         case('site_phase') ! Flag for site dependent phase
            read(ifile,*,iostat=i_err) site_phase
            if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err

         case('site_phase_file') ! Site dependent phase for the monochromatic microwave field
            read(ifile,'(a)',iostat=i_err) cache
            if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err
            site_phase_file=adjustl(trim(cache))

         case('prn_mwf') ! Flag for printing the monochromatic microwave field
            read(ifile,*,iostat=i_err) prn_mwf
            if(i_err/=0) write(*,*) 'ERROR:  Reading ', trim(keyword),' data',i_err

         case('mwf_step') ! Interval for sampling the monochromatic microwave field
            read(ifile,*,iostat=i_err) mwf_step
            if(i_err/=0) write(*,*) 'ERROR:  Reading ', trim(keyword),' data',i_err

         case('mwf_buff') ! Buffer size for the monochromatic microwave field
            read(ifile,*,iostat=i_err) mwf_buff
            if(i_err/=0) write(*,*) 'ERROR:  Reading ', trim(keyword),' data',i_err

         case('mwf_gauss') ! Flag for the frequency broadened microwave field
            read(ifile,*,iostat=i_err) mwf_gauss
            if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err

         case('mwf_gauss_dir') ! Direction of the frequency broadened microwave field
            read(ifile,*,iostat=i_err) mwf_gauss_dir
            if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err

         case('mwf_gauss_ampl') ! Amplitude of the frequency broadened microwave field
            read(ifile,*,iostat=i_err) mwf_gauss_ampl
            if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err

         case('mwf_gauss_freq') ! Frequency of the frequency broadened microwave field
            read(ifile,*,iostat=i_err) mwf_gauss_freq
            if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err

         case('mwf_gauss_time_sigma') ! Frequency sigma parameter for the frequency broadened microwave field
            read(ifile,*,iostat=i_err) mwf_gauss_time_sigma
            if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword), ' data',i_err

         case('mwf_gauss_pulse_time') ! Number of time steps in which the frequency broadened  microwave field is on
            read(ifile,*,iostat=i_err) mwf_gauss_pulse_time
            if(i_err/=0) write(*,*) 'ERROR:  Reading ', trim(keyword),' data',i_err

         case('mwf_gauss_site_file') ! Site dependent frequency broadended microwave field
            read(ifile,'(a)',iostat=i_err) cache
            if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err
            mwf_gauss_site_file=adjustl(trim(cache))

         case('prn_mwf_gauss') ! Flag for printing the static gaussian shaped field
            read(ifile,*,iostat=i_err) prn_mwf_gauss
            if(i_err/=0) write(*,*) 'ERROR:  Reading ', trim(keyword),' data',i_err

         case('mwf_gauss_step') ! Interval for sampling the frequency broadened microwave field
            read(ifile,*,iostat=i_err) mwf_gauss_step
            if(i_err/=0) write(*,*) 'ERROR:  Reading ', trim(keyword),' data',i_err

         case('mwf_gauss_buff') ! Buffer size for the frequency broadened microwave field
            read(ifile,*,iostat=i_err) mwf_gauss_buff
            if(i_err/=0) write(*,*) 'ERROR:  Reading ', trim(keyword),' data',i_err

         case('mwf_gauss_spatial') ! Flag for the frequency broadened gaussian shaped microwave field
            read(ifile,*,iostat=i_err) mwf_gauss_spatial
            if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err

         case('mwf_gauss_spatial_ampl') ! Amplitude for the frequency broadened gaussian shaped microwave field
            read(ifile,*,iostat=i_err) mwf_gauss_spatial_ampl
            if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err

         case('mwf_gauss_spatial_freq') ! Frequency for the frequency broadened gaussian shaped microwave field
            read(ifile,*,iostat=i_err) mwf_gauss_spatial_freq
            if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err

         case('mwf_gauss_spatial_pulse_time') ! Number of time steps in which the frequency broadened gaussian shaped microwave field is on
            read(ifile,*,iostat=i_err) mwf_gauss_spatial_pulse_time
            if(i_err/=0) write(*,*) 'ERROR:  Reading ', trim(keyword),' data',i_err

         case('mwf_gauss_spatial_time_sigma') ! Frequency sigma parameter for the frequency broadened gaussian shaped microwave field
            read(ifile,*,iostat=i_err) mwf_gauss_spatial_time_sigma
            if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword), ' data',i_err

         case('mwf_gauss_spatial_space_sigma') ! Spatial sigma parameter for the frequency broadened gaussian shaped microwave field
            read(ifile,*,iostat=i_err) mwf_gauss_spatial_space_sigma(1), mwf_gauss_spatial_space_sigma(2), mwf_gauss_spatial_space_sigma(3)
            if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword), ' data',i_err

         case('mwf_gauss_spatial_site_file') ! Site dependent gaussian shaped gaussian broadened microwave field
            read(ifile,'(a)',iostat=i_err) cache
            if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err
            mwf_gauss_spatial_site_file=adjustl(trim(cache))

         case('prn_mwf_gauss_spatial') ! Flag for printing the gaussian shaped frequency broadened microwave field
            read(ifile,*,iostat=i_err) prn_mwf_gauss_spatial
            if(i_err/=0) write(*,*) 'ERROR:  Reading ', trim(keyword),' data',i_err

         case('mwf_gauss_spatial_step') ! Interval for sampling the gaussian shaped frequency broadened microwave field
            read(ifile,*,iostat=i_err) mwf_gauss_spatial_step
            if(i_err/=0) write(*,*) 'ERROR:  Reading ', trim(keyword),' data',i_err

         case('mwf_gauss_spatial_buff') ! Buffer size for the gaussian shaped frequency broadened microwave field
            read(ifile,*,iostat=i_err) mwf_gauss_spatial_buff
            if(i_err/=0) write(*,*) 'ERROR:  Reading ', trim(keyword),' data',i_err

         case('mov_gauss') ! Flag for the moving static gaussian shaped pulse
            read(ifile,*,iostat=i_err) mov_gauss
            if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err

         case('mov_gauss_step') ! Number of time steps in which the positions of the moving static gaussian shaped pulse is updated
            read(ifile,*,iostat=i_err) mov_gauss_step
            if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err

         case('mov_gauss_ampl') ! Amplitude of the moving static gaussian shaped pulse
            read(ifile,*,iostat=i_err) mov_gauss_ampl
            if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err

         case('mov_gauss_pulse_time') ! Number of time steps in which the moving static gaussian shaped pulse is on
            read(ifile,*,iostat=i_err) mov_gauss_pulse_time
            if(i_err/=0) write(*,*) 'ERROR:  Reading ', trim(keyword),' data',i_err

         case('mov_gauss_space_sigma') ! Sigma parameter for the moving gaussian shaped static field
            read(ifile,*,iostat=i_err) mov_gauss_space_sigma(1), mov_gauss_space_sigma(2), mov_gauss_space_sigma(3)
            if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword), ' data',i_err

         case('mov_gauss_file') ! Moving static gaussian shaped field trajectory file
            read(ifile,'(a)',iostat=i_err) cache
            if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err
            mov_gauss_file=adjustl(trim(cache))

         case('prn_mov_gauss') ! Flag for printing the moving gaussian shaped static field
            read(ifile,*,iostat=i_err) prn_mov_gauss
            if(i_err/=0) write(*,*) 'ERROR:  Reading ', trim(keyword),' data',i_err

         case('mov_gauss_pstep') !  Interval for sampling the moving gaussian shaped static field
            read(ifile,*,iostat=i_err) mov_gauss_pstep
            if(i_err/=0) write(*,*) 'ERROR:  Reading ', trim(keyword),' data',i_err

         case('mov_gauss_buff') ! Buffer size for the moving gaussian shaped static field
            read(ifile,*,iostat=i_err) mov_gauss_buff
            if(i_err/=0) write(*,*) 'ERROR:  Reading ', trim(keyword),' data',i_err

         case('mov_circle') ! Flag for the moving static circular shaped pulse
            read(ifile,*,iostat=i_err) mov_circle
            if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err

         case('mov_circle_step') ! Number of time steps in which the positions of the moving static circular shaped pulse is updated
            read(ifile,*,iostat=i_err) mov_circle_step
            if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err

         case('mov_circle_ampl') ! Amplitude of the moving static circular shaped pulse
            read(ifile,*,iostat=i_err) mov_circle_ampl
            if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err

         case('mov_circle_pulse_time') ! Number of time steps in which the moving static circular shaped pulse is on
            read(ifile,*,iostat=i_err) mov_circle_pulse_time
            if(i_err/=0) write(*,*) 'ERROR:  Reading ', trim(keyword),' data',i_err

         case('mov_circle_radius') ! Radius for the moving circular shaped static field
            read(ifile,*,iostat=i_err) mov_circle_radius
            if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword), ' data',i_err

         case('mov_circle_file') ! Moving static circular shaped field trajectory file
            read(ifile,'(a)',iostat=i_err) cache
            if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err
            mov_circle_file=adjustl(trim(cache))

         case('prn_mov_circle') ! Flag for printing the moving gaussian shaped static field
            read(ifile,*,iostat=i_err) prn_mov_circle
            if(i_err/=0) write(*,*) 'ERROR:  Reading ', trim(keyword),' data',i_err

         case('mov_circle_pstep') !  Interval for sampling the moving circular shaped static field
            read(ifile,*,iostat=i_err) mov_circle_pstep
            if(i_err/=0) write(*,*) 'ERROR:  Reading ', trim(keyword),' data',i_err

         case('mov_circle_buff') ! Buffer size for the moving circular shaped static field
            read(ifile,*,iostat=i_err) mov_circle_buff
            if(i_err/=0) write(*,*) 'ERROR:  Reading ', trim(keyword),' data',i_err

         case('mov_square') ! Flag for the moving static cubic shaped pulse
            read(ifile,*,iostat=i_err) mov_square
            if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err

         case('mov_square_step') ! Number of time steps in which the positions of the moving static cubic shaped pulse is updated
            read(ifile,*,iostat=i_err) mov_square_step
            if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err

         case('mov_square_ampl') ! Amplitude of the moving static cubic shaped pulse
            read(ifile,*,iostat=i_err) mov_square_ampl
            if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err

         case('mov_square_pulse_time') ! Number of time steps in which the moving static cubic shaped pulse is on
            read(ifile,*,iostat=i_err) mov_square_pulse_time
            if(i_err/=0) write(*,*) 'ERROR:  Reading ', trim(keyword),' data',i_err

         case('mov_square_dimensions') ! Dimensions for the moving cubic shaped static field
            read(ifile,*,iostat=i_err) mov_square_dimensions(1), mov_square_dimensions(2), mov_square_dimensions(3)
            if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword), ' data',i_err

         case('mov_square_file') ! Moving static cubic shaped field trajectory file
            read(ifile,'(a)',iostat=i_err) cache
            if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err
            mov_square_file=adjustl(trim(cache))

         case('prn_mov_square') ! Flag for printing the moving cubic shaped static field
            read(ifile,*,iostat=i_err) prn_mov_square
            if(i_err/=0) write(*,*) 'ERROR:  Reading ', trim(keyword),' data',i_err

         case('mov_square_pstep') !  Interval for sampling the moving cubic shaped static field
            read(ifile,*,iostat=i_err) mov_square_pstep
            if(i_err/=0) write(*,*) 'ERROR:  Reading ', trim(keyword),' data',i_err

         case('mov_square_buff') ! Buffer size for the moving cubic shaped static field
            read(ifile,*,iostat=i_err) mov_square_buff
            if(i_err/=0) write(*,*) 'ERROR:  Reading ', trim(keyword),' data',i_err

         case('mwf_mov_gauss') ! Flag for the moving microwave gaussian shaped field
            read(ifile,*,iostat=i_err) mwf_mov_gauss
            if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err

         case('mwf_mov_gauss_step') ! Number of time steps in which the positions of the moving microwave gaussian shaped pulse is updated
            read(ifile,*,iostat=i_err) mwf_mov_gauss_step
            if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err

         case('mwf_mov_gauss_ampl') ! Amplitude of the moving microwave gaussian shaped pulse
            read(ifile,*,iostat=i_err) mwf_mov_gauss_ampl
            if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err

         case('mwf_mov_gauss_freq') ! Frequency of the moving microwave gaussian shaped pulse
            read(ifile,*,iostat=i_err) mwf_mov_gauss_freq
            if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err

         case('mwf_mov_gauss_pulse_time') ! Number of time steps in which the moving microwave gaussian shaped pulse is on
            read(ifile,*,iostat=i_err) mwf_mov_gauss_pulse_time
            if(i_err/=0) write(*,*) 'ERROR:  Reading ', trim(keyword),' data',i_err

         case('mwf_mov_gauss_space_sigma') ! Sigma parameter for the moving gaussian shaped microwave field
            read(ifile,*,iostat=i_err) mwf_mov_gauss_space_sigma(1), mwf_mov_gauss_space_sigma(2), mwf_mov_gauss_space_sigma(3)
            if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword), ' data',i_err

         case('mwf_mov_gauss_time_sigma') ! Sigma parameter for frequency gaussian in the moving gaussian shaped microwave field
            read(ifile,*,iostat=i_err) mwf_mov_gauss_time_sigma
            if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword), ' data',i_err

         case('mwf_mov_gauss_file') ! Moving frequency broadened gaussian shaped microwave filed trajectory file
            read(ifile,'(a)',iostat=i_err) cache
            if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err
            mwf_mov_circle_file=adjustl(trim(cache))

         case('prn_mwf_mov_gauss') ! Flag for printing the moving gaussian shaped microwave field
            read(ifile,*,iostat=i_err) prn_mwf_mov_gauss
            if(i_err/=0) write(*,*) 'ERROR:  Reading ', trim(keyword),' data',i_err

         case('mwf_mov_gauss_pstep') !  Interval for sampling the moving gaussian shaped microwave field
            read(ifile,*,iostat=i_err) mwf_mov_gauss_pstep
            if(i_err/=0) write(*,*) 'ERROR:  Reading ', trim(keyword),' data',i_err

         case('mwf_mov_gauss_buff') ! Buffer size for the moving gaussian shaped microwave field
            read(ifile,*,iostat=i_err) mwf_mov_gauss_buff
            if(i_err/=0) write(*,*) 'ERROR:  Reading ', trim(keyword),' data',i_err

         case('mwf_mov_circle') ! Flag for the moving microwave circular shaped field
            read(ifile,*,iostat=i_err) mwf_mov_circle
            if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err

         case('mwf_mov_circle_step') ! Number of time steps in which the positions of the moving microwave circle shaped pulse is updated
            read(ifile,*,iostat=i_err) mwf_mov_circle_step
            if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err

         case('mwf_mov_circle_ampl') ! Amplitude of the moving microwave circular shaped pulse
            read(ifile,*,iostat=i_err) mwf_mov_circle_ampl
            if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err

         case('mwf_mov_circle_freq') ! Frequency of the moving microwave circular shaped pulse
            read(ifile,*,iostat=i_err) mwf_mov_circle_freq
            if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err

         case('mwf_mov_circle_pulse_time') ! Number of time steps in which the moving microwave gaussian shaped pulse is on
            read(ifile,*,iostat=i_err) mwf_mov_gauss_pulse_time
            if(i_err/=0) write(*,*) 'ERROR:  Reading ', trim(keyword),' data',i_err

         case('mwf_mov_circle_radius') ! Radius for the moving circle shaped microwave field
            read(ifile,*,iostat=i_err) mwf_mov_circle_radius
            if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword), ' data',i_err

         case('mwf_mov_circle_time_sigma') ! Sigma parameter for frequency circle in the moving gaussian shaped microwave field
            read(ifile,*,iostat=i_err) mwf_mov_circle_time_sigma
            if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword), ' data',i_err

         case('mwf_mov_circle_file') ! Moving frequency broadened circular shaped microwave filed trajectory file
            read(ifile,'(a)',iostat=i_err) cache
            if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err
            mwf_mov_gauss_file=adjustl(trim(cache))

         case('prn_mwf_mov_circle') ! Flag for printing the moving circular shaped microwave field
            read(ifile,*,iostat=i_err) prn_mwf_mov_circle
            if(i_err/=0) write(*,*) 'ERROR:  Reading ', trim(keyword),' data',i_err

         case('mwf_mov_circle_pstep') !  Interval for sampling the moving circular shaped microwave field
            read(ifile,*,iostat=i_err) mwf_mov_circle_pstep
            if(i_err/=0) write(*,*) 'ERROR:  Reading ', trim(keyword),' data',i_err

         case('mwf_mov_circle_buff') ! Buffer size for the moving circular shaped microwave field
            read(ifile,*,iostat=i_err) mwf_mov_circle_buff
            if(i_err/=0) write(*,*) 'ERROR:  Reading ', trim(keyword),' data',i_err

         case('mwf_mov_square') ! Flag for the moving microwave cubic shaped field
            read(ifile,*,iostat=i_err) mwf_mov_square
            if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err

         case('mwf_mov_square_step') ! Number of time steps in which the positions of the moving microwave square shaped pulse is updated
            read(ifile,*,iostat=i_err) mwf_mov_square_step
            if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err

         case('mwf_mov_square_ampl') ! Amplitude of the moving microwave cubic shaped pulse
            read(ifile,*,iostat=i_err) mwf_mov_square_ampl
            if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err

         case('mwf_mov_square_freq') ! Frequency of the moving microwave cubic shaped pulse
            read(ifile,*,iostat=i_err) mwf_mov_square_freq
            if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err

         case('mwf_mov_square_pulse_time') ! Number of time steps in which the moving microwave cubic shaped pulse is on
            read(ifile,*,iostat=i_err) mwf_mov_square_pulse_time
            if(i_err/=0) write(*,*) 'ERROR:  Reading ', trim(keyword),' data',i_err

         case('mwf_mov_square_dimensions') ! Dimensions for the moving cubic shaped microwave field
            read(ifile,*,iostat=i_err) mwf_mov_square_dimensions(1), mwf_mov_square_dimensions(2),mwf_mov_square_dimensions(3)
            if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword), ' data',i_err

         case('mwf_mov_square_time_sigma') ! Sigma parameter for frequency cubic in the moving gaussian shaped microwave field
            read(ifile,*,iostat=i_err) mwf_mov_square_time_sigma
            if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword), ' data',i_err

         case('mwf_mov_square_file') ! Moving frequency broadened cubic shaped microwave filed trajectory file
            read(ifile,'(a)',iostat=i_err) cache
            if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err
            mwf_mov_square_file=adjustl(trim(cache))

         case('prn_mwf_mov_square') ! Flag for printing the moving cubic shaped microwave field
            read(ifile,*,iostat=i_err) prn_mwf_mov_square
            if(i_err/=0) write(*,*) 'ERROR:  Reading ', trim(keyword),' data',i_err

         case('mwf_mov_square_pstep') !  Interval for sampling the moving square shaped microwave field
            read(ifile,*,iostat=i_err) mwf_mov_square_pstep
            if(i_err/=0) write(*,*) 'ERROR:  Reading ', trim(keyword),' data',i_err

         case('mwf_mov_square_buff') ! Buffer size for the moving cubic shaped microwave field
            read(ifile,*,iostat=i_err) mwf_mov_square_buff
            if(i_err/=0) write(*,*) 'ERROR:  Reading ', trim(keyword),' data',i_err

         case('do_gauss') ! Static gaussian shaped field
            read(ifile,*,iostat=i_err) do_gauss
            if(i_err/=0) write(*,*) 'ERROR:  Reading ', trim(keyword),' data',i_err

         case('gauss_spatial_ampl') ! Amplitude of the gaussian shaped static field
            read(ifile,*,iostat=i_err) gauss_spatial_ampl
            if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err

         case('gauss_spatial_sigma') ! Sigma parameter for the gaussian shaped static field
            read(ifile,*,iostat=i_err) gauss_spatial_sigma(1), gauss_spatial_sigma(2), gauss_spatial_sigma(3)
            if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword), ' data',i_err

         case('gauss_pulse_time') ! Number of time steps in which the static gaussian shaped pulse is on
            read(ifile,*,iostat=i_err) gauss_pulse_time
            if(i_err/=0) write(*,*) 'ERROR:  Reading ', trim(keyword),' data',i_err

         case('gauss_site_file') ! Static gaussian shaped field
            read(ifile,'(a)',iostat=i_err) cache
            if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err
            gauss_site_file=adjustl(trim(cache))

         case('prn_gauss') ! Print gaussian shaped static field
            read(ifile,*,iostat=i_err) prn_gauss
            if(i_err/=0) write(*,*) 'ERROR:  Reading ', trim(keyword),' data',i_err

         case('gauss_step') ! Time interval to measure the microwave field
            read(ifile,*,iostat=i_err) gauss_step
            if(i_err/=0) write(*,*) 'ERROR:  Reading ', trim(keyword),' data',i_err

         case('gauss_buff') ! Buffer size for the static gaussian shaped field
            read(ifile,*,iostat=i_err) gauss_buff
            if(i_err/=0) write(*,*) 'ERROR:  Reading ', trim(keyword),' data',i_err

         end select
      end if

      ! End of file
      if (i_errb==20) goto 20
      ! End of row
      if (i_errb==10) goto 10
   end do

   20  continue

   rewind(ifile)
   return
end subroutine read_parameters_mwf

end module MicroWaveField
