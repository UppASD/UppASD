!------------------------------------------------------------------------------------
module Omegas
   use Parameters
   use Profiling
   use Correlation_type
   !
   implicit none
   !

   ! Working variables to perform the printing of the correlation
   !! integer :: nw  !< Number of frequencies to sample
   !! real(dblprec), dimension(:), allocatable :: w                     !< Frequencies
   !
   integer :: sc_window_fun  !< Choice of FFT window function (1=box, 2=Hann, 3=Hamming, 4=Blackman-Harris)
   logical :: print_real_w 
   
   public

contains

   !----------------------------------------------------------------------------------
   ! SUBROUTINE: set_w
   !> @brief Calculate suitable values of frequencies for \f$ \mathbf{S}\left(\mathbf{q},t\right) \rightarrow \mathbf{S}\left(\mathbf{q},\omega\right)\f$ transform
   !> @todo Change setup to smarter algorithm with respect to the exchange strength of the
   !> system
   !---------------------------------------------------------------------------------
   subroutine set_w(delta_t,cc)
      !
      use Constants, only : pi, hbar_mev, hbar, mry
      !
      implicit none
      !
      real(dblprec), intent(in) :: delta_t  !< Time step
      type(corr_t), intent(inout) :: cc
      !
      integer :: i_stat,j
      real(dblprec) :: dt !< Time step
      real(dblprec) :: dww
      real(dblprec) :: emin,emax
      !
      print *,'Debug:', 'henlo'
      ! Check if sc_eres and sc_emax are defined. 
      print '("Sampling intervals before",i6,i6)', cc%sc_nstep,cc%sc_step
      ! If so, set sc_step and sc_nstep accordingly
      if (cc%sc_emax>0.0_dblprec .and. cc%sc_eres>0.0_dblprec) then
         cc%sc_nstep = int(cc%sc_emax/cc%sc_eres)
         cc%sc_step = int ( 0.5_dblprec * (2.0_dblprec * pi * hbar) / (delta_t * (cc%sc_emax * mry)) ) 
      end if
      print '("Sampling intervals after",i6,i6)', cc%sc_nstep,cc%sc_step

      !
      dt=cc%sc_step*delta_t
      cc%nw = cc%sc_nstep
      print *,'dEb = ',cc%sc_nstep
      if (.not. allocated(cc%w)) then
         allocate(cc%w(cc%nw),stat=i_stat)
         call memocc(i_stat,product(shape(cc%w))*kind(cc%w),'cc%w','set_w')
      end if
      dww = 2.0_dblprec*pi/(cc%nw*dt)

      ! Previous convention was j*ddw
      do j=1,cc%nw
         cc%w(j)=(j-1)*dww
      end do

      emin=hbar_mev*(cc%w(2)-cc%w(1))
      emax=hbar_mev*cc%w(cc%nw)
      write(*,'(1x,a,f6.3,a,f6.1,a)') 'Spin wave sampling between ',emin,' meV and ',emax,' meV.'
   end subroutine set_w

   !---------------------------------------------------------------------------------
   ! SUBROUTINE: calc_corr_w
   !> @brief Calculate correction to the sampling frequencies, w
   !> Perform correction to all frequencies >= sc_tidx
   !---------------------------------------------------------------------------------
   subroutine calc_corr_w(deltat_upd,cc) !sc_step,sc_nstep,sc_tidx)

      use Constants, only : pi

      implicit none

      real(dblprec), intent(in) :: deltat_upd  !< Updated time step
      type(corr_t), intent(inout) :: cc
      !integer, intent(in) :: sc_step  !< Number of timesteps between samplings
      !integer, intent(in) :: sc_nstep !< Number of samplings
      !integer, intent(in) :: sc_tidx  !< Current sampling index

      integer :: j
      real(dblprec) :: dt
      real(dblprec) :: dww

      dt = cc%sc_step*deltat_upd ! Time interval between new samples
      cc%nw = cc%sc_nstep
      dww = 2*pi/(cc%nw*dt)

      ! Perform correction of the current/future sampling frequencies
      do j=cc%sc_tidx+1,cc%nw ! The addition of 1 is due to function calling at the end of the previous sampling period
         cc%w(j)=j*dww
      end do

   end subroutine calc_corr_w

!------------------------------------------------------------------------------------
! FUNCTION: sc_window_fac
!> @brief Window function factor for the different types of windowing
!> @author Anders Bergman
!------------------------------------------------------------------------------------
real(dblprec) function sc_window_fac(sc_window_fun,step,nstep)

   use Constants, only : pi

   !
   implicit none
   !
   integer, intent(in)  :: sc_window_fun
   integer, intent(in)  :: step
   integer, intent(in)  :: nstep
   !
   real(dblprec) :: dum
   !
   dum=1.0_dblprec
   select case(sc_window_fun)
      ! Hann
      case(2)
         dum= (0.50_dblprec-0.50_dblprec*cos(2.0_dblprec*pi*(step-1._dblprec)/(nstep-1._dblprec)))
      ! Hamming
      case(3)
         dum= (0.54_dblprec-0.46_dblprec*cos(2.0_dblprec*pi*(step-1._dblprec)/(nstep-1._dblprec)))
      ! Hamming v2
      case(32)
         dum= (0.53836_dblprec- 0.46164_dblprec*cos(2.0_dblprec*pi*(step-1._dblprec)/(nstep-1._dblprec)))
      ! Blackman-Harris
      case(4)
         dum=  &
            (0.35785_dblprec-0.48829_dblprec*cos(2.0_dblprec*pi*(step-1._dblprec)/(nstep-1.0_dblprec))+ &
            0.14128_dblprec*cos(4.0_dblprec*pi*(step-1._dblprec)/(nstep-1.0_dblprec))   &
            -0.01168_dblprec*cos(6.0_dblprec*pi*(step-1._dblprec)/(nstep-1.0_dblprec)))
      ! Nuttal
      case(5)
         dum=  &
            (0.355768_dblprec-0.478396_dblprec*cos(2.0_dblprec*pi*(step-1._dblprec)/(nstep-1.0_dblprec))+ &
            0.144232_dblprec*cos(4.0_dblprec*pi*(step-1._dblprec)/(nstep-1.0_dblprec))   &
            -0.012604_dblprec*cos(6.0_dblprec*pi*(step-1._dblprec)/(nstep-1.0_dblprec)))
      ! Nuttal
      !!!case(5)
      !!!   dum=  &
      !!!      (0.355768_dblprec-0.478396_dblprec*cos(2.0_dblprec*pi*(step)/(nstep))+ &
      !!!      0.144232_dblprec*cos(4.0_dblprec*pi*(step)/(nstep))   &
      !!!      -0.012604_dblprec*cos(6.0_dblprec*pi*(step)/(nstep)))
      !!!   ! Square windows
      case default
         dum=1.0_dblprec
      end select
      !
      sc_window_fac=dum
      return
      !
   end function sc_window_fac

   subroutine set_w_print_option(real_time_measure)
      !
      implicit none
      !
      character(len=1), intent(in)  :: real_time_measure
      !
      print_real_w = (real_time_measure=='Y')
      !
  end subroutine set_w_print_option


end module omegas
