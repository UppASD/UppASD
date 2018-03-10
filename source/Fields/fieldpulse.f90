!-------------------------------------------------------------------------------
! MODULE: FieldPulse
!> @brief Data and routines for calculating time dependent field pulses
!> @copyright
!! Copyright (C) 2008-2018 UppASD group
!! This file is distributed under the terms of the
!! GNU General Public License.
!! See http://www.gnu.org/copyleft/gpl.txt
!-------------------------------------------------------------------------------
module FieldPulse
   use Parameters
   use Profiling
   !
   implicit none
   !
   !Magnetic field pulse
   integer :: bpulse_npar !< Number of parameters used to describe pulse
   integer :: bpulse_step
   real(dblprec) :: bpulse_dt !< Time step for magnetic pulse
   real(dblprec) :: bpulse_time
   real(dblprec), dimension(3) :: bpulse_b0 !< Starting amplitude for magnetic field pulse
   real(dblprec), dimension(10) :: bpulse_par !< Parameters used to describe pulse
   real(dblprec), dimension(3) :: bpulsefield !< Current applied magnetic field from pulse
   real(dblprec) :: ba, bb
   !
   private

   public :: bpulse_time, bpulse_step, bpulsefield
   public :: bpulse, read_bpulse


contains

   !-----------------------------------------------------------------------------
   ! SUBROUTINE: bpulse
   !> @brief Evaluate field pulse at time t=tin
   !-----------------------------------------------------------------------------
   subroutine bpulse(do_bpulse,tin)

      !.. Implicit declarations
      implicit none

      integer, intent(in) :: do_bpulse  !< Add magnetic field pulse (0=no, 1-4 for different shapes)
      real(dblprec),intent(in) :: tin  !< Current time

      !.. Local arrays
      real(dblprec), dimension(3) :: befftemp
      real(dblprec) :: tpulse

      !.. Executable statements

      if (do_bpulse == 1) then
         call exppulse(tin,tpulse)
      else if (do_bpulse == 2) then
         call gaussianpulse(tin,tpulse)
      else if (do_bpulse == 3) then
         call polexppulse(tin,tpulse)
      else if (do_bpulse == 4) then
         call squarepulse(tin,tpulse)
      else
         tpulse=0.0_dblprec
      end if

      befftemp(1) = bpulse_b0(1)*tpulse
      befftemp(2) = bpulse_b0(2)*tpulse
      befftemp(3) = bpulse_b0(3)*tpulse
      bpulsefield(1) = befftemp(1)
      bpulsefield(2) = befftemp(2)
      bpulsefield(3) = befftemp(3)


   end subroutine bpulse

   !-----------------------------------------------------------------------------
   ! SUBROUTINE: squarepulse
   !> @brief A square pulse
   !-----------------------------------------------------------------------------
   subroutine squarepulse (t,tpulse)
      real(dblprec) :: t
      real(dblprec) :: tpulse
      if (t < bpulse_par(2) ) then
         tpulse = 0.0d0
      elseif (t >= bpulse_par(2) .and. t < bpulse_par(3) ) then
         tpulse = bpulse_par(6)
      else
         tpulse = 0.0d0
      end if
   end subroutine squarepulse

   !-----------------------------------------------------------------------------
   ! SUBROUTINE: exppulse
   !> An almost square pulse with exponential head and tail
   !-----------------------------------------------------------------------------
   subroutine exppulse (t,tpulse)
      real(dblprec) :: t
      real(dblprec) :: tpulse
      if (t <= bpulse_par(2) ) then
         tpulse = bpulse_par(6) * exp((t-bpulse_par(2))/ba)
      elseif (t > bpulse_par(2) .and. t < bpulse_par(3) ) then
         tpulse = bpulse_par(6)
      else
         tpulse = bpulse_par(6) * exp((bpulse_par(3)-t)/bb)
      end if
   end subroutine exppulse

   !-----------------------------------------------------------------------------
   ! SUBROUTINE: gaussianpulse
   !> A Gaussian shaped pulse
   !-----------------------------------------------------------------------------
   subroutine gaussianpulse(t,tpulse)
      real(dblprec) :: t
      real(dblprec) :: tpulse
      tpulse = bpulse_par(6)*ba*exp(bb*(t-bpulse_par(2))**2)
   end subroutine gaussianpulse

   !-----------------------------------------------------------------------------
   ! SUBROUTINE: polexppulse
   !> An exponentially decaying pulse
   !-----------------------------------------------------------------------------
   subroutine polexppulse(t,tpulse)
      real(dblprec) :: t
      real(dblprec) :: tpulse
      tpulse = bpulse_par(6)*bb*(t-bpulse_par(1))**bpulse_par(3)*exp(-ba*t)
   end subroutine polexppulse

   !-----------------------------------------------------------------------------
   ! SUBROUTINE: prn_bpulse0
   !> @brief Print header of magnetic pulse field
   !-----------------------------------------------------------------------------
   subroutine prn_bpulse0(simid,mstep)
      !

      !.. Implicit declarations
      implicit none

      character(len=8), intent(in) :: simid !< Name of simulation
      integer, intent(in) :: mstep !< Current simulation step

      character(len=30) :: filn

      !.. Executable statements
      write (filn,'(''bpulse.'',A8,''.out'')') simid
      open(ofileno, file=filn, position="append")
      write (ofileno,10002) mstep, bpulse_time, bpulsefield(1), bpulsefield(2), bpulsefield(3)
      close(ofileno)
      return

      write (*,*) "Error writing output file"

      10002 format (i8,2x,es16.8,2x,es16.8,es16.8,es16.8)

   end subroutine prn_bpulse0

   !-----------------------------------------------------------------------------
   ! SUBROUTINE: prn_bpulse
   !> @brief Print shape of magnetic pulse field
   !-----------------------------------------------------------------------------
   subroutine prn_bpulse(simid, mstep)
      !

      !.. Implicit declarations
      implicit none

      character(len=8), intent(in) :: simid !< Name of simulation
      integer, intent(in) :: mstep !< Current simulation step

      !integer :: mstep !< Current simulation step
      character(len=30) :: filn

      !.. Executable statements
      write (filn,'(''bpulse.'',A8,''.out'')') simid
      open(ofileno, file=filn, position="append")
      write (ofileno,10002) mstep, bpulse_time, bpulsefield(1), bpulsefield(2), bpulsefield(3)
      close(ofileno)
      return

      write (*,*) "Error writing output file"

      10002 format (i8,2x,es16.8,2x,es16.8,es16.8,es16.8)

   end subroutine prn_bpulse

   !-----------------------------------------------------------------------------
   ! SUBROUTINE: read_bpulse
   !> @brief Read in field pulse from file
   !-----------------------------------------------------------------------------
   subroutine read_bpulse(do_bpulse, bpulsefile)

      implicit none

      integer, intent(in) :: do_bpulse  !< Add magnetic field pulse (0=no, 1-4 for different shapes)
      character(LEN=30), intent(in) :: bpulsefile !< File containing data for magnetic field pulse

      integer :: i
      open (ifileno, file=adjustl(bpulsefile))
      read (ifileno,*)
      read (ifileno,*) bpulse_b0(1), bpulse_b0(2), bpulse_b0(3)
      read (ifileno,*) bpulse_dt
      read (ifileno,*) bpulse_step
      read (ifileno,*) bpulse_npar
      do i=1,bpulse_npar
         read (ifileno,*) bpulse_par(i)
      end do

      close(ifileno)
      if (do_bpulse == 1) then
         ba = (bpulse_par(1)-bpulse_par(2))/log(bpulse_par(5)/bpulse_par(6))
         bb = (bpulse_par(3)-bpulse_par(4))/log(bpulse_par(5)/bpulse_par(6))
      else if (do_bpulse == 2) then
         ba = 1
         bb = -1/(2*bpulse_par(3)**2)
      else if (do_bpulse == 3) then
         ba = bpulse_par(3)/(bpulse_par(2)-bpulse_par(1))
         bb = 1 / ( (bpulse_par(2)-bpulse_par(1))**bpulse_par(3) * exp(-ba*bpulse_par(2)) )
      else
         write (*,*) "No other pulse profiles implemented yet!"
      end if

   end subroutine read_bpulse


end module FieldPulse
