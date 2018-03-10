!> @brief
!! Fake calls for CUDA functionalities
!> @copyright
!! Copyright (C) 2008-2018 UppASD group
!! This file is distributed under the terms of the
!! GNU General Public License. 
!! See http://www.gnu.org/copyleft/gpl.txt
module nocuda
   use Parameters
   implicit none
   ! Initiate pointers for C/C++ implementation
   !> Calls functions in fortrandata.cpp

   public


contains
   subroutine fortrandata_setconstants()
      implicit none
      return
   end subroutine fortrandata_setconstants

   subroutine fortrandata_setmatrices()
      implicit none
      return
   end subroutine fortrandata_setmatrices

   subroutine fortrandata_setinputdata()
      implicit none
      return
   end subroutine fortrandata_setinputdata

   subroutine cudamdsim_initiateconstants()
      implicit none
      return
   end subroutine cudamdsim_initiateconstants

   subroutine cudamdsim_initiatematrices()
      implicit none
      return
   end subroutine cudamdsim_initiatematrices

   subroutine cudamdsim_measurementphase()
      implicit none
      return
   end subroutine cudamdsim_measurementphase


   subroutine cmdsim_initiateconstants()
      implicit none
      return
   end subroutine cmdsim_initiateconstants

   subroutine cmdsim_initiatefortran()
      implicit none
      return
   end subroutine cmdsim_initiatefortran

   subroutine cmdsim_measurementphase()
      implicit none
      return
   end subroutine cmdsim_measurementphase

   subroutine FortranData_Initiate(stt,btorque)
      use Parameters
      implicit none
      character(len=1), intent(in) :: STT       !< Treat spin transfer torque? (Y/N)
      real(dblprec), dimension(:,:,:), optional :: btorque !< Field from (m x dm/dr)
      return
   end subroutine FortranData_Initiate

end module nocuda
