module nocuda
   use Parameters
   implicit none
   ! Initiate pointers for C/C++ implementation
   !> Calls functions in fortrandata.cpp

   public


contains
   subroutine fortrandata_setflags()
      implicit none
      return
   end subroutine fortrandata_setflags

   subroutine fortrandata_setconstants()
      implicit none
      return
   end subroutine fortrandata_setconstants

   subroutine fortrandata_sethamiltonian()
      implicit none
      return
   end subroutine fortrandata_sethamiltonian

   subroutine fortrandata_setlattice()
      implicit none
      return
   end subroutine fortrandata_setlattice

   subroutine fortrandata_setmeasurables()
      implicit none
      return
   end subroutine fortrandata_setmeasurables

   subroutine fortrandata_setinputdata()
      implicit none
      return
   end subroutine fortrandata_setinputdata

   subroutine gpusim_initiateconstants()
      implicit none
      return
   end subroutine gpusim_initiateconstants

   subroutine gpusim_initiatematrices()
      implicit none
      return
   end subroutine gpusim_initiatematrices

   subroutine gpusim_gpurunsimulation(whichsim, whichphase)
      implicit none
      integer, intent(in) :: whichsim, whichphase
      character(len=1), intent(in) :: gpu_mc_bf
      return
  end subroutine gpusim_gpurunsimulation

   subroutine gpusim_release()
      implicit none
      return
   end subroutine gpusim_release

 !  subroutine cmdsim_initiateconstants()
 !     implicit none
 !     return
 !  end subroutine cmdsim_initiateconstants

 !  subroutine cmdsim_initiatefortran()
 !     implicit none
 !     return
 !  end subroutine cmdsim_initiatefortran

  ! subroutine cmdsim_measurementphase()
   !   implicit none
    !  return
  ! end subroutine cmdsim_measurementphase

   subroutine FortranData_Initiate(stt,btorque, cc)
      use Parameters
      use Correlation_type
      implicit none
      character(len=1), intent(in) :: STT       !< Treat spin transfer torque? (Y/N)
      real(dblprec), dimension(:,:,:), optional :: btorque !< Field from (m x dm/dr)
      type(corr_t), intent(inout) :: cc !< Derived type for correlation data
      return
   end subroutine FortranData_Initiate


end module nocuda

