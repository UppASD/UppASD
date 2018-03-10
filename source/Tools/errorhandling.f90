!====================================================================!
!> @brief
!> Intended use is container for any function dealing with error checking
!
!> @author
!> Thomas Nystrand
!> @copyright
!! Copyright (C) 2008-2018 UppASD group
!! This file is distributed under the terms of the
!! GNU General Public License. 
!! See http://www.gnu.org/copyleft/gpl.txt
!====================================================================!
module ErrorHandling
   use Parameters
   implicit none

   private

   public  ::                           &
      ErrorHandling_check_input_file,  &
      ErrorHandling_check_file_exists, &
      ErrorHandling_ERROR            , &
      ErrorHandling_missing

contains
   !---------------------------------------------------------------------------
   !> @brief
   !> Makes an initial check to make sure inpsd.dat exists and is not empty
   !
   !> @author
   !> Thomas Nystrand
   !---------------------------------------------------------------------------
   subroutine ErrorHandling_check_input_file()
      call ErrorHandling_check_file_exists('inpsd.dat','CD to appropriate folder or construct an inpsd.dat file')
   end subroutine ErrorHandling_check_input_file


   !---------------------------------------------------------------------------
   !> @brief
   !> Check if file exists, can be opened and prints
   !> an input error message if it does not exist
   !
   !> @author
   !> Thomas Nystrand
   !---------------------------------------------------------------------------
   subroutine ErrorHandling_check_file_exists(inputfile,message_missing) 
      use iso_fortran_env

      character(len=*),intent(in) :: inputfile
      character(len=*),optional,intent(in) :: message_missing 

      logical :: file_exist
      integer :: Read_Code, ifileno2
      character(10) :: teststring

      ifileno2=ifileno+1

      inquire(file=inputfile, exist=file_exist)
      if (file_exist) then
         open(ifileno2,file=inputfile)
         read(ifileno2,'(A)',advance='no',iostat=Read_Code) teststring
         if ( Read_Code /= 0 ) then
            if ( Read_Code == iostat_end ) then
               write(*,*)
               write(*,*) '*** ',inputfile,' file is empty ****'
               write(*,*) '*** execution terminated ***'
               stop
            else
               write ( *, '( / "read error: ", I0 )' )  Read_Code
               stop
            end if
         end if
         close(ifileno2)
      else
         write(*,*)
         write(*,*) '*** There is no file "',adjustl(trim(inputfile)),'" in current folder ****'
         if(present(message_missing)) write(*,*) '*** ',adjustl(trim(message_missing)),' ****'
         stop
      end if

   end subroutine ErrorHandling_check_file_exists



   !---------------------------------------------------------------------------
   ! Used to give a common error message structure
   !
   ! Written by: Thomas Nystrand
   !---------------------------------------------------------------------------
   subroutine ErrorHandling_ERROR(msg)
      character(len=*), intent(in) :: msg
      write(*,*) '                   **** ERROR ****                      '
      write(*,*) '========================================================'
      write(*,*) ' - ',msg 
      write(*,*) '========================================================'
      write(*,*) '             *** Terminating program ***                '
      stop
   end subroutine ErrorHandling_ERROR




   !---------------------------------------------------------------------------
   !> @brief
   !> Return string when calling unsupported features in this version    
   !!
   !---------------------------------------------------------------------------
   subroutine ErrorHandling_missing(feature)
      !
      implicit none
      !
      character(len=*),intent(in) :: feature

      write(*,999) feature
      STOP 
      999 format(1x,"This version does not support ",a,", stopping. Contact the developers for more info")
   end subroutine ErrorHandling_missing

end module ErrorHandling


