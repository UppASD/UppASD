!> Global parameters
!> /note Keep the contents of this module to a minimum
module Parameters

   implicit none

   integer, parameter :: snglprec = selected_real_kind(6, 37)  !< define precision for single reals
   integer, parameter :: dblprec = selected_real_kind(15, 307)  !< define precision for double reals
   !integer, parameter :: dblprec = selected_real_kind(6, 37)  !< define precision for single reals
   integer, parameter :: qdprec = selected_real_kind(33, 4931)  !< define precision for quad reals
   integer, parameter :: long = selected_int_kind(range(1)*2)   !< define long integer
   real(dblprec) :: dbl_tolerance=1e-14 !!! Note this is used for both ne and eq as > or < comparisons

   integer, parameter :: ifileno = 55 !< File handle number for input files
   integer, parameter :: ifileno2 = 56 !< File handle number for input files
   integer, parameter :: ofileno = 66 !< File handle number for output files
   integer, parameter :: mfileno = 77 !< File handle number for memory log file

   integer, parameter :: ofileno2 = 67 !< File handle number for output files
   integer, parameter :: ofileno3 = 68 !< File handle number for output files

   public
!!#ifdef ( __PATHSCALE__
!!!#if (defined __PATHSCALE__) || (defined __PGIF90__ ) || ( defined __GFORTRAN__ )
!!!#if (defined __PATHSCALE__) || (defined __PGIF90__ ) 
!!!contains
!!!
!!!   real(dblprec) function norm2(vec)
!!!      !
!!!      implicit none
!!!      !
!!!      real(dblprec), dimension(:), intent(in) :: vec
!!!      !
!!!      norm2=sqrt(sum(vec*vec))
!!!      return
!!!      !
!!!   end function norm2
!!!#endif
!!! #if (defined __PATHSCALE__) || (defined __PGIF90__ ) 
!!! contains
!!! 
!!!    function compiler_version() result(cname)
!!!       !
!!!       implicit none
!!!       !
!!!       character*1 :: cname
!!!       !
!!!       cname=' '
!!!       return
!!!       !
!!!    end function compiler_version
!!! #endif

end module
