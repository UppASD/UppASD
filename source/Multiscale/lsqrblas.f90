!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!     File lsqrblas.f90   (double precision)
!
!     This file contains the following BLAS routines
!        ddot, dnrm2, dscal
!     required by subroutines LSQR and Acheck.
!
! 24 Sep 2007: All parameters declared with correct intent
!              to avoid compiler warnings.
! 24 Oct 2007: Use real(8) instead of double precision or -r8.
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!! DDOT forms the dot product of two vectors.
!
!  Discussion:
!    This routine uses double precision real arithmetic.
!    This routine uses unrolled loops for increments equal to one.
!
!  Modified:
!    16 May 2005
!
!  Author:
!    Jack Dongarra
!    Fortran90 translation by John Burkardt.
!
!  Reference:
!    Jack Dongarra, Jim Bunch, Cleve Moler, Pete Stewart,
!    LINPACK User´s Guide,
!    SIAM, 1979,
!    ISBN13: 978-0-898711-72-1,
!    LC: QA214.L56.
!
!    Charles Lawson, Richard Hanson, David Kincaid, Fred Krogh,
!    Algorithm 539, 
!    Basic Linear Algebra Subprograms for Fortran Usage,
!    ACM Transactions on Mathematical Software, 
!    Volume 5, Number 3, September 1979, pages 308-323.
!
!  Parameters:
!
!    Input, integer N, the number of entries in the vectors.
!
!    Input, real ( kind = 8 ) DX(*), the first vector.
!
!    Input, integer INCX, the increment between successive entries in DX.
!
!    Input, real ( kind = 8 ) DY(*), the second vector.
!
!    Input, integer INCY, the increment between successive entries in DY.
!
!    Output, real ( kind = 8 ) DDOT, the sum of the product of the 
!    corresponding entries of DX and DY.


      real(8) function ddot(n,dx,incx,dy,incy)

      implicit none
      integer, intent(in) :: n,incx,incy
      real(8), intent(in) :: dx(*),dy(*)

      real(8)  dtemp
      integer  i,ix,iy,m

      ddot = 0.0d0
      dtemp = 0.0d0
      if ( n <= 0 ) then
         return
      end if

!  Code for unequal increments or equal increments
!  not equal to 1.

      if ( incx /= 1 .or. incy /= 1 ) then

         if ( 0 <= incx ) then
            ix = 1
         else
            ix = ( - n + 1 ) * incx + 1
         end if

         if ( 0 <= incy ) then
            iy = 1
         else
            iy = ( - n + 1 ) * incy + 1
         end if

         do i = 1, n
            dtemp = dtemp + dx(ix) * dy(iy)
            ix = ix + incx
            iy = iy + incy
         end do

!  Code for both increments equal to 1.

        else

           m = mod ( n, 5 )

           do i = 1, m
              dtemp = dtemp + dx(i) * dy(i)
           end do

           do i = m+1, n, 5
              dtemp = dtemp + dx(i)*dy(i) + dx(i+1)*dy(i+1) + dx(i+2)*dy(i+2) &
                                          + dx(i+3)*dy(i+3) + dx(i+4)*dy(i+4)
           end do

        end if

        ddot = dtemp
        return
end function ddot

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!*****************************************************************************80
!
!! DNRM2 returns the euclidean norm of a vector.
!
!  Discussion:
!    This routine uses real(8) real arithmetic.
!     DNRM2 ( X ) = sqrt ( X´ * X )
!
!  Modified:
!    16 May 2005
!
!  Author:
!    Sven Hammarling
!    Fortran90 translation by John Burkardt.
!
!  Reference:
!    Jack Dongarra, Jim Bunch, Cleve Moler, Pete Stewart,
!    LINPACK User´s Guide,
!    SIAM, 1979,
!    ISBN13: 978-0-898711-72-1,
!    LC: QA214.L56.
!
!    Charles Lawson, Richard Hanson, David Kincaid, Fred Krogh,
!    Algorithm 539, 
!    Basic Linear Algebra Subprograms for Fortran Usage,
!    ACM Transactions on Mathematical Software,
!    Volume 5, Number 3, September 1979, pages 308-323.
!
!  Parameters:
!
!    Input, integer N, the number of entries in the vector.
!
!    Input, real ( kind = 8 ) X(*), the vector whose norm is to be computed.
!
!    Input, integer INCX, the increment between successive entries of X.
!
!    Output, real ( kind = 8 ) DNRM2, the Euclidean norm of X.
!

      real(8) function dnrm2 ( n, x, incx)

      implicit none
      integer, intent(in) :: n,incx
      real(8), intent(in) :: x(*)

      integer  ix
      real(8)  ssq,absxi,norm,scale

      if ( n < 1 .or. incx < 1 ) then
         norm  = 0.d0
      else if ( n == 1 ) then
         norm  = abs ( x(1) )
      else
         scale = 0.d0
         ssq = 1.d0

         do ix = 1, 1 + ( n - 1 )*incx, incx
            if ( x(ix) /= 0.d0 ) then
               absxi = abs ( x(ix) )
               if ( scale < absxi ) then
                  ssq = 1.d0 + ssq * ( scale / absxi )**2
                  scale = absxi
               else
                  ssq = ssq + ( absxi / scale )**2
               end if
            end if
         end do
         norm  = scale * sqrt ( ssq )
      end if

      dnrm2 = norm
      return
end function dnrm2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! DSCAL scales a vector by a constant.
!
!  Discussion:
!    This routine uses real(8) arithmetic.
!
!  Modified:
!    08 April 1999
!
!  Author:
!    Jack Dongarra
!    Fortran90 translation by John Burkardt.
!
!  Reference:
!    Jack Dongarra, Jim Bunch, Cleve Moler, Pete Stewart,
!    LINPACK User´s Guide,
!    SIAM, 1979,
!    ISBN13: 978-0-898711-72-1,
!    LC: QA214.L56.
!
!    Charles Lawson, Richard Hanson, David Kincaid, Fred Krogh,
!    Algorithm 539, 
!    Basic Linear Algebra Subprograms for Fortran Usage,
!    ACM Transactions on Mathematical Software,
!    Volume 5, Number 3, September 1979, pages 308-323.
!
!  Parameters:
!
!    Input, integer N, the number of entries in the vector.
!
!    Input, real ( kind = 8 ) SA, the multiplier.
!
!    Input/output, real ( kind = 8 ) X(*), the vector to be scaled.
!
!    Input, integer INCX, the increment between successive entries of X.
!

      subroutine  dscal(n,sa,x,incx)

      implicit none
      integer, intent(in)    :: n,incx
      real(8), intent(in)    :: sa
      real(8), intent(inout) :: x(*)

      integer                   i,ix,m

      if ( n <= 0 ) then
         return
      else if ( incx == 1 ) then
         m = mod ( n, 5 )
         x(1:m) = sa * x(1:m)

         do i = m+1, n, 5
            x(i)   = sa * x(i)
            x(i+1) = sa * x(i+1)
            x(i+2) = sa * x(i+2)
            x(i+3) = sa * x(i+3)
            x(i+4) = sa * x(i+4)
         end do
      else
         if ( 0 <= incx ) then
            ix = 1
         else
            ix = ( - n + 1 ) * incx + 1
         end if

         do i = 1, n
            x(ix) = sa * x(ix)
            ix = ix + incx
         end do

      end if

      return
end subroutine dscal
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
