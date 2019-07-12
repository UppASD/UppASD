!-------------------------------------------------------------------------------
! MODULE: HermiteLib
!> @brief Routines for Hermite interpolation polynomials
!> @author Pavel Bessarab
!-------------------------------------------------------------------------------
module HermiteLib
   use Parameters

   implicit none

   private
   public :: hermite_fit,spline_hermite_val

contains

   !-----------------------------------------------------------------------------
   ! SUBROUTINE: spline_hermite_set
   !> @author Pavel Bessarab
   !-----------------------------------------------------------------------------
   subroutine spline_hermite_set(ndata,tdata,coef)
      !
      implicit none
      !
      integer, intent(in) :: ndata
      real(dblprec), dimension(ndata), intent(in) :: tdata
      ! .. In/Out variables
      real(dblprec), dimension(4,ndata), intent(inout) :: coef !< Coefficients of the piecewise Hermite polynomials
      ! .. Local variables
      integer :: ii
      real(dblprec) :: divdif1, divdif3,dt
      !
      do ii = 1, ndata-1
         dt = tdata(ii+1) - tdata(ii)
         divdif1 = ( coef(1,ii+1) - coef(1,ii) ) / dt
         divdif3 = coef(2,ii) + coef(2,ii+1) - 2.0_dblprec* divdif1
         coef(3,ii) = ( divdif1 - coef(2,ii) - divdif3 ) / dt
         coef(4,ii) = divdif3 / (dt*dt)
      end do

      coef(3,ndata) = 0.0_dblprec
      coef(4,ndata) = 0.0_dblprec

   end subroutine spline_hermite_set

   !-----------------------------------------------------------------------------
   ! SUBROUTINE: spline_hermite_val
   !> @author Pavel Bessarab
   !-----------------------------------------------------------------------------
   subroutine spline_hermite_val(ndata,tdata,coef,tval,sval,dsval)
      !
      implicit none
      !
      integer, intent(in) :: ndata
      !
      real(dblprec), intent(in) :: tval
      real(dblprec), dimension(ndata), intent(in) :: tdata
      real(dblprec), dimension(4,ndata), intent(in) :: coef !< Coefficients of the piecewise Hermite polynomials
      ! .. Output variables
      real(dblprec), intent(out) :: sval, dsval
      ! .. Local variables
      integer :: left
      real(dblprec) :: dt

      !
      !  Find the interval [ TDATA(LEFT), TDATA(RIGHT) ] that contains
      !  or is nearest to TVAL.
      !
      call rvec_bracket (ndata, tdata, tval, left)
      !
      !  Evaluate the cubic polynomial.
      !
      dt = tval - tdata(left)

      sval = coef(1,left) + dt * ( coef(2,left) + dt * ( coef(3,left) + dt * coef(4,left) ) )
      dsval = coef(2,left) + dt*(2.0_dblprec*coef(3,left) + dt*3.0_dblprec*coef(4,left))

   end subroutine spline_hermite_val

   !-----------------------------------------------------------------------------
   ! SUBROUTINE: rvec_bracket
   !> @author Pavel Bessarab
   !-----------------------------------------------------------------------------
   subroutine rvec_bracket(n,x,xval,left)
      !
      implicit none
      !
      integer, intent(in) :: n
      real(dblprec), intent(in) :: xval
      real(dblprec), dimension(n), intent(in) :: x
      ! .. Output variables
      integer, intent(out) :: left
      ! .. Local variables
      integer :: ii
      !
      do ii = 2, n - 1
         if ( xval < x(ii) ) then
            left = ii - 1
            return
         end if
      end do

      left = n - 1

   end subroutine rvec_bracket

   !-----------------------------------------------------------------------------
   ! SUBROUTINE: hermite_fit
   !> @author Pavel Bessarab
   !-----------------------------------------------------------------------------
   subroutine hermite_fit(n,nn,x,y,dydx,xx,yy,dyy,coef)

      implicit none
      ! .. Input variables
      integer, intent(in) :: n,nn
      real(dblprec), dimension(n),intent(in) :: x
      real(dblprec), dimension(n),intent(in) :: y
      real(dblprec), dimension(n),intent(in) :: dydx
      ! .. Output variables
      real(dblprec), dimension(nn), intent(out) :: xx
      real(dblprec), dimension(nn), intent(out) :: yy
      real(dblprec), dimension(nn), intent(out) :: dyy
      real(dblprec), dimension(4,n), intent(out) :: coef !< Coefficients of the piecewise Hermite polynomials
      ! .. Local variables
      integer :: ii
      real(dblprec) :: dx

      do ii =1,n
         coef(1,ii) = y(ii)
         coef(2,ii) = dydx(ii)
      end do

      dx = (x(n)-x(1))/(nn-1)

      call spline_hermite_set(n,x,coef)

      do ii=1,nn
         xx(ii) = x(1) + (ii-1)*dx
         call spline_hermite_val(n,x,coef,xx(ii),yy(ii),dyy(ii))
      end do

   end subroutine hermite_fit

end module HermiteLib
