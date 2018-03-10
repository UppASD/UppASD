module ziggurat

   use stdtypes
   use mtprng

   type(mtprng_state) :: state_z

   integer(INT32) :: kn(128)
   real(IEEE64) :: fn(128)
   real(IEEE64) :: wn(128)

   real(IEEE64), parameter :: r = 3.442620D+00


contains

   function r4_nor ( )

      !*****************************************************************************80
      !
      !! R4_NOR returns a normally distributed single precision real value.
      !
      !  Discussion:
      !
      !    The value returned is generated from a distribution with mean 0 and
      !    variance 1.
      !
      !    The underlying algorithm is the ziggurat method.
      !
      !    Before the first call to this function, the user must call R4_NOR_SETUP
      !    to determine the values of KN, FN and WN.
      !
      !  Licensing:
      !
      !    This code is distributed under the GNU LGPL license.
      !
      !  Modified:
      !
      !    04 May 2008
      !
      !  Author:
      !
      !    Original C version by George Marsaglia, Wai Wan Tsang
      !    FORTRAN90 version by John Burkardt
      !
      !  Reference:
      !
      !    George Marsaglia, Wai Wan Tsang,
      !    The Ziggurat Method for Generating Random Variables,
      !    Journal of Statistical Software,
      !    Volume 5, Number 8, October 2000, seven pages.
      !
      !  Parameters:
      !
      !    Input/output, integer(INT32) JSR, the seed.
      !
      !    Input, integer(INT32) KN(128), data computed by R4_NOR_SETUP.
      !
      !    Input, real ( kind = 8 ) FN(128), WN(128), data computed by R4_NOR_SETUP.
      !
      !    Output, real ( IEEE64 ) R4_NOR, a normally distributed random value.
      !
      implicit none

      integer(INT32) :: hz
      integer(INT32) :: iz
      real(IEEE64) :: r4_nor
      real(IEEE64) :: value
      real(IEEE64) :: x
      real(IEEE64) :: y

      hz = mtprng_rand(state_z)
      iz = iand ( hz, 127 )

      if ( abs ( hz ) < kn(iz+1) ) then

         value = real ( hz, IEEE64 ) * wn(iz+1)

      else

         do

            if ( iz == 0 ) then

               do
                  !         x = - 0.2904764D+00 * log (mtprng_rand_real1(state))
                  !         y = - log(mtprng_rand_real1(state))
                  x = - 0.2904764D+00 * log (real(mtprng_rand_real1(state_z)))
                  y = - log(real(mtprng_rand_real1(state_z)))
                  if ( x*x <= y+y ) then
                     exit
                  end if
               end do

               if ( hz <= 0 ) then
                  value = - r - x
               else
                  value = + r + x
               end if

               exit

            end if

            x = real ( hz, IEEE64 ) * wn(iz+1)

            if (fn(iz+1)+mtprng_rand_real1(state_z)*(fn(iz)-fn(iz+1)) < exp(-0.5E+00*x*x)) then
               value = x
               exit
            end if

            !      hz = mtprng_rand_real1( state_z )
            hz = mtprng_rand ( state_z )
            iz = iand ( hz, 127 )

            if ( abs ( hz ) < kn(iz+1) ) then
               value = real ( hz, IEEE64 ) * wn(iz+1)
               exit
            end if

         end do

      end if

      r4_nor = value

      return
   end function r4_nor


   subroutine r4_nor_setup (inseed)

      !*****************************************************************************80
      !
      !! R4_NOR_SETUP sets data needed by R4_NOR.
      !
      !  Licensing:
      !
      !    This code is distributed under the GNU LGPL license.
      !
      !  Modified:
      !
      !    04 May 2008
      !
      !  Author:
      !
      !    Original C version by George Marsaglia, Wai Wan Tsang
      !    FORTRAN90 version by John Burkardt
      !
      !  Reference:
      !
      !    George Marsaglia, Wai Wan Tsang,
      !    The Ziggurat Method for Generating Random Variables,
      !    Journal of Statistical Software,
      !    Volume 5, Number 8, October 2000, seven pages.
      !
      !  Parameters:
      !
      !    Output, integer(INT32) KN(128), data needed by R4_NOR.
      !
      !    Output, real ( IEEE64 ) FN(128), WN(128), data needed by R4_NOR.
      !
      implicit none

      integer, intent(in) :: inseed
      real(IEEE64) :: dn
      integer(INT32) i
      ! real(IEEE64) :: fn(128)
      ! integer(INT32) kn(128)
      ! real(IEEE64) :: wn(128)
      real(IEEE64), parameter :: m1 = 2147483648.0D+00
      real(IEEE64) :: q
      real(IEEE64) :: tn
      real(IEEE64), parameter :: vn = 9.91256303526217D-03


      ! state=instate
      call mt_ran_init_zg(inseed)

      dn = 3.442619855899D+00
      tn = 3.442619855899D+00

      q = vn / exp ( - 0.5D+00 * dn * dn )

      kn(1) = int ( ( dn / q ) * m1 )
      kn(2) = 0

      wn(1) = real ( q / m1, IEEE64 )
      wn(128) = real ( dn / m1, IEEE64 )

      fn(1) = 1.0E+00
      fn(128) = real ( exp ( - 0.5D+00 * dn * dn ), IEEE64 )

      do i = 127, 2, -1
         dn = sqrt ( - 2.0D+00 * log ( vn / dn + exp ( - 0.5D+00 * dn * dn ) ) )
         kn(i+1) = int ( ( dn / tn ) * m1 )
         tn = dn
         fn(i) = real ( exp ( - 0.5D+00 * dn * dn ), IEEE64 )
         wn(i) = real ( dn / m1, IEEE64 )
         !print *,kn(i+1),fn(i),wn(i)
      end do

      return
   end subroutine r4_nor_setup

   !> Initialize the Mersenne Twister generator #1
   subroutine mt_ran_init_zg(seed)
      use mtprng
      integer, intent(in) :: seed  !< Seed number for PRNG

      call mtprng_init(seed, state)
   end subroutine mt_ran_init_zg


end module ziggurat
