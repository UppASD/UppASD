!-------------------------------------------------------------------------------
! MODULE: Constants
!> @brief Physical constants
!> @details The magnitude of the parameters has been chosen so that it mtaches
!> the values given in the NIST database with the CODATA recommended 2014 values
!> https://physics.nist.gov/cuu/Constants/
!> @copyright
!> GNU Public License.
!-------------------------------------------------------------------------------
module Constants
   use Parameters
   implicit none
   !.. Scalar parameters
   real(dblprec) :: gama         = 1.760859644d11     ! s^(-1)*T^(-1)
   real(dblprec) :: k_bolt       = 1.38064852d-23     ! J/K
   real(dblprec) :: k_bolt_ev    = 8.6173303d-5       ! eV/K
   real(dblprec) :: mub          = 9.274009994d-24    ! J/T
   real(dblprec) :: mu0          = 1.2566370614d-6    ! N/A^2
   real(dblprec) :: mry          = 2.179872325d-21    ! J
   real(dblprec) :: hbar_mev     = 6.582119514d-13    ! meV*s
   real(dblprec) :: Joule_ev     = 6.241509126d18     ! eV
   real(dblprec) :: ry_ev        = 13.605693009_dblprec     ! eV
   real(dblprec) :: hbar         = 1.054571800e-34    ! J*s
   real(dblprec) :: ev           = 1.6021766208d-19   ! J
   real(dblprec) :: amu          = 1.660539040d-27    ! kg
   real(dblprec) :: angstrom     = 1.0d-10            ! m
   real(dblprec) :: a0           = 0.52917721067d-10  ! m
   real(dblprec) :: g_e_abs      = 2.00231930436182
   real(dblprec),parameter :: pi = 3.141592653589793_dblprec
end module Constants
