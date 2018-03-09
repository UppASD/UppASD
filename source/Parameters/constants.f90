!>Physical constants
!! @todo Change to double precision once regression tests are passed
!> @copyright
!! Copyright (C) 2008-2018 UppASD group
!! This file is distributed under the terms of the
!! GNU General Public License. 
!! See http://www.gnu.org/copyleft/gpl.txt
module Constants
  use Parameters
  implicit none
  !.. Scalar parameters
  real(dblprec) :: gama=1.76d11
  real(dblprec) :: k_bolt=1.38d-23
  real(dblprec) :: k_bolt_ev=8.61734d-5
  real(dblprec) :: mub=9.274d-24
  real(dblprec) :: mu0=1.2566d-6
  real(dblprec) :: mry=2.1799d-21
  real(dblprec) :: hbar=6.62606d-34
  real(dblprec) :: hbar_mev=6.58211899d-13
  real(dblprec) :: Joule_ev=6.24150934d18
  real(dblprec) :: ry_ev=13.605698066
  real(dblprec) :: ev=1.602176565d-19
  real(dblprec) :: amu=1.660539040d-27
  real(dblprec) :: angstrom=1.0d-10
  real(dblprec),parameter :: pi=3.141592653589793d0
end module Constants
