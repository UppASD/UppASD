!> Dipolar interactions using Ewald summation
  !> @author
  !> J. Chico etc
!> @copyright
!! Copyright (C) 2008-2018 UppASD group
!! This file is distributed under the terms of the
!! GNU General Public License. 
!! See http://www.gnu.org/copyleft/gpl.txt
module Ewaldmom
   use Parameters
   use Profiling
   use Constants

   implicit none

   real(dblprec), dimension(:,:), allocatable :: qcoord !< q-points
   real(dblprec), dimension(:,:), allocatable :: CFEwald
   real(dblprec), dimension(:,:), allocatable :: CCFEwald
   real(dblprec), dimension(:,:,:,:), allocatable :: R_EWALD
   real(dblprec), dimension(:,:,:,:), allocatable :: EQdip
   integer, dimension(:), allocatable :: ia
   real(dblprec), dimension(:), allocatable :: dq

   real(dblprec) :: EW_facself1, EW_facself2

   integer :: Nq
   integer,dimension(2) :: qmin  !< Index of smallest wave vector and gamma point

   private

   public :: SETUP_QCOORDS_EWALD, allocate_ewald, Ewald_Qdip, EQdip

   contains

!--------------------------------------------------------------------------
  !
  ! DESCRIPTION
  !> @brief
  !> Dipolar interactions in periodic systems using Ewald summation
  !---------------------------------------------------------------------------------
subroutine EWALD_QDIP(Natom,N1,N2,N3,C1,C2,C3,RMAX,KMAX,alat,Ewald_alpha,coord,EQdip)

   implicit none

   integer, intent(in) :: Natom !< Number of atoms in the system
   integer, intent(in) :: N1
   integer, intent(in) :: N2
   integer, intent(in) :: N3
   integer, dimension(3), intent(in) :: RMAX ! Maximum number of cell repetitions in the X-direction
   integer, dimension(3), intent(in) :: KMAX

   real(dblprec), intent(in) :: alat ! Lattice parameter
   real(dblprec), intent(in) :: Ewald_alpha ! Ewald parameter
   real(dblprec), dimension(3), intent(in) :: C1
   real(dblprec), dimension(3), intent(in) :: C2
   real(dblprec), dimension(3), intent(in) :: C3
   real(dblprec), dimension(3,Natom), intent(in) :: coord ! Atom coordinates

   real(dblprec), dimension(3,3,Natom,Natom), intent(out) :: EQdip ! The Ewald tensor

   real(dblprec) :: EW_facself1, EW_facself2, fac ! Factors for the self interaction and surface terms
   real(dblprec) :: cellx, celly, cellz, Vol, Vinv ! Volume variables

   integer :: iatom, jatom, mu, nu

            call ErrorHandling_missing('Kinetic MC')

end subroutine EWALD_QDIP

!> Real part of Ewald summation
subroutine REAL_EWALD(Natom,RMAX,N1,N2,N3,alat,EALPHA,coord)

   use Sorting, only : qsort

   implicit none

   integer, intent(in) :: N1 ! Number of repetitions in the X-direction of the unit cell
   integer, intent(in) :: N2 ! Number of repetitions in the Y-direction of the unit cell
   integer, intent(in) :: N3 ! Number of repetitions in the Z-direction of the unit cell
   integer, intent(in) :: Natom !< Number of atoms in the simulation cell
   integer, dimension(3), intent(in) :: RMAX ! Number of cells that will be repeated in the Ewald summation

   real(dblprec), intent(in) :: alat ! Lattice parameter
   real(dblprec), intent(in) :: EALPHA ! Ewald parameter
   real(dblprec), dimension(3,Natom), intent(in) :: coord ! Coordinates of the atoms

   integer, dimension(:), allocatable :: IR
   real(dblprec), dimension(:), allocatable :: RRSCALAR
   real(dblprec), dimension(:,:), allocatable :: RCOORD
   integer :: ITT, I1,I2,I3,RTOT, IREP

   real(dblprec) :: RSCALAR, R2, Rm5, fac
   real(dblprec), dimension(3) :: Rij
   real(dblprec) :: TESTING, ERROR

   logical :: stop_signal

   integer :: i_stat, i_all
   integer :: iatom, jatom, nu, mu

            call ErrorHandling_missing('Kinetic MC')

end subroutine REAL_EWALD

!> Error estimation
real(dblprec) function ERROR_DET(MATRIX)

   real(dblprec), dimension(3,3), intent(in) :: MATRIX

   ERROR_DET = MATRIX(1,1)*(MATRIX(2,2)*MATRIX(3,3) - MATRIX(3,2)*MATRIX(2,3)) &
             + MATRIX(1,2)*(MATRIX(3,1)*MATRIX(2,3) - MATRIX(2,1)*MATRIX(3,3)) &
             + MATRIX(1,3)*(MATRIX(2,1)*MATRIX(3,2) - MATRIX(3,1)*MATRIX(2,2))

end function ERROR_DET

!> Complex fourier part of the Ewald summation
subroutine COMPLEX_EWALD(Natom,Vinv,Ewald_alpha,coord)

   implicit none

   integer, intent(in) :: Natom !< Number of atoms in the system

   real(dblprec), intent(in) :: Vinv ! Inverse Volume of the simulation cell
   real(dblprec), intent(in) :: Ewald_alpha ! Ewald parameter
   real(dblprec), dimension(3,Natom), intent(in) :: coord ! Coordinates of the atoms

   real(dblprec) :: prefac, qrfac, EALP2, E2inv
   real(dblprec) :: fac2, expfacq, fac1, sqrtfac
   real(dblprec) :: TESTING_C, TESTING_CC, ERROR_C, ERROR_CC

   complex(dblprec) :: i, iqrfac, epowqr, Cepowqr
   logical :: stop_counter
   integer :: iq, jatom, nu, counter

            call ErrorHandling_missing('Kinetic MC')

end subroutine COMPLEX_EWALD

!> Function needed for the real part of the Ewald summation
real(dblprec) function B_funct(Rscalar,R2,Ewald_alpha)

   implicit none

   real(dblprec), intent(in) :: Ewald_alpha ! Ewald parameter
   real(dblprec), intent(in) :: Rscalar ! Distance between the moments
   real(dblprec), intent(in) :: R2

   real(dblprec) :: EALP2, Arg ! Local variables

   EALP2=(Ewald_alpha**2)
   Arg=Ewald_alpha*Rscalar

   B_funct = c_error_function(Arg) + (2.0d0*Arg/sqrt(pi))*exp(-EALP2*R2)

end function B_funct

!> Function needed for the real part of the Ewald summation
real(dblprec) function C_funct(Rscalar,R2,Ewald_alpha)

   implicit none

   real(dblprec), intent(in) :: Ewald_alpha
   real(dblprec), intent(in) :: Rscalar
   real(dblprec), intent(in) :: R2

   real(dblprec) ::  EALP2, Arg ! Local variables

   EALP2=(Ewald_alpha**2)
   Arg=(Ewald_alpha*Rscalar)

   C_funct = 3.0d0*c_error_function(Arg) + (2.0d0*Arg/sqrt(pi))*(3.0d0 + 2.0d0*EALP2*R2)*exp(-EALP2*R2)

end function C_funct


!> Complementary error function
real(dblprec) function c_error_function(x)
   ! Complementary Error Function obtained from Takuya OOURA
   ! Copyright(C) 1996 Takuya OOURA (email: ooura@mmm.t.u-tokyo.ac.jp).
   ! The copyright states that one may use, copy, modify this code for any purpose and without fee. You may distribute this ORIGINAL package.
   ! Modifications are done to the package to make it UppASD3.0 compliant
      implicit none

      real(dblprec), intent(in) :: x ! Input parameter for the function
      real(dblprec) :: t, u, y  ! Auxiliarly variables

      real(dblprec), parameter :: pa = 3.97886080735226000d+00, p0 = 2.75374741597376782d-01, p1 = 4.90165080585318424d-01, &
                  p2 = 7.74368199119538609d-01, p3 = 1.07925515155856677d+00, p4 = 1.31314653831023098d+00, &
                  p5 = 1.37040217682338167d+00, p6 = 1.18902982909273333d+00, p7 = 8.05276408752910567d-01, &
                  p8 = 3.57524274449531043d-01, p9 = 1.66207924969367356d-02, p10 = -1.19463959964325415d-01, &
                  p11 = -8.38864557023001992d-02

      real(dblprec), parameter :: p12 = 2.49367200053503304d-03, p13 = 3.90976845588484035d-02, p14 = 1.61315329733252248d-02, &
                  p15 = -1.33823644533460069d-02, p16 = -1.27223813782122755d-02, p17 = 3.83335126264887303d-03, &
                  p18 = 7.73672528313526668d-03, p19 = -8.70779635317295828d-04, p20 = -3.96385097360513500d-03, &
                  p21 = 1.19314022838340944d-04, p22 = 1.27109764952614092d-03

      t = pa / (pa + abs(x))
      u = t - 0.5d0

      y = (((((((((p22*u+p21)*u+p20)*u+p19)*u+p18)*u+p17)*u+p16)*u+p15)*u+p14)*u+p13)*u+p12

      y = ((((((((((((y*u+p11)*u+p10)*u+p9)*u+p8)*u+ p7)*u+p6)*u+p5)*u+p4)*u+p3)*u+p2)*u+p1)*u+p0)*t*exp(-x*x)

      if (x .lt. 0) y = 2 - y

      c_error_function = y

      return
end function c_error_function


!> Automatic setup of q-point mesh for calculating the reciprocal space Ewald summation
  subroutine SETUP_QCOORDS_EWALD(N1,N2,N3,KMAX,C1,C2,C3)

    use Sorting, only : qsort
    !
    !
    implicit none
    !
    integer, intent(in) :: N1  !< Number of cell repetitions in x direction
    integer, intent(in) :: N2  !< Number of cell repetitions in y direction
    integer, intent(in) :: N3  !< Number of cell repetitions in z direction

    integer, dimension(3), intent(in) :: KMAX

    real(dblprec), dimension(3), intent(in) :: C1 !< First lattice vector
    real(dblprec), dimension(3), intent(in) :: C2 !< Second lattice vector
    real(dblprec), dimension(3), intent(in) :: C3 !< Third lattice vector
    !
    integer :: iq,xq,yq,zq
    integer :: i_stat, i_all, BZ, BZMAX
    real(dblprec), dimension(3) :: b1,r1
    real(dblprec), dimension(3) :: b2,r2
    real(dblprec), dimension(3) :: b3,r3
    real(dblprec) :: c1r1, c2r2, c3r3

             call ErrorHandling_missing('Kinetic MC')


  end subroutine SETUP_QCOORDS_EWALD


!> Allocation of Ewald arrays
subroutine allocate_ewald(Natom,flag)

     implicit none

     integer, intent(in) :: flag
     integer, intent(in) :: Natom

     integer :: i_stat, i_all

  if (flag==1) then
     allocate(CFEwald(3,Natom),stat=i_stat)
     call memocc(i_stat,product(shape(CFEwald))*kind(CFEwald),'CFEwald','allocate_ewald')
     allocate(CCFEwald(3,Natom),stat=i_stat)
     call memocc(i_stat,product(shape(CCFEwald))*kind(CCFEwald),'CCFEwald','allocate_ewald')
     allocate(R_EWALD(3,3,Natom,Natom),stat=i_stat)
     call memocc(i_stat,product(shape(R_EWALD))*kind(R_EWALD),'REWALD','allocate_ewald')
     allocate(EQdip(3,3,Natom,Natom),stat=i_stat)
     call memocc(i_stat,product(shape(EQdip))*kind(EQdip),'EQdip','allocate_ewald')

  else if (flag==2) then
    i_all=-product(shape(CFEwald))*kind(CFEwald)
    deallocate(CFEwald,stat=i_stat)
    call memocc(i_stat,i_all,'CFEwald','allocate_ewald')

    i_all=-product(shape(R_EWALD))*kind(R_EWALD)
    deallocate(R_EWALD,stat=i_stat)
    call memocc(i_stat,i_all,'R_EWALD','allocate_ewald')

    i_all=-product(shape(CCFEwald))*kind(CCFEwald)
    deallocate(CCFEwald,stat=i_stat)
    call memocc(i_stat,i_all,'CCFEwald','allocate_ewald')

    if(allocated(qcoord)) then
       i_all=-product(shape(qcoord))*kind(qcoord)
       deallocate(qcoord,stat=i_stat)
       call memocc(i_stat,i_all,'qcoord','setup_Ewald_q_coord')
    end if

   if(allocated(dq)) then
       i_all=-product(shape(dq))*kind(dq)
       deallocate(dq,stat=i_stat)
       call memocc(i_stat,i_all,'dq','setup_Ewald_q_coord')
   endif

   if(allocated(ia)) then
       i_all=-product(shape(ia))*kind(ia)
       deallocate(ia,stat=i_stat)
       call memocc(i_stat,i_all,'ia','setup_Ewald_q_coord')
   endif

   else
     i_all=-product(shape(EQdip)*kind(EQdip))
     deallocate(EQdip,stat=i_stat)
     call memocc(i_stat,i_all,'EQdip','allocate_ewald')

  endif


end subroutine allocate_ewald


end module
