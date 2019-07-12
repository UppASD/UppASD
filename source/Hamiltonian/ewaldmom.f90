!> Dipolar interactions using Ewald summation
!> @author
!> J. Chico etc
!> @copyright
!> GNU Public License.
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

      ! Calculations to determine the volume of the cell
      cellx=maxval(coord(1,:))-minval(coord(1,:))
      celly=maxval(coord(2,:))-minval(coord(2,:))
      cellz=maxval(coord(3,:))-minval(coord(3,:))

      ! Volume of the cell
      Vol=cellx*celly*cellz

      Vinv=1.0_dblprec/(Vol)

      ! Factor of mu_0/4*pi
      fac=1.0d-7*mub/alat**3 ! This factor is setup for the sign of the field

      call SETUP_QCOORDS_EWALD(N1,N2,N3,KMAX,C1,C2,C3)

      EQdip=0.0_dblprec

      write(*,*) fac

      CFEwald=0.0_dblprec
      call COMPLEX_EWALD(Natom,Vinv,Ewald_alpha,coord)

      R_EWALD=0.0_dblprec

      call REAL_EWALD(Natom,RMAX,N1,N2,N3,alat,Ewald_alpha,coord)

      EW_facself1 = fac*(-4.0_dblprec*Ewald_alpha**3)/(3.0_dblprec*sqrt(pi))
      EW_facself2 = fac*(4.0_dblprec*pi*Vinv)/3.0_dblprec

      do iatom=1, Natom
         do jatom=1, Natom
            do nu=1,3
               do mu=1,3
                  EQdip(nu,mu,iatom,jatom) = R_EWALD(nu,mu,iatom,jatom)+fac*CFEwald(nu,iatom)*CCFEwald(mu,jatom)
               enddo
               EQdip(nu,nu,iatom,jatom) = EQdip(nu,nu,iatom,jatom) + EW_facself2
            enddo
         enddo
         EQdip(:,:,iatom,iatom) = EQdip(:,:,iatom,iatom) + EW_facself1
      enddo

      call allocate_Ewald(Natom,2)

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

      fac=1.0d-7/alat**3*mub

      R_EWALD=0.0_dblprec

      RTOT=RMAX(1)*RMAX(2)*RMAX(3)

      allocate(RCOORD(3,RTOT),stat=i_stat)
      call memocc(i_stat,product(shape(RCOORD))*kind(RCOORD),'RCOORD','setup_REWALD')
      allocate(IR(RTOT),stat=i_stat)
      call memocc(i_stat,product(shape(IR))*kind(IR),'IR','setup_IR')
      allocate(RRSCALAR(RTOT),stat=i_stat)
      call memocc(i_stat,product(shape(RRSCALAR))*kind(RRSCALAR),'RRSCALAR','setup_RRSCALAR')

      ITT=0
      open(ofileno,file='REWALD_COORD.out', position='append')
      ! Creation of coordinates for the cell images
      do I3=-(RMAX(3)-1)/2, RMAX(3)/2
         do I2=-(RMAX(2)-1)/2, RMAX(2)/2
            do I1=-(RMAX(1)-1)/2, RMAX(1)/2
               ITT=ITT+1
               RCOORD(1,ITT)=I1*N1
               RCOORD(2,ITT)=I2*N2
               RCOORD(3,ITT)=I3*N3
               RRSCALAR(ITT)=RCOORD(1,ITT)**2+RCOORD(2,ITT)**2+RCOORD(3,ITT)**2
               write(ofileno,'(i8,3es16.8,es16.8)') ITT,RCOORD(1:3,ITT),RRSCALAR(ITT)
            end do
         end do
      end do
      close(ofileno)

      ! Sort the coordinates of the cell images
      call qsort(RRSCALAR,IR,RTOT)

      ! Real part of the Ewald summation
      ERROR=1.0_dblprec
      do iatom=1, Natom
         do jatom=1, Natom
            IREP=1
            stop_signal=.false.
            ! Do loop over the images that stops when a treshold is reached in the tolerance
            do while ( IREP<=RTOT .and..not.stop_signal)
               if (RRSCALAR(IREP).ne.0) TESTING=ERROR_DET(R_EWALD(:,:,iatom,jatom))

               Rij(1)=coord(1,jatom)-coord(1,iatom)+RCOORD(1,IR(IREP))
               Rij(2)=coord(2,jatom)-coord(2,iatom)+RCOORD(2,IR(IREP))
               Rij(3)=coord(3,jatom)-coord(3,iatom)+RCOORD(3,IR(IREP))
               R2=Rij(1)**2+Rij(2)**2+Rij(3)**2
               RSCALAR=sqrt(R2)
               Rm5=R2**(-2.5)*fac

               do mu=1,3
                  do nu=1,3
                     R_EWALD(nu,mu,iatom,jatom) = R_EWALD(nu,mu,iatom,jatom) -1.0_dblprec*Rm5*Rij(mu)*Rij(nu)*C_funct(RSCALAR,R2,EALPHA)
                  enddo
                  R_EWALD(mu,mu,iatom,jatom) = R_EWALD(mu,mu,iatom,jatom) + Rm5*R2*B_funct(RSCALAR,R2,EALPHA)
               enddo

               if (RRSCALAR(IREP).ne.0) then
                  ERROR=TESTING-ERROR_DET(R_EWALD(:,:,iatom,jatom))
                  if (ERROR.lt.1e-10) then
                     stop_signal=.true.
                  endif
               endif
               if (RRSCALAR(IREP).eq.0) R_EWALD(:,:,iatom,iatom)=0.0_dblprec
               IREP=IREP+1
            enddo
         enddo
      enddo

      i_all=-product(shape(RRSCALAR))*kind(RRSCALAR)
      deallocate(RRSCALAR,stat=i_stat)
      call memocc(i_stat,i_all,'RRSCALAR','allocate_ewald')
      i_all=-product(shape(IR))*kind(IR)
      deallocate(IR,stat=i_stat)
      call memocc(i_stat,i_all,'IR','allocate_ewald')
      i_all=-product(shape(RCOORD))*kind(RCOORD)
      deallocate(RCOORD,stat=i_stat)
      call memocc(i_stat,i_all,'RCOORD','allocate_ewald')

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

      i=(0.0_dblprec,1.0_dblprec)

      fac1=4.0_dblprec*pi*Vinv

      EALP2=Ewald_alpha**2
      E2inv=1.0_dblprec/EALP2
      fac2=-1.0_dblprec*E2inv*(pi)**2

      ERROR_C=1.0_dblprec
      ERROR_CC=1.0_dblprec
      TESTING_C=0.0_dblprec
      TESTING_CC=0.0_dblprec
      stop_counter=.false.
      counter=1

      do jatom=1, Natom
         iq=1
         do while ( iq<=Nq .and. .not.stop_counter)
            if(dq(iq).ne.0.0_dblprec) then
               if (counter.gt.1) then
                  TESTING_C=CFEwald(1,jatom)**2+CFEwald(2,jatom)**2+CFEwald(3,jatom)**2
                  TESTING_CC=CCFEwald(2,jatom)**2+CCFEwald(2,jatom)**2+CFEwald(3,jatom)**2
               endif
               ! (exp(-(pi*q/Ewald_alpha)**2))**0.5
               expfacq=sqrt(exp(dq(iq)*fac2))
               ! sqrt(4*pi/(Vol*q**2))
               sqrtfac=sqrt(fac1/dq(iq))
               prefac=expfacq*sqrtfac
               !Dot product between the reciprocal and regular coordinates
               qrfac=qcoord(1,ia(iq))*coord(1,jatom)+qcoord(2,ia(iq))*coord(2,jatom)+qcoord(3,ia(iq))*coord(3,jatom)
               ! Setting the exponential factor as imaginary
               iqrfac=i*qrfac*2*pi
               ! sqrt(prefactor)*exp(i*2*pi*qrfac)
               epowqr=prefac*exp(iqrfac)
               Cepowqr=CONJG(epowqr)
               ! COMMENT: should CFE/CCFEwald really be real?
               do nu=1,3
                  CFEwald(nu,jatom)=real(CFEwald(nu,jatom)+qcoord(nu,ia(iq))*epowqr*expfacq)
                  CCFEwald(nu,jatom)=real(CCFEwald(nu,jatom)+qcoord(nu,ia(iq))*Cepowqr*expfacq)
               enddo
               counter=counter+1
               if (counter.eq.1) then
                  TESTING_C=CFEwald(1,jatom)**2+CFEwald(2,jatom)**2+CFEwald(3,jatom)**2
                  TESTING_CC=CCFEwald(2,jatom)**2+CCFEwald(2,jatom)**2+CFEwald(3,jatom)**2
               else
                  ERROR_C=abs(TESTING_C-(CFEwald(1,jatom)**2+CFEwald(2,jatom)**2+CFEwald(3,jatom)**2))
                  ERROR_CC=abs(TESTING_CC-(CCFEwald(2,jatom)**2+CCFEwald(2,jatom)**2+CFEwald(3,jatom)**2))
               endif
            endif

            if(ERROR_C.lt.1e-10.and. ERROR_CC.lt.1e-10) then
               stop_counter=.true.
            endif

            iq=iq+1

         enddo
      enddo

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

      B_funct = c_error_function(Arg) + (2.0_dblprec*Arg/sqrt(pi))*exp(-EALP2*R2)

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

      C_funct = 3.0_dblprec*c_error_function(Arg) + (2.0_dblprec*Arg/sqrt(pi))*(3.0_dblprec + 2.0_dblprec*EALP2*R2)*exp(-EALP2*R2)

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
      u = t - 0.5_dblprec

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

      ! Calculate reciprocal lattice vectors
      ! r1 = C2xC3
      r1(1)=C2(2)*C3(3)-C2(3)*C3(2)
      r1(2)=C2(3)*C3(1)-C2(1)*C3(3)
      r1(3)=C2(1)*C3(2)-C2(2)*C3(1)
      ! r2 = C3xC1
      r2(1)=C3(2)*C1(3)-C3(3)*C1(2)
      r2(2)=C3(3)*C1(1)-C3(1)*C1(3)
      r2(3)=C3(1)*C1(2)-C3(2)*C1(1)
      ! r3 = C1xC2
      r3(1)=C1(2)*C2(3)-C1(3)*C2(2)
      r3(2)=C1(3)*C2(1)-C1(1)*C2(3)
      r3(3)=C1(1)*C2(2)-C1(2)*C2(1)
      ! cell volume C1*(C2xC3)
      c1r1=C1(1)*r1(1)+C1(2)*r1(2)+C1(3)*r1(3)
      c2r2=C2(1)*r2(1)+C2(2)*r2(2)+C2(3)*r2(3)
      c3r3=C3(1)*r3(1)+C3(2)*r3(2)+C3(3)*r3(3)
      ! b1=(2pi)*r1/(C1*r1)
      b1(1)=r1(1)/c1r1
      b1(2)=r1(2)/c1r1
      b1(3)=r1(3)/c1r1
      ! b2=(2pi)*r2/(C1*r1)
      b2(1)=r2(1)/c2r2
      b2(2)=r2(2)/c2r2
      b2(3)=r2(3)/c2r2
      ! b3=(2pi)*r3/(C1*r1)
      b3(1)=r3(1)/c3r3
      b3(2)=r3(2)/c3r3
      b3(3)=r3(3)/c3r3
      !
      !print *,b3
      if(allocated(qcoord)) then
         i_all=-product(shape(qcoord))*kind(qcoord)
         deallocate(qcoord,stat=i_stat)
         call memocc(i_stat,i_all,'qcoord','setup_Ewald_q_coord')
      end if

      BZMAX=maxval(KMAX(:))

      Nq=N1*N2*N3*BZMAX
      ! write(*,*) N1*N2*N3*BZMAX
      allocate(ia(Nq),stat=i_stat)
      call memocc(i_stat,product(shape(ia))*kind(ia),'ia','setup_Ewald_q_coord')
      allocate(qcoord(3,Nq),stat=i_stat)
      call memocc(i_stat,product(shape(qcoord))*kind(qcoord),'qcoord','setup_Ewald_q_coord')
      allocate(dq(Nq),stat=i_stat)
      call memocc(i_stat,product(shape(dq))*kind(dq),'dq','setup_Ewald_q_coord')
      iq=0
      do BZ=1, BZMAX
         do zq=-(N3-1)/2,N3/2
            do yq=-(N2-1)/2,N2/2
               do xq=-(N1-1)/2,N1/2
                  iq=iq+1
                  qcoord(:,iq)=xq/(1.0_dblprec*N1/BZ)*b1+yq/(1.0_dblprec*N2/BZ)*b2+zq/(1.0_dblprec*N3/BZ)*b3
                  dq(iq)=qcoord(1,iq)**2+qcoord(2,iq)**2+qcoord(3,iq)**2
               end do
            end do
         end do
      enddo
      call qsort(dq,ia,Nq)
      qmin(1)=ia(1)
      qmin(2)=ia(2)

      open(ofileno,file='QPOINTS_EWALD.out',position='append')
      do iq=1, Nq

         write(ofileno,'(i8,i8,3es16.9,es16.8)')iq, ia(iq), qcoord(1:3,ia(iq)),dq(iq)

      enddo

      close(ofileno)

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
