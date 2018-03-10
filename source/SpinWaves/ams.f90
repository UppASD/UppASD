!> Routines for calculating adiabatic magnon spectra (AMS), and magnon density of states
!> @author
!> A. Bergman, L. Bergqvist, J. Chico, etc
!> @copyright
!! Copyright (C) 2008-2018 UppASD group
!! This file is distributed under the terms of the
!! GNU General Public License. 
!! See http://www.gnu.org/copyleft/gpl.txt
module AMS
   !
   ! AMS - Adiabatic Magnon Spectra calculation
   ! The adiabatic magnon spectra is written to ams.SIMID.out, and the
   ! energies are sorted by nsize on each row.
   !
   ! Theory taken from Kübler, Theory of Itinerant Electron Magnetism (2009)
   ! and Halilov et ad, Phys. Rev. B, 58:293, 1998. However the A-matrix proposed
   ! in Kübler has been modified to function with the j-exchange already implemented
   ! in UppASD.
   !
   ! The module contains the main subroutine calculate_ams(), using calc_j,
   ! eigenvalue_calculation and printEnergies.
   !
   use Parameters
   use Fileparser
   use Constants
   use InputData,          only : N1,N2,N2,N3,NA,NT,Natom,simid, do_dm, hfield,gsconf_num

   use Correlation,        only : q,nq
   use Hamiltoniandata,    only : nlist, nlistsize, ncoup, dmlist, dmlistsize, dm_vect, taniso, eaniso, kaniso
   use Systemdata,         only : anumb,coord, Landeg
   use Momentdata,         only : mmom, emom, emomM
   use ChemicalData,       only : achem_ch

   implicit none


   ! Adiabatic Magnon Spectra calculation flag
   character(LEN=1) :: do_ams                  !< Calculate AMS (N/Y)
   character(LEN=1) :: do_magdos               !< Calculate AMS DOS (N/Y/F)
   character(len=35) :: magdosfile             !< AMS DOS file
   real(dblprec)    :: magdos_sigma            !< Frequency broadening of AMS DOS (in meV)
   integer          :: magdos_freq             !< Number of frequencies of AMS DOS
   integer          :: magdos_lfreq             !< Low cut-off frequency of AMS DOS (in meV)
   integer          :: magdos_hfreq             !< High cut-off frequency of AMS DOS (in meV)

   integer          :: ams_hist_size=2000

   real(dblprec), allocatable, dimension(:,:) :: magdos
   real(dblprec) :: tcmfa,tcrpa,msat

   private

   public :: magdos,tcmfa,tcrpa,msat, magdos_freq, do_ams, do_magdos
   public :: magdosfile, magdos_sigma
   public :: calculate_ams,calc_j,eigenvalue_calculation_lapack,printEnergies,magdos_calc, read_parameters_ams, init_ams
   public :: calculate_random_ams

contains

   !--------------------------------------------------------------------------
   !
   ! DESCRIPTION
   !> @brief
   !> Main driver for adiabatic magnon spectra (AMS) in collinear magnets.
   !! Based on linear spin wave theory as described in ....
   !---------------------------------------------------------------------------------
   subroutine calculate_ams()
      implicit none
      integer :: mu
      integer :: nu
      integer :: lambda
      integer :: i,j
      integer :: I1,I2,I3 !< indices to find suitable cell to perform calculations

      integer :: countstart !< index for last atom "before" the cell where calculations are performed

      real(dblprec), dimension(3) :: zeroVector !(0,0,0) vector  to calculate J(0)
      real(dblprec), allocatable, dimension(:,:) :: wres,jqres
      real(dblprec) :: fc,fc2,mtmp

      real(dblprec) :: mag_sign_mu, mag_sign_nu, mag_sign_lambda
      real(dblprec), dimension(3) :: order_vector
      complex(dblprec), allocatable, dimension(:,:,:) :: A !< matrix whose eigenvalues are the sought energies fro AMS
      complex(dblprec), allocatable, dimension(:,:,:) :: eigv !< matrix whose eigenvalues are the sought energies fro AMS
      complex(dblprec), allocatable, dimension(:,:,:) :: B !< matrix whose eigenvalues are the sought energies for J(q)

      ! for output printing
      character(LEN = 16) :: output_file
      character(LEN = 18) :: output_file2
      character(LEN = 18) :: output_file3
      character(LEN = 19) :: output_file4
      integer :: iq

      ! Factors for mRy energy conversion
      fc = mry/mub
      fc2 = 2*mry/mub

      ! if no q allocated then no AMS calculations are performed
      if (.not. allocated(q)) then
         write(*,*) 'No q-file allocated. No Adiabatic Magnon Spectrum calculation is performed'
         goto 10 ! jumps to end of subroutine
      end if

      allocate(A(NA,NA,nq))
      allocate(eigv(NA,NA,nq))
      allocate(B(NA,NA,nq))
      allocate(wres(NA,nq))
      allocate(jqres(NA,nq))
      A = 0.0_dblprec
      B = 0.0_dblprec
      zeroVector = 0.0_dblprec

      output_file = 'ams.'//simid//'.out'
      output_file2 = 'jqams.'//simid//'.out'
      output_file3 = 'evams.'//simid//'.out'
      output_file4 = 'magdos.'//simid//'.out'

      ! since j(q) is depending on positions (scalar product with q-vector), it is desired
      ! to perform calculations on the lattice in the middle of the system, to avoid
      ! problems caused by periodicy. Calculations are done for the lattice starting with
      ! atom number countstart+1
      !
      I1 = N1/2
      I2 = N2/2
      I3 = N3/2

      countstart = 0+I1*NA+I2*N1*NA+I3*N2*N1*NA
      ! this routine might be redundant
      if (anumb(countstart+1)/=1) then
         do i = 1,NA
            if (anumb(countstart+i+1)==1) then
               countstart = countstart+i
               exit
            else if (anumb(countstart-i+1)==1) then
               countstart = countstart-i
               exit
            end if
         end do
      end if
      !
      do i=1,na
         mtmp=0.d0
         do j=1,nlistsize(i+countstart)
            mtmp=mtmp+mmom(nlist(j,i+countstart),1)
         enddo
         mtmp=mtmp/nlistsize(i+countstart)
         msat=msat+mtmp/na
      enddo

      order_vector(1)=1.0d0/sqrt(3.0d0)
      order_vector(2)=1.0d0/sqrt(3.0d0)
      order_vector(3)=1.0d0/sqrt(3.0d0)
      ! A-matrix calculation, one for each q
      write (*,'(1x,a)',advance="no") "Mounting A-matrix for AMS calculation"
      A=(0.0d0,0.0d0)
      do i=1,nq
         do mu=1,NA
            mag_sign_mu=sign(1.0d0,sum(emom(1:3,mu+countstart,1)*order_vector))
            do nu=1,NA
               mag_sign_nu=sign(1.0d0,sum(emom(1:3,nu+countstart,1)*order_vector))
               A(mu,nu,i) = -calc_j(mu,nu,q(:,i),countstart)*mmom(nu+countstart,1)*mag_sign_nu
               B(mu,nu,i) = -A(mu,nu,i)*mmom(nu+countstart,1)
               if (mu==nu) then
                  do lambda = 1,NA
                     mag_sign_lambda=sign(1.0d0,sum(emom(1:3,lambda+countstart,1)*order_vector))
                     A(mu,mu,i)=A(mu,mu,i) + calc_j(mu,lambda,zeroVector,countstart)* &
                        mmom(lambda+countstart,1)*mag_sign_lambda
                  end do
                  A(mu,mu,i)=A(mu,mu,i) + calc_ani(mu,mu,countstart)
                  A(mu,mu,i)=A(mu,mu,i) + sum(hfield*emomm(1:3,mu+countstart,1))*mag_sign_nu*0.5d0
               end if
               A(mu,nu,i)=A(mu,nu,i)*Landeg(mu)*Landeg(nu) !LB Needs to be fixed,also multiply with landeg_glob
               B(mu,nu,i)=B(mu,nu,i)*Landeg(mu)*Landeg(nu)
            end do
         end do
      end do


      write(*,'(a)') " done."
      ! unit conversion from mRy to mEv, by division by Rydberg's constant in (eV) (and an factor of 4 for some reason..)
      ! L.B. The factor is fine, a factor of 2 from the def of dynamical matrix and fc2 includes another factor of 2
      ! L.B Temporary fix before replacing 4 with Landeg_glob^2=4
      A = 4.0d0*A/fc2*ry_ev
      B = B/fc2*ry_ev

      ! eigenvalues of A (one set per NA) are written to wres  - AMS
      write (*,'(1x,a)',advance='yes') "Diagonalizing A-matrix for AMS calculation"
      ! The realness of A has been questioned.
      !According to Essenberger et. al PRB 84, 174425 (2011) Eqn. 22 it is correctly so.
      do iq=1,nq
         A(:,:,iq)=real(A(:,:,iq))
      end do
      call eigenvalue_calculation_lapack(A,wres,eigv,na,nq)


      !Take absolute values, in case of AF frequencies
      !Assumption: negative frequencies ok, still gives real energies (not complex frequencies/energies)
      !This is now taken care of in eigenvalue_calculation_lapack() where abs() is applied
      ! Sort eigenvalues and eigenvectors according to magnitude of eigenvalues
      do i=1,nq
         call sortEigenVVs(wres(1:na,i),eigv(1:na,1:na,i),na)
      end do

      ! print energies to output file
      call printEnergies(output_file,wres,msat,tcmfa,tcrpa,1)
      call printEigVects(output_file3,eigv,wres,q,nq,na)
      jqres=0.0d0
      eigv=0.0d0
      call eigenvalue_calculation_lapack(B,jqres,eigv,na,nq)
      call printEnergies(output_file2,jqres,msat,tcmfa,tcrpa,2)
      !
      call magdos_calc(output_file4,wres,magdos)

      !!
      deallocate(wres)
      deallocate(jqres)
      deallocate(A)
      deallocate(eigv)
      deallocate(B)
      write(*,'(a)') " done."
      write(*,'(1x,a)') 'Adiabatic Magnon Spectra Calculations done.'
      10 continue
      ! stop
   end subroutine calculate_ams

   !> Fourier transform of exchange interactions around central atom
   complex(dblprec) function calc_j(mu,nu,q_vect,countstart)
      !
      implicit none
      !
      integer,intent(in) :: mu
      integer,intent(in) :: nu
      integer,intent(in) :: countstart

      real(dblprec),dimension(:),intent(in) :: q_vect
      real(dblprec),dimension(3) :: q_vect2pi
      real(dblprec),dimension(3) :: dist
      integer :: j, mutemp, nutemp
      complex(dblprec) :: i
      real(dblprec) :: mdot

      real(dblprec) :: dmdot
      real(dblprec), dimension(3) :: q_hat

      calc_j=0.0_dblprec
      q_vect2pi=q_vect*2_dblprec*pi
      q_hat=abs(q_vect/sqrt(sum(q_vect*q_vect)+1.0d-15))

      mutemp = mu+countstart
      nutemp = nu+countstart

      i = (0.0d0,1.0d0)

      do j=1,nlistsize(mutemp)
         if (anumb(nutemp)==anumb(nlist(j,mutemp))) then
            dist(:)=coord(1:3,mutemp)-coord(1:3,nlist(j,mutemp))
            mdot=sum(emom(1:3,mutemp,1)*emom(1:3,nlist(j,mutemp),1))
            calc_j = calc_j+ncoup(j,mutemp,gsconf_num)*exp(i* &
               (q_vect2pi(1)*dist(1)+q_vect2pi(2)*dist(2)+q_vect2pi(3)*dist(3)))
         end if
      end do
      if(do_dm==1) then
         do j=1,dmlistsize(mutemp)
            if (anumb(nutemp)==anumb(dmlist(j,mutemp))) then
               dist(:)=coord(1:3,mutemp)-coord(1:3,nlist(j,mutemp))
               dmdot=sum(dm_vect(:,j,mutemp)*q_hat)
               calc_j = calc_j+dmdot*sin(1.0_dblprec*( q_vect2pi(1)*dist(1)+q_vect2pi(2)*dist(2)+q_vect2pi(3)*dist(3)))
            end if
         end do
      end if

   end function calc_j

   !> Fourier transform of anisotropy
   complex(dblprec) function calc_ani(mu,nu,countstart)
      !
      implicit none
      !
      integer,intent(in) :: mu
      integer,intent(in) :: nu
      integer,intent(in) :: countstart


      integer :: mutemp, nutemp

      real(dblprec) :: aedot

      mutemp = mu+countstart
      nutemp = nu+countstart

      calc_ani=0.0_dblprec
      !!!!Add constant shift for uniaxial anisotropy (ugly and wip...)
      if (mutemp==nutemp.and.taniso(nutemp)==1) then
         aedot = sum(eaniso(:,nutemp)*emomm(:,nutemp,1))
         calc_ani = 2.0d0*abs(kaniso(1,nutemp)*aedot**2) + abs(kaniso(2,nutemp)*(aedot**2)**2)
      end if
      return
   end function calc_ani

   !> Calculate eigenvalues and eigenvectors of the dynamical matrix
   subroutine eigenvalue_calculation_lapack(A,wres,eigv,na,nq)
      ! calculating eigenvalues for the matrix (for all k) A(:,:,k)
      implicit none
      integer, intent(in) :: na !< Number of atoms in cell
      integer, intent(in) :: nq !< Number of q-points
      complex(dblprec), dimension(na,na,nq),intent(in) :: A !< A-matrix whose eigenvalues are the sought energies
      complex(dblprec), dimension(na,na,nq),intent(inout) :: eigv !< Eigenvectors from J(q)
      real(dblprec),  dimension(na,nq) :: wres

      integer :: iq, ia
      integer :: lwork, info
      complex(dblprec),dimension(:),allocatable :: work
      complex(dblprec),dimension(:,:),allocatable :: ctemp
      real(dblprec),dimension(:),allocatable :: rwork
      real(dblprec),dimension(:),allocatable :: mineig
      complex(dblprec),dimension(:,:),allocatable :: A_inplace
      complex(dblprec), allocatable, dimension(:) :: cwres
      complex(dblprec), allocatable, dimension(:) :: eig_ave

      lwork=2*na
      allocate(work(lwork))
      allocate(ctemp(na,na))
      allocate(A_inplace(na,na))
      allocate(rwork(2*na))
      allocate(mineig(nq))
      allocate(cwres(NA))
      allocate(eig_ave(NA))

      ! eigenvalue calculations performed, energy = abs(real_part +i*imaginary part)
      do iq = 1,nq
         A_inplace=A(:,:,iq)
         call zgeev('V','V',NA, A(1:NA,1:NA,iq), NA, cwres(1:NA), ctemp, na, eigv(1:NA,1:NA,iq), NA, WORK, LWORK, RWORK, INFO)
         if(info.ne.0) then
            print '(2x,a,i4)', 'Problem in zgeev:',info
         end if
         ! Should we use the absolute value w=sqrt(re(w)**2+im(w)**2) ?
         wres(1:NA,iq)=(real(cwres(1:NA)))
         mineig(iq)=minval(wres(1:NA,iq))
      end do

      if(maxval(mineig)<1.0e-6_dblprec) then
         wres=abs(wres)
      else
         wres=wres-minval(mineig)
         wres=abs(wres)
      end if

      deallocate(work)
      deallocate(rwork)
      deallocate(mineig)
      deallocate(cwres)
      deallocate(ctemp)
      deallocate(A_inplace)
      deallocate(eig_ave)
   end subroutine eigenvalue_calculation_lapack

   !> Eigenvalue and eigenvector calculation
   subroutine eigenvalue_calculation_colpa(A,wres,eigv,na,nq)
      ! calculating eigenvalues for the matrix (for all k) A(:,:,k)
      implicit none
      integer, intent(in) :: na !< Number of atoms in cell
      integer, intent(in) :: nq !< Number of q-points
      complex(dblprec), dimension(na,na,nq),intent(in) :: A !< A-matrix whose eigenvalues are the sought energies
      complex(dblprec), dimension(na,na,nq),intent(inout) :: eigv !< Eigenvectors from J(q)
      real(dblprec),  dimension(na,nq) :: wres

      integer :: iq
      integer :: lwork, info
      complex(dblprec),dimension(:),allocatable :: work
      complex(dblprec),dimension(:,:),allocatable :: ctemp
      real(dblprec),dimension(:),allocatable :: rwork
      complex(dblprec),dimension(:,:),allocatable :: A_inplace
      complex(dblprec), allocatable, dimension(:) :: cwres

      lwork=2*na
      allocate(work(lwork))
      allocate(ctemp(na,na))
      allocate(A_inplace(na,na))
      allocate(rwork(2*na))
      allocate(cwres(NA))

      ! eigenvalue calculations performed, energy = abs(real_part +i*imaginary part)
      do iq = 1,nq
         call zgeev('V','V',NA, A(1:NA,1:NA,iq), NA, cwres(1:NA), ctemp, na, eigv(1:NA,1:NA,iq), NA, WORK, LWORK, RWORK, INFO)
         if(info.ne.0) then
            print '(2x,a,i4)', 'Problem in zgeev:',info
         end if
         ! Saving the real part of eigenvalues
         wres(1:NA,iq)=real(cwres(1:NA))
      end do

      deallocate(work)
      deallocate(rwork)
      deallocate(cwres)
      deallocate(ctemp)
      deallocate(A_inplace)
   end subroutine eigenvalue_calculation_colpa

   !> Print out magnon frequencies and calculate MFA and RPA Tc.
   subroutine printEnergies(filename,wres,msat,tcmfa,tcrpa,flag)
      implicit none

      integer :: i,k,r,flag, ia
      real(dblprec), allocatable, dimension(:,:) :: wres
      real(dblprec) :: tcmfa,tcrpa,msat
      character(LEN=*),intent(in) :: filename

      real(dblprec) :: maxdist, locdist

      open(ofileno,file=filename)

      if(flag==1) then
         tcmfa=0.d0 ; tcrpa=0.d0
      endif
      k=0 ; r=0
      !
      maxdist=get_maxdist(q,nq)
      locdist=0.0_dblprec
      ! before printing the eigenvalues of each matrix is sorted by size
      do i =1,nq
         ! Calculate "q-space distance" to allow for nicer plots
         if(i>1) then
            locdist=locdist+sqrt(sum((q(:,i-1)-q(:,i))**2))
         else
            locdist=0.0_dblprec
         end if

         write(ofileno,1001) i, wres(:,i), locdist/maxdist
         if (flag==1) then
            do ia=1,na
               if (wres(ia,i)>=1.0d-2) then
                  tcrpa=tcrpa+(1.d0/wres(ia,i))
                  r=r+1
               endif
            end do
            tcmfa=tcmfa+sum(wres(:,i))
            k=k+na
         endif
      end do
      if (flag==1) then
         tcmfa=(msat*tcmfa)/(6*k*k_bolt_ev*1000)
         write(*,1002) 'Tc-MFA from AMS:' , tcmfa
         tcrpa=((r*msat)/(6*k_bolt_ev*1000))*(1.d0/tcrpa)
         write(*,1002) 'Tc-RPA from AMS:' , tcrpa
      endif


      ! "safe" if someone calculates the spectra with many atoms in one cell.
      1001 format (2x,i7,2000f18.12)
      1002 format (2x,a,f10.1)
      close(ofileno)
   end subroutine printEnergies

   !> Print out magnon eigenvectors
   subroutine printEigVects(filename,eigv,wres,q,nq,na)
      implicit none

      integer,intent(in) :: nq
      integer,intent(in) :: na
      complex(dblprec),dimension(na,na,nq), intent(inout) :: eigv
      real(dblprec),dimension(na,nq), intent(inout) :: wres
      real(dblprec),dimension(3,nq), intent(in) :: q
      character(LEN=*),intent(in) :: filename

      integer :: i,j,k

      open(ofileno,file=filename)

      ! before printing the eigenvalues of each matrix is sorted by size
      do i =1,nq
         ! Ensuring that the "reference" direction is positive (to avoid fluctuations of theta along q)
         write(ofileno,1001) "# q-point index ",i, " vector: ", q(1:3,i)
         do j=1,na
            ! Compact format : Value only for each complex number
            write(ofileno,1003) wres(j,i),(( real(eigv(k,j,i))),k=1,na)
         end do
      end do

      close(ofileno)

      ! "safe" if someone calculates the spectra with many atoms in one cell.
      1001 format (a,2x,i5,a,3f12.6)
      1002 format (1x,2000f12.6)
      1003 format (1x,f12.6,4x,2000f12.6)
   end subroutine printEigVects

   subroutine sortEigenVVs(array,matrix,nsize)
      ! subroutine to sort the output from cg when calculating eigenvalues,
      implicit none
      integer,intent(in) :: nsize
      real(dblprec),dimension(nsize):: array
      complex(dblprec),dimension(nsize,nsize):: matrix
      integer :: min_index,i,k
      real(dblprec) :: min_value,val_inf
      real(dblprec),dimension(:),allocatable :: sorted_array
      complex(dblprec),dimension(:),allocatable :: temp_row
      complex(dblprec),dimension(:,:),allocatable :: sorted_matrix

      allocate(sorted_array(nsize))
      allocate(temp_row(nsize))
      allocate(sorted_matrix(nsize,nsize))
      sorted_matrix=0.0d0

      val_inf =  1.d10
      min_index = 1
      do k = 1,nsize
         min_value = val_inf
         temp_row=0.0d0
         do i = 1,nsize
            if (array(i)<min_value) then
               min_value = array(i)
               temp_row(1:nsize) = matrix(1:nsize,i)
               min_index = i
            end if
         end do
         array(min_index) = val_inf
         sorted_array(k) = min_value
         sorted_matrix(1:nsize,k) = temp_row(1:nsize)
      end do

      array = sorted_array
      matrix = sorted_matrix

      deallocate(sorted_matrix)
      deallocate(sorted_array)
      deallocate(temp_row)

   end subroutine sortEigenVVs

   !> Sorting of array
   subroutine sortRealArray(array,nsize)
      ! subroutine to sort the output from cg when calculating eigenvalues,
      implicit none
      real(dblprec),dimension(:):: array
      integer,intent(in) :: nsize
      integer :: max_index,i,k
      real(dblprec) :: max_value,minus_inf
      real(dblprec),dimension(:),allocatable :: sorted_array

      allocate(sorted_array(nsize))
      minus_inf = -1.d10
      max_index = 1
      do k = 1,nsize
         max_value = minus_inf
         do i = 1,nsize
            if (array(i)>max_value) then
               max_value = array(i)
               max_index = i
            end if
         end do
         array(max_index) = minus_inf
         sorted_array(k) = max_value
      end do

      array(:) = sorted_array

   end subroutine sortRealArray

   !> Print eigenvectors
   subroutine printEigenvectors()
      ! NOT FINISHED! for printing the eigenvectors
      implicit none
      integer :: j,i,k
      character(LEN=20) :: filename
      filename='ams.eigenvectors.out'


      open(ofileno,file=filename)
      do k=1,nq
         do i = 1,NA
            do j = 1,NA
               write(ofileno,*)
            end do
         end do
      end do
      close(ofileno)

   end subroutine printEigenvectors

   !> print dynamical matrix
   subroutine printA(A)
      ! obsolete, used for troubleshooting
      implicit none
      integer :: i,j
      character(LEN=11) :: nom = 'Amatrix.out'
      complex(dblprec), allocatable, dimension(:,:,:) :: A
      !
      open(ofileno,file=nom)
      do i=1,nq
         do j=1,NA
            write(ofileno,*) A(j,:,i)
         end do
         write(ofileno,*)
      end do
      write(*,*) 'A-matrix printed to Amatrix.out'
      close(ofileno)
   end subroutine printA

   !> Calculate AMS density of states
   subroutine magdos_calc(filename,wres,magdos)
      implicit none

      integer :: i, k, ia, flines
      real(dblprec), allocatable, dimension(:,:) :: wres
      real(dblprec), allocatable,dimension(:,:) :: magdos,mtemp
      real(dblprec) :: emin, emax, deltae, fac, dummy
      character(LEN=*),intent(in) :: filename

      if(do_magdos=='A') then

         ! Open input file
         open(ifileno, file=magdosfile)
         flines=0
         ! Pre-read file to get number of lines (frequencies)
         do
            read(ifileno,*,end=200)  dummy,dummy,dummy
            flines=flines+1
         end do
         200 continue
         rewind(ifileno)
         magdos_freq=flines
         allocate(magdos(magdos_freq,2))
         do k=1,magdos_freq
            read(ifileno,*) magdos(k,1),magdos(k,2)
         enddo
         close(ifileno)
         !
      elseif(do_magdos=='Y') then

         ! Calculate AMS m-DOS
         open(ofileno,file=filename)
         if (magdos_lfreq==0) then
            emin=minval(wres)
         else
            emin=magdos_lfreq
         end if
         if (magdos_hfreq==0) then
            emax=(maxval(wres)+magdos_sigma)
         else
            emax=magdos_hfreq
         endif
         deltae=(emax-emin)/(magdos_freq-1)
         allocate(magdos(magdos_freq,2))
         magdos=0.d0

         do i=1,magdos_freq
            magdos(i,1)=emin+(i-1)*deltae
         enddo

         do k=1,magdos_freq
            do i =1,nq
               do ia=1,na
                  fac=exp(-(wres(ia,i)-magdos(k,1))**2/magdos_sigma)
                  magdos(k,2)=magdos(k,2)+fac
               enddo
            enddo
         enddo

         !Normalize DOS

         magdos(:,2)=magdos(:,2)/(sum(magdos(:,2))*deltae)

         do k=1,magdos_freq
            write(ofileno,1001) magdos(k,1),magdos(k,2),sum(magdos(1:k,2))*deltae
         enddo
      elseif(do_magdos=='N') then
      else    !S/B
         ! Open input file
         open(ifileno, file=magdosfile)
         flines=0
         ! Pre-read file to get number of lines (frequencies)
         do
            read(ifileno,*,end=300)  dummy,dummy,dummy,dummy,dummy
            flines=flines+1
         end do
         300 continue
         rewind(ifileno)
         magdos_freq=flines
         allocate(mtemp(magdos_freq,2))
         allocate(magdos(magdos_freq,2))
         do k=1,magdos_freq
            read(ifileno,*) magdos(k,1),mtemp(k,1),mtemp(k,2)
         enddo
         close(ifileno)

         ! Transversal components (x and y assumed)
         do k=1,magdos_freq
            magdos(k,2)=sqrt(mtemp(k,1)**2+mtemp(k,2)**2)
         enddo
         ! Normalize
         deltae=(magdos(magdos_freq,1)-magdos(1,1))/(magdos_freq-1)
         magdos(:,2)=magdos(:,2)/(sum(magdos(:,2))*deltae)
         deallocate(mtemp)

      endif

      1001 format (2x,f18.12,2f18.12)
      1002 format (2x,a,f10.1)
      close(ofileno)
   end subroutine magdos_calc

   !> q-point distances for nicer plots
   real(dblprec) function get_maxdist(q,nq)
      !
      implicit none
      !
      integer, intent(in) :: nq
      real(dblprec), dimension(3,nq) :: q
      !
      integer :: iq
      real(dblprec) :: dist,maxdist
      !
      dist=0.0d0
      maxdist=0.0d0
      !
      do iq=2,nq
         dist=sqrt(sum((q(:,iq-1)-q(:,iq))**2))
         if (dist>maxdist) maxdist=dist
      end do
      !
      get_maxdist=maxdist
      !
      return
      !
   end function get_maxdist



   !---------------------------------------------------------------------------
   !> @brief
   !> Read input parameters.
   !
   !> @author
   !> Anders Bergman
   !---------------------------------------------------------------------------
   subroutine read_parameters_ams(ifile)
      use FileParser

      implicit none

      ! ... Formal Arguments ...
      integer, intent(in) :: ifile   !< File to read from
      !
      ! ... Local Variables ...
      character(len=50) :: keyword,cache
      integer :: rd_len, i_err, i_errb
      logical :: comment


      do
         10     continue
         ! Read file character for character until first whitespace
         keyword=""
         call bytereader(keyword,rd_len,ifile,i_errb)

         ! converting Capital letters
         call caps2small(keyword)

         ! check for comment markers (currently % and #)
         comment=(scan(trim(keyword),'%')==1).or.(scan(trim(keyword),'#')==1).or.&
            (scan(trim(keyword),'*')==1).or.(scan(trim(keyword),'=')==1.or.&
            (scan(trim(keyword),'!')==1))

         if (comment) then
            read(ifile,*)
         else
            ! Parse keyword
            keyword=trim(keyword)
            select case(keyword)
            !> - simid

         case('do_ams')
            read(ifile,*,iostat=i_err) do_ams
            if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err

         case('do_magdos')
            read(ifile,*,iostat=i_err) do_magdos
            if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err

         case('magdosfile')
            read(ifile,'(a)',iostat=i_err) cache
            if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err
            magdosfile=adjustl(trim(cache))

         case('magdos_sigma')
            read(ifile,*,iostat=i_err) magdos_sigma
            if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err

         case('magdos_freq')
            read(ifile,*,iostat=i_err) magdos_freq
            if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err

         case('magdos_lfreq')
            read(ifile,*,iostat=i_err) magdos_lfreq
            if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err

         case('magdos_hfreq')
            read(ifile,*,iostat=i_err) magdos_hfreq
            if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err

         end select
      end if

      ! End of file
      if (i_errb==20) goto 20
      ! End of row
      if (i_errb==10) goto 10
   end do

   20  continue

   rewind(ifile)
   return
end subroutine read_parameters_ams

subroutine init_ams()
   !
   !
   implicit none
   !
   !Adiabatic Magnon Spectra calculation
   do_ams = 'N'
   do_magdos = 'N'
   magdosfile= 'amsdosfile'
   magdos_sigma=30
   magdos_freq=200
   magdos_hfreq=0
   magdos_lfreq=0
   !
end subroutine init_ams

!--------------------------------------------------------------------------
!
! DESCRIPTION
!> @brief
!> Main driver for adiabatic magnon spectra (AMS) in collinear magnets.
!! Based on linear spin wave theory as described in ....
!---------------------------------------------------------------------------------
subroutine calculate_random_ams()
   implicit none
   integer :: mu
   integer :: nu
   integer :: lambda
   integer :: i,j,k, eidx
   integer :: I1,I2,I3 !< indices to find suitable cell to perform calculations

   integer :: countstart !< index for last atom "before" the cell where calculations are performed

   real(dblprec), dimension(3) :: zeroVector !(0,0,0) vector  to calculate J(0)
   real(dblprec), allocatable, dimension(:,:) :: wres,jqres
   real(dblprec), allocatable, dimension(:,:,:) :: wres_tot
   real(dblprec), allocatable, dimension(:,:) :: ams_hist
   real(dblprec) :: fc,fc2,mtmp, emax

   real(dblprec) :: mag_sign_mu, mag_sign_nu, mag_sign_lambda
   real(dblprec), dimension(3) :: order_vector
   complex(dblprec), allocatable, dimension(:,:,:) :: A !< matrix whose eigenvalues are the sought energies fro AMS
   complex(dblprec), allocatable, dimension(:,:,:) :: eigv !< matrix whose eigenvalues are the sought energies fro AMS
   complex(dblprec), allocatable, dimension(:,:,:) :: B !< matrix whose eigenvalues are the sought energies for J(q)

   ! for output printing
   character(LEN = 16) :: output_file
   character(LEN = 18) :: output_file2
   character(LEN = 18) :: output_file3
   character(LEN = 19) :: output_file4
   integer :: iq, idx

   ! Factors for mRy energy conversion
   fc = mry/mub
   fc2 = 2*mry/mub

   ! if no q allocated then no AMS calculations are performed
   if (.not. allocated(q)) then
      write(*,*) 'No q-file allocated. No Adiabatic Magnon Spectrum calculation is performed'
      goto 10 ! jumps to end of subroutine
   end if

   allocate(A(NA,NA,nq))
   allocate(eigv(NA,NA,nq))
   allocate(B(NA,NA,nq))
   allocate(wres(NA,nq))
   allocate(ams_hist(nq,ams_hist_size))
   allocate(wres_tot(NA,nq,N1*N2*N3))
   wres_tot=0.0d0
   allocate(jqres(NA,nq))
   A = 0.0_dblprec
   B = 0.0_dblprec
   zeroVector = 0.0_dblprec

   output_file4 = 'ra_ams.'//simid//'.out'

   ! since j(q) is depending on positions (scalar product with q-vector), it is desired
   ! to perform calculations on the lattice in the middle of the system, to avoid
   ! problems caused by periodicy. Calculations are done for the lattice starting with
   ! atom number countstart+1
   !

   write (*,'(1x,a)',advance="no") "Mounting A-matrices for random-AMS calculation"
   idx=0
   do I1 = 0,N1-1
      do I2 = 0,N2-1
         do I3 = 0,N3-1

            countstart = 0+I1*NA+I2*N1*NA+I3*N2*N1*NA
            idx=idx+1
            !
            order_vector(1)=1.0d0/sqrt(3.0d0)
            order_vector(2)=1.0d0/sqrt(3.0d0)
            order_vector(3)=1.0d0/sqrt(3.0d0)
            ! A-matrix calculation, one for each q
            A=(0.0d0,0.0d0)
            do i=1,nq
               do mu=1,NA
                  mag_sign_mu=sign(1.0d0,sum(emom(1:3,mu+countstart,1)*order_vector))
                  do nu=1,NA
                     mag_sign_nu=sign(1.0d0,sum(emom(1:3,nu+countstart,1)*order_vector))
                     A(mu,nu,i) = -calc_j(mu,nu,q(:,i),countstart)*mmom(nu+countstart,1)*mag_sign_nu
                     B(mu,nu,i) = -A(mu,nu,i)*mmom(nu+countstart,1)
                     if (mu==nu) then
                        do lambda = 1,NA
                           mag_sign_lambda=sign(1.0d0,sum(emom(1:3,lambda+countstart,1)*order_vector))
                           A(mu,mu,i)=A(mu,mu,i) + calc_j(mu,lambda,zeroVector,countstart)* &
                              mmom(lambda+countstart,1)*mag_sign_lambda
                        end do
                        A(mu,mu,i)=A(mu,mu,i) + calc_ani(mu,mu,countstart)
                        A(mu,mu,i)=A(mu,mu,i) + sum(hfield*emomm(1:3,mu+countstart,1))*mag_sign_nu*0.5d0
                     end if
                     A(mu,nu,i)=A(mu,nu,i)*Landeg(mu)*Landeg(nu) !LB Needs to be fixed,also multiply with landeg_glob
                     B(mu,nu,i)=B(mu,nu,i)*Landeg(mu)*Landeg(nu)
                  end do
               end do
            end do

            ! unit conversion from mRy to mEv, by division by Rydberg's constant in (eV) (and an factor of 4 for some reason..)
            ! L.B. The factor is fine, a factor of 2 from the def of dynamical matrix and fc2 includes another factor of 2
            ! L.B Temporary fix before replacing 4 with Landeg_glob^2=4
            A = 4.0d0*A/fc2*ry_ev !*sign(msat,1.d0)
            B = B/fc2*ry_ev


            ! eigenvalues of A (one set per NA) are written to wres  - AMS
            ! The realness of A has been questioned.
            !According to Essenberger et. al PRB 84, 174425 (2011) Eqn. 22 it is correctly so.
            do iq=1,nq
               A(:,:,iq)=real(A(:,:,iq))
            end do
            !
            call eigenvalue_calculation_lapack(A,wres,eigv,na,nq)
            wres_tot(:,:,idx)=wres(:,:)
         end do
      end do
   end do
   write(*,'(a)') " done."
   !
   emax=1.1d0*maxval(wres_tot)
   ams_hist=0.0d0
   do i=1,nq
      do j=1,na
         do k=1,idx
            eidx=nint(ams_hist_size*wres_tot(i,j,k)/emax)+1
            ams_hist(i,eidx)=ams_hist(i,eidx)+1.0d0
         end do
      end do
   end do
   open(ofileno,file=output_file4)
   do i=1,nq
      do eidx=1,ams_hist_size
         write(ofileno,'(1x,i4,g12.5,2x,g20.8)') i, emax*eidx/ams_hist_size, ams_hist(i,eidx)
      end do
   end do
   close(ofileno)

   !!
   deallocate(wres)
   deallocate(wres_tot)
   deallocate(ams_hist)
   deallocate(jqres)
   deallocate(A)
   deallocate(eigv)
   deallocate(B)
   write(*,'(1x,a)') 'Random-Adiabatic Magnon Spectra Calculations done.'
   10 continue
   ! stop
end subroutine calculate_random_ams


end module AMS
