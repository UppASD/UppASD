!-------------------------------------------------------------------------------
! MODULE: AMS
!> @brief
!> Routines for calculating adiabatic magnon spectra (AMS), and magnon density of states
!
!> @details
!> AMS - Adiabatic Magnon Spectra calculation
!> The adiabatic magnon spectra is written to ams.SIMID.out, and the
!> energies are sorted by nsize on each row.
!
!> Theory taken from Kübler, Theory of Itinerant Electron Magnetism (2009)
!> and Halilov et ad, Phys. Rev. B, 58:293, 1998. However the A-matrix proposed
!> in Kübler has been modified to function with the j-exchange already implemented
!> in UppASD.
!
!> The module contains the main subroutine calculate_ams(), using calc_j(),
!> eigenvalue_calculation() and printEnergies().
!
!> @author
!> A. Bergman, L. Bergqvist, J. Chico, etc
!> @copyright
!> GNU Public License.
!-------------------------------------------------------------------------------
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
   use InputDataType
   !use InputData,          only : N1,N2,N2,N3,NA,NT,Natom,simid,do_dm,hfield,gsconf_num,do_anisotropy,nchmax,ammom_inp,do_ralloy
   !!! use InputData,          only : N1,N2,N2,N3,NA,NT,Natom,simid,hfield,gsconf_num,nchmax,ammom_inp,do_ralloy, ham_inp
   use InputData, only : ham_inp
   use Qvectors,        only : q,nq,q_weight
   use Hamiltoniandata,    only : ham
   !use Systemdata,         only : anumb,coord, Landeg, atype
   use Systemdata,         only : coord, Landeg
   use Momentdata,         only : mmom, emom, emomM
   use ChemicalData,       only : achem_ch,asite_ch
   use Profiling
   use Math_functions, only : f_wrap_coord_diff

   implicit none


   ! Adiabatic Magnon Spectra calculation flag
   character(LEN=1) :: do_ams                  !< Calculate AMS (N/Y)
   character(LEN=1) :: do_magdos               !< Calculate AMS DOS (N/Y/F)
   character(len=35) :: magdosfile             !< AMS DOS file
   real(dblprec)    :: magdos_sigma            !< Frequency broadening of AMS DOS (in meV)
   integer          :: magdos_freq             !< Number of frequencies of AMS DOS
   integer          :: magdos_lfreq             !< Low cut-off frequency of AMS DOS (in meV)
   integer          :: magdos_hfreq             !< High cut-off frequency of AMS DOS (in meV)
   integer          :: magdos_rasamples        !< Number of samples for random alloy AMS

   integer          :: ams_hist_size=2000

   real(dblprec), allocatable, dimension(:,:) :: magdos
   real(dblprec) :: tcmfa,tcrpa,msat

   private

   public :: magdos,tcmfa,tcrpa,msat, magdos_freq, do_ams, do_magdos
   public :: magdosfile, magdos_sigma
   public :: calculate_ams,calc_j,eigenvalue_calculation_lapack,printEnergies
   public :: magdos_calc, read_parameters_ams, init_ams
   public :: calculate_random_ams

contains

   !--------------------------------------------------------------------------
   ! SUBROUTINE: calculate_ams
   ! DESCRIPTION
   !> @brief
   !> Main driver for adiabatic magnon spectra (AMS) in collinear magnets.
   !> Based on linear spin wave theory as described in ....
   !---------------------------------------------------------------------------------
   !use InputData,          only : N1,N2,N2,N3,NA,NT,Natom,simid,hfield,gsconf_num,nchmax,ammom_inp,do_ralloy, ham_inp
   subroutine calculate_ams(N1,N2,N3,NA,NT,Natom,anumb,NA_list,simid,hfield,do_ralloy)
      implicit none

      integer, intent(in) :: N1  !< Number of cell repetitions in x direction
      integer, intent(in) :: N2  !< Number of cell repetitions in y direction
      integer, intent(in) :: N3  !< Number of cell repetitions in z direction
      integer, intent(in) :: NA  !< Number of atoms in one cell
      integer, intent(in) :: NT  !< Number of types of atoms
      integer, intent(in) :: Natom     !< Number of atoms in system
      integer, dimension(Natom), intent(inout)   :: anumb       !< Atom number in cell
      integer, dimension(Natom), intent(inout)   :: NA_list     !< Reverse anumb list 
      character(len=8), intent(in) :: simid !< Name of simulation
      real(dblprec), dimension(3), intent(in) :: hfield !< Constant effective field
      integer, intent(in), optional :: do_ralloy     !< random alloy simulation (0/1)

      integer :: mu
      integer :: nu
      integer :: lambda
      integer :: i, naref
      integer :: I1,I2,I3 !< indices to find suitable cell to perform calculations
      integer :: mua, nua, lambdaa

      integer :: countstart !< index for last atom "before" the cell where calculations are performed

      real(dblprec), dimension(3) :: zeroVector !(0,0,0) vector  to calculate J(0)
      real(dblprec), allocatable, dimension(:,:) :: wres,jqres
      real(dblprec) :: fc, fc2

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
      integer :: i_stat, i_all

      ! Factors for mRy energy conversion
      fc = mry/mub
      fc2 = 2*mry/mub

      ! if no q allocated then no AMS calculations are performed
      if (.not. allocated(q)) then
         write(*,*) 'No q-file allocated. No Adiabatic Magnon Spectrum calculation is performed'
         goto 10 ! jumps to end of subroutine
      end if

      if (.not. allocated(q_weight)) then
         allocate(q_weight(nq),stat=i_stat)
         call memocc(i_stat,product(shape(q_weight))*kind(q_weight),'q_weight','calculate_ams')
         q_weight=1.0_dblprec
      end if

      allocate(A(NA,NA,nq),stat=i_stat)
      call memocc(i_stat,product(shape(A))*kind(A),'A','calculate_ams')
      allocate(eigv(NA,NA,nq),stat=i_stat)
      call memocc(i_stat,product(shape(eigv))*kind(eigv),'eigv','calculate_ams')
      allocate(B(NA,NA,nq),stat=i_stat)
      call memocc(i_stat,product(shape(B))*kind(B),'B','calculate_ams')
      allocate(wres(NA,nq),stat=i_stat)
      call memocc(i_stat,product(shape(wres))*kind(wres),'wres','calculate_ams')
      allocate(jqres(NA,nq),stat=i_stat)
      call memocc(i_stat,product(shape(jqres))*kind(jqres),'jqres','calculate_ams')
      
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

      naref=Natom/N1/N2/N3
      countstart = 0+I1*NAref+I2*N1*NAref+I3*N2*N1*NAref-1
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
      msat=sum(mmom(:,1),1)/natom

      order_vector(1)=1.0_dblprec/sqrt(3.0_dblprec)
      order_vector(2)=1.0_dblprec/sqrt(3.0_dblprec)
      order_vector(3)=1.0_dblprec/sqrt(3.0_dblprec)

      A=(0.0_dblprec,0.0_dblprec)
      B=(0.0_dblprec,0.0_dblprec)
      ! A-matrix calculation, one for each q
      write (*,'(1x,a)',advance="no") "Mounting A-matrix for AMS calculation"
      !$omp parallel do default(shared),private(i,mu,nu,lambda,mag_sign_mu,mag_sign_nu,mag_sign_lambda,mua,nua,lambdaa)
      do i=1,nq
         do mu=1,NA
            mua=NA_list(mu)
            mag_sign_mu=sign(1.0_dblprec,sum(emom(1:3,mua+countstart,1)*order_vector))
            do nu=1,NA
               nua=NA_list(nu)
               mag_sign_nu=sign(1.0_dblprec,sum(emom(1:3,nua+countstart,1)*order_vector))
               A(mu,nu,i) = -calc_j(mua,nua,q(:,i),countstart,anumb)*mmom(nua+countstart,1)*mag_sign_nu  !*mag_sign_mu
               B(mu,nu,i) = -A(mu,nu,i)*mmom(nua+countstart,1)
               if (mu==nu) then
                  do lambda = 1,NA
                     lambdaa = NA_list(lambda)
                     mag_sign_lambda=sign(1.0_dblprec,sum(emom(1:3,lambdaa+countstart,1)*order_vector))
                     A(mu,mu,i)=A(mu,mu,i) + calc_j(mua,lambdaa,zeroVector,countstart,anumb)* &
                        mmom(lambdaa+countstart,1)*mag_sign_lambda
                  end do
                  A(mu,mu,i)=A(mu,mu,i) + calc_ani(mua,mua,countstart,ham_inp%do_anisotropy)
                  A(mu,mu,i)=A(mu,mu,i) + sum(hfield*emomm(1:3,mua+countstart,1))*mag_sign_nu*0.50_dblprec
               end if
               A(mu,nu,i)=A(mu,nu,i)*Landeg(mua)*Landeg(nua) !LB Needs to be fixed,also multiply with landeg_glob
               B(mu,nu,i)=B(mu,nu,i)*Landeg(mua)*Landeg(nua)
            end do
         end do
      enddo
      !$omp end parallel do
      write(*,'(a)') " done."
      ! unit conversion from mRy to mEv, by division by Rydberg's constant in (eV) (and an factor of 4 for some reason..)
      ! L.B. The factor is fine, a factor of 2 from the def of dynamical matrix and fc2 includes another factor of 2
      ! L.B Temporary fix before replacing 4 with Landeg_glob^2=4
      A = 4.0_dblprec*A/fc2*ry_ev !*sign(msat,1._dblprec)
      B = B/fc2*ry_ev
      ! eigenvalues of A (one set per NA) are written to wres  - AMS
      ! write (*,'(1x,a)',advance='yes') "Diagonalizing A-matrix for AMS calculation"

      call eigenvalue_calculation_lapack(A,B,wres,jqres,eigv,na,nq)

      !Take absolute values, in case of AF frequencies
      !Assumption: negative frequencies ok, still gives real energies (not complex frequencies/energies)
      !This is now taken care of in eigenvalue_calculation_lapack() where abs() is applied
      ! Sort eigenvalues and eigenvectors according to magnitude of eigenvalues
      do i=1,nq
         call sortEigenVVs(wres(1:na,i),eigv(1:na,1:na,i),na)
      end do

      ! print energies to output file
      call printEnergies(output_file,wres,msat,tcmfa,tcrpa,na,1)
      call printEigVects(output_file3,eigv,wres,q,nq,na)
      !!

      call printEnergies(output_file2,jqres,msat,tcmfa,tcrpa,na,2)
      !
      call magdos_calc(output_file4,wres,na,nq)
      !call magdos_calc(output_file4,wres,magdos,na,nq)

      !!
      i_all=-product(shape(wres))*kind(wres)
      deallocate(wres,stat=i_stat)
      call memocc(i_stat,i_all,'wres','calculate_ams')
      i_all=-product(shape(jqres))*kind(jqres)
      deallocate(jqres,stat=i_stat)
      call memocc(i_stat,i_all,'jqres','calculate_ams')
      i_all=-product(shape(A))*kind(A)
      deallocate(A,stat=i_stat)
      call memocc(i_stat,i_all,'A','calculate_ams')
      i_all=-product(shape(eigv))*kind(eigv)
      deallocate(eigv,stat=i_stat)
      call memocc(i_stat,i_all,'eigv','calculate_ams')
      i_all=-product(shape(B))*kind(B)
      deallocate(B,stat=i_stat)
      call memocc(i_stat,i_all,'B','calculate_ams')
!      write(*,'(a)') " done."
      write(*,'(1x,a)') 'Adiabatic Magnon Spectra Calculations done.'
      10 continue
   end subroutine calculate_ams

   !-----------------------------------------------------------------------------
   ! FUNCTION: calc_j
   !> Fourier transform of exchange interactions around central atom
   !-----------------------------------------------------------------------------
   complex(dblprec) function calc_j(mu,nu,q_vect,countstart,anumb)
      !
      use InputData, only : gsconf_num, Natom

      implicit none
      !
      integer,intent(in) :: mu
      integer,intent(in) :: nu
      integer,intent(in) :: countstart
      integer, dimension(Natom), intent(in)   :: anumb       !< Atom number in cell

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

      i = (0.0_dblprec,1.0_dblprec)

      do j=1,ham%nlistsize(ham%aham(mutemp))
         if (anumb(nutemp)==anumb(ham%nlist(j,mutemp))) then
            !dist(:)=-redcoord(atype(mutemp),j,:)
            call f_wrap_coord_diff(Natom,coord,mutemp,ham%nlist(j,mutemp),dist)
            !dist(:)=coord(1:3,mutemp)-coord(1:3,ham%nlist(j,mutemp))
            mdot=sum(emom(1:3,mutemp,1)*emom(1:3,ham%nlist(j,mutemp),1))
            calc_j = calc_j+ham%ncoup(j,ham%aham(mutemp),gsconf_num)*exp(i* &
               (q_vect2pi(1)*dist(1)+q_vect2pi(2)*dist(2)+q_vect2pi(3)*dist(3)))
         end if
      end do
      if(ham_inp%do_dm==1) then
         do j=1,ham%dmlistsize(ham%aham(mutemp))
            if (anumb(nutemp)==anumb(ham%dmlist(j,mutemp))) then
               !dist(:)=-redcoord(atype(mutemp),j,:)
               call f_wrap_coord_diff(Natom,coord,mutemp,ham%nlist(j,mutemp),dist)
               !dist(:)=coord(1:3,mutemp)-coord(1:3,ham%nlist(j,mutemp))
               dmdot=sum(ham%dm_vect(:,j,ham%aham(mutemp))*q_hat)
               calc_j = calc_j+dmdot*sin(1.0_dblprec*( q_vect2pi(1)*dist(1)+q_vect2pi(2)*dist(2)+q_vect2pi(3)*dist(3)))
            end if
         end do
      end if

   end function calc_j

   !-----------------------------------------------------------------------------
   ! FUNCTION: calc_jRA
   !> Fourier transform of exchange interactions around central atom of random alloy
   !-----------------------------------------------------------------------------
   complex(dblprec) function calc_jRA(mu,nu,q_vect,iatom)
      !
      use InputData, only : gsconf_num, Natom
      !
      implicit none
      !
      integer,intent(in) :: mu
      integer,intent(in) :: nu
      integer,intent(in) :: iatom

      real(dblprec),dimension(:),intent(in) :: q_vect
      real(dblprec),dimension(3) :: q_vect2pi
      real(dblprec),dimension(3) :: dist
      integer :: j
      complex(dblprec) :: i
      real(dblprec) :: mdot

      real(dblprec) :: dmdot
      real(dblprec), dimension(3) :: q_hat

      calc_jRA=0.0_dblprec
      q_vect2pi=q_vect*2_dblprec*pi
      q_hat=abs(q_vect/sqrt(sum(q_vect*q_vect)+1.0d-15))

      i = (0.0_dblprec,1.0_dblprec)

      do j=1,ham%nlistsize(ham%aham(iatom))
         if (achem_ch(ham%nlist(j,iatom))==nu) then
            call f_wrap_coord_diff(Natom,coord,iatom,ham%nlist(j,iatom),dist)
            mdot=sum(emom(1:3,iatom,1)*emom(1:3,ham%nlist(j,iatom),1))
            calc_jRA = calc_jRA+ham%ncoup(j,ham%aham(iatom),gsconf_num)*exp(i* &
               (q_vect2pi(1)*dist(1)+q_vect2pi(2)*dist(2)+q_vect2pi(3)*dist(3)))
         end if
      end do
      if(ham_inp%do_dm==1) then
         do j=1,ham%dmlistsize(ham%aham(iatom))
            if (achem_ch(ham%dmlist(j,iatom))==nu) then
               call f_wrap_coord_diff(Natom,coord,iatom,ham%nlist(j,iatom),dist)
               dmdot=sum(ham%dm_vect(:,j,ham%aham(iatom))*q_hat)
               calc_jRA = calc_jRA+dmdot*sin(1.0_dblprec*( q_vect2pi(1)*dist(1)+q_vect2pi(2)*dist(2)+q_vect2pi(3)*dist(3)))
            end if
         end do
      end if

   end function calc_jRA

   !-----------------------------------------------------------------------------
   ! FUNCTION: calc_jDRA
   !> Fourier transform of exchange interactions around central atom of random alloy
   !-----------------------------------------------------------------------------
   complex(dblprec) function calc_jDRA(alfa,beta,q_vect,iatom)
      !
      use InputData, only : gsconf_num, Natom, ammom_inp
      implicit none
      !
      integer,intent(in) :: alfa
      integer,intent(in) :: beta
      integer,intent(in) :: iatom

      real(dblprec),dimension(:),intent(in) :: q_vect
      real(dblprec),dimension(3) :: q_vect2pi
      real(dblprec),dimension(3) :: dist
      integer :: j
      complex(dblprec) :: i
      real(dblprec) :: mdot

      real(dblprec) :: dmdot
      real(dblprec), dimension(3) :: q_hat

      calc_jDRA=0.0_dblprec
      q_vect2pi=q_vect*2_dblprec*pi
      q_hat=abs(q_vect/sqrt(sum(q_vect*q_vect)+1.0d-15))

      i = (0.0_dblprec,1.0_dblprec)

      do j=1,ham%nlistsize(ham%aham(iatom))
         if (asite_ch(ham%nlist(j,iatom))==beta) then
            call f_wrap_coord_diff(Natom,coord,iatom,ham%nlist(j,iatom),dist)
            mdot=sum(emom(1:3,iatom,1)*emom(1:3,ham%nlist(j,iatom),1))
            calc_jDRA = calc_jDRA+ham%ncoup(j,ham%aham(iatom),gsconf_num)*exp(i* &
               (q_vect2pi(1)*dist(1)+q_vect2pi(2)*dist(2)+q_vect2pi(3)*dist(3))) * &
               ammom_inp(asite_ch(ham%nlist(j,iatom)),achem_ch(ham%nlist(j,iatom)),gsconf_num)
         end if
      end do
      if(ham_inp%do_dm==1) then
         do j=1,ham%dmlistsize(ham%aham(iatom))
            if (asite_ch(ham%dmlist(j,iatom))==beta) then
               call f_wrap_coord_diff(Natom,coord,iatom,ham%nlist(j,iatom),dist)
               dmdot=sum(ham%dm_vect(:,j,ham%aham(iatom))*q_hat)
               calc_jDRA = calc_jDRA+dmdot*sin(1.0_dblprec*( q_vect2pi(1)*dist(1)+q_vect2pi(2)*dist(2)+q_vect2pi(3)*dist(3)))* &
                ammom_inp(asite_ch(ham%nlist(j,iatom)),achem_ch(ham%nlist(j,iatom)),gsconf_num)
            end if
         end do
      end if

   end function calc_jDRA
   !-----------------------------------------------------------------------------
   ! FUNCTION: calc_ani
   !> Fourier transform of anisotropy
   !-----------------------------------------------------------------------------
   complex(dblprec) function calc_ani(mu,nu,countstart,do_anisotropy)
      !
      implicit none
      !
      integer,intent(in) :: mu
      integer,intent(in) :: nu
      integer,intent(in) :: countstart
      integer,intent(in) :: do_anisotropy

      integer :: mutemp, nutemp

      real(dblprec) :: aedot

      mutemp = mu+countstart
      nutemp = nu+countstart
      !!
      calc_ani=0.0_dblprec
      !!!!Add constant shift for uniaxial anisotropy (ugly and wip...)
      if (do_anisotropy==1) then
         if (mutemp==nutemp.and.ham%taniso(nutemp)==1) then
            aedot = sum(ham%eaniso(:,nutemp)*emomm(:,nutemp,1))
            calc_ani = 2.0_dblprec*abs(ham%kaniso(1,nutemp)*aedot**2) + abs(ham%kaniso(2,nutemp)*(aedot**2)**2)
         end if
      endif
      return
   end function calc_ani

   !-----------------------------------------------------------------------------
   ! SUBROUTINE: eigenvalue_calculation_lapack
   !> @brief Calculate eigenvalues and eigenvectors of the dynamical matrix
   !-----------------------------------------------------------------------------
   subroutine eigenvalue_calculation_lapack(A,B,wres,jqres,eigv,na,nq)
      ! calculating eigenvalues for the matrix (for all k) A(:,:,k)
      implicit none
      integer, intent(in) :: na !< Number of atoms in cell
      integer, intent(in) :: nq !< Number of q-points
      complex(dblprec), dimension(na,na,nq),intent(in) :: A !< A-matrix whose eigenvalues are the sought energies
      complex(dblprec), dimension(na,na,nq),intent(in) :: B !< A-matrix whose eigenvalues are the sought energies
      complex(dblprec), dimension(na,na,nq),intent(inout) :: eigv !< Eigenvectors from J(q)
      real(dblprec),  dimension(na,nq), intent(out) :: wres
      real(dblprec),  dimension(na,nq), intent(out) :: jqres

      integer :: iq, i_stat,i_all
      integer :: lwork, info
      complex(dblprec),dimension(:),allocatable :: work
      complex(dblprec),dimension(:,:),allocatable :: ctemp
      real(dblprec),dimension(:),allocatable :: rwork
      real(dblprec),dimension(:),allocatable :: mineig
      !complex(dblprec),dimension(:,:),allocatable :: A_inplace
      complex(dblprec), allocatable, dimension(:) :: cwres
      complex(dblprec), allocatable, dimension(:) :: eig_ave
      complex(dblprec), dimension(na,na,nq) :: eigq !< Eigenvectors from J(q)

      lwork=2*na
      allocate(work(lwork),stat=i_stat)
      call memocc(i_stat,product(shape(work))*kind(work),'work','eigenvalue_calculation_lapack')
      allocate(ctemp(na,na),stat=i_stat)
      call memocc(i_stat,product(shape(ctemp))*kind(ctemp),'ctemp','eigenvalue_calculation_lapack')
      !allocate(A_inplace(na,na),stat=i_stat)
      !call memocc(i_stat,product(shape(A_inplace))*kind(A_inplace),'A_inplace','eigenvalue_calculation_lapack')
      allocate(rwork(2*na),stat=i_stat)
      call memocc(i_stat,product(shape(rwork))*kind(rwork),'rwork','eigenvalue_calculation_lapack')
      allocate(mineig(nq),stat=i_stat)
      call memocc(i_stat,product(shape(mineig))*kind(mineig),'mineig','eigenvalue_calculation_lapack')
      allocate(cwres(NA),stat=i_stat)
      call memocc(i_stat,product(shape(cwres))*kind(cwres),'cwres','eigenvalue_calculation_lapack')
      allocate(eig_ave(NA),stat=i_stat)
      call memocc(i_stat,product(shape(eig_ave))*kind(eig_ave),'eig_ave','eigenvalue_calculation_lapack')

      ! eigenvalue calculations performed, energy = abs(real_part +i*imaginary part)
      do iq = 1,nq
         !A_inplace=A(:,:,iq)
         call zgeev('V','V',NA, A(1:NA,1:NA,iq), NA, cwres(1:NA), ctemp, na, eigv(1:NA,1:NA,iq), NA, WORK, LWORK, RWORK, INFO)
         if(info.ne.0) then
            print '(2x,a,i4)', 'Problem in zgeev:',info
         end if
         wres(1:NA,iq)=abs(real(cwres(1:NA)))
         mineig(iq)=minval(wres(1:NA,iq))
      end do
      wres=abs(wres)

      ! eigenvalue calculations performed, energy = abs(real_part +i*imaginary part)
      do iq = 1,nq
         !A_inplace=B(:,:,iq)
         call zgeev('V','V',NA, B(1:NA,1:NA,iq), NA, cwres(1:NA), ctemp, na, eigq(1:NA,1:NA,iq), NA, WORK, LWORK, RWORK, INFO)
         if(info.ne.0) then
            print '(2x,a,i4)', 'Problem in zgeev:',info
         end if
         jqres(1:NA,iq)=real(cwres(1:NA))
      end do

      !!!if(maxval(mineig)<1.0e-6_dblprec) then
      !!!   wres=abs(wres)
      !!!else
      !!!   wres=wres-minval(mineig)
      !!!   wres=abs(wres)
      !!!end if
      !!!print *,maxval(wres(ia,:)),minval(wres(ia,:))

      !!!do ia=1,NA
      !!!   if(maxval(wres(ia,:))>0.0_dblprec.and.minval(wres(ia,:))<0.0_dblprec) then
      !!!      wres=wres-minval(wres(ia,:))
      !!!   end if
      !!!end do

      i_all=-product(shape(work))*kind(work)
      deallocate(work,stat=i_stat)
      call memocc(i_stat,i_all,'work','eigenvalue_calculation_lapack')
      i_all=-product(shape(rwork))*kind(rwork)
      deallocate(rwork,stat=i_stat)
      call memocc(i_stat,i_all,'rwork','eigenvalue_calculation_lapack')
      i_all=-product(shape(mineig))*kind(mineig)
      deallocate(mineig,stat=i_stat)
      call memocc(i_stat,i_all,'mineig','eigenvalue_calculation_lapack')
      i_all=-product(shape(cwres))*kind(cwres)
      deallocate(cwres,stat=i_stat)
      call memocc(i_stat,i_all,'cwres','eigenvalue_calculation_lapack')
      i_all=-product(shape(ctemp))*kind(ctemp)
      deallocate(ctemp,stat=i_stat)
      call memocc(i_stat,i_all,'ctemp','eigenvalue_calculation_lapack')
      !i_all=-product(shape(A_inplace))*kind(A_inplace)
      !deallocate(A_inplace,stat=i_stat)
      !call memocc(i_stat,i_all,'A_inplace','eigenvalue_calculation_lapack')
      i_all=-product(shape(eig_ave))*kind(eig_ave)
      deallocate(eig_ave,stat=i_stat)
      call memocc(i_stat,i_all,'eig_ave','eigenvalue_calculation_lapack')

   end subroutine eigenvalue_calculation_lapack


   !-----------------------------------------------------------------------------
   ! SUBROUTINE: eigenvalue_calculation_colpa
   !> @brief Eigenvalue and eigenvector calculation
   !-----------------------------------------------------------------------------
   subroutine eigenvalue_calculation_colpa(A,wres,eigv,NA,nq)
      ! calculating eigenvalues for the matrix (for all k) A(:,:,k)
      implicit none
      integer, intent(in) :: na !< Number of atoms in cell
      integer, intent(in) :: nq !< Number of q-points
      complex(dblprec), dimension(NA,NA,nq),intent(in) :: A !< A-matrix whose eigenvalues are the sought energies
      complex(dblprec), dimension(NA,NA,nq),intent(inout) :: eigv !< Eigenvectors from J(q)
      real(dblprec),  dimension(NA,nq) :: wres

      integer :: iq,i_stat,i_all
      integer :: lwork, info
      complex(dblprec),dimension(:),allocatable :: work
      complex(dblprec),dimension(:,:),allocatable :: ctemp
      real(dblprec),dimension(:),allocatable :: rwork
      complex(dblprec),dimension(:,:),allocatable :: A_inplace
      complex(dblprec), allocatable, dimension(:) :: cwres

      lwork=2*na
      allocate(work(lwork),stat=i_stat)
      call memocc(i_stat,product(shape(work))*kind(work),'work','eigenvalue_calculation_colpa')
      allocate(ctemp(NA,NA),stat=i_stat)
      call memocc(i_stat,product(shape(ctemp))*kind(ctemp),'ctemp','eigenvalue_calculation_colpa')
      allocate(A_inplace(NA,NA),stat=i_stat)
      call memocc(i_stat,product(shape(A_inplace))*kind(A_inplace),'A_inplace','eigenvalue_calculation_colpa')
      allocate(rwork(2*NA),stat=i_stat)
      call memocc(i_stat,product(shape(rwork))*kind(rwork),'rwork','eigenvalue_calculation_colpa')
      allocate(cwres(NA),stat=i_stat)
      call memocc(i_stat,product(shape(cwres))*kind(cwres),'cwres','eigenvalue_calculation_colpa')

      ! eigenvalue calculations performed, energy = abs(real_part +i*imaginary part)
      do iq = 1,nq
         call zgeev('V','V',NA, A(1:NA,1:NA,iq), NA, cwres(1:NA), ctemp, na, eigv(1:NA,1:NA,iq), NA, WORK, LWORK, RWORK, INFO)
         if(info.ne.0) then
            print '(2x,a,i4)', 'Problem in zgeev:',info
         end if
         ! Saving the real part of eigenvalues
         wres(1:NA,iq)=real(cwres(1:NA))
         ! Should we use the absolute value w=sqrt(re(w)**2+im(w)**2) ?
      end do

      i_all=-product(shape(work))*kind(work)
      deallocate(work,stat=i_stat)
      call memocc(i_stat,i_all,'work','eigenvalue_calculation_colpa')
      i_all=-product(shape(rwork))*kind(rwork)
      deallocate(rwork,stat=i_stat)
      call memocc(i_stat,i_all,'rwork','eigenvalue_calculation_colpa')
      i_all=-product(shape(cwres))*kind(cwres)
      deallocate(cwres,stat=i_stat)
      call memocc(i_stat,i_all,'cwres','eigenvalue_calculation_colpa')
      i_all=-product(shape(ctemp))*kind(ctemp)
      deallocate(ctemp,stat=i_stat)
      call memocc(i_stat,i_all,'ctemp','eigenvalue_calculation_colpa')
      i_all=-product(shape(A_inplace))*kind(A_inplace)
      deallocate(A_inplace,stat=i_stat)
      call memocc(i_stat,i_all,'A_inplace','eigenvalue_calculation_colpa')

   end subroutine eigenvalue_calculation_colpa

   !-----------------------------------------------------------------------------
   ! SUBROUTINE: printEnergies
   !> @brief Print out magnon frequencies and calculate MFA and RPA Tc.
   !-----------------------------------------------------------------------------
   subroutine printEnergies(filename,wres,msat,tcmfa,tcrpa,iat,flag)
      implicit none

      integer :: i,k,r,flag, ia,iat
      real(dblprec),intent(in), dimension(:,:) :: wres
      real(dblprec),intent(out) :: tcmfa,tcrpa
      real(dblprec), intent(in) :: msat
      character(LEN=*),intent(in) :: filename

      real(dblprec) :: maxdist, locdist

      open(ofileno,file=filename)

      if(flag==1.or.flag==3) then
         tcmfa=0.0_dblprec ; tcrpa=0.0_dblprec
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
         if (flag==1.or.flag==3) then
            do ia=1,iat
               if (wres(ia,i)>=1.0d-2) then
                  tcrpa=tcrpa+(1.0_dblprec/wres(ia,i))*q_weight(i)
                  r=r+q_weight(i)
               endif
               tcmfa=tcmfa+wres(ia,i)*q_weight(i)
               k=k+q_weight(i)
            end do
         endif
      end do
      if (flag==1) then
         tcmfa=(msat*tcmfa)/(6*k*k_bolt_ev*1000)
         write(*,1002) 'Tc-MFA from AMS:' , tcmfa
         tcrpa=((r*msat)/(6*k_bolt_ev*1000))*(1.0_dblprec/tcrpa)
         write(*,1002) 'Tc-RPA from AMS:' , tcrpa
      elseif (flag==3) then
         tcmfa=(msat*tcmfa)/(6*k*k_bolt_ev*1000)
         write(*,1002) 'Tc-MFA from nc-AMS:' , tcmfa
         tcrpa=((r*msat)/(6*k_bolt_ev*1000))*(1.0_dblprec/tcrpa)
         write(*,1002) 'Tc-RPA from nc-AMS:' , tcrpa
      endif

      ! "safe" if someone calculates the spectra with many atoms in one cell.
      1001 format (2x,i7,2000f18.12)
      1002 format (2x,a,f10.1)
      close(ofileno)
   end subroutine printEnergies

   !-----------------------------------------------------------------------------
   ! SUBROUTINE: printEigVects
   !> @brief Print out magnon eigenvectors
   !-----------------------------------------------------------------------------
   subroutine printEigVects(filename,eigv,wres,q,nq,NA)

      implicit none

      integer,intent(in) :: nq
      integer,intent(in) :: NA
      complex(dblprec),dimension(NA,NA,nq), intent(inout) :: eigv
      real(dblprec),dimension(NA,nq), intent(inout) :: wres
      real(dblprec),dimension(3,nq), intent(in) :: q
      character(LEN=*),intent(in) :: filename

      integer :: i,j,k

      open(ofileno,file=filename)

      ! before printing the eigenvalues of each matrix is sorted by size
      do i =1,nq
         ! Ensuring that the "reference" direction is positive (to avoid fluctuations of theta along q)
         write(ofileno,1001) "# q-point index ",i, " vector: ", q(1:3,i)
         do j=1,NA
            ! Compact format : Value only for each complex number
            write(ofileno,1003) wres(j,i),(( real(eigv(k,j,i))),k=1,na)
         end do
      end do

      close(ofileno)

      ! "safe" if someone calculates the spectra with many atoms in one cell.
      1001 format (a,2x,i5,a,3f12.6)
      1003 format (1x,f12.6,4x,2000f12.6)
   end subroutine printEigVects

   !-----------------------------------------------------------------------------
   ! SUBROUTINE: sortEigenVVs
   !> @brief Subroutine to sort the output from cg when calculating eigenvalues,
   !-----------------------------------------------------------------------------
   subroutine sortEigenVVs(array,matrix,nsize)
      ! subroutine to sort the output from cg when calculating eigenvalues,
      implicit none
      integer,intent(in) :: nsize
      real(dblprec),dimension(nsize):: array
      complex(dblprec),dimension(nsize,nsize):: matrix
      integer :: min_index,i,k,i_stat,i_all
      real(dblprec) :: min_value,val_inf
      real(dblprec),dimension(:),allocatable :: sorted_array
      complex(dblprec),dimension(:),allocatable :: temp_row
      complex(dblprec),dimension(:,:),allocatable :: sorted_matrix

      allocate(sorted_array(nsize),stat=i_stat)
      !call memocc(i_stat,product(shape(sorted_array))*kind(sorted_array),'sorted_array','sortEigenVVs')
      allocate(temp_row(nsize),stat=i_stat)
      !call memocc(i_stat,product(shape(temp_row))*kind(temp_row),'temp_row','sortEigenVVs')
      allocate(sorted_matrix(nsize,nsize),stat=i_stat)
      !call memocc(i_stat,product(shape(sorted_matrix))*kind(sorted_matrix),'sorted_matrix','sortEigenVVs')
      sorted_matrix=0.0_dblprec

      val_inf =  1.d10
      min_index = 1
      do k = 1,nsize
         min_value = val_inf
         temp_row=0.0_dblprec
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

      !i_all=-product(shape(sorted_matrix))*kind(sorted_matrix)
      deallocate(sorted_matrix,stat=i_stat)
      !call memocc(i_stat,i_all,'sorted_matrix','sortEigenVVs')
      !i_all=-product(shape(sorted_array))*kind(sorted_array)
      deallocate(sorted_array,stat=i_stat)
      !call memocc(i_stat,i_all,'sorted_array','sortEigenVVs')
      !i_all=-product(shape(temp_row))*kind(temp_row)
      deallocate(temp_row,stat=i_stat)
      !call memocc(i_stat,i_all,'temp_row','sortEigenVVs')

   end subroutine sortEigenVVs

   !-----------------------------------------------------------------------------
   ! SUBROUTINE: sortRealArray
   !> @brief Sorting of array
   !-----------------------------------------------------------------------------
   subroutine sortRealArray(array,nsize)
      ! subroutine to sort the output from cg when calculating eigenvalues,
      implicit none
      real(dblprec),dimension(:):: array
      integer,intent(in) :: nsize
      integer :: max_index,i,k,i_stat,i_all
      real(dblprec) :: max_value,minus_inf
      real(dblprec),dimension(:),allocatable :: sorted_array

      allocate(sorted_array(nsize),stat=i_stat)
      call memocc(i_stat,product(shape(sorted_array))*kind(sorted_array),'sorted_array','sortRealArray')
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

      i_all=-product(shape(sorted_array))*kind(sorted_array)
      deallocate(sorted_array,stat=i_stat)
      call memocc(i_stat,i_all,'sorted_array','sortRealArray')
   end subroutine sortRealArray

   !-----------------------------------------------------------------------------
   ! SUBROUTINE: printEigenvectors
   !> @brief Print eigenvectors
   !> @note NOT FINISHED! for printing the eigenvectors
   !-----------------------------------------------------------------------------
   subroutine printEigenvectors(nq,NA)
      ! NOT FINISHED! for printing the eigenvectors
      implicit none

      integer, intent(in) :: nq,NA
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

   !-----------------------------------------------------------------------------
   ! SUBROUTINE: printA
   !> @brief Print dynamical matrix
   !> @note obsolete, used for troubleshooting
   !-----------------------------------------------------------------------------
   subroutine printA(A,nq,NA)
      !
      implicit none
      integer :: i,j
      character(LEN=11) :: nom = 'Amatrix.out'
      complex(dblprec), allocatable, dimension(:,:,:) :: A
      integer, intent(in) :: nq,NA
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

   !-----------------------------------------------------------------------------
   ! SUBROUTINE: magdos_calc
   !> @brief Calculate AMS density of states
   !-----------------------------------------------------------------------------
   !subroutine magdos_calc(filename,wres,magdos,iat,nq)
   subroutine magdos_calc(filename,wres,iat,nq)
      implicit none

      integer, intent(in) :: iat, nq
      character(LEN=*),intent(in) :: filename
      real(dblprec), intent(in), dimension(iat,nq) :: wres
      !real(dblprec), allocatable, dimension(:,:) :: magdos
      !
      integer :: i, k, ia, flines,i_stat,i_all
      real(dblprec) :: emin, emax, deltae, fac, dummy
      real(dblprec), allocatable,dimension(:,:) :: mtemp

      if(do_magdos=='A') then

         ! Open input file
         open(ifileno, file=magdosfile)
         flines=-1
         ! Pre-read file to get number of lines (frequencies)
         do
            read(ifileno,*,end=200)  dummy,dummy,dummy
            flines=flines+1
         end do
         200 continue
         rewind(ifileno)
         magdos_freq=flines
         if(.not.allocated(magdos)) then 
            allocate(magdos(magdos_freq,2),stat=i_stat)
            call memocc(i_stat,product(shape(magdos))*kind(magdos),'magdos','magdos_calc')
         end if
         read(ifileno,*)
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
            emax=maxval(wres)
         else
            emax=magdos_hfreq
         endif
         deltae=(emax-emin)/(magdos_freq-1)
         if(.not.allocated(magdos)) then 
            allocate(magdos(magdos_freq,2),stat=i_stat)
            call memocc(i_stat,product(shape(magdos))*kind(magdos),'magdos','magdos_calc')
         end if
         magdos=0.0_dblprec

         do i=1,magdos_freq
            magdos(i,1)=emin+(i-1)*deltae
         enddo

         do k=1,magdos_freq
            do i =1,nq
               do ia=1,iat
                  fac=exp(-(wres(ia,i)-magdos(k,1))**2/magdos_sigma)*q_weight(i)
                  magdos(k,2)=magdos(k,2)+fac
               enddo
            enddo
         enddo

         !Normalize DOS

         magdos(:,2)=magdos(:,2)/(sum(magdos(:,2))*deltae)
        write(ofileno,'(a)') &
            "#          E(meV)         D(S(q,E))         Int(D(S))  "
         do k=1,magdos_freq
            write(ofileno,1001) magdos(k,1),magdos(k,2),sum(magdos(1:k,2))*deltae
         enddo
      elseif(do_magdos=='N') then
      else    !S/B
         ! Open input file
         open(ifileno, file=magdosfile)
         flines=-1
         ! Pre-read file to get number of lines (frequencies)
         do
            read(ifileno,*,end=300)  
            flines=flines+1
         end do
         300 continue
         rewind(ifileno)
         magdos_freq=flines
         allocate(mtemp(magdos_freq,2),stat=i_stat)
         call memocc(i_stat,product(shape(mtemp))*kind(mtemp),'mtemp','magdos_calc')
         if(.not.allocated(magdos)) then 
            allocate(magdos(magdos_freq,2),stat=i_stat)
            call memocc(i_stat,product(shape(magdos))*kind(magdos),'magdos','magdos_calc')
         end if
         read(ifileno,*) 
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
         i_all=-product(shape(mtemp))*kind(mtemp)
         deallocate(mtemp,stat=i_stat)
         call memocc(i_stat,i_all,'mtemp','magdos_calc')

      endif

      1001 format (2x,f18.12,2f18.12)
      close(ofileno)
   end subroutine magdos_calc

   !-----------------------------------------------------------------------------
   ! FUNCTION: get_maxdist
   !> @brief q-point distances for nicer plots
   !-----------------------------------------------------------------------------
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
      dist=0.0_dblprec
      maxdist=0.0_dblprec
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

         case('magdos_rasamples')
            read(ifile,*,iostat=i_err) magdos_rasamples
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

!-----------------------------------------------------------------------------
! SUBROUTINE: init_ams
!> @brief Set default values for the AMS parameters
!-----------------------------------------------------------------------------
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
   magdos_rasamples=0
   !
end subroutine init_ams

!-----------------------------------------------------------------------------
! SUBROUTINE: calculate_random_ams
! DESCRIPTION
!> @brief
!> Main driver for adiabatic magnon spectra (AMS) in collinear magnets.
!> Based on linear spin wave theory as described in ....
!-----------------------------------------------------------------------------
subroutine calculate_random_ams()
   use iso_fortran_env
   !use Montecarlo, only : choose_random_atom_x
   use InputData,          only : N1,N2,N3,NA,NT,Natom,simid,hfield,gsconf_num,nchmax,ammom_inp,do_ralloy, ham_inp

   implicit none
   integer :: mu
   integer :: nu
   integer :: lambda
   integer :: i, ia, alfa, beta, delta

   integer :: countstart !< index for last atom "before" the cell where calculations are performed

   real(dblprec), dimension(3) :: zeroVector !(0,0,0) vector  to calculate J(0)
   real(dblprec), allocatable, dimension(:,:) :: wres,jqres,wres_diag,jqres_diag
   real(dblprec), allocatable, dimension(:) :: wres_tmp,jqres_tmp
   real(dblprec), allocatable, dimension(:) :: xconc  ! Actual concentration of impurities in the supercell

   real(dblprec) :: fc, fc2

   real(dblprec) :: mag_sign_mu, mag_sign_nu, mag_sign_lambda
   real(dblprec), dimension(3) :: order_vector
   complex(dblprec), allocatable, dimension(:,:,:) :: A !< matrix whose eigenvalues are the sought energies fro AMS
   complex(dblprec), allocatable, dimension(:,:,:) :: Ad !< matrix whose eigenvalues are the sought energies fro AMS
   complex(dblprec), allocatable, dimension(:,:,:) :: eigv !< matrix whose eigenvalues are the sought energies fro AMS
   complex(dblprec), allocatable, dimension(:,:,:) :: eigv_diag !< matrix whose eigenvalues are the sought energies fro AMS
   complex(dblprec), allocatable, dimension(:,:,:) :: B !< matrix whose eigenvalues are the sought energies for J(q)
   complex(dblprec), allocatable, dimension(:,:,:) :: Bd !< matrix whose eigenvalues are the sought energies for J(q)

   complex(dblprec), allocatable, dimension(:,:,:) :: Asc !< matrix whose eigenvalues are the sought energies fro AMS
   complex(dblprec), allocatable, dimension(:,:,:) :: Bsc !< matrix whose eigenvalues are the sought energies fro J(q)
   complex(dblprec), allocatable, dimension(:,:,:) :: Adsc !< matrix whose eigenvalues are the sought energies fro AMS
   complex(dblprec), allocatable, dimension(:,:,:) :: Bdsc !< matrix whose eigenvalues are the sought energies fro J(q)

   ! for output printing
   character(LEN = 19) :: output_file4
   character(LEN = 19) :: output_file5
   character(LEN = 22) :: output_file6
   character(LEN = 22) :: output_file7
   character(LEN = 22) :: output_file8
   character(LEN = 22) :: output_file9
   character(LEN = 20) :: output_file15
   character(LEN = 23) :: output_file16
   integer :: iq, idx,i_stat,i_all, maxiter, iia
   integer, dimension(:), allocatable :: atomlist

   ! Factors for mRy energy conversion
   fc = mry/mub
   fc2 = 2*mry/mub

   ! if no q allocated then no AMS calculations are performed
   if (.not. allocated(q)) then
      write(*,*) 'No q-file allocated. No Adiabatic Magnon Spectrum calculation is performed'
      goto 10 ! jumps to end of subroutine
   end if

   if (.not. allocated(q_weight)) then
      allocate(q_weight(nq),stat=i_stat)
      call memocc(i_stat,product(shape(q_weight))*kind(q_weight),'q_weight','calculate_random_ams')
      q_weight=1.0_dblprec
   end if


   allocate(atomlist(natom),stat=i_stat)
   call memocc(i_stat,product(shape(atomlist))*kind(atomlist),'atomlist','calculate_random_ams')
   allocate(A(nchmax,nchmax,nq),stat=i_stat)
   call memocc(i_stat,product(shape(A))*kind(A),'A','calculate_random_ams')
   allocate(Ad(na,na,nq),stat=i_stat)
   call memocc(i_stat,product(shape(Ad))*kind(Ad),'Ad','calculate_random_ams')
   allocate(Asc(nchmax,nchmax,natom),stat=i_stat)
   call memocc(i_stat,product(shape(Asc))*kind(Asc),'Asc','calculate_random_ams')
   allocate(Adsc(na,na,natom),stat=i_stat)
   call memocc(i_stat,product(shape(Adsc))*kind(Adsc),'Adsc','calculate_random_ams')
   allocate(eigv(nchmax,nchmax,nq),stat=i_stat)
   call memocc(i_stat,product(shape(eigv))*kind(eigv),'eigv','calculate_random_ams')
   allocate(eigv_diag(na,na,nq),stat=i_stat)
   call memocc(i_stat,product(shape(eigv_diag))*kind(eigv_diag),'eigv_diag','calculate_random_ams')
   allocate(B(nchmax,nchmax,nq),stat=i_stat)
   call memocc(i_stat,product(shape(B))*kind(B),'B','calculate_random_ams')
   allocate(Bd(na,na,nq),stat=i_stat)
   call memocc(i_stat,product(shape(Bd))*kind(Bd),'Bd','calculate_random_ams')
   allocate(Bsc(nchmax,nchmax,natom),stat=i_stat)
   call memocc(i_stat,product(shape(Bsc))*kind(Bsc),'Bsc','calculate_random_ams')
   allocate(Bdsc(na,na,natom),stat=i_stat)
   call memocc(i_stat,product(shape(Bdsc))*kind(Bdsc),'Bdsc','calculate_random_ams')
   allocate(wres(nchmax,nq),stat=i_stat)
   call memocc(i_stat,product(shape(wres))*kind(wres),'wres','calculate_random_ams')
   allocate(wres_tmp(na),stat=i_stat)
   call memocc(i_stat,product(shape(wres_tmp))*kind(wres_tmp),'wres_tmp','calculate_random_ams')
   allocate(wres_diag(na,nq),stat=i_stat)
   call memocc(i_stat,product(shape(wres_diag))*kind(wres_diag),'wres_diag','calculate_random_ams')
   allocate(jqres(nchmax,nq),stat=i_stat)
   call memocc(i_stat,product(shape(jqres))*kind(jqres),'jqres','calculate_random_ams')
   allocate(jqres_tmp(na),stat=i_stat)
   call memocc(i_stat,product(shape(jqres_tmp))*kind(jqres_tmp),'jqres_tmp','calculate_random_ams')
   allocate(jqres_diag(na,nq),stat=i_stat)
   call memocc(i_stat,product(shape(jqres_diag))*kind(jqres_diag),'jqres_diag','calculate_random_ams')
   allocate(xconc(nchmax),stat=i_stat)
   call memocc(i_stat,product(shape(xconc))*kind(xconc),'xconc','calculate_random_ams')
   A = 0.0_dblprec
   Ad = 0.0_dblprec
   B = 0.0_dblprec
   Bd = 0.0_dblprec
   Asc = 0.0_dblprec
   Bsc = 0.0_dblprec
   Adsc = 0.0_dblprec
   Bdsc = 0.0_dblprec
   zeroVector = 0.0_dblprec
   xconc=0.0_dblprec

   output_file4 = 'ams.'//simid//'.out'
   output_file5 = 'ra_ams.'//simid//'.out'
   output_file6 = 'ra_magdos.'//simid//'.out'
   output_file7 = 'ra_jqams.'//simid//'.out'
   output_file8 = 'magdos.'//simid//'.out'
   output_file9 = 'jqams.'//simid//'.out'
   output_file15 = 'ara_ams.'//simid//'.out'
   output_file16 = 'ara_magdos.'//simid//'.out'

   !write (*,'(1x,a)',advance="no") "Mounting A-matrices for random-AMS calculation"
   write (*,'(1x,a)',advance="yes") "Mounting A-matrices for random-AMS calculation"
   idx=0
   maxiter=(N1)*(N2)*(N3)*NQ
   !write(output_unit,'(2x,1i3,1a1,2x,a1,40a1)', advance='no') 100*idx/maxiter,'%','[', ('#', k =1,20*idx/maxiter)
   !close(output_unit);open(output_unit)
   !
   order_vector(1)=1.0_dblprec/sqrt(3.0_dblprec)
   order_vector(2)=1.0_dblprec/sqrt(3.0_dblprec)
   order_vector(3)=1.0_dblprec/sqrt(3.0_dblprec)

   ! Calculate actual concentration of the impurities in the supercell

   if (magdos_rasamples==0) then
      magdos_rasamples=natom
   else
      magdos_rasamples = min(magdos_rasamples,natom)
   end if

   do ia=1,magdos_rasamples
      xconc(achem_ch(ia))=xconc(achem_ch(ia))+1.0_dblprec/(magdos_rasamples)
   enddo

   do ia=1,natom
      atomlist(ia)=ia
   enddo

   ! todo: change to summation over relevant atoms
   msat=sum(mmom(:,1),1)/(natom)
   
   !call choose_random_atom_x(Natom,atomlist)

   do i=1,nq
   !$omp parallel do default(shared),private(ia,iia,countstart,mu,nu,lambda,mag_sign_mu,mag_sign_nu,mag_sign_lambda,alfa,beta,delta)
      do ia=1,magdos_rasamples
         iia=atomlist(ia) !+magdos_rasamples
         !countstart=ia !(ia-1)
         idx=idx+1
         !if(mod(idx-1,maxiter/20)==0) then
         !  write(output_unit,'(47a1)', advance='no') (char(8), k =1,(20*idx/maxiter)+8)
         !  write(output_unit,'(2x,1i3,1a1,2x,a1,40a1)', advance='no') 100*idx/maxiter,'%','[', ('#', k =1,20*idx/maxiter)
         !  close(output_unit);open(output_unit)
         !end if

         mu=achem_ch(iia)
         alfa=asite_ch(iia)
!         mag_sign_mu=sign(1.0_dblprec,sum(emom(1:3,ia,1)*order_vector))   
         mag_sign_mu=1.0_dblprec
         mag_sign_nu=1.0_dblprec
         mag_sign_lambda=1.0_dblprec

         do nu=1,nchmax
!            mag_sign_nu=sign(1.0_dblprec,sum(emom(1:3,ia,1)*order_vector)) ! Bogus at the moment, only FM
            Asc(mu,nu,ia) = -calc_jRA(mu,nu,q(:,i),iia)*ammom_inp(alfa,nu,gsconf_num)*mag_sign_nu
            Bsc(mu,nu,ia) = -Asc(mu,nu,ia)*ammom_inp(alfa,mu,gsconf_num)
            if (mu==nu) then
               do lambda = 1,nchmax
!                  mag_sign_lambda=sign(1.0_dblprec,sum(emom(1:3,ia,1)*order_vector)) ! Bogus
                  Asc(mu,mu,ia)=Asc(mu,mu,ia) + calc_jRA(mu,lambda,zeroVector,iia)* &
                      ammom_inp(alfa,lambda,gsconf_num)*mag_sign_lambda
               end do
               Asc(mu,mu,ia)=Asc(mu,mu,ia) + calc_ani(mu,mu,iia,ham_inp%do_anisotropy) ! Probably wrong
               Asc(mu,mu,ia)=Asc(mu,mu,ia) + sum(hfield*emomm(1:3,iia,1))*mag_sign_nu*0.5_dblprec ! No
            end if
            Asc(mu,nu,ia)=(Asc(mu,nu,ia)*Landeg(alfa)*Landeg(alfa))/xconc(mu)
            Bsc(mu,nu,ia)=(Bsc(mu,nu,ia)*Landeg(alfa)*Landeg(alfa))/xconc(mu)
         end do

         do beta=1,na
!            mag_sign_nu=sign(1.0_dblprec,sum(emom(1:3,ia,1)*order_vector)) ! Bogus at the moment, only FM
            Adsc(alfa,beta,ia) = -calc_jDRA(alfa,beta,q(:,i),iia)*ammom_inp(alfa,mu,gsconf_num)*mag_sign_nu
            Bdsc(alfa,beta,ia) = -Adsc(alfa,beta,ia)
            if (alfa==beta) then
               do delta = 1,na
!                  mag_sign_lambda=sign(1.0_dblprec,sum(emom(1:3,ia,1)*order_vector)) ! Bogus
                  Adsc(alfa,alfa,ia)=Adsc(alfa,alfa,ia) + calc_jDRA(alfa,delta,zeroVector,iia)*mag_sign_lambda*&
                          ammom_inp(alfa,mu,gsconf_num)
               end do
               Adsc(alfa,alfa,ia)=Adsc(alfa,alfa,ia) + calc_ani(alfa,alfa,iia,ham_inp%do_anisotropy) ! Probably wrong
               Adsc(alfa,alfa,ia)=Adsc(alfa,alfa,ia) + sum(hfield*emomm(1:3,iia,1))*mag_sign_nu*0.5_dblprec ! No
            end if
            Adsc(alfa,beta,ia)=Adsc(alfa,beta,ia)*Landeg(alfa)*Landeg(beta)
            Bdsc(alfa,beta,ia)=Bdsc(alfa,beta,ia)*Landeg(alfa)*Landeg(beta)
         end do
      end do
      !$omp end parallel do

      ! Testing to solve eigen value problem for each site and average
      ! eigenvalues instead of Hamiltonian
      wres_diag(:,i) = 0.0_dblprec
      jqres_diag(:,i) = 0.0_dblprec

      do ia=1,magdos_rasamples
         Ad(:,:,i) = 4.0_dblprec*Adsc(:,:,ia)/fc2*ry_ev/msat
         Bd(:,:,i) = Bdsc(:,:,ia)/fc2*ry_ev
         print *,ia,atomlist(ia),achem_ch(atomlist(ia))
         print '(2f12.6)',real(Ad(:,:,i))
         print *,'------------------------------'
         call eigenvalue_calculation_lapack(Ad(:,:,i),Bd(:,:,i),wres_tmp,jqres_tmp,eigv(:,:,i),na,1)
         wres_diag(:,i) = wres_diag(:,i) + wres_tmp
         jqres_diag(:,i) = jqres_diag(:,i) + jqres_tmp
      end do

!!!      wres(:,i) = 0.0_dblprec
!!!      jqres(:,i) = 0.0_dblprec
!!!
!!!      do ia=1,magdos_rasamples
!!!         A(:,:,i) = 4.0_dblprec*Asc(:,:,ia)/fc2*ry_ev !*sign(msat,1._dblprec)
!!!         B(:,:,i) = Bsc(:,:,ia)/fc2*ry_ev
!!!         print *,ia,atomlist(ia),achem_ch(atomlist(ia))
!!!         print '(2f12.6)',real(A(:,:,i))
!!!         print *,'------------------------------'
!!!         call eigenvalue_calculation_lapack(A(:,:,i),B(:,:,i),wres_tmp,jqres_tmp,eigv(:,:,i),nchmax,1)
!!!         wres(:,i) = wres(:,i) + wres_tmp
!!!         jqres(:,i) = jqres(:,i) + jqres_tmp
!!!      end do

      A(:,:,i)=sum(Asc,3)/(magdos_rasamples)
      B(:,:,i)=sum(Bsc,3)/(magdos_rasamples)
      Ad(:,:,i)=sum(Adsc,3)/(magdos_rasamples)
      Bd(:,:,i)=sum(Bdsc,3)/(magdos_rasamples)
   end do
       
   write(*,'(a)') " done."
   wres_diag=wres_diag/magdos_rasamples
   jqres_diag=jqres_diag/magdos_rasamples
   call printEnergies(output_file15,wres_diag,msat,tcmfa,tcrpa,na,1)
   call magdos_calc(output_file16,wres_diag,na,nq)

   !!!wres=wres/magdos_rasamples
   !!!jqres=jqres/magdos_rasamples
   !!!call printEnergies(output_file15,wres,msat,tcmfa,tcrpa,nchmax,1)
   !!!call magdos_calc(output_file16,wres,nchmax,nq)
   ! unit conversion from mRy to mEv, by division by Rydberg's constant in (eV) (and an factor of 4 for some reason..)
   ! L.B. The factor is fine, a factor of 2 from the def of dynamical matrix and fc2 includes another factor of 2
   ! L.B Temporary fix before replacing 4 with Landeg_glob^2=4
   A = 4.0_dblprec*A/fc2*ry_ev !*sign(msat,1._dblprec)
   B = B/fc2*ry_ev
   Ad = (4.0_dblprec*Ad/fc2*ry_ev)/msat !*sign(msat,1._dblprec)
   Bd = Bd/fc2*ry_ev

   ! eigenvalues of A (one set per NA) are written to wres  - AMS
   write (*,'(1x,a)',advance='yes') "Diagonalizing A-matrix for AMS calculation"
   ! The realness of A has been questioned.
   !According to Essenberger et. al PRB 84, 174425 (2011) Eqn. 22 it is correctly so.
   !do iq=1,nq
   !   A(:,:,iq)=real(A(:,:,iq))
   !   Ad(:,:,iq)=real(Ad(:,:,iq))
   !end do
   !print *,A
   call eigenvalue_calculation_lapack(A,B,wres,jqres,eigv,nchmax,nq)
   call eigenvalue_calculation_lapack(Ad,Bd,wres_diag,jqres_diag,eigv_diag,na,nq)

   !Take absolute values, in case of AF frequencies
   !Assumption: negative frequencies ok, still gives real energies (not complex frequencies/energies)
   !This is now taken care of in eigenvalue_calculation_lapack() where abs() is applied
   ! Sort eigenvalues and eigenvectors according to magnitude of eigenvalues
 !  do i=1,nq
 !     call sortEigenVVs(wres(1:na,i),eigv(1:na,1:na,i),na)
 !  end do

   ! print energies to output file
   call printEnergies(output_file5,wres,msat,tcmfa,tcrpa,nchmax,1)
   call printEnergies(output_file7,jqres,msat,tcmfa,tcrpa,nchmax,2)
   call magdos_calc(output_file6,wres,nchmax,nq)
   !call magdos_calc(output_file6,wres,magdos,nchmax,nq)

   write (*,'(1x,a)',advance='yes') "Diagonal disorder approximation"
   call printEnergies(output_file4,wres_diag,msat,tcmfa,tcrpa,na,1)
   call printEnergies(output_file9,jqres_diag,msat,tcmfa,tcrpa,na,2)
   call magdos_calc(output_file8,wres_diag,na,nq)
   !call magdos_calc(output_file8,wres_diag,magdos,na,nq)
   !!
   i_all=-product(shape(wres))*kind(wres)
   deallocate(wres,stat=i_stat)
   call memocc(i_stat,i_all,'wres','calculate_random_ams')
   i_all=-product(shape(wres_tmp))*kind(wres_tmp)
   deallocate(wres_tmp,stat=i_stat)
   call memocc(i_stat,i_all,'wres_tmp','calculate_random_ams')
   i_all=-product(shape(wres_diag))*kind(wres_diag)
   deallocate(wres_diag,stat=i_stat)
   call memocc(i_stat,i_all,'wres_diag','calculate_random_ams')
   i_all=-product(shape(jqres))*kind(jqres)
   deallocate(jqres,stat=i_stat)
   call memocc(i_stat,i_all,'jqres','calculate_random_ams')
   i_all=-product(shape(jqres_tmp))*kind(jqres_tmp)
   deallocate(jqres_tmp,stat=i_stat)
   call memocc(i_stat,i_all,'jqres_tmp','calculate_random_ams')
   i_all=-product(shape(jqres_diag))*kind(jqres_diag)
   deallocate(jqres_diag,stat=i_stat)
   call memocc(i_stat,i_all,'jqres_diag','calculate_random_ams')

   i_all=-product(shape(atomlist))*kind(atomlist)
   deallocate(atomlist,stat=i_stat)
   call memocc(i_stat,i_all,'atomlist','calculate_random_ams')

   i_all=-product(shape(A))*kind(A)
   deallocate(A,stat=i_stat)
   call memocc(i_stat,i_all,'A','calculate_random_ams')
   i_all=-product(shape(Ad))*kind(Ad)
   deallocate(Ad,stat=i_stat)
   call memocc(i_stat,i_all,'Ad','calculate_random_ams')
   i_all=-product(shape(Asc))*kind(Asc)
   deallocate(Asc,stat=i_stat)
   call memocc(i_stat,i_all,'Asc','calculate_random_ams')
   i_all=-product(shape(Adsc))*kind(Adsc)
   deallocate(Adsc,stat=i_stat)
   call memocc(i_stat,i_all,'Adsc','calculate_random_ams')
   i_all=-product(shape(eigv))*kind(eigv)
   deallocate(eigv,stat=i_stat)
   call memocc(i_stat,i_all,'eigv','calculate_random_ams')
   i_all=-product(shape(eigv_diag))*kind(eigv_diag)
   deallocate(eigv_diag,stat=i_stat)
   call memocc(i_stat,i_all,'eigv_diag','calculate_random_ams')
   i_all=-product(shape(B))*kind(B)
   deallocate(B,stat=i_stat)
   call memocc(i_stat,i_all,'B','calculate_random_ams')
   i_all=-product(shape(Bd))*kind(Bd)
   deallocate(Bd,stat=i_stat)
   call memocc(i_stat,i_all,'Bd','calculate_random_ams')
   i_all=-product(shape(Bsc))*kind(Bsc)
   deallocate(Bsc,stat=i_stat)
   call memocc(i_stat,i_all,'Bsc','calculate_random_ams')
   i_all=-product(shape(Bdsc))*kind(Bdsc)
   deallocate(Bdsc,stat=i_stat)
   call memocc(i_stat,i_all,'Bdsc','calculate_random_ams')
   i_all=-product(shape(xconc))*kind(xconc)
   deallocate(xconc,stat=i_stat)
   call memocc(i_stat,i_all,'xconc','calculate_random_ams')
   write(*,'(1x,a)') 'Random-Adiabatic Magnon Spectra Calculations done.'
   10 continue
end subroutine calculate_random_ams

end module AMS
