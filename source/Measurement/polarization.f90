!> Data and routines for measuring polarization vector
!> @author
!> Anders Bergman, Johan Hellsvik, Jonathan Chico
!> @copyright
!! Copyright (C) 2008-2018 UppASD group
!! This file is distributed under the terms of the
!! GNU General Public License. 
!! See http://www.gnu.org/copyleft/gpl.txt
module Polarization
   use Parameters
   use Profiling
   !
   implicit none
   !
   ! Variables from input
   integer :: pol_buff        !< Buffer size for the polarization
   integer :: pol_step        !< Interval for sampling the polarization
   character(len=1) :: do_pol !< Do polarization
   character(len=1) :: do_loc_pol !< Measure polarization vector locally
   character(len=1) :: do_chiral !< Measure polarization vector locally
   real(dblprec), dimension(:,:,:), allocatable :: eij !< Normalized neighbour distances
   real(dblprec), dimension(:,:,:), allocatable :: local_pol !< Data for local polarization vector
   real(dblprec), dimension(:,:,:), allocatable :: s_cross_s !< Data for local polarization vector
   !
   integer :: max_pol_nn                               !< Maximum number of polarization neighbours
   integer, dimension(:), allocatable :: pollistsize     !< Size of neighbour list for polarization calculations

   ! Variables needed for printing
   integer :: bcount_pol     !< Counter of buffer for for polarization
   real(dblprec), dimension(:), allocatable :: indxb_pol     !< Step counter for polarization vector
   real(dblprec), dimension(:,:,:), allocatable :: spol_buff  !< Buffer data for polarization vector

   private :: indxb_pol,spol_buff,bcount_pol

contains



   !!!!---------------init_polarization----------------!!!!!
   !> Allocate counters for sampling
   !> Allocate arrays for sampling data
   subroutine init_polarization(Natom, Mensemble, max_no_neigh, &
         nlist, coord, flag)

      implicit none

      integer, intent(in) :: flag          !< Allocate or deallocate (1/-1)
      integer, intent(in) :: Natom         !< Number of atoms in system
      integer, intent(in) :: Mensemble     !< Number of ensembles

      integer, intent(in) :: max_no_neigh  !< Calculated maximum of neighbours for exchange
      integer, dimension(max_no_neigh,Natom), intent(in) :: nlist !< Neighbour list for Heisenberg exchange couplings
      real(dblprec), dimension(3,Natom) :: coord !< Coordinates for atoms
      !real(dblprec), dimension(nt,max_no_shells,3) :: redcoord !< Coordinates for Heisenberg exchange couplings


      integer :: i,j,n, i_stat, i_all
      real(dblprec) :: dx, dy, dz, dnorm

      if(flag>0) then
         allocate(eij(3,max_no_neigh,Natom),stat=i_stat)
         call memocc(i_stat,product(shape(eij))*kind(eij),'eij','init_polarization')
         eij=0.0d0

         if(do_loc_pol=='Y') then
            allocate(local_pol(3,Natom,Mensemble),stat=i_stat)
            call memocc(i_stat,product(shape(local_pol))*kind(local_pol),'local_pol','init_polarization')
         end if
         if(do_chiral=='Y') then
            allocate(s_cross_s(3,Natom,Mensemble),stat=i_stat)
            call memocc(i_stat,product(shape(s_cross_s))*kind(s_cross_s),'s_cross_s','init_polarization')
         end if
         allocate(pollistsize(Natom),stat=i_stat)
         call memocc(i_stat,product(shape(pollistsize))*kind(pollistsize),'pollistsize','allocate_hamiltoniandata')
         pollistsize=max_pol_nn
         ! Set up normalized neighbour vectors
         do i=1,Natom
            do n=1,pollistsize(i)
               j=nlist(n,i)
               dx=coord(1,j)-coord(1,i)
               dy=coord(2,j)-coord(2,i)
               dz=coord(3,j)-coord(3,i)
               dnorm=sqrt(dx*dx+dy*dy+dz*dz)
               eij(1,n,i)=dx/dnorm
               eij(2,n,i)=dy/dnorm
               eij(3,n,i)=dz/dnorm
            end do
         end do
      else
         i_all=-product(shape(eij))*kind(eij)
         deallocate(eij,stat=i_stat)
         call memocc(i_stat,i_all,'eij','init_polarization')

         if(do_loc_pol=='Y') then
            i_all=-product(shape(local_pol))*kind(local_pol)
            deallocate(local_pol,stat=i_stat)
            call memocc(i_stat,i_all,'local_pol','init_polarization')
         end if
         if(do_chiral=='Y') then
            i_all=-product(shape(s_cross_s))*kind(s_cross_s)
            deallocate(s_cross_s,stat=i_stat)
            call memocc(i_stat,i_all,'s_cross_s','init_polarization')
         end if
         deallocate(pollistsize,stat=i_stat)
         call memocc(i_stat,i_all,'pollistsize','allocate_hamiltoniandata')
      end if

   end subroutine init_polarization



   !!!!---------------allocate_prn_pol----------------!!!!!
   !> Allocation and initialization of the necesary arrays for printing
   subroutine allocate_prn_pol(Mensemble,flag)

      implicit none

      integer, intent(in) :: flag      !< Allocate or deallocate (1/-1)
      integer, intent(in) :: Mensemble !< Number of ensembles

      integer :: i_stat, i_all

      if (flag>0) then

         bcount_pol=1

         if (do_pol=='Y') then
            allocate(spol_buff(3,pol_buff,Mensemble),stat=i_stat)
            call memocc(i_stat,product(shape(spol_buff))*kind(spol_buff),'spol_buff','allocate_prn_pol')
            allocate(indxb_pol(pol_buff),stat=i_stat)
            call memocc(i_stat,product(shape(indxb_pol))*kind(indxb_pol),'indxb_pol','allocate_prn_pol')
         endif

      else

         if (do_pol=='Y') then
            i_all=-product(shape(spol_buff))*kind(spol_buff)
            deallocate(spol_buff,stat=i_stat)
            call memocc(i_stat,i_all,'spol_buff','allocate_prn_pol')
            i_all=-product(shape(indxb_pol))*kind(indxb_pol)
            deallocate(indxb_pol,stat=i_stat)
            call memocc(i_stat,i_all,'indxb_pol','allocate_prn_pol')
         endif

      endif

   end subroutine allocate_prn_pol



   !!!!---------------prin_pol_init----------------!!!!!
   !> Initialization of the default variables for the polarization printing
   subroutine prn_pol_init()

      implicit none

      do_pol            = 'N'
      pol_step          = 100
      pol_buff          = 10
      do_loc_pol        = 'N'
      do_chiral         = 'N'
      max_pol_nn        =  6

   end subroutine prn_pol_init



   !!!!---------------print_pol----------------!!!!!
   !> Wrapper subroutine for the polarization printing
   subroutine print_pol(sstep,mstep,Natom,Mensemble,max_no_neigh,nlist,nlistsize,&
         emomM,delta_t,simid,real_time_measure)

      implicit none

      integer, intent(in) :: sstep                                !< Current simulation step in logarithmic form
      integer, intent(in) :: mstep                                !< Current simulation step
      integer, intent(in) :: Natom                                !< Number of atoms in the system
      integer, intent(in) :: Mensemble                            !< Number of ensembles
      integer, intent(in) :: max_no_neigh                         !< Maximum number of neighbours for the neighbour lists
      integer, dimension(max_no_neigh,Natom), intent(in) :: nlist !< Neighbour list for Heisenberg exchange couplings
      integer, dimension(Natom),intent(in) :: nlistsize           !< Size of neighbour list for Heisenberg exchange couplings
      real(dblprec), intent(in) :: delta_t                        !< Time step for real time measurement
      real(dblprec), dimension(3,Natom, Mensemble), intent(in) :: emomM   !< Current magnetic moment vector
      character(len=8), intent(in) :: simid                       !< Name of simulation
      character(len=1), intent(in) :: real_time_measure           !< Measurements displayed in real time

      ! Polarization
      if (do_pol=='Y') then
         if ( mod(sstep-1,pol_step)==0) then

            ! Write step to buffer
            call buffer_pol(Natom, Mensemble,mstep-1, emomM, max_no_neigh, nlist, delta_t,real_time_measure)

            if (do_loc_pol=='Y') call measure_local_pol(Natom, Mensemble, emomM, max_no_neigh, nlist, nlistsize)
            if (do_chiral=='Y')  call measure_chirality(Natom, Mensemble, emomM, max_no_neigh, nlist, nlistsize)

            if (bcount_pol==pol_buff) then
               ! Write buffer to file
               call prn_pol(Natom, Mensemble, simid,real_time_measure)
               if (do_loc_pol=='Y') call print_local_polarity_chirality(Natom, Mensemble, mstep, simid,'P')
               if (do_chiral=='Y')  call print_local_polarity_chirality(Natom, Mensemble, mstep, simid,'C')

               ! Reset statistics buffer
               bcount_pol=1
            else
               bcount_pol=bcount_pol+1
            endif

         endif
      end if

   end subroutine print_pol



   !!!!---------------flush_polarization----------------!!!!!
   !> Flush the polarization measurements, i.e. print to file the polarization
   subroutine flush_polarization(mstep,Natom,Mensemble,simid,real_time_measure)

      implicit none

      integer, intent(in) :: mstep                      !< Current simulation step
      integer, intent(in) :: Natom                      !< Number of atoms in system
      integer, intent(in) :: Mensemble                  !< Number of ensembles
      character(len=8), intent(in) :: simid             !< Name of simulation
      character(len=1), intent(in) :: real_time_measure !< Measurements displayed in real time

      ! Polarization
      if (do_pol=='Y') then
         ! Write buffer to file
         bcount_pol=bcount_pol-1
         call prn_pol(Natom, Mensemble, simid,real_time_measure)

         if (do_loc_pol=='Y') call print_local_polarity_chirality(Natom, Mensemble, mstep, simid,'P')
         if (do_chiral=='Y')  call print_local_polarity_chirality(Natom, Mensemble, mstep, simid,'C')

      endif

   end subroutine flush_polarization




   !!!!---------------buffer_pol----------------!!!!!
   !> Buffer average magnetization
   subroutine buffer_pol(Natom, Mensemble, mstep, emomM, max_no_neigh, nlist, delta_t,real_time_measure)
      !
      !.. Implicit declarations
      implicit none

      integer, intent(in) :: mstep !< Current simulation step
      integer, intent(in) :: Natom !< Number of atoms in system
      integer, intent(in) :: Mensemble !< Number of ensembles
      integer, intent(in) :: max_no_neigh !< Calculated maximum of neighbours for exchange
      integer, dimension(max_no_neigh,Natom), intent(in) :: nlist !< Neighbour list for Heisenberg exchange couplings
      real(dblprec), intent(in) :: delta_t !< Current time step
      real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: emomM  !< Current magnetic moment vector
      character(len=1), intent(in) :: real_time_measure !< Measurements displayed in real time

      !.. Scalar variables
      integer :: i,j,k,n

      !.. Local arrays
      real(dblprec), dimension(Mensemble) ::  px, py, pz
      real(dblprec) ::  tx, ty, tz

      !.. Executable statements
      px(:)=0.0d0
      py(:)=0.0d0
      pz(:)=0.0d0

      !.. Sum over moments
      do k=1,Mensemble
         do i=1, Natom
            do n=1,pollistsize(i)
               j=nlist(n,i)
               tx = emomM(2,i,k)*emomM(3,j,k)-emomM(3,i,k)*emomM(2,j,k)
               ty = emomM(3,i,k)*emomM(1,j,k)-emomM(1,i,k)*emomM(3,j,k)
               tz = emomM(1,i,k)*emomM(2,j,k)-emomM(2,i,k)*emomM(1,j,k)
               px(k) = px(k) + eij(2,n,i)*tz-eij(3,n,i)*ty
               py(k) = py(k) + eij(3,n,i)*tx-eij(1,n,i)*tz
               pz(k) = pz(k) + eij(1,n,i)*ty-eij(2,n,i)*tx
            end do
         end do
      end do

      !.. Save in buffer
      do k=1,Mensemble
         spol_buff(1,bcount_pol,k) = px(k)
         spol_buff(2,bcount_pol,k) = py(k)
         spol_buff(3,bcount_pol,k) = pz(k)
      end do

      if (real_time_measure=='Y')  then
         indxb_pol(bcount_pol)=mstep*delta_t
      else
         indxb_pol(bcount_pol)=mstep
      endif

   end subroutine buffer_pol



   !!!!---------------measure_local_pol---------------!!!!!
   subroutine measure_local_pol(Natom, Mensemble, emomM, max_no_neigh, nlist, nlistsize)
      !
      !.. Implicit declarations
      implicit none

      integer, intent(in) :: Natom !< Number of atoms in system
      integer, intent(in) :: Mensemble !< Number of ensembles
      real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: emomM  !< Current magnetic moment vector
      integer, intent(in) :: max_no_neigh !< Calculated maximum of neighbours for exchange
      integer, dimension(max_no_neigh,Natom), intent(in) :: nlist !< Neighbour list for Heisenberg exchange couplings
      integer, dimension(Natom),intent(in) :: nlistsize !< Size of neighbour list for Heisenberg exchange couplings

      !.. Scalar variables
      integer :: i,j,k,n

      !.. Local arrays
      real(dblprec) ::  tx, ty, tz

      local_pol=0.0d0
      !.. Sum over moments
      do k=1,Mensemble
         do i=1, Natom
            do n=1,nlistsize(i)
               j=nlist(n,i)
               tx = emomM(2,i,k)*emomM(3,j,k)-emomM(3,i,k)*emomM(2,j,k)
               ty = emomM(3,i,k)*emomM(1,j,k)-emomM(1,i,k)*emomM(3,j,k)
               tz = emomM(1,i,k)*emomM(2,j,k)-emomM(2,i,k)*emomM(1,j,k)
               local_pol(1,i,k) = local_pol(1,i,k) + eij(2,n,i)*tz-eij(3,n,i)*ty
               local_pol(2,i,k) = local_pol(2,i,k) + eij(3,n,i)*tx-eij(1,n,i)*tz
               local_pol(3,i,k) = local_pol(3,i,k) + eij(1,n,i)*ty-eij(2,n,i)*tx
            end do
         end do
      end do

   end subroutine measure_local_pol



   !!!!---------------measure_chirality----------------!!!!!
   subroutine measure_chirality(Natom, Mensemble, emomM,  max_no_neigh, nlist, nlistsize)
      !
      !.. Implicit declarations
      implicit none

      integer, intent(in) :: Natom !< Number of atoms in system
      integer, intent(in) :: Mensemble !< Number of ensembles
      real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: emomM  !< Current magnetic moment vector
      integer, intent(in) :: max_no_neigh !< Calculated maximum of neighbours for exchange
      integer, dimension(max_no_neigh,Natom), intent(in) :: nlist !< Neighbour list for Heisenberg exchange couplings
      integer, dimension(Natom),intent(in) :: nlistsize !< Size of neighbour list for Heisenberg exchange couplings

      !.. Scalar variables
      integer :: i,j,k,n

      s_cross_s=0.0d0
      !.. Sum over moments
      do k=1,Mensemble
         do i=1, Natom
            do n=1,nlistsize(i)
               j=nlist(n,i)
               s_cross_s(1,i,k) = emomM(2,i,k)*emomM(3,j,k)-emomM(3,i,k)*emomM(2,j,k)
               s_cross_s(2,i,k) = emomM(3,i,k)*emomM(1,j,k)-emomM(1,i,k)*emomM(3,j,k)
               s_cross_s(3,i,k) = emomM(1,i,k)*emomM(2,j,k)-emomM(2,i,k)*emomM(1,j,k)
            end do
         end do
      end do

   end subroutine measure_chirality




   !!!!---------------prn_pol----------------!!!!!
   !> Print polarization vector
   subroutine prn_pol(Natom, Mensemble, simid,real_time_measure)
      !
      use Constants

      !.. Implicit declarations
      implicit none

      integer, intent(in) :: Natom     !< Number of atoms in system
      integer, intent(in) :: Mensemble !< Number of ensembles
      character(len=8), intent(in) :: simid !< Name of simulation
      character(len=1), intent(in) :: real_time_measure !< Measurements displayed in real time

      !.. Scalar variables
      integer :: i,k
      character(len=30) :: filn
      real(dblprec), dimension(Mensemble) :: avrgp
      real(dblprec), dimension(Mensemble) :: avrgpx, avrgpy, avrgpz


      real(dblprec) :: avrgpxm, avrgpym, avrgpzm, avrgpm
      real(dblprec) :: avrgps

      !.. Executable statements
      write (filn,'(''polarization.'',a8,''.out'')') simid
      open(ofileno, file=filn, position="append")

      ! Write header to output files for first iteration
      if(abs(indxb_pol(1))<dbl_tolerance) then
         write(ofileno,'(a)') "# Iter.     P_avg_x           P_avg_y         P_avg_z         P_avg  "
      end if

      do i=1, bcount_pol
         do k=1, Mensemble
            avrgpx(k) = spol_buff(1,i,k)/Natom
            avrgpy(k) = spol_buff(2,i,k)/Natom
            avrgpz(k) = spol_buff(3,i,k)/Natom
         end do

         do k=1,Mensemble
            avrgp(k) = avrgpx(k)**2+avrgpy(k)**2+avrgpz(k)**2
            avrgp(k) = sqrt(avrgp(k))
         end do

         ! mean polarizations
         avrgpxm=0d0
         avrgpym=0d0
         avrgpzm=0d0
         avrgpm=0d0

         !stdev polarization
         avrgps=0d0

         do k=1,Mensemble
            avrgpxm = avrgpxm + avrgpx(k)
            avrgpym = avrgpym + avrgpy(k)
            avrgpzm = avrgpzm + avrgpz(k)
            avrgpm = avrgpm + avrgp(k)
            avrgps = avrgps + avrgp(k)**2
         end do

         avrgpxm=avrgpxm/Mensemble
         avrgpym=avrgpym/Mensemble
         avrgpzm=avrgpzm/Mensemble
         avrgpm=avrgpm/Mensemble
         avrgps=avrgps/Mensemble - avrgpm**2

         if (real_time_measure=='Y')  then
            write (ofileno,10005) indxb_pol(i), avrgpxm, avrgpym, avrgpzm, avrgpm, avrgps
         else
            write (ofileno,10004) int(indxb_pol(i)), avrgpxm, avrgpym, avrgpzm, avrgpm, avrgps
         endif

      end do

      close(ofileno)

      return

      write (*,*) "Error writing the polarization file"

      10004 format (i8,6es16.8)
      10005 format (es16.4,6es16.8)

   end subroutine prn_pol



   !!!!---------------print_local_polarity_chirality----------------!!!!!
   !> Prints magnetic moment to restart file
   subroutine print_local_polarity_chirality(Natom, Mensemble, mstep, simid, printType)
      implicit none

      integer, intent(in) :: Natom      !< Number of atoms in system
      integer, intent(in) :: Mensemble  !< Number of ensembles
      integer, intent(in) :: mstep      !< Current simulation step
      character(len=8) :: simid         !< Name of simulation
      character(len=1) :: printType     !< Type of function to print

      !.. Scalar variables
      integer :: i,k
      character(len=30) :: filn

      real(dblprec) ::  tz

      !.. Executable statements

      if (printType=='C') then       ! Print chirality
         write (filn,'(''crestart.'',a8,''.out'')') simid
         open(ofileno, file=filn)
         do k=1,Mensemble
            do i=1, Natom
               tz=sum(s_cross_s(1:3,i,k)**2)**0.5d0+1.0d-12
               write (ofileno,10004) mstep, i, s_cross_s(1:3,i,k)/tz, sum(s_cross_s(1:3,i,k)**2)**0.5d0
            end do
         end do
      else if (printType=='P') then   ! Print local polarity
         write (filn,'(''prestart.'',a8,''.out'')') simid
         open(ofileno, file=filn)
         do k=1,Mensemble
            do i=1, Natom
               write (ofileno,10004) mstep, i, local_pol(1:3,i,k), sum(local_pol(1:3,i,k)**2)**0.5d0
            end do
         end do
      else
         write(*,*) '*** Unrecognized printType ****'
         stop
      end if

      close(ofileno)
      return

      10004 format (i8,i8,4es16.8)
   end subroutine print_local_polarity_chirality

end module polarization
