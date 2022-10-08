!------------------------------------------------------------------------------
!  MODULE: KMCData
!
!  DESCRIPTION:
!> This module contains the information about the necessary variables to perform
!> a KMC calculation. As well as a set of helper routines, such as array allocation
!> printing and reading of data.
!
!> @author
!> Jonathan Chico
!> @copyright
!> GNU Public License.
!------------------------------------------------------------------------------
module KMCData

   use Parameters
   use Profiling
   !
   implicit none
   !
   !From setup
   integer :: NA_KMC     !< Number of KMC particles found in the system
   integer :: kmc_method !< Which is the type of KMC algorithm tight now only kmc_method=1 (usual KMC) is used
   integer :: time_efield !< Time dependent electric field
   integer :: max_no_equiv_barriers  !< Calculated maximum of neighbours for the KMC barriers
   integer :: max_no_neigh_barriers  !< Calculated maximum of neighbours for KMC barriers
   integer :: max_no_shells_barriers !< Calculated maximum of shells for KMC barriers
   integer, dimension(:), allocatable :: ivisit_kmc !< Visit pattern for the neighbours
   integer, dimension(:), allocatable :: nn_barriers !< Number of neighbour shells for KMC barriers
   integer, dimension(:), allocatable :: achem_ch_KMC !< Chemical species for the KMC particles
   integer, dimension(:), allocatable :: atype_inp_KMC !< Atomic types for the KMC particles
   integer, dimension(:), allocatable :: kmc_index_prv !< Indices indicating positions of the KMC paritcles before and update
   integer, dimension(:), allocatable :: kmc_index_aft !< Indices indicating positions of the KMC paritcles after an update
   integer, dimension(:), allocatable :: kmc_time_steps !< Number of steps until next jump for all the KMC particles
   integer, dimension(:), allocatable :: nlistsize_barriers !< Size of neighbour list for KMC barriers
   integer, dimension(:,:), allocatable :: nlist_barriers   !< Neighbour list for KMC barriers

   real(dblprec) :: rate0 !< Attempt rate frequency, this is a mechanism dependent variable

   real(dblprec), dimension(3) :: efield !< Magnitude of the applied electric field (mostly for polarons)
   real(dblprec), dimension(:), allocatable :: rate, rate_full, rate_free
   real(dblprec), dimension(:,:), allocatable :: Bas_KMC !< Coordinates for basis atoms of the KMC particles
   real(dblprec), dimension(:,:), allocatable :: ncoup_barriers !< KMC barriers
   real(dblprec), dimension(:,:,:), allocatable ::redcoord_barriers
   real(dblprec), dimension(:,:,:,:,:), allocatable :: kmc_barriers !< KMC barriers input

   character(len=1) :: do_kmc !< Perform KMC calculation (Y/N)
   character(len=1) :: do_efield !< Include external electric field (Y/N)
   character(len=35) :: barrfile !< File name fo the KMC energy barriers
   character(len=35) :: kmc_posfile !< Position files for the KMC particles

   ! Counters for the measurements
   integer :: kmc_buff !< Size of the buffer for the KMC particles
   integer :: kmc_step !< Time interval for the printing of the KMC particles
   integer :: bcount_kmc !< Counter of buffer for KMC info
   integer, dimension(:,:), allocatable :: kmc_posb      !< Buffer for the polaron position
   integer, dimension(:,:), allocatable :: time_stepsb   !< Buffer for the time steps between jumps
   real(dblprec), dimension(:), allocatable :: indxb_kmc !< Step counter for the KMC info
   real(dblprec), dimension(:,:,:), allocatable :: KMC_ave ! < Buffer for the KMC prticles average magnetization
   real(dblprec), dimension(:,:,:), allocatable :: kmc_coordb !< Coordinates for the KMC particles
   character(len=1) :: do_prn_kmc !< Flag for printing the KMC info
   logical :: KMC_changed !< Flag to determine if there is going to be a printing

   integer, dimension(:), allocatable :: kmc_part_shuffle !< Array for shuffeling the update sequence to avoid biased movement

   public

contains

   !------------------------------------------------------------------------------
   ! SUBROUTINE: allocate_kmc
   !
   ! DESCRIPTION:
   !> In this routine the lists for the barriers are allocated
   !
   !> @author
   !> Jonathan Chico
   !------------------------------------------------------------------------------
   subroutine allocate_kmc(NA_KMC,Natom,max_no_neigh_barriers,flag)
      use MonteCarlo, only : choose_random_atom_x
      implicit none

      integer, intent(in) :: NA_KMC !< Number of KMC particles found in the system
      integer, intent(in), optional :: Natom                 !< Number of atoms in system
      integer, intent(in), optional :: max_no_neigh_barriers !< Calculated maximum of neighbours for exchange
      integer, intent(in) :: flag  !< Allocate or deallocate (1/-1)

      integer :: i_all, i_stat, i
      ! Energy_barriers
      if(flag>0) then
         allocate(nlistsize_barriers(Natom),stat=i_stat)
         call memocc(i_stat,product(shape(nlistsize_barriers))*kind(nlistsize_barriers),'nlistsize_barriers','allocate_kmc')
         nlistsize_barriers=0
         allocate(nlist_barriers(max_no_neigh_barriers,Natom),stat=i_stat)
         call memocc(i_stat,product(shape(nlist_barriers))*kind(nlist_barriers),'nlist_barriers','allocate_kmc')
         nlist_barriers=0
         allocate(ncoup_barriers(max_no_neigh_barriers,Natom),stat=i_stat)
         call memocc(i_stat,product(shape(ncoup_barriers))*kind(ncoup_barriers),'ncoup_barriers','allocate_kmc')
         ncoup_barriers=0.0_dblprec
         allocate(Rate(max_no_neigh_barriers),stat=i_stat)
         call memocc(i_stat,product(shape(Rate))*kind(Rate),'Rate','allocate_kmc')
         Rate=0.0_dblprec
         allocate(Rate_full(max_no_neigh_barriers),stat=i_stat)
         call memocc(i_stat,product(shape(Rate_full))*kind(Rate_full),'Rate_full','allocate_kmc')
         Rate_full=0.0_dblprec
         allocate(Rate_free(0:max_no_neigh_barriers),stat=i_stat)
         call memocc(i_stat,product(shape(Rate_free))*kind(Rate_free),'Rate_free','allocate_kmc')
         Rate_free=0.0_dblprec
         allocate(kmc_time_steps(NA_KMC),stat=i_stat)
         call memocc(i_stat,product(shape(kmc_time_steps))*kind(kmc_time_steps),'kmc_time_steps','allocate_kmc')
         kmc_time_steps=0
         allocate(kmc_index_prv(NA_KMC),stat=i_stat)
         call memocc(i_stat,product(shape(kmc_index_prv))*kind(kmc_index_prv),'kmc_index_prv','allocate_kmc')
         kmc_index_prv=0
         allocate(kmc_index_aft(NA_KMC),stat=i_stat)
         call memocc(i_stat,product(shape(kmc_index_aft))*kind(kmc_index_aft),'kmc_index_aft','allocate_kmc')
         kmc_index_aft=0
         allocate(kmc_part_shuffle(NA_KMC),stat=i_stat)
         call memocc(i_stat,product(shape(kmc_part_shuffle))*kind(kmc_part_shuffle),'kmc_part_shuffle','allocate_kmc')
         do i=1,NA_KMC
            kmc_part_shuffle(i)=i
         end do
         !     call choose_random_atom_x(NA_KMC,kmc_part_shuffle)

         allocate(ivisit_kmc(max_no_neigh_barriers),stat=i_stat)
         call memocc(i_stat,product(shape(ivisit_kmc))*kind(ivisit_kmc),'ivisit_kmc','allocate_kmc')
         ivisit_kmc=0

      else
         i_all=-product(shape(nlistsize_barriers))*kind(nlistsize_barriers)
         deallocate(nlistsize_barriers,stat=i_stat)
         call memocc(i_stat,i_all,'nlistsize_barriers','allocate_kmc')
         i_all=-product(shape(nlist_barriers))*kind(nlist_barriers)
         deallocate(nlist_barriers,stat=i_stat)
         call memocc(i_stat,i_all,'nlist_barriers','allocate_kmc')
         i_all=-product(shape(ncoup_barriers))*kind(ncoup_barriers)
         deallocate(ncoup_barriers,stat=i_stat)
         call memocc(i_stat,i_all,'ncoup_barriers','allocate_kmc')
         i_all=-product(shape(Rate))*kind(Rate)
         deallocate(Rate,stat=i_stat)
         call memocc(i_stat,i_all,'Rate','allocate_kmc')
         i_all=-product(shape(kmc_time_steps))*kind(kmc_time_steps)
         deallocate(kmc_time_steps,stat=i_stat)
         call memocc(i_stat,i_all,'kmc_time_steps','allocate_kmc')
         i_all=-product(shape(Bas_KMC))*kind(Bas_KMC)
         deallocate(Bas_KMC,stat=i_stat)
         call memocc(i_stat,i_all,'Bas_KMC','kmc_part_init')
         i_all=-product(shape(atype_inp_KMC))*kind(atype_inp_KMC)
         deallocate(atype_inp_KMC,stat=i_stat)
         call memocc(i_stat,i_all,'atype_inp_KMC','kmc_part_init')
         i_all=-product(shape(achem_ch_KMC))*kind(achem_ch_KMC)
         deallocate(achem_ch_KMC,stat=i_stat)
         call memocc(i_stat,i_all,'achem_ch_KMC','kmc_part_init')
         i_all=-product(shape(kmc_index_prv))*kind(kmc_index_prv)
         deallocate(kmc_index_prv,stat=i_stat)
         call memocc(i_stat,i_all,'kmc_index_prv','kmc_part_init')
         i_all=-product(shape(kmc_index_aft))*kind(kmc_index_aft)
         deallocate(kmc_index_aft,stat=i_stat)
         call memocc(i_stat,i_all,'kmc_index_aft','kmc_part_init')
         i_all=-product(shape(kmc_part_shuffle))*kind(kmc_part_shuffle)
         deallocate(kmc_part_shuffle,stat=i_stat)
         call memocc(i_stat,i_all,'kmc_part_shuffle','kmc_part_init')
         i_all=-product(shape(ivisit_kmc))*kind(ivisit_kmc)
         deallocate(ivisit_kmc,stat=i_stat)
         call memocc(i_stat,i_all,'ivisit_kmc','kmc_part_init')
      end if

   end subroutine allocate_kmc

   !-----------------------------------------------------------------------------
   ! SUBROUTINE: init_kmc
   !
   ! DESCRIPTION:
   !> Initialization of the kmc variables for algorithm and printing
   !
   !> @author
   !> Jonathan Chico
   !-----------------------------------------------------------------------------
   subroutine init_kmc()

      implicit none

      do_kmc           = 'N'
      do_efield        = 'N'
      do_prn_kmc       = 'N'
      kmc_posfile      = 'kmc_posfile'
      barrfile         = 'barrfile'
      kmc_step         = 100
      kmc_buff         = 10
      kmc_method       = 1
      rate0            = 0.0_dblprec
      efield           = 0.0_dblprec
      KMC_changed      =.true.
      time_efield      = 0

   end subroutine init_kmc

   !-----------------------------------------------------------------------------
   !  SUBROUTINE: read_positions_KMC
   !> @brief
   !> Read Positions for the KMC particles
   !
   !> @author
   !> Jonathan Chico, based on the routine by Anders Bergman
   !-----------------------------------------------------------------------------
   subroutine read_positions_KMC(C1,C2,C3,posfiletype)
      !
      implicit none

      real(dblprec), dimension(3), intent(in) :: C1 !< First lattice vector
      real(dblprec), dimension(3), intent(in) :: C2 !< Second lattice vector
      real(dblprec), dimension(3), intent(in) :: C3 !< Third lattice vector
      character(len=1), intent(in) :: posfiletype
      !
      integer :: flines,itype,mtype,iat,isite,msite,i_stat
      real(dblprec),dimension(3) :: tmp

      ! Open input file
      open(ifileno, file=kmc_posfile)
      ! Check if input file is for random alloy
      flines=0
      mtype=0
      msite=0
      ! Pre-read file to get max no. sites, types and chemical species
      do
         read(ifileno,*,end=200)  isite,itype
         flines=flines+1
         msite=max(msite,isite)
         mtype=max(mtype,itype)
      end do
      200 continue
      rewind(ifileno)

      ! Set no. sites, types and chemical species
      if(msite/=flines) write(*,*) "WARNING: Check input file ", kmc_posfile, &
         " for inconsistent information."
      NA_KMC=msite
      ! Really? I would expect nchmax=1 if no random alloy and if so it will be set later.
      ! Allocate input arrays
      allocate(Bas_KMC(3,NA_KMC),stat=i_stat)
      call memocc(i_stat,product(shape(Bas_KMC))*kind(Bas_KMC),'Bas_KMC','read_positions_KMC')
      Bas_KMC=0.0_dblprec
      allocate(atype_inp_KMC(NA_KMC),stat=i_stat)
      call memocc(i_stat,product(shape(atype_inp_KMC))*kind(atype_inp_KMC),'atype_inp_KMC','read_positions_KMC')
      atype_inp_KMC=0

      ! Read basis atoms and setup type array
      ! Site, Type, Rx, Ry, Rz
      if (posfiletype=='C') then
         do iat=1, NA_KMC
            read (ifileno,*) isite, itype,  Bas_KMC(1,isite), Bas_KMC(2,isite),&
               Bas_KMC(3,isite)
            atype_inp_KMC(isite)=itype
         enddo
      elseif (posfiletype=='D') then
         do iat=1, NA_KMC
            read (ifileno,*) isite, itype,  tmp(1),tmp(2),tmp(3)
            atype_inp_KMC(isite)=itype
            Bas_KMC(1,isite)=tmp(1)*C1(1)+tmp(2)*C2(1)+tmp(3)*C3(1)
            Bas_KMC(2,isite)=tmp(1)*C1(2)+tmp(2)*C2(2)+tmp(3)*C3(2)
            Bas_KMC(3,isite)=tmp(1)*C1(3)+tmp(2)*C2(3)+tmp(3)*C3(3)

         enddo
      endif
      close (ifileno)
   end subroutine read_positions_KMC

   !-----------------------------------------------------------------------------
   !  SUBROUTINE: read_positions_alloy_KMC
   !> @brief
   !> Read Positions for the KMC particles in the case of random alloy
   !
   !> @author
   !> Jonathan Chico, based on the routine by Anders Bergman
   !-----------------------------------------------------------------------------
   subroutine read_positions_alloy_KMC(C1,C2,C3,posfiletype)
      !
      implicit none

      real(dblprec), dimension(3), intent(in) :: C1 !< First lattice vector
      real(dblprec), dimension(3), intent(in) :: C2 !< Second lattice vector
      real(dblprec), dimension(3), intent(in) :: C3 !< Third lattice vector
      character(len=1), intent(in) :: posfiletype
      !
      integer :: flines,isite,itype,ichem
      integer :: msite,mtype,mchem,iat, i_stat
      real(dblprec) :: rconc
      real(dblprec),dimension(3) :: tmp

      ! Open input file
      open(ifileno, file=kmc_posfile)

      ! Check if input file is for random alloy
      flines=0
      msite=0
      mtype=0
      mchem=0

      ! Pre-read file to get max no. sites, types and chemical species
      do
         read(ifileno,*,end=200)  isite,itype,ichem
         flines=flines+1
         msite=max(msite,isite)
         mtype=max(mtype,itype)
         mchem=max(mchem,ichem)
      end do
      200 continue
      rewind(ifileno)

      ! Set no. sites, types and chemical species
      NA_KMC=msite

      ! Allocate input arrays
      allocate(Bas_KMC(3,msite),stat=i_stat)
      call memocc(i_stat,product(shape(Bas_KMC))*kind(Bas_KMC),'Bas_KMC','read_positions_alloy_KMC')
      Bas_KMC=0.0_dblprec
      allocate(atype_inp_KMC(NA_KMC),stat=i_stat)
      call memocc(i_stat,product(shape(atype_inp_KMC))*kind(atype_inp_KMC),'atype_inp_KMC','read_positions_alloy_KMC')
      atype_inp_KMC=0
      ! Chemical data
      allocate(achem_ch_KMC(NA_KMC),stat=i_stat)
      call memocc(i_stat,product(shape(achem_ch_kmc))*kind(achem_ch_kmc),'achem_ch_kmc','read_positions_alloy_KMC')
      achem_ch_KMC=0

      ! Read data
      ! Site  Type   Chem_type Conc, Rx, Ry, Rz
      if (posfiletype=='C') then
         do iat=1, flines
            read (ifileno,*) isite, itype, ichem, rconc, Bas_KMC(1,isite), Bas_KMC(2,isite), Bas_KMC(3,isite)
            atype_inp_KMC(isite)=itype
            achem_ch_KMC(isite)=ichem
         enddo
      elseif(posfiletype=='D') then
         do iat=1, flines
            read (ifileno,*) isite, itype,ichem,rconc,  tmp(1),tmp(2),tmp(3)
            Bas_KMC(1,isite)=tmp(1)*C1(1)+tmp(2)*C2(1)+tmp(3)*C3(1)
            Bas_KMC(2,isite)=tmp(1)*C1(2)+tmp(2)*C2(2)+tmp(3)*C3(2)
            Bas_KMC(3,isite)=tmp(1)*C1(3)+tmp(2)*C2(3)+tmp(3)*C3(3)
            atype_inp_KMC(isite)=itype
            achem_ch_KMC(isite)=ichem
         enddo
      endif
      close (ifileno)
   end subroutine read_positions_alloy_KMC

   !-----------------------------------------------------------------------------
   !  SUBROUTINE: print_kmc
   !> @brief
   !> Wrapper for printing the information from kinetic Monte Carlo
   !
   !> @author
   !> Jonathan Chico
   !-----------------------------------------------------------------------------
   subroutine print_kmc(sstep,mstep,Natom,Mensemble,do_ralloy,Natom_full,delta_t,&
         max_no_neigh,nlistsize,nlist,real_time_measure,emomM,simid,atype,&
         coord,achem_ch)

      implicit none

      integer, intent(in) :: mstep !< Current simulation step
      integer, intent(in) :: sstep !< Current simulation step in logarithmic scale
      integer, intent(in) :: Natom !< Number of atoms in system
      !  integer, intent(in) :: NA_KMC
      integer, intent(in) :: Mensemble !< Number of ensembles
      integer, intent(in) :: do_ralloy !< Random alloy simulation (0/1)
      integer, intent(in) :: Natom_full  !< Number of atoms for full system (=Natom if not dilute)
      integer, intent(in) :: max_no_neigh !< Calculated maximum of neighbours for exchange
      integer, dimension(Natom), intent(in) :: nlistsize !< Size of neighbour list for Heisenberg exchange couplings
      integer, dimension(Natom), intent(in) :: atype !< Type of atoms
      !  integer, dimension(NA_KMC), intent(in) :: kmc_index_prv
      integer, dimension(Natom_full), intent(in) :: achem_ch  !< Chemical type of atoms (reduced list)
      integer, dimension(max_no_neigh,Natom), intent(in) :: nlist !< Neighbour list for Heisenberg exchange couplings
      real(dblprec), intent(in) :: delta_t !< Time step
      real(dblprec), dimension(3,Natom), intent(in) :: coord !< Coordinates of atoms
      real(dblprec), dimension(3,Natom, Mensemble), intent(in) :: emomM  !< Current magnetic moment vector
      character(len=8), intent(in) :: simid !< Name of simulation
      character(len=1), intent(in) :: real_time_measure !< Measurements displayed in real time

      ! .. In/out variables
      !  logical, intent(inout) :: KMC_changed !< Flag to determine if there is going to be a printing

      ! Adaptive KMC printing
      if (do_prn_kmc=='D') then
         if (KMC_changed) then
            ! Wite the KMC information to the buffer
            call buffer_kmc(mstep-1,Natom,NA_KMC,do_ralloy,Mensemble,Natom_full,&
               bcount_kmc,max_no_neigh,atype,nlistsize,kmc_index_prv,achem_ch,nlist,&
               delta_t,coord,emomM,real_time_measure)
            KMC_changed=.false.
            if (bcount_kmc==kmc_buff) then
               ! Print the buffer the KMC buffer to file
               call prn_kmc(simid,real_time_measure)
               bcount_kmc=1
            else
               bcount_kmc=bcount_kmc+1
            endif
         endif

         ! Static KMC printing
      else if(do_prn_kmc=='Y') then
         if (mod(sstep-1,kmc_step)==0) then
            ! Wite the KMC information to the buffer
            call buffer_kmc(mstep-1,Natom,NA_KMC,do_ralloy,Mensemble,Natom_full,&
               bcount_kmc,max_no_neigh,atype,nlistsize,kmc_index_prv,achem_ch,nlist,&
               delta_t,coord,emomM,real_time_measure)
            KMC_changed=.false.
            if (bcount_kmc==kmc_buff) then
               ! Print the buffer the KMC buffer to file
               call prn_kmc(simid,real_time_measure)
               bcount_kmc=1
            else
               bcount_kmc=bcount_kmc+1
            endif
         endif
      end if

   end subroutine print_kmc

   !-----------------------------------------------------------------------------
   !  SUBROUTINE: allocate_prn_kmc
   !> @brief
   !> Allocate data for the printing of the KMC info
   !
   !> @author
   !> Anders Bergman
   !> Jonathan Chico --> Separated into its own individual module
   !-----------------------------------------------------------------------------
   subroutine allocate_prn_kmc(flag)

      implicit none

      integer, intent(in) :: flag       !< Allocate or deallocate (1/-1)

      !.. Local variables
      integer :: i_stat, i_all

      if (flag>0) then

         bcount_kmc=1
         if (do_prn_kmc.eq.'Y'.or.do_prn_kmc.eq.'D') then
            allocate(kmc_posb(NA_KMC,kmc_buff),stat=i_stat)
            call memocc(i_stat,product(shape(kmc_posb))*kind(kmc_posb),'kmc_posb','allocate_prn_mkc')
            kmc_posb=0
            allocate(kmc_coordb(3,NA_KMC,kmc_buff),stat=i_stat)
            call memocc(i_stat,product(shape(kmc_coordb))*kind(kmc_coordb),'kmc_coordb','allocate_prn_mkc')
            kmc_coordb=0.0_dblprec
            allocate(indxb_kmc(kmc_buff),stat=i_stat)
            call memocc(i_stat,product(shape(indxb_kmc))*kind(indxb_kmc),'indxb_kmc','allocate_prn_mkc')
            indxb_kmc=0
            allocate(KMC_ave(3,NA_KMC,kmc_buff),stat=i_stat)
            call memocc(i_stat,product(shape(KMC_ave))*kind(KMC_ave),'KMC_ave','allocate_prn_kmc')
            KMC_ave=0.0_dblprec
            allocate(time_stepsb(NA_KMC,kmc_buff),stat=i_stat)
            call memocc(i_stat,product(shape(time_stepsb))*kind(time_stepsb),'time_stepsb','allocate_prn_kmc')
            time_stepsb=0
         endif

      else

         if (do_prn_kmc.eq.'Y'.or.do_prn_kmc.eq.'D') then
            i_all=-product(shape(kmc_posb))*kind(kmc_posb)
            deallocate(kmc_posb,stat=i_stat)
            call memocc(i_stat,i_all,'kmc_posb','allocate_prn_kmc')
            i_all=-product(shape(indxb_kmc))*kind(indxb_kmc)
            deallocate(indxb_kmc,stat=i_stat)
            call memocc(i_stat,i_all,'indxb_kmc','allocate_prn_kmc')
            i_all=-product(shape(KMC_ave))*kind(KMC_ave)
            deallocate(KMC_ave,stat=i_stat)
            call memocc(i_stat,i_all,'KMC_ave','allocate_prn_kmc')
            i_all=-product(shape(time_stepsb))*kind(time_stepsb)
            deallocate(time_stepsb,stat=i_stat)
            call memocc(i_stat,i_all,'time_stepsb','allocate_prn_kmc')
            i_all=-product(shape(kmc_coordb))*kind(kmc_coordb)
            deallocate(kmc_coordb,stat=i_stat)
            call memocc(i_stat,i_all,'kmc_coordb','allocate_prn_kmc')
         endif

      end if

   end subroutine allocate_prn_kmc

   !-----------------------------------------------------------------------------
   !  SUBROUTINE: flush_prn_kmc
   !> @brief
   !> Subroutine to flush the KMC information for printing
   !
   !> @author
   !> Jonathan Chico
   !-----------------------------------------------------------------------------
   subroutine flush_prn_kmc(Natom,simid,real_time_measure)

      implicit none

      integer, intent(in) :: Natom !< Number of atoms in the system
      character(len=8), intent(in) :: simid              !< Name of simulation
      character(len=1), intent(in) :: real_time_measure  !< Measurements displayed in real time

      ! Print the KMC information
      if (do_prn_kmc.eq.'Y'.or.do_prn_kmc.eq.'D') then
         bcount_kmc=bcount_kmc-1
         ! Write the total skyrmion number buffer to file
         call prn_kmc(simid,real_time_measure)

      endif

   end subroutine flush_prn_kmc

   !-----------------------------------------------------------------------------
   !  SUBROUTINE: buffer_kmc
   !> @brief
   !> Subroutine to buff the KMC information
   !
   !> @author
   !> Jonathan Chico
   !-----------------------------------------------------------------------------
   subroutine buffer_kmc(mstep,Natom,NA_KMC,do_ralloy,Mensemble,Natom_full,&
         bcount_kmc,max_no_neigh,atype,nlistsize,kmc_index_prv,achem_ch,nlist,&
         delta_t,coord,emomM,real_time_measure)

      implicit none

      integer, intent(in) :: mstep  !< Current simulation step
      integer, intent(in) :: Natom  !< Number of atoms in system
      integer, intent(in) :: NA_KMC !< Number of KMC particles found in the system
      integer, intent(in) :: do_ralloy !< Random alloy simulation (0/1)
      integer, intent(in) :: Mensemble   !< Number of ensembles
      integer, intent(in) :: Natom_full  !< Number of atoms for full system (=Natom if not dilute)
      integer, intent(in) :: bcount_kmc  !< Counter of buffer for the KMC information
      integer, intent(in) :: max_no_neigh !< Calculated maximum of neighbours for exchange
      integer, dimension(Natom), intent(in) :: atype !< Type of atoms
      integer, dimension(Natom), intent(in) :: nlistsize !< Size of neighbour list for Heisenberg exchange couplings
      integer, dimension(NA_KMC), intent(in) :: kmc_index_prv !< Indices indicating positions of the KMC paritcles after an update
      integer, dimension(Natom_full), intent(in) :: achem_ch  !< Chemical type of atoms (reduced list)
      integer, dimension(max_no_neigh,Natom), intent(in) :: nlist !< Neighbour list for Heisenberg exchange couplings
      real(dblprec), intent(in) :: delta_t !< Current time step
      real(dblprec), dimension(3,Natom), intent(in) :: coord !< Coordinates of atoms
      real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: emomM  !< Current moment vector
      character(len=1), intent(in) :: real_time_measure !< Measurements displayed in real time

      ! .. Local variables
      integer :: ii,KMC_part_size,kmc_part,kmc_part_pos

      if (do_ralloy.eq.0) then
         ! Loop over the KMC particles
         do kmc_part=1, NA_KMC
            kmc_part_pos=kmc_index_prv(kmc_part)
            KMC_part_size=0
            ! Average for the current KMC particle position
            do ii=1, nlistsize(kmc_part_pos)
               ! Make sure that only atoms that belong to the actual polaron cloud are counted
               if (atype(kmc_part_pos).ne.atype(nlist(ii,kmc_part_pos))) write(*,*) ii, kmc_part_pos
               if (atype(kmc_part_pos).eq.atype(nlist(ii,kmc_part_pos))) then
                  KMC_ave(1,kmc_part,bcount_kmc)=KMC_ave(1,kmc_part,bcount_kmc)+emomM(1,nlist(ii,kmc_part_pos),1)
                  KMC_ave(2,kmc_part,bcount_kmc)=KMC_ave(2,kmc_part,bcount_kmc)+emomM(2,nlist(ii,kmc_part_pos),1)
                  KMC_ave(3,kmc_part,bcount_kmc)=KMC_ave(3,kmc_part,bcount_kmc)+emomM(3,nlist(ii,kmc_part_pos),1)
                  KMC_part_size=KMC_part_size+1
               endif
            enddo
            KMC_ave(1:3,kmc_part,bcount_kmc)=KMC_ave(1:3,kmc_part,bcount_kmc)/KMC_part_size
            kmc_posb(kmc_part,bcount_kmc)=kmc_part_pos
            kmc_coordb(1:3,kmc_part,bcount_kmc)=coord(1:3,kmc_part_pos)
            time_stepsb(kmc_part,bcount_kmc)=kmc_time_steps(kmc_part)
         enddo
      else
         ! Loop over the KMC particles
         do kmc_part=1, NA_KMC
            kmc_part_pos=kmc_index_prv(kmc_part)
            KMC_part_size=0
            ! Average for the current KMC particle position
            do ii=1, nlistsize(kmc_part_pos)
               ! Make sure that only atoms that belong to the actual polaron cloud are counted
               KMC_ave(1,kmc_part,bcount_kmc)=KMC_ave(1,kmc_part,bcount_kmc)+emomM(1,nlist(ii,kmc_part_pos),1)
               KMC_ave(2,kmc_part,bcount_kmc)=KMC_ave(2,kmc_part,bcount_kmc)+emomM(2,nlist(ii,kmc_part_pos),1)
               KMC_ave(3,kmc_part,bcount_kmc)=KMC_ave(3,kmc_part,bcount_kmc)+emomM(3,nlist(ii,kmc_part_pos),1)
               KMC_part_size=KMC_part_size+1
            enddo
            KMC_ave(1:3,kmc_part,bcount_kmc)=KMC_ave(1:3,kmc_part,bcount_kmc)/KMC_part_size
            kmc_posb(kmc_part,bcount_kmc)=kmc_part_pos
            kmc_coordb(1:3,kmc_part,bcount_kmc)=coord(1:3,kmc_part_pos)
            time_stepsb(kmc_part,bcount_kmc)=kmc_time_steps(kmc_part)
         enddo
      endif

      if (real_time_measure=='Y') then
         indxb_kmc(bcount_kmc)=nint(real(mstep,dblprec)*delta_t)
      else
         indxb_kmc(bcount_kmc)=mstep
      endif

   end subroutine buffer_kmc

   !-----------------------------------------------------------------------------
   !  SUBROUTINE: prn_kmc
   !> @brief
   !> Printing the type dependent projected skyrmion number
   !
   !> @author
   !> Anders Bergman
   !> Jonathan Chico --> Separated into its own individual module
   !-----------------------------------------------------------------------------
   subroutine prn_kmc(simid,real_time_measure)
      !
      !.. Implicit declarations
      implicit none

      character(len=8), intent(in) :: simid !< Simulation name
      character(len=1), intent(in) :: real_time_measure !< Measurements displayed in real time

      ! Local variables
      integer :: k,ii
      character(len=30) :: filn

      ! write info about KMC
      write(filn,'(''kmc_info.'',a,''.out'')') trim(simid)
      open(ofileno,file=filn, position="append")

      if (real_time_measure=='Y') then
         do k=1, bcount_kmc
            do ii=1,NA_KMC
               write (ofileno,2401) ii,indxb_kmc(k), kmc_posb(ii,k),kmc_coordb(1,ii,k),&
                  kmc_coordb(2,ii,k),kmc_coordb(3,ii,k),time_stepsb(ii,k),&
                  KMC_ave(1,ii,k),KMC_ave(2,ii,k),KMC_ave(3,ii,k)
            enddo
         enddo
      else
         do k=1, bcount_kmc
            do ii=1,NA_KMC
               write (ofileno,2400) ii,int(indxb_kmc(k)), kmc_posb(ii,k),kmc_coordb(1,ii,k),&
                  kmc_coordb(2,ii,k),kmc_coordb(3,ii,k),time_stepsb(ii,k),&
                  KMC_ave(1,ii,k),KMC_ave(2,ii,k),KMC_ave(3,ii,k)
            enddo
         enddo
      endif

      close(ofileno)

      return

      write(*,*) "Error writing the KMC info file"

      2400 format(i8,2x,i8,2x,i6,2x,f10.4,f10.4,f10.4,2x,i8,2x,f10.4,f10.4,f10.4)
      2401 format(es16.4,2x,i8,2x,i6,2x,f10.4,f10.4,f10.4,2x,i8,2x,f10.4,f10.4,f10.4)

   end subroutine prn_kmc

end module KMCData
