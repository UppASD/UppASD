!------------------------------------------------------------------------------------
! MODULE: prn_averages
!> @brief Module and data structures needed for printing the averages and cumulants
!> @todo Remove the ensemble energy arrays, as the energy module should take care of that
!> @copyright
!> GNU Public License
!------------------------------------------------------------------------------------
module prn_averages

   use AutoCorrelation
   use Energy, only : ene

   implicit none

   ! Printing definitions
   integer :: cumu_step !< Interval for sampling Binder cumulant
   integer :: cumu_buff !< Buffer size for Binder cumulant
   integer :: avrg_step !< Interval for sampling average magnetization
   integer :: avrg_buff !< Buffer size for average magnetization
   character(len=1) :: do_avrg                  !< Measure average magnetization (Y/N)
   character(len=1) :: do_cumu                  !< Measure Binder cumulant, susceptibility, and specific heat(Y/N)
   character(len=1) :: do_cuda_avrg             !< Measure average magnetization (Y/N) with CUDA
   character(len=1) :: do_cuda_cumu             !< Measure Binder cumulant, susceptibility, and specific heat(Y/N) with CUDA
   character(len=1) :: do_proj_avrg             !< Measure projected averages (Y/A/N)
   character(len=1) :: do_projch_avrg           !< Measure chemically projected averages (Y/N)
   character(len=1) :: do_cumu_proj             !< Measure Binder cumulant, susceptibility, and specific heat(Y/N)

   ! Local calculations for printing
   integer :: Navrgcum        !< Counter for number of cumulated averages
   real(dblprec) :: cumuw     !< Weight for current sample to cumulant
   real(dblprec) :: cumutotw  !< Sum of all cumulant weights 
   integer :: Navrgcum_proj   !< Counter for number of cumulated projected averages
   real(dblprec) :: mavg      !< Average magnetic moment
   real(dblprec) :: totene    !< Total energy
   real(dblprec) :: binderc   !< Binder cumulant
   real(dblprec) :: avrglcum  !< Cumulated average of l
   real(dblprec) :: avrgmcum  !< Cumulated average of m
   real(dblprec) :: avrgecum  !< Cumulated average of E
   real(dblprec) :: avrgetcum !< Cumulated average of E_xc
   real(dblprec) :: avrgelcum !< Cumulated average of E_LSF
   real(dblprec) :: avrgm4cum !< Cumulated average of m^4
   real(dblprec) :: avrgm2cum !< Cumulated average of m^2
   real(dblprec) :: avrgl4cum !< Cumulated average of l^4
   real(dblprec) :: avrgl2cum !< Cumulated average of l^2
   real(dblprec) :: avrge2cum !< Cumulated average of E^2

   real(dblprec), dimension(:), allocatable       :: indxb_avrg       !< Step counter for average magnetization
   real(dblprec), dimension(:,:,:), allocatable   :: mavg_buff        !< Buffer for average magnetizations
   real(dblprec), dimension(:,:,:), allocatable   :: mavg2_buff_proj  !< Buffer for squared projected averages
   real(dblprec), dimension(:,:,:,:), allocatable :: mavg_buff_proj   !< Buffer for projected averages
   real(dblprec), dimension(:,:,:,:), allocatable :: mavg_buff_projch !< Buffer for chemical projected averages

   real(dblprec), dimension(:), allocatable :: avrgm4cum_proj !< Cumulated average of projected m^4
   real(dblprec), dimension(:), allocatable :: avrgm2cum_proj !< Cumulated average of projected m^2
   real(dblprec), dimension(:), allocatable :: avrgmcum_proj  !< Cumulated average of projected m

   integer :: bcount_avrg    !< Counter of buffer for averages

   !public
   private

   ! Private variables
   public :: print_averages, flush_averages, averages_allocations, avrg_init, calc_and_print_cumulant
   public :: read_parameters_averages,zero_cumulant_counters
   public :: do_avrg, mavg, binderc,  avrg_step, do_cumu, cumu_step, cumu_buff, do_cuda_avrg, do_cuda_cumu
   public :: Navrgcum

contains

   !---------------------------------------------------------------------------------
   !> @brief
   !> Wrapper for printing the averages and cumulants
   !
   !> @author
   !> Anders Bergman
   !---------------------------------------------------------------------------------
   subroutine print_averages(NA,NT,N1,N2,N3,Natom,mstep,sstep,Nchmax,Mensemble,     &
      do_ralloy,Natom_full,plotenergy,atype,achem_ch,asite_ch,mmom,emom,emomM,Temp, &
      temprescale,temprescalegrad,delta_t,real_time_measure,simid,mode)

      use Restart, only : prn_mag_conf

      implicit none

      integer, intent(in) :: NA         !< Number of atoms in one cell
      integer, intent(in) :: NT         !< Number of types of atoms
      integer, intent(in) :: N1         !< Number of cell repetitions in x direction
      integer, intent(in) :: N2         !< Number of cell repetitions in y direction
      integer, intent(in) :: N3         !< Number of cell repetitions in z direction
      integer, intent(in) :: Natom      !< Number of atoms in the system
      integer, intent(in) :: mstep      !< Current simulation step
      integer, intent(in) :: sstep      !< Simulation step in logarithmic scale
      integer, intent(in) :: Nchmax     !< Max number of chemical components on each site in cell
      integer, intent(in) :: Mensemble  !< Number of ensembles
      integer, intent(in) :: do_ralloy  !< Random alloy simulation (0/1)
      integer, intent(in) :: Natom_full !< Number of atoms for full system (=Natom if not dilute)
      integer, intent(in) :: plotenergy !< Calculate and plot energy (0/1)
      real(dblprec), intent(in) :: Temp                      !< Temperature of the system (scalar)
      real(dblprec), intent(in) :: temprescale              !< Temperature rescaling if QHB
      real(dblprec), intent(in) :: temprescalegrad          !< Temperature rescaling gradient if QHB
      real(dblprec), intent(in) :: delta_t                   !< Current time step
      character(len=1), intent(in) :: mode   !< Simulation mode (S=SD, M=MC, H=MC Heat Bath, P=LD, C=SLD, G=GNEB)
      character(len=8), intent(in) :: simid             !< Simulation name
      character(len=1), intent(in) :: real_time_measure !< Measurement perfomred in real time
      integer, dimension(Natom), intent(in) :: atype         !< Type of atom
      integer, dimension(Natom), intent(in) :: achem_ch      !< Chemical type of atoms (reduced list) (achem_ch(i)=achtype(acellnumbrev(i)))
      integer, dimension(Natom_full), intent(in) :: asite_ch !< Actual site of atom for dilute system
      real(dblprec), dimension(Natom, Mensemble), intent(in) :: mmom     !< Magnitude of magnetic moments
      real(dblprec), dimension(3,Natom, Mensemble), intent(in) :: emom   !< Current unit moment vector
      real(dblprec), dimension(3,Natom, Mensemble), intent(in) :: emomM  !< Current magnetic moment vector
      ! .. Local variables
      integer :: k

      ! Calculate starting average magnetization
      if(mstep==0) then
         mavg=0.0_dblprec
         do k=1,Mensemble
            mavg = mavg + sqrt(sum(emomM(1,:,k))**2+sum(emomM(2,:,k))**2+sum(emomM(3,:,k))**2)/natom/Mensemble
         end do
      end if

      ! Averages
      if (do_avrg=='Y') then
         if ( mod(sstep-1,avrg_step)==0) then
            ! Write average step to buffer
            call buffer_avrg(Natom,Mensemble,mstep-1,emom,emomM,bcount_avrg,delta_t,real_time_measure)
            ! Write projected average step to buffer
            if (do_proj_avrg=='Y' .or. do_proj_avrg=='A') then
               call buffer_proj_avrg(Natom,Mensemble,NA,emomM,bcount_avrg,do_ralloy,&
                  Natom_full,asite_ch)
            endif
            ! Write chemically projected average step to buffer
            if (do_projch_avrg=='Y') then
               call buffer_projch_avrg(Natom,Mensemble,Nchmax,achem_ch,emomM,bcount_avrg)
            endif

            if (bcount_avrg==avrg_buff) then
               ! Write the total averages buffer to file
               call prn_avrg(Natom, Mensemble, simid,real_time_measure)
               if (do_proj_avrg=='Y' .or. do_proj_avrg=='A') then
                  call prn_proj_avrg(Natom,Mensemble,NA,N1,N2,N3,simid,NT,atype,    &
                     do_proj_avrg,real_time_measure)
               endif
               ! Write the chemically projected averages buffer to file
               if (do_projch_avrg=='Y') then
                  call prn_projch_avrg(Mensemble,N1,N2,N3,simid,Nchmax,real_time_measure)
               endif
               ! Create a restart file
               call prn_mag_conf(Natom,mstep,Mensemble,'R',simid,mmom,emom,'',mode)
               ! Reset statistics buffer
               bcount_avrg=1
            else
               bcount_avrg=bcount_avrg+1
            endif
         endif
      endif

      if (do_autocorr=='Y') then
         if ( mod(sstep-1,ac_step)==0.or.any(spinwaitt==sstep-1)) then
            call buffer_autocorr(Natom,Mensemble,do_autocorr,mstep,nspinwait,spinwait,emom,&
               emomM,bcount_ac,delta_t,real_time_measure)
            if (bcount_ac==ac_buff) then
               ! Write the autocorrelation buffer to file
               call prn_autocorr(Natom,simid,do_autocorr,nspinwait,real_time_measure)
               ! Reset statistics buffer
               bcount_ac=1
            else
               bcount_ac=bcount_ac+1
            endif
         endif
      end if
      ! Binder cumulant, susceptibility, and specific heat
      if(mod(sstep,cumu_step)==0) then
         if(do_cumu=='Y') then
            call calc_and_print_cumulant(Natom,Mensemble,emomM,simid,Temp,          &
               temprescale,temprescalegrad,plotenergy,cumu_buff,.true.)
         elseif(do_cumu=='A') then
            call calc_and_print_afm_cumulant(Natom,Mensemble,emomM,simid,Temp,      &
               temprescale,temprescalegrad,plotenergy,cumu_buff,NA,do_ralloy,Natom_full,asite_ch)
         end if

         if(do_cumu_proj=='Y') then
            call calc_and_print_cumulant_proj(Natom,Mensemble,emomM,simid,Temp,     &
               cumu_buff,NT,atype)
         end if
      endif

   end subroutine print_averages

   !---------------------------------------------------------------------------------
   ! SUBROUTINE: flush_averages
   !> Flush the averages, cumulants and autocorrelation measurements, i.e. print them to file in the last iteration
   !---------------------------------------------------------------------------------
   subroutine flush_averages(NT,NA,N1,N2,N3,mstep,Natom,Nchmax,Mensemble,atype,mmom,&
         emom,simid,real_time_measure,mode)

      use Restart, only : prn_mag_conf

      implicit none

      integer, intent(in) :: NT        !< Number of types of atoms
      integer, intent(in) :: NA        !< Number of atoms in one cell
      integer, intent(in) :: N1        !< Number of cell repetitions in x direction
      integer, intent(in) :: N2        !< Number of cell repetitions in y direction
      integer, intent(in) :: N3        !< Number of cell repetitions in z direction
      integer, intent(in) :: mstep     !< Current simulation step
      integer, intent(in) :: Natom     !< Number of atoms in the system
      integer, intent(in) :: Nchmax    !< Max number of chemical components on each site in cell
      integer, intent(in) :: Mensemble !< Number of ensembles
      character(len=1), intent(in) :: mode   !< Simulation mode (S=SD, M=MC, H=MC Heat Bath, P=LD, C=SLD, G=GNEB)
      character(len=8), intent(in) :: simid                !< Name of simulation
      character(len=1), intent(in) :: real_time_measure    !< Measurements displayed in real time
      integer, dimension(Natom), intent(in) :: atype       !< Type of atom
      real(dblprec), dimension(Natom, Mensemble), intent(in) :: mmom   !< Magnitude of magnetic moments
      real(dblprec), dimension(3,Natom, Mensemble), intent(in) :: emom !< Current unit moment vector

      ! Averages
      if (do_avrg=='Y') then
         ! Write buffer to file
         bcount_avrg=bcount_avrg-1
         call prn_avrg(Natom, Mensemble, simid,real_time_measure)


         if (do_proj_avrg=='Y' .or. do_proj_avrg=='A') then
            call prn_proj_avrg(Natom,Mensemble,NA,N1,N2,N3,simid,NT,atype,          &
               do_proj_avrg,real_time_measure)
         endif
         if (do_projch_avrg=='Y') then
            call prn_projch_avrg(Mensemble,N1,N2,N3,simid,Nchmax,real_time_measure)
         endif

         ! Create a restart file
         call prn_mag_conf(Natom,mstep,Mensemble,'R',simid,mmom,emom,'',mode)
         ! Reset statistics buffer
         bcount_avrg=1
      endif
      if (do_autocorr=='Y') then
         bcount_ac=bcount_ac+1
         call prn_autocorr(Natom, simid, do_autocorr, nspinwait, real_time_measure)
         bcount_ac=1
      end if

   end subroutine flush_averages

   !---------------------------------------------------------------------------------
   ! SUBROUTINE: avrg_init
   !> Initialization of the printing statements
   !---------------------------------------------------------------------------------
   subroutine avrg_init()

      implicit none

      !Averages variables
      do_avrg        = 'Y'
      do_proj_avrg   = 'N'
      do_projch_avrg = 'N'
      avrg_step      = 100
      avrg_buff      = 10
      !Cumulants variables
      do_cumu        = 'N'
      do_cumu_proj   = 'N'
      cumu_step      = 50
      cumu_buff      = 10
      !Autocorrelation variables
      ac_step        = 100
      ac_buff        = 10
      do_autocorr    = 'N'
      twfile         = 'twfile'

   end subroutine avrg_init

   !---------------------------------------------------------------------------------
   ! SUBROUTINE: zero_cumulant_counters
   !> Initialize cumulant counters
   !---------------------------------------------------------------------------------
   subroutine zero_cumulant_counters()
      !
      implicit none
      Navrgcum          = 0
      cumuw             = 0.0_dblprec
      cumutotw          = 0.0_dblprec
      Navrgcum_proj     = 0
      avrgmcum          = 0.0_dblprec
      avrgm2cum         = 0.0_dblprec
      avrgm4cum         = 0.0_dblprec
      avrglcum          = 0.0_dblprec
      avrgl2cum         = 0.0_dblprec
      avrgl4cum         = 0.0_dblprec
      avrgecum          = 0.0_dblprec
      avrgetcum         = 0.0_dblprec
      avrgelcum         = 0.0_dblprec
      avrge2cum         = 0.0_dblprec
      totene            = 0.0_dblprec
   end subroutine zero_cumulant_counters

   !---------------------------------------------------------------------------------
   ! SUBROUTINE: averages_allocations
   !> Allocation of the necessary arrays for the measurement of the averages, cumulants, and autocorrelation
   !---------------------------------------------------------------------------------
   subroutine averages_allocations(NA,NT,Nchmax,Mensemble,plotenergy,flag)

      implicit none

      integer, intent(in) :: NA         !< Number of atoms in one cell
      integer, intent(in) :: NT         !< Number of types of atoms
      integer, intent(in) :: flag       !< Allocate or deallocate (1/-1)
      integer, intent(in) :: Nchmax     !< Max number of chemical components on each site in cell
      integer, intent(in) :: Mensemble  !< Number of ensembles
      integer, intent(in) :: plotenergy !< Flag for plotting the energy

      !.. Local variables
      integer :: i_stat, i_all

      call allocate_projcumulants(nt,flag)

      if(flag>0) then
         bcount_avrg=1
         ! Allocations for the averages
         if (do_avrg=='Y'.or.do_proj_avrg=='Y' .or. do_proj_avrg=='A') then
            !allocate(mavg_buff(20,avrg_buff,Mensemble),stat=i_stat)
            allocate(mavg_buff(3,avrg_buff,Mensemble),stat=i_stat)
            call memocc(i_stat,product(shape(mavg_buff))*kind(mavg_buff),'mavg_buff','averages_allocations')
            mavg_buff=0.0_dblprec
            allocate(mavg_buff_proj(3,NA,avrg_buff,Mensemble),stat=i_stat)
            call memocc(i_stat,product(shape(mavg_buff_proj))*kind(mavg_buff_proj),'mavg_buff_proj','averages_allocations')
            mavg_buff_proj=0.0_dblprec
            allocate(mavg2_buff_proj(NA,avrg_buff,Mensemble),stat=i_stat)
            call memocc(i_stat,product(shape(mavg2_buff_proj))*kind(mavg2_buff_proj),'mavg2_buff_proj','allocate_measurements')
            mavg2_buff_proj=0.0_dblprec
         endif
         ! Index array should be allocated fpr all kind of possible measurements
         if (do_avrg=='Y'.or.do_proj_avrg=='Y' .or. do_proj_avrg=='A'.or.do_projch_avrg=='Y') then
            allocate(indxb_avrg(avrg_buff),stat=i_stat)
            call memocc(i_stat,product(shape(indxb_avrg))*kind(indxb_avrg),'indxb_avrg','averages_allocations')
            indxb_avrg=0
         endif
         ! Allocations for chemically projected averages
         if (do_projch_avrg=='Y') then
            allocate(mavg_buff_projch(4,Nchmax,avrg_buff,Mensemble),stat=i_stat)
            call memocc(i_stat,product(shape(mavg_buff_projch))*kind(mavg_buff_projch),'mavg_buff_projch','averages_allocations')
            mavg_buff_projch=0.0_dblprec
         end if
      else
         ! Deallocations for averages
         if (do_avrg=='Y'.or.do_proj_avrg=='Y' .or. do_proj_avrg=='A') then
            i_all=-product(shape(mavg_buff))*kind(mavg_buff)
            deallocate(mavg_buff,stat=i_stat)
            call memocc(i_stat,i_all,'mavg_buff','averages_allocations')
            i_all=-product(shape(mavg_buff_proj))*kind(mavg_buff_proj)
            deallocate(mavg_buff_proj,stat=i_stat)
            call memocc(i_stat,i_all,'mavg_buff_proj','averages_allocations')
            i_all=-product(shape(mavg2_buff_proj))*kind(mavg2_buff_proj)
            deallocate(mavg2_buff_proj,stat=i_stat)
            call memocc(i_stat,i_all,'mavg2_buff_proj','allocate_measurements')
         endif
         ! Index array should be allocated for all kind of possible measurements
         if (do_avrg=='Y'.or.do_proj_avrg=='Y' .or. do_proj_avrg=='A'.or.do_projch_avrg=='Y') then
            i_all=-product(shape(indxb_avrg))*kind(indxb_avrg)
            deallocate(indxb_avrg,stat=i_stat)
            call memocc(i_stat,i_all,'indxb_avrg','averages_allocations')
         endif
         ! Deallocations for chemically projected averages
         if(do_projch_avrg=='Y') then
            i_all=-product(shape(mavg_buff_projch))*kind(mavg_buff_projch)
            deallocate(mavg_buff_projch,stat=i_stat)
            call memocc(i_stat,i_all,'mavg_buff_projch','averages_allocations')
         end if

      endif

   end subroutine averages_allocations

   !---------------------------------------------------------------------------------
   ! SUBROUTINE: allocate_projcumulants
   !> Allocation of the necessary arrays for the projected cummulants 
   !---------------------------------------------------------------------------------
   subroutine allocate_projcumulants(nt,flag)

      implicit none

      integer, intent(in) :: NT  !< Number of atom types in the cell
      integer, intent(in) :: flag  !< Allocate or deallocate (1/-1)

      integer :: i_stat, i_all
      if (do_cumu_proj=='Y') then
         if (flag==1) then
            allocate(avrgmcum_proj(NT),stat=i_stat)
            call memocc(i_stat,product(shape(avrgmcum_proj))*kind(avrgmcum),'avrgmcum_proj','allocate_projcumulants')
            avrgmcum_proj=0.0_dblprec
            allocate(avrgm2cum_proj(NT),stat=i_stat)
            call memocc(i_stat,product(shape(avrgm2cum_proj))*kind(avrgm2cum_proj),'avrgm2cum_proj','allocate_projcumulants')
            avrgm2cum_proj=0.0_dblprec
            allocate(avrgm4cum_proj(NT),stat=i_stat)
            call memocc(i_stat,product(shape(avrgm4cum_proj))*kind(avrgm4cum_proj),'avrgm4cum_proj','allocate_projcumulants')
            avrgm4cum_proj=0.0_dblprec
         else
            i_all=-product(shape(avrgmcum_proj))*kind(avrgmcum_proj)
            deallocate(avrgmcum_proj,stat=i_stat)
            call memocc(i_stat,i_all,'avrgmcum_proj','allocate_projcumulants')
            i_all=-product(shape(avrgm2cum_proj))*kind(avrgm2cum_proj)
            deallocate(avrgm2cum_proj,stat=i_stat)
            call memocc(i_stat,i_all,'avrgm2cum_proj','allocate_projcumulants')
            i_all=-product(shape(avrgm4cum_proj))*kind(avrgm4cum_proj)
            deallocate(avrgm4cum_proj,stat=i_stat)
            call memocc(i_stat,i_all,'avrgm4cum_proj','allocate_projcumulants')
         end if
      end if
   end subroutine allocate_projcumulants

   !---------------------------------------------------------------------------------
   ! SUBROUTINE: buffer_avrg
   !> Buffer average magnetization
   !---------------------------------------------------------------------------------
   subroutine buffer_avrg(Natom,Mensemble,mstep,emom,emomM,bcount_avrg,delta_t,real_time_measure)
      !

      !.. Implicit declarations
      implicit none
      !
      integer, intent(in) :: mstep !< Current simulation step
      integer, intent(in) :: Natom !< Number of atoms in system
      integer, intent(in) :: Mensemble   !< Number of ensembles
      integer, intent(in) :: bcount_avrg !< Counter of buffer for averages
      real(dblprec), intent(in) :: delta_t               !< Current time step (used for real time measurements)
      character(len=1), intent(in) :: real_time_measure  !< Real time measurement flag
      real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: emom   !< Current unit moment vector
      real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: emomM  !< Current magnetic moment vector

      !.. Scalar variables
      integer :: i, k
      !real(dblprec) :: tm

      !.. Local arrays
      real(dblprec), dimension(3,Mensemble) ::  m

      !.. Executable statements
      m(:,:)=0.0_dblprec
      !.. Sum over moments
      do i=1, Natom
         do k=1,Mensemble
            m(:,k) = m(:,k) + emomM(:,i,k)
         end do
      end do

      !.. Save in buffer
      do k=1,Mensemble
         mavg_buff(1:3,bcount_avrg,k) = m(1:3,k)
      end do

      if (real_time_measure=='Y') then
         indxb_avrg(bcount_avrg)=mstep*delta_t
      else
         indxb_avrg(bcount_avrg)=mstep
      endif

   end subroutine buffer_avrg

   !---------------------------------------------------------------------------------
   ! SUBROUTINE: buffer_proj_avrg
   !> Buffer site projected average magnetizations
   !---------------------------------------------------------------------------------
   subroutine buffer_proj_avrg(Natom,Mensemble,NA,emomM,bcount_avrg,do_ralloy,      &
      Natom_full,asite_ch)

      !.. Implicit declarations
      implicit none

      integer, intent(in) :: NA  !< Number of atoms in one cell
      integer, intent(in) :: Natom !< Number of atoms in system
      integer, intent(in) :: Mensemble !< Number of ensembles
      integer, intent(in) :: do_ralloy  !< Random alloy simulation (0/1)
      integer, intent(in) :: Natom_full !< Number of atoms for full system (=Natom if not dilute)
      integer, intent(in) :: bcount_avrg !< Counter of buffer for averages
      integer, dimension(Natom_full), intent(in) :: asite_ch !< Actual site of atom for dilute system
      real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: emomM  !< Current magnetic moment vector

      !.. Locals
      real(dblprec), dimension(3,NA,Mensemble) ::  avg_mom !< Projected average magnetic moment

      !.. Scalar variables
      integer :: i,i_na,k

      !.. Executable statements
      avg_mom=0.0_dblprec

      !.. Sum over moments
      do k=1,Mensemble
         if(do_ralloy==0) then
            do i=1, Natom
               i_na=mod(i-1,NA)+1
               avg_mom(1,i_na,k) = avg_mom(1,i_na,k) + emomM(1,i,k)
               avg_mom(2,i_na,k) = avg_mom(2,i_na,k) + emomM(2,i,k)
               avg_mom(3,i_na,k) = avg_mom(3,i_na,k) + emomM(3,i,k)
            end do
         else
            do i=1, Natom
               i_na=asite_ch(i)
               avg_mom(1,i_na,k) = avg_mom(1,i_na,k) + emomM(1,i,k)
               avg_mom(2,i_na,k) = avg_mom(2,i_na,k) + emomM(2,i,k)
               avg_mom(3,i_na,k) = avg_mom(3,i_na,k) + emomM(3,i,k)
            end do
         end if

         !.. Save in buffer
         do i_na=1,NA
            mavg_buff_proj(1,i_na,bcount_avrg,k) = avg_mom(1,i_na,k)
            mavg_buff_proj(2,i_na,bcount_avrg,k) = avg_mom(2,i_na,k)
            mavg_buff_proj(3,i_na,bcount_avrg,k) = avg_mom(3,i_na,k)
            mavg2_buff_proj(i_na,bcount_avrg,k)  = (avg_mom(1,i_na,k)**2+avg_mom(3,i_na,k)**2+avg_mom(3,i_na,k)**2)
         end do

      end do
   end subroutine buffer_proj_avrg

   !---------------------------------------------------------------------------------
   ! SUBROUTINE: buffer_projch_avrg
   !> Buffer site and chemical projected average magnetizations
   !---------------------------------------------------------------------------------
   subroutine buffer_projch_avrg(Natom,Mensemble,Nchmax,achem_ch,emomM,bcount_avrg)

      !.. Implicit declarations
      implicit none

      integer, intent(in) :: Natom !< Number of atoms in system
      integer, intent(in) :: Nchmax !< Max number of chemical components on each site in cell
      integer, intent(in) :: Mensemble !< Number of ensembles
      integer, intent(in) :: bcount_avrg !< Counter of buffer for averages
      integer, dimension(Natom), intent(in) :: achem_ch !< Chemical type of atoms (reduced list) (achem_ch(i)=achtype(acellnumbrev(i)))
      real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: emomM  !< Current magnetic moment vector

      !.. Locals
      integer, dimension(Nchmax) :: counter
      real(dblprec), dimension(Nchmax) :: conc
      real(dblprec), dimension(4,Nchmax,Mensemble) ::  avg_mom !< Chemically projected average magnetic moment

      !.. Scalar variables
      integer :: i,ich,k

      !.. Executable statements
      avg_mom=0.0_dblprec

      !.. Sum over moments
      do k=1,Mensemble
         counter=0 ; conc=0.0_dblprec
         do i=1, Natom
            ich=achem_ch(i)
            counter(ich)=counter(ich)+1
            avg_mom(1,ich,k) = avg_mom(1,ich,k) + emomM(1,i,k)
            avg_mom(2,ich,k) = avg_mom(2,ich,k) + emomM(2,i,k)
            avg_mom(3,ich,k) = avg_mom(3,ich,k) + emomM(3,i,k)
            avg_mom(4,ich,k) = avg_mom(4,ich,k) + norm2(emomM(:,i,k)) 
         end do
         conc=(1.0_dblprec*counter/Natom)
         !.. Save in buffer
         do ich=1,Nchmax
            mavg_buff_projch(1,ich,bcount_avrg,k) = avg_mom(1,ich,k)/conc(ich)
            mavg_buff_projch(2,ich,bcount_avrg,k) = avg_mom(2,ich,k)/conc(ich)
            mavg_buff_projch(3,ich,bcount_avrg,k) = avg_mom(3,ich,k)/conc(ich)
            mavg_buff_projch(4,ich,bcount_avrg,k) = avg_mom(4,ich,k)/conc(ich)
         end do
      end do

   end subroutine buffer_projch_avrg

   !---------------------------------------------------------------------------------
   ! SUBROUTINE: prn_avrg
   !> Print instantaneous and cumulated average magnetizations
   !---------------------------------------------------------------------------------
   subroutine prn_avrg(Natom, Mensemble, simid,real_time_measure)

      !.. Implicit declarations
      implicit none

      integer, intent(in) :: Natom     !< Number of atoms in system
      integer, intent(in) :: Mensemble !< Number of ensembles
      character(len=8), intent(in) :: simid !< Name of simulation
      character(len=1), intent(in) :: real_time_measure !< Real time measurement flag

      !.. Scalar variables
      integer :: i,k
      character(len=30) :: filn
      real(dblprec), dimension(Mensemble) :: avrgm
      real(dblprec), dimension(3,Mensemble) :: avrgmv

      real(dblprec) :: avrgmm
      real(dblprec) :: avrgms
      real(dblprec), dimension(3) :: avrgmvm

      avrgmv=0.0_dblprec;avrgm=0.0_dblprec;avrgmvm=0.0_dblprec

      write (filn,'(''averages.'',a,''.out'')') trim(simid)
      open(ofileno, file=filn, position="append")

      ! Write header to output files for first iteration
      if(abs(indxb_avrg(1))==0.0e0_dblprec) then
         ! Averages
         if (real_time_measure=='Y') then
            write(ofileno,10007)"Time[s]","<M>_x","<M>_y","<M>_z","<M>","M_{stdv}"
         else
            write(ofileno,10006)"#Iter","<M>_x","<M>_y","<M>_z","<M>","M_{stdv}"
         endif
      end if

      !
      do i=1, bcount_avrg
         do k=1, Mensemble
            avrgmv(1:3,k) = mavg_buff(1:3,i,k)/Natom
         end do

         do k=1,Mensemble
            avrgm(k) = norm2(avrgmv(:,k))
         end do

         !mean magnetisations
         avrgmm   = 0.0_dblprec
         avrgmvm  = 0.0_dblprec
         !stdev magnetisation
         avrgms   = 0.0_dblprec

         do k=1,Mensemble
            avrgmvm(:) = avrgmvm(:) + avrgmv(:,k)
            avrgmm = avrgmm + avrgm(k)
            avrgms = avrgms + avrgm(k)**2
         end do

         avrgmvm=avrgmvm/Mensemble
         avrgmm=avrgmm/Mensemble

         !Do not use Mensemble-1 for numerical convenience
         avrgms=avrgms/Mensemble - avrgmm**2

         !Filter out negative errors
         if(avrgms<dbl_tolerance) then
            avrgms=0.0_dblprec
         else
            avrgms=sqrt(avrgms)
         end if

         if (real_time_measure=='Y') then
            write (ofileno,10005) indxb_avrg(i), avrgmvm(1), avrgmvm(2), avrgmvm(3), avrgmm, avrgms
         else
            write (ofileno,10004) int(indxb_avrg(i)), avrgmvm(1), avrgmvm(2), avrgmvm(3), avrgmm, avrgms
         endif
      end do

      close(ofileno)
      return

      write (*,*) "Error writing the averages file"

      10004 format (i8,6es16.8)
      10006 format (a8,6a16)
      10005 format (es12.4,6es16.8)
      10007 format (a12,6a16)

   end subroutine prn_avrg

   !---------------------------------------------------------------------------------
   ! SUBROUTINE: prn_proj_avrg
   !> Print site projected average magnetizations
   !---------------------------------------------------------------------------------
   subroutine prn_proj_avrg(Natom, Mensemble, NA, N1, N2, N3, simid, NT, atype,     &
         do_proj_avrg,real_time_measure)

      use Constants
      !.. Implicit declarations
      implicit none

      integer, intent(in) :: NA  !< Number of atoms in one cell
      integer, intent(in) :: N1  !< Number of cell repetitions in x direction
      integer, intent(in) :: N2  !< Number of cell repetitions in y direction
      integer, intent(in) :: N3  !< Number of cell repetitions in z direction
      integer, intent(in) :: NT  !< Number of types of atoms
      integer, intent(in) :: Natom     !< Number of atoms in system
      integer, intent(in) :: Mensemble !< Number of ensembles
      character(len=8), intent(in) :: simid          !< Name of simulation
      character(len=1), intent(in) :: do_proj_avrg   !< Measure projected averages (Y/A/N)
      character(len=1), intent(in) :: real_time_measure !< Real time measurement flag
      integer, dimension(Natom), intent(in) :: atype !< Type of atom

      !.. Scalar variables
      integer :: i,i_na,j,k,nproj
      character(len=30) :: filn

      real(dblprec), dimension(Mensemble,NA) :: avrgmx, avrgmy, avrgmz, avrgm, avrgm2, avrgm2t
      real(dblprec), dimension(NA) :: avrgmm, avrgms, avrgmxm, avrgmym, avrgmzm, avrgmm2

      i_na=-1;

      !.. Executable statements
      write (filn,'(''projavgs.'',a,''.out'')') trim(simid)
      open(ofileno, file=filn, position="append")

      if(abs(indxb_avrg(1))<=0.0e0_dblprec) then
         ! Averages
         if (real_time_measure=='Y') then
            write(ofileno,10007)"Time[s]","Proj","<M>","M_{stdv}","<M>_x","<M>_y","<M>_z"
         else
            write(ofileno,10006)"#Iter","Proj","<M>","M_{stdv}","<M>_x","<M>_y","<M>_z"
         endif
      end if

      do i=1, bcount_avrg
         if(do_proj_avrg=='Y') then
            do j=1, Mensemble
               do k=1, NT
                  avrgmx(j,k) =0.0_dblprec
                  avrgmy(j,k) =0.0_dblprec
                  avrgmz(j,k) =0.0_dblprec
                  avrgm2t(j,k)=0.0_dblprec

                  do i_na=1, NA
                     if(k==atype(i_na)) then
                        avrgmx(j,k) = avrgmx(j,k) + mavg_buff_proj(1,i_na,i,j)/(N1*N2*N3)
                        avrgmy(j,k) = avrgmy(j,k) + mavg_buff_proj(2,i_na,i,j)/(N1*N2*N3)
                        avrgmz(j,k) = avrgmz(j,k) + mavg_buff_proj(3,i_na,i,j)/(N1*N2*N3)
                        avrgm2t(j,k) = avrgm2t(j,k) + mavg2_buff_proj(i_na,i,j)/(N1*N2*N3)**2
                     end if
                  end do

                  !Type k projected mean magnetization for ensemble j
                  avrgm(j,k) = sqrt(avrgmx(j,k)**2+avrgmy(j,k)**2+avrgmz(j,k)**2)
                  avrgm2(j,k) = avrgm2t(j,k)
               end do
            end do
         end if
         if(do_proj_avrg=='A') then
            do j=1, Mensemble
               do k=1, NA
                  avrgmx(j,k)=mavg_buff_proj(1,k,i,j)/(N1*N2*N3)
                  avrgmy(j,k)=mavg_buff_proj(2,k,i,j)/(N1*N2*N3)
                  avrgmz(j,k)=mavg_buff_proj(3,k,i,j)/(N1*N2*N3)
                  avrgm2t(j,k) = mavg2_buff_proj(k,i,j)/(N1*N2*N3)**2
                  !Site k projected mean magnetization for ensemble j
                  avrgm(j,k) = sqrt(avrgmx(j,k)**2+avrgmy(j,k)**2+avrgmz(j,k)**2)
                  avrgm2(j,k) = avrgm2t(j,k)
               end do
            end do
         end if

         !Start with calculation of mean over the ensembles
         if(do_proj_avrg=='Y') then
            nproj=NT
         else if(do_proj_avrg=='A') then
            nproj=NA
         else
            nproj=NA
         end if

         do k=1, nproj

            !mean magnetisation over the ensembles
            avrgmm(k)  =0.0_dblprec
            avrgmm2(k) =0.0_dblprec

            !standard deviation of the mean magnetisation over the ensembles
            avrgms(k)=0.0_dblprec

            avrgmxm(k)=0.0_dblprec
            avrgmym(k)=0.0_dblprec
            avrgmzm(k)=0.0_dblprec
            do j=1, Mensemble
               avrgmxm(k) = avrgmxm(k) + avrgmx(j,k)
               avrgmym(k) = avrgmym(k) + avrgmy(j,k)
               avrgmzm(k) = avrgmzm(k) + avrgmz(j,k)
               avrgmm(k) = avrgmm(k) + avrgm(j,k)
               avrgms(k) = avrgms(k) + avrgm(j,k)**2
               avrgmm2(k) = avrgmm2(k) + avrgm2(j,k)/nproj
            end do
            avrgmxm(k)=avrgmxm(k)/Mensemble
            avrgmym(k)=avrgmym(k)/Mensemble
            avrgmzm(k)=avrgmzm(k)/Mensemble
            avrgmm(k)=avrgmm(k)/Mensemble
            avrgms(k)=avrgms(k)/Mensemble
            avrgms(k)=avrgms(k) - avrgmm(k)**2
            avrgmm2(k)=avrgmm2(k)/Mensemble

            !Filter out negative errors
            if(avrgms(k)<0) then
               avrgms(k)=0
            else
               avrgms(k)=sqrt(avrgms(k))
            end if

            !Write results to file
            if (real_time_measure=='Y') then
               write (ofileno,10005) indxb_avrg(i), k, avrgmm(k), avrgms(k),&
                  avrgmxm(k), avrgmym(k), avrgmzm(k)
            else
               write (ofileno,10004) int(indxb_avrg(i)), k, avrgmm(k), avrgms(k),&
                  avrgmxm(k), avrgmym(k), avrgmzm(k)
            endif
         end do
      end do
      close(ofileno)
      return

      write (*,*) "Error writing the projected averages file"
      10004 format (i8,i8,5es16.8)
      10006 format (a8,a,5a16)
      10005 format (es16.4,i8,5es16.8)
      10007 format (a16,a,5a16)

   end subroutine prn_proj_avrg

   !---------------------------------------------------------------------------------
   ! SUBROUTINE: prn_projch_avrg
   !> Print chemical projected average magnetizations
   !---------------------------------------------------------------------------------
   subroutine prn_projch_avrg(Mensemble, N1, N2, N3, simid, Nchmax,real_time_measure)
      !
      !.. Implicit declarations
      implicit none

      integer, intent(in) :: N1  !< Number of cell repetitions in x direction
      integer, intent(in) :: N2  !< Number of cell repetitions in y direction
      integer, intent(in) :: N3  !< Number of cell repetitions in z direction
      integer, intent(in) :: Nchmax    !< Max number of chemical components on each site in cell
      integer, intent(in) :: Mensemble !< Number of ensembles
      character(len=8), intent(in) :: simid !< Name of simulation
      character(len=1), intent(in) :: real_time_measure !< Real time measurement flag

      !.. Scalar variables
      integer :: i,j,k
      character(len=30) :: filn
      real(dblprec), dimension(Nchmax) :: avrgmm, avrgms, avrgmxm, avrgmym, avrgmzm, avrgmagm
      real(dblprec), dimension(Mensemble,Nchmax) :: avrgmx, avrgmy, avrgmz, avrgm, avrgmag

      avrgmx=0.0_dblprec;avrgmy=0.0_dblprec;avrgmz=0.0_dblprec;avrgm=0.0_dblprec
      avrgmag=0.0_dblprec
      !.. Executable statements
      write (filn,'(''projchavgs.'',a,''.out'')') trim(simid)
      open(ofileno, file=filn, position="append")

      if(abs(indxb_avrg(1))<=0.0e0_dblprec) then
         ! Averages
         if (real_time_measure=='Y') then
            write(ofileno,10007)"#Time[s]","Proj","<M>","M_{stdv}","<M>_x","<M>_y","<M>_z","\sum<M>"
         else
            write(ofileno,10006)"# Iter.","Proj","<M>","M_{stdv}","<M>_x","<M>_y","<M>_z","\sum<M>"
         endif
      end if

      do i=1, bcount_avrg
         do j=1, Mensemble
            do k=1, Nchmax
               avrgmx(j,k)    = mavg_buff_projch(1,k,i,j)/(N1*N2*N3)
               avrgmy(j,k)    = mavg_buff_projch(2,k,i,j)/(N1*N2*N3)
               avrgmz(j,k)    = mavg_buff_projch(3,k,i,j)/(N1*N2*N3)
               avrgmag(j,k)   = mavg_buff_projch(4,k,i,j)/(N1*N2*N3)
               !Calculate mean magnetisation sublatice for each ensemble
               ! average in ensemble j of sublattice k
               avrgm(j,k) = sqrt(avrgmx(j,k)**2+avrgmy(j,k)**2+avrgmz(j,k)**2)
            end do
         end do

         !Start with calculation of mean over the ensembles
         do k=1, Nchmax

            !mean magnetisation over the ensembles
            avrgmm(k)=0.0_dblprec

            !standard deviation of the mean magnetisation over the ensembles
            avrgms(k)   = 0.0_dblprec
            avrgmxm(k)  = 0.0_dblprec
            avrgmym(k)  = 0.0_dblprec
            avrgmzm(k)  = 0.0_dblprec
            avrgmagm(k) = 0.0_dblprec

            do j=1, Mensemble
               avrgmxm(k) = avrgmxm(k) + avrgmx(j,k)
               avrgmym(k) = avrgmym(k) + avrgmy(j,k)
               avrgmzm(k) = avrgmzm(k) + avrgmz(j,k)
               avrgmm(k) = avrgmm(k) + avrgm(j,k)
               avrgmagm(k) = avrgmagm(k) + avrgmag(j,k)
               avrgms(k) = avrgms(k) + avrgm(j,k)**2
            end do

            avrgmxm(k)=avrgmxm(k)/Mensemble
            avrgmym(k)=avrgmym(k)/Mensemble
            avrgmzm(k)=avrgmzm(k)/Mensemble
            avrgmm(k)=avrgmm(k)/Mensemble
            avrgms(k)=avrgms(k)/Mensemble
            avrgmagm(k)=avrgmagm(k)/Mensemble
            avrgms(k)=avrgms(k) - avrgmm(k)**2

            !Filter out negative errors
            if(avrgms(k)<0) then
               avrgms(k)=0
            else
               avrgms(k)=sqrt(avrgms(k))
            end if

            !Write results to file
            if (real_time_measure=='Y') then
               write (ofileno,10005) indxb_avrg(i), k, avrgmm(k), avrgms(k),&
                  avrgmxm(k), avrgmym(k), avrgmzm(k),avrgmagm(k)
            else
               write (ofileno,10004) int(indxb_avrg(i)), k, avrgmm(k), avrgms(k),&
                  avrgmxm(k), avrgmym(k), avrgmzm(k),avrgmagm(k)
            endif
         end do
      end do
      close(ofileno)
      return

      write (*,*) "Error writing the chemical projected averages file"
      10004 format (i8,i8,6es16.8)
      10006 format (a8,a,6a16)
      10005 format (es16.4,i8,6es16.8)
      10007 format (a16,a,6a16)

   end subroutine prn_projch_avrg

   !---------------------------------------------------------------------------------
   ! SUBROUTINE: calc_and_print_cumulant
   !> Calculate and print Binder cumulant, spin susceptibility, and specific heat
   !---------------------------------------------------------------------------------
   subroutine calc_and_print_cumulant(Natom,Mensemble,emomM,simid,Temp,temprescale, &
      temprescalegrad,plotenergy,cumu_buff,do_prn_flag)
      !
      use Constants
      use prn_topology
      use Topology, only : chi_cavg, kappa_cavg, n_chi_cavg, kappa_csum
      use Polarization, only : do_chiral

      !.. Implicit declarations
      implicit none

      integer, intent(in) :: Natom        !< Number of atoms in system
      integer, intent(in) :: Mensemble    !< Number of ensembles
      integer, intent(in) :: cumu_buff    !< Buffer size for Binder cumulant
      integer, intent(in) :: plotenergy   !< Calculate and plot energy (0/1)
      real(dblprec), intent(in) :: Temp   !< Temperature
      real(dblprec), intent(in) :: temprescale   !< Temperature rescaling if QHB
      real(dblprec), intent(in) :: temprescalegrad   !< Temperature rescaling gradient if QHB
      logical, intent(in) :: do_prn_flag  !< Print Yes or No (No to only run calc)
      character(len=8), intent(in) :: simid  !< Name of simulation
      real(dblprec), intent(in), dimension(3,Natom, Mensemble) :: emomM  !< Current magnetic moment vector

      !.. Scalar variables
      integer :: i,k
      character(len=30) :: filn
      real(dblprec) :: avrgme, avrgm2, avrgm4
      real(dblprec) :: cumulant, avrgmt, avrgmt2, avrgmt4
      real(dblprec) :: pmsusc, cv, avrgen, avrgen2, avrgent,avrgenl
      real(dblprec),dimension(3) :: m

      cumulant=0_dblprec;avrgmt=0_dblprec;avrgmt2=0_dblprec;avrgmt4=0_dblprec
      pmsusc=0_dblprec;cv=0_dblprec

      do k=1,Mensemble
         m=0.0_dblprec

#if _OPENMP >= 201307
         !$omp parallel do default(shared) private(i) reduction(+:m) schedule(static)
#endif
         do i=1, Natom
            m(:) = m(:) + emomM(:,i,k)
         end do
#if _OPENMP >= 201307
         !$omp end parallel do
#endif
         avrgme = norm2(m)/Natom
         avrgm2 = avrgme**2
         avrgm4 = avrgm2**2

         cumuw = cumuw + 1.0_dblprec

         avrgmt  = (avrgmcum*cumutotw+avrgme*cumuw )/(cumutotw+cumuw)
         avrgmt2 = (avrgm2cum*cumutotw+avrgm2*cumuw)/(cumutotw+cumuw)
         avrgmt4 = (avrgm4cum*cumutotw+avrgm4*cumuw)/(cumutotw+cumuw)

         !!!avrgmt = (Navrgcum*avrgmcum+avrgme)/(Navrgcum+1)
         !!!avrgmt2 = (Navrgcum*avrgm2cum+avrgm2)/(Navrgcum+1)
         !!!avrgmt4 = (Navrgcum*avrgm4cum+avrgm4)/(Navrgcum+1)

         cumulant = 1-(avrgmt4/3/avrgmt2**2)
         avrgmcum = avrgmt
         avrgm2cum = avrgmt2
         avrgm4cum = avrgmt4

         !Specific heat (cv) and susceptibility (pmsusc)
         if(Temp>0.0_dblprec) then
            pmsusc = (avrgmt2-avrgmt**2) * mub**2 * Natom / (k_bolt**2) / Temp     ! units of k_B
         else
            !For T=0, not the proper susceptibility
            pmsusc = (avrgmt2-avrgmt**2) * Natom * mub**2 / k_bolt
         end if

         if (plotenergy>0.and.Temp>0.0_dblprec) then
            !!! avrgen   = (Navrgcum*avrgecum+ene%energy(k))/(Navrgcum+1)
            !!! avrgen2  = (Navrgcum*avrge2cum+ene%energy(k)**2)/(Navrgcum+1)
            !!! avrgent  = (Navrgcum*avrgetcum+ene%ene_xc(k))/(Navrgcum+1)
            !!! avrgenl  = (Navrgcum*avrgelcum+ene%ene_lsf(k))/(Navrgcum+1)

            avrgen   = (avrgecum*cumutotw+ene%energy(k)*cumuw)/(cumutotw+cumuw)
            avrgen2  = (avrge2cum*cumutotw+ene%energy(k)**2*cumuw)/(cumutotw+cumuw)
            avrgent  = (avrgetcum*cumutotw+ene%ene_xc(k)*cumuw)/(cumutotw+cumuw)
            avrgenl  = (avrgelcum*cumutotw+ene%ene_lsf(k)*cumuw)/(cumutotw+cumuw)

            cv  =  ( avrgen2  -avrgen**2  ) * mry**2 * Natom  / (k_bolt**2) / Temp**2 * (temprescalegrad*Temp+temprescale) / temprescale**2 !* Natom   ! units of k_B
            avrgecum    = avrgen
            avrge2cum   = avrgen2
            avrgetcum   = avrgent
            avrgelcum   = avrgenl

         else if(plotenergy>0) then
            !!! avrgen   = (Navrgcum*avrgecum+ene%energy(k))/(Navrgcum+1)
            !!! avrgen2  = (Navrgcum*avrge2cum+ene%energy(k)**2)/(Navrgcum+1)
            !!! avrgent  = (Navrgcum*avrgetcum+ene%ene_xc(k))/(Navrgcum+1)
            !!! avrgenl  = (Navrgcum*avrgelcum+ene%ene_lsf(k))/(Navrgcum+1)
            avrgen   = (avrgecum*cumutotw+ene%energy(k)*cumuw)/(cumutotw+cumuw)
            avrgen2  = (avrge2cum*cumutotw+ene%energy(k)**2*cumuw)/(cumutotw+cumuw)
            avrgent  = (avrgetcum*cumutotw+ene%ene_xc(k)*cumuw)/(cumutotw+cumuw)
            avrgenl  = (avrgelcum*cumutotw+ene%ene_lsf(k)*cumuw)/(cumutotw+cumuw)

            cv = 0.0_dblprec
            avrgecum    = avrgen
            avrge2cum   = avrgen2
            avrgetcum   = avrgent
            avrgelcum   = avrgenl
         else
            cv = 0.0_dblprec
            avrgen      = 0.0_dblprec
            avrgen2     = 0.0_dblprec
            avrgent     = 0.0_dblprec
            avrgenl     = 0.0_dblprec
            avrgecum    = 0.0_dblprec
            avrge2cum   = 0.0_dblprec
            avrgetcum   = 0.0_dblprec
            avrgelcum   = 0.0_dblprec
         end if

         Navrgcum = Navrgcum+1
         cumutotw = cumutotw + cumuw
      end do

      ! Write output to file
      if(do_prn_flag .eqv. .true.) then
         if(mod(Navrgcum/Mensemble-1,cumu_buff)==0) then
            ! Open file
            write (filn,'(''cumulants.'',a,''.out'')') trim(simid)
            open(ofileno, file=filn, position="append")

            ! Write header the first time
            if(Navrgcum/Mensemble==1) then
               write(ofileno,10005)"#Iter","<M>","<M^2>","<M^4>","U_{Binder}",&
                  "\chi","C_v(tot)","<E>","<E_{exc}>","<E_{lsf}>"
            end if
            !if (real_time_measure=='Y') then
            !   write (ofileno,10014) Navrgcum/Mensemble,avrgmt,avrgmt2,avrgmt4,        &
            !      cumulant,pmsusc,cv,avrgen,avrgent,avrgenl
            !else
            write (ofileno,10004) Navrgcum/Mensemble,avrgmt,avrgmt2,avrgmt4,        &
               cumulant,pmsusc,cv,avrgen,avrgent,avrgenl
            !end if
            close(ofileno)

            ! Open file
            write (filn,'(''cumulants.'',a8,''.json'')') simid
            open(ofileno, file=filn)

            write(ofileno,'(a)') '{'
            write(ofileno,'(a,f16.8,a)') '    "temperature"     : ', Temp,' ,'
            write(ofileno,'(a,f16.8,a)') '    "magnetization"   : ', avrgmt,' ,'
            write(ofileno,'(a,f16.8,a)') '    "binder_cumulant" : ', cumulant,' ,'
            if (avrgen.ne.0.0_dblprec) then
               write(ofileno,'(a,f16.8,a)') '    "energy"          : ', avrgen, ' ,'
            else
               write(ofileno,'(a)') '    "energy"          :  null ,'
            end if
            if (skyno=='T'.or.skyno=='Y') then
               write(ofileno,'(a,f16.8,a)') '    "skyrmion_num"    : ', sk_avrg,' ,'
               write(ofileno,'(a,f16.8,a)') '    "skyrmion_std"    : ', sqrt(sk_var),' ,'
            else
               write(ofileno,'(a)') '    "skyrmion_num"    :  null ,'
               write(ofileno,'(a)') '    "skyrmion_std"    :  null ,'
            end if
            if (do_chiral == 'Y') then
               write(ofileno,'(a,f16.8,a)') '    "scalar_chirality"    : ', chi_cavg,' ,'
               kappa_cavg = kappa_csum / n_chi_cavg
               ! print *, 'kappa_csum', kappa_csum
               ! print *, 'n_chi_cavg', n_chi_cavg
               ! print *, 'kappa_cavg', kappa_cavg
               write(ofileno,'(a,a,f16.8,a,f16.8,a,f16.8,a)') '    "vector_chirality"    : ', '[ ', &
                  kappa_cavg(1), ', ', kappa_cavg(2), ', ', kappa_cavg(3), '],'
            end if
            write(ofileno,'(a,f16.8,a)') '    "susceptibility"  : ', pmsusc,' ,'
            write(ofileno,'(a,f16.8,a)') '    "specific_heat"   : ', cv
            write(ofileno,'(a)') '}'
            close(ofileno)

         end if
      end if

      binderc=cumulant

      return

      write(*,*) 'Error writing the cumulant file'
      10004 format (i8,10es16.8)
      10005 format (a8,9a16)
      10014 format (es16.8,10es16.8)

   end subroutine calc_and_print_cumulant

   !---------------------------------------------------------------------------------
   ! SUBROUTINE: calc_and_print_afm_cumulant
   !> Calculate and print Binder cumulant, spin susceptibility, and specific heat
   !! @todo Read AFM vector from file
   !---------------------------------------------------------------------------------
   subroutine calc_and_print_afm_cumulant(Natom,Mensemble,emomM,simid,Temp,temprescale,&
      temprescalegrad,plotenergy, cumu_buff,NA, do_ralloy, Natom_full, asite_ch)
      !
      use Constants

      !.. Implicit declarations
      implicit none

      integer, intent(in) :: NA         !< Number of atoms in one cell
      integer, intent(in) :: Natom      !< Number of atoms in system
      integer, intent(in) :: cumu_buff  !< Buffer size for Binder cumulant
      integer, intent(in) :: do_ralloy  !< Random alloy simulation (0/1)
      integer, intent(in) :: Mensemble  !< Number of ensembles
      integer, intent(in) :: Natom_full !< Number of atoms for full system (=Natom if not dilute)
      integer, intent(in) :: plotenergy !< Calculate and plot energy (0/1)
      real(dblprec), intent(in) :: Temp !< Temperature
      real(dblprec), intent(in) :: temprescale !< Temperature rescaling if QHB
      real(dblprec), intent(in) :: temprescalegrad !< Temperature rescaling gradient if QHB
      character(len=8), intent(in) :: simid !< Name of simulation
      integer, dimension(Natom_full), intent(in) :: asite_ch !< Actual site of atom for dilute system
      real(dblprec), dimension(3,Natom, Mensemble), intent(in) :: emomM  !< Current magnetic moment vector

      !.. Local arrays
      real(dblprec), dimension(4) :: lvec !< Antiferromagnetic order parameter mask
      real(dblprec), dimension(3,NA,Mensemble) ::  avg_mom !< Projected average magnetic moment

      !.. Scalar variables
      integer :: i,i_na,k
      character(len=30) :: filn
      real(dblprec) :: avrgle, avrgl2, avrgl4
      real(dblprec) :: lx, ly, lz
      real(dblprec) :: cumulant, avrglt, avrglt2, avrglt4
      real(dblprec) :: pmsusc, cv, avrgen, avrgen2

      !.. Executable statements
      avg_mom=0.0_dblprec
      lvec(1) =  1
      lvec(2) =  1
      lvec(3) = -1
      lvec(4) = -1

      !.. Sum over moments
      do k=1,Mensemble
         if(do_ralloy==0) then
            do i=1, Natom
               i_na=mod(i-1,NA)+1
               avg_mom(1,i_na,k) = avg_mom(1,i_na,k) + emomM(1,i,k)
               avg_mom(2,i_na,k) = avg_mom(2,i_na,k) + emomM(2,i,k)
               avg_mom(3,i_na,k) = avg_mom(3,i_na,k) + emomM(3,i,k)
            end do
         else
            do i=1, Natom
               i_na=asite_ch(i)
               avg_mom(1,i_na,k) = avg_mom(1,i_na,k) + emomM(1,i,k)
               avg_mom(2,i_na,k) = avg_mom(2,i_na,k) + emomM(2,i,k)
               avg_mom(3,i_na,k) = avg_mom(3,i_na,k) + emomM(3,i,k)
            end do
         end if
         lx=0.0_dblprec;ly=0.0_dblprec;lz=0.0_dblprec
         do i=1,NA
            lx = lx + lvec(i)*avg_mom(1,i,k)
            ly = ly + lvec(i)*avg_mom(2,i,k)
            lz = lz + lvec(i)*avg_mom(3,i,k)
         end do
         avrgle = sqrt(lx**2+ly**2+lz**2)/Natom
         avrgl2 = avrgle**2
         avrgl4 = avrgl2**2

         avrglt = (Navrgcum*avrglcum+avrgle)/(Navrgcum+1)
         avrglt2 = (Navrgcum*avrgl2cum+avrgl2)/(Navrgcum+1)
         avrglt4 = (Navrgcum*avrgl4cum+avrgl4)/(Navrgcum+1)

         cumulant = 1-(avrglt4/3/avrglt2**2)
         avrglcum = avrglt
         avrgl2cum = avrglt2
         avrgl4cum = avrglt4

         !Specific heat (cv) and susceptibility (pmsusc)
         if(Temp>0.0_dblprec) then
            pmsusc = (avrglt2-avrglt**2) * mub**2 * Natom / (k_bolt**2) / Temp       ! units of k_B
         else
            !For T=0, not the proper susceptibility
            pmsusc = (avrglt2-avrglt**2) * Natom * mub**2 / k_bolt
         end if

         if (plotenergy>0.and.Temp>0.0_dblprec) then
            avrgen = (Navrgcum*avrgecum+ene%energy(k))/(Navrgcum+1)
            avrgen2 = (Navrgcum*avrge2cum+ene%energy(k)**2)/(Navrgcum+1)
            cv = ( avrgen2 -avrgen**2 ) * mry**2 * Natom / (k_bolt**2) / Temp**2 * (temprescalegrad*Temp+temprescale) / temprescale**2 ! units of k_B
            avrgecum = avrgen
            avrge2cum = avrgen2
         else if(plotenergy>0) then
            avrgen = (Navrgcum*avrgecum+ene%energy(k))/(Navrgcum+1)
            avrgen2 = (Navrgcum*avrge2cum+ene%energy(k)**2)/(Navrgcum+1)
            cv = 0.0_dblprec
            avrgecum = avrgen
            avrge2cum = avrgen2
         else
            cv          = 0.0_dblprec
            avrgen      = 0.0_dblprec
            avrgen2     = 0.0_dblprec
            avrgecum    = 0.0_dblprec
            avrge2cum   = 0.0_dblprec
         end if
         Navrgcum = Navrgcum+1
      end do

      ! Write output to file
      if(mod(Navrgcum/Mensemble-1,cumu_buff)==0) then

         ! Open file
         write (filn,'(''afmcumulants.'',a,''.out'')') trim(simid)
         open(ofileno, file=filn, position="append")

         ! Write header the first time
         if(Navrgcum/Mensemble==1) then
            write (ofileno,10005)"#Iter.","<L>","<L^2>","<L^4>","U_{Binder}^L",&
               "Susc.","Cv"
         end if
         write (ofileno,10004) Navrgcum/Mensemble,avrglt,avrglt2,avrglt4,cumulant,  &
            pmsusc,cv
         close(ofileno)
      end if
      binderc=cumulant

      return

      write(*,*) 'Error writing the AFM cumulant file'
      10004 format (i8,6es16.8)
      10005 format (a8,6a16)

   end subroutine calc_and_print_afm_cumulant

   !---------------------------------------------------------------------------------
   ! SUBROUTINE: calc_and_print_cumulant_proj
   !> Calculate and print Binder cumulant, spin susceptibility, and specific heat
   !---------------------------------------------------------------------------------
   subroutine calc_and_print_cumulant_proj(Natom,Mensemble,emomM,simid,T,cumu_buff, &
      NT,atype)
      !
      use Constants

      !.. Implicit declarations
      implicit none

      integer, intent(in) :: NT           !< Number of atom types in system
      integer, intent(in) :: Natom        !< Number of atoms in system
      integer, intent(in) :: Mensemble    !< Number of ensembles
      integer, intent(in) :: cumu_buff    !< Buffer size for Binder cumulant
      real(dblprec),intent(in) :: T       !< Temperature
      character(len=8), intent(in) :: simid !< Name of simulation
      integer, dimension(Natom), intent(in) :: atype !< Type of atom
      real(dblprec), intent(in), dimension(3,Natom, Mensemble) :: emomM  !< Current magnetic moment vector

      !.. Scalar variables
      integer :: i,k,it
      character(len=30) :: filn
      integer,dimension(nt) :: ncounter
      real(dblprec),dimension(NT) :: mx, my, mz
      real(dblprec),dimension(NT) :: avrgme, avrgm2, avrgm4
      real(dblprec),dimension(NT) :: cumulant, avrgmt, avrgmt2, avrgmt4,pmsusc

      avrgme=0.0_dblprec;avrgm2=0.0_dblprec;avrgm4=0.0_dblprec
      mx=0.0_dblprec;my=0.0_dblprec;mz=0.0_dblprec
      cumulant=0.0_dblprec;avrgmt=0.0_dblprec;avrgmt2=0.0_dblprec;avrgmt4=0.0_dblprec
      pmsusc=0.0_dblprec
      ncounter=0
      do it=1,nt
         do i=1,Natom
            if (it==atype(i)) then
               ncounter(it)=ncounter(it)+1
            end if
         enddo
      enddo

      do k=1,Mensemble
         mx=0.0_dblprec;my=0.0_dblprec;mz=0.0_dblprec
         do i=1, Natom
            it=atype(i)
            mx(it) = mx(it) + emomM(1,i,k)
            my(it) = my(it) + emomM(2,i,k)
            mz(it) = mz(it) + emomM(3,i,k)
         end do

         do it=1,NT

            avrgme(it) = sqrt(mx(it)**2+my(it)**2+mz(it)**2)/ncounter(it)
            avrgm2(it) = avrgme(it)**2
            avrgm4(it) = avrgm2(it)**2

            avrgmt(it) = (Navrgcum_proj*avrgmcum_proj(it)+avrgme(it))/(Navrgcum_proj+1)
            avrgmt2(it) = (Navrgcum_proj*avrgm2cum_proj(it)+avrgm2(it))/(Navrgcum_proj+1)
            avrgmt4(it) = (Navrgcum_proj*avrgm4cum_proj(it)+avrgm4(it))/(Navrgcum_proj+1)

            cumulant(it) = 1-(avrgmt4(it)/3/avrgmt2(it)**2)
            avrgmcum_proj(it) = avrgmt(it)
            avrgm2cum_proj(it) = avrgmt2(it)
            avrgm4cum_proj(it) = avrgmt4(it)

            !Specific heat (cv) and susceptibility (pmsusc)
            if(T>0.0_dblprec) then
               pmsusc(it) = (avrgmt2(it)-avrgmt(it)**2) * mub**2 * ncounter(it) / (k_bolt**2) / T     ! units of k_B
            else
               !For T=0, not the proper susceptibility
               pmsusc(it) = (avrgmt2(it)-avrgmt(it)**2) * ncounter(it) * mub**2 / k_bolt
            end if
         end do
         Navrgcum_proj = Navrgcum_proj+1
      enddo

      ! Write output to file
      if(mod(Navrgcum_proj/Mensemble-1,cumu_buff)==0) then

         ! Open file
         write (filn,'(''projcumulants.'',a,''.out'')') trim(simid)
         open(ofileno, file=filn, position="append")

         ! Write header the first time
         if(Navrgcum_proj/Mensemble==1) then
            write (ofileno,10005) "# ter.","Type","<M>","<M^2>","<M^4>",     &
               "U_{Binder}","\chi"
         end if
         do it=1,NT
            write (ofileno,10004) Navrgcum_proj/Mensemble, it, avrgmt(it), avrgmt2(it), avrgmt4(it), cumulant(it), pmsusc(it)
         enddo
         close(ofileno)
      end if

      10004 format (i8,i6,5es16.8)
      10005 format (a8,a6,5a16)

   end subroutine calc_and_print_cumulant_proj

   !---------------------------------------------------------------------------
   !> @brief
   !> Read input parameters.
   !
   !> @author
   !> Anders Bergman
   !---------------------------------------------------------------------------
   subroutine read_parameters_averages(ifile)
      use FileParser

      implicit none

      ! ... Formal Arguments ...
      integer, intent(in) :: ifile   !< File to read from
      !
      ! ... Local Variables ...
      character(len=50) :: keyword
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

            case('do_avrg')
               read(ifile,*,iostat=i_err) do_avrg
               if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err

            case('do_proj_avrg')
               read(ifile,*,iostat=i_err) do_proj_avrg
               if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err

            case('do_projch_avrg')
               read(ifile,*,iostat=i_err) do_projch_avrg
               if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err

            case('avrg_step')
               read(ifile,*,iostat=i_err) avrg_step
               if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err

            case('avrg_buff')
               read(ifile,*,iostat=i_err) avrg_buff
               if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err

            case('do_cumu')
               read(ifile,*,iostat=i_err) do_cumu
               if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err

            case('do_cumu_proj')
               read(ifile,*,iostat=i_err) do_cumu_proj
               if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err

            case('cumu_step')
               read(ifile,*,iostat=i_err) cumu_step
               if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err

            case('cumu_buff')
               read(ifile,*,iostat=i_err) cumu_buff
               if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err

            case('do_cuda_avrg')
               read(ifile,*,iostat=i_err) do_cuda_avrg
               if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err

            case('do_cuda_cumu')
               read(ifile,*,iostat=i_err) do_cuda_cumu
               if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err

            !!! case('do_autocorr')
            !!!    read(ifile,*,iostat=i_err) do_autocorr
            !!!    if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err

            !!! case('acfile')
            !!!    read(ifile,'(a)',iostat=i_err) cache
            !!!    if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err
            !!!    twfile=trim(adjustl(cache))
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
   end subroutine read_parameters_averages

end module prn_averages
