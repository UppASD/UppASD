!-------------------------------------------------------------------------------
!  MODULE: Measurements
!> Data and routines for measuring averages and trajectories
!> @author
!> Anders Bergman, Lars Bergqvist, Johan Hellsvik, Jonathan Chico
!> @copyright
!> GNU Public License.
!-------------------------------------------------------------------------------
module Measurements
   use Parameters
   use Profiling
   use Math_functions, only : f_logstep
   !
   implicit none
   !
   private

   !public subroutines
   public :: measure, flush_measurements, do_measurements, allocate_measurementdata, calc_mavrg

contains

   !-----------------------------------------------------------------------------
   !  SUBROUTINE: do_measurements
   !> @brief
   !> Helper method for determining if measurements are neccessary
   !> used during gpu acceleration
   !>
   !> @author
   !> Niklas Fejes
   !-----------------------------------------------------------------------------
   subroutine do_measurements(mstep,do_avrg,do_tottraj,avrg_step,ntraj,tottraj_step,&
      traj_step,do_cumu,cumu_step,logsamp,do_copy,do_gpu_measurements)
      !
      implicit none
      !
      integer, intent(in)          :: mstep              !< Current simulation step
      character(len=1),intent(in)  :: do_avrg            !< Measure average magnetization (Y/N)
      character(len=1),intent(in)  :: do_tottraj         !< Measure magnetic moments
      integer, intent(in)          :: avrg_step          !< Interval for sampling average magnetization
      integer, intent(in)          :: ntraj              !< Number of trajectories to sample
      integer, intent(in)          :: tottraj_step       !< Interval for sampling magnetic moments
      integer, intent(in), dimension(ntraj) :: traj_step !< Interval for sampling individual trajectories
      character(len=1), intent(in) :: do_cumu            !< Measure Binder cumulant, susceptibility, and specific heat(Y/N)
      integer, intent(in)          :: cumu_step          !< Interval for sampling Binder cumulant
      character(len=1)             :: logsamp            !< Sample measurements logarithmically (Y/N)
      integer, intent(out)         :: do_copy            !< Flag if moment must be copied
      character(len=1),intent(in)  :: do_gpu_measurements  !< Measurements (Y/N) with CUDA

      integer :: i
      integer :: sstep

      sstep = f_logstep(mstep,logsamp)

      ! Averages
      if (do_avrg=='Y') then
         if ( mod(sstep-1,avrg_step)==0) then
            do_copy = 1
            return
         endif
      endif

      !  Total trajectory
      if (do_tottraj=='Y') then
         if (mod(sstep-1,tottraj_step)==0) then
            do_copy = 1
            return
         endif
      endif

      !  Atom trajectory
      if (ntraj>0) then
         do i=1,ntraj
            if (mod(sstep-1,traj_step(i))==0) then
               do_copy = 1
               return
            end if
         enddo
      endif

      ! Binder cumulant, susceptibility, and specific heat
      if(do_cumu=='Y') then
         if(mod(sstep,cumu_step)==0) then
            do_copy = 1
            return
         end if
      elseif(do_cumu=='A') then
         if(mod(sstep,cumu_step)==0) then
            do_copy = 1
            return
         end if
      end if

      ! If not returned yet, don't copy
      do_copy = 0

   end subroutine do_measurements

   !---------------------------------------------------------------------------------
   !  SUBROUTINE: measure
   !> @brief
   !> Wrapper routine for sampling and printing of several observables
   !> @detail 
   !-----------------------------------------------------------------------------
   subroutine measure(Natom,Mensemble,NT,NA,nHam,N1,N2,N3,simid,mstep,emom,emomM,   &
      mmom,Nchmax,do_ralloy,Natom_full,asite_ch,achem_ch,atype,plotenergy,Temp,     &
      temprescale,temprescalegrad,real_time_measure,delta_t,logsamp,max_no_neigh,nlist,ncoup,       &
      nlistsize,aham,thermal_field,beff,beff1,beff3,coord,ind_list_full,            &
      ind_nlistsize,ind_nlist,max_no_neigh_ind,sus_ind,do_mom_legacy,mode)
      !
      use prn_fields,       only : print_fields
      use prn_SpinIce,      only : print_vertices
      use Polarization,     only : print_pol
      use prn_averages,     only : print_averages
      use prn_topology,     only : print_topology
      use prn_currents,     only : print_currents
      use prn_microwaves,   only : print_mwf_fields
      use prn_trajectories, only : print_trajectories
      use prn_induced_info, only : print_ind_trajectories

      implicit none
      !
      integer, intent(in) :: NA            !< Number of atoms in one cell
      integer, intent(in) :: NT            !< Number of types of atoms
      integer, intent(in) :: N1            !< Number of cell repetitions in x direction
      integer, intent(in) :: N2            !< Number of cell repetitions in y direction
      integer, intent(in) :: N3            !< Number of cell repetitions in z direction
      integer, intent(in) :: mstep         !< Current simulation step
      integer, intent(in) :: Natom         !< Number of atoms in system
      integer, intent(in) :: nHam          !< Number of atoms in Hamiltonian
      integer, intent(in) :: Nchmax        !< Max number of chemical components on each site in cell
      integer, intent(in) :: Mensemble     !< Number of ensembles
      integer, intent(in) :: do_ralloy     !< Random alloy simulation (0/1)
      integer, intent(in) :: plotenergy    !< Calculate and plot energy (0/1)
      integer, intent(in) :: Natom_full    !< Number of atoms for full system (=Natom if not dilute)
      integer, intent(in) :: max_no_neigh  !< Maximum number of neighbours for the neighbour lists
      integer, intent(in) :: max_no_neigh_ind  !< Maximum number of neighbours for the induced moments
      real(dblprec), intent(in)     :: Temp                 !< Temperature
      real(dblprec), intent(in)     :: delta_t              !< Time step for real time measurement
      real(dblprec), intent(in)     :: temprescale             !< Temperature rescaling if QHB
      real(dblprec), intent(in)     :: temprescalegrad         !< Temperature rescaling gradient if QHB
      character(len=1), intent(in)  :: mode                 !< Simulation mode (S=SD, M=MC, H=MC Heat Bath, P=LD, C=SLD, G=GNEB)
      character(len=8), intent(in)  :: simid                !< Name of simulation
      character(len=1), intent(in)  :: logsamp              !< Sample measurements logarithmically (Y/N)
      character(len=1), intent(in)  :: do_mom_legacy         !< Flag to print/read moments in legacy output
      character(len=1), intent(in)  :: real_time_measure    !< Measurements displayed in real time
      integer, dimension(:), intent(in)         :: aham           !< Hamiltonian look-up table
      integer, dimension(:), intent(in)         :: atype          !< Type of atom
      integer, dimension(:), intent(in)         :: achem_ch       !< Chemical type of atoms (reduced list) (achem_ch(i)=achtype(acellnumbrev(i)))
      integer, dimension(:), intent(in)         :: nlistsize      !< Size of neighbour list for Heisenberg exchange couplings
      integer, dimension(:), intent(in)         :: asite_ch       !< Actual site of atom for dilute system
      integer, dimension(Natom),intent(in)      :: ind_nlistsize
      integer, dimension(Natom), intent(in)     :: ind_list_full
      integer, dimension(:,:), intent(in)       :: nlist          !< Neighbour list for Heisenberg exchange couplings
      integer, dimension(max_no_neigh_ind,Natom), intent(in) :: ind_nlist
      real(dblprec), dimension(Natom), intent(in) :: sus_ind
      real(dblprec), dimension(3,Natom)          :: coord !< Coordinates for all the atoms
      real(dblprec), dimension(:,:), intent(in)  :: mmom    !< Magnitude of magnetic moments
      real(dblprec), dimension(:,:,:), intent(in)   :: ncoup !< Heisenberg exchange couplings
      real(dblprec), dimension(:,:,:), intent(in)    :: emom   !< Current unit moment vector
      real(dblprec), dimension(:,:,:), intent(in)    :: emomM  !< Current magnetic moment vector
      real(dblprec), dimension(:,:,:), intent(inout) :: beff !< Current site dependent total effective field
      real(dblprec), dimension(:,:,:), intent(inout) :: beff1 !< Current site dependent internal field from spin Hamiltonian
      real(dblprec), dimension(:,:,:), intent(inout) :: beff3 !< Current site dependent internal field from mixed spin-lattice Hamiltonian
      real(dblprec), dimension(:,:,:), intent(inout) :: thermal_field !< Current site dependent stochastic field

      ! Local variables
      integer :: sstep

      sstep = f_logstep(mstep,logsamp)

      ! Print Microwave fields
      call print_mwf_fields(sstep,mstep,Natom,delta_t,real_time_measure,simid)
      ! Print trajectories of atoms
      call print_trajectories(Natom,sstep,mstep,Mensemble,emom,mmom,delta_t,        &
         real_time_measure,simid,do_mom_legacy,mode)
      ! Print vertex information for loop alogirthm
      call print_vertices(Natom,sstep,mstep,Mensemble,emom,simid)
      ! Print magnon and heat currents
      call print_currents(Natom,nHam, sstep,mstep,Mensemble,emomM,delta_t,          &
         real_time_measure,simid,ncoup,nlistsize,nlist,aham,max_no_neigh,coord)
      ! Print averages, cumulants and autocorrelation
      call print_averages(NA,NT,N1,N2,N3,Natom,mstep,sstep,Nchmax,Mensemble,        &
         do_ralloy,Natom_full,plotenergy,atype,achem_ch,asite_ch,mmom,emom,emomM,   &
         Temp,temprescale,temprescalegrad,delta_t,real_time_measure,simid,mode)
      ! Print the effective and thermal fields
      call print_fields(mstep,sstep,Natom,Mensemble,simid,real_time_measure,delta_t,&
         beff,thermal_field,beff1,beff3,emom)
      ! Print topological information of the sample
      call print_topology(NT,NA,N1,N2,sstep,mstep,Natom,Mensemble,delta_t,             &
         real_time_measure,emomM,emom,simid,atype)
      ! Print the polarization measurement
      call print_pol(sstep,mstep,Natom,Mensemble,max_no_neigh,nlist,nlistsize,emomM,&
         delta_t,simid,real_time_measure)

      ! Print information about the induced moments
      call print_ind_trajectories(Natom,Mensemble,sstep,mstep,max_no_neigh_ind,     &
         ind_nlistsize,ind_list_full,ind_nlist,delta_t,emomM,simid,                 &
         real_time_measure,sus_ind)

   end subroutine measure

   !-----------------------------------------------------------------------------
   !  SUBROUTINE: calc_mavrg
   !> @brief
   !> Calculate current average magnetization
   !> @todo Merge this routing with other buffer mavrg and prnavrg?
   !> @note this function is called by chelper from the C++ code
   !-----------------------------------------------------------------------------
   subroutine calc_mavrg(Natom, Mensemble, emomM, mavrg)
      !.. Implicit declarations
      implicit none

      integer, intent(in) :: Natom                                        !< Number of atoms in system
      integer, intent(in) :: Mensemble                                    !< Number of ensembles
      real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: emomM    !< Current unit  and magnitude moment vector
      real(dblprec), intent(out) :: mavrg                                 !< Current average magnetization

      !.. Scalar variables
      integer :: k

      !.. Local arrays
      real(dblprec), dimension(3,Mensemble) ::  m

      !.. Executable statements
      m(:,:) = 0.0_dblprec

      do k=1,Mensemble
            m(:,k) = sum(emomM(:,:,k),2)
      end do
      mavrg=sum(sum(m(:,:)**2,1)**0.5_dblprec)/Mensemble/Natom


   end subroutine calc_mavrg

   !---------------------------------------------------------------------------------
   !  SUBROUTINE: flush_measurements
   !> @brief
   !> Print remaining buffered measurements if simulation ends before all data has been written.
   !-----------------------------------------------------------------------------
   subroutine flush_measurements(Natom,Mensemble,NT,NA,N1,N2,N3,simid,mstep,emom,   &
      mmom,Nchmax,atype,real_time_measure,mcnstep,ind_list_full,do_mom_legacy,mode)

      use prn_fields,       only : flush_prn_fields
      use prn_SpinIce,      only : flush_vertices
      use Polarization,     only : flush_polarization
      use prn_averages,     only : flush_averages
      use prn_topology,     only : flush_topology
      use prn_microwaves,   only : flush_mwf_fields
      use prn_trajectories, only : flush_trajectories
      use prn_currents,     only : flush_currents
      use prn_induced_info, only : flush_ind_trajectories

      implicit none

      integer, intent(in) :: NT  !< Number of types of atoms
      integer, intent(in) :: NA  !< Number of atoms in one cell
      integer, intent(in) :: N1  !< Number of cell repetitions in x direction
      integer, intent(in) :: N2  !< Number of cell repetitions in y direction
      integer, intent(in) :: N3  !< Number of cell repetitions in z direction
      integer, intent(in) :: mstep !< Current simulation step
      integer, intent(in) :: Natom !< Number of atoms in system
      integer, intent(in) :: Nchmax    !< Max number of chemical components on each site in cell
      integer, intent(in) :: mcnstep        !< Total number of simulation step
      integer, intent(in) :: Mensemble !< Number of ensembles
      integer, dimension(Natom), intent(in) :: atype !< Type of atom
      integer, dimension(Natom), intent(in) :: ind_list_full
      real(dblprec), dimension(Natom, Mensemble) :: mmom   !< Magnitude of magnetic moments
      real(dblprec), dimension(3,Natom, Mensemble) :: emom !< Current unit moment vector
      character(len=1), intent(in) :: mode   !< Simulation mode (S=SD, M=MC, H=MC Heat Bath, P=LD, C=SLD, G=GNEB)
      character(len=8), intent(in) :: simid             !< Name of simulation
      character(len=1), intent(in) :: do_mom_legacy      !< Flag to print/read moments in legacy output
      character(len=1), intent(in) :: real_time_measure !< Measurements displayed in real time

      ! Flush the microwave field measurements
      call flush_mwf_fields(Natom,real_time_measure,simid)
      ! Flush the trajectories measurements
      call flush_trajectories(Natom,Mensemble,simid,real_time_measure,mmom,         &
         do_mom_legacy,mode)
      ! Flush the vertex measurements
      call flush_vertices(Natom,mcnstep,Mensemble,simid)
      ! Flush averages, cumulants and autocorrelation
      call flush_averages(NT,NA,N1,N2,N3,mstep,Natom,Nchmax,Mensemble,atype,mmom,emom,&
         simid,real_time_measure,mode)
      ! Flush thermal and effective fields
      call flush_prn_fields(Natom,Mensemble,simid,real_time_measure)
      ! Flush the polarization measurements
      call flush_polarization(mstep,Natom,Mensemble,simid,real_time_measure)
      ! Flush the topological information of the simulation
      call flush_topology(NT,Natom,simid,real_time_measure)
      ! Flush the current measurements
      call flush_currents(Natom,Mensemble,real_time_measure,simid)
      ! Flush the induced moments measurements
      call flush_ind_trajectories(Natom,Mensemble,ind_list_full,simid,real_time_measure)

   end subroutine flush_measurements

   !-----------------------------------------------------------------------------
   !  SUBROUTINE: allocate_measurementdata
   !> @brief
   !> Allocate arrays for sampling data
   !-----------------------------------------------------------------------------
   subroutine allocate_measurementdata(NA,NT,Natom,Mensemble,Nchmax,plotenergy,flag)

      use prn_fields,       only : allocate_field_print
      use prn_averages,     only : averages_allocations
      use Polarization,     only : allocate_prn_pol
      use prn_topology,     only : allocate_topology
      use prn_currents,     only : allocate_currents
      use prn_microwaves,   only : mwf_initializations
      use prn_trajectories, only : allocate_trajectories
      use prn_induced_info, only : allocate_ind_trajectories

      implicit none

      integer, intent(in) :: NA         !< Number of atoms in the unit cell
      integer, intent(in) :: NT         !< Number of types of atoms
      integer, intent(in) :: flag       !< Allocate or deallocate (1/-1)
      integer, intent(in) :: Natom      !< Number of atoms in system
      integer, intent(in) :: Nchmax     !< Max number of chemical components on each site in cell
      integer, intent(in) :: Mensemble  !< Number of ensembles
      integer, intent(in) :: plotenergy !< Plot the total energy

      ! Allocate field printing arrays
      call allocate_field_print(Natom,Mensemble,flag)
      ! Allocate and initialize the appropiated microwave field measurements
      call mwf_initializations(Natom,flag)
      ! Allocate and initialize the appropiated trajectory measurements arrays
      call allocate_trajectories(Natom,Mensemble,flag)
      ! Allocate and initialize the appropiated averages, cumulants and autocorrelation arrays
      call averages_allocations(NA,NT,Nchmax,Mensemble,plotenergy,flag)
      ! Allocate the measurement of polarization arrays
      call allocate_prn_pol(Mensemble,flag)
      ! Allocate the arrays for measuring topological quantities
      call allocate_topology(NT,Natom,Mensemble,flag)
      ! Allocate the arrays needed for measuring the magnon and heat currents
      call allocate_currents(Natom,Mensemble,flag)
      ! Allocate arrays for the measurement of induced moments
      call allocate_ind_trajectories(Natom,Mensemble,flag)

   end subroutine allocate_measurementdata

end module Measurements
