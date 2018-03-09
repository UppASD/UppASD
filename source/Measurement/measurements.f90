!-------------------------------------------------------------------------------
!> MODULE: Measurements
!> Data and routines for measuring averages and trajectories
!> @copyright
!! Copyright (C) 2008-2018 UppASD group
!! This file is distributed under the terms of the
!! GNU General Public License. 
!! See http://www.gnu.org/copyleft/gpl.txt
!-------------------------------------------------------------------------------
module Measurements
  use Parameters
  use Profiling
  !
  implicit none
  !
  private

  !public subroutines
  public :: measure, flush_measurements, do_measurements, allocate_measurementdata, calc_mavrg, logstep, logfun

contains

  !-----------------------------------------------------------------------------
  !> SUBROUTINE: do_measurements
  !> @Brief
  !> Helper method for determining if measurements are neccessary
  !> used during gpu acceleration
  !>
  !> @Author
  !> Niklas Fejes
  !-----------------------------------------------------------------------------
  subroutine do_measurements(mstep, do_avrg, do_tottraj, avrg_step, ntraj, tottraj_step, &
       traj_step, do_cumu, cumu_step, logsamp, do_copy)
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

    integer :: i
    integer :: sstep

    sstep = logstep(mstep,logsamp)

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



  !--------------------------------------------!
  ! @date 2014/09/01 - Thomas Nystrand
  ! - Moved to separate function
  !--------------------------------------------!
  function logstep(mstep,logsamp) result(sstep)
    character(len=1)   :: logsamp
    integer,intent(in) :: mstep

    integer            :: sstep

    if(logsamp=='Y') then
       sstep=logfun(mstep)
       if(sstep<0) return
    else
       sstep=mstep
    end if
  end function logstep

  !-----------------------------------------------------------------------------
  !> SUBROUTINE: measure
  !> @brief
  !> Wrapper routine for sampling and printing averages and trajectories
  !-----------------------------------------------------------------------------
  subroutine measure(Natom, Mensemble, NT, NA, N1, N2, N3, simid, mstep, emom, emomM, mmom, &
             Nchmax, do_ralloy, Natom_full, asite_ch, achem_ch, atype,  plotenergy, Temp, &
             real_time_measure, delta_t,logsamp, max_no_neigh, nlist,ncoup,nlistsize,&
             thermal_field,beff,beff1,beff3,coord,ind_mom,ind_nlistsize,ind_nlist,atype_ch)
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
    integer, intent(in) :: Nchmax        !< Max number of chemical components on each site in cell
    integer, intent(in) :: Mensemble     !< Number of ensembles
    integer, intent(in) :: do_ralloy     !< Random alloy simulation (0/1)
    integer, intent(in) :: plotenergy    !< Calculate and plot energy (0/1)
    integer, intent(in) :: Natom_full    !< Number of atoms for full system (=Natom if not dilute)
    integer, intent(in) :: max_no_neigh  !< Maximum number of neighbours for the neighbour lists
    integer, dimension(:), intent(in)      :: atype       !< Type of atom
    integer, dimension(:), intent(in)      :: atype_ch
    integer, dimension(:), intent(in)      :: achem_ch    !< Chemical type of atoms (reduced list) (achem_ch(i)=achtype(acellnumbrev(i)))
    integer, dimension(:), intent(in)      :: nlistsize   !< Size of neighbour list for Heisenberg exchange couplings
    integer, dimension(:), intent(in) :: asite_ch      !< Actual site of atom for dilute system
    integer, dimension(:,:), intent(in) :: nlist !< Neighbour list for Heisenberg exchange couplings
    integer, dimension(Natom),intent(in) :: ind_nlistsize
    integer, dimension(NA,Nchmax), intent(in) :: ind_mom
    integer, dimension(max_no_neigh,Natom), intent(in) :: ind_nlist
    character(len=8), intent(in)               :: simid                !< Name of simulation
    character(len=1), intent(in)               :: logsamp              !< Sample measurements logarithmically (Y/N)
    character(len=1), intent(in)               :: real_time_measure    !< Measurements displayed in real time
    real(dblprec), intent(in)                  :: Temp                    !< Temperature
    real(dblprec), intent(in)                  :: delta_t                 !< Time step for real time measurement
    real(dblprec), dimension(3,Natom)          :: coord !< Coordinates for all the atoms
    real(dblprec), dimension(:,:,:), intent(in)   :: ncoup !< Heisenberg exchange couplings
    real(dblprec), dimension(:,:), intent(in)     :: mmom    !< Magnitude of magnetic moments
    real(dblprec), dimension(:,:,:), intent(in)    :: emom   !< Current unit moment vector
    real(dblprec), dimension(:,:,:), intent(in)    :: emomM  !< Current magnetic moment vector
    real(dblprec), dimension(:,:,:), intent(inout) :: beff !< Current site dependent total effective field
    real(dblprec), dimension(:,:,:), intent(inout) :: beff1 !< Current site dependent internal field from spin Hamiltonian
    real(dblprec), dimension(:,:,:), intent(inout) :: beff3 !< Current site dependent internal field from mixed spin-lattice Hamiltonian
    real(dblprec), dimension(:,:,:), intent(inout) :: thermal_field !< Current site dependent stochastic field

    ! Local variables
    integer :: sstep

    sstep = logstep(mstep,logsamp)

    ! Print Microwave fields
    call print_mwf_fields(sstep,mstep,Natom,delta_t,real_time_measure,simid)
    ! Print trajectories of atoms
    call print_trajectories(Natom,sstep,mstep,Mensemble,emom,mmom,delta_t,real_time_measure,simid)
    ! Print vertex information for loop alogirthm
    call print_vertices(Natom,sstep,mstep,Mensemble,emom,simid)
    ! Print magnon and heat currents
    call print_currents(Natom,sstep,mstep,Mensemble,emomM,delta_t,real_time_measure,simid,&
         ncoup,nlistsize,nlist,max_no_neigh,coord)
    ! Print averages, cumulants and autocorrelation
    call print_averages(NA,NT,N1,N2,N3,Natom,mstep,sstep,Nchmax,Mensemble,do_ralloy,&
                        Natom_full,plotenergy,atype,achem_ch,asite_ch,mmom,emom,emomM,&
                        Temp,delta_t,real_time_measure,simid)
    ! Print the effective and thermal fields
    call print_fields(mstep, sstep, Natom, Mensemble, simid, real_time_measure, delta_t, &
         beff, thermal_field, beff1, beff3, emom)
    ! Print topological information of the sample
    call print_topology(NT,N1,N2,sstep,mstep,Natom,Mensemble,delta_t,real_time_measure,&
        emomM,emom,simid,atype)
    ! Print the polarization measurement
    call print_pol(sstep,mstep,Natom,Mensemble,max_no_neigh,nlist,nlistsize,&
         emomM,delta_t,simid,real_time_measure)

    ! Print information about the induced moments
    call print_ind_trajectories(NA,Natom,sstep,mstep,Nchmax,Mensemble,do_ralloy,Natom_full,max_no_neigh,&
         atype,ind_nlistsize,atype_ch,ind_mom,ind_nlist,delta_t,mmom,emom,simid,real_time_measure)

  end subroutine measure

  !-----------------------------------------------------------------------------
  !> SUBROUTINE: calc_mavrg
  !> @brief
  !> Calculate current average magnetization
  !> @TODO Merge this routing with other buffer mavrg and prnavrg?
  !> @Note this function is called by chelper from the C++ code
  !-----------------------------------------------------------------------------
  subroutine calc_mavrg(Natom, Mensemble, emomM, mavrg)
    !.. Implicit declarations
    implicit none

    integer, intent(in) :: Natom                                        !< Number of atoms in system
    integer, intent(in) :: Mensemble                                    !< Number of ensembles
    real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: emomM    !< Current unit  and magnitude moment vector
    real(dblprec), intent(out) :: mavrg                                 !< Current average magnetization

    !.. Scalar variables
    integer :: i,k

    !.. Local arrays
    real(dblprec), dimension(Mensemble) ::  mavrg_local
    real(dblprec), dimension(3,Mensemble) ::  m

    !.. Executable statements
    m(:,:)=0.0d0

#if _OPENMP >= 201307 && __INTEL_COMPILER < 1800
    !$omp parallel do default(shared),private(k,i),reduction(+:m),collapse(2),schedule(static)
#endif
    do k=1,Mensemble
       do i=1, Natom
          m(:,k) = m(:,k) + emomM(:,i,k)
       end do
    end do
#if _OPENMP >= 201307 && __INTEL_COMPILER < 1800
    !$omp end parallel do
#endif
    mavrg_local(:)=(m(1,:)**2+m(2,:)**2+m(3,:)**2)**0.5d0
    mavrg=sum(mavrg_local)/Mensemble/Natom

  end subroutine calc_mavrg

  !-----------------------------------------------------------------------------
  !> FUNCTION: logfun
  !> @brief
  !> Function to write the measurements in logarithmic scale
  !-----------------------------------------------------------------------------
  integer function logfun(mstep)
    implicit none
    !
    integer, intent(in) :: mstep !< Current simulation step
    !
    integer :: sstep,pstep,tstep
    !
    if(mstep==0) then
       sstep=0
    else
       tstep=int(log(real(mstep))*10)
       pstep=int(log(real(mstep-1))*10)
       if(tstep>pstep) then
          sstep=tstep
       else
          sstep=-1
       end if
    end if
    logfun=sstep
  end function logfun

  !-----------------------------------------------------------------------------
  !> SUBROUTINE: flush_measurements
  !> @brief
  !> Print remaining buffered measurements if simulation ends befor all data has been written.
  !-----------------------------------------------------------------------------
  subroutine flush_measurements(Natom, Mensemble, NT, NA, N1, N2, N3, simid, mstep,emom, mmom, &
             Nchmax,atype,real_time_measure,mcnstep,do_ralloy,Natom_full,atype_ch,ind_mom)

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
    integer, intent(in) :: do_ralloy
    integer, intent(in) :: Natom_full
    integer, dimension(Natom), intent(in) :: atype !< Type of atom
    integer, dimension(Natom_full), intent(in) :: atype_ch
    integer, dimension(NA,Nchmax), intent(in) :: ind_mom
    real(dblprec), dimension(Natom, Mensemble) :: mmom   !< Magnitude of magnetic moments
    real(dblprec), dimension(3,Natom, Mensemble) :: emom !< Current unit moment vector
    character(len=8), intent(in) :: simid             !< Name of simulation
    character(len=1), intent(in) :: real_time_measure !< Measurements displayed in real time

    ! Flush the microwave field measurements
    call flush_mwf_fields(Natom,real_time_measure,simid)
    ! Flush the trajectories measurements
    call flush_trajectories(Natom,Mensemble,simid,real_time_measure)
    ! Flush the vertex measurements
    call flush_vertices(Natom,mcnstep,Mensemble,simid)
    ! Flush averages, cumulants and autocorrelation
    call flush_averages(NT,NA,N1,N2,N3,mstep,Natom,Nchmax,Mensemble,atype,mmom,emom,&
         simid,real_time_measure)
    ! Flush thermal and effective fields
    call flush_prn_fields(Natom,Mensemble,simid,real_time_measure)
    ! Flush the polarization measurements
    call flush_polarization(mstep,Natom,Mensemble,simid,real_time_measure)
    ! Flush the topological information of the simulation
    call flush_topology(NT,Natom,simid,real_time_measure)
    ! Flush the current measurements
    call flush_currents(Natom,Mensemble,real_time_measure,simid)
    ! Flush the induced moments measurements
    call flush_ind_trajectories(NA,Natom,Nchmax,Mensemble,do_ralloy,Natom_full,atype,&
         atype_ch,ind_mom,simid,real_time_measure)

  end subroutine flush_measurements

  !-----------------------------------------------------------------------------
  !> SUBROUTINE: allocate_measurementdata
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
