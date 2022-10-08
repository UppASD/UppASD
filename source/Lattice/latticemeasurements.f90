!> Data and routines for measuring averages, trajectories, 
!> energies, and temperature for the displacement and velocites variables

module Latticemeasurements
   use Parameters
   use Profiling
   use Math_functions, only : f_logstep
   !
   implicit none
   !
   private

   !public subroutines
   public :: lattmeasure, flush_lattmeasurements, allocate_lattmeasurementdata, local_angular_momentum


contains


   !> Wrapper routine for sampling and printing averages and trajectories
   subroutine lattmeasure(simid, Natom, Mensemble, NT, NA, N1, N2, N3, &
         mstep, logsamp, coord, mion, uvec, vvec, lvec, eeff, eeff1, eeff3, ethermal_field, &
         Nchmax, do_ralloy, Natom_full, asite_ch, achem_ch, atype, Temp, &
         ldpot_energy, sdpot_energy, sldpot_energy, totpot_energy)
      !
      use prn_latticefields,   only : print_lattfields
      use prn_latticeaverages,     only : print_lattaverages
      use prn_latticetrajectories, only : print_latttrajectories

      implicit none
      !
      character(len=8), intent(in) :: simid    !< Name of simulation
      integer, intent(in) :: Natom         !< Number of atoms in system
      integer, intent(in) :: Mensemble     !< Number of ensembles
      integer, intent(in) :: NT            !< Number of types of atoms
      integer, intent(in) :: NA            !< Number of atoms in one cell
      integer, intent(in) :: N1            !< Number of cell repetitions in x direction
      integer, intent(in) :: N2            !< Number of cell repetitions in y direction
      integer, intent(in) :: N3            !< Number of cell repetitions in z direction

      integer, intent(in) :: mstep         !< Current simulation step
      character(len=1), intent(in) :: logsamp  !< Sample measurements logarithmically (Y/N)
      real(dblprec), dimension(3,Natom) :: coord !< Coordinates for all the atoms
      real(dblprec), dimension(Natom,Mensemble), intent(in) :: mion   !< Ionic mass
      real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: uvec   !< Current ionic displacement
      real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: vvec   !< Current ionic velocity
      real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: lvec   !< Current local angular momentum
      real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: eeff          !< Current site dependent total effective e-field
      real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: eeff1   !< Current site dependent internal effective e-field
      real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: eeff3   !< Current site dependent internal effective e-field from mixed spin-lattice Hamiltonian
      real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: ethermal_field !< Current site dependent stochastic e-field

      integer, intent(in) :: Nchmax        !< Max number of chemical components on each site in cell
      integer, intent(in) :: do_ralloy     !< Random alloy simulation (0/1)
      integer, intent(in) :: Natom_full    !< Number of atoms for full system (=Natom if not dilute)
      integer, dimension(:), intent(in) :: asite_ch      !< Actual site of atom for dilute system
      integer, dimension(:), intent(in) :: achem_ch    !< Chemical type of atoms (reduced list) (achem_ch(i)=achtype(acellnumbrev(i)))
      integer, dimension(:), intent(in) :: atype       !< Type of atom
      real(dblprec), intent(in) :: Temp   !< Temperature
      real(dblprec), dimension(Mensemble), intent(in) :: ldpot_energy     !< LD potential energy
      real(dblprec), dimension(Mensemble), intent(in) :: sdpot_energy        !< SD potential energy
      real(dblprec), dimension(Mensemble), intent(in) :: sldpot_energy    !< SLD potential energy (without pure LD or SD potential energies)
      real(dblprec), dimension(Mensemble), intent(in) :: totpot_energy    !< Total potential energy: LD + SD + SLD. No kinetic energy!

      ! Local variables
      integer :: sstep

      sstep = f_logstep(mstep,logsamp)

      ! SLDTODO: sample ionic DOF energies and temperature

      ! Print trajectories of atoms
      call print_latttrajectories(simid, Natom, Mensemble, mstep, sstep, mion, uvec, vvec, coord)

      ! Print averages, cumulants and autocorrelation
      call print_lattaverages(simid, Natom, Mensemble, NT, NA, N1, N2, N3, &
         mstep, sstep, logsamp, coord, mion, uvec, vvec, lvec, eeff, Temp, &
         Nchmax, do_ralloy, Natom_full, atype, achem_ch, asite_ch, &
         ldpot_energy, sdpot_energy, sldpot_energy, totpot_energy)

      ! Print the effective and thermal fields
      call print_lattfields(simid, Natom, Mensemble, mstep, sstep, eeff, eeff1, eeff3, ethermal_field)

   end subroutine lattmeasure


   !> Print remaining buffered measurements if simulation ends before all data has been written.
   subroutine flush_lattmeasurements(simid, Natom, Mensemble, NT, NA, N1, N2, N3, &
         mstep, mion, uvec, vvec, coord, Nchmax, atype, Temp)

      use prn_latticefields,       only : flush_lattfields
      use prn_latticeaverages,     only : flush_lattaverages
      use prn_latticetrajectories, only : flush_latttrajectories

      implicit none

      character(len=8), intent(in) :: simid             !< Name of simulation
      integer, intent(in) :: Natom !< Number of atoms in system
      integer, intent(in) :: Mensemble !< Number of ensembles
      integer, intent(in) :: NT  !< Number of types of atoms
      integer, intent(in) :: NA  !< Number of atoms in one cell
      integer, intent(in) :: N1  !< Number of cell repetitions in x direction
      integer, intent(in) :: N2  !< Number of cell repetitions in y direction
      integer, intent(in) :: N3  !< Number of cell repetitions in z direction

      integer, intent(in) :: mstep !< Current simulation step
      real(dblprec), dimension(Natom,Mensemble), intent(in) :: mion   !< Ionic mass
      real(dblprec), dimension(3, Natom, Mensemble), intent(in) :: uvec   !< Current ionic displacement
      real(dblprec), dimension(3, Natom, Mensemble), intent(in) :: vvec   !< Current ionic velocity
      real(dblprec), dimension(3, Natom), intent(in) :: coord !< Coordinates of atoms
      integer, intent(in) :: Nchmax    !< Max number of chemical components on each site in cell
      integer, dimension(Natom), intent(in) :: atype !< Type of atom
      real(dblprec), intent(in) :: Temp                      !< Temperature of the system (scalar)

      ! Flush the trajectories measurements
      call flush_latttrajectories(simid, Natom, Mensemble, coord)

      ! Flush averages, cumulants and autocorrelation
      call flush_lattaverages(simid, Natom, Mensemble, NT, NA, N1, N2, N3, &
         mstep, Nchmax, atype, mion, uvec, vvec, Temp)

      ! Flush thermal and effective fields
      call flush_lattfields(simid, Natom, Mensemble)

   end subroutine flush_lattmeasurements


   !> Allocate arrays for sampling data
   subroutine allocate_lattmeasurementdata(Natom, Mensemble, NA, NT, Nchmax, flag)

      use prn_latticefields,       only : allocate_lattfield_print
      use prn_latticeaverages,     only : allocate_lattaverages
      use prn_latticetrajectories, only : allocate_latttrajectories

      implicit none

      integer, intent(in) :: Natom      !< Number of atoms in system
      integer, intent(in) :: Mensemble  !< Number of ensembles
      integer, intent(in) :: NA         !< Number of atoms in the unit cell
      integer, intent(in) :: NT         !< Number of types of atoms
      integer, intent(in) :: Nchmax     !< Max number of chemical components on each site in cell
      integer, intent(in) :: flag       !< Allocate or deallocate (1/-1)


      ! Allocate field printing arrays
      call allocate_lattfield_print(Natom, Mensemble, flag)

      ! Allocate and initialize the appropiated trajectory measurements arrays
      call allocate_latttrajectories(Natom, Mensemble,flag)

      ! Allocate and initialize the appropiated averages, cumulants and autocorrelation arrays
      call allocate_lattaverages(Natom, Mensemble, NA, NT, Nchmax, flag)

   end subroutine allocate_lattmeasurementdata

   !> Calculate local angular momentum for lattice
   subroutine local_angular_momentum(Natom, Mensemble, mion,  uvec, vvec, lvec)

      use Math_functions, only : f_cross_product

      integer, intent(in) :: Natom !< Number of atoms in system
      integer, intent(in) :: Mensemble   !< Number of ensembles
      real(dblprec), dimension(Natom,Mensemble), intent(in) :: mion   !< Ionic mass
      real(dblprec), dimension(3, Natom, Mensemble), intent(in) :: uvec   !< Current ionic displacement
      real(dblprec), dimension(3, Natom, Mensemble), intent(in) :: vvec   !< Current ionic velocity
      real(dblprec), dimension(3, Natom, Mensemble), intent(out) :: lvec   !< Current local angular momentum

      integer :: i,k

      !$omp parallel do default(shared) private(k,i), collapse(2)
      do k=1, Mensemble
         do i=1,Natom
            lvec(:,i,k) = mion(i,1) * f_cross_product(uvec(:,i,k),vvec(:,i,k))
         end do
      end do

   end subroutine local_angular_momentum


end module Latticemeasurements
