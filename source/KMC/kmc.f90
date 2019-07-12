!-------------------------------------------------------------------------------
! MODULE: KMC
!
! DESCRIPTION:
!> This module contains the main routines to perform a standard Kinetic Monte Calo
!> (KMC) i.e. rare dynamics approach.
!> @details In this approach particles can go from one state to
!> another by overcoming an energy barrier, with a mean time for this to occur
!> dictated by an Arrhenius law.
!> This module contains all the different approaches to consider the dynamics
!> of KMC active particles in the sample.
!> In this approach a set of "particles" are considered to be able to move, overcoming
!> a certain energy barrier defined by the user.
!> As these particles move they locally change the magnitude of the magnetic moments
!> and exchage interactions present in the system.
!
!> @author
!> Jonathan Chico and Anders Bergman
!> @copyright
!> GNU Public License.
!> @note
!> Due to the inherent assumed symmetry breaking in this routine, the reduced
!> Hamiltonian approach (`do_reduced`) is not supported in this routine
!-------------------------------------------------------------------------------
module KMC
   use Parameters
   use Profiling
   use Constants
   use KMCData

   implicit none

   integer, dimension(:,:), allocatable :: nmdim_barriers !< Dimension of neighbour map for barriers
   integer, dimension(:,:,:), allocatable :: nm_barriers !< Neighbour map for barriers

   public

contains

   !-----------------------------------------------------------------------------
   ! SUBROUTINE: initialize_KMC
   !
   ! DESCRIPTION:
   !> Finds the location of the KMC particles in the starting configuration and calculates
   !> the neighbour map necessary for the energy barriers
   !
   !> @author
   !> Jonathan Chico
   !-----------------------------------------------------------------------------
   subroutine initialize_KMC(NT,NA,N1,N2,N3,sym,Natom,Nchmax,NA_KMC,do_ralloy,Mensemble,&
         Natom_full,max_no_neigh,do_prnstruct,max_no_shells_barriers,nn_barriers,atype,anumb,&
         nlistsize,atype_ch,acellnumb,nlist,C1,C2,C3,Bas,coord,ncoup,redcoord_barriers,&
         kmc_barriers,BC1,BC2,BC3,do_sortcoup,simid,max_no_neigh_barriers,asite_ch,&
         achem_ch,kmc_index_prv,kmc_index_aft,nlistsize_barriers,nlist_barriers,&
         ncoup_barriers,mmom,emom,emomM)

      use NeighbourMap, only : setup_nm

      implicit none

      ! .. Input variables
      integer, intent(in) :: NT  !< Number of types of atoms
      integer, intent(in) :: NA  !< Number of atoms in one cell
      integer, intent(in) :: N1  !< Number of cell repetitions in x direction
      integer, intent(in) :: N2  !< Number of cell repetitions in y direction
      integer, intent(in) :: N3  !< Number of cell repetitions in z direction
      integer, intent(in) :: sym !< Symmetry of system (0-3)
      integer, intent(in) :: Natom !< Number of atoms in system
      integer, intent(in) :: Nchmax !< Max number of chemical components on each site in cell
      integer, intent(in) :: NA_KMC !< Number of KMC particles found in the system
      integer, intent(in) :: do_ralloy !< Random alloy simulation (0/1)
      integer, intent(in) :: Mensemble !< Number of ensembles
      integer, intent(in) :: Natom_full !< Number of atoms for full system (=Natom if not dilute)
      integer, intent(in) :: max_no_neigh !< Calculated maximum of neighbours for exchange
      integer, intent(in) :: do_prnstruct !< Print Hamiltonian information (0-4)
      integer, intent(in) :: max_no_shells_barriers !< Calculated maximum of shells for KMC barriers
      integer, dimension(:), intent(in) :: nn_barriers !< Number of neighbour shells for barriers
      integer, dimension(Natom), intent(in) :: anumb !< Atom number in cell
      integer, dimension(Natom), intent(in) :: nlistsize !< Size of neighbour list for Heisenberg exchange couplings
      integer, dimension(Natom_full), intent(in) :: atype_ch  !< Actual type of atom for dilute system
      integer, dimension(Natom_full), intent(in) :: acellnumb !< List for translating atom no. in full cell to actual cell
      integer, dimension(max_no_neigh,Natom), intent(in) :: nlist !< Neighbour list for Heisenberg exchange couplings
      real(dblprec), dimension(3), intent(in) :: C1 !< First lattice vector
      real(dblprec), dimension(3), intent(in) :: C2 !< Second lattice vector
      real(dblprec), dimension(3), intent(in) :: C3 !< Third lattice vector
      real(dblprec), dimension(3,NA), intent(in) :: Bas !< Coordinates for basis atoms
      real(dblprec), dimension(3,Natom), intent(in) :: coord !< Coordinates of atoms
      real(dblprec), dimension(max_no_neigh,Natom), intent(inout) :: ncoup !< Heisenberg exchange couplings
      real(dblprec), dimension(NT,max_no_shells_barriers,3), intent(in) :: redcoord_barriers  !< Coordinates for Heisenberg exchange couplings
      real(dblprec), dimension(NT,max_no_shells_barriers,Nchmax,max(Nchmax,NT),1), intent(in) :: kmc_barriers !< Energy barriers for the KMC (input)
      character(len=1), intent(in) :: BC1 !< Boundary conditions in x-direction
      character(len=1), intent(in) :: BC2 !< Boundary conditions in y-direction
      character(len=1), intent(in) :: BC3 !< Boundary conditions in z-direction
      character(len=1), intent(in) :: do_sortcoup !< Sort the entries of ncoup arrays (Y/N)
      character(len=8), intent(in) :: simid !< Name of simulation
      ! .. In/out variables
      integer, intent(inout) :: max_no_neigh_barriers   !< Calculated maximum of neighbours for KMC barriers
      integer, dimension(Natom), intent(inout) :: atype !< Type of atom
      integer, dimension(Natom_full), intent(inout) :: asite_ch  !< Actual site of atom for dilute system
      integer, dimension(Natom_full), intent(inout) :: achem_ch  !< Chemical type of atoms (reduced list)
      integer, dimension(:), allocatable, intent(inout) :: kmc_index_aft !< Indices indicating positions of the KMC paritcles after an update
      integer, dimension(:), allocatable, intent(inout) :: kmc_index_prv !< Indices indicating positions of the KMC paritcles before an update
      integer, dimension(:), allocatable, intent(inout) :: nlistsize_barriers !< Size of neighbour list for KMC barriers
      integer, dimension(:,:), allocatable, intent(inout) :: nlist_barriers !< Neighbour list for Heisenberg exchange couplings
      real(dblprec), dimension(:,:), allocatable, intent(inout) :: ncoup_barriers !< KMC barriers
      real(dblprec), dimension(Natom,Mensemble), intent(inout) :: mmom    !< Magnitude of magnetic moments
      real(dblprec), dimension(3,Natom,Mensemble), intent(inout) :: emom  !< Current unit moment vector
      real(dblprec), dimension(3,Natom,Mensemble), intent(inout) :: emomM !< Current magnetic moment vector

      ! .. Local variables
      integer :: max_no_equiv_barriers

      ! Create the nieghbour list for the energy barriers
      ! This should be done by taking the neighbour map as in the case of the
      ! exchange
      write(*,'(2x,a)',advance='no') "Setting up KMC barriers "
      call setup_nm(Natom, NT, NA, N1, N2, N3, C1, C2, C3, BC1, BC2, BC3,1,&
         atype, Bas, max_no_neigh_barriers, max_no_shells_barriers, &
         max_no_equiv_barriers, sym, nn_barriers, redcoord_barriers, &
         nm_barriers, nmdim_barriers, do_ralloy, Natom_full, acellnumb, atype_ch)

      ! Now that the neighbour map is set, the barriers can now be mapped
      call allocate_kmc(NA_KMC,Natom,max_no_neigh_barriers,1)

      ! Map the energy barriers to ensure that they can be read in from a file
      ! similar to the jfile
      call setup_barriers(NT,Natom,Nchmax,do_ralloy,Natom_full,max_no_neigh_barriers,&
         max_no_equiv_barriers,max_no_shells_barriers,nn_barriers,anumb,atype,atype_ch,&
         asite_ch,achem_ch,nmdim_barriers,nm_barriers,do_sortcoup,kmc_barriers,&
         nlistsize_barriers,nlist_barriers,ncoup_barriers)
      write(*,'(a)') "done"

      if (do_prnstruct.eq.1.or.do_prnstruct.eq.4) then
         write(*,'(2x,a)',advance='no') "Printing KMC barriers "
         ! Printing the kmc barriers
         call printing_barriers(Natom,max_no_neigh_barriers,nlistsize_barriers,nlist_barriers,&
            ncoup_barriers,simid)
         write(*,'(a)') "done"
      endif

      write(*,'(2x,a)',advance='no') "Setting initial KMC particles positions "
      ! Finding the positions of the KMC particles in the simulation cell
      call kmc_part_init(Natom,NA_KMC,do_ralloy,Natom_full,coord,Bas_KMC,atype,&
         kmc_index_prv,achem_ch_kmc,atype_inp_KMC,achem_ch)
      write(*,'(a)') "done"

   end subroutine initialize_KMC

   !------------------------------------------------------------------------------
   ! SUBROUTINE: initial_wrapper_kmc
   !
   ! DESCRIPTION:
   !> Wraps up all the possible kmc routines that might be needed to run for the initial step
   !> kmc_mode 1 Standard fixed barrier kmc
   !
   !> @author
   !> Jonathan Chico
   !------------------------------------------------------------------------------
   subroutine initial_wrapper_kmc(mstep,Natom,NA_KMC,do_ralloy,Mensemble,kmc_method,&
         rate0,delta_t,Temp_array,coord,atype,kmc_index_prv,&
         kmc_index_aft,kmc_time_steps,time_efield)

      implicit none

      ! .. Input variables
      integer, intent(in) :: mstep !< Current simulation step
      integer, intent(in) :: Natom !< Number of atoms in system
      integer, intent(in) :: NA_KMC !< Number of KMC particles found in the system
      integer, intent(in) :: do_ralloy !< Random alloy simulation (0/1)
      integer, intent(in) :: Mensemble !< Number of ensembles
      integer, intent(in) :: kmc_method !< Which is the type of KMC algorithm tight now only kmc_method=1 (usual KMC) is used
      integer, intent(in) :: time_efield
      real(dblprec), intent(in) :: rate0 !< Attempt rate frequency, this is a mechanism dependent variable
      real(dblprec), intent(in) :: delta_t !< Time step
      real(dblprec), dimension(Natom), intent(in) :: Temp_array !< Temperature array
      real(dblprec), dimension(3,Natom), intent(in) :: coord !< Coordinates of atoms

      ! .. In/out variables
      integer, dimension(Natom), intent(inout) :: atype !< Type of atom
      integer, dimension(NA_KMC), intent(inout) :: kmc_index_aft  !< Indices indicating positions of the KMC paritcles after an update
      integer, dimension(NA_KMC), intent(inout) :: kmc_index_prv  !< Indices indicating positions of the KMC paritcles before and update
      integer, dimension(NA_KMC), intent(inout) :: kmc_time_steps !< Number of steps until next jump for all the KMC particles

      ! .. Local variables
      integer :: kmc_part

      ! It is important to notice that only the standard KMC is a "proper" KMC method, since in principle
      ! the number of states needed to be able to apply KMC is finite which is clearly not the case when one
      ! introduces thermal fluctuations, the spins can go anywhere they like, the question is whether the electron
      ! can or cannot do that and how does the motion of the spins affect this.

      ! Standard KMC
      if (kmc_method.eq.1) then
         ! This loop ensures that the rates and new positions are calculated for all
         ! KMC particles
         ! Probably need to do the trick done in normal MC in which an auxiliary list is created
         ! with the positions in the list being random as to not bias the choice of
         ! events
         ! The KMC index is updated individually, probably it is best to have it as an in/out
         ! variable with the checks of auto-update being done in-place
         ! The same is done for the kmc_time_steps which are sent as in/out variables
         do kmc_part=1, NA_KMC
            call trans_rate_calc(NA_KMC,Natom,max_no_neigh_barriers,rate0,delta_t,&
               coord,Temp_array,kmc_time_steps,kmc_index_prv,kmc_part,kmc_index_aft,&
               KMC_changed,mstep,time_efield)
         enddo
      endif

   end subroutine initial_wrapper_kmc

   !-----------------------------------------------------------------------------
   ! SUBROUTINE: updater_wrapper_kmc
   !
   ! DESCRIPTION:
   !> Wraps up all the possible kmc routines that might be needed to update the
   !> exchange and if necessary recalculate the KMC time
   !> kmc_mode 1 Standard fixed barrier kmc
   !
   !> @author
   !> Jonathan Chico
   !-----------------------------------------------------------------------------
   subroutine updater_wrapper_kmc(mstep,Natom,NA_KMC,do_ralloy,Mensemble,Natom_full,&
         kmc_method,max_no_neigh,nlistsize,nlist,rate0,delta_t,&
         Temp_array,coord,atype,kmc_index_prv,kmc_index_aft,kmc_time_steps,&
         achem_ch,mmom,ncoup,emom,emomM,time_efield)

      use MonteCarlo, only : choose_random_atom_x

      implicit none

      ! .. Input variables
      integer, intent(in) :: mstep !< Current simulation step
      integer, intent(in) :: Natom !< Number of atoms in system
      integer, intent(in) :: NA_KMC !< Number of KMC particles found in the system
      integer, intent(in) :: do_ralloy !< Random alloy simulation (0/1)
      integer, intent(in) :: Mensemble !< Number of ensembles
      integer, intent(in) :: Natom_full !< Number of atoms for full system (=Natom if not dilute)
      integer, intent(in) :: kmc_method !< Which is the type of KMC algorithm tight now only kmc_method=1 (usual KMC) is used
      integer, intent(in) :: time_efield
      integer, intent(in) :: max_no_neigh !< Calculated maximum of neighbours for exchange
      integer, dimension(Natom), intent(in) :: nlistsize !< Size of neighbour list for Heisenberg exchange couplings
      integer, dimension(max_no_neigh,Natom), intent(in) :: nlist !< Neighbour list for Heisenberg exchange couplings

      real(dblprec), intent(in) :: rate0 !< Attempt rate frequency, this is a mechanism dependent variable
      real(dblprec), intent(in) :: delta_t !< Time step
      real(dblprec), dimension(Natom), intent(in) :: Temp_array !< Temperature array
      real(dblprec), dimension(3,Natom), intent(in) :: coord !< Coordinates of atoms

      ! .. In/out variables
      integer, dimension(Natom), intent(inout) :: atype !< Type of atom
      integer, dimension(NA_KMC), intent(inout) :: kmc_index_aft  !< Indices indicating positions of the KMC paritcles after an update
      integer, dimension(NA_KMC), intent(inout) :: kmc_index_prv  !< Indices indicating positions of the KMC paritcles before and update
      integer, dimension(NA_KMC), intent(inout) :: kmc_time_steps !< Number of steps until next jump for all the KMC particles
      integer, dimension(Natom_full), intent(inout) :: achem_ch   !< Chemical type of atoms (reduced list)

      real(dblprec), dimension(Natom,Mensemble), intent(inout) :: mmom     !< Magnitude of magnetic moments
      real(dblprec), dimension(max_no_neigh,Natom), intent(inout) :: ncoup !< Heisenberg exchange couplings
      real(dblprec), dimension(3,Natom,Mensemble), intent(inout) :: emom   !< Current unit moment vector
      real(dblprec), dimension(3,Natom,Mensemble), intent(inout) :: emomM  !< Current magnetic moment vector

      ! .. Local variables
      integer :: kmc_part,kmc_part_pos,kmc_part_check, kmc_part_check_pos, ii
      logical*1 :: kmc_overlap

      ! Standard KMC
      if (kmc_method==1) then
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         ! The first thing that one needs to do is to find out which the KMC particle
         ! that must be updated
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         do kmc_part=1,NA_KMC
            if ( kmc_time_steps(kmc_part).eq.mstep ) then
               kmc_part_pos=kmc_index_aft(kmc_part)
               call update_hamiltonian(Natom,NA_KMC,kmc_part,do_ralloy,Mensemble,Natom_full,&
                  max_no_neigh,nlistsize,kmc_index_aft,kmc_index_prv,nlist,coord,emom,&
                  atype,achem_ch,mmom,ncoup,emomM)

               ! After the update of the Hamiltonian one must update the kmc particles
               ! indices so that they represent the new crystal and chemical structure

               kmc_index_prv(kmc_part)=kmc_index_aft(kmc_part)

               ! Probably moot but placed here for didactic reasons
               kmc_part_pos=kmc_index_prv(kmc_part)

               ! Re-calculate the rates and determine the new polaron candiate position
               call trans_rate_calc(NA_KMC,Natom,max_no_neigh_barriers,rate0,delta_t,&
                  coord,Temp_array,kmc_time_steps,kmc_index_prv,kmc_part,kmc_index_aft,&
                  KMC_changed,mstep,time_efield)

               !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
               ! Need to ensure that if a particle that was going to move to a site that is
               ! now occupied that rate is re-calculated
               !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
               do kmc_part_check=1,NA_KMC
                  if ((kmc_part_check.ne.kmc_part)) then
                     kmc_overlap=.false.
                     kmc_part_check_pos=kmc_index_prv(kmc_part_check)
                     do ii=1,nlistsize_barriers(kmc_part_pos)
                        kmc_overlap=kmc_overlap.or.(nlist_barriers(ii,kmc_part_pos)==kmc_index_prv(kmc_part_check))&
                           .or.(nlist_barriers(ii,kmc_part_pos)==kmc_index_aft(kmc_part_check))
                     end do
                     if(kmc_overlap) then
                        call trans_rate_calc(NA_KMC,Natom,max_no_neigh_barriers,rate0,delta_t,&
                           coord,Temp_array,kmc_time_steps,kmc_index_prv,kmc_part_check,kmc_index_aft,&
                           KMC_changed,mstep,time_efield)
                     end if
                  end if
               enddo
               !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
               ! End of check
               !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
               if (kmc_index_prv(kmc_part).eq.kmc_index_aft(kmc_part)) then
                  KMC_changed=.false.
               endif

            endif ! Selectively update
         end do ! Loop over KMC particles
      endif ! KMC method

   end subroutine updater_wrapper_kmc

   !------------------------------------------------------------------------------
   ! SUBROUTINE: kmc_part_init
   !
   ! DESCRIPTION:
   !> Subroutine to find the initial location of the kmc particles with respect
   !> to the positions of the host system
   !
   !> @author
   !> Jonathan Chico
   !------------------------------------------------------------------------------
   subroutine kmc_part_init(Natom,NA_KMC,do_ralloy,Natom_full,coord,Bas_KMC,atype,&
         kmc_index_prv,achem_ch_kmc,atype_inp_KMC,achem_ch)

      implicit none

      !.. Input variables
      integer, intent(in) :: Natom !< Number of atoms in system
      integer, intent(in) :: NA_KMC !< Number of KMC particles found in the system
      integer, intent(in) :: do_ralloy !< Random alloy simulation (0/1)
      integer, intent(in) :: Natom_full !< Number of atoms for full system (=Natom if not dilute)
      real(dblprec), dimension(3,Natom), intent(in) :: coord !< Coordinates of atoms
      real(dblprec), dimension(3,NA_KMC), intent(in) :: Bas_KMC !< Coordinates for basis atoms of the KMC particles

      !.. In/out variables
      integer, dimension(Natom), intent(inout) :: atype !< Type of atom
      integer, dimension(NA_KMC), intent(inout) :: kmc_index_prv !< Mapper of indices of the kmc particle positions
      integer, dimension(NA_KMC), intent(inout) :: achem_ch_kmc  !< Chemical species for the KMC particles
      integer, dimension(NA_KMC), intent(inout) :: atype_inp_KMC !< Atomic types for the KMC particles
      integer, dimension(Natom_full), intent(inout) :: achem_ch !< Chemical type of atoms (reduced list)

      !.. Local variables
      integer :: ii,iatom
      real(dblprec) :: tol,dist_coord

      tol=0.005_dblprec

      ! Do a loop over all the atoms in the system
      do iatom=1,Natom
         ! Do a loop over the KMC particles
         do ii=1, NA_KMC
            dist_coord=sqrt((coord(1,iatom)-Bas_KMC(1,ii))**2+&
               (coord(2,iatom)-Bas_KMC(2,ii))**2+&
               (coord(3,iatom)-Bas_KMC(3,ii))**2)
            if (do_ralloy.eq.0) then
               ! This will store the location of the KMC particles in the system
               if ((dist_coord.le.tol).and.(atype(iatom).eq.atype_inp_KMC(ii))) kmc_index_prv(ii)=iatom
            else
               ! For random alloys not only the position if not the chemical type must also match
               ! Must be careful with the posibility of having zeros stored here
               if ((dist_coord.le.tol).and.(achem_ch(iatom).eq.achem_ch_kmc(ii))) then
                  kmc_index_prv(ii)=iatom
               endif
            endif
         enddo
      enddo

   end subroutine kmc_part_init

   !------------------------------------------------------------------------------
   ! SUBROUTINE: Trans_rate_calc
   !
   ! DESCRIPTION:
   !> @brief This routine will calculate the tranisition rate between the possible states
   !> of this way the \f$J_{ij}\f$'s will be updated at the time calculates by this routine
   !> thus allowing the dynamics of magnetic polarons. This is done by a regular
   !> KMC method
   !
   !> This is known as kmc_method 1
   !
   !> @author
   !> Jonathan Chico
   !------------------------------------------------------------------------------
   subroutine trans_rate_calc(NA_KMC,Natom,max_no_neigh_barriers,rate0,delta_t,&
         coord,Temp_array,kmc_time_steps,kmc_index_prv,kmc_part,kmc_index_aft,&
         KMC_changed,mstep,time_efield)

      use RandomNumbers, only : rng_uniform

      implicit none

      integer, intent(in) :: mstep
      integer, intent(in) :: Natom !< Number of atoms in system
      integer, intent(in) :: NA_KMC !< Number of KMC particles found in the system
      integer, intent(in) :: kmc_part ! Index of the considered kmc particle
      integer, intent(in) :: time_efield
      integer, intent(in) :: max_no_neigh_barriers !< Calculated maximum of neighbours for KMC barriers
      integer, dimension(NA_KMC), intent(in) :: kmc_index_prv !< KMC particle indices before the selection of a next candidate

      real(dblprec), intent(in) :: rate0 !< Attempt rate frequency, this is a mechanism dependent variable
      real(dblprec), intent(in) :: delta_t !< Time step
      real(dblprec), dimension(Natom), intent(in) :: Temp_array !< Temperature array
      real(dblprec), dimension(3,Natom), intent(in) :: coord !< Coordinates of atoms

      ! .. In/out variables
      logical, intent(inout) :: KMC_changed
      integer, dimension(NA_KMC), intent(inout) :: kmc_index_aft  !< Indices indicating positions of the KMC paritcles before and update
      integer, dimension(NA_KMC), intent(inout) :: kmc_time_steps !< Number of steps until next jump for all the KMC particles

      ! .. Local variables
      integer :: ii, itemp,curr_neigh,kmc_part_pos,count,iprev
      integer :: nfree_neigh
      integer, dimension(max_no_neigh_barriers) :: visited_neigh
      integer, dimension(max_no_neigh_barriers) :: free_neigh

      real(dblprec) :: urand, trand,temp
      real(dblprec) :: dnorm,KMC_time
      real(dblprec) :: Efield_ene
      real(dblprec), dimension(2) :: rrand
      real(dblprec), dimension(3) :: dr
      real(dblprec), dimension(max_no_neigh_barriers) :: DeltaE
      logical :: accepted, full_occupied, trapped

      dr=0.0_dblprec
      Rate=0.0_dblprec
      dnorm=0.0_dblprec
      DeltaE=0.0_dblprec
      visited_neigh=0
      free_neigh=0
      count=0
      Efield_ene=0.0_dblprec
      !.. Logical
      accepted=.false.
      full_occupied=.false.
      trapped=.false.

      ! The energy barrier has to be calculated between all the possible states,
      ! in this case the states correspond to the position of the seed of the polaron
      ! that is one needs to calculate the probability of moving from the initial position
      ! to any other position

      ! The energy barriers are read from an special file which has the same
      ! structure than the jfile

      ! Calculate the total energy barrier

      kmc_part_pos=kmc_index_prv(kmc_part)

      do ii=1, nlistsize_barriers(kmc_part_pos)

         dr(1)= coord(1,nlist_barriers(ii,kmc_part_pos))-coord(1,kmc_part_pos)
         dr(2)= coord(2,nlist_barriers(ii,kmc_part_pos))-coord(2,kmc_part_pos)
         dr(3)= coord(3,nlist_barriers(ii,kmc_part_pos))-coord(3,kmc_part_pos)

         ! Norm of the vector
         dnorm=sqrt(dr(1)**2+dr(2)**2+dr(3)**2)

         ! The efield is in eV/m

         ! Check if the electric field is active
         if (do_efield.eq.'Y') then
            Efield_ene=efield(1)*dr(1)+efield(2)*dr(2)+efield(3)*dr(3)
         else if (do_efield.eq.'T'.and.(mstep.le.time_efield)) then
            Efield_ene=efield(1)*dr(1)+efield(2)*dr(2)+efield(3)*dr(3)
         else if (do_efield.eq.'T'.and.(mstep.le.time_efield)) then
            Efield_ene=0.0_dblprec
         endif

         ! Calculating the barrier for nieghbour ii when the Efield is considered
         DeltaE(ii)=ncoup_barriers(ii,kmc_part_pos)+Efield_ene

      enddo

      ! Do a loop over the four possible paths that the electron can follow
      ! Rate stores the cumulative sum, Rate_full each individual
      Rate(1)=rate0*exp(-DeltaE(1)/(k_bolt_ev*Temp_array(kmc_part_pos)))
      Rate_full(1)=rate0*exp(-DeltaE(1)/(k_bolt_ev*Temp_array(kmc_part_pos)))

      do ii=2, nlistsize_barriers(kmc_part_pos)
         ! Creating of the partial summation array
         Rate(ii)=Rate(ii-1)+rate0*exp(-DeltaE(ii)/(k_bolt_ev*Temp_array(kmc_part_pos)))
         Rate_full(ii)=rate0*exp(-DeltaE(ii)/(k_bolt_ev*Temp_array(kmc_part_pos)))
      enddo

      ! Storing all barriers
      Rate_free=0.0_dblprec

      ! See which positions are occupied/blocked
      ! Storing the cumulative rate sum in Rate_free
      visited_neigh=0
      free_neigh=0
      nfree_neigh=0
      do ii=1,nlistsize_barriers(kmc_part_pos)
         ! Store that this neighbour has been visited
         curr_neigh=nlist_barriers(ii,kmc_part_pos)
         ! Now one must check whether the neighbour is in the list of kmc particles
         if (any(kmc_index_prv.eq.curr_neigh)) then
            visited_neigh(ii)=1
         else
            nfree_neigh=nfree_neigh+1
            free_neigh(nfree_neigh)=ii
            if(nfree_neigh==1) then
               Rate_free(nfree_neigh)=Rate_full(ii)
            else
               Rate_free(nfree_neigh)=Rate_free(nfree_neigh-1)+Rate_full(ii)
            end if
         endif
      end do

      ! Check that the particle is not trapped
      if(nfree_neigh==0) then
         trapped=.true.
         accepted=.false.
         KMC_changed=.false.
         kmc_index_aft(kmc_part)=kmc_index_prv(kmc_part)
      else
         ! Choose among the available barriers to jump
         call rng_uniform(rrand,2)
         urand=rrand(1)*Rate_free(nfree_neigh)
         iprev=0
         do while(urand.gt.Rate_free(iprev).and.iprev<=nfree_neigh)
            iprev=iprev+1
         end do
         if(iprev<=nfree_neigh) then
            ii=free_neigh(iprev)
            ! Selected state from the partial summations
            kmc_index_aft(kmc_part)=nlist_barriers(ii,kmc_part_pos)
            accepted=.true.
            KMC_changed=.true.
            accepted=.false.
         else
            print *,'no dice'
         end if
      end if

      if (KMC_changed) then
         ! Calculate the actual time from the KMC
         trand=rrand(2)
         ! A point for debate. Should the time scale with all rates
         ! or only for those possible (i.e. not blocked paths)
         KMC_time=(1.0_dblprec/Rate_free(nfree_neigh))*log(1.0_dblprec/trand)
         if (KMC_time.ge.delta_t) then
            temp=KMC_time/delta_t
            itemp=nint(temp)+mstep
         else
            itemp=1+mstep
         endif
      else
         itemp=1+mstep
      endif

      ! Transform it to actual time step
      ! Typecast the time as integer ammount of time steps
      kmc_time_steps(kmc_part)=itemp

   end subroutine trans_rate_calc

   !------------------------------------------------------------------------------
   ! SUBROUTINE: update_hamiltonian
   !
   ! DESCRIPTION:
   !> This routine is a wrapper that will update the crystallographic, chemical and
   !> interactions of the system after the kmc particles move
   !
   !> @author
   !> Jonathan Chico
   !------------------------------------------------------------------------------
   subroutine update_hamiltonian(Natom,NA_KMC,kmc_part,do_ralloy,Mensemble,Natom_full,&
         max_no_neigh,nlistsize,kmc_index_aft,kmc_index_prv,nlist,coord,emom,&
         atype,achem_ch,mmom,ncoup,emomM)

      implicit none

      ! .. Input variables
      integer, intent(in) :: Natom !< Number of atoms in system
      integer, intent(in) :: NA_KMC !< Number of KMC particles found in the system
      integer, intent(in) :: kmc_part
      integer, intent(in) :: do_ralloy !< Random alloy simulation (0/1)
      integer, intent(in) :: Mensemble !< Number of ensembles
      integer, intent(in) :: Natom_full !< Number of atoms for full system (=Natom if not dilute)
      integer, intent(in) :: max_no_neigh !< Calculated maximum of neighbours for exchange

      integer, dimension(Natom), intent(in) :: nlistsize !< Size of neighbour list for Heisenberg exchange couplings
      integer, dimension(max_no_neigh,Natom), intent(in) :: nlist !< Neighbour list for Heisenberg exchange couplings

      real(dblprec), dimension(3,Natom), intent(in) :: coord !< Coordinates of atoms
      real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: emom !< Current unit moment vector

      ! .. In/out variables
      integer, dimension(Natom), intent(inout) :: atype !< Type of atom
      integer, dimension(NA_KMC), intent(inout) :: kmc_index_aft  !< Indices indicating positions of the KMC paritcles after an update
      integer, dimension(NA_KMC), intent(inout) :: kmc_index_prv  !< Indices indicating positions of the KMC paritcles before and update
      integer, dimension(Natom_full), intent(inout) :: achem_ch  !< Chemical type of atoms (reduced list)

      real(dblprec), dimension(Natom,Mensemble), intent(inout) :: mmom     !< Magnitude of magnetic moments
      real(dblprec), dimension(max_no_neigh,Natom), intent(inout) :: ncoup !< Heisenberg exchange couplings
      real(dblprec), dimension(3,Natom,Mensemble), intent(inout) :: emomM  !< Current magnetic moment vector

      ! Update the actual crystallographic, chemical information and moments of the
      ! previous and current positions of the KMC particles
      call update_cent_atom(Natom,NA_KMC,kmc_part,Mensemble,Natom_full,atype,&
         kmc_index_prv,kmc_index_aft,achem_ch,mmom,emom,emomM)

      ! Update the chemical, crystallographic and moment inromation of the neighbours
      ! and ther neighbours of the previous and current positions.
      ! The exchange interactions are also updated so that they represent KMC particles
      ! moving and affecting the interactions
      call update_neighbours(Natom,NA_KMC,kmc_part,do_ralloy,Mensemble,Natom_full,&
         max_no_neigh,nlistsize,kmc_index_aft,kmc_index_prv,nlist,coord,emom,&
         atype,achem_ch,mmom,emomM,ncoup)


   end subroutine update_hamiltonian

   !------------------------------------------------------------------------------
   ! SUBROUTINE: update_cent_atom
   !
   ! DESCRIPTION:
   !> This routine updates the information of the atoms at the position of the KMC
   !> particles. Both the positions before and after the motion are updated.
   !> @details The update is done by swapping the chemical data, crystallographic  data and
   !> moment magnitude of the sites. With the new setup information then a new
   !> Hamiltonian can be set
   !
   !> @author
   !> Jonathan Chico
   !------------------------------------------------------------------------------
   subroutine update_cent_atom(Natom,NA_KMC,kmc_part,Mensemble,Natom_full,atype,&
         kmc_index_prv,kmc_index_aft,achem_ch,mmom,emom,emomM)

      implicit none

      integer, intent(in) :: Natom !< Number of atoms in system
      integer, intent(in) :: NA_KMC !< Number of KMC particles found in the system
      integer, intent(in) :: kmc_part
      integer, intent(in) :: Mensemble !< Number of ensembles
      integer, intent(in) :: Natom_full !< Number of atoms for full system (=Natom if not dilute)
      real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: emom   !< Current unit moment vector

      integer, dimension(Natom), intent(inout) :: atype !< Type of atom
      integer, dimension(NA_KMC), intent(inout) :: kmc_index_prv !< Indices indicating positions of the KMC paritcles after an update
      integer, dimension(NA_KMC), intent(inout) :: kmc_index_aft !< Indices indicating positions of the KMC paritcles before and update
      integer, dimension(Natom_full), intent(inout) :: achem_ch  !< Chemical type of atoms (reduced list)

      real(dblprec), dimension(Natom,Mensemble), intent(inout) :: mmom     !< Magnitude of magnetic moments
      real(dblprec), dimension(3,Natom,Mensemble), intent(inout) :: emomM  !< Current magnetic moment vector

      ! .. Local variables
      integer :: temp_type, temp_chem
      real(dblprec) :: temp_mom

      ! Storing the data temporarily of the site which was not a kmc particle before
      ! the transfer rate is calculated
      temp_type=atype(kmc_index_aft(kmc_part))
      temp_chem=achem_ch(kmc_index_aft(kmc_part))
      temp_mom=mmom(kmc_index_aft(kmc_part),1)

      ! Updating the information of the new kmc particle position, with the one of the previous position
      atype(kmc_index_aft(kmc_part))=atype(kmc_index_prv(kmc_part))
      achem_ch(kmc_index_aft(kmc_part))=achem_ch(kmc_index_prv(kmc_part))
      mmom(kmc_index_aft(kmc_part),:)=mmom(kmc_index_prv(kmc_part),:)
      emomM(1,kmc_index_aft(kmc_part),:)=emom(1,kmc_index_aft(kmc_part),:)*mmom(kmc_index_aft(kmc_part),:)
      emomM(2,kmc_index_aft(kmc_part),:)=emom(2,kmc_index_aft(kmc_part),:)*mmom(kmc_index_aft(kmc_part),:)
      emomM(3,kmc_index_aft(kmc_part),:)=emom(3,kmc_index_aft(kmc_part),:)*mmom(kmc_index_aft(kmc_part),:)

      ! Now updating the properties of the previous position to be the ones of where the particle is located now
      atype(kmc_index_prv(kmc_part))=temp_type
      achem_ch(kmc_index_prv(kmc_part))=temp_chem
      mmom(kmc_index_prv(kmc_part),:)=temp_mom
      emomM(1,kmc_index_prv(kmc_part),:)=emom(1,kmc_index_prv(kmc_part),:)*mmom(kmc_index_prv(kmc_part),:)
      emomM(2,kmc_index_prv(kmc_part),:)=emom(2,kmc_index_prv(kmc_part),:)*mmom(kmc_index_prv(kmc_part),:)
      emomM(3,kmc_index_prv(kmc_part),:)=emom(3,kmc_index_prv(kmc_part),:)*mmom(kmc_index_prv(kmc_part),:)

      ! One still needs to do a loop over the neighbours and make sure that the
      ! chemical types and so on are properly updated

   end subroutine update_cent_atom

   !------------------------------------------------------------------------------
   ! SUBROUTINE: update_neighbours
   !
   ! DESCRIPTION:
   !> This routine updates the information of the neighbouring atoms to the KMC
   !> particles. Both the positions before and after the motion are updated.
   !> @details The update is done by swapping the chemical data, crystallographic  data and
   !> moment magnitude of the sites. With the new setup information then a new
   !> Hamiltonian can be set. It is important to update the neighbours to ensure
   !> that any new Hamiltonian has the correct local properties
   !
   !> @author
   !> Jonathan Chico
   !------------------------------------------------------------------------------
   subroutine update_neighbours(Natom,NA_KMC,kmc_part,do_ralloy,Mensemble,Natom_full,&
         max_no_neigh,nlistsize,kmc_index_aft,kmc_index_prv,nlist,coord,emom,&
         atype,achem_ch,mmom,emomM,ncoup)

      implicit none

      ! .. Input variables
      integer, intent(in) :: Natom !< Number of atoms in system
      integer, intent(in) :: NA_KMC !< Number of KMC particles found in the system
      integer, intent(in) :: kmc_part
      integer, intent(in) :: do_ralloy !< Random alloy simulation (0/1)
      integer, intent(in) :: Mensemble !< Number of ensembles
      integer, intent(in) :: Natom_full !< Number of atoms for full system (=Natom if not dilute)
      integer, intent(in) :: max_no_neigh !< Calculated maximum of neighbours for exchange

      integer, dimension(Natom), intent(in) :: nlistsize !< Size of neighbour list for Heisenberg exchange couplings
      integer, dimension(NA_KMC), intent(in) :: kmc_index_aft !< Indices indicating positions of the KMC paritcles after an update
      integer, dimension(NA_KMC), intent(in) :: kmc_index_prv !< Indices indicating positions of the KMC paritcles before and update
      integer, dimension(max_no_neigh,Natom), intent(in) :: nlist !< Neighbour list for Heisenberg exchange couplings

      real(dblprec), dimension(3,Natom), intent(in) :: coord !< Coordinates of atoms
      real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: emom   !< Current unit moment vector

      ! .. In/out variables
      integer, dimension(Natom), intent(inout) :: atype !< Type of atom
      integer, dimension(Natom_full), intent(inout) :: achem_ch  !< Chemical type of atoms (reduced list)
      real(dblprec), dimension(Natom,Mensemble), intent(inout) :: mmom     !< Magnitude of magnetic moments
      real(dblprec), dimension(max_no_neigh,Natom), intent(inout) :: ncoup !< Heisenberg exchange couplings
      real(dblprec), dimension(3,Natom,Mensemble), intent(inout) :: emomM  !< Current magnetic moment vector

      ! .. Local variables
      integer :: prev_type, prev_chem,curr_chem,curr_type
      integer :: prev_pos_neigh, curr_pos_neigh,i_neigh_prev
      integer :: icount, jcount, lcount, kcount,j_neigh_curr
      integer :: curr_pos, prev_pos, iatom, jatom
      real(dblprec) :: prev_mom, curr_mom,tol
      real(dblprec) :: temp_coup_prev, temp_coup_curr
      real(dblprec), dimension(3) :: rij_prev, rij_curr, rij_prev_neigh, rij_curr_neigh

      ! Value of the tolerance when comparing bonding vectors
      tol=0.005_dblprec

      ! This is the position of the KMC particle after the update
      ! The important thing is that this is used as the "HOST" to be copied over
      curr_pos=kmc_index_aft(kmc_part)
      ! Next position of the KMC particle before the update
      ! This contains the information of the "IMPURITY", that is the ones affected by the KMC particle
      prev_pos=kmc_index_prv(kmc_part)

      ! Check wheter one is dealing with a random alloy
      if (do_ralloy.eq.0) then

         ! Need to loop through all the neighbours of the previous position first
         ! That is this contains the info about the KMC particle
         do iatom=1,nlistsize(prev_pos)
            prev_pos_neigh=nlist(iatom,prev_pos)
            rij_prev(1)=coord(1,prev_pos_neigh)-coord(1,prev_pos)
            rij_prev(2)=coord(2,prev_pos_neigh)-coord(2,prev_pos)
            rij_prev(3)=coord(3,prev_pos_neigh)-coord(3,prev_pos)
            prev_type=atype(prev_pos)
            prev_mom=mmom(prev_pos,1)
            temp_coup_prev=ncoup(iatom,prev_pos)

            do jatom=1, nlistsize(curr_pos)
               curr_pos_neigh=nlist(iatom,prev_pos)
               rij_curr(1)=coord(1,curr_pos_neigh)-coord(1,curr_pos)
               rij_curr(2)=coord(2,curr_pos_neigh)-coord(2,curr_pos)
               rij_curr(3)=coord(3,curr_pos_neigh)-coord(3,curr_pos)
               curr_type=atype(curr_pos)
               curr_mom=mmom(curr_pos,1)
               temp_coup_curr=ncoup(jatom,curr_pos)

               ! If the neighbours have the same bonding vector
               if (abs(rij_prev(1)-rij_curr(1)).le.tol .and. abs(rij_prev(2)-rij_curr(2)).le. tol &
                  .and. abs(rij_prev(3)-rij_curr(3)).le. tol) then
                  atype(prev_pos)=curr_type
                  atype(curr_pos)=prev_type
                  mmom(prev_pos,:)=curr_mom
                  emomM(1,prev_pos,:)=mmom(prev_pos,:)*emom(1,prev_pos,:)
                  emomM(2,prev_pos,:)=mmom(prev_pos,:)*emom(2,prev_pos,:)
                  emomM(3,prev_pos,:)=mmom(prev_pos,:)*emom(3,prev_pos,:)
                  mmom(curr_pos,:)=prev_mom
                  emomM(1,curr_pos,:)=mmom(curr_pos,:)*emom(1,curr_pos,:)
                  emomM(2,curr_pos,:)=mmom(curr_pos,:)*emom(2,curr_pos,:)
                  emomM(3,curr_pos,:)=mmom(curr_pos,:)*emom(3,curr_pos,:)
                  ! The couplings between the atom with the previous KMC particle position
                  ! and its neighbours and the new KMC particle position and its
                  ! neighbours are swaped
                  ncoup(jatom,curr_pos)=temp_coup_prev
                  ncoup(iatom,prev_pos)=temp_coup_curr
               endif
            enddo
         enddo
         ! If the system is a random alloy
      else
         ! Need to loop through all the neighbours of the previous position first
         ! That is this contains the info about the KMC particle
         do iatom=1,nlistsize(prev_pos)
            prev_pos_neigh=nlist(iatom,prev_pos)
            rij_prev(1)=coord(1,prev_pos_neigh)-coord(1,prev_pos)
            rij_prev(2)=coord(2,prev_pos_neigh)-coord(2,prev_pos)
            rij_prev(3)=coord(3,prev_pos_neigh)-coord(3,prev_pos)
            prev_type=atype(prev_pos)
            prev_chem=achem_ch(prev_pos)
            prev_mom=mmom(prev_pos,1)
            temp_coup_prev=ncoup(iatom,prev_pos)

            ! This is the next position of the KMC particle
            do jatom=1, nlistsize(curr_pos)
               curr_pos_neigh=nlist(jatom,curr_pos)
               rij_curr(1)=coord(1,curr_pos_neigh)-coord(1,curr_pos)
               rij_curr(2)=coord(2,curr_pos_neigh)-coord(2,curr_pos)
               rij_curr(3)=coord(3,curr_pos_neigh)-coord(3,curr_pos)
               curr_type=atype(curr_pos)
               curr_chem=achem_ch(curr_pos)
               curr_mom=mmom(curr_pos,1)
               temp_coup_curr=ncoup(jatom,curr_pos)

               ! If the neighbours have the same bonding vector
               ! (this is not done with the modulus to let the algorithm be more general)
               if (abs(rij_prev(1)-rij_curr(1)).le. tol .and. &
                  abs(rij_prev(2)-rij_curr(2)).le. tol .and. &
                  abs(rij_prev(3)-rij_curr(3)).le. tol) then
                  ! Swapping the crystallographic types
                  atype(prev_pos)=curr_type
                  atype(curr_pos)=prev_type
                  ! Swapping the chemical types of the position
                  achem_ch(prev_pos)=curr_chem
                  achem_ch(curr_pos)=prev_chem
                  ! Swapping the magnetic moments between the previous and new positions
                  mmom(prev_pos,:)=curr_mom
                  emomM(1,prev_pos,:)=mmom(prev_pos,:)*emom(1,prev_pos,:)
                  emomM(2,prev_pos,:)=mmom(prev_pos,:)*emom(2,prev_pos,:)
                  emomM(3,prev_pos,:)=mmom(prev_pos,:)*emom(3,prev_pos,:)
                  mmom(curr_pos,:)=prev_mom
                  emomM(1,curr_pos,:)=mmom(curr_pos,:)*emom(1,curr_pos,:)
                  emomM(2,curr_pos,:)=mmom(curr_pos,:)*emom(2,curr_pos,:)
                  emomM(3,curr_pos,:)=mmom(curr_pos,:)*emom(3,curr_pos,:)
                  ! The couplings between the atom with the previous KMC particle position
                  ! and its neighbours and the new KMC particle position and its
                  ! neighbours are swaped
                  ncoup(jatom,curr_pos)=temp_coup_prev
                  ncoup(iatom,prev_pos)=temp_coup_curr
               endif
            enddo
         enddo
      endif

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Now the couplings are updated
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Loop over the neighbours of the previous position
      do icount=1, nlistsize(prev_pos)
         prev_pos_neigh=nlist(icount,prev_pos)
         ! Loop over the neighbours of the neighbours of the previous position
         do jcount=1, nlistsize(prev_pos_neigh)

            i_neigh_prev=nlist(jcount,prev_pos_neigh)
            ! Distance between the neighbours of the neighbours atoms of the kmc particle
            rij_prev_neigh(1)=coord(1,i_neigh_prev)-coord(1,prev_pos_neigh)
            rij_prev_neigh(2)=coord(2,i_neigh_prev)-coord(2,prev_pos_neigh)
            rij_prev_neigh(3)=coord(3,i_neigh_prev)-coord(3,prev_pos_neigh)
            ! Storing the coupling temporarily for the previous position neighbours
            temp_coup_prev=ncoup(jcount,prev_pos_neigh)

            ! Loop over the neighbours of the current position
            do kcount=1, nlistsize(curr_pos)
               curr_pos_neigh=nlist(kcount,curr_pos)
               ! Loop over the neighbours of the neighbours of the current position
               do lcount=1, nlistsize(curr_pos_neigh)

                  j_neigh_curr=nlist(lcount,curr_pos_neigh)
                  ! Distance between the neighbours of the neighbour atoms of the next kmc particle
                  rij_curr_neigh(1)=coord(1,j_neigh_curr)-coord(1,curr_pos_neigh)
                  rij_curr_neigh(2)=coord(2,j_neigh_curr)-coord(2,curr_pos_neigh)
                  rij_curr_neigh(3)=coord(3,j_neigh_curr)-coord(3,curr_pos_neigh)
                  ! Storing the coupling temporarily for the current position neighbours
                  temp_coup_curr=ncoup(lcount,curr_pos_neigh)

                  ! Checking if the bonding vectors are the same
                  if (abs(rij_prev_neigh(1)-rij_curr_neigh(1)).le.tol .and.&
                     abs(rij_prev_neigh(2)-rij_curr_neigh(2)).le.tol .and.&
                     abs(rij_prev_neigh(3)-rij_curr_neigh(3)).le.tol) then
                     ! Swapping the coupling constants
                     ncoup(lcount,curr_pos_neigh)=temp_coup_prev
                     ncoup(jcount,prev_pos_neigh)=temp_coup_curr

                  endif ! Check the bonding vector
               enddo ! Neighbour of the neighbours of the current kmc particle
            enddo ! Neighbour of the current kmc particle
         enddo ! Neighbour of the neighbour of the previous kmc particle
      enddo ! Neighbour of the previous kmc particle

   end subroutine update_neighbours

   !-------------------------------------------------------------------------------
   ! SUBROUTINE: setup_barriers
   !
   ! DESCRIPTION:
   !> In this the magnitude for the energy barriers for the neighbour lists for the
   !> KMC when dealing with the motion of an electron through the sample.
   !
   !> @author
   !  Jonathan Chico
   !-------------------------------------------------------------------------------
   subroutine setup_barriers(NT,Natom,Nchmax,do_ralloy,Natom_full,max_no_neigh_barriers,&
         max_no_equiv_barriers,max_no_shells_barriers,nn_barriers,anumb,atype,atype_ch,&
         asite_ch,achem_ch,nmdim_barriers,nm_barriers,do_sortcoup,kmc_barriers,&
         nlistsize_barriers,nlist_barriers,ncoup_barriers)
      !
      use Sorting, only : MergeSortIR
      !
      implicit none
      !
      integer, intent(in) :: NT    !< Number of types of atoms
      integer, intent(in) :: Natom !< Number of atoms in system
      integer, intent(in) :: Nchmax     !< Max number of chemical components on each site in cell
      integer, intent(in) :: do_ralloy  !< Random alloy simulation (0/1)
      integer, intent(in) :: Natom_full !< Number of atoms for full system (=Natom if not dilute)
      integer, intent(in) :: max_no_neigh_barriers  !< Calculated maximum of neighbours for KMC barriers
      integer, intent(in) :: max_no_equiv_barriers  !< Calculated maximum of neighbours in one shell for KMC barriers
      integer, intent(in) :: max_no_shells_barriers !< Calculated maximum of shells for KMC barriers
      integer, dimension(NT), intent(in):: nn_barriers !< Number of neighbour shells
      integer, dimension(Natom), intent(in) :: anumb   !< Atom number in cell
      integer, dimension(Natom), intent(in) :: atype   !< Type of atom
      integer, dimension(Natom_full), intent(in) :: atype_ch !< Actual site of atom for dilute system
      integer, dimension(Natom_full), intent(in) :: asite_ch !< Actual site of atom for dilute system
      integer, dimension(Natom_full), intent(in) :: achem_ch !< Chemical type of atoms (reduced list)
      integer, dimension(max_no_shells_barriers,Natom), intent(in) :: nmdim_barriers !< Dimension of neighbour map for barriers
      integer, dimension(Natom,max_no_shells_barriers,max_no_equiv_barriers), intent(in) :: nm_barriers !< Neighbour map

      character :: do_sortcoup !< Sort the entries of ncoup arrays (Y/N)
      real(dblprec), dimension(NT,max_no_shells_barriers,Nchmax,max(Nchmax,NT),1), intent(in) :: kmc_barriers !< KMC barriers input

      integer, dimension(Natom), intent(out):: nlistsize_barriers !< Size of neighbour list for KMC barriers
      integer, dimension(max_no_neigh_barriers,Natom), intent(out) :: nlist_barriers !< Neighbour list for KMC barriers
      real(dblprec), dimension(max_no_neigh_barriers,Natom), intent(out) :: ncoup_barriers !< KMC barriers
      !
      integer :: tempn
      integer :: i, j, k, l
      integer :: jatom, jchtype
      integer :: ncount,ncounter
      real(dblprec) :: tempcoup
      logical :: exis

      ncoup_barriers = 0.0_dblprec
      ! Factors for mRy energy conversion

      ! Loop over atoms
      !$omp parallel do default(shared) private(i,k,j,ncount,ncounter,exis,l,jatom,jchtype)
      do i=1, Natom
         ncount = 1
         ncounter = 1
         ! Shells
         do k=1, Nn_barriers(atype(i))
            ! Sites in shell
            do j=1, nmdim_barriers(k,i)
               ! Existing coupling
               if (nm_barriers(i,k,j)>0) then
                  exis=.false.
                  do l=1,ncount-1
                     if(nlist_barriers(l,i)==nm_barriers(i,k,j)) exis=.true.
                  end do
                  if(.not.exis) then
                     nlist_barriers(ncount,i) = nm_barriers(i,k,j)
                     if (do_ralloy==0) then
                        ncoup_barriers(ncount,i) = kmc_barriers(atype(i),k,1,atype(nm_barriers(i,k,j)),1)
                     else
                        jatom = nm_barriers(i,k,j)
                        jchtype = atype_ch(jatom)
                        ncoup_barriers(ncount,i) = kmc_barriers(atype_ch(i),k,achem_ch(i),achem_ch(jatom),1)
                     end if
                     ncount = ncount + 1
                  end if !If it does not exist
               end if ! If the coupling exists
            end do ! Sites in shell
         end do ! Shells
         nlistsize_barriers(i) = ncount-1
      end do ! Number of atoms
      !$omp end parallel do

      if(do_sortcoup == 'Y') then
         !$omp parallel do default(shared) private(i,j,k,tempn,tempcoup)
         do i=1,Natom
            do j=1,nlistsize_barriers(i)
               do k=1,nlistsize_barriers(i)-j
                  if ( nlist_barriers(k,i) .gt. nlist_barriers(k+1,i) ) then
                     tempn        = nlist_barriers(k,i)
                     nlist_barriers(k,i)   = nlist_barriers(k+1,i)
                     nlist_barriers(k+1,i) = tempn
                     tempcoup     = ncoup_barriers(k,i)
                     ncoup_barriers(k,i)   = ncoup_barriers(k+1,i)
                     ncoup_barriers(k+1,i) = tempcoup
                  end if ! Check if the position is the same
               end do ! j-th Site
            end do ! i-th Site
         end do ! Number of atoms
         !$omp end parallel do

      end if

   end subroutine setup_barriers

   !-------------------------------------------------------------------------------
   ! SUBROUTINE: printing_barriers
   !
   ! DESCRIPTION:
   !> In this subroutine the energy barriers will be printed, this might be useful for
   !> more complex systems that the ones treated previously
   !
   !> @author
   !  Jonathan Chico
   !-------------------------------------------------------------------------------
   subroutine printing_barriers(Natom,max_no_neigh_barriers,nlistsize_barriers,nlist_barriers,&
         ncoup_barriers,simid)
      !
      !.. Implicit declarations
      implicit none

      integer, intent(in) :: Natom !< Number of atoms in system
      integer, intent(in) :: max_no_neigh_barriers !< Calculated maximum of neighbours for the KMC barriers
      integer, dimension(Natom), intent(in) :: nlistsize_barriers !< Size of neighbour list for KMC barriers
      integer, dimension(max_no_neigh_barriers, Natom), intent(in) :: nlist_barriers !< Neighbour list for KMC barriers
      real(dblprec), dimension(max_no_neigh_barriers, Natom), intent(in) :: ncoup_barriers !< KMC barriers
      character(len=8),intent(in) :: simid !< Name of simulation

      !.. Local variables
      integer :: i
      character(len=30) :: filn

      !.. Executable statements
      write (filn,'(''struct_kmc.'',a8,''.out'')') simid
      open(ofileno, file=filn)

      ! print neighbor list - after sort
      write (ofileno,*) "Sorted data from KMC"
      do i=1,Natom
         write (ofileno,*) "----------------------------------"
         write (ofileno,70001) i,nlistsize_barriers(i)
         write (ofileno,70002) nlist_barriers(1:nlistsize_barriers(i),i)
         write (ofileno,70003) ncoup_barriers(1:nlistsize_barriers(i),i)
      end do
      close(ofileno)

      70001 format ("Atom=",i8,4x,"No neigh=",i7)
      70002 format ("            ",1X,5I6)
      70003 format (5es16.8)

   end subroutine printing_barriers

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!!!! DEPRECATED ROUTINES
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!------------------------------------------------------------------------------
   !! SUBROUTINE: dynamical_barrier
   !!
   !! DESCRIPTION:
   !!> This routine will calculate the magnetic barrier from moving the electron to
   !!> another site, the barrier is calculated as follows: The magnetic energy for the
   !!> polaron site is calculated with their current FM exchange, a "trial movement"
   !!> is made, the energy of the current configuration on the other site is calculated
   !!> by considering the current configuration (probably AFM) with the FM values.
   !!
   !!> This is referred as kmc_method 2
   !!
   !!> @author
   !!  Jonathan Chico
   !!------------------------------------------------------------------------------
   !  subroutine dynamical_barrier(Natom,Mensemble,max_no_neigh,prev_polaron,max_no_neigh_barriers,&
   !             atype,nlistsize,nlistsize_barriers,nlist,nlist_barriers,delta_t,lat_barr,Temp_array,&
   !             coord,mmom,ncoup,emom,emomM,next_polaron,time_jump)
   !
   !    use Constants
   !    use RandomNumbers, only : rng_uniform
   !
   !    implicit none
   !
   !    integer, intent(in) :: Natom
   !    integer, intent(in) :: Mensemble
   !    integer, intent(in) :: max_no_neigh
   !    integer, intent(in) :: prev_polaron
   !    integer, intent(in) :: max_no_neigh_barriers
   !    integer, dimension(Natom), intent(in) :: atype
   !    integer, dimension(Natom), intent(in) :: nlistsize
   !    integer, dimension(Natom), intent(in) :: nlistsize_barriers
   !    integer, dimension(max_no_neigh,Natom), intent(in) :: nlist
   !    integer, dimension(max_no_neigh_barriers,Natom), intent(in) :: nlist_barriers
   !
   !    real(dblprec), intent(in) :: lat_barr
   !    real(dblprec), intent(in) :: delta_t
   !    real(dblprec), dimension(Natom), intent(in) :: Temp_array
   !    real(dblprec), dimension(3,Natom), intent(in) :: coord
   !    real(dblprec), dimension(Natom,Mensemble), intent(in) :: mmom
   !    real(dblprec), dimension(max_no_neigh,Natom), intent(in) :: ncoup
   !    real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: emom
   !    real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: emomM
   !
   !    integer, intent(out) :: next_polaron
   !    integer, intent(out) :: time_jump
   !
   !    integer :: ii, kk, iatom,ene_barr,state, count_num,list_indx
   !    integer :: current_atom, current_neigh
   !
   !    real(dblprec) :: mdot, ene_curr,urand,trand,dnorm,Efield_ene,fcinv,mdot2,dummy_time
   !    real(dblprec), dimension(2) :: rrand
   !    real(dblprec), dimension(3) :: dr
   !    real(dblprec), dimension(max_no_neigh_barriers) :: cumu_barr,aux_list
   !    real(dblprec) :: fc2,KMC_time
   !
   !    fc2=2*mry/mub
   !
   !    fcinv=mub/mry
   !
   !    dym_rate=0.0_dblprec
   !    list_indx=0
   !    count_num=0
   !    ene_curr=0.0_dblprec
   !    cumu_barr=0.0_dblprec
   !    dym_ene_barr=0.0_dblprec
   !    Efield_ene=0.0_dblprec
   !    mdot=0.0_dblprec
   !
   !    ! Calculation of the magnetic energy for the polaron site
   !    do iatom=1, nlistsize(prev_polaron)
   !
   !       dnorm=sqrt((coord(1,nlist(iatom,prev_polaron))-coord(1,prev_polaron))**2+&
   !                  (coord(2,nlist(iatom,prev_polaron))-coord(2,prev_polaron))**2+&
   !                  (coord(3,nlist(iatom,prev_polaron))-coord(3,prev_polaron))**2)
   !       if (dnorm.le.1.0_dblprec.and.mod(dnorm,2.0_dblprec)==1.0_dblprec) then
   !       mdot=emomM(1,prev_polaron,1)*emomM(1,nlist(iatom,prev_polaron),1)+&
   !            emomM(2,prev_polaron,1)*emomM(2,nlist(iatom,prev_polaron),1)+&
   !            emomM(3,prev_polaron,1)*emomM(3,nlist(iatom,prev_polaron),1)
   !       ene_curr=ene_curr-ncoup(iatom,prev_polaron)*mdot
   !       endif
   !    enddo
   !
   !    ! Checking the neighbours of the current polaron site to see possible movement candidates
   !    do ene_barr=1, nlistsize_barriers(prev_polaron)
   !       current_atom=nlist_barriers(ene_barr,prev_polaron)
   !
   !       if ( atype(current_atom).ne. atype(prev_polaron) ) then
   !
   !         ! Looping through the neighbours of the possible movement candidates
   !         do iatom=1, nlistsize(current_atom)
   !            current_neigh=nlist(iatom,current_atom)
   !
   !            dnorm=sqrt((coord(1,current_neigh)-coord(1,current_atom))**2+&
   !                       (coord(2,current_neigh)-coord(2,current_atom))**2+&
   !                       (coord(3,current_neigh)-coord(3,current_atom))**2)
   !
   !            if (dnorm.le.1.0_dblprec.or.mod(dnorm,2.0_dblprec)==1.0_dblprec) then
   !            mdot=emomM(1,current_neigh,1)*emomM(1,current_atom,1)+&
   !                 emomM(2,current_neigh,1)*emomM(2,current_atom,1)+&
   !                 emomM(3,current_neigh,1)*emomM(3,current_atom,1)
   !            cumu_barr(ene_barr)=cumu_barr(ene_barr)-j_fm_fm*mdot
   !
   !            endif
   !         enddo
   !
   !         dr(1)= coord(1,current_atom)-coord(1,prev_polaron)
   !         dr(2)= coord(2,current_atom)-coord(2,prev_polaron)
   !         dr(3)= coord(3,current_atom)-coord(3,prev_polaron)
   !
   !         ! norm of the vector
   !         dnorm=sqrt(dr(1)**2+dr(2)**2+dr(3)**2)
   !
   !         ! check if the electric field is active$
   !         if (do_efield=='Y') then
   !            Efield_ene=efield(1)*dr(1)+efield(2)*dr(2)+efield(3)*dr(3)
   !         endif
   !
   !         ! Magnetic enerrgy difference, efield potential and lattice barrier
   !         dym_ene_barr(ene_barr)=((cumu_barr(ene_barr)-ene_curr)*fcinv*13.6/1000) + Lat_barr+Efield_ene
   !         ! Dynamical rates taking into account the magnetic configuration
   !
   !       endif
   !    enddo
   !
   !    if (dym_ene_barr(1).ne.0.0_dblprec) then
   !       dym_rate(1)=rate0*exp(-dym_ene_barr(1)/(k_bolt_ev*Temp_array(prev_polaron)))
   !    else
   !       dym_rate(1)=0.0_dblprec
   !    endif
   !
   !    do ene_barr=2,max_no_neigh/2
   !       if (dym_ene_barr(ene_barr).eq.0.0_dblprec) then
   !           dym_rate(ene_barr)=dym_rate(ene_barr-1)
   !       else
   !           dym_rate(ene_barr)=dym_rate(ene_barr-1)+rate0*exp(-dym_ene_barr(ene_barr)/(k_bolt_ev*Temp_array(polaron_seed)))
   !
   !       endif
   !     enddo
   !
   !    call rng_uniform(rrand,2)
   !
   !    ! Random number to check which state is chosen
   !    urand=(rrand(1))*maxval(dym_rate)
   !    state=maxloc(dym_rate,1)
   !    ! Check which one of the neighbouring states are the ones chosen
   !    do ii=1,max_no_neigh_barriers-1
   !
   !       if (count_num.lt.1) then
   !          if (dym_ene_barr(ii).gt.0) then
   !
   !             if ( (Dym_Rate(ii).lt.urand).and.(Dym_Rate(ii+1).ge.urand)) then
   !                ! Selected state from the partial summations
   !                next_polaron=nlist_barriers(ii,prev_polaron)
   !                count_num=count_num+1
   !              else
   !                 next_polaron=nlist_barriers(max_no_neigh_barriers,prev_polaron)
   !                 count_num=count_num+1
   !             endif
   !
   !          else
   !             next_polaron=prev_polaron
   !          endif
   !
   !       endif
   !    enddo
   !
   !    trand=(rrand(2)+1e-15)
   !
   !    new_polaron_pos=next_polaron
   !    KMC_time=(1.0_dblprec/Dym_rate(max_no_neigh_barriers))*log(1.0_dblprec/trand)
   !    ! Transform it to actual time step
   !    ! Typecast the time as integer ammount of time steps
   !    time_jump=floor(KMC_time/delta_t)
   !
   !  end subroutine dynamical_barrier


   !!------------------------------------------------------------------------------
   !! SUBROUTINE: linear_exchange_update
   !!
   !! DESCRIPTION:
   !!> This routine will linearly change the exchange interaction between FM and AFM
   !!> configurations, the objective behind this is to try to capture the coupled motion
   !!> of the magnetic polaron.
   !!> Ths time taken to perform such motion is calculated via kmc_method 1.
   !!
   !!> This is known as kmc_method 3
   !!
   !!> @author
   !!  Jonathan Chico
   !!------------------------------------------------------------------------------$
   !  subroutine linear_exchange_update(mstep,Nstep,Natom,Mensemble,max_no_neigh,KMC_STATIC_TIMES,&
   !             init_time,next_polaron,prev_polaron,atype,nlistsize,nlist,coord,mmom,&
   !             emom,emomM,ncoup)
   !
   !     implicit none
   !
   !     integer, intent(in) :: mstep
   !     integer, intent(in) :: Nstep
   !     integer, intent(in) :: Natom
   !     integer, intent(in) :: Mensemble
   !     integer, intent(in) :: init_time
   !     integer, intent(in) :: max_no_neigh
   !     integer, intent(in) :: next_polaron ! Next position of the polaron after jump
   !     integer, intent(in) :: prev_polaron ! Previous position of the polaron
   !     integer, intent(in) :: KMC_STATIC_TIMES
   !     integer, dimension(Natom), intent(inout) :: atype
   !     integer, dimension(Natom), intent(in) :: nlistsize
   !     integer, dimension(max_no_neigh,Natom), intent(in) :: nlist
   !
   !     real(dblprec), dimension(3,Natom), intent(in) :: coord
   !     real(dblprec), dimension(Natom,Mensemble), intent(inout) :: mmom
   !     real(dblprec), dimension(3,Natom,Mensemble), intent(inout) :: emom
   !     real(dblprec), dimension(3,Natom,Mensemble), intent(inout) :: emomM
   !     real(dblprec), dimension(max_no_neigh,Natom), intent(inout) :: ncoup
   !
   !     ! ... Local variables
   !     integer :: curr_time,rcoord,ii,curr_neigh
   !     integer :: iatom_prev, jatom_prev, neigh,neigh_cent,jatom_next
   !
   !     curr_time=Nstep-init_time
   !
   !!     The problem is then that the atom types need to be canged, and the atoms of type "one" when interacting
   !!     with type two need to be changed, after that the atom types themselves need to be changed
   !!
   !     ! Update the exchange interactions as a function of time
   !     ! Change the exchnage of the atoms that belong to the polaron, i.e. type 2 to being AFM
   !
   !     if (Nstep.eq.init_time) then
   !        ! Find the "link" between the polaron sites
   !        !$omp parallel do default(shared) private(iatom_prev,jatom_next,polaron_link)
   !        do iatom_prev=1, nlistsize(prev_polaron)
   !           do jatom_next=1, nlistsize(next_polaron)
   !
   !              if (nlist(iatom_prev,prev_polaron).eq.nlist(jatom_next,next_polaron)) then
   !                 polaron_link=nlist(iatom_prev,prev_polaron)
   !              endif
   !
   !           enddo
   !        enddo
   !       !$omp end parallel do
   !     endif
   !
   !    ! Also the magnetic moments need to be updated
   !    mmom(prev_polaron,1:Mensemble)=mom_fm+(mom_afm-mom_fm)*curr_time/KMC_STATIC_TIMES
   !
   !    ! This will update the center of the polaron
   !    do neigh=1, nlistsize(prev_polaron)
   !
   !       rcoord=(coord(1,nlist(neigh,prev_polaron))-coord(1,prev_polaron))**2+&
   !              (coord(2,nlist(neigh,prev_polaron))-coord(2,prev_polaron))**2+&
   !              (coord(3,nlist(neigh,prev_polaron))-coord(3,prev_polaron))**2
   !
   !       ! Update for all the neighbours, except for the link
   !       if (rcoord.le.1.0_dblprec.and.nlist(neigh,prev_polaron).ne.polaron_link) then
   !          ncoup(neigh,prev_polaron)=j_fm_fm+(j_afm_afm-j_fm_fm)*curr_time/KMC_STATIC_TIMES
   !       ! Update for the polaron link
   !       else if (rcoord.le.1.0_dblprec.and. nlist(neigh,prev_polaron).eq.polaron_link) then
   !          ! It goes from full FM to the interaction AMF with the FM region
   !          mmom(nlist(neigh,prev_polaron),1:Mensemble)=mom_fm+(mom_afm-mom_fm)*curr_time/KMC_STATIC_TIMES
   !          ncoup(neigh,prev_polaron)=j_fm_fm+(j_fm_afm-j_fm_fm)*curr_time/KMC_STATIC_TIMES
   !       else
   !          ncoup(neigh,prev_polaron)=0.0_dblprec
   !       endif
   !
   !    enddo
   !
   !    ! This will update the neighbours of the center of the polaron
   !    do neigh_cent=1, nlistsize(prev_polaron)
   !
   !       curr_neigh=nlist(neigh_cent,prev_polaron)
   !       ! Now update the exchange of the neighbours
   !        do neigh=1, nlistsize(curr_neigh)
   !
   !           rcoord=(coord(1,nlist(neigh,curr_neigh))-coord(1,curr_neigh))**2+&
   !                  (coord(2,nlist(neigh,curr_neigh))-coord(2,curr_neigh))**2+&
   !                  (coord(3,nlist(neigh,curr_neigh))-coord(3,curr_neigh))**2
   !
   !           ! Check that the current site is not the linking site
   !           if ( rcoord.le.1.0_dblprec .and.curr_neigh.ne.polaron_link) then
   !             ! Changing from FM to AFM interaction
   !             if (atype(nlist(neigh,curr_neigh)).eq.2)  then
   !                ncoup(neigh,curr_neigh)=j_fm_fm+(j_afm_afm-j_fm_fm)*curr_time/KMC_STATIC_TIMES
   !             ! Changing from border AFM to AFM
   !             else if (atype(nlist(neigh,curr_neigh)).eq.1) then
   !                ncoup(neigh,curr_neigh)=j_fm_afm+(j_afm_afm-j_fm_afm)*curr_time/KMC_STATIC_TIMES
   !             endif
   !           ! If the current site is the linking site
   !           else if ( rcoord.le.1.0_dblprec .and.curr_neigh.eq.polaron_link) then
   !               ! Changing from FM to border AFM
   !               if (nlist(neigh,curr_neigh).eq. prev_polaron) then
   !                  ncoup(neigh,curr_neigh)=j_fm_fm+(j_fm_afm-j_fm_fm)*curr_time/KMC_STATIC_TIMES
   !               ! Changing from border AFM to FM
   !             else if (nlist(neigh,curr_neigh).eq. next_polaron) then
   !                  ncoup(neigh,curr_neigh)=j_fm_afm+(j_fm_fm-j_fm_afm)*curr_time/KMC_STATIC_TIMES
   !              endif
   !           endif
   !
   !        enddo
   !
   !    enddo
   !
   !    ! Updating the moment for the new polaron center
   !    mmom(next_polaron,1:Mensemble)=mom_afm+(mom_fm-mom_afm)*curr_time/KMC_STATIC_TIMES
   !
   !    ! Change the exchange of the atoms that belong to the new polaron site to be FM
   !    do neigh=1, nlistsize(next_polaron)
   !
   !       rcoord=(coord(1,nlist(neigh,next_polaron))-coord(1,next_polaron))**2+&
   !              (coord(2,nlist(neigh,next_polaron))-coord(2,next_polaron))**2+&
   !              (coord(3,nlist(neigh,next_polaron))-coord(3,next_polaron))**2
   !
   !       ! Change the interaction between the new polaron center and the link site
   !       if (rcoord.le.1.0_dblprec.and.nlist(neigh,next_polaron).eq.polaron_link) then
   !           ncoup(neigh,next_polaron)=j_fm_afm+(j_fm_fm-j_fm_afm)*curr_time/KMC_STATIC_TIMES
   !       ! Change the interaction between the new polaron center and its other neighbours
   !     else if (rcoord.le.1.0_dblprec.and.nlist(neigh,next_polaron).ne.polaron_link) then
   !           ncoup(neigh,next_polaron)=j_afm_afm+(j_fm_fm-j_afm_afm)*curr_time/KMC_STATIC_TIMES
   !           mmom(nlist(neigh,next_polaron),1:Mensemble)=mom_afm+(mom_fm-mom_afm)*curr_time/KMC_STATIC_TIMES
   !       endif
   !
   !    enddo
   !
   !    ! Then must change the neighbour
   !    do neigh=1, nlistsize(next_polaron)
   !       curr_neigh=nlist(neigh,next_polaron)
   !
   !       ! Changing the interactions between the polaron link and its neighbours
   !       if (curr_neigh.eq.polaron_link) then
   !
   !          do ii=1, nlistsize(curr_neigh)
   !
   !             ! Bonding distance between the neighbours
   !             rcoord=(coord(1,nlist(ii,curr_neigh))-coord(1,curr_neigh))**2+&
   !                    (coord(2,nlist(ii,curr_neigh))-coord(2,curr_neigh))**2+&
   !                    (coord(3,nlist(ii,curr_neigh))-coord(3,curr_neigh))**2
   !             ! Change the interaction between the neighbours and the new polaron center
   !             if  (rcoord.le.1.0_dblprec.and.nlist(ii,curr_neigh).eq.next_polaron) then
   !                 ncoup(ii,curr_neigh)=j_fm_afm+(j_fm_fm-j_fm_afm)*curr_time/KMC_STATIC_TIMES
   !             ! Change the interactions between the neighbours and the background
   !             else if (rcoord.le.1.0_dblprec.and.nlist(ii,curr_neigh).ne.prev_polaron) then
   !                 ncoup(ii,curr_neigh)=j_fm_fm+(j_fm_afm-j_fm_fm)*curr_time/KMC_STATIC_TIMES
   !             endif
   !
   !          enddo
   !
   !       else
   !
   !          do ii=1, nlistsize(curr_neigh)
   !             ! Bonding distance between the neighbours
   !             rcoord=(coord(1,nlist(ii,curr_neigh))-coord(1,curr_neigh))**2+&
   !                    (coord(2,nlist(ii,curr_neigh))-coord(2,curr_neigh))**2+&
   !                    (coord(3,nlist(ii,curr_neigh))-coord(3,curr_neigh))**2
   !             ! Change the interaction between the neighbours and the new polaron center
   !             if (rcoord.le.1.0_dblprec.and.nlist(ii,curr_neigh).eq.next_polaron) then
   !                 ncoup(ii,curr_neigh)=j_afm_afm+(j_fm_fm-j_afm_afm)*curr_time/KMC_STATIC_TIMES
   !             ! Change the interactions between the neighbours and the background
   !           else if (rcoord.le.1.0_dblprec.and.nlist(ii,curr_neigh).ne.next_polaron) then
   !                 ncoup(ii,curr_neigh)=j_afm_afm+(j_fm_afm-j_afm_afm)*curr_time/KMC_STATIC_TIMES
   !             endif
   !          enddo
   !
   !       endif
   !
   !    enddo
   !
   !    if (curr_time.eq.KMC_STATIC_TIMES) then
   !       ! Change the type of the atom in which the polaron is localized now to type 2
   !       atype=1
   !       atype(next_polaron)=2
   !       do neigh=1, nlistsize(next_polaron)
   !          atype(neigh)=2
   !       enddo
   !    endif
   !
   !
   !  end subroutine linear_exchange_update
   !
   !! For Nina's angle approach
   !!------------------------------------------------------------------------------
   !! SUBROUTINE: anti_parallel_update
   !!
   !! DESCRIPTION:
   !!> The magnetic polaron is allowed to move if the spins are in-plane and
   !!> anti-parallel to each other. There is a quantization axis assumed to be in the
   !!> z-direction.
   !!
   !!> This is known as kmc_method 4
   !!
   !!> @author
   !!  Jonathan Chico
   !!------------------------------------------------------------------------------$
   !  subroutine anti_parallel_update(Natom,Mensemble,curr_polaron,max_no_neigh,max_no_neigh_barriers,&
   !             atype,nlistsize,nlistsize_barriers,nlist,nlist_barriers,coord,mmom,emom,emomM,&
   !             ncoup,do_ralloy,Natom_full,achem_ch)
   !
   !     use Constants
   !
   !     implicit none
   !
   !     integer, intent(in) :: Natom
   !     integer, intent(in) :: Natom_full
   !     integer, intent(in) :: Mensemble
   !     integer, intent(in) :: curr_polaron
   !     integer, intent(in) :: do_ralloy
   !     integer, intent(in) :: max_no_neigh
   !     integer, intent(in) :: max_no_neigh_barriers
   !     integer, dimension(Natom), intent(inout) :: atype
   !     integer, dimension(Natom_full), intent(inout) :: achem_ch
   !     integer, dimension(Natom), intent(in) :: nlistsize
   !     integer, dimension(Natom), intent(in) :: nlistsize_barriers
   !     integer, dimension(max_no_neigh,Natom), intent(in) :: nlist
   !     integer, dimension(max_no_neigh,Natom), intent(in) :: nlist_barriers
   !     real(dblprec), dimension(3,Natom), intent(in) :: coord
   !     real(dblprec), dimension(Natom,Mensemble), intent(inout) :: mmom
   !     real(dblprec), dimension(3,Natom,Mensemble), intent(inout) :: emom
   !     real(dblprec), dimension(3,Natom,Mensemble), intent(inout) :: emomM
   !     real(dblprec), dimension(max_no_neigh,Natom), intent(inout) :: ncoup
   !
   !     ! ... Local variables
   !     integer :: neigh,neigh_atom
   !     real(dblprec) :: neigh_angle, center_angle, diff_angle,rnorm,TOL
   !     real(dblprec), dimension(3) :: rcoord
   !
   !     rnorm=0.0_dblprec
   !     rcoord=0.0_dblprec
   !     diff_angle=0.0_dblprec
   !     neigh_angle=0.0_dblprec
   !     center_angle=0.0_dblprec
   !
   !     TOL=pi*0.10
   !
   !     do neigh=1, nlistsize_barriers(curr_polaron)
   !        ! Considered neighbour dummy
   !        neigh_atom=nlist_barriers(neigh,curr_polaron)
   !        ! Norm of the bonding vector
   !        rnorm=sqrt((coord(1,neigh_atom)-coord(1,curr_polaron))**2+&
   !                   (coord(2,neigh_atom)-coord(2,curr_polaron))**2+&
   !                   (coord(3,neigh_atom)-coord(3,curr_polaron))**2)
   !
   !        rcoord(1)=(coord(1,neigh_atom)-coord(1,curr_polaron))/rnorm
   !        rcoord(2)=(coord(2,neigh_atom)-coord(2,curr_polaron))/rnorm
   !        rcoord(3)=(coord(3,neigh_atom)-coord(3,curr_polaron))/rnorm
   !
   !        neigh_angle=acos(emom(1,neigh_atom,1)*rcoord(1)+emom(2,neigh_atom,1)*rcoord(2)+&
   !                         emom(3,neigh_atom,1)*rcoord(3))
   !
   !        center_angle=acos(emom(1,curr_polaron,1)*rcoord(1)+emom(2,curr_polaron,1)*rcoord(2)+&
   !                          emom(3,curr_polaron,1)*rcoord(3))
   !
   !        diff_angle=pi-neigh_angle
   !
   !        ! If the angle is below the tolerance accept the motion of the polaron
   !        if ( abs(neigh_angle).lt. TOL .and.  abs(diff_angle).lt. TOL ) then
   !
   !           new_polaron_pos=neigh_atom
   !           ! Update the neighbour list
   !           call update_neigh(Natom,Mensemble,max_no_neigh,new_polaron_pos,&
   !                nlist,nlistsize,atype,coord,mmom,emom,emomM,ncoup,&
   !                do_ralloy,achem_ch,Natom_full)
   !
   !        endif
   !
   !     enddo
   !
   !  end subroutine anti_parallel_update

end module KMC
