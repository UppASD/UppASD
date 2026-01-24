!******************************************************************
!*                     *** UppASD ***                             *
!*                                                                *
!*               Uppsala Atomic Spin Dynamics                     *
!*                                                                *
!*                   Version 5.0 Mar 2017                         *
!*                                                                *
!*       Anders Bergman                                           *
!*       Johan Hellsvik             Lars Bergqvist                *
!*       Jonathan Chico             Thomas Nystrand               *
!*       Bjorn Skubic               Olle Eriksson                 *
!*       Andrea Taroni              Lars Nordstrom                *
!*       + other contributors.                                    *
!*                                                                *
!*       Division of Materials Theory                             *
!*       Department of Physics and Materials Science              *
!*       Uppsala University                                       *
!*       Sweden                                                   *
!*                                                                *
!*       and                                                      *
!*                                                                *
!*       Department of Materials and Nano Science                 *
!*       KTH Royal Institute of Technology                        *
!*       Sweden                                                   *
!*                                                                *
!******************************************************************


!------------------------------------------------------------------------------------
! MODULE: pyasd
!
! DESCRIPTION:
!> Main program for atomistic spin dynamic simulations.
!> The program performs atomistic spin dynamics simulations
!> using the Landau-Lifshitz equation. The methodology and capabilties are described
!> in the textbook "Atomistic spin dynamics - Foundations and Applications" by
!> O. Eriksson, A. Bergman, L. Bergqvist and J. Hellsvik
!> Oxford University Press, Oxford, UK, 2017
!> ISBN: 9780198788669
!> https://global.oup.com/academic/product/atomistic-spin-dynamics-9780198788669
!
!------------------------------------------------------------------------------------
!module pyasd
!   !
!   use iso_c_binding
!   use omp_lib
!   use uppasd
!   use MomentData
!   use InputData
!
!   !.. Implicit declarations
!   implicit none
!
!   ! Variables for timing and number of processors
!
!   ! Local array for spins (for python reasons), low precision for now
!
!contains
!
   subroutine RunUppASD() bind(c, name='runuppasd_')
   use uppasd, only : main
      implicit none

      call main()
   end subroutine RunUppASD
   !==============================================================!
   ! Check if inpsd.dat exists and whether it contains anything
   !--------------------------------------------------------------!
   subroutine SanityCheck() bind(c, name='sanitycheck_')
           use ErrorHandling, only : ErrorHandling_check_input_file

      implicit none

      call ErrorHandling_check_input_file()
   end subroutine SanityCheck

   !==============================================================!
   ! Find number of OpenMP processors active
   !--------------------------------------------------------------!
   subroutine get_numprocs(nprocs) bind(c,name='get_numprocs_')
      use iso_c_binding
      use uppasd, only : number_of_active_processors

      implicit none
      !f2py intent(out) :: nprocs
      integer(c_int), intent(out) :: nprocs

      nprocs = number_of_active_processors()
   end subroutine get_numprocs

   !==============================================================!
   ! Initialize profiling routines and print the UppASD logo
   !--------------------------------------------------------------!
   subroutine PrintLogo() bind(c, name='printlogo_')
           use uppasd, only : print_logo
      implicit none

      call print_logo()
   end subroutine PrintLogo
   !==============================================================!
   ! Setup all things needed for a simulation
   !--------------------------------------------------------------!
   subroutine SetupAll(nat, men) bind(c, name='setupall_')
      use iso_c_binding
           use uppasd, only : setup_simulation
           use InputData, only : natom, mensemble
           use energy, only : allocate_energies
           use MonteCarlo, only : allocate_mcdata
             use RandomNumbers, only : allocate_randomwork
      use Depondt,               only : allocate_depondtfields
      use Midpoint,              only : allocate_aux_midpoint_fields
      implicit none
      !f2py intent(out) :: nat
      integer(c_int), intent(out) :: nat
      !f2py intent(out) :: men
      integer(c_int), intent(out) :: men
      
      call setup_simulation()
      call allocate_energies(1, Mensemble)
      call allocate_randomwork(Natom,Mensemble,1,'N')
      call allocate_depondtfields(Natom, Mensemble,1)
      call allocate_aux_midpoint_fields(1,Natom,Mensemble)
      call allocate_mcdata(Natom,1)
      men = mensemble
      nat = natom
   end subroutine SetupAll

   !==============================================================!
   ! Run an initial phase simulation 
   !--------------------------------------------------------------!
   subroutine InitialPhase() bind(c, name='initialphase_')
           use uppasd, only : run_initial_phase
      implicit none

      call run_initial_phase()
   end subroutine InitialPhase

   !==============================================================!
   ! Run an measurement phase simulation 
   !--------------------------------------------------------------!
   subroutine Measure() bind(c,name='measure_')
           use uppasd, only : run_measurement_phase
      implicit none
      call run_measurement_phase()
   end subroutine Measure


   !==============================================================!
   ! Clean up after simulation
   !--------------------------------------------------------------!
   subroutine CleanUp() bind(c,name='cleanup_')
           use uppasd, only : cleanup_simulation
           use energy, only : allocate_energies
           use MonteCarlo, only : allocate_mcdata
             use RandomNumbers, only : allocate_randomwork
      implicit none
      call cleanup_simulation()
      call allocate_energies(-1)
      call allocate_mcdata(1,-1)
      call allocate_randomwork(1,1,-1,'N')
   end subroutine CleanUp


!!!    !==============================================================!
!!!    ! Relaxation Runners
!!!    !--------------------------------------------------------------!
!!! 
    subroutine RelaxMonteCarlo(moments, natom, mensemble) bind(c, name='relaxmontecarlo_')
      use iso_c_binding
         use uppasd, only : mc_iphase
         use MomentData, only : emom
      use Profiling, only : timing
       implicit none
       !f2py intent(in) :: natom
       integer(c_int), intent(in) :: natom
       !f2py intent(in) :: mensemble
       integer(c_int), intent(in) :: mensemble
       !f2py intent(out) moments
       real(c_double), dimension(3,natom, mensemble), intent(out) :: moments
 
      call timing(0,'Initial       ','ON')
       call mc_iphase()
      call timing(0,'Initial       ','OF')
       moments = emom
    end subroutine RelaxMonteCarlo
 
    subroutine RelaxMetropolis(moments, natom, mensemble) bind(c, name='relaxmetropolis_')
      use iso_c_binding
         use InputData, only : ipmode
         use uppasd, only : mc_iphase
         use MomentData, only : emom
      use Profiling, only : timing
       implicit none
       !f2py intent(in) :: natom
       integer(c_int), intent(in) :: natom
       !f2py intent(in) :: mensemble
       integer(c_int), intent(in) :: mensemble
       !f2py intent(out) moments
       real(c_double), dimension(3,natom, mensemble), intent(out) :: moments
 
      call timing(0,'Initial       ','ON')
       ipmode='M'
       call mc_iphase()
      call timing(0,'Initial       ','OF')
       moments = emom
    end subroutine RelaxMetropolis
 
    subroutine RelaxHeatBath(moments, natom, mensemble) bind(c, name='relaxheatbath_')
      use iso_c_binding
       use InputData, only : ipmode
         use uppasd, only : mc_iphase
         use MomentData, only : emom
      use Profiling, only : timing
       implicit none
       !f2py intent(in) :: natom
       integer(c_int), intent(in) :: natom
       !f2py intent(in) :: mensemble
       integer(c_int), intent(in) :: mensemble
       !f2py intent(out) moments
       real(c_double), dimension(3,natom, mensemble), intent(out) :: moments

      call timing(0,'Initial       ','ON')
       ipmode='H'
       call mc_iphase()
      call timing(0,'Initial       ','OF')
       moments = emom
    end subroutine RelaxHeatBath

!!! 
!!!    subroutine RelaxMultiScale()
!!!       implicit none
!!! 
!!!       call ms_iphase()
!!!    end subroutine RelaxMultiScale
!!! 
    subroutine RelaxLLG(moments, natom, mensemble) bind(c, name='relaxllg_')
      use iso_c_binding
         use uppasd, only : sd_iphase
         use MomentData, only : emom
      use Profiling, only : timing
       implicit none
       !f2py intent(in) :: natom
       integer(c_int), intent(in) :: natom
       !f2py intent(in) :: mensemble
       integer(c_int), intent(in) :: mensemble
       !f2py intent(out) moments
       real(c_double), dimension(3,natom, mensemble), intent(out) :: moments
  
      call timing(0,'Initial       ','ON')
       call sd_iphase()
      call timing(0,'Initial       ','OF')
       moments = emom
    end subroutine RelaxLLG

    subroutine Relax(moments, natom, mensemble, imode, instep, itemperature, itimestep, idamping)  bind(c, name='relax_')
      use iso_c_binding
      use uppasd, only : sd_iphase, mc_iphase
      use sd_driver, only : sd_minimal
      use mc_driver, only : mc_minimal
      use Damping, only : damping1 => lambda1_array
      use MomentData, only : emomM, emom, mmom, emom2
      use Profiling, only : timing

      implicit none
      !f2py intent(out) moments
      real(c_double), dimension(3,natom, mensemble), intent(out) :: moments
      !f2py intent(in) :: natom
      integer(c_int), intent(in) :: natom
      !f2py intent(in) :: mensemble
      integer(c_int), intent(in) :: mensemble
      !f2py intent(in) :: imode
      character(kind=c_char), intent(in) :: imode
      !f2py intent(in) :: instep
      integer(c_int), intent(in) :: instep
      !f2py intent(in) :: itemperature
      real(c_double), intent(in) :: itemperature
      !f2py intent(in) :: itimestep
      real(c_double), intent(in) :: itimestep
      !f2py intent(in) :: idamping
      real(c_double), intent(in) :: idamping

      call timing(0,'Initial       ','ON')
      
      if(imode == 'M' .or. imode == 'H') then
         call mc_minimal(emomM, emom, mmom, instep, imode//" ", itemperature)
      else
         damping1 = idamping
         call sd_minimal(emomM, emom, mmom, instep, 1, itemperature)
      end if
      call timing(0,'Initial       ','OF')

      ! Copy data directly to output array
      ! F2PY allocates moments, we just fill it
      moments = emomM

      return
      
    end subroutine Relax
!!! 
!!!    !==============================================================!
!!!    ! Moment handling routines
!!!    !--------------------------------------------------------------!
!!! 
    subroutine get_coord(coords, natom) bind(c, name='get_coord_')
      use iso_c_binding
         use SystemData, only : icoord => coord
       implicit none
       !f2py intent(in) :: natom
       integer(c_int), intent(in) :: natom
       !f2py intent(out) coords
       real(c_double), dimension(3,natom), intent(out) :: coords

       coords = icoord
    end subroutine get_coord

    subroutine get_emom(moments, natom, mensemble) bind(c, name='get_emom_')
      use iso_c_binding
         use MomentData, only : emom
       implicit none
       !f2py intent(in) :: natom
       integer(c_int), intent(in) :: natom
       !f2py intent(in) :: mensemble
       integer(c_int), intent(in) :: mensemble
       !f2py intent(out) moments
       real(c_double), dimension(3,natom, mensemble), intent(out) :: moments

       moments = emom
    end subroutine get_emom

    subroutine put_emom(moments, natom, mensemble) bind(c, name='put_emom_')
      use iso_c_binding
         use MomentData, only : emom, emomM, mmom, emom2
       implicit none
       !f2py intent(in) :: natom
       integer(c_int), intent(in) :: natom
       !f2py intent(in) :: mensemble
       integer(c_int), intent(in) :: mensemble
       !f2py intent(in) moments
       real(c_double), dimension(3,natom, mensemble), intent(in) :: moments
       !
       integer(c_int) :: i_atom, i_ensemble

       emom = moments 
       emom2 = emom
       do i_ensemble = 1, mensemble
          do i_atom = 1, natom
            emomM(:,i_atom,i_ensemble) = moments(:,i_atom,i_ensemble) * mmom(i_atom,i_ensemble)
         end do
       end do
    end subroutine put_emom

!!!    !==============================================================!
!!!    ! Field handling routines
!!!    !--------------------------------------------------------------!
!!! 
    subroutine get_beff(fields, natom, mensemble) bind(c, name='get_beff_')
      use iso_c_binding
         use FieldData, only : beff
         use HamiltonianActions, only : effective_field
       implicit none
       !f2py intent(in) :: natom
       integer(c_int), intent(in) :: natom
       !f2py intent(in) :: mensemble
       integer(c_int), intent(in) :: mensemble
       !f2py intent(out) fields
       real(c_double), dimension(3,natom, mensemble), intent(out) :: fields

      call effective_field()
       fields = beff
    end subroutine get_beff

    subroutine put_beff(fields, natom, mensemble) bind(c, name='put_beff_')
      use iso_c_binding
         use FieldData, only : beff
       implicit none
       !f2py intent(in) :: natom
       integer(c_int), intent(in) :: natom
       !f2py intent(in) :: mensemble
       integer(c_int), intent(in) :: mensemble
       !f2py intent(in) moments
       real(c_double), dimension(3,natom, mensemble), intent(in) :: fields

       beff = fields
    end subroutine put_beff
!!!    !==============================================================!
!!!    ! Input data routines
!!!    !--------------------------------------------------------------!
!!! 
    subroutine get_nstep(nstep) bind(c, name='get_nstep_')
      use iso_c_binding
         use InputData, only : instep => nstep
       implicit none
       !f2py intent(out) nstep
       integer(c_int), intent(out) :: nstep

       nstep = instep
    end subroutine get_nstep

    subroutine put_nstep(nstep) bind(c, name='put_nstep_')
      use iso_c_binding
         use InputData, only : ostep => nstep
       implicit none
       !f2py intent(in) nstep
       integer(c_int), intent(in) :: nstep

       ostep = nstep
    end subroutine put_nstep

    subroutine get_hfield(hfield) bind(c, name='get_hfield_')
      use iso_c_binding
         use InputData, only : ihfield => hfield
       implicit none
       !f2py intent(out) hfield
       real(c_double), dimension(3), intent(out) :: hfield

       hfield = ihfield
    end subroutine get_hfield

    subroutine put_hfield(hfield) bind(c, name='put_hfield_')
      use iso_c_binding
         use InputData, only : ohfield => hfield
       implicit none
       !f2py intent(in) hfield
       real(c_double), dimension(3), intent(in) :: hfield

       ohfield = hfield
    end subroutine put_hfield

    subroutine get_iphfield(hfield) bind(c, name='get_iphfield_')
      use iso_c_binding
         use InputData, only : ihfield => iphfield
       implicit none
       !f2py intent(out) hfield
       real(c_double), dimension(3), intent(out) :: hfield

       hfield = ihfield
    end subroutine get_iphfield

    subroutine put_iphfield(hfield) bind(c, name='put_iphfield_')
      use iso_c_binding
         use InputData, only : ohfield => iphfield
       implicit none
       !f2py intent(in) hfield
       real(c_double), dimension(3), intent(in) :: hfield

       ohfield = hfield
    end subroutine put_iphfield


      subroutine get_ipmode(mode) bind(c, name='get_ipmode_')
         use iso_c_binding
             use InputData, only : ipmode
          implicit none
          !f2py intent(out) mode
          character(kind=c_char), dimension(2), intent(out) :: mode

          mode(1) = ipmode(1:1)
          mode(2) = ipmode(2:2)
      end subroutine get_ipmode


      subroutine put_ipmode(mode) bind(c, name='put_ipmode_')
         use iso_c_binding
             use InputData, only : ipmode
          implicit none
          !f2py intent(in) mode
          character(kind=c_char), dimension(2), intent(in) :: mode

          ipmode = mode(1)//mode(2)
      end subroutine put_ipmode


    subroutine get_mcnstep(nstep) bind(c, name='get_mcnstep_')
      use iso_c_binding
         use InputData, only : mcnstep
       implicit none
       !f2py intent(out) nstep
       integer(c_int), intent(out) :: nstep

       nstep = mcnstep
    end subroutine get_mcnstep

    subroutine put_mcnstep(nstep) bind(c, name='put_mcnstep_')
      use iso_c_binding
         use InputData, only : mcnstep
       implicit none
       !f2py intent(in) nstep
       integer(c_int), intent(in) :: nstep

       mcnstep = nstep
    end subroutine put_mcnstep

    subroutine get_temperature(temperature) bind(c, name='get_temperature_')
      use iso_c_binding
         use InputData, only : itemp => temp
       implicit none
       !f2py intent(out) temperature
       real(c_double), intent(out) :: temperature

       temperature = itemp
    end subroutine get_temperature

    subroutine put_temperature(temperature) bind(c, name='put_temperature_')
      use iso_c_binding
         use InputData, only : otemp => temp
       implicit none
       !f2py intent(in) temperature
       real(c_double), intent(in) :: temperature

       otemp = temperature
    end subroutine put_temperature

    subroutine get_iptemperature(temperature) bind(c, name='get_iptemperature_')
      use iso_c_binding
         use InputData, only : itemp => iptemp
       implicit none
       !f2py intent(out) temperature
       real(c_double), intent(out) :: temperature

       if (allocated(itemp)) then
          temperature = itemp(1)
       else
          temperature = 1.0d-6
       end if
    end subroutine get_iptemperature

    subroutine put_iptemperature(temperature) bind(c, name='put_iptemperature_')
      use iso_c_binding
         use InputData, only : otemp => iptemp
       implicit none
       !f2py intent(in) temperature
       real(c_double), intent(in) :: temperature

       otemp = temperature
    end subroutine put_iptemperature

    subroutine get_delta_t(timestep) bind(c, name='get_delta_t_')
      use iso_c_binding
         use InputData, only : delta_t
       implicit none
       !f2py intent(out) timestep
       real(c_double), intent(out) :: timestep

       timestep = delta_t
    end subroutine get_delta_t

    subroutine put_delta_t(timestep) bind(c, name='put_delta_t_')
      use iso_c_binding
         use InputData, only : delta_t
       implicit none
       !f2py intent(in) timestep
       real(c_double), intent(in) :: timestep

       delta_t = timestep
    end subroutine put_delta_t

   !!! 
   subroutine get_energy(energy) bind(c, name='get_energy_')
      use iso_c_binding
      use InputData, only : Natom, Mensemble
      use HamiltonianActions, only : effective_field
      implicit none

      !f2py intent(out) energy
      real(c_double), intent(out) :: energy

      call effective_field(energy)
      energy = energy  / (Natom * Mensemble)
   
   end subroutine get_energy
!!! 

!!! end module pyasd
