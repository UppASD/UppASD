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
   function NumProcs() result(nprocs) bind(c,name='numprocs_')
      use iso_c_binding
           use uppasd, only : number_of_active_processors
      implicit none
      integer :: nprocs

      nprocs = number_of_active_processors()
      return
   end function NumProcs



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
      use InputData, only : ipmode, nstep => ipnstep, temperature => iptemp, timestep => ipdelta_t, ipnphase
      use Damping, only : damping1 => lambda1_array
      use MomentData, only : emomM, emom, mmom
      use FieldData, only : beff
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
 

      !ipnphase = 1
      !   ipmode = imode
      !   if (allocated(damping)) deallocate(damping)
      !      allocate(damping(1))
      !    damping = idamping
      !   if (allocated(nstep)) deallocate(nstep)
      !      allocate(nstep(1))
      !    nstep = instep
      !    if (allocated(temperature)) deallocate(temperature)  
      !      allocate(temperature(1))
      !    temperature = itemperature
      !    if (allocated(timestep)) deallocate(timestep)
      !      allocate(timestep(1))
      !    timestep = itimestep

      call timing(0,'Initial       ','ON')
      if(imode == 'M' .or. imode == 'H') then
         !call mc_iphase()
         call mc_minimal(emomM,emom,mmom, 1000, 'M ', 10.0d0)
      else
         !call sd_iphase()
         damping1 = idamping
         call sd_minimal(emomM,emom,mmom, 1000, 1, 10.0d0)
      end if
      call timing(0,'Initial       ','OF')

      moments = emomM

      !!! nstep = nstep_old
      !!! temperature = temperature_old
      !!! timestep = timestep_old
      !!! ipmode = ipmode_old

      return

      
    end subroutine Relax
!!! 
!!!    subroutine RelaxMD()
!!!       implicit none
!!! 
!!!       call ld_iphase() 
!!!    end subroutine RelaxMD 
!!! 
!!!    subroutine RelaxSLDLLG()
!!!       implicit none 
!!! 
!!!       call sld_iphase() 
!!!    end subroutine RelaxSLDLLG 
!!! 
!!! 
!!!    !==============================================================!
!!!    ! Moment handling routines
!!!    !--------------------------------------------------------------!
!!! 
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
         use MomentData, only : emom
       implicit none
       !f2py intent(in) :: natom
       integer(c_int), intent(in) :: natom
       !f2py intent(in) :: mensemble
       integer(c_int), intent(in) :: mensemble
       !f2py intent(in) moments
       real(c_double), dimension(3,natom, mensemble), intent(in) :: moments

       emom = moments 
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

    subroutine get_mcnstep(nstep) bind(c, name='get_mcnstep_')
      use iso_c_binding
         use InputData, only : mcnstep
       implicit none
       !f2py intent(out) nstep
       integer(c_int), intent(out) :: nstep

       nstep = mcnstep
    end subroutine get_mcnstep

    subroutine get_temperature(temperature) bind(c, name='get_temperature_')
      use iso_c_binding
         use InputData, only : temp
       implicit none
       !f2py intent(out) nstep
       real(c_double), intent(out) :: temperature

       temperature = temp
    end subroutine get_temperature

    subroutine get_delta_t(timestep) bind(c, name='get_delta_t_')
      use iso_c_binding
         use InputData, only : delta_t
       implicit none
       !f2py intent(out) nstep
       real(c_double), intent(out) :: timestep

       timestep = delta_t
    end subroutine get_delta_t

!!!    !==============================================================!
!!!    ! Measurement Runners
!!!    !--------------------------------------------------------------!
!!! 
!!!    subroutine RelaxGNEB()
!!!       implicit none
!!! 
!!!       call em_iphase()
!!!    end subroutine RelaxGNEB
!!! 
!!!    subroutine RunMonteCarlo()
!!!       implicit none
!!! 
!!!       call mc_mphase() 
!!!    end subroutine RunMonteCarlo
!!! 
!!!    subroutine RunMultiScale()
!!!       implicit none
!!! 
!!!       call ms_mphase()
!!!    end subroutine RunMultiScale
!!! 
!!!    subroutine RunLLGlite()
!!!       implicit none
!!! 
!!!       call sd_mphase_lite() 
!!!    end subroutine RunLLGlite
!!! 
!!!    subroutine RunLLG()
!!!       implicit none
!!! 
!!!       call sd_mphase() 
!!!    end subroutine RunLLG
!!! 
!!!    subroutine RunLLGCUDA()
!!!       implicit none
!!! 
!!!       call sd_mphaseCUDA()
!!!    end subroutine RunLLGCUDA
!!! 
!!!    subroutine RunLD()
!!!       implicit none
!!! 
!!!       call ld_mphase()
!!!    end subroutine RunLD
!!! 
!!!    subroutine RunSLDLLG()
!!!       implicit none
!!! 
!!!       call sld_mphase()
!!!    end subroutine RunSLDLLG
!!! 
!!!    subroutine RunSLDLLGimplicit()
!!!       implicit none
!!! 
!!!       call sld_mphase3()
!!!    end subroutine RunSLDLLGimplicit
!!! 
!!!    subroutine RunGNEB()
!!!       implicit none
!!! 
!!!       call mep_mphase()
!!!    end subroutine RunGNEB
!!! 
function TotalEnergy() result(energy) bind(c, name='totalenergy_')
   use iso_c_binding
   use uppasd, only : calculate_energy
   implicit none

   real(c_double) :: energy

   call calculate_energy(energy)
   return

end function totalenergy
!!! 
!!!    !!! function getMoments() result(moments)
!!!    !!!    use MomentData
!!!    !!!    use InputData, only : Natom, Mensemble
!!!    !!!    implicit none
!!!    !!!
!!!    !!!    !f2py integer, parameter : Natom
!!!    !!!    !f2py integer, parameter : Mensemble
!!!    !!!    real(dblprec), dimension(3,Natom,Mensemble) :: moments
!!!    !!!
!!!    !!!    moments=emom
!!!    !!!    return 
!!!    !!! end function getMoments

!!! end module pyasd
