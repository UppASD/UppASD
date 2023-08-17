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


   !!==============================================================!
   !! Find number of OpenMP processors active
   !!--------------------------------------------------------------!
   !function NumProcs() result(nprocs) bind(c,name='numprocs_')
   !        use uppasd, only : number_of_active_processors
   !   implicit none
   !   integer :: nprocs

   !   nprocs = number_of_active_processors()
   !   return
   !end function NumProcs



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
   subroutine SetupAll() bind(c, name='setupall_')
           use uppasd, only : setup_simulation
      implicit none

      call setup_simulation()
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
      implicit none
      call cleanup_simulation()
   end subroutine CleanUp



!!!    !==============================================================!
!!!    ! Relaxation Runners
!!!    !--------------------------------------------------------------!
!!! 
!!!    subroutine RelaxMonteCarlo()
!!!       implicit none
!!! 
!!!       call mc_iphase()
!!!    end subroutine RelaxMonteCarlo
!!! 
!!!    subroutine RelaxMetropolis()
!!!       use InputData, only : ipmode
!!!       implicit none
!!! 
!!!       ipmode='M'
!!!       call mc_iphase()
!!!    end subroutine RelaxMetropolis
!!! 
!!!    subroutine RelaxHeatBath()
!!!       use InputData, only : ipmode
!!!       implicit none
!!! 
!!!       ipmode='H'
!!!       call mc_iphase()
!!!    end subroutine RelaxHeatBath
!!! 
!!!    subroutine RelaxMultiScale()
!!!       implicit none
!!! 
!!!       call ms_iphase()
!!!    end subroutine RelaxMultiScale
!!! 
!!!    subroutine RelaxLLG()
!!!       implicit none
!!! 
!!!       call sd_iphase()
!!!    end subroutine RelaxLLG
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
!!!    function TotalEnergy() result(energy)
!!!       implicit none
!!! 
!!!       real(dblprec) :: energy
!!! 
!!!       call calculate_energy(energy)
!!!       return
!!! 
!!!    end function totalenergy
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
