!******************************************************************
!*                     *** UppASD ***                             *
!*                                                                *
!*               Uppsala Atomic Spin Dynamics                     *
!*                                                                *
!*                   Version 6.0 Oct 2022                         *
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
! MODULE: uppasd
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
module uppasd
   !
   use omp_lib
   use iso_fortran_env
   use iso_c_binding
   use SD_driver
   use WL_driver
   use MC_driver
   use MS_driver
   use PT_driver
   use SX_driver
   use Profiling
   use SLD_driver
   use Parameters
   use InputData, only: Nstep
   use GNEB_driver
   use ErrorHandling

   !.. Implicit declarations
   implicit none

   ! Variables for timing and number of processors
   real(dblprec) :: time_a, time_b, time_c, time_d, time_e
   real(dblprec) :: time_I, time_II, time_III, time_IV, time_V
   integer :: nprocs

contains

   subroutine main()

      implicit none

   ! Executable statements
   !==============================================================!
   ! Check if inpsd.dat exists and whether it contains anything
   !--------------------------------------------------------------!
   call ErrorHandling_check_input_file()

   ! Find number of OpenMP processors active
   !--------------------------------------------------------------!
   nprocs = number_of_active_processors()

   ! Initialize profiling routines and print the UppASD logo
   !--------------------------------------------------------------!
   call print_logo()
   call memocc(0,0,'count','start')
   call timing(0,'serial        ','IN')
   call cpu_time(time_I)
   time_a=omp_get_wtime()

   ! Initialize arrays and read input
   !--------------------------------------------------------------!
   call timing(0,'Startup       ','ON')
   call setup_simulation()
   call timing(0,'Startup       ','OF')
   time_b=omp_get_wtime()
   call cpu_time(time_II)

   ! Initial phase
   !--------------------------------------------------------------!
   !call timing(0,'Initial       ','ON')
   call run_initial_phase()
   !call timing(0,'Initial       ','OF')
   time_c=omp_get_wtime()
   call cpu_time(time_III)

   ! Measurement Phase
   !--------------------------------------------------------------!
   call run_measurement_phase()
   time_d=omp_get_wtime()
   call cpu_time(time_IV)

   ! Deallocate left over data
   !--------------------------------------------------------------!
   call cleanup_simulation()

   ! End time analysis and stop memory analysis
   !--------------------------------------------------------------!
   call timing(0,'              ','RE')
   time_e=omp_get_wtime()
   call cpu_time(time_V)
   call sd_timing(Nstep)
   call memocc(0,0,'count','stop')

   end subroutine main

   !==============================================================!


   !---------------------------------------------------------------------------------
   !> @brief
   !> Find number of OpenMP processors active
   !
   !> @author
   !> Anders Bergman
   !
   !> @date 11/08/2014 - Thomas Nystrand
   !> - Moved to separate function
   !---------------------------------------------------------------------------------
   function number_of_active_processors() result(nprocs)
      use omp_lib
      integer :: nprocs

      !$omp parallel
      !$omp master
      nprocs=omp_get_num_threads()
      !write(*,'(1x,a18,i2,a16,i3,a10)') &
      !   "Using OpenMP with ",nprocs," threads out of",OMP_GET_NUM_PROCS(),"possible."
      !$omp end master
      !$omp end parallel
   end function number_of_active_processors


   !---------------------------------------------------------------------------------
   !> @brief
   !> Initialize simulation by setting up Hamiltonian and magnetic moments
   !> * Read input file
   !> * Print input data
   !> * Initialize random number generator
   !> * Set up Hamiltonian
   !> * Initialize magnetic moments
   !
   !> @author
   !> Anders Bergman
   !
   !> @date 11/08/2014 - Thomas Nystrand
   !> - Moved to separate function
   !---------------------------------------------------------------------------------
   subroutine run_initial_phase()
      use InputData, only : ipmode, iphfield
      use QMinimizer

      implicit none

      call timing(0,'Initial       ','ON')
      write (*,'(1x,a,2x,a)') "Enter initial phase:", ipmode

      if (ipmode=='M' .or. ipmode=='H'.or.ipmode=='I'.or.ipmode=='L'.or.ipmode=='W'.or.ipmode=='D') then
         ! Monte Carlo initial phase
         call mc_iphase()
      
      elseif (ipmode=='MS') then
        call ms_iphase()
        
      elseif (ipmode=='E') then
         ! Wang-Landau warmup
         call wl_iphase()
      elseif (ipmode=='B') then
         ! Spin-lattice Hamiltonian Monte Carlo initial phase
         call ErrorHandling_missing('Spin-lattice Monte Carlo')
      elseif (ipmode=='S') then
         ! Spin dynamics initial phase
         call sd_iphase()
      elseif (ipmode=='P') then
         write (*,'(1x,a)') "Calls ld_iphase"
         call ld_iphase() ! Lattice Dynamics initial phase
      elseif (ipmode=='R') then
         write (*,'(1x,a)') "Calls sld_iphase"
         call sld_iphase() ! Spin Lattice Dynamics initial phase
      elseif (ipmode=='X') then
         write (*,'(1x,a)') "Calls parallel tempering"
         call pt_iphase() ! Parallel tempering initial phase
      elseif (ipmode=='SX') then
         write (*,'(1x,a)') "Calls parallel tempering"
         ! 1Q init
         call qminimizer_wrapper('Q')
         call sx_copy_ens(1,iphfield)
         ! 3Q init
         call qminimizer_wrapper('Y')
         call sx_copy_ens(2,iphfield)
         ! Random init
         call sx_copy_ens(3,iphfield)
         ! FM init
         call sx_copy_ens(4,iphfield)
         call sx_iphase() ! Parallel tempering initial phase
      elseif (ipmode=='Q') then
         ! Spin spiral minimization initial phase
         call qminimizer_wrapper(ipmode)
      elseif (ipmode=='Z') then
         ! Spin spiral minimization initial phase
         call qminimizer_wrapper(ipmode)
      elseif (ipmode=='Y') then
         ! Spin spiral minimization initial phase
         call qminimizer_wrapper(ipmode)
      elseif (ipmode=='G') then
         !call ErrorHandling_missing('VPO minimization')
         call em_iphase()
      endif
      call timing(0,'Initial       ','OF')
   end subroutine run_initial_phase

   !---------------------------------------------------------------------------------
   !> @brief
   !> Run the main simulation loop and obtain measurements
   !> * Spin correlation
   !> * Magnetic averages
   !> * Spin transfer torque
   !
   !> @author
   !> Anders Bergman
   !
   !> @date 11/08/2014 - Thomas Nystrand
   !> - Moved to separate function
   !---------------------------------------------------------------------------------
   subroutine run_measurement_phase()

      use BLS,               only: calc_bls, deallocate_bls_data
      use Restart
      use InputData
      use macrocells
      use QMinimizer
      use MomentData
      use SystemData
      use LatticeData
      use Correlation
      use Correlation_print, only: print_gkw, print_gkt
      use ChemicalData
      use SimulationData
      use HamiltonianData
      use AutoCorrelation,   only: do_autocorr, allocate_autocorr, autocorr_init
      use OptimizationRoutines
      use DiaMag
      use ElkGeometry
      use MetaTypes
      
      integer :: cflag

      if(do_diamag=='Y') then
         call timing(0,'SpinCorr      ','ON')
         call setup_tensor_hamiltonian(NA,Natom,Mensemble,simid,emomM,mmom)
         call timing(0,'SpinCorr      ','OF')
      end if

      write (*,'(1x,a)') "Enter measurement phase:"

      ! Allocate and initialize spin arrays for autocorrelation
      call timing(0,'SpinCorr      ','ON')
      if(do_autocorr=='Y') call allocate_autocorr(Natom,Mensemble)
      if(do_autocorr=='Y') then 
         call autocorr_init(Natom,Mensemble,simid,rstep,initmag,emom)
      endif
      call timing(0,'SpinCorr      ','OF')

      if(mode=='M' .or. mode=='H'.or.mode=='I'.or.mode=='L'.or.mode=='W'.or.mode=='D') then
         call mc_mphase() ! Monte Carlo measurement phase
         
      elseif (mode=='MS') then
         call ms_mphase()

      elseif (mode=='B') then
         call print_siminfo()
         call sd_mphase_lite() ! Spin Dynamics measurement phase

      elseif (mode=='S') then
         call print_siminfo()
         if (gpu_mode==0) then !FORTRAN
            call sd_mphase() ! Spin Dynamics measurement phase
         else ! C++ or CUDA
            call sd_mphaseCUDA()
         endif

      elseif (mode=='E') then
         write (*,'(1x,a)') "Calls Wang-Landau measurements"
         call wl_mphase()

      elseif (mode=='B') then
         call ErrorHandling_missing('Spin-lattice Monte Carlo')

      elseif (mode=='P') then
         write (*,'(1x,a)') "Calls ld_mphase"
         call ld_mphase() ! Lattice Dynamics measurement phase

      elseif (mode=='R') then
         write (*,'(1x,a)') "Calls sld_mphase"
         if (SDEalgh .le. 20) then
            write (*,'(1x,a)') "velocity-Verlet combined with SD solver"
            call sld_mphase()
         elseif(SDEalgh==21 .or. SDEalgh==22) then
            write (*,'(1x,a)') "Fix point iteration Implicit midpoint solver"
            call sld_mphase3()
         end if

      elseif (mode=='Q') then
         ! Spin spiral minimization measurement phase
         call qminimizer_wrapper('S')

      elseif (mode=='G') then
         write (*,'(1x,a)') "GNEB measurement phase"
         call mep_mphase()
      endif

      ! Save moment information from final simulation step
      mstep = mstep-1
      call timing(0,'PrintRestart  ','ON')
      if (do_mom_legacy.ne.'Y') then
         call prn_mag_conf(Natom,mstep,Mensemble,'R',simid,mmom,emom,'',mode)
      else 
         call prnrestart(Natom,Mensemble,simid,mstep,emom,mmom)
      endif
      call timing(0,'PrintRestart  ','OF')

      if(do_prn_elk==1) then
         call prn_elk_geometry(Natom,Mensemble,simid,mstep,emom,mmom,NA,N1,N2,N3,C1,C2,& 
            C3,atype_inp,coord,jfile(1),maptype,posfiletype,Bas)
      end if

      call timing(0,'SpinCorr      ','ON')
      cflag = 2 


      ! Spin-correlation
      call correlation_wrapper(Natom,Mensemble,coord,simid,emomM,mstep,delta_t,NT_meta,  &
         atype_meta,Nchmax,achtype,sc,do_sc,do_sr,cflag)

      ! Lattice-correlation
      call correlation_wrapper(Natom,Mensemble,coord,simid,uvec,mstep,delta_t,NT,  &
         atype,Nchmax,achtype,uc,do_uc,do_ur,cflag)

      ! Velocity-correlation
      call correlation_wrapper(Natom,Mensemble,coord,simid,vvec,mstep,delta_t,NT,  &
         atype,Nchmax,achtype,vc,do_vc,do_vr,cflag)

      ! Angular momentum-correlation
      call correlation_wrapper(Natom,Mensemble,coord,simid,lvec,mstep,delta_t,NT,  &
         atype,Nchmax,achtype,lc,do_lc,do_lr,cflag)

      if (do_suc=='Y'.and.(do_sc=='Y'.or.do_sc=='Q').and.(do_uc=='Y'.or.do_uc=='Q')) then
         call print_gkw(NT, Nchmax, sc, uc, simid, 'su')
      end if

      if (do_suc=='Y'.and.(do_sc=='Y'.or.do_sc=='T').and.(do_uc=='Y'.or.do_uc=='T')) then
         call print_gkt(NT, Nchmax, sc, uc, simid, 'su')
      end if

      if (do_uvc=='Y'.and.(do_uc=='Y'.or.do_uc=='Q').and.(do_vc=='Y'.or.do_vc=='Q')) then
         call print_gkw(NT, Nchmax, uc, vc, simid, 'uv')
      end if

      if (do_uvc=='Y'.and.(do_uc=='Y'.or.do_uc=='T').and.(do_vc=='Y'.or.do_vc=='T')) then
         call print_gkt(NT, Nchmax, uc, vc, simid, 'uv')
      end if

      if (do_slc=='Y'.and.(do_sc=='Y'.or.do_sc=='T').and.(do_lc=='Y'.or.do_lc=='T')) then
         call print_gkt(NT, Nchmax, lc, sc, simid, 'sl')
      end if

      cflag = 3

      ! Spin-correlation deallocation
      call correlation_wrapper(Natom,Mensemble,coord,simid,emomM,mstep,delta_t,NT_meta,  &
         atype_meta,Nchmax,achtype,sc,do_sc,do_sr,cflag)

      ! Lattice-correlation deallocation
      call correlation_wrapper(Natom,Mensemble,coord,simid,uvec,mstep,delta_t,NT,  &
         atype,Nchmax,achtype,uc,do_uc,do_ur,cflag)

      ! Velocity-correlation deallocation
      call correlation_wrapper(Natom,Mensemble,coord,simid,vvec,mstep,delta_t,NT,  &
         atype,Nchmax,achtype,vc,do_vc,do_vr,cflag)

      ! Angular momentum-correlation deallocation
      call correlation_wrapper(Natom,Mensemble,coord,simid,lvec,mstep,delta_t,NT,  &
         atype,Nchmax,achtype,lc,do_lc,do_lr,cflag)

      !call lattcorrelation_wrapper(Natom,Mensemble,coord,simid,uvec,vvec,mstep,     &
      !   delta_t,NT,atype,Nchmax,achtype,cflag,cflag_p)

      ! Brillouin Light scattering function.
      cflag = 2 ; 
      call calc_bls(N1,N2,N3,C1,C2,C3,Natom,Mensemble,simid,coord,emomM,mstep,      &
      delta_t,cflag)

      call deallocate_bls_data()

      call timing(0,'SpinCorr      ','OF')

   end subroutine run_measurement_phase

   !---------------------------------------------------------------------------------
   !> @brief
   !> Finish simulation and deallocate leftovers
   !
   !> @author
   !> Anders Bergman
   !
   !> @date 11/08/2014 - Thomas Nystrand
   !> - Moved to separate function
   !---------------------------------------------------------------------------------
   subroutine cleanup_simulation()
      use LSF,          only : allocate_lsfdata
      use Tools
      use Energy,       only : allocate_energies
      use KMCData,      only : allocate_kmc,do_kmc,NA_KMC
      use Damping
      use clusters,     only : do_cluster, deallocate_cluster_info
      use InputData
      use macrocells,   only : Num_macro,allocate_macrocell,do_macro_cells
      use MomentData
      use SystemData
      use Correlation, only : do_sc, do_uc, do_vc, do_lc, sc, uc, vc, lc, do_sr
      use Correlation_utils, only : allocate_corr, corr_sr_init
      use SpinTorques
      use SpinIceData,  only : allocate_spinice_data, allocate_vertex_data
      use ChemicalData
      use Polarization
      use DipoleManager, only: dipole_cleanup
      use Qvectors, only : deallocate_q
      use MicroWaveField
      use HamiltonianData
      use Prn_Topology, only: skyno, do_proj_skyno
      use Gradients, only : deallocate_gradient_lists
      use LatticeInputData
      use OptimizationRoutines
      use LatticeHamiltonianData
      use MultiscaleInterpolation
      use MultiscaleSetupSystem
      use MultiscaleDampingBand

    if (do_multiscale) then
      call allocate_multiscale(flag=-1)
      call allocate_systemdata(flag=-1)
      call allocate_anisotropies(Natom,ham_inp%mult_axis,flag=-1)
      call allocate_mmoms(flag=-1)
      call allocate_emoms(Natom, 1, -1)
      call allocate_general(flag=-1)
      call allocate_damping(Natom,ipnphase,-1)
      call allocate_hamiltoniandata(Natom, 1, Natom,1,0, 0, 'N',-1, 'N','N')
      call deallocateDampingBand(dampingBand)
      call deleteLocalInterpolationInfo(interfaceInterpolation) 
   else

      write (*,'(1x,a)') "Simulation finished"
      call deallocate_q(do_sc) ! Deallocate spin correlation related arrays
      call allocate_mmoms(flag=-1)
      call deallocate_rest() ! Deallocate remaining arrays
      call allocate_initmag(flag=-1)
      call allocate_general(flag=-1)
      call allocate_damping(Natom,ipnphase,-1) ! Deallocate the damping
      call setup_mwf_fields(Natom,-1) ! Deallocate microwave fields
      call allocate_stt_data(Natom,Mensemble,flag=-1) ! Deallocate stt data
      call allocate_systemdata(flag=-1)
      call allocate_hamiltoniandata(do_jtensor=ham_inp%do_jtensor,do_lsf=do_lsf,flag=-1,    &
         lsf_field=lsf_field,exc_inter=ham_inp%exc_inter)
      ! Conditional deallocations
      if(ham_inp%do_dm==1) call allocate_dmhamiltoniandata(flag=-1)
      if(ham_inp%do_sa==1) call allocate_sahamiltoniandata(flag=-1)
      if(ham_inp%do_pd==1) call allocate_pdhamiltoniandata(flag=-1)
      if(ham_inp%do_bq==1) call allocate_bqhamiltoniandata(flag=-1)
      if(ham_inp%do_ring==1) call allocate_ringhamiltoniandata(flag=-1)      
      if(do_ll==1) call allocate_llhamiltoniandata(flag=-1)
      if(do_lll==1) call allocate_llhamiltoniandata(flag=-1)
      if(do_llll==1) call allocate_llhamiltoniandata(flag=-1)
      if(do_ml==1) call allocate_mlhamiltoniandata(flag=-1)
      if(do_mml>=1) call allocate_mmlhamiltoniandata(flag=-1)
      if(do_mmll==1) call allocate_mmllhamiltoniandata(flag=-1)
      if(do_ralloy==1.or.oldformat) call allocate_chemicaldata(flag=-1)
      if(do_ralloy==1.or.oldformat) call allocate_chemicalinput(flag=-1)
      if(ham_inp%do_anisotropy==1) call allocate_anisotropies(Natom,ham_inp%mult_axis,flag=-1)
      if(do_lsf=='Y') call allocate_lsfdata(Nchmax,conf_num,flag=-1)
      if(do_kmc.eq.'Y') call allocate_kmc(NA_KMC=NA_KMC,flag=-1)

      if(ham_inp%do_dip>0) call dipole_cleanup(ham_inp%do_dip,Num_macro)
      if(do_macro_cells=='Y'.or.ham_inp%do_dip==2) then
         call allocate_macrocell(-1,Natom,Mensemble)
      endif
      if(do_sc=='Y'.or.do_sc=='Q'.or.do_sc=='C') then
         call allocate_corr(Natom,Mensemble,sc,-1)
      end if
      if(do_uc=='Y'.or.do_uc=='Q'.or.do_uc=='C') then
         call allocate_corr(Natom,Mensemble,uc,-1)
      end if
      if(do_vc=='Y'.or.do_vc=='Q'.or.do_vc=='C') then
         call allocate_corr(Natom,Mensemble,vc,-1)
      end if
      if(do_lc=='Y'.or.do_lc=='Q'.or.do_lc=='C') then
         call allocate_corr(Natom,Mensemble,lc,-1)
      end if
      if(do_sr=='Y'.or.do_sr=='R'.or.do_sr=='T') then
         call corr_sr_init(Natom, sc, coord,-1)
      end if
      if(ind_mom_flag=='Y') call allocate_hamiltoniandata_ind(flag=-1)
      if(plotenergy>0) call allocate_energies(flag=-1)
      if(do_cluster=='Y') call deallocate_cluster_info()
      if (mode=='L') then
         call allocate_spinice_data(Natom,Mensemble,-1)
         call allocate_vertex_data(-1)
      end if
      ! --Optimization Region-- !
      if(OPT_ON) then
         call allocateOptimizationStuff(na,natom,Mensemble,.false.)
      end if
      if (stt/='N'.or.skyno=='Y'.or.do_proj_skyno=='Y') then
         call deallocate_gradient_lists()
      end if
    endif
   end subroutine cleanup_simulation

   !---------------------------------------------------------------------------------
   !> @brief
   !> Initialize simulation by setting up Hamiltonian and magnetic moments
   !
   !> @author
   !> Anders Bergman
   !---------------------------------------------------------------------------------
   subroutine setup_simulation()
      !
      !use Diamag_full
      use ams
      use KMC
      use BLS
      use LSF,             only : read_LSF,allocate_lsfdata
      use Sparse
      use Energy,          only: allocate_energies
      use Damping
      use KMCData
      use Damping
      use FixedMom
      use clusters
      use MC_Wolff
      use Topology
      use geometry,        only : setup_geometry, rescale_lattvec
      use gradients
      use Stiffness,       only : do_stiffness
      use InputData
      use SystemData
      use MomentData
      use fieldpulse
      use printinput
      use prn_fields
      use macrocells
      use latticeinit
      use Temperature
      use Correlation
      use Qvectors
      use SpinTorques
      use LatticeData
      use prn_topology
      use setupspinice
      use prn_currents
      use Polarization
      use ChemicalData
      use InputHandler
      use InputHandler_ext
      use RandomNumbers
      use clus_geometry
      use InducedMoments,  only : calculate_init_ind_sus!, renorm_ncoup_ind
      use Temperature_3TM
      use prn_microwaves
      use latticespectra
      use SimulationData
      use microwavefield
      use HamiltonianData
      use AutoCorrelation
      use hamiltonianinit
      use LatticeInputData
      use magnetizationinit
      use prn_latticefields
      use prn_micromagnetic
      use LatticeInputHandler
      use prn_latticeaverages
      use OptimizationRoutines
      use latticehamiltonianinit
      use LatticeHamiltonianData
      use prn_latticetrajectories
      use WangLandau
      use Diamag
      use ElkGeometry
      use Qminimizer
      use Correlation_core
      use Correlation_utils, only : allocate_corr, corr_sr_init
      use Omegas, only : set_w, set_w_print_option
      use MultiscaleInterpolation
      use MultiscaleSetupSystem,      only : setup_multiscale_system
      use MultiscaleDampingBand
      use prn_trajectories
      use prn_averages, only : read_parameters_averages,zero_cumulant_counters, avrg_init
      use MetaTypes
      
      implicit none

      integer :: i, i_stat,mconf

      ! Read inpsd.dat input file
      write(*,'(1x,a)',advance='no') 'Read input file'

      !Check which input format is used
      call check_format()

      write(*,'(1x,a)') 'using new format'
      ! Set defaults for KMC (ONLY FOR MAGNETIC POLARONS!)
      call init_kmc()
      ! Set defaults for Adiabatic Magnon Spectra calculations
      call init_ams()
      ! Set defaults for Spin Transfer Torques and related effects
      call init_stt()
      ! Set the defaults for the BLS analysis
      call setup_bls()
      ! Set the defaults for the trajectories
      call traj_init()
      ! Set the defaults for the averages, cumulants and autocorrelation
      call avrg_init()
      ! Set the default for the polarization measurements
      call prn_pol_init()
      ! Set defaults for current measurements
      call current_init()
      ! Set defaults for random number generation
      call rng_defaults()
      ! Set the defaults for the trajectories
      call latttraj_init()
      ! Set the defaults for the averages, cumulants and autocorrelation
      call lattavrg_init()
      ! Set defaults for vertex printing
      call spin_ice_init()
      ! Set defaults for fixed moments
      call init_fixed_mom()
      ! Set defaults for macrocells calculations
      call init_macrocell()
      ! Set defaults for ferromagnetic stiffness calculations
      call init_stiffness()
      ! Set defaults for the field printing
      call fields_prn_init()
      ! Set defaults for Optimization routines
      call set_opt_defaults(delta_t)
      ! Set input defaults for the correlation
      call correlation_init()
      ! Set input defaults for the correlation
      call correlation_init_loc(sc)
      call correlation_init_loc(uc)
      call correlation_init_loc(vc)
      call correlation_init_loc(lc)
      !call correlation_init_loc(suc)
      !call correlation_init_loc(uvc)
      ! Set input defaults for topological information
      call prn_topology_init()
      ! Set general input defaults
      call set_input_defaults()
      ! Set defaults for the lattice field printing
      call lattfields_prn_init()
      ! Set input defaults for the displacement correlation
      ! Set lattice dynamics input defaults
      call set_lattinput_defaults()
      ! Setup parameters for cluster embeding
      call set_input_defaults_clus()
      ! Set defaults for temperature routines
      call set_temperature_defaults()
      ! Set 3TM defaults
      call set_temperature_3tm_defaults()
      ! Set microwaves fields printing defaults
      call initialize_prn_microwaves()
      ! Set Wang-Landau defaults
      call wanglandau_init()
      ! Set q-sweep defaults
      call qminimizer_init()

      open(ifileno,file='inpsd.dat')
      call read_parameters(ifileno)
      rewind(ifileno)
      call read_parameters_sld(ifileno)
      rewind(ifileno)
      call read_parameters_bls(ifileno)
      rewind(ifileno)
      call read_parameters_mwf(ifileno)
      rewind(ifileno)
      call read_parameters_opt(ifileno)
      rewind(ifileno)
      call read_parameters_ams(ifileno)
      rewind(ifileno)
      call read_parameters_stt(ifileno)
      rewind(ifileno)
      call read_parameters_spinice(ifileno)
      rewind(ifileno)
      call read_parameters_cluster(ifileno,conf_num)
      rewind(ifileno)
      call read_parameters_averages(ifileno)
      rewind(ifileno)
      call read_parameters_trajectories(ifileno)
      rewind(ifileno)
      call read_parameters_autocorr(ifileno)
      rewind(ifileno)
      call read_parameters_correlation(ifileno)
      rewind(ifileno)
      call read_parameters_correlation_sc(ifileno)
      rewind(ifileno)
      call read_parameters_correlation_uc(ifileno)
      rewind(ifileno)
      call read_parameters_correlation_vc(ifileno)
      rewind(ifileno)
      call read_parameters_correlation_lc(ifileno)
      rewind(ifileno)
      call read_parameters_qvectors(ifileno)
      rewind(ifileno)
      call read_parameters_wanglandau(ifileno)
      rewind(ifileno)
      call read_parameters_diamag(ifileno)
      rewind(ifileno)
      call read_parameters_elkgeometry(ifileno)
      rewind(ifileno)
      call read_parameters_qminimizer(ifileno)
      rewind(ifileno)
      call read_parameters_3tm(ifileno)
      rewind(ifileno)
      call read_parameters_tempexp(ifileno)
      close(ifileno)
      
       
      !----------------------------!
      ! Print warning if warranted !
      !----------------------------!
      if(.not.sane_input) then
         call ErrorHandling_ERROR(&
            "Input data is inconsistent."//achar(10)// &
         " Please ensure consistency between simulation method and input parameters.")
      end if

      !----------------------------!
      ! Reading in data from files !
      !----------------------------!
    if (do_multiscale) then
     call initDampingBand(dampingBand)
     call newLocalInterpolationInfo(interfaceInterpolation)
     call setup_multiscale_system(multiscale_file_name)
     call setup_mwf_fields(Natom,1)
     if(SDEalgh>=1) &
     !call allocate_randomwork(Natom,Mensemble,1,'Y') 
      ipSDEalgh=SDEalgh
     ! Initialize random number generator
      !$omp parallel copyin(tseed)
      if(para_rng) tseed=tseed+omp_get_thread_num()
      call rng_init(tseed)
      !$omp end parallel 
#ifdef VSL
      if(.not.use_vsl) call setup_rng_hb(tseed,ziggurat,rngpol)
#else
      call setup_rng_hb(tseed,ziggurat,rngpol)
#endif
    else

      ! Read atom positions for unit cell or full input
      if(aunits=='Y') then
         write(*,'(1x,a)') 'Use atomic units'
         call change_constants()
      end if

      ! Rescale lattice vectors
      call rescale_lattvec(scalefac,C1,C2,C3)

      ! Read geometry and spin configuration on ELK format
      if(do_read_elk==1) then
         call read_elk_geometry()
      end if

      ! If the cluster embedding method is used
      if (do_cluster=='Y') then
         call reading_wrapper_clus(conf_num,do_ralloy,set_landeg,do_anisotropy_clus,&
            Landeg_glob,do_lsf,posfiletype,ind_mom_flag)
      endif

      if(.not.allocated(Bas) .and. do_read_elk==0) then
         if(do_ralloy==0) then
            call read_positions()
         else
            call read_positions_alloy()
         end if
      end if

      Natom_full=NA*N1*N2*N3
      !------------------------------------------------------------------------------
      ! Setting up the cluster geometry
      !------------------------------------------------------------------------------
      if (do_cluster=='Y') then
         write(*,'(1x,a)',advance='no') 'Setup cluster geometry'
         call setup_clus_geometry(NA,N1,N2,N3,N1_clus,N2_clus,N3_clus,NA_clus,      &
            do_ralloy,block_size,Natom_full,Nchmax_clus,Natom_full_clus,Nch_clus,   &
                  C1,C2,C3,C1_clus,C2_clus,C3_clus,Bas,chconc_clus,clus_expand,     &
            index_clus,Natom_clus,atype_inp_clus,anumb_inp_clus,Bas_clus,coord_clus)
         write(*,'(1x,a)') 'done'
      endif
      !------------------------------------------------------------------------------
      ! End of cluster geometry setup
      !------------------------------------------------------------------------------
      if (do_ralloy/=0) call allocate_chemicaldata(Natom_full,1)

      if(do_read_elk==0) then
         if (do_fixed_mom.eq.'Y') then
            call read_fixed_moments(Landeg_glob)
         else
            call read_moments(Landeg_glob)
         endif
      endif

      !Read ionic displacements for unit cell or full input
      if(do_ld == 'Y') then
         write(*,'(1x,a)') 'Calls read_phonons'
         call read_phonons()
      end if

      ! Read LSF related items
      if (do_lsf == 'Y') then
         mconf=gsconf_num
         call allocate_lsfdata(Nchmax,conf_num,1)
         call read_LSF(lsffile,nchmax,conf_num,mconf)
      else
         mconf=1
      endif

      ! Site dependent damping for the initial phase
      if (ipmode=='S'.and.do_site_ip_damping=='Y') then
         if (do_ralloy==0) then
            call read_ip_damping()
         else
            call read_ip_damping_alloy()
         endif
      end if
      ! Site dependent dampinf for the measurement phase
      if (mode=='S'.and.do_site_damping=='Y') then
         if(do_ralloy==0) then
            call read_damping()
         else
            call read_damping_alloy()
         end if
      endif


      ! Set initial phase solver to same by default
      if(ipSDEalgh==-1) ipSDEalgh=SDEalgh

      !----------------------------!
      ! Initialize and allocations !
      !----------------------------!

      ! Calculate the combined tensor after all reads
      allocate(atype(Natom_full),stat=i_stat)
      call memocc(i_stat,product(shape(atype))*kind(atype),'atype','read_initphase')
      allocate(anumb(Natom_full),stat=i_stat)
      call memocc(i_stat,product(shape(anumb))*kind(anumb),'anumb','read_initphase')
      allocate(Landeg(Natom_full),stat=i_stat)
      call memocc(i_stat,product(shape(Landeg))*kind(Landeg),'Landeg','read_initphase')
      call zero_cumulant_counters()

      write(*,'(a)') ' done.'

      ! Initialize random number generator
      !$omp parallel copyin(tseed)
      if(para_rng) tseed=tseed+omp_get_thread_num()
      call rng_init(tseed)
      !$omp end parallel
#ifdef VSL
      if(.not.use_vsl) call setup_rng_hb(tseed,ziggurat,rngpol)
#else
      call setup_rng_hb(tseed,ziggurat,rngpol)
#endif

      if (do_cluster=='Y') then
         ! Set up geometry containing chemical data with the cluster information
         call modify_global_coordinates(Natom,NA,N1,N2,N3,Bas,C1,C2,C3,coord,    &
            atype,anumb,atype_inp,anumb_inp,do_prnstruct,do_prn_poscar,tseed,simid, &
            do_ralloy,Natom_full,Natom_clus,Nchmax,Nch,acellnumb,acellnumbrev,      &
            achtype,chconc,atype_ch,asite_ch,achem_ch,block_size)
      else
         ! Set up geometry containing chemical data
         call setup_geometry(Natom,NA,N1,N2,N3,Bas,C1,C2,C3,coord,atype,anumb,      &
            atype_inp,anumb_inp,do_prnstruct,do_prn_poscar,tseed,simid,do_ralloy,   &
            Natom_full,Nchmax,Nch,acellnumb,acellnumbrev,achtype,chconc,atype_ch,   &
            asite_ch,achem_ch,block_size)
      endif

      ! Read exchange interaction from file, either normal or tensor format
      if(ham_inp%do_jtensor==0) then
         ! Normal exchange interaction read
         call read_exchange(ham_inp)
      elseif(ham_inp%do_jtensor==1 .and. ham_inp%calc_jtensor.eqv..false.) then
         ! Read tensor exchange format file for tensor coupling
         call read_exchange_tensor()
      elseif(ham_inp%do_jtensor==1 .and. ham_inp%calc_jtensor.eqv..true.) then
         ! Read normal exchange file which is thereafter used to build tensor coupling
         call read_exchange_build_tensor()
      else
         call ErrorHandling_ERROR("Unrecognized or unsupported combination of do_tensor &
            &and calc_tensor")
      end if
      ham_inp%max_no_shells=maxval(ham_inp%NN)

      ! Reading the positions for the KMC calculation
      if (do_kmc.eq.'Y') then
         if(.not.allocated(Bas_KMC)) then
            if(do_ralloy==0) then
               call read_positions_KMC(C1,C2,C3,posfiletype)
            else
               call read_positions_alloy_KMC(C1,C2,C3,posfiletype)
            end if
         end if
      endif

      if(ham_inp%do_dm==1)         call read_dmdata()
      if(ham_inp%do_sa==1)         call read_sadata()
      if(ham_inp%do_pd==1)         call read_pddata()
      if(ham_inp%do_chir==1)       call read_chirdata()
      if(ham_inp%do_bq==1)         call read_bqdata()
      if(ham_inp%do_ring==1)       call read_ringdata()      
      if(ham_inp%do_biqdm==1)      call read_biqdmdata()
      if(do_ll==1)         call read_lldata()
      if(do_ml==1)         call read_mldata()
      if(do_mml>=1)        call read_mmldata()
      if(do_lll==1)        call read_llldata()
      if(do_llll==1)       call read_lllldata()
      if(do_mmll==1)       call read_mmlldata()
      if(do_ll_phonopy==1) call read_llphonopydata()
      if(do_kmc=='Y')      call read_barriers()
      if(do_autocorr=='Y') call read_tw(twfile)
      if(do_bpulse>=1.and.do_bpulse<5) call read_bpulse(do_bpulse,bpulsefile)

      if(ham_inp%do_anisotropy==1) then
         if(do_ralloy==0) then
            call read_anisotropy()
         else
            call read_anisotropy_alloy()
         end if
      end if
      ! Allocate energy arrays
      if (plotenergy>0) call allocate_energies(1,Mensemble)
      ! Allocating the damping arrays
      call allocate_damping(Natom,ipnphase,0)

      write(*,'(1x,a)') ' Set up damping'
      call setup_damping(NA,Natom,Natom_full,ipmode,mode,ipnphase,do_ralloy,        &
         asite_ch,achem_ch,do_site_damping,do_site_ip_damping,iplambda1,iplambda2,  &
         mplambda1,mplambda2,iplambda1_array,iplambda2_array,lambda1_array,         &
         lambda2_array)

      write(*,'(1x,a)') ' done'

      if (do_macro_cells.eq.'Y'.or.ham_inp%do_dip.eq.2) then
         write(*,'(1x,a)',advance='no') ' Create macro cells '
         call create_macrocell(NA,N1,N2,N3,Natom,Mensemble,block_size,coord,        &
            Num_macro,max_num_atom_macro_cell,cell_index,macro_nlistsize,           &
            macro_atom_nlist,simid)
         write(*,'(a)') 'done'
      endif

      ! See if it is necesary to read the temperature file
      !if(grad.eq.'Y') call read_temperature_legacy()
      !if(grad.eq.'Y'.or.grad.eq.'F') call read_temperature()

      ! Allocating the temperature arrays
      call allocate_Temp(Natom,ipnphase)

      write(*,'(1x,a)') ' Set up temperature'
      if (ipnphase.ge.1) then
         do i=1, ipnphase
            if(grad=='Y') then
               !call SETUP_TEMP(NATOM,NT,NA,N1,N2,N3,NATOM_FULL,&
               !   DO_RALLOY,ATYPE,ACELLNUMB,ATYPE_CH,SIMID,iptemp(i),&
               !   C1,C2,C3,BC1,BC2,BC3,BAS,COORD,ipTemp_array(:,i))
               call Lparray(ipTemp_array(:,i),Natom,coord,iptemp(i),simid,.false.)
            else if(grad=='F') then
               call SETUP_TEMP(NATOM,NT,NA,N1,N2,N3,NATOM_FULL,DO_RALLOY,ATYPE,     &
                  ACELLNUMB,ATYPE_CH,SIMID,iptemp(i),C1,C2,C3,BC1,BC2,BC3,BAS,COORD,&
                  ipTemp_array(:,i))
            else
               ipTemp_array(:,i)=ipTemp(i)
            end if

         enddo
      endif
      if(grad=='Y') then
         !call SETUP_TEMP(NATOM,NT,NA,N1,N2,N3,NATOM_FULL,&
         !   DO_RALLOY,ATYPE,ACELLNUMB,ATYPE_CH,SIMID,TEMP,&
         !   C1,C2,C3,BC1,BC2,BC3,BAS,COORD,Temp_array)
         call Lparray(Temp_array,Natom,coord,Temp,simid,.true.)
      else if(grad=='F') then
         call SETUP_TEMP(NATOM,NT,NA,N1,N2,N3,NATOM_FULL,DO_RALLOY,ATYPE,ACELLNUMB, &
            ATYPE_CH,SIMID,TEMP,C1,C2,C3,BC1,BC2,BC3,BAS,COORD,Temp_array)
      else
         Temp_array=TEMP
      end if

      write(*,'(1x,a)') ' done'

      ! Spin-Ice setup
      if (mode=='L') then
         call read_vertices()
         write(*,'(a)') 'Setup vertex geometry'
         call setup_vertex_geometry(N1,N2,N3,do_prnstruct,simid)
         write(*,'(a)') 'Setup Spin Ice'
         call setup_ice_neighbours(Natom,Mensemble,NT,NA,N1,N2,N3,BC1,BC2,BC3,atype,&
            Bas,sym,simid,coord)
      endif

      !
      ! 
      i=max(nstep,mcnstep)
      ! Set up q-point grid if real space correlation function is to be calculated
         if (qpoints == 'C' .or. qpoints == 'P' .or. qpoints == 'A'.or.qpoints=='X'.or.qpoints=='H') then
            write (*,'(1x,a)') "Set up q-vectors"
            call setup_qcoord(N1, N2, N3, C1, C2, C3)
         else if(qpoints == 'F' .or. qpoints == 'D' .or. qpoints=='I' .or. qpoints=='B' .or. qpoints == 'G') then
            write (*,'(1x,a)') "Read q-vectors"
            call read_q(C1, C2, C3)
         end if
         call set_w_print_option(real_time_measure)
         if(do_sc=='Q'.or.do_sc=='Y') then
            write (*,'(1x,a)') "Set up spin-correlation frequencies"
            call set_w(delta_t,sc)
         end if
         if(do_uc=='Q'.or.do_uc=='Y') then
            write (*,'(1x,a)') "Set up lattice correlation frequencies"
            call set_w(delta_t,uc)
         end if
         if(do_vc=='Q'.or.do_vc=='Y') then
            write (*,'(1x,a)') "Set up velocity correlation frequencies"
            call set_w(delta_t,vc)
         end if
         if(do_lc=='Q'.or.do_lc=='Y') then
            write (*,'(1x,a)') "Set up angular momentum correlation frequencies"
            call set_w(delta_t,lc)
         end if

      !
      call set_correlation_average(sc,i)
      call set_correlation_average(uc,i)
      call set_correlation_average(vc,i)
      call set_correlation_average(lc,i)
      !

         if(do_sc=='Y'.or.do_sc=='Q'.or.do_sc=='C') then
            write (*,'(1x,a)') "Allocate spin-correlation data"
            call allocate_corr(Natom,Mensemble,sc,1)
         end if
         if(do_uc=='Y'.or.do_uc=='Q'.or.do_uc=='C') then
            write (*,'(1x,a)') "Allocate angular momentum correlation data"
            call allocate_corr(Natom,Mensemble,uc,1)
         end if
         if(do_vc=='Y'.or.do_vc=='Q'.or.do_vc=='C') then
            write (*,'(1x,a)') "Allocate angular momentum correlation data"
            call allocate_corr(Natom,Mensemble,vc,1)
         end if
         if(do_lc=='Y'.or.do_lc=='Q'.or.do_lc=='C') then
            write (*,'(1x,a)') "Allocate  angular momentum correlation data"
            call allocate_corr(Natom,Mensemble,lc,1)
         end if

         if(do_sr=='Y'.or.do_sr=='R'.or.do_sr=='T') then
            write(*,'(1x,a)',advance='no') "Set up neighbourmap for correlations"
            call corr_sr_init(Natom, sc, coord, 1)
            write(*,'(a)')" done."
         end if

      ! Print input
      write(*,'(1x,a)',advance='no') 'Write input data'
      call prninp()
      call print_yaml()
      write(*,'(a)') ' done.'
      write (*,'(1x,a)') "Set up Hamiltonian"

      ! Set up Hamiltonian, containing exchange, anisotropy, and optional terms
      ! like DM and dipolar interactions.
      !!! call setup_hamiltonian(NT,NA,N1,N2,N3,Nchmax,do_ralloy,Natom_full,Mensemble,  &
      !!!    nHam,Natom,achtype,atype_ch,asite_ch,achem_ch,acellnumb,acellnumbrev,atype,&
      !!!    anumb,alat,C1,C2,C3,Bas,ammom_inp,coord,BC1,BC2,BC3,sym,do_jtensor,        &
      !!!    max_no_shells,nn,redcoord,jc,jcD,jc_tens,do_dm,max_no_dmshells,dm_nn,      &
      !!!    dm_redcoord,dm_inpvect,do_anisotropy,random_anisotropy_density,            &
      !!!    anisotropytype,anisotropytype_diff,anisotropy,anisotropy_diff,             &
      !!!    random_anisotropy,mult_axis,mconf,conf_num,map_multiple,do_lsf,lsf_field,  &
      !!!    exc_inter,do_bq,max_no_bqshells,bq_nn,bq_redcoord,jc_bq,do_ring,           &
      !!!    max_no_ringshells,ring_nn,ring_redcoord_ij,ring_redcoord_ik,               &
      !!!    ring_redcoord_il,jc_ring,do_biqdm,                                         &
      !!!    max_no_biqdmshells,biqdm_nn,biqdm_redcoord,biqdm_inpvect,do_pd,            &
      !!!    max_no_pdshells,pd_nn,pd_redcoord,pd_inpvect,do_dip,                       &
      !!!    do_chir,max_no_chirshells,chir_nn,chir_redcoord,chir_inpval ,              &
      !!!    do_fourx,max_no_fourxshells,fourx_nn,fourx_redcoord,fourx_inpval ,         &
      !!!    ind_mom,ind_tol,ind_mom_flag,NA_clus,NT_clus,N1_clus,N2_clus,N3_clus,      &
      !!!    Natom_clus,Nchmax_clus,clus_expand,Natom_full_clus,max_no_shells_clus,     &
      !!!    max_no_dmshells_clus,NN_clus,dm_nn_clus,achtype_clus,atype_ch_clus,        &
      !!!    asite_ch_clus,achem_ch_clus,acellnumb_clus,acellnumbrev_clus,atype_clus,   &
      !!!    anumb_clus,atype_inp_clus,anumb_inp_clus,C1_clus,C2_clus,C3_clus,Bas_clus, &
      !!!    ammom_inp_clus,redcoord_clus,jc_clus,dm_redcoord_clus,dm_inpvect_clus,     &
      !!!    do_cluster,do_anisotropy_clus,random_anisotropy_density_clus,              &
      !!!    anisotropytype_clus,anisotropytype_diff_clus,anisotropy_clus,              &
      !!!    anisotropy_diff_clus,random_anisotropy_clus,mult_axis_clus,Num_macro,      &
      !!!    max_num_atom_macro_cell,macro_nlistsize,macro_atom_nlist,block_size,       &
      !!!    do_reduced,do_prnstruct,do_sortcoup,simid,print_dip_tensor,read_dipole,    &
      !!!    qdip_files)
      call setup_hamiltonian(NT,NA,N1,N2,N3,Nchmax,do_ralloy,Natom_full,Mensemble,  &
         nHam,Natom,achtype,atype_ch,asite_ch,achem_ch,acellnumb,acellnumbrev,atype,&
         anumb,alat,C1,C2,C3,Bas,ammom_inp,coord,BC1,BC2,BC3,sym,        &
         mconf,conf_num,do_lsf,lsf_field,  &
         ind_mom,ind_tol,ind_mom_flag,NA_clus,NT_clus,N1_clus,N2_clus,N3_clus,      &
         Natom_clus,Nchmax_clus,clus_expand,Natom_full_clus,max_no_shells_clus,     &
         max_no_dmshells_clus,NN_clus,dm_nn_clus,achtype_clus,atype_ch_clus,        &
         asite_ch_clus,achem_ch_clus,acellnumb_clus,acellnumbrev_clus,atype_clus,   &
         anumb_clus,atype_inp_clus,anumb_inp_clus,C1_clus,C2_clus,C3_clus,Bas_clus, &
         ammom_inp_clus,redcoord_clus,jc_clus,dm_redcoord_clus,dm_inpvect_clus,     &
         do_cluster,do_anisotropy_clus,random_anisotropy_density_clus,              &
         anisotropytype_clus,anisotropytype_diff_clus,anisotropy_clus,              &
         anisotropy_diff_clus,random_anisotropy_clus,mult_axis_clus,Num_macro,      &
         max_num_atom_macro_cell,macro_nlistsize,macro_atom_nlist,block_size,       &
         do_reduced,do_prnstruct,do_sortcoup,simid,ham_inp%print_dip_tensor,ham_inp%read_dipole,    &
         ham_inp%qdip_files)
      !call setup_reduced_hamiltonian(Natom,NA,conf_num)
      ! Allocate arrays for simulation and measurement
      call allocate_general(1)

      if(do_sparse=='Y') then
         if(ham_inp%do_jtensor==1) then
            write (*,'(1x,a)') "Setting up sparse Hamiltonian from tensor form"
            call setupSparseTensor(Natom,nHam,conf_num,ham%max_no_neigh,ham%nlist,   &
               ham%nlistsize,ham%j_tens, ham%aham)
         else if (ham_inp%do_dm==1) then
            write (*,'(1x,a)') "Setting up sparse Hamiltonian with DMI"
            call setupSparseBlock(Natom,nHam,conf_num,ham%max_no_neigh,ham%nlist,   &
               ham%nlistsize,ham%ncoup,ham%max_no_dmneigh,ham%dmlistsize,ham%dmlist,&
               ham%dm_vect,ham%aham)
         else
            write (*,'(1x,a)') "Setting up sparse Hamiltonian"
            call setupSparseScalar(Natom,nHam,conf_num,ham%max_no_neigh,ham%nlist,  &
               ham%nlistsize,ham%ncoup,ham%aham)
         end if
      end if

      write (*,'(1x,a)') "Initialize metatype "
      call find_metatype_coordination(Natom,atype,ham,metatype)
      call find_metanumb(Natom,NA,N1,N2,N3,Bas,C1,C2,C3,BC1,BC2,BC3,&
      coord,atype,anumb,do_ralloy,metanumb)

      write (*,'(1x,a)') "Initialize magnetic moments "

      ! Allocate arrays for magnetic moment magnitudes
      call allocate_mmoms(Natom,Mensemble,1)

      ! Fill arrays with magnetic moment magnitudes from input
      call setup_moment(Natom,Mensemble,NA,N1,N2,N3,ammom_inp,      &
         Landeg_ch,Landeg,mmom,mmom0,mmomi,do_ralloy,achtype,acellnumb,  &
         mconf)

      ! Set up magnetic moment vectors
      call magninit(Natom,Mensemble,NA,N1,N2,N3,initmag,Nchmax,aemom_inp,  &
         anumb,do_ralloy,Natom_full,achtype,acellnumb,emom,emom2,emomM,mmom,rstep,  &
         theta0,phi0,mavg0,restartfile,initrotang,initpropvec,initrotvec,coord,C1,C2,C3,  &
         do_fixed_mom,Nred,red_atom_list,ham%ind_list_full,ind_mom_flag,ind_mom,    &
         ham%fix_num,ham%fix_list,read_ovf,do_mom_legacy,relaxed_if)

      ! Treat the embedding of the cluster
      if(do_cluster=='Y') then
         if (Natom_clus.ge.Natom) write(*,'(a)') "WARNING: Number of atoms in the cluster is larger than in the host! Check input"
         write(*,'(1x,a)') "Embed cluster infromation into the host"
         call cluster_creation(NT,ham_inp%do_dm,Natom,initmag,conf_num,Mensemble,    &
            do_ralloy,Natom_full,ham_inp%do_jtensor,do_prnstruct,do_anisotropy_clus,        &
            index_clus,atype_clus,anumb_clus,coord,coord_clus,simid,mult_axis_clus, &
            atype,anumb,asite_ch,achem_ch,mmom,mmomi,emom,     &
            emomM,Landeg)
         call allocate_clusdata(flag=-1)
         if (do_ralloy/=0) call allocate_chemicaldata_clus(flag=-1)
         write(*,'(1x,a)')"Embedding done"
      endif

      ! Calculate the macrospin magnetic moments per macrocell if the dipolar interaction is considered
      ! with the macro spin model
      if (ham_inp%do_dip.eq.2) then
         call calc_macro_mom(Natom,Num_macro,Mensemble,max_num_atom_macro_cell,     &
            macro_nlistsize,macro_atom_nlist,mmom,emom,emomM,mmom_macro,emom_macro, &
            emomM_macro)
         if (prn_dip_subset.eq.'Y') then
            call read_dipole_subset()
            call calculate_dipole_subset(Natom,Mensemble,Num_macro,                 &
               max_num_atom_macro_cell,macro_nlistsize,macro_atom_nlist,emomM_macro,&
               ham%Qdip_macro,simid)
         endif
      endif

      ! If induced moments are considered, calculate the weight factor this needs to be calculated only once
      if (ind_mom_flag=='Y') then
         write(*,'(1x,a)',advance='no') "Setup induced moments information"
         call calculate_init_ind_sus(Natom,Mensemble,                &
            ham%max_no_neigh_ind, &
            ham%ind_list_full,ham%ind_nlistsize,&
            ham%ind_nlist,mmom,emom,emomM,ham%sus_ind)
         write(*,'(a)')" done"
      endif

      ! Rotation of the initial spin configuration
      if (roteul >= 1) then
         write(*,'(1x,a)',advance='no') 'Perform Euler rotation'
         call rotationeuler(Natom,Mensemble,roteul,rotang,emom,emomM,mmom,emom2)
         write (*,'(a)') ' done.'
      end if

      ! Excitation of the initial spin configuration
      if (initexc == 'I') then
         write(*,'(1x,a,f10.4,a)',advance='no') 'Introduces ', initconc*100, '% vacancies'
         call setinitexc(Natom,Mensemble,emom,emomM,mmom,mmom2,emom2,initexc,       &
            initconc,initneigh,initimp,ham%max_no_neigh,ham%nlist, mseed)
         write (*,'(a)') ' done.'
      else if (initexc == 'R') then
         write(*,'(1x,a,f10.4,a)',advance='no') 'Introduces ', initconc*100, '% two-magnon Raman spin flips'
         call setinitexc(Natom,Mensemble,emom,emomM,mmom,mmom2,emom2,initexc,       &
            initconc,initneigh,initimp,ham%max_no_neigh,ham%nlist, mseed)
         write (*,'(a)') ' done.'
      end if

      if (stt/='N'.or.do_she/='N'.or.do_sot/='N') then

         call read_jvecfile(Natom)

         if (stt=='A'.or.(stt/='A'.and.skyno=='Y')) then
            call setup_stencil_mesh(Natom,N1,N2,N3,C1,C2,C3,BC1,BC2,BC3,            &
               ham%max_no_neigh,ham%nlistsize,ham%nlist,coord)
         end if
         ! Call to allocate the needed stt data
         call allocate_stt_data(Natom,Mensemble,flag=1)
         if (stt/='N' .or. do_she/='N') then
            ! Call the printing of the current density in proper units
            call print_curr_density(NA,Nchmax,conf_num,alat,spin_pol,C1,C2,C3,jvec, &
               ammom_inp)
         endif

         if (do_sot/='N') call read_sot_pol_site(Natom)

      endif
      if (stt=='N'.and.(skyno=='Y'.or.do_proj_skyno=='Y')) then
         call setup_stencil_mesh(Natom,N1,N2,N3,C1,C2,C3,BC1,BC2,BC3,               &
            ham%max_no_neigh,ham%nlistsize,ham%nlist,coord)
      end if

      if (skyno=='T') then
         write(*,'(1x, a)') "Triangulating mesh"
         call delaunay_tri_tri(n1,n2,n3, NA)
      end if

      if (mode=='W') then
         call create_wolff_neigh_list(Natom,N1,N2,N3,C1,C2,C3,BC1,BC2,BC3,          &
            ham%max_no_neigh,ham%nlistsize,ham%nlist,coord,wolff_list_neigh_size,   &
            wolff_list_neigh)
      endif

      call setup_mwf_fields(Natom,1)

      if(do_bpulse==6) then
         call read_sitefield(Natom,sitenatomfld)
      endif

      ! --Optimization Region-- !
      if(OPT_flag) then
         call allocateOptimizationStuff(na,natom,Mensemble,.true.)
         call getCellPosNumNeigh(Natom, NA, ham%nlistsize, cellPosNumNeigh)
         OPT_ON=.true.
      end if

      ! Check AMS-flag
      if (do_ams =='Y'.and.ham_inp%do_jtensor.ne.1) then
         if(do_ralloy==0) call calculate_ams(N1,N2,N3,NA_meta,NT,Natom,anumb_meta,NA_metalist,simid,hfield,do_ralloy)
         if(do_ralloy>0)  call calculate_random_ams()
      end if

      ! Check if the micromagnetic information is calculated
      if (do_stiffness =='Y') then
         call stiffness_wrapper(NT,NA,N1,N2,N3,1,mconf,Natom,nHam,Nchmax,conf_num,  &
            do_ralloy,Natom_full,ham%max_no_neigh,Nch,anumb,atype,ham%aham,         &
            ham%nlistsize,atype_ch,asite_ch,achem_ch,ham%nlist,alat,C1,C2,C3,coord, &
            chconc,ammom_inp,ham%ncoup,ham%max_no_dmneigh,ham%dmlistsize,ham%dmlist,&
            ham%dm_vect,ham_inp%do_anisotropy,ham_inp%anisotropy,simid)
      end if

      ! Deallocate input data for Heisenberg Hamiltonian
      call allocate_hamiltonianinput(ham_inp,flag=-1)

      !------------------------------------------------------------------------------
      ! This is the initialization of the KMC particles
      !------------------------------------------------------------------------------
      ! Check if KMC is performed
      if (do_kmc=='Y') then
         write(*,'(1x,a)') "Setup KMC"
         call initialize_KMC(NT,NA,N1,N2,N3,sym,Natom,Nchmax,NA_KMC,do_ralloy,      &
            Mensemble,Natom_full,ham%max_no_neigh,do_prnstruct,                     &
            max_no_shells_barriers,nn_barriers,atype,anumb,ham%nlistsize,atype_ch,  &
            acellnumb,ham%nlist,C1,C2,C3,Bas,coord,ham%ncoup,redcoord_barriers,     &
            kmc_barriers,BC1,BC2,BC3,do_sortcoup,simid,max_no_neigh_barriers,       &
            asite_ch,achem_ch,kmc_index_prv,kmc_index_aft,nlistsize_barriers,       &
            nlist_barriers,ncoup_barriers,mmom,emom,emomM)
         write(*,'(1x,a)') "done"
      endif
      !------------------------------------------------------------------------------
      ! Ending of the KMC initialization
      !------------------------------------------------------------------------------

      if (mode=='P' .or. mode=='R') then

         write (*,'(1x,a)') "Set up lattice Hamiltonian"
         ! Set up lattice Hamiltonian, containing harmonic force constants
         ! higher order couplings, and spin-lattice coupling couplings
         call setup_latticehamiltonian(simid,Natom,NT,NA,N1,N2,N3,C1,C2,C3,BC1,BC2, &
            BC3,atype,anumb,coord,alat,Bas,do_prnstruct,do_sortcoup,sym,do_n3,do_ll,&
            ll_nn,nn_ll_tot,max_no_llshells,ll_listsize,ll_list,ll_tens,ll_redcoord,&
            ll_inptens,do_lll,lll_nn,nn_lll_tot,max_no_lllshells,lll_listsize,      &
            lll_list,lll_tens,lll_redcoord,lll_inptens,do_llll,llll_nn,nn_llll_tot, &
            max_no_llllshells,llll_listsize,llll_list,llll_tens,llll_redcoord,      &
            llll_inptens,do_ml,ml_nn,nn_ml_tot,max_no_mlshells,ml_listsize,ml_list, &
            ml_tens,ml_redcoord,ml_inptens,do_mml,mml_nn,nn_mml_tot,                &
            max_no_mmlshells,mml_listsize,mml_list,mml_tens,mml_redcoord,           &
            mml_inptens,mml_invsym,do_mmll,mmll_nn,nn_mmll_tot,max_no_mmllshells,   &
            mmll_listsize,mmll_list,mmll_tens,mmll_redcoord,mmll_inptens,do_ralloy, &
            acellnumb,acellnumbrev,achtype,nchmax,natom_full,ammom_inp,atype_ch,    &
            asite_ch,achem_ch,lm_listsize,lm_list,lm_tens,lmm_listsize,lmm_list,    &
            lmm_tens,llmm_listsize,llmm_list,llmm_tens,do_ll_phonopy,Natom_phonopy, &
            atomindex_phonopy,ll_inptens_phonopy)

         write (*,'(1x,a)') "Initialize ionic displacements and velocities"

         ! Allocate arrays for ionic displacements and velocities
         call allocate_latticedata(Natom,Mensemble, 1)

         ! Set up ionic displacements and velocities
         call lattinit(Natom,Mensemble,NA,N1,N2,N3,initlatt,Nchmax,mion_inp,        &
            uvec_inp,vvec_inp,anumb,do_ralloy,Natom_full,achtype,acellnumb,mion,    &
            mioninv,uvec,vvec,rstep,lattrestartfile)

         ! Rotation of the initial ionic displacement and velocity configuration
         if (lattroteul >= 1) then
            write(*,'(1x,a)',advance='no') 'Perform Euler rotation of ionic displacements and velocities'
            call lattrotationeuler(Natom, Mensemble, lattroteul, lattrotang, uvec, vvec)
            write (*,'(a)') ' done.'
         end if

         ! Calculate the center of mass
         call calc_masscenter(Natom, Mensemble, coord, mion, uvec)

         ! Check phonon spectra calculation flag
         if (do_phonspec =='Y') then
            call calculate_phondisp(simid,Natom,Mensemble,NT,NA,N1,N2,N3,C1,C2,C3,  &
               BC1,BC2,BC3,atype,anumb,coord,mioninv,Bas,nq,q,do_phondos,           &
               phondosfile,phondos_sigma,phondos_freq)
         end if
      end if
    endif

    !call get_rlist_corr(Natom,NA,coord,1)
   end subroutine setup_simulation

   !---------------------------------------------------------------------------------
   !> @brief
   !> Prints timing for the overall large program phases
   !> which is more complete than the profiling used otherwise
   !
   !> @author
   !> Anders Bergman
   !
   !> @date 2014/08/15 - Thomas Nystrand
   !> - Now both execution time and cpu_time is logged
   !---------------------------------------------------------------------------------
   subroutine sd_timing(Nstep)
      !
      implicit none
      !
      integer, intent(in) :: Nstep
      !
      character(30) :: format1
      character(30) :: format2

      format1 = "(1x,a,1x,a,4x,a)"
      format2 = "(2x,a,7x,f8.2,10x,f8.2)"

      write(*,'(1x,a)') '--------------MAIN PROGRAM PHASES TIME REPORT----------------'
      write (*,format1) 'CATEGORY                  ', 'EXECUTION TIME', 'TOTAL CPU TIME'
      write (*,format2)  'Time for initialization :', (time_b - time_a), (time_II  - time_I)
      write (*,format2)  'Time for initial phase  :', (time_c - time_b), (time_III - time_II)
      write (*,format2)  'Time for meas. phase    :', (time_d - time_c), (time_IV  - time_III)
      write (*,format2)  'Time for cleanup        :', (time_e - time_d), (time_V   - time_IV)

      write(*,'(1x,a)') '-------------------------------------------------------------'
      write (*,format2)  "Time total              :", (time_e - time_a), (time_V - time_I)
      write (*,format2)  "Time for one meas. iter :", (time_d - time_c)/Nstep, (time_IV - time_III)/Nstep

   end subroutine sd_timing


   !---------------------------------------------------------------------------------
   !> @brief
   !> Allocate arrays needed for the simulations
   !
   !> @author
   !> Anders Bergman
   !---------------------------------------------------------------------------------
   subroutine allocate_general(flag)
      !
      use LLGI,                  only : allocate_llgifields
      use Depondt,               only : allocate_depondtfields
      use Midpoint,              only : allocate_aux_midpoint_fields
      use InputData
      use FieldData,             only : allocate_fields, read_local_field, allocation_field_time
      use MonteCarlo,            only : mcmavg_buff, indxb_mcavrg
      use Heun_proper,           only : allocate_aux_heun_fields
      use Measurements,          only : allocate_measurementdata
      use RandomNumbers,         only : allocate_randomwork
      use LatticeInputData,      only : do_ld
      use LatticeFieldData,      only : latt_allocate_fields
      use LatticeMeasurements,   only : allocate_lattmeasurementdata
      !
      implicit none
      !
      integer, intent(in) :: flag !< Allocate or deallocate (1/-1)

      integer :: i_stat

      !
      call allocate_fields(Natom,Mensemble,flag)

      if(locfield=='Y'.and.flag>0)  call read_local_field(NA,locfieldfile)
      if(SDEalgh==5 .or. ipSDEalgh==5) then
         call allocate_depondtfields(Natom, Mensemble,flag)
      elseif(SDEalgh==11) then
         call allocate_llgifields(Natom, Mensemble,flag)
      end if

       if (SDEalgh==1 .or. ipSDEalgh==1) then
        if (mode.ne.'MS') then
        call allocate_aux_midpoint_fields(flag,Natom,Mensemble)
        endif
      endif
      if (SDEalgh==4 .or. ipSDEalgh==4) then
        if (ipmode.ne.'MS') then 
        call allocate_aux_heun_fields(flag,Natom,Mensemble)
        endif
      endif

      !if(SDEalgh>=1.and.SDEalgh<=4) call allocate_randomwork(Natom,Mensemble,flag)
      if(SDEalgh>=1.and.SDEalgh<=4.or.SDEalgh==6.or.SDEalgh==7.or.SDEalgh==21.or.SDEalgh==22) &
         call allocate_randomwork(Natom,Mensemble,flag,do_ld)

      call allocate_measurementdata(NA,NT,Natom,Mensemble,Nchmax,plotenergy,flag)

      if(allocated(mcmavg_buff)) then
         allocate(mcmavg_buff(3,mcavrg_buff),stat=i_stat)
         call memocc(i_stat,product(shape(mcmavg_buff))*kind(mcmavg_buff),'mcmavg_buff','allocate_general')
      end if

      if(allocated(indxb_mcavrg)) then
         allocate(indxb_mcavrg(mcavrg_buff),stat=i_stat)
         call memocc(i_stat,product(shape(indxb_mcavrg))*kind(indxb_mcavrg),'indxb_mcavrg','allocate_general')
      end if

      if(do_ld == 'Y') then
         call latt_allocate_fields(Natom,Mensemble,flag)
         call allocate_lattmeasurementdata(Natom, Mensemble, NA, NT, Nchmax, flag)
      end if

   end subroutine allocate_general

   !---------------------------------------------------------------------------------
   !> @brief
   !> Display information about the program
   !
   !> @author
   !> Anders Bergman
   !---------------------------------------------------------------------------------
   subroutine print_logo()

      implicit none

      integer, dimension(8) :: times

      call date_and_time(VALUES=times)

      ! Print logo
      write (*,'(1x, a)')    "--------------------------------------------------------------"
      write (*,'(1x, a)')    "            __  __          ___   _______    ____  ___        "
      write (*,'(1x, a)')    "           / / / /__  ___  / _ | / __/ _ \  / __/ / _ \       "
      write (*,'(1x, a)')    "          / /_/ / _ \/ _ \/ __ |_\ \/ // / / _ \_/ // /       "
      write (*,'(1x, a)')    "          \____/ .__/ .__/_/ |_/___/____/  \___(_)___/        "
      write (*,'(1x, a)')    "              /_/  /_/                                        "
      write (*,'(1x, a)')    "--------------------------------------------------------------"
      write (*,'(1x, a)')    "               https://github.com/UppASD/UppASD               "
      write (*,'(1x, a)')    "--------------------------------------------------------------"
      write (*,'(1x, a)')    "             Division of Materials Theory                     "
      write (*,'(1x, a)')    "             Department of Physics and Astronomy              "
      write (*,'(1x, a)')    "             Uppsala University                               "
      write (*,'(1x, a)')    "             Sweden                                           "
      write (*,'(1x, a)')    "---------------------Production-version-----------------------"
      ! Current logo using the Small Slant font from
      ! http://patorjk.com/software/taag/#p=display&f=Small%20Slant&t=UppASD%206.0
      ! Print git repo version
#if defined(VERSION)
      write (*,'(1x, a,a)')  "Git revision: ", VERSION
#else
      write (*,'(1x, a)')    "Git revision: unknown"
#endif
      write (*,'(1x, a)')    "--------------------------------------------------------------"

      ! Print compiler (deactivated due to pgf90 issues
!#if (!defined __PATHSCALE__) || (!defined __PGIF90__ )
!write (*,'(1x, a,a)')  "Fortran compiler: ", compiler_version()
!#endif
!      write (*,'(1x, a)')    "--------------------------------------------------------------"

      ! Print if MKL RNG is enabled
#ifdef VSL
      write (*,'(1x, a,a)')  "RNG: MKL Vector statistics library "
#else
      write (*,'(1x, a,a)')  "RNG: Mersenne-Twister (mtprng)"
#endif
      write (*,'(1x, a)')    "--------------------------------------------------------------"

      ! Print execution time
      write (*,10010)        "Execution Date: ",times(3),"/",times(2),"/",times(1),&
         "Execution Time: ",times(5),":",times(6),":",times(7),":",times(8)
      !write (*,'(1x, a,a)')  "  with the flags: ", compiler_options()
      write (*,'(1x, a)')    "--------------------------------------------------------------"

      ! Print OpenMP information
      !$omp parallel
      !$omp master
      !nprocs=omp_get_num_threads()
      write(*,'(1x,a18,i2,a16,i3,a10)') &
         "Using OpenMP with ",omp_get_num_threads()," threads out of",OMP_GET_NUM_PROCS(),"possible."
      !$omp end master
      !$omp end parallel
      write (*,'(1x, a)')    "--------------------------------------------------------------"

      10010 format(1x,a,i0.2,a,i0.2,a,i4,2x,a,i0.2,a,i0.2,a,i0.2,a,i0.3)
   end subroutine print_logo


   !---------------------------------------------------------------------------------
   !> @brief
   !> Print brief information about the simulation to stdout
   !
   !> @author
   !> Anders Bergman
   !---------------------------------------------------------------------------------
   subroutine print_siminfo()

      ! restruc later
      use Correlation, only : do_sc
      !use Correlation, only : do_sc, do_conv, sigma_q, sigma_w, LQfactor, LWfactor!, do_sc_proj
      use Temperature
      use SystemData
      use InputData
      use Damping
      use prn_averages
      use prn_trajectories

      implicit none

      integer :: i

      write (*,'(1x,a)')          "------------------------------------------"
      write (*,'(1x,a)')          "           Simulation data                "
      write (*,'(1x,a)')          "------------------------------------------"
      write (*,'(1x,a,i8)')       " Number of atoms", Natom
      write (*,'(1x,a,i4)')       " Number of ensembles", Mensemble
      write (*,'(1x,a,1x,a)')     " Gradient:", grad
      if(grad.eq.'N') then
         write (*,'(1x,a,f8.2,a)')   " Temperature:", Temp, " K"
      else
         write (*,'(1x,a,f8.2,f8.2,f8.2,f8.2,a)')   " Temperature boundaries:", Temp_low_x, Temp_high_x, Temp_low_y, Temp_high_y, " K"
      end if

      write (*,'(1x,a,i10)')      " Number of simulation steps:", nstep
      write (*,'(1x,a,G11.4,a)')  " Time step:", (delta_t*1.0d18), "as"
      write (*,'(1x,a,G11.4,a)')  " Simulation time:", (nstep*delta_t*1.0d12), "ps"
      if (do_site_damping=='Y') then
         write (*,'(2x,a,4x,a)') "Damping_1", "Site"
         do i=1, NA
            write (*,'(2x,G11.4,i3)') lambda1_array(i),i
         end do
      else
         write (*,'(1x,a,G11.4)')    " Damping parameter:", lambda1_array(1)
      endif
      write (*,'(1x,a,1x,a)')     " Sample averages:", do_avrg
      write (*,'(1x,a,1x,a)')     " Sample moments:", do_tottraj
      write (*,'(1x,a,1x,a)')     " Spin correlation:", do_sc
      write (*,'(1x,a)')          "------------------------------------------"
      write (*,'(1x,a)')          "        Progress of simulation:           "

   end subroutine print_siminfo

   !---------------------------------------------------------------------------------
   !> @brief
   !> Checks if the inpsd.dat file is using old legacy format
   !
   !> @author
   !> Anders Bergman
   !
   !> @date 11/08/2014 - Thomas Nystrand
   !> - Moved oldformat stop execution here
   !---------------------------------------------------------------------------------
   subroutine check_format()

      implicit none

      logical:: oldformat
      character(len=2) :: teststring
      open(ifileno, file='inpsd.dat')
      read(ifileno,'(a)') teststring
      rewind(ifileno)
      close(ifileno)
      oldformat=(teststring=="**")

      ! Old format
      if (oldformat) then
         write(*,'(1x,a)',advance='no') 'using old format - No longer supported'
         stop
      end if
      return
   end subroutine check_format


   subroutine calculate_energy(outenergy)
      use QHB,                   only : qhb_rescale
      use FixedMom
      use InputData
      use FieldData,             only : external_field,time_external_field
      use FieldPulse
      use SystemData
      use DemagField
      use MomentData
      use macrocells
      use Temperature
      use MicroWaveField
      use SimulationData,        only : bn
      use HamiltonianData
      use CalculateFields
      use HamiltonianActions
      use optimizationRoutines
      use AdaptiveTimestepping
      use Energy
      use Simulationdata, only : total_energy

      real(dblprec), optional :: outenergy

      call calc_energy(nHam,0,Natom,Nchmax, &
         conf_num,Mensemble,Natom,Num_macro,1,         &
         plotenergy,Temp,delta_t,do_lsf,        &
         lsf_field,lsf_interpolate,real_time_measure,simid,cell_index,            &
         macro_nlistsize,mmom,emom,emomM,emomM_macro,external_field,              &
         time_external_field,max_no_constellations,maxNoConstl,                   &
         unitCellType,constlNCoup,constellations,OPT_flag,                        &
         constellationsNeighType,total_energy,NA,N1,N2,N3)

      if(present(outenergy)) outenergy=total_energy
   end subroutine calculate_energy

   end module uppasd
