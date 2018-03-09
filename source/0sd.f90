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


!------------------------------------------------------------------------------
! MODULE: 0sd
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
!------------------------------------------------------------------------------
program SD
    !
    use Parameters
    use ErrorHandling
    use Profiling
    use omp_lib
    use InputData, only: Nstep
    !
    use SD_driver
    use MC_driver

   !.. Implicit declarations
   implicit none

   ! Variables for timing and number of processors
   real(dblprec) :: time_a, time_b, time_c, time_d, time_e
   real(dblprec) :: time_I, time_II, time_III, time_IV, time_V
   integer :: nprocs


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
   call timing(0,'Initial       ','ON')
   call run_initial_phase()
   call timing(0,'Initial       ','OF')
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

   !==============================================================!

contains


   !---------------------------------------------------------------------------
   !> @brief
   !> Find number of OpenMP processors active
   !
   !> @author
   !> Anders Bergman
   !
   !> @date 11/08/2014 - Thomas Nystrand
   !> - Moved to separate function
   !---------------------------------------------------------------------------
   function number_of_active_processors() result(nprocs)
      use omp_lib
      integer :: nprocs

      !$omp parallel
      !$omp master
      nprocs=omp_get_num_threads()
      write(*,'(1x,a18,i2,a16,i3,a10)') &
         "Using OpenMP with ",nprocs," threads out of",OMP_GET_NUM_PROCS(),"possible."
      !$omp end master
      !$omp end parallel
   end function number_of_active_processors


    !---------------------------------------------------------------------------
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
    !---------------------------------------------------------------------------
    subroutine run_initial_phase()

       use OptimizationRoutines
       use SystemData
       use InputData
       use HamiltonianData
       use MomentData
       use QMinimizer

       implicit none

       write (*,'(1x,a)') "Enter initial phase:"

       if (ipmode=='M' .or. ipmode=='H'.or.ipmode=='I'.or.ipmode=='L'.or.ipmode=='W'.or.ipmode=='D') then
          ! Monte Carlo initial phase
          call mc_iphase()
       elseif (ipmode=='S') then
          ! Spin dynamics initial phase
          call sd_iphase()
       elseif (ipmode=='Q') then
          ! Spin spiral minimization initial phase
                      call ErrorHandling_missing('Q-minimization')

       elseif (ipmode=='G') then
         call ErrorHandling_missing('VPO-minimization')
       endif
    end subroutine run_initial_phase



    !---------------------------------------------------------------------------
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
    !---------------------------------------------------------------------------
    subroutine run_measurement_phase()

       use BLS,               only: calc_bls, deallocate_bls_data
       use Correlation,       only: correlation_wrapper, sc_nstep
       use AutoCorrelation,   only: do_autocorr, allocate_autocorr, autocorr_init
       use OptimizationRoutines
       use SystemData
       use InputData
       use HamiltonianData
       use MomentData
       use ChemicalData
       use SimulationData
       use QMinimizer
       use Restart

       integer :: cflag, cflag_p


       write (*,'(1x,a)') "Enter measurement phase:"
       ! Allocate and initialize spin arrays for autocorrelation
       call timing(0,'SpinCorr      ','ON')
       if(do_autocorr=='Y') call allocate_autocorr(Natom)
       if(do_autocorr=='Y') call autocorr_init(Natom, Mensemble, simid, rstep, initmag, emom)
       call timing(0,'SpinCorr      ','OF')


       if(mode=='M' .or. mode=='H'.or.mode=='I'.or.mode=='L'.or.mode=='W'.or.mode=='D') then
          call mc_mphase() ! Monte Carlo measurement phase

       elseif (mode=='S') then
          call print_siminfo()
          if (gpu_mode==0) then !FORTRAN
             call sd_mphase() ! Spin Dynamics measurement phase
          else ! C++ or CUDA
             call sd_mphaseCUDA()
          endif

       elseif (mode=='Q') then
             call ErrorHandling_missing('Q-minimization')


       elseif (mode=='G') then
             call ErrorHandling_missing('GNEB')

       endif

       ! Save moment information from final simulation step
       mstep = mstep-1
       call timing(0,'PrintRestart  ','ON')
       call prnrestart(Natom, Mensemble, simid, mstep, emom, mmom)
       call timing(0,'PrintRestart  ','OF')

       call timing(0,'SpinCorr      ','ON')
       cflag = 2 ; cflag_p = 2

       ! Spin-correlation
       call correlation_wrapper(Natom, Mensemble, coord, simid, emomM, mstep, delta_t, NT, atype, Nchmax, achtype, cflag, cflag_p)

       ! Lattice-correlation
       call calc_bls(N1,N2,N3,C1,C2,C3,Natom, Mensemble, simid,coord, emomM, mstep,delta_t, cflag)
       call deallocate_bls_data()

       call timing(0,'SpinCorr      ','OF')

    end subroutine run_measurement_phase



    !---------------------------------------------------------------------------
    !> @brief
    !> Finish simulation and deallocate leftovers
    !
    !> @author
    !> Anders Bergman
    !
    !> @date 11/08/2014 - Thomas Nystrand
    !> - Moved to separate function
    !---------------------------------------------------------------------------
    subroutine cleanup_simulation()
      use Correlation
      use SpinIceData, only : allocate_spinice_data, allocate_vertex_data
      use InputData
      use SystemData
      use HamiltonianData
      use MomentData
      use ChemicalData
      use SpinTorques
      use OptimizationRoutines
      use MicroWaveField
      use Polarization
      use Damping
      use Tools
      use LSF, only : allocate_lsfdata

      write (*,'(1x,a)') "Simulation finished"
      call allocate_systemdata(flag=-1)
      call allocate_anisotropies(Natom,mult_axis,flag=-1)
      call allocate_hamiltoniandata(do_jtensor=do_jtensor,do_lsf=do_lsf,ind_mom_flag=ind_mom_flag,flag=-1,&
         lsf_field=lsf_field,exc_inter=exc_inter)
      if(do_dm==1) call allocate_dmhamiltoniandata(flag=-1)
      if(do_pd==1) call allocate_pdhamiltoniandata(flag=-1)
      if(do_bq==1) call allocate_bqhamiltoniandata(flag=-1)
      call allocate_mmoms(flag=-1)
      call allocate_general(flag=-1)
      if(do_ralloy==1.or.oldformat) call allocate_chemicaldata(flag=-1)
      if(do_ralloy==1.or.oldformat) call allocate_chemicalinput(flag=-1)
      if(do_lsf=='Y') call allocate_lsfdata(NA,Nchmax,conf_num,flag=-1)
      call allocate_initmag(flag=-1)
      call deallocate_q() ! Deallocate spin correlation related arrays
      call deallocate_rest() ! Deallocate remaining arrays

      call allocate_deltatcorr(.false.)

      if (mode=='L') then
        call allocate_spinice_data(Natom,Mensemble,-1)
        call allocate_vertex_data(-1)
      end if

      ! Deallocate dampings
      call allocate_damping(Natom,ipnphase,-1)
      ! Deallocate trajectory measurements
      ! Deallocate Microwave fields
      call setup_mwf_fields(Natom,-1)

      ! Deallocate the btorque if it is present
      call allocate_stt_data(Natom,Mensemble,flag=-1)
      ! --Optimization Region-- !
      if(OPT_ON) then
         call allocateOptimizationStuff(na,natom,Mensemble,.false.)
      end if
    end subroutine cleanup_simulation


   !---------------------------------------------------------------------------
   !> @brief
   !> Initialize simulation by setting up Hamiltonian and magnetic moments
   !
   !> @author
   !> Anders Bergman
   !---------------------------------------------------------------------------
   subroutine setup_simulation()
      !
      use SystemData
      use SimulationData
      use InputData
      use HamiltonianData
      use MomentData
      use ChemicalData
      use Damping
      use RandomNumbers
      use OptimizationRoutines
      use Correlation
      use Temperature
      use AutoCorrelation
      use SpinTorques
      use Topology
      use MC_Wolff
      use prn_micromagnetic
      use Stiffness,  only : do_stiffness
      use Polarization
      use BLS
      use Damping
      use InputHandler
      use ams
      use fieldpulse
      use geometry
      use gradients
      use hamiltonianinit
      use inputhandler
      use magnetizationinit
      use microwavefield
      use printinput
      use prn_fields
      use prn_microwaves
      use prn_topology
      use setupspinice
      use prn_currents
      use prn_averages
      use prn_trajectories
      use LSF, only : read_LSF,allocate_lsfdata
      use InducedMoments, only : calculate_ind_mom_linear

      implicit none

      integer :: i, i_stat,mconf

      ! Read inpsd.dat input file
      write(*,'(1x,a)',advance='no') 'Read input file'


      !Check which input format is used
      call check_format()


      write(*,'(1x,a)') 'using new format'
      ! Set the defaults for the BLS analysis
      call setup_bls()
      ! Set the defaults for the trajectories
      call traj_init()
      ! Set the defaults for the averages, cumulants and autocorrelation
      call avrg_init()
      ! Set the default for the polarization measurements
      call prn_pol_init()
      ! Set defaults for vertex printing
      call spin_ice_init()
      ! Set defaults for the field printing
      call fields_prn_init()
      ! Set input defaults for the correlation
      call correlation_init()
      ! Set input defaults for topological information
      call prn_topology_init()
      ! Set general input defaults
      call set_input_defaults()
      ! Set microwaves fields printing defaults
      call initialize_prn_microwaves()
      ! Set defaults for Adiabatic Magnon Spectra calculations
      call init_ams()
      ! Set defaults for Spin Transfer Torques and related effects
      call init_stt()
      ! Set defaults for ferromagnetic stiffness calculations
      call init_stiffness()
      ! Set defaults for current measurements
      call current_init()
      ! Set defaults for random number generation
      call rng_defaults()
      ! Set defaults for temperature routines
      call set_temperature_defaults()
      ! Set defaults for Optimization routines
      call set_opt_defaults(delta_t)

      open(ifileno,file='inpsd.dat')
      call read_parameters(ifileno)
      rewind(ifileno)
      call read_parameters_correlation(ifileno)
      rewind(ifileno)
      call read_parameters_bls(ifileno)
      rewind(ifileno)
      call read_parameters_mwf(ifileno)
      rewind(ifileno)
      call read_parameters_averages(ifileno)
      rewind(ifileno)
      call read_parameters_autocorr(ifileno)
      rewind(ifileno)
      call read_parameters_ams(ifileno)
      rewind(ifileno)
      call read_parameters_stt(ifileno)
      rewind(ifileno)
      close(ifileno)

      !----------------------------!
      ! Reading in data from files !
      !----------------------------!

      ! Read atom positions for unit cell or full input
      if(aunits=='Y') then
         write(*,*) 'Use atomic units'
         call change_constants()
      end if
      if(.not.allocated(Bas)) then
         if(do_ralloy==0) then
            call read_positions()
         else
            call read_positions_alloy()
            Natom_full=NA*N1*N2*N3
            call allocate_chemicaldata(Natom_full,1)
         end if
      end if
      Natom_full=NA*N1*N2*N3

      call read_moments(Landeg_glob)

      ! Read LSF related items
      if (do_lsf == 'Y') then
            call ErrorHandling_missing('Local spin fluctuations')
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

      ! Read exchange interaction from file, either normal or tensor format
      if(do_jtensor==0) then
         ! Normal exchange interaction read
         call read_exchange()
      elseif(do_jtensor==1 .and. calc_jtensor.eqv..false.) then
         ! Read tensor exchange format file for tensor coupling
         call read_exchange_tensor()
      elseif(do_jtensor==1 .and. calc_jtensor.eqv..true.) then
         ! Read normal exchange file which is thereafter used to build tensor coupling
         call read_exchange_build_tensor()
      else
         call ErrorHandling_ERROR("Unrecognized or unsupported combination of do_tensor &
            &and calc_tensor")
      end if
      max_no_shells=maxval(NN)


      if(do_dm==1)    call read_dmdata()
      if(do_pd==1)    call read_pddata()
      if(do_biqdm==1) call read_biqdmdata()
      if(do_bq==1)    call read_bqdata()

      if(do_autocorr=='Y') call read_tw(twfile)

      if(do_bpulse>=1.and.do_bpulse<5) call read_bpulse(do_bpulse,bpulsefile)

      if(do_anisotropy==1) then
         if(do_ralloy==0) then
            call read_anisotropy()
         else
            call read_anisotropy_alloy()
         end if
      end if

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

      ! Set up geometry containing chemical data
      call setup_geometry(Natom, NA, N1, N2, N3, Bas, C1, C2, C3,&
           coord, atype, anumb, atype_inp, anumb_inp, do_prnstruct, &
           tseed, simid, do_ralloy, Natom_full, Nchmax, Nch, acellnumb, acellnumbrev, achtype, &
           chconc, atype_ch, asite_ch, achem_ch)

      ! Allocating the damping arrays
      call allocate_damping(Natom,ipnphase,0)

      write(*,'(1x,a)') ' Setup up damping'
      call setup_damping(NA,Natom,Natom_full,ipmode,mode,ipnphase,&
           do_ralloy,asite_ch,achem_ch,do_site_damping,&
           do_site_ip_damping,iplambda1,iplambda2,mplambda1,mplambda2,&
           iplambda1_array,iplambda2_array,lambda1_array,lambda2_array)

      write(*,'(1x,a)') ' done'

      ! See if it is necesary to read the temperature file
      !if(grad.eq.'Y'.or.grad.eq.'F') call read_temperature()
      if(grad.eq.'Y'.or.grad.eq.'F') call read_temperature_legacy()

      ! Allocating the temperature arrays
      call allocate_temp(Natom,ipnphase)

      write(*,'(1x,a)') ' Setup up temperature'
      if (ipnphase.ge.1) then
         do i=1, ipnphase
            if(grad=='Y') then 
               !call SETUP_TEMP(NATOM,NT,NA,N1,N2,N3,NATOM_FULL,&
               !   DO_RALLOY,ATYPE,ACELLNUMB,ATYPE_CH,SIMID,iptemp(i),&
               !   C1,C2,C3,BC1,BC2,BC3,BAS,COORD,iptemp_array(:,i))
               call Lparray(ipTemp_array(:,i),Natom,coord,iptemp(i),simid,.false.)
            else
               iptemp_array(:,i)=iptemp(i)
            end if

         enddo
      endif
      if(grad=='Y') then
         !call SETUP_TEMP(NATOM,NT,NA,N1,N2,N3,NATOM_FULL,&
         !   DO_RALLOY,ATYPE,ACELLNUMB,ATYPE_CH,SIMID,TEMP,&
         !   C1,C2,C3,BC1,BC2,BC3,BAS,COORD,temp_array)
         call Lparray(Temp_array,Natom,coord,Temp,simid,.true.)

      else
         temp_array=TEMP
      end if

      write(*,'(1x,a)') ' done'
      print *,iptemp,TEMP

      if (mode=='L') then
            call ErrorHandling_missing('Spin-ice')

      endif

      ! Print input
      write(*,'(1x,a)',advance='no') 'Write input data'
      ! Modify for SKKR stuff later
      call prninp()
      write(*,'(a)') ' done.'
      !
      ! Set up q-point grid if real space correlation function is to be calculated
      if(do_sc=='C' .or. do_sc=='Q' .or. do_ams=='Y' ) then
         if (qpoints == 'C' .or. qpoints == 'P' .or. qpoints == 'A'.or.qpoints=='X'.or.qpoints=='H') then
            write (*,'(1x,a)') "Calls setup_qcoord"
            call setup_qcoord(N1, N2, N3, C1, C2, C3)
         else if(qpoints == 'F' .or. qpoints == 'D'.or.qpoints=='I'.or.qpoints=='B') then
            write (*,'(1x,a)') "Calls read_q"
            call read_q(C1, C2, C3)
         end if
         if(do_sc=='Q') then
            call set_w(delta_t)
         end if
         call allocate_deltatcorr(.true.)
      end if

      write (*,'(1x,a)') "Set up Hamiltonian"
      ! Set up Hamiltonian, containing exchange, anisotropy, and optional terms
      ! like DM and dipolar interactions.
      call setup_hamiltonian(NT,NA,N1,N2,N3,Nchmax,do_ralloy,Natom_full,Natom,acellnumb,acellnumbrev,&
           achtype,atype_ch,asite_ch,achem_ch,atype,anumb,alat,C1,C2,C3,Bas,ammom_inp,coord,BC1,BC2,BC3,&
           sym,do_jtensor,max_no_neigh,max_no_shells,nn,nlistsize,nlist,redcoord,jc,jcD,jc_tens,ncoup,&
           ncoupD,j_tens,do_dm,max_no_dmshells,max_no_dmneigh,dm_nn,dmlistsize,dmlist,dm_vect,dm_redcoord,&
           dm_inpvect,do_anisotropy,taniso,taniso_diff,random_anisotropy_density,anisotropytype,&
           anisotropytype_diff,anisotropy,anisotropy_diff,sb,sb_diff,eaniso,eaniso_diff,kaniso,&
           kaniso_diff,random_anisotropy,mult_axis,mconf,conf_num,fs_nlistsize,fs_nlist,nind,&
           map_multiple,do_lsf,lsf_field,exc_inter,do_bq,nn_bq_tot,max_no_bqshells,bq_nn,bqlistsize,&
           bqlist,bq_redcoord,jc_bq,j_bq,do_biqdm,nn_biqdm_tot,max_no_biqdmshells,biqdm_nn,&
           biqdmlistsize,biqdmlist,biqdm_redcoord,biqdm_inpvect,biqdm_vect,do_pd,nn_pd_tot,&
           max_no_pdshells,pd_nn,pdlistsize,pdlist,pd_redcoord,pd_inpvect,pd_vect,do_dip,qdip,&
           ind_mom,ind_nlistsize,ind_nlist,ind_tol,ind_mom_flag,&
           do_prnstruct,do_sortcoup,simid)

      ! Allocate arrays for simulation and measurement
      call allocate_general(1)


      write (*,'(1x,a)') "Initialize magnetic moments "

      ! Allocate arrays for magnetic moment magnitudes
      call allocate_mmoms(Natom, Mensemble,1)

      ! Fill arrays with magnetic moment magnitudes from input
      call setup_moment(Natom, conf_num, Mensemble, NA, N1, N2, N3, Nchmax, ammom_inp, Landeg_ch, Landeg, mmom, mmom0, mmomi,&
           do_ralloy, Natom_full, achtype, acellnumb,mconf)

      ! Set up magnetic moment vectors
      call magninit(Natom, conf_num, Mensemble, NA, N1, N2, N3, initmag, Nchmax, aemom_inp, anumb, do_ralloy, &
           Natom_full, achtype, acellnumb, emom, emom2, emomM, mmom, rstep, theta0, phi0, restartfile, &
           initrotang, initpropvec, initrotvec, coord, C1, C2, C3)!!!,do_fixed_mom,Nred,red_atom_list)

      ! Rotation of the initial spin configuration
      if (roteul >= 1) then
         write(*,'(1x,a)',advance='no') 'Perform Euler rotation'
         call rotationeuler(Natom, Mensemble,roteul, rotang,emom, emomM, mmom, emom2)
         write (*,'(a)') ' done.'
      end if

      ! Excitation of the initial spin configuration
      if (initexc == 'I') then
         write(*,'(1x,a,f10.4,a)',advance='no') 'Introduces ', initconc*100, '% vacancies'
         call setinitexc(Natom, Mensemble, emom, emomM, mmom, mmom2, emom2, &
              initexc, initconc, initneigh, initimp, max_no_neigh, nlist, mseed)
         write (*,'(a)') ' done.'
      else if (initexc == 'R') then
         write(*,'(1x,a,f10.4,a)',advance='no') 'Introduces ', initconc*100, '% two-magnon Raman spin flips'
         call setinitexc(Natom, Mensemble, emom, emomM, mmom, mmom2, emom2, &
              initexc, initconc, initneigh, initimp, max_no_neigh, nlist, mseed)
         write (*,'(a)') ' done.'
      end if

      if (stt/='N'.or.do_she/='N') then

         call read_jvecfile(Natom)

         if (stt=='A'.or.(stt/='A'.and.skyno=='Y')) call setup_stencil_mesh(Natom, N1, N2, N3, C1, C2, C3, BC1, BC2, BC3, max_no_neigh, nlistsize, nlist, coord)
         ! Call to allocate the needed stt data
         call allocate_stt_data(Natom,Mensemble,flag=1)
         ! Call the printing of the current density in proper units
         call print_curr_density(NA,Nchmax,conf_num,alat,spin_pol,C1,C2,C3,jvec,ammom_inp)

      endif
      if (stt=='N'.and.(skyno=='Y'.or.do_proj_skyno=='Y')) then
         call setup_stencil_mesh(Natom, N1, N2, N3, C1, C2, C3, BC1, BC2, BC3, max_no_neigh, nlistsize, nlist, coord)
      end if

      if (mode=='W') then
            call ErrorHandling_missing('Wolff algorithm')

         endif

      call setup_mwf_fields(Natom,1)

      if(do_bpulse==6) then
         call read_sitefield(Natom,sitenatomfld)
      endif

      ! --Optimization Region-- !
      if(OPT_flag) then
         call allocateOptimizationStuff(na,natom,Mensemble,.true.)
         call getCellPosNumNeigh(Natom, NA, nlistsize, cellPosNumNeigh)
         OPT_ON=.true.
      end if

      ! Check AMS-flag
      if (do_ams =='Y') then
         call calculate_ams()
         if(do_ralloy>0) call calculate_random_ams()
      end if

      ! Check if the micromagnetic information is calculated
      if (do_stiffness =='Y') then
        call stiffness_wrapper(NT,NA,N1,N2,N3,1,mconf,Natom,Nchmax,conf_num,do_ralloy,&
             Natom_full,max_no_neigh,Nch,anumb,atype,nlistsize,atype_ch,asite_ch,&
             achem_ch,nlist,alat,C1,C2,C3,coord,chconc,ammom_inp,ncoup,max_no_dmneigh,&
             dmlistsize,dmlist,dm_vect,do_anisotropy,anisotropy,simid)
      end if

      ! Deallocate input data for Heisenberg Hamiltonian
      call allocate_hamiltonianinput(flag=-1)


   end subroutine setup_simulation


   !---------------------------------------------------------------------------
   !> @brief
   !> Prints timing for the overall large program phases
   !> which is more complete than the profiling used otherwise
   !
   !> @author
   !> Anders Bergman
   !
   !> @date 2014/08/15 - Thomas Nystrand
   !> - Now both execution time and cpu_time is logged
   !---------------------------------------------------------------------------
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


    !---------------------------------------------------------------------------
    !> @brief
    !> Allocate arrays needed for the simulations
    !
    !> @author
    !> Anders Bergman
    !---------------------------------------------------------------------------
    subroutine allocate_general(flag)
       !
       use LLGI,          only : allocate_llgifields
       use Depondt,       only : allocate_depondtfields
       use InputData
       use FieldData,     only : allocate_fields, read_local_field, allocation_field_time
       use MonteCarlo,    only : mcmavg_buff, indxb_mcavrg
       use Measurements,  only : allocate_measurementdata
       use RandomNumbers, only : allocate_randomwork
       !
       implicit none
       !
       integer, intent(in) :: flag !< Allocate or deallocate (1/-1)

       integer :: i_stat

       !
       call allocate_fields(Natom,Mensemble,flag)

       if(locfield=='Y'.and.flag>0)  call read_local_field(NA,locfieldfile)
       if(SDEalgh==5) then
          call allocate_depondtfields(Natom, Mensemble,flag)
       elseif(SDEalgh==11) then
          call allocate_llgifields(Natom, Mensemble,flag)
       end if

       if(SDEalgh>=1.and.SDEalgh<=4.or.SDEalgh==6.or.SDEalgh==7.or.SDEalgh==21.or.SDEalgh==22) &
          call allocate_randomwork(Natom,Mensemble,flag,'N')

       call allocate_measurementdata(NA,NT,Natom,Mensemble,Nchmax,plotenergy,flag)

       if(allocated(mcmavg_buff)) then
          allocate(mcmavg_buff(3,mcavrg_buff),stat=i_stat)
          call memocc(i_stat,product(shape(mcmavg_buff))*kind(mcmavg_buff),'mcmavg_buff','allocate_general')
       end if

       if(allocated(indxb_mcavrg)) then
          allocate(indxb_mcavrg(mcavrg_buff),stat=i_stat)
          call memocc(i_stat,product(shape(indxb_mcavrg))*kind(indxb_mcavrg),'indxb_mcavrg','allocate_general')
       end if

    end subroutine allocate_general


   !---------------------------------------------------------------------------
   !> @brief
   !> Display information about the program
   !
   !> @author
   !> Anders Bergman
   !---------------------------------------------------------------------------
   subroutine print_logo()
      implicit none

      write (*,'(/,1x, a)')  "--------------------------------------------------------------"
      write (*,'(1x, a)')    "            __  __          ___   _______    ____  ___        "
      write (*,'(1x, a)')    "           / / / /__  ___  / _ | / __/ _ \  / __/ / _ \       "
      write (*,'(1x, a)')    "          / /_/ / _ \/ _ \/ __ |_\ \/ // / /__ \_/ // /       "
      write (*,'(1x, a)')    "          \____/ .__/ .__/_/ |_/___/____/ /____(_)___/        "
      write (*,'(1x, a)')    "              /_/  /_/                                        "
      write (*,'(1x, a)')    "--------------------------------------------------------------"
      write (*,'(1x, a)')    "             Division of Materials Theory                     "
      write (*,'(1x, a)')    "             Department of Physics and Astronomy              "
      write (*,'(1x, a)')    "             Uppsala University                               "
      write (*,'(1x, a)')    "             Sweden                                           "
      write (*,'(1x, a, /)') "----------------------Release-version-------------------------"
      ! Current logo using the Small Slant font from
      ! http://patorjk.com/software/taag/#p=display&f=Small%20Slant&t=UppASD%205.0

   end subroutine print_logo



   !---------------------------------------------------------------------------
   !> @brief
   !> Print brief information about the simulation to stdout
   !
   !> @author
   !> Anders Bergman
   !---------------------------------------------------------------------------
   !subroutine print_siminfo(Natom, Mensemble, nstep, delta_t, Temp, lambda1_array, do_avrg, do_tottraj, do_site_damping)
   subroutine print_siminfo()

      ! restruc later
      use Correlation, only : do_sc, do_conv, sigma_q, sigma_w, LQfactor, LWfactor, do_sc_proj
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
      write (*,'(1x,a,1x,a)')     " Sublattice projection of spin correlation:", do_sc_proj
      write (*,'(1x,a,1x,a)')     " Convolution for magnon DOS:", do_conv
      if (do_conv=='GW') then
         write(*,'(1x,a,f8.2)') " Sigma W:", sigma_w
      else if (do_conv=='GQ') then
         write(*,'(1x,a,f8.2)') " Sigma Q:", sigma_q
      else if (do_conv=='GY') then
         write(*,'(1x,a,f8.2,f8.2)') " Sigma W, Sigma Q:", sigma_w, sigma_q
      else if (do_conv=='LW') then
         write(*,'(1x,a,f8.2)') " Lorentz Factor W:", LWfactor
      else if (do_conv=='LQ') then
         write(*,'(1x,a,f8.2)') " Lorentz Factor Q:", LQfactor
      else if (do_conv=='LY') then
         write(*,'(1x,a,f8.2,f8.2)') " Lorentz Factor W, Factor Q:", LWfactor, LQfactor
      end if
      write (*,'(1x,a)')          "------------------------------------------"
      write (*,'(1x,a)')          "        Progress of simulation:           "

   end subroutine print_siminfo



   !---------------------------------------------------------------------------
   !> @brief
   !> Checks if the inpsd.dat file is using old legacy format
   !
   !> @author
   !> Anders Bergman
   !
   !> @date 11/08/2014 - Thomas Nystrand
   !> - Moved oldformat stop execution here
   !---------------------------------------------------------------------------
   subroutine check_format()
      implicit none
      logical:: oldformat
      character(len=2) :: teststring
      open(ifileno, file='inpsd.dat')
      read(ifileno,'(a2)') teststring
      close(ifileno)
      oldformat=(teststring=="**")

      ! Old format
      if (oldformat) then
         write(*,'(1x,a)',advance='no') 'using old format - No longer supported'
         stop
      end if
      return
   end subroutine check_format




end program SD
