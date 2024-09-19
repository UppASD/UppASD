!> authors
!> Edgar Mendez
!> Nikos Ntallis
!> Manuel Pereiro 
!> Nastaran Salehi

 module MultiscaleSetupSystem
  use Parameters, only : dblprec
  use Profiling, only : memocc
  
  use Multiscale, only : runMultiscaleSetup, MultiscaleSetup, MultiscaleRegions, &
     deallocateSetup, multiscaleBackbuffer, multiscaleBackbufferHead

  integer, parameter :: OUTPUT_FILE_LEN = 48
  
  private
  public &
       setupExchanges, setupInterpolation, setupAnisotropy, setupDampingBand, &
       setup_multiscale_system, initializeMomentData, allocate_multiscale

contains 

  pure function rescale_inputed_exchange(inputed_exchange, this_atom_moment_magnitude, other_atom_moment_magnitude, &
            spin_interaction_order) result(rescale_exchange)
    use Constants
  implicit none
    real(dblprec), dimension(:), intent(in) :: inputed_exchange
    real(dblprec), intent(in) :: this_atom_moment_magnitude
    real(dblprec), intent(in) :: other_atom_moment_magnitude
    integer, intent(in) :: spin_interaction_order
    
    real(dblprec), dimension(ubound(inputed_exchange, 1)) :: rescale_exchange

    real(dblprec) :: scale_factor
    scale_factor = 2*mry/mub
    if(this_atom_moment_magnitude*other_atom_moment_magnitude < 1e-6) then
          rescale_exchange = 0.d0
    else
          rescale_exchange  = inputed_exchange * scale_factor &
             / (this_atom_moment_magnitude*other_atom_moment_magnitude)**spin_interaction_order
    endif
  end function rescale_inputed_exchange



  subroutine setup_multiscale_system(config_file)
    use Configuration, only : readMultiscaleOptions, MultiscaleOptions, &
         deallocateOptions
    use MultiscaleInterpolation
    use MultiscaleDampingBand
    use MC_WOLFF,             only : create_wolff_neigh_list
    use Damping
    use Geometry,             only : setup_geometry
    use Ewaldmom,             only: allocate_Ewald
    use InputData            
    use PrintInput
    use prn_fields,           only : fields_prn_init
    use MomentData           
    use FieldPulse,           only : read_bpulse
    use SystemData,           only : atype, anumb, coord, Landeg
    !use SpinIceData,          only : spin_ice_init
    use Temperature
    !use Correlation        
    use Polarization,         only : prn_pol_init
    use SetupSpinIce
    use InputHandler_ext      
    use prn_averages,         only : avrg_init
    use ChemicalData          !only : acellnumb, acellnumbrev, achtype, atype_ch, asite_ch, achem_ch, chconceff
    use RandomNumbers         !only : mt_ran_init_b, mt_ran_init_c, setup_rng_hb, mt_ran_b
    use SimulationData,       only : rstep
    use MicroWaveField,       only : setup_mwf_fields
    use prn_microwaves,       only : initialize_prn_microwaves
    use HamiltonianInit,      only : setup_hamiltonian
    use AutoCorrelation,      only : read_tw
    use HamiltonianData
    use prn_trajectories,     only : traj_init
    use MagnetizationInit    !only : setup_moment, rotationeuler, magninit, setinitexc, loadrestart
    use Restart
    use optimizationRoutines
    use MagnetizationInit
    use Measurements
    !use Ams
    use Gradients
    use Spintorques
    !use energy_barriers
    use stiffness
    !use interpolation
    use LatticeInputData
    use LatticeInputHandler
    use LatticeData   !,          only : allocate_latticedata, uvec, uvec2, vvec, vvec2
    use LatticeInit    !,          only : lattrotationeuler, lattinit
    use LatticeHamiltonianData
    use LatticeHamiltonianInit     !,      only : setup_latticehamiltonian
    use omp_lib
    
    implicit none

    character(len=*), intent(in) :: config_file
    
    integer :: i_stat, i, j

    integer, dimension(:), allocatable :: nlistsizeTmp
    integer, dimension(:, :), allocatable :: nlistTmp
    real(dblprec), dimension(:, :, :), allocatable :: ncoupTmp

    integer, dimension(:), allocatable :: dmlistsizeTmp
    integer, dimension(:, :), allocatable :: dmlistTmp
    real(dblprec), dimension(:, :, :), allocatable :: dmTmp
    
    integer :: isize

    type(MultiscaleOptions) :: options
    type(MultiscaleSetup) :: setup
    type(MultiscaleRegions) :: regions

    real(dblprec) :: alpha, beta
    
    ! Output file names
    character(len=OUTPUT_FILE_LEN) :: intp_file_name
    character(len=OUTPUT_FILE_LEN) :: coord_file_name
    character(len=OUTPUT_FILE_LEN) :: dband_file_name
    character(len=OUTPUT_FILE_LEN) :: exchange_file_name
    character(len=OUTPUT_FILE_LEN) :: dm_file_name
    character(len=OUTPUT_FILE_LEN) :: grad_link_file_name
    character(len=OUTPUT_FILE_LEN) :: regions_file_name
    
    coord_file_name     = "coord."     // simid // ".out"
    intp_file_name      = "interface." // simid // ".out"
    dband_file_name     = "dband."     // simid // ".out"
    exchange_file_name  = "exchange."  // simid // ".out"
    dm_file_name        = "dm."        // simid // ".out"
    grad_link_file_name = "gradlink."  // simid // ".out"
    regions_file_name   = "msregions."   // simid // ".out"
    
    print *,"Reading multiscale options from '",trim(adjustl(config_file)),"'"
    
    call readMultiscaleOptions(trim(adjustl(config_file)),options)
    call runMultiscaleSetup(options,setup,regions)
    call deallocateOptions(options)

    if (do_prnmultiscale) then
       call printMultiscaleRegions(regions, regions_file_name)
    end if

    if(ham_inp%do_anisotropy /= 0 .and. .not. associated(setup%anisotropies%anisotropyTypes))then
       print *,"WARNING: do_anisotropy is set, but no anisotropies"
       print *,"         are generated by the setup."
       print *,"         The flag do_anisotropy will be ignored."
       ham_inp%do_anisotropy = 0
    end if

    !if (do_prnstruct /= 0 .or. prn_vertices /= 'N') then
    !   print *, "WARNING: do_prnstruct and prn_vertices are ignored with multiscale 1."
    !   print *, "         To see the structure of the space please check "
    !   print *, "         MUASD output files."
    !   prn_vertices = 'N'
    !   do_prnstruct = 0
    !end if

    if (allocated(jfile)) then
       isize = product(shape(jfile))*kind(jfile);
       deallocate(jfile, stat=i_stat)
       call memocc(i_stat, -isize, 'jfile', 'setup_multiscale_system')
    end if
    !*!!! Initialize random number generators
     call setup_rng_hb(tseed,ziggurat,rngpol) ! temperature
     
     !call mt_ran_init_b(mseed) ! magninit

    !* Note that the same seed value is used for
    !* Monte-Carlo as for the temperature to avoid an extra input variable
    !*!!! !$omp parallel
    !*!!!if(para_rng) tseed=tseed+omp_get_thread_num()
    !*!!!call mt_ran_init_c(tseed) ! Monte-Carlo
    !*!!!!!$omp end parallel

    coord = setup%positions

    ! Required values
    Natom = ubound(coord, 2)
    Natom_full = Natom
    NA = Natom
    NT = 1
    ! Cell reps: We are fooling UppASD into thinging there is one cell
    !  as big as the whole system. Doing otherwise might cause
    !  some module to think there is a regular lattice spacing, or other
    !  unexpected pitfalls that may arise from having both atoms and continuum.
    N1 = 1
    N2 = 1
    N3 = 1

    ! Lattice vectors
    C1 = (/options%space%universeSize(1),0.0_dblprec,0.0_dblprec/)
    C2 = (/0.0_dblprec,options%space%universeSize(2),0.0_dblprec/)
    C3 = (/0.0_dblprec,0.0_dblprec,options%space%universeSize(3)/)
    
    ! Periodic boundaries
    BC1 = '0'
    BC2 = '0'
    BC3 = '0'
    if (options%space%periodicBoundary(1)) then
       BC1 = 'P'
    endif
    if (options%space%periodicBoundary(2)) then
       BC2 = 'P'
    endif
    if (options%space%periodicBoundary(3)) then
       BC3 = 'P'
    endif

    ! Warn about maptype
    if (maptype/=1) then
       print*, "WARNING: Maptype does not work in multiscale, use the 'unitcell' keyword" &
            // NEW_LINE('A') // &
            "          in coupling lists inside the multiscale configuration file."       
    end if
    
    ! No random alloy
    nchmax = 1
    if(do_ralloy /= 0) then
       print *, "WARNING: Ralloy will not work in multiscaled configurations."
    end if
    
    ! Prevent LSF
    if(conf_num /= 1) then !< Number of configurations for LSF
       print *, "WARNING: conf_num other than 1 is ignored in multiscale mode."
       conf_num = 1
    end if
    !!

    ! Required for anisotropy, deallocatesd after setup
    allocate(anumb_inp(NA),stat=i_stat)
    call memocc(i_stat,NA*kind(anumb_inp),'anumb_inp','setup_multiscale_system')
    do i=1,NA
       anumb_inp(i) = i
    end do

    if(do_prnstruct /= 0) then
       open(unit=1234, file=coord_file_name)
       do i = 1, Natom
          write (1234,'(i12,2x,3F19.13,2x,a,2x,i12 )') i, coord(:, i), '1', i
       enddo
       close(1234)
    end if
    
    call allocate_mmoms(Natom, Mensemble, 1)
    call allocate_emoms(Natom, Mensemble, 1)

    call allocate_multiscale(Natom, Mensemble,1)

    if (initmag==4) then
       write (*,'(2x,a)',advance='no') "Read from restart file"

         !if(multiscale_old_format  =='Y') then      
         !call loadrestart_old(Natom,Mensemble,restartfile,rstep,mmom,emom,emomM)
         !else
         call read_mag_conf(Natom,Mensemble,do_mom_legacy,rstep,restartfile,  &
                  mmom,emom,emomM)  
        ! end if


       write (*,*) " done"
       mmom2 = mmom
       emom2 = emom
       mmom0 = mmom
       mmomi = 1.0_dblprec/mmom0
    else if(initmag == 1) then
       print *,"WARNING: Don't use random moments with continuum regions."
       print *,"I AM SETTING EVERYTHING ALONG Z"
       do i = 1, natom
          alpha = 0.0d0
          beta  = pi/2.0d0
          emom(1,i,:) = sin(alpha) * cos(beta)
          emom(2,i,:) = sin(alpha) * sin(beta)
          emom(3,i,:) = cos(alpha)
          mmom(i,: ) = sqrt(sum(setup%moments(:,i)**2))
       end do
       mmom2 = mmom
       emom2 = emom
       mmom0 = mmom
       mmomi = 1.0_dblprec/mmom0

    else       
       call initializeMomentData(setup%moments,Natom,Mensemble)
    end if

    ! Required for anisotropy
    allocate(ammom_inp(NA,Nchmax,conf_num),stat=i_stat)
    call memocc(i_stat,product(shape(ammom_inp))*kind(ammom_inp),'ammom_inp', &
         'setup_multiscale_system')
    do j=1,conf_num
       do i=1,NA     
          ammom_inp(i,:,:) = 0
          ammom_inp(i,1,j) = mmom(i,1)
       end do
    end do


    allocate(atype(Natom_full),stat=i_stat)
    call memocc(i_stat,product(shape(atype))*kind(atype),'atype','setup_multiscale_system')
    allocate(anumb(Natom_full),stat=i_stat)
    call memocc(i_stat,product(shape(anumb))*kind(anumb),'anumb','setup_multiscale_system')
    allocate(Landeg(Natom_full),stat=i_stat)
    call memocc(i_stat,product(shape(Landeg))*kind(Landeg),'Landeg','setup_multiscale_system')
    !call zero_cumulant_counters()

    atype = 1
    Landeg = 1.0_dblprec

    ! Allocating the damping arrays
    call allocate_damping(Natom,ipnphase,0)
    write(*,'(1x,a)') ' Setup up damping'
    call setup_damping(NA,Natom,Natom_full,ipmode,mode,ipnphase,&
         do_ralloy,asite_ch,achem_ch,do_site_damping,&
         do_site_ip_damping,iplambda1,iplambda2,mplambda1,mplambda2,&
         iplambda1_array,iplambda2_array,lambda1_array,lambda2_array)


    if (ipmode /= 'MS' .and. ipmode /= 'N') then
       ! Why?:
       ! - Damping band is not implemented for other methods than Midpoint and Depondt
       ! - mc_iphase does not implement interpolation
     
       write (*, *) 'ERROR: Initial phase must be iether F or N with multiscale.' 
       stop 1
    endif

    ! 1 is Midpoint, 5 is Depondt
    if(SDEalgh /= 1 .and.  SDEalgh /= 5) then
       ! This condition might be relaxed; If the mode is spin dynamics
       !  everything but the damping band should work.
       write (*,*) 'ERROR: The only algorithms supported with multiscale are:'
       write (*,*) "  - SDEalgh 1  Mentik's Semi-implicit midpoint method."
       write (*,*) "  - SDEalgh 5  Depondt method."
       stop 1
    end if

    lambda1_array = mplambda1
    lambda2_array = mplambda2

    ! Allocating the temperature arrays
    call allocate_temp(Natom,ipnphase)
    write(*,'(1x,a)') ' Setup up temperature'
    if (ipnphase.ge.1) then
       do i=1, ipnphase
       
         

          call setup_temp(Natom,NT,NA,N1,N2,N3,Natom_full,do_ralloy,atype,acellnumb,&
          atype_ch,simid,iptemp(i),C1,C2,C3,BC1,BC2,BC3,Bas,coord,iptemp_array(:,i) )
           
          print*, ipTemp(i) 
       enddo
    endif
    

         call setup_temp(Natom,NT,NA,N1,N2,N3,Natom_full,do_ralloy,atype,acellnumb,&
          atype_ch,simid,temp,C1,C2,C3,BC1,BC2,BC3,Bas,coord,temp_array )
    write(*,'(1x,a)') ' done'
    !temp_array = temp

    ! Print input

    write (*,'(1x,a)') "Set up Hamiltonian"

    call setupExchanges(setup%atomsExchange, Natom,conf_num, &
         ncoupTmp,nlistsizeTmp,nlistTmp, exchange_file_name, do_prnmultiscale)
          
        
         ham%max_no_neigh = ubound(nlistTmp,1);

    
     call allocate_hamiltoniandata(Natom, 1, Natom,1,ubound(ncoupTmp, 1), 0, 'N', 1, 'N','N')

    ham%ncoup = ncoupTmp
    ham%nlistsize = nlistsizeTmp
    ham%nlist = nlistTmp
    !ham%aham=nlistsizeTmp
     do i = 1, natom
      ham%aham(i)=i
     end do
    isize = product(shape(ncoupTmp))*kind(ncoupTmp)
    deallocate(ncoupTmp)
    call memocc(i_stat, -isize, 'ncoupTmp', 'setup_multiscale_system')

    isize = product(shape(nlistsizeTmp))*kind(nlistsizeTmp)
    deallocate(nlistsizeTmp)
    call memocc(i_stat, -isize, 'nlistsizeTmp', 'setup_multiscale_system')

    isize = product(shape(nlistTmp))*kind(nlistTmp)
    deallocate(nlistTmp)
    call memocc(i_stat, -isize, 'nlistTmp', 'setup_multiscale_system')

    call rescale_loaded_exchangeCouplings()
    
    call setupDm(setup%atomsDm, Natom, &
         dmTmp,dmlistsizeTmp,dmlistTmp,ham_inp%do_dm, &
         dm_file_name, do_prnmultiscale)
    
    ham%max_no_dmneigh = ubound(dmlistTmp,1)
    if (ham_inp%do_dm /= 0) then
       
       call allocate_dmhamiltoniandata(Natom,Natom, ubound(dmTmp, 1), 1)
       
       ham%dm_vect = dmTmp
       ham%dmlistsize = dmlistsizeTmp
       ham%dmlist = dmlistTmp
       call rescale_loaded_dmCouplings()
    end if
      
    
    call multiscale_setup_anisotropy(setup)  

    ! Allocate arrays for simulation and measurement
    call allocate_general(1)

    write (*,'(1x,a)') "Initialize damping band interpolation "

    call setupDampingBand(setup%dampingBandWeights, &
         setup%dampingBandAttenuation, &
         natom, mensemble, dampingBand, &
         dband_file_name, do_prnmultiscale)
    call setupInterpolation(setup%interpolationWeights, &
         natom,interfaceInterpolation, &
         intp_file_name, do_prnmultiscale)
    
    
    ! Release temporary arrays
    isize = product(shape(anumb_inp))*kind(anumb_inp)
    deallocate(anumb_inp, stat=i_stat)
    call memocc(i_stat, -isize, 'anumb_inp', 'setup_multiscale_system')

    isize = product(shape(ammom_inp))*kind(ammom_inp)
    deallocate(ammom_inp, stat=i_stat)
    call memocc(i_stat, -isize, 'ammom_inp', 'setup_multiscale_system')

    call setupGradients(setup%gradientLinks, grad_link_file_name, do_prnmultiscale)
    
    if (ipmode=='MS'.and.do_site_ip_damping=='Y') then
       call read_ip_damping()
    end if

    call deallocateSetup(setup)

  end subroutine setup_multiscale_system

    subroutine allocate_general(flag)
    !
    use LLGI,          only : allocate_llgifields
    !use Depondt,       only : allocate_depondtfields
    use InputData
    use FieldData,     only : allocate_fields, read_local_field, allocation_field_time
    use Measurements,  only : allocate_measurementdata
    use RandomNumbers, only : allocate_randomwork
    use Midpoint_ms,   only : allocate_midpointms_fields
    use Depondt_ms,    only : allocate_depondtms_fields

    !
    implicit none
    !
    integer, intent(in) :: flag !< Allocate or deallocate (1/-1)

    !
    call allocate_fields(Natom,Mensemble,flag)

    if(locfield=='Y'.and.flag>0)  call read_local_field(NA,locfieldfile)
    if(SDEalgh==5) then
       call allocate_depondtms_fields(flag,Natom, Mensemble)
    elseif(SDEalgh==11) then
       call allocate_llgifields(Natom, Mensemble,flag)
    end if

    if (SDEalgh==1 .or. ipSDEalgh==1) then
      call allocate_midpointms_fields(flag,Natom,Mensemble)
    endif


    if(SDEalgh>=1.and.SDEalgh<=5) call allocate_randomwork(Natom,Mensemble,flag,'N')
 

    call allocate_measurementdata(NA,NT,Natom,Mensemble,Nchmax,plotenergy,flag)
  end subroutine allocate_general

  
  subroutine multiscale_setup_anisotropy(setup)
    use InputHandler_ext, only : read_anisotropy
    use HamiltonianInit      !only : setup_hamiltonian, setup_anisotropies
    use Multiscale
    use InputData
    use HamiltonianData 
    use ChemicalData, only : achem_ch
    implicit none
    type(MultiscaleSetup), intent(in) :: setup

    integer :: i_stat, isize

    if (do_ralloy /= 0) then
       print *, "WARNING: do_ralloy is ignored with multiscale 1."
       do_ralloy = 0
    end if
    if (ham_inp%random_anisotropy) then
       print *, "WARNING: random_anisotropy is ignored with multiscale 1."     
       ham_inp%random_anisotropy = .false.
    end if
    ! ToDo : test multi-axial
    if (ham_inp%mult_axis/='N') then
       print *, "WARNING: Does multi-axial anisotropy work?."     
    end if


    call allocate_anisotropies(Natom, ham_inp%mult_axis, 1)

    if(ham_inp%do_anisotropy==1) then
       allocate(ham_inp%anisotropytype(NA,Nchmax),stat=i_stat)
       call memocc(i_stat,product(shape(ham_inp%anisotropytype))*kind(ham_inp%anisotropytype),'anisotropytype','multiscale_setup_anisotropy')
       ham_inp%anisotropytype=0

       allocate(ham_inp%anisotropy(NA,6,Nchmax),stat=i_stat)
       call memocc(i_stat,product(shape(ham_inp%anisotropy))*kind(ham_inp%anisotropy),'anisotropy','multiscale_setup_anisotropy')
       ham_inp%anisotropy=0.0d0

       if (ham_inp%mult_axis=='Y') then
          allocate(ham_inp%anisotropytype_diff(NA,Nchmax),stat=i_stat)
          call memocc(i_stat,product(shape(ham_inp%anisotropytype_diff))*kind(ham_inp%anisotropytype_diff),'anisotropytype_diff','multiscale_setup_anisotropy')
          allocate(ham_inp%anisotropy_diff(NA,6,Nchmax),stat=i_stat)
          call memocc(i_stat,product(shape(ham_inp%anisotropy_diff))*kind(ham_inp%anisotropy_diff),'anisotropy_diff','multiscale_setup_anisotropy')
          ham_inp%anisotropy_diff=0.0d0
       endif
    end if

    ham%taniso = 0;
    ham%eaniso = 0;
    ham%kaniso = 0;
    ham%sb = 0;
    if(ham_inp%mult_axis == 'Y') then
       ham%taniso_diff = 0;
       ham%eaniso_diff = 0;
       ham%kaniso_diff = 0;
       ham%sb_diff = 0;     
    end if

    if(ham_inp%do_anisotropy == 1) then
       
        call setupAnisotropy(setup%anisotropies,&
            ham_inp%anisotropytype,ham_inp%anisotropy, &
           ham_inp%anisotropytype_diff, ham_inp%anisotropy_diff)

       call setup_anisotropies(Natom,NA,anumb_inp,ham_inp%anisotropytype,ham%taniso,          &
            ham%eaniso,ham%kaniso,ham%sb,ham_inp%anisotropy,ham_inp%mult_axis,ham_inp%anisotropytype_diff,  &
            ham%taniso_diff,ham%eaniso_diff,ham%kaniso_diff,ham%sb_diff,            &
            ham_inp%anisotropy_diff,ham_inp%random_anisotropy,0.0d0,do_ralloy,  &
            Natom_full,Nchmax,achem_ch,ammom_inp,0)


        isize = product(shape(ham_inp%anisotropytype))*kind(ham_inp%anisotropytype)
        deallocate(ham_inp%anisotropytype,stat=i_stat)
        call memocc(i_stat,-isize,'anisotropytype','multiscale_setup_anisotropy')

       isize = product(shape(ham_inp%anisotropy))*kind(ham_inp%anisotropy)
       deallocate(ham_inp%anisotropy,stat=i_stat)
       call memocc(i_stat,-isize,'anisotropy','multiscale_setup_anisotropy')

       if (ham_inp%mult_axis=='Y') then
          isize = product(shape(ham_inp%anisotropytype_diff))*kind(ham_inp%anisotropytype_diff)
          deallocate(ham_inp%anisotropytype_diff,stat=i_stat)
          call memocc(i_stat,-isize,'anisotropytype_diff','multiscale_setup_anisotropy')

          isize = product(shape(ham_inp%anisotropy_diff))*kind(ham_inp%anisotropy_diff)
          deallocate(ham_inp%anisotropy_diff,stat=i_stat)
          call memocc(i_stat,-isize,'anisotropy_diff','multiscale_setup_anisotropy')
       endif
    end if


  end subroutine multiscale_setup_anisotropy

  subroutine rescale_loaded_exchangeCouplings()
    
    use InputData, only : Natom
    use HamiltonianData, only : ham
    use MomentData, only : mmom
  implicit none
    integer :: i, j
    real(dblprec), dimension(1) :: tmpNcoup

    do i = 1, Natom
       do j = 1,ham%nlistsize(i)
          tmpNcoup(1) = ham%ncoup(j, i, 1)
          tmpNcoup = rescale_inputed_exchange(tmpNcoup, mmom(i, 1), &
                                              mmom(ham%nlist(j, i), 1), 1)
          ham%ncoup(j, i, 1) = tmpNcoup(1)
       end do
    end do
  end subroutine rescale_loaded_exchangeCouplings

  !> Rescales Dzyaloshinsky-Moriya interaction vectors
  !! that are already loaded into dm_vect
  subroutine rescale_loaded_dmCouplings()
    
    use InputData, only : Natom
    use HamiltonianData, only :ham
    use MomentData, only : mmom
  implicit none
    integer :: i, j
    real(dblprec), dimension(3) :: tmpDmCoup

    do i = 1, Natom
       do j = 1,ham%dmlistsize(i)
          tmpDmCoup = ham%dm_vect(:, j, i)
          tmpDmCoup = rescale_inputed_exchange(tmpDmCoup, mmom(i, 1), &
                                               mmom(ham%dmlist(j, i), 1), 1)
          ham%dm_vect(:, j, i) = tmpDmCoup
       end do
    end do
  end subroutine rescale_loaded_dmCouplings


  !> Adapt anisotropy setup from libmuasd to UppASD's representation
  subroutine setupAnisotropy(anisotropies, &
       anisotropytype,anisotropyvalues, &
       axis2_anisotropytype,axis2_anisotropyvalues)
    use Multiscale
    implicit none      
    type(AtomAnisotropies), intent(in) :: anisotropies
     integer, dimension(:,:), intent(inout) :: anisotropyType
    real(dblprec), dimension(:,:,:), intent(inout) :: anisotropyvalues
    integer, dimension(:,:), intent(inout) :: axis2_anisotropyType
    real(dblprec), dimension(:,:,:), intent(inout) :: axis2_anisotropyvalues

    integer, parameter :: VALUE_K1 = 1
    integer, parameter :: VALUE_K2 = 2
    integer, parameter :: VALUE_E_X = 3
    integer, parameter :: VALUE_E_Y = 4
    integer, parameter :: VALUE_E_Z = 5
    integer, parameter :: VALUE_RATIO = 6


    anisotropyType(:,1) = anisotropies%anisotropyTypes(:,1)

    anisotropyValues(:,VALUE_K1,1) = anisotropies%anisotropyKs(1,:,1)
    anisotropyValues(:,VALUE_K2,1) = anisotropies%anisotropyKs(2,:,1)    
    anisotropyValues(:,VALUE_E_X,1) = anisotropies%anisotropyE(1,:,1)
    anisotropyValues(:,VALUE_E_Y,1) = anisotropies%anisotropyE(2,:,1)
    anisotropyValues(:,VALUE_E_Z,1) = anisotropies%anisotropyE(3,:,1)
    anisotropyValues(:,VALUE_RATIO,1) = anisotropies%anisotropyRatios(:,1)

    if(ubound(anisotropies%anisotropyKs,3) > 1) then
       axis2_anisotropyType(:,1) = anisotropies%anisotropyTypes(:,2)

       axis2_anisotropyValues(:,VALUE_K1,1) = anisotropies%anisotropyKs(1,:,2)
       axis2_anisotropyValues(:,VALUE_K2,1) = anisotropies%anisotropyKs(2,:,2)    
       axis2_anisotropyValues(:,VALUE_E_X,1) = anisotropies%anisotropyE(1,:,2)
       axis2_anisotropyValues(:,VALUE_E_Y,1) = anisotropies%anisotropyE(2,:,2)
       axis2_anisotropyValues(:,VALUE_E_Z,1) = anisotropies%anisotropyE(3,:,2)
       axis2_anisotropyValues(:,VALUE_RATIO,1) = anisotropies%anisotropyRatios(:,2)
    end if
  end subroutine setupAnisotropy

  !> Converts interpolation data from MUASD to UppASD
  !! @param[in] weights Interpolation weights in as a MUASD matrix. MUST be sorted by row.
  !! @param[in] natoms  Atoms in the system.
  !! @param[in,out] interpolation info in uppasd.
  !! @param[in] dump_name dump the interpolation data to a file.
  !! @param[in] dump if true, data will be written in the file dump_name
  subroutine setupInterpolation(weights, natoms, interpolation, dump_name, dump)
    use SparseMatrix
    use MultiscaleInterpolation
    use Profiling
    implicit none
    type(SpMatrix), intent(in) :: weights
    integer, intent(in) :: natoms
    type(LocalInterpolationInfo), intent(inout) :: interpolation
    character(len=OUTPUT_FILE_LEN), intent(in) :: dump_name
    logical, intent(in) :: dump

    integer :: i,j,pre_row, nrWeights, nRows, atom, row
    integer :: i_stat
    
    nrWeights = weights%row%length

    if (nrWeights /= 0) then
       
       nRows = 1
       do i = 2, nrWeights          
          if(weights%row%values(i-1) /= weights%row%values(i)) then
             nRows = nRows + 1
          end if
       end do

       allocate(interpolation%weights(nrWeights),stat=i_stat)
       call memocc(i_stat, &
            product(shape(interpolation%weights))*kind(interpolation%weights), &
            'LocalInterpolationInfo%weights','setupInterpolation')
       allocate(interpolation%firstNeighbour(nRows+1),stat=i_stat)
       call memocc(i_stat, &
            product(shape(interpolation%firstNeighbour))* &
            kind(interpolation%firstNeighbour), &
            'LocalInterpolationInfo%firstNeighbour','setupInterpolation')
       allocate(interpolation%neighbours(nrWeights),stat=i_stat)
       call memocc(i_stat, &
            product(shape(interpolation%neighbours))* &
            kind(interpolation%neighbours), &
            'LocalInterpolationInfo%neighbours','setupInterpolation')
       allocate(interpolation%indices(Natoms),stat=i_stat)
       call memocc(i_stat, &
            product(shape(interpolation%indices))*kind(interpolation%indices), &
            'LocalInterpolationInfo%indices','setupInterpolation')

       interpolation%weights    = weights%entries%values(1:nrWeights)
       interpolation%neighbours = weights%col%values(1:nrWeights)
       interpolation%indices = 0
       interpolation%firstNeighbour = 0
       interpolation%firstNeighbour(1) = 1

       pre_row = 1
       row = 1
       do i = 1, nrWeights
          atom = weights%row%values(i)
          if(interpolation%indices(atom) == 0) then             
             interpolation%indices(atom) = row
             do j=pre_row+1,row
                interpolation%firstNeighbour(j) = i
             end do
             pre_row = row
             row=row+1
          end if
       end do
       interpolation%firstNeighbour(nRows+1) = nrWeights+1       
       interpolation%nrInterpAtoms = nRows

       ! Dump to file
       if(dump) then
          open(unit=661, file=trim(dump_name))
          do atom = 1, natoms
             i = interpolation%indices(atom)
             if(i/=0) then
                do j = interpolation%firstNeighbour(i),&
                     interpolation%firstNeighbour(i+1)-1
                   write(661,*) atom, &
                        interpolation%neighbours(j), &
                        interpolation%weights(j)
                end do
             end if
          end do
          close(661)
       end if
    end if
    
  end subroutine setupInterpolation

  !> Converts a damping band interpolation from MUASD to UppASD representation.
  !! @param[in] muasd_weights Weights in as a MUASD matrix. MUST be sorted by row.
  !! @param[in] muasd_positional Positional coefficients from Muasd. The band strength already comes multiplied to them.
  !! @param[in] natoms  Atoms in the system.
  !! @param[in] mensemble  Number of ensembles, required to allocate the interpolation buffer.
  !! @param[in,out] dampingBand info in uppasd.
  !! @param[in] dump_name dump the interpolation data to a file.
  !! @param[in] dump dump the interpolation data to a file.
  subroutine setupDampingBand(muasd_weights, muasd_positional, &
       natom, mensemble, dampingBand, dump_file, dump)
    use SparseMatrix
    use MultiscaleDampingBand, only : DampingBandData
    implicit none
    type(SpMatrix), intent(in) :: muasd_weights
    real(dblprec),dimension(:), intent(in) :: muasd_positional
    integer, intent(in) :: natom, mensemble
    type(DampingBandData),intent(inout) :: dampingBand
    character(len=OUTPUT_FILE_LEN), intent(in) :: dump_file
    logical, intent(in) :: dump
    
    integer :: i_stat

    if(size(muasd_positional) > 0 .and. &
       muasd_weights%row%length > 0) then
       
       call setupInterpolation(muasd_weights,natom,&
            dampingBand%interpolation,dump_file, dump)
       dampingBand%enable = associated(dampingBand%interpolation%indices)
       allocate(dampingBand%coefficients(size(muasd_positional)),stat=i_stat)
       call memocc(i_stat,&
            product(shape(dampingBand%coefficients))*&
            kind(dampingBand%coefficients),&
            'DampingBandData%coefficients','setupDampingBand')       
       dampingBand%coefficients = muasd_positional
       allocate(dampingBand%preinterpolation(3,&
            dampingBand%interpolation%nrInterpAtoms,mensemble),&
            stat=i_stat)
       call memocc(i_stat,&
            product(shape(dampingBand%preinterpolation))*&
            kind(dampingBand%preinterpolation),&
            'DampingBandData%preinterpolation','setupDampingBand')
       dampingBand%enable = .true.

       ! Check for NaNs
       if (any(dampingBand%coefficients /= dampingBand%coefficients)) then
          print *,"NaN's found in dband coeffs"
          stop
       end if
       if (any(dampingBand%interpolation%weights /= dampingBand%interpolation%weights)) then
          print *,"NaN's found in dband weights. Is your window large enough?"
          stop
       end if
       
    else
       dampingBand%enable = .false.
    end if

  end subroutine setupDampingBand



  !> Converts exchanges from a MUASD sparse matrix to UppASD
  !!  internal representation
  !! @param[in] exchanges Exchanges in MUASD
  !! @param[in] natoms Number of atoms in the system
  !! @param[in] nconf Number of configurations, for LSF
  !! @param[in, out] ncoupTmp Values of the exchange, might be reallocated
  !! @param[in, out] nlistsizeTmp Number of entries per column, might be reallocated
  !! @param[in, out] nlistTmp Column index per element, might be reallocated
  !! @param[in] dump_file file to dump the exhange values.
  !! @param[in] dump dump the exhange values.
  !! ToDo: Profile, maybe use a common, more efficient format
  subroutine setupExchanges(exchange,natoms,nconf, &
       ncoupTmp,nlistsizeTmp,nlistTmp, dump_file, dump)
    use SparseMatrix
    implicit none
    type(SpMatrix), intent(in) :: exchange
    integer, intent(in) :: natoms
    integer, intent(in) :: nconf
    real(dblprec),dimension(:,:,:),allocatable,intent(inout) :: ncoupTmp
    integer, dimension(:), allocatable,intent(inout) :: nlistsizeTmp
    integer, dimension(:, :), allocatable, intent(inout) :: nlistTmp
    character(len=OUTPUT_FILE_LEN), intent(in) :: dump_file
    logical, intent(in) :: dump
    
    integer :: max_row_elems, nrows
    integer :: row, col
    integer :: columns, last_row
    integer :: ent, i, j
    real(dblprec) :: val

    nrows = natoms

    max_row_elems = 0
    ! Count row elems
    columns = 1
    do i=2,exchange%row%length
       if (exchange%row%values(i-1) == exchange%row%values(i))then
          columns = columns + 1
       else
          max_row_elems = max(max_row_elems,columns)
          columns = 1
       end if
    end do
    max_row_elems = max(max_row_elems,columns)

    ! Resize the matrix when needed
    call ensureNcoupSize(nrows, max_row_elems, nconf, ncoupTmp,nlistsizeTmp,nlistTmp)

    ! Write the elements
    ncoupTmp = 0d0 !
    nlistTmp = 0
    nlistsizeTmp = 0 ! Reset the number of elements per row, as the matrix might be larger than needed now.
    last_row = 0
    ent = 0
    do i = 1, exchange%row%length
       row = exchange%row%values(i)
       col = exchange%col%values(i)
       val = exchange%entries%values(i)
       if (row .ne. last_row) then
          ent = 0
          last_row = row
       end if
       ent = ent + 1
       nlistsizeTmp(row)   = ent
       ncoupTmp(ent,row,:) = val
       nlistTmp(ent,row)   = col
    end do

    if (dump) then 
       open(unit=112,file=trim(adjustl(dump_file)))
       do i = 1, ubound(nlistsizeTmp,1)
          do j=1,nlistsizeTmp(i)
             write(112,*) i,nlistTmp(j,i),ncoupTmp(j,i,1)
          end do
       end do
       close(112)
    end if
  end subroutine setupExchanges


  
  !> Converts dms from a MUASD sparse matrix to UppASD
  !!  internal representation
  !! @param[in] dm Dms in MUASD
  !! @param[in] natoms Number of atoms in the system
  !! @param[in] nconf Number of configurations, for LSF
  !! @param[in, out] dmTmp Values of the dm, might be reallocated
  !! @param[in, out] dmlistsize Number of entries per column, might be reallocated
  !! @param[in, out] dmlist Column index per element, might be reallocated
  !! @param[in] dump_file file to dump the DM vectors.
  !! @param[in] dump 
  !! ToDo: Profile, maybe use a common, more efficient format
  subroutine setupDm(dm,natoms, &
       dmTmp,dmlistsizeTmp,dmlistTmp,do_dm, &
       dump_file, dump)
    use SparseMatrix
    implicit none
    type(SpMatrix), intent(in) :: dm
    integer, intent(in), optional :: natoms
    real(dblprec),dimension(:,:,:),allocatable,intent(inout) :: dmTmp
    integer, dimension(:), allocatable,intent(inout) :: dmlistsizeTmp
    integer, dimension(:, :), allocatable, intent(inout) :: dmlistTmp
    integer, intent(inout) :: do_dm
    character(len=OUTPUT_FILE_LEN), intent(in) :: dump_file
    logical, intent(in) :: dump
    
    
    integer :: max_row_elems, nrows
    integer :: row, col, i_col
    integer :: columns, last_row
    integer :: ent, i, j, k
    real(dblprec) :: val

    
    nrows = natoms
        
    max_row_elems = 0
    ! Count row elems
    columns = 1
    do i=2,dm%row%length
       if (dm%row%values(i-1) == dm%row%values(i))then
          columns = columns + 1
       else
          max_row_elems = max(max_row_elems,columns)
          columns = 1
       end if
    end do
    max_row_elems = (max(max_row_elems,columns)-1) / 3 + 1
    
    ! Resize the matrix when needed
    call ensureDmSize(nrows, max_row_elems, dmTmp,dmlistsizeTmp,dmlistTmp)
    
    ! Write the elements 
    dmTmp = 0d0 !
    dmlistTmp = 0
    dmlistsizeTmp = 0 ! Reset the number of elements per row, as the matrix might be larger than needed now.
    last_row = 0
    ent = 0
    do_dm = 0
    if (dm%row%length > 0) then
       do_dm = 1
    end if
    do i = 1, dm%row%length
       row = dm%row%values(i)
       i_col = dm%col%values(i)
       col = (i_col - 1) / 3 + 1
       k = modulo(i_col-1,3)+1 
       val = dm%entries%values(i)
       if (row .ne. last_row) then
          ent = 0
          last_row = row
       end if
       ent = ent + 1
       dmlistsizeTmp(row) = (ent-1)/3+1
       dmTmp(k,(ent-1)/3+1,row) = val
       dmlistTmp((ent-1)/3+1,row) = col
    end do

    if (dump) then
       open(unit=112,file=trim(adjustl(dump_file)))
       do i = 1, ubound(dmlistsizeTmp,1)
          do j=1,dmlistsizeTmp(i)
             write(112,*) i,dmlistTmp(j,i),dmTmp(:,j,i)
          end do
       end do
       close(112)
    end if
    
  end subroutine setupDm

  
  ! [Re]Allocate, only if needed, the exchange matrix
  subroutine ensureNcoupSize(nrows,max_row_elems,nconf, &
       ncoupTmp,nlistsizeTmp,nlistTmp)
    use Profiling
    implicit none
    integer, intent(in) :: nrows
    integer, intent(in) :: max_row_elems
    integer, intent(in) :: nconf

    real(dblprec),dimension(:,:,:),allocatable,intent(inout) :: ncoupTmp  
    integer, dimension(:), allocatable,intent(inout) :: nlistsizeTmp
    integer, dimension(:, :), allocatable, intent(inout) :: nlistTmp

    integer :: i_stat, isize
    ! ncoup, of size (max_row_elems, nrows, nconf)
    if(allocated(ncoupTmp)) then
       if (any(size(ncoupTmp) .le. (/ max_row_elems, nrows, nconf /))) then

          isize = product(shape(ncoupTmp))*kind(ncoupTmp)
          deallocate(ncoupTmp, stat=i_stat)
          call memocc(i_stat, -isize, 'ncoupTmp', 'ensureNcoupSize')

          allocate(ncoupTmp(max_row_elems,nrows,nconf), stat=i_stat)
          call memocc(i_stat,product(shape(ncoupTmp))*kind(ncoupTmp),&
               'ncoupTmp','ensureNcoupSize')

       end if
    else
       allocate(ncoupTmp(max_row_elems,nrows,nconf), stat=i_stat)
       call memocc(i_stat,product(shape(ncoupTmp))*kind(ncoupTmp),&
            'ncoupTmp','ensureNcoupSize')
       
    end if
    
    ! nlistsizeTmp (nrows)
    if(allocated(nlistsizeTmp)) then
       if (any(size(nlistsizeTmp) .le. (/ nrows /))) then

          isize = product(shape(nlistsizeTmp))*kind(nlistsizeTmp)
          deallocate(nlistsizeTmp, stat=i_stat)
          call memocc(i_stat, -isize, 'nlistsizeTmp', 'ensureNcoupSize')

          allocate(nlistsizeTmp(nrows), stat=i_stat)
          call memocc(i_stat,product(shape(nlistsizeTmp))*kind(nlistsizeTmp),&
               'nlistsizeTmp','ensureNcoupSize')

       end if
    else
       allocate(nlistsizeTmp(nrows), stat=i_stat)
       call memocc(i_stat,product(shape(nlistsizeTmp))*kind(nlistsizeTmp),&
            'nlistsizeTmp','ensureNcoupSize')

    end if

    ! nlistTmp(max_row_elems, nrows)
    if(allocated(nlistTmp)) then
       if (any(size(nlistTmp) .le. (/ max_row_elems,nrows /))) then

          isize = product(shape(nlistTmp))*kind(nlistTmp)
          deallocate(nlistTmp, stat=i_stat)
          call memocc(i_stat, -isize, 'nlistTmp', 'ensureNcoupSize')

          allocate(nlistTmp(max_row_elems,nrows), stat=i_stat)
          call memocc(i_stat,product(shape(nlistTmp))*kind(nlistTmp),&
               'nlistTmp','ensureNcoupSize')

       end if
    else
       allocate(nlistTmp(max_row_elems,nrows), stat=i_stat)
       call memocc(i_stat,product(shape(nlistTmp))*kind(nlistTmp),&
            'nlistTmp','ensureNcoupSize')

    end if
  end subroutine ensureNcoupSize

  ![Re]Allocate dm matrix, almost the same as ensureNcoupSize
  ! but the names for memocc are different and the size of dmTmp works
  ! in another way.
  subroutine ensureDmSize(nrows,max_row_elems, &
       dmTmp,dmlistsizeTmp,dmlistTmp)
    use Profiling
    implicit none
    integer, intent(in) :: nrows
    integer, intent(in) :: max_row_elems

    real(dblprec),dimension(:,:,:),allocatable,intent(inout) :: dmTmp  
    integer, dimension(:), allocatable,intent(inout) :: dmlistsizeTmp
    integer, dimension(:, :), allocatable, intent(inout) :: dmlistTmp

    integer :: i_stat, isize

    ! dm, of size (max_row_elems, nrows, nconf)
    if(allocated(dmTmp)) then
       if (any(size(dmTmp) .le. (/ 3, max_row_elems, nrows /))) then

          isize = product(shape(dmTmp))*kind(dmTmp)
          deallocate(dmTmp, stat=i_stat)
          call memocc(i_stat, -isize, 'dmTmp', 'ensureDmSize')

          allocate(dmTmp(3, max_row_elems,nrows), stat=i_stat)
          call memocc(i_stat,product(shape(dmTmp))*kind(dmTmp),&
               'dmTmp','ensureDmSize')

       end if
    else
       allocate(dmTmp(3,max_row_elems,nrows), stat=i_stat)
       call memocc(i_stat,product(shape(dmTmp))*kind(dmTmp),&
            'dmTmp','ensureDmSize')

    end if

    ! dmlistsizeTmp (nrows)
    if(allocated(dmlistsizeTmp)) then
       if (any(size(dmlistsizeTmp) .le. (/ nrows /))) then          
          isize = product(shape(dmlistsizeTmp))*kind(dmlistsizeTmp)
          deallocate(dmlistsizeTmp, stat=i_stat)
          call memocc(i_stat, -isize, 'dmlistsizeTmp', 'ensureDmSize')

          allocate(dmlistsizeTmp(nrows), stat=i_stat)
          call memocc(i_stat,product(shape(dmlistsizeTmp))*kind(dmlistsizeTmp),&
               'dmlistsizeTmp','ensureDmSize')

       end if
    else
       allocate(dmlistsizeTmp(nrows), stat=i_stat)
       call memocc(i_stat,product(shape(dmlistsizeTmp))*kind(dmlistsizeTmp),&
            'dmlistsizeTmp','ensureDmSize')
    end if

    ! dmlistTmp(max_row_elems, nrows)
    if(allocated(dmlistTmp)) then
       if (any(size(dmlistTmp) .le. (/ max_row_elems,nrows /))) then

          isize = product(shape(dmlistTmp))*kind(dmlistTmp)
          deallocate(dmlistTmp, stat=i_stat)
          call memocc(i_stat, -isize, 'dmlistTmp', 'ensureDmSize')

          allocate(dmlistTmp(max_row_elems,nrows), stat=i_stat)
          call memocc(i_stat,product(shape(dmlistTmp))*kind(dmlistTmp),&
               'dmlistTmp','ensureDmSize')

       end if
    else
       allocate(dmlistTmp(max_row_elems,nrows), stat=i_stat)
       call memocc(i_stat,product(shape(dmlistTmp))*kind(dmlistTmp),&
            'dmlistTmp','ensureDmSize')

    end if
  end subroutine ensureDmSize
  
  subroutine setupGradients (links, dump_file, dump)   
    use InputData
    use SystemData, only : coord
    use HamiltonianData, only : ham       !only : max_no_neigh, nlistsize, ncoup, nlist
    use InputHandler !only: read_jvecfile
    use Prn_Topology, only: skyno, do_proj_skyno
    use Gradients
    use Spintorques
    use MultiscaleGradients
    use SparseMatrix, only : SpMatrix
  implicit none
    type(SpMatrix), intent(in) :: links
    character(len=OUTPUT_FILE_LEN), intent(in) :: dump_file
    logical, intent(in) :: dump

    integer :: i, row, col, i_stat
    real(dblprec) :: val

    ! Dump if needed
    if (dump) then 
       open(unit=112,file=trim(adjustl(dump_file)))
       do i = 1, links%row%length
          row = links%row%values(i)
          col = links%col%values(i)
          val = links%entries%values(i)
          write(112,*) row,col,val
       end do
       close(112)
    end if

    if (associated(links%row%values) .and. links%row%length > 0) then
       stt = 'A'
    else
       if (stt /= 'N') then
          print *,"WARNING: STT option is ignored in Multiscale mode, please check"
          print *,"         stt_vector and stt_window_size options on multiscale conf."
       end if
    end if
    
    ! site-dependent j is not supported, sorry
    ! It would be nice if someone implemented region-based specification in
    ! multiscale, like for moments, anisotropy or zones.
    if (jsite/='N') then
       print *,"Site-dependent stt vector is still not supported on multiscale."
       stop 1
    end if

    ! if STT is used initialize gradient interpolation
    if (stt/='N') then
       allocate(sitenatomjvec(3,Natom),stat=i_stat)
       call memocc(i_stat,product(shape(sitenatomjvec))*kind(sitenatomjvec),'sitenatomjvec','setupGradients')
       
       call initializeMultiscaleGradients(links, Natom)
        
       do i=1, Natom
          sitenatomjvec(1,i)=1.0_dblprec
          sitenatomjvec(2,i)=0.0_dblprec
          sitenatomjvec(3,i)=0.0_dblprec
       enddo
       !call allocate_stt_data(Natom,Mensemble,stt,do_she,flag=1)
        call allocate_stt_data(Natom,Mensemble,flag=1)

    endif
    if (skyno=='Y'.or.do_proj_skyno=='Y') then
       print *,"WARNING: In multiscale mode, skyno and proj_skyno"
       print *,"         are not tested."
       call setup_stencil_mesh( &
            Natom, N1, N2, N3, C1, C2, C3, BC1, BC2, BC3, &
            ham%max_no_neigh, ham%nlistsize, ham%nlist, coord)
    end if

  end subroutine setupGradients

  ! Writes the atom regions by indices
  ! Instead of writing all indices, secuential indices are clumped together in ranges, for the shake of brevity.
  ! Ranges are inlusive, meaning the file
  ! # fully coarse grained atoms
  ! 6 6
  ! Specifies only one fully coarse-grained atom of index 6
  ! Regions are guaranteed to be sorted, as the output from 'runMultiscaleSetup'
  ! is sorted
  subroutine printMultiscaleRegions(regions, dump_file)
    implicit none
    type(MultiscaleRegions), intent(in) :: regions
    character(len=OUTPUT_FILE_LEN), intent(in) :: dump_file
    

    open(unit=112,file=trim(adjustl(dump_file)))
    call write_ranges(112,'# fully coarse grained atoms',regions%coarseGrained)
    call write_ranges(112,'# partially coarse grained atoms',regions%partiallyCoarseGrained)
    call write_ranges(112,'# real atoms',regions%realAtoms)
    call write_ranges(112,'# no damped atoms',regions%nonDampingAtoms)
    close(112)

  contains
    
    subroutine write_ranges(file, tag, indices) 
      implicit none
      integer, intent(in) :: file
      character(len=*), intent(in) :: tag
      integer, dimension(:) :: indices
      integer :: i, lower, upper, index
      
      if (ubound(indices,1) > 0) then         
         lower = indices(1)
         upper = lower
         write (file, *) trim(tag)
         do i=2,ubound(indices,1)
            index = indices(i)
            if (index == upper+1) then
               upper = index
            else
               write (file,*) lower,upper
               lower = index
               upper = lower
            end if            
         end do
         write (file,*) lower,upper

      end if

         
    end subroutine write_ranges
  end subroutine printMultiscaleRegions


   subroutine initializeMomentData(moments,Natom,Mensemble)
    
    use MomentData           
    implicit none
    real(dblprec), dimension(:,:), intent(in) :: moments
    integer, intent(in) :: Natom, Mensemble

    integer :: atom,ens,i
    real(dblprec) :: mag

    do ens = 1,Mensemble
       emomM(:, :, ens) = moments
       do atom = 1,Natom
          mag = sqrt(sum(moments(:, atom)**2))
          mmom(atom, ens) = mag
          mmom2(atom, ens) = mag
          emom(1:3, atom, ens) = emomM(:, atom, ens)/mag
          emom2(1:3, atom, ens) = emom(1:3, atom, ens)
          do i = 1,ubound(multiscaleBackbuffer,4)
             multiscaleBackbuffer(1:3, atom, ens, i) = emom(1:3, atom, ens)
          end do
       end do
    end do
    mmom0 = mmom
    mmomi = 1.0_dblprec/mmom0
    
  end subroutine initializeMomentData


!> (De)Allocate backbuffer for multiscale
  subroutine allocate_multiscale(Natom,Mensemble,flag)
    implicit none
    integer, intent(in), optional :: Natom !< Number of atoms in system
    integer, intent(in), optional :: Mensemble !< Number of ensembles 
    integer, intent(in) :: flag  !< Allocate or deallocate (1/-1)
    integer :: i_all,i_stat
    if(flag>0) then
       allocate(multiscaleBackbuffer(3,Natom,Mensemble,3),stat=i_stat)
       call memocc(i_stat,&
            product(shape(multiscaleBackbuffer))*kind(multiscaleBackbuffer),&
            'multiscaleBackbuffer','allocate_multiscale')
       multiscaleBackbufferHead = 1
    else
       i_all=-product(shape(multiscaleBackbuffer))*kind(multiscaleBackbuffer)
       deallocate(multiscaleBackbuffer,stat=i_stat)
       call memocc(i_stat,i_all,'multiscaleBackbuffer','allocate_multiscale')       
    end if
    
  end subroutine allocate_multiscale


  
end module MultiscaleSetupSystem
