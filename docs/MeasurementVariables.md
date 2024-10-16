# Measurement variables
A listing of variables used for measurement of physical observables in UppASD.

## Magnetization and cumulants
#### source/Measurement/prn_averages.f90

   ! Printing definitions
   integer :: cumu_step !< Interval for sampling Binder cumulant
   integer :: cumu_buff !< Buffer size for Binder cumulant
   integer :: avrg_step !< Interval for sampling average magnetization
   integer :: avrg_buff !< Buffer size for average magnetization
   character(len=1) :: do_avrg                  !< Measure average magnetization (Y/N)
   character(len=1) :: do_cumu                  !< Measure Binder cumulant, susceptibility, and specific heat(Y/N)
   character(len=1) :: do_cuda_avrg             !< Measure average magnetization (Y/N) with CUDA
   character(len=1) :: do_cuda_cumu             !< Measure Binder cumulant, susceptibility, and specific heat(Y/N) with CUDA
   character(len=1) :: do_proj_avrg             !< Measure projected averages (Y/A/N)
   character(len=1) :: do_projch_avrg           !< Measure chemically projected averages (Y/N)
   character(len=1) :: do_cumu_proj             !< Measure Binder cumulant, susceptibility, and specific heat(Y/N)

   ! Local calculations for printing
   integer :: Navrgcum        !< Counter for number of cumulated averages
   real(dblprec) :: cumuw     !< Weight for current sample to cumulant
   real(dblprec) :: cumutotw  !< Sum of all cumulant weights
   integer :: Navrgcum_proj   !< Counter for number of cumulated projected averages
   real(dblprec) :: mavg      !< Average magnetic moment
   real(dblprec) :: totene    !< Total energy
   real(dblprec) :: binderc   !< Binder cumulant
   real(dblprec) :: avrglcum  !< Cumulated average of l
   real(dblprec) :: avrgmcum  !< Cumulated average of m
   real(dblprec) :: avrgecum  !< Cumulated average of E
   real(dblprec) :: avrgetcum !< Cumulated average of E_xc
   real(dblprec) :: avrgelcum !< Cumulated average of E_LSF
   real(dblprec) :: avrgm4cum !< Cumulated average of m^4
   real(dblprec) :: avrgm2cum !< Cumulated average of m^2
   real(dblprec) :: avrgl4cum !< Cumulated average of l^4
   real(dblprec) :: avrgl2cum !< Cumulated average of l^2
   real(dblprec) :: avrge2cum !< Cumulated average of E^2

   ! Definition of magnetization and cumulant arrays
   real(dblprec), dimension(:), allocatable       :: indxb_avrg       !< Step counter for average magnetization
   real(dblprec), dimension(:,:,:), allocatable   :: mavg_buff        !< Buffer for average magnetizations
   real(dblprec), dimension(:,:,:), allocatable   :: mavg2_buff_proj  !< Buffer for squared projected averages
   real(dblprec), dimension(:,:,:,:), allocatable :: mavg_buff_proj   !< Buffer for projected averages
   real(dblprec), dimension(:,:,:,:), allocatable :: mavg_buff_projch !< Buffer for chemical projected averages
   real(dblprec), dimension(:), allocatable :: avrgm4cum_proj !< Cumulated average of projected m^4
   real(dblprec), dimension(:), allocatable :: avrgm2cum_proj !< Cumulated average of projected m^2
   real(dblprec), dimension(:), allocatable :: avrgmcum_proj  !< Cumulated average of projected m
   integer :: bcount_avrg    !< Counter of buffer for averages

   ! Allocation of magnetization and cumulant arrays
   allocate(mavg_buff(3,avrg_buff,Mensemble),stat=i_stat)
   allocate(mavg_buff_proj(3,NA,avrg_buff,Mensemble),stat=i_stat)
   allocate(mavg2_buff_proj(NA,avrg_buff,Mensemble),stat=i_stat)
   allocate(indxb_avrg(avrg_buff),stat=i_stat)
   allocate(mavg_buff_projch(4,Nchmax,avrg_buff,Mensemble),stat=i_stat)
   allocate(avrgmcum_proj(NT),stat=i_stat)
   allocate(avrgm2cum_proj(NT),stat=i_stat)
   allocate(avrgm4cum_proj(NT),stat=i_stat)



## Trajectories
#### source/Measurements/prn_trajectories.f90

   ! Input parameters
   integer :: ntraj         !< Number of trajectories to sample
   integer :: tottraj_step  !< Interval for sampling magnetic moments
   integer :: tottraj_buff  !< Buffer size for magnetic moments
   character(len=1) :: do_tottraj !< Measure magnetic moments

   ! Definition of trajectories arrays
   integer, dimension(:), allocatable :: traj_step !< Interval for sampling individual trajectories
   integer, dimension(:), allocatable :: traj_buff !< Buffer size for individual trajectories
   integer, dimension(:), allocatable :: traj_atom !< List of atoms to sample trajectories for
   integer, dimension(:), allocatable :: bcount_traj   !< Counter of buffer for trajectories
   real(dblprec), dimension(:), allocatable :: indxb            !< Step counter for magnetic moments
   real(dblprec), dimension(:), allocatable :: scount_traj      !< Counter of sampling for trajectories
   real(dblprec), dimension(:,:), allocatable :: indxb_traj     !< Step counter for individual trajectories
   real(dblprec), dimension(:,:,:), allocatable :: mmomb        !< Buffer for all moment magnitudes
   real(dblprec), dimension(:,:,:), allocatable :: mmomb_traj   !< Buffer for selected moment magnitudes
   real(dblprec), dimension(:,:,:,:), allocatable :: emomb      !< Buffer for all individual trajectories
   real(dblprec), dimension(:,:,:,:), allocatable :: emomb_traj !< Buffer for selected individual trajectories

   ! Allocation of trajectories arrays
   allocate(scount_traj(ntraj),stat=i_stat)
   allocate(bcount_traj(ntraj),stat=i_stat)
   allocate(emomb_traj(3,maxval(traj_buff),ntraj,Mensemble),stat=i_stat)
   allocate(mmomb_traj(maxval(traj_buff),ntraj,Mensemble),stat=i_stat)
   allocate(indxb_traj(maxval(traj_buff),ntraj),stat=i_stat)
   allocate(mmomb(Natom,tottraj_buff,Mensemble),stat=i_stat)
   allocate(emomb(3,Natom,tottraj_buff,Mensemble),stat=i_stat)
   allocate(indxb(tottraj_buff),stat=i_stat)
   allocate(traj_atom(ntraj),stat=i_stat)
   allocate(traj_step(ntraj),stat=i_stat)
   allocate(traj_buff(ntraj),stat=i_stat)



## Autocorrelation
#### source/Correlation/autocorrelation.f90

   !From input
   integer :: nspinwait
   integer, dimension(:), allocatable :: spinwaitt
   character(len=35) :: twfile                  !< File name for autocorrelation waiting times
   character(len=1) :: do_autocorr              !< Perform autocorrelation (Y/N)

   ! Printing definitions
   integer :: ac_step !< Interval for sampling average magnetization
   integer :: ac_buff !< Buffer size for average magnetization

   ! Definition of autocorrelation arrays
   integer, dimension(:), allocatable :: spinwaitt
   real(dblprec), dimension(:,:,:), allocatable :: spinwait !< Data for autocorrelation analysis
   real(dblprec), dimension(:,:,:), allocatable   :: autocorr_buff        !< Buffer for average magnetizations
   real(dblprec), dimension(:), allocatable       :: indxb_ac       !< Step counter forautocorrelation
   real(dblprec), dimension(3,Mensemble) ::  m
   real(dblprec), dimension(:), allocatable :: autocorr
   real(dblprec), dimension(:), allocatable :: autocorr ! Note, in another subroutine

   ! Allocation of autocorrelation arrays
   allocate(spinwait(3,Natom,nspinwait),stat=i_stat)
   allocate(autocorr_buff(nspinwait,ac_buff,Mensemble),stat=i_stat)
   allocate(indxb_ac(ac_buff),stat=i_stat)
   allocate(autocorr(nspinwait),stat=i_stat)
   allocate(autocorr(nspinwait),stat=i_stat) ! Note, in another subroutine
   allocate(spinwaitt(nspinwait),stat=i_stat)



## Polarization
#### source/Measurement/polarization.f90

   ! Definition of polarization arrays
   real(dblprec), dimension(:,:,:), allocatable :: eij            !< Normalized neighbour distances
   real(dblprec), dimension(:,:,:), allocatable :: local_pol      !< Data for local polarization vector
   real(dblprec), dimension(:,:,:), allocatable :: s_cross_s      !< Data for local polarization vector
   integer, dimension(:), allocatable :: pollistsize              !< Size of neighbour list for polarization calculations
   real(dblprec), dimension(:), allocatable :: indxb_pol          !< Step counter for polarization vector
   real(dblprec), dimension(:,:,:), allocatable :: spol_buff      !< Buffer data for polarization vector

   ! Allocation of polarization arrays
   allocate(eij(3,max_no_neigh,Natom),stat=i_stat)
   allocate(local_pol(3,Natom,Mensemble),stat=i_stat)
   allocate(s_cross_s(3,Natom,Mensemble),stat=i_stat)
   allocate(pollistsize(Natom),stat=i_stat)
   allocate(spol_buff(3,pol_buff,Mensemble),stat=i_stat)
   allocate(indxb_pol(pol_buff),stat=i_stat)



## Currents
#### source/Measurement/prn_currents.f90

   ! Input parameters
   integer          :: quant_axis            !< Quantization axis for calculating psi
   integer          :: current_step          !< Interval for sampling currents
   integer          :: current_buff          !< Buffer size for currents
   character(len=1) :: do_currents           !< Measure magnon and heat currents

   ! Definition of current arrays
   real(dblprec), dimension(:), allocatable        :: indxb_currents   !< Step counter for the currents
   real(dblprec), dimension(:,:,:), allocatable    :: psi_amp_b        !< Buffer for the heat current
   real(dblprec), dimension(:,:,:), allocatable    :: psi_phase_b      !< Buffer for the magnon current
   real(dblprec), dimension(:,:,:), allocatable    :: heat_current_b   !< Buffer for the heat current
   real(dblprec), dimension(:,:,:), allocatable    :: heat_current2_b  !< Buffer for the heat current
   real(dblprec), dimension(:,:,:), allocatable    :: magnon_current_b !< Buffer for the magnon current
   complex(dblprec), dimension(:,:,:), allocatable :: psi_b
   complex(dblprec), dimension(:,:), allocatable   :: psi_curr           !< Psi variable for the current time step
   complex(dblprec), dimension(:,:), allocatable   :: psi_prev           !< Psi variable for the previous time step
   complex(dblprec), dimension(:,:), allocatable   :: psi_dot            !< Derivative of the psi variable
   complex(dblprec), dimension(:,:), allocatable   :: psi_dot2           !< Derivative of the psi variable
   complex(dblprec), dimension(:,:,:), allocatable :: psi_time
   real(dblprec), dimension(:,:), allocatable      :: heat_current       !< Heat current per site
   real(dblprec), dimension(:,:), allocatable      :: heat_current2      !< Heat current per site
   real(dblprec), dimension(:,:), allocatable      :: magnon_current     !< Magnon current per time
   real(dblprec), dimension(:,:), allocatable      :: ave_heat_current   !< Average heat current per site
   real(dblprec), dimension(:,:), allocatable      :: ave_magnon_current !< Average magnon current per time
   real(dblprec), dimension(:,:), allocatable      :: psi_amp            !< Amplitude of the complex number
   real(dblprec), dimension(:,:), allocatable      :: psi_phase          !< Phase of the complex number

   ! Allocation of current arrays
   allocate(heat_current(3,Natom),stat=i_stat)
   allocate(heat_current2(3,Natom),stat=i_stat)
   allocate(magnon_current(3,Natom),stat=i_stat)
   allocate(ave_heat_current(3,Natom),stat=i_stat)
   allocate(ave_magnon_current(3,Natom),stat=i_stat)
   allocate(psi_curr(Natom,Mensemble),stat=i_stat)
   allocate(psi_prev(Natom,Mensemble),stat=i_stat)
   allocate(psi_dot(Natom,Mensemble),stat=i_stat)
   allocate(psi_dot2(Natom,Mensemble),stat=i_stat)
   allocate(psi_time(3,Natom,Mensemble),stat=i_stat)
   allocate(psi_amp(Natom,Mensemble),stat=i_stat)
   allocate(psi_phase(Natom,Mensemble),stat=i_stat)
   allocate(psi_b(Natom,Mensemble,current_buff),stat=i_stat)
   allocate(heat_current_b(3,Natom,current_buff),stat=i_stat)
   allocate(psi_amp_b(Natom,Mensemble,current_buff),stat=i_stat)
   allocate(psi_phase_b(Natom,Mensemble,current_buff),stat=i_stat)
   allocate(heat_current2_b(3,Natom,current_buff),stat=i_stat)
   allocate(magnon_current_b(3,Natom,current_buff),stat=i_stat)
   allocate(indxb_currents(current_buff),stat=i_stat)



## Topology
#### source/Measurement/topology.f90

   ! Parameters for the printing
   integer :: skyno_step !< Interval for sampling the skyrmion number
   integer :: skyno_buff !< Buffer size for the sampling of the skyrmion number
   character(len=1) :: skyno      !< Perform skyrmion number measurement
   character(len=1) :: do_proj_skyno !< Perform type dependent skyrmion number measurement
   character(len=1) :: do_skyno_den  !< Perform site dependent skyrmion number measurement
   character(len=1) :: do_skyno_cmass  !< Perform center-of-mass skyrmion number measurement
   integer :: nsimp !< Number of simplices

   ! Definition of topology arrays
   integer, dimension(:,:), allocatable :: simp !< Array for storing Delaunay simplices

   ! Allocation of topology arrays
   allocate(simp(3,nsimp))



## Energy
#### source/Hamiltonian/energy.f90

   ! Input variables
   integer, intent(in) :: NA           !< Number of atoms in one cell
   integer, intent(in) :: N1           !< Number of cell repetitions in x direction
   integer, intent(in) :: N2           !< Number of cell repetitions in y direction
   integer, intent(in) :: N3           !< Number of cell repetitions in z direction
   integer, intent(in) :: nHam         !< Number of atoms in Hamiltonian
   integer, intent(in) :: mstep        !< Current simulation step
   integer, intent(in) :: Natom        !< Number of atoms in system
   integer, intent(in) :: Nchmax       !< Number of chemical type
   integer, intent(in) :: conf_num     !< number of configurations for LSF
   integer, intent(in) :: Mensemble    !< Number of ensembles
   integer, intent(in) :: stop_atom    !< Atom to end loop for
   integer, intent(in) :: Num_macro    !< Number of macrocells in the system
   integer, intent(in) :: start_atom   !< Atom to start loop for
   integer, intent(in) :: plotenergy   !< Calculate and plot energy (0/1)
   real(dblprec), intent(in) :: Temp         !< Temperature
   real(dblprec), intent(in) :: delta_t      !< Current time step
   character(len=1), intent(in) :: do_lsf    !< Including LSF energy
   character(len=1), intent(in) :: lsf_field          !< LSF field contribution (Local/Total)
   character(len=1), intent(in) :: lsf_interpolate    !< Interpolate LSF or not
   character(len=1), intent(in) :: real_time_measure  !< Display measurements in real time
   character(len=8), intent(in) :: simid              !< Name of simulation
   integer, dimension(Natom), intent(in) :: cell_index            !< Macrocell index for each atom
   integer, dimension(Num_macro), intent(in) :: macro_nlistsize   !< Number of atoms per macrocell
   real(dblprec), dimension(Natom,Mensemble), intent(in) :: mmom     !< Current magnetic moment
   real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: emom   !< Current unit moment vector
   real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: emomM  !< Current magnetic moment vector
   real(dblprec), dimension(3,Num_macro,Mensemble), intent(in) :: emomM_macro !< The full vector of the macrocell magnetic moment
   real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: external_field  !< External magnetic field
   real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: time_external_field !< External time-dependent magnetic field

   ! .. Output Variables
   real(dblprec), intent(out) :: totene !< Total energy

   ! Definition of energy arrays
   real(dblprec), dimension(:), allocatable :: ene_xc
   real(dblprec), dimension(:), allocatable :: ene_dm
   real(dblprec), dimension(:), allocatable :: ene_sa
   real(dblprec), dimension(:), allocatable :: ene_bq
   real(dblprec), dimension(:), allocatable :: ene_ring
   real(dblprec), dimension(:), allocatable :: ene_pd
   real(dblprec), dimension(:), allocatable :: ene_chir
   real(dblprec), dimension(:), allocatable :: ene_dip
   real(dblprec), dimension(:), allocatable :: ene_ani
   real(dblprec), dimension(:), allocatable :: ene_lsf
   real(dblprec), dimension(:), allocatable :: ene_ext
   real(dblprec), dimension(:), allocatable :: ene_pair
   real(dblprec), dimension(:), allocatable :: ene_bqdm
   real(dblprec), dimension(:), allocatable :: energy
   real(dblprec), dimension(:,:,:), allocatable :: bfield_dip
   real(dblprec), dimension(:,:,:), allocatable :: site_energy

   ! Allocation of energy arrays
   allocate(site_energy(11,Natom,Mensemble),stat=i_stat)
   allocate(bfield_dip(3,Natom,Mensemble),stat=i_stat)
   allocate(ene%energy(Mensemble),stat=i_stat)
   allocate(ene%ene_xc(Mensemble),stat=i_stat)
   allocate(ene%ene_dm(Mensemble),stat=i_stat)
   allocate(ene%ene_sa(Mensemble),stat=i_stat)
   allocate(ene%ene_bq(Mensemble),stat=i_stat)
   allocate(ene%ene_ring(Mensemble),stat=i_stat)
   allocate(ene%ene_pd(Mensemble),stat=i_stat)
   allocate(ene%ene_chir(Mensemble),stat=i_stat)
   allocate(ene%ene_ani(Mensemble),stat=i_stat)
   allocate(ene%ene_bqdm(Mensemble),stat=i_stat)
   allocate(ene%ene_ext(Mensemble),stat=i_stat)
   allocate(ene%ene_dip(Mensemble),stat=i_stat)
   allocate(ene%ene_lsf(Mensemble),stat=i_stat)
   allocate(ene%ene_pair(Mensemble),stat=i_stat)



## Magnetic field
#### source/Fields/prn_fields.f90

   ! Input parameters to be read
   integer :: beff_step               !< Interval between consecutive prints of the total effective field
   integer :: beff_buff               !< Buffer size for the total field
   integer :: binteff_step            !< Interval between consecutive prints of the internal effective field
   integer :: binteff_buff            !< Buffer size for the internal field
   integer :: thermfield_step         !< Interval between thermal field trajectories
   integer :: thermfield_buff         !< Buffer size for the stochastic field
   integer :: torques_step            !< Interval between consecutive prints of the resulting torques
   integer :: torques_buff            !< Buffer size for the resulting torques
   integer :: larm_step               !< Interval between consecutive prints of the larmor frequencies
   integer :: larm_buff               !< Buffer size for the larmor frequencies
   integer :: larm_dos_size           !< Number of windows for larmor dos histogram
   character(len=1) :: do_prn_beff    !< Flag governing file output of total effective fields (Y/N)
   character(len=1) :: do_prn_binteff !< Flag governing file output of internal effective fields (Y/N)
   character(len=1) :: do_prn_torques !< Flag governing file output of resulting torques (Y/N)
   character(len=1) :: do_thermfield  !< Thermal fields trajectory
   character(len=1) :: do_larmor_loc  !< Calculate local precession frequencies from local field (Y/N)
   character(len=1) :: do_larmor_dos  !< Calculate average precession frequencies from local field (Y/N)

   ! Local variables for buffering and indexing of fields
   integer :: bcount_beff    !< Counter of buffer for total field
   integer :: bcount_binteff !< Counter of buffer for internal field
   integer :: bcount_torques !< Counter of buffer for toruqes
   integer :: bcount_therm   !< Counter of buffer for stochastic field
   integer :: bcount_larm    !< Counter of buffer for Larmor frequencies

   ! Definition of magnetic field arrays
   real(dblprec), dimension(:), allocatable :: indxb_beff         !< Step counter for total field
   real(dblprec), dimension(:), allocatable :: indxb_binteff      !< Step counter for internal field
   real(dblprec), dimension(:), allocatable :: indxb_torques      !< Step counter for resulting torques
   real(dblprec), dimension(:), allocatable :: indxb_larm         !< Step counter for total field
   real(dblprec), dimension(:), allocatable :: indxb_therm        !< Step counter for stochastic field
   real(dblprec), dimension(:,:,:,:), allocatable :: beffb        !< Buffer the site dependent total field
   real(dblprec), dimension(:,:,:,:), allocatable :: binteffb     !< Buffer the site resulting torques
   real(dblprec), dimension(:,:,:,:), allocatable :: torquesb     !< Buffer the site dependent internal field
   real(dblprec), dimension(:,:,:), allocatable :: larmb          !< Buffer the site dependent larmor frequencies
   real(dblprec), dimension(:,:,:,:), allocatable :: therm_fieldb !< Buffer the site dependent stochastic field
   real(dblprec), dimension(:), allocatable :: larm_dos       !< Histogram array for Larmor frequencies

   ! Allocation of magnetic field arrays
   allocate(therm_fieldb(3,Natom,thermfield_buff,Mensemble),stat=i_stat)
   allocate(indxb_therm(thermfield_buff),stat=i_stat)
   allocate(beffb(3,Natom,beff_buff,Mensemble),stat=i_stat)
   allocate(indxb_beff(beff_buff),stat=i_stat)
   allocate(binteffb(6,Natom,binteff_buff,Mensemble),stat=i_stat)
   allocate(indxb_binteff(binteff_buff),stat=i_stat)
   allocate(torquesb(6,Natom,torques_buff,Mensemble),stat=i_stat)
   allocate(indxb_torques(torques_buff),stat=i_stat)
   allocate(larmb(Natom,larm_buff,Mensemble),stat=i_stat)
   allocate(indxb_larm(larm_buff),stat=i_stat)
   allocate(larm_dos(0:larm_dos_size))



## Ionic displacement and velocity
#### source/Lattice/prn_latticetrajectories.f90

   ! Input parameters
   integer :: lntraj         !< Number of displacement trajectories to sample
   integer :: ltottraj_step  !< Interval for sampling displacement trajectories
   integer :: ltottraj_buff  !< Buffer size for displacement trajectories
   character(len=1) :: do_ltottraj !< Measure displacements

   ! Definition of ionic displacement and velocity arrays
   integer, dimension(:), allocatable :: ltraj_step !< Interval for sampling individual displacement trajectories
   integer, dimension(:), allocatable :: ltraj_buff !< Buffer size for individual displacement trajectories
   integer, dimension(:), allocatable :: ltraj_atom !< List of atoms to sample displacement trajectories for
   integer, dimension(:), allocatable :: bcount_traj   !< Counter of buffer of displacements
   real(dblprec), dimension(:), allocatable :: scount_traj      !< Counter of sampling of displacements
   real(dblprec), dimension(:), allocatable :: indxb            !< Step counter for displacements
   real(dblprec), dimension(:,:), allocatable :: indxb_traj     !< Step counter for individual displacements
   real(dblprec), dimension(:,:,:,:), allocatable :: uvecb      !< Buffer for all individual displacements
   real(dblprec), dimension(:,:,:,:), allocatable :: uvecb_traj !< Buffer for selected individual displacements
   real(dblprec), dimension(:,:,:,:), allocatable :: vvecb      !< Buffer for all individual velocities
   real(dblprec), dimension(:,:,:,:), allocatable :: vvecb_traj !< Buffer for selected individual velocities
   real(dblprec), dimension(:,:,:,:), allocatable :: llatb      !< Buffer for all individual angular momenta
   real(dblprec), dimension(:,:,:,:), allocatable :: llatb_traj !< Buffer for selected individual angular momenta

   ! Allocation of displacement arrays
   allocate(scount_traj(lntraj),stat=i_stat)
   allocate(bcount_traj(lntraj),stat=i_stat)
   allocate(uvecb_traj(3,maxval(ltraj_buff),lntraj,Mensemble),stat=i_stat)
   allocate(vvecb_traj(3,maxval(ltraj_buff),lntraj,Mensemble),stat=i_stat)
   allocate(llatb_traj(3,maxval(ltraj_buff),lntraj,Mensemble),stat=i_stat)
   allocate(indxb_traj(maxval(ltraj_buff),lntraj),stat=i_stat)
   allocate(uvecb(3,Natom,ltottraj_buff,Mensemble),stat=i_stat)
   allocate(vvecb(3,Natom,ltottraj_buff,Mensemble),stat=i_stat)
   allocate(llatb(3,Natom,ltottraj_buff,Mensemble),stat=i_stat)
   allocate(indxb(ltottraj_buff),stat=i_stat)



## Ionic averages
#### source/Lattice/prn_latticeaverages.f90

   ! Printing definitions
   integer :: lavrg_step !< Interval for sampling average displacements
   integer :: lavrg_buff !< Buffer size for average displacements
   character(len=1) :: do_lavrg           !< Measure average displacments (Y/N)
   character(len=1) :: do_proj_lavrg      !< Measure projected displacements (Y/A/N)
   character(len=1) :: do_projch_lavrg    !< Measure chemically projected displacements (Y/N)
   character(len=1) :: do_lenergy         !< Measure susceptibility, and specific heat(Y/N)

   ! Definition of ionic averages
   real(dblprec), dimension(:,:), allocatable :: dcoord !< Coordinates of atoms
   real(dblprec), dimension(:), allocatable       :: indxb_uavrg       !< Step counter for average displacement
   real(dblprec), dimension(:,:,:), allocatable   :: uavrg_buff        !< Buffer for average displacements
   real(dblprec), dimension(:,:,:), allocatable   :: uavrg2_buff_proj  !< Buffer for squared projected displacements
   real(dblprec), dimension(:,:,:,:), allocatable :: uavrg_buff_proj   !< Buffer for projected displacements
   real(dblprec), dimension(:,:,:,:), allocatable :: uavrg_buff_projch !< Buffer for chemical projected displacements
   real(dblprec), dimension(:,:,:), allocatable   :: vavrg_buff        !< Buffer for average velocities
   real(dblprec), dimension(:,:,:), allocatable   :: vavrg2_buff_proj  !< Buffer for squared projected velocities
   real(dblprec), dimension(:,:,:,:), allocatable :: vavrg_buff_proj   !< Buffer for projected velocities
   real(dblprec), dimension(:,:,:,:), allocatable :: vavrg_buff_projch !< Buffer for chemical projected velocities
   real(dblprec), dimension(:,:,:), allocatable   :: pavrg_buff        !< Buffer for average ionic momentum
   real(dblprec), dimension(:,:,:), allocatable   :: langavrg_buff     !< Buffer for average ionic angular momentum
   real(dblprec), dimension(:,:), allocatable     :: ldpotenrg_buff    !< Buffer for average ionic potential energy
   real(dblprec), dimension(:,:), allocatable     :: sdpotenrg_buff    !< Buffer for average magnetic energy
   real(dblprec), dimension(:,:), allocatable     :: sldpotenrg_buff   !< Buffer for average spin-lattice potential energy
   real(dblprec), dimension(:,:), allocatable     :: totpotenrg_buff   !< Buffer for average total potential energy
   real(dblprec), dimension(:,:), allocatable     :: kinenrg_buff      !< Buffer for average ionic kinetic energy
   real(dblprec), dimension(:,:), allocatable     :: totenrg_buff      !< Buffer for average ionic total energy
   real(dblprec), dimension(:,:), allocatable     :: heatcap_buff      !< Buffer for ionic heat capacity
   real(dblprec), dimension(:,:), allocatable     :: iontemp_buff      !< Buffer for ionic temperature

   ! Allocation of ionic average arrays
   allocate(uavrg_buff(4,lavrg_buff,Mensemble),stat=i_stat)
   allocate(uavrg_buff_proj(3,NA,lavrg_buff,Mensemble),stat=i_stat)
   allocate(uavrg2_buff_proj(NA,lavrg_buff,Mensemble),stat=i_stat)
   allocate(vavrg_buff(4,lavrg_buff,Mensemble),stat=i_stat)
   allocate(vavrg_buff_proj(3,NA,lavrg_buff,Mensemble),stat=i_stat)
   allocate(vavrg2_buff_proj(NA,lavrg_buff,Mensemble),stat=i_stat)
   allocate(pavrg_buff(4,lavrg_buff,Mensemble),stat=i_stat)
   allocate(langavrg_buff(4,lavrg_buff,Mensemble),stat=i_stat)
   allocate(dcoord(3,Natom),stat=i_stat)
   allocate(ldpotenrg_buff(lavrg_buff,Mensemble),stat=i_stat)
   allocate(sdpotenrg_buff(lavrg_buff,Mensemble),stat=i_stat)
   allocate(sldpotenrg_buff(lavrg_buff,Mensemble),stat=i_stat)
   allocate(totpotenrg_buff(lavrg_buff,Mensemble),stat=i_stat)
   allocate(kinenrg_buff(lavrg_buff,Mensemble),stat=i_stat)
   allocate(totenrg_buff(lavrg_buff,Mensemble),stat=i_stat)
   allocate(heatcap_buff(lavrg_buff,Mensemble),stat=i_stat)
   allocate(iontemp_buff(lavrg_buff,Mensemble),stat=i_stat)
   allocate(indxb_uavrg(lavrg_buff),stat=i_stat)
   allocate(uavrg_buff_projch(4,Nchmax,lavrg_buff,Mensemble),stat=i_stat)
   allocate(vavrg_buff_projch(4,Nchmax,lavrg_buff,Mensemble),stat=i_stat)



## Ionic force field
#### source/Lattice/prn_latticefields.f90

   ! Input parameters to be read
   integer :: eeff_step               !< Interval between consecutive prints of the total effective field
   integer :: eeff_buff               !< Buffer size for the total field
   integer :: einteff_step               !< Interval between consecutive prints of the internal effective field
   integer :: einteff_buff               !< Buffer size for the internal field
   integer :: ethermfield_step         !< Interval between thermal field trajectories
   integer :: ethermfield_buff         !< Buffer size for the stochastic field
   character(len=1) :: do_prn_eeff    !< Flag governing file output of total effective fields (Y/N)
   character(len=1) :: do_prn_einteff    !< Flag governing file output of internal effective fields (Y/N)
   character(len=1) :: do_ethermfield  !< Thermal fields trajectory

   ! Definition of ionic force field arrays
   real(dblprec), dimension(:), allocatable :: indxb_eeff          !< Step counter for total field
   real(dblprec), dimension(:), allocatable :: indxb_einteff       !< Step counter for internal field
   real(dblprec), dimension(:), allocatable :: indxb_etherm        !< Step counter for stochastic field
   real(dblprec), dimension(:,:,:,:), allocatable :: eeffb         !< Buffer the site dependent total field
   real(dblprec), dimension(:,:,:,:), allocatable :: einteffb      !< Buffer the site dependent internal field
   real(dblprec), dimension(:,:,:,:), allocatable :: etherm_fieldb !< Buffer the site dependent stochastic field
   integer, dimension(:,:), allocatable :: simp !< Array for storing Delaunay simplices

   ! Allocation of ionic force field arrays
   real(dblprec), dimension(:), allocatable :: indxb_eeff          !< Step counter for total field
   real(dblprec), dimension(:), allocatable :: indxb_einteff       !< Step counter for internal field
   real(dblprec), dimension(:), allocatable :: indxb_etherm        !< Step counter for stochastic field
   real(dblprec), dimension(:,:,:,:), allocatable :: eeffb         !< Buffer the site dependent total field
   real(dblprec), dimension(:,:,:,:), allocatable :: einteffb      !< Buffer the site dependent internal field
   real(dblprec), dimension(:,:,:,:), allocatable :: etherm_fieldb !< Buffer the site dependent stochastic field
