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

   ! Definition of topology arrays
   integer, dimension(:,:), allocatable :: simp !< Array for storing Delaunay simplices

   ! Allocation of topology arrays
   allocate(simp(3,nsimp))



## Energy
#### source/Hamiltonian/energy.f90

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



## Ionic displacement



## Ionic velocity



## Ionic forces
