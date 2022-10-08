!------------------------------------------------------------------------------------
! MODULE: InputData
!> @brief New input data type for input from file
!> @copyright
!> GNU Public License.
!------------------------------------------------------------------------------------
module InputData

   use Parameters
   use Profiling

   implicit none
   !---------------------------------------------------------------------------------
   ! Geometry and composition
   !---------------------------------------------------------------------------------
   integer :: N1                                            !< Number of cell repetitions in x direction
   integer :: N2                                            !< Number of cell repetitions in y direction
   integer :: N3                                            !< Number of cell repetitions in z direction
   integer :: NA                                            !< Number of atoms in one cell
   integer :: NT                                            !< Number of types of atoms in the unit cell
   integer :: Sym                                           !< Assumed symmetry of the system
   integer :: nHam                                          !< Number of atoms for to be used in the Hamiltonian (Natom or NA)
   integer :: Natom                                         !< Number of atoms in system
   integer :: Natom_full                                    !< Number of atoms for full system (=Natom if not dilute)
   integer :: set_landeg                                    !< Set 'custom' value for gyromagnetic ration (1=yes,0=no)
   integer :: block_size                                    !< Size of the blocking parameter for the macro cell creation
   integer, dimension(:), allocatable :: acomp              !< Atom component
   integer, dimension(:), allocatable :: asite              !< Atom site
   integer, dimension(:), allocatable :: anumb_inp          !< Atom number in cell
   integer, dimension(:), allocatable :: atype_inp          !< Type of atom
   character :: posfiletype                                 !< posfile type (D)irect or (C)arteisian coordinates in posfile
   character :: do_sortcoup                                 !< Sort the entries of ncoup arrays (Y/N)
   character(len=1) :: BC1                                  !< Boundary conditions in x-direction
   character(len=1) :: BC2                                  !< Boundary conditions in y-direction
   character(len=1) :: BC3                                  !< Boundary conditions in z-direction
   character(len=35) :: posfile                             !< File name for coordinates
   real(dblprec) :: alat                                    !< Lattice parameter
   real(dblprec) :: scalefac                                !< Lattice vector rescaling factor (for debugging)
   real(dblprec) :: Landeg_glob                             !< Gyromagnetic ratio
   real(dblprec), dimension(3) :: C1                        !< First lattice vector
   real(dblprec), dimension(3) :: C2                        !< Second lattice vector
   real(dblprec), dimension(3) :: C3                        !< Third lattice vector
   real(dblprec), dimension(:,:), allocatable :: Bas        !< Coordinates for basis atoms
   real(dblprec), dimension(:,:), allocatable :: Bas0       !< Coordinates for basis atoms (raw)
   real(dblprec), dimension(:,:,:), allocatable :: Landeg_ch  !< Gyromagnetic ratio
   !---------------------------------------------------------------------------------
   ! Moment data
   !---------------------------------------------------------------------------------
   character(len=1) :: do_mom_legacy   !< Flag to print/read moments in legacy output
   character(len=1) :: renorm_coll                            !< Flag to force collienar calculation of the susceptibility for induced moments
   character(len=1) :: ind_mom_flag                           !< Flag for whether induced magnetic moments are present or not
   character(len=35) :: momfile                               !< File name for magnetic moments
   character(len=35) :: momfile_i                             !< File name for initial magn. configuration
   character(len=35) :: momfile_f                             !< File name for final magn. configuration
   real(dblprec) :: ind_tol                                   !< Tolerance for the induced moments
   real(dblprec) :: amp_rnd                                   !< Amplitude of random perturbation of the components of magnetic moments
   real(dblprec) :: amp_rnd_path                              !< Amplitue of the random perturbation for the path
   integer, dimension(:,:), allocatable :: ind_mom              !< Flag to decide whether a moment is induced or not
   real(dblprec), dimension(:,:,:), allocatable :: ammom_inp    !< Magnetic moment magnitudes from input (for alloys)
   real(dblprec), dimension(:,:,:,:), allocatable :: aemom_inp  !< Magnetic moment directions from input (for alloys)
   !---------------------------------------------------------------------------------
   ! Exchange data
   !---------------------------------------------------------------------------------
   integer :: maptype                        !< Format for input data (1=direct format,2=bgfm style)
   integer :: do_jtensor                     !< Use SKKR style exchange tensor (0=off, 1=on) --OLD: 2=with biquadratic exchange
   logical :: calc_jtensor                   !< Calculate or read tensor from input (T/F)
   logical :: map_multiple                   !< Allow for multiple couplings between atoms i and j (use only for small cells)
   integer :: max_no_shells                  !< Actual maximum number of shells
   integer, dimension(:), allocatable :: nn  !< Number of neighbour shells
   real(dblprec), dimension(:,:,:), allocatable :: redcoord       !< Coordinates for Heisenberg exchange couplings
   real(dblprec), dimension(:,:,:,:), allocatable :: jc_bq        !< Biquadratic exchange coupling
   real(dblprec), dimension(:,:,:,:,:), allocatable :: jc         !< Exchange couplings
   real(dblprec), dimension(:,:,:,:,:), allocatable :: jcD        !< Exchange couplings (DLM)
   real(dblprec), dimension(:,:,:,:,:,:), allocatable :: jc_tens  !< Tensorial exchange couplings (SKKR)
   character(len=1) :: exc_inter                         !< Interpolate Jij between FM and DLM states (Y(N)
   character(len=30) :: pre_jfile                        !< File name for exchange couplings
   character(len=30), dimension(:), allocatable :: jfile    !< File name for exchange couplings
   character(len=30), dimension(:), allocatable :: jfileD   !< File name for exchange couplings (DLM)
   real(dblprec) :: jij_scale                               !< Rescale Jij couplings manually
   !---------------------------------------------------------------------------------
   ! Anisotropy data
   !---------------------------------------------------------------------------------
   integer :: do_anisotropy                                          !< Read anisotropy data (1/0)
   real(dblprec) :: random_anisotropy_density                        !< Density for randomness in anisotropy
   logical :: random_anisotropy                                      !< Put random anisotropy in the sample (T/F)
   character(len=1) :: mult_axis                                     !< Flag to treat more than one anisotropy axis at the same time
   character(len=35) :: kfile                                        !< File name for anisotropy data
   integer, dimension(:,:), allocatable :: anisotropytype            !< Type of anisotropies (0-2)
   integer, dimension(:,:), allocatable :: anisotropytype_diff       !< Type of anisotropies when one is considering more than one anisotropy axis
   real(dblprec), dimension(:,:,:), allocatable :: anisotropy        !< Input data for anisotropies
   real(dblprec), dimension(:,:,:), allocatable :: anisotropy_diff   !< Input data for the second anisotropy axis
   !---------------------------------------------------------------------------------
   ! Parameters for energy minimization calculations
   !---------------------------------------------------------------------------------
   integer :: minalgo         !< Minimization algorithm (1=VPO,...)
   integer :: minitrmax       !< Maximum number of iterations in energy minimization procedure
   integer :: mintraj_step    !< Save configuration every 'mintraj_step' step during energy minimization
   real(dblprec) :: vpodt     !< VPO timestep
   real(dblprec) :: minftol   !< Convergence criterion in mRy
   real(dblprec) :: vpomass   !< VPO mass
   !---------------------------------------------------------------------------------
   ! Parameters for GNEB calculations
   !---------------------------------------------------------------------------------
   integer :: initpath              !< Path initialization method (1=geodesic, 2=read from file,...)
   integer :: mepitrmax             !< Maximum number of iterations in MEP finding procedure
   integer :: meptraj_step          !< Save configuration every 'meptraj_step' step during MEP finding procedure
   real(dblprec) :: spring          !< Magnitude of the spring constant in GNEB calculations
   real(dblprec) :: mepftol         !< Convergence criterion in the GNEB method, in mRy
   real(dblprec) :: mepftol_ci      !< Convergence criterion in the CI-GNEB method, in mRy
   character(len=1) :: do_gneb      !< Do GNEB calculations (Y/N)
   character(len=1) :: do_gneb_ci   !< Do CI-GNEB calculations (Y/N)
   character(len=1) :: do_norm_rx   !< Normalize reaction coordinate (Y/N)
   character(len=1) :: en_zero      !< Level of zero energy. 'I' - initial state; 'F' - final state; 'N' - 0.0
   character(len=1) :: relaxed_if   !< Use relaxed ini and fin states (Y/N)
   character(len=1) :: fixed_if     !< Fix endpoints of the path (Y) or allow them to move along energy isocontours (N)
   character(len=1) :: prn_gneb_fields !< Print the magnetic fields from GNEB
   !---------------------------------------------------------------------------------
   ! Parameters for Hessian calculations
   !---------------------------------------------------------------------------------
   character(len=1) :: do_hess_ini   !< Calculate Hessian at the initial state (Y/N)
   character(len=1) :: do_hess_fin   !< Calculate Hessian at the final state (Y/N)
   character(len=1) :: do_hess_sp    !< Calculate Hessian at the saddle point (Y/N)
   real(dblprec) :: eig_0            !< Minimal nonzero eigenvalue
   character(len=1) :: is_afm        !< is the texture AFM
   !---------------------------------------------------------------------------------
   ! Parameters for energy interpolation along the MEP
   !---------------------------------------------------------------------------------
   integer :: sample_num            !< Number of samples in the interpolated curve
   !---------------------------------------------------------------------------------
   ! Dzyaloshinskii-Moriya data
   !---------------------------------------------------------------------------------
   integer :: do_dm                                      !< Add Dzyaloshinskii-Moriya (DM) term to Hamiltonian (0/1)
   integer :: max_no_dmshells                            !< Actual maximum number of shells for DM interactions
   character(len=35) :: dmfile                           !< File name for Dzyaloshinskii-Moriya data
   real(dblprec) :: dm_scale                             !< Rescale Dij couplings manually
   integer, dimension(:), allocatable :: dm_nn           !< No. shells of neighbours for DM
   real(dblprec), dimension(:,:,:), allocatable :: dm_redcoord    !< Neighbour vectors for DM
   real(dblprec), dimension(:,:,:,:,:), allocatable :: dm_inpvect !< Neighbour vectors for DM
   !---------------------------------------------------------------------------------
   ! Scalar chirality data
   !---------------------------------------------------------------------------------
   integer :: do_chir                                       !< Add scalar chirality (CHIR) term to Hamiltonian (0/1)
   integer :: max_no_chirshells                             !< Actual maximum number of shells for CHIR interactions
   character(len=35) :: chirfile                            !< File name for Pseudo-Dipolar data
   integer, dimension(:), allocatable :: chir_nn            !< No. shells of neighbours for CHIR
   real(dblprec), dimension(:,:,:,:), allocatable :: chir_inpval     !< Interactions for CHIR
   real(dblprec), dimension(:,:,:,:), allocatable :: chir_redcoord   !< Neighbour vectors for CHIR
   !---------------------------------------------------------------------------------
   ! Pseudo-Dipolar data
   !---------------------------------------------------------------------------------
   integer :: do_pd                                      !< Add Pseudo-Dipolar (PD) term to Hamiltonian (0/1)
   integer :: max_no_pdshells                            !< Actual maximum number of shells for PD interactions
   character(len=35) :: pdfile                           !< File name for Pseudo-Dipolar data
   integer, dimension(:), allocatable :: pd_nn           !< No. shells of neighbours for PD
   real(dblprec), dimension(:,:,:), allocatable :: pd_redcoord       !< Neighbour vectors for PD
   real(dblprec), dimension(:,:,:,:,:), allocatable :: pd_inpvect    !< Neighbour vectors for PD
   !---------------------------------------------------------------------------------
   ! Biquadratic DM data (BIQDM)
   !---------------------------------------------------------------------------------
   integer :: do_biqdm                                      !< Add BIQDM term to Hamiltonian (0/1)
   integer :: max_no_biqdmshells                            !< Actual maximum number of shells for BIQDM interactions
   character(len=35) :: biqdmfile                           !< File name for BIQDM data
   integer, dimension(:), allocatable :: biqdm_nn           !< No. shells of neighbours for BIQDM
   real(dblprec), dimension(:,:,:), allocatable :: biqdm_redcoord       !< Neighbour vectors for BIQDM
   real(dblprec), dimension(:,:,:,:,:), allocatable :: biqdm_inpvect    !< Neighbour vectors for BIQDM
   !---------------------------------------------------------------------------------
   ! Biquadratic exchange data
   !---------------------------------------------------------------------------------
   integer :: do_bq                                   !< Add biquadratic (BQ) term to Hamiltonian (0/1)
   integer :: max_no_bqshells                         !< Actual maximum number of shells for BQ interactions
   character(len=35) :: bqfile                        !< File name for biquadratic data
   integer, dimension(:), allocatable :: bq_nn        !< No. shells of neighbours for BQ
   real(dblprec), dimension(:,:,:), allocatable :: bq_redcoord    !< Neighbour vectors for BQ
   !---------------------------------------------------------------------------------
   ! Dipole-dipole data
   !---------------------------------------------------------------------------------
   integer :: do_dip                      !< Calculate dipole-dipole contribution (0=Off, 1=Brute Force, 2=macrocell)
   character(len=1) :: read_dipole        !< Flag to read the dipole-dipole tensor from file
   character(len=1) :: print_dip_tensor   !< Flag to print the dipole tensor
   character(len=30) :: qdip_files        !< Input file that contains the dipole-dipole tensor
   !---------------------------------------------------------------------------------
   ! Ewald summation data
   !---------------------------------------------------------------------------------
   character(len=1) :: do_ewald  !< Perform Ewald summation
   integer, dimension(3) :: RMAX !< Maximum number of cells in real space taken into account
   integer, dimension(3) :: KMAX !< Maximum number of cells in Fourier soace taken into account
   real(dblprec) :: Ewald_alpha  !< Ewald parameter
   !---------------------------------------------------------------------------------
   ! Simulation paramters
   !---------------------------------------------------------------------------------
   character(len=8) :: simid  !< Name of simulation
   integer :: Mensemble       !< Number of independent simulations
   integer :: tseed           !< Temperature seed
   !$omp threadprivate(tseed)
   logical*1 :: para_rng      !< Use threaded random number generation (T/F)
   integer :: llg             !< Type of equation of motion (1=LLG)
   integer :: nstep           !< Number of steps in measurement phase
   integer :: SDEalgh         !< Solver for equations of motion (1-5)
   integer :: ipSDEalgh       !< Initial phase solver for equations of motion (1-5)
   character(len=1) :: aunits !< Atomic units to simulate model Hamiltonians (Y/N)
   character(len=1) :: perp   !< Remove the component of B_i parallel to m_i (Y/N)
   !---------------------------------------------------------------------------------
   ! Tasks
   !---------------------------------------------------------------------------------
   integer :: mompar           !< Parametrization of magnetic moment magnitudes (0=no)
   integer :: heisout          !< Print debug information from heisge
   integer :: evolveout        !< Print debug information from evolution
   integer :: plotenergy       !< Calculate and plot energy (0/1)
   integer :: do_hoc_debug     !< Print higher order couplings debug information (0/1)
   integer :: do_prnstruct     !< Print Hamiltonian information (0/1)
   integer :: do_prn_poscar    !< Print geometry on POSCAR format (0/1)
   integer :: do_prn_elk       !< Print geometry on ELK format (0/1)
   integer :: do_read_elk      !< Read geometry on ELK format (0/1)
   integer :: compensate_drift !< Correct for drift in RNG
   character(len=1) :: do_sparse    !< Use sparse linear algebra for effective field evaluation (T/F)
   character(len=1) :: do_reduced   !< Use reduced formulation of Hamiltonian (Y/N)
   !---------------------------------------------------------------------------------
   ! Measurement phase
   !---------------------------------------------------------------------------------
   real(dblprec) :: Temp                  !< Temperature
   real(dblprec) :: delta_t               !< Time step
   real(dblprec) :: relaxtime             !< Relaxation time for inertial regime (in LLG-I equation)
   real(dblprec) :: mplambda1             !< Damping parameter for measurement phase
   real(dblprec) :: mplambda2             !< Additional damping parameter measurement phase
   character(len=1) :: mode               !< Simulation mode (S=SD, M=MC, H=MC Heat Bath, P=LD, C=SLD, G=GNEB)
   real(dblprec) , dimension(3) :: hfield !< Applied magnetic field
   !---------------------------------------------------------------------------------
   ! Damping data
   !---------------------------------------------------------------------------------
   integer :: mpNlines                                           !< Number of lines of the damping file
   real(dblprec), dimension(:), allocatable :: mpdamping1        !< Damping parameter per atom in the unit cell for measurement phase
   real(dblprec), dimension(:), allocatable :: mpdamping2        !< Aditional damping parameter per atom in the unit cell for measurement phase
   real(dblprec), dimension(:,:), allocatable :: mpdampingalloy1 !< Damping parameter per atom in the unit cell for measurement phase (alloy)
   real(dblprec), dimension(:,:), allocatable :: mpdampingalloy2 !< Aditional damping parameter per atom in the unit cell for measurement phase (alloy)
   character(len=1) :: do_site_damping                           !< Flag for site dependent damping in measurement phase
   character(len=35) :: mp_dampfile                              !< File name for damping parameter of the atoms in the unit cell for measurement phase
   !---------------------------------------------------------------------------------
   ! Monte Carlo
   !---------------------------------------------------------------------------------
   integer :: do_mc        !< Do Monte Carlo instead of SD (Y/N)
   integer :: mcnstep      !< Number of Monte Carlo steps
   integer :: mcavrg_step  !< Sampling interval for averages
   integer :: mcavrg_buff  !< Buffer size for averages
   !---------------------------------------------------------------------------------
   ! Initial phase
   !---------------------------------------------------------------------------------
   integer :: ipnphase           !< Number of SD initial phases
   integer :: ipmcnphase         !< Number of MC initial phases
   character(len=1) :: ipmode    !< Simulation mode for initial phase (S=SD, M=MC, H=MC Heat Bath, G = EM)
   integer, dimension(:), allocatable :: ipnstep   !< Number of steps in initial phase
   integer, dimension(:), allocatable :: ipmcnstep !< Number of Monte Carlo steps for initial phase
   real(dblprec), dimension(3) :: iphfield         !< Applied magnetic field for initial phase
   real(dblprec), dimension(:), allocatable :: ipdelta_t !< Time step for initial phase
   real(dblprec), dimension(:), allocatable :: iplambda1 !< Damping parameter for initial phase
   real(dblprec), dimension(:), allocatable :: iplambda2 !< Additional damping parameter for initial phase
   real(dblprec), dimension(:), allocatable :: ipTemp    !< Temperature for initial phase
   !---------------------------------------------------------------------------------
   ! Initial phase Damping data
   !---------------------------------------------------------------------------------
   real(dblprec), dimension(:,:), allocatable :: ipdamping1        !< Damping parameter for initial phase per atom in the unit cell
   real(dblprec), dimension(:,:), allocatable :: ipdamping2        !< Additional damping parameter for initial phase per atom in the unit cell
   real(dblprec), dimension(:,:,:), allocatable :: ipdampingalloy1 !< Damping parameter for initial phase per atom in the unit cell (alloy)
   real(dblprec), dimension(:,:,:), allocatable :: ipdampingalloy2 !< Additional damping parameter for initial phase per atom in the unit cell (alloy)
   character(len=1) :: do_site_ip_damping !< Flag for site dependent dampin in the initial phase
   character(len=35) :: ip_dampfile       !< File name for damping parameter of the atoms in the unit cell for initial phase
   !---------------------------------------------------------------------------------
   ! Init mag
   !---------------------------------------------------------------------------------
   integer :: mseed            !< Seed for initialization of magnetic moments
   integer :: roteul           !< Global rotation of magnetic moments (0/1)
   integer :: initmag          !< Mode of initialization of magnetic moments (1-4)
   integer :: initneigh        !< Raman neighbour spin index
   real(dblprec) :: phi0       !< Cone angle phi
   real(dblprec) :: theta0     !< Cone angle theta
   real(dblprec) :: initimp    !< Size of impurity magnetic moment
   real(dblprec) :: initconc   !< Concentration of vacancies or two magnon Raman spin flips
   real(dblprec) :: initrotang !< Rotation angle phase for initial spin spiral
   real(dblprec), dimension(3) :: rotang        !< Euler angles for global rotation of magnetic moments
   real(dblprec), dimension(3) :: initrotvec    !< rotation vector for initial spin spiral
   real(dblprec), dimension(3) :: initpropvec   !< propagation vector for initial spin spiral
   character(len=1) :: initexc      !< Mode of excitation of initial magnetic moments (I=vacancies, R=two magnon Raman, F=no)
   character(len=35) :: restartfile !< File containing restart information
   !---------------------------------------------------------------------------------
   ! De-magnetization field
   !---------------------------------------------------------------------------------
   real(dblprec) :: demagvol  !< Volume for modelling demagnetization field
   character(len=1) :: demag  !< Model demagnetization field (Y/N)
   character(len=1) :: demag1 !< Model demagnetization field in x-direction (Y/N)
   character(len=1) :: demag2 !< Model demagnetization field in y-direction (Y/N)
   character(len=1) :: demag3 !< Model demagnetization field in z-direction (Y/N)
   !---------------------------------------------------------------------------------
   ! Localized magnetic fields
   !---------------------------------------------------------------------------------
   real(dblprec), dimension(:,:), allocatable :: sitenatomfld
   !---------------------------------------------------------------------------------
   ! Magnetic field pulse
   !---------------------------------------------------------------------------------
   integer :: do_bpulse              !< Add magnetic field pulse (1=no, 1-4 for different shapes and 5-7 with isite pulse and mwf)
   character :: locfield             !< File name for local applied field info
   character(len=35) :: bpulsefile   !< File name for magnetic field pulse info
   character(len=35) :: siteatomfile !< File name for local applied field info
   character(len=35) :: locfieldfile
   !---------------------------------------------------------------------------------
   ! LSF simulation
   !---------------------------------------------------------------------------------
   integer :: conf_num     !< Number of configurations for LSF
   integer :: gsconf_num   !< Ground state configuration in LSF
   integer :: lsf_metric   !< Metric in the phase space integration (1=Murata-Doniach, 2=Jacobian)
   real(dblprec) :: lsf_window   !< Range of moment variation in LSF
   character(len=1) :: do_lsf    !< (Y/N) Do LSF for MC. ASD LSF NOT IMPLEMENTED YET ASK FAN FOR DETAILS
   character(len=1) :: lsf_field
   character(len=1) :: lsf_interpolate
   character(len=35) :: lsffile !< File name for the configuration dependent LSF energy
   !---------------------------------------------------------------------------------
   ! Measurements
   !---------------------------------------------------------------------------------
   integer :: spintemp_step !< Interval for measuring spin temperature
   character(len=1) :: logsamp !< Sample measurements logarithmically
   character(len=1) :: do_spintemp !< Measure spin temperature
   character(len=1) :: real_time_measure !< Measurements displayed in real time
   !---------------------------------------------------------------------------------
   ! Random alloy data (chemical data)
   !---------------------------------------------------------------------------------
   integer :: Nchmax     !< Max number of chemical components on each site in cell
   integer :: do_ralloy  !< Random alloy simulation (0/1)
   integer, dimension(:), allocatable :: Nch            !< Number of chemical components on each site in cell
   integer, dimension(:,:), allocatable :: achtype_ch   !< Chemical type of atoms from input
   real(dblprec), dimension(:,:), allocatable :: chconc !< Chemical concentration on sites
   !---------------------------------------------------------------------------------
   ! File version info
   !---------------------------------------------------------------------------------
   logical :: oldformat
   !---------------------------------------------------------------------------------
   ! Random number related flag
   !---------------------------------------------------------------------------------
   character(len=1) :: rngpol   !< Use spherical coordinates for RNG generation (Y/N)
   character(len=1) :: ziggurat !< Use ziggurat transform for RNG instead of Box-Muller transform (Y/N)
   !---------------------------------------------------------------------------------
   ! GPU related flags
   !---------------------------------------------------------------------------------
   integer :: gpu_mode     !< What GPU mode to use (0 = FORTRAN, 1 = CUDA, 2 = C/C++)
   integer :: gpu_rng      !< What CURAND RNG to use (0 = DEFAULT, 1 = XORWOW, 2 = MRG32K3A, 3 = MTGP32)
   integer :: gpu_rng_seed !< Seed for RNG. If 0, the current time will be used.
   !---------------------------------------------------------------------------------
   ! I/O OVF related flags
   !---------------------------------------------------------------------------------
   character(len=1) :: prn_ovf  !< Print the magnetization data in the ovf format
   character(len=1) :: read_ovf !< Read the magnetization data in the ovf format

contains

   !---------------------------------------------------------------------------------
   ! SUBROUTINE: set_input_defaults
   !> @brief Sets default values to input variables
   !---------------------------------------------------------------------------------
   subroutine set_input_defaults()
      !
      implicit none

      real(dblprec) :: one=1.0_dblprec
      real(dblprec) :: zero=0.0_dblprec

      !Geometry and composition
      N1                = 0
      N2                = 0
      N3                = 0
      NA                = 0
      BC1               = '0'
      BC2               = '0'
      BC3               = '0'
      C1                = (/one,zero,zero/)
      C2                = (/zero,one,zero/)
      C3                = (/zero,zero,one/)
      Landeg_glob       = 2.0_dblprec
      set_landeg        = 0
      NT                = 0
      Sym               = 0
      do_sortcoup       = 'N'
      posfile           = ''
      posfiletype       = 'C'
      alat              = 1.0_dblprec
      scalefac          = 1.0_dblprec
      momfile           = ''
      momfile_i         = 'momfile_i'
      momfile_f         = 'momfile_f'
      amp_rnd           = 0.0_dblprec
      amp_rnd_path      = 0.0_dblprec
      block_size        = 1
      relaxed_if        = 'Y'
      fixed_if          = 'Y'

      !Induced moment data
      ind_mom_flag      = 'N'
      renorm_coll       = 'N'
      ind_tol           = 0.0010_dblprec

      !Exchange data
      maptype           = 1
      exc_inter         = 'N'
      map_multiple      = .false.
      jij_scale         = 1.0_dblprec

      !Anisotropy data
      kfile             = 'kfile'
      do_anisotropy     = 0
      mult_axis         = 'N'

      !Dzyaloshinskii-Moriya data
      dmfile            = 'dmfile'
      do_dm             = 0
      dm_scale          = 1.0_dblprec

      !Pseudo-Dipolar data
      pdfile            = 'pdfile'
      do_pd             = 0

      !Biquadratic DM data
      biqdmfile         = 'biqdmfile'
      do_biqdm          = 0

      !Biquadratic exchange data
      bqfile            = 'bqfile'
      do_bq             = 0

      !Tensorial exchange (SKKR) data
      do_jtensor        = 0
      calc_jtensor      = .true.

      !Dipole-dipole data
      do_dip            = 0
      print_dip_tensor  = 'N'
      read_dipole       = 'N'
      qdip_files        = 'qdip_file'

      !Ewald summation data
      do_ewald          = 'N'
      Ewald_alpha       = 0.0_dblprec
      KMAX              = (/0,0,0/)
      RMAX              = (/0,0,0/)

      !Parameters for energy minimization calculations
      minalgo           = 1
      minftol           = 0.000000001_dblprec
      mintraj_step      = 100
      vpodt             = 0.010_dblprec
      vpomass           = 1.0_dblprec
      minitrmax         = 10000000

      !Parameters for GNEB calculations
      initpath          = 1
      spring            = 0.50_dblprec
      mepftol           = 0.0010_dblprec
      mepftol_ci        = 0.000010_dblprec
      mepitrmax         = 10000000
      meptraj_step      = 100
      do_gneb           = 'Y'
      do_gneb_ci        = 'N'
      do_norm_rx        = 'N'
      en_zero           = 'N'
      prn_gneb_fields   = 'N'

      !Parameters for Hessian calculations
      do_hess_ini       = 'N'
      do_hess_fin       = 'N'
      do_hess_sp        = 'N'
      eig_0             = 0.0000010_dblprec
      is_afm            = 'N'

      !Parameters for energy interpolation along the MEP
      sample_num        = 500

      !Simulation parameters
      simid             = "_UppASD_"
      Mensemble         = 1
      tseed             = 1
      para_rng          = .false.
      llg               = 1
      nstep             = 1
      SDEalgh           = 1
      aunits            = 'N'
      perp              = 'N'

      !Tasks
      compensate_drift  = 0
      do_prnstruct      = 0
      do_prn_poscar     = 0
      do_prn_elk        = 0
      do_read_elk       = 0
      do_hoc_debug      = 0
      evolveout         = 0
      heisout           = 0
      mompar            = 0
      plotenergy        = 0
      do_sparse         = 'N'
      do_reduced        = 'N'

      !Measurement phase
      mode              = 'S'
      hfield            = (/zero,zero,zero/)
      mplambda1         = 0.050_dblprec
      mplambda2         = zero
      Temp              = zero
      delta_t           = 1.0e-16
      relaxtime         = 0.0e-16

      !Monte Carlo
      mcnstep           = 0
      mcavrg_step       = 0
      mcavrg_buff       = 0

      !Initial phase
      ipmode            = 'N'
      ipnphase          = 0
      ipmcnphase        = 1
      iphfield          = (/zero,zero,zero/)

      !Init mag
      mseed             = 1
      restartfile       = 'restart'
      initmag           = 4
      initexc           = 'N'
      initconc          = 0.0_dblprec
      initneigh         = 1
      initimp           = 0.0_dblprec
      theta0            = zero
      phi0              = zero
      roteul            = 0
      rotang            = (/zero,zero,zero/)
      initrotang        = 0.0_dblprec
      initpropvec       = (/zero,zero,zero/)
      initrotvec        = (/one,zero,zero/)
      do_mom_legacy     = 'N'

      !De-magnetization field
      demag             = 'N'
      demag1            = 'N'
      demag2            = 'N'
      demag3            = 'N'
      demagvol          = zero

      !Magnetic field pulse
      do_bpulse         = 0
      bpulsefile        = 'bpulsefile'
      locfield          = 'N'

      ! LSF
      conf_num          = 1
      gsconf_num        = 1
      lsf_metric        = 1
      do_lsf            = 'N'
      lsffile           = 'lsffile'
      lsf_interpolate   = 'Y'
      lsf_field         = 'T'
      lsf_window        = 0.050_dblprec

      !Measurements
      logsamp           = 'N'
      real_time_measure = 'N'
      do_spintemp       = 'N'
      spintemp_step     = 100

      !Random alloy data
      do_ralloy         = 0
      nchmax            = 1

      ! Random number transform
      ziggurat          = 'Y'
      rngpol            = 'N'

      ! GPU
      gpu_mode          = 0
      gpu_rng           = 0
      gpu_rng_seed      = 0

      ! I/O OVF
      prn_ovf           = 'N'
      read_ovf          = 'N'

   end subroutine set_input_defaults

   !---------------------------------------------------------------------------------
   ! SUBROUTINE: allocate_initmag
   !> @brief Allocate arrays for storing input moments
   !---------------------------------------------------------------------------------
   subroutine allocate_initmag(NA, Nchmax,conf_num,flag)

      implicit none

      integer,optional, intent(in) :: NA  !< Number of atoms in one cell
      integer,optional, intent(in) :: conf_num !< Number of configurations for LSF
      integer,optional, intent(in) :: Nchmax !< Max number of chemical components on each site in cell
      integer, intent(in) :: flag  !< Allocate or deallocate (1/-1)

      integer :: i_stat, i_all

      if(flag>0) then
         allocate(ammom_inp(NA,Nchmax,conf_num),stat=i_stat)
         call memocc(i_stat,product(shape(ammom_inp))*kind(ammom_inp),'ammom_inp','read_initmag')
         allocate(aemom_inp(3,NA,Nchmax,conf_num),stat=i_stat)
         call memocc(i_stat,product(shape(aemom_inp))*kind(aemom_inp),'aemom_inp','read_initmag')
         allocate(Landeg_ch(NA,Nchmax,conf_num),stat=i_stat)
         call memocc(i_stat,product(shape(Landeg_ch))*kind(Landeg_ch),'Landeg_ch','read_initmag')
         if (ind_mom_flag=='Y') then
            allocate(ind_mom(NA,Nchmax),stat=i_stat)
            call memocc(i_stat,product(shape(ind_mom))*kind(ind_mom),'ind_mom','read_initmag')
         endif
      else
         i_all=-product(shape(ammom_inp))*kind(ammom_inp)
         deallocate(ammom_inp,stat=i_stat)
         call memocc(i_stat,i_all,'ammom_inp','allocate_initmag')
         i_all=-product(shape(aemom_inp))*kind(aemom_inp)
         deallocate(aemom_inp,stat=i_stat)
         call memocc(i_stat,i_all,'aemom_inp','allocate_initmag')
         i_all=-product(shape(Landeg_ch))*kind(Landeg_ch)
         deallocate(Landeg_ch,stat=i_stat)
         call memocc(i_stat,i_all,'Landeg_ch','allocate_initmag')
         if (ind_mom_flag=='Y') then
            i_all=-product(shape(ind_mom))*kind(ind_mom)
            deallocate(ind_mom,stat=i_stat)
            call memocc(i_stat,i_all,'ind_mom','allocate_initmag')
         endif
      end if

   end subroutine allocate_initmag

   !---------------------------------------------------------------------------------
   ! SUBROUTINE: allocate_chemicalinput
   !> @brief Allocate arrays for random alloys
   !---------------------------------------------------------------------------------
   subroutine allocate_chemicalinput(NA,Nchmax,flag)
      implicit none

      integer, optional, intent(in) :: NA  !< Number of atoms in one cell
      integer, optional, intent(in) :: Nchmax !< Max number of chemical components on each site in cell
      integer, intent(in) :: flag !< Allocate or deallocate (1/-1)

      integer :: i_all, i_stat

      if(flag>0) then
         allocate(Nch(NA),stat=i_stat)
         call memocc(i_stat,product(shape(Nch))*kind(Nch),'Nch','allocate_chemicalinput')
         allocate(achtype_ch(NA,Nchmax),stat=i_stat)
         call memocc(i_stat,product(shape(achtype_ch))*kind(achtype_ch),'achtype_ch','allocate_chemicalinput')
         allocate(chconc(NA,Nchmax),stat=i_stat)
         call memocc(i_stat,product(shape(chconc))*kind(chconc),'chconc','allocate_chemicalinput')
      else
         i_all=-product(shape(Nch))*kind(Nch)
         deallocate(Nch,stat=i_stat)
         call memocc(i_stat,i_all,'Nch','allocate_chemicalinput')
         if(allocated(achtype_ch)) then
            i_all=-product(shape(achtype_ch))*kind(achtype_ch)
            deallocate(achtype_ch,stat=i_stat)
            call memocc(i_stat,i_all,'achtype_ch','allocate_chemicalinput')
         end if
         i_all=-product(shape(chconc))*kind(chconc)
         deallocate(chconc,stat=i_stat)
         call memocc(i_stat,i_all,'chconc','allocate_chemicalinput')
      end if

   end subroutine allocate_chemicalinput

   !---------------------------------------------------------------------------------
   ! SUBROUTINE: reshape_hamiltonianinput
   !> @brief Container routine to allocate/deallocate temporal cechange interactions
   !>  from the input
   !---------------------------------------------------------------------------------
   subroutine reshape_hamiltonianinput()

      implicit none

      real(dblprec), dimension(:,:,:,:), allocatable :: jc_tmp !< Exchange couplings
      real(dblprec), dimension(:,:,:), allocatable :: redcoord_tmp !< Coordinates for Heisenberg exchange couplings

      integer :: i,j,k,l,i_stat, i_all

      allocate(jc_tmp(NT,max_no_shells,Nchmax,Nchmax),stat=i_stat)
      call memocc(i_stat,product(shape(jc_tmp))*kind(jc_tmp),'jc_tmp','reshape_hamiltonian')
      allocate(redcoord_tmp(NT,max_no_shells,3),stat=i_stat)
      call memocc(i_stat,product(shape(redcoord_tmp))*kind(redcoord_tmp),'redcoord_tmp','reshape_hamiltonian')

      do l=1,Nchmax
         do k=1,Nchmax
            do j=1,max_no_shells
               do i=1,NT
                  jc_tmp(i,j,k,l)=jc(i,j,k,l,1)
               end do
            end do
         end do
      end do

      do k=1,3
         do j=1,max_no_shells
            do i=1,NT
               redcoord_tmp(i,j,k)=redcoord(i,j,k)
            end do
         end do
      end do

      i_all=-product(shape(jc))*kind(jc)
      deallocate(jc,stat=i_stat)
      call memocc(i_stat,i_all,'jc','reshape_hamiltonian')
      i_all=-product(shape(redcoord))*kind(redcoord)
      deallocate(redcoord,stat=i_stat)
      call memocc(i_stat,i_all,'redcoord','reshape_hamiltonian')

      allocate(jc(NT,max_no_shells,Nchmax,Nchmax,1),stat=i_stat)
      call memocc(i_stat,product(shape(jc))*kind(jc),'jc','reshape_hamiltonian')
      allocate(redcoord(NT,max_no_shells,3),stat=i_stat)
      call memocc(i_stat,product(shape(redcoord))*kind(redcoord),'redcoord','reshape_hamiltonian')

      jc(:,:,:,:,1)=jc_tmp
      redcoord=redcoord_tmp

      i_all=-product(shape(jc_tmp))*kind(jc_tmp)
      deallocate(jc_tmp,stat=i_stat)
      call memocc(i_stat,i_all,'jc_tmp','reshape_hamiltonian')
      i_all=-product(shape(redcoord_tmp))*kind(redcoord_tmp)
      deallocate(redcoord_tmp,stat=i_stat)
      call memocc(i_stat,i_all,'redcoord_tmp','reshape_hamiltonian')

   end subroutine reshape_hamiltonianinput


end module InputData
