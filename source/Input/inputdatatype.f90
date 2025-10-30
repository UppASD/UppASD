!------------------------------------------------------------------------------------
! MODULE: InputDataType
!> @brief New input data type for input from file
!> @copyright
!> GNU Public License.
!------------------------------------------------------------------------------------
module InputDataType

   use Parameters

   implicit none
   
   ! Custom type for Hamiltonian input data
   type ham_inp_t
      !sequence
      !---------------------------------------------------------------------------------
      ! Exchange data
      !---------------------------------------------------------------------------------
      integer :: do_jtensor                     !< Use SKKR style exchange tensor (0=off, 1=on) --OLD: 2=with biquadratic exchange
      logical :: calc_jtensor                   !< Calculate or read tensor from input (T/F)
      logical :: map_multiple                   !< Allow for multiple couplings between atoms i and j (use only for small cells)
      integer :: max_no_shells                  !< Actual maximum number of shells
      integer, dimension(:), allocatable :: nn  !< Number of neighbour shells
      integer, dimension(:,:), allocatable :: nntype                  !< Type for exchange target atom
      real(dblprec), dimension(:,:,:), allocatable :: redcoord       !< Coordinates for Heisenberg exchange couplings
      real(dblprec), dimension(:,:,:,:,:), allocatable :: jc         !< Exchange couplings
      real(dblprec), dimension(:,:,:,:,:), allocatable :: jcD        !< Exchange couplings (DLM)
      real(dblprec), dimension(:,:,:,:,:,:), allocatable :: jc_tens  !< Tensorial exchange couplings (SKKR)
      character(len=1) :: exc_inter                         !< Interpolate Jij between FM and DLM states (Y(N)
      !character(len=30) :: pre_jfile                        !< File name for exchange couplings
      !character(len=30), dimension(:), allocatable :: jfile    !< File name for exchange couplings
      !character(len=30), dimension(:), allocatable :: jfileD   !< File name for exchange couplings (DLM)
      real(dblprec) :: jij_scale                               !< Rescale Jij couplings manually
      logical ::  ea_model                                     !< Randomize exchange couplings for Edwards-Anderson model T,(F)
      real(dblprec) ::  ea_sigma                               !< Standard deviation for Edwards-Anderson randomization
      character(len=1) :: ea_algo                              !< Algoritm for Edwards-Anderson randomization
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
      ! Dzyaloshinskii-Moriya data
      !---------------------------------------------------------------------------------
      integer :: do_dm                                      !< Add Dzyaloshinskii-Moriya (DM) term to Hamiltonian (0/1)
      integer :: max_no_dmshells                            !< Actual maximum number of shells for DM interactions
      character(len=35) :: dmfile                           !< File name for Dzyaloshinskii-Moriya data
      real(dblprec) :: dm_scale                             !< Rescale Dij couplings manually
      logical ::  rdm_model                                 !< Randomize Dij couplings
      real(dblprec) ::  rdm_sigma                           !< Standard deviation for random Dij couplings
      character(len=1) :: rdm_algo                          !< Algoritm for random Dij couplings
      integer, dimension(:), allocatable :: dm_nn           !< No. shells of neighbours for DM
      real(dblprec), dimension(:,:,:), allocatable :: dm_redcoord    !< Neighbour vectors for DM
      real(dblprec), dimension(:,:,:,:,:), allocatable :: dm_inpvect !< Neighbour vectors for DM
      !---------------------------------------------------------------------------------
      ! Symmetric anisotropic data
      !---------------------------------------------------------------------------------
      integer :: do_sa                                      !< Add Symmetric anisotropic (SA) term to Hamiltonian (0/1)
      integer :: max_no_sashells                            !< Actual maximum number of shells for SA interactions
      character(len=35) :: safile                           !< File name for Symmetric anisotropic data
      real(dblprec) :: sa_scale                             !< Rescale Cij couplings manually
      integer, dimension(:), allocatable :: sa_nn           !< No. shells of neighbours for SA
      real(dblprec), dimension(:,:,:), allocatable :: sa_redcoord    !< Neighbour vectors for SA
      real(dblprec), dimension(:,:,:,:,:), allocatable :: sa_inpvect !< Neighbour vectors for SA
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
      ! Scalar four-site data
      !---------------------------------------------------------------------------------
      integer :: do_fourx                                       !< Add scalar four-site exchange (fourx) term to Hamiltonian (0/1)
      integer :: max_no_fourxshells                             !< Actual maximum number of shells for fourx interactions
      character(len=35) :: fourxfile                            !< File name for Pseudo-Dipolar data
      integer, dimension(:), allocatable :: fourx_nn            !< No. shells of neighbours forfourx 
      real(dblprec), dimension(:,:,:,:), allocatable :: fourx_inpval     !< Interactions for fourx 
      real(dblprec), dimension(:,:,:,:), allocatable :: fourx_redcoord   !< Neighbour vectors for fourx 
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
      real(dblprec), dimension(:,:,:,:,:), allocatable :: biqdm_inpvect    !< BIQDM interaction vectors
      !---------------------------------------------------------------------------------
      ! Biquadratic exchange data
      !---------------------------------------------------------------------------------
      integer :: do_bq                                   !< Add biquadratic (BQ) term to Hamiltonian (0/1)
      integer :: max_no_bqshells                         !< Actual maximum number of shells for BQ interactions
      character(len=35) :: bqfile                        !< File name for biquadratic data
      integer, dimension(:), allocatable :: bq_nn        !< No. shells of neighbours for BQ
      real(dblprec), dimension(:,:,:,:), allocatable :: jc_bq        !< Biquadratic exchange coupling
      real(dblprec), dimension(:,:,:), allocatable :: bq_redcoord    !< Neighbour vectors for BQ
      !---------------------------------------------------------------------------------
      ! Four-spin ring exchange data
      !---------------------------------------------------------------------------------
      integer :: do_ring                                 !< Add four-spin ring (4SR) term to Hamiltonian (0/1)
      integer :: max_no_ringshells                         !< Actual maximum number of shells for 4SR interactions
      character(len=35) :: ringfile                        !< File name for four-spin ring data
      integer, dimension(:), allocatable :: ring_nn        !< No. shells of neighbours for 4SR
      real(dblprec), dimension(:,:,:), allocatable :: ring_redcoord_ij    !< Neighbour vector rij for 4SR
      real(dblprec), dimension(:,:,:), allocatable :: ring_redcoord_ik    !< Neighbour vector rik for 4SR
      real(dblprec), dimension(:,:,:), allocatable :: ring_redcoord_il    !< Neighbour vector ril for 4SR
      real(dblprec), dimension(:), allocatable :: jc_ring  !<Four-spin ring couplings
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
   end type ham_inp_t

   public

end module InputDataType
