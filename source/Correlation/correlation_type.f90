module Correlation_type
   use Parameters
   use Profiling
   !
   implicit none
   !

   type corr_t
      !!! !Spin-spin correlation parameters from input file
      !!! character(len=1) :: do_sc        !< Perform spin-correlation sampling (Y/N/C)
      !!! character(len=1) :: do_sr        !< Perform spin-correlation sampling in real space direclty (Y/N)
      character(len=1) :: do_proj   !< Measure sublattice projection of S(q,w) (Y/N/C)
      character(len=1) :: do_projch   !< Measure chemical sublattice projection of S(q,w) (Y/N/C)
      character(len=1) :: do_sc_local_axis   !< Perform SQW along local quantization axis (Y/N)
      character(len=1) :: do_sc_dosonly !< Do not print s(q,w), only Magnon DOS (Y/N)
      character(len=1) :: do_sc_complex !< Print the complex values s(q,w) (Y/N)
      character(len=1) :: do_sc_tens   !< Print the tensorial values s(q,w) (Y/N)
      character(len=1) :: do_sc_list   !< Print the max locations of s(q,w) as list (Y/N)
      character(len=1) :: do_qt_traj   !< Measure time trajectory of S(q,t) (Y/N)
      character(len=1) :: do_connected !< Perform the connected part S(q,w)
      character(len=1) :: sc_average   !< Averaging of S(q,w): (F)ull, (E)ven, or (N)one

      character(len=2) :: label        !< Label for the correlation type (sc, uc)
      integer :: sc_sep                !< Separation between sampling of C(q)
      integer :: sc_nsamp              !< Number of C(q) samples
      integer :: sc_step               !< Separation between sampling steps of C(q,t)
      integer :: sc_nstep              !< Number of steps to sample C(q,t)
      integer :: sc_max_nstep          !< Max number of sampling opportunities of C(q,t)
      integer :: sc_naverages          !< Number of averages for S(q,w) averaging 

      real(dblprec) :: sc_emax         !< Energy range for S(q,w) in mRy, replaces `sc_step`
      real(dblprec) :: sc_eres         !< Energy resolution for S(q,w) in mRy, replaces `sc_nstep`

      !!! ! Working variables to perform the printing of the correlation
      integer :: sc_tidx
      integer :: sc_samp_done !< Flag to keep track of if S(q) sampling is done (for do_sc='C')
      integer :: sc_samp_done_sr !< Flag to keep track of if S(r) sampling is done (for do_sr='Y')
      
      complex(dblprec), dimension(:,:), allocatable :: m_k                ! Correlation in q G(k)
      complex(dblprec), dimension(:,:,:), allocatable :: m_k_proj       ! Correlation for G(k) (projected)
      complex(dblprec), dimension(:,:,:), allocatable :: m_k_projch     ! Correlation for G(k) (chemical)

      complex(dblprec), dimension(:,:), allocatable :: m_kt0           ! Temporary for g(q,t0)
      real(dblprec), dimension(:,:), allocatable :: corr_srt           ! Correlation in r G(r,t)
      real(dblprec), dimension(:,:), allocatable :: corr_sr            ! Correlation in r G(r)
      real(dblprec), dimension(:,:), allocatable :: corr_r             ! Correlation in r G(r) calculated directly

      complex(dblprec), dimension(:,:), allocatable :: corr_sA              ! Correlation in r G(r) for variable A
      complex(dblprec), dimension(:,:), allocatable :: corr_sB              ! Correlation in r G(r) for variable B

      real(dblprec), dimension(:,:,:,:), allocatable :: corr_k_proj     ! Correlation in q G(k) (sublattice)
      
      ! Change to a_kt and add b_kt
      complex(dblprec), dimension(:,:,:), allocatable   :: m_kt            ! Correlation for G(k,t)
      complex(dblprec), dimension(:,:,:,:), allocatable :: m_kt_proj       ! Correlation for G(k,t)
      complex(dblprec), dimension(:,:,:,:), allocatable :: m_kt_projch     ! Correlation for G(k,t)

      ! Change to a_kt and add b_kt
      complex(dblprec), dimension(:,:,:), allocatable   :: m_kw            ! Correlation for G(k,t)
      complex(dblprec), dimension(:,:,:,:), allocatable :: m_kw_proj       ! Correlation for G(k,t)
      complex(dblprec), dimension(:,:,:,:), allocatable :: m_kw_projch     ! Correlation for G(k,t)



      real(dblprec), dimension(:), allocatable :: deltat_corr ! Array storing delta_t for each sample
      real(dblprec), dimension(:), allocatable :: scstep_arr  ! Array storing sc_step for each sample

      real(dblprec), dimension(:,:), allocatable :: SA_axis      !< Internal array for finding the global direction of magnetization
      real(dblprec), dimension(:,:), allocatable :: SB_axis      !< Internal array for finding the global direction of magnetization

      real(dblprec), dimension(3) :: SA_avrg      !< Internal array for finding the global direction of magnetization
      real(dblprec), dimension(3) :: SB_avrg      !< Internal array for finding the global direction of magnetization

      real(dblprec), dimension(:,:,:), allocatable :: mavg_local_axis!< Internal array for finding the local average direction of magnetization
      real(dblprec), dimension(:,:,:,:), allocatable :: mavg_local_rotmat !< Internal array for finding rotatons to project the local average direction of magnetization
      real(dblprec), dimension(:,:,:), allocatable :: mort_axis    !< Orthogonal (and parallel ) components of magnetization
      real(dblprec), dimension(:,:,:), allocatable ::  m_loc      !< Array containing moments projected to parallel/perpendicular directions
      !!! !
      !!! !
      real(dblprec) :: sc_local_axis_mix  !< Mixing parameter for updating the local quantization axis. Rec. value: 0.1-0.5
      !!! !
      !!! real(dblprec), dimension(3) :: r_mid !< Local variable containing the center coordinate for the system. Used for FT in calc_gk
      integer :: gkt_flag    !< Flag to keep track of dynamic sampling
      integer :: gk_flag    !< Flag to keep track of dynamic sampling

      character :: prefix    !< Prefix for output files. s=spin, u=displacement, v=velocities
      integer :: cmode       !< Flag to determine correlation function mode. 1=single, 2=dual

      real(dblprec), dimension(:), allocatable :: time    !< Accumulated time 
      real(dblprec), dimension(:), allocatable :: dt      !< Timesteps per iterations (for adaptive stepping)

      integer :: nw                                       !< Number of frequencies to sample
      real(dblprec), dimension(:), allocatable :: w       !< Frequencies

      ! Real-space correlation
      integer, dimension(:,:), allocatable :: cr_list
      real(dblprec), dimension(:,:,:), allocatable :: cr_vect
      real(dblprec) :: cr_cut
      integer :: cr_nnmax
      integer :: cr_nhist
      integer :: cr_flag
      real(dblprec), dimension(:), allocatable :: cr_hist
      real(dblprec), dimension(:,:), allocatable :: cr_uniq_vect
      integer, dimension(:,:), allocatable :: cr_lut



      !!! private

   end type corr_t

end module Correlation_type
