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

## Autocorrelation

## Energy

## Polarization

## Magnetic field

## Ionic displacement

## Ionic velocity
