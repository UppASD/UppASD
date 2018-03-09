!> Optimiztion stuff ....
!> @author
!! Johan Venemalm, Thomas Nystrand, Anders Bergman
!> @copyright
!! Copyright (C) 2008-2018 UppASD group
!! This file is distributed under the terms of the
!! GNU General Public License. 
!! See http://www.gnu.org/copyleft/gpl.txt
module optimizationRoutines
   use Parameters
   use Profiling
   use NeighbourMap, only : setup_nm
   use ErrorHandling

   implicit none

   !Optimization routines variables
   logical :: OPT_flag               !< Optimization flag
   integer :: OPT_rebuild_time       !< Time between consecutive rebuilds for optimization data structures
   logical ::  print_atom_class_flag !< Print cores flag - necessary in order to run the Python viewrestart2 script
   real(dblprec) :: cos_thr    !< Cosine similarity threshold
   integer :: larmor_numrev    !< Number of time steps per revolution
   real(dblprec) :: larmor_thr !< Larmor frequency threshold
   logical :: adaptive_time_flag, ip_adaptive_time_flag   !< Adaptive time stepping flag
   integer :: adapt_time_interval, ip_adapt_time_interval !< Interval between adaptive time step checkings

   ! To be used in conjunction with spin correlation only
   integer :: adapt_to_sc_ratio !< Number of spin correlation sampling periods between consecutive recomputations of delta_t

   logical :: totalsimtime_flag

   ! Define total simulation time for the adaptive scheme (to be used in inputHandler; it is not defaulted here in inputData)
   real(dblprec) :: totalsimtime
   real(dblprec), dimension(:), allocatable :: ip_totalsimtime

   integer, dimension(:,:), allocatable :: unitCellType    !< Three types exists:
   !< : 0  for noise cell
   !< : > 0 for core cell (index corresponds to position in constellation array)
   !< : < 0 for boundary cell

   integer, dimension(:),allocatable :: cellPosNumNeigh    !< Maximum number of neighbours at each position in the unit cell

   real(dblprec), dimension(:,:,:), allocatable :: constellations      !< Saved fixed unit cell configurations,
   !< These represent a configuration present
   !< in the domain with same configuration
   !< of unit cells in the neighborhood

   integer, dimension(:,:,:), allocatable :: constellationsNeighType   !< Every constellation atom's neighbour
   !< type atoms.  This will tell which values
   !< in the constellation to look for
   !< during the applyhamiltionian step

   real(dblprec), dimension(:,:,:), allocatable :: constellationsUnitVec  !< Normalized constellation matrix
   !< It is used for efficient cosinus comparisons in evolve step
   real(dblprec), dimension(:,:,:), allocatable :: constellationsUnitVec2
   real(dblprec), dimension(:,:), allocatable :: constellationsMag       !< Magnitude of magnetic moments
   !< for the constellation array
   real(dblprec), dimension(:,:,:), allocatable :: constlNCoup           !< Couplings for saved constellations

   integer, dimension(:), allocatable :: maxNoConstl ! Number of existing entries
   ! in for each ensemble in the constellation matrix
   integer :: max_no_constellations                  ! Max number of constellations for all ensembles
   logical :: OPT_ON = .false.                       ! Stop or start optimization during execution

   !    private :: printMatrix
   public  :: buildOptimizationRegions
   public  :: allocateOptimizationStuff
   public  :: smodeulermpt_constl
   public  :: modeulermpf_constl

   ! @todo THESE should be put in a type (struct)
   ! integer, intent(in) :: max_no_constellations ! The maximum (global) length of the constellation matrix
   ! integer, dimension(:), intent(in) :: maxNoConstl
   ! integer, dimension(:,:), intent(in) :: unitCellType ! Array of constellation id and classification (core, boundary, or noise) per atom
   ! real(dblprec), dimension(:,:,:), intent(in) :: constlNCoup
   ! real(dblprec), dimension(:,:,:), intent(in) :: constellations
   ! logical, intent(in) :: OPT_flag
   ! integer, dimension(:,:,:), intent(in) :: constellationsNeighType
contains

   !------------------------------------------------------------------------------------!
   subroutine buildOptimizationRegions(na,natom,Mensemble,max_no_neigh,atype,emom,mmom, &
         ncoup,nlist,nlistsize,cellPosNumNeigh,cos_thr)

      implicit none

      integer, intent(in) :: Na                !< Number of atoms in one cell
      integer, intent(in) :: Natom             !< Number of atoms in system
      integer, intent(in) :: Mensemble         !< Number of ensembles
      integer, intent(in) :: max_no_neigh      !< Calculated maximum of neighbours for exchange

      integer, dimension(Natom), intent(in) :: atype                      !< Type of atom
      real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: emom     !< Current unit moment vector
      real(dblprec), dimension(Natom,Mensemble), intent(in) :: mmom       !< Magnitude of magnetic moments
      real(dblprec), dimension(max_no_neigh,Natom), intent(in) :: ncoup   !< Heisenberg exchange couplings
      integer, dimension(max_no_neigh,Natom), intent(in) :: nlist         !< Neighbour list for Heisenberg exchange couplings
      integer, dimension(Natom), intent(in) :: nlistsize                  !< Size of neighbour list for Heisenberg exchange couplings
      integer, dimension(:),intent(in) :: cellPosNumNeigh                 !< Maximum number of neighbours at each position in the unit cell
      real(dblprec), intent(in) :: cos_thr                                            !< Cosine similarity threshold

      ! Private stuff
      real(dblprec), dimension(3,natom/na,Mensemble) :: tmp_constellations   ! Temporary holder for constellations
      integer, dimension(natom/na,Mensemble) :: atomConstlEx                 ! Holder for example atom index using constellation
      logical, dimension(na) :: tmp_diff                                     ! Temporary difference used in calculations
      logical :: fix_bool                                                    ! Used for internal checks
      integer :: neigh_index,h,i,j,k,naE,found_neigh,const_exists            ! Internal variables
      integer :: i_all,i_stat                                                ! Profiling related
      integer, dimension(Mensemble) :: nOfCoresFound

        call ErrorHandling_missing('Coarse-graining')
  
   end subroutine buildOptimizationRegions


   !DEC$ ATTRIBUTES FORCEINLINE :: invalidationCheck
   subroutine invalidationCheck(i, j, nlist, nlistsize, constellationsUnitVec, unitCellType, cos_thr, emom2)

      ! Optimization used variables
      integer, intent(in) :: i,j
      integer, dimension(:,:), intent(in) :: nlist ! Neighbour list for Heisenberg exchange couplings
      integer, dimension(:), intent(in) :: nlistsize ! Size of neighbour list for Heisenberg exchange couplings
      real(dblprec), dimension(:,:,:), intent(in) :: constellationsUnitVec ! Normalized magnetic moments for each atom and each constellation
      integer, dimension(:,:), intent(inout) :: unitCellType ! Type of each atom (core, boundary, or noise)
      real(dblprec), intent(in) :: cos_thr ! Cosine threshold value
      real(dblprec), dimension(:,:,:), intent(inout) :: emom2  !< Final (or temporary) unit moment vector

      ! Local optimization variables
      integer :: temp_index,k

        call ErrorHandling_missing('Coarse-graining')
  
   end subroutine invalidationCheck

   ! Allocation is not done for unitCellType since this must
   !   be done after builtOptimizationRegions has been called
   !------------------------------------------------------------------------------------!
   subroutine allocateOptimizationStuff(na,natom,Mensemble,ALLOC_flag)
      implicit none

      integer, intent(in) :: Na               !< Number of atoms in one cell
      integer, intent(in) :: Natom            !< Number of atoms in system
      integer, intent(in) :: Mensemble        !< Number of ensembles
      logical,intent(in)  :: ALLOC_flag       !< true for allocation, false for deallocation
      integer :: i_all,i_stat

      if (ALLOC_flag) then
         allocate(cellPosNumNeigh(Na),stat=i_stat)
         call memocc(i_stat,product(shape(cellPosNumNeigh))*kind(cellPosNumNeigh),'cellPosNumNeigh','allocateOptimizationStuff')

         allocate(maxNoConstl(Mensemble),stat=i_stat)
         call memocc(i_stat,product(shape(maxNoConstl))*kind(maxNoConstl),'maxNoConstl','allocateOptimizationStuff')

         allocate(unitCellType(natom,Mensemble),stat=i_stat)
         call memocc(i_stat,product(shape(unitCellType))*kind(unitCellType),'unitCellType','allocateOptimizationStuff')
      else
         i_all=-product(shape(cellPosNumNeigh))*kind(cellPosNumNeigh)
         deallocate(cellPosNumNeigh,stat=i_stat)
         call memocc(i_stat,i_all,'cellPosNumNeigh','allocateOptimizationStuff')

         i_all=-product(shape(constellations))*kind(constellations)
         deallocate(constellations,stat=i_stat)
         call memocc(i_stat,i_all,'constellations','allocateOptimizationStuff')

         i_all=-product(shape(constellationsNeighType))*kind(constellationsNeighType)
         deallocate(constellationsNeighType,stat=i_stat)
         call memocc(i_stat,i_all,'constellationsNeighType','allocateOptimizationStuff')

         i_all=-product(shape(constellationsUnitVec))*kind(constellationsUnitVec)
         deallocate(constellationsUnitVec,stat=i_stat)
         call memocc(i_stat,i_all,'constellationsUnitVec','allocateOptimizationStuff')

         i_all=-product(shape(constellationsUnitVec2))*kind(constellationsUnitVec2)
         deallocate(constellationsUnitVec2,stat=i_stat)
         call memocc(i_stat,i_all,'constellationsUnitVec2','buildOptRegions')

         i_all=-product(shape(unitCellType))*kind(unitCellType)
         deallocate(unitCellType,stat=i_stat)
         call memocc(i_stat,i_all,'unitCellType','allocateOptimizationStuff')

         i_all=-product(shape(maxNoConstl))*kind(maxNoConstl)
         deallocate(maxNoConstl,stat=i_stat)
         call memocc(i_stat,i_all,'maxNoConstl','allocateOptimizationStuff')

         i_all=-product(shape(constlNCoup))*kind(constlNCoup)
         deallocate(constlNCoup,stat=i_stat)
         call memocc(i_stat,i_all,'constlNCoup','allocateOptimizationStuff')

         i_all=-product(shape(constellationsMag))*kind(constellationsMag)
         deallocate(constellationsMag,stat=i_stat)
         call memocc(i_stat,i_all,'constellationsMag','buildOptRegions')

      end if

   end subroutine allocateOptimizationStuff



   ! Function used while debugging
   !------------------------------------------------------------------------------------!
   !subroutine printMatrix(matrix,n,m,o)
   !    integer,intent(in) :: n,m,o
   !    real(dblprec), dimension(n,m,o),intent(in) :: matrix
   !    integer :: i
   !    do i=1,m
   !        write(*,*) matrix(:,i,:)
   !    end do
   !end subroutine printMatrix


   ! Routine used for finding the maximum number of neighbours for each position in the unit cell
   subroutine getCellPosNumNeigh(Natom, NA, nlistsize, cellPosNumNeigh)

      implicit none

      integer, intent(in)                   :: Natom            !< Number of atoms in system
      integer, intent(in)                   :: NA               !< Number of atoms in one cell
      integer, dimension(NA), intent(out)   :: cellPosNumNeigh  !< Maximum number of neighbours for each position in the unit cell
      integer, dimension(Natom), intent(in) :: nlistsize        !< Size of neighbour list for Heisenberg exchange couplings

      ! Local variables
      integer :: i

        cellPosNumNeigh=0
        call ErrorHandling_missing('Coarse-graining')
  
   end subroutine getCellPosNumNeigh

   ! First evolve step for constellations
   !------------------------------------------------------------------------------------!
   subroutine smodeulermpt_constl(Natom, Mensemble, Landeg,bn, lambda1_array, beff, emom, emom2, emomM, mmom, deltat)
      use Constants

      implicit none

      integer, intent(in) :: Natom !< Number of atoms in system
      integer, intent(in) :: Mensemble !< Number of ensembles
      real(dblprec), dimension(Natom), intent(in) :: Landeg  !< Gyromagnetic ratio
      real(dblprec), intent(in) :: bn !< Scaling factor for LLG equation (parameter=1.0d0)
      real(dblprec), dimension(Natom), intent(in) :: lambda1_array !< Damping parameter
      real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: beff !< Total effective field from application of Hamiltonian
      real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: emom   !< Current unit moment vector
      real(dblprec), dimension(3,Natom,Mensemble), intent(inout) :: emom2  !< Final (or temporary) unit moment vector
      real(dblprec), dimension(3,Natom,Mensemble), intent(inout) :: emomM  !< Current magnetic moment vector
      real(dblprec), dimension(Natom,Mensemble), intent(in) :: mmom !< Magnitude of magnetic moments
      real(dblprec), intent(in) :: deltat !< Time step

      ! deterministic variables
      real(dblprec) :: a1x, a1y, a1z
      real(dblprec) :: b1x, b1y, b1z

      ! auxilary variables
      real(dblprec) :: Ax, Ay, Az, detAi
      real(dblprec) :: a2x, a2y, a2z
      !
      ! time steps
      real(dblprec) :: dt,sqrtdt
      real(dblprec) :: dtg, sqrtdtg

      !
      ! spins
      real(dblprec)  :: etx, ety, etz
      real(dblprec)  :: e1x, e1y, e1z

      ! ... Local variables ...
      integer :: i, j, ij
      !

        call ErrorHandling_missing('Coarse-graining')
  
   end subroutine smodeulermpt_constl




   !> Second evolve step of constellations
   !------------------------------------------------------------------------------------!
   subroutine modeulermpf_constl(Natom, Mensemble, Landeg,bn, lambda1_array, beff, emom, emom2, deltat)
      use Constants

      implicit none

      integer, intent(in) :: Natom !< Number of atoms in system
      integer, intent(in) :: Mensemble !< Number of ensembles
      real(dblprec), dimension(Natom), intent(in) :: Landeg  !< Gyromagnetic ratio
      real(dblprec), intent(in) :: bn !< Scaling factor for LLG equation (parameter=1.0d0)
      real(dblprec), dimension(Natom), intent(in) :: lambda1_array !< Damping parameter
      real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: beff !< Total effective field from application of Hamiltonian
      real(dblprec), dimension(3,Natom,Mensemble), intent(inout) :: emom   !< Current unit moment vector
      real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: emom2  !< Final (or temporary) unit moment vector
      real(dblprec), intent(in) :: deltat !< Time step

      ! deterministic variables
      real(dblprec) :: a1x, a1y, a1z
      real(dblprec) :: b1x, b1y, b1z
      !

      !
      ! auxilary variables
      real(dblprec) :: Ax, Ay, Az, detAi
      real(dblprec) :: a2x, a2y, a2z
      !
      ! time steps
      real(dblprec) :: dt,sqrtdt
      real(dblprec) :: dtg, sqrtdtg
      !
      ! spins
      real(dblprec) :: etpx, etpy, etpz
      real(dblprec) :: e1x, e1y, e1z
      !
      ! ... Local variables ...
      integer :: i, j, ij
      !
        call ErrorHandling_missing('Coarse-graining')
  
   end subroutine modeulermpf_constl


   !---------------------------------------------------------------------------
   !> @brief
   !> Read input parameters.
   !
   !> @author
   !> Anders Bergman
   !---------------------------------------------------------------------------
   subroutine read_parameters_opt(ifile)
      use FileParser

      implicit none

      ! ... Formal Arguments ...
      integer, intent(in) :: ifile   !< File to read from
      !
      ! ... Local Variables ...
      character(len=50) :: keyword
      integer :: rd_len, i_err, i_errb
      logical :: comment


      character(len=1) :: OPT_flag_str, adapt_flag_str, ip_adapt_flag_str, OPT_printcores_flag_str

        call ErrorHandling_missing('Coarse-graining')
  
   end subroutine read_parameters_opt

   subroutine set_opt_defaults(delta_t)
      !
      implicit none
      !
      real(dblprec), intent(in)  :: delta_t               !< Time step
      !
      ! Optimization routines
      OPT_flag               = .false.
      OPT_rebuild_time       = 500
      print_atom_class_flag  = .false.
      cos_thr                = 0.99d0
      adaptive_time_flag     = .false.
      ip_adaptive_time_flag  = .false.
      adapt_time_interval    = 100
      ip_adapt_time_interval = 100
      adapt_to_sc_ratio      = 1
      larmor_numrev          = 20
      larmor_thr             = 0.15d0
      totalsimtime           = delta_t
      !
   end subroutine set_opt_defaults

end module optimizationRoutines
