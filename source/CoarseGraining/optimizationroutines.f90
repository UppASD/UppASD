!> Optimiztion stuff ....
module optimizationRoutines
   use Parameters
   use Profiling
   use NeighbourMap, only : setup_nm

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

   integer, dimension(:,:,:), allocatable :: constellationsNeighType   !< Every constellation atom´s neighbour
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

contains

   !------------------------------------------------------------------------------------!
   subroutine buildOptimizationRegions(na,natom,nHam,Mensemble,max_no_neigh,atype,emom,mmom, &
         ncoup,nlist,nlistsize,aham, cellPosNumNeigh,cos_thr)

      implicit none

      integer, intent(in) :: Na                !< Number of atoms in one cell
      integer, intent(in) :: Natom             !< Number of atoms in system
      integer, intent(in) :: nHam              !< Number of atoms in Hamiltonian
      integer, intent(in) :: Mensemble         !< Number of ensembles
      integer, intent(in) :: max_no_neigh      !< Calculated maximum of neighbours for exchange

      integer, dimension(Natom), intent(in) :: atype                      !< Type of atom
      real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: emom     !< Current unit moment vector
      real(dblprec), dimension(Natom,Mensemble), intent(in) :: mmom       !< Magnitude of magnetic moments
      real(dblprec), dimension(max_no_neigh,nHam), intent(in) :: ncoup   !< Heisenberg exchange couplings
      integer, dimension(max_no_neigh,Natom), intent(in) :: nlist         !< Neighbour list for Heisenberg exchange couplings
      integer, dimension(nHam), intent(in) :: nlistsize                  !< Size of neighbour list for Heisenberg exchange couplings
      integer, dimension(Natom), intent(in) :: aham                       !< Hamiltonian look-up table
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

      ! Initialize parameter values
      naE                    = na-1
      tmp_constellations     = 0
      unitCellType           = 0
      nOfCoresFound          = 0
      maxNoConstl            = 1
      max_no_constellations  = 0

      ensamble: do h=1,Mensemble
         ! First pass: find core points and mark the neighbours as boundary points
         ! Note that the structure unitCellType is initialized with zero entries
         ! Thus, all points are initially marked as noise points


         do i=1,natom,na
            fix_bool = .true.
            const_exists = 0
            found_neigh = 0

            ! Check if a constellation is previously saved with the same current configuration
            do j=1,maxNoConstl(h),na
               do k=0,naE
                  tmp_diff(k+1) = dot_product(tmp_constellations(:,j+k,h),emom(:,i+k,h)).ge.cos_thr
               end do
               if (ALL(tmp_diff(:))) then
                  const_exists = j
                  exit
               end if
            end do

            ! Loop over the unit cell
            do j=0,naE
               do k=1,nlistsize(i+j)
                  neigh_index = nlist(k,i+j)
                  ! Only atoms of same type must spin in same dir
                  if (atype(i+j).eq.atype(neigh_index)) then
                     ! Potential core neighborhood should only have neighboors corresponding to same constellation
                     found_neigh = abs(unitCellType(neigh_index,h))
                     ! Any (previous or current) found non-zero neighbors (boundaries or cores)
                     ! must belong to same constellation
                     if (found_neigh.gt.0) then
                        if (found_neigh.ne.const_exists) then
                           fix_bool = .false.
                           exit
                        end if
                     else
                        fix_bool = fix_bool.and.dot_product(emom(:,i+j,h),emom(:,neigh_index,h)).ge.cos_thr
                     end if
                  end if
               end do
               if (.not.fix_bool) then
                  exit
               end if
            end do

            ! If the found point is a fixed core point
            if (fix_bool) then
               nOfCoresFound(h) = nOfCoresFound(h)+na
               ! If there exists a previously saved constellation
               if (const_exists.gt.0) then
                  unitCellType(i:i+naE,h) = const_exists
                  ! If there did not exist a previous constellation (found point is an isolated fixed core point)
               else
                  tmp_constellations(:,maxNoConstl(h):maxNoConstl(h)+naE,h) = emom(:,i:i+naE,h)
                  unitCellType(i:i+naE,h) = maxNoConstl(h)
                  const_exists = maxNoConstl(h)
                  maxNoConstl(h) = maxNoConstl(h)+na
               end if

               ! Change all previous noise-marked atoms to neighbors for the found core point
               do j=0,naE
                  do k=1,nlistsize(i+j)
                     neigh_index = nlist(k,i+j)
                     if(unitCellType(neigh_index,h).eq.0) then
                        unitCellType(neigh_index,h) = -const_exists
                     end if
                  end do
               end do
            end if
         end do

         !call cpu_time(timeb)


         ! Update highest constellation numbers
         maxNoConstl(h)=maxNoConstl(h)-1
         if (maxNoConstl(h).gt.max_no_constellations) then
            max_no_constellations = maxNoConstl(h)
         end if

         ! Second pass: find all domain boundary points and mark these as boundary points.
         ! This will cause a few extra calculations in next time step since some boundary
         ! points might in reality be noise points.
         do i=1,natom
            if (nlistsize(i).lt.cellPosNumNeigh(mod(i-1,Na)+1)) then
               ! If the considered point is a core point, degrade it to a boundary point
               ! Core points are not allowed on the boundaries
               if (unitCellType(i,h).gt.0) then
                  nOfCoresFound(h) = nOfCoresFound(h)-1
                  unitCellType(i,h)=-unitCellType(i,h)
               end if
            end if
         end do

         ! Third pass: find example/representative core atoms
         ! (Used to figure out couplings and neighbours to be calculated in applyhamiltonian)
         tmp_constellations(:,:,h) = 0
         do k=1,max_no_constellations,na
            i=1
            do while(unitCellType(i,h).ne.k)
               i=i+1
               if(i.ge.natom) then
                  i=1
                  exit
               end if
            end do

            do j=0,naE
               atomConstlEx(k+j,h) = i+j
            end do

            ! Taking average constellation and normalizing it (avoid large changes in moments when rebuilding)
            i=0
            do j=1,natom,na
               if (unitCellType(j,h).eq.k) then
                  i=i+1
                  tmp_constellations(:,k:k+naE,h) = tmp_constellations(:,k:k+naE,h)+emom(:,j:j+naE,h)
               end if
            end do
            if (i.gt.1) then
               tmp_constellations(:,k:k+naE,h)=tmp_constellations(:,k:k+naE,h)/i
            end if
            do j=1,na
               tmp_constellations(:,j,h) = tmp_constellations(:,j,h)/norm2(tmp_constellations(:,j,h))
            end do
         end do
      end do ensamble

      ! Save memory by using a smaller constellation and ncoupling matrix
      ! First deallocate any potential previous allocation
      if (allocated(constellations)) then
         i_all=-product(shape(constellations))*kind(constellations)
         deallocate(constellations,stat=i_stat)
         call memocc(i_stat,i_all,'constellations','buildOptRegions')

         i_all=-product(shape(constellationsNeighType))*kind(constellationsNeighType)
         deallocate(constellationsNeighType,stat=i_stat)
         call memocc(i_stat,i_all,'constellationsNeighType','buildOptRegions')

         i_all=-product(shape(constlNCoup))*kind(constlNCoup)
         deallocate(constlNCoup,stat=i_stat)
         call memocc(i_stat,i_all,'constlNCoup','buildOptRegions')

         i_all=-product(shape(constellationsUnitVec))*kind(constellationsUnitVec)
         deallocate(constellationsUnitVec,stat=i_stat)
         call memocc(i_stat,i_all,'constellationsUnitVec','buildOptRegions')

         i_all=-product(shape(constellationsUnitVec2))*kind(constellationsUnitVec2)
         deallocate(constellationsUnitVec2,stat=i_stat)
         call memocc(i_stat,i_all,'constellationsUnitVec2','buildOptRegions')

         i_all=-product(shape(constellationsMag))*kind(constellationsMag)
         deallocate(constellationsMag,stat=i_stat)
         call memocc(i_stat,i_all,'constellationsMag','buildOptRegions')
      end if

      allocate(constellations(3,max_no_constellations,Mensemble),stat=i_stat)
      call memocc(i_stat,product(shape(constellations))*kind(constellations),'constellations','buildOptRegions')

      allocate(constellationsNeighType(max_no_neigh,max_no_constellations,Mensemble),stat=i_stat)
      call memocc(i_stat,product(shape(constellationsNeighType))*kind(constellationsNeighType),'constellationsNeighType','buildOptRegions')

      allocate(constlNCoup(max_no_neigh,max_no_constellations,Mensemble),stat=i_stat)
      call memocc(i_stat,product(shape(constlNCoup))*kind(constlNCoup),'constlNCoup','buildOptRegions')

      allocate(constellationsUnitVec(3,max_no_constellations,Mensemble),stat=i_stat)
      call memocc(i_stat,product(shape(constellationsUnitVec))*kind(constellationsUnitVec),'constellationsUnitVec','buildOptRegions')

      allocate(constellationsUnitVec2(3,max_no_constellations,Mensemble),stat=i_stat)
      call memocc(i_stat,product(shape(constellationsUnitVec2))*kind(constellationsUnitVec2),'constellationsUnitVec2','buildOptRegions')

      allocate(constellationsMag(max_no_constellations,Mensemble),stat=i_stat)
      call memocc(i_stat,product(shape(constellationsMag))*kind(constellationsMag),'constellationsMag','buildOptRegions')

      ! Setting up couplings and neighbours needed during applyhamiltonian for evaluation of constellations
      constellationsNeighType = 1
      constlNCoup = 0 ! Needed for applyhamiltonian
      do k=1,Mensemble
         do i=1,maxNoConstl(k)
            do j=1,max_no_neigh
               if (nlist(j,atomConstlEx(i,k)).ne.0) then
                  constellationsNeighType(j,i,k) = abs(unitCellType(nlist(j,atomConstlEx(i,k)),k))

                  ! Ugly setting the basic neigh type (just to prevent indexing errors in applyhamiltionian)
                  if (constellationsNeighType(j,i,k).eq.0) constellationsNeighType(j,i,k) = 1
               end if
            end do
            constellationsMag(i,k) = mmom(atomConstlEx(i,k),k)
            constlNCoup(:,i,k) = ncoup(:,aham(atomConstlEx(i,k)))
         end do
      end do

      ! Setting up constellations
      constellationsUnitVec(:,:,:) = tmp_constellations(:,1:max_no_constellations,:)

      do i=1,3
         constellations(i,:,:) = tmp_constellations(i,1:max_no_constellations,:)*constellationsMag(:,:)
      end do

      write(*,*) 'Number of core atoms: ',nOfCoresFound
      if (any(nOfCoresFound.lt.natom/3)) then
         if (OPT_flag) then
            write(*,*) 'Optimization Off'
            OPT_flag = .false.
         end if
      else if (.not.OPT_flag) then
         write(*,*) 'Optimization On'
         OPT_flag = .true.
      end if

      if(print_atom_class_flag) then
         open(ofileno,file='unitcelldata.out')
         do j=1,Mensemble
            do i=1,natom
               write(ofileno,110) i,unitCellType(i,j),nlistsize(i),emom(:,i,j)
            end do
         end do
         110 format(i8,i8,i8,8x,es16.8,es16.8,es16.8)
         close(ofileno)
      end if

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

      ! Validity testing to see if the moments due to randomness have spun out of control
      if (unitCellType(i,j).ne.0) then
         temp_index = abs(unitCellType(i,j))
         if(dot_product(constellationsUnitVec(:,temp_index,j), emom2(:,i,j)).lt.cos_thr) then
            unitCellType(i,j) = 0
            do k=1,nlistsize(i)
               if (unitCellType(nlist(k,i),j).gt.0) then
                  unitCellType(nlist(k,i),j) = -temp_index
               end if
            end do
         end if
      end if

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

      cellPosNumNeigh = 0 ! Initialize the array to zero

      do i=1,Natom
         if(cellPosNumNeigh(mod(i-1,Na)+1) < nlistsize(i) ) then
            cellPosNumNeigh(mod(i-1,Na)+1)= nlistsize(i)
         end if
      end do

   end subroutine getCellPosNumNeigh

   ! First evolve step for constellations
   !------------------------------------------------------------------------------------!
   subroutine smodeulermpt_constl(Natom, Mensemble, Landeg,bn, lambda1_array, beff, emom, emom2, emomM, mmom, deltat)
      use Constants

      implicit none

      integer, intent(in) :: Natom !< Number of atoms in system
      integer, intent(in) :: Mensemble !< Number of ensembles
      real(dblprec), dimension(Natom), intent(in) :: Landeg  !< Gyromagnetic ratio
      real(dblprec), intent(in) :: bn !< Scaling factor for LLG equation (parameter=1.0_dblprec)
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

      dt=deltat*bn*gama !dimm. less time
      sqrtdt=sqrt(dt)

      !$omp parallel do default(shared) schedule(static) &
      !$omp private(ij,i,j,e1x,e1y,e1z,etx,ety,etz,&
      !$omp a1x,a1y,a1z,b1x,b1y,b1z,Ax,Ay,Az,detAi,a2x,a2y,a2z,dtg,sqrtdtg)
      do ij=1,Natom*Mensemble

         i=mod(ij-1,Natom)+1
         j=int((ij-1)/Natom)+1
         dtg=dt*Landeg(i)
         sqrtdtg=sqrtdt*Landeg(i)
         e1x=emom(1,i,j) !emom is value previous timestep
         e1y=emom(2,i,j) !in step t emom=emom2 due to copym after step f
         e1z=emom(3,i,j)
         b1x=beff(1,i,j) ! effective field approximated with et
         b1y=beff(2,i,j)
         b1z=beff(3,i,j)

         ! a1 = -b1 - lambda*(e1 cross b1)
         a1x=-b1x-lambda1_array(i)*(e1y*b1z-e1z*b1y)
         a1y=-b1y-lambda1_array(i)*(e1z*b1x-e1x*b1z)
         a1z=-b1z-lambda1_array(i)*(e1x*b1y-e1y*b1x)

         !
         ! semi-implicitness midpoint requires solution of linear system:
         ! A*e2 = At*e1, At=transpose(A) => e2=inv(A)*At*e1
         ! A = I + skew(dt*a1/2 + sqrt(dt)*s1/2)
         ! write A*e2=a2, a2=At*e1 => e2=inv(A)*a2
         ! Ax,Ay,Az off-diagonal components of A
         ! solve with Cramers´ rule => define detAi=1/determinant(A)
         !
         Ax=0.5_dblprec*dtg*a1x
         Ay=0.5_dblprec*dtg*a1y
         Az=0.5_dblprec*dtg*a1z

         detAi=1.0_dblprec/(1.0_dblprec+Ax*Ax+Ay*Ay+Az*Az)
         !
         a2x=e1x+e1y*Az-e1z*Ay
         a2y=e1y+e1z*Ax-e1x*Az
         a2z=e1z+e1x*Ay-e1y*Ax
         !
         etx=a2x*(1+Ax*Ax)+a2y*(Ax*Ay+Az)+a2z*(Ax*Az-Ay)
         ety=a2x*(Ay*Ax-Az)+a2y*(1+Ay*Ay)+a2z*(Ay*Az+Ax)
         etz=a2x*(Az*Ax+Ay)+a2y*(Az*Ay-Ax)+a2z*(1+Az*Az)

         ! now use et for writing et´=(e1+et)/2 in emom2
         etx=0.5_dblprec*(e1x+etx*detAi)
         ety=0.5_dblprec*(e1y+ety*detAi)
         etz=0.5_dblprec*(e1z+etz*detAi)
         !
         ! write et´=(e1+et)/2 in emom2
         emom2(1,i,j)=etx
         emom2(2,i,j)=ety
         emom2(3,i,j)=etz

         ! write new emomM
         ! effective field for step f approximated with et
         ! no new variables (e.g emomt) to reduce memory requirements
         emomM(1,i,j)=etx*mmom(i,j)
         emomM(2,i,j)=ety*mmom(i,j)
         emomM(3,i,j)=etz*mmom(i,j)
      end do
      !$omp end parallel do
      return
   end subroutine smodeulermpt_constl




   !> Second evolve step of constellations
   !------------------------------------------------------------------------------------!
   subroutine modeulermpf_constl(Natom, Mensemble, Landeg,bn, lambda1_array, beff, emom, emom2, deltat)
      use Constants

      implicit none

      integer, intent(in) :: Natom !< Number of atoms in system
      integer, intent(in) :: Mensemble !< Number of ensembles
      real(dblprec), dimension(Natom), intent(in) :: Landeg  !< Gyromagnetic ratio
      real(dblprec), intent(in) :: bn !< Scaling factor for LLG equation (parameter=1.0_dblprec)
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
      ! scale dt (from evolve.f90)
      dt=deltat*bn*gama !dimm. less time
      sqrtdt=sqrt(dt)
      !
      !$omp parallel do default(shared) schedule(static) &
      !$omp private(ij,i,j,e1x,e1y,e1z,etpx,etpy,etpz,&
      !$omp a1z,b1x,b1y,b1z,Ax,Ay,Az,detAi,a2x,a2y,a2z,dtg,sqrtdtg)
      do ij=1,Natom*Mensemble
         i=mod(ij-1,Natom)+1
         j=int((ij-1)/Natom)+1
         dtg=dt*Landeg(i)
         sqrtdtg=sqrtdt*Landeg(i)
         e1x=emom(1,i,j) ! load e1 back from emom
         e1y=emom(2,i,j) ! emom unchanged by step t
         e1z=emom(3,i,j)
         etpx=emom2(1,i,j) ! load etp=et´ back from emom2
         etpy=emom2(2,i,j)
         etpz=emom2(3,i,j)
         b1x=beff(1,i,j) ! effective field approximated with e1
         b1y=beff(2,i,j)
         b1z=beff(3,i,j)

         !
         ! a1 = -b1 - lambda*(et cross b1)
         a1x=-b1x-lambda1_array(i)*(etpy*b1z-etpz*b1y)
         a1y=-b1y-lambda1_array(i)*(etpz*b1x-etpx*b1z)
         a1z=-b1z-lambda1_array(i)*(etpx*b1y-etpy*b1x)
         !
         !
         ! semi-implicitness midpoint requires solution of linear system:
         ! A*e2 = At*e1, At=transpose(A) => e2=inv(A)*At*e1
         ! A = I + skew(dt*a1/2 + sqrt(dt)*s1/2)
         ! write A*e2=a2, a2=At*e1 => e2=inv(A)*a2
         ! Ax,Ay,Az off-diagonal components of A
         ! solve with Cramers´ rule => define detAi=1/determinant(A)
         !
         Ax=0.5_dblprec*dtg*a1x
         Ay=0.5_dblprec*dtg*a1y
         Az=0.5_dblprec*dtg*a1z
         detAi=1.0_dblprec/(1+Ax*Ax+Ay*Ay+Az*Az)
         !
         a2x=e1x+e1y*Az-e1z*Ay
         a2y=e1y+e1z*Ax-e1x*Az
         a2z=e1z+e1x*Ay-e1y*Ax

         ! now use etp simply to write emom2
         etpx=a2x*(1+Ax*Ax)+a2y*(Ax*Ay+Az)+a2z*(Ax*Az-Ay)
         etpy=a2x*(Ay*Ax-Az)+a2y*(1+Ay*Ay)+a2z*(Ay*Az+Ax)
         etpz=a2x*(Az*Ax+Ay)+a2y*(Az*Ay-Ax)+a2z*(1+Az*Az)
         !
         emom(1,i,j)=etpx*detAi
         emom(2,i,j)=etpy*detAi
         emom(3,i,j)=etpz*detAi

      end do
      !$omp end parallel do
      return
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


      do
         10     continue
         ! Read file character for character until first whitespace
         keyword=""
         call bytereader(keyword,rd_len,ifile,i_errb)

         ! converting Capital letters
         call caps2small(keyword)

         ! check for comment markers (currently % and #)
         comment=(scan(trim(keyword),'%')==1).or.(scan(trim(keyword),'#')==1).or.&
            (scan(trim(keyword),'*')==1).or.(scan(trim(keyword),'=')==1.or.&
            (scan(trim(keyword),'!')==1))

         if (comment) then
            read(ifile,*)
         else
            ! Parse keyword
            keyword=trim(keyword)
            select case(keyword)

            case('do_region_optimization')
               read(ifile,*,iostat=i_err) OPT_flag_str
               if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err
               if(OPT_flag_str=='Y') then
                  OPT_flag = .true.
               else
                  OPT_flag = .false.
               end if

            case('region_optimization_rebuild_step')
               read(ifile,*,iostat=i_err) OPT_rebuild_time
               if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err

            case('glob_cos_thr')
               read(ifile,*,iostat=i_err) cos_thr
               if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err

            case('do_printatomclass')
               read(ifile,*,iostat=i_err) OPT_printcores_flag_str
               if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err
               if(OPT_printcores_flag_str=='Y') then
                  print_atom_class_flag = .true.
               else
                  print_atom_class_flag = .false.
               end if

            case('do_adaptive_timestepping')
               read(ifile,*,iostat=i_err) adapt_flag_str
               if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err
               if(adapt_flag_str=='Y') then
                  adaptive_time_flag = .true.
               else
                  adaptive_time_flag = .false.
               end if

            case('do_ip_adaptive_timestepping')
               read(ifile,*,iostat=i_err) ip_adapt_flag_str
               if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err
               if(ip_adapt_flag_str=='Y') then
                  ip_adaptive_time_flag = .true.
               else
                  ip_adaptive_time_flag = .false.
               end if

            case('adaptive_timestep_interval')
               read(ifile,*,iostat=i_err) adapt_time_interval
               if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err

            case('ip_adaptive_timestep_interval')
               read(ifile,*,iostat=i_err) ip_adapt_time_interval
               if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err

            case('adapt_to_sc_ratio')
               read(ifile,*,iostat=i_err) adapt_to_sc_ratio
               if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err

            case('larmor_numrev')
               read(ifile,*,iostat=i_err) larmor_numrev
               if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err

            case('larmor_thr')
               read(ifile,*,iostat=i_err) larmor_thr
               if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err

            case('totalsimtime')
               totalsimtime_flag = .true.
               read(ifile,*,iostat=i_err) totalsimtime
               if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err

            end select
         end if

         ! End of file
         if (i_errb==20) goto 20
         ! End of row
         if (i_errb==10) goto 10
      end do

      20  continue

      rewind(ifile)
      return
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
      cos_thr                = 0.99_dblprec
      adaptive_time_flag     = .false.
      ip_adaptive_time_flag  = .false.
      adapt_time_interval    = 100
      ip_adapt_time_interval = 100
      adapt_to_sc_ratio      = 1
      larmor_numrev          = 20
      larmor_thr             = 0.15_dblprec
      totalsimtime           = delta_t
      !
   end subroutine set_opt_defaults

end module optimizationRoutines
