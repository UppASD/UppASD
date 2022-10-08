!-------------------------------------------------------------------------------
! MODULE: temperature
!> @brief module to implement different temperature profiles in the sample
!> @author Jonathan Chico, Lars Bergqvist, Anders Bergman
!> @copyright
!> GNU Public License.
!> @todo generalize to all kind of system shapes
!-------------------------------------------------------------------------------
module temperature
   use Parameters
   use Profiling
   use Constants
   use temperature_common

   implicit none

   real(dblprec), dimension(:), allocatable   :: temp_array      !< temperature as an array
   real(dblprec), dimension(:,:),allocatable  :: iptemp_array    !< temperature (array) for initial phase

   real(dblprec), dimension(:), allocatable   :: temp_s_array      !< spin temperature as an array
   real(dblprec), dimension(:,:),allocatable  :: iptemp_s_array    !< spin temperature (array) for initial phase
   real(dblprec), dimension(:), allocatable   :: temp_l_array      !< lattice temperature as an array
   real(dblprec), dimension(:,:),allocatable  :: iptemp_l_array    !< lattice temperature (array) for initial phase
   real(dblprec), dimension(:), allocatable   :: temp_e_array      !< electron temperature as an array
   real(dblprec), dimension(:,:),allocatable  :: iptemp_e_array    !< electron temperature (array) for initial phase

   real(dblprec),dimension(:),allocatable     :: spintemp_sample  !< temperature (array) for initial phase
   integer :: spintemp_samp_done

   ! temperature gradient
   character(len=1) :: grad                 !< if temperature gradient is active (y/n)
   character(len=20) :: tempfile            !< temperature file
   real(dblprec) :: temp_max                !< amplitude of the gaussian temperature profile
   integer :: temp_solver                   !< solver used for the temperature routine
   integer :: crys_symm                     !< symmetry in the crystallographic cell
   integer :: dim_sys                       !< system simensionality
   integer :: no_shells_num                 !< number of shells considered
   integer :: num_tot_neigh                 !< number of total neighbours
   integer :: temp_max_no_shells            !< actual maximum number of shells
   character(len=1) :: grid_type            !< type of gird being considered
   character(len=1) :: eq_type              !< type of equation being solved (poisson or laplace)
   character(len=1) :: source_type          !< source type for the poisson equation
   character(len=35) :: loadtemp            !< file name for temperature "restartfile"
   character(len=10) :: I1_max_border       !< type of boundary conditions in the i1-max border
   character(len=10) :: I1_min_border       !< type of boundary conditions in the i1-min border
   character(len=10) :: I2_max_border       !< type of boundary conditions in the i2-max border
   character(len=10) :: I2_min_border       !< type of boundary conditions in the i2-min border
   character(len=10) :: I3_max_border       !< type of boundary conditions in the i3-max border
   character(len=10) :: I3_min_border       !< type of boundary conditions in the i3-min border
   real(dblprec) :: temp_I1_min             !< temperature of border i1=0 for constant temp
   real(dblprec) :: temp_I1_max             !< temperature of border i1=n1 for constant temp
   real(dblprec) :: temp_I2_min             !< temperature of border i2=0 for constant temp
   real(dblprec) :: temp_I2_max             !< temperature of border i2=n2 for constant temp
   real(dblprec) :: temp_I3_min             !< temperature of border i3=0 for constant temp
   real(dblprec) :: temp_I3_max             !< temperature of border i3=n3 for constant temp
   real(dblprec) :: temp_I1_min_low         !< minimum tempearture of border i1=0 for linear profile
   real(dblprec) :: temp_I1_min_high        !< maximum tempearture of border i1=0 for linear profile
   real(dblprec) :: temp_I1_max_low         !< minimum tempearture of border i1=n1 for linear profile
   real(dblprec) :: temp_I1_max_high        !< maximum tempearture of border i1=n1 for linear profile
   real(dblprec) :: temp_I2_min_low         !< minimum tempearture of border i2=0 for linear profile
   real(dblprec) :: temp_I2_min_high        !< maximum tempearture of border i2=0 for linear profile
   real(dblprec) :: temp_I2_max_low         !< minimum tempearture of border i2=n2 for linear profile
   real(dblprec) :: temp_I2_max_high        !< maximum tempearture of border i2=n2 for linear profile
   real(dblprec) :: temp_I3_min_low         !< minimum tempearture of border i3=0 for linear profile
   real(dblprec) :: temp_I3_min_high        !< maximum tempearture of border i3=0 for linear profile
   real(dblprec) :: temp_I3_max_low         !< minimum tempearture of border i3=n3 for linear profile
   real(dblprec) :: temp_I3_max_high        !< maximum tempearture of border i3=n3 for linear profile
   real(dblprec), dimension(3) :: r_center  !< center of the gaussian temperature profile
   real(dblprec), dimension(3) :: sigmatemp !< sigma parameter for the gaussian temperature profile
   integer :: init_temp
   integer, dimension(:), allocatable :: temp_nn
   real(dblprec), dimension(:,:,:), allocatable :: temp_neigh_dist

   ! legacy variables
   ! temperature gradient
   character(len=20), dimension(6) :: bounds !< array for the boundary conditions of the temperature profile
   real(dblprec) :: temp_low_x,temp_high_x, temp_high_y, temp_low_y, temp_high_z,temp_low_z !< boundary conditions for the 1d laplace equation
   integer :: barr_size ! size of the stpe for the tprof=1dstep

   ! Three-temperature model

   character(len=1) :: do_3tm               !< Use the 3TM for spin-lattice simulations (Y)/N/E 
   real(dblprec) :: temp_spin_init          !< Initial spin temperature in 3TM (default Temp)
   real(dblprec) :: temp_latt_init          !< Initial spin temperature in 3TM (default Temp)
   real(dblprec) :: temp_elec_init          !< Initial spin temperature in 3TM (default Temp)


   private

   public :: temp_array, iptemp_array, spintemperature, do_3tm
   public :: grad, tempfile,temp_low_x,temp_high_x, temp_high_y, temp_low_y
   public :: allocate_temp, setup_temp, read_temperature, deallocate_temp, set_temperature_defaults
   public :: read_temperature_legacy, lparray
   public :: f_spintemp

contains

   !-----------------------------------------------------------------------------
   ! SUBROUTINE: spintemperature
   !> @brief Calculation of the spin temperature
   !> @author Lars Bergqvist
   !-----------------------------------------------------------------------------
   subroutine spintemperature(natom, mensemble, mstep,ntemp,simid,emomm, beff,flag)

      implicit none

      integer, intent(in) :: natom !< number of atoms in system
      integer, intent(in) :: mensemble !< number of ensembles
      integer, intent(in) :: mstep !< present step
      integer, intent(in) :: ntemp !< total steps
      character(len=8), intent(in) :: simid !< name of simulation
      real(dblprec), dimension(3,natom,mensemble), intent(in) :: emomm  !< current magnetic moment vector
      real(dblprec), dimension(3,natom,mensemble), intent(in) :: beff !< effective field
      integer,intent(in) :: flag !< flag for init,meas and print

      integer :: l,r,i_stat, i_all
      character(len=30) :: filn
      real(dblprec) :: spintemp, spintempnom, spintempdenom

      if (flag==0) then
         allocate(spintemp_sample(ntemp),stat=i_stat)
         call memocc(i_stat,product(shape(spintemp_sample))*kind(spintemp_sample),'spintemp_sample','spintemperature')
         spintemp_sample=0.0_dblprec
         spintemp_samp_done=1
      endif
      if (flag==1) then
         spintemp=0.0_dblprec
         spintempnom=0_dblprec
         spintempdenom=0_dblprec
         do l=1,mensemble
            do r=1,natom
               spintempnom = spintempnom + &
                  ((emomm(2,r,l)*beff(3,r,l)-emomm(3,r,l)*beff(2,r,l))**2 + &
                  (emomm(3,r,l)*beff(1,r,l)-emomm(1,r,l)*beff(3,r,l))**2 + &
                  (emomm(1,r,l)*beff(2,r,l)-emomm(2,r,l)*beff(1,r,l))**2)
               spintempdenom = spintempdenom + &
                  (emomm(1,r,l)*beff(1,r,l)+emomm(2,r,l)*beff(2,r,l)+ emomm(3,r,l)*beff(3,r,l))
            end do
         enddo
         spintemp = spintempnom / spintempdenom

         spintemp=spintemp/mensemble
         spintemp=(mub*spintemp)/(2*k_bolt)

         spintemp_sample(spintemp_samp_done)=spintemp

         write (filn,'(''spintemp.'',a,''.out'')') trim(simid)
         open(ofileno,file=filn, position='append')

         write(ofileno,20002) mstep,spintemp_sample(spintemp_samp_done), &
            sum(spintemp_sample(1:spintemp_samp_done))/spintemp_samp_done, &
            sum(abs(spintemp_sample(1:spintemp_samp_done)))/spintemp_samp_done

         close(ofileno)
         spintemp_samp_done=spintemp_samp_done+1

      endif
      if (flag==2) then
         i_all=-product(shape(spintemp_sample))*kind(spintemp_sample)
         deallocate(spintemp_sample,stat=i_stat)
         call memocc(i_stat,i_all,'spintemp_sample','spintemperature')
      endif

      20002 format (i8,2x,es16.8,es16.8,es16.8)

   end subroutine spintemperature

   !-----------------------------------------------------------------------------
   ! FUNCTION: f_spintemp
   !> @brief Calculation of the spin temperature 
   !> @author Anders Bergman (based on spintemperature())
   !-----------------------------------------------------------------------------
   function f_spintemp(natom, mensemble, emomm, beff) result(spintemp)

      implicit none

      integer, intent(in) :: natom !< number of atoms in system
      integer, intent(in) :: mensemble !< number of ensembles
      real(dblprec), dimension(3,natom,mensemble), intent(in) :: emomm  !< current magnetic moment vector
      real(dblprec), dimension(3,natom,mensemble), intent(in) :: beff !< effective field

      integer :: l,r
      real(dblprec) :: spintemp
      real(dblprec) ::  spintempnom, spintempdenom

         spintemp=0.0_dblprec
         spintempnom=0_dblprec
         spintempdenom=0_dblprec
         do l=1,mensemble
            do r=1,natom
               spintempnom = spintempnom + &
                  ((emomm(2,r,l)*beff(3,r,l)-emomm(3,r,l)*beff(2,r,l))**2 + &
                  (emomm(3,r,l)*beff(1,r,l)-emomm(1,r,l)*beff(3,r,l))**2 + &
                  (emomm(1,r,l)*beff(2,r,l)-emomm(2,r,l)*beff(1,r,l))**2)
               spintempdenom = spintempdenom + &
                  (emomm(1,r,l)*beff(1,r,l)+emomm(2,r,l)*beff(2,r,l)+ emomm(3,r,l)*beff(3,r,l))
            end do
         enddo
         spintempdenom = spintempdenom + 1.0e-12_dblprec
         spintemp = spintempnom / spintempdenom

         spintemp=spintemp/mensemble
         spintemp=(mub*spintemp)/(2*k_bolt)

   end function f_spintemp

   !-----------------------------------------------------------------------------
   ! SUBROUTINE: allocate_temp
   !> @brief Subroutine to allocate the arrays that contain the temperature
   !> @author Jonathan Chico
   !-----------------------------------------------------------------------------
   subroutine allocate_temp(natom, ip_nphase)

      implicit none

      integer, intent(in) :: natom
      integer, intent(in) :: ip_nphase

      integer :: i_stat

      allocate(temp_array(natom),stat=i_stat)
      call memocc(i_stat,product(shape(temp_array))*kind(temp_array),'temp_array','allocate_temp')
      if(ip_nphase>0) then
         allocate(iptemp_array(natom,ip_nphase),stat=i_stat)
         call memocc(i_stat,product(shape(iptemp_array))*kind(iptemp_array),'iptemp_array','allocate_temp')
      end if

   end subroutine allocate_temp

   !-----------------------------------------------------------------------------
   ! SUBROUTINE: deallocate_temp
   !> @brief Subroutine to deallocate the arrays that contain the temperature and
   !> all the temperature gradients arrays
   !> @author Jonathan Chico
   !-----------------------------------------------------------------------------
   subroutine deallocate_temp()

      integer :: i_stat,i_all

      i_all=-product(shape(temp_array))*kind(temp_array)
      deallocate(temp_array,stat=i_stat)
      call memocc(i_stat,i_all,'temp_array','deallocate_temp')
      if(allocated(iptemp_array)) then
         i_all=-product(shape(iptemp_array))*kind(iptemp_array)
         deallocate(iptemp_array,stat=i_stat)
         call memocc(i_stat,i_all,'iptemp_array','deallocate_temp')
      end if

      if (allocated(temp_nlist)) then
         i_all=-product(shape(temp_nlist))*kind(temp_nlist)
         deallocate(temp_nlist,stat=i_stat)
         call memocc(i_stat,i_all,'temp_nlist','deallocate_temp')
      endif

      if (allocated(temp_nm)) then
         i_all=-product(shape(temp_nm))*kind(temp_nm)
         deallocate(temp_nm,stat=i_stat)
         call memocc(i_stat,i_all,'temp_nm','deallocate_temp')
      endif

      if (allocated(temp_nmdim)) then
         i_all=-product(shape(temp_nmdim))*kind(temp_nmdim)
         deallocate(temp_nmdim,stat=i_stat)
         call memocc(i_stat,i_all,'temp_nmdim','deallocate_temp')
      endif

      if (allocated(temp_nlistsize)) then
         i_all=-product(shape(temp_nlistsize))*kind(temp_nlistsize)
         deallocate(temp_nlistsize,stat=i_stat)
         call memocc(i_stat,i_all,'temp_nlistsize','deallocate_temp')
      endif

   end subroutine deallocate_temp

   !-----------------------------------------------------------------------------
   ! SUBROUTINE setup_temp
   !> @brief Setup the temperature array of the system
   !> @author Jonathan Chico
   !-----------------------------------------------------------------------------
   subroutine setup_temp(Natom,NT,NA,N1,N2,N3,Natom_full,do_ralloy,atype,acellnumb,&
         atype_ch,simid,temp,C1,C2,C3,BC1,BC2,BC3,Bas,coord,temp_array)

      implicit none

      integer, intent(in) :: NT  !< number of types of atoms
      integer, intent(in) :: NA  !< number of atoms in one cell
      integer, intent(in) :: N1  !< number of cell repetitions in x direction
      integer, intent(in) :: N2  !< number of cell repetitions in y direction
      integer, intent(in) :: N3  !< number of cell repetitions in z direction
      integer, intent(in) :: Natom      !< number of atoms in system
      integer, intent(in) :: Natom_full !< number of atoms for full system (=natom if not dilute)
      integer, intent(in), optional :: do_ralloy  !< random alloy simulation (0/1)
      integer, dimension(Natom), intent(in) :: atype !< type of atom
      integer, dimension(:), optional ,intent(in) :: acellnumb !< list for translating atom no. in full cell to actual cell
      integer, dimension(:), optional ,intent(in) :: atype_ch  !< actual type of atom for dilute system
      character(len=1), intent(in) :: BC1  !< boundary conditions in x-direction
      character(len=1), intent(in) :: BC2  !< boundary conditions in y-direction
      character(len=1), intent(in) :: BC3  !< boundary conditions in z-direction
      character(len=8), intent(in) :: simid       !< simulation name
      real(dblprec), intent(in) :: temp              !< constant temperature from input
      real(dblprec), dimension(3), intent(in) :: C1       !< first lattice vector
      real(dblprec), dimension(3), intent(in) :: C2       !< second lattice vector
      real(dblprec), dimension(3), intent(in) :: C3       !< third lattice vector
      real(dblprec), dimension(3,NA), intent(in) :: Bas   !< coordinates for basis atoms
      real(dblprec), dimension(3,Natom), intent(in) :: coord !< coordinates of the atoms in the system
      real(dblprec), dimension(Natom), intent(out) :: temp_array !< temperature array

      !.. local variables
      integer :: temp_max_no_neigh !< calculated maximum of neighbours for temperature grid
      integer :: count_I1_min  !< number of atoms belonging to the I1=0 border
      integer :: count_I1_max  !< number of atoms belonging to the I1=n1 border
      integer :: count_I2_min  !< number of atoms belonging to the I2=0 border
      integer :: count_I2_max  !< number of atoms belonging to the I1=n2 border
      integer :: count_I3_min  !< number of atoms belonging to the I3=0 border
      integer :: count_I3_max  !< number of atoms belonging to the I3=n3 border
      integer, dimension(:), allocatable :: temp_nlistsize !< size of neighbour list for temperature grid
      integer, dimension(:), allocatable :: borders_I1_min !< list of atoms belonging to the I1=0 border
      integer, dimension(:), allocatable :: borders_I1_max !< list of atoms belonging to the I1=n1 border
      integer, dimension(:), allocatable :: borders_I2_min !< list of atoms belonging to the I2=0 border
      integer, dimension(:), allocatable :: borders_I2_max !< list of atoms belonging to the I2=n2 border
      integer, dimension(:), allocatable :: borders_I3_min !< list of atoms belonging to the I3=0 border
      integer, dimension(:), allocatable :: borders_I3_max !< list of atoms belonging to the I3=n3 border
      integer, dimension(:,:), allocatable :: temp_nlist !< neighbour list for temperature grid
      logical, dimension(Natom) :: border_belong !< array to identify wether an atom belongs or not in a border

      real(dblprec), dimension(Natom) :: temp_init !< initial temperature + source terms

      integer :: i, i_stat,i_all
      real(dblprec) :: tol

      real(dblprec) :: c_fac, rmax, err
      character(len=30) :: filn
      integer :: weight_type, dim_base_poly, nmax,itmax, iter, itol

      ! if there is no gradient the array is filled with a constant temperature
      if (grad.eq.'N') then
         temp_array=temp
         write(*,'(2x,a)') ' Homogeneous Temperature'
         ! if the gradient is set to f one reads the temperature from a file
      else if (grad.eq.'F') then
         write(*,'(2x,a)') ' Read Temperature from file'
         temp_init=0.0_dblprec
         !call readloadtemp(natom,loadtemp,temp_array)
         call readloadtemp(natom,tempfile,temp_array)
         ! if the gradient is set to y the gradient is calculated
      else if (grad.eq.'Y') then
         write(*,'(2x,a)') ' Calculate Gradient'
         ! setup the temperature grid to find border atoms
         call setup_temp_grid(NT,NA,N1,N2,N3,Natom,crys_symm,Natom_full,temp_max_no_shells,&
            do_ralloy,temp_nn,atype,atype_ch,acellnumb,C1,C2,C3,BC1,BC2,BC3,Bas,&
            temp_neigh_dist,temp_max_no_neigh,temp_max_no_equiv,temp_nm,temp_nmdim,temp_nlistsize,&
            temp_nlist,count_I1_min,count_I1_max,count_I2_min,count_I2_max,count_I3_min,count_I3_max,&
            borders_i1_min,borders_i1_max,borders_i2_min,borders_i2_max,borders_i3_min,borders_i3_max,&
            border_belong,simid)

         ! setup the initial temperature state, boundaries and source terms
         call temperature_init(Natom,init_temp,count_I1_min,count_I1_max,count_I2_min,count_I2_max,&
            count_I3_min,count_I3_max,borders_I1_min,borders_I1_max,borders_I2_min,borders_I2_max,&
            borders_I3_min,borders_I3_max,eq_type,source_type,loadtemp,I1_min_border,I1_max_border,&
            I2_min_border,I2_max_border,I3_min_border,I3_max_border,temp_max,temp_I1_min,&
            temp_I1_max,temp_I2_min,temp_I2_max,temp_I3_min,temp_I3_max,temp_I1_min_low,temp_I1_min_high,&
            temp_I1_max_low,temp_I1_max_high,temp_I2_min_low,temp_I2_min_high,temp_I2_max_low,&
            temp_I2_max_high,temp_I3_min_low,temp_I3_min_high,temp_I3_max_low,temp_I3_max_high,&
            r_center,sigmatemp,coord,temp_init)

         ! third the algorithm used must be selected (meshfree, gauss-seidel, sor methods) this is important for one can solve the poisson or laplace equation
         ! meshfree solver
         if (temp_solver.eq.2) then
            ! once the intial temperature (really more boundary conditions and source terms) one calls the mls function creation
            c_fac=0.0_dblprec
            weight_type=2
            rmax=15.0_dblprec
            dim_base_poly=10
            ! call the mls shape function creator
            call mean_least_squares_method(natom,weight_type,dim_base_poly,c_fac,&
               rmax,temp_max_no_neigh,temp_nlistsize,temp_nlist,coord,border_belong)

            ! store the stress matrix in a sparse matrix way (must allocate the a_matrix and see if nmax can be different from n_size)
            call dsprsin(s_matrix,Natom,nmax,a_matrix,ija)

            deallocate(s_matrix,stat=i_stat)
            call memocc(i_stat,-product(shape(s_matrix))*kind(s_matrix),'s_matrix','setup_temp')

            ! this solves a linear equation where the matrix is solved in a sparse way
            itmax=1000
            itol=1
            call linbcg(Natom,nmax,a_matrix,temp_init,temp_array,itol,tol,itmax,iter,err,ija)

         else if (temp_solver.eq.1) then
            ! these are the finite difference methods solvers
            call finite_differences(Natom,dim_sys,eq_type,temp_solver,temp_max_no_neigh,&
               temp_nlistsize,temp_nlist,grid_type,temp_init,coord,border_belong,simid,temp_array)
         endif

      endif

      ! printing the final temperature file and for comparison the intial temperature file
      if (grad.eq.'Y'.or.grad.eq.'F') then
         write(filn,'(''temperature.'',a,''.out'')') trim(simid)
         open(ofileno,file=filn,position="append")
         do i=1, Natom
            write(ofileno,'(i8,3f16.8,f16.8,2x,f16.8)') i, coord(1:3,i),temp_array(i), temp_init(i)
         enddo
         close(ofileno)
      endif

      if (allocated(borders_I1_min)) then
         i_all=-product(shape(borders_I1_min))*kind(borders_I1_min)
         deallocate(borders_I1_min,stat=i_stat)
         call memocc(i_stat,i_all,'borders_i1_min','setup_temp')
      endif
      if (allocated(borders_i2_min)) then
         i_all=-product(shape(borders_I2_min))*kind(borders_I2_min)
         deallocate(borders_I2_min,stat=i_stat)
         call memocc(i_stat,i_all,'borders_I2_min','setup_temp')
      endif
      if (allocated(borders_I3_min)) then
         i_all=-product(shape(borders_I3_min))*kind(borders_I3_min)
         deallocate(borders_I3_min,stat=i_stat)
         call memocc(i_stat,i_all,'borders_I3_min','setup_temp')
      endif
      if (allocated(borders_I1_max)) then
         i_all=-product(shape(borders_I1_max))*kind(borders_I1_max)
         deallocate(borders_I1_max,stat=i_stat)
         call memocc(i_stat,i_all,'borders_I1_max','setup_temp')
      endif
      if (allocated(borders_I2_max)) then
         i_all=-product(shape(borders_I2_max))*kind(borders_I2_max)
         deallocate(borders_I2_max,stat=i_stat)
         call memocc(i_stat,i_all,'borders_I2_max','setup_temp')
      endif
      if (allocated(borders_I3_max)) then
         i_all=-product(shape(borders_I3_max))*kind(borders_I3_max)
         deallocate(borders_I3_max,stat=i_stat)
         call memocc(i_stat,i_all,'borders_I3_max','setup_temp')
      endif

      if (allocated(temp_nlist)) then
         i_all=-product(shape(temp_nlist))*kind(temp_nlist)
         deallocate(temp_nlist,stat=i_stat)
         call memocc(i_stat,i_all,'temp_nlist','setup_temp')
      endif
      if (allocated(temp_nlistsize)) then
         i_all=-product(shape(temp_nlistsize))*kind(temp_nlistsize)
         deallocate(temp_nlistsize,stat=i_stat)
         call memocc(i_stat,i_all,'temp_nlistsize','setup_temp')
      endif

      return

   end subroutine setup_temp

   !-----------------------------------------------------------------------------
   ! SUBROUTINE: read_temperature
   !> @brief Routine to read the temperature file for thermal gradients
   !> @author Anders Bergman, Jonathan Chico
   !-----------------------------------------------------------------------------
   subroutine read_temperature()

      use FileParser
      use InputData, only: maptype, c1, c2, c3, bas, nt, posfiletype, atype_inp

      implicit none

      integer :: i_err,rd_len,i_errb
      integer :: num_tot_neigh,isite,i_stat,jsite
      integer :: itype,jtype,iline,ishell,i_all
      logical :: unique
      logical :: comment
      real(dblprec):: tol, norm,t1,t2
      real(dblprec), dimension(3) :: r_tmp, r_red
      integer :: no_shells_num
      character(len=50) :: keyword
      real(dblprec), dimension(:,:,:), allocatable :: temp_neigh_dist_tmp

      ! set tolerance for neighbour shells
      tol=1.0d-5

      open(ifileno,file=trim(tempfile))
      do
         10 continue

         ! read file character for character until first whitespace
         keyword=""
         call bytereader(keyword,rd_len,ifileno,i_errb)

         ! converting capital letters
         call caps2small(keyword)

         ! check for comment markers (currently % and #)
         comment=(scan(trim(keyword),'%')==1).or.(scan(trim(keyword),'#')==1).or.&
            (scan(trim(keyword),'*')==1).or.(scan(trim(keyword),'=')==1)

         ! parse keyword
         keyword=trim(keyword)

         select case(keyword)
         case('shells_nums')
            read(ifileno,*,iostat=i_err) no_shells_num
            if(i_err/=0) write(*,*) 'error: reading ',trim(keyword),' data',i_err
         case('num_tot_neigh')
            read(ifileno,*,iostat=i_err) num_tot_neigh
            if(i_err/=0) write(*,*) 'error: reading ',trim(keyword),' data',i_err
         case('temp_solver')
            read(ifileno,*,iostat=i_err) temp_solver
            if(i_err/=0) write(*,*) 'error: reading ',trim(keyword),' data',i_err
         case('crys_symm')
            read(ifileno,*,iostat=i_err) crys_symm
            if(i_err/=0) write(*,*) 'error: reading ',trim(keyword),' data',i_err
         case('grid_type')
            read(ifileno,*,iostat=i_err) grid_type
            if(i_err/=0) write(*,*) 'error: reading ',trim(keyword),' data',i_err
         case('eq_type')
            read(ifileno,*,iostat=i_err) eq_type
            if(i_err/=0) write(*,*) 'error: reading ',trim(keyword),' data',i_err
         case('init_temp')
            read(ifileno,*,iostat=i_err) init_temp
            if(i_err/=0) write(*,*) 'error: reading ',trim(keyword),' data',i_err
         case('dim_sys')
            read(ifileno,*,iostat=i_err) dim_sys
            if(i_err/=0) write(*,*) 'error: reading ',trim(keyword),' data',i_err
         case('x_min_border')
            read(ifileno,*,iostat=i_err) I1_min_border
            if(i_err/=0) write(*,*) 'error: reading ',trim(keyword),' data',i_err
            if (I1_min_border.eq.'constant') then
               read(ifileno,*,iostat=i_err) temp_I1_min
            else if(I1_min_border.eq.'linear') then
               read(ifileno,*,iostat=i_err) t1,t2
               temp_I1_min_low=min(t1,t2)
               temp_I1_min_high=max(t1,t2)
            endif
         case('x_max_border')
            read(ifileno,*,iostat=i_err) i1_max_border
            if(i_err/=0) write(*,*) 'error: reading ',trim(keyword),' data',i_err
            if (I1_max_border.eq.'constant') then
               read(ifileno,*,iostat=i_err) temp_I1_max
            else if(I1_max_border.eq.'linear') then
               read(ifileno,*,iostat=i_err) t1,t2
               temp_I1_max_low=min(t1,t2)
               temp_I1_max_high=max(t1,t2)
            endif
         case('y_min_border')
            read(ifileno,*,iostat=i_err) I2_min_border
            if(i_err/=0) write(*,*) 'error: reading ',trim(keyword),' data',i_err
            if (I2_min_border.eq.'constant') then
               read(ifileno,*,iostat=i_err) temp_I2_min
            else if(I2_min_border.eq.'linear') then
               read(ifileno,*,iostat=i_err) t1,t2
               temp_I2_min_low=min(t1,t2)
               temp_I2_min_high=max(t1,t2)
            endif
         case('y_max_border')
            read(ifileno,*,iostat=i_err) I2_max_border
            if(i_err/=0) write(*,*) 'error: reading ',trim(keyword),' data',i_err
            if (I2_max_border.eq.'constant') then
               read(ifileno,*,iostat=i_err) temp_i2_max
            else if(I2_max_border.eq.'linear') then
               read(ifileno,*,iostat=i_err) t1,t2
               temp_I2_max_low=min(t1,t2)
               temp_I2_max_high=max(t1,t2)
            endif
         case('z_min_border')
            read(ifileno,*,iostat=i_err) i3_min_border
            if(i_err/=0) write(*,*) 'error: reading ',trim(keyword),' data',i_err
            if (I3_min_border.eq.'constant') then
               read(ifileno,*,iostat=i_err) temp_I3_min
            else if(I3_min_border.eq.'linear') then
               read(ifileno,*,iostat=i_err) t1,t2
               temp_I3_min_low=min(t1,t2)
               temp_I3_min_high=max(t1,t2)
            endif
         case('z_max_border')
            read(ifileno,*,iostat=i_err) i3_max_border
            if(i_err/=0) write(*,*) 'error: reading ',trim(keyword),' data',i_err
            if (I3_max_border.eq.'constant') then
               read(ifileno,*,iostat=i_err) temp_I3_max
            else if(I3_max_border.eq.'linear') then
               read(ifileno,*,iostat=i_err) t1,t2
               temp_I3_max_low=min(t1,t2)
               temp_I3_max_high=max(t1,t2)
            endif
         case('source_type')
            read(ifileno,*,iostat=i_err) source_type
            if(i_err/=0) write(*,*) 'error: reading ',trim(keyword),' data',i_err
            if (source_type.eq.'G') then
               read(ifileno,*,iostat=i_err) temp_max
               read(ifileno,*,iostat=i_err) sigmatemp(1), sigmatemp(2), sigmatemp(3)
               read(ifileno,*,iostat=i_err) r_center(1), r_center(2), r_center(3)
            endif
         case('load_temp')
            read(ifileno,*,iostat=i_err) loadtemp
            if(i_err/=0) write(*,*) 'error: reading ',trim(keyword),' data',i_err
         case('temp_neigh')
            allocate(temp_neigh_dist_tmp(nt,no_shells_num,3),stat=i_stat)
            call memocc(i_stat,product(shape(temp_neigh_dist_tmp))*kind(temp_neigh_dist_tmp),'temp_neigh_dist_tmp','read_temperature')
            allocate(temp_nn(nt),stat=i_stat)
            call memocc(i_stat,product(shape(temp_nn))*kind(temp_nn),'temp_nn','read_temperature')
            temp_nn(1:nt)=0
            do iline=1, num_tot_neigh
               ! read indices and coordinates
               read (ifileno,*) isite, jsite, r_tmp
               ! find type of site
               itype=atype_inp(isite)
               jtype=atype_inp(jsite)
               if(maptype==2) then
                  ! calculate proper neighbour vector (from "bgfm")
                  r_red(1)=bas(1,jsite)-bas(1,isite)+c1(1)*r_tmp(1)+c2(1)*r_tmp(2)+c3(1)*r_tmp(3)
                  r_red(2)=bas(2,jsite)-bas(2,isite)+c1(2)*r_tmp(1)+c2(2)*r_tmp(2)+c3(2)*r_tmp(3)
                  r_red(3)=bas(3,jsite)-bas(3,isite)+c1(3)*r_tmp(1)+c2(3)*r_tmp(2)+c3(3)*r_tmp(3)
               else
                  ! calculates neighbour vectors from direct coordinates or cartesian
                  ! coordinates, corresponding to how the atomic positions are entered
                  if (posfiletype=='C') then
                     r_red=r_tmp
                  elseif (posfiletype=='D') then
                     r_red(1)=r_tmp(1)*c1(1)+r_tmp(2)*c2(1)+r_tmp(3)*c3(1)
                     r_red(2)=r_tmp(1)*c1(2)+r_tmp(2)*c2(2)+r_tmp(3)*c3(2)
                     r_red(3)=r_tmp(1)*c1(3)+r_tmp(2)*c2(3)+r_tmp(3)*c3(3)
                  else
                     stop 'only posfiletype= c or d is currently supported'
                  endif
               end if
               ! loop through earlier vectors to find equivalent shells
               unique=.true.
               do ishell=1,temp_nn(itype)
                  norm=(r_red(1)-temp_neigh_dist_tmp(itype,ishell,1))**2+ &
                     (r_red(2)-temp_neigh_dist_tmp(itype,ishell,2))**2+ &
                     (r_red(3)-temp_neigh_dist_tmp(itype,ishell,3))**2
                  if(norm<tol) then
                     unique=.false.
                  end if
               end do
               if (unique) then
                  temp_nn(itype)=temp_nn(itype)+1
                  temp_neigh_dist_tmp(itype,temp_nn(itype),1:3)=r_red(1:3)
               endif
            enddo

            temp_max_no_shells=maxval(temp_nn)

            allocate(temp_neigh_dist(nt,temp_max_no_shells,3),stat=i_stat)
            call memocc(i_stat,product(shape(temp_neigh_dist))*kind(temp_neigh_dist),'temp_neigh_dist','read_temperature')
            temp_neigh_dist=0
            !
            do ishell=1,temp_max_no_shells
               do itype=1,nt
                  temp_neigh_dist(itype,ishell,:)=temp_neigh_dist_tmp(itype,ishell,:)
               end do
            end do
            !
            i_all=-product(shape(temp_neigh_dist_tmp))*kind(temp_neigh_dist_tmp)
            deallocate(temp_neigh_dist_tmp,stat=i_stat)
            call memocc(i_stat,i_all,'temp_neigh_dist_tmp','read_temperature')

         case default
            if(comment.or.len(trim(keyword))==0) then
            else
               print *,"keyword '",trim(keyword),"' is not recognized"
            end if

         end select

         ! end of file
         if (i_errb==20) goto 20
         ! end of row
         if (i_errb==10) goto 10

      end do
      goto 30

      20  continue
      30  continue

      return

      close(ofileno)

   end subroutine read_temperature

   !-----------------------------------------------------------------------------
   ! SUBROUTINE: read_temperature_legacy
   !> @brief Routine to read the needed info for the legacy version of the temperature
   !> gradient
   !-----------------------------------------------------------------------------
   subroutine read_temperature_legacy()

      use FileParser

      implicit none

      integer :: i_err, rd_len, i_errb
      character(len=50) :: keyword
      logical :: comment

      open(500,file=trim(tempfile))
      do
         10   continue
         ! read file character for character until first whitespace
         keyword=""
         call bytereader(keyword,rd_len,500,i_errb)
         ! converting capital letters
         call caps2small(keyword)
         ! check for comment markers (currently % and #)
         comment=(scan(trim(keyword),'%')==1).or.(scan(trim(keyword),'#')==1).or.&
            (scan(trim(keyword),'*')==1).or.(scan(trim(keyword),'=')==1)
         ! parse keyword
         keyword=trim(keyword)
         select case(keyword)

         case('tbound_x')

            read(500,*,iostat=i_err) bounds(1),bounds(2)
            read(500,*,iostat=i_err) temp_high_x, temp_low_x
            if((bounds(1).eq.'step').or.(bounds(2).eq.'step')) then
               read(500,*,iostat=i_err) barr_size
            else if ((bounds(1).eq.'step').and.(bounds(2).eq.'step')) then
               read(500,*,iostat=i_err) barr_size
            end if
            if(i_err/=0) write(*,*) 'error: reading ',trim(keyword),' data',i_err

         case('tbound_y')

            read(500,*,iostat=i_err) bounds(3),bounds(4)
            if((bounds(3).ne.'n').and.(bounds(4).ne.'n')) then
               read(500,*,iostat=i_err) temp_high_y, temp_low_y
            end if
            if(i_err/=0) write(*,*) 'error: reading ',trim(keyword),' data',i_err

         case('tbound_z')

            read(500,*,iostat=i_err) bounds(5),bounds(6)
            if((bounds(5).ne.'n').and.(bounds(6).ne.'n')) then
               read(500,*,iostat=i_err) temp_high_z, temp_low_z
            end if
            if(i_err/=0) write(*,*) 'error: reading ',trim(keyword),' data',i_err

         case default
            if(comment.or.len(trim(keyword))==0) then
            else
               print *,"keyword '",trim(keyword),"' is not recognized"
            end if

         end select
         ! end of file
         if (i_errb==20) goto 20
         ! end of row
         if (i_errb==10) goto 10
      end do
      goto 30
      20  continue
      30  continue
      return

      close(500)

   end subroutine read_temperature_legacy

   !-----------------------------------------------------------------------------
   ! SUBROUTINE: set_temperature_defaults
   !> @brief Sets the default values for the temperature gradient
   !> @author Jonathan Chico
   !-----------------------------------------------------------------------------
   subroutine set_temperature_defaults()
      !
      implicit none
      !
      !temperature gradient
      tempfile         = 'tempfile'
      grad             = 'N'
      no_shells_num    = 0
      num_tot_neigh    = 0
      temp_solver      = 1
      crys_symm        = 1
      grid_type        = 'N'
      eq_type          = 'N'
      init_temp        = 1
      dim_sys          = 3
      I1_min_border    = 'N'
      I1_max_border    = 'N'
      I2_min_border    = 'N'
      I2_max_border    = 'N'
      I3_min_border    = 'N'
      I3_max_border    = 'N'
      temp_max         = 0.0_dblprec
      temp_I1_max      = 0.0_dblprec
      temp_I2_max      = 0.0_dblprec
      temp_I3_max      = 0.0_dblprec
      temp_I1_min      = 0.0_dblprec
      temp_I2_min      = 0.0_dblprec
      temp_I2_min      = 0.0_dblprec
      temp_I1_max_low  = 0.0_dblprec
      temp_I2_max_low  = 0.0_dblprec
      temp_I3_max_low  = 0.0_dblprec
      temp_I1_min_low  = 0.0_dblprec
      temp_I2_min_low  = 0.0_dblprec
      temp_I3_min_low  = 0.0_dblprec
      temp_I1_max_high = 0.0_dblprec
      temp_I2_max_high = 0.0_dblprec
      temp_I3_max_high = 0.0_dblprec
      temp_I1_min_high = 0.0_dblprec
      temp_I2_min_high = 0.0_dblprec
      temp_I3_min_high = 0.0_dblprec
      sigmatemp        = 0.0_dblprec
      r_center         = 0.0_dblprec

   end subroutine set_temperature_defaults

   !-----------------------------------------------------------------------------
   ! SUBROUTINE: lparray
   !> @brief subroutine to create a site dependant temperature that follows the
   !> laplace equation for a time independet configuration
   !> @author Jonathan Chico
   !-----------------------------------------------------------------------------
   subroutine lparray(tarray,natom,coord,temp,simid,printme)

      real(dblprec), intent(in) :: temp  !< temperature
      integer, intent(in) :: natom !< array size
      real(dblprec), dimension(3,natom), intent(in) :: coord !< coordinates of the atoms in the crystal
      real(dblprec), dimension(natom), intent(out) :: tarray !< temperature array
      character(len=8), intent(in) :: simid  !< name of simulation
      logical, intent(in) :: printme  !< flag to determine printing of temperature file

      integer :: i,n
      character(len=30) :: filn

      n=1

      !print '(4f10.4)',temp_low_y,temp_high_y,temp_low_x,temp_high_x
      do i=1, natom
         ! fill up the array with the results of the chosen function given by the boundary conditions
         tarray(i) = lpfunction(temp_low_y,temp_high_y,temp_low_x,temp_high_x,coord,natom,i,n,grad,bounds,temp)
      end do

      ! printing the final temperature file and for comparison the intial temperature file
      if (printme) then
         write(filn,'(''temperature.'',a,''.out'')') trim(simid)
         open(45,file=filn)
         do i=1, natom
            write(45,'(i8,3f16.8,f16.8)') i, coord(1:3,i),tarray(i)
         enddo
         close(45)
      endif


   end subroutine lparray

   !-----------------------------------------------------------------------------
   !> @brief this creates a linear function to fill up a 1d lattice
   !> @author Jonathan Chico
   !-----------------------------------------------------------------------------
   function linear_1d(natom,temp_high_x,temp_low_x,coord,i)

      real(dblprec), intent(in) :: temp_high_x, temp_low_x ! boundaries
      real(dblprec), dimension(3,natom), intent(in) :: coord !> coordinates of the  atoms in the lattice
      real(dblprec) :: linear_1d
      integer, intent(in) :: Natom
      integer :: i, x_size

      x_size = maxval(coord(1,:)) ! this works only if the person puts the direction along x maybe put in a logical variable to check

      linear_1d = temp_low_x + (temp_high_x-temp_low_x)*coord(1,i)/x_size


   end function linear_1d

   ! function to include a step function for a 1d array
   function step_x(barr_size,natom,temp_high_x,temp_low_x,coord,i)

      real(dblprec), intent(in) :: temp_high_x, temp_low_x ! temperatures for the step
      real(dblprec), dimension(3,natom), intent(in) :: coord ! coordinates of the atoms in the lattice
      real(dblprec) :: step_x
      integer, intent(in) :: natom, barr_size ! number of atoms and size of the hot part of the lattice
      integer :: i

      if(int(coord(1,i)).le.barr_size) then
         step_x = temp_high_x
      else
         step_x = temp_low_x
      end if

   end function step_x

   ! this deals with the boundary condition f(x,ysize)=tymax
   function cts_2d_x_max(ty_max,n,coord,natom,i)

      integer, intent(in) :: n, i, natom
      real(dblprec), intent(in) :: ty_max ! boundary condition
      real(dblprec), dimension(3,natom),intent(in) :: coord ! coordinates of the atoms in the lattice
      real(dblprec) :: cts_2d_x_max
      real(dblprec) x_size, pi, arg_x, y_size, x , y

      parameter(pi=4._dblprec*datan(1._dblprec)) ! the best way to define pi

      y_size = maxval(coord(2,:)) ! to define the size of the integration region in the y direction
      x_size = maxval(coord(1,:)) ! to define the size of the integration region in the x direction

      x = coord(1,i)
      y = coord(2,i)

      arg_x=n*pi/y_size

      cts_2d_x_max = 2*ty_max*(1-(-1)**n)*sinh(arg_x*x)*sin(arg_x*y)
      cts_2d_x_max = cts_2d_x_max/(n*pi*sinh(arg_x*x_size))

   end function cts_2d_x_max

   ! this deals with the boundary condition f(x,ysize)=tymin
   function cts_2d_x_min(ty_min,n,coord,natom,i)

      integer, intent(in) :: n, i, natom
      real(dblprec), intent(in) :: ty_min ! boundary condition
      real(dblprec), dimension(3,natom),intent(in) :: coord ! coordinates for the atoms in the lattice
      real(dblprec) :: cts_2d_x_min
      real(dblprec) x_size, pi, arg_x, y_size, x, y

      parameter(pi=4._dblprec*datan(1._dblprec)) ! the best way to define pi

      y_size = maxval(coord(2,:)) ! size of the lattice in the y direction
      x_size = maxval(coord(1,:)) ! size of the lattice in the x direction

      x = coord(1,i)
      y = coord(2,i)

      arg_x=n*pi/y_size

      cts_2d_x_min = -2*ty_min*(1-(-1)**n)*sinh(arg_x*(x-x_size))*sin(arg_x*y)
      cts_2d_x_min = cts_2d_x_min/(n*pi*sinh(arg_x*x_size))

   end function cts_2d_x_min

   ! this deals with the boundary condition f(x,ysize)=txmax
   function cts_2d_y_max(tx_max,n,coord,natom,i)

      integer, intent(in) :: n, i, natom
      real(dblprec), intent(in) :: tx_max ! boundary condition
      real(dblprec), dimension(3,natom),intent(in) :: coord ! coordinates of the atoms in the lattice
      real(dblprec) :: cts_2d_y_max
      real(dblprec) x_size, pi, arg_y, y_size, x ,y

      parameter(pi=4._dblprec*datan(1._dblprec)) ! the best way to define pi

      y_size = maxval(coord(2,:)) ! size of the lattice in the y direction
      x_size = maxval(coord(1,:)) ! size of the lattice in the x direction

      x = coord(1,i)
      y = coord(2,i)

      arg_y=n*pi/x_size

      cts_2d_y_max = 2*tx_max*(1-(-1)**n)*sinh(arg_y*y)*sin(arg_y*x)
      cts_2d_y_max = cts_2d_y_max/(n*pi*sinh(arg_y*y_size))

   end function cts_2d_y_max

   !! this deals with the boundary condition f(x,ysize)=txmin
   function cts_2d_y_min(tx_min,n,coord,natom,i)

      integer, intent(in) :: n, i, natom
      real(dblprec), intent(in) :: tx_min ! boundary condition
      real(dblprec), dimension(3,natom),intent(in) :: coord ! coordinates for the atoms in the lattice
      real(dblprec) :: cts_2d_y_min
      real(dblprec) x_size, pi, arg_y, y_size,x ,y

      parameter(pi=4._dblprec*datan(1._dblprec)) ! the best way to define pi

      y_size = maxval(coord(2,:)) ! size of the lattice in the y direction
      x_size = maxval(coord(1,:)) ! size of the lattice in the x direction

      x = coord(1,i)
      y = coord(2,i)

      arg_y=n*pi/x_size

      ! f_y_min, deals with f(x,0)=tx_min and all the other boundary conditions set to zero
      cts_2d_y_min = -2*tx_min*(1-(-1)**n)*sinh(arg_y*(y-y_size))*sin(arg_y*x)
      cts_2d_y_min = cts_2d_y_min/(n*pi*sinh(arg_y*y_size))

   end function cts_2d_y_min

   ! this deals with the boundary condition f(x,ymax)= x(txmax-txmin)/xsize + txmin
   function linear_2d_y_max(tx_min,tx_max,ty_min,ty_max,n,coord,natom,i)

      integer, intent(in) :: n, i, natom
      real(dblprec), intent(in) :: tx_min, tx_max, ty_min, ty_max ! boundary conditions
      real(dblprec), dimension(3,natom), intent(in) :: coord ! coordinates for the atoms in the lattice
      real(dblprec) :: linear_2d_y_max
      real(dblprec) :: pi, arg_y, x_size, y_size, x, y

      parameter(pi=4._dblprec*datan(1._dblprec)) ! the best way to define pi

      y_size = maxval(coord(2,:)) ! size of the lattice in the y direction
      x_size = maxval(coord(1,:)) ! size of the lattice in the x direction

      x = coord(1,i)
      y = coord(2,i)

      arg_y=n*pi/x_size

      linear_2d_y_max = 2*(tx_min*(1-2*(-1)**n)-tx_max*(-1)**n)
      linear_2d_y_max = linear_2d_y_max*sinh(arg_y*y)*sin(arg_y*x)
      linear_2d_y_max = linear_2d_y_max/(n*pi*sinh(arg_y*y_size))

   end function linear_2d_y_max

   ! this deals with the boundary condition f(x,0)= x(txmax-txmin)/xsize + txmin
   function linear_2d_y_min(tx_min,tx_max,ty_min,ty_max,n,coord,natom,i)

      integer, intent(in) :: n, i, natom
      real(dblprec), intent(in) :: tx_min, tx_max, ty_min, ty_max ! boundary conditions
      real(dblprec), dimension(3,natom), intent(in) :: coord ! coordinates for the atoms in the lattice
      real(dblprec) :: linear_2d_y_min
      real(dblprec) :: pi, arg_y, x_size, y_size, x, y

      parameter(pi=4._dblprec*datan(1._dblprec)) ! the best way to define pi

      y_size = maxval(coord(2,:)) ! size of the lattice in the y direction
      x_size = maxval(coord(1,:)) ! size of the lattice in the x direction

      x = coord(1,i)
      y = coord(2,i)

      arg_y=n*pi/x_size

      linear_2d_y_min = -2*(tx_min*(1-2*(-1)**n)-tx_max*(-1)**n)
      linear_2d_y_min = linear_2d_y_min*sinh(arg_y*(y-y_size))*sin(arg_y*x)
      linear_2d_y_min = linear_2d_y_min/(n*pi*sinh(arg_y*y_size))

   end function linear_2d_y_min

   ! this deals with the boundary condition f(0,y)= y(tymax-tymin)/ysize + tymin
   function linear_2d_x_min(tx_min,tx_max,ty_min,ty_max,n,coord,natom,i)

      integer, intent(in) :: n, i, natom
      real(dblprec), intent(in) :: tx_min, tx_max, ty_min, ty_max ! boundary conditions
      real(dblprec), dimension(3,natom), intent(in) :: coord ! coordinates for the atoms in the lattice
      real(dblprec) :: linear_2d_x_min
      real(dblprec) :: pi, arg_x, x_size, y_size, x, y

      parameter(pi=4._dblprec*datan(1._dblprec)) ! the best way to define pi

      y_size = maxval(coord(2,:)) ! size of the lattice in the y direction
      x_size = maxval(coord(1,:)) ! size of the lattice in the x direction

      x = coord(1,i)
      y = coord(2,i)

      arg_x=n*pi/y_size

      linear_2d_x_min = -2*(ty_min*(1-2*(-1)**n)-ty_max*(-1)**n)
      linear_2d_x_min = linear_2d_x_min*sinh(arg_x*(x-x_size))*sin(arg_x*y)
      linear_2d_x_min = linear_2d_x_min/(n*pi*sinh(arg_x*x_size))

   end function linear_2d_x_min

   ! this deals with the boundary condition f(x_max,y)= y(tymax-tymin)/ysize + tymin
   function linear_2d_x_max(tx_min,tx_max,ty_min,ty_max,n,coord,natom,i)

      integer, intent(in) :: n, i, natom
      real(dblprec), intent(in) :: tx_min, tx_max, ty_min, ty_max ! boundary conditions
      real(dblprec), dimension(3,natom), intent(in) :: coord ! coordinates for the atoms in the lattice
      real(dblprec) :: linear_2d_x_max
      real(dblprec) :: pi, arg_x, x_size, y_size, x, y

      parameter(pi=4._dblprec*datan(1._dblprec)) ! the best way to define pi

      y_size = maxval(coord(2,:)) ! size of the lattice in the y direction
      x_size = maxval(coord(1,:)) ! size of the lattice in the x direction

      x = coord(1,i)
      y = coord(2,i)

      arg_x=n*pi/y_size

      linear_2d_x_max = 2*(ty_min*(1-2*(-1)**n)-ty_max*(-1)**n)
      linear_2d_x_max = linear_2d_x_max*sinh(arg_x*x)*sin(arg_x*y)
      linear_2d_x_max = linear_2d_x_max/(n*pi*sinh(arg_x*x_size))

   end function linear_2d_x_max

   ! subroutine to choose the function which will fill up the lattice
   function lpfunction(ty_min,ty_max,tx_min,tx_max,coord,natom,i,n,grad,bounds,temp)

      integer, intent(in) :: n, i, natom
      integer ::  dimn, j
      real(dblprec), intent(in) :: tx_min, tx_max, ty_min, ty_max, temp ! boundary conditions
      real(dblprec), dimension(3,natom), intent(in) :: coord ! coordinates for the atoms in the lattice
      real(dblprec) :: f_x_max, f_x_min, f_y_min, f_y_max, f_z_min, f_z_max, lpfunction ! functions relating to each of the 6 boundary conditions
      character(len=1) :: grad ! varibale for the gradient of the system
      character(len=20), dimension(6), intent(in) :: bounds ! array that states the boundary conditions of the sample

      ! if there is no gradient the function is just a constant
      if(grad.eq.'n') then
         lpfunction = temp
         ! if there is a gradient the value of the function will be given by the values stored at the bounds array
      else
         do j=1,6
            ! counting how many boundaries are 'on' n states that the boudary is off therefore indicating the dimension of the system
            if(bounds(j).ne.'n') dimn=dimn+1
         end do

         ! if there are two or less values of boundaries different from n 1d solutions are chosen
         if(dimn.le.2) then
            if(bounds(1).eq.'constant') then
               f_x_min = linear_1d(natom,temp_high_x,temp_low_x,coord,i)
            else if(bounds(2).eq.'constant') then
               f_x_min = linear_1d(natom,temp_high_x,temp_low_x,coord,i)
            else if((bounds(1).eq.'constant').and.(bounds(2).eq.'constant')) then
               f_x_min = linear_1d(natom,temp_high_x,temp_low_x,coord,i)
            else if(bounds(1).eq.'step') then
               f_x_min = step_x(barr_size,natom,temp_high_x,temp_low_x,coord,i)
            else if(bounds(2).eq.'step') then
               f_x_min = step_x(barr_size,natom,temp_high_x,temp_low_x,coord,i)
            else if((bounds(1).eq.'step').and.(bounds(2).eq.'step')) then
               f_x_min = step_x(barr_size,natom,temp_high_x,temp_low_x,coord,i)

               f_x_max = 0.0_dblprec
               f_y_min = 0.0_dblprec
               f_y_max = 0.0_dblprec
               f_z_max = 0.0_dblprec
               f_z_min = 0.0_dblprec

            end if
         end if

         ! if there are only two values of bounds which values equal to n then the temperature array is for 2d systems
         if(dimn.eq.4) then
            select case(bounds(1))
            case('constant')
               f_x_min = cts_2d_x_min(ty_min,n,coord,natom,i)
            case('linear')
               f_x_min = linear_2d_x_min(tx_min,tx_max,ty_min,ty_max,n,coord,natom,i)
            end select
            select case(bounds(2))
            case('constant')
               f_x_max =  cts_2d_x_max(ty_max,n,coord,natom,i)
            case('linear')
               f_x_max = linear_2d_x_max(tx_min,tx_max,ty_min,ty_max,n,coord,natom,i)
            end select
            select case(bounds(3))
            case('constant')
               f_y_min = cts_2d_y_min(tx_min,n,coord,natom,i)
            case('linear')
               f_y_min = linear_2d_y_min(tx_min,tx_max,ty_min,ty_max,n,coord,natom,i)
            end select
            select case(bounds(4))
            case('constant')
               f_y_max = cts_2d_y_max(tx_max,n,coord,natom,i)
            case('linear')
               f_y_max = linear_2d_y_max(tx_min,tx_max,ty_min,ty_max,n,coord,natom,i)
            end select
            f_z_max = 0.0_dblprec
            f_z_min = 0.0_dblprec
         end if

         ! if all the boundaries are different from n the temperature array will be filled with solutions of the 3d laplace equation
         !	if(dimn.eq.6)
         !
         !		select case(bounds(1))
         !
         !		case('constant')
         !			f_x_min = cts_3d_x_min(ty_min,n,coord,natom,i)
         !		case('linear')
         !			f_x_min = linear_3d_x_min(tx_min,tx_max,ty_min,ty_max,n,coord,natom,i)
         !
         !		end select case
         !
         !		select case(bounds(2))
         !
         !		case('constant')
         !			f_x_max =  cts_3d_x_max(ty_max,n,coord,natom,i)
         !		case('linear')
         !			f_x_max = linear_3d_x_max(tx_min,tx_max,ty_min,ty_max,n,coord,natom,i)
         !
         !		end select case
         !
         !		select case(bounds(3))
         !
         !		case('constant')
         !			f_y_min = cts_3d_y_min(tx_min,n,coord,natom,i)
         !		case('linear')
         !			f_y_min = linear_3d_y_min(tx_min,tx_max,ty_min,ty_max,n,coord,natom,i)
         !
         !		end select case
         !
         !		select case(bounds(4))
         !
         !		case('constant')
         !			f_y_max = cts_3d_y_max(tx_max,n,coord,natom,i)
         !		case('linear')
         !			f_y_max = linear_3d_y_max(tx_min,tx_max,ty_min,ty_max,n,coord,natom,i)
         !
         !		end select case
         !
         !		select case(bounds(5))
         !
         !		case('constant')
         !			f_z_min = cts_3d_z_min(tx_max,n,coord,natom,i)
         !		case('linear')
         !			f_z_min = linear_3d_z_min(tx_min,tx_max,ty_min,ty_max,n,coord,natom,i)
         !
         !		end select case
         !
         !		select case(bounds(6))
         !
         !		case('constant')
         !			f_z_max = cts_3d_z_max(tx_max,n,coord,natom,i)
         !		case('linear')
         !			f_z_max = linear_3d_z_max(tx_min,tx_max,ty_min,ty_max,n,coord,natom,i)
         !
         !		end select case
         !
         !	end if
         lpfunction = f_x_max + f_x_min + f_y_max + f_y_min + f_z_min + f_z_max

      end if

   end function lpfunction

end module temperature
