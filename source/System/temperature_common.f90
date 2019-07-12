!-------------------------------------------------------------------------------
! MODULE temperature_common
!> @brief Module with helper routines for the tempearture gradients calculations
!> @author Jonathan Chico, Anders Bergman
!> @copyright
!> GNU Public License.
!> @todo generalize to all kind of system shapes
!-------------------------------------------------------------------------------
module temperature_common
   use Parameters
   use Profiling
   use Geometry
   use Constants

   implicit none

   ! fdm variables
   integer :: temp_max_no_neigh
   integer :: temp_max_no_equiv
   integer, dimension(:), allocatable         :: temp_nlistsize !< size of neighbour list for temperature grid
   integer, dimension(:,:), allocatable       :: temp_nlist !< list of neighbours for temperature grid
   integer, dimension(:,:), allocatable       :: temp_nmdim !< dimension of neighbour map for temperature grid
   integer, dimension(:,:,:), allocatable     :: temp_nm  !< neighbour map for temeprature grid
   ! meshless
   real(dblprec), dimension(:,:), allocatable :: s_matrix
   real(dblprec), dimension(:), allocatable   :: a_matrix
   integer, dimension(:), allocatable         :: indx
   integer, dimension(:), allocatable         :: ija
   real(dblprec), dimension(:), allocatable   :: p_vec
   real(dblprec), dimension(:), allocatable   :: c_vec
   real(dblprec), dimension(:), allocatable   :: w_vector
   real(dblprec), dimension(:,:), allocatable :: ptp
   real(dblprec), dimension(:,:), allocatable :: u_matrix
   real(dblprec), dimension(:,:), allocatable :: v_matrix
   real(dblprec), dimension(:,:), allocatable :: mom_matrix
   real(dblprec), dimension(:,:), allocatable :: div_mom_matrix
   real(dblprec), dimension(:,:), allocatable :: div2_mom_matrix

   public

contains
   !-----------------------------------------------------------------------------
   ! SUBROUTINE setup_temp_grid
   !> @brief subroutine to create the grid information necessary for the finite
   !> difference method calculation of the temperature gradient
   !> @author Jonathan Chico
   !-----------------------------------------------------------------------------
   subroutine setup_temp_grid(NT,NA,N1,N2,N3,Natom,crys_symm,natom_full,temp_max_no_shells,&
         do_ralloy,temp_nn,atype,atype_ch,acellnumb,C1,C2,C3,BC1,BC2,BC3,Bas,&
         temp_neigh_dist,temp_max_no_neigh,temp_max_no_equiv,temp_nm,temp_nmdim,temp_nlistsize,&
         temp_nlist,count_I1_min,count_I1_max,count_I2_min,count_I2_max,count_I3_min,count_I3_max,&
         borders_I1_min,borders_I1_max,borders_I2_min,borders_I2_max,borders_I3_min,borders_I3_max,&
         border_belong,simid)

      use neighbourmap

      implicit none

      !.. input needed
      integer, intent(in) :: NT  !< number of types of atoms
      integer, intent(in) :: NA  !< number of atoms in one cell
      integer, intent(in) :: N1  !< number of cell repetitions in x direction
      integer, intent(in) :: N2  !< number of cell repetitions in y direction
      integer, intent(in) :: N3  !< number of cell repetitions in z direction
      integer, intent(in) :: Natom      !< number of atoms in system
      integer, intent(in) :: crys_symm  !< symmetry of system (0-3)
      integer, intent(in) :: Natom_full !< number of atoms for full system (=natom if not dilute)
      integer, intent(in) :: temp_max_no_shells !< calculated maximum of neighbours for exchange
      integer, intent(in), optional :: do_ralloy     !< random alloy simulation (0/1)
      integer, dimension(NT), intent(in) :: temp_nn  !< number of neighbour shells
      integer, dimension(Natom), intent(in) :: atype !< type of atom
      integer, dimension(Natom_full), optional ,intent(in) :: atype_ch  !< actual type of atom for dilute system
      integer, dimension(Natom_full), optional ,intent(in) :: acellnumb !< list for translating atom no. in full cell to actual cell
      real(dblprec), dimension(3), intent(in) :: C1 !< first lattice vector
      real(dblprec), dimension(3), intent(in) :: C2 !< second lattice vector
      real(dblprec), dimension(3), intent(in) :: C3 !< third lattice vector
      real(dblprec), dimension(3,NA), intent(in) :: Bas !< coordinates for basis atoms
      real(dblprec), dimension(NT,temp_max_no_shells,3), intent(in) :: temp_neigh_dist !< coordinates for temperature gird neighbours
      character(len=1), intent(in) :: BC1 !< boundary conditions in x-direction
      character(len=1), intent(in) :: BC2 !< boundary conditions in y-direction
      character(len=1), intent(in) :: BC3 !< boundary conditions in z-direction
      character(len=8), intent(in) :: simid !< simulation name

      !.. outputs from the subroutine
      integer, intent(out) :: count_I1_min  !< Number of atoms belonging to the I1=0 border
      integer, intent(out) :: count_I1_max  !< Number of atoms belonging to the I1=N1 border
      integer, intent(out) :: count_I2_min  !< Number of atoms belonging to the I2=0 border
      integer, intent(out) :: count_I2_max  !< Number of atoms belonging to the I2=N2 border
      integer, intent(out) :: count_I3_min  !< Number of atoms belonging to the I3=0 border
      integer, intent(out) :: count_I3_max  !< Number of atoms belonging to the I3=N3 border
      integer, intent(out) :: temp_max_no_neigh !< calculated maximum of neighbours for temperature grid
      integer, intent(out) :: temp_max_no_equiv !< calculated maximum of neighbours in one shell for temperature grid
      integer, dimension(:), allocatable, intent(out) :: borders_I1_max !< list of atoms belonging to the I1=N1 border
      integer, dimension(:), allocatable, intent(out) :: borders_I1_min !< list of atoms belonging to the I1=0 border
      integer, dimension(:), allocatable, intent(out) :: borders_I2_max !< list of atoms belonging to the I2=N2 border
      integer, dimension(:), allocatable, intent(out) :: borders_I2_min !< list of atoms belonging to the I2=0 border
      integer, dimension(:), allocatable, intent(out) :: borders_I3_max !< list of atoms belonging to the I3=N3 border
      integer, dimension(:), allocatable, intent(out) :: borders_I3_min !< list of atoms belonging to the I3=0 border
      integer, dimension(:), allocatable, intent(out) :: temp_nlistsize !< size of neighbour list for temperature grid
      integer, dimension(:,:), allocatable, intent(out) :: temp_nlist   !< neighbour list for the temperature grid
      integer, dimension(:,:), allocatable, intent(out) :: temp_nmdim   !< dimension of neighbour map for the temprature grid
      integer, dimension(:,:,:), allocatable, intent(out) :: temp_nm    !< neighbour map for the temperature grid
      logical, dimension(natom), intent(out) :: border_belong  !< array to identify wether an atom belongs or not in a border

      !.. local scalar variables
      integer :: ncount
      integer :: count_num_border
      integer :: I1,I2,I3,I0
      integer :: i, j, k, l ,i_stat
      logical :: exis
      real(dblprec) :: xmax,xmin,ymax,ymin,zmax,zmin
      real(dblprec) :: tol

      tol =0.0050000_dblprec
      border_belong=.false.

      count_I1_min=0
      count_I1_max=0
      count_I2_min=0
      count_I2_max=0
      count_I3_min=0
      count_I3_max=0
      count_num_border=0

      ! need to introduce a variable to differentiate which atoms from the unit cell
      xmax=maxval(Bas(1,:))
      xmin=minval(Bas(1,:))
      ymax=maxval(Bas(2,:))
      ymin=minval(Bas(2,:))
      zmax=maxval(Bas(3,:))
      zmin=minval(Bas(3,:))

      ! one can use the already defined neighbour lists from the exchange to deal with the fact that the nearest neighbours are needed to calculate the temperature distribution
      call setup_nm(Natom, NT, NA,N1, N2, N3, C1, C2, C3, BC1, BC2, BC3, 1, atype, Bas, &
         temp_max_no_neigh, temp_max_no_shells, temp_max_no_equiv, crys_symm, &
         temp_nn, temp_neigh_dist, temp_nm, temp_nmdim, do_ralloy, Natom_full, acellnumb, atype_ch)

      ! the full list of neighbours for the solution of the laplace of poisson equation
      allocate(temp_nlist(temp_max_no_neigh,Natom))
      call memocc(i_stat,product(shape(temp_nlist))*kind(temp_nlist),'temp_nlist','setup_temp_grid')
      allocate(temp_nlistsize(Natom))
      call memocc(i_stat,product(shape(temp_nlistsize))*kind(temp_nlistsize),'temp_nlistsize','setup_temp_grid')

      do I3=0, N3-1
         do I2=0, N2-1
            do I1=0, N1-1
               do I0=1, NA
                  i=I0+I1*NA+I2*N1*Na+I3*N2*N1*NA
                  ncount = 1
                  ! shells
                  do k=1, temp_nn(atype(i))
                     ! sites in shell
                     do j=1, temp_nmdim(k,i)
                        ! existing coupling
                        if (temp_nm(i,k,j)>0) then
                           exis=.false.
                           do l=1,ncount-1
                              if (temp_nlist(l,i)==temp_nm(i,k,j)) exis=.true.
                           end do
                           temp_nlist(ncount,i) = temp_nm(i,k,j)
                           ncount = ncount + 1
                        end if
                     end do
                  end do
                  temp_nlistsize(i) = ncount-1
                  ! the border points in the direction given by c1 (lower part)
                  if ((I1==0.and.N1-1.ne.0).and.Bas(1,I0).eq.xmin) then
                     count_I1_min=count_I1_min+1
                     count_num_border=count_num_border+1
                     border_belong(i)=.true.
                     ! the border points in the direction given by c1 (higher part)
                  else if ((I1==N1-1.and.N1-1.ne.0).and.Bas(1,I0).eq.xmax) then
                     count_I1_max=count_I1_max+1
                     count_num_border=count_num_border+1
                     border_belong(i)=.true.
                     ! the border points in the direction given by c2 (lower part)
                  else if ((I2==0.and.N2-1.ne.0).and.Bas(2,I0).eq.ymin) then
                     count_I2_min=count_I2_min+1
                     count_num_border=count_num_border+1
                     border_belong(i)=.true.
                     ! the border points in the direction given by c2 (higher part)
                  else if ((I2==N2-1.and.N2-1.ne.0).and.Bas(2,I0).eq.ymax) then
                     count_I2_max=count_I2_max+1
                     count_num_border=count_num_border+1
                     border_belong(i)=.true.
                     ! the border points in the direction given by c3 (lower part)
                  else if ((I3==0.and.N3-1.ne.0).and.Bas(3,I0).eq.zmin) then
                     count_I3_min=count_I3_min+1
                     count_num_border=count_num_border+1
                     border_belong(i)=.true.
                     ! the border points in the direction given by c3 (higher part)
                  else if ((I3==N3-1.and.N3-1.ne.0).and.Bas(3,I0).eq.zmax) then
                     count_I3_max=count_I3_max+1
                     count_num_border=count_num_border+1
                     border_belong(i)=.true.
                  endif
               enddo
            enddo
         enddo
      enddo

      ! allocation of arrays that deal with the borders of the system
      if (count_I1_max.ne.0) then
         allocate(borders_I1_max(count_I1_max))
         call memocc(i_stat,product(shape(borders_I1_max))*kind(borders_I1_max),'borders_I1_max','temperature_init')
      endif
      if (count_I1_min.ne.0) then
         allocate(borders_I1_min(count_I1_min))
         call memocc(i_stat,product(shape(borders_I1_min))*kind(borders_I1_min),'borders_I1_min','temperature_init')
      endif
      if (count_I2_max.ne.0) then
         allocate(borders_I2_max(count_I2_max))
         call memocc(i_stat,product(shape(borders_I2_max))*kind(borders_I2_max),'borders_I2_max','temperature_init')
      endif
      if (count_I2_min.ne.0) then
         allocate(borders_I2_min(count_I2_min))
         call memocc(i_stat,product(shape(borders_I2_min))*kind(borders_I2_min),'borders_I2_min','temperature_init')
      endif
      if (count_I3_max.ne.0) then
         allocate(borders_I3_max(count_I3_max))
         call memocc(i_stat,product(shape(borders_I3_max))*kind(borders_I3_max),'borders_I3_max','temperature_init')
      endif
      if (count_i3_min.ne.0) then
         allocate(borders_I3_min(count_I3_min))
         call memocc(i_stat,product(shape(borders_I3_min))*kind(borders_I3_min),'borders_I3_min','temperature_init')
      endif

      ! initialize the variables to have edges arrays
      count_I1_min=0
      count_I1_max=0
      count_I2_min=0
      count_I2_max=0
      count_I3_min=0
      count_I3_max=0

      do I3=0, N3-1
         do I2=0, N2-1
            do I1=0, N1-1
               do I0=1, NA
                  i=I0+I1*NA+I2*N1*NA+I3*N2*N1*NA
                  ! the border points in the direction given by c1 (lower part)
                  if ((I1==0.and.N1-1.ne.0).and.(allocated(borders_I1_min)).and.Bas(1,I0).eq.xmin) then
                     count_I1_min=count_I1_min+1
                     borders_I1_min(count_I1_min)=i
                     ! the border points in the direction given by c1 (higher part)
                  else if ((I1==N1-1.and.N1-1.ne.0).and.(allocated(borders_I1_max)).and.Bas(1,I0).eq.xmax) then
                     count_I1_max=count_I1_max+1
                     borders_I1_max(count_I1_max)=i
                     ! the border points in the direction given by c2 (lower part)
                  else if ((I2==0.and.N2-1.ne.0).and.(allocated(borders_I2_min)).and.Bas(2,I0).eq.ymin) then
                     count_I2_min=count_I2_min+1
                     borders_I2_min(count_I2_min)=i
                     ! the border points in the direction given by c2 (higher part)
                  else if ((I2==N2-1.and.N2-1.ne.0).and.(allocated(borders_I2_max)).and.Bas(2,I0).eq.ymax) then
                     count_I2_max=count_I2_max+1
                     borders_I2_max(count_I2_max)=i
                     ! the border points in the direction given by c3 (lower part)
                  else if ((I3==0.and.N3-1.ne.0).and.(allocated(borders_I3_min)).and.Bas(3,I0).eq.zmin) then
                     count_I3_min=count_I3_min+1
                     borders_I3_min(count_I3_min)=i
                     ! the border points in the direction given by c3 (higher part)
                  else if ((I3==N3-1.and.N3-1.ne.0).and.(allocated(borders_I3_max)).and.Bas(3,I0).eq.zmax) then
                     count_i3_max=count_i3_max+1
                     borders_i3_max(count_i3_max)=i
                  endif
               enddo
            enddo
         enddo
      enddo

      ! print the neighbour map to be sure it is correct
      call prn_temp_neigh(natom,temp_max_no_neigh, temp_nlistsize, temp_nlist, simid)

      return

   end subroutine setup_temp_grid

   !-----------------------------------------------------------------------------
   ! SUBROUTINE: temperature_init
   !> @brief Setup the initial temperature distribution, i.e. the boundaries or source terms
   !> either via file or functions activated by input variables
   !> @author Jonathan Chico
   !-----------------------------------------------------------------------------
   subroutine temperature_init(natom,init_temp,count_i1_min,count_i1_max,count_i2_min,count_i2_max,&
         count_i3_min,count_i3_max,borders_i1_min,borders_i1_max,borders_i2_min,borders_i2_max,&
         borders_i3_min,borders_i3_max,eq_type,source_type,loadtemp,i1_min_border,i1_max_border,&
         i2_min_border,i2_max_border,i3_min_border,i3_max_border,temp_max,temp_i1_min,&
         temp_i1_max,temp_i2_min,temp_i2_max,temp_i3_min,temp_i3_max,temp_i1_min_low,temp_i1_min_high,&
         temp_i1_max_low,temp_i1_max_high,temp_i2_min_low,temp_i2_min_high,temp_i2_max_low,&
         temp_i2_max_high,temp_i3_min_low,temp_i3_min_high,temp_i3_max_low,temp_i3_max_high,&
         r_center,sigmatemp,coord,temp_init)

      implicit none

      integer, intent(in) :: Natom !< number of atoms in the system
      integer, intent(in) :: init_temp    !< flag for the type of initial temperature
      integer, intent(in) :: count_I1_min !< number of atoms belonging to the i1=0 border
      integer, intent(in) :: count_I1_max !< number of atoms belonging to the i1=n1 border
      integer, intent(in) :: count_I2_min !< number of atoms belonging to the i2=0 border
      integer, intent(in) :: count_I2_max !< number of atoms belonging to the i1=n2 border
      integer, intent(in) :: count_I3_min !< number of atoms belonging to the i3=0 border
      integer, intent(in) :: count_I3_max !< number of atoms belonging to the i3=n3 border
      integer, dimension(count_I1_min), intent(in) :: borders_I1_min !< list of atoms belonging to the i1=0 border
      integer, dimension(count_I1_max), intent(in) :: borders_I1_max !< list of atoms belonging to the i1=n1 border
      integer, dimension(count_I2_min), intent(in) :: borders_I2_min !< list of atoms belonging to the i2=0 border
      integer, dimension(count_I2_max), intent(in) :: borders_I2_max !< list of atoms belonging to the i2=n2 border
      integer, dimension(count_I3_min), intent(in) :: borders_I3_min !< list of atoms belonging to the i3=0 border
      integer, dimension(count_I3_max), intent(in) :: borders_I3_max !< list of atoms belonging to the i3=n3 border

      character(len=1), intent(in) :: eq_type     !< type of equaition poisson (p) or laplace (l)
      character(len=1), intent(in) :: source_type !< type of source for the poisson equation (currently only gaussian)
      character(len=35), intent(in) :: loadtemp        !< "restart file" for the temperature
      character(len=10), intent(in) :: I1_min_border   !< type of border in the i1=0 border
      character(len=10), intent(in) :: I1_max_border   !< type of border in the i1=n1 border
      character(len=10), intent(in) :: I2_min_border   !< type of border in the i2=0 border
      character(len=10), intent(in) :: I2_max_border   !< type of border in the i2=n2 border
      character(len=10), intent(in) :: I3_min_border   !< type of border in the i3=0 border
      character(len=10), intent(in) :: I3_max_border   !< type of border in the i3=n3 border

      real(dblprec), intent(in) :: temp_max         !< amplitude of the gaussian pulse
      real(dblprec), intent(in) :: temp_I1_min      !< temperature of border i1=0 for constant temp
      real(dblprec), intent(in) :: temp_I1_max      !< temperature of border i1=n1 for constant temp
      real(dblprec), intent(in) :: temp_I2_min      !< temperature of border i2=0 for constant temp
      real(dblprec), intent(in) :: temp_I2_max      !< temperature of border i2=n2 for constant temp
      real(dblprec), intent(in) :: temp_I3_min      !< temperature of border i3=0 for constant temp
      real(dblprec), intent(in) :: temp_I3_max      !< temperature of border i3=n3 for constant temp
      real(dblprec), intent(in) :: temp_I1_min_low  !< minimum tempearture of border i1=0 for linear profile
      real(dblprec), intent(in) :: temp_I1_min_high !< maximum tempearture of border i1=0 for linear profile
      real(dblprec), intent(in) :: temp_I1_max_low  !< minimum tempearture of border i1=n1 for linear profile
      real(dblprec), intent(in) :: temp_I1_max_high !< maximum tempearture of border i1=n1 for linear profile
      real(dblprec), intent(in) :: temp_I2_min_low  !< minimum tempearture of border i2=0 for linear profile
      real(dblprec), intent(in) :: temp_I2_min_high !< maximum tempearture of border i2=0 for linear profile
      real(dblprec), intent(in) :: temp_I2_max_low  !< minimum tempearture of border i2=n2 for linear profile
      real(dblprec), intent(in) :: temp_I2_max_high !< maximum tempearture of border i2=n2 for linear profile
      real(dblprec), intent(in) :: temp_I3_min_low  !< minimum tempearture of border i3=0 for linear profile
      real(dblprec), intent(in) :: temp_I3_min_high !< maximum tempearture of border i3=0 for linear profile
      real(dblprec), intent(in) :: temp_I3_max_low  !< minimum tempearture of border i3=n3 for linear profile
      real(dblprec), intent(in) :: temp_I3_max_high !< maximum tempearture of border i3=n3 for linear profile

      real(dblprec), dimension(3), intent(in) :: r_center    !< position of the center of the gaussian profile
      real(dblprec), dimension(3), intent(in) :: sigmatemp   !< sigma parameters for the gaussian profile
      real(dblprec), dimension(3,Natom), intent(in) :: coord !< atom coordinates

      real(dblprec), dimension(Natom), intent(out) :: temp_init !< final temperature array

      !.. local variables
      integer :: i
      real(dblprec) :: xmax, ymax, zmax
      real(dblprec) :: I1_max, I2_max, I3_max
      real(dblprec), dimension(3) :: gauss_factor_r

      ! the temperature is initialized to zero
      temp_init=0.0_dblprec

      ! initial temperature values, boundaries and source terms
      ! from input file
      if (init_temp.eq.1) then

         write(*,'(2x,a)',advance='no') "read temperature from file"
         call readloadtemp(natom,loadtemp,temp_init)
         ! initial temperature values, boundaries and source terms
         ! calculated from the routine
      else if (init_temp.eq.2) then

         ! fill up source terms for poisson eq
         if (eq_type.eq.'P') then
            ! several type of functions can be chosen as source terms
            if (source_type.eq.'G') then
               ! gaussian source term
               xmax=maxval(coord(1,:))
               ymax=maxval(coord(2,:))
               zmax=maxval(coord(3,:))

               do i=1,Natom
                  gauss_factor_r(1)=((coord(1,i)-r_center(1))**2)/(2*xmax*sigmatemp(1)**2)
                  if (xmax.eq.0.0_dblprec) gauss_factor_r(1)=0.0_dblprec
                  gauss_factor_r(2)=((coord(2,i)-r_center(2))**2)/(2*xmax*sigmatemp(2)**2)
                  if (ymax.eq.0.0_dblprec) gauss_factor_r(2)=0.0_dblprec
                  gauss_factor_r(3)=((coord(3,i)-r_center(3))**2)/(2*xmax*sigmatemp(3)**2)
                  if (zmax.eq.0.0_dblprec) gauss_factor_r(3)=0.0_dblprec
                  temp_init(i)=temp_init(i)+temp_max*exp(-(gauss_factor_r(1)+gauss_factor_r(2)+gauss_factor_r(3)))
               enddo

            endif

         endif

         ! fill boundary conditions
         ! boundary conditions for the direction given by c1 (lower part)
         if (count_I1_min.ne.0) then
            if (I1_min_border.eq.'constant') then
               ! fill up the border with a constat temperature
               do i=1, count_I1_min
                  temp_init(borders_i1_min(i))=temp_i1_min
               enddo
            else if (I1_min_border.eq.'linear') then
               ! fill up the border with a linear temperature profile
               ! calculate the maximum distance in the x-direction
               ! one should be careful this means the i1 is constrained to be the x-direction!
               I1_max=maxval(coord(2,borders_I1_min(:)))
               do i=1, count_I1_min
                  temp_init(borders_I1_min(i))=temp_I1_min_low+&
                     (temp_I1_min_high-temp_I1_min_low)*coord(2,borders_I1_min(i))/I1_max
               enddo
            endif
         endif
         ! boundary conditions for the direction given by c1 (higher part)
         if (count_I1_max.ne.0) then
            if (I1_max_border.eq.'constant') then
               ! fill up the border with a constat temperature
               do i=1, count_I1_max
                  temp_init(borders_I1_max(i))=temp_I1_max
               enddo
            else if (I1_max_border.eq.'linear') then
               ! fill up the border with a linear temperature profile
               ! calculate the maximum distance in the x-direction
               ! one should be careful this means the i1 is constrained to be the x-direction!
               I1_max=maxval(coord(2,borders_I1_max(:)))
               do i=1, count_I1_max
                  temp_init(borders_I1_max(i))=temp_I1_max_low+&
                     (temp_I1_max_high-temp_I1_max_low)*coord(2,borders_I1_min(i))/I1_max
               enddo
            endif
         endif
         ! boundary conditions for the direction given by c2 (lower part)
         if (count_I2_min.ne.0) then
            if (I2_min_border.eq.'constant') then
               ! fill up the border with a constat temperature
               do i=1, count_I2_min
                  temp_init(borders_I2_min(i))=temp_I2_min
               enddo
            else if (I2_min_border.eq.'linear') then
               ! fill up the border with a linear temperature profile
               ! calculate the maximum distance in the x-direction
               ! one should be careful this means the i2 is constrained to be the y-direction!
               I2_max=maxval(coord(2,borders_I2_min(:)))
               do i=1, count_I2_min
                  temp_init(borders_I2_min(i))=temp_I2_min_low+&
                     (temp_I2_min_high-temp_I2_min_low)*coord(2,borders_I2_min(i))/I2_max
               enddo
            endif
         endif
         ! boundary conditions for the direction given by c2 (higher part)
         if (count_I2_max.ne.0) then
            if (I2_max_border.eq.'constant') then
               ! fill up the border with a constat temperature
               do i=1, count_I2_max
                  temp_init(borders_I2_max(i))=temp_I2_max
               enddo
            else if (I2_max_border.eq.'linear') then
               ! fill up the border with a linear temperature profile
               ! calculate the maximum distance in the x-direction
               ! one should be careful this means the i2 is constrained to be the y-direction!
               I2_max=maxval(coord(2,borders_I2_max(:)))
               do i=1, count_I2_max
                  temp_init(borders_I2_max(i))=temp_I2_max_low+&
                     (temp_I2_max_high-temp_I2_max_low)*coord(2,borders_I2_max(i))/I2_max
               enddo
            endif
         endif
         ! boundary conditions for the direction given by c3 (lower part)
         if (count_I3_min.ne.0) then
            if (I3_min_border.eq.'constant') then
               ! fill up the border with a constat temperature
               do i=1, count_I3_min
                  temp_init(borders_I3_min(i))=temp_I3_min
               enddo
            else if (I3_min_border.eq.'linear') then
               ! fill up the border with a linear temperature profile
               ! calculate the maximum distance in the x-direction
               ! one should be careful this means the i3 is constrained to be the z-direction!
               I3_max=maxval(coord(3,borders_I3_min(:)))
               do i=1, count_I3_min
                  temp_init(borders_I3_min(i))=temp_I3_min_low+&
                     (temp_I3_min_high-temp_I3_min_low)*coord(1,borders_I3_min(i))/I3_max
               enddo
            endif
         endif
         ! boundary conditions for the direction given by c3 (higher part)
         if (count_I3_max.ne.0) then
            if (I3_max_border.eq.'constant') then
               ! fill up the border with a constat temperature
               do i=1, count_I3_max
                  temp_init(borders_I3_max(i))=temp_I3_max
               enddo
            else if (I3_max_border.eq.'linear') then
               ! fill up the border with a linear temperature profile
               ! calculate the maximum distance in the x-direction
               ! one should be careful this means the i3 is constrained to be the z-direction!
               I3_max=maxval(coord(3,borders_I3_max(:)))
               do i=1, count_I3_max
                  temp_init(borders_I3_max(i))=temp_I3_max_low+&
                     (temp_I3_max_high-temp_I3_max_low)*coord(3,borders_I3_max(i))/I3_max
               enddo
            endif
         endif
      end if

      return

   end subroutine temperature_init

   !-----------------------------------------------------------------------------
   ! SUBROUTINE: finite_differences
   !> @brief This is the finite diference methods code only use if one has a cubic or semi-cubic geometry
   !> @details This is a finite diference solver to deal with both the laplace and poisson equation,
   !> two types of situations can be handled by the routine, a simple cubic(square) lattice or
   !> interpolations to obtain information about the missing neighbours
   !> @author Jonathan Chico
   !-----------------------------------------------------------------------------
   subroutine finite_differences(Natom,dim_sys,eq_type,temp_solver,temp_max_no_neigh,temp_nlistsize,&
         temp_nlist,grid_type,temp_init,coord,border_belong,simid,temp_array)

      implicit none

      integer, intent(in) :: Natom
      integer, intent(in) :: dim_sys
      integer, intent(in) :: temp_solver
      integer, intent(in) :: temp_max_no_neigh
      integer, dimension(natom), intent(in) :: temp_nlistsize
      integer, dimension(temp_max_no_neigh,Natom), intent(in) :: temp_nlist
      character(len=1), intent(in) :: eq_type
      character(len=1), intent(in) :: grid_type
      character(len=8), intent(in) :: simid
      real(dblprec), dimension(Natom), intent(in) :: temp_init
      real(dblprec), dimension(3,Natom), intent(in) :: coord

      logical, dimension(Natom), intent(in) :: border_belong

      real(dblprec), dimension(Natom), intent(out) :: temp_array

      ! this is the regular grid setup (sc or square lattice)
      if (grid_type.eq.'R') then

         if (temp_solver.eq.1) then
            ! succesive over-relaxation solver
            call SOR_solver(natom,dim_sys,temp_max_no_neigh,temp_nlistsize,temp_nlist,&
               eq_type,temp_init,coord,border_belong,simid,temp_array)
         else if (temp_solver.eq.2) then
            ! gauss-seidel solver
            !        call gauss_seidel()
         endif

      endif

      return

   end subroutine finite_differences

   !-----------------------------------------------------------------------------
   ! SUBROUTINE: SOR_solver
   !> @brief subroutine to use a succesive over relaxation scheme based in the numerical recipes
   !> this solver only works for a regular cubic(square) grid
   !> @author Jonathan Chico
   !-----------------------------------------------------------------------------
   subroutine SOR_solver(Natom,dim_sys,temp_max_no_neigh,temp_nlistsize,temp_nlist,&
         eq_type,temp_init,coord,border_belong,simid,temp_array)

      implicit none

      integer, intent(in) :: Natom   !< number of atoms in the system
      integer, intent(in) :: dim_sys !< dimension of the system
      integer, intent(in) :: temp_max_no_neigh !< calculated maximum of neighbours for temperature grid
      integer, dimension(natom), intent(in) :: temp_nlistsize !< size of neighbour list for temperature grid
      integer, dimension(temp_max_no_neigh,natom), intent(in) :: temp_nlist !< neighbour list for the temperature grid
      character(len=1), intent(in) :: eq_type !< type of equation to solve poisson (p) or laplace (l)
      character(len=8), intent(in) :: simid   !< simulation name
      real(dblprec), dimension(Natom), intent(in) :: temp_init !< initial value of the temperature + source term
      real(dblprec), dimension(3,Natom), intent(in) :: coord !< coordinates of the atoms
      logical, dimension(Natom), intent(in) :: border_belong !< logical array to see if an atom is in the border

      real(dblprec), dimension(Natom), intent(inout) :: temp_array  !< temperature array

      integer :: max_its
      real(dblprec) :: eps    !< the tolerance level of the solution
      real(dblprec) :: rjac   !< the jacobbi radius of the matrix
      real(dblprec) :: xmax   !< maximum value of the coordinates is x-direction
      real(dblprec) :: ymax   !< maximum value of the coordinates is y-direction
      real(dblprec) :: zmax   !< maximum value of the coordinates is z-direction
      real(dblprec) :: omega  !< this is the relaxation parameter
      real(dblprec) :: errmax !< maximum error
      real(dblprec) :: sum_temp    ! temporal summation of the temperature
      real(dblprec) :: source_tot  ! measurement of the source terms
      real(dblprec), dimension(natom) :: res ! residue

      parameter (eps=1.0e-4)
      parameter (max_its=1e5)

      character(len=30) :: filn
      integer :: i,j,k,fac

      k=1
      fac=0
      res=0.0_dblprec
      rjac=0.0_dblprec
      omega=1.0_dblprec
      errmax=1.0_dblprec
      sum_temp=0.0_dblprec
      temp_array=0.0_dblprec
      xmax=maxval(coord(1,:))
      ymax=maxval(coord(2,:))
      zmax=maxval(coord(3,:))

      ! save in temp_array the values from the borders
      do i=1, Natom
         if (border_belong(i)) temp_array(i)=temp_init(i)
      enddo

      ! calculate the jacobi radius of the solution in accordance to the numerical recipes
      if (xmax.gt.0_dblprec) then
         rjac=cos(pi/xmax)
         fac=fac+1
      else if (ymax.gt.0_dblprec) then
         rjac=rjac+cos(pi/ymax)
         fac=fac+1
      else if (zmax.gt.0_dblprec) then
         rjac=rjac+cos(pi/zmax)
         fac=fac+1
      endif
      rjac=rjac/fac

      ! file for the convergence of the sor method
      write(filn,'(''sor.'',a8,''.out'')') simid
      open(ofileno,file=filn,position="append")

      do while ((k.le.max_its).and.(errmax.gt.eps))
         do i=1, natom
            ! one only considers the atoms which are not in the border
            if (.not.border_belong(i)) then
               sum_temp=0.0_dblprec
               ! obtain the temperature in atom i by averaging over its neighbours
               do j=1,temp_nlistsize(i)
                  sum_temp=sum_temp+temp_array(temp_nlist(j,i))
               enddo
               ! calculate the residue by taking into account the source terms
               res(i)=abs(sum_temp-temp_array(i)*dim_sys*2.0_dblprec+temp_init(i))
               ! calculate the temperature using the relaxation parameter to increase cpnvergence
               temp_array(i)=(1.0_dblprec-omega)*temp_array(i)+omega*(1.0_dblprec/(dim_sys*2.0_dblprec))*(sum_temp+temp_init(i))
               if (k.eq.1) then
                  source_tot=source_tot+abs(temp_init(i))
               endif
            endif
         enddo

         ! change the relaxation parameter to increase convergence speed
         if (k.eq.1) then
            omega=1.0_dblprec/(1.0_dblprec-0.5_dblprec*rjac**2)
         else
            omega=1.0_dblprec/(1.0_dblprec-(0.25_dblprec*rjac**2)*omega)
         endif

         ! calculation of the error depends on wether one solves the poisson or laplace eq.
         if (eq_type.eq.'P') then
            errmax=maxval(res)/source_tot
         else if (eq_type.eq.'L') then
            errmax=maxval(res)
         endif

         ! write the convergence info to file
         write(ofileno,'(i8,es16.8,es16.8)') k, omega, errmax

         k=k+1

         ! if the maximum number of iterations has been done and the error is too high the simulation stops
         if (k.eq.max_its.and.errmax.gt.eps) then
            write(*,*) 'ERROR: No convergence in SOR solver for temperature'
            stop
         endif
      enddo

      close(ofileno)

      return

   end subroutine sor_solver

   !-----------------------------------------------------------------------------
   ! SUBROUTINE: prn_temp_neigh
   !> Print the neighbours for the generation of the temperature gradient
   !-----------------------------------------------------------------------------
   subroutine prn_temp_neigh(natom, temp_max_no_neigh, temp_nlistsize, temp_nlist,simid)
      !
      !.. implicit declarations
      implicit none

      integer, intent(in) :: natom !< number of atoms in system
      integer, intent(in) :: temp_max_no_neigh !< calculated maximum of neighbours for exchange
      integer, dimension(natom), intent(in) :: temp_nlistsize !< size of neighbour list for heisenberg exchange couplings
      integer, dimension(temp_max_no_neigh, natom), intent(in) :: temp_nlist !< neighbour list for heisenberg exchange couplings
      character(len=8), intent(in) :: simid !< simulation name

      !.. local variables
      integer :: i
      character(len=30) :: filn

      !.. executable statements
      write(filn,'(''temp_neigh.'',a8,''.out'')') simid
      open(45, file=filn,position="append")

      ! print neighbor list - after sort
      write (45,*) "sorted data from heisge0"
      do i=1,natom
         write (45,*) "----------------------------------"
         write (45,10001) i,temp_nlistsize(i)
         write (45,10002) temp_nlist(1:temp_nlistsize(i),i)
      end do
      close(45)

      10001 format ("atom=",i8,4x,"no neigh=",i7)
      10002 format ("            ",1x,5i6)
      10003 format (5es16.8)
      10004 format (9es16.8)

   end subroutine prn_temp_neigh

   !-----------------------------------------------------------------------------
   ! SUBROUTINE: readloadtemp
   !> @brief Read the temperature from a file (either for boundaries or complete temperature)
   !> @author Jonathan Chico
   !-----------------------------------------------------------------------------
   subroutine readloadtemp(natom,loadtemp,temp_array)

      implicit none

      integer, intent(in) :: Natom
      character(len=20), intent(in) :: loadtemp
      !character(len=35), intent(in) :: loadtemp
      real(dblprec), dimension(Natom) :: temp_array

      integer :: i,k,ios
      logical :: exists

      print *,loadtemp
      inquire(file=loadtemp,exist=exists)
      if (exists) then
         open(ifileno, iostat=ios,file=loadtemp,status="old")
         do i=1, natom
            read(ifileno,*) k, temp_array(i)
         enddo
         close(ifileno)
      else
         write(*,*) 'ERROR: load_temp file ', trim(adjustl(loadtemp)), ' does not exist.'
         stop
      endif

      return

   end subroutine readloadtemp

   subroutine mean_least_squares_method(natom,weight_type,dim_base_poly,c_fac,&
         rmax,temp_max_no_neigh,temp_nlistsize,temp_nlist,coord,border_belong)

      implicit none

      integer, intent(in) :: natom  ! number of atoms of the system
      integer, intent(in) :: weight_type   ! weight function type
      integer, intent(in) :: dim_base_poly ! dimension of the polynomial basis
      integer, intent(in) :: temp_max_no_neigh ! calculated maximum of neighbours for exchange
      integer, dimension(natom), intent(in) :: temp_nlistsize ! size of neighbour list for heisenberg exchange couplings
      integer, dimension(temp_max_no_neigh, natom), intent(in) :: temp_nlist ! neighbour list for heisenberg exchange couplings

      logical, dimension(natom), intent(in) :: border_belong

      real(dblprec), intent(in) :: rmax  ! radius of the support of the weight, maximum neigh_dist
      real(dblprec), intent(in) :: c_fac ! factor for the gaussian weight
      real(dblprec), dimension(3,natom), intent(in) :: coord

      real(dblprec) :: phi_mls ! the mls shape functions
      real(dblprec) :: div_phi ! the laplacian of the mls shape functions
      real(dblprec) :: tot_phi
      real(dblprec) :: tot_div2_phi
      real(dblprec), dimension(3) :: tot_div2_phi_temp

      integer :: d_par, index_num, coord_index, shape_index
      integer :: i_stat, i_all, node, neigh_nodes, cart_coord

      ! allocate the mxm matrix that contains the polynomial basis
      allocate(ptp(dim_base_poly,dim_base_poly),stat=i_stat)
      call memocc(i_stat,product(shape(ptp))*kind(ptp),'ptp','mean_least_squares_methods')

      ! allocate the mxm moment matrix
      allocate(mom_matrix(dim_base_poly,dim_base_poly),stat=i_stat)
      call memocc(i_stat,product(shape(mom_matrix))*kind(mom_matrix),'mom_matrix','mean_least_squares_methods')

      ! allocate the mxm moment matrix
      allocate(u_matrix(dim_base_poly,dim_base_poly),stat=i_stat)
      call memocc(i_stat,product(shape(u_matrix))*kind(u_matrix),'u_matrix','mean_least_squares_methods')

      ! allocate the mxm moment matrix
      allocate(v_matrix(dim_base_poly,dim_base_poly),stat=i_stat)
      call memocc(i_stat,product(shape(v_matrix))*kind(v_matrix),'v_matrix','mean_least_squares_methods')

      ! allocate the mxm moment matrix
      allocate(w_vector(dim_base_poly),stat=i_stat)
      call memocc(i_stat,product(shape(w_vector))*kind(w_vector),'w_vector','mean_least_squares_methods')

      ! allocate the mxm derivative of the moment matrix
      allocate(div_mom_matrix(dim_base_poly,dim_base_poly),stat=i_stat)
      call memocc(i_stat,product(shape(div_mom_matrix))*kind(div_mom_matrix),'div_mom_matrix','mean_least_squares_methods')

      ! allocate the mxm second order derivative of the moment matrix
      allocate(div2_mom_matrix(dim_base_poly,dim_base_poly),stat=i_stat)
      call memocc(i_stat,product(shape(div2_mom_matrix))*kind(div2_mom_matrix),'div2_mom_matrix','mean_least_squares_methods')

      ! allocate the array for the permutation matrix for the lu decomposition
      allocate(indx(dim_base_poly),stat=i_stat)
      call memocc(i_stat,product(shape(indx))*kind(indx),'indx','mean_least_squares_methods')

      ! allocate the array for the polynomial basis
      allocate(p_vec(dim_base_poly),stat=i_stat)
      call memocc(i_stat,product(shape(p_vec))*kind(p_vec),'p_vec','mean_least_squares_methods')

      ! allocate the array for the bacwards substituion
      allocate(c_vec(dim_base_poly),stat=i_stat)
      call memocc(i_stat,product(shape(c_vec))*kind(c_vec),'c_vec','mean_least_squares_methods')

      open(75,file='stressmatrix.out',position='append')

      ! allocate the stiffness matrix
      allocate(s_matrix(natom,natom),stat=i_stat)
      call memocc(i_stat,product(shape(s_matrix))*kind(s_matrix),'s_matrix','mean_least_squares_methods')

      do index_num=1, natom
         ! loop over all the nodes in the system and calculate their shape functions
         do shape_index=1, natom

            write(*,*) shape_index, index_num
            ! assembly of the stiffness matrix
            if ( border_belong(index_num) ) then
               ! for the atoms at the border one considers only the shape function
               do coord_index=1, natom
                  ! the entries in the border for the stiffness matrix are given by
                  ! kij=sum_k ni(xk)*nj(xk)
                  call calc_phi_mls(natom,shape_index,coord_index,weight_type,dim_base_poly,&
                     temp_max_no_neigh,temp_nlistsize,temp_nlist,rmax,c_fac,coord,phi_mls,c_vec,w_vector,u_matrix,v_matrix)

                  tot_phi=phi_mls

                  call calc_phi_mls(natom,index_num,coord_index,weight_type,dim_base_poly,&
                     temp_max_no_neigh,temp_nlistsize,temp_nlist,rmax,c_fac,coord,phi_mls,c_vec,w_vector,u_matrix,v_matrix)

                  tot_phi=tot_phi*phi_mls

                  s_matrix(index_num,shape_index) = s_matrix(index_num,shape_index)+tot_phi

               enddo

            else if ( .not.border_belong(index_num) ) then
               ! for the atoms inside the domain one considers the derivatives of the shape function

               tot_div2_phi=0.0_dblprec
               tot_div2_phi_temp=0.0_dblprec
               do coord_index=1, natom
                  ! the entried in the domain for the stiffness matrix are given by
                  ! kij=sum_k d(ni(xk))*d(nj(xk)) where d is the differential operator

                  call calc_phi_mls(natom,shape_index,coord_index,weight_type,dim_base_poly,&
                     temp_max_no_neigh,temp_nlistsize,temp_nlist,rmax,c_fac,coord,phi_mls,c_vec,w_vector,u_matrix,v_matrix)

                  do cart_coord=1,3

                     tot_div2_phi_temp(cart_coord)=tot_div2_phi_temp(cart_coord)+div2_phi_mls(natom,cart_coord,shape_index,coord_index,weight_type,&
                        dim_base_poly,temp_max_no_neigh,temp_nlistsize,temp_nlist,rmax,c_fac,coord)

                  enddo

                  call calc_phi_mls(natom,index_num,coord_index,weight_type,dim_base_poly,&
                     temp_max_no_neigh,temp_nlistsize,temp_nlist,rmax,c_fac,coord,phi_mls,c_vec,w_vector,u_matrix,v_matrix)

                  do cart_coord=1,3

                     tot_div2_phi=tot_div2_phi+tot_div2_phi_temp(cart_coord)*div2_phi_mls(natom,cart_coord,index_num,coord_index,weight_type,&
                        dim_base_poly,temp_max_no_neigh,temp_nlistsize,temp_nlist,rmax,c_fac,coord)
                  enddo

               enddo
               s_matrix(index_num,shape_index) = s_matrix(index_num,shape_index)+ tot_div2_phi
            endif
            write(75,'(i8,i8,es16.8)') index_num,shape_index,s_matrix(index_num,shape_index)
         enddo
      enddo

      close(75)

      i_all=-product(shape(ptp))*kind(ptp)
      deallocate(ptp,stat=i_stat)
      call memocc(i_stat,i_all,'ptp','mean_least_squares_methods')

      i_all=-product(shape(mom_matrix))*kind(mom_matrix)
      deallocate(mom_matrix,stat=i_stat)
      call memocc(i_stat,i_all,'mom_matrix','mean_least_squares_methods')

      i_all=-product(shape(div_mom_matrix))*kind(div_mom_matrix)
      deallocate(div_mom_matrix,stat=i_stat)
      call memocc(i_stat,i_all,'div_mom_matrix','mean_least_squares_methods')

      i_all=-product(shape(div2_mom_matrix))*kind(div2_mom_matrix)
      deallocate(div2_mom_matrix,stat=i_stat)
      call memocc(i_stat,i_all,'div2_mom_matrix','mean_least_squares_methods')

      i_all=-product(shape(indx))*kind(indx)
      deallocate(indx,stat=i_stat)
      call memocc(i_stat,i_all,'indx','mean_least_squares_methods')

      i_all=-product(shape(c_vec))*kind(c_vec)
      deallocate(c_vec,stat=i_stat)
      call memocc(i_stat,i_all,'c_vec','mean_least_squares_methods')

      i_all=-product(shape(p_vec))*kind(p_vec)
      deallocate(p_vec,stat=i_stat)
      call memocc(i_stat,i_all,'p_vec','mean_least_squares_methods')


   end subroutine mean_least_squares_method

   ! subroutine to calculate the shape functions for the mls
   subroutine calc_phi_mls(natom,shape_index,coord_index,weight_type,dim_base_poly,&
         temp_max_no_neigh,temp_nlistsize,temp_nlist,rmax,c_fac,coord,phi_mls,c_vec,w_vector,u_matrix,v_matrix)

      implicit none

      integer, intent(in) :: natom
      integer, intent(in) :: shape_index ! the shape index i ni(xk)
      integer, intent(in) :: coord_index ! the coordinate index k in ni(xk)
      integer, intent(in) :: weight_type   ! weight function type
      integer, intent(in) :: dim_base_poly ! dimension of the polynomial basis
      integer, intent(in) :: temp_max_no_neigh ! calculated maximum of neighbours for exchange
      integer, dimension(natom), intent(in) :: temp_nlistsize ! size of neighbour list for heisenberg exchange couplings
      integer, dimension(temp_max_no_neigh, natom), intent(in) :: temp_nlist ! neighbour list for heisenberg exchange couplings

      real(dblprec), intent(in) :: rmax  ! radius of the support of the weight, maximum neigh_dist
      real(dblprec), intent(in) :: c_fac ! factor for the gaussian weight
      real(dblprec), dimension(3,natom), intent(in) :: coord

      real(dblprec), intent(out) :: phi_mls
      real(dblprec), dimension(dim_base_poly), intent(inout) :: c_vec

      real(dblprec), dimension(dim_base_poly), intent(inout) :: w_vector
      real(dblprec), dimension(dim_base_poly,dim_base_poly), intent(inout) :: u_matrix
      real(dblprec), dimension(dim_base_poly,dim_base_poly), intent(inout) :: v_matrix
      integer :: neigh_nodes, i,j, d_par
      real(dblprec), dimension(3) :: rij
      real(dblprec) :: wmax, wmin

      mom_matrix=0.0_dblprec
      c_vec=0.0_dblprec
      p_vec=0.0_dblprec
      w_vector=0.0_dblprec
      v_matrix=0.0_dblprec
      u_matrix=0.0_dblprec

      ! for a certain node calculate the shapse function
      ! for this one need to consider all the nodes inside the support of the current node
      do neigh_nodes=1, temp_nlistsize(coord_index)

         ! distance between the current node and its neighbours
         rij(1)=coord(1,temp_nlist(neigh_nodes,shape_index))-coord(1,shape_index)
         rij(2)=coord(2,temp_nlist(neigh_nodes,shape_index))-coord(2,shape_index)
         rij(3)=coord(3,temp_nlist(neigh_nodes,shape_index))-coord(3,shape_index)

         ptp=0.0_dblprec
         do i=1, dim_base_poly
            p_vec(i)=poly_base(coord(1:3,temp_nlist(neigh_nodes,shape_index)),i)
            do j=1, dim_base_poly
               ptp(i,j)=poly_base(coord(1:3,temp_nlist(neigh_nodes,shape_index)),i)*poly_base(coord(1:3,temp_nlist(neigh_nodes,shape_index)),j)
            enddo
         enddo

         ! calculate the moment matrix
         mom_matrix=mom_matrix+weight_function(weight_type,rmax,c_fac,rij)*ptp

      enddo

      call svd(mom_matrix,u_matrix,w_vector,v_matrix,dim_base_poly,dim_base_poly)

      wmax=maxval(w_vector(:))
      wmin=wmax*1.0d-6

      ! evaluate the monomial base at the evaluation point
      do i=1,dim_base_poly
         p_vec(i)=poly_base(coord(1:3,shape_index),i)
      enddo

      ! perform the backwards substitution for the shape functions
      call svbksb(u_matrix,w_vector,v_matrix,dim_base_poly,dim_base_poly,dim_base_poly,dim_base_poly,dim_base_poly,p_vec,c_vec)

      ! distance between the current node (shape_index) and the evaluation point (coord_index)
      rij(1)=coord(1,shape_index)-coord(1,coord_index)
      rij(2)=coord(2,shape_index)-coord(2,coord_index)
      rij(3)=coord(3,shape_index)-coord(3,coord_index)

      ! output shape function
      phi_mls=dot_product(c_vec,p_vec)*weight_function(weight_type,rmax,c_fac,rij)

      return

   end subroutine calc_phi_mls

   real(dblprec) function div2_phi_mls(natom,cart_coord,shape_index,coord_index,weight_type,&
         dim_base_poly,temp_max_no_neigh,temp_nlistsize,temp_nlist,rmax,c_fac,coord)

      implicit none

      integer, intent(in) :: natom
      integer, intent(in) :: cart_coord
      integer, intent(in) :: shape_index ! the shape index i ni(xk)
      integer, intent(in) :: coord_index ! the coordinate index k in ni(xk)
      integer, intent(in) :: weight_type   ! weight function type
      integer, intent(in) :: dim_base_poly ! dimension of the polynomial basis
      integer, intent(in) :: temp_max_no_neigh ! calculated maximum of neighbours for exchange
      integer, dimension(natom), intent(in) :: temp_nlistsize ! size of neighbour list for heisenberg exchange couplings
      integer, dimension(temp_max_no_neigh, natom), intent(in) :: temp_nlist ! neighbour list for heisenberg exchange couplings

      real(dblprec), intent(in) :: rmax  ! radius of the support of the weight, maximum neigh_dist
      real(dblprec), intent(in) :: c_fac ! factor for the gaussian weight
      real(dblprec), dimension(3,natom), intent(in) :: coord

      integer :: neigh_nodes, i,j
      real(dblprec), dimension(3) :: rij
      real(dblprec), dimension(dim_base_poly) :: b_vector
      real(dblprec), dimension(dim_base_poly) :: b2_vector
      real(dblprec), dimension(dim_base_poly) :: div_p_vec
      real(dblprec), dimension(dim_base_poly) :: div_c_vec
      real(dblprec), dimension(dim_base_poly) :: div2_p_vec
      real(dblprec), dimension(dim_base_poly) :: div2_c_vec

      div_mom_matrix=0.0_dblprec

      ! for a certain node calculate the laplacian of the shape function
      ! for that we first need to calculate the derivative of the moment matrix
      ! for this one need to consider all the nodes inside the support of the current node
      do neigh_nodes=1, temp_nlistsize(shape_index)

         ! distance between the current node and its neighbours
         rij(1)=coord(1,temp_nlist(neigh_nodes,shape_index))-coord(1,shape_index)
         rij(2)=coord(2,temp_nlist(neigh_nodes,shape_index))-coord(2,shape_index)
         rij(3)=coord(3,temp_nlist(neigh_nodes,shape_index))-coord(3,shape_index)

         ptp=0.0_dblprec
         do i=1, dim_base_poly
            p_vec(i)=poly_base(coord(1:3,temp_nlist(neigh_nodes,shape_index)),i)
            do j=1, dim_base_poly
               ptp(i,j)=poly_base(coord(1:3,temp_nlist(neigh_nodes,shape_index)),i)*poly_base(coord(1:3,temp_nlist(neigh_nodes,shape_index)),j)
            enddo
         enddo

         ! calculate the derivative of the moment matrix
         div_mom_matrix=div_mom_matrix+div_weight_function(weight_type,rmax,c_fac,cart_coord,rij)*ptp
         ! calculate the second order derivative of the moment matrix
         div2_mom_matrix=div2_mom_matrix+div2_weight_function(weight_type,rmax,c_fac,cart_coord,rij)*ptp

      enddo

      ! evaluate the derivative of the monomial base at the evaluation point
      do i=1,dim_base_poly
         div_p_vec(i)=div_poly_base(coord(1:3,shape_index),i,cart_coord)
         div2_p_vec(i)=div2_poly_base(coord(1:3,shape_index),i,cart_coord)
      enddo

      ! to perform the backward substitution one must use the following vector instead of just p_vec
      ! -div_mom_matrix*c_vec+div_p_vec

      b_vector=matmul(div_mom_matrix,c_vec)+div_p_vec

      ! there is no need to do the lu decomposition of mom_matrix as it has been done for phi_mls
      ! perform the backward substituion to obtain div_c_vec
      call svbksb(u_matrix,w_vector,v_matrix,dim_base_poly,dim_base_poly,dim_base_poly,dim_base_poly,dim_base_poly,b_vector,div_c_vec)

      ! to perform the backward substitution one must use the following vector instead of just p_vec
      ! div2_p_vec-2*div_mom_matrix*div_c_vec-div2_mom_matrix

      b2_vector=div2_p_vec-2.0_dblprec*matmul(div_mom_matrix,div_c_vec)-matmul(div2_mom_matrix,c_vec)

      ! there is no need to do the lu decomposition of mom_matrix as it has been done for phi_mls
      ! perform the backward substituion to obtain div2_c_vec
      call svbksb(u_matrix,w_vector,v_matrix,dim_base_poly,dim_base_poly,dim_base_poly,dim_base_poly,dim_base_poly,b2_vector,div2_c_vec)

      ! second order derivative of the shape function
      div2_phi_mls=dot_product(div2_c_vec,p_vec)*weight_function(weight_type,rmax,c_fac,rij)+&
         dot_product(div_c_vec,p_vec)*div_weight_function(weight_type,rmax,c_fac,cart_coord,rij)+&
         dot_product(c_vec,p_vec)*div2_weight_function(weight_type,rmax,c_fac,cart_coord,rij)


      return

   end function div2_phi_mls


   ! vector that contains the basis function of the polynomials
   real(dblprec) function poly_base(eval_coord,array_ind)

      implicit none

      integer, intent(in) :: array_ind
      real(dblprec), dimension(3), intent(in) :: eval_coord

      if (array_ind.eq.1) then
         poly_base = 1.0_dblprec
      else if (array_ind.eq.2) then
         poly_base = eval_coord(1)
      else if (array_ind.eq.3) then
         poly_base = eval_coord(2)
      else if (array_ind.eq.4) then
         poly_base = eval_coord(3)
      else if (array_ind.eq.5) then
         poly_base = eval_coord(1)**2
      else if (array_ind.eq.6) then
         poly_base = eval_coord(2)**2
      else if (array_ind.eq.7) then
         poly_base = eval_coord(3)**2
      else if (array_ind.eq.8) then
         poly_base = eval_coord(1)*eval_coord(2)
      else if (array_ind.eq.9) then
         poly_base = eval_coord(2)*eval_coord(3)
      else if (array_ind.eq.10) then
         poly_base = eval_coord(1)*eval_coord(3)
      endif

      return

   end function poly_base


   ! derivatives of the polynomial basis
   real(dblprec) function div_poly_base(eval_coord,array_ind,cart_coord)

      implicit none

      integer, intent(in) :: array_ind
      integer, intent(in) :: cart_coord
      real(dblprec), dimension(3), intent(in) :: eval_coord

      ! derivative with respect to the x-coordinate
      if (cart_coord.eq.1) then

         if (array_ind.eq.1) then
            div_poly_base = 0.0_dblprec
         else if (array_ind.eq.2) then
            div_poly_base = 1.0_dblprec
         else if (array_ind.eq.3) then
            div_poly_base = 0.0_dblprec
         else if (array_ind.eq.4) then
            div_poly_base = 2*eval_coord(1)
         else if (array_ind.eq.5) then
            div_poly_base = 0.0_dblprec
         else if (array_ind.eq.6) then
            div_poly_base = 0.0_dblprec
         else if (array_ind.eq.7) then
            div_poly_base = 0.0_dblprec
         else if (array_ind.eq.8) then
            div_poly_base = eval_coord(2)
         else if (array_ind.eq.9) then
            div_poly_base = 0.0_dblprec
         else if (array_ind.eq.10) then
            div_poly_base = eval_coord(3)
         endif
         ! derivative with respect to the y-coordinate
      else if (cart_coord.eq.2) then

         if (array_ind.eq.1) then
            div_poly_base = 0.0_dblprec
         else if (array_ind.eq.2) then
            div_poly_base = 0.0_dblprec
         else if (array_ind.eq.3) then
            div_poly_base = 1.0_dblprec
         else if (array_ind.eq.4) then
            div_poly_base = 0.0_dblprec
         else if (array_ind.eq.5) then
            div_poly_base = 0.0_dblprec
         else if (array_ind.eq.6) then
            div_poly_base = 2*eval_coord(2)
         else if (array_ind.eq.7) then
            div_poly_base = 0.0_dblprec
         else if (array_ind.eq.8) then
            div_poly_base = eval_coord(1)
         else if (array_ind.eq.9) then
            div_poly_base = eval_coord(3)
         else if (array_ind.eq.10) then
            div_poly_base = 0.0_dblprec
         endif
         ! derivative with respect to the z-coordinate
      else if (cart_coord.eq.3) then

         if (array_ind.eq.1) then
            div_poly_base = 0.0_dblprec
         else if (array_ind.eq.2) then
            div_poly_base = 0.0_dblprec
         else if (array_ind.eq.3) then
            div_poly_base = 0.0_dblprec
         else if (array_ind.eq.4) then
            div_poly_base = 1.0_dblprec
         else if (array_ind.eq.5) then
            div_poly_base = 0.0_dblprec
         else if (array_ind.eq.6) then
            div_poly_base = 0.0_dblprec
         else if (array_ind.eq.7) then
            div_poly_base = 2*eval_coord(3)
         else if (array_ind.eq.8) then
            div_poly_base = 0.0_dblprec
         else if (array_ind.eq.9) then
            div_poly_base = eval_coord(2)
         else if (array_ind.eq.10) then
            div_poly_base = eval_coord(1)
         endif

      endif

      return

   end function div_poly_base


   ! 2nd order derivatives of the polynomial basis
   real(dblprec) function div2_poly_base(eval_coord,array_ind,cart_coord)

      implicit none

      integer, intent(in) :: array_ind
      integer, intent(in) :: cart_coord
      real(dblprec), dimension(3), intent(in) :: eval_coord

      ! derivative with respect to the x-coordinate
      if (cart_coord.eq.1) then

         if (array_ind.eq.1) then
            div2_poly_base = 0.0_dblprec
         else if (array_ind.eq.2) then
            div2_poly_base = 0.0_dblprec
         else if (array_ind.eq.3) then
            div2_poly_base = 0.0_dblprec
         else if (array_ind.eq.4) then
            div2_poly_base = 0.0_dblprec
         else if (array_ind.eq.5) then
            div2_poly_base = 2.0_dblprec
         else if (array_ind.eq.6) then
            div2_poly_base = 0.0_dblprec
         else if (array_ind.eq.7) then
            div2_poly_base = 0.0_dblprec
         else if (array_ind.eq.8) then
            div2_poly_base = 0.0_dblprec
         else if (array_ind.eq.9) then
            div2_poly_base = 0.0_dblprec
         else if (array_ind.eq.10) then
            div2_poly_base = 0.0_dblprec
         endif
         ! derivative with respect to the y-coordinate
      else if (cart_coord.eq.2) then

         if (array_ind.eq.1) then
            div2_poly_base = 0.0_dblprec
         else if (array_ind.eq.2) then
            div2_poly_base = 0.0_dblprec
         else if (array_ind.eq.3) then
            div2_poly_base = 0.0_dblprec
         else if (array_ind.eq.4) then
            div2_poly_base = 0.0_dblprec
         else if (array_ind.eq.5) then
            div2_poly_base = 0.0_dblprec
         else if (array_ind.eq.6) then
            div2_poly_base = 2.0_dblprec
         else if (array_ind.eq.7) then
            div2_poly_base = 0.0_dblprec
         else if (array_ind.eq.8) then
            div2_poly_base = 0.0_dblprec
         else if (array_ind.eq.9) then
            div2_poly_base = 0.0_dblprec
         else if (array_ind.eq.10) then
            div2_poly_base = 0.0_dblprec
         endif
         ! derivative with respect to the z-coordinate
      else if (cart_coord.eq.3) then

         if (array_ind.eq.1) then
            div2_poly_base = 0.0_dblprec
         else if (array_ind.eq.2) then
            div2_poly_base = 0.0_dblprec
         else if (array_ind.eq.3) then
            div2_poly_base = 0.0_dblprec
         else if (array_ind.eq.4) then
            div2_poly_base = 0.0_dblprec
         else if (array_ind.eq.5) then
            div2_poly_base = 0.0_dblprec
         else if (array_ind.eq.6) then
            div2_poly_base = 0.0_dblprec
         else if (array_ind.eq.7) then
            div2_poly_base = 2.0_dblprec
         else if (array_ind.eq.8) then
            div2_poly_base = 0.0_dblprec
         else if (array_ind.eq.9) then
            div2_poly_base = 0.0_dblprec
         else if (array_ind.eq.10) then
            div2_poly_base = 0.0_dblprec
         endif

      endif

      return

   end function div2_poly_base


   ! weight function for the moving least squares method
   real(dblprec) function weight_function(weight_type,rmax,c_fac,rij)

      implicit none

      integer, intent(in) :: weight_type     ! type of weight function

      real(dblprec), intent(in) :: c_fac ! factor for the gaussian function
      real(dblprec), intent(in) :: rmax  ! radius of the support of the weight function
      real(dblprec), dimension(3), intent(in) :: rij ! distance between the ceenter and meassurement point

      real(dblprec) :: inv_rmax
      real(dblprec) :: norm_rij
      real(dblprec) :: rscalar

      inv_rmax=1.0_dblprec/rmax
      rscalar=sqrt(rij(1)**2+rij(2)**2+rij(3)**2)
      norm_rij=rscalar*inv_rmax

      ! decide which kind of weight function will be used
      if (weight_type.eq.1) then
         ! do the cubic spline weight
         if (norm_rij.le.0.5) then
            weight_function = (2.0_dblprec/3.0_dblprec) -4.0_dblprec*(rscalar*inv_rmax)**2
            weight_function = weight_function + 4.0_dblprec*(rscalar*inv_rmax)**3
         else if ((norm_rij.gt.0.5).and.(norm_rij.le.1.0_dblprec)) then
            weight_function = (4.0_dblprec/3.0_dblprec) -4.0_dblprec*rscalar*inv_rmax +4.0_dblprec*(rscalar*inv_rmax)**2
            weight_function = weight_function -(4.0_dblprec/3.0_dblprec)*(rscalar*inv_rmax)**3
         else
            weight_function = 0.0_dblprec
         endif
      else if (weight_type.eq.2) then
         ! do the quadratic spline weight
         if (norm_rij.le.1.0_dblprec) then
            weight_function = 1.0_dblprec -6.0_dblprec*(rscalar*inv_rmax)**2 +8.0_dblprec*(rscalar*inv_rmax)**3
            weight_function = weight_function -3.0_dblprec*(rscalar*inv_rmax)**4
         else
            weight_function = 0.0_dblprec
         endif
      else if (weight_type.eq.3) then
         ! do the gaussian weight
         if (norm_rij.le.1.0_dblprec) then
            weight_function = exp(-(rscalar/c_fac)**2)-exp(-(rmax/c_fac)**2)
            weight_function = weight_function/(1.0_dblprec-exp(-(rmax/c_fac)**2))
         else
            weight_function = 0.0_dblprec
         endif
      endif

      return

   end function weight_function


   ! derivative of the weight function for the mls
   real(dblprec) function div_weight_function(weight_type,rmax,c_fac,cart_coord,rij)

      implicit none

      integer, intent(in) :: cart_coord
      integer, intent(in) :: weight_type

      real(dblprec), intent(in) :: rmax
      real(dblprec), intent(in) :: c_fac
      real(dblprec), dimension(3), intent(in) :: rij

      real(dblprec) :: inv_rmax
      real(dblprec) :: norm_rij
      real(dblprec) :: rscalar

      inv_rmax=1.0_dblprec/rmax
      rscalar=sqrt(rij(1)**2+rij(2)**2+rij(3)**2)
      norm_rij=rscalar*inv_rmax

      div_weight_function=0.0_dblprec

      ! decide which weight function is used
      if (weight_type.eq.1) then
         ! do the cubic spline weight
         if (norm_rij.le.0.5) then
            ! first order derivative with respect to the x-coordinate
            if (cart_coord.eq.1) then
               div_weight_function = 4.0_dblprec*rij(1)*(3.0_dblprec*rscalar*inv_rmax**3)
               div_weight_function = div_weight_function -16.0_dblprec*rij(1)*inv_rmax**2/3.0_dblprec
               ! first order derivative with respect to the y-coordinate
            else if (cart_coord.eq.2) then
               div_weight_function = 4.0_dblprec*rij(2)*(3.0_dblprec*rscalar*inv_rmax**3)
               div_weight_function = div_weight_function -16.0_dblprec*rij(2)*inv_rmax**2/3.0_dblprec
               ! first order derivative with respect to the z-coordinate
            else if (cart_coord.eq.3) then
               div_weight_function = 4.0_dblprec*rij(3)*(3.0_dblprec*rscalar*inv_rmax**3)
               div_weight_function = div_weight_function -16.0_dblprec*rij(3)*inv_rmax**2/3.0_dblprec
            endif
         else if ((norm_rij.gt.0.5).and.(norm_rij.le.1.0_dblprec)) then
            ! first order derivative with respect to the x-coordinate
            if (cart_coord.eq.1) then
               div_weight_function = -4.0_dblprec*rij(1)*(rscalar-rmax)**2
               div_weight_function = div_weight_function*inv_rmax**3/rscalar
               ! first order derivative with respect to the y-coordinate
            else if (cart_coord.eq.2) then
               div_weight_function = -4.0_dblprec*rij(2)*(rscalar-rmax)**2
               div_weight_function = div_weight_function*inv_rmax**3/rscalar
               ! first order derivative with respect to the z-coordinate
            else if (cart_coord.eq.3) then
               div_weight_function = -4.0_dblprec*rij(3)*(rscalar-rmax)**2
               div_weight_function = div_weight_function*inv_rmax**3/rscalar
            endif
         else
            div_weight_function=0.0_dblprec
         endif
      else if (weight_type.eq.2) then
         ! do the quartic spline weight
         if (norm_rij.le.1) then
            ! first order derivative with respect to the x-coordinate
            if (cart_coord.eq.1) then
               div_weight_function = 12.0_dblprec*rij(1)*(2.0_dblprec*rscalar*inv_rmax**3-rscalar**2*inv_rmax**4)
               div_weight_function = div_weight_function -12.0_dblprec*rij(1)*inv_rmax**2
               ! first order derivative with respect to the y-coordinate
            else if (cart_coord.eq.2) then
               div_weight_function = 12.0_dblprec*rij(2)*(2.0_dblprec*rscalar*inv_rmax**3-rscalar**2*inv_rmax**4)
               div_weight_function = div_weight_function -12.0_dblprec*rij(2)*inv_rmax**2
               ! first order derivative with respect to the z-coordinate
            else if (cart_coord.eq.3) then
               div_weight_function = 12.0_dblprec*rij(3)*(2.0_dblprec*rscalar*inv_rmax**3-rscalar**2*inv_rmax**4)
               div_weight_function = div_weight_function -12.0_dblprec*rij(3)*inv_rmax**2
            endif
         else
            div_weight_function=0.0_dblprec
         endif
      else if (weight_type.eq.3) then
         ! do the gaussian weight
         if(weight_type.le.1) then
            ! first order derivative with respect to the x-coordinate
            if (cart_coord.eq.1) then
               div_weight_function = -4.0_dblprec*rij(1)*(rscalar**2)*exp(-(rscalar/c_fac)**2)
               div_weight_function = div_weight_function/((c_fac**2)*(1-exp(-(rmax/c_fac)**2)))
               ! first order derivative with respect to the y-coordinate
            else if (cart_coord.eq.2) then
               div_weight_function = -4.0_dblprec*rij(1)*(rscalar**2)*exp(-(rscalar/c_fac)**2)
               div_weight_function = div_weight_function/((c_fac**2)*(1-exp(-(rmax/c_fac)**2)))
               ! first order derivative with respect to the z-coordinate
            else if (cart_coord.eq.3) then
               div_weight_function = -4.0_dblprec*rij(1)*(rscalar**2)*exp(-(rscalar/c_fac)**2)
               div_weight_function = div_weight_function/((c_fac**2)*(1-exp(-(rmax/c_fac)**2)))
            endif
         else
            div_weight_function=0.0_dblprec
         endif
      endif

      return

   end function div_weight_function


   real(dblprec) function div2_weight_function(weight_type,rmax,c_fac,cart_coord,rij)

      implicit none

      integer, intent(in) :: weight_type
      integer, intent(in) :: cart_coord

      real(dblprec), intent(in) :: rmax
      real(dblprec), intent(in) :: c_fac
      real(dblprec), dimension(3), intent(in) :: rij

      real(dblprec) :: inv_rmax
      real(dblprec) :: norm_rij
      real(dblprec) :: rscalar

      inv_rmax=1.0_dblprec/rmax
      rscalar=sqrt(rij(1)**2+rij(2)**2+rij(3)**2)
      norm_rij=rscalar*inv_rmax

      div2_weight_function=0.0_dblprec

      ! decide which weight function is used
      if (weight_type.eq.1) then
         ! do the cubic spline weight
         if (norm_rij.le.0.5) then
            ! second order derivative with respect to the x-coordinate
            if (cart_coord.eq.1) then
               div2_weight_function = 12.0_dblprec*(inv_rmax**3)*(rij(1)**2)/rscalar
               div2_weight_function = div2_weight_function+12.0_dblprec*(inv_rmax**3)*rscalar
               div2_weight_function = div2_weight_function-16.0_dblprec*inv_rmax**2/3.0_dblprec
               ! second order derivative with respect to the y-coordinate
            else if (cart_coord.eq.2) then
               div2_weight_function = 12.0_dblprec*(inv_rmax**3)*(rij(2)**2)/rscalar
               div2_weight_function = div2_weight_function+12.0_dblprec*(inv_rmax**3)*rscalar
               div2_weight_function = div2_weight_function-16.0_dblprec*inv_rmax**2/3.0_dblprec
               ! second order derivative with respect to the z-coordinate
            else if (cart_coord.eq.3) then
               div2_weight_function = 12.0_dblprec*(inv_rmax**3)*(rij(3)**2)/rscalar
               div2_weight_function = div2_weight_function+12.0_dblprec*(inv_rmax**3)*rscalar
               div2_weight_function = div2_weight_function-16.0_dblprec*inv_rmax**2/3.0_dblprec
            endif
         else if ((norm_rij.gt.0.5).and.(norm_rij.le.1.0_dblprec)) then
            ! second order derivative with respect to the x-coordinate
            if (cart_coord.eq.1) then
               div2_weight_function = -4.0_dblprec*rij(1)**2*inv_rmax**3/rscalar-4.0_dblprec*inv_rmax**3*rscalar
               div2_weight_function = div2_weight_function -4.0_dblprec*inv_rmax/rscalar +4.0_dblprec*rij(1)**2*inv_rmax +8.0_dblprec*inv_rmax**2
               ! second order derivative with respect to the y-coordinate
            else if (cart_coord.eq.2) then
               div2_weight_function = -4.0_dblprec*rij(2)**2*inv_rmax**3/rscalar-4.0_dblprec*inv_rmax**3*rscalar
               div2_weight_function = div2_weight_function -4.0_dblprec*inv_rmax/rscalar +4.0_dblprec*rij(2)**2*inv_rmax +8.0_dblprec*inv_rmax**2
               ! second order derivative with respect to the z-coordinate
            else if (cart_coord.eq.3) then
               div2_weight_function = -4.0_dblprec*rij(3)**2*inv_rmax**3/rscalar-4.0_dblprec*inv_rmax**3*rscalar
               div2_weight_function = div2_weight_function -4.0_dblprec*inv_rmax/rscalar +4.0_dblprec*rij(3)**2*inv_rmax +8.0_dblprec*inv_rmax**2
            endif
         else
            div2_weight_function=0.0_dblprec
         endif
      else if (weight_type.eq.2) then
         if(norm_rij.le.1) then
            ! second order derivative with respect to the x-coordinate
            if (cart_coord.eq.1) then
               div2_weight_function = -12.0_dblprec*inv_rmax**4*(rscalar**2+2.0_dblprec*rij(1)**2)
               div2_weight_function = div2_weight_function+24.0_dblprec*inv_rmax**3*rij(1)**2/rscalar
               div2_weight_function = div2_weight_function-12.0_dblprec*inv_rmax**2
               ! second order derivative with respect to the y-coordinate
            else if (cart_coord.eq.2) then
               div2_weight_function = -12.0_dblprec*inv_rmax**4*(rscalar**2+2.0_dblprec*rij(1)**2)
               div2_weight_function = div2_weight_function+24.0_dblprec*inv_rmax**3*rij(1)**2/rscalar
               div2_weight_function = div2_weight_function-12.0_dblprec*inv_rmax**2
               ! second order derivative with respect to the z-coordinate
            else if (cart_coord.eq.3) then
               div2_weight_function = -12.0_dblprec*inv_rmax**4*(rscalar**2+2.0_dblprec*rij(1)**2)
               div2_weight_function = div2_weight_function+24.0_dblprec*inv_rmax**3*rij(1)**2/rscalar
               div2_weight_function = div2_weight_function-12.0_dblprec*inv_rmax**2
            endif
         else
            div2_weight_function=0.0_dblprec
         endif
      else if (weight_type.eq.3) then
         if(norm_rij.le.1) then
            ! second order derivative with respect to the x-coordinate
            if (cart_coord.eq.1) then
               div2_weight_function = exp(-(rscalar/c_fac)**2)/(1.0_dblprec-exp((rmax/c_fac)**2))*(-8.0_dblprec*(rij(1)/c_fac)**2-4.0_dblprec*(rscalar/c_fac)**2)
               div2_weight_function = div2_weight_function +exp(-(rscalar/c_fac)**2)/(1.0_dblprec-exp((rmax/c_fac)**2))*(16.0_dblprec*rij(1)**2*rscalar**2/c_fac**4)
               ! second order derivative with respect to the y-coordinate
            else if (cart_coord.eq.2) then
               div2_weight_function = exp(-(rscalar/c_fac)**2)/(1.0_dblprec-exp((rmax/c_fac)**2))*(-8.0_dblprec*(rij(2)/c_fac)**2-4.0_dblprec*(rscalar/c_fac)**2)
               div2_weight_function = div2_weight_function +exp(-(rscalar/c_fac)**2)/(1.0_dblprec-exp((rmax/c_fac)**2))*(16.0_dblprec*rij(2)**2*rscalar**2/c_fac**4)
               ! second order derivative with respect to the z-coordinate
            else if (cart_coord.eq.3) then
               div2_weight_function = exp(-(rscalar/c_fac)**2)/(1.0_dblprec-exp((rmax/c_fac)**2))*(-8.0_dblprec*(rij(3)/c_fac)**2-4.0_dblprec*(rscalar/c_fac)**2)
               div2_weight_function = div2_weight_function +exp(-(rscalar/c_fac)**2)/(1.0_dblprec-exp((rmax/c_fac)**2))*(16.0_dblprec*rij(3)**2*rscalar**2/c_fac**4)
            endif
         else
            div2_weight_function=0.0_dblprec
         endif

      endif

      return

   end function div2_weight_function

   ! subroutine to perform the lu decomposition of matrices
   ! this routine is taken from the numerical recipes for fortran 90 the art of programming
   subroutine ludcmp(a_matrix,n_size,indx,d_par)

      implicit none

      integer, intent(in) :: n_size
      real(dblprec), dimension(n_size,n_size), intent(inout) :: a_matrix

      integer, intent(out) :: d_par
      integer, dimension(n_size), intent(out) :: indx

      integer :: imax,i,j,k

      real(dblprec) :: dum
      real(dblprec) :: aamax
      real(dblprec) :: tiny_par
      real(dblprec) :: sum_a_mat
      real(dblprec), dimension(n_size) :: vv

      parameter(tiny_par=1.0e-20)

      d_par=1

      do i=1, n_size
         aamax=0.0_dblprec
         do j=1, n_size
            if ( abs(a_matrix(i,j)).gt.aamax ) aamax=abs(a_matrix(i,j))
         enddo
         !if (aamax.eq.0.0_dblprec) pause 'singular matrix in ludcmp'
         if (aamax.eq.0.0_dblprec) write(*,*) '*** singular matrix in ludcmp ***'
         vv(i)=1.0_dblprec/aamax
      enddo

      do j=1, n_size
         do i=1, j-1
            sum_a_mat=a_matrix(i,j)
            do k=1, i-1
               sum_a_mat=sum_a_mat-a_matrix(i,k)*a_matrix(k,j)
            enddo
            a_matrix(i,j)=sum_a_mat
         enddo
         aamax=0.0_dblprec
         do i=j,n_size
            sum_a_mat=a_matrix(i,j)
            do k=1,j-1
               sum_a_mat=sum_a_mat-a_matrix(i,k)*a_matrix(k,j)
            enddo
            a_matrix(i,j)=sum_a_mat
            dum=vv(i)*abs(sum_a_mat)
            if (dum.ge.aamax) then
               imax=i
               aamax=dum
            endif
         enddo
         if (j.ne.imax) then
            do k=1,n_size
               dum=a_matrix(imax,k)
               a_matrix(imax,k)=a_matrix(j,k)
               a_matrix(j,k)=dum
            enddo
            d_par=-d_par
            vv(imax)=vv(j)
         endif
         indx(j)=imax
         if (a_matrix(j,j).eq.0.0_dblprec) a_matrix(j,j)=tiny_par
         if (j.ne.n_size) then
            dum=1.0_dblprec/a_matrix(j,j)
            do i=j+1,n_size
               a_matrix(i,j)=a_matrix(i,j)*dum
            enddo
         endif
      enddo

      return

   end subroutine ludcmp

   ! subroutine for backsubstitution to solve linear equations a*x=b
   subroutine lubksb(a_matrix,n_size,b_vector,indx,x_vector)

      implicit none

      integer, intent(in) :: n_size
      integer, dimension(n_size), intent(in) :: indx
      real(dblprec), dimension(n_size), intent(inout) :: b_vector
      real(dblprec), dimension(n_size,n_size), intent(in) :: a_matrix

      real(dblprec), dimension(n_size), intent(out) :: x_vector

      integer :: i,ii,j,ll
      real(dblprec) :: sum_a_mat

      ii=0

      do i=1, n_size
         ll=indx(i)
         sum_a_mat=b_vector(ll)
         b_vector(ll)=b_vector(i)
         if (ii.ne.0) then
            do j=ii,i-1
               sum_a_mat=sum_a_mat-a_matrix(i,j)*b_vector(j)
            enddo
         else if (sum_a_mat.ne.0) then
            ii=i
         endif
         b_vector(i)=sum_a_mat
      enddo

      do i=n_size,1,-1
         sum_a_mat=b_vector(i)
         do j=i+1,n_size
            sum_a_mat=sum_a_mat-a_matrix(i,j)*b_vector(j)
         enddo
         x_vector(i)=sum_a_mat/a_matrix(i,i)
      enddo

      return

   end subroutine lubksb


   subroutine linbcg(n_size,nmax,a_matrix,b_vector,x_vector,itol,tol,itmax,iter,err,ija)
      !extracted from the numerical recepies the art of scientific computing for fortran 90
      !solves a·x=b for x, given b of the same length, by the iterative biconjugate gradient method.
      !on input x should be set to an initial guess of the solution (or all zeros);
      !itol is 1,2,3, or 4, specifying which convergence test is applied (see text);
      !itmax is the maximum number of allowed iterations; and tol is the desired convergence tolerance.
      !on output, x is reset to the improved solution, iter is the number of iterations actually taken,
      !and err is the estimated error. the matrix a is referenced only through the user-supplied routines atimes,
      !which computes the product of either a or its transpose on a vector; and asolve, which solves a·x=b or a t·x=b
      !for some preconditioner matrix a(possibly the trivial diagonal part of a).

      implicit none

      integer, intent(in) :: nmax   ! size of the matrix that has the equation
      integer, intent(in) :: itol   ! is the type of convergence test used (1-4)
      integer, intent(in) :: itmax  ! maximum number of iterations allowed
      integer, intent(in) :: n_size ! size of the vectors that define the system of equations
      integer, dimension(nmax), intent(in) :: ija
      real(dblprec), intent(in) :: tol ! the desired tolerance
      real(dblprec), dimension(nmax), intent(in) :: a_matrix   ! matrix that defines the system of equations
      real(dblprec), dimension(n_size), intent(in) :: b_vector ! vector of the right hand side of the equation

      integer, intent(out) :: iter ! the iteration in which the system converges
      real(dblprec), intent(out) :: err ! the error at the end of the convergence
      real(dblprec), dimension(n_size), intent(inout) :: x_vector ! vector that one needs to find

      integer :: j
      real(dblprec) :: eps
      real(dblprec) :: ak,akden,bk,bkden,bknum,bnrm,dxnrm,xnrm,zm1nrm,znrm
      real(dblprec), dimension(nmax) :: p_vector, r_vector, z_vector
      real(dblprec), dimension(nmax) :: pp_vector, rr_vector, zz_vector
      parameter (eps=1.d-14)

      iter=0

      ! this calculates b=a.x and returns b or b=a^t.x returning r_vector
      call atimes(n_size,x_vector,r_vector,nmax,a_matrix,0,ija)

      do j=1,n_size
         r_vector(j)=b_vector(j)-r_vector(j)
         rr_vector(j)=r_vector(j)
      enddo
      !    call atimes(n,r,rr,0)
      if (itol.eq.1) then
         ! calculate the norm of b_vector
         bnrm=snrm(n_size,b_vector,itol)
         ! solves the system z_vector = r_vector.a_matrix^-1
         call asolve(n_size,nmax,a_matrix,r_vector,z_vector)
      else if(itol.eq.2) then
         ! solves the system z_vector = b_vector.a_matrix^-1
         call asolve(n_size,nmax,a_matrix,b_vector,z_vector)
         ! calculate the norm of z_vector
         bnrm=snrm(n_size,z_vector,itol)
         ! solves the system z_vector = r_vector.a_matrix^-1
         call asolve(n_size,nmax,a_matrix,r_vector,z_vector)
      else if(itol.eq.3.or.itol.eq.4) then
         ! solves the system z_vector = b_vector.a_matrix^-1
         call asolve(n_size,nmax,a_matrix,b_vector,z_vector)
         ! calculate the norm of z_vector
         bnrm=snrm(n_size,z_vector,itol)
         ! solves the system z_vector = r_vector.a_matrix^-1
         call asolve(n_size,nmax,a_matrix,r_vector,z_vector)
         ! calculate the norm of z_vector
         znrm=snrm(n_size,z_vector,itol)
      else
         !pause 'illegal itol in linbcg'
         write(*,*) '*** illegal itol in linbcg ***'
      endif

      100  if (iter.le.itmax) then
         iter=iter+1
         ! solves the system zz_vector = rr_vector.a_matrix^-1
         call asolve(n_size,nmax,a_matrix,rr_vector,zz_vector)

         bknum=0.0_dblprec
         do j=1,n_size
            bknum=bknum+z_vector(j)*rr_vector(j)
         enddo

         if(iter.eq.1) then
            do j=1,n_size
               p_vector(j)=z_vector(j)
               pp_vector(j)=zz_vector(j)
            enddo
         else
            bk=bknum/bkden
            do j=1,n_size
               p_vector(j)=bk*p_vector(j)+z_vector(j)
               pp_vector(j)=bk*pp_vector(j)+zz_vector(j)
            enddo
         endif

         bkden=bknum
         ! this calculates b=a.x and returns b or b=a^t.x returning z_vector
         call atimes(n_size,p_vector,z_vector,nmax,a_matrix,0,ija)

         akden=0.0_dblprec
         do j=1, n_size
            akden=akden+z_vector(j)*pp_vector(j)
         enddo
         ak=bknum/akden
         ! this calculates b=a.x and returns b or b=a^t.x returning zz_vector
         call atimes(n_size,pp_vector,zz_vector,nmax,a_matrix,1,ija)

         do j=1,n_size
            x_vector(j)=x_vector(j)+ak*p_vector(j)
            r_vector(j)=r_vector(j)-ak*z_vector(j)
            rr_vector(j)=rr_vector(j)-ak*zz_vector(j)
         enddo

         ! solves the system z_vector = r_vector.a_matrix^-1
         call asolve(n_size,nmax,a_matrix,r_vector,z_vector)

         if (itol.eq.1) then
            ! calculate the norm of r_vector and divide it by the norm of b_vector
            err=snrm(n_size,r_vector,itol)/bnrm
         else if (itol.eq.2) then
            ! calculate the norm of z_vector and divide it by the norm of b_vector
            err=snrm(n_size,z_vector,itol)/bnrm
         else if (itol.eq.3.or.itol.eq.4) then
            zm1nrm=znrm
            ! calculate the norm of z_vector
            znrm=snrm(n_size,z_vector,itol)
            if (abs(zm1nrm-znrm).gt.eps*znrm) then
               ! one calculates the norm of p_vector
               dxnrm=abs(ak)*snrm(n_size,p_vector,itol)
               err=znrm/abs(zm1nrm-znrm)*dxnrm
            else
               err=znrm/bnrm
               goto 100
            endif
            ! calculate the norm of the x_vector
            xnrm=snrm(n_size,x_vector,itol)
            if (err.le.0.5_dblprec*xnrm) then
               err=err/xnrm
            else
               err=znrm/bnrm
               goto 100
            endif
         endif

         write(*,*) 'iter=',iter,'err=',err
         if(err.gt.tol) goto 100
      endif

      return

   end subroutine linbcg

   ! this function calculates the norm of a vector that will be used in linbcg
   ! this function is also extracted from the numerical recipes for fortran 90
   real(dblprec) function snrm(n_size,s_vector,itol)

      implicit none

      integer, intent(in) :: itol ! flag to see which type of calculation will be made
      integer, intent(in) :: n_size ! size of the vector

      real(dblprec), dimension(n_size), intent(in) :: s_vector ! vector which will be calculated

      integer :: isamax

      integer :: i

      ! compute the vector magnitude norm
      if (itol.le.3) then
         snrm=0
         do i=1,n_size
            snrm=snrm+s_vector(i)**2
         enddo
         snrm=sqrt(snrm)
      else
         isamax=1
         ! largest component norm
         do i=1,n_size
            if (abs(s_vector(i)).gt.abs(s_vector(isamax))) isamax=i
         enddo
         snrm=abs(s_vector(isamax))
      endif

      return

   end function snrm

   ! multiplication of a matrix by a vector a.x=b that will be used in linbcg
   ! this function is also extracted from the numerical recipes for fortran 90
   subroutine atimes(n_size,x_vector,r_solution,nmax,a_matrix,itrnsp,ija)

      implicit none

      integer, intent(in) :: nmax ! size of the matrix
      integer, intent(in) :: n_size ! size of the vector in x which multiplies the matrix
      integer, intent(in) :: itrnsp
      integer, dimension(nmax), intent(in) :: ija ! index vector array storage
      real(dblprec), dimension(nmax), intent(in) :: a_matrix
      real(dblprec), dimension(n_size), intent(in) :: x_vector

      real(dblprec), dimension(n_size), intent(out) :: r_solution

      if (itrnsp.eq.0) then
         ! multiply a_matrix.x_vector=r_solution
         call dsprsax(a_matrix,ija,x_vector,r_solution,n_size,nmax)
      else
         ! multiply (a_matrix)^t.x_vector=r_solution
         call dsprtstx(a_matrix,ija,x_vector,r_solution,n_size,nmax)
      endif

      return

   end subroutine atimes

   ! solves the system x=b.a^-1 for the diagonal part of the matrix
   ! this function is also extracted from the numerical recipes for fortran 90
   subroutine asolve(n_size,nmax,a_matrix,b_solution,x_vector)

      implicit none

      integer, intent(in) :: nmax
      integer, intent(in) :: n_size ! size of the vector in x which multiplies the matrix

      real(dblprec), dimension(nmax), intent(in) :: a_matrix ! matrix ordered in a single array
      real(dblprec), dimension(n_size), intent(in) :: b_solution ! the solution of multiplying a.x=b

      real(dblprec), dimension(n_size), intent(out) :: x_vector ! vector which multiplies the matrix

      integer :: i

      do i=1, n_size
         x_vector(i)=b_solution(i)/a_matrix(i)
      enddo

      return

   end subroutine asolve

   ! multiply a matrix in a row index sparse storage arrays sa and ija by a vector x(1:n), giving a vector b(1:n)
   ! this function is also extracted from the numerical recipes for fortran 90
   subroutine dsprsax(a_matrix,ija,x_vector,b_solution,n_size,nmax)

      implicit none

      integer, intent(in) :: nmax
      integer, intent(in) :: n_size ! size of the vector in x which multiplies the matrix
      integer, dimension(nmax), intent(in) :: ija ! index vector array storage
      real(dblprec), dimension(nmax), intent(in) :: a_matrix ! matrix ordered in a single array
      real(dblprec), dimension(n_size), intent(in) :: x_vector ! vector which multiplies the matrix

      real(dblprec), dimension(n_size), intent(out) :: b_solution ! the solution of multiplying a.x=b

      integer :: i,k

      !if (ija(1).ne.n_size+2) pause 'mismatched vector and matrix dsprsax'
      if (ija(1).ne.n_size+2) write(*,*) '*** mismatched vector and matrix dsprsax ***'

      do i=1,n_size
         ! loop over diagonal terms
         b_solution(i)=a_matrix(i)*x_vector(i)
         do k=ija(i),ija(i+1)-1
            ! loop over off-diagonal terms
            b_solution(i)=b_solution(i)+a_matrix(k)*x_vector(ija(k))
         enddo
      enddo

      return

   end subroutine dsprsax

   ! multiply the transpose matrix in a row-index sparse storage arrays sa and ija by a vector x(1:n), giving vector b(1:n)
   ! this function is also extracted from the numerical recipes for fortran 90
   subroutine dsprtstx(a_matrix,ija,x_vector,b_solution,n_size,nmax)

      implicit none

      integer, intent(in) :: nmax
      integer, intent(in) :: n_size ! size of the vector in x which multiplies the matrix
      integer, dimension(nmax), intent(in) :: ija ! index vector array storage
      real(dblprec), dimension(nmax), intent(in) :: a_matrix ! matrix ordered in a single array
      real(dblprec), dimension(n_size), intent(in) :: x_vector ! vector which multiplies the matrix

      real(dblprec), dimension(n_size), intent(out) :: b_solution ! the solution of multiplying a.x=b

      integer :: i,j,k

      !     if (ija(1).ne.n_size+2) pause 'mismatched vector and matrix in dsprstx'
      if (ija(1).ne.n_size+2) write(*,*) '*** mismatched vector and matrix in dsprstx ***'

      ! start with diagonal terms
      do i=1,n_size
         b_solution(i)=a_matrix(i)*x_vector(i)
      enddo

      ! loop over off diagonal terms
      do i=1,n_size
         do k=ija(i),ija(i+1)-1
            j=ija(k)
            b_solution(j)=b_solution(j)+a_matrix(k)*x_vector(i)
         enddo
      enddo

      return

   end subroutine dsprtstx


   subroutine dsprsin(amm,n_size,nmax,a_matrix,ija)

      implicit none

      integer, intent(in) :: n_size ! number of rows(columns) of the matrix
      real(dblprec), dimension(n_size,n_size), intent(in) :: amm ! matrix that is going to be saved in a vector arrangement

      integer, intent(out) :: nmax ! size of the a_matrix
      integer, dimension(:), allocatable, intent(out) :: ija ! saving the indexes of the matrix in a vector arrangement
      real(dblprec), dimension(:), allocatable, intent(out) :: a_matrix ! matrix ordered as a vector arrangement

      integer :: i,j,k
      integer :: i_stat

      nmax=(n_size+1)**2

      allocate(ija(nmax),stat=i_stat)
      call memocc(i_stat,product(shape(ija))*kind(ija),'ija','dsprsin')

      allocate(a_matrix(nmax),stat=i_stat)
      call memocc(i_stat,product(shape(a_matrix))*kind(a_matrix),'a_matrix','dsprsin')

      ! saving the diagonal elements of the a(1:n,1:n) matrix in the a_matrix
      do j=1, n_size
         a_matrix(j)=amm(j,j)
      enddo
      ija(1)=n_size+2
      k=n_size+1

      ! saving the off-diagonal terms of the a(1:n,1:n) matrix
      do i=1, n_size
         do j=1,n_size
            if (i.ne.j) then
               k=k+1
               a_matrix(k)=amm(i,j)
               ija(k)=j
            endif
         enddo
         ija(i+1)=k+1
      enddo

      return

   end subroutine dsprsin


   ! svd decomposition of a matrix to deal with singular systems
   subroutine svdcmp(a_matrix,v_matrix,w_vector,nmax,m_size,n_size,mp_size,np_size)
      ! given a matrix a(1:m,1:n), with physical dimensions mp by np, this routine computes
      ! its singular value decomposition, a=uwv^t. the matrix u replaces a on output.
      ! the diagonal matrix of singulat values w is output as a vector w(1:n).
      ! the matrix v (not the transpose v^t) is output as v(1:n,1:n).

      implicit none

      integer, intent(in) :: nmax
      integer, intent(in) :: m_size
      integer, intent(in) :: n_size
      integer, intent(in) :: mp_size
      integer, intent(in) :: np_size
      real(dblprec), dimension(n_size), intent(out) :: w_vector
      real(dblprec), dimension(n_size,n_size), intent(out) :: v_matrix
      real(dblprec), dimension(m_size,n_size), intent(inout) :: a_matrix

      integer :: i, its, j, jj, k, l, nm, max_its
      real(dblprec) :: anorm, c_val,f_val,g_val,h_val,s_val,scale_val,x_val,y_val,z_val
      real(dblprec), dimension(nmax) :: rv1

      g_val=0.0_dblprec
      scale_val=0.0_dblprec
      anorm=0.0_dblprec
      max_its=1e6

      do i=1, n_size
         l=i+1
         rv1(i)=scale_val*g_val
         g_val=0.0_dblprec
         s_val=0.0_dblprec
         scale_val=0.0_dblprec

         if (i.le.m_size) then
            do k=i, m_size
               scale_val=scale_val+abs(a_matrix(k,i))
            enddo
            if (scale_val.ne.0.0_dblprec) then
               do k=i,m_size
                  a_matrix(k,i)=a_matrix(k,i)/scale_val
                  s_val=s_val+a_matrix(k,i)*a_matrix(k,i)
               enddo
               f_val=a_matrix(i,i)
               g_val=-sign(sqrt(s_val),f_val)
               h_val=f_val*g_val-s_val
               a_matrix(i,i)=f_val-g_val
               do j=l,n_size
                  s_val=0.0_dblprec
                  do k=i,m_size
                     s_val=s_val+a_matrix(k,i)*a_matrix(k,j)
                  enddo
                  f_val=s_val/h_val
                  do k=i,m_size
                     a_matrix(k,j)=a_matrix(k,j)+f_val*a_matrix(k,i)
                  enddo
               enddo
               do k=i,m_size
                  a_matrix(k,i)=scale_val*a_matrix(k,i)
               enddo
            endif
         endif

         w_vector(i)=scale_val*g_val
         g_val=0.0_dblprec
         s_val=0.0_dblprec
         scale_val=0.0_dblprec

         if ((i.le.m_size).and.(i.ne.n_size)) then
            do k=l,n_size
               scale_val=scale_val+abs(a_matrix(i,k))
            enddo
            if (scale_val.ne.0.0_dblprec) then
               do k=l,n_size
                  a_matrix(i,k) = a_matrix(i,k)/scale_val
                  s_val=s_val+a_matrix(i,k)*a_matrix(i,k)
               enddo
               f_val=a_matrix(i,l)
               g_val=-sign(sqrt(s_val),f_val)
               h_val=f_val*g_val-s_val
               a_matrix(i,l)=f_val-g_val
               do k=l,n_size
                  rv1(k)=a_matrix(i,k)/h_val
               enddo
               do j=l,m_size
                  s_val=0.0_dblprec
                  do k=l,n_size
                     s_val=s_val+a_matrix(j,k)*a_matrix(i,k)
                  enddo
                  do k=l,n_size
                     a_matrix(j,k)=a_matrix(j,k)+s_val*rv1(k)
                  enddo
               enddo
               do k=1, n_size
                  a_matrix(i,k)=scale_val*a_matrix(i,k)
               enddo
            endif
         endif
         anorm=max(anorm,(abs(w_vector(i))+abs(rv1(i))))
      enddo

      do i=n_size,1,-1
         if (i.lt.n_size) then
            if (g_val.ne.0.0_dblprec)  then
               do j=l, n_size
                  v_matrix(j,i)=(a_matrix(i,j)/a_matrix(i,l))/g_val
               enddo
               do j=l,n_size
                  s_val=0.0_dblprec
                  do k=l, n_size
                     s_val=s_val+a_matrix(i,k)*v_matrix(k,j)
                  enddo
               enddo
            endif
            do j=l, n_size
               v_matrix(i,j)=0.0_dblprec
               v_matrix(j,i)=0.0_dblprec
            enddo
         endif
         v_matrix(i,i)=1.0_dblprec
         g_val=rv1(i)
         l=i
      enddo

      do i=min(m_size,n_size),1,-1
         l=i+1
         g_val=w_vector(i)
         do j=l, n_size
            a_matrix(i,j)=0.0_dblprec
         enddo
         if (g_val.ne.0.0_dblprec) then
            g_val=1.0_dblprec/g_val
            do j=l, n_size
               s_val=0.0_dblprec
               do k=l,m_size
                  s_val=s_val+a_matrix(k,i)*a_matrix(k,j)
               enddo
               f_val=(s_val/a_matrix(i,i))*g_val
               do k=i, m_size
                  a_matrix(k,j)=a_matrix(k,j)+f_val*a_matrix(k,i)
               enddo
            enddo
            do j=i, m_size
               a_matrix(j,i)=a_matrix(j,i)*g_val
            enddo
         else
            do j=i, m_size
               a_matrix(j,i)=0.0_dblprec
            enddo
         endif
         a_matrix(i,i)=a_matrix(i,i)+1.0
      enddo

      do k=n_size,1,-1
         do its=1,max_its
            do l=k,1,-1
               nm=l-1
               !write(*,*) k,l,nm, rv1(l)
               if ((abs(rv1(l))+anorm).eq.anorm) goto 2
               if ((abs(w_vector(nm))+anorm).eq.anorm) goto 1
            enddo
            1         c_val=0.0_dblprec
            s_val=0.0_dblprec
            do i=l,k
               f_val=s_val*rv1(i)
               rv1(i)=c_val*rv1(i)
               if ((abs(f_val)+anorm).eq.anorm) goto 2
               g_val=w_vector(i)
               h_val=pythag(f_val,g_val)
               w_vector(i)=h_val
               h_val=1.0_dblprec/h_val
               c_val=(g_val*h_val)
               s_val=-(f_val*h_val)
               do j=1, m_size
                  y_val=a_matrix(j,nm)
                  z_val=a_matrix(j,i)
                  a_matrix(j,nm)=(y_val*c_val)+(z_val*s_val)
                  a_matrix(j,i)=-(y_val*s_val)+(z_val*c_val)
               enddo
            enddo
            2         z_val=w_vector(k)
            if (l.eq.k) then
               if (z_val.lt.0.0_dblprec) then
                  w_vector(k)=-z_val
                  do j=1, n_size
                     v_matrix(j,k)=-v_matrix(j,k)
                  enddo
               endif
               goto 3
            endif
            !          if (its.eq.max_its) pause 'no convergence in svdcmp'
            if (its.eq.max_its) write(*,*) '*** no convergence in svdcmp ***'
            x_val=w_vector(l)
            nm=k-1
            y_val=w_vector(nm)
            g_val=rv1(nm)
            h_val=rv1(k)
            f_val=((y_val-z_val)*(y_val+z_val)+(g_val-h_val)*(g_val+h_val))/(2.0_dblprec*h_val*y_val)
            g_val=pythag(f_val,1.0_dblprec)
            f_val=((x_val-z_val)*(x_val+z_val)+h_val*((y_val/(f_val+sign(g_val,f_val)))-h_val))/x_val
            c_val=1.0_dblprec
            s_val=1.0_dblprec
            do j=l, nm
               i=j+1
               g_val=rv1(i)
               y_val=w_vector(i)
               h_val=s_val*g_val
               g_val=c_val*g_val
               z_val=pythag(f_val,h_val)
               rv1(j)=z_val
               c_val=f_val/z_val
               s_val=h_val/z_val
               f_val=(x_val*c_val)+(g_val*s_val)
               g_val=-(x_val*s_val)+(g_val*c_val)
               h_val=y_val*s_val
               y_val=y_val*c_val
               do jj=1, n_size
                  x_val=v_matrix(jj,j)
                  z_val=v_matrix(jj,i)
                  v_matrix(jj,j)=(x_val*c_val)+(z_val*s_val)
                  v_matrix(jj,i)=-(x_val*s_val)+(z_val*c_val)
               enddo
               z_val=pythag(f_val,h_val)
               w_vector(j)=z_val
               if (z_val.ne.0.0_dblprec) then
                  z_val=1.0_dblprec/z_val
                  c_val=f_val*z_val
                  s_val=h_val*z_val
               endif
               f_val=(c_val*g_val)+(s_val*s_val)
               x_val=-(s_val*g_val)+(c_val*y_val)
               do jj=1, m_size
                  y_val=a_matrix(jj,j)
                  z_val=a_matrix(jj,i)
                  a_matrix(jj,j)=(y_val*c_val)+(z_val*s_val)
                  a_matrix(jj,i)=-(y_val*s_val)+(z_val*c_val)
               enddo
            enddo
            rv1(l)=0.0_dblprec
            rv1(k)=f_val
            w_vector(k)=x_val
         enddo
         3      continue
      enddo

      return

   end subroutine svdcmp


   ! backsubstitution for the svd scheme several transformations must be done before
   subroutine svbksb(u_matrix,w_vector,v_matrix,m_size,n_size,mp_size,np_size,nmax,b_vector,x_vector)
      ! solves a.x=b for a vector x, where a is specified by the arrays u, w, v as returned by svdcmp.
      ! m and n are the logical dimensions of a, and will be equal for square matrices.
      ! mp and np are the physical dimensions of a. b(1:m) is the input right-hand side.
      ! x(1:n) is the output solution vector. no input quantities are destroyed, so the routine may be called
      ! sequentially with differen b's.

      implicit none

      integer, intent(in) :: n_size
      integer, intent(in) :: m_size
      integer, intent(in) :: mp_size
      integer, intent(in) :: np_size
      integer, intent(in) :: nmax

      real(dblprec), dimension(mp_size), intent(in) :: b_vector ! rhs of the equation to solve
      ! obtained from the svdcmp
      real(dblprec), dimension(np_size), intent(in) :: w_vector
      real(dblprec), dimension(np_size,np_size), intent(in) :: v_matrix
      real(dblprec), dimension(mp_size,np_size), intent(in) :: u_matrix !a_matrix as output from the sdvcmp

      real(dblprec), dimension(np_size), intent(out) :: x_vector

      integer :: i,j,jj
      real(dblprec) :: s_val
      real(dblprec), dimension(nmax) :: tmp

      do j=1, n_size
         s_val=0.0_dblprec
         if (w_vector(j).ne.0.0_dblprec) then
            do i=1, m_size
               s_val=s_val+u_matrix(i,j)*b_vector(i)
            enddo
            s_val=s_val/w_vector(j)
         endif
         tmp(j)=s_val
      enddo
      do j=1,n_size
         s_val=0.0_dblprec
         do jj=1, n_size
            s_val=s_val+v_matrix(j,jj)*tmp(jj)
         enddo
         x_vector(j)=s_val
      enddo

      return

   end subroutine svbksb

   ! calculate the pythagoras rule
   real(dblprec) function pythag(a_val,b_val)

      implicit none

      real(dblprec), intent(in) :: a_val
      real(dblprec), intent(in) :: b_val

      real(dblprec) :: absa,absb

      absa=abs(a_val)
      absb=abs(b_val)

      if (absa.gt.absb) then
         pythag=absa*sqrt(1.0_dblprec+(absb/absa)**2)
      else
         if (absb.eq.0.0_dblprec) then
            pythag=0.0_dblprec
         else
            pythag=absb*sqrt(1.0_dblprec+(absa/absb)**2)
         endif
      endif

      return
   end function

   ! singular value decomposition caller
   subroutine svd(a_matrix,u_matrix,s_vector,v_matrix,m_size,n_size)

      implicit none

      integer, intent(in) :: n_size ! number of columns of the matrix
      integer, intent(in) :: m_size ! number of rows of the matrix

      real(dblprec), dimension(n_size,m_size), intent(in) :: a_matrix ! input matrix that will be subjected to svd

      real(dblprec), dimension(n_size), intent(out) :: s_vector  ! singular values vector
      real(dblprec), dimension(m_size,m_size), intent(out) :: u_matrix  ! output unitary matrix
      real(dblprec), dimension(n_size,n_size), intent(out) :: v_matrix  ! output unitary matrix
      real(dblprec), dimension(n_size,n_size) :: vt_matrix

      !
      ! program computes the matrix singular value decomposition.
      ! using lapack library.
      !
      ! programmed by sukhbinder singh
      ! 14th january 2011
      !


      real(dblprec), dimension(:), allocatable :: work
      integer :: lda,m,n,lwork,ldvt,info,ldu,i,j
      character :: jobu, jobvt

      jobu='a'
      jobvt='a'
      lda=m_size
      ldu=m_size
      ldvt=n_size

      lwork=max(1,3*min(m_size,n_size)+max(m_size,n_size),5*min(m_size,n_size))

      allocate(work(lwork))

      call dgesvd(jobu, jobvt, m_size, n_size, a_matrix, lda, s_vector, u_matrix, ldu, vt_matrix, ldvt, &
         work, lwork, info )

      do i=1,2
         do j=1,2
            v_matrix(j,i)=vt_matrix(i,j)
         end do
      end do

   end subroutine svd

end module temperature_common
