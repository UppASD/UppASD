!< Module to implement different temperature profiles in the sample
!@TODO generaliza to all kind of system shapes
!> @copyright
!! Copyright (C) 2008-2018 UppASD group
!! This file is distributed under the terms of the
!! GNU General Public License. 
!! See http://www.gnu.org/copyleft/gpl.txt
module temperature
   use Parameters
   use Profiling
   USE CONSTANTS

   implicit none

   ! FDM variables
   integer :: temp_max_no_neigh
   integer :: temp_max_no_equiv
   integer, dimension(:), allocatable         :: temp_nlistsize  !< Size of neighbour list for temperature grid
   integer, dimension(:,:), allocatable       :: temp_nlist      !< List of neighbours for temperature grid
   integer, dimension(:,:), allocatable       :: temp_nmdim      !< Dimension of neighbour map for temperature grid
   integer, dimension(:,:,:), allocatable     :: temp_nm         !< Neighbour map for temeprature grid

   !!!$ ! Meshless
   !!!$ real(dblprec), dimension(:,:), allocatable :: S_Matrix
   !!!$ real(dblprec), dimension(:), allocatable   :: A_Matrix
   !!!$ integer, dimension(:), allocatable         :: indx
   !!!$ integer, dimension(:), allocatable         :: ija
   !!!$ real(dblprec), dimension(:), allocatable   :: p_vec
   !!!$ real(dblprec), dimension(:), allocatable   :: c_vec
   !!!$ real(dblprec), dimension(:), allocatable   :: W_vector
   !!!$ real(dblprec), dimension(:,:), allocatable :: pTp
   !!!$ real(dblprec), dimension(:,:), allocatable :: U_Matrix
   !!!$ real(dblprec), dimension(:,:), allocatable :: V_Matrix
   !!!$ real(dblprec), dimension(:,:), allocatable :: Mom_Matrix
   !!!$ real(dblprec), dimension(:,:), allocatable :: Div_Mom_Matrix
   !!!$ real(dblprec), dimension(:,:), allocatable :: Div2_Mom_Matrix

   real(dblprec), dimension(:), allocatable   :: Temp_array      !< Temperature as an array
   real(dblprec), dimension(:,:),allocatable  :: ipTemp_array    !< Temperature (array) for initial phase

   real(dblprec),dimension(:),allocatable     :: spinTemp_array  !< Temperature (array) for initial phase
   integer :: spintemp_samp_done

   ! Temperature gradient
   character(len=1) :: grad                 !< If temperature gradient is active (Y/N)
   character(len=20) :: tempfile            !< Temperature file
   real(dblprec) :: TEMP_MAX                !< Amplitude of the gaussian temperature profile
   integer :: temp_solver                   !< Solver used for the temperature routine
   integer :: crys_symm                     !< Symmetry in the Crystallographic cell
   integer :: dim_sys                       !< System simensionality
   integer :: no_shells_num                 !< Number of shells considered
   integer :: num_tot_neigh                 !< Number of total neighbours
   integer :: temp_max_no_shells            !< Actual maximum number of shells
   character(len=1) :: grid_type            !< Type of gird being considered
   character(len=1) :: eq_type              !< Type of equation being solved (Poisson or Laplace)
   character(len=1) :: source_type          !< Source type for the Poisson equation
   CHARACTER(LEN=35) :: LOADTEMP            !< File name for temperature "restartfile"
   character(len=10) :: I1_MAX_BORDER       !< Type of boundary conditions in the I1-Max border
   character(len=10) :: I1_MIN_BORDER       !< Type of boundary conditions in the I1-Min border
   character(len=10) :: I2_MAX_BORDER       !< Type of boundary conditions in the I2-Max border
   character(len=10) :: I2_MIN_BORDER       !< Type of boundary conditions in the I2-Min border
   character(len=10) :: I3_MAX_BORDER       !< Type of boundary conditions in the I3-Max border
   character(len=10) :: I3_MIN_BORDER       !< Type of boundary conditions in the I3-Min border
   REAL(DBLPREC) :: TEMP_I1_MIN             !< Temperature of border I1=0 for constant temp
   REAL(DBLPREC) :: TEMP_I1_MAX             !< Temperature of border I1=N1 for constant temp
   REAL(DBLPREC) :: TEMP_I2_MIN             !< Temperature of border I2=0 for constant temp
   REAL(DBLPREC) :: TEMP_I2_MAX             !< Temperature of border I2=N2 for constant temp
   REAL(DBLPREC) :: TEMP_I3_MIN             !< Temperature of border I3=0 for constant temp
   REAL(DBLPREC) :: TEMP_I3_MAX             !< Temperature of border I3=N3 for constant temp
   REAL(DBLPREC) :: TEMP_I1_MIN_LOW         !< Minimum tempearture of border I1=0 for linear profile
   REAL(DBLPREC) :: TEMP_I1_MIN_HIGH        !< Maximum tempearture of border I1=0 for linear profile
   REAL(DBLPREC) :: TEMP_I1_MAX_LOW         !< Minimum tempearture of border I1=N1 for linear profile
   REAL(DBLPREC) :: TEMP_I1_MAX_HIGH        !< Maximum tempearture of border I1=N1 for linear profile
   REAL(DBLPREC) :: TEMP_I2_MIN_LOW         !< Minimum tempearture of border I2=0 for linear profile
   REAL(DBLPREC) :: TEMP_I2_MIN_HIGH        !< Maximum tempearture of border I2=0 for linear profile
   REAL(DBLPREC) :: TEMP_I2_MAX_LOW         !< Minimum tempearture of border I2=N2 for linear profile
   REAL(DBLPREC) :: TEMP_I2_MAX_HIGH        !< Maximum tempearture of border I2=N2 for linear profile
   REAL(DBLPREC) :: TEMP_I3_MIN_LOW         !< Minimum tempearture of border I3=0 for linear profile
   REAL(DBLPREC) :: TEMP_I3_MIN_HIGH        !< Maximum tempearture of border I3=0 for linear profile
   REAL(DBLPREC) :: TEMP_I3_MAX_LOW         !< Minimum tempearture of border I3=N3 for linear profile
   REAL(DBLPREC) :: TEMP_I3_MAX_HIGH        !< Maximum tempearture of border I3=N3 for linear profile
   REAL(DBLPREC), DIMENSION(3) :: R_CENTER  !< Center of the Gaussian temperature profile
   REAL(DBLPREC), DIMENSION(3) :: SIGMATEMP !< Sigma parameter for the Gaussian temperature profile
   integer :: init_temp
   integer, dimension(:), allocatable :: temp_NN
   real(dblprec), dimension(:,:,:), allocatable :: temp_neigh_dist

   ! Legacy variables
   ! Temperature gradient
   character(len=20), dimension(6) :: bounds !< array for the boundary conditions of the temperature profile
   real(dblprec) :: Temp_low_x,Temp_high_x, Temp_high_y, Temp_low_y, Temp_high_z,Temp_low_z !< Boundary conditions for the 1D laplace equation
   integer :: barr_size ! Size of the stpe for the tprof=1Dstep


   private

   public :: Temp_array, ipTemp_array, spintemperature
   public :: grad, tempfile,Temp_low_x,Temp_high_x, Temp_high_y, Temp_low_y
   public :: allocate_temp, setup_temp, read_temperature, deallocate_temp, set_temperature_defaults
   public :: read_temperature_legacy, Lparray

contains

   !!!!!++ spintemp
   !!-------------Spin Temperature (Lars)----------------------!!
   subroutine spintemperature(Natom, Mensemble, mstep,ntemp,simid,emomM, beff,flag)

      implicit none

      integer, intent(in) :: Natom !< Number of atoms in system
      integer, intent(in) :: Mensemble !< Number of ensembles
      integer, intent(in) :: mstep !< Present Step
      integer, intent(in) :: ntemp !< Total Steps
      character(len=8), intent(in) :: simid !< Name of simulation
      real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: emomM  !< Current magnetic moment vector
      real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: beff !< Effective field
      integer,intent(in) :: flag !< Flag for init,meas and print

      integer :: l,r,i_stat
      character(len=30) :: filn
      real(dblprec) :: spintemp, spintempnom, spintempdenom

      if (flag==0) then
         allocate(spinTemp_array(ntemp),stat=i_stat)
         call memocc(i_stat,product(shape(spinTemp_array))*kind(spinTemp_array),'spinTemp_array','spintemperature')
         spintemp_array=0.d0
         spintemp_samp_done=1
      endif
      if (flag==1) then
         spintemp=0.d0
         spintempnom=0d0
         spintempdenom=0d0
         do l=1,Mensemble
            do r=1,Natom
               spintempnom = spintempnom + &
                  ((emomM(2,r,l)*beff(3,r,l)-emomM(3,r,l)*beff(2,r,l))**2 + &
                  (emomM(3,r,l)*beff(1,r,l)-emomM(1,r,l)*beff(3,r,l))**2 + &
                  (emomM(1,r,l)*beff(2,r,l)-emomM(2,r,l)*beff(1,r,l))**2)
               spintempdenom = spintempdenom + &
                  (emomM(1,r,l)*beff(1,r,l)+emomM(2,r,l)*beff(2,r,l)+ emomM(3,r,l)*beff(3,r,l))
            end do
         enddo
         spintemp = spintempnom / spintempdenom

         spintemp=spintemp/Mensemble
         spintemp=(mub*spintemp)/(2*k_bolt)

         spintemp_array(spintemp_samp_done)=spintemp

         write (filn,'(''spintemp.'',a8,''.out'')') simid
         open(ofileno,file=filn, position='append')

         write(ofileno,20002) mstep,spintemp_array(spintemp_samp_done), &
            sum(spintemp_array(1:spintemp_samp_done))/spintemp_samp_done, &
            sum(abs(spintemp_array(1:spintemp_samp_done)))/spintemp_samp_done

         close(ofileno)
         spintemp_samp_done=spintemp_samp_done+1

      endif
      if (flag==2) then
         deallocate(spinTemp_array,stat=i_stat)
         call memocc(i_stat,-product(shape(spinTemp_array))*kind(spinTemp_array),'spinTemp_array','spintemperature')
      endif

      20002 format (i8,2x,es16.8,es16.8,es16.8)

   end subroutine spintemperature




   subroutine allocate_temp (Natom, ip_nphase)

      integer, intent(in) :: Natom
      integer, intent(in) :: ip_nphase

      integer :: i_stat

      allocate(Temp_array(Natom),stat=i_stat)
      call memocc(i_stat,product(shape(Temp_array))*kind(Temp_array),'Temp_array','allocate_temp')
      if(ip_nphase>0) then
         allocate(ipTemp_array(Natom,ip_nphase),stat=i_stat)
         call memocc(i_stat,product(shape(ipTemp_array))*kind(ipTemp_array),'ipTemp_array','allocate_temp')
      end if

   end subroutine allocate_temp


   subroutine deallocate_temp ()

      integer :: i_stat,i_all

      i_all=-product(shape(Temp_array))*kind(Temp_array)
      deallocate(Temp_array,stat=i_stat)
      call memocc(i_stat,i_all,'Temp_array','deallocate_temp')
      if(allocated(ipTemp_array)) then
         i_all=-product(shape(ipTemp_array))*kind(ipTemp_array)
         deallocate(ipTemp_array,stat=i_stat)
         call memocc(i_stat,i_all,'ipTemp_array','deallocate_temp')
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

   !> Setup the temperature array of the system
   subroutine SETUP_TEMP(NATOM,NT,NA,N1,N2,N3,NATOM_FULL,&
         DO_RALLOY,ATYPE,ACELLNUMB,ATYPE_CH,SIMID,TEMP,&
         C1,C2,C3,BC1,BC2,BC3,BAS,COORD,TEMP_ARRAY)

      !
      IMPLICIT NONE

      INTEGER, INTENT(IN) :: NT  !< Number of types of atoms
      INTEGER, INTENT(IN) :: NA  !< Number of atoms in one cell
      INTEGER, INTENT(IN) :: N1  !< Number of cell repetitions in x direction
      INTEGER, INTENT(IN) :: N2  !< Number of cell repetitions in y direction
      INTEGER, INTENT(IN) :: N3  !< Number of cell repetitions in z direction
      INTEGER, INTENT(IN) :: Natom      !< Number of atoms in system
      INTEGER, INTENT(IN) :: NATOM_FULL !< Number of atoms for full system (=Natom if not dilute)
      INTEGER, INTENT(IN), OPTIONAL :: DO_RALLOY  !< Random alloy simulation (0/1)
      INTEGER, DIMENSION(NATOM), INTENT(IN) :: ATYPE !< Type of atom
      INTEGER, DIMENSION(:), OPTIONAL ,INTENT(IN) :: ACELLNUMB !< List for translating atom no. in full cell to actual cell
      INTEGER, DIMENSION(:), OPTIONAL ,INTENT(IN) :: ATYPE_CH  !< Actual type of atom for dilute system
      CHARACTER(LEN=1), INTENT(IN) :: BC1  !< Boundary conditions in x-direction
      CHARACTER(LEN=1), INTENT(IN) :: BC2  !< Boundary conditions in y-direction
      CHARACTER(LEN=1), INTENT(IN) :: BC3  !< Boundary conditions in z-direction
      CHARACTER(LEN=8), INTENT(IN) :: SIMID       !< Simulation name
      REAL(DBLPREC), INTENT(IN) :: TEMP              !< Constant temperature from input

      REAL(DBLPREC), DIMENSION(3), INTENT(IN) :: C1       !< First lattice vector
      REAL(DBLPREC), DIMENSION(3), INTENT(IN) :: C2       !< Second lattice vector
      REAL(DBLPREC), DIMENSION(3), INTENT(IN) :: C3       !< Third lattice vector
      REAL(DBLPREC), DIMENSION(3,NA), INTENT(IN) :: BAS   !< Coordinates for basis atoms
      REAL(DBLPREC), DIMENSION(3,NATOM), INTENT(IN) :: COORD !< Coordinates of the atoms in the system
      REAL(DBLPREC), DIMENSION(NATOM), INTENT(OUT) :: TEMP_ARRAY !< Temperature array

      !.. Local variables

      ! If there is no gradient the array is filled with a constant temperature
      IF (GRAD.EQ.'N') THEN
         TEMP_ARRAY=TEMP
         WRITE(*,'(2x,a)') ' Homogeneous Temperature'
         ! If the gradient is set to F one reads the temperature from a file
      ENDIF


      ! PRINTING THE FINAL TEMPERATURE FILE AND FOR COMPARISON THE INTIAL TEMPERATURE FILE

      RETURN

   END SUBROUTINE SETUP_TEMP


   ! Routine to read the temperature file for thermal gradients
   subroutine  read_temperature()

      use FileParser
      use InputData, only: maptype, C1, C2, C3, Bas, nt, posfiletype, atype_inp

      implicit none

      integer :: i_err,rd_len,i_errb
      integer :: num_tot_neigh,isite,i_stat,jsite
      integer :: itype,jtype,iline,ishell,i_all
      logical :: unique
      logical :: comment
      real(dblprec):: tol, norm,T1,T2
      real(dblprec), dimension(3) :: r_tmp, r_red
      integer :: no_shells_num
      character(len=50) :: keyword
      real(dblprec), dimension(:,:,:), allocatable :: temp_neigh_dist_tmp

      ! Set tolerance for neighbour shells
      tol=1.0d-5

      open(ifileno,file=trim(tempfile))
      do
         10     continue

         ! Read file character for character until first whitespace
         keyword=""
         call bytereader(keyword,rd_len,ifileno,i_errb)

         ! converting Capital letters
         call caps2small(keyword)

         ! check for comment markers (currently % and #)
         comment=(scan(trim(keyword),'%')==1).or.(scan(trim(keyword),'#')==1).or.&
            (scan(trim(keyword),'*')==1).or.(scan(trim(keyword),'=')==1)

         ! Parse keyword
         keyword=trim(keyword)

         select case(keyword)
         case('shells_nums')
            read(ifileno,*,iostat=i_err) no_shells_num
            if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err
         case('num_tot_neigh')
            read(ifileno,*,iostat=i_err) num_tot_neigh
            if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err
         case('temp_solver')
            read(ifileno,*,iostat=i_err) temp_solver
            if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err
         case('crys_symm')
            read(ifileno,*,iostat=i_err) crys_symm
            if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err
         case('grid_type')
            read(ifileno,*,iostat=i_err) grid_type
            if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err
         case('eq_type')
            read(ifileno,*,iostat=i_err) eq_type
            if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err
         case('init_temp')
            read(ifileno,*,iostat=i_err) init_temp
            if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err
         case('dim_sys')
            read(ifileno,*,iostat=i_err) dim_sys
            if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err
         case('x_min_border')
            read(ifileno,*,iostat=i_err) I1_MIN_BORDER
            if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err
            if (I1_MIN_BORDER.eq.'constant') then
               read(ifileno,*,iostat=i_err) TEMP_I1_MIN
            else if(I1_MIN_BORDER.eq.'linear') then
               read(ifileno,*,iostat=i_err) T1,T2
               TEMP_I1_MIN_LOW=min(T1,T2)
               TEMP_I1_MIN_HIGH=max(T1,T2)
            endif
         case('x_max_border')
            read(ifileno,*,iostat=i_err) I1_MAX_BORDER
            if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err
            if (I1_MAX_BORDER.eq.'constant') then
               read(ifileno,*,iostat=i_err) TEMP_I1_MAX
            else if(I1_MAX_BORDER.eq.'linear') then
               read(ifileno,*,iostat=i_err) T1,T2
               TEMP_I1_MAX_LOW=min(T1,T2)
               TEMP_I1_MAX_HIGH=max(T1,T2)
            endif
         case('y_min_border')
            read(ifileno,*,iostat=i_err) I2_MIN_BORDER
            if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err
            if (I2_MIN_BORDER.eq.'constant') then
               read(ifileno,*,iostat=i_err) TEMP_I2_MIN
            else if(I2_MIN_BORDER.eq.'linear') then
               read(ifileno,*,iostat=i_err) T1,T2
               TEMP_I2_MIN_LOW=min(T1,T2)
               TEMP_I2_MIN_HIGH=max(T1,T2)
            endif
         case('y_max_border')
            read(ifileno,*,iostat=i_err) I2_MAX_BORDER
            if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err
            if (I2_MAX_BORDER.eq.'constant') then
               read(ifileno,*,iostat=i_err) TEMP_I2_MAX
            else if(I2_MAX_BORDER.eq.'linear') then
               read(ifileno,*,iostat=i_err) T1,T2
               TEMP_I2_MAX_LOW=min(T1,T2)
               TEMP_I2_MAX_HIGH=max(T1,T2)
            endif
         case('z_min_border')
            read(ifileno,*,iostat=i_err) I3_MIN_BORDER
            if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err
            if (I3_MIN_BORDER.eq.'constant') then
               read(ifileno,*,iostat=i_err) TEMP_I3_MIN
            else if(I3_MIN_BORDER.eq.'linear') then
               read(ifileno,*,iostat=i_err) T1,T2
               TEMP_I3_MIN_LOW=min(T1,T2)
               TEMP_I3_MIN_HIGH=max(T1,T2)
            endif
         case('z_max_border')
            read(ifileno,*,iostat=i_err) I3_MAX_BORDER
            if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err
            if (I3_MAX_BORDER.eq.'constant') then
               read(ifileno,*,iostat=i_err) TEMP_I3_MAX
            else if(I3_MAX_BORDER.eq.'linear') then
               read(ifileno,*,iostat=i_err) T1,T2
               TEMP_I3_MAX_LOW=min(T1,T2)
               TEMP_I3_MAX_HIGH=max(T1,T2)
            endif
         case('source_type')
            read(ifileno,*,iostat=i_err) source_type
            if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err
            if (source_type.eq.'G') then
               read(ifileno,*,iostat=i_err) TEMP_MAX
               read(ifileno,*,iostat=i_err) sigmatemp(1), sigmatemp(2), sigmatemp(3)
               read(ifileno,*,iostat=i_err) R_CENTER(1), R_CENTER(2), R_CENTER(3)
            endif
         case('load_temp')
            read(ifileno,*,iostat=i_err) LOADTEMP
            if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err
         case('temp_neigh')
            allocate(temp_neigh_dist_tmp(NT,no_shells_num,3),stat=i_stat)
            call memocc(i_stat,product(shape(temp_neigh_dist_tmp))*kind(temp_neigh_dist_tmp),'temp_neigh_dist_tmp','read_temperature')
            allocate(temp_NN(NT),stat=i_stat)
            call memocc(i_stat,product(shape(temp_NN))*kind(temp_NN),'temp_NN','read_temperature')
            temp_NN(1:NT)=0
            do iline=1, num_tot_neigh
               ! Read indices and coordinates
               read (ifileno,*) isite, jsite, r_tmp
               ! Find type of site
               itype=atype_inp(isite)
               jtype=atype_inp(jsite)
               if(maptype==2) then
                  ! Calculate proper neighbour vector (from "bgfm")
                  r_red(1)=Bas(1,jsite)-Bas(1,isite)+C1(1)*r_tmp(1)+C2(1)*r_tmp(2)+C3(1)*r_tmp(3)
                  r_red(2)=Bas(2,jsite)-Bas(2,isite)+C1(2)*r_tmp(1)+C2(2)*r_tmp(2)+C3(2)*r_tmp(3)
                  r_red(3)=Bas(3,jsite)-Bas(3,isite)+C1(3)*r_tmp(1)+C2(3)*r_tmp(2)+C3(3)*r_tmp(3)
               else
                  ! Calculates neighbour vectors from direct coordinates or Cartesian
                  ! coordinates, corresponding to how the atomic positions are entered
                  if (posfiletype=='C') then
                     r_red=r_tmp
                  elseif (posfiletype=='D') then
                     r_red(1)=r_tmp(1)*C1(1)+r_tmp(2)*C2(1)+r_tmp(3)*C3(1)
                     r_red(2)=r_tmp(1)*C1(2)+r_tmp(2)*C2(2)+r_tmp(3)*C3(2)
                     r_red(3)=r_tmp(1)*C1(3)+r_tmp(2)*C2(3)+r_tmp(3)*C3(3)
                  else
                     stop 'Only posfiletype= C or D is currently supported'
                  endif
               end if
               ! Loop through earlier vectors to find equivalent shells
               unique=.true.
               do ishell=1,temp_NN(itype)
                  norm=(r_red(1)-temp_neigh_dist_tmp(itype,ishell,1))**2+ &
                     (r_red(2)-temp_neigh_dist_tmp(itype,ishell,2))**2+ &
                     (r_red(3)-temp_neigh_dist_tmp(itype,ishell,3))**2
                  if(norm<tol) then
                     unique=.false.
                  end if
               end do
               if (unique) then
                  temp_NN(itype)=temp_NN(itype)+1
                  temp_neigh_dist_tmp(itype,temp_NN(itype),1:3)=r_red(1:3)
               endif
            enddo

            temp_max_no_shells=maxval(temp_NN)

            allocate(temp_neigh_dist(NT,temp_max_no_shells,3),stat=i_stat)
            call memocc(i_stat,product(shape(temp_neigh_dist))*kind(temp_neigh_dist),'temp_neigh_dist','read_temperature')
            temp_neigh_dist=0
            !
            do ishell=1,temp_max_no_shells
               do itype=1,NT
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
               print *,"Keyword '",trim(keyword),"' is not recognized"
            end if

         end select

         ! End of file
         if (i_errb==20) goto 20

         ! End of row
         if (i_errb==10) goto 10

      end do
      goto 30

      20  continue

      30  continue

      return

      close(ofileno)

   end subroutine read_temperature

   subroutine read_temperature_legacy()
      use FileParser

      implicit none

      integer :: i_err,i_stat,rd_len,i_errb
      character(len=50) :: keyword, cache
      logical :: comment

      open(500,file=trim(tempfile))
      do
         10   continue
         ! Read file character for character until first whitespace
         keyword=""
         call bytereader(keyword,rd_len,500,i_errb)
         ! converting Capital letters
         call caps2small(keyword)
         ! check for comment markers (currently % and #)
         comment=(scan(trim(keyword),'%')==1).or.(scan(trim(keyword),'#')==1).or.&
            (scan(trim(keyword),'*')==1).or.(scan(trim(keyword),'=')==1)
         !print *,'-->','|',trim(keyword),'|'
         ! Parse keyword
         keyword=trim(keyword)
         !print *,keyword
         select case(keyword)

         case('tbound_x')

            read(500,*,iostat=i_err) bounds(1),bounds(2)
            read(500,*,iostat=i_err) Temp_high_x, Temp_low_x
            if((bounds(1).eq.'step').or.(bounds(2).eq.'step')) then
               read(500,*,iostat=i_err) barr_size
            else if ((bounds(1).eq.'step').and.(bounds(2).eq.'step')) then
               read(500,*,iostat=i_err) barr_size
            end if
            if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err

         case('tbound_y')

            read(500,*,iostat=i_err) bounds(3),bounds(4)
            if((bounds(3).ne.'N').and.(bounds(4).ne.'N')) then
               read(500,*,iostat=i_err) Temp_high_y, Temp_low_y
            end if
            if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err

         case('tbound_z')

            read(500,*,iostat=i_err) bounds(5),bounds(6)
            if((bounds(5).ne.'N').and.(bounds(6).ne.'N')) then
               read(500,*,iostat=i_err) Temp_high_z, Temp_low_z
            end if
            if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err

         case default
            if(comment.or.len(trim(keyword))==0) then
               !        print *,"Comment"
               !read(ifile,*,iostat=i_err) cache
            else
               print *,"Keyword '",trim(keyword),"' is not recognized"
            end if
            !       read(ifile,*)

         end select
         ! End of file
         if (i_errb==20) goto 20
         ! End of row
         if (i_errb==10) goto 10
      end do
      goto 30
      20  continue
      !   print *,'End of file ', keyword
      30  continue
      return

      close(500)

   end subroutine read_temperature_legacy

   subroutine set_temperature_defaults()
      !
      implicit none
      !
      !Temperature gradient
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
      I1_MIN_BORDER    = 'N'
      I1_MAX_BORDER    = 'N'
      I2_MIN_BORDER    = 'N'
      I2_MAX_BORDER    = 'N'
      I3_MIN_BORDER    = 'N'
      I3_MAX_BORDER    = 'N'
      TEMP_MAX         = 0.0D0
      TEMP_I1_MAX      = 0.0D0
      TEMP_I2_MAX      = 0.0D0
      TEMP_I3_MAX      = 0.0D0
      TEMP_I1_MIN      = 0.0D0
      TEMP_I2_MIN      = 0.0D0
      TEMP_I2_MIN      = 0.0D0
      TEMP_I1_MAX_LOW  = 0.0D0
      TEMP_I2_MAX_LOW  = 0.0D0
      TEMP_I3_MAX_LOW  = 0.0D0
      TEMP_I1_MIN_LOW  = 0.0D0
      TEMP_I2_MIN_LOW  = 0.0D0
      TEMP_I3_MIN_LOW  = 0.0D0
      TEMP_I1_MAX_HIGH = 0.0D0
      TEMP_I2_MAX_HIGH = 0.0D0
      TEMP_I3_MAX_HIGH = 0.0D0
      TEMP_I1_MIN_HIGH = 0.0D0
      TEMP_I2_MIN_HIGH = 0.0D0
      TEMP_I3_MIN_HIGH = 0.0D0
      SIGMATEMP        = 0.0D0
      R_CENTER         = 0.0D0

   end subroutine set_temperature_defaults

   subroutine Lparray(Tarray,Natom,coord,Temp,simid,printme)
      !< Subroutine to create a site dependant temperature that follows the laplace equation for a time independet configuration

      real(dblprec), intent(in) :: Temp  !< Temperature
      integer, intent(in) :: Natom !< Array size
      real(dblprec), dimension(3,Natom), intent(in) :: coord !< Coordinates of the atoms in the crystal
      real(dblprec), dimension(Natom), intent(out) :: Tarray !< Temperature array
      character(len=8), intent(in) :: simid  !< Name of simulation
      logical, intent(in) :: printme  !< Flag to determine printing of temperature file

      integer :: i,n
      character(len=30) :: filn

      n=1

      print '(4f10.4)',Temp_low_y,Temp_high_y,Temp_low_x,Temp_high_x
      do i=1, Natom
         ! Fill up the array with the results of the chosen function given by the boundary conditions

         Tarray(i) = lpfunction(Temp_low_y,Temp_high_y,Temp_low_x,Temp_high_x,coord,Natom,i,n,grad,bounds,Temp)

      end do

      ! PRINTING THE FINAL TEMPERATURE FILE AND FOR COMPARISON THE INTIAL TEMPERATURE FILE
      IF (printme) then
         WRITE(filn,'(''temperature.'',a8,''.out'')') simid
         OPEN(45,FILE=filn)
         do I=1, NATOM
            WRITE(45,'(i8,3f16.8,f16.8)') I, COORD(1:3,I),Tarray(I)
         ENDDO
         CLOSE(45)
      ENDIF


   end subroutine Lparray

   function linear_1D(Natom,Temp_high_x,Temp_low_x,coord,i) !< This creates a linear function to fill up a 1D lattice

      real(dblprec), intent(in) :: Temp_high_x, Temp_low_x ! Boundaries
      real(dblprec), dimension(3,Natom), intent(in) :: coord ! Coordinates of the  atoms in the lattice
      real(dblprec) :: linear_1D
      integer, intent(in) :: Natom
      integer :: i, x_size

      x_size = maxval(coord(1,:)) ! This works only if the person puts the direction along x maybe put in a logical variable to check

      linear_1D = Temp_low_x + (Temp_high_x-Temp_low_x)*coord(1,i)/x_size


   end function linear_1D

   ! Function to include a step function for a 1D array
   function step_x(barr_size,Natom,Temp_high_x,Temp_low_x,coord,i)

      real(dblprec), intent(in) :: Temp_high_x, Temp_low_x ! Temperatures for the step
      real(dblprec), dimension(3,Natom), intent(in) :: coord ! coordinates of the atoms in the lattice
      real(dblprec) :: step_x
      integer, intent(in) :: Natom, barr_size ! Number of atoms and size of the hot part of the lattice
      integer :: i

      if(int(coord(1,i)).le.barr_size) then

         step_x = Temp_high_x

      else

         step_x = Temp_low_x

      end if

   end function step_x

   ! this deals with the boundary condition f(x,ysize)=Tymax
   function cts_2D_x_max(Ty_max,n,coord,Natom,i)

      integer, intent(in) :: n, i, Natom
      real(dblprec), intent(in) :: Ty_max ! Boundary condition
      real(dblprec), dimension(3,Natom),intent(in) :: coord ! Coordinates of the atoms in the lattice
      real(dblprec) :: cts_2D_x_max
      real(dblprec) x_size, pi, arg_x, y_size, x , y

      parameter(pi=4.D0*DATAN(1.D0)) ! The best way to define pi

      y_size = maxval(coord(2,:)) ! To define the size of the integration region in the y direction
      x_size = maxval(coord(1,:)) ! To define the size of the integration region in the x direction

      x = coord(1,i)
      y = coord(2,i)

      arg_x=n*pi/y_size

      cts_2D_x_max = 2*Ty_max*(1-(-1)**n)*sinh(arg_x*x)*sin(arg_x*y)
      cts_2D_x_max = cts_2D_x_max/(n*pi*sinh(arg_x*x_size))

   end function cts_2D_x_max

   ! this deals with the boundary condition f(x,ysize)=Tymin
   function cts_2D_x_min(Ty_min,n,coord,Natom,i)

      integer, intent(in) :: n, i, Natom
      real(dblprec), intent(in) :: Ty_min ! Boundary condition
      real(dblprec), dimension(3,Natom),intent(in) :: coord ! Coordinates for the atoms in the lattice
      real(dblprec) :: cts_2D_x_min
      real(dblprec) x_size, pi, arg_x, y_size, x, y

      parameter(pi=4.D0*DATAN(1.D0)) ! The best way to define pi

      y_size = maxval(coord(2,:)) ! Size of the lattice in the y direction
      x_size = maxval(coord(1,:)) ! Size of the lattice in the x direction

      x = coord(1,i)
      y = coord(2,i)

      arg_x=n*pi/y_size

      cts_2D_x_min = -2*Ty_min*(1-(-1)**n)*sinh(arg_x*(x-x_size))*sin(arg_x*y)
      cts_2D_x_min = cts_2D_x_min/(n*pi*sinh(arg_x*x_size))

   end function cts_2D_x_min

   ! this deals with the boundary condition f(x,ysize)=Txmax
   function cts_2D_y_max(Tx_max,n,coord,Natom,i)

      integer, intent(in) :: n, i, Natom
      real(dblprec), intent(in) :: Tx_max ! Boundary condition
      real(dblprec), dimension(3,Natom),intent(in) :: coord ! Coordinates of the atoms in the lattice
      real(dblprec) :: cts_2D_y_max
      real(dblprec) x_size, pi, arg_y, y_size, x ,y

      parameter(pi=4.D0*DATAN(1.D0)) ! The best way to define pi

      y_size = maxval(coord(2,:)) ! Size of the lattice in the y direction
      x_size = maxval(coord(1,:)) ! Size of the lattice in the x direction

      x = coord(1,i)
      y = coord(2,i)

      arg_y=n*pi/x_size

      cts_2D_y_max = 2*Tx_max*(1-(-1)**n)*sinh(arg_y*y)*sin(arg_y*x)
      cts_2D_y_max = cts_2D_y_max/(n*pi*sinh(arg_y*y_size))

   end function cts_2D_y_max

   !! this deals with the boundary condition f(x,ysize)=Txmin
   function cts_2D_y_min(Tx_min,n,coord,Natom,i)

      integer, intent(in) :: n, i, Natom
      real(dblprec), intent(in) :: Tx_min ! Boundary condition
      real(dblprec), dimension(3,Natom),intent(in) :: coord ! Coordinates for the atoms in the lattice
      real(dblprec) :: cts_2D_y_min
      real(dblprec) x_size, pi, arg_y, y_size,x ,y

      parameter(pi=4.D0*DATAN(1.D0)) ! The best way to define pi

      y_size = maxval(coord(2,:)) ! Size of the lattice in the y direction
      x_size = maxval(coord(1,:)) ! Size of the lattice in the x direction

      x = coord(1,i)
      y = coord(2,i)

      arg_y=n*pi/x_size

      ! f_y_min, deals with f(x,0)=Tx_min and all the other boundary conditions set to zero
      cts_2D_y_min = -2*Tx_min*(1-(-1)**n)*sinh(arg_y*(y-y_size))*sin(arg_y*x)
      cts_2D_y_min = cts_2D_y_min/(n*pi*sinh(arg_y*y_size))

   end function cts_2D_y_min

   ! This deals with the boundary condition f(x,ymax)= x(Txmax-Txmin)/xsize + Txmin
   function linear_2D_y_max(Tx_min,Tx_max,Ty_min,Ty_max,n,coord,Natom,i)

      integer, intent(in) :: n, i, Natom
      real(dblprec), intent(in) :: Tx_min, Tx_max, Ty_min, Ty_max ! Boundary conditions
      real(dblprec), dimension(3,Natom), intent(in) :: coord ! Coordinates for the atoms in the lattice
      real(dblprec) :: linear_2D_y_max
      real(dblprec) :: pi, arg_y, x_size, y_size, x, y

      parameter(pi=4.D0*DATAN(1.D0)) ! The best way to define pi

      y_size = maxval(coord(2,:)) ! Size of the lattice in the y direction
      x_size = maxval(coord(1,:)) ! Size of the lattice in the x direction

      x = coord(1,i)
      y = coord(2,i)

      arg_y=n*pi/x_size

      linear_2D_y_max = 2*(Tx_min*(1-2*(-1)**n)-Tx_max*(-1)**n)
      linear_2D_y_max = linear_2D_y_max*sinh(arg_y*y)*sin(arg_y*x)
      linear_2D_y_max = linear_2D_y_max/(n*pi*sinh(arg_y*y_size))

   end function linear_2D_y_max

   ! This deals with the boundary condition f(x,0)= x(Txmax-Txmin)/xsize + Txmin
   function linear_2D_y_min(Tx_min,Tx_max,Ty_min,Ty_max,n,coord,Natom,i)

      integer, intent(in) :: n, i, Natom
      real(dblprec), intent(in) :: Tx_min, Tx_max, Ty_min, Ty_max ! Boundary conditions
      real(dblprec), dimension(3,Natom), intent(in) :: coord ! Coordinates for the atoms in the lattice
      real(dblprec) :: linear_2D_y_min
      real(dblprec) :: pi, arg_y, x_size, y_size, x, y

      parameter(pi=4.D0*DATAN(1.D0)) ! The best way to define pi

      y_size = maxval(coord(2,:)) ! Size of the lattice in the y direction
      x_size = maxval(coord(1,:)) ! Size of the lattice in the x direction

      x = coord(1,i)
      y = coord(2,i)

      arg_y=n*pi/x_size

      linear_2D_y_min = -2*(Tx_min*(1-2*(-1)**n)-Tx_max*(-1)**n)
      linear_2D_y_min = linear_2D_y_min*sinh(arg_y*(y-y_size))*sin(arg_y*x)
      linear_2D_y_min = linear_2D_y_min/(n*pi*sinh(arg_y*y_size))

   end function linear_2D_y_min

   ! This deals with the boundary condition f(0,y)= y(Tymax-Tymin)/ysize + Tymin
   function linear_2D_x_min(Tx_min,Tx_max,Ty_min,Ty_max,n,coord,Natom,i)

      integer, intent(in) :: n, i, Natom
      real(dblprec), intent(in) :: Tx_min, Tx_max, Ty_min, Ty_max ! Boundary conditions
      real(dblprec), dimension(3,Natom), intent(in) :: coord ! Coordinates for the atoms in the lattice
      real(dblprec) :: linear_2D_x_min
      real(dblprec) :: pi, arg_x, x_size, y_size, x, y

      parameter(pi=4.D0*DATAN(1.D0)) ! The best way to define pi

      y_size = maxval(coord(2,:)) ! Size of the lattice in the y direction
      x_size = maxval(coord(1,:)) ! Size of the lattice in the x direction

      x = coord(1,i)
      y = coord(2,i)

      arg_x=n*pi/y_size

      linear_2D_x_min = -2*(Ty_min*(1-2*(-1)**n)-Ty_max*(-1)**n)
      linear_2D_x_min = linear_2D_x_min*sinh(arg_x*(x-x_size))*sin(arg_x*y)
      linear_2D_x_min = linear_2D_x_min/(n*pi*sinh(arg_x*x_size))

   end function linear_2D_x_min

   ! This deals with the boundary condition f(x_max,y)= y(Tymax-Tymin)/ysize + Tymin
   function linear_2D_x_max(Tx_min,Tx_max,Ty_min,Ty_max,n,coord,Natom,i)

      integer, intent(in) :: n, i, Natom
      real(dblprec), intent(in) :: Tx_min, Tx_max, Ty_min, Ty_max ! Boundary conditions
      real(dblprec), dimension(3,Natom), intent(in) :: coord ! Coordinates for the atoms in the lattice
      real(dblprec) :: linear_2D_x_max
      real(dblprec) :: pi, arg_x, x_size, y_size, x, y

      parameter(pi=4.D0*DATAN(1.D0)) ! The best way to define pi

      y_size = maxval(coord(2,:)) ! Size of the lattice in the y direction
      x_size = maxval(coord(1,:)) ! Size of the lattice in the x direction

      x = coord(1,i)
      y = coord(2,i)

      arg_x=n*pi/y_size

      linear_2D_x_max = 2*(Ty_min*(1-2*(-1)**n)-Ty_max*(-1)**n)
      linear_2D_x_max = linear_2D_x_max*sinh(arg_x*x)*sin(arg_x*y)
      linear_2D_x_max = linear_2D_x_max/(n*pi*sinh(arg_x*x_size))

   end function linear_2D_x_max

   ! Subroutine to choose the function which will fill up the lattice
   function lpfunction(Ty_min,Ty_max,Tx_min,Tx_max,coord,Natom,i,n,grad,bounds,Temp)

      integer, intent(in) :: n, i, Natom
      integer ::  dimn, j
      real(dblprec), intent(in) :: Tx_min, Tx_max, Ty_min, Ty_max, Temp ! Boundary conditions
      real(dblprec), dimension(3,Natom), intent(in) :: coord ! Coordinates for the atoms in the lattice
      real(dblprec) :: f_x_max, f_x_min, f_y_min, f_y_max, f_z_min, f_z_max, lpfunction ! Functions relating to each of the 6 boundary conditions
      character(len=1) :: grad ! Varibale for the gradient of the system
      character(len=20), dimension(6), intent(in) :: bounds ! Array that states the boundary conditions of the sample

      ! If there is no gradient the function is just a constant
      if(grad.eq.'N') then

         lpfunction = Temp

         ! If there is a gradient the value of the function will be given by the values stored at the bounds array
      else

         do j=1,6

            ! Counting how many boundaries are 'on' N states that the boudary is off therefore indicating the dimension of the system
            if(bounds(j).ne.'N') dimn=dimn+1

         end do

         ! If there are two or less values of boundaries different from N 1D solutions are chosen
         if(dimn.le.2) then

            if(bounds(1).eq.'constant') then
               f_x_min = linear_1D(Natom,Temp_high_x,Temp_low_x,coord,i)

            else if(bounds(2).eq.'constant') then
               f_x_min = linear_1D(Natom,Temp_high_x,Temp_low_x,coord,i)

            else if((bounds(1).eq.'constant').and.(bounds(2).eq.'constant')) then
               f_x_min = linear_1D(Natom,Temp_high_x,Temp_low_x,coord,i)

            else if(bounds(1).eq.'step') then
               f_x_min = step_x(barr_size,Natom,Temp_high_x,Temp_low_x,coord,i)

            else if(bounds(2).eq.'step') then
               f_x_min = step_x(barr_size,Natom,Temp_high_x,Temp_low_x,coord,i)

            else if((bounds(1).eq.'step').and.(bounds(2).eq.'step')) then
               f_x_min = step_x(barr_size,Natom,Temp_high_x,Temp_low_x,coord,i)

               f_x_max = 0.D0
               f_y_min = 0.D0
               f_y_max = 0.D0
               f_z_max = 0.D0
               f_z_min = 0.D0

            end if

         end if

         ! If there are only two values of bounds which values equal to N then the temperature array is for 2D systems
         if(dimn.eq.4) then

            select case(bounds(1))

            case('constant')
               f_x_min = cts_2D_x_min(Ty_min,n,coord,Natom,i)
            case('linear')
               f_x_min = linear_2D_x_min(Tx_min,Tx_max,Ty_min,Ty_max,n,coord,Natom,i)

            end select

            select case(bounds(2))

            case('constant')
               f_x_max =  cts_2D_x_max(Ty_max,n,coord,Natom,i)
            case('linear')
               f_x_max = linear_2D_x_max(Tx_min,Tx_max,Ty_min,Ty_max,n,coord,Natom,i)

            end select

            select case(bounds(3))

            case('constant')
               f_y_min = cts_2D_y_min(Tx_min,n,coord,Natom,i)
            case('linear')
               f_y_min = linear_2D_y_min(Tx_min,Tx_max,Ty_min,Ty_max,n,coord,Natom,i)

            end select

            select case(bounds(4))

            case('constant')
               f_y_max = cts_2D_y_max(Tx_max,n,coord,Natom,i)
            case('linear')
               f_y_max = linear_2D_y_max(Tx_min,Tx_max,Ty_min,Ty_max,n,coord,Natom,i)

            end select

            f_z_max = 0.D0
            f_z_min = 0.D0

         end if

         ! If all the boundaries are different from N the temperature array will be filled with solutions of the 3D laplace equation
         !	if(dimn.eq.6)
         !
         !		select case(bounds(1))
         !
         !		case('constant')
         !			f_x_min = cts_3D_x_min(Ty_min,n,coord,Natom,i)
         !		case('linear')
         !			f_x_min = linear_3D_x_min(Tx_min,Tx_max,Ty_min,Ty_max,n,coord,Natom,i)
         !
         !		end select case
         !
         !		select case(bounds(2))
         !
         !		case('constant')
         !			f_x_max =  cts_3D_x_max(Ty_max,n,coord,Natom,i)
         !		case('linear')
         !			f_x_max = linear_3D_x_max(Tx_min,Tx_max,Ty_min,Ty_max,n,coord,Natom,i)
         !
         !		end select case
         !
         !		select case(bounds(3))
         !
         !		case('constant')
         !			f_y_min = cts_3D_y_min(Tx_min,n,coord,Natom,i)
         !		case('linear')
         !			f_y_min = linear_3D_y_min(Tx_min,Tx_max,Ty_min,Ty_max,n,coord,Natom,i)
         !
         !		end select case
         !
         !		select case(bounds(4))
         !
         !		case('constant')
         !			f_y_max = cts_3D_y_max(Tx_max,n,coord,Natom,i)
         !		case('linear')
         !			f_y_max = linear_3D_y_max(Tx_min,Tx_max,Ty_min,Ty_max,n,coord,Natom,i)
         !
         !		end select case
         !
         !		select case(bounds(5))
         !
         !		case('constant')
         !			f_z_min = cts_3D_z_min(Tx_max,n,coord,Natom,i)
         !		case('linear')
         !			f_z_min = linear_3D_z_min(Tx_min,Tx_max,Ty_min,Ty_max,n,coord,Natom,i)
         !
         !		end select case
         !
         !		select case(bounds(6))
         !
         !		case('constant')
         !			f_z_max = cts_3D_z_max(Tx_max,n,coord,Natom,i)
         !		case('linear')
         !			f_z_max = linear_3D_z_max(Tx_min,Tx_max,Ty_min,Ty_max,n,coord,Natom,i)
         !
         !		end select case
         !
         !	end if
         lpfunction = f_x_max + f_x_min + f_y_max + f_y_min + f_z_min + f_z_max

      end if

   end function lpfunction

end module temperature
