!------------------------------------------------------------------------------------
! MODULE: prn_trajectories 
!> Data and routines necessary for printing the trajectories of inidividual atomic moments
!> @copyright
!> GNU Public License
!------------------------------------------------------------------------------------
module prn_trajectories

   use Parameters
   use Profiling

   implicit none

   ! Input parameters
   integer :: ntraj         !< Number of trajectories to sample
   integer :: tottraj_step  !< Interval for sampling magnetic moments
   integer :: tottraj_buff  !< Buffer size for magnetic moments
   character(len=1) :: do_tottraj !< Measure magnetic moments
   integer, dimension(:), allocatable :: traj_step !< Interval for sampling individual trajectories
   integer, dimension(:), allocatable :: traj_buff !< Buffer size for individual trajectories
   integer, dimension(:), allocatable :: traj_atom !< List of atoms to sample trajectories for

   ! Measurement variables
   integer :: bcount_tottraj !< Counter of buffer for moments
   integer :: scount_tottraj !< Counter of sampling for moments
   integer, dimension(:), allocatable :: bcount_traj   !< Counter of buffer for trajectories
   real(dblprec), dimension(:), allocatable :: indxb            !< Step counter for magnetic moments
   real(dblprec), dimension(:), allocatable :: scount_traj      !< Counter of sampling for trajectories
   real(dblprec), dimension(:,:), allocatable :: indxb_traj     !< Step counter for individual trajectories
   real(dblprec), dimension(:,:,:), allocatable :: mmomb        !< Buffer for all moment magnitudes
   real(dblprec), dimension(:,:,:), allocatable :: mmomb_traj   !< Buffer for selected moment magnitudes
   real(dblprec), dimension(:,:,:,:), allocatable :: emomb      !< Buffer for all individual trajectories
   real(dblprec), dimension(:,:,:,:), allocatable :: emomb_traj !< Buffer for selected individual trajectories

   private

   public :: print_trajectories, flush_trajectories, allocate_trajectories, traj_init
   public :: read_parameters_trajectories
   public :: do_tottraj, tottraj_step, tottraj_buff, ntraj, traj_atom, traj_step, traj_buff
   public :: mmomb, mmomb_traj, emomb, emomb_traj

contains

   !---------------------------------------------------------------------------------
   ! SUBROUTINE: print_trajectories
   !> @brief Wrapper routine to print all the spin trajectories 
   !---------------------------------------------------------------------------------
   subroutine print_trajectories(Natom,sstep,mstep,Mensemble,emom,mmom,delta_t,     &
      real_time_measure,simid,do_mom_legacy,mode)

      implicit none

      integer, intent(in) :: Natom         !< Number of atoms in the system
      integer, intent(in) :: sstep         ! Simulation step in logarithmic scale
      integer, intent(in) :: mstep         !< Current simulation step
      integer, intent(in) :: Mensemble     !< Number of ensembles
      character(len=1), intent(in) :: mode   !< Simulation mode (S=SD, M=MC, H=MC Heat Bath, P=LD, C=SLD, G=GNEB)
      character(len=8), intent(in) :: simid             !< Simulation ID
      character(len=1), intent(in) :: do_mom_legacy      !< Flag to print/read moments in legacy output
      character(len=1), intent(in) :: real_time_measure !< Measurements displayed in real time   
      real(dblprec), intent(in) :: delta_t !< Time step for real time measurement
      real(dblprec), dimension(Natom, Mensemble), intent(in) :: mmom     !< Magnitude of magnetic moments
      real(dblprec), dimension(3,Natom, Mensemble), intent(in) :: emom   !< Current unit moment vector

      !.. Local variables
      integer :: i

      !  Total trajectory
      if (do_tottraj=='Y') then

         if (mod(sstep-1,tottraj_step)==0) then
            ! Write step to buffer
            call buffer_tottraj(Natom,Mensemble,mstep-1,emom,mmom,bcount_tottraj,   &
               delta_t,real_time_measure)
            if (bcount_tottraj==tottraj_buff) then
               !Write buffer to file
               call prn_tottraj(Natom, Mensemble, simid,real_time_measure,mmom,     &
                  do_mom_legacy,mode)
               bcount_tottraj=1
            else
               bcount_tottraj=bcount_tottraj+1
            endif
         else

         endif
      endif

      !  Atom trajectory
      if(ntraj>0) then
         do i=1,ntraj
            if (mod(sstep-1,traj_step(i))==0) then
               ! Write step to buffer
               call buffer_traj(Natom,Mensemble,mstep-1,emom,mmom,i,ntraj,traj_atom,&
                  bcount_traj(i),delta_t,real_time_measure)
               scount_traj(i)=1
               if (bcount_traj(i)==traj_buff(i)) then
                  call prn_traj(Mensemble,simid,ntraj,traj_atom,i,real_time_measure)
                  bcount_traj(i)=1
               else
                  bcount_traj(i)=bcount_traj(i)+1
               endif
            else
               scount_traj(i)=scount_traj(i)+1
            end if
         enddo
      endif

   end subroutine print_trajectories

   !---------------------------------------------------------------------------------
   ! SUBROUTINE: traj_init
   !> Initialization of variables with default variables for the trajectories
   !---------------------------------------------------------------------------------
   subroutine traj_init()

      implicit none

      ntraj        = 0
      do_tottraj   = 'N'
      tottraj_step = 1000
      tottraj_buff = 10

   end subroutine traj_init

   !---------------------------------------------------------------------------------
   ! SUBROUTINE: flush_trajectories
   !> Flush the trajectory measurements, i.e. print to file in the last iteration
   !---------------------------------------------------------------------------------
   subroutine flush_trajectories(Natom,Mensemble,simid,real_time_measure,mmom,      &
      do_mom_legacy,mode)

      implicit none

      integer, intent(in) :: Natom     !< Number of atoms in the system
      integer, intent(in) :: Mensemble !< Number of ensembles
      character(len=1), intent(in) :: mode   !< Simulation mode (S=SD, M=MC, H=MC Heat Bath, P=LD, C=SLD, G=GNEB)
      character(len=8), intent(in) :: simid             !< Simulation name ID
      character(len=1), intent(in) :: do_mom_legacy      !< Flag to print/read moments in legacy output
      character(len=1), intent(in) :: real_time_measure !< Perfomr measurements in real time
      real(dblprec), dimension(Natom, Mensemble), intent(in) :: mmom     !< Magnitude of magnetic moments

      !.. Local variables
      integer :: i

      ! All trajectories are printed to file in the last iteration
      if (do_tottraj=='Y') then
         !Write buffer to file  
         bcount_tottraj=bcount_tottraj-1
         call prn_tottraj(Natom, Mensemble, simid,real_time_measure,mmom,           &
            do_mom_legacy,mode)
      endif

      !  Selected trajectories are printed to file in the last iteration
      if(ntraj>0) then
         ! Write buffer to file
         do i=1,ntraj
            bcount_traj(i)=bcount_traj(i)-1
            call prn_traj(Mensemble, simid, ntraj, traj_atom, i,real_time_measure)
            bcount_traj(i)=1
         enddo
      endif

   end subroutine flush_trajectories

   !---------------------------------------------------------------------------------
   ! SUBROUTINES: allocate_trajectories
   !> Allocate and initialize the variables needed for printing the trajectories
   !---------------------------------------------------------------------------------
   subroutine allocate_trajectories(Natom,Mensemble,flag)

      implicit none
      !
      integer, intent(in) :: flag      !< Allocate or deallocate (1/-1)
      integer, intent(in) :: Natom     !< Number of atoms in the system
      integer, intent(in) :: Mensemble !< Number of ensembles
      !
      integer :: i, i_stat, i_all
      !.. Allocate measurement counters
      if(flag>0) then

         if (ntraj>0) then
            allocate(scount_traj(ntraj),stat=i_stat)
            call memocc(i_stat,product(shape(scount_traj))*kind(scount_traj),'scount_traj','allocate_trajectories')
            scount_traj=0
            allocate(bcount_traj(ntraj),stat=i_stat)
            call memocc(i_stat,product(shape(bcount_traj))*kind(bcount_traj),'bcount_traj','allocate_trajectories')
            bcount_traj=0
            allocate(emomb_traj(3,maxval(traj_buff),ntraj,Mensemble),stat=i_stat)
            call memocc(i_stat,product(shape(emomb_traj))*kind(emomb_traj),'emomb_traj','allocate_trajectories')
            emomb_traj=0.0_dblprec
            allocate(mmomb_traj(maxval(traj_buff),ntraj,Mensemble),stat=i_stat)
            call memocc(i_stat,product(shape(mmomb_traj))*kind(mmomb_traj),'mmomb_traj','allocate_trajectories')
            mmomb_traj=0.0_dblprec
            allocate(indxb_traj(maxval(traj_buff),ntraj),stat=i_stat)
            call memocc(i_stat,product(shape(indxb_traj))*kind(indxb_traj),'indxb_traj','allocate_trajectories')
            indxb_traj=0
         endif

         if (do_tottraj=='Y') then
            allocate(mmomb(Natom,tottraj_buff,Mensemble),stat=i_stat)
            call memocc(i_stat,product(shape(mmomb))*kind(mmomb),'mmomb','allocate_trajectories')
            mmomb=0.0_dblprec
            allocate(emomb(3,Natom,tottraj_buff,Mensemble),stat=i_stat)
            call memocc(i_stat,product(shape(emomb))*kind(emomb),'emomb','allocate_trajectories')
            emomb=0.0_dblprec
            allocate(indxb(tottraj_buff),stat=i_stat)
            call memocc(i_stat,product(shape(indxb))*kind(indxb),'indxb','allocate_trajectories')
            indxb=0
         endif

         !.. Initiate trajectory measurements counters
         do i=1,ntraj
            scount_traj(i)=1
         enddo

         !.. Initiate trajectory buffer counters
         bcount_tottraj=1
         do i=1,ntraj
            bcount_traj(i)=1
         enddo
      else

         if (do_tottraj=='Y'.or. ntraj>0) then

            if (do_tottraj=='Y') then
               i_all=-product(shape(mmomb))*kind(mmomb)
               deallocate(mmomb,stat=i_stat)
               call memocc(i_stat,i_all,'mmomb','allocate_trajectories')
               i_all=-product(shape(emomb))*kind(emomb)
               deallocate(emomb,stat=i_stat)
               call memocc(i_stat,i_all,'emomb','allocate_trajectories')
               i_all=-product(shape(indxb))*kind(indxb)
               deallocate(indxb,stat=i_stat)
               call memocc(i_stat,i_all,'indxb','allocate_trajectories')
            endif

            if (ntraj>0)  then
               i_all=-product(shape(mmomb_traj))*kind(mmomb_traj)
               deallocate(mmomb_traj,stat=i_stat)
               call memocc(i_stat,i_all,'mmomb_traj','allocate_trajectories')
               i_all=-product(shape(emomb_traj))*kind(emomb_traj)
               deallocate(emomb_traj,stat=i_stat)
               call memocc(i_stat,i_all,'emomb_traj','allocate_trajectories')
               i_all=-product(shape(bcount_traj))*kind(bcount_traj)
               deallocate(bcount_traj,stat=i_stat)
               call memocc(i_stat,i_all,'bcount_traj','allocate_trajectories')
               i_all=-product(shape(scount_traj))*kind(scount_traj)
               deallocate(scount_traj,stat=i_stat)
               call memocc(i_stat,i_all,'scount_traj','allocate_trajectories')
               i_all=-product(shape(indxb_traj))*kind(indxb)
               deallocate(indxb_traj,stat=i_stat)
               call memocc(i_stat,i_all,'indxb_traj','allocate_trajectories')
            endif
         endif
      endif

   end subroutine allocate_trajectories

   !---------------------------------------------------------------------------------
   ! SUBROUTINE: buffer_tottraj
   !> Buffer all moments
   !---------------------------------------------------------------------------------
   subroutine buffer_tottraj(Natom,Mensemble,mstep,emom,mmom,bcount_tottraj,delta_t,&
      real_time_measure)

      implicit none

      integer, intent(in) :: mstep     !< Current simulation step
      integer, intent(in) :: Natom     !< Number of atoms in system
      integer, intent(in) :: Mensemble !< Number of ensembles 
      integer, intent(in) :: bcount_tottraj !< Counter of buffer for moments
      real(dblprec), intent(in) :: delta_t  !< Current time step (used for real time measurements)
      character(len=1), intent(in) :: real_time_measure   !< Real time measurement flag
      real(dblprec), dimension(Natom,Mensemble), intent(in) :: mmom   !< Magnitude of magnetic moments
      real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: emom !< Current unit moment vector

      !.. Local variables
      integer :: i,k

      do k=1, Mensemble
         do i=1, Natom
            emomb(1:3,i,bcount_tottraj,k)=emom(1:3,i,k)
            mmomb(i,bcount_tottraj,k)=mmom(i,k)
         end do
      end do

      if (real_time_measure=='Y') then
         indxb(bcount_tottraj)=mstep*delta_t
      else
         indxb(bcount_tottraj)=mstep
      endif

   end subroutine buffer_tottraj

   !---------------------------------------------------------------------------------
   ! SUBROUTINE: prn_tottraj
   !> Print all moments to file
   !---------------------------------------------------------------------------------
   subroutine prn_tottraj(Natom,Mensemble,simid,real_time_measure,mmom,             &
      do_mom_legacy,mode)
      !
      !.. Implicit declarations
      use restart, only : prn_mag_conf, prn_tottraj_legacy
      implicit none

      integer, intent(in) :: Natom     !< Number of atoms in system
      integer, intent(in) :: Mensemble !< Number of ensembles 
      character(len=1), intent(in) :: mode   !< Simulation mode (S=SD, M=MC, H=MC Heat Bath, P=LD, C=SLD, G=GNEB)
      character(len=8), intent(in) :: simid             !< Name of simulation 
      character(len=1), intent(in) :: do_mom_legacy      !< Flag to print/read moments in legacy output
      character(len=1), intent(in) :: real_time_measure !< Measurements displayed in real time
      real(dblprec), dimension(Natom,Mensemble), intent(in) :: mmom   !< Magnitude of magnetic moments

      !.. Local variables
      integer :: i

      !.. Executable statements

      if (do_mom_legacy.ne.'Y') then 
         if (real_time_measure=='Y') then
            do i=1, bcount_tottraj
               call prn_mag_conf(Natom,indxb(i),Mensemble,'M',simid,mmom,           &
               emomb(:,:,i,:),'',mode)
            enddo
         else
            do i=1, bcount_tottraj
               call prn_mag_conf(Natom,int(indxb(i)),Mensemble,'M',simid,mmom,      &
               emomb(:,:,i,:),'',mode)
            enddo
         endif
      else
         call prn_tottraj_legacy(Natom,Mensemble,tottraj_buff,bcount_tottraj,simid, &
            real_time_measure,indxb,emomb)
      endif
      return
      write (*,*) "Error writing the total trajectories file"

   end subroutine prn_tottraj

   !---------------------------------------------------------------------------------
   ! SUBROUTINE: buffer_traj
   !> Buffer selected trajectories
   !---------------------------------------------------------------------------------
   subroutine buffer_traj(Natom,Mensemble,mstep,emom,mmom,traj,ntraj,traj_atom,     &
      bcount,delta_t,real_time_measure)

      implicit none

      integer, intent(in) :: traj   !< Trajectory number
      integer, intent(in) :: ntraj  !< Number of trajectories to sample
      integer, intent(in) :: Natom  !< Number of atoms in system
      integer, intent(in) :: mstep  !< Current simulation step
      integer, intent(in) :: bcount !< Counter of buffer
      integer, intent(in) :: Mensemble !< Number of ensembles
      real(dblprec), intent(in) :: delta_t                !< Current time step (used for real time measurements)
      character(len=1), intent(in) :: real_time_measure   !< Real time measurement flag
      integer, intent(in), dimension(ntraj) :: traj_atom  !< List of atoms to sample trajectories for    
      real(dblprec), dimension(Natom,Mensemble) :: mmom   !< Magnitude of magnetic moments
      real(dblprec), dimension(3,Natom,Mensemble) :: emom !< Current unit moment vector

      integer :: k !< Current ensemble 

      do k=1, Mensemble
         emomb_traj(1:3,bcount,traj,k)=emom(1:3,traj_atom(traj),k)
         mmomb_traj(bcount,traj,k)=mmom(traj_atom(traj),k)
      end do

      if (real_time_measure=='Y') then
         indxb_traj(bcount,traj)=mstep*delta_t
      else
         indxb_traj(bcount,traj)=mstep
      endif

   end subroutine buffer_traj

   !---------------------------------------------------------------------------------
   ! SUBROUTINE: prn_traj
   !> Print selected trajectories
   !---------------------------------------------------------------------------------
   subroutine prn_traj(Mensemble, simid, ntraj, traj_atom, traj,real_time_measure)
      !
      !.. Implicit declarations
      implicit none

      integer, intent(in) :: traj      !< Trajectory number
      integer, intent(in) :: ntraj     !< Number of trajectories to sample
      integer, intent(in) :: Mensemble !< Number of ensembles 
      integer, intent(in), dimension(ntraj) :: traj_atom !< List of atoms to sample trajectories for
      character(len=8), intent(in) :: simid              !< Name of simulation 
      character(len=1), intent(in) :: real_time_measure  !< Measurements displayed in real time

      !.. Local variables
      integer :: i, j
      real(dblprec) ::  tempm, tempmx, tempmy, tempmz
      character(len=30) :: filn

      !.. Executable statements
      do j=1, Mensemble
         write (filn,'(''trajectory.'',a,''.'',i0.3,''.'',i1,''.out'')') trim(simid), traj, j

         open(ofileno, file=filn, position="append")
         do i=1, bcount_traj(traj)
            tempmx = emomb_traj(1,i,traj,j)
            tempmy = emomb_traj(2,i,traj,j)
            tempmz = emomb_traj(3,i,traj,j)
            tempm = mmomb_traj(i,traj,j)
            if (real_time_measure=='Y') then
               write (ofileno,10003) indxb_traj(i,traj), traj_atom(traj), &
                  tempmx, tempmy, tempmz, tempm
            else
               write (ofileno,10002) int(indxb_traj(i,traj)), traj_atom(traj), &
                  tempmx, tempmy, tempmz, tempm
            endif
         end do
         close(ofileno)
      end do
      return
      write (*,*) "Error writing the selected trajectories file"

      10002 format (i8,2x,i8,2x,2x, es16.8,es16.8,es16.8,es16.8)
      10003 format (es16.8,2x,i8,2x,2x, es16.8,es16.8,es16.8,es16.8)

   end subroutine prn_traj

   !---------------------------------------------------------------------------
   !> @brief
   !> Read input parameters.
   !
   !> @author
   !> Anders Bergman
   !---------------------------------------------------------------------------
   subroutine read_parameters_trajectories(ifile)
      use FileParser

      implicit none

      ! ... Formal Arguments ...
      integer, intent(in) :: ifile   !< File to read from
      !
      ! ... Local Variables ...
      character(len=50) :: keyword
      integer :: rd_len, i_err, i_errb, i_stat, i
      logical :: comment

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
            !------------------------------------------------------------------------
            ! START OF VARIABLES FOR PRINTING MOMENTS TRAJECTORIES
            !------------------------------------------------------------------------

            case('do_tottraj')
               read(ifile,*,iostat=i_err) do_tottraj
               if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err

            case('tottraj_step')
               read(ifile,*,iostat=i_err) tottraj_step
               if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err

            case('tottraj_buff')
               read(ifile,*,iostat=i_err) tottraj_buff
               if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err

            case('ntraj')
               read(ifile,*,iostat=i_err) ntraj
               if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err
               if(ntraj>0) then
                  allocate(traj_atom(ntraj),stat=i_stat)
                  call memocc(i_stat,product(shape(traj_atom))*kind(traj_atom),'traj_atom','read_parameters')
                  allocate(traj_step(ntraj),stat=i_stat)
                  call memocc(i_stat,product(shape(traj_step))*kind(traj_step),'traj_step','read_parameters')
                  allocate(traj_buff(ntraj),stat=i_stat)
                  call memocc(i_stat,product(shape(traj_buff))*kind(traj_buff),'traj_buff','read_parameters')
                  do i=1,ntraj
                     read(ifile,*,iostat=i_err) traj_atom(i), traj_step(i), traj_buff(i)
                  end do
               else
                  read(ifile,*)
               end if

            !------------------------------------------------------------------------
            ! END OF VARIABLES FOR PRINTING MOMENTS TRAJECTORIES
            !------------------------------------------------------------------------
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
   end subroutine read_parameters_trajectories

end module prn_trajectories
