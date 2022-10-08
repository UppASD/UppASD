!> Data and routines necessary for printing the trajectories of individual atomic displacements and velocities
module prn_latticetrajectories

   use Parameters
   use Profiling

   implicit none

   ! SLDTODO Here some variables are passed as arguments in functions calls even though
   ! SLDTODO they are already within the scope. Make consistent.

   ! Input parameters
   integer :: lntraj         !< Number of displacement trajectories to sample
   integer :: ltottraj_step  !< Interval for sampling displacement trajectories
   integer :: ltottraj_buff  !< Buffer size for displacement trajectories
   integer, dimension(:), allocatable :: ltraj_step !< Interval for sampling individual displacement trajectories
   integer, dimension(:), allocatable :: ltraj_buff !< Buffer size for individual displacement trajectories
   integer, dimension(:), allocatable :: ltraj_atom !< List of atoms to sample displacement trajectories for
   character(len=1) :: do_ltottraj !< Measure displacements

   ! Measurement variables
   integer :: bcount_tottraj !< Counter of buffer for displacements
   integer, dimension(:), allocatable :: bcount_traj   !< Counter of buffer of displacements
   real(dblprec), dimension(:), allocatable :: scount_traj      !< Counter of sampling of displacements
   real(dblprec), dimension(:), allocatable :: indxb            !< Step counter for displacements
   real(dblprec), dimension(:,:), allocatable :: indxb_traj     !< Step counter for individual displacements
   real(dblprec), dimension(:,:,:,:), allocatable :: uvecb      !< Buffer for all individual displacements
   real(dblprec), dimension(:,:,:,:), allocatable :: uvecb_traj !< Buffer for selected individual displacements
   real(dblprec), dimension(:,:,:,:), allocatable :: vvecb      !< Buffer for all individual velocities
   real(dblprec), dimension(:,:,:,:), allocatable :: vvecb_traj !< Buffer for selected individual velocities
   real(dblprec), dimension(:,:,:,:), allocatable :: llatb      !< Buffer for all individual angular momenta
   real(dblprec), dimension(:,:,:,:), allocatable :: llatb_traj !< Buffer for selected individual angular momenta

   public

   private :: bcount_tottraj, bcount_traj, scount_traj, indxb, indxb_traj, uvecb, uvecb_traj, vvecb, vvecb_traj
   private :: llatb, llatb_traj


contains


   !> Wrapper routine to print all the displacement trajectories 
   subroutine print_latttrajectories(simid, Natom, Mensemble, mstep, sstep, mion, uvec, vvec, coord)

      implicit none

      character(len=8), intent(in) :: simid             !< Simulation ID 
      integer, intent(in) :: Natom         !< Number of atoms in the system
      integer, intent(in) :: Mensemble     !< Number of ensembles
      integer, intent(in) :: mstep         !< Current simulation step
      integer, intent(in) :: sstep         ! Simulation step in logarithmic scale
      real(dblprec), dimension(Natom,Mensemble), intent(in) :: mion   !< Ionic mass
      real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: uvec   !< Current ionic displacement
      real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: vvec   !< Current ionic velocity
      real(dblprec), dimension(3,Natom), intent(in) :: coord !< Coordinates of atoms

      !.. Local variables
      integer :: i

      !  Total trajectory
      if (do_ltottraj=='Y') then

         if (mod(sstep-1,ltottraj_step)==0) then

            ! Write step to buffer
            call buffer_tottraj(Natom, Mensemble, mstep-1, bcount_tottraj, uvec, vvec)

            if (bcount_tottraj==ltottraj_buff) then
               !Write buffer to file
               call prn_tottraj(simid, Natom, Mensemble,coord)
               bcount_tottraj=1
            else
               bcount_tottraj=bcount_tottraj+1
            endif

         endif
      endif

      !  Atom trajectory
      if(lntraj>0) then
         do i=1,lntraj
            if (mod(sstep-1,ltraj_step(i))==0) then

               ! Write step to buffer
               call buffer_traj(Natom, Mensemble, mstep-1, uvec, vvec, i, lntraj, ltraj_atom, bcount_traj(i))
               scount_traj(i)=1

               if (bcount_traj(i)==ltraj_buff(i)) then
                  call prn_traj(simid, Mensemble, lntraj, ltraj_atom, i)
                  bcount_traj(i)=1
               else
                  bcount_traj(i)=bcount_traj(i)+1
               endif
            else
               scount_traj(i)=scount_traj(i)+1
            end if
         enddo
      endif

   end subroutine print_latttrajectories


   !> Initialization of variables with default variables for the trajectories
   subroutine latttraj_init()

      implicit none

      lntraj        = 0
      do_ltottraj   = 'N'
      ltottraj_step = 1000
      ltottraj_buff = 10

   end subroutine latttraj_init


   !> Flush the trajectory measurements, i.e. print to file in the last iteration
   subroutine flush_latttrajectories(simid, Natom, Mensemble,coord)

      implicit none

      character(len=8), intent(in) :: simid             !< Simulation name ID
      integer, intent(in) :: Natom     !< Number of atoms in the system
      integer, intent(in) :: Mensemble !< Number of ensembles
      real(dblprec), dimension(3,Natom), intent(in) :: coord !< Coordinates of atoms

      !.. Local variables
      integer :: i

      ! All trajectories are printed to file in the last iteration
      if (do_ltottraj=='Y') then
         !Write buffer to file  
         bcount_tottraj=bcount_tottraj-1
         call prn_tottraj(simid, Natom, Mensemble,coord)
      endif

      !  Selected trajectories are printed to file in the last iteration
      if(lntraj>0) then
         ! Write buffer to file
         do i=1,lntraj
            bcount_traj(i)=bcount_traj(i)-1
            call prn_traj(simid, Mensemble, lntraj, ltraj_atom, i)
         enddo
      endif

   end subroutine flush_latttrajectories


   !> Allocate and initialize the variables needed for printing the displacement trajectories
   subroutine allocate_latttrajectories(Natom, Mensemble, flag)

      implicit none
      !
      integer, intent(in) :: Natom     !< Number of atoms in the system
      integer, intent(in) :: Mensemble !< Number of ensembles
      integer, intent(in) :: flag      !< Allocate or deallocate (1/-1)
      !
      integer :: i, i_stat, i_all
      !.. Allocate measurement counters
      if(flag>0) then

         if (lntraj>0) then
            allocate(scount_traj(lntraj),stat=i_stat)
            call memocc(i_stat,product(shape(scount_traj))*kind(scount_traj),'scount_traj','allocate_trajectories')
            allocate(bcount_traj(lntraj),stat=i_stat)
            call memocc(i_stat,product(shape(bcount_traj))*kind(bcount_traj),'bcount_traj','allocate_trajectories')
            allocate(uvecb_traj(3,maxval(ltraj_buff),lntraj,Mensemble),stat=i_stat)
            call memocc(i_stat,product(shape(uvecb_traj))*kind(uvecb_traj),'uvecb_traj','allocate_trajectories')
            allocate(vvecb_traj(3,maxval(ltraj_buff),lntraj,Mensemble),stat=i_stat)
            call memocc(i_stat,product(shape(vvecb_traj))*kind(vvecb_traj),'vvecb_traj','allocate_trajectories')
            allocate(llatb_traj(3,maxval(ltraj_buff),lntraj,Mensemble),stat=i_stat)
            call memocc(i_stat,product(shape(llatb_traj))*kind(llatb_traj),'llatb_traj','allocate_trajectories')
            allocate(indxb_traj(maxval(ltraj_buff),lntraj),stat=i_stat)
            call memocc(i_stat,product(shape(indxb_traj))*kind(indxb_traj),'indxb_traj','allocate_trajectories')
         endif

         if (do_ltottraj=='Y') then
            allocate(uvecb(3,Natom,ltottraj_buff,Mensemble),stat=i_stat)
            call memocc(i_stat,product(shape(uvecb))*kind(uvecb),'emomb','allocate_trajectories')
            allocate(vvecb(3,Natom,ltottraj_buff,Mensemble),stat=i_stat)
            call memocc(i_stat,product(shape(vvecb))*kind(vvecb),'vvecb','allocate_trajectories')
            allocate(llatb(3,Natom,ltottraj_buff,Mensemble),stat=i_stat)
            call memocc(i_stat,product(shape(llatb))*kind(llatb),'llatb','allocate_trajectories')
            allocate(indxb(ltottraj_buff),stat=i_stat)
            call memocc(i_stat,product(shape(indxb))*kind(indxb),'indxb','allocate_trajectories')
         endif

         !.. Initiate trajectory measurements counters
         do i=1,lntraj
            scount_traj(i)=1
         enddo

         !.. Initiate trajectory buffer counters
         bcount_tottraj=1
         do i=1,lntraj
            bcount_traj(i)=1
         enddo

      else

         if (do_ltottraj=='Y'.or. lntraj>0) then

            if (do_ltottraj=='Y') then
               i_all=-product(shape(uvecb))*kind(uvecb)
               deallocate(uvecb,stat=i_stat)
               call memocc(i_stat,i_all,'uvecb','allocate_trajectories')
               i_all=-product(shape(vvecb))*kind(vvecb)
               deallocate(vvecb,stat=i_stat)
               call memocc(i_stat,i_all,'vvecb','allocate_trajectories')
               i_all=-product(shape(llatb))*kind(llatb)
               deallocate(llatb,stat=i_stat)
               call memocc(i_stat,i_all,'llatb','allocate_trajectories')
               i_all=-product(shape(indxb))*kind(indxb)
               deallocate(indxb,stat=i_stat)
               call memocc(i_stat,i_all,'indxb','allocate_trajectories')
            endif

            if (lntraj>0)  then
               i_all=-product(shape(uvecb_traj))*kind(uvecb_traj)
               deallocate(uvecb_traj,stat=i_stat)
               call memocc(i_stat,i_all,'uvecb_traj','allocate_trajectories')
               i_all=-product(shape(vvecb_traj))*kind(vvecb_traj)
               deallocate(vvecb_traj,stat=i_stat)
               call memocc(i_stat,i_all,'vvecb_traj','allocate_trajectories')
               i_all=-product(shape(llatb_traj))*kind(llatb_traj)
               deallocate(llatb_traj,stat=i_stat)
               call memocc(i_stat,i_all,'llatb_traj','allocate_trajectories')
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


   end subroutine allocate_latttrajectories


   !> Buffer all displacements
   subroutine buffer_tottraj(Natom, Mensemble, mstep, bcount_tottraj, uvec, vvec)

      use math_functions, only : f_cross_product
      use prn_latticeaverages,only : dcoord
      use latticedata, only : mion

      implicit none


      integer, intent(in) :: Natom     !< Number of atoms in system
      integer, intent(in) :: Mensemble !< Number of ensembles 
      integer, intent(in) :: mstep     !< Current simulation step
      integer, intent(in) :: bcount_tottraj !< Counter of buffer for moments
      real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: uvec   !< Current ionic displacement
      real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: vvec   !< Current ionic velocity

      !.. Local variables
      integer :: i, k
      real(dblprec), dimension(3) :: cdum

      do k=1, Mensemble
         do i=1, Natom
            uvecb(1:3,i,bcount_tottraj,k)=uvec(1:3,i,k)
            vvecb(1:3,i,bcount_tottraj,k)=vvec(1:3,i,k)
            cdum=0*dcoord(1:3,i)+uvec(1:3,i,k)
            llatb(1:3,i,bcount_tottraj,k)=f_cross_product(cdum,vvec(1:3,i,k))*mion(i,k)
         end do
      end do

      indxb(bcount_tottraj)=mstep

   end subroutine buffer_tottraj


   !> Print all displacements to file
   subroutine prn_tottraj(simid, Natom, Mensemble,coord)
      use prn_latticeaverages,only : dcoord
      use Constants, only: amu,hbar,angstrom
      !
      !.. Implicit declarations
      implicit none

      character(len=8), intent(in) :: simid             !< Name of simulation 
      integer, intent(in) :: Natom     !< Number of atoms in system
      integer, intent(in) :: Mensemble !< Number of ensembles 
      real(dblprec), dimension(3,Natom), intent(in) :: coord !< Coordinates of atoms

      !.. Local variables
      integer :: i, j
      character(len=30) :: filn, filn2
      real(dblprec) :: spinfac

      !.. Executable statements
      write (filn,'(''disp.'',a,''.out'')') trim(simid)
      open(ofileno, file=filn, position="append")
      do i=1, bcount_tottraj
         do j=1, Natom
            write (ofileno,10002) int(indxb(i)), j, uvecb(1,j,i,Mensemble), uvecb(2,j,i,Mensemble), uvecb(3,j,i,Mensemble), &
               vvecb(1,j,i,Mensemble), vvecb(2,j,i,Mensemble), vvecb(3,j,i,Mensemble), llatb(1:3,j,i,Mensemble)
         end do
      end do
      close(ofileno)

      spinfac=0.5_dblprec*amu*angstrom/hbar*angstrom
      write (filn,'(''gdisp.'',a,''.out'')') trim(simid)
      open(ofileno, file=filn, position="append")
      do i=1, bcount_tottraj
         do j=1, Natom
            write (ofileno,10002) int(indxb(i)), j, dcoord(:,j)+uvecb(:,j,i,Mensemble), vvecb(:,j,i,Mensemble), &
                llatb(1:3,j,i,Mensemble)*spinfac
         end do
      end do
      close(ofileno)

      write (filn2,'(''rcoord.'',a,''.out'')') trim(simid)
      open(ofileno, file=filn2, position="append")
      do i=1, bcount_tottraj
         do j=1, Natom
            write (ofileno,10003) j, coord(1,j) + uvecb(1,j,i,Mensemble), &
               coord(2,j) + uvecb(2,j,i,Mensemble), coord(3,j) +  uvecb(3,j,i,Mensemble)
         end do
      end do
      close(ofileno)

      return

      write (*,*) "Error writing the total trajectories file"

      !10002 format (i8,2x,i8,2x,2x,es16.8E3,2x,es16.8E3,2x,es16.8E3,2x,es16.8E3,2x,es16.8E3,2x,es16.8E3,2x,es16.8E3,2x,es16.8E3)
      10002 format (i8,2x,i8,4x,3es16.6,4x,3es16.6,4x,3es16.6)
      10003 format (i8,4x,3es16.6,4x,3es16.6,4x,3es16.6)
      !10003 format (i8,2x,es16.8E3,2x,es16.8E3,2x,es16.8E3,2x,es16.8E3,2x,es16.8E3,2x,es16.8E3,2x,es16.8E3,2x,es16.8E3)

   end subroutine prn_tottraj


   !> Buffer selected trajectories
   subroutine buffer_traj(Natom, Mensemble, mstep, uvec, vvec, traj, lntraj, ltraj_atom, bcount)

      use math_functions, only : f_cross_product
      use prn_latticeaverages,only : dcoord
      use latticedata, only : mion

      implicit none

      integer, intent(in) :: Natom  !< Number of atoms in system
      integer, intent(in) :: Mensemble !< Number of ensembles 
      integer, intent(in) :: mstep  !< Current simulation step
      real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: uvec   !< Current ionic displacement
      real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: vvec   !< Current ionic velocity
      integer, intent(in) :: traj   !< Trajectory number
      integer, intent(in) :: lntraj  !< Number of trajectories to sample
      integer, intent(in), dimension(lntraj) :: ltraj_atom  !< List of atoms to sample trajectories for    
      integer, intent(in) :: bcount !< Counter of buffer

      !.. Local variables
      integer :: k !< Current ensemble 
      real(dblprec), dimension(3) :: cdum

      do k=1, Mensemble
         uvecb_traj(1:3,bcount,traj,k)=uvec(1:3,ltraj_atom(traj),k)
         vvecb_traj(1:3,bcount,traj,k)=vvec(1:3,ltraj_atom(traj),k)
         ! Notice the zero below. Currently local L only.
         cdum=0*dcoord(:,ltraj_atom(traj))+uvec(1:3,ltraj_atom(traj),k)
         llatb_traj(1:3,bcount,traj,k)=f_cross_product(cdum,vvec(1:3,ltraj_atom(traj),k))*mion(ltraj_atom(traj),k)
      end do

      indxb_traj(bcount,traj)=mstep

   end subroutine buffer_traj


   !> Print selected trajectories
   subroutine prn_traj(simid, Mensemble, lntraj, ltraj_atom, traj)
      !
      !.. Implicit declarations
      implicit none

      character(len=8), intent(in) :: simid    !< Name of simulation 
      integer, intent(in) :: Mensemble !< Number of ensembles 
      integer, intent(in) :: lntraj     !< Number of trajectories to sample
      integer, intent(in), dimension(lntraj) :: ltraj_atom !< List of atoms to sample trajectories for
      integer, intent(in) :: traj      !< Trajectory number

      !.. Local variables
      integer :: i, j
      real(dblprec) ::  tempux, tempuy, tempuz, tempvx, tempvy, tempvz, templx, temply, templz
      character(len=30) :: filn

      !.. Executable statements
      do j=1, Mensemble
         write (filn,'(''disptraj.'',a,''.'',i1,''.'',i1,''.out'')') trim(simid), traj, j
         open(ofileno, file=filn, position="append")
         do i=1, bcount_traj(traj)
            tempux = uvecb_traj(1,i,traj,j)
            tempuy = uvecb_traj(2,i,traj,j)
            tempuz = uvecb_traj(3,i,traj,j)
            tempvx = vvecb_traj(1,i,traj,j)
            tempvy = vvecb_traj(2,i,traj,j)
            tempvz = vvecb_traj(3,i,traj,j)
            templx = llatb_traj(1,i,traj,j)
            temply = llatb_traj(2,i,traj,j)
            templz = llatb_traj(3,i,traj,j)
            write (ofileno,10002) int(indxb_traj(i,traj)), ltraj_atom(traj), &
               tempux, tempuy, tempuz, tempvx, tempvy, tempvz, templx,temply,templz
         end do
         close(ofileno)
      end do
      return
      write (*,*) "Error writing the selected trajectories file"

      10002 format (i8,2x,i8,4x,3es16.8,4x,3es16.8,4x,3es16.8)
      !10002 format (i8,2x,i8,2x,2x,8es16.8)

   end subroutine prn_traj


end module prn_latticetrajectories
