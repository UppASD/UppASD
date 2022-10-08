!-------------------------------------------------------------------------------
!> @brief
!> Module to treat the prinitng of extra information regarding the induced moments
!> @details This module is very similar to what is obtained from prn_trajectories
!> The main difference is that some extra information for the induced moments
!> will be provided, all this info is also encoded in the trajectories, and
!> an avid user should be able to obtain it via clever scripting
!> @author
!> Jonathan Chico
!> @copyright
!> GNU Public License.
!-------------------------------------------------------------------------------
module prn_induced_info

   use Parameters
   use Profiling

   implicit none

   integer :: ind_step  !< Interval for sampling the induced magnetic moments
   integer :: ind_buff  !< Buffer size for the induced magnetic moments
   character(len=1) :: do_prn_induced !< Print extra information of the induced moments

   integer :: bcount_indtraj !< Counter of buffer for moments

   real(dblprec), dimension(:), allocatable :: ind_indxb            !< Step counter for the induced magnetic moments
   real(dblprec), dimension(:,:), allocatable :: sus_indb           !< Buffer for the susceptibility
   real(dblprec), dimension(:,:,:), allocatable :: neigh_modb       !< Buffer for the average moment of the neighbouring moment for the induced moments
   real(dblprec), dimension(:,:,:,:), allocatable :: ind_emomMb     !< Buffer for all individual trajectories
   real(dblprec), dimension(:,:,:,:), allocatable :: neigh_aveb     !< Buffer for the average magnitude of the nearest neighbours

   private

   public :: print_ind_trajectories,flush_ind_trajectories,allocate_ind_trajectories
   public :: do_prn_induced, ind_step, ind_buff,ind_traj_init

contains

   !----------------------------------------------------------------------------
   ! SUBROUTINE: print_ind_trajectories
   !> @brief Wrapper routine to print the trajectories of the induced moments
   !> @author Jonathan Chico
   !----------------------------------------------------------------------------
   subroutine print_ind_trajectories(Natom,Mensemble,sstep,mstep,max_no_neigh_ind,  &
         ind_nlistsize,ind_list_full,ind_nlist,delta_t,emomM,simid,                 &
         real_time_measure,sus_ind)

      implicit none

      integer, intent(in) :: Natom         !< Number of atoms in the system
      integer, intent(in) :: sstep         ! Simulation step in logarithmic scale
      integer, intent(in) :: mstep         !< Current simulation step
      integer, intent(in) :: Mensemble     !< Number of ensembles
      integer, intent(in) :: max_no_neigh_ind   !< Calculated maximum of neighbours for induced moments
      integer, dimension(Natom),intent(in) :: ind_nlistsize !< Size of neighbour list for induced moments
      integer, dimension(Natom), intent(in) :: ind_list_full   !< Indication of whether a given moment is induced/fixed (1/0) for all atoms
      integer, dimension(max_no_neigh_ind,Natom), intent(in) :: ind_nlist  !< Neighbour list for iduced moments
      real(dblprec), intent(in) :: delta_t !< Time step for real time measurement
      real(dblprec), dimension(Natom), intent(in) :: sus_ind
      real(dblprec), dimension(3,Natom, Mensemble), intent(in) :: emomM   !< Current magnetic moment vector
      character(len=8), intent(in) :: simid             !< Simulation ID
      character(len=1), intent(in) :: real_time_measure !< Measurements displayed in real time

      !  Total trajectory
      if (do_prn_induced=='Y') then

         if (mod(sstep-1,ind_step)==0) then

            ! Write step to buffer
            call buffer_indtraj(mstep-1,Natom,Mensemble,max_no_neigh_ind,bcount_indtraj,ind_nlistsize,ind_nlist,&
               delta_t,emomM,real_time_measure,sus_ind)
            if (bcount_indtraj==ind_buff) then

               !Write buffer to file
               call prn_indtraj(Natom, Mensemble,ind_list_full,simid,real_time_measure)
               bcount_indtraj=1
            else
               bcount_indtraj=bcount_indtraj+1
            endif
         else

         endif
      endif

   end subroutine print_ind_trajectories

   !----------------------------------------------------------------------------
   ! SUBROUTINE: ind_traj_init
   !> @brief Initialization of variables with default variables for the trajectories
   !> @author Jonathan Chico
   !----------------------------------------------------------------------------
   subroutine ind_traj_init()

      implicit none

      do_prn_induced   = 'N'
      ind_step = 1000
      ind_buff = 10

   end subroutine ind_traj_init

   !----------------------------------------------------------------------------
   ! SUBROUTINE: flush_ind_trajectories
   !> @brief Flush the trajectory measurements, i.e. print to file in the last iteration
   !> @author Jonathan Chico
   !----------------------------------------------------------------------------
   subroutine flush_ind_trajectories(Natom,Mensemble,ind_list_full,simid,           &
      real_time_measure)

      implicit none

      integer, intent(in) :: Natom     !< Number of atoms in the system
      integer, intent(in) :: Mensemble !< Number of ensembles
      integer, dimension(Natom), intent(in) :: ind_list_full   !< Indication of whether a given moment is induced/fixed (1/0) for all atoms
      character(len=8), intent(in) :: simid             !< Simulation name ID
      character(len=1), intent(in) :: real_time_measure !< Perfomr measurements in real time

      ! All trajectories are printed to file in the last iteration
      if (do_prn_induced=='Y') then
         !Write buffer to file
         bcount_indtraj=bcount_indtraj-1
         call prn_indtraj(Natom, Mensemble,ind_list_full,simid,real_time_measure)
      endif

   end subroutine flush_ind_trajectories

   !----------------------------------------------------------------------------
   ! SUBROUTINE: allocate_ind_trajectories
   !> @brief Allocate and initialize the variables needed for printing the trajectories
   !> @author Jonathan Chico
   !----------------------------------------------------------------------------
   subroutine allocate_ind_trajectories(Natom,Mensemble,flag)

      implicit none
      !
      integer, intent(in) :: flag      !< Allocate or deallocate (1/-1)
      integer, intent(in) :: Natom     !< Number of atoms in the system
      integer, intent(in) :: Mensemble !< Number of ensembles
      !
      integer :: i_stat, i_all
      !.. Allocate measurement counters
      if(flag>0) then

         if (do_prn_induced=='Y') then
            allocate(neigh_modb(Natom,ind_buff,Mensemble),stat=i_stat)
            call memocc(i_stat,product(shape(neigh_modb))*kind(neigh_modb),'neigh_modb','allocate_ind_trajectories')
            allocate(ind_emomMb(3,Natom,ind_buff,Mensemble),stat=i_stat)
            call memocc(i_stat,product(shape(ind_emomMb))*kind(ind_emomMb),'ind_emomMb','allocate_ind_trajectories')
            allocate(neigh_aveb(3,Natom,ind_buff,Mensemble),stat=i_stat)
            call memocc(i_stat,product(shape(neigh_aveb))*kind(neigh_aveb),'neigh_aveb','allocate_ind_trajectories')
            allocate(ind_indxb(ind_buff),stat=i_stat)
            call memocc(i_stat,product(shape(ind_indxb))*kind(ind_indxb),'ind_indxb','allocate_ind_trajectories')
            allocate(sus_indb(Natom,ind_buff),stat=i_stat)
            call memocc(i_stat,product(shape(sus_indb))*kind(sus_indb),'sus_indb','allocate_ind_trajectories')

            ! Start the counter
            bcount_indtraj=1
         endif
      else
         if (do_prn_induced=='Y') then

            i_all=-product(shape(neigh_modb))*kind(neigh_modb)
            deallocate(neigh_modb,stat=i_stat)
            call memocc(i_stat,i_all,'neigh_modb','allocate_ind_trajectories')
            i_all=-product(shape(ind_emomMb))*kind(ind_emomMb)
            deallocate(ind_emomMb,stat=i_stat)
            call memocc(i_stat,i_all,'ind_emomMb','allocate_ind_trajectories')
            i_all=-product(shape(neigh_aveb))*kind(neigh_aveb)
            deallocate(neigh_aveb,stat=i_stat)
            call memocc(i_stat,i_all,'neigh_aveb','allocate_ind_trajectories')
            i_all=-product(shape(ind_indxb))*kind(ind_indxb)
            deallocate(ind_indxb,stat=i_stat)
            call memocc(i_stat,i_all,'ind_indxb','allocate_ind_trajectories')
            i_all=-product(shape(sus_indb))*kind(sus_indb)
            deallocate(sus_indb,stat=i_stat)
            call memocc(i_stat,i_all,'sus_indb','allocate_ind_trajectories')

         endif
      endif

   end subroutine allocate_ind_trajectories

   !----------------------------------------------------------------------------
   ! SUBROUTINE: buffer_indtraj
   !> @brief Buffer all moments
   !> @author Jonathan Chico
   !----------------------------------------------------------------------------
   subroutine buffer_indtraj(mstep,Natom,Mensemble,max_no_neigh_ind,bcount_indtraj, &
      ind_nlistsize,ind_nlist,delta_t,emomM,real_time_measure,sus_ind)

      implicit none

      integer, intent(in) :: mstep     !< Current simulation step
      integer, intent(in) :: Natom     !< Number of atoms in system
      integer, intent(in) :: Mensemble !< Number of ensembles
      integer, intent(in) :: max_no_neigh_ind   !< Calculated maximum of neighbours for induced moments
      integer, intent(in) :: bcount_indtraj !< Counter of buffer for the induced moments
      integer, dimension(Natom), intent(in):: ind_nlistsize !< Size of neighbour list for induced moments
      integer, dimension(max_no_neigh_ind,Natom), intent(in) :: ind_nlist !< Neighbour list for induced moments
      real(dblprec), intent(in) :: delta_t  !< Current time step (used for real time measurements)
      real(dblprec), dimension(Natom), intent(in) :: sus_ind
      real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: emomM !< Current magnetic moment vector
      character(len=1), intent(in) :: real_time_measure   !< Real time measurement flag

      !.. Local variables
      integer :: ii,kk,neigh,fix
      real(dblprec) :: mod_ave
      real(dblprec), dimension(3) :: ave_mom

      do kk=1, Mensemble
         do ii=1, Natom
            ind_emomMb(1:3,ii,bcount_indtraj,kk)=emomM(1:3,ii,kk)
            ave_mom=0.0_dblprec
            do neigh=1, ind_nlistsize(ii)
               fix=ind_nlist(neigh,ii)
               ave_mom(1)=ave_mom(1)+emomM(1,fix,kk)
               ave_mom(2)=ave_mom(2)+emomM(2,fix,kk)
               ave_mom(3)=ave_mom(3)+emomM(3,fix,kk)
            enddo
            mod_ave=norm2(ave_mom)
            neigh_aveb(1:3,ii,bcount_indtraj,kk)=ave_mom(1:3)
            neigh_modb(ii,bcount_indtraj,kk)=mod_ave/ind_nlistsize(ii)
            sus_indb(ii,bcount_indtraj)=sus_ind(ii)
         end do
      end do

      if (real_time_measure=='Y') then
         ind_indxb(bcount_indtraj)=mstep*delta_t
      else
         ind_indxb(bcount_indtraj)=mstep
      endif

   end subroutine buffer_indtraj

   !----------------------------------------------------------------------------
   ! SUBROUTINE: prn_indtraj
   !> @brief Print the induced magnetic moments to file
   !> @author Jonathan Chico
   !----------------------------------------------------------------------------
   subroutine prn_indtraj(Natom, Mensemble,ind_list_full,simid,real_time_measure)
      !
      !.. Implicit declarations
      implicit none

      integer, intent(in) :: Natom     !< Number of atoms in system
      integer, intent(in) :: Mensemble !< Number of ensembles
      integer, dimension(Natom),intent(in) :: ind_list_full !< Indication of whether a given moment is induced/fixed (1/0) for all atoms
      character(len=8), intent(in) :: simid             !< Name of simulation
      character(len=1), intent(in) :: real_time_measure !< Measurements displayed in real time

      !.. Local variables
      integer :: ii, jj
      character(len=30) :: filn

      !.. Executable statements
      write (filn,'(''ind_moment.'',a,''.out'')') trim(simid)
      open(ofileno, file=filn, position="append")
      do ii=1, bcount_indtraj
         do jj=1, Natom
            if (real_time_measure=='Y') then
               if (ind_list_full(jj)==1) then
                  write (ofileno,10043) ind_indxb(ii), jj, norm2(ind_emomMb(:,jj,ii,Mensemble)),&
                  neigh_aveb(1,jj,ii,Mensemble),neigh_aveb(2,jj,ii,Mensemble),neigh_aveb(3,jj,ii,Mensemble),&
                  neigh_modb(jj,ii,Mensemble),sus_indb(jj,ii)
               endif
            else
               if (ind_list_full(jj)==1) then
                  write (ofileno,10042) int(ind_indxb(ii)), jj,norm2(ind_emomMb(:,jj,ii,Mensemble)), &
                  neigh_aveb(1,jj,ii,Mensemble),neigh_aveb(2,jj,ii,Mensemble),neigh_aveb(3,jj,ii,Mensemble),&
                  neigh_modb(jj,ii,Mensemble),sus_indb(jj,ii)
               endif
            endif
         end do
      end do
      close(ofileno)
      return
      write (*,*) "Error writing the induced moment trajectories file"
      10042 format (i8,2x,i8,2x,2x, es16.8E3,2x,es16.8E3,2x,es16.8E3,2x,es16.8E3,2x,es16.8E3,2x,es16.8E3)
      10043 format (es16.4,2x,i8,2x,2x, es16.8E3,2x,es16.8E3,2x,es16.8E3,2x,es16.8E3,2x,es16.8E3)

   end subroutine prn_indtraj

end module prn_induced_info
