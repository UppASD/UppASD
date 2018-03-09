!====================================================================!
!> @brief
!> Module to treat the prinitng of extra information regarding the induced moments
!
!> @author
!> Jonathan Chico
!> @copyright
!! Copyright (C) 2008-2018 UppASD group
!! This file is distributed under the terms of the
!! GNU General Public License. 
!! See http://www.gnu.org/copyleft/gpl.txt
!====================================================================!
! This module is very similar to what is obtained from prn_trajectories
! The main difference is that some extra information for the induced moments
! will be provided, all this info is also encoded in the trajectories, and
! an avid user should be able to obtain it via clever scripting
module prn_induced_info

  use Parameters
  use Profiling

  implicit none

  integer :: ind_step  !< Interval for sampling the induced magnetic moments
  integer :: ind_buff  !< Buffer size for the induced magnetic moments
  character(len=1) :: do_prn_induced !< Print extra information of the induced moments

  integer :: bcount_indtraj !< Counter of buffer for moments
  integer :: scount_tottraj !< Counter of sampling for moments

  real(dblprec), dimension(:), allocatable :: ind_indxb            !< Step counter for the induced magnetic moments
  real(dblprec), dimension(:,:,:,:), allocatable :: ind_emomb      !< Buffer for all individual trajectories
  real(dblprec), dimension(:,:,:,:), allocatable :: neigh_aveb     !< Buffer for the average magnitude of the nearest neighbours
  real(dblprec), dimension(:,:,:), allocatable :: ind_mmomb        !< Buffer for all moment magnitudes for induced moments
  real(dblprec), dimension(:,:,:), allocatable :: neigh_modb       !< Buffer for the average moment of the neighbouring moment for the induced moments

  private

  public :: print_ind_trajectories, flush_ind_trajectories, allocate_ind_trajectories, ind_traj_init
  public :: do_prn_induced, ind_step, ind_buff

contains

  !> Wrapper routine to print all the trajectories
  subroutine print_ind_trajectories(NA,Natom,sstep,mstep,Nchmax,Mensemble,do_ralloy,Natom_full,max_no_neigh,&
             atype,ind_nlistsize,atype_ch,ind_mom,ind_nlist,delta_t,mmom,emom,simid,real_time_measure)

    implicit none

    integer, intent(in) :: NA
    integer, intent(in) :: Natom         !< Number of atoms in the system
    integer, intent(in) :: sstep         ! Simulation step in logarithmic scale
    integer, intent(in) :: mstep         !< Current simulation step
    integer, intent(in) :: Nchmax
    integer, intent(in) :: Mensemble     !< Number of ensembles
    integer, intent(in) :: do_ralloy
    integer, intent(in) :: Natom_full
    integer, intent(in) :: max_no_neigh
    integer, dimension(Natom), intent(in) :: atype
    integer, dimension(Natom),intent(in) :: ind_nlistsize
    integer, dimension(Natom_full), intent(in) :: atype_ch
    integer, dimension(NA,Nchmax), intent(in) :: ind_mom
    integer, dimension(max_no_neigh,Natom), intent(in) :: ind_nlist
    real(dblprec), intent(in) :: delta_t !< Time step for real time measurement
    real(dblprec), dimension(Natom, Mensemble), intent(in) :: mmom     !< Magnitude of magnetic moments
    real(dblprec), dimension(3,Natom, Mensemble), intent(in) :: emom   !< Current unit moment vector
    character(len=8), intent(in) :: simid             !< Simulation ID
    character(len=1), intent(in) :: real_time_measure !< Measurements displayed in real time

    !  Total trajectory
    if (do_prn_induced=='Y') then

       if (mod(sstep-1,ind_step)==0) then

          ! Write step to buffer
          call buffer_indtraj(mstep-1,Natom,Mensemble,max_no_neigh,bcount_indtraj,ind_nlistsize,ind_nlist,&
               delta_t,mmom,emom,real_time_measure)
          if (bcount_indtraj==ind_buff) then

             !Write buffer to file
             call prn_indtraj(NA,Natom, Nchmax,Mensemble,do_ralloy,Natom_full,atype,atype_ch,ind_mom,&
                  simid,real_time_measure)
             bcount_indtraj=1
          else
             bcount_indtraj=bcount_indtraj+1
          endif
       else

       endif
    endif

  end subroutine print_ind_trajectories

  !> Initialization of variables with default variables for the trajectories
  subroutine ind_traj_init()

     implicit none

     do_prn_induced   = 'N'
     ind_step = 1000
     ind_buff = 10


  end subroutine ind_traj_init

  !> Flush the trajectory measurements, i.e. print to file in the last iteration
  subroutine flush_ind_trajectories(NA,Natom,Nchmax,Mensemble,do_ralloy,Natom_full,atype,&
             atype_ch,ind_mom,simid,real_time_measure)

    implicit none

    integer, intent(in) :: NA
    integer, intent(in) :: Natom     !< Number of atoms in the system
    integer, intent(in) :: Nchmax
    integer, intent(in) :: Mensemble !< Number of ensembles
    integer, intent(in) :: do_ralloy
    integer, intent(in) :: Natom_full
    integer, dimension(Natom), intent(in) :: atype
    integer, dimension(Natom_full), intent(in) :: atype_ch
    integer, dimension(NA,Nchmax), intent(in) :: ind_mom
    character(len=8), intent(in) :: simid             !< Simulation name ID
    character(len=1), intent(in) :: real_time_measure !< Perfomr measurements in real time

    ! All trajectories are printed to file in the last iteration
    if (do_prn_induced=='Y') then
       !Write buffer to file
       bcount_indtraj=bcount_indtraj-1
       call prn_indtraj(NA,Natom, Nchmax,Mensemble,do_ralloy,Natom_full,atype,atype_ch,ind_mom,&
            simid,real_time_measure)
    endif


  end subroutine flush_ind_trajectories

  !> Allocate and initialize the variables needed for printing the trajectories
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
          allocate(ind_mmomb(Natom,ind_buff,Mensemble),stat=i_stat)
          call memocc(i_stat,product(shape(ind_mmomb))*kind(ind_mmomb),'ind_mmomb','allocate_ind_trajectories')
          allocate(neigh_modb(Natom,ind_buff,Mensemble),stat=i_stat)
          call memocc(i_stat,product(shape(neigh_modb))*kind(neigh_modb),'neigh_modb','allocate_ind_trajectories')
          allocate(ind_emomb(3,Natom,ind_buff,Mensemble),stat=i_stat)
          call memocc(i_stat,product(shape(ind_emomb))*kind(ind_emomb),'ind_emomb','allocate_ind_trajectories')
          allocate(neigh_aveb(3,Natom,ind_buff,Mensemble),stat=i_stat)
          call memocc(i_stat,product(shape(neigh_aveb))*kind(neigh_aveb),'neigh_aveb','allocate_ind_trajectories')
          allocate(ind_indxb(ind_buff),stat=i_stat)
          call memocc(i_stat,product(shape(ind_indxb))*kind(ind_indxb),'ind_indxb','allocate_ind_trajectories')

          ! Start the counter
          bcount_indtraj=1
       endif

    else
       if (do_prn_induced=='Y') then

          i_all=-product(shape(ind_mmomb))*kind(ind_mmomb)
          deallocate(ind_mmomb,stat=i_stat)
          call memocc(i_stat,i_all,'ind_mmomb','allocate_ind_trajectories')
          i_all=-product(shape(neigh_modb))*kind(neigh_modb)
          deallocate(neigh_modb,stat=i_stat)
          call memocc(i_stat,i_all,'neigh_modb','allocate_ind_trajectories')
          i_all=-product(shape(ind_emomb))*kind(ind_emomb)
          deallocate(ind_emomb,stat=i_stat)
          call memocc(i_stat,i_all,'ind_emomb','allocate_ind_trajectories')
          i_all=-product(shape(neigh_aveb))*kind(neigh_aveb)
          deallocate(neigh_aveb,stat=i_stat)
          call memocc(i_stat,i_all,'neigh_aveb','allocate_ind_trajectories')
          i_all=-product(shape(ind_indxb))*kind(ind_indxb)
          deallocate(ind_indxb,stat=i_stat)
          call memocc(i_stat,i_all,'ind_indxb','allocate_ind_trajectories')

       endif

    endif


  end subroutine allocate_ind_trajectories

  !> Buffer all moments
  subroutine buffer_indtraj(mstep,Natom,Mensemble,max_no_neigh,bcount_indtraj,ind_nlistsize,ind_nlist,&
             delta_t,mmom,emom,real_time_measure)

    implicit none

    integer, intent(in) :: mstep     !< Current simulation step
    integer, intent(in) :: Natom     !< Number of atoms in system
    integer, intent(in) :: Mensemble !< Number of ensembles
    integer, intent(in) :: max_no_neigh
    integer, intent(in) :: bcount_indtraj !< Counter of buffer for the induced moments
    integer, dimension(Natom), intent(in):: ind_nlistsize !< Size of neighbour list for Heisenberg exchange couplings
    integer, dimension(max_no_neigh,Natom), intent(in) :: ind_nlist !< Neighbour list for Heisenberg exchange couplings
    real(dblprec), intent(in) :: delta_t  !< Current time step (used for real time measurements)
    real(dblprec), dimension(Natom,Mensemble), intent(in) :: mmom   !< Magnitude of magnetic moments
    real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: emom !< Current unit moment vector
    character(len=1), intent(in) :: real_time_measure   !< Real time measurement flag

    !.. Local variables
    integer :: i,k,neigh
    real(dblprec) :: mod_ave
    real(dblprec), dimension(3) :: ave_mom

    do k=1, Mensemble
       do i=1, Natom
          ind_emomb(1:3,i,bcount_indtraj,k)=emom(1:3,i,k)
          ind_mmomb(i,bcount_indtraj,k)=mmom(i,k)
          ave_mom=0.0d0
          do neigh=1, ind_nlistsize(i)
            ave_mom(1)=ave_mom(1)+emom(1,ind_nlist(neigh,i),k)*mmom(ind_nlist(neigh,i),k)
            ave_mom(2)=ave_mom(2)+emom(2,ind_nlist(neigh,i),k)*mmom(ind_nlist(neigh,i),k)
            ave_mom(3)=ave_mom(3)+emom(3,ind_nlist(neigh,i),k)*mmom(ind_nlist(neigh,i),k)
          enddo
          mod_ave=sqrt(ave_mom(1)**2+ave_mom(2)**2+ave_mom(3)**2)
          neigh_aveb(1:3,i,bcount_indtraj,k)=ave_mom(1:3)/mod_ave
          neigh_modb(i,bcount_indtraj,k)=mod_ave/ind_nlistsize(i)
       end do
    end do

    if (real_time_measure=='Y') then
       ind_indxb(bcount_indtraj)=mstep*delta_t
    else
       ind_indxb(bcount_indtraj)=mstep
    endif

  end subroutine buffer_indtraj

  !> Print all moments to file
  subroutine prn_indtraj(NA,Natom, Nchmax,Mensemble,do_ralloy,Natom_full,atype,atype_ch,ind_mom,&
             simid,real_time_measure)
    !
    !.. Implicit declarations
    implicit none

    integer, intent(in) :: NA
    integer, intent(in) :: Natom     !< Number of atoms in system
    integer, intent(in) :: Nchmax
    integer, intent(in) :: Mensemble !< Number of ensembles
    integer, intent(in) :: do_ralloy
    integer, intent(in) :: Natom_full
    integer, dimension(Natom), intent(in) :: atype !< Type of atom
    integer, dimension(Natom_full), intent(in) :: atype_ch !< Actual site of atom for dilute system
    integer, dimension(NA,Nchmax),intent(in) :: ind_mom
    character(len=8), intent(in) :: simid             !< Name of simulation
    character(len=1), intent(in) :: real_time_measure !< Measurements displayed in real time

    !.. Local variables
    integer :: i, j,ichem
    character(len=30) :: filn

    !.. Executable statements
    write (filn,'(''ind_moment.'',a8,''.out'')') simid
    open(ofileno, file=filn, position="append")
    do i=1, bcount_indtraj
       do j=1, Natom
          if (do_ralloy==0) then
            ichem=1
          else
            ichem=atype_ch(j)
          endif
          if (real_time_measure=='Y') then
             write (ofileno,10043) ind_indxb(i), j, ind_mom(atype(j),ichem),ind_emomb(1,j,i,Mensemble), ind_emomb(2,j,i,Mensemble), ind_emomb(3,j,i,Mensemble), &
                   neigh_aveb(1,j,i,Mensemble),neigh_aveb(2,j,i,Mensemble),neigh_aveb(3,j,i,Mensemble),neigh_modb(j,i,Mensemble)
          else
             write (ofileno,10042) int(ind_indxb(i)), j, ind_mom(atype(j),ichem),ind_emomb(1,j,i,Mensemble), ind_emomb(2,j,i,Mensemble), ind_emomb(3,j,i,Mensemble), &
                   neigh_aveb(1,j,i,Mensemble),neigh_aveb(2,j,i,Mensemble),neigh_aveb(3,j,i,Mensemble),neigh_modb(j,i,Mensemble)
          endif
       end do
    end do
    close(ofileno)
    return
    write (*,*) "Error writing the induced moment trajectories file"
10042 format (i8,2x,i8,2x,i8,2x,2x, es16.8E3,2x,es16.8E3,2x,es16.8E3,2x,es16.8E3,2x,es16.8E3,2x,es16.8E3,2x,es16.8E3)
10043 format (es16.4,2x,i8,2x,i8,2x,2x, es16.8E3,es16.8E3,es16.8E3,es16.8E3,2x,es16.8E3,2x,es16.8E3,2x,es16.8E3)

  end subroutine prn_indtraj

end module prn_induced_info
