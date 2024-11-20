!------------------------------------------------------------------------------------
! MODULE: prn_cudameasurements
!> @brief Data structures and subroutines for printing measurements calculated in CUDA
!> @author Johan Hellsvik
!> @copyright
!> GNU Public License
!-------------------------------------------------------------------------------------
module prn_cudameasurements

   use Parameters
   use ErrorHandling

   implicit none

   ! Printing definitions

   ! Local scalars

   ! Local arrays

   ! Private variables
   private

   ! Public variables

contains

   !---------------------------------------------------------------------------------
   !> @brief
   !> Wrapper for printing measurements calculated in CUDA
   !
   !> @author Johan Hellsvik
   !---------------------------------------------------------------------------------
  subroutine print_observable(simid, Mensemble, obs_name, obs_step, obs_buff, &
       obs_dim, indxb_obs, obs_buffer, obs_label, real_time_measure, delta_t, mstep)

    implicit none

    ! Printing definitions
    character(len=8), intent(in) :: simid !< Simulation name
    integer, intent(in) :: Mensemble !< Number of ensembles
    character(len=16), intent(in) :: obs_name !< Observable name
    integer, intent(in) :: obs_step !< Interval for sampling the observable
    integer, intent(in) :: obs_buff !< Buffer size for storing the observable
    integer, intent(in) :: obs_dim !< Number of columns for the observable
    real(dblprec), dimension(:), allocatable      :: indxb_obs   !< Step counter for the buffer of the observable
    real(dblprec), dimension(:,:,:), allocatable  :: obs_buffer  !< Buffer for the observable
    character(len=16), dimension(:), allocatable  :: obs_label   !< Labels for the components of the observable
    character(len=1), intent(in) :: real_time_measure !< Measurement performed in real time
    real(dblprec), intent(in) :: delta_t !< Current time step
    integer, intent(in) :: mstep !< Current simulation step

    ! Local scalars
    character(len=30) :: filn
    integer :: bcount_obs !< Counter of buffer for observables

    write (filn,'(a,a,''.out'')') obs_name, trim(simid)
    open(ofileno, file=filn, position="append")

    ! Write header to output files for first iteration
    if(abs(indxb_obs(1))==0.0e0_dblprec) then
       ! Averages
       if (real_time_measure=='Y') then
          !write(ofileno,10002)"Time[s]","<M>_x","<M>_y","<M>_z","<M>","M_{stdv}"
          write(ofileno,10002) "Time[s]", obs_label
       else
          !write(ofileno,10002)"#Iter","<M>_x","<M>_y","<M>_z","<M>","M_{stdv}"
          write(ofileno,10002) "#Iter", obs_label
       endif
    end if

    ! Write the contents of the buffer
    do bcount_obs=1, obs_buff
       if (real_time_measure=='Y') then
          indxb_obs(bcount_obs)=mstep*delta_t
       else
          indxb_obs(bcount_obs)=mstep
       endif
       ! Is index for Mensemble needed here?
       write(ofileno,10001) indxb_obs, obs_buffer(2:obs_dim,bcount_obs,1)
    end do

    close(ofileno)
    return

    write (*,*) "Error writing the averages file"

10001 format (i8,12es16.8)
10002 format (a8,12a16)

  end subroutine print_observable

  ! SUBROUTINE: cudameasurements_allocations
  !> Allocation of the necessary arrays for the printing of observables sampled within CUDA code
  !---------------------------------------------------------------------------------------------
  subroutine cudameasurements_allocations(Mensemble, obs_step, obs_buff, obs_dim, &
      indxb_obs, obs_buffer, obs_label, flag)

    implicit none

    integer, intent(in) :: Mensemble  !< Number of ensembles
    integer, intent(in) :: obs_step !< Interval for sampling the observable
    integer, intent(in) :: obs_buff !< Buffer size for storing the observable
    integer, intent(in) :: obs_dim !< Number of columns for the observable
    real(dblprec), dimension(:), allocatable      :: indxb_obs   !< Step counter for the buffer of the observable
    real(dblprec), dimension(:,:,:), allocatable  :: obs_buffer  !< Buffer for the observable
    character(len=16), dimension(:), allocatable  :: obs_label   !< Labels for the components of the observable
    integer, intent(in) :: flag       !< Allocate or deallocate (1/-1)
    character(len=30) :: filn

    !.. Local variables
    integer :: i_stat, i_all

    if(flag>0) then
       ! Allocations for the averages
       !allocate(mavg_buff(20,avrg_buff,Mensemble),stat=i_stat)
       allocate(obs_buffer(obs_dim,obs_buff,Mensemble),stat=i_stat)
       call memocc(i_stat,product(shape(obs_buffer))*kind(obs_buffer),'obs_buffer','cudameasurements_allocations')
       obs_buffer=0.0_dblprec
       allocate(obs_label(20),stat=i_stat)
       call memocc(i_stat,product(shape(obs_label))*kind(obs_label),'obs_label','cudameasurements_allocations')
       obs_buffer=0.0_dblprec
       allocate(indxb_obs(obs_buff),stat=i_stat)
       call memocc(i_stat,product(shape(indxb_obs))*kind(indxb_obs),'indxb_obs','cudameasurements_allocations')
       indxb_obs=0
    else
       ! Deallocations for averages
       i_all=-product(shape(obs_buffer))*kind(obs_buffer)
       deallocate(obs_buffer,stat=i_stat)
       call memocc(i_stat,i_all,'obs_buffer','cudameasurements_allocations')
       i_all=-product(shape(indxb_obs))*kind(indxb_obs)
       deallocate(indxb_obs,stat=i_stat)
       call memocc(i_stat,i_all,'indxb_obs','cudameasurements_allocations')
    end if
  end subroutine cudameasurements_allocations

end module prn_cudameasurements
