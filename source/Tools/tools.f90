!> @copyright
!! Copyright (C) 2008-2018 UppASD group
!! This file is distributed under the terms of the
!! GNU General Public License. 
!! See http://www.gnu.org/copyleft/gpl.txt
module tools

contains
   !> Clean-up routine to deallocate left over arrays
   subroutine deallocate_rest()
      use Parameters
      use Profiling
      !
      use SystemData
      use Measurements
      use MomentData
      use Profiling
      use InputData
      use AutoCorrelation
      use Temperature
      use Gradients
      use prn_trajectories, only : traj_atom, traj_step, traj_buff, ntraj
      !
      implicit none
      !
      !
      integer :: i_all,i_stat
      !
      !
      i_all=-product(shape(coord))*kind(coord)
      deallocate(coord,stat=i_stat)
      call memocc(i_stat,i_all,'coord','deallocate_rest')

      i_all=-product(shape(redcoord))*kind(redcoord)
      deallocate(redcoord,stat=i_stat)
      call memocc(i_stat,i_all,'redcoord','deallocate_rest')

      i_all=-product(shape(bas))*kind(bas)
      deallocate(bas,stat=i_stat)
      call memocc(i_stat,i_all,'bas','deallocate_rest')
      !
      if(ntraj>0) then
         i_all=-product(shape(traj_atom))*kind(traj_atom)
         deallocate(traj_atom,stat=i_stat)
         call memocc(i_stat,i_all,'traj_atom','allocate_measurements')
         i_all=-product(shape(traj_step))*kind(traj_step)
         deallocate(traj_step,stat=i_stat)
         call memocc(i_stat,i_all,'traj_step','allocate_measurements')
         i_all=-product(shape(traj_buff))*kind(traj_buff)
         deallocate(traj_buff,stat=i_stat)
         call memocc(i_stat,i_all,'traj_buff','allocate_measurements')
      end if
      i_all=-product(shape(emomM))*kind(emomM)
      deallocate(emomM,stat=i_stat)
      call memocc(i_stat,i_all,'emomM','deallocate_rest')
      !
      i_all=-product(shape(emom2))*kind(emom2)
      deallocate(emom2,stat=i_stat)
      call memocc(i_stat,i_all,'emom2','deallocate_rest')
      !
      i_all=-product(shape(emom))*kind(emom)
      deallocate(emom,stat=i_stat)
      call memocc(i_stat,i_all,'emom','deallocate_rest')
      !
      if(allocated(spinwait)) then
         i_all=-product(shape(spinwait))*kind(spinwait)
         deallocate(spinwait,stat=i_stat)
         call memocc(i_stat,i_all,'spinwait','deallocate_rest')
      end if
      !
      if(allocated(spinwaitt)) then
         i_all=-product(shape(spinwaitt))*kind(spinwaitt)
         deallocate(spinwaitt,stat=i_stat)
         call memocc(i_stat,i_all,'spinwaitt','deallocate_rest')
      end if
      !
      if(allocated(ipdelta_t)) then
         i_all=-product(shape(ipdelta_t))*kind(ipdelta_t)
         deallocate(ipdelta_t,stat=i_stat)
         call memocc(i_stat,i_all,'ipdelta_t','deallocate_rest')
      end if
      !
      if(allocated(ipTemp)) then
         i_all=-product(shape(ipTemp))*kind(ipTemp)
         deallocate(ipTemp,stat=i_stat)
         call memocc(i_stat,i_all,'ipTemp','deallocate_rest')
      end if
      !
      if(allocated(ipnstep)) then
         i_all=-product(shape(ipnstep))*kind(ipnstep)
         deallocate(ipnstep,stat=i_stat)
         call memocc(i_stat,i_all,'ipnstep','deallocate_rest')
      end if
      if(allocated(ipmcnstep)) then
         i_all=-product(shape(ipmcnstep))*kind(ipmcnstep)
         deallocate(ipmcnstep,stat=i_stat)
         call memocc(i_stat,i_all,'ipmcnstep','deallocate_rest')
      end if
      if(allocated(iplambda1)) then
         i_all=-product(shape(iplambda1))*kind(iplambda1)
         deallocate(iplambda1,stat=i_stat)
         call memocc(i_stat,i_all,'iplambda1','deallocate_rest')
      end if
      if(allocated(iplambda2)) then
         i_all=-product(shape(iplambda2))*kind(iplambda2)
         deallocate(iplambda2,stat=i_stat)
         call memocc(i_stat,i_all,'iplambda2','deallocate_rest')
      end if

      call deallocate_temp()
      call deallocate_gradient_lists()

   end subroutine deallocate_rest

end module tools
