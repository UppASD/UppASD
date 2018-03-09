!< Module to implement different temperature profiles in the sample
!@TODO generaliza to all kind of system shapes
!> @copyright
!! Copyright (C) 2008-2018 UppASD group
!! This file is distributed under the terms of the
!! GNU General Public License. 
!! See http://www.gnu.org/copyleft/gpl.txt
module Damping
  use Parameters
  use Profiling
  use ErrorHandling
  use InputData, only : mpdamping1,mpdamping2,mpdampingalloy1,&
      mpdampingalloy2,ipdamping1,ipdamping2,ipdampingalloy1,ipdampingalloy2

  implicit none

  real(dblprec), dimension(:),   allocatable :: lambda1_array
  real(dblprec), dimension(:),   allocatable :: lambda2_array
  real(dblprec), dimension(:,:), allocatable :: iplambda1_array
  real(dblprec), dimension(:,:), allocatable :: iplambda2_array

  public :: lambda1_array, lambda2_array, iplambda1_array, iplambda2_array

contains

  !> Allocate the dampings arrays for the site dependent damping
  subroutine allocate_damping(Natom,ipnphase,flag)

     implicit none

     integer, intent(in) :: flag
     integer, intent(in) :: Natom
     integer, intent(in) :: ipnphase

     !.. Local varaibles
     integer :: i_stat, i_all

     !.. Executable parts
     if (flag==0) then
        ! Allocating the damping parameter for the measurement phase
        allocate(lambda1_array(Natom),stat=i_stat)
        call memocc(i_stat,product(shape(lambda1_array))*kind(lambda1_array),'lambda1_array','allocate_damping')
        allocate(lambda2_array(Natom),stat=i_stat)
        call memocc(i_stat,product(shape(lambda2_array))*kind(lambda2_array),'lambda2_array','allocate_damping')
        ! Allocating the site dependent damping for the initial phase
        if(ipnphase>0) then
           allocate(iplambda1_array(ipnphase,Natom),stat=i_stat)
           call memocc(i_stat,product(shape(iplambda1_array))*kind(iplambda1_array),'iplambda1_array','allocate_damping')
           iplambda1_array=0.0d0
           allocate(iplambda2_array(ipnphase,Natom),stat=i_stat)
           call memocc(i_stat,product(shape(iplambda2_array))*kind(iplambda2_array),'iplambda2_array','allocate_damping')
           iplambda2_array=0.0d0
        end if
     else
        ! Deallocating the damping parameter for the measurement phase
        i_all=-product(shape(lambda1_array))*kind(lambda1_array)
        deallocate(lambda1_array,stat=i_stat)
        call memocc(i_stat,i_all,'lambda1_array','allocate_damping')
        i_all=-product(shape(lambda2_array))*kind(lambda2_array)
        deallocate(lambda2_array,stat=i_stat)
        call memocc(i_stat,i_all,'lambda2_array','allocate_damping')

        ! Allocating the site dependent damping for the initial phase
        if(ipnphase>0) then
           i_all=-product(shape(iplambda1_array))*kind(iplambda1_array)
           deallocate(iplambda1_array,stat=i_stat)
           call memocc(i_stat,i_all,'iplambda1_array','allocate_damping')
           i_all=-product(shape(iplambda2_array))*kind(iplambda2_array)
           deallocate(iplambda2_array,stat=i_stat)
           call memocc(i_stat,i_all,'iplambda2_array','allocate_damping')

        end if
     endif

  end subroutine allocate_damping

  !> Filling up the damping arrays
  subroutine setup_damping(NA,Natom,Natom_full,ip_mode,mode,ipnphase,do_ralloy,asite_ch,achem_ch,do_site_damping,&
             do_site_ip_damping,iplambda1,iplambda2,mplambda1,mplambda2,iplambda1_array,iplambda2_array,lambda1_array,lambda2_array)

     implicit none

     integer, intent(in) :: NA
     integer, intent(in) :: Natom
     integer, intent(in) :: ipnphase
     integer, intent(in) :: do_ralloy
     integer, intent(in) :: Natom_full
     integer, dimension(:), intent(in) :: asite_ch
     integer, dimension(:), intent(in) :: achem_ch

     character(len=1), intent(in) :: mode
     character(len=1), intent(in) :: ip_mode
     character(len=1), intent(in) :: do_site_damping
     character(len=1), intent(in) :: do_site_ip_damping

     real(dblprec), intent(in) :: mplambda1
     real(dblprec), intent(in) :: mplambda2
     real(dblprec), dimension(:), intent(in) :: iplambda1
     real(dblprec), dimension(:), intent(in) :: iplambda2
     real(dblprec), dimension(:), intent(out) :: lambda1_array
     real(dblprec), dimension(:), intent(out) :: lambda2_array
     real(dblprec), dimension(:,:), intent(out) :: iplambda1_array
     real(dblprec), dimension(:,:), intent(out) :: iplambda2_array

     integer :: i,j,k

     ! If the site dependent damping is turned on  for the initial phase
     ! Fill up the corresponding arrays

     if (ip_mode=='S') then
        if (do_site_ip_damping=='Y') then
           do i=1, ipnphase
           ! Check if the random alloy flag is on
              if (do_ralloy==0) then
                 ! If there is no random alloy
                 do k=1, Natom
                   if (k<=NA) then
                      j=k
                   else if (k>NA .AND. mod(k,NA)==0) then
                      j=mod(k,NA)+1
                   else if (k>NA .AND. mod(k,NA)/=0) then
                      j=mod(k,NA)
                   end if
                   iplambda1_array(i,k)=ipdamping1(i,j)
                   iplambda2_array(i,k)=ipdamping2(i,j)
                 end do
              else
                 do j=1, Natom_full
                   iplambda1_array(i,j)=ipdampingalloy1(i,asite_ch(j),achem_ch(j))
                   iplambda2_array(i,j)=ipdampingalloy2(i,asite_ch(j),achem_ch(j))
                 end do
              endif
           enddo
        else
           do i=1, ipnphase
              iplambda1_array(i,:)=iplambda1(i)
           end do
        endif
     ! If there is no site dependent damping
     elseif (ip_mode=='H' .or. ip_mode=='M'.or.ip_mode=='I'.or.ip_mode=='D'.or.ip_mode=='R'.or.ip_mode=='L') then
         do i=1, ipnphase
            iplambda1_array(i,:)=iplambda1(i)
            iplambda2_array(i,:)=iplambda2(1)
         enddo
     elseif (ip_mode=='N') then
     elseif (ip_mode=='O') then
     elseif (ip_mode=='G') then
     elseif (ip_mode=='B') then
     elseif (ip_mode=='Q') then
     else
         call ErrorHandling_ERROR('Unrecognized ip_mode: '//ip_mode)
     endif

     ! Fill up the site dependent arrays for the measurement phase
     ! Mapping the damping parameter over the number of atoms Natom
     if (mode=='S' .or. mode=='R') then
        if (do_site_damping=='Y') then
           ! Check if the system is a random alloy
           if (do_ralloy==0) then
              do i=1,Natom
                 if (i<=NA) then
                   j=i
                 else if (i>NA .AND. mod(i,NA)==0) then
                   j=mod(i,NA)+1
                 else if (i>NA .AND. mod(i,NA)/=0) then
                   j=mod(i,NA)
                 end if
                   lambda1_array(i)=mpdamping1(j)
                   lambda2_array(i)=mpdamping2(j)
              enddo
           else
              do i=1, Natom_full
                 lambda1_array(i)=mpdampingalloy1(asite_ch(i),achem_ch(i))
                 lambda2_array(i)=mpdampingalloy2(asite_ch(i),achem_ch(i))
              enddo
           endif
        else

           lambda1_array=mplambda1
           lambda2_array=mplambda2

        endif
     endif

  end subroutine setup_damping

end module Damping
