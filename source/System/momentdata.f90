!> Module for storing magnetic moment information
!> @copyright
!! Copyright (C) 2008-2018 UppASD group
!! This file is distributed under the terms of the
!! GNU General Public License. 
!! See http://www.gnu.org/copyleft/gpl.txt
module MomentData
  use Parameters
  use Profiling
  !
  implicit none
  ! 
  real(dblprec), dimension(:,:,:), allocatable :: emom   !< Current unit moment vector
  real(dblprec), dimension(:,:), allocatable :: mmom !< Magnitude of magnetic moments
  real(dblprec), dimension(:,:), allocatable :: mmomi !< Inverse of magnitude of magnetic moments
  real(dblprec), dimension(:,:,:), allocatable :: emom2  !< Final (or temporary) unit moment vector
  real(dblprec), dimension(:,:,:), allocatable :: emomM  !< Current magnetic moment vector
  real(dblprec), dimension(:,:), allocatable :: mmom2 !< Temporary value of magnitude of magnetic moments
  real(dblprec), dimension(:,:), allocatable :: mmom0 !< Starting magnitude of magnetic moments
  !
  public


contains


  !> Allocate magnitude of magnetic moments
  subroutine allocate_mmoms(Natom,Mensemble,flag)

    implicit none

    integer, intent(in), optional :: Natom !< Number of atoms in system
    integer, intent(in), optional :: Mensemble !< Number of ensembles 
    integer, intent(in) :: flag  !< Allocate or deallocate (1/-1)

    integer :: i_all, i_stat

    if(flag>0) then
       allocate(mmom(Natom,Mensemble),stat=i_stat)
       call memocc(i_stat,product(shape(mmom))*kind(mmom),'mmom','allocate_mmoms')
       allocate(mmom0(Natom,Mensemble),stat=i_stat)
       call memocc(i_stat,product(shape(mmom0))*kind(mmom0),'mmom0','allocate_mmoms')
       allocate(mmom2(Natom,Mensemble),stat=i_stat)
       call memocc(i_stat,product(shape(mmom2))*kind(mmom2),'mmom2','allocate_mmoms')
       allocate(mmomi(Natom,Mensemble),stat=i_stat)
       call memocc(i_stat,product(shape(mmomi))*kind(mmomi),'mmomi','allocate_mmoms')
    else
       i_all=-product(shape(mmom))*kind(mmom)
       deallocate(mmom,stat=i_stat)
       call memocc(i_stat,i_all,'mmom','allocate_moms')
       i_all=-product(shape(mmom0))*kind(mmom0)
       deallocate(mmom0,stat=i_stat)
       call memocc(i_stat,i_all,'mmom0','allocate_moms')
       i_all=-product(shape(mmom2))*kind(mmom2)
       deallocate(mmom2,stat=i_stat)
       call memocc(i_stat,i_all,'mmom2','allocate_moms')
       i_all=-product(shape(mmomi))*kind(mmomi)
       deallocate(mmomi,stat=i_stat)
       call memocc(i_stat,i_all,'mmomi','allocate_moms')
    end if

  end subroutine allocate_mmoms


  !> Allocate directions of magnetic moments
  subroutine allocate_emoms(Natom,Mensemble,flag)

    implicit none

    integer, intent(in) :: Natom !< Number of atoms in system
    integer, intent(in) :: Mensemble !< Number of ensembles 
    integer, intent(in) :: flag  !< Allocate or deallocate (1/-1)

    integer :: i_all, i_stat

    if(flag>0) then
       allocate(emom(3,Natom,Mensemble),stat=i_stat)
       call memocc(i_stat,product(shape(emom))*kind(emom),'emom','allocate_emoms')
       allocate(emom2(3,Natom,Mensemble),stat=i_stat)
       call memocc(i_stat,product(shape(emom2))*kind(emom2),'emom2','allocate_emoms')
       allocate(emomM(3,Natom,Mensemble),stat=i_stat)
       call memocc(i_stat,product(shape(emomM))*kind(emomM),'emomM','allocate_emoms')
    else
       i_all=-product(shape(emom))*kind(emom)
       deallocate(emom,stat=i_stat)
       call memocc(i_stat,i_all,'emom','allocate_eoms')
       i_all=-product(shape(emomM))*kind(emomM)
       deallocate(emomM,stat=i_stat)
       call memocc(i_stat,i_all,'emomM','allocate_eoms')
       i_all=-product(shape(emom2))*kind(emom2)
       deallocate(emom2,stat=i_stat)
       call memocc(i_stat,i_all,'emom2','allocate_eoms')
    end if
  end subroutine allocate_emoms


end module MomentData


