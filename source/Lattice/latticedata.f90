!> Module for storing ionic displacement information
!> @author
!> Johan Hellsvik
module LatticeData
   use Parameters
   use Profiling
   !
   implicit none
   !
   real(dblprec), dimension(:,:), allocatable :: mion   !< Ionic mass
   real(dblprec), dimension(:,:), allocatable :: mioninv   !< Inverse ionic mass
   !
   real(dblprec), dimension(:,:,:), allocatable :: uvec   !< Current ionic displacement
   real(dblprec), dimension(:,:,:), allocatable :: uvec2  !< Final (or temporary) ionic displacement
   !
   real(dblprec), dimension(:,:,:), allocatable :: vvec   !< Current ionic velocity
   real(dblprec), dimension(:,:,:), allocatable :: vvec2  !< Final (or temporary) ionic velocity
   !
   real(dblprec), dimension(:,:,:), allocatable :: lvec   !< Current local angular momentum
   real(dblprec), dimension(:,:,:), allocatable :: lvec2  !< Final (or temporary) local angular momentum


   public


contains


   !> Allocate ionic displacement arrays
   subroutine allocate_latticedata(Natom, Mensemble, flag)

      implicit none

      integer, intent(in), optional :: Natom !< Number of atoms in system
      integer, intent(in), optional :: Mensemble !< Number of ensembles 
      integer, intent(in) :: flag  !< Allocate or deallocate (1/-1)

      integer :: i_all, i_stat

      if(flag>0) then
         allocate(mion(Natom,Mensemble),stat=i_stat)
         call memocc(i_stat,product(shape(mion))*kind(mion),'mion','allocate_latticedata')
         allocate(mioninv(Natom,Mensemble),stat=i_stat)
         call memocc(i_stat,product(shape(mioninv))*kind(mioninv),'mioninv','allocate_latticedata')
         allocate(uvec(3,Natom,Mensemble),stat=i_stat)
         call memocc(i_stat,product(shape(uvec))*kind(uvec),'uvec','allocate_latticedata')
         allocate(uvec2(3,Natom,Mensemble),stat=i_stat)
         call memocc(i_stat,product(shape(uvec2))*kind(uvec2),'uvec2','allocate_latticedata')
         allocate(vvec(3,Natom,Mensemble),stat=i_stat)
         call memocc(i_stat,product(shape(vvec))*kind(vvec),'vvec','allocate_latticedata')
         allocate(vvec2(3,Natom,Mensemble),stat=i_stat)
         call memocc(i_stat,product(shape(vvec2))*kind(vvec2),'vvec2','allocate_latticedata')
         allocate(lvec(3,Natom,Mensemble),stat=i_stat)
         call memocc(i_stat,product(shape(lvec))*kind(lvec),'lvec','allocate_latticedata')
      else
         i_all=-product(shape(mion))*kind(mion)
         deallocate(mion,stat=i_stat)
         call memocc(i_stat,i_all,'mion','allocate_latticedata')
         i_all=-product(shape(mioninv))*kind(mioninv)
         deallocate(mioninv,stat=i_stat)
         call memocc(i_stat,i_all,'mioninv','allocate_latticedata')
         i_all=-product(shape(uvec))*kind(uvec)
         deallocate(uvec,stat=i_stat)
         call memocc(i_stat,i_all,'uvec','allocate_latticedata')
         i_all=-product(shape(uvec2))*kind(uvec2)
         deallocate(uvec2,stat=i_stat)
         i_all=-product(shape(vvec))*kind(vvec)
         deallocate(vvec,stat=i_stat)
         call memocc(i_stat,i_all,'vvec','allocate_latticedata')
         i_all=-product(shape(vvec2))*kind(vvec2)
         deallocate(vvec2,stat=i_stat)
         call memocc(i_stat,i_all,'vvec2','allocate_latticedata')
         i_all=-product(shape(lvec))*kind(lvec)
         deallocate(lvec,stat=i_stat)
         call memocc(i_stat,i_all,'lvec','allocate_latticedata')
      end if

   end subroutine allocate_latticedata


end module LatticeData


