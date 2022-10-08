program decodeS
   !
   implicit none
   !
   integer :: nq,nw
   integer :: iq,iw,i,j,k
   !
   real*8, dimension(:,:), allocatable :: sqw,slist,sqw_s
   real*8, dimension(:,:,:), allocatable :: sqw_all
   real*8, dimension(4) :: s_max
   real*8 :: fac1,fac2,t,zres,pi,sigma,sfac,tmp
   integer, dimension(4) :: s_place
   real*8 , dimension(4) :: s_pl_int
   !
   read(*,*) nq,nw
   allocate(sqw(4,nw))
   allocate(sqw_s(4,nw))
   allocate(sqw_all(4,nq,nw))
   allocate(slist(4,nq))
   open(10,file='sqw.dat')
   open(11,file='sqw.norm.dat')
   open(12,file='sqw.list.dat')
   open(14,file='sqw.plot.dat')
   !
   !sigma=10
   sigma=1
!   sigma=5
   zres = 1.0/nw
   pi = 4.0*atan(1.0)
   fac1 = -1/(2*sigma**2)
   fac2 = 1/(sqrt(2*pi)*sigma*nw)
   !
   do iq=1,nq
      s_max=0.0d0
      s_place=0
      sqw=0.0d0
      do iw=1,nw
         read(10,*) i,tmp,tmp,tmp,j,(sqw(k,iw),k=1,4)
      end do
      sqw_s=0.0d0
      do i = 1, nw/2
        do j = 1, nw-1
           t = abs((i-j))
           sfac=exp(fac1*t**2)*fac2
           sqw_s(1,i) = sqw_s(1,i) + sqw(1,j)*sfac
           sqw_s(2,i) = sqw_s(2,i) + sqw(2,j)*sfac
           sqw_s(3,i) = sqw_s(3,i) + sqw(3,j)*sfac
           sqw_s(4,i) = sqw_s(4,i) + sqw(4,j)*sfac
         end do
      end do
      do iw = 1, nw/2
         s_max(1)=max(s_max(1),sqw_s(1,iw))
         if(s_max(1)==sqw_s(1,iw)) s_place(1)=iw
         s_max(2)=max(s_max(2),sqw_s(2,iw))
         if(s_max(2)==sqw_s(2,iw)) s_place(2)=iw
         s_max(3)=max(s_max(3),sqw_s(3,iw))
         if(s_max(3)==sqw_s(3,iw)) s_place(3)=iw
         s_max(4)=max(s_max(4),sqw_s(4,iw))
         if(s_max(4)==sqw_s(4,iw)) s_place(4)=iw
      end do
      if(s_max(1)>0.and.s_max(2)>0.and.s_max(3)>0.and.s_max(4)>0) then
         do iw=1,nw
            sqw_s(1,iw)=sqw_s(1,iw)/s_max(1)
            sqw_s(2,iw)=sqw_s(2,iw)/s_max(2)
            sqw_s(3,iw)=sqw_s(3,iw)/s_max(3)
            sqw_s(4,iw)=sqw_s(4,iw)/s_max(4)
            sqw_all(1,iq,iw)=sqw_s(1,iw)
            sqw_all(2,iq,iw)=sqw_s(2,iw)
            sqw_all(3,iq,iw)=sqw_s(3,iw)
            sqw_all(4,iq,iw)=sqw_s(4,iw)
        end do
      end if
      do iw = 1, nw/2
        write(13,'(2i7,4e14.6)') iq,iw,(sqw(k,iw)**1,k=1,4)
        write(11,'(2i7,4e14.6)') iq,iw,(sqw_s(k,iw)**1,k=1,4)
        
        !write(11,'(2i7,4f12.8)') iq,iw,(sqw_s(k,iw)**1,k=1,4)
!         write(11,'(2i7,4f12.8)') iq,iw,(sqw_s(k,iw)**10,k=1,4)
      end do
      do k=1,4
         iw=s_place(k)
         s_pl_int(k)=iw+0.5*sqw_s(k,iw+1)-0.5*sqw_s(k,iw-1)
      end do
      write(12,'(i7,4f8.2)') iq,s_pl_int(1:4)
      !write(12,'(i7,4i7)') iq,s_place(1:4)
!     write(12,'(i7,4f12.8,4i7)') iq,s_max,s_place(1:4)
   end do
      do iw = 1, nw/2
        write(14,'(i7,400e14.6)') iw,(sqw_all(4,iq,iw)**1+0.2d0*iq,iq=1,nq)
      end do
   close(10)
   close(11)
   close(12)
   close(14)
   deallocate(sqw)
   deallocate(slist)
end program decodeS
