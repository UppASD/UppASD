module Qvectors
   use Parameters
   use Profiling
  

   !
   implicit none
   !

   character(len=1) :: qpoints      !< Flag for q-point generation (F=file,A=automatic,C=full cell)
   character(LEN=35) :: qfile       !< File name for q-points
   integer :: nq  !< Number of q-points to sample
   real(dblprec), dimension(:), allocatable :: q_weight              !< Weights of the q points for DOS calculations
   real(dblprec), dimension(:), allocatable :: qnumber               !< Distance from BZ Center ! diamag full reads qfile that has a forth column indicating the distance form the first point of in BZ.  
   real(dblprec), dimension(:,:), allocatable :: q                   !< q-points
   real(dblprec), dimension(:,:), allocatable :: qdir                !< q-points expressed in the basis of the reciprocal lattice vectors
   integer,dimension(2) :: qmin  !< Index of smallest wave vector and gamma point
   real(dblprec), dimension(3) :: r_mid   !< Center of cell (used for q-transforms)
   

   public 


contains
   !---------------------------------------------------------------------------------
   ! SUBROUTINE: deallocate
   !> @brief Deallocate arrays for the q-points
   !---------------------------------------------------------------------------------
   subroutine deallocate_q(do_sc)
      !
      implicit none
      !
      character(len=1), intent(in) :: do_sc
      !
      integer :: i_all,i_stat
      !
      if(allocated(q)) then
         i_all=-product(shape(q))*kind(q)
         deallocate(q,stat=i_stat)
         call memocc(i_stat,i_all,'q','deallocate_q')
      end if
      if(allocated(qdir)) then
         i_all=-product(shape(qdir))*kind(qdir)
         deallocate(qdir,stat=i_stat)
         call memocc(i_stat,i_all,'qdir','deallocate_q')
      end if
      if(allocated(qnumber)) then
         i_all=-product(shape(qnumber))*kind(qnumber)
         deallocate(qnumber,stat=i_stat)
         call memocc(i_stat,i_all,'qnumber','deallocate_q')
      end if

      if(do_sc=='C'.or.do_sc=='D') then

         if (allocated(q_weight)) then
            i_all=-product(shape(q_weight))*kind(q_weight)
            deallocate(q_weight,stat=i_stat)
            call memocc(i_stat,i_all,'q_weight','deallocate_q')
         endif

         !if (do_sc=='D') then
         !   i_all=-product(shape(w))*kind(w)
         !   deallocate(w,stat=i_stat)
         !   call memocc(i_stat,i_all,'w','deallocate_q')
         !end if
      end if
   end subroutine deallocate_q

   !---------------------------------------------------------------------------------
   ! SUBROUTINE: read_q
   !> @brief Read q-points from file for \f$ \mathbf{S}\left(\mathbf{r},t\right) \rightarrow \mathbf{S}\left(\mathbf{q},t\right)\f$ transform
   !---------------------------------------------------------------------------------
   subroutine read_q(C1,C2,C3)
      !
      use Sorting, only : qsort
      use math_functions
      !use Diamag_full, only : do_diamg_full
      implicit none
      !
      real(dblprec), dimension(3), intent(in) :: C1 !< First lattice vector
      real(dblprec), dimension(3), intent(in) :: C2 !< Second lattice vector
      real(dblprec), dimension(3), intent(in) :: C3 !< Third lattice vector
      !integer, intent(out) :: nq                     !< Number of q-vectors
      integer :: iq, idum
      integer :: i_stat,i_all
      integer, dimension(:), allocatable :: ia
      real(dblprec), dimension(3) :: b1,r1
      real(dblprec), dimension(3) :: b2,r2
      real(dblprec), dimension(3) :: b3,r3
      real(dblprec), dimension(3) :: q_tmp
      real(dblprec), dimension(:), allocatable :: dq
      real(dblprec) :: cv
      

      

      open(ifileno, file=adjustl(qfile))
      read (ifileno,*) nq
           
      allocate(q(3,nq),stat=i_stat)
      call memocc(i_stat,product(shape(q))*kind(q),'q','read_q')
      allocate(dq(nq),stat=i_stat)
      call memocc(i_stat,product(shape(dq))*kind(dq),'dq','read_q')
      allocate(ia(nq),stat=i_stat)
      call memocc(i_stat,product(shape(ia))*kind(ia),'ia','read_q')
      if (qpoints=='F') then
         do iq=1,nq
            read (ifileno,*) q(1,iq), q(2,iq), q(3,iq)
         enddo
      else if (qpoints=='G') then
         close(ifileno)
         open(ifileno, file='qpoints.dat')
         do iq=1,nq
            read (ifileno,*) idum, q(1,iq), q(2,iq), q(3,iq)
         enddo
         close(ifileno)
      else if (qpoints=='I') then
         allocate(q_weight(nq),stat=i_stat)
         call memocc(i_stat,product(shape(q_weight))*kind(q_weight),'q_weight','read_q')
         do iq=1,nq
            read(ifileno,*) q(1,iq), q(2,iq), q(3,iq), q_weight(iq)
         enddo
      else  ! qpoints==D
         ! Calculate reciprocal lattice vectors
         print '(2x,a)','-- Real-space lattice --'
         print '(2x,3f12.5)',C1
         print '(2x,3f12.5)',C2
         print '(2x,3f12.5)',C3

         ! r1 = C2xC3
         r1=f_cross_product(C2,C3)
         !r1(1)=C2(2)*C3(3)-C2(3)*C3(2)
         !r1(2)=C2(3)*C3(1)-C2(1)*C3(3)
         !r1(3)=C2(1)*C3(2)-C2(2)*C3(1)
         ! r2 = C3xC1
         r2=f_cross_product(C3,C1)
         !r2(1)=C3(2)*C1(3)-C3(3)*C1(2)
         !r2(2)=C3(3)*C1(1)-C3(1)*C1(3)
         !r2(3)=C3(1)*C1(2)-C3(2)*C1(1)
         ! r3 = C1xC2
         r3=f_cross_product(C1,C2)
         !r3(1)=C1(2)*C2(3)-C1(3)*C2(2)
         !r3(2)=C1(3)*C2(1)-C1(1)*C2(3)
         !r3(3)=C1(1)*C2(2)-C1(2)*C2(1)
         ! cell volume C1*(C2xC3)
         cv=f_volume(C1,C2,C3)
         !c1r1=C1(1)*r1(1)+C1(2)*r1(2)+C1(3)*r1(3)
         !c2r2=C2(1)*r2(1)+C2(2)*r2(2)+C2(3)*r2(3)
         !c3r3=C3(1)*r3(1)+C3(2)*r3(2)+C3(3)*r3(3)
         ! b1=(2pi)*r1/(C1*r1)
         b1=r1/cv
         !b1(1)=r1(1)/c1r1
         !b1(2)=r1(2)/c1r1
         !b1(3)=r1(3)/c1r1
         ! b2=(2pi)*r2/(C1*r1)
         b2=r2/cv
         !b2(1)=r2(1)/c2r2
         !b2(2)=r2(2)/c2r2
         !b2(3)=r2(3)/c2r2
         b3=r3/cv
         ! b3=(2pi)*r3/(C1*r1)
         !b3(1)=r3(1)/c3r3
         !b3(2)=r3(2)/c3r3
         !b3(3)=r3(3)/c3r3
         print '(2x,a)','-- Reciprocal lattice --'
         print '(2x,3f12.5)',b1
         print '(2x,3f12.5)',b2
         print '(2x,3f12.5)',b3

         if (qpoints=='B') then
            allocate(q_weight(nq),stat=i_stat)
            call memocc(i_stat,product(shape(q_weight))*kind(q_weight),'q_weight','read_q')
            do iq=1,nq
               read (ifileno,*) q_tmp(1), q_tmp(2) , q_tmp(3), q_weight(iq)
               q(1,iq)=q_tmp(1)*b1(1)+q_tmp(2)*b2(1)+q_tmp(3)*b3(1)
               q(2,iq)=q_tmp(1)*b1(2)+q_tmp(2)*b2(2)+q_tmp(3)*b3(2)
               q(3,iq)=q_tmp(1)*b1(3)+q_tmp(2)*b2(3)+q_tmp(3)*b3(3)
            enddo

         else
            do iq=1,nq
               read (ifileno,*) q_tmp(1), q_tmp(2) , q_tmp(3)
               q(:,iq)=q_tmp(1)*b1+q_tmp(2)*b2+q_tmp(3)*b3

               !!! !q(1,iq)=q_tmp(1)*b1(1)+q_tmp(2)*b1(2)+q_tmp(3)*b1(3)
               !!! !q(2,iq)=q_tmp(1)*b2(1)+q_tmp(2)*b2(2)+q_tmp(3)*b2(3)
               !!! !q(3,iq)=q_tmp(1)*b3(1)+q_tmp(2)*b3(2)+q_tmp(3)*b3(3)
               !!! q(1,iq)=q_tmp(1)*b1(1)+q_tmp(2)*b1(2)+q_tmp(3)*b1(3)
               !!! q(2,iq)=q_tmp(1)*b2(1)+q_tmp(2)*b2(2)+q_tmp(3)*b2(3)
               !!! q(3,iq)=q_tmp(1)*b3(1)+q_tmp(2)*b3(2)+q_tmp(3)*b3(3)
            enddo
         endif
      endif

      do iq=1,nq
         dq(iq)=(q(1,iq)**2+q(2,iq)**2+q(3,iq)**2)
      enddo
      qmin(1)=minloc(dq,1)
      qmin(2)=minloc(dq,1,abs(dq-dq(qmin(1)))> 0.0_dblprec)
      if(qmin(2)==0) qmin(2)=qmin(1)

      close(ifileno)
      open(ofileno,file='qpoints.out',status='replace')
      do iq=1,nq
         write(ofileno,'(i6,3f14.6)') iq,q(1,iq),q(2,iq),q(3,iq)
      end do
      close(ofileno)

      i_all=-product(shape(dq))*kind(dq)
      deallocate(dq,stat=i_stat)
      call memocc(i_stat,i_all,'dq','read_q')
      i_all=-product(shape(ia))*kind(ia)
      deallocate(ia,stat=i_stat)
      call memocc(i_stat,i_all,'ia','read_q')

   end subroutine read_q


   !---------------------------------------------------------------------------------
   ! SUBROUTINE: setup_qcoord
   !> @brief Automatic setup of q-point mesh for calculating \f$\mathbf{S}\left(\mathbf{q}\right)\rightarrow\mathbf{S}\left(\mathbf{r}\right)\f$
   !---------------------------------------------------------------------------------
   subroutine setup_qcoord(N1,N2,N3,C1,C2,C3)

      use Sorting, only : qsort
      !
      implicit none
      !
      integer, intent(in) :: N1  !< Number of cell repetitions in x direction
      integer, intent(in) :: N2  !< Number of cell repetitions in y direction
      integer, intent(in) :: N3  !< Number of cell repetitions in z direction
      real(dblprec), dimension(3), intent(in) :: C1 !< First lattice vector
      real(dblprec), dimension(3), intent(in) :: C2 !< Second lattice vector
      real(dblprec), dimension(3), intent(in) :: C3 !< Third lattice vector
      !integer, intent(out) :: nq !< Number of q-vectors

      !or (P)lane 2D extended cell
      !
      integer :: iq,xq,yq,zq
      integer :: i_stat, i_all
      integer,dimension(:), allocatable :: ia
      real(dblprec), dimension(3) :: b1,r1
      real(dblprec), dimension(3) :: b2,r2
      real(dblprec), dimension(3) :: b3,r3
      real(dblprec), dimension(:), allocatable :: dq
      real(dblprec) :: c1r1, c2r2, c3r3

      ! Calculate reciprocal lattice vectors
      ! r1 = C2xC3
      r1(1)=C2(2)*C3(3)-C2(3)*C3(2)
      r1(2)=C2(3)*C3(1)-C2(1)*C3(3)
      r1(3)=C2(1)*C3(2)-C2(2)*C3(1)
      ! r2 = C3xC1
      r2(1)=C3(2)*C1(3)-C3(3)*C1(2)
      r2(2)=C3(3)*C1(1)-C3(1)*C1(3)
      r2(3)=C3(1)*C1(2)-C3(2)*C1(1)
      ! r3 = C1xC2
      r3(1)=C1(2)*C2(3)-C1(3)*C2(2)
      r3(2)=C1(3)*C2(1)-C1(1)*C2(3)
      r3(3)=C1(1)*C2(2)-C1(2)*C2(1)
      ! cell volume C1*(C2xC3)
      c1r1=C1(1)*r1(1)+C1(2)*r1(2)+C1(3)*r1(3)
      c2r2=C2(1)*r2(1)+C2(2)*r2(2)+C2(3)*r2(3)
      c3r3=C3(1)*r3(1)+C3(2)*r3(2)+C3(3)*r3(3)
      ! b1=(2pi)*r1/(C1*r1)
      b1(1)=r1(1)/c1r1
      b1(2)=r1(2)/c1r1
      b1(3)=r1(3)/c1r1
      ! b2=(2pi)*r2/(C1*r1)
      b2(1)=r2(1)/c2r2
      b2(2)=r2(2)/c2r2
      b2(3)=r2(3)/c2r2
      ! b3=(2pi)*r3/(C1*r1)
      b3(1)=r3(1)/c3r3
      b3(2)=r3(2)/c3r3
      b3(3)=r3(3)/c3r3
      !
      if(allocated(q)) then
         i_all=-product(shape(q))*kind(q)
         deallocate(q,stat=i_stat)
         call memocc(i_stat,i_all,'q','setup_qcoord')
      end if

      if(qpoints=='C') then
         nq=0
         do zq=-(N3-1)/2,(N3)/2
            do yq=-(N2-1)/2,(N2)/2
               do xq=-(N1-1)/2,(N1)/2
                  nq=nq+1
               end do
            end do
         end do
         allocate(ia(nq),stat=i_stat)
         call memocc(i_stat,product(shape(ia))*kind(ia),'ia','setup_qcoord')
         allocate(q(3,nq),stat=i_stat)
         call memocc(i_stat,product(shape(q))*kind(q),'q','setup_qcoord')
         allocate(qdir(3,nq),stat=i_stat)
         call memocc(i_stat,product(shape(qdir))*kind(qdir),'qdir','setup_qcoord')
         allocate(dq(nq),stat=i_stat)
         call memocc(i_stat,product(shape(dq))*kind(dq),'dq','setup_qcoord')
         iq=0

         ! In qpointsdir.out are stored the q-points expressed
         ! in the bases of the reciprocal lattice vectors
         open(ofileno,file='qpointsdir.out',status='replace')
         do zq=-(N3-1)/2,(N3)/2
            do yq=-(N2-1)/2,(N2)/2
               do xq=-(N1-1)/2,(N1)/2
                  iq=iq+1
                  q(:,iq)=xq/(1.0_dblprec*N1)*b1+yq/(1.0_dblprec*N2)*b2+zq/(1.0_dblprec*N3)*b3
                  qdir(1,iq)=xq/(1.0_dblprec*N1)
                  qdir(2,iq)=yq/(1.0_dblprec*N2)
                  qdir(3,iq)=zq/(1.0_dblprec*N3)
                  write(ofileno,'(i6,6f14.6)') iq,q(1,iq),q(2,iq),q(3,iq), &
                     qdir(1,iq), qdir(2,iq), qdir(3,iq)
                  !xq/(1.0_dblprec*N1), yq/(1.0_dblprec*N2), zq/(1.0_dblprec*N3)
               end do
            end do
         end do
         do iq=1,nq
            dq(iq)=q(1,iq)**2+q(2,iq)**2+q(3,iq)**2
         enddo
         close(ofileno)
         qmin(1)=minloc(dq,1)
         qmin(2)=minloc(dq,1,abs(dq-dq(qmin(1)))> 0.0_dblprec)

      else if(qpoints=='X') then    ! Like 'C' but for double size of q-grid
         nq=0
         do zq=-(N3),(N3)
            do yq=-(N2),(N2)
               do xq=-(N1),(N1)
                  nq=nq+1
               end do
            end do
         end do

         allocate(ia(nq),stat=i_stat)
         call memocc(i_stat,product(shape(ia))*kind(ia),'ia','setup_qcoord')
         allocate(q(3,nq),stat=i_stat)
         call memocc(i_stat,product(shape(q))*kind(q),'q','setup_qcoord')
         allocate(dq(nq),stat=i_stat)
         call memocc(i_stat,product(shape(dq))*kind(dq),'dq','setup_qcoord')
         iq=0

         do zq=-(N3),(N3)
            do yq=-(N2),(N2)
               do xq=-(N1),(N1)
                  iq=iq+1
                  q(:,iq)=xq/(1.0_dblprec*N1)*b1+yq/(1.0_dblprec*N2)*b2+zq/(1.0_dblprec*N3)*b3
                  dq(iq)=q(1,iq)**2+q(2,iq)**2+q(3,iq)**2
               end do
            end do
         end do
         qmin(1)=minloc(dq,1)
         qmin(2)=minloc(dq,1,abs(dq-dq(qmin(1)))> 0.0_dblprec)

      else if(qpoints=='H') then    ! Like 'C' but for double size of q-grid
         nq=0
         do zq=0,N3-1 !-(N3),(N3)
            do yq=0,N2-1 !-(N2),(N2)
               do xq=0,N1-1 !-(N1),(N1)
                  nq=nq+1
               end do
            end do
         end do

         allocate(ia(nq),stat=i_stat)
         call memocc(i_stat,product(shape(ia))*kind(ia),'ia','setup_qcoord')
         allocate(q(3,nq),stat=i_stat)
         call memocc(i_stat,product(shape(q))*kind(q),'q','setup_qcoord')
         allocate(dq(nq),stat=i_stat)
         call memocc(i_stat,product(shape(dq))*kind(dq),'dq','setup_qcoord')
         iq=0

         do zq=0,N3-1 !-(N3),(N3)
            do yq=0,N2-1 !-(N2),(N2)
               do xq=0,N1-1 !-(N1),(N1)
                  iq=iq+1
                  !q(:,iq)=xq/(2.0_dblprec*N1)*b1+yq/(2.0_dblprec*N2)*b2+zq/(2.0_dblprec*N3)*b3
                  q(:,iq)=xq/(1.0_dblprec*N1)*b1+yq/(1.0_dblprec*N2)*b2+zq/(1.0_dblprec*N3)*b3
                  dq(iq)=q(1,iq)**2+q(2,iq)**2+q(3,iq)**2
               end do
            end do
         end do
         qmin(1)=minloc(dq,1)
         qmin(2)=minloc(dq,1,abs(dq-dq(qmin(1)))> 0.0_dblprec)

      else if(qpoints=='P') then
         nq=(4*N1+1)*(4*N3+1)
         allocate(ia(nq),stat=i_stat)
         call memocc(i_stat,product(shape(ia))*kind(ia),'ia','setup_qcoord')
         allocate(q(3,nq),stat=i_stat)
         call memocc(i_stat,product(shape(q))*kind(q),'q','setup_qcoord')
         allocate(dq(nq),stat=i_stat)
         call memocc(i_stat,product(shape(dq))*kind(dq),'dq','setup_qcoord')
         iq=0
         do zq=-2*N3,2*N3
            yq=0
            do xq=-2*N3,2*N3
               iq=iq+1
               q(:,iq)=xq/(1.0_dblprec*N1)*b1+yq/(1.0_dblprec*N2)*b2+zq/(1.0_dblprec*N3)*b3
               dq(iq)=q(1,iq)**2+q(2,iq)**2+q(3,iq)**2
            end do
         end do
         qmin(1)=minloc(dq,1)
         qmin(2)=minloc(dq,1,abs(dq-dq(qmin(1)))> 0.0_dblprec)

      else if(qpoints=='A') then
         ! Currently G-x-y-G-z-y
         ! where x,y,z are the directions of the reciprocal lattice vectors
         nq=N1+(N2-1)+(N2-1)+(N3-1)+(N2-1)+(N3-1)
         allocate(q(3,nq),stat=i_stat)
         call memocc(i_stat,product(shape(q))*kind(q),'q','setup_qcoord')
         allocate(dq(nq),stat=i_stat)
         call memocc(i_stat,product(shape(dq))*kind(dq),'dq','setup_qcoord')
         q=0.0_dblprec
         iq=0
         xq=0
         yq=0
         zq=0
         qmin(1)=1
         ! G->x
         do xq=0,N1-1
            iq=iq+1
            q(:,iq)=xq/(1.0_dblprec*N1)*b1+yq/(1.0_dblprec*N2)*b2+zq/(1.0_dblprec*N3)*b3
         end do
         ! x-> y
         do yq=1,N2-1
            !       xq=(N1-1)-yq*(N1-1)/(N2-1)
            iq=iq+1
            q(:,iq)=xq/(1.0_dblprec*N1)*b1+yq/(1.0_dblprec*N2)*b2+zq/(1.0_dblprec*N3)*b3
         end do
         ! xy->G
         do yq=N2-2,0,-1
            xq=yq/(N2-1)
            iq=iq+1
            q(:,iq)=xq/(1.0_dblprec*N1)*b1+yq/(1.0_dblprec*N2)*b2+zq/(1.0_dblprec*N3)*b3
         end do
         ! G->y
         do yq=1,N2-1
            iq=iq+1
            q(:,iq)=xq/(1.0_dblprec*N1)*b1+yq/(1.0_dblprec*N2)*b2+zq/(1.0_dblprec*N3)*b3
         end do
         ! y->z
         do zq=0,N3-1
            yq=N2-1-zq*(N2-1)/(N3-1)
            iq=iq+1
            q(:,iq)=xq/(1.0_dblprec*N1)*b1+yq/(1.0_dblprec*N2)*b2+zq/(1.0_dblprec*N3)*b3
         end do
         ! z-G
         do zq=N3-2,1,-1
            iq=iq+1
            q(:,iq)=xq/(1.0_dblprec*N1)*b1+yq/(1.0_dblprec*N2)*b2+zq/(1.0_dblprec*N3)*b3
         end do
      end if
      open(ofileno,file='qpoints.out',status='replace')
      do iq=1,nq
         write(ofileno,'(i6,3f14.6)') iq,q(1,iq),q(2,iq),q(3,iq)
      end do
      close(ofileno)
      if (qpoints=='A') then
         do iq=1,nq
            dq(iq)=q(1,iq)**2+q(2,iq)**2+q(3,iq)**2
         enddo
         qmin(1)=minloc(dq,1)
         qmin(2)=minloc(dq,1,abs(dq-dq(qmin(1)))> 0.0_dblprec)
      endif

      if (qpoints=='P' .or. qpoints=='C') then
         i_all=-product(shape(ia))*kind(ia)
         deallocate(ia,stat=i_stat)
         call memocc(i_stat,i_all,'ia','setup_qcoord')

         i_all=-product(shape(dq))*kind(dq)
         deallocate(dq,stat=i_stat)
         call memocc(i_stat,i_all,'dq','setup_qcoord')
      end if
   end subroutine setup_qcoord

   subroutine qvectors_init()

      implicit none

      qfile        = 'qfile'
      qpoints      = 'F'

   end subroutine qvectors_init

   !---------------------------------------------------------------------------------
   !> @brief
   !> Read input parameters.
   !
   !> @author
   !> Anders Bergman
   !---------------------------------------------------------------------------------
   subroutine read_parameters_qvectors(ifile)

      use FileParser

      implicit none

      ! ... Formal Arguments ...
      integer, intent(in) :: ifile   !< File to read from
      !
      ! ... Local Variables ...
      character(len=50) :: keyword, cache
      integer :: rd_len, i_err, i_errb
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
            ! This is the flags for the S(q,w)
            !> - qpoints
            !! Specify format of qfile for correlation (F=file in cart. coord,
            !! D=file in direct coord, C=full BZ + .... )
         case('qpoints')
            read(ifile,*,iostat=i_err) qpoints
            if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err
            !> - qfile
            !! Name of qfile for correlation
         case('qfile')
            read(ifile,'(a)',iostat=i_err) cache
            if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err
            qfile=trim(adjustl(cache))
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
end subroutine read_parameters_qvectors

end module Qvectors
