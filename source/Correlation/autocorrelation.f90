!-------------------------------------------------------------------------------
! MODULE: AUTOCORRELATION
!> Routines for performing autocorrelation sampling
!! @todo Consider averaging over multiple ensembles
!> @authors
!> Johan Hellsvik, Anders Bergman
!> @copyright
!> GNU Public License.
!-------------------------------------------------------------------------------
module AutoCorrelation
   use Parameters
   use Profiling
   !
   implicit none
   !
   !From input
   integer :: nspinwait
   integer, dimension(:), allocatable :: spinwaitt
   character(len=35) :: twfile                  !< File name for autocorrelation waiting times
   character(len=1) :: do_autocorr              !< Perform autocorrelation (Y/N)
   ! Printing definitions
   integer :: ac_step !< Interval for sampling average magnetization
   integer :: ac_buff !< Buffer size for average magnetization
   !From simulation
   integer :: n0spinwait
   real(dblprec), dimension(:,:,:), allocatable :: spinwait !< Data for autocorrelation analysis
   !
   real(dblprec), dimension(:,:,:), allocatable   :: autocorr_buff        !< Buffer for average magnetizations
   real(dblprec), dimension(:), allocatable       :: indxb_ac       !< Step counter forautocorrelation
   integer :: bcount_ac    !< Counter of buffer for autocorrelation
   !
   public

contains

   !> Allocates the array spinwait that hold sampled moments for different waiting times
   subroutine allocate_autocorr(Natom,Mensemble)

      !.. Implicit declarations
      implicit none

      integer, intent(in) :: Natom !< Number of atoms in system
      integer, intent(in) :: Mensemble !< Number of ensembles

      !.. Local scalars
      integer :: i_stat

      ! Allocate spin_ arrays

      allocate(spinwait(3,Natom,nspinwait),stat=i_stat)
      call memocc(i_stat,product(shape(spinwait))*kind(spinwait),'spinwait','autocorr_allocate')
      if (do_autocorr=='Y') then
         allocate(autocorr_buff(nspinwait,ac_buff,Mensemble),stat=i_stat)
         call memocc(i_stat,product(shape(autocorr_buff))*kind(autocorr_buff),'autocorr_buff','autocorr_allocate')
         autocorr_buff=0.0_dblprec
         allocate(indxb_ac(ac_buff),stat=i_stat)
         call memocc(i_stat,product(shape(indxb_ac))*kind(indxb_ac),'indxb_ac','autocorr_allocate')
         indxb_ac=0
      end if
      bcount_ac=1

     !   i_all=-product(shape(autocorr_buff))*kind(autocorr_buff)
     !   deallocate(autocorr_buff,stat=i_stat)
     !   call memocc(i_stat,i_all,'autocorr_buff','prn_autocorr')

   end subroutine allocate_autocorr


   !> Initializes the array spinwait for sampling moments during different waiting times
   subroutine autocorr_init(Natom, Mensemble, simid, rstep, initmag, emom)

      !.. Implicit declarations
      implicit none

      integer, intent(in) :: Natom !< Number of atoms in system
      integer, intent(in) :: Mensemble !< Number of ensembles
      character(len=8), intent(in) :: simid !< Name of simulation
      integer, intent(in) :: rstep !< Starting simulation step
      integer, intent(in) :: initmag !< Mode of initialization of magnetic moments (1-4)
      real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: emom   !< Current unit moment vector

      !.. Local scalars
      character(len=30) :: filn
      integer :: i,j

      ! If new simulation, reset arrays
      if(initmag==1.or.initmag==2.or.initmag==3.or.initmag==5.or.initmag==4) then
         j=1
         write (filn,'(''spin.'',a,''.'',i1,''.out'')') trim(simid), j
         !open(ofileno, file=filn, status="new")
         open(ofileno, file=filn)
         do i=1,Natom
            spinwait(1:3,i,1) = emom(1:3,i,1)
            write(ofileno,10001) emom(1:3,i,1)
            do j=2,nspinwait
               spinwait(1,i,j)=0.0_dblprec
               spinwait(2,i,j)=0.0_dblprec
               spinwait(3,i,j)=0.0_dblprec
            end do
         enddo
         close(ofileno)
         n0spinwait = 1
      else
         ! else continue from last waiting time
         i=1
         do
            if(spinwaitt(i).ge.rstep) exit
            if(i.gt.nspinwait) exit
            i=i+1
         end do
         n0spinwait = i-1
         write(*,*) "n0spinwait ", n0spinwait
         do j=1,n0spinwait
            if (j.lt.10) then
               write (filn,'(''spin.'',a,''.'',i1,''.out'')') trim(simid), j
            else
               write (filn,'(''spin.'',a,''.'',i2,''.out'')') trim(simid), j
            end if
            open(ofileno, file=filn, status="old")
            do i=1,Natom
               read(ofileno,10001) spinwait(1:3,i,j)
            enddo
            close(ofileno)
         end do
         do j=n0spinwait+1,nspinwait
            do i=1,Natom
               spinwait(1,i,j)=0.0_dblprec
               spinwait(2,i,j)=0.0_dblprec
               spinwait(3,i,j)=0.0_dblprec
            end do
         end do
      end if
      10001 format (3es16.8)
   end subroutine autocorr_init


   !> Samples the magnetic moments for given waiting times
   subroutine autocorr_sample(Natom, Mensemble, simid, mstep, emom)

      !.. Implicit declarations
      implicit none

      integer, intent(in) :: Natom !< Number of atoms in system
      integer, intent(in) :: Mensemble !< Number of ensembles
      character(len=8), intent(in) :: simid !< Name of simulation
      integer, intent(in) :: mstep !< Current simulation step
      real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: emom   !< Current unit moment vector

      !.. Local scalars
      character(len=30) :: filn
      integer :: i,j

      ! Loop over waiting times
      do j=n0spinwait,nspinwait
         if(mstep.eq.spinwaitt(j)) then
            if (j.lt.10) then
               write (filn,'(''spin.'',a,''.'',i1,''.out'')') trim(simid), j
            else
               write (filn,'(''spin.'',a,''.'',i2,''.out'')') trim(simid), j
            end if

            !open(ofileno, file=filn, status="new")
            open(ofileno, file=filn)
            do i=1,Natom
               spinwait(1:3,i,j) = emom(1:3,i,1)
               write(ofileno,10001) emom(1:3,i,1)
            enddo
            close(ofileno)
         end if
      end do
      10001 format (3es16.8)
   end subroutine autocorr_sample

   !---------------------------------------------------------------------------------
   ! SUBROUTINE: buffer_autocorr
   !> Buffer average magnetization
   !---------------------------------------------------------------------------------
   subroutine buffer_autocorr(Natom,Mensemble,do_autocorr,mstep,nspinwait,spinwait,emom,&
      emomM,bcount_ac,delta_t,real_time_measure)
      !

      !.. Implicit declarations
      implicit none
      !
      integer, intent(in) :: mstep !< Current simulation step
      integer, intent(in) :: Natom !< Number of atoms in system
      integer, intent(in) :: Mensemble   !< Number of ensembles
      integer, intent(in) :: nspinwait   !< Number of waiting times for autocorrelatio
      integer, intent(in) :: bcount_ac !< Counter of buffer for averages
      real(dblprec), intent(in) :: delta_t               !< Current time step (used for real time measurements)
      character(len=1), intent(in) :: do_autocorr        !< Perform autocorrelation (Y/N)
      character(len=1), intent(in) :: real_time_measure  !< Real time measurement flag
      real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: emom   !< Current unit moment vector
      real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: emomM  !< Current magnetic moment vector
      real(dblprec), dimension(:,:,:), intent(in) :: spinwait !< Data for autocorrelation analysis

      !.. Scalar variables
      integer :: i,j,k, i_stat, i_all
      !real(dblprec) :: tm

      !.. Local arrays
      real(dblprec), dimension(3,Mensemble) ::  m
      real(dblprec), dimension(:), allocatable :: autocorr

      !.. Executable statements
      if(do_autocorr=='Y') then
         allocate(autocorr(nspinwait),stat=i_stat)
         call memocc(i_stat,product(shape(autocorr))*kind(autocorr),'autocorr','buffer_autocorr')
         autocorr=0.0_dblprec

         !.. Sum over moments
         do i=1, Natom
            do k=1,Mensemble
               m(:,k) = m(:,k) + emomM(:,i,k)
            end do
            !Autocorr only over sample 1
            do j=1, nspinwait
               autocorr(j) = autocorr(j)+spinwait(1,i,j)*emom(1,i,1)+&
                  spinwait(2,i,j)*emom(2,i,1)+spinwait(3,i,j)*emom(3,i,1)
            end do
         end do

         do i=1, Natom
            do k=1,Mensemble
               m(:,k) = m(:,k) + emomM(:,i,k)
            end do
         end do

         do j=1, nspinwait
            autocorr_buff(j,bcount_ac,1) = autocorr(j)
         end do

         i_all=-product(shape(autocorr))*kind(autocorr)
         deallocate(autocorr,stat=i_stat)
         call memocc(i_stat,i_all,'autocorr','buffer_ac')

         if (real_time_measure=='Y') then
            indxb_ac(bcount_ac)=mstep*delta_t
         else
            indxb_ac(bcount_ac)=mstep
         endif

      end if

   end subroutine buffer_autocorr

   !---------------------------------------------------------------------------------
   ! SUBROUTINE: prn_autocorr
   !> Print autocorrelation data
   !---------------------------------------------------------------------------------
   subroutine prn_autocorr(Natom, simid, do_autocorr, nspinwait,real_time_measure)

      !.. Implicit declarations
      implicit none

      integer, intent(in) :: Natom     !< Number of atoms in system
      integer, intent(in) :: nspinwait !< Number of waiting times for autocorrelation
      character(len=8), intent(in) :: simid       !< Name of simulation
      character(len=1), intent(in) :: do_autocorr !< Perform autocorrelation (Y/N)
      character(len=1), intent(in) :: real_time_measure !< Real time measurement flag
      integer :: i,j, i_all, i_stat
      real(dblprec), dimension(:), allocatable :: autocorr
      character(len=30) :: filn4

      !.. Executable statements
      if(do_autocorr=='Y') then

         allocate(autocorr(nspinwait),stat=i_stat)
         call memocc(i_stat,product(shape(autocorr))*kind(autocorr),'autocorr','prn_autocorr')

         write (filn4,'(''autocorr.'',a,''.out'')') trim(simid)
         open(ofileno, file=filn4, position="append")

         do i=1, bcount_ac
            do j=1, nspinwait
               autocorr(j) = autocorr_buff(j,i,1)/Natom
            end do
            if (real_time_measure=='Y') then
               write (ofileno,10005) indxb_ac(i), autocorr(1:nspinwait)
            else
               write (ofileno,10004) int(indxb_ac(i)), autocorr(1:nspinwait)
            endif
         end do

         i_all=-product(shape(autocorr))*kind(autocorr)
         deallocate(autocorr,stat=i_stat)
         call memocc(i_stat,i_all,'autocorr','prn_autocorr')


         close(ofileno)
      end if

      return

      write(*,*) 'Eror writing the autocorrelation file'
      10004 format (i8,300es16.8)
      10005 format (es16.8,300es16.8)

   end subroutine prn_autocorr

   !> Reads the waiting times for autocorrelation measurements
   subroutine read_tw(twfile)
      !
      !
      implicit none
      !
      character(LEN=30), intent(in) :: twfile !< Name of input file
      !
      integer :: iw
      integer :: i_stat

      open(ifileno, file=adjustl(twfile))
      !read (ifileno,'(10x,i10)') nspinwait
      read (ifileno,*) nspinwait

      allocate(spinwaitt(nspinwait),stat=i_stat)
      call memocc(i_stat,product(shape(spinwaitt))*kind(spinwaitt),'spinwaitt','read_tw')
      do iw=1,nspinwait
         read (ifileno,*) spinwaitt(iw)
      enddo

      close(ifileno)
   end subroutine read_tw

   subroutine read_parameters_autocorr(ifile)
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

            case('do_autocorr')
               read(ifile,*,iostat=i_err) do_autocorr
               if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err

            case('ac_buff')
               read(ifile,*,iostat=i_err) ac_buff
               if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err

            case('ac_step')
               read(ifile,*,iostat=i_err) ac_step
               if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err

            case('acfile')
               read(ifile,'(a)',iostat=i_err) cache
               if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err
               twfile=trim(adjustl(cache))
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
   end subroutine read_parameters_autocorr


end module AutoCorrelation
