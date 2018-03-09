!> Routines for performing autocorrelation sampling
!! @todo Consider averaging over multiple ensembles
!> @author
!! Johan Hellsvik, Anders Bergman
!> @copyright
!! Copyright (C) 2008-2018 UppASD group
!! This file is distributed under the terms of the
!! GNU General Public License.
!! See http://www.gnu.org/copyleft/gpl.txt
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
  !From simulation
  integer :: n0spinwait
  real(dblprec), dimension(:,:,:), allocatable :: spinwait !< Data for autocorrelation analysis
  !
  public


contains


  !> Allocates the array spinwait that hold sampled moments for different waiting times
  subroutine allocate_autocorr(Natom)

    !.. Implicit declarations
    implicit none

    integer, intent(in) :: Natom !< Number of atoms in system

    !.. Local scalars
    integer :: i_stat

    ! Allocate spin_ arrays
    allocate(spinwait(3,Natom,nspinwait),stat=i_stat)
    call memocc(i_stat,product(shape(spinwait))*kind(spinwait),'spinwait','autocorr_allocate')
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
       write (filn,'(''spin.'',a8,''.'',i1,''.out'')') simid, j
       open(ofileno, file=filn, status="new")
       do i=1,Natom
          spinwait(1:3,i,1) = emom(1:3,i,1)
          write(ofileno,10001) emom(1:3,i,1)
          do j=2,nspinwait
             spinwait(1,i,j)=0d0
             spinwait(2,i,j)=0d0
             spinwait(3,i,j)=0d0
          end do
       enddo
       close(ofileno)
       n0spinwait = 1
    else ! if(initmag==4) then
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
             write (filn,'(''spin.'',a8,''.'',i1,''.out'')') simid, j
          else
             write (filn,'(''spin.'',a8,''.'',i2,''.out'')') simid, j
          end if
          open(ofileno, file=filn, status="old")
          do i=1,Natom
             read(ofileno,10001) spinwait(1:3,i,j)
          enddo
          close(ofileno)
       end do
       do j=n0spinwait+1,nspinwait
          do i=1,Natom
             spinwait(1,i,j)=0d0
             spinwait(2,i,j)=0d0
             spinwait(3,i,j)=0d0
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
             write (filn,'(''spin.'',a8,''.'',i1,''.out'')') simid, j
          else
             write (filn,'(''spin.'',a8,''.'',i2,''.out'')') simid, j
          end if

          open(ofileno, file=filn, status="new")
          do i=1,Natom
             spinwait(1:3,i,j) = emom(1:3,i,1)
             write(ofileno,10001) emom(1:3,i,1)
          enddo
          close(ofileno)
       end if
    end do
10001 format (3es16.8)
  end subroutine autocorr_sample


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
    read (ifileno,'(10x,i10)') nspinwait

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
