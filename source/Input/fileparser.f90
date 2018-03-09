!> Routines for parsing input files
!> @copyright
!! Copyright (C) 2008-2018 UppASD group
!! This file is distributed under the terms of the
!! GNU General Public License. 
!! See http://www.gnu.org/copyleft/gpl.txt
module fileparser
  use Parameters

  implicit none


contains


  !> Routine for reading and isolating keywords
  subroutine bytereader(keyword,rd_len,ifile,i_err)
    ! 
    implicit none
    !
    ! ... Formal Arguments ...
    character(len=*), intent(out) :: keyword  !< Parsed keyword
    integer, intent(out) :: rd_len !< Length of parsed keyword
    integer, intent(in) :: ifile  !< File to read from
    integer, intent(out) :: i_err  !< Error status of reading
    !
    ! ... Local Variables ...
    logical :: rd_done,rd_start
    !
    rd_done=.false.
    rd_start=.false.
    rd_len=0
    do while(.not.rd_done.and.rd_len<len(keyword))
       rd_len=rd_len+1
       read(ifile,'(a1)',advance='no',end=20,eor=10) keyword(rd_len:rd_len)
       rd_start=rd_start.or.(keyword(rd_len:rd_len)/=" ")
       rd_done=rd_start.and.(keyword(rd_len:rd_len)==" ")
    end do
    ! happy ending
    i_err=0
    keyword=adjustl(keyword(1:rd_len)//'')
    return
    ! final word
10  continue
    i_err=10
    keyword=adjustl(keyword(1:rd_len)//'')
    !        print *,'<<<<',trim(keyword),'>>>>>'
    return
    ! end of file
20  continue
    i_err=20
    keyword=adjustl(keyword(1:rd_len)//'')
    return
    !
  end subroutine bytereader


  !> Convert lower case characters to upper case
  subroutine small2caps(str)
    implicit none
    character(len=*),intent(inout):: str  !< string to convert
    !
    integer i
    !
    do i=1,len(str)
       if(str(i:i)>="a" .and. str(i:i)<= "z") str(i:i)=achar(iachar(str(i:i))-32)
    end do
    !
  end subroutine small2caps


  !> Convert upper case characters to lower case
  subroutine caps2small(str)
    implicit none
    character(len=*),intent(inout):: str  !< string to convert
    !
    integer i
    !
    do i=1,len(str)
       if(str(i:i)>="A" .and. str(i:i)<= "Z") str(i:i)=achar(iachar(str(i:i))+32)
    end do
    !
  end subroutine caps2small


end module fileparser
