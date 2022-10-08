!! Data formatting.
!> authors
!> Edgar Mendez
!> Nikos  Ntallis
!> Manuel Pereiro

module Formats
  use Parameters
  implicit none

  !! Identifies a file as a valid binary
  byte,dimension(4),parameter :: BIN_MAGIC = (/77_1,85_1,80_1,0_1/) !reads 'MUP\0'
  character(len=*), parameter :: MATLAB_HEADER = '% Multiscale auto-generated matlab file.' // new_line('a')
  
  !! Size in bytes of the genType structure.
  integer, parameter :: GEN_TYPE_SZ = 16
  !! Identifiers for types accepted by genType
  integer, parameter :: TYPE_NONE    = 0
  integer, parameter :: TYPE_INTEGER = 2
  integer, parameter :: TYPE_REAL_SP = 4
  integer, parameter :: TYPE_REAL_DP = 5
  integer, parameter :: TYPE_LOGICAL = 6
  integer, parameter :: TYPE_STR     = 7

  !! Maximum length for a text or binary sequence representing a number
  integer, parameter :: NUMBER_BUFFER_LEN = 31  

  !! Format strings for plain text:
  !> Consider the output is trimmed, and MUST fit NUMBER_BUFFER_LEN-1
  character(len=*),parameter :: FORMAT_DBLPREC = '(ES27.20)'
  character(len=*),parameter :: FORMAT_REAL    = '(ES27.20)'
  character(len=*),parameter :: FORMAT_INTEGER = '(I25)'
  character(len=*),parameter :: FORMAT_LOGICAL = '(L1)'

  !! Format value types
  integer, parameter :: FORMAT_TYPE_SCALAR        = 0
  integer, parameter :: FORMAT_TYPE_ARRAY         = 1
  integer, parameter :: FORMAT_TYPE_SPARSE_MATRIX = 2

  !! Generic type, useful to have a single formatter for many types
  type genType
     integer :: type
     byte,dimension(GEN_TYPE_SZ) :: data
  end type genType

  !! Interface for any function that dumps bytes somewhere.  
  abstract interface
     subroutine writerIf(bytes)
       implicit none
       byte,dimension(:),intent(in) :: bytes
     end subroutine writerIf
  end interface

  !! Structure that knows how to format data for a specific file type.
  !> magic is called when the file is created, to add a magic number.
  !> valueStart is called before adding data from an scalar, matrix or array
  !> valueData is called for every entry of the matrix or array, or for the value of the scalar
  !> valueEnd is called when no more entries are to be written.
  !> Keep in mind that these three functions can be called more than once per file,
  !> but always all three for each variable on the file and always in that order.
  type formatter
     procedure(writeMagicIf), pointer, nopass :: magic
     procedure(writeValueStartIf), pointer, nopass :: valueStart
     procedure(writeValueIf), pointer, nopass :: valueData
     procedure(writeValueEndIf), pointer, nopass :: valueEnd
     character(len=3) :: fileExtension
  end type formatter

  !! Interfaces for formatter
  abstract interface
     !! Write a magic number for the format, in case it is needed.
     !> @param writer
     subroutine writeMagicIf(writer)
       import :: writerIf
       procedure(writerIf) :: writer
     end subroutine writeMagicIf
     !! In case is required for the format, write data indicating the beginning of
     !! a value.
     !> @param name Name of the value.
     !> @param type Type of the value, one of the FORMAT_TYPE_* parameters.
     !> @param shape Dimensions of the value
     !> @param writer
     subroutine writeValueStartIf(name,type,shape,writer)
       import :: writerIf
       character(len=*), intent(in) :: name
       integer, intent(in) :: type
       integer, dimension(3), intent(in) :: shape
       procedure(writerIf) :: writer
     end subroutine writeValueStartIf
     !! Write the contents of value according to the format
     !> @param name Name of the value.
     !> @param type Type of the value, one of the FORMAT_TYPE_* parameters.
     !> @param shape Dimensions of the value
     !> @param position coordinates of the value written
     !> @param val GenType containing the data to write
     !> @param writer
     subroutine writeValueIf(name,type,shape,position,val,writer)
       import :: writerIf
       import :: genType
       character(len=*), intent(in) :: name
       integer, intent(in) :: type
       integer, dimension(3), intent(in) :: shape, position
       type(genType), intent(in) :: val
       procedure(writerIf) :: writer
     end subroutine writeValueIf
     !! When required by the format, write data indicating the end of a value.
     !> @param name Name of the value.
     !> @param type Type of the value, one of the FORMAT_TYPE_* parameters.
     !> @param shape Dimensions of the value
     !> @param writer
     subroutine writeValueEndIf(name,type,shape,writer)
       import :: writerIf
       character(len=*), intent(in) :: name
       integer, intent(in) :: type
       integer, dimension(3), intent(in) :: shape
       procedure(writerIf) :: writer
     end subroutine writeValueEndIf

     !! Function that builds a formatter structure for a specific format.
     function formatterIf() result (r)
       import :: formatter
       type(formatter) :: r
     end function formatterIf

  end interface


  !! Create a generic type struct holding a value:
  interface gen
     procedure gen_type_none
     procedure gen_type_integer
     procedure gen_type_real_sp
     procedure gen_type_real_dp
     procedure gen_type_logical
     procedure gen_type_str
  end interface gen

  !! Get a value from a generic type struct.
  !> @returns .false. if the type does not match.:
  interface getGen
     procedure  get_gen_type_none
     procedure  get_gen_type_integer
     procedure  get_gen_type_real_sp
     procedure  get_gen_type_real_dp
     procedure  get_gen_type_logical
     procedure  get_gen_type_str
  end interface getGen

  !! Formats a Fortran variable as binary or text.
  interface fformat
     procedure fformat_d
     procedure fformat_r
     procedure fformat_i
     procedure fformat_l
  end interface fformat

  !! Cast anything to a byte array, used because transfer is slow
  interface toByte
     procedure toByte_s
     procedure toByte_d
     procedure toByte_r
     procedure toByte_i
     procedure toByte_l
  end interface toByte
     
  
  !! Given a value or variable, finds the size in bytes required to hold it.
  interface sizeOf
     procedure sizeOf_d
     procedure sizeOf_r
     procedure sizeOf_i
     procedure sizeOf_l
  end interface sizeOf

  private

  public formatGen,fformat,gen,getGen, sizeof
  public writerIf
  public genType

  public formatter
  public plaintextFormat, binaryFormat, matlabFormat

  public FORMAT_TYPE_SCALAR, FORMAT_TYPE_ARRAY, FORMAT_TYPE_SPARSE_MATRIX

contains


  subroutine toByte_s(v, bytearray, length)
    !character(len=*),intent(in) :: v
    character, dimension(:),intent(in) :: v
    byte, dimension(:), intent(inout)::bytearray
    integer, intent(out):: length
    
    length = sizeof(v)
    
  end subroutine toByte_s
  subroutine toByte_d(v, bytearray, length)
    real(dblprec),intent(in) :: v
    byte, dimension(:), intent(inout)::bytearray
    integer, intent(out):: length
    
    length = sizeof(v)
    
  end subroutine toByte_d
  subroutine toByte_r(v, bytearray, length)
    real,intent(in) :: v
    byte, dimension(:), intent(inout)::bytearray
    integer, intent(out):: length
    
    length = sizeof(v)
    
  end subroutine toByte_r
  subroutine toByte_i(v, bytearray, length)
    integer,intent(in) :: v
    byte, dimension(:), intent(inout)::bytearray
    integer, intent(out):: length
    
    length = sizeof(v)
    
  end subroutine toByte_i
  subroutine toByte_l(v, bytearray, length)
    logical,intent(in) :: v
    byte, dimension(:), intent(inout)::bytearray
    integer, intent(out):: length
    
    length = sizeof(v)
    
  end subroutine toByte_l
    
    
    
  
  subroutine formatGen(gt,binary,writer)
    implicit none
    type(genType), intent(in) :: gt
    logical, intent(in) :: binary
    procedure(writerIf) :: writer

    integer :: int
    real :: rsp
    real(dblprec) :: rdp
    logical :: bool
    character(len=GEN_TYPE_SZ) :: str

    integer :: stat
    byte,dimension(NUMBER_BUFFER_LEN) :: buffer
    integer :: length
    equivalence(int,rsp,rdp,bool,str)

    if(getGen(gt)) then
       if(.not. binary) then
          length = 6
          buffer(1:length) = transfer('(NONE)',buffer)
       else
          stat = 0
          return
       end if
    elseif (getGen(gt,int)) then
       call fformat(int,binary,buffer,length,stat)
    elseif (getGen(gt,rsp)) then
       call fformat(rsp,binary,buffer,length,stat)
    elseif (getGen(gt,rdp)) then
       call fformat(rdp,binary,buffer,length,stat)
    elseif (getGen(gt,bool)) then
       call fformat(bool,binary,buffer,length,stat)
    elseif (getGen(gt,str)) then
       length = min(ubound(buffer,1),len(trim(str)))    
       call strToBytes(str,buffer,length)
       stat = 0
    else
       print *,"Attempt to format unkonw type: ",gt%type
       stat = -1
       return
    end if
    call writer(buffer(1:length))
  end subroutine formatGen

  ! Transfer is slow in comparison.
  ! Using comprehensions causes mallocs/frees for some reason
  subroutine strToBytes(str, buffer, length)
    character(len=*),intent(in) :: str
    byte,dimension(:),intent(inout) :: buffer
    integer, intent(in) :: length
    integer :: i
    do i=1,length
       buffer(i) = ichar(str(i:i))
    end do
  end subroutine strToBytes
  
  subroutine fformat_d(v, binary, raw, length, stat)
    implicit none
    real(dblprec), intent(in) :: v
    logical, intent(in) :: binary
    byte,dimension(:), intent(inout) :: raw
    integer, intent(inout) :: length
    integer, intent(out):: stat
    character(len=NUMBER_BUFFER_LEN) :: str
    if (binary) then
       length = sizeOf(v)
       raw = transfer(v,raw)
       stat = 0
    else
       write (str,FORMAT_DBLPREC,iostat=stat) v
       str = adjustl(str)
       length = min(ubound(raw,1),len(trim(str))+1)
       call strToBytes(str,raw,length)
    end if
  end subroutine fformat_d
  subroutine fformat_r(v, binary, raw, length, stat)
    implicit none
    real, intent(in) :: v
    logical, intent(in) :: binary
    byte,dimension(:), intent(inout) :: raw
    integer, intent(inout) :: length
    integer, intent(out):: stat
    character(len=NUMBER_BUFFER_LEN) :: str
    if (binary) then
       length = sizeOf(v)
       raw = transfer(v,raw)
       stat = 0
    else
       write (str,FORMAT_REAL,iostat=stat) v
       str = adjustl(str)
       length = min(ubound(raw,1),len(trim(str))+1)
       call strToBytes(str,raw,length)
    end if
  end subroutine fformat_r
  subroutine fformat_i(v, binary, raw, length, stat)
    integer, intent(in) :: v
    logical, intent(in) :: binary
    byte,dimension(:), intent(inout) :: raw
    integer, intent(inout) :: length
    integer, intent(out):: stat
    character(len=NUMBER_BUFFER_LEN) :: str
    if (binary) then
       length = sizeOf(v)
       raw = transfer(v,raw)
       stat = 0
    else
       write (str,FORMAT_INTEGER,iostat=stat) v
       str = adjustl(str)
       length = min(ubound(raw,1),len(trim(str))+1)
       call strToBytes(str,raw,length)
    end if
  end subroutine fformat_i
  subroutine fformat_l(v, binary, raw, length, stat)
    logical, intent(in) :: v
    logical, intent(in) :: binary
    byte,dimension(:), intent(inout) :: raw
    integer, intent(inout) :: length
    integer, intent(out):: stat
    character(len=NUMBER_BUFFER_LEN) :: str
    if (binary) then
       length = sizeOf(v)
       raw = transfer(v,raw)
       stat = 0
    else
       write (str,FORMAT_LOGICAL,iostat=stat) v
       str = adjustl(str)
       length = 2
       call strToBytes(str,raw,length)
    end if
  end subroutine fformat_l


  function gen_type_none() result(gt)
    implicit none
    type(genType) :: gt
    gt%type = TYPE_NONE
    gt%data = 0
  end function gen_type_none
  function gen_type_integer(v) result(gt)
    implicit none
    integer, intent(in) :: v
    type(genType) :: gt
    gt%type = TYPE_INTEGER
    gt%data = 0
    gt%data = transfer(v,gt%data)
  end function gen_type_integer
  function gen_type_real_sp(v) result(gt)
    implicit none
    real, intent(in) :: v
    type(genType) :: gt
    gt%type = TYPE_REAL_SP
    gt%data = transfer(v,gt%data)
  end function gen_type_real_sp
  function gen_type_real_dp(v) result(gt)
    implicit none
    real(dblprec), intent(in) :: v
    type(genType) :: gt
    gt%type = TYPE_REAL_DP
    gt%data = transfer(v,gt%data)
  end function gen_type_real_dp
  function gen_type_logical(v) result(gt)
    implicit none
    logical, intent(in) :: v
    type(genType) :: gt
    gt%type = TYPE_LOGICAL
    gt%data = transfer(v,gt%data)
  end function gen_type_logical
  function gen_type_str(v) result(gt)
    implicit none
    character(len=*), intent(in) :: v
    type(genType) :: gt
    gt%type = TYPE_STR
    if(len(v) > GEN_TYPE_SZ) then
       stop "String too large for gen_type"
    end if
    gt%data = ichar(' ')
    gt%data = transfer(v,gt%data)
  end function gen_type_str



  !! Unpacking gen type  
  function get_gen_type_none(gt) result(success)
    implicit none
    type(genType), intent(in) :: gt
    logical :: success
    success = .false.
    if(gt%type == TYPE_NONE) then
       success = .true.
    end if
  end function get_gen_type_none
  function get_gen_type_integer(gt, v) result(success)
    implicit none
    integer, intent(out) :: v
    type(genType), intent(in) :: gt
    logical :: success
    success = .false.
    if(gt%type == TYPE_INTEGER) then
       v = transfer(gt%data,v) 
       success = .true.
    end if
  end function get_gen_type_integer
  function get_gen_type_real_sp(gt, v) result(success)
    implicit none
    real, intent(out) :: v
    type(genType), intent(in) :: gt
    logical :: success
    success = .false.
    if(gt%type == TYPE_REAL_SP) then
       v = transfer(gt%data,v) 
       success = .true.
    end if
  end function get_gen_type_real_sp
  function get_gen_type_real_dp(gt, v) result(success)
    implicit none
    real(dblprec), intent(out) :: v
    type(genType), intent(in) :: gt
    logical :: success
    success = .false.
    if(gt%type == TYPE_REAL_DP) then
       v = transfer(gt%data,v) 
       success = .true.
    end if

  end function get_gen_type_real_dp
  function get_gen_type_logical(gt, v) result(success)
    implicit none
    logical, intent(out) :: v
    type(genType), intent(in) :: gt
    logical :: success
    success = .false.
    if(gt%type == TYPE_LOGICAL) then
       v = transfer(gt%data,v) 
       success = .true.
    end if

  end function get_gen_type_logical
  function get_gen_type_str(gt, v) result(success)
    implicit none
    character(len=*), intent(inout) :: v
    type(genType), intent(in) :: gt
    logical :: success
    integer :: length
    success = .false.
    if(gt%type == TYPE_STR) then
       length = min(ubound(gt%data,1),len(v))
       v = " "
       v(1:length) = transfer(gt%data(1:length),v)
       success = .true.
    end if

  end function get_gen_type_str




  !! Given a value, finds the size in bytes required to hold it.
  integer pure function sizeOf_d(r)
    implicit none
    real(dblprec), intent(in) :: r
    real(dblprec),parameter :: r2 = 0
    byte, parameter, dimension(1) :: c = (/0_1/)
    integer, parameter :: sizeOf_d_v = size(transfer(r2,c))
    sizeOf_d = sizeOf_d_v
    return
  end function sizeOf_d
  integer pure function sizeOf_r(r)
    implicit none
    real, intent(in):: r
    real, parameter :: r2 = 0
    byte, parameter, dimension(1) :: c = (/0_1/)
    integer, parameter :: sizeOf_r_v = size(transfer(r2,c))
    sizeOf_r = sizeOf_r_v
    return
  end function sizeOf_r
  integer pure function sizeOf_i(r)
    implicit none
    Integer, intent(in):: r
    Integer, parameter :: r2 = 0
    byte, parameter, dimension(1) :: c = (/0_1/)
    integer, parameter :: sizeOf_i_v = size(transfer(r2,c))
    sizeOf_i = sizeOf_i_v
    return
  end function sizeOf_i
  integer pure function sizeOf_l(r)
    implicit none
    logical, intent(in):: r
    logical, parameter :: r2 = .false.
    byte, parameter, dimension(1) :: c = (/0_1/)
    integer, parameter :: sizeOf_l_v = size(transfer(r2,c))
    sizeOf_l = sizeOf_l_v
    return
  end function sizeOf_l
  integer pure function sizeOf_gt(gt)
    implicit none
    type(genType), intent(in) :: gt

    if (gt%type == TYPE_NONE) then
       sizeOf_gt = 0;
    else if (gt%type == TYPE_INTEGER) then
       sizeOf_gt = sizeOf(1)
    else if (gt%type == TYPE_REAL_SP) then
       sizeOf_gt = sizeOf(1.0)
    else if (gt%type == TYPE_REAL_DP) then
       sizeOf_gt = sizeOf(1.0_dblprec)
    else if (gt%type == TYPE_LOGICAL) then
       sizeOf_gt = sizeOf(.true.)
    else if (gt%type == TYPE_STR) then
       sizeOf_gt = GEN_TYPE_SZ
    end if

    return
  end function sizeOf_gt



  integer function posToIdx(pos,vshape)
    integer, dimension(3), intent(in) :: pos,vshape
    posToIdx = pos(1) + pos(2)*vshape(1) + pos(3)*vshape(1)*vshape(2) - 1
  end function posToIdx

  subroutine dummy_writeMagic(writer)
    procedure(writerIf) :: writer
  end subroutine dummy_writeMagic
  subroutine dummy_writeValueStart(name,type,shape,writer)
    character(len=*), intent(in) :: name
    integer, intent(in) :: type
    integer, dimension(3), intent(in) :: shape
    procedure(writerIf) :: writer
  end subroutine dummy_writeValueStart
  subroutine dummy_writeValueEnd(name,type,shape,writer)
    character(len=*), intent(in) :: name
    integer, intent(in) :: type
    integer, dimension(3), intent(in) :: shape
    procedure(writerIf) :: writer
  end subroutine dummy_writeValueEnd


  subroutine binary_writeMagic(writer)
    procedure(writerIf) :: writer
    call writer(BIN_MAGIC);
  end subroutine binary_writeMagic
  subroutine binary_writeValueStart(name,type,vshape,writer)
    character(len=*), intent(in) :: name
    integer, intent(in) :: type
    integer, dimension(3), intent(in) :: vshape
    procedure(writerIf) :: writer

    integer :: i,rank
    rank=1
    if(vshape(2) .ne. 0) then
       rank=2
       if(vshape(3) .ne. 0) then
          rank=3
       end if
    end if

    call formatGen(gen(type),.true.,writer)
    if (type == FORMAT_TYPE_ARRAY) then
       call formatGen(gen(rank),.true.,writer)
       do i=1,rank
          call formatGen(gen(vshape(i)),.true.,writer)
       end do
    elseif (type == FORMAT_TYPE_SPARSE_MATRIX) then
       do i=1,3
          call formatGen(gen(vshape(i)),.true.,writer)
       end do
    end if
  end subroutine binary_writeValueStart
  subroutine binary_writeValue(name,type,shape,position,val,writer)
    character(len=*), intent(in) :: name
    integer, intent(in) :: type
    integer, dimension(3), intent(in) :: shape, position
    type(genType), intent(in) :: val
    procedure(writerIf) :: writer

    if(type == FORMAT_TYPE_SPARSE_MATRIX) then
       call formatGen(gen(position(1)),.true., writer)
       call formatGen(gen(position(2)),.true., writer)
    end if
    call formatGen(val,.true.,writer)

  end subroutine binary_writeValue
  subroutine binary_writeValueEnd(name,type,shape,writer)
    implicit none
    character(len=*), intent(in) :: name
    integer, intent(in) :: type
    integer, dimension(3), intent(in) :: shape
    procedure(writerIf) :: writer
    if(type == FORMAT_TYPE_SPARSE_MATRIX) then
       call formatGen(gen(-1),.true., writer)
       call formatGen(gen(-1),.true., writer)
    end if
  end subroutine binary_writeValueEnd

  subroutine plaintext_writeValueStart(name,type,vshape,writer)
    character(len=*), intent(in) :: name
    integer, intent(in) :: type
    integer, dimension(3), intent(in) :: vshape
    procedure(writerIf) :: writer
    byte,dimension(1) :: mold
    integer :: i,rank
    rank=1
    if(vshape(2) .ne. 0) then
       rank=2
       if(vshape(3) .ne. 0) then
          rank=3
       end if
    end if

    if(type == FORMAT_TYPE_ARRAY) then
       call writer(transfer("array ",mold))
       do i=1,rank
          call formatGen(gen(vshape(i)),.false.,writer)
       end do
    else if (type == FORMAT_TYPE_SPARSE_MATRIX) then
       call writer(transfer("matrix ",mold))
       do i=1,3
          call formatGen(gen(vshape(i)),.false.,writer)
       end do
    end if

    call formatGen(gen(new_line('a')),.false., writer)

  end subroutine plaintext_writeValueStart
  subroutine plaintext_writeValue(name,type,shape,position,val,writer)
    implicit none
    character(len=*), intent(in) :: name
    integer, intent(in) :: type
    integer, dimension(3), intent(in) :: shape, position
    type(genType), intent(in) :: val
    procedure(writerIf) :: writer

    if(type == FORMAT_TYPE_SPARSE_MATRIX) then
       call formatGen(gen(position(1)),.false., writer)
       call formatGen(gen(position(2)),.false., writer)
    end if

    call formatGen(val,.false.,writer)

    if(type == FORMAT_TYPE_SPARSE_MATRIX .or. & 
         type == FORMAT_TYPE_ARRAY) then
       call formatGen(gen(new_line('a')),.false., writer)
    end if
  end subroutine plaintext_writeValue
  subroutine plaintext_writeValueEnd(name,type,shape,writer)
    character(len=*), intent(in) :: name
    integer, intent(in) :: type
    integer, dimension(3), intent(in) :: shape
    procedure(writerIf) :: writer
    if(type == FORMAT_TYPE_SPARSE_MATRIX) then
       call formatGen(gen("end"),.false., writer)
    end if
    call formatGen(gen(new_line('a')),.false., writer)
  end subroutine plaintext_writeValueEnd


  ! MATLAB
  
  subroutine matlab_writeMagic(writer)
    procedure(writerIf) :: writer
    byte,dimension(1) :: mold
    call writer(transfer(trim(MATLAB_HEADER),mold));
  end subroutine matlab_writeMagic


  subroutine matlab_writeValueStart(name,type,vshape,writer)
    character(len=*), intent(in) :: name
    integer, intent(in) :: type
    integer, dimension(3), intent(in) :: vshape
    procedure(writerIf) :: writer
    byte,dimension(300) :: as_bytes
    character(len=300) :: buffer
    equivalence (buffer,as_bytes)
    
    if (type == FORMAT_TYPE_SPARSE_MATRIX) then
       buffer = name // " = sparse([],[],[]);" // new_line('a')
    elseif (type == FORMAT_TYPE_ARRAY) then
       buffer = name // " = [ "
    else
       buffer = name // " = "
    end if
    call writer(as_bytes(1:len(trim(buffer))))
  end subroutine matlab_writeValueStart

  subroutine matlab_writeValue(name,type,shape,position,val,writer)
    implicit none
    character(len=*), intent(in) :: name
    integer, intent(in) :: type
    integer, dimension(3), intent(in) :: shape, position
    type(genType), intent(in) :: val
    procedure(writerIf) :: writer
    integer :: i,rank
    byte,dimension(300) :: as_bytes
    character(len=300) :: buffer
    equivalence (buffer,as_bytes)

    rank=1
    if(shape(2) .ne. 0) then
       rank=2
       if(shape(3) .ne. 0) then
          rank=3
       end if
    end if

    if(type == FORMAT_TYPE_SPARSE_MATRIX) then
       write(buffer,*) trim(name),'(',position(1),',',position(2),') = '
       call writer(as_bytes(1:len(trim(buffer))))
       call formatGen(val,.false.,writer)
       call formatGen(gen(';' // new_line('a')),.false.,writer)       
    elseif(type == FORMAT_TYPE_ARRAY) then
       buffer = ''
       do i=1,rank-1
          if (position(i) == shape(i)) then
             buffer = trim(buffer) // ";" // new_line('a')
          end if
       end do
       call formatGen(val,.false.,writer)
       call writer(as_bytes(1:len(trim(buffer))))
    else
       call formatGen(val,.false.,writer)
    end if
  end subroutine matlab_writeValue
  subroutine matlab_writeValueEnd(name,type,shape,writer)
    character(len=*), intent(in) :: name
    integer, intent(in) :: type
    integer, dimension(3), intent(in) :: shape
    procedure(writerIf) :: writer
    if(type == FORMAT_TYPE_SPARSE_MATRIX) then
        ! Do nothing
    elseif(type == FORMAT_TYPE_ARRAY) then
       call formatGen(gen("];"),.false., writer)
    else
       call formatGen(gen(";"),.false.,writer)
    end if
    call formatGen(gen(new_line('a')),.false., writer)
  end subroutine matlab_writeValueEnd


  ! Format generators

  function plaintextFormat() result(f)
  implicit none
    type(formatter) :: f
    f%magic => dummy_writeMagic
    f%valueStart => plaintext_writeValueStart
    f%valueData => plaintext_writeValue
    f%valueEnd => plaintext_writeValueEnd
    f%fileExtension = 'dat'
  end function plaintextFormat

  function matlabFormat() result(f)
    type(formatter) :: f
    f%magic => matlab_writeMagic
    f%valueStart => matlab_writeValueStart
    f%valueData => matlab_writeValue
    f%valueEnd => matlab_writeValueEnd
    f%fileExtension = 'm'
  end function matlabFormat


  function binaryFormat() result(f)
    type(formatter) :: f
    f%magic => binary_writeMagic
    f%valueStart => binary_writeValueStart
    f%valueData => binary_writeValue
    f%valueEnd => binary_writeValueEnd
    f%fileExtension = 'dat'
  end function binaryFormat



end module Formats
