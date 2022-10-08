!! Reading and writing of files for all formats supported by multiscale.
!> authors
!> Edgar Mendez
!> Nikos  Ntallis
!> Manuel Pereiro

module SpinDynamicsIO
    use Parameters
    use DynamicArray
    use SparseMatrix
    use Formats
implicit none

    integer, parameter :: CLOSED_MODE = 0 !< Indicates that the file is not open
    integer, parameter :: READ_MODE   = 1 !< Indicates that the file is open for reading
    integer, parameter :: WRITE_MODE  = 2 !< Indicates that the file is open for writing

    integer, parameter :: MAX_PATH = 260 !< Maximum length of file names and paths 
    integer, parameter :: TAG_LEN = 6 !< length of type tags


    
    type OutFile
       integer :: unit     !< File unit.
       character(len=MAX_PATH) :: filename !< File name, useful for error reports.
       type(formatter) :: format !< format to use when writing
       integer :: position !< record if binary==.true., line when reading text
       integer :: mode     !< One of CLOSED_MODE, READ_MODE or WRITE_MODE
       integer :: estat    !< Error status for the last IO operation on this file.
    end type OutFile

        
    !! Writes a Fortran variable to a file.
    !> Works only for a file containing a single value.
    !> @param[in] filename Name of the file.
    !> @param[in] format Formatter function 
    !> @param[in] name   Name of the value to put into the file.
    !> @param[in] value  
    interface putFile
       ! rank 0
       procedure putFile_d 
       procedure putFile_r
       procedure putFile_i
       procedure putFile_l
       ! rank 1
       procedure putFile_d1 
       procedure putFile_r1
       procedure putFile_i1
       procedure putFile_l1
       ! rank 2
       procedure putFile_d2 
       procedure putFile_r2
       procedure putFile_i2
       procedure putFile_l2
       !! dynamic arrays
       procedure putFile_dyn_i
       procedure putFile_dyn_d
       !! sparse matrices (given iterators)
       procedure putFile_sp
    end interface putFile

    
    !! Writes a Fortran variable to an OutFile.
    interface put
       ! rank 0
       procedure put_d 
       procedure put_r
       procedure put_i
       procedure put_l
       ! rank 1
       procedure put_d1 
       procedure put_r1
       procedure put_i1
       procedure put_l1
       ! rank 2
       procedure put_d2
       procedure put_r2
       procedure put_i2
       procedure put_l2
       !! dynamic arrays
       procedure put_dyn_i
       procedure put_dyn_d
       !! sparse matrices (given iterators)
       procedure put_sp
    end interface put

    private

    public OutFile, createFile, closeFile, put, putFile
    public testDir, remove, getUnit

contains

  !> Determines whether it is possible to create a file in the given directory
  !! @param[in] dir Directory to test.
  !! @return .true. if it is possible to create and write a file in \p dir, .false. otherwise.
  function testDir(dir) result (canWrite)
    character(len=*),intent(in) :: dir
    logical :: canWrite
    character(len=MAX_PATH) :: filename
    integer :: stat,unit
    filename = dir // "/.probe_directory"
    canWrite = .false.
    unit = getUnit()
    open(unit=unit, iostat=stat, file=filename, status='replace')
    if (stat == 0) then
       close(unit, status='delete')
       canWrite = .true.
    end if
    
  end function testDir

  !> Removes a file
  !! @param[in] filename Path of the file to delete
  !! @return .true. on success, .false. otherwise.
  function remove(filename) result(success)
    character(len=*),intent(in) :: filename
    logical :: success
    integer :: unit, stat
    
    success = .false.
    unit = getUnit()
    open(unit=unit, iostat=stat, file=filename, status='old')
    if (stat == 0) then
       close(unit, status='delete')
       success = .true.
    end if
  end function remove
  
  !> Retrieves an unused file unit
  function getUnit() result(unit)
    integer :: unit
    integer, save :: next_unit = 1
    integer :: iostat
    logical :: opened
    
    do 
       inquire (unit=next_unit, opened=opened, iostat=iostat)
       if (iostat.ne.0) cycle
       if (.not.opened) exit
       next_unit = next_unit+1
    end do
    unit = next_unit    
    return
  end function getUnit
  
  !> Opens a file for writing.
  !! If the file exists, the old data is discarded.
  !! @param[in] filename File to write
  !! @param     format   Formatting function used for this file.
  !! @return An OutFile structure describing the open file. Check the field \p estat for errors.
  function createFile(filename, fileformat) result(file)
    character(len=*), intent(in) :: filename
    type(formatter) :: fileformat
    
    type(OutFile) :: file
    integer :: funit
    procedure(writerIf), pointer :: w
    character(len=len(filename)+5) :: finalName
    w => writer

    funit = getUnit()

    finalName = filename // '.' // trim(fileFormat%fileExtension)
    file%unit     = funit
    file%filename = finalName
    file%mode     = CLOSED_MODE
    file%position = 1
    file%format = fileformat

    open(unit=file%unit, file=file%filename, iostat=file%estat, &
         status='replace', access='stream', form='unformatted')
    
    if (file%estat .eq. 0) then
       file%mode = WRITE_MODE
       call file%format%magic(w)
    end if

  contains
    subroutine writer(bytes)
      implicit none
      byte,dimension(:),intent(in) :: bytes
      call putRaw(bytes,file)
    end subroutine writer
  end function createFile

  !> Closes a file if it is open.
  !! The file is reset.
  !! @param[in,out] file OutFile structure describing the file.
  subroutine closeFile(file)
    type(OutFile), intent(inout) :: file
    if (file%unit .ge. 0) then
       close (file%unit,iostat=file%estat)
    end if
    file%unit = -1
    file%filename = "CLOSED"
    file%mode = CLOSED_MODE
    file%position = 0
  end subroutine closeFile

  
  !> Writes a sequence of bytes to a file.
  !! To be used in put_* functions.
  !! @param[in,out] file OutFile structure describing the file.
  !! @param[in] data Byte array to write.
  subroutine putRaw(data, file)
    implicit none
    byte,dimension(:),intent(in) :: data
    type(OutFile),intent(inout) :: file
    write(file%unit,iostat=file%estat) data
    file%position = file%position + size(data,1)
  end subroutine putRaw

  subroutine put_i(name, v, file)
    character(len=*),intent(in) :: name
    integer,intent(in) :: v
    type(OutFile),intent(inout) :: file

    procedure(writerIf), pointer :: w
    w => writer
    call file%format%valueStart(name,FORMAT_TYPE_SCALAR,(/0,0,0/),w)
    call file%format%valueData(name,FORMAT_TYPE_SCALAR,(/0,0,0/),(/1,0,0/),gen(v),w)
    call file%format%valueEnd(name,FORMAT_TYPE_SCALAR,(/0,0,0/),w)
  contains
    subroutine writer(bytes)
      implicit none
      byte,dimension(:), intent(in) :: bytes
      call putRaw(bytes,file)
    end subroutine writer
  end subroutine put_i
  subroutine put_r(name, v, file)
    character(len=*),intent(in) :: name
    real,intent(in) :: v
    type(OutFile),intent(inout) :: file

    procedure(writerIf), pointer :: w
    w => writer
    call file%format%valueStart(name,FORMAT_TYPE_SCALAR,(/0,0,0/),w)
    call file%format%valueData(name,FORMAT_TYPE_SCALAR,(/0,0,0/),(/1,0,0/),gen(v),w)
    call file%format%valueEnd(name,FORMAT_TYPE_SCALAR,(/0,0,0/),w)
  contains
    subroutine writer(bytes)
      implicit none
      byte,dimension(:), intent(in) :: bytes
      call putRaw(bytes,file)
    end subroutine writer
  end subroutine put_r
  subroutine put_d(name, v, file)
    character(len=*),intent(in) :: name
    real(dblprec),intent(in) :: v
    type(OutFile),intent(inout) :: file
    
    procedure(writerIf), pointer :: w
    w => writer
    call file%format%valueStart(name,FORMAT_TYPE_SCALAR,(/0,0,0/),w)
    call file%format%valueData(name,FORMAT_TYPE_SCALAR,(/0,0,0/),(/1,0,0/),gen(v),w)
    call file%format%valueEnd(name,FORMAT_TYPE_SCALAR,(/0,0,0/),w)
  contains
    subroutine writer(bytes)
      implicit none
      byte,dimension(:), intent(in) :: bytes
      call putRaw(bytes,file)
    end subroutine writer
  end subroutine put_d
  subroutine put_l(name, v, file)
    character(len=*),intent(in) :: name
    logical,intent(in) :: v
    type(OutFile),intent(inout) :: file

    procedure(writerIf), pointer :: w
    w => writer
    call file%format%valueStart(name,FORMAT_TYPE_SCALAR,(/0,0,0/),w)
    call file%format%valueData(name,FORMAT_TYPE_SCALAR,(/0,0,0/),(/1,0,0/),gen(v),w)
    call file%format%valueEnd(name,FORMAT_TYPE_SCALAR,(/0,0,0/),w)
  contains
    subroutine writer(bytes)
      implicit none
      byte,dimension(:), intent(in) :: bytes
      call putRaw(bytes,file)
    end subroutine writer
  end subroutine put_l


  subroutine put_i1(name, v, file)
    character(len=*),intent(in) :: name
    integer,dimension(:),intent(in) :: v
    type(OutFile),intent(inout) :: file
    integer :: i
    integer, dimension(3) :: length
    procedure(writerIf), pointer :: w
    length = 0
    length(1) = ubound(v,1)
    w => writer
    call file%format%valueStart(name,FORMAT_TYPE_ARRAY,length,w)
    do i=1,length(1)
       call file%format%valueData(name,FORMAT_TYPE_ARRAY,&
            length,(/i,0,0/),gen(v(i)),w)
    end do
    call file%format%valueEnd(name,FORMAT_TYPE_ARRAY,length,w)
  contains
    subroutine writer(bytes)
      implicit none
      byte,dimension(:), intent(in) :: bytes
      call putRaw(bytes,file)
    end subroutine writer
  end subroutine put_i1
  subroutine put_r1(name, v, file)
    character(len=*),intent(in) :: name
    real,dimension(:),intent(in) :: v
    type(OutFile),intent(inout) :: file
    integer :: i
    integer, dimension(3) :: length
    procedure(writerIf), pointer :: w
    length = 0
    length(1) = ubound(v,1)
    w => writer
    call file%format%valueStart(name,FORMAT_TYPE_ARRAY,length,w)
    do i=1,length(1)
       call file%format%valueData(name,FORMAT_TYPE_ARRAY,&
            length,(/i,0,0/),gen(v(i)),w)
    end do
    call file%format%valueEnd(name,FORMAT_TYPE_ARRAY,length,w)
  contains
    subroutine writer(bytes)
      implicit none
      byte,dimension(:), intent(in) :: bytes
      call putRaw(bytes,file)
    end subroutine writer
  end subroutine put_r1
  subroutine put_d1(name, v, file)
    character(len=*),intent(in) :: name
    real(dblprec),dimension(:),intent(in) :: v
    type(OutFile),intent(inout) :: file
    integer :: i
    integer, dimension(3) :: length
    procedure(writerIf), pointer :: w
    length = 0
    length(1) = ubound(v,1)
    w => writer
    call file%format%valueStart(name,FORMAT_TYPE_ARRAY,length,w)
    do i=1,length(1)
       call file%format%valueData(name,FORMAT_TYPE_ARRAY,&
            length,(/i,0,0/),gen(v(i)),w)
    end do
    call file%format%valueEnd(name,FORMAT_TYPE_ARRAY,length,w)
  contains
    subroutine writer(bytes)
      implicit none
      byte,dimension(:), intent(in) :: bytes
      call putRaw(bytes,file)
    end subroutine writer
  end subroutine put_d1
  subroutine put_l1(name, v, file)
    character(len=*),intent(in) :: name
    logical,dimension(:),intent(in) :: v
    type(OutFile),intent(inout) :: file
    integer :: i
    integer, dimension(3) :: length
    procedure(writerIf), pointer :: w
    length = 0
    length(1) = ubound(v,1)
    w => writer
    call file%format%valueStart(name,FORMAT_TYPE_ARRAY,length,w)
    do i=1,length(1)
       call file%format%valueData(name,FORMAT_TYPE_ARRAY,&
            length,(/i,0,0/),gen(v(i)),w)
    end do
    call file%format%valueEnd(name,FORMAT_TYPE_ARRAY,length,w)
  contains
    subroutine writer(bytes)
      implicit none
      byte,dimension(:), intent(in) :: bytes
      call putRaw(bytes,file)
    end subroutine writer
  end subroutine put_l1

  subroutine put_i2(name, v, file)
    implicit none
    character(len=*),intent(in) :: name
    integer,dimension(:,:),intent(in) :: v
    type(OutFile),intent(inout) :: file
    integer :: j,i
    integer, dimension(3) :: length
    procedure(writerIf), pointer :: w
    length(1:2) = ubound(v)
    length(3) = 0
    w => writer
    call file%format%valueStart(name,FORMAT_TYPE_ARRAY,length,w)
    do j=1,length(2)
       do i=1,length(1)
          call file%format%valueData(name,FORMAT_TYPE_ARRAY,length,&
               (/i,j,0/),gen(v(i,j)),w)
       end do
    end do
    call file%format%valueEnd(name,FORMAT_TYPE_ARRAY,length,w)
  contains
    subroutine writer(bytes)
      implicit none
      byte,dimension(:), intent(in) :: bytes
      call putRaw(bytes,file)
    end subroutine writer
  end subroutine put_i2
  subroutine put_r2(name, v, file)
    character(len=*),intent(in) :: name
    real,dimension(:,:),intent(in) :: v
    type(OutFile),intent(inout) :: file
    integer :: j,i
    integer, dimension(3) :: length
    procedure(writerIf), pointer :: w
    length(1:2) = ubound(v)
    length(3) = 0
    w => writer
    call file%format%valueStart(name,FORMAT_TYPE_ARRAY,length,w)
    do j=1,length(2)
       do i=1,length(1)
          call file%format%valueData(name,FORMAT_TYPE_ARRAY,length,&
               (/i,j,0/),gen(v(i,j)),w)
       end do
    end do
    call file%format%valueEnd(name,FORMAT_TYPE_ARRAY,length,w)
  contains
    subroutine writer(bytes)
      implicit none
      byte,dimension(:), intent(in) :: bytes
      call putRaw(bytes,file)
    end subroutine writer
  end subroutine put_r2
  subroutine put_d2(name, v, file)
    character(len=*),intent(in) :: name
    real(dblprec),dimension(:,:),intent(in) :: v
    type(OutFile),intent(inout) :: file
    integer :: j,i
    integer, dimension(3) :: length
    procedure(writerIf), pointer :: w
    length(1:2) = ubound(v)
    length(3) = 0
    w => writer
    call file%format%valueStart(name,FORMAT_TYPE_ARRAY,length,w)
    do j=1,length(2)
       do i=1,length(1)
          call file%format%valueData(name,FORMAT_TYPE_ARRAY,length,&
               (/i,j,0/),gen(v(i,j)),w)
       end do
    end do
    call file%format%valueEnd(name,FORMAT_TYPE_ARRAY,length,w)
  contains
    subroutine writer(bytes)
      implicit none
      byte,dimension(:), intent(in) :: bytes
      call putRaw(bytes,file)
    end subroutine writer
  end subroutine put_d2
  subroutine put_l2(name, v, file)
    character(len=*),intent(in) :: name
    logical,dimension(:,:),intent(in) :: v
    type(OutFile),intent(inout) :: file
    integer :: j,i
    integer, dimension(3) :: length
    procedure(writerIf), pointer :: w
    length(1:2) = ubound(v)
    length(3) = 0
    w => writer
    call file%format%valueStart(name,FORMAT_TYPE_ARRAY,length,w)
    do j=1,length(2)
       do i=1,length(1)
          call file%format%valueData(name,FORMAT_TYPE_ARRAY,&
               length,(/i,j,0/),gen(v(i,j)),w)
       end do
    end do
    call file%format%valueEnd(name,FORMAT_TYPE_ARRAY,length,w)
  contains
    subroutine writer(bytes)
      implicit none
      byte,dimension(:), intent(in) :: bytes
      call putRaw(bytes,file)
    end subroutine writer
  end subroutine put_l2


  subroutine put_dyn_i(tag,v,file)
    implicit none
    character(len=*), intent(in) :: tag
    type(DynArrayInt), intent(in) :: v
    type(OutFile),intent(inout) :: file

    call put(tag,v%values(1:v%length), file)
    
  end subroutine put_dyn_i
                
  subroutine put_dyn_d(tag,v,file)
    implicit none
    character(len=*), intent(in) :: tag
    type(DynArrayReal), intent(in) :: v
    type(OutFile),intent(inout) :: file

    call put(tag,v%values(1:v%length), file)
    
  end subroutine put_dyn_d

  subroutine put_sp(name, matrix, file, nrows)
    character(len=*), intent(in) :: name    
    type(SpMatrix), intent(inout) :: matrix
    type(OutFile),intent(inout) :: file
    integer, optional, intent(in) :: nrows

    integer :: rows, max_row_elems, elems
    integer :: row, col, elem, count
    real(dblprec) :: val

    if (matrix%entries%length == 0) then
        call file%format%valueStart(name,FORMAT_TYPE_SPARSE_MATRIX,&
             (/0,0,0/),writer)
        call file%format%valueEnd(name,FORMAT_TYPE_SPARSE_MATRIX,&
             (/0,0,0/),writer)
        return
    end if
    call sortMatrixByRow(matrix)
    rows = matrix%row%values(matrix%row%length)
    elems = matrix%entries%length
    max_row_elems = 1
    count = 1
    do elem=2,elems
       if (matrix%row%values(elem-1) == matrix%row%values(elem)) then
          count = count + 1
       else
          max_row_elems = max(max_row_elems,count)
          count = 1
       end if
    end do
    max_row_elems = max(max_row_elems,count)

    if (present(nrows)) then
       rows = nrows
    end if
    call file%format%valueStart(name,FORMAT_TYPE_SPARSE_MATRIX,&
         (/rows,max_row_elems,elems/),writer)
    do elem=1,elems
       row = matrix%row%values(elem)
       col = matrix%col%values(elem)
       val = matrix%entries%values(elem)
       call file%format%valueData(name,FORMAT_TYPE_SPARSE_MATRIX,&
            (/rows,max_row_elems,elems/),(/row,col,elem/),gen(val),writer)
    end do
    call file%format%valueEnd(name,FORMAT_TYPE_SPARSE_MATRIX,&
         (/rows,max_row_elems,elems/),writer)
  contains
    subroutine writer(bytes)
      implicit none
      byte,dimension(:), intent(in) :: bytes
      call putRaw(bytes,file)
    end subroutine writer
  end subroutine put_sp


  subroutine putFile_d(file, format, name, v)
    implicit none
    character(len=*), intent(in) :: file
    type(formatter), intent(in) :: format
    character(len=*), intent(in) :: name
    real(dblprec), intent(in) :: v

    type(OutFile) :: out_file
    out_file = createFile(file,format)
    call put(name,v,out_file)
    call closeFile(out_file)    
  end subroutine putFile_d
  subroutine putFile_r(file, format, name, v)
    implicit none
    character(len=*), intent(in) :: file
    type(formatter), intent(in) :: format
    character(len=*), intent(in) :: name
    real, intent(in) :: v

    type(OutFile) :: out_file
    out_file = createFile(file,format)
    call put(name,v,out_file)
    call closeFile(out_file)    
  end subroutine putFile_r
  subroutine putFile_i(file, format, name, v)
    implicit none
    character(len=*), intent(in) :: file
    type(formatter), intent(in) :: format
    character(len=*), intent(in) :: name
    integer, intent(in) :: v

    type(OutFile) :: out_file
    out_file = createFile(file,format)
    call put(name,v,out_file)
    call closeFile(out_file)    
  end subroutine putFile_i
  subroutine putFile_l(file, format, name, v)
    implicit none
    character(len=*), intent(in) :: file
    type(formatter), intent(in) :: format
    character(len=*), intent(in) :: name
    logical, intent(in) :: v

    type(OutFile) :: out_file
    out_file = createFile(file,format)
    call put(name,v,out_file)
    call closeFile(out_file)    
  end subroutine putFile_l

  subroutine putFile_d1(file, format, name, v)
    implicit none
    character(len=*), intent(in) :: file
    type(formatter), intent(in) :: format
    character(len=*), intent(in) :: name
    real(dblprec), dimension(:), intent(in) :: v

    type(OutFile) :: out_file
    out_file = createFile(file,format)
    call put(name,v,out_file)
    call closeFile(out_file)    
  end subroutine putFile_d1
  subroutine putFile_r1(file, format, name, v)
    implicit none
    character(len=*), intent(in) :: file
    type(formatter), intent(in) :: format
    character(len=*), intent(in) :: name
    real, dimension(:), intent(in) :: v

    type(OutFile) :: out_file
    out_file = createFile(file,format)
    call put(name,v,out_file)
    call closeFile(out_file)    
  end subroutine putFile_r1
  subroutine putFile_i1(file, format, name, v)
    implicit none
    character(len=*), intent(in) :: file
    type(formatter), intent(in) :: format
    character(len=*), intent(in) :: name
    integer, dimension(:), intent(in) :: v

    type(OutFile) :: out_file
    out_file = createFile(file,format)
    call put(name,v,out_file)
    call closeFile(out_file)    
  end subroutine putFile_i1
  subroutine putFile_l1(file, format, name, v)
    implicit none
    character(len=*), intent(in) :: file
    type(formatter), intent(in) :: format
    character(len=*), intent(in) :: name
    logical, dimension(:), intent(in) :: v

    type(OutFile) :: out_file
    out_file = createFile(file,format)
    call put(name,v,out_file)
    call closeFile(out_file)    
  end subroutine putFile_l1


  subroutine putFile_d2(file, format, name, v)
    implicit none
    character(len=*), intent(in) :: file
    type(formatter), intent(in) :: format
    character(len=*), intent(in) :: name
    real(dblprec), dimension(:,:), intent(in) :: v

    type(OutFile) :: out_file
    out_file = createFile(file,format)
    call put(name,v,out_file)
    call closeFile(out_file)    
  end subroutine putFile_d2
  subroutine putFile_r2(file, format, name, v)
    implicit none
    character(len=*), intent(in) :: file
    type(formatter), intent(in) :: format
    character(len=*), intent(in) :: name
    real, dimension(:,:), intent(in) :: v

    type(OutFile) :: out_file
    out_file = createFile(file,format)
    call put(name,v,out_file)
    call closeFile(out_file)    
  end subroutine putFile_r2
  subroutine putFile_i2(file, format, name, v)
    implicit none
    character(len=*), intent(in) :: file
    type(formatter), intent(in) :: format
    character(len=*), intent(in) :: name
    integer, dimension(:,:), intent(in) :: v

    type(OutFile) :: out_file
    out_file = createFile(file,format)
    call put(name,v,out_file)
    call closeFile(out_file)    
  end subroutine putFile_i2
  subroutine putFile_l2(file, format, name, v)
    implicit none
    character(len=*), intent(in) :: file
    type(formatter), intent(in) :: format
    character(len=*), intent(in) :: name
    logical, dimension(:,:), intent(in) :: v

    type(OutFile) :: out_file
    out_file = createFile(file,format)
    call put(name,v,out_file)
    call closeFile(out_file)    
  end subroutine putFile_l2
  
  subroutine putFile_dyn_i(file, format, name, v)
    implicit none
    character(len=*), intent(in) :: file
    type(formatter), intent(in) :: format
    character(len=*), intent(in) :: name
    type(DynArrayInt), intent(in) :: v

    type(OutFile) :: out_file
    out_file = createFile(file,format)
    call put(name,v,out_file)
    call closeFile(out_file)    
  end subroutine putFile_dyn_i
  subroutine putFile_dyn_d(file, format, name, v)
    implicit none
    character(len=*), intent(in) :: file
    type(formatter), intent(in) :: format
    character(len=*), intent(in) :: name
    type(DynArrayReal), intent(in) :: v

    type(OutFile) :: out_file
    out_file = createFile(file,format)
    call put(name,v,out_file)
    call closeFile(out_file)    
  end subroutine putFile_dyn_d

  subroutine putFile_sp(file, format, name, v, nrows)
    implicit none
    character(len=*), intent(in) :: file
    type(formatter), intent(in) :: format
    character(len=*), intent(in) :: name
    type(SpMatrix), intent(inout) :: v
    integer,intent(in),optional :: nrows
    
    type(OutFile) :: out_file
    out_file = createFile(file,format)
    call put(name,v,out_file,nrows)
    call closeFile(out_file)    
  end subroutine putFile_sp

end module SpinDynamicsIO
