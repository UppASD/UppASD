!> authors
!> Edgar Mendez
!> Nikos  Ntallis
!> Manuel Pereiro

module MultiscaleFileParser
    use Parameters
implicit none

    !> Structure describing the parsing state of an open file.
    type FileData
       character(len=128) :: filename
       integer :: fileId
       character(len=256) :: line
       integer :: lineNr
       integer :: startWord, endWord
       character(len=128) :: word
       integer :: ierr
       character(len=1024) :: errMsg
       logical :: atEndOfLine
    end type FileData

    !! Predefined parsing error messages
    character(len=*), parameter :: ERROR_PARSING_INT = "Expected an integer."
    character(len=*), parameter :: ERROR_PARSING_REAL = "Expected a real number."
    character(len=*), parameter :: ERROR_PARSING_LOGICAL = "Expected a logical, " // & 
         "use either T (true) or F (false)."
    character(len=*), parameter :: ERROR_PARSING_CHARACTER = "Expected a character."

    character(len=*), parameter :: DELIMITER = "/" ! Optional delimiter between fields

    interface
       subroutine filePeeker(fData)
         import FileData
         implicit none
         type(FileData), intent(inout) :: fData
       end subroutine filePeeker
    end interface

    
private

public openFile, FileData, isAtEndOfLine, parseCharacter, &
       parseLogical, parseReal, parseInt, countLinesStartWithInteger, &
       readNextLine, createSimpleErrorMsg, copyErrorMsg, &
       createErrorMsg, readNextWord, createErrorLineMsg, &
       skipOptionalDelimiter, peekFile, skipToNextLine
contains

!< Ensure the parser is reading at the end of the line,
!! or generate an error otherwise
subroutine verifyEndOfLine(fData)
implicit none
    type(FileData), intent(inout) :: fData
    
    if (.not. isAtEndOfLine(fData)) &
        call createErrorMsg(fData, 1, 'Some junk at the end of the line.')
end subroutine

!< Copy the error state and error message from one parser to another
subroutine copyErrorMsg(src, dest)
implicit none
    type(FileData), intent(in) :: src
    type(FileData), intent(inout) :: dest
    
    dest%ierr = src%ierr
    dest%errMsg = src%errMsg
end subroutine

!< Puts together an error message
!! @param fData The FileData structure where the error message will be placed
!! @param ierr The error code that the error message is connected with
!! @param msg The error message that should be presented
subroutine createErrorMsg(fData, ierr, msg)
implicit none
    type(FileData), intent(inout) :: fData
    integer, intent(in) :: ierr
    character(len=*), intent(in) :: msg
    
    call createErrorLineMsg(fData)

    write (fData%errMsg, '(a)') trim(fData%errMsg) // NEW_LINE('A') // trim(msg)
    fData%ierr = ierr
end subroutine createErrorMsg

!< Creates an error message that contains the file, line and column where the
!! parser is reading
subroutine createErrorLineMsg(fData)
    type(FileData), intent(inout) :: fData
    
    character(len=10) :: lineNrString
    character(len=10) :: columnNrString
    character(len=512) :: pointerLine
    
    integer :: errorCursor
        
    errorCursor = min(len(trim(fData%line)) + 1, fData%startWord)
    pointerLine = ''
    pointerLine(errorCursor:errorCursor) = '^'

    write(lineNrString, '(i8)') fData%lineNr
    write(columnNrString, '(i8)') errorCursor
    
    write (fData%errMsg, '(4a)') 'Error in file "', trim(fData%filename), '", line: ' // trim(adjustl(lineNrString)) &
        // ', column: ' // trim(adjustl(columnNrString)) // NEW_LINE('A') // &
        trim(adjustl(fData%line)) // NEW_LINE('A') // trim(pointerLine)
end subroutine

!< Create an error message with given error code and text
subroutine createSimpleErrorMsg(fData, ierr, msg)
implicit none
    type(FileData), intent(inout) :: fData
    integer, intent(in) :: ierr
    character(len=*), intent(in) :: msg
    
    write (fData%errMsg, '(4a)') 'Error in file "', trim(fData%filename), '": ', trim(msg)
    fData%ierr = ierr
end subroutine

!< Returns a file identifier that is not being used
integer function getUnusedFileId()
implicit none
    integer, save :: next_unit = 1
    integer :: ierr
    logical :: opened
    
    do
        inquire(unit=next_unit, opened=opened, iostat=ierr)
        if (ierr /= 0) cycle
        if (.not. opened) exit
        next_unit = next_unit + 1
    end do
    getUnusedFileId = next_unit
end function

!< Initiate the parsing data structure from a file name
subroutine openFile(fData, filename)
implicit none
    type(FileData), intent(inout) :: fData
    character(len=*), intent(in) :: filename
    
    fData%fileId = getUnusedFileId()
    fData%filename = filename
    open(unit=fData%fileId, file=fData%filename, status='old', iostat=fData%ierr)
    if (fData%ierr /= 0) then
        fData%errMsg = 'Cant open the file: ' // trim(fData%filename)
        return
    endif
    fData%lineNr = 0
    fData%startWord = 1
    fData%endWord = 0
    fData%line = ''
    fData%word = ''
    fData%ierr = 0
    fData%atEndOfLine = .true.
    call skipToNextLine(fData)
end subroutine openFile

!< Place the parser at the beginning of the giving line
subroutine gotoLine(fData, lineNr)
implicit none
    type(FileData), intent(inout) :: fData
    integer, intent(in) :: lineNr
    
    close(fData%fileId)
    call openFile(fData, fData%filename)
    do while (fData%lineNr < lineNr)
        call skipToNextLine(fData)
    enddo
end subroutine gotoLine

!< Checks if the parser is reading at the end of line
logical pure function isAtEndOfLine(fData)
implicit none
    type(FileData), intent(in) :: fData
    
    isAtEndOfLine = fData%atEndOfLine
end function 

!< Counts lines from the current parser position that starts with an integer
!! This is used to count the number of entries in some lists
!! Todo: This function would be more useful if it used a predicate instead of
!!       being limited to `lines starting with integers'
integer function countLinesStartWithInteger(fData)
implicit none
    type(FileData), intent(inout) :: fData
    integer :: nrOfLines
    integer :: lineStart
    integer :: testRead
    integer :: errConvert

    lineStart = fData%lineNr
    nrOfLines = 0
    do while (.true.)
        read(fData%word, *, iostat=errConvert) testRead
        if (errConvert /= 0) then
            exit
        else
            nrOfLines = nrOfLines + 1
            call skipToNextLine(fData)
            if (fData%ierr /= 0) then
                fData%ierr = 0
                exit
            endif
        endif
    enddo
    call gotoLine(fData, lineStart)
    countLinesStartWithInteger = nrOfLines
end function

!< Apply a function on the current line of the parser.
!! On return the parser is restored to its previous position.
subroutine peekFile(fData, f)
  implicit none
    type(FileData), intent(inout) :: fData
    procedure(filePeeker) :: f
    
    integer :: startWord, endWord, lineNr
    startWord = fData%startWord
    endWord = fData%endWord
    lineNr = fData%lineNr
    
    call f(fData)    
    
    close(fData%fileId)
    call openFile(fData, fData%filename)
    do while (fData%lineNr < lineNr)
        call skipToNextLine(fData)
    enddo
    do while (fData%startWord < startWord)
        call readNextWord(fData)
    enddo

    if(fData%endWord .ne. endWord) then
       print *, "Something went wrong when reading the configuration file."
       print *, trim(fData%filename)
       print *, "Unlikely, but possibly, it was modified while being read."
       fData%ierr = 1
    end if    
    
end subroutine peekFile
  



!< Advance the parser to the next line.
!! If there are unread symbols generates an error message.
subroutine readNextLine(fData)
implicit none
    type(FileData), intent(inout) :: fData
    
    call readNextWord(fData)
    call verifyEndOfLine(fData)
    if (fData%ierr /= 0) return
    call skipToNextLine(fData)
end subroutine readNextLine

!< Makes the parser ignore symbols until it reaches the beginning of the next line.
subroutine skipToNextLine(fData)
    type(FileData), intent(inout) :: fData

    fData%line = ''
    fData%atEndOfLine = .true.
    if (is_iostat_end(fData%ierr)) then
        fData%startWord = 1
        fData%endWord = 0
        fData%atEndOfLine = .true.
        return
    endif
    do while (isAtEndOfLine(fData))
        read(fData%fileId, '(a)', iostat=fData%ierr) fData%line
        fData%lineNr = fData%lineNr + 1
        if (is_iostat_end(fData%ierr)) then
            return
        elseif (fData%ierr == 0) then
            fData%line = trim(adjustl(fData%line))
            fData%startWord = 1
            fData%endWord = 0
            fData%atEndOfLine = .false.
            call readNextWord(fData)
        endif
    enddo
end subroutine skipToNextLine

!< If the next word is the DELIMITER, advances the parser one word
subroutine skipOptionalDelimiter(fData)
implicit none
    type(FileData), intent(inout) :: fData

    if(.not. fData%atEndOfLine .and. fData%ierr == 0) then
       if (fData%word == DELIMITER) then
          call readNextWord (fData)
       end if
    end if
end subroutine skipOptionalDelimiter
  
!< Reads the next word in the line and updates the state of the parser
subroutine readNextWord(fData)
implicit none
    type(FileData), intent(inout) :: fData
    
    character(len=1) :: currentChar
    
    if (is_iostat_end(fData%ierr) .or. isAtEndOfLine(fData)) then
        return
    endif
    fData%startWord = fData%endWord + 1
    do while (fData%startWord <= len(fData%line))
        currentChar = fData%line(fData%startWord:fData%startWord)
        if (isWhitespace(currentChar)) then
            fData%startWord = fData%startWord + 1
        else
            exit
        endif
    enddo
    if (fData%startWord > len(fData%line) .or. currentChar == "#") then
        fData%atEndOfLine = .true.
        fData%word = ''
    else
        fData%atEndOfLine = .false.
    
        fData%endWord = fData%startWord
        do while (fData%endWord < len(fData%line))
            currentChar = fData%line((fData%endWord + 1):(fData%endWord + 1))
            if (isWhitespace(currentChar) .or. currentChar == "#") then
                exit
            else
                fData%endWord = fData%endWord + 1
            endif
        enddo
        
        fData%word = fData%line(fData%startWord:fData%endWord)
    endif
end subroutine readNextWord

!< Checks if c is a whitespace or horizontal tab
logical function isWhitespace(c)
implicit none
    character(len=1), intent(in) :: c
    
    isWhitespace = c == ' ' .or. iachar(c) == 9 !9 is Horizontal tab
end function

!< If the current word in the parser is a valid integer, parses it and writes the
!! value in val.
!! Otherwise it creates an error message.
subroutine parseInt(fData, val)
implicit none
    type(FileData), intent(inout) :: fData
    integer, intent(out) :: val
    ! Sometimes fortran's read does neither fail nor read a value
    !  when the input is incorrect. We detect that with a fixed initial value.
    ! Let's hope no parameter can be realistically HUGE(val)
    val = HUGE(val)
    
    if (fData%ierr /= 0) then
        call createErrorMsg(fData, 1, ERROR_PARSING_INT)
        return
    endif
    read(fData%word, *, iostat=fData%ierr) val
    if (val == HUGE(val) .or. fData%ierr /= 0) call createErrorMsg(fData, 1, ERROR_PARSING_INT)
end subroutine parseInt

!< If the current word in the parser is a valid real, parses it and writes the
!! value in val.
!! Otherwise it creates an error message. 
subroutine parseReal(fData, val)
implicit none
    type(FileData), intent(inout) :: fData
    real(dblprec), intent(out) :: val
    ! Sometimes fortran's read does neither fail nor read a value
    !  when the input is incorrect. We detect that with a fixed initial value.
    ! Let's hope no parameter can be realistically HUGE(val)
    val = HUGE(val)
    if (fData%ierr /= 0) then
        call createErrorMsg(fData, 1, ERROR_PARSING_REAL)
        return
     end if
     read(fData%word, *, iostat=fData%ierr) val
     if (val == HUGE(val) .or. fData%ierr /= 0) then
        call createErrorMsg(fData, 1, ERROR_PARSING_REAL)
     end if
end subroutine parseReal

!< If the current word in the parser is a logical-parseable character
!! ('tT1' or 'fF0'), parses it and writes the value in val.
!! Otherwise it creates an error message.
subroutine parseLogical(fData, val)
implicit none
    type(FileData), intent(inout) :: fData
    logical, intent(out) :: val
    
    if (fData%ierr /= 0) then
        call createErrorMsg(fData, 1, ERROR_PARSING_LOGICAL)
        return
    end if
    if(fData%word == "t" .or. fData%word == "T" .or. fData%word == "1") then
       val = .true.
    elseif (fData%word == "f" .or. fData%word == "F" .or. fData%word == "0") then
       val = .false.
    else
       call createErrorMsg(fData, 1, ERROR_PARSING_LOGICAL)
    end if
end subroutine parseLogical

!< If the current word in the parser is a single character,
!! it reads it copies it into val.
!! Otherwise it creates an error message.
subroutine parseCharacter(fData, val)
implicit none
    type(FileData), intent(inout) :: fData
    character, intent(out) :: val
    
    if (fData%ierr /= 0) then
        call createErrorMsg(fData, 1, ERROR_PARSING_CHARACTER)
        return
    end if
    read(fData%word, *, iostat=fData%ierr) val
    if (fData%ierr /= 0) call createErrorMsg(fData, 1, ERROR_PARSING_CHARACTER)
end subroutine parseCharacter

end module MultiscaleFileParser
