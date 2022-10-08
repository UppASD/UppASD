!> authors
!> Edgar Mendez
!> Nikos  Ntallis
!> Manuel Pereiro

module FileInput
    use MultiscaleFileParser
    use Parameters
implicit none

    !> Procedure type for a parser function. 
    interface
        subroutine parserFunction(fData)
            use MultiscaleFileParser
        implicit none
            type(FileData), intent(inout) :: fData
        end subroutine parserFunction
    end interface

    !> Structure that may hold a pointer to any of the data types supported
    !!  by the fields. A 'one-fits-all' parser function can be used if none
    !!  of the other types fit the use.
    type GenericPtr
        integer, pointer :: intValue
        real(dblprec), pointer :: realValue
        logical, pointer :: logicalValue
        character, pointer :: charValue
        integer, dimension(:), pointer :: intArray
        real(dblprec), dimension(:), pointer :: realArray
        logical, dimension(:), pointer :: logicalArray
        procedure(parserFunction), pointer, nopass :: prsFunc
    end type

    !> Linked list of fields
    type InputFields
        type(InputField), pointer :: startField
        type(InputField), pointer :: endField
    end type

    !> Input field descriptor. Acts as a linked list node.
    type InputField
        character(len=64) :: fieldName
        logical :: isRequired
        logical :: isPresent
        logical :: canRepeat
        integer :: presentOnLine
        character(len=512) :: description
        type(GenericPtr) :: dataValue
        type(InputField), pointer :: nextField
    end type

    !> Add a field that will be parsed depending on the data type
    interface addField
        procedure addFieldInt
        procedure addFieldReal
        procedure addFieldLogical
        procedure addFieldCharacter
        procedure addFieldIntArray
        procedure addFieldRealArray
        procedure addFieldLogicalArray
        procedure addFieldParserFunc
    end interface


public loadInputFile, addField, InputFields, getDescription, allocateInputFields, &
        deallocateInputFields, getClosestMatchingField, parseExternalFile, parserFunction, &
        tryToParseFromExternalFile, toLowerCase

private 

contains

!< Initialises an InputFields structure
subroutine allocateInputFields(fields)
implicit none
    type(InputFields), intent(inout) :: fields
    
    nullify(fields%startField)
    nullify(fields%endField)
end subroutine allocateInputFields

!< Deallocates an InputFields structre
subroutine deallocateInputFields(fields)
implicit none
    type(InputFields), intent(inout) :: fields
    
    type(InputField), pointer :: currentField
    type(InputField), pointer :: nextField
    currentField => fields%startField
    do while (associated(currentField))
        nextField => currentField%nextField
        deallocate(currentField)
        currentField => nextField
    enddo
    nullify(fields%startField)
    nullify(fields%endField)
end subroutine deallocateInputFields

!< Appends an int field to the fields structure
subroutine addFieldInt(fields, fieldName, intValue, required, description)
implicit none
    type(InputFields), intent(inout) :: fields
    character(len=*), intent(in) :: fieldName
    integer, target, intent(in) :: intValue
    logical, intent(in) :: required
    character(len=*) :: description

    type(InputField), pointer :: field
    allocate(field)
    call initializeGenericPtr(field%dataValue)
    
    field%dataValue%intValue => intValue
    call genericAddField(fields, field, fieldName, required, .false., description)
end subroutine addFieldInt

!< Appends a real field to the fields structure
subroutine addFieldReal(fields, fieldName, realValue, required, description)
implicit none
    type(InputFields), intent(inout) :: fields
    character(len=*), intent(in) :: fieldName
    real(dblprec), target, intent(in) :: realValue
    logical, intent(in) :: required
    character(len=*) :: description

    type(InputField), pointer :: field
    allocate(field)
    call initializeGenericPtr(field%dataValue)
    
    field%dataValue%realValue => realValue
    call genericAddField(fields, field, fieldName, required, .false., description)
end subroutine addFieldReal

!< Appends a logical field to the fields structure
subroutine addFieldLogical(fields, fieldName, logicalValue, required, description)
implicit none
    type(InputFields), intent(inout) :: fields
    character(len=*), intent(in) :: fieldName
    logical, target, intent(in) :: logicalValue
    logical, intent(in) :: required
    character(len=*) :: description

    type(InputField), pointer :: field
    allocate(field)
    call initializeGenericPtr(field%dataValue)
    
    field%dataValue%logicalValue => logicalValue
    call genericAddField(fields, field, fieldName, required, .false., description)
end subroutine addFieldLogical

!< Appends a (single) character field to the fields structure
subroutine addFieldCharacter(fields, fieldName, charValue, required, description)
implicit none
    type(InputFields), intent(inout) :: fields
    character(len=*), intent(in) :: fieldName
    character, target, intent(in) :: charValue
    logical, intent(in) :: required
    character(len=*) :: description

    type(InputField), pointer :: field
    allocate(field)
    call initializeGenericPtr(field%dataValue)
    
    field%dataValue%charValue => charValue
    call genericAddField(fields, field, fieldName, required, .false., description)
end subroutine addFieldCharacter

!< Appends an int array field to the fields structure
subroutine addFieldIntArray(fields, fieldName, intArray, required, description)
implicit none
    type(InputFields), intent(inout) :: fields
    character(len=*), intent(in) :: fieldName
    integer, dimension(:), target, intent(in) :: intArray
    logical, intent(in) :: required
    character(len=*) :: description

    type(InputField), pointer :: field
    allocate(field)
    call initializeGenericPtr(field%dataValue)
    
    field%dataValue%intArray => intArray
    call genericAddField(fields, field, fieldName, required, .false., description)
end subroutine addFieldIntArray

!< Appends a real array field to the fields structure
subroutine addFieldRealArray(fields, fieldName, realArray, required, description)
implicit none
    type(InputFields), intent(inout) :: fields
    character(len=*), intent(in) :: fieldName
    real(dblprec), dimension(:), target, intent(in) :: realArray
    logical, intent(in) :: required
    character(len=*) :: description

    type(InputField), pointer :: field
    allocate(field)
    call initializeGenericPtr(field%dataValue)
    
    field%dataValue%realArray => realArray
    call genericAddField(fields, field, fieldName, required, .false., description)
end subroutine addFieldRealArray

!< Appends a logical array field to the fields structure
subroutine addFieldLogicalArray(fields, fieldName, logicalArray, required, description)
implicit none
    type(InputFields), intent(inout) :: fields
    character(len=*), intent(in) :: fieldName
    logical, dimension(:), target, intent(in) :: logicalArray
    logical, intent(in) :: required
    character(len=*) :: description

    type(InputField), pointer :: field
    allocate(field)
    call initializeGenericPtr(field%dataValue)
    
    field%dataValue%logicalArray => logicalArray
    call genericAddField(fields, field, fieldName, required, .false., description)
end subroutine addFieldLogicalArray

!< Appends a field with an specific parser function to the fields structur.e
subroutine addFieldParserFunc(fields, fieldName, prsFunc, required, repeatable, description)
implicit none
    type(InputFields), intent(inout) :: fields
    character(len=*), intent(in) :: fieldName
    procedure(parserFunction) :: prsFunc
    logical, intent(in) :: required
    logical, intent(in) :: repeatable
    character(len=*) :: description

    type(InputField), pointer :: field
    allocate(field)
    call initializeGenericPtr(field%dataValue)

    field%dataValue%prsFunc => prsFunc
    call genericAddField(fields, field, fieldName, required, repeatable, description)
end subroutine addFieldParserFunc

!< Appends a field to an InputFields structure.
subroutine genericAddField(fields, field, fieldName, required, repeatable, description)
    type(InputFields), intent(inout) :: fields
    type(InputField), target, intent(inout) :: field
    character(len=*), intent(in) :: fieldName
    logical, intent(in) :: required
    logical, intent(in) :: repeatable
    character(len=*) :: description

    field%fieldName = fieldName
    call toLowerCase(field%fieldName)
    field%isRequired = required
    field%isPresent = .false.
    field%canRepeat = repeatable
    field%description = description
    nullify(field%nextField)
    if (.not. associated(fields%endField)) then
        fields%endField => field
        fields%startField => field
    else
        fields%endField%nextField => field
        fields%endField => field
    end if
end subroutine genericAddField

!< Returns the field in fields with a given field name
!! @param fields The fields that is search among
!! @param fieldName The field name that is looked for
!! @param[out] field The field that is found, if it is not found this pointer is not associated
subroutine getField(fields, fieldName, field)
implicit none
    type(InputFields), intent(in) :: fields
    character(len=*), intent(in) :: fieldName
    type(InputField), pointer, intent(out) :: field
    
    type(InputField), pointer :: currentField
    currentField => fields%startField
    
    nullify(field)
    do while (associated(currentField))
        if (trim(currentField%fieldName) == trim(fieldName)) then
            field => currentField
            return
        endif
        currentField => currentField%nextField
    enddo
end subroutine getField

!< Gets the description of a field
!! @param fields The fields that is being search among
!! @param fieldName The name of the field 
!! @param[out] description The description of the field
!! @param[out] found If the field was found or not
subroutine getDescription(fields, fieldName, description, found)
implicit none
    type(InputFields), intent(in) :: fields
    character(len=*), intent(in) :: fieldName
    character(len=512), intent(out) :: description
    logical, intent(out) :: found

    type(InputField), pointer :: field
    
    call getField(fields, fieldName, field)
    if (associated(field)) then
        found = .true.
        description = field%description
    else
        found = .false.
        description = ''
    endif
end subroutine getDescription

!> Find closest matching field given a field name.
!! @param[in] fields Descriptors for known fields
!! @param[in] fieldName name of the field to be searched
!! @param[out] closestField name of the best-matching field
!! @param[out] distance distance value between the given and the found field names.
!! Distance is guaranteed to be 0 if and only if the field names are the same
subroutine getClosestMatchingField(fields, fieldName, closestField, distance)
implicit none
    type(InputFields), intent(in) :: fields
    character(len=*), intent(in) :: fieldName
    character(len=64), intent(out) :: closestField
    integer, intent(out) :: distance
    
    type(InputField), pointer :: currentField
    integer :: currentDistance
    distance = huge(distance)
    currentField => fields%startField
    do while (associated(currentField))
        currentDistance = damerauLevenshteinDistance(trim(currentField%fieldName), trim(fieldName))
        if (currentDistance < distance) then
            distance = currentDistance
            closestField = currentField%fieldName
        end if
        currentField => currentField%nextField
    end do
end subroutine

!> Initialize a generic pointer (do not use one before initialising!)
subroutine initializeGenericPtr(genericPointer)
implicit none
    type(GenericPtr), intent(inout) :: genericPointer

    nullify(genericPointer%intValue)
    nullify(genericPointer%realValue)
    nullify(genericPointer%logicalValue)
    nullify(genericPointer%charValue)
    nullify(genericPointer%realArray)
    nullify(genericPointer%logicalArray)
    nullify(genericPointer%intArray)
    nullify(genericPointer%realArray)
    nullify(genericPointer%logicalArray)
    nullify(genericPointer%prsFunc)
end subroutine initializeGenericPtr

!> Open a file and parse the fields
!! @param[in] filename The file path
!! @param[in] fields Descriptors of the fields and pointers to store parsed values
!! @param[out] ierr An error code if something went wrong, anthing other than 0 indicates an error
!! @param[out] errMsg A textual error message, in case ierr /= 0
subroutine loadInputFile(filename, fields, ierr, errMsg)
implicit none
    character(len=*), intent(in) :: filename
    type(InputFields), intent(in) :: fields
    integer, intent(out) :: ierr
    character(len=*), intent(out) :: errMsg
    
    type(FileData) :: fData
    call openFile(fData, filename)
    if (fData%ierr /= 0) then
        ierr = fData%ierr
        errMsg = 'The file ' // trim(fData%filename) // ' cant be opened'
        return
    endif
    
    call parseInputFile(fData, fields)
    close(fData%fileId)
    ierr = fData%ierr
    errMsg = fData%errMsg
    if (.not. is_iostat_end(ierr)) then !Indicates end of file
        if (ierr == 0) then
            errMsg = 'Something is at the end of the multiscale file'
            ierr = 1
        endif
        return
    else
        ierr = 0
    endif
end subroutine loadInputFile

!> Error message for fields that are missing in the fields structure.
!! @param[inout] fData FileData used to create the error message
!! @param[in] ierr error code (nonzero to cause the parsing to stop, specific codes are not being used now)
!! @param[in] fields Fields structure used to find the closest-named field
!! @param[in] unknownField name of the missing field
subroutine createUnknownFieldMsg(fData, ierr, fields, unknownField)
implicit none
    type(FileData), intent(inout) :: fData
    integer, intent(in) :: ierr
    type(InputFields), intent(in) :: fields
    character(len=*), intent(in) :: unknownField
    
    character(len=64) :: closestMatch
    integer :: matchDistance

    call getClosestMatchingField(fields, unknownField, closestMatch, matchDistance)
    call createErrorLineMsg(fData)
    write (fData%errMsg, '(a)') trim(fData%errMsg) // NEW_LINE('A') // 'Unkown field'
    if (len(trim(closestMatch)) >= 3*matchDistance) then
        write (fData%errMsg, '(a)') trim(fData%errMsg) // ', did you mean ' // trim(closestMatch) // '?'
    else
        write (fData%errMsg, '(a)') trim(fData%errMsg) // '.'
    end if
    fData%ierr = ierr
end subroutine

!< Attempts to parse all fields in an input file.
!! @param fData The file data object that is associated with the input file which will be parsed
!! @param setup Description of the fields that will be parsed (includes the destination for the parsed data)
subroutine parseInputFile(fData, fields)
implicit none
    type(FileData), intent(inout) :: fData
    type(InputFields), intent(in) :: fields
    
    type(InputField), pointer :: currentField

    do while (.not. is_iostat_end(fData%ierr))
        call toLowerCase(fData%word)
        call getField(fields, fData%word, currentField)
        if (.not. associated(currentField)) then
            call createUnknownFieldMsg(fData, 1, fields, fData%word)
            return
        else
            call parseField(fData, currentField)
            if (fData%ierr /= 0 .and. .not. is_iostat_end(fData%ierr)) return
        endif
    enddo
    if (fData%ierr == 0 .or. is_iostat_end(fData%ierr)) &
        call verifyInputFields(fields, fData)
end subroutine parseInputFile

!< Parses fields from an additional file
!! Used when a field in the configuration is specified to take its values from
!! another file. (e.g: lists that contain the 'file' keyword)
subroutine parseExternalFile(fData, filename, prsFunc)
implicit none
    type(FileData), intent(inout) :: fData
    character(len=*), intent(in) :: filename
    procedure(parserFunction) :: prsFunc

    type(FileData) :: parseFile
    
    call openFile(parseFile, filename)
    
    if (parseFile%ierr /= 0) then
        call copyErrorMsg(parseFile, fData)
        return
    end if
    call prsFunc(parseFile)
    
    close(parseFile%fileId)
    if (.not. is_iostat_end(parseFile%ierr)) then
        if (parseFile%ierr == 0) then
            call createErrorMsg(parseFile, 1, 'Some junk at the end of the file')
        end if
        call copyErrorMsg(parseFile, fData)
        return
    end if
end subroutine

!< Parses the 'file' keyword from a configuration file (fData), if the keyword
!! is present, parses the next word as a filename and applies the parser function
!! on its content.
subroutine tryToParseFromExternalFile(fData, parser)
implicit none
    type(FileData), intent(inout) :: fData
    procedure(parserFunction) :: parser

    if (trim(fData%word) == 'file') then
        call readNextWord(fData)
        if (isAtEndOfLine(fData)) then
            call createErrorMsg(fData, 1, 'A filename must be provided')
            return
        end if
        call parseExternalFile(fData, trim(fData%word), parser)
        if (fData%ierr /= 0) return
        call readNextLine(fData)
        if (fData%ierr /= 0) return
    else
        call createErrorMsg(fData, 1, 'Unknown argument ' // fData%word)
    end if
end subroutine


!> Parses a given field 
subroutine parseField(fData, field)
implicit none
    type(FileData), intent(inout) :: fData
    type(InputField), intent(inout) :: field
    character(len=256) :: errMsg
    character(len=10) :: lineNrString
    integer :: i
    
    if (field%isPresent .and. .not. field%canRepeat) then
        write (lineNrString, '(i10)') field%presentOnLine
        write (errMsg, '(4a)') 'The field ', trim(field%fieldName), ' is already defined on line ', trim(adjustl(lineNrString))
        call createErrorMsg(fData, 1, errMsg)
        return
    endif
    
    call readNextWord(fData)
    
    field%isPresent = .true.
    field%presentOnLine = fData%lineNr
    if (associated(field%dataValue%intValue)) then
       call parseInt(fData, field%dataValue%intValue)
       if (fData%ierr /= 0) return
       call readNextLine(fData)
       if (fData%ierr /= 0) return
    elseif (associated(field%dataValue%realValue)) then
       call parseReal(fData, field%dataValue%realValue)
       if (fData%ierr /= 0) return
       call readNextLine(fData)
       if (fData%ierr /= 0) return
    elseif (associated(field%dataValue%logicalValue)) then
       call parseLogical(fData, field%dataValue%logicalValue)
       if (fData%ierr /= 0) return
       call readNextLine(fData)
       if (fData%ierr /= 0) return
    elseif (associated(field%dataValue%charValue)) then
       call parseCharacter(fData, field%dataValue%charValue)
       if (fData%ierr /= 0) return
       call readNextLine(fData)
       if (fData%ierr /= 0) return
    elseif (associated(field%dataValue%intArray)) then
       do i = 1, ubound(field%dataValue%intArray, 1)
          call parseInt(fData, field%dataValue%intArray(i))
          if (fData%ierr /= 0) return
          if (i < ubound(field%dataValue%intArray, 1)) call readNextWord(fData)
       enddo
       call readNextLine(fData)
       if (fData%ierr /= 0) return
    elseif (associated(field%dataValue%realArray)) then
       do i = 1, ubound(field%dataValue%realArray, 1)
          call parseReal(fData, field%dataValue%realArray(i))
          if (fData%ierr /= 0) return
          if (i < ubound(field%dataValue%realArray, 1)) call readNextWord(fData)
       enddo
       call readNextLine(fData)
       if (fData%ierr /= 0) return
    elseif (associated(field%dataValue%logicalArray)) then
       do i = 1, ubound(field%dataValue%logicalArray, 1)
          call parseLogical(fData, field%dataValue%logicalArray(i))
          if (fData%ierr /= 0) return
          if (i < ubound(field%dataValue%logicalArray, 1)) call readNextWord(fData)
       enddo
       call readNextLine(fData)
       if (fData%ierr /= 0) return
    elseif (associated(field%dataValue%prsFunc)) then
       call field%dataValue%prsFunc(fData)
       if (fData%ierr /= 0) return
    else    
       call createErrorMsg(fData, 1,&
           'Something is wrong in the code, the field: ' &
           // trim(field%fieldName) // &
           ' is not associated with a type. Check parseField in module FileInput.')
       return
    endif
end subroutine parseField

!! Validates that the required parameters are present.
!! @param[in] fields fields structure to validate.
!! @param[inout] fData filedata used to generate an error if a field is missing.
subroutine verifyInputFields(fields, fData)
  implicit none
  type(InputFields), intent(in) :: fields
  type(FileData), intent(inout) :: fData
  
  type(InputField), pointer :: currentField
  
  currentField => fields%startField
  do while (associated(currentField))
     if (currentField%isRequired .and. .not. currentField%isPresent) then
        call createSimpleErrorMsg(fData, 1, &
             'The field ' // trim(currentField%fieldName) // &
             ' must be present in the input file')
        return
     endif
     currentField => currentField%nextField
  enddo
end subroutine verifyInputFields
  
!! Converts (inplace) a string to lower case
!! @param[inout] str 
subroutine toLowerCase(str)
implicit none
    character(len=*), intent(inout) :: str
    
    integer :: i
    do i = 1, len(str)
        if (str(i:i) >= 'A' .and. str(i:i) <= 'Z') then
            str(i:i) = achar(iachar(str(i:i)) + 32)
        endif
    enddo
end subroutine toLowerCase

!> string-to-string distance used for parameter suggestion when an
!! unknown parameter is found.
!! @param[in] string1 First string to compare 
!! @param[in] string2 Second string to compare
!! @return nonnegative integer, zero only if the strings are equal.
!! corresponds to the number of editions or permutations required to transform
!! one string into the other
integer function damerauLevenshteinDistance(string1, string2) result(d)
implicit none
    character(len=*), intent(in) :: string1
    character(len=*), intent(in) :: string2
    
    integer, dimension(-1:len(string1), -1:len(string2)) :: distance
    integer, dimension(256) :: da
    
    integer :: i, j, k, l, cost, maxDistance, db
    
    da = 0
    maxDistance = len(string1) + len(string2)
    distance(-1, -1) = maxDistance
    do i = 0, len(string1)
        distance(i, -1) = maxDistance
        distance(i, 0) = i
    end do
    do j = 0, len(string2)
        distance(-1, j) = maxDistance
        distance(0, j) = j
    end do
    
    do i = 1, len(string1)
        db = 0
        do j = 1, len(string2)
            k = da(iachar(string2(j:j)))
            l = db
            if (string1(i:i) == string2(j:j)) then
                cost = 0
                db = j
            else
                cost = 1
            end if
            
            distance(i, j) = minval((/distance(i - 1, j - 1) + cost, &
                                 distance(i, j - 1) + 1, &
                                 distance(i - 1, j) + 1, &
                                 distance(k - 1, l - 1) + (i - k - 1) + 1 + (j - l - 1) /))
                                     
        end do
        da(iachar(string1(i:i))) = i
    end do
    
    d = distance(len(string1), len(string2))
end function damerauLevenshteinDistance

end module FileInput
