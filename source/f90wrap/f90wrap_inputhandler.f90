! Module inputhandler defined in file Input/inputhandler.f90

subroutine f90wrap_read_parameters(ifile)
    use inputhandler, only: read_parameters
    implicit none
    
    integer, intent(in) :: ifile
    call read_parameters(ifile=ifile)
end subroutine f90wrap_read_parameters

subroutine f90wrap_change_constants
    use inputhandler, only: change_constants
    implicit none
    
    call change_constants()
end subroutine f90wrap_change_constants

subroutine f90wrap_inputhandler__get__sane_input(f90wrap_sane_input)
    use inputhandler, only: inputhandler_sane_input => sane_input
    implicit none
    logical, intent(out) :: f90wrap_sane_input
    
    f90wrap_sane_input = inputhandler_sane_input
end subroutine f90wrap_inputhandler__get__sane_input

subroutine f90wrap_inputhandler__set__sane_input(f90wrap_sane_input)
    use inputhandler, only: inputhandler_sane_input => sane_input
    implicit none
    logical, intent(in) :: f90wrap_sane_input
    
    inputhandler_sane_input = f90wrap_sane_input
end subroutine f90wrap_inputhandler__set__sane_input

! End of module inputhandler defined in file Input/inputhandler.f90

