# Short hand for converting paths to absolute.
macro(PATH_TO_ABSOLUTE var_name)
    get_filename_component(${var_name} "${${var_name}}" ABSOLUTE)
endmacro()

# Check that a required variable is set.
macro(CHECK_REQUIRED_VARIABLE var_name)
    if(NOT DEFINED ${var_name})
        message(FATAL_ERROR "The \"${var_name}\" variable must be defined.")
    endif()
    PATH_TO_ABSOLUTE(${var_name})
endmacro()

# Check that an optional variable is set, or, set it to a default value.
macro(CHECK_OPTIONAL_VARIABLE var_name default_value)
    if(NOT DEFINED ${var_name})
        set(${var_name} ${default_value})
    endif()
    PATH_TO_ABSOLUTE(${var_name})
endmacro()

# Check the optional git variable.
# If it's not set, we'll try to find it using the CMake packaging system.
if(NOT DEFINED GIT_EXECUTABLE)
    find_package(Git QUIET REQUIRED)
endif()
CHECK_REQUIRED_VARIABLE(GIT_EXECUTABLE)


# Description: gets the current state of the git repo.
# Args:
#   _working_dir (in)  string; the directory from which git commands will be executed.
#   _state       (out) list; a collection of variables representing the state of the
#                            repository (e.g. commit SHA).
function(GetGitState _working_dir _gitlabel)
    # Get whether or not the working tree is dirty.
    execute_process(COMMAND
       "${GIT_EXECUTABLE}" describe --abbrev=4 --dirty --always --tags
        WORKING_DIRECTORY "${_working_dir}"
        OUTPUT_VARIABLE git_output
        ERROR_QUIET
        OUTPUT_STRIP_TRAILING_WHITESPACE)
    # Return a list of our variables to the parent scope.
    set(${_gitlabel} VERSION=\"${git_output}\" PARENT_SCOPE)
endfunction()

function(GetGitVersion _working_dir _gitversion)
    # Get whether or not the working tree is dirty.
    execute_process(COMMAND
       "${GIT_EXECUTABLE}" describe --abbrev=0 --tags
        WORKING_DIRECTORY "${_working_dir}"
        OUTPUT_VARIABLE GIT_VERSION
        ERROR_QUIET
        OUTPUT_STRIP_TRAILING_WHITESPACE)
    # If git describe fails, fallback to a default version
    if(NOT GIT_VERSION)
        set(GIT_VERSION "0.0.0")
    endif()
    
    # Remove the leading 'v' if present
    string(REGEX REPLACE "^v" "" GIT_VERSION ${GIT_VERSION})    
    set(${_gitversion} "${GIT_VERSION}" PARENT_SCOPE)
endfunction()

GetGitState("${_working_dir}" gitlabel)
GetGitVersion("${_working_dir}" gitversion)

# Ensure that preprocessor flags are invoked
add_compile_definitions(${gitlabel})

message(STATUS "Git tag found: ${gitlabel}.")
