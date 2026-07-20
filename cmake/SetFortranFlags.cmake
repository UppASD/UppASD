######################################################
# Determine and set the Fortran compiler flags we want 
######################################################

message(STATUS "Setting compiler flags.")
####################################################################
# Make sure that the default build type is RELEASE if not specified.
####################################################################
INCLUDE(${CMAKE_MODULE_PATH}/SetCompileFlag.cmake)

# Accept the standard CMake configuration names, case-insensitively, and
# normalise to their canonical spelling.  CMake looks up
# CMAKE_<LANG>_FLAGS_<CONFIG> with <CONFIG> upper-cased, so the historical
# all-caps spellings (RELEASE, DEBUG) keep selecting exactly the same flags --
# existing build trees, CMakeLists.txt.local and CI are unaffected.
#
# Coverage is this project's own configuration (-O2 plus gcov instrumentation).
# It was called TESTING; the old name is still accepted.
set(UPPASD_BUILD_TYPES Debug Release RelWithDebInfo MinSizeRel Coverage)

STRING(TOUPPER "${CMAKE_BUILD_TYPE}" BT)

IF(NOT BT)
    SET(CMAKE_BUILD_TYPE Release CACHE STRING
      "Choose the type of build: ${UPPASD_BUILD_TYPES}." FORCE)
    MESSAGE(STATUS "CMAKE_BUILD_TYPE not given, defaulting to Release")
ELSE()
    SET(_bt_canonical "")
    FOREACH(_bt IN LISTS UPPASD_BUILD_TYPES)
        STRING(TOUPPER "${_bt}" _bt_upper)
        IF(BT STREQUAL _bt_upper)
            SET(_bt_canonical "${_bt}")
        ENDIF()
    ENDFOREACH()
    # Deprecated spelling of Coverage.
    IF(BT STREQUAL "TESTING")
        SET(_bt_canonical "Coverage")
        MESSAGE(DEPRECATION "CMAKE_BUILD_TYPE=TESTING is deprecated; use Coverage")
    ENDIF()
    IF(NOT _bt_canonical)
        MESSAGE(FATAL_ERROR
          "CMAKE_BUILD_TYPE '${CMAKE_BUILD_TYPE}' not valid, choices are: ${UPPASD_BUILD_TYPES}")
    ENDIF()
    SET(CMAKE_BUILD_TYPE "${_bt_canonical}" CACHE STRING
      "Choose the type of build: ${UPPASD_BUILD_TYPES}." FORCE)
ENDIF()

# Offer the list as a drop-down in cmake-gui / ccmake.
SET_PROPERTY(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS ${UPPASD_BUILD_TYPES})
MESSAGE(STATUS "Build type: ${CMAKE_BUILD_TYPE}")

# Re-derive the upper-case form; the canonical spelling may differ from what
# the user passed, and the flag blocks below key off the upper-case name.
STRING(TOUPPER "${CMAKE_BUILD_TYPE}" BT)

#########################################################
# If the compiler flags have already been set, return now
#########################################################
#
IF(CMAKE_Fortran_FLAGS_RELEASE AND CMAKE_Fortran_FLAGS_COVERAGE AND CMAKE_Fortran_FLAGS_DEBUG
   AND CMAKE_Fortran_FLAGS_RELWITHDEBINFO AND CMAKE_Fortran_FLAGS_MINSIZEREL)
   #unset(CMAKE_Fortran_FLAGS CACHE)
    RETURN ()
ENDIF()

########################################################################
# Determine the appropriate flags for this compiler for each build type.
# For each option type, a list of possible flags is given that work
# for various compilers.  The first flag that works is chosen.
# If none of the flags work, nothing is added (unless the REQUIRED 
# flag is given in the call).  This way unknown compiles are supported.
#######################################################################

#####################
### GENERAL FLAGS ###
#####################


## Don't add underscores in symbols for C-compatability
#SET_COMPILE_FLAG(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS}"
#                 Fortran "-fno-underscoring")

# No limits on line-lengths with GNU
SET_COMPILE_FLAG(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS}"
                 Fortran "-ffree-line-length-0")
SET_COMPILE_FLAG(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS}"
                 Fortran "-std=legacy")

if(USE_VSL)
   SET_COMPILE_FLAG(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS}"
                 Fortran "-fno-range-check")
endif()

# Ensure that preprocessor flags are invoked
SET_COMPILE_FLAG(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS}"
                 Fortran "-cpp"
                         "-fpp")

# There is some bug where -march=native doesn't work on Mac
IF(APPLE)
    SET(GNUNATIVE "-mtune=native")
ELSE()
    SET(GNUNATIVE "-march=native")
ENDIF()
# Optimize for the host's architecture
SET_COMPILE_FLAG(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS}"
                 Fortran "-xHost"        # Intel
                         ${GNUNATIVE}    # GNU
                         "-ta=host"      # Portland Group
                         "/QxHost"       # Intel Windows
                )

###################
### DEBUG FLAGS ###
###################

# NOTE: debugging symbols (-g or /debug:full) are already on by default

# Disable optimizations
SET_COMPILE_FLAG(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG}"
                           Fortran REQUIRED "-O0" # All compilers not on Windows
                                            "/Od" # Intel Windows
                )

# Turn on all warnings 
SET_COMPILE_FLAG(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG}"
                 Fortran "-warn all" # Intel
                         "-Minform=warn" #Portland
                         "-Wall"     # GNU
                          "/warn:all" # Intel Windows
                )

# Traceback
SET_COMPILE_FLAG(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG}"
                 Fortran "-traceback"   # Intel/Portland Group
                         "-fbacktrace"  # GNU (gfortran)
                         "-ftrace=full" # GNU (g95)
                          "/traceback"   # Intel Windows
                )

# Check array bounds
SET_COMPILE_FLAG(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG}"
                 Fortran "-check bounds"  # Intel
                         "-fcheck=bounds" # GNU (New style)
                         "-fbounds-check" # GNU (Old style)
                         "-Mbounds"       # Portland Group
                          "/check:bounds"  # Intel Windows
                )

######################
### COVERAGE FLAGS ###
######################

# Optimizations
SET_COMPILE_FLAG(CMAKE_Fortran_FLAGS_COVERAGE "${CMAKE_Fortran_FLAGS_COVERAGE}"
                 Fortran REQUIRED "-O2" # All compilers not on Windows
                                  "/O2" # Intel Windows
                )
# gcov instrumentation.  This deliberately does not go through
# SET_COMPILE_FLAG: that macro probes a flag by compiling *and linking* a test
# program, and -fprofile-arcs fails to link without a matching --coverage on the
# link line.  The old TESTING configuration used the probe and so silently ended
# up with -ftest-coverage but no -fprofile-arcs, emitting .gcno notes and never
# producing the .gcda data that gcov actually reads.
IF(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
    SET(CMAKE_Fortran_FLAGS_COVERAGE "${CMAKE_Fortran_FLAGS_COVERAGE} --coverage"
        CACHE STRING "Fortran flags used for Coverage builds" FORCE)
ELSE()
    MESSAGE(WARNING
      "Coverage instrumentation is only wired up for gfortran; "
      "${CMAKE_Fortran_COMPILER_ID} will build optimised but uninstrumented.")
ENDIF()

#####################
### RELEASE FLAGS ###
#####################

# NOTE: agressive optimizations (-O3) are already turned on by default
# Optimizations
SET_COMPILE_FLAG(CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS_RELEASE}"
                 Fortran REQUIRED "-fastsse" # All compilers not on Windows
                                  "-Ofast" # All compilers not on Windows
                                  "-O3" # All compilers not on Windows
                                  "/O3" # Intel Windows
                )


# Unroll loops
SET_COMPILE_FLAG(CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS_RELEASE}"
                 Fortran "-funroll-loops" # GNU
                         "-unroll"        # Intel
                         "-Munroll"       # Portland Group
                         "/unroll"        # Intel Windows
                )

# Inline functions
SET_COMPILE_FLAG(CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS_RELEASE}"
                 Fortran "-inline"            # Intel
                         "-finline-functions" # GNU
                         "-Minline"           # Portland Group
                         "/Qinline"           # Intel Windows
                )

             ## Interprocedural (link-time) optimizations
             #SET_COMPILE_FLAG(CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS_RELEASE}"
             #                 Fortran "-ipo"     # Intel
             #                         "/Qipo"    # Intel Windows
             #                         "-Mipa=fast"    # Portland Group
             #)

# Single-file optimizations.  ifx accepts -ip only as a warning and ignores it;
# retain it for classic ifort but do not pass it to the LLVM-based compiler.
if(NOT CMAKE_Fortran_COMPILER_ID STREQUAL "IntelLLVM")
   SET_COMPILE_FLAG(CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS_RELEASE}"
                    Fortran "-ip"  # Intel classic
                            "-Mnoipa"    # Portland
                            "/Qip" # Intel Windows
                   )
endif()

# Vectorize code
SET_COMPILE_FLAG(CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS_RELEASE}"
                 Fortran "-qopt-report0"  # Intel
                         "-Mvect"        # Portland Group
                         "/Qvec-report0" # Intel Windows
                )

             
############################
### RELWITHDEBINFO FLAGS ###
############################

# Optimised, but with symbols and without -Ofast's fast-math relaxations, so a
# profile or a core dump maps back to source.  This is the configuration to use
# for Nsight/rocprof work.
SET_COMPILE_FLAG(CMAKE_Fortran_FLAGS_RELWITHDEBINFO "${CMAKE_Fortran_FLAGS_RELWITHDEBINFO}"
                 Fortran REQUIRED "-O2" # All compilers not on Windows
                                  "/O2" # Intel Windows
                )

SET_COMPILE_FLAG(CMAKE_Fortran_FLAGS_RELWITHDEBINFO "${CMAKE_Fortran_FLAGS_RELWITHDEBINFO}"
                 Fortran "-g"          # GNU/Intel/Portland
                         "/debug:full" # Intel Windows
                )

SET_COMPILE_FLAG(CMAKE_Fortran_FLAGS_RELWITHDEBINFO "${CMAKE_Fortran_FLAGS_RELWITHDEBINFO}"
                 Fortran "-funroll-loops" # GNU
                         "-unroll"        # Intel
                         "-Munroll"       # Portland Group
                         "/unroll"        # Intel Windows
                )

#########################
### MINSIZEREL FLAGS ###
#########################

SET_COMPILE_FLAG(CMAKE_Fortran_FLAGS_MINSIZEREL "${CMAKE_Fortran_FLAGS_MINSIZEREL}"
                 Fortran REQUIRED "-Os" # All compilers not on Windows
                                  "/O1" # Intel Windows
                )

##############################
### COVERAGE LINKER FLAGS ###
##############################

# gcov instrumentation has to reach the link line too, or the build fails on
# undefined __gcov_* symbols.  The old TESTING configuration set compile flags
# only, which is one reason it was not usable.
SET(CMAKE_EXE_LINKER_FLAGS_COVERAGE "--coverage" CACHE STRING
    "Linker flags used for Coverage builds")
SET(CMAKE_SHARED_LINKER_FLAGS_COVERAGE "--coverage" CACHE STRING
    "Shared-linker flags used for Coverage builds")
mark_as_advanced(CMAKE_EXE_LINKER_FLAGS_COVERAGE CMAKE_SHARED_LINKER_FLAGS_COVERAGE)

# Coverage is not one of CMake's built-in configurations, so the C/C++/GPU flag
# variables for it start out empty and those sources would compile unoptimised
# and uninstrumented.  Fortran is what the coverage numbers are about, so the
# GPU languages just get an optimisation level.
SET(CMAKE_C_FLAGS_COVERAGE "-O2 --coverage" CACHE STRING "C flags for Coverage builds")
SET(CMAKE_CXX_FLAGS_COVERAGE "-O2 --coverage" CACHE STRING "C++ flags for Coverage builds")
SET(CMAKE_CUDA_FLAGS_COVERAGE "-O2" CACHE STRING "CUDA flags for Coverage builds")
SET(CMAKE_HIP_FLAGS_COVERAGE "-O2" CACHE STRING "HIP flags for Coverage builds")
mark_as_advanced(CMAKE_C_FLAGS_COVERAGE CMAKE_CXX_FLAGS_COVERAGE
                 CMAKE_CUDA_FLAGS_COVERAGE CMAKE_HIP_FLAGS_COVERAGE)

mark_as_advanced(CMAKE_Fortran_FLAGS_COVERAGE)
list(REMOVE_DUPLICATES CMAKE_Fortran_FLAGS_RELEASE)
list(REMOVE_DUPLICATES CMAKE_Fortran_FLAGS_COVERAGE)
list(REMOVE_DUPLICATES CMAKE_Fortran_FLAGS_DEBUG)
list(REMOVE_DUPLICATES CMAKE_Fortran_FLAGS_RELWITHDEBINFO)
list(REMOVE_DUPLICATES CMAKE_Fortran_FLAGS_MINSIZEREL)
list(REMOVE_DUPLICATES CMAKE_Fortran_FLAGS)
