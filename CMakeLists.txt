##############################################################################################
# CMAKE Configurations for UppASD
# (C) Anders Bergman 2019
#
# Acknowledgements: 
# https://github.com/SethMMorton/cmake_fortran_template  (for starting template)
# https://github.com/openmeeg/findmkl_cmake/blob/master/cmake/FindMKL.cmake (for FindMKL)
# https://github.com/libigl/eigen/blob/master/cmake/FindFFTW.cmake (for FindFFTW)
# 
##############################################################################################

##############################################################################################
# Define the project and the depencies that it has
##############################################################################################
CMAKE_MINIMUM_REQUIRED(VERSION 2.8.5)
PROJECT(UppASD Fortran)
# Set the project version
SET(VERSION 1.2.3)
##############################################################################################

include(CMakeDetermineFortranCompiler)
message(STATUS "CMAKE_Fortran_COMPILER_ID ${CMAKE_Fortran_COMPILER_ID}.")
include(FindPython)
find_package(Python)
message(STATUS "Python_STDLIB ${Python_STDLIB}.")
message(STATUS "Python_SITELIB ${Python_SITELIB}")
message(STATUS "Python_LIBRARIES ${Python_LIBRARIES}.")
message(STATUS "Python_LIBRARY_DIRS ${Python_LIBRARY_DIRS}.")

##############################################################################################
# There is an error in CMAKE with this flag for pgf90.  Unset it
##############################################################################################
GET_FILENAME_COMPONENT(FCNAME ${CMAKE_Fortran_COMPILER} NAME)
set(CMAKE_Fortran_COMPILE_OPTIONS_PIE "")
IF(FCNAME STREQUAL "pgf90")
   UNSET(CMAKE_SHARED_LIBRARY_LINK_Fortran_FLAGS)
ENDIF(FCNAME STREQUAL "pgf90")
##############################################################################################


##############################################################################################
# Add our local modlues to the module path
##############################################################################################
SET(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake/")
##############################################################################################


##############################################################################################
# Uncomment if it is required that Fortran 90 is supported
##############################################################################################
IF(NOT CMAKE_Fortran_COMPILER_SUPPORTS_F90)
   MESSAGE(FATAL_ERROR "Fortran compiler does not support F90")
ENDIF(NOT CMAKE_Fortran_COMPILER_SUPPORTS_F90)
##############################################################################################


##############################################################################################
# Macro used to add sources one-by-one
##############################################################################################
macro (add_sources)
   file (RELATIVE_PATH _relPath "${PROJECT_SOURCE_DIR}/source" "${CMAKE_CURRENT_SOURCE_DIR}")
   foreach (_src ${ARGN})
      if (_relPath)
         list (APPEND SRCS "${_relPath}/${_src}")
      else()
         list (APPEND SRCS "${_src}")
      endif()
   endforeach()
   if (_relPath)
      # propagate SRCS to parent directory
      set (SRCS ${SRCS} PARENT_SCOPE)
   endif()
   #message(STATUS "${SRCS}")
endmacro()
#############################################################################################


##############################################################################################
# Set relevant user-controlled options for the compilation
##############################################################################################
include(CMakeDependentOption)
OPTION(USE_OPENMP "Use OpenMP for parallelization" ON)
OPTION(USE_CUDA "Compile with CUDA support" OFF)
OPTION(USE_VSL "Use Intel VSL (needs MKL)" OFF)
OPTION(USE_OVF "Use OVF format" OFF)
OPTION(BUILD_MANUAL "Build manual"  OFF)
OPTION(BUILD_DOXYGEN "Build doxygen"  OFF)
CMAKE_DEPENDENT_OPTION(USE_MKL "Use Intel MKL for BLAS/LAPACK" ON "NOT USE_VSL" OFF)
CMAKE_DEPENDENT_OPTION(USE_MKL_FFT "Use Intel MKL for FFT" OFF "NOT USE_VSL" ON)
CMAKE_DEPENDENT_OPTION(USE_FFTW "Use FFTW for FFT" OFF "NOT USE_VSL" ON)
CMAKE_DEPENDENT_OPTION(BUILD_CUDA_DOCS "Build LaTeX notes of the CUDA implementation" OFF "BUILD_MANUAL" OFF) 
CMAKE_DEPENDENT_OPTION(BUILD_GNEB_DOCS "Build LaTeX notes of the GNEB implementation" OFF "BUILD_MANUAL" OFF) 
CMAKE_DEPENDENT_OPTION(USE_FFTW "Use FFTW for FFT" OFF "NOT USE_VSL" ON)
##############################################################################################

##############################################################################################
# Build LaTeX manual (and extra notes if requested)
##############################################################################################
if(BUILD_MANUAL)
   include(${CMAKE_MODULE_PATH}/UseLATEX.cmake)
   add_subdirectory(docs/Manual)
   if(BUILD_GNEB_DOCS)
      add_subdirectory(docs/GNEB)
   endif()
   if(BUILD_CUDA_DOCS)
      add_subdirectory(docs/CUDA)
   endif()
endif()
##############################################################################################

##############################################################################################
# Build Doxygen documentation
##############################################################################################
if(BUILD_DOXYGEN)
   # check if Doxygen is installed
   find_package(Doxygen)
   if (DOXYGEN_FOUND)
      ## set input and output files
      set(DOXYGEN_IN ${CMAKE_CURRENT_SOURCE_DIR}/docs/Doxyfile)
      #set(DOXYGEN_OUT ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile)

      ## request to configure the file
      #configure_file(${DOXYGEN_IN} ${DOXYGEN_OUT} @ONLY)
      message(STATUS "Doxygen build started")

      # note the option ALL which allows to build the docs together with the application
    add_custom_target( doc_doxygen ALL
       COMMAND ${DOXYGEN_EXECUTABLE} ${DOXYGEN_IN}
       WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/docs
         COMMENT "Generating API documentation with Doxygen"
         VERBATIM )
   else (DOXYGEN_FOUND)
      message("Doxygen need to be installed to generate the doxygen documentation")
   endif (DOXYGEN_FOUND)
endif()
##############################################################################################



##############################################################################################
# If MKL/VSL is requested. First checks here to set the linear algebra packages
##############################################################################################
if(USE_VSL OR USE_MKL OR USE_MKL_FFT)
#if(USE_VSL OR USE_MKL)
   set(BLA_DEFINITIONS)

   #if ("$ENV{MKLROOT}" STREQUAL "")
   #   message(FATAL_ERROR "MKLROOT is not set. Please source the Intel MKL mklvars.sh file.")
   #endif()

   # user defined options for MKL
   #option(MKL_USE_parallel "Use MKL parallel" True)
   #option(MKL_USE_sdl "Single Dynamic Library or static/dynamic" False)
   #set(MKL_USE_interface "lp64" CACHE STRING "for Intel(R)64 compatible arch: ilp64/lp64 or for ia32 arch: cdecl/stdcall")
   set(MKL_USE_parallel  True)
   set(MKL_USE_sdl False)
   set(MKL_USE_interface "lp64")

   if (BLA_VENDOR MATCHES "_seq")
      set(MKL_USE_parallel OFF)
   else()
      set(MKL_USE_parallel ON)
   endif()

   if(USE_VSL)
      find_package(MKL REQUIRED)
   else()
      find_package(MKL)
   endif(USE_VSL)

   if (MKL_FOUND)
      set(BLA_INCLUDE_DIR ${MKL_INCLUDE_DIR})
      set(BLAS_LIBRARIES ${MKL_LIBRARIES})
      set(LAPACK_LIBRARIES ${MKL_LIBRARIES})
      set(BLA_DEFINITIONS USE_MKL)
   else()
      include(FindBLAS)
      include(FindLAPACK)
      find_package(BLAS)
      find_package(LAPACK REQUIRED)
   endif()

   #message(STATUS "BLAS Libraries: ${BLAS_LIBRARIES}")
   #include(${CMAKE_MODULE_PATH}/FindMKL)
   #find_package(MKL)
   #message(STATUS " Mkl: ${MKL_FOUND} ${LAPACK_FOUND} ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES}")
   message(STATUS "mkl-location: ${MKL_LIBRARY_DIR}")
   link_directories(${MKL_LIBRARY_DIR})
else()
   include(FindBLAS)
   include(FindLAPACK)
   find_package(BLAS)
   find_package(LAPACK REQUIRED)
   ###message(STATUS "MAth: ${BLAS_FOUND} ${LAPACK_FOUND} ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES}")
endif()
message(STATUS "Math libs: 
       BLAS:${BLAS_FOUND} ${BLAS_LIBRARIES} 
       LAPACK:${LAPACK_FOUND} ${LAPACK_LIBRARIES} 
       MKL:${MKL_FOUND} ${MKL_LIBRARIES}")
##############################################################################################


##############################################################################################
# Determine and set compiler-dependent flags
##############################################################################################
INCLUDE(${CMAKE_MODULE_PATH}/SetFortranFlags.cmake) 
##############################################################################################

##############################################################################################
# assume built-in pthreads on MacOS (Is this needed?)
##############################################################################################
IF(APPLE)
   set(CMAKE_THREAD_LIBS_INIT "-lpthread")
   set(CMAKE_HAVE_THREADS_LIBRARY 1)
   set(CMAKE_USE_WIN32_THREADS_INIT 0)
   set(CMAKE_USE_PTHREADS_INIT 1)
   set(THREADS_PREFER_PTHREAD_FLAG ON)
ENDIF()
##############################################################################################


##############################################################################################
# Extract git-stamp
##############################################################################################
INCLUDE(${CMAKE_MODULE_PATH}/git-watcher.cmake)
##############################################################################################

##############################################################################################
# Define the actual files and folders that make up the build
##############################################################################################

##############################################################################################
# Defines the directories for source, modules and object filed
SET(SRC ${CMAKE_SOURCE_DIR}/source)
#SET(MOD ${CMAKE_SOURCE_DIR}/source/mod)
SET(MOD ${CMAKE_CACHEFILE_DIR}/source/mod)
#SET(BIN ${CMAKE_SOURCE_DIR}/bin)
# Have the .mod files placed in the lib folder
SET(CMAKE_Fortran_MODULE_DIRECTORY ${MOD})
# Have the .o files placed in the bin folder
#set(CMAKE_BUILD_FILES_DIRECTORY ${CMAKE_CACHEFILE_DIR}/bin)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/bin)
##############################################################################################


##############################################################################################
# Defines that the sources and following CMakeLists.txt are indeed in the source dir
##############################################################################################
#add_subdirectory(${SRC} ${BIN})
add_subdirectory(${SRC})
##############################################################################################


##############################################################################################
# Add a distclean target to the Makefile
##############################################################################################
ADD_CUSTOM_TARGET(distclean 
   COMMAND ${CMAKE_COMMAND} -P ${CMAKE_SOURCE_DIR}/distclean.cmake
   )
##############################################################################################


##############################################################################################
# Print ALL variables (only for testing)
##############################################################################################
#get_cmake_property(_variableNames VARIABLES)
#list (SORT _variableNames)
#foreach (_variableName ${_variableNames})
#    message(STATUS "${_variableName}=${${_variableName}}")
#endforeach()
##############################################################################################

