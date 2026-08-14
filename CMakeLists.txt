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
#CMAKE_MINIMUM_REQUIRED(VERSION 2.8.12)
CMAKE_MINIMUM_REQUIRED(VERSION 3.13)
CMAKE_POLICY(VERSION 3.13)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CUDA_STANDARD 20)
set(CMAKE_HIP_STANDARD 20)
set(CMAKE_CUDA_STANDARD_REQUIRED TRUE)
set(CMAKE_HIP_STANDARD_REQUIRED TRUE)
set(CMAKE_CXX_STANDARD_REQUIRED TRUE)

PROJECT(UppASD Fortran)
# Set the project version
SET(VERSION 5.1.0)
SET(CMAKE_CXX_STANDARD 20)
##############################################################################################

##############################################################################################
# Check OS and compiler status
##############################################################################################
message(STATUS "Operating system= ${CMAKE_SYSTEM_NAME}")
include(CMakeDetermineFortranCompiler)
message(STATUS "CMAKE_Fortran_COMPILER_ID ${CMAKE_Fortran_COMPILER_ID} ${CMAKE_Fortran_COMPILER}.")
include(CMakeDetermineCCompiler)
message(STATUS "CMAKE_C_COMPILER_ID ${CMAKE_C_COMPILER_ID} ${CMAKE_C_COMPILER}.")
enable_language(C)
##############################################################################################


##############################################################################################
# Add our local modlues to the module path
##############################################################################################
SET(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake/")
##############################################################################################

##############################################################################################
# Create lib directory if missing
##############################################################################################
set(LIB ${CMAKE_SOURCE_DIR}/lib)
#message(STATUS "The lib is ${LIB}")
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${LIB})
message(STATUS "The lib is ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}")
file(MAKE_DIRECTORY ${CMAKE_LIBRARY_OUTPUT_DIRECTORY})
set(LEGACY_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/source)
message(STATUS "The legacy output directory is ${LEGACY_OUTPUT_DIRECTORY}")
##############################################################################################


##############################################################################################
# Defines the directories for source, modules and object filed
SET(SRC ${CMAKE_SOURCE_DIR}/source)
SET(MOD ${CMAKE_BINARY_DIR}/mod)
SET(BIN ${CMAKE_BINARY_DIR}/bin)
# Keep generated Fortran modules inside the active build tree.  Writing these
# files to source/mod makes concurrent builds and different compiler builds
# corrupt one another's module cache.
file(MAKE_DIRECTORY ${MOD})
SET(CMAKE_Fortran_MODULE_DIRECTORY ${MOD})
# Have the .o files placed in the bin folder
file(MAKE_DIRECTORY ${BIN})
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${BIN})
##############################################################################################


##############################################################################################
# Uncomment if it is required that Fortran 90 is supported
##############################################################################################
IF(NOT CMAKE_Fortran_COMPILER_SUPPORTS_F90)
   MESSAGE(FATAL_ERROR "Fortran compiler does not support F90")
ENDIF(NOT CMAKE_Fortran_COMPILER_SUPPORTS_F90)
##############################################################################################

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
# GPU backend and precision selection
#
# UPPASD_GPU_BACKEND and UPPASD_PRECISION are the canonical knobs.  USE_CUDA and
# USE_HIP are kept as deprecated aliases: existing build trees, CI and
# CMakeLists.txt.local still pass them.  Both directions are resolved here, so
# the rest of the file (and the gpu_files subdirectories) can keep testing
# USE_CUDA/USE_HIP until CM2-CM4 convert them.
##############################################################################################
set(UPPASD_GPU_BACKEND "OFF" CACHE STRING "GPU backend to build for: OFF, CUDA or HIP")
set_property(CACHE UPPASD_GPU_BACKEND PROPERTY STRINGS OFF CUDA HIP)
set(UPPASD_PRECISION "DOUBLE" CACHE STRING "Floating point precision: DOUBLE, SINGLE or MIXED")
set_property(CACHE UPPASD_PRECISION PROPERTY STRINGS DOUBLE SINGLE MIXED)
OPTION(USE_CUDA "Deprecated alias for -DUPPASD_GPU_BACKEND=CUDA" OFF)
OPTION(USE_HIP "Deprecated alias for -DUPPASD_GPU_BACKEND=HIP" OFF)
OPTION(USE_FAST_COPY "Use the pinned fast-copy path for Fortran<->GPU transfers" OFF)

if(USE_CUDA AND USE_HIP)
   message(FATAL_ERROR "Cannot use both CUDA and HIP simultaneously.")
endif()

string(TOUPPER "${UPPASD_GPU_BACKEND}" UPPASD_GPU_BACKEND)
if(UPPASD_GPU_BACKEND STREQUAL "NONE")
   set(UPPASD_GPU_BACKEND "OFF")
endif()

if(UPPASD_GPU_BACKEND STREQUAL "OFF")
   # Nothing requested directly; fall back to the deprecated booleans.
   if(USE_CUDA)
      message(DEPRECATION "USE_CUDA is deprecated; use -DUPPASD_GPU_BACKEND=CUDA")
      set(UPPASD_GPU_BACKEND "CUDA")
   elseif(USE_HIP)
      message(DEPRECATION "USE_HIP is deprecated; use -DUPPASD_GPU_BACKEND=HIP")
      set(UPPASD_GPU_BACKEND "HIP")
   endif()
elseif(NOT UPPASD_GPU_BACKEND STREQUAL "CUDA" AND NOT UPPASD_GPU_BACKEND STREQUAL "HIP")
   message(FATAL_ERROR
      "UPPASD_GPU_BACKEND must be OFF, CUDA or HIP (got '${UPPASD_GPU_BACKEND}')")
endif()

# Derive the legacy booleans from the resolved backend, and reject the case
# where the two spellings disagree rather than silently picking one.
if(UPPASD_GPU_BACKEND STREQUAL "CUDA")
   if(USE_HIP)
      message(FATAL_ERROR "UPPASD_GPU_BACKEND=CUDA conflicts with USE_HIP=ON")
   endif()
   set(USE_CUDA ON)
   set(USE_HIP OFF)
elseif(UPPASD_GPU_BACKEND STREQUAL "HIP")
   if(USE_CUDA)
      message(FATAL_ERROR "UPPASD_GPU_BACKEND=HIP conflicts with USE_CUDA=ON")
   endif()
   set(USE_CUDA OFF)
   set(USE_HIP ON)
else()
   set(USE_CUDA OFF)
   set(USE_HIP OFF)
endif()

string(TOUPPER "${UPPASD_PRECISION}" UPPASD_PRECISION)
if(UPPASD_PRECISION STREQUAL "MIXED")
   # MIXED is specified against storage_t/accum_t typedefs that do not exist
   # yet (Phase 4, PR1).  Refuse rather than quietly building DOUBLE.
   message(FATAL_ERROR
      "UPPASD_PRECISION=MIXED is not implemented yet; use DOUBLE or SINGLE")
elseif(NOT UPPASD_PRECISION STREQUAL "DOUBLE" AND NOT UPPASD_PRECISION STREQUAL "SINGLE")
   message(FATAL_ERROR
      "UPPASD_PRECISION must be DOUBLE, SINGLE or MIXED (got '${UPPASD_PRECISION}')")
endif()

message(STATUS "GPU backend: ${UPPASD_GPU_BACKEND}, precision: ${UPPASD_PRECISION}")

# The language kernel-bearing sources compile as, empty for a CPU-only build.
if(UPPASD_GPU_BACKEND STREQUAL "CUDA")
   set(UPPASD_GPU_LANGUAGE CUDA)
elseif(UPPASD_GPU_BACKEND STREQUAL "HIP")
   set(UPPASD_GPU_LANGUAGE HIP)
else()
   set(UPPASD_GPU_LANGUAGE "")
endif()

# Mark sources containing device code so they compile as CUDA/HIP rather than
# plain C++.  Source-file properties are directory-scoped and ${UppASD_LIB} is
# defined in this directory, so a bare set_source_files_properties() down in
# source/gpu_files/ would be silently ignored; TARGET_DIRECTORY (CMake >= 3.18)
# aims the property at the directory that owns the target.  This lets each
# gpu_files subdirectory own its own source list.
function(uppasd_mark_gpu_sources)
   if(NOT UPPASD_GPU_LANGUAGE)
      return()
   endif()
   set(_abs "")
   foreach(_src IN LISTS ARGN)
      if(IS_ABSOLUTE "${_src}")
         list(APPEND _abs "${_src}")
      else()
         list(APPEND _abs "${CMAKE_CURRENT_SOURCE_DIR}/${_src}")
      endif()
   endforeach()
   set_source_files_properties(${_abs}
      TARGET_DIRECTORY ${UppASD_LIB}
      PROPERTIES LANGUAGE ${UPPASD_GPU_LANGUAGE})
endfunction()
##############################################################################################


##############################################################################################
# Define the executable name
##############################################################################################
set(ASD_binary sd)
#if (${CMAKE_SYSTEM_NAME} MATCHES "Windows")
#   string(APPEND CMAKE_Fortran_COMPILER ".exe")
#   #message(STATUS "---> ${CMAKE_Fortran_COMPILER}")
#endif()
# Find the vernacular compiler name
string(REGEX MATCH gfortran|g95|ifort|pgf90|f95|f77|nv current_compiler  ${CMAKE_Fortran_COMPILER})
message(STATUS "Matched compiler: ${current_compiler}.  ${CMAKE_Fortran_COMPILER} ${FCNAME}")
if(current_compiler)
   string(APPEND ASD_binary ".${current_compiler}")
endif()
# add a `cuda` tag to the binary if used.
if(USE_CUDA)
   string(APPEND ASD_binary ".cuda")
elseif(USE_HIP)
   string(APPEND ASD_binary ".hip")
endif()
SET(UppASD_EXE ${ASD_binary})
message(STATUS "ASD_binary ${ASD_binary}.")
##############################################################################################


##############################################################################################
# Define the executable in terms of the source files
##############################################################################################
ADD_EXECUTABLE(${UppASD_EXE} ${SRC}/0sd.f90 )
##############################################################################################


##############################################################################################
# Define the library name
##############################################################################################
set(ASD_library asdlib)
SET(UppASD_LIB ${ASD_library})
##############################################################################################


##############################################################################################
# Define the executable in terms of the source files
##############################################################################################
set(lib_type STATIC)
if (BUILD_PYTHON)
   set(CMAKE_POSITION_INDEPENDENT_CODE ON)
endif(BUILD_PYTHON)
message(STATUS " LIB INFO: ${UppASD_LIB} ${lib_type} ${SRC}")
ADD_LIBRARY(${UppASD_LIB} ${lib_type})
set_target_properties(${UppASD_LIB} PROPERTIES Fortran_MODULE_DIRECTORY "${MOD}")
set_target_properties(${UppASD_EXE} PROPERTIES Fortran_MODULE_DIRECTORY "${MOD}")
target_include_directories(${UppASD_LIB} PUBLIC "${MOD}")
target_include_directories(${UppASD_EXE} PRIVATE "${MOD}")
##############################################################################################


##############################################################################################
include(${CMAKE_MODULE_PATH}/FindFFTW.cmake)
if (FFTW_FOUND)
   message(STATUS "FFTW found ${FFTW_FOUND} ${FFTW_INCLUDES} ${FFTW_INCLUDE_DIR}")
else (FFTW_FOUND)
   set(FFTW_LIB "")
   set(FFTW_LIBRARIES "")
endif (FFTW_FOUND)
##############################################################################################

##############################################################################################
# Set relevant user-controlled options for the compilation
##############################################################################################
include(CMakeDependentOption)
OPTION(USE_OPENMP "Use OpenMP for parallelization" ON)
# USE_CUDA/USE_HIP/USE_FAST_COPY and the UPPASD_GPU_BACKEND / UPPASD_PRECISION
# knobs they alias are declared and resolved above, before the binary is named.
# Keep option-dependent Fortran objects from being reused after a build tree is
# reconfigured.  In particular this invalidates objects compiled with the old
# __INTEL_MKL__ sparse path, and distinguishes CUDA/HIP/no-GPU variants.
target_compile_definitions(${UppASD_LIB} PRIVATE
   UPPASD_CMAKE_CONFIG_STAMP=1
   UPPASD_CMAKE_USE_CUDA=${USE_CUDA}
   UPPASD_CMAKE_USE_HIP=${USE_HIP}
   UPPASD_MKL_LEGACY_SPARSE=0)
if(UPPASD_PRECISION STREQUAL "SINGLE")
   target_compile_definitions(${UppASD_LIB} PRIVATE SINGLE_PREC)
endif()
if(USE_FAST_COPY)
   target_compile_definitions(${UppASD_LIB} PRIVATE USE_FAST_COPY)
endif()
CMAKE_DEPENDENT_OPTION(USE_CUCOMPLEX_CORRELATIONS "Use native cuDoubleComplex for GPU correlations (Thrust by default)" OFF "USE_CUDA" OFF)
# The cuDoubleComplex correlations path is unsupported: its only source file,
# correlations/gpuCorrelations_cucomplex.cu, is not attached to any target and is
# not stream-correct (it launches on the default stream). Turning this option ON
# used to silently build the Thrust path anyway; fail instead of misleading.
if(USE_CUCOMPLEX_CORRELATIONS)
  message(FATAL_ERROR
    "USE_CUCOMPLEX_CORRELATIONS is not supported: gpuCorrelations_cucomplex.cu is "
    "not built and is not stream-correct. Use the Thrust correlations path "
    "(USE_CUCOMPLEX_CORRELATIONS=OFF). See the header comment in that file to revive it.")
endif()
CMAKE_DEPENDENT_OPTION(ON_LUMI "Use Thrust complex for HIP correlations on LUMI (hipComplex by default)" OFF "USE_HIP" OFF)
# Linear algebra providers are selected explicitly; merely having oneAPI in the
# environment must not make a non-MKL build compile or link against MKL.
OPTION(USE_MKL "Use Intel MKL for BLAS/LAPACK" OFF)
CMAKE_DEPENDENT_OPTION(USE_VSL "Use Intel VSL (needs MKL)" OFF "USE_MKL" OFF)
CMAKE_DEPENDENT_OPTION(USE_MKL_FFT "Use Intel MKL for FFT" OFF "USE_MKL" OFF)
CMAKE_DEPENDENT_OPTION(USE_FFTW "Use FFTW for FFT" OFF "FFTW_FOUND" OFF)
CMAKE_DEPENDENT_OPTION(BUILD_CUDA_DOCS "Build LaTeX notes of the CUDA implementation" OFF "BUILD_MANUAL" OFF) 
####### for hip here???
CMAKE_DEPENDENT_OPTION(BUILD_GNEB_DOCS "Build LaTeX notes of the GNEB implementation" OFF "BUILD_MANUAL" OFF) 
OPTION(USE_OVF "Use OVF format" OFF)
OPTION(BUILD_MANUAL "Build manual"  OFF)
OPTION(BUILD_DOXYGEN "Build doxygen"  OFF)
OPTION(BUILD_PYTHON "Build Python interface"  OFF)

##############################################################################################
# Enable and define tests
##############################################################################################
include(CTest)
if(BUILD_TESTING)
   enable_language(CXX)
endif()
CMAKE_DEPENDENT_OPTION(RUN_ASD_TESTS "Run asd-tests"  ON "BUILD_TESTING" OFF)
CMAKE_DEPENDENT_OPTION(RUN_SLD_TESTS "Run sld-tests"  OFF "BUILD_TESTING" OFF)
CMAKE_DEPENDENT_OPTION(RUN_GNEB_TESTS "Run gneb-tests"  OFF "BUILD_TESTING" OFF)
CMAKE_DEPENDENT_OPTION(RUN_REG_TESTS "Run regression-test"  ON "BUILD_TESTING" OFF)
CMAKE_DEPENDENT_OPTION(RUN_CUDA_TESTS "Run cuda-test"  ON "BUILD_TESTING;USE_CUDA" OFF)
#######hip here????
##############################################################################################

##############################################################################################
# Determine and set compiler-dependent flags
##############################################################################################
INCLUDE(${CMAKE_MODULE_PATH}/SetFortranFlags.cmake) 
##############################################################################################

##############################################################################################
# Add OpenMP compiler/linker flags
##############################################################################################
IF(USE_OPENMP)
   FIND_PACKAGE(OpenMP)
   message(STATUS "Did we find OpenMP? ${OpenMP_Fortran_FOUND}")
   message(STATUS "Which are the OpenMP flags?  ${OpenMP_Fortran_FLAGS}")


   if(OpenMP_Fortran_FOUND)
        set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${OpenMP_Fortran_FLAGS}")
        target_link_libraries(${UppASD_LIB} PRIVATE OpenMP::OpenMP_Fortran)

   endif()

ENDIF(USE_OPENMP)
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
      #######hip here????
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
# Select exactly one linear algebra provider.  The LP64 interface matches
# UppASD's default Fortran INTEGER kind on every supported compiler.
##############################################################################################
set(UPPASD_LINALG_LIBRARIES)
if(USE_MKL)
   if(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
      set(_uppasd_mkl_interface gf_lp64)
   else()
      set(_uppasd_mkl_interface intel_lp64)
   endif()
   # Let UppASD own OpenMP parallelism.  This avoids loading a second
   # libiomp5 from an MKL installation that differs from the compiler's
   # oneAPI release.  Advanced users can select intel_thread explicitly.
   set(UPPASD_MKL_THREADING "sequential" CACHE STRING
       "MKL threading layer (sequential, intel_thread, or gnu_thread)")
   set_property(CACHE UPPASD_MKL_THREADING PROPERTY STRINGS sequential intel_thread gnu_thread)
   set(MKL_THREADING "${UPPASD_MKL_THREADING}" CACHE STRING "MKL threading layer" FORCE)
   set(MKL_INTERFACE_FULL "${_uppasd_mkl_interface}" CACHE STRING "MKL Fortran interface" FORCE)
   find_package(MKL CONFIG REQUIRED)
   if(NOT TARGET MKL::MKL)
      message(FATAL_ERROR "USE_MKL=ON requires the MKL::MKL imported target")
   endif()
   list(APPEND UPPASD_LINALG_LIBRARIES MKL::MKL)
   message(STATUS "Linear algebra provider: Intel MKL (${MKL_INTERFACE_FULL})")
else()
   if(NOT DEFINED BLA_VENDOR OR BLA_VENDOR STREQUAL "")
      set(BLA_VENDOR Generic CACHE STRING "BLAS/LAPACK implementation when USE_MKL=OFF")
   endif()
   find_package(BLAS REQUIRED)
   find_package(LAPACK REQUIRED)
   list(APPEND UPPASD_LINALG_LIBRARIES ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})
   message(STATUS "Linear algebra provider: ${BLA_VENDOR}")
endif()
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
# Add CUDA flags and preprocessor commands
##############################################################################################
if (USE_CUDA) 

   # CMAKE_CUDA_ARCHITECTURES=native, used as the default below, needs 3.24.
   if(CMAKE_VERSION VERSION_LESS "3.24")
      message(FATAL_ERROR
         "GPU builds require CMake >= 3.24 (found ${CMAKE_VERSION}); "
         "CPU-only builds still work with the project minimum.")
   endif()

   ENABLE_LANGUAGE(C CUDA CXX)
   find_package(CUDAToolkit REQUIRED)
   target_sources(${UppASD_LIB} PRIVATE ${SRC}/chelper.f90)
   # The GPU source lists, and which of them carry device code, live in
   # source/gpu_files/ and its subdirectories -- see uppasd_mark_gpu_sources().
   add_subdirectory(${SRC}/gpu_files)
   target_compile_definitions(${UppASD_LIB} PRIVATE CUDA_V)

   # Default to the architecture of the GPU in this machine.  Machines with no
   # visible GPU -- compile-only CI, cross-compilation -- must pass one
   # explicitly, e.g. -DCMAKE_CUDA_ARCHITECTURES=80 or "80;90".
   if(NOT CMAKE_CUDA_ARCHITECTURES)
      set(CMAKE_CUDA_ARCHITECTURES native)
   endif()
   SET_TARGET_PROPERTIES(${UppASD_LIB} PROPERTIES CUDA_ARCHITECTURES "${CMAKE_CUDA_ARCHITECTURES}")
   message(STATUS "CMAKE_CUDA_ARCHITECTURES: ${CMAKE_CUDA_ARCHITECTURES}")

   # Link the toolkit libraries that are actually used, rather than the blanket
   # CUDA::toolkit.  NVTX is a genuine dependency the blanket link was hiding:
   # gpuEventPool.cpp and measurement/fortranMeasurement.cpp include
   # nvToolsExtCuda.h.  Its imported target was renamed across toolkit versions.
   target_link_libraries(${UppASD_LIB} PRIVATE CUDA::cudart CUDA::curand CUDA::cufft)
   if(TARGET CUDA::nvtx3)
      target_link_libraries(${UppASD_LIB} PRIVATE CUDA::nvtx3)
   elseif(TARGET CUDA::nvToolsExt)
      target_link_libraries(${UppASD_LIB} PRIVATE CUDA::nvToolsExt)
   endif()

   # The imported targets carry the toolkit include directories, but the plain
   # CXX sources are compiled without them, so keep an explicit include path.
   target_include_directories(${UppASD_LIB} PRIVATE ${CUDAToolkit_INCLUDE_DIRS})
   if(DEFINED CUDAToolkit_MATH_INCLUDE_DIR)
      target_include_directories(${UppASD_LIB} PRIVATE ${CUDAToolkit_MATH_INCLUDE_DIR})
   endif()
   foreach(_cuda_inc IN LISTS CUDAToolkit_INCLUDE_DIRS)
       if(EXISTS "${_cuda_inc}/cccl/thrust/complex.h")
           target_include_directories(${UppASD_LIB} PRIVATE "${_cuda_inc}/cccl")
           message(STATUS "Added CUDA CCCL include dir: ${_cuda_inc}/cccl")
       endif()
   endforeach()

   if(EXISTS "${CUDAToolkit_TARGET_DIR}/include/cccl/thrust/complex.h")
       target_include_directories(${UppASD_LIB} PRIVATE "${CUDAToolkit_TARGET_DIR}/include/cccl")
       message(STATUS "Added CUDA CCCL include dir: ${CUDAToolkit_TARGET_DIR}/include/cccl")
   endif()
elseif(USE_HIP)
   find_package(HIP REQUIRED)
   target_compile_definitions(${UppASD_LIB} PRIVATE HIP_V)
   if(ON_LUMI)
     target_compile_definitions(${UppASD_LIB} PRIVATE ON_LUMI)
     message(STATUS "ON_LUMI: using Thrust complex for HIP correlations")
   endif()

   ENABLE_LANGUAGE(HIP C CXX)
    set(CMAKE_HIP_COMPILER hipcc)
    set(CMAKE_CXX_COMPILER hipcc)
   #set(CMAKE_CXX_COMPILER "/opt/rocm-6.3.3/bin/hipcc")
   #set(CMAKE_HIP_COMPILER "/opt/rocm-6.3.3/bin/hipcc")
   message(STATUS "HIP CXX compiler: ${CMAKE_CXX_COMPILER}")
   target_sources(${UppASD_LIB} PRIVATE ${SRC}/chelper.f90)
   add_subdirectory(${SRC}/gpu_files)
  
    # HIP include and compile flags
   set(PP_HIP_COMPILER_FLAG "${CMAKE_Fortran_COMPILER_ID}")
   set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -D__${PP_HIP_COMPILER_FLAG}__")
   set(CMAKE_HIP_FLAGS "${CMAKE_HIP_FLAGS} -D__${PP_HIP_COMPILER_FLAG}__")
   message(STATUS "HIP PPFLAG: __${PP_HIP_COMPILER_FLAG}__")
   
   
    set(ROCM_ROOT $ENV{ROCM_PATH})
    #link_directories("/opt/rocm/lib")
    
    target_include_directories(${UppASD_LIB} PRIVATE
        ${ROCM_ROOT}/include
    )
    

    # hiprand
    find_library(HIPRAND_LIBRARY
        NAMES hiprand
        PATHS /opt/rocm/lib /opt/rocm/lib64
        REQUIRED
    )

    target_link_libraries(${UppASD_LIB} PRIVATE ${HIPRAND_LIBRARY})

   if(NOT CMAKE_HIP_ARCHITECTURES)
      set(CMAKE_HIP_ARCHITECTURES "gfx90a")
   endif()
   SET_TARGET_PROPERTIES(${UppASD_LIB} PROPERTIES HIP_ARCHITECTURES "${CMAKE_HIP_ARCHITECTURES}")
   set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} --offload-arch=${CMAKE_HIP_ARCHITECTURES}")
   message(STATUS "HIP GPU architecture: ${CMAKE_HIP_ARCHITECTURES}")

else ()
   target_sources(${UppASD_LIB} PRIVATE ${SRC}/Tools/nocuda.f90)
endif ()
##############################################################################################


##############################################################################################
# Link with OVF library if available 
##############################################################################################
IF(USE_OVF)
   find_library(OVF_LIB ovf /lib/python3.7/site-packages/ovf/ /usr/lib64/)
   message(STATUS "Looking for OVF library, found in ${OVF_LIB}.")
   add_compile_definitions(USE_OVF)
   target_sources(${UppASD_LIB} PRIVATE ${SRC}/Third_party/ovf/fortran/ovf.f90)
ELSE()
   unset(OVF_LIB)
ENDIF(USE_OVF)
##############################################################################################

##############################################################################################
# MKL support for BLAS/LAPACK
##############################################################################################
if(USE_MKL)
   target_include_directories(${UppASD_LIB} PRIVATE ${MKL_INCLUDE})
   # The legacy mkl_dbsrmv/mkl_dcsrmm entry points used by sparse.f90 are
   # not exported by current oneMKL releases.  Do not enable that code path:
   # the portable OpenMP sparse implementation remains available while MKL
   # continues to provide the BLAS/LAPACK routines selected above.
endif()
##############################################################################################


##############################################################################################
# Include VSL preprocessor flag (recommended)
##############################################################################################
if(USE_VSL AND MKL_FOUND)
   add_compile_definitions(VSL)
   #add_compile_definitions(USE_VSL)
      include_directories($ENV{MKLROOT}/include)
endif()
##############################################################################################


##############################################################################################
# Check for FFT treatment (MKL default, FFTW available)
##############################################################################################
if(USE_MKL_FFT AND MKL_FOUND)
   target_sources(${UppASD_LIB} PRIVATE ${SRC}/Hamiltonian/fftdipole_mkl.f90)
   add_compile_definitions(USE_MKL_FFT)
   target_include_directories(${UppASD_LIB} PRIVATE ${MKL_INCLUDE})
   #message(STATUS "FFT-> ${USE_MKL_FFT} ")
elseif(USE_FFTW)
   target_sources(${UppASD_LIB} PRIVATE ${SRC}/Hamiltonian/fftdipole_fftw.f90)
   add_compile_definitions(USE_FFTW)
   message(STATUS "FFT-> ${USE_FFTW} ")
   #include(${CMAKE_MODULE_PATH}/FindFFTW.cmake)
   include_directories(${FFTW_INCLUDE_DIR})
   include_directories(${FFTW_INCLUDES})
   message(STATUS "FFTW: -I${FFTW_INCLUDES}  -L${FFTW_LIBRARIES} ${FFTW_LIB}")
endif()
##############################################################################################


##############################################################################################
# Add a distclean target to the Makefile
##############################################################################################
ADD_CUSTOM_TARGET(distclean 
   COMMAND ${CMAKE_COMMAND} -P ${CMAKE_SOURCE_DIR}/distclean.cmake
   )
##############################################################################################


##############################################################################################
# Build Python interface
###############################################################################################
if(BUILD_PYTHON)
   # check if Python is installed
   include(FindPython)
   find_package(Python REQUIRED COMPONENTS Interpreter Development.Module NumPy)
   #find_package(Python)
   #message(STATUS "Python_STDLIB ${Python_STDLIB}.")
   #message(STATUS "Python_SITELIB ${Python_SITELIB}")
   #message(STATUS "Python_LIBRARIES ${Python_LIBRARIES}.")
   #message(STATUS "Python_LIBRARY_DIRS ${Python_LIBRARY_DIRS}.")

   if (Python_FOUND)
      # Grab the variables from a local Python installation
      # F2PY headers
      execute_process(
         COMMAND "${Python_EXECUTABLE}"
         -c "import numpy.f2py; print(numpy.f2py.get_include())"
         OUTPUT_VARIABLE F2PY_INCLUDE_DIR
         OUTPUT_STRIP_TRAILING_WHITESPACE
         )

      # Common variables
      set(f2py_module_name "uppasd")
      SET(UppASD_F90w f90_wraps)
      set(F90wrap_kind
         ${CMAKE_SOURCE_DIR}/cmake/.f2py_f2cmapPRIVATE
         )

      # Uncomment if f90wrap needs to be invoked
      #set(F90wrap_sources 
      #   ${SRC}/pyasd.f90 
      #   ${SRC}/uppasd.f90 
      #   ${SRC}/Input/inputdatatype.f90
      #   ${SRC}/Input/inputdata.f90
      #   ${SRC}/Input/inputhandler.f90
      #   ${SRC}/Input/inputhandler_ext.f90
      #   ${SRC}/System/momentdata.f90 
      #   ${SRC}/System/simulationdata.f90
      #   )
      set(F2PY_sources 
         ${SRC}/f90wrap/f90wrap_pyasd.f90
         ${SRC}/f90wrap/f90wrap_uppasd.f90
         ${SRC}/f90wrap/f90wrap_inputdatatype.f90
         ${SRC}/f90wrap/f90wrap_inputdata.f90
         ${SRC}/f90wrap/f90wrap_inputhandler.f90
         ${SRC}/f90wrap/f90wrap_inputhandler_ext.f90
         ${SRC}/f90wrap/f90wrap_simulationdata.f90
         ${SRC}/f90wrap/f90wrap_momentdata.f90
         )

      # Bild library list for f2py
      ########hip here??
      foreach ( elem ${OpenMP_Fortran_LIBRARIES} ${OVF_LIB} ${FFTW_LIB} ${CUDA_curand_LIBRARY})
         message(STATUS "Element ${elem}")
         string(SUBSTRING ${elem} 0 2 head)
         if(${head} MATCHES "/")
            if(NOT ${elem} MATCHES "\\.dylib" AND NOT ${elem} MATCHES "\\.so" AND NOT ${elem} MATCHES "\\.a")
               #if(NOT ${elem} MATCHES ".dylib" )
               string(PREPEND elem "-L")
            endif()
            if(${elem} MATCHES "libfftw3.dylib")
               string(LENGTH ${elem} str_len)
               string(LENGTH "libfftw3.dylib" str_lib)
               math(EXPR str_sub "${str_len}-${str_lib}")
               string(SUBSTRING ${elem} 0  ${str_sub}  leader)

               list(APPEND F2PY_libs -L${leader})
               set(elem "-lfftw3")
            endif()
            if(${elem} MATCHES "libgomp.dylib")
               set(elem "-lgomp")
            endif()
         endif()
         list(APPEND F2PY_libs ${elem})
      endforeach()
      list(APPEND F2PY_libs "-lgfortran")
      list(REMOVE_DUPLICATES F2PY_libs)
      foreach (elem ${F2PY_libs})
         message(STATUS "------------------ ${elem}")
      endforeach()


      if(current_compiler MATCHES "gfortran" OR current_compiler MATCHES "g95" OR current_compiler MATCHES "g95")
         set(F2PY_compiler --fcompiler=gnu95)
         set(F2PY_f90flags "--f90flags='${fflags}'")
      elseif (current_compiler MATCHES "ifort")
         set(F2PY_compiler --fcompiler=intelem)
         string(STRIP ${CMAKE_Fortran_FLAGS} fflags)
         set(F2PY_f90flags "--f90flags='${fflags}'")
      elseif (current_compiler MATCHES "nv")
         set(F2PY_compiler --fcompiler=nv)
      endif()
      #message(STATUS "FORTRAN FLAGS= ${CMAKE_Fortran_FLAGS}")
      #message(STATUS "F2PY_compiler= ${F2PY_compiler}")
      #message(STATUS "F2PY_flags   = ${F2PY_f90flags}")
      #message(STATUS "F2PY_flags   = --f90flags='${F2PY_f90flags} -O2'")


      add_custom_target ( ${UppASD_F90w}  ALL
         DEPENDS ${UppASD_LIB} 
         # Uncomment if f90wrap needs to be invoked
         #DEPENDS ${UppASD_LIB} ${F2PY_sources}
         WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
         COMMAND f2py -c ${F2PY_compiler}  ${F2PY_sources} -m _${f2py_module_name} -I${CMAKE_Fortran_MODULE_DIRECTORY} -L${CMAKE_BINARY_DIR} -l${UppASD_LIB} ${F2PY_libs} --link-blas_opt --link-lapack_opt
         COMMENT "Generating Python modules"
         VERBATIM 
         )

      # Not used by default (invoke only if new header files need to be generated)
      add_custom_command( 
         OUTPUT ${F2PY_sources}
         WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
         COMMAND f90wrap  -P -k ${F90wrap_kind} -m ${f2py_module_name}  ${F90wrap_sources}
         COMMENT "Writing F90wrappers for F2PY"
         DEPENDS ${F90wrap_sources}
         )

      add_dependencies(${UppASD_F90w} ${UppASD_LIB})

      add_custom_command(
         TARGET ${UppASD_F90w} POST_BUILD
         COMMAND ${CMAKE_COMMAND} -E copy
         ${CMAKE_BINARY_DIR}/_*.so
         ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}
         )

      else (Python_FOUND)
      message("A proper Python environment needs to be installed to generate the python libraries")
      endif (Python_FOUND)
endif(BUILD_PYTHON)
##############################################################################################


##############################################################################################
mark_as_advanced( ENABLE_BLACS ENABLE_BLAS95 ENABLE_CDFT ENABLE_CPARDISO ENABLE_LAPACK95 ENABLE_OMP_OFFLOAD ENABLE_SCALAPACK MKL_DIR MKL_H MKL_LINK MKL_MPI mkl_gf_ilp64_file mkl_intel_thread_file mkl_core_file mkl_intel_ilp64_file mkl_gnu_thread_file mkl_intel_lp64_file mkl_gf_lp64_file OMP_LIBRARY)

##############################################################################################
# Add the needed libraries and special compiler flags
##############################################################################################
########hip here?????
TARGET_LINK_LIBRARIES(${UppASD_LIB} PRIVATE
   ${UPPASD_LINALG_LIBRARIES}
   ${CMAKE_THREAD_LIBS_INIT} 
   ${OVF_LIB} 
   ${FFTW_LIB}  
   ${CUDA_curand_LIBRARY}
   ${OpenMP_Fortran_LIBRARIES}
   #    ${hiprand_LIBRARIES}
#    ${HIP_LIBRARIES}
   )

TARGET_LINK_LIBRARIES(${UppASD_EXE} PRIVATE
   ${UppASD_LIB}
   ${UPPASD_LINALG_LIBRARIES}
   ${CMAKE_THREAD_LIBS_INIT} 
   ${OVF_LIB} 
   ${FFTW_LIB}  
   ${CUDA_curand_LIBRARY}
   ${OpenMP_Fortran_LIBRARIES}
   #    ${hiprand_LIBRARIES}
#    ${HIP_LIBRARIES}
   )


##############################################################################################
# Tell where to install this executable
##############################################################################################
IF(WIN32)
   SET(CMAKE_INSTALL_PREFIX "C:\\Program Files")
ELSE()
   SET(CMAKE_INSTALL_PREFIX /usr/local)
ENDIF(WIN32)


SET_TARGET_PROPERTIES(${UppASD_EXE} PROPERTIES LINKER_LANGUAGE "Fortran")
if (NOT BUILD_PYTHON)
   INSTALL(TARGETS ${UppASD_EXE} RUNTIME DESTINATION bin)
   INSTALL(TARGETS ${UppASD_LIB} RUNTIME DESTINATION lib)
   #INSTALL(TARGETS UppASD_F90w RUNTIME DESTINATION bin)
endif(NOT BUILD_PYTHON)
##############################################################################################

##############################################################################################
# Copy static library to ./lib
##############################################################################################
add_custom_command(
   TARGET ${UppASD_LIB} POST_BUILD
   COMMAND ${CMAKE_COMMAND} -E copy
   ${CMAKE_BINARY_DIR}/lib${UppASD_LIB}.a
   ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}
   )

#################################################################################################
#### Copy binary to ./source/sd for legacy purposes
#### currently disabled..
#################################################################################################
###add_custom_command(
###   TARGET ${UppASD_EXE} POST_BUILD
###   COMMAND ${CMAKE_COMMAND} -E copy
###   ${CMAKE_SOURCE_DIR}/bin/${UppASD_EXE}
###   ${LEGACY_OUTPUT_DIRECTORY}/sd
###   )


###############################################################################################
## Add Sources in sub-directories
###############################################################################################
add_subdirectory(${SRC})
##############################################################################################

###############################################################################################
## Add tests directory
###############################################################################################
set(TestBinary ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/${UppASD_EXE})
set(TestDIR ${CMAKE_SOURCE_DIR}/tests)
if (RUN_REG_TESTS)
add_test(NAME regression-test  COMMAND python3 -u bergtest.py --file regressionResaro.yaml --binary ${TestBinary}  WORKING_DIRECTORY ${TestDIR})
endif (RUN_REG_TESTS)
if (RUN_ASD_TESTS)
add_test(NAME asd-tests  COMMAND python3 -u bergtest.py --file regulartests.yaml --binary ${TestBinary}  WORKING_DIRECTORY ${TestDIR})
set_property(TEST asd-tests PROPERTY ENVIRONMENT "OMP_NUM_THREADS=1")
endif (RUN_ASD_TESTS)
if (RUN_SLD_TESTS)
add_test(NAME sld-tests  COMMAND python3 -u ./bergtest.py --file sldtests.yaml --binary ${TestBinary}  WORKING_DIRECTORY ${TestDIR})
endif (RUN_SLD_TESTS)
if (RUN_GNEB_TESTS)
add_test(NAME gneb-tests  COMMAND python3 -u ./bergtest.py --file regressionGNEB_15d7169b.yaml --binary ${TestBinary}  WORKING_DIRECTORY ${TestDIR})
endif (RUN_GNEB_TESTS)
#######hip here????
if (RUN_CUDA_TESTS)
add_test(NAME cuda-tests  COMMAND python3 -u bergtest.py --file cudatests.yaml --binary ${TestBinary}  WORKING_DIRECTORY ${TestDIR})
endif (RUN_CUDA_TESTS)
if(BUILD_TESTING)
   # Canonical adaptive coarse-graining topology tests are deliberately
   # independent of the dipole numerical kernels and the UppASD executable.
   add_executable(block_topology_tests
      ${CMAKE_SOURCE_DIR}/source/CoarseGraining/blocktopology.f90
      ${CMAKE_SOURCE_DIR}/tests/coarse_graining/test_block_topology.f90)
   set_target_properties(block_topology_tests PROPERTIES
      Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/modules/block_topology_tests)
   add_test(NAME coarse-graining-block-topology COMMAND block_topology_tests)

   # CG-03 material extraction fixtures exercise the production stiffness
   # module while keeping all calculation and validation I/O inside the test.
   add_executable(stiffness_material_tests
      ${CMAKE_SOURCE_DIR}/tests/coarse_graining/test_stiffness_material.f90)
   set_target_properties(stiffness_material_tests PROPERTIES
      Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/modules/stiffness_material_tests)
   target_link_libraries(stiffness_material_tests PRIVATE ${UppASD_LIB})
   add_test(NAME coarse-graining-stiffness-material COMMAND stiffness_material_tests)

   # RCG-02 analytic source-of-truth for the directed-pair DMI energy and
   # both dimer fields.  It calls the independent CPU field entrypoints.
   add_executable(dmi_dimer_energy_tests
      ${CMAKE_SOURCE_DIR}/source/Hamiltonian/applyhamiltonian.f90
      ${CMAKE_SOURCE_DIR}/tests/coarse_graining/test_dmi_dimer_energy.f90)
   set_target_properties(dmi_dimer_energy_tests PROPERTIES
      Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/modules/dmi_dimer_energy_tests)
   target_link_libraries(dmi_dimer_energy_tests PRIVATE ${UppASD_LIB})
   add_test(NAME coarse-graining-dmi-dimer-energy COMMAND dmi_dimer_energy_tests)

   # RCG-09A.1: exact adaptive-setup predicates for reciprocal directed
   # exchange/DMI lists, periodic-image aliasing, and ensemble moments.
   add_executable(adaptive_hamiltonian_contract_tests
      ${CMAKE_SOURCE_DIR}/tests/coarse_graining/test_adaptive_hamiltonian_contract.f90)
   set_target_properties(adaptive_hamiltonian_contract_tests PROPERTIES
      Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/modules/adaptive_hamiltonian_contract_tests)
   target_link_libraries(adaptive_hamiltonian_contract_tests PRIVATE ${UppASD_LIB})
   add_test(NAME coarse-graining-adaptive-hamiltonian-contract
      COMMAND adaptive_hamiltonian_contract_tests)

   # CG-04 CPU reference tensor operator: explicit energy/field fixtures on
   # cubic and skew canonical grids, independent of adaptive and GPU paths.
   add_executable(coarse_tensor_operator_tests
      ${CMAKE_SOURCE_DIR}/tests/coarse_graining/test_coarse_tensor_operator.f90)
   set_target_properties(coarse_tensor_operator_tests PROPERTIES
      Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/modules/coarse_tensor_operator_tests)
   target_link_libraries(coarse_tensor_operator_tests PRIVATE ${UppASD_LIB})
   add_test(NAME coarse-graining-tensor-operator COMMAND coarse_tensor_operator_tests)

   # CG-11 controlled two-sublattice CPU reference.  This target intentionally
   # has no GPU dependency: GPU enablement remains blocked on this diagnostic
   # gate and an independent physics review.
   add_executable(multichannel_coarse_tensor_operator_tests
      ${CMAKE_SOURCE_DIR}/tests/coarse_graining/test_multichannel_coarse_tensor_operator.f90)
   set_target_properties(multichannel_coarse_tensor_operator_tests PROPERTIES
      Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/modules/multichannel_coarse_tensor_operator_tests)
   target_link_libraries(multichannel_coarse_tensor_operator_tests PRIVATE ${UppASD_LIB})
   add_test(NAME coarse-graining-multichannel-tensor-operator
      COMMAND multichannel_coarse_tensor_operator_tests)

   # CG-05 CPU smooth projection: normalized trilinear prolongation, exact
   # adjoint restriction, read-only atomistic bilinear bonds, and spiral checks.
   add_executable(smooth_projected_operator_tests
      ${CMAKE_SOURCE_DIR}/tests/coarse_graining/test_smooth_projected_operator.f90)
   set_target_properties(smooth_projected_operator_tests PROPERTIES
      Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/modules/smooth_projected_operator_tests)
   target_link_libraries(smooth_projected_operator_tests PRIVATE ${UppASD_LIB})
   add_test(NAME coarse-graining-smooth-projected-operator
      COMMAND smooth_projected_operator_tests)

   # CG-06 static CPU hybrid: limiting masks, conservative buffer ownership,
   # normalized ghost reactions, patch/texture tests, and mesh refinement.
   add_executable(static_hybrid_operator_tests
      ${CMAKE_SOURCE_DIR}/tests/coarse_graining/test_static_hybrid_operator.f90)
   set_target_properties(static_hybrid_operator_tests PROPERTIES
      Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/modules/static_hybrid_operator_tests)
   target_link_libraries(static_hybrid_operator_tests PRIVATE ${UppASD_LIB})
   add_test(NAME coarse-graining-static-hybrid-operator
      COMMAND static_hybrid_operator_tests)

   # RCG-07 (F-08): local directional-stencil dilation correctness (brute-
   # force reference cross-check) and a wall-clock scaling negative control
   # discriminating the fixed local algorithm from the previous
   # all-block-by-all-seed quadratic scan.
   add_executable(static_hybrid_dilation_scaling_tests
      ${CMAKE_SOURCE_DIR}/tests/coarse_graining/test_static_hybrid_dilation_scaling.f90)
   set_target_properties(static_hybrid_dilation_scaling_tests PROPERTIES
      Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/modules/static_hybrid_dilation_scaling_tests)
   target_link_libraries(static_hybrid_dilation_scaling_tests PRIVATE ${UppASD_LIB})
   add_test(NAME coarse-graining-static-hybrid-dilation-scaling
      COMMAND static_hybrid_dilation_scaling_tests)

   # CG-07 is deliberately independent of field operators and reconstruction:
   # criteria only produce requests, while the explicit state-machine sync
   # point is solely responsible for changing block resolution.
   add_executable(block_selector_tests
      ${CMAKE_SOURCE_DIR}/source/CoarseGraining/blockselector.f90
      ${CMAKE_SOURCE_DIR}/tests/coarse_graining/test_block_selector.f90)
   set_target_properties(block_selector_tests PROPERTIES
      Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/modules/block_selector_tests)
   add_test(NAME coarse-graining-block-selector COMMAND block_selector_tests)

   # CG-08 transactional adaptive CPU orchestration: channel restriction,
   # deterministic reconstruction, energy-gated rollback, safe compaction,
   # complete-step timing, and representative texture transitions.
   add_executable(adaptive_hybrid_solver_tests
      ${CMAKE_SOURCE_DIR}/tests/coarse_graining/test_adaptive_hybrid_solver.f90)
   set_target_properties(adaptive_hybrid_solver_tests PROPERTIES
      Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/modules/adaptive_hybrid_solver_tests)
   target_link_libraries(adaptive_hybrid_solver_tests PRIVATE ${UppASD_LIB})
   add_test(NAME coarse-graining-adaptive-hybrid-solver
      COMMAND adaptive_hybrid_solver_tests)

   # RCG-03 polarization safety gate (F-14): a single-channel block whose
   # channel-resolved resultant/moment-sum ratio is below the accepted
   # threshold, or whose resultant is undefined (near-zero, never
   # normalized), must never be judged safe to coarsen.
   add_executable(polarization_gate_tests
      ${CMAKE_SOURCE_DIR}/tests/coarse_graining/test_polarization_gate.f90)
   set_target_properties(polarization_gate_tests PROPERTIES
      Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/modules/polarization_gate_tests)
   target_link_libraries(polarization_gate_tests PRIVATE ${UppASD_LIB})
   add_test(NAME coarse-graining-polarization-gate
      COMMAND polarization_gate_tests)

   # RCG-04B deterministic moving-state generators (conical spiral, +-q
   # chiral partner, periodic domain-wall pair). Host-only Python, no
   # UppASD/FFT/CUDA/HIP dependency -- generator-level unit tests only, no
   # production dynamics is run or claimed (see
   # docs/RCG-04_MOVING_E2E_EVIDENCE.md, RCG-04B section).
   add_test(NAME coarse-graining-moving-state-generator
      COMMAND python3
         ${CMAKE_SOURCE_DIR}/tests/coarse_graining/test_moving_state_generator.py -v)

   # RCG-04C reusable trajectory parsing/validation/comparison infrastructure
   # (strict moment/restart-file parser, component/angular error, spin
   # displacement, named energy/field series comparison, adaptive transition
   # and resolution-state history parsing, conical-mode amplitude/phase/
   # frequency, signed chirality, domain-wall centre/crossing tracking).
   # Host-only Python, no UppASD/FFT/CUDA/HIP dependency -- parser/metric
   # unit tests against synthetic records only; see
   # docs/RCG-04_MOVING_E2E_EVIDENCE.md (RCG-04C section) for the separate
   # manual production-output parsing evidence.
   add_test(NAME coarse-graining-trajectory-evidence
      COMMAND python3
         ${CMAKE_SOURCE_DIR}/tests/coarse_graining/test_trajectory_evidence.py -v)

   # RCG-04D independent initial effective-field/torque oracle, calibrated
   # against two real production runs (see torque_oracle.py's module
   # docstring) but not deriving the RCG-04D moving-state acceptance gate
   # from the diagnostic it gates. Host-only Python, no UppASD/CUDA/HIP
   # dependency; see docs/RCG-04_MOVING_E2E_EVIDENCE.md (RCG-04D section).
   add_test(NAME coarse-graining-torque-oracle
      COMMAND python3
         ${CMAKE_SOURCE_DIR}/tests/coarse_graining/test_torque_oracle.py -v)

   # RCG-04E independent structural nontriviality oracle for the all-coarse
   # block state (uncalibrated against production's regularized stiffness
   # fit -- see coarse_torque_oracle.py's module docstring for scope). Used
   # only as a pre-acceptance gate; the RCG-04E accuracy oracle is the
   # accepted all-fine atomistic trajectory. Host-only Python, no
   # UppASD/CUDA/HIP dependency.
   add_test(NAME coarse-graining-coarse-torque-oracle
      COMMAND python3
         ${CMAKE_SOURCE_DIR}/tests/coarse_graining/test_coarse_torque_oracle.py -v)

   # RCG-04F independent static fine/interface(buffer)/coarse block-ownership
   # oracle, re-deriving StaticHybridOperator's periodic dilation algorithm
   # from source (never reading production's own ownership arrays back).
   # Host-only Python, no UppASD/CUDA/HIP dependency.
   add_test(NAME coarse-graining-static-topology-oracle
      COMMAND python3
         ${CMAKE_SOURCE_DIR}/tests/coarse_graining/test_static_topology_oracle.py -v)

   # RCG-05C independent unit tests for the CPU/GPU ownership-map comparator
   # itself (parsing, index-set comparison, self-consistency oracle
   # cross-checks, bond coverage) against real captured stdout excerpts.
   # Host-only Python, no UppASD/CUDA/HIP dependency; the e2e comparator run
   # against the real executable is registered separately below
   # (adaptive-cg-ownership-map-comparator, GPU-configured builds only).
   add_test(NAME coarse-graining-ownership-map-comparator
      COMMAND python3
         ${CMAKE_SOURCE_DIR}/tests/coarse_graining/test_ownership_map_comparator.py -v)

   # RCG-06A MEM-LARGE-HOST (F-13/F-17): a minimal standalone reproduction of
   # the automatic-array pattern eliminated from
   # evaluate_static_hybrid_operator, deliberately NOT linked against
   # asdlib -- it must remain valid fault-injection evidence independent of
   # further production source changes. Driven only by
   # run_mem_large_host.py, never registered as its own passing ctest target
   # (it is expected to crash under the driver's constrained ulimit).
   add_executable(stack_overflow_fault_injection_tests
      ${CMAKE_SOURCE_DIR}/tests/coarse_graining/test_stack_overflow_fault_injection.f90)
   set_target_properties(stack_overflow_fault_injection_tests PROPERTIES
      Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/modules/stack_overflow_fault_injection_tests)

   # RCG-06A MEM-LARGE-HOST: the ordinary UppASD executable on a Natom=8000
   # adaptive-CG host under a constrained process stack limit, paired with
   # the fault-injection negative control above. See
   # tests/coarse_graining/e2e/mem_large_host/README.md.
   add_test(NAME adaptive-cg-mem-large-host
      COMMAND python3
         ${CMAKE_SOURCE_DIR}/tests/coarse_graining/run_mem_large_host.py
         --binary $<TARGET_FILE:${UppASD_EXE}>
         --fault-injection-binary $<TARGET_FILE:stack_overflow_fault_injection_tests>)
   set_tests_properties(adaptive-cg-mem-large-host
      PROPERTIES LABELS "coarse-graining;cg13;cg13-cpu"
                 ENVIRONMENT "OMP_NUM_THREADS=1")

   # RCG-06D RECONSTRUCTION-RNG-SPATIAL-STATS (F-20): standalone driver that
   # emits raw draws from the literal tuple-seeded MINSTD reconstruction RNG
   # (deterministic_reconstruction_seed + next_uniform, AdaptiveHybridSolver)
   # over a sweep of block indices, so run_reconstruction_rng_spatial_stats.py
   # can measure the spatial-correlation statistics F-20 asks for. See
   # docs/RCG-06_MEMORY_TIMING_PRECISION_EVIDENCE.md (RCG-06D section) for
   # the measured numbers and the zero-cone-vs-nonzero-cone scope decision.
   add_executable(reconstruction_rng_spatial_stats_tests
      ${CMAKE_SOURCE_DIR}/tests/coarse_graining/test_reconstruction_rng_spatial_stats.f90)
   set_target_properties(reconstruction_rng_spatial_stats_tests PROPERTIES
      Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/modules/reconstruction_rng_spatial_stats_tests)
   target_link_libraries(reconstruction_rng_spatial_stats_tests PRIVATE ${UppASD_LIB})

   set(RCG06D_GPU_BINARY_ARGS)
   if(USE_CUDA OR USE_HIP)
      list(APPEND RCG06D_GPU_BINARY_ARGS
         --gpu-binary $<TARGET_FILE:reconstruction_rng_gpu_parity_tests>)
   endif()
   add_test(NAME adaptive-cg-reconstruction-rng-spatial-stats
      COMMAND python3
         ${CMAKE_SOURCE_DIR}/tests/coarse_graining/run_reconstruction_rng_spatial_stats.py
         --cpu-binary $<TARGET_FILE:reconstruction_rng_spatial_stats_tests>
         ${RCG06D_GPU_BINARY_ARGS})
   set_tests_properties(adaptive-cg-reconstruction-rng-spatial-stats
      PROPERTIES LABELS "coarse-graining;cg13;cg13-cpu")
   if(USE_CUDA)
      set_property(TEST adaptive-cg-reconstruction-rng-spatial-stats APPEND PROPERTY LABELS cg13-cuda)
   elseif(USE_HIP)
      set_property(TEST adaptive-cg-reconstruction-rng-spatial-stats APPEND PROPERTY LABELS cg13-hip)
   endif()

   # CG-13 labels are the repeatable suite interface.  Keep the focused
   # analytic/reference tests in the CPU suite even in a GPU build; the
   # production and backend tests below add the matching CUDA/HIP label.
   set_tests_properties(
      coarse-graining-block-topology
      coarse-graining-dmi-dimer-energy
      coarse-graining-stiffness-material
      coarse-graining-tensor-operator
      coarse-graining-multichannel-tensor-operator
      coarse-graining-smooth-projected-operator
      coarse-graining-static-hybrid-operator
      coarse-graining-static-hybrid-dilation-scaling
      coarse-graining-block-selector
      coarse-graining-adaptive-hybrid-solver
      coarse-graining-polarization-gate
      coarse-graining-moving-state-generator
      coarse-graining-trajectory-evidence
      coarse-graining-torque-oracle
      coarse-graining-coarse-torque-oracle
      coarse-graining-static-topology-oracle
      coarse-graining-ownership-map-comparator
      PROPERTIES LABELS "coarse-graining;cg13;cg13-cpu;reference")
   if(USE_CUDA)
      set_property(TEST
         coarse-graining-block-topology
         coarse-graining-stiffness-material
         coarse-graining-tensor-operator
         coarse-graining-multichannel-tensor-operator
         coarse-graining-smooth-projected-operator
         coarse-graining-static-hybrid-operator
         coarse-graining-static-hybrid-dilation-scaling
         coarse-graining-block-selector
         coarse-graining-adaptive-hybrid-solver
         coarse-graining-polarization-gate
         coarse-graining-moving-state-generator
         coarse-graining-trajectory-evidence
         coarse-graining-torque-oracle
         coarse-graining-coarse-torque-oracle
         coarse-graining-static-topology-oracle
         coarse-graining-ownership-map-comparator
         APPEND PROPERTY LABELS cg13-cuda)
   elseif(USE_HIP)
      set_property(TEST
         coarse-graining-block-topology
         coarse-graining-stiffness-material
         coarse-graining-tensor-operator
         coarse-graining-multichannel-tensor-operator
         coarse-graining-smooth-projected-operator
         coarse-graining-static-hybrid-operator
         coarse-graining-static-hybrid-dilation-scaling
         coarse-graining-block-selector
         coarse-graining-adaptive-hybrid-solver
         coarse-graining-polarization-gate
         coarse-graining-moving-state-generator
         coarse-graining-trajectory-evidence
         coarse-graining-torque-oracle
         coarse-graining-coarse-torque-oracle
         coarse-graining-static-topology-oracle
         coarse-graining-ownership-map-comparator
         APPEND PROPERTY LABELS cg13-hip)
   endif()

   # CG-10.5 must be exercised through the ordinary executable and ordinary
   # inpsd.dat files; focused operator calls are insufficient evidence.
   set(ADAPTIVE_CG_GPU_TEST_ARGS)
   if(USE_CUDA OR USE_HIP)
      list(APPEND ADAPTIVE_CG_GPU_TEST_ARGS
         --gpu-binary $<TARGET_FILE:${UppASD_EXE}>)
   endif()
   add_test(NAME adaptive-cg-production-e2e
      COMMAND python3
         ${CMAKE_SOURCE_DIR}/tests/coarse_graining/run_production_e2e.py
         --binary $<TARGET_FILE:${UppASD_EXE}>
         ${ADAPTIVE_CG_GPU_TEST_ARGS})
   add_test(NAME adaptive-cg-setup-rejection-matrix
      COMMAND python3
         ${CMAKE_SOURCE_DIR}/tests/coarse_graining/run_setup_rejection_matrix.py
         --binary $<TARGET_FILE:${UppASD_EXE}>)
   # RCG-04D E2E-MOVING-OFF-FINE: genuinely nonstationary feature-off/all-fine
   # parity, independent-oracle nontriviality gate, and defect-sensitivity
   # negative controls. CPU-only in this slice; CPU/GPU backend parity for
   # this fixture is RCG-04I's responsibility, not claimed here.
   add_test(NAME adaptive-cg-moving-off-fine
      COMMAND python3
         ${CMAKE_SOURCE_DIR}/tests/coarse_graining/run_moving_off_fine.py
         --binary $<TARGET_FILE:${UppASD_EXE}>)
   set_tests_properties(adaptive-cg-moving-off-fine
      PROPERTIES LABELS "coarse-graining;cg13;cg13-cpu;moving-parity"
                 ENVIRONMENT "OMP_NUM_THREADS=1")
   # RCG-04E E2E-MOVING-ALL-COARSE: all-coarse long-wavelength dynamics at
   # multiple physical block sizes, validated against the accepted all-fine
   # atomistic reference (re-verified at this fixture's own wide geometry),
   # plus independent nontriviality gates and defect-sensitivity negative
   # controls. CPU-only in this slice; CPU/GPU backend parity is RCG-04I's
   # responsibility, not claimed here.
   add_test(NAME adaptive-cg-moving-all-coarse
      COMMAND python3
         ${CMAKE_SOURCE_DIR}/tests/coarse_graining/run_moving_all_coarse.py
         --binary $<TARGET_FILE:${UppASD_EXE}>)
   set_tests_properties(adaptive-cg-moving-all-coarse
      PROPERTIES LABELS "coarse-graining;cg13;cg13-cpu;moving-parity"
                 ENVIRONMENT "OMP_NUM_THREADS=1")
   # RCG-04F E2E-MOVING-STATIC: static fine/interface(buffer)/coarse mixed-mask
   # decomposition of the RCG-04E wide conical spiral, with independent
   # ownership-topology verification, a fixed-physical-partition spatial-
   # refinement pair, and shifted-interface sensitivity. CPU-only in this
   # slice; CPU/GPU backend parity is RCG-04I's responsibility, not claimed
   # here.
   add_test(NAME adaptive-cg-moving-static-mixed
      COMMAND python3
         ${CMAKE_SOURCE_DIR}/tests/coarse_graining/run_moving_static_mixed.py
         --binary $<TARGET_FILE:${UppASD_EXE}>)
   set_tests_properties(adaptive-cg-moving-static-mixed
      PROPERTIES LABELS "coarse-graining;cg13;cg13-cpu;moving-parity"
                 ENVIRONMENT "OMP_NUM_THREADS=1")
   # RCG-04G E2E-MOVING-ADAPTIVE: accepted adaptive-resolution transitions
   # during genuine domain-wall motion (a deterministic periodic Neel
   # wall-pair, RCG-04B's domain_wall_pair_state), with independent
   # nontriviality gates, wall-centre/block-boundary-crossing tracking, full
   # transition-history and RCG-03 polarization-safety assertions, and
   # defect-sensitivity negative controls. CPU-only in this slice; CPU/GPU
   # backend parity is RCG-04I's responsibility, not claimed here.
   add_test(NAME adaptive-cg-moving-adaptive-wall
      COMMAND python3
         ${CMAKE_SOURCE_DIR}/tests/coarse_graining/run_moving_adaptive_wall.py
         --binary $<TARGET_FILE:${UppASD_EXE}>)
   set_tests_properties(adaptive-cg-moving-adaptive-wall
      PROPERTIES LABELS "coarse-graining;cg13;cg13-cpu;moving-parity"
                 ENVIRONMENT "OMP_NUM_THREADS=1")
   # RCG-04H DMI-HYBRID-CROSSING: RCG-04B's deterministic +q/-q chiral-partner
   # pair on the RCG-04F fine/interface/coarse partition, with DMI and
   # uniaxial anisotropy both enabled, asserting signed chirality, +q/-q DMI
   # energy ordering, named DMI/anisotropy energy series, a complete
   # mixed-vs-all-fine trajectory comparison, a spatial-refinement pair, and
   # DMI-specific interface ownership -- plus a sign-reversed-DMI negative
   # control against a fixed handedness oracle anchored to RCG-02. CPU-only
   # in this slice; CPU/GPU backend parity is RCG-04I's responsibility, not
   # claimed here.
   add_test(NAME adaptive-cg-moving-dmi-chiral
      COMMAND python3
         ${CMAKE_SOURCE_DIR}/tests/coarse_graining/run_moving_dmi_chiral.py
         --binary $<TARGET_FILE:${UppASD_EXE}>)
   set_tests_properties(adaptive-cg-moving-dmi-chiral
      PROPERTIES LABELS "coarse-graining;cg13;cg13-cpu;moving-parity"
                 ENVIRONMENT "OMP_NUM_THREADS=1")
   # RCG-04I: CPU/GPU backend precision parity for every RCG-04D-H accepted
   # moving fixture (this build's own executable run twice per fixture --
   # once as-is for the CPU path, once with an appended do_gpu/do_gpu_llg
   # block for the GPU path -- see run_moving_backend_parity.py's module
   # docstring for why one CUDA/HIP-enabled binary suffices for both legs).
   # Only meaningful, and only registered, in a GPU-enabled configuration;
   # a plain CPU build has no second backend to compare against and does
   # not register this test (matching adaptive-cg-production-e2e's
   # existing --gpu-binary-only-if-GPU-configured pattern).
   set(RCG04I_GPU_BINARY_ARGS)
   if(USE_CUDA)
      set(RCG04I_BACKEND_LABEL cg13-cuda)
      if(UPPASD_PRECISION STREQUAL "SINGLE")
         list(APPEND RCG04I_GPU_BINARY_ARGS --cuda-fp32-binary $<TARGET_FILE:${UppASD_EXE}>)
      else()
         list(APPEND RCG04I_GPU_BINARY_ARGS --cuda-fp64-binary $<TARGET_FILE:${UppASD_EXE}>)
      endif()
   elseif(USE_HIP)
      set(RCG04I_BACKEND_LABEL cg13-hip)
      if(UPPASD_PRECISION STREQUAL "SINGLE")
         list(APPEND RCG04I_GPU_BINARY_ARGS --hip-fp32-binary $<TARGET_FILE:${UppASD_EXE}>)
      else()
         list(APPEND RCG04I_GPU_BINARY_ARGS --hip-fp64-binary $<TARGET_FILE:${UppASD_EXE}>)
      endif()
   endif()
   if(USE_CUDA OR USE_HIP)
      add_test(NAME adaptive-cg-moving-backend-parity
         COMMAND python3
            ${CMAKE_SOURCE_DIR}/tests/coarse_graining/run_moving_backend_parity.py
            --cpu-binary $<TARGET_FILE:${UppASD_EXE}>
            ${RCG04I_GPU_BINARY_ARGS}
            --workspace-root ${CMAKE_CURRENT_BINARY_DIR}/rcg04i_backend_parity_workspace)
      set_tests_properties(adaptive-cg-moving-backend-parity
         PROPERTIES LABELS "coarse-graining;cg13;${RCG04I_BACKEND_LABEL};moving-parity;backend-parity"
                    ENVIRONMENT "OMP_NUM_THREADS=1")
   endif()
   # RCG-06C (F-18): CPU/GPU adaptive-CG phase-timing reconciliation.  Reuses
   # RCG-04I's fixture list, workspace preparation, and single-fixture
   # runner (run_moving_backend_parity.py) rather than a new mechanism.  Runs
   # moving_wall_adaptive (the only tracked fixture exercising every phase,
   # including the selector/transition path) and asserts the printed
   # wall/phase_sum/unaccounted reconciliation lines are internally
   # consistent and phase_sum never materially exceeds wall time.  The GPU
   # leg only runs, and is only registered, in a GPU-enabled configuration
   # (the same UppASD_EXE binary runs both legs, matching RCG-04I's
   # mechanism); a plain CPU build still runs and asserts the CPU leg.
   set(RCG06C_GPU_BINARY_ARGS)
   if(USE_CUDA OR USE_HIP)
      list(APPEND RCG06C_GPU_BINARY_ARGS --gpu-binary $<TARGET_FILE:${UppASD_EXE}>)
   endif()
   add_test(NAME adaptive-cg-timing-reconciliation
      COMMAND python3
         ${CMAKE_SOURCE_DIR}/tests/coarse_graining/run_timing_reconciliation.py
         --cpu-binary $<TARGET_FILE:${UppASD_EXE}>
         ${RCG06C_GPU_BINARY_ARGS}
         --workspace-root ${CMAKE_CURRENT_BINARY_DIR}/rcg06c_timing_reconciliation_workspace)
   set_tests_properties(adaptive-cg-timing-reconciliation
      PROPERTIES LABELS "coarse-graining;cg13;cg13-cpu;moving-parity"
                 ENVIRONMENT "OMP_NUM_THREADS=1")
   if(USE_CUDA)
      set_property(TEST adaptive-cg-timing-reconciliation APPEND PROPERTY LABELS cg13-cuda)
   elseif(USE_HIP)
      set_property(TEST adaptive-cg-timing-reconciliation APPEND PROPERTY LABELS cg13-hip)
   endif()
   # RCG-05C: CPU/GPU ownership-map comparator run against the real
   # executable -- every RCG-04I adaptive_cg=True moving fixture (regression,
   # must match exactly) plus the new e2e/ownership_aniso_buffer fixture
   # (defect demonstration: CPU matches the correct per-axis dilation
   # oracle, GPU matches the isotropic-collapsed one). Same
   # GPU-binary-only-if-GPU-configured pattern as RCG-04I immediately above.
   set(RCG05C_GPU_BINARY_ARGS)
   if(USE_CUDA)
      set(RCG05C_BACKEND_LABEL cg13-cuda)
      if(UPPASD_PRECISION STREQUAL "SINGLE")
         list(APPEND RCG05C_GPU_BINARY_ARGS --cuda-fp32-binary $<TARGET_FILE:${UppASD_EXE}>)
      else()
         list(APPEND RCG05C_GPU_BINARY_ARGS --cuda-fp64-binary $<TARGET_FILE:${UppASD_EXE}>)
      endif()
   elseif(USE_HIP)
      set(RCG05C_BACKEND_LABEL cg13-hip)
      if(UPPASD_PRECISION STREQUAL "SINGLE")
         list(APPEND RCG05C_GPU_BINARY_ARGS --hip-fp32-binary $<TARGET_FILE:${UppASD_EXE}>)
      else()
         list(APPEND RCG05C_GPU_BINARY_ARGS --hip-fp64-binary $<TARGET_FILE:${UppASD_EXE}>)
      endif()
   endif()
   # RCG-05D: the buffer-width scalarization defect RCG-05C demonstrated
   # above is now fixed (per-axis dilation restored end to end). GPU's own
   # reported map, re-dilated from GPU's own FINE seed set, must match the
   # correct per-axis oracle and no longer the isotropic-collapsed one --
   # --expect-buffer-width-fixed isolates that claim from a separate,
   # already-documented seed-set gap (GPU's hardAtomisticBlockMask is
   # sourced only from the polarization gate, not cg_static_mask_file) that
   # keeps the full CPU/GPU maps from agreeing exactly on this fixture and
   # is out of RCG-05D's authorized scope. This registration is also the
   # permanent regression against the buffer-width defect being
   # reintroduced: if bufferDilationBlocks were ever scalarized again, GPU
   # would once more match the isotropic oracle and this assertion would
   # fail.
   if(USE_CUDA OR USE_HIP)
      add_test(NAME adaptive-cg-ownership-map-comparator
         COMMAND python3
            ${CMAKE_SOURCE_DIR}/tests/coarse_graining/run_ownership_map_comparator.py
            --cpu-binary $<TARGET_FILE:${UppASD_EXE}>
            ${RCG05C_GPU_BINARY_ARGS}
            --expect-buffer-width-fixed
            --workspace-root ${CMAKE_CURRENT_BINARY_DIR}/rcg05c_ownership_comparator_workspace)
      set_tests_properties(adaptive-cg-ownership-map-comparator
         PROPERTIES LABELS "coarse-graining;cg13;${RCG05C_BACKEND_LABEL};moving-parity;backend-parity;ownership"
                    ENVIRONMENT "OMP_NUM_THREADS=1")
   endif()
   # RCG-05F: extends RCG-05C's ownership-map comparator (identical
   # --cpu-binary/GPU-binary-args pattern, GPU-configured builds only, since
   # this exercises gpu_dipole_mode EWALD3D_FFT which requires do_gpu=Y) to
   # confirm the regular-grid FFT dipole does not perturb fine/buffer/coarse
   # ownership: a genuinely anisotropic RCG-05B-generator fixture
   # (ownership_dipole_unequal_width) plus the pre-existing RCG-04-era
   # gpu_fft_static_mixed/parity_static_cpu dipole pair, both compared by
   # RCG-05C's full per-block identity comparator rather than an aggregate
   # count.
   if(USE_CUDA OR USE_HIP)
      add_test(NAME adaptive-cg-dipole-ownership-check
         COMMAND python3
            ${CMAKE_SOURCE_DIR}/tests/coarse_graining/run_dipole_ownership_check.py
            --cpu-binary $<TARGET_FILE:${UppASD_EXE}>
            ${RCG05C_GPU_BINARY_ARGS}
            --workspace-root ${CMAKE_CURRENT_BINARY_DIR}/rcg05f_dipole_ownership_workspace)
      set_tests_properties(adaptive-cg-dipole-ownership-check
         PROPERTIES LABELS "coarse-graining;cg13;${RCG05C_BACKEND_LABEL};backend-parity;ownership"
                    ENVIRONMENT "OMP_NUM_THREADS=1")
   endif()
   add_test(NAME adaptive-cg-fixture-dependencies
      COMMAND python3
         ${CMAKE_SOURCE_DIR}/tests/coarse_graining/audit_fixture_dependencies.py)
   set_tests_properties(adaptive-cg-fixture-dependencies
      PROPERTIES LABELS "coarse-graining;cg13;packaging")
   set_tests_properties(adaptive-cg-setup-rejection-matrix
      PROPERTIES LABELS "coarse-graining;cg13;cg13-cpu;rejection"
                 ENVIRONMENT "OMP_NUM_THREADS=1")
   if(USE_CUDA)
      set(CG13_PRODUCTION_LABEL cg13-cuda)
   elseif(USE_HIP)
      set(CG13_PRODUCTION_LABEL cg13-hip)
   else()
      set(CG13_PRODUCTION_LABEL cg13-cpu)
   endif()
   set_tests_properties(adaptive-cg-setup-rejection-matrix
      PROPERTIES LABELS "coarse-graining;cg13;${CG13_PRODUCTION_LABEL};rejection")
   set_tests_properties(adaptive-cg-production-e2e
      PROPERTIES LABELS "coarse-graining;cg13;${CG13_PRODUCTION_LABEL};production"
                 ENVIRONMENT "OMP_NUM_THREADS=1")

   # RCG-05F: reuses RCG-04G's transition-history infrastructure
   # (trajectory_evidence.parse_transition_events/parse_resolution_state_history)
   # to confirm ownership invariants (a complete, well-formed fine/buffer/
   # coarse partition; no accepted transition silently dropped; no rejected
   # transition silently applied) hold across every recorded mask-rebuild
   # transition on ownership_aniso_buffer, not only at the final snapshot.
   # CPU-only: trajectory_evidence.py's own module docstring already
   # documents GPU stdout has no per-event transition log.
   add_test(NAME adaptive-cg-transition-ownership-invariants
      COMMAND python3
         ${CMAKE_SOURCE_DIR}/tests/coarse_graining/check_transition_ownership_invariants.py
         --cpu-binary $<TARGET_FILE:${UppASD_EXE}>
         --workspace-root ${CMAKE_CURRENT_BINARY_DIR}/rcg05f_transition_invariants_workspace)
   set_tests_properties(adaptive-cg-transition-ownership-invariants
      PROPERTIES LABELS "coarse-graining;cg13;${CG13_PRODUCTION_LABEL};ownership"
                 ENVIRONMENT "OMP_NUM_THREADS=1")

   # Luna's host-only OPEN_FFT oracle has no UppASD, FFT, CUDA, or HIP
   # dependency.  Keep it runnable before the production open builder exists.
   add_test(NAME dipole-open-fft-oracle
      COMMAND python3 ${CMAKE_SOURCE_DIR}/tests/dipole_validation/test_open_fft_oracle.py -v)

   # Backend-neutral host-fp64 OPEN_FFT builder tests.  These targets construct
   # and apply only the finite padded tensor; they add no production selector,
   # Fortran bridge, FFT backend, CUDA/HIP context, or periodic expected value.
   add_executable(dipole_open_host_tests
      ${CMAKE_SOURCE_DIR}/source/gpu_files/dipoleOpenKernel.cpp
      ${CMAKE_SOURCE_DIR}/tests/dipole_validation/test_host_open_dipole_kernel.cpp)
   target_include_directories(dipole_open_host_tests PRIVATE ${CMAKE_SOURCE_DIR}/source/gpu_files)
   add_test(NAME dipole-open-host-builder COMMAND dipole_open_host_tests)

   add_executable(dipole_open_kernel_host_driver
      ${CMAKE_SOURCE_DIR}/source/gpu_files/dipoleOpenKernel.cpp
      ${CMAKE_SOURCE_DIR}/tests/dipole_validation/open_kernel_host_driver.cpp)
   target_include_directories(dipole_open_kernel_host_driver PRIVATE
      ${CMAKE_SOURCE_DIR}/source/gpu_files)
   add_test(NAME dipole-open-host-goldens
      COMMAND python3 ${CMAKE_SOURCE_DIR}/tests/dipole_validation/test_open_kernel_goldens.py
              --driver $<TARGET_FILE:dipole_open_kernel_host_driver>)

   add_executable(dipole_ewald_host_tests
      ${CMAKE_SOURCE_DIR}/source/gpu_files/dipoleEwaldKernel.cpp
      ${CMAKE_SOURCE_DIR}/tests/dipole_validation/test_host_periodic_ewald_kernel.cpp)
   target_include_directories(dipole_ewald_host_tests PRIVATE ${CMAKE_SOURCE_DIR}/source/gpu_files)
   add_test(NAME dipole-ewald-host-builder COMMAND dipole_ewald_host_tests)

   # This target is intentionally built by CMake with the active GPU backend;
   # do not recompile a driver from Python/nvcc.  It exercises the complete
   # Builder A tensor upload, batched kernel R2C, and stored-spectrum contract.
   if(USE_CUDA OR USE_HIP)
      # RCG-02 device half of DMI-DIMER-CPU-GPU.  The test runs the actual
      # hamdev::dm_field device helper, rather than reproducing it on host.
      add_executable(gpu_dmi_dimer_tests
         ${CMAKE_SOURCE_DIR}/tests/coarse_graining/test_gpu_dmi_dimer.cpp)
      target_include_directories(gpu_dmi_dimer_tests PRIVATE
         ${CMAKE_SOURCE_DIR}/source/gpu_files)
      if(USE_CUDA)
         set_source_files_properties(${CMAKE_SOURCE_DIR}/tests/coarse_graining/test_gpu_dmi_dimer.cpp
            PROPERTIES LANGUAGE CUDA)
         target_compile_definitions(gpu_dmi_dimer_tests PRIVATE CUDA_V)
         target_link_libraries(gpu_dmi_dimer_tests PRIVATE CUDA::cudart)
      else()
         set_source_files_properties(${CMAKE_SOURCE_DIR}/tests/coarse_graining/test_gpu_dmi_dimer.cpp
            PROPERTIES LANGUAGE HIP)
         target_compile_definitions(gpu_dmi_dimer_tests PRIVATE HIP_V)
      endif()
      add_test(NAME coarse-graining-gpu-dmi-dimer COMMAND gpu_dmi_dimer_tests)
      set_tests_properties(coarse-graining-gpu-dmi-dimer PROPERTIES SKIP_RETURN_CODE 77)
      if(USE_CUDA)
         set_tests_properties(coarse-graining-gpu-dmi-dimer
            PROPERTIES LABELS "coarse-graining;cg13;cg13-cuda;reference")
      else()
         set_tests_properties(coarse-graining-gpu-dmi-dimer
            PROPERTIES LABELS "coarse-graining;cg13;cg13-hip;reference")
      endif()

      # RCG-06B ENERGY-FP32-ACCUM (F-11): standalone microbenchmark isolating
      # the atomicAdd accumulation-order defect class from the full adaptive
      # CG physics pipeline. Deliberately NOT linked against asdlib, matching
      # RCG-06A's stack_overflow_fault_injection precedent. Shares the exact
      # atomicAddEnergyTerm CAS-loop helper (gpuAtomicDouble.hpp) that the
      # production kernels in gpuAdaptiveRuntime.cpp use.
      add_executable(energy_fp32_accum_tests
         ${CMAKE_SOURCE_DIR}/tests/coarse_graining/test_energy_fp32_accum.cpp)
      target_include_directories(energy_fp32_accum_tests PRIVATE
         ${CMAKE_SOURCE_DIR}/source/gpu_files)
      if(USE_CUDA)
         set_source_files_properties(${CMAKE_SOURCE_DIR}/tests/coarse_graining/test_energy_fp32_accum.cpp
            PROPERTIES LANGUAGE CUDA)
         target_compile_definitions(energy_fp32_accum_tests PRIVATE CUDA_V)
         target_link_libraries(energy_fp32_accum_tests PRIVATE CUDA::cudart)
      else()
         set_source_files_properties(${CMAKE_SOURCE_DIR}/tests/coarse_graining/test_energy_fp32_accum.cpp
            PROPERTIES LANGUAGE HIP)
         target_compile_definitions(energy_fp32_accum_tests PRIVATE HIP_V)
      endif()
      add_test(NAME adaptive-cg-energy-fp32-accum COMMAND energy_fp32_accum_tests)
      set_tests_properties(adaptive-cg-energy-fp32-accum PROPERTIES SKIP_RETURN_CODE 77)
      if(USE_CUDA)
         set_tests_properties(adaptive-cg-energy-fp32-accum
            PROPERTIES LABELS "coarse-graining;cg13;cg13-cuda")
      else()
         set_tests_properties(adaptive-cg-energy-fp32-accum
            PROPERTIES LABELS "coarse-graining;cg13;cg13-hip")
      endif()

      # RCG-06D RECONSTRUCTION-RNG-SPATIAL-STATS GPU companion (F-20):
      # standalone, deliberately NOT linked against asdlib, matching
      # RCG-06B's energy_fp32_accum precedent. Shares the literal device
      # generator (adaptiveTupleSeed/adaptiveNextUniform,
      # gpuAdaptiveReconstructionRng.hpp) that gpuAdaptiveRuntime.cpp's
      # production kernels use. Not registered as its own ctest target --
      # only used as --gpu-binary by
      # adaptive-cg-reconstruction-rng-spatial-stats (registered in the
      # host-side block above), matching stack_overflow_fault_injection_tests'
      # pairing pattern.
      add_executable(reconstruction_rng_gpu_parity_tests
         ${CMAKE_SOURCE_DIR}/tests/coarse_graining/test_reconstruction_rng_gpu_parity.cpp)
      target_include_directories(reconstruction_rng_gpu_parity_tests PRIVATE
         ${CMAKE_SOURCE_DIR}/source/gpu_files)
      if(USE_CUDA)
         set_source_files_properties(${CMAKE_SOURCE_DIR}/tests/coarse_graining/test_reconstruction_rng_gpu_parity.cpp
            PROPERTIES LANGUAGE CUDA)
         target_compile_definitions(reconstruction_rng_gpu_parity_tests PRIVATE CUDA_V)
         target_link_libraries(reconstruction_rng_gpu_parity_tests PRIVATE CUDA::cudart)
      else()
         set_source_files_properties(${CMAKE_SOURCE_DIR}/tests/coarse_graining/test_reconstruction_rng_gpu_parity.cpp
            PROPERTIES LANGUAGE HIP)
         target_compile_definitions(reconstruction_rng_gpu_parity_tests PRIVATE HIP_V)
      endif()
      if(UPPASD_PRECISION STREQUAL "SINGLE")
         target_compile_definitions(reconstruction_rng_gpu_parity_tests PRIVATE SINGLE_PREC)
      endif()

      # CG-09/CG-10 immutable topology transfer, mutable block/channel
      # storage, compact work scheduling, tensor/interface/adaptive parity,
      # measured phase boundaries, and cleanup accounting on the active
      # CUDA/HIP backend.
      add_executable(gpu_adaptive_runtime_tests
         ${CMAKE_SOURCE_DIR}/tests/coarse_graining/test_gpu_adaptive_runtime.cpp)
      target_include_directories(gpu_adaptive_runtime_tests PRIVATE
         ${CMAKE_SOURCE_DIR}/source/gpu_files
         ${CMAKE_SOURCE_DIR}/source/gpu_files/measurement)
      target_link_libraries(gpu_adaptive_runtime_tests PRIVATE ${UppASD_LIB})
      if(USE_CUDA)
         set_source_files_properties(${CMAKE_SOURCE_DIR}/tests/coarse_graining/test_gpu_adaptive_runtime.cpp
            PROPERTIES LANGUAGE CUDA)
         target_compile_definitions(gpu_adaptive_runtime_tests PRIVATE CUDA_V)
         target_link_libraries(gpu_adaptive_runtime_tests PRIVATE CUDA::cudart)
      else()
         set_source_files_properties(${CMAKE_SOURCE_DIR}/tests/coarse_graining/test_gpu_adaptive_runtime.cpp
            PROPERTIES LANGUAGE HIP)
         target_compile_definitions(gpu_adaptive_runtime_tests PRIVATE HIP_V)
      endif()
      if(UPPASD_PRECISION STREQUAL "SINGLE")
         target_compile_definitions(gpu_adaptive_runtime_tests PRIVATE SINGLE_PREC)
      endif()
      add_test(NAME coarse-graining-gpu-adaptive-runtime COMMAND gpu_adaptive_runtime_tests)
      if(USE_CUDA)
         set(CG13_GPU_LABEL cg13-cuda)
      else()
         set(CG13_GPU_LABEL cg13-hip)
      endif()
      set_tests_properties(coarse-graining-gpu-adaptive-runtime
         PROPERTIES LABELS "coarse-graining;cg13;${CG13_GPU_LABEL};reference")

      # RCG-05E: compute-sanitizer racecheck/memcheck/synccheck/initcheck
      # against dilateAdaptiveState specifically (--kernel-name), run through
      # gpu_adaptive_runtime_tests above, which already launches it with a
      # nonzero per-axis buffer width. Kept out of cg13-cuda/cg13-hip (a
      # dedicated "sanitizer" label instead): sanitized runs carry real
      # instrumentation overhead and an extra tool dependency neither of the
      # baseline regression labels should require. compute-sanitizer ships
      # alongside nvcc in every CUDA toolkit that can build USE_CUDA, so
      # find_program only fails to register this test on a CUDA build whose
      # toolkit installation was stripped of it. HIP's equivalent
      # (rocprof/rocgdb-based race detection) is not wired here: this
      # repository has no HIP toolchain to build or verify it against (see
      # RCG-04-FU1's still-open deferral); recorded as a deferral, not
      # attempted blind.
      find_program(UPPASD_COMPUTE_SANITIZER NAMES compute-sanitizer)
      if(USE_CUDA AND UPPASD_COMPUTE_SANITIZER)
         add_test(NAME adaptive-cg-dilation-sanitizer
            COMMAND python3
               ${CMAKE_SOURCE_DIR}/tests/gpu_regression/sanitize.py
               --dilation-kernel
               --binary $<TARGET_FILE:gpu_adaptive_runtime_tests>
               --sanitizer ${UPPASD_COMPUTE_SANITIZER}
               --workdir ${CMAKE_CURRENT_BINARY_DIR}/rcg05e_dilation_sanitizer_workspace)
         set_tests_properties(adaptive-cg-dilation-sanitizer
            PROPERTIES LABELS "coarse-graining;cg13;sanitizer")
      endif()

      # RCG-09 PERF-ATOMISTIC-PROD / PERF-CG-SWEEP.  Hardware-dependent
      # acceptance harness, intentionally not a CTest: it checks that the
      # adaptive and production atomistic paths evaluate the same field, then
      # times UppASD's real GpuHamiltonianCalculations + GpuDepondtIntegrator
      # step against the adaptive step over a fine-fraction sweep, reporting
      # active atom/block/interface/bond fractions, per-phase device times,
      # launch counts, selector/compaction overhead, memory, and raw samples.
      add_executable(gpu_adaptive_runtime_benchmark
         ${CMAKE_SOURCE_DIR}/tests/coarse_graining/benchmark_gpu_adaptive_runtime.cpp)
      target_include_directories(gpu_adaptive_runtime_benchmark PRIVATE
         ${CMAKE_SOURCE_DIR}/source/gpu_files
         ${CMAKE_SOURCE_DIR}/source/gpu_files/measurement)
      target_link_libraries(gpu_adaptive_runtime_benchmark PRIVATE ${UppASD_LIB})
      if(USE_CUDA)
         set_source_files_properties(
            ${CMAKE_SOURCE_DIR}/tests/coarse_graining/benchmark_gpu_adaptive_runtime.cpp
            PROPERTIES LANGUAGE CUDA)
         target_compile_definitions(gpu_adaptive_runtime_benchmark PRIVATE CUDA_V)
         target_link_libraries(gpu_adaptive_runtime_benchmark PRIVATE CUDA::cudart)
      else()
         set_source_files_properties(
            ${CMAKE_SOURCE_DIR}/tests/coarse_graining/benchmark_gpu_adaptive_runtime.cpp
            PROPERTIES LANGUAGE HIP)
         target_compile_definitions(gpu_adaptive_runtime_benchmark PRIVATE HIP_V)
      endif()
      if(UPPASD_PRECISION STREQUAL "SINGLE")
         target_compile_definitions(gpu_adaptive_runtime_benchmark PRIVATE SINGLE_PREC)
      endif()

      # RCG-09 (F-10): the synthetic fused-multiply-add control, separated out
      # of the benchmark above so it can never again be read as a production
      # baseline.  It measures only whether an inactive adaptive runtime costs
      # anything; it contains no UppASD physics.
      add_executable(gpu_adaptive_inactive_overhead_microbenchmark
         ${CMAKE_SOURCE_DIR}/tests/coarse_graining/benchmark_gpu_inactive_runtime_overhead.cpp)
      target_include_directories(gpu_adaptive_inactive_overhead_microbenchmark PRIVATE
         ${CMAKE_SOURCE_DIR}/source/gpu_files
         ${CMAKE_SOURCE_DIR}/source/gpu_files/measurement)
      target_link_libraries(gpu_adaptive_inactive_overhead_microbenchmark PRIVATE ${UppASD_LIB})
      if(USE_CUDA)
         set_source_files_properties(
            ${CMAKE_SOURCE_DIR}/tests/coarse_graining/benchmark_gpu_inactive_runtime_overhead.cpp
            PROPERTIES LANGUAGE CUDA)
         target_compile_definitions(gpu_adaptive_inactive_overhead_microbenchmark PRIVATE CUDA_V)
         target_link_libraries(gpu_adaptive_inactive_overhead_microbenchmark PRIVATE CUDA::cudart)
      else()
         set_source_files_properties(
            ${CMAKE_SOURCE_DIR}/tests/coarse_graining/benchmark_gpu_inactive_runtime_overhead.cpp
            PROPERTIES LANGUAGE HIP)
         target_compile_definitions(gpu_adaptive_inactive_overhead_microbenchmark PRIVATE HIP_V)
      endif()
      if(UPPASD_PRECISION STREQUAL "SINGLE")
         target_compile_definitions(gpu_adaptive_inactive_overhead_microbenchmark PRIVATE SINGLE_PREC)
      endif()

      add_executable(dipole_gpu_fft_tests
         ${CMAKE_SOURCE_DIR}/tests/dipole_validation/test_gpu_dipole_fft.cpp)
      target_include_directories(dipole_gpu_fft_tests PRIVATE
         ${CMAKE_SOURCE_DIR}/source/gpu_files
         ${CMAKE_SOURCE_DIR}/source/gpu_files/measurement)
      target_link_libraries(dipole_gpu_fft_tests PRIVATE ${UppASD_LIB})
      if(USE_CUDA)
         set_source_files_properties(${CMAKE_SOURCE_DIR}/tests/dipole_validation/test_gpu_dipole_fft.cpp
            PROPERTIES LANGUAGE CUDA)
         target_compile_definitions(dipole_gpu_fft_tests PRIVATE CUDA_V)
         target_link_libraries(dipole_gpu_fft_tests PRIVATE CUDA::cudart CUDA::cufft)
      else()
         set_source_files_properties(${CMAKE_SOURCE_DIR}/tests/dipole_validation/test_gpu_dipole_fft.cpp
            PROPERTIES LANGUAGE HIP)
         target_compile_definitions(dipole_gpu_fft_tests PRIVATE HIP_V)
         target_link_libraries(dipole_gpu_fft_tests PRIVATE ${HIPFFT_LIBRARY})
      endif()
      add_test(NAME dipole-gpu-fft-convolution COMMAND dipole_gpu_fft_tests)
      set_tests_properties(dipole-gpu-fft-convolution
         PROPERTIES LABELS "coarse-graining;cg13;${CG13_GPU_LABEL};dipole")
      # Test-only OPEN_FFT E2E seam.  It has no production selector, but
      # exercises both the explicit-tensor diagnostics and the production
      # finite-kernel construction/upload through active-to-padded field and
      # energy evaluation.
      add_executable(dipole_open_fft_layout_tests
         ${CMAKE_SOURCE_DIR}/tests/dipole_validation/open_fft_test_seam.cpp
         ${CMAKE_SOURCE_DIR}/tests/dipole_validation/test_open_fft_layout_red.cpp)
      target_include_directories(dipole_open_fft_layout_tests PRIVATE
         ${CMAKE_SOURCE_DIR}/source/gpu_files
         ${CMAKE_SOURCE_DIR}/source/gpu_files/measurement
         ${CMAKE_SOURCE_DIR}/tests/dipole_validation)
      target_link_libraries(dipole_open_fft_layout_tests PRIVATE ${UppASD_LIB})
      if(USE_CUDA)
         set_source_files_properties(${CMAKE_SOURCE_DIR}/tests/dipole_validation/open_fft_test_seam.cpp
            PROPERTIES LANGUAGE CUDA)
         target_compile_definitions(dipole_open_fft_layout_tests PRIVATE CUDA_V)
         target_link_libraries(dipole_open_fft_layout_tests PRIVATE CUDA::cudart CUDA::cufft)
      else()
         set_source_files_properties(${CMAKE_SOURCE_DIR}/tests/dipole_validation/open_fft_test_seam.cpp
            PROPERTIES LANGUAGE HIP)
         target_compile_definitions(dipole_open_fft_layout_tests PRIVATE HIP_V)
         target_link_libraries(dipole_open_fft_layout_tests PRIVATE ${HIPFFT_LIBRARY})
      endif()
      if(UPPASD_PRECISION STREQUAL "SINGLE")
         target_compile_definitions(dipole_open_fft_layout_tests PRIVATE SINGLE_PREC)
      endif()
      add_test(NAME dipole-open-fft-layout COMMAND dipole_open_fft_layout_tests)
      set_tests_properties(dipole-open-fft-layout
         PROPERTIES LABELS "coarse-graining;cg13;${CG13_GPU_LABEL};dipole")
      # Explicit performance harness.  It is a build target rather than a
      # CTest because timings are hardware-dependent and are recorded by WP9
      # sweep jobs, not used as pass/fail correctness expectations.
      add_executable(dipole_gpu_fft_benchmark
         ${CMAKE_SOURCE_DIR}/tests/dipole_validation/benchmark_gpu_dipole_fft.cpp)
      target_include_directories(dipole_gpu_fft_benchmark PRIVATE
         ${CMAKE_SOURCE_DIR}/source/gpu_files
         ${CMAKE_SOURCE_DIR}/source/gpu_files/measurement)
      target_link_libraries(dipole_gpu_fft_benchmark PRIVATE ${UppASD_LIB})
      if(USE_CUDA)
         set_source_files_properties(${CMAKE_SOURCE_DIR}/tests/dipole_validation/benchmark_gpu_dipole_fft.cpp
            PROPERTIES LANGUAGE CUDA)
         target_compile_definitions(dipole_gpu_fft_benchmark PRIVATE CUDA_V)
         target_link_libraries(dipole_gpu_fft_benchmark PRIVATE CUDA::cudart CUDA::cufft)
      else()
         set_source_files_properties(${CMAKE_SOURCE_DIR}/tests/dipole_validation/benchmark_gpu_dipole_fft.cpp
            PROPERTIES LANGUAGE HIP)
         target_compile_definitions(dipole_gpu_fft_benchmark PRIVATE HIP_V)
         target_link_libraries(dipole_gpu_fft_benchmark PRIVATE ${HIPFFT_LIBRARY})
      endif()
      # real_type.h is consumed directly by these standalone targets, while
      # SINGLE_PREC is private to asdlib.  Mirror the selected storage type so
      # an fp32 validation build cannot accidentally link a fp64 test object.
      if(UPPASD_PRECISION STREQUAL "SINGLE")
         target_compile_definitions(dipole_gpu_fft_tests PRIVATE SINGLE_PREC)
         target_compile_definitions(dipole_gpu_fft_benchmark PRIVATE SINGLE_PREC)
      endif()
      # This invokes the real UppASD GPU SD executable: descriptor/kernel
      # construction, macro-moment refresh, predictor/corrector dispatch, and
      # the GPU measurement writer are checked against the independent oracle.
      set(WP5_E2E_PRECISION_ARGUMENT)
      if(UPPASD_PRECISION STREQUAL "SINGLE")
         set(WP5_E2E_PRECISION_ARGUMENT --fp32)
      endif()
      add_test(NAME dipole-gpu-wp5-e2e
         COMMAND python3 ${CMAKE_SOURCE_DIR}/tests/dipole_validation/run_wp5_e2e.py
                 --binary ${TestBinary} ${WP5_E2E_PRECISION_ARGUMENT})
      set_property(TEST dipole-gpu-wp5-e2e PROPERTY ENVIRONMENT "OMP_NUM_THREADS=1")
      set(OPEN_FFT_E2E_PRECISION_ARGUMENT)
      if(UPPASD_PRECISION STREQUAL "SINGLE")
         set(OPEN_FFT_E2E_PRECISION_ARGUMENT --fp32)
      endif()
      add_test(NAME dipole-open-fft-coarse-e2e
         COMMAND python3 ${CMAKE_SOURCE_DIR}/tests/dipole_validation/run_wp10_open_coarse_e2e.py
                 --binary ${TestBinary} ${OPEN_FFT_E2E_PRECISION_ARGUMENT})
      set_property(TEST dipole-open-fft-coarse-e2e PROPERTY ENVIRONMENT "OMP_NUM_THREADS=1")
   endif()
endif()
message(STATUS "Output binary:  ${TestBinary}")
##############################################################################################
