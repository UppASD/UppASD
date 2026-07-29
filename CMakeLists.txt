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

   # CG-04 CPU reference tensor operator: explicit energy/field fixtures on
   # cubic and skew canonical grids, independent of adaptive and GPU paths.
   add_executable(coarse_tensor_operator_tests
      ${CMAKE_SOURCE_DIR}/tests/coarse_graining/test_coarse_tensor_operator.f90)
   set_target_properties(coarse_tensor_operator_tests PROPERTIES
      Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/modules/coarse_tensor_operator_tests)
   target_link_libraries(coarse_tensor_operator_tests PRIVATE ${UppASD_LIB})
   add_test(NAME coarse-graining-tensor-operator COMMAND coarse_tensor_operator_tests)

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
      # CG-09 immutable topology transfer, mutable block/channel storage,
      # compact work scheduling, measured block-mask synchronization, and
      # cleanup accounting on the active CUDA/HIP backend.
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
