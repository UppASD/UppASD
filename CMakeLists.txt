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

PROJECT(uppasd LANGUAGES Fortran)
# Set the project version
SET(VERSION 5.1.0)
##############################################################################################
unset(CMAKE_Fortran_IMPLICIT_LINK_LIBRARIES)
##############################################################################################
# Check OS and compiler status
##############################################################################################
message(STATUS "Operating system= ${CMAKE_SYSTEM_NAME}")
include(CMakeDetermineFortranCompiler)
message(STATUS "CMAKE_Fortran_COMPILER_ID ${CMAKE_Fortran_COMPILER_ID} ${CMAKE_Fortran_COMPILER}.")
message(STATUS "CMAKE_Fortran_IMPLICIT_LINK_LIBRARIES ${CMAKE_Fortran_IMPLICIT_LINK_LIBRARIES}.")
message(STATUS "CMAKE_Fortran_IMPLICIT_LINK_DIRECTORIES ${CMAKE_Fortran_IMPLICIT_LINK_DIRECTORIES}.")
include(CMakeDetermineCCompiler)
message(STATUS "CMAKE_C_COMPILER_ID ${CMAKE_C_COMPILER_ID} ${CMAKE_C_COMPILER}.")
enable_language(C)
##############################################################################################

##############################################################################################
# Set the fortran compiler as the linker
##############################################################################################
#set(CMAKE_LINKER ${CMAKE_Fortran_COMPILER})
#############################################################################################

##############################################################################################
# Add our local modlues to the module path
##############################################################################################
message(STATUS "Initial module path: ${CMAKE_MODULE_PATH}")
SET(LOCAL_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake/")
SET(CMAKE_MODULE_PATH "${CMAKE_MODULE_PATH};${LOCAL_MODULE_PATH}")
message(STATUS "Final module path  : ${CMAKE_MODULE_PATH}")
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
##############################################################################################
SET(SRC ${CMAKE_SOURCE_DIR}/source)
#SET(MOD ${CMAKE_SOURCE_DIR}/source/mod)
SET(MOD ${CMAKE_BINARY_DIR}/mod)
SET(BIN ${CMAKE_SOURCE_DIR}/bin)
# Have the .mod files placed in the lib folder
SET(CMAKE_Fortran_MODULE_DIRECTORY ${MOD})
# Have the .o files placed in the bin folder
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/bin)
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
UNSET(CMAKE_SHARED_LIBRARY_LINK_Fortran_FLAGS)
GET_FILENAME_COMPONENT(FCNAME ${CMAKE_Fortran_COMPILER} NAME)
set(CMAKE_Fortran_COMPILE_OPTIONS_PIE "")
IF(FCNAME STREQUAL "pgf90")
   UNSET(CMAKE_SHARED_LIBRARY_LINK_Fortran_FLAGS)
ENDIF(FCNAME STREQUAL "pgf90")
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
string(REGEX MATCH gfortran|g95|ifort|pgf90|f95|f77|nv|ifx current_compiler  ${CMAKE_Fortran_COMPILER})
message(STATUS "Matched compiler: ${current_compiler}.  ${CMAKE_Fortran_COMPILER} ${FCNAME}")
message(STATUS "Matched compiler: ${current_compiler}.  ${CMAKE_Fortran_COMPILER} ${FCNAME}")
if(current_compiler)
   string(APPEND ASD_binary ".${current_compiler}")
endif()
# add a `cuda` tag to the binary if used.
if(USE_CUDA)
   string(APPEND ASD_binary ".cuda")
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
SET(UppASD_LIB asdlib)
#SET(UppASD_pyLIB pyasdlib)
##############################################################################################


##############################################################################################
# Define the executable in terms of the source files
##############################################################################################
set(lib_type STATIC)
ADD_LIBRARY(${UppASD_LIB} ${lib_type})
message(STATUS " LIB INFO: ${UppASD_LIB} ${lib_type} ${SRC} ${CMAKE_POSITION_INDEPENDENT_CODE}")
if (BUILD_PYTHON)
   set(CMAKE_POSITION_INDEPENDENT_CODE ON)
   set_target_properties(${UppASD_LIB} PROPERTIES POSITION_INDEPENDENT_CODE TRUE)
   #   ADD_LIBRARY(${UppASD_pyLIB} SHARED)
endif(BUILD_PYTHON)
add_dependencies(${UppASD_EXE} ${UppASD_LIB})
##############################################################################################


##############################################################################################
# MKL / FFTW actions
##############################################################################################
#if (USE_MKL)
set(MKL_LINK static)
set(MKL_INTERFACE lp64)
set(MKL_THREADING sequential)
find_package(MKL CONFIG)
message(STATUS "MKL_IMPORTED_TARGETS:   ${MKL_IMPORTED_TARGETS}") 
#endif ()

set(FFTW_USE_STATIC_LIBS ON)
if (USE_FFTW)
find_package(FFTW)
#include(${CMAKE_MODULE_PATH}/FindFFTW.cmake)
if (FFTW_FOUND)
	message(STATUS "FFTW found       : ${FFTW_FOUND}, ${FFTW_INCLUDES}, ${FFTW_INCLUDE_DIRS}")
	message(STATUS "FFTW_LIBRARIES   : ${FFTW_LIBRARIES}")
	message(STATUS "FFTW_INCLUDE_DIRS: ${FFTW_INCLUDE_DIRS}")
	message(STATUS "FFTW_DOUBLE_LIB  : ${FFTW_DOUBLE_LIB}")
	message(STATUS "FFTW_DOUBLE_OPENMP_LIB  : ${FFTW_DOUBLE_OPENMP_LIB}")
   set(FFTW_LIB ${FFTW_LIBRARIES})
   #set(FFTW_LIB "${FFTW_DOUBLE_LIB};${FFTW_DOUBLE_OPENMP_LIB}")
   # It sets the following variables:
#   FFTW_FOUND                  ... true if fftw is found on the system
#   FFTW_[component]_LIB_FOUND  ... true if the component is found on the system (see components below)
#   FFTW_LIBRARIES              ... full paths to all found fftw libraries
#   FFTW_[component]_LIB        ... full path to one of the components (see below)
#   FFTW_INCLUDE_DIRS           ... fftw include directory paths

else (FFTW_FOUND)
   set(FFTW_LIB "")
   set(FFTW_LIBRARIES "")
endif (FFTW_FOUND)
endif (USE_FFTW)
##############################################################################################

##############################################################################################
# Set relevant user-controlled options for the compilation
##############################################################################################
include(CMakeDependentOption)
OPTION(USE_OPENMP "Use OpenMP for parallelization" ON)
OPTION(USE_CUDA "Compile with CUDA support" OFF)
# Uncomment below to make gnu defaults without MKL
if(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
   CMAKE_DEPENDENT_OPTION(USE_MKL "Use Intel MKL for BLAS/LAPACK" OFF "MKL_FOUND" OFF)
   CMAKE_DEPENDENT_OPTION(USE_VSL "Use Intel VSL (needs MKL)" OFF "MKL_FOUND" OFF)
   CMAKE_DEPENDENT_OPTION(USE_MKL_FFT "Use Intel MKL for FFT" OFF "MKL_FOUND"  OFF)
   CMAKE_DEPENDENT_OPTION(USE_FFTW "Use FFTW for FFT" OFF "FFTW_FOUND" OFF)
else()
   CMAKE_DEPENDENT_OPTION(USE_MKL "Use Intel MKL for BLAS/LAPACK" ON "MKL_FOUND" OFF)
   CMAKE_DEPENDENT_OPTION(USE_VSL "Use Intel VSL (needs MKL)" ON "MKL_FOUND" OFF)
   CMAKE_DEPENDENT_OPTION(USE_MKL_FFT "Use Intel MKL for FFT" ON "MKL_FOUND"  OFF)
   CMAKE_DEPENDENT_OPTION(USE_FFTW "Use FFTW for FFT" OFF "FFTW_FOUND" OFF)
endif()
CMAKE_DEPENDENT_OPTION(BUILD_CUDA_DOCS "Build LaTeX notes of the CUDA implementation" OFF "BUILD_MANUAL" OFF) 
CMAKE_DEPENDENT_OPTION(BUILD_GNEB_DOCS "Build LaTeX notes of the GNEB implementation" OFF "BUILD_MANUAL" OFF) 
OPTION(USE_OVF "Use OVF format" OFF)
OPTION(BUILD_MANUAL "Build manual"  OFF)
OPTION(BUILD_DOXYGEN "Build doxygen"  OFF)
OPTION(BUILD_PYTHON "Build Python interface"  OFF)

##############################################################################################
# Enable and define tests
##############################################################################################
include(CTest)
CMAKE_DEPENDENT_OPTION(RUN_ASD_TESTS "Run asd-tests"  ON "BUILD_TESTING" OFF)
CMAKE_DEPENDENT_OPTION(RUN_SLD_TESTS "Run sld-tests"  OFF "BUILD_TESTING" OFF)
CMAKE_DEPENDENT_OPTION(RUN_GNEB_TESTS "Run gneb-tests"  OFF "BUILD_TESTING" OFF)
CMAKE_DEPENDENT_OPTION(RUN_REG_TESTS "Run regression-test"  ON "BUILD_TESTING" OFF)
CMAKE_DEPENDENT_OPTION(RUN_CUDA_TESTS "Run cuda-test"  ON "BUILD_TESTING;USE_CUDA" OFF)
##############################################################################################

##############################################################################################
# Determine and set compiler-dependent flags
##############################################################################################
INCLUDE(${LOCAL_MODULE_PATH}/SetFortranFlags.cmake) 
##############################################################################################

##############################################################################################
# Add OpenMP compiler/linker flags
##############################################################################################
IF(USE_OPENMP)
   FIND_PACKAGE(OpenMP)
   message(STATUS "Did we find OpenMP? ${OpenMP_Fortran_FOUND}")
   message(STATUS "            OpenMP flags: ${OpenMP_Fortran_FLAGS}")
   message(STATUS "            OpenMP libraries: ${OpenMP_Fortran_LIBRARIES}")

   if(OpenMP_Fortran_FOUND)
      set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${OpenMP_Fortran_FLAGS}")
   endif()

ENDIF(USE_OPENMP)
##############################################################################################

#set(CMAKE_Fortran_FLAGS "-assume underscore")

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
set(BLA_STATIC ON)

message( STATUS "USING::::: ${USE_VSL} ${USE_MKL} ${USE_MKL_FFT}")
if(USE_VSL OR USE_MKL OR USE_MKL_FFT)

   unset(BLA_DEFINITIONS)

   if (BLA_VENDOR MATCHES "_seq")
      set(MKL_USE_parallel OFF)
   else()
      set(MKL_USE_parallel ON)
   endif()

   if (MKL_FOUND)
      message( STATUS "MKL_FOUND condition true")
      #set(BLA_INCLUDE_DIR ${MKL_INCLUDE_DIR})
      #set(BLAS_LIBRARIES ${MKL_LIBRARIES})
      #set(LAPACK_LIBRARIES ${MKL_LIBRARIES})
      #set(BLA_DEFINITIONS USE_MKL)
      #set(MKL_INTERFACE_FULL "lp64")
      #set(BLA_VENDOR Intel)
      #find_package(BLAS REQUIRED)
      #find_package(LAPACK REQUIRED)
   else()
	   #include(FindBLAS)
      	   #include(FindLAPACK)
      find_package(BLAS REQUIRED)
      find_package(LAPACK REQUIRED)
   endif()
   message(STATUS "mkl-location: ${MKL_LIBRARY_DIR}")
   message(STATUS "mkl-includes: ${MKL_INCLUDE}")
   #link_directories(${MKL_LIBRARY_DIR})

   target_compile_options(${UppASD_LIB} PUBLIC $<TARGET_PROPERTY:MKL::MKL,INTERFACE_COMPILE_OPTIONS>)
   target_include_directories(${UppASD_LIB} PUBLIC $<TARGET_PROPERTY:MKL::MKL,INTERFACE_INCLUDE_DIRECTORIES>)

else()
   find_package(BLAS REQUIRED)
   find_package(LAPACK REQUIRED)
endif()
#find_package(BLAS REQUIRED)
#find_package(LAPACK REQUIRED)
message(STATUS "Math libs: 
 BLAS:    ${BLAS_FOUND}, ${BLAS_LIBRARIES},  ${BLAS_LINKER_FLAGS}
 LAPACK:  ${LAPACK_FOUND}, ${LAPACK_LIBRARIES},  ${LAPACK_LINKER_FLAGS}
 MKL:     ${MKL_FOUND}, ${MKL_LIBRARIES}, ${MKL_LINKER_FLAGS}"
)
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
INCLUDE(${LOCAL_MODULE_PATH}/git-watcher.cmake)
##############################################################################################

##############################################################################################
# Add CUDA flags and preprocessor commands
##############################################################################################
if (USE_CUDA) 

   ENABLE_LANGUAGE(C CUDA CXX)
   if (CMAKE_VERSION VERSION_LESS "3.17.0")
	   find_package(CUDA)
   else ()
   	find_package(CUDAToolkit)
   endif ()
   target_sources(${UppASD_LIB} PRIVATE ${SRC}/chelper.f90)
   add_subdirectory(${SRC}/gpu_files)
   add_compile_definitions(CUDA)
   set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} -I${CMAKE_CUDA_TOOLKIT_INCLUDE_DIRECTORIES} -DDEBUG")
   set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -I${CMAKE_CUDA_TOOLKIT_INCLUDE_DIRECTORIES}")
   if (DEFINED CUDAToolkit_MATH_INCLUDE_DIR )
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -I${CUDAToolkit_MATH_INCLUDE_DIR}")
   endif()

   # Set preprocessing flags indicating compiler
   #string(TOUPPER ${current_compiler} PP_CUDA_COMPILER_FLAG)
   set(PP_CUDA_COMPILER_FLAG "${CMAKE_Fortran_COMPILER_ID}")
   set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -D__${PP_CUDA_COMPILER_FLAG}__")
   set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} -D__${PP_CUDA_COMPILER_FLAG}__")
   message(STATUS "CUDA PPFLAG: __${PP_CUDA_COMPILER_FLAG}__, ${CMAKE_Fortran_COMPILER_ID}")

   # Find CUDA arch on local machine
   include(FindCUDA/select_compute_arch)
   CUDA_DETECT_INSTALLED_GPUS(INSTALLED_GPU_CCS_1)
   string(STRIP "${INSTALLED_GPU_CCS_1}" INSTALLED_GPU_CCS_2)
   string(REPLACE " " ";" INSTALLED_GPU_CCS_3 "${INSTALLED_GPU_CCS_2}")
   string(REPLACE "." "" CUDA_ARCH_LIST "${INSTALLED_GPU_CCS_3}")
   SET(CMAKE_CUDA_ARCHITECTURES ${CUDA_ARCH_LIST})
   SET_TARGET_PROPERTIES(${UppASD_LIB} PROPERTIES CUDA_ARCHITECTURES "${CUDA_ARCH_LIST}")

else ()
   target_sources(${UppASD_LIB} PRIVATE ${SRC}/Tools/nocuda.f90)
endif (USE_CUDA)
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
# MKL edits in case of GNU compiler
##############################################################################################
if(USE_MKL)
   if(CMAKE_Fortran_COMPILER_ID MATCHES "GNU" OR CMAKE_Fortran_COMPILER_ID MATCHES "NVHPC")
      target_include_directories(${UppASD_LIB} PRIVATE ${MKL_INCLUDE})
   endif()
   set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -D__INTEL_MKL__")

endif()
##############################################################################################


##############################################################################################
# Include VSL preprocessor flag (recommended)
##############################################################################################
if(USE_VSL AND MKL_FOUND)
   add_compile_definitions(VSL)
   #add_compile_definitions(USE_VSL)
endif()
##############################################################################################


##############################################################################################
# Check for FFT treatment (MKL default, FFTW available)
##############################################################################################
if(USE_MKL_FFT AND MKL_FOUND)
   target_sources(${UppASD_LIB} PRIVATE ${SRC}/Hamiltonian/fftdipole_mkl.f90)
   add_compile_definitions(USE_MKL_FFT)
   include_directories(${BLA_INCLUDE_DIR})
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
   #set(Python_VERSION $ENV{PYTHON_VERSION})
   #message(STATUS "Python target version: ${Python_VERSION}")
   message(STATUS "Python target version: $ENV{PYTHON_VERSION}")

   # Fix since FindPython is a bit flawed
   # https://github.com/scikit-build/scikit-build/issues/506#issuecomment-706791791
   find_package(Python3 COMPONENTS Interpreter Development.Module NumPy)
   
   message(STATUS "============================================")
   message(STATUS "Python3_EXECUTABLE         ${Python3_EXECUTABLE}")
   message(STATUS "Python3_SOABI              ${Python3_SOABI}")
   message(STATUS "Python3_STDLIB             ${Python3_STDLIB}.")
   message(STATUS "Python3_SITELIB            ${Python3_SITELIB}")
   message(STATUS "Python3_LIBRARIES          ${Python3_LIBRARIES}.")
   message(STATUS "Python3_LIBRARY_DIRS       ${Python3_LIBRARY_DIRS}.")
   message(STATUS "Python3_INCLUDE_DIRS       ${Python3_INCLUDE_DIRS}.")
   message(STATUS "Python3_NumPy_INCLUDE_DIRS ${Python3_NumPy_INCLUDE_DIRS}.")
   message(STATUS "============================================")

   if (Python3_FOUND)
      # Grab the variables from a local Python installation
      # F2PY headers
      execute_process(
         COMMAND "${Python3_EXECUTABLE}" -c "import numpy.f2py; print(numpy.f2py.get_include())"
         OUTPUT_VARIABLE F2PY_INCLUDE_DIR
         OUTPUT_STRIP_TRAILING_WHITESPACE
         )

      execute_process(
         COMMAND "${Python3_EXECUTABLE}" -c "import numpy; print(numpy.version.version)"
         OUTPUT_VARIABLE F2PY_version
         OUTPUT_STRIP_TRAILING_WHITESPACE
         )

      # Common variables
      set(f2py_module_name "uppasd")
      set(f2py_module_c "_${f2py_module_name}module.c")
      set(f2py_wrapper_f "_${f2py_module_name}-f2pywrappers.f")
      set(f2py_wrapper_f90 "_${f2py_module_name}-f2pywrappers2.f90")
      set(f2py_wrapper_o "_${f2py_module_name}-f2pywrappers2.o")
      set(UppASD_F90w f90_wraps)
      set(F90wrap_kind
         ${CMAKE_SOURCE_DIR}/cmake/.f2py_f2cmap
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
         ${SRC}/pyasd.f90
      #   ${SRC}/f90wrap/f90wrap_pyasd.f90
      #   #${SRC}/f90wrap/f90wrap_uppasd.f90
      #   #${SRC}/f90wrap/f90wrap_inputdatatype.f90
      #   #${SRC}/f90wrap/f90wrap_inputdata.f90
      #   #${SRC}/f90wrap/f90wrap_inputhandler.f90
      #   #${SRC}/f90wrap/f90wrap_inputhandler_ext.f90
      #   #${SRC}/f90wrap/f90wrap_simulationdata.f90
      #   #${SRC}/f90wrap/f90wrap_momentdata.f90
         )


##############################################################################################
# Transform F2PY include-string to acceptable slash convention 
##############################################################################################
      # Input file name string (with double backslashes)
      if (F2PY_INCLUDE_DIR)
         #message(STATUS "Original F2PY-include: ${F2PY_INCLUDE_DIR}")
         file(TO_CMAKE_PATH ${F2PY_INCLUDE_DIR} F2PY_INCLUDE_DIR)
         #message(STATUS "Transformed F2PY-include: ${F2PY_INCLUDE_DIR}")
      endif ()
##############################################################################################

##############################################################################################
# New Python generation
##############################################################################################
      add_custom_target ( ${UppASD_F90w}
         DEPENDS "${CMAKE_CURRENT_BINARY_DIR}/${f2py_module_c}"
         VERBATIM
         )

      #add_custom_command( 
      #   OUTPUT "${CMAKE_CURRENT_BINARY_DIR}/${f2py_wrapper_o}"
      #   COMMAND ${CMAKE_Fortran_COMPILER} -c  ${CMAKE_CURRENT_BINARY_DIR}/${f2py_wrapper_f}  -I${MOD}
      #   COMMENT "Writing F90wrappers for F2PY"
      #   DEPENDS "${CMAKE_CURRENT_BINARY_DIR}/${f2py_wrapper_f}"
      #   )

      add_custom_command(
         OUTPUT "${CMAKE_CURRENT_BINARY_DIR}/${f2py_module_c}" "${CMAKE_CURRENT_BINARY_DIR}/${f2py_wrapper_f}"
         "${CMAKE_CURRENT_BINARY_DIR}/${f2py_wrapper_f90}"  
         COMMAND ${Python3_EXECUTABLE}  -m "numpy.f2py"  ${F2PY_sources} -m _${f2py_module_name}  --lower --quiet
         #COMMAND f2py   ${F2PY_sources} -m _${f2py_module_name}  --lower --quiet
         COMMENT "Generating Python modules"
         DEPENDS ${UppASD_LIB}  ${F2PY_sources}
         VERBATIM
         )

      message(STATUS "F2PY_INCLUDE_DIR" ${F2PY_INCLUDE_DIR})
      #message(STATUS "Python3_NumPy_INCLUDE_DIRS" ${Python3_NumPy_INCLUDE_DIRS})

      # Set up target
      #message(STATUS "F2PY_INCLUDE_DIR:        ${F2PY_INCLUDE_DIR}")
      #message(STATUS "=============>  ${CMAKE_CURRENT_BINARY_DIR}/${f2py_module_c}")

      # Not used by default (invoke only if new header files need to be generated)
      # add_custom_command( 
      #    OUTPUT "${F2PY_sources}"
      #    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
      #    COMMAND f90wrap  -P -k ${F90wrap_kind} -m ${f2py_module_name}  ${F90wrap_sources}
      #    COMMENT "Writing F90wrappers for F2PY"
      #    DEPENDS ${F90wrap_sources}
      #    )

      ##message(STATUS " F2PY_sauce: ${F2PY_sources}")
      set(CMAKE_C_FLAGS "-std=c99")
      Python3_add_library(_${CMAKE_PROJECT_NAME} MODULE WITH_SOABI
         "${CMAKE_CURRENT_BINARY_DIR}/${f2py_module_c}"
         #"${F2PY_sources}"
         "${F2PY_INCLUDE_DIR}/fortranobject.c" # From NumPy
         "${CMAKE_CURRENT_BINARY_DIR}/${f2py_wrapper_f}"
         #"${CMAKE_CURRENT_BINARY_DIR}/${f2py_wrapper_f90}"
         )

      add_dependencies(_${CMAKE_PROJECT_NAME} ${UppASD_F90w})

      target_link_libraries(_${CMAKE_PROJECT_NAME}  PRIVATE Python3::NumPy  ${UppASD_LIB})

      #if(USE_VSL OR USE_MKL OR USE_MKL_FFT)
      #   TARGET_LINK_LIBRARIES(_${CMAKE_PROJECT_NAME} PUBLIC  $<LINK_ONLY:MKL::MKL>)
      #endif(USE_VSL OR USE_MKL OR USE_MKL_FFT)
      message(STATUS " Python linking check: 
            ${UppASD_LIB},  ${BLAS_LIBRARIES}, ${OpenMP_Fortran_LIBRARIES}")
      message(STATUS " CMAKE_MODULE_LINKER_FLAGS : ${CMAKE_MODULE_LINKER_FLAGS}")
      target_include_directories(_${CMAKE_PROJECT_NAME} PRIVATE "${F2PY_INCLUDE_DIR}" )
##############################################################################################

##############################################################################################
# Set static linking for different compilers (for Python compilation)
##############################################################################################
if (${CMAKE_SYSTEM_NAME} MATCHES "Windows")
   if( CMAKE_Fortran_COMPILER_ID STREQUAL Intel )
      #set_target_properties( _${CMAKE_PROJECT_NAME} PROPERTIES LINK_FLAGS -static -static-intel )
      #set(CMAKE_MODULE_LINKER_FLAGS "-static-intel -lifcoremt")
      set(CMAKE_MODULE_LINKER_FLAGS "-static -static-intel  -liomp5 /libs:static")
   elseif( CMAKE_Fortran_COMPILER_ID STREQUAL GNU )
      #set_target_properties( _${CMAKE_PROJECT_NAME} PROPERTIES LINK_FLAGS -static -lgfortran -libquadmath)
      #set(CMAKE_MODULE_LINKER_FLAGS "-lgfortran -lquadmath")
   endif()
endif()
##############################################################################################

#############################################################################################
# Copy the library to the designated place
#############################################################################################
###add_custom_command(
###   TARGET ${UppASD_F90w} POST_BUILD
###   COMMAND ${CMAKE_COMMAND} -E 
###   #${CMAKE_BINARY_DIR}/_*.so
###   copy -t ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}
###   ${CMAKE_BINARY_DIR}/_*.*
###   )
#message(STATUS "COPY:::" ${CMAKE_BINARY_DIR} )
#############################################################################################

   else (Python3_FOUND)
      message("A proper Python environment needs to be installed to generate the python libraries")
   endif (Python3_FOUND)
endif(BUILD_PYTHON)
##############################################################################################


#if (${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
#   set(OpenMP_Fortran_LIBRARIES "${OpenMP_Fortran_FLAGS}")
#endif()



##############################################################################################
mark_as_advanced( ENABLE_BLACS ENABLE_BLAS95 ENABLE_CDFT ENABLE_CPARDISO ENABLE_LAPACK95 ENABLE_OMP_OFFLOAD ENABLE_SCALAPACK MKL_DIR MKL_H MKL_LINK MKL_MPI mkl_gf_ilp64_file mkl_intel_thread_file mkl_core_file mkl_intel_ilp64_file mkl_gnu_thread_file mkl_intel_lp64_file mkl_gf_lp64_file OMP_LIBRARY)

##############################################################################################
# Add the needed libraries and special compiler flags
##############################################################################################
# First fix for msys/mingw having semi-shared .dll.a libraries
##############################################################################################
# message(STATUS "***********************************************")
# message(STATUS "--> MSYS ENVIROnMENT BEFORE PURGE:")
# message(STATUS "     BLAS_LIBRARIES: ${BLAS_LIBRARIES}")
# message(STATUS "     LAPAC_LIBRARIES: ${LAPACK_LIBRARIES}")
# message(STATUS "     OpenMP_Fortran_FLAGS: ${OpenMP_Fortran_FLAGS}")
# message(STATUS "     OpenMP_Fortran_LIB_NAMES: ${OpenMP_Fortran_LIB_NAMES}")
# message(STATUS "     OpenMP_Fortran_LIBRARIES: ${OpenMP_Fortran_LIBRARIES}")
# message(STATUS "     CMAKE_THREAD_LIBS_INIT: ${CMAKE_THREAD_LIBS_INIT}")
# message(STATUS "     FFTW_LIBRARIES: ${FFTW_LIBRARIES}")
# message(STATUS "     FFTW_LIB: ${FFTW_LIB}")

# message(STATUS " MSYS ENVIROnMENT $ENV{MSYSTEM} ${BLAS_LIBRARIES}")
# message(STATUS "-> MSYS ENVIROnMENT $ENV{MSYSTEM} ${BLAS_LIBRARIES}, ${OpenMP_Fortran_LIB_NAMES}, ${OpenMP_Fortran_LIBRARIES}.")
if (${CMAKE_SYSTEM_NAME} MATCHES "Windows")
#if ( $ENV{MSYSTEM})
   message(STATUS " Performing MSYS library pruning.")
   if(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
   #if ( ($ENV{MSYSTEM} STREQUAL MINGW) OR  ($ENV{MSYSTEM} STREQUAL UCRT) )
      message(STATUS " Performing MSYS library pruning for $ENV{MSYSTEM}.")
      # BLAS
      unset(UppASD_links)
      foreach ( elem ${BLAS_LIBRARIES})
         string(REPLACE ".dll.a" ".a" filtered_elem ${elem})
         set(UppASD_links "${filtered_elem}")
      endforeach()
      unset(BLAS_LIBRARIES)
      set(BLAS_LIBRARIES ${UppASD_links})
      message(STATUS "BLAS_LIBRARIES pruned: ${BLAS_LIBRARIES} ")
      # LAPACK
      unset(UppASD_links)
      foreach ( elem ${LAPACK_LIBRARIES})
         string(REPLACE ".dll.a" ".a" filtered_elem ${elem})
	      list(APPEND UppASD_links "${filtered_elem}")
	      #set(UppASD_links "${filtered_elem}")
	      #message(STATUS " PRUNING LAPACK: ", ${filtered_elem})
	      #message(STATUS " PRUNING LAPACK::: ", ${UppASD_links})
      endforeach()
      unset(LAPACK_LIBRARIES)
      #message(STATUS " LAPACK PRUNED: ", ${filtered_elem})
      set(LAPACK_LIBRARIES ${UppASD_links})
      # FFTW
      unset(UppASD_links)
      foreach ( elem ${FFTW_LIBRARIES})
         string(REPLACE ".dll.a" ".a" filtered_elem ${elem})
         set(UppASD_links "${filtered_elem}")
      endforeach()
      unset(FFTW_LIBRARIES)
      set(FFTW_LIBRARIES ${UppASD_links})
      # OpenMP
      unset(UppASD_links)
         foreach ( elem ${OpenMP_Fortran_LIBRARIES})
         string(REPLACE ".dll.a" ".a" filtered_elem ${elem})
      set(UppASD_links "${filtered_elem}")
      endforeach()
      unset(OpenMP_Fortran_LIBRARIES)
      set(OpenMP_Fortran_LIBRARIES ${UppASD_links})
      unset(UppASD_links)
   #endif(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
   endif ()
endif ()

# message(STATUS "--> MSYS ENVIROnMENT $ENV{MSYSTEM} ${BLAS_LIBRARIES}, ${OpenMP_Fortran_LIB_NAMES}, ${OpenMP_Fortran_LIBRARIES}.")
# message(STATUS "***********************************************")
# message(STATUS "--> MSYS ENVIROnMENT AFTER PURGE:")
# message(STATUS "     BLAS_LIBRARIES: ${BLAS_LIBRARIES}")
# message(STATUS "     LAPAC_LIBRARIES: ${LAPACK_LIBRARIES}")
# message(STATUS "     OpenMP_Fortran_FLAGS: ${OpenMP_Fortran_FLAGS}")
# message(STATUS "     OpenMP_Fortran_LIB_NAMES: ${OpenMP_Fortran_LIB_NAMES}")
# message(STATUS "     OpenMP_Fortran_LIBRARIES: ${OpenMP_Fortran_LIBRARIES}")
# message(STATUS "     CMAKE_THREAD_LIBS_INIT: ${CMAKE_THREAD_LIBS_INIT}")
# message(STATUS "     FFTW_LIBRARIES: ${FFTW_LIBRARIES}")
# message(STATUS "     FFTW_LIB: ${FFTW_LIB}")
# message(STATUS "***********************************************")

# Detect the platform and configure linker flags and libraries
unset(UppASD_links)
#if (CMAKE_SYSTEM_NAME STREQUAL "Linux" OR CMAKE_SYSTEM_NAME STREQUAL "Windows")
if (CMAKE_SYSTEM_NAME STREQUAL "Windows")
    # On Linux, use static linking and include libquadmath
    set(UppASD_links "-static")
#elseif (CMAKE_SYSTEM_NAME STREQUAL "Darwin")
#    # On macOS, avoid static linking
#    set(CMAKE_STATIC_LINKER_FLAGS "${CMAKE_STATIC_LINKER_FLAGS} -dynamic")
endif()

TARGET_LINK_LIBRARIES(${UppASD_LIB} PUBLIC
   # -static
   # -fopenmp
   # -lgomp
   # -lgfortran
   # -lquadmath
         ${UppASD_links}
         ${BLAS_LIBRARIES}
         ${LAPACK_LIBRARIES}
         ${CMAKE_THREAD_LIBS_INIT} 
         ${OVF_LIB} 
         ${FFTW_LIB}  
         ${CUDA_curand_LIBRARY}
	 # ${OpenMP_Fortran_LIB_NAMES}
         ${OpenMP_Fortran_FLAGS}
         ${OpenMP_Fortran_LIBRARIES} 
      )
      if(MKL_FOUND AND USE_MKL)
         TARGET_LINK_LIBRARIES(${UppASD_LIB} PUBLIC $<LINK_ONLY:MKL::MKL>)
      endif(MKL_FOUND AND USE_MKL)
      if(USE_VSL OR USE_MKL OR USE_MKL_FFT)
         TARGET_LINK_LIBRARIES(${UppASD_LIB} PUBLIC $<LINK_ONLY:MKL::MKL>)
      endif(USE_VSL OR USE_MKL OR USE_MKL_FFT)
   

TARGET_LINK_LIBRARIES(${UppASD_EXE} PRIVATE ${UppASD_LIB})


##############################################################################################
# Tell where to install this executable
##############################################################################################
IF(WIN32)
   SET(CMAKE_INSTALL_PREFIX "C:\\Program Files")
ELSE()
   SET(CMAKE_INSTALL_PREFIX /usr/local)
ENDIF(WIN32)


SET_TARGET_PROPERTIES(${UppASD_EXE} PROPERTIES LINKER_LANGUAGE "Fortran")
SET_TARGET_PROPERTIES(${UppASD_LIB} PROPERTIES LINKER_LANGUAGE "Fortran")
if (NOT BUILD_PYTHON)
   #INSTALL(TARGETS ${UppASD_LIB} RUNTIME DESTINATION lib)
   INSTALL(TARGETS ${UppASD_EXE} RUNTIME DESTINATION bin)
   #INSTALL(TARGETS UppASD_F90w RUNTIME DESTINATION bin)
else()
   #   INSTALL(TARGETS ${UppASD_pyLIB} RUNTIME DESTINATION lib)
   INSTALL(TARGETS ${UppASD_LIB} RUNTIME DESTINATION lib)
   INSTALL(TARGETS ${UppASD_EXE} RUNTIME DESTINATION bin)
endif(NOT BUILD_PYTHON)
install(FILES "${CMAKE_SOURCE_DIR}/LICENSE" DESTINATION .)
install(FILES "${CMAKE_SOURCE_DIR}/README.md" DESTINATION .)

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
set(TestBinary ${CMAKE_SOURCE_DIR}/bin/${UppASD_EXE})
set(TestDIR ${CMAKE_SOURCE_DIR}/tests)
if (RUN_REG_TESTS)
add_test(NAME regression-test  COMMAND python3 -u bergtest.py --file regressionResaro.yaml --binary ${TestBinary}  WORKING_DIRECTORY ${TestDIR})
endif (RUN_REG_TESTS)
if (RUN_ASD_TESTS)
   add_test(NAME asd-tests  COMMAND python3 -u bergtest.py --file regulartests.yaml --binary ${TestBinary}  WORKING_DIRECTORY ${TestDIR})
endif (RUN_ASD_TESTS)
if (RUN_SLD_TESTS)
   add_test(NAME sld-tests  COMMAND python3 -u ./bergtest.py --file sldtests.yaml --binary ${TestBinary}  WORKING_DIRECTORY ${TestDIR})
endif (RUN_SLD_TESTS)
if (RUN_GNEB_TESTS)
   add_test(NAME gneb-tests  COMMAND python3 -u ./bergtest.py --file regressionGNEB_15d7169b.yaml --binary ${TestBinary}  WORKING_DIRECTORY ${TestDIR})
endif (RUN_GNEB_TESTS)
if (RUN_CUDA_TESTS)
add_test(NAME cuda-tests  COMMAND python3 -u bergtest.py --file cudatests.yaml --binary ${TestBinary}  WORKING_DIRECTORY ${TestDIR})
endif (RUN_CUDA_TESTS)
set_property(TEST asd-tests PROPERTY ENVIRONMENT "OMP_NUM_THREADS=1")
message(STATUS "Output binary:  ${TestBinary}")
##############################################################################################
message(STATUS "--> ${CMAKE_SYSTEM_NAME}")


##############################################################################################
## ADD Cpack
##############################################################################################
# Include CPack

# Set package properties
set(CPACK_PACKAGE_NAME "UppASD")
set(CPACK_PACKAGE_VERSION "1.0.1")
set(CPACK_PACKAGE_CONTACT "uppasd@physics.uu.se")
set(CPACK_PACKAGE_VENDOR "The UppASD group")
set(CPACK_PACKAGE_ICON "${CMAKE_SOURCE_DIR}/cpack/icon.jpg")
set(CPACK_RESOURCE_FILE_LICENSE "${CMAKE_SOURCE_DIR}/LICENSE")
set(CPACK_RESOURCE_FILE_README "${CMAKE_SOURCE_DIR}/README.md")
set(CPACK_RESOURCE_FILE_WELCOME "${CMAKE_SOURCE_DIR}/cpack/WELCOME.txt")
set(CPACK_PACKAGE_DESCRIPTION_FILE "${CMAKE_SOURCE_DIR}/cpack/README.txt")
set(CPACK_PACKAGE_DESCRIPTION_SUMMARY "UppASD is a package for atomistic spin dynamics and Monte Carlo simulations.")
set(CPACK_PACKAGE_HOMEPAGE_URL "https://github.com/UppASD/UppASD")
#set(CPACK_GENERATOR "STGZ;DEB;RPM")

if(WIN32)
   include(InstallRequiredSystemLibraries)

    # Include NSIS for creating Windows installers
    # set(CPACK_GENERATOR "NSIS64")

    # Include the custom NSIS script
    set(CPACK_NSIS_INCLUDE "${CMAKE_SOURCE_DIR}/cpack/CustomInstallNSIS.nsh")

    # Additional NSIS settings if needed
    set(CPACK_NSIS_MODIFY_PATH ON)
    set(CPACK_NSIS_ENABLE_UNINSTALL_BEFORE_INSTALL ON)
    set(CPACK_NSIS_EXTRA_INSTALL_COMMANDS "Call AddToPath")

    # WIX installer
    # set(CPACK_WIX_UPGRADE_GUID "E25F58C1-EC45-422A-9700-7254DB19547B")
    # WIX installer
    set(CPACK_WIX_LICENSE_RTF "${CMAKE_SOURCE_DIR}/cpack/LICENSE.rtf")
    set(CPACK_WIX_UI_BANNER "${CMAKE_SOURCE_DIR}/cpack/banner.jpg")
    set(CPACK_WIX_UI_DIALOG "${CMAKE_SOURCE_DIR}/cpack/dialog.jpg")
    set(CPACK_WIX_UPGRADE_GUID "02cf4cf7-fa5b-43de-a8ff-b39098d35db3")
    set(CPACK_WIX_PATCH_FILE "${CMAKE_SOURCE_DIR}/cpack/wixpath.xml")
    # set(CPACK_WIX_EXTRA_SOURCES "${CMAKE_SOURCE_DIR}/PathEnvironment.wxs")
   set(CPACK_WIX_VERSION "4")
   endif()

include(CPack)
