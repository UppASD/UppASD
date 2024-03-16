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

PROJECT(UppASD Fortran)
# Set the project version
SET(VERSION 5.1.0)
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
SET(MOD ${CMAKE_SOURCE_DIR}/source/mod)
#SET(MOD ${CMAKE_BINARY_DIR}/mod)
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
string(REGEX MATCH gfortran|g95|ifort|pgf90|f95|f77|nv current_compiler  ${CMAKE_Fortran_COMPILER})
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
##############################################################################################


##############################################################################################
set(MKL_LINK static)
find_package(MKL CONFIG)

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
CMAKE_DEPENDENT_OPTION(RUN_CUDA_TESTS "Run cuda-test"  OFF "BUILD_TESTING" OFF)
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

   if(OpenMP_Fortran_FOUND)
      set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${OpenMP_Fortran_FLAGS}")
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
#if(USE_VSL OR USE_MKL OR USE_MKL_FFT)
#
#   set(BLA_DEFINITIONS)
#
#   if (BLA_VENDOR MATCHES "_seq")
#      set(MKL_USE_parallel OFF)
#   else()
#      set(MKL_USE_parallel ON)
#   endif()
#
#   if (MKL_FOUND)
#      message( STATUS "MKL_FOUND condition true")
#      set(BLA_INCLUDE_DIR ${MKL_INCLUDE_DIR})
#      set(BLAS_LIBRARIES ${MKL_LIBRARIES})
#      set(LAPACK_LIBRARIES ${MKL_LIBRARIES})
#      set(BLA_DEFINITIONS USE_MKL)
#      set(MKL_INTERFACE_FULL "lp64")
#   else()
#      #include(FindBLAS)
#      #include(FindLAPACK)
#      find_package(BLAS)
#      find_package(LAPACK REQUIRED)
#   endif()
#   message(STATUS "mkl-location: ${MKL_LIBRARY_DIR}")
#   message(STATUS "mkl-includes: ${MKL_INCLUDE}")
#   link_directories(${MKL_LIBRARY_DIR})
#
#else()
#   find_package(BLAS)
#   find_package(LAPACK REQUIRED)
#endif()
#set(BLA_STATIC ON)
find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)
message(STATUS "Math libs: 
BLAS:${BLAS_FOUND} ${BLAS_LIBRARIES}  ${BLAS_LINKER_FLAGS}
LAPACK:${LAPACK_FOUND} ${LAPACK_LIBRARIES}  ${LAPACK_LINKER_FLAGS}
MKL:${MKL_FOUND} ${MKL_LIBRARIES} ${MKL_LINKER_FLAGS}"
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
INCLUDE(${CMAKE_MODULE_PATH}/git-watcher.cmake)
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
TARGET_LINK_LIBRARIES(${UppASD_LIB}
   -lgfortran
   ${BLAS_LIBRARIES}
   ${LAPACK_LIBRARIES}
   ${CMAKE_THREAD_LIBS_INIT} 
   ${OVF_LIB} 
   ${FFTW_LIB}  
   ${CUDA_curand_LIBRARY}
   ${OpenMP_Fortran_LIBRARIES}
   )

TARGET_LINK_LIBRARIES(${UppASD_EXE}
   -lgfortran
   ${UppASD_LIB}
   ${BLAS_LIBRARIES}
   ${LAPACK_LIBRARIES}
   ${CMAKE_THREAD_LIBS_INIT} 
   ${OVF_LIB} 
   ${FFTW_LIB}  
   ${CUDA_curand_LIBRARY}
   ${OpenMP_Fortran_LIBRARIES}
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

##############################################################################################
# Copy binary to ./source/sd for legacy purposes
##############################################################################################
add_custom_command(
   TARGET ${UppASD_EXE} POST_BUILD
   COMMAND ${CMAKE_COMMAND} -E copy
   ${CMAKE_SOURCE_DIR}/bin/${UppASD_EXE}
   ${LEGACY_OUTPUT_DIRECTORY}/sd
   )


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
