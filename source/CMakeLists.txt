##############################################################################################
# CMAKE Configurations for UppASD
# (C) Anders Bergman 2019
#
# Acknowledgements:
# https://github.com/SethMMorton/cmake_fortran_template  (for starting template)
# https://github.com/openmeeg/findmkl_cmake/blob/master/cmake/FindMKL.cmake (for FindMKL)
# https://github.com/libigl/eigen/blob/master/cmake/FindFFTW.cmake (for FindFFTW)
# https://github.com/andrew-hardin/cmake-git-version-tracking (for git tracking)
#
##############################################################################################


##############################################################################################
# Add CUDA flags and preprocessor commands
##############################################################################################
if (USE_CUDA) 
   ENABLE_LANGUAGE(C CUDA CXX)
   find_package(CUDA)
   add_sources(chelper.f90)
   add_subdirectory(gpu_files)
   add_compile_definitions(CUDA)
   set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} -I${CMAKE_CUDA_TOOLKIT_INCLUDE_DIRECTORIES}")
   #set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -I${CMAKE_CUDA_TOOLKIT_INCLUDE_DIRECTORIES}")
   set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -I${CMAKE_CUDA_TOOLKIT_INCLUDE_DIRECTORIES}")
   #message(STATUS "CUDA_LIBRARIES ${CUDA_LIBRARIES}.")
   #message(STATUS "CUDA_NVCC_FLAGS ${CUDA_NVCC_FLAGS}")
   #message(STATUS "CMAKE_CUDA_TOOLKIT_INCLUDE_DIRECTORIES ${CMAKE_CUDA_TOOLKIT_INCLUDE_DIRECTORIES}")
   #message(STATUS "CUDA_CURAND ${CUDA_CURAND}.")
   #message(STATUS "CUDA_curand_LIBRARY ${CUDA_curand_LIBRARY}.")
else ()
   add_sources(Tools/nocuda.f90)
endif (USE_CUDA)
##############################################################################################


##############################################################################################
# Link with OVF library if available 
##############################################################################################
IF(USE_OVF)
   find_library(OVF_LIB ovf /lib/python3.7/site-packages/ovf/ /usr/lib64/)
   message(STATUS "Looking for OVF library, found in ${OVF_LIB}.")
   add_compile_definitions(USE_OVF)
   add_sources(Third_party/ovf/fortran/ovf.f90)
   ### ENABLE_LANGUAGE(Fortran C CXX)
   ### set(CMAKE_CXX_STANDARD 11)

   ### IF(APPLE)
   ###    SET(CMAKE_CXX_COMPILER gcc-8)
   ### ENDIF(APPLE)

   ### add_compile_definitions(USE_OVF)
   ### #add_subdirectory (Third_party/ovf/src Third_party/ovf/fortran)

   ### add_library(ovf STATIC  Third_party/ovf/src/ovf.cpp Third_party/ovf/fortran/ovf.f90)
   ### target_include_directories( ovf PUBLIC Third_party/ovf/include)
   ### set_target_properties( ovf PROPERTIES COMPILE_FLAGS "-lstdc++" LINK_FLAGS "-lstdc++")
   ### target_link_libraries (${UppASD_EXE} ovf)
   ### # Don't add underscores in symbols for C-compatability
   ### #set(OVF_LIB "-L. -lovf")
   ### #message(STATUS "OVF_LIB ${OVF_LIB}")
ELSE()
   unset(OVF_LIB)
ENDIF(USE_OVF)
##############################################################################################


##############################################################################################
# Include VSL preprocessor flag (recommended)
##############################################################################################
if(USE_VSL AND MKL_FOUND)
   add_compile_definitions(USE_VSL)
endif()
##############################################################################################


##############################################################################################
# Check for FFT treatment (MKL default, FFTW available)
##############################################################################################
if(USE_MKL_FFT AND MKL_FOUND)
   add_sources(Hamiltonian/fftdipole_mkl.f90)
   add_compile_definitions(USE_MKL_FFT)
   include_directories(${BLA_INCLUDE_DIR})
   message(STATUS "FFT-> ${USE_MKL_FFT} ")
elseif(USE_FFTW)
   add_sources(Hamiltonian/fftdipole_fftw.f90)
   add_compile_definitions(USE_FFTW)
   message(STATUS "FFT-> ${USE_FFTW} ")
   include(${CMAKE_MODULE_PATH}/FindFFTW.cmake)
   include_directories(${FFTW_INCLUDE_DIR})
   message(STATUS "FFTW: -I${FFTW_INCLUDE_DIR}  -L${FFTW_LIBRARIES} ${FFTW_LIB}")
endif()
##############################################################################################


##############################################################################################
# Add remaining sources
##############################################################################################
# Source files in root directory
add_sources(0sd.f90	gneb_driver.f90	mc_driver.f90	pt_driver.f90	sd_driver.f90	wl_driver.f90 sld_driver.f90 sld_mc_driver.f90 uppasd.f90)
# Sources in sub-directories
add_subdirectory(Clusters)
add_subdirectory(CoarseGraining)
add_subdirectory(Correlation)
add_subdirectory(Evolution)
add_subdirectory(Fields)
add_subdirectory(GNEB)
add_subdirectory(Hamiltonian)
add_subdirectory(Input)
add_subdirectory(KMC)
add_subdirectory(LSF)
add_subdirectory(Measurement)
add_subdirectory(MonteCarlo)
add_subdirectory(Parameters)
add_subdirectory(RNG)
add_subdirectory(SpinWaves)
add_subdirectory(System)
add_subdirectory(Tools)
add_subdirectory(Lattice)
add_subdirectory(Third_party/Gridding)
##############################################################################################


##############################################################################################
# Define the executable name
##############################################################################################
set(ASD_binary sd)
# Find the vernacular compiler name
string(REGEX MATCH gfortran|g95|ifort|pgf90 current_compiler  ${CMAKE_Fortran_COMPILER})
message(STATUS "Matched compiler:${current_compiler}.  ${CMAKE_Fortran_COMPILER}")
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
ADD_EXECUTABLE(${UppASD_EXE} ${SRCS} ${HDRS})
##############################################################################################


##############################################################################################
# Add OpenMP compiler/linker flags
##############################################################################################
IF(USE_OPENMP)
   include(FindOpenMP)
   if(OpenMP_Fortran_FOUND)
      set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${OpenMP_Fortran_FLAGS}")
      message(STATUS "OMP found, ${OpenMP_Fortran_FLAGS}.")
   endif()
   #SET_TARGET_PROPERTIES(${UppASD_EXE} PROPERTIES
   #                       COMPILE_FLAGS "${OpenMP_Fortran_FLAGS}"
   #                       LINK_FLAGS "${OpenMP_Fortran_FLAGS}")
ENDIF(USE_OPENMP)
##############################################################################################


##############################################################################################
# Add the needed libraries and special compiler flags
##############################################################################################
TARGET_LINK_LIBRARIES(${UppASD_EXE} ${BLAS_LIBRARIES}
                                ${LAPACK_LIBRARIES}
                                ${CMAKE_THREAD_LIBS_INIT} 
                                ${OVF_LIB} 
                                ${FFTW_LIB}  
                                ${CUDA_curand_LIBRARY}
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
INSTALL(TARGETS ${UppASD_EXE} RUNTIME DESTINATION bin)
