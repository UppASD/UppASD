# Sources containing device code.  These compile as CUDA or HIP; everything
# else here is plain C++.
set(GPU_KERNEL_SOURCES
gpuThermfield.cpp
gpuDepondtIntegrator.cpp
gpuEventPool.cpp
gpuHamiltonianCalculations.cpp
gpuDipoleConvolution.cpp
gpuLatticeConvolutionHamiltonian.cpp
gpuMomentUpdater.cpp
gpuSDSimulation.cpp
gpuMCSimulation.cpp
gpuSimulation.cpp
gpuParallelizationHelper.cpp
gpuMetropolis.cpp
gpuMetropolis_bruteforce.cpp
)
uppasd_mark_gpu_sources(${GPU_KERNEL_SOURCES})

target_sources(${UppASD_LIB} PRIVATE ${GPU_KERNEL_SOURCES}
dipoleEwaldKernel.cpp
c_helper.h
base.hpp
gpuParallelizationHelper.tpp
gpuCommon.hpp
gpuDepondtIntegrator.hpp
gpuEventPool.hpp
gpuHamiltonianBackend.hpp
gpuHamiltonianCalculations.hpp
gpuDipoleConvolution.hpp
dipoleEwaldKernel.hpp
gpuHamiltonianDeviceOps.hpp
gpuFftWrapper.hpp
gpuLatticeConvolutionHamiltonian.hpp
gpuSimulation.hpp
gpuMetropolis.hpp
gpuMetropolis_bruteforce.hpp
gpuMomentUpdater.hpp
gpuMultiscaleField.cpp
gpuMultiscaleField.hpp
gpuParallelizationHelper.hpp
fort_helper.cpp
fortranData.cpp
fortranData.hpp
gridHelper.hpp
gpu_wrappers.h
printDebug.hpp
real_type.h
stopwatch.hpp
stopwatchDeviceSync.hpp
stopwatchPool.cpp
stopwatchPool.hpp
stopwatch_fortran.cpp
tensor.hpp
gpuThermfield.hpp
)

target_include_directories(${UppASD_LIB}
        PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}
)

add_subdirectory(measurement)
add_subdirectory(correlations)

if(USE_CUDA)
        target_link_libraries(${UppASD_LIB} PRIVATE CUDA::cufft)
elseif(USE_HIP)
        find_library(HIPFFT_LIBRARY
                NAMES hipfft
                PATHS /opt/rocm/lib /opt/rocm/lib64
                REQUIRED
        )
        target_link_libraries(${UppASD_LIB} PRIVATE ${HIPFFT_LIBRARY})
endif()
