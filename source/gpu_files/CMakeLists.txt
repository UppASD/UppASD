if(USE_HIP OR USE_CUDA)
set(GPU_TOOLS gpuThermfield.cpp gpuThermfield.hpp)
endif()
set(GPU_PARALLEL gpuParallelizationHelper.tpp)
target_sources(${UppASD_LIB} PRIVATE c_helper.h
base.hpp
${GPU_PARALLEL}
${GPU_TOOLS}
gpuCommon.hpp
gpuDepondtIntegrator.cpp
gpuDepondtIntegrator.hpp
gpuEventPool.cpp
gpuEventPool.hpp
gpuHamiltonianBackend.hpp
gpuHamiltonianCalculations.cpp
gpuHamiltonianCalculations.hpp
gpuHamiltonianDeviceOps.hpp
gpuFftWrapper.hpp
gpuLatticeConvolutionHamiltonian.cpp
gpuLatticeConvolutionHamiltonian.hpp
gpuSimulation.cpp
gpuSimulation.hpp
gpuSDSimulation.cpp
gpuMCSimulation.cpp
#gpuMeasurement.cpp
gpuMetropolis.hpp
gpuMetropolis.cpp
gpuMetropolis_bruteforce.hpp
gpuMetropolis_bruteforce.cpp
gpuMomentUpdater.cpp
gpuMomentUpdater.hpp
gpuMultiscaleField.cpp
gpuMultiscaleField.hpp
gpuParallelizationHelper.hpp 
gpuParallelizationHelper.cpp 
fort_helper.cpp
fortranData.cpp
fortranData.hpp
gridHelper.hpp
gpu_wrappers.h
#measurementQueue.cpp
#measurementQueue.hpp
printDebug.hpp
real_type.h
stopwatch.hpp
stopwatchDeviceSync.hpp
stopwatchPool.cpp
stopwatchPool.hpp
stopwatch_fortran.cpp
tensor.hpp
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

#common.cuh
#cudaMatrix.hpp
#depondtIntegrator.cpp
#depondtIntegrator.hpp
#fortMatrix.hpp
#hamiltonianCalculations.cpp
#hamiltonianCalculations.hpp
#hostMatrix.hpp
#matrix.hpp
#mdSimulation.cpp
#mdSimulation.hpp

#momentUpdater.cpp
#momentUpdater.hpp

#randomnum.cpp
#randomnum.hpp

#thermfield.cpp
#thermfield.hpp
