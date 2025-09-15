if(USE_HIP)
set(GPU_PARALLEL hipParallelizationHelper.tpp)
set(GPU_TOOLS hipThermfield.cpp hipThermfield.hpp)
elseif(USE_CUDA)
set(GPU_PARALLEL cudaParallelizationHelper.tpp)
set(GPU_TOOLS cudaThermfield.cu cudaThermfield.hpp)
endif()
target_sources(${UppASD_LIB} PRIVATE c_helper.h
base.hpp
${GPU_PARALLEL}
${GPU_TOOLS}
gpuCommon.hpp
gpuDepondtIntegrator.cpp
gpuDepondtIntegrator.hpp
gpuEventPool.cpp
gpuEventPool.hpp
gpuHamiltonianCalculations.cpp
gpuHamiltonianCalculations.hpp
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

#common.cuh
#cudaMatrix.hpp
#depondtIntegrator.cpp
#depondtIntegrator.hpp
#fortMatrix.hpp
#hamiltonianCalculations.cpp
#hamiltonianCalculations.hpp
#hipParallelizationHelper.hip.cpp
#hipParallelizationHelper.hpp
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




