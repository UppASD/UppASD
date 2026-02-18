if(USE_HIP)
set(GPU_TOOLS gpuCorrelations.hpp gpuCorrelation.hip.cpp)
elseif(USE_CUDA)
  if(USE_CUCOMPLEX_CORRELATIONS)
    set(GPU_TOOLS gpuCorrelations.cuh gpuCorrelations_cucomplex.cu)
    message(STATUS "GPU Correlations: Using native cuDoubleComplex kernels (gpuCorrelations_cucomplex.cu)")
  else()
    set(GPU_TOOLS gpuCorrelations.cuh gpuCorrelations.cu)
    message(STATUS "GPU Correlations: Using Thrust complex kernels (gpuCorrelations.cu)")
  endif()
endif()
target_sources(${UppASD_LIB} PRIVATE 
${GPU_TOOLS}
correlation.hpp
fortranCorrelation.hpp
fortranCorrelation.cpp
#correlationQueue.hpp
#correlationQueue.cpp
correlationFactory.hpp
correlation_kernels.cuh
correlation_kernels.cu


#kernels.hpp
#kernels.cpp
)

target_include_directories(${UppASD_LIB}
        PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}
)


