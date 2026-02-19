if(USE_HIP)
set(GPU_TOOLS gpuCorrelations.hpp gpuCorrelation.hip.cpp)
elseif(USE_CUDA)
  set(GPU_TOOLS gpuCorrelations.cuh gpuCorrelations.cu)
  if(USE_CUCOMPLEX_CORRELATIONS)
    message(STATUS "GPU Correlations: Using native cuDoubleComplex kernels (conditional compilation)")
    target_compile_definitions(${UppASD_LIB} PRIVATE USE_CUCOMPLEX_CORRELATIONS)
  else()
    message(STATUS "GPU Correlations: Using Thrust complex kernels (default)")
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


