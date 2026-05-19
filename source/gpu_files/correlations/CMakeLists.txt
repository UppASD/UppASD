#if(USE_HIP)
#set(GPU_TOOLS gpuCorrelations.hpp gpuCorrelation.hip.cpp)
#elseif(USE_CUDA)
#  if(USE_CUCOMPLEX_CORRELATIONS)
#    set(GPU_TOOLS gpuCorrelations.hpp gpuCorrelations_cucomplex.cu)
#    message(STATUS "GPU Correlations: Using native cuDoubleComplex kernels (gpuCorrelations_cucomplex.cu)")
#  else()
#    set(GPU_TOOLS gpuCorrelations.hpp gpuCorrelations.cu)
#    message(STATUS "GPU Correlations: Using Thrust complex kernels (gpuCorrelations.cu)")
#  endif()
#endif()
if(USE_HIP)
  if(ON_LUMI)
    message(STATUS "GPU Correlations: HIP with Thrust complex (ON_LUMI)")
  else()
    message(STATUS "GPU Correlations: HIP with native hipComplex")
  endif()
  set(GPU_TOOLS correlation_kernels.hip.cpp)
elseif(USE_CUDA)
  message(STATUS "GPU Correlations: CUDA with Thrust complex")
  set(GPU_TOOLS  correlation_kernels.cpp)
endif()
target_sources(${UppASD_LIB} PRIVATE
${GPU_TOOLS}
gpuCorrelations.hpp
correlation_kernels.hpp 
correlation.hpp
fortranCorrelation.hpp
fortranCorrelation.cpp
correlationFactory.hpp
gpuCorrelations.cpp
correlation_types.h


#kernels.hpp
#kernels.cpp
)

target_include_directories(${UppASD_LIB}
        PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}
)
