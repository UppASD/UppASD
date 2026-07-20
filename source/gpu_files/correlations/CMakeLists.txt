if(USE_HIP)
  if(ON_LUMI)
    message(STATUS "GPU Correlations: HIP with Thrust complex (ON_LUMI)")
  else()
    message(STATUS "GPU Correlations: HIP with native hipComplex")
  endif()
elseif(USE_CUDA)
  message(STATUS "GPU Correlations: CUDA with Thrust complex")
endif()

# Sources containing device code; everything else here is plain C++.
set(GPU_KERNEL_SOURCES
correlation_kernels.cpp
fortranCorrelation.cpp
gpuCorrelations.cpp
)
uppasd_mark_gpu_sources(${GPU_KERNEL_SOURCES})

target_sources(${UppASD_LIB} PRIVATE ${GPU_KERNEL_SOURCES}
gpuCorrelations.hpp
correlation_kernels.hpp
correlation.hpp
fortranCorrelation.hpp
correlationFactory.hpp
correlation_types.h
)

target_include_directories(${UppASD_LIB}
        PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}
)

