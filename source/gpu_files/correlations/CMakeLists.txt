if(USE_HIP)
set(GPU_TOOLS gpuCorrelations.hpp gpuCorrelation.hip.cpp)
elseif(USE_CUDA)
set(GPU_TOOLS gpuCorrelations.cuh gpuCorrelations.cu)
endif()
target_sources(${UppASD_LIB} PRIVATE 
${GPU_TOOLS}
correlation.hpp
fortranCorrelation.hpp
fortranCorrelation.cpp
correlationQueue.hpp
correlationQueue.cpp
correlationFactory.hpp

#kernels.hpp
#kernels.cpp
)

target_include_directories(${UppASD_LIB}
        PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}
)


