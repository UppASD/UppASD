target_sources(${UppASD_LIB} PRIVATE
        gpuDepondtIntegrator_MS.hpp
        gpuDepondtIntegrator_MS.cpp
        gpuInterpolation.hpp
        gpuInterpolation.cpp
        gpuDampingBand.hpp
        gpuDampingBand.hpp
        interpolation_kernels.hpp
        interpolation_kernels.cpp
)

target_include_directories(${UppASD_LIB}
        PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}
)
