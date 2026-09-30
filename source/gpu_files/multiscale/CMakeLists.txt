target_sources(${UppASD_LIB} PRIVATE
        gpuDepondtIntegrator_MS.hpp
        gpuDepondtIntegrator_MS.cpp
        gpuInterpolation.hpp
        gpuInterpolation.cpp
        gpuDampingBand.hpp
        gpuDampingBand.hpp
        multiscale_kernels.hpp
        multiscale_kernels.cpp
)

target_include_directories(${UppASD_LIB}
        PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}
)
