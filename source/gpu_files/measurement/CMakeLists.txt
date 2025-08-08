target_sources(${UppASD_LIB} PRIVATE
        measurable.hpp
        fortranMeasurement.cuh
        fortranMeasurement.cu
        measurementQueue.hpp
        measurementQueue.cpp
        cudaMeasurement.cuh
        cudaMeasurement.cu
        measurementFactory.cuh
        measurementWriter.cuh
        measurementWriter.cu
        kernels.cuh
        kernels.cu
        measurementData.h
)

target_include_directories(${UppASD_LIB}
        PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}
)