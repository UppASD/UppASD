target_sources(${UppASD_LIB} PRIVATE
        cudaAverageMagnetization.cuh
        cudaAverageMagnetization.cu
        measurable.hpp
        fortranMeasurement.cuh
        fortranMeasurement.cu
        measurementQueue.hpp
        measurementQueue.cpp
        cudaMeasurement.cuh
        cudaMeasurement.cu
        measurementFactory.cuh
        cudaBinderCumulant.cuh
        cudaBinderCumulant.cu
        cudaSkyrmionNumber.cuh
        cudaSkyrmionNumber.cu
        cudaEnergy.cuh
        cudaEnergy.cu
        MeasurementWriter.h
        MeasurementWriter.cpp
)

target_include_directories(${UppASD_LIB}
        PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}
)