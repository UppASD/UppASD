# Sources containing device code; everything else here is plain C++.
set(GPU_KERNEL_SOURCES
        cpuRestMeasurement.cpp
        fortranMeasurement.cpp
        gpuMeasurement.cpp
        kernels.cpp
        autocorrelation_kernels.cpp
)
uppasd_mark_gpu_sources(${GPU_KERNEL_SOURCES})

target_sources(${UppASD_LIB} PRIVATE ${GPU_KERNEL_SOURCES}
        measurable.hpp
        cpuRestMeasurement.hpp
        fortranMeasurement.hpp
        measurementQueue.hpp
        measurementQueue.cpp
        gpuMeasurement.hpp
        measurementFactory.hpp
        measurementWriter.h
        measurementWriter.cpp
        measurementWriterHelper.h
        kernels.hpp
        measurementData.h
        memoryMeasurement.h
        memoryMeasurement.cpp
        autocorrelation_kernels.hpp
)

target_include_directories(${UppASD_LIB}
        PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}
)
