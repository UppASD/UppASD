target_sources(${UppASD_LIB} PRIVATE
        measurable.hpp
        fortranMeasurement.hpp
        fortranMeasurement.cpp
        measurementQueue.hpp
        measurementQueue.cpp
        gpuMeasurement.hpp
        gpuMeasurement.cpp
        measurementFactory.hpp
        measurementWriter.h
        measurementWriter.cpp
        measurementWriterHelper.h
        kernels.hpp
        kernels.cpp
        measurementData.h
        memoryMeasurement.h
        memoryMeasurement.cpp
)

target_include_directories(${UppASD_LIB}
        PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}
)
