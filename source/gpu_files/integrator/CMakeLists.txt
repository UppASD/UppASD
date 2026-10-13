target_sources(${UppASD_LIB} PRIVATE
        depondt_kernels.hpp
        depondt_kernels.cpp
)

target_include_directories(${UppASD_LIB}
        PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}
)
