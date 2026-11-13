# SPDX-FileCopyrightText: 2026 CSC – IT Center for Science
#
# SPDX-License-Identifier: Apache-2.0

set(SOURCES
    src/errhandler.c
    src/debug.c
    src/stopwatch.c
    src/index.c
    src/math-toolkit.c
    src/verify.c
)

set(HPC_REAL_TYPE "double" CACHE STRING "Floating-point type of real_t")
set_property(CACHE HPC_REAL_TYPE PROPERTY STRINGS float double)

add_library(hpc-common ${SOURCES})
add_library(hpc::common ALIAS hpc-common)
target_include_directories(hpc-common PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_compile_definitions(hpc-common PUBLIC HPC_REAL_TYPE=${HPC_REAL_TYPE})
target_link_libraries(hpc-common PRIVATE m)
