# SPDX-FileCopyrightText: 2026 CSC – IT Center for Science
#
# SPDX-License-Identifier: Apache-2.0

# CMake settings
cmake_minimum_required(VERSION 3.31.0)

# Project settings
project(hpc-toolkit)

# Output binaries in bin/, libraries in lib/
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)

# Standards
set(CMAKE_C_STANDARD 23)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CUDA_STANDARD 17)
set(CMAKE_HIP_STANDARD 17)

set(CMAKE_C_STANDARD_REQUIRED ON)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CUDA_STANDARD_REQUIRED ON)
set(CMAKE_HIP_STANDARD_REQUIRED ON)

# Build type
if (NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE RelWithDebInfo)
endif()

# Linting
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# Project-wide compilation flags
add_compile_options(# C/CXX
                    # $<$<COMPILE_LANGUAGE:C,CXX>:-Werror>
                    $<$<COMPILE_LANGUAGE:C,CXX>:-Wfatal-errors>
                    $<$<COMPILE_LANGUAGE:C,CXX>:-Wall>
                    $<$<COMPILE_LANGUAGE:C,CXX>:-Wextra>
                    $<$<COMPILE_LANGUAGE:C,CXX>:-Wdouble-promotion>
                    $<$<COMPILE_LANGUAGE:C,CXX>:-Wfloat-conversion>
                    $<$<COMPILE_LANGUAGE:C,CXX>:-Wshadow>
                    $<$<COMPILE_LANGUAGE:C,CXX>:-Wno-unused-result>
                    $<$<COMPILE_LANGUAGE:C,CXX>:-Wconversion>
                    $<$<COMPILE_LANGUAGE:C,CXX>:-Wsign-conversion>
                    $<$<COMPILE_LANGUAGE:C,CXX>:-Wvla>
                    $<$<COMPILE_LANGUAGE:C,CXX>:-Wcast-align>
                    $<$<COMPILE_LANGUAGE:C,CXX>:-Wcast-qual>
                    $<$<COMPILE_LANGUAGE:C,CXX>:-Wfloat-equal>
                    $<$<COMPILE_LANGUAGE:C,CXX>:-Wformat=2>
                    $<$<COMPILE_LANGUAGE:C,CXX>:-Wmissing-declarations>
                    $<$<COMPILE_LANGUAGE:C,CXX>:-Wmissing-include-dirs>
                    $<$<COMPILE_LANGUAGE:C,CXX>:-Wpointer-arith>
                    $<$<COMPILE_LANGUAGE:C,CXX>:-Wredundant-decls>
                    $<$<COMPILE_LANGUAGE:C,CXX>:-Wsequence-point>
                    $<$<COMPILE_LANGUAGE:C,CXX>:-Wswitch>
                    $<$<COMPILE_LANGUAGE:C,CXX>:-Wunreachable-code>
                    $<$<COMPILE_LANGUAGE:C,CXX>:-Wunused-but-set-parameter>
                    $<$<COMPILE_LANGUAGE:C,CXX>:-Wwrite-strings>
                    $<$<COMPILE_LANGUAGE:C,CXX>:-Wpedantic>
                    $<$<COMPILE_LANGUAGE:C,CXX>:-pedantic-errors>
                    $<$<COMPILE_LANGUAGE:C,CXX>:-Werror=implicit-fallthrough>
                    $<$<COMPILE_LANGUAGE:C,CXX>:-Wundef>
                    # C
                    $<$<COMPILE_LANGUAGE:C>:-Wimplicit-function-declaration>
                    $<$<COMPILE_LANGUAGE:C>:-Wimplicit-int>
                    $<$<COMPILE_LANGUAGE:C>:-Wbad-function-cast>
                    $<$<COMPILE_LANGUAGE:C>:-Wmissing-prototypes>
                    $<$<COMPILE_LANGUAGE:C>:-Wnested-externs>
)

# Options
option(HPC_BUILD_TESTS "Build the tests" ON)

# Subdirectories
add_subdirectory(common)
add_subdirectory(host-toolkit)
add_subdirectory(device-toolkit)
add_subdirectory(mpi-toolkit)

# Testing
if (HPC_BUILD_TESTS)
    enable_testing()
    add_subdirectory(tests)
endif()
