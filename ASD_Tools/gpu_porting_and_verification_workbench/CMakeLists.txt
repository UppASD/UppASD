# SPDX-FileCopyrightText: 2026 CSC – IT Center for Science
#
# SPDX-License-Identifier: Apache-2.0

cmake_minimum_required(VERSION 3.31.0)

project(uppasd-multiscale LANGUAGES C Fortran)

set(CMAKE_C_STANDARD 23)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CUDA_STANDARD 17)
set(CMAKE_HIP_STANDARD 17)

set(CMAKE_C_STANDARD_REQUIRED ON)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CUDA_STANDARD_REQUIRED ON)
set(CMAKE_HIP_STANDARD_REQUIRED ON)

if (NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE RelWithDebInfo)
endif()

set(CMAKE_HIP_ARCHITECTURES gfx90a)

set(HPC_BUILD_TESTS OFF) # Disable HPC toolkit tests
add_subdirectory(thirdparty/hpc-toolkit)

# Must match the backend of hpc::device-toolkit, which prefers HIP over CUDA.
# The toolkit's check_language() caches CMAKE_HIP_COMPILER.
if (CMAKE_HIP_COMPILER)
    enable_language(HIP)
elseif (CMAKE_CUDA_COMPILER)
    enable_language(CUDA)
endif()

add_subdirectory(src)
