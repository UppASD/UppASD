# SPDX-FileCopyrightText: 2026 CSC – IT Center for Science
#
# SPDX-License-Identifier: Apache-2.0

include(CheckLanguage)
check_language(HIP)
check_language(CUDA)

set(DEVICE_SOURCES
    src/memory.cu
    src/rand.cu
    src/reduce.cu
)

if (CMAKE_HIP_COMPILER)
    set(DEVICE_TOOLKIT_AMDGPU_TARGETS "gfx90a;gfx942")
    set(CMAKE_HIP_ARCHITECTURES ${DEVICE_TOOLKIT_AMDGPU_TARGETS})

    # Workaround until CMake's HIP support handles GPU targets out of the box
    set(AMDGPU_TARGETS ${DEVICE_TOOLKIT_AMDGPU_TARGETS})
    set(GPU_TARGETS ${DEVICE_TOOLKIT_AMDGPU_TARGETS})

    enable_language(HIP)
    find_package(hip REQUIRED)
    find_package(hiprand REQUIRED)
    find_package(hipcub REQUIRED)

    add_library(hpc-device-toolkit STATIC ${DEVICE_SOURCES})
    set_source_files_properties(${DEVICE_SOURCES} PROPERTIES LANGUAGE HIP)
    target_link_libraries(hpc-device-toolkit PUBLIC hip::host hip::device hip::hiprand hip::hipcub)
    target_compile_definitions(hpc-device-toolkit PUBLIC HIP_ENABLED)
elseif (CMAKE_CUDA_COMPILER)
    set(CMAKE_CUDA_ARCHITECTURES 90)

    enable_language(CUDA)
    find_package(CUDAToolkit REQUIRED)

    add_library(hpc-device-toolkit STATIC ${DEVICE_SOURCES})
    target_link_libraries(hpc-device-toolkit PUBLIC CUDA::cudart CUDA::curand)
else()
    message(STATUS "No HIP or CUDA compiler found: device-toolkit will not be built")
endif()

if (TARGET hpc-device-toolkit)
    add_library(hpc::device-toolkit ALIAS hpc-device-toolkit)
    target_include_directories(hpc-device-toolkit PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
    target_link_libraries(hpc-device-toolkit PUBLIC hpc::common)
endif()
