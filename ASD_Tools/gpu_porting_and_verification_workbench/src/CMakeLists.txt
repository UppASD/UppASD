# SPDX-FileCopyrightText: 2026 CSC – IT Center for Science
#
# SPDX-License-Identifier: Apache-2.0

# Catches out-of-range indices in the test data
set_source_files_properties(model.f90 PROPERTIES COMPILE_OPTIONS $<$<Fortran_COMPILER_ID:GNU>:-fcheck=bounds>)

add_executable(test-host test-host.c testcase.c testcase-verify.c host-implementation.c model.f90)
target_link_libraries(test-host PRIVATE hpc::common hpc::host-toolkit m)

if (TARGET hpc::device-toolkit)
    add_executable(test-device test-device.c testcase.c testcase-device.c testcase-verify.c device-implementation.cu model.f90)
    if (CMAKE_HIP_COMPILER)
        set_source_files_properties(device-implementation.cu PROPERTIES LANGUAGE HIP)
    endif()
    target_link_libraries(test-device PRIVATE hpc::common hpc::host-toolkit hpc::device-toolkit m)
endif()
