# SPDX-FileCopyrightText: 2026 CSC – IT Center for Science
#
# SPDX-License-Identifier: Apache-2.0

add_executable(basic basic.c)
target_link_libraries(basic PRIVATE hpc::common hpc::host-toolkit)
add_test(NAME basic COMMAND basic)

add_executable(indextest indextest.c)
target_link_libraries(indextest PRIVATE hpc::common)
add_test(NAME indextest COMMAND indextest)

add_executable(isclosetest isclosetest.c)
target_link_libraries(isclosetest PRIVATE hpc::common)
add_test(NAME isclosetest COMMAND isclosetest)

add_executable(verifytest verifytest.c)
target_link_libraries(verifytest PRIVATE hpc::common)
add_test(NAME verifytest COMMAND verifytest)

add_executable(xcorrtest xcorrtest.c)
target_link_libraries(xcorrtest PRIVATE hpc::common hpc::host-toolkit)
add_test(NAME xcorrtest COMMAND xcorrtest)

if (TARGET hpc::device-toolkit)
    add_executable(devicetest devicetest.c)
    target_link_libraries(devicetest PRIVATE hpc::common hpc::host-toolkit hpc::device-toolkit)
    add_test(NAME devicetest COMMAND devicetest)
endif()

if (TARGET hpc::mpi-toolkit AND TARGET hpc::device-toolkit)
    add_executable(mpitest mpitest.c)
    target_link_libraries(mpitest PRIVATE hpc::common hpc::host-toolkit hpc::device-toolkit
                                          hpc::mpi-toolkit)
    add_test(NAME mpitest COMMAND mpitest)
endif()

include(CheckLanguage)
check_language(Fortran)

if (CMAKE_Fortran_COMPILER)
    enable_language(Fortran)

    add_executable(fortrantest fortrantest.f90)
    target_link_libraries(fortrantest PRIVATE hpc::host-toolkit)
    add_test(NAME fortrantest COMMAND fortrantest)

    if (TARGET hpc::device-toolkit)
        add_executable(fortrantest-device fortrantest-device.f90)
        target_link_libraries(fortrantest-device PRIVATE hpc::device-toolkit)
        add_test(NAME fortrantest-device COMMAND fortrantest-device)
    endif()
else()
    message(STATUS "No Fortran compiler found: fortrantest will not be built")
endif()
