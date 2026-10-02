# SPDX-FileCopyrightText: 2026 CSC – IT Center for Science
#
# SPDX-License-Identifier: Apache-2.0

find_package(MPI)

if (MPI_FOUND)
    add_library(hpc-mpi-toolkit INTERFACE)
    add_library(hpc::mpi-toolkit ALIAS hpc-mpi-toolkit)
    target_include_directories(hpc-mpi-toolkit INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/include)
    target_link_libraries(hpc-mpi-toolkit INTERFACE hpc::common MPI::MPI_C)
else()
    message(STATUS "MPI not found: mpi-toolkit will not be built")
endif()
