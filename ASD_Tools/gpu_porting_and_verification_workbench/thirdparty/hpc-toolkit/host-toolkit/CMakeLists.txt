# SPDX-FileCopyrightText: 2026 CSC – IT Center for Science
#
# SPDX-License-Identifier: Apache-2.0

set(SOURCES
    src/memory.c
    src/rand.c
    src/reduce.c
    src/xcorr.c
)

add_library(hpc-host-toolkit ${SOURCES})
add_library(hpc::host-toolkit ALIAS hpc-host-toolkit)
target_include_directories(hpc-host-toolkit PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_link_libraries(hpc-host-toolkit PUBLIC hpc::common)
