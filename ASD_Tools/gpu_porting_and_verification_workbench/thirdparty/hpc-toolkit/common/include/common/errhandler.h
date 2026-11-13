/*
 * SPDX-FileCopyrightText: 2026 CSC – IT Center for Science
 *
 * SPDX-License-Identifier: Apache-2.0
 */

#pragma once

#include "common/datatypes.h"

#ifdef __cplusplus
extern "C" {
#endif

void errhandler_log(const char* file, const int line, const char* expr, const int result,
                    const char* name, const char* description);
void errhandler_warn(const char* file, const int line, const char* expr, const int result,
                     const char* name, const char* description);
void errhandler_err(const char* file, const int line, const char* expr, const int result,
                    const char* name, const char* description);
void errhandler_atexit(void (*callback)(void));

#ifdef __cplusplus
}
#endif

#define ERRCHK(x)                                                                                  \
    do {                                                                                           \
        if (!(x)) {                                                                                \
            errhandler_err(__FILE__, __LINE__, #x, (x), "Generic error", "N/A");                   \
        }                                                                                          \
    } while (0)

#define ERRCHK_HPC(x)                                                                              \
    do {                                                                                           \
        const hpc_status_t status = (x);                                                           \
        if (status != HPC_SUCCESS) {                                                               \
            errhandler_err(__FILE__, __LINE__, #x, (int)status, "HPC toolkit error", "N/A");       \
        }                                                                                          \
    } while (0)
