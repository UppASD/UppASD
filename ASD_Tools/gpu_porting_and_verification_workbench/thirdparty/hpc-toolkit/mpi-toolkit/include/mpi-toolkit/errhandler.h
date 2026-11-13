/*
 * SPDX-FileCopyrightText: 2026 CSC – IT Center for Science
 *
 * SPDX-License-Identifier: Apache-2.0
 */

#pragma once

#include <mpi.h>

#include "common/errhandler.h"

#define ERRCHK_MPI(x)                                                                              \
    do {                                                                                           \
        const int result = (x);                                                                    \
        if (result != MPI_SUCCESS) {                                                               \
            char error_string[MPI_MAX_ERROR_STRING];                                               \
            int  error_string_len;                                                                 \
            MPI_Error_string(result, error_string, &error_string_len);                             \
            errhandler_err(__FILE__, __LINE__, #x, result, "MPI error", error_string);             \
        }                                                                                          \
    } while (0)
