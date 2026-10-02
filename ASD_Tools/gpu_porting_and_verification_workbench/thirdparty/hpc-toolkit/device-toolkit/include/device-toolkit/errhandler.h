/*
 * SPDX-FileCopyrightText: 2026 CSC – IT Center for Science
 *
 * SPDX-License-Identifier: Apache-2.0
 */

#pragma once

#include "common/errhandler.h"

#if defined(HIP_ENABLED)
#include "device-toolkit/hip.h"
#else
#include <cuda_runtime.h>
#endif

#define ERRCHK_CUDA(x)                                                                             \
    do {                                                                                           \
        const cudaError_t error = (x);                                                             \
        if (error != cudaSuccess) {                                                                \
            errhandler_err(__FILE__,                                                               \
                           __LINE__,                                                               \
                           #x,                                                                     \
                           (int)error,                                                             \
                           cudaGetErrorName(error),                                                \
                           cudaGetErrorString(error));                                             \
        }                                                                                          \
    } while (0)

#define ERRCHK_CUDA_KERNEL()                                                                       \
    do {                                                                                           \
        ERRCHK_CUDA(cudaGetLastError());                                                           \
        ERRCHK_CUDA(cudaDeviceSynchronize());                                                      \
    } while (0)
