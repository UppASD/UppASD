/*
 * SPDX-FileCopyrightText: 2026 CSC – IT Center for Science
 *
 * SPDX-License-Identifier: Apache-2.0
 */

#pragma once

// Placeholder CUDA -> HIP runtime API translation. Swap this out for CSC's
// own translation header once available.

#include <hip/hip_runtime.h>

typedef hipError_t cudaError_t;

#define cudaSuccess hipSuccess
#define cudaGetErrorName hipGetErrorName
#define cudaGetErrorString hipGetErrorString

#define cudaMalloc hipMalloc
#define cudaFree hipFree

// hipHostMalloc() takes a required "flags" argument (hip_runtime_api.h is a
// plain-C header, so it has no C++ default argument to fall back on), unlike
// cudaMallocHost(). NOT verified against a real ROCm install - double check
// hipHostMallocDefault is the right flag before relying on this.
#define cudaMallocHost(ptr, size) hipHostMalloc((ptr), (size), hipHostMallocDefault)
#define cudaFreeHost hipHostFree

#define cudaMemcpy hipMemcpy
#define cudaMemcpyHostToDevice hipMemcpyHostToDevice
#define cudaMemcpyDeviceToHost hipMemcpyDeviceToHost
#define cudaMemcpyHostToHost hipMemcpyHostToHost
#define cudaMemcpyDeviceToDevice hipMemcpyDeviceToDevice

#define cudaGetLastError hipGetLastError
#define cudaDeviceSynchronize hipDeviceSynchronize
