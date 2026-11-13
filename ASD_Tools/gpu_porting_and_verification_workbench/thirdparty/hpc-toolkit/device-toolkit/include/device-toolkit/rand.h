/*
 * SPDX-FileCopyrightText: 2026 CSC – IT Center for Science
 *
 * SPDX-License-Identifier: Apache-2.0
 */

#pragma once

#include <stddef.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

#include "common/datatypes.h"

#if defined(HIP_ENABLED)
#include <hip/hip_fp16.h> // Workaround for ROCm 6.3.4 issue
#include <hiprand/hiprand_kernel.h>
typedef hiprandState_t curandState;
#define curand_init hiprand_init
#define curand_uniform hiprand_uniform
#define curand_uniform_double hiprand_uniform_double
#else
#include <curand_kernel.h>
#endif

#ifdef __cplusplus
extern "C" {
#endif

// Device RNG states
typedef struct {
    size_t       count;
    curandState* states;
} device_rng_t;

device_rng_t device_rng_create(const uint64_t seed, const size_t count);
void         device_rng_destroy(device_rng_t* rng);

void device_randomize(const size_t count, real_t* data, device_rng_t* rng);

#ifdef __cplusplus
}
#endif
