// SPDX-FileCopyrightText: 2026 CSC – IT Center for Science
//
// SPDX-License-Identifier: Apache-2.0

#include "device-toolkit/rand.h"

#include <stdint.h>
#include <stdio.h>
#include <type_traits>

#include "common/errhandler.h"
#include "device-toolkit/errhandler.h"

__global__ void
device_rng_init_kernel(const uint64_t seed, const size_t count, curandState* states)
{
    const size_t idx = (size_t)threadIdx.x + (size_t)blockIdx.x * (size_t)blockDim.x;
    if (idx >= count)
        return;

    curand_init(seed, idx, 0, &states[idx]);
}

device_rng_t
device_rng_create(const uint64_t seed, const size_t count)
{
    ERRCHK(count > 0);

    device_rng_t rng = {
        .count  = count,
        .states = nullptr,
    };
    ERRCHK_CUDA(cudaMalloc((void**)&rng.states, count * sizeof(rng.states[0])));
    ERRCHK(rng.states != nullptr);

    constexpr uint32_t tpb = 256;
    const uint32_t     bpg = (uint32_t)((count + tpb - 1) / tpb);
    device_rng_init_kernel<<<bpg, tpb>>>(seed, count, rng.states);
    ERRCHK_CUDA_KERNEL();

    return rng;
}

void
device_rng_destroy(device_rng_t* rng)
{
    ERRCHK(rng != nullptr);
    ERRCHK(rng->states != nullptr);

    ERRCHK_CUDA(cudaFree(rng->states));
    rng->states = nullptr;
    rng->count  = 0;
}

__global__ void
device_randomize_kernel(const size_t count, real_t* data, curandState* states)
{
    const size_t idx = (size_t)threadIdx.x + (size_t)blockIdx.x * (size_t)blockDim.x;
    if (idx >= count)
        return;

    if constexpr (std::is_same_v<real_t, double>)
        data[idx] = curand_uniform_double(&states[idx]);
    else
        data[idx] = curand_uniform(&states[idx]);
}

void
device_randomize(const size_t count, real_t* data, device_rng_t* rng)
{
    ERRCHK(count > 0);
    ERRCHK(data != nullptr);
    ERRCHK(rng != nullptr);
    ERRCHK(rng->states != nullptr);
    ERRCHK(count <= rng->count);

    constexpr uint32_t tpb = 256;
    const uint32_t     bpg = (uint32_t)((count + tpb - 1) / tpb);
    device_randomize_kernel<<<bpg, tpb>>>(count, data, rng->states);
    ERRCHK_CUDA_KERNEL();
}
