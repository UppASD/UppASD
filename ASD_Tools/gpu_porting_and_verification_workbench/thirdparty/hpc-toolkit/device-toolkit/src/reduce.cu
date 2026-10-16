// SPDX-FileCopyrightText: 2026 CSC – IT Center for Science
//
// SPDX-License-Identifier: Apache-2.0

#include "device-toolkit/reduce.h"

// Placeholder CUDA -> HIP CUB translation. Swap this out for CSC's own
// translation header once available.
#if defined(HIP_ENABLED)
#include <hipcub/hipcub.hpp>
namespace cub = hipcub;
#else
#include <cub/cub.cuh>
#endif

#include "common/errhandler.h"
#include "device-toolkit/errhandler.h"

device_reduce_workspace_t
device_reduce_workspace_create(const size_t bytes)
{
    ERRCHK(bytes > 0);

    device_reduce_workspace_t workspace = {
        .bytes = bytes,
        .data  = nullptr,
    };
    ERRCHK_CUDA(cudaMalloc(&workspace.data, bytes));
    ERRCHK(workspace.data != nullptr);

    return workspace;
}

void
device_reduce_workspace_destroy(device_reduce_workspace_t* workspace)
{
    ERRCHK(workspace != nullptr);
    ERRCHK(workspace->data != nullptr);

    ERRCHK_CUDA(cudaFree(workspace->data));
    workspace->data  = nullptr;
    workspace->bytes = 0;
}

size_t
device_reduce_workspace_required_bytes(const size_t count, const real_t* data, real_t* result)
{
    ERRCHK(count > 0);
    ERRCHK(data != nullptr);
    ERRCHK(result != nullptr);

    size_t bytes = 0;
    ERRCHK_CUDA(cub::DeviceReduce::Sum(nullptr, bytes, data, result, count));

    return bytes;
}

void
device_reduce_sum(const size_t count, const real_t* data, real_t* result,
                  device_reduce_workspace_t* workspace)
{
    ERRCHK(count > 0);
    ERRCHK(data != nullptr);
    ERRCHK(result != nullptr);
    ERRCHK(workspace != nullptr);
    ERRCHK(workspace->data != nullptr);

    size_t required_bytes = 0;
    ERRCHK_CUDA(cub::DeviceReduce::Sum(nullptr, required_bytes, data, result, count));
    ERRCHK(required_bytes <= workspace->bytes);

    ERRCHK_CUDA(cub::DeviceReduce::Sum(workspace->data, workspace->bytes, data, result, count));
}
