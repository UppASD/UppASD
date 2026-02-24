#pragma once

#include "tensor.hpp"
#include "real_type.h"
#include <cmath>
#include "gpu_wrappers.h"
#if defined(HIP_V)
#include <hip/hip_runtime.h>
#elif defined(CUDA_V)
#include <cuda_runtime.h>
#endif

__global__ void fill_spinwait(GpuTensor<real, 4> spinwait, const GpuTensor<real, 3> emom, const int taskNum, const int swIdx);

__global__ void calc_autocorr_block(GpuTensor<real, 2> ac_block, const GpuTensor<real, 4> spinwait, const GpuTensor<real, 3> emom, const int taskNum, const int swIdx);


__global__ void calc_autocorr_final(GpuTensor<real, 2> ac_block, GpuTensor<real, 3> ac, const int taskNum, const int ac_count);
