#pragma once

#if defined (CUDA_V)
#include <curand.h>
#elif defined (HIP_V)
 #include <hiprand/hiprand.h>
#endif
#include "c_headers.hpp"
#include "tensor.hpp"
#include "gpuStructures.hpp"
#include "real_type.h"
#include "gpuInterpolation.hpp"
#include "gpu_wrappers.h"
#include "interpolation_kernels.hpp"



GpuInterpolation::GpuInterpolation(){


}

GpuInterpolation::~GpuInterpolation(){

}

void GpuInterpolation::interpolate(){


}

