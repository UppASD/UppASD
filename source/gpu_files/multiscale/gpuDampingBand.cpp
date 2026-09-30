#pragma once

#if defined(CUDA_V)
#include <curand.h>
#elif defined(HIP_V)
#include <hiprand/hiprand.h>
#endif
#include "c_headers.hpp"
#include "gpuDampingBand.hpp"
#include "gpuStructures.hpp"
#include "gpu_wrappers.h"
#include "multiscale_kernels.hpp"
#include "real_type.h"
#include "tensor.hpp"


/*
deviceMultiscaleRest
   GpuTensor<real, 4> backbuffer;
   int backbufferHead;


deviceDampingBand 
   deviceInterpolationInfo interpolation;
   GpuTensor<real, 1> coefficients;
   GpuTensor<real, 3> preinterpolation;
   bool  enable;

*/

GpuDampingBand::GpuDampingBand(const unsigned int p_N, const unsigned int p_M, deviceDampingBand& p_gpuDampingBand, 
                     deviceLattice& p_gpuLattice, GpuTensor<real, 4> p_backbuffer, unsigned int p_backbufferHead)
    : N(p_N)
    , M(p_M)
    , gpuDampingBand(p_gpuDampingBand)
    , gpuLattice(p_gpuLattice)
    , backbuffer(p_backbuffer)
    , backbufferHead(p_backbufferHead)
{
    blocks  = {1, 1, 1};
    threads = {1, 1, 1};
}

GpuDampingBand::~GpuDampingBand() {}

void GpuDampingBand::preInterpolate(){
}

void GpuDampingBand::corrPreInterpolation(){

}  

void GpuDampingBand::corrPreInterpolationAvg(){
    
}