#pragma once

#include "gpu_wrappers.h"
#include <stdexcept>

#if defined(CUDA_V)
#include <cufft.h>
#elif defined(HIP_V)
#include <hipfft/hipfft.h>
#endif

#if defined(CUDA_V)
using GpuFftHandle = cufftHandle;
using GpuFftResult = cufftResult;
#ifdef SINGLE_PREC
using GpuFftComplex = cufftComplex;
#define GPUFFT_C2C CUFFT_C2C
#define GPUFFT_EXEC_C2C cufftExecC2C
#else
using GpuFftComplex = cufftDoubleComplex;
#define GPUFFT_C2C CUFFT_Z2Z
#define GPUFFT_EXEC_C2C cufftExecZ2Z
#endif
#define GPUFFT_SUCCESS CUFFT_SUCCESS
#define GPUFFT_FORWARD CUFFT_FORWARD
#define GPUFFT_BACKWARD CUFFT_INVERSE
#define GPUFFT_PLAN_MANY cufftPlanMany
#define GPUFFT_DESTROY cufftDestroy
#define GPUFFT_SET_STREAM cufftSetStream
#elif defined(HIP_V)
using GpuFftHandle = hipfftHandle;
using GpuFftResult = hipfftResult;
#ifdef SINGLE_PREC
using GpuFftComplex = hipfftComplex;
#define GPUFFT_C2C HIPFFT_C2C
#define GPUFFT_EXEC_C2C hipfftExecC2C
#else
using GpuFftComplex = hipfftDoubleComplex;
#define GPUFFT_C2C HIPFFT_Z2Z
#define GPUFFT_EXEC_C2C hipfftExecZ2Z
#endif
#define GPUFFT_SUCCESS HIPFFT_SUCCESS
#define GPUFFT_FORWARD HIPFFT_FORWARD
#define GPUFFT_BACKWARD HIPFFT_BACKWARD
#define GPUFFT_PLAN_MANY hipfftPlanMany
#define GPUFFT_DESTROY hipfftDestroy
#define GPUFFT_SET_STREAM hipfftSetStream
#endif

inline void assertGpuFft(GpuFftResult result) {
   if(result != GPUFFT_SUCCESS) {
      throw std::runtime_error("GPU FFT call failed");
   }
}
