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
#define GPUFFT_R2C CUFFT_R2C
#define GPUFFT_C2R CUFFT_C2R
#define GPUFFT_EXEC_R2C cufftExecR2C
#define GPUFFT_EXEC_C2R cufftExecC2R
#else
using GpuFftComplex = cufftDoubleComplex;
#define GPUFFT_R2C CUFFT_D2Z
#define GPUFFT_C2R CUFFT_Z2D
#define GPUFFT_EXEC_R2C cufftExecD2Z
#define GPUFFT_EXEC_C2R cufftExecZ2D
#endif
#define GPUFFT_SUCCESS CUFFT_SUCCESS
#define GPUFFT_FORWARD CUFFT_FORWARD
#define GPUFFT_BACKWARD CUFFT_INVERSE
#define GPUFFT_PLAN_MANY cufftPlanMany
#define GPUFFT_CREATE cufftCreate
#define GPUFFT_SET_AUTO_ALLOCATION cufftSetAutoAllocation
#define GPUFFT_MAKE_PLAN_MANY cufftMakePlanMany
#define GPUFFT_ESTIMATE_MANY cufftEstimateMany
#define GPUFFT_SET_WORK_AREA cufftSetWorkArea
#define GPUFFT_DESTROY cufftDestroy
#define GPUFFT_SET_STREAM cufftSetStream
#elif defined(HIP_V)
using GpuFftHandle = hipfftHandle;
using GpuFftResult = hipfftResult;
#ifdef SINGLE_PREC
using GpuFftComplex = hipfftComplex;
#define GPUFFT_R2C HIPFFT_R2C
#define GPUFFT_C2R HIPFFT_C2R
#define GPUFFT_EXEC_R2C hipfftExecR2C
#define GPUFFT_EXEC_C2R hipfftExecC2R
#else
using GpuFftComplex = hipfftDoubleComplex;
#define GPUFFT_R2C HIPFFT_D2Z
#define GPUFFT_C2R HIPFFT_Z2D
#define GPUFFT_EXEC_R2C hipfftExecD2Z
#define GPUFFT_EXEC_C2R hipfftExecZ2D
#endif
#define GPUFFT_SUCCESS HIPFFT_SUCCESS
#define GPUFFT_FORWARD HIPFFT_FORWARD
#define GPUFFT_BACKWARD HIPFFT_BACKWARD
#define GPUFFT_PLAN_MANY hipfftPlanMany
#define GPUFFT_CREATE hipfftCreate
#define GPUFFT_SET_AUTO_ALLOCATION hipfftSetAutoAllocation
#define GPUFFT_MAKE_PLAN_MANY hipfftMakePlanMany
#define GPUFFT_ESTIMATE_MANY hipfftEstimateMany
#define GPUFFT_SET_WORK_AREA hipfftSetWorkArea
#define GPUFFT_DESTROY hipfftDestroy
#define GPUFFT_SET_STREAM hipfftSetStream
#endif

inline void assertGpuFft(GpuFftResult result) {
   if(result != GPUFFT_SUCCESS) {
      throw std::runtime_error("GPU FFT call failed");
   }
}
