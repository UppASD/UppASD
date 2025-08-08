#pragma once

#include <cuda_runtime.h>
#include <cooperative_groups.h>
#include "real_type.h"
#include "tensor.cuh"
#include "measurementData.h"

#ifdef SINGLE_PREC
typedef float3 real3;
#else
typedef double3 real3;
#endif


namespace kernels
{
    __inline__ __device__ real warpReduceSum(real val, cooperative_groups::thread_block_tile<32> warp);


    __global__ void sumOverAtoms(const CudaTensor<real, 3> in_tensor,
                                 CudaTensor<real, 2> out_tensor);


    __global__ void averageMagnetization(const CudaTensor<real, 2> emomMSum,
                                          uint atoms,
                                          uint ensembles,
                                          AverageMagnetizationData* out);

    __global__ void binderCumulantNoEnergy(const CudaTensor<real, 2> emomMSum,
                                           uint atoms,
                                           uint ensembles,
                                           real temp,
                                           real mub,
                                           real k_bolt,
                                           BinderCumulantData* out);

    __global__ void grad_moments(const CudaTensor<real, 3> emomM,
                                 const CudaTensor<real, 3> dxyz_vec,
                                 const CudaTensor<int, 2> dxyz_atom,
                                 const CudaTensor<int, 1> dxyz_list,
                                 CudaTensor<real, 4> grad_mom_out);

    __global__ void pontryagin_number(const CudaTensor<real, 3> emomM,
                                      const CudaTensor<real, 4> grad_mom,
                                      uint sk_num_count,
                                      SkyrmionNumberData* out);
}




