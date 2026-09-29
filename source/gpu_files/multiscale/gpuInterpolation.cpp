#pragma once

#if defined(CUDA_V)
#include <curand.h>
#elif defined(HIP_V)
#include <hiprand/hiprand.h>
#endif
#include "c_headers.hpp"
#include "gpuInterpolation.hpp"
#include "gpuStructures.hpp"
#include "gpu_wrappers.h"
#include "interpolation_kernels.hpp"
#include "real_type.h"
#include "tensor.hpp"

/*
gpuLattice:

   GpuTensor<real, 3> emomM;
   GpuTensor<real, 3> emom;
   GpuTensor<real, 3>  emom2;
   GpuTensor<real, 2>  mmom;
   GpuTensor<real, 2>  mmom0;
   GpuTensor<real, 2>  mmom2;
   GpuTensor<real, 2>  mmomi;

deviceInterpolationInfo:
    unsigned int nrInterpAtoms; // Number of atoms affected
    GpuTensor<int, 1> indices;   //index on firstNeighbour corresponding to each atom. indices(i)
contains 0 if the atom is not affected by the interpolation. GpuTensor<int, 1>  firstNeighbour; //
index of the first neighbour in weights and neighbours GpuTensor<real, 1>  weights; // Per-neighbour
coefficient GpuTensor<int, 1>  neighbours; // Atom indices for neighbours participating in the
interpolation

*/

GpuInterpolation::GpuInterpolation(const unsigned int p_N, const unsigned int p_M,
                                   deviceInterpolationInfo& p_gpuInterpolationInfo,
                                   deviceLattice&           p_gpuLattice)
    : N(p_N), M(p_M), gpuInterpolationInfo(p_gpuInterpolationInfo), gpuLattice(p_gpuLattice)
{
    blocks  = {1, 1, 1};
    threads = {1, 1, 1};
}

GpuInterpolation::~GpuInterpolation() {}

void
GpuInterpolation::interpolate()
{
}
