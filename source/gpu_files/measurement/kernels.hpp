#pragma once

#include "tensor.hpp"
#include "real_type.h"
#include "measurementData.h"
#if defined(HIP_V)
#include <hip/hip_runtime.h>
#elif defined(CUDA_V)
#include <cuda_runtime.h>
//#include "device_launch_parameters.h"
#endif

namespace kernels::measurement
{

    // ---------------- helper partial structs (for two-phase kernels) ----------------
    struct AvgMPart   { real mx, my, mz, m, m2; };
    struct BinderPart { real s1, s2, s4; };
    struct SumPart    { real s; };

    // ----------------------------- kernels: declarations ----------------------------

    // Sum over atoms: out_tensor has shape (3, M).  **IMPORTANT**: out_tensor must be
    // ZEROED before launch (we do atomicAdd of per-block partials).
    // Launch suggestion (good when N ≫ M):
    //   dim3 block(256,1);
    //   dim3 grid(ceil_div(N, block.x) /*tiles over atoms*/, 3*M /*i,k pairs*/);
    //   shared = ((block.x*block.y + 31)/32) * sizeof(real)
    __global__ void sumOverAtoms(const GpuTensor<real, 3>& in_tensor,
                                 GpuTensor<real, 2>& out_tensor);

    // Average magnetization over ensembles (two-phase)
    __global__ void averageMagnetization_partial(const GpuTensor<real, 2>& emomMSum,
                                                 uint atoms,
                                                 uint ensembles,
                                                 AvgMPart* __restrict__ block_parts);

    __global__ void averageMagnetization_finalize(const AvgMPart* __restrict__ block_parts,
                                                  uint nblocks,
                                                  uint ensembles,
                                                  AverageMagnetizationData& out);

    // Binder cumulant (no energy) over ensembles (two-phase)
    __global__ void binderCumulantNoEnergy_partial(const GpuTensor<real, 2>& emomMSum,
                                                   uint atoms,
                                                   uint ensembles,
                                                   BinderPart* __restrict__ block_parts);

    __global__ void binderCumulantNoEnergy_finalize(const BinderPart* __restrict__ block_parts,
                                                    uint nblocks,
                                                    uint atoms,
                                                    uint ensembles,
                                                    real temp,
                                                    real mub,
                                                    real k_bolt,
                                                    BinderCumulantData& d);

    // Gradient of magnetization field
    __global__ void grad_moments(const GpuTensor<real, 3>& emomM,
                                 const GpuTensor<real, 3>& dxyz_vec,
                                 const GpuTensor<int, 2>& dxyz_atom,
                                 const GpuTensor<int, 1>& dxyz_list,
                                 GpuTensor<real, 4>& grad_mom);

    // Pontryagin (brute force) over (atom, ensemble): two-phase reduce
    __global__ void pontryagin_no_partial(const GpuTensor<real, 3>& emomM,
                                          const GpuTensor<real, 4>& grad_mom,
                                          SumPart* __restrict__ block_parts);

    __global__ void pontryagin_no_finalize(const SumPart* __restrict__ block_parts,
                                           uint nblocks_total,
                                           uint ensembles,
                                           uint sk_num_count,
                                           SkyrmionNumberData& d);

    // Delaunay triangulation (0-based vertex ids), simp has extents (3, nsimp)
    __global__ void delaunay_tri_tri(uint nx, uint ny, uint nz, uint nt,
                                     GpuTensor<uint, 2>& simp);

    // Triangulation Pontryagin (sum over triangles × ensembles)
    // emom: (3, Natom, Mensemble)
    // simp: (3, nsimp) with 0-based vertex ids from delaunay_tri_tri
    __global__ void pontryagin_tri_partial(const GpuTensor<real, 3>& emom,
                                           const GpuTensor<uint, 2>& simp,
                                           SumPart* __restrict__ block_parts);

    __global__ void pontryagin_tri_finalize(const SumPart* __restrict__ block_parts,
                                            uint nblocks_total,
                                            uint ensembles,
                                            uint sk_num_count,
                                            SkyrmionNumberData& d);

    // ----------------------------- small helpers (host) -----------------------------
    template<class T>
    inline uint ceil_div(T a, T b) { return (a + b - T(1)) / b; }


    inline uint nwarps(const dim3& threads)
    {
        return ceil_div(threads.x * threads.y * threads.z, 32u);
    }

} // namespace kernels::measurement
