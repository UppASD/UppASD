#include "kernels.hpp"
#include <cmath>
#include "gpu_wrappers.h"



// ============================================================================
// Utilities (device)
// First index is contiguous in CudaTensor, so reading (i=0..2, j fixed) is contiguous.
// ============================================================================
namespace
{
    __device__ __forceinline__ real warp_reduce_sum(real v)
    {
        unsigned mask = 0xffffffffu;
 
        #pragma unroll
    for (int offset = WARPSIZE / 2; offset > 0; offset /= 2) {
 #if defined(HIP_V)
        v+= __shfl_down(v, offset);
 #elif defined (CUDA_V)
 #if CUDA_VERSION < 9000
        real shfl_v = __shfl_down(v, offset);
 #else
        real shfl_v = __shfl_down_sync(mask, v, offset);
 #endif
       v += shfl_v;
 #endif
    }        
        return v;
    }

    // 1D block-wide reduce; shared must have >= number of warps elements
    __device__ __forceinline__ real block_reduce_sum_1d(real v, real* shared)
    {
        v = warp_reduce_sum(v);
        int lane  = threadIdx.x & (WARPSIZE - 1);       
        int wid   = threadIdx.x >> (WARPSIZE == 32 ? 5 : 6); 
        int nwarp = (blockDim.x + WARPSIZE - 1) / WARPSIZE;
        if (lane == 0) shared[wid] = v;
        __syncthreads();
        real out = (threadIdx.x < nwarp) ? shared[threadIdx.x] : real(0);
        if (wid == 0) out = warp_reduce_sum(out);
        return out; // valid when threadIdx.x==0
    }

    __device__ __forceinline__ real vnorm3(real x, real y, real z) { return sqrt(x*x + y*y + z*z); }

} // anon


// ============================================================================
// sumOverAtoms (two-phase across ensembles)
// ============================================================================
__global__ void kernels::measurement::sumOverAtoms_partial(const GpuTensor<real, 3> in_tensor,
        GpuTensor<real, 3> block_parts)
{
    const uint N_atoms     = static_cast<uint>(in_tensor.extent(1));
    const uint M_ensembles = static_cast<uint>(in_tensor.extent(2));

    const uint q = blockIdx.y;      // 0..(3*M-1)
    const uint i = q % 3u;          // component
    const uint k = q / 3u;          // ensemble
    if (k >= M_ensembles) return;

    real local{};
    for (uint j = blockIdx.x * blockDim.x + threadIdx.x; j < N_atoms; j += blockDim.x * gridDim.x)
        local += in_tensor(i, j, k);


    extern __shared__ real s[];
    real block_sum = block_reduce_sum_1d(local, s);

    if (threadIdx.x == 0)
        block_parts(blockIdx.x, i, k) = block_sum;
}


__global__ void kernels::measurement::sumOverAtoms_finalize(const GpuTensor<real, 3> block_parts,
                                                            uint nblocks,
                                                            GpuTensor<real, 2> emomMEnsembleSums)
{
    const uint M_ensembles = static_cast<uint>(emomMEnsembleSums.extent(1));

    const uint q = blockIdx.x;      // 0 .. (3*M_ensembles - 1)
    const uint i = q % 3u;          // component
    const uint k = q / 3u;          // ensemble
    if (k >= M_ensembles) return;

    extern __shared__ real s[];

    real local = 0;
    for (uint bid = threadIdx.x; bid < nblocks; bid += blockDim.x)
        local += block_parts(bid, i, k);

    real sum = block_reduce_sum_1d(local, s);

    if (threadIdx.x == 0)
        emomMEnsembleSums(i, k) = sum;
}


// ============================================================================
// averageMagnetization (two-phase across ensembles)
// ============================================================================
__global__ void kernels::measurement::averageMagnetization_partial(const GpuTensor<real, 2> emomMSum,
                                                      uint atoms,
                                                      uint ensembles,
                                                      AvgMPart* __restrict__ block_parts)
{
    const uint k0 = blockIdx.x * blockDim.x + threadIdx.x;
    const uint stride = blockDim.x * gridDim.x;

    real sx=0, sy=0, sz=0, sm=0, sm2=0;

    for (uint k = k0; k < ensembles; k += stride)
    {
        // emomMSum(:,k) holds sums over atoms -> divide by atoms for averages
        real mx = emomMSum(0, k) / static_cast<real>(atoms);
        real my = emomMSum(1, k) / static_cast<real>(atoms);
        real mz = emomMSum(2, k) / static_cast<real>(atoms);
        real mn = vnorm3(mx, my, mz);
        sx  += mx;  sy += my;  sz += mz;
        sm  += mn;  sm2 += mn*mn;
    }

    extern __shared__ real s[];
    real r;
    r = block_reduce_sum_1d(sx, s);  if (threadIdx.x==0) block_parts[blockIdx.x].mx = r; __syncthreads();
    r = block_reduce_sum_1d(sy, s);  if (threadIdx.x==0) block_parts[blockIdx.x].my = r; __syncthreads();
    r = block_reduce_sum_1d(sz, s);  if (threadIdx.x==0) block_parts[blockIdx.x].mz = r; __syncthreads();
    r = block_reduce_sum_1d(sm, s);  if (threadIdx.x==0) block_parts[blockIdx.x].m  = r; __syncthreads();
    r = block_reduce_sum_1d(sm2,s);  if (threadIdx.x==0) block_parts[blockIdx.x].m2 = r;
}


__global__ void kernels::measurement::averageMagnetization_finalize(const AvgMPart* __restrict__ block_parts,
                                                       uint nblocks,
                                                       uint ensembles,
                                                       AverageMagnetizationData& out)
{
    real sx=0, sy=0, sz=0, sm=0, sm2=0;
    for (uint i = threadIdx.x; i < nblocks; i += blockDim.x) {
        sx  += block_parts[i].mx;
        sy  += block_parts[i].my;
        sz  += block_parts[i].mz;
        sm  += block_parts[i].m;
        sm2 += block_parts[i].m2;
    }

    extern __shared__ real s[];
    real r;
    r = block_reduce_sum_1d(sx, s);  if (threadIdx.x==0) out.m_x = r / ensembles; __syncthreads();
    r = block_reduce_sum_1d(sy, s);  if (threadIdx.x==0) out.m_y = r / ensembles; __syncthreads();
    r = block_reduce_sum_1d(sz, s);  if (threadIdx.x==0) out.m_z = r / ensembles; __syncthreads();
    r = block_reduce_sum_1d(sm, s);  if (threadIdx.x==0) out.m   = r / ensembles; __syncthreads();
    r = block_reduce_sum_1d(sm2,s);
    if (threadIdx.x==0) {
        real Em2 = r / ensembles;
        real Em  = out.m;
        real var = fmax(Em2 - Em*Em, real(0));
        out.m_stdv = sqrt(var);
    }
}


// ============================================================================
// Binder cumulant (no energy) (two-phase across ensembles)
// ============================================================================
__global__ void kernels::measurement::binderCumulantNoEnergy_partial(const GpuTensor<real, 2> emomMSum,
                                                        uint atoms,
                                                        uint ensembles,
                                                        BinderPart* __restrict__ block_parts)
{
    const uint k0 = blockIdx.x * blockDim.x + threadIdx.x;
    const uint stride = blockDim.x * gridDim.x;

    const real atoms_inv = 1.0 / static_cast<real>(atoms);

    real s1=0, s2=0, s4=0;

    for (uint k = k0; k < ensembles; k += stride)
    {
        const real mx = emomMSum(0, k) * atoms_inv;
        const real my = emomMSum(1, k) * atoms_inv;
        const real mz = emomMSum(2, k) * atoms_inv;
        const real m  = vnorm3(mx, my, mz);
        const real m2 = m*m;
        s1 += m;  s2 += m2;  s4 += m2*m2;
    }

    extern __shared__ real s[];
    real r;
    r = block_reduce_sum_1d(s1, s); if (threadIdx.x==0) block_parts[blockIdx.x].s1 = r; __syncthreads();
    r = block_reduce_sum_1d(s2, s); if (threadIdx.x==0) block_parts[blockIdx.x].s2 = r; __syncthreads();
    r = block_reduce_sum_1d(s4, s); if (threadIdx.x==0) block_parts[blockIdx.x].s4 = r;
}


__global__ void kernels::measurement::binderCumulantNoEnergy_finalize(const BinderPart* __restrict__ block_parts,
                                                         uint nblocks,
                                                         uint atoms,
                                                         uint ensembles,
                                                         real temp,
                                                         real mub,
                                                         real k_bolt,
                                                         BinderCumulantData& d)
{
    // Strided accumulation over blocks
    real s1 = 0, s2 = 0, s4 = 0;
    for (uint i = threadIdx.x; i < nblocks; i += blockDim.x) {
        s1 += block_parts[i].s1;
        s2 += block_parts[i].s2;
        s4 += block_parts[i].s4;
    }

    // Shared memory must be at least (#warps) * sizeof(real)
    extern __shared__ real sh[];

    // Reduce each scalar with the block helper.
    // IMPORTANT: barrier between calls since the helper reuses 'sh' and has no trailing sync.
    real S1 = 0, S2 = 0, S4 = 0;

    real r = block_reduce_sum_1d(s1, sh);
    if (threadIdx.x == 0) S1 = r;
    __syncthreads();

    r = block_reduce_sum_1d(s2, sh);
    if (threadIdx.x == 0) S2 = r;
    __syncthreads();

    r = block_reduce_sum_1d(s4, sh);
    if (threadIdx.x == 0) S4 = r;
    __syncthreads();

    if (threadIdx.x == 0) {
        const real w_new = static_cast<real>(ensembles);
        const real tot   = d.cumutotw + w_new;

        const real avrgme   = S1 / w_new; // batch mean |m|
        const real avrgm2_b = S2 / w_new;
        const real avrgm4_b = S4 / w_new;

        d.avrgmcum  = (d.avrgmcum  * d.cumutotw + avrgme   * w_new) / tot;
        d.avrgm2cum = (d.avrgm2cum * d.cumutotw + avrgm2_b * w_new) / tot;
        d.avrgm4cum = (d.avrgm4cum * d.cumutotw + avrgm4_b * w_new) / tot;

        const real eps = (sizeof(real) == sizeof(double)) ? real(1e-15) : real(1e-7);
        const real m2c = (fabs(d.avrgm2cum) < eps) ? eps : d.avrgm2cum;

        d.binderc = real(1) - (d.avrgm4cum / (real(3) * m2c * m2c));

        if (temp > real(0))
            d.pmsusc = (d.avrgm2cum - d.avrgmcum * d.avrgmcum) * mub*mub * atoms / (k_bolt*k_bolt) / temp;
        else
            d.pmsusc = (d.avrgm2cum - d.avrgmcum * d.avrgmcum) * mub*mub * atoms / k_bolt;

        d.cumuw    += w_new;
        d.cumutotw  =  tot;
    }
}

// ============================================================================
// Binder cumulant (with energy) (two-phase across ensembles)
// ============================================================================
__global__ void kernels::measurement::binderCumulantEnergy_partial(const GpuTensor<real, 2> emomMSum,
                                                        const GpuTensor<real, 2> energyM,
                                                        uint atoms,
                                                        uint ensembles,
                                                        BinderEnePart* __restrict__ block_parts)
{
    const uint k0 = blockIdx.x * blockDim.x + threadIdx.x;
    const uint stride = blockDim.x * gridDim.x;
    const uint emtype = blockIdx.y;
    int tid = threadIdx.x;
    int lane   = tid & (WARPSIZE - 1);
    int warp   = tid / WARPSIZE;

    const real atoms_inv = 1.0 / static_cast<real>(atoms);

    real s1=0, s2=0, s4=0;
    //real exch = 0, tt, total = 0, total2 = 0;
   real mx, my, mz, m, m2;
   real ene, ene2;

    for (uint k = k0; k < ensembles; k += stride)
    {
        //0 - magnetizarion, 1 - exchange energy, 2 - total energy
        switch(emtype){
        case(0):
        {
            mx = emomMSum(0, k) * atoms_inv;
            my = emomMSum(1, k) * atoms_inv;
            mz = emomMSum(2, k) * atoms_inv;
            m  = vnorm3(mx, my, mz);
            m2 = m*m;
            s1 += m;  s2 += m2;  s4 += m2*m2; 
            break;

        }
  
        case(1):
        {
            ene = energyM(k, 0);
            ene2 = ene*ene;
            s1 += ene;  s2 += ene2; 
            break;
        }

        case(2):
        {
            ene = energyM(k, 5);
            ene2 = ene*ene;
            s1 += ene;  s2 += ene2; 
            break;

        }

        }
    }

    __shared__ real sum1[32];
    __shared__ real sum2[32];
    __shared__ real sum4[32];

     #pragma unroll
    for (int offset = WARPSIZE / 2; offset > 0; offset >>= 1)
    {
        s1  += SHFL_DOWN(s1, offset);
        s2  += SHFL_DOWN(s2, offset);
        s4  += SHFL_DOWN(s4, offset);

    }

    if (lane == 0)
    {
        sum1[warp]  = s1;
        sum2[warp]  = s2;
        sum4[warp]  = s4;

    }

    __syncthreads();

    int num_warps = (blockDim.x + WARPSIZE - 1) / WARPSIZE;

    if (warp == 0)
    {
        s1  = (lane < num_warps) ? sum1[lane]  : (real)0;
        s2  = (lane < num_warps) ? sum2[lane]  : (real)0;
        s4  = (lane < num_warps) ? sum4[lane]  : (real)0;


        #pragma unroll
        for (int offset = WARPSIZE / 2; offset > 0; offset >>= 1)
        {
            s1  += SHFL_DOWN(s1, offset);
            s2  += SHFL_DOWN(s2, offset);
            s4  += SHFL_DOWN(s4, offset);
        }
    }

    if (tid == 0)
    {
        BinderEnePart& part = block_parts[blockIdx.y * gridDim.x + blockIdx.x];

        part.sum[emtype]  = s1;
        part.sum2[emtype] = s2;
        part.sum4[emtype] = s4;
    }


   //if (threadIdx.x==0) {
    //    block_parts[blockIdx.x].s1 = r1; 
    //    block_parts[blockIdx.x].s2 = r2;
    //    block_parts[blockIdx.x].s4 = r4;
    //}
}


__global__ void kernels::measurement::binderCumulantEnergy_finalize(
                                                         const BinderEnePart* __restrict__ block_parts,
                                                         uint nblocks,
                                                         uint atoms,
                                                         uint ensembles,
                                                         real temp,
                                                         real mub,
                                                         real k_bolt,
                                                         real mry,
                                                         //real temprescale,
                                                         //real temprescalegrad,
                                                         BinderCumulantData& d)
{
     constexpr int N = N_BINDER_TYPES;

    __shared__ real sums[N][32];
    __shared__ real sums2[N][32];
    __shared__ real sums4[N][32];

    real local[N]  = {};
    real local2[N] = {};
    real local4[N] = {};

    int tid  = threadIdx.x;
    int lane = tid & (WARPSIZE - 1);
    int warp = tid / WARPSIZE;


    if (tid < nblocks)
    {
        #pragma unroll
        for (int t = 0; t < N; ++t)
        {
            const BinderEnePart& p =
                block_parts[t * nblocks + tid];

            local[t]  = p.sum[t];
            local2[t] = p.sum2[t];
            local4[t] = p.sum4[t];
        }
    }

    #pragma unroll
    for (int offset = WARPSIZE / 2; offset > 0; offset >>= 1)
    {
        #pragma unroll
        for (int t = 0; t < N; ++t)
        {
            local[t]  += SHFL_DOWN(local[t], offset);
            local2[t] += SHFL_DOWN(local2[t], offset);
            local4[t] += SHFL_DOWN(local4[t], offset);
        }
    }


    if (lane == 0)
    {
        #pragma unroll
        for (int t = 0; t < N; ++t)
        {
            sums[t][warp]  = local[t];
            sums2[t][warp] = local2[t];
            sums4[t][warp] = local4[t];
        }
    }

    __syncthreads();

    int num_warps = (blockDim.x + WARPSIZE - 1) / WARPSIZE;

    if (warp == 0)
    {
        #pragma unroll
        for (int t = 0; t < N; ++t)
        {
            local[t]  = (lane < num_warps) ? sums[t][lane]  : (real)0;
            local2[t] = (lane < num_warps) ? sums2[t][lane] : (real)0;
            local4[t] = (lane < num_warps) ? sums4[t][lane] : (real)0;
        }

        #pragma unroll
        for (int offset = WARPSIZE / 2; offset > 0; offset >>= 1)
        {
            #pragma unroll
            for (int t = 0; t < N; ++t)
            {
                local[t]  += SHFL_DOWN(local[t], offset);
                local2[t] += SHFL_DOWN(local2[t], offset);
                local4[t] += SHFL_DOWN(local4[t], offset);
            }
        }
    }

    if (tid == 0)    
    {
        real m1     = local[MAGCUM];
        real m2     = local2[MAGCUM];
        real m4     = local4[MAGCUM];

        real exch1  = local[EXCHCUM];
        real total1 = local[TOTALCUM];
        real total2 = local2[TOTALCUM];

        const real w_new = static_cast<real>(ensembles);
        const real tot   = d.cumutotw + w_new;

        const real avrgme   = m1 / w_new;
        const real avrgm2_b = m2 / w_new;
        const real avrgm4_b = m4 / w_new;

        const real exch_mean  = exch1 / w_new;
        const real total_mean  = total1 / w_new;
        const real total2_mean = total2 / w_new;

        d.avrgmcum  = (d.avrgmcum  * d.cumutotw + avrgme   * w_new) / tot;
        d.avrgm2cum = (d.avrgm2cum * d.cumutotw + avrgm2_b * w_new) / tot;
        d.avrgm4cum = (d.avrgm4cum * d.cumutotw + avrgm4_b * w_new) / tot;

        const real eps = (sizeof(real) == sizeof(double)) ? real(1e-15) : real(1e-7);
        const real m2c = (fabs(d.avrgm2cum) < eps) ? eps : d.avrgm2cum;

        d.binderc = real(1) - (d.avrgm4cum / (real(3) * m2c * m2c));

        d.avrgecum  = (d.avrgecum  * d.cumutotw + total_mean  * w_new) / tot;
        d.avrge2cum = (d.avrge2cum * d.cumutotw + total2_mean * w_new) / tot;
        d.avrgetcum = (d.avrgetcum * d.cumutotw + exch_mean  * w_new) / tot;

        if (temp > real(0))
        {
            d.pmsusc = (d.avrgm2cum - d.avrgmcum * d.avrgmcum)* mub * mub * atoms / (k_bolt * k_bolt) / temp;

            d.cv = (d.avrge2cum - d.avrgecum * d.avrgecum) * mry * mry * atoms  / (k_bolt * k_bolt) / (temp * temp);
                // * (temprescalegrad * temp + temprescale) / (temprescale * temprescale);
        }
        else
        {
            d.pmsusc = (d.avrgm2cum - d.avrgmcum * d.avrgmcum) * mub * mub * atoms / k_bolt;
            d.cv = real(0);
        }

        d.cumuw    += w_new;
        d.cumutotw  = tot;
        d.Navrgcum += 1;
    }
}

// ============================================================================
// grad_moments (first index contiguous; N ≫ M → set block.y=1 typically)
// ============================================================================
__global__ void kernels::measurement::grad_moments(const GpuTensor<real, 3> emomM,
                                      const GpuTensor<real, 3> dxyz_vec,
                                      const GpuTensor<int, 2> dxyz_atom,
                                      const GpuTensor<int, 1> dxyz_list,
                                      GpuTensor<real, 4> grad_mom)
{
    const uint N  = static_cast<uint>(emomM.extent(1));
    const uint M  = static_cast<uint>(emomM.extent(2));
    const uint i  = blockIdx.x * blockDim.x + threadIdx.x; // atom
    const uint kk = blockIdx.y * blockDim.y + threadIdx.y; // ensemble

    if (i >= N || kk >= M) return;

    // zero 3x3
    #pragma unroll
    for (int r=0; r<3; ++r)
        #pragma unroll
        for (int c=0; c<3; ++c)
            grad_mom(r, c, i, kk) = real(0);

    const uint nnei = static_cast<uint>(dxyz_list(i));

    for (uint jn = 0; jn < nnei; ++jn) {
        const uint j = static_cast<uint>(dxyz_atom(jn, i)) - 1u; // 1-based neighbor index
        if (j >= N) continue;

        const real dx = dxyz_vec(0, jn, i);
        const real dy = dxyz_vec(1, jn, i);
        const real dz = dxyz_vec(2, jn, i);

        const real dmx = emomM(0, j, kk) - emomM(0, i, kk);
        const real dmy = emomM(1, j, kk) - emomM(1, i, kk);
        const real dmz = emomM(2, j, kk) - emomM(2, i, kk);

        if (fabs(dx) > real(1e-12)) { const real inv = real(1)/dx;
            grad_mom(0,0,i,kk) += dmx*inv; grad_mom(1,0,i,kk) += dmy*inv; grad_mom(2,0,i,kk) += dmz*inv; }
        if (fabs(dy) > real(1e-12)) { const real inv = real(1)/dy;
            grad_mom(0,1,i,kk) += dmx*inv; grad_mom(1,1,i,kk) += dmy*inv; grad_mom(2,1,i,kk) += dmz*inv; }
        if (fabs(dz) > real(1e-12)) { const real inv = real(1)/dz;
            grad_mom(0,2,i,kk) += dmx*inv; grad_mom(1,2,i,kk) += dmy*inv; grad_mom(2,2,i,kk) += dmz*inv; }
    }

    const real inv_n = (nnei > 0) ? (real(1)/real(nnei)) : real(0);
    #pragma unroll
    for (int c=0; c<3; ++c) {
        grad_mom(0,c,i,kk) *= inv_n;
        grad_mom(1,c,i,kk) *= inv_n;
        grad_mom(2,c,i,kk) *= inv_n;
    }
}


// ============================================================================
// pontryagin_no (brute-force): two-phase across (atom, ensemble)
// Good for N ≫ M → use block=(256,1), grid=(ceil_div(N,256), M)
// ============================================================================
__global__ void kernels::measurement::pontryagin_no_partial(const GpuTensor<real, 3> emomM,
                                               const GpuTensor<real, 4> grad_mom,
                                               SumPart* __restrict__ block_parts)
{
    const uint N = static_cast<uint>(emomM.extent(1));
    const uint M = static_cast<uint>(emomM.extent(2));

    const uint i = blockIdx.x * blockDim.x + threadIdx.x; // atom
    const uint k = blockIdx.y * blockDim.y + threadIdx.y; // ensemble (usually 1)

    real local = real(0);

    if (i < N && k < M) {
        const real cvec_x = grad_mom(1,0,i,k) * grad_mom(2,1,i,k)
                            - grad_mom(2,0,i,k) * grad_mom(1,1,i,k);
        const real cvec_y = grad_mom(2,0,i,k) * grad_mom(0,1,i,k)
                            - grad_mom(0,0,i,k) * grad_mom(2,1,i,k);
        const real cvec_z = grad_mom(0,0,i,k) * grad_mom(1,1,i,k)
                            - grad_mom(1,0,i,k) * grad_mom(0,1,i,k);

        local = emomM(0,i,k)*cvec_x + emomM(1,i,k)*cvec_y + emomM(2,i,k)*cvec_z;
    }

    extern __shared__ real s[];
    real block_sum = block_reduce_sum_1d(local, s);

    if (threadIdx.x == 0) {
        const uint bid = blockIdx.y * gridDim.x + blockIdx.x;
        block_parts[bid].s = block_sum;
    }
}


__global__ void kernels::measurement::pontryagin_no_finalize(const SumPart* __restrict__ block_parts,
                                                uint nblocks_total,
                                                uint ensembles,
                                                uint sk_num_count,
                                                SkyrmionNumberData& d)
{
    real sum = real(0);
    for (uint i = threadIdx.x; i < nblocks_total; i += blockDim.x) sum += block_parts[i].s;

    extern __shared__ real s[];
    s[threadIdx.x] = sum; __syncthreads();
    for (int off = blockDim.x>>1; off>0; off>>=1) {
        if (threadIdx.x < off) s[threadIdx.x] += s[threadIdx.x+off];
        __syncthreads();
    }

    if (threadIdx.x == 0) {
        d.skyno = s[0] / (M_PI * static_cast<real>(ensembles));
        const real prev = d.skyno_avg;
        d.skyno_avg  = prev + (d.skyno - prev) / static_cast<real>(sk_num_count);
        d.skyno_stdv += (d.skyno - prev) * (d.skyno - d.skyno_avg) / static_cast<real>(sk_num_count);
    }
}


// ============================================================================
// Delaunay triangulation (namespace kernels; uint/ULL types; 0-based ids)
// ============================================================================
__global__ void kernels::measurement::delaunay_tri_tri(uint nx, uint ny, uint nz, uint nt,
                                          GpuTensor<uint, 2> simp)
{
    const uint pairs = nx * ny * nz * nt;
    const uint tid = blockIdx.x * blockDim.x + threadIdx.x;
    if (tid >= pairs) return;

    const uint it = tid % nt;
    uint t = tid / nt;
    const uint x = t % nx; t /= nx;
    const uint y = t % ny; t /= ny;
    const uint z = t;

    const uint xp1 = (x + 1u) % nx;
    const uint yp1 = (y + 1u) % ny;
    const uint xm1 = (x + nx - 1u) % nx;

    const auto vid = [nx, ny, nt](uint cx, uint cy, uint cz, uint it0) -> uint {
        const uint lin = cz * ny + cy * nx + cx;
        return nt * lin + it0;
    };

    const uint tri1 = tid;
    const uint tri2 = tid + pairs;

    // Triangle 1: [0 -> +x -> +y]
    simp(0, tri1) = vid(x,   y,   z, it);
    simp(1, tri1) = vid(xp1, y,   z, it);
    simp(2, tri1) = vid(x,   yp1, z, it);

    // Triangle 2: [0 -> +y -> +y - x]
    simp(0, tri2) = vid(x,   y,   z, it);
    simp(1, tri2) = vid(x,   yp1, z, it);
    simp(2, tri2) = vid(xm1, yp1, z, it);
}



// ======================= Triangulation Pontryagin (two-phase) =======================
// Phase A: per-block partial sums over triangles (and all ensembles in a small loop)
__global__ void kernels::measurement::pontryagin_tri_partial(const GpuTensor<real, 3> emom,
                                                const GpuTensor<uint, 2>  simp,
                                                SumPart* __restrict__ block_parts)
{
    const uint nsimp = static_cast<uint>(simp.extent(1));   // triangles
    const uint M     = static_cast<uint>(emom.extent(2));   // ensembles (usually small, often 1)

    const uint tid     = blockIdx.x * blockDim.x + threadIdx.x;
    const uint stride  = blockDim.x * gridDim.x;

    real local = real(0);

    // Grid-stride over triangles; inner loop over ensembles (M small)
    for (uint isimp = tid; isimp < nsimp; isimp += stride) {
        const uint i1 = simp(0, isimp);
        const uint i2 = simp(1, isimp);
        const uint i3 = simp(2, isimp);

        for (uint k = 0; k < M; ++k) {
            // Load m1, m2, m3 (first index contiguous → three coalesced loads per vector)
            const real m1x = emom(0, i1, k), m1y = emom(1, i1, k), m1z = emom(2, i1, k);
            const real m2x = emom(0, i2, k), m2y = emom(1, i2, k), m2z = emom(2, i2, k);
            const real m3x = emom(0, i3, k), m3y = emom(1, i3, k), m3z = emom(2, i3, k);

            // Cross and dot products
            const real cx = m2y*m3z - m2z*m3y;
            const real cy = m2z*m3x - m2x*m3z;
            const real cz = m2x*m3y - m2y*m3x;

            const real triple = m1x*cx + m1y*cy + m1z*cz;

            const real m1m2 = m1x*m2x + m1y*m2y + m1z*m2z;
            const real m1m3 = m1x*m3x + m1y*m3y + m1z*m3z;
            const real m2m3 = m2x*m3x + m2y*m3y + m2z*m3z;

            const real denom = real(1) + m1m2 + m1m3 + m2m3;
            // Optional tiny guard if needed:
            // const real eps = (sizeof(real)==8)? 1e-15 : 1e-7f;
            // const real q = real(2) * atan(triple / ((fabs(denom) < eps) ? copysign(eps, denom) : denom));

            const real q = real(2) * atan(triple / denom);
            local += q;
        }
    }

    // Block reduction → one partial per block
    extern __shared__ real s[];
    real block_sum = block_reduce_sum_1d(local, s);

    if (threadIdx.x == 0) {
        block_parts[blockIdx.x].s = block_sum;
    }
}


// Phase B: reduce partials, scale by (4π·M), update SkyrmionNumberData (Welford-like)
__global__ void kernels::measurement::pontryagin_tri_finalize(const SumPart* __restrict__ block_parts,
                                                 uint nblocks_total,
                                                 uint ensembles,
                                                 uint sk_num_count,
                                                 SkyrmionNumberData& d)
{
    real sum = real(0);
    for (uint i = threadIdx.x; i < nblocks_total; i += blockDim.x)
        sum += block_parts[i].s;

    extern __shared__ real s[];
    s[threadIdx.x] = sum;
    __syncthreads();
    for (int off = blockDim.x >> 1; off > 0; off >>= 1) {
        if (threadIdx.x < off) s[threadIdx.x] += s[threadIdx.x + off];
        __syncthreads();
    }

    if (threadIdx.x == 0) {
        // Fortran: pontryagin_tri = thesum / (4*pi) / Mensemble
        d.skyno = s[0] / (real(4) * real(M_PI) * static_cast<real>(ensembles));

        const real prev = d.skyno_avg;
        d.skyno_avg  = prev + (d.skyno - prev) / static_cast<real>(sk_num_count);
        d.skyno_stdv += (d.skyno - prev) * (d.skyno - d.skyno_avg) / static_cast<real>(sk_num_count);
    }
}

__global__ void kernels::measurement::averageEnergy_partial(const GpuTensor<real, 2> energyM,
                           uint M,
                           EnePart* __restrict__ block_parts)
{
    const int etype = blockIdx.y;

    __shared__ real ene_sums[32];
    __shared__ real ene2_sums[32];

    int tid    = threadIdx.x;
    int lane   = tid & (WARPSIZE - 1);
    int warp   = tid / WARPSIZE;

    int i0 = blockIdx.x * blockDim.x + tid;
    int stride = blockDim.x * gridDim.x;

    real ene  = 0;
    real ene2 = 0;

    // energyM is stored as {exchange, anisotropy, DM, tensor, external,
    // total}; EnergyType is ordered for the output format instead.
    int energy_column = etype;
    switch (etype) {
        case DM:     energy_column = 2; break;
        case ANI:    energy_column = 1; break;
        case EXT:    energy_column = 4; break;
        case TENSOR: energy_column = 3; break;
        default: break; // EXCH and TOTAL have matching indices.
    }

    for (int i = i0; i < M; i += stride)
    {
        real e = energyM(i, energy_column);

        ene  += e;
        ene2 += e * e;
    }

    #pragma unroll
    for (int offset = WARPSIZE / 2; offset > 0; offset >>= 1)
    {
        ene  += SHFL_DOWN(ene, offset);
        ene2 += SHFL_DOWN(ene2, offset);
    }

    if (lane == 0)
    {
        ene_sums[warp]  = ene;
        ene2_sums[warp] = ene2;
    }

    __syncthreads();

    int num_warps = (blockDim.x + WARPSIZE - 1) / WARPSIZE;

    if (warp == 0)
    {
        ene  = (lane < num_warps) ? ene_sums[lane]  : (real)0;
        ene2 = (lane < num_warps) ? ene2_sums[lane] : (real)0;

        #pragma unroll
        for (int offset = WARPSIZE / 2; offset > 0; offset >>= 1)
        {
            ene  += SHFL_DOWN(ene, offset);
            ene2 += SHFL_DOWN(ene2, offset);
        }
    }

    if (tid == 0)
    {
        EnePart& part = block_parts[blockIdx.y * gridDim.x + blockIdx.x];

        part.sum[etype]  = ene;
        part.sum2[etype] = ene2;
    }
}
    

__global__ void kernels::measurement::averageEnergy_final(const EnePart* __restrict__ block_parts,
                                                            uint nblocks_x,
                                                            uint M,
                                                            real fcinv,
                                                            EnergyData& ene)
{
    constexpr int N = N_ENERGY_TYPES;

    __shared__ real sums[N][32];
    __shared__ real sums2[N][32];

    real local[N]  = {};
    real local2[N] = {};

    int tid  = threadIdx.x;
    int lane = tid & (WARPSIZE - 1);
    int warp = tid / WARPSIZE;


    if (tid < nblocks_x)
    {
        #pragma unroll
        for (int t = 0; t < N; ++t)
        {
            // averageEnergy_partial stores one component per grid.y row.
            // Reading block_parts[tid] for every component reads only the
            // exchange row, which made Tot and all non-exchange outputs zero.
            const EnePart& p = block_parts[t * nblocks_x + tid];
            local[t]  = p.sum[t];
            local2[t] = p.sum2[t];
        }
    }

    #pragma unroll
    for (int offset = WARPSIZE / 2; offset > 0; offset >>= 1)
    {
        #pragma unroll
        for (int t = 0; t < N; ++t)
        {
            local[t]  += SHFL_DOWN(local[t], offset);
            local2[t] += SHFL_DOWN(local2[t], offset);
        }
    }


    if (lane == 0)
    {
        #pragma unroll
        for (int t = 0; t < N; ++t)
        {
            sums[t][warp]  = local[t];
            sums2[t][warp] = local2[t];
        }
    }

    __syncthreads();

    int num_warps = (blockDim.x + WARPSIZE - 1) / WARPSIZE;

    if (warp == 0)
    {
        #pragma unroll
        for (int t = 0; t < N; ++t)
        {
            local[t]  = (lane < num_warps) ? sums[t][lane]  : (real)0;
            local2[t] = (lane < num_warps) ? sums2[t][lane] : (real)0;
        }

        #pragma unroll
        for (int offset = WARPSIZE / 2; offset > 0; offset >>= 1)
        {
            #pragma unroll
            for (int t = 0; t < N; ++t)
            {
                local[t]  += SHFL_DOWN(local[t], offset);
                local2[t] += SHFL_DOWN(local2[t], offset);
            }
        }
    }

    if (tid == 0)
    {
        #pragma unroll
        for (int t = 0; t < N; ++t)
        {
            local[t]  /= static_cast<real>(M);
            local2[t] /= static_cast<real>(M);
        }

        ene.exchange  = local[EXCH]*fcinv;
        ene.DM        = local[DM]*fcinv;
        ene.anisotropy= local[ANI]*fcinv;
        ene.Zeeman    = local[EXT]*fcinv;
        ene.pair      = local[TENSOR]*fcinv;
        ene.total     = local[TOTAL]*fcinv;

        ene.std_exchange =
            sqrt(fmax(local2[EXCH] - local[EXCH]*local[EXCH], real(0)))*fcinv;

        ene.std_DM =
            sqrt(fmax(local2[DM] - local[DM]*local[DM], real(0)))*fcinv;

        ene.std_anisotropy =
            sqrt(fmax(local2[ANI] - local[ANI]*local[ANI], real(0)))*fcinv;

        ene.std_Zeeman =
            sqrt(fmax(local2[EXT] - local[EXT]*local[EXT], real(0)))*fcinv;

        ene.std_pair =
            sqrt(fmax(local2[TENSOR] - local[TENSOR]*local[TENSOR], real(0)))*fcinv;

        ene.std_total =
            sqrt(fmax(local2[TOTAL] - local[TOTAL]*local[TOTAL], real(0)))*fcinv;
    }
}

// ============================================================================
// Projected stuff
// ============================================================================
__global__ void kernels::measurement::sumOverAtomsProj_partial(const GpuTensor<real, 3> in_tensor,
        GpuTensor<real, 3> block_parts)
{
    const uint N_atoms     = static_cast<uint>(in_tensor.extent(1));
    const uint M_ensembles = static_cast<uint>(in_tensor.extent(2));

    const uint q = blockIdx.y;      // 0..(3*M-1)
    const uint i = q % 3u;          // component
    const uint k = q / 3u;          // ensemble
    if (k >= M_ensembles) return;

    real local{};
    for (uint j = blockIdx.x * blockDim.x + threadIdx.x; j < N_atoms; j += blockDim.x * gridDim.x)
        local += in_tensor(i, j, k);


    extern __shared__ real s[];
    real block_sum = block_reduce_sum_1d(local, s);

    if (threadIdx.x == 0)
        block_parts(blockIdx.x, i, k) = block_sum;
}


__global__ void kernels::measurement::sumOverAtomsProj_finalize(const GpuTensor<real, 3> block_parts,
                                                            uint nblocks,
                                                            GpuTensor<real, 2> emomMEnsembleSums)
{
    const uint M_ensembles = static_cast<uint>(emomMEnsembleSums.extent(1));

    const uint q = blockIdx.x;      // 0 .. (3*M_ensembles - 1)
    const uint i = q % 3u;          // component
    const uint k = q / 3u;          // ensemble
    if (k >= M_ensembles) return;

    extern __shared__ real s[];

    real local = 0;
    for (uint bid = threadIdx.x; bid < nblocks; bid += blockDim.x)
        local += block_parts(bid, i, k);

    real sum = block_reduce_sum_1d(local, s);

    if (threadIdx.x == 0)
        emomMEnsembleSums(i, k) = sum;
}
