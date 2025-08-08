#include "kernels.cuh"

namespace cg = cooperative_groups;


__inline__ __device__ real kernels::warpReduceSum(real val, cg::thread_block_tile<32> warp)
{
    for (int offset = 16; offset > 0; offset /= 2)
        val += warp.shfl_down(val, offset);
    return val;
}


__global__ void kernels::sumOverAtoms(const CudaTensor<real, 3> in_tensor,
                                        CudaTensor<real, 2> out_tensor)
{
    const uint N_atoms = in_tensor.extent(1);
    const uint M_ensembles = in_tensor.extent(2);

    assert(out_tensor.extent(0) == 3);
    assert(out_tensor.extent(1) == M_ensembles);

    const uint i = blockIdx.x; // 0 <= i < 3
    const uint k = blockIdx.y; // 0 <= k < M_ensembles

    const uint tid = threadIdx.x;

    real local_sum = 0.0;

    for (uint j = tid; j < N_atoms; j += blockDim.x)
    {
        local_sum += in_tensor(i, j, k);
    }


    cg::thread_block cta = cg::this_thread_block();
    auto warp = cg::tiled_partition<32>(cta);
    local_sum = kernels::warpReduceSum(local_sum, warp);

    __shared__ real warp_sums[32];
    if (warp.thread_rank() == 0)
        warp_sums[warp.meta_group_rank()] = local_sum;
    cg::sync(cta);

    real block_sum = 0.0;
    if (warp.meta_group_rank() == 0)
    {
        if (tid < (blockDim.x / 32))
            block_sum = warp_sums[tid];
        block_sum = kernels::warpReduceSum(block_sum, warp);
        if (tid == 0)
            out_tensor(i, k) = block_sum;
    }
}


__global__ void kernels::averageMagnetization(const CudaTensor<real, 2> emomMSum,
                                                 uint atoms,
                                                 uint ensembles,
                                                 AverageMagnetizationData* out)
{
    const uint k = blockIdx.x * blockDim.x + threadIdx.x;

    if (k >= ensembles)
        return;

    real m[3] = {
        emomMSum(0, k), 
        emomMSum(1, k), 
        emomMSum(2, k)
    };

    m[0] /= atoms;
    m[1] /= atoms;
    m[2] /= atoms;

    const real m_norm = norm(3, m);
    const real avrgms = pow(m_norm, 2);

    atomicAdd(&(out->m_x), m[0]);
    atomicAdd(&(out->m_y), m[1]);
    atomicAdd(&(out->m_z), m[2]);
    atomicAdd(&(out->m), m_norm);
    atomicAdd(&(out->m_stdv), avrgms);

    __syncthreads();

    if (k == 0)
    {
        out->m_x /= ensembles;
        out->m_y /= ensembles;
        out->m_z /= ensembles;
        out->m /= ensembles;
        out->m_stdv = out->m_stdv / ensembles - pow(out->m, 2);
        out->m_stdv = sqrt(max(out->m_stdv, 0.0));
    }
}


__global__ void kernels::binderCumulantNoEnergy(const CudaTensor<real, 2> emomMSum,
                                                   uint atoms,
                                                   uint ensembles,
                                                   real temp,
                                                   real mub,
                                                   real k_bolt,
                                                   BinderCumulantData* d)
{
    const uint tid = blockIdx.x * blockDim.x + threadIdx.x;

    if (tid >= 1)
        return;

    for (uint k = 0; k < ensembles; ++k)
    {
        real m[3] = {
            emomMSum(0, k),
            emomMSum(1, k),
            emomMSum(2, k)
        };

        const real avrgme = norm(3, m) / atoms;
        const real avrgm2 = pow(avrgme, 2);
        const real avrgm4 = pow(avrgm2, 2);
        d->cumuw += 1;

        d->avrgmcum = (d->avrgmcum * d->cumutotw + avrgme * d->cumuw ) / (d->cumutotw + d->cumuw);
        d->avrgm2cum = (d->avrgm2cum * d->cumutotw + avrgm2 * d->cumuw ) / (d->cumutotw + d->cumuw);
        d->avrgm4cum = (d->avrgm4cum * d->cumutotw + avrgm4 * d->cumuw ) / (d->cumutotw + d->cumuw);

        assert(d->avrgm2cum != 0.0); // this is not checked in Fortran
        d->binderc = 1.0 - (d->avrgm4cum / 3.0 / pow(d->avrgm2cum, 2));

        // For T=0, not the proper susceptibility
    //        const real extra = (temp == 0.0)? 1.0 : k_bolt * temp;
    //        d->pmsusc = (d->avrgm2cum - pow(d->avrgmcum, 2)) * pow(mub, 2) * atoms / (k_bolt * extra);

        if (temp > 0.0)
            d->pmsusc = (d->avrgm2cum - pow(d->avrgmcum, 2)) * pow(mub, 2) * atoms / pow(k_bolt, 2) / temp;
        else
            d->pmsusc = (d->avrgm2cum - pow(d->avrgmcum, 2)) * pow(mub, 2) * atoms / k_bolt;

        d->cumutotw += d->cumuw;
    }
}



__global__ void kernels::grad_moments(const CudaTensor<real, 3> emomM,
                                        const CudaTensor<real, 3> dxyz_vec,
                                        const CudaTensor<int, 2> dxyz_atom,
                                        const CudaTensor<int, 1> dxyz_list,
                                        CudaTensor<real, 4> grad_mom)
{
    const uint N = emomM.extent(1);
    const uint M = emomM.extent(2);
    const uint iatom = blockDim.x * blockIdx.x + threadIdx.x;
    const uint kk = blockDim.y * blockIdx.y + threadIdx.y;

    if (iatom >= N || kk >= M)
        return;

    assert(dxyz_list(iatom) < N);

    for (uint i = 0; i < 3; ++i)
        for (uint j = 0; j < 3; ++j)
            grad_mom(i, j, iatom, kk) = 0;

    for (uint jneigh = 0; jneigh < dxyz_list(iatom); ++jneigh)
    {
        assert(jneigh < 26);
        const uint jatom = dxyz_atom(jneigh, iatom) - 1; // needs -1 here since it gives the index of a neighboring atom

        assert(jatom < N);
        const real d_mom[3] = {
                emomM(0, jatom, kk) - emomM(0, iatom, kk),
                emomM(1, jatom, kk) - emomM(1, iatom, kk),
                emomM(2, jatom, kk) - emomM(2, iatom, kk)
        };

        const real dv[3] = { // dv = {dx, dy, dz}
                dxyz_vec(0, jneigh, iatom),
                dxyz_vec(1, jneigh, iatom),
                dxyz_vec(2, jneigh, iatom)
        };

        for (uint coord = 0; coord < 3; ++coord)
        {
            if (abs( dv[coord] ) > 1e-7)
            {
                grad_mom(0, coord, iatom, kk) += d_mom[0] / dv[coord];
                grad_mom(1, coord, iatom, kk) += d_mom[1] / dv[coord];
                grad_mom(2, coord, iatom, kk) += d_mom[2] / dv[coord];
            }
        }
    }

    for (uint coord = 0; coord < 3; ++coord)
    {
        grad_mom(coord, 0, iatom, kk) /= dxyz_list(iatom);
        grad_mom(coord, 1, iatom, kk) /= dxyz_list(iatom);
        grad_mom(coord, 2, iatom, kk) /= dxyz_list(iatom);
    }
}


__global__ void kernels::pontryagin_number(const CudaTensor<real, 3> emomM,
                                     const CudaTensor<real, 4> grad_mom,
                                     uint sk_num_count,
                                     SkyrmionNumberData* d)
{
    const uint N = emomM.extent(1);
    const uint M = emomM.extent(2);

    const uint iatom = blockDim.x * blockIdx.x + threadIdx.x;
    const uint k = blockDim.y * blockIdx.y + threadIdx.y;

    if (iatom >= N || k >= M)
        return;

    const real cvec_x = grad_mom(1,0,iatom,k) * grad_mom(2,1,iatom,k)
                        - grad_mom(2,0,iatom,k) * grad_mom(1,1,iatom,k);

    const real cvec_y = grad_mom(2,0,iatom,k) * grad_mom(0,1,iatom,k)
                        - grad_mom(0,0,iatom,k) * grad_mom(2,1,iatom,k);

    const real cvec_z = grad_mom(0,0,iatom,k) * grad_mom(1,1,iatom,k)
                        - grad_mom(1,0,iatom,k) * grad_mom(0,1,iatom,k);


    const real partial_sum = emomM(0,iatom,k) * cvec_x
                             + emomM(1,iatom,k) * cvec_y
                             + emomM(2,iatom,k) * cvec_z;

    atomicAdd(&d->skyno, partial_sum);

    __syncthreads();

    if (iatom == 0 && k == 0)
    {
        d->skyno /= (M_PI * M);

        const real skyno_avg_prev = d->skyno_avg;
        d->skyno_avg = skyno_avg_prev + (d->skyno - skyno_avg_prev) / sk_num_count;
        d->skyno_stdv += (d->skyno - skyno_avg_prev) * (d->skyno - d->skyno_avg) / sk_num_count;
    }
}
