#include "cudaAverageMagnetization.cuh"
#include <cuda.h>
#include "fortranData.hpp"
#include "cudaParallelizationHelper.hpp"
#include <cooperative_groups.h>
#include <cooperative_groups/reduce.h>

// #include "MeasurementWriter.h"

namespace cg = cooperative_groups;

// should spawn one thread per ensemble M
__global__ void naiveAverageMagnetization_kernel(const CudaTensor<real, 3> emomM,
                                          AverageMagnetizationData* d);

__inline__ __device__ real warpReduceSum_kernel(real val);

__global__ void GPUEnergySum_kernel(const CudaTensor<real, 3> spin,
                             const CudaTensor<real, 3> beff,
                             CudaTensor<real, 2> eblock,
                             uint tasks);

__global__ void GPUMagnSum_kernel(const CudaTensor<real, 3> spin,
                           CudaTensor<real, 2> mblock,
                           uint tasks);

__global__ void GPUMagnFinalSum_kernel(CudaTensor<real, 2> mblock,
                                CudaTensor<real, 3> msum,
                                uint numBlocks,
                                uint curstep);

__global__ void GPUScalarFinalSum_kernel(CudaTensor<real, 2> block_value,
                                  CudaTensor<real, 2> scalar_sum,
                                  uint numBlocks,
                                  uint curstep);


AverageMagnetization::AverageMagnetization(const CudaTensor<real, 3>& emomM, dim3 cudaThreads)
: emomM(emomM)
, mavg_buff_size(*FortranData::avrg_buff)
, N(*FortranData::Natom)
, M(*FortranData::Mensemble)
, avrg_step(*FortranData::avrg_step)
, threads(cudaThreads)
, blocks( (M + threads.x - 1) / threads.x ) // could use GridHelper for this
{
    // mavg_buff_cpu.set(FortranData::mavg_buff, 3, mavg_buff_size, M);
    // mavg_buff_cpu.AllocateHost(mavg_buff_size);

    mavg_buff.Allocate(mavg_buff_size);
    mavg_buff_cpu.AllocateHost(mavg_buff_size);
    indxb_mavg.AllocateHost(mavg_buff_size);

    mavg_block_buff.Allocate(3 * blocks.x, M);

    mavg_buff.zeros();
    mavg_buff_cpu.zeros();
    indxb_mavg.zeros();
    mavg_block_buff.zeros();
}


AverageMagnetization::~AverageMagnetization()
{
    mavg_buff.Free();
    mavg_buff_cpu.FreeHost();
    indxb_mavg.FreeHost();
    mavg_block_buff.Free();
}


void AverageMagnetization::measure(std::size_t mstep)
{
    if (mstep % avrg_step != 0)
        return;

    cudaStream_t workStream = CudaParallelizationHelper::def.getWorkStream();
    mavg_block_buff.zeros();

    // Mariia's kernels
//    GPUMagnSum<<<blocks, threads, 0, workStream>>>(emomM, mavg_block_buff, 3 * N);
//    GPUMagnFinalSum<<<M, threads, 0, workStream>>> (mavg_block_buff, mavg_buff, blocks.x, mstep);

    // My simple kernel for testing
    naiveAverageMagnetization_kernel<<<blocks, threads, 0, workStream>>>(
            emomM,
            mavg_buff.data() + bcount_mavrg
    );

    indxb_mavg(bcount_mavrg++) = static_cast<uint>(mstep);
    // cudaDeviceSynchronize(); // for printing


    if (bcount_mavrg >= mavg_buff_size)
        flushMeasurements(mstep);

}


void AverageMagnetization::flushMeasurements(std::size_t mstep)
{
    // copy to fortran since buffer is full, and reset buffer
    // mavg_buff_fortran.copy_sync(mavg_buff);
    mavg_buff_cpu.copy_sync(mavg_buff);

    for (uint i = 0; i < bcount_mavrg; ++i)
        measurementWriter.write(MeasurementType::AverageMagnetization, indxb_mavg(i), (real*)&mavg_buff_cpu[i], sizeof(AverageMagnetizationData) / sizeof(real));


    mavg_buff.zeros();
    mavg_buff_cpu.zeros();
    bcount_mavrg = 0;
    // TODO: add function call to Fortran for writing the buffer to file
}


__global__ void naiveAverageMagnetization_kernel(const CudaTensor<real, 3> emomM,
                                                AverageMagnetizationData* d)
{
    const uint k = blockIdx.x * blockDim.x + threadIdx.x;
    const uint atoms = emomM.extent(1);
    const uint ensembles = emomM.extent(2);

    if (k >= ensembles)
        return;

    real m[3] = {0};

    for (uint i = 0; i < atoms; ++i)
    {
        m[0] += emomM(0, i, k);
        m[1] += emomM(1, i, k);
        m[2] += emomM(2, i, k);
    }

    m[0] /= atoms;
    m[1] /= atoms;
    m[2] /= atoms;

    const real m_norm = norm(3, m);
    const real avrgms = pow(m_norm, 2);

    atomicAdd(&(d->m_x), m[0]);
    atomicAdd(&(d->m_y), m[1]);
    atomicAdd(&(d->m_z), m[2]);
    atomicAdd(&(d->m), m_norm);
    atomicAdd(&(d->m_stdv), avrgms);

    __syncthreads();

    if (k == 0)
    {
        d->m_x /= ensembles;
        d->m_y /= ensembles;
        d->m_z /= ensembles;
        d->m /= ensembles;
        d->m_stdv = d->m_stdv / ensembles - pow(d->m, 2);
        d->m_stdv = sqrt(max(d->m_stdv, real(0.0)));

        // for debug
//        printf("[naiveAverageMagnetization_kernel] %e, %e, %e, %e, %e\n",
//               d->m_x,
//               d->m_y,
//               d->m_z,
//               d->m,
//               d->m_stdv
//        );
    }
}



__inline__ __device__
real warpReduceSum_kernel(real val) {
//    for (int offset = warpSize / 2; offset > 0; offset /= 2)
//        val += __shfl_down(val, offset);
    return val;
}


__global__ void GPUMagnSum_kernel(const CudaTensor<real, 3> spin, CudaTensor<real, 2> mblock, uint tasks) {
    auto grid = cg::this_grid();
    auto block = cg::this_thread_block();
    auto warp = cg::tiled_partition<32>(block);

    uint lane = warp.thread_rank();
    uint wid = warp.meta_group_rank();
    uint wSize = warp.size();
    uint wNum = warp.meta_group_size();
    uint tid = grid.thread_rank();
    uint tid_in_block = block.thread_rank();

    uint mInd = grid.block_index().y;
    uint offsetM = mInd * tasks;
    uint tid_in_M = grid.block_index().x * block.num_threads() + tid_in_block;
    uint stride = grid.dim_blocks().x * block.num_threads();

    real mySum[3] = { 0.0, 0.0, 0.0 };
    static __shared__ real shared0[32];
    static __shared__ real shared1[32];
    static __shared__ real shared2[32];


    for (uint id = tid_in_M; id < tasks; id += stride) {
        mySum[id % 3] += spin[id + offsetM];
        //  printf("tid = %i, mInd = %i, stride = %i, data_id = %i, mySum = %.3f\n", tid, mInd , stride, id + offsetM, mySum[id % 3]);
    }
    warp.sync();

    mySum[0] = warpReduceSum_kernel(mySum[0]);
    mySum[1] = warpReduceSum_kernel(mySum[1]);
    mySum[2] = warpReduceSum_kernel(mySum[2]);

    if (lane == 0) {
        shared0[wid] = mySum[0];
        shared1[wid] = mySum[1];
        shared2[wid] = mySum[2];
    }

    __syncthreads();              // Wait for all partial reductions
    mySum[0] = (tid_in_block < wNum) ? shared0[lane] : 0;
    mySum[1] = (tid_in_block < wNum) ? shared1[lane] : 0;
    mySum[2] = (tid_in_block < wNum) ? shared2[lane] : 0;

    if (wid == 0) {
        mySum[0] = warpReduceSum_kernel(mySum[0]); //Final reduce within first warp
        mySum[1] = warpReduceSum_kernel(mySum[1]);
        mySum[2] = warpReduceSum_kernel(mySum[2]);
    }

    if (tid_in_block == 0) {
        mblock(block.group_index().x, mInd) = mySum[0];
        mblock(block.group_index().x + grid.group_dim().x, mInd) = mySum[1];
        mblock(block.group_index().x + 2 * grid.group_dim().x, mInd) = mySum[2];
        // printf("mInd = %i, bid = %i, mblock0 = %.3f\n", mInd, block.group_index().x, mblock(block.group_index().x, mInd));
        //printf("tid = %i, mInd = %i, mblock0 = %lf, mblock1 = %lf, mblock2 = %lf\n", tid, mInd, mblock(block.group_index().x, mInd), mblock(block.group_index().x + grid.group_dim().x, mInd), mblock(block.group_index().x + 2 * grid.group_dim().x, mInd));
    }
}


__global__ void GPUMagnFinalSum_kernel(CudaTensor<real, 2> mblock, CudaTensor<real, 3> msum, uint numBlocks, uint curstep)
{
    auto grid = cg::this_grid();
    auto block = cg::this_thread_block();
    auto warp = cg::tiled_partition<32>(block);

    uint lane = warp.thread_rank();
    uint wid = warp.meta_group_rank();
    uint wSize = warp.size();
    uint wNum = warp.meta_group_size();
    uint tid = grid.thread_rank();
    uint tNum = block.size();
    uint tid_in_block = block.thread_rank();

    uint mInd = grid.block_index().x;
    uint offsetM = mInd * numBlocks;
    uint tid_in_M = tid_in_block;
    //printf("numblocks = %i\n", numBlocks);

    real mySum[3] = { 0.0, 0.0, 0.0 };
    static __shared__ real shared0[32];
    static __shared__ real shared1[32];
    static __shared__ real shared2[32];

    if (tid_in_M < numBlocks) {
        mySum[0] += mblock(tid_in_M, mInd);
        mySum[1] += mblock(tid_in_M + numBlocks, mInd);
        mySum[2] += mblock(tid_in_M + 2 * numBlocks, mInd);
        //printf("tid_in_m = %i, mInd = %i, mblock = %.3f\n", tid_in_M, mInd, mblock(tid_in_M, mInd));
    }

    warp.sync();

    mySum[0] = warpReduceSum_kernel(mySum[0]); //Final reduce within first warp
    mySum[1] = warpReduceSum_kernel(mySum[1]); //Final reduce within first warp
    mySum[2] = warpReduceSum_kernel(mySum[2]); //Final reduce within first warp

    if (lane == 0) {
        shared0[wid] = mySum[0];
        shared1[wid] = mySum[1];
        shared2[wid] = mySum[2];
    }

    __syncthreads();              // Wait for all partial reductions
    mySum[0] = (tid_in_block < wNum) ? shared0[lane] : 0;
    mySum[1] = (tid_in_block < wNum) ? shared1[lane] : 0;
    mySum[2] = (tid_in_block < wNum) ? shared2[lane] : 0;
    if (wid == 0) mySum[0] = warpReduceSum_kernel(mySum[0]); //Final reduce within first warp
    if (wid == 0) mySum[1] = warpReduceSum_kernel(mySum[1]); //Final reduce within first warp
    if (wid == 0) mySum[2] = warpReduceSum_kernel(mySum[2]); //Final reduce within first warp

    if (tid_in_block == 0) {
        msum(0, curstep, mInd) = mySum[0];
        msum(1, curstep, mInd) = mySum[1];
        msum(2, curstep, mInd) = mySum[2];
        /*mblock_gpu[block.group_index().x] += mySum[0];
        mblock_gpu[block.group_index().x + grid.group_dim().x] += mySum[1];
        mblock_gpu[block.group_index().x + 2 * grid.group_dim().x] += mySum[2];*/
        printf("mInd = %i, mblock0 = %lf, mblock1 = %lf, mblock2 = %lf\n", mInd, msum(0, curstep, mInd), msum(1, curstep, mInd), msum(2, curstep, mInd));
    }
}

__global__ void GPUEnergySum(const CudaTensor<real, 3> spin, const CudaTensor<real, 3> beff, CudaTensor<real, 2> eblock, uint tasks) {
    auto grid = cg::this_grid();
    auto block = cg::this_thread_block();
    auto warp = cg::tiled_partition<32>(block);

    int lane = warp.thread_rank();
    int wid = warp.meta_group_rank();
    int wSize = warp.size();
    int wNum = warp.meta_group_size();
    int tid = grid.thread_rank();
    int tid_in_block = block.thread_rank();

    int mInd = grid.block_index().y;
    int offsetM = mInd * tasks;
    int tid_in_M = grid.block_index().x * block.num_threads() + tid_in_block;
    int stride = grid.dim_blocks().x * block.num_threads();

    real mySum = 0.0;
    //real curSum = 0.0;
    static __shared__ real shared0[32];

    for (int id = tid_in_M; id < tasks; id += stride) {
        mySum += spin[id + offsetM] * beff[id + offsetM];
        //curSum = spin[id + offsetM]*beff[id + offsetM];
        //printf("tid = %i, site = %i, comp = %i, beff = %lf, spin = %lf, curE = %lf, sumE = %lf\n", id, id/3, id%3, beff[id + offsetM], spin[id + offsetM], curSum, mySum);
    }
    warp.sync();
    //if (tid < tasks) printf("before warp, tid = %i, sumE = %lf\n", tid, mySum);
    mySum = warpReduceSum_kernel(mySum); //Final reduce within first warp
    //if (tid < tasks) printf("after warp, tid = %i, sumE = %lf\n", tid, mySum);
    if (lane == 0) {
        shared0[wid] = mySum;
    }
    __syncthreads();              // Wait for all partial reductions
    mySum = (tid_in_block < wNum) ? shared0[lane] : 0;
    if (wid == 0) mySum = warpReduceSum_kernel(mySum); //Final reduce within first warp
    if (tid_in_block == 0) {
        eblock(block.group_index().x, mInd) = mySum;
        // printf("mInd = %i, block = %i,eblock = %lf\n", mInd, block.group_index().x, eblock(block.group_index().x, mInd));
    }
}

__global__ void GPUScalarFinalSum_kernel(CudaTensor<real, 2> block_value, CudaTensor<real, 2> scalar_sum, uint numBlocks, uint curstep)
{
    auto grid = cg::this_grid();
    auto block = cg::this_thread_block();
    auto warp = cg::tiled_partition<32>(block);

    int lane = warp.thread_rank();
    int wid = warp.meta_group_rank();
    int wSize = warp.size();
    int wNum = warp.meta_group_size();
    int tid = grid.thread_rank();
    int tNum = block.size();
    int tid_in_block = block.thread_rank();

    int mInd = grid.block_index().x;
    int offsetM = mInd * numBlocks;
    int tid_in_M = tid_in_block;

    real mySum = 0.0;
    static __shared__ real shared0[32];

    if (tid_in_M < numBlocks) { mySum += block_value(tid_in_M, mInd); }

    warp.sync();
    mySum = warpReduceSum_kernel(mySum); //Final reduce within first warp
    if (lane == 0) { shared0[wid] = mySum; }

    __syncthreads();              // Wait for all partial reductions

    mySum = (tid_in_block < wNum) ? shared0[lane] : 0;
    if (wid == 0) mySum = warpReduceSum_kernel(mySum); //Final reduce within first warp
    if (tid_in_block == 0) {

        scalar_sum(curstep, mInd) = mySum;
        printf("mInd = %i, e = %.3f\n", mInd, scalar_sum(curstep, mInd));
    }
}



//__global__ void meas::naiveAverageMagnetization(const CudaTensor<real, 3> emomM,
//                                                CudaTensor<real, 3> mavg_buff,
//                                                uint mavg_buff_i)
//{
//    const uint atom = blockIdx.x * blockDim.x + threadIdx.x;
//    if (atom >= emomM.extent(1))
//        return;
//
//    const uint ensemble = blockIdx.y * blockDim.y + threadIdx.y;
//    if (ensemble >= emomM.extent(2))
//        return;
//
//    const real px = emomM(0, atom, ensemble);
//    const real py = emomM(1, atom, ensemble);
//    const real pz = emomM(2, atom, ensemble);
//
//    atomicAdd(&mavg_buff(0, mavg_buff_i, ensemble), px);
//    atomicAdd(&mavg_buff(1, mavg_buff_i, ensemble), py);
//    atomicAdd(&mavg_buff(2, mavg_buff_i, ensemble), pz);
//
//
//    // for debug
//    __syncthreads();
//
//    if (atom == 0 && ensemble == 0)
//    {
//        printf("[cuda_kernels::measurement::averageMagnetization] mavg_buff(:, %d, %d) = {%e, %e, %e}\n",
//               mavg_buff_i, ensemble,
//               mavg_buff(0, mavg_buff_i, ensemble),
//               mavg_buff(1, mavg_buff_i, ensemble),
//               mavg_buff(2, mavg_buff_i, ensemble)
//        );
//    }
//
//}

//__global__ void CudaMeasurement::Kernels::compute_mavg_per_block(const real* emomM, uint N, uint M, real* mavg_out)
//{
//    const uint atom = blockIdx.x * blockDim.x + threadIdx.x;
//    if (atom >= N) return;
//
//    extern __shared__ real smem[];
//    real* sx = smem;
//    real* sy = &sx[blockDim.x];
//    real* sz = &sy[blockDim.x];
//
//    for (int k = 0; k < M; ++k)
//    {
//        sx[threadIdx.x] = emomM[idx3(0, atom, k, N, M)];
//        sy[threadIdx.x] = emomM[idx3(1, atom, k, N, M)];
//        sz[threadIdx.x] = emomM[idx3(2, atom, k, N, M)];
//        __syncthreads();
//
//        for (uint s = blockDim.x / 2; s > 0; s >>= 1)
//        {
//            if (threadIdx.x < s)
//            {
//                sx[threadIdx.x] += sx[threadIdx.x + s];
//                sy[threadIdx.x] += sy[threadIdx.x + s];
//                sz[threadIdx.x] += sz[threadIdx.x + s];
//            }
//            __syncthreads();
//        }
//
//        if (threadIdx.x == 0)
//        {
//            mavg_block_out[blockIdx.x].x += sx[0];
//            mavg_block_out[blockIdx.x].y += sy[0];
//            mavg_block_out[blockIdx.x].z += sz[0];
//        }
//        __syncthreads();
//
////        for (uint k = 0; k < M; ++k) {
////
////        }
//    }
//}