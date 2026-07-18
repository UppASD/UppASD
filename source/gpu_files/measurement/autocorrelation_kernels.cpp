#include "autocorrelation_kernels.hpp"

__global__ void fill_spinwait(GpuTensor<real, 4> spinwait, const GpuTensor<real, 3> emom, const int taskNum, const int swIdx){
    int idx = blockIdx.x*blockDim.x + threadIdx.x;
    if(idx < taskNum){
    const uint N = static_cast<uint>(spinwait.extent(1));
    const uint M = static_cast<uint>(spinwait.extent(2));

    int sp = 3 * N;
    int mInd = idx / sp;
    int nn = (idx % sp);
    int nInd = nn / 3;
    int cInd = nn % 3;

    spinwait(cInd, nInd, mInd, swIdx) = emom(cInd, nInd, mInd);
    } 
}

__global__ void calc_autocorr_block(GpuTensor<real, 2> ac_block, const GpuTensor<real, 4> spinwait, const GpuTensor<real, 3> emom){
    auto grid = cg::this_grid();
    auto block = cg::this_thread_block();
    auto warp = cg::tiled_partition<WARPSIZE>(block);

    int lane = warp.thread_rank();
    int wid = warp.meta_group_rank();
    int wNum = warp.meta_group_size();
    int tid = grid.thread_rank();
    int tid_in_block = threadIdx.x;

    int swInd = blockIdx.y;

    int tid_in_X = blockIdx.x * blockDim.x + tid_in_block;

    int stride = gridDim.x * blockDim.x;

    // Register-based accumulation: store real and imaginary parts separately (reduced register pressure)
    real sum_nm = 0.0;
 
    unsigned int mInd, cInd, nInd, ii;
 
    __shared__ real shared_nm[32];

    unsigned int N = static_cast<unsigned int>(emom.extent(1));
    unsigned int M = static_cast<unsigned int>(emom.extent(2));
    unsigned int tasks = 3*M*N;
   
    #pragma unroll 4
    for (int id = tid_in_X; id < tasks; id += stride) {
        ii = id / 3;
        cInd = id % 3;
        nInd = ii % N;
        mInd = ii / N; 
        
        sum_nm += spinwait(cInd, nInd, mInd, swInd)*emom(cInd, nInd, mInd);
    }

    
    warp.sync();

    // Warp-level reduction using register-based functions
    warpReduceScalar(sum_nm);

    // Store warp results to shared memory
    if (lane == 0) {
        shared_nm[wid] = sum_nm;
    }

    __syncthreads();              // Wait for all partial reductions
    
    // Load results from shared memory for final warp reduction
    sum_nm = (tid_in_block < wNum) ? shared_nm[lane] : 0;


    // Final reduction in first warp
    if (wid == 0) {
        warpReduceScalar(sum_nm);
    }

    if (tid_in_block == 0) {
        ac_block(block.group_index().x, swInd) = sum_nm;
    }

}

__global__ void calc_autocorr_final(GpuTensor<real, 2> ac_block, GpuTensor<real, 2> ac, const real norm, const int ac_count, const int numBlocks){

    auto grid = cg::this_grid();
    auto block = cg::this_thread_block();
    auto warp = cg::tiled_partition<WARPSIZE>(block);

    int lane = warp.thread_rank();
    int wid = warp.meta_group_rank();
    int wNum = warp.meta_group_size();
    int tid = grid.thread_rank();
    int tNum = block.size();
    int tid_in_block = block.thread_rank();

    int swInd = blockIdx.x;
    int tid_in_SW = tid_in_block;

    // Register-based accumulators
    real sum_nm = 0.0; 
    __shared__ real shared_nm[32];

    if (tid_in_SW < numBlocks) {
        sum_nm += ac_block(tid_in_SW, swInd);     
    }

    warp.sync();

    // Warp-level reduction
    warpReduceScalar(sum_nm);

    if (lane == 0) {
        shared_nm[wid] = sum_nm;

    }
    __syncthreads();
    sum_nm = (tid_in_block < wNum) ? shared_nm[lane] : 0;

    
    if (wid == 0) {
        warpReduceScalar(sum_nm);
    }

    if (tid_in_block == 0) {
        ac(swInd, ac_count) += (sum_nm*norm);
    }
};
