#pragma once


#include "c_headers.hpp"
#include "tensor.hpp"
#include "real_type.h"
#include <thrust/complex.h>

struct blocksQ{
    unsigned int tasks;
    unsigned int x;
    unsigned int y;
    unsigned int z;
    dim3 blocks;
};

struct blocksW{
    unsigned int tasks;
    unsigned int x;
    unsigned int y;
    unsigned int z;
    dim3 blocks;
};


struct SC{
    //Proj correlations
    GpuTensor<thrust::complex<real>, 2> sc_block_gpu;
    GpuTensor<thrust::complex<real>, 3> sc_block_w_gpu;
    GpuTensor<thrust::complex<real>, 2> sc_q_gpu;
    GpuTensor<thrust::complex<real>, 3> sc_qt_gpu;
    GpuTensor<thrust::complex<real>, 3> sc_qw_gpu;  
    char do_sc;


};

struct SC_proj{
    //Proj correlations
    GpuTensor<thrust::complex<real>, 3> sc_block_gpu;
    GpuTensor<thrust::complex<real>, 4> sc_block_w_gpu;
    GpuTensor<thrust::complex<real>, 3> sc_q_gpu;
    GpuTensor<thrust::complex<real>, 4> sc_qt_gpu;
    GpuTensor<thrust::complex<real>, 4> sc_qw_gpu;
    GpuVector<int> atype;    
    //size_t NT;
    //char do_sc;

    SC(char do_sc, std::size_t p_nw, std::size_t p_nq, std::size_t p_nt, std::size_t p_nproj, blocksQ bq, blocksW bw){

        int bl;
        long int nw = static_cast <long int> p_nw;
        long int nq = static_cast <long int> p_nq;
        long int nt = static_cast <long int> p_nt;
        long int nproj = static_cast <long int> p_nproj;

        sc_block_gpu.Allocate(static_cast <long int>(3 * numBlocksX_q), static_cast <long int>(nq));
        setZero<2> << <bl, numThreads >> > (sc_block_gpu, 3 * numBlocksX_q * nq);
        atype.Allocate(nproj);           

        if ((do_sc == 'C') || (do_sc == 'Y')) {
            sc_q_gpu.Allocate(3, nproj, nq);           
            bl = (3 * nq + numThreads - 1) / numThreads;
            setZero<2> << <bl, numThreads >> > (sc_q_gpu, 3 * nproj * nq);

        }
        if ((do_sc == 'Q') || (do_sc == 'Y')) {
            sc_qt_gpu.Allocate(3, nproj, nq, nt);
            sc_qw_gpu.Allocate(3, nproj, nq, nw);
            sc_block_w_gpu.Allocate(static_cast <long int>(3 * numBlocksX_w), static_cast <long int>(nq), static_cast <long int>(nw));

            bl = (3 * nq* sc_max_nstep + numThreads - 1) / numThreads;
            setZero<3> << <bl, numThreads >> > (sc_qt_gpu, 3 * nq* sc_max_nstep);
            bl = (3 * nq * nw + numThreads - 1) / numThreads;
            setZero<3> << <bl, numThreads >> > (sc_qw_gpu, 3 * nq*nw);

        }

    }
};

