#pragma once


#include "c_headers.hpp"
#include "tensor.hpp"
#include "real_type.h"
#include <thrust/complex.h>
#if defined(HIP_V)
#include <correlation_kernels.hpp>
#elif defined(CUDA_V)
#include <correlation_kernels.cuh>
#endif

struct blocksQW{
    unsigned int tasks;
    unsigned int x;
    unsigned int y;
    unsigned int z;
    dim3 blocks;

    blocksQW(unsigned int N, unsigned int M, unsigned int nq, unsigned int numThreads, unsigned int maxBlocks){
        tasks = 3 * N * M;
        x = std::min(((tasks + numThreads - 1) / numThreads), maxBlocks);
        y = nq;
        z = 1;
        blocks = {x, y, z};
    }

    blocksQW(unsigned int N, unsigned int M, unsigned int nq, unsigned int sc_max_nstep, unsigned int nw, unsigned int numThreads, unsigned int maxBlocks){
        tasks = 3 * sc_max_nstep;;
        x = std::min(((tasks + numThreads - 1) / numThreads), maxBlocks);
        y = nq ;
        z = nw;
        blocks = {x, y, z};
    }
};

struct blocksQWproj{
    unsigned int tasks;
    unsigned int x;
    unsigned int y;
    unsigned int z;
    unsigned int blocksNum;
    dim3 blocks;
    dim3 blocksFin;

    blocksQWproj(unsigned int N, unsigned int M, unsigned int nq, unsigned int nt, unsigned int numThreads, unsigned int maxBlocks){
        tasks = 3 * N * M;
        blocksNum =  std::min(((tasks + numThreads - 1) / numThreads), maxBlocks);
        x = blocksNum;
        y = nt;
        z = nq;
        blocks = {x, y, z};
        blocksFin = {nt, nq, 1};
    }

    blocksQWproj(unsigned int N, unsigned int M, unsigned int nq, unsigned int sc_max_nstep, unsigned int nw, unsigned int nt, unsigned int numThreads, unsigned int maxBlocks){
        tasks = 3 * sc_max_nstep;
        blocksNum = std::min(((tasks + numThreads - 1) / numThreads), maxBlocks);
        x = blocksNum * nt;
        y = nq;
        z = nw;
        blocks = {x, y, z};
        blocksFin = {nt, nq, nw};
    }
};




struct SC{
    GpuTensor<thrust::complex<real>, 2> q_block;
    GpuTensor<thrust::complex<real>, 3> w_block;
    GpuTensor<thrust::complex<real>, 2> q;
    GpuTensor<thrust::complex<real>, 3> qt;
    GpuTensor<thrust::complex<real>, 3> qw;  
    //char do_sc;

    SC(char do_sc, std::size_t p_nw, std::size_t p_nq, std::size_t p_sc_max_nstep, unsigned int numThreads, blocksQW blq, blocksQW blw){

        long int bl;
        long int nw = static_cast <long int>(p_nw);
        long int nq = static_cast <long int> (p_nq);
        long int sc_max_nstep = static_cast <long int>(p_sc_max_nstep);
        

        q_block.Allocate(static_cast <long int>(3 * blq.x), nq);
        bl = (3 * nq + numThreads - 1) / numThreads;
        setZero<2> <<<bl, numThreads>>> (q_block, 3 * blq.x * nq);


    
        if ((do_sc == 'C') || (do_sc == 'Y')) {
            q.Allocate(3, nq);           
            bl = (3 * nq + numThreads - 1) / numThreads;
            setZero<2> << <bl, numThreads >> > (q, 3 * nq);
        }
 


        if ((do_sc == 'Q') || (do_sc == 'Y')) {
            qt.Allocate(3, nq, sc_max_nstep);
            qw.Allocate(3, nq, nw);
            w_block.Allocate(static_cast <long int>(3 * blw.x), nq, nw);

            bl = (3 * nq* sc_max_nstep + numThreads - 1) / numThreads;
            setZero<3> <<<bl, numThreads>>> (qt, 3 * nq * sc_max_nstep);
            bl = (3 * nq * nw + numThreads - 1) / numThreads;
            setZero<3> <<<bl, numThreads>>> (qw, 3 * nq * nw);

            bl = (3 * blw.x * nq * nw + numThreads - 1) / numThreads;
            setZero<3> << <bl, numThreads >> > (w_block, 3 * blw.x * nq * nw);

        }
 

        

    }

    void free(char do_sc){
        q_block.Free();
            if ((do_sc == 'C') || (do_sc == 'Y')) {
                q.Free();
            }
            if ((do_sc == 'Q') || (do_sc == 'Y')) {
                qt.Free();
                qw.Free();
                w_block.Free();
            }
    }


};

struct SC_proj{
    //Proj correlations
    GpuTensor<thrust::complex<real>, 3> q_block;
    GpuTensor<thrust::complex<real>, 4> w_block;
    GpuTensor<thrust::complex<real>, 3> q;
    GpuTensor<thrust::complex<real>, 4> qt;
    GpuTensor<thrust::complex<real>, 4> qw;
    GpuVector<int> aproj;    
    //size_t NT;
    //char do_sc;

    SC_proj(char do_sc, std::size_t p_nw, std::size_t p_nq, std::size_t p_sc_max_nstep, std::size_t p_nproj, Vector<int>aproj_cpu, std::size_t p_N, std::size_t numThreads, blocksQWproj blq, blocksQWproj blw){

        int bl;
        long int nw = static_cast <long int> (p_nw);
        long int nq = static_cast <long int> (p_nq);
        long int sc_max_nstep = static_cast <long int> (p_sc_max_nstep);
        long int nproj = static_cast <long int> (p_nproj);
        long int N = static_cast <long int> (p_N);

        q_block.Allocate(static_cast <long int>(3 * blq.x), nproj, nq);
        bl = (3 * nq *nproj + numThreads - 1) / numThreads;
        setZero<3> << <bl, numThreads >> > (q_block, 3 * blq.x * nq * nproj);
        aproj.Allocate(aproj_cpu.extent(0)); 
        printf("\n EXTENTS: gpu = %i, cpu = %i\n", aproj.extent(0), aproj_cpu.extent(0));

        aproj.copy_sync(aproj_cpu);         

        if ((do_sc == 'C') || (do_sc == 'Y')) {
            q.Allocate(3, nproj, nq);           
            bl = (3 * nq * nproj + numThreads - 1) / numThreads;
            setZero<3> << <bl, numThreads >> > (q, 3 * nproj * nq);
        }
        if ((do_sc == 'Q') || (do_sc == 'Y')) {
            qt.Allocate(3, nproj, nq, sc_max_nstep);
            qw.Allocate(3, nproj, nq, nw);
            w_block.Allocate(static_cast <long int>(3 * blw.x), nproj, nq, nw);

            bl = (3 * nq* sc_max_nstep* nproj + numThreads - 1) / numThreads;
            setZero<4> << <bl, numThreads >> > (qt, 3 * nq * sc_max_nstep * nproj);
            bl = (3 * nq * nw* nproj + numThreads - 1) / numThreads;
            setZero<4> << <bl, numThreads >> > (qw, 3 * nq * nw * nproj);
            bl = (3 * blw.x * nq * nw* nproj + numThreads - 1) / numThreads;
            setZero<4> << <bl, numThreads >> > (w_block, 3 * blw.x * nq * nw* nproj);

        }
    }

        void free(char do_sc){
            q_block.Free();
            aproj.Free();
                if ((do_sc == 'C') || (do_sc == 'Y')) {
                    q.Free();
                }
                if ((do_sc == 'Q') || (do_sc == 'Y')) {
                    qt.Free();
                    qw.Free();
                    w_block.Free();
                }
        }
    
};

