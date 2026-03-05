#pragma once


#include "c_headers.hpp"
#include "tensor.hpp"
#include "real_type.h"
#include <thrust/complex.h>

struct blocksQW{
    unsigned int tasks;
    unsigned int x;
    unsigned int y;
    unsigned int z;
    dim3 blocks;

    blocksQW(insigned int N, unsigned int M, unsigned int nq, unsigned int numThreads, unsigned int maxBlocks){
        tasks = 3 * N * M;
        x = std::min(((tasks + numThreads - 1) / numThreads), maxBlocks);
        y = nq;
        z = 1;
        blocks = {x, y, z};
    }

    blocksQW(insigned int N, unsigned int M, unsigned int nq, unsigned int sc_max_nstep, unsigned int nw, unsigned int numThreads, unsigned int maxBlocks){
        tasks = 3 * sc_max_nstep;;
        x = std::min(((tasks + numThreads - 1) / numThreads), maxBlocks);
        y = nq * nw;
        z = 1;
        blocks = {x, y, z};
    }
};

struct blocksQWproj{
    unsigned int tasks;
    unsigned int x;
    unsigned int y;
    unsigned int z;
    dim3 blocks;

    blocksQWproj(insigned int N, unsigned int M, unsigned int nq, unsigned int nt, unsigned int numThreads, unsigned int maxBlocks){
        tasks = 3 * N * M;
        x = std::min(((tasks + numThreads - 1) / numThreads), maxBlocks);
        y = nq;
        z = 1;
        blocks = {x, y, z};
    }

    blocksQW(insigned int N, unsigned int M, unsigned int nq, unsigned int sc_max_nstep, unsigned int nw, unsigned int nt, unsigned int numThreads, unsigned int maxBlocks){
        tasks = 3 * sc_max_nstep;;
        x = std::min(((tasks + numThreads - 1) / numThreads), maxBlocks);
        y = nq * nw;
        z = 1;
        blocks = {x, y, z};
    }
};




struct SC{
    GpuTensor<thrust::complex<real>, 2> q_block;
    GpuTensor<thrust::complex<real>, 3> w_block;
    GpuTensor<thrust::complex<real>, 2> q;
    GpuTensor<thrust::complex<real>, 3> qt;
    GpuTensor<thrust::complex<real>, 3> qw;  
    //char do_sc;

    SC(char do_sc, std::size_t p_nw, std::size_t p_nq, std::size_t p_sc_max_nstep, blocksQW blq, blocksQW blw){

        int bl;
        long int nw = static_cast <long int> p_nw;
        long int nq = static_cast <long int> p_nq;
        long int sc_max_nstep = static_cast <long int> p_sc_max_nstep;

        q_block.Allocate(static_cast <long int>(3 * blq.x), nq);
        setZero<2> <<<bl, numThreads>>> (q_block, 3 * blq.x * nq);
    
        if ((do_sc == 'C') || (do_sc == 'Y')) {
            q.Allocate(3, nq);           
            bl = (3 * nq + numThreads - 1) / numThreads;
            setZero<2> << <bl, numThreads >> > (q, 3 * nq);

            qt.Allocate(3, nq, sc_max_nstep);
            qw.Allocate(3, nq, nw);
            w_block.Allocate(static_cast <long int>(3 * blw.x), nq, nw);

            bl = (3 * nq* sc_max_nstep + numThreads - 1) / numThreads;
            setZero<3> <<<bl, numThreads>>> (qt, 3 * nq * sc_max_nstep);
            bl = (3 * nq * nw + numThreads - 1) / numThreads;
            setZero<3> <<<bl, numThreads>>> (qw, 3 * nq * nw);

            bl = (3 * blw.x * nq * nw + numThreads - 1) / numThreads;
            setZero<3> << <bl, numThreads >> > (w_block, 3 * blw.x * nq * nw);
            bl = (3 * nq * nw + numThreads - 1) / numThreads;
            setZero<3> << <bl, numThreads >> > (qw, 3 * nq * nw);

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
    GpuVector<int> atype;    
    //size_t NT;
    //char do_sc;

    /*SC(char do_sc, std::size_t p_nw, std::size_t p_nq, std::size_t p_sc_max_nstep, std::size_t p_nproj, blocksQW blq, blocksQW blw){

        int bl;
        long int nw = static_cast <long int> p_nw;
        long int nq = static_cast <long int> p_nq;
        long int sc_max_nstep = static_cast <long int> p_sc_max_nstep;
        long int nproj = static_cast <long int> p_nproj;

        q_block.Allocate(static_cast <long int>(3 * blq.x), nq);
        setZero<2> << <bl, numThreads >> > (q_block, 3 * blq.x * nq);
        atype.Allocate(nproj);           

        if ((do_sc == 'C') || (do_sc == 'Y')) {
            q.Allocate(3, nproj, nq);           
            bl = (3 * nq + numThreads - 1) / numThreads;
            setZero<2> << <bl, numThreads >> > (q, 3 * nproj * nq);

        }
        if ((do_sc == 'Q') || (do_sc == 'Y')) {
            qt.Allocate(3, nproj, nq, sc_max_nstep);
            qw.Allocate(3, nproj, nq, nw);
            w_block.Allocate(static_cast <long int>(3 * blw.x), nq, nw);

            bl = (3 * nq* sc_max_nstep + numThreads - 1) / numThreads;
            setZero<3> << <bl, numThreads >> > (qt, 3 * nq * sc_max_nstep * nproj);
            bl = (3 * nq * nw + numThreads - 1) / numThreads;
            setZero<3> << <bl, numThreads >> > (qw, 3 * nq * nw * n_proj);

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
    */
};

