#pragma once

#include "c_headers.hpp"
#include "tensor.hpp"
#include "real_type.h"
#include "gpuStructures.hpp"
#include "gpuCorrelations.cuh"
#include "fortranData.hpp"
#include <numeric>
#include <curand.h>
#include <cuda.h>
#include <correlation_kernels.cuh>


// Constructor
GpuCorrelations::GpuCorrelations(const Flag Flags, const SimulationParameters SimParam, const deviceLattice& gpuLattice, const hostCorrelations& cpuCorrelations)
: emomM(gpuLattice.emomM)
, emom(gpuLattice.emom)
, mmom(gpuLattice.mmom) {

    isallocated = 0; 
    if(!initiate(Flags, SimParam, cpuCorrelations)) {  
      std::fprintf(stderr, "GpuCorrelations: correlations failed to initiate!\n");
      return;
   }
}
// Destructor
GpuCorrelations::~GpuCorrelations() {
    release();
}
// Initiator
bool GpuCorrelations::initiate(const Flag Flags, const SimulationParameters SimParam, const hostCorrelations& cpuCorrelations) {
    // Assert that we're not already initialized
    //release();

    // Parameters
    if(Flags.do_gpu_correlations){

        N = SimParam.N;
        M = SimParam.M;
        nq = SimParam.nq;
        sc_max_nstep = SimParam.sc_max_nstep;
        sc_window_fun = SimParam.sc_window_fun;
        nw = SimParam.nw;
        delta_t = SimParam.delta_t;
        t_cur = 0;
        n_samples = 0;
        do_sc = Flags.do_sc;
        sc_sep = SimParam.sc_sep;
        sc_step = SimParam.sc_step;
        // nainv = 1 / N;
        // Blocks and threads
        maxThreads = 512;
        maxBlocks = 1024; 
        tasksTot_q = 3 * N * M;
        tasksTot_w = 3 * sc_max_nstep;
        // maxBlocks = 1023; //must be devidable by 3, less than 1024
        numThreads = maxThreads;
        //numBlocks = std::min(((3 * ((spinTot + 2) / 3) + numThreads - 1) / numThreads), maxBlocks);
        numBlocksX_q = std::min(((tasksTot_q + numThreads - 1) / numThreads), maxBlocks);
        numBlocksX_w = std::min(((tasksTot_w + numThreads - 1) / numThreads), maxBlocks);
        numBlocksY_q = nq;
        numBlocksY_w = nq*nw;//TODO
        blocks_q = { numBlocksX_q, numBlocksY_q, 1 };
        blocks_w = { numBlocksX_w, numBlocksY_w, 1 };//TODO
        threads = { numThreads, 1, 1 };
        //printf("numBlocks = %i\n", numBlocksX_q);

        //iqfac = thrust::complex<real>(0, 2 * M_PI);
        r_mid.Allocate(static_cast <long int>(3));
        q.Allocate(static_cast <long int>(3), static_cast <long int>(nq));
        coord.Allocate(static_cast <long int>(3), static_cast <long int>(N));

        r_mid.copy_sync(cpuCorrelations.r_mid);
        q.copy_sync(cpuCorrelations.q);
        coord.copy_sync(cpuCorrelations.coord);
        int bl;

        sc_block_gpu.Allocate(static_cast <long int>(3 * numBlocksX_q), static_cast <long int>(nq));
        if ((do_sc == 'C') || (do_sc == 'Y')) {
            sc_q_gpu.Allocate(static_cast <long int>(3), static_cast <long int>(nq));
            sc_q_cpu.AllocateHost(static_cast <long int>(3), static_cast <long int>(nq));
            bl = (3 * nq + numThreads - 1) / numThreads;
            setZero<2> << <bl, numThreads >> > (sc_q_gpu, 3 * nq);

        }
        if ((do_sc == 'Q') || (do_sc == 'Y')) {
            // CRITICAL: Match Fortran memory layout: (component, q, time) NOT (component, time, q)
            sc_qt_gpu.Allocate(static_cast <long int>(3), static_cast <long int>(nq), static_cast <long int>(sc_max_nstep));
            sc_qw_gpu.Allocate(static_cast <long int>(3), static_cast <long int>(nq), static_cast <long int>(nw));
            sc_block_w_gpu.Allocate(static_cast <long int>(3 * numBlocksX_w), static_cast <long int>(nq), static_cast <long int>(nw));
            dt.Allocate(static_cast <long int>(sc_max_nstep));
            dt_cpu.AllocateHost(static_cast <long int>(sc_max_nstep));
            sc_step_arr_cpu.AllocateHost(static_cast <long int>(sc_max_nstep));
            w.Allocate(static_cast <long int>(nw));
            //dt.copy_sync(cpuCorrelations.dt);const deviceLattice& gpuLattice, const int curstep
            w.copy_sync(cpuCorrelations.w);


            bl = (3 * nq* sc_max_nstep + numThreads - 1) / numThreads;
            setZero<3> << <bl, numThreads >> > (sc_qt_gpu, 3 * nq* sc_max_nstep);
            bl = (3 * nq * nw + numThreads - 1) / numThreads;
            setZero<3> << <bl, numThreads >> > (sc_qw_gpu, 3 * nq*nw);
            //bl = (3 * nq * sc_max_nstep + numThreads - 1) / numThreads;
        }

        //mbuff_gpu.Allocate(static_cast <long int>(3), static_cast <long int>(avrg_buff), static_cast <long int>(M));
        isallocated = 1;
        bl = (3 * numBlocksX_q * nq + numThreads - 1) / numThreads;
        setZero<2> << <bl, numThreads >> > (sc_block_gpu, 3 * numBlocksX_q * nq);
        //sc_block_gpu.zeros();
        //sc_gpu.zeros();
    }

    // All initialized?
    if (cudaDeviceSynchronize() != cudaSuccess) {
        release();
        return false;
    }

    return true;
}
void GpuCorrelations::release() {
    if (isallocated) {
        r_mid.Free();
        coord.Free();
        q.Free();
        sc_block_gpu.Free();
        if ((do_sc == 'C') || (do_sc == 'Y')) {
            sc_q_gpu.Free();
            sc_q_cpu.FreeHost();
        }
        if ((do_sc == 'Q') || (do_sc == 'Y')) {
            sc_qt_gpu.Free();
            sc_qw_gpu.Free();
            sc_block_w_gpu.Free();
            w.Free();
            dt.Free();
            dt_cpu.FreeHost();
            sc_step_arr_cpu.FreeHost();
        }
        isallocated = 0;
    }

}

void GpuCorrelations::measure(std::size_t mstep) {
    
    std::size_t curstep = mstep;
    switch (do_sc) {
    case 'C':
        if ((curstep % sc_sep) == 0) {
            GPUSqSum << <blocks_q, threads >> > (emomM, coord, q, r_mid, sc_block_gpu, tasksTot_q, N);
            GPUSqFinalSum_stat << <nq, 1024 >> > (sc_block_gpu, sc_q_gpu, numBlocksX_q);
            cudaDeviceSynchronize();
            n_samples++;
        }
        break;

    case 'Q':
        if ((curstep % sc_step) == 0 && t_cur < static_cast<unsigned int>(sc_max_nstep)) {
            // printf("[GPU-SAMPLE Q] mstep=%zu, curstep=%zu, condition: (curstep %% sc_step)==0, t_cur=%u (max=%u)\n", 
            //        mstep, curstep, t_cur, sc_max_nstep);
            
            // Record metadata BEFORE kernel call (mirrors Fortran: increment then record)
            if (t_cur < static_cast<unsigned int>(dt_cpu.extent(0))) {
                dt_cpu(int(t_cur)) = delta_t * sc_step;  // Record time step at current index

            }
            if (t_cur < static_cast<unsigned int>(sc_step_arr_cpu.extent(0))) {
                sc_step_arr_cpu(int(t_cur)) = static_cast<real>(sc_step);  // Record step width

            }
            
            // Kernel writes to m_kt(:,:,t_cur)
            GPUSqSum << <blocks_q, threads >> > (emomM, coord, q, r_mid, sc_block_gpu, tasksTot_q, N);
            GPUSqFinalSum_dyn << <nq, 1024 >> > (sc_block_gpu, sc_qt_gpu, numBlocksX_q, t_cur);
            cudaDeviceSynchronize();
            t_cur++;  // Increment AFTER writing to that time slice
        } else {
            if ((curstep % sc_step) == 0) {
            }
        }
        break;

    case 'Y':
        if (((curstep % sc_step) == 0) && ((curstep % sc_sep) == 0) && t_cur < static_cast<unsigned int>(sc_max_nstep)) {
            both_flag = 2;

            // Record metadata for time-domain sample
            if (t_cur < static_cast<unsigned int>(dt_cpu.extent(0))) {
                dt_cpu(int(t_cur)) = delta_t * sc_step;
            }
            if (t_cur < static_cast<unsigned int>(sc_step_arr_cpu.extent(0))) {
                sc_step_arr_cpu(int(t_cur)) = static_cast<real>(sc_step);
            }
            
            GPUSqSum << <blocks_q, threads >> > (emomM, coord, q, r_mid, sc_block_gpu, tasksTot_q, N);
            GPUSqFinalSum_both << <nq, 1024 >> > (sc_block_gpu, sc_q_gpu, sc_qt_gpu, numBlocksX_q, t_cur, both_flag);
            cudaDeviceSynchronize();
            t_cur++;
            n_samples++;

        }
        else if ((curstep % sc_step) == 0 && t_cur < static_cast<unsigned int>(sc_max_nstep)) {
            both_flag = 1;

            // Record metadata for time-domain sample
            if (t_cur < static_cast<unsigned int>(dt_cpu.extent(0))) {
                dt_cpu(int(t_cur)) = delta_t * sc_step;
            }
            if (t_cur < static_cast<unsigned int>(sc_step_arr_cpu.extent(0))) {
                sc_step_arr_cpu(int(t_cur)) = static_cast<real>(sc_step);
            }
            
            GPUSqSum << <blocks_q, threads >> > (emomM, coord, q, r_mid, sc_block_gpu, tasksTot_q, N);
            GPUSqFinalSum_both << <nq, 1024 >> > (sc_block_gpu, sc_q_gpu, sc_qt_gpu, numBlocksX_q, t_cur, both_flag);
            cudaDeviceSynchronize();
            t_cur++;
        }
        else if ((curstep % sc_sep) == 0) {
            both_flag = 0;

            GPUSqSum << <blocks_q, threads >> > (emomM, coord, q, r_mid, sc_block_gpu, tasksTot_q, N);
            GPUSqFinalSum_both << <nq, 1024 >> > (sc_block_gpu, sc_q_gpu, sc_qt_gpu, numBlocksX_q, t_cur, both_flag);
            cudaDeviceSynchronize();
            n_samples++;
        }
        break;

    }    

}

void GpuCorrelations::flushCorrelations(hostCorrelations& cpuCorrelations, std::size_t mstep) {
    int tasks; int bl;
    switch (do_sc) {
    case 'C': {
        //tasks = 3 * nq;
        //bl = (tasks + maxThreads - 1) / maxThreads;
        //GPUSqAvrg << <bl, maxThreads >> > (sc_q_gpu, n_samples, tasks, M);  
        //cudaDeviceSynchronize();
        
        // Transfer to host CPU tensor first
        //sc_q_cpu.copy_sync(sc_q_gpu);
        
        // Transfer to cpuCorrelations for Fortran
        cpuCorrelations.m_k.copy_sync(sc_q_gpu);
        break;
    }

    case 'Q': {
        // Copy time step data to GPU
        dt.copy_sync(dt_cpu);
        
        // Zero intermediate and output buffers
        int bl_w = (3 * numBlocksX_w * nq * nw + numThreads - 1) / numThreads;
        setZero<3> << <bl_w, numThreads >> > (sc_block_w_gpu, 3 * numBlocksX_w * nq * nw);
        bl_w = (3 * nq * nw + numThreads - 1) / numThreads;
        setZero<3> << <bl_w, numThreads >> > (sc_qw_gpu, 3 * nq * nw);
        cudaDeviceSynchronize();
        
        // Compute partial S(q,ω) from S(q,t) using Fourier transform
        GPUSwSum << <blocks_w, threads >> > (sc_qt_gpu, dt, w, sc_block_w_gpu, tasksTot_w, sc_max_nstep, nq, sc_max_nstep, sc_window_fun);
        cudaDeviceSynchronize();
        
        // Reduce block results to final S(q,ω)
        GPUSwFinalSum << <nq * nw, 1024 >> > (sc_block_w_gpu, sc_qw_gpu, numBlocksX_w, nq);
        cudaDeviceSynchronize();
        
        // Transfer time-domain correlations for Fortran reference
        if (sc_qt_gpu.extent(0) == cpuCorrelations.m_kt.extent(0) &&
            sc_qt_gpu.extent(1) == cpuCorrelations.m_kt.extent(1) &&
            sc_qt_gpu.extent(2) == cpuCorrelations.m_kt.extent(2)) {
            cpuCorrelations.m_kt.copy_sync(sc_qt_gpu);
            cpuCorrelations.sc_tidx = static_cast<int>(sc_qt_gpu.extent(2));
        }
        
        // Transfer frequency-domain correlations (m_kw)
        if (sc_qw_gpu.extent(0) == cpuCorrelations.m_kw.extent(0) &&
            sc_qw_gpu.extent(1) == cpuCorrelations.m_kw.extent(1) &&
            sc_qw_gpu.extent(2) == cpuCorrelations.m_kw.extent(2)) {
            cpuCorrelations.m_kw.copy_sync(sc_qw_gpu);
        }
        break;
    }

    case 'Y': {
        // Copy time step data to GPU
        dt.copy_sync(dt_cpu);
        
        // Zero intermediate and output FFT buffers
        int bl_w = (3 * numBlocksX_w * nq * nw + numThreads - 1) / numThreads;
        setZero<3> << <bl_w, numThreads >> > (sc_block_w_gpu, 3 * numBlocksX_w * nq * nw);
        bl_w = (3 * nq * nw + numThreads - 1) / numThreads;
        setZero<3> << <bl_w, numThreads >> > (sc_qw_gpu, 3 * nq * nw);
        cudaDeviceSynchronize();
        
        // Average static S(q) correlations
        //tasks = 3 * nq;
        //bl = (tasks + maxThreads - 1) / maxThreads;
        //GPUSqAvrg << <bl, maxThreads >> > (sc_q_gpu, n_samples, tasks, M);
        //cudaDeviceSynchronize();
        
        // Compute partial S(q,ω) from S(q,t) using Fourier transform
        GPUSwSum << <blocks_w, threads >> > (sc_qt_gpu, dt, w, sc_block_w_gpu, tasksTot_w, sc_max_nstep, nq, sc_max_nstep, sc_window_fun);
        cudaDeviceSynchronize();
        
        // Reduce block results to final S(q,ω)
        GPUSwFinalSum << <nq * nw, 1024 >> > (sc_block_w_gpu, sc_qw_gpu, numBlocksX_w, nq);
        cudaDeviceSynchronize();
        
        // Transfer static S(q)
        cpuCorrelations.m_k.copy_sync(sc_q_gpu);
        
        // Transfer time-domain S(q,t)
        if (sc_qt_gpu.extent(0) == cpuCorrelations.m_kt.extent(0) &&
            sc_qt_gpu.extent(1) == cpuCorrelations.m_kt.extent(1) &&
            sc_qt_gpu.extent(2) == cpuCorrelations.m_kt.extent(2)) {
            cpuCorrelations.m_kt.copy_sync(sc_qt_gpu);
            cpuCorrelations.sc_tidx = static_cast<int>(sc_qt_gpu.extent(2));
        }
        
        // Transfer frequency-domain S(q,ω)
        if (sc_qw_gpu.extent(0) == cpuCorrelations.m_kw.extent(0) &&
            sc_qw_gpu.extent(1) == cpuCorrelations.m_kw.extent(1) &&
            sc_qw_gpu.extent(2) == cpuCorrelations.m_kw.extent(2)) {
            cpuCorrelations.m_kw.copy_sync(sc_qw_gpu);
        }
        break;
    }
    
    default: {
        // Case C: S(q) static correlations only
        //tasks = 3 * nq;
        //bl = (tasks + maxThreads - 1) / maxThreads;
        //GPUSqAvrg << <bl, maxThreads >> > (sc_q_gpu, n_samples, tasks, M);
        //cudaDeviceSynchronize();
        
        // Transfer to cpuCorrelations for Fortran
        cpuCorrelations.m_k.copy_sync(sc_q_gpu);
        break;
    }
    }  // end switch
    
    cudaDeviceSynchronize();
    
    // Publish sampling info to Fortran (CRITICAL: updates sc_nsamp and sc_tidx)
    publishSamplingInfo(cpuCorrelations);

}

void GpuCorrelations::recordSample() {
    // Record current time and step size when a sample is taken
    if (n_samples < static_cast<std::size_t>(dt_cpu.extent(0))) {
        dt_cpu(int(n_samples)) = delta_t;
    }
    // Track which step this sample corresponds to
    if (t_cur < static_cast<unsigned int>(sc_step_arr_cpu.extent(0))) {
        sc_step_arr_cpu(int(t_cur)) = static_cast<real>(sc_step);
    }
    t_cur++;
    n_samples++;
}

void GpuCorrelations::publishSamplingInfo(hostCorrelations& cpuCorrelations) {
    // Copy the recorded delta_t values back to Fortran's deltat_corr array
    if (FortranData::deltat_corr != nullptr && dt_cpu.extent(0) > 0) {
        std::size_t n_copy = std::min(static_cast<std::size_t>(dt_cpu.extent(0)), n_samples);
        for (std::size_t i = 0; i < n_copy; i++) {
            FortranData::deltat_corr[i] = dt_cpu(int(i));
        }
    }
    
    // Copy the recorded sc_step values back to Fortran's scstep_arr array
    if (FortranData::scstep_arr != nullptr && sc_step_arr_cpu.extent(0) > 0) {
        std::size_t n_copy = std::min(static_cast<std::size_t>(sc_step_arr_cpu.extent(0)), n_samples);
        for (std::size_t i = 0; i < n_copy; i++) {
            FortranData::scstep_arr[i] = sc_step_arr_cpu(int(i));
        }
    }
    
    // Update cpuCorrelations with the sample count and time index
    // GPU does NOT normalize - Fortran will handle normalization using sc_nsamp
    cpuCorrelations.sc_nsamp = static_cast<int>(n_samples);
    cpuCorrelations.sc_tidx = static_cast<int>(t_cur);
    
    // Write sample count and time index back to Fortran through FortranData pointers
    if (FortranData::sc_nsamp_ptr != nullptr) {
        *FortranData::sc_nsamp_ptr = static_cast<int>(n_samples);
    }
    if (FortranData::sc_tidx_ptr != nullptr) {
        *FortranData::sc_tidx_ptr = static_cast<int>(t_cur);
    }
}






