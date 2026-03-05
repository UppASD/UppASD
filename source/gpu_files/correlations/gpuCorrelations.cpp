#pragma once

#include "c_headers.hpp"
#include "tensor.hpp"
#include "real_type.h"
#include "gpuStructures.hpp"
#include "fortranData.hpp"
#include <numeric>

#include "gpu_wrappers.h"

#if defined(HIP_V)
#include <correlation_kernels.hpp>
#include "gpuCorrelations.hpp"
#elif defined(CUDA_V)
#include <correlation_kernels.cuh>
#include "gpuCorrelations.cuh"
#endif


// Constructor
GpuCorrelations::GpuCorrelations(const Flag Flags, const SimulationParameters SimParam, const deviceLattice& gpuLattice, const hostCorrelations& cpuCorrelations)
: emomM(gpuLattice.emomM)
, emom(gpuLattice.emom)
, mmom(gpuLattice.mmom)
, do_sc(Flags.do_sc)
, do_proj(Flags.do_sc_proj)
, do_projch(Flags.do_sc_projch)
, sc_step(SimParam.sc_step)
, sc_sep(SimParam.sc_sep)
, N(SimParam.N)
, M(SimParam.M)
, nq(SimParam.nq)
, sc_max_nstep(SimParam.sc_max_nstep)
, sc_window_fun(SimParam.sc_window_fun)
, nw(SimParam.nw)
, NT(SimParam.NT)
, Nchmax(SimParam.Nchmax)
, delta_t(SimParam.delta_t)
, t_cur(0)
, n_samples(0)
, maxThreads(256)
, maxBlocks(1024)//should be equal to current threads per block limit of GPU
, numThreads(maxThreads)
, blQ(N, M, nq, numThreads, maxBlocks)
, blW(N, M, nq, sc_max_nstep, nw, numThreads, maxBlocks)
, blQproj(N, M, nq, NT, numThreads, maxBlocks)
, blWproj(N, M, nq, sc_max_nstep, nw, NT, numThreads, maxBlocks)
, blQprojch(N, M, nq, Nchmax, numThreads, maxBlocks)
, blWprojch(N, M, nq, sc_max_nstep, nw, Nchmax, numThreads, maxBlocks)
, sc(do_sc, nw, nq, sc_max_nstep, blQ, blW)
, sc_proj(do_sc_proj, nw, nq, sc_max_nstep, NT, blQproj, blWproj)
, sc_projch(do_sc_projch, nw, nq, sc_max_nstep, Nchmax, blQprojch, blWprojch)

{
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

        threads = { numThreads, 1, 1 };

        r_mid.Allocate(static_cast <long int>(3));
        q.Allocate(static_cast <long int>(3), static_cast <long int>(nq));
        coord.Allocate(static_cast <long int>(3), static_cast <long int>(N));

        r_mid.copy_sync(cpuCorrelations.r_mid);
        q.copy_sync(cpuCorrelations.q);
        coord.copy_sync(cpuCorrelations.coord);


        if ((do_sc == 'Q') || (do_sc == 'Y')) {
            dt.Allocate(static_cast <long int>(sc_max_nstep));
            dt_cpu.AllocateHost(static_cast <long int>(sc_max_nstep));
            sc_step_arr_cpu.AllocateHost(static_cast <long int>(sc_max_nstep));
            w.Allocate(static_cast <long int>(nw));
            w.copy_sync(cpuCorrelations.w);
        }

        isallocated = 1; 
    }

    // All initialized?
    if (GPU_DEVICE_SYNCHRONIZE() != GPU_SUCCESS) {
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
        if ((do_sc == 'Q') || (do_sc == 'Y')) {
            w.Free();
            dt.Free();
            dt_cpu.FreeHost();
            sc_step_arr_cpu.FreeHost();
        }
        sc.free(do_sc);
        sc_proj.free(do_sc_proj);
        sc_projch.free(do_sc_projch);

        isallocated = 0;
    }

}

void GpuCorrelations::measure(std::size_r mstep){
    measure_SC(std::size_t mstep);
    

}

void GpuCorrelations::measure_SC(std::size_t mstep) {
    
    std::size_t curstep = mstep;
    switch (do_sc) {
    case 'C':
        if ((curstep % sc_sep) == 0) {
            GPUSqSum << <blQ.blocks, threads >> > (emomM, coord, q, r_mid, sc.q_block, blQ.tasks, N);
            GPUSqFinalSum_stat << <nq, maxBlocks>> > (sc.q_block, sc.q, blQ.x);
            GPU_DEVICE_SYNCHRONIZE();
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
            GPUSqSum << <blQ.blocks, threads >> > (emomM, coord, q, r_mid, sc.q_block, blQ.tasks, N);
            GPUSqFinalSum_dyn << <nq, maxBlocks>> > (sc.q_block, sc.qt, blQ.x, t_cur);
            GPU_DEVICE_SYNCHRONIZE();
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
            
            GPUSqSum << <blQ.blocks, threads >> > (emomM, coord, q, r_mid, sc.q_block, blQ.tasks, N);
            GPUSqFinalSum_both << <nq, maxBlocks >> > (sc.q_block, sc.q, sc.qt, blQ.x, t_cur, both_flag);
            GPU_DEVICE_SYNCHRONIZE();
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
            
            GPUSqSum << <blQ.block, threads >> > (emomM, coord, q, r_mid, sc.q_block, blq.tasks, N);
            GPUSqFinalSum_both << <nq, maxBlocks >> > (sc.q_block, sc.q, sc.qt, blQ.x, t_cur, both_flag);
            GPU_DEVICE_SYNCHRONIZE();
            t_cur++;
        }
        else if ((curstep % sc_sep) == 0) {
            both_flag = 0;

            GPUSqSum << <blQ.blocks, threads >> > (emomM, coord, q, r_mid, sc_block_gpu, tasksTot_q, N);
            GPUSqFinalSum_both << <nq, maxBlocks >> > (sc.q_block, sc.q, sc.qt, blQ.x, t_cur, both_flag);
            GPU_DEVICE_SYNCHRONIZE();
            n_samples++;
        }
        break;

    }    

}

void GpuCorrelations::measure_SC_proj(std::size_t mstep) {
    
   /* std::size_t curstep = mstep;
    switch (do_sc) {
    case 'C':
        if ((curstep % sc_sep) == 0) {
            GPUSqSum << <blQ.blocks, threads >> > (emomM, coord, q, r_mid, sc.q_block, blQ.tasks, N);
            GPUSqFinalSum_stat << <nq, maxBlocks>> > (sc.q_block, sc.q, blQ.x);
            GPU_DEVICE_SYNCHRONIZE();
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
            GPUSqSum << <blQ.blocks, threads >> > (emomM, coord, q, r_mid, sc.q_block, blQ.tasks, N);
            GPUSqFinalSum_dyn << <nq, maxBlocks>> > (sc.q_block, sc.qt, blQ.x, t_cur);
            GPU_DEVICE_SYNCHRONIZE();
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
            
            GPUSqSum << <blQ.blocks, threads >> > (emomM, coord, q, r_mid, sc.q_block, blQ.tasks, N);
            GPUSqFinalSum_both << <nq, maxBlocks >> > (sc.q_block, sc.q, sc.qt, blQ.x, t_cur, both_flag);
            GPU_DEVICE_SYNCHRONIZE();
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
            
            GPUSqSum << <blQ.block, threads >> > (emomM, coord, q, r_mid, sc.q_block, blq.tasks, N);
            GPUSqFinalSum_both << <nq, maxBlocks >> > (sc.q_block, sc.q, sc.qt, blQ.x, t_cur, both_flag);
            GPU_DEVICE_SYNCHRONIZE();
            t_cur++;
        }
        else if ((curstep % sc_sep) == 0) {
            both_flag = 0;

            GPUSqSum << <blQ.blocks, threads >> > (emomM, coord, q, r_mid, sc_block_gpu, tasksTot_q, N);
            GPUSqFinalSum_both << <nq, maxBlocks >> > (sc.q_block, sc.q, sc.qt, blQ.x, t_cur, both_flag);
            GPU_DEVICE_SYNCHRONIZE();
            n_samples++;
        }
        break;

    }  */  

}


void GpuCorrelations::flushCorrelations(hostCorrelations& cpuCorrelations, std::size_t mstep) {
    int tasks; int bl;
    switch (do_sc) {
    case 'C': {
        cpuCorrelations.m_k.copy_sync(sc.q);
        break;
    }

    case 'Q': {
        // Copy time step data to GPU
        dt.copy_sync(dt_cpu);
        
        // Compute partial S(q,ω) from S(q,t) using Fourier transform
        GPUSwSum << <blW.blocks, threads >> > (sc.qt, dt, w, sc.w_block, blW.tasks, sc_max_nstep, nq, sc_max_nstep, sc_window_fun);
        GPUSwFinalSum << <nq * nw, maxBlocks >> > (sc.w_block, sc.qw, blW.x, nq);
        
        // Transfer time-domain correlations for Fortran reference
        if (sc.qt.extent(0) == cpuCorrelations.m_kt.extent(0) &&
            sc.qt.extent(1) == cpuCorrelations.m_kt.extent(1) &&
            sc.qt.extent(2) == cpuCorrelations.m_kt.extent(2)) {
            cpuCorrelations.m_kt.copy_sync(sc.qt);
            cpuCorrelations.sc_tidx = static_cast<int>(sc.qt.extent(2));
        }
        
        // Transfer frequency-domain correlations (m_kw)
        if (sc.qw.extent(0) == cpuCorrelations.m_kw.extent(0) &&
            sc.qw.extent(1) == cpuCorrelations.m_kw.extent(1) &&
            sc.qw.extent(2) == cpuCorrelations.m_kw.extent(2)) {
            cpuCorrelations.m_kw.copy_sync(sc.qw);
        }
        break;
    }

    case 'Y': {
        // Copy time step data to GPU
        dt.copy_sync(dt_cpu);
        
        // Compute partial S(q,ω) from S(q,t) using Fourier transform
        GPUSwSum << <blW.blocks, threads >> > (sc.qt, dt, w, sc.w_block, blW.tasks, sc_max_nstep, nq, sc_max_nstep, sc_window_fun);
        GPUSwFinalSum << <nq * nw, maxBlocks >> > (sc.w_block, sc.qw, blW.x, nq);
        
        // Transfer static S(q)
        cpuCorrelations.m_k.copy_sync(sc.q);
        
        // Transfer time-domain S(q,t)
        if (sc.qt.extent(0) == cpuCorrelations.m_kt.extent(0) &&
            sc.qt.extent(1) == cpuCorrelations.m_kt.extent(1) &&
            sc.qt.extent(2) == cpuCorrelations.m_kt.extent(2)) {
            cpuCorrelations.m_kt.copy_sync(sc.qt);
            cpuCorrelations.sc_tidx = static_cast<int>(sc.qt.extent(2));
        }
        
        // Transfer frequency-domain S(q,ω)
        if (sc.qw.extent(0) == cpuCorrelations.m_kw.extent(0) &&
            sc.qw.extent(1) == cpuCorrelations.m_kw.extent(1) &&
            sc.qw.extent(2) == cpuCorrelations.m_kw.extent(2)) {
            cpuCorrelations.m_kw.copy_sync(sc.qw);
        }
        break;
    }
    
    default: {
        // Case C: S(q) static correlations only
        cpuCorrelations.m_k.copy_sync(sc.q);
        break;
    }
    }  // end switch
    
    GPU_DEVICE_SYNCHRONIZE();
    
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






