#pragma once

// Ivan Zivkovic, ivanzi@kth.se
// Requires C++ and CUDA 20 support

#include "gpuStructures.hpp"
#include "measurable.hpp"
#include "tensor.hpp"
#include "real_type.h"
#include "stopwatchDeviceSync.hpp"
#include "fortranData.hpp"
#include "measurementWriter.h"
#include "measurementData.h"
#include "kernels.hpp"
#include "gpu_wrappers.h"
#include "gpuParallelizationHelper.hpp"
#include "gpuCommon.hpp"
#include "measurementQueue.hpp"
#include "cpuRestMeasurement.hpp"
#include "autocorrelation_kernels.hpp"


class GpuMeasurement : public Measurable
{
public:
   GpuMeasurement(const deviceLattice& gpuLattice,
                    const deviceEnergies& gpuEnergies,
                    Tensor<real, 3>& f_emomM, 
                    Tensor<real, 3>& f_emom,
                    Tensor<real, 2>& f_mmom,
                    Tensor<real, 3>& f_beff,
                    MeasurementQueue& mq,
                    bool p_do_jtensor,
                    bool alwaysCopy = false
                    );
    ~GpuMeasurement() override;
    void measure(size_t mstep) override;
    void updateAC(size_t mstep) override;
    void flushMeasurements(size_t mstep) override;

private:
    bool timeToMeasure(MeasurementType mtype, size_t mstep) const;
    template<typename T> inline void fill_index(Vector<real>& iter, T step, size_t& count);
    void saveToFile(MeasurementType mtype);
    void calculateEmomMSum();
    void measureAverageMagnetization(size_t mstep);
    void measureBinderCumulant(size_t mstep);
    void measureSkyrmionNumber(size_t mstep);
    void measureEnergy(size_t mstep);
    static dim3 skyrmionKernelNumBlocks(SkyrmionMethod method, uint N, uint M, uint nsimp, uint kernel_threads);

    void measureAutocorrelation(size_t mstep);
    void release();
    MeasurementQueue& mqueue;
    CpuRestMeasurement cpuMeas;
    
private:
    bool isAllocated;
    const bool do_gpu_measurements;
    const deviceLattice& gpuLattice;
    const deviceEnergies& gpuEnergies;
    const uint N;
    const uint M;
    const uint NX, NY, NZ, NT;
    const uint NA, Natom_full, Nchmax;
    const char do_avrg_proj;
    const bool do_avrg_projch;
    const char do_cumu_proj;
    const bool do_cumu_projch;

    bool do_jtensor;
    const int do_ralloy;
    GPU_STREAM_T workStream;
    StopwatchDeviceSync stopwatch;
    MeasurementWriter measurementWriter;


    bool asitealloc;

    // Average magnetization
    const bool do_avrg;
    GpuVector<AverageMagnetizationData> mavg_buff_gpu;
    GpuTensor<AverageMagnetizationData, 2> mavg_proj_buff_gpu;
    GpuTensor<AverageMagnetizationData, 2> mavg_projch_buff_gpu;
    Vector<AverageMagnetizationData> mavg_buff_cpu;
    Tensor<AverageMagnetizationData, 2> mavg_proj_buff_cpu;
    Tensor<AverageMagnetizationData, 2> mavg_projch_buff_cpu;
    GpuVector<kernels::measurement::AvgMPart> mavg_partial_buff;
    GpuTensor<kernels::measurement::AvgMPart, 2> mavg_proj_partial_buff;
    GpuTensor<kernels::measurement::AvgMPart, 2> mavg_projch_partial_buff;
    Vector<real> mavg_iter;
    const dim3 mavg_kernel_threads;
    const dim3 mavg_kernel_blocks;
    dim3 mavg_proj_kernel_threads;
    dim3 mavg_proj_kernel_blocks;
    dim3 mavg_projch_kernel_threads;
    dim3 mavg_projch_kernel_blocks;
    size_t mavg_count = 0;
    GpuVector<int> achem_ch_gpu;
    GpuVector<int> asite_ch_gpu;
    GpuVector<int> atype_gpu;
    Vector<int> achem_ch_cpu;
    Vector<int> asite_ch_cpu;
    Vector<int> atype_cpu;


    // Energy
    const int do_ene;
    const int ene_types;
    dim3 ene_kernel_threads;
    dim3 ene_kernel_blocks;
    unsigned int ene_maxBlocks;
    unsigned int ene_maxThreads;
    const unsigned int ene_step;
    const unsigned int ene_buff;
    GpuVector<EnergyData> energy_buff_gpu;
    Vector<EnergyData> energy_buff_cpu;
    GpuTensor<kernels::measurement::EnePart, 2> energy_partial_buff;

    // Binder cumulant
    const bool do_cumu;
    GpuVector<BinderCumulantData> cumu_buff_gpu; // scalar but tensor of rank 0 is not allowed, so rank 1 is size 1
    Vector<BinderCumulantData> cumu_buff_cpu;
    GpuVector<kernels::measurement::BinderPart> cumu_partial_buff;
    //proj
    GpuVector<BinderCumulantData> cumu_proj_buff_gpu; // scalar but tensor of rank 0 is not allowed, so rank 1 is size 1
    Vector<BinderCumulantData> cumu_proj_buff_cpu;
    GpuTensor<kernels::measurement::BinderPart, 2> cumu_proj_partial_buff;
    //projch
    GpuVector<BinderCumulantData> cumu_projch_buff_gpu; // scalar but tensor of rank 0 is not allowed, so rank 1 is size 1
    Vector<BinderCumulantData> cumu_projch_buff_cpu;
    GpuTensor<kernels::measurement::BinderPart, 2> cumu_projch_partial_buff;
    const dim3 cumu_kernel_threads;
    const dim3 cumu_kernel_blocks;
    size_t cumu_count = 0;
    GpuVector<kernels::measurement::BinderEnePart> cumu_ene_partial_buff;
    const unsigned int cumu_ene_maxBlocks;
    const unsigned int cumu_ene_maxThreads;
 
    const dim3 cumu_ene_kernel_threads;
    const dim3 cumu_ene_kernel_blocks;

    dim3 cumu_proj_kernel_threads;
    dim3 cumu_proj_kernel_blocks;

    dim3 cumu_projch_kernel_threads;
    dim3 cumu_projch_kernel_blocks;


    // Used for both Average magnetization and Binder cumulant
    GpuTensor<real, 3> emomMEnsembleSums_partial;
    GpuTensor<real, 2> emomMEnsembleSums; // tensor of dim = 3 x M
    GpuTensor<real, 4> emomMEnsembleNTSums_partial;
    GpuTensor<real, 3> emomMEnsembleNTSums; 
    GpuTensor<real, 4> emomMEnsembleNASums_partial;
    GpuTensor<real, 3> emomMEnsembleNASums; 
    GpuTensor<real, 4> emomMEnsembleNCSums_partial;
    GpuTensor<real, 3> emomMEnsembleNCSums; 
    const dim3 sumOverAtoms_kernel_threads;
    const dim3 sumOverAtoms_kernel_blocks;

    dim3 sumOverAtoms_NT_kernel_threads;
    dim3 sumOverAtoms_NT_kernel_blocks;
    dim3 sumOverAtoms_NA_kernel_threads;
    dim3 sumOverAtoms_NA_kernel_blocks;
    dim3 sumOverAtoms_NC_kernel_threads;
    dim3 sumOverAtoms_NC_kernel_blocks;


    // Skyrmion number
    const SkyrmionMethod do_skyno;
    GpuVector<SkyrmionNumberData> skyno_buff_gpu;
    Vector<SkyrmionNumberData> skyno_buff_cpu;
    GpuVector<kernels::measurement::SumPart> skyno_partial_buff;
    Vector<real> skyno_iter;
    size_t skyno_count = 0;

    // Skyrmion number method brute force
    GpuTensor<real, 3> dxyz_vec;
    GpuTensor<int, 2> dxyz_atom;
    GpuVector<int> dxyz_list;
    GpuTensor<real, 4> grad_mom; // 3 x 3 x N x M, is initialized to zeros in Fortran, no need to copy over

    // Skyrmion number method triangulation
    GpuTensor<uint, 2> simp;
    const uint nsimp;
    const dim3 skyno_kernel_threads;
    const dim3 skyno_kernel_blocks;

    const real mub;
    const real mry;
    const real fcinv;
    //GpuVector<EnergyStdData> energy_std_buff_gpu;
    //Vector<EnergyData> energy_std_buff_cpu;

    Vector<real> energy_iter;
    size_t energy_count = 0;


    // Autocorrelations 
    const char do_autocorr;
   // GpuVector<unsigned int> spinwaittable_gpu;
    GpuTensor<real, 2> autocorr_buff_gpu;
    Tensor<real, 2> autocorr_buff_cpu;
    GpuTensor<real, 2> ac_block_gpu;
    GpuTensor<real, 4> spinwait_gpu;

    Vector<unsigned int> spinwaittable_cpu;
    //Tensor<real, 3> spinwait_cpu;

    Vector<real> indxb_ac;

    dim3 ac_threads;
    dim3 ac_blocks;
    unsigned int ac_threadsX;
    unsigned int ac_tasksX;
    unsigned int ac_blocksX;
    const unsigned int ac_maxThreads;
    const unsigned int ac_maxBlocks;

    size_t ac_count;
    unsigned int sw_next;
    int sw_curr;
    unsigned int sw_curIdx;
    const unsigned int ac_buff;
    const unsigned int ac_step;
    const unsigned int nspinwait;
    const unsigned int n0spinwait;
    unsigned int sw_threads;
    unsigned int sw_tasks;
    unsigned int sw_blocks;

};

