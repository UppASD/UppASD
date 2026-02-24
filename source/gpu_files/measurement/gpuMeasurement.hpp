#pragma once

// Ivan Zivkovic, ivanzi@kth.se
// Requires C++ and CUDA 20 support

#include "measurable.hpp"
#include "tensor.hpp"
#include "real_type.h"
#include "stopwatchDeviceSync.hpp"
#include "fortranData.hpp"
#include "measurementWriter.h"
#include "measurementData.h"
#include "kernels.hpp"
#include "gpu_wrappers.h"
#include "measurementQueue.hpp"
#include "cpuRestMeasurement.hpp"
#include "autocorrelation_kernels.hpp"


class GpuMeasurement : public Measurable
{
public:
    GpuMeasurement(const GpuTensor<real, 3>& emomM,
                    const GpuTensor<real, 3>& emom,
                    const GpuTensor<real, 2>& mmom,
                    const GpuTensor<real, 3>& beff,
                    Tensor<real, 3>& f_emomM, 
                    Tensor<real, 3>& f_emom,
                    Tensor<real, 2>& f_mmom,
                    Tensor<real, 3>& f_beff,
                    MeasurementQueue& mq,
                    bool alwaysCopy = false);
    ~GpuMeasurement() override;
    void measure(size_t mstep) override;
    void flushMeasurements(size_t mstep) override;


private:
    bool isAllocated;
    bool timeToMeasure(MeasurementType mtype, size_t mstep) const;
    void saveToFile(MeasurementType mtype);
    void calculateEmomMSum();
    void measureAverageMagnetization(size_t mstep);
    void measureBinderCumulant(size_t mstep);
    void measureSkyrmionNumber(size_t mstep);
    void measureAutocorrelation(size_t mstep);
    void release();
    MeasurementQueue& mqueue;
    CpuRestMeasurement cpuMeas;
    
private:
    const GpuTensor<real, 3>& emomM;
    const GpuTensor<real, 3>& emom;
    const GpuTensor<real, 2>& mmom;
    const uint N;
    const uint M;
    const uint NX, NY, NZ, NT;
    GPU_STREAM_T workStream;
    StopwatchDeviceSync stopwatch;
    MeasurementWriter measurementWriter;


    // Average magnetization
    const bool do_avrg;
    GpuVector<AverageMagnetizationData> mavg_buff_gpu;
    Vector<AverageMagnetizationData> mavg_buff_cpu;
    GpuVector<kernels::measurement::AvgMPart> mavg_partial_buff;
    Vector<size_t> mavg_iter;
    const dim3 mavg_kernel_threads;
    const dim3 mavg_kernel_blocks;
    size_t mavg_count = 0;

    // Binder cumulant
    const bool do_cumu;
    GpuVector<BinderCumulantData> cumu_buff_gpu; // scalar but tensor of rank 0 is not allowed, so rank 1 is size 1
    Vector<BinderCumulantData> cumu_buff_cpu;
    GpuVector<kernels::measurement::BinderPart> cumu_partial_buff;
    const dim3 cumu_kernel_threads;
    const dim3 cumu_kernel_blocks;
    size_t cumu_count = 0;

    // Used for both Average magnetization and Binder cumulant
    GpuTensor<real, 3> emomMEnsembleSums_partial;
    GpuTensor<real, 2> emomMEnsembleSums; // tensor of dim = 3 x M
    const dim3 sumOverAtoms_kernel_threads;
    const dim3 sumOverAtoms_kernel_blocks;

    // Skyrmion number
    const SkyrmionMethod do_skyno;
    GpuVector<SkyrmionNumberData> skyno_buff_gpu;
    Vector<SkyrmionNumberData> skyno_buff_cpu;
    GpuVector<kernels::measurement::SumPart> skyno_partial_buff;
    Vector<size_t> skyno_iter;
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

    // Autocorrelations 
    const bool do_autocorr;
   // GpuVector<unsigned int> spinwaittable_gpu;
    GpuTensor<real, 2> autocorr_buff;
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

    unsigned int ac_count;
    unsigned int sw_next;
    unsigned int sw_curIdx;
    const unsigned int ac_buff;
    const unsigned int ac_step;
    const unsigned int nspinwait;
    const unsigned int n0spinwait;
    unsigned int sw_threads;
    unsigned int sw_tasks;
    unsigned int sw_blocks;

};

