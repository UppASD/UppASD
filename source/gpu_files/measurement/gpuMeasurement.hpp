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


class GpuMeasurement : public Measurable
{
public:
    GpuMeasurement(const GpuTensor<real, 3>& emomM,
                    const GpuTensor<real, 3>& emom,
                    const GpuTensor<real, 2>& mmom,
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
    void release();
    
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
    GpuTensor<real, 2> emomMEnsembleSums; // tensor of dim = 3 x M

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

};

