#pragma once

// Ivan Zivkovic, ivanzi@kth.se
// Requires C++ and CUDA 20 support

#include "measurable.hpp"
#include "tensor.cuh"
#include "real_type.h"
#include "stopwatchDeviceSync.hpp"
#include "fortranData.hpp"
#include "measurementWriter.h"
#include "measurementData.h"
#include "kernels.cuh"


class CudaMeasurement : public Measurable
{
public:
    CudaMeasurement(const CudaTensor<real, 3>& emomM,
                    const CudaTensor<real, 3>& emom,
                    const CudaTensor<real, 2>& mmom,
                    bool alwaysCopy = false);
    ~CudaMeasurement() override;
    void measure(size_t mstep) override;
    void flushMeasurements(size_t mstep) override;

private:
    bool timeToMeasure(MeasurementType mtype, size_t mstep) const;
    void saveToFile(MeasurementType mtype);
    void calculateEmomMSum();
    void measureAverageMagnetization(size_t mstep);
    void measureBinderCumulant(size_t mstep);
    void measureSkyrmionNumber(size_t mstep);
    
private:
    const CudaTensor<real, 3>& emomM;
    const CudaTensor<real, 3>& emom;
    const CudaTensor<real, 2>& mmom;
    const uint N;
    const uint M;
    const uint NX, NY, NZ, NT;
    cudaStream_t workStream;
    StopwatchDeviceSync stopwatch;
    MeasurementWriter measurementWriter;


    // Average magnetization
    const bool do_avrg;
    CudaVector<AverageMagnetizationData> mavg_buff_gpu;
    Vector<AverageMagnetizationData> mavg_buff_cpu;
    CudaVector<kernels::measurement::AvgMPart> mavg_partial_buff;
    Vector<size_t> mavg_iter;
    const dim3 mavg_kernel_threads;
    const dim3 mavg_kernel_blocks;
    size_t mavg_count = 0;

    // Binder cumulant
    const bool do_cumu;
    CudaVector<BinderCumulantData> cumu_buff_gpu; // scalar but tensor of rank 0 is not allowed, so rank 1 is size 1
    Vector<BinderCumulantData> cumu_buff_cpu;
    CudaVector<kernels::measurement::BinderPart> cumu_partial_buff;
    const dim3 cumu_kernel_threads;
    const dim3 cumu_kernel_blocks;
    size_t cumu_count = 0;

    // Used for both Average magnetization and Binder cumulant
    CudaTensor<real, 2> emomMEnsembleSums; // tensor of dim = 3 x M

    // Skyrmion number
    const SkyrmionMethod do_skyno;
    CudaVector<SkyrmionNumberData> skyno_buff_gpu;
    Vector<SkyrmionNumberData> skyno_buff_cpu;
    CudaVector<kernels::measurement::SumPart> skyno_partial_buff;
    Vector<size_t> skyno_iter;
    size_t skyno_count = 0;

    // Skyrmion number method brute force
    CudaTensor<real, 3> dxyz_vec;
    CudaTensor<int, 2> dxyz_atom;
    CudaVector<int> dxyz_list;
    CudaTensor<real, 4> grad_mom; // 3 x 3 x N x M, is initialized to zeros in Fortran, no need to copy over

    // Skyrmion number method triangulation
    CudaTensor<uint, 2> simp;
    const uint nsimp;
    const dim3 skyno_kernel_threads;
    const dim3 skyno_kernel_blocks;

};

