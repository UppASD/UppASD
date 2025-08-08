#pragma once

#include "measurable.hpp"
#include "tensor.cuh"
#include "real_type.h"
#include "stopwatchDeviceSync.hpp"
#include "fortranData.hpp"
#include "measurementWriter.cuh"
#include "measurementData.h"


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
    cudaStream_t workStream;
    StopwatchDeviceSync stopwatch;
    MeasurementWriter measurementWriter;

    const bool do_avrg;
    const bool do_cumu;
    const bool do_skyno;

    // Average magnetization
    CudaVector<AverageMagnetizationData> mavg_buff_gpu;
    Vector<AverageMagnetizationData> mavg_buff_cpu;
    Vector<size_t> mavg_iter;
    uint mavg_count = 0;

    // Binder cumulant
    CudaVector<BinderCumulantData> cumu_buff_gpu; // scalar but tensor of rank 0 is not allowed, so rank 1 is size 1
    Vector<BinderCumulantData> cumu_buff_cpu;
    size_t cumu_count = 0;

    // Used for both Average magnetization and Binder cumulant
    CudaTensor<real, 2> emomMEnsembleSums; // tensor of dim = 3 x M

    // Skyrmion number
    CudaTensor<real, 3> dxyz_vec;
    CudaTensor<int, 2> dxyz_atom;
    CudaTensor<int, 1> dxyz_list;
    CudaTensor<real, 4> grad_mom; // 3 x 3 x N x M, is initialized to zeros in Fortran, no need to copy over
    CudaVector<SkyrmionNumberData> skyno_buff_gpu;
    Vector<SkyrmionNumberData> skyno_buff_cpu;
    Vector<size_t> skyno_iter;
    uint skyno_count = 0;
};


