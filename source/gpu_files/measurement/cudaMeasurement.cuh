#pragma once

#include "measurable.hpp"
#include "tensor.cuh"
#include "real_type.h"
#include "stopwatchDeviceSync.hpp"
// #include "cudaStructures.hpp"
#include "MeasurementWriter.h"

class CudaMeasurement : public Measurable
{
public:
    CudaMeasurement(const CudaTensor<real, 3>& emomM,
                    const CudaTensor<real, 3>& emom,
                    const CudaTensor<real, 2>& mmom,
                    bool alwaysCopy = false);
    ~CudaMeasurement() override;
    void measure(std::size_t mstep) override;
    void flushMeasurements(std::size_t mstep) override;

    struct AverageMagnetizationData { real m_x, m_y, m_z, m, m_stdv; };

    struct BinderCumulantData
    {
        real avrgmcum;      // Cumulated average of m
        real avrgm2cum;     // Cumulated average of m^2
        real avrgm4cum;     // Cumulated average of m^4
        real binderc;       // Binder cumulant
        real pmsusc;        // Susceptibility
        real cv;            // Specific heat
        real avrgecum;      // Cumulated average of E
        // real avrge2cum;     // Cumulated average of E^2
        real avrgetcum;     // Cumulated average of E_xc
        real avrgelcum;     // Cumulated average of E_LSF

        real cumuw;         // Weight for current sample to cumulant
        real cumutotw;      // Sum of all cumulant weights
        uint Navrgcum;      // Counter for number of cumulated averages
    };

private:
    void measureAverageMagnetization(std::size_t mstep);
    void measureBinderCumulant(std::size_t mstep);
    void measureSkyrmionNumber(std::size_t mstep);


private:
    const CudaTensor<real, 3>& emomM;
    const CudaTensor<real, 3>& emom;
    const CudaTensor<real, 2>& mmom;
    StopwatchDeviceSync stopwatch;
    MeasurementWriter measurementWriter;

    // Average magnetization
    CudaTensor<AverageMagnetizationData, 1> mavg_buff_gpu;
    Tensor<AverageMagnetizationData, 1> mavg_buff_cpu;
    Tensor<std::size_t, 1> mavg_iter;
    uint mavg_count = 0;

    // Binder cumulant
    CudaTensor<BinderCumulantData, 1> cumu_buff_gpu; // constant but tensor of rank 0 is not allowed, so rank 1 is size 1
    Tensor<BinderCumulantData, 1> cumu_buff_cpu;
    uint cumu_count = 0;

    // Skyrmion number
    CudaTensor<real, 3> dxyz_vec;
    CudaTensor<int, 2> dxyz_atom;
    CudaTensor<int, 1> dxyz_list;
    CudaTensor<real, 4> grad_mom; // 3 x 3 x N x M, is initialized to zeros in Fortran, no need to copy over
    CudaTensor<real, 2> skyno_buff_gpu;
    Tensor<real, 2> skyno_buff_cpu;
    Tensor<std::size_t, 1> skyno_iter;
    uint skyno_count = 0;
};

__global__ void naiveAverageMagnetization_kernel(const CudaTensor<real, 3> emomM,
                                                 CudaMeasurement::AverageMagnetizationData* d);