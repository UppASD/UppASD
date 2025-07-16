#pragma once
#include "measurable.hpp"
#include "MeasurementWriter.h"

struct AverageMagnetizationData
{
    real m_x;
    real m_y;
    real m_z;
    real m;
    real m_stdv;
};

class AverageMagnetization : public Measurable
{
public:
    explicit AverageMagnetization(const CudaTensor<real, 3>& emomM, dim3 cudaThreads = 1024);
    ~AverageMagnetization() override;
    void measure(std::size_t mstep) override;
    void flushMeasurements(std::size_t mstep) override;

    // void copyToFortran() override;

private:
    const CudaTensor<real, 3>& emomM;
    const uint mavg_buff_size;
    const uint N;
    const uint M;
    const uint avrg_step;
    const dim3 threads;
    const dim3 blocks;

    MeasurementWriter measurementWriter;

    CudaTensor<AverageMagnetizationData, 1> mavg_buff;
    Tensor<AverageMagnetizationData, 1> mavg_buff_cpu;
    Tensor<size_t, 1> indxb_mavg;

    CudaTensor<real, 2> mavg_block_buff;
//    Tensor<real, 3> mavg_buff_cpu;
    uint bcount_mavrg = 0;
    //Block variables

};