#pragma once

#include "measurable.hpp"
#include "tensor.cuh"


class SkyrmionNumber : public Measurable
{
public:

    SkyrmionNumber(const CudaTensor<real, 3>& emomM, const CudaTensor<real, 3>& emom);
    ~SkyrmionNumber() override;
    void measure(std::size_t mstep) override;
    void flushMeasurements(std::size_t mstep) override;

    // fortran globals used:
    // indxb_skyno, skynob, sk_num_cum, sk_avrg, sk_var, sk_num_num?

    // values for printing:
    // "# Iter","Skx num", "Skx avg", "Skx std"
    // indxb_skyno(k), skynob(k), sk_avrg, sk_var/sk_num_num


private:
    const CudaTensor<real, 3>& emomM;
    const CudaTensor<real, 3>& emom;
    const uint skyno_step;
    const uint buffer_size;
    CudaTensor<real, 3> dxyz_vec;
    CudaTensor<uint, 2> dxyz_atom;
    CudaTensor<uint, 1> dxyz_list;

    CudaTensor<real, 4> grad_mom; // 3 x 3 x N x M
    CudaTensor<real, 1> skynob;
    CudaTensor<real, 1> sk_avrg;
    CudaTensor<real, 1> sk_var;
    Tensor<real, 1> indxb_skyno;
    uint buffer_count = 0;
    real sk_num_num = 0;
};



