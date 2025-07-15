#pragma once

#include "measurable.hpp"


class Energy : public Measurable
{
public:
    Energy();
    ~Energy() override;
    void measure(std::size_t mstep) override;
    void flushMeasurements(std::size_t mstep) override;

private:
    const uint avrg_step;
    const uint buffer_size;
    CudaTensor<real, 3> eavg_buff;
};
