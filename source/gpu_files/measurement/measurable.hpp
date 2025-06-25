#pragma once

#include <memory>
#include "c_headers.hpp"
#include "tensor.cuh"
#include "real_type.h"

class Measurable
{
public:
    virtual void measure(std::size_t mstep) = 0;
    virtual void flushMeasurements(std::size_t mstep) = 0;
    // virtual void copyToFortran() = 0;
    virtual ~Measurable() = default;

protected:
    Measurable() = default;
};










