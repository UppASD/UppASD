#pragma once
#include "gpuStructures.hpp"


class Correlation
{
public:
    virtual void measure(std::size_t mstep) = 0;
    virtual void flushCorrelations(hostCorrelations& cpuCorrelations, std::size_t mstep) = 0;
    // virtual void copyToFortran() = 0;
    virtual ~Correlation() = default;

protected:
    Correlation() = default;
};










