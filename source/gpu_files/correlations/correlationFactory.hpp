#pragma once

#include "c_headers.hpp"
#include "tensor.hpp"
#include "real_type.h"
#include "fortranData.hpp"
#include "gpuStructures.hpp"
#include "correlation.hpp"

#include "fortranCorrelation.hpp"

#if defined(HIP_V)
#include "gpuCorrelations.hpp"
#elif defined(CUDA_V)
#include "gpuCorrelations.cuh"
#endif

#include <iostream>

class CorrelationFactory
{
public:
    // could be moved to a .cu file, but the function was so short, so I implemented it
    // directly in the header
    static std::unique_ptr<Correlation> create(const deviceLattice& gpuLattice, hostLattice& cpuLattice, const Flag Flags, const SimulationParameters SimParam, const hostCorrelations& cpuCorrelations)
    {
        if (*FortranData::do_gpu_correlations == 'Y')
        {
            std::cout << "GpuCorrelation used" << std::endl;
            return std::make_unique<GpuCorrelations>(Flags, SimParam, gpuLattice, cpuCorrelations);
        }
        else
        {
            std::cout << "FortranCorrelation used" << std::endl;
            return std::make_unique<FortranCorrelation>(
                gpuLattice.emomM,
                gpuLattice.emom,
                gpuLattice.mmom,
                cpuLattice.emomM,
                cpuLattice.emom,
                cpuLattice.mmom
            );
        }
    }
};
