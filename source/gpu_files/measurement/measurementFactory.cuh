#pragma once

#include "c_headers.hpp"
#include "tensor.cuh"
#include "real_type.h"
#include "fortranData.hpp"
#include "cudaStructures.hpp"
#include "measurable.hpp"
#include "cudaMeasurement.cuh"
#include "fortranMeasurement.cuh"

#include <iostream>

class MeasurementFactory
{
public:
    // could be moved to a .cu file, but the function was so short, so I implemented it
    // directly in the header
    static std::unique_ptr<Measurable> create(const cudaLattice& gpuLattice, hostLattice& cpuLattice)
    {
        if (*FortranData::do_cuda_measurements == 'Y')
        {
            std::cout << "CudaMeasurement used" << std::endl;
            return std::make_unique<CudaMeasurement>(
                gpuLattice.emomM,
                gpuLattice.emom,
                gpuLattice.mmom
            );
        }
        else
        {
            std::cout << "FortranMeasurement used" << std::endl;
            return std::make_unique<FortranMeasurement>(
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
