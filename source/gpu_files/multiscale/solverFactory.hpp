#pragma once

#include "c_headers.hpp"
#include "tensor.hpp"
#include "real_type.h"
#include "fortranData.hpp"
#include "gpuStructures.hpp"
#include "solver.hpp"
//#include "gpuMeasurement.hpp"
//#include "fortranMeasurement.hpp"
//#include "measurementQueue.hpp"

#include <iostream>

class SolverFactory
{
public:
    // could be moved to a .cu file, but the function was so short, so I implemented it
    // directly in the header
    static std::unique_ptr<Solver> create(const deviceLattice& gpuLattice, hostLattice& cpuLattice, 
                                              const deviceEnergies& gpuEnergies, MeasurementQueue& mq,
                                              bool do_jtensor)
    {
        if (*FortranData::SDEalgh == 5)
        {
            std::cout << "Depondt integrator used" << std::endl;
            return std::make_unique<DepondtIntegratorMS>();
        }
        else if(*FortranData::SDEalgh == 1)
        {
            std::cout << "Midpoint Solver used" << std::endl;
            return std::make_unique<MidpointSolverMS>(
                gpuLattice.emomM,
                gpuLattice.emom,
                gpuLattice.mmom,
                cpuLattice.emomM,
                cpuLattice.emom,
                cpuLattice.mmom,
                mq
            );
        }
        else{
            std::cout << "Unavailable solver type on GPU" << std::endl;

        }
    }
};
