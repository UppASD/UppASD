#pragma once

#include "c_headers.hpp"
#include "tensor.hpp"
#include "real_type.h"
#include "fortranData.hpp"
#include "gpuStructures.hpp"
#include "hamiltonian.hpp"
#include "gpuMacroHamiltonianCalculations.hpp"
#include "gpuHamiltonianCalculations.hpp"
#include "measurementQueue.hpp"
#include <iostream>

class HamiltonianFactory
{
public:
    // could be moved to a .cu file, but the function was so short, so I implemented it
    // directly in the header
    static std::unique_ptr<Hamiltonian> create(const Flag Flags, const SimulationParameters SimParam, 
                                        const deviceHamiltonian& gpuHamiltonian, const deviceMacrocell& gpuMacro)
    {
        if (*FortranData::do_macro_cells == 'Y')
        {
            std::cout << "Macro cells used" << std::endl;
            return std::make_unique<GpuMacroHamiltonianCalculations>(Flags, SimParam, gpuHamiltonian, gpuMacro);
        }
        else
        {
            std::cout << "Macro cells NOT used" << std::endl;
            return std::make_unique<GpuHamiltonianCalculations>(Flags, SimParam, gpuHamiltonian,);
        }
    }
};
