#pragma once
#include "gpuStructures.hpp"
#include "tensor.hpp"



class Solver
{
public:
    virtual bool initiate(const SimulationParameters SimParam) = 0;
    virtual bool initiateConstants(const SimulationParameters SimParam, const Tensor<real, 1>temperature) = 0;
    virtual void evolveFirst(deviceLattice& gpuLattice) = 0;
    virtual void evolveSecond(deviceLattice& gpuLattice) = 0;

    virtual ~Solver() = default;

protected:
    Solver() = default;
};










