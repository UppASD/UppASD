#pragma once

#include <cuda_runtime.h>
#include <vector>
#include <memory>
#include "measurable.hpp"
#include "tensor.cuh"
#include "real_type.h"
#include "stopwatchDeviceSync.hpp"
#include "cudaStructures.hpp"

class CudaMeasurement : public Measurable
{
public:
    CudaMeasurement(const CudaTensor<real, 3>& emomM,
                    const CudaTensor<real, 3>& emom,
                    const CudaTensor<real, 2>& mmom,
                    bool alwaysCopy = false);
    ~CudaMeasurement() override = default;
    void measure(std::size_t mstep) override;
    void flushMeasurements(std::size_t mstep) override;

private:
    std::vector<std::unique_ptr<Measurable>> measurables;
    StopwatchDeviceSync stopwatch;
    bool alwaysCopy;
};

