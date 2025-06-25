#pragma once

#include <cuda_runtime.h>
#include "c_headers.hpp"

#include "measurable.hpp"
#include "cudaEventPool.hpp"
#include "tensor.cuh"
#include "cudaParallelizationHelper.hpp"
#include "measurementQueue.hpp"
#include "real_type.h"
#include "stopwatch.hpp"
#include "stopwatchDeviceSync.hpp"

#if defined(USE_FAST_COPY)
#define DEFAULT_FAST_COPY true
#else
#define DEFAULT_FAST_COPY false
#endif

class FortranMeasurement : public Measurable {
    // Queue callback data struct
    struct queue_callback_data {
        queue_callback_data(FortranMeasurement* m, std::size_t s) : me(m), step(s) {
        }

        FortranMeasurement* me;
        std::size_t step;
    };

    // Queue callback
    static void queue_callback(cudaStream_t, cudaError_t, void* data);

    // Temporary device storage vectors
    CudaTensor<real, 3> tmp_emomM;
    CudaTensor<real, 3> tmp_emom;
    CudaTensor<real, 2> tmp_mmom;

    // Temporary host storage (pinned memory)
    Tensor<real, 3> pinned_emomM;
    Tensor<real, 3> pinned_emom;
    Tensor<real, 2> pinned_mmom;

    // Vectors to copy
    const CudaTensor<real, 3>& emomM;
    const CudaTensor<real, 3>& emom;
    const CudaTensor<real, 2>& mmom;

    Tensor<real, 3>& fortran_emomM;
    Tensor<real, 3>& fortran_emom;
    Tensor<real, 2>& fortran_mmom;

    // Event stack
    CudaEventPool eventPool;

    // Measure queue
    MeasurementQueue measurementQueue;

    // Timer
    StopwatchDeviceSync stopwatch;

    // Parallelization helper
    CudaParallelizationHelper& parallel;

    // Control flags
    bool alwaysCopy;
    bool fastCopy;

    // Helpers
    void queueMeasurement(std::size_t mstep);
    void copyQueueFast(std::size_t mstep);
    void copyQueueSlow(std::size_t mstep);

public:
    // TODO add flag for fast_copy
    FortranMeasurement(const CudaTensor<real, 3>& emomM, const CudaTensor<real, 3>& emom,
                    const CudaTensor<real, 2>& mmom, Tensor<real, 3>& f_emomM, Tensor<real, 3>& f_emom,
                    Tensor<real, 2>& f_mmom, bool fastCopy = DEFAULT_FAST_COPY,
                    bool alwaysCopy = false);
    ~FortranMeasurement() override;

    // Access methods
    void measure(std::size_t mstep) override;
    void flushMeasurements(std::size_t mstep) override;
};

