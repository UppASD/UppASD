#pragma once

#include "gpu_wrappers.h"
#include "real_type.h"
#include "tensor.hpp"

#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

// GPU contract for the canonical BlockTopology Fortran type.  Integer maps
// retain their Fortran meaning (one-based ids, except for the zero-based CSR
// offsets).  The same plain descriptor is compiled by the CUDA and HIP paths.
struct GpuAdaptiveTopologyInput {
   int geometryMode = 0;
   std::size_t atoms = 0;
   std::size_t blocks = 0;
   std::size_t basis = 0;
   std::size_t fftChannelsPerBlock = 0;
   std::size_t fftGridChannels = 0;
   std::size_t dynamicChannels = 0;
   std::size_t ensembles = 0;
   const int* repetitionShape = nullptr;
   const int* blockShape = nullptr;
   const int* blockGrid = nullptr;
   const double* cellVectors = nullptr;
   const double* blockVectors = nullptr;

   const int* atomToBlock = nullptr;
   const int* atomToBasis = nullptr;
   const int* atomToDynamicChannel = nullptr;
   const int* atomToFftChannel = nullptr;
   const int* atomToFftGridIndex = nullptr;
   const int* basisToDynamicChannel = nullptr;
   const int* basisToFftChannel = nullptr;
   const int* blockAtomCount = nullptr;
   const int* blockAtomOffset = nullptr;
   const int* blockAtoms = nullptr;
   const int* blockGridCoordinate = nullptr;
   const int* blockBasisPopulation = nullptr;
   const int* blockFftChannelPopulation = nullptr;
   const int* blockDynamicChannelPopulation = nullptr;
   const double* blockCenter = nullptr;
   const double* blockVolume = nullptr;
};

struct GpuAdaptiveRuntimeInput {
   // CG-06 ownership states: 0 coarse, 1 buffer/interface, 2 fine.
   const int* blockState = nullptr;
   const int* pendingState = nullptr;
   const unsigned int* stateAge = nullptr;
   const unsigned int* transitionEpoch = nullptr;
   std::size_t selectorCriteria = 0;
   const double* selectorScores = nullptr;

   // Fortran order is (3, channel, block, ensemble).
   const double* coarseMoment = nullptr;
   const double* coarseDirection = nullptr;
   const double* coarseField = nullptr;
   // Fortran order is (channel, block, ensemble).
   const double* channelMomentSum = nullptr;
};

struct GpuAdaptiveDeviceTopology {
   int geometryMode = 0;
   std::size_t atoms = 0;
   std::size_t blocks = 0;
   std::size_t basis = 0;
   std::size_t fftChannelsPerBlock = 0;
   std::size_t fftGridChannels = 0;
   std::size_t dynamicChannels = 0;
   std::size_t ensembles = 0;
   int repetitionShape[3] = {};
   int blockShape[3] = {};
   int blockGrid[3] = {};
   real cellVectors[9] = {};
   real blockVectors[9] = {};
   const int* atomToBlock = nullptr;
   const int* atomToBasis = nullptr;
   const int* atomToDynamicChannel = nullptr;
   const int* atomToFftChannel = nullptr;
   const int* atomToFftGridIndex = nullptr;
   const int* basisToDynamicChannel = nullptr;
   const int* basisToFftChannel = nullptr;
   const int* blockAtomCount = nullptr;
   const int* blockAtomOffset = nullptr;
   const int* blockAtoms = nullptr;
   const int* blockGridCoordinate = nullptr;
   const int* blockBasisPopulation = nullptr;
   const int* blockFftChannelPopulation = nullptr;
   const int* blockDynamicChannelPopulation = nullptr;
   const real* blockCenter = nullptr;
   const real* blockVolume = nullptr;
};

struct GpuAdaptiveDeviceRuntime {
   std::size_t selectorCriteria = 0;
   int* blockState = nullptr;
   int* pendingState = nullptr;
   unsigned int* stateAge = nullptr;
   unsigned int* transitionEpoch = nullptr;
   unsigned char* atomisticBlockMask = nullptr;
   unsigned char* coarseBlockMask = nullptr;
   unsigned char* atomisticAtomMask = nullptr;
   unsigned char* interfaceAtomMask = nullptr;
   int* activeAtomList = nullptr;
   int* activeBlockList = nullptr;
   int* interfaceAtomList = nullptr;
   unsigned int* workCounts = nullptr; // active atoms, active blocks, interface atoms
   real* selectorScores = nullptr;
   real* coarseMoment = nullptr;
   real* coarseDirection = nullptr;
   real* coarseField = nullptr;
   real* channelMomentSum = nullptr;
};

struct GpuAdaptiveWorkSnapshot {
   std::vector<int> activeAtoms;
   std::vector<int> activeBlocks;
   std::vector<int> interfaceAtoms;
   std::vector<unsigned char> atomisticBlockMask;
   std::vector<unsigned char> coarseBlockMask;
   std::vector<unsigned char> atomisticAtomMask;
   std::vector<unsigned char> interfaceAtomMask;
};

struct GpuAdaptiveCompactionMetrics {
   std::uint64_t hostSynchronizations = 0;
   std::uint64_t blockBytesUploaded = 0;
   double elapsedMilliseconds = 0.0;       // complete staged update wall time
   double hostWaitMilliseconds = 0.0;      // time blocked in end-event synchronization
   double deviceMilliseconds = 0.0;        // GPU event interval
};

// Owns immutable topology and mutable adaptive state for one GPU simulation
// lifetime.  It deliberately contains no operator polymorphism: later CG-10
// kernels consume the compact POD descriptors above.
//
// Stream ordering: the owner creates one stream before allocation; topology
// upload, runtime staging, mask rebuilds, and future adaptive kernels are
// ordered on that stream.  initialize() synchronizes it once before publishing
// ready=true, updateBlockState() synchronizes only its measured end event so
// the pinned block staging buffer can be reused, and release() synchronizes
// before freeing.  Consumers must launch on stream() (or explicitly establish
// an event dependency) rather than assuming legacy-default-stream ordering.
class GpuAdaptiveRuntime {
public:
   GpuAdaptiveRuntime() = default;
   ~GpuAdaptiveRuntime();
   GpuAdaptiveRuntime(const GpuAdaptiveRuntime&) = delete;
   GpuAdaptiveRuntime& operator=(const GpuAdaptiveRuntime&) = delete;

   static bool validate(const GpuAdaptiveTopologyInput& topology,
                        const GpuAdaptiveRuntimeInput& runtime,
                        std::size_t expectedAtoms,
                        std::size_t expectedEnsembles,
                        std::string& diagnostic);
   static std::size_t estimateBytes(const GpuAdaptiveTopologyInput& topology,
                                    const GpuAdaptiveRuntimeInput& runtime);

   void initialize(const GpuAdaptiveTopologyInput& topology,
                   const GpuAdaptiveRuntimeInput& runtime,
                   std::size_t expectedAtoms,
                   std::size_t expectedEnsembles);
   void release();

   // This is the only initial CPU-selector synchronization path.  It stages
   // and uploads exactly blocks*sizeof(int), launches compaction on the owned
   // stream, and measures the localized event synchronization.
   void updateBlockState(const int* blockState, std::size_t count);

   bool ready() const { return ready_; }
   std::size_t allocatedBytes() const { return allocatedBytes_; }
   GPU_STREAM_T stream() const { return stream_; }
   const GpuAdaptiveDeviceTopology& deviceTopology() const { return deviceTopology_; }
   const GpuAdaptiveDeviceRuntime& deviceRuntime() const { return deviceRuntime_; }
   const GpuAdaptiveCompactionMetrics& compactionMetrics() const { return metrics_; }

   // Explicit diagnostics/test path.  It is not used by selector updates.
   GpuAdaptiveWorkSnapshot downloadWorkSnapshot();

private:
   void allocate(const GpuAdaptiveTopologyInput&, const GpuAdaptiveRuntimeInput&);
   void uploadTopology(const GpuAdaptiveTopologyInput&);
   void uploadRuntime(const GpuAdaptiveTopologyInput&, const GpuAdaptiveRuntimeInput&);
   void launchCompaction();
   void refreshDeviceDescriptors();

   bool ready_ = false;
   bool streamCreated_ = false;
   bool startEventCreated_ = false;
   bool endEventCreated_ = false;
   std::size_t allocatedBytes_ = 0;
   std::size_t atoms_ = 0;
   std::size_t blocks_ = 0;
   std::size_t dynamicChannels_ = 0;
   std::size_t ensembles_ = 0;
   GPU_STREAM_T stream_{};
   GPU_EVENT_T updateStart_{};
   GPU_EVENT_T updateEnd_{};
   GpuAdaptiveCompactionMetrics metrics_{};
   std::vector<std::vector<real>> convertedStaging_;

   Tensor<int, 1> stagedBlockState_;

   GpuTensor<int, 1> atomToBlock_, atomToBasis_, atomToDynamicChannel_;
   GpuTensor<int, 1> atomToFftChannel_, atomToFftGridIndex_;
   GpuTensor<int, 1> basisToDynamicChannel_, basisToFftChannel_;
   GpuTensor<int, 1> blockAtomCount_, blockAtomOffset_, blockAtoms_;
   GpuTensor<int, 1> blockGridCoordinate_;
   GpuTensor<int, 1> blockBasisPopulation_, blockFftChannelPopulation_;
   GpuTensor<int, 1> blockDynamicChannelPopulation_;
   GpuTensor<real, 1> blockCenter_, blockVolume_;

   GpuTensor<int, 1> blockState_, pendingState_;
   GpuTensor<unsigned int, 1> stateAge_, transitionEpoch_;
   GpuTensor<unsigned char, 1> atomisticBlockMask_, coarseBlockMask_;
   GpuTensor<unsigned char, 1> atomisticAtomMask_, interfaceAtomMask_;
   GpuTensor<real, 1> selectorScores_;
   GpuTensor<real, 1> coarseMoment_, coarseDirection_, coarseField_;
   GpuTensor<real, 1> channelMomentSum_;
   GpuTensor<int, 1> activeAtomList_, activeBlockList_, interfaceAtomList_;
   GpuTensor<unsigned int, 1> workCounts_;

   GpuAdaptiveDeviceTopology deviceTopology_{};
   GpuAdaptiveDeviceRuntime deviceRuntime_{};
};
