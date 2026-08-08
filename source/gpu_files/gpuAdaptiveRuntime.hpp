#pragma once

#include "gpu_wrappers.h"
#include "real_type.h"
#include "tensor.hpp"

#include <cstddef>
#include <cstdint>
#include <cmath>
#include <functional>
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

   // Optional CG-10 short-range operator inventory.  A null atomMoment keeps
   // the CG-09 storage-only path byte-for-byte unchanged.  Every array uses
   // the accepted Fortran ordering documented below and is copied during
   // initialize(), so construction and scratch storage are part of the same
   // preflight as the topology/runtime arrays.
   struct KernelInput {
      const double* atomMoment = nullptr;       // (atom)
      const int* atomAnisotropyAxisCount = nullptr; // (atom), 0..2
      const double* atomAnisotropyAxis = nullptr;   // (3,2,atom)
      const double* atomAnisotropyK1 = nullptr;     // (2,atom), joule
      const double* atomAnisotropyK2 = nullptr;     // (2,atom), joule
      const int* projectionBlock = nullptr;    // (8,atom), one-based
      const double* projectionWeight = nullptr;// (8,atom)
      std::size_t bonds = 0;
      const int* bondAtom = nullptr;            // (2,bond), one-based
      const double* bondMatrix = nullptr;       // (3,3,bond), joule
      std::size_t selectorEdges = 0;
      const int* selectorEdge = nullptr;        // (2,edge), one-based
      const double* inverseBlockTranspose = nullptr; // (3,3), m^-1
      const double* exchangeStiffness = nullptr;     // (3,3), J/m
      const double* spiralization = nullptr;         // (3,3), J/m^2
      const int* anisotropyAxisCount = nullptr;      // (block), 0..2
      const double* anisotropyAxis = nullptr;        // (3,2,block)
      const double* anisotropyK1 = nullptr;          // (2,block), J/m^3
      const double* anisotropyK2 = nullptr;          // (2,block), J/m^3
      double normalizationFloor = 1.0e-12;
      double magneticMomentSi = 9.2740100783e-24;
      double gammaPerTs = 0.0;
      double damping = 0.0;
   } kernels;
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
   unsigned char* polarizationUnsafeMask = nullptr;
   // RCG-03 diagnostic: worst (minimum) resultant/moment-sum ratio observed
   // over every channel/ensemble at each block, mirroring the Fortran
   // evaluate_polarization_gate block_ratio output; see its header comment
   // for the 0.0 undefined-direction floor convention.
   real* polarizationRatio = nullptr;
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

struct GpuAdaptiveEnergy {
   double atomisticBilinearJ = 0.0;
   double atomisticOnsiteJ = 0.0;
   double coarseExchangeJ = 0.0;
   double coarseSpiralizationJ = 0.0;
   double coarseAnisotropyJ = 0.0;
   double coarseExternalJ = 0.0;
   double dipoleJ = 0.0;
   double totalJ = 0.0;
};

struct GpuAdaptivePhaseMetrics {
   double atomisticMilliseconds = 0.0;
   double coarseMilliseconds = 0.0;
   double interfaceMilliseconds = 0.0;
   double selectorMilliseconds = 0.0;
   double polarizationMilliseconds = 0.0;
   double compactionMilliseconds = 0.0;
   double fftMilliseconds = 0.0;
   double integrationMilliseconds = 0.0;
};

struct GpuAdaptiveDiagnosticSnapshot {
   std::vector<int> blockState;
   std::vector<unsigned int> stateAge;
   std::vector<unsigned int> transitionEpoch;
   std::vector<real> selectorScores;
   std::vector<real> polarizationRatio;
   GpuAdaptiveEnergy energy{};
   double atomFieldSumT = 0.0;
   double atomFieldNorm2T2 = 0.0;
   double coarseFieldSumT = 0.0;
   double coarseFieldNorm2T2 = 0.0;
   double directionSum = 0.0;
   double directionNorm2 = 0.0;
};

// Basis-resolved uniform FFT field view. The convolution owner retains the
// padded storage; adaptive kernels consume it on the same device stream.
struct GpuAdaptiveUniformFftField {
   const real* paddedField = nullptr;
   std::size_t activeN1 = 0;
   std::size_t activeN2 = 0;
   std::size_t activeN3 = 0;
   std::size_t fftN1 = 0;
   std::size_t fftN2 = 0;
   std::size_t fftCells = 0;
   unsigned int basis = 0;
   real prefactorT = real(0);

   bool valid() const {
      return paddedField && activeN1 && activeN2 && activeN3 &&
             fftN1 && fftN2 && fftCells && basis &&
             std::isfinite(static_cast<double>(prefactorT));
   }
};

using GpuAdaptiveFftEvaluator =
   std::function<GpuAdaptiveUniformFftField(const real*)>;

struct GpuAdaptiveSelectorPolicy {
   real refineThreshold = real(0.25);
   real coarsenThreshold = real(0.10);
   unsigned int minimumDwellUpdates = 0;
   unsigned int bufferDilationBlocks = 0;
};

enum class GpuAdaptiveReconstruction {
   Aligned = 1,
   ConstrainedCone = 2
};

struct GpuAdaptiveReconstructionPolicy {
   GpuAdaptiveReconstruction scheme = GpuAdaptiveReconstruction::Aligned;
   real coneAngleRadians = real(0);
   real resultantTolerance = sizeof(real) == sizeof(double) ? real(1.0e-12) : real(2.0e-5);
   std::uint64_t globalSeed = 1;
};

// Owns immutable topology and mutable adaptive state for one GPU simulation
// lifetime.  It deliberately contains no operator polymorphism: CG-10
// kernels consume compact POD descriptors and this owner's tracked storage.
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

   // CG-10 device workflow.  Arrays use (3,atom,ensemble) and
   // (3,channel,block,ensemble) Fortran ordering.  The FFT field is an
   // optional already-dispatched all-grid field: adaptive state never changes
   // its resolution or ownership.
   bool kernelsReady() const { return kernelsReady_; }
   void restrictMoments(const real* atomDirection);
   void evaluateSelectorScores(const real* atomDirection);
   void evaluatePolarizationGate(real polarizationThreshold);
   const unsigned char* polarizationUnsafeBlockMask() const {
      return polarizationUnsafeBlockMask_.data();
   }
   // RCG-03 diagnostic: device pointer to the per-block ratio computed by
   // the same evaluatePolarizationGate() call; see GpuAdaptiveDeviceRuntime.
   const real* polarizationRatioBlock() const {
      return polarizationRatioBlock_.data();
   }
   void proposeSelectorState(const GpuAdaptiveSelectorPolicy& policy,
                             const unsigned char* hardAtomisticBlockMask = nullptr);
   void publishProposedState(real* atomDirection,
                             const GpuAdaptiveReconstructionPolicy& policy,
                             bool completeIntegrationStep,
                             const unsigned char* acceptedBlockMask = nullptr);
   GpuAdaptiveEnergy evaluateHybrid(const real* atomDirection,
                                    const real* externalCoarseField,
                                    const real* uniformFftDipoleField,
                                    real* atomField,
                                    real* coarseField,
                                    const GpuAdaptiveUniformFftField* basisResolvedFftField = nullptr);
   void integrateHeun(real timeStepSeconds, real* atomDirection,
                      const real* externalCoarseField = nullptr,
                      const real* uniformFftDipoleField = nullptr,
                      const GpuAdaptiveFftEvaluator& basisResolvedFftEvaluator = {});
   void synchronizeAtomicState(
      real* atomDirection,
      const GpuAdaptiveReconstructionPolicy& policy);
   void recordFftMilliseconds(double elapsed);
   const GpuAdaptivePhaseMetrics& phaseMetrics() const { return phaseMetrics_; }
   void resetPhaseMetrics() { phaseMetrics_ = {}; }

   bool ready() const { return ready_; }
   std::size_t allocatedBytes() const { return allocatedBytes_; }
   GPU_STREAM_T stream() const { return stream_; }
   const GpuAdaptiveDeviceTopology& deviceTopology() const { return deviceTopology_; }
   const GpuAdaptiveDeviceRuntime& deviceRuntime() const { return deviceRuntime_; }
   const GpuAdaptiveCompactionMetrics& compactionMetrics() const { return metrics_; }

   // Explicit diagnostics/test path.  It is not used by selector updates.
   GpuAdaptiveWorkSnapshot downloadWorkSnapshot();
   GpuAdaptiveDiagnosticSnapshot diagnosticSnapshot(const real* atomDirection);

private:
   void allocate(const GpuAdaptiveTopologyInput&, const GpuAdaptiveRuntimeInput&);
   void uploadTopology(const GpuAdaptiveTopologyInput&);
   void uploadRuntime(const GpuAdaptiveTopologyInput&, const GpuAdaptiveRuntimeInput&);
   void launchCompaction();
   void refreshDeviceDescriptors();
   double finishPhase(double& accumulator);

   bool ready_ = false;
   bool streamCreated_ = false;
   bool startEventCreated_ = false;
   bool endEventCreated_ = false;
   bool phaseEventsCreated_ = false;
   bool kernelsReady_ = false;
   std::size_t allocatedBytes_ = 0;
   std::size_t atoms_ = 0;
   std::size_t blocks_ = 0;
   std::size_t dynamicChannels_ = 0;
   std::size_t ensembles_ = 0;
   GPU_STREAM_T stream_{};
   GPU_EVENT_T updateStart_{};
   GPU_EVENT_T updateEnd_{};
   GPU_EVENT_T phaseStart_{};
   GPU_EVENT_T phaseEnd_{};
   GpuAdaptiveCompactionMetrics metrics_{};
   GpuAdaptivePhaseMetrics phaseMetrics_{};
   GpuAdaptiveEnergy lastEnergy_{};
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
   GpuTensor<unsigned int, 1> workCounts_, compactionScanA_, compactionScanB_;

   // CG-10 immutable operator data and preflight-visible scratch.
   std::size_t bonds_ = 0;
   std::size_t selectorEdges_ = 0;
   real normalizationFloor_ = real(0);
   real magneticMomentSi_ = real(0);
   real gammaPerTs_ = real(0);
   real damping_ = real(0);
   GpuTensor<real, 1> atomMoment_;
   GpuTensor<int, 1> atomAnisotropyAxisCount_;
   GpuTensor<real, 1> atomAnisotropyAxis_, atomAnisotropyK1_, atomAnisotropyK2_;
   GpuTensor<int, 1> projectionBlock_, bondAtom_, selectorEdge_;
   GpuTensor<real, 1> projectionWeight_, bondMatrix_;
   GpuTensor<real, 1> inverseBlockTranspose_, exchangeStiffness_, spiralization_;
   GpuTensor<int, 1> anisotropyAxisCount_;
   GpuTensor<real, 1> anisotropyAxis_, anisotropyK1_, anisotropyK2_;
   GpuTensor<real, 1> ghostDirection_, projectionNorm_, atomFieldScratch_;
   GpuTensor<real, 1> coarseFieldScratch_, predictorAtom_, predictorCoarse_;
   GpuTensor<real, 1> initialAtomField_, initialCoarseField_, predictorCoarseField_;
   GpuTensor<real, 1> energyTerms_;
   GpuTensor<unsigned char, 1> acceptedBlockMask_;
   GpuTensor<unsigned char, 1> polarizationUnsafeBlockMask_;
   GpuTensor<real, 1> polarizationRatioBlock_;

   GpuAdaptiveDeviceTopology deviceTopology_{};
   GpuAdaptiveDeviceRuntime deviceRuntime_{};
};
