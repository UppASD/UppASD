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
   // RCG-09C: non-atomistic atoms that are an endpoint of at least one live
   // bond -- the only coarse atoms the atomistic reaction field can reach.
   unsigned char* ghostAtomMask = nullptr;
   int* activeAtomList = nullptr;
   int* activeBlockList = nullptr;
   int* interfaceAtomList = nullptr;
   int* ghostAtomList = nullptr;
   // RCG-09C: the compact live-bond list -- unique bonds that remain
   // atomistically owned, in ascending bond order.
   int* activeBondList = nullptr;
   // active atoms, active (coarse) blocks, interface atoms, ghost atoms,
   // live bonds
   unsigned int* workCounts = nullptr;
   real* selectorScores = nullptr;
   real* coarseMoment = nullptr;
   real* coarseDirection = nullptr;
   real* coarseField = nullptr;
   real* channelMomentSum = nullptr;
   // CGP-03: one-element device counter, atomicAdd'd by publishAdaptiveState
   // for every block whose blockState actually changes (not merely
   // proposed/rolled-back). Cleared before each publishAdaptiveState launch;
   // read back by GpuAdaptiveRuntime::transitionsPublishedLastCall().
   unsigned int* transitionsPublished = nullptr;
};

struct GpuAdaptiveWorkSnapshot {
   std::vector<int> activeAtoms;
   std::vector<int> activeBlocks;
   std::vector<int> interfaceAtoms;
   std::vector<int> ghostAtoms;
   std::vector<int> activeBonds;
   std::vector<unsigned char> atomisticBlockMask;
   std::vector<unsigned char> coarseBlockMask;
   std::vector<unsigned char> atomisticAtomMask;
   std::vector<unsigned char> interfaceAtomMask;
   std::vector<unsigned char> ghostAtomMask;
};

struct GpuAdaptiveCompactionMetrics {
   std::uint64_t hostSynchronizations = 0;
   std::uint64_t blockBytesUploaded = 0;
   double elapsedMilliseconds = 0.0;       // complete staged update wall time
   double hostWaitMilliseconds = 0.0;      // time blocked in end-event synchronization
   double deviceMilliseconds = 0.0;        // GPU event interval
   // RCG-09C: compaction now publishes its work counts to the host so launch
   // sizes can be cut to the live lists.  That readback is one stream
   // synchronization plus a five-word copy per rebuild -- not per timestep --
   // and is measured here rather than folded into the numbers above, so the
   // price of live-list launching is visible next to its benefit.
   std::uint64_t rebuilds = 0;
   std::uint64_t workCountReadbacks = 0;
   double workCountReadbackMilliseconds = 0.0;
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
   // RCG-06C (F-18): host wall-clock duration of the complete adaptive step
   // (GpuSimulation::advanceAdaptiveStep, entry to exit), accumulated
   // independently of the device-event phase timers above. The difference
   // between this and the sum of the other fields is genuinely unaccounted
   // time, not assumed to be ~0; see recordStepWallMilliseconds().
   double stepWallMilliseconds = 0.0;

   // RCG-08: kernel launches issued per phase, so a reader can distinguish
   // "this phase is slow" from "this phase is launch-bound", and can see the
   // compaction launch count stop growing with system size once the
   // Hillis--Steele sweep (F-12) was replaced by the hierarchical scan.
   unsigned long long atomisticLaunches = 0;
   unsigned long long coarseLaunches = 0;
   unsigned long long interfaceLaunches = 0;
   unsigned long long selectorLaunches = 0;
   unsigned long long polarizationLaunches = 0;
   unsigned long long compactionLaunches = 0;
   unsigned long long integrationLaunches = 0;
   // Device synchronizations performed by the phase timers themselves. Each
   // finishPhase() blocks the host on an event, so this counts real host
   // waits inside the adaptive step, not just instrumentation bookkeeping.
   unsigned long long phaseSynchronizations = 0;
   // RCG-09 (RCG-08-FU2): false when the per-phase device timers were disabled
   // for this measurement, in which case every *Milliseconds field above except
   // stepWallMilliseconds is zero because it was never measured -- not because
   // the phase was free.  Consumers must report "not measured", never 0.0.
   bool phaseTimingEnabled = true;
};

// CGP-00 evidence-only diagnostic: the real per-block FP64 partials each
// production energy kernel writes for the most recent evaluateHybrid() call,
// sliced to only the blocks that call actually wrote (stale entries from a
// previous, larger step are never included). Not part of the production
// timestep's data path -- this exists so the CGP-00 hierarchical-precision
// evidence harness can be fed the real magnitude/sign distribution
// production physics produces, without building parallel precision-variant
// copies of the physics kernels themselves. See
// GpuAdaptiveRuntime::downloadEnergyPartialsSnapshot().
struct GpuAdaptiveEnergyPartialsSnapshot {
   std::vector<double> atomisticBilinearJ;   // term 0: evaluateAdaptiveAtomisticBonds
   std::vector<double> atomisticOnsiteJ;     // term 1: evaluateAdaptiveAtomisticOnsite
   std::vector<double> coarseExchangeJ;      // term 2: evaluateAdaptiveCoarseTensor
   std::vector<double> coarseSpiralizationJ; // term 3: evaluateAdaptiveCoarseTensor
   std::vector<double> coarseAnisotropyJ;    // term 4: finalizeAdaptiveCoarseLocal
   std::vector<double> coarseExternalJ;      // term 5: finalizeAdaptiveCoarseLocal
   std::vector<double> dipoleJ;              // term 6: addAdaptiveDipole / basis-resolved
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
   double atomisticDirectionSum = 0.0;
   double atomisticDirectionNorm2 = 0.0;
   double coarseDirectionSum = 0.0;
   double coarseDirectionNorm2 = 0.0;
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
   // Per-axis (x,y,z) buffer dilation width in blocks, matching the CPU's
   // directional statichybridoperator%buffer_width_blocks(3) exactly
   // (RCG-05D). Index order matches blockGridCoordinate/blockGrid: 0=x,
   // 1=y, 2=z.
   unsigned int bufferDilationBlocks[3] = {0, 0, 0};
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

// CGP-03: named reasons for calling materializeCoarseAtoms(), so every
// full-lattice atom-direction commit in production code is explicit and
// countable rather than an unnamed call to synchronizeAtomicState() buried
// among several call sites (see docs/CGP-03_LAZY_MATERIALIZATION_EVIDENCE.md
// Part 3). This task keeps the OrdinaryStep commit unconditional -- the
// evidence doc explains the GpuMomentUpdater/measurement-schedule coupling
// that makes skipping it unsafe without further, larger changes -- and makes
// only the Transition commit conditional on a transition having actually
// been published.
enum class GpuAdaptiveMaterializationReason {
   OrdinaryStep,
   Transition
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
   // CGP-04: fullFieldClear selects which clear kernels evaluateHybridImpl()
   // launches (see gpuAdaptiveRuntime.cpp and
   // docs/CGP-04_COMPACT_FIELD_CLEAR_EVIDENCE.md for the validity contract).
   // Default true preserves the exact pre-CGP-04 full-system clear for every
   // existing caller (tests, benchmarks, and any call site that does not
   // explicitly opt in); production's ordinary-step call sites in
   // gpuSimulation.cpp pass false except when a full materialization is
   // already required for another reason (diagnostics>=3, or the next
   // iteration's measurement/print cadence).
   GpuAdaptiveEnergy evaluateHybrid(const real* atomDirection,
                                    const real* externalCoarseField,
                                    const real* uniformFftDipoleField,
                                    real* atomField,
                                    real* coarseField,
                                    const GpuAdaptiveUniformFftField* basisResolvedFftField = nullptr,
                                    bool fullFieldClear = true);
   // CGP-01: the field-only sibling of evaluateHybrid(), sharing the same
   // kernels (see evaluateHybridImpl<MeasureEnergy>) but compiled without any
   // energy contribution arithmetic, block/partial reduction, finalization
   // launch, or energy D2H readback.  lastEnergy_ is left untouched -- callers
   // that discard evaluateHybrid()'s return value (the ordinary Depondt
   // predictor/corrector field evaluations) should call this instead.
   void evaluateField(const real* atomDirection,
                      const real* externalCoarseField,
                      const real* uniformFftDipoleField,
                      real* atomField,
                      real* coarseField,
                      const GpuAdaptiveUniformFftField* basisResolvedFftField = nullptr,
                      bool fullFieldClear = true);
   // CGP-01 negative control: incremented exactly once per evaluateHybrid()
   // energy reduction/finalization/D2H readback, never by evaluateField().
   // Tests use this to prove the field-only path performs zero energy work.
   std::size_t energyEvaluationCount() const { return energyEvaluationCount_; }
   // Adaptive field assembly is separate from spin integration.  The caller
   // evaluates the total field at the initial and predictor states, while the
   // shared production Depondt integrator owns the fine-atom updates.
   real* coarseFieldData() { return coarseField_.data(); }
   const int* activeAtomList() const { return activeAtomList_.data(); }
   // RCG-09C: served from the host mirror refreshed by every compaction.  The
   // counts can only change when compaction runs, so this no longer costs a
   // stream synchronization and a device read per timestep.
   std::size_t activeAtomCount() const { return hostWorkCounts_[0]; }
   std::size_t activeBlockCount() const { return hostWorkCounts_[1]; }
   std::size_t interfaceAtomCount() const { return hostWorkCounts_[2]; }
   std::size_t ghostAtomCount() const { return hostWorkCounts_[3]; }
   std::size_t liveBondCount() const { return hostWorkCounts_[4]; }
   std::size_t bondCount() const { return bonds_; }
   void prepareCoarsePredictor(real timeStepSeconds,
                               const real* initialCoarseField);
   void correctCoarse(real timeStepSeconds,
                      const real* predictorCoarseField);
   void synchronize();
   void synchronizeAtomicState(
      real* atomDirection,
      const GpuAdaptiveReconstructionPolicy& policy);
   // CGP-03 Part 3: the one named entry point production code should call to
   // commit coarse-block state into atomDirection, instead of calling
   // synchronizeAtomicState() directly from scattered call sites. Internally
   // still synchronizeAtomicState() -- no reconstruction formula changes --
   // but every call now carries an explicit reason and is countable via
   // materializationCount().
   void materializeCoarseAtoms(
      GpuAdaptiveMaterializationReason reason, real* atomDirection,
      const GpuAdaptiveReconstructionPolicy& policy);
   std::size_t materializationCount(
      GpuAdaptiveMaterializationReason reason) const {
      return materializationCounts_[static_cast<std::size_t>(reason)];
   }
   // CGP-03: exact count of blocks whose blockState actually changed in the
   // most recent publishProposedState() call (0 if none did, e.g. a
   // selector-update step where every score stayed within its dwell/
   // threshold band). Read back as part of publishProposedState()'s existing
   // refreshHostWorkCounts() synchronization -- no extra host wait.
   std::size_t transitionsPublishedLastCall() const {
      return transitionsPublishedLastCall_;
   }
   void recordFftMilliseconds(double elapsed);
   // RCG-06C (F-18): called once per complete adaptive step with an
   // independent host wall-clock measurement of the whole step (not derived
   // from the device-event phase timers), so unaccounted time can be
   // reported rather than assumed zero.
   void recordStepWallMilliseconds(double elapsed);
   const GpuAdaptivePhaseMetrics& phaseMetrics() const { return phaseMetrics_; }
   const GpuAdaptiveEnergy& lastEnergy() const { return lastEnergy_; }
   void resetPhaseMetrics() {
      const bool timing = phaseMetrics_.phaseTimingEnabled;
      phaseMetrics_ = {};
      phaseMetrics_.phaseTimingEnabled = timing;
   }

   // RCG-09 (RCG-08-FU2).  Per-phase timing costs one blocking host
   // synchronization per phase boundary -- ten per step, measured at ~38% of
   // step wall time in RCG-08 SS6.5.  That cost is instrumentation, not
   // physics, so a wall-time measurement that leaves it enabled measures the
   // instrument.  Disabling it removes only the event synchronization; kernel
   // order, stream, and every launch are unchanged, and launch counters keep
   // working.  Default is enabled, so RCG-06C's accepted per-phase evidence
   // and every existing caller are bit-for-bit unaffected.
   void setPhaseTimingEnabled(bool enabled) {
      phaseMetrics_.phaseTimingEnabled = enabled;
   }
   bool phaseTimingEnabled() const { return phaseMetrics_.phaseTimingEnabled; }

   bool ready() const { return ready_; }
   std::size_t allocatedBytes() const { return allocatedBytes_; }
   GPU_STREAM_T stream() const { return stream_; }
   const GpuAdaptiveDeviceTopology& deviceTopology() const { return deviceTopology_; }
   const GpuAdaptiveDeviceRuntime& deviceRuntime() const { return deviceRuntime_; }
   const GpuAdaptiveCompactionMetrics& compactionMetrics() const { return metrics_; }

   // Explicit diagnostics/test path.  It is not used by selector updates.
   GpuAdaptiveWorkSnapshot downloadWorkSnapshot();
   GpuAdaptiveDiagnosticSnapshot diagnosticSnapshot(const real* atomDirection);
   // CGP-00 evidence-only: must be called immediately after evaluateHybrid()
   // (same contract as diagnosticSnapshot()'s dependence on lastEnergy_) so
   // the downloaded partials reflect that call, not a stale earlier one.
   GpuAdaptiveEnergyPartialsSnapshot downloadEnergyPartialsSnapshot(
      bool dipoleActive) const;

private:
   // CGP-01: evaluateHybrid()/evaluateField() are thin wrappers selecting the
   // MeasureEnergy=true/false instantiation of this template, which is the
   // single field implementation the design requirement asks for -- compile-
   // time specialized rather than branching inside every thread.  Defined and
   // instantiated only in gpuAdaptiveRuntime.cpp.
   template<bool MeasureEnergy>
   GpuAdaptiveEnergy evaluateHybridImpl(const real* atomDirection,
                                        const real* externalCoarseField,
                                        const real* uniformFftDipoleField,
                                        real* atomField,
                                        real* coarseField,
                                        const GpuAdaptiveUniformFftField* basisResolvedFftField,
                                        bool fullFieldClear);
   void allocate(const GpuAdaptiveTopologyInput&, const GpuAdaptiveRuntimeInput&);
   void uploadTopology(const GpuAdaptiveTopologyInput&);
   void uploadRuntime(const GpuAdaptiveTopologyInput&, const GpuAdaptiveRuntimeInput&);
   void launchCompaction();
   // RCG-09C: publish the freshly compacted work counts to hostWorkCounts_.
   // Called once per rebuild, never per timestep.
   void refreshHostWorkCounts();
   void refreshDeviceDescriptors();
   void beginPhase();
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
   // CGP-01 negative control; see energyEvaluationCount() above.
   std::size_t energyEvaluationCount_ = 0;
   // CGP-03: indexed by GpuAdaptiveMaterializationReason; see
   // materializationCount() above.
   std::size_t materializationCounts_[2] = {};
   // CGP-03: host mirror of transitionsPublished_, refreshed inside
   // publishProposedState()'s existing refreshHostWorkCounts() call.
   std::size_t transitionsPublishedLastCall_ = 0;
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
   // RCG-05E: dilateAdaptiveState's write target.  pendingState_ is read-only
   // for the duration of that kernel launch (see the invariant comment above
   // dilateAdaptiveState in gpuAdaptiveRuntime.cpp); the dilated result is
   // written here instead, then copied back into pendingState_ once the
   // kernel has completed, so no thread ever reads a location another thread
   // concurrently writes.
   GpuTensor<int, 1> dilatedState_;
   GpuTensor<unsigned int, 1> stateAge_, transitionEpoch_;
   GpuTensor<unsigned char, 1> atomisticBlockMask_, coarseBlockMask_;
   GpuTensor<unsigned char, 1> atomisticAtomMask_, interfaceAtomMask_;
   GpuTensor<unsigned char, 1> ghostAtomMask_;
   GpuTensor<real, 1> selectorScores_;
   GpuTensor<real, 1> coarseMoment_, coarseDirection_, coarseField_;
   GpuTensor<real, 1> channelMomentSum_;
   GpuTensor<int, 1> activeAtomList_, activeBlockList_, interfaceAtomList_;
   GpuTensor<int, 1> ghostAtomList_, activeBondList_;
   GpuTensor<unsigned int, 1> workCounts_, compactionScanA_, compactionScanB_;
   // RCG-09C: the live-bond scan is a separate, single-component sweep over
   // the unique-bond list.  Bonds outnumber atoms by the coordination number,
   // so folding them into the atom/block scan would make every atom and block
   // component that many times longer for no benefit.
   GpuTensor<unsigned int, 1> bondScanA_, bondScanB_, bondScanLevels_;
   std::vector<std::size_t> bondScanLevelItems_, bondScanLevelOffset_;
   // Host mirror of workCounts_, refreshed by refreshHostWorkCounts().
   unsigned int hostWorkCounts_[5] = {};
   // RCG-08 (F-12): tile-total storage for the hierarchical compaction scan.
   // One flat buffer holding levels 1..L back to back; scanLevelItems_[i] is
   // level i+1's per-component item count and scanLevelOffset_[i] its element
   // offset into compactionScanLevels_. Both are fixed at allocate() time and
   // are included in estimateBytes(), so the scan's temporary storage is
   // visible to the ordinary GPU memory preflight rather than allocated on
   // the hot path.
   GpuTensor<unsigned int, 1> compactionScanLevels_;
   std::vector<std::size_t> scanLevelItems_, scanLevelOffset_;

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
   GpuTensor<real, 1> coarseFieldScratch_, transitionBackup_, predictorCoarse_;
   GpuTensor<real, 1> initialCoarseField_, predictorCoarseField_;
   // RCG-06B (F-11): global energy accumulation is FP64 unconditionally,
   // independent of `real`'s build precision (SINGLE_PREC/DOUBLE_PREC), per
   // the parent blueprint's precision contract (section 6.6). Field/output
   // values remain explicit `real` everywhere else; only this 8-slot
   // reduction target is widened.
   GpuTensor<double, 1> energyTerms_;
   // RCG-09B: each energy-producing kernel writes one deterministic FP64
   // partial per launch block into this compact term-major buffer.  The
   // final ordered reduction consumes these partials; no adaptive energy
   // kernel writes a contended global scalar.
   GpuTensor<double, 1> energyPartials_;
   std::size_t energyPartialBlocks_ = 0;
   GpuTensor<unsigned char, 1> acceptedBlockMask_;
   GpuTensor<unsigned char, 1> polarizationUnsafeBlockMask_;
   GpuTensor<real, 1> polarizationRatioBlock_;
   // CGP-03: one-element device counter backing
   // GpuAdaptiveDeviceRuntime::transitionsPublished.
   GpuTensor<unsigned int, 1> transitionsPublished_;

   GpuAdaptiveDeviceTopology deviceTopology_{};
   GpuAdaptiveDeviceRuntime deviceRuntime_{};
};
