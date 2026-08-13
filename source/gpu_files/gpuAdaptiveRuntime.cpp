#include "gpuAdaptiveRuntime.hpp"

#include "base.hpp"
#include "gpuAdaptiveReconstructionRng.hpp"
#include "gpuAtomicDouble.hpp"
#include "measurement/memoryMeasurement.h"

#include <algorithm>
#include <chrono>
#include <climits>
#include <cmath>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <vector>

namespace {

constexpr int regularReplicatedCell = 1;
constexpr int coarseState = 0;
constexpr int bufferState = 1;
constexpr int fineState = 2;
constexpr real piValue = real(3.141592653589793238462643383279502884L);
constexpr unsigned int adaptiveThreads = 256;

inline bool anyBufferDilation(const unsigned int width[3]) {
   return width[0] > 0 || width[1] > 0 || width[2] > 0;
}

inline unsigned int adaptiveGrid(std::size_t workItems) {
   return static_cast<unsigned int>(
      std::max<std::size_t>(1, (workItems + adaptiveThreads - 1) /
                               adaptiveThreads));
}

__device__ inline std::size_t adaptiveThreadIndex() {
   return static_cast<std::size_t>(blockIdx.x) * blockDim.x + threadIdx.x;
}

// RCG-08: every launch site this task introduced or rewrote goes through this
// macro.  No HIP toolchain or device exists in any environment used across
// RCG-02..RCG-08, so the HIP spelling of a new launch cannot be compiled, let
// alone executed, here; writing each new launch twice by hand would make a
// silent CUDA/HIP divergence (different grid, different argument order) both
// easy to introduce and impossible to detect locally.  Routing them through
// one macro makes the two backends structurally identical by construction,
// which is the property section 5.5 of the remediation blueprint requires and
// the only one that can honestly be claimed without HIP hardware.
//
// Pre-existing launch sites are deliberately left in their explicit
// #if/#elif form: converting them would enlarge this patch's diff without
// changing behaviour, and section 5.2 asks for one defect class per patch.
#if defined(CUDA_V)
#define ADAPTIVE_LAUNCH(kernel, grid, stream, ...) \
   kernel<<<(grid), adaptiveThreads, 0, (stream)>>>(__VA_ARGS__)
#elif defined(HIP_V)
#define ADAPTIVE_LAUNCH(kernel, grid, stream, ...) \
   hipLaunchKernelGGL(kernel, dim3(grid), dim3(adaptiveThreads), 0, (stream), \
                      __VA_ARGS__)
#endif

// Level sizes for the hierarchical compaction scan: n_0 = scanItems and
// n_{k+1} = ceil(n_k / adaptiveThreads), stopping at the first level that fits
// in a single tile.  Returns n_1..n_L (empty when level 0 already fits).
// Shared by allocate() and estimateBytes() so the preflight cannot disagree
// with the actual allocation.
inline std::vector<std::size_t> compactionScanLevelItems(std::size_t scanItems) {
   std::vector<std::size_t> levels;
   std::size_t items = scanItems;
   while(items > adaptiveThreads) {
      items = (items + adaptiveThreads - 1) / adaptiveThreads;
      levels.push_back(items);
   }
   return levels;
}

struct AdaptiveKernelDevice {
   std::size_t bonds;
   std::size_t selectorEdges;
   real normalizationFloor;
   real magneticMomentSi;
   real gammaPerTs;
   real damping;
   const real* atomMoment;
   const int* atomAnisotropyAxisCount;
   const real* atomAnisotropyAxis;
   const real* atomAnisotropyK1;
   const real* atomAnisotropyK2;
   const int* projectionBlock;
   const real* projectionWeight;
   const int* bondAtom;
   const real* bondMatrix;
   const int* selectorEdge;
   const real* inverseBlockTranspose;
   const real* exchangeStiffness;
   const real* spiralization;
   const int* anisotropyAxisCount;
   const real* anisotropyAxis;
   const real* anisotropyK1;
   const real* anisotropyK2;
   real* ghostDirection;
   real* projectionNorm;
   real* atomFieldScratch;
   real* coarseFieldScratch;
   // RCG-06B (F-11): FP64 regardless of `real`'s build precision; see the
   // header comment on GpuAdaptiveRuntime::energyTerms_.
   double* energyTerms;
   real* transitionBackup;
};

__host__ __device__ inline std::size_t atomVectorIndex(
   int xyz, std::size_t atom, std::size_t ensemble, std::size_t atoms) {
   return static_cast<std::size_t>(xyz) + 3 * (atom + atoms * ensemble);
}

__host__ __device__ inline std::size_t coarseScalarIndex(
   std::size_t channel, std::size_t block, std::size_t ensemble,
   std::size_t channels, std::size_t blocks) {
   return channel + channels * (block + blocks * ensemble);
}

__host__ __device__ inline std::size_t coarseVectorIndex(
   int xyz, std::size_t channel, std::size_t block, std::size_t ensemble,
   std::size_t channels, std::size_t blocks) {
   return static_cast<std::size_t>(xyz) +
          3 * coarseScalarIndex(channel, block, ensemble, channels, blocks);
}

__device__ inline void crossDevice(const real a[3], const real b[3], real result[3]) {
   result[0] = a[1] * b[2] - a[2] * b[1];
   result[1] = a[2] * b[0] - a[0] * b[2];
   result[2] = a[0] * b[1] - a[1] * b[0];
}

__device__ inline real dotDevice(const real a[3], const real b[3]) {
   return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

__device__ inline real normDevice(const real a[3]) {
   return sqrt(dotDevice(a, a));
}

// Shared "never normalize noise" epsilon, matching restrict_channel_moments
// on the Fortran side.  A resultant is only trusted once its length exceeds
// 64*epsilon*scale; both restrictAdaptiveMoments and evaluatePolarizationGate
// must agree on this constant.
__device__ inline real epsilonDevice() {
#ifdef SINGLE_PREC
   return real(1.1920928955078125e-7);
#else
   return real(2.2204460492503131e-16);
#endif
}

__device__ inline std::size_t plusBlock(const GpuAdaptiveDeviceTopology& topology,
                                        std::size_t block, int direction) {
   int coordinate[3] = {
      topology.blockGridCoordinate[0 + 3 * block],
      topology.blockGridCoordinate[1 + 3 * block],
      topology.blockGridCoordinate[2 + 3 * block]
   };
   coordinate[direction] = (coordinate[direction] + 1) %
                           topology.blockGrid[direction];
   return static_cast<std::size_t>(
      coordinate[0] + topology.blockGrid[0] *
      (coordinate[1] + topology.blockGrid[1] * coordinate[2]));
}

__device__ inline void loadAtomVector(const real* data,
                                      const GpuAdaptiveDeviceTopology& topology,
                                      std::size_t atom, std::size_t ensemble,
                                      real value[3]) {
   for(int xyz = 0; xyz < 3; ++xyz)
      value[xyz] = data[atomVectorIndex(xyz, atom, ensemble, topology.atoms)];
}

__device__ inline void loadCoarseVector(const real* data,
                                        const GpuAdaptiveDeviceTopology& topology,
                                        std::size_t channel, std::size_t block,
                                        std::size_t ensemble, real value[3]) {
   for(int xyz = 0; xyz < 3; ++xyz)
      value[xyz] = data[coarseVectorIndex(xyz, channel, block, ensemble,
                                         topology.dynamicChannels, topology.blocks)];
}

// RCG-08 (F-09): restriction used to be a single-lane kernel that zeroed the
// whole coarse state, scattered every atom of every ensemble into it, and then
// finalized every channel/block/ensemble scalar -- O(atoms * ensembles)
// serialized onto thread (0,0).
//
// It is now one thread per (channel, block, ensemble) scalar.  Each thread
// owns its output scalar exclusively and gathers over that block's own CSR
// atom range, so the kernel needs no atomic and no separate clear pass: the
// accumulate and finalize halves both act on the single scalar the thread
// owns, and are therefore fused into one launch.
//
// The gather order is not merely race-free, it is bitwise faithful to the
// serial reference it replaces.  BlockTopology fills block_atoms by scanning
// atoms 1..Natom with a per-block cursor
// (source/CoarseGraining/blocktopology.f90 build loop), so each block's CSR
// slice is strictly ascending in atom index.  The old global atom loop
// contributed to a given (channel, block, ensemble) exactly the atoms of that
// block carrying that channel, in ascending atom order; summing the CSR slice
// visits the identical summands in the identical order.  Floating-point
// addition is order-sensitive but not set-sensitive, so the result is
// bit-identical rather than merely within tolerance -- which is what keeps the
// polarization gate's decisions (and hence adaptive transitions) deterministic
// after parallelization.
__global__ void restrictAdaptiveMoments(GpuAdaptiveDeviceTopology topology,
                                        GpuAdaptiveDeviceRuntime runtime,
                                        AdaptiveKernelDevice kernels,
                                        const real* atomDirection) {
   const std::size_t scalarCount =
      topology.dynamicChannels * topology.blocks * topology.ensembles;
   const std::size_t scalar = adaptiveThreadIndex();
   if(scalar >= scalarCount) return;
   const std::size_t channel = scalar % topology.dynamicChannels;
   const std::size_t block =
      (scalar / topology.dynamicChannels) % topology.blocks;
   const std::size_t ensemble =
      scalar / (topology.dynamicChannels * topology.blocks);

   const std::size_t begin =
      static_cast<std::size_t>(topology.blockAtomOffset[block]);
   const std::size_t end =
      static_cast<std::size_t>(topology.blockAtomOffset[block + 1]);
   real momentSum = real(0);
   real resultant[3] = {real(0), real(0), real(0)};
   for(std::size_t position = begin; position < end; ++position) {
      const std::size_t atom =
         static_cast<std::size_t>(topology.blockAtoms[position] - 1);
      const int rawChannel = topology.atomToDynamicChannel[atom];
      if(rawChannel <= 0 ||
         static_cast<std::size_t>(rawChannel - 1) != channel) continue;
      const real moment = kernels.atomMoment[atom];
      momentSum += moment;
      for(int xyz = 0; xyz < 3; ++xyz)
         resultant[xyz] += moment * atomDirection[
            atomVectorIndex(xyz, atom, ensemble, topology.atoms)];
   }
   runtime.channelMomentSum[scalar] = momentSum;

   if(runtime.blockState[block] == coarseState) {
      for(int xyz = 0; xyz < 3; ++xyz)
         runtime.coarseMoment[3 * scalar + xyz] =
            momentSum * runtime.coarseDirection[3 * scalar + xyz];
      return;
   }
   for(int xyz = 0; xyz < 3; ++xyz)
      runtime.coarseMoment[3 * scalar + xyz] = resultant[xyz];
   const real length = normDevice(resultant);
   const real scale = momentSum > real(1) ? momentSum : real(1);
   if(length > real(64) * epsilonDevice() * scale) {
      runtime.coarseDirection[3 * scalar] = resultant[0] / length;
      runtime.coarseDirection[3 * scalar + 1] = resultant[1] / length;
      runtime.coarseDirection[3 * scalar + 2] = resultant[2] / length;
   } else {
      runtime.coarseDirection[3 * scalar] = real(0);
      runtime.coarseDirection[3 * scalar + 1] = real(0);
      runtime.coarseDirection[3 * scalar + 2] = real(0);
   }
}

__device__ inline void atomicMaxSelector(real* address, real value) {
#ifdef SINGLE_PREC
   auto* bits = reinterpret_cast<unsigned int*>(address);
   unsigned int observed = *bits;
   while(value > __uint_as_float(observed)) {
      const unsigned int assumed = observed;
      observed = atomicCAS(bits, assumed, __float_as_uint(value));
      if(observed == assumed) break;
   }
#else
   auto* bits = reinterpret_cast<unsigned long long*>(address);
   unsigned long long observed = *bits;
   while(value > __longlong_as_double(static_cast<long long>(observed))) {
      const unsigned long long assumed = observed;
      observed = atomicCAS(
         bits, assumed,
         static_cast<unsigned long long>(__double_as_longlong(value)));
      if(observed == assumed) break;
   }
#endif
}

// RCG-06B (F-11): FP64 energy-accumulator atomicAdd. See
// gpuAtomicDouble.hpp (included above) for the CAS-loop implementation and
// why it does not use the native intrinsic unconditionally.

__global__ void clearSelectorAdaptiveScores(
   GpuAdaptiveDeviceTopology topology, GpuAdaptiveDeviceRuntime runtime) {
   const std::size_t block = adaptiveThreadIndex();
   if(block < topology.blocks)
      runtime.selectorScores[runtime.selectorCriteria * block] = real(0);
}

__global__ void selectorAdaptiveScores(GpuAdaptiveDeviceTopology topology,
                                       GpuAdaptiveDeviceRuntime runtime,
                                       AdaptiveKernelDevice kernels,
                                       const real* atomDirection) {
   const std::size_t work = adaptiveThreadIndex();
   const std::size_t count = kernels.selectorEdges * topology.ensembles;
   if(work >= count) return;
   const std::size_t edge = work % kernels.selectorEdges;
   const std::size_t ensemble = work / kernels.selectorEdges;
   const std::size_t atomI =
      static_cast<std::size_t>(kernels.selectorEdge[2 * edge] - 1);
   const std::size_t atomJ =
      static_cast<std::size_t>(kernels.selectorEdge[2 * edge + 1] - 1);
   real si[3], sj[3];
   loadAtomVector(atomDirection, topology, atomI, ensemble, si);
   loadAtomVector(atomDirection, topology, atomJ, ensemble, sj);
   const real q = real(1) - dotDevice(si, sj);
   const real score = q > real(0) ? q : real(0);
   const std::size_t blockI =
      static_cast<std::size_t>(topology.atomToBlock[atomI] - 1);
   const std::size_t blockJ =
      static_cast<std::size_t>(topology.atomToBlock[atomJ] - 1);
   atomicMaxSelector(
      &runtime.selectorScores[runtime.selectorCriteria * blockI], score);
   atomicMaxSelector(
      &runtime.selectorScores[runtime.selectorCriteria * blockJ], score);
}

// RCG-03 (F-14): a block is unsafe to coarsen whenever any dynamical
// channel/ensemble has no defined resultant direction (near-zero moment,
// never normalized -- see epsilonDevice()) or a resultant/moment-sum ratio
// below the accepted threshold.  One thread per block, looping internally
// over channel/ensemble, so there is no multi-writer race on the output
// byte (mirrors proposeAdaptiveState's one-thread-per-block style rather
// than selectorAdaptiveScores' atomic-reduction style).  Pure function of
// restrictAdaptiveMoments' own outputs, exactly like the Fortran
// evaluate_polarization_gate it mirrors.
__global__ void evaluateAdaptivePolarizationGate(GpuAdaptiveDeviceTopology topology,
                                                 GpuAdaptiveDeviceRuntime runtime,
                                                 real polarizationThreshold) {
   const std::size_t block = adaptiveThreadIndex();
   if(block >= topology.blocks) return;
   unsigned char unsafe = 0;
   // RCG-03 diagnostic: track the worst (minimum) ratio across every
   // channel/ensemble, mirroring the Fortran block_ratio output.  This
   // cannot short-circuit on the first unsafe hit the way the mask-only
   // computation could, or later (safer) channels would never be visited.
   real worst = real(1);
   for(std::size_t channel = 0; channel < topology.dynamicChannels; ++channel) {
      for(std::size_t ensemble = 0; ensemble < topology.ensembles; ++ensemble) {
         const std::size_t scalar = coarseScalarIndex(
            channel, block, ensemble, topology.dynamicChannels, topology.blocks);
         const real resultant[3] = {
            runtime.coarseMoment[3 * scalar],
            runtime.coarseMoment[3 * scalar + 1],
            runtime.coarseMoment[3 * scalar + 2]
         };
         const real length = normDevice(resultant);
         const real total = runtime.channelMomentSum[scalar];
         const real scale = total > real(1) ? total : real(1);
         const bool directionDefined = length > real(64) * epsilonDevice() * scale;
         if(!directionDefined) {
            unsafe = 1;
            worst = real(0);
            continue;
         }
         const real ratio = length / total;
         if(ratio < polarizationThreshold) unsafe = 1;
         if(ratio < worst) worst = ratio;
      }
   }
   runtime.polarizationUnsafeMask[block] = unsafe;
   if(runtime.polarizationRatio) runtime.polarizationRatio[block] = worst;
}

__global__ void proposeAdaptiveState(GpuAdaptiveDeviceTopology topology,
                                     GpuAdaptiveDeviceRuntime runtime,
                                     GpuAdaptiveSelectorPolicy policy,
                                     const unsigned char* hardMask) {
   const std::size_t block = adaptiveThreadIndex();
   if(block >= topology.blocks) return;
   int proposed = runtime.blockState[block];
   const real score =
      runtime.selectorScores[runtime.selectorCriteria * block];
   const bool dwellReady = runtime.stateAge[block] >= policy.minimumDwellUpdates;
   if(hardMask && hardMask[block]) proposed = fineState;
   else if(runtime.blockState[block] == coarseState &&
           score >= policy.refineThreshold && dwellReady)
      proposed = fineState;
   else if(runtime.blockState[block] != coarseState &&
           score <= policy.coarsenThreshold && dwellReady)
      proposed = coarseState;
   runtime.pendingState[block] = proposed;
}

// RCG-05E: monotonic pending-state invariant this kernel depends on.
// dilateAdaptiveState reads runtime.pendingState[source] for neighbouring
// blocks while (conceptually) deciding whether to promote its own target
// block to bufferState. runtime.pendingState itself is populated by the
// separate, prior proposeAdaptiveState kernel launch (gpuAdaptiveRuntime.cpp,
// proposeSelectorState's launch site) and by construction holds only
// coarseState or fineState values when this kernel begins -- proposeAdaptive
// State never writes bufferState (gpuAdaptiveRuntime.cpp:302-324). Two
// kernels on the same stream execute in launch order with the first fully
// complete before the second starts (ordinary CUDA/HIP stream semantics), so
// that precondition holds for every thread of this launch, not just some.
// This kernel itself never writes runtime.pendingState (see `dilatedState`
// below) -- the set of blocks holding fineState is therefore invariant
// (monotonically unchanged) for the entire lifetime of this launch, and is
// exactly why a thread's read of a neighbour's runtime.pendingState is safe
// to interpret as "genuinely fineState" regardless of what any other thread
// is concurrently doing. Any future change that made this kernel write
// runtime.pendingState directly (as it did before RCG-05E) would violate
// this invariant and reintroduce the read/write race RCG-05E removed; see
// the RCG-05E section of docs/RCG-05_GEOMETRY_OWNERSHIP_EVIDENCE.md for the
// sanitizer evidence and reasoning this comment summarizes.
//
// The actual write target is the separate `dilatedState` buffer
// (GpuAdaptiveRuntime::dilatedState_, gpuAdaptiveRuntime.hpp), not
// runtime.pendingState: each thread still writes only its own unique
// `target` index (one thread per block, so there is no writer-writer
// collision either), but that index now lives in an array no other thread
// in this launch ever reads, so the read (runtime.pendingState) and write
// (dilatedState) sets are disjoint by construction -- a real double buffer,
// not merely a value-domain argument for why a same-array race would have
// been benign. proposeSelectorState copies pendingState_ into dilatedState_
// before this launch (so untouched blocks keep their proposed state) and
// copies dilatedState_ back into pendingState_ after it, on the same stream,
// so every other consumer of deviceRuntime().pendingState still sees exactly
// the merged, fully-dilated result this kernel used to write in place.
__global__ void dilateAdaptiveState(GpuAdaptiveDeviceTopology topology,
                                    GpuAdaptiveDeviceRuntime runtime,
                                    GpuAdaptiveSelectorPolicy policy,
                                    int* dilatedState) {
   const std::size_t target = adaptiveThreadIndex();
   if(target >= topology.blocks ||
      runtime.pendingState[target] != coarseState) return;
   const int widthX = static_cast<int>(policy.bufferDilationBlocks[0]);
   const int widthY = static_cast<int>(policy.bufferDilationBlocks[1]);
   const int widthZ = static_cast<int>(policy.bufferDilationBlocks[2]);
   if(widthX == 0 && widthY == 0 && widthZ == 0) return;
   const int tx = topology.blockGridCoordinate[3 * target];
   const int ty = topology.blockGridCoordinate[3 * target + 1];
   const int tz = topology.blockGridCoordinate[3 * target + 2];
   for(int dz = -widthZ; dz <= widthZ; ++dz)
      for(int dy = -widthY; dy <= widthY; ++dy)
         for(int dx = -widthX; dx <= widthX; ++dx) {
            const int x = (tx + dx % topology.blockGrid[0] +
                           topology.blockGrid[0]) % topology.blockGrid[0];
            const int y = (ty + dy % topology.blockGrid[1] +
                           topology.blockGrid[1]) % topology.blockGrid[1];
            const int z = (tz + dz % topology.blockGrid[2] +
                           topology.blockGrid[2]) % topology.blockGrid[2];
            const std::size_t source = static_cast<std::size_t>(
               x + topology.blockGrid[0] *
               (y + topology.blockGrid[1] * z));
            if(runtime.pendingState[source] == fineState) {
               dilatedState[target] = bufferState;
               return;
            }
         }
}

__device__ inline void transverseBasis(const real mean[3], real e1[3], real e2[3]) {
   real reference[3] = {real(1), real(0), real(0)};
   if(fabs(mean[0]) > real(0.8)) {
      reference[0] = real(0);
      reference[1] = real(1);
   }
   crossDevice(mean, reference, e1);
   const real n = normDevice(e1);
   for(int xyz = 0; xyz < 3; ++xyz) e1[xyz] /= n;
   crossDevice(mean, e1, e2);
}

__device__ real coneLongitudinal(const GpuAdaptiveDeviceTopology& topology,
                                 const AdaptiveKernelDevice& kernels,
                                 std::size_t begin, std::size_t end, real amplitude) {
   real sum = real(0);
   for(std::size_t position = begin; position < end; ++position) {
      const std::size_t atom = static_cast<std::size_t>(topology.blockAtoms[position] - 1);
      const real tx = kernels.ghostDirection[3 * atom];
      const real ty = kernels.ghostDirection[3 * atom + 1];
      const real radius2 = amplitude * amplitude * (tx * tx + ty * ty);
      sum += kernels.atomMoment[atom] * sqrt(radius2 < real(1) ? real(1) - radius2 : real(0));
   }
   return sum;
}

// RCG-08 (F-09): publication used to walk every block on thread (0,0),
// serializing backup, reconstruction, rollback and the state/age/epoch update.
// It is now one thread per spatial block.
//
// Every block's work is independent by construction, so this needs no atomic
// and no inter-thread ordering:
//   * reads are own-block only -- coarseMoment/channelMomentSum/coarseDirection
//     at this block's scalars, blockState/pendingState/stateAge/transitionEpoch
//     at this block, and acceptedMask[block];
//   * writes to atomDirection, kernels.transitionBackup and
//     kernels.ghostDirection touch only atoms in this block's CSR range, and
//     validate() has already proved the CSR membership is a bijection (each
//     atom belongs to exactly one block), so no two threads can address the
//     same atom;
//   * writes to pendingState/blockState/stateAge/transitionEpoch are at this
//     thread's own block index.
// Determinism is unchanged: the reconstruction RNG is still seeded from the
// (globalSeed, block, channel, ensemble, epoch) tuple, so each block draws the
// identical stream it drew when the loop was serial, independently of the
// order the blocks now execute in.
__global__ void publishAdaptiveState(GpuAdaptiveDeviceTopology topology,
                                     GpuAdaptiveDeviceRuntime runtime,
                                     AdaptiveKernelDevice kernels,
                                     real* atomDirection,
                                     GpuAdaptiveReconstructionPolicy policy,
                                     const unsigned char* acceptedMask) {
   const std::size_t block = adaptiveThreadIndex();
   if(block >= topology.blocks) return;
   {
      if(acceptedMask && !acceptedMask[block]) {
         runtime.pendingState[block] = runtime.blockState[block];
         if(runtime.stateAge[block] != UINT_MAX)
            ++runtime.stateAge[block];
         return;
      }
      const int oldState = runtime.blockState[block];
      const int newState = runtime.pendingState[block];
      bool reconstructionOk = true;
      if(oldState != newState && oldState == coarseState && newState != coarseState) {
         const std::size_t begin = static_cast<std::size_t>(topology.blockAtomOffset[block]);
         const std::size_t end = static_cast<std::size_t>(topology.blockAtomOffset[block + 1]);
         for(std::size_t ensemble = 0; ensemble < topology.ensembles; ++ensemble)
            for(std::size_t position = begin; position < end; ++position) {
               const std::size_t atom =
                  static_cast<std::size_t>(topology.blockAtoms[position] - 1);
               for(int xyz = 0; xyz < 3; ++xyz)
                  kernels.transitionBackup[atomVectorIndex(
                     xyz, atom, ensemble, topology.atoms)] =
                     atomDirection[atomVectorIndex(
                        xyz, atom, ensemble, topology.atoms)];
            }
         for(std::size_t channel = 0; channel < topology.dynamicChannels; ++channel) {
            for(std::size_t ensemble = 0; ensemble < topology.ensembles; ++ensemble) {
               const std::size_t scalar = coarseScalarIndex(channel, block, ensemble,
                                                            topology.dynamicChannels,
                                                            topology.blocks);
               real resultant[3] = {
                  runtime.coarseMoment[3 * scalar],
                  runtime.coarseMoment[3 * scalar + 1],
                  runtime.coarseMoment[3 * scalar + 2]
               };
               const real requested = normDevice(resultant);
               const real total = runtime.channelMomentSum[scalar];
               const real scale = total > real(1) ? total : real(1);
               if(requested <= policy.resultantTolerance * scale ||
                  requested > total + policy.resultantTolerance * scale) {
                  reconstructionOk = false;
                  break;
               }
               if(policy.scheme == GpuAdaptiveReconstruction::Aligned &&
                  fabs(requested - total) > policy.resultantTolerance * scale) {
                  reconstructionOk = false;
                  break;
               }
               real mean[3] = {
                  resultant[0] / requested,
                  resultant[1] / requested,
                  resultant[2] / requested
               };
               if(policy.scheme == GpuAdaptiveReconstruction::Aligned) {
                  for(std::size_t position = begin; position < end; ++position) {
                     const std::size_t atom =
                        static_cast<std::size_t>(topology.blockAtoms[position] - 1);
                     if(topology.atomToDynamicChannel[atom] != static_cast<int>(channel + 1)) continue;
                     for(int xyz = 0; xyz < 3; ++xyz)
                        atomDirection[atomVectorIndex(xyz, atom, ensemble, topology.atoms)] =
                           mean[xyz];
                  }
               } else {
                  real e1[3], e2[3];
                  transverseBasis(mean, e1, e2);
                  std::uint64_t rng = adaptiveTupleSeed(policy.globalSeed, block, channel, ensemble,
                                                runtime.transitionEpoch[block]);
                  real sampledTotal = real(0), meanX = real(0), meanY = real(0);
                  for(std::size_t position = begin; position < end; ++position) {
                     const std::size_t atom =
                        static_cast<std::size_t>(topology.blockAtoms[position] - 1);
                     if(topology.atomToDynamicChannel[atom] != static_cast<int>(channel + 1)) continue;
                     const real angle = real(2) * piValue * adaptiveNextUniform(rng);
                     const real radius = real(0.5) + real(0.5) * adaptiveNextUniform(rng);
                     const real tx = radius * cos(angle);
                     const real ty = radius * sin(angle);
                     kernels.ghostDirection[3 * atom] = tx;
                     kernels.ghostDirection[3 * atom + 1] = ty;
                     sampledTotal += kernels.atomMoment[atom];
                     meanX += kernels.atomMoment[atom] * tx;
                     meanY += kernels.atomMoment[atom] * ty;
                  }
                  meanX /= sampledTotal;
                  meanY /= sampledTotal;
                  real maximumRadius = real(0);
                  for(std::size_t position = begin; position < end; ++position) {
                     const std::size_t atom =
                        static_cast<std::size_t>(topology.blockAtoms[position] - 1);
                     if(topology.atomToDynamicChannel[atom] != static_cast<int>(channel + 1)) continue;
                     real& tx = kernels.ghostDirection[3 * atom];
                     real& ty = kernels.ghostDirection[3 * atom + 1];
                     tx -= meanX;
                     ty -= meanY;
                     const real radius = sqrt(tx * tx + ty * ty);
                     if(radius > maximumRadius) maximumRadius = radius;
                  }
                  if(maximumRadius > real(0)) {
                     for(std::size_t position = begin; position < end; ++position) {
                        const std::size_t atom =
                           static_cast<std::size_t>(topology.blockAtoms[position] - 1);
                        if(topology.atomToDynamicChannel[atom] != static_cast<int>(channel + 1)) continue;
                        kernels.ghostDirection[3 * atom] /= maximumRadius;
                        kernels.ghostDirection[3 * atom + 1] /= maximumRadius;
                     }
                  }
                  real lower = real(0);
                  real upper = sin(policy.coneAngleRadians);
                  const real minimumResultant =
                     coneLongitudinal(topology, kernels, begin, end, upper);
                  if(requested < minimumResultant -
                     policy.resultantTolerance * scale) {
                     reconstructionOk = false;
                     break;
                  }
                  for(int iteration = 0; iteration < 100; ++iteration) {
                     const real amplitude = real(0.5) * (lower + upper);
                     const real longitudinal =
                        coneLongitudinal(topology, kernels, begin, end, amplitude);
                     if(longitudinal > requested) lower = amplitude;
                     else upper = amplitude;
                     if(fabs(longitudinal - requested) <=
                        real(0.1) * policy.resultantTolerance * scale) break;
                  }
                  const real amplitude = real(0.5) * (lower + upper);
                  for(std::size_t position = begin; position < end; ++position) {
                     const std::size_t atom =
                        static_cast<std::size_t>(topology.blockAtoms[position] - 1);
                     if(topology.atomToDynamicChannel[atom] != static_cast<int>(channel + 1)) continue;
                     const real tx = kernels.ghostDirection[3 * atom];
                     const real ty = kernels.ghostDirection[3 * atom + 1];
                     const real radius2 = amplitude * amplitude * (tx * tx + ty * ty);
                     const real longitudinal = sqrt(radius2 < real(1) ?
                                                    real(1) - radius2 : real(0));
                     for(int xyz = 0; xyz < 3; ++xyz)
                        atomDirection[atomVectorIndex(xyz, atom, ensemble, topology.atoms)] =
                           longitudinal * mean[xyz] + amplitude * tx * e1[xyz] +
                           amplitude * ty * e2[xyz];
                  }
               }
            }
            if(!reconstructionOk) break;
         }
         if(!reconstructionOk) {
            for(std::size_t ensemble = 0; ensemble < topology.ensembles; ++ensemble)
               for(std::size_t position = begin; position < end; ++position) {
                  const std::size_t atom =
                     static_cast<std::size_t>(topology.blockAtoms[position] - 1);
                  for(int xyz = 0; xyz < 3; ++xyz)
                     atomDirection[atomVectorIndex(
                        xyz, atom, ensemble, topology.atoms)] =
                        kernels.transitionBackup[atomVectorIndex(
                           xyz, atom, ensemble, topology.atoms)];
               }
            runtime.pendingState[block] = oldState;
            if(runtime.stateAge[block] != UINT_MAX) ++runtime.stateAge[block];
            return;
         }
      }
      if(oldState != newState) {
         runtime.blockState[block] = newState;
         runtime.stateAge[block] = 0;
         ++runtime.transitionEpoch[block];
      } else if(runtime.stateAge[block] != UINT_MAX) {
         ++runtime.stateAge[block];
      }
   }
}

__global__ void prolongateAdaptiveGhosts(GpuAdaptiveDeviceTopology topology,
                                         GpuAdaptiveDeviceRuntime runtime,
                                         AdaptiveKernelDevice kernels) {
   const std::size_t work = adaptiveThreadIndex();
   const std::size_t count = topology.atoms * topology.ensembles;
   if(work >= count) return;
   const std::size_t atom = work % topology.atoms;
   if(runtime.atomisticAtomMask[atom]) return;
   const std::size_t ensemble = work / topology.atoms;
   const int rawChannel = topology.atomToDynamicChannel[atom];
   if(rawChannel <= 0) return;
   const std::size_t channel = static_cast<std::size_t>(rawChannel - 1);
   real interpolated[3] = {real(0), real(0), real(0)};
   for(int corner = 0; corner < 8; ++corner) {
      const std::size_t block = static_cast<std::size_t>(
         kernels.projectionBlock[corner + 8 * atom] - 1);
      const real weight = kernels.projectionWeight[corner + 8 * atom];
      for(int xyz = 0; xyz < 3; ++xyz)
         interpolated[xyz] += weight *
            runtime.coarseDirection[coarseVectorIndex(
               xyz, channel, block, ensemble, topology.dynamicChannels,
               topology.blocks)];
   }
   const real length = normDevice(interpolated);
   kernels.projectionNorm[atom + topology.atoms * ensemble] = length;
   for(int xyz = 0; xyz < 3; ++xyz)
      kernels.ghostDirection[atomVectorIndex(xyz, atom, ensemble, topology.atoms)] =
         length > kernels.normalizationFloor ? interpolated[xyz] / length : real(0);
}

__global__ void commitAdaptiveGhosts(
   GpuAdaptiveDeviceTopology topology, GpuAdaptiveDeviceRuntime runtime,
   AdaptiveKernelDevice kernels, real* atomDirection,
   GpuAdaptiveReconstructionPolicy policy) {
   const std::size_t work = adaptiveThreadIndex();
   const std::size_t count = topology.atoms * topology.ensembles;
   if(work >= count) return;
   const std::size_t atom = work % topology.atoms;
   if(runtime.atomisticAtomMask[atom]) return;
   const std::size_t ensemble = work / topology.atoms;
   const std::size_t block =
      static_cast<std::size_t>(topology.atomToBlock[atom] - 1);
   const int rawChannel = topology.atomToDynamicChannel[atom];
   if(rawChannel <= 0) return;
   const std::size_t channel = static_cast<std::size_t>(rawChannel - 1);
   for(int xyz = 0; xyz < 3; ++xyz) {
      real value = kernels.ghostDirection[atomVectorIndex(
         xyz, atom, ensemble, topology.atoms)];
      if(policy.scheme == GpuAdaptiveReconstruction::Aligned)
         value = runtime.coarseDirection[coarseVectorIndex(
            xyz, channel, block, ensemble, topology.dynamicChannels,
            topology.blocks)];
      atomDirection[atomVectorIndex(
         xyz, atom, ensemble, topology.atoms)] = value;
   }
}

__device__ inline void effectiveAtomDirection(
   const GpuAdaptiveDeviceTopology& topology,
   const GpuAdaptiveDeviceRuntime& runtime,
   const AdaptiveKernelDevice& kernels,
   const real* atomDirection, std::size_t atom, std::size_t ensemble,
   real value[3]) {
   const real* source = runtime.atomisticAtomMask[atom] ?
                        atomDirection : kernels.ghostDirection;
   loadAtomVector(source, topology, atom, ensemble, value);
}

// RCG-08 (F-09): the atomistic Hamiltonian used to run entirely on thread
// (0,0) -- the clear pass, the unique-pair bond loop, the on-site anisotropy
// loop and the compact-list writeback, for every ensemble.  On the 8192-atom
// /24576-bond benchmark fixture that single lane accounted for 94.6% of the
// all-fine field-evaluation wall time.  It is now four parallel launches,
// ordered on the runtime stream so each sees the previous one complete.

__global__ void clearAdaptiveAtomistic(
   GpuAdaptiveDeviceTopology topology, AdaptiveKernelDevice kernels,
   real* atomField) {
   const std::size_t atomVectors = 3 * topology.atoms * topology.ensembles;
   const std::size_t index = adaptiveThreadIndex();
   if(index < atomVectors) {
      kernels.atomFieldScratch[index] = real(0);
      atomField[index] = real(0);
   }
   if(index < 2) kernels.energyTerms[index] = 0.0;
}

// One thread per (bond, ensemble).  The unique-pair ownership contract is
// unchanged -- every bond is still visited exactly once and still scatters
// the reaction field onto both of its endpoints -- but two bonds can now share
// an endpoint atom concurrently, so the endpoint accumulation must be atomic.
// This mirrors the choice RCG-07 accepted for the CPU analogue of the same
// loop (per-component atomic updates rather than an atom-sized per-thread
// reduction array, which would multiply an O(natoms) array by the thread
// count and undo RCG-06A).  The energy accumulates through the RCG-06B FP64
// accumulator, so the term stays FP64 regardless of build precision; its
// summation order becomes atomic-arrival order, exactly as the already
// accepted coarse energy terms (energyTerms[2..6]) have been since CG-10.
__global__ void evaluateAdaptiveAtomisticBonds(
   GpuAdaptiveDeviceTopology topology, GpuAdaptiveDeviceRuntime runtime,
   AdaptiveKernelDevice kernels, const real* atomDirection) {
   const std::size_t work = adaptiveThreadIndex();
   const std::size_t count = kernels.bonds * topology.ensembles;
   if(work >= count) return;
   const std::size_t bond = work % kernels.bonds;
   const std::size_t ensemble = work / kernels.bonds;
   const std::size_t atomI =
      static_cast<std::size_t>(kernels.bondAtom[2 * bond] - 1);
   const std::size_t atomJ =
      static_cast<std::size_t>(kernels.bondAtom[2 * bond + 1] - 1);
   if(!runtime.atomisticAtomMask[atomI] && !runtime.atomisticAtomMask[atomJ]) return;
   real si[3], sj[3], ksiJ[3] = {}, ktSi[3] = {};
   effectiveAtomDirection(topology, runtime, kernels, atomDirection, atomI, ensemble, si);
   effectiveAtomDirection(topology, runtime, kernels, atomDirection, atomJ, ensemble, sj);
   for(int row = 0; row < 3; ++row) {
      for(int column = 0; column < 3; ++column) {
         const real matrix = kernels.bondMatrix[row + 3 * (column + 3 * bond)];
         ksiJ[row] += matrix * sj[column];
         ktSi[column] += matrix * si[row];
      }
   }
   atomicAddEnergyTerm(&kernels.energyTerms[0],
                       -static_cast<double>(dotDevice(si, ksiJ)));
   for(int xyz = 0; xyz < 3; ++xyz) {
      atomicAdd(&kernels.atomFieldScratch[atomVectorIndex(
                   xyz, atomI, ensemble, topology.atoms)],
                ksiJ[xyz] / (kernels.magneticMomentSi * kernels.atomMoment[atomI]));
      atomicAdd(&kernels.atomFieldScratch[atomVectorIndex(
                   xyz, atomJ, ensemble, topology.atoms)],
                ktSi[xyz] / (kernels.magneticMomentSi * kernels.atomMoment[atomJ]));
   }
}

// One thread per (compact active-atom slot, ensemble).  activeAtomList holds
// each atom at most once, so (atom, ensemble) is unique per thread and the
// field update needs no atomic; only the shared energy term does.
__global__ void evaluateAdaptiveAtomisticOnsite(
   GpuAdaptiveDeviceTopology topology, GpuAdaptiveDeviceRuntime runtime,
   AdaptiveKernelDevice kernels, const real* atomDirection) {
   const std::size_t index = adaptiveThreadIndex();
   const std::size_t work = index % topology.atoms;
   const std::size_t ensemble = index / topology.atoms;
   if(ensemble >= topology.ensembles || work >= runtime.workCounts[0]) return;
   const std::size_t atom =
      static_cast<std::size_t>(runtime.activeAtomList[work] - 1);
   real direction[3];
   loadAtomVector(atomDirection, topology, atom, ensemble, direction);
   double onsite = 0.0;
   for(int axisIndex = 0;
       axisIndex < kernels.atomAnisotropyAxisCount[atom]; ++axisIndex) {
      real axis[3];
      for(int xyz = 0; xyz < 3; ++xyz)
         axis[xyz] = kernels.atomAnisotropyAxis[
            xyz + 3 * (axisIndex + 2 * atom)];
      const real c = dotDevice(direction, axis);
      const real k1 = kernels.atomAnisotropyK1[axisIndex + 2 * atom];
      const real k2 = kernels.atomAnisotropyK2[axisIndex + 2 * atom];
      onsite += static_cast<double>(
         k1 * c * c + k2 * (real(2) * c * c - c * c * c * c));
      const real derivative = real(2) * c *
         (k1 + real(2) * k2 * (real(1) - c * c));
      for(int xyz = 0; xyz < 3; ++xyz)
         kernels.atomFieldScratch[atomVectorIndex(
            xyz, atom, ensemble, topology.atoms)] -=
            derivative * axis[xyz] /
            (kernels.magneticMomentSi * kernels.atomMoment[atom]);
   }
   if(onsite != 0.0) atomicAddEnergyTerm(&kernels.energyTerms[1], onsite);
}

// The writeback traverses the compact active-atom list.  Coarse ghosts keep
// their reaction field private for the exact projection adjoint.
__global__ void writebackAdaptiveAtomistic(
   GpuAdaptiveDeviceTopology topology, GpuAdaptiveDeviceRuntime runtime,
   AdaptiveKernelDevice kernels, real* atomField) {
   const std::size_t index = adaptiveThreadIndex();
   const std::size_t work = index % topology.atoms;
   const std::size_t ensemble = index / topology.atoms;
   if(ensemble >= topology.ensembles || work >= runtime.workCounts[0]) return;
   const std::size_t atom =
      static_cast<std::size_t>(runtime.activeAtomList[work] - 1);
   for(int xyz = 0; xyz < 3; ++xyz)
      atomField[atomVectorIndex(xyz, atom, ensemble, topology.atoms)] =
         kernels.atomFieldScratch[atomVectorIndex(
            xyz, atom, ensemble, topology.atoms)];
}

__global__ void clearAdaptiveInterface(
   GpuAdaptiveDeviceTopology topology, GpuAdaptiveDeviceRuntime runtime,
   AdaptiveKernelDevice kernels) {
   const std::size_t coarseVectors =
      3 * topology.dynamicChannels * topology.blocks * topology.ensembles;
   const std::size_t index = adaptiveThreadIndex();
   if(index < coarseVectors) kernels.coarseFieldScratch[index] = real(0);
}

__global__ void restrictAdaptiveInterface(
   GpuAdaptiveDeviceTopology topology, GpuAdaptiveDeviceRuntime runtime,
   AdaptiveKernelDevice kernels) {
   const std::size_t work = adaptiveThreadIndex();
   const std::size_t count = topology.atoms * topology.ensembles;
   if(work >= count) return;
   const std::size_t atom = work % topology.atoms;
   if(runtime.atomisticAtomMask[atom]) return;
   const std::size_t ensemble = work / topology.atoms;
   const int rawChannel = topology.atomToDynamicChannel[atom];
   if(rawChannel <= 0) return;
   const std::size_t channel = static_cast<std::size_t>(rawChannel - 1);
   real direction[3], field[3], tangent[3];
   loadAtomVector(kernels.ghostDirection, topology, atom, ensemble, direction);
   loadAtomVector(kernels.atomFieldScratch, topology, atom, ensemble, field);
   const real projection = dotDevice(direction, field);
   for(int xyz = 0; xyz < 3; ++xyz)
      tangent[xyz] = field[xyz] - direction[xyz] * projection;
   const real rawNorm = kernels.projectionNorm[atom + topology.atoms * ensemble];
   for(int corner = 0; corner < 8; ++corner) {
      const std::size_t block = static_cast<std::size_t>(
         kernels.projectionBlock[corner + 8 * atom] - 1);
      const std::size_t scalar = coarseScalarIndex(
         channel, block, ensemble, topology.dynamicChannels, topology.blocks);
      const real coarseMoment = runtime.channelMomentSum[scalar];
      if(rawNorm <= kernels.normalizationFloor || coarseMoment <= real(0)) continue;
      const real factor = kernels.atomMoment[atom] *
         kernels.projectionWeight[corner + 8 * atom] /
         (coarseMoment * rawNorm);
      for(int xyz = 0; xyz < 3; ++xyz)
         atomicAdd(&kernels.coarseFieldScratch[3 * scalar + xyz],
                   factor * tangent[xyz]);
   }
}

__device__ inline void physicalGradient(
   const GpuAdaptiveDeviceTopology& topology, const AdaptiveKernelDevice& kernels,
   const real* coarseDirection, std::size_t block, std::size_t ensemble,
   int physical, real gradient[3]) {
   gradient[0] = gradient[1] = gradient[2] = real(0);
   for(int direction = 0; direction < 3; ++direction) {
      const real coefficient = kernels.inverseBlockTranspose[physical + 3 * direction];
      const std::size_t plus = plusBlock(topology, block, direction);
      for(int xyz = 0; xyz < 3; ++xyz)
         gradient[xyz] += coefficient *
            (coarseDirection[coarseVectorIndex(
                xyz, 0, plus, ensemble, topology.dynamicChannels, topology.blocks)] -
             coarseDirection[coarseVectorIndex(
                xyz, 0, block, ensemble, topology.dynamicChannels, topology.blocks)]);
   }
}

__device__ inline bool tensorTermOwned(
   const GpuAdaptiveDeviceTopology& topology,
   const GpuAdaptiveDeviceRuntime& runtime,
   const AdaptiveKernelDevice& kernels, std::size_t block, int p, int q) {
   if(!runtime.coarseBlockMask[block]) return false;
   for(int direction = 0; direction < 3; ++direction) {
      if(kernels.inverseBlockTranspose[p + 3 * direction] == real(0) &&
         (q < 0 || kernels.inverseBlockTranspose[q + 3 * direction] == real(0)))
         continue;
      if(!runtime.coarseBlockMask[plusBlock(topology, block, direction)]) return false;
   }
   return true;
}

__device__ inline void atomicAddDerivativeStencil(
   const GpuAdaptiveDeviceTopology& topology, const AdaptiveKernelDevice& kernels,
   real* derivative, std::size_t block, std::size_t ensemble,
   int physical, const real value[3], real scale) {
   for(int direction = 0; direction < 3; ++direction) {
      const real coefficient =
         scale * kernels.inverseBlockTranspose[physical + 3 * direction];
      const std::size_t plus = plusBlock(topology, block, direction);
      for(int xyz = 0; xyz < 3; ++xyz) {
         atomicAdd(&derivative[coarseVectorIndex(
                      xyz, 0, plus, ensemble, topology.dynamicChannels,
                      topology.blocks)],
                   coefficient * value[xyz]);
         atomicAdd(&derivative[coarseVectorIndex(
                      xyz, 0, block, ensemble, topology.dynamicChannels,
                      topology.blocks)],
                   -coefficient * value[xyz]);
      }
   }
}

__global__ void clearAdaptiveCoarse(
   GpuAdaptiveDeviceTopology topology, GpuAdaptiveDeviceRuntime runtime,
   AdaptiveKernelDevice kernels, real* coarseField) {
   const std::size_t vectorCount =
      3 * topology.dynamicChannels * topology.blocks * topology.ensembles;
   const std::size_t index = adaptiveThreadIndex();
   if(index < vectorCount) {
      runtime.coarseField[index] = real(0);
      coarseField[index] = real(0);
   }
   if(index < 6) kernels.energyTerms[index + 2] = 0.0;
}

__global__ void evaluateAdaptiveCoarseTensor(
   GpuAdaptiveDeviceTopology topology, GpuAdaptiveDeviceRuntime runtime,
   AdaptiveKernelDevice kernels) {
   const std::size_t index = adaptiveThreadIndex();
   const std::size_t work = index % topology.blocks;
   const std::size_t ensemble = index / topology.blocks;
   if(ensemble >= topology.ensembles || work >= runtime.workCounts[1]) return;
   const std::size_t block =
      static_cast<std::size_t>(runtime.activeBlockList[work] - 1);

   for(int p = 0; p < 3; ++p) {
      for(int q = 0; q < 3; ++q) {
         if(!tensorTermOwned(topology, runtime, kernels, block, p, q)) continue;
         real gradientP[3], gradientQ[3], derivative[3];
         physicalGradient(topology, kernels, runtime.coarseDirection,
                          block, ensemble, p, gradientP);
         physicalGradient(topology, kernels, runtime.coarseDirection,
                          block, ensemble, q, gradientQ);
         const real volume = topology.blockVolume[block];
         const real stiffness = kernels.exchangeStiffness[p + 3 * q];
         atomicAddEnergyTerm(&kernels.energyTerms[2], static_cast<double>(
                   volume * stiffness * dotDevice(gradientP, gradientQ)));
         for(int xyz = 0; xyz < 3; ++xyz)
            derivative[xyz] = volume * stiffness * gradientQ[xyz];
         atomicAddDerivativeStencil(topology, kernels, runtime.coarseField,
                                    block, ensemble, p, derivative, real(2));
      }
   }

   real direction[3];
   loadCoarseVector(runtime.coarseDirection, topology, 0, block,
                    ensemble, direction);
   for(int k = 0; k < 3; ++k) {
      real basis[3] = {real(0), real(0), real(0)};
      basis[k] = real(1);
      for(int p = 0; p < 3; ++p) {
         if(!tensorTermOwned(topology, runtime, kernels, block, p, -1)) continue;
         real gradient[3], crossMG[3], crossBasisGradient[3];
         real crossBasisDirection[3], derivative[3];
         physicalGradient(topology, kernels, runtime.coarseDirection,
                          block, ensemble, p, gradient);
         crossDevice(direction, gradient, crossMG);
         crossDevice(basis, gradient, crossBasisGradient);
         crossDevice(basis, direction, crossBasisDirection);
         const real volume = topology.blockVolume[block];
         const real spiral = kernels.spiralization[k + 3 * p];
         atomicAddEnergyTerm(&kernels.energyTerms[3],
                   static_cast<double>(volume * spiral * crossMG[k]));
         for(int xyz = 0; xyz < 3; ++xyz) {
            atomicAdd(&runtime.coarseField[coarseVectorIndex(
                         xyz, 0, block, ensemble, topology.dynamicChannels,
                         topology.blocks)],
                      -volume * spiral * crossBasisGradient[xyz]);
            derivative[xyz] = volume * spiral * crossBasisDirection[xyz];
         }
         atomicAddDerivativeStencil(topology, kernels, runtime.coarseField,
                                    block, ensemble, p, derivative, real(1));
      }
   }
}

__global__ void finalizeAdaptiveCoarseLocal(
   GpuAdaptiveDeviceTopology topology, GpuAdaptiveDeviceRuntime runtime,
   AdaptiveKernelDevice kernels, const real* externalField) {
   const std::size_t index = adaptiveThreadIndex();
   const std::size_t work = index % topology.blocks;
   const std::size_t ensemble = index / topology.blocks;
   if(ensemble >= topology.ensembles || work >= runtime.workCounts[1]) return;
   const std::size_t block =
      static_cast<std::size_t>(runtime.activeBlockList[work] - 1);
   real direction[3];
   loadCoarseVector(runtime.coarseDirection, topology, 0, block,
                    ensemble, direction);
   for(int axisIndex = 0;
       axisIndex < kernels.anisotropyAxisCount[block]; ++axisIndex) {
      real axis[3];
      for(int xyz = 0; xyz < 3; ++xyz)
         axis[xyz] = kernels.anisotropyAxis[
            xyz + 3 * (axisIndex + 2 * block)];
      const real c = dotDevice(direction, axis);
      const real k1 = kernels.anisotropyK1[axisIndex + 2 * block];
      const real k2 = kernels.anisotropyK2[axisIndex + 2 * block];
      const real volume = topology.blockVolume[block];
      atomicAddEnergyTerm(&kernels.energyTerms[4], static_cast<double>(volume *
         (k1 * c * c + real(2) * k2 * c * c -
          k2 * c * c * c * c)));
      const real derivative = volume * real(2) * c *
         (k1 + real(2) * k2 * (real(1) - c * c));
      for(int xyz = 0; xyz < 3; ++xyz)
         runtime.coarseField[coarseVectorIndex(
            xyz, 0, block, ensemble, topology.dynamicChannels,
            topology.blocks)] += derivative * axis[xyz];
   }
   const std::size_t scalar = coarseScalarIndex(
      0, block, ensemble, topology.dynamicChannels, topology.blocks);
   const real moment = runtime.channelMomentSum[scalar];
   for(int xyz = 0; xyz < 3; ++xyz) {
      const std::size_t vector = 3 * scalar + xyz;
      runtime.coarseField[vector] *=
         -real(1) / (kernels.magneticMomentSi * moment);
      if(externalField) runtime.coarseField[vector] += externalField[vector];
   }
   if(externalField) {
      real external[3];
      for(int xyz = 0; xyz < 3; ++xyz)
         external[xyz] = externalField[3 * scalar + xyz];
      atomicAddEnergyTerm(&kernels.energyTerms[5],
                -static_cast<double>(kernels.magneticMomentSi * moment *
                 dotDevice(external, direction)));
   }
}

__global__ void addAdaptiveDipole(
   GpuAdaptiveDeviceTopology topology, GpuAdaptiveDeviceRuntime runtime,
   AdaptiveKernelDevice kernels, const real* atomDirection,
   const real* fftDipoleField, real* atomField) {
   const std::size_t index = adaptiveThreadIndex();
   const std::size_t count = topology.blocks * topology.ensembles;
   if(index >= count) return;
   const std::size_t block = index % topology.blocks;
   const std::size_t ensemble = index / topology.blocks;
   const std::size_t scalar = coarseScalarIndex(
      0, block, ensemble, topology.dynamicChannels, topology.blocks);
   real direction[3] = {}, dipole[3], source[3] = {};
   for(int xyz = 0; xyz < 3; ++xyz)
      dipole[xyz] = fftDipoleField[3 * scalar + xyz];
   if(runtime.coarseBlockMask[block]) {
      loadCoarseVector(runtime.coarseDirection, topology, 0, block,
                       ensemble, direction);
      for(int xyz = 0; xyz < 3; ++xyz)
         source[xyz] = runtime.channelMomentSum[scalar] * direction[xyz];
   } else {
      const int begin = topology.blockAtomOffset[block];
      const int end = topology.blockAtomOffset[block + 1];
      for(int position = begin; position < end; ++position) {
         const std::size_t atom =
            static_cast<std::size_t>(topology.blockAtoms[position] - 1);
         real atomVector[3];
         loadAtomVector(atomDirection, topology, atom, ensemble, atomVector);
         for(int xyz = 0; xyz < 3; ++xyz) {
            source[xyz] += kernels.atomMoment[atom] * atomVector[xyz];
            atomField[atomVectorIndex(xyz, atom, ensemble, topology.atoms)] +=
               dipole[xyz];
         }
      }
   }
   atomicAddEnergyTerm(&kernels.energyTerms[6],
             -static_cast<double>(real(0.5) * kernels.magneticMomentSi *
              dotDevice(dipole, source)));
   if(runtime.coarseBlockMask[block]) {
      for(int xyz = 0; xyz < 3; ++xyz)
         runtime.coarseField[3 * scalar + xyz] += dipole[xyz];
   }
}

__device__ inline void loadBasisResolvedFftField(
   const GpuAdaptiveUniformFftField& fft, std::size_t block,
   unsigned int basis, unsigned int ensemble, real field[3]) {
   const std::size_t q1 = block % fft.activeN1;
   const std::size_t q2 = (block / fft.activeN1) % fft.activeN2;
   const std::size_t q3 = block / (fft.activeN1 * fft.activeN2);
   const std::size_t fftCell =
      q1 + fft.fftN1 * (q2 + fft.fftN2 * q3);
   for(int xyz = 0; xyz < 3; ++xyz) {
      const std::size_t fieldIndex =
         fftCell + fft.fftCells *
            (static_cast<std::size_t>(xyz) +
             3 * (basis + fft.basis * ensemble));
      field[xyz] = fft.prefactorT * fft.paddedField[fieldIndex];
   }
}

__global__ void addAdaptiveBasisResolvedDipole(
   GpuAdaptiveDeviceTopology topology, GpuAdaptiveDeviceRuntime runtime,
   AdaptiveKernelDevice kernels, const real* atomDirection,
   GpuAdaptiveUniformFftField fft, real* atomField) {
   const std::size_t index = adaptiveThreadIndex();
   const std::size_t count = topology.blocks * topology.ensembles;
   if(index >= count) return;
   const std::size_t block = index % topology.blocks;
   const std::size_t ensemble = index / topology.blocks;
   const std::size_t scalar = coarseScalarIndex(
      0, block, ensemble, topology.dynamicChannels, topology.blocks);
   real coarseDirection[3] = {};
   if(runtime.coarseBlockMask[block])
      loadCoarseVector(runtime.coarseDirection, topology, 0, block,
                       ensemble, coarseDirection);
   real coarseWeightedField[3] = {};
   // RCG-06B (F-11): this per-thread partial sum feeds the FP64 global
   // energy accumulator below; keep it double even though the field terms
   // it is built from stay `real`.
   double dipoleEnergy = 0.0;
   const int begin = topology.blockAtomOffset[block];
   const int end = topology.blockAtomOffset[block + 1];
   for(int position = begin; position < end; ++position) {
      const std::size_t atom =
         static_cast<std::size_t>(topology.blockAtoms[position] - 1);
      const int oneBasedFftChannel = topology.atomToFftChannel[atom];
      if(oneBasedFftChannel <= 0 ||
         static_cast<unsigned int>(oneBasedFftChannel) > fft.basis) continue;
      real field[3], sourceDirection[3];
      loadBasisResolvedFftField(
         fft, block, static_cast<unsigned int>(oneBasedFftChannel - 1),
         static_cast<unsigned int>(ensemble), field);
      if(runtime.coarseBlockMask[block]) {
         for(int xyz = 0; xyz < 3; ++xyz)
            sourceDirection[xyz] = coarseDirection[xyz];
      } else {
         loadAtomVector(atomDirection, topology, atom, ensemble,
                        sourceDirection);
         for(int xyz = 0; xyz < 3; ++xyz)
            atomField[atomVectorIndex(
               xyz, atom, ensemble, topology.atoms)] += field[xyz];
      }
      const real moment = kernels.atomMoment[atom];
      dipoleEnergy -= static_cast<double>(real(0.5) * kernels.magneticMomentSi *
                      moment * dotDevice(field, sourceDirection));
      for(int xyz = 0; xyz < 3; ++xyz)
         coarseWeightedField[xyz] += moment * field[xyz];
   }
   atomicAddEnergyTerm(&kernels.energyTerms[6], dipoleEnergy);
   if(runtime.coarseBlockMask[block]) {
      const real inverseMoment = real(1) / runtime.channelMomentSum[scalar];
      for(int xyz = 0; xyz < 3; ++xyz)
         runtime.coarseField[3 * scalar + xyz] +=
            inverseMoment * coarseWeightedField[xyz];
   }
}

__global__ void writeAdaptiveCoarse(
   GpuAdaptiveDeviceTopology topology, GpuAdaptiveDeviceRuntime runtime,
   AdaptiveKernelDevice kernels, real* coarseField) {
   const std::size_t index = adaptiveThreadIndex();
   const std::size_t work = index % topology.blocks;
   const std::size_t ensemble = index / topology.blocks;
   if(ensemble >= topology.ensembles || work >= runtime.workCounts[1]) return;
   const std::size_t block =
      static_cast<std::size_t>(runtime.activeBlockList[work] - 1);
   const std::size_t scalar = coarseScalarIndex(
      0, block, ensemble, topology.dynamicChannels, topology.blocks);
   for(int xyz = 0; xyz < 3; ++xyz) {
      const std::size_t vector = 3 * scalar + xyz;
      coarseField[vector] = runtime.coarseField[vector] +
                            kernels.coarseFieldScratch[vector];
   }
}

// RCG-08 (F-09) documented exception: this kernel is a genuine O(1) reduction
// -- seven FP64 loads, six adds, one store -- of accumulators that every other
// energy-producing kernel has already reduced in parallel.  It is launched
// <<<1,1>>> because there is exactly one element of work, not because O(N)
// work has been serialized onto one lane, which is what F-09 objected to and
// what the other five kernels in that finding actually did.  It is left
// single-threaded deliberately; parallelizing a seven-term sum would add a
// launch and a reduction buffer to save nanoseconds.  The guard is written
// against the shared thread-index helper so it stays correct if the launch
// geometry is ever widened.
__global__ void finalizeAdaptiveEnergy(AdaptiveKernelDevice kernels) {
   if(adaptiveThreadIndex() != 0) return;
   kernels.energyTerms[7] = kernels.energyTerms[0] + kernels.energyTerms[1] +
      kernels.energyTerms[2] + kernels.energyTerms[3] +
      kernels.energyTerms[4] + kernels.energyTerms[5] +
      kernels.energyTerms[6];
}

__device__ inline void llgRhs(const real direction[3], const real field[3],
                              real gamma, real damping, real rhs[3]) {
   real first[3], second[3];
   crossDevice(direction, field, first);
   crossDevice(direction, first, second);
   const real prefactor = -gamma / (real(1) + damping * damping);
   for(int xyz = 0; xyz < 3; ++xyz)
      rhs[xyz] = prefactor * (first[xyz] + damping * second[xyz]);
}

// RCG-08 (F-09): both Heun stages ran on thread (0,0), including the two
// full-state save copies at the top of the predictor.  They are now one thread
// per compact-list slot per ensemble, split into an atom kernel and a block
// kernel so each launch has a single uniform work definition.  The saves moved
// out of the kernel entirely and are issued as stream-ordered device-to-device
// copies by integrateHeun().
//
// Each thread owns one (atom, ensemble) or (block, ensemble) entirely: the
// compact lists contain each atom/block at most once, so there is no atomic,
// no reduction, and no cross-thread read of anything this launch writes.  The
// accepted Heun update, its normalization, and the compact-list work
// definition are unchanged, so results are bitwise identical to the serial
// reference.

__global__ void predictorAdaptiveHeunAtoms(
   GpuAdaptiveDeviceTopology topology, GpuAdaptiveDeviceRuntime runtime,
   AdaptiveKernelDevice kernels, real dt, real* atomDirection,
   const real* savedAtom, const real* atomField) {
   const std::size_t index = adaptiveThreadIndex();
   const std::size_t work = index % topology.atoms;
   const std::size_t ensemble = index / topology.atoms;
   if(ensemble >= topology.ensembles || work >= runtime.workCounts[0]) return;
   const std::size_t atom =
      static_cast<std::size_t>(runtime.activeAtomList[work] - 1);
   real direction[3], field[3], rhs[3], candidate[3];
   loadAtomVector(savedAtom, topology, atom, ensemble, direction);
   loadAtomVector(atomField, topology, atom, ensemble, field);
   llgRhs(direction, field, kernels.gammaPerTs, kernels.damping, rhs);
   for(int xyz = 0; xyz < 3; ++xyz) candidate[xyz] = direction[xyz] + dt * rhs[xyz];
   const real n = normDevice(candidate);
   for(int xyz = 0; xyz < 3; ++xyz)
      atomDirection[atomVectorIndex(xyz, atom, ensemble, topology.atoms)] =
         candidate[xyz] / n;
}

__global__ void predictorAdaptiveHeunBlocks(
   GpuAdaptiveDeviceTopology topology, GpuAdaptiveDeviceRuntime runtime,
   AdaptiveKernelDevice kernels, real dt, const real* savedCoarse,
   const real* coarseField) {
   const std::size_t index = adaptiveThreadIndex();
   const std::size_t work = index % topology.blocks;
   const std::size_t ensemble = index / topology.blocks;
   if(ensemble >= topology.ensembles || work >= runtime.workCounts[1]) return;
   const std::size_t block =
      static_cast<std::size_t>(runtime.activeBlockList[work] - 1);
   real direction[3], field[3], rhs[3], candidate[3];
   loadCoarseVector(savedCoarse, topology, 0, block, ensemble, direction);
   loadCoarseVector(coarseField, topology, 0, block, ensemble, field);
   llgRhs(direction, field, kernels.gammaPerTs, kernels.damping, rhs);
   for(int xyz = 0; xyz < 3; ++xyz) candidate[xyz] = direction[xyz] + dt * rhs[xyz];
   const real n = normDevice(candidate);
   for(int xyz = 0; xyz < 3; ++xyz)
      runtime.coarseDirection[coarseVectorIndex(
         xyz, 0, block, ensemble, topology.dynamicChannels,
         topology.blocks)] = candidate[xyz] / n;
}

__global__ void correctorAdaptiveHeunAtoms(
   GpuAdaptiveDeviceTopology topology, GpuAdaptiveDeviceRuntime runtime,
   AdaptiveKernelDevice kernels, real dt, real* atomDirection,
   const real* savedAtom, const real* initialAtomField,
   const real* predictorAtomField) {
   const std::size_t index = adaptiveThreadIndex();
   const std::size_t work = index % topology.atoms;
   const std::size_t ensemble = index / topology.atoms;
   if(ensemble >= topology.ensembles || work >= runtime.workCounts[0]) return;
   const std::size_t atom =
      static_cast<std::size_t>(runtime.activeAtomList[work] - 1);
   real initial[3], predictor[3], field0[3], field1[3], rhs0[3], rhs1[3], value[3];
   loadAtomVector(savedAtom, topology, atom, ensemble, initial);
   loadAtomVector(atomDirection, topology, atom, ensemble, predictor);
   loadAtomVector(initialAtomField, topology, atom, ensemble, field0);
   loadAtomVector(predictorAtomField, topology, atom, ensemble, field1);
   llgRhs(initial, field0, kernels.gammaPerTs, kernels.damping, rhs0);
   llgRhs(predictor, field1, kernels.gammaPerTs, kernels.damping, rhs1);
   for(int xyz = 0; xyz < 3; ++xyz)
      value[xyz] = initial[xyz] + real(0.5) * dt * (rhs0[xyz] + rhs1[xyz]);
   const real n = normDevice(value);
   for(int xyz = 0; xyz < 3; ++xyz)
      atomDirection[atomVectorIndex(xyz, atom, ensemble, topology.atoms)] =
         value[xyz] / n;
}

__global__ void correctorAdaptiveHeunBlocks(
   GpuAdaptiveDeviceTopology topology, GpuAdaptiveDeviceRuntime runtime,
   AdaptiveKernelDevice kernels, real dt, const real* savedCoarse,
   const real* initialCoarseField, const real* predictorCoarseField) {
   const std::size_t index = adaptiveThreadIndex();
   const std::size_t work = index % topology.blocks;
   const std::size_t ensemble = index / topology.blocks;
   if(ensemble >= topology.ensembles || work >= runtime.workCounts[1]) return;
   const std::size_t block =
      static_cast<std::size_t>(runtime.activeBlockList[work] - 1);
   real initial[3], predictor[3], field0[3], field1[3], rhs0[3], rhs1[3], value[3];
   loadCoarseVector(savedCoarse, topology, 0, block, ensemble, initial);
   loadCoarseVector(runtime.coarseDirection, topology, 0, block, ensemble, predictor);
   loadCoarseVector(initialCoarseField, topology, 0, block, ensemble, field0);
   loadCoarseVector(predictorCoarseField, topology, 0, block, ensemble, field1);
   llgRhs(initial, field0, kernels.gammaPerTs, kernels.damping, rhs0);
   llgRhs(predictor, field1, kernels.gammaPerTs, kernels.damping, rhs1);
   for(int xyz = 0; xyz < 3; ++xyz)
      value[xyz] = initial[xyz] + real(0.5) * dt * (rhs0[xyz] + rhs1[xyz]);
   const real n = normDevice(value);
   for(int xyz = 0; xyz < 3; ++xyz)
      runtime.coarseDirection[coarseVectorIndex(
         xyz, 0, block, ensemble, topology.dynamicChannels,
         topology.blocks)] = value[xyz] / n;
}

bool checkedAdd(std::size_t& total, std::size_t count, std::size_t elementBytes) {
   if(count > std::numeric_limits<std::size_t>::max() / elementBytes) return false;
   const std::size_t bytes = count * elementBytes;
   if(bytes > std::numeric_limits<std::size_t>::max() - total) return false;
   total += bytes;
   return true;
}

bool checkedProduct(std::initializer_list<std::size_t> factors, std::size_t& result) {
   result = 1;
   for(const std::size_t factor : factors) {
      if(factor != 0 && result > std::numeric_limits<std::size_t>::max() / factor) return false;
      result *= factor;
   }
   return true;
}

template <typename T>
bool required(const T* pointer) {
   return pointer != nullptr;
}

template <typename T>
void uploadNative(GpuTensor<T, 1>& destination, const T* source, std::size_t count,
                  GPU_STREAM_T stream) {
   Tensor<T, 1> view(const_cast<T*>(source), static_cast<index_t>(count));
   destination.copy_async(view, stream);
}

void uploadReal(GpuTensor<real, 1>& destination, const double* source, std::size_t count,
                GPU_STREAM_T stream, std::vector<std::vector<real>>& staging) {
   if(!source) {
      destination.zeros_async(stream);
      return;
   }
#ifndef SINGLE_PREC
      (void)staging;
      Tensor<real, 1> view(const_cast<real*>(source), static_cast<index_t>(count));
      destination.copy_async(view, stream);
#else
      staging.emplace_back(count);
      auto& converted = staging.back();
      std::transform(source, source + count, converted.begin(),
                     [](double value) { return static_cast<real>(value); });
      Tensor<real, 1> view(converted.data(), static_cast<index_t>(count));
      destination.copy_async(view, stream);
#endif
}

void freeIfAllocated(GpuTensor<int, 1>& tensor) {
   if(!tensor.empty()) tensor.Free();
}
void freeIfAllocated(GpuTensor<unsigned int, 1>& tensor) {
   if(!tensor.empty()) tensor.Free();
}
void freeIfAllocated(GpuTensor<unsigned char, 1>& tensor) {
   if(!tensor.empty()) tensor.Free();
}
void freeIfAllocated(GpuTensor<real, 1>& tensor) {
   if(!tensor.empty()) tensor.Free();
}
#ifdef SINGLE_PREC
// RCG-06B (F-11): energyTerms_ is GpuTensor<double, 1> unconditionally
// (independent of `real`'s build precision). In a DOUBLE_PREC build `real`
// already is `double`, so the overload above already covers it; a second
// overload would be a duplicate-signature redefinition, not a distinct
// overload. Only SINGLE_PREC builds need this one.
void freeIfAllocated(GpuTensor<double, 1>& tensor) {
   if(!tensor.empty()) tensor.Free();
}
#endif

__global__ void initializeAdaptiveWorkScan(
   GpuAdaptiveDeviceTopology topology, GpuAdaptiveDeviceRuntime runtime,
   unsigned int* scan, std::size_t scanItems) {
   const std::size_t index = adaptiveThreadIndex();
   if(index >= scanItems) return;
   unsigned int activeAtom = 0;
   unsigned int activeBlock = 0;
   unsigned int interfaceAtom = 0;
   if(index < topology.blocks) {
      const std::size_t block = index;
      const int state = runtime.blockState[block];
      const bool atomistic = state != coarseState;
      const bool coarse = state == coarseState;
      runtime.atomisticBlockMask[block] = static_cast<unsigned char>(atomistic);
      runtime.coarseBlockMask[block] = static_cast<unsigned char>(coarse);
      activeBlock = static_cast<unsigned int>(coarse);
   }
   if(index < topology.atoms) {
      const std::size_t atom = index;
      const int block = topology.atomToBlock[atom] - 1;
      const int state = runtime.blockState[block];
      const bool atomistic = state != coarseState;
      const bool interface = state == bufferState;
      runtime.atomisticAtomMask[atom] = static_cast<unsigned char>(atomistic);
      runtime.interfaceAtomMask[atom] = static_cast<unsigned char>(interface);
      activeAtom = static_cast<unsigned int>(atomistic);
      interfaceAtom = static_cast<unsigned int>(interface);
   }
   scan[index] = activeAtom;
   scan[scanItems + index] = activeBlock;
   scan[2 * scanItems + index] = interfaceAtom;
}

// RCG-08 (F-12): the compaction scan used to be a global Hillis--Steele
// sweep -- ceil(log2(N)) kernel launches, each reading and writing all 3N
// flags, for O(N log N) element work with an N-dependent launch count.
//
// It is replaced by the standard three-phase hierarchical scan below, which
// is linear-work and portable (plain shared memory, no CUB/rocPRIM, so CUDA
// and HIP compile the identical algorithm rather than two different vendor
// primitives that would have to be validated separately -- and HIP has no
// device here to validate on).
//
//   1. scanAdaptiveTiles     -- each thread block inclusively scans one
//                               adaptiveThreads-sized tile of one component in
//                               shared memory and publishes that tile's total.
//   2. scanAdaptiveTiles     -- applied again to the (much shorter) array of
//                               tile totals, recursively, until one tile
//                               covers a level.
//   3. addAdaptiveTileOffsets - adds each level's scanned tile totals back
//                               down into the level below.
//
// Work accounting: the tile-local scan is Kogge--Stone, so it costs
// T*log2(T) operations per T=adaptiveThreads elements.  T is a compile-time
// constant, so that is a fixed factor (log2(256)=8) per element, not a growing
// one: total element work is O(N), and the launch count is O(log_T N) -- three
// launches for the 2048-block benchmark fixture and five up to ~16M items,
// against the eleven the old sweep needed at 2048 alone.  The tile-local
// factor of 8 is stated rather than hidden: a work-efficient Blelloch tile
// scan would reduce it to ~2, but the scan is already far from the hot path
// and the simpler formulation is easier to keep identical across backends.
//
// `values` is scanned in place when it is a tile-sum level, and out-of-place
// for level 0; both are safe because every thread block reads and writes only
// its own tile.
__global__ void scanAdaptiveTiles(
   const unsigned int* input, unsigned int* output, unsigned int* tileTotals,
   std::size_t itemsPerComponent, std::size_t tilesPerComponent) {
   __shared__ unsigned int shared[adaptiveThreads];
   const std::size_t component =
      static_cast<std::size_t>(blockIdx.x) / tilesPerComponent;
   const std::size_t tile =
      static_cast<std::size_t>(blockIdx.x) % tilesPerComponent;
   const std::size_t base = component * itemsPerComponent + tile * adaptiveThreads;
   const std::size_t local = threadIdx.x;
   const std::size_t index = tile * adaptiveThreads + local;
   shared[local] = index < itemsPerComponent ? input[base + local] : 0u;
   __syncthreads();
   for(unsigned int offset = 1; offset < adaptiveThreads; offset <<= 1) {
      unsigned int addend = 0u;
      if(local >= offset) addend = shared[local - offset];
      __syncthreads();
      if(local >= offset) shared[local] += addend;
      __syncthreads();
   }
   if(index < itemsPerComponent) output[base + local] = shared[local];
   if(local == adaptiveThreads - 1 && tileTotals)
      tileTotals[component * tilesPerComponent + tile] = shared[local];
}

__global__ void addAdaptiveTileOffsets(
   unsigned int* values, const unsigned int* scannedTileTotals,
   std::size_t itemsPerComponent, std::size_t tilesPerComponent) {
   const std::size_t flat = adaptiveThreadIndex();
   const std::size_t count = 3 * itemsPerComponent;
   if(flat >= count) return;
   const std::size_t component = flat / itemsPerComponent;
   const std::size_t index = flat - component * itemsPerComponent;
   const std::size_t tile = index / adaptiveThreads;
   if(tile == 0) return;
   values[flat] +=
      scannedTileTotals[component * tilesPerComponent + tile - 1];
}

__global__ void scatterAdaptiveWork(
   GpuAdaptiveDeviceTopology topology, GpuAdaptiveDeviceRuntime runtime,
   const unsigned int* inclusiveScan, std::size_t scanItems) {
   const std::size_t index = adaptiveThreadIndex();
   if(index >= scanItems) return;
   if(index < topology.atoms && runtime.atomisticAtomMask[index]) {
      const unsigned int position = inclusiveScan[index] - 1;
      runtime.activeAtomList[position] = static_cast<int>(index + 1);
   }
   if(index < topology.blocks && runtime.coarseBlockMask[index]) {
      const unsigned int position = inclusiveScan[scanItems + index] - 1;
      runtime.activeBlockList[position] = static_cast<int>(index + 1);
   }
   if(index < topology.atoms && runtime.interfaceAtomMask[index]) {
      const unsigned int position = inclusiveScan[2 * scanItems + index] - 1;
      runtime.interfaceAtomList[position] = static_cast<int>(index + 1);
   }
   if(index == 0) {
      runtime.workCounts[0] = inclusiveScan[scanItems - 1];
      runtime.workCounts[1] = inclusiveScan[2 * scanItems - 1];
      runtime.workCounts[2] = inclusiveScan[3 * scanItems - 1];
   }
}

} // namespace

GpuAdaptiveRuntime::~GpuAdaptiveRuntime() {
   if(ready_ || streamCreated_) {
      try {
         release();
      } catch(...) {
         // Destructors must not propagate driver-shutdown errors.
      }
   }
}

bool GpuAdaptiveRuntime::validate(const GpuAdaptiveTopologyInput& t,
                                  const GpuAdaptiveRuntimeInput& r,
                                  std::size_t expectedAtoms,
                                  std::size_t expectedEnsembles,
                                  std::string& diagnostic) {
   diagnostic.clear();
   if(t.geometryMode != regularReplicatedCell) {
      diagnostic = "GPU adaptive runtime supports REGULAR_REPLICATED_CELL topology only";
      return false;
   }
   if(t.atoms == 0 || t.blocks == 0 || t.basis == 0 ||
      t.fftChannelsPerBlock == 0 || t.fftGridChannels == 0 ||
      t.dynamicChannels == 0 || t.ensembles == 0) {
      diagnostic = "GPU adaptive topology counts must all be positive";
      return false;
   }
   if(t.atoms != expectedAtoms || t.ensembles != expectedEnsembles) {
      diagnostic = "Fortran and GPU adaptive topology atom/ensemble counts do not match";
      return false;
   }
   if(t.atoms > static_cast<std::size_t>(std::numeric_limits<int>::max()) ||
      t.blocks > static_cast<std::size_t>(std::numeric_limits<int>::max())) {
      diagnostic = "GPU adaptive topology exceeds the compact-list integer range";
      return false;
   }
   if(!required(t.repetitionShape) || !required(t.blockShape) ||
      !required(t.blockGrid) || !required(t.cellVectors) ||
      !required(t.blockVectors) || !required(t.atomToBlock) || !required(t.atomToBasis) ||
      !required(t.atomToDynamicChannel) || !required(t.atomToFftChannel) ||
      !required(t.atomToFftGridIndex) || !required(t.basisToDynamicChannel) ||
      !required(t.basisToFftChannel) || !required(t.blockAtomCount) ||
      !required(t.blockAtomOffset) || !required(t.blockAtoms) ||
      !required(t.blockGridCoordinate) || !required(t.blockBasisPopulation) ||
      !required(t.blockFftChannelPopulation) ||
      !required(t.blockDynamicChannelPopulation) || !required(t.blockCenter) ||
      !required(t.blockVolume) || !required(r.blockState)) {
      diagnostic = "GPU adaptive topology staging is missing a mandatory array";
      return false;
   }
   std::size_t expectedBlocks = 1;
   std::size_t expectedAtomsFromShape = t.basis;
   for(int d = 0; d < 3; ++d) {
      if(t.repetitionShape[d] <= 0 || t.blockShape[d] <= 0 || t.blockGrid[d] <= 0 ||
         static_cast<long long>(t.blockShape[d]) * t.blockGrid[d] != t.repetitionShape[d]) {
         diagnostic = "GPU adaptive repetition, block-shape, and block-grid counts disagree";
         return false;
      }
      if(expectedBlocks > std::numeric_limits<std::size_t>::max() /
                             static_cast<std::size_t>(t.blockGrid[d]) ||
         expectedAtomsFromShape > std::numeric_limits<std::size_t>::max() /
                                     static_cast<std::size_t>(t.repetitionShape[d])) {
         diagnostic = "GPU adaptive topology shape counts overflow";
         return false;
      }
      expectedBlocks *= static_cast<std::size_t>(t.blockGrid[d]);
      expectedAtomsFromShape *= static_cast<std::size_t>(t.repetitionShape[d]);
   }
   std::size_t expectedFftGrid = 0;
   if(!checkedProduct({t.fftChannelsPerBlock, t.blocks}, expectedFftGrid) ||
      expectedBlocks != t.blocks || expectedAtomsFromShape != t.atoms ||
      t.fftGridChannels != expectedFftGrid) {
      diagnostic = "GPU adaptive scalar topology counts are internally inconsistent";
      return false;
   }
   for(int i = 0; i < 9; ++i) {
      if(!std::isfinite(t.cellVectors[i]) || !std::isfinite(t.blockVectors[i])) {
         diagnostic = "GPU adaptive topology vectors must be finite";
         return false;
      }
   }
   if(r.selectorCriteria > 0 && !r.selectorScores) {
      diagnostic = "GPU adaptive selector score count is nonzero but storage was not staged";
      return false;
   }
   if(t.blockAtomOffset[0] != 0 ||
      t.blockAtomOffset[t.blocks] != static_cast<int>(t.atoms)) {
      diagnostic = "GPU adaptive block CSR offsets must span exactly all atoms";
      return false;
   }

   std::vector<unsigned char> seen(t.atoms, 0);
   std::vector<std::size_t> observedBlock(t.blocks, 0);
   std::vector<std::size_t> observedBasis(t.basis * t.blocks, 0);
   std::vector<std::size_t> observedFft(t.fftChannelsPerBlock * t.blocks, 0);
   std::vector<std::size_t> observedDynamic(t.dynamicChannels * t.blocks, 0);
   for(std::size_t block = 0; block < t.blocks; ++block) {
      if(t.blockAtomOffset[block] > t.blockAtomOffset[block + 1] ||
         t.blockAtomCount[block] != t.blockAtomOffset[block + 1] - t.blockAtomOffset[block]) {
         diagnostic = "GPU adaptive block CSR offsets/counts are inconsistent";
         return false;
      }
      if(r.blockState[block] < coarseState || r.blockState[block] > fineState) {
         diagnostic = "GPU adaptive block state must be coarse(0), buffer(1), or fine(2)";
         return false;
      }
      if(!std::isfinite(t.blockVolume[block]) || t.blockVolume[block] <= 0.0) {
         diagnostic = "GPU adaptive block volumes must be finite and positive";
         return false;
      }
      for(int d = 0; d < 3; ++d) {
         if(!std::isfinite(t.blockCenter[d + 3 * block]) ||
            t.blockGridCoordinate[d + 3 * block] < 0) {
            diagnostic = "GPU adaptive block centers must be finite";
            return false;
         }
      }
   }
   for(std::size_t basis = 0; basis < t.basis; ++basis) {
      const int dynamic = t.basisToDynamicChannel[basis];
      const int fft = t.basisToFftChannel[basis];
      if((dynamic != -1 && (dynamic < 1 || dynamic > static_cast<int>(t.dynamicChannels))) ||
         fft < 1 || fft > static_cast<int>(t.fftChannelsPerBlock)) {
         diagnostic = "GPU adaptive basis topology map contains an out-of-range id";
         return false;
      }
   }
   for(std::size_t atom = 0; atom < t.atoms; ++atom) {
      const int block = t.atomToBlock[atom] - 1;
      const int basis = t.atomToBasis[atom] - 1;
      const int fft = t.atomToFftChannel[atom] - 1;
      const int rawDynamic = t.atomToDynamicChannel[atom];
      const int dynamic = rawDynamic == -1 ? -1 : rawDynamic - 1;
      if(block < 0 || block >= static_cast<int>(t.blocks) ||
         basis < 0 || basis >= static_cast<int>(t.basis) ||
         fft < 0 || fft >= static_cast<int>(t.fftChannelsPerBlock) ||
         dynamic < -1 || dynamic >= static_cast<int>(t.dynamicChannels) ||
         t.atomToFftGridIndex[atom] < 1) {
         diagnostic = "GPU adaptive atom topology map contains an out-of-range id";
         return false;
      }
      if(t.basisToDynamicChannel[basis] != rawDynamic ||
         t.basisToFftChannel[basis] != t.atomToFftChannel[atom]) {
         diagnostic = "GPU adaptive atom and basis topology maps disagree";
         return false;
      }
      ++observedBlock[block];
      ++observedBasis[basis + t.basis * block];
      ++observedFft[fft + t.fftChannelsPerBlock * block];
      if(dynamic >= 0) ++observedDynamic[dynamic + t.dynamicChannels * block];
   }
   for(std::size_t position = 0; position < t.atoms; ++position) {
      const int atom = t.blockAtoms[position] - 1;
      if(atom < 0 || atom >= static_cast<int>(t.atoms) || seen[atom]) {
         diagnostic = "GPU adaptive block atom CSR is not a permutation";
         return false;
      }
      seen[atom] = 1;
      const auto upper = std::upper_bound(t.blockAtomOffset, t.blockAtomOffset + t.blocks + 1,
                                          static_cast<int>(position));
      const std::size_t block = static_cast<std::size_t>(upper - t.blockAtomOffset - 1);
      if(t.atomToBlock[atom] != static_cast<int>(block + 1)) {
         diagnostic = "GPU adaptive block atom CSR disagrees with atom_to_block";
         return false;
      }
   }
   for(std::size_t block = 0; block < t.blocks; ++block) {
      if(observedBlock[block] != static_cast<std::size_t>(t.blockAtomCount[block])) {
         diagnostic = "GPU adaptive block populations disagree with atom maps";
         return false;
      }
      for(std::size_t basis = 0; basis < t.basis; ++basis) {
         const std::size_t i = basis + t.basis * block;
         if(observedBasis[i] != static_cast<std::size_t>(t.blockBasisPopulation[i])) {
            diagnostic = "GPU adaptive basis populations disagree with atom maps";
            return false;
         }
      }
      for(std::size_t fft = 0; fft < t.fftChannelsPerBlock; ++fft) {
         const std::size_t i = fft + t.fftChannelsPerBlock * block;
         if(observedFft[i] != static_cast<std::size_t>(t.blockFftChannelPopulation[i])) {
            diagnostic = "GPU adaptive FFT populations disagree with atom maps";
            return false;
         }
      }
      for(std::size_t channel = 0; channel < t.dynamicChannels; ++channel) {
         const std::size_t i = channel + t.dynamicChannels * block;
         if(observedDynamic[i] != static_cast<std::size_t>(t.blockDynamicChannelPopulation[i])) {
            diagnostic = "GPU adaptive dynamic populations disagree with atom maps";
            return false;
         }
      }
   }
   const auto& k = r.kernels;
   if(k.atomMoment) {
      if(t.dynamicChannels != 1 || r.selectorCriteria == 0 ||
         !required(k.atomAnisotropyAxisCount) || !required(k.atomAnisotropyAxis) ||
         !required(k.atomAnisotropyK1) || !required(k.atomAnisotropyK2) ||
         !required(k.projectionBlock) || !required(k.projectionWeight) ||
         !required(k.inverseBlockTranspose) || !required(k.exchangeStiffness) ||
         !required(k.spiralization) || !required(k.anisotropyAxisCount) ||
         !required(k.anisotropyAxis) || !required(k.anisotropyK1) ||
         !required(k.anisotropyK2) ||
         (k.bonds > 0 && (!required(k.bondAtom) || !required(k.bondMatrix))) ||
         (k.selectorEdges > 0 && !required(k.selectorEdge))) {
         diagnostic = "GPU adaptive kernels require the complete accepted single-FM operator inventory";
         return false;
      }
      if(!std::isfinite(k.normalizationFloor) || k.normalizationFloor <= 0.0 ||
         !std::isfinite(k.magneticMomentSi) || k.magneticMomentSi <= 0.0 ||
         !std::isfinite(k.gammaPerTs) || k.gammaPerTs <= 0.0 ||
         !std::isfinite(k.damping) || k.damping < 0.0) {
         char values[256];
         std::snprintf(values, sizeof(values),
            "GPU adaptive kernel scalars are invalid (normalization=%g, moment_SI=%g, gamma=%g, damping=%g)",
            k.normalizationFloor, k.magneticMomentSi, k.gammaPerTs, k.damping);
         diagnostic = values;
         return false;
      }
      for(std::size_t atom = 0; atom < t.atoms; ++atom) {
         if(!std::isfinite(k.atomMoment[atom]) ||
            (t.atomToDynamicChannel[atom] > 0 && k.atomMoment[atom] <= 0.0)) {
            diagnostic = "GPU adaptive kernel atom moments must be finite and positive for magnetic atoms";
            return false;
         }
         if(k.atomAnisotropyAxisCount[atom] < 0 ||
            k.atomAnisotropyAxisCount[atom] > 2) {
            diagnostic = "GPU adaptive atomistic anisotropy supports zero, one, or two axes";
            return false;
         }
         for(int axis = 0; axis < 2; ++axis) {
            double norm2 = 0.0;
            for(int xyz = 0; xyz < 3; ++xyz) {
               const double value =
                  k.atomAnisotropyAxis[xyz + 3 * (axis + 2 * atom)];
               if(!std::isfinite(value)) {
                  diagnostic = "GPU adaptive atomistic anisotropy axes must be finite";
                  return false;
               }
               norm2 += value * value;
            }
            if(axis < k.atomAnisotropyAxisCount[atom] &&
               std::abs(norm2 - 1.0) > 1.0e-10) {
               diagnostic = "GPU adaptive active atomistic anisotropy axes must be normalized";
               return false;
            }
            if(!std::isfinite(k.atomAnisotropyK1[axis + 2 * atom]) ||
               !std::isfinite(k.atomAnisotropyK2[axis + 2 * atom])) {
               diagnostic = "GPU adaptive atomistic anisotropy coefficients must be finite";
               return false;
            }
         }
         double weightSum = 0.0;
         for(int corner = 0; corner < 8; ++corner) {
            const std::size_t index = static_cast<std::size_t>(corner) + 8 * atom;
            if(k.projectionBlock[index] < 1 ||
               k.projectionBlock[index] > static_cast<int>(t.blocks) ||
               !std::isfinite(k.projectionWeight[index])) {
               diagnostic = "GPU adaptive projection stencil contains an invalid block or weight";
               return false;
            }
            weightSum += k.projectionWeight[index];
         }
         if(std::abs(weightSum - 1.0) > 1.0e-12) {
            diagnostic = "GPU adaptive projection weights do not form a partition of unity";
            return false;
         }
      }
      for(std::size_t bond = 0; bond < k.bonds; ++bond) {
         const int atomI = k.bondAtom[2 * bond];
         const int atomJ = k.bondAtom[2 * bond + 1];
         if(atomI < 1 || atomI > static_cast<int>(t.atoms) ||
            atomJ < 1 || atomJ > static_cast<int>(t.atoms) || atomI == atomJ) {
            diagnostic = "GPU adaptive bond list contains an invalid unique pair";
            return false;
         }
         for(int component = 0; component < 9; ++component) {
            if(!std::isfinite(k.bondMatrix[component + 9 * bond])) {
               diagnostic = "GPU adaptive bond tensor must be finite";
               return false;
            }
         }
      }
      for(std::size_t edge = 0; edge < k.selectorEdges; ++edge) {
         const int atomI = k.selectorEdge[2 * edge];
         const int atomJ = k.selectorEdge[2 * edge + 1];
         if(atomI < 1 || atomI > static_cast<int>(t.atoms) ||
            atomJ < 1 || atomJ > static_cast<int>(t.atoms) || atomI == atomJ) {
            diagnostic = "GPU adaptive selector edge list contains an invalid pair";
            return false;
         }
      }
      for(int component = 0; component < 9; ++component) {
         if(!std::isfinite(k.inverseBlockTranspose[component]) ||
            !std::isfinite(k.exchangeStiffness[component]) ||
            !std::isfinite(k.spiralization[component])) {
            diagnostic = "GPU adaptive coarse tensor coefficients must be finite";
            return false;
         }
      }
      double exchangeScale = std::numeric_limits<double>::min();
      for(int component = 0; component < 9; ++component)
         exchangeScale = std::max(exchangeScale,
                                  std::abs(k.exchangeStiffness[component]));
      for(int row = 0; row < 3; ++row)
         for(int column = 0; column < 3; ++column) {
            if(std::abs(k.exchangeStiffness[row + 3 * column] -
                        k.exchangeStiffness[column + 3 * row]) >
               1.0e-12 * exchangeScale) {
               diagnostic = "GPU adaptive exchange stiffness must be Cartesian-symmetric";
               return false;
            }
         }
      for(std::size_t block = 1; block < t.blocks; ++block) {
         const double scale = std::max(
            std::abs(t.blockVolume[0]), std::abs(t.blockVolume[block]));
         if(std::abs(t.blockVolume[block] - t.blockVolume[0]) >
            1.0e-12 * scale) {
            diagnostic = "GPU adaptive tensor kernels require a uniform block volume";
            return false;
         }
      }
      for(std::size_t block = 0; block < t.blocks; ++block) {
         if(k.anisotropyAxisCount[block] < 0 || k.anisotropyAxisCount[block] > 2) {
            diagnostic = "GPU adaptive anisotropy supports zero, one, or two axes";
            return false;
         }
         for(int axis = 0; axis < 2; ++axis) {
            real norm2 = real(0);
            for(int xyz = 0; xyz < 3; ++xyz) {
               const double value = k.anisotropyAxis[
                  xyz + 3 * (axis + 2 * block)];
               if(!std::isfinite(value)) {
                  diagnostic = "GPU adaptive anisotropy axes must be finite";
                  return false;
               }
               norm2 += static_cast<real>(value * value);
            }
            if(axis < k.anisotropyAxisCount[block] &&
               std::abs(static_cast<double>(norm2) - 1.0) > 1.0e-10) {
               diagnostic = "GPU adaptive active anisotropy axes must be normalized";
               return false;
            }
            if(!std::isfinite(k.anisotropyK1[axis + 2 * block]) ||
               !std::isfinite(k.anisotropyK2[axis + 2 * block])) {
               diagnostic = "GPU adaptive anisotropy coefficients must be finite";
               return false;
            }
         }
      }
   }
   return true;
}

std::size_t GpuAdaptiveRuntime::estimateBytes(const GpuAdaptiveTopologyInput& t,
                                              const GpuAdaptiveRuntimeInput& r) {
   if(t.atoms == 0 || t.blocks == 0) return 0;
   std::size_t total = 0;
   std::size_t vectorState = 0;
   std::size_t scalarState = 0;
   if(!checkedProduct({3, t.dynamicChannels, t.blocks, t.ensembles}, vectorState) ||
      !checkedProduct({t.dynamicChannels, t.blocks, t.ensembles}, scalarState)) {
      throw std::overflow_error("GPU adaptive runtime memory estimate overflow");
   }
   const bool ok =
      checkedAdd(total, 6 * t.atoms, sizeof(int)) &&
      checkedAdd(total, 2 * t.basis, sizeof(int)) &&
      checkedAdd(total, 2 * t.blocks + 1, sizeof(int)) &&
      checkedAdd(total, 3 * t.blocks, sizeof(int)) &&
      checkedAdd(total, t.basis * t.blocks, sizeof(int)) &&
      checkedAdd(total, t.fftChannelsPerBlock * t.blocks, sizeof(int)) &&
      checkedAdd(total, t.dynamicChannels * t.blocks, sizeof(int)) &&
      checkedAdd(total, 4 * t.blocks, sizeof(real)) &&
      checkedAdd(total, 2 * t.blocks + 2 * t.atoms, sizeof(unsigned char)) &&
      checkedAdd(total, 2 * t.blocks, sizeof(int)) &&
      checkedAdd(total, 2 * t.blocks, sizeof(unsigned int)) &&
      checkedAdd(total, r.selectorCriteria * t.blocks, sizeof(real)) &&
      checkedAdd(total, 3 * vectorState + scalarState, sizeof(real)) &&
      checkedAdd(total, 2 * t.atoms + t.blocks, sizeof(int)) &&
      // RCG-05E: dilatedState_, the dilateAdaptiveState double buffer.
      checkedAdd(total, t.blocks, sizeof(int)) &&
      checkedAdd(total, 3, sizeof(unsigned int)) &&
      checkedAdd(total, 6 * std::max(t.atoms, t.blocks),
                 sizeof(unsigned int));
   if(!ok) throw std::overflow_error("GPU adaptive runtime memory estimate overflow");
   // RCG-08 (F-12): tile-total levels for the hierarchical compaction scan.
   // Derived from the same helper allocate() uses, so preflight and
   // allocation cannot drift.  The series is geometric in adaptiveThreads, so
   // this adds well under 1% of the level-0 scan buffers.
   for(const std::size_t levelItems :
       compactionScanLevelItems(std::max(t.atoms, t.blocks))) {
      if(!checkedAdd(total, 3 * levelItems, sizeof(unsigned int)))
         throw std::overflow_error("GPU adaptive runtime memory estimate overflow");
   }
   if(r.kernels.atomMoment) {
      std::size_t atomEnsembles = 0;
      if(!checkedProduct({t.atoms, t.ensembles}, atomEnsembles))
         throw std::overflow_error("GPU adaptive kernel memory estimate overflow");
      const bool kernelOk =
         checkedAdd(total, 9 * t.atoms + 2 * r.kernels.bonds +
                           2 * r.kernels.selectorEdges + t.blocks, sizeof(int)) &&
         checkedAdd(total, 11 * t.atoms + 8 * t.atoms + 9 * r.kernels.bonds + 27 +
                           10 * t.blocks + 13 * atomEnsembles +
                           4 * vectorState + t.blocks, sizeof(real)) &&
         // RCG-06B (F-11): energyTerms_'s 8 slots are FP64 unconditionally,
         // independent of `real`'s build precision, so they are no longer
         // part of the sizeof(real) bucket above.
         checkedAdd(total, 8, sizeof(double)) &&
         checkedAdd(total, 2 * t.blocks, sizeof(unsigned char));
      if(!kernelOk)
         throw std::overflow_error("GPU adaptive kernel memory estimate overflow");
   }
   return total;
}

void GpuAdaptiveRuntime::initialize(const GpuAdaptiveTopologyInput& topology,
                                    const GpuAdaptiveRuntimeInput& runtime,
                                    std::size_t expectedAtoms,
                                    std::size_t expectedEnsembles) {
   if(ready_ || streamCreated_) throw std::logic_error("GPU adaptive runtime is already initialized");
   metrics_ = {};
   phaseMetrics_ = {};
   std::string diagnostic;
   if(!validate(topology, runtime, expectedAtoms, expectedEnsembles, diagnostic))
      throw std::invalid_argument(diagnostic);

   atoms_ = topology.atoms;
   blocks_ = topology.blocks;
   dynamicChannels_ = topology.dynamicChannels;
   ensembles_ = topology.ensembles;
   allocatedBytes_ = estimateBytes(topology, runtime);
   try {
      ASSERT_GPU(GPU_STREAM_CREATE(&stream_));
      streamCreated_ = true;
      ASSERT_GPU(GPU_EVENT_CREATE(&updateStart_));
      startEventCreated_ = true;
      ASSERT_GPU(GPU_EVENT_CREATE(&updateEnd_));
      endEventCreated_ = true;
      ASSERT_GPU(GPU_EVENT_CREATE(&phaseStart_));
      ASSERT_GPU(GPU_EVENT_CREATE(&phaseEnd_));
      phaseEventsCreated_ = true;
      allocate(topology, runtime);
      deviceTopology_.geometryMode = topology.geometryMode;
      deviceTopology_.fftGridChannels = topology.fftGridChannels;
      for(int d = 0; d < 3; ++d) {
         deviceTopology_.repetitionShape[d] = topology.repetitionShape[d];
         deviceTopology_.blockShape[d] = topology.blockShape[d];
         deviceTopology_.blockGrid[d] = topology.blockGrid[d];
      }
      for(int i = 0; i < 9; ++i) {
         deviceTopology_.cellVectors[i] = static_cast<real>(topology.cellVectors[i]);
         deviceTopology_.blockVectors[i] = static_cast<real>(topology.blockVectors[i]);
      }
      uploadTopology(topology);
      uploadRuntime(topology, runtime);
      launchCompaction();
      ASSERT_GPU(GPU_STREAM_SYNC(stream_));
      convertedStaging_.clear();
      kernelsReady_ = runtime.kernels.atomMoment != nullptr;
      ready_ = true;
   } catch(...) {
      release();
      throw;
   }
}

void GpuAdaptiveRuntime::allocate(const GpuAdaptiveTopologyInput& t,
                                  const GpuAdaptiveRuntimeInput& r) {
   const auto n = static_cast<index_t>(t.atoms);
   const auto b = static_cast<index_t>(t.blocks);
   const auto basis = static_cast<index_t>(t.basis);
   const auto fft = static_cast<index_t>(t.fftChannelsPerBlock);
   const auto channels = static_cast<index_t>(t.dynamicChannels);
   const auto ensembles = static_cast<index_t>(t.ensembles);
   const auto criteria = static_cast<index_t>(r.selectorCriteria);

   atomToBlock_.Allocate(n);
   atomToBasis_.Allocate(n);
   atomToDynamicChannel_.Allocate(n);
   atomToFftChannel_.Allocate(n);
   atomToFftGridIndex_.Allocate(n);
   basisToDynamicChannel_.Allocate(basis);
   basisToFftChannel_.Allocate(basis);
   blockAtomCount_.Allocate(b);
   blockAtomOffset_.Allocate(b + 1);
   blockAtoms_.Allocate(n);
   blockGridCoordinate_.Allocate(3 * b);
   blockBasisPopulation_.Allocate(basis * b);
   blockFftChannelPopulation_.Allocate(fft * b);
   blockDynamicChannelPopulation_.Allocate(channels * b);
   blockCenter_.Allocate(3 * b);
   blockVolume_.Allocate(b);

   blockState_.Allocate(b);
   pendingState_.Allocate(b);
   dilatedState_.Allocate(b);
   stateAge_.Allocate(b);
   transitionEpoch_.Allocate(b);
   atomisticBlockMask_.Allocate(b);
   coarseBlockMask_.Allocate(b);
   atomisticAtomMask_.Allocate(n);
   interfaceAtomMask_.Allocate(n);
   if(criteria > 0) selectorScores_.Allocate(criteria * b);
   coarseMoment_.Allocate(3 * channels * b * ensembles);
   coarseDirection_.Allocate(3 * channels * b * ensembles);
   coarseField_.Allocate(3 * channels * b * ensembles);
   channelMomentSum_.Allocate(channels * b * ensembles);
   activeAtomList_.Allocate(n);
   activeBlockList_.Allocate(b);
   interfaceAtomList_.Allocate(n);
   workCounts_.Allocate(3);
   const auto scanItems = std::max(n, b);
   compactionScanA_.Allocate(3 * scanItems);
   compactionScanB_.Allocate(3 * scanItems);
   scanLevelItems_ =
      compactionScanLevelItems(static_cast<std::size_t>(scanItems));
   scanLevelOffset_.assign(scanLevelItems_.size(), 0);
   std::size_t scanLevelTotal = 0;
   for(std::size_t level = 0; level < scanLevelItems_.size(); ++level) {
      scanLevelOffset_[level] = scanLevelTotal;
      scanLevelTotal += 3 * scanLevelItems_[level];
   }
   if(scanLevelTotal > 0)
      compactionScanLevels_.Allocate(static_cast<index_t>(scanLevelTotal));
   if(r.kernels.atomMoment) {
      const auto atomEnsembles = n * ensembles;
      const auto vectorState = 3 * channels * b * ensembles;
      bonds_ = r.kernels.bonds;
      selectorEdges_ = r.kernels.selectorEdges;
      atomMoment_.Allocate(n);
      atomAnisotropyAxisCount_.Allocate(n);
      atomAnisotropyAxis_.Allocate(6 * n);
      atomAnisotropyK1_.Allocate(2 * n);
      atomAnisotropyK2_.Allocate(2 * n);
      projectionBlock_.Allocate(8 * n);
      projectionWeight_.Allocate(8 * n);
      if(bonds_ > 0) {
         bondAtom_.Allocate(static_cast<index_t>(2 * bonds_));
         bondMatrix_.Allocate(static_cast<index_t>(9 * bonds_));
      }
      if(selectorEdges_ > 0)
         selectorEdge_.Allocate(static_cast<index_t>(2 * selectorEdges_));
      inverseBlockTranspose_.Allocate(9);
      exchangeStiffness_.Allocate(9);
      spiralization_.Allocate(9);
      anisotropyAxisCount_.Allocate(b);
      anisotropyAxis_.Allocate(6 * b);
      anisotropyK1_.Allocate(2 * b);
      anisotropyK2_.Allocate(2 * b);
      ghostDirection_.Allocate(3 * atomEnsembles);
      projectionNorm_.Allocate(atomEnsembles);
      atomFieldScratch_.Allocate(3 * atomEnsembles);
      coarseFieldScratch_.Allocate(vectorState);
      predictorAtom_.Allocate(3 * atomEnsembles);
      predictorCoarse_.Allocate(vectorState);
      initialAtomField_.Allocate(3 * atomEnsembles);
      initialCoarseField_.Allocate(vectorState);
      predictorCoarseField_.Allocate(vectorState);
      energyTerms_.Allocate(8);
      acceptedBlockMask_.Allocate(b);
      polarizationUnsafeBlockMask_.Allocate(b);
      polarizationRatioBlock_.Allocate(b);
      // Make zero-step diagnostics deterministic before the first field
      // evaluation; normal kernels overwrite these buffers thereafter.
      atomFieldScratch_.zeros_async(stream_);
      predictorCoarseField_.zeros_async(stream_);
      energyTerms_.zeros_async(stream_);
   }
   stagedBlockState_.AllocateHost(b);
   refreshDeviceDescriptors();
}

void GpuAdaptiveRuntime::uploadTopology(const GpuAdaptiveTopologyInput& t) {
   uploadNative(atomToBlock_, t.atomToBlock, t.atoms, stream_);
   uploadNative(atomToBasis_, t.atomToBasis, t.atoms, stream_);
   uploadNative(atomToDynamicChannel_, t.atomToDynamicChannel, t.atoms, stream_);
   uploadNative(atomToFftChannel_, t.atomToFftChannel, t.atoms, stream_);
   uploadNative(atomToFftGridIndex_, t.atomToFftGridIndex, t.atoms, stream_);
   uploadNative(basisToDynamicChannel_, t.basisToDynamicChannel, t.basis, stream_);
   uploadNative(basisToFftChannel_, t.basisToFftChannel, t.basis, stream_);
   uploadNative(blockAtomCount_, t.blockAtomCount, t.blocks, stream_);
   uploadNative(blockAtomOffset_, t.blockAtomOffset, t.blocks + 1, stream_);
   uploadNative(blockAtoms_, t.blockAtoms, t.atoms, stream_);
   uploadNative(blockGridCoordinate_, t.blockGridCoordinate, 3 * t.blocks, stream_);
   uploadNative(blockBasisPopulation_, t.blockBasisPopulation, t.basis * t.blocks, stream_);
   uploadNative(blockFftChannelPopulation_, t.blockFftChannelPopulation,
                t.fftChannelsPerBlock * t.blocks, stream_);
   uploadNative(blockDynamicChannelPopulation_, t.blockDynamicChannelPopulation,
                t.dynamicChannels * t.blocks, stream_);
   uploadReal(blockCenter_, t.blockCenter, 3 * t.blocks, stream_, convertedStaging_);
   uploadReal(blockVolume_, t.blockVolume, t.blocks, stream_, convertedStaging_);
}

void GpuAdaptiveRuntime::uploadRuntime(const GpuAdaptiveTopologyInput& t,
                                      const GpuAdaptiveRuntimeInput& r) {
   std::copy(r.blockState, r.blockState + t.blocks, stagedBlockState_.data());
   blockState_.copy_async(stagedBlockState_, stream_);
   if(r.pendingState) uploadNative(pendingState_, r.pendingState, t.blocks, stream_);
   else pendingState_.copy_async(blockState_, stream_);
   if(r.stateAge) uploadNative(stateAge_, r.stateAge, t.blocks, stream_);
   else stateAge_.zeros_async(stream_);
   if(r.transitionEpoch) uploadNative(transitionEpoch_, r.transitionEpoch, t.blocks, stream_);
   else transitionEpoch_.zeros_async(stream_);
   if(r.selectorCriteria > 0)
      uploadReal(selectorScores_, r.selectorScores, r.selectorCriteria * t.blocks, stream_,
                 convertedStaging_);

   const std::size_t vectorState = 3 * t.dynamicChannels * t.blocks * t.ensembles;
   const std::size_t scalarState = t.dynamicChannels * t.blocks * t.ensembles;
   uploadReal(coarseMoment_, r.coarseMoment, vectorState, stream_, convertedStaging_);
   uploadReal(coarseDirection_, r.coarseDirection, vectorState, stream_, convertedStaging_);
   uploadReal(coarseField_, r.coarseField, vectorState, stream_, convertedStaging_);
   uploadReal(channelMomentSum_, r.channelMomentSum, scalarState, stream_, convertedStaging_);
   if(r.kernels.atomMoment) {
      const auto& k = r.kernels;
      uploadReal(atomMoment_, k.atomMoment, t.atoms, stream_, convertedStaging_);
      uploadNative(atomAnisotropyAxisCount_, k.atomAnisotropyAxisCount, t.atoms, stream_);
      uploadReal(atomAnisotropyAxis_, k.atomAnisotropyAxis, 6 * t.atoms, stream_,
                 convertedStaging_);
      uploadReal(atomAnisotropyK1_, k.atomAnisotropyK1, 2 * t.atoms, stream_,
                 convertedStaging_);
      uploadReal(atomAnisotropyK2_, k.atomAnisotropyK2, 2 * t.atoms, stream_,
                 convertedStaging_);
      uploadNative(projectionBlock_, k.projectionBlock, 8 * t.atoms, stream_);
      uploadReal(projectionWeight_, k.projectionWeight, 8 * t.atoms, stream_,
                 convertedStaging_);
      if(k.bonds > 0) {
         uploadNative(bondAtom_, k.bondAtom, 2 * k.bonds, stream_);
         uploadReal(bondMatrix_, k.bondMatrix, 9 * k.bonds, stream_,
                    convertedStaging_);
      }
      if(k.selectorEdges > 0)
         uploadNative(selectorEdge_, k.selectorEdge, 2 * k.selectorEdges, stream_);
      uploadReal(inverseBlockTranspose_, k.inverseBlockTranspose, 9, stream_,
                 convertedStaging_);
      uploadReal(exchangeStiffness_, k.exchangeStiffness, 9, stream_,
                 convertedStaging_);
      uploadReal(spiralization_, k.spiralization, 9, stream_, convertedStaging_);
      uploadNative(anisotropyAxisCount_, k.anisotropyAxisCount, t.blocks, stream_);
      uploadReal(anisotropyAxis_, k.anisotropyAxis, 6 * t.blocks, stream_,
                 convertedStaging_);
      uploadReal(anisotropyK1_, k.anisotropyK1, 2 * t.blocks, stream_,
                 convertedStaging_);
      uploadReal(anisotropyK2_, k.anisotropyK2, 2 * t.blocks, stream_,
                 convertedStaging_);
      normalizationFloor_ = static_cast<real>(k.normalizationFloor);
      magneticMomentSi_ = static_cast<real>(k.magneticMomentSi);
      gammaPerTs_ = static_cast<real>(k.gammaPerTs);
      damping_ = static_cast<real>(k.damping);
   }
}

void GpuAdaptiveRuntime::refreshDeviceDescriptors() {
   deviceTopology_.atoms = atoms_;
   deviceTopology_.blocks = blocks_;
   deviceTopology_.basis = static_cast<std::size_t>(basisToDynamicChannel_.size());
   deviceTopology_.fftChannelsPerBlock =
      blocks_ == 0 ? 0 :
      static_cast<std::size_t>(blockFftChannelPopulation_.size()) / blocks_;
   deviceTopology_.fftGridChannels =
      deviceTopology_.fftChannelsPerBlock * deviceTopology_.blocks;
   deviceTopology_.dynamicChannels = dynamicChannels_;
   deviceTopology_.ensembles = ensembles_;
   deviceTopology_.atomToBlock = atomToBlock_.data();
   deviceTopology_.atomToBasis = atomToBasis_.data();
   deviceTopology_.atomToDynamicChannel = atomToDynamicChannel_.data();
   deviceTopology_.atomToFftChannel = atomToFftChannel_.data();
   deviceTopology_.atomToFftGridIndex = atomToFftGridIndex_.data();
   deviceTopology_.basisToDynamicChannel = basisToDynamicChannel_.data();
   deviceTopology_.basisToFftChannel = basisToFftChannel_.data();
   deviceTopology_.blockAtomCount = blockAtomCount_.data();
   deviceTopology_.blockAtomOffset = blockAtomOffset_.data();
   deviceTopology_.blockAtoms = blockAtoms_.data();
   deviceTopology_.blockGridCoordinate = blockGridCoordinate_.data();
   deviceTopology_.blockBasisPopulation = blockBasisPopulation_.data();
   deviceTopology_.blockFftChannelPopulation = blockFftChannelPopulation_.data();
   deviceTopology_.blockDynamicChannelPopulation = blockDynamicChannelPopulation_.data();
   deviceTopology_.blockCenter = blockCenter_.data();
   deviceTopology_.blockVolume = blockVolume_.data();

   deviceRuntime_.selectorCriteria =
      selectorScores_.empty() ? 0 : static_cast<std::size_t>(selectorScores_.size() / blocks_);
   deviceRuntime_.blockState = blockState_.data();
   deviceRuntime_.pendingState = pendingState_.data();
   deviceRuntime_.stateAge = stateAge_.data();
   deviceRuntime_.transitionEpoch = transitionEpoch_.data();
   deviceRuntime_.atomisticBlockMask = atomisticBlockMask_.data();
   deviceRuntime_.coarseBlockMask = coarseBlockMask_.data();
   deviceRuntime_.polarizationUnsafeMask = polarizationUnsafeBlockMask_.data();
   deviceRuntime_.polarizationRatio = polarizationRatioBlock_.data();
   deviceRuntime_.atomisticAtomMask = atomisticAtomMask_.data();
   deviceRuntime_.interfaceAtomMask = interfaceAtomMask_.data();
   deviceRuntime_.activeAtomList = activeAtomList_.data();
   deviceRuntime_.activeBlockList = activeBlockList_.data();
   deviceRuntime_.interfaceAtomList = interfaceAtomList_.data();
   deviceRuntime_.workCounts = workCounts_.data();
   deviceRuntime_.selectorScores = selectorScores_.data();
   deviceRuntime_.coarseMoment = coarseMoment_.data();
   deviceRuntime_.coarseDirection = coarseDirection_.data();
   deviceRuntime_.coarseField = coarseField_.data();
   deviceRuntime_.channelMomentSum = channelMomentSum_.data();
}

// RCG-08 (F-12): linear-work hierarchical compaction.  See the comment above
// scanAdaptiveTiles for the algorithm and its work/launch accounting.  The
// level structure is fixed at allocate() time (scanLevelItems_/scanLevelOffset_)
// so this hot path performs no allocation and no host synchronization.
void GpuAdaptiveRuntime::launchCompaction() {
   const std::size_t n0 = std::max(atoms_, blocks_);
   unsigned int* flags = compactionScanA_.data();
   unsigned int* scanned = compactionScanB_.data();
   unsigned int* levels =
      scanLevelItems_.empty() ? nullptr : compactionScanLevels_.data();
   const std::size_t levelCount = scanLevelItems_.size();
   const auto tilesOf = [](std::size_t items) {
      return (items + adaptiveThreads - 1) / adaptiveThreads;
   };

   ADAPTIVE_LAUNCH(initializeAdaptiveWorkScan, adaptiveGrid(n0), stream_,
                   deviceTopology_, deviceRuntime_, flags, n0);
   ++phaseMetrics_.compactionLaunches;

   // Level 0: scan every tile of all three components, publishing tile totals
   // into level 1 when there is more than one tile.
   const std::size_t tiles0 = tilesOf(n0);
   ADAPTIVE_LAUNCH(scanAdaptiveTiles, 3 * tiles0, stream_,
                   flags, scanned,
                   levelCount ? levels + scanLevelOffset_[0] : nullptr,
                   n0, tiles0);
   ++phaseMetrics_.compactionLaunches;

   // Levels 1..L: scan the tile totals in place, each publishing the level
   // above it.  The last level fits in one tile and publishes nothing.
   for(std::size_t level = 0; level < levelCount; ++level) {
      const std::size_t items = scanLevelItems_[level];
      unsigned int* values = levels + scanLevelOffset_[level];
      unsigned int* totals = level + 1 < levelCount ?
         levels + scanLevelOffset_[level + 1] : nullptr;
      ADAPTIVE_LAUNCH(scanAdaptiveTiles, 3 * tilesOf(items), stream_,
                      values, values, totals, items, tilesOf(items));
      ++phaseMetrics_.compactionLaunches;
   }

   // Propagate each level's scanned tile totals back down, finishing with
   // level 0's own output.
   for(std::size_t level = levelCount; level-- > 1;) {
      const std::size_t items = scanLevelItems_[level - 1];
      ADAPTIVE_LAUNCH(addAdaptiveTileOffsets, adaptiveGrid(3 * items), stream_,
                      levels + scanLevelOffset_[level - 1],
                      levels + scanLevelOffset_[level],
                      items, tilesOf(items));
      ++phaseMetrics_.compactionLaunches;
   }
   if(levelCount) {
      ADAPTIVE_LAUNCH(addAdaptiveTileOffsets, adaptiveGrid(3 * n0), stream_,
                      scanned, levels + scanLevelOffset_[0], n0, tiles0);
      ++phaseMetrics_.compactionLaunches;
   }

   ADAPTIVE_LAUNCH(scatterAdaptiveWork, adaptiveGrid(n0), stream_,
                   deviceTopology_, deviceRuntime_, scanned, n0);
   ++phaseMetrics_.compactionLaunches;
   ASSERT_GPU(GPU_GET_LAST_ERROR());
}

void GpuAdaptiveRuntime::updateBlockState(const int* blockState, std::size_t count) {
   if(!ready_) throw std::logic_error("GPU adaptive runtime is not initialized");
   if(!blockState || count != blocks_)
      throw std::invalid_argument("GPU adaptive block-state update must contain exactly nblocks entries");
   for(std::size_t block = 0; block < blocks_; ++block) {
      if(blockState[block] < coarseState || blockState[block] > fineState)
         throw std::invalid_argument("GPU adaptive block state must be coarse(0), buffer(1), or fine(2)");
   }
   std::copy(blockState, blockState + blocks_, stagedBlockState_.data());
   const auto wallStart = std::chrono::steady_clock::now();
   ASSERT_GPU(GPU_EVENT_RECORD(updateStart_, stream_));
   blockState_.copy_async(stagedBlockState_, stream_);
   pendingState_.copy_async(blockState_, stream_);
   launchCompaction();
   ASSERT_GPU(GPU_EVENT_RECORD(updateEnd_, stream_));
   const auto waitStart = std::chrono::steady_clock::now();
   ASSERT_GPU(GPU_EVENT_SYNCHRONIZE(updateEnd_));
   const auto wallEnd = std::chrono::steady_clock::now();
   float elapsed = 0.0f;
   ASSERT_GPU(GPU_EVENT_ELAPSED_TIME(&elapsed, updateStart_, updateEnd_));
   ++metrics_.hostSynchronizations;
   metrics_.blockBytesUploaded += blocks_ * sizeof(int);
   metrics_.elapsedMilliseconds +=
      std::chrono::duration<double, std::milli>(wallEnd - wallStart).count();
   metrics_.hostWaitMilliseconds +=
      std::chrono::duration<double, std::milli>(wallEnd - waitStart).count();
   metrics_.deviceMilliseconds += static_cast<double>(elapsed);
   phaseMetrics_.compactionMilliseconds += static_cast<double>(elapsed);
}

double GpuAdaptiveRuntime::finishPhase(double& accumulator) {
   ASSERT_GPU(GPU_EVENT_RECORD(phaseEnd_, stream_));
   // RCG-08: this event synchronization blocks the host until the phase's
   // kernels retire.  It is how RCG-06C obtains per-phase device times, but it
   // also means each phase boundary is a real host wait that prevents adjacent
   // phases from overlapping -- the largest remaining serial section in the
   // adaptive step.  Counted here so the benchmark can report it rather than
   // leaving it to be inferred from unaccounted time.
   ++phaseMetrics_.phaseSynchronizations;
   ASSERT_GPU(GPU_EVENT_SYNCHRONIZE(phaseEnd_));
   float elapsed = 0.0f;
   ASSERT_GPU(GPU_EVENT_ELAPSED_TIME(&elapsed, phaseStart_, phaseEnd_));
   accumulator += static_cast<double>(elapsed);
   return static_cast<double>(elapsed);
}

void GpuAdaptiveRuntime::restrictMoments(const real* atomDirection) {
   if(!ready_ || !kernelsReady_)
      throw std::logic_error("GPU adaptive restriction requires initialized CG-10 kernels");
   if(!atomDirection) throw std::invalid_argument("GPU adaptive restriction requires device directions");
   AdaptiveKernelDevice kernels{
      bonds_, selectorEdges_, normalizationFloor_, magneticMomentSi_,
      gammaPerTs_, damping_, atomMoment_.data(), atomAnisotropyAxisCount_.data(),
      atomAnisotropyAxis_.data(), atomAnisotropyK1_.data(), atomAnisotropyK2_.data(),
      projectionBlock_.data(),
      projectionWeight_.data(), bondAtom_.data(), bondMatrix_.data(),
      selectorEdge_.data(), inverseBlockTranspose_.data(),
      exchangeStiffness_.data(), spiralization_.data(),
      anisotropyAxisCount_.data(), anisotropyAxis_.data(),
      anisotropyK1_.data(), anisotropyK2_.data(), ghostDirection_.data(),
      projectionNorm_.data(), atomFieldScratch_.data(),
      coarseFieldScratch_.data(), energyTerms_.data(), predictorAtom_.data()
   };
   ASSERT_GPU(GPU_EVENT_RECORD(phaseStart_, stream_));
   ADAPTIVE_LAUNCH(restrictAdaptiveMoments,
                   adaptiveGrid(dynamicChannels_ * blocks_ * ensembles_),
                   stream_, deviceTopology_, deviceRuntime_, kernels,
                   atomDirection);
   ++phaseMetrics_.integrationLaunches;
   ASSERT_GPU(GPU_GET_LAST_ERROR());
   finishPhase(phaseMetrics_.integrationMilliseconds);
}

void GpuAdaptiveRuntime::evaluateSelectorScores(const real* atomDirection) {
   if(!ready_ || !kernelsReady_)
      throw std::logic_error("GPU adaptive selector requires initialized CG-10 kernels");
   if(!atomDirection) throw std::invalid_argument("GPU adaptive selector requires device directions");
   AdaptiveKernelDevice kernels{
      bonds_, selectorEdges_, normalizationFloor_, magneticMomentSi_,
      gammaPerTs_, damping_, atomMoment_.data(), atomAnisotropyAxisCount_.data(),
      atomAnisotropyAxis_.data(), atomAnisotropyK1_.data(), atomAnisotropyK2_.data(),
      projectionBlock_.data(),
      projectionWeight_.data(), bondAtom_.data(), bondMatrix_.data(),
      selectorEdge_.data(), inverseBlockTranspose_.data(),
      exchangeStiffness_.data(), spiralization_.data(),
      anisotropyAxisCount_.data(), anisotropyAxis_.data(),
      anisotropyK1_.data(), anisotropyK2_.data(), ghostDirection_.data(),
      projectionNorm_.data(), atomFieldScratch_.data(),
      coarseFieldScratch_.data(), energyTerms_.data(), predictorAtom_.data()
   };
   ASSERT_GPU(GPU_EVENT_RECORD(phaseStart_, stream_));
#if defined(CUDA_V)
   clearSelectorAdaptiveScores<<<
      adaptiveGrid(blocks_), adaptiveThreads, 0, stream_>>>(
      deviceTopology_, deviceRuntime_);
   selectorAdaptiveScores<<<
      adaptiveGrid(selectorEdges_ * ensembles_), adaptiveThreads, 0, stream_>>>(
      deviceTopology_, deviceRuntime_, kernels, atomDirection);
#else
   hipLaunchKernelGGL(
      clearSelectorAdaptiveScores, dim3(adaptiveGrid(blocks_)),
      dim3(adaptiveThreads), 0, stream_, deviceTopology_, deviceRuntime_);
   hipLaunchKernelGGL(
      selectorAdaptiveScores,
      dim3(adaptiveGrid(selectorEdges_ * ensembles_)),
      dim3(adaptiveThreads), 0, stream_, deviceTopology_, deviceRuntime_,
      kernels, atomDirection);
#endif
   phaseMetrics_.selectorLaunches += 2;
   ASSERT_GPU(GPU_GET_LAST_ERROR());
   finishPhase(phaseMetrics_.selectorMilliseconds);
}

void GpuAdaptiveRuntime::evaluatePolarizationGate(real polarizationThreshold) {
   if(!ready_ || !kernelsReady_)
      throw std::logic_error("GPU adaptive polarization gate requires initialized CG-10 kernels");
   ASSERT_GPU(GPU_EVENT_RECORD(phaseStart_, stream_));
#if defined(CUDA_V)
   evaluateAdaptivePolarizationGate<<<
      adaptiveGrid(blocks_), adaptiveThreads, 0, stream_>>>(
      deviceTopology_, deviceRuntime_, polarizationThreshold);
#else
   hipLaunchKernelGGL(
      evaluateAdaptivePolarizationGate, dim3(adaptiveGrid(blocks_)),
      dim3(adaptiveThreads), 0, stream_, deviceTopology_, deviceRuntime_,
      polarizationThreshold);
#endif
   ++phaseMetrics_.polarizationLaunches;
   ASSERT_GPU(GPU_GET_LAST_ERROR());
   finishPhase(phaseMetrics_.polarizationMilliseconds);
}

void GpuAdaptiveRuntime::proposeSelectorState(
   const GpuAdaptiveSelectorPolicy& policy,
   const unsigned char* hardAtomisticBlockMask) {
   if(!ready_ || !kernelsReady_)
      throw std::logic_error("GPU adaptive state proposal requires initialized CG-10 kernels");
   if(policy.coarsenThreshold > policy.refineThreshold)
      throw std::invalid_argument("GPU adaptive selector coarsen threshold exceeds refine threshold");
   ASSERT_GPU(GPU_EVENT_RECORD(phaseStart_, stream_));
   const bool dilate = anyBufferDilation(policy.bufferDilationBlocks);
#if defined(CUDA_V)
   proposeAdaptiveState<<<
      adaptiveGrid(blocks_), adaptiveThreads, 0, stream_>>>(
      deviceTopology_, deviceRuntime_, policy, hardAtomisticBlockMask);
   if(dilate) {
      // RCG-05E: dilatedState_ starts as a copy of pendingState_ so blocks
      // the kernel never visits (already fine, or coarse with no fine
      // neighbour) keep their proposeAdaptiveState value; see the invariant
      // comment above dilateAdaptiveState for why this makes the read/write
      // sets of that kernel disjoint.
      dilatedState_.copy_async(pendingState_, stream_);
      dilateAdaptiveState<<<
         adaptiveGrid(blocks_), adaptiveThreads, 0, stream_>>>(
         deviceTopology_, deviceRuntime_, policy, dilatedState_.data());
      pendingState_.copy_async(dilatedState_, stream_);
   }
#else
   hipLaunchKernelGGL(
      proposeAdaptiveState, dim3(adaptiveGrid(blocks_)),
      dim3(adaptiveThreads), 0, stream_, deviceTopology_, deviceRuntime_,
      policy, hardAtomisticBlockMask);
   if(dilate) {
      dilatedState_.copy_async(pendingState_, stream_);
      hipLaunchKernelGGL(
         dilateAdaptiveState, dim3(adaptiveGrid(blocks_)),
         dim3(adaptiveThreads), 0, stream_, deviceTopology_, deviceRuntime_,
         policy, dilatedState_.data());
      pendingState_.copy_async(dilatedState_, stream_);
   }
#endif
   phaseMetrics_.selectorLaunches += dilate ? 2 : 1;
   ASSERT_GPU(GPU_GET_LAST_ERROR());
   finishPhase(phaseMetrics_.selectorMilliseconds);
}

void GpuAdaptiveRuntime::publishProposedState(
   real* atomDirection, const GpuAdaptiveReconstructionPolicy& policy,
   bool completeIntegrationStep, const unsigned char* acceptedBlockMask) {
   if(!ready_ || !kernelsReady_)
      throw std::logic_error("GPU adaptive publication requires initialized CG-10 kernels");
   if(!completeIntegrationStep)
      throw std::invalid_argument("GPU adaptive state may be published only after a complete integration step");
   if(!atomDirection)
      throw std::invalid_argument("GPU adaptive reconstruction requires device directions");
   if(policy.scheme != GpuAdaptiveReconstruction::Aligned &&
      policy.scheme != GpuAdaptiveReconstruction::ConstrainedCone)
      throw std::invalid_argument("GPU adaptive reconstruction scheme is invalid");
   if(policy.coneAngleRadians < real(0) ||
      policy.coneAngleRadians > real(0.5) * piValue ||
      policy.resultantTolerance <= real(0))
      throw std::invalid_argument("GPU adaptive reconstruction tolerance or cone angle is invalid");
   AdaptiveKernelDevice kernels{
      bonds_, selectorEdges_, normalizationFloor_, magneticMomentSi_,
      gammaPerTs_, damping_, atomMoment_.data(), atomAnisotropyAxisCount_.data(),
      atomAnisotropyAxis_.data(), atomAnisotropyK1_.data(), atomAnisotropyK2_.data(),
      projectionBlock_.data(),
      projectionWeight_.data(), bondAtom_.data(), bondMatrix_.data(),
      selectorEdge_.data(), inverseBlockTranspose_.data(),
      exchangeStiffness_.data(), spiralization_.data(),
      anisotropyAxisCount_.data(), anisotropyAxis_.data(),
      anisotropyK1_.data(), anisotropyK2_.data(), ghostDirection_.data(),
      projectionNorm_.data(), atomFieldScratch_.data(),
      coarseFieldScratch_.data(), energyTerms_.data(), predictorAtom_.data()
   };
   ASSERT_GPU(GPU_EVENT_RECORD(phaseStart_, stream_));
   ADAPTIVE_LAUNCH(publishAdaptiveState, adaptiveGrid(blocks_), stream_,
                   deviceTopology_, deviceRuntime_, kernels, atomDirection,
                   policy, acceptedBlockMask);
   ++phaseMetrics_.integrationLaunches;
   ASSERT_GPU(GPU_GET_LAST_ERROR());
   finishPhase(phaseMetrics_.integrationMilliseconds);
   ASSERT_GPU(GPU_EVENT_RECORD(phaseStart_, stream_));
   launchCompaction();
   finishPhase(phaseMetrics_.compactionMilliseconds);
}

GpuAdaptiveEnergy GpuAdaptiveRuntime::evaluateHybrid(
   const real* atomDirection, const real* externalCoarseField,
   const real* uniformFftDipoleField, real* atomField, real* coarseField,
   const GpuAdaptiveUniformFftField* basisResolvedFftField) {
   if(!ready_ || !kernelsReady_)
      throw std::logic_error("GPU adaptive evaluation requires initialized CG-10 kernels");
   if(!atomDirection || !atomField || !coarseField)
      throw std::invalid_argument("GPU adaptive evaluation requires device direction/field arrays");
   if(basisResolvedFftField) {
      if(!basisResolvedFftField->valid())
         throw std::invalid_argument("GPU adaptive FFT field view is incomplete");
      if(basisResolvedFftField->activeN1 !=
            static_cast<std::size_t>(deviceTopology_.blockGrid[0]) ||
         basisResolvedFftField->activeN2 !=
            static_cast<std::size_t>(deviceTopology_.blockGrid[1]) ||
         basisResolvedFftField->activeN3 !=
            static_cast<std::size_t>(deviceTopology_.blockGrid[2]) ||
         basisResolvedFftField->basis != deviceTopology_.fftChannelsPerBlock)
         throw std::invalid_argument(
            "GPU adaptive FFT field grid/channel mapping does not match BlockTopology");
   }
   AdaptiveKernelDevice kernels{
      bonds_, selectorEdges_, normalizationFloor_, magneticMomentSi_,
      gammaPerTs_, damping_, atomMoment_.data(), atomAnisotropyAxisCount_.data(),
      atomAnisotropyAxis_.data(), atomAnisotropyK1_.data(), atomAnisotropyK2_.data(),
      projectionBlock_.data(),
      projectionWeight_.data(), bondAtom_.data(), bondMatrix_.data(),
      selectorEdge_.data(), inverseBlockTranspose_.data(),
      exchangeStiffness_.data(), spiralization_.data(),
      anisotropyAxisCount_.data(), anisotropyAxis_.data(),
      anisotropyK1_.data(), anisotropyK2_.data(), ghostDirection_.data(),
      projectionNorm_.data(), atomFieldScratch_.data(),
      coarseFieldScratch_.data(), energyTerms_.data(), predictorAtom_.data()
   };
   ASSERT_GPU(GPU_EVENT_RECORD(phaseStart_, stream_));
#if defined(CUDA_V)
   prolongateAdaptiveGhosts<<<
      adaptiveGrid(atoms_ * ensembles_), adaptiveThreads, 0, stream_>>>(
      deviceTopology_, deviceRuntime_, kernels);
#else
   hipLaunchKernelGGL(
      prolongateAdaptiveGhosts, dim3(adaptiveGrid(atoms_ * ensembles_)),
      dim3(adaptiveThreads), 0, stream_, deviceTopology_, deviceRuntime_,
      kernels);
#endif
   ++phaseMetrics_.interfaceLaunches;
   phaseMetrics_.interfaceLaunches += 2;
   ASSERT_GPU(GPU_GET_LAST_ERROR());
   finishPhase(phaseMetrics_.interfaceMilliseconds);

   ASSERT_GPU(GPU_EVENT_RECORD(phaseStart_, stream_));
   // RCG-08 (F-09): four ordered launches replacing one serial kernel.  The
   // clear must precede the bond scatter (which accumulates atomically into
   // the scratch it zeroes), the bond scatter must precede the on-site pass
   // (which subtracts into the same scratch non-atomically at atoms it owns
   // exclusively), and the writeback must be last.  Stream order supplies
   // exactly that: consecutive launches on one stream do not overlap.
   ADAPTIVE_LAUNCH(clearAdaptiveAtomistic,
                   adaptiveGrid(3 * atoms_ * ensembles_), stream_,
                   deviceTopology_, kernels, atomField);
   ADAPTIVE_LAUNCH(evaluateAdaptiveAtomisticBonds,
                   adaptiveGrid(bonds_ * ensembles_), stream_,
                   deviceTopology_, deviceRuntime_, kernels, atomDirection);
   ADAPTIVE_LAUNCH(evaluateAdaptiveAtomisticOnsite,
                   adaptiveGrid(atoms_ * ensembles_), stream_,
                   deviceTopology_, deviceRuntime_, kernels, atomDirection);
   ADAPTIVE_LAUNCH(writebackAdaptiveAtomistic,
                   adaptiveGrid(atoms_ * ensembles_), stream_,
                   deviceTopology_, deviceRuntime_, kernels, atomField);
   phaseMetrics_.atomisticLaunches += 4;
   ASSERT_GPU(GPU_GET_LAST_ERROR());
   finishPhase(phaseMetrics_.atomisticMilliseconds);

   ASSERT_GPU(GPU_EVENT_RECORD(phaseStart_, stream_));
#if defined(CUDA_V)
   clearAdaptiveInterface<<<
      adaptiveGrid(3 * dynamicChannels_ * blocks_ * ensembles_),
      adaptiveThreads, 0, stream_>>>(
      deviceTopology_, deviceRuntime_, kernels);
   restrictAdaptiveInterface<<<
      adaptiveGrid(atoms_ * ensembles_), adaptiveThreads, 0, stream_>>>(
      deviceTopology_, deviceRuntime_, kernels);
#else
   hipLaunchKernelGGL(
      clearAdaptiveInterface,
      dim3(adaptiveGrid(3 * dynamicChannels_ * blocks_ * ensembles_)),
      dim3(adaptiveThreads), 0, stream_, deviceTopology_, deviceRuntime_,
      kernels);
   hipLaunchKernelGGL(
      restrictAdaptiveInterface, dim3(adaptiveGrid(atoms_ * ensembles_)),
      dim3(adaptiveThreads), 0, stream_, deviceTopology_, deviceRuntime_,
      kernels);
#endif
   ASSERT_GPU(GPU_GET_LAST_ERROR());
   finishPhase(phaseMetrics_.interfaceMilliseconds);

   ASSERT_GPU(GPU_EVENT_RECORD(phaseStart_, stream_));
#if defined(CUDA_V)
   clearAdaptiveCoarse<<<
      adaptiveGrid(3 * dynamicChannels_ * blocks_ * ensembles_),
      adaptiveThreads, 0, stream_>>>(
      deviceTopology_, deviceRuntime_, kernels, coarseField);
   evaluateAdaptiveCoarseTensor<<<
      adaptiveGrid(blocks_ * ensembles_), adaptiveThreads, 0, stream_>>>(
      deviceTopology_, deviceRuntime_, kernels);
   finalizeAdaptiveCoarseLocal<<<
      adaptiveGrid(blocks_ * ensembles_), adaptiveThreads, 0, stream_>>>(
      deviceTopology_, deviceRuntime_, kernels, externalCoarseField);
   if(uniformFftDipoleField)
      addAdaptiveDipole<<<
         adaptiveGrid(blocks_ * ensembles_), adaptiveThreads, 0, stream_>>>(
         deviceTopology_, deviceRuntime_, kernels, atomDirection,
         uniformFftDipoleField, atomField);
   if(basisResolvedFftField)
      addAdaptiveBasisResolvedDipole<<<
         adaptiveGrid(blocks_ * ensembles_), adaptiveThreads, 0, stream_>>>(
         deviceTopology_, deviceRuntime_, kernels, atomDirection,
         *basisResolvedFftField, atomField);
   writeAdaptiveCoarse<<<
      adaptiveGrid(blocks_ * ensembles_), adaptiveThreads, 0, stream_>>>(
      deviceTopology_, deviceRuntime_, kernels, coarseField);
   finalizeAdaptiveEnergy<<<1, 1, 0, stream_>>>(kernels);
#else
   hipLaunchKernelGGL(
      clearAdaptiveCoarse,
      dim3(adaptiveGrid(3 * dynamicChannels_ * blocks_ * ensembles_)),
      dim3(adaptiveThreads), 0, stream_, deviceTopology_, deviceRuntime_,
      kernels, coarseField);
   hipLaunchKernelGGL(
      evaluateAdaptiveCoarseTensor,
      dim3(adaptiveGrid(blocks_ * ensembles_)), dim3(adaptiveThreads), 0,
      stream_, deviceTopology_, deviceRuntime_, kernels);
   hipLaunchKernelGGL(
      finalizeAdaptiveCoarseLocal,
      dim3(adaptiveGrid(blocks_ * ensembles_)), dim3(adaptiveThreads), 0,
      stream_, deviceTopology_, deviceRuntime_, kernels,
      externalCoarseField);
   if(uniformFftDipoleField)
      hipLaunchKernelGGL(
         addAdaptiveDipole, dim3(adaptiveGrid(blocks_ * ensembles_)),
         dim3(adaptiveThreads), 0, stream_, deviceTopology_, deviceRuntime_,
         kernels, atomDirection, uniformFftDipoleField, atomField);
   if(basisResolvedFftField)
      hipLaunchKernelGGL(
         addAdaptiveBasisResolvedDipole,
         dim3(adaptiveGrid(blocks_ * ensembles_)),
         dim3(adaptiveThreads), 0, stream_, deviceTopology_, deviceRuntime_,
         kernels, atomDirection, *basisResolvedFftField, atomField);
   hipLaunchKernelGGL(
      writeAdaptiveCoarse, dim3(adaptiveGrid(blocks_ * ensembles_)),
      dim3(adaptiveThreads), 0, stream_, deviceTopology_, deviceRuntime_,
      kernels, coarseField);
   hipLaunchKernelGGL(finalizeAdaptiveEnergy, dim3(1), dim3(1), 0, stream_,
                      kernels);
#endif
   phaseMetrics_.coarseLaunches += 5 +
      (uniformFftDipoleField ? 1 : 0) + (basisResolvedFftField ? 1 : 0);
   ASSERT_GPU(GPU_GET_LAST_ERROR());
   finishPhase(phaseMetrics_.coarseMilliseconds);

   // RCG-06B (F-11): energyTerms_ is FP64 device storage now, so this
   // readback is a plain double copy -- no per-term precision loss remains
   // between the device reduction and this struct.
   double terms[8] = {};
   TensorDataMovementTracker::add_d2h(sizeof(terms));
   ASSERT_GPU(GPU_MEMCPY(terms, energyTerms_.data(), sizeof(terms),
                         GPU_MEMCPY_DEVICE_TO_HOST));
   GpuAdaptiveEnergy result;
   result.atomisticBilinearJ = terms[0];
   result.atomisticOnsiteJ = terms[1];
   result.coarseExchangeJ = terms[2];
   result.coarseSpiralizationJ = terms[3];
   result.coarseAnisotropyJ = terms[4];
   result.coarseExternalJ = terms[5];
   result.dipoleJ = terms[6];
   result.totalJ = terms[7];
   lastEnergy_ = result;
   return result;
}

void GpuAdaptiveRuntime::integrateHeun(
   real timeStepSeconds, real* atomDirection, const real* externalCoarseField,
   const real* uniformFftDipoleField,
   const GpuAdaptiveFftEvaluator& basisResolvedFftEvaluator) {
   if(!ready_ || !kernelsReady_)
      throw std::logic_error("GPU adaptive integration requires initialized CG-10 kernels");
   if(!atomDirection || !std::isfinite(static_cast<double>(timeStepSeconds)) ||
      timeStepSeconds <= real(0))
      throw std::invalid_argument("GPU adaptive Heun step requires positive dt and device directions");
   AdaptiveKernelDevice kernels{
      bonds_, selectorEdges_, normalizationFloor_, magneticMomentSi_,
      gammaPerTs_, damping_, atomMoment_.data(), atomAnisotropyAxisCount_.data(),
      atomAnisotropyAxis_.data(), atomAnisotropyK1_.data(), atomAnisotropyK2_.data(),
      projectionBlock_.data(),
      projectionWeight_.data(), bondAtom_.data(), bondMatrix_.data(),
      selectorEdge_.data(), inverseBlockTranspose_.data(),
      exchangeStiffness_.data(), spiralization_.data(),
      anisotropyAxisCount_.data(), anisotropyAxis_.data(),
      anisotropyK1_.data(), anisotropyK2_.data(), ghostDirection_.data(),
      projectionNorm_.data(), atomFieldScratch_.data(),
      coarseFieldScratch_.data(), energyTerms_.data(), predictorAtom_.data()
   };
   GpuAdaptiveUniformFftField initialFftField{};
   const GpuAdaptiveUniformFftField* initialFftFieldPtr = nullptr;
   if(basisResolvedFftEvaluator) {
      initialFftField = basisResolvedFftEvaluator(atomDirection);
      initialFftFieldPtr = &initialFftField;
   }
   (void)evaluateHybrid(atomDirection, externalCoarseField, uniformFftDipoleField,
                        initialAtomField_.data(), initialCoarseField_.data(),
                        initialFftFieldPtr);
   ASSERT_GPU(GPU_EVENT_RECORD(phaseStart_, stream_));
   // RCG-08 (F-09): the two full-state saves were the serial kernel's first
   // two loops.  They are pure copies, so they are issued as stream-ordered
   // device-to-device transfers rather than as kernel work -- the copy engine
   // does them at memory bandwidth and the following kernels still see them
   // complete, because everything here is ordered on stream_.  They remain an
   // O(atoms + blocks) pass per step regardless of the active fraction; that
   // is recorded as a known remaining full-system cost rather than hidden.
   ASSERT_GPU(GPU_MEMCPY_ASYNC(
      predictorAtom_.data(), atomDirection,
      3 * atoms_ * ensembles_ * sizeof(real),
      GPU_MEMCPY_DEVICE_TO_DEVICE, stream_));
   predictorCoarse_.copy_async(coarseDirection_, stream_);
   ADAPTIVE_LAUNCH(predictorAdaptiveHeunAtoms,
                   adaptiveGrid(atoms_ * ensembles_), stream_,
                   deviceTopology_, deviceRuntime_, kernels, timeStepSeconds,
                   atomDirection, predictorAtom_.data(),
                   initialAtomField_.data());
   ADAPTIVE_LAUNCH(predictorAdaptiveHeunBlocks,
                   adaptiveGrid(blocks_ * ensembles_), stream_,
                   deviceTopology_, deviceRuntime_, kernels, timeStepSeconds,
                   predictorCoarse_.data(), initialCoarseField_.data());
   phaseMetrics_.integrationLaunches += 2;
   ASSERT_GPU(GPU_GET_LAST_ERROR());
   finishPhase(phaseMetrics_.integrationMilliseconds);
   GpuAdaptiveUniformFftField predictorFftField{};
   const GpuAdaptiveUniformFftField* predictorFftFieldPtr = nullptr;
   if(basisResolvedFftEvaluator) {
      predictorFftField = basisResolvedFftEvaluator(atomDirection);
      predictorFftFieldPtr = &predictorFftField;
   }
   (void)evaluateHybrid(atomDirection, externalCoarseField, uniformFftDipoleField,
                        atomFieldScratch_.data(), predictorCoarseField_.data(),
                        predictorFftFieldPtr);
   ASSERT_GPU(GPU_EVENT_RECORD(phaseStart_, stream_));
   ADAPTIVE_LAUNCH(correctorAdaptiveHeunAtoms,
                   adaptiveGrid(atoms_ * ensembles_), stream_,
                   deviceTopology_, deviceRuntime_, kernels, timeStepSeconds,
                   atomDirection, predictorAtom_.data(),
                   initialAtomField_.data(), atomFieldScratch_.data());
   ADAPTIVE_LAUNCH(correctorAdaptiveHeunBlocks,
                   adaptiveGrid(blocks_ * ensembles_), stream_,
                   deviceTopology_, deviceRuntime_, kernels, timeStepSeconds,
                   predictorCoarse_.data(), initialCoarseField_.data(),
                   predictorCoarseField_.data());
   phaseMetrics_.integrationLaunches += 2;
   ASSERT_GPU(GPU_GET_LAST_ERROR());
   finishPhase(phaseMetrics_.integrationMilliseconds);
}

void GpuAdaptiveRuntime::synchronizeAtomicState(
   real* atomDirection,
   const GpuAdaptiveReconstructionPolicy& policy) {
   if(!ready_ || !kernelsReady_ || !atomDirection)
      throw std::logic_error("GPU adaptive reconstruction requires initialized state");
   AdaptiveKernelDevice kernels{
      bonds_, selectorEdges_, normalizationFloor_, magneticMomentSi_,
      gammaPerTs_, damping_, atomMoment_.data(), atomAnisotropyAxisCount_.data(),
      atomAnisotropyAxis_.data(), atomAnisotropyK1_.data(), atomAnisotropyK2_.data(),
      projectionBlock_.data(),
      projectionWeight_.data(), bondAtom_.data(), bondMatrix_.data(),
      selectorEdge_.data(), inverseBlockTranspose_.data(),
      exchangeStiffness_.data(), spiralization_.data(),
      anisotropyAxisCount_.data(), anisotropyAxis_.data(),
      anisotropyK1_.data(), anisotropyK2_.data(), ghostDirection_.data(),
      projectionNorm_.data(), atomFieldScratch_.data(),
      coarseFieldScratch_.data(), energyTerms_.data(), predictorAtom_.data()
   };
   ASSERT_GPU(GPU_EVENT_RECORD(phaseStart_, stream_));
#if defined(CUDA_V)
   prolongateAdaptiveGhosts<<<
      adaptiveGrid(atoms_ * ensembles_), adaptiveThreads, 0, stream_>>>(
      deviceTopology_, deviceRuntime_, kernels);
   commitAdaptiveGhosts<<<
      adaptiveGrid(atoms_ * ensembles_), adaptiveThreads, 0, stream_>>>(
      deviceTopology_, deviceRuntime_, kernels, atomDirection, policy);
#else
   hipLaunchKernelGGL(
      prolongateAdaptiveGhosts, dim3(adaptiveGrid(atoms_ * ensembles_)),
      dim3(adaptiveThreads), 0, stream_, deviceTopology_, deviceRuntime_,
      kernels);
   hipLaunchKernelGGL(
      commitAdaptiveGhosts, dim3(adaptiveGrid(atoms_ * ensembles_)),
      dim3(adaptiveThreads), 0, stream_, deviceTopology_, deviceRuntime_,
      kernels, atomDirection, policy);
#endif
   phaseMetrics_.integrationLaunches += 2;
   ASSERT_GPU(GPU_GET_LAST_ERROR());
   finishPhase(phaseMetrics_.integrationMilliseconds);
}

void GpuAdaptiveRuntime::recordFftMilliseconds(double elapsed) {
   if(!std::isfinite(elapsed) || elapsed < 0.0)
      throw std::invalid_argument("GPU adaptive FFT phase time must be finite and nonnegative");
   phaseMetrics_.fftMilliseconds += elapsed;
}

void GpuAdaptiveRuntime::recordStepWallMilliseconds(double elapsed) {
   if(!std::isfinite(elapsed) || elapsed < 0.0)
      throw std::invalid_argument("GPU adaptive step wall time must be finite and nonnegative");
   phaseMetrics_.stepWallMilliseconds += elapsed;
}

GpuAdaptiveWorkSnapshot GpuAdaptiveRuntime::downloadWorkSnapshot() {
   if(!ready_) throw std::logic_error("GPU adaptive runtime is not initialized");
   ASSERT_GPU(GPU_STREAM_SYNC(stream_));
   unsigned int counts[3] = {};
   TensorDataMovementTracker::add_d2h(sizeof(counts));
   ASSERT_GPU(GPU_MEMCPY(counts, workCounts_.data(), sizeof(counts), GPU_MEMCPY_DEVICE_TO_HOST));

   GpuAdaptiveWorkSnapshot result;
   result.activeAtoms.resize(counts[0]);
   result.activeBlocks.resize(counts[1]);
   result.interfaceAtoms.resize(counts[2]);
   result.atomisticBlockMask.resize(blocks_);
   result.coarseBlockMask.resize(blocks_);
   result.atomisticAtomMask.resize(atoms_);
   result.interfaceAtomMask.resize(atoms_);
   auto download = [](void* destination, const void* source, std::size_t bytes) {
      if(bytes == 0) return;
      TensorDataMovementTracker::add_d2h(bytes);
      ASSERT_GPU(GPU_MEMCPY(destination, source, bytes, GPU_MEMCPY_DEVICE_TO_HOST));
   };
   download(result.activeAtoms.data(), activeAtomList_.data(), counts[0] * sizeof(int));
   download(result.activeBlocks.data(), activeBlockList_.data(), counts[1] * sizeof(int));
   download(result.interfaceAtoms.data(), interfaceAtomList_.data(), counts[2] * sizeof(int));
   download(result.atomisticBlockMask.data(), atomisticBlockMask_.data(), blocks_);
   download(result.coarseBlockMask.data(), coarseBlockMask_.data(), blocks_);
   download(result.atomisticAtomMask.data(), atomisticAtomMask_.data(), atoms_);
   download(result.interfaceAtomMask.data(), interfaceAtomMask_.data(), atoms_);
   return result;
}

GpuAdaptiveDiagnosticSnapshot GpuAdaptiveRuntime::diagnosticSnapshot(
   const real* atomDirection) {
   if(!ready_ || !atomDirection)
      throw std::logic_error("GPU adaptive diagnostics require initialized state");
   ASSERT_GPU(GPU_STREAM_SYNC(stream_));
   GpuAdaptiveDiagnosticSnapshot result;
   result.blockState.resize(blocks_);
   result.stateAge.resize(blocks_);
   result.transitionEpoch.resize(blocks_);
   result.selectorScores.resize(deviceRuntime_.selectorCriteria * blocks_);
   result.polarizationRatio.resize(blocks_);
   const auto download = [](void* destination, const void* source,
                            std::size_t bytes) {
      if(bytes == 0) return;
      TensorDataMovementTracker::add_d2h(bytes);
      ASSERT_GPU(GPU_MEMCPY(destination, source, bytes,
                            GPU_MEMCPY_DEVICE_TO_HOST));
   };
   download(result.blockState.data(), blockState_.data(),
            blocks_ * sizeof(int));
   download(result.stateAge.data(), stateAge_.data(),
            blocks_ * sizeof(unsigned int));
   download(result.transitionEpoch.data(), transitionEpoch_.data(),
            blocks_ * sizeof(unsigned int));
   download(result.selectorScores.data(), selectorScores_.data(),
            result.selectorScores.size() * sizeof(real));
   download(result.polarizationRatio.data(), polarizationRatioBlock_.data(),
            result.polarizationRatio.size() * sizeof(real));
   std::vector<real> atomField(3 * atoms_ * ensembles_);
   std::vector<real> coarseField(3 * dynamicChannels_ * blocks_ * ensembles_);
   std::vector<real> direction(3 * atoms_ * ensembles_);
   std::vector<unsigned char> atomisticMask(atoms_);
   download(atomField.data(), atomFieldScratch_.data(),
            atomField.size() * sizeof(real));
   download(coarseField.data(), predictorCoarseField_.data(),
            coarseField.size() * sizeof(real));
   download(direction.data(), atomDirection,
            direction.size() * sizeof(real));
   download(atomisticMask.data(), atomisticAtomMask_.data(),
            atomisticMask.size() * sizeof(unsigned char));
   for(std::size_t ensemble = 0; ensemble < ensembles_; ++ensemble)
      for(std::size_t atom = 0; atom < atoms_; ++atom) {
         if(!atomisticMask[atom]) continue;
         for(int xyz = 0; xyz < 3; ++xyz) {
            const real value = atomField[
               xyz + 3 * (atom + atoms_ * ensemble)];
            result.atomFieldSumT += static_cast<double>(value);
            result.atomFieldNorm2T2 +=
               static_cast<double>(value) * static_cast<double>(value);
         }
      }
   for(const real value : coarseField) {
      result.coarseFieldSumT += static_cast<double>(value);
      result.coarseFieldNorm2T2 +=
         static_cast<double>(value) * static_cast<double>(value);
   }
   for(const real value : direction) {
      result.directionSum += static_cast<double>(value);
      result.directionNorm2 +=
         static_cast<double>(value) * static_cast<double>(value);
   }
   result.energy = lastEnergy_;
   return result;
}

void GpuAdaptiveRuntime::release() {
   if(streamCreated_) ASSERT_GPU(GPU_STREAM_SYNC(stream_));
   freeIfAllocated(acceptedBlockMask_);
   freeIfAllocated(polarizationUnsafeBlockMask_);
   freeIfAllocated(polarizationRatioBlock_);
   freeIfAllocated(energyTerms_);
   freeIfAllocated(predictorCoarseField_);
   freeIfAllocated(initialCoarseField_);
   freeIfAllocated(initialAtomField_);
   freeIfAllocated(predictorCoarse_);
   freeIfAllocated(predictorAtom_);
   freeIfAllocated(coarseFieldScratch_);
   freeIfAllocated(atomFieldScratch_);
   freeIfAllocated(projectionNorm_);
   freeIfAllocated(ghostDirection_);
   freeIfAllocated(anisotropyK2_);
   freeIfAllocated(anisotropyK1_);
   freeIfAllocated(anisotropyAxis_);
   freeIfAllocated(anisotropyAxisCount_);
   freeIfAllocated(spiralization_);
   freeIfAllocated(exchangeStiffness_);
   freeIfAllocated(inverseBlockTranspose_);
   freeIfAllocated(selectorEdge_);
   freeIfAllocated(bondMatrix_);
   freeIfAllocated(bondAtom_);
   freeIfAllocated(projectionWeight_);
   freeIfAllocated(projectionBlock_);
   freeIfAllocated(atomAnisotropyK2_);
   freeIfAllocated(atomAnisotropyK1_);
   freeIfAllocated(atomAnisotropyAxis_);
   freeIfAllocated(atomAnisotropyAxisCount_);
   freeIfAllocated(atomMoment_);
   freeIfAllocated(compactionScanLevels_);
   scanLevelItems_.clear();
   scanLevelOffset_.clear();
   freeIfAllocated(compactionScanB_);
   freeIfAllocated(compactionScanA_);
   freeIfAllocated(workCounts_);
   freeIfAllocated(interfaceAtomList_);
   freeIfAllocated(activeBlockList_);
   freeIfAllocated(activeAtomList_);
   freeIfAllocated(channelMomentSum_);
   freeIfAllocated(coarseField_);
   freeIfAllocated(coarseDirection_);
   freeIfAllocated(coarseMoment_);
   freeIfAllocated(selectorScores_);
   freeIfAllocated(interfaceAtomMask_);
   freeIfAllocated(atomisticAtomMask_);
   freeIfAllocated(coarseBlockMask_);
   freeIfAllocated(atomisticBlockMask_);
   freeIfAllocated(transitionEpoch_);
   freeIfAllocated(stateAge_);
   freeIfAllocated(dilatedState_);
   freeIfAllocated(pendingState_);
   freeIfAllocated(blockState_);
   freeIfAllocated(blockVolume_);
   freeIfAllocated(blockCenter_);
   freeIfAllocated(blockDynamicChannelPopulation_);
   freeIfAllocated(blockFftChannelPopulation_);
   freeIfAllocated(blockBasisPopulation_);
   freeIfAllocated(blockGridCoordinate_);
   freeIfAllocated(blockAtoms_);
   freeIfAllocated(blockAtomOffset_);
   freeIfAllocated(blockAtomCount_);
   freeIfAllocated(basisToFftChannel_);
   freeIfAllocated(basisToDynamicChannel_);
   freeIfAllocated(atomToFftGridIndex_);
   freeIfAllocated(atomToFftChannel_);
   freeIfAllocated(atomToDynamicChannel_);
   freeIfAllocated(atomToBasis_);
   freeIfAllocated(atomToBlock_);
   if(!stagedBlockState_.empty()) stagedBlockState_.FreeHost();
   if(startEventCreated_) {
      ASSERT_GPU(GPU_EVENT_DESTROY(updateStart_));
      startEventCreated_ = false;
   }
   if(endEventCreated_) {
      ASSERT_GPU(GPU_EVENT_DESTROY(updateEnd_));
      endEventCreated_ = false;
   }
   if(phaseEventsCreated_) {
      ASSERT_GPU(GPU_EVENT_DESTROY(phaseStart_));
      ASSERT_GPU(GPU_EVENT_DESTROY(phaseEnd_));
      phaseEventsCreated_ = false;
   }
   if(streamCreated_) {
      ASSERT_GPU(GPU_STREAM_DESTROY(stream_));
      streamCreated_ = false;
   }
   ready_ = false;
   kernelsReady_ = false;
   allocatedBytes_ = 0;
   atoms_ = blocks_ = dynamicChannels_ = ensembles_ = 0;
   bonds_ = selectorEdges_ = 0;
   normalizationFloor_ = magneticMomentSi_ = gammaPerTs_ = damping_ = real(0);
   metrics_ = {};
   phaseMetrics_ = {};
   lastEnergy_ = {};
   deviceTopology_ = {};
   deviceRuntime_ = {};
   convertedStaging_.clear();
}
