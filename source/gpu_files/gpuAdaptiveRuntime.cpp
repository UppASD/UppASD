#include "gpuAdaptiveRuntime.hpp"

#include "base.hpp"
#include "measurement/memoryMeasurement.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <limits>
#include <numeric>
#include <stdexcept>

namespace {

constexpr int regularReplicatedCell = 1;
constexpr int coarseState = 0;
constexpr int bufferState = 1;
constexpr int fineState = 2;

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

__global__ void compactAdaptiveWork(GpuAdaptiveDeviceTopology topology,
                                    GpuAdaptiveDeviceRuntime runtime) {
   // CG-09 intentionally chooses the accepted simple device scan.  One thread
   // gives stable Fortran ordering without atomics or a host mask round trip.
   if(blockIdx.x != 0 || threadIdx.x != 0) return;

   unsigned int activeAtoms = 0;
   unsigned int activeBlocks = 0;
   unsigned int interfaceAtoms = 0;
   for(std::size_t block = 0; block < topology.blocks; ++block) {
      const int state = runtime.blockState[block];
      const bool atomistic = state != coarseState;
      const bool coarse = state == coarseState;
      runtime.atomisticBlockMask[block] = static_cast<unsigned char>(atomistic);
      runtime.coarseBlockMask[block] = static_cast<unsigned char>(coarse);
      if(coarse) runtime.activeBlockList[activeBlocks++] = static_cast<int>(block + 1);
   }
   for(std::size_t atom = 0; atom < topology.atoms; ++atom) {
      const int block = topology.atomToBlock[atom] - 1;
      const int state = runtime.blockState[block];
      const bool atomistic = state != coarseState;
      const bool interfaceAtom = state == bufferState;
      runtime.atomisticAtomMask[atom] = static_cast<unsigned char>(atomistic);
      runtime.interfaceAtomMask[atom] = static_cast<unsigned char>(interfaceAtom);
      if(atomistic) runtime.activeAtomList[activeAtoms++] = static_cast<int>(atom + 1);
      if(interfaceAtom) runtime.interfaceAtomList[interfaceAtoms++] = static_cast<int>(atom + 1);
   }
   runtime.workCounts[0] = activeAtoms;
   runtime.workCounts[1] = activeBlocks;
   runtime.workCounts[2] = interfaceAtoms;
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
      checkedAdd(total, 3, sizeof(unsigned int));
   if(!ok) throw std::overflow_error("GPU adaptive runtime memory estimate overflow");
   return total;
}

void GpuAdaptiveRuntime::initialize(const GpuAdaptiveTopologyInput& topology,
                                    const GpuAdaptiveRuntimeInput& runtime,
                                    std::size_t expectedAtoms,
                                    std::size_t expectedEnsembles) {
   if(ready_ || streamCreated_) throw std::logic_error("GPU adaptive runtime is already initialized");
   metrics_ = {};
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

void GpuAdaptiveRuntime::launchCompaction() {
#if defined(CUDA_V)
   compactAdaptiveWork<<<1, 1, 0, stream_>>>(deviceTopology_, deviceRuntime_);
#elif defined(HIP_V)
   hipLaunchKernelGGL(compactAdaptiveWork, dim3(1), dim3(1), 0, stream_,
                      deviceTopology_, deviceRuntime_);
#endif
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

void GpuAdaptiveRuntime::release() {
   if(streamCreated_) ASSERT_GPU(GPU_STREAM_SYNC(stream_));
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
   if(streamCreated_) {
      ASSERT_GPU(GPU_STREAM_DESTROY(stream_));
      streamCreated_ = false;
   }
   ready_ = false;
   allocatedBytes_ = 0;
   atoms_ = blocks_ = dynamicChannels_ = ensembles_ = 0;
   deviceTopology_ = {};
   deviceRuntime_ = {};
   convertedStaging_.clear();
}
