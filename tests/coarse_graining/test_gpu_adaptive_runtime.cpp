#include "gpuAdaptiveRuntime.hpp"
#include "measurement/memoryMeasurement.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

struct Fixture {
   static constexpr std::size_t atoms = 8;
   static constexpr std::size_t blocks = 4;
   static constexpr std::size_t basis = 2;
   static constexpr std::size_t channels = 2;
   static constexpr std::size_t ensembles = 3;
   int repetitionShape[3] = {4, 1, 1};
   int blockShape[3] = {1, 1, 1};
   int blockGrid[3] = {4, 1, 1};
   double cellVectors[9] = {1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0};
   double blockVectors[9] = {1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0};
   int atomToBlock[atoms] = {1, 1, 2, 2, 3, 3, 4, 4};
   int atomToBasis[atoms] = {1, 2, 1, 2, 1, 2, 1, 2};
   int atomToDynamic[atoms] = {1, 2, 1, 2, 1, 2, 1, 2};
   int atomToFft[atoms] = {1, 2, 1, 2, 1, 2, 1, 2};
   int atomToFftGrid[atoms] = {1, 2, 3, 4, 5, 6, 7, 8};
   int basisToDynamic[basis] = {1, 2};
   int basisToFft[basis] = {1, 2};
   int blockCount[blocks] = {2, 2, 2, 2};
   int blockOffset[blocks + 1] = {0, 2, 4, 6, 8};
   int blockAtoms[atoms] = {1, 2, 3, 4, 5, 6, 7, 8};
   int blockCoordinate[3 * blocks] = {0, 0, 0, 1, 0, 0, 2, 0, 0, 3, 0, 0};
   int basisPopulation[basis * blocks] = {1, 1, 1, 1, 1, 1, 1, 1};
   int fftPopulation[basis * blocks] = {1, 1, 1, 1, 1, 1, 1, 1};
   int dynamicPopulation[channels * blocks] = {1, 1, 1, 1, 1, 1, 1, 1};
   double center[3 * blocks] = {0.5, 0.0, 0.0, 1.5, 0.0, 0.0,
                                2.5, 0.0, 0.0, 3.5, 0.0, 0.0};
   double volume[blocks] = {1.0, 1.0, 1.0, 1.0};
   int state[blocks] = {0, 1, 2, 0};
   std::vector<double> moment;
   std::vector<double> direction;
   std::vector<double> field;
   std::vector<double> momentSum;

   Fixture()
      : moment(3 * channels * blocks * ensembles),
        direction(moment.size()),
        field(moment.size()),
        momentSum(channels * blocks * ensembles) {
      for(std::size_t ensemble = 0; ensemble < ensembles; ++ensemble) {
         for(std::size_t block = 0; block < blocks; ++block) {
            for(std::size_t channel = 0; channel < channels; ++channel) {
               const std::size_t scalar = channel + channels * (block + blocks * ensemble);
               momentSum[scalar] = 100.0 * ensemble + 10.0 * block + channel;
               for(std::size_t xyz = 0; xyz < 3; ++xyz) {
                  const std::size_t vector = xyz + 3 * scalar;
                  moment[vector] = 1000.0 * ensemble + 100.0 * block + 10.0 * channel + xyz;
                  direction[vector] = moment[vector] + 0.25;
                  field[vector] = moment[vector] - 0.5;
               }
            }
         }
      }
   }

   GpuAdaptiveTopologyInput topology() const {
      GpuAdaptiveTopologyInput t;
      t.geometryMode = 1;
      t.atoms = atoms;
      t.blocks = blocks;
      t.basis = basis;
      t.fftChannelsPerBlock = basis;
      t.fftGridChannels = basis * blocks;
      t.dynamicChannels = channels;
      t.ensembles = ensembles;
      t.repetitionShape = repetitionShape;
      t.blockShape = blockShape;
      t.blockGrid = blockGrid;
      t.cellVectors = cellVectors;
      t.blockVectors = blockVectors;
      t.atomToBlock = atomToBlock;
      t.atomToBasis = atomToBasis;
      t.atomToDynamicChannel = atomToDynamic;
      t.atomToFftChannel = atomToFft;
      t.atomToFftGridIndex = atomToFftGrid;
      t.basisToDynamicChannel = basisToDynamic;
      t.basisToFftChannel = basisToFft;
      t.blockAtomCount = blockCount;
      t.blockAtomOffset = blockOffset;
      t.blockAtoms = blockAtoms;
      t.blockGridCoordinate = blockCoordinate;
      t.blockBasisPopulation = basisPopulation;
      t.blockFftChannelPopulation = fftPopulation;
      t.blockDynamicChannelPopulation = dynamicPopulation;
      t.blockCenter = center;
      t.blockVolume = volume;
      return t;
   }

   GpuAdaptiveRuntimeInput runtime() const {
      GpuAdaptiveRuntimeInput r;
      r.blockState = state;
      r.selectorCriteria = 2;
      static double scores[2 * blocks] = {};
      r.selectorScores = scores;
      r.coarseMoment = moment.data();
      r.coarseDirection = direction.data();
      r.coarseField = field.data();
      r.channelMomentSum = momentSum.data();
      return r;
   }
};

void require(bool condition, const char* message) {
   if(!condition) throw std::runtime_error(message);
}

void expectEqual(const std::vector<int>& actual, std::initializer_list<int> expected,
                 const char* message) {
   require(actual == std::vector<int>(expected), message);
}

} // namespace

int main() {
   Fixture fixture;
   const auto topology = fixture.topology();
   const auto input = fixture.runtime();
   std::string diagnostic;

   require(GpuAdaptiveRuntime::estimateBytes({}, {}) == 0,
           "feature-off estimate is not zero");
   require(GpuAdaptiveRuntime::validate(topology, input, Fixture::atoms,
                                        Fixture::ensembles, diagnostic),
           "valid topology was rejected");
   require(!GpuAdaptiveRuntime::validate(topology, input, Fixture::atoms + 1,
                                         Fixture::ensembles, diagnostic),
           "mismatched Fortran/GPU atom count was accepted");
   require(diagnostic.find("do not match") != std::string::npos,
           "count mismatch diagnostic is not specific");

   auto badTopology = topology;
   int badOffset[Fixture::blocks + 1] = {0, 2, 4, 6, 7};
   badTopology.blockAtomOffset = badOffset;
   require(!GpuAdaptiveRuntime::validate(badTopology, input, Fixture::atoms,
                                         Fixture::ensembles, diagnostic),
           "invalid topology CSR was accepted");

   const auto baseline = TensorMemoryTracker::current_device();
   {
      GpuAdaptiveRuntime featureOff;
      require(!featureOff.ready() && featureOff.allocatedBytes() == 0,
              "feature-off runtime changed its inventory");
      require(TensorMemoryTracker::current_device() == baseline,
              "feature-off runtime allocated device storage");
   }
   GpuAdaptiveRuntime runtime;
   runtime.initialize(topology, input, Fixture::atoms, Fixture::ensembles);
   require(runtime.ready(), "runtime did not become ready");
   require(runtime.allocatedBytes() == GpuAdaptiveRuntime::estimateBytes(topology, input),
           "runtime allocation and preflight estimates differ");
   require(TensorMemoryTracker::current_device() - baseline ==
              static_cast<std::int64_t>(runtime.allocatedBytes()),
           "device memory tracker and runtime inventory differ");

   auto snapshot = runtime.downloadWorkSnapshot();
   expectEqual(snapshot.activeAtoms, {3, 4, 5, 6}, "initial active atom list is wrong");
   expectEqual(snapshot.activeBlocks, {1, 4}, "initial active block list is wrong");
   expectEqual(snapshot.interfaceAtoms, {3, 4}, "initial interface atom list is wrong");
   require(snapshot.atomisticBlockMask == std::vector<unsigned char>({0, 1, 1, 0}),
           "initial atomistic block mask is wrong");
   require(snapshot.coarseBlockMask == std::vector<unsigned char>({1, 0, 0, 1}),
           "initial coarse block mask is wrong");

   int changed[Fixture::blocks] = {2, 0, 1, 2};
   runtime.updateBlockState(changed, Fixture::blocks);
   snapshot = runtime.downloadWorkSnapshot();
   expectEqual(snapshot.activeAtoms, {1, 2, 5, 6, 7, 8},
               "rebuilt active atom list is wrong");
   expectEqual(snapshot.activeBlocks, {2}, "rebuilt active block list is wrong");
   expectEqual(snapshot.interfaceAtoms, {5, 6}, "rebuilt interface atom list is wrong");
   const auto metrics = runtime.compactionMetrics();
   require(metrics.hostSynchronizations == 1, "selector synchronization was not localized");
   require(metrics.blockBytesUploaded == Fixture::blocks * sizeof(int),
           "selector update transferred data beyond the block mask");
   require(metrics.elapsedMilliseconds >= 0.0, "selector synchronization was not timed");
   require(metrics.hostWaitMilliseconds >= 0.0 && metrics.deviceMilliseconds >= 0.0,
           "selector host wait and device interval were not measured");

   int invalid[Fixture::blocks] = {2, 3, 1, 2};
   bool rejected = false;
   try {
      runtime.updateBlockState(invalid, Fixture::blocks);
   } catch(const std::invalid_argument&) {
      rejected = true;
   }
   require(rejected, "invalid mask update was accepted");
   snapshot = runtime.downloadWorkSnapshot();
   expectEqual(snapshot.activeBlocks, {2}, "rejected mask update changed compact work");

   // Verify the explicit (3,channel,block,ensemble) storage contract at a
   // nonzero ensemble/channel/block index.
   const std::size_t scalar = 1 + Fixture::channels * (2 + Fixture::blocks * 2);
   const std::size_t vector = 2 + 3 * scalar;
   real stagedMoment = 0.0;
   ASSERT_GPU(GPU_MEMCPY(&stagedMoment,
                         runtime.deviceRuntime().coarseMoment + vector,
                         sizeof(stagedMoment), GPU_MEMCPY_DEVICE_TO_HOST));
   const double tolerance = sizeof(real) == sizeof(double) ? 1.0e-12 : 1.0e-5;
   require(std::abs(static_cast<double>(stagedMoment) - fixture.moment[vector]) < tolerance,
           "multiple-ensemble block/channel indexing is wrong");

   runtime.release();
   require(!runtime.ready(), "runtime remained ready after release");
   require(TensorMemoryTracker::current_device() == baseline,
           "release did not restore device memory accounting");
   std::cout << "CG-09 GPU adaptive runtime tests passed\n";
   return 0;
}
