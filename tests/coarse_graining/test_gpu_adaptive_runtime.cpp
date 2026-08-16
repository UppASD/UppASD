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

struct KernelFixture {
   static constexpr std::size_t atoms = 8;
   static constexpr std::size_t blocks = 4;
   static constexpr std::size_t basis = 2;
   static constexpr std::size_t channels = 1;
   static constexpr std::size_t ensembles = 1;
   static constexpr std::size_t bonds = 8;
   int repetitionShape[3] = {4, 1, 1};
   int blockShape[3] = {1, 1, 1};
   int blockGrid[3] = {4, 1, 1};
   double cellVectors[9] = {1, 0, 0, 0, 1, 0, 0, 0, 1};
   double blockVectors[9] = {1, 0, 0, 0, 1, 0, 0, 0, 1};
   int atomToBlock[atoms] = {1, 1, 2, 2, 3, 3, 4, 4};
   int atomToBasis[atoms] = {1, 2, 1, 2, 1, 2, 1, 2};
   int atomToDynamic[atoms] = {1, 1, 1, 1, 1, 1, 1, 1};
   int atomToFft[atoms] = {1, 2, 1, 2, 1, 2, 1, 2};
   int atomToFftGrid[atoms] = {1, 2, 3, 4, 5, 6, 7, 8};
   int basisToDynamic[basis] = {1, 1};
   int basisToFft[basis] = {1, 2};
   int blockCount[blocks] = {2, 2, 2, 2};
   int blockOffset[blocks + 1] = {0, 2, 4, 6, 8};
   int blockAtoms[atoms] = {1, 2, 3, 4, 5, 6, 7, 8};
   int blockCoordinate[3 * blocks] = {0, 0, 0, 1, 0, 0, 2, 0, 0, 3, 0, 0};
   int basisPopulation[basis * blocks] = {1, 1, 1, 1, 1, 1, 1, 1};
   int fftPopulation[basis * blocks] = {1, 1, 1, 1, 1, 1, 1, 1};
   int dynamicPopulation[channels * blocks] = {2, 2, 2, 2};
   double center[3 * blocks] = {0.5, 0, 0, 1.5, 0, 0, 2.5, 0, 0, 3.5, 0, 0};
   double volume[blocks] = {1, 1, 1, 1};
   int state[blocks] = {0, 1, 2, 0};
   double scores[blocks] = {};
   double atomMoment[atoms] = {1, 1, 1, 1, 1, 1, 1, 1};
   int atomAxisCount[atoms] = {};
   double atomAxis[6 * atoms] = {};
   double atomK1[2 * atoms] = {};
   double atomK2[2 * atoms] = {};
   int projectionBlock[8 * atoms] = {};
   double projectionWeight[8 * atoms] = {};
   int bondAtom[2 * bonds] = {
      1, 2, 3, 4, 5, 6, 7, 8,
      2, 3, 4, 5, 6, 7, 8, 1
   };
   double bondMatrix[9 * bonds] = {};
   int selectorEdge[2 * bonds] = {
      1, 2, 3, 4, 5, 6, 7, 8,
      2, 3, 4, 5, 6, 7, 8, 1
   };
   double inverseBlockTranspose[9] = {1, 0, 0, 0, 1, 0, 0, 0, 1};
   double exchange[9] = {};
   double spiralization[9] = {};
   int axisCount[blocks] = {1, 1, 1, 1};
   double axis[6 * blocks] = {};
   double k1[2 * blocks] = {};
   double k2[2 * blocks] = {};
   std::vector<double> coarseMoment =
      std::vector<double>(3 * channels * blocks * ensembles, 0.0);
   std::vector<double> coarseDirection =
      std::vector<double>(3 * channels * blocks * ensembles, 0.0);
   std::vector<double> coarseField =
      std::vector<double>(3 * channels * blocks * ensembles, 0.0);
   std::vector<double> momentSum =
      std::vector<double>(channels * blocks * ensembles, 2.0);

   KernelFixture() {
      for(std::size_t atom = 0; atom < atoms; ++atom) {
         projectionBlock[8 * atom] = atomToBlock[atom];
         projectionWeight[8 * atom] = 1.0;
         for(int corner = 1; corner < 8; ++corner)
            projectionBlock[corner + 8 * atom] = atomToBlock[atom];
      }
      for(std::size_t bond = 0; bond < bonds; ++bond)
         for(int xyz = 0; xyz < 3; ++xyz)
            bondMatrix[xyz + 3 * xyz + 9 * bond] = 0.4;
      for(std::size_t block = 0; block < blocks; ++block) {
         axis[2 + 3 * (2 * block)] = 1.0;
         k1[2 * block] = 2.0;
         k2[2 * block] = 0.5;
         coarseMoment[2 + 3 * block] = 2.0;
         coarseDirection[2 + 3 * block] = 1.0;
      }
      atomAxisCount[2] = 1;
      atomAxis[2 + 3 * (2 * 2)] = 1.0;
      atomK1[2 * 2] = 0.25;
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
      r.selectorCriteria = 1;
      r.selectorScores = scores;
      r.coarseMoment = coarseMoment.data();
      r.coarseDirection = coarseDirection.data();
      r.coarseField = coarseField.data();
      r.channelMomentSum = momentSum.data();
      r.kernels.atomMoment = atomMoment;
      r.kernels.atomAnisotropyAxisCount = atomAxisCount;
      r.kernels.atomAnisotropyAxis = atomAxis;
      r.kernels.atomAnisotropyK1 = atomK1;
      r.kernels.atomAnisotropyK2 = atomK2;
      r.kernels.projectionBlock = projectionBlock;
      r.kernels.projectionWeight = projectionWeight;
      r.kernels.bonds = bonds;
      r.kernels.bondAtom = bondAtom;
      r.kernels.bondMatrix = bondMatrix;
      r.kernels.selectorEdges = bonds;
      r.kernels.selectorEdge = selectorEdge;
      r.kernels.inverseBlockTranspose = inverseBlockTranspose;
      r.kernels.exchangeStiffness = exchange;
      r.kernels.spiralization = spiralization;
      r.kernels.anisotropyAxisCount = axisCount;
      r.kernels.anisotropyAxis = axis;
      r.kernels.anisotropyK1 = k1;
      r.kernels.anisotropyK2 = k2;
      r.kernels.magneticMomentSi = 1.0;
      r.kernels.gammaPerTs = 2.0;
      r.kernels.damping = 0.1;
      return r;
   }
};

void require(bool condition, const char* message) {
   if(!condition) throw std::runtime_error(message);
}

void require(bool condition, const std::string& message) {
   if(!condition) throw std::runtime_error(message);
}
void expectEqual(const std::vector<int>& actual, std::initializer_list<int> expected,
                 const char* message) {
   require(actual == std::vector<int>(expected), message);
}
// RCG-08: the compaction fixture's reference lists are computed, not literal,
// so they arrive as a vector rather than an initializer list.
void expectEqual(const std::vector<int>& actual, const std::vector<int>& expected,
                 const char* message) {
   require(actual == expected, message);
}

std::vector<real> deviceVector(const std::vector<double>& source) {
   std::vector<real> result(source.size());
   std::transform(source.begin(), source.end(), result.begin(),
                  [](double value) { return static_cast<real>(value); });
   return result;
}

void upload(real* destination, const std::vector<real>& source) {
   ASSERT_GPU(GPU_MEMCPY(destination, source.data(), source.size() * sizeof(real),
                         GPU_MEMCPY_HOST_TO_DEVICE));
}

std::vector<real> download(const real* source, std::size_t count) {
   std::vector<real> result(count);
   ASSERT_GPU(GPU_MEMCPY(result.data(), source, count * sizeof(real),
                         GPU_MEMCPY_DEVICE_TO_HOST));
   return result;
}

// RCG-09C: the compact live-bond list and the ghost shell it induces.
//
// KernelFixture's state {coarse, buffer, fine, coarse} over a periodic chain
// exercises every ownership case the accepted hybrid rule distinguishes, so
// the expected lists below are written out rather than recomputed with the
// same predicate the implementation uses:
//
//   atoms 1,2 -> block 1 (coarse)   atoms 5,6 -> block 3 (fine)
//   atoms 3,4 -> block 2 (buffer)   atoms 7,8 -> block 4 (coarse)
//
//   bond 1 (1,2) coarse-coarse   -> owned by the continuum operator, excluded
//   bond 2 (3,4) buffer-buffer   -> atomistic
//   bond 3 (5,6) fine-fine       -> atomistic
//   bond 4 (7,8) coarse-coarse   -> excluded
//   bond 5 (2,3) coarse-buffer   -> atomistic, ghost endpoint atom 2
//   bond 6 (4,5) buffer-fine     -> atomistic
//   bond 7 (6,7) fine-coarse     -> atomistic, ghost endpoint atom 7
//   bond 8 (8,1) coarse-coarse   -> excluded
void testLiveBondCompaction() {
   KernelFixture fixture;
   const auto topology = fixture.topology();
   const auto input = fixture.runtime();
   GpuAdaptiveRuntime runtime;
   runtime.initialize(topology, input, KernelFixture::atoms,
                      KernelFixture::ensembles);

   const std::vector<int> expectedBonds{2, 3, 5, 6, 7};
   const std::vector<int> expectedGhosts{2, 7};
   const auto snapshot = runtime.downloadWorkSnapshot();
   expectEqual(snapshot.activeBonds, expectedBonds,
               "live-bond compaction changed the hybrid bond ownership rule");
   expectEqual(snapshot.ghostAtoms, expectedGhosts,
               "ghost shell does not contain exactly the coarse endpoints of "
               "live bonds");
   require(runtime.liveBondCount() == expectedBonds.size(),
           "host live-bond mirror disagrees with the device list");
   require(runtime.ghostAtomCount() == expectedGhosts.size(),
           "host ghost-shell mirror disagrees with the device list");
   require(runtime.liveBondCount() < runtime.bondCount(),
           "fixture does not actually exclude any coarse-coarse bond, so it "
           "cannot demonstrate compaction");

   // Ordering is a function of state alone: the list is strictly ascending,
   // and rebuilding from the identical state reproduces it exactly.  Field
   // scatter and the block-partial energy reduction both depend on that.
   for(std::size_t slot = 1; slot < snapshot.activeBonds.size(); ++slot)
      require(snapshot.activeBonds[slot - 1] < snapshot.activeBonds[slot],
              "live-bond list is not in ascending bond order");
   runtime.updateBlockState(fixture.state, KernelFixture::blocks);
   const auto rebuilt = runtime.downloadWorkSnapshot();
   expectEqual(rebuilt.activeBonds, expectedBonds,
               "recompaction from an unchanged state produced a different "
               "live-bond ordering");
   expectEqual(rebuilt.ghostAtoms, expectedGhosts,
               "recompaction from an unchanged state produced a different "
               "ghost shell");

   // All-fine is the limit the RCG-09A parity gates run in: every bond stays
   // live and the compact list must be the identity permutation, so the
   // all-fine launch is the same work in the same order it was before.
   std::vector<int> allFine(KernelFixture::blocks, 2);
   runtime.updateBlockState(allFine.data(), allFine.size());
   const auto fine = runtime.downloadWorkSnapshot();
   require(fine.activeBonds.size() == KernelFixture::bonds,
           "all-fine state dropped a bond from the live list");
   for(std::size_t bond = 0; bond < KernelFixture::bonds; ++bond)
      require(fine.activeBonds[bond] == static_cast<int>(bond + 1),
              "all-fine live-bond list is not the identity permutation");
   require(fine.ghostAtoms.empty(),
           "all-fine state produced a ghost shell");

   // All-coarse is the other limit: the continuum operator owns every bond,
   // so the atomistic evaluator has nothing to launch at all.
   std::vector<int> allCoarse(KernelFixture::blocks, 0);
   runtime.updateBlockState(allCoarse.data(), allCoarse.size());
   const auto coarse = runtime.downloadWorkSnapshot();
   require(coarse.activeBonds.empty() && coarse.ghostAtoms.empty() &&
           runtime.liveBondCount() == 0 && runtime.activeAtomCount() == 0,
           "all-coarse state left atomistic work in the live lists");

   // CGP-02: an empty ghost shell must issue no prolongation kernel at all,
   // at either call site -- not merely a cheap one-block no-op launch.
   const std::size_t atomVectors = 3 * KernelFixture::atoms;
   const std::size_t coarseVectors = 3 * KernelFixture::blocks;
   GpuTensor<real, 1> atomDirection, atomField, coarseField;
   atomDirection.Allocate(static_cast<index_t>(atomVectors));
   atomField.Allocate(static_cast<index_t>(atomVectors));
   coarseField.Allocate(static_cast<index_t>(coarseVectors));
   std::vector<double> hostAtom(atomVectors, 0.0);
   for(std::size_t atom = 0; atom < KernelFixture::atoms; ++atom)
      hostAtom[2 + 3 * atom] = 1.0;
   upload(atomDirection.data(), deviceVector(hostAtom));

   const auto beforeHybrid = runtime.phaseMetrics();
   (void)runtime.evaluateHybrid(atomDirection.data(), nullptr, nullptr,
                                atomField.data(), coarseField.data());
   const auto afterHybrid = runtime.phaseMetrics();
   // clearAdaptiveInterface always launches (it zeroes coarseFieldScratch
   // unconditionally); the ghost-shell-gated launches -- the compact
   // prolongation here, and restrictAdaptiveInterface further down the same
   // call -- must contribute nothing on top of it.
   require(afterHybrid.interfaceLaunches == beforeHybrid.interfaceLaunches + 1,
           "an empty ghost shell issued more than the unconditional "
           "clearAdaptiveInterface launch from evaluateHybrid");

   const auto beforeSync = runtime.phaseMetrics();
   runtime.synchronizeAtomicState(atomDirection.data(), {});
   const auto afterSync = runtime.phaseMetrics();
   require(afterSync.integrationLaunches == beforeSync.integrationLaunches + 1,
           "an empty ghost shell issued more than the required commit launch "
           "from synchronizeAtomicState under the Aligned policy");

   coarseField.Free();
   atomField.Free();
   atomDirection.Free();
   runtime.release();
}

void testKernelParityAndWorkflow() {
   KernelFixture fixture;
   const auto topology = fixture.topology();
   const auto input = fixture.runtime();
   std::string diagnostic;
   require(GpuAdaptiveRuntime::validate(topology, input, KernelFixture::atoms,
                                        KernelFixture::ensembles, diagnostic),
           "valid CG-10 operator inventory was rejected");

   const auto baseline = TensorMemoryTracker::current_device();
   GpuAdaptiveRuntime runtime;
   runtime.initialize(topology, input, KernelFixture::atoms,
                      KernelFixture::ensembles);
   require(runtime.kernelsReady(), "CG-10 kernels were not published ready");
   require(TensorMemoryTracker::current_device() - baseline ==
              static_cast<std::int64_t>(GpuAdaptiveRuntime::estimateBytes(topology, input)),
           "CG-10 scratch or construction storage bypassed memory preflight");

   const std::size_t atomVectors = 3 * KernelFixture::atoms;
   const std::size_t coarseVectors = 3 * KernelFixture::blocks;
   GpuTensor<real, 1> atomDirection, atomField, coarseField, externalField, dipoleField;
   atomDirection.Allocate(static_cast<index_t>(atomVectors));
   atomField.Allocate(static_cast<index_t>(atomVectors));
   coarseField.Allocate(static_cast<index_t>(coarseVectors));
   externalField.Allocate(static_cast<index_t>(coarseVectors));
   dipoleField.Allocate(static_cast<index_t>(coarseVectors));
   std::vector<double> hostAtom(atomVectors, 0.0);
   std::vector<double> hostExternal(coarseVectors, 0.0);
   std::vector<double> hostDipole(coarseVectors, 0.0);
   for(std::size_t atom = 0; atom < KernelFixture::atoms; ++atom)
      hostAtom[2 + 3 * atom] = 1.0;
   for(std::size_t block = 0; block < KernelFixture::blocks; ++block) {
      hostExternal[3 * block] = 0.25;
      hostDipole[2 + 3 * block] = 0.1;
   }
   upload(atomDirection.data(), deviceVector(hostAtom));
   upload(externalField.data(), deviceVector(hostExternal));
   upload(dipoleField.data(), deviceVector(hostDipole));

   runtime.restrictMoments(atomDirection.data());
   const auto restricted = download(runtime.deviceRuntime().coarseMoment, coarseVectors);
   for(std::size_t block = 0; block < KernelFixture::blocks; ++block) {
      require(std::abs(static_cast<double>(restricted[2 + 3 * block]) - 2.0) < 1.0e-12,
              "GPU restriction does not match the CPU moment sum");
   }

   const auto energy = runtime.evaluateHybrid(
      atomDirection.data(), externalField.data(), dipoleField.data(),
      atomField.data(), coarseField.data());
   const double tolerance = sizeof(real) == sizeof(double) ? 2.0e-12 : 2.0e-5;
   require(std::abs(energy.atomisticBilinearJ + 2.0) < tolerance,
           "GPU atomistic unique-pair energy violates the CPU sign contract");
   require(std::abs(energy.atomisticOnsiteJ - 0.25) < tolerance,
           "GPU atomistic anisotropy energy violates the CPU onsite contract");
   require(std::abs(energy.coarseExchangeJ) < tolerance &&
           std::abs(energy.coarseSpiralizationJ) < tolerance,
           "uniform state produced a coarse gradient energy");
   require(std::abs(energy.coarseAnisotropyJ - 5.0) < tolerance,
           "GPU anisotropy energy violates the CPU tensor-index contract");
   require(std::abs(energy.coarseExternalJ) < tolerance,
           "orthogonal external field produced energy");
   require(std::abs(energy.dipoleJ + 0.4) < tolerance,
           "all-grid FFT dipole energy was adaptively masked");
   require(std::abs(energy.totalJ - 2.85) < 5.0 * tolerance,
           "GPU hybrid total energy does not equal the CPU reference");

   // CGP-00: downloadEnergyPartialsSnapshot() must reproduce the same
   // per-term totals evaluateHybrid() just computed -- it downloads the
   // identical per-block FP64 partials from that same call, sliced to the
   // blocks that call actually wrote, so an ordered sum over each term's
   // partials must equal that term's entry in `energy`. This doubles as real
   // production-kernel evidence for CGP-00 Part C: the per-block
   // magnitude/sign distribution the actual atomistic exchange+DMI-folded
   // bond, onsite anisotropy, coarse exchange/spiralization/anisotropy/
   // external, and dipole kernels produce, not a synthetic proxy.
   const auto partials = runtime.downloadEnergyPartialsSnapshot(true);
   auto orderedSum = [](const std::vector<double>& values) {
      double total = 0.0;
      for(const double value : values) total += value;
      return total;
   };
   require(std::abs(orderedSum(partials.atomisticBilinearJ) - energy.atomisticBilinearJ) < tolerance,
           "downloaded atomistic bilinear partials do not reproduce evaluateHybrid's own total");
   require(std::abs(orderedSum(partials.atomisticOnsiteJ) - energy.atomisticOnsiteJ) < tolerance,
           "downloaded atomistic onsite partials do not reproduce evaluateHybrid's own total");
   require(std::abs(orderedSum(partials.coarseExchangeJ) - energy.coarseExchangeJ) < tolerance,
           "downloaded coarse exchange partials do not reproduce evaluateHybrid's own total");
   require(std::abs(orderedSum(partials.coarseSpiralizationJ) - energy.coarseSpiralizationJ) < tolerance,
           "downloaded coarse spiralization partials do not reproduce evaluateHybrid's own total");
   require(std::abs(orderedSum(partials.coarseAnisotropyJ) - energy.coarseAnisotropyJ) < tolerance,
           "downloaded coarse anisotropy partials do not reproduce evaluateHybrid's own total");
   require(std::abs(orderedSum(partials.coarseExternalJ) - energy.coarseExternalJ) < tolerance,
           "downloaded coarse external partials do not reproduce evaluateHybrid's own total");
   require(std::abs(orderedSum(partials.dipoleJ) - energy.dipoleJ) < tolerance,
           "downloaded dipole partials do not reproduce evaluateHybrid's own total");

   const auto fineField = download(atomField.data(), atomVectors);
   for(std::size_t atom = 2; atom <= 5; ++atom)
      require(std::abs(static_cast<double>(fineField[2 + 3 * atom]) -
                       (atom == 2 ? 0.4 : 0.9)) < tolerance,
              "GPU atomistic/interface field does not contain the FFT dipole exactly once");
   const auto blockField = download(coarseField.data(), coarseVectors);
   for(const std::size_t block : {std::size_t(0), std::size_t(3)}) {
      require(std::abs(static_cast<double>(blockField[3 * block]) - 0.25) < tolerance,
              "GPU external-field derivative sign differs from CPU");
      require(std::abs(static_cast<double>(blockField[2 + 3 * block]) + 1.9) <
                 5.0 * tolerance,
              "GPU anisotropy/dipole field differs from CPU");
      // Directional derivative for m(theta)=(sin theta,0,cos theta):
      // dE/dtheta|0 = -mu_block Bx = -2*0.25.
      require(std::abs(-2.0 * static_cast<double>(blockField[3 * block]) + 0.5) <
                 tolerance,
              "GPU field fails the CPU energy directional-derivative fixture");
   }

   GpuTensor<real, 1> basisResolvedDipole;
   basisResolvedDipole.Allocate(
      3 * KernelFixture::basis * KernelFixture::blocks);
   std::vector<double> hostBasisDipole(
      3 * KernelFixture::basis * KernelFixture::blocks, 0.0);
   for(std::size_t block = 0; block < KernelFixture::blocks; ++block) {
      hostBasisDipole[block + KernelFixture::blocks * 2] = 0.1;
      hostBasisDipole[
         block + KernelFixture::blocks * (2 + 3)] = 0.3;
   }
   upload(basisResolvedDipole.data(), deviceVector(hostBasisDipole));
   GpuAdaptiveUniformFftField basisView{};
   basisView.paddedField = basisResolvedDipole.data();
   basisView.activeN1 = KernelFixture::blocks;
   basisView.activeN2 = 1;
   basisView.activeN3 = 1;
   basisView.fftN1 = KernelFixture::blocks;
   basisView.fftN2 = 1;
   basisView.fftCells = KernelFixture::blocks;
   basisView.basis = KernelFixture::basis;
   basisView.prefactorT = real(2);
   const auto basisEnergy = runtime.evaluateHybrid(
      atomDirection.data(), externalField.data(), nullptr,
      atomField.data(), coarseField.data(), &basisView);
   require(std::abs(basisEnergy.dipoleJ + 1.6) < 5.0 * tolerance,
           "basis-resolved FFT dipole energy was not counted exactly once");
   const auto basisAtomField = download(atomField.data(), atomVectors);
   for(std::size_t atom = 2; atom <= 5; ++atom) {
      const double expected = atom == 2 ? 0.5 : (atom % 2 == 0 ? 1.0 : 1.4);
      require(std::abs(
                 static_cast<double>(basisAtomField[2 + 3 * atom]) -
                 expected) < 5.0 * tolerance,
              "basis-resolved FFT field did not enter an atomistic/interface equation exactly once");
   }
   const auto basisCoarseField = download(coarseField.data(), coarseVectors);
   for(const std::size_t block : {std::size_t(0), std::size_t(3)})
      require(std::abs(
                 static_cast<double>(basisCoarseField[2 + 3 * block]) +
                 1.6) < 5.0 * tolerance,
              "basis-resolved FFT field did not enter the coarse equation exactly once");
   basisResolvedDipole.Free();

   hostAtom[0] = 1.0;
   hostAtom[2] = 0.0;
   upload(atomDirection.data(), deviceVector(hostAtom));
   runtime.evaluateSelectorScores(atomDirection.data());
   const auto selector = download(runtime.deviceRuntime().selectorScores,
                                  KernelFixture::blocks);
   require(std::abs(static_cast<double>(selector[0]) - 1.0) < tolerance &&
           std::abs(static_cast<double>(selector[3]) - 1.0) < tolerance &&
           std::abs(static_cast<double>(selector[1])) < tolerance &&
           std::abs(static_cast<double>(selector[2])) < tolerance,
           "GPU selector scores differ from the CPU maximum-misalignment reference");
   GpuAdaptiveSelectorPolicy selectorPolicy;
   selectorPolicy.refineThreshold = real(0.5);
   selectorPolicy.coarsenThreshold = real(0.1);
   auto dilatedPolicy = selectorPolicy;
   // RCG-05D: this fixture's block grid is 1D ({4,1,1}), so only axis 0 (x)
   // can move a block from coarse to buffer; y/z are deliberately left at 0
   // (rather than mirroring the old isotropic scalar into all three) to
   // confirm the per-axis field is actually read component-wise by the real
   // CUDA/HIP kernel launch, not broadcast from a single value.
   dilatedPolicy.bufferDilationBlocks[0] = 1;
   dilatedPolicy.bufferDilationBlocks[1] = 0;
   dilatedPolicy.bufferDilationBlocks[2] = 0;
   runtime.proposeSelectorState(dilatedPolicy);
   std::vector<int> pending(KernelFixture::blocks);
   ASSERT_GPU(GPU_MEMCPY(pending.data(), runtime.deviceRuntime().pendingState,
                         pending.size() * sizeof(int), GPU_MEMCPY_DEVICE_TO_HOST));
   require(pending == std::vector<int>({2, 1, 1, 2}),
           "GPU selector buffer dilation differs from the periodic CPU state map");
   runtime.proposeSelectorState(selectorPolicy);
   bool stageRejected = false;
   try {
      runtime.publishProposedState(atomDirection.data(), {}, false);
   } catch(const std::invalid_argument&) {
      stageRejected = true;
   }
   require(stageRejected, "GPU state transition was allowed inside an integration stage");
   GpuTensor<unsigned char, 1> accepted;
   accepted.Allocate(KernelFixture::blocks);
   const std::vector<unsigned char> hostAccepted = {1, 1, 0, 0};
   ASSERT_GPU(GPU_MEMCPY(accepted.data(), hostAccepted.data(), hostAccepted.size(),
                         GPU_MEMCPY_HOST_TO_DEVICE));
   runtime.publishProposedState(atomDirection.data(), {}, true, accepted.data());
   const auto snapshot = runtime.downloadWorkSnapshot();
   expectEqual(snapshot.activeAtoms, {1, 2, 5, 6},
               "accepted/rejected GPU transitions produced the wrong active atoms");
   expectEqual(snapshot.activeBlocks, {2, 4},
               "accepted/rejected GPU transitions produced the wrong coarse blocks");
   const auto reconstructed = download(atomDirection.data(), atomVectors);
   require(std::abs(static_cast<double>(reconstructed[2]) - 1.0) < tolerance &&
           std::abs(static_cast<double>(reconstructed[5]) - 1.0) < tolerance,
           "aligned GPU reconstruction is not exact");

   // The constrained-cone path uses the CPU tuple seed
   // (global,block,channel,ensemble,prospective epoch).  Two fresh-equivalent
   // runtimes must therefore publish bitwise-identical directions.
   real reducedResultant = real(1.8);
   ASSERT_GPU(GPU_MEMCPY(runtime.deviceRuntime().coarseMoment + 2 + 3 * 3,
                         &reducedResultant, sizeof(real), GPU_MEMCPY_HOST_TO_DEVICE));
   const std::vector<unsigned char> acceptLast = {0, 0, 0, 1};
   ASSERT_GPU(GPU_MEMCPY(accepted.data(), acceptLast.data(), acceptLast.size(),
                         GPU_MEMCPY_HOST_TO_DEVICE));
   runtime.proposeSelectorState(selectorPolicy);
   GpuAdaptiveReconstructionPolicy cone;
   cone.scheme = GpuAdaptiveReconstruction::ConstrainedCone;
   cone.coneAngleRadians = real(0.8);
   cone.globalSeed = 99173;
   runtime.publishProposedState(atomDirection.data(), cone, true, accepted.data());
   const auto firstCone = download(atomDirection.data(), atomVectors);
   require(std::abs(static_cast<double>(firstCone[0 + 3 * 6] +
                                        firstCone[0 + 3 * 7])) < 2.0e-10 &&
           std::abs(static_cast<double>(firstCone[1 + 3 * 6] +
                                        firstCone[1 + 3 * 7])) < 2.0e-10 &&
           std::abs(static_cast<double>(firstCone[2 + 3 * 6] +
                                        firstCone[2 + 3 * 7]) - 1.8) <
              (sizeof(real) == sizeof(double) ? 2.0e-10 : 2.0e-4),
           "constrained-cone GPU reconstruction missed the requested resultant");

   {
      GpuAdaptiveRuntime repeat;
      repeat.initialize(topology, input, KernelFixture::atoms,
                        KernelFixture::ensembles);
      GpuTensor<real, 1> repeatDirection;
      GpuTensor<unsigned char, 1> repeatAccepted;
      repeatDirection.Allocate(static_cast<index_t>(atomVectors));
      repeatAccepted.Allocate(KernelFixture::blocks);
      upload(repeatDirection.data(), deviceVector(hostAtom));
      repeat.evaluateSelectorScores(repeatDirection.data());
      repeat.proposeSelectorState(selectorPolicy);
      ASSERT_GPU(GPU_MEMCPY(repeat.deviceRuntime().coarseMoment + 2 + 3 * 3,
                            &reducedResultant, sizeof(real),
                            GPU_MEMCPY_HOST_TO_DEVICE));
      ASSERT_GPU(GPU_MEMCPY(repeatAccepted.data(), acceptLast.data(),
                            acceptLast.size(), GPU_MEMCPY_HOST_TO_DEVICE));
      repeat.publishProposedState(repeatDirection.data(), cone, true,
                                  repeatAccepted.data());
      const auto secondCone = download(repeatDirection.data(), atomVectors);
      for(std::size_t atom = 6; atom < 8; ++atom)
         for(int xyz = 0; xyz < 3; ++xyz)
            require(firstCone[xyz + 3 * atom] == secondCone[xyz + 3 * atom],
                    "tuple-seeded cone reconstruction is not deterministic");
      repeatAccepted.Free();
      repeatDirection.Free();
      repeat.release();
   }

   const auto beforeStep = runtime.downloadWorkSnapshot();
   (void)runtime.evaluateHybrid(atomDirection.data(), externalField.data(),
                                dipoleField.data(), atomField.data(),
                                coarseField.data());
   runtime.prepareCoarsePredictor(real(1.0e-3), coarseField.data());
   runtime.synchronize();
   (void)runtime.evaluateHybrid(atomDirection.data(), externalField.data(),
                                dipoleField.data(), atomField.data(),
                                coarseField.data());
   runtime.correctCoarse(real(1.0e-3), coarseField.data());
   runtime.synchronize();
   const auto afterStep = runtime.downloadWorkSnapshot();
   require(beforeStep.activeAtoms == afterStep.activeAtoms &&
           beforeStep.activeBlocks == afterStep.activeBlocks,
           "coarse predictor/corrector changed adaptive ownership inside a stage");

   runtime.recordFftMilliseconds(0.125);
   const auto phases = runtime.phaseMetrics();
   require(phases.atomisticMilliseconds >= 0.0 &&
           phases.coarseMilliseconds >= 0.0 &&
           phases.interfaceMilliseconds >= 0.0 &&
           phases.selectorMilliseconds >= 0.0 &&
           phases.compactionMilliseconds >= 0.0 &&
           phases.integrationMilliseconds >= 0.0 &&
           std::abs(phases.fftMilliseconds - 0.125) < 1.0e-15,
           "CG-10 phases are not independently accounted");

   accepted.Free();
   dipoleField.Free();
   externalField.Free();
   coarseField.Free();
   atomField.Free();
   atomDirection.Free();
   runtime.release();
   require(TensorMemoryTracker::current_device() == baseline,
           "CG-10 runtime cleanup did not restore tracked device memory");
}

// CGP-01: evaluateField() must be a genuine field-only sibling of
// evaluateHybrid() -- same field arithmetic (evaluateHybridImpl<MeasureEnergy>
// is the single shared implementation), but zero energy contribution
// arithmetic, zero block reduction, zero finalization launch, and zero D2H
// energy readback.  This is both the negative control the CGP-01 design
// requirement asks for (energyEvaluationCount() must stay exactly zero across
// field-only calls) and the field-parity check (field-only output must equal
// the field evaluateHybrid() produces for the identical inputs).
void testFieldOnlyEvaluation() {
   KernelFixture fixture;
   const auto topology = fixture.topology();
   const auto input = fixture.runtime();
   const auto baseline = TensorMemoryTracker::current_device();
   GpuAdaptiveRuntime runtime;
   runtime.initialize(topology, input, KernelFixture::atoms,
                      KernelFixture::ensembles);
   require(runtime.kernelsReady(), "CG-10 kernels were not published ready");

   const std::size_t atomVectors = 3 * KernelFixture::atoms;
   const std::size_t coarseVectors = 3 * KernelFixture::blocks;
   GpuTensor<real, 1> atomDirection, atomFieldOnly, atomWithEnergy,
      coarseFieldOnly, coarseWithEnergy, externalField, dipoleField;
   atomDirection.Allocate(static_cast<index_t>(atomVectors));
   atomFieldOnly.Allocate(static_cast<index_t>(atomVectors));
   atomWithEnergy.Allocate(static_cast<index_t>(atomVectors));
   coarseFieldOnly.Allocate(static_cast<index_t>(coarseVectors));
   coarseWithEnergy.Allocate(static_cast<index_t>(coarseVectors));
   externalField.Allocate(static_cast<index_t>(coarseVectors));
   dipoleField.Allocate(static_cast<index_t>(coarseVectors));
   std::vector<double> hostAtom(atomVectors, 0.0);
   std::vector<double> hostExternal(coarseVectors, 0.0);
   std::vector<double> hostDipole(coarseVectors, 0.0);
   for(std::size_t atom = 0; atom < KernelFixture::atoms; ++atom)
      hostAtom[2 + 3 * atom] = 1.0;
   for(std::size_t block = 0; block < KernelFixture::blocks; ++block) {
      hostExternal[3 * block] = 0.25;
      hostDipole[2 + 3 * block] = 0.1;
   }
   upload(atomDirection.data(), deviceVector(hostAtom));
   upload(externalField.data(), deviceVector(hostExternal));
   upload(dipoleField.data(), deviceVector(hostDipole));
   runtime.restrictMoments(atomDirection.data());

   require(runtime.energyEvaluationCount() == 0,
           "a freshly initialized runtime already recorded an energy evaluation");
   require(runtime.lastEnergy().totalJ == 0.0,
           "a freshly initialized runtime already has a cached energy");

   // Negative control: repeated field-only calls must never touch the energy
   // reduction/finalization/D2H readback path.
   for(int repeat = 0; repeat < 3; ++repeat) {
      runtime.evaluateField(atomDirection.data(), externalField.data(),
                            dipoleField.data(), atomFieldOnly.data(),
                            coarseFieldOnly.data());
      require(runtime.energyEvaluationCount() == 0,
              "evaluateField() incremented the CGP-01 energy negative control");
      require(runtime.lastEnergy().totalJ == 0.0,
              "evaluateField() updated lastEnergy_ with an invented/stale value");
   }

   const auto energy = runtime.evaluateHybrid(
      atomDirection.data(), externalField.data(), dipoleField.data(),
      atomWithEnergy.data(), coarseWithEnergy.data());
   require(runtime.energyEvaluationCount() == 1,
           "evaluateHybrid() did not increment the CGP-01 energy negative control");
   require(energy.totalJ != 0.0,
           "evaluateHybrid() produced a trivial energy on a non-uniform fixture");
   require(runtime.lastEnergy().totalJ == energy.totalJ,
           "evaluateHybrid() did not update lastEnergy_ with its own result");

   // Field parity: the field-only pass must reproduce exactly what
   // evaluateHybrid() writes for the identical inputs -- one field
   // implementation, only the energy path differs.
   const auto fieldOnlyAtom = download(atomFieldOnly.data(), atomVectors);
   const auto energyAtom = download(atomWithEnergy.data(), atomVectors);
   for(std::size_t index = 0; index < atomVectors; ++index)
      require(fieldOnlyAtom[index] == energyAtom[index],
              "evaluateField()'s atomistic field differs from evaluateHybrid()'s");
   const auto fieldOnlyCoarse = download(coarseFieldOnly.data(), coarseVectors);
   const auto energyCoarse = download(coarseWithEnergy.data(), coarseVectors);
   for(std::size_t index = 0; index < coarseVectors; ++index)
      require(fieldOnlyCoarse[index] == energyCoarse[index],
              "evaluateField()'s coarse field differs from evaluateHybrid()'s");

   // A field-only call issued after an energy call must still leave the
   // negative control and the cached energy untouched.
   runtime.evaluateField(atomDirection.data(), externalField.data(),
                         dipoleField.data(), atomFieldOnly.data(),
                         coarseFieldOnly.data());
   require(runtime.energyEvaluationCount() == 1,
           "evaluateField() incremented the negative control after a prior energy call");
   require(runtime.lastEnergy().totalJ == energy.totalJ,
           "evaluateField() overwrote lastEnergy_ after a prior energy call");

   dipoleField.Free();
   externalField.Free();
   coarseWithEnergy.Free();
   coarseFieldOnly.Free();
   atomWithEnergy.Free();
   atomFieldOnly.Free();
   atomDirection.Free();
   runtime.release();
   require(TensorMemoryTracker::current_device() == baseline,
           "field-only test cleanup did not restore tracked device memory");
}

// RCG-03 (F-14): the GPU polarization gate must mirror the Fortran
// evaluate_polarization_gate contract -- a block with a well-aligned
// channel resultant stays eligible to coarsen, while a block whose
// resultant is exactly cancelled (never normalized) is unconditionally
// unsafe. KernelFixture atoms are one-based-pair-per-block: atoms {3,4}
// (zero-based indices 2,3) are block 2 (one-based), atoms {5,6}
// (zero-based indices 4,5) are block 3 (one-based).
void testPolarizationGate() {
   KernelFixture fixture;
   const auto topology = fixture.topology();
   const auto input = fixture.runtime();
   GpuAdaptiveRuntime runtime;
   runtime.initialize(topology, input, KernelFixture::atoms, KernelFixture::ensembles);
   require(runtime.kernelsReady(), "polarization gate fixture kernels were not published ready");

   const std::size_t atomVectors = 3 * KernelFixture::atoms;
   GpuTensor<real, 1> atomDirection;
   atomDirection.Allocate(static_cast<index_t>(atomVectors));
   std::vector<double> hostAtom(atomVectors, 0.0);
   hostAtom[2 + 3 * 2] = 1.0;   // block 2, atom 3 (zero-based 2): +z
   hostAtom[2 + 3 * 3] = 1.0;   // block 2, atom 4 (zero-based 3): +z, aligned
   hostAtom[2 + 3 * 4] = 1.0;   // block 3, atom 5 (zero-based 4): +z
   hostAtom[2 + 3 * 5] = -1.0;  // block 3, atom 6 (zero-based 5): -z, exactly cancelled
   upload(atomDirection.data(), deviceVector(hostAtom));

   runtime.restrictMoments(atomDirection.data());
   runtime.evaluatePolarizationGate(real(0.9));

   std::vector<unsigned char> unsafe(KernelFixture::blocks, 0);
   ASSERT_GPU(GPU_MEMCPY(unsafe.data(), runtime.polarizationUnsafeBlockMask(),
                         unsafe.size(), GPU_MEMCPY_DEVICE_TO_HOST));
   require(unsafe[1] == 0, "aligned block was flagged unsafe to coarsen");
   require(unsafe[2] == 1, "exactly-cancelled block was not flagged unsafe");

   // RCG-03 diagnostic: the same call also reports the ratio that produced
   // each verdict above, not just the pass/fail bit.
   std::vector<real> ratio(KernelFixture::blocks, real(-1));
   ASSERT_GPU(GPU_MEMCPY(ratio.data(), runtime.polarizationRatioBlock(),
                         ratio.size() * sizeof(real), GPU_MEMCPY_DEVICE_TO_HOST));
   require(std::abs(static_cast<double>(ratio[1]) - 1.0) < 1.0e-12,
           "aligned block ratio was not measured as 1.0");
   require(static_cast<double>(ratio[2]) == 0.0,
           "exactly-cancelled block ratio was not reported at the undefined-direction floor");

   atomDirection.Free();
   runtime.release();
}

// RCG-05D descriptor layout check: GpuAdaptiveSelectorPolicy::
// bufferDilationBlocks changed from a single scalar to a 3-element per-axis
// array. This is a host-only, GPU-independent proof that the other fields
// (refineThreshold, coarsenThreshold, minimumDwellUpdates) were not shifted
// or aliased by that change, and that all three new components are
// independently addressable, by round-tripping distinct sentinel values
// through every field (including a copy, since call sites copy-construct a
// dilated variant from a base policy) and confirming none of them observe
// any other field's value.
void testSelectorPolicyDescriptorLayout() {
   GpuAdaptiveSelectorPolicy policy;
   policy.refineThreshold = real(0.123456);
   policy.coarsenThreshold = real(0.654321);
   policy.minimumDwellUpdates = 777u;
   policy.bufferDilationBlocks[0] = 11u;
   policy.bufferDilationBlocks[1] = 22u;
   policy.bufferDilationBlocks[2] = 33u;

   const auto copied = policy;
   // RCG-05-FU4: compare against the pre-copy `policy` value, not against the
   // hardcoded double literal.  `real` is `float` in a SINGLE-precision build,
   // and float(0.123456) does not round-trip to the double 0.123456, so the
   // old exact-equality-against-a-literal form failed every fp32 build and
   // aborted this binary before any later test ran -- masking, among other
   // things, RCG-08's fp32 evidence.  Comparing the two `real` values is the
   // check this test actually wants (did the field survive a copy unshifted)
   // and is precision-independent.  Confirmed by RCG-05G not to be a real
   // descriptor-layout defect.
   require(copied.refineThreshold == policy.refineThreshold,
           "descriptor layout: refineThreshold shifted across a copy");
   require(copied.coarsenThreshold == policy.coarsenThreshold,
           "descriptor layout: coarsenThreshold shifted across a copy");
   // The sentinels must still be distinguishable from each other after the
   // round trip, or the comparison above would pass on an aliased field.
   require(copied.refineThreshold != copied.coarsenThreshold,
           "descriptor layout: refineThreshold and coarsenThreshold alias");
   require(copied.minimumDwellUpdates == 777u,
           "descriptor layout: minimumDwellUpdates shifted across a copy");
   require(copied.bufferDilationBlocks[0] == 11u &&
           copied.bufferDilationBlocks[1] == 22u &&
           copied.bufferDilationBlocks[2] == 33u,
           "descriptor layout: bufferDilationBlocks components shifted, "
           "aliased each other, or aliased another field across a copy");

   // A default-constructed policy must not observe any of the sentinels
   // above (rules out a static/shared-storage aliasing bug).
   const GpuAdaptiveSelectorPolicy fresh;
   require(fresh.bufferDilationBlocks[0] == 0 &&
           fresh.bufferDilationBlocks[1] == 0 &&
           fresh.bufferDilationBlocks[2] == 0,
           "descriptor layout: default bufferDilationBlocks is not zero-initialized");
   // RCG-05-FU4: same reason as above -- compare `real` against `real`.  This
   // direction happened to pass at fp32 (float(0.123456) != double 0.123456 is
   // true), but for the wrong reason, so it is corrected too rather than left
   // as a check that only works by accident.
   require(fresh.refineThreshold != policy.refineThreshold,
           "descriptor layout: default-constructed policy observed a sentinel from another instance");
}

// RCG-08 (F-12): the compaction scan was a global Hillis--Steele sweep --
// ceil(log2(N)) launches over all 3N flags -- and is now a hierarchical
// tile scan.  Nothing in the suite above reaches past a single 256-element
// tile: the fixtures have 8 atoms and 4 blocks, so they exercise neither the
// tile-total level nor the offset propagation between levels, and would keep
// passing if either were wrong.
//
// This fixture is sized so scanItems = 76800 forces two tile-total levels
// (76800 -> 300 -> 2), exercising the complete down-propagation chain, and
// asserts two independent things:
//
//   1. Correctness -- the three compact work lists match a host reference
//      computed directly from the block states.  A wrong tile offset would
//      scatter entries to wrong positions and change the lists.
//   2. Complexity (the actual finding) -- the compaction launch count stays
//      near log_256(N) rather than log2(N).  This is the negative control for
//      F-12: the pre-fix sweep needed 19 launches at this size (17 doubling
//      steps plus initialize and scatter) and fails the bound asserted here,
//      while the hierarchical scan needs 7.  The old scan was correct but
//      slow, so a correctness-only test could not have distinguished them.
void testHierarchicalCompaction() {
   constexpr std::size_t blocks = 300;
   constexpr std::size_t basis = 256;
   constexpr std::size_t atoms = basis * blocks;

   int repetitionShape[3] = {static_cast<int>(blocks), 1, 1};
   int blockShape[3] = {1, 1, 1};
   int blockGrid[3] = {static_cast<int>(blocks), 1, 1};
   double cellVectors[9] = {1, 0, 0, 0, 1, 0, 0, 0, 1};
   double blockVectors[9] = {1, 0, 0, 0, 1, 0, 0, 0, 1};

   std::vector<int> atomToBlock(atoms), atomToBasis(atoms),
      atomToDynamic(atoms, 1), atomToFft(atoms), atomToFftGrid(atoms),
      blockAtoms(atoms), blockCount(blocks, static_cast<int>(basis)),
      blockOffset(blocks + 1), blockCoordinate(3 * blocks),
      basisPopulation(basis * blocks, 1), fftPopulation(basis * blocks, 1),
      dynamicPopulation(blocks, static_cast<int>(basis)),
      basisToDynamic(basis, 1), basisToFft(basis), state(blocks);
   std::vector<double> center(3 * blocks), volume(blocks, 1.0);
   for(std::size_t member = 0; member < basis; ++member)
      basisToFft[member] = static_cast<int>(member + 1);
   for(std::size_t block = 0; block < blocks; ++block) {
      blockOffset[block] = static_cast<int>(block * basis);
      blockCoordinate[3 * block] = static_cast<int>(block);
      center[3 * block] = static_cast<double>(block) + 0.5;
      for(std::size_t member = 0; member < basis; ++member) {
         const std::size_t atom = member + basis * block;
         atomToBlock[atom] = static_cast<int>(block + 1);
         atomToBasis[atom] = static_cast<int>(member + 1);
         atomToFft[atom] = static_cast<int>(member + 1);
         atomToFftGrid[atom] = static_cast<int>(atom + 1);
         blockAtoms[atom] = static_cast<int>(atom + 1);
      }
      // Deliberately irregular, and coprime-strided so runs of fine/coarse
      // blocks straddle tile boundaries instead of aligning with them.
      state[block] = (block * 7) % 11 < 4 ? 2 : ((block % 13) == 0 ? 1 : 0);
   }
   blockOffset[blocks] = static_cast<int>(atoms);

   GpuAdaptiveTopologyInput t;
   t.geometryMode = 1;
   t.atoms = atoms;
   t.blocks = blocks;
   t.basis = basis;
   t.fftChannelsPerBlock = basis;
   t.fftGridChannels = basis * blocks;
   t.dynamicChannels = 1;
   t.ensembles = 1;
   t.repetitionShape = repetitionShape;
   t.blockShape = blockShape;
   t.blockGrid = blockGrid;
   t.cellVectors = cellVectors;
   t.blockVectors = blockVectors;
   t.atomToBlock = atomToBlock.data();
   t.atomToBasis = atomToBasis.data();
   t.atomToDynamicChannel = atomToDynamic.data();
   t.atomToFftChannel = atomToFft.data();
   t.atomToFftGridIndex = atomToFftGrid.data();
   t.basisToDynamicChannel = basisToDynamic.data();
   t.basisToFftChannel = basisToFft.data();
   t.blockAtomCount = blockCount.data();
   t.blockAtomOffset = blockOffset.data();
   t.blockAtoms = blockAtoms.data();
   t.blockGridCoordinate = blockCoordinate.data();
   t.blockBasisPopulation = basisPopulation.data();
   t.blockFftChannelPopulation = fftPopulation.data();
   t.blockDynamicChannelPopulation = dynamicPopulation.data();
   t.blockCenter = center.data();
   t.blockVolume = volume.data();

   std::vector<double> moment(3 * blocks, 0.0), direction(3 * blocks, 0.0),
      field(3 * blocks, 0.0), momentSum(blocks, 1.0), scores(blocks, 0.0);
   GpuAdaptiveRuntimeInput r;
   r.blockState = state.data();
   r.selectorCriteria = 1;
   r.selectorScores = scores.data();
   r.coarseMoment = moment.data();
   r.coarseDirection = direction.data();
   r.coarseField = field.data();
   r.channelMomentSum = momentSum.data();

   std::string diagnostic;
   require(GpuAdaptiveRuntime::validate(t, r, atoms, 1, diagnostic),
           "hierarchical-compaction fixture was rejected: " + diagnostic);

   GpuAdaptiveRuntime runtime;
   runtime.initialize(t, r, atoms, 1);
   runtime.resetPhaseMetrics();
   runtime.updateBlockState(state.data(), blocks);

   std::vector<int> expectedAtoms, expectedBlocks, expectedInterface;
   for(std::size_t atom = 0; atom < atoms; ++atom) {
      const int blockState = state[atomToBlock[atom] - 1];
      if(blockState != 0) expectedAtoms.push_back(static_cast<int>(atom + 1));
      if(blockState == 1) expectedInterface.push_back(static_cast<int>(atom + 1));
   }
   for(std::size_t block = 0; block < blocks; ++block)
      if(state[block] == 0) expectedBlocks.push_back(static_cast<int>(block + 1));

   const auto snapshot = runtime.downloadWorkSnapshot();
   expectEqual(snapshot.activeAtoms, expectedAtoms,
               "hierarchical scan produced the wrong active atom list");
   expectEqual(snapshot.activeBlocks, expectedBlocks,
               "hierarchical scan produced the wrong coarse block list");
   expectEqual(snapshot.interfaceAtoms, expectedInterface,
               "hierarchical scan produced the wrong interface atom list");
   // RCG-09C: this fixture carries no operator inventory, so it has no bonds
   // and therefore no ghost shell.  The scan still carries the fourth
   // component; it must come back empty rather than uninitialized.
   require(snapshot.ghostAtoms.empty() && snapshot.activeBonds.empty(),
           "bondless fixture produced a nonempty ghost shell or live-bond list");

   // 76800 items -> tile levels {300, 2}: initialize, three tile scans, two
   // offset additions, scatter = 7.  RCG-09C widened the scan from three
   // components to four (the ghost shell) and added a separate single-
   // component bond scan, but neither changes the launch count here: the tile
   // scans are one launch per level regardless of component count, and this
   // fixture has no bonds.  Bound left at 8 to tolerate one extra level if
   // adaptiveThreads is ever retuned, while still failing hard on the
   // pre-fix sweep's 19.
   const auto metrics = runtime.phaseMetrics();
   require(metrics.compactionLaunches > 0,
           "compaction launch counter did not record the rebuild");
   require(metrics.compactionLaunches <= 8,
           "compaction still scales its launch count like log2(N) (F-12): got " +
              std::to_string(metrics.compactionLaunches) + " launches for " +
              std::to_string(atoms) + " scan items");
   runtime.release();
}

// CG-14 CONTINUUM-OPERATOR-HOST-REFERENCE.
//
// Why this fixture exists.  CG-14 rewrote evaluateAdaptiveCoarseTensor, and
// four deliberate defects were injected into the rewrite to find out which
// oracle would catch them.  Only two were caught, and by the whole-run
// fixtures rather than by any unit test:
//
//   dropped stencil neighbour along direction 0   caught (production e2e)
//   inverted basis-cross chirality                caught (production e2e)
//   dropped stencil neighbour along direction 2   NOT caught by anything
//   transposed inverseBlockTranspose index        NOT caught by anything
//   exchange derivative index p/q swapped         NOT caught by anything
//
// The three misses are structural, not accidental.  Every geometry that runs
// end to end today is orthogonal (RCG-05-FU3: neighbourmap.f90 rejects the
// skew fixture at setup), so inverseBlockTranspose is diagonal and transposing
// it is the identity; the tracked materials are cubic, so exchangeStiffness is
// diagonal and Sum_q C_pq grad_q collapses to a multiple of grad_p; and every
// tracked texture varies along one axis, so the transverse stencil branches
// carry exact zero.
//
// This fixture removes all three degeneracies at the unit level: a 3x3x3
// periodic block grid so each of the three forward neighbours is a distinct
// block, a dense non-symmetric inverseBlockTranspose, a dense *symmetric*
// exchangeStiffness (validate() requires symmetry, and symmetric-but-not-
// diagonal is what discriminates the p/q contract), a dense spiralization, and
// a direction field that varies in all three block directions.  It compares
// the device against a host reference written straight from the operator's
// definition -- per-term gradients, per-term scatter, the pre-CG-14
// formulation -- rather than against a restatement of the device code.
struct ContinuumFixture {
   static constexpr int side = 3;
   static constexpr std::size_t blocks = side * side * side;
   static constexpr std::size_t atoms = blocks;
   static constexpr std::size_t basis = 1;
   static constexpr std::size_t channels = 1;
   static constexpr std::size_t ensembles = 1;
   static constexpr std::size_t bonds = atoms;
   static constexpr double blockVolumeValue = 1.75;
   static constexpr double moment = 1.5;

   int repetitionShape[3] = {side, side, side};
   int blockShape[3] = {1, 1, 1};
   int blockGrid[3] = {side, side, side};
   double cellVectors[9] = {1, 0, 0, 0, 1, 0, 0, 0, 1};
   double blockVectors[9] = {1, 0, 0, 0, 1, 0, 0, 0, 1};

   // Dense and deliberately non-symmetric: B[physical + 3*direction].  A
   // transposed index therefore selects a different number.
   double inverseBlockTranspose[9] = {
      1.30, 0.40, -0.20,
      0.25, 0.90, 0.35,
      -0.15, 0.50, 1.10
   };
   // Dense and symmetric (validate() enforces symmetry).  Off-diagonal entries
   // are what make Sum_q C_pq grad_q differ from a multiple of grad_p.
   double exchange[9] = {
      0.90, 0.30, -0.20,
      0.30, 0.70, 0.15,
      -0.20, 0.15, 1.10
   };
   // Dense: D[k + 3*physical].  Every (k, p) spiralization term is live.
   double spiralization[9] = {
      0.21, -0.13, 0.34,
      0.17, 0.29, -0.11,
      -0.23, 0.19, 0.26
   };

   int atomToBlock[atoms] = {};
   int atomToBasis[atoms] = {};
   int atomToDynamic[atoms] = {};
   int atomToFft[atoms] = {};
   int atomToFftGrid[atoms] = {};
   int basisToDynamic[basis] = {1};
   int basisToFft[basis] = {1};
   int blockCount[blocks] = {};
   int blockOffset[blocks + 1] = {};
   int blockAtoms[atoms] = {};
   int blockCoordinate[3 * blocks] = {};
   int basisPopulation[basis * blocks] = {};
   int fftPopulation[basis * blocks] = {};
   int dynamicPopulation[channels * blocks] = {};
   double center[3 * blocks] = {};
   double volume[blocks] = {};
   int state[blocks] = {};
   double scores[blocks] = {};
   double atomMoment[atoms] = {};
   int atomAxisCount[atoms] = {};
   double atomAxis[6 * atoms] = {};
   double atomK1[2 * atoms] = {};
   double atomK2[2 * atoms] = {};
   int projectionBlock[8 * atoms] = {};
   double projectionWeight[8 * atoms] = {};
   int bondAtom[2 * bonds] = {};
   double bondMatrix[9 * bonds] = {};
   int selectorEdge[2 * bonds] = {};
   // No coarse anisotropy: this fixture is about the tensor operator alone.
   int axisCount[blocks] = {};
   double axis[6 * blocks] = {};
   double k1[2 * blocks] = {};
   double k2[2 * blocks] = {};

   std::vector<double> coarseMoment =
      std::vector<double>(3 * blocks, 0.0);
   std::vector<double> coarseDirection =
      std::vector<double>(3 * blocks, 0.0);
   std::vector<double> coarseField =
      std::vector<double>(3 * blocks, 0.0);
   std::vector<double> momentSum =
      std::vector<double>(blocks, moment);

   ContinuumFixture() {
      for(std::size_t block = 0; block < blocks; ++block) {
         const int ix = static_cast<int>(block % side);
         const int iy = static_cast<int>((block / side) % side);
         const int iz = static_cast<int>(block / (side * side));
         blockCoordinate[0 + 3 * block] = ix;
         blockCoordinate[1 + 3 * block] = iy;
         blockCoordinate[2 + 3 * block] = iz;
         center[0 + 3 * block] = ix + 0.5;
         center[1 + 3 * block] = iy + 0.5;
         center[2 + 3 * block] = iz + 0.5;
         volume[block] = blockVolumeValue;
         blockCount[block] = 1;
         blockOffset[block] = static_cast<int>(block);
         blockAtoms[block] = static_cast<int>(block) + 1;
         basisPopulation[block] = 1;
         fftPopulation[block] = 1;
         dynamicPopulation[block] = 1;
         state[block] = 0;   // every block coarse
         atomToBlock[block] = static_cast<int>(block) + 1;
         atomToBasis[block] = 1;
         atomToDynamic[block] = 1;
         atomToFft[block] = 1;
         atomToFftGrid[block] = static_cast<int>(block) + 1;
         atomMoment[block] = moment;
         for(int corner = 0; corner < 8; ++corner)
            projectionBlock[corner + 8 * block] = static_cast<int>(block) + 1;
         projectionWeight[8 * block] = 1.0;
         // A texture that varies along all three block directions, so no
         // stencil branch and no gradient component is identically zero.
         const double theta = 0.35 + 0.41 * ix + 0.23 * iy + 0.17 * iz;
         const double phi = 0.19 + 0.37 * ix - 0.29 * iy + 0.53 * iz;
         coarseDirection[0 + 3 * block] = std::sin(theta) * std::cos(phi);
         coarseDirection[1 + 3 * block] = std::sin(theta) * std::sin(phi);
         coarseDirection[2 + 3 * block] = std::cos(theta);
         for(int xyz = 0; xyz < 3; ++xyz)
            coarseMoment[xyz + 3 * block] =
               moment * coarseDirection[xyz + 3 * block];
      }
      blockOffset[blocks] = static_cast<int>(atoms);
      // A periodic bond ring.  Every bond is coarse-coarse, so the continuum
      // operator owns all of them and no atomistic work is launched.
      for(std::size_t bond = 0; bond < bonds; ++bond) {
         bondAtom[bond] = static_cast<int>(bond) + 1;
         bondAtom[bonds + bond] = static_cast<int>((bond + 1) % bonds) + 1;
         selectorEdge[bond] = bondAtom[bond];
         selectorEdge[bonds + bond] = bondAtom[bonds + bond];
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

   GpuAdaptiveRuntimeInput runtime() {
      GpuAdaptiveRuntimeInput r;
      r.blockState = state;
      r.selectorCriteria = 1;
      r.selectorScores = scores;
      r.coarseMoment = coarseMoment.data();
      r.coarseDirection = coarseDirection.data();
      r.coarseField = coarseField.data();
      r.channelMomentSum = momentSum.data();
      r.kernels.atomMoment = atomMoment;
      r.kernels.atomAnisotropyAxisCount = atomAxisCount;
      r.kernels.atomAnisotropyAxis = atomAxis;
      r.kernels.atomAnisotropyK1 = atomK1;
      r.kernels.atomAnisotropyK2 = atomK2;
      r.kernels.projectionBlock = projectionBlock;
      r.kernels.projectionWeight = projectionWeight;
      r.kernels.bonds = bonds;
      r.kernels.bondAtom = bondAtom;
      r.kernels.bondMatrix = bondMatrix;
      r.kernels.selectorEdges = bonds;
      r.kernels.selectorEdge = selectorEdge;
      r.kernels.inverseBlockTranspose = inverseBlockTranspose;
      r.kernels.exchangeStiffness = exchange;
      r.kernels.spiralization = spiralization;
      r.kernels.anisotropyAxisCount = axisCount;
      r.kernels.anisotropyAxis = axis;
      r.kernels.anisotropyK1 = k1;
      r.kernels.anisotropyK2 = k2;
      r.kernels.magneticMomentSi = 1.0;
      r.kernels.gammaPerTs = 2.0;
      r.kernels.damping = 0.1;
      return r;
   }
};

// The host reference.  Written from the operator's definition in the same
// shape the pre-CG-14 device kernel had: neighbour indices, gradients and the
// stencil are all evaluated per tensor term, without any of the hoisting or
// factoring CG-14 introduced.  Its agreement with the device is therefore
// evidence about the rewrite, not a restatement of it.
struct ContinuumReference {
   std::vector<double> derivative;   // dE/dm, before the -1/(mu*M) scaling
   double exchangeEnergy = 0.0;
   double spiralEnergy = 0.0;
};

ContinuumReference continuumHostReference(const ContinuumFixture& fixture) {
   constexpr int side = ContinuumFixture::side;
   constexpr std::size_t blocks = ContinuumFixture::blocks;
   const double* B = fixture.inverseBlockTranspose;
   const double* C = fixture.exchange;
   const double* D = fixture.spiralization;
   const std::vector<double>& m = fixture.coarseDirection;

   ContinuumReference result;
   result.derivative.assign(3 * blocks, 0.0);

   auto plusBlock = [&](std::size_t block, int direction) {
      int c[3] = {fixture.blockCoordinate[0 + 3 * block],
                  fixture.blockCoordinate[1 + 3 * block],
                  fixture.blockCoordinate[2 + 3 * block]};
      c[direction] = (c[direction] + 1) % side;
      return static_cast<std::size_t>(c[0] + side * (c[1] + side * c[2]));
   };
   auto gradient = [&](std::size_t block, int physical, double g[3]) {
      g[0] = g[1] = g[2] = 0.0;
      for(int direction = 0; direction < 3; ++direction) {
         const double coefficient = B[physical + 3 * direction];
         const std::size_t plus = plusBlock(block, direction);
         for(int xyz = 0; xyz < 3; ++xyz)
            g[xyz] += coefficient * (m[xyz + 3 * plus] - m[xyz + 3 * block]);
      }
   };
   auto scatter = [&](std::size_t block, int physical, const double value[3],
                      double scale) {
      for(int direction = 0; direction < 3; ++direction) {
         const double coefficient = scale * B[physical + 3 * direction];
         const std::size_t plus = plusBlock(block, direction);
         for(int xyz = 0; xyz < 3; ++xyz) {
            result.derivative[xyz + 3 * plus] += coefficient * value[xyz];
            result.derivative[xyz + 3 * block] -= coefficient * value[xyz];
         }
      }
   };
   auto cross = [](const double a[3], const double b[3], double out[3]) {
      out[0] = a[1] * b[2] - a[2] * b[1];
      out[1] = a[2] * b[0] - a[0] * b[2];
      out[2] = a[0] * b[1] - a[1] * b[0];
   };

   for(std::size_t block = 0; block < blocks; ++block) {
      const double volume = fixture.volume[block];
      for(int p = 0; p < 3; ++p)
         for(int q = 0; q < 3; ++q) {
            double gp[3], gq[3], derivative[3];
            gradient(block, p, gp);
            gradient(block, q, gq);
            const double stiffness = C[p + 3 * q];
            result.exchangeEnergy += volume * stiffness *
               (gp[0] * gq[0] + gp[1] * gq[1] + gp[2] * gq[2]);
            for(int xyz = 0; xyz < 3; ++xyz)
               derivative[xyz] = volume * stiffness * gq[xyz];
            scatter(block, p, derivative, 2.0);
         }
      double direction[3] = {m[0 + 3 * block], m[1 + 3 * block],
                             m[2 + 3 * block]};
      for(int k = 0; k < 3; ++k) {
         double basis[3] = {0.0, 0.0, 0.0};
         basis[k] = 1.0;
         for(int p = 0; p < 3; ++p) {
            double g[3], crossMG[3], crossBG[3], crossBD[3], derivative[3];
            gradient(block, p, g);
            cross(direction, g, crossMG);
            cross(basis, g, crossBG);
            cross(basis, direction, crossBD);
            const double spiral = D[k + 3 * p];
            result.spiralEnergy += volume * spiral * crossMG[k];
            for(int xyz = 0; xyz < 3; ++xyz) {
               result.derivative[xyz + 3 * block] -=
                  volume * spiral * crossBG[xyz];
               derivative[xyz] = volume * spiral * crossBD[xyz];
            }
            scatter(block, p, derivative, 1.0);
         }
      }
   }
   return result;
}

void testContinuumOperatorAgainstHostReference() {
   ContinuumFixture fixture;
   const auto topology = fixture.topology();
   auto input = fixture.runtime();
   std::string diagnostic;
   require(GpuAdaptiveRuntime::validate(topology, input, ContinuumFixture::atoms,
                                        ContinuumFixture::ensembles, diagnostic),
           "CG-14 continuum reference fixture was rejected: " + diagnostic);

   GpuAdaptiveRuntime runtime;
   runtime.initialize(topology, input, ContinuumFixture::atoms,
                      ContinuumFixture::ensembles);
   require(runtime.activeBlockCount() == ContinuumFixture::blocks,
           "CG-14 reference fixture did not put every block in the coarse list");
   require(runtime.liveBondCount() == 0,
           "CG-14 reference fixture launched atomistic work, so its coarse "
           "field would not isolate the continuum operator");

   const std::size_t atomVectors = 3 * ContinuumFixture::atoms;
   const std::size_t coarseVectors = 3 * ContinuumFixture::blocks;
   GpuTensor<real, 1> atomDirection, atomField, coarseField;
   atomDirection.Allocate(static_cast<index_t>(atomVectors));
   atomField.Allocate(static_cast<index_t>(atomVectors));
   coarseField.Allocate(static_cast<index_t>(coarseVectors));
   // Every block is coarse, so restriction keeps the fixture's own coarse
   // directions and the atom directions never enter the coarse equation.
   upload(atomDirection.data(),
          deviceVector(std::vector<double>(atomVectors, 0.0)));
   runtime.restrictMoments(atomDirection.data());
   const auto energy = runtime.evaluateHybrid(
      atomDirection.data(), nullptr, nullptr, atomField.data(),
      coarseField.data());
   const auto deviceField = download(coarseField.data(), coarseVectors);

   const auto reference = continuumHostReference(fixture);

   // A gate on the fixture itself: if the reference operator produced a
   // near-zero field or near-zero energies, agreement would be vacuous.
   double referenceScale = 0.0;
   for(const double value : reference.derivative)
      referenceScale = std::max(referenceScale, std::abs(value));
   require(referenceScale > 1.0e-3,
           "CG-14 reference fixture produced a vanishing continuum derivative");
   require(std::abs(reference.exchangeEnergy) > 1.0e-3 &&
           std::abs(reference.spiralEnergy) > 1.0e-3,
           "CG-14 reference fixture produced vanishing coarse energies");

   const double tolerance = sizeof(real) == sizeof(double) ? 1.0e-11 : 5.0e-4;
   require(std::abs(energy.coarseExchangeJ - reference.exchangeEnergy) <=
              tolerance * std::abs(reference.exchangeEnergy),
           "GPU coarse exchange energy disagrees with the host reference "
           "operator");
   require(std::abs(energy.coarseSpiralizationJ - reference.spiralEnergy) <=
              tolerance * std::abs(reference.spiralEnergy),
           "GPU coarse spiralization energy disagrees with the host reference "
           "operator");

   // The published coarse field is -(dE/dm)/(mu*M); there is no external
   // field, no dipole and no coarse anisotropy in this fixture, and the
   // interface scratch is zero because no atom is atomistically owned.
   const double scale = -1.0 / (1.0 * ContinuumFixture::moment);
   double worst = 0.0;
   for(std::size_t component = 0; component < coarseVectors; ++component) {
      const double expected = scale * reference.derivative[component];
      worst = std::max(worst,
         std::abs(static_cast<double>(deviceField[component]) - expected));
   }
   require(worst <= tolerance * referenceScale * std::abs(scale),
           "GPU continuum coarse field disagrees with the host reference "
           "operator: worst component deviation " + std::to_string(worst));
   runtime.release();
}

// CGP-02: synchronizeAtomicState() was previously exercised only through the
// Fortran GPU e2e fixtures, never this unit-test binary, so the new
// policy-dependent branch (Aligned launches the compact
// prolongateAdaptiveGhostsCompact; ConstrainedCone keeps the untouched
// full-range prolongateAdaptiveGhosts, since commitAdaptiveGhosts reads
// ghostDirection for every non-atomistic atom under that policy -- see
// docs/rcg09/rcg09c_live_bond_compaction.txt section 6) is exercised here
// directly. KernelFixture's projection is degenerate (every corner maps to
// the atom's own block with corner-0 weight 1, see its constructor above),
// so prolongation reduces to an exact copy of the owning block's
// coarseDirection for every non-atomistic atom -- meaning Aligned (which
// overwrites with coarseDirection directly) and ConstrainedCone (which
// commits the prolongated ghostDirection unchanged) must agree exactly.
// That equality cross-checks the new compact kernel's math against the
// unmodified full-range kernel it replaces at this call site.
void testSynchronizeAtomicStatePolicies() {
   KernelFixture fixture;
   const auto topology = fixture.topology();
   const auto input = fixture.runtime();
   GpuAdaptiveRuntime runtime;
   runtime.initialize(topology, input, KernelFixture::atoms,
                      KernelFixture::ensembles);
   require(runtime.ghostAtomCount() == 2,
           "fixture ghost shell drifted from the documented {2,7}");

   const std::size_t atomVectors = 3 * KernelFixture::atoms;
   GpuTensor<real, 1> atomDirection;
   atomDirection.Allocate(static_cast<index_t>(atomVectors));
   std::vector<double> hostAtom(atomVectors, 0.0);
   // Sentinel x-direction so an atomistic atom the commit is not supposed to
   // touch, or a non-atomistic atom the commit forgets to touch, is
   // detectable.
   for(std::size_t atom = 0; atom < KernelFixture::atoms; ++atom)
      hostAtom[0 + 3 * atom] = 1.0;
   upload(atomDirection.data(), deviceVector(hostAtom));

   const double tolerance = sizeof(real) == sizeof(double) ? 1.0e-12 : 1.0e-5;
   const auto beforeAligned = runtime.phaseMetrics();
   runtime.synchronizeAtomicState(atomDirection.data(), {});
   const auto afterAligned = runtime.phaseMetrics();
   require(afterAligned.integrationLaunches ==
              beforeAligned.integrationLaunches + 2,
           "Aligned synchronizeAtomicState with a nonempty ghost shell "
           "issued an unexpected launch count (expected exactly the compact "
           "prolongation plus the commit)");
   const auto aligned = download(atomDirection.data(), atomVectors);
   // Non-atomistic atoms 1,2,7,8 (0-indexed 0,1,6,7) must align to their
   // block's +z direction; atomistic atoms 3-6 (0-indexed 2-5) must be
   // untouched (still the +x sentinel).
   for(const std::size_t atom : {std::size_t(0), std::size_t(1),
                                  std::size_t(6), std::size_t(7)})
      require(std::abs(static_cast<double>(aligned[2 + 3 * atom]) - 1.0) <
                 tolerance &&
              std::abs(static_cast<double>(aligned[0 + 3 * atom])) <
                 tolerance,
              "Aligned reconstruction did not align a non-atomistic atom to "
              "its block direction");
   for(const std::size_t atom : {std::size_t(2), std::size_t(3),
                                  std::size_t(4), std::size_t(5)})
      require(std::abs(static_cast<double>(aligned[0 + 3 * atom]) - 1.0) <
                 tolerance,
              "Aligned reconstruction touched an atomistic atom's direction");

   upload(atomDirection.data(), deviceVector(hostAtom));
   GpuAdaptiveReconstructionPolicy cone;
   cone.scheme = GpuAdaptiveReconstruction::ConstrainedCone;
   const auto beforeCone = runtime.phaseMetrics();
   runtime.synchronizeAtomicState(atomDirection.data(), cone);
   const auto afterCone = runtime.phaseMetrics();
   require(afterCone.integrationLaunches == beforeCone.integrationLaunches + 2,
           "ConstrainedCone synchronizeAtomicState issued an unexpected "
           "launch count (expected exactly the full-range prolongation plus "
           "the commit)");
   const auto cone_result = download(atomDirection.data(), atomVectors);
   for(std::size_t atom = 0; atom < KernelFixture::atoms; ++atom)
      for(int xyz = 0; xyz < 3; ++xyz)
         require(std::abs(static_cast<double>(cone_result[xyz + 3 * atom]) -
                          static_cast<double>(aligned[xyz + 3 * atom])) <
                    tolerance,
                 "ConstrainedCone and Aligned reconstruction disagree under "
                 "a degenerate single-corner projection, where they must "
                 "match exactly");

   atomDirection.Free();
   runtime.release();
}

// CGP-02: no existing GPU adaptive-runtime test exercises ghost
// prolongation's ensemble indexing (`work / ghostAtoms` in
// prolongateAdaptiveGhostsCompact, and the analogous division in
// synchronizeAtomicState's full-range fallback) because every fixture above
// hardcodes ensembles=1, so every atom happens to occupy ensemble 0 and a
// swapped or aliased ensemble index would be invisible. This fixture is
// KernelFixture's exact eight-atom/four-block mixed chain (ghost shell
// {2,7}, documented above testLiveBondCompaction()) at a configurable
// ensemble count, with a different uniform coarse direction/moment in each
// ensemble so a cross-ensemble mixup produces a detectably wrong
// reconstructed direction rather than silently reusing another ensemble's
// answer.
struct GhostEnsembleFixture {
   static constexpr std::size_t atoms = 8;
   static constexpr std::size_t blocks = 4;
   static constexpr std::size_t basis = 2;
   static constexpr std::size_t channels = 1;
   static constexpr std::size_t bonds = 8;
   std::size_t ensembles;
   int repetitionShape[3] = {4, 1, 1};
   int blockShape[3] = {1, 1, 1};
   int blockGrid[3] = {4, 1, 1};
   double cellVectors[9] = {1, 0, 0, 0, 1, 0, 0, 0, 1};
   double blockVectors[9] = {1, 0, 0, 0, 1, 0, 0, 0, 1};
   int atomToBlock[atoms] = {1, 1, 2, 2, 3, 3, 4, 4};
   int atomToBasis[atoms] = {1, 2, 1, 2, 1, 2, 1, 2};
   int atomToDynamic[atoms] = {1, 1, 1, 1, 1, 1, 1, 1};
   int atomToFft[atoms] = {1, 2, 1, 2, 1, 2, 1, 2};
   int atomToFftGrid[atoms] = {1, 2, 3, 4, 5, 6, 7, 8};
   int basisToDynamic[basis] = {1, 1};
   int basisToFft[basis] = {1, 2};
   int blockCount[blocks] = {2, 2, 2, 2};
   int blockOffset[blocks + 1] = {0, 2, 4, 6, 8};
   int blockAtoms[atoms] = {1, 2, 3, 4, 5, 6, 7, 8};
   int blockCoordinate[3 * blocks] = {0, 0, 0, 1, 0, 0, 2, 0, 0, 3, 0, 0};
   int basisPopulation[basis * blocks] = {1, 1, 1, 1, 1, 1, 1, 1};
   int fftPopulation[basis * blocks] = {1, 1, 1, 1, 1, 1, 1, 1};
   int dynamicPopulation[channels * blocks] = {2, 2, 2, 2};
   double center[3 * blocks] = {0.5, 0, 0, 1.5, 0, 0, 2.5, 0, 0, 3.5, 0, 0};
   double volume[blocks] = {1, 1, 1, 1};
   int state[blocks] = {0, 1, 2, 0};
   double scores[blocks] = {};
   double atomMoment[atoms] = {1, 1, 1, 1, 1, 1, 1, 1};
   int atomAxisCount[atoms] = {};
   double atomAxis[6 * atoms] = {};
   double atomK1[2 * atoms] = {};
   double atomK2[2 * atoms] = {};
   int projectionBlock[8 * atoms] = {};
   double projectionWeight[8 * atoms] = {};
   int bondAtom[2 * bonds] = {
      1, 2, 3, 4, 5, 6, 7, 8,
      2, 3, 4, 5, 6, 7, 8, 1
   };
   double bondMatrix[9 * bonds] = {};
   int selectorEdge[2 * bonds] = {
      1, 2, 3, 4, 5, 6, 7, 8,
      2, 3, 4, 5, 6, 7, 8, 1
   };
   double inverseBlockTranspose[9] = {1, 0, 0, 0, 1, 0, 0, 0, 1};
   double exchange[9] = {};
   double spiralization[9] = {};
   int axisCount[blocks] = {};
   double axis[6 * blocks] = {};
   double k1[2 * blocks] = {};
   double k2[2 * blocks] = {};
   std::vector<double> coarseMoment;
   std::vector<double> coarseDirection;
   std::vector<double> coarseField;
   std::vector<double> momentSum;
   // Field-evaluation input: a distinct uniform direction per ensemble,
   // matching coarseDirection below.
   std::vector<double> atomDirectionHost;

   explicit GhostEnsembleFixture(std::size_t ensembleCount)
      : ensembles(ensembleCount),
        coarseMoment(3 * channels * blocks * ensembles, 0.0),
        coarseDirection(3 * channels * blocks * ensembles, 0.0),
        coarseField(3 * channels * blocks * ensembles, 0.0),
        momentSum(channels * blocks * ensembles, 2.0),
        atomDirectionHost(3 * atoms * ensembles, 0.0) {
      for(std::size_t atom = 0; atom < atoms; ++atom) {
         projectionBlock[8 * atom] = atomToBlock[atom];
         projectionWeight[8 * atom] = 1.0;
         for(int corner = 1; corner < 8; ++corner)
            projectionBlock[corner + 8 * atom] = atomToBlock[atom];
      }
      for(std::size_t bond = 0; bond < bonds; ++bond)
         for(int xyz = 0; xyz < 3; ++xyz)
            bondMatrix[xyz + 3 * xyz + 9 * bond] = 0.4;
      // Ensemble e gets a uniform direction along axis (e % 3) with moment
      // (2.0 + e) -- e=0 reproduces KernelFixture's own +z/2.0 exactly, e=1
      // is +x/3.0, so a two-ensemble run has two axis-disjoint references.
      for(std::size_t ensemble = 0; ensemble < ensembles; ++ensemble) {
         const int component = static_cast<int>(ensemble % 3);
         const double moment = 2.0 + static_cast<double>(ensemble);
         for(std::size_t block = 0; block < blocks; ++block) {
            coarseMoment[component + 3 * (block + blocks * ensemble)] = moment;
            coarseDirection[component + 3 * (block + blocks * ensemble)] = 1.0;
         }
         for(std::size_t atom = 0; atom < atoms; ++atom)
            atomDirectionHost[component + 3 * (atom + atoms * ensemble)] = 1.0;
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
      r.selectorCriteria = 1;
      r.selectorScores = scores;
      r.coarseMoment = coarseMoment.data();
      r.coarseDirection = coarseDirection.data();
      r.coarseField = coarseField.data();
      r.channelMomentSum = momentSum.data();
      r.kernels.atomMoment = atomMoment;
      r.kernels.atomAnisotropyAxisCount = atomAxisCount;
      r.kernels.atomAnisotropyAxis = atomAxis;
      r.kernels.atomAnisotropyK1 = atomK1;
      r.kernels.atomAnisotropyK2 = atomK2;
      r.kernels.projectionBlock = projectionBlock;
      r.kernels.projectionWeight = projectionWeight;
      r.kernels.bonds = bonds;
      r.kernels.bondAtom = bondAtom;
      r.kernels.bondMatrix = bondMatrix;
      r.kernels.selectorEdges = bonds;
      r.kernels.selectorEdge = selectorEdge;
      r.kernels.inverseBlockTranspose = inverseBlockTranspose;
      r.kernels.exchangeStiffness = exchange;
      r.kernels.spiralization = spiralization;
      r.kernels.anisotropyAxisCount = axisCount;
      r.kernels.anisotropyAxis = axis;
      r.kernels.anisotropyK1 = k1;
      r.kernels.anisotropyK2 = k2;
      r.kernels.magneticMomentSi = 1.0;
      r.kernels.gammaPerTs = 2.0;
      r.kernels.damping = 0.1;
      return r;
   }
};

void testMultiEnsembleGhostProlongation() {
   GhostEnsembleFixture combined(2);
   const auto topology = combined.topology();
   const auto input = combined.runtime();
   GpuAdaptiveRuntime runtime;
   runtime.initialize(topology, input, GhostEnsembleFixture::atoms,
                      combined.ensembles);
   const auto snapshot = runtime.downloadWorkSnapshot();
   expectEqual(snapshot.ghostAtoms, {2, 7},
               "multi-ensemble ghost shell differs from the single-ensemble "
               "reference -- the ghost list is supposed to be a physical-atom "
               "set shared across ensembles");

   const std::size_t atomVectors = 3 * GhostEnsembleFixture::atoms * combined.ensembles;
   const std::size_t coarseVectors = 3 * GhostEnsembleFixture::blocks * combined.ensembles;
   GpuTensor<real, 1> atomDirection, atomField, coarseField;
   atomDirection.Allocate(static_cast<index_t>(atomVectors));
   atomField.Allocate(static_cast<index_t>(atomVectors));
   coarseField.Allocate(static_cast<index_t>(coarseVectors));
   upload(atomDirection.data(), deviceVector(combined.atomDirectionHost));

   const auto combinedEnergy = runtime.evaluateHybrid(
      atomDirection.data(), nullptr, nullptr, atomField.data(), coarseField.data());
   const auto combinedField = download(atomField.data(), atomVectors);
   const auto combinedCoarse = download(coarseField.data(), coarseVectors);
   runtime.synchronizeAtomicState(atomDirection.data(), {});
   const auto combinedReconstructed = download(atomDirection.data(), atomVectors);
   runtime.release();

   const double tolerance = sizeof(real) == sizeof(double) ? 1.0e-11 : 5.0e-5;
   double soloEnergySum = 0.0;
   double combinedEnergySum = combinedEnergy.totalJ;
   for(std::size_t ensemble = 0; ensemble < combined.ensembles; ++ensemble) {
      GhostEnsembleFixture solo(1);
      // Overwrite the solo fixture's single-ensemble slice with this
      // ensemble's slice from the combined fixture, so the solo run is an
      // independent ensembles=1 oracle for exactly the data the combined
      // run assigned to this ensemble.
      for(std::size_t block = 0; block < GhostEnsembleFixture::blocks; ++block)
         for(int xyz = 0; xyz < 3; ++xyz) {
            solo.coarseMoment[xyz + 3 * block] =
               combined.coarseMoment[xyz + 3 * (block + GhostEnsembleFixture::blocks * ensemble)];
            solo.coarseDirection[xyz + 3 * block] =
               combined.coarseDirection[xyz + 3 * (block + GhostEnsembleFixture::blocks * ensemble)];
         }
      for(std::size_t atom = 0; atom < GhostEnsembleFixture::atoms; ++atom)
         for(int xyz = 0; xyz < 3; ++xyz)
            solo.atomDirectionHost[xyz + 3 * atom] =
               combined.atomDirectionHost[xyz + 3 * (atom + GhostEnsembleFixture::atoms * ensemble)];

      const auto soloTopology = solo.topology();
      const auto soloInput = solo.runtime();
      GpuAdaptiveRuntime soloRuntime;
      soloRuntime.initialize(soloTopology, soloInput, GhostEnsembleFixture::atoms, 1);
      GpuTensor<real, 1> soloDirection, soloField, soloCoarse;
      soloDirection.Allocate(static_cast<index_t>(3 * GhostEnsembleFixture::atoms));
      soloField.Allocate(static_cast<index_t>(3 * GhostEnsembleFixture::atoms));
      soloCoarse.Allocate(static_cast<index_t>(3 * GhostEnsembleFixture::blocks));
      upload(soloDirection.data(), deviceVector(solo.atomDirectionHost));
      const auto soloEnergy = soloRuntime.evaluateHybrid(
         soloDirection.data(), nullptr, nullptr, soloField.data(), soloCoarse.data());
      const auto soloField_ = download(soloField.data(), 3 * GhostEnsembleFixture::atoms);
      const auto soloCoarse_ = download(soloCoarse.data(), 3 * GhostEnsembleFixture::blocks);
      soloRuntime.synchronizeAtomicState(soloDirection.data(), {});
      const auto soloReconstructed = download(soloDirection.data(), 3 * GhostEnsembleFixture::atoms);
      soloEnergySum += soloEnergy.totalJ;

      for(std::size_t atom = 0; atom < GhostEnsembleFixture::atoms; ++atom)
         for(int xyz = 0; xyz < 3; ++xyz) {
            const std::size_t combinedIndex =
               xyz + 3 * (atom + GhostEnsembleFixture::atoms * ensemble);
            require(std::abs(static_cast<double>(combinedField[combinedIndex]) -
                             static_cast<double>(soloField_[xyz + 3 * atom])) < tolerance,
                    "combined multi-ensemble field does not match the solo "
                    "single-ensemble oracle -- ensemble indexing in the "
                    "compact prolongation launch is suspect");
            require(std::abs(static_cast<double>(combinedReconstructed[combinedIndex]) -
                             static_cast<double>(soloReconstructed[xyz + 3 * atom])) < tolerance,
                    "combined multi-ensemble synchronizeAtomicState does not "
                    "match the solo single-ensemble oracle");
         }
      for(std::size_t block = 0; block < GhostEnsembleFixture::blocks; ++block)
         for(int xyz = 0; xyz < 3; ++xyz) {
            const std::size_t combinedIndex =
               xyz + 3 * (block + GhostEnsembleFixture::blocks * ensemble);
            require(std::abs(static_cast<double>(combinedCoarse[combinedIndex]) -
                             static_cast<double>(soloCoarse_[xyz + 3 * block])) < tolerance,
                    "combined multi-ensemble coarse field does not match the "
                    "solo single-ensemble oracle");
         }
      soloCoarse.Free();
      soloField.Free();
      soloDirection.Free();
      soloRuntime.release();
   }
   require(std::abs(combinedEnergySum - soloEnergySum) <
              tolerance * std::max(1.0, std::abs(soloEnergySum)),
           "combined multi-ensemble total energy is not the sum of the "
           "independent per-ensemble energies");

   atomDirection.Free();
   atomField.Free();
   coarseField.Free();
}

// CGP-02 negative control: drop a required ghost atom from the compact list
// -- without changing ghostAtomCount(), i.e. only the *set* of processed
// atoms is wrong, not the launch geometry -- and show that the field the
// atomistic reaction reaches at the interface is now wrong. This proves the
// interface-restriction/prolongation pair genuinely depends on
// ghostAtomList's contents, not merely its length.
void testGhostListNegativeControl() {
   KernelFixture fixture;
   const auto topology = fixture.topology();
   const auto input = fixture.runtime();
   GpuAdaptiveRuntime runtime;
   runtime.initialize(topology, input, KernelFixture::atoms,
                      KernelFixture::ensembles);
   const auto snapshot = runtime.downloadWorkSnapshot();
   expectEqual(snapshot.ghostAtoms, {2, 7},
               "fixture ghost shell drifted from the documented {2,7}");

   const std::size_t atomVectors = 3 * KernelFixture::atoms;
   const std::size_t coarseVectors = 3 * KernelFixture::blocks;
   GpuTensor<real, 1> atomDirection, atomField, coarseField;
   atomDirection.Allocate(static_cast<index_t>(atomVectors));
   atomField.Allocate(static_cast<index_t>(atomVectors));
   coarseField.Allocate(static_cast<index_t>(coarseVectors));
   // The fixture's coarseDirection is uniformly +z (see KernelFixture's
   // constructor above); an atom direction that is also +z would make the
   // ghost-shell restriction's tangential component (field minus its
   // projection onto ghostDirection) identically zero regardless of which
   // atom is processed, which would make this negative control vacuous. +x
   // is orthogonal to that, so the restricted tangential field is the whole
   // field and genuinely depends on which atom contributed it.
   std::vector<double> hostAtom(atomVectors, 0.0);
   for(std::size_t atom = 0; atom < KernelFixture::atoms; ++atom)
      hostAtom[0 + 3 * atom] = 1.0;
   upload(atomDirection.data(), deviceVector(hostAtom));

   (void)runtime.evaluateHybrid(atomDirection.data(), nullptr, nullptr,
                                atomField.data(), coarseField.data());
   const auto goodCoarse = download(coarseField.data(), coarseVectors);

   // Overwrite ghost atom 7's compact slot (1-based index 7, 0-indexed
   // block 3) with a duplicate of ghost atom 2's slot; ghostAtomCount() is
   // unchanged so the launch still issues exactly the same number of
   // threads, but atom 7 is no longer among the atoms those threads visit.
   const std::vector<int> corrupted = {2, 2};
   ASSERT_GPU(GPU_MEMCPY(runtime.deviceRuntime().ghostAtomList, corrupted.data(),
                         corrupted.size() * sizeof(int), GPU_MEMCPY_HOST_TO_DEVICE));

   (void)runtime.evaluateHybrid(atomDirection.data(), nullptr, nullptr,
                                atomField.data(), coarseField.data());
   const auto corruptedCoarse = download(coarseField.data(), coarseVectors);

   const double tolerance = sizeof(real) == sizeof(double) ? 1.0e-12 : 1.0e-5;
   // Block 3 (0-indexed, owned by the now-dropped ghost atom 7) must lose
   // its interface restriction and revert toward the unrestricted value;
   // block 0 (owned by ghost atom 2, now double-scattered) must change too.
   bool block3Differs = false;
   bool block0Differs = false;
   for(int xyz = 0; xyz < 3; ++xyz) {
      if(std::abs(static_cast<double>(goodCoarse[xyz + 3 * 3]) -
                  static_cast<double>(corruptedCoarse[xyz + 3 * 3])) > tolerance)
         block3Differs = true;
      if(std::abs(static_cast<double>(goodCoarse[xyz + 3 * 0]) -
                  static_cast<double>(corruptedCoarse[xyz + 3 * 0])) > tolerance)
         block0Differs = true;
   }
   require(block3Differs,
           "dropping ghost atom 7 from the compact list did not change its "
           "block's coarse field -- the negative control is not "
           "discriminating");
   require(block0Differs,
           "duplicating ghost atom 2 into atom 7's compact slot did not "
           "change atom 2's block's coarse field -- the negative control is "
           "not discriminating");

   coarseField.Free();
   atomField.Free();
   atomDirection.Free();
   runtime.release();
}

// CGP-03: transitionsPublishedLastCall() must be an exact count of blocks
// whose blockState actually changed -- zero when every proposal is rejected,
// and equal to an independently-computed before/after blockState diff when
// some are accepted. materializeCoarseAtoms() must reproduce
// synchronizeAtomicState() bitwise and credit only the reason it was called
// with.
void testTransitionCounterAndMaterialization() {
   KernelFixture fixture;
   const auto topology = fixture.topology();
   const auto input = fixture.runtime();
   GpuAdaptiveRuntime runtime;
   runtime.initialize(topology, input, KernelFixture::atoms,
                      KernelFixture::ensembles);

   const std::size_t atomVectors = 3 * KernelFixture::atoms;
   GpuTensor<real, 1> atomDirection;
   atomDirection.Allocate(static_cast<index_t>(atomVectors));
   std::vector<double> hostAtom(atomVectors, 0.0);
   for(std::size_t atom = 0; atom < KernelFixture::atoms; ++atom)
      hostAtom[0 + 3 * atom] = 1.0;
   upload(atomDirection.data(), deviceVector(hostAtom));

   auto downloadBlockState = [&]() {
      std::vector<int> state(KernelFixture::blocks);
      ASSERT_GPU(GPU_MEMCPY(state.data(), runtime.deviceRuntime().blockState,
                            state.size() * sizeof(int), GPU_MEMCPY_DEVICE_TO_HOST));
      return state;
   };
   const auto blockStateBefore = downloadBlockState();
   require(blockStateBefore == std::vector<int>({0, 1, 2, 0}),
           "fixture initial block state drifted from the documented {0,1,2,0}");

   require(runtime.transitionsPublishedLastCall() == 0,
           "transitionsPublishedLastCall() must start at zero, before any publish");

   // Rejecting every proposal must publish exactly zero transitions and
   // leave blockState untouched.
   runtime.evaluateSelectorScores(atomDirection.data());
   GpuAdaptiveSelectorPolicy policy;
   policy.refineThreshold = real(0.5);
   policy.coarsenThreshold = real(0.1);
   runtime.proposeSelectorState(policy);
   GpuTensor<unsigned char, 1> none;
   none.Allocate(KernelFixture::blocks);
   const std::vector<unsigned char> hostNone(KernelFixture::blocks, 0);
   ASSERT_GPU(GPU_MEMCPY(none.data(), hostNone.data(), hostNone.size(),
                         GPU_MEMCPY_HOST_TO_DEVICE));
   runtime.publishProposedState(atomDirection.data(), {}, true, none.data());
   require(runtime.transitionsPublishedLastCall() == 0,
           "rejecting every proposed transition must publish exactly zero transitions");
   require(downloadBlockState() == blockStateBefore,
           "rejected proposals must not change blockState");

   // Accepting every proposal must publish exactly as many transitions as an
   // independent before/after blockState diff shows, and that diff must be
   // nonzero or this scenario is not exercising the counter at all.
   runtime.proposeSelectorState(policy);
   GpuTensor<unsigned char, 1> all;
   all.Allocate(KernelFixture::blocks);
   const std::vector<unsigned char> hostAll(KernelFixture::blocks, 1);
   ASSERT_GPU(GPU_MEMCPY(all.data(), hostAll.data(), hostAll.size(),
                         GPU_MEMCPY_HOST_TO_DEVICE));
   runtime.publishProposedState(atomDirection.data(), {}, true, all.data());
   const auto blockStateAfterAccept = downloadBlockState();
   std::size_t actualDiffs = 0;
   for(std::size_t block = 0; block < KernelFixture::blocks; ++block)
      if(blockStateAfterAccept[block] != blockStateBefore[block]) ++actualDiffs;
   require(actualDiffs > 0,
           "test scenario must exercise at least one real transition");
   require(runtime.transitionsPublishedLastCall() == actualDiffs,
           "transitionsPublishedLastCall() must exactly match the number of "
           "blocks whose state actually changed");

   // materializeCoarseAtoms() must reproduce synchronizeAtomicState()
   // bitwise, and credit only the reason it was called with.
   require(runtime.materializationCount(GpuAdaptiveMaterializationReason::OrdinaryStep) == 0 &&
           runtime.materializationCount(GpuAdaptiveMaterializationReason::Transition) == 0,
           "materialization counters must start at zero");
   GpuTensor<real, 1> viaHelper, viaDirect;
   viaHelper.Allocate(static_cast<index_t>(atomVectors));
   viaDirect.Allocate(static_cast<index_t>(atomVectors));
   viaHelper.copy_async(atomDirection, runtime.stream());
   viaDirect.copy_async(atomDirection, runtime.stream());
   ASSERT_GPU(GPU_STREAM_SYNC(runtime.stream()));
   runtime.materializeCoarseAtoms(GpuAdaptiveMaterializationReason::OrdinaryStep,
                                  viaHelper.data(), {});
   runtime.synchronizeAtomicState(viaDirect.data(), {});
   const auto helperResult = download(viaHelper.data(), atomVectors);
   const auto directResult = download(viaDirect.data(), atomVectors);
   for(std::size_t i = 0; i < atomVectors; ++i)
      require(helperResult[i] == directResult[i],
              "materializeCoarseAtoms must reproduce synchronizeAtomicState bitwise");
   require(runtime.materializationCount(GpuAdaptiveMaterializationReason::OrdinaryStep) == 1 &&
           runtime.materializationCount(GpuAdaptiveMaterializationReason::Transition) == 0,
           "materializeCoarseAtoms(OrdinaryStep, ...) must credit only OrdinaryStep");
   runtime.materializeCoarseAtoms(GpuAdaptiveMaterializationReason::Transition,
                                  viaHelper.data(), {});
   require(runtime.materializationCount(GpuAdaptiveMaterializationReason::Transition) == 1,
           "materializeCoarseAtoms(Transition, ...) must credit Transition");

   viaDirect.Free();
   viaHelper.Free();
   all.Free();
   none.Free();
   atomDirection.Free();
   runtime.release();
}

// CGP-03 Negative control (Part 6): a fine/buffer -> coarse transition is
// reclassified by publishAdaptiveState (blockState flips, transitionsPublished
// is incremented) but that kernel never touches atomDirection in this
// direction -- only the coarse -> fine branch writes atoms directly (see
// gpuAdaptiveRuntime.cpp, the `oldState == coarseState` guard around the
// reconstruction block). This proves gpuSimulation.cpp's
// `if(transitionsPublishedLastCall() > 0)` guard around the post-selector
// materializeCoarseAtoms(Transition, ...) call is load-bearing, not merely an
// optimization: skipping it after a real fine->coarse transition leaves
// stale, unaligned atom directions for the newly-coarse block, exactly the
// parity failure CGP-03's own Part 6 negative control asks to demonstrate.
void testTransitionMaterializationNegativeControl() {
   KernelFixture fixture;
   const auto topology = fixture.topology();
   const auto input = fixture.runtime();
   GpuAdaptiveRuntime runtime;
   runtime.initialize(topology, input, KernelFixture::atoms,
                      KernelFixture::ensembles);

   const std::size_t atomVectors = 3 * KernelFixture::atoms;
   GpuTensor<real, 1> atomDirection;
   atomDirection.Allocate(static_cast<index_t>(atomVectors));
   // +x, orthogonal to the fixture's uniformly +z coarseDirection (see
   // KernelFixture's constructor), so "still at its pre-transition value"
   // and "aligned to coarseDirection" are numerically distinguishable.
   std::vector<double> hostAtom(atomVectors, 0.0);
   for(std::size_t atom = 0; atom < KernelFixture::atoms; ++atom)
      hostAtom[0 + 3 * atom] = 1.0;
   upload(atomDirection.data(), deviceVector(hostAtom));

   // coarsenThreshold above every attainable score and refineThreshold
   // likewise unreachable: every already-atomistic block (1: buffer, 2: fine)
   // proposes coarseState unconditionally; every already-coarse block (0, 3)
   // proposes to stay coarse. Accepting only block 2 (0-indexed, the fine
   // block) isolates a single fine->coarse transition.
   runtime.evaluateSelectorScores(atomDirection.data());
   GpuAdaptiveSelectorPolicy policy;
   policy.refineThreshold = real(2.0);
   policy.coarsenThreshold = real(2.0);
   runtime.proposeSelectorState(policy);
   GpuTensor<unsigned char, 1> onlyBlock2;
   onlyBlock2.Allocate(KernelFixture::blocks);
   const std::vector<unsigned char> hostMask = {0, 0, 1, 0};
   ASSERT_GPU(GPU_MEMCPY(onlyBlock2.data(), hostMask.data(), hostMask.size(),
                         GPU_MEMCPY_HOST_TO_DEVICE));
   runtime.publishProposedState(atomDirection.data(), {}, true, onlyBlock2.data());
   require(runtime.transitionsPublishedLastCall() == 1,
           "isolated single-block acceptance must publish exactly one transition");

   // Block 2 (0-indexed) owns atoms 5,6 (1-based) = indices 4,5 (0-based).
   const auto afterPublish = download(atomDirection.data(), atomVectors);
   const double tolerance = sizeof(real) == sizeof(double) ? 1.0e-12 : 1.0e-5;
   for(const std::size_t atom : {std::size_t(4), std::size_t(5)}) {
      require(std::abs(static_cast<double>(afterPublish[0 + 3 * atom]) - 1.0) < tolerance &&
              std::abs(static_cast<double>(afterPublish[2 + 3 * atom])) < tolerance,
              "publishAdaptiveState must NOT reconstruct a fine->coarse "
              "block's atoms directly -- if this fails, the follow-up "
              "materializeCoarseAtoms(Transition, ...) call in "
              "gpuSimulation.cpp would be provably redundant instead of "
              "load-bearing, and this negative control would need "
              "rewriting");
   }

   // Skipping materialization here is exactly what gpuSimulation.cpp's
   // transitionsPublishedLastCall()==0 fast path would (correctly) never do
   // when a real transition occurred; this control proves why: without it,
   // the newly-coarse block's atoms stay at their stale pre-transition
   // value, not aligned to coarseDirection ((0,0,1) per the fixture).
   runtime.materializeCoarseAtoms(GpuAdaptiveMaterializationReason::Transition,
                                  atomDirection.data(), {});
   const auto afterMaterialize = download(atomDirection.data(), atomVectors);
   for(const std::size_t atom : {std::size_t(4), std::size_t(5)}) {
      require(std::abs(static_cast<double>(afterMaterialize[0 + 3 * atom])) < tolerance &&
              std::abs(static_cast<double>(afterMaterialize[2 + 3 * atom]) - 1.0) < tolerance,
              "materializeCoarseAtoms(Transition, ...) must align the "
              "newly-coarse block's atoms to coarseDirection");
   }

   onlyBlock2.Free();
   atomDirection.Free();
   runtime.release();
}

} // namespace

int main() {
   testSelectorPolicyDescriptorLayout();

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
   // RCG-09C: the live lists are a precondition of every field launch, so the
   // structural check runs before the numerical one; a compaction defect then
   // reports itself as a wrong list rather than as a wrong energy.
   testLiveBondCompaction();
   testKernelParityAndWorkflow();
   testFieldOnlyEvaluation();
   testContinuumOperatorAgainstHostReference();
   testPolarizationGate();
   testHierarchicalCompaction();
   testSynchronizeAtomicStatePolicies();
   testMultiEnsembleGhostProlongation();
   testGhostListNegativeControl();
   testTransitionCounterAndMaterialization();
   testTransitionMaterializationNegativeControl();
   std::cout << "CG-09/CG-10 GPU adaptive runtime tests passed\n";
   return 0;
}
