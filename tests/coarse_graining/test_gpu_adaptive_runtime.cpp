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

void expectEqual(const std::vector<int>& actual, std::initializer_list<int> expected,
                 const char* message) {
   require(actual == std::vector<int>(expected), message);
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
   dilatedPolicy.bufferDilationBlocks = 1;
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
   runtime.integrateHeun(real(1.0e-3), atomDirection.data(),
                         externalField.data(), dipoleField.data());
   const auto afterStep = runtime.downloadWorkSnapshot();
   require(beforeStep.activeAtoms == afterStep.activeAtoms &&
           beforeStep.activeBlocks == afterStep.activeBlocks,
           "Heun integration changed adaptive state inside a stage");
   const auto stepped = download(atomDirection.data(), atomVectors);
   for(const int oneBasedAtom : afterStep.activeAtoms) {
      const std::size_t atom = static_cast<std::size_t>(oneBasedAtom - 1);
      double norm2 = 0.0;
      for(int xyz = 0; xyz < 3; ++xyz)
         norm2 += static_cast<double>(stepped[xyz + 3 * atom]) *
                  static_cast<double>(stepped[xyz + 3 * atom]);
      require(std::abs(norm2 - 1.0) < (sizeof(real) == sizeof(double) ? 2.0e-12 : 2.0e-5),
              "GPU Heun integration did not preserve unit directions");
   }

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

   atomDirection.Free();
   runtime.release();
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
   testKernelParityAndWorkflow();
   testPolarizationGate();
   std::cout << "CG-09/CG-10 GPU adaptive runtime tests passed\n";
   return 0;
}
