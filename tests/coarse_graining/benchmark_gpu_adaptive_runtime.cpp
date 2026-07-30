#include "gpuAdaptiveRuntime.hpp"
#include "measurement/memoryMeasurement.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <numeric>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

constexpr int coarseState = 0;
constexpr int bufferState = 1;
constexpr int fineState = 2;

struct Options {
   std::size_t blocks = 2048;
   std::size_t atomsPerBlock = 4;
   unsigned int warmup = 2;
   unsigned int iterations = 10;
   unsigned int repetitions = 7;
   unsigned int featureIterations = 100;
   unsigned int featureRounds = 64;
   double featureNoiseFraction = 0.03;
   double crossoverMargin = 0.02;
   bool requireAcceptance = false;
};

std::size_t parseSize(const char* value, const char* option) {
   char* end = nullptr;
   const unsigned long long parsed = std::strtoull(value, &end, 10);
   if(!value[0] || !end || *end || parsed == 0)
      throw std::runtime_error(std::string("invalid ") + option);
   return static_cast<std::size_t>(parsed);
}

double parsePercent(const char* value, const char* option) {
   char* end = nullptr;
   const double parsed = std::strtod(value, &end);
   if(!value[0] || !end || *end || !std::isfinite(parsed) ||
      parsed <= 0.0 || parsed >= 100.0)
      throw std::runtime_error(std::string("invalid ") + option);
   return parsed / 100.0;
}

Options parse(int argc, char** argv) {
   Options result;
   for(int index = 1; index < argc; ++index) {
      const std::string option = argv[index];
      if(option == "--blocks" && index + 1 < argc)
         result.blocks = parseSize(argv[++index], "--blocks");
      else if(option == "--atoms-per-block" && index + 1 < argc)
         result.atomsPerBlock = parseSize(argv[++index], "--atoms-per-block");
      else if(option == "--warmup" && index + 1 < argc)
         result.warmup = static_cast<unsigned int>(parseSize(argv[++index], "--warmup"));
      else if(option == "--iterations" && index + 1 < argc)
         result.iterations = static_cast<unsigned int>(parseSize(argv[++index], "--iterations"));
      else if(option == "--repetitions" && index + 1 < argc)
         result.repetitions = static_cast<unsigned int>(parseSize(argv[++index], "--repetitions"));
      else if(option == "--feature-iterations" && index + 1 < argc)
         result.featureIterations =
            static_cast<unsigned int>(parseSize(argv[++index], "--feature-iterations"));
      else if(option == "--feature-rounds" && index + 1 < argc)
         result.featureRounds =
            static_cast<unsigned int>(parseSize(argv[++index], "--feature-rounds"));
      else if(option == "--feature-noise-percent" && index + 1 < argc)
         result.featureNoiseFraction =
            parsePercent(argv[++index], "--feature-noise-percent");
      else if(option == "--crossover-margin-percent" && index + 1 < argc)
         result.crossoverMargin =
            parsePercent(argv[++index], "--crossover-margin-percent");
      else if(option == "--require-acceptance")
         result.requireAcceptance = true;
      else if(option == "--help") {
         std::puts(
            "usage: gpu_adaptive_runtime_benchmark [--blocks N] "
            "[--atoms-per-block N] [--warmup N] [--iterations N] "
            "[--repetitions N] [--feature-iterations N] "
            "[--feature-rounds N] [--feature-noise-percent P] "
            "[--crossover-margin-percent P] [--require-acceptance]");
         std::exit(0);
      } else {
         throw std::runtime_error("unknown or incomplete benchmark option: " + option);
      }
   }
   if(result.blocks < 8 || result.atomsPerBlock < 2 || result.repetitions < 3)
      throw std::runtime_error(
         "benchmark requires at least 8 blocks, 2 atoms/block, and 3 repetitions");
   return result;
}

void gpuCheck(GPU_ERROR_T status, const char* operation) {
   if(status != GPU_SUCCESS)
      throw std::runtime_error(std::string(operation) + ": " +
                               GPU_GET_ERROR_STRING(status));
}

double median(std::vector<double> values) {
   if(values.empty()) return 0.0;
   std::sort(values.begin(), values.end());
   const std::size_t middle = values.size() / 2;
   return values.size() % 2 ? values[middle] :
      0.5 * (values[middle - 1] + values[middle]);
}

double medianAbsoluteDeviation(const std::vector<double>& values) {
   const double center = median(values);
   std::vector<double> deviations(values.size());
   std::transform(values.begin(), values.end(), deviations.begin(),
                  [center](double value) { return std::abs(value - center); });
   return median(std::move(deviations));
}

__global__ void featureOffAtomisticWork(real* values, std::size_t count,
                                        unsigned int rounds) {
   const std::size_t index =
      static_cast<std::size_t>(blockIdx.x) * blockDim.x + threadIdx.x;
   if(index >= count) return;
   real value = values[index];
   for(unsigned int round = 0; round < rounds; ++round) {
      value = value * real(0.9999997) +
              real(0.0000001) * static_cast<real>((index + round) % 17);
   }
   values[index] = value;
}

void launchFeatureOffWork(real* values, std::size_t count, unsigned int rounds,
                          GPU_STREAM_T stream,
                          const GpuAdaptiveRuntime* inactiveRuntime) {
   if(inactiveRuntime && (inactiveRuntime->ready() ||
                          inactiveRuntime->allocatedBytes() != 0))
      throw std::logic_error("feature-off control unexpectedly owns GPU state");
   const unsigned int threads = 256;
   const unsigned int blocks =
      static_cast<unsigned int>((count + threads - 1) / threads);
#if defined(CUDA_V)
   featureOffAtomisticWork<<<blocks, threads, 0, stream>>>(values, count, rounds);
#else
   hipLaunchKernelGGL(featureOffAtomisticWork, dim3(blocks), dim3(threads), 0,
                      stream, values, count, rounds);
#endif
   gpuCheck(GPU_GET_LAST_ERROR(), "feature-off kernel launch");
}

double timeFeatureOffBatch(real* values, std::size_t count,
                           const Options& options, GPU_STREAM_T stream,
                           const GpuAdaptiveRuntime* inactiveRuntime) {
   GPU_EVENT_T begin{}, end{};
   gpuCheck(GPU_EVENT_CREATE(&begin), "feature-off begin-event creation");
   try {
      gpuCheck(GPU_EVENT_CREATE(&end), "feature-off end-event creation");
      gpuCheck(GPU_EVENT_RECORD(begin, stream), "feature-off begin-event record");
      for(unsigned int iteration = 0; iteration < options.featureIterations; ++iteration)
         launchFeatureOffWork(values, count, options.featureRounds, stream,
                              inactiveRuntime);
      gpuCheck(GPU_EVENT_RECORD(end, stream), "feature-off end-event record");
      gpuCheck(GPU_EVENT_SYNCHRONIZE(end), "feature-off event synchronization");
      float elapsed = 0.0f;
      gpuCheck(GPU_EVENT_ELAPSED_TIME(&elapsed, begin, end),
               "feature-off elapsed time");
      GPU_EVENT_DESTROY(begin);
      GPU_EVENT_DESTROY(end);
      return 1000.0 * static_cast<double>(elapsed) /
             options.featureIterations;
   } catch(...) {
      if(begin) GPU_EVENT_DESTROY(begin);
      if(end) GPU_EVENT_DESTROY(end);
      throw;
   }
}

struct AdaptiveFixture {
   std::size_t atoms = 0;
   std::size_t blocks = 0;
   std::size_t basis = 0;
   std::size_t bonds = 0;
   int repetitionShape[3] = {};
   int blockShape[3] = {1, 1, 1};
   int blockGrid[3] = {};
   double cellVectors[9] = {1, 0, 0, 0, 1, 0, 0, 0, 1};
   double blockVectors[9] = {1, 0, 0, 0, 1, 0, 0, 0, 1};
   std::vector<int> atomToBlock, atomToBasis, atomToDynamic, atomToFft;
   std::vector<int> atomToFftGrid, basisToDynamic, basisToFft;
   std::vector<int> blockCount, blockOffset, blockAtoms, blockCoordinate;
   std::vector<int> basisPopulation, fftPopulation, dynamicPopulation;
   std::vector<double> center, volume;
   std::vector<int> state;
   std::vector<double> scores, atomMoment;
   std::vector<int> atomAxisCount;
   std::vector<double> atomAxis, atomK1, atomK2;
   std::vector<int> projectionBlock, bondAtom, selectorEdge;
   std::vector<double> projectionWeight, bondMatrix;
   double inverseBlockTranspose[9] = {1, 0, 0, 0, 1, 0, 0, 0, 1};
   double exchange[9] = {};
   double spiralization[9] = {};
   std::vector<int> axisCount;
   std::vector<double> axis, k1, k2;
   std::vector<double> coarseMoment, coarseDirection, coarseField, momentSum;
   std::vector<real> atomDirection;

   AdaptiveFixture(std::size_t blockCountValue, std::size_t atomsPerBlock)
      : atoms(blockCountValue * atomsPerBlock),
        blocks(blockCountValue),
        basis(atomsPerBlock),
        atomToBlock(atoms), atomToBasis(atoms), atomToDynamic(atoms, 1),
        atomToFft(atoms), atomToFftGrid(atoms),
        basisToDynamic(basis, 1), basisToFft(basis),
        blockCount(blocks, static_cast<int>(atomsPerBlock)),
        blockOffset(blocks + 1), blockAtoms(atoms),
        blockCoordinate(3 * blocks), basisPopulation(basis * blocks, 1),
        fftPopulation(basis * blocks, 1),
        dynamicPopulation(blocks, static_cast<int>(atomsPerBlock)),
        center(3 * blocks), volume(blocks, 1.0), state(blocks, fineState),
        scores(blocks), atomMoment(atoms, 1.0), atomAxisCount(atoms),
        atomAxis(6 * atoms), atomK1(2 * atoms), atomK2(2 * atoms),
        projectionBlock(8 * atoms), projectionWeight(8 * atoms),
        axisCount(blocks), axis(6 * blocks), k1(2 * blocks),
        k2(2 * blocks), coarseMoment(3 * blocks),
        coarseDirection(3 * blocks), coarseField(3 * blocks),
        momentSum(blocks, static_cast<double>(atomsPerBlock)),
        atomDirection(3 * atoms) {
      repetitionShape[0] = static_cast<int>(blocks);
      repetitionShape[1] = repetitionShape[2] = 1;
      blockGrid[0] = static_cast<int>(blocks);
      blockGrid[1] = blockGrid[2] = 1;
      exchange[0] = 0.05;
      for(std::size_t member = 0; member < basis; ++member)
         basisToFft[member] = static_cast<int>(member + 1);
      for(std::size_t block = 0; block < blocks; ++block) {
         blockOffset[block] = static_cast<int>(block * basis);
         blockCoordinate[3 * block] = static_cast<int>(block);
         center[3 * block] = static_cast<double>(block) + 0.5;
         const double phase =
            2.0 * 3.14159265358979323846 * static_cast<double>(block) /
            static_cast<double>(blocks);
         double direction[3] = {std::cos(phase), std::sin(phase), 0.2};
         const double norm = std::sqrt(
            direction[0] * direction[0] + direction[1] * direction[1] +
            direction[2] * direction[2]);
         for(int xyz = 0; xyz < 3; ++xyz) {
            direction[xyz] /= norm;
            coarseDirection[xyz + 3 * block] = direction[xyz];
            coarseMoment[xyz + 3 * block] =
               static_cast<double>(basis) * direction[xyz];
         }
         for(std::size_t member = 0; member < basis; ++member) {
            const std::size_t atom = member + basis * block;
            atomToBlock[atom] = static_cast<int>(block + 1);
            atomToBasis[atom] = static_cast<int>(member + 1);
            atomToFft[atom] = static_cast<int>(member + 1);
            atomToFftGrid[atom] = static_cast<int>(atom + 1);
            blockAtoms[atom] = static_cast<int>(atom + 1);
            projectionBlock[8 * atom] = static_cast<int>(block + 1);
            projectionWeight[8 * atom] = 1.0;
            for(int corner = 1; corner < 8; ++corner)
               projectionBlock[corner + 8 * atom] =
                  static_cast<int>(block + 1);
            for(int xyz = 0; xyz < 3; ++xyz)
               atomDirection[xyz + 3 * atom] =
                  static_cast<real>(direction[xyz]);
         }
      }
      blockOffset[blocks] = static_cast<int>(atoms);

      // One physical within-block chain plus one periodic same-basis bond.
      for(std::size_t block = 0; block < blocks; ++block) {
         for(std::size_t member = 0; member + 1 < basis; ++member) {
            bondAtom.push_back(static_cast<int>(member + basis * block + 1));
            bondAtom.push_back(static_cast<int>(member + 1 + basis * block + 1));
         }
         const std::size_t next = (block + 1) % blocks;
         for(std::size_t member = 0; member < basis; ++member) {
            bondAtom.push_back(static_cast<int>(member + basis * block + 1));
            bondAtom.push_back(static_cast<int>(member + basis * next + 1));
         }
      }
      bonds = bondAtom.size() / 2;
      selectorEdge = bondAtom;
      bondMatrix.assign(9 * bonds, 0.0);
      for(std::size_t bond = 0; bond < bonds; ++bond)
         for(int xyz = 0; xyz < 3; ++xyz)
            bondMatrix[xyz + 3 * xyz + 9 * bond] = 0.1;
   }

   GpuAdaptiveTopologyInput topology() const {
      GpuAdaptiveTopologyInput result;
      result.geometryMode = 1;
      result.atoms = atoms;
      result.blocks = blocks;
      result.basis = basis;
      result.fftChannelsPerBlock = basis;
      result.fftGridChannels = basis * blocks;
      result.dynamicChannels = 1;
      result.ensembles = 1;
      result.repetitionShape = repetitionShape;
      result.blockShape = blockShape;
      result.blockGrid = blockGrid;
      result.cellVectors = cellVectors;
      result.blockVectors = blockVectors;
      result.atomToBlock = atomToBlock.data();
      result.atomToBasis = atomToBasis.data();
      result.atomToDynamicChannel = atomToDynamic.data();
      result.atomToFftChannel = atomToFft.data();
      result.atomToFftGridIndex = atomToFftGrid.data();
      result.basisToDynamicChannel = basisToDynamic.data();
      result.basisToFftChannel = basisToFft.data();
      result.blockAtomCount = blockCount.data();
      result.blockAtomOffset = blockOffset.data();
      result.blockAtoms = blockAtoms.data();
      result.blockGridCoordinate = blockCoordinate.data();
      result.blockBasisPopulation = basisPopulation.data();
      result.blockFftChannelPopulation = fftPopulation.data();
      result.blockDynamicChannelPopulation = dynamicPopulation.data();
      result.blockCenter = center.data();
      result.blockVolume = volume.data();
      return result;
   }

   GpuAdaptiveRuntimeInput runtime() const {
      GpuAdaptiveRuntimeInput result;
      result.blockState = state.data();
      result.selectorCriteria = 1;
      result.selectorScores = scores.data();
      result.coarseMoment = coarseMoment.data();
      result.coarseDirection = coarseDirection.data();
      result.coarseField = coarseField.data();
      result.channelMomentSum = momentSum.data();
      result.kernels.atomMoment = atomMoment.data();
      result.kernels.atomAnisotropyAxisCount = atomAxisCount.data();
      result.kernels.atomAnisotropyAxis = atomAxis.data();
      result.kernels.atomAnisotropyK1 = atomK1.data();
      result.kernels.atomAnisotropyK2 = atomK2.data();
      result.kernels.projectionBlock = projectionBlock.data();
      result.kernels.projectionWeight = projectionWeight.data();
      result.kernels.bonds = bonds;
      result.kernels.bondAtom = bondAtom.data();
      result.kernels.bondMatrix = bondMatrix.data();
      result.kernels.selectorEdges = bonds;
      result.kernels.selectorEdge = selectorEdge.data();
      result.kernels.inverseBlockTranspose = inverseBlockTranspose;
      result.kernels.exchangeStiffness = exchange;
      result.kernels.spiralization = spiralization;
      result.kernels.anisotropyAxisCount = axisCount.data();
      result.kernels.anisotropyAxis = axis.data();
      result.kernels.anisotropyK1 = k1.data();
      result.kernels.anisotropyK2 = k2.data();
      result.kernels.magneticMomentSi = 1.0;
      result.kernels.gammaPerTs = 1.0;
      result.kernels.damping = 0.1;
      return result;
   }

   void setFraction(double fraction) {
      std::fill(state.begin(), state.end(), coarseState);
      const std::size_t fineBlocks = static_cast<std::size_t>(
         std::llround(fraction * static_cast<double>(blocks)));
      for(std::size_t block = 0; block < std::min(fineBlocks, blocks); ++block)
         state[block] = fineState;
      if(fineBlocks > 0 && fineBlocks < blocks) {
         state[fineBlocks] = bufferState;
         state[blocks - 1] = bufferState;
      }
   }
};

struct SweepSample {
   double requestedFraction = 0.0;
   double activeDofRatio = 0.0;
   double activeAtomFraction = 0.0;
   double activeBlockFraction = 0.0;
   double interfaceFraction = 0.0;
   std::size_t activeAtoms = 0;
   std::size_t activeBlocks = 0;
   std::size_t interfaceAtoms = 0;
   std::size_t allocatedBytes = 0;
   std::vector<double> wallUs;
   std::vector<double> atomisticUs;
   std::vector<double> coarseUs;
   std::vector<double> interfaceUs;
   std::vector<double> selectorUs;
   std::vector<double> compactionUs;
   std::vector<double> fftUs;
   std::vector<double> integrationUs;
   std::vector<double> unaccountedUs;
   double wallMedian = 0.0;
   double wallMad = 0.0;
};

SweepSample measureFraction(GpuAdaptiveRuntime& runtime,
                            AdaptiveFixture& fixture, double fraction,
                            const Options& options, const real* atomDirection,
                            real* atomField, real* coarseField) {
   fixture.setFraction(fraction);
   runtime.updateBlockState(fixture.state.data(), fixture.blocks);
   for(unsigned int iteration = 0; iteration < options.warmup; ++iteration)
      (void)runtime.evaluateHybrid(atomDirection, nullptr, nullptr,
                                   atomField, coarseField);
   SweepSample result;
   result.requestedFraction = fraction;
   const auto snapshot = runtime.downloadWorkSnapshot();
   result.activeAtoms = snapshot.activeAtoms.size();
   result.activeBlocks = snapshot.activeBlocks.size();
   result.interfaceAtoms = snapshot.interfaceAtoms.size();
   result.activeAtomFraction = static_cast<double>(result.activeAtoms) /
      static_cast<double>(fixture.atoms);
   result.activeBlockFraction = static_cast<double>(result.activeBlocks) /
      static_cast<double>(fixture.blocks);
   result.interfaceFraction = static_cast<double>(result.interfaceAtoms) /
      static_cast<double>(fixture.atoms);
   result.activeDofRatio =
      static_cast<double>(result.activeAtoms + result.activeBlocks) /
      static_cast<double>(fixture.atoms);
   result.allocatedBytes = runtime.allocatedBytes();
   for(unsigned int repetition = 0; repetition < options.repetitions; ++repetition) {
      runtime.resetPhaseMetrics();
      const auto begin = std::chrono::steady_clock::now();
      for(unsigned int iteration = 0; iteration < options.iterations; ++iteration)
         (void)runtime.evaluateHybrid(atomDirection, nullptr, nullptr,
                                      atomField, coarseField);
      const double wallUs = 1.0e6 *
         std::chrono::duration<double>(std::chrono::steady_clock::now() - begin).count() /
         options.iterations;
      const auto phases = runtime.phaseMetrics();
      result.wallUs.push_back(wallUs);
      result.atomisticUs.push_back(
         1000.0 * phases.atomisticMilliseconds / options.iterations);
      result.coarseUs.push_back(
         1000.0 * phases.coarseMilliseconds / options.iterations);
      result.interfaceUs.push_back(
         1000.0 * phases.interfaceMilliseconds / options.iterations);
      result.selectorUs.push_back(
         1000.0 * phases.selectorMilliseconds / options.iterations);
      result.compactionUs.push_back(
         1000.0 * phases.compactionMilliseconds / options.iterations);
      result.fftUs.push_back(
         1000.0 * phases.fftMilliseconds / options.iterations);
      result.integrationUs.push_back(
         1000.0 * phases.integrationMilliseconds / options.iterations);
      result.unaccountedUs.push_back(
         wallUs - 1000.0 * (phases.atomisticMilliseconds +
            phases.coarseMilliseconds + phases.interfaceMilliseconds +
            phases.selectorMilliseconds + phases.compactionMilliseconds +
            phases.fftMilliseconds + phases.integrationMilliseconds) /
            options.iterations);
   }
   result.wallMedian = median(result.wallUs);
   result.wallMad = medianAbsoluteDeviation(result.wallUs);
   return result;
}

void printSweep(const SweepSample& sample) {
   std::printf(
      "adaptive-sweep requested_fraction=%.6f active_atoms=%zu "
      "active_blocks=%zu interface_atoms=%zu active_atom_fraction=%.6f "
      "active_block_fraction=%.6f interface_fraction=%.6f active_dof_ratio=%.6f "
      "wall_us=%.6f wall_mad_us=%.6f atomistic_us=%.6f coarse_us=%.6f "
      "interface_us=%.6f selector_us=%.6f compaction_us=%.6f fft_us=%.6f "
      "integration_us=%.6f unaccounted_us=%.6f device_bytes=%zu\n",
      sample.requestedFraction, sample.activeAtoms, sample.activeBlocks,
      sample.interfaceAtoms, sample.activeAtomFraction,
      sample.activeBlockFraction, sample.interfaceFraction,
      sample.activeDofRatio, sample.wallMedian, sample.wallMad,
      median(sample.atomisticUs), median(sample.coarseUs),
      median(sample.interfaceUs), median(sample.selectorUs),
      median(sample.compactionUs), median(sample.fftUs),
      median(sample.integrationUs), median(sample.unaccountedUs),
      sample.allocatedBytes);
}

} // namespace

int main(int argc, char** argv) {
   try {
      const Options options = parse(argc, argv);
      std::printf(
         "adaptive-benchmark precision=%s backend=%s blocks=%zu "
         "atoms_per_block=%zu warmup=%u iterations=%u repetitions=%u\n",
         sizeof(real) == sizeof(double) ? "fp64" : "fp32",
#if defined(CUDA_V)
         "CUDA",
#else
         "HIP",
#endif
         options.blocks, options.atomsPerBlock, options.warmup,
         options.iterations, options.repetitions);

      // Paired feature-off control.  Both samples execute the same atomistic
      // kernel; the candidate sample merely keeps the default, uninitialized
      // adaptive owner in scope.  Alternating order suppresses thermal drift.
      const std::size_t featureCount = 3 * options.blocks * options.atomsPerBlock;
      GpuTensor<real, 1> featureValues;
      featureValues.Allocate(static_cast<index_t>(featureCount));
      featureValues.zeros();
      GPU_STREAM_T featureStream{};
      gpuCheck(GPU_STREAM_CREATE(&featureStream), "feature-off stream creation");
      std::vector<double> featureBaseline, featureCandidate;
      const auto featureInventory = TensorMemoryTracker::current_device();
      {
         GpuAdaptiveRuntime featureOff;
         if(featureOff.ready() || featureOff.allocatedBytes() != 0 ||
            TensorMemoryTracker::current_device() != featureInventory)
            throw std::runtime_error("feature-off runtime changed the device inventory");
         for(unsigned int repetition = 0; repetition < options.repetitions; ++repetition) {
            if(repetition % 2 == 0) {
               featureBaseline.push_back(timeFeatureOffBatch(
                  featureValues.data(), featureCount, options, featureStream, nullptr));
               featureCandidate.push_back(timeFeatureOffBatch(
                  featureValues.data(), featureCount, options, featureStream, &featureOff));
            } else {
               featureCandidate.push_back(timeFeatureOffBatch(
                  featureValues.data(), featureCount, options, featureStream, &featureOff));
               featureBaseline.push_back(timeFeatureOffBatch(
                  featureValues.data(), featureCount, options, featureStream, nullptr));
            }
         }
      }
      gpuCheck(GPU_STREAM_DESTROY(featureStream), "feature-off stream destruction");
      featureValues.Free();
      const double featureBaselineMedian = median(featureBaseline);
      const double featureCandidateMedian = median(featureCandidate);
      const double featureDelta =
         (featureCandidateMedian - featureBaselineMedian) / featureBaselineMedian;
      const double featureNoise = 3.0 *
         (medianAbsoluteDeviation(featureBaseline) +
          medianAbsoluteDeviation(featureCandidate));
      const double featureBudget = std::max(
         options.featureNoiseFraction * featureBaselineMedian, featureNoise);
      const bool featureAccepted =
         std::abs(featureCandidateMedian - featureBaselineMedian) <= featureBudget;
      std::printf(
         "feature-off baseline_us=%.6f candidate_us=%.6f delta_percent=%.6f "
         "baseline_mad_us=%.6f candidate_mad_us=%.6f budget_us=%.6f "
         "inventory_delta_bytes=0 result=%s\n",
         featureBaselineMedian, featureCandidateMedian, 100.0 * featureDelta,
         medianAbsoluteDeviation(featureBaseline),
         medianAbsoluteDeviation(featureCandidate), featureBudget,
         featureAccepted ? "PASS" : "NOISY_OR_REGRESSED");

      AdaptiveFixture fixture(options.blocks, options.atomsPerBlock);
      const auto topology = fixture.topology();
      const auto input = fixture.runtime();
      std::string diagnostic;
      if(!GpuAdaptiveRuntime::validate(topology, input, fixture.atoms, 1, diagnostic))
         throw std::runtime_error("benchmark fixture rejected: " + diagnostic);
      GpuAdaptiveRuntime runtime;
      runtime.initialize(topology, input, fixture.atoms, 1);
      GpuTensor<real, 1> atomDirection, atomField, coarseField;
      atomDirection.Allocate(static_cast<index_t>(3 * fixture.atoms));
      atomField.Allocate(static_cast<index_t>(3 * fixture.atoms));
      coarseField.Allocate(static_cast<index_t>(3 * fixture.blocks));
      gpuCheck(GPU_MEMCPY(atomDirection.data(), fixture.atomDirection.data(),
                          fixture.atomDirection.size() * sizeof(real),
                          GPU_MEMCPY_HOST_TO_DEVICE),
               "benchmark direction upload");

      const std::vector<double> fractions = {
         1.0, 0.75, 0.5, 0.25, 0.125, 0.0625, 0.0
      };
      std::vector<SweepSample> sweep;
      for(const double fraction : fractions) {
         sweep.push_back(measureFraction(
            runtime, fixture, fraction, options, atomDirection.data(),
            atomField.data(), coarseField.data()));
         printSweep(sweep.back());
      }

      // Selector and compaction are measured at the mixed 50% point.  Report
      // both device-event time and the mandatory localized host wait.
      fixture.setFraction(0.5);
      runtime.updateBlockState(fixture.state.data(), fixture.blocks);
      GpuAdaptiveSelectorPolicy selectorPolicy;
      selectorPolicy.refineThreshold = real(0.25);
      selectorPolicy.coarsenThreshold = real(0.10);
      runtime.resetPhaseMetrics();
      const auto selectorBegin = std::chrono::steady_clock::now();
      for(unsigned int iteration = 0; iteration < options.iterations; ++iteration) {
         runtime.evaluateSelectorScores(atomDirection.data());
         runtime.proposeSelectorState(selectorPolicy);
      }
      const double selectorWallUs = 1.0e6 *
         std::chrono::duration<double>(
            std::chrono::steady_clock::now() - selectorBegin).count() /
         options.iterations;
      const double selectorDeviceUs =
         1000.0 * runtime.phaseMetrics().selectorMilliseconds /
         options.iterations;

      const auto compactionBefore = runtime.compactionMetrics();
      const auto compactionWallBegin = std::chrono::steady_clock::now();
      for(unsigned int iteration = 0; iteration < options.iterations; ++iteration) {
         fixture.setFraction(iteration % 2 ? 0.5 : 0.25);
         runtime.updateBlockState(fixture.state.data(), fixture.blocks);
      }
      const double compactionWallUs = 1.0e6 *
         std::chrono::duration<double>(
            std::chrono::steady_clock::now() - compactionWallBegin).count() /
         options.iterations;
      const auto compactionAfter = runtime.compactionMetrics();
      const double compactionDeviceUs = 1000.0 *
         (compactionAfter.deviceMilliseconds - compactionBefore.deviceMilliseconds) /
         options.iterations;
      const double compactionHostWaitUs = 1000.0 *
         (compactionAfter.hostWaitMilliseconds -
          compactionBefore.hostWaitMilliseconds) / options.iterations;
      const double mixedWallUs = sweep[2].wallMedian;
      std::printf(
         "adaptive-overhead selector_device_us=%.6f selector_wall_us=%.6f "
         "compaction_device_us=%.6f compaction_host_wait_us=%.6f "
         "compaction_wall_us=%.6f block_bytes_per_update=%zu "
         "device_bytes=%zu mixed_unaccounted_us=%.6f "
         "overhead_percent_of_mixed_step=%.6f\n",
         selectorDeviceUs, selectorWallUs, compactionDeviceUs,
         compactionHostWaitUs, compactionWallUs,
         fixture.blocks * sizeof(int),
         runtime.allocatedBytes(), sweep[2].unaccountedUs.empty() ? 0.0 :
            median(sweep[2].unaccountedUs),
         100.0 * (selectorWallUs + compactionWallUs) / mixedWallUs);

      // The fraction sweep is a field-evaluation benchmark.  Also measure a
      // complete two-stage Heun call so integration/reconstruction overhead is
      // visible in the same artifact instead of being inferred from speedup.
      for(unsigned int iteration = 0; iteration < options.warmup; ++iteration)
         runtime.integrateHeun(real(1.0e-6), atomDirection.data());
      runtime.resetPhaseMetrics();
      const auto stepBegin = std::chrono::steady_clock::now();
      for(unsigned int iteration = 0; iteration < options.iterations; ++iteration)
         runtime.integrateHeun(real(1.0e-6), atomDirection.data());
      const double stepWallUs = 1.0e6 *
         std::chrono::duration<double>(
            std::chrono::steady_clock::now() - stepBegin).count() /
         options.iterations;
      const auto stepPhase = runtime.phaseMetrics();
      const double stepDeviceUs = 1000.0 *
         (stepPhase.atomisticMilliseconds + stepPhase.coarseMilliseconds +
          stepPhase.interfaceMilliseconds + stepPhase.selectorMilliseconds +
          stepPhase.compactionMilliseconds + stepPhase.fftMilliseconds +
          stepPhase.integrationMilliseconds) / options.iterations;
      std::printf(
         "adaptive-step wall_us=%.6f phase_us(atomistic=%.6f coarse=%.6f "
         "interface=%.6f selector=%.6f compaction=%.6f fft=%.6f "
         "integration=%.6f) device_phase_sum_us=%.6f unaccounted_us=%.6f\n",
         stepWallUs, 1000.0 * stepPhase.atomisticMilliseconds / options.iterations,
         1000.0 * stepPhase.coarseMilliseconds / options.iterations,
         1000.0 * stepPhase.interfaceMilliseconds / options.iterations,
         1000.0 * stepPhase.selectorMilliseconds / options.iterations,
         1000.0 * stepPhase.compactionMilliseconds / options.iterations,
         1000.0 * stepPhase.fftMilliseconds / options.iterations,
         1000.0 * stepPhase.integrationMilliseconds / options.iterations,
         stepDeviceUs, stepWallUs - stepDeviceUs);

      const SweepSample& atomistic = sweep.front();
      const SweepSample* crossover = nullptr;
      for(std::size_t index = 1; index < sweep.size(); ++index) {
         const double uncertainty = 3.0 *
            (atomistic.wallMad + sweep[index].wallMad);
         const double target =
            atomistic.wallMedian * (1.0 - options.crossoverMargin);
         if(sweep[index].wallMedian + uncertainty < target) {
            crossover = &sweep[index];
            break;
         }
      }
      if(crossover) {
         std::printf(
            "active-dof-crossover result=PASS active_dof_ratio=%.6f "
            "requested_fraction=%.6f wall_us=%.6f atomistic_wall_us=%.6f "
            "speedup=%.6f acceptance_margin_percent=%.6f\n",
            crossover->activeDofRatio, crossover->requestedFraction,
            crossover->wallMedian, atomistic.wallMedian,
            atomistic.wallMedian / crossover->wallMedian,
            100.0 * options.crossoverMargin);
      } else {
         std::printf(
            "active-dof-crossover result=NOT_OBSERVED atomistic_wall_us=%.6f "
            "acceptance_margin_percent=%.6f\n",
            atomistic.wallMedian, 100.0 * options.crossoverMargin);
      }

      coarseField.Free();
      atomField.Free();
      atomDirection.Free();
      runtime.release();
      if(options.requireAcceptance && (!featureAccepted || !crossover))
         return 2;
      return 0;
   } catch(const std::exception& error) {
      std::fprintf(stderr, "FAIL GPU adaptive runtime benchmark: %s\n",
                   error.what());
      return 1;
   }
}
