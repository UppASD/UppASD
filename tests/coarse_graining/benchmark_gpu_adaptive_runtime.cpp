// RCG-09: PERF-ATOMISTIC-PROD and PERF-CG-SWEEP.
//
// This benchmark compares adaptive coarse graining against UppASD's real
// production atomistic GPU path -- GpuHamiltonianCalculations (the same
// `heisge` the feature-off timestep loop calls) plus GpuDepondtIntegrator (the
// same SDEalgh=1 evolveFirst/evolveSecond pair) -- on one shared geometry, at
// one precision, with one ensemble count and one timestep.
//
// Finding F-10 was that the previous version of this file compared adaptive
// coarse graining against a synthetic fused-multiply-add loop, and that its
// "crossover" compared the adaptive runtime at a low fine fraction against the
// adaptive runtime at 100% fine -- two modes of the same code, never against
// UppASD.  Both defects are corrected here:
//
//   * the synthetic loop moved to
//     tests/coarse_graining/benchmark_gpu_inactive_runtime_overhead.cpp, where
//     it is labelled as an inactive-runtime overhead control and cannot be
//     mistaken for a baseline;
//   * the crossover below is computed against the production step, and the
//     all-fine adaptive point is reported alongside it as a diagnostic only.
//
// Before any timing, the two paths are required to agree numerically on the
// same state (`atomistic-parity` below).  That check is what licenses the
// comparison: it is how "identical supported physics" is established as a
// measurement rather than asserted as a claim.  If it fails, the benchmark
// fails and reports no speedup.

#include "gpuAdaptiveRuntime.hpp"
#include "gpuHamiltonianCalculations.hpp"
#include "gpuDepondtIntegrator.hpp"
#include "gpuParallelizationHelper.hpp"
#include "gpuStructures.hpp"
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

// One isotropic exchange constant shared by the production neighbour list and
// the adaptive bond matrix.  Identical physics starts here.
constexpr double exchangeJ = 0.1;

enum class Texture { Spiral, DomainWall, Uniform };

struct Options {
   std::size_t blocks = 2048;
   std::size_t atomsPerBlock = 4;
   unsigned int warmup = 2;
   unsigned int iterations = 10;
   unsigned int repetitions = 7;
   double crossoverMargin = 0.02;
   double parityTolerance = 0.0;  // 0 => precision-derived default
   Texture texture = Texture::Spiral;
   bool requireAcceptance = false;
   bool requireCrossover = false;
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

double parsePositive(const char* value, const char* option) {
   char* end = nullptr;
   const double parsed = std::strtod(value, &end);
   if(!value[0] || !end || *end || !std::isfinite(parsed) || parsed <= 0.0)
      throw std::runtime_error(std::string("invalid ") + option);
   return parsed;
}

const char* textureName(Texture texture) {
   switch(texture) {
      case Texture::Spiral: return "spiral";
      case Texture::DomainWall: return "domain_wall";
      case Texture::Uniform: return "uniform";
   }
   return "unknown";
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
      else if(option == "--crossover-margin-percent" && index + 1 < argc)
         result.crossoverMargin = parsePercent(argv[++index], "--crossover-margin-percent");
      else if(option == "--parity-tolerance" && index + 1 < argc)
         result.parityTolerance = parsePositive(argv[++index], "--parity-tolerance");
      else if(option == "--texture" && index + 1 < argc) {
         const std::string value = argv[++index];
         if(value == "spiral") result.texture = Texture::Spiral;
         else if(value == "domain-wall") result.texture = Texture::DomainWall;
         else if(value == "uniform") result.texture = Texture::Uniform;
         else throw std::runtime_error("unknown --texture: " + value);
      }
      else if(option == "--require-acceptance")
         result.requireAcceptance = true;
      else if(option == "--require-crossover")
         result.requireCrossover = true;
      else if(option == "--help") {
         std::puts(
            "usage: gpu_adaptive_runtime_benchmark [--blocks N] "
            "[--atoms-per-block N] [--warmup N] [--iterations N] "
            "[--repetitions N] [--crossover-margin-percent P] "
            "[--parity-tolerance T] [--texture spiral|domain-wall|uniform] "
            "[--require-acceptance] [--require-crossover]\n"
            "\n"
            "Measures adaptive coarse graining against UppASD's production\n"
            "atomistic GPU Hamiltonian and integrator on one shared geometry.\n"
            "--require-acceptance fails on a parity violation; "
            "--require-crossover additionally fails when no crossover exists.");
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

void printSamples(const char* tag, const std::vector<double>& values) {
   for(std::size_t index = 0; index < values.size(); ++index)
      std::printf("%s index=%zu wall_us=%.6f\n", tag, index, values[index]);
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

   AdaptiveFixture(std::size_t blockCountValue, std::size_t atomsPerBlock,
                   Texture texture)
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
         double direction[3] = {};
         switch(texture) {
            case Texture::Spiral: {
               const double phase =
                  2.0 * 3.14159265358979323846 * static_cast<double>(block) /
                  static_cast<double>(blocks);
               direction[0] = std::cos(phase);
               direction[1] = std::sin(phase);
               direction[2] = 0.2;
               break;
            }
            case Texture::DomainWall: {
               // Two domains joined by a tanh wall a few blocks wide, so the
               // selector sees a genuinely localized refinement region rather
               // than uniformly distributed curvature.
               const double width = std::max(2.0, static_cast<double>(blocks) / 32.0);
               const double position =
                  (static_cast<double>(block) - 0.5 * static_cast<double>(blocks)) / width;
               direction[0] = 1.0 / std::cosh(position);
               direction[1] = 0.0;
               direction[2] = std::tanh(position);
               break;
            }
            case Texture::Uniform:
               direction[0] = 0.0;
               direction[1] = 0.0;
               direction[2] = 1.0;
               break;
         }
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
            bondMatrix[xyz + 3 * xyz + 9 * bond] = exchangeJ;
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
      // Unit moment and unit gyromagnetic scale: with mmom == 1 the adaptive
      // field b_i = sum_j (K s_j)/(mu m_i) and the production field
      // beff_i = sum_j ncoup_ij emomM_j are the same quantity, which is what
      // the parity check below verifies rather than assumes.
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

   // Short-range bonds whose field the atomistic kernel must actually evaluate:
   // a bond is live when either endpoint is atomistically owned.  Reported so a
   // reader can tell a reduction in active atoms from a reduction in real
   // Hamiltonian work (blueprint SS8.9: "the short-range Hamiltonian work must
   // also be reduced").
   std::size_t liveBonds(const std::vector<int>& activeAtoms,
                         const std::vector<int>& interfaceAtoms) const {
      std::vector<char> owned(atoms, 0);
      for(const int atom : activeAtoms)
         if(atom >= 0 && static_cast<std::size_t>(atom) < atoms) owned[atom] = 1;
      for(const int atom : interfaceAtoms)
         if(atom >= 0 && static_cast<std::size_t>(atom) < atoms) owned[atom] = 1;
      std::size_t live = 0;
      for(std::size_t bond = 0; bond < bonds; ++bond) {
         const std::size_t atomI = static_cast<std::size_t>(bondAtom[2 * bond] - 1);
         const std::size_t atomJ = static_cast<std::size_t>(bondAtom[2 * bond + 1] - 1);
         if(owned[atomI] || owned[atomJ]) ++live;
      }
      return live;
   }
};

// UppASD's real atomistic GPU path, constructed on the fixture geometry.
// Nothing here is a reimplementation: the neighbour lists are built in the
// layout GpuHamiltonianCalculations::initiate expects, and every field and
// integration call below is the same call gpuSDSimulation.cpp makes on a
// feature-off timestep.
struct ProductionAtomisticBaseline {
   Flag flags{};
   SimulationParameters simParam{};
   deviceHamiltonian hamiltonian;
   deviceLattice lattice;
   deviceEnergies energies;
   GpuHamiltonianCalculations hamCalc;
   GpuDepondtIntegrator integrator;
   unsigned int maxNeighbours = 0;
   double setupWallUs = 0.0;

   void build(const AdaptiveFixture& fixture, real timestep) {
      const auto setupBegin = std::chrono::steady_clock::now();
      const std::size_t atoms = fixture.atoms;
      const std::size_t ensembles = 1;

      // Bond list -> symmetric neighbour list.  The adaptive kernels store each
      // pair once and add the field to both endpoints; the production kernel
      // walks a per-atom list, so the same pair appears in both atoms' lists.
      std::vector<std::vector<unsigned int>> neighbours(atoms);
      for(std::size_t bond = 0; bond < fixture.bonds; ++bond) {
         const auto atomI = static_cast<unsigned int>(fixture.bondAtom[2 * bond] - 1);
         const auto atomJ = static_cast<unsigned int>(fixture.bondAtom[2 * bond + 1] - 1);
         neighbours[atomI].push_back(atomJ);
         neighbours[atomJ].push_back(atomI);
      }
      maxNeighbours = 0;
      for(const auto& list : neighbours)
         maxNeighbours = std::max(maxNeighbours,
                                  static_cast<unsigned int>(list.size()));
      if(maxNeighbours == 0)
         throw std::runtime_error("production baseline geometry has no bonds");

      const auto N = static_cast<long int>(atoms);
      const auto M = static_cast<long int>(ensembles);
      const auto mnn = static_cast<long int>(maxNeighbours);

      // nlist is indexed [atom + neighbour * N] and ncoup [site + neighbour *
      // NH]; SetupNeighbourList converts the Fortran 1-based entries below to
      // 0-based and zeroes the padding, exactly as in a production run.
      std::vector<unsigned int> hostNlist(atoms * maxNeighbours, 1u);
      std::vector<unsigned int> hostNlistSize(atoms, 0u);
      std::vector<unsigned int> hostAham(atoms, 0u);
      std::vector<real> hostNcoup(atoms * maxNeighbours, real(0));
      for(std::size_t atom = 0; atom < atoms; ++atom) {
         hostAham[atom] = static_cast<unsigned int>(atom + 1);
         hostNlistSize[atom] = static_cast<unsigned int>(neighbours[atom].size());
         for(std::size_t slot = 0; slot < neighbours[atom].size(); ++slot) {
            hostNlist[atom + slot * atoms] = neighbours[atom][slot] + 1u;
            hostNcoup[atom + slot * atoms] = static_cast<real>(exchangeJ);
         }
      }

      hamiltonian.aHam.Allocate(N);
      hamiltonian.nlist.Allocate(N, mnn);
      hamiltonian.nlistsize.Allocate(N);
      hamiltonian.ncoup.Allocate(N, mnn);
      hamiltonian.extfield.Allocate(3, N, M);
      hamiltonian.extfield.zeros();
      gpuCheck(GPU_MEMCPY(hamiltonian.aHam.data(), hostAham.data(),
                          hostAham.size() * sizeof(unsigned int),
                          GPU_MEMCPY_HOST_TO_DEVICE), "aHam upload");
      gpuCheck(GPU_MEMCPY(hamiltonian.nlist.data(), hostNlist.data(),
                          hostNlist.size() * sizeof(unsigned int),
                          GPU_MEMCPY_HOST_TO_DEVICE), "nlist upload");
      gpuCheck(GPU_MEMCPY(hamiltonian.nlistsize.data(), hostNlistSize.data(),
                          hostNlistSize.size() * sizeof(unsigned int),
                          GPU_MEMCPY_HOST_TO_DEVICE), "nlistsize upload");
      gpuCheck(GPU_MEMCPY(hamiltonian.ncoup.data(), hostNcoup.data(),
                          hostNcoup.size() * sizeof(real),
                          GPU_MEMCPY_HOST_TO_DEVICE), "ncoup upload");

      lattice.beff.Allocate(3, N, M);
      lattice.b2eff.Allocate(3, N, M);
      lattice.eneff.Allocate(3, N, M);
      lattice.emom.Allocate(3, N, M);
      lattice.emom2.Allocate(3, N, M);
      lattice.emomM.Allocate(3, N, M);
      lattice.mmom.Allocate(N, M);
      lattice.mmom0.Allocate(N, M);
      lattice.mmom2.Allocate(N, M);
      lattice.mmomi.Allocate(N, M);
      lattice.beff.zeros();
      lattice.b2eff.zeros();
      lattice.eneff.zeros();
      // Unit moments: emomM == emom, so the production and adaptive fields are
      // directly comparable without a moment-weighting conversion.
      std::vector<real> unit(atoms * ensembles, real(1));
      gpuCheck(GPU_MEMCPY(lattice.mmom.data(), unit.data(),
                          unit.size() * sizeof(real), GPU_MEMCPY_HOST_TO_DEVICE),
               "mmom upload");
      gpuCheck(GPU_MEMCPY(lattice.mmom0.data(), unit.data(),
                          unit.size() * sizeof(real), GPU_MEMCPY_HOST_TO_DEVICE),
               "mmom0 upload");
      gpuCheck(GPU_MEMCPY(lattice.mmom2.data(), unit.data(),
                          unit.size() * sizeof(real), GPU_MEMCPY_HOST_TO_DEVICE),
               "mmom2 upload");
      gpuCheck(GPU_MEMCPY(lattice.mmomi.data(), unit.data(),
                          unit.size() * sizeof(real), GPU_MEMCPY_HOST_TO_DEVICE),
               "mmomi upload");
      uploadState(fixture.atomDirection);
      energies.energyM.Allocate(M, static_cast<long int>(7));
      energies.energyM.zeros();

      // Per-atom temperature, matching hostLattice::temperature in a real run:
      // GpuThermfield::resetConstants copies it straight into an N-element
      // sigma factor, so an ensemble-sized array is the wrong shape.
      Tensor<real, 1> temperature;
      temperature.AllocateHost(N);
      for(long int atom = 0; atom < N; ++atom)
         temperature(atom) = real(0);
      lattice.temperature.Allocate(N);
      lattice.temperature.zeros();

      flags.do_dm = false;
      flags.do_jtensor = false;
      flags.do_aniso = 0;
      flags.do_ene = 0;
      flags.do_gpu_measurements = false;
      flags.do_gpu_correlations = false;

      simParam.stt = 'N';
      simParam.SDEalgh = 1;
      simParam.N = atoms;
      simParam.NH = atoms;
      simParam.M = ensembles;
      simParam.mnn = maxNeighbours;
      simParam.mnndm = 0;
      simParam.do_gpu_convolution = false;
      simParam.delta_t = timestep;
      simParam.gamma = real(1);
      simParam.k_bolt = real(1);
      simParam.mub = real(1);
      simParam.damping = real(0.1);
      simParam.Temp = real(0);
      simParam.rngType = defaultRngType();
      simParam.randomSeed = 1234ull;

      ParallelizationHelperInstance.initiate(
         static_cast<unsigned int>(atoms), static_cast<unsigned int>(ensembles),
         static_cast<unsigned int>(atoms));

      hamiltonian.neighbourListsPrepared = false;
      if(!hamCalc.initiate(flags, simParam, hamiltonian))
         throw std::runtime_error("production atomistic Hamiltonian failed to initiate");
      if(!integrator.initiate(simParam))
         throw std::runtime_error("production Depondt integrator failed to initiate");
      if(!integrator.initiateConstants(simParam, temperature))
         throw std::runtime_error("production Depondt constants failed to initiate");
      gpuCheck(GPU_DEVICE_SYNCHRONIZE(), "production baseline setup synchronization");
      setupWallUs = 1.0e6 * std::chrono::duration<double>(
         std::chrono::steady_clock::now() - setupBegin).count();
   }

   static decltype(SimulationParameters::rngType) defaultRngType() {
#if defined(HIP_V)
      return HIPRAND_RNG_PSEUDO_DEFAULT;
#else
      return CURAND_RNG_PSEUDO_DEFAULT;
#endif
   }

   void uploadState(const std::vector<real>& direction) {
      gpuCheck(GPU_MEMCPY(lattice.emom.data(), direction.data(),
                          direction.size() * sizeof(real),
                          GPU_MEMCPY_HOST_TO_DEVICE), "emom upload");
      gpuCheck(GPU_MEMCPY(lattice.emom2.data(), direction.data(),
                          direction.size() * sizeof(real),
                          GPU_MEMCPY_HOST_TO_DEVICE), "emom2 upload");
      gpuCheck(GPU_MEMCPY(lattice.emomM.data(), direction.data(),
                          direction.size() * sizeof(real),
                          GPU_MEMCPY_HOST_TO_DEVICE), "emomM upload");
   }

   void evaluateField() { hamCalc.heisge(lattice, energies, false); }

   // The complete feature-off timestep of gpuSDSimulation.cpp: measure field,
   // predictor, predictor field, corrector.  This -- not the field alone -- is
   // what adaptive coarse graining has to beat.
   void step() {
      hamCalc.heisge(lattice, energies, false);
      integrator.evolveFirst(lattice);
      hamCalc.heisge(lattice, energies, false);
      integrator.evolveSecond(lattice);
   }

   std::vector<real> downloadField() {
      std::vector<real> field(static_cast<std::size_t>(lattice.beff.size()));
      gpuCheck(GPU_MEMCPY(field.data(), lattice.beff.data(),
                          field.size() * sizeof(real),
                          GPU_MEMCPY_DEVICE_TO_HOST), "beff download");
      return field;
   }

   void release() {
      // GpuDepondtIntegrator::~GpuDepondtIntegrator already calls release(),
      // and GpuTensor::Free() leaves the pointer dangling, so releasing it here
      // as well would double-free.  Only the storage this struct allocated is
      // freed here.
      hamiltonian.aHam.Free();
      hamiltonian.nlist.Free();
      hamiltonian.nlistsize.Free();
      hamiltonian.ncoup.Free();
      hamiltonian.extfield.Free();
      lattice.beff.Free();
      lattice.b2eff.Free();
      lattice.eneff.Free();
      lattice.emom.Free();
      lattice.emom2.Free();
      lattice.emomM.Free();
      lattice.mmom.Free();
      lattice.mmom0.Free();
      lattice.mmom2.Free();
      lattice.mmomi.Free();
      lattice.temperature.Free();
      energies.energyM.Free();
   }
};

struct TimingSample {
   std::vector<double> wallUs;
   double medianUs = 0.0;
   double madUs = 0.0;

   void finish() {
      medianUs = median(wallUs);
      madUs = medianAbsoluteDeviation(wallUs);
   }
};

template<typename Body>
TimingSample timeSteadyState(const Options& options, Body body) {
   for(unsigned int iteration = 0; iteration < options.warmup; ++iteration)
      body();
   gpuCheck(GPU_DEVICE_SYNCHRONIZE(), "warm-up synchronization");
   TimingSample sample;
   for(unsigned int repetition = 0; repetition < options.repetitions; ++repetition) {
      const auto begin = std::chrono::steady_clock::now();
      for(unsigned int iteration = 0; iteration < options.iterations; ++iteration)
         body();
      gpuCheck(GPU_DEVICE_SYNCHRONIZE(), "steady-state synchronization");
      sample.wallUs.push_back(
         1.0e6 * std::chrono::duration<double>(
            std::chrono::steady_clock::now() - begin).count() / options.iterations);
   }
   sample.finish();
   return sample;
}

struct SweepPoint {
   double requestedFraction = 0.0;
   double activeDofRatio = 0.0;
   double activeAtomFraction = 0.0;
   double activeBlockFraction = 0.0;
   double interfaceFraction = 0.0;
   double liveBondFraction = 0.0;
   std::size_t activeAtoms = 0;
   std::size_t activeBlocks = 0;
   std::size_t interfaceAtoms = 0;
   std::size_t liveBonds = 0;
   std::size_t allocatedBytes = 0;
   // Headline timing, per-phase instrumentation disabled (RCG-08-FU2).
   TimingSample untimedStep;
   // Same work with per-phase device timers enabled: the breakdown, and the
   // cost of obtaining it.
   TimingSample instrumentedStep;
   GpuAdaptivePhaseMetrics phases{};
   double phaseIterations = 0.0;
};

} // namespace

int main(int argc, char** argv) {
   try {
      const Options options = parse(argc, argv);
      const bool fp64 = sizeof(real) == sizeof(double);
      const double parityTolerance = options.parityTolerance > 0.0 ?
         options.parityTolerance : (fp64 ? 1.0e-12 : 5.0e-5);
      std::printf(
         "adaptive-benchmark precision=%s backend=%s blocks=%zu "
         "atoms_per_block=%zu atoms=%zu ensembles=1 texture=%s warmup=%u "
         "iterations=%u repetitions=%u crossover_margin_percent=%.6f "
         "parity_tolerance=%.3e\n",
         fp64 ? "fp64" : "fp32",
#if defined(CUDA_V)
         "CUDA",
#else
         "HIP",
#endif
         options.blocks, options.atomsPerBlock,
         options.blocks * options.atomsPerBlock, textureName(options.texture),
         options.warmup, options.iterations, options.repetitions,
         100.0 * options.crossoverMargin, parityTolerance);

      AdaptiveFixture fixture(options.blocks, options.atomsPerBlock, options.texture);
      std::printf("adaptive-benchmark-geometry atoms=%zu blocks=%zu basis=%zu "
                  "bonds=%zu exchange_j=%.6f\n",
                  fixture.atoms, fixture.blocks, fixture.basis, fixture.bonds,
                  exchangeJ);

      const real timestep = real(1.0e-6);

      // ---- PERF-ATOMISTIC-PROD: UppASD's real atomistic GPU path ----------
      ProductionAtomisticBaseline baseline;
      baseline.build(fixture, timestep);
      std::printf("production-baseline-setup wall_us=%.6f max_neighbours=%u "
                  "sdealgh=1 integrator=depondt temperature=0 damping=%.6f\n",
                  baseline.setupWallUs, baseline.maxNeighbours,
                  static_cast<double>(baseline.simParam.damping));

      // ---- Adaptive runtime on the identical geometry ----------------------
      const auto topology = fixture.topology();
      const auto input = fixture.runtime();
      std::string diagnostic;
      if(!GpuAdaptiveRuntime::validate(topology, input, fixture.atoms, 1, diagnostic))
         throw std::runtime_error("benchmark fixture rejected: " + diagnostic);
      GpuAdaptiveRuntime runtime;
      const auto adaptiveSetupBegin = std::chrono::steady_clock::now();
      runtime.initialize(topology, input, fixture.atoms, 1);
      gpuCheck(GPU_DEVICE_SYNCHRONIZE(), "adaptive setup synchronization");
      const double adaptiveSetupUs = 1.0e6 * std::chrono::duration<double>(
         std::chrono::steady_clock::now() - adaptiveSetupBegin).count();
      std::printf("adaptive-setup wall_us=%.6f device_bytes=%zu\n",
                  adaptiveSetupUs, runtime.allocatedBytes());

      GpuTensor<real, 1> atomDirection, atomField, coarseField;
      atomDirection.Allocate(static_cast<index_t>(3 * fixture.atoms));
      atomField.Allocate(static_cast<index_t>(3 * fixture.atoms));
      coarseField.Allocate(static_cast<index_t>(3 * fixture.blocks));
      gpuCheck(GPU_MEMCPY(atomDirection.data(), fixture.atomDirection.data(),
                          fixture.atomDirection.size() * sizeof(real),
                          GPU_MEMCPY_HOST_TO_DEVICE),
               "benchmark direction upload");

      // ---- Identical-physics check ----------------------------------------
      // All blocks fine: the adaptive path then owns every atom atomistically
      // and must reproduce the production field on the same state.  Reported
      // before any timing so a failure can never be mistaken for a slow result.
      fixture.setFraction(1.0);
      runtime.updateBlockState(fixture.state.data(), fixture.blocks);
      (void)runtime.evaluateHybrid(atomDirection.data(), nullptr, nullptr,
                                   atomField.data(), coarseField.data());
      gpuCheck(GPU_DEVICE_SYNCHRONIZE(), "parity synchronization");
      std::vector<real> adaptiveFieldHost(3 * fixture.atoms);
      gpuCheck(GPU_MEMCPY(adaptiveFieldHost.data(), atomField.data(),
                          adaptiveFieldHost.size() * sizeof(real),
                          GPU_MEMCPY_DEVICE_TO_HOST), "adaptive field download");
      baseline.uploadState(fixture.atomDirection);
      baseline.evaluateField();
      gpuCheck(GPU_DEVICE_SYNCHRONIZE(), "production parity synchronization");
      const std::vector<real> productionFieldHost = baseline.downloadField();
      double maxAbsolute = 0.0, referenceScale = 0.0;
      for(std::size_t component = 0; component < adaptiveFieldHost.size(); ++component) {
         const double adaptiveValue = static_cast<double>(adaptiveFieldHost[component]);
         const double productionValue = static_cast<double>(productionFieldHost[component]);
         maxAbsolute = std::max(maxAbsolute, std::abs(adaptiveValue - productionValue));
         referenceScale = std::max(referenceScale, std::abs(productionValue));
      }
      // Normalized against the largest reference component rather than
      // component-wise: a spiral has components that pass through zero, and a
      // component-wise ratio there reports a large "relative error" for a
      // difference at the last bit of the field's actual scale.
      const double maxRelative = referenceScale > 0.0 ?
         maxAbsolute / referenceScale : maxAbsolute;
      const bool parityAccepted = maxRelative <= parityTolerance;
      std::printf(
         "atomistic-parity max_abs_diff=%.6e max_rel_diff=%.6e "
         "reference_scale=%.6e tolerance=%.3e components=%zu result=%s\n",
         maxAbsolute, maxRelative, referenceScale, parityTolerance,
         adaptiveFieldHost.size(), parityAccepted ? "PASS" : "FAIL");
      if(!parityAccepted)
         std::puts("atomistic-parity note=baseline and adaptive do not evaluate "
                   "the same physics on this geometry; no speedup is reported");

      // ---- PERF-ATOMISTIC-PROD steady state -------------------------------
      TimingSample productionField, productionStep;
      TimingSample adaptiveAllFineStep;
      std::vector<SweepPoint> sweep;
      if(parityAccepted) {
         productionField = timeSteadyState(options, [&]() { baseline.evaluateField(); });
         productionStep = timeSteadyState(options, [&]() { baseline.step(); });
         std::printf(
            "production-atomistic field_wall_us=%.6f field_mad_us=%.6f "
            "step_wall_us=%.6f step_mad_us=%.6f atoms=%zu bonds=%zu "
            "repetitions=%u iterations=%u\n",
            productionField.medianUs, productionField.madUs,
            productionStep.medianUs, productionStep.madUs,
            fixture.atoms, fixture.bonds, options.repetitions, options.iterations);
         // Disclosure, not a correction.  GpuDepondtIntegrator::evolveFirst
         // calls thermfield.randomize() unconditionally, so the production step
         // generates 3*N*M random numbers every step even at temp 0.  The
         // adaptive Heun path is deterministic and does no RNG.  That is a real
         // difference between the two production paths and part of what a user
         // actually pays, so the accepted crossover below is computed against
         // the complete step -- but a speedup that came only from skipping the
         // thermal field would not be a coarse-graining result, so the
         // Hamiltonian-only floor (two heisge calls, the work adaptive really
         // replaces) is reported alongside it as the conservative bound.
         std::printf(
            "production-atomistic-composition thermfield_rng_per_step=%zu "
            "heisge_calls_per_step=2 hamiltonian_only_floor_us=%.6f "
            "rng_and_integration_residual_us=%.6f note=%s\n",
            3 * fixture.atoms, 2.0 * productionField.medianUs,
            productionStep.medianUs - 2.0 * productionField.medianUs,
            "adaptive Heun is deterministic; production Depondt randomizes even at T=0");
         printSamples("production-atomistic-step-sample", productionStep.wallUs);
         printSamples("production-atomistic-field-sample", productionField.wallUs);

         // ---- PERF-CG-SWEEP ------------------------------------------------
         const std::vector<double> fractions = {
            1.0, 0.75, 0.5, 0.25, 0.125, 0.0625, 0.0
         };
         for(const double fraction : fractions) {
            SweepPoint point;
            point.requestedFraction = fraction;
            fixture.setFraction(fraction);
            runtime.updateBlockState(fixture.state.data(), fixture.blocks);
            const auto snapshot = runtime.downloadWorkSnapshot();
            point.activeAtoms = snapshot.activeAtoms.size();
            point.activeBlocks = snapshot.activeBlocks.size();
            point.interfaceAtoms = snapshot.interfaceAtoms.size();
            point.liveBonds = fixture.liveBonds(snapshot.activeAtoms,
                                                snapshot.interfaceAtoms);
            point.activeAtomFraction = static_cast<double>(point.activeAtoms) /
               static_cast<double>(fixture.atoms);
            point.activeBlockFraction = static_cast<double>(point.activeBlocks) /
               static_cast<double>(fixture.blocks);
            point.interfaceFraction = static_cast<double>(point.interfaceAtoms) /
               static_cast<double>(fixture.atoms);
            point.liveBondFraction = static_cast<double>(point.liveBonds) /
               static_cast<double>(fixture.bonds);
            point.activeDofRatio =
               static_cast<double>(point.activeAtoms + point.activeBlocks) /
               static_cast<double>(fixture.atoms);
            point.allocatedBytes = runtime.allocatedBytes();

            // Headline: complete Heun step, phase instrumentation off.
            runtime.setPhaseTimingEnabled(false);
            point.untimedStep = timeSteadyState(options, [&]() {
               runtime.integrateHeun(timestep, atomDirection.data());
            });

            // Breakdown: identical work, per-phase device timers on.
            runtime.setPhaseTimingEnabled(true);
            runtime.resetPhaseMetrics();
            point.instrumentedStep = timeSteadyState(options, [&]() {
               runtime.integrateHeun(timestep, atomDirection.data());
            });
            point.phases = runtime.phaseMetrics();
            point.phaseIterations = static_cast<double>(options.iterations) *
               static_cast<double>(options.repetitions);
            sweep.push_back(point);

            const double scale = 1000.0 / point.phaseIterations;
            std::printf(
               "adaptive-sweep requested_fraction=%.6f active_atoms=%zu "
               "active_blocks=%zu interface_atoms=%zu live_bonds=%zu "
               "active_atom_fraction=%.6f active_block_fraction=%.6f "
               "interface_fraction=%.6f live_bond_fraction=%.6f "
               "active_dof_ratio=%.6f step_wall_us=%.6f step_mad_us=%.6f "
               "instrumented_wall_us=%.6f instrumented_mad_us=%.6f "
               "instrumentation_overhead_percent=%.6f device_bytes=%zu\n",
               point.requestedFraction, point.activeAtoms, point.activeBlocks,
               point.interfaceAtoms, point.liveBonds, point.activeAtomFraction,
               point.activeBlockFraction, point.interfaceFraction,
               point.liveBondFraction, point.activeDofRatio,
               point.untimedStep.medianUs, point.untimedStep.madUs,
               point.instrumentedStep.medianUs, point.instrumentedStep.madUs,
               100.0 * (point.instrumentedStep.medianUs - point.untimedStep.medianUs) /
                  point.untimedStep.medianUs,
               point.allocatedBytes);
            std::printf(
               "adaptive-sweep-phases requested_fraction=%.6f atomistic_us=%.6f "
               "coarse_us=%.6f interface_us=%.6f selector_us=%.6f "
               "polarization_us=%.6f compaction_us=%.6f fft_us=%.6f "
               "integration_us=%.6f phase_sum_us=%.6f unaccounted_us=%.6f "
               "phase_syncs_per_step=%.3f\n",
               point.requestedFraction,
               scale * point.phases.atomisticMilliseconds,
               scale * point.phases.coarseMilliseconds,
               scale * point.phases.interfaceMilliseconds,
               scale * point.phases.selectorMilliseconds,
               scale * point.phases.polarizationMilliseconds,
               scale * point.phases.compactionMilliseconds,
               scale * point.phases.fftMilliseconds,
               scale * point.phases.integrationMilliseconds,
               scale * (point.phases.atomisticMilliseconds +
                        point.phases.coarseMilliseconds +
                        point.phases.interfaceMilliseconds +
                        point.phases.selectorMilliseconds +
                        point.phases.polarizationMilliseconds +
                        point.phases.compactionMilliseconds +
                        point.phases.fftMilliseconds +
                        point.phases.integrationMilliseconds),
               point.instrumentedStep.medianUs -
                  scale * (point.phases.atomisticMilliseconds +
                           point.phases.coarseMilliseconds +
                           point.phases.interfaceMilliseconds +
                           point.phases.selectorMilliseconds +
                           point.phases.polarizationMilliseconds +
                           point.phases.compactionMilliseconds +
                           point.phases.fftMilliseconds +
                           point.phases.integrationMilliseconds),
               static_cast<double>(point.phases.phaseSynchronizations) /
                  point.phaseIterations);
            std::printf(
               "adaptive-sweep-launches requested_fraction=%.6f atomistic=%.3f "
               "coarse=%.3f interface=%.3f selector=%.3f polarization=%.3f "
               "compaction=%.3f integration=%.3f\n",
               point.requestedFraction,
               static_cast<double>(point.phases.atomisticLaunches) / point.phaseIterations,
               static_cast<double>(point.phases.coarseLaunches) / point.phaseIterations,
               static_cast<double>(point.phases.interfaceLaunches) / point.phaseIterations,
               static_cast<double>(point.phases.selectorLaunches) / point.phaseIterations,
               static_cast<double>(point.phases.polarizationLaunches) / point.phaseIterations,
               static_cast<double>(point.phases.compactionLaunches) / point.phaseIterations,
               static_cast<double>(point.phases.integrationLaunches) / point.phaseIterations);
            char tag[96];
            std::snprintf(tag, sizeof(tag),
                          "adaptive-sweep-sample fraction=%.6f", fraction);
            printSamples(tag, point.untimedStep.wallUs);
         }
         adaptiveAllFineStep = sweep.front().untimedStep;
      }

      // ---- Selector and compaction cost at a mixed point -------------------
      double selectorWallUs = 0.0, selectorDeviceUs = 0.0;
      double compactionWallUs = 0.0, compactionDeviceUs = 0.0, compactionHostWaitUs = 0.0;
      if(parityAccepted) {
         fixture.setFraction(0.5);
         runtime.setPhaseTimingEnabled(true);
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
         selectorWallUs = 1.0e6 * std::chrono::duration<double>(
            std::chrono::steady_clock::now() - selectorBegin).count() / options.iterations;
         selectorDeviceUs = 1000.0 * runtime.phaseMetrics().selectorMilliseconds /
            options.iterations;

         const auto compactionBefore = runtime.compactionMetrics();
         const auto compactionWallBegin = std::chrono::steady_clock::now();
         for(unsigned int iteration = 0; iteration < options.iterations; ++iteration) {
            fixture.setFraction(iteration % 2 ? 0.5 : 0.25);
            runtime.updateBlockState(fixture.state.data(), fixture.blocks);
         }
         compactionWallUs = 1.0e6 * std::chrono::duration<double>(
            std::chrono::steady_clock::now() - compactionWallBegin).count() /
            options.iterations;
         const auto compactionAfter = runtime.compactionMetrics();
         compactionDeviceUs = 1000.0 *
            (compactionAfter.deviceMilliseconds - compactionBefore.deviceMilliseconds) /
            options.iterations;
         compactionHostWaitUs = 1000.0 *
            (compactionAfter.hostWaitMilliseconds - compactionBefore.hostWaitMilliseconds) /
            options.iterations;
         std::printf(
            "adaptive-overhead selector_device_us=%.6f selector_wall_us=%.6f "
            "compaction_device_us=%.6f compaction_host_wait_us=%.6f "
            "compaction_wall_us=%.6f block_bytes_per_update=%zu device_bytes=%zu\n",
            selectorDeviceUs, selectorWallUs, compactionDeviceUs,
            compactionHostWaitUs, compactionWallUs,
            fixture.blocks * sizeof(int), runtime.allocatedBytes());
      }

      // ---- Crossover against the production baseline -----------------------
      const SweepPoint* crossover = nullptr;
      if(parityAccepted) {
         for(const SweepPoint& point : sweep) {
            const double uncertainty =
               3.0 * (productionStep.madUs + point.untimedStep.madUs);
            const double target = productionStep.medianUs * (1.0 - options.crossoverMargin);
            if(point.untimedStep.medianUs + uncertainty < target) {
               crossover = &point;
               break;
            }
         }
         if(crossover) {
            std::printf(
               "production-crossover result=PASS requested_fraction=%.6f "
               "active_dof_ratio=%.6f active_atom_fraction=%.6f "
               "live_bond_fraction=%.6f adaptive_step_wall_us=%.6f "
               "production_step_wall_us=%.6f speedup=%.6f "
               "uncertainty_us=%.6f acceptance_margin_percent=%.6f\n",
               crossover->requestedFraction, crossover->activeDofRatio,
               crossover->activeAtomFraction, crossover->liveBondFraction,
               crossover->untimedStep.medianUs, productionStep.medianUs,
               productionStep.medianUs / crossover->untimedStep.medianUs,
               3.0 * (productionStep.madUs + crossover->untimedStep.madUs),
               100.0 * options.crossoverMargin);
            std::printf(
               "production-crossover-conservative hamiltonian_only_floor_us=%.6f "
               "adaptive_step_wall_us=%.6f conservative_speedup=%.6f note=%s\n",
               2.0 * productionField.medianUs, crossover->untimedStep.medianUs,
               2.0 * productionField.medianUs / crossover->untimedStep.medianUs,
               "excludes the production thermfield RNG the adaptive path never runs");
         } else {
            std::printf(
               "production-crossover result=NOT_OBSERVED "
               "production_step_wall_us=%.6f best_adaptive_step_wall_us=%.6f "
               "acceptance_margin_percent=%.6f\n",
               productionStep.medianUs,
               sweep.empty() ? 0.0 : sweep.back().untimedStep.medianUs,
               100.0 * options.crossoverMargin);
         }
         // Diagnostic only.  This is the number F-10 objected to: adaptive
         // against itself.  It is never the accepted crossover.
         if(!sweep.empty())
            std::printf(
               "adaptive-self-reference note=DIAGNOSTIC_NOT_A_BASELINE "
               "all_fine_step_wall_us=%.6f zero_fine_step_wall_us=%.6f "
               "self_ratio=%.6f\n",
               adaptiveAllFineStep.medianUs, sweep.back().untimedStep.medianUs,
               adaptiveAllFineStep.medianUs / sweep.back().untimedStep.medianUs);
      }

      coarseField.Free();
      atomField.Free();
      atomDirection.Free();
      runtime.release();
      baseline.release();

      if(options.requireAcceptance && !parityAccepted) return 2;
      if(options.requireCrossover && !crossover) return 3;
      return 0;
   } catch(const std::exception& error) {
      std::fprintf(stderr, "FAIL GPU adaptive runtime benchmark: %s\n",
                   error.what());
      return 1;
   }
}
