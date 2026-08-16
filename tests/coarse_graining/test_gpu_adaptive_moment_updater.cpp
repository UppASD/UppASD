// CGP-03B: GpuAdaptiveMomentUpdater must (1) reproduce GpuMomentUpdater
// bitwise when asked to touch every atom (updateFull(), and
// updateActiveOnly() given the complete atom list -- the "all-fine control"
// the task text asks for), and (2) leave every atom outside the active list
// exactly as it was before the call when given a genuine subset -- proving
// the touched-only scoping is real, not merely untested. See
// gpuAdaptiveMomentUpdater.hpp and
// docs/CGP-03B_ADAPTIVE_MOMENT_UPDATER_EVIDENCE.md.

#include "gpuMomentUpdater.hpp"
#include "gpuAdaptiveMomentUpdater.hpp"
#include "gpuParallelizationHelper.hpp"
#include "gpuStructures.hpp"
#include "tensor.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

constexpr std::size_t atoms = 6;
constexpr std::size_t ensembles = 2;
constexpr std::size_t vectorCount = 3 * atoms * ensembles;
constexpr std::size_t scalarCount = atoms * ensembles;

bool backendAvailable() {
   int count = 0;
#if defined(CUDA_V)
   const auto status = cudaGetDeviceCount(&count);
   return status == cudaSuccess && count > 0 && cudaFree(nullptr) == cudaSuccess;
#elif defined(HIP_V)
   const auto status = hipGetDeviceCount(&count);
   return status == hipSuccess && count > 0 && hipFree(nullptr) == hipSuccess;
#endif
}

void check(GPU_ERROR_T status, const char* context) {
   if(status != GPU_SUCCESS)
      throw std::runtime_error(std::string(context) + ": " + GPU_GET_ERROR_STRING(status));
}

struct State {
   deviceLattice lattice;

   State() {
      lattice.mmom.Allocate(atoms, ensembles);
      lattice.mmom0.Allocate(atoms, ensembles);
      lattice.mmom2.Allocate(atoms, ensembles);
      lattice.mmomi.Allocate(atoms, ensembles);
      lattice.emom.Allocate(3, atoms, ensembles);
      lattice.emom2.Allocate(3, atoms, ensembles);
      lattice.emomM.Allocate(3, atoms, ensembles);
   }

   ~State() {
      lattice.mmom.Free();
      lattice.mmom0.Free();
      lattice.mmom2.Free();
      lattice.mmomi.Free();
      lattice.emom.Free();
      lattice.emom2.Free();
      lattice.emomM.Free();
   }

   State(const State&) = delete;
   State& operator=(const State&) = delete;
};

// A fixed, physically-plausible direction field for emom2 (this step's
// freshly-integrated/committed corrector result -- the buffer both updaters
// read mz from under mompar 1/2).
std::vector<real> directions() {
   std::vector<real> result(vectorCount);
   for(std::size_t ensemble = 0; ensemble < ensembles; ++ensemble) {
      for(std::size_t atom = 0; atom < atoms; ++atom) {
         const real x = real(0.11) * real(atom + 1) - real(0.02) * real(ensemble);
         const real y = real(0.05) * real(ensemble + 1);
         const real z = real(1) - real(0.04) * real(atom) + real(0.01) * real(ensemble);
         const real norm = std::sqrt(x * x + y * y + z * z);
         const std::size_t element = 3 * (atom + atoms * ensemble);
         result[element] = x / norm;
         result[element + 1] = y / norm;
         result[element + 2] = z / norm;
      }
   }
   return result;
}

std::vector<real> baseMoments() {
   std::vector<real> result(scalarCount);
   for(std::size_t i = 0; i < scalarCount; ++i) result[i] = real(1.6) + real(0.1) * real(i);
   return result;
}

// Sentinel values distinguishable from anything either updater's formulas
// would ever produce from the inputs above, so "untouched" is unambiguous.
constexpr real sentinelScalar = real(-777.0);
constexpr real sentinelVector = real(-999.0);

void uploadScalar(GpuTensor<real, 2>& tensor, const std::vector<real>& host) {
   check(GPU_MEMCPY(tensor.data(), host.data(), host.size() * sizeof(real),
                    GPU_MEMCPY_HOST_TO_DEVICE), "scalar upload");
}

void uploadVector(GpuTensor<real, 3>& tensor, const std::vector<real>& host) {
   check(GPU_MEMCPY(tensor.data(), host.data(), host.size() * sizeof(real),
                    GPU_MEMCPY_HOST_TO_DEVICE), "vector upload");
}

std::vector<real> downloadScalar(const GpuTensor<real, 2>& tensor) {
   std::vector<real> result(scalarCount);
   check(GPU_MEMCPY(result.data(), tensor.data(), result.size() * sizeof(real),
                    GPU_MEMCPY_DEVICE_TO_HOST), "scalar download");
   return result;
}

std::vector<real> downloadVector(const GpuTensor<real, 3>& tensor) {
   std::vector<real> result(vectorCount);
   check(GPU_MEMCPY(result.data(), tensor.data(), result.size() * sizeof(real),
                    GPU_MEMCPY_DEVICE_TO_HOST), "vector download");
   return result;
}

// Common, deliberately-not-yet-updated starting state shared by every
// fixture below: emom2 holds this step's real direction; mmom holds a
// distinct "previous step" value; mmom2/mmomi/emomM hold sentinels so a
// skipped write is unambiguous; emom/mmom0 are the remaining required
// inputs.
void seed(State& state, const std::vector<real>& direction, const std::vector<real>& base) {
   const std::vector<real> previousMmom(scalarCount, real(1.23));
   const std::vector<real> sentinelMmom2(scalarCount, sentinelScalar);
   const std::vector<real> sentinelMmomi(scalarCount, sentinelScalar);
   const std::vector<real> sentinelEmomM(vectorCount, sentinelVector);
   const std::vector<real> staleEmom(vectorCount, real(0.0));

   uploadVector(state.lattice.emom, staleEmom);
   uploadVector(state.lattice.emom2, direction);
   uploadScalar(state.lattice.mmom, previousMmom);
   uploadScalar(state.lattice.mmom0, base);
   uploadScalar(state.lattice.mmom2, sentinelMmom2);
   uploadScalar(state.lattice.mmomi, sentinelMmomi);
   uploadVector(state.lattice.emomM, sentinelEmomM);
}

void requireBitwiseEqual(const std::vector<real>& actual, const std::vector<real>& expected,
                         const char* label) {
   if(actual.size() != expected.size() ||
      std::memcmp(actual.data(), expected.data(), actual.size() * sizeof(real)) != 0)
      throw std::runtime_error(std::string("mismatch in ") + label);
}

struct DeviceList {
   GpuTensor<int, 1> values;

   explicit DeviceList(const std::vector<int>& host) {
      values.Allocate(static_cast<index_t>(host.size()));
      Tensor<int, 1> view(const_cast<int*>(host.data()), static_cast<index_t>(host.size()));
      values.copy_sync(view);
   }

   ~DeviceList() { values.Free(); }

   DeviceList(const DeviceList&) = delete;
   DeviceList& operator=(const DeviceList&) = delete;
};

// Part 1: updateFull() must be byte-for-byte identical to
// GpuMomentUpdater::update() -- same buffers, same formulas, only the
// object issuing the launches differs.
void runFullParity(int mompar, char initexc, const char* label) {
   const auto direction = directions();
   const auto base = baseMoments();

   State reference;
   State adaptive;
   seed(reference, direction, base);
   seed(adaptive, direction, base);

   GpuMomentUpdater referenceUpdater(reference.lattice, mompar, initexc);
   GpuAdaptiveMomentUpdater adaptiveUpdater(adaptive.lattice, mompar, initexc);
   referenceUpdater.update();
   adaptiveUpdater.updateFull();
   check(GPU_STREAM_SYNC(ParallelizationHelperInstance.getWorkStream()), "full-parity sync");

   requireBitwiseEqual(downloadScalar(adaptive.lattice.mmom), downloadScalar(reference.lattice.mmom),
                       (std::string(label) + " mmom").c_str());
   requireBitwiseEqual(downloadScalar(adaptive.lattice.mmomi), downloadScalar(reference.lattice.mmomi),
                       (std::string(label) + " mmomi").c_str());
   requireBitwiseEqual(downloadVector(adaptive.lattice.emom), downloadVector(reference.lattice.emom),
                       (std::string(label) + " emom").c_str());
   requireBitwiseEqual(downloadVector(adaptive.lattice.emomM), downloadVector(reference.lattice.emomM),
                       (std::string(label) + " emomM").c_str());
}

// Part 2 ("all-fine control"): updateActiveOnly() given every atom must
// also match GpuMomentUpdater::update() bitwise -- the active-list launch
// geometry must not itself change any answer.
void runActiveAllParity(int mompar, char initexc, const char* label) {
   const auto direction = directions();
   const auto base = baseMoments();
   std::vector<int> all;
   for(std::size_t i = 1; i <= atoms; ++i) all.push_back(static_cast<int>(i));
   DeviceList allList(all);

   State reference;
   State adaptive;
   seed(reference, direction, base);
   seed(adaptive, direction, base);

   GpuMomentUpdater referenceUpdater(reference.lattice, mompar, initexc);
   GpuAdaptiveMomentUpdater adaptiveUpdater(adaptive.lattice, mompar, initexc);
   referenceUpdater.update();
   // active_atom_kernel launches activeCount*M threads and maps each
   // (slot, ensemble) pair to the ordinary flattened atom index itself, so
   // supplying every physical atom id here covers every ensemble too.
   adaptiveUpdater.updateActiveOnly(allList.values.data(), all.size());
   check(GPU_STREAM_SYNC(ParallelizationHelperInstance.getWorkStream()), "active-all-parity sync");

   requireBitwiseEqual(downloadScalar(adaptive.lattice.mmom), downloadScalar(reference.lattice.mmom),
                       (std::string(label) + " mmom").c_str());
   requireBitwiseEqual(downloadScalar(adaptive.lattice.mmomi), downloadScalar(reference.lattice.mmomi),
                       (std::string(label) + " mmomi").c_str());
   requireBitwiseEqual(downloadVector(adaptive.lattice.emom), downloadVector(reference.lattice.emom),
                       (std::string(label) + " emom").c_str());
   requireBitwiseEqual(downloadVector(adaptive.lattice.emomM), downloadVector(reference.lattice.emomM),
                       (std::string(label) + " emomM").c_str());
}

// Part 3 (negative control): a genuine subset must leave every atom outside
// it exactly as it was before the call (mmom/mmomi/emomM pinned at the
// sentinel seeded above), while atoms inside it match the reference full
// update, and every atom's emom (the whole-tensor swap) tracks the
// reference regardless of active-list membership -- proving CGP-03's own
// atom-direction contract is untouched by this task.
void runActiveSubsetNegativeControl(int mompar, char initexc, const char* label) {
   const auto direction = directions();
   const auto base = baseMoments();
   const std::vector<int> subset = {2, 4, 6};
   DeviceList subsetList(subset);

   State reference;
   State adaptive;
   seed(reference, direction, base);
   seed(adaptive, direction, base);

   GpuMomentUpdater referenceUpdater(reference.lattice, mompar, initexc);
   GpuAdaptiveMomentUpdater adaptiveUpdater(adaptive.lattice, mompar, initexc);
   referenceUpdater.update();
   adaptiveUpdater.updateActiveOnly(subsetList.values.data(), subset.size());
   check(GPU_STREAM_SYNC(ParallelizationHelperInstance.getWorkStream()), "active-subset sync");

   const auto refMmom = downloadScalar(reference.lattice.mmom);
   const auto refMmomi = downloadScalar(reference.lattice.mmomi);
   const auto refEmomM = downloadVector(reference.lattice.emomM);
   const auto refEmom = downloadVector(reference.lattice.emom);
   const auto gotMmom = downloadScalar(adaptive.lattice.mmom);
   const auto gotMmomi = downloadScalar(adaptive.lattice.mmomi);
   const auto gotEmomM = downloadVector(adaptive.lattice.emomM);
   const auto gotEmom = downloadVector(adaptive.lattice.emom);

   auto isSelected = [&](std::size_t atom) {
      for(int s : subset)
         if(static_cast<std::size_t>(s) == atom + 1) return true;
      return false;
   };

   for(std::size_t ensemble = 0; ensemble < ensembles; ++ensemble) {
      for(std::size_t atom = 0; atom < atoms; ++atom) {
         const std::size_t scalarIndex = atom + atoms * ensemble;
         const std::size_t vecBase = 3 * scalarIndex;
         if(isSelected(atom)) {
            if(gotMmom[scalarIndex] != refMmom[scalarIndex] ||
               gotMmomi[scalarIndex] != refMmomi[scalarIndex] ||
               std::memcmp(&gotEmomM[vecBase], &refEmomM[vecBase], 3 * sizeof(real)) != 0)
               throw std::runtime_error(std::string(label) + ": active atom diverged from reference");
         } else {
            // Not in the active list: mmomi must be pinned exactly at the
            // sentinel seeded pre-call (proving the Copy1/Copy2 kernel was
            // skipped, not merely coincidentally correct), while emom still
            // tracks the reference exactly (proving the whole-tensor swap
            // stayed global, matching the CGP-03 atom-direction contract
            // this task does not change).
            //
            // mmom is a special case for mompar==0 only: that branch is a
            // double whole-tensor swap (mmom2.swap(mmom) then, later,
            // mmom.swap(mmom2)) with no per-atom kernel at all -- see
            // GpuAdaptiveMomentUpdater::updateActiveOnly's case-0 comment --
            // so mmom nets out unchanged (and therefore still equal to the
            // reference) for every atom regardless of active-list
            // membership. Only mompar 1/2 actually replace that first swap
            // with a real per-atom write to mmom2, which is what makes mmom
            // pinnable at the sentinel for an untouched atom.
            if(mompar == 0) {
               if(gotMmom[scalarIndex] != refMmom[scalarIndex])
                  throw std::runtime_error(std::string(label) +
                                           ": mompar=0 inactive atom's mmom diverged from the "
                                           "double-swap identity reference");
            } else if(gotMmom[scalarIndex] != sentinelScalar) {
               throw std::runtime_error(std::string(label) + ": inactive atom's mmom was not pinned");
            }
            if(gotMmomi[scalarIndex] != sentinelScalar)
               throw std::runtime_error(std::string(label) + ": inactive atom's mmomi was not pinned");
            for(int c = 0; c < 3; ++c)
               if(gotEmomM[vecBase + c] != sentinelVector)
                  throw std::runtime_error(std::string(label) + ": inactive atom's emomM was not pinned");
            if(std::memcmp(&gotEmom[vecBase], &refEmom[vecBase], 3 * sizeof(real)) != 0)
               throw std::runtime_error(std::string(label) + ": inactive atom's emom diverged from the "
                                                              "whole-tensor swap reference");
         }
      }
   }
}

// Not a ctest (hardware-dependent numbers, no pass/fail): a standalone
// --benchmark mode, following the project's accepted ABBA protocol (one
// discarded warm-up, then alternating arms, median+MAD reported) --
// see docs/CGP-03B_ADAPTIVE_MOMENT_UPDATER_EVIDENCE.md for the readings
// this produced and project memory ("GPU benchmark cold start",
// "Shared GPU host contention") for why the protocol looks like this.
double median(std::vector<double> samples) {
   std::sort(samples.begin(), samples.end());
   const std::size_t n = samples.size();
   return (n % 2 == 0) ? (samples[n / 2 - 1] + samples[n / 2]) / 2.0 : samples[n / 2];
}

double mad(std::vector<double> samples, double centre) {
   for(double& s : samples) s = std::fabs(s - centre);
   return median(samples);
}

void runBenchmark(std::size_t benchAtoms, std::size_t benchEnsembles, double activeFraction,
                  int mompar) {
   ParallelizationHelperInstance.free();
   ParallelizationHelperInstance.initiate(static_cast<unsigned int>(benchAtoms),
                                          static_cast<unsigned int>(benchEnsembles), 0);

   deviceLattice lattice;
   lattice.mmom.Allocate(benchAtoms, benchEnsembles);
   lattice.mmom0.Allocate(benchAtoms, benchEnsembles);
   lattice.mmom2.Allocate(benchAtoms, benchEnsembles);
   lattice.mmomi.Allocate(benchAtoms, benchEnsembles);
   lattice.emom.Allocate(3, benchAtoms, benchEnsembles);
   lattice.emom2.Allocate(3, benchAtoms, benchEnsembles);
   lattice.emomM.Allocate(3, benchAtoms, benchEnsembles);

   const std::size_t scalarN = benchAtoms * benchEnsembles;
   const std::size_t vectorN = 3 * scalarN;
   std::vector<real> ones(scalarN, real(1.2));
   std::vector<real> dir(vectorN);
   for(std::size_t i = 0; i < scalarN; ++i) {
      dir[3 * i] = real(0.1);
      dir[3 * i + 1] = real(0.2);
      dir[3 * i + 2] = real(0.97);
   }
   uploadScalar(lattice.mmom, ones);
   uploadScalar(lattice.mmom0, ones);
   uploadScalar(lattice.mmom2, ones);
   uploadScalar(lattice.mmomi, ones);
   uploadVector(lattice.emom, dir);
   uploadVector(lattice.emom2, dir);
   uploadVector(lattice.emomM, dir);

   const std::size_t activeCount =
      std::max<std::size_t>(1, static_cast<std::size_t>(activeFraction * double(benchAtoms)));
   std::vector<int> activeHost(activeCount);
   for(std::size_t i = 0; i < activeCount; ++i) activeHost[i] = static_cast<int>(i + 1);
   DeviceList activeList(activeHost);

   GpuAdaptiveMomentUpdater updater(lattice, mompar, 'N');

   auto timeCall = [&](bool full) {
      const auto started = std::chrono::steady_clock::now();
      if(full)
         updater.updateFull();
      else
         updater.updateActiveOnly(activeList.values.data(), activeCount);
      check(GPU_STREAM_SYNC(ParallelizationHelperInstance.getWorkStream()), "benchmark sync");
      const auto stopped = std::chrono::steady_clock::now();
      return std::chrono::duration<double, std::micro>(stopped - started).count();
   };

   // Discarded warm-up (project memory: first call at a given size is
   // systematically slow from cold SM clocks).
   timeCall(true);

   std::vector<double> fullSamples, activeSamples;
   for(int i = 0; i < 8; ++i) {
      // ABBA order per sample pair.
      fullSamples.push_back(timeCall(true));
      activeSamples.push_back(timeCall(false));
      activeSamples.push_back(timeCall(false));
      fullSamples.push_back(timeCall(true));
   }

   const double fullMedian = median(fullSamples);
   const double activeMedian = median(activeSamples);
   std::printf("benchmark atoms=%zu ensembles=%zu active_fraction=%.4f mompar=%d "
               "full_us_median=%.3f full_us_mad=%.3f active_us_median=%.3f active_us_mad=%.3f "
               "speedup=%.3f\n",
               benchAtoms, benchEnsembles, activeFraction, mompar, fullMedian,
               mad(fullSamples, fullMedian), activeMedian, mad(activeSamples, activeMedian),
               fullMedian / activeMedian);

   lattice.mmom.Free();
   lattice.mmom0.Free();
   lattice.mmom2.Free();
   lattice.mmomi.Free();
   lattice.emom.Free();
   lattice.emom2.Free();
   lattice.emomM.Free();
   ParallelizationHelperInstance.free();
}

} // namespace

int main(int argc, char** argv) {
   if(!backendAvailable()) {
      std::puts("GPU-ADAPTIVE-MOMENT-UPDATER unavailable: no backend device");
      return 77;
   }
   if(argc > 1 && std::strcmp(argv[1], "--benchmark") == 0) {
      try {
         // 65536 atoms matches CGP-07's historical reference size (and the
         // one CGP-03's own materialization-overhead benchmark used), so
         // this is directly comparable to that evidence.
         for(double fraction : {0.0125, 0.0625, 0.25, 1.0})
            for(int mompar : {0, 1, 2})
               runBenchmark(65536, 1, fraction, mompar);
         return 0;
      } catch(const std::exception& error) {
         std::fprintf(stderr, "%s\n", error.what());
         return 1;
      }
   }
   try {
      ParallelizationHelperInstance.initiate(atoms, ensembles, 0);

      for(int mompar : {0, 1, 2}) {
         for(char initexc : {'N', 'I'}) {
            const std::string label =
               "mompar=" + std::to_string(mompar) + " initexc=" + std::string(1, initexc);
            runFullParity(mompar, initexc, label.c_str());
            runActiveAllParity(mompar, initexc, label.c_str());
            runActiveSubsetNegativeControl(mompar, initexc, label.c_str());
         }
      }

      ParallelizationHelperInstance.free();
      std::puts("GPU-ADAPTIVE-MOMENT-UPDATER passed");
      return 0;
   } catch(const std::exception& error) {
      std::fprintf(stderr, "%s\n", error.what());
      ParallelizationHelperInstance.free();
      return 1;
   }
}
