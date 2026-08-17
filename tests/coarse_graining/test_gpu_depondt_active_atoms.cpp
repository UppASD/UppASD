// RCG-09A.2: the active-list path must be the production Depondt algorithm,
// not an adaptive approximation.  This fixture compares the observable
// predictor/corrector buffers bitwise against the ordinary full-range calls.

#include "gpuDepondtIntegrator.hpp"
#include "gpuParallelizationHelper.hpp"
#include "gpuStructures.hpp"
#include "gpuThermfield.hpp"
#include "tensor.hpp"

#include <cmath>
#include <cstdio>
#include <cstring>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

constexpr std::size_t atoms = 5;
constexpr std::size_t ensembles = 2;
constexpr std::size_t vectorCount = 3 * atoms * ensembles;

bool backendAvailable() {
   int count = 0;
#if defined(CUDA_V)
   const auto status = cudaGetDeviceCount(&count);
   // Some scheduler/container combinations report a stale visible count until
   // the first runtime call.  Require a usable context, not merely a count,
   // so this fixture follows the project's SKIP_RETURN_CODE contract instead
   // of turning an unavailable accelerator into a false regression failure.
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

SimulationParameters parameters() {
   SimulationParameters result{};
   result.stt = 'N';
   result.N = atoms;
   result.M = ensembles;
   result.delta_t = real(0.017);
   result.gamma = real(1.3);
   result.k_bolt = real(1);
   result.mub = real(1);
   result.damping = real(0.11);
   result.randomSeed = 918273ULL;
#if defined(CUDA_V)
   result.rngType = CURAND_RNG_PSEUDO_DEFAULT;
#else
   result.rngType = HIPRAND_RNG_PSEUDO_DEFAULT;
#endif
   return result;
}

struct State {
   deviceLattice lattice;

   State() {
      lattice.beff.Allocate(3, atoms, ensembles);
      lattice.b2eff.Allocate(3, atoms, ensembles);
      lattice.eneff.Allocate(3, atoms, ensembles);
      lattice.emomM.Allocate(3, atoms, ensembles);
      lattice.emom.Allocate(3, atoms, ensembles);
      lattice.emom2.Allocate(3, atoms, ensembles);
      lattice.mmom.Allocate(atoms, ensembles);
      lattice.btorque.Allocate(3, atoms, ensembles);
      lattice.b2eff.zeros();
      lattice.eneff.zeros();
      lattice.btorque.zeros();
   }

   ~State() {
      lattice.beff.Free();
      lattice.b2eff.Free();
      lattice.eneff.Free();
      lattice.emomM.Free();
      lattice.emom.Free();
      lattice.emom2.Free();
      lattice.mmom.Free();
      lattice.btorque.Free();
   }

   State(const State&) = delete;
   State& operator=(const State&) = delete;
};

std::vector<real> directions() {
   std::vector<real> result(vectorCount);
   for(std::size_t ensemble = 0; ensemble < ensembles; ++ensemble) {
      for(std::size_t atom = 0; atom < atoms; ++atom) {
         const real x = real(0.13) * real(atom + 1);
         const real y = real(-0.07) * real(ensemble + 1);
         const real z = real(1) - real(0.03) * real(atom);
         const real norm = std::sqrt(x * x + y * y + z * z);
         const std::size_t element = 3 * (atom + atoms * ensemble);
         result[element] = x / norm;
         result[element + 1] = y / norm;
         result[element + 2] = z / norm;
      }
   }
   return result;
}

std::vector<real> field(real scale) {
   std::vector<real> result(vectorCount);
   for(std::size_t element = 0; element < result.size(); ++element)
      result[element] = scale * real(0.01) * real(1 + element);
   return result;
}

void upload(State& state, const std::vector<real>& direction,
            const std::vector<real>& effectiveField) {
   const std::vector<real> moments(atoms * ensembles, real(1));
   check(GPU_MEMCPY(state.lattice.emom.data(), direction.data(), direction.size() * sizeof(real),
                    GPU_MEMCPY_HOST_TO_DEVICE), "emom upload");
   check(GPU_MEMCPY(state.lattice.emom2.data(), direction.data(), direction.size() * sizeof(real),
                    GPU_MEMCPY_HOST_TO_DEVICE), "emom2 upload");
   check(GPU_MEMCPY(state.lattice.emomM.data(), direction.data(), direction.size() * sizeof(real),
                    GPU_MEMCPY_HOST_TO_DEVICE), "emomM upload");
   check(GPU_MEMCPY(state.lattice.mmom.data(), moments.data(), moments.size() * sizeof(real),
                    GPU_MEMCPY_HOST_TO_DEVICE), "mmom upload");
   check(GPU_MEMCPY(state.lattice.beff.data(), effectiveField.data(),
                    effectiveField.size() * sizeof(real), GPU_MEMCPY_HOST_TO_DEVICE),
         "beff upload");
}

void uploadField(State& state, const std::vector<real>& effectiveField) {
   check(GPU_MEMCPY(state.lattice.beff.data(), effectiveField.data(),
                    effectiveField.size() * sizeof(real), GPU_MEMCPY_HOST_TO_DEVICE),
         "predictor beff upload");
}

std::vector<real> download(const GpuTensor<real, 3>& tensor) {
   std::vector<real> result(static_cast<std::size_t>(tensor.size()));
   check(GPU_MEMCPY(result.data(), tensor.data(), result.size() * sizeof(real),
                    GPU_MEMCPY_DEVICE_TO_HOST), "tensor download");
   return result;
}

void requireBitwiseEqual(const std::vector<real>& actual, const std::vector<real>& expected,
                         const char* label) {
   if(actual.size() != expected.size() ||
      std::memcmp(actual.data(), expected.data(), actual.size() * sizeof(real)) != 0)
      throw std::runtime_error(std::string("Depondt active-list mismatch in ") + label);
}

void requireSelectedUnchanged(const std::vector<real>& before, const std::vector<real>& after,
                              const std::vector<int>& oneBasedActive) {
   for(std::size_t ensemble = 0; ensemble < ensembles; ++ensemble) {
      for(std::size_t atom = 0; atom < atoms; ++atom) {
         bool selected = false;
         for(const int active : oneBasedActive)
            selected = selected || active == static_cast<int>(atom + 1);
         if(selected) continue;
         const std::size_t element = 3 * (atom + atoms * ensemble);
         if(std::memcmp(before.data() + element, after.data() + element,
                        3 * sizeof(real)) != 0)
            throw std::runtime_error("Depondt active-list modified an inactive atom");
      }
   }
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

void runAllActiveParity(real temperatureValue, const char* label) {
   const auto inputDirection = directions();
   const auto initialField = field(real(1));
   const auto predictorField = field(real(-0.6));
   const std::vector<int> all = {1, 2, 3, 4, 5};
   DeviceList deviceAll(all);
   State ordinary;
   State active;
   upload(ordinary, inputDirection, initialField);
   upload(active, inputDirection, initialField);

   Tensor<real, 1> temperature;
   temperature.AllocateHost(atoms);
   for(std::size_t atom = 0; atom < atoms; ++atom) temperature(atom) = temperatureValue;
   const auto params = parameters();
   GpuDepondtIntegrator ordinaryIntegrator;
   GpuDepondtIntegrator activeIntegrator;
   if(!ordinaryIntegrator.initiate(params) || !activeIntegrator.initiate(params) ||
      !ordinaryIntegrator.initiateConstants(params, temperature) ||
      !activeIntegrator.initiateConstants(params, temperature))
      throw std::runtime_error("Depondt fixture initialization failed");

   ordinaryIntegrator.evolveFirst(ordinary.lattice);
   activeIntegrator.evolveFirst(active.lattice, deviceAll.values.data(), all.size());
   check(GPU_STREAM_SYNC(ParallelizationHelperInstance.getWorkStream()), "predictor synchronization");
   requireBitwiseEqual(download(active.lattice.emom), download(ordinary.lattice.emom), label);
   requireBitwiseEqual(download(active.lattice.emom2), download(ordinary.lattice.emom2), label);
   requireBitwiseEqual(download(active.lattice.emomM), download(ordinary.lattice.emomM), label);
   requireBitwiseEqual(download(active.lattice.b2eff), download(ordinary.lattice.b2eff), label);

   uploadField(ordinary, predictorField);
   uploadField(active, predictorField);
   ordinaryIntegrator.evolveSecond(ordinary.lattice);
   activeIntegrator.evolveSecond(active.lattice, deviceAll.values.data(), all.size());
   check(GPU_STREAM_SYNC(ParallelizationHelperInstance.getWorkStream()), "corrector synchronization");
   requireBitwiseEqual(download(active.lattice.emom), download(ordinary.lattice.emom), label);
   requireBitwiseEqual(download(active.lattice.emom2), download(ordinary.lattice.emom2), label);
   requireBitwiseEqual(download(active.lattice.emomM), download(ordinary.lattice.emomM), label);
   requireBitwiseEqual(download(active.lattice.b2eff), download(ordinary.lattice.b2eff), label);
   temperature.FreeHost();
}

void runReorderedSubsetIdentity() {
   const auto inputDirection = directions();
   const auto initialField = field(real(1));
   const auto predictorField = field(real(-0.6));
   const std::vector<int> forward = {2, 5};
   const std::vector<int> reverse = {5, 2};
   DeviceList deviceForward(forward);
   DeviceList deviceReverse(reverse);
   State first;
   State second;
   upload(first, inputDirection, initialField);
   upload(second, inputDirection, initialField);
   Tensor<real, 1> temperature;
   temperature.AllocateHost(atoms);
   for(std::size_t atom = 0; atom < atoms; ++atom) temperature(atom) = real(0.75);
   const auto params = parameters();
   GpuDepondtIntegrator firstIntegrator;
   GpuDepondtIntegrator secondIntegrator;
   if(!firstIntegrator.initiate(params) || !secondIntegrator.initiate(params) ||
      !firstIntegrator.initiateConstants(params, temperature) ||
      !secondIntegrator.initiateConstants(params, temperature))
      throw std::runtime_error("Depondt subset fixture initialization failed");

   firstIntegrator.evolveFirst(first.lattice, deviceForward.values.data(), forward.size());
   secondIntegrator.evolveFirst(second.lattice, deviceReverse.values.data(), reverse.size());
   uploadField(first, predictorField);
   uploadField(second, predictorField);
   firstIntegrator.evolveSecond(first.lattice, deviceForward.values.data(), forward.size());
   secondIntegrator.evolveSecond(second.lattice, deviceReverse.values.data(), reverse.size());
   check(GPU_STREAM_SYNC(ParallelizationHelperInstance.getWorkStream()), "subset synchronization");
   requireBitwiseEqual(download(first.lattice.emom2), download(second.lattice.emom2),
                       "reordered active list");
   requireBitwiseEqual(download(first.lattice.emomM), download(second.lattice.emomM),
                       "reordered active list");
   requireSelectedUnchanged(inputDirection, download(first.lattice.emom2), forward);
   temperature.FreeHost();
}

// CGP-06B negative control (task's own suggested control #1): a careless
// active-scoped write might index the drawn random value by compact-list
// slot instead of by physical site. The production kernel
// (active_atom_site_kernel, gpuParallelizationHelper.hpp) never does this --
// it derives site from oneBasedAtoms[slot] before indexing -- but this toy
// kernel deliberately reproduces the bug to prove the permutation-invariance
// test above is actually discriminating, not vacuously passing.
__global__ void buggySlotIndexedWrite(real* field, const real* randomValues,
                                      const int* oneBasedAtoms, unsigned int activeCount) {
   const unsigned int slot = blockIdx.x * blockDim.x + threadIdx.x;
   if(slot >= activeCount) return;
   const int oneBasedAtom = oneBasedAtoms[slot];
   const unsigned int atom = static_cast<unsigned int>(oneBasedAtom - 1);
   field[atom] = randomValues[slot];  // BUG: should be randomValues[atom]
}

void runNegativeControlSlotIndexing() {
   const std::vector<real> randomHost = {real(11), real(22), real(33), real(44), real(55)};
   GpuTensor<real, 1> randomValues;
   randomValues.Allocate(static_cast<index_t>(randomHost.size()));
   {
      Tensor<real, 1> view(const_cast<real*>(randomHost.data()),
                           static_cast<index_t>(randomHost.size()));
      randomValues.copy_sync(view);
   }
   GpuTensor<real, 1> fieldForward, fieldReverse;
   fieldForward.Allocate(static_cast<index_t>(atoms));
   fieldReverse.Allocate(static_cast<index_t>(atoms));
   const std::vector<real> poison(atoms, real(-9999));
   check(GPU_MEMCPY(fieldForward.data(), poison.data(), atoms * sizeof(real),
                    GPU_MEMCPY_HOST_TO_DEVICE), "forward poison upload");
   check(GPU_MEMCPY(fieldReverse.data(), poison.data(), atoms * sizeof(real),
                    GPU_MEMCPY_HOST_TO_DEVICE), "reverse poison upload");

   const std::vector<int> forward = {2, 5};
   const std::vector<int> reverse = {5, 2};
   DeviceList deviceForward(forward);
   DeviceList deviceReverse(reverse);

   buggySlotIndexedWrite<<<1, 8, 0, ParallelizationHelperInstance.getWorkStream()>>>(
      fieldForward.data(), randomValues.data(), deviceForward.values.data(), 2);
   buggySlotIndexedWrite<<<1, 8, 0, ParallelizationHelperInstance.getWorkStream()>>>(
      fieldReverse.data(), randomValues.data(), deviceReverse.values.data(), 2);
   check(GPU_STREAM_SYNC(ParallelizationHelperInstance.getWorkStream()),
        "slot-indexing negative control synchronization");

   std::vector<real> hostForward(atoms), hostReverse(atoms);
   check(GPU_MEMCPY(hostForward.data(), fieldForward.data(), atoms * sizeof(real),
                    GPU_MEMCPY_DEVICE_TO_HOST), "forward download");
   check(GPU_MEMCPY(hostReverse.data(), fieldReverse.data(), atoms * sizeof(real),
                    GPU_MEMCPY_DEVICE_TO_HOST), "reverse download");

   // Site 2 (index 1) and site 5 (index 4) must disagree between forward and
   // reverse compact-list order under the slot-indexed bug; if they agree,
   // this control is not discriminating.
   if(hostForward[1] == hostReverse[1] && hostForward[4] == hostReverse[4])
      throw std::runtime_error(
         "negative control not discriminating: slot-indexed bug did not manifest");

   fieldForward.Free();
   fieldReverse.Free();
   randomValues.Free();
}

// CGP-06B: the active-scoped randomize() overload must write exactly the
// active (site,ensemble) slots and leave every other slot untouched. Proven
// by poisoning the field buffer before the call (GPU_MALLOC does not
// zero-initialize) and checking active-list entries match a same-seeded,
// full-write reference bit-for-bit (same generator sequence -- CGP-06A
// invariant 1's mechanism) while every other entry still holds the poison.
void runPartialActiveFieldScoping() {
   const auto params = parameters();
   Tensor<real, 1> temperature;
   temperature.AllocateHost(atoms);
   for(std::size_t atom = 0; atom < atoms; ++atom) temperature(atom) = real(0.6);

   GpuThermfield fieldFull;
   GpuThermfield fieldScoped;
   if(!fieldFull.initiate(atoms, ensembles, params.rngType, params.randomSeed) ||
      !fieldScoped.initiate(atoms, ensembles, params.rngType, params.randomSeed) ||
      !fieldFull.initiateConstants(temperature, params.delta_t, params.gamma, params.k_bolt,
                                   params.mub, params.damping) ||
      !fieldScoped.initiateConstants(temperature, params.delta_t, params.gamma, params.k_bolt,
                                     params.mub, params.damping))
      throw std::runtime_error("Thermfield scoping fixture initialization failed");

   GpuTensor<real, 2> mmom;
   mmom.Allocate(static_cast<index_t>(atoms), static_cast<index_t>(ensembles));
   {
      const std::vector<real> moments(atoms * ensembles, real(1));
      check(GPU_MEMCPY(mmom.data(), moments.data(), moments.size() * sizeof(real),
                       GPU_MEMCPY_HOST_TO_DEVICE), "scoping mmom upload");
   }

   const std::vector<real> poison(vectorCount, real(-9999.0));
   check(GPU_MEMCPY(const_cast<real*>(fieldScoped.getField().data()), poison.data(),
                    vectorCount * sizeof(real), GPU_MEMCPY_HOST_TO_DEVICE), "poison upload");

   const std::vector<int> active = {2, 5};
   DeviceList deviceActive(active);

   fieldFull.randomize(mmom);
   fieldScoped.randomize(mmom, deviceActive.values.data(),
                         static_cast<unsigned int>(active.size()));
   check(GPU_STREAM_SYNC(ParallelizationHelperInstance.getWorkStream()),
        "scoping synchronization");

   const auto fullHost = download(fieldFull.getField());
   const auto scopedHost = download(fieldScoped.getField());

   for(std::size_t ensemble = 0; ensemble < ensembles; ++ensemble) {
      for(std::size_t atom = 0; atom < atoms; ++atom) {
         const std::size_t oneBased = atom + 1;
         const bool isActive = oneBased == 2 || oneBased == 5;
         const std::size_t element = 3 * (atom + atoms * ensemble);
         if(isActive) {
            if(std::memcmp(scopedHost.data() + element, fullHost.data() + element,
                           3 * sizeof(real)) != 0)
               throw std::runtime_error("Active-scoped thermfield write diverged from full write");
         } else {
            for(int component = 0; component < 3; ++component)
               if(scopedHost[element + static_cast<std::size_t>(component)] != real(-9999.0))
                  throw std::runtime_error("Active-scoped thermfield wrote an inactive site");
         }
      }
   }

   mmom.Free();
   temperature.FreeHost();
}

// CGP-06B negative control (task's own suggested control #2): Strategy 1
// keeps the generator's full-N,M draw every step specifically so a physical
// atom's stochastic continuation across a coarse<->fine transition never
// depends on how long it spent coarse. This shows why: skipping the
// generate call entirely on a hypothetical all-coarse step -- the naive
// optimization this task's own text warns against -- desynchronizes the
// generator from what an always-advancing reference produces on the very
// next call.
void runNegativeControlSkippedAdvancement() {
   const auto params = parameters();
   Tensor<real, 1> temperature;
   temperature.AllocateHost(atoms);
   for(std::size_t atom = 0; atom < atoms; ++atom) temperature(atom) = real(0.6);

   GpuThermfield reference;
   GpuThermfield skipped;
   if(!reference.initiate(atoms, ensembles, params.rngType, params.randomSeed) ||
      !skipped.initiate(atoms, ensembles, params.rngType, params.randomSeed) ||
      !reference.initiateConstants(temperature, params.delta_t, params.gamma, params.k_bolt,
                                   params.mub, params.damping) ||
      !skipped.initiateConstants(temperature, params.delta_t, params.gamma, params.k_bolt,
                                 params.mub, params.damping))
      throw std::runtime_error("Thermfield advancement fixture initialization failed");

   GpuTensor<real, 2> mmom;
   mmom.Allocate(static_cast<index_t>(atoms), static_cast<index_t>(ensembles));
   {
      const std::vector<real> moments(atoms * ensembles, real(1));
      check(GPU_MEMCPY(mmom.data(), moments.data(), moments.size() * sizeof(real),
                       GPU_MEMCPY_HOST_TO_DEVICE), "advancement mmom upload");
   }

   // "step 1": reference advances, as production always does. skipped
   // deliberately omits the call, simulating a hypothetical
   // skip-generate-during-coarse-residence optimization.
   reference.randomize(mmom);
   // "step 2": both draw again -- skipped's call is its first-ever draw.
   reference.randomize(mmom);
   skipped.randomize(mmom);
   check(GPU_STREAM_SYNC(ParallelizationHelperInstance.getWorkStream()),
        "advancement synchronization");

   const auto referenceHost = download(reference.getField());
   const auto skippedHost = download(skipped.getField());
   if(std::memcmp(referenceHost.data(), skippedHost.data(),
                 referenceHost.size() * sizeof(real)) == 0)
      throw std::runtime_error(
         "negative control not discriminating: skipped advancement did not desynchronize");

   mmom.Free();
   temperature.FreeHost();
}

} // namespace

int main() {
   if(!backendAvailable()) {
      std::puts("DEPONDT-ACTIVE-ATOMS unavailable: no backend device");
      return 77;
   }
   try {
      ParallelizationHelperInstance.initiate(atoms, ensembles, 0);
      runAllActiveParity(real(0), "all-active T=0");
      runAllActiveParity(real(0.75), "all-active finite temperature");
      runReorderedSubsetIdentity();
      runPartialActiveFieldScoping();
      runNegativeControlSlotIndexing();
      runNegativeControlSkippedAdvancement();
      ParallelizationHelperInstance.free();
      std::puts("DEPONDT-ACTIVE-ATOMS passed");
      return 0;
   } catch(const std::exception& error) {
      std::fprintf(stderr, "%s\n", error.what());
      ParallelizationHelperInstance.free();
      return 1;
   }
}
