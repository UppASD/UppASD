// RCG-09A.2: the active-list path must be the production Depondt algorithm,
// not an adaptive approximation.  This fixture compares the observable
// predictor/corrector buffers bitwise against the ordinary full-range calls.

#include "gpuDepondtIntegrator.hpp"
#include "gpuParallelizationHelper.hpp"
#include "gpuStructures.hpp"
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
      ParallelizationHelperInstance.free();
      std::puts("DEPONDT-ACTIVE-ATOMS passed");
      return 0;
   } catch(const std::exception& error) {
      std::fprintf(stderr, "%s\n", error.what());
      ParallelizationHelperInstance.free();
      return 1;
   }
}
