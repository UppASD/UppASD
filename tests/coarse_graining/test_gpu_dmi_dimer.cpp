#include "gpuHamiltonianDeviceOps.hpp"

#include <cmath>
#include <cstdio>
#include <cstdlib>

#if defined(CUDA_V)
#include <cuda_runtime.h>
#define GPU_CHECK(call) do { const auto error = (call); if (error != cudaSuccess) { \
   std::fprintf(stderr, "%s\\n", cudaGetErrorString(error)); return 1; } } while (false)
#define GPU_NO_DEVICE cudaErrorNoDevice
#define GPU_GET_COUNT cudaGetDeviceCount
#define GPU_MALLOC cudaMalloc
#define GPU_MEMCPY cudaMemcpy
#define GPU_MEMCPY_D2H cudaMemcpyDeviceToHost
#define GPU_FREE cudaFree
#elif defined(HIP_V)
#include <hip/hip_runtime.h>
#define GPU_CHECK(call) do { const auto error = (call); if (error != hipSuccess) { \
   std::fprintf(stderr, "%s\\n", hipGetErrorString(error)); return 1; } } while (false)
#define GPU_NO_DEVICE hipErrorNoDevice
#define GPU_GET_COUNT hipGetDeviceCount
#define GPU_MALLOC hipMalloc
#define GPU_MEMCPY hipMemcpy
#define GPU_MEMCPY_D2H hipMemcpyDeviceToHost
#define GPU_FREE hipFree
#else
#error "This fixture requires the configured CUDA or HIP backend"
#endif

using gpu::hamiltonian::Vec3;
using gpu::hamiltonian::dm_field;
using gpu::hamiltonian::make_vec3;

__global__ void dmi_dimer_kernel(Vec3* fields) {
   // D_12=+D z, M_2=+y; D_21=-D z, M_1=+x.
   fields[0] = dm_field(make_vec3(0, 0, 2.5), make_vec3(0, 1, 0));
   fields[1] = dm_field(make_vec3(0, 0, -2.5), make_vec3(1, 0, 0));
}

int main() {
   int count = 0;
   const auto count_status = GPU_GET_COUNT(&count);
   if (count_status == GPU_NO_DEVICE || count == 0) {
      std::puts("DMI-DIMER-CPU-GPU unavailable: no backend device");
      return 77;
   }
   GPU_CHECK(count_status);

   Vec3* device_fields = nullptr;
   Vec3 fields[2] = {};
   GPU_CHECK(GPU_MALLOC(reinterpret_cast<void**>(&device_fields), sizeof(fields)));
   dmi_dimer_kernel<<<1, 1>>>(device_fields);
#if defined(CUDA_V)
   GPU_CHECK(cudaGetLastError());
   GPU_CHECK(cudaDeviceSynchronize());
#else
   GPU_CHECK(hipGetLastError());
   GPU_CHECK(hipDeviceSynchronize());
#endif
   GPU_CHECK(GPU_MEMCPY(fields, device_fields, sizeof(fields), GPU_MEMCPY_D2H));
   GPU_CHECK(GPU_FREE(device_fields));

   const auto close = [](double actual, double expected) { return std::abs(actual - expected) < 1.0e-12; };
   if (!close(fields[0].x, -2.5) || !close(fields[0].y, 0.0) || !close(fields[0].z, 0.0) ||
       !close(fields[1].x, 0.0) || !close(fields[1].y, -2.5) || !close(fields[1].z, 0.0)) {
      std::fprintf(stderr, "GPU dm_field does not implement D_ij cross M_j\\n");
      return 1;
   }
   std::puts("DMI-DIMER-CPU-GPU passed");
   return 0;
}
