#include <cstdlib>
#include <memory>

#include "tensor.cuh"


// CudaTensors must be passed by value(uses shallow copy)
template <typename T, index_t dim>
__global__ void Compute(CudaTensor<T, dim> A) {
   int i = blockDim.x * blockIdx.x + threadIdx.x;
   for(; i < A.size(); i += gridDim.x * blockDim.x) {
      A[i] += 1000;
   }
}


// initialize receives shallow copy
template <typename T, index_t dim>
void initialize(Tensor<T, dim> A) {
   for(int i = 0; i < A.size(); ++i) {
      A[i] = i;
   }
}


int main() {
   Extents<3> ext(3, 3, 3);
   Tensor<double, 3> A;
   A.AllocateHost(ext);
   initialize(A);

   CudaTensor<double, 3> A_gpu;
   A_gpu.Allocate(ext);

   A_gpu.copy_async(A);
   Compute<<<2, 2, 0, 0>>>(A_gpu);
   A.copy_async(A_gpu);
   cudaDeviceSynchronize();

   std::cout << A << std::endl;

   A_gpu.Free();
   A.FreeHost();
   return EXIT_SUCCESS;
}
