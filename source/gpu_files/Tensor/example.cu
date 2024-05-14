#include <cstdlib>
#include <memory>

#include "tensor.cuh"


// CudaTensors must be passed by value(uses shallow copy)
template <typename T, index_t dim>
__global__ void Compute(const CudaTensor<T, dim> A, CudaTensor<T, dim> B) {
   int i = blockDim.x * blockIdx.x + threadIdx.x;
   for(; i < A.size(); i += gridDim.x * blockDim.x) {
      B[i] = 1000 * A[i];
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
   std::unique_ptr<double[]> data{new double[ext.size()]};

   Tensor<double, 3> A;
   A.set(data.get(), ext);
   initialize(A);

   CudaTensor<double, 3> A_gpu, B_gpu;
   A_gpu.Allocate(ext);
   B_gpu.Allocate(ext);
   A_gpu.copy_sync(A);

   Compute<<<2, 2>>>(A_gpu, B_gpu);  // multiply B = A * 1000
   A.copy_sync(B_gpu);
   std::cout << A << std::endl;

   A_gpu.Free();
   B_gpu.Free();

   return EXIT_SUCCESS;
}
