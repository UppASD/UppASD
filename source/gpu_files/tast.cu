#include "tensor.cuh"

#include <cstdlib>

int main() {
   double* data = new double[27];
   for(int i = 0; i < 27; ++i) {
      data[i] = i;
   }
   Tensor<double, 3> z;
   z.set(data, Extents<3>{3, 3, 3});
   std::cout << z << std::endl;
   z.resize(Extents<3>{8, 2, 2}, false);
   std::cout << z << std::endl;

   delete[] data;

   return EXIT_SUCCESS;
}

