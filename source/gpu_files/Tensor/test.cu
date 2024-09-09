#include "measurements.cuh"
#include <cstdlib>

int main() {
   Matrix<double> x;
   x.AllocateHost(6, 3);

   x(0, 0) = 0.1; x(0, 1) = 0.2; x(0, 2) = 0.3;
   x(1, 0) = 1.1; x(1, 1) = 1.2; x(1, 2) = 1.3;
   x(2, 0) = 2.1; x(2, 1) = 2.2; x(2, 2) = 2.3;
   x(3, 0) = 3.1; x(3, 1) = 3.2; x(3, 2) = 3.3;
   x(4, 0) = 4.1; x(4, 1) = 4.2; x(4, 2) = 4.3;
   x(5, 0) = 5.1; x(5, 1) = 5.2; x(5, 2) = 5.3;

   std::cout << Measurements(x) << std::endl;

   return EXIT_SUCCESS;
}
