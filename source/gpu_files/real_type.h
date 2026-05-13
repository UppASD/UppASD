#pragma once

#ifdef HIP_V
    #include <hip/hip_complex.h>
    #ifdef SINGLE_PREC
        typedef float real;
        typedef hipFloatComplex gpu_complex;
        typedef std::complex<float> cpu_complex;

        #define make_hipRealComplex  make_hipFloatComplex
        #define hipCadd_complex     hipCaddf
        #define hipCmul_complex     hipCmulf
        #define hipCreal_complex   hipCrealf
        #define hipCimag_complex    hipCimagf
    #else
        typedef double real;
        typedef hipDoubleComplex gpu_complex;
        typedef std::complex<double> cpu_complex;

        #define make_hipRealComplex  make_hipDoubleComplex
        #define hipCadd_complex     hipCadd
        #define hipCmul_complex     hipCmul
        #define hipCreal_complex    hipCreal
        #define hipCimag_complex    hipCimag
    #endif

#elif defined(CUDA_V)
    #include <thrust/complex.h>
    #ifdef SINGLE_PREC
        typedef float real;
        typedef thrust::complex<float> gpu_complex;
        typedef thrust::complex<float> cpu_complex;

    #else
        typedef double real;
        typedef thrust::complex<double> gpu_complex;
        typedef thrust::complex<double> cpu_complex;

    #endif
#endif