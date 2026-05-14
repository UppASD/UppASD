#pragma once

#ifdef HIP_V
    #ifdef ON_LUMI
        // LUMI (AMD GPU): use thrust::complex for improved interoperability with CUDA code
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
        #define GPU_CREAL(x)             ((x).real())
        #define GPU_CIMAG(x)             ((x).imag())
        #define MAKE_GPU_COMPLEX(re, im) gpu_complex(re, im)
    #else
        // Standard HIP: use native hipComplex types
        #include <hip/hip_complex.h>
        #ifdef SINGLE_PREC
            typedef float real;
            typedef hipFloatComplex gpu_complex;
            typedef std::complex<float> cpu_complex;

            #define make_hipRealComplex  make_hipFloatComplex
            #define hipCadd_complex     hipCaddf
            #define hipCmul_complex     hipCmulf
            #define hipCreal_complex    hipCrealf
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
        #define GPU_CREAL(x)             hipCreal_complex(x)
        #define GPU_CIMAG(x)             hipCimag_complex(x)
        #define MAKE_GPU_COMPLEX(re, im) make_hipRealComplex(re, im)
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
    #define GPU_CREAL(x)             ((x).real())
    #define GPU_CIMAG(x)             ((x).imag())
    #define MAKE_GPU_COMPLEX(re, im) gpu_complex(re, im)
#endif
