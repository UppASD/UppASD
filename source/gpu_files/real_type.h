#pragma once

#if defined(HIP_V) 
    #include <hip/hip_complex.h>
    #ifdef SINGLE_PREC
        typedef float real;
        typedef hipFloatComplex hipComplex_t;

        #define make_hipComplex  make_hipFloatComplex
        #define hipCadd_complex     hipCaddf
        #define hipCmul_complex     hipCmulf
        #define hipCreal_complex   hipCrealf
        #define hipCimag_complex    hipCimagf
    #else
        typedef double real;
        typedef hipDoubleComplex hipComplex_t;

        #define make_hipComplex  make_hipDoubleComplex
        #define hipCadd_complex     hipCadd
        #define hipCmul_complex     hipCmul
        #define hipCreal_complex    hipCreal
        #define hipCimag_complex    hipCimag
    #endif
#elif defined(CUDA_V)
    #ifdef SINGLE_PREC
        typedef float real;

    #else
        typedef double real;
    #endif
#endif