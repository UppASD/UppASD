##########################################
## Makefile Settings for ifort profile  ##
##########################################

#------------------------------------------------#
# The different compilers
#------------------------------------------------#

# Fortran compiler
FC = ifort

# C compiler
CC = icc

# C++ compiler
CXX = icpc

# CUDA compiler
CUDA = nvcc


#------------------------------------------------#
# Flags for FORTRAN compilation 
#------------------------------------------------#
# Basic optimization settings explained
# -ip         Inline function, substantioal speed up
# -O3         Optimization, faster execution, slow make times
# -ipo        Inline between files
# -xP         Intel processor specific optimizations
# -fast       Uses -ipo -O3 -xP  -static
FCFLAGS = -O3 -ip -xHost

# Basic debug settings explained
# -g           Debug with gdb
# -CB          Array bound checks (same as check)
# -traceback   Gdb traceback
# -p           Function profiling with gprof
# -warn all    Show all warning messages
# -check all   Check all
# -xT          Optimization for intel(R) core(TM)2 Duo
FCDEBUG = #-g -traceback

# OpenMp related flags (-mp on PGI, -openmp on ifort)
FCOMPFLAGS = -qopenmp -qopenmp-simd

# Library flags
# -lblas       Basic Linear Algebra Subprograms
# -llapack     Linear Algebra Package (for eigenvalue, cholesky etc...)
# -lmkl        Includes lapack and blas
#
MKLROOT=/pdc/vol/i-compilers/16.0.3/mkl
FLIBFLAGS=-Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_intel_lp64.a ${MKLROOT}/lib/intel64/libmkl_sequential.a ${MKLROOT}/lib/intel64/libmkl_core.a -Wl,--end-group -lpthread -lm -ldl


# ifort mod folder flag (used to put .mods in separate files)
FCMODFLAG = -module

# Declare what fortran compiler is used (for C/C++/CUDA code)
C_FCFLAG = -D__IFORT__


#------------------------------------------------#
# Flags for C compilation 
#------------------------------------------------#
CCFLAGS = -O3 -ip -xHost
CCLIBFLAGS = -lstdc++ -qopenmp 

#------------------------------------------------#
# Flags for C++ compilation
#------------------------------------------------#
CXXFLAGS = -O3 -ip -xHost
CXXLIBFLAGS = -qopenmp  


#------------------------------------------------#
# Flags for CUDA compilation
#------------------------------------------------#
# -g	Debug with cuda-gdb
# -G	Debug with cuda-gdb
NVCCFLAGS = -O3 -arch=sm_37

# Sepcific setting for each NIVIDA graphic card
# Other common alternatives:
GENCODE_ARCH  = -gencode=arch=compute_37,code=\"sm_37,compute_37\" 

# CUDA install, include and library paths is not matching default
# This is computer specific. Change if needed.
CUDA_INSTALL_PATH = /pdc/vol/cuda/cuda-7.5/
CUDA_INCLUDE_PATH = /pdc/vol/cuda/cuda-7.5/include/
CUDA_LIBRARY_PATH = /pdc/vol/cuda/cuda-7.5/lib64/

# Extra libraries needed for compilation
#  - lcublas         BLAS for CUDA
#  - use_fast_math
CUDA_LIB = -lcudart -lcuda -lcurand -lnvToolsExt

# Enable CUDA support
USE_CUDA = YES

# # Enable Intel Vector Statistical Library support for RNG
USE_VSL = YES


#------------------------------------------------#
# Common parameters for C/C++/CUDA code (T/F)
#------------------------------------------------#

# Enable error checking in matrix.hpp
DDEBUG                   = F

# Interrupts program if error occurs in matrix.hpp
DMATRIX_ERROR_INTERRUPT  = F

# Disable stopwatch
DDUMMY_STOPWATCH         = T

# Dont synchronize with device before timing
DASYNC_STOPWATCH         = F

# Big grid (parallelization indexing method)
DUSE_BIG_GRID            = T

# Fast copy (allows the cpu to work ahead of the gpu, requires more memory)
DUSE_FAST_COPY           = T

# Use single preision in C/C++/CUDA code (must also change preisison in parameters.f90)
DSINGLE_PREC             = F

# Profiling (adds measurement events for nvvp)
DNVPROF                  = F

PREPROC = -cpp
