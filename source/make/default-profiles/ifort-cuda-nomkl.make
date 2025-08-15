#####################################################################################
## Makefile Settings for ifort profile  														  ##
#####################################################################################

#------------------------------------------------------------------------------------
# The different compilers
#------------------------------------------------------------------------------------

# Fortran compiler
FC = ifort

# C compiler
CC = icc

# C++ compiler
CXX = icc

# CUDA compiler
CUDA = nvcc

#------------------------------------------------------------------------------------
# Flags for FORTRAN compilation
#------------------------------------------------------------------------------------
# Basic optimization settings explained
# -ip         Inline function, substantioal speed up
# -O3         Optimization, faster execution, slow make times
# -ipo        Inline between files
# -xP         Intel processor specific optimizations
# -fast       Uses -ipo -O3 -xP  -static
FCFLAGS = -O3 -ip -xHost -align array16byte

# Basic debug settings explained
# -g           Debug with gdb
# -CB          Array bound checks (same as check)
# -traceback   Gdb traceback
# -p           Function profiling with gprof
# -warn all    Show all warning messages
# -check all   Check all
# -xT          Optimization for intel(R) core(TM)2 Duo
#FCDEBUG = -g -traceback
FCDEBUG =
#------------------------------------------------------------------------------------
# Flags for C compilation
#------------------------------------------------------------------------------------
CCFLAGS = -O3 -g 
CCLIBFLAGS =-qopenmp
# Declare what fortran compiler is used (for C/C++/CUDA code)
C_FCFLAG = -D__Intel__
#------------------------------------------------------------------------------------
# Flags for C++ compilation
#------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------
# Define the path for the GCC
#------------------------------------------------------------------------------------
GCC_ID :=$(shell which icc | sed 's/bin/lib/g')
GCCPATH :=$(dir $(GCC_ID))
CXXFLAGS = -O3 -g
CXXLIBFLAGS = -L${GCCPATH} -lstdc++ -qopenmp 
#------------------------------------------------------------------------------------
# OpenMp related flags (-mp on PGI, -openmp on ifort)
#------------------------------------------------------------------------------------
# Special treatment for ifort compiler to ensure correct [q]openmp flags.
# First find compiler version
IFORTVER := $(shell $(FC) --version | sed 's/\./ /;s/ //;s/[^ ]* *//;s/ .*//;q')
# Then use -openmp if IFORTVER <= 14.x; otherwise use -qopenmp
ifneq ($(shell test $(IFORTVER) -gt 14; echo $$?),0)
FCOMPFLAGS = -openmp -openmp-simd
else
FCOMPFLAGS = -qopenmp
endif

# Library flags
# -lblas       Basic Linear Algebra Subprograms
# -llapack     Linear Algebra Package (for eigenvalue, cholesky etc...)
# -lmkl        Includes lapack and blas
FLIBFLAGS = -lblas -llapack

# ifort mod folder flag (used to put .mods in separate files)
FCMODFLAG = -module
#------------------------------------------------#
# Flags for CUDA compilation
#------------------------------------------------#
# -g	Debug with cuda-gdb
# -G	Debug with cuda-gdb
NVCCFLAGS = -O3

# Sepcific setting for each NIVIDA graphic card
# Other common alternatives:
#  - gencode=arch=compute_70,code=\"sm_70,compute_70\"
GENCODE_ARCH  = -gencode=arch=compute_70,code=\"sm_70,compute_70\"

# CUDA install, include and library paths is not matching default
# This is computer specific. Change if needed.
#CUDA_INSTALL_PATH = /usr/lib/nvidia-cuda-toolkit
#CUDA_INCLUDE_PATH = /usr/lib/nvidia-cuda-toolkit/include
#CUDA_LIBRARY_PATH = /usr/lib/nvidia-304-updates/
#CUDA_INSTALL_PATH = /opt/cuda-7.5
#CUDA_INCLUDE_PATH = /opt/cuda-7.5/include
#CUDA_LIBRARY_PATH = /opt/cuda-7.5/lib64

# Extra libraries needed for compilation
#  - lcublas         BLAS for CUDA
#  - use_fast_math
CUDA_LIB = -lcudart -lcuda -lcurand -lnvToolsExt

# Enable CUDA support
USE_CUDA = YES

# Enable Intel Vector Statistical Library support for RNG
USE_VSL = NO

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

PREPROC = -fpp

# Enable FFTW Support
USE_FFTW = NO
# Enable MKL FFTW Support
USE_MKL_FFT = NO

# Enable OVF support
USE_OVF = NO


