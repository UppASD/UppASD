#####################################################################################
## Makefile Settings for ifort profile  														  ##
#####################################################################################

#------------------------------------------------------------------------------------
# The different compilers
#------------------------------------------------------------------------------------

# Fortran compiler
FC = gfortran

# C compiler
CC = gcc

# C++ compiler
CXX = g++

#------------------------------------------------------------------------------------
# Flags for FORTRAN compilation
#------------------------------------------------------------------------------------
# First check compiler version for version dependent flag setting
GFORTVER := $(shell $(FC) --version | head -1 | awk '{ print $$NF}' | sed 's/\./ /g' | awk '{ print $$1}')
$(info $$GFORTVER is [$(GFORTVER)] [$(FC)])
# Then use -fallow-argument-mismatch if GFORTVER > 8.x; otherwise use flag is not needed
ifneq ($(shell test $(GFORTVER) -gt 9; echo $$?),0)
FCARGCHECK =
else
FCARGCHECK = -fallow-argument-mismatch
endif

# Basic optimization settings explained
# -ip         Inline function, substantioal speed up
# -O3         Optimization, faster execution, slow make times
# -ipo        Inline between files
# -xP         Intel processor specific optimizations
# -fast       Uses -ipo -O3 -xP  -static
FCFLAGS =  -O3 -ffree-line-length-0 $(FCARGCHECK)

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
CCFLAGS = -O3 -g -pthread
CCLIBFLAGS =-fopenmp -llapack
# Declare what fortran compiler is used (for C/C++/CUDA code)
C_FCFLAG = -D__GFORTRAN__
#------------------------------------------------------------------------------------
# Flags for C++ compilation
#------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------
# Define the path for the GCC
#------------------------------------------------------------------------------------
GCC_ID :=$(shell which g++ | sed 's/bin/lib/g')
GCCPATH :=$(dir $(GCC_ID))
CXXFLAGS = -O3 -g -pthread
CXXLIBFLAGS = -L${GCCPATH} -lstdc++ -fopenmp -Wl,-no_compact_unwind 
# OpenMp related flags (-mp on PGI, -openmp on ifort)
FCOMPFLAGS = -fopenmp 
#------------------------------------------------------------------------------------
# Library flags
#------------------------------------------------------------------------------------
# -lblas       Basic Linear Algebra Subprograms
# -llapack     Linear Algebra Package (for eigenvalue, cholesky etc...)
# -lmkl        Includes lapack and blas
FLIBFLAGS = -framework Accelerate
# ifort mod folder flag (used to put .mods in separate files)
FCMODFLAG = -J

PREPROC = -cpp

# Enable FFTW Support
USE_FFTW = NO
# Enable MKL FFT Support
USE_MKL_FFT = NO 

# Enable CUDA support
USE_CUDA = NO

# Enable Intel Vector Statistical Library support for RNG
USE_VSL = NO
# Enable OVF support
USE_OVF = NO


