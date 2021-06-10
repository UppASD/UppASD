#############################################
## Makefile Settings for gfortran profile  ##
#############################################

#------------------------------------------------#
# The different compilers
#------------------------------------------------#

# Fortran compiler
FC = gfortran

# C compiler
CC = gcc

# C++ compiler
CXX = g++
#------------------------------------------------#
# Flags for FORTRAN compilation
#------------------------------------------------#
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
# -O3                      Optimization, faster exeuction, slower make times
# -ffree-line-length-200   Allow long lines
FCFLAGS = -O3 -ffree-line-length-0 $(FCARGCHECK)

#------------------------------------------------#
# Flags for C compilation
#------------------------------------------------#
CCFLAGS = -O3 -g -pthread
CCLIBFLAGS = -fopenmp

#------------------------------------------------------------------------------------
# Flags for C++ compilation
#------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------
# Define the path for the GCC
#------------------------------------------------------------------------------------
GCC_ID :=$(shell which g++ | sed 's/bin/lib/g')
GCCPATH :=$(dir $(GCC_ID))
CXXFLAGS = -O3 -g -pthread
CXXLIBFLAGS = -L${GCCPATH} -lstdc++ -fopenmp 
# Basic debug settings explained
# -pedantic
# -p
# -pg
# -fbacktrace
# -Wall
# -Wextra
# -fbounds-check
# -fprofile-arcs
# -ftest-coverage
#FCDEBUG = -g -fbacktrace -Wall -fcheck=bounds
#FCDEBUG = -g -fbacktrace -pg
FCDEBUG =

# OpenMp related flags (-mp on PGI, -openmp on ifort)
FCOMPFLAGS = -fopenmp

# Link flags
# -llapack_atlas lf77blas
# -lcblas
# -lf77blas
# -latlas
# -lblas
# -llapack
# -lmkl
FLIBFLAGS = -llapack -lblas

# gfortran mod flag (used to put .mods in separate files)
FCMODFLAG = -J

# Enable CUDA support
USE_CUDA = NO

# Enable Intel Vector Statistical Library support for RNG
USE_VSL = NO

PREPROC = -cpp

# Enable FFTW support
USE_FFTW = NO
# Enable MKL FFT Support
USE_MKL_FFT = NO

# Enable OVF support
USE_OVF = NO


