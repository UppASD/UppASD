##########################################
## Makefile Settings for ifort profile  ##
##########################################

#------------------------------------------------#
# The different compilers
#------------------------------------------------#

# Fortran compiler
FC = ftn

#------------------------------------------------#
# Flags for FORTRAN compilation
#------------------------------------------------#
# Basic optimization settings explained
# -ip         Inline function, substantioal speed up
# -O3         Optimization, faster execution, slow make times
# -ipo        Inline between files
# -xP         Intel processor specific optimizations
# -fast       Uses -ipo -O3 -xP  -static
#FCFLAGS = -O3 -ip -openmp -openmp-simd
FCFLAGS = -O3 -h omp

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

# Library flags
# -lblas       Basic Linear Algebra Subprograms
# -llapack     Linear Algebra Package (for eigenvalue, cholesky etc...)
# -lmkl        Includes lapack and blas
FLIBFLAGS = #-mkl=sequential

# ifort mod folder flag (used to put .mods in separate files)
FCMODFLAG = -J

# Enable CUDA support
USE_CUDA = NO

# Enable Intel Vector Statistical Library support for RNG
USE_VSL = NO

PREPROC = -F

# Enable FFTW Support
USE_FFTW = NO
# Enable MKL FFTW Support
USE_MKL_FFT = NO

# Enable OVF support
USE_OVF = NO


