# The fortran compiler and its options
#FC = pathf95 -O2 -ff90 -fexceptions -zerouv
#FC = pathf95 -O2
#FC = pathf95 -Ofast -OPT:Olimit=0
FC = pathf95

#FCFLAGS = -O3 -ffast-math -fno-math-errno -OPT:Olimit=0:ro=1:div_split=ON:alias=typed
#FCFLAGS = -O1  -OPT:IEEE_arithmetic=1
FCFLAGS = -O1 -inline -msse3
# O3 too hard..
#FCFLAGS = -O3 -OPT:Ofast -msse3 -inline

#FCDEBUG = -g -pg -C -trapuv -std=c89 -pedantic-errors -Wall -ffortran-bounds-check
FCDEBUG =

# OpenMp related flags (-mp on PGI, -openmp on ifort)
FCOMPFLAGS = -mp

# Link flags
FLIBFLAGS = -lblas -llapack
#
FCMODFLAG = -module

# Enable CUDA support
USE_CUDA = NO

# # Enable Intel Vector Statistical Library support for RNG
USE_VSL = NO

# Enable FFTW Support
USE_FFTW = NO
# Enable MKL FFTW Support
USE_MKL_FFT = NO
