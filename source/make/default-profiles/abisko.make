# The fortran compiler and its options
FC = mpif90 -O3 -mp -ffast-math -fno-math-errno -OPT:Olimit=0:ro=1:div_split=on:alias=typed:early_mp=on:Ofast



DEBUG = -g

# OpenMp related flags (-mp on PGI, -openmp on ifort)
OMPFLAGS = -mp

# Enable CUDA support
USE_CUDA = NO

# # Enable Intel Vector Statistical Library support for RNG
USE_VSL = NO

PREPROC = -cpp 
