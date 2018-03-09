#############################################
## Makefile Settings for gfortran profile  ##
#############################################

#------------------------------------------------#
# The different compilers
#------------------------------------------------#

# Fortran compiler
FC = gfortran

#------------------------------------------------#
# Flags for FORTRAN compilation 
#------------------------------------------------#
# Basic optimization settings explained
# -O3                      Optimization, faster exeuction, slower make times
# -ffree-line-length-200   Allow long lines
FCFLAGS = -O3 -ffree-line-length-0

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
#FCDEBUG = -g -fbacktrace -Wall
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
FLIBFLAGS = -framework Accelerate

# gfortran mod flag (used to put .mods in separate files)
FCMODFLAG = -J

# Enable CUDA support
USE_CUDA = NO

# Enable Intel Vector Statistical Library support for RNG
USE_VSL = NO

PREPROC = -cpp
