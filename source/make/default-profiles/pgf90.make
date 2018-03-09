# The fortran compiler and its options
#FC = pgf90 -fastsse -Minfo=ccff -Mprof=lines -mp -Mipa -g
FC = pgf90 

FCFLAGS = -fastsse 

#FCDEBUG = -g -pg
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


PREPROC = -Mpreprocess -D__PGIF90__
