FC = ftn -O3 -hfp3                   #Cray Fortran (recommended)
#FC = ftn -fast -openmp          #Intel (do not use !)
#FC = ftn -fastsse -mp -Mipa   #PGI (forbidden .... really slow)

# MPI related flags (usually -lfmpich -lmpich -lpthread if needed)
MPIFLAGS =

# OpenMp related flags (-mp on PGI, -openmp on ifort)
OMPFLAGS =

DEBUG = -g

# Enable CUDA support
USE_CUDA = NO

# Enable Intel Vector Statistical Library support for RNG
USE_VSL = NO

PREPROC = -cpp

# Enable FFTW Support
USE_FFTW = NO
# Enable MKL FFTW Support
USE_MKL_FFT = NO

# Enable OVF support
USE_OVF = NO


