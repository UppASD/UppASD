#!/bin/bash
#SBATCH -N 1
#SBATCH -t 0:59:00
#SBATCH -J bccFe     
#SBATCH -U SNIC001-10-192
#SBATCH --mail-type=ALL  
#SBATCH --mail-user=andrea.taroni@fysik.uu.se

: "${SD_BINARY:?Set SD_BINARY to the full path of the UppASD executable}"

export OMP_NUM_THREADS=8
time ${SD_BINARY} > out.dm
gzip --best *.out
