#!/bin/bash
#SBATCH -N 1
#SBATCH -t 0:59:00
#SBATCH -J bccFe     
#SBATCH -U SNIC001-10-192
#SBATCH --mail-type=ALL  
#SBATCH --mail-user=andrea.taroni@fysik.uu.se

export OMP_NUM_THREADS=8
time /home/x_andta/UppASD_3.0/source/sd > out.dm
gzip --best *.out
