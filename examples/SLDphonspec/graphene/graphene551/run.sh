#!/bin/bash
#SBATCH -A snic2021-5-304
#SBATCH -t 30:00:00
#SBATCH -n 32 
#SBATCH -J fcc-Pd


module load VASP/5.4.4.16052018-nsc2-intel-2018a-eb
mpprun vasp_ncl
