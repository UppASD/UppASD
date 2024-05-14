#!/bin/bash -l
# Job script for VASP

#SBATCH -A pdc.staff
#SBATCH --reservation=sbtest-s11
#SBATCH -p mtn
#SBATCH -J vaspjob
#SBATCH -t 1-00:00:00
# number of nodes
#SBATCH --nodes=4
# Number of MPI processes per node
#SBATCH --ntasks-per-node=128

# Ignore the LMOD cache
export LMOD_IGNORE_CACHE=1

ml PDCTEST/22.06
ml vasp/6.3.2-vanilla

# Different version of libfabric on login node uan02 and on the compute nodes
export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:/opt/cray/libfabric/1.15.0.0/lib64

echo "Script initiated at `date` on `hostname`"
time srun -n 512 vasp > out.log
echo "Script finished at `date` on `hostname`"
