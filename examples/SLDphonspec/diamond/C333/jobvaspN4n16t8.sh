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
#SBATCH --ntasks-per-node=16
# Note that cpus-per-task is set as 2x OMP_NUM_THREADS
#SBATCH --cpus-per-task=16

# Ignore the LMOD cache
export LMOD_IGNORE_CACHE=1

# Number and placement of OpenMP threads
export OMP_NUM_THREADS=8
export OMP_PLACES=cores
export OMP_PROC_BIND=false

ml PDCTEST/22.06
ml vasp/6.3.2-vanilla

# Different version of libfabric on login node uan02 and on the compute nodes
export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:/opt/cray/libfabric/1.15.0.0/lib64

echo "Script initiated at `date` on `hostname`"
time srun -n 64 vasp > out.log
echo "Script finished at `date` on `hostname`"
