#!/bin/bash
#SBATCH -A project_465002525  # Set the allocation to be charged for this job
#SBATCH -J uppasd             # The name of the script is myjob
#SBATCH -t 24:00:00           # 24 hours wall-clock time
#SBATCH -p standard-g         # The partition
#SBATCH -n 1                  # Number of ranks
#SBATCH -c 16                 # Number of cpus per task

export SD_BINARY=/pfs/lustrep2/users/johellsv/GitProjects/UppASD/UppASD_gpu_hip_cu/bin/sd.hip.gpu_gnu
export OMP_NUM_THREADS=16

echo "Script initiated at `date` on `hostname`"

time -p ./runme.sh
#time -p $SD_BINARY > out.log

echo "Script finished at `date` on `hostname`"
