#!/bin/bash
#SBATCH -A pdc.staff          # Set the allocation to be charged for this job
#SBATCH -J uppasd             # The name of the script is myjob
#SBATCH -t 24:00:00           # 2 hours wall-clock time
#SBATCH -p shared             # The partition
#SBATCH -N 1                  # Number of nodes
#SBATCH -c 16                 # Number of cpus per task

export SD_BINARY=/cfs/klemming/home/h/hellsvik/Thora/GitProjects/UppASD-cuda/bin/sd.amdcpu
export OMP_NUM_THREADS=16

echo "Script initiated at `date` on `hostname`"

$SD_BINARY > out.log

echo "Script finished at `date` on `hostname`"
