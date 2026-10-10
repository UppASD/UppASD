#!/bin/bash
#SBATCH -A naiss2026-3-418-gpu          # Set the allocation to be charged for this job
#SBATCH -J uppasd             # The name of the script is myjob
#SBATCH -t 24:00:00           # 24 hours wall-clock time
#SBATCH -p gpu                # The partition
#SBATCH -G 1                  # Number of GPUs
#SBATCH -n 1                  # Number of tasks
#SBATCH -c 16                 # Number of cpus per task

export SD_BINARY=/home/hellsvik/GitProjects/UppASD/UppASD/bin/sd.arm.cuda
export OMP_NUM_THREADS=16

echo "Script initiated at `date` on `hostname`"

time -p $SD_BINARY > out.log

echo "Script finished at `date` on `hostname`"
