#!/bin/bash
#SBATCH -A <project>          # Set the allocation to be charged for this job
#SBATCH -J uppasd             # The name of the script is myjob
#SBATCH -t 24:00:00           # 24 hours wall-clock time
#SBATCH -p small-g            # The partition
#SBATCH -G 1                  # Number of GPUs
#SBATCH -n 1                  # Number of ranks
#SBATCH -c 14                 # Number of cpus per task

# The executable is in the bin directory of the build directory, e.g. <UppASD dir>/build_GPU_GNU/bin/sd.hip
export SD_BINARY=<UppASD executable>
export OMP_NUM_THREADS=14

echo "Script initiated at `date` on `hostname`"

time -p ./runme.sh
#time -p $SD_BINARY > out.log

echo "Script finished at `date` on `hostname`"
