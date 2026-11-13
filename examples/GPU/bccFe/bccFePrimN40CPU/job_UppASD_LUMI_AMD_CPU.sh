#!/bin/bash
#SBATCH -A <project>          # Set the allocation to be charged for this job
#SBATCH -J uppasd             # The name of the script is myjob
#SBATCH -t 24:00:00           # 24 hours wall-clock time
#SBATCH -p small              # The partition
#SBATCH -n 1                  # Number of tasks
#SBATCH -c 16                 # Number of cpus per task

# The executable is in the bin directory of the build directory, e.g. <UppASD dir>/build_AMD_CPU/bin/sd
export SD_BINARY=<UppASD executable>
export OMP_NUM_THREADS=16

echo "Script initiated at `date` on `hostname`"

time -p ./runme.sh
#time -p $SD_BINARY > out.log

echo "Script finished at `date` on `hostname`"
