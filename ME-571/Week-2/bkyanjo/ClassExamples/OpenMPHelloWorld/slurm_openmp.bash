#!/bin/bash

#SBATCH -J OpenMP_TEST     # job name
#SBATCH -o log_slurm.o%j  # output and error file name (%j expands to jobID)
#SBATCH -n 4           # total number of tasks requested
#SBATCH -N 1 		  # number of nodes you want to run on	
#SBATCH -p classroom  # queue (partition) -- defq, eduq, gpuq, shortq
#SBATCH -t 12:00:00       # run time (hh:mm:ss) - 12.0 hours in this example.
# Generally needed modules:
module load slurm

#module load matlab/r2020a
# Execute the program
export OMP_NUM_THREADS=4
./hellow

## Some examples:
# mpirun python3 script.py

# Exit if mpirun errored:
status=$?
if [ $status -ne 0 ]; then
    exit $status
fi

# Do some post processing.
