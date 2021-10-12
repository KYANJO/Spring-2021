#!/bin/bash
###
###
#SBATCH --time=01:00:00
#SBATCH --tasks=64
#SBATCH --partition=classroom
#SBATCH --job-name=dudx
#SBATCH --output=dudx.o%j

module load gcc
module load openmpi/gcc
#source ~/.bashrc

echo "nproc, npts, time_v1, error_v1, time_v2, error_v2, time_v3, error_v3"
for p in 1 2 4 8 16 32 64
do
    mpirun -np $p ./derivative_mpi $((2**24)) 
done
