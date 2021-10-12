#!/bin/bash
###
###
#SBATCH --time=01:00:00
#SBATCH --tasks=4
#SBATCH --partition=classroom
#SBATCH --job-name=trapez
#SBATCH --output=trapez.o%j

module load gcc
module load openmpi/gcc

mpirun -np 4 ./trapez_integration_mpi 100000000
