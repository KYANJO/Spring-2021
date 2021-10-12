#!/bin/bash
###
###
#SBATCH --time=12:00:00
#SBATCH --tasks=4
#SBATCH --partition=classroom
#SBATCH --job-name=fdiff
#SBATCH --output=fdiff.o%j

module load gcc
module load openmpi/gcc

mpirun -np 4 ./finite_difference 100000000
