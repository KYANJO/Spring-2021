#!/bin/bash
###
###
#SBATCH --time=01:00:00
#SBATCH --tasks=4
#SBATCH --partition=classroom
#SBATCH --job-name=hello
#SBATCH --output=hello.o%j

module load gcc
module load openmpi

mpiexec -np 4 ./hello_mpi
