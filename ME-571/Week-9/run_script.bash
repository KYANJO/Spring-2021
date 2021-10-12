#!/bin/bash
###
###
#SBATCH --time=01:00:00
#SBATCH --tasks=64
#SBATCH --partition=classroom
#SBATCH --job-name=poisson
#SBATCH --output=poisson.o%j

module load gcc
module load openmpi/gcc
#source ~/.bashrc

./poisson 100
#mpirun -np 4 ./poisson_mpi_sr 100
#mpirun -np 4 ./poisson_mpi_nb 100 

