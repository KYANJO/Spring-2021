#!/bin/bash
###
###
#SBATCH --time=01:00:00
#SBATCH --tasks=64
#SBATCH --partition=shortq
#SBATCH --job-name=barkley_mpi
#SBATCH --output=barkley_mpi.o%j

module load gcc
module load openmpi/gcc
#source ~/.bashrc

#./barkley 100
mpirun -np 4 ./barkley_mpi 100
#mpirun -np 4 ./poisson_mpi_nb 100 

