#!/bin/bash
###
###
#SBATCH --time=01:00:00
#SBATCH --tasks=64
#SBATCH --partition=classroom
#SBATCH --job-name=wave
#SBATCH --output=wave.o%j

module load gcc
module load openmpi/gcc
#source ~/.bashrc

mpirun -np 64 ./wave_mpi 1024 0.0005 1 

