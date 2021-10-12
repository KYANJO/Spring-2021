#!/bin/bash
###
###
#SBATCH --time=01:00:00
#SBATCH --tasks=4
#SBATCH --partition=classroom
#SBATCH --job-name=p2phello
#SBATCH --output=p2phello.o%j

module load gcc
module load openmpi/gcc
#source ~/.bashrc

mpirun -np 4 ./p2p_hello_world

