#!/bin/bash
###
###
#SBATCH --time=01:00:00
#SBATCH --tasks=4
#SBATCH --partition=classroom
#SBATCH --job-name=std_dev
#SBATCH --output=std_dev.o%j

module load gcc
module load openmpi/gcc
#source ~/.bashrc

mpirun -np 1 ./std_dev_mpi data_test.csv

