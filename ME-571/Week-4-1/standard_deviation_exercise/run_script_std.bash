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

./std_run men_data.csv
#mpirun -np 1 ./std_run_mpi women_data.csv

