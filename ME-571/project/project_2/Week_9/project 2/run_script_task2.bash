#!/bin/bash
###
###
#SBATCH --time=01:00:00
#SBATCH --tasks=64
#SBATCH --partition=classroom
#SBATCH --job-name=t2_strong
#SBATCH --output=t2_strong.o%j

module load gcc
module load openmpi/gcc
#source ~/.bashrc

echo "nproc,N,gtime" >> task2_strong.csv
for p in 1 4 9 16 25 36 49
do
mpirun -np $p ./barkley_mpi 840 >> task2_strong.csv
done


