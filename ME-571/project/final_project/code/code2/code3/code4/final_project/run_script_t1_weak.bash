#!/bin/bash
###
###
#SBATCH --time=01:00:00
#SBATCH --tasks=64
#SBATCH --partition=classroom
#SBATCH --job-name=monte
#SBATCH --output=monte.o%j

module load gcc
module load openmpi/gcc

echo "N,nproc,uexact,umc,time" >> task1_weak.csv

p=(1 2 4 8 16 32 64)
for i in 0 1 2 3 4 5 6
do
  mpirun -np ${p[i]} ./mc_task1_1 $[2**(22+$i)] >> task1_weak.csv
done	 

