#!/bin/bash
###
###
#SBATCH --time=01:00:00
#SBATCH --tasks=64
#SBATCH --partition=classroom
#SBATCH --job-name=t3_strong
#SBATCH --output=t3_strong.o%j

module load gcc
module load openmpi/gcc
#source ~/.bashrc

echo "nproc,N,gtime" >> task3_strong.csv
for p in 1 4 9 16 25 36 49
do
mpirun -np $p ./barkley_task3 840 >> task3_strong.csv
done


