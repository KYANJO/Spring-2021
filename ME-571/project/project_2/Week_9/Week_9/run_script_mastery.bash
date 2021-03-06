#!/bin/bash
###
###
#SBATCH --time=01:00:00
#SBATCH --tasks=64
#SBATCH --partition=classroom
#SBATCH --job-name=mastery_strong
#SBATCH --output=mastery_strong.o%j

module load gcc
module load openmpi/gcc
#source ~/.bashrc

echo "nproc,N,gtime" >> mastery_strong.csv
for p in 1 4 9 16 25 36 49
do
mpirun -np $p ./barkley_mastery 840 >> mastery_strong.csv
done


