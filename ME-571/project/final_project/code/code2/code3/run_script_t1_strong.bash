#!/bin/bash
###
###
#SBATCH --time=01:00:00
#SBATCH --tasks=64
#SBATCH --partition=classroom
#SBATCH --job-name=mc
#SBATCH --output=mc.o%j

module load gcc
module load openmpi/gcc

echo "N,nproc,uexact,umc,time" >> task1_strong.csv
for p in 1 2 4 8 16 32 64
do
  mpirun -np $p ./mc_task1_1 $[2**(28)] >> task1_strong.csv
done	 

