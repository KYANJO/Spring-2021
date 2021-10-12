#!/bin/bash
###
###
#SBATCH --time=01:00:00
#SBATCH --tasks=64
#SBATCH --partition=shortq
#SBATCH --job-name=monte
#SBATCH --output=output.o%j

module load gcc
module load openmpi/gcc
#source ~/.bashrc

#echo "N,uexact,umc,time" >> task1.csv
#for q in 3 4 5 6 7 8 9 10
#do
mpirun -np 64 ./mc_task1_2 $((10**10)) #$((10**$q)) >> task1.csv
#done
