#!/bin/bash
###
###
#SBATCH --time=01:00:00
#SBATCH --tasks=64
#SBATCH --partition=classroom
#SBATCH --job-name=mastery_4
#SBATCH --output=mastery_4.o%j

module load gcc
module load openmpi/gcc
#source ~/.bashrc

echo "nproc,Nloc,tcomm,tcomp" >> mastery_4.csv
for p in 1 4 9 16 25 36 49
do
mpirun -np $p ./barkley_mastery_4 840 >> mastery_4.csv
done


