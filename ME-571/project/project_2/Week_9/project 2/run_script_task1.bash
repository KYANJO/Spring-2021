#!/bin/bash
###
###
#SBATCH --time=01:00:00
#SBATCH --tasks=64
#SBATCH --partition=classroom
#SBATCH --job-name=task1
#SBATCH --output=t1.o%j

module load gcc
module load openmpi/gcc
#source ~/.bashrc

./barkley 100 



