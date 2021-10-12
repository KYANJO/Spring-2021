#!/bin/bash
###
###
#SBATCH --time=01:00:00
#SBATCH --tasks=1
#SBATCH --partition=classroom
#SBATCH --job-name=trapez
#SBATCH --output=trapez.o%j

module load gcc
module load openmpi/gcc

./trapez 100
