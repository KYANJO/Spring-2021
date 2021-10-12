#!/bin/bash
###
###
#SBATCH --time=01:00:00
#SBATCH --tasks=64
#SBATCH --partition=classroom
#SBATCH --job-name=std_dev_weak
#SBATCH --output=std_dev_weak.o%j

module load gcc
module load openmpi/gcc

echo "nproc,sigma_1,elapsed_time,time_comput,time_comm,time_bcast,time_allreduce,time_reduce" >> task1_weak_mastery.csv

p=(1 2 4 8 16 32 64)
for i in 0 1 2 3 4 5 6
do
  mpirun -np ${p[i]} ./std_dev_task1 $[2**(50+$i)] >> task1_weak_mastery.csv
done	 

