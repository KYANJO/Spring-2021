#!/bin/bash
###
###
#SBATCH --time=01:00:00
#SBATCH --tasks=64
#SBATCH --partition=classroom
#SBATCH --job-name=std_dev
#SBATCH --output=std_dev.o%j

module load gcc
module load openmpi/gcc

echo "nproc,sigma_1,elapsed_time,time_comput,time_comm,time_bcast,time_allreduce,time_reduce"
for p in 1 2 4 8 16 32 64
do
  mpirun -np $p ./std_dev_mpi 268435456
done	 

