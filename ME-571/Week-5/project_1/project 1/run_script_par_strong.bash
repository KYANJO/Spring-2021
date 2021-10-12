#!/bin/bash
###
###
#SBATCH --time=01:00:00
#SBATCH --tasks=64
#SBATCH --partition=classroom
#SBATCH --job-name=std_dev_par_st
#SBATCH --output=std_dev_par_st.o%j

module load gcc
module load openmpi/gcc

echo "nproc,mu,sigma,elapsed_time,time_comput,time_comm,time_gather,time_reduce"
for p in 1 2 4 8 16 32 64
do
  mpirun -np $p ./std_dev_parallel 268435456
done	 

