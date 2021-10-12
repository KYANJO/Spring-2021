#!/bin/bash
###
###
#SBATCH --time=01:00:00
#SBATCH --tasks=64
#SBATCH --partition=classroom
#SBATCH --job-name=mpi_task3
#SBATCH --output=mpi_task3.o%j

module load gcc
module load openmpi/gcc

echo "N,time" >> result3_mpi.csv
#for p in 1 2 4 8 16 32 64
#do
mpirun -np 28 ./finite_difference_mpi 1000 >> result3_mpi.csv
mpirun -np 28 ./finite_difference_mpi 100000 >> result3_mpi.csv
mpirun -np 28 ./finite_difference_mpi 1000000 >> result3_mpi.csv
mpirun -np 28 ./finite_difference_mpi 10000000 >> result3_mpi.csv
mpirun -np 28 ./finite_difference_mpi 1000000000 >> result3_mpi.csv



#done	 

