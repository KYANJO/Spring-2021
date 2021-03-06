#!/bin/bash                                                                          
###                                                                                  
###                                                                                  
#SBATCH --time=01:00:00                                                              
#SBATCH --tasks=64                                                                   
#SBATCH --partition=classroom                                                        
#SBATCH --job-name=wave                                                              
#SBATCH --output=wave.o%j                                                            

module load gcc
module load openmpi/gcc
#source ~/.bashrc                                                                    

echo "nproc, N, u_error, gtime" >> strong.csv
for p in 1 2 4 8 16 32 64
do
    mpirun -np $p ./wave_mpi 1024 0.0005 1 >> strong.csv
done



