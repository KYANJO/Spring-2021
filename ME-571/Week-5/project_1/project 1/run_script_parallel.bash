#!/bin/bash                                                                   
###                                                                            
###                                                                            
#SBATCH --time=01:00:00                                                        
#SBATCH --tasks=64                                                             
#SBATCH --partition=classroom
                                                  
#SBATCH --job-name=std_dev_parallel                                              
#SBATCH --output=std_dev_parallel.o%j                                            

module load gcc
module load openmpi/gcc

echo "mean,sigma"
for p in 1 2 4 8 16 32 64
do
    mpirun -np $p ./std_dev_parallel $[2**28]
  # printf>>output.txt
done
