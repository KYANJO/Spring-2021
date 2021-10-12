#!/bin/bash                                                                   
###                                                                            
###                                                                            
#SBATCH --time=01:00:00                                                        
#SBATCH --tasks=64                                                             
#SBATCH --partition=classroom
                                                  
#SBATCH --job-name=std_dev_serial                                              
#SBATCH --output=std_dev_serial.o%j                                            

module load gcc
module load openmpi/gcc

./std_dev_serial 268435456


#mpirun -np 1 ./std_dev_serial  268435456
