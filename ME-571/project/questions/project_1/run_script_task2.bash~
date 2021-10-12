#!/bin/bash                                                                   
###                                                                            
###                                                                            
#SBATCH --time=01:00:00                                                        
#SBATCH --tasks=64                                                             
#SBATCH --partition=classroom
                                                  
#SBATCH --job-name=std_dev_task2                                              
#SBATCH --output=std_dev_task2.o%j                                            

module load gcc
module load openmpi/gcc

./std_dev_task2 268435456 >> task2.csv

