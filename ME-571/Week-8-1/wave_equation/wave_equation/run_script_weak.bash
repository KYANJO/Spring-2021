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
  
echo "nproc, N, u_error, gtime" >> weak.csv
p=(1 2 4 8 16 32 64)
for i in 0 1 2 3 4 5 6
do
   tfinal=0.5
   dt=0.0005
   var=$[2**$i]
   d=$(echo "scale=10; $dt/$var" | bc -l)
   mpirun -np ${p[i]} ./wave_mpi $[2**(10+$i)] $d $tfinal >> weak.csv 
done
