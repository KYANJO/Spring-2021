ls
passwd
logout
ls
cd scratch/
ls
logout
cp -r /cm/shared/examples/ClassExamples/ .
ls
ClassExamples/
ls
cd ClassExamples/
ls
cd gcc
ls
vim main.c 
ls
vim slurm_gcc.bash 
ls
sbatch slurm_gcc.bash 
ls
cat log_slurm.o304104 
ls
vim main.c 
ls
cd ..
ls
cd MPIHelloWorld/
ls
vim helloworld.c 
make
vim makefile 
make
module avail
module load mpich/ge/gcc/64/3.2.1 
module load slurm
module load slurm_MPI.bash 
ls
sbatch slurm_MPI.bash 
ls
cat log_slurm.o304111 
ls
make
ls
cd ..
ls
cd CUDA-Mandelbrot/
ls
vim cudaMandy.cu 
make
module load cuda10.0/toolkit/10.0.130 
make
ls
vim slurm_cuda.bash 
sh slurm_cuda.bash 
ls
sbatch slurm_cuda.bash 
ls
cat log_slurm.o304120 
squeue | grep classtester
cat log_slurm.o304120 
squeue
cd
ls
mkdir week-2
ls
pwd
cd week-2/
pwd
