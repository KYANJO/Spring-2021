#module load mpich/ge/gcc/64/3.2.1
hellow: helloworld.c
	mpicc -o $@ helloworld.c

clean: 
	rm hellow log_slurm.*
