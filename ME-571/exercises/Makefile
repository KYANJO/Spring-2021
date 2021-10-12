finite_difference: finite_difference.c
	gcc -o finite_difference -lm finite_difference.c

finite_difference_mpi: finite_difference_mpi.c
	mpicc -o finite_difference_mpi -lm finite_difference_mpi.c

finite_difference_cuda: finite_difference_cuda.cu
	nvcc -o finite_difference_cuda finite_difference_cuda.cu

task3_1: task3_1.cu
	nvcc -o task3_1 task3_1.cu

task3_2: task3_2.cu
	nvcc -o task3_2 task3_2.cu
