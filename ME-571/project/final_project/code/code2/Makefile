mc_task1_1: mc_task1_1.c
	mpicc -o mc_task1_1 -lm mc_task1_1.c

mc_task2: mc_task2.cu
	nvcc -o mc_task2 -lm mc_task2.cu

mc_task3: mc_task3.cu
	nvcc -o mc_task3 -lm mc_task3.cu

trial: trial.cu
	nvcc -o trial -lm trial.cu

clean:
	rm mc_task1_1 mc_task2 mc_task3
