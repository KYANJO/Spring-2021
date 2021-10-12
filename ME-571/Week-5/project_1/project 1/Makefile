std_dev: std_dev.c
	mpicc -o std_dev -lm std_dev.c

std_dev_mpi: std_dev_mpi.c
	mpicc -o std_dev_mpi -lm std_dev_mpi.c

std_dev_serial: std_dev_serial.c
	gcc -o std_dev_serial -lm std_dev_serial.c

std_dev_parallel: std_dev_parallel.c
	mpicc -o std_dev_parallel -lm std_dev_parallel.c

clean:
	rm std_dev std_dev_mpi std_dev_seial std_dev_parallel
