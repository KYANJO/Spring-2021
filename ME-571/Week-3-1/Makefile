finite_difference: finite_difference.c
	gcc -o finite_difference -lm finite_difference.c

finite_difference_mpi: finite_difference_mpi.c
	mpicc -o finite_difference_mpi -lm finite_difference_mpi.c

trapez_integration: trapez_integration.c
	gcc -o trapez_integration -lm trapez_integration.c

clean:
	rm finite_difference finite_difference_mpi
