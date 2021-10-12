CC =gcc
CFLAGS= -Wall -lm
MPICC=mpicc

poisson: poisson.c
	$(MPICC) $(CFLAGS) -o poisson poisson.c

poisson_mpi: poisson_mpi.c
	$(MPICC) $(CFLAGS) -o poisson_mpi poisson_mpi.c

poisson_debug: poisson.c
	$(MPICC) $(CFLAGS) -DDEBUG=1 -o poisson poisson.c

poisson_mpi_debug: poisson_mpi.c
	$(MPICC) $(CFLAGS) -DDEBUG=1 -o poisson_mpi poisson_mpi.c


clean:
	rm -f poisson
	rm -f poisson_mpi
