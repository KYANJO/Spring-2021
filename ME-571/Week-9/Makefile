CC =gcc
CFLAGS= -Wall -lm
MPICC=mpicc

poisson: poisson.c
	$(MPICC) $(CFLAGS) -o poisson poisson.c

poisson_mpi_sr: poisson_mpi_sr.c
	$(MPICC) $(CFLAGS) -o poisson_mpi_sr poisson_mpi_sr.c

poisson_mpi_nb: poisson_mpi_nb.c
	$(MPICC) $(CFLAGS) -o poisson_mpi_nb poisson_mpi_nb.c

poisson_debug: poisson.c
	$(MPICC) $(CFLAGS) -DDEBUG=1 -o poisson poisson.c

poisson_mpi_sr_debug: poisson_mpi_sr.c
	$(MPICC) $(CFLAGS) -DDEBUG=1 -o poisson_mpi_sr poisson_mpi_sr.c

poisson_mpi_nb_debug: poisson_mpi_nb.c
	$(MPICC) $(CFLAGS) -DDEBUG=1 -o poisson_mpi_nb poisson_mpi_nb.c


clean:
	rm -f poisson
	rm -f poisson_mpi_sr
	rm -f poisson_mpi_nb
