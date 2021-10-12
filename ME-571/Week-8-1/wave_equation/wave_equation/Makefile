CC =gcc
CFLAGS= -Wall -lm
MPICC=mpicc

wave: wave.c
	$(MPICC) $(CFLAGS) -o wave wave.c

wave_mpi: wave_mpi.c
	$(MPICC) $(CFLAGS) -o wave_mpi wave_mpi.c

wave_debug: wave.c
	$(MPICC) $(CFLAGS) -DDEBUG=1 -o wave wave.c

wave_mpi_debug: wave_mpi.c
	$(MPICC) $(CFLAGS) -DDEBUG=1 -o wave_mpi wave_mpi.c


clean:
	rm -f wave
	rm -f wave_mpi
