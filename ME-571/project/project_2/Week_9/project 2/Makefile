CC =gcc
CFLAGS= -Wall -lm
MPICC=mpicc

barkley: barkley.c
	$(MPICC) $(CFLAGS) -o barkley barkley.c

barkley_mpi: barkley_mpi.c
	$(MPICC) $(CFLAGS) -o barkley_mpi barkley_mpi.c

barkley_task3: barkley_task3.c
	$(MPICC) $(CFLAGS) -o barkley_task3 barkley_task3.c

barkley_mastery: barkley_mastery.c
	$(MPICC) $(CFLAGS) -o barkley_mastery barkley_mastery.c

barkley_mastery_4: barkley_mastery_4.c
	$(MPICC) $(CFLAGS) -o barkley_mastery_4 barkley_mastery_4.c
clean:
	rm -f barkley
	rm -f barkley_mpi
	rm -f barkley_task3
	rm -f barkley_mastery
	rm -f barkley_mastery_4
