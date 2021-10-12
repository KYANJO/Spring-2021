PROG = steady

HEADERS = driver.h operators.h solver.h ezplot.h
SRC = driver.c operators.c solver.c realft.c 
OBJ= driver.o solver.o operators.o realft.o 

# SGI
#CFLAGS = -O2 
#LIBRARIES = -lcomplib.sgimath -lm 
#LDFLAGS = -L$(HOME)/lib -L/usr/local/lib $(LIBRARIES)

# Linux: set LIBRARIES to use ezplot graphics or not
# Atlas with graphics
# LIBRARIES = -L/usr/local/lib -llapack -lblas -latlas -lf2c -lezplot -lX11 -lm
# Atlas without graphics
# LIBRARIES = -L/usr/local/lib -llapack -lblas -latlas -lf2c -lm 
# standard blas/lapack without graphics

LIBRARIES = -llapack -lblas -lf2c -lm 
CFLAGS = -g -O2 -Wall 
LDFLAGS = -L$(HOME)/lib -L/usr/X11R6/lib $(LIBRARIES)

# Fortran (not necessary with direct solves) 
F77 = fort77
FFLAGS = -g -C -Wall 

.c.o :
	$(CC) -c $(CFLAGS) $(NOPTION) $*.c

.f.o :
	$(F77) -c $(FFLAGS) $(NOPTION) $*.f

$(PROG): $(OBJ) 
	$(CC) $(CFLAGS) $(OBJ) -o $(PROG) $(LDFLAGS)
#	$(F77) $(CFLAGS) $(OBJ) -o $(PROG) $(LDFLAGS)


# Dependencies
# ------------

driver.o     :driver.h
solver.o     :solver.h driver.h
operators.o  :operators.h driver.h

clean:	tidy
	rm -f $(PROG)  *.o *~
	rm -f core *~ .*~ 	
	rm -f *.log *.dvi *.aux

tidy:
	rm -f $(PROG) *.o 

tar:
	tar -cvf steady.tar $(SRC) $(HEADERS) Makefile dgbtr_manpage \
	ic.dat.example.fife.gz  ic.dat.example.db.gz 


