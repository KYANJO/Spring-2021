#
# This is the Makefile for EZ-Spiral.
#
# You must use an ANSI C compiler for EZ-Spiral.
# This Makefile should work with gcc (GNU C) on almost any machine.
# The -Wall option may give you warnings, but this is probably ok.
# If you don't have gcc you should try something like:
# CC=cc -O2 
# or whatever C compiler you want.  You may have to tell cc where 
# the Xlib libraries are, e.g:
#	$(CC) $(OBJ) -o ezspiral -lm -L/usr/X11/lib -lX11

OBJ = ezspiral.o ezstep.o eztip.o ezgraph.o 
CC  = gcc -O2 -Wall 

ezspiral: $(OBJ) 
	$(CC) $(OBJ) -o ezspiral -lm -lX11

clean:	tidy
	rm -f ezspiral *.o *.u *~
	rm -f core *~ .*~ 	
	rm -f *.log *.dvi *.aux

tidy:
	rm -f *.o *.u

# Dependencies
# ------------

ezspiral.o :ezspiral.h 
ezstep.o   :ezspiral.h ezstep.h 
ezgraph.o  :ezspiral.h ezgraph.h 
eztip.o    :ezspiral.h 


