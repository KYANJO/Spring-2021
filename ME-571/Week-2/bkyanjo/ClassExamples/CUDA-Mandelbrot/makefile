#module load cuda10.0
EXEC=cudaMandy 

all: $(EXEC)

cudaMandy: cudaMandy.cu 
	nvcc -o $@ cudaMandy.cu $(OBJS)



clean: 
	rm -f *.o $(OBJS) $(EXEC)
	rm -f mandelbrotCuda.log
	rm -f Mandelbrot.pgm
 
