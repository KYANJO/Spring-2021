APPS=poisson 

all: ${APPS}

%: %.cu
	nvcc -o $@ $<
clean:
	rm -f ${APPS}
