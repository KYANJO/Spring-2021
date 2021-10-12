APPS=poisson kernels

all: ${APPS}

%: %.cu
	nvcc -o $@ $<
clean:
	rm -f ${APPS}
