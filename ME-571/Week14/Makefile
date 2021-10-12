APPS=sumMatrixOnGPU-2D-grid-2D-block sumMatrix

all: ${APPS}

sumMatrix: sumMatrix.cu kernel.cu
	nvcc -o $@ sumMatrix.cu
%: %.cu
	nvcc -o $@ $<
clean:
	rm -f ${APPS}
