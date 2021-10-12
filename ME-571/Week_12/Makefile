APPS=finite_difference finite_difference1 checkDeviceInfor sumArraysOnGPU-timer sumMatrixOnGPU-1D-grid-1D-block sumMatrixOnGPU-2D-grid-1D-block sumMatrixOnGPU-2D-grid-2D-block sumMatrixOnGPU

all: ${APPS}

%: %.cu
	nvcc -o $@ $<
clean:
	rm -f ${APPS}
