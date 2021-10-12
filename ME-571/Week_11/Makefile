APPS=add sumArraysOnHost sumArraysOnGPU sumArraysOnGPU-timer checkDeviceInfor

all: ${APPS}

%: %.cu
	nvcc -o $@ $<
clean:
	rm -f ${APPS}
