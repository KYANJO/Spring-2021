CU_APPS=dot_product


all: ${CU_APPS}

%: %.cu
	nvcc -o $@ $<
clean:
	rm -f ${CU_APPS}
