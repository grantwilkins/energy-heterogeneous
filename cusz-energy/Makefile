NVCC	:=nvcc --cudart=static -ccbin g++ --extended-lambda
CFLAGS	:=-O3 -std=c++14
ARCHES	:=-gencode arch=compute_61,code=\"compute_61,sm_61\"

LIB_PAPI := -L /home/grant/software/papi/src/install/lib -lpapi
INC_PAPI := -I /home/grant/software/papi/src/install/include
LIB_CUSZ := -L /home/grant/software/cusz-latest/build -lcusz
INC_CUSZ := -I /home/grant/software/cusz-latest/include 


LIBS	:= $(LIB_PAPI) $(LIB_CUSZ)
INC     := $(INC_PAPI) $(INC_CUSZ)

SOURCES := cusz-energy \

all: $(SOURCES)
.PHONY: all

cusz-energy: cusz-energy.cu
	$(NVCC) $(CFLAGS) ${ARCHES} $^ -o $@ $(INC) $(LIBS)
	
clean:
	rm -f $(SOURCES)
