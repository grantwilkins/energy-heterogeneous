CC := gcc
CFLAGS	:=-O3 

LIB_SZX := -L /home/grant/software/SZx-main/install/lib -lSZx
INC_SZX := -I /home/grant/software/SZx-main/install/include
LIB_PAPI := -L /home/grant/software/papi/src/install/lib -lpapi
INC_PAPI := -I /home/grant/software/papi/src/install/include

LIBS	:= $(LIB_PAPI) $(LIB_SZX) -lm
INC     := $(INC_PAPI) $(INC_SZX)

SOURCES := szx-energy \

all: $(SOURCES)
.PHONY: all

szx-energy: szx-energy.c
	$(CC) $(CFLAGS) $^ -o $@ $(INC) $(LIBS)
	
clean:
	rm -f $(SOURCES)