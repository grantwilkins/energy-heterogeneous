CC := g++
CFLAGS	:=-O3 -std=c++14

LIB_SZ3 := -L /home/grant/software/SZ3/install/lib 
INC_SZ3 := -I /home/grant/software/SZ3/install/include
LIB_PAPI := -L /home/grant/software/papi/src/install/lib -lpapi
INC_PAPI := -I /home/grant/software/papi/src/install/include

LIBS	:= $(LIB_PAPI) $(LIB_SZ3)
INC     := $(INC_PAPI) $(INC_SZ3)

SOURCES := sz3-energy \

all: $(SOURCES)
.PHONY: all

sz3-energy: sz3-energy.cpp
	$(CC) $(CFLAGS) $^ -o $@ $(INC) $(LIBS)
	
clean:
	rm -f $(SOURCES)