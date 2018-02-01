
TARGET = lf30bus
CC = icc
CFLAGS = -std=c++11 -fopenmp -xCORE-AVX2 -axMIC-AVX512 -O3  

DEP = LF_include.h LF_functions.h
OBJ = LF30.o functions.o

%.o: %.cpp $(DEP)
	$(CC) -o $(CFLAGS)
lf30bus: $(OBJ)
	$(CC) -o $@ $^ $(CFLAGS)

all: $(TARGET)

run: all
	- ./lf30bus
clean:
	rm -rf *.o lf30bus	


