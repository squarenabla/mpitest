CC = mpicc
CXX = mpicxx

SRC_CPP = $(wildcard *.cpp)
SRC_C = $(wildcard *.c)

OBJ=$(SRC_C:.c=.o) $(SRC_CPP:.cpp=.o)

INCPATH =  -I.
CXXFLAGS =  -march=native -mtune=native -std=c++11 -Ofast
#LIBS = `pkg-config --cflags --libs mpich`

%.o: %.cpp
	$(CXX)  $(INCPATH)  $^ -o $@ -c

all: $(OBJ)
	$(CXX) $(CXXFLAGS)  $(OBJ) -o main
