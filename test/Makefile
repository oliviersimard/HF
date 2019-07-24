SRC=$(PWD)/../src

CXX=g++ -std=c++14 -Wall 
INC=${HOME}/include
LIB=${HOME}/lib
CXXFLAGS= -L$(LIB) -ljson_spirit -lfftw3 -larmadillo

PROG=HubbardUnitTests

all: $(PROG)

$(PROG): $(OBJ)/main.o $(OBJ)/json_utils.o $(OBJ)/fft.o $(OBJ)/param.o $(OBJ)/Hubbard.o $(OBJ)/integral_utils.o $(OBJ)/HubbardSelfCon.o
	$(CXX) $(OBJ)/main.o $(OBJ)/json_utils.o $(OBJ)/fft.o $(OBJ)/param.o $(OBJ)/Hubbard.o $(OBJ)/integral_utils.o $(OBJ)/HubbardSelfCon.o $(CXXFLAGS) -o main.out

$(OBJ)/main.o: $(SRC)/Main.cpp
	$(CXX) -c -I$(INC) $(SRC)/Main.cpp -o $(OBJ)/main.o

$(OBJ)/json_utils.o: $(SRC)/json_utils.cpp
	$(CXX) -c -I$(INC) $(SRC)/json_utils.cpp -o $(OBJ)/json_utils.o
