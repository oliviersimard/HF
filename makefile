SRC=$(PWD)/src
OBJ=$(PWD)/obj

CXX=g++ -std=c++14 -Wall 
INC=${HOME}/include
LIB=${HOME}/lib
CXXFLAGS= -L$(LIB) -ljson_spirit -lfftw3 -larmadillo

PROG=PiTon

all: $(PROG)

$(PROG): $(OBJ)/main.o $(OBJ)/json_utils.o $(OBJ)/fft.o $(OBJ)/param.o $(OBJ)/Hubbard.o $(OBJ)/integral_utils.o
	$(CXX) $(OBJ)/main.o $(OBJ)/json_utils.o $(OBJ)/fft.o $(OBJ)/param.o $(OBJ)/Hubbard.o $(OBJ)/integral_utils.o $(CXXFLAGS) -o main.out

$(OBJ)/main.o: $(SRC)/Main.cpp
	$(CXX) -c -I$(INC) $(SRC)/Main.cpp -o $(OBJ)/main.o

$(OBJ)/json_utils.o: $(SRC)/json_utils.cpp
	$(CXX) -c -I$(INC) $(SRC)/json_utils.cpp -o $(OBJ)/json_utils.o

$(OBJ)/fft.o: $(SRC)/fft.cpp
	$(CXX) -c -I$(INC) $(SRC)/fft.cpp -o $(OBJ)/fft.o

$(OBJ)/param.o: $(SRC)/param.cpp
	$(CXX) -c -I$(INC) $(SRC)/param.cpp -o $(OBJ)/param.o

$(OBJ)/Hubbard.o: $(SRC)/Hubbard.cpp
	$(CXX) -c -I$(INC) $(SRC)/Hubbard.cpp -o $(OBJ)/Hubbard.o

$(OBJ)/integral_utils.o: $(SRC)/integral_utils.cpp
	$(CXX) -c -I$(INC) $(SRC)/integral_utils.cpp -o $(OBJ)/integral_utils.o
clean:
	rm $(OBJ)/* $(PWD)/main.out

