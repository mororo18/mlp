CXX      ?= g++
CXXFLAGS ?= -std=c++17 -O3 -march=native -Wall

SRC := $(wildcard src/*.cpp)
OBJ := $(SRC:.cpp=.o)
BIN := solve

all: $(BIN)

$(BIN): $(OBJ)
	$(CXX) $(CXXFLAGS) -o $@ $(OBJ)

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	rm -f $(OBJ) $(BIN)

.PHONY: all clean
