CXX = g++
CXXFLAGS = -std=c++17 -O3 -Wall

all: main

main: main.cpp readData.cpp readData.h
	$(CXX) $(CXXFLAGS) main.cpp readData.cpp -o main

clean:
	rm -f main
