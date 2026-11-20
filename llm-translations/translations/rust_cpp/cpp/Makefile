CXX ?= g++
CXXFLAGS ?= -std=c++17 -O3 -g -march=native
TARGET := solve

SRCS := main.cpp data.cpp
OBJS := $(SRCS:.cpp=.o)

ifdef FLAT
CXXFLAGS += -DFLAT
endif

$(TARGET): $(OBJS)
	$(CXX) $(CXXFLAGS) -o $@ $(OBJS)

%.o: %.cpp subseq.hpp data.hpp
	$(CXX) $(CXXFLAGS) -c -o $@ $<

clean:
	rm -f $(OBJS) $(TARGET)

.PHONY: clean
