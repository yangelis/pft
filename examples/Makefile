ROOTCXXFLAGS := $(shell root-config --cflags)
CXXFLAGS :=-I../ -Wall -Wextra -Werror -pedantic -ggdb -O0
LIBS := $(shell root-config --ldflags --libs)
#LIBS := -L/home/john/Programs/root6.22/lib -lCore -lMathCore
EXAMPLES := file2tree matrix_test argparse

.PHONY: all
all: $(EXAMPLES)

file2tree: file2tree.cpp
	$(CXX) $(CXXFLAGS) $(ROOTCXXFLAGS) file2tree.cpp -o file2tree $(LIBS)

matrix_test: matrix_test.cpp
	$(CXX) $(CXXFLAGS) $(ROOTCXXFLAGS) matrix_test.cpp -o matrix_test $(LIBS)

argparse: argparse.cpp
	$(CXX) $(CXXFLAGS) --std=c++17 argparse.cpp -o argparse

.PHONY: clean
clean:
	rm -f *.o $(EXAMPLES) 
