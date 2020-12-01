CXXFLAGS      := -std=c++17 -I../ -Wall -Wextra -Werror -pedantic -ggdb -O0
LIBS          :=
ROOTCXXFLAGS  := $(shell root-config --cflags)
ROOTLIBS      := $(shell root-config --ldflags --libs) -lCore -lMathCore
EXAMPLES      := file2tree matrix_test argparse

.PHONY: all
all: $(EXAMPLES)

file2tree: file2tree.cpp
	@printf "# compiling $<\n"
	$(CXX) $(CXXFLAGS) $(ROOTCXXFLAGS) file2tree.cpp -o file2tree $(ROOTLIBS)

matrix_test: matrix_test.cpp
	@printf "# compiling $<\n"
	$(CXX) $(CXXFLAGS) matrix_test.cpp -o matrix_test

argparse: argparse.cpp
	@printf "# compiling $<\n"
	$(CXX) $(CXXFLAGS) argparse.cpp -o argparse

.PHONY: clean
clean:
	@printf "Cleaning executables\n"
	rm -f *.o $(EXAMPLES) 
	rm -f *.root 
