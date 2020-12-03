CXXFLAGS      := -std=c++17 -I../ -Wall -Wextra -Werror -pedantic -ggdb -O0
LIBS          :=
ROOTCXXFLAGS  := $(shell root-config --cflags)
ROOTLIBS      := $(shell root-config --ldflags --libs) -lCore -lMathCore
ROOTEXAMPLES  := file2tree
EXAMPLES      := matrix_test argparse readfile

.PHONY: all
all: $(EXAMPLES) $(ROOTEXAMPLES)

$(EXAMPLES): %: %.cpp
	@printf "# compiling $<\n"
	$(CXX) $(CXXFLAGS)  $^ -o $@

$(ROOTEXAMPLES): %: %.cpp
	@printf "# compiling $<\n"
	$(CXX) $(CXXFLAGS) $(ROOTCXXFLAGS) $^ -o $@ $(ROOTLIBS)

.PHONY: clean
clean:
	@printf "Cleaning executables\n"
	rm -f *.o $(EXAMPLES)
	rm -f *.o $(ROOTEXAMPLES)
	rm -f *.root 
