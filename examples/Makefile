CXXFLAGS      := -std=c++17 -I../ -fno-exceptions -Wall -Wextra -Wshadow \
               -Wpedantic -Wcast-align -ggdb -O0
LIBS          :=
ROOTCXXFLAGS  := $(shell root-config --cflags)
ROOTLIBS      := $(shell root-config --ldflags --libs) -lCore -lMathCore
ROOTEXAMPLES  := file2tree
EXAMPLES_FFTW := peak_finding sg_filter ffts
EXAMPLES      := matrix_test argparse readfile

.PHONY: all
all: $(EXAMPLES) $(ROOTEXAMPLES) $(EXAMPLES_FFTW)

.PHONY: examples
examples: $(EXAMPLES)

$(EXAMPLES): %:%.cpp
	@printf "# compiling $<\n"
	$(CXX) $(CXXFLAGS) $< -o $@

$(ROOTEXAMPLES): %:%.cpp
	@printf "# compiling $<\n"
	$(CXX) $(CXXFLAGS) $(ROOTCXXFLAGS) $< -o $@ $(ROOTLIBS)

$(EXAMPLES_FFTW): %: %.cpp
	@printf "# compiling $<\n"
	$(CXX) $(CXXFLAGS) $(ROOTCXXFLAGS) $< -o $@ $(ROOTLIBS) -lfftw3

.PHONY: clean
clean:
	@printf "Cleaning executables\n"
	rm -f *.o $(EXAMPLES)
	rm -f *.o $(ROOTEXAMPLES)
	rm -f *.o $(EXAMPLES_FFTW)
	rm -f *.root 
