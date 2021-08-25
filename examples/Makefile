CXXFLAGS      := -std=c++17 -I../ -Wall -Wextra -Wshadow \
               -Wpedantic -Wcast-align -ggdb -O0
SANITIZE      := -fsanitize=address -fsanitize=leak -fsanitize=undefined
SANITIZE_LIBS := -static-libasan
LIBS          :=
ROOTCXXFLAGS  := $(CXXFLAGS) $(shell root-config --cflags)
ROOTLIBS      := $(LIBS) $(shell root-config --ldflags --libs) -lCore -lMathCore
ROOTEXAMPLES  := file2tree
EXAMPLES_FFTW := peak_finding sg_filter ffts fftw_wrapper
EXAMPLES      := matrix_test readfile gol

.PHONY: all
all: $(EXAMPLES) $(ROOTEXAMPLES) $(EXAMPLES_FFTW)

.PHONY: examples
examples: $(EXAMPLES)

$(EXAMPLES): %:%.cpp
	@printf "# compiling $<\n"
	$(CXX) $(CXXFLAGS) -fno-omit-frame-pointer -fsanitize=address -lrt $< -o $@

$(ROOTEXAMPLES): %:%.cpp
	@printf "# compiling $<\n"
	$(CXX) $(ROOTCXXFLAGS) $< -o $@ $(ROOTLIBS)

$(EXAMPLES_FFTW): %: %.cpp
	@printf "# compiling $<\n"
	$(CXX) $(ROOTCXXFLAGS) $< -o $@ $(ROOTLIBS) -lfftw3

convolutions: convolutions.cpp
	@printf "# compiling $<\n"
	$(CXX) $(CXXFLAGS) $< -o $@ -lfftw3

placeholder: placeholder.cpp
	@printf "# compiling $<\n"
	$(CXX) $(CXXFLAGS) -fno-omit-frame-pointer $(SANITIZE) -lrt $< -o $@

.PHONY: clean
clean:
	@printf "Cleaning executables\n"
	rm -f *.o $(EXAMPLES)
	rm -f *.o $(ROOTEXAMPLES)
	rm -f *.o $(EXAMPLES_FFTW)
	rm -f *.root
