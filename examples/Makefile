ROOTCXXFLAGS	:= `root-config --cflags`
CXXFLAGS			:= ${ROOTCXXFLAGS} -Wall -ggdb -O0
LIBS					:= `root-config --ldflags --libs`
# LIBS 					:= -L/home/john/Programs/root6.22/lib -lCore -lMathCore

.PHONY: all clean
all: file2tree matrix_test

file2tree: file2tree.cpp
	$(CXX) $(CXXFLAGS) file2tree.cpp -o file2tree $(LIBS)

matrix_test: matrix_test.cpp
	$(CXX) $(CXXFLAGS) matrix_test.cpp -o matrix_test $(LIBS)

clean:
	rm -f *.o file2tree matrix_test
