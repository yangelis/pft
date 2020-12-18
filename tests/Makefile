CXXFLAGS = -std=c++17 -I../ -Wall -Wextra -pedantic -O0

.PHONY: test
test: testing
	./testing

testing: testing.cpp ../pft.hpp
	$(CXX) $(CXXFLAGS) $< -o $@
