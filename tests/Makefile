CXXFLAGS      := -std=c++17 -I../ -fno-exceptions -Wall -Wextra -Wshadow \
                 -Wpedantic -Wcast-align -ggdb -O0

.PHONY: test
test: testing
	./testing

testing: testing.cpp ../pft.hpp
	$(CXX) $(CXXFLAGS) $< -o $@
