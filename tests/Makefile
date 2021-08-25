CXXFLAGS      := -std=c++17 -I../ -Wall -Wextra -Wshadow \
                 -Wpedantic -Wcast-align -ggdb -O3
SANITIZE      := -fsanitize=address -fsanitize=leak -fsanitize=undefined

.PHONY: test
test: testing
	./testing

testing: testing.cpp ../pft.hpp
	$(CXX) $(CXXFLAGS) -fno-omit-frame-pointer $(SANITIZE) -lrt $< -o $@
