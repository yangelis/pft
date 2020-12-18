#define PFT_IMPLEMENTATION
#include "../pft.hpp"
#include <iostream>
#include <list>

using namespace std;

int main(int argc, char* argv[]) {

  pft::AParse args(argc, argv);
  args.Add({"-h", "--help", "print this help message", false});
  args.Add({"-f", "--file", "-f which file to read", true});

  args.Parse();

  auto t = args.value_of("-f");
  pft::println(stdout, t);

  return 0;
}
