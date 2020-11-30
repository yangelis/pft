#include "../pft.hpp"
#include <iostream>
#include <list>
using namespace std;

int main(int argc, char* argv[]) {

  pft::AParse args(argc, argv);
  args.Add({"-h", "--help", "print this help message", false});
  args.Add({"-f", "--file", "-f which file to read", true});

  args.Parse();

  auto t = args.arg_table.find("-f");

  if (t != args.arg_table.end()) {
    cout << "AYY LMAO\n";
    cout << args.arg_table["-f"].second.c_str() << '\n';
  } else {
    cout << "PepeHands\n";
  }

  return 0;
}
