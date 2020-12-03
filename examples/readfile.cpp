#define PFT_IMPLEMENTATION
#include "../pft.hpp"

int main(int argc, char* argv[]) {
  pft::Maybe<pft::StringView> filename;
  if (argc > 1) {
    filename = {true, argv[1]};
  } else {
    fprintf(stderr, "No filename given!\n");
    exit(1);
  }

  auto lines = pft::readlines(filename.unwrap.data);

  pft::pop(lines, 4);
  pft::drop(lines, 4);

  return 0;
}
