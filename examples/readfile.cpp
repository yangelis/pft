#include "../pft.hpp"
#include <iostream>

auto main(i32 argc, c8* argv[]) -> i32
{
  pft::Maybe<pft::StringView> filename;
  if (argc > 1) {
    filename = {true, argv[1]};
  } else {
    pft::panic("No filename given!\n");
  }

  // auto file = pft::read_file_as_string_view(filename.unwrap.data()).unwrap;
  // pft::println(stdout, file);
  auto lines = pft::readlines(filename.unwrap.data());

  pft::println(stdout, lines);
  // pft::pop(lines, 4);
  // pft::drop(lines, 4);

  return 0;
}
