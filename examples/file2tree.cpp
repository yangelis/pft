#include "../pft.hpp"
#include <iostream>

#include <TF1.h>

using namespace std;

void generate_file(const char* filename) {
  FILE *f = fopen(filename, "w");
  fprintf(f, "This is the header\n with some text\n to simulate a real world "
             "case\n\n");
  auto sigmoid =
      new TF1("sigmoid", "1 / (1.0 + exp(-[1] * (x - [0])))", 0.0, 30.0);
  sigmoid->SetParameters(2, 0.25);

  for (float x = 0.0; x <= 30.0; x += 0.3) {
    fprintf(f, "%f\n", (*sigmoid)(x));
  }

  fclose(f);
}

int main(int argc, char *argv[]) {
  const char * filename = "data.txt";
  generate_file(filename);
  auto buffer = pft::read_file_as_string_view(filename);
  if (!buffer.has_value) {
    cerr << "Could not read file" << '\n';
    exit(1);
  }
  auto vec = split_by(buffer.unwrap, '\n');
  for (int i = 0; i < 10; i++) {
    cout << vec[i] << '\n';
  }
  cout << "\nafter deleting header\n";
  pft::ignore_header_lines(vec, 4);
  auto buf = pft::as_float(vec);
  for (int i = 0; i < 10; i++) {
    cout << buf[i] << '\n';
  }
  return 0;
}
