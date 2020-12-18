#define PFT_IMPLEMENTATION
#include "../pft.hpp"
#include <iostream>
#include <memory>

#include <TF1.h>
#include <TFile.h>
#include <TTree.h>

using namespace std;

void generate_file(const char* filename) {
  FILE* f = fopen(filename, "w");
  fprintf(f, "This is the header\nwith some text\nto simulate a real world "
             "case\n\n");
  auto sigmoid =
      new TF1("sigmoid", "1 / (1.0 + exp(-[1] * (x - [0])))", 0.0, 30.0);
  sigmoid->SetParameters(2, 0.25);

  for (float x = 0.0; x <= 30.0; x += 0.3) {
    fprintf(f, "%f\n", (*sigmoid)(x));
  }

  fclose(f);
}

template <typename T>
shared_ptr<TTree> createTree(vector<T>& vec) {
  shared_ptr<TTree> tree = make_shared<TTree>("tree", "Tree from a file");
  tree->Branch("vec", &vec);
  tree->Fill();
  return tree;
}

int main(int argc, char* argv[]) {
  pft::Maybe<pft::StringView> filename;
  if (argc > 1) {
    filename = {true, argv[1]};
  } else {
    fprintf(stderr, "No filename given!\n");
    exit(1);
  }
  generate_file(filename.unwrap.data());
  auto buffer = pft::read_file_as_string_view(filename.unwrap.data());
  if (!buffer.has_value) {
    fprintf(stderr, "Could not read file\n");
    exit(1);
  }

  auto vec = split_by(buffer.unwrap, '\n');
  pft::ignore_header_lines(vec, 4);

  auto buf = pft::as_floats(vec);

  TFile hfile("file.root", "RECREATE");
  auto myTree = createTree(buf);
  myTree->Write();
  hfile.Close();

  return 0;
}
