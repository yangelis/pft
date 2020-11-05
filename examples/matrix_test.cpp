// #define PFT_USE_ROOT
#include "../pft.hpp"
#include <iostream>
#include <random>

// #include <TUnixSystem.h>

using namespace std;

int main() {
  // gSystem->IgnoreSignal(kSigSegmentationViolation, true);
  constexpr size_t N = 5000;
  constexpr size_t M = 5000;

  random_device rd;
  mt19937 gen(rd());
  uniform_real_distribution<> dis(-1.0, 1.0);


  pft::Matrix2v<float> matrix2(N, M);

  for (size_t i = 0; i < N; ++i) {
    for (size_t j = 0; j < M; ++j) {
      matrix2(i, j) = dis(gen);
    }
  }

  return 0;
}
