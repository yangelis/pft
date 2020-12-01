#define PFT_IMPLEMENTATION
#include "../pft.hpp"
#include <iostream>
#include <random>

using namespace std;

int main() {
  constexpr size_t N = 15;
  constexpr size_t M = 15;

  random_device rd;
  mt19937 gen(rd());
  uniform_int_distribution<> dis(-1, 1);

  pft::Matrix<i32> matrix1(N, M);

  for (size_t i = 0; i < N; ++i) {
    for (size_t j = 0; j < M; ++j) {
      matrix1(i, j) = dis(gen);
    }
  }
  pft::println(stdout, "Transpose matrix\n", matrix1);

  auto mT = matrix1.transpose();
  pft::println(stdout, "Transpose matrix\n", mT);

  constexpr size_t N2 = 10;
  constexpr size_t M2 = 10;
  pft::Matrix<i64> matrix2(N2, M2);
  matrix2.diagonal();

  pft::println(stdout, "Diagonal matrix\n", matrix2);

  return 0;
}
