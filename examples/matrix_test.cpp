#include "../pft.hpp"
#include "../utils.hpp"
#include <iostream>
#include <random>

using namespace std;

int main() {
  constexpr size_t N = 3;
  constexpr size_t M = 3;

  random_device rd;
  mt19937 gen(rd());
  uniform_int_distribution<> dis(-1, 1);

  pft::Matrix<i32> matrix1(N, M);

  for (size_t i = 0; i < N; ++i) {
    for (size_t j = 0; j < M; ++j) {
      matrix1(i, j) = dis(gen);
    }
  }
  pft::println(stdout, "Random matrix\n", matrix1);

  auto mT = matrix1.transpose();
  pft::println(stdout, "Transpose matrix\n", mT);

  pft::println(stdout, matrix1.mult(mT));

  constexpr size_t N2 = 10;
  constexpr size_t M2 = 10;
  pft::Matrix<i64> matrix2(N2, M2);
  matrix2.diagonal();

  pft::println(stdout, "Diagonal matrix\n", matrix2);

  double in[][3] = {
      {12, -51, 4}, {6, 167, -68}, {-4, 24, -41}, {-1, 1, 0}, {2, 0, 3},
  };

  pft::Matrix<f64> mat_from_array(in);

  pft::println(stdout, "From 2d array\n", mat_from_array);

  auto [Q, R] = QRDecomposition(mat_from_array);

  pft::println(stdout, "Q:\n", Q);
  pft::println(stdout, "R:\n", R);
  pft::println(stdout, "A=Q*R:\n", Q.mult(R));

  double arr2[][2] = {{4, 3}, {6, 3}};

  pft::Matrix<f64> mat_for_lu(arr2);
  LUdecomposition lumat(mat_for_lu);

  pft::println(stdout, lumat.lu);

  return 0;
}
