#include "../pft.hpp"
#include "../utils.hpp"
#include <numeric>

int main() {
  ////////////////////////////////////////////////
  // Maybe
  ////////////////////////////////////////////////
  constexpr auto m = pft::Some<i32>(123459);
  static_assert(m.unwrap == 123459 && m.has_value);

  pft::println(stdout, "[PASSED] ", "Maybe tests");

  ////////////////////////////////////////////////
  // Matrix
  ////////////////////////////////////////////////

  double in[][3] = {
      {12, -51, 4}, {6, 167, -68}, {-4, 24, -41}, {-1, 1, 0}, {2, 0, 3},
  };

  pft::Matrix<f64> mat_from_array(in);

  assert(mat_from_array(0, 0) == 12);
  assert(mat_from_array(0, 1) == -51);
  assert(mat_from_array(0, 2) == 4);
  assert(mat_from_array(1, 0) == 6);
  assert(mat_from_array(1, 1) == 167);
  assert(mat_from_array(1, 2) == -68);
  assert(mat_from_array(2, 0) == -4);
  assert(mat_from_array(2, 1) == 24);
  assert(mat_from_array(2, 2) == -41);
  assert(mat_from_array(3, 0) == -1);
  assert(mat_from_array(3, 1) == 1);
  assert(mat_from_array(3, 2) == 0);
  assert(mat_from_array(4, 0) == 2);
  assert(mat_from_array(4, 1) == 0);
  assert(mat_from_array(4, 2) == 3);
  ////////////////////////////////////////////////
  // filter
  ////////////////////////////////////////////////
  auto xx       = pft::arange(0, 10);
  auto filtered = pft::filter([](const auto x) { return x > 5; }, xx);
  for (std::size_t i = 0; i < filtered.size(); ++i) {
    assert(filtered[i] == (i32)i + 6);
  }

  pft::println(stdout, "[PASSED] ", "filter tests");

  ////////////////////////////////////////////////
  // map
  ////////////////////////////////////////////////
  auto to_int    = [](const pft::StringView x) { return stoi(std::string(x)); };
  auto to_double = [](const pft::StringView x) { return stod(std::string(x)); };

  std::vector<pft::StringView> num_test = {"123.456", "9999.9", "69.420"};
  std::vector<f64> nums_double          = {123.456, 9999.9, 69.420};
  std::vector<i32> nums_ints            = {123, 9999, 69};

  assert(nums_double == pft::map(to_double, num_test));
  assert(nums_ints == pft::map(to_int, num_test));

  auto ones = std::vector<i32>(10, 1);
  auto twos = std::vector<i32>(10, 2);
  auto d    = pft::map([](const auto&... args) { return (args + ...); }, ones,
                    ones, ones, twos);

  for (std::size_t i = 0; i < d.size(); ++i) {
    assert(d[i] == 3 * ones[i] + twos[i]);
  }

  d = pft::map([](const auto&... args) { return (args + ... + 1); }, ones);

  for (std::size_t i = 0; i < d.size(); ++i) {
    assert(d[i] == ones[i] + 1);
  }

  pft::println(stdout, "[PASSED] ", "map tests");

  ////////////////////////////////////////////////
  // StringView
  ////////////////////////////////////////////////
  pft::StringView text = {"  just a nice line with some spaces \n \t"};
  assert(pft::triml(text) == "just a nice line with some spaces \n \t");
  assert(pft::trimr(text) == "just a nice line with some spaces");
  assert(text.chop(4) == " a nice line with some spaces");

  pft::println(stdout, "[PASSED] ", "StringView tests");

  ////////////////////////////////////////////////
  // factorial
  ////////////////////////////////////////////////
  assert(utils::factorial(10) == 3628800.0);
  assert(utils::factorial(11) == 39916800.0);
  assert(utils::factorial(12) == 479001600.0);
  assert(utils::factorial(13) == 6227020800.0);
  assert(utils::factorial(14) == 87178291200.0);
  assert(utils::factorial(15) == 1307674368000.0);
  assert(utils::factorial(16) == 20922789888000.0);
  assert(utils::factorial(17) == 355687428096000.0);
  assert(utils::factorial(18) == 6402373705728000.0);
  assert(utils::factorial(19) == 121645100408832000.0);
  assert(utils::factorial(20) == 2432902008176640000.0);
  assert(utils::factorial(21) == 51090942171709440000.0);
  pft::println(stdout, "[PASSED] ", "factorial tests");

  ////////////////////////////////////////////////
  // arange
  ////////////////////////////////////////////////
  auto big_vec = pft::arange(0, 1000000);
  std::vector<i32> vec_iota(big_vec.size());
  std::iota(vec_iota.begin(), vec_iota.end(), 0);
  for (i32 i = 0; i < static_cast<i32>(big_vec.size()); ++i) {
    assert(big_vec[i] == i);
  }
  pft::println(stdout, "[PASSED] ", "arange tests");

  ////////////////////////////////////////////////
  // zip_to_pair, zip_with
  ////////////////////////////////////////////////
  auto xx_reversed = pft::arange(9, -1, -1);
  for (const auto& ii : pft::zip_to_pair(xx, xx_reversed)) {
    assert(ii.first + ii.second == 9);
  }

  for (const auto& ii : pft::zip_to_pair(xx, xx)) {
    assert(ii.first == ii.second);
  }

  for (const auto& [i1, i2] : pft::zip_to_pair(xx, xx)) {
    assert(i1 == i2);
  }

  for (const auto &x : pft::zip_with(
           [](const auto&a, const auto&b) { return a + b; }, xx, xx_reversed)) {
    assert(x == 9);
  }
  pft::println(stdout, "[PASSED] ", "zip tests");

  ////////////////////////////////////////////////
  // pad_right, pad_left, pad
  ////////////////////////////////////////////////
  auto paded_xx_right = pft::pad_right(xx, 10, 69);
  for (size_t i = xx.size(); i < 10 + xx.size(); ++i) {
    assert(paded_xx_right[i] == 69);
  }

  auto paded_xx_left = pft::pad_left(xx, 10, 69);
  for (size_t i = 0; i < 10; ++i) {
    assert(paded_xx_left[i] == 69);
  }

  auto paded_xx = pft::pad(xx, 10, 69);
  for (size_t i = 0; i < 10; ++i) {
    assert(paded_xx[i] == 69);
    assert(paded_xx[i + 10 + xx.size()] == 69);
  }

  pft::println(stdout, "[PASSED] ", "pad tests");

  ////////////////////////////////////////////////
  // chuks
  ////////////////////////////////////////////////
  return 0;
}
