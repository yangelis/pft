#include "../pft.hpp"
#include "../utils.hpp"

int main()
{
  ////////////////////////////////////////////////
  // Maybe
  ////////////////////////////////////////////////
  {
    constexpr auto m = pft::Some<i32>(123459);
    static_assert(m.unwrap == 123459 && m.has_value);

    pft::println(stdout, "[PASSED] ", "Maybe tests");
  }

  ////////////////////////////////////////////////
  // Matrix
  ////////////////////////////////////////////////
  {
    const f64 in[][3] = {
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

    pft::println(stdout, "[PASSED] ", "Matrix tests");
  }

  ////////////////////////////////////////////////
  // filter
  ////////////////////////////////////////////////
  {
    const auto xx       = pft::arange(0, 10);
    const auto filtered = pft::filter([](auto x) { return x > 5; }, xx);
    for (std::size_t i = 0; i < filtered.size(); ++i) {
      assert(filtered[i] == static_cast<i32>(i) + 6);
    }

    auto xx2 = pft::arange(-2, 10000000);
    const auto filtered2 =
        pft::filter([](auto x) { return x > 0; }, std::move(xx2));
    // const auto filtered2 =
    //     pft::filter([](auto x) { return x > 0; }, xx2);
    //    pft::println(stdout, xx2);
    pft::println(stdout, "[PASSED] ", "filter tests");
  }

  ////////////////////////////////////////////////
  // map
  ////////////////////////////////////////////////
  {
    auto to_int    = [](pft::StringView x) { return std::stoi(x.as_string()); };
    auto to_double = [](pft::StringView x) { return std::stod(x.as_string()); };

    const std::vector<pft::StringView> num_test = {"123.456", "9999.9",
                                                   "69.420"};
    const std::vector<f64> nums_double          = {123.456, 9999.9, 69.420};
    const std::vector<i32> nums_ints            = {123, 9999, 69};

    assert(nums_double == pft::map(to_double, num_test));
    assert(nums_ints == pft::map(to_int, num_test));

    const auto ones = std::vector<i32>(10, 1);
    const auto twos = std::vector<i32>(10, 2);
    auto d = pft::map([](const auto&... args) { return (args + ...); }, ones,
                      ones, ones, twos);

    for (std::size_t i = 0; i < d.size(); ++i) {
      assert(d[i] == 3 * ones[i] + twos[i]);
    }

    d = pft::map([](const auto&... args) { return (args + ... + 1); }, ones);

    for (std::size_t i = 0; i < d.size(); ++i) {
      assert(d[i] == ones[i] + 1);
    }

    pft::println(stdout, "[PASSED] ", "map tests");
  }

  ////////////////////////////////////////////////
  // StringView
  ////////////////////////////////////////////////
  {
    pft::StringView text = {"  just a nice line with some spaces \n \t"};
    assert(pft::triml(text) == "just a nice line with some spaces \n \t");
    assert(pft::trimr(text) == "just a nice line with some spaces");
    assert(text.chop(4) == " a nice line with some spaces");

    pft::println(stdout, "[PASSED] ", "StringView tests");
  }

  ////////////////////////////////////////////////
  // factorial
  ////////////////////////////////////////////////
  {
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
  }

  ////////////////////////////////////////////////
  // arange
  ////////////////////////////////////////////////
  {
    const auto big_vec = pft::arange(0, 1000000);
    std::vector<i32> vec_iota(big_vec.size());
    std::iota(std::begin(vec_iota), std::end(vec_iota), 0);
    for (i32 i = 0; i < static_cast<i32>(big_vec.size()); ++i) {
      assert(big_vec[i] == i);
    }
    pft::println(stdout, "[PASSED] ", "arange tests");
  }

  ////////////////////////////////////////////////
  // zip_to_pair, zip_with
  ////////////////////////////////////////////////
  {
    const auto xx          = pft::arange(0, 10);
    const auto xx_reversed = pft::arange(9, -1, -1);
    for (const auto& ii : pft::zip_to_pair(xx, xx_reversed)) {
      assert(ii.first + ii.second == 9);
    }

    for (const auto& ii : pft::zip_to_pair(xx, xx)) {
      assert(ii.first == ii.second);
    }

    for (const auto& [i1, i2] : pft::zip_to_pair(xx, xx)) {
      assert(i1 == i2);
    }

    for (const auto &x :
         pft::zip_with([](const auto&a, const auto&b) { return a + b; }, xx,
                       xx_reversed)) {
      assert(x == 9);
    }
    pft::println(stdout, "[PASSED] ", "zip tests");
  }

  ////////////////////////////////////////////////
  // pad_right, pad_left, pad
  ////////////////////////////////////////////////
  {
    const auto xx             = pft::arange(0, 10);
    const auto paded_xx_right = pft::pad_right(xx, 10, 69);
    for (size_t i = 0; i < xx.size(); ++i) {
      assert(paded_xx_right[i] == xx[i]);
    }
    for (size_t i = xx.size(); i < 10 + xx.size(); ++i) {
      assert(paded_xx_right[i] == 69);
    }

    const auto paded_xx_left = pft::pad_left(xx, 10, 69);
    for (size_t i = 0; i < 10; ++i) {
      assert(paded_xx_left[i] == 69);
    }
    for (size_t i = 0; i < xx.size(); ++i) {
      assert(xx[0] == paded_xx_left[0 + 10]);
    }

    const auto paded_xx = pft::pad(xx, 10, 69);
    for (size_t i = 0; i < 10; ++i) {
      assert(paded_xx[i] == 69);
      assert(paded_xx[i + 10 + xx.size()] == 69);
    }

    const auto paded_xx_until = pft::pad_right_until(xx, pft::Slice{1, 5}, 33);
    assert(paded_xx_until.size() == 33);
    for (size_t i = 0; i < xx.size(); ++i) {
      assert(paded_xx_until[i] == xx[i]);
    }

    pft::println(stdout, "[PASSED] ", "pad tests");
  }

  ////////////////////////////////////////////////
  // chunks
  ////////////////////////////////////////////////
  {
    const auto vec_to_be_chunked = pft::arange(0, 8);
    const auto chunked_vec       = pft::chunks(3, vec_to_be_chunked);
    const auto check_chunk1      = chunked_vec[0] == std::vector<i32>{0, 1, 2};
    const auto check_chunk2      = chunked_vec[1] == std::vector<i32>{3, 4, 5};
    const auto check_chunk3      = chunked_vec[2] == std::vector<i32>{6, 7};
    assert(check_chunk1);
    assert(check_chunk2);
    assert(check_chunk3);

    pft::println(stdout, "[PASSED] ", "chunking tests");
  }

  ////////////////////////////////////////////////
  // adjacent_*
  ////////////////////////////////////////////////
  {
    const std::vector<i32> v1{0, 1, 2, 3, 40, 40, 41, 41, 5};
    std::vector<i32> res;
    std::vector<i32> res2;
    std::adjacent_difference(std::begin(v1), std::end(v1),
                             std::back_inserter(res), std::plus{});
    pft::adjacent_transformN<2>(std::begin(v1), std::end(v1),
                                std::back_inserter(res2), std::plus{});

    for (size_t i = 0; i < v1.size(); ++i) {
      assert(res[i] == res[i]);
    }
    pft::println(stdout, "[PASSED] ", "adjacent_* tests");
  }

  return 0;
}
