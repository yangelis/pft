#include "../pft.hpp"
#include "../utils.hpp"

using namespace std;

int main(int argc, char* argv[]) {
  auto to_int    = [](const pft::StringView x) { return stoi(string(x)); };
  auto to_double = [](const pft::StringView x) { return stod(string(x)); };

  vector<pft::StringView> num_test = {"123.456", "9999.9", "69.420"};
  vector<double> nums_double       = {123.456, 9999.9, 69.420};
  vector<int> nums_ints            = {123, 9999, 69};

  assert(nums_double == pft::map(to_double, num_test));
  assert(nums_ints == pft::map(to_int, num_test));

  pft::StringView text = {"  just a nice line with some spaces \n \t"};
  assert(pft::triml(text) == "just a nice line with some spaces \n \t");
  assert(pft::trimr(text) == "just a nice line with some spaces");
  assert(text.chop(4) == " a nice line with some spaces");

  assert(factorial(10) == 3628800.0);
  assert(factorial(11) == 39916800.0);
  assert(factorial(12) == 479001600.0);
  assert(factorial(13) == 6227020800.0);
  assert(factorial(14) == 87178291200.0);
  assert(factorial(15) == 1307674368000.0);
  assert(factorial(16) == 20922789888000.0);
  assert(factorial(17) == 355687428096000.0);
  assert(factorial(18) == 6402373705728000.0);
  assert(factorial(19) == 121645100408832000.0);
  assert(factorial(20) == 2432902008176640000.0);
  assert(factorial(21) == 51090942171709440000.0);
  return 0;
}
