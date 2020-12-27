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

  assert(factorial(10) == 3628800);
  assert(factorial(11) == 39916800);
  assert(factorial(12) == 479001600);
  assert(factorial(13L) == 6227020800);
  assert(factorial(14L) == 87178291200);
  assert(factorial(15L) == 1307674368000);
  assert(factorial(16L) == 20922789888000);
  assert(factorial(17L) == 355687428096000);
  assert(factorial(18L) == 6402373705728000);
  assert(factorial(19L) == 121645100408832000);
  assert(factorial(20L) == 2432902008176640000);

  return 0;
}
