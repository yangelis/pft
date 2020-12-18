#include "../pft.hpp"

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

  return 0;
}
