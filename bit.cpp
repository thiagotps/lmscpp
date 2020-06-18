#include <bits/stdint-uintn.h>
#include <bitsery/bitsery.h>
#include <bitsery/adapter/buffer.h>
#include <bitsery/brief_syntax.h>
#include <bitsery/brief_syntax/vector.h>
#include <bitsery/brief_syntax/string.h>
#include <initializer_list>
#include <iostream>

using namespace bitsery;
using namespace std;
using namespace string_literals;

using buffer = vector<uint8_t>;
using output_adapter = OutputBufferAdapter<buffer>;
using input_adapter = InputBufferAdapter<buffer>;


class my_struct{
  uint32_t i;
  string str;
  vector<float> fs;
public:
  my_struct(): i{}, str{}, fs{{}} {}
  my_struct(uint32_t a, string s, initializer_list<float> l): i{a}, str{s}, fs{l} {}
  bool operator==(const my_struct& m){
    return i == m.i and str == m.str and fs == m.fs;
  }
  template<typename S>
  void serialize(S& s){
    s(i, str, fs);
  }
};


int main() {
  my_struct a{1,"Hello"s, {1.0,2.0,3.0}};
  buffer b;
  auto written = quickSerialization(output_adapter{b}, a);

  my_struct m{};
  auto state = quickDeserialization(input_adapter{b.begin(), written}, m);

  assert(state.first == ReaderError::NoError and state.second);
  assert(a == m);
 

  return 0;
}
