#include <bitsery/bitsery.h>
#include <bitsery/adapter/buffer.h>
#include <bitsery/brief_syntax.h>
#include <bitsery/brief_syntax/vector.h>
#include <bitsery/brief_syntax/string.h>
#include <bitsery/brief_syntax/memory.h>
#include <bits/stdc++.h>

using namespace bitsery;
using namespace std;

using buffer = vector<uint8_t>;
using output_adapter = OutputBufferAdapter<buffer>;
using input_adapter = InputBufferAdapter<buffer>;


class bank{
  vector<string> names{};
  vector<uint32_t> ids{};
public:
  bank() = default;
  bank(initializer_list<string> n, initializer_list<uint32_t> i): names{n}, ids{i} {}
  bool operator==(const bank& b) const {
    return names == b.names and ids == b.ids;
  }

  template<typename S>
  void serialize(S& s){
    s(names, ids);
  }

};

class contry{
  string contryname{};
  bank centralbank{};
public:
  contry() = default;
  contry(string&& name, bank& cb): contryname{name}, centralbank{cb} {}
  bool operator==(const contry& b) const {
    return contryname == b.contryname and centralbank == b.centralbank;
  }

  template<typename S>
  void serialize(S &s){
    s(contryname, centralbank);
  }
};

int main() {
  buffer b;
  bank ba{{"Thiago", "Diego"}, {1, 2}};
  contry ca{"brasil", ba}, c2{};

  auto written = quickSerialization(output_adapter{b}, ca);
  cout << "Serialized " << written << " bytes." << endl;
  auto state = quickDeserialization(input_adapter{b.begin(), written}, c2);


  assert(state.first == ReaderError::NoError and state.second);
  assert(ca == c2);

  return 0;
}

//    class teste{
//   vector<string> v{};
//   uint32_t level{};
// public:
//   teste() = default;
//   teste(initializer_list<string> l, uint32_t lv): v{l}, level{lv} {}
//   bool operator==(const teste &lhs) const {
//     return v == lhs.v and level == lhs.level;
//   }
//   template<typename S>
//   void serialize(S &s){
//     s(v, level);
//   }
// };

// class teste2{
//   vector<teste> v{};
//   string teste_name{};
// public:
//   teste2() = default;
//   teste2(initializer_list<teste> lt, string s): v{lt}, teste_name{s} {}
//   bool operator==(const teste2 &other) const {
//     return v == other.v and teste_name == other.teste_name;
//   }

//   template<typename S>
//   void serialize(S &s){
//     s(v, teste_name);
//   }
// };

// int main(){
//   buffer b;
//   teste2 t{{{{"Hello", "World"}, 10},{{"Ola", "Mundo"}, 20}}, "João"};

//   auto written = quickSerialization(output_adapter{b}, t);
//   cout << "Serialized " << written << " bytes." << endl;

//   teste2 t2{};
//   auto state = quickDeserialization(input_adapter{b.begin(), written}, t2);

//   assert(state.first == ReaderError::NoError and state.second);
//   assert(t == t2);
// }
