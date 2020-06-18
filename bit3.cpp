#include <bitsery/bitsery.h>
#include <bitsery/traits/vector.h>
#include <bitsery/traits/string.h>
#include <bitsery/adapter/buffer.h>
#include <bitsery/ext/pointer.h>
#include <bitsery/ext/inheritance.h>
#include <bitsery/ext/std_smart_ptr.h>
#include <bits/stdc++.h>

using namespace std;
using namespace bitsery;
using bitsery::ext::StdSmartPtr;
using Buffer = vector<uint8_t>;
using Writer = OutputBufferAdapter<Buffer>;
using Reader = InputBufferAdapter<Buffer>;
// using TContext = std::tuple<ext::PointerLinkingContext>;
using TContext = ext::PointerLinkingContext;
using MySerializer = Serializer<Writer, TContext>;
using MyDeserializer = Deserializer<Reader, TContext>;

struct teste{
  string x{};
  bool operator==(const teste &o) const {
    return x == o.x;
  }

  template<typename S>
  void serialize(S &s){
    // s.value4b(x);
    s.text1b(x,100);
  }
};

template<typename S>
void serialize(S &s, unique_ptr<teste> &p) {
  s.ext(p, StdSmartPtr{});
};


int main() {
  unique_ptr<teste> p1{new teste{"Ola Mundo"}}, p2{};
  Buffer buffer{};


  //  TContext ctx{};
  // MySerializer ser{ctx, buffer};
  // ser.object(p1);
  // ser.adapter().flush();
  // auto written = ser.adapter().writtenBytesCount();

  // MyDeserializer des{ctx, buffer.begin(), written};
  // des.object(p2);
  // assert(*p1 == *p2);

  TContext ctx{};
  auto written = quickSerialization(ctx, Writer{buffer}, p1);
  auto state = quickDeserialization(ctx, Reader{buffer.begin(), written}, p2);

  assert(*p1 == *p2);

  return 0;
}

