#include <bits/stdc++.h>

#include <symengine/expression.h>
#include <symengine/functions.h>
#include <symengine/symbol.h>
#include <symengine/add.h>
#include <symengine/mul.h>
#include <symengine/subs.h>
#include <symengine/logic.h>

#include <bitsery/bitsery.h>
#include <bitsery/adapter/buffer.h>
#include <bitsery/brief_syntax.h>
#include <bitsery/brief_syntax/vector.h>
#include <bitsery/brief_syntax/string.h>

#include <lmscpp/operators.hpp>
#include <lmscpp/utils.hpp>
#include <lmscpp/nodevisitor.hpp>

using namespace stochastic;
using namespace SymEngine::OverloadedOperators;
using namespace SymEngine;
using namespace NodeVisitor;
using namespace bitsery;
using namespace std;

using buffer = vector<uint8_t>;
using output_adapter = OutputBufferAdapter<buffer>;
using input_adapter = InputBufferAdapter<buffer>;

int main() {
  const auto gamma{symbol("γ")}, alpha{symbol("α")}, k{symbol("k")};
  auto fk{function_symbol("f", k)};
  auto a = [](const auto & n) {return function_symbol("a", n);};
  RCP<const Basic> expr1{gamma*a(zero)*a(one) + pow(fk, 2_i) * alpha};

  Node h1;
  expr1->accept(h1);


  buffer b;
  auto written = quickSerialization(output_adapter{b}, h1);

  Node h2;
  auto state = quickDeserialization(input_adapter{b.begin(), written}, h2);

  assert(state.first == ReaderError::NoError and state.second);
  assert(h1 == h2);

  auto expr2{h2.to_basic()};

  assert(eq(*expr1,*expr2));
  return 0;
}
