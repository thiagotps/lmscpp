#include <bits/stdc++.h>
#include <symengine/expression.h>
#include <symengine/functions.h>
#include <symengine/symbol.h>
#include <symengine/ntheory.h>

#include "stochastic.hpp"
#include "utils.hpp"

using namespace SymEngine;
using namespace std;
using namespace stochastic;

bool eq(const vec_basic &a, const vec_basic &b){
  if (a.size() != b.size())
    return false;

  for (size_t i = 0; i < a.size(); i++)
    if (not eq(*a[i],*b[i]))
      return false;

  return true;
}

bool even(const Integer &i){
  static RCP<const Integer> two{integer(2)};
  return mod(i,*two)->is_zero();
}

RCP<const Integer> operator""_i(unsigned long long i) {return integer(i);}

// TODO: Currently we only test the case with one variable

void test_StochasticProcess(){
  RCP<const Basic> b{symbol("γ")};
  StochasticProcess u{"u", [b](const vec_basic& args, const RCP<const Integer> &n) -> RCP<const Basic> {
                             if (even(*n))
                               return zero;

                             return pow(b,n);
                           }};

  auto inst = u(symbol("x"));
  assert(StochasticProcess::is_random(inst));

  auto expr1 = pow(inst,10_i);
  auto expr2 = pow(inst,23_i);

  assert(eq(*StochasticProcess::moment(expr1), *zero) == true);
  assert(eq(*StochasticProcess::moment(expr2), *pow(b,23_i)) == true);

  auto v = make_rcp<const FunctionSymbol>("ν", symbol("k"));
  assert(StochasticProcess::is_random(v) == false);

  cout << "StochasticProcess [OK]" << endl;
}


void test_get_multuple(){
  auto x = symbol("x");
  auto y = symbol("y");
  auto z = symbol("z");

  auto two = 2_i;
  auto expr1 = mul({x,y,z});
  auto expr2 = mul({two,x,y,z});

  assert(eq(get_multuple(expr1), {one, x, y, z}));
  assert(eq(get_multuple(expr2) , {two,x,y,z}));
  assert(eq(get_multuple(x) , {one, x}));

  cout << "test_get_multuple [OK]" << endl;
}

void test_get_addtuple(){
  auto x = symbol("x");
  auto y = symbol("y");
  auto z = symbol("z");

  auto two = 2_i;
  auto expr1 = add({x,y,z});
  auto expr2 = add({two,x,y,z});


  // NOTE: For some reason get_args() returns {z,y,x} instead of {x,y,z}
  // If there is some change in symengine affecting this ordering, this test will fail.
  assert(eq(get_addtuple(expr1), {zero, z, y, x}));
  assert(eq(get_addtuple(expr2) , {two, z, y, x}));
  assert(eq(get_addtuple(x) , {zero, x}));

  cout << "test_get_addtuple [OK]" << endl;
}

void test_split_independents(){
  ExpectedOperator E{[](const FunctionSymbol& x, const FunctionSymbol& y){
                       auto xy = [](const FunctionSymbol &x, const FunctionSymbol &y)
                                 {
                                   if (x.get_name() == "v")
                                     return true;
                                   if (x.get_name() == "u" and y.get_name() == "u")
                                     return not eq(x,y);

                                   return false;
                                 };
                       return xy(x,y) or xy(y,x);
                     }};
  StochasticProcess u{"u", [](const vec_basic&, const RCP<const Integer>& n){
                             static const auto gamma{symbol("γ")};
                             return pow(gamma, n);
                           }};
  StochasticProcess v{"v", [](const vec_basic&, const RCP<const Integer>&){
                             return zero;
                           }};
  auto k = symbol("k");
  auto uk{u(k)}, vk{v(k)};
  RCP<const Basic> expr1 = vk;
  expr1 = add(expr1,pow(uk,2_i));
  expr1  = add(expr1, uk);
  auto e1  = rcp_dynamic_cast<const FunctionSymbol>(E(expr1));

  cout << *e1 << endl;
  cout << *E.expand(e1) << endl;
}

int main() {

  test_StochasticProcess();
  test_get_multuple();
  test_get_addtuple();
  test_split_independents();
  return 0;
}
