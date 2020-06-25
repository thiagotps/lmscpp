#include <bits/stdc++.h>
#include <symengine/expression.h>
#include <symengine/functions.h>
#include <symengine/symbol.h>
#include <symengine/ntheory.h>

#include "stochastic.hpp"
#include "utils.hpp"
#include "operators.hpp"

using namespace SymEngine;
using namespace SymEngine::OverloadedOperators;
using namespace std;
using namespace stochastic;




//    TODO: Currently we only test the case with one variable
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


void test_ExptectOperator_and_EquationSet(){
  // u is independent of v but v is not independent of w. w and u are independent.

  const ExpectedOperator E{[](const FunctionSymbol& x, const FunctionSymbol& y){
                       auto xy = [](const FunctionSymbol &x, const FunctionSymbol &y)
                                 {
                                   auto xname{x.get_name()};
                                   auto yname{y.get_name()};

                                   if (xname == "u" and yname == "v")
                                     return true;

                                   if (xname == "w" and yname == "u")
                                     return true;

                                   return false;
                                 };
                       return xy(x,y) or xy(y,x);
                     }};

  const auto gamma{symbol("γ")}, sigma{symbol("σ")}, alpha{symbol("α")};
  StochasticProcess u{"u", [&gamma](const vec_basic&, const RCP<const Integer>& n){
                             return pow(gamma, n);
                           }};
  StochasticProcess v{"v", [](const vec_basic&, const RCP<const Integer>&){
                             return zero;
                           }};
  StochasticProcess w{"w", [&](const vec_basic&, const RCP<const Integer>& n) -> RCP<const Basic> {
                             if (even(*n))
                               return zero;

                             return pow(sigma, n);
                           }};
  auto k = symbol("k");
  auto uk{u(k)}, vk{v(k)}, wk(w(k));

  // Testing if E[u(k) + u(k)**2 + v(k) + alpha] = gamma + gamma**2 + alpha
  RCP<const Basic> expr1 = add(vk, alpha);
  expr1 = add(expr1,pow(uk,2_i));
  expr1  = add(expr1, uk);
  auto e1  = rcp_dynamic_cast<const FunctionSymbol>(E(expr1));

  RCP<const Basic> r1{pow(gamma,2_i)};
  r1 = add(r1, gamma);
  r1 = add(r1, alpha);
  assert(eq(*E.expand(e1),*r1));

  // Testing if E[alpha*u(k)**2v(k)w(k)] = alpha*E[u(k)**2]E[v(k)w(k)] = alpha*gammma**2*E[v(k)*w(k)]
  RCP<const Basic> expr2{mul({pow(uk,2_i), vk, wk, alpha})};
  auto e2 = rcp_dynamic_cast<const FunctionSymbol>(E(expr2));
  auto r2{mul(pow(gamma,2_i),alpha)};
  r2 = mul(r2, E(mul(vk,wk)));
  assert(eq(*E.expand(e2), *r2));
  cout << "test_ExptectOperator [OK]" << endl;

  auto mulexpr = mul({pow(u(k), 4_i),
                      pow(u(add(k,1_i)), 2_i),
                      pow(v(add(k, 2_i)), 2_i)
    });

  auto mulexpr_res = mul({pow(u(sub(k, 2_i)), 4_i),
                      pow(u(sub(k,one)), 2_i),
                      pow(v(k), 2_i)
    });

  auto gcvf = get_cnt_var_func(rcp_dynamic_cast<const FunctionSymbol>(E(mulexpr)));
  assert(eq(*gcvf.first, *2_i));
  assert(eq(*gcvf.second, *k));

  cout << "test get_cnt_var_func [OK]" << endl;

  RCP<const Basic> lhsb = pow(u(k), 2_i) * pow(u(k - one), 4_i) * pow(w(k + one), 2_i);
  RCP<const Basic> lhsb2 = pow(u(k), 4_i) * pow(u(k + one), 2_i) * pow(w(k + 2_i), 2_i);
  auto lhs = rcp_dynamic_cast<const FunctionSymbol>(E(lhsb));
  auto lhs2 = rcp_dynamic_cast<const FunctionSymbol>(E(lhsb2));

  RCP<const Basic> rhs = pow(gamma, 2_i)*E(u(k)*w(k)) + E(pow(u(k - 1_i),2_i)*w(k));

  EquationSet eqs(k);
  eqs.setitem(lhs, rhs);

  auto res = pow(gamma, 2_i)*E(u(k + 1_i)*w(k + 1_i)) + E(pow(u(k),2_i)*w(k + 1_i));
  assert(eq(*eqs.getitem(lhs2), *res));

  cout << "EquationSet [OK]" << endl;

  auto t1 = E(v(k) * w(k));
  auto t2 = E(pow(v(k), 2_i) * w(k));
  auto t3 = E(pow(v(k), 3_i) * w(k));

  auto stv = states_vars(alpha*t1 + 2_i*pow(sigma, 10_i)*t2*t3);
  assert(eq(*stv[0], *t3));
  assert(eq(*stv[1], *t2));
  assert(eq(*stv[2], *t1));

  cout << "state_vars [OK]" << endl;
}

const auto SYMALPHA{symbol("α")}, SYMBETA{symbol("β")};

void test_expand()
{
  auto x = [](const auto & k) {return make_rcp<const FunctionSymbol>("x", k); };
  auto y = [](const auto & k) {return make_rcp<const FunctionSymbol>("y", k); };

  const auto k{symbol("k")};

  EquationSet eqs(k);


  // x(k + 1) = αx(k) + 10
  eqs.setitem(x(k + 1_i) , SYMALPHA * x(k) + 10_i);
  // y(k + 1) = 2βy(k)**2
  eqs.setitem(y(k + 1_i), 2_i * SYMBETA * pow(y(k), 2_i));
  // 2*x(k)*y(k+2)**2
  auto expr  = rcp_dynamic_cast<const Mul>(2_i * x(k) * pow(y(k + 2_i), 2_i));
  // The expected result is 8β**2 * y(k + 1)**4 * ( αx(k - 1) + 10)
  auto expected_res = 8_i * pow(SYMBETA, 2_i) * pow(y(k + 1_i), 4_i) * (SYMALPHA * x(k - 1_i) + 10_i);
  // What we really get
  auto res = expand(expr, eqs);
  assert(eq(*res, *expected_res));

  cout << "test_expand[OK]" << endl;
}


int main() {

  test_StochasticProcess();
  test_get_multuple();
  test_get_addtuple();
  test_ExptectOperator_and_EquationSet();
  test_expand();
  return 0;
}
