#include <bits/stdc++.h>
#include <symengine/expression.h>
#include <symengine/functions.h>
#include <symengine/symbol.h>
#include <symengine/add.h>
#include <symengine/mul.h>
#include <symengine/subs.h>

using namespace SymEngine;
using namespace std;

RCP<const Basic>
operator^(const RCP<const Basic> &l, const RCP<const Basic> &r)
{
  return pow(l,r);
}

RCP<const Basic>
operator+(const RCP<const Basic> &l, const RCP<const Basic> &r)
{
  return add(l,r);
}

RCP<const Basic>
operator-(const RCP<const Basic> &l, const RCP<const Basic> &r)
{
  return sub(l,r);
}

RCP<const Basic>
operator*(const RCP<const Basic> &l, const RCP<const Basic> &r)
{
  return mul(l,r);
}


int main() {
  auto x = symbol("x");
  auto y = symbol("y");
  auto f = function_symbol("f", x);
  auto g = function_symbol("f", add(x, one));

  auto expr = mul(f,g);

  auto expr2 = xreplace(expr,{{x, add(y, one)}});

  cout << *expr2 << endl;

  cout << *(x ^ integer(2)) << endl;

  return 0;
}
