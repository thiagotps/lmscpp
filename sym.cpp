#include <bits/stdc++.h>
#include <symengine/expression.h>
#include <symengine/functions.h>
#include <symengine/symbol.h>
#include <symengine/add.h>
#include <symengine/mul.h>

using namespace SymEngine;
using namespace std;




int main() {
  auto x{symbol("x")}, y{symbol("y")};
  auto expr1{add(x,y)};
  expr1 = mul(expr1, expr1);

  cout << *expr1 << endl;
  cout << *expand(expr1) << endl;


  return 0;
}
