#ifndef __UTILS_H_
#define __UTILS_H_

#include <symengine/expression.h>
#include <symengine/functions.h>
#include <symengine/symbol.h>
#include <tuple>

namespace stochastic{
  using namespace SymEngine;
  using namespace std;
 
  tuple<RCP<const Basic>, RCP<const Basic>> get_base_exp(const RCP<const Basic> &);

  vec_basic get_multuple(const RCP<const Basic>&);
  vec_basic get_addtuple(const RCP<const Basic>&);

  pair<RCP<const Basic>, RCP<const Number>>
  canonical_form(const RCP<const FunctionSymbol>&);

  pair<RCP<const Number>, RCP<const Basic>> max_cnt_var(const vec_basic &v);
  pair<RCP<const Number>, RCP<const Basic>> get_cnt_var_func(const RCP<const FunctionSymbol> &expr);

  bool eq(const vec_basic &, const vec_basic &);
  bool even(const Integer &);

  RCP<const Integer> operator""_i(unsigned long long);
}

#endif // __UTILS_H_
