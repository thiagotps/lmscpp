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
}

#endif // __UTILS_H_
