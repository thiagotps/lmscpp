#include "utils.hpp"
#include <symengine/number.h>
#include <symengine/basic-inl.h>

namespace stochastic{
  tuple<RCP<const Basic>, RCP<const Basic>> get_base_exp(const RCP<const Basic> &symbol_or_pow){
    if (is_a<const Pow>(*symbol_or_pow)){
      auto p = rcp_static_cast<const Pow>(symbol_or_pow);
      return {p->get_base(), p->get_exp()};
    }
    return {symbol_or_pow, one};
  }

  vec_basic get_multuple(const RCP<const Basic>& mulexpr){
    if (not is_a<Mul>(*mulexpr))
      return {one,mulexpr};

    auto m = rcp_static_cast<const Mul>(mulexpr);
    auto mulvec = m->get_args();

    if (m->get_coef()->is_one())
      mulvec.insert(begin(mulvec), one);

    return mulvec;
  }

  vec_basic get_addtuple(const RCP<const Basic>& addexpr){
    if (not is_a<Add>(*addexpr))
      return {zero, addexpr};

    auto a = rcp_static_cast<const Add>(addexpr);
    auto addtuple = a->get_args();

    if (a->get_coef()->is_zero())
      addtuple.insert(begin(addtuple), zero);

    return addtuple;
  }
}
