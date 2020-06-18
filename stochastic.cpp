#include "stochastic.hpp"
#include "utils.hpp"
#include <symengine/basic-inl.h>

#include <numeric>

namespace stochastic{
  unordered_map<string,moment_type> StochasticProcess::moment_{};


  StochasticProcess::StochasticProcess(string name, moment_type moment): name_{name}
  {
    moment_[name_] = moment;
  }
  RCP<const FunctionSymbol> StochasticProcess::operator()(const vec_basic &arg){
    return make_rcp<const FunctionSymbol>(name_, arg);
  }

  RCP<const FunctionSymbol> StochasticProcess::operator()(const RCP<const Basic> &arg){
    return make_rcp<const FunctionSymbol>(name_, arg);
  }

   RCP<const Basic> StochasticProcess::moment(const RCP<const Basic> &f) {
    auto [b,n] = get_base_exp(f);
    auto bf = rcp_dynamic_cast<const FunctionSymbol>(b);
    auto m = moment_[bf->get_name()];
    //  return m(bf->get_args(), down_cast<Integer>(n));
    return m(bf->get_args(), rcp_dynamic_cast<const Integer>(n));
  }

  bool StochasticProcess::is_random(const RCP<const Basic> &sym){
    auto [s,b] = get_base_exp(sym);
    if (is_a<const FunctionSymbol>(*s)){
      auto f = rcp_dynamic_cast<const FunctionSymbol>(s);
      return moment_.find(f->get_name()) != moment_.end();
    }

    return false;
  }

  ExpectedOperator::ExpectedOperator(function<bool (const FunctionSymbol&,const FunctionSymbol&)> is_ind_func):
    is_independent_{is_ind_func} {}

  RCP<const Basic> ExpectedOperator::operator()(const RCP<const Basic> &arg) const{
    if (is_a_Number(*arg))
      return arg;

    if (StochasticProcess::is_random(arg))
      return StochasticProcess::moment(arg);

   return make_rcp<const FunctionSymbol>("E", arg);
  }

  bool ExpectedOperator::is_independent(const Basic& b1, const Basic& b2) const {
    // NOTE: This will throw a exception if something diferent of a Pow and FunctionSymbol
    // is passed in. If the Pow's base is not a FunctionSymbol, this will raise a exception too.
    auto convert =  [](const Basic &arg) -> const FunctionSymbol& {
                      if (is_a<Pow>(arg)) {
                        auto base = static_cast<const Pow&>(arg).get_base();
                        return dynamic_cast<const FunctionSymbol&>(*base);
                      }
                      return dynamic_cast<const FunctionSymbol&>(arg);
                    };

    const auto& f1{convert(b1)};
    const auto& f2{convert(b2)};
    return is_independent_(f1, f2);
  };

  vec_basic ExpectedOperator::split_independents(const RCP<const Basic> &expr) const
  {
    if (not is_a<Mul>(*expr))
      return {expr};

    auto m = rcp_static_cast<const Mul>(expr);
    RCP<const Basic> not_ind{one};

    auto rvs{get_multuple(m)};
    rvs.erase(begin(rvs));


    vec_basic res{};
    for (auto &x : rvs)
      {
        bool ind = true;
        for (auto &y : rvs)
          {
            if (x != y and not is_independent(dynamic_cast<const FunctionSymbol&>(*x),dynamic_cast<const FunctionSymbol&>(*y)))
              {
                not_ind = mul(not_ind, x);
                ind = false;
                break;
              }
          }

        if (ind)
          res.push_back(x);
      }
    res.push_back(not_ind);

    return res;
  }

  RCP<const Basic> ExpectedOperator::expand(const RCP<const FunctionSymbol> &term) const
  {
    vec_basic addtuple{get_addtuple(SymEngine::expand(term->get_args()[0]))};
    auto s{addtuple[0]};

    addtuple.erase(begin(addtuple));
    for (auto &addterm : addtuple){
      vec_basic multuple{get_multuple(addterm)};
      auto cnts = multuple[0];
      RCP<const Basic> not_cnts{one};

      multuple.erase(begin(multuple));
      for (auto &multerm : multuple){
        if (StochasticProcess::is_random(multerm))
          not_cnts = mul(not_cnts, multerm);
        else
          cnts = mul(cnts, multerm);
      }

      vec_basic e{};
      for (auto &t : split_independents(not_cnts))
        e.emplace_back((*this)(t));

      auto mul_e = reduce(begin(e), end(e), rcp_dynamic_cast<const Basic>(one), [](auto &a1, auto &a2)
                                                        {
                                                          return mul(a1,a2);
                                                        });

      s = add(s,mul(cnts, mul_e));
    }

    return s;
  }

}
