#include "stochastic.hpp"
#include "utils.hpp"
#include <symengine/basic-inl.h>
#include <symengine/subs.h>

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

    if (StochasticProcess::is_random(arg)){
      auto r =  StochasticProcess::moment(arg);
      if (not r.is_null())
        return r;
    }

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
            if (x != y and not is_independent(*x,*y))
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

      auto mul_e = reduce(begin(e), end(e), rcp_static_cast<const Basic>(one), [](auto &a1, auto &a2)
                                                        {
                                                          return mul(a1,a2);
                                                        });

      s = add(s,mul(cnts, mul_e));
    }

    return s;
  }

  EquationSet::EquationSet(const RCP<const Basic> &var): var_{var}, eqs_{} {};

  RCP<const Basic> EquationSet::getitem(const RCP<const FunctionSymbol> &lhs) const
  {
    auto [canon_lhs, cnt] = canonical_form(lhs);
    return xreplace(eqs_.at(canon_lhs->hash()), {{var_, add(var_, cnt)}});
  }

  void EquationSet::setitem(const RCP<const FunctionSymbol> &lhs,const RCP<const Basic> &rhs)
  {
    auto [canon_lhs, cnt] = canonical_form(lhs);
    auto x_rhs = xreplace(rhs, {{var_, sub(var_, cnt)}});
    eqs_[canon_lhs->hash()] = x_rhs;
  }
  size_t EquationSet::size() const {return eqs_.size();}

  bool EquationSet::contains(const RCP<const FunctionSymbol> & key) const
  {
    auto [canon, cnt] = canonical_form(key);
    return eqs_.find(canon->hash()) != eqs_.end();
  }

  RCP<const Basic> expand(const RCP<const Basic>& expr, const EquationSet & eqs)
  {
    auto ml{get_multuple(expr)};
    RCP<const Basic> result{ml[0]};

    for (const auto & x : prefix(ml, 1))
      {
        auto [b, n] = get_base_exp(x);
        if (is_a<FunctionSymbol>(*b))
          if (auto bf =  rcp_static_cast<const FunctionSymbol>(b); eqs.contains(bf))
            {
              result = mul(result, expand(pow(eqs.getitem(bf), n)));
              continue;
            }

        result = mul(result, x);
      }

    return result;
  }

  size_t coumpute_eqs(const vec_func & seed, const EquationSet& inieqs, const ExpectedOperator & E)
  {
    EquationSet eqs{inieqs.get_var()};
    list<RCP<const FunctionSymbol>> s{seed.begin(), seed.end()};
    while (not s.empty())
      {
        auto t = s.front();
        s.pop_front();
        if (not eqs.contains(t))
          {
            // cout << *t << endl;
            // NOTE: The result of this expectation operator should be a FunctionSymbol
            // This should not result in a number. Maybe I should do the expansion inside
            // the operator () of E instead of calling it separately.
            auto eee = E(expand(t->get_args()[0], inieqs));
            auto u = E.expand(rcp_dynamic_cast<const FunctionSymbol>(eee));
            // NOTE: Test
            // u = expand(u);

            auto stvec = states_vars(u);
            move(stvec.begin(), stvec.end(), back_inserter(s));
            eqs.setitem(t, u);
          }
      }
    return eqs.size();
  }
}
