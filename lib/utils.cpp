#include <lmscpp/utils.hpp>
#include <symengine/number.h>
#include <symengine/basic-inl.h>
#include <symengine/subs.h>

namespace stochastic{
  inline RCP<const FunctionSymbol>
  get_if_expected(const RCP<const Basic> &b)
  {
    if (is_a<FunctionSymbol>(*b))
      {
        auto ptr =  rcp_static_cast<const FunctionSymbol>(b);
        if (ptr->get_name() == "E")
          return ptr;
      }

    return null;
  }

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

  pair<RCP<const Number>, RCP<const Basic>> max_cnt_var(const vec_basic &v){
    pair<RCP<const Number>, RCP<const Basic>> coef{NegInf,NegInf};

    for (auto & rv: v){
      auto [b,n] = get_base_exp(rv);
      auto cp = get_addtuple(b->get_args()[0]);

      auto cpn = rcp_dynamic_cast<const Number>(cp[0]);
      if (cpn->sub(*coef.first)->is_positive()){
        coef.first = cpn;
        coef.second = cp[1];
      }
    }

    return coef;
  }

  pair<RCP<const Number>, RCP<const Basic>>
  get_cnt_var_func(RCP<const Basic> expr){
    if (auto ptr = get_if_expected(expr); not ptr.is_null())
      expr = ptr->get_args()[0];

    auto rvvec = get_multuple(expr);
    rvvec.erase(begin(rvvec));

    return max_cnt_var(rvvec);
  }

  pair<RCP<const Basic>, RCP<const Number>>
  canonical_form(const RCP<const FunctionSymbol>& expr)
  {
    auto [cnt, var] = get_cnt_var_func(expr);
    return {xreplace(expr, {{var, sub(var,cnt)}}), cnt};
  }

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


  vec_func states_vars(const RCP<const Basic> & expr)
  {
    vec_func v{};
    auto addlist = get_addtuple(expr);
    for (const auto& addterm : prefix(addlist, 1))
      {
        auto mullist = get_multuple(addterm);
        for (const auto& multerm : prefix(mullist, 1))
          {
            auto m = get<0>(get_base_exp(multerm));
            auto s = get_if_expected(m);
            if (not s.is_null())
              v.push_back(s);
          }
      }
    return v;
  }

  // Give a expression like 2*gamma*E(v(k+1)*v(k+2)) + 3*alpha + 4, this function
  // gives a vector of tuples in which the first term of each tuple is the multiplier constant
  // and the second is the state variable's hash. The sum constant, in this example 3*alpha + 4, is always
  // placed in the first slot of the tuple whereas its second element is the actual vector.
  tuple<RCP<const Basic>,vector<tuple<RCP<const Basic>,size_t>>>
  cnt_st_terms(const RCP<const Basic> & expr)
  {
    auto addlist = get_addtuple(expr);
    auto expr_cnt = addlist[0];

    vector<tuple<RCP<const Basic>,size_t>> result{};
    for (const auto & addterm : prefix(addlist, 1))
      {
        auto mullist = get_multuple(addterm);
        auto mcnt = mullist[0];

        RCP<const FunctionSymbol> st{null};
        for (const auto & multerm : prefix(mullist, 1))
          {
            st = get_if_expected(multerm);
            if (st.is_null())
              mcnt = mul(mcnt, multerm);
          }

        // No state variable was found. We have only a constant.
        if (st.is_null())
          expr_cnt = add(expr_cnt, mcnt);
        else
          result.push_back({mcnt, canonical_form(st).first->hash()});
      }

    return {expr_cnt, result};
  }




}
