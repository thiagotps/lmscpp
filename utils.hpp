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

  pair<RCP<const Number>, RCP<const Basic>> max_cnt_var(const vec_basic &);
  pair<RCP<const Number>, RCP<const Basic>> get_cnt_var_func(RCP<const Basic>);

  bool eq(const vec_basic &, const vec_basic &);
  bool even(const Integer &);

  template<typename T>
  using Iterator = typename T::iterator;

  template<typename T>
  class prefix{
    T v_{};
    Iterator<T> b_,e_;
  public:
    prefix(T&& v, int start): v_{v}, b_{v_.begin()+start}, e_{v_.end()}  {}
    prefix(T& v, int start): b_{v.begin()+start}, e_{v.end()} {}
    Iterator<T>& begin() {return b_;}
    Iterator<T>& end() {return e_;}
    const Iterator<T>& begin() const {return b_;}
    const Iterator<T>& end() const {return e_;}
  };

  using vec_func = vector<RCP<const Function>>;
  vec_func states_vars(const RCP<const Basic> &);


  RCP<const Integer> operator""_i(unsigned long long);
}

#endif // __UTILS_H_
