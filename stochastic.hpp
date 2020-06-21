#ifndef __STOCHASTIC_H_
#define __STOCHASTIC_H_

#include <symengine/expression.h>
#include <symengine/functions.h>
#include <symengine/symbol.h>
#include <functional>
#include <initializer_list>
#include <cstdint>
#include <string>

namespace stochastic{
  using namespace std;
  using namespace SymEngine;
  using moment_type = function<RCP<const Basic> (const vec_basic&,const RCP<const Integer>&)>;

  class StochasticProcess{
    static unordered_map<string,moment_type> moment_;
    const string name_;
  public:
    StochasticProcess(string name, moment_type moment);
    RCP<const FunctionSymbol> operator()(const vec_basic &arg);
    RCP<const FunctionSymbol> operator()(const RCP<const Basic> &arg);

    static RCP<const Basic> moment(const RCP<const Basic> &f);
    static bool is_random(const RCP<const Basic> &sym);
  };


  class ExpectedOperator{
    function<bool (const FunctionSymbol&,const FunctionSymbol&)> is_independent_;
  public:
    ExpectedOperator(function<bool (const FunctionSymbol&,const FunctionSymbol&)> is_ind_func);

    RCP<const Basic> operator()(const RCP<const Basic> &arg) const;
    // NOTE: Basic can carry a FunctionSymbol or a Pow
    bool is_independent(const Basic& b1, const Basic& b2) const;

    vec_basic split_independents(const RCP<const Basic> &expr) const ;

    RCP<const Basic> expand(const RCP<const FunctionSymbol> &term) const;
  };

  class EquationSet{
    unordered_map<size_t, RCP<const Basic>> eqs_;
    RCP<const Basic> var_;
  public:
    EquationSet(const RCP<const Basic> &var);
    RCP<const Basic> getitem(const RCP<const FunctionSymbol> &lhs);
    void setitem(const RCP<const FunctionSymbol> &lhs,const RCP<const Basic> &rhs);
    size_t size() const;
  };

}


#endif // __STOCHASTIC_H_
