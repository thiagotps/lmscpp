#ifndef __STOCHASTIC_H_
#define __STOCHASTIC_H_

#include <symengine/expression.h>
#include <symengine/functions.h>
#include <symengine/symbol.h>
#include <symengine/matrix.h>
#include <functional>
#include <initializer_list>
#include <cstdint>
#include <string>
#include "utils.hpp"

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

    // TODO: This is only a quick WorkAround
    static inline void clear() {moment_.clear();}
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
    RCP<const Basic> getitem(const RCP<const FunctionSymbol> &lhs) const;

    //    getitem by its canonical form's hash.
    // If you are unsure if your expression is canonical or not then you shouldn't use
    // this method.
    RCP<const Basic> getitem(size_t hash) const;
    void setitem(const RCP<const FunctionSymbol> &lhs,const RCP<const Basic> &rhs);
    size_t size() const;

    auto get_var() const {return var_;};
    bool contains(const RCP<const FunctionSymbol> & key) const;
  };

  // NOTE: Obsolete
  size_t coumpute_eqs(const vec_func &, const EquationSet&, const ExpectedOperator &);
  RCP<const Basic> expand(const RCP<const Basic>&, const EquationSet &);


  class Experiment
  {
    DenseMatrix A_,B_,Yk_;
    size_t number_of_eqs_;
    const EquationSet& inieqs_;
    const ExpectedOperator& E_;
  public:
    Experiment(const EquationSet & inieqs, const ExpectedOperator & E): inieqs_{inieqs}, E_{E} {};
    void save(ofstream&) const;
    void load(ifstream&);
    void compute(const vec_func &);
    inline const DenseMatrix & get_A() const {return A_;};
    inline const DenseMatrix & get_Yk() const {return Yk_;};
    inline const DenseMatrix & get_B() const {return B_;};
    inline size_t get_number_of_eqs() const {return number_of_eqs_;};
  };
}


#endif // __STOCHASTIC_H_
