/*
 * This header declares functions and classes of variable uses, such as get_multuple, Prefix,
 * even, state_vars, etc...
  */
#ifndef __UTILS_H_
#define __UTILS_H_

#include <symengine/expression.h>
#include <symengine/functions.h>
#include <symengine/symbol.h>
#include <symengine/matrix.h>
#include <tuple>

namespace stochastic{
  using namespace SymEngine;
  using namespace std;
 
  // This function returns a tuple where the first element is the base and the second is the
  // exponent. If the argument does not contain a Pow, then the base is just the same expression
  // and the exponent will be one.
  tuple<RCP<const Basic>, RCP<const Basic>> get_base_exp(const RCP<const Basic> &);

  // NOTE: This should be converted into a generator when c++20 is available.
  //  If the argument is a Mul, then returns each element of the multiplication with the first
  //  element being the constant of the multiplication.
  //  If the argument is not a Mul than it will be returned as the second element while the first will
  //  be simply a 'one'.
  vec_basic get_multuple(const RCP<const Basic>&);
  // NOTE: This should be converted into a generator when c++20 is available.
  //  If the argument is a Add, then returns each element of the summation with the first
  //  element being the constant of the summation.
  //  If the argument is not a Add than it will be returned as the second element while the first will
  //  be simply a 'zero'.
  vec_basic get_addtuple(const RCP<const Basic>&);

  // The argument is a vector containing powers of FunctionSymbol
  // (i.e: It can be a FunctionSymbol itself or a Pow where the base is a FunctonSymbol).
  // The returned value is pair containing the maximum constant of the expression as its first value.
  // The second value is the independent variable.
  // Ex: The argument {u(k + 1)**3, v(k - 1), u(k+2)**2, n(k-2)} results in {2, k}
  pair<RCP<const Number>, RCP<const Basic>> max_cnt_var(const vec_basic &);
  // The argument can be an FunctionSymbol describing a expected operator  or it can be a multiplication made by powers of FunctionSymbols.
  // If it is a expected operator, the function will act only on its arguments, which should be multiplication made by powers of FunctionSymbols.
  // The return is a pair where the first element is the maximum constant of the expression and the
  // second is the independent variable.
  // Ex: With the argument u(k + 1)*v(k - 2)*u(k + 2) the result will be {2, k}
  // Ex: With the argument E(u(k + 1)*v(k)*u(k - 2)) the result will be {1, k}
  pair<RCP<const Number>, RCP<const Basic>> get_cnt_var_func(RCP<const Basic>);

  // NOTE: Maybe this should be moved to stochastic.hpp
  // The source code of this function is small (two lines) and self explanatory.
  pair<RCP<const Basic>, RCP<const Number>>
  canonical_form(const RCP<const FunctionSymbol>&);


  // Check if the two vec_basic are equal.
  bool eq(const vec_basic &, const vec_basic &);
  // Check if the integer is even.
  bool even(const Integer &);

  template<typename T>
  using Iterator = typename T::iterator;

  // This class can be used in situations where we need to loop through a vector but skiping its first
  // n elements. Ex: for (const auto & t: prefix{v,1}) loop through 'v' skiping its first element.
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

  using vec_func = vector<RCP<const FunctionSymbol>>;
  // Take an expression and return a vector of states vars found in this expression.
  vec_func states_vars(const RCP<const Basic> &);

  // Give a expression like 2*gamma*E(v(k+1)*v(k+2)) + 3*alpha + 4, this function
  // gives a vector of tuples in which the first term of each tuple is the multiplier constant
  // and the second is the state variable's hash. The sum constant, in this example 3*alpha + 4, is always
  // placed in the first slot of the tuple whereas its second element is the actual vector.
  tuple<RCP<const Basic>,vector<tuple<RCP<const Basic>,size_t>>>
  cnt_st_terms(const RCP<const Basic> &);


  inline RCP<const Integer> operator""_i(unsigned long long i) {return integer(i);}
}

#endif // __UTILS_H_
