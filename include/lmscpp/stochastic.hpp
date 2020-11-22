/*
 * This is the start point of the library. Here we define basic classes like
 * StochasticProcess, ExpectedOperator, EquationSet and Experiment. In order to make it easier
 * to calculate the evolution of the state variables matrix, here is included the iterator num_stv_iter.
 */


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

#include <lmscpp/utils.hpp>
#include <lmscpp/nummatrix.hpp>

namespace stochastic{
  using namespace std;
  using namespace SymEngine;
  using moment_type = function<RCP<const Basic> (const vec_basic&,const RCP<const Integer>&)>;

  /*  This is the class representing a stochastic process. It takes the name of the process and
      a function used to take the moments of random variables taken from that process.

      The function must receive as the first parameter
      the list of arguments passed to a particular random variable taken from this process. The second
      argument is the exponent of that random variable. As an example, for a process named 'u' and random
      variable 'u(x)**n' the first argument will be '{x}' and the second will be 'n'.

      The function can return null to indicate that the variables from this process don't have a moment
      defined for a give argument or exponent.
  */
  class StochasticProcess{
    // A hash table linking a process' name to its moment function.
    static unordered_map<string,moment_type> moment_;
    // The name of the process.
    const string name_;
  public:
    // The constructor
    StochasticProcess(string name, moment_type moment);

    // Return a function symbol with the same name of this process and with 'arg' as the argument list.
    RCP<const FunctionSymbol> operator()(const vec_basic &arg) const;
    // The same from above, but with only one argument.
    RCP<const FunctionSymbol> operator()(const RCP<const Basic> &arg) const;

    // Receive a Basic object 'f' representing a random variable. It can be the FunctionSymbol returned
    // by the operator () or a Pow in which the base is a FunctionSymbol returned from that operator.
    static RCP<const Basic> moment(const RCP<const Basic> &f);
    // Check if 'sym' is a random variable. That is, if it can be passed to the method named 'moment'.
    static bool is_random(const RCP<const Basic> &sym);

    // TODO: This is only a quick WorkAround
    static inline void clear() {moment_.clear();}
  };


  /*  This class implements the notion of a expected operator. It is a functional whose constructor
      takes only a function that check if two random variables are independent. */
  class ExpectedOperator{
    // The function responsible for checking if two random variables are independent.
    function<bool (const FunctionSymbol&,const FunctionSymbol&)> is_independent_;
  public:
    ExpectedOperator(function<bool (const FunctionSymbol&,const FunctionSymbol&)> is_ind_func);

    // Returns E(arg) with the following simplifications:
    // 1. If 'arg' is a number, then just return 'arg'.
    // 2. If 'arg' is random variable, then return its moment if it has one.
    // 3. In the case that 'arg' does not have a moment or it is not a random variable, just return
    // a FunctionSymbol whose name is 'E' and whose argument is 'arg'.
    RCP<const Basic> operator()(const RCP<const Basic> &arg) const;
    //  Basic can carry a FunctionSymbol or a Pow. It returns whether the two random variables are independent or not.
    bool is_independent(const Basic& b1, const Basic& b2) const;

    // NOTE: This should be turned into a generator when C++20 becomes available.
    //  Takes a expression and breaks it into a list of independent random variables.
    // The last element is the part of the expression that cannot be breaked into independent parts.
    vec_basic split_independents(const RCP<const Basic> &expr) const ;

    // This method takes a FunctionSymbol in the form E(...) and try to expand it accordingly with
    // the properties of the expected operator. Example: if it receives E(u(x)**2 + 2*v(x+1) + n(x-1)) and
    // the second moment of u(x) is γ2 , the first moment of v(x+1) is β and n(x-1) doesn't have a moment then the result of this
    // method will be γ2 + 2*β + E(n(x-1)).
    RCP<const Basic> expand(const RCP<const FunctionSymbol> &term) const;
  };

  /*  This class is a set in which each element represents a state variable's update equation.
      For example, the update equation can be something like E[v(k+1)**2] = γE[v(k)**2] + μE[v(k)**2 * u(k-1)**4] + 4*β**2.
      This class can be used for equations that are not necessarily an update for  some state variable(e.g: v(k+1) = γv(k) + n(k))

      To accomplish this we take the canonical_form of both sides of an equation, get the hash of the left side and
      save both this hash and the canonical_form of the right hand side in an unordered_map.
      */
  class EquationSet{
    // Map the hash of the left hand side of an equation against its right hand side.
    unordered_map<size_t, RCP<const Basic>> eqs_;
    // The independent variable in the equation.
    RCP<const Basic> var_;
  public:
    // The constructor takes the independent variable of the update equation.
    EquationSet(const RCP<const Basic> &var);
    // Get the right hand side given its left hand side.
    RCP<const Basic> getitem(const RCP<const FunctionSymbol> &lhs) const;

    // getitem directly by its canonical form's hash.
    // If you are unsure if your expression is canonical or not then you shouldn't use
    // this method.
    RCP<const Basic> getitem(size_t hash) const;
    // Save an equation of the form lhs = rhs in this set.
    void setitem(const RCP<const FunctionSymbol> &lhs,const RCP<const Basic> &rhs);
    // How many equations are actually stored.
    size_t size() const;

    // Get the independent variable.
    auto get_var() const {return var_;};
    // Check if there is any equation whose left hand side has the form of 'key'.
    bool contains(const RCP<const FunctionSymbol> & key) const;
  };

  // NOTE: Obsolete
  // size_t coumpute_eqs(const vec_func &, const EquationSet&, const ExpectedOperator &);
  // The first argument is any expression of the form u0(k+a0)**n0 * u1(k+a1)**n1 * ...
  // The second is an set of equations containing expansions for each of the elements passed in the first argument.
  RCP<const Basic> expand(const RCP<const Basic>&, const EquationSet &);

  using nmatrix = NumMatrix::NumMatrix<double>;

  // This iterator takes into its constructor numerical values for the matrix A, Y0 and B.
  // Each time that the operator ++ is called an iteration of the form Yk+1 = A*Yk + B is realized..
  // The operator * returns the actual value for Yk.
  class num_stv_iter
  {
    const nmatrix a_,b_;
    nmatrix y_;
  public:
    using difference_type = ptrdiff_t;
    using value_type = nmatrix;
    using pointer = nmatrix*;
    using reference = nmatrix&;
    using iterator_category = forward_iterator_tag;

    num_stv_iter(const nmatrix& a, const nmatrix &b, const nmatrix &y0): a_{a}, b_{b}, y_{y0} {};
    inline const nmatrix& operator*() const {return y_;};

    // Perform Yk+1 = A*Yk + B
    inline const nmatrix& operator++()
    {
      y_ = a_ * y_ + b_;
      return y_;
    };
  };

  /* This is a help class that can compute the equations, mount both symbolic and numerical matrices, save/load
     these matrices into/from a cache, etc ...*/
  using fallback_func_type = function<RCP<const Basic> (const Symbol &x)>;
  class Experiment
  {
    DenseMatrix A_,B_,Yk_;
    size_t number_of_eqs_;
    //TODO: These references are dangerous. Maybe I should use a smart pointer instead.
    const EquationSet& inieqs_;
    const ExpectedOperator& E_;
    const map_basic_basic& inivalsmap_;
    const fallback_func_type fallback_; // Why I can't put a & here ?

  public:
    // Convert the symbolic DenseMatrix in the numerical nmatrix.
    nmatrix sym2num(const DenseMatrix&) const;

    /* 'inieqs' is the initial equations from the model, 'E' is the expected operator,
       'inivalsmap' is a map associating each constant with its value (which can be symbolic or numeric).
       When, for a symbol, a initial value cannot be found, the non null value returned by 'fallback' will be used.
    */
    Experiment(const EquationSet & inieqs, const ExpectedOperator & E,
               const map_basic_basic& inivalsmap = map_basic_basic{},
               fallback_func_type fallback = fallback_func_type{}): inieqs_{inieqs}, E_{E}, inivalsmap_{inivalsmap},
                                                                       fallback_{fallback}
    {};

    // Save the symbolical matrices and the number of equations in a file.
    void save(fstream&) const;
    // Load the symbolical matrices and the number of equations from a file.
    void load(fstream&);
    // Compute all update equations beginning from this vector of state variables.
    void compute(const vec_func &);
    inline const DenseMatrix & get_A() const {return A_;};
    inline const DenseMatrix & get_Yk() const {return Yk_;};
    inline const DenseMatrix & get_B() const {return B_;};
    nmatrix get_num_A() const {return sym2num(A_);};
    DenseMatrix get_sym_Y0() const;
    nmatrix get_num_Y0() const {return sym2num(get_sym_Y0());}
    nmatrix get_num_B() const {return sym2num(B_);}
    inline size_t get_number_of_eqs() const {return number_of_eqs_;};

    void write_skewness(int niter,RCP<const Basic> randexpr, ofstream & os) const;
    pair<int, double> skewness_steady_state(RCP<const Basic> randexpr) const;
    void write_expression(int niter,RCP<const Basic> randexpr, ofstream & os) const;
  };

}


#endif // __STOCHASTIC_H_
