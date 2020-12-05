#include <iostream>
#include <fstream>
#include <chrono>

#include <stdexcept>
#include <symengine/expression.h>
#include <symengine/functions.h>
#include <symengine/symbol.h>
#include <symengine/ntheory.h>
#include <symengine/logic.h>
#include <symengine/matrix.h>
#include <symengine/printers/latex.h>

#include <argparse.hpp>

#include <lmscpp/stochastic.hpp>
#include <lmscpp/utils.hpp>
#include <lmscpp/operators.hpp>
#include <lmscpp/nummatrix.hpp>

#include "common.hpp"

using namespace argparse;
using namespace SymEngine;
using namespace SymEngine::OverloadedOperators;
using namespace std;
using namespace stochastic;

using namespace chrono;


// NOTE: This a redefininition. Maybe I should create a header file with definitions.
using vec_func = vector<RCP<const FunctionSymbol>>;

enum cachemode
  {
   IGNORE=0,
   READ,
   WRITE,
  };
//  NOTE: This should be shared with classical.cpp
void read_write_compute(const string& cachename, cachemode cmode, Experiment& todo, const vec_func & seeds)
{
  MeasureDuration duration{};

  if (cmode == cachemode::READ)
    {
      fstream fs{cachename, fstream::binary | fstream::in};
      if (fs.is_open())
        {
          cout << "Using cached data..." << endl;
          duration.reset();
          todo.load(fs);
          duration.show("todo.load()");
        }
      else
        throw runtime_error{"It was not possible to open cache " + cachename + " for reading."};
    }
  else
    {
      duration.reset();
      todo.compute(seeds);
      duration.show("todo.compute()");

      if (cmode == cachemode::WRITE)
        {
          cout << "Writing cache..." << endl;
          fstream fs{cachename, fstream::binary | fstream::out};

          duration.reset();
          todo.save(fs);
          duration.show("todo.save()");
        }
    }
}

enum indmode
  {
   IA = 0,
   EEA,
  };

enum outmode
  {
   SK = 0,
   MSE,
  };

enum distmode
  {
   UNSPECIFIED = 0,
   GAUSS,
   LAP,
  };




string latex(const DenseMatrix& d) {
  string out{};

  auto n{d.nrows()}, m{d.ncols()};
  for (auto i = 0; i < n; i++) {
    for (auto j = 0; j < m; j++){
      out += latex(*d.get(i, j)) + " &";
    }
    out += " \\\\\n";
  }

  return out;
}

const double binary_precision{pow(10,-3)};
double binary(function<bool (double)> f, double lo, double hi, double ml) {
  double m{(lo + hi)/2.0};
  if (abs(m - ml) < binary_precision*abs(ml))
    return m;

  return f(m) ? binary(f, m, hi, m) : binary(f, lo, m, m);
}

const auto sigmav{ symbol("σ_v") }, step_size{ symbol("β") }, k{ symbol("k") };

double get_max_beta(double lo, double hi, const DenseMatrix & m,
                map_basic_basic ini,
                const fallback_func_type fallback) {

  return binary([&](double beta) {
    ini[step_size] = number(beta);
    return NumMatrix::max_eigen_value(sym2num(m, ini, fallback)) < 1;
  }, lo, hi, lo);
}

int main(int argc, char ** argv)
{
  argparse::ArgumentParser program(argv[0]);

  program.add_argument("--readcache")
    .help("The cache file to use.")
    .default_value(""s)
    .action([](const string &val){return val;});

  program.add_argument("--writecache")
    .help("The cache file to write.")
    .default_value(""s)
    .action([](const string &val){return val;});

  program.add_argument("-N").help("filter length").required().action([](const string &val){return stoi(val);});

  program.add_argument("-M").help("data length").required().action([](const string &val){return stoi(val);});

  program.add_argument("--indmode").help("ia or eea").required()
    .action([](const string &val)
            {
              static const map<string, indmode> choices {{"ia", indmode::IA}, {"eea", indmode::EEA}};
              auto it = choices.find(val);
              if (it == choices.end())
                throw runtime_error{"Invalid value for --iamode"};

              return it->second;
            });

  program.add_argument("--outmode").help("sk or mse").required()
    .action([](const string &val)
            {
              static const map<string, outmode> choices {{"sk", outmode::SK}, {"mse", outmode::MSE}};
              auto it = choices.find(val);
              if (it == choices.end())
                throw runtime_error{"Invalid value for --outmode"};

              return it->second;
            });

  program.add_argument("-b","--beta").help("β").default_value(-1.0)
    .action([](const string &val){return stod(val);});

  program.add_argument("--sv2","--sigmav2").help("variance (σᵥ²)").default_value(-1.0)
    .action([](const string &val){return stod(val);});

  program.add_argument("-n","--niter").help("Number of iterations.").default_value(-1)
    .action([](const string &val){return stoi(val);});

  program.add_argument("-o","--output").help("The file where the output will be stored.").default_value(""s)
    .action([](const string &val){return val;});

  program.add_argument("-d", "--dist").help("gauss or lap")
    .default_value(distmode::UNSPECIFIED)
    .action([](const string &val)
            {
              static const map<string, distmode> choices {{"gauss", distmode::GAUSS}, {"lap", distmode::LAP}};
              auto it = choices.find(val);
              if (it == choices.end())
                throw runtime_error{"Invalid value for --dist"};

              return it->second;
            });

  program.add_argument("--steady-state")
    .help("The steady-state value. Currently only support the skewness.")
    .implicit_value(true)
    .default_value(false);

  program.add_argument("--latex")
    .help("print the matrices to the standard output ")
    .implicit_value(true)
    .default_value(false);

  program.add_argument("--max-beta").help("get the maximum beta").default_value(0.0)
    .action([](const string &val){return stod(val);});

  // TODO: Maybe instead of an int we should return a uintmax_t ?
  program.add_argument("--max-high").help("Max high/depth of the recursion!").default_value(numeric_limits<int>::max())
    .action([](const string &val){return stoi(val);});

  try {
    program.parse_args(argc, argv);
  } catch (const runtime_error &err) {
    cout << err.what() << endl;
    cout << program;
    exit(1);
  }

  const auto readcache{program.get<string>("--readcache")};
  const auto writecache{program.get<string>("--writecache")};
  const auto N{program.get<int>("-N")};
  const auto M{program.get<int>("-M")};
  const auto beta{program.get<double>("--beta")};
  const auto sigmav2{program.get<double>("--sigmav2")};
  const auto niter{program.get<int>("--niter")};
  const auto ofilename {program.get<string>("--output")};
  const auto ind_mode = program.get<indmode>("--indmode");
  const auto out_mode = program.get<outmode>("--outmode");
  const auto dist_mode = program.get<distmode>("--dist");
  const auto steady_state = program.get<bool>("--steady-state");
  const auto print_latex = program.get<bool>("--latex");
  const auto max_beta = program.get<double>("--max-beta");
  const auto max_high = program.get<int>("--max-high");


  if (not ofilename.empty())
    if (beta < 0 or sigmav2 < 0 or niter < 0 or dist_mode == distmode::UNSPECIFIED)
      throw runtime_error{"Missing some of the following: --beta, --sv2, --niter, --dist "};



  const string cachename{ (not readcache.empty()) ? readcache : writecache };
  const auto cache_mode = [&]()
                          {
                            if (not readcache.empty())
                              return cachemode::READ;
                            if (not writecache.empty())
                              return cachemode::WRITE;

                            return cachemode::IGNORE;
                          }();

  map_basic_basic cache{ {sigmav, number(sqrt(sigmav2))}, {step_size, number(beta)} };

  StochasticProcess::clear();

  StochasticProcess u{"u",  [](const vec_basic&, const RCP<const Integer>& n) -> RCP<const Basic>
                      {
                       if (n->as_int() % 2 == 0)
                         return symbol("γ_" + to_string(n->as_int()));

                       return zero;
                      }};

  StochasticProcess v{"v", [](const vec_basic&, const RCP<const Integer>& n) -> RCP<const Basic>
                           {
                             if (even(*n))
                               return pow(sigmav, n);

                             return zero;
                           }};

  vector<StochasticProcess> wtil;
  for (auto i = 0; i < N; i++)
    wtil.emplace_back("ẅ_" + to_string(i), [](const vec_basic& vec, const RCP<const Integer>& n) -> RCP<const Basic>
                           {
                             if (eq(*vec[0], *zero))
                               return one;

                            return null;
                           });

  int recursion_depth = 0;
  ExpectedOperator E{[ind_mode, max_high, &recursion_depth](const FunctionSymbol& x, const FunctionSymbol& y, int high){
    recursion_depth = max(recursion_depth, high);
    auto local_ind_mode = high <= max_high ? ind_mode : indmode::IA;
                             auto xy = [local_ind_mode, max_high](const FunctionSymbol &x, const FunctionSymbol &y)
                                       {
                                         auto xname{x.get_name()};
                                         auto yname{y.get_name()};

                                         if (xname == "u" and yname == "u")
                                           return not eq(x,y);

                                         // NOTE: This is extremely ugly. Maybe there is some better way to do this
                                         // comparasion.
                                         if (xname.find("ẅ") != string::npos and yname == "u")
                                           {
                                             if (local_ind_mode == indmode::EEA)
                                               return not rcp_dynamic_cast<const Integer>(expand(x.get_args()[0] - y.get_args()[0]))->is_positive();
                                             else if (local_ind_mode == indmode::IA)
                                               return true;
                                           }

                                         if (xname == "v")
                                           return true;

                                         return false;
                                       };
                             return xy(x,y) or xy(y,x);
                           }};

  array bvalues{1., 0.8, 0.8, -0.7, 0.6, -0.5, 0.4, -0.3, 0.2, -0.1};
  if (bvalues.size() < M)
    {
      cerr << "There size of 'bvalues' is less than " << M << endl;
      return 1;
    }

  vector<RCP<const FunctionSymbol>> b;
  for (auto i = 0; i < M; i++)
    {
      auto bi { function_symbol("b", integer(i)) };
      b.push_back( rcp_static_cast<const FunctionSymbol>(bi) );
      cache[bi] = number(bvalues[i]);
    }

  auto x  = [&](const auto & k)
            {
              RCP<const Basic> sum{zero};
              // NOTE: This integer could be memorised.
              for (auto i = 0; i < M; i++)
                sum = sum + b[i] * u(k - integer(i));

              return sum;
            };

  EquationSet eqs{k};

  RCP<const Basic> inn{zero};
  for (auto j = 0; j < N; j++)
    inn = inn + x(k - integer(j))*wtil[j](k);

  for (auto i = 0; i < N; i++)
    eqs.setitem(wtil[i](k+one), wtil[i](k) - step_size*x(k - integer(i)) * inn - step_size*v(k)*x(k - integer(i)));

  const auto fallback = [dist_mode](const Symbol & x) -> RCP<const Basic>
                                 {
                                  auto name = x.__str__();
                                  if (name.size() >= 4 and name.substr(0, 2) == "γ")
                                    {
                                      auto p{ stoll(name.substr(3)) };
                                      if (dist_mode == distmode::LAP)
                                        {
                                          static const double scale{ 1.0/sqrt(2) };
                                          auto res = numfactorial(p)*pow(scale, p);
                                          return number(res);
                                        }
                                      else if (dist_mode == distmode::GAUSS)
                                        {
                                          static constexpr double sd{1.0};
                                          auto res = pow(sd, p) * num_doublefactorial(p - 1);
                                          return number(res);
                                        }
                                    }

                                  return null;
                                 };

  Experiment todo{eqs, E, cache, fallback};

  RCP<const Basic> mse_expanded{null};
  vec_func seeds{};

  if (out_mode == outmode::SK)
    {
      seeds.push_back(rcp_dynamic_cast<const FunctionSymbol>(E(wtil[0](k))));
      seeds.push_back(rcp_dynamic_cast<const FunctionSymbol>(E(pow(wtil[0](k), 2_i))));
      seeds.push_back(rcp_dynamic_cast<const FunctionSymbol>(E(pow(wtil[0](k), 3_i))));
    }
  else if (out_mode == outmode::MSE)
    {
      auto ee = rcp_dynamic_cast<const FunctionSymbol>( E(expand(pow(inn + v(k), 2_i ))) );
      // TODO: The expand should have a default value of 0!
      mse_expanded = E.expand(ee, 0);
      seeds = states_vars(mse_expanded);
    }

  read_write_compute(cachename, cache_mode, todo, seeds);
  cout << "NUMBER_OF_EQUATIONS: " << todo.get_number_of_eqs() << endl;
  cout << "RECURSION_DEPTH: " << recursion_depth << endl;

  MeasureDuration duration{};
  if (not ofilename.empty())
    {
      ofstream os{ofilename, ofstream::out | ofstream::trunc};
      if (!os)
        {
          cerr << "It was not possible to open file " << ofilename << endl;
          return 1;
        }

      if (out_mode == outmode::SK)
        {
          duration.reset();
          todo.write_skewness(niter, wtil[0](k), os);
          duration.show("todo.write_skewness()");
        }
      else if (out_mode == outmode::MSE)
        {
          duration.reset();
          todo.write_expression(niter, mse_expanded, os);
          duration.show("todo.write_expression()");
        }
    }

  if (print_latex)
    {
      cout << "Printing latex...\n" << endl;
      cout << "Matrix A: " << endl;
      cout << latex(todo.get_A()) << endl;

      cout << "Matrix Yk: " << endl;
      cout << latex(todo.get_Yk()) << endl;

      cout << "Matrix B: " << endl;
      cout << latex(todo.get_B()) << endl;

      cout << "Printing numeric matrix...\n" << endl;

      cout << "Matrix A: " << endl;
      cout << todo.get_num_A() << endl;
      cout << "Matrix B: " << endl;
      cout << todo.get_num_B() << endl;
    }

  if (steady_state)
    {
      if (out_mode == outmode::SK) {
        auto p = todo.skewness_steady_state(wtil[0](k));
        cout << "steady-state: " << p.second << endl;
        cerr << "Steady-state found with " << p.first << " iterations." << endl;
      }
      else
        throw runtime_error{"The steady-state value of MSE is currently not supported!"};
    }

  if (max_beta) {
    cout << "beta max = " << get_max_beta(0.0, max_beta, todo.get_A(), cache, fallback) << endl;
  }

  return 0;
}
