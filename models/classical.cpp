#include <cerrno>
#include <iostream>
#include <fstream>
#include <chrono>

#include <symengine/expression.h>
#include <symengine/functions.h>
#include <symengine/symbol.h>
#include <symengine/ntheory.h>
#include <symengine/logic.h>
#include <symengine/matrix.h>
#include <argparse/argparse.hpp>


#include <lmscpp/stochastic.hpp>
#include <lmscpp/utils.hpp>
#include <lmscpp/operators.hpp>

#include "common.hpp"

using namespace argparse;
using namespace SymEngine;
using namespace SymEngine::OverloadedOperators;
using namespace std;
using namespace stochastic;

using namespace chrono;


// classical -L -M --symmatrix file --nummatrix file --evolution file --force
// NUMBER_OF_EQUATIONS: 1000
int main(int argc, char ** argv)
{
  argparse::ArgumentParser program(argv[0]);
  program.add_argument("-L").help("Filter length").required().action([](const string &val){return stoi(val);});

  program.add_argument("-M").help("Data length").required().action([](const string &val){return stoi(val);});

  program.add_argument("--symmatrix").help("The file where symbolic matrices will be printed.").default_value(""s)
    .action([](const string &val){return val;});

  program.add_argument("--nummatrix").help("The file where numeric matrices will be printed.").default_value(""s)
    .action([](const string &val){return val;});

  program.add_argument("--evolution").help("The file where numeric matrices will be printed.").default_value(""s)
    .action([](const string &val){return val;});

  program.add_argument("--emse").help("The file where Excess MSE evolution will be printed.").default_value(""s)
    .action([](const string &val){return val;});

  program.add_argument("--niter").help("Number of iterations. Only needed when --evolution is passed.").default_value(0)
    .action([](const string &val){return stoi(val);});

  program.add_argument("--readcache")
    .help("Read cache if available.")
    .default_value(false)
    .implicit_value(true);

  program.add_argument("--writecache")
    .help("If the matrices was computed, instead of read from the cache, write them.")
    .default_value(false)
    .implicit_value(true);

  program.add_argument("-b","--beta").help("β").default_value(-1.0)
    .action([](const string &val){return stod(val);});

  program.add_argument("--sv2","--sigmav2").help("variance (σᵥ²)").default_value(-1.0)
    .action([](const string &val){return stod(val);});




  try {
    program.parse_args(argc, argv);
  } catch (const runtime_error &err) {
    cout << err.what() << endl;
    cout << program;
    exit(1);
  }

  const int L{program.get<int>("-L")};
  const int M{program.get<int>("-M")};
  const string symmatrix_filename{program.get<string>("--symmatrix")};
  const string nummatrix_filename{program.get<string>("--nummatrix")};
  const string evolution_filename{program.get<string>("--evolution")};
  const string emse_filename{program.get<string>("--emse")};
  const int niter{program.get<int>("--niter")};
  const auto beta{program.get<double>("--beta")};
  const auto sigmav2{program.get<double>("--sigmav2")};
  const bool readcache{program.get<bool>("--readcache")};
  const bool writecache{program.get<bool>("--writecache")};

  // TODO: It should not be mandatory in all ocasions.
  if (beta < 0 or sigmav2 < 0 or niter < 0)
    throw runtime_error{"Missing some of the following: --beta, --sv2, --niter"};


  const auto sigma{ symbol("σ_n") };
  map_basic_basic cache{ {sigma, number(sqrt(sigmav2))} };

  StochasticProcess::clear();

  StochasticProcess u{"u",  [](const vec_basic&, const RCP<const Integer>& n) -> RCP<const Basic>
                      {
                       if (n->as_int() % 2 == 0)
                         return symbol("γ_" + to_string(n->as_int()));

                       return zero;
                      }};

  StochasticProcess n{"n", [&sigma](const vec_basic&, const RCP<const Integer>& n) -> RCP<const Basic>
                           {
                             if (even(*n))
                               return pow(sigma, n);

                             return zero;
                           }};

  vector<StochasticProcess> V;
  for (auto i = 0; i < L; i++)
    V.emplace_back("v_" + to_string(i), [](const vec_basic& vec, const RCP<const Integer>& n) -> RCP<const Basic>
                           {
                             if (eq(*vec[0], *zero))
                               return one;

                            return null;
                           });

  ExpectedOperator E{[](const FunctionSymbol& x, const FunctionSymbol& y){
                             auto xy = [](const FunctionSymbol &x, const FunctionSymbol &y)
                                       {
                                         auto xname{x.get_name()};
                                         auto yname{y.get_name()};

                                         if (xname == "u" and yname == "u")
                                           return not eq(x,y);

                                         // NOTE: This is extremely ugly. Maybe there is some better way to do this
                                         // comparasion.
                                         if (xname[0] == 'v' and xname[1] == '_' and yname == "u")
                                           return not rcp_dynamic_cast<const Integer>(expand(x.get_args()[0] - y.get_args()[0]))->is_positive();

                                         if (xname == "n")
                                           return true;


                                         return false;
                                       };
                             return xy(x,y) or xy(y,x);
                           }};




  vector<RCP<const FunctionSymbol>> a;
  for (auto i = 0; i < M; i++)
    {
      auto ai { function_symbol("a", integer(i)) };
      a.push_back( rcp_static_cast<const FunctionSymbol>(ai) );
      cache[ai] = number(1.0/(i+1));
    }

  auto x  = [&](const auto & k)
            {
              RCP<const Basic> sum{zero};
              // NOTE: This integer could be memorised.
              for (auto i = 0; i < M; i++)
                sum = sum + a[i] * u(k - integer(i));

              return sum;
            };


  auto step_size{symbol("μ")}, k{symbol("k")};
  cache[step_size] = number(beta);

  EquationSet eqs{k};

  RCP<const Basic> inn{zero};
  for (auto j = 0; j < L; j++)
    inn = inn + x(k - integer(j))*V[j](k);

  for (auto i = 0; i < L; i++)
    eqs.setitem(V[i](k+one), V[i](k) - step_size*x(k - integer(i)) * inn + step_size*n(k)*x(k - integer(i)));

  auto epislonk = E.expand(rcp_dynamic_cast<const FunctionSymbol>(E(expand(pow(inn, 2_i)))));

  MeasureDuration duration;
  Experiment todo{eqs, E, cache, [](const Symbol & x) -> RCP<const Basic>
                                 {
                                  static const double scale = 0.5;
                                  auto name = x.__str__();
                                  if (name.size() >= 4 and name.substr(0, 2) == "γ")
                                    {
                                      auto p{ stoll(name.substr(3)) };
                                      auto res = numfactorial(p)*pow(scale, p);
                                      return number(res);
                                    }

                                  return null;
                                 }};
  // Experiment todo{eqs, E, cache};

  {
    string filename{"cache"+to_string(L)+to_string(M)+".bin"};
    fstream fs{filename, fstream::binary | fstream::in};
    if (readcache and fs.is_open())
      {
        cout << "Using cached data..." << endl;
        duration.reset();
        todo.load(fs);
        duration.show("todo.load()");
      }
    else
      {
        if (readcache)
          cout << "It was not possible to open the cache. Computing equations..." << endl;
        else
          cout << "Computing equations..." << endl;

        duration.reset();
        todo.compute(states_vars(epislonk));
        duration.show("todo.compute()");

        if (writecache)
          {
            cout << "Writing cache..." << endl;
            fs.open(filename, fstream::binary | fstream::out);

            duration.reset();
            todo.save(fs);
            duration.show("todo.save()");
          }
      }
  }

  if (symmatrix_filename != "")
    {
      ofstream os{symmatrix_filename, ofstream::out | ofstream::trunc};
      os << todo.get_A() << endl;
      os << todo.get_Yk() << endl;
      os << todo.get_B() << endl;
    }

  if (nummatrix_filename != "")
    {
      ofstream os{nummatrix_filename, ofstream::out | ofstream::trunc};
      duration.reset();
      os << todo.get_num_A() << "\n" << endl;
      os << todo.get_num_Y0() << "\n" << endl;
      os << todo.get_num_B() << "\n" << endl;
      duration.show("Numerical evaluation of the matrices");
    }

  if (evolution_filename != "")
    {
      ofstream os{evolution_filename, ofstream::out | ofstream::trunc};

      duration.reset();
      num_stv_iter iter{todo.get_num_A(), todo.get_num_B(), todo.get_num_Y0()};
      for (auto i = 0; i < niter; i++)
        {
          os << *iter << "\n" << endl;
          ++iter;
        }
      duration.show("Matrix evolution");
    }

  if (emse_filename != "") {
    ofstream os{emse_filename, ofstream::out | ofstream::trunc};
    todo.write_expression(niter, epislonk, os);
  }

  cout << "NUMBER_OF_EQUATIONS: " << todo.get_number_of_eqs() << endl;

  return 0;
}
