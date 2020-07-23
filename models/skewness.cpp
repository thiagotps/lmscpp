#include <iostream>
#include <fstream>
#include <chrono>

#include <symengine/expression.h>
#include <symengine/functions.h>
#include <symengine/symbol.h>
#include <symengine/ntheory.h>
#include <symengine/logic.h>
#include <symengine/matrix.h>
#include <argparse.hpp>

#include <lmscpp/stochastic.hpp>
#include <lmscpp/utils.hpp>
#include <lmscpp/operators.hpp>

using namespace argparse;
using namespace SymEngine;
using namespace SymEngine::OverloadedOperators;
using namespace std;
using namespace stochastic;

using namespace chrono;

// NOTE: This was copied from classical.cpp, maybe I should put this in its own library
class MeasureDuration
{
  time_point<high_resolution_clock> t0{high_resolution_clock::now()};
public:
  void reset() {t0 = high_resolution_clock::now();}
  void show(string fname) const
  {
    auto t1{ high_resolution_clock::now() };
    cout << fname
         << " taked " <<  duration_cast<seconds>(t1 - t0).count() << "s." << endl;
  }
};

// NOTE: This was copied from classical.cpp, maybe I should put this in its own library
uintmax_t numfactorial(uintmax_t n)
{
  uintmax_t f{1};
  while (n)
    {
      f *= n;
      n--;
    }
  return f;
}

// NOTE: This a redefininition. Maybe I should create a header file with definitions.
using vec_func = vector<RCP<const FunctionSymbol>>;

// NOTE: This should be shared with classical.cpp
void read_write_compute(string filename, bool readcache, bool writecache, Experiment& todo, const vec_func & seeds)
{
  MeasureDuration duration{};
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
      todo.compute(seeds);
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

enum indmode
  {
   IA = 1,
   EEA,
  };

int main(int argc, char ** argv)
{
  argparse::ArgumentParser program(argv[0]);
  program.add_argument("--readcache")
    .help("Read cache if available.")
    .default_value(false)
    .implicit_value(true);

  program.add_argument("--writecache")
    .help("If the matrices was computed, instead of read from the cache, write them.")
    .default_value(false)
    .implicit_value(true);

  program.add_argument("-N").help("filter length").required().action([](const string &val){return stoi(val);});

  program.add_argument("-M").help("data length").required().action([](const string &val){return stoi(val);});

  program.add_argument("--beta").help("β").required().default_value(0.0)
    .action([](const string &val){return stod(val);});

  program.add_argument("--sigmav2").help("variance (σᵥ²)").required().default_value(0.0)
    .action([](const string &val){return stod(val);});

  program.add_argument("--niter").help("Number of iterations.").required().default_value(0)
    .action([](const string &val){return stoi(val);});

  program.add_argument("--sk").help("The file where the skewness' evolution will be stored.").default_value(""s)
    .action([](const string &val){return val;});

  program.add_argument("--mse").help("The file where the MSE's evolution will be stored.").default_value(""s)
    .action([](const string &val){return val;});

  program.add_argument("--mode").help("ia or eea").required()
    .action([](const string &val)
            {
              static const map<string, indmode> choices {{"ia", indmode::IA}, {"eea", indmode::EEA}};
              auto it = choices.find(val);
              if (it == choices.end())
                throw runtime_error{"Invalid value for --mode"};

              return it->second;
            });



  try {
    program.parse_args(argc, argv);
  } catch (const runtime_error &err) {
    cout << err.what() << endl;
    cout << program;
    exit(1);
  }

  const int N{program.get<int>("-N")};
  const int M{program.get<int>("-M")};
  const double beta{program.get<double>("--beta")};
  const double sigmav2{program.get<double>("--sigmav2")};
  const int niter{program.get<int>("--niter")};
  const string sk_filename{program.get<string>("--sk")};
  const string mse_filename{program.get<string>("--mse")};
  const bool readcache{program.get<bool>("--readcache")};
  const bool writecache{program.get<bool>("--writecache")};
  const auto mode = program.get<indmode>("--mode");

  const auto sigmav{ symbol("σ_v") }, step_size{ symbol("β") }, k{ symbol("k") };
  map_basic_basic cache{ {sigmav, number(sqrt(sigmav2))}, {step_size, number(beta)} };

  StochasticProcess::clear();

  StochasticProcess u{"u",  [](const vec_basic&, const RCP<const Integer>& n) -> RCP<const Basic>
                      {
                       if (n->as_int() % 2 == 0)
                         return symbol("γ_" + to_string(n->as_int()));

                       return zero;
                      }};

  StochasticProcess v{"v", [&sigmav](const vec_basic&, const RCP<const Integer>& n) -> RCP<const Basic>
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

  ExpectedOperator E{[mode](const FunctionSymbol& x, const FunctionSymbol& y){
                             auto xy = [mode](const FunctionSymbol &x, const FunctionSymbol &y)
                                       {
                                         auto xname{x.get_name()};
                                         auto yname{y.get_name()};

                                         if (xname == "u" and yname == "u")
                                           return not eq(x,y);

                                         // NOTE: This is extremely ugly. Maybe there is some better way to do this
                                         // comparasion.
                                         if (xname.find("ẅ") != string::npos and yname == "u")
                                           {
                                             if (mode == indmode::EEA)
                                               return not rcp_dynamic_cast<const Integer>(expand(x.get_args()[0] - y.get_args()[0]))->is_positive();
                                             else if (mode == indmode::IA)
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
    eqs.setitem(wtil[i](k+one), wtil[i](k) - step_size*x(k - integer(i)) * inn + step_size*v(k)*x(k - integer(i)));

  Experiment todo{eqs, E, cache, [](const Symbol & x) -> RCP<const Basic>
                                 {
                                  static const double scale = 1.0/sqrt(2);
                                  auto name = x.__str__();
                                  if (name.size() >= 4 and name.substr(0, 2) == "γ")
                                    {
                                      auto p{ stoll(name.substr(3)) };
                                      auto res = numfactorial(p)*pow(scale, p);
                                      return number(res);
                                    }

                                  return null;
                                 }};

  vec_func seeds{};
  seeds.push_back(rcp_dynamic_cast<const FunctionSymbol>(E(wtil[0](k))));
  seeds.push_back(rcp_dynamic_cast<const FunctionSymbol>(E(pow(wtil[0](k), 2_i))));
  seeds.push_back(rcp_dynamic_cast<const FunctionSymbol>(E(pow(wtil[0](k), 3_i))));

  read_write_compute("cache" + to_string(N) + to_string(M) + ".bin", readcache, writecache, todo, seeds);
  cout << "NUMBER_OF_EQUATIONS: " << todo.get_number_of_eqs() << endl;

  MeasureDuration duration{};
  if (not sk_filename.empty())
    {
      ofstream os{sk_filename, ofstream::out | ofstream::trunc};

      duration.reset();
      todo.write_skewness(niter, wtil[0](k), os);
      duration.show("todo.write_skewness()");
    }
  return 0;
}
