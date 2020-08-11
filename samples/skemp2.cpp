#include <iostream>
#include <fstream>
#include <valarray>
#include <random>
#include <string>
#include <future>
#include <algorithm>

#include <argparse.hpp>

using namespace std;
using namespace chrono;

struct result_struct
{
  // m1, m2 and m3 are the raw moments.
  valarray<double> m1, m2, m3;
  inline void operator += (const result_struct & o )
  {
    m1 += o.m1;
    m2 += o.m2;
    m3 += o.m3;
  }

  template<typename T>
  inline void operator /= (T n)
  {
    m1 /= n;
    m2 /= n;
    m3 /= n;
  }
};

template <typename T>
class neg_idx_vector
{
  const size_t offset_;
  vector<T> vec_;
public:
  neg_idx_vector(size_t n, size_t offset): offset_{offset}, vec_(n) {}

  inline T& at(int idx) { return vec_.at(idx + offset_); }
  inline const T& at(int idx) const { return vec_.at(idx + offset_); }

  inline auto begin() { return vec_.begin(); }
  inline auto end() { return vec_.end(); }
  inline const auto begin() const { return vec_.begin(); }
  inline const auto end() const { return vec_.end(); }
};

struct experiment
{
  const int N, M, niter;
  const double beta, sigmav;

  template<typename T>
  auto compute(uintmax_t nsamples, T seed) const
  {
    static constexpr double sigmau {1};
    static constexpr array b {1., 0.8, 0.8, -0.7, 0.6, -0.5, 0.4, -0.3, 0.2, -0.1};

    default_random_engine gen{seed};
    normal_distribution<double> u_dist_gauss{0.0, sigmau};
    normal_distribution<double> nu_dist_gauss{0.0, sigmav};
    auto u_dist = [&gen,&u_dist_gauss]() {return u_dist_gauss(gen);};
    auto nu_dist = [&gen,&nu_dist_gauss]() {return nu_dist_gauss(gen);};

    vector<double> v(niter + 1);
    neg_idx_vector<double> u(niter + N + M - 1, N + M - 2 );

    result_struct rs;
    rs.m1.resize(niter + 1, 0);
    rs.m2.resize(niter + 1, 0);
    rs.m3.resize(niter + 1, 0);

    vector<double> wtilk(N);
    vector<double> wtilkplus1(N);

    for (uintmax_t iter = 0; iter < nsamples; iter++)
      {
        // wtil(0) = [1 ... 1].transpose()
        fill(wtilk.begin(), wtilk.end(), 1);

        for (auto &val : v)
          val = nu_dist();

        for (auto &val : u)
          val = u_dist();

        neg_idx_vector<double> x(niter + N, N - 1);
        for (auto k = 1 - N; k <= niter; k++)
          {
            x.at(k) = 0;
            for (auto m = 0; m < M; m++)
              x.at(k) += b.at(m) * u.at(k - m);
          }

        // Execute one trial of the experiment.
        for (auto k = 0; k <= niter; k++)
          {
            // Update the raw moments
            rs.m1[k] += wtilk[0];
            rs.m2[k] += wtilk[0]*wtilk[0];
            rs.m3[k] += wtilk[0]*wtilk[0]*wtilk[0];

            double inner_prod{0};
            for (auto j = 0; j < N; j++)
              inner_prod += x.at(k - j) * wtilk.at(j);

            for (auto i = 0; i < N; i++)
              wtilkplus1.at(i) = wtilk.at(i) - beta*x.at(k - i)*inner_prod
                - beta*x.at(k - i)*v.at(k);

            wtilk = wtilkplus1;
          }

      }

    return rs;
  }
};


bool is_random_device_random()
{
  random_device rd1,rd2;
  // default_random_engine rd1,rd2;
  uniform_int_distribution<int> uni{1,100};

  array<int,10> ar1, ar2;
  for (auto &n : ar1)
    n = uni(rd1);

  for (auto &n : ar2)
    n = uni(rd2);

  return ar1 != ar2;
}

int main(int argc, char ** argv) {

  argparse::ArgumentParser program(argv[0]);
  program.add_argument("--filter-length","-N")
    .help("filter length ")
    .required()
    .action([](const string& val){return stoi(val);});

  program.add_argument("--data-length", "-M")
    .help("data model length ")
    .required()
    .action([](const string& val){return stoi(val);});

  program.add_argument("--niter")
    .help("number of iterations ")
    .required()
    .action([](const string& val){return stoi(val);});

  program.add_argument("--exp")
    .help("the number of Monte Carlos trials will be 10^exp ")
    .required()
    .action([](const string& val){return stoi(val);});

  program.add_argument("--beta")
    .help("step-size parameter β ")
    .required()
    .action([](const string& val){return stod(val);});

  program.add_argument("--sigmav2","--sv2")
    .help("variance (σᵥ²) ")
    .required()
    .action([](const string& val){return stod(val);});

  program.add_argument("--dist")
    .help("The distribution used for the input signal. Possible values are gauss and lap(for laplacian).")
    .required()
    .action([](const string& val)
            {
              if (val != "gauss" and val != "lap")
                throw runtime_error("--dist must be gauss or lap");

              return val;
            });

  program.add_argument("--skewness-file","--sk-file")
    .help("the file where the evolution of the skewness of the first filter's coefficient will be stored")
    .default_value(string(""))
    .action([](const string& val){return val;});

  try {
    program.parse_args(argc, argv);
  } catch (const runtime_error &err) {
    cout << err.what() << endl;
    cout << program;
    exit(0);
  }

  const auto N = program.get<int>("-N");
  const auto M = program.get<int>("-M");
  const auto niter = program.get<int>("--niter");
  const auto exp = program.get<int>("--exp");
  const uintmax_t nsamples { static_cast<uintmax_t>(pow(10, exp)) };
  const auto beta = program.get<double>("--beta");
  const auto sigmav = sqrt(program.get<double>("--sigmav2"));
  const auto skewness_file = program.get<string>("--skewness-file");
  const auto dist = program.get<string>("--dist");
  const auto NCPU { thread::hardware_concurrency() };

  if (not is_random_device_random())
    {
      cerr << "This device does not support true random numbers" << endl;
      return 1;
    }

  random_device d;
  uniform_int_distribution<uintmax_t> intdist{1,numeric_limits<uintmax_t>::max()};

  const experiment e{.N = N, .M = M, .niter = niter, .beta = beta,
    .sigmav = sigmav};

  auto q{nsamples/NCPU};
  auto r{nsamples % NCPU};
  auto func = [&e](auto n, auto s) { return e.compute(n, s); };

  auto as1 = async(func, q + r, intdist(d));
  vector<decltype(as1)> v_async {};
  v_async.emplace_back(move(as1));

  for (auto i = 0; i < NCPU - 1; i++)
    v_async.emplace_back(async(func, q, intdist(d)));

  auto res = v_async[0].get();
  for (auto it = v_async.begin()+1; it != v_async.end(); it++)
    res += it->get();

  res /= nsamples;

  if (not skewness_file.empty())
    {
      ofstream skfile{skewness_file};
      for (auto k = 0; k <= niter; k++)
        {
          double u {res.m1[k]};
          double sd { sqrt(res.m2[k] - u*u) };
          double sk { (res.m3[k] - 3*u*sd*sd - u*u*u)/(sd*sd*sd) };
          skfile << k << " " << sk << endl;
        }
    }

  return 0;
}
