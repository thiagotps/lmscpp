#include <iostream>
#include <fstream>
#include <valarray>
#include <random>
#include <string>
#include <future>
#include <algorithm>
#include <set>

#include <argparse.hpp>
#include <pcg_random.hpp>

using namespace std;
using namespace chrono;

struct result_struct
{
  // m1, m2 and m3 are the raw moments.
  valarray<double> m1, m2, m3, m4, error_square, pdf, fourth_moment, msd;
  inline void operator += (const result_struct & o )
  {
    m1 += o.m1; m2 += o.m2; m3 += o.m3; m4 += o.m4; error_square += o.error_square; pdf += o.pdf, msd += o.msd;
  }

  template<typename T>
  inline void operator /= (T n)
  {
    m1 /= n; m2 /= n; m3 /= n; m4 /= n; error_square /= n; pdf /= n, msd /= n;
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

template <typename A>
A sgn(A a)
{
  if (a < 0)
    return -1;

  return 1;
}

class laplace_distribution
{
  uniform_real_distribution<> _uni{-1.0/2.0,1.0/2.0};
  double _location, _scale;
public:

  laplace_distribution(double location, double scale): _location{location}, _scale{scale} {}

  template<class Generator>
  double operator()(Generator &g)
  {
    double u = _uni(g);
    return _location - _scale*sgn(u)*log(1 - 2*abs(u));
  }

};

enum
  {
   SKEWNESS = 1 << 0,
   MSE = 1 << 2,
   PDF = 1 << 3,
   FOUR = 1 << 4,
   MSD = 1 << 5,
  };


struct experiment
{
  const int N, M, niter;
  const double beta, sigmav;
  const uint16_t modes;
  const string dist;

  const int pdf_instant;
  const double pdf_start, pdf_end;
  const int pdf_samples, kernel_exp;
  const bool pdf_square;

  inline double normal_pdf(double x, double mean) const
  {
    const double sigma {pow(10, kernel_exp)};
    return  exp(-(x - mean)*(x - mean)/(2.0*sigma*sigma))/(sigma * sqrt(2.0 * M_PI));
  }

  template<typename T>
  auto compute(uintmax_t nsamples, T seed) const
  {
    static constexpr double sigmau {1};
    static constexpr array b {1., 0.8, 0.8, -0.7, 0.6, -0.5, 0.4, -0.3, 0.2, -0.1};

    pcg64 gen{seed};
    normal_distribution<double> u_dist_gauss{0.0, sigmau};
    laplace_distribution u_dist_lap(0.0, 1.0/sqrt(2.0));
    normal_distribution<double> nu_dist_gauss{0.0, sigmav};

    function<double()> u_dist;
    if (dist == "gauss")
      u_dist = [&gen,&u_dist_gauss]() {return u_dist_gauss(gen);};
    else if (dist == "lap")
        u_dist = [&gen,&u_dist_lap]() {return u_dist_lap(gen);};
    else
      throw runtime_error("Unknow dist: " + dist);

    auto nu_dist = [&gen,&nu_dist_gauss]() {return nu_dist_gauss(gen);};

    vector<double> v(niter + 1);
    neg_idx_vector<double> u(niter + N + M - 1, N + M - 2 );

    result_struct rs;
    if (modes & (SKEWNESS | PDF))
      {
        rs.m1.resize(niter + 1, 0);
        rs.m2.resize(niter + 1, 0);
        rs.m3.resize(niter + 1, 0);
      }

    if (modes & MSE)
      rs.error_square.resize(niter + 1, 0);

    if (modes & MSD)
      rs.msd.resize(niter + 1, 0);

    if (modes & PDF)
      rs.pdf.resize(pdf_samples+1, 0);

    if (modes & FOUR)
      rs.m4.resize(niter + 1, 0);

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
            if (modes & (SKEWNESS | PDF))
              {
                rs.m1[k] += wtilk[0];
                rs.m2[k] += wtilk[0]*wtilk[0];
                rs.m3[k] += wtilk[0]*wtilk[0]*wtilk[0];
              }

            if (modes & FOUR)
              rs.m4[k] += wtilk[0]*wtilk[0]*wtilk[0]*wtilk[0];

            if (modes & MSE)
              {
                double error = v.at(k);
                for (auto i = 0; i < N; i++)
                  error += wtilk.at(i)*x.at(k - i);

                rs.error_square[k] += error*error;
              }

              if (modes & MSD) {
                rs.msd[k] = 0.0;
                for (int i = 0; i < N; i++)
                  rs.msd[k] += wtilk.at(i)*wtilk.at(i);
              }

            if (modes & PDF and k == pdf_instant)
              {
                double delta { (pdf_end - pdf_start)/pdf_samples };
                for (auto n = 0; n <= pdf_samples; n++)
                  {
                    double x {pdf_start + n*delta};
                    rs.pdf[n] += normal_pdf(x,  pdf_square ? wtilk.at(0)*wtilk.at(0) : wtilk.at(0) );
                  }
              }

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

// Return a vector with n distinct random uintmax_t randomly choosen
// from 1 to max value holded by this type.
auto gen_random_seeds(int n)
{
  set<uintmax_t> s;
  mt19937_64 gen{chrono::high_resolution_clock::now().time_since_epoch().count()};
  uniform_int_distribution<uintmax_t> uri{1, numeric_limits<uintmax_t>::max()};

  while (s.size() != n)
    s.insert(uri(gen));

  return vector(s.begin(), s.end());
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

  program.add_argument("--ncpu")
    .help("The number of cpus to be used. Default is the number of cpus in the system.")
    .default_value(thread::hardware_concurrency())
    .action([](const string& val){return (unsigned int)stoi(val);});

  program.add_argument("--skewness-file","--sk-file")
    .help("the file where the evolution of the skewness of the first filter's coefficient will be stored")
    .default_value(string(""))
    .action([](const string& val){return val;});

  program.add_argument("--mse-file")
    .help("the file where the evolution of the MSE of the first filter's coefficient will be stored")
    .default_value(string(""))
    .action([](const string& val){return val;});

  program.add_argument("--msd-file")
    .help("the file where the evolution of the MSD will be stored")
    .default_value(string(""))
    .action([](const string& val){return val;});


  program.add_argument("--pdf-instant")
    .help("The instant k to which the estimated PDF refers to.")
    .default_value(0)
    .action([](const string& val){return stoi(val);});

  program.add_argument("--pdf-start")
    .help("The first sample of the PDF to generate.")
    .default_value(0.0)
    .action([](const string& val){return stod(val);});

  program.add_argument("--pdf-end")
    .help("The last sample of the PDF to generate.")
    .default_value(0.0)
    .action([](const string& val){return stod(val);});

  program.add_argument("--pdf-samples")
    .help("The number of samples to describe the PDF in the interval [pdf-start, pdf-end]")
    .default_value(0)
    .action([](const string& val){return stoi(val);});

  program.add_argument("--kernel-exp")
    .help("The value of the standard-deviation of the gaussian kernel will be pow(10, kernel_exp).")
    .default_value(-7)
    .action([](const string& val){return stoi(val);});

  program.add_argument("--pdf-square")
    .help("Return the PDF of the square.")
    .default_value(false)
    .implicit_value(true);

  program.add_argument("--pdf-file")
    .help("The file where the PDF's evolution of the first filter coefficient at a determined instant k will be stored.")
    .default_value(string(""))
    .action([](const string& val){return val;});

  program.add_argument("--fourth-file")
    .help("The file where the fourth moment of the first desviation coefficient will be stored.")
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
  const auto mse_file = program.get<string>("--mse-file");
  const auto msd_file = program.get<string>("--msd-file");
  const auto pdf_file = program.get<string>("--pdf-file");
  const auto fourth_file = program.get<string>("--fourth-file");
  const auto dist = program.get<string>("--dist");
  const auto NCPU = program.get<unsigned int>("--ncpu");
  const auto pdf_instant = program.get<int>("--pdf-instant");
  const auto pdf_start = program.get<double>("--pdf-start");
  const auto pdf_end = program.get<double>("--pdf-end");
  const auto pdf_samples = program.get<int>("--pdf-samples");
  const auto kernel_exp = program.get<int>("--kernel-exp");
  const auto pdf_square = program.get<bool>("--pdf-square");

  uint16_t modes{0};
  if (not skewness_file.empty())
    modes |= SKEWNESS;

  if (not mse_file.empty())
    modes |= MSE;

  if (not msd_file.empty())
    modes |= MSD;

  if (not pdf_file.empty())
    modes |= PDF;

  if (not fourth_file.empty())
    modes |= FOUR;

  auto random_seeds = gen_random_seeds(NCPU);

  const experiment e{.N = N, .M = M, .niter = niter, .beta = beta,
    .sigmav = sigmav, .modes = modes, .dist = dist,
    .pdf_instant = pdf_instant, .pdf_start = pdf_start, .pdf_end = pdf_end,
    .pdf_samples = pdf_samples, .kernel_exp = kernel_exp,
    .pdf_square = pdf_square};

  auto q{nsamples/NCPU};
  auto r{nsamples % NCPU};
  auto func = [&e](auto n, auto seed) { return e.compute(n, seed); };


  auto as1 = async(func, q + r, random_seeds[0]);
  vector<decltype(as1)> v_async {};
  v_async.emplace_back(move(as1));

  for (auto i = 1; i <= NCPU - 1; i++)
    v_async.emplace_back(async(func, q, random_seeds.at(i)));

  auto res = v_async[0].get();
  for (auto it = v_async.begin()+1; it != v_async.end(); it++)
    res += it->get();

  res /= nsamples;

  if (modes & SKEWNESS)
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

  if (modes & MSE)
    {
      ofstream msefile{mse_file};
      for (auto k = 0; k <= niter; k++)
        msefile << k << " " << res.error_square[k] << endl;
    }

  if (modes & MSD)
    {
      ofstream msdfile{msd_file};
      for (auto k = 0; k <= niter; k++)
        msdfile << k << " " << res.msd[k] << endl;
    }


  if (modes & PDF)
    {
      ofstream pdf {pdf_file};
      double u { res.m1[pdf_instant] };
      double sd { sqrt(res.m2[e.pdf_instant] - u*u) };

      pdf << "#Mean,Standard-Deviation" << endl;
      pdf << "#" << u << "," << sd << endl;

      double delta {(pdf_end - pdf_start)/pdf_samples};

      for (auto n = 0; n <= pdf_samples; n++)
        pdf << pdf_start + delta*n << " " << res.pdf[n] << endl;
    }

  if (modes & FOUR) {
    ofstream fourfile {fourth_file};

    for (auto k = 0; k <= niter; k++)
      fourfile << k << " " << res.m4[k] << endl;
  }

  return 0;
}
