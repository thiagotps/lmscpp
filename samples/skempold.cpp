#include <fstream>
#include <functional>
#include <iostream>
#include <array>
#include <cmath>
#include <math.h>
#include <random>
#include <regex>
#include <sstream>
#include <stdexcept>
#include <string>
#include <tuple>
#include <vector>
#include <queue>
#include <future>
#include <limits>
#include <chrono>

#include <argparse/argparse.hpp>

using namespace std;
using namespace chrono;

typedef unsigned long long ull;

template <typename A>
A sgn(A a)
{
  if (a < 0)
    return -1;

  return 1;
}

class laplace_distribution
{
  uniform_real_distribution<> _uni{-1.0/2,1.0/2};
  double _location, _scale;
public:

  laplace_distribution(double location, double scale): _location{location}, _scale{scale} 
  {

  }

  template<class Generator>
  double operator()(Generator &g)
  {
    double u = _uni(g);
    return _location - _scale*sgn(u)*log(1 - 2*abs(u));
  }

};


template <typename T>
struct neg_idx_vector: public vector<T>
{
  int offset = 0; 		// Default no offset
  using vector<T>::vector;
  T& operator[](int idx)
  {
    return vector<T>::operator[](idx+offset);
  }
};

struct result_struct
{
  // NOTE: This could be replaced by a valarray
  vector<double> cm1, cm2, cm3, error_square, pdf;
};

enum
  {
   SKEWNESS = 1 << 0,
   MSE = 1 << 2,
   PDF = 1 << 3,
   WTILDE = 1 << 4,
  };

struct experiment
{
  int N, M, niter;
  double beta, sigmav;
  int modes;

  int pdf_instant;
  double pdf_start, pdf_end;
  int pdf_samples, kernel_exp;
  string dist;

  // TODO: Make this function only callable when the struct is instantiated with const.
  template<typename T>
  auto compute(ull nsamples, T seed) const {
    double sigmau = 1;
    // TODO: The parameters b are still passed manually in the code.
    array<double, 10> b = {1., 0.8, 0.8, -0.7, 0.6, -0.5, 0.4, -0.3, 0.2, -0.1};

    default_random_engine generator(seed);
    normal_distribution<double> u_dist_gauss(0.0, sigmau);
    laplace_distribution u_dist_lap(0.0, 1/sqrt(2));
    normal_distribution<double> nu_dist(0.0, sigmav);

    function<double(default_random_engine &)> u_dist;
    if (dist == "gauss")
      u_dist = [&u_dist_gauss](default_random_engine &g) {return u_dist_gauss(g);};
    else if (dist == "lap")
      u_dist = [&u_dist_lap](default_random_engine &g) {return u_dist_lap(g);};
    else
      throw runtime_error("Unknow dist: " + dist);

    vector<double> v(niter + 1);
    neg_idx_vector<double> u(niter + M);
    u.offset = M - 1;

    // auto x = [&u,&b,M](int k)
    //     {
    // 	double xk = 0;
    // 	for (auto j = 0; j < M; j++)
    // 	  xk += b[j] * u[k - j];
    // 	return xk;
    //     };

    result_struct rs;
    if (modes & (SKEWNESS | WTILDE | PDF) ) {
      rs.cm1.resize(niter + 1, 0);
      rs.cm2.resize(niter + 1, 0);
      rs.cm3.resize(niter + 1, 0);
    }

    if (modes & MSE)
      rs.error_square.resize(niter + 1, 0);

    if (modes & PDF)
      rs.pdf.resize(pdf_samples+1, 0);

    vector<double> wtildek(N);
    vector<double> wtilde_kplus1(N);

    for (ull iter = 0; iter < nsamples; iter++) {
      wtildek.assign(N, 1);

      for (auto &val : v)
        val = nu_dist(generator);
      for (auto &val : u)
        val = u_dist(generator);

      // This must be updated after the u vector.
      neg_idx_vector<double> x(niter + N);
      x.offset = N - 1;
      for (auto k = 1 - N; k <= niter; k++) {
        x[k] = 0;
        for (auto j = 0; j < M; j++)
          x[k] += b[j] * u[k - j];
      }

      for (auto k = 0; k <= niter; k++) {

        // Central Moments of W₀(k).
        if (modes & (SKEWNESS | WTILDE | PDF) ) {
          rs.cm1[k] += wtildek[0];
          rs.cm2[k] += (wtildek[0] * wtildek[0]);
          rs.cm3[k] += (wtildek[0] * wtildek[0] * wtildek[0]);
        }

        if (modes & MSE) {
          double error = v[k];
          for (auto i = 0; i < N; i++)
            error += wtildek[i] * x[k - i];

          // TODO: This can cause overflow.
          rs.error_square[k] += error * error;
        }

        if (modes & PDF and k == pdf_instant){
          const double sigma = pow(10,kernel_exp);

          double delta = (pdf_end - pdf_start)/pdf_samples;
          for (auto n = 0; n <= pdf_samples; n++){
            double x = pdf_start + n*delta;
            double t = exp(-((x - wtildek[0]) * (x - wtildek[0])) /
                           (2 * sigma * sigma));
            t /= sqrt(2 * M_PI) * sigma;
            rs.pdf[n] += t;
          }
        }

        for (auto i = 0; i < N; i++) {
          double s = 0;
          for (auto t = 0; t < N; t++)
            s += x[k - t] * wtildek[t];

          wtilde_kplus1[i] =
            wtildek[i] - beta * x[k - i] * s + beta * v[k] * x[k - i];
        }
        wtildek = wtilde_kplus1;

      }
    }

    return rs;
  }
};

// Test if the random_device is truly random. If the hardware doesn't support random numbers
// the random_device will fallback to a pseudo random generator and thus the arrays ar1 and ar2
// will be necessarily equal. If the hardware supports true random numbers, then
// the probability of this function returns false (i.e: A false negative) is 10**-20.
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

const auto INVALID_INT = numeric_limits<int>::max();
const auto INVALID_DOUBLE = numeric_limits<double>::quiet_NaN();

bool is_invalid(double val) {return isinf(val);}
bool is_invalid(int val) {return val == INVALID_INT;}

int main(int argc, char **argv)
{
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

  // TODO: The program should verify the correct relation between pdf-instant and the other
  // pdf parameters.
  program.add_argument("--pdf-instant")
    .help("The instant k to which the estimated PDF refers to.")
    .default_value(INVALID_INT)
    .action([](const string& val){return stoi(val);});

  program.add_argument("--pdf-start")
    .help("The first sample of the PDF to generate.")
    .default_value(INVALID_DOUBLE)
    .action([](const string& val){return stod(val);});

  program.add_argument("--pdf-end")
    .help("The last sample of the PDF to generate.")
    .default_value(INVALID_DOUBLE)
    .action([](const string& val){return stod(val);});

  program.add_argument("--pdf-samples")
    .help("The number of samples to describe the PDF in the interval [pdf-start, pdf-end]")
    .default_value(INVALID_INT)
    .action([](const string& val){return stoi(val);});

  program.add_argument("--kernel-exp")
    .help("The value of the standard-deviation of the gaussian kernel will be pow(10, kernel_exp).")
    .default_value(INVALID_INT)
    .action([](const string& val){return stoi(val);});

  program.add_argument("--pdf-file")
    .help("The file where the PDF's evolution of the first filter coefficient at a determined instant k will be stored.")
    .default_value(string(""))
    .action([](const string& val){return val;});

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

  program.add_argument("--mse-file")
    .help("the file where the evolution of the MSE of the first filter's coefficient will be stored")
    .default_value(string(""))
    .action([](const string& val){return val;});

  program.add_argument("--wfile")
    .help("The file where the evolution of the the first filter coefficient")
    .default_value(string(""))
    .action([](const string& val){return val;});

  program.add_argument("--ncpu")
    .help("The number of cpus to be used. Default is the number of cpus in the system.")
    .default_value(thread::hardware_concurrency())
    .action([](const string& val){return stoi(val);});

  try {
    program.parse_args(argc, argv);
  } catch (const runtime_error &err) {
    cout << err.what() << endl;
    cout << program;
    exit(0);
  }

  int N = program.get<int>("--filter-length");
  int M = program.get<int>("--data-length");
  int niter = program.get<int>("--niter"); 
  int exp = program.get<int>("--exp");
  ull nsamples = pow(10, exp);
  double beta = program.get<double>("--beta");
  double sigmav = sqrt(program.get<double>("--sigmav2"));
  string skewness_file = program.get<string>("--skewness-file");
  string mse_file = program.get<string>("--mse-file");
  const auto NCPU = program.get<unsigned int>("--ncpu");
  int pdf_instant = program.get<int>("--pdf-instant");
  double pdf_start = program.get<double>("--pdf-start");
  double pdf_end = program.get<double>("--pdf-end");
  int pdf_samples = program.get<int>("--pdf-samples");
  int kernel_exp = program.get<int>("--kernel-exp");
  string pdf_file = program.get<string>("--pdf-file");
  string dist = program.get<string>("--dist");
  string w_file = program.get<string>("--wfile");

  if (not pdf_file.empty())
    if (is_invalid(pdf_instant)  or is_invalid(pdf_start)  or is_invalid(pdf_end)  or is_invalid(pdf_samples)  or is_invalid(kernel_exp))
      throw runtime_error{"The following options are needed when using --pdf-file: --pdf-instant, --pdf-start, "s +
          "--pdf-end, --pdf-samples, --kernel-exp"s};

  if (not is_random_device_random())
    {
      cerr << "This device does not support true random numbers" << endl;
      return 1;
    }

  // The initializer {} avoids implicit conversion.
  ull indv_samples{ nsamples/NCPU };
  ull rest_samples{ nsamples % NCPU };

  int modes = 0;
  if (not skewness_file.empty())
    modes |= SKEWNESS;

  if (not mse_file.empty())
    modes |= MSE;

  if (not pdf_file.empty())
    modes |= PDF;

  if (not w_file.empty())
    modes |= WTILDE;


  const experiment e{.N = N, .M = M, .niter = niter, .beta = beta,
               .sigmav = sigmav, .modes = modes, .pdf_instant = pdf_instant,
               .pdf_start = pdf_start, .pdf_end = pdf_end, .pdf_samples = pdf_samples,
               .kernel_exp = kernel_exp, .dist = dist};

  random_device rd;
  auto current_time = high_resolution_clock::now().time_since_epoch().count();
  uniform_int_distribution<decltype(current_time)> uni{0, current_time};

  auto compute = [e](auto nsamples, auto seed) { return e.compute(nsamples, seed); };
  vector<decltype(async(compute, nsamples, current_time))> r;

  if (NCPU > 1)
    for (auto i = 0; i < NCPU-1; i++)
      r.push_back(async(compute,indv_samples, uni(rd)));

  auto rs = compute(indv_samples+rest_samples, uni(rd));
  for (auto &ft:r)
    {
      auto rs1 = ft.get();
      for (auto k = 0; k <= niter; k++)
        {
          if (modes & (SKEWNESS | WTILDE | PDF) ){
            rs.cm1[k] += rs1.cm1[k];
            rs.cm2[k] += rs1.cm2[k];
            rs.cm3[k] += rs1.cm3[k];
          }

          if (modes & MSE)
            rs.error_square[k] += rs1.error_square[k];
        }

      if (modes & PDF)
        for (auto n = 0; n <= e.pdf_samples; n++)
          rs.pdf[n] += rs1.pdf[n];
    }

  for (auto k = 0; k <= niter; k++) {
    if (modes & (SKEWNESS | WTILDE | PDF) ){
      rs.cm1[k] /= nsamples;
      rs.cm2[k] /= nsamples;
      rs.cm3[k] /= nsamples;
    }

    if (modes & MSE)
      rs.error_square[k] /= nsamples;
  }


  if (modes & PDF) {
    ofstream pdf(pdf_file);
    pdf.precision(numeric_limits<double>::max_digits10);
    double u = rs.cm1[e.pdf_instant];
    double sd = sqrt(rs.cm2[e.pdf_instant] - u * u);
    pdf << "#Mean,Standard-Deviation" << endl;
    pdf << "#" << u << "," << sd << endl;

    double delta = (e.pdf_end - e.pdf_start)/e.pdf_samples;
    for (auto n = 0; n <= e.pdf_samples; n++) {
      rs.pdf[n] /= nsamples;
      pdf << e.pdf_start + delta * n << " " << rs.pdf[n] << endl;
    }
  }

  if (modes & SKEWNESS) {
    ofstream skfile(skewness_file);

    // We print the SKEWNESS' evolution for W₀(k+1)
    for (auto k = 0; k <= niter; k++) {
      double u = rs.cm1[k+1];
      double sd = sqrt(rs.cm2[k+1] - u * u);
      double skewness = (rs.cm3[k+1] - 3 * u * sd * sd - u * u * u) / (sd * sd * sd);
      skfile << k << " " << skewness << endl;
    }
  }

  if (modes & MSE){
    ofstream msefile(mse_file);
    for (auto k = 0; k <= niter; k++)
      msefile << k << " " << rs.error_square[k] << endl;
  }

  if (modes & WTILDE){
    ofstream wtilde_file(w_file);
    for (auto k = 0; k <= niter; k++)
      wtilde_file << k << " " << rs.cm1[k] << endl;
  }

  return 0;
}
