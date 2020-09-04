#include <iostream>
#include <cmath>
#include <random>
#include <chrono>
#include <valarray>

#include <argparse.hpp>

using namespace::std;

inline double kernel_pdf(double x, double mean, double sigma)
{
  return exp(-(x - mean)*(x - mean)/(2.0*sigma*sigma))/(sigma * sqrt(2.0 * M_PI));
}


int main(int argc, char ** argv) {

  argparse::ArgumentParser program(argv[0]);

  program.add_argument("--mean")
    .help("Mean")
    .required()
    .action([](const string& val){return stod(val);});

  program.add_argument("--sd")
    .help("Standard Deviation")
    .required()
    .action([](const string& val){return stod(val);});

  program.add_argument("--start")
    .help("Start")
    .required()
    .action([](const string& val){return stod(val);});

  program.add_argument("--end")
    .help("End")
    .required()
    .action([](const string& val){return stod(val);});

  program.add_argument("--samples")
    .help("Number of samples to take in the interval [start, end]")
    .required()
    .action([](const string& val){return stoi(val);});

  program.add_argument("--exp")
    .help("The number of Monte Carlo trials will be 10**exp")
    .required()
    .action([](const string& val){return stoi(val);});

  program.add_argument("--kernel-exp")
    .help("The gaussian kernel used")
    .required()
    .action([](const string& val){return stod(val);});

  program.add_argument("--square")
    .help("Standard Deviation")
    .default_value(false)
    .implicit_value(true);

  try {
    program.parse_args(argc, argv);
  } catch (const runtime_error &err) {
    cout << err.what() << endl;
    cout << program;
    exit(0);
  }

  const auto mean = program.get<double>("--mean");
  const auto sd = program.get<double>("--sd");
  const auto start = program.get<double>("--start");
  const auto end = program.get<double>("--end");
  const auto samples = program.get<int>("--samples");
  const auto exp = program.get<int>("--exp");
  const auto kernel_exp = program.get<double>("--kernel-exp");
  const auto square = program.get<bool>("--square");

  const double delta{(end - start)/samples};
  const uintmax_t monte_carlo_samples{static_cast<uintmax_t>(pow(10.0, exp))};
  const double sigma {pow(10.0, kernel_exp)};

  valarray<double> pdf(0.0, samples + 1);

  mt19937_64 gen{static_cast<uintmax_t>(chrono::high_resolution_clock::now().time_since_epoch().count())};
  normal_distribution<double> gauss{mean, sd};

  for (uintmax_t it = 0; it < monte_carlo_samples; it++)
    {
      auto g = gauss(gen);
      if (square)
        g *= g;

      for (auto n = 0; n <= samples; n++)
        {
          double x{start + n*delta};
          pdf[n] += kernel_pdf(x, g, sigma);
        }
    }

  pdf /= monte_carlo_samples;

  for (auto n = 0; n <= samples; n++)
    {
      double x{start + n*delta};
      cout << x << " " << pdf[n] << endl;
    }


  return 0;
}
