#ifndef __COMMON_H_
#define __COMMON_H_

#include <chrono>
#include <iostream>
#include <string>

using namespace std;
using namespace chrono;

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

double numfactorial(uintmax_t);
double num_doublefactorial(uintmax_t);

#endif // __COMMON_H_
