#include "common.hpp"

double numfactorial(uintmax_t n)
{
  double f{1};
  while (n)
    {
      f *= n;
      n--;
    }
  return f;
}

double num_doublefactorial(uintmax_t n)
{
  return n <= 1 ? 1 : n*num_doublefactorial(n - 2);
}
