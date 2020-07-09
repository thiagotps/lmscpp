#include <numeric>
#include <stdexcept>
#include <algorithm>

#include <lmscpp/nummatrix.hpp>

namespace NumMatrix
{
  using namespace std;

  template<typename T>
  NumMatrix<T>& NumMatrix<T>::operator+=(const NumMatrix<T>& a)
  {
    throws_if_dimension_mismatch(a);
    auto n{a.nrows()}, m{a.ncols()};

    // NOTE: This loop could be parallelised
    for (auto i = 0; i < n; i++)
      for (auto j = 0; j < m; j++)
        m_[i][j] += a[i][j];

    return *this;
  }

  template<typename T>
  NumMatrix<T>& NumMatrix<T>::operator*=(const NumMatrix<T>& a)
  {
    throws_if_dimension_mismatch(a);
    auto n{a.nrows()}, m{a.ncols()};

    NumMatrix<T> c(n,m);
    // NOTE: This loop could be parallelised
    for (auto i = 0; i < n; i++)
      for (auto j = 0; j < m; j++)
        c[i][j] = inner_product(begin(a[j]), end(a[j]), (*this)[i]);

    return *this;
  }

  template<typename T>
  void NumMatrix<T>::throws_if_dimension_mismatch(const NumMatrix<T> &a)
  {
    if (not (this->nrows() == a.nrows() and this->ncols() == a.ncols()))
      throw runtime_error{"Matrix dimension mismatch"};
  }
};
