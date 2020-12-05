/*
 * This header declares our numerical matrix library.
 * It is small and self explanatory.
  */

#ifndef __NUMMATRIX_H_
#define __NUMMATRIX_H_

#include <vector>
#include <ostream>
#include <cmath>

namespace NumMatrix
{
  using namespace std;

  template <typename T>
  class NumMatrix
  {
    vector<vector<T>> m_;
  public:
    NumMatrix(unsigned int rows, unsigned int cols): m_{rows,vector<T>(cols, 0)} {}
    NumMatrix(initializer_list<vector<T>> l): m_{l} {};


    inline unsigned int  nrows() const {return m_.size();}
    inline unsigned int  ncols() const {return m_[0].size();}
    vector<T>& operator[](unsigned int  idx) {return m_.at(idx);}
    const vector<T>& operator[](unsigned int idx) const {return m_.at(idx);}

    T get(unsigned int i, unsigned int j) const {return (*this)[i][j];}
    void set(unsigned int i, unsigned int j, T t) {(*this)[i][j] = t;}
    auto begin() {return m_.begin();};
    auto end() {return m_.end();};
    auto begin() const {return m_.begin();};
    auto end() const {return m_.end();};
  };

  template<typename T>
  NumMatrix<T> operator+(const NumMatrix<T>& a, const NumMatrix<T>& b)
  {
    if (not (b.nrows() == a.nrows() and b.ncols() == a.ncols()))
      throw runtime_error{"Matrix dimension mismatch while adding."};

    auto n{a.nrows()}, m{a.ncols()};
    NumMatrix<T> tmp{n, m};

    // NOTE: This loop could be parallelised
    for (auto i = 0; i < n; i++)
      for (auto j = 0; j < m; j++)
        tmp[i][j] = a[i][j] + b[i][j];

    return tmp;
  }

  template<typename T>
  NumMatrix<T> operator*(const NumMatrix<T>& a, const NumMatrix<T>& b)
  {
    auto n{a.nrows()}, p{a.ncols()}, q{b.nrows()}, m{b.ncols()};
    if (p != q)
      throw runtime_error{"Matrix dimension mismatch while multiplying."};


    NumMatrix<T> tmp(n,m);
    // NOTE: This loop could be parallelised
    for (auto i = 0; i < n; i++)
      for (auto j = 0; j < m; j++)
        {
          tmp[i][j] = 0;
          for (auto k = 0; k < q; k++)
            tmp[i][j] += a[i][k]*b[k][j];
        }

    return tmp;
  }



  template<typename T>
  ostream& operator << (ostream & os, const NumMatrix<T>& m)
  {
    bool f1 = true;
    for (const auto & l : m)
      {
        if (!f1)
          os << endl;
        os << "[";

        bool f2 = true;
        for (const auto & e : l)
          {
            if (!f2)
              os << ", ";
            os << e;
            f2 = false;
          }

        os << "]";
        f1 = false;
      }

    return os;
  }
  using nmatrix = NumMatrix<double>;

  nmatrix operator*(double, const nmatrix &);
  double max_eigen_value(const nmatrix &);

};

#endif // __NUMMATRIX_H_
