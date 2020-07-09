#ifndef __NUMMATRIX_H_
#define __NUMMATRIX_H_

#include <vector>
namespace NumMatrix
{
  using namespace std;

  template <typename T>
  class NumMatrix
  {
    vector<vector<T>> m_;
    void throws_if_dimension_mismatch(const NumMatrix<T> &);
  public:
    NumMatrix(int nrows, int cols): m_{nrows,vector<T>(ncols, 0)} {}
    inline int nrows() const {return m_.size();}
    inline int ncols() const {return m_[0].size();}
    vector<T>& operator[](int idx) {return m_[idx];}

    auto begin() {return m_.begin();};
    auto end() {return m_.end();};
    auto begin() const {return m_.begin();};
    auto end() const {return m_.end();};

    NumMatrix<T>& operator*=(const NumMatrix<T> &);
    NumMatrix<T>& operator+=(const NumMatrix<T> &);
  };
};

#endif // __NUMMATRIX_H_
