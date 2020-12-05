#include <iostream>
#include <lmscpp/nummatrix.hpp>


using namespace std;
using nmatrix = NumMatrix::NumMatrix<double>;

int main() {
  nmatrix a({{3,0,1},{2,2,2},{4,2,5}});
  nmatrix b({{1,0,0},{0,-2,0},{0,0,1.5}});

  cout << NumMatrix::max_eigen_value(a) << endl;
  cout << NumMatrix::max_eigen_value(b) << endl;
  return 0;
}
