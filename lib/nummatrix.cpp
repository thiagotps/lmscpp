#include <lmscpp/nummatrix.hpp>


namespace NumMatrix {

nmatrix operator*(double v, const nmatrix &a) {
  nmatrix m{a.nrows(), a.ncols()};
  for (auto i = 0; i < a.nrows(); i++)
    for (auto j = 0; j < a.ncols(); j++)
      m.set(i, j, a.get(i, j) * v);

  return m;
}

// TODO: Move it to a separate file!
double max_eigen_value(const nmatrix &A) {
  if (A.ncols() != A.nrows())
    throw "This is not a square matrix";

  auto mnorm = [](const auto &z) {
    double res{-1.0};
    for (auto r = 0; r < z.nrows(); r++)
      res = max(res, abs(z.get(r, 0)));
    return res;
  };

  auto inner = [](const auto &a, const auto &b) {
    double res{0.0};
    auto n{a.nrows()};
    for (auto i = 0; i < n; i++)
        res += a.get(i,0) * b.get(i,0);
    return res;
  };

  constexpr double precision{pow(10.0, -5.0)};
  const int n = A.ncols();
  nmatrix y(n, 1);
  for (auto it = 0; it < n; it++)
    y.set(it, 0, 1);

  double l{0.0}, ll{0.0};
  bool flag;

  do {

    auto z = A * y;

    ll = l;
    l = inner(y, z)/inner(y,y);

    y = (1/mnorm(z))*z;

  } while (abs(l - ll) >= precision*abs(ll));


  return l;
}

} // namespace NumMatrix
