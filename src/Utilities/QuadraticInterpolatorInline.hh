#include "Utilities/DBC.hh"
#include <algorithm>

#include <Eigen/Dense>

namespace Spheral {
  
//------------------------------------------------------------------------------
// Construct to fit the given function
//------------------------------------------------------------------------------
template<typename Func>
inline
QuadraticInterpolator::QuadraticInterpolator(const double xmin,
                                             const double xmax,
                                             const size_t n,
                                             const Func& F) {
  //mN1(n - 1),
  //mXmin(xmin),
  //mXmax(xmax),
  //mXstep((xmax - xmin)/n),
  //mcoeffs(3u*n)
  mN1=n - 1;
  mXmin=xmin;
  mXmax=xmax;
  mXstep= (xmax - xmin)/n;
  mcoeffs = CoeffsType(3u*n);

  // Preconditions
  VERIFY2(n > 0, "QuadraticInterpolator requires n > 1 : n=" << n);
  VERIFY2(xmax > xmin, "QuadraticInterpolator requires a positive domain: [" << xmin << " " << xmax << "]");

  typedef Eigen::Matrix<double, 3, 3, Eigen::RowMajor> EMatrix;
  typedef Eigen::Matrix<double, 3, 1> EVector;

  // We use simple least squares fitting for each 3-point interval, giving us an
  // exact parabolic fit for those three points.
  // Find the coefficient fits.
  double x0, x1, x2;
  EMatrix A;
  EVector X, B;
  for (auto i0 = 0u; i0 < n; ++i0) {
    x0 = xmin + i0*mXstep;
    x1 = x0 + 0.5*mXstep;
    x2 = x0 + mXstep;
    A << 1.0, x0, x0*x0,
         1.0, x1, x1*x1,
         1.0, x2, x2*x2;
    B << F(x0), F(x1), F(x2);
    X = A.inverse()*B;
    mcoeffs[3*i0     ] = X(0);
    mcoeffs[3*i0 + 1u] = X(1);
    mcoeffs[3*i0 + 2u] = X(2);
  }
}

}
