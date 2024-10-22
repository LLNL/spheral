//---------------------------------Spheral++----------------------------------//
// QuadraticInterpolator
//
// Encapsulates the algorithm and data for parabolic interpolation in 1D
// Assumes the results is interpolated as y_interp = a + b*x + c*x^2
//
// Created by JMO, Fri Dec  4 14:28:08 PST 2020
//----------------------------------------------------------------------------//
#include "QuadraticInterpolator.hh"

namespace Spheral {

VVI_IMPL_BEGIN

//------------------------------------------------------------------------------
// Initialize the interpolation to fit the given data
//------------------------------------------------------------------------------
void
QuadraticInterpolator::initialize(const double xmin,
                                  const double xmax,
                                  const std::vector<double>& yvals) {
  const auto n = yvals.size();
  VERIFY2(n > 2, "QuadraticInterpolator::initialize requires at least 3 unique values to fit");
  VERIFY2(n % 2 == 1, "QuadraticInterpolator::initialize requires an odd number of tabulated values");
  VERIFY2(xmax > xmin, "QuadraticInterpolator::initialize requires a positive domain: [" << xmin << " " << xmax << "]");

  mN1 = (n - 1u)/2u - 1u;  // Maximum index into arrays
  mXmin = xmin;
  mXmax = xmax;
  mXstep = (xmax - xmin)/(mN1 + 1u);
  mcoeffs.resize(3*(mN1 + 1u));

  typedef Eigen::Matrix<double, 3, 3, Eigen::RowMajor> EMatrix;
  typedef Eigen::Matrix<double, 3, 1> EVector;

  // We use simple least squares fitting for each 3-point interval, giving us an
  // exact parabolic fit for those three points.
  // Find the coefficient fits.
  double x0, x1, x2;
  EMatrix A;
  EVector X, B;
  for (auto i0 = 0u; i0 <= mN1; ++i0) {
    x0 = xmin + i0*mXstep;
    x1 = x0 + 0.5*mXstep;
    x2 = x0 + mXstep;
    A << 1.0, x0, x0*x0,
         1.0, x1, x1*x1,
         1.0, x2, x2*x2;
    B << yvals[2u*i0], yvals[2u*i0 + 1u], yvals[2u*i0 + 2u];
    X = A.inverse()*B;
    mcoeffs[3*i0     ] = X(0);
    mcoeffs[3*i0 + 1u] = X(1);
    mcoeffs[3*i0 + 2u] = X(2);
  }
}

VVI_IMPL_END

} // namespace Spheral
