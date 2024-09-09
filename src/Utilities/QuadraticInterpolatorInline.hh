#include "Utilities/DBC.hh"
#include <algorithm>

#include <Eigen/Dense>

namespace Spheral {

VVI_IMPL_BEGIN

//------------------------------------------------------------------------------
// Construct to fit the given function
//------------------------------------------------------------------------------
template<typename Func>
inline
QuadraticInterpolator::QuadraticInterpolator(const double xmin,
                                             const double xmax,
                                             const size_t n,
                                             const Func& F):
  mN1(n - 1),
  mXmin(xmin),
  mXmax(xmax),
  mXstep((xmax - xmin)/n),
  mcoeffs(3u*n) {

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

//------------------------------------------------------------------------------
// Equivalence
//------------------------------------------------------------------------------
bool
inline
QuadraticInterpolator::operator==(const QuadraticInterpolator& rhs) const {
  return ((mN1 == rhs.mN1) and
          (mXmin == rhs.mXmin) and
          (mXmax == rhs.mXmax) and
          (mcoeffs == rhs.mcoeffs));
}

//------------------------------------------------------------------------------
// Interpolate for the given x value.
//------------------------------------------------------------------------------
inline
double
QuadraticInterpolator::operator()(const double x) const {
  const auto i0 = lowerBound(x);
  return mcoeffs[i0] + (mcoeffs[i0 + 1] + mcoeffs[i0 + 2]*x)*x;
}

inline
double
QuadraticInterpolator::operator()(const double x,
                                  const size_t i0) const {
  REQUIRE(i0 <= 3u*mN1);
  return mcoeffs[i0] + (mcoeffs[i0 + 1] + mcoeffs[i0 + 2]*x)*x;
}

//------------------------------------------------------------------------------
// Interpolate the first derivative the given x value.
//------------------------------------------------------------------------------
inline
double
QuadraticInterpolator::prime(const double x) const {
  const auto i0 = lowerBound(x);
  return mcoeffs[i0 + 1] + 2.0*mcoeffs[i0 + 2]*x;
}

inline
double
QuadraticInterpolator::prime(const double x,
                             const size_t i0) const {
  REQUIRE(i0 <= 3u*mN1);
  return mcoeffs[i0 + 1] + 2.0*mcoeffs[i0 + 2]*x;
}

//------------------------------------------------------------------------------
// Interpolate the second derivative for the given x value.
// Just a constant value, so not a great fit.
//------------------------------------------------------------------------------
inline
double
QuadraticInterpolator::prime2(const double x) const {
  const auto i0 = lowerBound(x);
  return 2.0*mcoeffs[i0 + 2];
}

inline
double
QuadraticInterpolator::prime2(const double /*x*/,
                              const size_t i0) const {
  REQUIRE(i0 <= 3u*mN1);
  return 2.0*mcoeffs[i0 + 2];
}

//------------------------------------------------------------------------------
// Return the lower bound entry in the table for the given x coordinate
//------------------------------------------------------------------------------
inline
size_t
QuadraticInterpolator::lowerBound(const double x) const {
  const auto result = 3u*RAJA_MIN(mN1, size_t(RAJA_MAX(0.0, x - mXmin)/mXstep));
  ENSURE(result <= 3u*mN1);
  return result;
}

//------------------------------------------------------------------------------
// Data accessors
//------------------------------------------------------------------------------
inline
size_t
QuadraticInterpolator::size() const {
  return mcoeffs.size();
}

inline
double
QuadraticInterpolator::xmin() const {
  return mXmin;
}

inline
double
QuadraticInterpolator::xmax() const {
  return mXmax;
}

inline
double
QuadraticInterpolator::xstep() const {
  return mXstep;
}

inline
const vvi::vector<double>&
QuadraticInterpolator::coeffs() const {
  return mcoeffs;
}

VVI_IMPL_END

}
