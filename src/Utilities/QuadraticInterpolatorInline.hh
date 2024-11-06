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
QuadraticInterpolator::QuadraticInterpolator(double xmin,
                                             double xmax,
                                             size_t n,
                                             const Func& F):
  mN1(),
  mXmin(),
  mXmax(),
  mXstep(),
  mcoeffs() {
  this->initialize(xmin, xmax, n, F);
}

//------------------------------------------------------------------------------
// Initialize to fit the given function
//------------------------------------------------------------------------------
template<typename Func>
inline
void
QuadraticInterpolator::initialize(double xmin,
                                  double xmax,
                                  size_t n,
                                  const Func& F) {
  // Preconditions
  VERIFY2(n > 1, "QuadraticInterpolator requires n > 1 : n=" << n);
  VERIFY2(xmax > xmin, "QuadraticInterpolator requires a positive domain: [" << xmin << " " << xmax << "]");

  // Build up an array of the function values and use the array based initialization.
  if (n % 2 == 0) ++n;  // Need odd number of samples to hit both endpoints of the range
  mXstep = (xmax - xmin)/(n - 1u);
  std::vector<double> yvals(n);
  for (auto i = 0u; i < n; ++i) yvals[i] = F(xmin + i*mXstep);
  this->initialize(xmin, xmax, yvals);
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
