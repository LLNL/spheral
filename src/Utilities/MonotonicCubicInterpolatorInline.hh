#include "Utilities/DBC.hh"

#include <Eigen/Dense>

namespace Spheral {

//------------------------------------------------------------------------------
// Construct to fit the given function -- we have to estimate the gradient at
// each point.
//------------------------------------------------------------------------------
template<typename Func>
inline
MonotonicCubicInterpolator::MonotonicCubicInterpolator(const double xmin,
                                                       const double xmax,
                                                       const size_t n,
                                                       const Func& F):
  mN(n),
  mXmin(xmin),
  mXmax(xmax),
  mXstep((xmax - xmin)/(n - 1u)),
  mVals(2u*n) {

  // Preconditions
  VERIFY2(n > 2u, "MonotonicCubicInterpolator requires n >= 3 without a gradient function : n=" << n);
  VERIFY2(xmax > xmin, "MonotonicCubicInterpolator requires a positive domain: [" << xmin << " " << xmax << "]");

  // Compute the function values
  for (auto i = 0u; i < mN; ++i) mVals[i] = F(xmin + i*mXstep);

  // Estimate the gradients at each interpolation node
  const auto dx = 0.05*mXstep;
  for (auto i = 0u; i < mN; ++i) {
    const auto xi = xmin + i*mXstep;
    mVals[mN + i] = (F(xi + dx) - F(xi - dx))/(2.0*dx);
  }
}

//------------------------------------------------------------------------------
// Construct to fit the given function with it's gradient
//------------------------------------------------------------------------------
template<typename Func, typename GradFunc>
inline
MonotonicCubicInterpolator::MonotonicCubicInterpolator(const double xmin,
                                                       const double xmax,
                                                       const size_t n,
                                                       const Func& F,
                                                       const GradFunc& Fgrad):
  mN(n),
  mXmin(xmin),
  mXmax(xmax),
  mXstep((xmax - xmin)/(n - 1u)),
  mVals(2u*n) {

  // Preconditions
  VERIFY2(n > 1u, "MonotonicCubicInterpolator requires n >= 2 : n=" << n);
  VERIFY2(xmax > xmin, "MonotonicCubicInterpolator requires a positive domain: [" << xmin << " " << xmax << "]");

  // Compute the function and gradient values
  for (auto i = 0u; i < mN; ++i) {
    const auto xi = xmin + i*mXstep;
    mVals[i] = F(xi);
    mVals[mN + i] = Fgrad(xi);
  }
}

//------------------------------------------------------------------------------
// Interpolate for the given x value.
//------------------------------------------------------------------------------
inline
double
MonotonicCubicInterpolator::operator()(const double x) const {
  const auto i0 = lowerBound(x);
  return this->operator()(x, i0);
}

inline
double
MonotonicCubicInterpolator::operator()(const double x,
                                       const size_t i0) const {
  REQUIRE(i0 <= mN - 2u);
  const auto t = std::max(0.0, std::min(1.0, (x - mXmin - i0*mXstep)/mXstep));
  return (h00(t)*mVals[i0] + h01(t)*mVals[i0 + 1u] + 
          mXstep*(h10(t)*mVals[mN + i0] + h11(t)*mVals[mN + i0 + 1u]));
}

//------------------------------------------------------------------------------
// Return the lower bound entry in the table for the given x coordinate
//------------------------------------------------------------------------------
inline
size_t
MonotonicCubicInterpolator::lowerBound(const double x) const {
  const auto result = std::min(mN - 2u, size_t(std::max(0.0, x - mXmin)/mXstep));
  ENSURE(result <= mN - 2u);
  return result;
}

//------------------------------------------------------------------------------
// Hermite basis functions
//------------------------------------------------------------------------------
inline
double
MonotonicCubicInterpolator::h00(const double t) const {
  return (2.0*t - 3.0)*t*t + 1.0;
}

inline
double
MonotonicCubicInterpolator::h10(const double t) const {
  return (t - 2.0)*t*t + t;
}

inline
double
MonotonicCubicInterpolator::h01(const double t) const {
  return (3.0 - 2.0*t)*t*t;
}

inline
double
MonotonicCubicInterpolator::h11(const double t) const {
  return (t - 1.0)*t*t;
}

//------------------------------------------------------------------------------
// Data accessors
//------------------------------------------------------------------------------
inline
size_t
MonotonicCubicInterpolator::N() const {
  return mN;
}

inline
double
MonotonicCubicInterpolator::xmin() const {
  return mXmin;
}

inline
double
MonotonicCubicInterpolator::xmax() const {
  return mXmax;
}

inline
double
MonotonicCubicInterpolator::xstep() const {
  return mXstep;
}

inline
const std::vector<double>&
MonotonicCubicInterpolator::vals() const {
  return mVals;
}

}
