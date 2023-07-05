#include "Utilities/DBC.hh"

#include <Eigen/Dense>

namespace Spheral {

//------------------------------------------------------------------------------
// Construct to fit the given function -- we have to estimate the gradient at
// each point.
//------------------------------------------------------------------------------
template<typename Func>
inline
CubicHermiteInterpolator::CubicHermiteInterpolator(const double xmin,
                                                   const double xmax,
                                                   const size_t n,
                                                   const Func& F):
  mN(n),
  mXmin(xmin),
  mXmax(xmax),
  mXstep((xmax - xmin)/(n - 1u)),
  mVals(2u*n) {

  // Preconditions
  VERIFY2(n > 2u, "CubicHermiteInterpolator requires n >= 3 without a gradient function : n=" << n);
  VERIFY2(xmax > xmin, "CubicHermiteInterpolator requires a positive domain: [" << xmin << " " << xmax << "]");

  // Compute the function values
  for (auto i = 0u; i < mN; ++i) mVals[i] = F(xmin + i*mXstep);

  // Estimate the gradients at each interpolation node
  const auto dx = 0.001*mXstep;
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
CubicHermiteInterpolator::CubicHermiteInterpolator(const double xmin,
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
  VERIFY2(n > 1u, "CubicHermiteInterpolator requires n >= 2 : n=" << n);
  VERIFY2(xmax > xmin, "CubicHermiteInterpolator requires a positive domain: [" << xmin << " " << xmax << "]");

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
CubicHermiteInterpolator::operator()(const double x) const {
  const auto i0 = lowerBound(x);
  return this->operator()(x, i0);
}

inline
double
CubicHermiteInterpolator::operator()(const double x,
                                     const size_t i0) const {
  REQUIRE(i0 <= mN - 2u);
  const auto t = std::max(0.0, std::min(1.0, (x - mXmin - i0*mXstep)/mXstep));
  const auto t2 = t*t;
  const auto t3 = t*t2;
  return ((2.0*t3 - 3.0*t2 + 1.0)*mVals[i0] +          // h00
          (-2.0*t3 + 3.0*t2)*mVals[i0 + 1] +           // h01
          mXstep*((t3 - 2.0*t2 + t)*mVals[mN + i0] +   // h10
                  (t3 - t2)*mVals[mN + i0 + 1u]));     // h11
  // return (h00(t)*mVals[i0] + h01(t)*mVals[i0 + 1u] + 
  //         mXstep*(h10(t)*mVals[mN + i0] + h11(t)*mVals[mN + i0 + 1u]));
}

//------------------------------------------------------------------------------
// Return the lower bound entry in the table for the given x coordinate
//------------------------------------------------------------------------------
inline
size_t
CubicHermiteInterpolator::lowerBound(const double x) const {
  const auto result = std::min(mN - 2u, size_t(std::max(0.0, x - mXmin)/mXstep));
  ENSURE(result <= mN - 2u);
  return result;
}

//------------------------------------------------------------------------------
// Hermite basis functions
//------------------------------------------------------------------------------
inline
double
CubicHermiteInterpolator::h00(const double t) const {
  return (2.0*t - 3.0)*t*t + 1.0;
}

inline
double
CubicHermiteInterpolator::h10(const double t) const {
  return (t - 2.0)*t*t + t;
}

inline
double
CubicHermiteInterpolator::h01(const double t) const {
  return (3.0 - 2.0*t)*t*t;
}

inline
double
CubicHermiteInterpolator::h11(const double t) const {
  return (t - 1.0)*t*t;
}

//------------------------------------------------------------------------------
// Data accessors
//------------------------------------------------------------------------------
inline
size_t
CubicHermiteInterpolator::N() const {
  return mN;
}

inline
double
CubicHermiteInterpolator::xmin() const {
  return mXmin;
}

inline
double
CubicHermiteInterpolator::xmax() const {
  return mXmax;
}

inline
double
CubicHermiteInterpolator::xstep() const {
  return mXstep;
}

inline
const std::vector<double>&
CubicHermiteInterpolator::vals() const {
  return mVals;
}

}
