#include "Utilities/DBC.hh"

namespace Spheral {

//------------------------------------------------------------------------------
// Interpolate for the given x value.
//------------------------------------------------------------------------------
inline
double
ParabolicInterpolator::operator()(const double x) const {
  const auto i0 = lowerBound(x);
  return mA[i0] + mB[i0]*x + mC[i0]*x*x;
}

//------------------------------------------------------------------------------
// Interpolate the first derivative the given x value.
//------------------------------------------------------------------------------
inline
double
ParabolicInterpolator::prime(const double x) const {
  const auto i0 = lowerBound(x);
  return mB[i0] + 2.0*mC[i0]*x;
}

//------------------------------------------------------------------------------
// Interpolate the second derivative for the given x value.
// Just a constant value, so not a great fit.
//------------------------------------------------------------------------------
inline
double
ParabolicInterpolator::prime2(const double x) const {
  const auto i0 = lowerBound(x);
  return 2.0*mC[i0];
}

//------------------------------------------------------------------------------
// Return the lower bound entry in the table for the given x coordinate
//------------------------------------------------------------------------------
inline
size_t
ParabolicInterpolator::lowerBound(const double x) const {
  const auto n = mA.size();
  CHECK(n > 0);
  const auto result = std::min(n - 1u, size_t(std::max(0.0, x - mXmin)/mXstep));
  ENSURE(result >= 0 && result < n);
  return result;
}

//------------------------------------------------------------------------------
// Data accessors
//------------------------------------------------------------------------------
inline
size_t
ParabolicInterpolator::size() const {
  return mA.size();
}

inline
double
ParabolicInterpolator::xmin() const {
  return mXmin;
}

inline
double
ParabolicInterpolator::xmax() const {
  return mXmax;
}

inline
double
ParabolicInterpolator::xstep() const {
  return mXstep;
}

inline
const std::vector<double>&
ParabolicInterpolator::a() const {
  return mA;
}

inline
const std::vector<double>&
ParabolicInterpolator::b() const {
  return mB;
}

inline
const std::vector<double>&
ParabolicInterpolator::c() const {
  return mC;
}



}
