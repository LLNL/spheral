#include "Utilities/DBC.hh"

namespace Spheral {

//------------------------------------------------------------------------------
// Interpolate for the given x value.
//------------------------------------------------------------------------------
inline
double
ParabolicInterpolator::operator()(const double x) const {
  const auto n = mA.size();
  const auto i0 = std::min(n - 3, lowerBound(x));
  const auto i1 = i0 + 1;
  CHECK(i1 >= 1 and i1 <= n - 2);
  const auto dx = x/mXstep - i0;
  CHECK(dx >= 0.0);
  return mA[i1] + mB[i1]*dx + mC[i1]*dx*dx;
}

//------------------------------------------------------------------------------
// Interpolate the first derivative the given x value.
//------------------------------------------------------------------------------
inline
double
ParabolicInterpolator::prime(const double x) const {
  const auto n = mA.size();
  const auto i0 = std::min(n - 3, lowerBound(x));
  const auto i1 = i0 + 1;
  CHECK(i1 >= 1 and i1 <= n - 2);
  const auto dx = x/mXstep - i0;
  CHECK(dx >= 0.0);
  return mB[i1] + mC[i1]*dx;
}

//------------------------------------------------------------------------------
// Interpolate the second derivative for the given x value.
//------------------------------------------------------------------------------
inline
double
ParabolicInterpolator::prime2(const double x) const {
  const auto n = mA.size();
  const auto i0 = std::min(n - 3, lowerBound(x));
  const auto i1 = i0 + 1;
  CHECK(i1 >= 1 and i1 <= n - 2);
  return mC[i1];
}

//------------------------------------------------------------------------------
// Return the lower bound entry in the table for the given x coordinate
//------------------------------------------------------------------------------
inline
size_t
ParabolicInterpolator::lowerBound(const double x) const {
  const auto n = mA.size();
  const auto result = std::min(n - 1, size_t(std::max(x, mXmin)/mXstep));
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
