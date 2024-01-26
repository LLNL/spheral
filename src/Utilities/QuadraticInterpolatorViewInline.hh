#include "Utilities/DBC.hh"
#include <algorithm>

#include <Eigen/Dense>

namespace Spheral {

//------------------------------------------------------------------------------
// Default constructor
//------------------------------------------------------------------------------
inline
QuadraticInterpolatorView::QuadraticInterpolatorView():
  mN1(),
  mXmin(),
  mXmax(),
  mXstep(),
  mcoeffs() {
}

//------------------------------------------------------------------------------
// Interpolate for the given x value.
//------------------------------------------------------------------------------
inline
double
QuadraticInterpolatorView::operator()(const double x) const {
  const auto i0 = lowerBound(x);
  return mcoeffs[i0] + (mcoeffs[i0 + 1] + mcoeffs[i0 + 2]*x)*x;
}

inline
double
QuadraticInterpolatorView::operator()(const double x,
                                  const size_t i0) const {
  REQUIRE(i0 <= 3u*mN1);
  return mcoeffs[i0] + (mcoeffs[i0 + 1] + mcoeffs[i0 + 2]*x)*x;
}

//------------------------------------------------------------------------------
// Interpolate the first derivative the given x value.
//------------------------------------------------------------------------------
inline
double
QuadraticInterpolatorView::prime(const double x) const {
  const auto i0 = lowerBound(x);
  return mcoeffs[i0 + 1] + 2.0*mcoeffs[i0 + 2]*x;
}

inline
double
QuadraticInterpolatorView::prime(const double x,
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
QuadraticInterpolatorView::prime2(const double x) const {
  const auto i0 = lowerBound(x);
  return 2.0*mcoeffs[i0 + 2];
}

inline
double
QuadraticInterpolatorView::prime2(const double /*x*/,
                              const size_t i0) const {
  REQUIRE(i0 <= 3u*mN1);
  return 2.0*mcoeffs[i0 + 2];
}

//------------------------------------------------------------------------------
// Return the lower bound entry in the table for the given x coordinate
//------------------------------------------------------------------------------
inline
size_t
QuadraticInterpolatorView::lowerBound(const double x) const {
  const auto result = 3u*std::min(mN1, size_t(std::max(0.0, x - mXmin)/mXstep));
  ENSURE(result <= 3u*mN1);
  return result;
}

//------------------------------------------------------------------------------
// Data accessors
//------------------------------------------------------------------------------
inline
size_t
QuadraticInterpolatorView::size() const {
  return mcoeffs.size();
}

inline
double
QuadraticInterpolatorView::xmin() const {
  return mXmin;
}

inline
double
QuadraticInterpolatorView::xmax() const {
  return mXmax;
}

inline
double
QuadraticInterpolatorView::xstep() const {
  return mXstep;
}

inline
//const std::vector<double>&
const typename QuadraticInterpolatorView::CoeffsType&
QuadraticInterpolatorView::coeffs() const {
  return mcoeffs;
}

//------------------------------------------------------------------------------
// Destructor
//------------------------------------------------------------------------------
inline
QuadraticInterpolatorView::~QuadraticInterpolatorView() {
}

//------------------------------------------------------------------------------
// Equivalence
//------------------------------------------------------------------------------
inline
bool
QuadraticInterpolatorView::
operator==(const QuadraticInterpolatorView& rhs) const {
  return ((mN1 == rhs.mN1) and
          (mXmin == rhs.mXmin) and
          (mXmax == rhs.mXmax) and
          (mcoeffs == rhs.mcoeffs));
}

}
