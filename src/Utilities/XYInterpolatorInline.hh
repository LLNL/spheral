#include "Utilities/DBC.hh"
#include <iostream>
#include <cmath>

namespace Spheral {

//------------------------------------------------------------------------------
// Default constructor
//------------------------------------------------------------------------------
inline
XYInterpolator::XYInterpolator():
  mnx1(),
  mny1(),
  mxmin(),
  mxmax(),
  mymin(),
  mymax(),
  mxstep(),
  mystep() {
}

//------------------------------------------------------------------------------
// Construct with limits
//------------------------------------------------------------------------------
inline
XYInterpolator::XYInterpolator(const double xmin,
                               const double xmax,
                               const double ymin,
                               const double ymax,
                               const size_t nx,
                               const size_t ny,
                               const bool xlog,
                               const bool ylog):
  mxlog(xlog),
  mylog(ylog),
  mnx1(std::max(size_t(2u), nx) - 1u),
  mny1(std::max(size_t(2u), ny) - 1u),
  mxmin(xmin),
  mxmax(xmax),
  mymin(ymin),
  mymax(ymax),
  mxstep(),
  mystep() {
  // Figure out the sampling steps.
  mxstep = (xmax - xmin)/(xlog ? 1.0 : mnx1);
  mystep = (ymax - ymin)/(ylog ? 1.0 : mny1);
}

//------------------------------------------------------------------------------
// Destructor
//------------------------------------------------------------------------------
inline
XYInterpolator::~XYInterpolator() {
}

//------------------------------------------------------------------------------
// Return the lower bound entry in the table for the given (x,y) position
//------------------------------------------------------------------------------
inline
void
XYInterpolator::lowerBound(const double x, const double y,
                                    size_t& ix, size_t& iy, size_t& i) const {
  ix = (mxlog ?
        std::min(mnx1 - 1u, size_t(std::max(0.0, mnx1 + log(std::min(1.0, std::max(1e-10, (x - mxmin)/mxstep)))))) :
        std::min(mnx1 - 1u, size_t(std::max(0.0, x - mxmin)/mxstep)));
  iy = (mylog ?
        std::min(mny1 - 1u, size_t(std::max(0.0, mny1 + log(std::min(1.0, std::max(1e-10, (y - mymin)/mystep)))))) :
        std::min(mny1 - 1u, size_t(std::max(0.0, y - mymin)/mystep)));
  CHECK(ix < mnx1 and iy < mny1);
  i = 9u*(mnx1*iy + ix);
  ENSURE(i <= 9u*mnx1*mny1);
}

//------------------------------------------------------------------------------
// Data accessors
//------------------------------------------------------------------------------
inline
double
XYInterpolator::xmin() const {
  return mxmin;
}

inline
double
XYInterpolator::xmax() const {
  return mxmax;
}

inline
double
XYInterpolator::ymin() const {
  return mymin;
}

inline
double
XYInterpolator::ymax() const {
  return mymax;
}

inline
double
XYInterpolator::xstep() const {
  return mxstep;
}

inline
double
XYInterpolator::ystep() const {
  return mystep;
}

inline
bool
XYInterpolator::xlog() const {
  return mxlog;
}

inline
bool
XYInterpolator::ylog() const {
  return mylog;
}

//------------------------------------------------------------------------------
// operator==
//------------------------------------------------------------------------------
inline
bool
XYInterpolator::operator==(const XYInterpolator& rhs) const {
  return ((mnx1 == rhs.mnx1) and
          (mny1 == rhs.mny1) and
          (mxmin == rhs.mxmin) and
          (mxmax == rhs.mxmax) and
          (mymin == rhs.mymin) and
          (mymax == rhs.mymax) and
          (mxstep == rhs.mxstep) and
          (mystep == rhs.mystep));
}

//------------------------------------------------------------------------------
// Compute the coordinate value for a grid point
//------------------------------------------------------------------------------
inline
double
XYInterpolator::coord(const double xmin,
                      const double dx,
                      const size_t ix,
                      const size_t nx,
                      const bool xlog) const {
  REQUIRE(ix <= nx);
  if (xlog) {
    return xmin + dx*exp(int(ix) - int(nx));
  } else {
    return xmin + ix*dx;
  }
}

//------------------------------------------------------------------------------
// Similar to above, but compute the relative normalized coordinate inside
// a grid patch for the fit (range [0,1]).
//------------------------------------------------------------------------------
inline
void
XYInterpolator::eta_coords(const double xi,
                           const double yi,
                           double& etax,
                           double& etay,
                           size_t& ix,
                           size_t& iy,
                           size_t& i0) const {
  const auto xb = std::max(mxmin, std::min(mxmax, xi));
  const auto yb = std::max(mymin, std::min(mymax, yi));
  lowerBound(xb, yb, ix, iy, i0);
  const auto x0 = coord(mxmin, mxstep, ix, mnx1, mxlog);
  const auto y0 = coord(mymin, mystep, iy, mny1, mylog);
  const auto dx = (mxlog ? coord(mxmin, mxstep, ix + 1u, mnx1, mxlog) - x0 : mxstep);
  const auto dy = (mylog ? coord(mymin, mystep, iy + 1u, mny1, mylog) - y0 : mystep);
  etax = std::max(0.0, std::min(1.0, (xb - x0)/dx));
  etay = std::max(0.0, std::min(1.0, (yb - y0)/dy));
}

}
