#include "Utilities/DBC.hh"
#include <iostream>
#include <cmath>

namespace Spheral {

//------------------------------------------------------------------------------
// Default constructor
//------------------------------------------------------------------------------
inline
XYInterpolator::XYInterpolator():
  mxlog(),
  mylog(),
  mnx1(),
  mny1(),
  mncoeffs(),
  mxmin(),
  mxmax(),
  mymin(),
  mymax(),
  mxstep(),
  mystep(),
  mcoeffs() {
}

//------------------------------------------------------------------------------
// Construct with limits
//------------------------------------------------------------------------------
inline
XYInterpolator::XYInterpolator(const size_t order,
                               const double xmin,
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
  mncoeffs((order + 1u)*(order + 1u)),
  mxmin(xmin),
  mxmax(xmax),
  mymin(ymin),
  mymax(ymax),
  mxstep(),
  mystep(),
  mcoeffs() {

  // Figure out the sampling steps.
  mxstep = (xmax - xmin)/(xlog ? 1.0 : mnx1);
  mystep = (ymax - ymin)/(ylog ? 1.0 : mny1);

  // Size coefficients array
  mcoeffs.resize(mncoeffs*mnx1*mny1);
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
  i = mncoeffs*(mnx1*iy + ix);
  ENSURE(i <= mncoeffs*mnx1*mny1);
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

inline
size_t
XYInterpolator::size() const {
  return mcoeffs.size();
}

inline
const std::vector<double>&
XYInterpolator::coeffs() const {
  return mcoeffs;
}

//------------------------------------------------------------------------------
// operator==
//------------------------------------------------------------------------------
inline
bool
XYInterpolator::operator==(const XYInterpolator& rhs) const {
  return ((mxlog == rhs.mxlog) and
          (mylog == rhs.mylog) and
          (mnx1 == rhs.mnx1) and
          (mny1 == rhs.mny1) and
          (mncoeffs == rhs.mncoeffs) and
          (mxmin == rhs.mxmin) and
          (mxmax == rhs.mxmax) and
          (mymin == rhs.mymin) and
          (mymax == rhs.mymax) and
          (mxstep == rhs.mxstep) and
          (mystep == rhs.mystep) and
          (mcoeffs == rhs.mcoeffs));
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
// Compute the x-coordinate for a grid point
//------------------------------------------------------------------------------
inline
double
XYInterpolator::xcoord(const size_t ix) const {
  return coord(mxmin, mxstep, ix, mnx1, mxlog);
}

//------------------------------------------------------------------------------
// Compute the y-coordinate for a grid point
//------------------------------------------------------------------------------
inline
double
XYInterpolator::ycoord(const size_t iy) const {
  return coord(mymin, mystep, iy, mny1, mylog);
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
  const auto x0 = xcoord(ix);
  const auto y0 = ycoord(iy);
  const auto dx = (mxlog ? coord(mxmin, mxstep, ix + 1u, mnx1, mxlog) - x0 : mxstep);
  const auto dy = (mylog ? coord(mymin, mystep, iy + 1u, mny1, mylog) - y0 : mystep);
  etax = std::max(0.0, std::min(1.0, (xb - x0)/dx));
  etay = std::max(0.0, std::min(1.0, (yb - y0)/dy));
  // std::cerr << "(" << xi << " " << yi << ") : (" << x0 << " " << y0 << ") : " << dx << " " << dy << " : (" << etax << " " << etay << ")" << std::endl;
}

}
