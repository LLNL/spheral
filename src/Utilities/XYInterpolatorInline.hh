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
  mAx(),
  mBx(),
  mAy(),
  mBy(),
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
  mxstep(0.0),
  mystep(0.0),
  mAx(0.0),
  mBx(0.0),
  mAy(0.0),
  mBy(0.0),
  mcoeffs() {

  // Figure out the sampling steps.  Note A and B are for logarithmic stepping.
  mxstep = (xmax - xmin)/mnx1;
  mystep = (ymax - ymin)/mny1;
  mBx = (xmax - xmin)/(1.0 - exp(-double(mnx1)));
  mAx = xmax - mBx;
  mBy = (ymax - ymin)/(1.0 - exp(-double(mny1)));
  mAy = ymax - mBy;

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
        std::min(mnx1 - 1u, size_t(std::max(0.0, double(mnx1) + log((x - mAx)/mBx)))) :
        std::min(mnx1 - 1u, size_t(std::max(0.0, x - mxmin)/mxstep)));
  iy = (mylog ?
        std::min(mny1 - 1u, size_t(std::max(0.0, double(mny1) + log((y - mAy)/mBy)))) :
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
double
XYInterpolator::Ax() const {
  return mAx;
}  

inline
double
XYInterpolator::Bx() const {
  return mBx;
}  

inline
double
XYInterpolator::Ay() const {
  return mAy;
}  

inline
double
XYInterpolator::By() const {
  return mBy;
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
          (mAx == rhs.mAx) and 
          (mBx == rhs.mBx) and 
          (mAy == rhs.mAy) and 
          (mBy == rhs.mBy) and 
          (mcoeffs == rhs.mcoeffs));
}

//------------------------------------------------------------------------------
// Compute the coordinate value for a grid point
//------------------------------------------------------------------------------
inline
double
XYInterpolator::coord(const size_t ix,
                      const size_t nx,
                      const double xmin,
                      const double dx,
                      const double A,
                      const double B,
                      const bool xlog) const {
  REQUIRE(ix <= nx);
  return (xlog ?
          A + B*exp(double(ix) - double(nx)) :
          xmin + ix*dx);
}

//------------------------------------------------------------------------------
// Compute the x-coordinate for a grid point
//------------------------------------------------------------------------------
inline
double
XYInterpolator::xcoord(const size_t ix) const {
  return coord(ix, mnx1, mxmin, mxstep, mAx, mBx, mxlog);
}

//------------------------------------------------------------------------------
// Compute the y-coordinate for a grid point
//------------------------------------------------------------------------------
inline
double
XYInterpolator::ycoord(const size_t iy) const {
  return coord(iy, mny1, mymin, mystep, mAy, mBy, mylog);
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
  const auto dx = (mxlog ? xcoord(ix + 1u) - x0 : mxstep);
  const auto dy = (mylog ? ycoord(iy + 1u) - y0 : mystep);
  etax = std::max(0.0, std::min(1.0, (xb - x0)/dx));
  etay = std::max(0.0, std::min(1.0, (yb - y0)/dy));
  // std::cerr << "XYInterpolator::eta_coords: (" << xi << " " << yi << ") : (" << x0 << " " << y0 << ") : " << dx << " " << dy << " : (" << etax << " " << etay << ")" << std::endl;
}

}
