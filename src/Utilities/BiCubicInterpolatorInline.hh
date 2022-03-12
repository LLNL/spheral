#include "Utilities/DBC.hh"
#include <Eigen/Dense>
#include <iostream>
#include <cmath>

namespace Spheral {

//------------------------------------------------------------------------------
// Default constructor
//------------------------------------------------------------------------------
inline
BiCubicInterpolator::BiCubicInterpolator():
  mnx1(),
  mny1(),
  mxmin(),
  mxmax(),
  mymin(),
  mymax(),
  mxstep(),
  mystep(),
  mcoeffs() {
}

//------------------------------------------------------------------------------
// Construct with tabulated data
//------------------------------------------------------------------------------
template<typename Func, typename GradFunc>
BiCubicInterpolator::BiCubicInterpolator(const double xmin,
                                         const double xmax,
                                         const double ymin,
                                         const double ymax,
                                         const size_t nx,
                                         const size_t ny,
                                         const Func& F,
                                         const GradFunc& gradF,
                                         const bool xlog,
                                         const bool ylog):
  mxlog(),
  mylog(),
  mnx1(),
  mny1(),
  mxmin(),
  mxmax(),
  mxstep(),
  mcoeffs() {
  this->initialize(xmin, xmax, ymin, ymax, nx, ny, F, gradF, xlog, ylog);
}

//------------------------------------------------------------------------------
// Destructor
//------------------------------------------------------------------------------
inline
BiCubicInterpolator::~BiCubicInterpolator() {
}

//------------------------------------------------------------------------------
// Initialize the interpolation to fit the given data
//------------------------------------------------------------------------------
template<typename Func, typename GradFunc>
inline
void
BiCubicInterpolator::initialize(const double xmin,
                                const double xmax,
                                const double ymin,
                                const double ymax,
                                const size_t nx,
                                const size_t ny,
                                const Func& F,
                                const GradFunc& gradF,
                                const bool xlog,
                                const bool ylog) {

  // Size stuff up.
  REQUIRE(nx > 1u);
  REQUIRE(ny > 1u);
  mnx1 = nx - 1u;
  mny1 = ny - 1u;
  mxlog = xlog;
  mylog = ylog;
  mcoeffs.resize(16*mnx1*mny1);

  // Figure out the sampling steps.
  mxmin = xmin;
  mxmax = xmax;
  mymin = ymin;
  mymax = ymax;
  mxstep = (xmax - xmin)/(xlog ? 1.0 : mnx1);
  mystep = (ymax - ymin)/(ylog ? 1.0 : mny1);

  // Fit the coefficients
  double x0, x1, y0, y1;
  Dim<2>::SymTensor gradF00, gradF01, gradF10, gradF11;
  Eigen::MatrixXd Ainv(16, 16);
  Eigen::VectorXd b(16), c(16);
  Ainv <<  1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
           0,  0,  0,  0,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
          -3,  3,  0,  0, -2, -1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
           2, -2,  0,  0,  1,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  
           0,  0,  0,  0,  0,  0,  0,  0,  1,  0,  0,  0,  0,  0,  0,  0,  
           0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  0,  0,  0,  
           0,  0,  0,  0,  0,  0,  0,  0, -3,  3,  0,  0, -2, -1,  0,  0,  
           0,  0,  0,  0,  0,  0,  0,  0,  2, -2,  0,  0,  1,  1,  0,  0,  
          -3,  0,  3,  0,  0,  0,  0,  0, -2,  0, -1,  0,  0,  0,  0,  0,  
           0,  0,  0,  0, -3,  0,  3,  0,  0,  0,  0,  0, -2,  0, -1,  0,  
           9, -9, -9,  9,  6,  3, -6, -3,  6, -6,  3, -3,  4,  2,  2,  1,
          -6,  6,  6, -6, -3, -3,  3,  3, -4,  4, -2,  2, -2, -2, -1, -1, 
           2,  0, -2,  0,  0,  0,  0,  0,  1,  0,  1,  0,  0,  0,  0,  0, 
           0,  0,  0,  0,  2,  0, -2,  0,  0,  0,  0,  0,  1,  0,  1,  0,  
          -6,  6,  6, -6, -4, -2,  4,  2, -3,  3, -3,  3, -2, -1, -2, -1, 
           4, -4, -4,  4,  2,  2, -2, -2,  2, -2,  2, -2,  1,  1,  1,  1;

  for (auto i = 0u; i < mnx1; ++i) {
    for (auto j = 0u; j < mny1; ++j) {
      x0 = coord(xmin, mxstep, i,      mnx1, xlog);
      x1 = coord(xmin, mxstep, i + 1u, mnx1, xlog);
      y0 = coord(ymin, mystep, j,      mny1, ylog);
      y1 = coord(ymin, mystep, j + 1u, mny1, ylog);
      gradF00 = gradF(x0, y0);
      gradF01 = gradF(x0, y1);
      gradF10 = gradF(x1, y0);
      gradF11 = gradF(x1, y1);
      b << F(x0, y0),
           F(x1, y0),
           F(x0, y1),
           F(x1, x1),
           gradF00.xx(),  // partial_x
           gradF10.xx(),  // partial_x
           gradF01.xx(),  // partial_x
           gradF11.xx(),  // partial_x
           gradF00.yy(),  // partial_y
           gradF10.yy(),  // partial_y
           gradF01.yy(),  // partial_y
           gradF11.yy(),  // partial_y
           gradF00.xy(),  // partial_xy
           gradF10.xy(),  // partial_xy
           gradF01.xy(),  // partial_xy
           gradF11.xy();  // partial_xy
      CHECK(b == b);
      c = Ainv*b;
      const auto k = 16*(i + j*mnx1);
      mcoeffs[k     ] = c(0);
      mcoeffs[k +  1] = c(1);
      mcoeffs[k +  2] = c(2);
      mcoeffs[k +  3] = c(3);
      mcoeffs[k +  4] = c(4);
      mcoeffs[k +  5] = c(5);
      mcoeffs[k +  6] = c(6);
      mcoeffs[k +  7] = c(7);
      mcoeffs[k +  8] = c(8);
      mcoeffs[k +  9] = c(9);
      mcoeffs[k + 10] = c(10);
      mcoeffs[k + 11] = c(11);
      mcoeffs[k + 12] = c(12);
      mcoeffs[k + 13] = c(13);
      mcoeffs[k + 14] = c(14);
      mcoeffs[k + 15] = c(15);
    }
  }
}

//------------------------------------------------------------------------------
// Interpolate for the given coordinate.
//------------------------------------------------------------------------------
inline
double
BiCubicInterpolator::operator()(const double xi,
                                const double yi) const {
  const auto xb = std::max(mxmin, std::min(mxmax, xi));
  const auto yb = std::max(mymin, std::min(mymax, yi));
  size_t ix, iy, i0;
  lowerBound(xb, yb, ix, iy, i0);
  const auto x = xb - coord(mxmin, mxstep, ix, mnx1, mxlog);
  const auto y = yb - coord(mymin, mystep, iy, mny1, mylog);
  const auto x2 = x*x, x3 = x*x*x;
  const auto y2 = y*y, y3 = y*y*y;
  return (mcoeffs[i0]               +          // c00
          mcoeffs[i0 +  1]*x        +          // c10
          mcoeffs[i0 +  2]*x2       +          // c20
          mcoeffs[i0 +  3]*x3       +          // c30
          mcoeffs[i0 +  4]    * y   +          // c01
          mcoeffs[i0 +  5]*x  * y   +          // c11
          mcoeffs[i0 +  6]*x2 * y   +          // c21
          mcoeffs[i0 +  7]*x3 * y   +          // c31
          mcoeffs[i0 +  8]    * y2  +          // c02
          mcoeffs[i0 +  9]*x  * y2  +          // c12
          mcoeffs[i0 + 10]*x2 * y2  +          // c22
          mcoeffs[i0 + 11]*x3 * y2  +          // c32
          mcoeffs[i0 + 12]    * y3  +          // c03
          mcoeffs[i0 + 13]*x  * y3  +          // c13
          mcoeffs[i0 + 14]*x2 * y3  +          // c23
          mcoeffs[i0 + 15]*x3 * y3);           // c33
}

//------------------------------------------------------------------------------
// Interpolate for the gradient (x)
//------------------------------------------------------------------------------
inline
double
BiCubicInterpolator::prime_x(const double xi, const double yi) const {
  return 0.0;
}

//------------------------------------------------------------------------------
// Interpolate for the gradient (y)
//------------------------------------------------------------------------------
inline
double
BiCubicInterpolator::prime_y(const double xi, const double yi) const {
  return 0.0;
}

//------------------------------------------------------------------------------
// Interpolate for the gradient2 (xx)
//------------------------------------------------------------------------------
inline
double
BiCubicInterpolator::prime2_xx(const double xi, const double yi) const {
  return 0.0;
}

//------------------------------------------------------------------------------
// Interpolate for the gradient2 (xy)
//------------------------------------------------------------------------------
inline
double
BiCubicInterpolator::prime2_xy(const double xi, const double yi) const {
  return 0.0;
}

//------------------------------------------------------------------------------
// Interpolate for the gradient2 (yy)
//------------------------------------------------------------------------------
inline
double
BiCubicInterpolator::prime2_yy(const double xi, const double yi) const {
  return 0.0;
}

//------------------------------------------------------------------------------
// Return the lower bound entry in the table for the given (x,y) position
//------------------------------------------------------------------------------
inline
void
BiCubicInterpolator::lowerBound(const double x, const double y,
                                size_t& ix, size_t& iy, size_t& i) const {
  ix = (mxlog ?
        std::min(mnx1 - 1u, size_t(mnx1 + log(std::min(1.0, std::max(0.0, (x - mxmin)/mxstep))))) :
        std::min(mnx1 - 1u, size_t(std::max(0.0, x - mxmin)/mxstep)));
  iy = (mylog ?
        std::min(mny1 - 1u, size_t(mny1 + log(std::min(1.0, std::max(0.0, (y - mymin)/mystep))))) :
        std::min(mny1 - 1u, size_t(std::max(0.0, y - mymin)/mystep)));
  CHECK(ix < mnx1 and iy < mny1);
  i = 16u*(mnx1*iy + ix);
  ENSURE(i <= 16u*mnx1*mny1);
}

//------------------------------------------------------------------------------
// Data accessors
//------------------------------------------------------------------------------
inline
size_t
BiCubicInterpolator::size() const {
  return mcoeffs.size();
}

inline
double
BiCubicInterpolator::xmin() const {
  return mxmin;
}

inline
double
BiCubicInterpolator::xmax() const {
  return mxmax;
}

inline
double
BiCubicInterpolator::ymin() const {
  return mymin;
}

inline
double
BiCubicInterpolator::ymax() const {
  return mymax;
}

inline
double
BiCubicInterpolator::xstep() const {
  return mxstep;
}

inline
double
BiCubicInterpolator::ystep() const {
  return mystep;
}

inline
bool
BiCubicInterpolator::xlog() const {
  return mxlog;
}

inline
bool
BiCubicInterpolator::ylog() const {
  return mylog;
}

inline
const std::vector<double>&
BiCubicInterpolator::coeffs() const {
  return mcoeffs;
}

//------------------------------------------------------------------------------
// operator==
//------------------------------------------------------------------------------
inline
bool
BiCubicInterpolator::operator==(const BiCubicInterpolator& rhs) const {
  return ((mnx1 == rhs.mnx1) and
          (mny1 == rhs.mny1) and
          (mxmin == rhs.mxmin) and
          (mxmax == rhs.mxmax) and
          (mymin == rhs.mymin) and
          (mymax == rhs.mymax) and
          (mxstep == rhs.mxstep) and
          (mystep == rhs.mystep) and
          (mcoeffs == rhs.mcoeffs));
}

//------------------------------------------------------------------------------
// Compute the coordinate value
//------------------------------------------------------------------------------
inline
double
BiCubicInterpolator::coord(const double xmin,
                           const double dx,
                           const size_t ix,
                           const size_t nx,
                           const bool xlog) const {
  REQUIRE(ix <= nx);
  if (xlog) {
    return xmin + dx*exp(ix - nx);
  } else {
    return xmin + ix*dx;
  }
}

}
