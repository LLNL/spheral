#include "Utilities/DBC.hh"
#include <Eigen/Dense>
#include <iostream>
#include <cmath>

namespace Spheral {

//------------------------------------------------------------------------------
// Default constructor
//------------------------------------------------------------------------------
inline
BiQuadraticInterpolator::BiQuadraticInterpolator():
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
template<typename Func>
BiQuadraticInterpolator::BiQuadraticInterpolator(const double xmin,
                                                 const double xmax,
                                                 const double ymin,
                                                 const double ymax,
                                                 const size_t nx,
                                                 const size_t ny,
                                                 const Func& F,
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
  this->initialize(xmin, xmax, ymin, ymax, nx, ny, F);
}

//------------------------------------------------------------------------------
// Destructor
//------------------------------------------------------------------------------
inline
BiQuadraticInterpolator::~BiQuadraticInterpolator() {
}

//------------------------------------------------------------------------------
// Initialize the interpolation to fit the given data
//------------------------------------------------------------------------------
template<typename Func>
inline
void
BiQuadraticInterpolator::initialize(const double xmin,
                                    const double xmax,
                                    const double ymin,
                                    const double ymax,
                                    const size_t nx,
                                    const size_t ny,
                                    const Func& F,
                                    const bool xlog,
                                    const bool ylog) {

  // Size stuff up.
  REQUIRE(nx > 2u);
  REQUIRE(ny > 2u);
  mnx1 = nx - 1u;
  mny1 = ny - 1u;
  mxlog = xlog;
  mylog = ylog;
  mcoeffs.resize(9*mnx1*mny1);

  // Figure out the sampling steps.
  mxmin = xmin;
  mxmax = xmax;
  mymin = ymin;
  mymax = ymax;
  mxstep = (xmax - xmin)/mnx1;
  mystep = (ymax - ymin)/mny1;
  // mxstep = (xmax - xmin)/(xlog ? 1.0 : mnx1);
  // mystep = (ymax - ymin)/(ylog ? 1.0 : mny1);

  // Fit the coefficients
  double x1, x2, x3, y1, y2, y3;
  Eigen::MatrixXd A(9, 9);
  Eigen::VectorXd b(9), c(9);
  for (auto i = 0u; i < mnx1; ++i) {
    for (auto j = 0u; j < mny1; ++j) {
      x1 = xmin + i*mxstep;
      x2 = x1 + 0.5*mxstep;
      x3 = x1 +     mxstep;
      y1 = ymin + j*mystep;
      y2 = y1 + 0.5*mystep;
      y3 = y1 +     mystep;
      // x1 = coord(xmin, mxstep, i,      mnx1, xlog);
      // x3 = coord(xmin, mxstep, i + 1u, mnx1, xlog);
      // x2 = 0.5*(x1 + x3);
      // y1 = coord(ymin, mystep, j,      mny1, ylog);
      // y3 = coord(ymin, mystep, j + 1u, mny1, ylog);
      // y2 = 0.5*(y1 + y3);
      A << 1.0, x1, y1, x1*y1, x1*x1, x1*x1*y1, y1*y1, x1*y1*y1, x1*x1*y1*y1,
           1.0, x2, y1, x2*y1, x2*x2, x2*x2*y1, y1*y1, x2*y1*y1, x2*x2*y1*y1,
           1.0, x3, y1, x3*y1, x3*x3, x3*x3*y1, y1*y1, x3*y1*y1, x3*x3*y1*y1,
           1.0, x1, y2, x1*y2, x1*x1, x1*x1*y2, y2*y2, x1*y2*y2, x1*x1*y2*y2,
           1.0, x2, y2, x2*y2, x2*x2, x2*x2*y2, y2*y2, x2*y2*y2, x2*x2*y2*y2,
           1.0, x3, y2, x3*y2, x3*x3, x3*x3*y2, y2*y2, x3*y2*y2, x3*x3*y2*y2,
           1.0, x1, y3, x1*y3, x1*x1, x1*x1*y3, y3*y3, x1*y3*y3, x1*x1*y3*y3,
           1.0, x2, y3, x2*y3, x2*x2, x2*x2*y3, y3*y3, x2*y3*y3, x2*x2*y3*y3,
           1.0, x3, y3, x3*y3, x3*x3, x3*x3*y3, y3*y3, x3*y3*y3, x3*x3*y3*y3;
      b << F(x1, y1),
           F(x2, y1),
           F(x3, y1),
           F(x1, y2),
           F(x2, y2),
           F(x3, y2),
           F(x1, y3),
           F(x2, y3),
           F(x3, y3);
      CHECK2(b == b, "BiQuadraticInterpoolator function return error: \n"
             << "(" << x1 << " " << y1 << ") : " << b[0] << "\n"
             << "(" << x2 << " " << y1 << ") : " << b[1] << "\n"
             << "(" << x3 << " " << y1 << ") : " << b[2] << "\n"
             << "(" << x1 << " " << y2 << ") : " << b[3] << "\n"
             << "(" << x2 << " " << y2 << ") : " << b[4] << "\n"
             << "(" << x3 << " " << y2 << ") : " << b[5] << "\n"
             << "(" << x1 << " " << y3 << ") : " << b[6] << "\n"
             << "(" << x2 << " " << y3 << ") : " << b[7] << "\n"
             << "(" << x3 << " " << y3 << ") : " << b[8] << "\n");
      c = A.inverse()*b;
      // c = A.fullPivHouseholderQr().solve(b);
      // c = A.bdcSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(b);
      // std::cerr << "------------------------------------------------------------------------------\n"
      //           << "x00: " << x00 << "\n"
      //           << "x10: " << x10 << "\n"
      //           << "x20: " << x20 << "\n"
      //           << "x01: " << x01 << "\n"
      //           << "x11: " << x11 << "\n"
      //           << "x21: " << x21 << "\n"
      //           << "x02: " << x02 << "\n"
      //           << "x12: " << x12 << "\n"
      //           << "x22: " << x22 << "\n"
      //           << "A:\n" << A << "\n"
      //           << "b:\n" << b << "\n"
      //           << "c:\n" << c << "\n";
      const auto k = 9*(i + j*mnx1);
      mcoeffs[k    ] = c(0);
      mcoeffs[k + 1] = c(1);
      mcoeffs[k + 2] = c(2);
      mcoeffs[k + 3] = c(3);
      mcoeffs[k + 4] = c(4);
      mcoeffs[k + 5] = c(5);
      mcoeffs[k + 6] = c(6);
      mcoeffs[k + 7] = c(7);
      mcoeffs[k + 8] = c(8);
    }
  }
}

//------------------------------------------------------------------------------
// Interpolate for the given coordinate.
//------------------------------------------------------------------------------
inline
double
BiQuadraticInterpolator::operator()(const double xi,
                                    const double yi) const {
  const auto x = std::max(mxmin, std::min(mxmax, xi));
  const auto y = std::max(mymin, std::min(mymax, yi));
  const auto i0 = lowerBound(x, y);
  // std::cerr << "================================================================================\n"
  //           << "mxlog, mylog : " << mxlog << " " << mylog << "\n"
  //           << "pos   : " << pos << "\n"
  //           << "(x,y) : " << x << " " << y << "\n"
  //           << "i0    : " << i0 << "\n"
  //           << "coeffs: " << mcoeffs[i0] << " " << mcoeffs[i0 + 1] << " " << mcoeffs[i0 + 2] << " " << mcoeffs[i0 + 3] << " " << mcoeffs[i0 + 4] << " " << mcoeffs[i0 + 5] << "\n"
  //           << "F(x,y): " << mcoeffs[i0] + mcoeffs[i0 + 1]*x + mcoeffs[i0 + 2]*y + mcoeffs[i0 + 3]*x*y + mcoeffs[i0 + 4]*x*x + mcoeffs[i0 + 5]*y*y << "\n";
  return mcoeffs[i0] + mcoeffs[i0 + 1]*x + mcoeffs[i0 + 2]*y + mcoeffs[i0 + 3]*x*y + mcoeffs[i0 + 4]*x*x  + mcoeffs[i0 + 5]*x*x*y + mcoeffs[i0 + 6]*y*y + mcoeffs[i0 + 7]*x*y*y + mcoeffs[i0 + 8]*x*x*y*y;
}

//------------------------------------------------------------------------------
// Interpolate for the gradient (x)
//------------------------------------------------------------------------------
inline
double
BiQuadraticInterpolator::prime_x(const double xi, const double yi) const {
  const auto x = std::max(mxmin, std::min(mxmax, xi));
  const auto y = std::max(mymin, std::min(mymax, yi));
  const auto i0 = lowerBound(x, y);
  return mcoeffs[i0 + 1] + mcoeffs[i0 + 3]*y + 2.0*mcoeffs[i0 + 4]*x + 2.0*mcoeffs[i0 + 5]*x*y + mcoeffs[i0 + 7]*y*y + 2.0*mcoeffs[i0 + 8]*x*y*y;
}

//------------------------------------------------------------------------------
// Interpolate for the gradient (y)
//------------------------------------------------------------------------------
inline
double
BiQuadraticInterpolator::prime_y(const double xi, const double yi) const {
  const auto x = std::max(mxmin, std::min(mxmax, xi));
  const auto y = std::max(mymin, std::min(mymax, yi));
  const auto i0 = lowerBound(x, y);
  return mcoeffs[i0 + 2] + mcoeffs[i0 + 3]*x + mcoeffs[i0 + 5]*x*x + 2.0*mcoeffs[i0 + 6]*y + 2.0*mcoeffs[i0 + 7]*x*y + 2.0*mcoeffs[i0 + 8]*x*x*y;
}

//------------------------------------------------------------------------------
// Interpolate for the gradient2 (xx)
//------------------------------------------------------------------------------
inline
double
BiQuadraticInterpolator::prime2_xx(const double xi, const double yi) const {
  const auto x = std::max(mxmin, std::min(mxmax, xi));
  const auto y = std::max(mymin, std::min(mymax, yi));
  const auto i0 = lowerBound(x, y);
  return 2.0*(mcoeffs[i0 + 4] + mcoeffs[i0 + 5]*y + mcoeffs[i0 + 8]*y*y);
}

//------------------------------------------------------------------------------
// Interpolate for the gradient2 (xy)
//------------------------------------------------------------------------------
inline
double
BiQuadraticInterpolator::prime2_xy(const double xi, const double yi) const {
  const auto x = std::max(mxmin, std::min(mxmax, xi));
  const auto y = std::max(mymin, std::min(mymax, yi));
  const auto i0 = lowerBound(x, y);
  return mcoeffs[i0 + 3] + 2.0*mcoeffs[i0 + 5]*x + 2.0*mcoeffs[i0 + 7]*y + 4.0*mcoeffs[i0 + 8]*x*y;
}

//------------------------------------------------------------------------------
// Interpolate for the gradient2 (yx)
//------------------------------------------------------------------------------
inline
double
BiQuadraticInterpolator::prime2_yx(const double xi, const double yi) const {
  const auto x = std::max(mxmin, std::min(mxmax, xi));
  const auto y = std::max(mymin, std::min(mymax, yi));
  const auto i0 = lowerBound(x, y);
  return mcoeffs[i0 + 3] + 2.0*mcoeffs[i0 + 5]*x + 2.0*mcoeffs[i0 + 7]*y + 4.0*mcoeffs[i0 + 8]*x*y;
}

//------------------------------------------------------------------------------
// Interpolate for the gradient2 (yy)
//------------------------------------------------------------------------------
inline
double
BiQuadraticInterpolator::prime2_yy(const double xi, const double yi) const {
  const auto x = std::max(mxmin, std::min(mxmax, xi));
  const auto y = std::max(mymin, std::min(mymax, yi));
  const auto i0 = lowerBound(x, y);
  return 2.0*(mcoeffs[i0 + 6] + mcoeffs[i0 + 7]*x + mcoeffs[i0 + 8]*x*x);
}

//------------------------------------------------------------------------------
// Return the lower bound entry in the table for the given x coordinate
//------------------------------------------------------------------------------
inline
size_t
BiQuadraticInterpolator::lowerBound(const double x, const double y) const {
  const auto result = 9u*(mnx1*std::min(mny1 - 1u, size_t(std::max(0.0, y - mymin)/mystep)) +
                               std::min(mnx1 - 1u, size_t(std::max(0.0, x - mxmin)/mxstep)));
  // size_t ix, iy;
  // if (mxlog) {
  //   ix = std::min(mnx1 - 1u, size_t(mnx1 + log(std::min(1.0, std::max(0.0, (x - mxmin)/mxstep)))));
  // } else {
  //   ix = std::min(mnx1 - 1u, size_t(std::max(0.0, x - mxmin)/mxstep));
  // }
  // if (mylog) {
  //   iy = std::min(mny1 - 1u, size_t(mny1 + log(std::min(1.0, std::max(0.0, (y - mymin)/mystep)))));
  // } else {
  //   iy = std::min(mny1 - 1u, size_t(std::max(0.0, y - mymin)/mystep));
  // }
  // CHECK(ix < mnx1 and iy < mny1);
  // const auto result = 9u*(mnx1*iy + ix);
  ENSURE(result <= 9u*mnx1*mny1);
  return result;
}

//------------------------------------------------------------------------------
// Data accessors
//------------------------------------------------------------------------------
inline
size_t
BiQuadraticInterpolator::size() const {
  return mcoeffs.size();
}

inline
double
BiQuadraticInterpolator::xmin() const {
  return mxmin;
}

inline
double
BiQuadraticInterpolator::xmax() const {
  return mxmax;
}

inline
double
BiQuadraticInterpolator::ymin() const {
  return mymin;
}

inline
double
BiQuadraticInterpolator::ymax() const {
  return mymax;
}

inline
double
BiQuadraticInterpolator::xstep() const {
  return mxstep;
}

inline
double
BiQuadraticInterpolator::ystep() const {
  return mystep;
}

inline
bool
BiQuadraticInterpolator::xlog() const {
  return mxlog;
}

inline
bool
BiQuadraticInterpolator::ylog() const {
  return mylog;
}

inline
const std::vector<double>&
BiQuadraticInterpolator::coeffs() const {
  return mcoeffs;
}

//------------------------------------------------------------------------------
// operator==
//------------------------------------------------------------------------------
inline
bool
BiQuadraticInterpolator::operator==(const BiQuadraticInterpolator& rhs) const {
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
BiQuadraticInterpolator::coord(const double xmin,
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
