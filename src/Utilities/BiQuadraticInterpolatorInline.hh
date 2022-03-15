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
  XYInterpolator() {
}

//------------------------------------------------------------------------------
// Construct by sampling the given functor
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
  XYInterpolator(2u, xmin, xmax, ymin, ymax, nx, ny, xlog, ylog) {

  // We can predetermine A based on a unit square
  double x1 = 0.0, x2 = 0.5, x3 = 1.0;
  double y1 = 0.0, y2 = 0.5, y3 = 1.0;
  Eigen::MatrixXd A(9, 9);
  A << 1.0, x1, y1, x1*y1, x1*x1, x1*x1*y1, y1*y1, x1*y1*y1, x1*x1*y1*y1,
       1.0, x2, y1, x2*y1, x2*x2, x2*x2*y1, y1*y1, x2*y1*y1, x2*x2*y1*y1,
       1.0, x3, y1, x3*y1, x3*x3, x3*x3*y1, y1*y1, x3*y1*y1, x3*x3*y1*y1,
       1.0, x1, y2, x1*y2, x1*x1, x1*x1*y2, y2*y2, x1*y2*y2, x1*x1*y2*y2,
       1.0, x2, y2, x2*y2, x2*x2, x2*x2*y2, y2*y2, x2*y2*y2, x2*x2*y2*y2,
       1.0, x3, y2, x3*y2, x3*x3, x3*x3*y2, y2*y2, x3*y2*y2, x3*x3*y2*y2,
       1.0, x1, y3, x1*y3, x1*x1, x1*x1*y3, y3*y3, x1*y3*y3, x1*x1*y3*y3,
       1.0, x2, y3, x2*y3, x2*x2, x2*x2*y3, y3*y3, x2*y3*y3, x2*x2*y3*y3,
       1.0, x3, y3, x3*y3, x3*x3, x3*x3*y3, y3*y3, x3*y3*y3, x3*x3*y3*y3;
  const auto Ainv = A.inverse();
  
  // Fit the coefficients
  Eigen::VectorXd b(9), c(9);
  for (auto i = 0u; i < mnx1; ++i) {
    for (auto j = 0u; j < mny1; ++j) {
      x1 = coord(xmin, mxstep, i,      mnx1, xlog);
      x3 = coord(xmin, mxstep, i + 1u, mnx1, xlog);
      x2 = 0.5*(x1 + x3);
      y1 = coord(ymin, mystep, j,      mny1, ylog);
      y3 = coord(ymin, mystep, j + 1u, mny1, ylog);
      y2 = 0.5*(y1 + y3);
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
      c = Ainv*b;
      // c = A.bdcSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(b);
      // c = A.inverse()*b;
      // c = A.fullPivHouseholderQr().solve(b);
      // std::cerr << "------------------------------------------------------------------------------\n"
      //           << "x1 : " << x1 << "\n"
      //           << "x2 : " << x2 << "\n"
      //           << "x3 : " << x3 << "\n"
      //           << "y1 : " << y1 << "\n"
      //           << "y2 : " << y2 << "\n"
      //           << "y3 : " << y3 << "\n"
      //           << "A:\n" << A << "\n"
      //           << "b:\n" << b << "\n"
      //           << "c:\n" << c << "\n";
      const auto k = mncoeffs*(i + j*mnx1);
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
// Destructor
//------------------------------------------------------------------------------
inline
BiQuadraticInterpolator::~BiQuadraticInterpolator() {
}

//------------------------------------------------------------------------------
// Interpolate for the given coordinate.
//------------------------------------------------------------------------------
inline
double
BiQuadraticInterpolator::operator()(const double xi,
                                    const double yi) const {
  double x, y;
  size_t ix, iy, i0;
  eta_coords(xi, yi, x, y, ix, iy, i0);
  const auto x2 = x*x;
  const auto y2 = y*y;
  // std::cerr << "================================================================================\n"
  //           << "mxlog, mylog : " << mxlog << " " << mylog << "\n"
  //           << "pos   : " << pos << "\n"
  //           << "(x,y) : " << x << " " << y << "\n"
  //           << "i0    : " << i0 << "\n"
  //           << "coeffs: " << mcoeffs[i0] << " " << mcoeffs[i0 + 1] << " " << mcoeffs[i0 + 2] << " " << mcoeffs[i0 + 3] << " " << mcoeffs[i0 + 4] << " " << mcoeffs[i0 + 5] << "\n"
  //           << "F(x,y): " << mcoeffs[i0] + mcoeffs[i0 + 1]*x + mcoeffs[i0 + 2]*y + mcoeffs[i0 + 3]*x*y + mcoeffs[i0 + 4]*x*x + mcoeffs[i0 + 5]*y*y << "\n";
  return (mcoeffs[i0    ] +
          mcoeffs[i0 + 1]*x +
          mcoeffs[i0 + 2]*y +
          mcoeffs[i0 + 3]*x*y +
          mcoeffs[i0 + 4]*x2 +
          mcoeffs[i0 + 5]*x2*y +
          mcoeffs[i0 + 6]*y2 +
          mcoeffs[i0 + 7]*x*y2 +
          mcoeffs[i0 + 8]*x2*y2);
}

//------------------------------------------------------------------------------
// Interpolate for the gradient (x)
//------------------------------------------------------------------------------
inline
double
BiQuadraticInterpolator::prime_x(const double xi, const double yi) const {
  double x, y;
  size_t ix, iy, i0;
  eta_coords(xi, yi, x, y, ix, iy, i0);
  return (mcoeffs[i0 + 1] +
          mcoeffs[i0 + 3]*y +
          mcoeffs[i0 + 4]*2.0*x +
          mcoeffs[i0 + 5]*2.0*x*y +
          mcoeffs[i0 + 7]*y*y +
          mcoeffs[i0 + 8]*2.0*x*y*y);
}

//------------------------------------------------------------------------------
// Interpolate for the gradient (y)
//------------------------------------------------------------------------------
inline
double
BiQuadraticInterpolator::prime_y(const double xi, const double yi) const {
  double x, y;
  size_t ix, iy, i0;
  eta_coords(xi, yi, x, y, ix, iy, i0);
  return (mcoeffs[i0 + 2] +
          mcoeffs[i0 + 3]*x +
          mcoeffs[i0 + 5]*x*x +
          mcoeffs[i0 + 6]*2.0*y +
          mcoeffs[i0 + 7]*2.0*x*y +
          mcoeffs[i0 + 8]*2.0*x*x*y);
}

//------------------------------------------------------------------------------
// Interpolate for the gradient2 (xx)
//------------------------------------------------------------------------------
inline
double
BiQuadraticInterpolator::prime2_xx(const double xi, const double yi) const {
  double x, y;
  size_t ix, iy, i0;
  eta_coords(xi, yi, x, y, ix, iy, i0);
  return 2.0*(mcoeffs[i0 + 4] +
              mcoeffs[i0 + 5]*y +
              mcoeffs[i0 + 8]*y*y);
}

//------------------------------------------------------------------------------
// Interpolate for the gradient2 (xy)
//------------------------------------------------------------------------------
inline
double
BiQuadraticInterpolator::prime2_xy(const double xi, const double yi) const {
  double x, y;
  size_t ix, iy, i0;
  eta_coords(xi, yi, x, y, ix, iy, i0);
  return (mcoeffs[i0 + 3] +
          mcoeffs[i0 + 5]*2.0*x +
          mcoeffs[i0 + 7]*2.0*y +
          mcoeffs[i0 + 8]*4.0*x*y);
}

//------------------------------------------------------------------------------
// Interpolate for the gradient2 (yy)
//------------------------------------------------------------------------------
inline
double
BiQuadraticInterpolator::prime2_yy(const double xi, const double yi) const {
  double x, y;
  size_t ix, iy, i0;
  eta_coords(xi, yi, x, y, ix, iy, i0);
  return 2.0*(mcoeffs[i0 + 6] +
              mcoeffs[i0 + 7]*x +
              mcoeffs[i0 + 8]*x*x);
}

}
